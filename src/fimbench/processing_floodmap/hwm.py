"""
Authors: Supath Dhital (sdhital@ua.edu), Dipsikha Devi (ddevi@ua.edu)
Updated: June 2026

High Water Mark (HWM) flood map processing.

Standardizes a single USGS High Water Mark flood raster into the FIM database
format: renamed GeoTIFF, metadata.json, and AOI.gpkg, using a user-supplied
flood start and end date.
"""

import re
import json
import shutil
import tempfile
from pathlib import Path
import numpy as np
import rasterio
from rasterio.features import shapes
from rasterio.warp import calculate_default_transform, reproject, Resampling
from shapely.geometry import shape
from shapely.ops import unary_union
from shapely.validation import make_valid
import geopandas as gpd

from .._log import log as _log
from .utils import get_intersected_huc8, list_input_tifs, flood_polygon, write_extent_geojson
from .recurrence_interval import (
    get_event_recurrence_interval_nwm_range,
    get_event_recurrence_interval_usgs_range,
)


class HwmProcessor:
    """Standardize a High Water Mark flood map for the FIM database."""

    SENSOR_CODE = "HWM"
    Full_Form = "High Water Mark"
    SOURCE = "Dr. Dinuke Munasinghe,The University of Alabama, Email id: dsmunasinghe@ua.edu"

    NODATA_VAL = -9999
    BLOCK_SIZE = 512
    SIMPLIFY_TOL = 0.0001  # ~11 m at the equator

    # Extent GeoJSON is simplified by this many metres (EPSG:5070) to stay
    # light; the metadata/AOI geometry is not affected. None disables it.
    GEOJSON_SIMPLIFY_TOL = 10.0

    # metadata.json fields copied onto the extent GeoJSON feature.
    GEOJSON_FIELDS = (
        "File_Name",
        "Full form of the sensor code",
        "HUC8",
        "Resolution in meter",
        "Extent",
        "State",
        "River Basin Name",
        "Start Date of the Flood",
        "End Date of the Flood",
        "Event Recurrence Interval (USGS)",
        "Event Discharge (USGS) [cms]",
        "USGS Gauge ID",
        "Event Recurrence Interval (NWM)",
        "Event Discharge (NWM) [cms]",
    )

    # Reference data for the Event Recurrence Interval lookups (see
    # recurrence_interval.py). None disables the corresponding lookup.
    GAGES_SHP = None
    FLOW_SHP = None
    RECURR_FOLDER = None
    DOWNSTREAM_GAGE_MAX_KM = 100.0

    def __init__(
        self,
        *,
        sensor_code=None,
        full_form=None,
        source=None,
        nodata_val=None,
        block_size=None,
        simplify_tol=None,
        gages_shp=None,
        flow_shp=None,
        recurr_folder=None,
        downstream_gage_max_km=None,
    ):
        # Class attributes are the defaults.
        if sensor_code is not None:
            self.SENSOR_CODE = sensor_code
        if full_form is not None:
            self.Full_Form = full_form
        if source is not None:
            self.SOURCE = source
        if nodata_val is not None:
            self.NODATA_VAL = nodata_val
        if block_size is not None:
            self.BLOCK_SIZE = block_size
        if simplify_tol is not None:
            self.SIMPLIFY_TOL = simplify_tol
        if gages_shp is not None:
            self.GAGES_SHP = gages_shp
        if flow_shp is not None:
            self.FLOW_SHP = flow_shp
        if recurr_folder is not None:
            self.RECURR_FOLDER = recurr_folder
        if downstream_gage_max_km is not None:
            self.DOWNSTREAM_GAGE_MAX_KM = downstream_gage_max_km

        # Folder for intermediate rasters; set per run in process().
        self._work_dir = None

    def log(self, msg):
        _log(msg)

    def _work_path(self, name):
        # Place an intermediate file in the work dir (else next to the input).
        if self._work_dir is not None:
            return Path(self._work_dir) / name
        return Path(name)

    # Two consecutive YYYYMMDD dates joined by "_" or "-": the flood window
    # in an HWM name, e.g. HWM_10_0m_20100329_20100401_710139W415659N_BM.tif
    _YMD = r"(?:19|20)\d{2}(?:0[1-9]|1[0-2])(?:0[1-9]|[12]\d|3[01])"
    _WINDOW_RE = re.compile(rf"(?<!\d)({_YMD})[_-]({_YMD})(?!\d)")

    def infer_flood_window_from_path(self, path: Path):
        """
        Extract (start_date, end_date) as "YYYYMMDD" strings from the file
        name, else from its parent folders. Returns (None, None) if no
        YYYYMMDD_YYYYMMDD pair is found.
        """
        path = Path(path)
        for text in (path.name, path.as_posix()):
            m = self._WINDOW_RE.search(text)
            if m:
                return m.group(1), m.group(2)
        return None, None

    # Assign nodata, datatype and tile size
    def process_tif(self, input_tif, nodata_val=None, block_size=None):
        nodata_val = self.NODATA_VAL if nodata_val is None else nodata_val
        block_size = self.BLOCK_SIZE if block_size is None else block_size
        input_tif = Path(input_tif)
        output_tif = self._work_path(input_tif.stem + "_nodata.tif")
        output_tif.parent.mkdir(parents=True, exist_ok=True)

        with rasterio.open(input_tif) as src:
            m = src.read(masked=True)
            src_dtype = np.dtype(src.dtypes[0])
            if np.issubdtype(src_dtype, np.floating):
                out_dtype = "float32"
            else:
                out_dtype = "int32"

            out = np.full(m.shape, nodata_val, dtype=out_dtype)
            valid = ~m.mask
            out[valid] = m.data[valid].astype(out_dtype, copy=False)

            profile = src.profile.copy()
            profile.update(
                dtype=out_dtype,
                nodata=nodata_val,
                tiled=True,
                blockxsize=block_size,
                blockysize=block_size,
                compress="lzw",
            )

            with rasterio.open(output_tif, "w", **profile) as dst:
                dst.write(out)
                try:
                    dst.update_tags(**src.tags())
                except Exception:
                    pass

        return output_tif

    # Reproject to the target coordinate system (default EPSG:5070)
    def reproject_raster(self, input_file, target_epsg):
        input_file = Path(input_file)
        output_file = self._work_path(input_file.stem + f"_reprojected_{target_epsg}.tif")
        output_file.parent.mkdir(parents=True, exist_ok=True)

        with rasterio.open(input_file) as src:
            transform, width, height = calculate_default_transform(
                src.crs, f"EPSG:{target_epsg}", src.width, src.height, *src.bounds
            )

            metadata = src.meta.copy()
            metadata.update(
                {
                    "crs": f"EPSG:{target_epsg}",
                    "transform": transform,
                    "width": width,
                    "height": height,
                    "compress": "lzw",
                }
            )

            with rasterio.open(output_file, "w", **metadata) as dst:
                for i in range(1, src.count + 1):
                    reproject(
                        source=rasterio.band(src, i),
                        destination=rasterio.band(dst, i),
                        src_transform=src.transform,
                        src_crs=src.crs,
                        dst_transform=transform,
                        dst_crs=f"EPSG:{target_epsg}",
                        resampling=Resampling.nearest,
                    )
        return output_file

    # Resample the raster to a target resolution (default 10m) for visualization
    def resample_raster(self, in_file, target_res=10.0):
        in_file = Path(in_file)
        out_file = self._work_path(in_file.stem + f"_{int(target_res)}m.tif")

        with rasterio.open(in_file) as src:
            # New dimensions from the target resolution (scale = current_res / target_res)
            scale_factor_x = src.res[0] / target_res
            scale_factor_y = src.res[1] / target_res

            new_width = int(src.width * scale_factor_x)
            new_height = int(src.height * scale_factor_y)

            # Update the transform to reflect the new resolution
            new_transform = src.transform * src.transform.scale(
                (src.width / new_width), (src.height / new_height)
            )

            metadata = src.meta.copy()
            metadata.update(
                {
                    "driver": "GTiff",
                    "height": new_height,
                    "width": new_width,
                    "transform": new_transform,
                    "compress": "lzw",
                }
            )

            with rasterio.open(out_file, "w", **metadata) as dst:
                for i in range(1, src.count + 1):
                    reproject(
                        source=rasterio.band(src, i),
                        destination=rasterio.band(dst, i),
                        src_transform=src.transform,
                        src_crs=src.crs,
                        dst_transform=new_transform,
                        dst_crs=src.crs,
                        resampling=Resampling.bilinear,  # bilinear for continuous data, nearest for discrete
                    )
                final_crs = dst.crs
        return out_file, final_crs

    # Extract the resolution from the projected file
    def get_resolution_string(self, geo_tif_5070):
        with rasterio.open(geo_tif_5070) as src:
            res_x, res_y = src.res
            avg_res = (res_x + res_y) / 2.0
            rounded_res = round(avg_res, 1)
            res_str = f"{str(rounded_res).replace('.', '_')}m"  # e.g., 10.0 -> "10_0m"
            res_txt = f"{rounded_res:.2f} m"
        return rounded_res, res_str, res_txt

    # Extract the centroid from the raster
    def get_raster_centroid(self, out_4326):
        with rasterio.open(out_4326) as src:
            b = src.bounds
            centroid_x = (b.left + b.right) / 2
            centroid_y = (b.top + b.bottom) / 2
            crs = src.crs
            data = src.read(1)

        return centroid_x, centroid_y, crs, b

    # Decimal degrees to DMS, e.g. 95 31'44" W and 31 04'36" N -> 953144W310436N
    def decimal_to_dms_str(self, dec, is_lat=True):
        direction = ""
        if is_lat:
            direction = "N" if dec >= 0 else "S"
        else:
            direction = "E" if dec >= 0 else "W"

        dec = abs(dec)
        degrees = int(dec)
        minutes = int((dec - degrees) * 60)
        seconds = int(((dec - degrees) * 60 - minutes) * 60)

        return f"{degrees:02d}{minutes:02d}{seconds:02d}{direction}"

    def build_metadata_dict(
        self,
        src_4326,
        src_5070,
        new_filename,
        res_m_float,
        res_txt,
        dms_code,
        x,
        y,
        state_list,
        huc8_list,
        name_list,
        start_date,
        end_date,
    ):
        with rasterio.open(src_4326) as src:
            meta = {
                "File_Name": new_filename,
                "Full form of the sensor code": self.Full_Form,
                "HUC8": huc8_list,
                "Format": src.driver,
                "Nodata_Value": src.nodata,
                "Resolution in meter": res_m_float,
                "Datatype": src.dtypes[0],
                "Compression Type": src.tags(ns="IMAGE_STRUCTURE"),
                "Extent": {
                    "xmin": src.bounds.left,
                    "ymin": src.bounds.bottom,
                    "xmax": src.bounds.right,
                    "ymax": src.bounds.top,
                },
                "DMS_Code_centroid": dms_code,
                "Projection": src.crs.to_string() if src.crs else "Unknown",
                "Bands": src.count,
                "Rows": src.height,
                "Columns": src.width,
                "State": f"{state_list}, USA" if state_list else "USA",
                "Description": f"",
                "River Basin Name": name_list,
                "Source": self.SOURCE,
                "Location of the centroid of the flood map": (x, y),
                "Start Date of the Flood": start_date,
                "End Date of the Flood": end_date,
                "Event Recurrence Interval (USGS)": None,
                "Event Discharge (USGS) [cms]": None,
                "USGS Gauge ID": None,
                "Event Recurrence Interval (NWM)": None,
                "Event Discharge (NWM) [cms]": None,
                "Keywords": ["flood", "hazard", "GIS"],
                "Access_Rights": "Public",
                "References": [
                    "1. Cohen, S., Baruah, A., Nikrou, P., Tian, D., & Liu, H. (2025). Toward robust evaluations of flood inundation predictions using remote sensing derived benchmark maps. Water Resources Research, 61(8), e2024WR039574.",
                    "2. Munasinghe, D., Cohen, S., Tian, D., Liu, H., Baruah, A., and Devi, D. (2025). A Database of Flood Maps using high-resolution Airborne Imagery and Machine Learning Models. In: CIROH Developers' Conference, May 28-30, 2025.",
                    "3. Devi, D., Dhital, S., Munasinghe, D., Cohen, S., Baruah, A., Chen, Y., ... & Pruitt, C. (2025). A Framework for the Evaluation of Flood Inundation Predictions Over Extensive Benchmark Databases.",
                    "4. Tian.D, Liu.H,  Wang.L, Cohen.S, Thapa.P (2024).Enhancing satellite image-derived flood maps with hydrologically guided region growing method and high-resolution DEMs.Chapman Conference on Remote Sensing of the Water Cycle",
                ],
            }

        # Flood polygon from the 5070 raster: every pixel > 0, dissolved
        with rasterio.open(src_5070) as src_geom:
            geom_crs = src_geom.crs
            unified_geom = flood_polygon(
                src_geom.read(1, masked=True), src_geom.transform, self.SIMPLIFY_TOL
            )
        if unified_geom is None:
            self.log("No flooded areas found in 5070 raster.")
            meta["FIM_Geometry"] = None
            return meta, None

        # Event Recurrence Interval, computed against the unified flood
        # polygon: USGS (max observed gage discharge over [start_date,
        # end_date] vs. published FFA) and NWM (max modeled outlet discharge
        # over the same window vs. NWM v2.1 RI tables), plus the discharge
        # (cms) that produced each. Unlike Tier 1/2/3 (a single event
        # instant), HWM uses the MAXIMUM discharge anywhere in the
        # user-supplied flood window -- see recurrence_interval's *_range
        # functions. USGS RI is "No USGS Gauge" when no gage was found and
        # "No usable USGS Gauge" when none could be used; a field is None (logged, not raised) when its reference
        # data isn't configured, or the lookup couldn't be attempted at all.
        ri_usgs, discharge_usgs, _scenario_usgs, gauge_id_usgs = (
            get_event_recurrence_interval_usgs_range(
                unified_geom,
                geom_crs,
                start_date,
                end_date,
                gages_shp=self.GAGES_SHP,
                flow_shp=self.FLOW_SHP,
                downstream_max_km=self.DOWNSTREAM_GAGE_MAX_KM,
                log=self.log,
            )
        )
        ri_nwm, discharge_nwm, _scenario_nwm = get_event_recurrence_interval_nwm_range(
            unified_geom,
            geom_crs,
            start_date,
            end_date,
            flow_shp=self.FLOW_SHP,
            recurr_folder=self.RECURR_FOLDER,
            log=self.log,
        )
        meta["Event Recurrence Interval (USGS)"] = ri_usgs
        meta["Event Discharge (USGS) [cms]"] = (
            round(discharge_usgs, 2) if discharge_usgs is not None else None
        )
        meta["USGS Gauge ID"] = gauge_id_usgs
        meta["Event Recurrence Interval (NWM)"] = ri_nwm
        meta["Event Discharge (NWM) [cms]"] = (
            round(discharge_nwm, 2) if discharge_nwm is not None else None
        )

        return meta, unified_geom

    def raster_to_polygon(
        self, in_raster, gpkg_path, metadata_row, valid_classes=(0, 1), nodata=None, connectivity=8
    ):
        nodata = self.NODATA_VAL if nodata is None else nodata
        in_raster = Path(in_raster)
        with rasterio.open(in_raster) as src:
            arr = src.read(1)
            tfm = src.transform
            crs = src.crs
            mask = arr != nodata
            polygons = []
            for geom, val in shapes(arr, mask=mask, transform=tfm, connectivity=connectivity):
                if val in valid_classes:
                    poly = shape(geom)
                    if not poly.is_valid:
                        poly = poly.buffer(0)
                    if not poly.is_empty:
                        polygons.append(poly)

        if not polygons:
            raise ValueError("No polygons found to build AOI.")

        unified = unary_union(polygons)
        gdf = gpd.GeoDataFrame([metadata_row], geometry=[unified], crs=crs)
        gdf.to_file(gpkg_path, layer="AOI", driver="GPKG")

    def raster_data_bbox(self, in_raster, gpkg_path, metadata_row, nodata=None):
        nodata = self.NODATA_VAL if nodata is None else nodata
        with rasterio.open(in_raster) as src:
            arr = src.read(1)
            crs = src.crs
            valid = arr != nodata
            if np.issubdtype(arr.dtype, np.floating):
                valid &= ~np.isnan(arr)
            valid_mask = valid.astype(np.uint8)
            geoms = [
                shape(geom)
                for geom, _ in shapes(valid_mask, mask=valid_mask, transform=src.transform)
            ]

        if not geoms:
            raise ValueError("No valid data found to build AOI.")

        unified = unary_union(geoms)
        if self.SIMPLIFY_TOL:
            unified = unified.simplify(self.SIMPLIFY_TOL, preserve_topology=True)
            unified = make_valid(unified)

        gdf = gpd.GeoDataFrame([metadata_row], geometry=[unified], crs=crs)
        gdf.to_file(gpkg_path, layer="AOI", driver="GPKG")

    # Unified flood polygon + selected metadata -> one-feature GeoJSON
    def write_extent_geojson(self, unified_geom, geom_crs, meta_dict, geojson_dest, bm_filename):
        out_path = Path(geojson_dest) / bm_filename.replace("BM.tif", "extent.geojson")
        return write_extent_geojson(
            unified_geom,
            geom_crs,
            meta_dict,
            self.GEOJSON_FIELDS,
            out_path,
            simplify_tol=self.GEOJSON_SIMPLIFY_TOL,
            log=self.log,
        )

    def process(
        self,
        input_path,
        base_dest,
        start_date=None,
        end_date=None,
        intermediate_folder=None,
        geojson_dest=None,
    ):
        """
        Process a folder of rasters or a single .tif.

        start_date / end_date: the flood window as "YYYYMMDD" (e.g.
        "20160928" / "20161009"). Pass both to use them for every file;
        leave both None to read each file's window from its name
        (..._YYYYMMDD_YYYYMMDD_...). If it cannot be read, a single file
        stops with an error; a folder skips that file.
        intermediate_folder: where temporary rasters are written; defaults to a
        fresh temp dir. Removed (best-effort) when processing finishes.
        geojson_dest: folder for the per-map extent GeoJSON (unified flood
        polygon + key metadata); must differ from base_dest. None skips it.
        """
        BASE_DEST = Path(base_dest)
        if geojson_dest is not None and Path(geojson_dest).resolve() == BASE_DEST.resolve():
            raise ValueError("geojson_dest must be a different folder from base_dest.")
        if (start_date is None) != (end_date is None):
            raise ValueError("Pass both start_date and end_date, or neither to read them per file.")

        tif_list = list_input_tifs(input_path)
        if not tif_list:
            self.log(f"No .tif files found under {input_path}")
            return

        is_single = len(tif_list) == 1
        self.log(f"Found {len(tif_list)} .tif file(s) under {input_path}")

        # Set up the intermediate work dir (temp if not provided).
        if intermediate_folder is None:
            self._work_dir = Path(tempfile.mkdtemp(prefix="fimbench_"))
        else:
            self._work_dir = Path(intermediate_folder)
        self._work_dir.mkdir(parents=True, exist_ok=True)

        try:
            for tif in tif_list:
                # Resolve the flood window (outside the per-file try so a single
                # un-datable file stops, while a folder run skips and continues).
                if start_date is not None:
                    file_start, file_end = start_date, end_date
                else:
                    file_start, file_end = self.infer_flood_window_from_path(tif)
                if file_start is None:
                    msg = (
                        f"Could not read the flood window from {tif.name}; expected "
                        "..._YYYYMMDD_YYYYMMDD_... or pass start_date/end_date."
                    )
                    if is_single:
                        raise ValueError(msg)
                    self.log(f"[notice] {msg} Skipping.")
                    continue

                try:
                    self.log(f"Processing: {tif} (flood window {file_start}-{file_end})")

                    # Pre-processing steps
                    nodata_tif = self.process_tif(tif)
                    out_5070 = self.reproject_raster(nodata_tif, 5070)
                    res_5070, final_crs = self.resample_raster(out_5070, target_res=10.0)
                    out_4326 = self.reproject_raster(nodata_tif, 4326)
                    res_m_float, res_str, res_txt = self.get_resolution_string(out_5070)

                    # Centroid (4326) + DMS
                    x, y, raster_crs, b = self.get_raster_centroid(out_4326)
                    lon_str = self.decimal_to_dms_str(x, is_lat=False)
                    lat_str = self.decimal_to_dms_str(y, is_lat=True)
                    dms_code = lon_str + lat_str

                    # Metadata and geometry generation (called first to get the merged geometry)
                    temp_filename = (
                        f"{self.SENSOR_CODE}_{res_str}_{file_start}_{file_end}_{dms_code}_BM.tif"
                    )

                    meta_dict, unified_geom = self.build_metadata_dict(
                        out_4326,
                        res_5070,
                        temp_filename,
                        res_m_float,
                        res_txt,
                        dms_code,
                        x,
                        y,
                        "",
                        [],
                        [],
                        file_start,
                        file_end,
                    )

                    # HUC overlay via the ArcGIS REST service (geometry is in final_crs)
                    huc8_list, name_list, state_list = get_intersected_huc8(unified_geom, final_crs)

                    # Update the metadata dictionary with the actual HUC results
                    meta_dict["HUC8"] = huc8_list
                    meta_dict["River Basin Name"] = name_list
                    meta_dict["State"] = f"{state_list}, USA" if state_list else "USA"
                    meta_dict["Description"] = (
                        f"The Flood Inundation Map (FIM) was generated using USGS High Water Marks acquired "
                        f"for the {file_start}-{file_end} flood with a spatial resolution of {res_txt}. "
                        f"The corresponding HUC IDs are {huc8_list}"
                    )

                    # File management
                    folder_name = f"{self.SENSOR_CODE}_{file_start}_{file_end}_{dms_code}"
                    destination_folder = BASE_DEST / folder_name
                    destination_folder.mkdir(parents=True, exist_ok=True)

                    new_filename = (
                        f"{self.SENSOR_CODE}_{res_str}_{file_start}_{file_end}_{dms_code}_BM.tif"
                    )
                    new_path = destination_folder / new_filename

                    shutil.copy2(out_4326, new_path)
                    self.log(f"Saved GeoTIFF: {new_path.name}")

                    # Save JSON
                    metadata_filename = new_path.name.replace("BM.tif", "metadata.json")
                    metadata_path = new_path.parent / metadata_filename
                    with open(metadata_path, "w") as jf:
                        json.dump(meta_dict, jf, indent=4)
                    self.log(f"Wrote metadata: {metadata_path.name}")

                    # AOI.gpkg
                    gpkg_path = new_path.parent / new_path.name.replace("BM.tif", "AOI.gpkg")
                    self.raster_data_bbox(out_4326, gpkg_path, meta_dict, nodata=self.NODATA_VAL)
                    self.log(f"Wrote AOI: {gpkg_path.name}")

                    # Extent GeoJSON (separate destination)
                    if geojson_dest is not None:
                        self.write_extent_geojson(
                            unified_geom, final_crs, meta_dict, geojson_dest, new_filename
                        )

                except Exception as e:
                    self.log(f"Error processing {tif.name}: {e}")
                    import traceback

                    traceback.print_exc()

            self.log("Done.")
        finally:
            # Remove the intermediate work dir so the input folder stays clean
            # (best-effort; ignore if it cannot be removed).
            shutil.rmtree(self._work_dir, ignore_errors=True)
            self._work_dir = None


def process(
    input_path,
    base_dest,
    start_date=None,
    end_date=None,
    intermediate_folder=None,
    geojson_dest=None,
    **overrides,
):
    """
    Convenience wrapper: run HwmProcessor over a folder or a single .tif.

    start_date / end_date: the flood window as "YYYYMMDD" (e.g. "20160928" /
    "20161009"); leave both None to read each file's window from its name.
    intermediate_folder: temp raster dir (defaults to a temp dir, removed after).
    geojson_dest: folder for the extent GeoJSON (None skips it).
    Config overrides (e.g. source=...) are passed to the processor.
    """
    HwmProcessor(**overrides).process(
        input_path,
        base_dest,
        start_date,
        end_date,
        intermediate_folder=intermediate_folder,
        geojson_dest=geojson_dest,
    )
