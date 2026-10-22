"""
Authors: Supath Dhital (sdhital@ua.edu), Dipsikha Devi (ddevi@ua.edu)
Updated: June 2026

Tier 1 flood map processing (Aerial Imagery).

Standardizes raw Aerial Imagery flood rasters into the FIM database format:
renamed GeoTIFF, metadata.json, and AOI.gpkg per flood map.
"""

import re
import copy
import json
import shutil
import tempfile
from pathlib import Path
from typing import ClassVar
import numpy as np
import rasterio
from rasterio.features import shapes
from rasterio.warp import calculate_default_transform, reproject, Resampling
from shapely.geometry import shape
from shapely.ops import unary_union
import geopandas as gpd

from .._log import log as _log
from .utils import (
    data_footprint,
    flood_polygon,
    get_intersected_huc8,
    list_input_tifs,
    write_extent_geojson,
)
from .recurrence_interval import (
    get_event_recurrence_interval_nwm,
    get_event_recurrence_interval_usgs,
    ri_years,
)


class Tier1Processor:
    """Standardize Tier 1 (Aerial Imagery) flood maps for the FIM database."""

    SENSOR_CODE = "AI"
    # Fixed metadata.json values, keyed by their field names (as in
    # GEOJSON_FIELDS). full_form=, data_provider= and quality= override them.
    METADATA_FIELDS: ClassVar[dict] = {
        "Full sensor name": "Aerial Imagery",
        "Data provider": "Dr. Dinuke Munasinghe, The University of Alabama. Email: dsmunasinghe@ua.edu",
        "Recurrence Interval Methodology Documentation": "https://github.com/sdmlua/fimbench/blob/main/docs/recurrence-interval_documentation.md",
        "Keywords": ["flood inundation", "benchmark dataset", "geospatial data", "aerial imagery"],
        "Access Rights": "Public",
        "Product Category": "Tier 1",
        "Metadata Documentation": "https://github.com/sdmlua/fimbench/blob/main/docs/metadata_documentation.md",
        "References": [
            "1. Cohen, S., Baruah, A., Nikrou, P., Tian, D., & Liu, H. (2025). Toward robust evaluations of flood inundation predictions using remote sensing derived benchmark maps. Water Resources Research, 61(8), e2024WR039574. https://doi.org/10.1029/2024WR039574",
            "2. Devi, D., Dhital, S., Munasinghe, D., Cohen, S., Baruah, A., Chen, Y., ... & Pruitt, C. (2025). A Framework for the Evaluation of Flood Inundation Predictions Over Extensive Benchmark Databases. https://doi.org/10.1016/j.envsoft.2025.106786",
        ],
    }

    # METADATA_FIELDS written at a fixed position near the top of
    # metadata.json; the rest follow the recurrence-interval fields.
    _PLACED_FIELDS = ("Full sensor name", "Data provider")

    NODATA_VAL = -9999
    BLOCK_SIZE = 512
    SIMPLIFY_TOL = 0.0001  # ~11 m at the equator

    # Extent GeoJSON is simplified by this many metres (EPSG:5070) to stay
    # light; the metadata/AOI geometry is not affected. None disables it.
    GEOJSON_SIMPLIFY_TOL = 10.0

    # metadata.json fields copied onto the extent GeoJSON feature.
    GEOJSON_FIELDS = (
        "File Name",
        "Full sensor name",
        "HUC8 ID(s)",
        "Spatial resolution (m)",
        "Bounding box",
        "State",
        "River Basin Name(s)",
        "Image acquisition date",
        "Event Recurrence Interval (USGS) [yr]",
        "Recurrence Interval Discharge (USGS) [cms]",
        "Reference USGS Streamgauge ID",
        "Event Recurrence Interval (NWM) [yr]",
        "Recurrence Interval Discharge (NWM) [cms]",
        "Product Category",
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
        data_provider=None,
        source=None,
        quality=None,
        nodata_val=None,
        block_size=None,
        simplify_tol=None,
        gages_shp=None,
        flow_shp=None,
        recurr_folder=None,
        downstream_gage_max_km=None,
    ):
        # Class attributes are the defaults; each processor gets its own copy
        # of METADATA_FIELDS so overrides never change the class defaults.
        self.METADATA_FIELDS = copy.deepcopy(self.METADATA_FIELDS)
        if sensor_code is not None:
            self.SENSOR_CODE = sensor_code
        if full_form is not None:
            self.METADATA_FIELDS["Full sensor name"] = full_form
        # source= is the older name for data_provider=
        if data_provider is not None or source is not None:
            self.METADATA_FIELDS["Data provider"] = (
                data_provider if data_provider is not None else source
            )
        if quality is not None:
            self.METADATA_FIELDS["Product Category"] = quality
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

    def infer_flood_date_from_path(self, path: Path, default=None):
        """
        Extracts flood date/time from filename.
        Returns:
          - "YYYYMMDD"
          - "YYYYMMDDT<time>"  (time can be any number of digits)
        """
        s = path.as_posix()

        # YYYYMMDD with optional T/_ separator and any length time digits
        m_iso = re.search(
            r"(?P<date>"
            r"20\d{2}"  # Year
            r"(?:0[1-9]|1[0-2])"  # Month
            r"(?:0[1-9]|[12]\d|3[01])"  # Day
            r")"
            r"(?:[Tt_](?P<time>\d+)(?![\dNSEWnsew]))?",  # Time: any digits, not a DMS code
            s,
        )

        if m_iso:
            date_part = m_iso.group("date")
            time_part = m_iso.group("time")
            return f"{date_part}T{time_part}" if time_part else date_part

        # YYYY-MM-DD or YYYY_MM_DD
        m_sep = re.search(
            r"(?P<y>20\d{2})[-_]" r"(?P<m>0[1-9]|1[0-2])[-_]" r"(?P<d>0[1-9]|[12]\d|3[01])", s
        )
        if m_sep:
            return f"{m_sep.group('y')}{m_sep.group('m')}{m_sep.group('d')}"

        # DD-MM-YYYY or MM-DD-YYYY
        m_rev = re.search(
            r"(?P<d>0[1-9]|[12]\d|3[01])[-_]" r"(?P<m>0[1-9]|1[0-2])[-_]" r"(?P<y>20\d{2})", s
        )
        if m_rev:
            return f"{m_rev.group('y')}{m_rev.group('m')}{m_rev.group('d')}"

        return default

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
        flood_date,
    ):
        fixed = self.METADATA_FIELDS
        with rasterio.open(src_4326) as src:
            meta = {
                "File Name": new_filename,
                "Full sensor name": fixed["Full sensor name"],
                "HUC8 ID(s)": huc8_list,
                "Format": src.driver,
                "Nodata_Value": src.nodata,
                "Spatial resolution (m)": res_m_float,
                "Datatype": src.dtypes[0],
                "Compression type": src.tags(ns="IMAGE_STRUCTURE"),
                "Bounding box": {
                    "xmin": src.bounds.left,
                    "ymin": src.bounds.bottom,
                    "xmax": src.bounds.right,
                    "ymax": src.bounds.top,
                },
                "Centroid DMS code": dms_code,
                "Coordinate reference system": src.crs.to_string() if src.crs else "Unknown",
                "Number of bands": src.count,
                "Rows": src.height,
                "Columns": src.width,
                "State": state_list or None,
                "Description": "",
                "River Basin Name(s)": name_list,
                "Data provider": fixed["Data provider"],
                "Map centroid coordinates (longitude, latitude)": (x, y),
                "Image acquisition date": flood_date,
                "Event Recurrence Interval (USGS) [yr]": None,
                "Recurrence Interval Discharge (USGS) [cms]": None,
                "Reference USGS Streamgauge ID": None,
                "Event Recurrence Interval (NWM) [yr]": None,
                "Recurrence Interval Discharge (NWM) [cms]": None,
                # The remaining METADATA_FIELDS, under their own names and in
                # table order, so renaming a key there needs no change here.
                **{
                    k: list(v) if isinstance(v, (list, tuple)) else v
                    for k, v in fixed.items()
                    if k not in self._PLACED_FIELDS
                },
            }

        # Flood polygon from the 5070 raster: every pixel > 0, dissolved
        with rasterio.open(src_5070) as src_geom:
            geom_crs = src_geom.crs
            unified_geom = flood_polygon(
                src_geom.read(1, masked=True), src_geom.transform, self.SIMPLIFY_TOL
            )
        if unified_geom is None:
            self.log("No flooded areas found in 5070 raster.")
            return meta, None

        # Event Recurrence Interval, computed against the unified flood
        # polygon: USGS (observed gage discharge vs. published FFA) and NWM
        # (modeled outlet discharge vs. NWM v2.1 RI tables), plus the
        # discharge (cms) that produced each. USGS RI is "No USGS Gauge" when
        # no gage was found and "No usable USGS Gauge" when none could be
        # used; a field is None (logged, not raised) when
        # its reference data isn't configured, or the lookup couldn't be
        # attempted/resolved at all.
        ri_usgs, discharge_usgs, _scenario_usgs, gauge_id_usgs = get_event_recurrence_interval_usgs(
            unified_geom,
            geom_crs,
            flood_date,
            gages_shp=self.GAGES_SHP,
            flow_shp=self.FLOW_SHP,
            downstream_max_km=self.DOWNSTREAM_GAGE_MAX_KM,
            log=self.log,
        )
        ri_nwm, discharge_nwm, _scenario_nwm = get_event_recurrence_interval_nwm(
            unified_geom,
            geom_crs,
            flood_date,
            flow_shp=self.FLOW_SHP,
            recurr_folder=self.RECURR_FOLDER,
            log=self.log,
        )
        # Years as a number ("2 yr" -> 2); "> 100" / "< 2" and no-estimate text
        # (e.g. "No usable USGS Gauge") stay as text.
        meta["Event Recurrence Interval (USGS) [yr]"] = ri_years(ri_usgs)
        meta["Recurrence Interval Discharge (USGS) [cms]"] = (
            round(discharge_usgs, 2) if discharge_usgs is not None else None
        )
        meta["Reference USGS Streamgauge ID"] = gauge_id_usgs
        meta["Event Recurrence Interval (NWM) [yr]"] = ri_years(ri_nwm)
        meta["Recurrence Interval Discharge (NWM) [cms]"] = (
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
        # Valid-data footprint, read in chunks on a coarse grid: the 0.3 m
        # aerial rasters are far too large to read and polygonize at full
        # resolution (see utils.data_footprint).
        nodata = self.NODATA_VAL if nodata is None else nodata
        unified, crs = data_footprint(in_raster, nodata, self.SIMPLIFY_TOL)
        if unified is None:
            raise ValueError("No valid data found to build AOI.")

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
        self, input_path, base_dest, flood_date=None, intermediate_folder=None, geojson_dest=None
    ):
        """
        Process a folder of rasters or a single .tif.

        flood_date: pass it (e.g. for a single raster) to use it for every file;
        leave None to infer it per file from the filename. If it cannot be
        inferred, a single file stops with a notice; a folder skips that file.
        intermediate_folder: where temporary rasters are written; defaults to a
        fresh temp dir. Removed (best-effort) when processing finishes.
        geojson_dest: folder for the per-map extent GeoJSON (unified flood
        polygon + key metadata); must differ from base_dest. None skips it.
        """
        BASE_DEST = Path(base_dest)
        if geojson_dest is not None and Path(geojson_dest).resolve() == BASE_DEST.resolve():
            raise ValueError("geojson_dest must be a different folder from base_dest.")

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
                # Resolve the flood date (outside the per-file try so a single
                # un-datable file stops, while a folder run skips and continues).
                file_flood_date = (
                    flood_date if flood_date is not None else self.infer_flood_date_from_path(tif)
                )
                if file_flood_date is None:
                    msg = f"Could not infer flood_date from {tif.name}; pass flood_date=YYYYMMDD."
                    if is_single:
                        raise ValueError(msg)
                    self.log(f"[notice] {msg} Skipping.")
                    continue

                try:
                    self.log(f"Processing: {tif}")

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
                        f"{self.SENSOR_CODE}_{res_str}_{file_flood_date}_{dms_code}_BM.tif"
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
                        file_flood_date,
                    )

                    # HUC overlay via the ArcGIS REST service (geometry is in final_crs)
                    huc8_list, name_list, state_list = get_intersected_huc8(unified_geom, final_crs)

                    # Update the metadata dictionary with the actual HUC results
                    meta_dict["HUC8 ID(s)"] = huc8_list
                    meta_dict["River Basin Name(s)"] = name_list
                    meta_dict["State"] = state_list or None
                    meta_dict["Description"] = (
                        f"The Flood Inundation Map (FIM) was generated using NOAA Emergency Response Imagery "
                        f"acquired on {file_flood_date} with a spatial resolution of {res_txt}. "
                        f"The corresponding HUC IDs are {huc8_list}"
                    )

                    # File management
                    folder_name = f"{self.SENSOR_CODE}_{file_flood_date}_{dms_code}"
                    destination_folder = BASE_DEST / folder_name
                    destination_folder.mkdir(parents=True, exist_ok=True)

                    new_filename = (
                        f"{self.SENSOR_CODE}_{res_str}_{file_flood_date}_{dms_code}_BM.tif"
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
                    self.log(f"Skipped {tif.name} due to error: {e}")
                    import traceback

                    traceback.print_exc()  # Useful for debugging the exact line of failure

            self.log("Done.")
        finally:
            # Remove the intermediate work dir so the input folder stays clean
            # (best-effort; ignore if it cannot be removed).
            shutil.rmtree(self._work_dir, ignore_errors=True)
            self._work_dir = None


def process(
    input_path, base_dest, flood_date=None, intermediate_folder=None, geojson_dest=None, **overrides
):
    """
    Convenience wrapper: run Tier1Processor over a folder or a single .tif.

    flood_date: pass it for a single raster, else inferred per file.
    intermediate_folder: temp raster dir (defaults to a temp dir, removed after).
    geojson_dest: folder for the extent GeoJSON (None skips it).
    Config overrides (e.g. data_provider=...) are passed to the processor.
    """
    Tier1Processor(**overrides).process(
        input_path,
        base_dest,
        flood_date=flood_date,
        intermediate_folder=intermediate_folder,
        geojson_dest=geojson_dest,
    )
