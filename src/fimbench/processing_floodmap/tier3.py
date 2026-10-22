"""
Authors: Supath Dhital (sdhital@ua.edu), Dipsikha Devi (ddevi@ua.edu)
Updated: June 2026

Tier 3 flood map processing (Sentinel 1A).

Standardizes raw Sentinel 1A flood rasters into the FIM database format:
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
from shapely.validation import make_valid
import geopandas as gpd

from .._log import log as _log
from .utils import (
    get_intersected_huc8,
    list_input_tifs,
    flood_polygon,
    split_flood_datetime,
    write_extent_geojson,
)
from .recurrence_interval import (
    get_event_recurrence_interval_nwm,
    get_event_recurrence_interval_usgs,
    ri_years,
)


class Tier3Processor:
    """Standardize Tier 3 (Sentinel 1A) flood maps for the FIM database."""

    SENSOR_CODE = "S1A"
    # Fixed metadata.json values, keyed by their field names (as in
    # GEOJSON_FIELDS). full_form=, data_provider= and quality= override them.
    METADATA_FIELDS: ClassVar[dict] = {
        "Full sensor name": "Sentinel-1A",
        # Nominal Sentinel-1 resolution, used everywhere instead of the
        # measured pixel size: file names (S1A_10_0m_...), metadata, description.
        "Spatial resolution (m)": 10,
        "Data provider": "Dr. Dan Tian, The University of Alabama. Email ID: dtian1@ua.edu",
        "Recurrence Interval Methodology Documentation": "https://github.com/sdmlua/fimbench/blob/main/docs/recurrence-interval_documentation.md",
        "Keywords": ["flood inundation", "benchmark dataset", "geospatial data", "sentinel-1"],
        "Access Rights": "Public",
        "Product Category": "Tier 3",
        "Metadata Documentation": "https://github.com/sdmlua/fimbench/blob/main/docs/metadata_documentation.md",
        "References": [
            "1. Tian, D., Liu, H., Wang, L., Cohen, S., & Mandal, T. (2026). RS-FloodXDepth: Enhancing remote sensing-derived flood extent and estimating flood depth using a hydrologically guided region-growing method and high-resolution DEMs. Water Resources Research, 62, e2025WR042384. https://doi.org/10.1029/2025WR042384",
            "2. Cohen, S., Baruah, A., Nikrou, P., Tian, D., & Liu, H. (2025). Toward robust evaluations of flood inundation predictions using remote sensing derived benchmark maps. Water Resources Research, 61(8), e2024WR039574. https://doi.org/10.1029/2024WR039574",
        ],
    }

    # METADATA_FIELDS written at a fixed position near the top of
    # metadata.json; the rest follow the recurrence-interval fields.
    _PLACED_FIELDS = ("Full sensor name", "Spatial resolution (m)", "Data provider")

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
        "Image acquisition datetime (UTC)",
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

    # Extract the resolution from the projected file
    def get_resolution_string(self, geo_tif_5070=None):
        # Tier 3 always uses the nominal Sentinel-1 resolution from
        # METADATA_FIELDS (10 m), not the measured pixel size, so file names,
        # metadata and description agree.
        rounded_res = float(self.METADATA_FIELDS["Spatial resolution (m)"])
        res_str = f"{str(rounded_res).replace('.', '_')}m"  # e.g., 10.0 -> "10_0m"
        res_txt = f"{rounded_res:.2f} m"
        return rounded_res, res_str, res_txt

    def get_raster_centroid(self, out_4326):
        with rasterio.open(out_4326) as src:
            b = src.bounds
            centroid_x = (b.left + b.right) / 2
            centroid_y = (b.top + b.bottom) / 2
            crs = src.crs
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
        Returns: "YYYYMMDD" or "YYYYMMDDTYYYYMM" (String)
        """
        s = path.as_posix()

        # YYYYMMDD with optional T/separator and time. Captures: 20170502, 20170502T161000, 20170502_161000
        m_iso = re.search(
            r"(?P<date>20\d{2}(?:0[1-9]|1[0-2])(?:0[1-9]|[12]\d|3[01]))"  # YYYYMMDD
            r"(?:[Tt_](?P<time>\d{6})(?![\dNSEWnsew]))?",  # Optional Time, not a DMS code
            s,
        )
        if m_iso:
            date_part = m_iso.group("date")
            time_part = m_iso.group("time")
            return f"{date_part}T{time_part}" if time_part else date_part

        # YYYY-MM-DD or YYYY_MM_DD (hyphenated/underscore)
        m_sep = re.search(
            r"(?P<y>20\d{2})[-_](?P<m>0[1-9]|1[0-2])[-_](?P<d>0[1-9]|[12]\d|3[01])", s
        )
        if m_sep:
            return f"{m_sep.group('y')}{m_sep.group('m')}{m_sep.group('d')}"

        # DD-MM-YYYY or MM-DD-YYYY (defaults to DD-MM-YYYY logic)
        m_rev = re.search(
            r"(?P<d>0[1-9]|[12]\d|3[01])[-_](?P<m>0[1-9]|1[0-2])[-_](?P<y>20\d{2})", s
        )
        if m_rev:
            return f"{m_rev.group('y')}{m_rev.group('m')}{m_rev.group('d')}"

        return default

    @staticmethod
    def describe(flood_date, res_m, huc8_list):
        """Description text, e.g. "Flood Inundation Map (FIM) generated using
        Sentinel-1 Imagery acquired on 2017-05-04 at 23:55:30 UTC, ...".
        The time is left out when flood_date has none."""
        day, time_utc = split_flood_datetime(flood_date)
        when = f"{day[:4]}-{day[4:6]}-{day[6:8]}" if len(day) >= 8 else day
        if time_utc:
            when += f" at {time_utc} UTC"
        hucs = [str(h) for h in huc8_list or []]
        if len(hucs) > 1:
            covers = (
                f" and covers portions of HUC8 watersheds {', '.join(hucs[:-1])} and {hucs[-1]}"
            )
        elif hucs:
            covers = f" and covers portions of HUC8 watershed {hucs[0]}"
        else:
            covers = ""
        return (
            f"Flood Inundation Map (FIM) generated using Sentinel-1 Imagery acquired on {when}, "
            "with a gap-filling algorithm applied to reconstruct inundation in unmapped areas. "
            f"The map has a nominal spatial resolution of {res_m:g} m{covers}."
        )

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
        # flood_date may carry the acquisition time, e.g. 20170504T235530
        flood_day, acquisition_time_utc = split_flood_datetime(flood_date)
        fixed = self.METADATA_FIELDS
        with rasterio.open(src_4326) as src:
            meta = {
                "File Name": new_filename,
                "Full sensor name": fixed["Full sensor name"],
                "HUC8 ID(s)": huc8_list,
                "Format": src.driver,
                "Nodata_Value": src.nodata,
                "Spatial resolution (m)": fixed["Spatial resolution (m)"],
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
                "Image acquisition datetime (UTC)": flood_date,
                "Image acquisition date": flood_day,
                "Image acquisition time (UTC)": acquisition_time_utc,
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
        # discharge (cms) that produced each. flood_date here typically
        # carries a UTC time (Tier 2/3), so the discharge used is the
        # observation/model timestep nearest that exact instant rather than
        # the day's max -- see recurrence_interval. USGS RI is "No USGS
        # Gauge" when no gage was found and "No usable USGS Gauge" when none
        # could be used; a field is None (logged, not
        # raised) when its reference data isn't configured, or the lookup
        # couldn't be attempted/resolved at all.
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
        self,
        in_raster,
        gpkg_path,
        metadata_row,
        valid_classes=(0, 1, 2),
        nodata=None,
        connectivity=8,
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
                    out_4326 = self.reproject_raster(nodata_tif, 4326)
                    res_m_float, res_str, res_txt = self.get_resolution_string(out_5070)

                    # Centroid (4326) + DMS
                    x, y, raster_crs, b = self.get_raster_centroid(out_4326)
                    lon_str = self.decimal_to_dms_str(x, is_lat=False)
                    lat_str = self.decimal_to_dms_str(y, is_lat=True)
                    dms_code = lon_str + lat_str

                    # Metadata and geometry generation
                    temp_filename = (
                        f"{self.SENSOR_CODE}_{res_str}_{file_flood_date}_{dms_code}_BM.tif"
                    )

                    meta_dict, unified_geom = self.build_metadata_dict(
                        out_4326,
                        out_5070,
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

                    # HUC overlay via the ArcGIS REST service (unified_geom is in EPSG:5070)
                    huc8_list, name_list, state_list = get_intersected_huc8(unified_geom, 5070)

                    # Update the metadata dictionary with the actual HUC results
                    meta_dict["HUC8 ID(s)"] = huc8_list
                    meta_dict["River Basin Name(s)"] = name_list
                    meta_dict["State"] = state_list or None
                    meta_dict["Description"] = self.describe(
                        file_flood_date, self.METADATA_FIELDS["Spatial resolution (m)"], huc8_list
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
                            unified_geom, 5070, meta_dict, geojson_dest, new_filename
                        )

                except Exception as e:
                    self.log(f"Skipped {tif.name} due to error: {e}")
                    import traceback

                    traceback.print_exc()

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
    Convenience wrapper: run Tier3Processor over a folder or a single .tif.

    flood_date: pass it for a single raster, else inferred per file.
    intermediate_folder: temp raster dir (defaults to a temp dir, removed after).
    geojson_dest: folder for the extent GeoJSON (None skips it).
    Config overrides (e.g. source=...) are passed to the processor.
    """
    Tier3Processor(**overrides).process(
        input_path,
        base_dest,
        flood_date=flood_date,
        intermediate_folder=intermediate_folder,
        geojson_dest=geojson_dest,
    )
