"""
Test: standardize raw flood maps (fimbench.processing_floodmap).

Each writes the renamed GeoTIFF + metadata.json + AOI.gpkg into out_dir, and
(if geojson_dest is set) an extent GeoJSON into geojson_dest.

input_path can be a FOLDER (every .tif under it) or a single .tif.
  - tier1/2/3: flood_date is inferred per file from the filename when processing
    a folder; pass flood_date=... for a single raster.
  - fema_ble: pass event (return period).
  - hwm: pass start_date and end_date.
Config defaults (source, sensor_code, ...) can be overridden as keyword args.
Run with pytest; comment out any test function you do not want to run.
"""

import fimbench

input_folder = "path/to/input_folder"
# input_tif = "path/to/input.tif"
out_dir = "path/to/output_folder"
# Extent GeoJSON (flood polygon + key metadata), one file per map. Must be a
# different folder from out_dir. None -> no GeoJSON is written.
geojson_dest = "path/to/geojson_folder"
event = "500"
start_date = "20100329"  # HWM flood window start (YYYYMMDD) -- edit to your event
end_date = "20100401"  # HWM flood window end (YYYYMMDD) -- edit to your event
# flood_date = "20160928"
intermediate_folder = None  # None -> a temp dir is used and removed afterward

# Config overrides. None -> keep the module default.
# source = "Dr. Dinuke Munasinghe,The University of Alabama, Email id: dsmunasinghe@ua.edu"  # e.g. "Dr. X, ..." Creator name and email
nodata_val = None  # default -9999
block_size = None  # default 512 (GeoTIFF tile size)
simplify_tol = None  # default 0.0001 deg (~11 m) geometry simplify tolerance

# Tier 1/2/3 and HWM: reference data for the "Event Recurrence Interval (USGS)" and
# "Event Recurrence Interval (NWM)" metadata fields. None -> that lookup is
# skipped (logged, not an error) and the metadata field is written as None.
gages_shp = "path/to/USGSgages.shp"  # USGS gage points: SOURCE_FEA, FLComID
flow_shp = "path/to/flowline_network.shp"  # NHDPlus flowlines: COMID, FromNode, ToNode
recurr_folder = "path/to/NWM_recurrence_folder"  # folder of nwm21_17C_recurr_{yr}_0_cms.csv
downstream_gage_max_km = None  # default 100.0 km downstream USGS gage search radius

# def test_processing_tier1_folder():
#     # Folder: flood_date inferred per file
# fimbench.Tier1Processor(nodata_val=nodata_val, block_size=block_size, simplify_tol=simplify_tol,gages_shp=gages_shp,flow_shp=flow_shp,recurr_folder=recurr_folder,downstream_gage_max_km=downstream_gage_max_km).process(
# input_folder, out_dir, intermediate_folder=intermediate_folder, geojson_dest=geojson_dest
# )


# def test_processing_tier1_single():
# Single raster: pass the flood_date; override source if you want
# fimbench.Tier1Processor(
#      source=source,
#     nodata_val=nodata_val,
#     block_size=block_size,
#      simplify_tol=simplify_tol,
#     gages_shp=gages_shp,
#     flow_shp=flow_shp,
#     recurr_folder=recurr_folder,
#     downstream_gage_max_km=downstream_gage_max_km,
# ).process(input_tif, out_dir, flood_date=flood_date, intermediate_folder=intermediate_folder, geojson_dest=geojson_dest)


# def test_processing_tier2_folder():
# fimbench.Tier2Processor(nodata_val=nodata_val, block_size=block_size, simplify_tol=simplify_tol,gages_shp=gages_shp,flow_shp=flow_shp,recurr_folder=recurr_folder,downstream_gage_max_km=downstream_gage_max_km).process(
#  input_folder, out_dir, intermediate_folder=intermediate_folder, geojson_dest=geojson_dest
# )


# def test_processing_tier2_single():
#     fimbench.Tier2Processor(source=source, nodata_val=nodata_val, block_size=block_size, simplify_tol=simplify_tol).process(
#         input_tif, out_dir, flood_date=flood_date, intermediate_folder=intermediate_folder, geojson_dest=geojson_dest
#     )


# def test_processing_tier3_folder():
# fimbench.Tier3Processor(nodata_val=nodata_val, block_size=block_size, simplify_tol=simplify_tol,gages_shp=gages_shp,flow_shp=flow_shp,recurr_folder=recurr_folder,downstream_gage_max_km=downstream_gage_max_km).process(
# input_folder, out_dir, intermediate_folder=intermediate_folder, geojson_dest=geojson_dest
# )


# def test_processing_tier3_single():
# fimbench.Tier3Processor(
#    source=source, nodata_val=nodata_val, block_size=block_size, simplify_tol=simplify_tol
# ).process(input_tif, out_dir, flood_date=flood_date, intermediate_folder=intermediate_folder, geojson_dest=geojson_dest)


def test_processing_fema_ble():
    fimbench.FemaBleProcessor(
        nodata_val=nodata_val, block_size=block_size, simplify_tol=simplify_tol
    ).process(
        input_folder,
        out_dir,
        event,
        intermediate_folder=intermediate_folder,
        geojson_dest=geojson_dest,
    )


# def test_processing_hwm_folder():
# Folder: each file's start/end date is read from its name
# (..._YYYYMMDD_YYYYMMDD_...); a file without one is skipped
# fimbench.HwmProcessor(
# nodata_val=nodata_val,
# block_size=block_size,
# simplify_tol=simplify_tol,
# gages_shp=gages_shp,
# flow_shp=flow_shp,
# recurr_folder=recurr_folder,
# downstream_gage_max_km=downstream_gage_max_km,
# ).process(
# input_folder,
# out_dir,
# intermediate_folder=intermediate_folder,
# geojson_dest=geojson_dest,
# )


# def test_processing_hwm_single():
# Single raster: pass the flood window (start_date/end_date) for this map
# fimbench.HwmProcessor(
# nodata_val=nodata_val,
# block_size=block_size,
# simplify_tol=simplify_tol,
# gages_shp=gages_shp,
# flow_shp=flow_shp,
# recurr_folder=recurr_folder,
# downstream_gage_max_km=downstream_gage_max_km,
# ).process(
# input_tif,
# out_dir,
# start_date,
# end_date,
# intermediate_folder=intermediate_folder,
# geojson_dest=geojson_dest,
# )
