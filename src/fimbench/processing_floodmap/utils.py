"""
Authors: Supath Dhital (sdhital@ua.edu), Dipsikha Devi (ddevi@ua.edu)
Updated: June 2026

Resolve the HUC8 watersheds a flood map overlaps via the public ArcGIS REST
service, so no HUC8 shapefile has to be loaded. get_intersected_huc8(geom, crs)
returns (huc8_list, name_list, state_list) ready for the metadata, matching the
tiers' previous shapefile-intersection output. get_dominant_huc8 returns the
same triple but with only the HUC8 that covers most of the flood extent, which
is what the per-watershed FEMA BLE maps need.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import requests
import geopandas as gpd
import pandas as pd
import rasterio
import shapely
from affine import Affine
from rasterio.features import shapes
from rasterio.windows import Window
from shapely.geometry import shape
from shapely.prepared import prep


def list_input_tifs(input_path):
    """
    Resolve input_path to a list of .tif files.

    Accepts a single .tif file or a folder (searched recursively). Raises
    FileNotFoundError if the path does not exist.
    """
    p = Path(input_path)
    if p.is_file():
        return [p]
    if p.is_dir():
        return sorted(f for f in p.rglob("*.tif") if f.is_file())
    raise FileNotFoundError(f"Input path not found: {p}")


# Public ArcGIS REST service for WBD HUC8 boundaries.
HUC8_REST_URL = (
    "https://services.arcgis.com/ts4gk3YgS68yLGFl/arcgis/rest/services/"
    "HUC8_Boundaries/FeatureServer/0"
)

# One response can only carry so many features, so results are paged.
_PAGE_SIZE = 1000


def _bounds_4326(geom, crs):
    """Return the (minx, miny, maxx, maxy) of ``geom`` in EPSG:4326."""
    minx, miny, maxx, maxy = gpd.GeoSeries([geom], crs=crs).to_crs(4326).total_bounds
    return float(minx), float(miny), float(maxx), float(maxy)


def query_huc8(geom, crs):
    """
    Query the REST service for HUC8 features intersecting ``geom``.

    geom: a shapely geometry. crs: its CRS (EPSG int or anything pyproj accepts).
    Returns a GeoDataFrame (EPSG:4326) with columns HUC8, NAME, STATES; empty if
    nothing intersects.

    Only the bounding box of ``geom`` is sent, so the request stays the same
    small size no matter how many vertices the flood polygons have (a fine
    resolution raster easily vectorizes into hundreds of thousands of vertices,
    which the service rejects with HTTP 413). The bbox over-selects, and
    get_intersected_huc8 narrows the result down locally against the real
    geometry.
    """
    empty = gpd.GeoDataFrame({"HUC8": [], "NAME": [], "STATES": []}, geometry=[], crs=4326)
    if geom is None or geom.is_empty:
        return empty

    minx, miny, maxx, maxy = _bounds_4326(geom, crs)
    if not all(v == v for v in (minx, miny, maxx, maxy)):  # NaN bounds
        return empty

    params = {
        "f": "geojson",
        "where": "1=1",
        "geometry": f"{minx},{miny},{maxx},{maxy}",
        "geometryType": "esriGeometryEnvelope",
        "spatialRel": "esriSpatialRelIntersects",
        "inSR": 4326,
        "outFields": "HUC8,name,states",
        "returnGeometry": "true",
        "outSR": 4326,
        "resultRecordCount": _PAGE_SIZE,
    }

    # Page until the service stops reporting more features, so an arbitrarily
    # large extent still returns every HUC8 it touches.
    frames = []
    offset = 0
    while True:
        res = requests.post(
            f"{HUC8_REST_URL}/query", data={**params, "resultOffset": offset}, timeout=60
        )
        res.raise_for_status()
        payload = json.loads(res.content)
        features = payload.get("features") or []
        if features:
            frames.append(gpd.GeoDataFrame.from_features(features, crs=4326))
        if not features or not payload.get("exceededTransferLimit"):
            break
        offset += len(features)

    if not frames:
        return empty

    gdf = gpd.GeoDataFrame(pd.concat(frames, ignore_index=True), geometry="geometry", crs=4326)

    # Match the column names the tiers expect.
    gdf = gdf.rename(columns={"name": "NAME", "states": "STATES"})
    return gdf


def get_intersected_huc8(geom, crs):
    """
    Return (huc8_list, name_list, state_list) for the HUC8s overlapping ``geom``.

    Reproduces the tiers' previous output: sorted unique HUC8 and NAME, and a
    comma-split/sorted/de-duplicated STATES string. Falls back to
    ([], [], "USA") when nothing intersects.
    """
    gdf = query_huc8(geom, crs)
    if gdf.empty:
        return [], [], "USA"

    # The REST query used the bounding box, so narrow it down here against the
    # real geometry. prep() keeps this cheap for very detailed flood polygons.
    g = gpd.GeoSeries([geom], crs=crs).to_crs(4326).iloc[0]
    hit = gdf[gdf.geometry.apply(prep(g).intersects)].copy()
    if hit.empty:
        # Sliver geometries can miss on exact predicates; retry with a ~0.1 m
        # tolerance before giving up.
        hit = gdf[gdf.geometry.apply(prep(g.buffer(1e-6)).intersects)].copy()
    if hit.empty:
        return [], [], "USA"

    return _format_huc8_lists(hit)


def _format_huc8_lists(hit):
    """Build the (huc8_list, name_list, state_list) triple the tiers write out."""
    huc8_list = sorted(hit["HUC8"].dropna().unique().tolist())
    name_list = sorted(hit["NAME"].dropna().unique().tolist())
    raw_states = hit["STATES"].dropna().unique().tolist()
    state_set = set()
    for s in raw_states:
        for part in str(s).split(","):
            state_set.add(part.strip())
    state_list = ", ".join(sorted(list(state_set)))
    return huc8_list, name_list, state_list


# FEMA BLE maps are produced per watershed, so a single HUC8 should own almost
# all of the flood extent; anything else is edge overlap from neighbours.
DOMINANT_HUC8_FRACTION = 0.70

# Equal-area CRS (CONUS Albers) so the area shares are not latitude-distorted.
_AREA_CRS = 5070


def get_dominant_huc8(geom, crs, min_fraction=DOMINANT_HUC8_FRACTION, log=None):
    """
    Return (huc8_list, name_list, state_list) for the single HUC8 that holds
    most of ``geom``'s area.

    Same output shape as get_intersected_huc8 -- huc8_list and name_list are
    lists (of one entry), state_list a comma-joined string -- so metadata keys
    and value formats stay identical. The HUC8 covering at least ``min_fraction``
    of the flood area wins; if none reaches that share the largest one is used
    and the shares are logged. Falls back to ([], [], "USA") when nothing
    intersects.
    """
    gdf = query_huc8(geom, crs)
    if gdf.empty:
        return [], [], "USA"

    # Measure the overlap in an equal-area projection.
    g = gpd.GeoSeries([geom], crs=crs).to_crs(_AREA_CRS).iloc[0]
    total = g.area
    if total <= 0:
        return [], [], "USA"

    areas = gdf.to_crs(_AREA_CRS).geometry.apply(
        lambda h: g.intersection(h).area if g.intersects(h) else 0.0
    )
    shares = areas / total
    if shares.max() <= 0:
        return [], [], "USA"

    top = shares.idxmax()
    if log is not None and shares[top] < min_fraction:
        ranked = ", ".join(
            f"{gdf.loc[i, 'HUC8']}={shares[i]:.0%}"
            for i in shares.sort_values(ascending=False).index[:4]
            if shares[i] > 0
        )
        log(
            f"No HUC8 covers {min_fraction:.0%} of the flood extent "
            f"({ranked}); using the largest, {gdf.loc[top, 'HUC8']}."
        )

    return _format_huc8_lists(gdf.loc[[top]])


def _json_default(o):
    # numpy scalars/arrays (e.g. resolution, discharge) are not JSON-native.
    if hasattr(o, "item"):
        return o.item()
    if hasattr(o, "tolist"):
        return o.tolist()
    return str(o)


def flood_polygon(arr, transform, simplify_tol=None):
    """
    Dissolved polygon of every flooded pixel (value > 0) of ``arr`` (a masked
    array; masked pixels are dry), in the raster's CRS. None if nothing floods.

    The binary flood mask is polygonized directly, so each connected flood
    patch comes out as one polygon -- no union of per-class pixel polygons
    (the slow step for large rasters) is needed.
    """
    flood = np.ma.filled(np.ma.asarray(arr) > 0, False)
    if not flood.any():
        return None

    parts = np.array(
        [shape(g) for g, _ in shapes(flood.view(np.uint8), mask=flood, transform=transform)]
    )
    return _as_multipolygon(parts, simplify_tol)


# Pixels read per chunk by data_footprint (~128 MB of int32), so a
# fine-resolution raster never has to fit in memory at once.
_FOOTPRINT_CHUNK_PIXELS = 1 << 25


def data_footprint(in_raster, nodata, simplify_tol=None):
    """
    Polygon of the valid-data (non-nodata, non-NaN) area of ``in_raster``, in
    its CRS, and that CRS. The polygon is None when the raster has no valid data.

    The raster is read in row chunks and reduced onto a coarse grid (cell =
    half of ``simplify_tol``) where a cell counts as valid if any of its pixels
    is. Only that small grid is polygonized, so a 0.3 m aerial map (tens of
    GB of pixels) is handled in bounded memory. Nodata gaps smaller than a
    cell close up, and specks smaller than ``simplify_tol`` drop out. These
    are details below the AOI's simplify tolerance anyway.
    """
    with rasterio.open(in_raster) as src:
        res = max(abs(src.res[0]), abs(src.res[1]))
        f = max(1, int(simplify_tol / 2 // res)) if simplify_tol else 1
        H, W = src.height, src.width
        h, w = -(-H // f), -(-W // f)
        coarse = np.zeros((h, w), dtype=bool)
        is_float = np.issubdtype(np.dtype(src.dtypes[0]), np.floating)

        # A multiple of f, so every chunk starts on a coarse-cell boundary.
        rows = f * max(1, _FOOTPRINT_CHUNK_PIXELS // (W * f))
        for r0 in range(0, H, rows):
            n = min(rows, H - r0)
            arr = src.read(1, window=Window(0, r0, W, n))
            valid = arr != nodata
            if is_float:
                valid &= ~np.isnan(arr)
            ch = -(-n // f)
            padded = np.zeros((ch * f, w * f), dtype=bool)
            padded[:n, :W] = valid
            coarse[r0 // f : r0 // f + ch] |= padded.reshape(ch, f, w, f).any(axis=(1, 3))

        transform = src.transform * Affine.scale(f)
        crs = src.crs

    if not coarse.any():
        return None, crs

    parts = np.array(
        [shape(g) for g, _ in shapes(coarse.view(np.uint8), mask=coarse, transform=transform)]
    )
    footprint = _as_multipolygon(parts, simplify_tol, preserve_topology=False)
    if footprint.is_empty:  # data made only of specks keeps its unsimplified footprint
        footprint = _as_multipolygon(parts)
    return footprint, crs


def _as_multipolygon(parts, simplify_tol=None, preserve_topology=True):
    # Simplify (per part, so it stays fast), repair, and gather the polygon
    # parts into one MultiPolygon. preserve_topology=False also lets specks
    # smaller than the tolerance collapse and drop out.
    if simplify_tol:
        parts = shapely.simplify(parts, simplify_tol, preserve_topology=preserve_topology)
    parts = shapely.get_parts(shapely.make_valid(parts))
    parts = parts[(shapely.get_type_id(parts) == 3) & ~shapely.is_empty(parts)]
    return shapely.multipolygons(parts)


def split_flood_datetime(flood_date):
    """
    Split a flood_date like "20170504T235530" into ("20170504", "23:55:30").

    The time is formatted HH:MM:SS (6 digits) or HH:MM (4 digits); any other
    length is returned as written. Time is None when flood_date has none.
    """
    s = str(flood_date)
    date_part, sep, time_part = s.replace("_", "T").partition("T")
    if not sep or not time_part:
        return s, None
    if len(time_part) == 6:
        return date_part, f"{time_part[:2]}:{time_part[2:4]}:{time_part[4:]}"
    if len(time_part) == 4:
        return date_part, f"{time_part[:2]}:{time_part[2:]}"
    return date_part, time_part


def _light_geometry(geom, crs, simplify_tol, precision):
    # Per-part simplify in the source (projected) CRS, dropping specks below
    # the tolerance, then EPSG:4326 with rounded coordinates -- the bulk of a
    # GeoJSON's size is coordinates.
    if simplify_tol:
        light = _as_multipolygon(shapely.get_parts(geom), simplify_tol, preserve_topology=False)
        if not light.is_empty:  # a map made only of specks keeps its full geometry
            geom = light
    g = gpd.GeoSeries([geom], crs=crs).to_crs(4326).iloc[0]
    if precision is not None:
        g = shapely.transform(g, lambda c: np.round(c, precision))
    # Simplified patches can overlap slightly; repair the whole so the
    # GeoJSON stays valid.
    return _as_multipolygon(shapely.get_parts(shapely.make_valid(g)))


def write_extent_geojson(
    geom, crs, meta, fields, out_path, simplify_tol=None, precision=6, log=None
):
    """
    Write ``geom`` (the unified flood polygon, in ``crs``) as a one-feature
    GeoJSON in EPSG:4326, with the ``fields`` of ``meta`` as its properties.

    To keep the file light, the geometry is simplified by ``simplify_tol``
    (``crs`` units, e.g. metres for EPSG:5070) -- patches smaller than that
    drop out -- and coordinates are rounded to ``precision`` decimals
    (6 ~ 0.1 m). Only the GeoJSON copy is lightened.

    Fields missing from ``meta`` are left out. Returns the written path, or
    None when there is no geometry to write.
    """
    if geom is None or geom.is_empty:
        if log is not None:
            log(f"No flood geometry; skipped extent GeoJSON {Path(out_path).name}")
        return None

    g = _light_geometry(geom, crs, simplify_tol, precision)
    feature = {
        "type": "Feature",
        "properties": {k: meta[k] for k in fields if k in meta},
        "geometry": g.__geo_interface__,
    }

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(
            {"type": "FeatureCollection", "features": [feature]},
            f,
            default=_json_default,
        )
    if log is not None:
        log(f"Wrote extent GeoJSON: {out_path}")
    return out_path
