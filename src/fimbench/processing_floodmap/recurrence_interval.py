"""
Authors: Reihaneh Zarrabi (rzarrabi@crimson.ua.edu),Dipsikha Devi (ddevi@ua.edu),Supath Dhital (sdhital@ua.edu),
Updated: September 2026

Event Recurrence Interval (RI) for a single flood map, ported from the
FIMbench_RI notebook pipeline (fim_recurrence_pipeline_T1.ipynb) to run against
one flood polygon/date instead of a batched GeoJSON of many events.

Two independent RI estimates are produced for the same flood polygon and
event date/time:
  - USGS: resolve the polygon to a USGS streamgage (in-polygon, closest to the
    outlet first; else the nearest usable gage downstream along the flow
    network), fetch its observed discharge for the event, and classify it
    against that gage's published USGS StreamStats flood-frequency
    (Bulletin 17C) quantiles.
  - NWM: resolve the polygon to its outlet NHDPlus reach (COMID), fetch NWM
    streamflow at that reach for the event, and classify it against the
    NWM v2.1 recurrence-interval tables.

flood_date accepts either form the notebook does, and the discharge timing
follows suit exactly as in fim_recurrence_pipeline_T1.ipynb:
  - "YYYYMMDD" (Tier 1, date-only): the daily maximum hourly discharge for
    that day (USGS 00060 instantaneous -> hourly max -> daily max; NWM hourly
    max for the day).
  - "YYYYMMDDTHHMMSS" (Tier 2/3, UTC date+time): the discharge/streamflow at
    the observation/model timestep nearest that exact UTC instant.
A "T..." suffix that doesn't parse as a full HHMMSS falls back to the
date-only path (logged, not raised).

Both need local reference data that is not part of the FIMbench package:
  gages_shp:     USGS gage points shapefile with SOURCE_FEA (site number) and
                 FLComID (the NHDPlus COMID it sits on).
  flow_shp:      NHDPlus flowline network shapefile with COMID, FromNode,
                 ToNode.
  recurr_folder: folder of NWM v2.1 RI CSVs (nwm21_17C_recurr_{yr}_0_cms.csv).
If a required path is not configured or the lookup can't be resolved, the
corresponding function returns None and logs why, instead of raising -- a
missing/failed RI lookup should not stop the rest of the metadata from being
written.
"""

from __future__ import annotations

import os
import re
from collections import deque
from datetime import datetime, timedelta
from functools import lru_cache

import numpy as np
import pandas as pd
import geopandas as gpd
import pyogrio
import requests
from pyproj import Transformer

CFS_TO_CMS = 0.0283168
RI_YEARS = [2, 5, 10, 25, 50, 100]
RETRO_MAX_YEAR = 2023
DOWNSTREAM_GAGE_MAX_KM = 100.0
GAP_TOL = 200  # meters; bridges multipart flood polygons for the spatial join
REGION_SLACK_M = 10_000  # meters of extra margin when reading flowlines/gages around the flood
COMID_COL = "COMID"

NWM_RETRO_ZARR = "s3://noaa-nwm-retrospective-3-0-pds/CONUS/zarr/chrtout.zarr"
NWM_OPS_BUCKET = "s3://noaa-nwm-pds"

GAGESTATS_BASE = "https://streamstats.usgs.gov/gagestatsservices"

USGS_AT_TIME_WINDOW_HOURS = 1  # +/- search window for the nearest USGS IV observation
USGS_CATALOG_RETRIES = 4  # attempts for the NWIS site-catalog request (diagnostics only)

# get_event_recurrence_interval_usgs's label when it definitively resolved
# that no usable gage exists (as opposed to being unable to look at all,
# e.g. missing config/geometry/date -- those still return None):
# NO_GAGE_LABEL when no gage was found in the polygon or downstream, and
# NO_USABLE_GAGE_LABEL when gages were found but none could be used. The
# per-gage reasons (not recording discharge, discontinued, data gap, no
# published flood-frequency statistics, ...) are logged.
NO_GAGE_LABEL = "No USGS Gauge"
NO_USABLE_GAGE_LABEL = "No usable USGS Gauge"


def _noop_log(msg):
    pass


def _parse_flood_datetime(flood_date, log=None):
    """Parse a flood_date string the way fim_recurrence_pipeline_T1.ipynb's
    parse_event_datetime does, dispatching Tier 1 (date-only) vs. Tier 2/3
    (date+UTC time) by whether a time is present.

    "YYYYMMDD"        -> (naive midnight datetime, has_time=False): use the
                          day-based (daily max) fetch path.
    "YYYYMMDDTHHMMSS" -> (UTC-aware pd.Timestamp, has_time=True): use the
                          at-time (nearest observation/timestep) fetch path.

    A "T..." suffix that isn't a full 6-digit HHMMSS falls back to the
    date-only path (logged, not raised). Returns (None, None) if even the
    date portion doesn't parse.
    """
    log = log or _noop_log
    text = str(flood_date).strip()
    if "T" in text:
        try:
            dt = datetime.strptime(text[:15], "%Y%m%dT%H%M%S")
            return pd.Timestamp(dt, tz="UTC"), True
        except ValueError:
            log(
                f"Event Recurrence Interval: could not parse time from {flood_date!r} as "
                "YYYYMMDDTHHMMSS; falling back to date-only."
            )
    try:
        return datetime.strptime(text[:8], "%Y%m%d"), False
    except ValueError:
        return None, None


def _parse_date_only(date_str, log=None):
    """Parse a plain "YYYYMMDD" date string (HWM start/end dates -- no time
    component). Returns a naive midnight datetime, or None if it doesn't
    parse (logged, not raised)."""
    log = log or _noop_log
    text = str(date_str).strip()
    try:
        return datetime.strptime(text[:8], "%Y%m%d")
    except ValueError:
        log(f"Event Recurrence Interval: could not parse date {date_str!r} as YYYYMMDD.")
        return None


# ---------------------------------------------------------------------------
# Shared: load reference layers, network topology, RI classification
# ---------------------------------------------------------------------------


def _read_region(path, columns, region, proj):
    """
    Read ``columns`` of a vector file (no CRS -> EPSG:4326), keeping only
    features within ``region`` = (geom_proj, margin_m) when given. The
    national flowline/gage files are huge, so a bbox read around the flood
    is far faster and lighter than loading everything.
    """
    bbox = None
    if region is not None:
        geom_proj, margin_m = region
        file_crs = pyogrio.read_info(path)["crs"] or "EPSG:4326"
        xmin, ymin, xmax, ymax = geom_proj.bounds
        bbox = Transformer.from_crs(proj, file_crs, always_xy=True).transform_bounds(
            xmin - margin_m, ymin - margin_m, xmax + margin_m, ymax + margin_m
        )
    gdf = gpd.read_file(path, columns=columns, bbox=bbox)
    if gdf.crs is None:
        gdf = gdf.set_crs(4326)
    return gdf.to_crs(proj)


def _load_gages(gages_shp, proj, region=None):
    gages = _read_region(gages_shp, ["SOURCE_FEA", "FLComID", "geometry"], region, proj)
    gages["FLComID"] = pd.to_numeric(gages["FLComID"], errors="coerce").astype("Int64")
    return gages


def _load_flow(flow_shp, proj, region=None):
    flow = _read_region(
        flow_shp, ["COMID", "FromNode", "ToNode", "TotDASqKm", "geometry"], region, proj
    )
    for col in ["FromNode", "ToNode", COMID_COL]:
        flow[col] = pd.to_numeric(flow[col], errors="coerce").astype("Int64")
    flow = flow.dropna(subset=["FromNode", "ToNode", COMID_COL]).copy()
    for col in ["FromNode", "ToNode", COMID_COL]:
        flow[col] = flow[col].astype("int64")
    return flow


def _clip_to_geom(gdf, geom_proj, proj):
    # Dilate by GAP_TOL so gages/reaches sitting in gaps between parts of a
    # fragmented (multipart) flood polygon are still captured. This affects
    # the join only -- the flood geometry itself is never modified.
    # "dwithin" equals intersecting a GAP_TOL buffer, without buffering a
    # (possibly million-vertex) flood polygon.
    flood = gpd.GeoDataFrame(geometry=[geom_proj], crs=proj)
    return gpd.sjoin(gdf, flood, how="inner", predicate="dwithin", distance=GAP_TOL).drop(
        columns="index_right"
    )


def _find_outlets(reaches):
    """OUTLET = a reach whose ToNode is not any other reach's FromNode."""
    from_nodes = set(reaches["FromNode"])
    return reaches[~reaches["ToNode"].isin(from_nodes)]


def _hops_from_outlet(reaches, outlet_comid, comid_col=COMID_COL):
    """BFS upstream from the outlet, counting reaches (outlet = hop 0)."""
    upstream_of, fromnode = {}, {}
    for _, r in reaches.iterrows():
        upstream_of.setdefault(r["ToNode"], []).append(r[comid_col])
        fromnode[r[comid_col]] = r["FromNode"]
    hops = {outlet_comid: 0}
    q = deque([outlet_comid])
    while q:
        cur = q.popleft()
        for up in upstream_of.get(fromnode.get(cur), []):
            if up not in hops:
                hops[up] = hops[cur] + 1
                q.append(up)
    return hops


def _pick_outlet_reach(reaches, comid_col=COMID_COL):
    """The most-downstream reach, chosen by topology only: among the
    topological sinks, the one with the most reaches upstream of it."""
    if len(reaches) == 0:
        return None
    outs = _find_outlets(reaches)
    if len(outs) <= 1:
        return outs.iloc[0] if len(outs) else reaches.iloc[0]
    best = max(outs[comid_col], key=lambda cid: len(_hops_from_outlet(reaches, cid, comid_col)))
    return reaches[reaches[comid_col] == best].iloc[0]


def _rank_gages_from_outlet(gages_in, reaches_in, outlet, outlet_comid, comid_col=COMID_COL):
    """Order gages by proximity to the outlet: connected gages by fewest
    reaches upstream (hops), ties broken by distance, then disconnected gages
    by distance. Returns [(gage_id, FLComID), ...]."""
    if len(gages_in) == 0:
        return []
    if outlet is None or outlet_comid is None:
        return [(r["SOURCE_FEA"], int(r["FLComID"])) for _, r in gages_in.iterrows()]
    hops = _hops_from_outlet(reaches_in, outlet_comid, comid_col)
    cand = gages_in.copy()
    cand["_hops"] = cand["FLComID"].map(hops)
    cand["_dist"] = cand.geometry.distance(outlet.geometry)
    connected = cand.dropna(subset=["_hops"]).sort_values(["_hops", "_dist"])
    disconnected = cand[cand["_hops"].isna()].sort_values("_dist")
    ordered = pd.concat([connected, disconnected])
    return [(r["SOURCE_FEA"], int(r["FLComID"])) for _, r in ordered.iterrows()]


def _downstream_gages_from_outlet(
    outlet_comid, flow_all, gages_all, max_distance_km=DOWNSTREAM_GAGE_MAX_KM, comid_col=COMID_COL
):
    """USGS gages downstream of the outlet within max_distance_km of
    river-network distance. Returns [(gage_id, FLComID, distance_km), ...],
    closest first."""
    if outlet_comid is None:
        return []

    max_dist_m = float(max_distance_km) * 1000.0
    reaches = flow_all.set_index(comid_col, drop=False)
    if int(outlet_comid) not in reaches.index:
        return []

    fromnode_to_comids = {}
    for _, r in flow_all.iterrows():
        fromnode_to_comids.setdefault(int(r["FromNode"]), []).append(int(r[comid_col]))

    gages_by_comid = {}
    for _, g in gages_all.dropna(subset=["FLComID"]).iterrows():
        cid = int(g["FLComID"])
        gages_by_comid.setdefault(cid, []).append(g)

    outlet = reaches.loc[int(outlet_comid)]
    if isinstance(outlet, pd.DataFrame):
        outlet = outlet.iloc[0]
    start_node = int(outlet["ToNode"])

    q = deque([(start_node, 0.0)])
    best_node_distance = {start_node: 0.0}
    candidates = []

    while q:
        node, dist_to_node = q.popleft()
        for cid in fromnode_to_comids.get(node, []):
            r = reaches.loc[cid]
            if isinstance(r, pd.DataFrame):
                r = r.iloc[0]
            reach_length = float(r.geometry.length)

            for g in gages_by_comid.get(cid, []):
                try:
                    along_reach = float(r.geometry.project(g.geometry))
                except Exception:
                    along_reach = 0.0
                gage_dist = dist_to_node + along_reach
                if gage_dist <= max_dist_m:
                    candidates.append((g["SOURCE_FEA"], cid, gage_dist / 1000.0))

            new_dist = dist_to_node + reach_length
            if new_dist > max_dist_m:
                continue
            next_node = int(r["ToNode"])
            if new_dist < best_node_distance.get(next_node, float("inf")):
                best_node_distance[next_node] = new_dist
                q.append((next_node, new_dist))

    candidates.sort(key=lambda x: x[2])
    out, seen = [], set()
    for gid, cid, dist_km in candidates:
        key = str(gid)
        if key in seen:
            continue
        seen.add(key)
        out.append((gid, cid, dist_km))
    return out


def ri_years(label):
    """
    A recurrence-interval label as years, for fields whose name already says
    [yr]: "25 yr" -> 25, "> 100 yr" -> "> 100", "< 2 yr" -> "< 2". Anything
    else (None, "No USGS Gauge", an NWM failure reason) is returned unchanged.
    """
    if not isinstance(label, str):
        return label
    m = re.fullmatch(r"\s*([<>]?)\s*(\d+)\s*yr\s*", label)
    if not m:
        return label
    sign, years = m.groups()
    return f"{sign} {years}" if sign else int(years)


def _classify_ri(discharge_cms, ri_map):
    """Match a discharge to a recurrence interval. ri_map = {year: discharge_cms}.
    Returns (label, year_or_None, ri_discharge_cms)."""
    years = sorted(ri_map)
    lo, hi = years[0], years[-1]
    if discharge_cms > ri_map[hi]:
        return f"> {hi} yr", None, float(ri_map[hi])
    if discharge_cms < ri_map[lo]:
        return f"< {lo} yr", None, float(ri_map[lo])
    best = min(years, key=lambda yr: abs(ri_map[yr] - discharge_cms))
    return f"{best} yr", best, float(ri_map[best])


# ---------------------------------------------------------------------------
# USGS: observed discharge + published flood-frequency (StreamStats)
# ---------------------------------------------------------------------------


def _usgs_discharge_series_cfs(site_id, start, end):
    from dataretrieval import nwis

    site = str(site_id).strip()
    if site.isdigit() and len(site) < 8:
        site = site.zfill(8)
    try:
        df, _ = nwis.get_iv(sites=site, start=start, end=end, parameterCd="00060")
    except Exception as e:
        return None, site, f"USGS request failed for site {site}: {e}"
    if df is None or df.empty:
        return None, site, f"no instantaneous discharge at site {site}"
    cols = [c for c in df.columns if "00060" in c and not str(c).endswith("_cd")]
    if not cols:
        return None, site, f"site {site}: no 00060 discharge column returned"
    series = pd.to_numeric(df[cols[0]], errors="coerce").dropna()
    if series.empty:
        return None, site, f"site {site}: discharge empty/NaN"
    return series, site, None


def fetch_usgs_max_over_range_cms(site_id, start_date, end_date):
    """HWM: USGS 00060 instantaneous discharge -> hourly maxima -> the max
    over [start_date, end_date) (end_date exclusive -- pass end_date + 1 day
    to include a full calendar day). Returns (value_cms, time_used_utc, reason)."""
    start = start_date.strftime("%Y-%m-%d")
    end = end_date.strftime("%Y-%m-%d")
    series, site, reason = _usgs_discharge_series_cfs(site_id, start, end)
    if series is None:
        return None, None, reason
    hourly = series.resample("1h").max().dropna()
    if hourly.empty:
        return None, None, f"site {site}: no hourly discharge between {start} and {end}"
    peak_time = hourly.idxmax()
    peak_time = (
        peak_time.tz_localize("UTC") if peak_time.tzinfo is None else peak_time.tz_convert("UTC")
    )
    return float(hourly.max()) * CFS_TO_CMS, peak_time, None


def fetch_usgs_max_hourly_cms(site_id, date):
    """Tier 1: USGS 00060 instantaneous discharge -> hourly maxima -> daily max.
    Returns (value_cms, time_used_utc, reason)."""
    return fetch_usgs_max_over_range_cms(site_id, date, date + timedelta(days=1))


def fetch_usgs_at_time_cms(site_id, event_time, search_window_hours=USGS_AT_TIME_WINDOW_HOURS):
    """Tier 2/3: USGS 00060 observation nearest the exact UTC event_time.
    A +/- search_window_hours request window is used because USGS IV
    observations are usually discrete (often 5- or 15-minute), so an exact
    HH:MM:SS match may not exist. Returns (value_cms, time_used_utc, reason)."""
    event_time = pd.Timestamp(event_time)
    event_time = (
        event_time.tz_localize("UTC") if event_time.tzinfo is None else event_time.tz_convert("UTC")
    )
    start = (event_time - pd.Timedelta(hours=search_window_hours)).strftime("%Y-%m-%dT%H:%M:%SZ")
    end = (event_time + pd.Timedelta(hours=search_window_hours)).strftime("%Y-%m-%dT%H:%M:%SZ")
    series, _site, reason = _usgs_discharge_series_cfs(site_id, start, end)
    if series is None:
        return None, None, reason
    idx = pd.DatetimeIndex(series.index)
    idx = idx.tz_localize("UTC") if idx.tz is None else idx.tz_convert("UTC")
    series = series.copy()
    series.index = idx
    pos = int(np.argmin(np.abs((series.index - event_time).asi8)))
    used_time = series.index[pos]
    value_cms = float(series.iloc[pos]) * CFS_TO_CMS
    return value_cms, used_time, None


def fetch_usgs_event_discharge(site_id, event_dt, at_time):
    """Dispatch USGS discharge retrieval: the Tier 1 daily max (at_time=False,
    event_dt a naive midnight datetime) or the Tier 2/3 nearest-observation
    lookup (at_time=True, event_dt a UTC-aware pd.Timestamp)."""
    if at_time:
        return fetch_usgs_at_time_cms(site_id, event_dt)
    return fetch_usgs_max_hourly_cms(site_id, event_dt)


def _is_peakflow_group(grp):
    """True for the Peak-Flow Statistics group (guards against Flow-Duration
    '% duration' stats, which are not flood quantiles)."""
    if not grp:
        return False
    code = str(grp.get("code", "")).strip().upper()
    name = str(grp.get("name", "")).lower()
    return code == "PFS" or ("peak" in name and "flow" in name)


def _aep_percent(rt):
    """Annual-exceedance-probability percent for a peak-flow regressionType
    blob, e.g. '2-percent AEP flood' -> 2.0; code 'PK0_2AEP' -> 0.2."""
    name = str(rt.get("name", "")).lower()
    if "aep" in name:
        m = re.search(r"([\d.]+)\s*-?\s*percent", name)
        if m:
            return float(m.group(1))
    code = str(rt.get("code", "")).strip().upper()
    m = re.match(r"^PK([0-9_]+)AEP$", code)
    if m:
        return float(m.group(1).replace("_", "."))
    return None


def _is_cfs(unit):
    """USGS publishes peak-flow stats in cfs, so assume cfs when unspecified."""
    if not unit:
        return True
    text = " ".join(str(unit.get(k, "")) for k in ("abbreviation", "unit", "name")).lower()
    if "m^3" in text or "cubic meter" in text or "cms" in text:
        return False
    return True


@lru_cache(maxsize=512)
def _get_gagestats_stats(site):
    """Published statistic records for a station from GageStatsServices."""
    headers = {"Accept": "application/json"}
    try:
        r = requests.get(f"{GAGESTATS_BASE}/stations/{site}", headers=headers, timeout=60)
        if r.status_code == 200:
            data = r.json()
            if isinstance(data, dict) and data.get("statistics"):
                return data["statistics"]
    except Exception:
        pass
    r = requests.get(
        f"{GAGESTATS_BASE}/statistics",
        params={"stationIDOrCode": site},
        headers=headers,
        timeout=60,
    )
    r.raise_for_status()
    data = r.json()
    return data if isinstance(data, list) else []


def fetch_usgs_ffa_ri(site_id, ri_years=RI_YEARS):
    """Published peak-flow flood-frequency quantiles (cms) for a USGS gage,
    from the USGS StreamStats GageStatsServices API (published B17C values).
    Returns ({year: cms}, None) or (None, reason)."""
    site = str(site_id).strip()
    if site.isdigit() and len(site) < 8:
        site = site.zfill(8)

    try:
        stats_list = _get_gagestats_stats(site)
    except Exception as e:
        return None, f"flood-frequency request to USGS StreamStats failed ({e})"
    if not stats_list:
        return None, "no published flood-frequency statistics (USGS StreamStats)"

    wanted = set(ri_years)
    ri_cms = {}
    for s in stats_list:
        if not _is_peakflow_group(s.get("statisticGroupType")):
            continue
        aep = _aep_percent(s.get("regressionType") or {})
        if not aep or aep <= 0:
            continue
        ri = 100.0 / aep
        year = int(round(ri))
        if abs(ri - year) > 1e-6 or year not in wanted:
            continue
        val = s.get("value")
        if val is None:
            continue
        q = float(val)
        if _is_cfs(s.get("unitType")):
            q *= CFS_TO_CMS
        ri_cms.setdefault(year, q)

    if len(ri_cms) < 2:
        years = "/".join(map(str, sorted(wanted)))
        if not ri_cms:
            return None, f"no published {years}-yr peak-flow quantiles (USGS StreamStats)"
        return None, (
            f"only {len(ri_cms)} of the {years}-yr peak-flow quantiles published "
            "(USGS StreamStats; need >= 2)"
        )
    return ri_cms, None


@lru_cache(maxsize=512)
def _usgs_discharge_catalog(site):
    """USGS NWIS series catalog for discharge (00060) at a site:
    [(data_type_cd, begin, end)] with begin/end as pd.Timestamp (NaT if
    unparseable). data_type_cd is "uv"/"iv" (continuous instantaneous),
    "dv" (daily values) or "qw"/"sv"/"pk"/... (discrete measurements)."""
    import time

    from dataretrieval import nwis

    for attempt in range(USGS_CATALOG_RETRIES):
        try:
            df, _ = nwis.get_info(sites=site, seriesCatalogOutput="true")
            break
        except Exception:
            if attempt == USGS_CATALOG_RETRIES - 1:
                raise
            time.sleep(2**attempt)  # NWIS intermittently returns HTTP 503
    if df is None or df.empty or "parm_cd" not in df:
        return ()
    df = df[df["parm_cd"].astype(str) == "00060"]
    return tuple(
        (
            str(r["data_type_cd"]).lower(),
            pd.to_datetime(r["begin_date"], errors="coerce"),
            pd.to_datetime(r["end_date"], errors="coerce"),
        )
        for _, r in df.iterrows()
    )


def _fmt_date(ts):
    return "unknown" if pd.isna(ts) else ts.strftime("%Y-%m-%d")


def diagnose_usgs_no_discharge(site_id, window_start, window_end, window_desc):
    """Explain why a gage returned no instantaneous discharge for the event
    window, using its NWIS series catalog: the gage doesn't measure discharge,
    only has discrete field measurements, was discontinued / not yet installed
    (continuous record ends before / begins after the event), or had a data
    gap during the event. window_start/window_end are naive-UTC datetimes;
    window_desc is a human phrase like "within +/-1 h of ...". Returns a
    short reason string."""
    site = str(site_id).strip()
    if site.isdigit() and len(site) < 8:
        site = site.zfill(8)
    try:
        catalog = _usgs_discharge_catalog(site)
    except Exception:
        return f"no discharge data {window_desc} (USGS site catalog unavailable to tell why)"

    w0 = pd.Timestamp(window_start).tz_localize(None).normalize()
    w1 = pd.Timestamp(window_end).tz_localize(None).normalize()

    if not catalog:
        return "gage does not record discharge"

    continuous = [(b, e) for cd, b, e in catalog if cd in ("uv", "iv")]
    daily = [(b, e) for cd, b, e in catalog if cd == "dv"]

    if not continuous:
        if daily:
            b = min(b for b, _ in daily)
            e = max(e for _, e in daily)
            if not pd.isna(e) and e < w0:
                return (
                    f"gage discontinued (daily discharge record {_fmt_date(b)} to {_fmt_date(e)})"
                )
            return (
                f"no continuous (instantaneous) discharge record; daily values only "
                f"({_fmt_date(b)} to {_fmt_date(e)})"
            )
        return "no continuous discharge record (discrete field measurements only)"

    b = min(b for b, _ in continuous)
    e = max(e for _, e in continuous)
    if not pd.isna(e) and e < w0:
        return f"gage discontinued (continuous discharge record ended {_fmt_date(e)})"
    if not pd.isna(b) and b > w1:
        return f"gage not yet operating (continuous discharge record begins {_fmt_date(b)})"
    return f"gage not reporting {window_desc} (data gap in continuous discharge record)"


def _resolve_usgs_ri(
    geom,
    geom_crs,
    gages_shp,
    flow_shp,
    fetch_discharge,
    window,
    proj=5070,
    downstream_max_km=DOWNSTREAM_GAGE_MAX_KM,
    log=None,
):
    """Shared USGS gage resolution + FFA classification, used by both the
    single-event (get_event_recurrence_interval_usgs) and date-range
    (get_event_recurrence_interval_usgs_range) lookups. Everything about
    *which* gage to use (in-polygon, closest to the outlet first, else the
    nearest usable gage downstream; Scenario 1/2/3 labeling) is identical
    either way -- only how a candidate gage's discharge is fetched differs,
    supplied as ``fetch_discharge(gage_id) -> (q_cms, used_time, reason)``.
    ``window`` is (start, end, description) of the period discharge is
    requested for; it is used to explain why a gage had no discharge.

    Returns (label, discharge_cms, scenario, gage_id) -- see
    get_event_recurrence_interval_usgs for the exact contract.
    """
    log = log or _noop_log

    if not gages_shp or not flow_shp:
        log("Event Recurrence Interval (USGS): gages_shp/flow_shp not configured; skipping.")
        return None, None, None, None
    if geom is None or geom.is_empty:
        return None, None, None, None

    try:
        geom_proj = gpd.GeoSeries([geom], crs=geom_crs).to_crs(proj).iloc[0]
        # The downstream walk stays within downstream_max_km of river from the
        # outlet, so features beyond that (plus slack for the outlet reach
        # itself) can never be used.
        region = (geom_proj, float(downstream_max_km) * 1000.0 + REGION_SLACK_M)
        gages = _load_gages(gages_shp, proj, region)
        flow = _load_flow(flow_shp, proj, region)

        gages_in = _clip_to_geom(gages, geom_proj, proj)
        flow_in = _clip_to_geom(flow, geom_proj, proj)

        outlet = _pick_outlet_reach(flow_in)
        outlet_comid = None if outlet is None else int(outlet[COMID_COL])

        try_order, seen = [], set()
        for gid, flc in _rank_gages_from_outlet(gages_in, flow_in, outlet, outlet_comid):
            if gid is None or str(gid) in seen:
                continue
            seen.add(str(gid))
            try_order.append((gid, flc))
        n_in_polygon = len(try_order)

        skip_reasons = []

        def _why_no_discharge(gid, q_reason):
            if q_reason and q_reason.startswith("USGS request failed"):
                return q_reason
            return diagnose_usgs_no_discharge(gid, *window)

        # 1) Try every gage inside the polygon/corridor, closest to the outlet first.
        for rank_idx, (gid, _flc) in enumerate(try_order):
            q_cms, _used_time, q_reason = fetch_discharge(gid)
            if q_cms is None:
                skip_reasons.append(f"{gid}: {_why_no_discharge(gid, q_reason)}")
                continue
            ffa, ffa_reason = fetch_usgs_ffa_ri(gid)
            if ffa is None:
                skip_reasons.append(f"{gid}: {ffa_reason}")
                continue
            label, _yr, _q_ri = _classify_ri(q_cms, ffa)
            log(f"Event Recurrence Interval (USGS): gage {gid} -> {q_cms:.2f} cms -> {label}")
            if rank_idx == 0:
                scenario = (
                    f"Scenario 1: single in-polygon gage {gid}"
                    if n_in_polygon == 1
                    else f"Scenario 2: {n_in_polygon} in-polygon gages, fewest reaches from "
                    f"outlet -> gage {gid}"
                )
            else:
                scenario = (
                    f"Scenario 2: primary in-polygon gage unusable; used fallback gage {gid} "
                    f"(proximity rank {rank_idx})"
                )
            return label, float(q_cms), scenario, str(gid)

        # 2) If no in-polygon gage is usable, search downstream from the outlet.
        if outlet_comid is not None:
            for gid, _flc, distance_km in _downstream_gages_from_outlet(
                outlet_comid, flow, gages, max_distance_km=downstream_max_km
            ):
                if str(gid) in seen:
                    continue
                seen.add(str(gid))
                q_cms, _used_time, q_reason = fetch_discharge(gid)
                if q_cms is None:
                    skip_reasons.append(
                        f"{gid} ({distance_km:.1f} km downstream): "
                        f"{_why_no_discharge(gid, q_reason)}"
                    )
                    continue
                ffa, ffa_reason = fetch_usgs_ffa_ri(gid)
                if ffa is None:
                    skip_reasons.append(f"{gid} ({distance_km:.1f} km downstream): {ffa_reason}")
                    continue
                label, _yr, _q_ri = _classify_ri(q_cms, ffa)
                log(
                    f"Event Recurrence Interval (USGS): downstream gage {gid} "
                    f"({distance_km:.1f} km) -> {q_cms:.2f} cms -> {label}"
                )
                scenario = (
                    f"Scenario 2: no usable gage in polygon ({n_in_polygon} found); used "
                    f"downstream gage {gid} ({distance_km:.1f} km downstream of outlet)"
                )
                return label, float(q_cms), scenario, str(gid)

        label = NO_USABLE_GAGE_LABEL if skip_reasons else NO_GAGE_LABEL
        detail = "; ".join(skip_reasons) if skip_reasons else "no gages found"
        log(
            f"Event Recurrence Interval (USGS): {label} in polygon or within "
            f"{downstream_max_km:.0f} km downstream ({detail})"
        )
        scenario = (
            f"Scenario 3: no usable USGS gage in polygon ({n_in_polygon} found) or within "
            f"{downstream_max_km:.0f} km downstream -- use NWM"
        )
        return label, None, scenario, None
    except Exception as e:
        log(f"Event Recurrence Interval (USGS): lookup failed: {e}")
        return None, None, None, None


def get_event_recurrence_interval_usgs(
    geom,
    geom_crs,
    flood_date,
    gages_shp,
    flow_shp,
    proj=5070,
    downstream_max_km=DOWNSTREAM_GAGE_MAX_KM,
    log=None,
):
    """
    Event Recurrence Interval (USGS): resolve ``geom`` to a USGS streamgage
    (in-polygon, closest to the outlet first; else the nearest usable gage
    downstream along the flow network within ``downstream_max_km``), fetch its
    observed discharge for ``flood_date``, and classify it against that
    gage's published USGS flood-frequency quantiles.

    ``flood_date`` is "YYYYMMDD" (Tier 1: daily max hourly discharge) or
    "YYYYMMDDTHHMMSS" UTC (Tier 2/3: observation nearest that exact instant)
    -- see _parse_flood_datetime.

    Returns (label, discharge_cms, scenario, gage_id):
      - label is "25 yr" / "> 100 yr" / "< 2 yr"; ``NO_GAGE_LABEL`` ("No
        USGS Gauge") if no gage was found; ``NO_USABLE_GAGE_LABEL`` ("No
        usable USGS Gauge") if gages were found but none had usable
        discharge/FFA (why each was rejected is logged); or
        None if the lookup could not even be attempted (missing
        gages_shp/flow_shp, empty geometry, unparseable flood_date, or an
        unexpected error).
      - discharge_cms is the gage's observed discharge (cms) that produced
        ``label``, or None whenever label is not a classified return period.
      - scenario describes which resolution path was used, mirroring the
        Scenario 1/2/3 terminology from fim_recurrence_pipeline_T1.ipynb:
        Scenario 1 (single in-polygon gage), Scenario 2 (multiple in-polygon
        gages / an in-polygon or downstream fallback), or Scenario 3 (no
        usable gage anywhere -- NWM should be used instead). None when the
        lookup could not be attempted at all.
      - gage_id is the USGS site number (str) of the gage that produced
        ``label``, or None whenever label is not a classified return period.
    All cases are logged via ``log``.
    """
    log = log or _noop_log

    if not gages_shp or not flow_shp:
        log("Event Recurrence Interval (USGS): gages_shp/flow_shp not configured; skipping.")
        return None, None, None, None
    if geom is None or geom.is_empty:
        return None, None, None, None

    event_dt, has_time = _parse_flood_datetime(flood_date, log=log)
    if event_dt is None:
        log(f"Event Recurrence Interval (USGS): could not parse flood_date {flood_date!r}.")
        return None, None, None, None

    if has_time:
        t = pd.Timestamp(event_dt).tz_convert(None)
        h = pd.Timedelta(hours=USGS_AT_TIME_WINDOW_HOURS)
        window = (
            t - h,
            t + h,
            f"within +/-{USGS_AT_TIME_WINDOW_HOURS} h of {t:%Y-%m-%d %H:%M} UTC",
        )
    else:
        window = (event_dt, event_dt, f"on {event_dt:%Y-%m-%d}")

    return _resolve_usgs_ri(
        geom,
        geom_crs,
        gages_shp,
        flow_shp,
        fetch_discharge=lambda gid: fetch_usgs_event_discharge(gid, event_dt, has_time),
        window=window,
        proj=proj,
        downstream_max_km=downstream_max_km,
        log=log,
    )


def get_event_recurrence_interval_usgs_range(
    geom,
    geom_crs,
    start_date,
    end_date,
    gages_shp,
    flow_shp,
    proj=5070,
    downstream_max_km=DOWNSTREAM_GAGE_MAX_KM,
    log=None,
):
    """
    Event Recurrence Interval (USGS) over a date range, for High Water Marks:
    resolve ``geom`` to a USGS streamgage exactly as
    get_event_recurrence_interval_usgs does, but instead of a single event's
    discharge, use the MAXIMUM observed discharge anywhere in
    [start_date, end_date] inclusive (USGS 00060 instantaneous -> hourly
    maxima -> max over the whole window), then classify that peak against the
    gage's published USGS flood-frequency quantiles.

    start_date / end_date are "YYYYMMDD" strings (only the first 8 characters
    are used). Swapped automatically if end_date is earlier than start_date.

    Returns (label, discharge_cms, scenario, gage_id) -- same contract as
    get_event_recurrence_interval_usgs.
    """
    log = log or _noop_log

    if not gages_shp or not flow_shp:
        log("Event Recurrence Interval (USGS): gages_shp/flow_shp not configured; skipping.")
        return None, None, None, None
    if geom is None or geom.is_empty:
        return None, None, None, None

    start_dt = _parse_date_only(start_date, log=log)
    end_dt = _parse_date_only(end_date, log=log)
    if start_dt is None or end_dt is None:
        log(
            "Event Recurrence Interval (USGS): could not parse start_date/end_date "
            f"({start_date!r}, {end_date!r})."
        )
        return None, None, None, None
    if end_dt < start_dt:
        start_dt, end_dt = end_dt, start_dt
    range_end = end_dt + timedelta(days=1)  # include the full end_date calendar day

    return _resolve_usgs_ri(
        geom,
        geom_crs,
        gages_shp,
        flow_shp,
        fetch_discharge=lambda gid: fetch_usgs_max_over_range_cms(gid, start_dt, range_end),
        window=(start_dt, end_dt, f"between {start_dt:%Y-%m-%d} and {end_dt:%Y-%m-%d}"),
        proj=proj,
        downstream_max_km=downstream_max_km,
        log=log,
    )


# ---------------------------------------------------------------------------
# NWM: modeled discharge at the outlet reach + NWM v2.1 RI tables
# ---------------------------------------------------------------------------


def fetch_nwm_retro_max_over_range(comid, start_date, end_date, zarr_url=NWM_RETRO_ZARR):
    """HWM: max hourly NWM retrospective streamflow (cms) over [start_date,
    end_date] (inclusive of end_date's instant -- pass end_date at 23:59:59
    to include a full calendar day). Returns (value_cms, time_used_utc, reason)."""
    import fsspec
    import xarray as xr

    try:
        ds = xr.open_zarr(fsspec.get_mapper(zarr_url, anon=True), consolidated=True)
    except Exception as e:
        return None, None, f"could not open retrospective zarr: {e}"
    try:
        window = ds["streamflow"].sel(feature_id=int(comid), time=slice(start_date, end_date))
    except KeyError:
        return None, None, f"COMID {comid} not in retrospective"
    if int(window.time.size) == 0:
        return (
            None,
            None,
            (
                f"no retrospective steps between {start_date:%Y-%m-%d} and {end_date:%Y-%m-%d} "
                f"at {comid}"
            ),
        )
    vals = np.asarray(window.values)
    if vals.size == 0 or np.all(np.isnan(vals)):
        return (
            None,
            None,
            (
                f"retrospective discharge empty/NaN between {start_date:%Y-%m-%d} and "
                f"{end_date:%Y-%m-%d} at {comid}"
            ),
        )
    pos = int(np.nanargmax(vals))
    used_time = pd.Timestamp(window.time.values[pos], tz="UTC")
    return float(vals[pos]), used_time, None


def fetch_nwm_retro_max_hourly(comid, date, zarr_url=NWM_RETRO_ZARR):
    """Tier 1: max hourly NWM retrospective streamflow (cms) for the event day.
    Returns (value_cms, time_used_utc, reason)."""
    return fetch_nwm_retro_max_over_range(
        comid, date, date + timedelta(hours=23, minutes=59, seconds=59), zarr_url=zarr_url
    )


def fetch_nwm_shortrange_max_over_range(comid, start_date, end_date, bucket=NWM_OPS_BUCKET):
    """HWM: max hourly NWM short-range streamflow (cms) valid over
    [start_date, end_date] (inclusive, by calendar day). Returns (value_cms,
    time_used_utc, reason)."""
    import fsspec
    import xarray as xr

    fs = fsspec.filesystem("s3", anon=True)
    day = datetime(start_date.year, start_date.month, start_date.day)
    last_day = datetime(end_date.year, end_date.month, end_date.day)
    peak, peak_time, found = None, None, False
    while day <= last_day:
        ymd = day.strftime("%Y%m%d")
        for hh in range(24):
            key = f"{bucket}/nwm.{ymd}/short_range/nwm.t{hh:02d}z.short_range.channel_rt.f001.conus.nc"
            if not fs.exists(key):
                continue
            try:
                with fs.open(key) as fobj:
                    ds = xr.open_dataset(fobj)
                    v = float(ds["streamflow"].sel(feature_id=int(comid)).values)
                found = True
                # f001 valid time = initialization + 1 hour
                valid_time = pd.Timestamp(
                    datetime.strptime(ymd + f"{hh:02d}", "%Y%m%d%H"), tz="UTC"
                ) + pd.Timedelta(hours=1)
                if peak is None or v > peak:
                    peak, peak_time = v, valid_time
            except Exception:
                continue
        day += timedelta(days=1)
    if not found:
        return (
            None,
            None,
            (
                f"no short-range files between {start_date:%Y-%m-%d} and {end_date:%Y-%m-%d} "
                f"at {comid} (archive expired?)"
            ),
        )
    return peak, peak_time, None


def fetch_nwm_shortrange_max_hourly(comid, date, bucket=NWM_OPS_BUCKET):
    """Tier 1: max hourly NWM short-range streamflow (cms) valid on the event
    day. Returns (value_cms, time_used_utc, reason)."""
    return fetch_nwm_shortrange_max_over_range(comid, date, date, bucket=bucket)


def fetch_nwm_retro_at_time(comid, event_time, zarr_url=NWM_RETRO_ZARR):
    """Tier 2/3: NWM retrospective streamflow at the model timestep nearest
    the exact UTC event_time. Returns (value_cms, time_used_utc, reason)."""
    import fsspec
    import xarray as xr

    event_time = pd.Timestamp(event_time)
    event_time = (
        event_time.tz_localize("UTC") if event_time.tzinfo is None else event_time.tz_convert("UTC")
    )
    target = event_time.tz_localize(None)  # NWM xarray time coordinate is UTC-naive

    try:
        ds = xr.open_zarr(fsspec.get_mapper(zarr_url, anon=True), consolidated=True)
    except Exception as e:
        return None, None, f"could not open retrospective zarr: {e}"
    try:
        q = ds["streamflow"].sel(feature_id=int(comid), time=target, method="nearest")
    except KeyError:
        return None, None, f"COMID {comid} not in retrospective"
    except Exception as e:
        return None, None, f"could not retrieve retrospective flow for COMID {comid}: {e}"
    used_time = pd.Timestamp(q["time"].values, tz="UTC")
    return float(q.values), used_time, None


def fetch_nwm_shortrange_at_time(comid, event_time, bucket=NWM_OPS_BUCKET):
    """Tier 2/3: NWM short-range streamflow nearest the exact UTC event_time.
    NWM channel output is hourly; the event is rounded to the nearest valid
    hour, using the f001 convention (initialization one hour before valid
    time). Returns (value_cms, time_used_utc, reason)."""
    import fsspec
    import xarray as xr

    event_time = pd.Timestamp(event_time)
    event_time = (
        event_time.tz_localize("UTC") if event_time.tzinfo is None else event_time.tz_convert("UTC")
    )
    valid_time = event_time.round("1h")
    init_time = valid_time - pd.Timedelta(hours=1)
    ymd = init_time.strftime("%Y%m%d")
    hh = init_time.hour
    fs = fsspec.filesystem("s3", anon=True)
    key = f"{bucket}/nwm.{ymd}/short_range/nwm.t{hh:02d}z.short_range.channel_rt.f001.conus.nc"
    if not fs.exists(key):
        return (
            None,
            None,
            (
                f"no short-range file for valid time {valid_time} at COMID {comid} (archive expired?)"
            ),
        )
    try:
        with fs.open(key) as fobj:
            ds = xr.open_dataset(fobj)
            value = float(ds["streamflow"].sel(feature_id=int(comid)).values)
    except Exception as e:
        return None, None, f"could not retrieve NWM short-range flow for COMID {comid}: {e}"
    return value, valid_time, None


def fetch_nwm_event_discharge(comid, event_dt, retro_max_year=RETRO_MAX_YEAR, at_time=False):
    """Dispatch NWM discharge retrieval: the Tier 1 daily max (at_time=False,
    event_dt a naive midnight datetime) or the Tier 2/3 nearest-timestep
    lookup (at_time=True, event_dt a UTC-aware pd.Timestamp). Either way,
    events through retro_max_year use NWM Retrospective v3.0; later events use
    the NWM short-range archive. Returns (value_cms, time_used_utc, reason, source)."""
    if at_time:
        if event_dt.year <= retro_max_year:
            return fetch_nwm_retro_at_time(comid, event_dt) + ("NWM retrospective v3.0",)
        return fetch_nwm_shortrange_at_time(comid, event_dt) + ("NWM short-range",)
    if event_dt.year <= retro_max_year:
        return fetch_nwm_retro_max_hourly(comid, event_dt) + ("NWM retrospective v3.0",)
    return fetch_nwm_shortrange_max_hourly(comid, event_dt) + ("NWM short-range",)


def fetch_nwm_max_over_range_discharge(comid, start_date, end_date, retro_max_year=RETRO_MAX_YEAR):
    """HWM: max NWM streamflow over [start_date, end_date]. Uses NWM
    Retrospective v3.0 when end_date.year <= retro_max_year, else the
    short-range archive for the whole range -- a range straddling the
    boundary is resolved by its later (and likelier to be short-range-only)
    end. Returns (value_cms, time_used_utc, reason, source)."""
    if end_date.year <= retro_max_year:
        full_end = datetime(end_date.year, end_date.month, end_date.day, 23, 59, 59)
        return fetch_nwm_retro_max_over_range(comid, start_date, full_end) + (
            "NWM retrospective v3.0",
        )
    return fetch_nwm_shortrange_max_over_range(comid, start_date, end_date) + ("NWM short-range",)


def load_ri_tables(folder, feature_ids, ri_years=RI_YEARS, log=None):
    """Read the NWM RI CSVs (one per return period), filtered to the given
    COMIDs. Returns {year: Series(discharge_cms indexed by feature_id)}."""
    log = log or _noop_log
    ids = {int(x) for x in feature_ids if x is not None}
    tables, missing_files = {}, []
    for yr in ri_years:
        path = os.path.join(folder, f"nwm21_17C_recurr_{yr}_0_cms.csv")
        if not os.path.exists(path):
            missing_files.append(os.path.basename(path))
            continue
        df = pd.read_csv(path, usecols=["feature_id", "discharge"])
        s = (
            df[df["feature_id"].isin(ids)]
            .dropna(subset=["discharge"])
            .drop_duplicates(subset="feature_id", keep="last")
            .set_index("feature_id")["discharge"]
        )
        tables[yr] = s
    if missing_files:
        log(
            "Event Recurrence Interval (NWM): missing RI CSV(s), those return "
            f"periods unavailable: {', '.join(missing_files)}"
        )
    return tables


def ri_from_nwm(discharge_cms, comid, ri_tables):
    """Classify a NWM discharge against every RI table the COMID appears in.
    Returns (label, reason): label is None with a reason when the COMID is
    absent from every RI table."""
    ri_map = {}
    for yr, s in ri_tables.items():
        if comid in s.index:
            ri_map[yr] = float(s.loc[comid])
    if not ri_map:
        return None, f"COMID {comid} not in any RI table"
    label, _yr, _q_ri = _classify_ri(discharge_cms, ri_map)
    return label, None


def _resolve_nwm_ri(
    geom,
    geom_crs,
    flow_shp,
    recurr_folder,
    fetch_discharge,
    proj=5070,
    ri_years=RI_YEARS,
    log=None,
):
    """Shared NWM outlet resolution + RI classification, used by both the
    single-event (get_event_recurrence_interval_nwm) and date-range
    (get_event_recurrence_interval_nwm_range) lookups. Only how the outlet
    reach's discharge is fetched differs, supplied as
    ``fetch_discharge(comid) -> (q_cms, used_time, reason, source)``.

    Returns (label, discharge_cms, scenario) -- see
    get_event_recurrence_interval_nwm for the exact contract.
    """
    log = log or _noop_log

    if not flow_shp or not recurr_folder:
        log("Event Recurrence Interval (NWM): flow_shp/recurr_folder not configured; skipping.")
        reason = "NWM not configured"
        return reason, None, reason
    if geom is None or geom.is_empty:
        reason = "Invalid flood geometry"
        return reason, None, reason

    try:
        geom_proj = gpd.GeoSeries([geom], crs=geom_crs).to_crs(proj).iloc[0]
        flow = _load_flow(flow_shp, proj, (geom_proj, GAP_TOL + REGION_SLACK_M))
        flow_in = _clip_to_geom(flow, geom_proj, proj)

        outlet = _pick_outlet_reach(flow_in)
        if outlet is None:
            log("Event Recurrence Interval (NWM): could not resolve an outlet reach.")
            reason = "No outlet reach found"
            return reason, None, reason
        comid = int(outlet[COMID_COL])

        q_cms, _used_time, reason, src = fetch_discharge(comid)
        if q_cms is None:
            log(
                f"Event Recurrence Interval (NWM): discharge fetch failed for COMID {comid}: {reason}"
            )
            short = f"NWM discharge unavailable: {reason}"
            return short, None, short

        ri_tables = load_ri_tables(recurr_folder, [comid], ri_years=ri_years, log=log)
        label, ri_reason = ri_from_nwm(q_cms, comid, ri_tables)
        if label is None:
            log(f"Event Recurrence Interval (NWM): {ri_reason}")
            short = f"NWM RI unavailable: {ri_reason}"
            return short, float(q_cms), short

        log(f"Event Recurrence Interval (NWM): {src} COMID {comid} -> {q_cms:.2f} cms -> {label}")
        scenario = f"{src} discharge at outlet reach COMID {comid}"
        return label, float(q_cms), scenario
    except Exception as e:
        log(f"Event Recurrence Interval (NWM): lookup failed: {e}")
        reason = f"NWM lookup failed: {e}"
        return reason, None, reason


def get_event_recurrence_interval_nwm(
    geom,
    geom_crs,
    flood_date,
    flow_shp,
    recurr_folder,
    proj=5070,
    retro_max_year=RETRO_MAX_YEAR,
    ri_years=RI_YEARS,
    log=None,
):
    """
    Event Recurrence Interval (NWM): resolve ``geom``'s outlet NHDPlus reach
    (COMID), fetch NWM streamflow there for ``flood_date``, and classify it
    against the NWM v2.1 recurrence-interval tables in ``recurr_folder``.

    ``flood_date`` is "YYYYMMDD" (Tier 1: daily max hourly streamflow) or
    "YYYYMMDDTHHMMSS" UTC (Tier 2/3: streamflow nearest that exact instant)
    -- see _parse_flood_datetime.

    Returns (label, discharge_cms, scenario, gage_id):
      - label is "25 yr" / "> 100 yr" / "< 2 yr" on success. When the RI
        can't be resolved, label carries a short reason instead of None (e.g.
        "NWM not configured", "No outlet reach found", "NWM discharge
        unavailable: ..."), so the metadata field always says *why* rather
        than being a bare null.
      - discharge_cms is the NWM discharge (cms) fetched at the outlet reach,
        or None whenever it could not be fetched. It can be set even when
        label carries a reason (discharge fetched fine, but the COMID has no
        RI table match).
      - scenario describes the resolution path (outlet COMID and NWM source
        used -- "NWM retrospective v3.0" or "NWM short-range", per
        fim_recurrence_pipeline_T1.ipynb) on success, or the same reason as
        label on failure.
    All cases are logged via ``log``.
    """
    log = log or _noop_log

    if not flow_shp or not recurr_folder:
        log("Event Recurrence Interval (NWM): flow_shp/recurr_folder not configured; skipping.")
        reason = "NWM not configured"
        return reason, None, reason
    if geom is None or geom.is_empty:
        reason = "Invalid flood geometry"
        return reason, None, reason

    event_dt, has_time = _parse_flood_datetime(flood_date, log=log)
    if event_dt is None:
        log(f"Event Recurrence Interval (NWM): could not parse flood_date {flood_date!r}.")
        reason = "Could not parse flood_date"
        return reason, None, reason

    return _resolve_nwm_ri(
        geom,
        geom_crs,
        flow_shp,
        recurr_folder,
        fetch_discharge=lambda comid: fetch_nwm_event_discharge(
            comid, event_dt, retro_max_year, at_time=has_time
        ),
        proj=proj,
        ri_years=ri_years,
        log=log,
    )


def get_event_recurrence_interval_nwm_range(
    geom,
    geom_crs,
    start_date,
    end_date,
    flow_shp,
    recurr_folder,
    proj=5070,
    retro_max_year=RETRO_MAX_YEAR,
    ri_years=RI_YEARS,
    log=None,
):
    """
    Event Recurrence Interval (NWM) over a date range, for High Water Marks:
    resolve ``geom``'s outlet NHDPlus reach exactly as
    get_event_recurrence_interval_nwm does, but instead of a single event's
    streamflow, use the MAXIMUM NWM streamflow anywhere in
    [start_date, end_date] inclusive, then classify that peak against the
    NWM v2.1 recurrence-interval tables.

    start_date / end_date are "YYYYMMDD" strings (only the first 8 characters
    are used). Swapped automatically if end_date is earlier than start_date.
    If the range straddles retro_max_year, the source is chosen by
    end_date's year (see fetch_nwm_max_over_range_discharge).

    Returns (label, discharge_cms, scenario) -- same contract as
    get_event_recurrence_interval_nwm.
    """
    log = log or _noop_log

    if not flow_shp or not recurr_folder:
        log("Event Recurrence Interval (NWM): flow_shp/recurr_folder not configured; skipping.")
        reason = "NWM not configured"
        return reason, None, reason
    if geom is None or geom.is_empty:
        reason = "Invalid flood geometry"
        return reason, None, reason

    start_dt = _parse_date_only(start_date, log=log)
    end_dt = _parse_date_only(end_date, log=log)
    if start_dt is None or end_dt is None:
        log(
            "Event Recurrence Interval (NWM): could not parse start_date/end_date "
            f"({start_date!r}, {end_date!r})."
        )
        reason = "Could not parse start_date/end_date"
        return reason, None, reason
    if end_dt < start_dt:
        start_dt, end_dt = end_dt, start_dt

    return _resolve_nwm_ri(
        geom,
        geom_crs,
        flow_shp,
        recurr_folder,
        fetch_discharge=lambda comid: fetch_nwm_max_over_range_discharge(
            comid, start_dt, end_dt, retro_max_year
        ),
        proj=proj,
        ri_years=ri_years,
        log=log,
    )
