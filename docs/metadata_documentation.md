# Flood map metadata

Every flood map in the FIMbench database has a `metadata.json` file stored
next to its raster. The file describes the raster itself, where it is, which
flood it shows, where it came from, and how severe the event was. This page
explains every metadata field, how its value is produced, and how it differs
between FIM categories. It also covers the files built from the metadata: the
AOI GeoPackage, the extent GeoJSON, and the web catalog.

The metadata is written by the processors in
[`fimbench.processing_floodmap`](../src/fimbench/processing_floodmap/), one per
category.

---

## Contents

1. [What each flood map consists of](#1-what-each-flood-map-consists-of)
2. [FIM categories at a glance](#2-fim-categories-at-a-glance)
3. [Which fields each category has](#3-which-fields-each-category-has)
4. [Field reference](#4-field-reference)
5. [Example metadata.json](#5-example-metadatajson)
6. [Companion files: AOI.gpkg and extent GeoJSON](#6-companion-files-aoigpkg-and-extent-geojson)
7. [From metadata to the web catalog](#7-from-metadata-to-the-web-catalog)
8. [Notes](#8-notes)

---

## 1. What each flood map consists of

Processing one raw flood raster produces a folder containing three files, plus
an extent GeoJSON written to a separate folder:

```
{folder}/
├── {name}_BM.tif           # the benchmark flood map (GeoTIFF, EPSG:4326)
├── {name}_metadata.json    # the metadata described on this page
└── {name}_AOI.gpkg         # area of interest: valid-data footprint + metadata
{geojson_dest}/
└── {name}_extent.geojson   # flooded area + key metadata
```

FEMA BLE folders also hold the flow table used for the map,
`BLE_HUC_{HUC8}_FLOWS_{return period}YR_{DMS}.csv`.

### Naming convention

| Category | Folder | File stem `{name}` |
|---|---|---|
| Tier 1, 2, 3 | `{SENSOR}_{date}_{DMS}` | `{SENSOR}_{res}_{date}_{DMS}` |
| HWM | `HWM_{start}_{end}_{DMS}` | `HWM_{res}_{start}_{end}_{DMS}` |
| FEMA BLE | `BLE_{return period}_{DMS}` | `BLE_{res}_{return period}_{DMS}` |

- `{SENSOR}` is the sensor code: `AI`, `PSS`, `REO`, `S1A`, `HWM` or `BLE`.
- `{res}` is the value of `Spatial resolution (m)` with `_` in place of the
  decimal point, for example `0_3m` or `10_0m`.
- `{date}` is the acquisition date as given, either `YYYYMMDD` or `YYYYMMDDTHHMMSS` (UTC).
- `{start}` and `{end}` are the HWM flood window, each as `YYYYMMDD`.
- `{DMS}` is the `Centroid DMS code`, for example `891934W365313N`.

Examples:

```
Tier_1/AI_20110412_965145W472341N/AI_0_5m_20110412_965145W472341N_metadata.json
Tier_2/PSS_20161009T150712_775851W352159N/PSS_3_0m_20161009T150712_775851W352159N_metadata.json
Tier_3/S1A_20170504T235530_891934W365313N/S1A_10_0m_20170504T235530_891934W365313N_metadata.json
HWM/HWM_20100329_20100401_710139W415659N/HWM_10_0m_20100329_20100401_710139W415659N_metadata.json
FEMA_BLE/BLE_100_1000404W323634N/BLE_10_0m_100_1000404W323634N_metadata.json
```

In S3, these folders sit under `FIM_Database/{Tier_1 | Tier_2 | Tier_3 | HWM | FEMA_BLE}/`.

---

## 2. FIM categories at a glance

| Category | Processor | Sensor code | `Full sensor name` | `Product Category` | `Spatial resolution (m)` | Event timing |
|---|---|---|---|---|---|---|
| **Tier 1**: Aerial imagery | `Tier1Processor` | `AI` | Aerial Imagery | `Tier 1` | Measured | Acquisition date |
| **Tier 2**: PlanetScope | `Tier2Processor` | `PSS` | Planet Scope Scene | `Tier 2` | Measured | Acquisition date + UTC time |
| **Tier 2**: RapidEye | `Tier2Processor` | `REO` | RapidEye Ortho Tile | `Tier 2` | Measured | Acquisition date + UTC time |
| **Tier 3**: Sentinel-1 | `Tier3Processor` | `S1A` | Sentinel-1A | `Tier 3` | Fixed `10` | Acquisition date + UTC time |
| **HWM**: High-water marks | `HwmProcessor` | `HWM` | High-water mark | `High-water mark` | Fixed `10` | Flood window (start–end) |
| **FEMA BLE**: Base Level Engineering | `FemaBleProcessor` | `BLE` | Base Level Engineering | `FEMA_BLE` | Fixed `10` | Return period (synthetic) |

`Data provider` per category:

| Category | `Data provider` |
|---|---|
| Tier 1 | Dr. Dinuke Munasinghe, The University of Alabama. Email: dsmunasinghe@ua.edu |
| Tier 2 | Dr. Dan Tian, The University of Alabama. Email: dtian1@ua.edu |
| Tier 3 | Dr. Dan Tian, The University of Alabama. Email ID: dtian1@ua.edu |
| HWM | Dr. Dinuke Munasinghe, The University of Alabama. Email ID: dsmunasinghe@ua.edu |
| FEMA BLE | NOAA/NWS Office of Water Prediction (OWP) |

The fixed values (`Full sensor name`, `Data provider`, `Product Category`,
`Keywords`, `Access Rights`, documentation links, `References`, and the fixed
resolution) live in each processor's `METADATA_FIELDS` table, keyed by their
field names. `full_form=`, `data_provider=` (or `source=`) and `quality=`
override `Full sensor name`, `Data provider` and `Product Category` when
creating a processor. For RapidEye maps, pass `sensor_code="REO"` and
`full_form="RapidEye Ortho Tile"` to `Tier2Processor`.

---

## 3. Which fields each category has

Fields appear in `metadata.json` in the order below. ✓ = present, — = not present.

| Field | Tier 1 | Tier 2 | Tier 3 | HWM | FEMA BLE |
|---|:-:|:-:|:-:|:-:|:-:|
| `File Name` | ✓ | ✓ | ✓ | ✓ | ✓ |
| `Full sensor name` | ✓ | ✓ | ✓ | ✓ | ✓ |
| `HUC8 ID(s)` | ✓ | ✓ | ✓ | ✓ | ✓ (dominant only) |
| `Format` | ✓ | ✓ | ✓ | ✓ | ✓ |
| `Nodata_Value` | ✓ | ✓ | ✓ | ✓ | ✓ |
| `Spatial resolution (m)` | ✓ | ✓ | ✓ | ✓ | ✓ |
| `Datatype` | ✓ | ✓ | ✓ | ✓ | ✓ |
| `Compression type` | ✓ | ✓ | ✓ | ✓ | ✓ |
| `Bounding box` | ✓ | ✓ | ✓ | ✓ | ✓ |
| `Centroid DMS code` | ✓ | ✓ | ✓ | ✓ | ✓ |
| `Coordinate reference system` | ✓ | ✓ | ✓ | ✓ | ✓ |
| `Number of bands` | ✓ | ✓ | ✓ | ✓ | ✓ |
| `Rows`, `Columns` | ✓ | ✓ | ✓ | ✓ | ✓ |
| `State` | ✓ | ✓ | ✓ | ✓ | ✓ |
| `Description` | ✓ | ✓ | ✓ | ✓ | ✓ |
| `River Basin Name(s)` | ✓ | ✓ | ✓ | ✓ | ✓ (dominant only) |
| `Data provider` | ✓ | ✓ | ✓ | ✓ | ✓ |
| `Map centroid coordinates (longitude, latitude)` | ✓ | ✓ | ✓ | ✓ | ✓ |
| `Image acquisition datetime (UTC)` | — | ✓ | ✓ | — | — |
| `Image acquisition date` | ✓ | ✓ | ✓ | — | — |
| `Image acquisition time (UTC)` | — | ✓ | ✓ | — | — |
| `Start Date of Flood`, `End Date of Flood` | — | — | — | ✓ | — |
| `Return period [yr]` | — | — | — | — | ✓ |
| `Event Recurrence Interval (USGS) [yr]` | ✓ | ✓ | ✓ | ✓ | — |
| `Recurrence Interval Discharge (USGS) [cms]` | ✓ | ✓ | ✓ | ✓ | — |
| `Reference USGS Streamgauge ID` | ✓ | ✓ | ✓ | ✓ | — |
| `Event Recurrence Interval (NWM) [yr]` | ✓ | ✓ | ✓ | ✓ | — |
| `Recurrence Interval Discharge (NWM) [cms]` | ✓ | ✓ | ✓ | ✓ | — |
| `Recurrence Interval Methodology Documentation` | ✓ | ✓ | ✓ | ✓ | — |
| `Keywords` | ✓ | ✓ | ✓ | ✓ | ✓ |
| `Access Rights` | ✓ | ✓ | ✓ | ✓ | ✓ |
| `Product Category` | ✓ | ✓ | ✓ | ✓ | ✓ |
| `Metadata Documentation` | ✓ | ✓ | ✓ | ✓ | ✓ |
| `References` | ✓ | ✓ | ✓ | ✓ | ✓ |

---

## 4. Field reference

### 4.1 Identification

#### File Name
The file name of the benchmark raster, for example
`S1A_10_0m_20170504T235530_891934W365313N_BM.tif`. It follows the
[naming convention](#naming-convention) and is the key used to find the
matching metadata, AOI and extent files.

#### Full sensor name
The data source written out in full, for example `"Aerial Imagery"` or
`"Sentinel-1A"`. See the table in [Section 2](#2-fim-categories-at-a-glance).
The web catalog copies it into `platform`.

#### Product Category
The product tier: `"Tier 1"`, `"Tier 2"`, `"Tier 3"`, `"High-water mark"` or
`"FEMA_BLE"`.

#### Data provider
Who produced the flood map, with an institution and contact email. See the
table in [Section 2](#2-fim-categories-at-a-glance).

#### Description
One or two sentences describing the map. It is filled in after the HUC8
lookup and follows a fixed template per category:

| Category | Template |
|---|---|
| Tier 1 | The Flood Inundation Map (FIM) was generated using NOAA Emergency Response Imagery acquired on `{date}` with a spatial resolution of `{res}`. The corresponding HUC IDs are `{HUC8}` |
| Tier 2 | Flood Inundation Map (FIM) generated using PlanetScope imagery acquired on `{YYYY-MM-DD}` at `{HH:MM:SS}` UTC, with a gap-filling algorithm applied to reconstruct inundation in unmapped areas. The map has a nominal spatial resolution of `{res}` m and covers portions of HUC8 watersheds `{a}, {b} and {c}`. |
| Tier 3 | As Tier 2, with "Sentinel-1 Imagery" and a resolution of `10` m |
| HWM | The Flood Inundation Map (FIM) was generated using USGS High Water Marks acquired for the `{start}-{end}` flood with a spatial resolution of 10 m. The corresponding HUC IDs are `{HUC8}` |
| FEMA BLE | The Flood Inundation Map (FIM) corresponds to FEMA Base Level Engineering for the `{return period}` flood with a spatial resolution of 10 m. The corresponding HUC IDs are `{HUC8}` |

For Tier 2 and Tier 3, "at … UTC" is left out when there is no acquisition
time, and a single HUC8 is written as "HUC8 watershed `{a}`". The Tier 2
processor always writes "PlanetScope imagery"; the RapidEye map in the database
(`REO_20161010T175532_775759W352159N`) says "RapidEye Ortho Tile imagery".

#### Keywords
`["flood inundation", "benchmark dataset", "geospatial data", …]`, followed by
the category's keyword(s):

| Category | Last keyword(s) |
|---|---|
| Tier 1 | `"aerial imagery"` |
| Tier 2 | `"planet"` |
| Tier 3 | `"sentinel-1"` |
| HWM | `"high-water mark"` |
| FEMA BLE | `"fema"`, `"base level engineering"` |

#### Access Rights
Always `"Public"`.

#### Recurrence Interval Methodology Documentation
Link to the [recurrence-interval documentation](recurrence-interval_documentation.md).
Not present for FEMA BLE.

#### Metadata Documentation
Link to this page.

#### References
Citations for the data and methods behind each category:

| Category | References |
|---|---|
| Tier 1 | Cohen et al. (2025), *Water Resources Research*, https://doi.org/10.1029/2024WR039574; Devi et al. (2025), *Environmental Modelling & Software*, https://doi.org/10.1016/j.envsoft.2025.106786 |
| Tier 2 | Tian et al. (2026), RS-FloodXDepth, *Water Resources Research*, https://doi.org/10.1029/2025WR042384; Cohen et al. (2025), https://doi.org/10.1029/2024WR039574 |
| Tier 3 | Tian et al. (2026), https://doi.org/10.1029/2025WR042384; Cohen et al. (2025), https://doi.org/10.1029/2024WR039574 |
| HWM | Baruah et al. (2026), *Natural Hazards and Earth System Sciences*, https://doi.org/10.5194/nhess-26-4753-2026 |
| FEMA BLE | Aristizabal et al. (2023), *Water Resources Research*, https://doi.org/10.1029/2022WR032039; Devi et al. (2025), https://doi.org/10.1016/j.envsoft.2025.106786 |

### 4.2 Raster properties

All raster properties describe the final `_BM.tif`, which is in EPSG:4326.

| Field | Meaning | Typical value |
|---|---|---|
| `Format` | GDAL driver of the raster | `"GTiff"` |
| `Nodata_Value` | Value used for pixels with no data | `-9999.0` |
| `Datatype` | Pixel data type. Integer sources are stored as `int32`, floating-point sources as `float32` | `"int32"` |
| `Compression type` | GeoTIFF compression and interleaving | `{"COMPRESSION": "LZW", "INTERLEAVE": "BAND"}` |
| `Coordinate reference system` | Coordinate reference system | `"EPSG:4326"` |
| `Number of bands` | Number of bands | `1` |
| `Rows`, `Columns` | Raster height and width in pixels | `10351`, `7736` |

#### Spatial resolution (m)
The pixel size in metres:

- **Tier 1 and Tier 2:** measured. The raster is reprojected to EPSG:5070
  (CONUS Albers, metres), the x and y pixel sizes are averaged, and the result
  is rounded to 0.1 m, for example `0.5` or `3.0`.
- **Tier 3, HWM and FEMA BLE:** always `10`, the nominal resolution of the
  product, whatever the measured pixel size.

The same value appears in the file name, for example `0_5m` or `10_0m`.

### 4.3 Location

#### Bounding box
The bounding box of the `_BM.tif`, in decimal degrees (EPSG:4326):

```json
"Bounding box": {"xmin": -89.684, "ymin": 36.408, "xmax": -88.968, "ymax": 37.366}
```

This box covers the whole raster, including any nodata margin. It is not the
box around the flooded pixels.

#### Map centroid coordinates (longitude, latitude)
`[longitude, latitude]` of the centre of the `Bounding box`, in decimal
degrees. It is the centre of the raster's bounding box, not the centre of the
flooded area.

#### Centroid DMS code
The centroid written as degrees-minutes-seconds without separators: longitude
first, then latitude, each with its hemisphere letter. For example,
`891934W365313N` means 89°19′34″ W, 36°53′13″ N. Values are truncated, not
rounded, to whole seconds. Longitudes of 100° or more give a longer code, for
example `1091442W451109N`.

The code identifies the map's location in its folder and file names.

#### HUC8 ID(s)
A list of the 8-digit USGS Watershed Boundary Dataset (WBD) HUC codes for the
watersheds the flooded area falls in, for example `["03020201", "03020202"]`.

They are found by querying the public ArcGIS REST service for WBD HUC8
boundaries with the flood polygon (every pixel with value > 0, from the
EPSG:5070 raster):

- **Tier 1, 2, 3 and HWM:** every HUC8 the flooded area intersects, sorted.
- **FEMA BLE:** only the single HUC8 that holds most of the flooded area,
  because BLE maps are produced per watershed. If no HUC8 holds at least 70 %
  of the area, the largest one is still used, and the shares are logged.

#### River Basin Name(s)
The names of the HUC8 watersheds in `HUC8 ID(s)`, sorted, for example
`["Middle Neuse", "Upper Neuse"]`.

#### State
The two-letter codes of the states the HUC8 watersheds cover, sorted and
comma-separated, for example `"AR, IL, KY, MO, MS, TN"`. Territories appear by
their code, for example `"PR"`. It is `null` when no HUC8 was found.

These are the states of the watersheds, so a map can list a state it doesn't
itself reach.

### 4.4 Event timing

The event fields depend on the category.

#### Tier 1: Image acquisition date
The acquisition date as given, usually `YYYYMMDD`, for example `"20110412"`.

#### Tier 2 and Tier 3: Image acquisition datetime (UTC), Image acquisition date, Image acquisition time (UTC)

| Field | Example |
|---|---|
| `Image acquisition datetime (UTC)` | `"20161009T150712"`, the date and time as given |
| `Image acquisition date` | `"20161009"`, the date part |
| `Image acquisition time (UTC)` | `"15:07:12"` (or `HH:MM` for a 4-digit time); `null` when no time was given |

For Tier 1–3, the date (and time) is passed to the processor or read from the
file name (`YYYYMMDD`, `YYYY-MM-DD` or `DD-MM-YYYY`, optionally followed by
`T<time>`). A location code after the date, such as `_953036W…`, is not taken
as a time. Times are in UTC.

#### HWM: Start Date of Flood, End Date of Flood
The flood window, each as `YYYYMMDD`, for example `"20100329"` and
`"20100401"`. High water marks record the highest water reached during an
event, so a window is used instead of a single date. The dates are read from
the file name (`…_YYYYMMDD_YYYYMMDD_…`), or passed to the processor.

#### FEMA BLE: Return period [yr]
The return period of the design flood in years, as a number, for example
`100` or `500`. BLE maps model a hypothetical flood of that size, not an
observed event, so they have no date.

### 4.5 Event magnitude (recurrence interval)

Observed-flood maps (Tier 1, 2, 3 and HWM) carry two independent estimates of
how rare the flood was. One uses observed USGS streamflow, the other the
National Water Model (NWM). The full procedure, including every scenario and
failure reason, is in
**[recurrence-interval_documentation.md](recurrence-interval_documentation.md)**.

| Field | Type | Values |
|---|---|---|
| `Event Recurrence Interval (USGS) [yr]` | number / string / null | Return period in years: `2`, `5`, `10`, `25`, `50`, `100`; `"< 2"` or `"> 100"` outside the published range; `"No USGS Gauge"`; `"No usable USGS Gauge"`; `null` if the lookup couldn't be attempted |
| `Recurrence Interval Discharge (USGS) [cms]` | number / null | Observed discharge at the gage in m³/s, 2 decimals; `null` unless a return period was assigned |
| `Reference USGS Streamgauge ID` | string / null | USGS site number of the gage used, for example `"02089000"` |
| `Event Recurrence Interval (NWM) [yr]` | number / string | A return period as above, or a short reason why none could be assigned, for example `"No outlet reach found"` |
| `Recurrence Interval Discharge (NWM) [cms]` | number / null | Modeled discharge at the flood's outlet reach in m³/s, 2 decimals |

The discharge behind each estimate depends on the category:

| Category | Discharge used |
|---|---|
| Tier 1 (date only) | Daily maximum of hourly discharge on the acquisition date |
| Tier 2, Tier 3 (date + time) | Discharge nearest the acquisition time |
| HWM | Maximum discharge over the flood window |

FEMA BLE maps have no recurrence-interval fields, because the return period is
already part of their definition.

---

## 5. Example metadata.json

### Tier 3 (Sentinel-1)

`Tier_3/S1A_20170504T235530_891934W365313N/S1A_10_0m_20170504T235530_891934W365313N_metadata.json`

```json
{
    "File Name": "S1A_10_0m_20170504T235530_891934W365313N_BM.tif",
    "Full sensor name": "Sentinel-1A",
    "HUC8 ID(s)": ["05140206", "07140105", "07140107", "07140108", "08010100",
                   "08010201", "08010202", "08010300", "08020204"],
    "Format": "GTiff",
    "Nodata_Value": -9999.0,
    "Spatial resolution (m)": 10,
    "Datatype": "int32",
    "Compression type": {"COMPRESSION": "LZW", "INTERLEAVE": "BAND"},
    "Bounding box": {"xmin": -89.68436, "ymin": 36.40796, "xmax": -88.96806, "ymax": 37.36639},
    "Centroid DMS code": "891934W365313N",
    "Coordinate reference system": "EPSG:4326",
    "Number of bands": 1,
    "Rows": 10351,
    "Columns": 7736,
    "State": "AR, IL, KY, MO, MS, TN",
    "Description": "Flood Inundation Map (FIM) generated using Sentinel-1 Imagery acquired on 2017-05-04 at 23:55:30 UTC, with a gap-filling algorithm applied to reconstruct inundation in unmapped areas. The map has a nominal spatial resolution of 10 m and covers portions of HUC8 watersheds 05140206, 07140105, 07140107, 07140108, 08010100, 08010201, 08010202, 08010300 and 08020204.",
    "River Basin Name(s)": ["Bayou De Chien-Mayfield", "Cache", "Little River Ditches",
                            "Lower Mississippi-Memphis", "Lower Ohio", "Obion",
                            "Town of Madrid-Saint Johns Bayou",
                            "Upper Mississippi-Cape Girardeau", "Whitewater"],
    "Data provider": "Dr. Dan Tian, The University of Alabama. Email ID: dtian1@ua.edu",
    "Map centroid coordinates (longitude, latitude)": [-89.32621, 36.88718],
    "Image acquisition datetime (UTC)": "20170504T235530",
    "Image acquisition date": "20170504",
    "Image acquisition time (UTC)": "23:55:30",
    "Event Recurrence Interval (USGS) [yr]": "No usable USGS Gauge",
    "Recurrence Interval Discharge (USGS) [cms]": null,
    "Reference USGS Streamgauge ID": null,
    "Event Recurrence Interval (NWM) [yr]": "< 2",
    "Recurrence Interval Discharge (NWM) [cms]": 27866.35,
    "Recurrence Interval Methodology Documentation": "https://github.com/sdmlua/fimbench/blob/main/docs/recurrence-interval_documentation.md",
    "Keywords": ["flood inundation", "benchmark dataset", "geospatial data", "sentinel-1"],
    "Access Rights": "Public",
    "Product Category": "Tier 3",
    "Metadata Documentation": "https://github.com/sdmlua/fimbench/blob/main/docs/metadata_documentation.md",
    "References": ["1. Tian, D., Liu, H., Wang, L., Cohen, S., & Mandal, T. (2026). ...", "2. ..."]
}
```

(Coordinates are shortened and the references abbreviated here.)

### How the other categories differ

**Tier 1** has a single date field:

```json
"Image acquisition date": "20110412",
"Event Recurrence Interval (USGS) [yr]": 50,
"Recurrence Interval Discharge (USGS) [cms]": 1718.83,
"Reference USGS Streamgauge ID": "05064500",
"Event Recurrence Interval (NWM) [yr]": "< 2",
"Recurrence Interval Discharge (NWM) [cms]": 362.74,
"Product Category": "Tier 1"
```

**Tier 2** has the same structure as Tier 3, with its measured resolution:

```json
"Spatial resolution (m)": 3.0,
"Image acquisition datetime (UTC)": "20161009T150712",
"Image acquisition date": "20161009",
"Image acquisition time (UTC)": "15:07:12",
"Event Recurrence Interval (USGS) [yr]": 5,
"Reference USGS Streamgauge ID": "02089000",
"Event Recurrence Interval (NWM) [yr]": "< 2",
"Product Category": "Tier 2"
```

**HWM** replaces the acquisition date with the flood window:

```json
"Spatial resolution (m)": 10,
"Start Date of Flood": "20100329",
"End Date of Flood": "20100401",
"Event Recurrence Interval (USGS) [yr]": "> 100",
"Event Recurrence Interval (NWM) [yr]": "> 100",
"Product Category": "High-water mark"
```

**FEMA BLE** replaces all event and magnitude fields with the return period:

```json
"Spatial resolution (m)": 10,
"Data provider": "NOAA/NWS Office of Water Prediction (OWP)",
"Return period [yr]": 100,
"Keywords": ["flood inundation", "benchmark dataset", "geospatial data", "fema", "base level engineering"],
"Product Category": "FEMA_BLE"
```

---

## 6. Companion files: AOI.gpkg and extent GeoJSON

### `{name}_AOI.gpkg`
The area of interest for the map: the footprint of the raster's valid data.

- **Layer:** `AOI`, one feature.
- **Geometry:** the area where the raster has data (not nodata or NaN),
  whether flooded or dry, in EPSG:4326. It is simplified with a tolerance of
  0.0001° (about 11 m).
  - For Tier 1 aerial maps (0.2–0.5 m pixels), the footprint is built on a
    coarse grid of about 5 m and read in chunks. This keeps very large rasters
    fast and the file small. Nodata gaps and specks smaller than the tolerance
    are dropped.
- **Attributes:** all fields of `metadata.json`.

It is used to check whether a model's prediction overlaps a benchmark map, for
example by the `query` module.

### `{name}_extent.geojson`
The flooded area with a selection of metadata. It is written when
`geojson_dest` is passed to `process()`, into that separate folder.

- **Geometry:** every flooded pixel (value > 0) of the EPSG:5070 raster,
  merged into one MultiPolygon. It is simplified by 10 m, which also drops
  patches smaller than that, then converted to EPSG:4326 with coordinates
  rounded to 6 decimals (about 0.1 m).
- **Properties:** a subset of the metadata fields (`GEOJSON_FIELDS` in each
  processor):

| Category | Properties |
|---|---|
| Tier 1 | `File Name`, `Full sensor name`, `HUC8 ID(s)`, `Spatial resolution (m)`, `Bounding box`, `State`, `River Basin Name(s)`, `Image acquisition date`, the five event-magnitude fields, `Product Category` |
| Tier 2 | As Tier 1, plus `Data provider`, and `Image acquisition datetime (UTC)` instead of `Image acquisition date` |
| Tier 3 | As Tier 2, without `Data provider` |
| HWM | As Tier 3, but `Start Date of Flood` and `End Date of Flood` instead of `Image acquisition datetime (UTC)` |
| FEMA BLE | `File Name`, `Full sensor name`, `HUC8 ID(s)`, `Spatial resolution (m)`, `Bounding box`, `State`, `River Basin Name(s)`, `Return period [yr]`, `Product Category` |

The extent GeoJSONs of all maps are merged and smoothed to build the web
map's `FIM_extents.geojson` and vector tiles (see the next section).

---

## 7. From metadata to the web catalog

`FIMCatalogBuilder` ([`webcontent_utils/build_catalog.py`](../src/fimbench/webcontent_utils/build_catalog.py))
reads every `*_metadata.json` under `FIM_Database/` in S3 and turns each one
into a record in `catalog_core.json`. The web map, the vector tiles and the
`query` module all use these records.

| Catalog field | Comes from |
|---|---|
| `id` | `{tier}/{folder}/{File Name without .tif}` |
| `site_id` | The map's folder name |
| `tier` | The S3 prefix: `Tier_1`, `Tier_2`, `Tier_3`, `HWM` or `FEMA_BLE` |
| `date_ymd` | The acquisition date as `YYYY-MM-DD` (Tier 1, 2, 3) |
| `event_ts` | `Image acquisition datetime (UTC)`, or `Image acquisition date` for Tier 1, as written, for example `20161009T150712` (Tier 1, 2, 3) |
| `start_date_ymd`, `end_date_ymd` | `Start Date of Flood`, `End Date of Flood` as `YYYY-MM-DD` (HWM) |
| `return_period` | `Return period [yr]` (FEMA BLE) |
| `json_url`, `tif_url`, `gpkg_url` | S3 URLs of the metadata, raster and AOI |
| `s3_prefix`, `s3_key` | The map's folder and metadata file in S3 |
| `centroid` | `Map centroid coordinates (longitude, latitude)` |
| `bbox` | `Bounding box` as `[xmin, ymin, xmax, ymax]` |
| `resolution_m` | `Spatial resolution (m)` |
| `huc8` | `HUC8 ID(s)` |
| `state`, `basin` | `State`, `River Basin Name(s)` |
| `source` | `Data provider` |
| `access_rights` | `Access Rights` |
| `description`, `references` | `Description`, `References` |
| `file_name` | `File Name` |
| `quality` | `Product Category` for Tier 1, 2, 3; `"High Water FIM"` for HWM; `"FEMA Base Level Engineering"` for FEMA BLE |
| `platform` | `Full sensor name`, for example `"Sentinel-1A"` |

The event-magnitude fields and the acquisition time are also copied, with
their values unchanged (including text such as `"No usable USGS Gauge"`), so
the web map's popup can show them without fetching each `metadata.json`. A key
is added only when the map's metadata has that field:

| Catalog field | Comes from | Present for |
|---|---|---|
| `recurrence_usgs` | `Event Recurrence Interval (USGS) [yr]` | Tier 1, 2, 3, HWM |
| `recurrence_nwm` | `Event Recurrence Interval (NWM) [yr]` | Tier 1, 2, 3, HWM |
| `discharge_usgs_cms` | `Recurrence Interval Discharge (USGS) [cms]` | Tier 1, 2, 3, HWM |
| `discharge_nwm_cms` | `Recurrence Interval Discharge (NWM) [cms]` | Tier 1, 2, 3, HWM |
| `usgs_gauge_id` | `Reference USGS Streamgauge ID` | Tier 1, 2, 3, HWM |
| `acquisition_time_utc` | `Image acquisition time (UTC)` | Tier 2, 3 |

`CatalogandTileManager` ([`webcontent_utils/tiling.py`](../src/fimbench/webcontent_utils/tiling.py))
then joins the catalog fields onto the merged extents to build
`FIM_extents.geojson` and the vector tiles in `FIM_Database/FIM_Viz/`.

---

## 8. Notes

- **`Bounding box` and the centroid describe the whole raster**, including
  nodata margins, not just the flooded area.
- **`Spatial resolution (m)`** is measured (after reprojection to EPSG:5070)
  for Tier 1 and Tier 2, and fixed at 10 m for Tier 3, HWM and FEMA BLE.
- **`State` comes from the HUC8 watersheds**, so it can include a state the
  flooded area doesn't reach.
- **Dates and times** are strings in `YYYYMMDD` or `YYYYMMDDTHHMMSS` form, and
  times are UTC. The catalog converts dates to `YYYY-MM-DD`.
- **`event_ts` in `FIM_extents.geojson` and the tiles** is converted to a
  number when the extents are built: a date-only value such as `20110412` is
  kept, while a value with a time such as `20161009T150712` becomes empty. The
  full value is in the catalog.
