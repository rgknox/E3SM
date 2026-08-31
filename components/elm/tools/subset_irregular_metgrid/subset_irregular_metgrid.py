# Regrid large DATM datasets for ELM runs that require a small number of sites
# or a regional domain. Supports nearest-neighbor selection by explicit
# coordinate list ("points" mode) or by reading coordinates from a target
# domain NetCDF file ("domain" mode).
#
# Usage:
#   python subset_irregular_metgrid.py --config my_config.json [--dry-run] [--yes]
#
# Questions: Ryan Knox rgknox@lbl.gov

import sys
import json
import argparse
import re
import numpy as np
import xarray as xr
import os
from pathlib import Path
import xml.etree.ElementTree as ET

_SCRIPT_DIR = Path(__file__).resolve().parent
_DEFAULT_STREAM_LIST_FILE = (
    _SCRIPT_DIR / "../../../data_comps/datm/cime_config/namelist_definition_datm.xml"
).resolve()

def GetXmlVals(root, group_str, tag_id_name, tag_name):
    entry = root.find(f"./entry[@id='{group_str}']")
    if entry is None:
        raise KeyError(f"No <entry id='{group_str}'> in stream list file")

    if tag_name == 'list':
        for val in entry.findall('.//value'):
            if val.attrib.get(tag_id_name):
                print(val.attrib[tag_id_name])
        return []

    match, fallback = None, None
    for val in entry.findall('.//value'):
        text = (val.text or '').strip()
        if not text:
            continue
        pattern = val.attrib.get(tag_id_name)
        if pattern is None:
            fallback = text
        elif re.search(pattern, tag_name):
            match = text          # last match wins, as CIME does

    result = match if match is not None else fallback
    if result is None:
        patterns = sorted({p for v in entry.findall('.//value')
                           if (p := v.attrib.get(tag_id_name))})
        raise KeyError(f"No value in '{group_str}' matches "
                       f"{tag_id_name}='{tag_name}'. Patterns: {patterns}")
    return [s.strip() for s in result.split(',')]


def GetXmlValsOld(root, group_str, tag_id_name, tag_name):

    text_list = []
    entry = f"./entry[@id='{group_str}']"
    streams_entry = root.find(entry)
    if streams_entry is None:
        print('big problems, no strealist?')
    else:
        value_tags = streams_entry.findall('.//value')
        for val in value_tags:
            mode_attr = val.attrib.get(tag_id_name, '')
            if tag_name == 'list':
                if mode_attr:
                    print(mode_attr)
            else:
                if mode_attr == tag_name and val.text:
                    text_list = [stream.strip() for stream in val.text.split(',')]

    return text_list


def RegridToPoints(ds: xr.Dataset, target_lats: list, target_lons: list):
    """
    Identify spatial indices in a dataset nearest to the provided lat-lon list.
    Works on any dataset with LATIXY and LONGXY coordinate variables.

    Returns integer vectors of (lat_indices, lon_indices) pairs.
    """
    latixy = ds["LATIXY"].values
    longxy = ds["LONGXY"].values

    domain_max_lon = np.nanmax(longxy)
    if domain_max_lon > 180.0:
        print(f"--> Detected Surface Convention: 0 to 360 (Max lon: {domain_max_lon:.2f})")
        target_lons_adjusted = np.array(target_lons) % 360
    else:
        print(f"--> Detected Surface Convention: -180 to 180 (Max lon: {domain_max_lon:.2f})")
        target_lons_adjusted = ((np.array(target_lons) + 180) % 360) - 180
    target_lats = np.array(target_lats)

    lat_indices = []
    lon_indices = []

    for tlat, tlon in zip(target_lats, target_lons_adjusted):
        dist = np.sqrt((latixy - tlat) ** 2 + (longxy - tlon) ** 2)
        idx = np.unravel_index(np.argmin(dist), dist.shape)
        lat_indices.append(idx[0])
        lon_indices.append(idx[1])

    print(f"base lon: {ds.LONGXY.isel(lat=lat_indices[0], lon=lon_indices[0]).values}")

    return lat_indices, lon_indices


def RegridSurface(base_ds: xr.Dataset,
                  target_lats: list,
                  target_lons: list,
                  lat_name: str = "lsmlat",
                  lon_name: str = "lsmlon") -> xr.Dataset:
    """
    Subset a surface dataset by nearest-neighbor matching to a list of
    (lat, lon) target points, collapsing the two spatial dimensions
    into a single 'gridcell' dimension (ELM vector layout).
    """
    lat_mesh = base_ds["LATIXY"].values
    lon_mesh = base_ds["LONGXY"].values

    print(f"base lon range:{np.min(lon_mesh)} {np.max(lon_mesh)}")

    domain_max_lon = np.nanmax(lon_mesh)
    if domain_max_lon > 180.0:
        print(f"--> Detected Surface Convention: 0 to 360 (Max lon: {domain_max_lon:.2f})")
        target_lons_adjusted = np.array(target_lons) % 360
    else:
        print(f"--> Detected Surface Convention: -180 to 180 (Max lon: {domain_max_lon:.2f})")
        target_lons_adjusted = ((np.array(target_lons) + 180) % 360) - 180

    target_lats = np.array(target_lats)

    matched_lat_indices = []
    matched_lon_indices = []

    for t_lat, t_lon in zip(target_lats, target_lons_adjusted):
        diff_squared = (lat_mesh - t_lat) ** 2 + (lon_mesh - t_lon) ** 2
        lat_idx, lon_idx = np.unravel_index(np.argmin(diff_squared), diff_squared.shape)
        print(f"pulling lat: {lat_mesh[lat_idx, lon_idx]:.4f}, "
              f"lon: {lon_mesh[lat_idx, lon_idx]:.4f}")
        matched_lat_indices.append(lat_idx)
        matched_lon_indices.append(lon_idx)

    lat_indices = np.array(matched_lat_indices)
    lon_indices = np.array(matched_lon_indices)

    lat_da = xr.DataArray(lat_indices, dims="gridcell")
    lon_da = xr.DataArray(lon_indices, dims="gridcell")

    new_ds = base_ds.isel({lat_name: lat_da, lon_name: lon_da})

    print(f"new surface lats: {new_ds['LATIXY'].values}, lons:{new_ds['LONGXY'].values}")

    return new_ds


def RegridMet(ds: xr.Dataset, lat_indices: list, lon_indices: list) -> xr.Dataset:
    """
    Subset a meteorological dataset at pre-identified grid point indices.

    Variables with "lat" and "lon" dimensions are extracted at the specified
    index pairs and stacked into a new dataset with lat=1 and lon=N_points.
    All other variables are passed through unchanged.
    """
    data_vars = {}

    for var in ds.data_vars:
        da = ds[var]

        if "lat" in da.dims and "lon" in da.dims:
            slices = [
                da.isel(lat=li, lon=loi).values
                for li, loi in zip(lat_indices, lon_indices)
            ]
            stacked = np.stack(slices, axis=-1)

            other_dims = [d for d in da.dims if d not in ("lat", "lon")]
            new_dims = other_dims + ["lat", "lon"]

            stacked = stacked[..., np.newaxis, :]

            data_vars[var] = xr.DataArray(stacked, dims=new_dims)

        elif "lat" in da.dims or "lon" in da.dims:
            data_vars[var] = da

        else:
            data_vars[var] = da

    new_coords = {}
    if "time" in ds.coords:
        new_coords["time"] = ds.coords["time"]

    return xr.Dataset(data_vars, coords=new_coords)


def RegridDomain(base_ds: xr.Dataset, lat_points: list, lon_points: list) -> xr.Dataset:
    """
    Subset a domain dataset (nj/ni conventions, xc/yc coordinate variables)
    by nearest-neighbor matching.
    """
    lat_mesh = base_ds["yc"].values
    lon_mesh = base_ds["xc"].values

    domain_max_lon = np.nanmax(lon_mesh)
    if domain_max_lon > 180.0:
        print(f"--> Detected Domain Convention: 0 to 360 (Max lon found: {domain_max_lon:.2f})")
        target_lons_adjusted = np.array(lon_points) % 360
    else:
        print(f"--> Detected Domain Convention: -180 to 180 (Max lon found: {domain_max_lon:.2f})")
        target_lons_adjusted = ((np.array(lon_points) + 180) % 360) - 180

    lat_points = np.array(lat_points)

    matched_lat_indices = []
    matched_lon_indices = []

    for t_lat, t_lon in zip(lat_points, target_lons_adjusted):
        diff_squared = (lat_mesh - t_lat) ** 2 + (lon_mesh - t_lon) ** 2
        lat_idx, lon_idx = np.unravel_index(np.argmin(diff_squared), diff_squared.shape)
        matched_lat_indices.append(lat_idx)
        matched_lon_indices.append(lon_idx)

    lat_indices = np.array(matched_lat_indices)
    lon_indices = np.array(matched_lon_indices)

    nj_indices_2d = lat_indices.reshape(1, -1)
    ni_indices_2d = lon_indices.reshape(1, -1)

    nj_da = xr.DataArray(nj_indices_2d, dims=("nj", "ni"))
    ni_da = xr.DataArray(ni_indices_2d, dims=("nj", "ni"))

    return base_ds.isel(nj=nj_da, ni=ni_da)


def GetCoordsFromDomain(domain_file: str):
    """
    Extract (lat, lon) coordinate pairs from a domain NetCDF file.

    Reads yc/xc coordinate arrays and filters to cells where mask == 1
    (if a mask variable is present).

    Returns (lat_list, lon_list) as Python lists.
    """
    ds = xr.open_dataset(domain_file, engine='netcdf4')
    yc = ds['yc'].values.flatten()
    xc = ds['xc'].values.flatten()
    mask = ds.get('mask')
    if mask is not None:
        valid = mask.values.flatten().astype(bool)
        yc = yc[valid]
        xc = xc[valid]
    ds.close()
    return yc.tolist(), xc.tolist()


def TrimDinLoc(text_str: str,
               din_loc_root: str,
               din_loc_root_clmforc: str,
               new_loc_root: str,
               new_loc_root_clmforc: str):
    """
    Split a path string at a DIN_LOC token and return the path suffix,
    the original base path, and the replacement base path.
    """
    if '$DIN_LOC_ROOT_CLMFORC' in text_str:
        path_end = text_str.split('$DIN_LOC_ROOT_CLMFORC')[1]
        path_base = din_loc_root_clmforc
        new_path_base = new_loc_root_clmforc
    elif '$DIN_LOC_ROOT' in text_str:
        path_end = text_str.split('$DIN_LOC_ROOT')[1]
        path_base = din_loc_root
        new_path_base = new_loc_root
    else:
        print('stream dir NOT found')
        sys.exit(2)

    return path_end, path_base, new_path_base


def validate_config(config: dict):
    required_always = [
        "din_loc_root", "din_loc_root_clmforc",
        "new_loc_root", "new_loc_root_clmforc",
        "target_mode", "target_coords_mode",
    ]
    for key in required_always:
        if key not in config:
            raise ValueError(f"Config missing required field: '{key}'")

    if not config.get("base_lnd_domain"):
        raise ValueError("'base_lnd_domain' is required in all modes")

    mode = config["target_coords_mode"]
    if mode == "points":
        if not config.get("lat_list") or not config.get("lon_list"):
            raise ValueError("'lat_list' and 'lon_list' are required when target_coords_mode is 'points'")
        if len(config["lat_list"]) != len(config["lon_list"]):
            raise ValueError("'lat_list' and 'lon_list' must have the same length")
    elif mode == "domain":
        if not config.get("target_domain_file"):
            raise ValueError("'target_domain_file' is required when target_coords_mode is 'domain'")
    else:
        raise ValueError(f"Unknown target_coords_mode: '{mode}'. Must be 'points' or 'domain'")


def _confirm(prompt: str, yes: bool) -> bool:
    if yes:
        print(f"{prompt} [auto-yes]")
        return True
    response = input(f"{prompt} (Y/N): ").strip().lower()
    return response == "y"


def main():
    parser = argparse.ArgumentParser(
        description="Regrid DATM meteorological data to a subset of grid points for ELM runs."
    )
    parser.add_argument(
        "--config", required=True,
        help="Path to JSON configuration file"
    )
    parser.add_argument(
        "--dry-run", action="store_true", default=False,
        help="Print what would be done without writing any output files"
    )
    parser.add_argument(
        "--yes", "-y", action="store_true", default=False,
        help="Skip all interactive confirmation prompts"
    )
    args = parser.parse_args()

    with open(args.config, "r") as f:
        config = json.load(f)

    validate_config(config)

    dry_run = args.dry_run or config.get("dry_run", False)
    yes = args.yes

    din_loc_root           = config["din_loc_root"]
    din_loc_root_clmforc   = config["din_loc_root_clmforc"]
    new_loc_root           = config["new_loc_root"]
    new_loc_root_clmforc   = config["new_loc_root_clmforc"]
    target_mode            = config["target_mode"]
    target_coords_mode     = config["target_coords_mode"]
    base_surf_files        = config.get("base_surf_files", [])

    stream_list_file = config.get("stream_list_file") or str(_DEFAULT_STREAM_LIST_FILE)
    print(f"Using stream list file: {stream_list_file}")

    # --- Determine target coordinates ---
    if target_coords_mode == "points":
        lat_list = config["lat_list"]
        lon_list = config["lon_list"]
        print(f"Mode: points ({len(lat_list)} point(s))")
    else:
        print(f"Mode: domain — reading coordinates from {config['target_domain_file']}")
        lat_list, lon_list = GetCoordsFromDomain(config["target_domain_file"])
        print(f"Loaded {len(lat_list)} target point(s) from domain file")

    # --- Parse stream XML ---
    tree = ET.parse(stream_list_file)
    stream_root = tree.getroot()

    stream_list = GetXmlVals(stream_root, "streamslist", "datm_mode", target_mode)
    print(f"Streams for mode '{target_mode}': {stream_list}")

    domain_dir  = GetXmlVals(stream_root, "strm_domdir", "stream", stream_list[0])
    domain_file = GetXmlVals(stream_root, "strm_domfil", "stream", stream_list[0])
    
    domain_end, domain_base, new_domain_base = TrimDinLoc(
        domain_dir[0] + '/' + domain_file[0],
        din_loc_root, din_loc_root_clmforc, new_loc_root, new_loc_root_clmforc
    )
    new_domain_file = os.path.join(new_domain_base, domain_end.lstrip('/'))
    new_domain_path = Path(new_domain_file)
    new_domain_dir  = str(new_domain_path.parent)

    print(f"\nWill create directory: '{new_domain_dir}'")
    print(f"  for modified domain file: {new_domain_path.name}")
    print(f"  full path: {new_domain_file}")
    if not _confirm("Continue?", yes):
        print("Exiting.")
        sys.exit()

    if not dry_run:
        new_domain_path.parent.mkdir(parents=True, exist_ok=True)

    # --- Create output domain file ---
    # Always subset base_lnd_domain to the target points so the output domain
    # has the same nj=1/ni=N vector layout as the output met data.
    base_lnd_domain = config["base_lnd_domain"]
    print(f"\nWill subset domain: {base_lnd_domain}")
    print(f"  -> {new_domain_file}")
    if not _confirm("Continue?", yes):
        print("Exiting.")
        sys.exit()
    if not dry_run:
        base_ds = xr.open_dataset(base_lnd_domain, engine='netcdf4')
        new_ds = RegridDomain(base_ds, lat_list, lon_list)
        new_ds.to_netcdf(new_domain_file)
        base_ds.close()
        new_ds.close()

    # --- Convert surface files ---
    for base_surf_file in base_surf_files:
        base_surf_path = Path(base_surf_file)
        new_surf_file = new_domain_dir + '/' + base_surf_path.name

        print(f"\nWill subset surface: {base_surf_file}")
        print(f"  -> {new_surf_file}")
        if not _confirm("Continue?", yes):
            print("Exiting.")
            sys.exit()
        if not dry_run:
            base_ds = xr.open_dataset(base_surf_file, engine='netcdf4')
            new_ds = RegridSurface(base_ds, lat_list, lon_list)
            new_ds.to_netcdf(new_surf_file)

    # --- Regrid met data streams ---
    stream_list = GetXmlVals(stream_root, "streamslist", "datm_mode", target_mode)

    print(f"Streams for mode '{target_mode}': {stream_list}")

    domain_dir  = GetXmlVals(stream_root, "strm_domdir", "stream", stream_list[0])
    domain_file = GetXmlVals(stream_root, "strm_domfil", "stream", stream_list[0])
    
    for stream in stream_list:
        stream_dir = GetXmlVals(stream_root, "strm_datdir", "stream", stream)
        stream_path_end, stream_path_base, new_path_base = TrimDinLoc(
            stream_dir[0],
            din_loc_root, din_loc_root_clmforc, new_loc_root, new_loc_root_clmforc
        )

        stream_path_dir = os.path.join(stream_path_base, stream_path_end.lstrip('/'))
        new_path_dir    = os.path.join(new_path_base, stream_path_end.lstrip('/'))
        new_path        = Path(new_path_dir)

        if os.path.isdir(new_path_dir):
            print(f"Folder exists: {new_path_dir}")
        else:
            print(f"{'[dry-run] Would create' if dry_run else 'Creating'} folder: {new_path_dir}")
            if not dry_run:
                new_path.mkdir(parents=True)

        if not os.path.isdir(stream_path_dir):
            if dry_run:
                print(f"  [dry-run] Source stream path not found (expected on non-NERSC): {stream_path_dir}")
                continue
            print(f"Critical error, cannot find stream path: {stream_path_dir}")
            sys.exit(2)

        if dry_run:
            nc_files = sorted(Path(stream_path_dir).glob("clmforc*.nc"))
            print(f"  Would convert {len(nc_files)} file(s) from {stream_path_dir}")
            for fp in nc_files:
                print(f"    {fp.name} -> {new_path_dir}/{fp.name}")
            continue

        lat_indices = None
        lon_indices = None


        failed = []
        for ifp, file_path in enumerate(sorted(Path(stream_path_dir).glob("clmforc*.nc"))):
            base_file_path = str(file_path)
            new_file_path = new_path_dir + '/' + file_path.name

            try:
                base_ds = xr.open_dataset(base_file_path, engine='netcdf4')
            except OSError as err:
                print(f"  SKIPPING unreadable file {file_path.name}: {err}")
                failed.append(file_path.name)
                continue

            with base_ds:
                if lat_indices is None:
                    lat_indices, lon_indices = RegridToPoints(base_ds, lat_list, lon_list)
                    print(f"Setting up grid translation")
                print(f"Converting: {base_file_path} -> {new_file_path}")
                new_ds = RegridMet(base_ds, lat_indices, lon_indices)
                new_ds.to_netcdf(new_file_path)
                new_ds.close()

        if failed:
            print(f"\n{len(failed)} file(s) could not be read: {failed}")

    print("\nDone.")


if __name__ == "__main__":
    main()
