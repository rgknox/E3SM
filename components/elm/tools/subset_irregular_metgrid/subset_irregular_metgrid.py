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
import shutil
import os
from pathlib import Path
import xml.etree.ElementTree as ET

_SCRIPT_DIR = Path(__file__).resolve().parent
_DEFAULT_STREAM_LIST_FILE = (
    _SCRIPT_DIR / "../../../data_comps/datm/cime_config/namelist_definition_datm.xml"
).resolve()


def _is_good(path):
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return False
    try:
        with xr.open_dataset(path, engine='netcdf4'):
            return True
    except OSError:
        return False

    
def PrintCaseSetup(config, new_stream_domain, new_lnd_domain, stream_dirs,stream_list, new_surf_files):
    """
    Print a copy-pasteable block of case-setup commands for the ELM run.

    stream_dirs   : dict of {stream_name: new forcing directory}
    new_surf_files: list of subsetted surface dataset paths
    """
    bar = "=" * 78

    print(f"\n{bar}")
    print("THE FOLLOWING IS FOR YOUR ELM CASE SETUP SCRIPT")
    print("Paste into the bash script you use to configure this case, after")
    print("create_newcase and before ./case.build.")
    print(f"{bar}\n")

    print("# ---- domains (model domain: stock masked land domain) ----")
    print(f"./xmlchange ATM_DOMAIN_PATH={new_lnd_domain.parent},ATM_DOMAIN_FILE={new_lnd_domain.name}")
    print(f"./xmlchange LND_DOMAIN_PATH={new_lnd_domain.parent},LND_DOMAIN_FILE={new_lnd_domain.name}\n")

    print(f'./xmlchange DIN_LOC_ROOT_CLMFORC={config["new_loc_root_clmforc"]}\n')
    
    print("# ---- generate namelists so the streams files exist ----")
    print("./preview_namelists\n")


    print("RUNDIR=$(./xmlquery -value RUNDIR)\n")

    print("# ---- redirect the stream domain (DIN_LOC_ROOT_CLMFORC does not")
    print("#      reach it: strm_domdir is often rooted at DIN_LOC_ROOT) ----")
    print(f'MYDOM_DIR={new_stream_domain.parent}\n')

    print("streams=(" + " ".join(stream_list) + ")")
    print("flddirs=(" + " ".join(stream_dirs[s] for s in stream_list) + ")\n")

    print('for i in "${!streams[@]}"; do')
    print("  s=${streams[$i]}; fld=${flddirs[$i]}")
    print("  f=$RUNDIR/datm.streams.txt.$s")
    print('  [ -f "$f" ] || { echo "MISSING: $f"; exit 1; }')
    print('  awk -v dom="$MYDOM_DIR" -v fld="$fld" \'')
    print('    /<filePath>/ { n++; print; print "     " (n==1 ? dom : fld); skip=1; next }')
    print('    skip && /<\\/filePath>/ { skip=0 }')
    print('    skip { next }')
    print('    { print }')
    print("  ' \"$f\" > ./user_datm.streams.txt.$s")
    print('  np=$(grep -c "<filePath>" ./user_datm.streams.txt.$s)')
    print('  [ "$np" = "2" ] || { echo "WARNING: $s has $np filePath entries"; }')
    print('  echo "wrote user_datm.streams.txt.$s"')
    print("done\n")

    print("./preview_namelists")
    print('for s in "${streams[@]}"; do')
    print('  echo "--- $s"; grep -A1 "<filePath>" $RUNDIR/datm.streams.txt.$s')
    print("done")

    print("# ---- mapalgo: 'copy' for subsetted streams, 'nn' for the rest ----")
    print("# Order follows the 'streams' array in datm_in")
    print("grep -A8 '^ *streams' $RUNDIR/datm_in   # verify before trusting this")
    print("cat >> user_nl_datm <<'EOF'")
    print(f"mapalgo = {','.join(chr(39)+'copy'+chr(39) for _ in stream_dirs)},'nn','nn'")
    print("EOF\n")

    if new_surf_files:
        print("# ---- surface dataset ----")
        print("cat >> user_nl_clm <<EOF")
        for s in new_surf_files:
            print(f"fsurdat = '{s}'")
        print("EOF\n")


    print(f"\n{bar}")
    print("END OF CASE SETUP BLOCK")
    print(f"{bar}\n")

    
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

def _unit_vectors(lat_deg, lon_deg):
    """Convert degrees to unit vectors on the sphere (broadcasts over any shape)."""
    la = np.deg2rad(np.asarray(lat_deg, dtype=np.float64))
    lo = np.deg2rad(np.asarray(lon_deg, dtype=np.float64))
    return np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)


def NearestGridIndex(mesh_lat, mesh_lon, t_lat, t_lon):
    """
    Index of the mesh point nearest to (t_lat, t_lon) by great-circle distance.

    Returns (index_tuple, distance_km). Minimizing chord length is equivalent
    to minimizing arc length, so no arccos is needed for the search itself.
    """
    mx, my, mz = _unit_vectors(mesh_lat, mesh_lon)
    tx, ty, tz = _unit_vectors(t_lat, t_lon)
    chord2 = (mx - tx) ** 2 + (my - ty) ** 2 + (mz - tz) ** 2
    idx = np.unravel_index(np.nanargmin(chord2), chord2.shape)
    dist_km = 2.0 * 6371.0 * np.arcsin(np.clip(np.sqrt(chord2[idx]) / 2.0, 0.0, 1.0))
    return idx, dist_km


def RegridToPoints(ds: xr.Dataset, target_lats: list, target_lons: list):
    """
    Identify spatial indices in a dataset nearest to the provided lat-lon list.
    Works on any dataset with LATIXY and LONGXY coordinate variables.

    Returns integer vectors of (lat_indices, lon_indices) pairs.
    """
    latixy = ds["LATIXY"].values
    longxy = ds["LONGXY"].values

    lat_indices, lon_indices, distances = [], [], []
    for tlat, tlon in zip(target_lats, target_lons):
        idx, dkm = NearestGridIndex(latixy, longxy, tlat, tlon)
        lat_indices.append(idx[0])
        lon_indices.append(idx[1])
        distances.append(dkm)

    distances = np.array(distances)
    print(f"--> matched {len(distances)} point(s); "
          f"nearest-neighbor distance min/mean/max = "
          f"{distances.min():.1f}/{distances.mean():.1f}/{distances.max():.1f} km")
    worst = np.argsort(distances)[-5:][::-1]
    for w in worst:
        print(f"    target ({target_lats[w]:.2f}, {target_lons[w]:.2f}) -> "
              f"({latixy[lat_indices[w], lon_indices[w]]:.2f}, "
              f"{longxy[lat_indices[w], lon_indices[w]]:.2f})  {distances[w]:.1f} km")

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
    lat_mesh = np.atleast_2d(base_ds["LATIXY"].values)
    lon_mesh = np.atleast_2d(base_ds["LONGXY"].values)

    matched_lat_indices, matched_lon_indices, distances = [], [], []
    for t_lat, t_lon in zip(target_lats, target_lons):
        idx, dkm = NearestGridIndex(lat_mesh, lon_mesh, t_lat, t_lon)
        matched_lat_indices.append(idx[0])
        matched_lon_indices.append(idx[1])
        distances.append(dkm)

    distances = np.array(distances)
    print(f"--> surface: matched {len(distances)} point(s); "
          f"distance min/mean/max = {distances.min():.1f}/"
          f"{distances.mean():.1f}/{distances.max():.1f} km")

    lat_da = xr.DataArray(np.array(matched_lat_indices), dims="gridcell")
    lon_da = xr.DataArray(np.array(matched_lon_indices), dims="gridcell")

    new_ds = base_ds.isel({lat_name: lat_da, lon_name: lon_da})
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

def CheckDomainVars(ds: xr.Dataset, filepath: str):

    DOMAIN_VARS_REQUIRED = ("xc", "yc", "xv", "yv", "mask", "area", "frac")
    missing = [v for v in DOMAIN_VARS_REQUIRED if v not in ds.variables]
    if missing:
        raise ValueError(
            f"Domain File:{filepath}: missing {missing}. This workflow requires a land domain "
            f"file (share/domains/domain.lnd.*.nc). Stream domains such as "
            f"domain.T62.050609.nc lack frac and are not supported as input."
        )
    
def RegridDomain(base_ds: xr.Dataset, lat_points: list, lon_points: list, filepath: str) -> xr.Dataset:
    """
    Subset a domain dataset (nj/ni conventions, xc/yc coordinate variables)
    by nearest-neighbor matching, returning an nj=1, ni=N vector layout.
    Note that we only want to work with Land domains, that have the frac variable.
    """
    CheckDomainVars(base_ds,filepath)
        
    lat_mesh = np.atleast_2d(base_ds["yc"].values)
    lon_mesh = np.atleast_2d(base_ds["xc"].values)

    matched_lat_indices, matched_lon_indices, distances = [], [], []
    for t_lat, t_lon in zip(lat_points, lon_points):
        idx, dkm = NearestGridIndex(lat_mesh, lon_mesh, t_lat, t_lon)
        matched_lat_indices.append(idx[0])
        matched_lon_indices.append(idx[1])
        distances.append(dkm)

    distances = np.array(distances)
    print(f"--> domain: matched {len(distances)} point(s); "
          f"distance min/mean/max = {distances.min():.1f}/"
          f"{distances.mean():.1f}/{distances.max():.1f} km")

    nj_da = xr.DataArray(np.array(matched_lat_indices).reshape(1, -1), dims=("nj", "ni"))
    ni_da = xr.DataArray(np.array(matched_lon_indices).reshape(1, -1), dims=("nj", "ni"))

    return base_ds.isel(nj=nj_da, ni=ni_da)

def GetCoordsFromDomain(domain_file: str, apply_mask: bool = False):
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
    if apply_mask and 'mask' in ds:
        valid = ds['mask'].values.flatten().astype(bool)
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
        "datm_mode_basis", "coords_mode",
    ]
    for key in required_always:
        if key not in config:
            raise ValueError(f"Config missing required field: '{key}'")

    mode = config["coords_mode"]
    if mode == "points":
        if not config.get("lat_list") or not config.get("lon_list"):
            raise ValueError("'lat_list' and 'lon_list' are required in 'points' mode")
        if len(config["lat_list"]) != len(config["lon_list"]):
            raise ValueError("'lat_list' and 'lon_list' must have the same length")
    elif mode == "domain":
        if not config.get("target_lnd_domain"):
            raise ValueError("'target_lnd_domain' is required in 'domain' mode")

    else:
        raise ValueError(f"Unknown coords_mode: '{mode}'")
    

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
    parser.add_argument(
        "--skip-existing", action="store_true", default=False,
        help="Skip regridding any output file that already exists and is non-empty"
    )
    args = parser.parse_args()

    with open(args.config, "r") as f:
        config = json.load(f)

    validate_config(config)

    dry_run = args.dry_run or config.get("dry_run", False)
    skip_existing = args.skip_existing or config.get("skip_existing", False)
    yes = args.yes

    din_loc_root           = config["din_loc_root"]
    din_loc_root_clmforc   = config["din_loc_root_clmforc"]
    new_loc_root           = config["new_loc_root"]
    new_loc_root_clmforc   = config["new_loc_root_clmforc"]
    target_lnd_domain      = config["target_lnd_domain"]
    datm_mode_basis        = config["datm_mode_basis"]
    coords_mode            = config["coords_mode"]
    base_surf_files        = config.get("base_surf_files", [])
    stream_list_file       = config.get("stream_list_file") or str(_DEFAULT_STREAM_LIST_FILE)
    
    print(f"Using stream list file: {stream_list_file}")

    # --- Determine target coordinates ---
    if coords_mode == "points":
        lat_list = config["lat_list"]
        lon_list = config["lon_list"]
        print(f"Mode: points ({len(lat_list)} point(s))")
    else:
        print(f"Mode: domain — reading coordinates from {config['target_lnd_domain']}")
        lat_list, lon_list = GetCoordsFromDomain(config["target_lnd_domain"])
        print(f"Loaded {len(lat_list)} target point(s) from domain file")

    # --- Parse stream XML ---
    tree = ET.parse(stream_list_file)
    stream_root = tree.getroot()

    stream_list = GetXmlVals(stream_root, "streamslist", "datm_mode", datm_mode_basis)
    print(f"Streams for mode '{datm_mode_basis}': {stream_list}")


    # Create the new domain file for the new stream files
    stream_domain_dir  = GetXmlVals(stream_root, "strm_domdir", "stream", stream_list[0])
    stream_domain_file = GetXmlVals(stream_root, "strm_domfil", "stream", stream_list[0])
    
    stream_domain_end, stream_domain_base, new_domain_base = TrimDinLoc(
        stream_domain_dir[0] + '/' + stream_domain_file[0],
        din_loc_root, din_loc_root_clmforc, new_loc_root, new_loc_root_clmforc
    )
    
    new_lnd_domain_name = Path(target_lnd_domain).name.removesuffix(".nc")+"_lnd.nc"
    new_stream_domain = Path(os.path.join(new_domain_base,stream_domain_file[0]))
    new_lnd_domain = Path(os.path.join(new_domain_base,new_lnd_domain_name))
    
    #p.name      # 'domain.ne4pg2.datm.nc'
    #p.parent    # PosixPath('/global/cfs/cdirs/m2420/rgknox')

    print(f"\nWill create directory: '{new_domain_base}'")
    print(f"  for new stream domain file: {new_stream_domain.name}")
    print(f"  and new land domain file: {new_lnd_domain.name}")
    if not _confirm("Continue?", yes):
        print("Exiting.")
        sys.exit()

    if not dry_run:
        new_stream_domain.parent.mkdir(parents=True, exist_ok=True)

    # --- Create the new domain files ---
    # Always subset base_lnd_domain to the target points so the output domain
    # has the same nj=1/ni=N vector layout as the output met data.

    if coords_mode == "points":
        print(f"\nWill subset domain: {target_lnd_domain}")
        print(f"  -> {new_domain_file}")
        if not _confirm("Continue?", yes):
            print("Exiting.")
            sys.exit()
        if not dry_run:
            with xr.open_dataset(target_lnd_domain, engine='netcdf4') as base_ds:
                RegridDomain(base_ds, lat_list, lon_list, target_lnd_domain).to_netcdf(new_lnd_domain)
    else:
        # If in domain mode, the new land domain is an exact copy of the target
        # And the stream domain is an unmasked copy
        shutil.copyfile(target_lnd_domain, new_lnd_domain)


    # In both cases, create the stream domain file
    ds = xr.open_dataset(new_lnd_domain)
    CheckDomainVars(ds,new_lnd_domain)
    ds["mask"] = xr.ones_like(ds["mask"])
    ds.to_netcdf(new_stream_domain)

        

    # --- Convert surface files ---
    new_surf_files = []
    for base_surf_file in base_surf_files:
        base_surf_path = Path(base_surf_file)
        new_surf_file = new_domain_dir + '/' + base_surf_path.name
        new_surf_files.append(new_surf_file)
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
    stream_dirs = {}
    for stream in stream_list:
        stream_dir = GetXmlVals(stream_root, "strm_datdir", "stream", stream)
        stream_path_end, stream_path_base, new_path_base = TrimDinLoc(
            stream_dir[0],
            din_loc_root, din_loc_root_clmforc, new_loc_root, new_loc_root_clmforc
        )
        stream_path_dir = os.path.join(stream_path_base, stream_path_end.lstrip('/'))
        new_path_dir    = os.path.join(new_path_base, stream_path_end.lstrip('/'))
        stream_dirs[stream] = new_path_dir
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

        n_skipped = 0
        failed = []
        for ifp, file_path in enumerate(sorted(Path(stream_path_dir).glob("clmforc*.nc"))):
            base_file_path = str(file_path)
            new_file_path = new_path_dir + '/' + file_path.name

            if skip_existing and _is_good(new_file_path):
                print(f"  exists, skipping: {file_path.name}")
                n_skipped += 1
                continue
            
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

    PrintCaseSetup(config, new_stream_domain, new_lnd_domain, stream_dirs,stream_list, new_surf_files)
            
    print("\nDone.")


if __name__ == "__main__":
    main()
