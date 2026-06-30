# This notebook is intended to re-grid large DATM 
# datasets for use in ELM runs that either only require
# a small number of sites, or are production runs.
# The intent here is to speed up the model
# Questions: Ryan Knox rgknox@lbl.gov

import sys
import json
import runpy
import numpy as np
import xarray as xr
import sys
import os
import warnings
import matplotlib.pyplot as plt # Import matplotlib for inline plotting
from pathlib import Path
import xml.etree.ElementTree as ET

# Here are the components as listed in a datm_in from a pre-industrial run:
#streams = "datm.streams.txt.CLMGSWP3v1.Solar 1 1948 1972",
#      "datm.streams.txt.CLMGSWP3v1.Precip 1 1948 1972",
#      "datm.streams.txt.CLMGSWP3v1.TPQW 1 1948 1972",
#      "datm.streams.txt.presaero.clim_1850 1 1 1",
#      "datm.streams.txt.topo.observed 1 1 1"


def RegridToPoints(ds: xr.Dataset, target_lats: list, target_lons: list):
    """
    Identify spatial indices in a dataset that are nearest to the lat-lon list
    provided.  This can work on any dataset that has geographic coordinates
    specified by LATIXY and LONGXY.
    
    Parameters
    ----------
    ds : xr.Dataset
        Input dataset with LATIXY and LONGXY variables defining grid coordinates.
    target_lats : list
        List of target latitude coordinates.
    target_lons : list
        List of target longitude coordinates.
    
    Returns integer vectors of lat/lon indices pairs
    -------

    """
    # Extract the 2D coordinate arrays
    latixy = ds["LATIXY"].values  # shape (lat, lon)
    longxy = ds["LONGXY"].values  # shape (lat, lon)

    lat_indices = []
    lon_indices = []

    for tlat, tlon in zip(target_lats, target_lons):
        # Compute distance (simple Euclidean in lat/lon space)
        # For more accuracy over large distances, use haversine instead
        dist = np.sqrt((latixy - tlat) ** 2 + (longxy - tlon) ** 2)
        idx = np.unravel_index(np.argmin(dist), dist.shape)
        lat_indices.append(idx[0])
        lon_indices.append(idx[1])

    return lat_indices,lon_indices

def RegridSurface(base_ds: xr.Dataset,
                  target_lats: list,
                  lon_points: list,
                  lat_name: str = "lsmlat",
                  lon_name: str = "lsmlon") -> xr.Dataset:
    """
    Subset a surface dataset by nearest-neighbor matching to a list of
    (lat, lon) target points, collapsing the two spatial dimensions
    (lsmlat, lsmlon) into a single 'gridcell' dimension of length N
    (ELM vector layout).

    Coordinate fields are LATIXY (degrees north) and LONGXY
    (degrees east), both shaped (lsmlat, lsmlon). lsmlat/lsmlon are the
    trailing dimensions on all spatial variables, so higher-dimensional
    variables (e.g. CV_IMPROAD(nlevurb, lsmlat, lsmlon)) keep their
    leading dims and gain a single trailing 'gridcell' dim.
    """

    # --- coordinate meshes --------------------------------------------
    lat_mesh = base_ds["LATIXY"].values   # Shape: (lsmlat, lsmlon)
    lon_mesh = base_ds["LONGXY"].values   # Shape: (lsmlat, lsmlon)

    # --- normalize longitude convention -------------------------------
    domain_max_lon = np.nanmax(lon_mesh)
    if domain_max_lon > 180.0:
        print(f"--> Detected Surface Convention: 0 to 360 (Max lon: {domain_max_lon:.2f})")
        target_lons_adjusted = np.array(lon_points) % 360
    else:
        print(f"--> Detected Surface Convention: -180 to 180 (Max lon: {domain_max_lon:.2f})")
        target_lons_adjusted = ((np.array(lon_points) + 180) % 360) - 180

    target_lats = np.array(target_lats)

    # --- nearest-neighbor index search --------------------------------
    matched_lat_indices = []
    matched_lon_indices = []

    for t_lat, t_lon in zip(lat_points, target_lons_adjusted):
        diff_squared = (lat_mesh - t_lat) ** 2 + (lon_mesh - t_lon) ** 2
        lat_idx, lon_idx = np.unravel_index(np.argmin(diff_squared),
                                            diff_squared.shape)
        print(f"pulling lat: {lat_mesh[lat_idx, lon_idx]:.4f}, "
              f"lon: {lon_mesh[lat_idx, lon_idx]:.4f}")
        matched_lat_indices.append(lat_idx)
        matched_lon_indices.append(lon_idx)

    lat_indices = np.array(matched_lat_indices)
    lon_indices = np.array(matched_lon_indices)

    # --- pointwise selection -> single 'gridcell' dim -----------------
    # Sharing one dim name on both index arrays triggers vectorized
    # (pointwise) selection, collapsing lsmlat/lsmlon to 'gridcell'.
    lat_da = xr.DataArray(lat_indices, dims="gridcell")
    lon_da = xr.DataArray(lon_indices, dims="gridcell")

    new_ds = base_ds.isel({lat_name: lat_da, lon_name: lon_da})

    return new_ds

    
    

 def RegridMet(ds: xr.Dataset, lat_indices: list, lon_indices: list) -> xr.Dataset:
    """
    Subset a dataset. This assumes that the desired points have already
    been identified. Data arrays with dimensions "lat" and "lon" will be subset using the
    index lists provided, all other data arrays will simply be copies
    of the original data array provided. Note that the new lat dim will
    be of LENGTH 1, and the lon dim will be: len(target_lats).
    
    Parameters
    ----------
    ds : xr.Dataset
        Input dataset with LATIXY and LONGXY variables defining grid coordinates.
    target_lats : list
        List of target latitude coordinates.
    target_lons : list
        List of target longitude coordinates.
    
    Returns
    -------
    xr.Dataset
        New dataset with lat dim = 1 and lon dim = len(target_lats),
        containing only the nearest neighbor grid points.
    """
    
    # Extract data at the selected grid points
    data_vars = {}

    for var in ds.data_vars:
        da = ds[var]

        if "lat" in da.dims and "lon" in da.dims:
            # Gather slices at each selected point
            slices = [
                da.isel(lat=li, lon=loi).values
                for li, loi in zip(lat_indices, lon_indices)
            ]
            # Stack along a new axis to form the new "lon" dimension
            # Each slice has shape (time,) or () depending on the variable
            stacked = np.stack(slices, axis=-1)  # shape: (time, n_points) or (n_points,)

            # Build new dims: replace lat/lon with lat=1, lon=n_points
            other_dims = [d for d in da.dims if d not in ("lat", "lon")]
            new_dims = other_dims + ["lat", "lon"]

            # Add the lat=1 dimension
            stacked = stacked[..., np.newaxis, :]  # shape: (..., 1, n_points)

            data_vars[var] = xr.DataArray(stacked, dims=new_dims)

        elif "lat" in da.dims or "lon" in da.dims:
            # Edge case: variable has only one of lat or lon — just pass through
            data_vars[var] = da

        else:
            # No spatial dimensions (e.g., scalar or time-only vars)
            data_vars[var] = da

    # Build coordinates for the new dataset
    new_coords = {}
    if "time" in ds.coords:
        new_coords["time"] = ds.coords["time"]

    new_ds = xr.Dataset(data_vars, coords=new_coords)

    return new_ds

def RegridDomain(base_ds : xr.Dataset, lat_points: list, lon_points: list) -> xr.Dataset:

    """
    Subset a domain dataset. Similar to Regrid(), but for domain files
    which have different conventions for spatial indexing and
    coordinate names (ie nj,ni, etc)
    """
    
    #<xarray.Dataset> Size: 24MB
    #Dimensions:  (nj: 360, ni: 720, nv: 4)
    #Dimensions without coordinates: nj, ni, nv
    #Data variables:
    #    xc       (nj, ni) float64 2MB ...
    #    yc       (nj, ni) float64 2MB ...
    #    xv       (nv, nj, ni) float64 8MB ...
    #    yv       (nv, nj, ni) float64 8MB ...
    #    mask     (nj, ni) int32 1MB ...
    #    area     (nj, ni) float64 2MB ...
    #Attributes:
    #    case_title:  GSWP3 3-Hourly Atmospheric Forcing
    #lon_2d = ds["xc"].values
    #lat_2d = ds["yc"].values
    
    # Extract the coordinate meshes for index mapping
    lat_mesh = base_ds["yc"].values  # Shape: (nj,ni)
    lon_mesh = base_ds["xc"].values  # Shape: (nj,ni)

    domain_max_lon = np.nanmax(lon_mesh)
    if domain_max_lon > 180.0:
        print(f"--> Detected Domain Convention: 0 to 360 (Max lon found: {domain_max_lon:.2f})")
        # Convert input targets to 0-360
        target_lons_adjusted = np.array(lon_points) % 360
    else:
        print(f"--> Detected Domain Convention: -180 to 180 (Max lon found: {domain_max_lon:.2f})")
        # Convert input targets to -180 to 180
        target_lons_adjusted = ((np.array(lon_points) + 180) % 360) - 180

    lat_points = np.array(lat_points)
    
    # -----------------------------------------------------------------
    # STEP A: Map target coordinates to raw 2D index pairs
    # -----------------------------------------------------------------
    matched_lat_indices = []
    matched_lon_indices = []

    for t_lat, t_lon in zip(lat_points, lon_points):

        diff_squared = (lat_mesh - t_lat) ** 2 + (lon_mesh - t_lon) ** 2
        lat_idx, lon_idx = np.unravel_index(np.argmin(diff_squared), diff_squared.shape)

        print(f"pulling lat: {lat_mesh[lat_idx,lon_idx]}, lon: {lon_mesh[lat_idx,lon_idx]}")
        matched_lat_indices.append(lat_idx)
        matched_lon_indices.append(lon_idx)
    
    # Convert lists to NumPy arrays for advanced indexing
    lat_indices = np.array(matched_lat_indices)
    lon_indices = np.array(matched_lon_indices)

    # 1. Reshape the 1D index arrays into 2D arrays with shape (1, number_of_points)
    nj_indices_2d = lat_indices.reshape(1, -1)  # Shape: (1, N)
    ni_indices_2d = lon_indices.reshape(1, -1)  # Shape: (1, N)

    # 2. Convert them to xarray DataArrays, explicitly naming the target dimensions.
    # This forces xarray to map the extracted data to 'nj' and 'ni'.
    nj_da = xr.DataArray(nj_indices_2d, dims=("nj", "ni"))
    ni_da = xr.DataArray(ni_indices_2d, dims=("nj", "ni"))

    # 3. Use advanced indexing to slice the original domain file
    # 2D variables will be (1, ni); 3D variables will automatically become (nv, 1, ni)
    new_ds = base_ds.isel(nj=nj_da, ni=ni_da)

    return new_ds

def TrimDinLoc(text_str):

    # Remove the DIN_LOC strings...
    
    if '$DIN_LOC_ROOT_CLMFORC' in text_str:
        path_end = text_str.split('$DIN_LOC_ROOT_CLMFORC')[1]
        path_base = din_loc_root_clmforc
        new_path_base = new_loc_root_clmforc
    else:
        if '$DIN_LOC_ROOT' in text_str:
            path_end = text_str.split('$DIN_LOC_ROOT')[1]
            path_base = din_loc_root
            new_path_base = new_loc_root
        else:
            print('stream dir NOT found')
            exit(2)

    return path_end,path_base,new_path_base
    


# 1) Read in user input
# 3) Create file structure
# 4) Create atm domain file
# 5) Perform regridding


din_loc_root = '/dvs_ro/cfs/cdirs/e3sm/inputdata/' 
din_loc_root_clmforc = '/dvs_ro/cfs/cdirs/e3sm/inputdata/atm/datm7/' 

new_loc_root = '/global/cfs/cdirs/m2420/rgknox/site_drivers/ZF2.2/inputdata/'
new_loc_root_clmforc = '/global/cfs/cdirs/m2420/rgknox/site_drivers/ZF2.2/inputdata/atm/datm7'

stream_list_file = '/global/homes/r/rgknox/E3SM/components/data_comps/datm/cime_config/namelist_definition_datm.xml'

# Surface Files
# It might be useful to convert several surface files,
# like if you are running both pre and post industrial
base_surf_files = ['/dvs_ro/cfs/cdirs/e3sm/inputdata/lnd/clm2/surfdata_map/surfdata_360x720cru_simyr2000_c180216.nc',
                   '/dvs_ro/cfs/cdirs/e3sm/inputdata/lnd/clm2/surfdata_map/surfdata_360x720cru_simyr1850_c180216.nc']

#target_mode = "ELMGSWP3w5e5"
target_mode = "CLMGSWP3v1"


# 2) Read in namelist definitions for the domain file
# Why do we also make a copy of the met domain file, instead of
# re-using the existing one? because I don't want to mess
# with the stream file definitions, this keeps them simpler
# by simply updating the DIN_LOC etc...


tree = ET.parse(stream_list_file)
stream_root = tree.getroot()
domain_dir = GetXmlVals(stream_root,"strm_domdir","stream",target_mode)
domain_file = GetXmlVals(stream_root,"strm_domfil","stream",target_mode)
domain_end,domain_base,new_domain_base = TrimDinLoc(domain_dir[0]+'/'+domain_file[0])
new_domain_file = new_domain_base+domain_end
new_domain_path = Path(new_domain_file)
new_domain_dir = new_domain_path.parent

print(f"Will create directory: '{new_domain_dir}'")
print(f" for modified version of domain file: {domain_end}")
response = input("Continue? (Y/N): ").strip().lower()

if response != "y":
    print("Exiting.")
    sys.exit()

new_domain_path.parent.mkdir(parents=True, exist_ok=True)

# Lets convert that domain file
base_ds = xr.open_dataset(domain_base+domain_end,engine='netcdf4')
new_ds = RegridDomain(base_ds,lat_flat,lon_flat)
new_ds.to_netcdf(new_domain_file)

print(f"Created new domain dataset:\n")
print(new_ds)
response = input("Continue? (Y/N): ").strip().lower()
if response != "y":
    print("Exiting.")
    sys.exit()

base_ds.close()
new_ds.close()

# lets convert the surface files and put them with the domain file
for base_surf_file in base_surf_files:
    base_surf_path = Path(base_surf_file)
    surf_file_name = base_surf_path.name
    # Put the new surface files next to the domain file
    new_surf_file = new_domain_dir+'/'+surf_file_name
    
    base_ds = xr.open_dataset(base_surf_file,engine='netcdf4')
    new_ds = RegridSurface(base_ds,lat_flat,lon_flat)
    new_ds.to_netcdf(new_surf_file)



# Now, for the met data,
# 1) Identify the streams needed from the configuration
# 2) go through those streams and create the folders for the new output
# 3) find the indices we will pull from by querying the first file
# 4) go through and convert all files in those directories and put result in new folders


stream_list = GetXmlVals(stream_root,"streamslist","datm_mode",target_mode)
#GetXmlVals(root,"strm_datdir","stream","list")

stream_paths = []
for stream in stream_list:

    stream_dir = GetXmlVals(stream_root,"strm_datdir","stream",stream)

    stream_path_end,stream_path_base,new_path_base = TrimDinLoc(stream_dir[0])

    stream_path_dir = stream_path_base+stream_path_end
    stream_path = Path(stream_path_dir)
    new_path_dir = new_path_base+stream_path_end
    new_path = Path(new_path_dir)

    if os.path.isdir(new_path_dir):
        print('folder: '+new_path.parent" already exists")
    else:
        new_path.mkdir(parents=True)
        print('creating folder: '+new_path_dir)

    if not os.path.isdir(stream_path_dir):
        print(f"Critical error, cannot find stream path: {stream_path_dir}")
        exit(2)

    # Loop only through NetCDF files
    
    for ifp,file_path in enumerate(stream_path.glob("clmforc*.nc")):

        base_file_path = stream_path_dir+'/'+file_path.name
        new_file_path = new_path_dir+'/'+file_path.name
        
        base_ds = xr.open_dataset(base_file_path,engine='netcdf4')
        print(f"Converting: {file_path.name} to {new_file_path}")
        if(ifp==0):
            lat_indices,lon_indices = regrid_to_points(base_ds, lat_flat, lon_flat) 

        new_ds = Regrid(base_ds,lat_indices,lon_indices)

        new_ds.to_netcdf(new_file_path)
        base_ds.close()
        new_ds.close()    

