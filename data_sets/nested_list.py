# Scrpt that creates list of gridcells containing woSIS and SRDB V5 data for optimization
# Importing modules
import numpy as np
import xarray as xr
import time
from scipy.interpolate import interp1d
import pickle
import matplotlib.pyplot as plt

# Importing CLASSIC outputs for gridcell dimensions and properties
CLASSIC = xr.open_dataset('/home/charlesgauthier/project/global_opt_files/rsFile_modified.nc')  # Loading netcdf file
actlyr = xr.open_dataset('/media/charlesgauthier/Seagate Backup Plus Drive/global_opt_files_latest/actlyrmax_annually.nc')
classic_lon = np.array(CLASSIC['lon'].data)     # Assigning longitudes [0, 360 degrees]
classic_lat = np.array(CLASSIC['lat'].data)     # Assigning latitudes [-90, 90 degrees]
classic_ipeat = CLASSIC['ipeatland'].data[0,:,:]

# Importing wosis observation data
woSIS = xr.open_dataset('/media/charlesgauthier/Seagate Backup Plus Drive/datasets/WoSIS/woSIS_soilc_sgbd_5p.nc') # Loading netcdf
woSIS_lon = woSIS['longitude'].data                 # Assigning longitude
woSIS_lat = woSIS['latitude'].data                  # Assigning latitude

# Importing srdb observation data
srdb = xr.open_dataset('/media/charlesgauthier/Seagate Backup Plus Drive/datasets/SRDB_V5/srdb_formatted.nc')
srdb_lon = srdb['longitude'].data       # Assigning longitude
srdb_lat = srdb['latitude'].data        # Assigning latitude
srdb_resp = srdb['rs_annual'].data      # Annual soil C resp. [kg C m^-2]
srdb_resp_std = srdb['rs_annual'].data  # Annual soil C resp std [kg C m^-2]
srdb_start = srdb['start_year'].data    # Start year
srdb_stop = srdb['stop_year'].data      # Stop year

# Comparison Classic/woSIS/srdb
classic_lonw = np.zeros(len(classic_lon)+1)  # Creating array like classic_lon with one more value for a wraparound longitude
classic_lonw[-1] = classic_lon[-1] + np.abs(classic_lon[0] - classic_lon[1])  # Adding value for wrap around longitude
classic_lonw[:-1] = classic_lon

'''
Building nested list in this fashion: Every element of the list corresponds to a CLASSIC gridcell. Instead of keeping 
gridcells in a 2D array format, they are flatten and only gridcells that have woSIS/srdb data within them are kept for 
optimization. Corresponding flatten lons & lats arrays are saved in .npy files to be read-in alongside the nested list
when optimizing. In the nested list, every element (or gridcell) is formatted like so:
0- woSIS organic carbon content [kg C m^-2]
1- woSIS lower depth of measurment [m]
2- srbd annual soil C respiration [kg C m^-2 yr^-1]
3- srdb annual soil C resp. std [kg C m^-2 yr^-1]
4- srdb start year
5- srdb stop year
6- Dataset flag: 1=Only woSIS data, 2=woSIS & srdb, 3= Only srdb

NaNs indicate that there is no data for the gridcell
'''
grid = []   # Empty nested list
grid_lats = np.array([])    # Empty array for latitudes
grid_lons = np.array([])    # Empty array for longitudes
data_wosis_lats = np.array([])
data_wosis_lons = np.array([])
data_srdb_lats = np.array([])
data_srdb_lons = np.array([])
# Looping on each grid cell and assigning each wosis/srdb datapoints that fits in the gridcell
for i in range(0, len(classic_lonw) - 1):   # Looping on longitude
    for j in range(0, len(classic_lat)):    # Looping on latitude
        print(i)
        if classic_ipeat[j,i] == 0: # Non peatlands only
            # Get index of woSIS datapoints that fit in current gridcell
            iwosis = np.argwhere((classic_lonw[i] < woSIS_lon) & (woSIS_lon <= classic_lonw[i + 1]) &
                                 (classic_lat[j - 1] < woSIS_lat) & (woSIS_lat <= classic_lat[j]))
            # Get index of srdb datapoints that fit in the gridcell
            isrdb = np.argwhere((classic_lonw[i] < srdb_lon) & (srdb_lon <= classic_lonw[i + 1]) &
                                 (classic_lat[j - 1] < srdb_lat) & (srdb_lat <= classic_lat[j]))

            if len(iwosis) != 0 and len(isrdb) == 0:    # If there is wosis data and no srdb
                # Creating list element and appending to nested list
                nested_list = [[woSIS['orgc_value'].data[iwosis]],[woSIS['lower_depth'].data[iwosis]],[np.nan],[np.nan],[np.nan],[np.nan],[1]]
                grid.append(nested_list)
                grid_lats = np.append(grid_lats, classic_lat[j])
                grid_lons = np.append(grid_lons, classic_lonw[i])
                data_wosis_lats = np.append(data_wosis_lats,woSIS_lat[iwosis][:,0])
                data__wosislons = np.append(data_wosis_lons,woSIS_lon[iwosis][:,0])
            elif len(iwosis) != 0 and len(isrdb) != 0:  # It there is wosis AND srdb data
                # Creating list element and appending to nested list
                nested_list = [[woSIS['orgc_value'].data[iwosis]],[woSIS['lower_depth'].data[iwosis]],[srdb_resp[isrdb]],[srdb_resp_std[isrdb]],[srdb_start[isrdb]],[srdb_stop[isrdb]],[2]]
                grid.append(nested_list)
                grid_lats = np.append(grid_lats, classic_lat[j])
                grid_lons = np.append(grid_lons, classic_lonw[i])
                data_wosis_lats = np.append(data_wosis_lats,woSIS_lat[iwosis][:,0])
                data_wosis_lons = np.append(data_wosis_lons,woSIS_lon[iwosis][:,0])
                data_srdb_lats = np.append(data_srdb_lats,srdb_lat[isrdb][:,0])
                data_srdb_lons = np.append(data_srdb_lons,srdb_lon[isrdb][:,0])
            elif len(iwosis) == 0 and len(isrdb) != 0:  # If there is srdb data and no wosis
                # Creating list element and appending to nested list
                nested_list = [[np.nan],[np.nan],[srdb_resp[isrdb]],[srdb_resp_std[isrdb]],[srdb_start[isrdb]],[srdb_stop[isrdb]],[3]]
                grid.append(nested_list)
                grid_lats = np.append(grid_lats, classic_lat[j])
                grid_lons = np.append(grid_lons, classic_lonw[i])
                data_srdb_lats = np.append(data_srdb_lats,srdb_lat[isrdb][:,0])
                data_srdb_lons = np.append(data_srdb_lons,srdb_lon[isrdb][:,0])

# Creating array of dataset flags
dataset_flag = np.zeros(len(grid))

# Looping on the grid and adding gridcell data flag to the dataset_flag array
for f in range(0,len(grid)):
    flag = grid[f][-1][0]
    dataset_flag[f] = flag

# Saving dataset_flag array
np.save('/home/charlesgauthier/project/global_opt_files/dataset_flag_5p', dataset_flag)
# Saving nested list
with open('wosis_srdb_grid_5p.txt','wb') as fp:
    pickle.dump(grid, fp)
#with open('dataset_spatial_grid.txt','wb') as fp:
#    pickle.dump(grid2, fp)
# Saving lons & lats associated with nested list
np.save('/home/charlesgauthier/project/global_opt_files/grid_opt_lats_5p', grid_lats)
np.save('/home/charlesgauthier/project/global_opt_files/grid_opt_lons_5p', grid_lons)




