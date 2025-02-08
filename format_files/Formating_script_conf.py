import numpy as np
import xarray as xr
import time
import pickle
# Loading optimization lats and lons
lats = np.load('opt_lats.npy')
lons = np.load('opt_lons.npy')

# Loading input files to be formatted
cSoilperlay_annually = xr.open_dataset('/home/chagau1/JAMES/Data/cSoilperlay_annually.nc').isel(time=slice(261,318))
rSoil_monthly = xr.open_dataset('/home/chagau1/JAMES/Data/rSoil_monthly.nc').isel(time=slice(12,696))

# Global variable
ncell = len(lats)
nyear = 57
nmonths = 684
ignd = 20  # init.dims['layer']


# Creating empty arrays to store formatted variables
cSoilperlay_a = np.zeros([nyear,ignd,ncell])
rSoil_monthly_a = np.zeros([nmonths,ncell])

# Looping on all gridcell and adding corresponding variable to their empty arrays
for i in range(0,ncell):
    t1 = time.time()
    cSoilperlay_a[:,:,i] = cSoilperlay_annually.cSoilperlay.sel(latitude=lats[i], longitude=lons[i])
    rSoil_monthly_a[:,i] = rSoil_monthly.rSoil.sel(latitude=lats[i], longitude=lons[i])

    tf = time.time() - t1
    print('Time spent on current iteration: ', int(tf//60), 'minute ',round(tf%60,0), 'secondes')
    print(i+1, '/', ncell)

np.save('cSoilperlay_complete',cSoilperlay_a)
np.save('rSoil_monthly_complete', rSoil_monthly_a)
