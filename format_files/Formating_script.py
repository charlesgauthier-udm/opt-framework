import numpy as np
import xarray as xr
import time
import frmatctem

# Loading lats & lons of gridcell that contain data
lats = np.load('/home/charlesgauthier/project/global_opt_files/grid_opt_lats_5p.npy')
lons = np.load('/home/charlesgauthier/project/global_opt_files/grid_opt_lons_5p.npy')

# START OF BLOCK OF CODE TO COMPUTE FILTRD INDEX
# We then use actlyr because it is a small file and it allows us to see which gridcells containing observations are not simulated by the model
actlyr = xr.open_dataset('/home/charlesgauthier/project/global_opt_files/actlyrmax_annually.nc')
init = xr.open_dataset('/home/charlesgauthier/project/global_opt_files/rsFile_modified.nc')
ncell = len(lats)
nyears = len(actlyr.time.data)

actlyr_tab = np.zeros([nyears,ncell])
for i in range(0,ncell):
    #print(i)
    if lats[i] < -57.20663:
        actlyr_tab[:, i] = np.nan
        print(lats[i])
    else:
        actlyr_tab[:, i] = actlyr['actlyrmax'].sel(latitude=lats[i], longitude=lons[i])


# When then identify index of gridcells that contain model values and we save those indices
np.save('/home/charlesgauthier/project/global_opt_files/filtrd_index_5p',np.where(~np.isnan(actlyr_tab[0])))
# END OF BLOCK OF CODE TO FORMAT FILTRDINDEX


# BLOCK OF CODE TO REMOVE CLASSIC GRIDCELL WITH OBS. IN THEM THAT ARE NOT SIMULATED BY THE MODEL
#Removing gridcells where the model doesn't simulate anything
ind = np.load('/home/charlesgauthier/project/global_opt_files/filtrd_index_5p.npy')
lats = lats[ind[0]]
lons = lons[ind[0]]
np.save('/home/charlesgauthier/project/global_opt_files/opt_lats.npy',lats)
np.save('/home/charlesgauthier/project/global_opt_files/opt_lons.npy',lons)
# END OF BLOCK OF CODE TO REMOVE GRIDCELLS WITH OBS NOT SIMULATED BY THE MODEL

# Formating starts here. Here is some precisions on the few lat/lons files
# grid_lats_5p = lats/lon of all gridcells containing observations regardless of if the cell is simulated by the model
# opt_lats = lat/lon of gridcells containing observations and that are simulated by the model <- this one is used to format files

# Loading optimization lats and lons

lats = np.load('/home/charlesgauthier/project/global_opt_files/opt_lats.npy')
lons = np.load('/home/charlesgauthier/project/global_opt_files/opt_lons.npy')

# Loading input files to be formatted
init = xr.open_dataset('/home/charlesgauthier/project/global_opt_files/rsFile_modified.nc')
rmrveg = xr.open_dataset('/home/charlesgauthier/project/global_opt_files/rmrveg_annually.nc')
sftlf = xr.open_dataset('/home/charlesgauthier/project/global_opt_files/sftlf.nc')
actlyrmax = xr.open_dataset('/home/charlesgauthier/project/global_opt_files/actlyrmax_annually.nc')
npp = xr.open_dataset('/home/charlesgauthier/project/global_opt_files/npp_daily_perpft.nc')
fRootLitter = xr.open_dataset('/home/charlesgauthier/project/global_opt_files/fRootLitter_daily_perpft.nc')
fStemLitter = xr.open_dataset('/home/charlesgauthier/project/global_opt_files/fStemLitter_daily_perpft.nc')
fLeafLitter = xr.open_dataset('/home/charlesgauthier/project/global_opt_files/fLeafLitter_daily_perpft.nc')
cRoot = xr.open_dataset('/home/charlesgauthier/project/global_opt_files/cRoot_daily_perpft.nc')
mrsfl = xr.open_dataset('/home/charlesgauthier/project/global_opt_files/mrsfl_daily.nc')
mrsll = xr.open_dataset('/home/charlesgauthier/project/global_opt_files/mrsll_daily.nc')
tsl = xr.open_dataset('/home/charlesgauthier/project/global_opt_files/tsl_daily.nc')

# Global variable
ncell = len(lats)
ndays = 44165  # tsl.dims['time']
nyear = int(ndays / 365)
ignd = 20  # init.dims['layer']
iccp2 = 11  # init.dims['iccp2']

# Creating empty arrays to store formatted variables
grclarea = np.zeros(ncell)
fcancmx = np.zeros([iccp2,ncell])
litrmass = np.zeros([ignd,iccp2,ncell])
soilcmass = np.zeros([ignd,iccp2,ncell])
CLAY = np.zeros([ignd,ncell])
SAND = np.zeros([ignd,ncell])
sdepth = np.zeros(ncell)
sftlf_a = np.zeros([ncell])
rmrveg_a = np.zeros([321, ncell])
actlyrmax_a = np.zeros([321,ncell])
fLeafLitter_a = np.zeros([ndays,iccp2,ncell])
fStemLitter_a = np.zeros([ndays,iccp2,ncell])
fRootLitter_a = np.zeros([ndays,iccp2,ncell])
npp_a = np.zeros([ndays,iccp2,ncell])
cRoot_a = np.zeros([ndays,iccp2,ncell])
mice = np.zeros([ndays, ignd, ncell])
mliq = np.zeros([ndays, ignd, ncell])
tbar = np.zeros([ndays,ignd,ncell])

# Looping on all gridcell and adding corresponding variable to their empty arrays
for i in range(0,ncell):
    t1 = time.time()
    rmrveg_a[:,i] = rmrveg.rmrveg.sel(latitude=lats[0], longitude=lons[0])
    grclarea[i] = init['grclarea'].sel(lat=lats[i],lon=lons[i])
    fcancmx[:9,i] = init['fcancmx'].sel(lat=lats[i], lon=lons[i])
    litrmass[:,:,i] = init['litrmass'].sel(lat=lats[i],lon=lons[i])[0]
    soilcmass[:,:,i] - init['soilcmas'].sel(lat=lats[i],lon=lons[i])[0]
    CLAY[:,i] = init['CLAY'].sel(lat=lats[i], lon=lons[i])[0]
    SAND[:,i] = init['SAND'].sel(lat=lats[i], lon=lons[i])[0]
    sdepth[i] = init['SDEP'].sel(lat=lats[i], lon=lons[i])[0]
    sftlf_a[i] = sftlf['sftlf'].sel(latitude=lats[i],longitude=lons[i])
    actlyrmax_a[:,i] = actlyrmax['actlyrmax'].sel(latitude=lats[i],longitude=lons[i])
    fLeafLitter_a[:,:9,i] = fLeafLitter['fLeafLitter'].sel(latitude=lats[i],longitude=lons[i])
    fStemLitter_a[:,:9,i] = fStemLitter['fStemLitter'].sel(latitude=lats[i],longitude=lons[i])
    fRootLitter_a[:,:9,i] = fRootLitter['fRootLitter'].sel(latitude=lats[i],longitude=lons[i])
    npp_a[:,:9,i] = npp.npp.sel(latitude=lats[i],longitude=lons[i])
    cRoot_a[:,:9,i] = cRoot['cRoot'].sel(latitude=lats[i],longitude=lons[i])
    mice[:,:,i] = mrsfl['mrsfl'].sel(latitude=lats[i],longitude=lons[i])
    mliq[:,:,i] = mrsll['mrsll'].sel(latitude=lats[i],longitude=lons[i])
    tbar[:,:,i] = tsl.tsl.sel(latitude=lats[i],longitude=lons[i])

    tf = time.time() - t1
    print('Time spent on current iteration: ', int(tf//60), 'minute ',round(tf%60,0), 'secondes')
    print(i+1, '/', ncell)

np.save('rmrveg_complete',rmrveg_a)
np.save('grclarea_complete',grclarea)
np.save('fcancmx_complete',fcancmx)
np.save('litrmass_complete',litrmass)
np.save('CLAY_complete',CLAY)
np.save('/SAND_complete',SAND)
np.save('sdepth_complete',sdepth)
np.save('sftlf_complete',sftlf_a)
np.save('actlyrmax_monthly_complete',actlyrmax_a)
np.save('fStemLitter_complete',fStemLitter_a)
np.save('fRootLitter_complete',fRootLitter_a)
np.save('npp_complete.npy',npp_a)
np.save('cRoot_complete',cRoot_a)
np.save('mice_complete',mice)
np.save('mliq_complete',mliq)
np.save('tbar_complete',tbar)

# Computing maximum active layer depth daily from actlyrmax_monthly
actlyrmax_a = np.load('actlyrmax_monthly_complete.npy')[200:] # ensures we start at year 1900
maxactlyr_daily = np.zeros([ndays, ncell])
count = 0
while count < ndays:
    f =count // 365
    maxactlyr_daily[count,:] = actlyrmax_a[f,:]
    count +=1
    print(f)
np.save('maxactlyr_daily_complete',maxactlyr_daily)


# Computing intermediate variables zbotw and delzw
SAND = np.load('SAND_complete.npy')
sdepth = np.load("sdepth_complete.npy")
isand = np.where((SAND == -4) | (SAND == -3))       # Flag for soil type -4 = ice, -3 = rock
delz = init['DELZ'].data[:]                         # Ground layer thickness [m]#
zbot = np.zeros_like(delz)                          # Dept of the bottom of each soil layer [m]
for i in range(0, len(delz)):
    zbot[i] = np.sum(delz[:i+1])
delz = np.tile(delz,(ncell,1)).transpose()
zbot = np.tile(zbot,(ncell,1)).transpose()
delzw = np.zeros_like(delz)                         # thickness of permeable part of soil layer
isice = np.where(SAND == -4)                        # Soil layers that are ice#
isrock = np.where(SAND == -3)                       # Soil layers that are rock
isbelow = np.where(sdepth >= zbot)                  #
isnear = np.where(sdepth < (zbot - delz + 0.025))

delzw = np.maximum(0.05, (sdepth - (zbot - delz)))

delzw[isnear] = 0.0
SAND[isnear] = -3

delzw[isbelow] = delz[isbelow]
delzw[isrock] = 0
for i in range(0,ncell):
    if SAND[0,i] == -4:
        delzw[:,i] = delz[:,i]
        SAND[:,i] = -4

zbotw = np.maximum(0, (zbot - delz)) + delzw          # Bottom of soil layers [m]

np.save('delzw_complete',delzw)
np.save('zbotw_complete', zbotw)

# Computing and formating rmatctem, distribution of carbon from roots
cRoot = np.load('cRoot_complete.npy')
maxannualactlyr = np.load('maxactlyr_daily_complete.npy')
zbotw = np.load('zbotw_complete.npy')
delzw = np.load('delzw_complete.npy')
sdepth = np.load('sdepth_complete.npy')

# PFT-specific parameters provided in the config file
alpha = np.array([0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0,0])
avertmas = np.array([1.85, 1.45, 2.45, 2.10, 2.10,0.1, 0.1, 0.7,0.7,0,0])
abar = np.array([4.70,5.86,3.87,3.46,3.97,3.97,3.97,5.86,4.92,4,4])
mxrtdpth = np.array([3,3,5,5,3,2,2,1,1,0,0])

# Here we split the arrays to compute subset of rmatctem because it is a too large array to handle at once
ncell = 10
tot_ncell = 1507
split = tot_ncell // ncell +1

cRoot_split = np.array_split(cRoot, split, axis=2)
maxannualactlyr_split = np.array_split(maxannualactlyr, split, axis=1)
zbotw_split = np.array_split(zbotw,split, axis=1)
sdepth_split = np.array_split(sdepth, split)

final = np.empty([44165,11,20,0],dtype='float32')
for i in range(0,split):
    t1 = time.time()
    rmatctem = frmatctem.rmatctemf2(zbotw_split[i],alpha,cRoot_split[i],avertmas,abar,sdepth_split[i],maxannualactlyr,mxrtdpth)
    rmatctem.astype('float32')

    np.save('/home/charlesgauthier/project/global_opt_files/opt_files/rmatctem_complete'+str(i)+'.npy',rmatctem)
    tf = time.time() - t1
    print('Time spent on current iteration: ', int(tf//60), 'minute ',round(tf%60,0), 'secondes')
    print(i,' / ',split)
