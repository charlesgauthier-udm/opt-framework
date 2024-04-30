import numpy as np

size = 48 # CHANGE TO NUMBER OF CORES USED

# Each _complete.npy file is split along the gridcell axis which should be the last axis

### litrmass ###
litrmass = np.load('litrmass_complete.npy')
split = np.array_split(litrmass,size,axis=2)
for i in range(0,len(split)):
    np.save('litrmass'+str(i)+'.npy',split[i])
del litrmass
del split

### soilcmass ###
soilcmass = np.load('soilcmass_complete.npy')
split = np.array_split(soilcmass,size,axis=2)
for i in range(0,len(split)):
    np.save('soilcmass'+str(i)+'.npy',split[i])
del soilcmass
del split

### actlyr ###
actlyr = np.load('maxactlyr_daily.npy')
split = np.array_split(actlyr,size,axis=1)
for i in range(0,len(split)):
    np.save('actlyr'+str(i)+'.npy', split[i])
del actlyr
del split

### fcancmx ###
fcancmx = np.load('fcancmx_complete.npy')
split = np.array_split(fcancmx,size,axis=1)
for i in range(0,len(split)):
    np.save('fcancmx'+str(i)+'.npy', split[i])
del fcancmx
del split

### SAND ###
SAND = np.load('SAND_complete.npy')
split = np.array_split(SAND,size,axis=1)
for i in range(0,len(split)):
    np.save('SAND'+str(i)+'.npy', split[i])
del SAND
del split

### CLAY ###
CLAY = np.load('CLAY_complete.npy')
split = np.array_split(CLAY,size,axis=1)
for i in range(0,len(split)):
    np.save('CLAY'+str(i)+'.npy', split[i])
del CLAY
del split

### zbotw ###
zbotw = np.load('zbotw_complete.npy')
split = np.array_split(zbotw,size,axis=1)
for i in range(0,len(split)):
    np.save('zbotw'+str(i)+'.npy', split[i])
del zbotw
del split

### delzw ###
delzw = np.load('delzw_complete.npy')
split = np.array_split(delzw,size,axis=1)
for i in range(0,len(split)):
    np.save('delzw'+str(i)+'.npy', split[i])
del delzw
del split

### tbar ###
tbar = np.load('tbar_complete.npy')
split = np.array_split(tbar,size,axis=2)
for i in range(0,len(split)):
    np.save('tbar'+str(i)+'.npy', split[i])
del tbar
del split

### mliq ###
mliq = np.load('mliq_complete.npy')
split = np.array_split(mliq,size,axis=2)
for i in range(0,len(split)):
    np.save('mliq'+str(i)+'.npy', split[i])
del mliq
del split

### mice ###
mice = np.load('mice_complete.npy')
split = np.array_split(mice,size,axis=2)
for i in range(0,len(split)):
    np.save('mice'+str(i)+'.npy', split[i])
del mice
del split

### Cinput ###
Cinput = np.load('Cinput.npy')
split = np.array_split(Cinput,size,axis=3)
for i in range(0,len(split)):
    np.save('Cinput'+str(i)+'.npy', split[i])
del Cinput
del split
