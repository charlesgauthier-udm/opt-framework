import numpy as np
# This scripts compute the Carbon input to the soil carbon module from rmatctem,fVegLitter,fStemLitter, fLeafLitter and
# npp

##### SECTION BELOW DEALS WITH SUBSETS OF RMATCTEM, UNCOMMENT IF SUBSETS OF RMATCTEM WERE USED ########
# # We concatenate the rmatctem files into one array
# final = np.empty([44165,11,20,0],dtype='float32') # Array to receive the subsets
# for i in range(0,34): # CHANGE TO NUMBER OF SUBSETS OF RMATCTEM
#    rmatctem = np.load('/home/chargaut/scratch/formating_global_files/formated/rmatctem'+str(i)+'.npy')
#    final = np.concatenate([final,rmatctem],axis=3)
# np.save('rmatctem_complete.npy',final) # Save the final complete array
########################################################################################################

# Once we have rmatctem we can compute that Carbon input coming to the soil carbon module
rmatctem = np.load('/home/chagau1/opt-framework/rmatctem_complete.npy')
fLeafLitter = np.load('/home/chagau1/opt-framework/fLeafLitter_complete.npy')
fRootLitter = np.load('/home/chagau1/opt-framework/fRootLitter_complete.npy')
fStemLitter = np.load('/home/chagau1/opt-framework/fStemLitter_complete.npy')
npp = np.load('/home/chagau1/opt-framework/npp_complete.npy')


root_Cinput = (fRootLitter * rmatctem.transpose(2,0,1,3)).transpose(1,2,0,3) # Carbon input through roots
del rmatctem
del fRootLitter

Cinput = np.zeros_like(root_Cinput) # Empty array to receive
Cinput[:,:,0,:] = fLeafLitter+fStemLitter+npp # Inputs to the first soil layer
del fLeafLitter
del fStemLitter

Cinput += root_Cinput
Cinput = Cinput * 3.1536e7 # Conversion from [kg C/m2/s] --> [kg C/m2/yr]

# Saving complete carbon input file
np.save('Cinput_complete.npy',Cinput)
