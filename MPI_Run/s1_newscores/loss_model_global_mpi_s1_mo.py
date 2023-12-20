import numpy as np
import xarray as xr
import func_model_s1
import pickle
from scipy.interpolate import interp1d
import MPI_array_allocation
import fscore
from mpi4py import MPI
from hyperopt import STATUS_OK


def loss(params):
    # Starting MPI comm world
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    eo_score = None
    mo_score = None

    # Loading files on core 0
    if rank == 0:
        init = xr.open_dataset('rsFile_modified.nc') # Restart file
        ignd = init.dims['layer']
        iccp2 = init.dims['iccp2']
        ncell = 1507

        # Now create output split and disp
        arg1_output = np.zeros([ncell,ignd,20440])
        arg1_output = np.ascontiguousarray(arg1_output)
        split_arg1_output, split_sizes_arg1_output, disp_arg1_output = MPI_array_allocation.fallocate(arg1_output, size)
        arg2_output = np.zeros([ncell,ignd,20440])
        split_arg2_output, split_sizes_arg2_output, disp_arg2_output = MPI_array_allocation.fallocate(arg2_output,size)
        arg3_output = np.zeros([ncell,20440])
        split_arg3_output, split_sizes_arg3_output, disp_arg3_output = MPI_array_allocation.fallocate(arg3_output, size)

    else:
        ignd = None
        iccp2 = None

        arg1_output = None
        split_arg1_output = None
        split_sizes_arg1_output = None
        disp_arg1_output = None
        arg2_output = None
        split_arg2_output = None
        split_sizes_arg2_output = None
        disp_arg2_output = None
        arg3_output = None
        split_arg3_output = None
        split_sizes_arg3_output = None
        disp_arg3_output = None

    ignd = comm.bcast(ignd, root=0)
    iccp2 = comm.bcast(iccp2, root=0)

    # Same for outputs
    split_sizes_arg1_output = comm.bcast(split_sizes_arg1_output, root=0)
    disp_arg1_output = comm.bcast(disp_arg1_output, root=0)
    split_sizes_arg2_output = comm.bcast(split_sizes_arg2_output, root=0)
    disp_arg2_output = comm.bcast(disp_arg2_output, root=0)
    split_sizes_arg3_output = comm.bcast(split_sizes_arg3_output, root=0)
    disp_arg3_output = comm.bcast(disp_arg3_output, root=0)

    # Loading input files
    litrmass = np.load('litrmass_eq'+str(rank)+'.npy') # Initial litrmass
    soilcmass = np.load('soilcmass_eq'+str(rank)+'.npy') # Initial soilcmass
    actlyr = np.load('actlyr'+str(rank)+'.npy') # Active layer
    fcancmx = np.load('fcancmx'+str(rank)+'.npy') # pft distribution of gridcells
    SAND = np.load('SAND'+str(rank)+'.npy') # SAND content
    CLAY = np.load('CLAY'+str(rank)+'.npy') # CLAY content
    zbotw = np.load('zbotw'+str(rank)+'.npy') # zbotw
    delzw = np.load('delzw'+str(rank)+'.npy') # thickness of permeable soil layer
    tbar = np.load('tbar'+str(rank)+'.npy') # Soil temperature
    mliq = np.load('mliq'+str(rank)+'.npy') # Liquid water content
    mice = np.load('mice'+str(rank)+'.npy') # Initial litrmass
    Cinput = np.load('Cinput'+str(rank)+'.npy') # Initial litrmass
    
# Computing the rest of variables needed for model run
    isand = np.where((SAND == -4) | (SAND == -3))  # Flag ,-4 = ice, -3 = ice
    ncell = CLAY.shape[1]

    thliq = mliq / 1000 / delzw  # Conversion to m^3/m^3
    thice = mice / 1000 / delzw  # Conversion to m^3/m^3

    thpor = (-0.126 * SAND + 48.9) / 100  # Soil total porority [m^3 m^-3]
    thpor[isand] = 0  # No value where soil is rock or ice
    b = 0.159 * CLAY + 2.91  # Clapp and Hornberger empirical param []
    b[isand] = 0  # No value where soil is rock or ice
    psisat = 0.01 * np.exp(-0.0302 * SAND + 4.33)  # Soil moisture suction at saturation [m]
    psisat[isand] = 0

    # Global variables
    ndays = tbar.shape[0]
    tf = int((ndays)/365)

    itermax = 500000
    delt = 1/365
    spinup = 50000

    # START MODEL RUN
    arg1, arg2, arg3 = func_model_s1.model(params,ncell,iccp2,ignd,litrmass,soilcmass,tf,delt,ndays, thice, thpor, thliq,
                                    psisat, tbar, b, itermax,spinup,zbotw,actlyr,SAND,fcancmx,delzw,
                                    Cinput,rank)

    # Sending model output back to main core
    # Transpose output arrays to gather them
    arg1 = np.ascontiguousarray(arg1.transpose(1,0,2))
    arg2 = np.ascontiguousarray(arg2.transpose(1,0,2))
    arg3 = np.ascontiguousarray(arg3)

    # Gather output data together on core 0
    comm.Gatherv(arg1, [arg1_output, split_sizes_arg1_output, disp_arg1_output, MPI.DOUBLE],root=0)
    comm.Gatherv(arg2, [arg2_output, split_sizes_arg2_output, disp_arg2_output, MPI.DOUBLE], root=0)
    comm.Gatherv(arg3, [arg3_output, split_sizes_arg3_output, disp_arg3_output, MPI.DOUBLE], root=0)

    # Once outputs are gathered, we compute the score
    if rank == 0:
        # Transpose output arrays back to original shape
        arg1_output = arg1_output.transpose(1,0,2)
        arg2_output = arg2_output.transpose(1,0,2)
        arg_output = arg2_output + arg1_output

        init = xr.open_dataset('rsFile_modified.nc') # init file
        filtrd = np.load('filtrd_index_5p.npy')[0] # init file
        rmrveg = np.load('rmrveg_complete.npy') # loading autotrophic respiration
        rmrveg = rmrveg[261:317]
        #rmrveg = rmrveg[:,filtrd]
        # Reading in formated woSIS/srdb dataset
        with open('wosis_srdb_grid_5p.txt', 'rb') as fp:
            grid = pickle.load(fp)
        grid = np.array(grid, dtype=object)
        grid = grid[filtrd]

        # Reading in dataset_flag array
        dataset_flag = np.load('dataset_flag_5p.npy')
        dataset_flag = dataset_flag[filtrd]

        # Depth of soil layers
        delz = init['DELZ'].data[:]                         # Ground layer thickness [m]
        zbot = np.zeros_like(delz)                          # Dept of the bottom of each soil layer [m]
        for i in range(0, len(delz)):
            zbot[i] = np.sum(delz[:i+1])
        classic_zbot = np.zeros(len(zbot)+1)
        classic_zbot[1:] = zbot

        # Wosis score
        iwosis = np.where((dataset_flag == 1) | (dataset_flag== 2))
        wosis_score = fscore.custom_score_wosis(grid[iwosis], arg_output[:,iwosis,:][:,0,:,:],classic_zbot, 0.7995)

        # Srdb score
        isrdb = np.where((dataset_flag == 2) | (dataset_flag == 3))
        srdb_score = fscore.custom_score_srdb(grid[isrdb], arg3_output[isrdb],rmrveg[:,isrdb[0]])
        # Computing score
        alpha = 0.5
        #score = alpha*wosis_score[1] + (1-alpha)*srdb_score[1]
        eo_score = alpha*wosis_score[0] + (1-alpha)*srdb_score[0]
        mo_score = alpha*wosis_score[1] + (1-alpha)*srdb_score[1]

    eo_score = comm.bcast(eo_score,root=0)
    mo_score = comm.bcast(mo_score,root=0)
    return {
        'loss': mo_score,
        'status': STATUS_OK,
        'other_score': eo_score
    }

