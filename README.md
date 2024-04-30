# opt-framework
Optimisation framework to optimize the soil carbon parameters of CLASSIC.

`MPI_run` folder contains the script needed to perform the optimization runs.

`data_set` contains the scripts related to the dataset that the framework uses to optimize the model against.

`format_files` folder contains the script needed to format the typicall netcdf files produced by CLASSIC to the numpy arrays needed by the framework.

`input_files` is an empty folder that has to be filled with the `.npy` files

`lib_fixed_turbation` folder contains the python wrapping of the modified fortran turbation routine used by the model.

## Getting started

### Input files
First ensure that all `.npy` files needed by the framework are located in the `opt-framework/input_files/` directory as the launch scripts will look for the files there.


### Virtual Environment
Create a virtual environment and ensure that all the packages listed in the `requirements.txt` file are correctly installed. The `lib_fixed_turbation` package has to be manually installed as it is a python wrapper of a modified version of the `soilCProcesses.f90` fortran subroutine. To install it, ensure that your virtual environment is activated, ensure that you are in the `/opt-framework/lib_fixed_turbation/` directory where the `setup.py` file is located and install the package using pip:


    python -m pip install .

### Running a job
Once the input file are available and once the environment has been set up, the framework should be ready to go. To launch an optimization job, go to the `/opt-framework/MPI_Run/runs/` directory and launch the `.sh` script corresponding to the optimization job you want to run. If you are starting the optimization from the start, you need to change the value of the `--iter` argument in the `mpiexex` call to the number of optimization iteration you want to launch. Once the optimization reaches the specified number of optimization steps, a `.pkl` file is generated and moved to the `/opt-framework/MPI_Run/outputs` directory. This file is a dictionnary containing information on the optimization run and can also be used to restart the optimization from that point. To restart the optimization from a `.pkl` file, a few things have to be changed in the `.sh` job file. First, the `cp` command line moving the `.pkl` file to the tmp-dir has to be uncommented so that the file is available to the framework. Then, the `--outfile` argument has to be added to the `mpiexec` call with only the name of the `.pkl` file e.g. `--outfile=s1_eo_opt_out_1200.pkl`

## Generating new input files
To generate new input files usable by the framework, the first step is to generate the nested list. In order to save on computation, the framework is set up to run on a vector of gridcells containing WoSIS/SRDB data. The nested list calculates where the observations fall within the gridcells.
To generate the list, execute the `data_sets/nested_list.py` script. The script needs three files:
1. A `rsFile_modified.nc` which gives information on the CLASSIC grid used
2. The `woSIS_soilc_sgbd_5p.nc` NetCDF file containing the WoSIS soil profiles
3. The `srdb_formatted.nc` NetCDF file containing the SRDB observations

Once completed, the scripts saves a `wosis_srdb_grid_5p.txt` file which is the vector of CLASSIC grid cells containing observation. The scripts also generates a pair of `grid_opt_lats`/`grid_opt_lons` which are the coordinates of every grid cell in the vector. Information on how the "grid" is constructed is available in the `nested_list.py` script.

Once the wosis_srdb_grid file is generated, it is possible to format new input files. To do so, execute the `/format_files/Formating_script.py` script.
This script takes in all the netCDF files of the CLASSIC variables needed in the soil carbon module and converts them to numpy arrays with a flatten gridcell axis similar to the `wosis_srdb_grid` vector.

Due to the original set up in which the framework was created, some extra steps had to be taken to deal with computational limitations. This is why the rmatctem variable can be split and computed in many smaller arrays (NOT TO BE CONFUSED WITH PRE_SPLITING THE FILES FOR MPI).

With either the split or whole rmatctem files in hand, run the `/format_files/rmatctem_2_Cinput.py` script which computes the C input to the soil carbon module using rmatctem,fLeaf,fStem and fRootLitter and npp.

Once the Cinput file has been computed, everything is in place to pre-splt the files for the optimization in parallel. To do so, execute the `format_files/pre_splitting_files.py` script.