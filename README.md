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
