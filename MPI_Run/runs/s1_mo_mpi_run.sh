#!/bin/bash
#SBATCH --time=20:00:00
#SBATCH --account=def-sonol
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem=249G
#SBATCH --mail-user=charles.gauthier.1@umontreal.ca
#SBATCH --mail-type=ALL


# Get parent directories name
runs=$(pwd)

# Get the parent directory path
MPI_Run=$(dirname "$runs")
optfdir=$(dirname "$MPI_Run")


module load python/3.8
module load scipy-stack
module load mpi4py

cd $SLURM_TMPDIR/

source ~/ENV_2/bin/activate # activate environment

# copy files to tmpdir
cp -v $MPI_Run/s1_newscores/main_s1_mo.py ./
cp -v $MPI_Run/fenv.py ./
cp -v $MPI_Run/fscore.py ./
cp -v $MPI_Run/fturbation.py ./
cp -v $MPI_Run/s1_newscores/algo_s1_mo.py ./
cp -v $MPI_Run/s1_newscores/loss_model_global_mpi_s1_mo.py ./
cp -v $MPI_Run/s1_newscores/func_model_s1.py ./
cp -v $MPI_Run/MPI_array_allocation.py ./
cp -v $optfdir/input_files/* ./

# If we start from a checkpoint we uncomment the cp command below and change the directory accordingly and we
# add --outfile [name_of_the_file.py] to the mpiexec call to the main script, e.g. --outfile=s1_mo_opt_out_2900.pkl
# cp -v $MPI_Run/outputs/s1_mo_opt_out_2900.pkl ./

mpiexec -n 32 python3 main_s1_mo.py --prior=$MPI_Run/prior_classic_scenario1w.txt --loss=$MPI_Run/s1_newscores/loss_model_global_mpi_s1_mo.py --algo=tpe --iter=100


cp -v ./s1_mo_opt_out_* $MPI_Run/outputs
