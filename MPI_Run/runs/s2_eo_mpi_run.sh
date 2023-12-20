#!/bin/bash
#SBATCH --time=20:00:00
#SBATCH --account=def-sonol
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem=249G
#SBATCH --mail-user=charles.gauthier.1@umontreal.ca
#SBATCH --mail-type=ALL


module load python/3.8
module load scipy-stack
module load mpi4py

cd $SLURM_TMPDIR/

source ~/ENV_2/bin/activate

cp -v /home/chargaut/scratch/MPI_Run/s2_newscores/main_s2_eo.py ./
cp -v /home/chargaut/scratch/MPI_Run/fenv.py ./
cp -v /home/chargaut/scratch/MPI_Run/fscore.py ./
cp -v /home/chargaut/scratch/MPI_Run/fturbation.py ./
cp -v /home/chargaut/scratch/MPI_Run/s2_newscores/algo_s2_eo.py ./
cp -v /home/chargaut/scratch/MPI_Run/s2_newscores/loss_model_global_mpi_s2_eo.py ./
cp -v /home/chargaut/scratch/MPI_Run/s2_newscores/func_model_s2.py ./
cp -v /home/chargaut/scratch/MPI_Run/MPI_array_allocation.py ./
cp -v /home/chargaut/scratch/formating_global_files/pre_split_narval/* ./
cp -v /home/chargaut/scratch/MPI_Run/outputs/article_results/s2_eo_opt_out_2900.pkl ./

#mpiexec -n 48 python3 pre_spliting_files.py
mpiexec -n 32 python3 main_s2_eo.py --prior=/home/chargaut/scratch/MPI_Run/prior_classic_scenario2w.txt --loss=/home/chargaut/scratch/MPI_Run/s2_newscores/loss_model_global_mpi_s2_eo.py --algo=tpe


cp -v ./s2_eo_opt_out_* /home/chargaut/scratch/MPI_Run/outputs/article_results

