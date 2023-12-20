#!/bin/bash
#SBATCH --time=24:00:00
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

cp -v /home/chargaut/scratch/MPI_Run/s3_newscores/main_scores_mo.py ./
cp -v /home/chargaut/scratch/MPI_Run/fenv.py ./
cp -v /home/chargaut/scratch/MPI_Run/fscore.py ./
cp -v /home/chargaut/scratch/MPI_Run/fturbation.py ./
cp -v /home/chargaut/scratch/MPI_Run/s3_newscores/algo_scores_mo.py ./
cp -v /home/chargaut/scratch/MPI_Run/s3_newscores/loss_model_global_mpi_scores_mo.py ./
cp -v /home/chargaut/scratch/MPI_Run/s3_newscores/func_model_scores.py ./
cp -v /home/chargaut/scratch/MPI_Run/MPI_array_allocation.py ./
cp -v /home/chargaut/scratch/formating_global_files/pre_split_narval/* ./
cp -v /home/chargaut/scratch/MPI_Run/outputs/article_results/s3_mo_opt_out_2900.pkl ./

#mpiexec -n 64 python3 pre_spliting_files.py
mpiexec -n 32 python3 main_scores_mo.py --prior=/home/chargaut/scratch/MPI_Run/prior_classic_scenario3.txt --loss=/home/chargaut/scratch/MPI_Run/loss_model_global_mpi_scores_mo.py --algo=tpe


cp -v ./s3_mo_opt_out* /home/chargaut/scratch/MPI_Run/outputs/article_results
