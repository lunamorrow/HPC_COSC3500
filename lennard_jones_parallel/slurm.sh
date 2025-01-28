#!/bin/bash
#SBATCH --job-name=luna_ljmodel       # Job name
#SBATCH --output=make.txt             # Slurm with make commands
#SBATCH --error=error.txt             # Standard error file
#SBATCH --partition=cosc              # Partition or queue name
#SBATCH --nodes=4                     # Number of nodes
#SBATCH --ntasks-per-node=1           # Number of tasks per node
#SBATCH --cpus-per-task=4             # Number of CPU cores per task
#SBATCH --mem-per-cpu=4G
#SBATCH --time=0:30:00                # Maximum runtime (D-HH:MM:SS)
#SBATCH --nodelist=smp-9-16           # Run on same node each time for consistency

module add mpi/openmpi-x86_64
export PATH=/usr/lib64/openmpi/bin:$PATH

# make clean
# make
# mpiexec -n 4 ./ljmodel 16 10 > slurms/mpi16.txt
# mpiexec -n 4 ./ljmodel 25 12.5 > slurms/mpi25.txt
# mpiexec -n 4 ./ljmodel 36 15 > slurms/mpi36.txt
# mpiexec -n 4 ./ljmodel 49 17.5 > slurms/mpi49.txt
# mpiexec -n 4 ./ljmodel 64 20 > slurms/mpi64.txt
# mpiexec -n 4 ./ljmodel 81 22.5 > slurms/mpi81.txt
# mpiexec -n 4 ./ljmodel 100 25 > slurms/mpi100.txt
# mpiexec -n 4 ./ljmodel 225 37.5 > slurms/mpi225.txt
# mpiexec -n 4 ./ljmodel 400 50 > slurms/mpi400.txt
# mpiexec -n 4 ./ljmodel 625 62.5 > slurms/mpi625.txt
# mpiexec -n 4 ./ljmodel 900 75 > slurms/mpi900.txt
# mpiexec -n 4 ./ljmodel 1225 87.5 > slurms/mpi1225.txt
# mpiexec -n 4 ./ljmodel 1600 100 > slurms/mpi1600.txt

make clean
make
mpiexec -n 4 -map-by node -bind-to none ./ljmodel 225 37.5 > slurms/mpi225.txt
mpiexec -n 4 -map-by node -bind-to none ./ljmodel 400 50 > slurms/mpi400.txt
mpiexec -n 4 -map-by node -bind-to none ./ljmodel 625 62.5 > slurms/mpi625.txt
mpiexec -n 4 -map-by node -bind-to none ./ljmodel 900 75 > slurms/mpi900.txt
mpiexec -n 4 -map-by node -bind-to none ./ljmodel 1225 87.5 > slurms/mpi1225.txt
mpiexec -n 4 -map-by node -bind-to none ./ljmodel 1600 100 > slurms/mpi1600.txt
mpiexec -n 4 -map-by node -bind-to none ./ljmodel 2500 125 > slurms/mpi2500.txt
mpiexec -n 4 -map-by node -bind-to none ./ljmodel 3600 150 > slurms/mpi3600.txt