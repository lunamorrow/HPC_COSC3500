#!/bin/bash
#SBATCH --job-name=luna_ljmodel       # Job name
#SBATCH --output=make.txt
#SBATCH --error=error.txt             # Standard error file
#SBATCH --partition=cosc              # Partition or queue name
#SBATCH --nodes=1                     # Number of nodes
#SBATCH --ntasks-per-node=1           # Number of tasks per node
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem-per-cpu=10G
#SBATCH --time=1:00:00                # Maximum runtime (D-HH:MM:SS)
#SBATCH --nodelist=smp-9-17           # Run on same node each time for consistency

# --- Run tests on Lennard-Jones model ----
make clean
make
./ljmodel 9 7.5 > optimised1.txt
./ljmodel 16 10 > optimised2.txt
./ljmodel 25 12.5 > optimised3.txt
./ljmodel 36 15 > optimised4.txt
./ljmodel 49 17.5 > optimised5.txt
./ljmodel 64 20 > optimised6.txt
./ljmodel 81 22.5 > optimised7.txt
./ljmodel 100 25 > optimised8.txt
./ljmodel 60 20 > optimised60.txt

# --- Analyse breakdown of time taken with gprof ----
# gprof ./ljmodel gmon.out > analysis.txt
