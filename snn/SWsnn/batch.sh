#!/bin/bash
#SBATCH -J SWsnn
#SBATCH -p dongsheng
#SBATCH -N 4
#SBATCH --ntasks-per-node 4
#SBATCH --gres=dcu:4
date
mpirun -n 16 ./swtest
date
