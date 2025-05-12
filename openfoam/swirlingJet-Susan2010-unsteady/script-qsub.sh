#!/bin/bash

#SBATCH --mail-user=artur.gesla@epfl.ch

#SBATCH --mail-type=ALL

#SBATCH --job-name="test"

#SBATCH --mem 4000

#SBATCH --chdir /scratch/gesla/openfoam-work/swirlingJet-Susan2010-unsteady/

#SBATCH --ntasks-per-node 1

#SBATCH --cpus-per-task 1

#SBATCH --time 10:0:0

#SBATCH -A head

 

module purge

module load gcc

module load openmpi

module load python

module load openfoam-org
 

echo $1 > SESSION.NAME

echo `pwd`'/' >> SESSION.NAME

 

echo "SESSION.NAME --------------------------"

cat SESSION.NAME

echo "---------------------------------------"

 

. clean.sh ; 
blockMesh; 
srun pimpleFoam;

#srun ./nek5000 $1 > output.txt
#srun simpleFoam $1 > output.txt



