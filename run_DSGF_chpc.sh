#!/bin/bash
#SBATCH -A francoeur
#SBATCH -p notchpeak
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -t 24:00:00
#SBATCH -C rom

# declare STRING variable
TITLE="Near-field radiative heat transfer between thermal objects"
METHOD="DSGF: Discrete System Green Function method" 
VERSION="Developed by RETL Lab at the Department of Mechanical Engineering, The University of Utah, UT, USA"

# Use this command if permission is denied: sudo chmod u+x run_dsgf.sh
# Use this command if bad interpreter appears: sed -i -e 's/\r$//' run_dsgf.sh

#print variable on a screen
echo $TITLE
echo $METHOD
echo $VERSION
echo "Compilation in process..."
#Cuma  

# Insert these command lines when open a new terminal to use icc compiler
module load intel-oneapi-compilers/2022.2.1 intel-oneapi-mkl/2022.2.1
# source /opt/intel/oneapi/mkl/2022.0.2/env/vars.sh 
# source /opt/intel/oneapi/compiler/2022.0.2/env/vars.sh 
# The compilers and mkl were installed from #https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit-download.html?operatingsystem=linux&distributions=webdownload&options=offline

# export OMP_THREADS=4 # use multiple cores during the calculation
export OMP_THREADS=$SLURM_NTASKS

icc -Ofast -qopenmp -Wall -std=c99 -I include/ src/computational/GreensFunction.c src/array_functions.c src/geometry/sample.c src/geometry/user_defined.c src/material/u_SiC.c src/material/SiN.c src/material/SiO2.c src/material/SiC.c src/DSGF_main.c src/file_utils.c src/functions_DSGF.c src/computational/solvers/iterative_solver.c src/computational/solvers/direct_solver.c src/computational/ThermalPower.c src/geometry/shared.c -o DSGF -mcmodel=medium  -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl  -I"${MKLROOT}/include" # -openmp 
 
# old geometry files src/geometry/sphere.c src/geometry/thin_film.c

echo "Compilation complete"
#  MALLOC_sgf-sphere-lapacke-linux.c  sgf-flat-surfaces-lapacke-linux.c  sgf-lapacke_original.c or sgf-lapacke-linux.c sgf-sphere-lapacke-linux.c sgf-flat-surfaces-lapacke-linux_v2.c  MALLOC_sgf-sphere-lapacke-linux    

# Improvement MALLOC:
# https://stackoverflow.com/questions/13534966/how-to-dynamically-allocate-a-contiguous-block-of-memory-for-a-2d-array
# https://stackoverflow.com/questions/39108092/allocating-contiguous-memory-for-a-3d-array-in-c

# OpenMP: https://cvw.cac.cornell.edu/OpenMP/default
# Intel advisor: https://www.intel.com/content/www/us/en/develop/documentation/get-started-with-advisor/top/prototype-threading-designs.html#prototype-threading-designs

echo "Starting the run on node" `hostname`
./DSGF

