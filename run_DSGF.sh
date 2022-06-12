#!/bin/bash
# declare STRING variable
TITLE="Near-field radiative heat transfer between thermal objects"
METHOD="DSGF: Discrete System Green Function method" 
VERSION="Developed by RETL Lab at the Department of Mechanical Engineering, The University of Utah, UT, USA"

#print variable on a screen
echo $TITLE
echo $METHOD
echo $VERSION
echo "Compilation in process..."
#Cuma  

# Insert these command lines when open a new terminal to use icc compiler
# source /opt/intel/oneapi/mkl/2022.0.2/env/vars.sh 
# source /opt/intel/oneapi/compiler/2022.0.2/env/vars.sh 
# The compilers and mkl were installed from #https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit-download.html?operatingsystem=linux&distributions=webdownload&options=offline

export OMP_THREADS=4 # use multiple cores during the calculation

icc -O3  -std=c99 DSGF_main.c -o DSGF -mcmodel=medium  -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl  -I"${MKLROOT}/include" # -openmp 

echo "Compilation complete"
#  MALLOC_sgf-sphere-lapacke-linux.c  sgf-flat-surfaces-lapacke-linux.c  sgf-lapacke_original.c or sgf-lapacke-linux.c sgf-sphere-lapacke-linux.c sgf-flat-surfaces-lapacke-linux_v2.c  MALLOC_sgf-sphere-lapacke-linux    

time ./DSGF  #time function displays the real, user and system time for running the code;    /.dsgf_LAPACKE_9N2 runs the compiled code. 

# Improvement MALLOC:
# https://stackoverflow.com/questions/13534966/how-to-dynamically-allocate-a-contiguous-block-of-memory-for-a-2d-array
# https://stackoverflow.com/questions/39108092/allocating-contiguous-memory-for-a-3d-array-in-c

# OpenMP: https://cvw.cac.cornell.edu/OpenMP/default
# Intel advisor: https://www.intel.com/content/www/us/en/develop/documentation/get-started-with-advisor/top/prototype-threading-designs.html#prototype-threading-designs

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#old and slow compilation:

#gcc -O3 ../lapack-3.10.0/LAPACKE/include -c -std=c99 dsgf.c -o dsgf_LAPACKE_9N2.o -mcmodel=medium # sgf-flat-surfaces-lapacke-linux.c  sgf-lapacke_original.c or sgf-lapacke-linux.c sgf-sphere-lapacke-linux.c sgf-flat-surfaces-lapacke-linux_v2.c

#gfortran -O2 -frecursive -mcmodel=medium -o dsgf_LAPACKE_9N2 dsgf_LAPACKE_9N2.o ../lapack-3.10.0/liblapacke.a ../lapack-3.10.0/liblapack.a ../lapack-3.10.0/librefblas.a

#./dsgf_LAPACKE_9N2
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#https://community.intel.com/t5/Intel-Fortran-Compiler/Fortran-additional-relocation-overflows-error-while-reading-big/m-p/1162453 
#https://community.intel.com/t5/Intel-Fortran-Compiler/relocation-truncated-to-fit-R-X86-64-PC32/td-p/762814
#-mcmodel=medium

#https://stackoverflow.com/questions/12916176/gfortran-for-dummies-what-does-mcmodel-medium-do-exactly  do not use -mcmodel=large

