#!/bin/bash


icc -debug -g DSGF_main.c file_utils.c functions_DSGF.c iterative_solver.c indexing_util.c debugging/debugging_utils.c -o DSGF_debug -mcmodel=medium  -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl  -I"${MKLROOT}/include"

gdb-oneapi DSGF_debug