#!/bin/sh -x
ifort -O2 -qopenmp -xHost -c gen-def-corr.f90  
ifort -O2 -qopenmp -xHost -c diagV.F90  
ifort -O2 -qopenmp -xHost -o a.out gen-def-corr.o diagV.o -lmkl_intel_lp64 -Wl,--start-group -lmkl_intel_thread -lmkl_core -Wl,--end-group -liomp5 -lpthread 
