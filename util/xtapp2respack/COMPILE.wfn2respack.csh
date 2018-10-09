#!/bin/sh -x
#intel
ifort -qopenmp -O3 -o wfn2respack wfn2respack.F90 
#GNU 
#gfortran -fopenmp -O2 -o wfn2respack wfn2respack.F90 

