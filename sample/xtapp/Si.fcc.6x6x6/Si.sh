#!/bin/sh 

set -x 

#xtapp
#rm -f fort.* Si.lpt Si.rho Si.str Si.wfn *.log
#rm -f fort.10 ; ln -s Si.cg fort.10
#mpirun -np 1 inipot > log.Si-inipot 
#mpirun -np 1 cgmrpt > log.Si-cgmrpt 
#rm -f fort.10 ; ln -s Si.vb fort.10
#mpirun -np 1 inipot > log.Si-inipot-vb 
#mpirun -np 1 vbpef  > log.Si-vbpef 
#vbpef2gp-lsda Si.band 

#interface 
#./xtapp2respack.sh -b ./wfn2respack -s ./strconv Si 

#respack 
#mpirun -np 1 calc_wannier < input.in > log.Si-wannier
#mpirun -np 1 calc_chiqw < input.in   > log.Si-chiqw
#mpirun -np 1 calc_w3d < input.in     > log.Si-calc_w3d
#mpirun -np 1 calc_j3d < input.in     > log.Si-calc_j3d

