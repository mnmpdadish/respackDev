#!/bin/sh 

set -x 

#xtapp
#rm -f fort.* Al.lpt Al.rho Al.str Al.wfn *.log
#rm -f fort.10 ; ln -s Al.cg fort.10
#mpirun -np 1 inipot > log.Al-inipot 
#mpirun -np 1 cgmrpt > log.Al-cgmrpt 
#rm -f fort.10 ; ln -s Al.vb fort.10
#mpirun -np 1 inipot > log.Al-inipot-vb 
#mpirun -np 1 vbpef > log.Al-vbpef 
#vbpef2gp-lsda Al.band 

#interface 
#./xtapp2respack.sh -b ./wfn2respack -s ./strconv Al 

#respack 
#mpirun -np 1 calc_wannier < input.in > log.Al-wannier
#mpirun -np 1 calc_chiqw < input.in > log.Al-chiqw
#mpirun -np 1 calc_w3d < input.in > log.Al-calc_w3d
#mpirun -np 1 calc_j3d < input.in > log.Al-calc_j3d

