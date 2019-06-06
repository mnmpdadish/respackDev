#!/bin/sh 

set -x 

#xtapp
#rm -f fort.* SrVO3.lpt SrVO3.rho SrVO3.str SrVO3.wfn *.log
#rm -f fort.10 ; ln -s SrVO3.cg fort.10
#mpirun -np 1 inipot > log.SrVO3-inipot 
#mpirun -np 1 cgmrpt > log.SrVO3-cgmrpt 
#rm -f fort.10 ; ln -s SrVO3.vb fort.10
#mpirun -np 1 inipot > log.SrVO3-inipot-vb 
#mpirun -np 1 vbpef  > log.SrVO3-vbpef 
#vbpef2gp-lsda SrVO3.band 

#interface 
#./xtapp2respack.sh -b ./wfn2respack -s ./strconv SrVO3

#respack 
#mpirun -np 1 calc_wannier < input.in > log.SrVO3-wannier
#mpirun -np 1 calc_chiqw < input.in   > log.SrVO3-chiqw
#mpirun -np 1 calc_w3d < input.in     > log.SrVO3-calc_w3d
#mpirun -np 1 calc_j3d < input.in     > log.SrVO3-calc_j3d

