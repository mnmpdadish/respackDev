#!/bin/sh 

set -x 

#xtapp
#rm -f fort.* La2CuO4.lpt La2CuO4.rho La2CuO4.str La2CuO4.wfn *.log
#rm -f fort.10 ; ln -s La2CuO4.cg fort.10
#mpirun -np 1 inipot > log.La2CuO4-inipot 
#mpirun -np 1 cgmrpt > log.La2CuO4-cgmrpt 
#rm -f fort.10 ; ln -s La2CuO4.vb fort.10
#mpirun -np 1 inipot > log.La2CuO4-inipot-vb 
#mpirun -np 1 vbpef  > log.La2CuO4-vbpef 
#vbpef2gp-lsda La2CuO4.band 

#interface 
#./xtapp2respack.sh -b ./wfn2respack -s ./strconv La2CuO4

#respack 
#mpirun -np 1 calc_wannier < input.in > log.La2CuO4-wannier
#mpirun -np 1 calc_chiqw < input.in   > log.La2CuO4-chiqw
#mpirun -np 1 calc_w3d < input.in     > log.La2CuO4-calc_w3d
#mpirun -np 1 calc_j3d < input.in     > log.La2CuO4-calc_j3d

