#!/bin/sh 

set -x 

#xtapp
#rm -f fort.* TiO2.lpt TiO2.rho TiO2.str TiO2.wfn *.log
#rm -f fort.10 ; ln -s TiO2.cg fort.10
#mpirun -np 1 inipot > log.TiO2-inipot 
#mpirun -np 1 cgmrpt > log.TiO2-cgmrpt 
#rm -f fort.10 ; ln -s TiO2.vb fort.10
#mpirun -np 1 inipot > log.TiO2-inipot-vb 
#mpirun -np 1 vbpef  > log.TiO2-vbpef 
#vbpef2gp-lsda TiO2.band 

#interface 
#./xtapp2respack.sh -b ./wfn2respack -s ./strconv TiO2 

#respack 
#mpirun -np 1 calc_wannier < input.in > log.TiO2-wannier

