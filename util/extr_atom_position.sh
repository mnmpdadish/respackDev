#!/bin/sh
natom=`grep 'nat' ./*.scf.in | awk '{print $3}'` 
echo $natom > ./dat.atom_position 
grep -A $natom 'ATOMIC_POSITIONS' ./*.scf.in | tail -$natom >> ./dat.atom_position 

