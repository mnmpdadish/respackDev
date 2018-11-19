#!/bin/sh
natom=`grep 'nat' ./*.scf.in | awk -F'[=,]' '{print $2}'`
echo $natom > ./dat.atom_position
grep -A $natom 'ATOMIC_POSITIONS' $1 | tail -n$natom >> ./dat.atom_position
mv dat.atom_position ./dir-wfn/.

