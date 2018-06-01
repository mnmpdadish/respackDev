#!/bin/bash 
pw.x < GaAs.scf.in 
pw.x < GaAs.band.in 
bands.x < GaAs.bands.in 
./a.out GaAs.band
mkdir dir-wfn
qe2respack work/GaAs.save/
calc_wannier < input.in 
