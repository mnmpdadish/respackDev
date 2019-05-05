#coding: utf-8
#(c) Kazuma Nakamura
import sys
import os
import subprocess 
import numpy as np
import argparse
#
parser=argparse.ArgumentParser() 
parser.add_argument('--dos',action='store_true',help='DOS calc for zvo data') 
parser.add_argument('--bnd',action='store_true',help='Band dispersion calc for zvo data') 
parser.add_argument('--frm',action='store_true',help='Fermi surface calc for fine k-mesh') 
parser.add_argument('--his',action='store_true',help='Histgram analysis for ztrans.def') 
parser.add_argument('--delt',type=float,default=0.01,help='<< 0.01eV: DEFAULT >> Greens function delt in eV')
parser.add_argument('--elec',type=float,default=0.00,help='<< 0.00: DEFAULT >> Total number of electrons in unitcell') 
parser.add_argument('--ecut',type=float,default=0.00,help='<< 0.00eV: DEFAULT >> Energy cutoff for transfer in eV')
parser.add_argument('--rcut',type=float,default=100.00,help='<< 100.0AA: DEFAULT >> Distance cutoff for transfer in AA')
parser.add_argument('--diff',type=float,default=0.01,help='<< 0.01eV: DEFAULT>> Match threshold for two transfers in eV')
parser.add_argument('--kgrd',default='0 0 0',help='<< The same as RESPACK calculation: DEFAULT >> k grid') 
#
args = parser.parse_args()
#
#print args.dos  
#print args.bnd
#print args.frm 
#print args.his  
#print args.delt
#print args.elec 
#print args.ecut
#print args.rcut
#print args.diff 
#print args.kgrd 
#
kgrd=args.kgrd 
kgrd=kgrd.split() 
kgrd=[int(i) for i in kgrd] 
#
cmd=["./calc_tr", str(args.delt), str(args.elec), str(args.ecut), str(args.rcut), str(args.diff)\
                , str(kgrd[0]), str(kgrd[1]), str(kgrd[2])\
                , str(args.dos), str(args.bnd), str(args.frm), str(args.his)]
subprocess.Popen(cmd)
