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
parser.add_argument('--bnd',action='store_true',help='BaND dispersion calc for zvo data') 
parser.add_argument('--frm',action='store_true',help='FeRMi surface calc for fine k-mesh') 
parser.add_argument('--his',action='store_true',help='HIStgram analysis for ztrans.def') 
parser.add_argument('--gd',type=float,default=0.01,help='Greens function Delt in eV')
parser.add_argument('--we',type=int,default=0,help='Flg whether calculate Weighted transfers (0:not calc, 1:calc)') 
parser.add_argument('--th',type=float,default=0.0,help='THreshold for transfer integral')
parser.add_argument('--el',type=float,default=0.0,help='total number of ELectrons in unitcell') 
parser.add_argument('--kgd',default='0 0 0',help='k grid') 
#
args = parser.parse_args()
#
#print args.dos  
#print args.bnd
#print args.frm 
#print args.his  
#print args.gd
#print args.we
#print args.th
#print args.el 
#print args.kgd 
#
kgd=args.kgd 
kgd=kgd.split() 
kgd=[int(i) for i in kgd] 
#
cmd=["./transfer_analysis", str(args.gd), str(args.we), str(args.th), str(args.el) \
                          , str(kgd[0]), str(kgd[1]), str(kgd[2]) \
                          , str(args.dos), str(args.bnd), str(args.frm), str(args.his)]
subprocess.Popen(cmd)
