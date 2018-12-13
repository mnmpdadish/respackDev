#!/usr/bin/env python
#coding: utf-8
#(c) Kazuma Nakamura
import sys
import os
import subprocess 
import numpy as np
import argparse
#
parser=argparse.ArgumentParser() 
parser.add_argument('--zvo',action='store_true',help='transfer analysis for ZVO data') 
parser.add_argument('--ztr',action='store_true',help='transfer analysis for ZTRans.def') 
parser.add_argument('--gd',type=float,default=0.01,help='Greens function Delt in eV')
parser.add_argument('--ta',type=float,default=0.001,help='Ttrhdrn parameter dmnA in eV')
parser.add_argument('--tr',type=float,default=0.001,help='Ttrhdrn parameter dmnR in eV')
parser.add_argument('--gs',type=float,default=0.02,help='Grid Spacing in eV')
parser.add_argument('--we',type=int,default=0,help='Flg whether calculate Weighted transfers (0:not calc, 1:calc)') 
parser.add_argument('--th',type=float,default=0.0,help='THreshold for transfer integral')
parser.add_argument('--el',type=float,default=0.0,help='total number of ELectrons in unitcell') 
#
args = parser.parse_args()
#
#print args.zvo  
#print args.ztr  
#print args.gd
#print args.ta
#print args.tr
#print args.gs
#print args.we
#print args.th
#print args.el 
#
zvo=str(args.zvo)
ztr=str(args.ztr)
delt=str(args.gd)
ttra=str(args.ta)
ttrr=str(args.tr)
#
if(args.zvo): 
    delw=str(args.gs)
if(args.ztr):
    delw=str(10.0*args.gs) 
#
flwe=str(args.we)
thtr=str(args.th)
elnm=str(args.el)
#
cmd=["./transfer_analysis",delt,ttra,ttrr,delw,flwe,thtr,elnm,zvo,ztr]
subprocess.Popen(cmd)
