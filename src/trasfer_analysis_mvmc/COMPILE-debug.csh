#!/bin/sh -x
ifort -O2 -qopenmp -xHost -debug full -check all -warn all -traceback -g -c m_rd_dat_mvmc.f90  
ifort -O2 -qopenmp -xHost -c m_tetrahedron_20170325.F 
ifort -O2 -qopenmp -xHost -debug full -check all -warn all -traceback -g -c calc_dos.f90  
ifort -O2 -qopenmp -xHost -debug full -check all -warn all -traceback -g -c make_grd.f90  
ifort -O2 -qopenmp -xHost -debug full -check all -warn all -traceback -g -c sub_search_Rmin.f90  
ifort -O2 -qopenmp -xHost -debug full -check all -warn all -traceback -g -c est_nkbi.f90 
ifort -O2 -qopenmp -xHost -debug full -check all -warn all -traceback -g -c truncation.f90  
ifort -O2 -qopenmp -xHost -debug full -check all -warn all -traceback -g -c make_eig.f90  
ifort -O2 -qopenmp -xHost -debug full -check all -warn all -traceback -g -c make_kdata.f90  
ifort -O2 -qopenmp -xHost -debug full -check all -warn all -traceback -g -c sub_eigenvalue.f90  
ifort -O2 -qopenmp -xHost -debug full -check all -warn all -traceback -g -c sub_mkidx.f90  
ifort -O2 -qopenmp -xHost -debug full -check all -warn all -traceback -g -c util.f90  
ifort -O2 -qopenmp -xHost -debug full -check all -warn all -traceback -g -c est_ef.f90  
ifort -O2 -qopenmp -xHost -debug full -check all -warn all -traceback -g -c wrt.f90  
ifort -O2 -qopenmp -xHost -debug full -check all -warn all -traceback -g -c transfer_analysis.f90  
ifort -O2 -qopenmp -xHost -debug full -check all -warn all -traceback -g -o transfer_analysis transfer_analysis.o m_rd_dat_mvmc.o est_nkbi.o truncation.o make_eig.o make_kdata.o sub_eigenvalue.o util.o make_grd.o sub_search_Rmin.o m_tetrahedron_20170325.o calc_dos.o sub_mkidx.o est_ef.o wrt.o -lmkl_intel_lp64 -Wl,--start-group -lmkl_intel_thread -lmkl_core -Wl,--end-group -liomp5 -lpthread 
