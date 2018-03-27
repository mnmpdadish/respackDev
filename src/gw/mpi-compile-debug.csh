#!/bin/sh -x
mpif90 -O2 -openmp -debug full -check all -CB -warn all -traceback -g -c m_rd_input.f90
mpif90 -O2 -openmp -debug full -check all -CB -warn all -traceback -g -c m_rd_dat_wfn.f90
mpif90 -O2 -openmp -debug full -check all -CB -warn all -traceback -g -c m_rd_dat_wan.f90
mpif90 -O2 -openmp -debug full -check all -CB -warn all -traceback -g -c m_rd_dat_eps.f90
mpif90 -O2 -openmp -debug full -check all -CB -warn all -traceback -g -c m_fft3d_20150826.f90 
mpif90 -O2 -openmp -debug full -check all -CB -warn all -traceback -g -c sub_band_disp.f90 
mpif90 -O2 -openmp -debug full -check all -CB -warn all -traceback -g -c util.f90  
mpif90 -O2 -openmp -debug full -check all -CB -warn all -traceback -g -c makekpts.f90 
mpif90 -O2 -openmp -debug full -check all -CB -warn all -traceback -g -c search_Rmin.f90 
mpif90 -O2 -openmp -debug full -check all -CB -warn all -traceback -g -c sub_wfn.f90 
mpif90 -O2 -openmp -debug full -check all -CB -warn all -traceback -g -c sub_eps.f90 
mpif90 -O2 -openmp -debug full -check all -CB -warn all -traceback -g -c sub_gw.f90 
mpif90 -O2 -openmp -debug full -check all -CB -warn all -traceback -g -c gw.f90 
mpif90 -O2 -openmp -debug full -check all -CB -warn all -traceback -g -o calc_gw gw.o sub_gw.o sub_wfn.o sub_eps.o sub_band_disp.o makekpts.o search_Rmin.o m_fft3d_20150826.o m_rd_input.o m_rd_dat_wfn.o m_rd_dat_wan.o m_rd_dat_eps.o util.o -lmkl_intel_lp64 -Wl,--start-group -lmkl_intel_thread -lmkl_core -Wl,--end-group -liomp5 -lpthread  
