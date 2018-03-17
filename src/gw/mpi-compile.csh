#!/bin/sh -x
mpiifort -O2 -openmp  -c m_rdinput.f90
mpiifort -O2 -openmp  -c m_rd_dat_wfn.f90
mpiifort -O2 -openmp  -c m_rd_dat_wan.f90
mpiifort -O2 -openmp  -c m_rd_dat_eps.f90
mpiifort -O2 -openmp  -c m_fft3d_20150826.f90 
mpiifort -O2 -openmp  -c sub_band_disp.f90 
mpiifort -O2 -openmp  -c est_latparam.f90
mpiifort -O2 -openmp  -c est_nwx2.f90 
mpiifort -O2 -openmp  -c inv.f90  
mpiifort -O2 -openmp  -c invmat_complex.f90
mpiifort -O2 -openmp  -c diagV.f90  
mpiifort -O2 -openmp  -c makekpts.f90 
mpiifort -O2 -openmp  -c search_Rmin.f90 
mpiifort -O2 -openmp  -c sub_wfn.f90 
mpiifort -O2 -openmp  -c sub_eps.f90 
mpiifort -O2 -openmp  -c sub_vxc.f90 
mpiifort -O2 -openmp  -c sub_sx.f90 
mpiifort -O2 -openmp  -c gw.f90 
mpiifort -O2 -openmp  -o calc_gw gw.o sub_vxc.o sub_wfn.o sub_eps.o sub_sx.o sub_band_disp.o est_nwx2.o est_latparam.o diagV.o makekpts.o search_Rmin.o m_fft3d_20150826.o m_rdinput.o m_rd_dat_wfn.o m_rd_dat_wan.o m_rd_dat_eps.o inv.o invmat_complex.o -lmkl_intel_lp64 -Wl,--start-group -lmkl_intel_thread -lmkl_core -Wl,--end-group -liomp5 -lpthread  
