#!/bin/sh -x
mpiifort -O2 -qopenmp -c m_rd_input.f90
mpiifort -O2 -qopenmp -c m_rd_dat_wfn.f90
mpiifort -O2 -qopenmp -c m_rd_dat_wan.f90
mpiifort -O2 -qopenmp -c m_rd_dat_eps.f90
mpiifort -O2 -qopenmp -c m_fft3d_20150826.f90 
mpiifort -O2 -qopenmp -c m_tetrahedron_20170325.F 
mpiifort -O2 -qopenmp -c sub_wfn.f90 
mpiifort -O2 -qopenmp -c sub_eps.f90 
mpiifort -O2 -qopenmp -c sub_gw.f90 
mpiifort -O2 -qopenmp -c util.f90  
mpiifort -O2 -qopenmp -c sub_band_disp.f90 
mpiifort -O2 -qopenmp -c make_kpts.f90 
mpiifort -O2 -qopenmp -c make_sgmw.f90 
mpiifort -O2 -qopenmp -c calc_gwdos.f90 
mpiifort -O2 -qopenmp -c calc_gwakw.f90 
mpiifort -O2 -qopenmp -c search_Rmin.f90 
mpiifort -O2 -qopenmp -c sub_ksdos.f90  
mpiifort -O2 -qopenmp -c sub_gwdos.f90  
mpiifort -O2 -qopenmp -c sub_eigenvalue.f90  
mpiifort -O2 -qopenmp -c sub_mkidx.f90  
mpiifort -O2 -qopenmp -c sub_det_shift.f90  
mpiifort -O2 -qopenmp -c gw.f90 
mpiifort -O2 -qopenmp -o calc_gw gw.o sub_gw.o sub_wfn.o sub_eps.o util.o sub_band_disp.o make_kpts.o search_Rmin.o make_sgmw.o calc_gwdos.o calc_gwakw.o m_fft3d_20150826.o m_rd_input.o m_rd_dat_wfn.o m_rd_dat_wan.o m_rd_dat_eps.o sub_gwdos.o sub_ksdos.o sub_det_shift.o sub_mkidx.o sub_eigenvalue.o m_tetrahedron_20170325.o -lmkl_intel_lp64 -Wl,--start-group -lmkl_intel_thread -lmkl_core -Wl,--end-group -liomp5 -lpthread  


