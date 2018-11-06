PROGRAM main 
  use m_rd_dat_mvmc 
  include "config.h" 
  !
  !default 
  !
  delt=0.01d0 !ttrhdrn Green's function delt (eV)
  dmna=0.001d0 !ttrhdrn (eV)
  dmnr=0.001d0 !ttrhdrn (eV)
  deltw=2.0d0*delt !grid_spacing (eV)
  flg_weight=.false. !flg to calculate weights (eV)
  threshold_transfer=0.0d0 !threshold for transfer (eV)
  electron_number=0.0d0 !total number of electrons in unitcell
  !
  !read input from command line 
  !
  ncount=iargc() 
  allocate(real_arg(ncount)); real_arg=0.0d0 
  do i=1,ncount 
   call getarg(i,arg) 
   read(arg,*) real_arg(ncount)  
  enddo 
  threshold_transfer=real_arg(1)
  electron_number=real_arg(2)  
  !
  delt =delt/au !au <- eV
  dmna =dmna/au !au <- eV
  dmnr =dmnr/au !au <- eV
  deltw=deltw/au!au <- eV
  threshold_transfer=threshold_transfer/au!au <- eV
  write(6,'(a33,x,f10.5)')'Greens function delt (eV)=',delt*au
  write(6,'(a33,x,f10.5)')'terahedon parameter dmna (eV)=',dmna*au
  write(6,'(a33,x,f10.5)')'terahedon parameter dmnr (eV)=',dmnr*au
  write(6,'(a33,x,f10.5)')'Grid spacing deltw (eV)=',deltw*au
  write(6,'(a33,x,f10.5)')'threshold transfer (eV)=',threshold_transfer*au 
  write(6,'(a33,x,f10.5)')'electron numbers in unit cell=',electron_number  
  write(6,*) 
  !
  !read zvo(mvmc files) 
  !
  call rd_dat_hr 
  call rd_dat_geom 
  call rd_dat_bandkpts 
  call rd_dat_mkkpts 
  !
  !###########################
  !  WANNIER-DOS CALCULATION
  !###########################
  !
  !create dir-tr
  !
  call system('rm -rf dir-tr') 
  call system('mkdir dir-tr') 
  !
  !truncate HMATR
  !
  call truncation(NWF,Na1,Na2,Na3,threshold_transfer,HR(1,1,-Na1,-Na2,-Na3))
  !
  !estimate nbb1,nkb2,nkb3 
  !
  call est_nkbi(NTK,SK0,nkb1,nkb2,nkb3)  
  !
  !make EKS
  !
  allocate(EKS(NWF,NTK)); EKS=0.0d0 
  allocate(VKS(NWF,NWF,NTK)); VKS=0.0d0 
  call make_eig(NWF,NTK,Na1,Na2,Na3,nkb1,nkb2,nkb3,flg_weight,a1(1),a2(1),a3(1),SK0(1,1),HR(1,1,-Na1,-Na2,-Na3),EKS(1,1),VKS(1,1,1)) 
  !
  !make dos-grid
  !
  call est_ndosgrd(NWF,NTK,EKS(1,1),deltw,emin,emax,ndosgrd)  
  allocate(dosgrd(ndosgrd));dosgrd=0.0d0 
  call make_dosgrd(emin,deltw,ndosgrd,dosgrd(1))  
  !
  !calc dos 
  !
  allocate(dos(ndosgrd));dos=0.0d0 
  call calc_dos(NWF,NTK,nkb1,nkb2,nkb3,ndosgrd,dosgrd(1),EKS(1,1),SK0(1,1),delt,dmnr,dmna,b1(1),b2(1),b3(1),dos(1)) 
  !
  !estimate ef 
  !
  call est_ef(ndosgrd,deltw,electron_number,dosgrd(1),dos(1),FermiEnergy) 
  allocate(efline(ndosgrd));efline=0.0d0 
  call make_efline(ndosgrd,FermiEnergy,deltw,dosgrd(1),dos(1),efline(1)) 
  !
  !wrt dat.dos 
  !
  call wrt_dos(ndosgrd,dosgrd(1),dos(1),efline(1)) 
  !
  deallocate(EKS,VKS) 
  !
  !###############################
  !  BAND DISPERSION CALCULATION
  !###############################
  !
  !make EKS
  !
  allocate(EKS(NWF,NSK_BAND_DISP)); EKS=0.0d0 
  allocate(VKS(NWF,NWF,NSK_BAND_DISP)); VKS=0.0d0 
  call make_eig(NWF,NSK_BAND_DISP,Na1,Na2,Na3,nkb1,nkb2,nkb3,flg_weight,a1(1),a2(1),a3(1),SK_BAND_DISP(1,1),HR(1,1,-Na1,-Na2,-Na3),&
  EKS(1,1),VKS(1,1,1)) 
  !
  !make kdata 
  !
  allocate(kdata(NSK_BAND_DISP)); kdata=0.0d0 
  call make_kdata(Ndiv,N_sym_points,NSK_BAND_DISP,SK_BAND_DISP(1,1),b1(1),b2(1),b3(1),kdata(1)) 
  !
  !wrt dat.iband 
  !
  call wrt_iband(NWF,NSK_BAND_DISP,kdata(1),EKS(1,1)) 
  !
  deallocate(EKS,VKS) 
  !
  stop
end     
