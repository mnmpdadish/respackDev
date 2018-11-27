PROGRAM main 
  use m_rd_dat_zvo 
  use m_rd_transdef  
  include "config.h" 
  !
  !read input from command line 
  !
  !(default for zvo) 
  !delt=0.01d0     !Greens function delt in eV
  !dmna=0.001d0    !Ttrhdrn parameter dmna in eV
  !dmnr=0.001d0    !Ttrhdrn parameter dmnr in eV
  !delw=2.0d0*delt !Grid width in eV
  !flwe=0          !Flg whether calculate weighted transfers (0:not calc, 1:calc)
  !thtr=0.0d0      !Threshold for transfer integral
  !elnm=0.0d0      !Total number of electrons in unitcell
  !
  !(default for ztr) 
  !delw=10.0d0*delt !Grid width in eV
  !
  ncount=iargc() 
  allocate(real_arg(ncount)); real_arg=0.0d0 
  do i=1,ncount-2 
   call getarg(i,arg) 
   read(arg,*) real_arg(i)  
  enddo 
  !zvo
  call getarg(ncount-1,arg)
  read(arg,*) zvo 
  !ztr
  call getarg(ncount,arg)
  read(arg,*) ztr 
  !
  !write(6,*) ncount
  !write(6,*) real_arg(1)
  !write(6,*) real_arg(2)
  !write(6,*) real_arg(3)
  !write(6,*) real_arg(4)
  !write(6,*) real_arg(5)
  !write(6,*) real_arg(6)
  !write(6,*) real_arg(7)
  !
  delt=real_arg(1)
  dmna=real_arg(2)  
  dmnr=real_arg(3)  
  delw=real_arg(4)  
  flwe=real_arg(5)  
  thtr=real_arg(6)  
  elnm=real_arg(7)  
  !
  delt=delt/au !au <- eV
  dmna=dmna/au !au <- eV
  dmnr=dmnr/au !au <- eV
  delw=delw/au !au <- eV
  flg_weight=nint(flwe)  
  if(flg_weight/=0) flg_weight=1
  threshold_transfer=thtr/au !au <- eV
  electron_number=elnm 
  !
  if(ztr)then
    write(6,*) 
    write(6,'(a50)')'##### TRANSFER ANALYSIS (trans.def) #####'
    write(6,*) 
    !
    !read trans.def 
    !
    call rd_transdef  
    !
    !diagonalize TR 
    !
    allocate(EIG_TR(Ndim_TR)); EIG_TR=0.0d0 
    call diagN(Ndim_TR,TR(1,1),EIG_TR(1)) 
    !
    !make dos-grid
    !
    call est_ndosgrd(Ndim_TR,1,EIG_TR(1),delw,emin,emax,ndosgrd)  
    !
    !calc histgram
    !
    emin_grd=nint(emin/delw) 
    emax_grd=nint(emax/delw)
    allocate(hist(emin_grd:emax_grd)); hist=0.0d0 
    call calc_hist(Ndim_TR,emin_grd,emax_grd,delw,EIG_TR(1),hist(emin_grd)) 
    !
    !wrt histgram
    ! 
    call wrt_hist(emin_grd,emax_grd,delw,hist(emin_grd)) 
    !
    write(6,'(a50)')'##### FINISH TRANSFER ANALYSIS (trans.def) #####'
    !
  endif 
  if(zvo)then 
    write(6,*) 
    write(6,'(a50)')'##### TRANSFER ANALYSIS (ZVO) #####'
    write(6,*) 
    write(6,'(a50,x,f10.5)')'Greens function delt (eV)=',delt*au
    write(6,'(a50,x,f10.5)')'Terahedon parameter dmna (eV)=',dmna*au
    write(6,'(a50,x,f10.5)')'Terahedon parameter dmnr (eV)=',dmnr*au
    write(6,'(a50,x,f10.5)')'Grid spacing delw (eV)=',delw*au
    write(6,'(a50,x,i10)')'Use weigted transfer (0:not)=',flg_weight 
    write(6,'(a50,x,f10.5)')'Threshold for transfer (eV)=',threshold_transfer*au
    write(6,'(a50,x,f10.5)')'Electron numbers in unit cell=',electron_number  
    write(6,*) 
    !
    !read zvo(mvmc files) 
    !
    call rd_dat_mkkpts 
    call rd_dat_hr 
    call rd_dat_geom 
    call rd_dat_bandkpts 
    !
    !create dir-tr
    !
    call system('rm -rf dir-tr') 
    call system('mkdir dir-tr') 
    !
    !truncate HR
    !
    call truncation(NWF,Na1,Na2,Na3,threshold_transfer,HR(1,1,-Na1,-Na2,-Na3))
    !
    !make EKS
    !
    allocate(EKS(NWF,NTK)); EKS=0.0d0 
    allocate(VKS(NWF,NWF,NTK)); VKS=0.0d0 
    call make_eig(NWF,NTK,Na1,Na2,Na3,nkb1,nkb2,nkb3,flg_weight,a1(1),a2(1),a3(1),SK0(1,1),HR(1,1,-Na1,-Na2,-Na3),EKS(1,1),VKS(1,1,1)) 
    !
    !make dos-grid
    !
    call est_ndosgrd(NWF,NTK,EKS(1,1),delw,emin,emax,ndosgrd)  
    allocate(dosgrd(ndosgrd));dosgrd=0.0d0 
    call make_dosgrd(emin,delw,ndosgrd,dosgrd(1))  
    !
    !calc dos 
    !
    allocate(dos(ndosgrd));dos=0.0d0 
    call calc_dos(NWF,NTK,nkb1,nkb2,nkb3,ndosgrd,dosgrd(1),EKS(1,1),SK0(1,1),delt,dmnr,dmna,b1(1),b2(1),b3(1),dos(1)) 
    !
    !estimate ef 
    !
    call est_ef(ndosgrd,delw,electron_number,dosgrd(1),dos(1),FermiEnergy) 
    allocate(efline(ndosgrd));efline=0.0d0 
    call make_efline(ndosgrd,FermiEnergy,delw,dosgrd(1),dos(1),efline(1)) 
    !
    !wrt dat.dos 
    !
    call wrt_dos(threshold_transfer,ndosgrd,dosgrd(1),dos(1),efline(1)) 
    !
    deallocate(EKS,VKS) 
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
    call wrt_iband(threshold_transfer,NWF,NSK_BAND_DISP,kdata(1),EKS(1,1)) 
    !
    deallocate(EKS,VKS) 
    !
    write(6,'(a50)')'##### FINISH TRANSFER ANALYSIS (ZVO) #####'
  endif!zvo 
  stop
end     
