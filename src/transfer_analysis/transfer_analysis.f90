!
!(default parameters) 
!
!delt=0.01d0          !Greens function delt in eV
!dmna=0.001d0         !Ttrhdrn parameter dmna in eV
!dmnr=0.001d0         !Ttrhdrn parameter dmnr in eV
!delw=2.0d0*delt      !Grid width in eV [dos]
!delw=10.0d0*delt     !Grid width in eV [hist]
!ecut=0.0d0           !Energy cutoff for transfer integral in eV
!rcut=100.0d0         !Distance cutoff for transfer integral in AA
!diff=0.01d0          !Match threshold for two transfer integral in eV 
!elnm=0.0d0           !Total number of electrons in unitcell
!kgrd='nkb1 nkb2 nkb3'!k grid 
!
PROGRAM main 
  use m_rd_dat_zvo 
  use m_rd_transdef  
  use m_truncation 
  use m_dos, only: calculate_dos 
  use m_eigenstate, only: calculate_eigenstate
  use m_hist, only: calculate_hist 
  use m_band, only: calculate_banddisp 
  use m_frmsf, only: wrt_frmsf 
  include "config.h" 
  !
  !read zvo files in dir-model 
  !
  call rd_dat_mkkpts 
  call rd_dat_hr 
  call rd_dat_geom 
  call rd_dat_bandkpts 
  call rd_dat_ef 
  !
  !read input from command line 
  !
  ncount=iargc() 
  allocate(real_arg(ncount)); real_arg=0.0d0 
  do i=1,ncount-4 
   call getarg(i,arg) 
   read(arg,*) real_arg(i)  
  enddo 
  !dos
  call getarg(ncount-3,arg)
  read(arg,*) dos 
  !bnd
  call getarg(ncount-2,arg)
  read(arg,*) bnd
  !frm 
  call getarg(ncount-1,arg)
  read(arg,*) frm 
  !his 
  call getarg(ncount-0,arg)
  read(arg,*) his 
  !
  !write(6,*) ncount
  !write(6,*) real_arg(1)
  !write(6,*) real_arg(2)
  !write(6,*) real_arg(3)
  !write(6,*) real_arg(4)
  !write(6,*) real_arg(5)
  !write(6,*) real_arg(6)
  !write(6,*) real_arg(7)
  !write(6,*) real_arg(8)
  !
  delt=real_arg(1)
  elnm=real_arg(2)  
  ecut=real_arg(3)  
  rcut=real_arg(4)  
  diff=real_arg(5)  
  !
  kgd(1)=nint(real_arg(6)); kgd(2)=nint(real_arg(7)); kgd(3)=nint(real_arg(8)) 
  !
  delt=delt/au !au <- eV
  electron_number=elnm 
  threshold_e=ecut/au !au <- eV
  threshold_r=rcut !AA 
  diff_transfers=diff/au !au <- eV  
  !
  write(6,*) 
  write(6,'(a50)')'##### TRANSFER ANALYSIS #####'
  write(6,*) 
  write(6,'(a50,x,f10.5)')'Greens function delt (eV):',delt*au
  write(6,'(a50,x,f10.5)')'Grid spacing delw (eV) [dos]:',2.0d0*delt*au
  write(6,'(a50,x,f10.5)')'Grid spacing delw (eV) [hist]:',10.0d0*delt*au
  write(6,'(a50,x,f10.5)')'Electron numbers in unit cell:',electron_number  
  write(6,'(a50,x,f10.5)')'Energy cutoff for transfer (eV):',threshold_e*au
  write(6,'(a50,x,f10.5)')'Distance cutoff for transfer (AA):',threshold_r 
  write(6,'(a50,x,f10.5)')'Match threshold for transfers (eV):',diff_transfers*au
  write(6,'(a50,x,f10.5)')'FermiEnergy in band calculation (eV):',FermiEnergy_bandcalc*au  
  if(kgd(1)==0.and.kgd(2)==0.and.kgd(3)==0)then 
      write(6,'(a50,x,3i10 )')'k grid:',nkb1,nkb2,nkb3  
  else
      write(6,'(a50,x,3i10 )')'k grid:',kgd 
  endif  
  write(6,*) 
  !
  !truncate H(R) on threshold 
  !
  call truncation(NWF,Na1,Na2,Na3,threshold_e,threshold_r,diff_transfers,a1(1),a2(1),a3(1),wcenter_lat(1,1),HR(1,1,-Na1,-Na2,-Na3))
  !
  if(dos)then 
    !
    !set kvec
    !
    if(kgd(1)==0.and.kgd(2)==0.and.kgd(3)==0)then 
        kgd(1)=nkb1; kgd(2)=nkb2; kgd(3)=nkb3 
    endif 
    ncalck=PRODUCT(kgd(1:3))
    allocate(kvec(3,ncalck)); kvec=0.0d0 
    call set_kgrid(kgd(1),kvec(1,1)) 
    write(6,*) 
    write(6,'(a7)')'kvec:' 
    do i=1,ncalck
     write(6,'(3f15.10)') kvec(:,i) 
    enddo 
    !
    !make EKS
    !
    allocate(EKS(NWF,ncalck)); EKS=0.0d0; allocate(VKS(NWF,NWF,ncalck)); VKS=0.0d0 
    flg_weight=0 !Flg whether calculate weighted transfers (0:not calc, 1:calc)
    call calculate_eigenstate(NWF,ncalck,Na1,Na2,Na3,nkb1,nkb2,nkb3,flg_weight,a1(1),a2(1),a3(1),kvec(1,1),HR(1,1,-Na1,-Na2,-Na3),EKS(1,1),VKS(1,1,1)) 
    !
    !make DOS 
    !
    delw=2.0d0*delt !Grid spacing delw [dos] 
    call calculate_dos(NWF,ncalck,kgd(1),kgd(2),kgd(3),electron_number,threshold_e,delt,delw,b1(1),b2(1),b3(1),kvec(1,1),EKS(1,1),VKS(1,1,1))
    deallocate(EKS,VKS) 
    deallocate(kvec) 
    !
    write(6,'(a50)')'##### FINISH TRANSFER ANALYSIS (DOS) #####'
  endif!dos 
  !
  if(bnd)then 
    !
    !set F(R) 
    !
    if(kgd(1)==0.and.kgd(2)==0.and.kgd(3)==0)then 
        kgd(1)=nkb1; kgd(2)=nkb2; kgd(3)=nkb3 
    endif 
    La1=kgd(1)/2; La2=kgd(2)/2; La3=kgd(3)/2 
    allocate(FR(NWF,NWF,-La1:La1,-La2:La2,-La3:La3)); FR(:,:,:,:,:)=0.0d0  
    call set_FR(NWF,nkb1,nkb2,nkb3,kgd(1),kgd(2),kgd(3),Na1,Na2,Na3,La1,La2,La3,HR(1,1,-Na1,-Na2,-Na3),FR(1,1,-La1,-La2,-La3)) 
    !
    !make EKS
    !
    allocate(EKS(NWF,NSK_BAND_DISP)); EKS=0.0d0; allocate(VKS(NWF,NWF,NSK_BAND_DISP)); VKS=0.0d0 
    flg_weight=0 !Flg whether calculate weighted transfers (0:not calc, 1:calc)
    call calculate_eigenstate(NWF,NSK_BAND_DISP,La1,La2,La3,kgd(1),kgd(2),kgd(3),flg_weight,a1(1),a2(1),a3(1),SK_BAND_DISP(1,1),FR(1,1,-La1,-La2,-La3),EKS(1,1),VKS(1,1,1)) 
    !
    !calc band disp
    !
    call calculate_banddisp(NWF,NSK_BAND_DISP,Ndiv,N_sym_points,threshold_e,b1(1),b2(1),b3(1),SK_BAND_DISP(1,1),EKS(1,1)) 
    deallocate(EKS,VKS) 
    !
    write(6,'(a50)')'##### FINISH TRANSFER ANALYSIS (BND) #####'
  endif!bnd 
  !
  if(frm)then
    !
    !wrt fermi surface 
    !
    call wrt_frmsf(NWF,kgd(1),Na1,Na2,Na3,nkb1,nkb2,nkb3,a1(1),a2(1),a3(1),b1(1),b2(1),b3(1),FermiEnergy_bandcalc,HR(1,1,-Na1,-Na2,-Na3)) 
    !
    write(6,'(a50)')'##### FINISH TRANSFER ANALYSIS (FRM) #####'
  endif!frm 
  !
  if(his)then
    write(6,*) 
    write(6,'(a50)')'##### TRANSFER ANALYSIS (trans.def) #####'
    write(6,*) 
    !
    !read trans.def 
    !
    call rd_transdef  
    !
    !calc histgram
    ! 
    delw=10.0d0*delt !Grid spacing delw [hist] 
    call calculate_hist(Ndim_TR,delw,TR(1,1)) 
    !
    write(6,'(a50)')'##### FINISH TRANSFER ANALYSIS (HIS) #####'
  endif!his
  !
  stop
end     
