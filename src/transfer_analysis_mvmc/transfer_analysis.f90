!
!(default parameters) 
!
!delt=0.01d0           !Greens function delt in eV
!dmna=0.001d0          !Ttrhdrn parameter dmna in eV
!dmnr=0.001d0          !Ttrhdrn parameter dmnr in eV
!delw=2.0d0*delt       !Grid width in eV [dos]
!delw=10.0d0*delt      ! Grid width in eV [hist]
!flwe=0                !Flg whether calculate weighted transfers (0:not calc, 1:calc)
!thtr=0.0d0            !Threshold for transfer integral
!elnm=0.0d0            !Total number of electrons in unitcell
!kgd='nkb1 nkb2 nkb3'  !k grid 
!
PROGRAM main 
  use m_rd_dat_zvo 
  use m_rd_transdef  
  use m_dos, only: calculate_dos 
  use m_eigenstate, only: calculate_eigenstate
  use m_hist, only: calculate_hist 
  use m_truncation, only: truncation 
  use m_band, only: calculate_banddisp 
  use m_frmsf, only: wrt_frmsf 
  include "config.h" 
  !
  !read zvo(mvmc files) 
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
  !
  delt=real_arg(1)
  flwe=real_arg(2)  
  thtr=real_arg(3)  
  elnm=real_arg(4)  
  !
  kgd(1)=nint(real_arg(5)); kgd(2)=nint(real_arg(6)); kgd(3)=nint(real_arg(7)) 
  if(kgd(1)==0) kgd(1)=nkb1
  if(kgd(2)==0) kgd(2)=nkb2
  if(kgd(3)==0) kgd(3)=nkb3 
  !
  flg_weight=nint(flwe)  
  if(flg_weight/=0) flg_weight=1
  delt=delt/au !au <- eV
  threshold_transfer=thtr/au !au <- eV
  electron_number=elnm 
  !
  write(6,*) 
  write(6,'(a50)')'##### TRANSFER ANALYSIS #####'
  write(6,*) 
  write(6,'(a50,x,f10.5)')'Greens function delt (eV)=',delt*au
  write(6,'(a50,x,f10.5)')'Grid spacing delw (eV) [dos]=',2.0d0*delt*au
  write(6,'(a50,x,f10.5)')'Grid spacing delw (eV) [hist]=',10.0d0*delt*au
  write(6,'(a50,x,i10  )')'Use weigted transfer (0:not)=',flg_weight 
  write(6,'(a50,x,f10.5)')'Threshold for transfer (eV)=',threshold_transfer*au
  write(6,'(a50,x,f10.5)')'Electron numbers in unit cell=',electron_number  
  write(6,'(a50,x,3i10 )')'k grid=',kgd 
  write(6,'(a50,x,f10.5 )')'FermiEnergy in band calculation (eV)=',FermiEnergy_bandcalc*au  
  write(6,*) 
  !
  if(dos)then 
    !
    !truncate HR
    !
    if(threshold_transfer/=0.0d0)then 
     call truncation(NWF,Na1,Na2,Na3,threshold_transfer,HR(1,1,-Na1,-Na2,-Na3))
    endif 
    !
    !make EKS
    !
    allocate(EKS(NWF,NTK)); EKS=0.0d0; allocate(VKS(NWF,NWF,NTK)); VKS=0.0d0 
    call calculate_eigenstate(NWF,NTK,Na1,Na2,Na3,nkb1,nkb2,nkb3,flg_weight,a1(1),a2(1),a3(1),SK0(1,1),HR(1,1,-Na1,-Na2,-Na3),EKS(1,1),VKS(1,1,1)) 
    !
    !make DOS 
    !
    delw=2.0d0*delt !Grid spacing delw [dos] 
    call calculate_dos(NWF,NTK,nkb1,nkb2,nkb3,electron_number,threshold_transfer,delt,delw,b1(1),b2(1),b3(1),SK0(1,1),EKS(1,1))
    deallocate(EKS,VKS) 
    !
    write(6,'(a50,f10.5)')'threshold_transfer [eV]=',threshold_transfer*au 
    write(6,'(a50)')'##### FINISH TRANSFER ANALYSIS (DOS) #####'
  endif!dos 
  !
  if(bnd)then 
    !
    !truncate HR
    !
    call truncation(NWF,Na1,Na2,Na3,threshold_transfer,HR(1,1,-Na1,-Na2,-Na3))
    !
    !make EKS
    !
    allocate(EKS(NWF,NSK_BAND_DISP)); EKS=0.0d0; allocate(VKS(NWF,NWF,NSK_BAND_DISP)); VKS=0.0d0 
    call calculate_eigenstate(NWF,NSK_BAND_DISP,Na1,Na2,Na3,nkb1,nkb2,nkb3,flg_weight,a1(1),a2(1),a3(1),SK_BAND_DISP(1,1),HR(1,1,-Na1,-Na2,-Na3),EKS(1,1),VKS(1,1,1)) 
    !
    !calc band disp
    !
    call calculate_banddisp(NWF,NSK_BAND_DISP,Ndiv,N_sym_points,threshold_transfer,b1(1),b2(1),b3(1),SK_BAND_DISP(1,1),EKS(1,1)) 
    deallocate(EKS,VKS) 
    !
    write(6,'(a50,f10.5)')'threshold_transfer [eV]=',threshold_transfer*au 
    write(6,'(a50)')'##### FINISH TRANSFER ANALYSIS (BND) #####'
  endif!bnd 
  !
  if(frm)then
    !
    write(6,'(a50)')'##### STILL NOT SURPORTED (FRM) #####'
    !
    !truncate HR
    !
    call truncation(NWF,Na1,Na2,Na3,threshold_transfer,HR(1,1,-Na1,-Na2,-Na3))
    !
    !wrt fermi surface 
    !
    call wrt_frmsf(NWF,kgd(1),Na1,Na2,Na3,nkb1,nkb2,nkb3,a1(1),a2(1),a3(1),b1(1),b2(1),b3(1),FermiEnergy_bandcalc,HR(1,1,-Na1,-Na2,-Na3)) 
    !
    write(6,'(a50,f10.5)')'threshold_transfer [eV]=',threshold_transfer*au 
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
