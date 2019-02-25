subroutine calc_gwdos(NWF,NTK,nsgm,Na1,Na2,Na3,nkb1,nkb2,nkb3,idlt,dmna,dmnr,FermiEnergy,a1,a2,a3,b1,b2,b3,sgmw,SK0,&
  KS_R,XC_R,SX_R,SC_R,shift_value,gwdos,gw_sigma_dos) 
  !
  implicit none 
  integer::NWF,NTK,nsgm,Na1,Na2,Na3,nkb1,nkb2,nkb3 
  real(8)::idlt,dmna,dmnr
  real(8)::a1(3),a2(3),a3(3)
  real(8)::b1(3),b2(3),b3(3)
  real(8)::sgmw(nsgm)  
  real(8)::SK0(3,NTK) 
  complex(8)::KS_R(NWF,NWF,-Na1:Na1,-Na2:Na2,-Na3:Na3)   
  complex(8)::XC_R(NWF,NWF,-Na1:Na1,-Na2:Na2,-Na3:Na3)   
  complex(8)::SX_R(NWF,NWF,-Na1:Na1,-Na2:Na2,-Na3:Na3)   
  complex(4)::SC_R(NWF,NWF,-Na1:Na1,-Na2:Na2,-Na3:Na3,nsgm)   
  !
  real(8),allocatable::WEIGHT_R(:,:,:)!WEIGHT_R(-Na1:Na1,-Na2:Na2,-Na3:Na3)
  real(8),allocatable::EKS(:,:)!EKS(NWF,NTK)           
  complex(8),allocatable::VKS(:,:,:)!VKS(NWF,NWF,NTK)   
  complex(8),allocatable::pf(:,:,:,:)!pf(-Na1:Na1,-Na2:Na2,-Na3:Na3,NTK)   
  complex(8),allocatable::EMK(:,:,:)!EMK(NWF,NTK,nsgm)           
  complex(4),allocatable::GW_R(:,:,:,:,:,:)!GW_R(NWF,NWF,-Na1:Na1,-Na2:Na2,-Na3:Na3,nsgm)   
  ! 
  integer::ia1min,ia2min,ia3min 
  integer::ia1,ia2,ia3,ik,ib,jb,ie  
  real(8)::PHASE,FermiEnergy 
  real(8)::SUM_REAL
  complex(8)::SUM_CMPX 
  !
  real(8),parameter::au=27.21151d0
  real(8),parameter::tpi=2.0d0*dacos(-1.0d0)
  complex(8),parameter::ci=(0.0D0,1.0D0) 
  !
  real(8),intent(out)::shift_value 
  real(8),intent(out)::gwdos(nsgm) 
  complex(8),intent(out)::gw_sigma_dos(nsgm) 
  !
  !1. make GR_R 
  !
  allocate(GW_R(NWF,NWF,-Na1:Na1,-Na2:Na2,-Na3:Na3,nsgm));GW_R=0.0d0   
  do ie=1,nsgm
   do ia1=-Na1,Na1
    do ia2=-Na2,Na2
     do ia3=-Na3,Na3
      do ib=1,NWF          
       do jb=1,NWF          
        GW_R(ib,jb,ia1,ia2,ia3,ie)&
       =KS_R(ib,jb,ia1,ia2,ia3)&
       -XC_R(ib,jb,ia1,ia2,ia3)&
       -SX_R(ib,jb,ia1,ia2,ia3)&
       +SC_R(ib,jb,ia1,ia2,ia3,ie) 
       enddo 
      enddo 
     enddo 
    enddo 
   enddo 
  enddo 
  write(6,*)'# finish make GW_R'
  !
  !2. WEIGHT_R BY Y.Nomura NOMURA 
  !
  allocate(WEIGHT_R(-Na1:Na1,-Na2:Na2,-Na3:Na3)); WEIGHT_R=1.0d0
  SUM_REAL=0.0d0 
  do ia1=-Na1,Na1
   do ia2=-Na2,Na2
    do ia3=-Na3,Na3
     if(abs(ia1)==Na1.and.mod(nkb1,2)==0.and.Na1/=0) WEIGHT_R(ia1,ia2,ia3)=WEIGHT_R(ia1,ia2,ia3)*0.5d0 
     if(abs(ia2)==Na2.and.mod(nkb2,2)==0.and.Na2/=0) WEIGHT_R(ia1,ia2,ia3)=WEIGHT_R(ia1,ia2,ia3)*0.5d0 
     if(abs(ia3)==Na3.and.mod(nkb3,2)==0.and.Na3/=0) WEIGHT_R(ia1,ia2,ia3)=WEIGHT_R(ia1,ia2,ia3)*0.5d0 
     SUM_REAL=SUM_REAL+WEIGHT_R(ia1,ia2,ia3)
    enddo!ia3
   enddo!ia2
  enddo!ia1
  write(6,'(a20,f15.8,i8)')'SUM_WEIGHT,NTK',SUM_REAL,NTK  
  if(abs(SUM_REAL-dble(NTK))>1.0d-6)then 
   stop 'SUM_WEIGHT/=NTK'
  endif 
  !
  allocate(pf(-Na1:Na1,-Na2:Na2,-Na3:Na3,NTK)); pf=0.0d0 
  do ik=1,NTK 
   do ia3=-Na3,Na3 
    do ia2=-Na2,Na2 
     do ia1=-Na1,Na1 
      !
      !NEAREST R SEARCH BY Y.Yoshimoto 
      !
      call search_Rmin(ia1,ia2,ia3,nkb1,nkb2,nkb3,a1(1),a2(1),a3(1),ia1min,ia2min,ia3min)
      PHASE=tpi*(SK0(1,ik)*DBLE(ia1min)+SK0(2,ik)*DBLE(ia2min)+SK0(3,ik)*DBLE(ia3min))           
      pf(ia1,ia2,ia3,ik)=EXP(ci*PHASE)*WEIGHT_R(ia1,ia2,ia3)          
     enddo!ia1 
    enddo!ia2 
   enddo!ia3 
  enddo!ik  
  write(6,*)'# finish make pf'
  !
  !3. HKS IN WANNIER BASIS. AND DIAGONALIZE
  !
  allocate(EKS(NWF,NTK)); EKS=0.0d0 
  allocate(VKS(NWF,NWF,NTK)); VKS=0.0d0 
  call make_eks(NTK,NWF,Na1,Na2,Na3,KS_R(1,1,-Na1,-Na2,-Na3),pf(-Na1,-Na2,-Na3,1),EKS(1,1),VKS(1,1,1))  
  !
  !4. HGW IN WANNIER BASIS. AND DIAGONALIZE
  !
  allocate(EMK(NWF,NTK,nsgm)); EMK=0.0d0 
  call make_emk(NTK,NWF,nsgm,Na1,Na2,Na3,GW_R(1,1,-Na1,-Na2,-Na3,1),pf(-Na1,-Na2,-Na3,1),EMK(1,1,1))  
  !
  !5. DETERMINE SHIFT 
  !
  call det_shift(NTK,NWF,nsgm,FermiEnergy,sgmw(1),EKS(1,1),EMK(1,1,1),shift_value) 
  !
  !6. GW-DOS CALC
  !
  gwdos=0.0d0  
  call calc_dos_GW(nkb1,nkb2,nkb3,NTK,NWF,nsgm,idlt,dmna,dmnr,shift_value,b1(1),b2(1),b3(1),sgmw(1),SK0(1,1),EMK(1,1,1),gwdos(1)) 
  ! 
  !7. GW-SIGMA-DOS CALC
  !
  gw_sigma_dos=0.0d0  
  do ie=1,nsgm 
   SUM_CMPX=0.0d0 
   do ik=1,NTK 
    do ib=1,NWF
     SUM_CMPX=SUM_CMPX+(EMK(ib,ik,ie)-cmplx(EKS(ib,ik)))   
    enddo!ib
   enddo!ik
   gw_sigma_dos(ie)=SUM_CMPX/dble(NTK) 
  enddo!ie 
  !
  deallocate(WEIGHT_R,EKS,VKS,pf,EMK,GW_R) 
  !
return 
end
