subroutine make_eig(NWF,NTK,Na1,Na2,Na3,nkb1,nkb2,nkb3,flg_weight,a1,a2,a3,SK0,KS_R,EKS,VKS) 
  implicit none 
  integer,intent(in)::NWF,NTK,Na1,Na2,Na3,nkb1,nkb2,nkb3 
  integer,intent(in)::flg_weight 
  real(8),intent(in)::a1(3),a2(3),a3(3)
  real(8),intent(in)::SK0(3,NTK) 
  complex(8),intent(in)::KS_R(NWF,NWF,-Na1:Na1,-Na2:Na2,-Na3:Na3)   
  real(8),intent(out)::EKS(NWF,NTK)           
  complex(8),intent(out)::VKS(NWF,NWF,NTK)   
  real(8),allocatable::WEIGHT_R(:,:,:)!WEIGHT_R(-Na1:Na1,-Na2:Na2,-Na3:Na3)
  complex(8),allocatable::pf(:,:,:,:)!pf(-Na1:Na1,-Na2:Na2,-Na3:Na3,NTK)   
  integer::ia1min,ia2min,ia3min 
  integer::ia1,ia2,ia3,ik 
  real(8)::PHASE 
  real(8)::SUM_REAL
  real(8),parameter::au=27.21151d0
  real(8),parameter::tpi=2.0d0*dacos(-1.0d0)
  complex(8),parameter::ci=(0.0D0,1.0D0) 
  !
  !1. WEIGHT_R BY Y.Nomura NOMURA 
  !
  allocate(WEIGHT_R(-Na1:Na1,-Na2:Na2,-Na3:Na3)); WEIGHT_R=1.0d0
  if(flg_weight.eq.1)then
   write(6,'(a20)')'WEIGHT CALCULATED' 
   SUM_REAL=0.0d0 
   do ia1=-Na1,Na1
    do ia2=-Na2,Na2
     do ia3=-Na3,Na3
      if(abs(ia1)==Na1.and.mod(NTK,2)==0.and.Na1/=0) WEIGHT_R(ia1,ia2,ia3)=WEIGHT_R(ia1,ia2,ia3)*0.5d0 
      if(abs(ia2)==Na2.and.mod(NTK,2)==0.and.Na2/=0) WEIGHT_R(ia1,ia2,ia3)=WEIGHT_R(ia1,ia2,ia3)*0.5d0 
      if(abs(ia3)==Na3.and.mod(NTK,2)==0.and.Na3/=0) WEIGHT_R(ia1,ia2,ia3)=WEIGHT_R(ia1,ia2,ia3)*0.5d0 
      SUM_REAL=SUM_REAL+WEIGHT_R(ia1,ia2,ia3)
     enddo!ia3
    enddo!ia2
   enddo!ia1
   write(6,'(a20,f15.8,i8)')'SUM_WEIGHT,NTK',SUM_REAL,NTK  
   if(abs(SUM_REAL-dble(NTK))>1.0d-6)then 
    stop 'SUM_WEIGHT/=NTK'
   endif 
  else
   write(6,'(a20)')'WEIGHT NOT CALCULATED' 
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
  !write(6,*)'# finish make pf'
  !
  !2. HKS IN WANNIER BASIS. AND DIAGONALIZE
  !
  EKS=0.0d0 
  VKS=0.0d0 
  call eigenvalue(NTK,NWF,Na1,Na2,Na3,KS_R(1,1,-Na1,-Na2,-Na3),pf(-Na1,-Na2,-Na3,1),EKS(1,1),VKS(1,1,1))  
  !
  deallocate(WEIGHT_R,pf) 
  !
  return 
end subroutine 
