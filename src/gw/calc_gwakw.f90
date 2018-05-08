subroutine calc_gwakw(NWF,NTK,nsgm,Na1,Na2,Na3,nkb1,nkb2,nkb3,Ndiv,N_sym_points,NSK_BAND_DISP,idlt,shift_ef,&
  SK_BAND_DISP,a1,a2,a3,b1,b2,b3,sgmw,KS_R,XC_R,SX_R,SC_R)
  !
  implicit none 
  integer::NWF,NTK,nsgm,Na1,Na2,Na3,nkb1,nkb2,nkb3,Ndiv,N_sym_points,NSK_BAND_DISP
  real(8)::idlt,shift_ef 
  real(8)::SK_BAND_DISP(3,NSK_BAND_DISP)
  real(8)::a1(3),a2(3),a3(3)
  real(8)::b1(3),b2(3),b3(3)
  real(8)::sgmw(nsgm)  
  complex(8)::KS_R(NWF,NWF,-Na1:Na1,-Na2:Na2,-Na3:Na3)   
  complex(8)::XC_R(NWF,NWF,-Na1:Na1,-Na2:Na2,-Na3:Na3)   
  complex(8)::SX_R(NWF,NWF,-Na1:Na1,-Na2:Na2,-Na3:Na3)   
  complex(4)::SC_R(NWF,NWF,-Na1:Na1,-Na2:Na2,-Na3:Na3,nsgm)   
  !
  real(8),allocatable::AKW(:,:)!AKW(NSK_BAND_DISP,nsgm)           
  real(8),allocatable::WEIGHT_R(:,:,:)!WEIGHT_R(-Na1:Na1,-Na2:Na2,-Na3:Na3)
  complex(8),allocatable::pf(:,:,:,:)!pf(-Na1:Na1,-Na2:Na2,-Na3:Na3,NSK_BAND_DISP)   
  complex(8),allocatable::EMK(:,:,:)!EMK(NWF,NSK_BAND_DISP,nsgm)           
  complex(4),allocatable::GW_R(:,:,:,:,:,:)!GW_R(NWF,NWF,-Na1:Na1,-Na2:Na2,-Na3:Na3,nsgm)   
  !
  real(8),allocatable::DIST_K(:)!DIST_K(NSK_BAND_DISP)
  real(8),allocatable::dist(:)!dist(0:Nblk-1)
  real(8)::DIST_B(3),DIST_KSPACE
  !
  integer::ik,ib,jb,ix,ia1,ia2,ia3,ie,ks,ke,ia1min,ia2min,ia3min 
  real(8)::PHASE           
  real(8)::omega,delta 
  real(8)::SUM_REAL 
  complex(8)::w,d,en 
  !
  real(8),parameter::au=27.21151d0
  real(8),parameter::tpi=2.0d0*dacos(-1.0d0)
  complex(8),parameter::ci=(0.0D0,1.0D0) 
  ! 
  !GW_R 
  !
  allocate(GW_R(NWF,NWF,-Na1:Na1,-Na2:Na2,-Na3:Na3,nsgm));GW_R=0.0d0   
  do ia1=-Na1,Na1 
   do ia2=-Na2,Na2 
    do ia3=-Na3,Na3 
     do ib=1,NWF          
      do jb=1,NWF          
       do ie=1,nsgm
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
  !
  write(6,*)'# finish make H(R)'
  !
  allocate(WEIGHT_R(-Na1:Na1,-Na2:Na2,-Na3:Na3)); WEIGHT_R=1.0d0
  SUM_REAL=0.0d0 
  do ia1=-Na1,Na1
   do ia2=-Na2,Na2
    do ia3=-Na3,Na3
     if(abs(ia1)==Na1.and.mod(NTK,2)==0.and.Na1/=0) WEIGHT_R(ia1,ia2,ia3)=WEIGHT_R(ia1,ia2,ia3)*0.5d0 
     if(abs(ia2)==Na2.and.mod(NTK,2)==0.and.Na2/=0) WEIGHT_R(ia1,ia2,ia3)=WEIGHT_R(ia1,ia2,ia3)*0.5d0 
     if(abs(ia3)==Na3.and.mod(NTK,2)==0.and.Na3/=0) WEIGHT_R(ia1,ia2,ia3)=WEIGHT_R(ia1,ia2,ia3)*0.5d0 
     SUM_REAL=SUM_REAL+WEIGHT_R(ia1,ia2,ia3)
    enddo
   enddo
  enddo 
  write(6,'(a20,f15.8,i8)')'SUM_WEIGHT,NTK',SUM_REAL,NTK  
  if(abs(SUM_REAL-dble(NTK))>1.0d-6)then 
   stop 'SUM_WEIGHT/=NTK'
  endif 
  !
  allocate(pf(-Na1:Na1,-Na2:Na2,-Na3:Na3,NSK_BAND_DISP)); pf=0.0d0 
  do ik=1,NSK_BAND_DISP          
   do ia3=-Na3,Na3 
    do ia2=-Na2,Na2 
     do ia1=-Na1,Na1 
      !
      !NEAREST R SEARCH BY Y.Yoshimoto 
      !
      call search_Rmin(ia1,ia2,ia3,nkb1,nkb2,nkb3,a1(1),a2(1),a3(1),ia1min,ia2min,ia3min)
      PHASE=tpi*(SK_BAND_DISP(1,ik)*DBLE(ia1min)+SK_BAND_DISP(2,ik)*DBLE(ia2min)+SK_BAND_DISP(3,ik)*DBLE(ia3min))           
      pf(ia1,ia2,ia3,ik)=EXP(ci*PHASE)*WEIGHT_R(ia1,ia2,ia3)          
     enddo 
    enddo 
   enddo 
  enddo 
  write(6,*)'# finish make pf'
  !
  !H(k) IN WANNIER BASIS. AND DIAGONALIZE
  !
  allocate(EMK(NWF,NSK_BAND_DISP,nsgm)); EMK=0.0d0 
  call make_emk(NSK_BAND_DISP,NWF,nsgm,Na1,Na2,Na3,GW_R(1,1,-Na1,-Na2,-Na3,1),pf(-Na1,-Na2,-Na3,1),EMK(1,1,1))  
  !
  !CALC SPECTRAL FUNCTION 
  !
  allocate(DIST_K(NSK_BAND_DISP));DIST_K(:)=0.0d0 
  allocate(dist(0:N_sym_points-1));dist(:)=0.0d0 
  !--
  dist(0)=0.0d0 
  do ix=1,N_sym_points-1 
   ks=(ix-1)*Ndiv+1
   ke=(ix)*Ndiv+1
   do ik=ks,ke 
    DIST_B(:)=(SK_BAND_DISP(1,ik)-SK_BAND_DISP(1,ks))*b1(:)&
             +(SK_BAND_DISP(2,ik)-SK_BAND_DISP(2,ks))*b2(:)&
             +(SK_BAND_DISP(3,ik)-SK_BAND_DISP(3,ks))*b3(:)
    DIST_KSPACE=DSQRT(DIST_B(1)**2+DIST_B(2)**2+DIST_B(3)**2)
    DIST_K(ik)=DIST_KSPACE+dist(ix-1) 
   enddo!ik 
   dist(ix)=dist(ix-1)+DIST_KSPACE
  enddo!ix 
  !
  allocate(AKW(NSK_BAND_DISP,nsgm)); AKW=0.0d0
  do ik=1,NSK_BAND_DISP 
   do jb=1,NWF 
    do ie=1,nsgm 
     w=cmplx(sgmw(ie)+shift_ef,-idlt) 
     en=EMK(jb,ik,ie)
     d=1.0d0/(w-en) 
     AKW(ik,ie)=AKW(ik,ie)+abs(imag(d)) 
    enddo 
   enddo 
  enddo 
  !
  rewind(161) 
  do ie=1,nsgm 
   omega=sgmw(ie) 
   do ik=1,NSK_BAND_DISP 
    delta=DIST_K(ik)/DIST_K(NSK_BAND_DISP)
    if(AKW(ik,ie)>=500.0d0) AKW(ik,ie)=500.0d0 
    write(161,*) delta,omega*au,AKW(ik,ie) 
   enddo 
   write(161,*) 
  enddo 
  !
  deallocate(AKW,WEIGHT_R,pf,EMK,GW_R,DIST_K,dist) 
  !
return
end
