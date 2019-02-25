subroutine calc_band_disp(Ndiv,N_sym_points,NTK,NSK_BAND_DISP,Na1,Na2,Na3,n_occ,SK_BAND_DISP,H_MAT_R,nkb1,nkb2,nkb3,a1,a2,a3,b1,b2,b3,&
  kdata,E_BAND_DISP)  
  implicit none 
  !
  integer::Ndiv,N_sym_points,NTK,NSK_BAND_DISP,Na1,Na2,Na3,n_occ
  integer::nkb1,nkb2,nkb3 
  real(8)::SK_BAND_DISP(3,NSK_BAND_DISP)
  real(8)::a1(3),a2(3),a3(3),b1(3),b2(3),b3(3) 
  complex(8)::H_MAT_R(n_occ,n_occ,-Na1:Na1,-Na2:Na2,-Na3:Na3)
  !
  real(8),allocatable::E_TMP(:)!E_TMP(n_occ)  
  real(8),allocatable::WEIGHT_R(:,:,:)!WEIGHT_R(-Na1:Na1,-Na2:Na2,-Na3:Na3)
  complex(8),allocatable::H_TMP_IN(:,:)!H_TMP_IN(n_occ,n_occ)
  complex(8),allocatable::H_TMP_OUT(:,:)!H_TMP_OUT(n_occ,n_occ)
  real(8),allocatable::DIST_K(:)!DIST_K(NTK)
  real(8),allocatable::dist(:)!dist(0:Nblk-1)
  !
  integer::ik,ib,jb,ix,ia1,ia2,ia3,ks,ke,ia1min,ia2min,ia3min 
  real(8)::DIST_B(3),DIST_KSPACE
  real(8)::PHASE 
  real(8)::SUM_REAL 
  complex(8)::PHASE_FACTOR     
  complex(8)::SUM_CMPX 
  !
  real(8),parameter::au=27.21151d0
  real(8),parameter::tpi=2.0d0*dacos(-1.0d0)
  complex(8),parameter::ci=(0.0D0,1.0D0) 
  !
  real(8),intent(out)::kdata(NSK_BAND_DISP)
  real(8),intent(out)::E_BAND_DISP(n_occ,NSK_BAND_DISP)
  !
  !1. SAMPLE k-POINTS FOR INTERPOLATED BAND DISPERSION 
  !
  write(6,*) 
  write(6,*)'==================================='
  write(6,*)'CALCULATED kpts FOR BAND DISPERSION'
  write(6,*)'==================================='
  write(6,*) 
  write(6,'(a20,i8)')'NSK_BAND_DISP=',NSK_BAND_DISP
  do ik=1,NSK_BAND_DISP 
   write(6,*) SK_BAND_DISP(:,ik)
  enddo!ik  
  !
  !2. WEIGHT_R BY YUSUKE NOMURA 
  !
  allocate(WEIGHT_R(-Na1:Na1,-Na2:Na2,-Na3:Na3)) 
  WEIGHT_R=1.0d0; SUM_REAL=0.0d0 
  do ia1=-Na1,Na1!-1         
   do ia2=-Na2,Na2!-1         
    do ia3=-Na3,Na3!-1         
     if(abs(ia1)==Na1.and.mod(nkb1,2)==0.and.Na1/=0) WEIGHT_R(ia1,ia2,ia3)=WEIGHT_R(ia1,ia2,ia3)*0.5d0 
     if(abs(ia2)==Na2.and.mod(nkb2,2)==0.and.Na2/=0) WEIGHT_R(ia1,ia2,ia3)=WEIGHT_R(ia1,ia2,ia3)*0.5d0 
     if(abs(ia3)==Na3.and.mod(nkb3,2)==0.and.Na3/=0) WEIGHT_R(ia1,ia2,ia3)=WEIGHT_R(ia1,ia2,ia3)*0.5d0 
     SUM_REAL=SUM_REAL+WEIGHT_R(ia1,ia2,ia3)
    enddo
   enddo
  enddo 
  write(6,'(a20,f15.8,i8)')'SUM_WEIGHT,NTK',SUM_REAL,NTK  
  if(abs(SUM_REAL-dble(NTK))>1.0d-6)then 
   stop 'SUM_WEIGHT/=NTK'
  endif 
  !
  !3. H_MAT(k') IN WANNIER BASIS AND DIAGONALIZE 
  !
  allocate(E_TMP(n_occ));E_TMP(:)=0.0d0            
  allocate(H_TMP_IN(n_occ,n_occ));H_TMP_IN(:,:)=0.0d0           
  allocate(H_TMP_OUT(n_occ,n_occ));H_TMP_OUT(:,:)=0.0d0  
  !allocate(E_BAND_DISP(n_occ,NSK_BAND_DISP))
  E_BAND_DISP=0.0d0 
  !
  do ik=1,NSK_BAND_DISP          
   H_TMP_IN(:,:)=0.0D0               
   do ib=1,n_occ      
    do jb=1,n_occ      
     SUM_CMPX=0.0D0                    
     do ia1=-Na1,Na1!-1         
      do ia2=-Na2,Na2!-1         
       do ia3=-Na3,Na3!-1         
        !
        !NEAREST R SEARCH BY YOSHIHIDE YOSHIMOTO 
        !
        call search_Rmin(ia1,ia2,ia3,nkb1,nkb2,nkb3,a1(1),a2(1),a3(1),ia1min,ia2min,ia3min)
        PHASE=tpi*(SK_BAND_DISP(1,ik)*DBLE(ia1min)+SK_BAND_DISP(2,ik)*DBLE(ia2min)+SK_BAND_DISP(3,ik)*DBLE(ia3min))                
        PHASE_FACTOR=EXP(ci*PHASE)*WEIGHT_R(ia1,ia2,ia3) 
        SUM_CMPX=SUM_CMPX+H_MAT_R(ib,jb,ia1,ia2,ia3)*PHASE_FACTOR 
       enddo!ia3            
      enddo!ia2                       
     enddo!ia1
     H_TMP_IN(ib,jb)=SUM_CMPX           
    enddo!jb
   enddo!ib
   !
   !4. diag H_TMP_IN
   !
   E_TMP(:)=0.0D0                
   call diagV(n_occ,H_TMP_IN(1,1),E_TMP(1)) 
   H_TMP_OUT(:,:)=H_TMP_IN(:,:)
   E_BAND_DISP(:,ik)=E_TMP(:)            
  enddo!ik           
  !
  !5. GENERATE KDATA (DISPERSION LINE) 
  !
  allocate(DIST_K(NSK_BAND_DISP));DIST_K(:)=0.0d0 
  allocate(dist(0:N_sym_points-1));dist(:)=0.0d0 
  !
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
  kdata=0.0d0 
  kdata(:)=DIST_K(:)/DIST_K(NSK_BAND_DISP)
  !
  deallocate(E_TMP,WEIGHT_R,H_TMP_IN,H_TMP_OUT,DIST_K,dist) 
  !
return
end
