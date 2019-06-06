module m_band 
  implicit none
contains
  !
  subroutine calculate_banddisp(NWF,NSK_BAND_DISP,Ndiv,N_sym_points,threshold_transfer,b1,b2,b3,SK_BAND_DISP,EKS) 
    implicit none 
    integer,intent(in)::NWF,NSK_BAND_DISP,Ndiv,N_sym_points
    real(8),intent(in)::threshold_transfer 
    real(8),intent(in)::b1(3),b2(3),b3(3)
    real(8),intent(in)::SK_BAND_DISP(3,NSK_BAND_DISP) 
    real(8),intent(in)::EKS(NWF,NSK_BAND_DISP) 
    !
    real(8),allocatable::kdata(:)!kdata(NSK_BAND_DISP) 
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
    return
  end subroutine      
  !
  subroutine make_kdata(Ndiv,N_sym_points,NSK_BAND_DISP,SK_BAND_DISP,b1,b2,b3,kdata) 
    implicit none 
    integer,intent(in)::Ndiv,N_sym_points,NSK_BAND_DISP 
    real(8),intent(in)::SK_BAND_DISP(3,NSK_BAND_DISP)
    real(8),intent(in)::b1(3),b2(3),b3(3) 
    real(8),intent(out)::kdata(NSK_BAND_DISP)
    integer::ik,ix,ks,ke 
    real(8)::DIST_B(3),DIST_KSPACE
    real(8),allocatable::DIST_K(:)!DIST_K(NTK)
    real(8),allocatable::dist(:)!dist(0:Nblk-1)
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
    deallocate(DIST_K,dist) 
    return
  end subroutine 
  !
  subroutine wrt_iband(threshold_transfer,NWF,NSK_BAND_DISP,kdata,EKS) 
    implicit none
    integer,intent(in)::NWF,NSK_BAND_DISP
    real(8),intent(in)::kdata(NSK_BAND_DISP)
    real(8),intent(in)::EKS(NWF,NSK_BAND_DISP) 
    real(8),intent(in)::threshold_transfer 
    integer::ib,ik 
    LOGICAL::REVERSE 
    real(8),parameter::au=27.21151d0 
    !
    !OPEN(114,W,FILE='./dat.iband') 
    !
    OPEN(114,FILE='./dat.iband') 
    write(114,'(a)')'#Wannier interpolaed band'
    write(114,'(a,x,f10.5)')'#Energy cutoff for transfer (eV)=',threshold_transfer*au
    write(114,'(a)')'#1:k, 2:Energy [eV]' 
    REVERSE=.TRUE.        
    do ib=1,NWF 
     if(REVERSE)then 
      do ik=1,NSK_BAND_DISP                     
       write(114,'(2f20.10)') kdata(ik),EKS(ib,ik)*au
      enddo!ik        
      REVERSE=.FALSE.        
     else         
      do ik=NSK_BAND_DISP,1,-1          
       write(114,'(2f20.10)') kdata(ik),EKS(ib,ik)*au
      enddo!ik        
      REVERSE=.TRUE.        
     endif!REVERSE                   
    enddo!ib 
    close(114) 
    !
    return 
  end subroutine   
  !
end module m_band 
