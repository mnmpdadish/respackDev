!
!GENERATE KDATA (DISPERSION LINE) 
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
