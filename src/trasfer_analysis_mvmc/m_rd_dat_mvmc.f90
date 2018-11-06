module m_rd_dat_mvmc
  implicit none
  private 
  public::rd_dat_hr
  public::rd_dat_geom 
  public::rd_dat_bandkpts
  public::rd_dat_mkkpts 
  !h_mat_r(300) 
  integer,public::Na1,Na2,Na3,NWF,NTK 
  complex(8),public,allocatable::HR(:,:,:,:,:)!HR(NWF,NWF,-Na1:Na1,-Na2:Na2,-Na3:Na3) 
  !geom(303)
  integer,public::N_wannier_center 
  real(8),public::a1(3),a2(3),a3(3) 
  real(8),public::b1(3),b2(3),b3(3)
  real(8),public::VOLUME 
  real(8),public,allocatable::wcenter_lat(:,:)!wcenter_lat(3,N_wannier_center) 
  !bandkpts(304)
  integer,public::Ndiv,N_sym_points,NSK_BAND_DISP
  real(8),public,allocatable::SK_BAND_DISP(:,:)!SK_BAND_DISP(3,NSK_BAND_DISP) 
  !mkkpts(305)  
  integer,public::N_mkmesh 
  real(8),public,allocatable::SK0(:,:)!SK0(3,N_mkmesh) 
  contains
!--
  subroutine rd_dat_hr 
  implicit none 
  integer::dum1,dum2,dum3,dumi,dumj 
  character(len=80)::dum_ch
  integer,allocatable::unit_vec(:)!unit_vec(NTK) 
  integer::i,ia1,ia2,ia3,ib,jb 
  real(8),parameter::au=27.21151d0
  !
  !OPEN(300,R,FILE='./dir-mvmc/zvo_hr.dat') 
  !
  OPEN(300,FILE='./dir-mvmc/zvo_hr.dat') 
  rewind(300) 
  read(300,'(a)') dum_ch 
  read(300,'(i10)') NWF 
  read(300,'(i10)') NTK,Na1,Na2,Na3  
  allocate(unit_vec(NTK)); unit_vec=1
  read(300,'(15i5)')(unit_vec(i),i=1,NTK) 
  allocate(HR(NWF,NWF,-Na1:Na1,-Na2:Na2,-Na3:Na3));HR(:,:,:,:,:)=0.0d0  
  do ia1=-Na1,Na1
   do ia2=-Na2,Na2
    do ia3=-Na3,Na3 
     do ib=1,NWF!n_occ
      do jb=1,NWF!n_occ
       read(300,'(i5,i5,i5,i5,i5,2f15.10)') dum1,dum2,dum3,dumi,dumj,HR(ib,jb,ia1,ia2,ia3) 
      enddo!jb
     enddo!ib
    enddo!ia3
   enddo!ia2
  enddo!ia1 
  HR=HR/au !au<-eV
  end subroutine
  
  subroutine rd_dat_geom 
  implicit none 
  real(8),parameter::tpi=2.0d0*dacos(-1.0d0)
  integer::ib,i 
  !
  !OPEN(303,R,FILE='./dir-mvmc/zvo_geom.dat') 
  !
  OPEN(303,FILE='./dir-mvmc/zvo_geom.dat') 
  read(303,'(3f15.10)')(a1(i),i=1,3) 
  read(303,'(3f15.10)')(a2(i),i=1,3) 
  read(303,'(3f15.10)')(a3(i),i=1,3) 
  read(303,'(i10)') N_wannier_center 
  allocate(wcenter_lat(3,N_wannier_center)) 
  do ib=1,N_wannier_center 
   read(303,*)(wcenter_lat(i,ib),i=1,3) 
  enddo 
  !
  call OUTER_PRODUCT(a2(1),a3(1),b1(1))
  VOLUME=a1(1)*b1(1)+a1(2)*b1(2)+a1(3)*b1(3)
  b1(:)=b1(:)*tpi/VOLUME 
  call OUTER_PRODUCT(a3(1),a1(1),b2(1))
  b2(:)=b2(:)*tpi/VOLUME 
  call OUTER_PRODUCT(a1(1),a2(1),b3(1))
  b3(:)=b3(:)*tpi/VOLUME 
  end subroutine

  subroutine rd_dat_bandkpts 
  implicit none 
  integer::ik,i 
  !
  !OPEN(304,R,FILE='zvo_bandkpts.dat') 
  !
  OPEN(304,FILE='./dir-mvmc/zvo_bandkpts.dat') 
  read(304,'(3i10)') NSK_BAND_DISP,Ndiv,N_sym_points 
  allocate(SK_BAND_DISP(3,NSK_BAND_DISP));SK_BAND_DISP=0.0d0 
  do ik=1,NSK_BAND_DISP
   read(304,*)(SK_BAND_DISP(i,ik),i=1,3) 
  enddo 
  end subroutine

  subroutine rd_dat_mkkpts 
  implicit none 
  integer::ik,i 
  !
  !OPEN(305,R,FILE='zvo_mkkpts.dat') 
  !
  OPEN(305,FILE='./dir-mvmc/zvo_mkkpts.dat') 
  read(305,'(i10)') N_mkmesh 
  allocate(SK0(3,N_mkmesh));SK0=0.0d0 
  do ik=1,N_mkmesh 
   read(305,*)(SK0(i,ik),i=1,3) 
  enddo 
  end subroutine
!--
end module 
