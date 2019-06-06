module m_rd_dat_zvo 
  implicit none
  private 
  public::rd_dat_hr
  public::rd_dat_geom 
  public::rd_dat_bandkpts
  public::rd_dat_mkkpts 
  public::rd_dat_ef 
  !h_mat_r(300) 
  integer,public::NWF 
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
  integer,public::NTK 
  integer,public::nkb1,nkb2,nkb3 
  integer,public::Na1,Na2,Na3 
  real(8),public,allocatable::SK0(:,:)!SK0(3,NTK) 
  !ef(306) 
  real(8),public::FermiEnergy_bandcalc 
  contains
  !--
  subroutine rd_dat_hr 
    implicit none 
    integer::dum1,dum2,dum3,dumi,dumj
    integer::N_element 
    character(len=80)::dum_ch
    integer,allocatable::unit_vec(:)!unit_vec(N_element) 
    integer::i,ia1,ia2,ia3,ib,jb 
    real(8),parameter::au=27.21151d0
    !
    !OPEN(300,R,FILE='./dir-model/zvo_hr.dat') 
    !
    OPEN(300,FILE='./dir-model/zvo_hr.dat') 
    rewind(300) 
    read(300,'(a)') dum_ch 
    read(300,'(i10)') NWF 
    read(300,'(i10)') N_element 
    allocate(unit_vec(N_element)); unit_vec=1
    read(300,'(15i5)')(unit_vec(i),i=1,N_element) 
    deallocate(unit_vec) 
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
  ! 
  subroutine rd_dat_geom 
    implicit none 
    real(8),parameter::tpi=2.0d0*dacos(-1.0d0)
    integer::ib,i 
    !
    !OPEN(303,R,FILE='./dir-model/zvo_geom.dat') 
    !
    OPEN(303,FILE='./dir-model/zvo_geom.dat') 
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
  !
  subroutine rd_dat_bandkpts 
    implicit none 
    integer::ik,i 
    !
    !OPEN(304,R,FILE='zvo_bandkpts.dat') 
    !
    OPEN(304,FILE='./dir-model/zvo_bandkpts.dat') 
    read(304,'(3i10)') NSK_BAND_DISP,Ndiv,N_sym_points 
    allocate(SK_BAND_DISP(3,NSK_BAND_DISP));SK_BAND_DISP=0.0d0 
    do ik=1,NSK_BAND_DISP
     read(304,*)(SK_BAND_DISP(i,ik),i=1,3) 
    enddo 
  end subroutine
  !
  subroutine rd_dat_mkkpts 
    implicit none 
    integer::ik,i 
    !
    !OPEN(305,R,FILE='zvo_mkkpts.dat') 
    !
    OPEN(305,FILE='./dir-model/zvo_mkkpts.dat') 
    read(305,'(i10)') NTK 
    allocate(SK0(3,NTK));SK0=0.0d0 
    do ik=1,NTK 
     read(305,*)(SK0(i,ik),i=1,3) 
    enddo 
    call est_nkbi(NTK,SK0,nkb1,nkb2,nkb3)  
    Na1=nkb1/2; Na2=nkb2/2; Na3=nkb3/2
  end subroutine
  !
  subroutine rd_dat_ef 
    implicit none 
    !
    !OPEN(306,R,FILE='zvo_ef.dat') 
    !
    OPEN(306,FILE='./dir-model/zvo_ef.dat') 
    read(306,'(f15.10)') FermiEnergy_bandcalc 
  end subroutine rd_dat_ef 
  !
  subroutine OUTER_PRODUCT(vec_x,vec_y,vec_z)
    implicit none 
    real(8)::vec_x(3),vec_y(3),vec_z(3) 
    !
    vec_z(1)=vec_x(2)*vec_y(3)-vec_x(3)*vec_y(2)
    vec_z(2)=vec_x(3)*vec_y(1)-vec_x(1)*vec_y(3) 
    vec_z(3)=vec_x(1)*vec_y(2)-vec_x(2)*vec_y(1)
    !
    return
  end subroutine 
  !
  subroutine est_nkbi(N,SK,nkb1,nkb2,nkb3)  
    implicit none 
    integer,intent(in)::N
    real(8),intent(in)::SK(3,N) 
    integer,intent(out)::nkb1,nkb2,nkb3
    integer::NTK  
    integer::i 
    real(8)::x 
    real(8),parameter::dlt_BZ=1.0d-6 
    x=1.0d0 
    do i=1,N
     if(abs(SK(1,i))<dlt_BZ) cycle 
     if(abs(SK(1,i))<x) then 
      x=abs(SK(1,i))  
     endif 
    enddo    
    nkb1=nint(1.0d0/x)  
    x=1.0d0 
    do i=1,N
     if(abs(SK(2,i))<dlt_BZ) cycle 
     if(abs(SK(2,i))<x) then 
      x=abs(SK(2,i))  
     endif 
    enddo    
    nkb2=nint(1.0d0/x)  
    x=1.0d0 
    do i=1,N
     if(abs(SK(3,i))<dlt_BZ) cycle 
     if(abs(SK(3,i))<x) then 
      x=abs(SK(3,i))  
     endif 
    enddo    
    nkb3=nint(1.0d0/x)  
    NTK=nkb1*nkb2*nkb3 
    write(6,'(a50,4i10)')'nkb1, nkb2, nkb3, NTK:',nkb1,nkb2,nkb3,NTK  
    return 
  end subroutine 
  !
end module 
