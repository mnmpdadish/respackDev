subroutine wrt_mvmc(NTK,Na1,Na2,Na3,n_occ,HR,WR) 
  implicit none 
  integer::NTK,Na1,Na2,Na3,n_occ 
  integer::L1,L2,L3,Nsite 
  integer::ia1,ia2,ia3,ib,jb,i 
  integer::unit_vec(NTK) 
  complex(8)::HR(n_occ,n_occ,-Na1:Na1,-Na2:Na2,-Na3:Na3)
  real(8)::WR(-Na1:Na1,-Na2:Na2,-Na3:Na3)
  real(8),parameter::au=27.21151d0
  !
  !OPEN(300,W,FILE='zvo_hr.dat') 
  !
  OPEN(300,FILE='./dir-mvmc/zvo_hr.dat') 
  write(300,'(a)')'wannier90 format for mvmcdry'
  write(300,'(i10)') n_occ
  write(300,'(i10)') NTK 
  unit_vec=1
  write(300,'(15i5)')(unit_vec(i),i=1,NTK) 
  do ia1=-Na1,Na1
   do ia2=-Na2,Na2
    do ia3=-Na3,Na3 
     do ib=1,n_occ
      do jb=1,n_occ
       write(300,'(i5,i5,i5,i5,i5,2f15.10)') ia1,ia2,ia3,ib,jb,HR(ib,jb,ia1,ia2,ia3)*WR(ia1,ia2,ia3)*au   
      enddo!jb
     enddo!ib
    enddo!ia3
   enddo!ia2
  enddo!ia1 
return
end 
!
subroutine wrt_mvmc_wcenter(n_occ,a1,a2,a3,wcenter)
  implicit none 
  integer::n_occ 
  real(8)::a1(3),a2(3),a3(3),wcenter(3,n_occ) 
  real(8)::ainv(3,3)
  integer::ib,i,j
  real(8)::wcenter_lat(3) 
  real(8)::SUM_REAL
  real(8),parameter::bohr=0.529177249d0 
  !
  ainv(:,1)=a1(:)
  ainv(:,2)=a2(:)
  ainv(:,3)=a3(:)
  !
  call invmat(3,ainv(1,1)) 
  !
  !
  !OPEN(303,W,FILE='zvo_geom.dat') 
  !
  OPEN(303,FILE='./dir-mvmc/zvo_geom.dat') 
  write(303,'(3f15.10)')(a1(i)*bohr,i=1,3) 
  write(303,'(3f15.10)')(a2(i)*bohr,i=1,3) 
  write(303,'(3f15.10)')(a3(i)*bohr,i=1,3) 
  write(303,'(i10)') n_occ
  do ib=1,n_occ 
   wcenter_lat=0.0d0 
   do i=1,3
    SUM_REAL=0.0d0 
    do j=1,3
     SUM_REAL=SUM_REAL+ainv(i,j)*wcenter(j,ib) 
    enddo!j 
    wcenter_lat(i)=SUM_REAL
   enddo!i
   write(303,*)(wcenter_lat(i),i=1,3) 
  enddo 
  !
return
end
