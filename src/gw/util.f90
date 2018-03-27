subroutine invmat(nm,mat)
  implicit none 
  integer , intent(in) :: nm
  real(8) , intent(inout) :: mat(nm,nm)
  integer :: ipiv(nm)
  integer :: Lwork 
  real(8) , allocatable :: work(:)
  integer :: info 

  Lwork = 10*nm
  allocate (work(Lwork))
  info = 0
  call dgetrf(nm,nm,mat,nm,ipiv,info)
  call dgetri(nm,mat,nm,ipiv,work,Lwork,info)

  if(info /= 0) then
    write(6,*) 'info (subrouitine inv):' , info
    stop
  end if 
  deallocate(work)
  return 
end subroutine
!
subroutine invmat_complex(nm,mat)
  implicit none 
  integer,intent(in)::nm
  complex(8),intent(inout)::mat(nm,nm)
  integer::ipiv(nm)
  integer::Lwork 
  complex(8),allocatable::work(:)
  integer::info 
  !
  Lwork = 10*nm
  allocate (work(Lwork))
  info = 0
  call zgetrf(nm,nm,mat,nm,ipiv,info)
  call zgetri(nm,mat,nm,ipiv,work,Lwork,info)
  !
  if(info /= 0) then
    write(6,*) 'info (subrouitine inv):' , info
    stop
  end if 
  deallocate(work)
  return 
end subroutine
!
subroutine diagV(nm,mat,eig)
  implicit none 
  integer,intent(in)::nm
  complex(8),intent(inout)::mat(nm,nm)
  real(8),intent(out)::eig(nm)
  integer::LWORK,LRWORK,LIWORK  
  integer,allocatable::iwork_zheevd(:)
  real(8),allocatable::rwork_zheevd(:)
  complex(8),allocatable::work_zheevd(:)
  integer::ind
  real(8)::eps 
  !
  LWORK= 2*nm+nm**2
  LRWORK=1+12*nm+3*nm**2
  LIWORK=3+10*nm 
  allocate(work_zheevd(LWORK));work_zheevd(:)=0.0d0
  allocate(rwork_zheevd(LRWORK));rwork_zheevd(:)=0.0d0
  allocate(iwork_zheevd(LIWORK));iwork_zheevd(:)=0
  eps=1.0d-18
  ind=0                 
  !
  call zheevd("V","U",nm,mat,nm,eig,work_zheevd,LWORK,rwork_zheevd,LRWORK,iwork_zheevd,LIWORK,ind)
  !
  if(ind/=0)then 
   write(6,*)'ind=',ind 
   stop
  endif 
  !
  deallocate(work_zheevd,rwork_zheevd,iwork_zheevd) 
  return 
end subroutine
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
