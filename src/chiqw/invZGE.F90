subroutine invZGE(nm,mat)
  implicit none 
  integer,intent(in)::nm
  complex(8),intent(inout)::mat(nm,nm)
  integer::ipiv(nm)
  integer::Lwork 
  complex(8),ALLOCATABLE::work(:)
  integer::info 
  Lwork=10*nm
  allocate(work(Lwork))
  info=0
  call zgetrf(nm,nm,mat,nm,ipiv,info)
!  write(6,*)'finish zgetrf'
  call zgetri(nm,mat,nm,ipiv,work,Lwork,info)
!  write(6,*)'finish zgetri'
  if(info/=0) then
    write(6,*)'info (subrouitine inv):',info
    stop
  end if 
  deallocate(work)
  return 
end subroutine

