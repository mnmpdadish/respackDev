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
!---
! L_wzgetri=10*nen 
! ALLOCATE(ipiv(nen)) 
! ALLOCATE(wzgetri(L_wzgetri))
! info=0
! call zgetrf(nen,nen,mat_b,nen,ipiv,info)
! call zgetri(nen,mat_b,nen,ipiv,wzgetri,L_wzgetri,info) 
! DEALLOCATE(ipiv,wzgetri)
!---
