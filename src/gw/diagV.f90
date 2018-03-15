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
!---
LWORK= 2*nm+nm**2
LRWORK=1+12*nm+3*nm**2
LIWORK=3+10*nm 
allocate(work_zheevd(LWORK));work_zheevd(:)=0.0d0
allocate(rwork_zheevd(LRWORK));rwork_zheevd(:)=0.0d0
allocate(iwork_zheevd(LIWORK));iwork_zheevd(:)=0
eps=1.0d-18
ind=0                 
!---
call zheevd("V","U",nm,mat,nm,eig,work_zheevd,LWORK, &
        rwork_zheevd,LRWORK,iwork_zheevd,LIWORK,ind)
!---
if(ind/=0) then 
write(6,*) 'ind=',ind 
stop
endif 
!---
deallocate(work_zheevd,rwork_zheevd,iwork_zheevd) 
return 
end subroutine
