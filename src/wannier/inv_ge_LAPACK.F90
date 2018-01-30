subroutine inv_ge_lapack(L,NBh,mat) 
implicit none 
integer,intent(in)::L 
integer,intent(in)::NBh 
real(8),intent(inout)::mat(L,L)!notice 
integer,allocatable::ipiv(:)!ipiv(NBh)
real(8),allocatable::work(:)
integer::Lwork,info,m,n,lda  
m=NBh  
n=L 
lda=L   
Lwork=10*m
allocate(ipiv(NBh));ipiv=0 
allocate(work(Lwork));work=0.0d0 
info=0
call dgetrf(n,m,mat,lda,ipiv,info)
if(info/=0)then
 write(6,*)'info(subrouitine dgetrf):',info
!stop
endif 
call dgetri(m,mat,lda,ipiv,work,Lwork,info)
if(info/=0)then
 write(6,*)'info(subrouitine dgetri):',info
!stop
endif 
deallocate(work)
!write(6,*) mat 
return 
end subroutine
