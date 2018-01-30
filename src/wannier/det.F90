subroutine calcdet(N,mat,det) 
implicit none
complex(8)::det
integer,intent(in)::N 
complex(8),intent(inout)::mat(N,N) 
integer::ipiv(N)
integer::i,info 
!---
ipiv=0
info=0 
call zgetrf(N,N,mat,N,ipiv,info)
if(info/=0) then 
 write(6,*) 'info=',info 
endif 
det=1.0d0 
do i=1,N!-1
 if(ipiv(i)==i) then
  det=det*mat(i,i)
 else 
  det=-det*mat(i,i)
 endif  
enddo
return 
end
