subroutine est_NTK(N,SK,NTK,nkb1,nkb2,nkb3)  
implicit none 
integer::N,nkb1,nkb2,nkb3,NTK  
real(8)::SK(3,N) 
integer::i 
real(8)::x 
!---
x=1.0d0 
do i=1,N
 if(abs(SK(1,i))<1.0d-7) cycle 
 if(abs(SK(1,i))<x) then 
  x=abs(SK(1,i))  
 endif 
enddo    
nkb1=nint(1.0d0/x)  
!---
x=1.0d0 
do i=1,N
 if(abs(SK(2,i))<1.0d-7) cycle 
 if(abs(SK(2,i))<x) then 
  x=abs(SK(2,i))  
 endif 
enddo    
nkb2=nint(1.0d0/x)  
!---
x=1.0d0 
do i=1,N
 if(abs(SK(3,i))<1.0d-7) cycle 
 if(abs(SK(3,i))<x) then 
  x=abs(SK(3,i))  
 endif 
enddo    
nkb3=nint(1.0d0/x)  
!---
NTK=nkb1*nkb2*nkb3 
!write(6,'(a24,4i10)') 'nkb1,nkb2,nkb3,NTK=',nkb1,nkb2,nkb3,NTK  
!---
return 
end 
