!
!ESTIMATE nkb1,nkb2,nkb3 
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
  write(6,'(a24,4i10)') 'nkb1,nkb2,nkb3,NTK=',nkb1,nkb2,nkb3,NTK  
return 
end 
