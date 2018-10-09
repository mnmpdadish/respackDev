subroutine gen_grid(Nt,E,em)  
implicit none 
real(8),parameter::au=27.21151d0
real(8),parameter::dlt=0.05d0/au!default  
real(8),parameter::ratio=0.1d0!default  
real(8),parameter::eps=1.0d-6!default 
real(8)::alp,x,x1,dl,a,b,c,d,E,f,g  
integer::Nt,N,M,i,k 
complex(8)::em(Nt) 
if(Nt==1)then 
 em(1)=0.0d0 
else
 M=int(ratio*dble(Nt)) 
 N=Nt-M
 alp=E/dlt
 x=1.5d0 
 dl=1.0d0 
 k=0
 do while(dl>eps) 
 f=x**(N-1)-alp*x+alp-1.0d0
 g=dble(N-1)*x**(N-2)-alp 
 x1=x-f/g
 dl=dabs(x1-x)
 k=k+1
 x=x1 
 enddo 
 a=dlt/(x-1.0d0)
 b=log(x)
 do i=1,N
 em(i)=dcmplx(a*exp(b*dble(i-1))-a,0.0d0) 
 enddo 
 !--
 !20180906 
 !c=E*exp(-log(2.0d0)*dble(N)/dble(M)) 
 !d=log(2.0d0)/dble(M) 
 c=E*exp(-log(3.0d0)*dble(N)/dble(M)) 
 d=log(3.0d0)/dble(M) 
 !--
 do i=N+1,N+M 
 em(i)=dcmplx(c*exp(d*dble(i)),0.0d0)  
 enddo 
endif 
return 
end
