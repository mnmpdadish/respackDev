subroutine  est_latparam(a1,a2,a3,a,b,c,alp,bet,gmm)   
implicit none 
real(8)::a1(3),a2(3),a3(3) 
real(8)::a,b,c,alp,bet,gmm  
real(8),parameter::pi=DACOS(-1.0d0)
!---
!read(100,*)aa(1,1),aa(2,1),aa(3,1)
!read(100,*)aa(1,2),aa(2,2),aa(3,2)
!read(100,*)aa(1,3),aa(2,3),aa(3,3)
!---
a=0.0d0 
b=0.0d0 
c=0.0d0 
a=dsqrt(a1(1)**2+a1(2)**2+a1(3)**2) 
b=dsqrt(a2(1)**2+a2(2)**2+a2(3)**2) 
c=dsqrt(a3(1)**2+a3(2)**2+a3(3)**2) 
alp=(a2(1)*a3(1)+a2(2)*a3(2)+a2(3)*a3(3))/b/c 
bet=(a3(1)*a1(1)+a3(2)*a1(2)+a3(3)*a1(3))/c/a
gmm=(a1(1)*a2(1)+a1(2)*a2(2)+a1(3)*a2(3))/a/b 
alp=dacos(alp)*180.0d0/pi  
bet=dacos(bet)*180.0d0/pi  
gmm=dacos(gmm)*180.0d0/pi  
!write(6,'(6f15.10)') a,b,c,alp,bet,gmm 
return 
end 
