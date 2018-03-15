subroutine make_LG0(NTG,b1,b2,b3,Gcut_for_eps,Gcut_for_psi,q1,q2,q3,LG0,NG_for_eps,NG_for_psi)
implicit none 
integer::NTG,igL,igL1,igL2,igL3
integer::NG_for_eps,NG_for_psi            
integer::LG0(3,NTG)    
real(8)::b1(3),b2(3),b3(3) 
real(8)::Gcut_for_eps,Gcut_for_psi  
real(8)::q1,q2,q3,qgL2,qgL(3)
integer,parameter::NGL1=100 
integer,parameter::NGL2=100 
integer,parameter::NGL3=100 
!--
igL=0
do igL1=-NGL1,NGL1 
 do igL2=-NGL2,NGL2 
  do igL3=-NGL3,NGL3 
   qgL(:)=(q1+dble(igL1))*b1(:)+(q2+dble(igL2))*b2(:)+(q3+dble(igL3))*b3(:)    
   qgL2=qgL(1)**2+qgL(2)**2+qgL(3)**2
   if(qgL2<=Gcut_for_eps)then 
    igL=igL+1 
    LG0(1,igL)=igL1
    LG0(2,igL)=igL2
    LG0(3,igL)=igL3 
    !write(6,*) igL
   endif  
  enddo 
 enddo 
enddo 
NG_for_eps=igL 
!--
do igL1=-NGL1,NGL1 
 do igL2=-NGL2,NGL2 
  do igL3=-NGL3,NGL3 
   qgL(:)=(q1+dble(igL1))*b1(:)+(q2+dble(igL2))*b2(:)+(q3+dble(igL3))*b3(:)    
   qgL2=qgL(1)**2+qgL(2)**2+qgL(3)**2
   if(qgL2>Gcut_for_eps.and.qgL2<=Gcut_for_psi)then 
    igL=igL+1 
    LG0(1,igL)=igL1
    LG0(2,igL)=igL2
    LG0(3,igL)=igL3 
    !write(6,*) igL
   endif  
  enddo 
 enddo 
enddo 
NG_for_psi=igL 
!--
!Do igL=1,NG_for_eps 
!write(6,*) igL,LG0(:,igL)
!enddo 
!--
RETURN 
END 
