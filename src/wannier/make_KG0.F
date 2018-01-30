!--
      call make_KG0(NTG,b1(1),b2(1),b3(1),Gcut_for_psi,ktmp(1),ktmp(2),ktmp(3),KGtmp(1,1),NG_for_psi)
!--
      subroutine make_KG0(NTG,b1,b2,b3,Gcut,q1,q2,q3,KG0,NG) 
      implicit none 
      integer,intent(in)::NTG
      real(8),intent(in)::b1(3),b2(3),b3(3) 
      real(8),intent(in)::Gcut,q1,q2,q3        
      integer,intent(out)::KG0(3,NTG),NG 
      integer::ig,igL,igL1,igL2,igL3
      real(8)::qgL(3),qgL2  
      integer,parameter::NGL1=50
      integer,parameter::NGL2=50 
      integer,parameter::NGL3=50

      igL=0
      do igL1=-NGL1,NGL1 
       do igL2=-NGL2,NGL2 
        do igL3=-NGL3,NGL3 
         qgL(:)=(q1+dble(igL1))*b1(:)+(q2+dble(igL2))*b2(:)+(q3+dble(igL3))*b3(:)    
         qgL2=qgL(1)**2+qgL(2)**2+qgL(3)**2
         if(qgL2<=Gcut) then 
          igL=igL+1 
          KG0(1,igL)=igL1
          KG0(2,igL)=igL2
          KG0(3,igL)=igL3 
         endif  
        enddo 
       enddo 
      enddo 
      NG=igL 
      RETURN 
      END 
!---
