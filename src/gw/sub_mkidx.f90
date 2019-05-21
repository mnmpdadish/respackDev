SUBROUTINE make_index_kpt(NTK,nkb1,nkb2,nkb3,SK0,index_kpt)      
  implicit none 
  integer::NTK,nkb1,nkb2,nkb3
  real(8)::SK0(3,NTK)  
  integer::ik,ix,iy,iz
  real(8)::x,y,z
  integer::index_kpt(nkb1,nkb2,nkb3)    
  ! 
  !if(MOD(NTK,2)/=0)then 
  ! !write(6,*)'i am in make_index for odd'
  ! do ik=1,NTK 
  !  x=SK0(1,ik)*dble(nkb1) 
  !  y=SK0(2,ik)*dble(nkb2)
  !  z=SK0(3,ik)*dble(nkb3)  
  !  x=x+(dble(nkb1)-1.0d0)/2.0d0 
  !  y=y+(dble(nkb2)-1.0d0)/2.0d0
  !  z=z+(dble(nkb3)-1.0d0)/2.0d0 
  !  ix=idnint(x)+1
  !  iy=idnint(y)+1
  !  iz=idnint(z)+1
  !  index_kpt(ix,iy,iz)=ik
  ! enddo 
  !else!20170316 
  ! !write(6,*)'i am in make_index for even'
  ! do ik=1,NTK 
  !  x=SK0(1,ik)*dble(nkb1) 
  !  y=SK0(2,ik)*dble(nkb2)
  !  z=SK0(3,ik)*dble(nkb3)  
  !  x=x+dble(nkb1)/2.0d0 
  !  y=y+dble(nkb2)/2.0d0
  !  z=z+dble(nkb3)/2.0d0 
  !  ix=idnint(x)
  !  iy=idnint(y)
  !  iz=idnint(z)
  !  index_kpt(ix,iy,iz)=ik
  ! enddo 
  !endif 
  !--
  !
  !20190520 Kazuma Nakamura
  !
  do ik=1,NTK 
   x=SK0(1,ik)*dble(nkb1) 
   y=SK0(2,ik)*dble(nkb2)
   z=SK0(3,ik)*dble(nkb3)  
   x=x+(dble(nkb1)-dble(mod(nkb1,2)))/2.0d0 
   y=y+(dble(nkb2)-dble(mod(nkb2,2)))/2.0d0 
   z=z+(dble(nkb3)-dble(mod(nkb3,2)))/2.0d0 
   ix=idnint(x)+mod(nkb1,2)
   iy=idnint(y)+mod(nkb2,2)
   iz=idnint(z)+mod(nkb3,2)
   index_kpt(ix,iy,iz)=ik
  enddo 
  ! 
  RETURN  
END SUBROUTINE make_index_kpt  
