module m_wrt_frmsf 
  implicit none
contains
subroutine wrt_frmsf(NTB,NTK,nkb1,nkb2,nkb3,E_EIG,SK0,FermiEnergy,b1,b2,b3)                     
  implicit none 
  integer::NTB,NTK,nkb1,nkb2,nkb3 
  real(8)::E_EIG(NTB,NTK)           
  real(8)::SK0(3,NTK)           
  integer::index_kpt(nkb1,nkb2,nkb3) 
  integer::ib,ik,ikb1,ikb2,ikb3
  real(8)::FermiEnergy,b1(3),b2(3),b3(3)
  real(8)::E_3D(nkb3,nkb2,nkb1,NTB) 
  integer::fo,ishift  
  index_kpt=0     
  E_3D=0.0d0 
  call make_index_kpt(NTK,nkb1,nkb2,nkb3,SK0(1,1),index_kpt(1,1,1)) 
  do ib=1,NTB 
   do ikb1=1,nkb1
    do ikb2=1,nkb2
     do ikb3=1,nkb3 
      ik=index_kpt(ikb1,ikb2,ikb3)
      E_3D(ikb3,ikb2,ikb1,ib)=E_EIG(ib,ik)-FermiEnergy 
     enddo 
    enddo 
   enddo 
  enddo 
  ! 
  ishift=0 
  fo=301
  OPEN(fo,FILE='./dir-wan/dat.frmsf') 
  rewind(fo) 
  write(fo,*) nkb1,nkb2,nkb3
  write(fo,*) ishift!not grid shift
  write(fo,*) NTB 
  write(fo,*) real(b1(1:3))
  write(fo,*) real(b2(1:3))
  write(fo,*) real(b3(1:3))
  do ib=1,NTB 
   do ikb1=1,nkb1
    do ikb2=1,nkb2
     do ikb3=1,nkb3
      write(fo,*) real(E_3D(ikb3,ikb2,ikb1,ib))
     end do
    end do
   end do
  end do
  do ib=1,NTB
   do ikb1=1,nkb1
    do ikb2=1,nkb2
     do ikb3=1,nkb3
      write(fo,*) real(ib)!real(phys(ikb3,ikb2,ikb1,ib))
     end do
    end do
   end do
  end do
  close(fo)
  return
end subroutine 
!
SUBROUTINE make_index_kpt(NTK,nkb1,nkb2,nkb3,SK0,index_kpt)      
  implicit none 
  integer::NTK,nkb1,nkb2,nkb3
  real(8)::SK0(3,NTK)  
  integer::ik,ix,iy,iz
  real(8)::x,y,z
  integer::index_kpt(nkb1,nkb2,nkb3)    
 
  if(MOD(NTK,2)/=0) then 
   do ik=1,NTK 
    x=SK0(1,ik)*dble(nkb1) 
    y=SK0(2,ik)*dble(nkb2)
    z=SK0(3,ik)*dble(nkb3)  
    x=x+(dble(nkb1)-1.0d0)/2.0d0 
    y=y+(dble(nkb2)-1.0d0)/2.0d0
    z=z+(dble(nkb3)-1.0d0)/2.0d0 
    ix=idnint(x)+1
    iy=idnint(y)+1
    iz=idnint(z)+1
    index_kpt(ix,iy,iz)=ik
   enddo 
  else!20170316 
   do ik=1,NTK 
    x=SK0(1,ik)*dble(nkb1) 
    y=SK0(2,ik)*dble(nkb2)
    z=SK0(3,ik)*dble(nkb3)  
    x=x+dble(nkb1)/2.0d0 
    y=y+dble(nkb2)/2.0d0
    z=z+dble(nkb3)/2.0d0 
    ix=idnint(x)
    iy=idnint(y)
    iz=idnint(z)
    index_kpt(ix,iy,iz)=ik
   enddo 
  endif 
  RETURN  
end subroutine 
  !
end module m_wrt_frmsf 
