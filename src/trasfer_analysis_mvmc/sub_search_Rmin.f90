subroutine search_Rmin(i,j,k,nkb1,nkb2,nkb3,a1,a2,a3,imin,jmin,kmin)
implicit none
integer::i,j,k,nkb1,nkb2,nkb3
integer::imin,jmin,kmin
real(8)::a1(3),a2(3),a3(3) 
integer::nmin,mmin,lmin
integer::n,m,l 
real(8)::R_pos(3),R_abs,R_min,R_bfr
R_pos(:)=dble(i)*a1(:)+dble(j)*a2(:)+dble(k)*a3(:)
nmin=0;mmin=0;lmin=0
R_bfr=dsqrt(R_pos(1)**2+R_pos(2)**2+R_pos(3)**2) 
R_min=R_bfr 
do n=-3,3
 do m=-3,3
  do l=-3,3
   R_pos(:)=dble(i+n*nkb1)*a1(:)+dble(j+m*nkb2)*a2(:)+dble(k+l*nkb3)*a3(:)
   R_abs=dsqrt(R_pos(1)**2+R_pos(2)**2+R_pos(3)**2) 
   if(R_min>R_abs)then 
    R_min=R_abs
    nmin=n
    mmin=m
    lmin=l
   endif 
  enddo
 enddo
enddo
imin=i+nmin*nkb1
jmin=j+mmin*nkb2
kmin=k+lmin*nkb3
!if(nmin/=0.or.mmin/=0.or.lmin/=0)then 
! write(6,'(a15,f15.8)')'R_before      ',R_bfr 
! write(6,'(a15,f15.8)')'R_after       ',R_min 
! write(6,'(a15,3i5)')  'i   ,j   ,k   ',i,j,k
! write(6,'(a15,3i5)')  'nmin,mmin,lmin',nmin,mmin,lmin
! write(6,'(a15,3i5)')  'imin,jmin,kmin',imin,jmin,kmin
!endif 
return
end
