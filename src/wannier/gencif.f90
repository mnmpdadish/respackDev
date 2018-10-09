!
! Copyright (c) 2013-2018 Yoshihide Yoshimoto and Kazuma Nakamura 
!
subroutine gencif(a1,a2,a3,amin,amax,bmin,bmax,cmin,cmax,nsi,cs,asi)    
use m_gencif_sub!subr_fmtconv: original name in xtapp 
implicit none 
real(8),intent(in)::a1(3),a2(3),a3(3)
integer,intent(in)::amin,amax,bmin,bmax,cmin,cmax 
!
integer::nsi 
character(len=2)::cs(nsi) 
real(8)::asi(3,nsi) 
!
integer::nkd 
character(len=2)::csr 
character(len=2),allocatable::csref(:)!csref(nkd) 
integer,allocatable::kd(:)!kd(nsi) 
real(8),allocatable::zn(:)!zn(nkd)
!
integer,allocatable::kd_sc(:) 
real(8),allocatable::asi_sc(:,:) 
real(8)::aa(3,3),aa_sc(3,3)  
integer::i,is,js,ks,ip,j,ia  
integer::na,nb,nc,nsi_sc 
character(len=256)::cbuf
! 
allocate(kd(nsi));kd=0 
csr=cs(1)
is=1
kd(1)=is 
do ia=2,nsi
 if(csr/=cs(ia))then
  csr=cs(ia)
  is=is+1
  kd(ia)=is 
 else
  kd(ia)=is 
 endif 
enddo 
nkd=is 
!
allocate(csref(nkd))
csr=cs(1)
is=1
csref(is)=csr 
do ia=2,nsi
 if(csr/=cs(ia))then
  csr=cs(ia)
  is=is+1
  csref(is)=csr  
 endif 
enddo 
!
allocate(zn(nkd));zn=0.0d0  
do i=1,nkd 
 do is=0,102
  if(csref(i).eq.atom_name(is))then
   js=is 
  endif 
 enddo 
 zn(i)=dble(js) 
enddo 
!
aa(:,1)=a1(:)
aa(:,2)=a2(:)
aa(:,3)=a3(:)
!
na=amax-amin
nb=bmax-bmin
nc=cmax-cmin
nsi_sc=nsi*na*nb*nc 
!
allocate(kd_sc(nsi_sc));kd_sc=0 
allocate(asi_sc(3,nsi_sc));asi_sc=0.0d0
! 
asi(1,:)=asi(1,:)/dble(na) 
asi(2,:)=asi(2,:)/dble(nb) 
asi(3,:)=asi(3,:)/dble(nc) 
! 
do is=1,na
 do js=1,nb
  do ks=1,nc 
   do i=1,nsi 
    ip=i+nsi*(ks-1)+nsi*nc*(js-1)+nsi*nc*nb*(is-1)
    asi_sc(1,ip)=asi(1,i)+dble((is-1))/dble(na) 
    asi_sc(2,ip)=asi(2,i)+dble((js-1))/dble(nb) 
    asi_sc(3,ip)=asi(3,i)+dble((ks-1))/dble(nc) 
    kd_sc(ip)=kd(i)
   enddo!i
  enddo!ks
 enddo!js
enddo!is 
!
aa_sc(:,1)=aa(:,1)*dble(na)
aa_sc(:,2)=aa(:,2)*dble(nb)
aa_sc(:,3)=aa(:,3)*dble(nc)
!
!do j=1,3
! write(6,'(3f15.10)')(aa_sc(i,j),i=1,3) 
!enddo 
!
!do i=1,nsi_sc
! write(6,'(i5,3f15.10)')kd_sc(i),(asi_sc(j,i),j=1,3) 
!enddo 
!
cbuf='cif'
call printcif(cbuf,aa_sc(1,1),nkd,zn(1),nsi_sc,kd_sc(1),asi_sc(1,1),na,nb,nc)
!
return 
end 

