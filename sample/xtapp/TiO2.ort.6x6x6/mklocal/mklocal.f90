program mklocal !by Kazum NAKAMURA 
implicit none 
integer::NTA 
character(2),allocatable::kd(:)
real,allocatable::lat(:,:)
integer::i,ix 
real(8)::a1(3),a2(3),a3(3),alat 
real(8)::locx(3),locy(3),locz(3),length,inp  
!--
read(100,*) alat 
read(100,*) a1(1),a1(2),a1(3)
read(100,*) a2(1),a2(2),a2(3)
read(100,*) a3(1),a3(2),a3(3)
read(100,*) NTA
allocate(kd(NTA))
allocate(lat(3,NTA)) 
do i=1,NTA
 read(100,*) kd(i),(lat(ix,i),ix=1,3) 
enddo 
!--
!local for TiO2 
 call gen_loc_1_TiO2(NTA,kd(1),lat(1,1),alat,a1(1),a2(1),a3(1),locx(1),locy(1),locz(1)) 
 write(6,'(9f10.5)')locx,locy,locz 
!--
!local for TiO2 
 call gen_loc_2_TiO2(NTA,kd(1),lat(1,1),alat,a1(1),a2(1),a3(1),locx(1),locy(1),locz(1)) 
 write(6,'(9f10.5)')locx,locy,locz 
!--
stop
end
!--
SUBROUTINE OUTER_PRODUCT(vec_x,vec_y,vec_z)
implicit none 
real(8)::vec_x(3),vec_y(3),vec_z(3) 
vec_z(1)=vec_x(2)*vec_y(3)-vec_x(3)*vec_y(2)
vec_z(2)=vec_x(3)*vec_y(1)-vec_x(1)*vec_y(3) 
vec_z(3)=vec_x(1)*vec_y(2)-vec_x(2)*vec_y(1)
RETURN
END
!--
subroutine gen_loc_1_TiO2(NTA,kd,lat,alat,a1,a2,a3,locx,locy,locz) 
implicit none 
integer::NTA 
character(2)::kd(3)
real::lat(3,NTA)
integer::i,ix 
real(8)::a1(3),a2(3),a3(3),alat 
real(8)::x1(3),x2(3),x3(3) 
real(8)::locx(3),locy(3),locz(3),length  
x1=0.0d0;x2=0.0d0;x3=0.0d0 
locx=0.0d0;locy=0.0d0;locz=0.0d0 
!--
do i=1,NTA 
 if(i==1)then
  !write(6,*) kd(i),(lat(ix,i),ix=1,3) 
  x1(:)=lat(1,i)*alat*a1(:)+lat(2,i)*alat*a2(:)+lat(3,i)*alat*a3(:)  
 endif
 if(i==6)then 
  !write(6,*) kd(i),(lat(ix,i),ix=1,3) 
  x2(:)=lat(1,i)*alat*a1(:)+(lat(2,i)-1.0d0)*alat*a2(:)+lat(3,i)*alat*a3(:)  
  x3(:)=lat(1,i)*alat*a1(:)+(lat(2,i)-1.0d0)*alat*a2(:)+(lat(3,i)-1.0d0)*alat*a3(:)  
 endif 
enddo 
!--
locx=x2-x1 
length=dsqrt(locx(1)**2+locx(2)**2+locx(3)**2)
locx=locx/length 
!--
locy(:)=x3-x1
length=dsqrt(locy(1)**2+locy(2)**2+locy(3)**2)
locy=locy/length 
!--
call OUTER_PRODUCT(locx(1),locy(1),locz(1))
length=dsqrt(locz(1)**2+locz(2)**2+locz(3)**2)
locz=locz/length 
return
end
!--
subroutine gen_loc_2_TiO2(NTA,kd,lat,alat,a1,a2,a3,locx,locy,locz) 
implicit none 
integer::NTA 
character(2)::kd(3)
real::lat(3,NTA)
integer::i,ix 
real(8)::a1(3),a2(3),a3(3),alat 
real(8)::x1(3),x2(3),x3(3) 
real(8)::locx(3),locy(3),locz(3),length  
x1=0.0d0;x2=0.0d0;x3=0.0d0 
locx=0.0d0;locy=0.0d0;locz=0.0d0 
!--
do i=1,NTA 
 if(i==2)then
  !write(6,*) kd(i),(lat(ix,i),ix=1,3) 
  x1(:)=lat(1,i)*alat*a1(:)+lat(2,i)*alat*a2(:)+lat(3,i)*alat*a3(:)  
 endif
 if(i==3)then 
  !write(6,*) kd(i),(lat(ix,i),ix=1,3) 
  x2(:)=lat(1,i)*alat*a1(:)+lat(2,i)*alat*a2(:)+(lat(3,i)+1.0d0)*alat*a3(:)  
 endif 
 if(i==4)then 
  x3(:)=lat(1,i)*alat*a1(:)+lat(2,i)*alat*a2(:)+(lat(3,i)+1.0d0)*alat*a3(:)  
 endif 
enddo 
!--
locx=x2-x1 
length=dsqrt(locx(1)**2+locx(2)**2+locx(3)**2)
locx=locx/length 
!--
locy(:)=x3-x1
length=dsqrt(locy(1)**2+locy(2)**2+locy(3)**2)
locy=locy/length 
!--
call OUTER_PRODUCT(locx(1),locy(1),locz(1))
length=dsqrt(locz(1)**2+locz(2)**2+locz(3)**2)
locz=locz/length 
return
end
