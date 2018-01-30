program mksupercell !by Kazum NAKAMURA 
implicit none 
real(8),parameter::pi=DACOS(-1.0d0)
integer,parameter::nx=2 
integer,parameter::ny=2 
integer,parameter::nz=2 
integer,parameter::nxyz=nx*ny*nz 
character(2),allocatable::kdp(:)!kd(NTA)
character(2),allocatable::kd(:)!kd(nxyz*NTA)
real(8),allocatable::xp(:),yp(:),zp(:)!xp(NTA),yp(NTA),zp(NTA)  
real(8),allocatable::x(:),y(:),z(:)!x(nxyz*NTA),y(nxyz*NTA),z(nxyz*NTA)  
integer::NTA
real(8)::a1(3),a2(3),a3(3),alat,talat 
integer::i,j,n,ip,il,jl,kl  
real(8)::R_MIN,R,a,b,c,beta  
!--
read(100,*) alat 
read(100,*) a1(1),a1(2),a1(3)
read(100,*) a2(1),a2(2),a2(3)
read(100,*) a3(1),a3(2),a3(3)
read(100,*) NTA
allocate(kdp(NTA))
allocate(xp(NTA)) 
allocate(yp(NTA)) 
allocate(zp(NTA)) 
!--
do i=1,NTA
 read(100,*) kdp(i),xp(i),yp(i),zp(i) 
enddo 
!--
xp(:)=xp(:)/dble(nx) 
yp(:)=yp(:)/dble(ny) 
zp(:)=zp(:)/dble(nz) 
!--
allocate(kd(nxyz*NTA))
allocate(x(nxyz*NTA)) 
allocate(y(nxyz*NTA)) 
allocate(z(nxyz*NTA)) 
!--
do kl=1,nz
 do jl=1,ny
  do il=1,nx
   do i=1,NTA
    ip=i+(il-1)*NTA+(jl-1)*NTA*nx+(kl-1)*NTA*nx*ny 
    x(ip)=xp(i)+dble(il-1)/dble(nx) 
    y(ip)=yp(i)+dble(jl-1)/dble(ny) 
    z(ip)=zp(i)+dble(kl-1)/dble(nz) 
    kd(ip)=kdp(i)
   enddo 
  enddo 
 enddo 
enddo 
!--
write(6,'(f15.10)') alat 
write(6,*) a1(1)*dble(nx),a1(2)*dble(nx),a1(3)*dble(nx)
write(6,*) a2(1)*dble(ny),a2(2)*dble(ny),a2(3)*dble(ny)
write(6,*) a3(1)*dble(nz),a3(2)*dble(nz),a3(3)*dble(nz) 
write(6,*) NTA*nxyz 
do i=1,NTA*nxyz 
 write(6,'(a5,3f15.10)') kd(i),x(i),y(i),z(i) 
enddo 
!--
STOP
END  
