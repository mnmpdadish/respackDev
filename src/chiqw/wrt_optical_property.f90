subroutine wrt_optical_property(ngrd,grd,func) 
  implicit none
  integer,intent(in)::ngrd
  complex(8),intent(in)::grd(ngrd)
  complex(8),intent(in)::func(ngrd,3)
  integer::ix,ie 
  real(8)::x,y,r  
  complex(8)::zsqrt 
  complex(8)::z 
  integer::file_num 
  character(99)::filename 
  real(8),parameter::au=27.21151D0!hartree 
  real(8),parameter::pi=dacos(-1.0d0)
  real(8),parameter::tpi=2.0d0*pi 
  real(8),parameter::unit_for_optical_conductivity=11.12650d0/2.41004d0![10^6Ohm^-1m^-1] 
  real(8),parameter::unit_for_optical_conductivity_secinv=41.341379d0![10^15sec-1] 
  integer,parameter::file_num_eels=9000 
  integer,parameter::file_num_optical=9100 
  integer,parameter::file_num_reflectivity=9200 
  integer,parameter::file_num_macroscopic_epsilon=9300 
  !
  !eels
  !
  do ix=1,3
   file_num=file_num_eels+ix
   !--
   !OPEN(900%,W,FILE='dat.eels-x,y,z')
   if(ix==1)write(filename,"('dat.eels-x')") 
   if(ix==2)write(filename,"('dat.eels-y')") 
   if(ix==3)write(filename,"('dat.eels-z')") 
   OPEN(file_num,FILE=filename) 
   REWIND(file_num)
   do ie=1,ngrd
    x=dble(func(ie,ix))
    y=imag(func(ie,ix)) 
    write(file_num,'(f15.8,2f25.12)') dble(grd(ie))*au,x,-y 
    !write(file_num,'(f15.8,2f25.12)') dble(grd(ie))*au,-y 
   enddo!ie
  enddo!ix 
  !
  !optical conductivity 
  !
  do ix=1,3
   file_num=file_num_optical+ix
   !--
   !OPEN(910%,W,FILE='dat.optical_conductivity-x,y,z')
   if(ix==1)write(filename,"('dat.optical_conductivity-x')") 
   if(ix==2)write(filename,"('dat.optical_conductivity-y')") 
   if(ix==3)write(filename,"('dat.optical_conductivity-z')") 
   OPEN(file_num,FILE=filename) 
   REWIND(file_num)
   do ie=1,ngrd
    !
    !20191021 Kazuma Nakamura 
    !
    !x=dble(func(ie,ix))
    !y=imag(func(ie,ix)) 
    !write(file_num,'(f15.8,2f25.12)') dble(grd(ie))*au,x/(x**2+y**2),-y/(x**2+y**2) 
    z=func(ie,ix) 
    !
    !unit: 10^6 Ohm-1 m-1
    !
    !write(file_num,'(f15.8,2f25.12)') dble(grd(ie))*au,& 
    !dble(grd(ie))/(2.0d0*tpi)*dble(1.0d0/z)*unit_for_optical_conductivity,&   
    !dble(grd(ie))/(2.0d0*tpi)*imag(1.0d0/z)*unit_for_optical_conductivity  
    !
    !unit: 10^15 sec-1
    !
    write(file_num,'(f15.8,3f25.12)') dble(grd(ie))*au,&
    dble(grd(ie))/(2.0d0*tpi)*dble(1.0d0/z)*unit_for_optical_conductivity_secinv,&   
    dble(grd(ie))/(2.0d0*tpi)*imag(1.0d0/z)*unit_for_optical_conductivity_secinv
    !
   enddo!ie
  enddo!ix 
  !
  !reflectivity 
  !
  do ix=1,3
   file_num=file_num_reflectivity+ix
   !--
   !OPEN(920%,W,FILE='dat.reflectivity-x,y,z')
   if(ix==1)write(filename,"('dat.reflectivity-x')") 
   if(ix==2)write(filename,"('dat.reflectivity-y')") 
   if(ix==3)write(filename,"('dat.reflectivity-z')") 
   OPEN(file_num,FILE=filename) 
   REWIND(file_num)
   do ie=1,ngrd
    zsqrt=func(ie,ix)**0.5d0
    r=(abs((1.0d0-zsqrt)/(1.0d0+zsqrt)))**2 
    write(file_num,'(2f15.8)') dble(grd(ie))*au,r 
   enddo!ie 
  enddo!ix  
  !
  !macroscpic dielectric function 
  !
  do ix=1,3
   file_num=file_num_macroscopic_epsilon+ix
   !--
   !OPEN(930%,W,FILE='dat.macroscopic_epsilon-x,y,z')
   if(ix==1)write(filename,"('dat.macroscopic_epsilon-x')") 
   if(ix==2)write(filename,"('dat.macroscopic_epsilon-y')") 
   if(ix==3)write(filename,"('dat.macroscopic_epsilon-z')") 
   OPEN(file_num,FILE=filename) 
   REWIND(file_num)
   do ie=1,ngrd
    !
    !20191021 Kazuma Nakamura 
    !
    !x=dble(func(ie,ix))
    !y=imag(func(ie,ix)) 
    !write(file_num,'(f15.8,2f25.12)') dble(grd(ie))*au,x/(x**2+y**2),-y/(x**2+y**2) 
    z=func(ie,ix) 
    write(file_num,'(f15.8,2f25.12)') dble(grd(ie))*au,dble(1.0d0/z),imag(1.0d0/z) 
    !
   enddo!ie
  enddo!ix 
  !
return
end
