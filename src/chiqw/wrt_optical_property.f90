subroutine wrt_optical_property(ngrd,grd,func) 
  implicit none
  integer,intent(in)::ngrd
  complex(8),intent(in)::grd(ngrd)
  complex(8),intent(in)::func(ngrd,3)
  integer::ix,ie 
  real(8)::x,y,r  
  complex(8)::zsqrt 
  integer::file_num 
  character(99)::filename 
  real(8),parameter::au=27.21151D0!hartree 
  integer,parameter::file_num_eels=9000 
  integer,parameter::file_num_optical=9100 
  integer,parameter::file_num_reflectivity=9200 
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
    x=dble(func(ie,ix))
    y=imag(func(ie,ix)) 
    write(file_num,'(f15.8,2f25.12)') dble(grd(ie))*au,x/(x**2+y**2),-y/(x**2+y**2) 
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
return
end
