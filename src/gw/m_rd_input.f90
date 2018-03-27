module m_rdinput 
implicit none
private
public::read_input 
!&param_calc_gw 
integer,public::calc_ifreq!    
integer,public::ix_intJ_min!
integer,public::ix_intJ_max!
integer,public::iy_intJ_min!
integer,public::iy_intJ_max!
integer,public::iz_intJ_min!
integer,public::iz_intJ_max!
real(8),public::wcut_mvmc!Cutoff of U' for mvmc 
real(8),public::jcut_mvmc!Cutoff of J  for mvmc 
real(8),public::Green_func_delt!ttrhdrn Green's function delt (eV)
real(8),public::idlt!ttrhdrn Green's function delt (eV)
integer,public::Rc_range_spacing!Range of attenuation potential cutoff 
!&param_interpolation   
integer,public::N_sym_points!The number of k-point points in symmetry line
integer,public::Ndiv!Separation between symmetry points  
integer,public::reading_sk_format!20170709 
real(8),public,allocatable::SK_sym_pts(:,:) 
integer,public,allocatable::dense(:)!dense(3)!Dense k-grid for the Wnnier-interpolated FS
namelist/param_interpolation/N_sym_points,Ndiv,reading_sk_format,dense
namelist/param_calc_gw/calc_ifreq,ix_intJ_min,ix_intJ_max,iy_intJ_min,iy_intJ_max,&
iz_intJ_min,iz_intJ_max,wcut_mvmc,jcut_mvmc,Green_func_delt,Rc_range_spacing 
contains
subroutine read_input 
integer::ix,ik  
!--
!&param_interpolation
!--
!default
Ndiv=40!Separation between symmetry points  
allocate(dense(3));dense(1:3)=0
!--
rewind(5) 
!read(999,nml=param_interpolation)
read(5,nml=param_interpolation)
write(6,param_interpolation) 
!k for band dispersion
allocate(SK_sym_pts(3,N_sym_points))  
if(reading_sk_format==0)then 
 write(6,*)'READING SK_sym_pts format=0: respack'
 do ik=1,N_sym_points 
  !read(999,*)(SK_sym_pts(ix,ik),ix=1,3) 
  read(5,*)(SK_sym_pts(ix,ik),ix=1,3) 
 enddo 
endif 
if(reading_sk_format==1)then 
 write(6,*)'READING SK_sym_pts format=1: xtapp'
 do ix=1,3 
  !read(999,*)(SK_sym_pts(ix,ik),ik=1,N_sym_points) 
  read(5,*)(SK_sym_pts(ix,ik),ik=1,N_sym_points) 
 enddo 
endif 
!--
write(6,param_interpolation) 
do ik=1,N_sym_points 
 write(6,'(3f10.5)')(SK_sym_pts(ix,ik),ix=1,3) 
enddo 
write(6,*) 
!--
!&param_calc_gw 
!--
!default
CALC_IFREQ=1!omega=0
IX_INTJ_MIN=0
IX_INTJ_MAX=0
IY_INTJ_MIN=0
IY_INTJ_MAX=0
IZ_INTJ_MIN=0
IZ_INTJ_MAX=0 
WCUT_MVMC=2.0d0!eV
JCUT_MVMC=0.3d0!eV   
GREEN_FUNC_DELT=0.1d0!eV
Rc_range_spacing=2!3  
!--
!open(999,file='input.in')
!read(999,nml=param_calc_gw)
!rewind(5) 
read(5,nml=param_calc_gw)
write(6,param_calc_gw) 
idlt=Green_func_delt 
!close(5)
!--
end subroutine
end module
