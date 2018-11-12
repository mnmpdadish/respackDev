module m_rd_input 
implicit none
private
public::read_input 
!&param_calc_gw 
integer,public::Rc_range_spacing!Range of attenuation potential cutoff 
integer,public::Ncalc!The number of bands considered in the GW calculation  
logical,public::calc_SC!flag to calc SC or not 
real(8),public::gw_grid_separation!minimum separation of GW grid (eV) 
!real(8),public::Green_func_delt!ttrhdrn Green's function delt (eV)
!real(8),public::idlt!ttrhdrn Green's function delt (eV)
!&param_interpolation   
integer,public::N_sym_points!The number of k-point points in symmetry line
integer,public::Ndiv!Separation between symmetry points  
integer,public::reading_sk_format!20170709 
real(8),public,allocatable::SK_sym_pts(:,:) 
integer,public,allocatable::dense(:)!dense(3)!Dense k-grid for the Wnnier-interpolated FS
namelist/param_interpolation/N_sym_points,Ndiv,reading_sk_format,dense
namelist/param_calc_gw/Rc_range_spacing,Ncalc,calc_SC,gw_grid_separation   
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
!
!GREEN_FUNC_DELT=0.1d0!eV
!
Rc_range_spacing=2!3  
Ncalc=0!Ncalc is set to NTB after
calc_SC=.true.!calc SC: .true., not calc SC: .false.
gw_grid_separation=0.05d0!(eV) 
!--
!open(999,file='input.in')
!read(999,nml=param_calc_gw)
read(5,nml=param_calc_gw)
write(6,param_calc_gw) 
!
!idlt=Green_func_delt 
!
!--
end subroutine
end module
