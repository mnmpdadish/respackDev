module m_rdinput 
implicit none
private
public::read_input 
!--
integer,public::calc_ifreq!    
integer,public::ix_intJ_min!
integer,public::ix_intJ_max!
integer,public::iy_intJ_min!
integer,public::iy_intJ_max!
integer,public::iz_intJ_min!
integer,public::iz_intJ_max!
!real(8),public::wcut_model!Cutoff of U' for model 
!real(8),public::jcut_model!Cutoff of J  for model 
namelist/param_calc_int/calc_ifreq,ix_intJ_min,ix_intJ_max,iy_intJ_min,iy_intJ_max,&
iz_intJ_min,iz_intJ_max!,wcut_model,jcut_model   
contains
subroutine read_input 
!--
!default
CALC_IFREQ=1!omega=0
IX_INTJ_MIN=0
IX_INTJ_MAX=0
IY_INTJ_MIN=0
IY_INTJ_MAX=0
IZ_INTJ_MIN=0
IZ_INTJ_MAX=0 
!WCUT_MVMC=2.0d0!eV
!JCUT_MVMC=0.3d0!eV   
!--
!open(999,file='input.in')
!read(999,nml=param_calc_int)
read(5,nml=param_calc_int)
write(6,param_calc_int) 
end subroutine
end module
