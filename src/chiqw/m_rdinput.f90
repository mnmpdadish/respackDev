module m_rdinput 
implicit none
private
public::read_input 
!--
integer,public::N_CALC_BAND!TOTAL NUMBER OF CALCULATION BANDS 
real(8),public::Ecut_for_eps!cutoff energy for dielectric matrix (Ry)
real(8),public::shift_ef!artificial shift of Ef (eV) 
real(8),public::Max_excitation_energy!delta_ex!excitation constraint (eV)
real(8),public::delta_ex!excitation constraint (eV)
real(8),public::Green_func_delt!delt!ttrhdrn Green's function delt (eV)
real(8),public::delt!ttrhdrn Green's function delt (eV)
real(8),public::ttrhdrn_dmna!dmna!ttrhdrn
real(8),public::dmna!ttrhdrn
real(8),public::ttrhdrn_dmnr!dmnr!ttrhdrn
real(8),public::dmnr!ttrhdrn
integer,public::MPI_io_rank!io_rank !mpi 
integer,public::io_rank !mpi 
integer,public::MPI_num_proc_per_qcomm!num_of_mpi_per_comm !mpi 
integer,public::num_of_mpi_per_comm !mpi 
integer,public::MPI_num_qcomm!num_of_comm !mpi 
integer,public::num_of_comm !mpi 
integer,public::Num_freq_grid!nen !frequency grid 
integer,public::nen !frequency grid 
real(8),public::Lower_bound_energy_window!LOWER BOUND OF ENERGY WINDOW (eV)
real(8),public::E_LOWER!LOWER BOUND OF ENERGY WINDOW (eV)
real(8),public::Upper_bound_energy_window!UPPER BOUND OF ENERGY WINDOW (eV)
real(8),public::E_UPPER!UPPER BOUND OF ENERGY WINDOW (eV)
integer,public::flg_cRPA 
integer,public::flg_cRPA_band 
integer,public::flg_cRPA_ewin
integer,public::flg_calc_type 
integer,public::file_num_log_start 
integer,public::file_num_chi_start 
integer,public::file_num_eps_start 
integer,public::file_num_chi_base_start 
integer,public::file_num_eps_base_start
!calc_num_k 
integer,public::n_calc_q!total number of calculated q pts
integer,public,allocatable::calc_num_k(:)!calc_num_k(n_calc_q) 
namelist/param_chiqw/N_CALC_BAND,Ecut_for_eps,shift_ef,Max_excitation_energy,& 
Green_func_delt,ttrhdrn_dmna,ttrhdrn_dmnr,MPI_io_rank,MPI_num_proc_per_qcomm,&
MPI_num_qcomm,Num_freq_grid,Lower_bound_energy_window,Upper_bound_energy_window,&
flg_cRPA,flg_calc_type,file_num_log_start,file_num_chi_start,file_num_eps_start,&
file_num_chi_base_start,file_num_eps_base_start,n_calc_q 
contains
subroutine read_input(nproc)  
integer::nproc 
integer::ix 
!--
!default
Ecut_for_eps=0.0d0!1/10 of Ecut_for_psi 
Num_freq_grid=70!log grid 
N_CALC_BAND=0!NTB
SHIFT_EF=0.0d0!eV
MAX_EXCITATION_ENERGY=200.0d0!eV 
GREEN_FUNC_DELT=0.1d0!eV
TTRHDRN_DMNA=0.001d0!eV
TTRHDRN_DMNR=0.001d0!eV
MPI_IO_RANK=0!master
MPI_NUM_QCOMM=0!see below
MPI_NUM_PROC_PER_QCOMM=0!see below
LOWER_BOUND_ENERGY_WINDOW=0.0d0!not active
UPPER_BOUND_ENERGY_WINDOW=0.0d0!not active
FLG_CRPA=0!fRPA
FLG_CALC_TYPE=0!all-q-calc
FILE_NUM_LOG_START=400
FILE_NUM_CHI_START=500
FILE_NUM_EPS_START=600
FILE_NUM_CHI_BASE_START=500000
FILE_NUM_EPS_BASE_START=600000
N_CALC_Q=0
!--
!open(999,file='input.in')
!read(999,nml=param_chiqw)
read(5,nml=param_chiqw)
!--
!default mpi setting
if(MPI_num_qcomm==0)then
 MPI_num_qcomm=1
 MPI_num_proc_per_qcomm=nproc 
endif 
if(MPI_num_qcomm/=0)then
 if(mod(nproc,MPI_num_qcomm)/=0)then 
  write(6,*)'ERROR: mpi setting wrong'
  stop
 endif 
 MPI_num_proc_per_qcomm=nproc/MPI_num_qcomm 
endif 
!--
if(flg_cRPA==0)then
 write(6,*) 
 write(6,*)'========================'
 write(6,*)'THIS CALCULATION IS fRPA'
 write(6,*)'========================'
 write(6,*) 
elseif(flg_cRPA==1)then  
 write(6,*) 
 write(6,*)'========================'
 write(6,*)'THIS CALCULATION IS cRPA'
 write(6,*)'========================'
 write(6,*) 
endif 
!--
if(flg_calc_type==2.and.n_calc_q==0)then!manual mode 
 write(6,*)'ERROR: n_calc_q should not be zero'
 stop
endif 
!--
if(n_calc_q/=0)then!manual mode 
 allocate(calc_num_k(n_calc_q));calc_num_k=0
 !read(999,*)(calc_num_k(ix),ix=1,n_calc_q)
 read(5,*)(calc_num_k(ix),ix=1,n_calc_q)
else
 allocate(calc_num_k(n_calc_q));calc_num_k=0!not use
endif 
!
write(6,nml=param_chiqw) 
delta_ex=Max_excitation_energy 
delt=Green_func_delt 
dmna=ttrhdrn_dmna 
dmnr=ttrhdrn_dmnr  
io_rank=MPI_io_rank 
num_of_mpi_per_comm=MPI_num_proc_per_qcomm 
num_of_comm=MPI_num_qcomm
nen=Num_freq_grid
E_LOWER=Lower_bound_energy_window
E_UPPER=Upper_bound_energy_window
end subroutine
end module
