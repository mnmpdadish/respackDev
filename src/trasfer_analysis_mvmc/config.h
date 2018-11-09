!
!config.h 
!
implicit none 
real(8),parameter::au=27.21151d0
real(8),parameter::pi=dacos(-1.0d0)
real(8),parameter::tpi=2.0d0*pi 
complex(8),parameter::ci=(0.0D0,1.0D0) 
complex(8),parameter::tci=(0.0D0,2.0D0) 
!
!command line
!
integer::iargc 
integer::ncount  
character(len=8)::arg
real(8),allocatable::real_arg(:)!real_arg(ncount)  
!
!&param_transfer_analysis 
!
logical::zvo  !transfer analysis for ZVO data 
logical::ztr  !transfer analysis for ZTRans.def 
real(8)::delt !Greens function delt in eV
real(8)::dmna !Ttrhdrn parameter dmna in eV
real(8)::dmnr !Ttrhdrn parameter dmnr in eV
real(8)::delw !Grid width in eV
real(8)::flwe !Flg whether calculate weighted transfers  
real(8)::thtr !Threshold for transfer integral in eV
real(8)::elnm !Total number of electrons in unitcell
!
integer::flg_weight         !Flg whether calculate weighted transfers (0:not calc, 1:calc)
real(8)::threshold_transfer !Threshold for tranfer integral (eV)
real(8)::electron_number    !Total number of electron in unitcell 
!
!eigenvalue and eigenstates
!
integer::nkb1,nkb2,nkb3 
real(8),allocatable::EKS(:,:) !EKS(NWF,NTK) 
complex(8),allocatable::VKS(:,:,:) !VKS(NWF,NWF,NTK)   
!
!dos 
!
real(8)::emax !=maxval(E_EIG)
real(8)::emin !=minval(E_EIG)
real(8)::FermiEnergy 
integer::ndosgrd !=int(2.0d0*diff/dlt)+1
real(8),allocatable::dosgrd(:) !dosgrd(ndosgrd) 
real(8),allocatable::dos(:) !dos(ndosgrd) 
real(8),allocatable::efline(:) !efline(ndosgrd)   
!
!interpolated band disp
!
real(8),allocatable::kdata(:) !kdata(NSK_BAND_DISP) 
!
!local
!
integer::i !ie,je,ik,ikir,ib,i    
!real(8)::SUM_REAL
!
!end config.h 
!
