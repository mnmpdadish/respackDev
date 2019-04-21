!---
!config.h 
!---
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
integer::i 
!
!transfer_analysis 
!
logical::dos !DOS calc for zvo data 
logical::bnd !BaND dispersion calc for zvo data 
logical::frm !FeRMi surface calc for fine k-mesh 
logical::his !HIStgram analysis for ztrans.def 
real(8)::delt !Greens function delt in eV
real(8)::elnm !Total number of electrons in unitcell
real(8)::ecut !Energy cutoff for transfer integral in eV
real(8)::rcut !Distance cutoff for transfer integral in AA 
real(8)::diff !Match threshold for two transfer integral in eV 
real(8)::electron_number !Total number of electron in unitcell 
real(8)::threshold_e !Energy cutoff for tranfer integral in eV
real(8)::threshold_r !Discance cutoff for tranfer integral in AA
real(8)::diff_transfers !Match threshold for two transfer integral in eV 
real(8)::delw !Grid width in eV
integer::kgd(3) !k grid  
integer::flg_weight !Flg whether calculate weighted transfers (0:not calc, 1:calc)
!
!eigenvalue and eigenstates
!
integer::ncalck,La1,La2,La3 
real(8),allocatable::kvec(:,:)!kvec(3,ncalck) 
real(8),allocatable::EKS(:,:)!EKS(NWF,ncalck) 
complex(8),allocatable::VKS(:,:,:)!VKS(NWF,NWF,ncalck)   
complex(8),allocatable::FR(:,:,:,:,:)!FR(NWF,NWF,-La1:La1,-La2:La2,-La3:La3) 
!
!--
!end config.h 
!--
