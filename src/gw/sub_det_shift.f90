subroutine det_shift(NTK,NWF,nsgm,FermiEnergy,sgmw,EKS,EMK,shift_value) 
  implicit none
  integer::NTK,NWF,nsgm
  real(8)::FermiEnergy 
  real(8)::sgmw(nsgm)  
  complex(8)::EMK(NWF,NTK,nsgm)           
  real(8)::EKS(NWF,NTK)           
  integer::min_ib,min_ik,min_ie 
  real(8)::min_diff,min_eak,diff_ie,diff_ef 
  integer::ik,ib,ie
  real(8)::shift_value 
  real(8),parameter::au=27.21151d0
  !
  min_diff=1.0d0
  do ik=1,NTK 
   do ib=1,NWF 
    diff_ef=dabs(EKS(ib,ik)-FermiEnergy)
    if(diff_ef<min_diff)then 
     min_diff=diff_ef  
     min_eak=EKS(ib,ik) 
     min_ib=ib
     min_ik=ik
    endif 
   enddo!ib
  enddo!ik 
  !
  min_diff=0.1d0
  do ie=1,nsgm 
   diff_ie=abs(sgmw(ie)-min_eak)
   if(diff_ie<min_diff)then 
    min_diff=diff_ie 
    min_ie=ie
   endif 
   !write(6,'(a10,f15.10,a10,f15.10)')'sgmw(ie)',sgmw(ie),'min_diff',min_diff 
  enddo!ie 
  shift_value=dble(EMK(min_ib,min_ik,min_ie))-EKS(min_ib,min_ik) 
  !
  write(6,*)'min_ib= ',min_ib
  write(6,*)'min_ik= ',min_ik
  write(6,*)'min_ie= ',min_ie 
  write(6,*)'closest EKS(eV)=',EKS(min_ib,min_ik)*au 
  write(6,*)'closest sgmw(eV)=',sgmw(min_ie)*au 
  write(6,*)'FermiEnergy(eV)=',FermiEnergy*au 
  write(6,*)'shift(eV)=',shift_value*au 
  ! 
return
end
