subroutine est_ef(ndosgrd,deltw,electron_number,dosgrd,dos,FermiEnergy)
  implicit none
  integer,intent(in)::ndosgrd
  real(8),intent(in)::deltw,electron_number 
  real(8),intent(in)::dosgrd(ndosgrd) 
  real(8),intent(in)::dos(ndosgrd) 
  real(8),intent(out)::FermiEnergy 
  real(8)::SUM_REAL
  integer::ie 
  real(8),parameter::au=27.21151d0 
  SUM_REAL=0.0d0 
  do ie=1,ndosgrd
   SUM_REAL=SUM_REAL+deltw*dos(ie) 
  enddo 
  write(6,*)'SUM of DOS',SUM_REAL
  SUM_REAL=0.0d0 
  do ie=1,ndosgrd
   if(SUM_REAL>=electron_number)goto 3000 
   SUM_REAL=SUM_REAL+deltw*dos(ie) 
  enddo 
3000 FermiEnergy=dosgrd(ie) 
  write(6,*)'electron_number',SUM_REAL
  write(6,*)'FermiEnergy=',dosgrd(ie)*au  
  return
end subroutine 

subroutine make_efline(ndosgrd,FermiEnergy,deltw,dosgrd,dos,efline) 
  implicit none 
  integer,intent(in)::ndosgrd
  real(8),intent(in)::FermiEnergy,deltw
  real(8),intent(in)::dosgrd(ndosgrd) 
  real(8),intent(in)::dos(ndosgrd) 
  real(8),intent(out)::efline(ndosgrd) 
  real(8)::diff,diff_old,dosmax,N0ttr 
  integer::ie,ief 
  real(8)::SUM_REAL 
  real(8),parameter::au=27.21151d0
  !
  !efline(ndosgrd) 
  !
  diff_old=1.0d0 
  do ie=1,ndosgrd
   diff=dabs(dosgrd(ie)-FermiEnergy) 
   if(diff_old>diff) then 
    ief=ie
    diff_old=diff 
   endif  
  enddo!ie  
  dosmax=maxval(abs(dos(:)))/au  
  dosmax=dble(nint(dosmax*10.0d0))/10.0d0 
  do ie=1,ndosgrd 
   if(ie<=ief) then 
    efline(ie)=dosmax
   else 
    efline(ie)=0.0d0 
   endif 
  enddo!ie 
  !
  !Integral dos(w) dw
  !
  SUM_REAL=0.0d0 
  do ie=1,ndosgrd 
   SUM_REAL=SUM_REAL+dos(ie)*deltw 
   !write(6,'(2f15.10)') dosgrd(ie)*au,dos(ie)/au  
  enddo 
  write(6,*)
  write(6,'(a40,f15.7)')'Integral dos(w) dw',SUM_REAL 
  N0ttr=dos(ief)/2.0d0 !2 is spin
  write(6,'(a40,f15.7)')'N(0) in au per calc cell',N0ttr 
  write(6,'(a40,f15.7)')'N(0) in eV per calc cell',N0ttr/au   
  write(6,*)
  return
end subroutine 
