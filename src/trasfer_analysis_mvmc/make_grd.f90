subroutine est_ndosgrd(NTB,NTK,E_EIG,deltw,emin,emax,ndosgrd)  
  implicit none
  integer,intent(in)::NTB,NTK
  real(8),intent(in)::E_EIG(NTB,NTK)
  real(8),intent(in)::deltw 
  real(8),intent(out)::emax,emin
  integer,intent(out)::ndosgrd
  real(8)::diff 
  real(8),parameter::au=27.21151d0
  !
  emax=maxval(E_EIG)
  emin=minval(E_EIG)
  diff=emax-emin 
  !
  !define grid range
  !
  ndosgrd=int(2.0d0*diff/deltw)+1
  emax=emax+0.5d0*diff
  emin=emin-0.5d0*diff
  !
  write(6,*)
  write(6,'(a)')'GRID DATA FOR DOS'
  write(6,'(a10,f12.7)')'emin(eV)',emin*au 
  write(6,'(a10,f12.7)')'emax(eV)',emax*au 
  write(6,'(a10,f12.7)')'deltw(eV)',deltw*au 
  write(6,'(a10,i12)')'ndosgrd',ndosgrd 
  write(6,*)
  return
end subroutine

subroutine make_dosgrd(emin,deltw,ndosgrd,dosgrd)  
  implicit none
  real(8),intent(in)::emin,deltw 
  integer,intent(in)::ndosgrd
  real(8),intent(out)::dosgrd(ndosgrd) 
  integer::ie 
  !
  dosgrd=0.0d0 
  do ie=1,ndosgrd
   dosgrd(ie)=emin+dble(ie)*deltw 
  enddo 
  return
end subroutine 
