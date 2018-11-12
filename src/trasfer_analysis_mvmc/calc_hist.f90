subroutine calc_hist(N_eig,N_ini,N_end,del,EIG,hist) 
  implicit none
  integer,intent(in)::N_eig,N_ini,N_end 
  real(8),intent(in)::del 
  real(8),intent(in)::EIG(N_eig) 
  real(8),intent(out)::hist(N_ini:N_end) 
  integer::ie,ix 
  real(8),parameter::au=27.21151d0
  !
  hist=0.0d0 
  do ie=1,N_eig 
   ix=nint(EIG(ie)/del) 
   hist(ix)=hist(ix)+1.0d0 
  enddo
  !
  return 
end subroutine 
