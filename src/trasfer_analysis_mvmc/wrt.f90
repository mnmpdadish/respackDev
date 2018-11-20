subroutine wrt_dos(ndosgrd,dosgrd,dos,efline) 
  implicit none
  integer,intent(in)::ndosgrd
  real(8),intent(in)::dosgrd(ndosgrd)
  real(8),intent(in)::dos(ndosgrd)
  real(8),intent(in)::efline(ndosgrd)
  integer::ie 
  real(8),parameter::au=27.21151d0 
  !
  !OPEN(300,W,file='./dir-tr/dat.dos')
  !OPEN(301,W,file='./dir-tr/dat.efline')
  !
  OPEN(300,file='./dir-tr/dat.dos') 
  OPEN(301,file='./dir-tr/dat.efline') 
  rewind(300);rewind(301) 
  do ie=1,ndosgrd
   write(300,'(2f15.10)') dosgrd(ie)*au,dos(ie)/au 
   write(301,'(2f15.10)') dosgrd(ie)*au,efline(ie)  
  enddo!ie 
  close(300)
  close(301) 
  return
end subroutine 

subroutine wrt_iband(NWF,NSK_BAND_DISP,kdata,EKS) 
  implicit none
  integer,intent(in)::NWF,NSK_BAND_DISP
  real(8),intent(in)::kdata(NSK_BAND_DISP)
  real(8),intent(in)::EKS(NWF,NSK_BAND_DISP) 
  integer::ib,ik 
  LOGICAL::REVERSE 
  real(8),parameter::au=27.21151d0 
  !
  !OPEN(114,W,FILE='./dir-tr/dat.iband') 
  !
  OPEN(114,FILE='./dir-tr/dat.iband') 
  write(114,'(a)')'#Wannier interpolaed band'
  write(114,'(a)')'#1:k, 2:Energy [eV]' 
  REVERSE=.TRUE.        
  do ib=1,NWF 
   if(REVERSE)then 
    do ik=1,NSK_BAND_DISP                     
     write(114,'(2f20.10)') kdata(ik),EKS(ib,ik)*au
    enddo!ik        
    REVERSE=.FALSE.        
   else         
    do ik=NSK_BAND_DISP,1,-1          
     write(114,'(2f20.10)') kdata(ik),EKS(ib,ik)*au
    enddo!ik        
    REVERSE=.TRUE.        
   endif!REVERSE                   
  enddo!ib 
  close(114) 
  !
  return 
end subroutine   
  
subroutine wrt_hist(emin_grd,emax_grd,delw,hist) 
  implicit none 
  integer,intent(in)::emin_grd,emax_grd
  real(8),intent(in)::delw 
  real(8),intent(in)::hist(emin_grd:emax_grd) 
  integer::ix 
  real(8),parameter::au=27.21151d0 
  !
  !OPEN(302,W,file='./dir-tr/dat.hist')
  !
  OPEN(302,file='./dir-tr/dat.hist') 
  rewind(302) 
  do ix=emin_grd,emax_grd 
   write(302,'(2f15.10)') dble(ix)*delw*au,hist(ix)/au  
  enddo 
  close(302) 
  !
  return 
end subroutine 
