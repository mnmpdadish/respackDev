module m_hist 
  implicit none
contains
  !
  subroutine calculate_hist(Ndim_TR,delw,TR) 
    implicit none 
    integer,intent(in)::Ndim_TR 
    real(8),intent(in)::delw!Grid width in histgram 
    complex(8),intent(inout)::TR(Ndim_TR,Ndim_TR) 
    !
    real(8),allocatable::EIG_TR(:)!EIG_TR(Ndim_TR) 
    real(8),allocatable::hist(:)!hist(ndosgrd) 
    real(8)::emin,emax
    integer::ndosgrd 
    integer::emin_grd,emax_grd 
    !
    !diagonalize TR 
    !
    allocate(EIG_TR(Ndim_TR)); EIG_TR=0.0d0 
    call diagN(Ndim_TR,TR(1,1),EIG_TR(1)) 
    !
    !make dos-grid
    !
    call est_ndosgrd(Ndim_TR,1,EIG_TR(1),delw,emin,emax,ndosgrd)  
    !
    !calc histgram
    !
    emin_grd=nint(emin/delw) 
    emax_grd=nint(emax/delw)
    allocate(hist(emin_grd:emax_grd)); hist=0.0d0 
    call calc_hist(Ndim_TR,emin_grd,emax_grd,delw,EIG_TR(1),hist(emin_grd)) 
    !
    !wrt histgram
    ! 
    call wrt_hist(emin_grd,emax_grd,delw,hist(emin_grd)) 
    deallocate(EIG_TR,hist) 
    !
    return
  end subroutine 
  !
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
  !
  subroutine diagN(nm,mat,eig)
    implicit none 
    integer,intent(in)::nm
    complex(8),intent(inout)::mat(nm,nm)
    real(8),intent(out)::eig(nm)
    integer::LWORK,LRWORK,LIWORK  
    integer,allocatable::iwork_zheevd(:)
    real(8),allocatable::rwork_zheevd(:)
    complex(8),allocatable::work_zheevd(:)
    integer::ind
    real(8)::eps 
    !
    LWORK= 2*nm+nm**2
    LRWORK=1+12*nm+3*nm**2
    LIWORK=3+10*nm 
    allocate(work_zheevd(LWORK));work_zheevd(:)=0.0d0
    allocate(rwork_zheevd(LRWORK));rwork_zheevd(:)=0.0d0
    allocate(iwork_zheevd(LIWORK));iwork_zheevd(:)=0
    eps=1.0d-18
    ind=0                 
    !
    call zheevd("N","U",nm,mat,nm,eig,work_zheevd,LWORK,rwork_zheevd,LRWORK,iwork_zheevd,LIWORK,ind)
    !
    if(ind/=0)then 
     write(6,*)'ind=',ind 
     stop
    endif 
    !
    deallocate(work_zheevd,rwork_zheevd,iwork_zheevd) 
    return 
  end subroutine
  !
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
  !
  subroutine wrt_hist(emin_grd,emax_grd,delw,hist) 
    implicit none 
    integer,intent(in)::emin_grd,emax_grd
    real(8),intent(in)::delw 
    real(8),intent(in)::hist(emin_grd:emax_grd) 
    integer::ix 
    real(8),parameter::au=27.21151d0 
    !
    !OPEN(302,W,file='./dat.hist')
    !
    OPEN(302,file='./dat.hist') 
    rewind(302) 
    do ix=emin_grd,emax_grd 
     write(302,'(2f15.10)') dble(ix)*delw*au,hist(ix)/au  
    enddo 
    close(302) 
    !
    return 
  end subroutine 
  !
end module m_hist 
