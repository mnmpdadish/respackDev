module m_rd_transdef 
  implicit none
  private 
  public::rd_transdef 
  !trans.def(400) 
  integer,public::Ndim_TR 
  complex(8),public,allocatable::TR(:,:) 
  contains
  !--
  subroutine rd_transdef 
    integer,allocatable::isite(:)
    integer,allocatable::ispin(:)
    integer,allocatable::jsite(:)
    integer,allocatable::jspin(:)
    complex(8),allocatable::TR1dim(:) 
    character(len=80)::dum_ch
    integer::NTransfer 
    integer::Nlat_times_Norb 
    integer::i,n,m
    real(8),parameter::au=27.21151d0
    !
    !OPEN(400,R,FILE='./dir-model/trans.def') 
    !
    OPEN(400,FILE='./dir-model/trans.def') 
    rewind(400) 
    read(400,*) dum_ch
    read(400,*) dum_ch,NTransfer 
    read(400,*) dum_ch
    read(400,*) dum_ch
    read(400,*) dum_ch
    allocate(isite(NTransfer));isite=0
    allocate(ispin(NTransfer));ispin=0
    allocate(jsite(NTransfer));jsite=0
    allocate(jspin(NTransfer));ispin=0
    allocate(TR1dim(NTransfer));TR1dim=0.0d0
    do i=1,NTransfer
     read(400,*) isite(i),ispin(i),jsite(i),jspin(i),TR1dim(i) 
    enddo 
    close(400) 
    !
    Nlat_times_Norb=maxval(isite)+1
    !
    if(Nlat_times_Norb/=(maxval(jsite)+1))then
     write(6,*) 'error; stop'
     stop
    endif 
    !
    Ndim_TR=Nlat_times_Norb*2 
    !
    write(6,'(i10)') NTransfer  
    write(6,'(i10)') Nlat_times_Norb 
    write(6,'(i10)') Ndim_TR 
    !
    allocate(TR(Ndim_TR,Ndim_TR));TR=0.0d0 
    do i=1,NTransfer
     n=isite(i)+ispin(i)*Nlat_times_Norb+1
     m=jsite(i)+jspin(i)*Nlat_times_Norb+1
     TR(n,m)=TR1dim(i) 
    enddo  
    TR=-TR/au !au<-eV
    ! 
    deallocate(isite,ispin,jsite,jspin,TR1dim) 
  end subroutine
  !--
end module 
