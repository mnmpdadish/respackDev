program main 
  implicit none 
  integer,parameter::NWF=5!TOTAL NUMBER OF OCCUPIED BAND
  real(8),parameter::s=1.0d0 
  real(8),parameter::t=0.0d0 
  integer::i,j,k,l,is,js,ks,ls,dum1,dum2,dum3,dum4,m     
  real(8)::dataT0(NWF,NWF),dataT(NWF,NWF),dataW4(NWF,NWF,NWF,NWF) 
  !--
  real(8)::eps 
  real(8),ALLOCATABLE::Eout(:),eig(:) 
  complex(8),ALLOCATABLE::Hin(:,:),coef(:,:) 
  !-- 
  real(8)::occ(NWF),dlt(NWF,NWF)
  complex(8)::DM(NWF,NWF),SUM_CMPX  
  !--
  do i=1,NWF
   read(101,*)(dataT0(i,j),j=1,NWF)
  enddo 
  do i=1,NWF
   do j=1,NWF
    do k=1,NWF
     do l=1,NWF
      read(102,*) dum1,dum2,dum3,dum4,dataW4(i,j,k,l) 
     enddo 
    enddo 
   enddo 
  enddo 
  dataW4=dataW4*s 
  !---
  !call two_center_only(NWF,dataW4(1,1,1,1)) 
  !---
  allocate(Hin(NWF,NWF));Hin(:,:)=0.0D0               
  allocate(coef(NWF,NWF));coef(:,:)=0.0D0               
  allocate(Eout(NWF));Eout(:)=0.0d0
  allocate(eig(NWF));eig(:)=0.0d0
  Hin=dataT0  
  call diagV(NWF,Hin(1,1),Eout(1)) 
  coef=Hin 
  eig=Eout 
  write(6,'(a20)')'EIGENVALUS of dataT0:' 
  write(6,*) 
  do i=1,NWF
   write(6,*) Eout(i)
  enddo 
  write(6,*) 
  write(6,'(a20,f15.10)')'shift value eV:',t 
  write(6,*) 
  do i=1,NWF
  if(i>3) eig(i)=eig(i)+t 
  write(6,*) eig(i)
  enddo 
  deallocate(Hin,Eout) 
  !---
  dataT=0.0d0 
  do i=1,NWF
   do j=1,NWF
    SUM_CMPX=0.0d0 
    do m=1,NWF
     SUM_CMPX=SUM_CMPX+coef(i,m)*eig(m)*conjg(coef(j,m))
    enddo 
    dataT(i,j)=dble(SUM_CMPX) 
   enddo 
  enddo 
  write(6,*) 
  write(6,'(a20)')'dataT0:'
  write(6,*) 
  do i=1,NWF
   write(6,'(100f15.5)')(dataT0(i,j),j=1,NWF)
  enddo 
  write(6,*) 
  write(6,'(a20)')'shifted dataT:'
  write(6,*) 
  do i=1,NWF
   write(6,'(100f15.5)')(dataT(i,j),j=1,NWF)
  enddo 
  !---
  !occ for type-A config
  !occ(1)=2.0d0
  !occ(2)=0.5d0
  !occ(3)=0.5d0
  !occ(4)=0.0d0
  !occ(5)=0.0d0
  !---
  !occ for type-B config
  !occ(1)=1.5d0
  !occ(2)=1.5d0
  !occ(3)=0.0d0
  !occ(4)=0.0d0
  !occ(5)=0.0d0
  !---
  !occ xtapp 
  occ(1)=1.0d0
  occ(2)=1.0d0
  occ(3)=1.0d0
  occ(4)=0.0d0
  occ(5)=0.0d0
  !---
  do i=1,NWF
   do j=1,NWF
    SUM_CMPX=0.0d0
    do m=1,NWF
     SUM_CMPX=SUM_CMPX+occ(m)*coef(i,m)*CONJG(coef(j,m))
    enddo
    DM(i,j)=SUM_CMPX 
   enddo
  enddo
  write(6,*) 
  write(6,'(a20)')'DENSITY MATRIX:'
  write(6,*) 
  SUM_CMPX=0.0d0
  do i=1,NWF
   SUM_CMPX=SUM_CMPX+DM(i,i) 
   write(6,'(10f15.10)')(dble(DM(i,j)),j=1,NWF)
  enddo 
  write(6,*) 
  write(6,'(a20,2f15.10)')'Ntot:',SUM_CMPX 
  !---
  !No 
  !dlt=0.0d0 
  !---
  !HA(diag) 
  !dlt=0.0d0 
  !do i=1,NWF
  ! do j=1,NWF
  !  do k=1,NWF
  !   dlt(i,j)=dlt(i,j)+dble(DM(k,k))*dataW4(i,j,k,k)
  !  enddo 
  ! enddo 
  !enddo 
  !---
  !HA
  !dlt=0.0d0 
  !do i=1,NWF
  ! do j=1,NWF
  !  do k=1,NWF
  !   do l=1,NWF
  !    dlt(i,j)=dlt(i,j)+dble(DM(k,l))*dataW4(i,j,l,k)
  !   enddo 
  !  enddo 
  ! enddo 
  !enddo 
  !---
  !HF(diag)
  dlt=0.0d0 
  do i=1,NWF
   do j=1,NWF
    do k=1,NWF
     dlt(i,j)=dlt(i,j)+dble(DM(k,k))*(dataW4(i,j,k,k)-0.5d0*dataW4(i,k,k,j)) 
    enddo 
   enddo 
  enddo 
  !---
  !HF
  !dlt=0.0d0 
  !do i=1,NWF
  ! do j=1,NWF
  !  do k=1,NWF
  !   do l=1,NWF
  !    dlt(i,j)=dlt(i,j)+dble(DM(k,l))*(dataW4(i,j,l,k)-0.5d0*dataW4(i,k,l,j)) 
  !   enddo 
  !  enddo 
  ! enddo 
  !enddo 
  !---
  write(6,*)  
  write(6,'(a20)')'ONE-BODY CORRECTION:'
  write(6,*) 
  do i=1,NWF
   write(6,'(10f15.10)')(dlt(i,j),j=1,NWF)
  enddo 
  !transfer  
  write(6,*)  
  write(6,*)'#transfer'
  k=0
  do i=1,NWF
   do j=1,NWF
    do is=1,2
     k=k+1
     write(6,'(i0,3i3,2f10.5)') i-1,is-1,j-1,is-1,-dataT(i,j)+dlt(i,j),0.0d0 
    enddo 
   enddo 
  enddo 
  write(6,*) k
  !W4 
  write(6,*)'#NinterAll'
  m=0
  do i=1,NWF
   do j=1,NWF 
    do k=1,NWF
     do l=1,NWF 
      do is=1,2
       do ks=1,2
        m=m+1
        write(6,'(i0,7i3,2f10.5)') i-1,is-1,j-1,is-1,k-1,ks-1,l-1,ks-1,0.5d0*dataW4(i,j,k,l),0.0d0
       enddo
      enddo
     enddo
    enddo
   enddo
  enddo
  write(6,*) m 
  !---
  stop
end program 

subroutine two_center_only(NWF,dataW4)
  implicit none 
  integer::NWF
  real(8),intent(inout)::dataW4(NWF,NWF,NWF,NWF) 
  real(8),allocatable::dataW4_tmp(:,:,:,:)!(NWF,NWF,NWF,NWF) 
  integer::i,j,k,l 
  ! 
  allocate(dataW4_tmp(NWF,NWF,NWF,NWF)); dataW4_tmp=0.0d0 
  do i=1,NWF
   do j=1,NWF
    do k=1,NWF
     do l=1,NWF 
      if((i==j).and.(k==l))then
       dataW4_tmp(i,j,k,l)=dataW4(i,j,k,l) 
      endif 
      if((i==k).and.(j==l))then
       dataW4_tmp(i,j,k,l)=dataW4(i,j,k,l) 
      endif 
      if((i==l).and.(j==k))then
       dataW4_tmp(i,j,k,l)=dataW4(i,j,k,l) 
      endif 
     enddo 
    enddo 
   enddo 
  enddo 
  dataW4=dataW4_tmp 
  deallocate(dataW4_tmp) 
  ! 
  return
end subroutine 
