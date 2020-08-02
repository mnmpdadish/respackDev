!
!Kazuma Nakamura 20200729
!
!
module m_spinvector 
  implicit none
contains
  !
  subroutine calc_spinvector(NTG,ncomp,NWF,NTK,NG0,C0WN,a1,a2,a3)
    implicit none 
    integer,intent(in)::NTG,ncomp,NWF,NTK 
    real(8),intent(in)::a1(3),a2(3),a3(3) 
    integer,intent(in)::NG0(NTK) 
    complex(8),intent(in)::C0WN(NTG,ncomp,NWF,NTK) 
    !
    complex(8),allocatable::rho(:,:,:)!rho(ncomp,ncomp,NWF) 
    real(8),allocatable::amat(:,:)!amat(3,3) 
    real(8),allocatable::amatinv(:,:)!amatinv(3,3)
    complex(8),allocatable::svec(:)!svec(3)
    complex(8),allocatable::sveclat(:,:)!sveclat(3,NWF)
    integer::iw,ic,jc,ik,ig,iL
    complex(8)::SUM_CMPX 
    complex(8),parameter::ci=(0.0D0,1.0D0) 
    !
    if(ncomp==2)then 
     write(6,*)'.....' 
     allocate(rho(ncomp,ncomp,NWF)); rho=0.0d0 
     do iw=1,NWF      
      do ic=1,ncomp 
       do jc=1,ncomp
        SUM_CMPX=0.0d0
        do ik=1,NTK   
         do ig=1,NG0(ik)  
          SUM_CMPX=SUM_CMPX+CONJG(C0WN(ig,ic,iw,ik))*C0WN(ig,jc,iw,ik)
         enddo!ig
        enddo!ik
        rho(ic,jc,iw)=SUM_CMPX/dble(NTK)  
       enddo!jc
      enddo!ic
     enddo!iw
     !
     allocate(amat(3,3)); amat=0.0d0
     allocate(amatinv(3,3)); amatinv=0.0d0  
     allocate(svec(3)); svec=0.0d0 
     allocate(sveclat(3,NWF)); sveclat=0.0d0  
     !
     amat(:,1)=a1(:); amat(:,2)=a2(:); amat(:,3)=a3(:) 
     amatinv=amat 
     call invmat_real(3,amatinv(1,1)) 
     !
     write(6,*) 
     do iw=1,NWF 
      !       write(6,fmt='(2(1x,F20.2,SP,F20.2,"i   ",1x))')((rho(ic,jc,iw),jc=1,2),ic=1,2)
      !       write(6,'(a8,x,i5,x,6f15.8)')'iw=',iw,
      !    +   rho(2,1,iw)+rho(1,2,iw),
      !    +  (rho(2,1,iw)-rho(1,2,iw))*(-ci),
      !    +   rho(1,1,iw)-rho(2,2,iw) 
      !       write(6,*) 
      svec(1)= rho(1,2,iw)+rho(2,1,iw)
      svec(2)=(rho(1,2,iw)-rho(2,1,iw))*(-ci)
      svec(3)= rho(1,1,iw)-rho(2,2,iw) 
      !
      do iL=1,3
       SUM_CMPX=0.0d0 
       do jc=1,3
        SUM_CMPX=SUM_CMPX+amatinv(iL,jc)*svec(jc)  
       enddo!jc 
       sveclat(iL,iw)=SUM_CMPX 
      enddo!iL
      !
      !write(6,'(a8,x,i5,x,3f15.8)')'iw=',iw,(dble(sveclat(iL,iw)),iL=1,3)
      !
     enddo!iw
     !
     !write spinvector 
     !
     call wrt_spinvector(NWF,sveclat(1,1)) 
     !
    else
     write(6,'(a40)')'ncomp=1: Skip spin-vector calculation.'
    endif 
    !
    deallocate(rho,amat,amatinv,svec,sveclat) 
    !
    return
  end subroutine calc_spinvector 
  !
  subroutine wrt_spinvector(NWF,sveclat)
    implicit none 
    integer::NWF 
    complex(8)::sveclat(3,NWF)
    integer::iw,iL 
    ! 
    OPEN(194,FILE='./dir-wan/dat.spinvector') 
    write(194,'(a)')'#Spin vectror for Wannier spinor in lattice coordinates'
    write(194,'(a)')'#1:iw, 2:S_a1, 3:S_a2, 4:S_a3' 
    do iw=1,NWF
     write(194,'(i5,x,3f15.8)') iw,(dble(sveclat(iL,iw)),iL=1,3)
    enddo 
    !
    return 
  end subroutine wrt_spinvector  
  !
  subroutine invmat_real(nm,mat)
    implicit none 
    integer,intent(in)::nm
    real(8),intent(inout)::mat(nm,nm)
    integer::ipiv(nm)
    integer::Lwork 
    real(8),allocatable::work(:)
    integer::info 
    Lwork=10*nm
    allocate(work(Lwork))
    info=0
    call dgetrf(nm,nm,mat,nm,ipiv,info)
    call dgetri(nm,mat,nm,ipiv,work,Lwork,info)
    if(info/=0) then
     write(6,*) 'info (subrouitine inv):',info
     stop
    endif 
    deallocate(work)
    return 
  end subroutine invmat_real
  !
end module m_spinvector 
