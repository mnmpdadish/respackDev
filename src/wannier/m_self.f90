module m_self 
  implicit none
contains
  !
  subroutine calculate_self(ncomp,NTB,NTK,nkb1,nkb2,nkb3,NK_irr,numirr,numMK,FermiEnergy,filename,b1,b2,b3,SK0,E_EIGI)
    implicit none 
    integer,intent(in)::ncomp,NTB,NTK,nkb1,nkb2,nkb3,NK_irr
    integer,intent(in)::numirr(NTK)
    integer,intent(in)::numMK(NK_irr)
    real(8),intent(inout)::FermiEnergy
    character(99),intent(in)::filename 
    real(8),intent(in)::b1(3),b2(3),b3(3)
    real(8),intent(in)::SK0(3,NTK)
    real(8),intent(in)::E_EIGI(NTB,NK_irr) 
    !
    real(8)::emax!=maxval(E_EIGI)
    real(8)::emin!=minval(E_EIGI)
    integer::ndosgrd!=int(2.0d0*diff/dlt)+1
    real(8),allocatable::dosgrd(:)!dosgrd(ndosgrd) 
    real(8),allocatable::dos_i(:)!dos_i(ndosgrd) 
    real(8),allocatable::dos_r(:)!dos_r(ndosgrd) 
    integer::flg_causal 
    !
    real(8),parameter::au=27.21151d0
    !real(8),parameter::delt=0.00001d0/au!Greens function delt in au 
    !real(8),parameter::delt=0.0001d0/au!Greens function delt in au 
    !real(8),parameter::delt=0.001d0/au!Greens function delt in au 
    !real(8),parameter::delt=0.005d0/au!Greens function delt in au 
    real(8),parameter::delt=0.01d0/au!Greens function delt in au 
    !real(8),parameter::delt=0.02d0/au!Greens function delt in au 
    !real(8),parameter::delt=0.05d0/au!Greens function delt in au 
    !real(8),parameter::delt=0.1d0/au!Greens function delt in au 
    !
    real(8),parameter::dmna=1.0d-3!Ttrhdrn parameter dmna in au 
    real(8),parameter::dmnr=1.0d-3!Ttrhdrn parameter dmnr in au 
    !
    !real(8),parameter::delw=0.001d0/au!2.0d0*delt!Grid width in au 
    !real(8),parameter::delw=0.005d0/au!2.0d0*delt!Grid width in au 
    !real(8),parameter::delw=0.01d0/au!2.0d0*delt!Grid width in au 
    !real(8),parameter::delw=0.05d0/au!2.0d0*delt!Grid width in au 
    real(8),parameter::delw=0.1d0/au!2.0d0*delt!Grid width in au 
    !real(8),parameter::delw=2.0d0*delt!Grid width in au 
    ! 
    !grid
    !
    call est_ndosgrd(NTB,NK_irr,E_EIGI(1,1),delw,emin,emax,ndosgrd)  
    allocate(dosgrd(ndosgrd));dosgrd=0.0d0 
    call make_dosgrd(emin,delw,ndosgrd,dosgrd(1))  
    !
    !calc self 
    !
    allocate(dos_r(ndosgrd));dos_r=0.0d0 
    allocate(dos_i(ndosgrd));dos_i=0.0d0 
    call calc_dos_self(ncomp,NTB,NTK,nkb1,nkb2,nkb3,NK_irr,numirr(1),numMK(1),FermiEnergy,ndosgrd,dosgrd(1),E_EIGI(1,1),SK0(1,1),delt,dmnr,dmna,b1(1),b2(1),b3(1),dos_r(1),dos_i(1));flg_causal=1 
    !
    !wrt dat.self.imag 
    !
    call wrt_dos_imag(filename,ndosgrd,dosgrd(1),dos_i(1),FermiEnergy) 
    !
    !wrt dat.dos.real  
    !
    if(flg_causal==1)then
     call wrt_dos_real(filename,ndosgrd,dosgrd(1),dos_r(1),FermiEnergy) 
    endif 
    !
    return
  end subroutine calculate_self 
  ! 
  subroutine est_ndosgrd(NTB,NK_irr,E_EIGI,deltw,emin,emax,ndosgrd)  
    implicit none
    integer,intent(in)::NTB,NK_irr
    real(8),intent(in)::E_EIGI(NTB,NK_irr)
    real(8),intent(in)::deltw 
    real(8),intent(out)::emax,emin
    integer,intent(out)::ndosgrd
    real(8)::diff 
    real(8),parameter::au=27.21151d0
    !
    emax=maxval(E_EIGI)
    emin=minval(E_EIGI)
    diff=emax-emin 
    !
    !define grid range
    !
    ndosgrd=int(2.0d0*diff/deltw)+1
    emax=emax+0.50d0*diff
    emin=emin-0.50d0*diff
    !
    write(6,*)
    write(6,'(a40)')'+++ m_dos: est_ndosgrd +++'
    write(6,'(a40)')'GRID DATA FOR DOS'
    write(6,'(a40,f15.8)')'emin(eV)',emin*au 
    write(6,'(a40,f15.8)')'emax(eV)',emax*au 
    write(6,'(a40,f15.8)')'deltw(eV)',deltw*au 
    write(6,'(a40,i15)')'ndosgrd',ndosgrd 
    write(6,*)
    return
  end subroutine
  !
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
  !
  subroutine make_index_kpt(NTK,nkb1,nkb2,nkb3,SK0,index_kpt)      
    implicit none 
    integer::NTK,nkb1,nkb2,nkb3
    real(8)::SK0(3,NTK)  
    integer::ik,ix,iy,iz
    real(8)::x,y,z
    integer::index_kpt(nkb1,nkb2,nkb3)    
    ! 
    !if(MOD(NTK,2)/=0)then 
    ! !write(6,*)'i am in make_index for odd'
    ! do ik=1,NTK 
    !  x=SK0(1,ik)*dble(nkb1) 
    !  y=SK0(2,ik)*dble(nkb2)
    !  z=SK0(3,ik)*dble(nkb3)  
    !  x=x+(dble(nkb1)-1.0d0)/2.0d0 
    !  y=y+(dble(nkb2)-1.0d0)/2.0d0
    !  z=z+(dble(nkb3)-1.0d0)/2.0d0 
    !  ix=idnint(x)+1
    !  iy=idnint(y)+1
    !  iz=idnint(z)+1
    !  index_kpt(ix,iy,iz)=ik
    ! enddo 
    !else!20170316 
    ! !write(6,*)'i am in make_index for even'
    ! do ik=1,NTK 
    !  x=SK0(1,ik)*dble(nkb1) 
    !  y=SK0(2,ik)*dble(nkb2)
    !  z=SK0(3,ik)*dble(nkb3)  
    !  x=x+dble(nkb1)/2.0d0 
    !  y=y+dble(nkb2)/2.0d0
    !  z=z+dble(nkb3)/2.0d0 
    !  ix=idnint(x)
    !  iy=idnint(y)
    !  iz=idnint(z)
    !  index_kpt(ix,iy,iz)=ik
    ! enddo 
    !endif 
    !--
    !
    !20190520 Kazuma Nakamura
    !
    do ik=1,NTK 
     x=SK0(1,ik)*dble(nkb1) 
     y=SK0(2,ik)*dble(nkb2)
     z=SK0(3,ik)*dble(nkb3)  
     x=x+(dble(nkb1)-dble(mod(nkb1,2)))/2.0d0 
     y=y+(dble(nkb2)-dble(mod(nkb2,2)))/2.0d0 
     z=z+(dble(nkb3)-dble(mod(nkb3,2)))/2.0d0 
     ix=idnint(x)+mod(nkb1,2)
     iy=idnint(y)+mod(nkb2,2)
     iz=idnint(z)+mod(nkb3,2)
     index_kpt(ix,iy,iz)=ik
    enddo 
    ! 
    return 
  end subroutine  
  !
  subroutine wrt_dos_imag(filename,ndosgrd,dosgrd,dos,FermiEnergy) 
    implicit none
    integer,intent(in)::ndosgrd
    real(8),intent(in)::dosgrd(ndosgrd)
    real(8),intent(in)::dos(ndosgrd)
    real(8),intent(in)::FermiEnergy 
    character(99),intent(in)::filename
    character(99)::fname,fname_arg 
    integer::ie 
    real(8),parameter::au=27.21151d0 
    !
    !OPEN(300,W,file='./dir-wan/dat.dos.xxx-total')
    !
    fname_arg='-total'
    !
    write(fname,'(a,a)') trim(filename),trim(fname_arg) 
    OPEN(300,file=trim(fname)) 
    rewind(300) 
    write(300,'(a,f15.10)')'# FermiEnergy (eV):',FermiEnergy*au 
    do ie=1,ndosgrd
     !write(300,'(2f20.10)') dosgrd(ie)*au,dos(ie)/au 
     write(300,'(2f20.10)')(dosgrd(ie)-FermiEnergy)*au,dos(ie)/au 
    enddo!ie 
    close(300)
    return
  end subroutine wrt_dos_imag 
  !
  subroutine wrt_dos_real(filename,ndosgrd,dosgrd,dos,FermiEnergy) 
    implicit none
    integer,intent(in)::ndosgrd
    real(8),intent(in)::dosgrd(ndosgrd)
    real(8),intent(in)::dos(ndosgrd)
    real(8),intent(in)::FermiEnergy 
    character(99),intent(in)::filename
    character(99)::fname,fname_arg 
    integer::ie 
    real(8),parameter::au=27.21151d0 
    !
    !OPEN(500,W,file='./dir-wan/dat.dos.xxx-total-real')
    !
    fname_arg='-total-real'
    !
    write(fname,'(a,a)') trim(filename),trim(fname_arg) 
    OPEN(500,file=trim(fname)) 
    rewind(500) 
    write(500,'(a,f15.10)')'# FermiEnergy (eV):',FermiEnergy*au 
    do ie=1,ndosgrd
     !write(500,'(2f20.10)') dosgrd(ie)*au,dos(ie)/au 
     write(500,'(3f20.10)')(dosgrd(ie)-FermiEnergy)*au,dos(ie)/au 
    enddo!ie 
    close(500)
    return
  end subroutine wrt_dos_real
  !
  subroutine calc_dos_self(ncomp,NTB,NTK,nkb1,nkb2,nkb3,NK_irr,numirr,numMK,FermiEnergy,ndosgrd,dosgrd,E_EIGI,SK0,delt,dmnr,dmna,b1,b2,b3,self_r,self_i) 
    use m_tetrahedron
    implicit none
    integer,intent(in)::ncomp,NTB,NTK,nkb1,nkb2,nkb3,NK_irr,ndosgrd
    integer,intent(in)::numirr(NTK) 
    integer,intent(in)::numMK(NK_irr) 
    real(8),intent(in)::FermiEnergy 
    real(8),intent(in)::dosgrd(ndosgrd) 
    real(8),intent(in)::E_EIGI(NTB,NK_irr)           
    real(8),intent(in)::SK0(3,NTK)           
    real(8),intent(in)::delt,dmnr,dmna 
    real(8),intent(in)::b1(3),b2(3),b3(3)
    real(8),intent(out)::self_r(ndosgrd) 
    real(8),intent(out)::self_i(ndosgrd) 
    !
    integer::index_kpt(nkb1,nkb2,nkb3) 
    integer::imt1(4*nkb1*nkb2*nkb3*6)  
    integer::ie,jb,ik,ikb1,ikb2,ikb3
    integer::iomp,omp_get_thread_num  
    real(8)::SUM_REAL1,SUM_REAL2  
    !
    integer::NTQ!NTQ=NTK 
    integer::nqb1,nqb2,nqb3,iqb1,iqb2,iqb3 
    real(8)::SQ(3,NTK)!SQ=SK0            
    real(8)::q1,q2,q3 
    integer::ne,je,iq,ikq,ikir,ikqir   
    integer::shift_G(3) 
    real(8)::Rezj,Imzj,sgn 
    complex(8),allocatable::pole_of_chi(:)!pole_of_chi(ne)
    complex(8)::zj 
    complex(8)::SUM_CMPX 
    !
    complex(8),allocatable::ekq_1D(:)!ekq_1D(NTQ)
    complex(8),allocatable::ekq_3D(:,:,:)!ekq_3D(nqb1,nqb2,nqb3) 
    complex(8),allocatable::self(:,:)!self(ndosgrd,NK_irr) 
    complex(8),allocatable::pself(:,:)!pself(ndosgrd,NK_irr) 
    complex(8),allocatable::xow(:)!xow(ndosgrd) 
    complex(8),allocatable::xowt(:,:,:,:)!xowt(ndosgrd,nqb1,nqb2,nqb3) 
    complex(8),allocatable::ca1(:)!ca1(4*ndosgrd) 
    complex(8),allocatable::dosgrd_cmplx(:)!dosgrd_cmplx(ndosgrd) 
    real(8),parameter::pi=dacos(-1.0d0)
    real(8),parameter::au=27.21151d0 
    complex(8),parameter::ci=(0.0D0,1.0D0) 
    !
    OPEN(100,FILE='./dat.pole') 
    rewind(100) 
    read(100,*) ne 
    allocate(pole_of_chi(ne)); pole_of_chi(:)=0.0d0 
    do je=1,ne 
     read(100,'(2f15.10)') pole_of_chi(je) 
    enddo 
    !
    allocate(dosgrd_cmplx(ndosgrd));dosgrd_cmplx=0.0d0   
    do ie=1,ndosgrd 
     dosgrd_cmplx(ie)=dcmplx(dosgrd(ie),0.0d0) 
    enddo
    !
    NTQ=NTK  
    SQ=SK0 
    nqb1=nkb1!important
    nqb2=nkb2!important
    nqb3=nkb3!important
    !
    call make_index_kpt(NTQ,nqb1,nqb2,nqb3,SQ(1,1),index_kpt(1,1,1)) 
    call ttrhdrn_mkidx(nqb1,nqb2,nqb3,imt1(1),b1(1),b2(1),b3(1))
    !
    allocate(self(ndosgrd,NK_irr)); self=0.0d0  
!$OMP PARALLEL PRIVATE(ca1,pself,xow,xowt,ekq_1D,ekq_3D,jb,iomp,ikir,ik,je,zj,Rezj,ImZj,sgn,iq,q1,q2,q3,shift_G,ikq,ikqir,iqb3,iqb2,iqb1,ie,SUM_CMPX) 
    allocate(ca1(4*ndosgrd)); ca1=0.d0   
    allocate(pself(ndosgrd,NK_irr)); pself=0.0d0  
    allocate(xow(ndosgrd)); xow=0.0d0  
    allocate(xowt(ndosgrd,nqb1,nqb2,nqb3)); xowt=0.0d0  
    allocate(ekq_1D(NTQ)); ekq_1D=0.0d0
    allocate(ekq_3D(nqb1,nqb2,nqb3)); ekq_3D=0.0d0  
!$OMP DO 
    do jb=1,NTB
     iomp=omp_get_thread_num() 
     do ikir=1,NK_irr 
      if(iomp.eq.0) write(6,*)'#ikir=',ikir 
      ik=numMK(ikir) 
      xow=0.0d0 
      do je=1,ne 
       zj=pole_of_chi(je) 
       !
       Rezj=dble(zj)
       Imzj=dimag(zj)*1.0d0!0.1d0!0.01d0 
       !
       zj=cmplx(Rezj,Imzj)
       !
       ekq_1D=0.0d0
       do iq=1,NTQ 
        q1=SQ(1,iq)
        q2=SQ(2,iq)
        q3=SQ(3,iq)
        shift_G(:)=0
        call search_kq(NTK,SK0(1,1),-q1,-q2,-q3,ik,ikq,shift_G(1))
        ikqir=numirr(ikq) 
        ekq_1D(iq)=cmplx(E_EIGI(jb,ikqir),0.0d0) 
       enddo!iq 
       ekq_3D=0.0d0 
       do iqb3=1,nqb3
        do iqb2=1,nqb2
         do iqb1=1,nqb1 
          iq=index_kpt(iqb1,iqb2,iqb3) 
          ekq_3D(iqb1,iqb2,iqb3)=ekq_1D(iq)
         enddo!iqb1 
        enddo!iqb2 
       enddo!iqb3 
       !--
       !simple sum
       !
       !xowt=0.0d0 
       !do iqb3=1,nqb3
       ! do iqb2=1,nqb2
       !  do iqb1=1,nqb1 
       !   if(dble(ekq_3D(iqb1,iqb2,iqb3))>FermiEnergy)then 
       !    sgn=1.0 !unocc
       !   else
       !    sgn=-1.0 !occ
       !   endif 
       !   xowt(:,iqb1,iqb2,iqb3)=1.0d0/(dosgrd_cmplx(:)-ekq_3D(iqb1,iqb2,iqb3)-(zj-ci*delt)*sgn)
       !  enddo!iqb1
       ! enddo!iqb2
       !enddo!iqb3
       !--
       !tetrahedron 
       !
       ca1=0.0d0
       xowt=0.0d0 
       call ttrhdrn_causal(dmna,dmnr,nqb1,nqb2,nqb3,imt1(1),ekq_3D(1,1,1),FermiEnergy,delt,ndosgrd,dosgrd_cmplx(1),zj,ca1(1),xowt(1,1,1,1))
       !--
       do ie=1,ndosgrd 
        SUM_CMPX=0.0d0 
        do iqb3=1,nqb3
         do iqb2=1,nqb2 
          do iqb1=1,nqb1 
           SUM_CMPX=SUM_CMPX+xowt(ie,iqb1,iqb2,iqb3) 
          enddo!iqb1
         enddo!iqb2
        enddo!iqb3 
        xow(ie)=xow(ie)+SUM_CMPX/dble(NTQ) 
       enddo!ie 
      enddo!je 
      pself(:,ikir)=pself(:,ikir)+xow(:) 
     enddo!ik
     if(iomp.eq.0) write(6,*)'#',jb   
    enddo!jb 
!$OMP END DO 
!$OMP CRITICAL
    self=self+pself  
!$OMP END CRITICAL
    deallocate(ca1,pself,xow,xowt,ekq_1D,ekq_3D)  
!$OMP END PARALLEL 
    !
    self_i=0.0d0 
    self_r=0.0d0 
    do ie=1,ndosgrd 
     SUM_REAL1=0.0d0 
     SUM_REAL2=0.0d0 
     do ik=1,NTK 
       ikir=numirr(ik) 
       SUM_REAL1=SUM_REAL1+dble(self(ie,ikir))/pi  
       SUM_REAL2=SUM_REAL2+dimag(self(ie,ikir))/pi  
     enddo 
     self_r(ie)=(2.0d0/dble(ncomp))*SUM_REAL1/dble(NTK)!2 is spin 
     self_i(ie)=(2.0d0/dble(ncomp))*SUM_REAL2/dble(NTK)!2 is spin 
    enddo 
    !
    deallocate(pole_of_chi,dosgrd_cmplx,self) 
    !    
    return
  end subroutine calc_dos_self 
  !
  subroutine search_kq(NTK,SK0,q1,q2,q3,ik,ikq,shift_G)
    implicit none 
    integer::NTK
    real(8)::SK0(3,NTK)  
    real(8)::q1,q2,q3
    integer::ik
    integer::jk
    real(8)::SKQ(3)
    integer::ikq,shift_G(3)
    real(8),parameter::dlt_BZ=1.0d-6!20170322  
    !
    SKQ(1)=SK0(1,ik)+q1
    SKQ(2)=SK0(2,ik)+q2 
    SKQ(3)=SK0(3,ik)+q3 
    !
    if(SKQ(1)>1.50d0+dlt_BZ)then 
     SKQ(1)=SKQ(1)-2.0D0 
     shift_G(1)=+2 
    endif 
    if(SKQ(1)>0.5D0+dlt_BZ)then 
     SKQ(1)=SKQ(1)-1.0D0 
     shift_G(1)=+1
    endif 
    if(SKQ(1)<=-1.5D0+dlt_BZ)then 
     SKQ(1)=SKQ(1)+2.0D0 
     shift_G(1)=-2 
    endif 
    if(SKQ(1)<=-0.5D0+dlt_BZ)then 
     SKQ(1)=SKQ(1)+1.0D0 
     shift_G(1)=-1
    endif 
    !
    if(SKQ(2)>1.5D0+dlt_BZ)then 
     SKQ(2)=SKQ(2)-2.0D0 
     shift_G(2)=+2 
    endif 
    if(SKQ(2)>0.5D0+dlt_BZ)then 
     SKQ(2)=SKQ(2)-1.0D0 
     shift_G(2)=+1
    endif 
    if(SKQ(2)<=-1.5D0+dlt_BZ)then 
     SKQ(2)=SKQ(2)+2.0D0 
     shift_G(2)=-2 
    endif 
    if(SKQ(2)<=-0.5D0+dlt_BZ)then 
     SKQ(2)=SKQ(2)+1.0D0 
     shift_G(2)=-1
    endif 
    !
    if(SKQ(3)>1.5D0+dlt_BZ)then 
     SKQ(3)=SKQ(3)-2.0D0 
     shift_G(3)=+2 
    endif 
    if(SKQ(3)>0.5D0+dlt_BZ)then 
     SKQ(3)=SKQ(3)-1.0D0 
     shift_G(3)=+1
    endif 
    if(SKQ(3)<=-1.5D0+dlt_BZ)then 
     SKQ(3)=SKQ(3)+2.0D0 
     shift_G(3)=-2 
    endif 
    if(SKQ(3)<=-0.5D0+dlt_BZ)then 
     SKQ(3)=SKQ(3)+1.0D0 
     shift_G(3)=-1
    endif 
    !
    do jk=1,NTK
     if(ABS(SK0(1,jk)-SKQ(1))<1.D-6.and.ABS(SK0(2,jk)-SKQ(2))<1.D-6.and.ABS(SK0(3,jk)-SKQ(3))<1.D-6)then 
      ikq=jk 
     endif 
    enddo 
    RETURN  
  end subroutine search_kq 
  !
end module m_self 
