module m_dos 
  implicit none
contains
  !
  subroutine calculate_dos(NTB,NTK,nkb1,nkb2,nkb3,electron_number,FermiEnergy,filename,b1,b2,b3,SK0,EIG,VKS)
    implicit none 
    integer,intent(in)::NTB,NTK,nkb1,nkb2,nkb3
    real(8),intent(in)::electron_number
    real(8),intent(in)::b1(3),b2(3),b3(3)
    real(8),intent(in)::SK0(3,NTK)
    real(8),intent(in)::EIG(NTB,NTK) 
    real(8),intent(inout)::FermiEnergy
    character(99),intent(in)::filename 
    !
    complex(8),intent(in)::VKS(NTB,NTB,NTK)
    real(8),allocatable::pdos(:,:)!pdos(ndosgrd,NTB) 
    real(8)::flg_pdos
    !
    real(8)::emax!=maxval(EIG)
    real(8)::emin!=minval(EIG)
    integer::ndosgrd!=int(2.0d0*diff/dlt)+1
    real(8),allocatable::dosgrd(:)!dosgrd(ndosgrd) 
    real(8),allocatable::dos(:)!dos(ndosgrd) 
    real(8),allocatable::efline(:)!efline(ndosgrd)   
    !
    real(8),parameter::au=27.21151d0
    !real(8),parameter::delt=0.005d0/au!Greens function delt in au 
    real(8),parameter::delt=0.01d0/au!Greens function delt in au 
    real(8),parameter::dmna=1.0d-3!Ttrhdrn parameter dmna in au 
    real(8),parameter::dmnr=1.0d-3!Ttrhdrn parameter dmnr in au 
    real(8),parameter::delw=2.0d0*delt!Grid width in au 
    ! 
    !dos-grid
    !
    call est_ndosgrd(NTB,NTK,EIG(1,1),delw,emin,emax,ndosgrd)  
    allocate(dosgrd(ndosgrd));dosgrd=0.0d0 
    call make_dosgrd(emin,delw,ndosgrd,dosgrd(1))  
    !
    !calc dos 
    !
    allocate(dos(ndosgrd));dos=0.0d0 
    call calc_dos(NTB,NTK,nkb1,nkb2,nkb3,ndosgrd,dosgrd(1),EIG(1,1),SK0(1,1),delt,dmnr,dmna,b1(1),b2(1),b3(1),dos(1)) 
    !
    !calc pdos 
    !
    flg_pdos=maxval(abs(VKS)) 
    if(flg_pdos/=0.0d0)then
      write(6,'(a40,f15.8,a)')'flg_pdos/=0.0:',flg_pdos,': PDOS calculate:' 
      allocate(pdos(ndosgrd,NTB));pdos=0.0d0 
      call calc_pdos(NTB,NTK,nkb1,nkb2,nkb3,ndosgrd,dosgrd(1),EIG(1,1),VKS(1,1,1),SK0(1,1),delt,dmnr,dmna,b1(1),b2(1),b3(1),pdos(1,1)) 
    endif 
    !
    !estimate ef 
    !
    if(electron_number/=0.0d0)then
     call est_ef(ndosgrd,delw,electron_number,dosgrd(1),dos(1),FermiEnergy) 
    endif 
    allocate(efline(ndosgrd));efline=0.0d0 
    call make_efline(ndosgrd,FermiEnergy,delw,dosgrd(1),dos(1),efline(1)) 
    !
    !wrt dat.dos 
    !
    call wrt_dos(filename,ndosgrd,dosgrd(1),dos(1),efline(1)) 
    !
    !wrt dat.pdos 
    !
    if(flg_pdos/=0.0d0)then
      call wrt_pdos(filename,NTB,ndosgrd,dosgrd(1),pdos(1,1),efline(1)) 
    endif 
    !
    return
  end subroutine calculate_dos 
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
    !
    SUM_REAL=0.0d0 
    do ie=1,ndosgrd-1
     !--
     !SUM_REAL=SUM_REAL+deltw*dos(ie) 
     !--
     !Trapezoidal 
     SUM_REAL=SUM_REAL+deltw*(dos(ie)+dos(ie+1))/2.0d0  
    enddo 
    !
    write(6,'(a40)')'+++ m_dos: est_ef +++'
    write(6,'(a40,f15.8)')'SUM of DOS',SUM_REAL
    write(6,'(a40,f15.8)')'FermiEnergy(before)=',FermiEnergy*au  
    SUM_REAL=0.0d0 
    do ie=1,ndosgrd-1
     if(SUM_REAL>=electron_number)goto 3000 
       !--
       !SUM_REAL=SUM_REAL+deltw*dos(ie) 
       !--
       !Trapezoidal 
       SUM_REAL=SUM_REAL+deltw*(dos(ie)+dos(ie+1))/2.0d0  
    enddo 
  3000 FermiEnergy=dosgrd(ie) 
    write(6,'(a40,f15.8)')'FermiEnergy(after)=',dosgrd(ie)*au  
    write(6,'(a40,f15.8)')'electron_number',SUM_REAL
    return
  end subroutine 
  !
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
    dosmax=1.2d0*maxval(abs(dos(:)))/au  
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
    do ie=1,ndosgrd-1 
     !--
     !SUM_REAL=SUM_REAL+dos(ie)*deltw 
     !--
     !Trapezoidal 
     SUM_REAL=SUM_REAL+deltw*(dos(ie)+dos(ie+1))/2.0d0  
     !write(6,'(2f15.10)') dosgrd(ie)*au,dos(ie)/au  
    enddo 
    N0ttr=dos(ief)/2.0d0 !2 is spin
    write(6,*)
    write(6,'(a40)')'+++ m_dos: est_efline +++'
    write(6,'(a40,f15.8)')'Integral dos(w) dw',SUM_REAL 
    write(6,'(a40,f15.8)')'N(0) in au per calc cell',N0ttr 
    write(6,'(a40,f15.8)')'N(0) in eV per calc cell',N0ttr/au   
    write(6,*)
    return
  end subroutine 
  !
  subroutine calc_dos(NTB,NTK,nkb1,nkb2,nkb3,ndosgrd,dosgrd,E_EIG,SK0,delt,dmnr,dmna,b1,b2,b3,dos) 
    use m_tetrahedron
    implicit none
    integer,intent(in)::NTB,NTK,nkb1,nkb2,nkb3,ndosgrd
    real(8),intent(in)::dosgrd(ndosgrd) 
    real(8),intent(in)::E_EIG(NTB,NTK)           
    real(8),intent(in)::SK0(3,NTK)           
    real(8),intent(in)::delt,dmnr,dmna 
    real(8),intent(in)::b1(3),b2(3),b3(3)
    real(8),intent(out)::dos(ndosgrd) 
    integer::index_kpt(nkb1,nkb2,nkb3) 
    integer::imt1(4*nkb1*nkb2*nkb3*6)  
    integer::ie,jb,ik,ikb1,ikb2,ikb3
    integer::iomp,omp_get_thread_num  
    real(8)::SUM_REAL 
    complex(8)::fk_1D(NTK)
    complex(8)::fk_3D(nkb1,nkb2,nkb3) 
    complex(8)::gk_1D(NTK)
    complex(8)::gk_3D(nkb1,nkb2,nkb3)
    complex(8)::xo(nkb1,nkb2,nkb3) 
    complex(8)::xow(ndosgrd,nkb1,nkb2,nkb3)
    real(8),parameter::pi=dacos(-1.0d0)
    real(8),parameter::au=27.21151d0 
    !
    call make_index_kpt(NTK,nkb1,nkb2,nkb3,SK0(1,1),index_kpt(1,1,1)) 
    call ttrhdrn_mkidx(nkb1,nkb2,nkb3,imt1(1),b1(1),b2(1),b3(1))
    xow=0.0d0 
!$OMP PARALLEL PRIVATE(ie,jb,fk_1D,gk_1D,ik,fk_3D,gk_3D,&
!$OMP ikb3,ikb2,ikb1,xo,iomp) 
!$OMP DO 
    do ie=1,ndosgrd 
     do jb=1,NTB
      fk_1D=0.0d0
      gk_1D=0.0d0
      do ik=1,NTK 
       fk_1D(ik)=1.0d0 
       gk_1D(ik)=cmplx(dosgrd(ie)-E_EIG(jb,ik),-delt) 
      enddo!ik 
      fk_3D=0.0d0 
      gk_3D=0.0d0 
      do ikb3=1,nkb3
       do ikb2=1,nkb2
        do ikb1=1,nkb1 
         ik=index_kpt(ikb1,ikb2,ikb3) 
         fk_3D(ikb1,ikb2,ikb3)=fk_1D(ik)
         gk_3D(ikb1,ikb2,ikb3)=gk_1D(ik)
        enddo!ikb1 
       enddo!ikb2 
      enddo!ikb3 
      xo=0.0d0 
      call ttrhdrn_simple(dmna,dmnr,nkb1,nkb2,nkb3,imt1(1),& 
       fk_3D(1,1,1),gk_3D(1,1,1),xo(1,1,1))
      xow(ie,:,:,:)=xow(ie,:,:,:)+xo(:,:,:) 
     enddo!jb 
     iomp=omp_get_thread_num() 
     !if(iomp.eq.0) write(6,*) '#',ie  
    enddo!ie
!$OMP END DO 
!$OMP END PARALLEL 
    dos=0.0d0 
    do ie=1,ndosgrd 
     SUM_REAL=0.0d0 
     do ikb3=1,nkb3
      do ikb2=1,nkb2
       do ikb1=1,nkb1 
        SUM_REAL=SUM_REAL+dabs(dimag(xow(ie,ikb1,ikb2,ikb3)))/pi  
       enddo 
      enddo 
     enddo 
     dos(ie)=2.0d0*SUM_REAL/dble(NTK)!2 is spin 
    enddo 
    !
    return
  end subroutine calc_dos
  !
  subroutine calc_pdos(NTB,NTK,nkb1,nkb2,nkb3,ndosgrd,dosgrd,E_EIG,VKS,SK0,delt,dmnr,dmna,b1,b2,b3,pdos) 
    use m_tetrahedron
    implicit none
    integer,intent(in)::NTB,NTK,nkb1,nkb2,nkb3,ndosgrd
    real(8),intent(in)::dosgrd(ndosgrd) 
    real(8),intent(in)::E_EIG(NTB,NTK)           
    real(8),intent(in)::SK0(3,NTK)           
    real(8),intent(in)::delt,dmnr,dmna 
    real(8),intent(in)::b1(3),b2(3),b3(3)
    !
    complex(8),intent(in)::VKS(NTB,NTB,NTK)
    complex(8)::pxow(ndosgrd,NTB,nkb1,nkb2,nkb3)
    real(8),intent(out)::pdos(ndosgrd,NTB) 
    !
    integer::index_kpt(nkb1,nkb2,nkb3) 
    integer::imt1(4*nkb1*nkb2*nkb3*6)  
    integer::ie,ib,jb,ik,ikb1,ikb2,ikb3
    integer::iomp,omp_get_thread_num  
    real(8)::SUM_REAL 
    complex(8)::fk_1D(NTK)
    complex(8)::fk_3D(nkb1,nkb2,nkb3) 
    complex(8)::gk_1D(NTK)
    complex(8)::gk_3D(nkb1,nkb2,nkb3)
    complex(8)::xo(nkb1,nkb2,nkb3) 
    real(8),parameter::pi=dacos(-1.0d0)
    real(8),parameter::au=27.21151d0 
    !
    call make_index_kpt(NTK,nkb1,nkb2,nkb3,SK0(1,1),index_kpt(1,1,1)) 
    call ttrhdrn_mkidx(nkb1,nkb2,nkb3,imt1(1),b1(1),b2(1),b3(1))
    pxow=0.0d0 
!$OMP PARALLEL PRIVATE(ie,jb,fk_1D,gk_1D,ik,fk_3D,gk_3D,&
!$OMP ikb3,ikb2,ikb1,xo,iomp) 
!$OMP DO 
    do ie=1,ndosgrd 
     do jb=1,NTB
      fk_1D=0.0d0
      gk_1D=0.0d0
      do ik=1,NTK 
       fk_1D(ik)=1.0d0 
       gk_1D(ik)=cmplx(dosgrd(ie)-E_EIG(jb,ik),-delt) 
      enddo!ik 
      fk_3D=0.0d0 
      gk_3D=0.0d0 
      do ikb3=1,nkb3
       do ikb2=1,nkb2
        do ikb1=1,nkb1 
         ik=index_kpt(ikb1,ikb2,ikb3) 
         fk_3D(ikb1,ikb2,ikb3)=fk_1D(ik)
         gk_3D(ikb1,ikb2,ikb3)=gk_1D(ik)
        enddo!ikb1 
       enddo!ikb2 
      enddo!ikb3 
      xo=0.0d0 
      call ttrhdrn_simple(dmna,dmnr,nkb1,nkb2,nkb3,imt1(1),& 
       fk_3D(1,1,1),gk_3D(1,1,1),xo(1,1,1))
      pxow(ie,jb,:,:,:)=xo(:,:,:) 
     enddo!jb 
     iomp=omp_get_thread_num() 
     !if(iomp.eq.0) write(6,*) '#',ie  
    enddo!ie
!$OMP END DO 
!$OMP END PARALLEL 
    pdos=0.0d0 
    do ie=1,ndosgrd 
     do ib=1,NTB
      SUM_REAL=0.0d0 
      do jb=1,NTB 
       do ikb3=1,nkb3
        do ikb2=1,nkb2
         do ikb1=1,nkb1 
          ik=index_kpt(ikb1,ikb2,ikb3) 
          SUM_REAL=SUM_REAL+(abs(VKS(jb,ib,ik))**2)*dabs(dimag(pxow(ie,jb,ikb1,ikb2,ikb3)))/pi  
         enddo!ikb1 
        enddo!ikb2 
       enddo!ikb3 
      enddo!jb  
      pdos(ie,ib)=2.0d0*SUM_REAL/dble(NTK)!2 is spin 
     enddo!ib 
    enddo!ie 
    !
    return
  end subroutine calc_pdos 
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
  subroutine wrt_dos(filename,ndosgrd,dosgrd,dos,efline) 
    implicit none
    integer,intent(in)::ndosgrd
    real(8),intent(in)::dosgrd(ndosgrd)
    real(8),intent(in)::dos(ndosgrd)
    real(8),intent(in)::efline(ndosgrd)
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
    do ie=1,ndosgrd
     write(300,'(3f15.10)') dosgrd(ie)*au,dos(ie)/au,efline(ie)  
    enddo!ie 
    close(300)
    return
  end subroutine 
  !
  subroutine wrt_pdos(filename,NTB,ndosgrd,dosgrd,pdos,efline) 
    implicit none
    integer,intent(in)::NTB
    integer,intent(in)::ndosgrd
    real(8),intent(in)::dosgrd(ndosgrd)
    real(8),intent(in)::pdos(ndosgrd,NTB)
    real(8),intent(in)::efline(ndosgrd)
    character(99),intent(in)::filename
    character(99)::fname,fname_arg 
    integer::ie,ib 
    real(8),parameter::au=27.21151d0 
    !
    !OPEN(400,W,file='./dir-wan/dat.dos.wannier-.ib')
    !
    do ib=1,NTB 
     fname_arg='-'
     write(fname,'(a,a,i3.3)') trim(filename),trim(fname_arg),ib
     OPEN(400,FILE=TRIM(fname)) 
     rewind(400)  
     do ie=1,ndosgrd
      write(400,'(3f15.10)') dosgrd(ie)*au,pdos(ie,ib)/au,efline(ie)  
     enddo!ie 
     close(400)
    enddo!ib 
    return
  end subroutine wrt_pdos 
  !
end module m_dos 
