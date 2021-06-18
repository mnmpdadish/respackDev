module m_self 
  implicit none
contains
  !
  subroutine calculate_self(ncomp,NTB,NTK,NTQ,nkb1,nkb2,nkb3,NK_irr,numirr,numMK,FermiEnergy,nsgm,sgmw,ne,pole_of_chi,E_EIGI,SK0,SQ,b1,b2,b3,&
  filename,nproc,pnq,bnq,enq,myrank,rec_len) 
    implicit none 
    integer,intent(in)::ncomp,NTB,NTK,NTQ,nkb1,nkb2,nkb3,NK_irr,nsgm,ne,nproc,pnq,bnq,enq,myrank,rec_len
    integer,intent(in)::numirr(NTK)
    integer,intent(in)::numMK(NK_irr)
    real(8),intent(inout)::FermiEnergy
    character(99),intent(in)::filename 
    real(8),intent(in)::b1(3),b2(3),b3(3)
    real(8),intent(in)::SK0(3,NTK),SQ(3,NTQ) 
    real(8),intent(in)::E_EIGI(NTB,NK_irr) 
    real(8),intent(in)::sgmw(nsgm) 
    complex(8),intent(in)::pole_of_chi(ne)!pole_of_chi(ne)
    !
    real(8),allocatable::dos_i(:)!dos_i(nsgm) 
    real(8),allocatable::dos_r(:)!dos_r(nsgm) 
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
    !calc self 
    !
    allocate(dos_r(nsgm));dos_r=0.0d0 
    allocate(dos_i(nsgm));dos_i=0.0d0 
    call calc_dos_self(ncomp,NTB,NTK,NTQ,nkb1,nkb2,nkb3,NK_irr,numirr(1),numMK(1),FermiEnergy,nsgm,sgmw(1),ne,pole_of_chi(1),E_EIGI(1,1),SK0(1,1),SQ(1,1),&
    delt,dmnr,dmna,b1(1),b2(1),b3(1),nproc,pnq,bnq,enq,myrank,rec_len,dos_r(1),dos_i(1));flg_causal=1
    !
    !wrt dat.self.imag 
    !
    if(myrank.eq.0)then 
     call wrt_dos_imag(filename,nsgm,sgmw(1),dos_i(1),FermiEnergy) 
     !
     !wrt dat.dos.real  
     !
     if(flg_causal==1)then
      call wrt_dos_real(filename,nsgm,sgmw(1),dos_r(1),FermiEnergy) 
     endif 
    endif 
    !
    return
  end subroutine calculate_self 
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
  subroutine calc_dos_self(ncomp,NTB,NTK,NTQ,nkb1,nkb2,nkb3,NK_irr,numirr,numMK,FermiEnergy,nsgm,sgmw,ne,pole_of_chi,E_EIGI,SK0,SQ,delt,dmnr,dmna,b1,b2,b3,&
  nproc,pnq,bnq,enq,myrank,rec_len,self_r,self_i)
    use m_tetrahedron
    implicit none
    integer,intent(in)::ncomp,NTB,NTK,NTQ,nkb1,nkb2,nkb3,NK_irr,nsgm,ne
    integer,intent(in)::numirr(NTK) 
    integer,intent(in)::numMK(NK_irr) 
    integer,intent(in)::nproc,pnq,bnq,enq,myrank 
    integer,intent(in)::rec_len
    real(8),intent(in)::FermiEnergy 
    real(8),intent(in)::sgmw(nsgm) 
    real(8),intent(in)::E_EIGI(NTB,NK_irr)           
    real(8),intent(in)::SK0(3,NTK),SQ(3,NTQ) 
    real(8),intent(in)::delt,dmnr,dmna 
    real(8),intent(in)::b1(3),b2(3),b3(3)
    complex(8),intent(in)::pole_of_chi(ne)!pole_of_chi(ne)
    real(8),intent(out)::self_r(nsgm) 
    real(8),intent(out)::self_i(nsgm) 
    !
    integer::index_kpt(nkb1,nkb2,nkb3) 
    integer::imt1(4*nkb1*nkb2*nkb3*6)  
    integer::ie,ib,ik,ikb1,ikb2,ikb3,iqb1,iqb2,iqb3 
    integer::je,iq,ikq,ikir,ikqir   
    integer::iomp,omp_get_thread_num,impi,file_num   
    integer::nqb1,nqb2,nqb3
    integer::rec_num 
    real(8)::SUM_REAL1,SUM_REAL2  
    real(8)::q1,q2,q3 
    integer::shift_G(3) 
    real(8)::Rezj,Imzj,sgn 
    complex(8)::zj 
    complex(8)::SUM_CMPX 
    character(99)::filename,command 
    !
    complex(8),allocatable::ekq_1D(:)!ekq_1D(NTQ)
    complex(8),allocatable::ekq_3D(:,:,:)!ekq_3D(nqb1,nqb2,nqb3) 
    complex(8),allocatable::self(:,:)!self(nsgm,NK_irr) 
    complex(8),allocatable::pself(:,:)!pself(nsgm,NK_irr) 
    complex(8),allocatable::xow(:)!xow(nsgm) 
    complex(8),allocatable::xowt(:,:,:,:)!xowt(nsgm,nqb1,nqb2,nqb3) 
    complex(8),allocatable::xowt_1D(:,:)!xowt_1D(nsgm,NTQ)
    complex(4),allocatable::xowtj(:,:,:)!xowtj(ne,nsgm,NTQ) 
    complex(8),allocatable::ca1(:)!ca1(4*nsgm) 
    complex(8),allocatable::sgmw_cmplx(:)!sgmw_cmplx(nsgm) 
    real(8),parameter::pi=dacos(-1.0d0)
    real(8),parameter::au=27.21151d0 
    complex(8),parameter::ci=(0.0D0,1.0D0) 
    !--
    nqb1=nkb1!important
    nqb2=nkb2!important
    nqb3=nkb3!important
    !--
    call make_index_kpt(NTQ,nqb1,nqb2,nqb3,SQ(1,1),index_kpt(1,1,1)) 
    call ttrhdrn_mkidx(nqb1,nqb2,nqb3,imt1(1),b1(1),b2(1),b3(1))
    !--
    allocate(sgmw_cmplx(nsgm)); sgmw_cmplx=0.0d0   
    do ie=1,nsgm
     sgmw_cmplx(ie)=dcmplx(sgmw(ie),0.0d0) 
    enddo
    !--
    !do impi=0,nproc-1
    !write(command,"('mkdir -p ./xqdata',i3.3)")myrank!impi  
    !call system(command) 
    !enddo!impi  
    !--
    !rec_len=ne*nsgm*2 
    !allocate(xowtj(ne,nsgm,1)); xowtj=0.0d0 
    !inquire(iolength=rec_len)((xowtj(je,ie,1),je=1,ne),ie=1,nsgm) 
    !write(6,'(a,i30)')'rec_len (byte):',rec_len 
    !deallocate(xowtj) 
    !--
    allocate(self(nsgm,NK_irr)); self=0.0d0  
!$OMP PARALLEL PRIVATE(ca1,pself,xow,xowt,ekq_1D,ekq_3D,ib,iomp,ikir,ik,je,zj,Rezj,ImZj,sgn,iq,q1,q2,q3,shift_G,ikq,ikqir,iqb3,iqb2,iqb1,ie,SUM_CMPX,&
!$OMP&                 xowt_1D,xowtj,filename,file_num,rec_num)
    allocate(ca1(4*nsgm)); ca1=0.d0   
    allocate(pself(nsgm,NK_irr)); pself=0.0d0  
    allocate(xow(nsgm)); xow=0.0d0  
    allocate(xowt(nsgm,nqb1,nqb2,nqb3)); xowt=0.0d0  
    allocate(ekq_1D(NTQ)); ekq_1D=0.0d0
    allocate(ekq_3D(nqb1,nqb2,nqb3)); ekq_3D=0.0d0  
    allocate(xowt_1D(nsgm,NTQ)); xowt_1D=0.0d0  
    allocate(xowtj(ne,nsgm,pnq)); xowtj=0.0d0 
!$OMP DO 
    do ib=1,NTB
     iomp=omp_get_thread_num() 
     !--
     !file open
     !--
     !do impi=0,nproc-1 
     !write(filename,"('/var/tmp/xqdata',i3.3,'/dat.ib',i4.4)")myrank,ib 
     write(filename,"('./dir-gw/xqdata',i3.3,'/dat.ib',i4.4)")myrank,ib 
     file_num=(myrank+1)*10000+ib  
     open(file_num,FILE=filename,FORM='unformatted',access='direct',recl=rec_len) 
     !open(file_num,FILE=filename,FORM='unformatted')
     !enddo!impi 
     !
     do ikir=1,NK_irr 
      if(iomp.eq.0.and.myrank.eq.0) write(6,*)'#ikir=',ikir 
      ik=numMK(ikir) 
      xowtj=0.0d0 
      xow=0.0d0 
      do je=1,ne 
       zj=pole_of_chi(je) 
       !
       Rezj=dble(zj)
       Imzj=dimag(zj)*1.0d0 !0.1d0!0.01d0 
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
        ekq_1D(iq)=cmplx(E_EIGI(ib,ikqir),0.0d0) 
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
       !--
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
       !--
       ca1=0.0d0
       xowt=0.0d0 
       call ttrhdrn_causal(dmna,dmnr,nqb1,nqb2,nqb3,imt1(1),ekq_3D(1,1,1),FermiEnergy,delt,nsgm,sgmw_cmplx(1),zj,ca1(1),xowt(1,1,1,1))
       !--
       do ie=1,nsgm 
        SUM_CMPX=0.0d0 
        do iqb3=1,nqb3
         do iqb2=1,nqb2 
          do iqb1=1,nqb1 
           SUM_CMPX=SUM_CMPX+xowt(ie,iqb1,iqb2,iqb3) 
           iq=index_kpt(iqb1,iqb2,iqb3) 
           xowt_1D(:,iq)=xowt(:,iqb1,iqb2,iqb3) 
          enddo!iqb1
         enddo!iqb2
        enddo!iqb3 
        xow(ie)=xow(ie)+SUM_CMPX/dble(NTQ) 
       enddo!ie 
       !--
       !make xowtj
       !--
       do ie=1,nsgm
        do iq=1,pnq !NTQ 
         xowtj(je,ie,iq)=xowt_1D(ie,bnq+iq-1) 
         !xowtj(je,ie,iq)=xowt_1D(ie,iq) 
        enddo!iq
       enddo!ie 
       !--
      enddo!je
      !--
      !wrt xowtj via fileIO 
      !--
      !do impi=0,nproc-1
      ! pnqIO=NTQ/nproc 
      ! if(impi.lt.mod(NTQ,nproc))then
      !  pnqIO=pnqIO+1
      !  bnqIO=pnqIO*impi+1
      !  enqIO=bnqIO+pnqIO-1
      ! else
      !  bnqIO=(pnqIO+1)*mod(NTQ,nproc)+pnqIO*(impi-mod(NTQ,nproc))+1
      !  enqIO=bnqIO+pnqIO-1
      ! endif 
      ! file_num=(impi+1)*10000+ib 
      do iq=1,pnq
       rec_num=iq+(ikir-1)*pnq
       !write(file_num,rec=rec_num)((xowtj(je,ie,bnq+iq-1),je=1,ne),ie=1,nsgm) 
       write(file_num,rec=rec_num)((xowtj(je,ie,iq),je=1,ne),ie=1,nsgm) 
       !write(file_num)((xowtj(je,ie,iq),je=1,ne),ie=1,nsgm) 
      enddo!iq
      !enddo!impi 
      !
      pself(:,ikir)=pself(:,ikir)+xow(:) 
      !
     enddo!ikir 
     !--
     !file close 
     !--
     !do impi=0,nproc-1 
     ! file_num=(impi+1)*10000+ib 
     close(file_num) 
     !enddo!impi 
     !
     if(iomp.eq.0) write(6,*)'#',ib   
     !
    enddo!ib 
!$OMP END DO 
!$OMP CRITICAL
    self=self+pself  
!$OMP END CRITICAL
    deallocate(ca1,pself,xow,xowt,ekq_1D,ekq_3D,xowt_1D,xowtj)  
!$OMP END PARALLEL 
    !
    self_i=0.0d0 
    self_r=0.0d0 
    do ie=1,nsgm 
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
    deallocate(sgmw_cmplx,self) 
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
