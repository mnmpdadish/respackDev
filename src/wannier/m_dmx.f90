module m_dmx 
  implicit none
contains
  !
  subroutine calculate_dmx(NTB,NTK,NWF,nkb1,nkb2,nkb3,Na1,Na2,Na3,FermiEnergy,a1,a2,a3,b1,b2,b3,SK0,EIG,UNT)
    implicit none 
    integer,intent(in)::NTB,NWF,NTK,nkb1,nkb2,nkb3,Na1,Na2,Na3 
    real(8),intent(in)::FermiEnergy 
    real(8),intent(in)::a1(3),a2(3),a3(3)
    real(8),intent(in)::b1(3),b2(3),b3(3)
    real(8),intent(in)::SK0(3,NTK) 
    real(8),intent(in)::EIG(NTB,NTK)
    complex(8),intent(in)::UNT(NTB,NWF,NTK) 
    !
    integer,allocatable::lat_num_a1(:),lat_num_a2(:),lat_num_a3(:) 
    real(8),allocatable::WEIGHT_R(:,:,:) 
    real(8),allocatable::xow(:,:,:,:)
    complex(8),allocatable::pf(:,:,:,:) 
    complex(8),allocatable::dmx(:,:,:,:,:) 
    integer::NR 
    integer,parameter::flg_weight=1 
    !
    !1. lat_num_a%
    !
    NR=(2*Na1+1)*(2*Na2+1)*(2*Na3+1) 
    allocate(lat_num_a1(NR)); lat_num_a1=0  
    allocate(lat_num_a2(NR)); lat_num_a2=0  
    allocate(lat_num_a3(NR)); lat_num_a3=0  
    call make_lat_num(Na1,Na2,Na3,NR,lat_num_a1(1),lat_num_a2(1),lat_num_a3(1)) 
    !
    !2. xow(ib,ik) 
    !
    allocate(xow(NTB,nkb1,nkb2,nkb3)); xow=0.0d0 
    call make_xow(nkb1,nkb2,nkb3,NTB,NTK,FermiEnergy,b1(1),b2(1),b3(1),SK0(1,1),EIG(1,1),xow(1,1,1,1)) 
    !
    !3. WEIGHT_R
    !
    allocate(WEIGHT_R(-Na1:Na1,-Na2:Na2,-Na3:Na3)); WEIGHT_R=1.0d0
    call make_WEIGHT_R(nkb1,nkb2,nkb3,Na1,Na2,Na3,NTK,flg_weight,WEIGHT_R(-Na1,-Na2,-Na3)) 
    !
    !4. phase factor
    !
    allocate(pf(-Na1:Na1,-Na2:Na2,-Na3:Na3,NTK)); pf=0.0d0 
    call make_pf(Na1,Na2,Na3,nkb1,nkb2,nkb3,NTK,a1(1),a2(1),a3(1),SK0(1,1),pf(-Na1,-Na2,-Na3,1)) 
    !
    !5. density matrix 
    !
    allocate(dmx(NWF,NWF,-Na1:Na1,-Na2:Na2,-Na3:Na3)); dmx=0.0d0   
    call make_dmx(Na1,Na2,Na3,NR,NTB,NTK,NWF,nkb1,nkb2,nkb3,lat_num_a1(1),lat_num_a2(1),lat_num_a3(1),&
    SK0(1,1),UNT(1,1,1),xow(1,1,1,1),pf(-Na1,-Na2,-Na3,1),WEIGHT_R(-Na1,-Na2,-Na3),dmx(1,1,-Na1,-Na2,-Na3)) 
    !
    !6. write density matrix
    call wrt_model_dr(NWF,Na1,Na2,Na3,dmx(1,1,-Na1,-Na2,-Na3))
    ! 
    deallocate(WEIGHT_R,pf,xow,lat_num_a1,lat_num_a2,lat_num_a3,dmx) 
    !
    return 
  end subroutine calculate_dmx 

  subroutine make_lat_num(Na1,Na2,Na3,NR,lat_num_a1,lat_num_a2,lat_num_a3) 
    implicit none
    integer,intent(in)::Na1,Na2,Na3,NR 
    integer,intent(out)::lat_num_a1(NR)
    integer,intent(out)::lat_num_a2(NR)
    integer,intent(out)::lat_num_a3(NR)
    integer::ia1,ia2,ia3,iR 
    !
    do ia1=-Na1,Na1
     do ia2=-Na2,Na2
      do ia3=-Na3,Na3
       iR=(ia3+Na3+1)+(ia2+Na2)*(2*Na3+1)+(ia1+Na1)*(2*Na3+1)*(2*Na2+1) 
       lat_num_a1(iR)=ia1
       lat_num_a2(iR)=ia2
       lat_num_a3(iR)=ia3 
      enddo
     enddo
    enddo 
    !
    return 
  end subroutine make_lat_num 

  subroutine make_xow(nkb1,nkb2,nkb3,NTB,NTK,FermiEnergy,b1,b2,b3,SK0,E_EIG,xow) 
    use m_tetrainteg
    implicit none 
    integer,intent(in)::nkb1,nkb2,nkb3,NTB,NTK
    real(8),intent(in)::FermiEnergy 
    real(8),intent(in)::b1(3),b2(3),b3(3)
    real(8),intent(in)::SK0(3,NTK) 
    real(8),intent(in)::E_EIG(NTB,NTK) 
    real(8),intent(out)::xow(NTB,nkb1,nkb2,nkb3) 
    !
    integer::index_kpt(nkb1,nkb2,nkb3) 
    integer::imt1(4*nkb1*nkb2*nkb3*6)  
    complex(8)::fk_3D(nkb1,nkb2,nkb3) 
    complex(8)::gk_1D(NTK)
    complex(8)::gk_3D(nkb1,nkb2,nkb3)
    complex(8)::xo(nkb1,nkb2,nkb3)  
    complex(8)::ca1(4)
    complex(8)::omega(1)  
    integer::ib,ikb1,ikb2,ikb3,ik 
    real(8)::SUM_REAL 
    real(8),parameter::au=27.21151d0
    !real(8),parameter::delt=0.005d0/au!Greens function delt in au 
    real(8),parameter::delt=0.01d0/au!Greens function delt in au 
    real(8),parameter::dmna=1.0d-3!Ttrhdrn parameter dmna in au 
    real(8),parameter::dmnr=1.0d-3!Ttrhdrn parameter dmnr in au 
    real(8),parameter::pi=dacos(-1.0d0)
    !
    index_kpt=0 
    call make_index_kpt(NTK,nkb1,nkb2,nkb3,SK0(1,1),index_kpt(1,1,1)) 
    !
    imt1=0 
    call ttritg_mkidx(nkb1,nkb2,nkb3,imt1(1),b1(1),b2(1),b3(1)) 
    !
    omega(1)=cmplx(FermiEnergy,0.0d0)
    !
    !   !$OMP PARALLEL PRIVATE(ik,ib,gk_1D,fk_3D,gk_3D,ikb3,ikb2,ikb1,xo,ca1) 
    !   !$OMP DO 
    do ib=1,NTB
     do ik=1,NTK 
      gk_1D(ik)=cmplx(-E_EIG(ib,ik),-delt) 
     enddo 
     gk_3D=0.0d0 
     do ikb3=1,nkb3
      do ikb2=1,nkb2
       do ikb1=1,nkb1 
        ik=index_kpt(ikb1,ikb2,ikb3)
        gk_3D(ikb1,ikb2,ikb3)=gk_1D(ik) 
       enddo
      enddo
     enddo 
     fk_3D=1.0d0 
     xo=0.0d0 
     ca1=0.0d0 
     call ttritg_sum(dmna,dmnr,nkb1,nkb2,nkb3,imt1(1),fk_3D(1,1,1),gk_3D(1,1,1),1,omega(1),ca1(1),xo(1,1,1)) 
     xow(ib,:,:,:)=dabs(dimag(xo(:,:,:)))/pi  
    enddo!ib 
    !    !$OMP END DO 
    !    !$OMP END PARALLEL 
    !
    SUM_REAL=0.0d0 
    do ikb3=1,nkb3
     do ikb2=1,nkb2
      do ikb1=1,nkb1
       do ib=1,NTB 
        SUM_REAL=SUM_REAL+xow(ib,ikb1,ikb2,ikb3)  
       enddo 
      enddo 
     enddo 
    enddo 
    write(6,'(a40)')'+++ m_dmx: make_xow +++' 
    write(6,'(a40,f15.8)')'(2/N)sum_{a,k}f(a,k)',2.0d0*SUM_REAL/dble(NTK) !2 is spin   
    write(6,*) 
    return 
  end subroutine make_xow 

  subroutine make_index_kpt(NTK,nkb1,nkb2,nkb3,SK0,index_kpt)      
    implicit none 
    integer,intent(in)::NTK,nkb1,nkb2,nkb3
    real(8),intent(in)::SK0(3,NTK)  
    integer,intent(out)::index_kpt(nkb1,nkb2,nkb3)    
    integer::ik,ix,iy,iz
    real(8)::x,y,z
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

  subroutine make_WEIGHT_R(nkb1,nkb2,nkb3,Na1,Na2,Na3,NTK,flg_weight,WEIGHT_R) 
    implicit none 
    integer,intent(in)::nkb1,nkb2,nkb3,Na1,Na2,Na3,NTK,flg_weight  
    real(8),intent(out)::WEIGHT_R(-Na1:Na1,-Na2:Na2,-Na3:Na3)
    integer::ia1,ia2,ia3
    real(8)::SUM_REAL 
    real(8),parameter::dlt_eps=1.0d-6 
    !
    if(flg_weight.eq.1)then
     write(6,'(a40)')'WEIGHT CALCULATED' 
     SUM_REAL=0.0d0 
     do ia1=-Na1,Na1
      do ia2=-Na2,Na2
       do ia3=-Na3,Na3
        if(abs(ia1)==Na1.and.mod(nkb1,2)==0.and.Na1/=0) WEIGHT_R(ia1,ia2,ia3)=WEIGHT_R(ia1,ia2,ia3)*0.5d0 
        if(abs(ia2)==Na2.and.mod(nkb2,2)==0.and.Na2/=0) WEIGHT_R(ia1,ia2,ia3)=WEIGHT_R(ia1,ia2,ia3)*0.5d0 
        if(abs(ia3)==Na3.and.mod(nkb3,2)==0.and.Na3/=0) WEIGHT_R(ia1,ia2,ia3)=WEIGHT_R(ia1,ia2,ia3)*0.5d0 
        SUM_REAL=SUM_REAL+WEIGHT_R(ia1,ia2,ia3)
       enddo!ia3
      enddo!ia2
     enddo!ia1
     write(6,'(a40,f15.8,i8)')'SUM_WEIGHT,NTK',SUM_REAL,NTK  
     if(abs(SUM_REAL-dble(NTK))>dlt_eps)then 
      stop 'SUM_WEIGHT/=NTK'
     endif 
    else
     write(6,'(a40)')'WEIGHT NOT CALCULATED' 
    endif 
    !
    return 
  end subroutine make_WEIGHT_R 

  subroutine make_pf(Na1,Na2,Na3,nkb1,nkb2,nkb3,NTK,a1,a2,a3,SK0,pf) 
    implicit none 
    integer,intent(in)::Na1,Na2,Na3,nkb1,nkb2,nkb3,NTK 
    real(8),intent(in)::a1(3),a2(3),a3(3) 
    real(8),intent(in)::SK0(3,NTK) 
    complex(8),intent(out)::pf(-Na1:Na1,-Na2:Na2,-Na3:Na3,NTK)   
    integer::ia1,ia2,ia3,ik 
    integer::ia1min,ia2min,ia3min 
    real(8)::PHASE 
    complex(8),parameter::ci=(0.0D0,1.0D0) 
    real(8),parameter::tpi=2.0d0*dacos(-1.0d0)
    !
    do ik=1,NTK 
     do ia3=-Na3,Na3 
      do ia2=-Na2,Na2 
       do ia1=-Na1,Na1 
        call search_Rmin(ia1,ia2,ia3,nkb1,nkb2,nkb3,a1(1),a2(1),a3(1),ia1min,ia2min,ia3min)
        PHASE=tpi*(SK0(1,ik)*DBLE(ia1min)+SK0(2,ik)*DBLE(ia2min)+SK0(3,ik)*DBLE(ia3min))           
        pf(ia1,ia2,ia3,ik)=EXP(-ci*PHASE) !notice minus sign
       enddo!ia1 
      enddo!ia2 
     enddo!ia3 
    enddo!ik  
    !
    return 
  end subroutine make_pf 

  subroutine search_Rmin(i,j,k,nkb1,nkb2,nkb3,a1,a2,a3,imin,jmin,kmin)
    implicit none
    integer,intent(in)::i,j,k,nkb1,nkb2,nkb3
    real(8),intent(in)::a1(3),a2(3),a3(3) 
    integer,intent(out)::imin,jmin,kmin
    integer::nmin,mmin,lmin
    integer::n,m,l 
    real(8)::R_pos(3),R_abs,R_min,R_bfr
    R_pos(:)=dble(i)*a1(:)+dble(j)*a2(:)+dble(k)*a3(:)
    nmin=0;mmin=0;lmin=0
    R_bfr=dsqrt(R_pos(1)**2+R_pos(2)**2+R_pos(3)**2) 
    R_min=R_bfr 
    do n=-3,3
     do m=-3,3
      do l=-3,3
       R_pos(:)=dble(i+n*nkb1)*a1(:)+dble(j+m*nkb2)*a2(:)+dble(k+l*nkb3)*a3(:)
       R_abs=dsqrt(R_pos(1)**2+R_pos(2)**2+R_pos(3)**2) 
       if(R_min>R_abs)then 
        R_min=R_abs
        nmin=n
        mmin=m
        lmin=l
       endif 
      enddo
     enddo
    enddo
    imin=i+nmin*nkb1
    jmin=j+mmin*nkb2
    kmin=k+lmin*nkb3
    return
  end subroutine search_Rmin
  
  subroutine make_dmx(Na1,Na2,Na3,NR,NTB,NTK,NWF,nkb1,nkb2,nkb3,lat_num_a1,lat_num_a2,lat_num_a3,SK0,UNT,xow,pf,WEIGHT_R,dmx) 
    implicit none
    integer,intent(in)::Na1,Na2,Na3,NR,NTB,NTK,NWF,nkb1,nkb2,nkb3   
    integer,intent(in)::lat_num_a1(NR),lat_num_a2(NR),lat_num_a3(NR)
    real(8),intent(in)::SK0(3,NTK) 
    real(8),intent(in)::WEIGHT_R(-Na1:Na1,-Na2:Na2,-Na3:Na3)
    real(8),intent(in)::xow(NTB,nkb1,nkb2,nkb3) 
    complex(8),intent(in)::UNT(NTB,NWF,NTK) 
    complex(8),intent(in)::pf(-Na1:Na1,-Na2:Na2,-Na3:Na3,NTK) 
    complex(8),intent(out)::dmx(NWF,NWF,-Na1:Na1,-Na2:Na2,-Na3:Na3) 
    integer::iR,ia1,ia2,ia3,iw,jw,ib,ikb1,ikb2,ikb3,ik 
    integer::index_kpt(nkb1,nkb2,nkb3) 
    real(8)::DMX_REAL 
    complex(8)::SUM_CMPX,DMX_CMPX
    !
    index_kpt=0 
    call make_index_kpt(NTK,nkb1,nkb2,nkb3,SK0(1,1),index_kpt(1,1,1)) 
    !
    !write(6,*) 
    !write(6,'(a40)')'+++ m_dmx: UNT check +++' 
    !do iw=1,NWF
    ! SUM_CMPX=0.0d0 
    ! do ib=1,NTB 
    !  do ikb1=1,nkb1
    !   do ikb2=1,nkb2
    !    do ikb3=1,nkb3 
    !     ik=index_kpt(ikb1,ikb2,ikb3)
    !     write(6,'(a,i8,2f15.8)')'pf000ik',ik,pf(0,0,0,ik)  
    !     SUM_CMPX=SUM_CMPX+CONJG(UNT(ib,iw,ik))*UNT(ib,iw,ik)*pf(0,0,0,ik)*WEIGHT_R(0,0,0)*xow(ib,ikb1,ikb2,ikb3) 
    !    enddo
    !   enddo
    !  enddo
    ! enddo 
    ! write(6,'(a40,i8,x,2f15.10)')'(2/N)sum_{a,k}|<ak|i0>|^2',iw,2.0d0*SUM_CMPX/dble(NTK)  
    !enddo 
    !
    dmx=0.0d0 
    do iR=1,NR 
     ia1=lat_num_a1(iR)
     ia2=lat_num_a2(iR)
     ia3=lat_num_a3(iR)
     !write(6,'(3i5,f15.8)')ia1,ia2,ia3,WEIGHT_R(ia1,ia2,ia3) 
     do iw=1,NWF
      do jw=1,NWF 
       SUM_CMPX=0.0d0 
       do ib=1,NTB 
        do ikb1=1,nkb1
         do ikb2=1,nkb2
          do ikb3=1,nkb3 
           ik=index_kpt(ikb1,ikb2,ikb3)
           SUM_CMPX&
          =SUM_CMPX&
          +CONJG(UNT(ib,iw,ik))*UNT(ib,jw,ik)*xow(ib,ikb1,ikb2,ikb3)*pf(ia1,ia2,ia3,ik)*WEIGHT_R(ia1,ia2,ia3)  
          enddo!ikb1
         enddo!ikb2
        enddo!ikb3
       enddo!ib  
       dmx(iw,jw,ia1,ia2,ia3)=2.0d0*SUM_CMPX/dble(NTK) !2 is spin   
      enddo!jw 
     enddo!iw 
    enddo!iR 
    ! 
    !impose sum rule
    !
    SUM_CMPX=0.0d0 
    do iw=1,NWF
     SUM_CMPX=SUM_CMPX+dmx(iw,iw,0,0,0) 
    enddo 
    write(6,*) 
    write(6,'(a40)')'+++ m_dmx: Tr(DMX) +++' 
    DMX_REAL=dble(nint(dble(SUM_CMPX))) 
    DMX_CMPX=cmplx(DMX_REAL,0.0d0)    
    write(6,'(a40,2f15.10)')'before correction',SUM_CMPX 
    write(6,'(a40,2f15.10)')'after  correction',DMX_CMPX 
    do iw=1,NWF
     dmx(iw,iw,0,0,0)=dmx(iw,iw,0,0,0)*(DMX_CMPX/SUM_CMPX) 
    enddo 
    !
    return 
  end subroutine make_dmx 

  subroutine wrt_model_dr(NWF,Na1,Na2,Na3,DR) 
    implicit none 
    integer,intent(in)::NWF,Na1,Na2,Na3 
    complex(8),intent(in)::DR(NWF,NWF,-Na1:Na1,-Na2:Na2,-Na3:Na3)
    integer::N_element 
    integer::ia1,ia2,ia3,ib,jb,i 
    integer,allocatable::unit_vec(:)!unit_vec(N_element) 
    real(8),parameter::au=27.21151d0
    !
    N_element=(2*Na1+1)*(2*Na2+1)*(2*Na3+1)  
    !
    !OPEN(307,W,FILE='zvo_dr.dat') 
    !
    OPEN(307,FILE='./dir-model/zvo_dr.dat') 
    write(307,'(a)')'wannier90 format for vmcdry.out or HPhi -sdry'
    write(307,'(i10)') NWF 
    write(307,'(i10)') N_element
    !
    allocate(unit_vec(N_element)); unit_vec=1
    write(307,'(15i5)')(unit_vec(i),i=1,N_element) 
    deallocate(unit_vec) 
    !
    do ia1=-Na1,Na1
     do ia2=-Na2,Na2
      do ia3=-Na3,Na3 
       do ib=1,NWF 
        do jb=1,NWF 
         write(307,'(i5,i5,i5,i5,i5,2f15.10)') ia1,ia2,ia3,ib,jb,DR(ib,jb,ia1,ia2,ia3) 
        enddo!jb
       enddo!ib
      enddo!ia3
     enddo!ia2
    enddo!ia1 
    return
    !
  end subroutine wrt_model_dr 
  !
end module m_dmx 
