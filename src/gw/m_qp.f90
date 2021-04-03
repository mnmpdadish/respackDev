module m_qp  
  implicit none
contains
  
  subroutine calculate_quasi_particle_band(& 
    Nk_irr,NTK,Mb,NTB,NWF,La1,La2,La3,nsgm,nsgmqp,FermiEnergy,shift_value,&
    numirr,numMK,Nb,Ns,E_EIGI,sgmw,sgmwqp,& 
    SCirr,SXirr,VXCirr,UNT,& 
    NSK_BAND_DISP,nkb1,nkb2,nkb3,Ndiv,N_sym_points,& 
    a1,a2,a3,b1,b2,b3,SK0,SK_BAND_DISP)  
    !
    use m_FR, only: calculate_FR_from_FKdiag 
    use m_eigenstate, only: calculate_qpvalue
    use m_band, only: calculate_banddisp 
    !
    implicit none 
    !
    !in
    !
    integer,intent(in)::Nk_irr,NTK,Mb,NTB,NWF,La1,La2,La3,nsgm,nsgmqp 
    integer,intent(in)::numirr(NTK),numMK(Nk_irr),Nb(NTK),Ns(NTK)
    real(8),intent(in)::FermiEnergy,shift_value
    real(8),intent(in)::E_EIGI(NTB,NK_irr),sgmw(nsgm),sgmwqp(nsgmqp) 
    complex(4),intent(in)::SCirr(nsgm,Mb,Mb,Nk_irr)
    complex(8),intent(in)::SXirr(Mb,Mb,Nk_irr)
    complex(8),intent(in)::VXCirr(Mb,Mb,Nk_irr)
    complex(8),intent(in)::UNT(Mb,NWF,NTK) 
    integer,intent(in)::NSK_BAND_DISP,nkb1,nkb2,nkb3,Ndiv,N_sym_points 
    real(8),intent(in)::a1(3),a2(3),a3(3),b1(3),b2(3),b3(3) 
    real(8),intent(in)::SK0(3,NTK),SK_BAND_DISP(3,NSK_BAND_DISP) 
    !
    !local 
    !
    complex(8),allocatable::eqpMK(:,:)!eqpMK(Mb,NTK) 
    complex(8),allocatable::eqpr(:,:,:,:,:)!eqpr(NWF,NWF,-La1:La1,-La2:La2,-La3:La3)
    complex(8),allocatable::EQP(:,:)!EQP(NWF,ncalck) 
    real(8),allocatable::ReEQP(:,:)!ReEQP(NWF,ncalck) 
    integer::flg_weight
    real(8)::threshold_e 
    real(8),parameter::au=27.21151d0
    !
    !1. set Eqp(k) 
    !
    allocate(eqpMK(Mb,NTK));eqpMK(:,:)=0.0d0  
    call calculate_eqp(Nk_irr,NTK,Mb,NTB,NWF,nsgm,nsgmqp,FermiEnergy,shift_value,numirr(1),numMK(1),Nb(1),Ns(1),E_EIGI(1,1),&
    sgmw(1),sgmwqp(1),SCirr(1,1,1,1),SXirr(1,1,1),VXCirr(1,1,1),eqpMK(1,1),SK0(1,1))  
    !
    !2. set Eqp(R) m_FR
    !
    allocate(eqpr(NWF,NWF,-La1:La1,-La2:La2,-La3:La3));eqpr(:,:,:,:,:)=0.0d0  
    call calculate_FR_from_FKdiag(Mb,NWF,NTK,La1,La2,La3,Nb(1),SK0(1,1),UNT(1,1,1),eqpMK(1,1),eqpr(1,1,-La1,-La2,-La3)) 
    deallocate(eqpMK)  
    !
    !3. make Eqp(k') m_eigenstate 
    !
    allocate(EQP(NWF,NSK_BAND_DISP));EQP=0.0d0 
    flg_weight=1 !Flg whether calculate weighted transfers (0:not calc, 1:calc)
    call calculate_qpvalue(NWF,NSK_BAND_DISP,La1,La2,La3,nkb1,nkb2,nkb3,flg_weight,a1(1),a2(1),a3(1),SK_BAND_DISP(1,1),eqpr(1,1,-La1,-La2,-La3),EQP(1,1)) 
    deallocate(eqpr) 
    !
    !4. calc band disp m_band 
    !  
    threshold_e=0.0d0 !Variable for transfer cutoff 
    allocate(ReEQP(NWF,NSK_BAND_DISP));ReEQP=0.0d0 
    ReEQP(:,:)=dble(EQP(:,:)) 
    call calculate_banddisp(NWF,NSK_BAND_DISP,Ndiv,N_sym_points,threshold_e,b1(1),b2(1),b3(1),SK_BAND_DISP(1,1),ReEQP(1,1)) 
    deallocate(EQP,ReEQP) 
    return
  end subroutine calculate_quasi_particle_band       
  
  subroutine calculate_eqp(Nk_irr,NTK,Mb,NTB,NWF,nsgm,nsgmqp,FermiEnergy,shift_value,numirr,numMK,Nb,Ns,E_EIGI,&
  sgmw,sgmwqp,SCirr,SXirr,VXCirr,eqpMK,SK0)
    !
    use bspline
    !
    implicit none 
    !
    !in
    !    
    integer,intent(in)::Nk_irr,NTK,Mb,NTB,NWF,nsgm,nsgmqp  
    integer,intent(in)::numirr(NTK),numMK(Nk_irr),Nb(NTK),Ns(NTK)
    real(8),intent(in)::FermiEnergy,shift_value 
    real(8),intent(in)::E_EIGI(NTB,Nk_irr),sgmw(nsgm) 
    real(8),intent(in)::SK0(3,NTK) 
    complex(4),intent(in)::SCirr(nsgm,Mb,Mb,Nk_irr)
    complex(8),intent(in)::SXirr(Mb,Mb,Nk_irr)
    complex(8),intent(in)::VXCirr(Mb,Mb,Nk_irr)
    !
    !out
    !
    complex(8),intent(out)::eqpMK(Mb,NTK) 
    !
    !local 
    !
    integer::ie,ik,jb,je,ia1,ia2,ia3,iw,jw,ikir 
    real(8)::dscdw,zak,dlt,dlt_min,en  
    complex(8)::dsc
    complex(8)::SUM_CMPX 
    complex(8),allocatable::eqpirr(:,:)!eqp(Mb,Nk_irr)  
    !
    !spline 
    !
    integer,parameter::irnk=2  
    integer::nknot
    integer::ierr,igz
    real(8),intent(in)::sgmwqp(nsgmqp)!sgmwqp(nknot+irnk) 
    real(8)::x(nsgmqp),c(nsgmqp)!x(nknot+irnk),c(nknot+irnk)
    real(8)::xx(nsgm),yy(nsgm),s(nsgm)
    real(8)::gzai(1-irnk:nsgmqp)!gzai(1-irnk:nknot+irnk) 
    real(8)::rn(0:irnk)
    real(8)::rnz
    !
    !parameter 
    !
    real(8),parameter::au=27.21151d0
    !
    nknot=nsgmqp-irnk 
    !
    do ie=1,nsgmqp!nknot+irnk
     x(ie)=sgmwqp(ie) 
    enddo!ie
    call formgzai(nknot+irnk,irnk,x,nknot,gzai,ierr)
    !
    allocate(eqpirr(Mb,Nk_irr));eqpirr(:,:)=0.0d0 
    !
    write(6,*) 
    write(6,'(a50)')'#omega, Z(ak), Re[Z(ak)*DSC], Im[Z(ak)*DSC]'
    write(6,*) 
    do ikir=1,Nk_irr 
     ik=numMK(ikir) 
     do jb=1,Nb(ik)!Mb!NWF 
      do ie=1,nsgm 
       xx(ie)=sgmw(ie) 
       yy(ie)=dble(SCirr(ie,jb,jb,ikir))   
      enddo!ie  
      do ie=1,nsgm 
       call bspl(xx(ie),nknot,irnk,gzai,irnk,igz,rn)
      enddo!ie
      call bsmooth(nsgm,xx,yy,nknot,gzai,irnk,c,ierr)
      call evspline(nknot,gzai,irnk,0,c,xx,nsgm,s)
      !---
      !write(6,*)'# SK=',SK0(1,ik),SK0(2,ik),SK0(3,ik) 
      !write(6,*)'# jb=',jb 
      !do ie=1,nsgm 
      ! write(6,'(3f12.6)') xx(ie)*au,yy(ie)*au,s(ie)*au
      !end do
      !write(6,*)'# diff point=',E_EIGI(jb+Ns(ik),ikir)*au  
      !stop
      !---
      dlt_min=1.0d0 
      do ie=1,nsgm 
       dlt=dabs(sgmw(ie)-E_EIGI(jb+Ns(ik),ikir)) 
       if(dlt<dlt_min)then 
        dlt_min=dlt 
        je=ie 
       endif 
      enddo!ie 
      dscdw=(s(je+1)-s(je-1))/(sgmw(je+1)-sgmw(je-1)) 
      zak=1.0d0/(1.0d0-dscdw) 
      dsc=-VXCirr(jb,jb,ikir)-SXirr(jb,jb,ikir)+SCirr(je,jb,jb,ikir)   
      !
      !20191205
      !
      eqpirr(jb,ikir)=E_EIGI(jb+Ns(ik),ikir)+(dsc-shift_value)*zak  
      !eqpirr(jb,ikir)=E_EIGI(jb+Ns(ik),ikir)+dsc*zak-shift_value  
      !
      en=E_EIGI(jb+Ns(ik),ikir)-FermiEnergy 
      write(6,'(4F20.10)') en*au, zak, zak*dsc*au    
      !--
      if(zak>1.0d0)then 
       do ie=1,nsgm 
        write(6,'(3f12.6)')(xx(ie)-FermiEnergy)*au,yy(ie)*au,s(ie)*au
       enddo!ie 
       !stop
      endif  
      !--
      !do ie=1,nsgm 
      ! SC(ie,jb,ik)=cmplx(s(ie),dimag(SC(ie,jb,ik)))      
      !enddo 
      !do ie=1,nsgm  
      ! en=wini+wstep*dble(ie-1)
      ! write(6,'(3F15.10)') en*au,SC(ie,jb,ik)*au 
      !enddo 
      !stop
      !--
     enddo!jb  
    enddo!ikir   
    !
    eqpMK(:,:)=0.0d0 
    do ik=1,NTK 
     ikir=numirr(ik) 
     eqpMK(:,ik)=eqpirr(:,ikir) 
    enddo!ik 
    deallocate(eqpirr) 
    !
    return 
  end subroutine calculate_eqp 
  
end module m_qp 
