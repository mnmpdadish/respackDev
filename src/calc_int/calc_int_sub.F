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
!--
      SKQ(1)=SK0(1,ik)+q1
      SKQ(2)=SK0(2,ik)+q2 
      SKQ(3)=SK0(3,ik)+q3 
!--
!     if(SKQ(1)<-0.5D0)then 
!      SKQ(1)=SKQ(1)+1.0D0 
!      shift_G(1)=-1
!      endif 
!      if(SKQ(2)<-0.5D0)then 
!      SKQ(2)=SKQ(2)+1.0D0 
!      shift_G(2)=-1
!      endif 
!      if(SKQ(3)<-0.5D0)then 
!      SKQ(3)=SKQ(3)+1.0D0 
!      shift_G(3)=-1
!      endif 
!      if(SKQ(1)>0.5D0)then 
!      SKQ(1)=SKQ(1)-1.0D0 
!      shift_G(1)=+1
!      endif 
!      if(SKQ(2)>0.5D0)then 
!      SKQ(2)=SKQ(2)-1.0D0 
!      shift_G(2)=+1
!      endif 
!      if(SKQ(3)>0.5D0)then 
!      SKQ(3)=SKQ(3)-1.0D0 
!      shift_G(3)=+1
!      endif 
!--
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
!--
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
!--
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
!--
      do jk=1,NTK
       if(ABS(SK0(1,jk)-SKQ(1))<1.D-6.and.
     +    ABS(SK0(2,jk)-SKQ(2))<1.D-6.and. 
     +    ABS(SK0(3,jk)-SKQ(3))<1.D-6)then 
        ikq=jk 
       endif 
      enddo 
      RETURN  
      END 
!--
      subroutine calc_InterStateMatrix(NTK,NTG,NG0,KG0,C0_K,C0_KQ,
     +                                 ik,ikq,nwx2,nwy2,nwz2,
     +                                 nfft1,nfft2,Nl123,
     +                                 wfunc,fftwk,fs,LG0,NG_for_eps,
     +                                 shift_G,m_tmp)
      use fft_3d 
      implicit none 
      type(fft3_struct)::fs 
      integer::NTK,NTG 
      integer::NG0(NTK),KG0(3,NTG,NTK) 
      complex(8)::C0_K(NTG),C0_KQ(NTG) 
      integer::ik,ikq 
      integer::nwx2,nwy2,nwz2,nfft1,nfft2,Nl123 
      real(8)::wfunc(Nl123*2),fftwk(Nl123*2) 
      integer::NG_for_eps 
      integer::LG0(3,NTG) 
      integer::shift_G(3)
      integer::ig,igb1,igb2,igb3,ind,igL,igL1,igL2,igL3  
      integer::ir,ir1,ir2,ir3
      real(8)::dG1,dG2,dG3,phase 
      complex(8),allocatable::cell_periodic_k(:)
      complex(8),allocatable::cell_periodic_kq(:)
      complex(8)::f,pf
      complex(kind=8)::phx(nwx2),phy(nwy2),phz(nwz2)
      complex(kind=8)::m_tmp(NG_for_eps) 
      real(8),parameter::pi=dacos(-1.0d0)
      real(8),parameter::tpi=2.0d0*pi 
      complex(8),parameter::ci=(0.0D0,1.0D0) 
!--
      allocate(cell_periodic_k(Nl123))
      allocate(cell_periodic_kq(Nl123))
      wfunc(:)=0.0D0
      fftwk(:)=0.0D0
      do ig=1,NG0(ik) 
      igb1=KG0(1,ig,ik) 
      igb2=KG0(2,ig,ik) 
      igb3=KG0(3,ig,ik) 
      igb1=MOD(nwx2+igb1,nwx2)+1
      igb2=MOD(nwy2+igb2,nwy2)+1
      igb3=MOD(nwz2+igb3,nwz2)+1
      ind=igb1+(igb2-1)*nfft1+(igb3-1)*nfft1*nfft2 
      wfunc(ind)=dble(C0_K(ig))
      wfunc(ind+Nl123)=dimag(C0_K(ig))
      enddo 
      call fft3_bw(fs,wfunc,fftwk) 
      do ir=1,Nl123 
      cell_periodic_k(ir)=cmplx(wfunc(ir),wfunc(ir+Nl123)) 
      enddo 
!--
      wfunc(:)=0.0D0
      fftwk(:)=0.0D0
      do ig=1,NG0(ikq)
      igb1=KG0(1,ig,ikq)
      igb2=KG0(2,ig,ikq) 
      igb3=KG0(3,ig,ikq) 
      igb1=MOD(nwx2+igb1,nwx2)+1
      igb2=MOD(nwy2+igb2,nwy2)+1      
      igb3=MOD(nwz2+igb3,nwz2)+1 
      ind=igb1+(igb2-1)*nfft1+(igb3-1)*nfft1*nfft2 
      wfunc(ind)=dble(C0_KQ(ig))
      wfunc(ind+Nl123)=dimag(C0_KQ(ig)) 
      enddo 
      call fft3_bw(fs,wfunc,fftwk)    
      do ir=1,Nl123
      cell_periodic_kq(ir)=cmplx(wfunc(ir),wfunc(ir+Nl123)) 
      enddo 
!--
      dG1=dble(shift_G(1))
      dG2=dble(shift_G(2))
      dG3=dble(shift_G(3))
!YOSHIHIDE YOSHIMOTO 20080620 
      do ir1=1,nwx2
        phase = tpi*(dble(ir1-1)*dG1/dble(nwx2))
        phx(ir1) = exp(-ci*phase)
      end do
      do ir2=1,nwy2
        phase = tpi*(dble(ir2-1)*dG2/dble(nwy2))
        phy(ir2) = exp(-ci*phase)
      end do
      do ir3=1,nwz2
        phase = tpi*(dble(ir3-1)*dG3/dble(nwz2))
        phz(ir3) = exp(-ci*phase)
      end do
!--
      do ir3=1,nwz2
      do ir2=1,nwy2
      do ir1=1,nwx2
       ir=ir1+(ir2-1)*nfft1+(ir3-1)*nfft1*nfft2 
!YOSHIHIDE YOSHIMOTO 20080620 
         pf=phx(ir1)*phy(ir2)*phz(ir3)
         cell_periodic_kq(ir)
     + = cell_periodic_kq(ir) * pf 
      enddo 
      enddo 
      enddo 
!--
      wfunc(:)=0.0D0
      fftwk(:)=0.0D0
      do ir=1,Nl123 
      f=CONJG(cell_periodic_kq(ir))*cell_periodic_k(ir)
      wfunc(ir)=dble(f)
      wfunc(ir+Nl123)=dimag(f) 
      enddo 
      call fft3_fw(fs,wfunc,fftwk) 
      do igL=1,NG_for_eps
      igL1=-LG0(1,igL)
      igL2=-LG0(2,igL)
      igL3=-LG0(3,igL) 
      igL1=MOD(nwx2+igL1,nwx2)+1
      igL2=MOD(nwy2+igL2,nwy2)+1
      igL3=MOD(nwz2+igL3,nwz2)+1
      ind=igL1+(igL2-1)*nfft1+(igL3-1)*nfft1*nfft2
      m_tmp(igL)=cmplx(wfunc(ind),wfunc(ind+Nl123))
      enddo 
!--
      deallocate(cell_periodic_k)
      deallocate(cell_periodic_kq)
      RETURN  
      END 
!--
      subroutine make_eps(NTG,NTGQ,ne,itrs,NG,LGtmp,RWtmp,rginvtmp,
     + pgtmp,nnp,L1,L2,L3,packtmp,epsirr,epsmk) 
      implicit none 
      integer::NTG,NTGQ,itrs,NG,ne,L1,L2,L3,nnp 
      integer::LGtmp(3,NTG) 
      integer::RWtmp(3) 
      real(8)::rginvtmp(3,3) 
      integer::pgtmp(3) 
      integer::packtmp(-L1:L1,-L2:L2,-L3:L3) 
      complex(4)::epsirr(NTGQ,NTGQ,ne) 
      integer::ig,jg,i1,i2,i3,j1,j2,j3,k1,k2,k3,iig,jjg 
      real(8)::phase 
      complex(8)::pf1,pf2 
      complex(4)::epsmk(NTGQ,NTGQ,ne) 
      real(8),parameter::pi=dacos(-1.0d0)
      real(8),parameter::tpi=2.0d0*pi 
      complex(8),parameter::ci=(0.0D0,1.0D0) 
!--
      epsmk(:,:,:)=0.0d0 
      select case(itrs) 
      case(1)!=== not time-reversal ===*      
      do ig=1,NG 
       i1=LGtmp(1,ig); j1=LGtmp(2,ig); k1=LGtmp(3,ig) 
       i2=i1+RWtmp(1); j2=j1+RWtmp(2); k2=k1+RWtmp(3) 
       i3=int(rginvtmp(1,1))*i2+int(rginvtmp(1,2))*j2 
     +   +int(rginvtmp(1,3))*k2 
       j3=int(rginvtmp(2,1))*i2+int(rginvtmp(2,2))*j2 
     +   +int(rginvtmp(2,3))*k2 
       k3=int(rginvtmp(3,1))*i2+int(rginvtmp(3,2))*j2 
     +   +int(rginvtmp(3,3))*k2 
       iig=packtmp(i3,j3,k3) 
       phase=tpi*(dble(i1)*dble(pgtmp(1))
     +           +dble(j1)*dble(pgtmp(2))
     +           +dble(k1)*dble(pgtmp(3))) 
       pf1=exp(-ci*phase/dble(nnp)) 
       do jg=1,NG 
        i1=LGtmp(1,jg); j1=LGtmp(2,jg); k1=LGtmp(3,jg) 
        i2=i1+RWtmp(1); j2=j1+RWtmp(2); k2=k1+RWtmp(3) 
        i3=int(rginvtmp(1,1))*i2+int(rginvtmp(1,2))*j2 
     +    +int(rginvtmp(1,3))*k2 
        j3=int(rginvtmp(2,1))*i2+int(rginvtmp(2,2))*j2 
     +    +int(rginvtmp(2,3))*k2 
        k3=int(rginvtmp(3,1))*i2+int(rginvtmp(3,2))*j2 
     +    +int(rginvtmp(3,3))*k2 
        jjg=packtmp(i3,j3,k3) 
        phase=tpi*(dble(i1)*dble(pgtmp(1))
     +            +dble(j1)*dble(pgtmp(2))
     +            +dble(k1)*dble(pgtmp(3))) 
        pf2=exp(ci*phase/dble(nnp)) 
!       write(6,'(a,x,2i5)') 'iig jjg',iig,jjg  
        epsmk(ig,jg,:)=epsirr(iig,jjg,:)*pf1*pf2 
       enddo!jg 
      enddo!ig 
      case(-1)!=== time-reversal ===*      
      do ig=1,NG 
       i1=-LGtmp(1,ig); j1=-LGtmp(2,ig); k1=-LGtmp(3,ig) 
       i2=i1+RWtmp(1); j2=j1+RWtmp(2); k2=k1+RWtmp(3) 
       i3=int(rginvtmp(1,1))*i2+int(rginvtmp(1,2))*j2 
     +   +int(rginvtmp(1,3))*k2 
       j3=int(rginvtmp(2,1))*i2+int(rginvtmp(2,2))*j2 
     +   +int(rginvtmp(2,3))*k2 
       k3=int(rginvtmp(3,1))*i2+int(rginvtmp(3,2))*j2 
     +   +int(rginvtmp(3,3))*k2 
       iig=packtmp(i3,j3,k3) 
       phase=tpi*(dble(i1)*dble(pgtmp(1))
     +           +dble(j1)*dble(pgtmp(2))
     +           +dble(k1)*dble(pgtmp(3))) 
       pf1=exp(-ci*phase/dble(nnp)) 
       do jg=1,NG 
        i1=-LGtmp(1,jg); j1=-LGtmp(2,jg); k1=-LGtmp(3,jg) 
        i2=i1+RWtmp(1); j2=j1+RWtmp(2); k2=k1+RWtmp(3) 
        i3=int(rginvtmp(1,1))*i2+int(rginvtmp(1,2))*j2 
     +    +int(rginvtmp(1,3))*k2 
        j3=int(rginvtmp(2,1))*i2+int(rginvtmp(2,2))*j2 
     +    +int(rginvtmp(2,3))*k2 
        k3=int(rginvtmp(3,1))*i2+int(rginvtmp(3,2))*j2 
     +    +int(rginvtmp(3,3))*k2 
        jjg=packtmp(i3,j3,k3) 
        phase=tpi*(dble(i1)*dble(pgtmp(1))
     +            +dble(j1)*dble(pgtmp(2))
     +            +dble(k1)*dble(pgtmp(3))) 
        pf2=exp(ci*phase/dble(nnp)) 
        epsmk(ig,jg,:)=epsirr(jjg,iig,:)*pf1*pf2 
       enddo!jg 
      enddo!ig 
      end select 
!     write(6,*)'finish make_eps'
!--
      return 
      end 
!--
      SUBROUTINE OUTER_PRODUCT(vec_x,vec_y,vec_z)
      implicit none 
      real(8)::vec_x(3),vec_y(3),vec_z(3) 
      vec_z(1)=vec_x(2)*vec_y(3)-vec_x(3)*vec_y(2)
      vec_z(2)=vec_x(3)*vec_y(1)-vec_x(1)*vec_y(3) 
      vec_z(3)=vec_x(1)*vec_y(2)-vec_x(2)*vec_y(1)
      RETURN
      END
!--
      subroutine search_q(SQI,SQ)
      implicit none 
      real(8)::SQI(3),SQ(3)
      real(8),parameter::dlt_BZ=1.0d-6!20170322  
      SQ=SQI 
      if(SQ(1)>  1.50d0+dlt_BZ) SQ(1)=SQ(1)-2.0d0 
      if(SQ(1)>  0.50d0+dlt_BZ) SQ(1)=SQ(1)-1.0d0 
      if(SQ(1)<=-1.50d0+dlt_BZ) SQ(1)=SQ(1)+2.0d0 
      if(SQ(1)<=-0.50d0+dlt_BZ) SQ(1)=SQ(1)+1.0d0 
      if(SQ(2)>  1.50d0+dlt_BZ) SQ(2)=SQ(2)-2.0d0 
      if(SQ(2)>  0.50d0+dlt_BZ) SQ(2)=SQ(2)-1.0d0 
      if(SQ(2)<=-1.50d0+dlt_BZ) SQ(2)=SQ(2)+2.0d0 
      if(SQ(2)<=-0.50d0+dlt_BZ) SQ(2)=SQ(2)+1.0d0 
      if(SQ(3)>  1.50d0+dlt_BZ) SQ(3)=SQ(3)-2.0d0 
      if(SQ(3)>  0.50d0+dlt_BZ) SQ(3)=SQ(3)-1.0d0 
      if(SQ(3)<=-1.50d0+dlt_BZ) SQ(3)=SQ(3)+2.0d0 
      if(SQ(3)<=-0.50d0+dlt_BZ) SQ(3)=SQ(3)+1.0d0 
      RETURN  
      END 
!--
      subroutine kcheck(ktmp,RWtmp) 
      implicit none 
      real(8),intent(inout)::ktmp(3)
      integer,intent(out)::RWtmp(3) 
      real(8),parameter::dlt_BZ=1.0d-6!20170321 
!--
!      write(6,'(a8,3f15.10)')'before',ktmp(1),ktmp(2),ktmp(3)
!--
      if(ktmp(1)>1.50d0+dlt_BZ)then 
       ktmp(1)=ktmp(1)-2.0d0;RWtmp(1)=-2
      endif 
      if(ktmp(1)>0.50d0+dlt_BZ)then 
       ktmp(1)=ktmp(1)-1.0d0;RWtmp(1)=-1
      endif 
      if(ktmp(1)<=-1.50d0+dlt_BZ)then 
       ktmp(1)=ktmp(1)+2.0d0;RWtmp(1)=2 
      endif 
      if(ktmp(1)<=-0.50d0+dlt_BZ)then 
       ktmp(1)=ktmp(1)+1.0d0;RWtmp(1)=1
      endif 
!---
      if(ktmp(2)>1.50d0+dlt_BZ)then 
       ktmp(2)=ktmp(2)-2.0d0;RWtmp(2)=-2 
      endif 
      if(ktmp(2)>0.50d0+dlt_BZ)then 
       ktmp(2)=ktmp(2)-1.0d0;RWtmp(2)=-1
      endif 
      if(ktmp(2)<=-1.50d0+dlt_BZ)then 
       ktmp(2)=ktmp(2)+2.0d0;RWtmp(2)=2 
      endif 
      if(ktmp(2)<=-0.50d0+dlt_BZ)then 
       ktmp(2)=ktmp(2)+1.0d0;RWtmp(2)=1
      endif 
!--
      if(ktmp(3)>1.50d0+dlt_BZ)then 
       ktmp(3)=ktmp(3)-2.0d0;RWtmp(3)=-2 
      endif 
      if(ktmp(3)>0.50d0+dlt_BZ)then 
       ktmp(3)=ktmp(3)-1.0d0;RWtmp(3)=-1
      endif 
      if(ktmp(3)<=-1.50d0+dlt_BZ)then 
       ktmp(3)=ktmp(3)+2.0d0;RWtmp(3)=2 
      endif 
      if(ktmp(3)<=-0.50d0+dlt_BZ)then 
       ktmp(3)=ktmp(3)+1.0d0;RWtmp(3)=1
      endif 
!--
!     write(6,'(a8,3f15.10)')'after',ktmp(1),ktmp(2),ktmp(3)
!--
      return
      end 
!--
!20170322  
      subroutine kcheck_trs(ktmp,RWtmp)!20170322 
      implicit none 
      real(8),intent(inout)::ktmp(3)
      integer,intent(out)::RWtmp(3) 
      real(8),parameter::dlt_BZ=-1.0d-6 
!--
!      write(6,'(a8,3f15.10)') 'before',ktmp(1),ktmp(2),ktmp(3)
!--
      if(ktmp(1)>=1.50d0+dlt_BZ)then 
       ktmp(1)=ktmp(1)-2.0d0;RWtmp(1)=-2
      endif 
      if(ktmp(1)>=0.50d0+dlt_BZ)then 
       ktmp(1)=ktmp(1)-1.0d0;RWtmp(1)=-1
      endif 
      if(ktmp(1)<-1.50d0+dlt_BZ)then 
       ktmp(1)=ktmp(1)+2.0d0;RWtmp(1)=2 
      endif 
      if(ktmp(1)<-0.50d0+dlt_BZ)then 
       ktmp(1)=ktmp(1)+1.0d0;RWtmp(1)=1
      endif 
!--
      if(ktmp(2)>=1.50d0+dlt_BZ)then 
       ktmp(2)=ktmp(2)-2.0d0;RWtmp(2)=-2 
      endif 
      if(ktmp(2)>=0.50d0+dlt_BZ)then 
       ktmp(2)=ktmp(2)-1.0d0;RWtmp(2)=-1
      endif 
      if(ktmp(2)<-1.50d0+dlt_BZ)then 
       ktmp(2)=ktmp(2)+2.0d0;RWtmp(2)=2 
      endif 
      if(ktmp(2)<-0.50d0+dlt_BZ)then 
       ktmp(2)=ktmp(2)+1.0d0;RWtmp(2)=1
      endif 
!--
      if(ktmp(3)>=1.50d0+dlt_BZ)then 
       ktmp(3)=ktmp(3)-2.0d0;RWtmp(3)=-2 
      endif 
      if(ktmp(3)>=0.50d0+dlt_BZ)then 
       ktmp(3)=ktmp(3)-1.0d0;RWtmp(3)=-1
      endif 
      if(ktmp(3)<-1.50d0+dlt_BZ)then 
       ktmp(3)=ktmp(3)+2.0d0;RWtmp(3)=2 
      endif 
      if(ktmp(3)<-0.50d0+dlt_BZ)then 
       ktmp(3)=ktmp(3)+1.0d0;RWtmp(3)=1
      endif 
!--
!      write(6,'(a8,3f15.10)') 'after',ktmp(1),ktmp(2),ktmp(3)
!--
      return
      end 
!--
      subroutine make_LG0(NTG,b1,b2,b3,Gcut_for_eps,Gcut_for_psi, 
     + q1,q2,q3,LG0,NG_for_eps,NG_for_psi)
      implicit none 
      integer::NTG,igL,igL1,igL2,igL3
      integer::NG_for_eps,NG_for_psi            
      integer::LG0(3,NTG)    
      real(8)::b1(3),b2(3),b3(3) 
      real(8)::Gcut_for_eps,Gcut_for_psi  
      real(8)::q1,q2,q3,qgL2,qgL(3)
      integer,parameter::NGL1=150!100!300 
      integer,parameter::NGL2=150!100!300 
      integer,parameter::NGL3=150!100!300 
!--
      igL=0
      do igL1=-NGL1,NGL1 
      do igL2=-NGL2,NGL2 
      do igL3=-NGL3,NGL3 
      qgL(:)=(q1+dble(igL1))*b1(:)
     +      +(q2+dble(igL2))*b2(:)
     +      +(q3+dble(igL3))*b3(:)    
      qgL2=qgL(1)**2+qgL(2)**2+qgL(3)**2
      if(qgL2 <= Gcut_for_eps) then 
       igL=igL+1 
       LG0(1,igL)=igL1
       LG0(2,igL)=igL2
       LG0(3,igL)=igL3 
!      write(6,*) igL
      endif  
      enddo 
      enddo 
      enddo 
      NG_for_eps=igL 
!--
      do igL1=-NGL1,NGL1 
      do igL2=-NGL2,NGL2 
      do igL3=-NGL3,NGL3 
      qgL(:)=(q1+dble(igL1))*b1(:)
     +      +(q2+dble(igL2))*b2(:)
     +      +(q3+dble(igL3))*b3(:)    
      qgL2=qgL(1)**2+qgL(2)**2+qgL(3)**2
      if(qgL2>Gcut_for_eps.and.qgL2<=Gcut_for_psi) then 
       igL=igL+1 
       LG0(1,igL)=igL1
       LG0(2,igL)=igL2
       LG0(3,igL)=igL3 
!      write(6,*) igL
      endif  
      enddo 
      enddo 
      enddo 
      NG_for_psi=igL 
!--
!Do igL=1,NG_for_eps 
!write(6,*) igL,LG0(:,igL)
!enddo 
!--
      RETURN 
      END 
!--
      subroutine make_KG0(NTG,b1,b2,b3,Gcut,q1,q2,q3,KG0,NG) 
      implicit none 
      integer,intent(in)::NTG
      real(8),intent(in)::b1(3),b2(3),b3(3) 
      real(8),intent(in)::Gcut 
      real(8),intent(in)::q1,q2,q3        
      integer,intent(out)::KG0(3,NTG)    
      integer,intent(out)::NG 
      integer::igL,igL1,igL2,igL3
      real(8)::qgL(3),qgL2  
      integer,parameter::NGL1=150!100!300
      integer,parameter::NGL2=150!100!300 
      integer,parameter::NGL3=150!100!300
      igL=0
      do igL1=-NGL1,NGL1 
      do igL2=-NGL2,NGL2 
      do igL3=-NGL3,NGL3 
      qgL(:)=(q1+dble(igL1))*b1(:)+(q2+dble(igL2))*b2(:)
     +      +(q3+dble(igL3))*b3(:)    
      qgL2=qgL(1)**2+qgL(2)**2+qgL(3)**2
      if(qgL2<=Gcut) then 
       igL=igL+1 
       KG0(1,igL)=igL1
       KG0(2,igL)=igL2
       KG0(3,igL)=igL3 
      endif  
      enddo 
      enddo 
      enddo 
      NG=igL 
!--
!write(6,*) "maxabs KG0_1=",maxval(abs(KG0(1,:))) 
!write(6,*) "maxabs KG0_2=",maxval(abs(KG0(2,:))) 
!write(6,*) "maxabs KG0_3=",maxval(abs(KG0(3,:))) 
!--
!do igL=1,NG 
! write(6,*) igL,KG0(:,igL)
!enddo 
!--
      RETURN 
      END 
