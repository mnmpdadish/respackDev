      subroutine search_kq(NTK,SK0,q1,q2,q3,ik,ikq,shift_G)
      implicit none 
      integer::NTK
      real(8)::SK0(3,NTK)  
      real(8)::q1,q2,q3
      integer::ik,jk
      real(8)::SKQ(3)
      integer::ikq,shift_G(3) 
      real(8),parameter::dlt_BZ=1.0d-6!20170316 
!--
      SKQ(1)=SK0(1,ik)+q1
      SKQ(2)=SK0(2,ik)+q2 
      SKQ(3)=SK0(3,ik)+q3 
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
!      write(6,'(a6,3f15.7)')'SKQ',SKQ
      enddo 
      RETURN  
      END 
!--
      subroutine search_list(SQI,SQ)
      implicit none 
      real(8)::SQI(3),SQ(3)
      real(8),parameter::dlt_BZ=1.0d-6!20170316 
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
      if(SQ(3)<=-0.50d0+dlt_BZ) SQ(3)=SQ(3)+1.0d0 
      if(SQ(3)<=-1.50d0+dlt_BZ) SQ(3)=SQ(3)+2.0d0 
      RETURN  
      END 
!--
      subroutine calc_InterStateMatrix(NTK,NTG,ncomp,NG0,KG0,
     +                                 C0_K,C0_KQ,
     +                                 ik,ikq,nwx2,nwy2,nwz2,
     +                                 nfft1,nfft2,Nl123,
     +                                 wfunc,fftwk,fs,LG0,NG_for_eps,
     +                                 shift_G,m_tmp)
!--
      use fft_3d 
      implicit none 
      type(fft3_struct)::fs 
      integer::NTK,NTG,ncomp
      integer::NG0(NTK),KG0(3,NTG,NTK) 
      complex(8)::C0_K(NTG,ncomp),C0_KQ(NTG,ncomp) 
      integer::ik,ikq 
      integer::nwx2,nwy2,nwz2,nfft1,nfft2,Nl123 
      real(8)::wfunc(Nl123*2),fftwk(Nl123*2)
      integer::NG_for_eps 
      integer::LG0(3,NTG) 
      integer::shift_G(3)
      integer::ig,igb1,igb2,igb3,ind,igL,igL1,igL2,igL3, ic
      integer::ir,ir1,ir2,ir3 
      complex(8)::cell_periodic_k(Nl123)  
      complex(8)::cell_periodic_kq(Nl123)  
      complex(8)::f
      real(8)::dG1,dG2,dG3
      real(8)::phase 
      complex(8)::pf 
      complex(8)::m_tmp(NG_for_eps) 
      real(8),parameter::pi=dacos(-1.0d0)
      real(8),parameter::tpi=2.0d0*pi 
      complex(8),parameter::ci=(0.0D0,1.0D0) 
!--

      m_tmp(:)=0.0d0
      do ic=1,ncomp ! we sum spins: m_tmp = <up|up> + <down|down>
      
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
       wfunc(ind)=dble(C0_K(ig,ic))
       wfunc(ind+Nl123)=dimag(C0_K(ig,ic))
      enddo 
      call fft3_bw(fs,wfunc(1),fftwk(1)) 
      do ir=1,Nl123 
       cell_periodic_k(ir)=cmplx(wfunc(ir),wfunc(ir+Nl123)) 
      enddo 
!     write(6,*) 'cell_k',cell_periodic_k  
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
       wfunc(ind)=dble(C0_KQ(ig,ic))
       wfunc(ind+Nl123)=dimag(C0_KQ(ig,ic)) 
      enddo 
      call fft3_bw(fs,wfunc,fftwk)    
      do ir=1,Nl123
       cell_periodic_kq(ir)=cmplx(wfunc(ir),wfunc(ir+Nl123)) 
      enddo 
!    write(6,*) 'cell_kq',cell_periodic_kq  
!--
      dG1=dble(shift_G(1))
      dG2=dble(shift_G(2))
      dG3=dble(shift_G(3))
!     write(6,*) 'dG=',dG1,dG2,dG3
!--
      do ir3=1,nwz2
       do ir2=1,nwy2
        do ir1=1,nwx2
         ir=ir1+(ir2-1)*nfft1+(ir3-1)*nfft1*nfft2 
         phase=tpi*(dble(ir1-1)*dG1/dble(nwx2)
     +             +dble(ir2-1)*dG2/dble(nwy2) 
     +             +dble(ir3-1)*dG3/dble(nwz2)) 
         pf=exp(-ci*phase) 
         cell_periodic_kq(ir)
     +  =cell_periodic_kq(ir)*pf 
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
       m_tmp(igL)=m_tmp(igL)
     +   +cmplx(wfunc(ind),wfunc(ind+Nl123)) !adding
      enddo
      
      enddo !ic=1,ncomp 
!--
      RETURN  
      END 
!--
      subroutine calc_VMab(NTK,NTG,ncomp,NG0,KG0,C0_K,C0_KQ,
     +                                        ik,b1,b2,b3,vm)  
      implicit none 
      integer::NTK,NTG,ik, ncomp
      integer::NG0(NTK),KG0(3,NTG,NTK) 
      real(8)::b1(3),b2(3),b3(3) 
      complex(8)::C0_K(NTG,ncomp),C0_KQ(NTG,ncomp) 
      integer::ig,igb1,igb2,igb3, ic
      real(8)::G_VEC(3)
      complex(8)::SUM_VEC(3)
      complex(8)::vm(3) 
      complex(8),parameter::ci=(0.0D0,1.0D0) 
!--
      SUM_VEC(:)=0.0D0 
      do ic=1,ncomp
       do ig=1,NG0(ik) 
        igb1=KG0(1,ig,ik) 
        igb2=KG0(2,ig,ik) 
        igb3=KG0(3,ig,ik) 
        G_VEC(:)=dble(igb1)*b1(:)+dble(igb2)*b2(:)+dble(igb3)*b3(:) 
        SUM_VEC(:)=SUM_VEC(:)+G_VEC(:)*
     +      CONJG(C0_KQ(ig,ic))*C0_K(ig,ic)
       enddo!ig
      enddo!ic
!"vm" depends on the direction of q approaching the zero limit
!We employ q=[001],[010],[100] therefore
!Kazuma Nakamura bugfix 2010 10 12
      vm(:)=SUM_VEC(:)*ci!ci is multiplied
      RETURN  
      END 
!--
      subroutine calc_VMaa(NTK,NTG,ncomp,NG0,KG0,C0_K,C0_KQ,
     +                                    ik,b1,b2,b3,SKT,vm)  
      implicit none 
      integer::NTK,NTG,ik, ncomp
      integer::NG0(NTK),KG0(3,NTG,NTK) 
      real(8)::b1(3),b2(3),b3(3),SKT(3) 
      complex(8)::C0_K(NTG,ncomp),C0_KQ(NTG,ncomp) 
      integer::ig,igb1,igb2,igb3,ic
      real(8)::G_VEC(3)
      complex(8)::SUM_VEC(3)
      complex(8) :: vm(3) 
      complex(8),parameter::ci=(0.0D0,1.0D0) 
!--
      SUM_VEC(:)=0.0D0 
      do ic=1,ncomp
       do ig=1,NG0(ik) 
        igb1=KG0(1,ig,ik) 
        igb2=KG0(2,ig,ik) 
        igb3=KG0(3,ig,ik) 
        G_VEC(:)=(SKT(1)+dble(igb1))*b1(:)
     +          +(SKT(2)+dble(igb2))*b2(:)
     +          +(SKT(3)+dble(igb3))*b3(:) 
        SUM_VEC(:)=SUM_VEC(:)+G_VEC(:)*
     +       CONJG(C0_KQ(ig,ic))*C0_K(ig,ic)
       enddo!ig
      enddo!ic
!"vm" depends on the direction of q approaching the zero limit 
!We employ q=[001],[010],[100] therefore
!Kazuma Nakamura bugfix 2010 10 12
      vm(:)=SUM_VEC(:)*ci!ci is multiplied  
      RETURN
      END
!--
      subroutine make_C0(NTG,ncomp,rinv_SO,itrs,NG,KGtmp,RWtmp,
     +           rginvtmp,pgtmp,nnp,L1,L2,L3,packtmp,OCCtmp,C0_K) 
      implicit none 
      integer::NTG,ncomp,L1,L2,L3,nnp  
      integer::itrs 
      integer::NG
      integer::KGtmp(3,NTG) 
      integer::RWtmp(3) 
      real(8)::rginvtmp(3,3) 
      integer::pgtmp(3) 
      integer::packtmp(-L1:L1,-L2:L2,-L3:L3) 
      complex(8)::OCCtmp(NTG,ncomp) 
      integer::ig,jg,i1,i2,i3,j1,j2,j3,k1,k2,k3 
      real(8)::phase 
      complex(8)::pf,rinv_SO(2,2)
      complex(8)::C0_K(NTG,ncomp) 
      real(8),parameter::pi=dacos(-1.0d0)
      real(8),parameter::tpi=2.0d0*pi 
      complex(8),parameter::ci=(0.0D0,1.0D0) 
!--
      C0_K(:,:)=0.0d0 
      if(itrs==1) then 
       do ig=1,NG 
        i1=KGtmp(1,ig); j1=KGtmp(2,ig); k1=KGtmp(3,ig) 
        i2=i1+RWtmp(1); j2=j1+RWtmp(2); k2=k1+RWtmp(3) 
        i3=int(rginvtmp(1,1))*i2
     +    +int(rginvtmp(1,2))*j2 
     +    +int(rginvtmp(1,3))*k2 
        j3=int(rginvtmp(2,1))*i2
     +    +int(rginvtmp(2,2))*j2 
     +    +int(rginvtmp(2,3))*k2 
        k3=int(rginvtmp(3,1))*i2
     +    +int(rginvtmp(3,2))*j2 
     +    +int(rginvtmp(3,3))*k2 
        jg=packtmp(i3,j3,k3) 
        phase=tpi*(dble(i1)*dble(pgtmp(1))
     +            +dble(j1)*dble(pgtmp(2))
     +            +dble(k1)*dble(pgtmp(3))) 
        pf=exp(-ci*phase/dble(nnp)) 
!        C0_K(ig)=OCCtmp(jg)*pf 

        if(ncomp==2) then
         C0_K(ig,1)=rinv_SO(1,1)*OCCtmp(jg,1)*pf
     +             +rinv_SO(1,2)*OCCtmp(jg,2)*pf
         C0_K(ig,2)=rinv_SO(2,1)*OCCtmp(jg,1)*pf
     +             +rinv_SO(2,2)*OCCtmp(jg,2)*pf
        else
         C0_K(ig,:)=OCCtmp(jg,:)*pf
        endif

       enddo!ig 
      elseif(itrs==-1) then  
       do ig=1,NG 
        i1=-KGtmp(1,ig); j1=-KGtmp(2,ig); k1=-KGtmp(3,ig) 
        i2=i1+RWtmp(1); j2=j1+RWtmp(2); k2=k1+RWtmp(3) 
        i3=int(rginvtmp(1,1))*i2
     +    +int(rginvtmp(1,2))*j2 
     +    +int(rginvtmp(1,3))*k2 
        j3=int(rginvtmp(2,1))*i2
     +    +int(rginvtmp(2,2))*j2 
     +    +int(rginvtmp(2,3))*k2 
        k3=int(rginvtmp(3,1))*i2
     +    +int(rginvtmp(3,2))*j2 
     +    +int(rginvtmp(3,3))*k2 
        jg=packtmp(i3,j3,k3) 
        phase=tpi*(dble(i1)*dble(pgtmp(1))
     +            +dble(j1)*dble(pgtmp(2))
     +            +dble(k1)*dble(pgtmp(3))) 
        pf=exp(-ci*phase/dble(nnp)) 

        !C0_K(ig)=OCCtmp(jg)*pf 
        if(ncomp==2) then
         C0_K(ig,2)= conjg(rinv_SO(1,1)*OCCtmp(jg,1)*pf
     +                    +rinv_SO(1,2)*OCCtmp(jg,2)*pf)
         C0_K(ig,1)=-conjg(rinv_SO(2,1)*OCCtmp(jg,1)*pf
     +                    +rinv_SO(2,2)*OCCtmp(jg,2)*pf)
        else
         C0_K(ig,:)=conjg(OCCtmp(jg,:)*pf)
        endif 

       enddo!ig 
!       C0_K(:)=conjg(C0_K(:))  
      endif 
      return 
      end 
!---
      SUBROUTINE make_index_kpt(NTK,nkb1,nkb2,nkb3,SK0,index_kpt)      
      implicit none 
      integer::NTK,nkb1,nkb2,nkb3
      real(8)::SK0(3,NTK)  
      integer::ik,ix,iy,iz
      real(8)::x,y,z
      integer::index_kpt(nkb1,nkb2,nkb3)    
 
      if(MOD(NTK,2)/=0) then 
!      write(6,*)'i am in make_index for odd'
       do ik=1,NTK 
        x=SK0(1,ik)*dble(nkb1) 
        y=SK0(2,ik)*dble(nkb2)
        z=SK0(3,ik)*dble(nkb3)  
        x=x+(dble(nkb1)-1.0d0)/2.0d0 
        y=y+(dble(nkb2)-1.0d0)/2.0d0
        z=z+(dble(nkb3)-1.0d0)/2.0d0 
        ix=idnint(x)+1
        iy=idnint(y)+1
        iz=idnint(z)+1
        index_kpt(ix,iy,iz)=ik
!        write(6,*) ix,iy,iz
       enddo 
      else!20170316 
!      write(6,*)'i am in make_index for even'
       do ik=1,NTK 
!       x=SK0(1,ik)*2.0D0*dble(nkb1)+1.0D0  
!       y=SK0(2,ik)*2.0D0*dble(nkb2)+1.0D0  
!       z=SK0(3,ik)*2.0D0*dble(nkb3)+1.0D0  
!       x=(x+dble(nkb1))/2.0D0 
!       y=(y+dble(nkb2))/2.0D0 
!       z=(z+dble(nkb3))/2.0D0 
!       ix=idnint(x)
!       iy=idnint(y)
!       iz=idnint(z)
        x=SK0(1,ik)*dble(nkb1) 
        y=SK0(2,ik)*dble(nkb2)
        z=SK0(3,ik)*dble(nkb3)  
        x=x+dble(nkb1)/2.0d0 
        y=y+dble(nkb2)/2.0d0
        z=z+dble(nkb3)/2.0d0 
        ix=idnint(x)
        iy=idnint(y)
        iz=idnint(z)
        index_kpt(ix,iy,iz)=ik
!       write(6,*) ix,iy,iz
       enddo 
      endif 
!--
!      do iz=1,nkb3 
!       do iy=1,nkb2
!        do ix=1,nkb1
!         ik=index_kpt(ix,iy,iz)
!         write(6,'(i8,3f15.10)') ik,SK0(:,ik)
!        enddo
!       enddo
!      enddo
!      write(6,*)'finish make_index_kpt'
!--
      RETURN  
      END 
!---
      subroutine calc_eps_rpa(NTG,NG_for_eps,LG0,
     +                   q1,q2,q3,b1,b2,b3,
     +                   chi0,tpi,eps_rpa, 
     +                   file_num_chi,file_num_eps)  
      implicit none 
      integer::NTG,NG_for_eps  
      integer::LG0(3,NTG)    
      real(8)::q1,q2,q3,b1(3),b2(3),b3(3),tpi 
      integer::file_num_chi,file_num_eps 
      complex(8)::chi0(NG_for_eps,NG_for_eps)
      integer::ig,jg,igL1,igL2,igL3,jgL1,jgL2,jgL3 
      real(8)::qgL(3),qgL1,qgL2,qgLi,qgLj   
      complex(8)::f
      complex(8)::eps_rpa(NG_for_eps,NG_for_eps) 
!
      if(file_num_chi==0) write(6,*)'not write chi'
!Factor of 2 is needed for only LDA-polarization calculation, pointed by Yohiro Nohara, 2009 5 27 ---* 
      chi0(:,:)=2.0d0*chi0(:,:) 
!cal_eps_rpa 
      write(6,*)'calc eps_rpa'
      write(6,*)' ' 
      eps_rpa(:,:)=0.0D0 
      do ig=1,NG_for_eps 
       igL1=LG0(1,ig)
       igL2=LG0(2,ig)
       igL3=LG0(3,ig)
       qgL(:)=(q1+dble(igL1))*b1(:)
     +       +(q2+dble(igL2))*b2(:)
     +       +(q3+dble(igL3))*b3(:)
       qgL2=qgL(1)**2+qgL(2)**2+qgL(3)**2
       qgLi=dsqrt(qgL2) 
      do jg=1,NG_for_eps   
       jgL1=LG0(1,jg)
       jgL2=LG0(2,jg)
       jgL3=LG0(3,jg)
       qgL(:)=(q1+dble(jgL1))*b1(:)
     +       +(q2+dble(jgL2))*b2(:)
     +       +(q3+dble(jgL3))*b3(:)
       qgL2=qgL(1)**2+qgL(2)**2+qgL(3)**2
       qgLj=dsqrt(qgL2) 
!diag 
       if(ig==jg)then 
!       write(file_num_chi,'(3F20.10)') qgLi,chi0(ig,jg)
        eps_rpa(ig,jg)=1.0D0-chi0(ig,jg)*2.0D0*tpi/qgLi/qgLj 
       endif 
!off diag
       if(ig/=jg)then 
        eps_rpa(ig,jg)=-chi0(ig,jg)*2.0D0*tpi/qgLi/qgLj 
       endif 
      enddo 
      enddo 
!---
      call invZGE(NG_for_eps,eps_rpa(1,1)) 
!---
      write(6,*)'calc eps_rpa_inv'
      write(6,*)' ' 
      do ig=1,NG_for_eps 
       igL1=LG0(1,ig)
       igL2=LG0(2,ig)
       igL3=LG0(3,ig)
       qgL(:)=(q1+dble(igL1))*b1(:)
     +       +(q2+dble(igL2))*b2(:)
     +       +(q3+dble(igL3))*b3(:)
       qgL2=qgL(1)**2+qgL(2)**2+qgL(3)**2
       qgL1=dsqrt(qgL2) 
       f=eps_rpa(ig,ig)
!      write(file_num_eps,'(3F20.10)') qgL1,dble(1.0D0/f),dimag(1.0d0/f)  
       write(file_num_eps,'(3F20.10)') qgL1,dble(f),dimag(f)  
      enddo 
!---
      return 
      end 
!---
      subroutine calc_eps_rpa_q_0(NTG,NG_for_eps,LG0,
     +                   q1,q2,q3,b1,b2,b3,
     +                   chi0,tpi,No_G_0,eps_rpa,
     +                   wd,delt,wcmplx,   
     +                   file_num_chi,file_num_eps)  
      implicit none 
      integer::NTG,NG_for_eps,No_G_0,file_num_chi,file_num_eps 
      integer::LG0(3,NTG)    
      real(8)::q1,q2,q3,b1(3),b2(3),b3(3),tpi,wd,delt,w  
      complex(8)::chi0(NG_for_eps,NG_for_eps)
      complex(8)::wcmplx,f
      integer::ig,jg,igL1,igL2,igL3,jgL1,jgL2,jgL3 
      real(8)::qgL(3),qgL1,qgL2,qgLi,qgLj   
      complex(8)::eps_rpa(NG_for_eps,NG_for_eps) 
!
      if(file_num_chi==0) write(6,*)'not write chi'
!Factor of 2 is needed for only LDA-polarization calculation, pointed by Yohiro Nohara, 2009 5 27 ---* 
      chi0(:,:)=2.0d0*chi0(:,:) 
      w=dble(wcmplx)  
      if(w==0.0d0)then 
       w=1.0d-5 
       write(6,*) wcmplx,w 
      endif 
!     wd=0.0d0 
!calc_eps_rpa_q=0 
      write(6,*)'calc eps_rpa_q=0'
      write(6,*)' ' 
      eps_rpa(:,:)=0.0D0 
      do ig=1,NG_for_eps 
       igL1=LG0(1,ig)
       igL2=LG0(2,ig)
       igL3=LG0(3,ig)
       qgL(:)=(q1+dble(igL1))*b1(:)
     +       +(q2+dble(igL2))*b2(:)
     +       +(q3+dble(igL3))*b3(:)
       qgL2=qgL(1)**2+qgL(2)**2+qgL(3)**2
       qgLi=dsqrt(qgL2) 
      do jg=1,NG_for_eps   
       jgL1=LG0(1,jg)
       jgL2=LG0(2,jg)
       jgL3=LG0(3,jg)
       qgL(:)=(q1+dble(jgL1))*b1(:)
     +       +(q2+dble(jgL2))*b2(:)
     +       +(q3+dble(jgL3))*b3(:)
       qgL2=qgL(1)**2+qgL(2)**2+qgL(3)**2
       qgLj=dsqrt(qgL2) 
!head 
       if(ig==jg.and.ig==No_G_0)then 
!       write(file_num_chi,'(3F20.10)') qgLi,chi0(ig,jg)
        eps_rpa(ig,jg)
     + =1.0D0-chi0(ig,jg)*2.0D0*tpi 
     + +cmplx(-wd**2/(w**2+delt**2),delt*wd**2/(w**3+w*delt**2))
       endif 
!body(diag)
       if(ig==jg.and.ig/=No_G_0)then 
!       write(file_num_chi,'(3F20.10)') qgLi,chi0(ig,jg)
        eps_rpa(ig,jg)=1.0D0-chi0(ig,jg)*2.0D0*tpi/qgLi/qgLj  
       endif 
!wing(column) 
       if(ig/=jg.and.ig==No_G_0)then 
        eps_rpa(ig,jg)=-chi0(ig,jg)*2.0D0*tpi/qgLj  
       endif 
!wing(raw) 
       if(ig/=jg.and.jg==No_G_0)then 
        eps_rpa(ig,jg)=-chi0(ig,jg)*2.0D0*tpi/qgLi
       endif      
!body(off diag)
       if(ig/=jg.and.ig/=No_G_0.and.jg/=No_G_0) then 
        eps_rpa(ig,jg)=-chi0(ig,jg)*2.0D0*tpi/qgLi/qgLj 
       endif 
      enddo 
      enddo 
!---
!      write(6,*)' '
!      write(6,*)'eps(i,G0)'
!      do ig=1,NG_for_eps
!       write(6,*) ig,eps_rpa(ig,No_G_0)
!      enddo 
!      write(6,*)' '
!      write(6,*)'eps(G0,j)'
!      do ig=1,NG_for_eps
!       write(6,*) ig,eps_rpa(No_G_0,ig)
!      enddo 
!      write(6,*)' '
!---
!      do ig=1,NG_for_eps
!       do jg=1,NG_for_eps
!        write(6,'(2i8,2f15.8)') ig,jg,eps_rpa(ig,jg)
!       enddo 
!      enddo 
!---
      call invZGE(NG_for_eps,eps_rpa(1,1)) 
!---
      write(6,*)'calc eps_rpa_inv_q=0'
      write(6,*)' ' 
      do ig=1,NG_for_eps 
       igL1=LG0(1,ig)
       igL2=LG0(2,ig)
       igL3=LG0(3,ig)
       qgL(:)=(q1+dble(igL1))*b1(:)
     +       +(q2+dble(igL2))*b2(:)
     +       +(q3+dble(igL3))*b3(:)
       qgL2=qgL(1)**2+qgL(2)**2+qgL(3)**2
       qgL1=dsqrt(qgL2) 
       f=eps_rpa(ig,ig)
!     write(file_num_eps,'(3F20.10)') qgL1,dble(1.0D0/f),dimag(1.0d0/f) 
      write(file_num_eps,'(3F20.10)') qgL1,dble(f),dimag(f)  
      enddo 
      return 
      end 
!--
      subroutine judge_FermiInside(NTK,NTB,E_EIG,FermiEnergy,
     + WindowInside,Tindx,E_AVE,band_max,band_min) 
      implicit none 
      integer::NTK,NTB
      real(8)::E_EIG(NTB,NTK)  
      real(8)::FermiEnergy
      integer::ik,ib,fs 
      real(8)::e,top,btm 
      integer::WindowInside(NTB)    
      integer::Tindx(NTB)    
      real(8)::E_AVE(NTB)    
      real(8)::band_max(NTB)
      real(8)::band_min(NTB)
!--
      WindowInside(:)=0
      Tindx(:)=0
      do ib=1,NTB
       do ik=1,NTK 
        e=E_EIG(ib,ik) 
        if(ik==1)then 
         top=e
         btm=e
        endif 
        if(e<btm) btm=e
        if(e>top) top=e
       enddo!ik 
!uocc
       if(btm>FermiEnergy) fs=2
!pocc
       if(btm<=FermiEnergy.and.top>=FermiEnergy) then
        fs=1
        Tindx(ib)=ib 
       endif 
!docc
       if(top<FermiEnergy) fs=0
!
       WindowInside(ib)=fs 
       E_AVE(ib)=(top+btm)/2.0d0 
       band_max(ib)=top
       band_min(ib)=btm
      enddo!ib 
!--
      RETURN 
      END 
!--
      subroutine judge_WindowInside(NTK,NTB,E_EIG,EL,EU,WindowInside,
     + Tindx,E_AVE,band_max,band_min)
      implicit none 
      integer::NTK,NTB
      real(8)::E_EIG(NTB,NTK)  
      real(8)::EL,EU
      integer::ik,ib,fs 
      real(8)::e,top,btm 
      integer::WindowInside(NTB)    
      integer::Tindx(NTB)    
      real(8)::E_AVE(NTB)    
      real(8)::band_max(NTB)
      real(8)::band_min(NTB)
!--
      WindowInside(:)=0
      Tindx(:)=0
      do ib=1,NTB
       do ik=1,NTK 
        e=E_EIG(ib,ik) 
        if(ik==1)then 
         top=e
         btm=e
        endif 
        if(e<btm) btm=e
        if(e>top) top=e
       enddo!ik 
!case1:docc
       if(top<EL) fs=0
!case2:docc
       if(EL<=top.and.top<EU.and.btm<EL) fs=0 
!case3:pocc
       if(EL<=top.and.top<EU.and.EL<=btm.and.btm<EU)then 
        fs=1
        Tindx(ib)=ib
       endif   
!case4:pocc
       if(btm<EL.and.top>EU)then   
        fs=1 
        Tindx(ib)=ib
       endif   
!case5:uocc
       if(EL<btm.and.btm<EU.and.top>EU) fs=2 
!case6:uocc 
       if(btm>EU) fs=2 
!--
       WindowInside(ib)=fs 
       E_AVE(ib)=(top+btm)/2.0d0 
       band_max(ib)=top
       band_min(ib)=btm
      enddo!ib 
      RETURN 
      END 
!--
      subroutine kcheck(ktmp,RWtmp) 
      implicit none 
      real(8),intent(inout)::ktmp(3)
      integer,intent(out)::RWtmp(3) 
      real(8),parameter::dlt_BZ=1.0d-6!20170316 
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
!--
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
!20170316 
      subroutine kcheck_trs(ktmp,RWtmp)!20170316 
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
      subroutine make_KG0(NTG,b1,b2,b3,Gcut,q1,q2,q3,KG0,NG) 
      implicit none 
      integer,intent(in)::NTG
      real(8),intent(in)::b1(3),b2(3),b3(3) 
      real(8),intent(in)::q1,q2,q3,Gcut 
      integer::igL,igL1,igL2,igL3
      real(8)::qgL(3),qgL2  
      integer,intent(out)::KG0(3,NTG)    
      integer,intent(out)::NG 
      integer,parameter::NGL1=100
      integer,parameter::NGL2=100 
      integer,parameter::NGL3=100
      igL=0
      do igL1=-NGL1,NGL1 
       do igL2=-NGL2,NGL2 
        do igL3=-NGL3,NGL3 
         qgL(:)=(q1+dble(igL1))*b1(:)
     +         +(q2+dble(igL2))*b2(:)
     +         +(q3+dble(igL3))*b3(:)    
         qgL2=qgL(1)**2+qgL(2)**2+qgL(3)**2
         if(qgL2<=Gcut)then 
          igL=igL+1 
          KG0(1,igL)=igL1
          KG0(2,igL)=igL2
          KG0(3,igL)=igL3 
         endif  
        enddo 
       enddo 
      enddo 
      NG=igL 
!---
!write(6,*)"maxabs KG0_1=",maxval(abs(KG0(1,:))) 
!write(6,*)"maxabs KG0_2=",maxval(abs(KG0(2,:))) 
!write(6,*)"maxabs KG0_3=",maxval(abs(KG0(3,:))) 
!---
!do igL=1,NG 
! write(6,*) igL,KG0(:,igL)
!enddo 
      RETURN 
      END 
!---
      SUBROUTINE OUTER_PRODUCT(vec_x,vec_y,vec_z)
      implicit none 
      real(8)::vec_x(3),vec_y(3),vec_z(3) 
      vec_z(1)=vec_x(2)*vec_y(3)-vec_x(3)*vec_y(2)
      vec_z(2)=vec_x(3)*vec_y(1)-vec_x(1)*vec_y(3) 
      vec_z(3)=vec_x(1)*vec_y(2)-vec_x(2)*vec_y(1)
      RETURN
      END
