      PROGRAM WANNIER          
      use m_rdinput       
      use fft_3d
      use m_bvector, only: generate_bvectors
      use m_frmsf_wan, only: wrt_frmsf_wan
      use m_dos, only: calculate_dos 
      use m_eigenstate, only: calculate_eigenstate
      use m_wrt_frmsf, only: wrt_frmsf 
      use m_wrt_mvmc  
      use m_dmx, only: calculate_dmx 
      include "config.h"
      !--
      call read_input
      NGAUSS=n_occ!TOTAL NUMBER OF WANNIER ORBITAL  
      E_LOWER=E_LOWER/au 
      E_UPPER=E_UPPER/au 
      E_LOWER_inner=E_LOWER_inner/au 
      E_UPPER_inner=E_UPPER_inner/au 
      write(6,*)'E_LOWER in au',E_LOWER 
      write(6,*)'E_UPPER in au',E_UPPER 
      write(6,*)'E_LOWER_inner in au',E_LOWER_inner 
      write(6,*)'E_UPPER_inner in au',E_UPPER_inner 
!--
!      do ig=1,nigs 
!       write(6,*) vec_ini(ig) 
!      enddo 
!--
!      do ix=1,3 
!       write(6,'(100f10.5)')(SK_sym_pts(ix,ik),ik=1,N_sym_points)
!      enddo 
!--
      WRITE(6,*) 
      WRITE(6,*)'=========================================='
      WRITE(6,*)'=== WANNIER FUNCTION CALCULATION START ==='
      WRITE(6,*)'=========================================='
      WRITE(6,*) 
!--
!OPEN(117,R,FILE='dat.bandcalc') 
      OPEN(117,FILE='./dir-wfn/dat.bandcalc') 
      rewind(117) 
      read(117,*) Ecut_for_psi 
      read(117,*) FermiEnergy  
      read(117,*) Etot
      write(6,*)'Ecut_for_psi=',Ecut_for_psi 
      write(6,*)'FermiEnergy=',FermiEnergy  
      write(6,*)'Etot=',Etot
!--
!OPEN(105,R,FILE='dat.lattice')
      OPEN(105,FILE='./dir-wfn/dat.lattice') 
      REWIND(105)
      READ(105,*) a1(1),a1(2),a1(3)!a1 vector
      READ(105,*) a2(1),a2(2),a2(3)!a2 vector
      READ(105,*) a3(1),a3(2),a3(3)!a3 vector
      CLOSE(105)
      !--
      call OUTER_PRODUCT(a2(1),a3(1),b1(1))
      VOLUME=a1(1)*b1(1)+a1(2)*b1(2)+a1(3)*b1(3)
      b1(:)=b1(:)*tpi/VOLUME 
      call OUTER_PRODUCT(a3(1),a1(1),b2(1))
      b2(:)=b2(:)*tpi/VOLUME 
      call OUTER_PRODUCT(a1(1),a2(1),b3(1))
      b3(:)=b3(:)*tpi/VOLUME 
      !--
      write(6,*) 
      write(6,*)'=============='
      write(6,*)'LATTICE VECTOR'
      write(6,*)'=============='
      write(6,*) 
      write(6,*) a1(1), a1(2), a1(3)
      write(6,*) a2(1), a2(2), a2(3)
      write(6,*) a3(1), a3(2), a3(3)
      write(6,*) 
      write(6,*)'VOLUME OF UNIT CELL=',VOLUME  
      write(6,*) 
      write(6,*)'========================='
      write(6,*)'RECIPROCAL LATTICE VECTOR'
      write(6,*)'========================='
      write(6,*) 
      write(6,*) b1(1), b1(2), b1(3)
      write(6,*) b2(1), b2(2), b2(3)
      write(6,*) b3(1), b3(2), b3(3)
      write(6,*) 
      write(6,*)'======================'
      write(6,*)'CHECK OF (A_VEC*B_VEC)'
      write(6,*)'======================'
      write(6,*) 
      write(6,*)'a1*b1=',b1(1)*a1(1)+b1(2)*a1(2)+b1(3)*a1(3)
      write(6,*)'a2*b2=',b2(1)*a2(1)+b2(2)*a2(2)+b2(3)*a2(3)
      write(6,*)'a3*b3=',b3(1)*a3(1)+b3(2)*a3(2)+b3(3)*a3(3)
      call est_latparam(a1(1),a2(1),a3(1),a,b,c,alp,bet,gmm) 
      write(6,'(6f15.10)') a,b,c,alp,bet,gmm 
!--
!OPEN(107,R,FILE='dat.atom_kind')
!      OPEN(107,FILE='./dir-wfn/dat.atom_kind') 
!      rewind(107)
!      read(107,*) nkd 
!      allocate(zo(nkd));zo=0.0d0
!      allocate(zn(nkd));zn=0.0d0
!      do i=1,nkd 
!       read(107,*) zo(i),zn(i)
!      enddo
!      close(107)
!--
!OPEN(108,R,FILE='dat.atom_position')
!      OPEN(108,FILE='./dir-wfn/dat.atom_position') 
!      rewind(108)
!      read(108,*) nsi 
!      allocate(kd(nsi));kd=0 
!      allocate(asi(3,nsi));asi=0.0d0
!      do i=1,nsi  
!       read(108,*) kd(i),(asi(j,i),j=1,3) 
!      enddo 
!      close(108)
!--
!OPEN(108,R,FILE='dat.atom_position')
      OPEN(108,FILE='./dir-wfn/dat.atom_position') 
      rewind(108)
      read(108,*) nsi 
      allocate(chemical_species(nsi))
      allocate(asi(3,nsi));asi=0.0d0
      do i=1,nsi  
       read(108,*) chemical_species(i),(asi(j,i),j=1,3) 
      enddo 
      close(108)
!--
!OPEN(100,R,FILE='dat.symmetry') 
      OPEN(100,FILE='./dir-wfn/dat.symmetry') 
      rewind(100) 
      read(100,*) nsymq 
      read(100,*) nnp 
      write(6,*)'nsymq=',nsymq 
      write(6,*)'nnp=',nnp 
      allocate(rg(3,3,nsymq));rg=0
      allocate(pg(3,nsymq));pg=0
      allocate(rginv(3,3,nsymq));rginv=0.0d0 
      do iop=1,nsymq
       read(100,*)((rg(i,j,iop),i=1,3),j=1,3) 
       read(100,*)(pg(i,iop),i=1,3)   
      enddo 
      close(100) 
      !--
      write(6,*) 
      write(6,*)'========'
      write(6,*)'SYMMETRY'
      write(6,*)'========'
      write(6,*) 
      do iop=1,nsymq  
       write(6,'(a5,i4)')'sym=',iop
       write(6,'(a10)')'rg and pg' 
       do i=1,3
        write(6,'(3i5,5x,i5)')(rg(i,j,iop),j=1,3),pg(i,iop) 
       enddo
       write(6,*) 
      enddo
      !--
      rginv=rg 
      do iop=1,nsymq
       call invmat(3,rginv(1,1,iop)) 
      enddo 
      !--
      do iop=1,nsymq  
       write(6,'(a5,i4)')'sym=',iop
       write(6,'(a10)')'rginv' 
       do i=1,3
        write(6,'(3f10.5)')(rginv(i,j,iop),j=1,3) 
       enddo
       write(6,*) 
      enddo
      write(6,*)'finish read dat.symmetry'
!--
!OPEN(101,R,FILE='dat.sample-k') 
      OPEN(101,FILE='./dir-wfn/dat.sample-k') 
      rewind(101) 
      read(101,*) Nk_irr 
      allocate(SKI(3,Nk_irr));SKI(:,:)=0.0D0 
      do ik=1,Nk_irr 
       read(101,*)(SKI(i,IK),i=1,3) 
      enddo 
      close(101) 
      !call est_NTK(NK_irr,SKI(1,1),NTK,nkb1,nkb2,nkb3)  
      !write(6,'(a24,4i10)')'nkb1,nkb2,nkb3,NTK=',nkb1,nkb2,nkb3,NTK  
      !Na1=nkb1/2;Na2=nkb2/2;Na3=nkb3/2
      !write(6,'(a24,3i10)')'Na1,Na2,Na3=',Na1,Na2,Na3 
!--
!OPEN(132,R,FILE='dat.nkm') 
      OPEN(132,FILE='./dir-wfn/dat.nkm') 
      allocate(NGI(Nk_irr));NGI(:)=0
      rewind(132)
      do ik=1,Nk_irr 
       read(132,*) NGI(ik) 
      enddo 
      close(132) 
      NTG=maxval(abs(NGI(:))) 
      write(6,'(a10,i10)')'NTG=',NTG  
!--
!OPEN(104,R,FILE='dat.kg') 
      OPEN(104,FILE='./dir-wfn/dat.kg') 
      rewind(104) 
      allocate(KGI(3,NTG,Nk_irr));KGI=0 
      do ik=1,Nk_irr 
       read(104,*) NKG_tmp
       if(NKG_tmp/=NGI(ik)) then 
        write(6,*)'ERROR; STOP; NKG_tmp should be NGI(ik)'   
        write(6,*)'NKG_tmp=',NKG_tmp,'NGI(ik)=',NGI(ik)
        stop 
       endif 
       do ig=1,NKG_tmp
        read(104,*)(KGI(i,ig,ik),i=1,3) 
       enddo!ig 
      enddo!ik  
      close(104)
      !--
      L1=maxval(abs(KGI(1,:,:)))+1; write(6,*)'L1=',L1 
      L2=maxval(abs(KGI(2,:,:)))+1; write(6,*)'L2=',L2 
      L3=maxval(abs(KGI(3,:,:)))+1; write(6,*)'L3=',L3 
      allocate(packing(-L1:L1,-L2:L2,-L3:L3,Nk_irr))   
      packing(:,:,:,:)=0 
      do ik=1,Nk_irr 
       do ig=1,NGI(ik) 
        i1=KGI(1,ig,ik); j1=KGI(2,ig,ik); k1=KGI(3,ig,ik) 
        packing(i1,j1,k1,ik)=ig 
       enddo 
      enddo 
      !--
      allocate(LKGI(NTG,Nk_irr));LKGI=0.0d0 
      do ik=1,Nk_irr 
       do ig=1,NGI(ik) 
       ktmp(1)=(SKI(1,ik)+dble(KGI(1,ig,ik)))*b1(1)
     +        +(SKI(2,ik)+dble(KGI(2,ig,ik)))*b2(1) 
     +        +(SKI(3,ik)+dble(KGI(3,ig,ik)))*b3(1) 

       ktmp(2)=(SKI(1,ik)+dble(KGI(1,ig,ik)))*b1(2)
     +        +(SKI(2,ik)+dble(KGI(2,ig,ik)))*b2(2) 
     +        +(SKI(3,ik)+dble(KGI(3,ig,ik)))*b3(2) 

       ktmp(3)=(SKI(1,ik)+dble(KGI(1,ig,ik)))*b1(3)
     +        +(SKI(2,ik)+dble(KGI(2,ig,ik)))*b2(3) 
     +        +(SKI(3,ik)+dble(KGI(3,ig,ik)))*b3(3) 
       LKGI(ig,ik)=ktmp(1)**2+ktmp(2)**2+ktmp(3)**2
       enddo!ig 
       write(6,*) maxval(LKGI(:,ik))
      enddo!ik 
      write(6,*) 
      write(6,*) maxval(LKGI(:,:))
      Ecut_for_psi=maxval(LKGI(:,:))+1.0d-8
      write(6,*)'Ecut_for_psi=',Ecut_for_psi 
!--
!fft grid
      !--
      !m1=maxval(abs(KGI(1,:,:))) 
      !m2=maxval(abs(KGI(2,:,:))) 
      !m3=maxval(abs(KGI(3,:,:))) 
      !nwx2=algn235(2*m1) 
      !nwy2=algn235(2*m2) 
      !nwz2=algn235(2*m3) 
      !--
      htmp=dsqrt(dot_product(a1,a1))
      h1(:)=a1(:)/htmp 
      htmp=dsqrt(dot_product(a2,a2))
      h2(:)=a2(:)/htmp 
      htmp=dsqrt(dot_product(a3,a3)) 
      h3(:)=a3(:)/htmp 
      d1=abs(dot_product(b1,h1)) 
      d2=abs(dot_product(b2,h2)) 
      d3=abs(dot_product(b3,h3)) 
      qwf=2.0d0*dsqrt(Ecut_for_psi) 
      !--
      !nwx2=algn235(int(qwf/d1)+1,1)!wrong 
      !nwy2=algn235(int(qwf/d2)+1,1)!wrong 
      !nwz2=algn235(int(qwf/d3)+1,1)!wrong 
      !--
      nwx2=algn235(int(qwf/d1)+1)!20180412 by Y.Nomura 
      nwy2=algn235(int(qwf/d2)+1)!20180412 by Y.Nomura 
      nwz2=algn235(int(qwf/d3)+1)!20180412 by Y.Nomura 
      !--
      nfft1=nwx2+1;nfft2=nwy2+1;nfft3=nwz2+1
      Nl123=nfft1*nfft2*nfft3 
      write(6,'(a20,i10)')'nwx2:',nwx2 
      write(6,'(a20,i10)')'nwy2:',nwy2 
      write(6,'(a20,i10)')'nwz2:',nwz2 
      write(6,'(a20,i10)')'nfft1:',nfft1 
      write(6,'(a20,i10)')'nfft2:',nfft2 
      write(6,'(a20,i10)')'nfft3:',nfft3 
      write(6,'(a20,i10)')'NL123:',Nl123  
!--
!call fft3_init(nwx2,nwy2,nwz2,nfft1,nfft2,nfft3,fs) 
!--
!OPEN(111,R,FILE='dat.eigenvalue') 
      OPEN(111,FILE='./dir-wfn/dat.eigenvalue') 
      rewind(111)
      read(111,*) NBAND 
      allocate(E_EIGI(NBAND,Nk_irr));E_EIGI=0.0d0
      do ik=1,Nk_irr 
       do ib=1,NBAND 
        read(111,*) E_EIGI(ib,ik)
       enddo!ib
      enddo!ik          
      !--
      call calc_NB_start_end(NK_irr,NBAND,E_EIGI(1,1),E_LOWER,E_UPPER,
     +     NB_start,NB_end) 
      !--
      NTB_glb=NBAND 
      NTB=NB_end-NB_start+1!TOTAL NUMBER OF CONSIDERED BAND in wannier  
      write(6,'(a20,i10)')'NTB_glb:',NTB_glb  
      write(6,'(a20,i10)')'NB_start:',NB_start 
      write(6,'(a20,i10)')'NB_end:',NB_end 
      write(6,'(a20,i10)')'NTB:',NTB 
      deallocate(E_EIGI) 
!--
!OPEN(111,R,FILE='dat.eigenvalue') 
      OPEN(111,FILE='./dir-wfn/dat.eigenvalue') 
      rewind(111)
      read(111,*) NBAND 
      allocate(E_EIGI(NTB,Nk_irr));E_EIGI=0.0d0
      do ik=1,Nk_irr 
       do ib=1,NBAND!NBI(ik) 
        if((ib.ge.NB_start).and.(ib.le.NB_end)) then  
         read(111,*) E_EIGI(ib-NB_start+1,ik)
        else 
         read(111,*) 
        endif 
       enddo!ib
      enddo!ik          
!--
      call est_NTK(Nk_irr,Nsymq,SKI(1,1),rg(1,1,1),NTK)
      write(6,'(a20,i10)')'Estimated NTK:',NTK 
!--
!SK0,numirr,numrot,trs,RW
      allocate(SK0(3,NTK));SK0(:,:)=0.0d0
      allocate(numirr(NTK));numirr(:)=0
      allocate(numrot(NTK));numrot(:)=0
      allocate(trs(NTK));trs(:)=0
      allocate(RW(3,NTK));RW(:,:)=0
!--
!      do ik=1,Nk_irr 
!       SK0(:,ik)=SKI(:,ik) 
!       numirr(ik)=ik;numrot(ik)=1;trs(ik)=1;RW(1:3,ik)=0
!      enddo 
!      jk=Nk_irr 
!--
      jk=0
      do ik=1,Nk_irr
       !write(6,*)'--' 
       !write(6,'(a8,i8,3f15.10)')'ik=',ik,SKI(:,ik)  
       !write(6,*)'--' 
      do iop=1,Nsymq
!sym
       ktmp(:)=0.0d0;RWtmp(:)=0  
       ktmp(1)=rg(1,1,iop)*SKI(1,ik)+rg(1,2,iop)*SKI(2,ik)
     +        +rg(1,3,iop)*SKI(3,ik)
       ktmp(2)=rg(2,1,iop)*SKI(1,ik)+rg(2,2,iop)*SKI(2,ik)
     +        +rg(2,3,iop)*SKI(3,ik)
       ktmp(3)=rg(3,1,iop)*SKI(1,ik)+rg(3,2,iop)*SKI(2,ik)
     +        +rg(3,3,iop)*SKI(3,ik)
       call kcheck(ktmp(1),RWtmp(1))!rewind check 
       do iik=1,jk
        if(abs(SK0(1,iik)-ktmp(1))<1.0d-4.and. 
     +     abs(SK0(2,iik)-ktmp(2))<1.0d-4.and. 
     +     abs(SK0(3,iik)-ktmp(3))<1.0d-4)goto 1000
       enddo!iik
       jk=jk+1
       !write(6,'(a8,i8,3f15.10)')'jk=',jk,ktmp(:)  
       SK0(:,jk)=ktmp(:)
       numirr(jk)=ik;numrot(jk)=iop;trs(jk)=1;RW(:,jk)=RWtmp(:)
       !if(.true.)goto 2000
!time-reversal
1000   ktmp(:)=0.0d0;RWtmp(:)=0  
       ktmp(1)=rg(1,1,iop)*SKI(1,ik)+rg(1,2,iop)*SKI(2,ik)
     +        +rg(1,3,iop)*SKI(3,ik)
       ktmp(2)=rg(2,1,iop)*SKI(1,ik)+rg(2,2,iop)*SKI(2,ik)
     +        +rg(2,3,iop)*SKI(3,ik)
       ktmp(3)=rg(3,1,iop)*SKI(1,ik)+rg(3,2,iop)*SKI(2,ik)
     +        +rg(3,3,iop)*SKI(3,ik)
       call kcheck_trs(ktmp(1),RWtmp(1))!rewind check modified 20170316  
       do iik=1,jk
        if(abs(SK0(1,iik)-(-ktmp(1)))<1.0d-4.and. 
     +     abs(SK0(2,iik)-(-ktmp(2)))<1.0d-4.and. 
     +     abs(SK0(3,iik)-(-ktmp(3)))<1.0d-4)goto 2000
       enddo!iik
       jk=jk+1
       !write(6,'(a8,i8,3f15.10)')'jk=',jk,ktmp(:)  
       SK0(:,jk)=-ktmp(:) 
       numirr(jk)=ik;numrot(jk)=iop;trs(jk)=-1;RW(:,jk)=RWtmp(:) 
2000  enddo 
      enddo 
      !--
      if(NTK/=jk) then 
       write(6,*)'ERROR;STOP;NTK should be jk'   
       write(6,*)'NTK=',NTK,'jk=',jk;STOP
      endif 
      !--
      write(6,*) 
      write(6,*)'====================='
      write(6,*)'WRITE SAMPLE K-POINTS'
      write(6,*)'====================='
      do ik=1,NTK
       WRITE(6,'(I5,3F15.10)') ik,SK0(:,ik) 
      enddo  
!--
!      do ik=1,NTK 
!       write(6,*) ik,numrot(ik) 
!      enddo 
!--
      call est_nkbi(NTK,SK0(1,1),nkb1,nkb2,nkb3)  
      write(6,*) 
      write(6,'(a24,4i10)')'nkb1,nkb2,nkb3,NTK=',nkb1,nkb2,nkb3,NTK  
      Na1=nkb1/2;Na2=nkb2/2;Na3=nkb3/2
      write(6,'(a24,3i10)')'Na1,Na2,Na3=',Na1,Na2,Na3 
!--
      call system('rm -rf dir-wan') 
      call system('mkdir dir-wan') 
!--
!write Fermi surface 
      if(CALC_FERMISURFACE)then 
        write(6,*) 
        write(6,'(a30)')'write Fermi surface'
        write(6,*) 
        allocate(E_EIG(NTB,NTK));E_EIG(:,:)=0.0d0
        do jk=1,NTK 
         if(trs(jk)==1) then 
          ik=numirr(jk) 
          E_EIG(:,jk)=E_EIGI(:,ik) 
         elseif(trs(jk)==-1)then  
          ik=numirr(jk) 
          E_EIG(:,jk)=E_EIGI(:,ik) 
         endif 
        enddo!jk  
        filename='./dir-wan/dat.frmsf'
        call wrt_frmsf(NTB,NTK,nkb1,nkb2,nkb3,E_EIG(1,1),SK0(1,1),
     +       FermiEnergy,filename,b1(1),b2(1),b3(1)) 
        deallocate(E_EIG)
      endif 
!-- 
!write global DOS
      if(CALC_GLOBAL_DOS)then 
        write(6,*) 
        write(6,'(a30)')'write global DOS'
        write(6,*) 
        allocate(E_EIG(NTB,NTK));E_EIG(:,:)=0.0d0
        do jk=1,NTK 
         if(trs(jk)==1) then 
          ik=numirr(jk) 
          E_EIG(:,jk)=E_EIGI(:,ik) 
         elseif(trs(jk)==-1)then  
          ik=numirr(jk) 
          E_EIG(:,jk)=E_EIGI(:,ik) 
         endif 
        enddo!jk  
        filename='./dir-wan/dat.dos.global'
        call calculate_dos(NTB,NTK,nkb1,nkb2,nkb3,0.0d0,FermiEnergy,
     +       filename,b1(1),b2(1),b3(1),SK0(1,1),E_EIG(1,1))
        deallocate(E_EIG)
      endif 
!--
!OPEN(102,R,FILE='dat.wfn',FORM='unformatted') 
      OPEN(102,FILE='./dir-wfn/dat.wfn',FORM='unformatted') 
!old--
!      allocate(CIR(NTG,NTB,Nk_irr));CIR=0.0d0  
!      rewind(102)
!      do ik=1,Nk_irr 
!       do ib=1,NBAND!NBI(ik) 
!        if((ib.ge.NB_start).and.(ib.le.NB_end)) then  
!         read(102)(CIR(IG,ib-NB_start+1,ik),IG=1,NGI(ik))
!        else 
!         read(102) 
!        endif 
!       enddo!ib 
!      enddo!ik          
!      close(102) 
!new--
      rewind(102)
      read(102) ncomp 
      if(ncomp/=1)then 
       write(6,*)'This program not suport ncomp/=1; then stop'
       stop
      endif 
      allocate(CIR(NTG,NTB,Nk_irr));CIR=0.0d0  
      do ik=1,Nk_irr 
       do ib=1,NBAND 
        if((ib.ge.NB_start).and.(ib.le.NB_end)) then  
         read(102)(CIR(IG,ib-NB_start+1,ik),IG=1,NGI(ik))
        else 
         read(102) 
        endif 
       enddo!ib 
      enddo!ik          
      close(102) 
!--
!gen(C0,E_EIG,NG0,NB0,KG0) 
      allocate(C0(NTG,NTB,NTK));C0(:,:,:)=0.0d0
      allocate(E_EIG(NTB,NTK));E_EIG(:,:)=0.0d0
      allocate(NB0(NTK));NB0(:)=0
      allocate(NG0(NTK));NG0(:)=0
      allocate(KG0(3,NTG,NTK));KG0(:,:,:)=0 
      allocate(KGtmp(3,NTG));KGtmp(:,:)=0 
!20161207
!      do jk=1,Nk_irr 
!       C0(:,:,jk)=CIR(:,:,jk) 
!       E_EIG(:,jk)=E_EIGI(:,jk) 
!       NB0(jk)=NBI(jk)
!       NG0(jk)=NGI(jk)
!       KG0(:,:,jk)=KGI(:,:,jk)
!      enddo 
!      do jk=Nk_irr+1,NTK 
!--
       do jk=1,NTK 
       if(trs(jk)==1) then 
        ik=numirr(jk);iop=numrot(jk) 
        ktmp(1)=rg(1,1,iop)*SKI(1,ik)+rg(1,2,iop)*SKI(2,ik)
     +         +rg(1,3,iop)*SKI(3,ik)+dble(RW(1,jk)) 
        ktmp(2)=rg(2,1,iop)*SKI(1,ik)+rg(2,2,iop)*SKI(2,ik)
     +         +rg(2,3,iop)*SKI(3,ik)+dble(RW(2,jk))  
        ktmp(3)=rg(3,1,iop)*SKI(1,ik)+rg(3,2,iop)*SKI(2,ik)
     +         +rg(3,3,iop)*SKI(3,ik)+dble(RW(3,jk))  
        E_EIG(:,jk)=E_EIGI(:,ik) 
        NB0(jk)=NBAND!NBI(ik)
!--
!20170316 KG0s are generated by original KGI
!        NG0(jk)=NGI(ik) 
!        do ig=1,NG0(jk)  
!         i1=rg(1,1,iop)*KGI(1,ig,ik)
!     +     +rg(1,2,iop)*KGI(2,ig,ik)
!     +     +rg(1,3,iop)*KGI(3,ig,ik)-RW(1,jk) 
!         i2=rg(2,1,iop)*KGI(1,ig,ik)
!     +     +rg(2,2,iop)*KGI(2,ig,ik)
!     +     +rg(2,3,iop)*KGI(3,ig,ik)-RW(2,jk) 
!         i3=rg(3,1,iop)*KGI(1,ig,ik)
!     +     +rg(3,2,iop)*KGI(2,ig,ik)
!     +     +rg(3,3,iop)*KGI(3,ig,ik)-RW(3,jk) 
!         KG0(1,ig,jk)=i1 
!         KG0(2,ig,jk)=i2 
!         KG0(3,ig,jk)=i3 
!         phase=tpi*(dble(i1)*dble(pg(1,iop))
!     +             +dble(i2)*dble(pg(2,iop))
!     +             +dble(i3)*dble(pg(3,iop))) 
!         pf=exp(-ci*phase/dble(nnp)) 
!         C0(ig,:,jk)=CIR(ig,:,ik)*pf 
!        enddo!ig  
!       write(6,'(i5,3f15.10,i5)') jk,ktmp(:),trs(jk)  
!--
        call make_KG0(NTG,b1(1),b2(1),b3(1),Ecut_for_psi,
     +       ktmp(1),ktmp(2),ktmp(3),KG0(1,1,jk),NG_for_psi)
        if(NG_for_psi/=NGI(ik)) then 
         write(6,*)'ERROR; STOP; NG_for_psi should be NGI(ik)'   
         write(6,*)'NG_for_psi=',NG_for_psi,'NG0(ik)=',NGI(ik)
         write(6,*)'ik,jk',ik,jk;STOP
        endif 
        NG0(jk)=NG_for_psi  
        call make_C0(NTG,NTB,trs(jk),NG0(jk),KG0(1,1,jk),RW(1,jk),
     +  rginv(1,1,iop),pg(1,iop),nnp,L1,L2,L3,packing(-L1,-L2,-L3,ik),
     +  CIR(1,1,ik),C0(1,1,jk)) 
       !--
       elseif(trs(jk)==-1)then  
        ik=numirr(jk);iop=numrot(jk) 
        ktmp(1)=rg(1,1,iop)*SKI(1,ik)+rg(1,2,iop)*SKI(2,ik)
     +         +rg(1,3,iop)*SKI(3,ik)+dble(RW(1,jk)) 
        ktmp(2)=rg(2,1,iop)*SKI(1,ik)+rg(2,2,iop)*SKI(2,ik)
     +         +rg(2,3,iop)*SKI(3,ik)+dble(RW(2,jk))  
        ktmp(3)=rg(3,1,iop)*SKI(1,ik)+rg(3,2,iop)*SKI(2,ik)
     +         +rg(3,3,iop)*SKI(3,ik)+dble(RW(3,jk))  
        !write(6,'(i5,3f15.10,i5)') jk,ktmp(:),trs(jk)  
        E_EIG(:,jk)=E_EIGI(:,ik) 
        NB0(jk)=NBAND!NBI(ik)
        KGtmp(:,:)=0 
        call make_KG0(NTG,b1(1),b2(1),b3(1),Ecut_for_psi,
     +       ktmp(1),ktmp(2),ktmp(3),KGtmp(1,1),NG_for_psi)
        if(NG_for_psi/=NGI(ik))then 
         write(6,*)'ERROR; STOP; NG_for_psi should be NGI(ik)'   
         write(6,*)'NG_for_psi=',NG_for_psi,'NGI(ik)=',NGI(ik);STOP
        endif 
        NG0(jk)=NG_for_psi  
        KG0(:,:,jk)=-KGtmp(:,:)!notice on '-' sign 
        call make_C0(NTG,NTB,trs(jk),NG0(jk),KGtmp(1,1),RW(1,jk), 
     +  rginv(1,1,iop),pg(1,iop),nnp,L1,L2,L3,packing(-L1,-L2,-L3,ik),
     +  CIR(1,1,ik),C0(1,1,jk))!modified 20170316 KG0->KGtmp 
       endif 
      enddo 
      !--
      WRITE(6,*) 
      write(6,*)'==================='
      write(6,*)'CALCULATED K-POINTS'
      write(6,*)'==================='
      WRITE(6,*) 
      do ik=1,NTK
       write(6,'(i5,3f15.10,i5)') ik,SK0(:,ik),trs(ik)  
      enddo 
!--
!     WRITE(6,*) 
!     write(6,*)'====================='
!     write(6,*)'WRITE ALL EIGENVALUES'
!     write(6,*)'====================='
!     WRITE(6,*) 
!     do ik=1,NTK
!      write(6,'(5F15.10)')(E_EIG(ib,ik),ib=1,NTB) 
!     enddo 
!--
!
!=========================================================
!=== CALCULATE THE NUMBER OF BANDS INSIDE OUTER WINDOW ===
!=========================================================

      allocate(N_BAND(NTK));N_BAND(:)=0
      allocate(N_BAND_BTM(NTK));N_BAND_BTM(:)=0
      do ik=1,NTK       
       SUM_INT=0        
       do ib=1,NTB        
        IF(E_LOWER<=E_EIG(ib,ik).AND.E_EIG(ib,ik)<=E_UPPER)THEN    
         SUM_INT=SUM_INT+1         
         IF(SUM_INT==1) N_BAND_BTM(ik)=ib-1  
        ENDIF        
       enddo!ib 
       N_BAND(ik)=SUM_INT 
      enddo!ik       
      !
      write(6,*) 
      write(6,*)'===================================='
      write(6,*)'BAND INFO INSIDE OUTER ENERGY WINDOW'
      write(6,*)'===================================='
      write(6,*) 
      !
      if(maxval(N_BAND)/=NTB)then
       write(6,'(a)')'maxval(N_BAND) SHOULD BE THE SAME AS NTB; STOP'
      endif 
      !
      do ik=1,NTK
!check N_BAND(ik) >= N_wannier 
       if(N_BAND(ik)<N_wannier)then 
        write(6,'(a)')'energy window setting: WRONG'
        write(6,'(a)')'N_BAND(ik) SHOULD BE LARGER THAN N_wannier: STOP'
        write(6,'(a,3i8)')'ik, N_BAND(ik), N_wannier',
     +                     ik, N_BAND(ik), N_wannier 
        stop 
       endif  
       write(6,*) ik,N_BAND_BTM(ik),N_BAND(ik) 
      enddo!ik 
      write(6,*) 
!--
!OPEN(149,W,FILE='dat.ns-nb')
      OPEN(149,FILE='./dir-wan/dat.ns-nb') 
      rewind(149)
      do ik=1,NTK
       write(149,*) N_BAND_BTM(ik)+NB_start-1,N_BAND(ik) 
      enddo!ik 
!--
!20180906 
      if(CALC_REAL_SPACE_BLOCH)then 
       write(6,*) 
       write(6,*)'================================'
       write(6,*)'=== VISUALIZE BLOCH FUNCTION ==='
       write(6,*)'================================'
       write(6,*) 
       !--
       nfft1=nwx2+1;nfft2=nwy2+1;nfft3=nwz2+1
       Nl123=nfft1*nfft2*nfft3 
       write(6,'(a20,3i10)')'nwx2,nwy2,nwz2:',nwx2,nwy2,nwz2 
       write(6,'(a20,3i10)')'nfft1,nfft2,nfft2:',nfft1,nfft2,nfft3 
       write(6,'(a20,i10)')'Nl123:',Nl123 
       call fft3_init(nwx2,nwy2,nwz2,nfft1,nfft2,nfft3,fs) 
       allocate(fftwk(Nl123*2),stat=err) 
       allocate(wfunc(Nl123*2),stat=err) 
       !--
       write(6,'(a20,3f15.10)')'calc_k:',calc_k
       call search_ik(NTK,SK0(1,1),calc_k(1),ik) 
       write(6,'(a20,i10)')'ik:',ik 
       !
       amax=1; bmax=1; cmax=1 
       if(CALC_BLOCH_VIS_RANGE)then 
        write(6,*) 
        write(6,'(a40)')'CALCULATE BLOCH VISUALIZATION RANGE'
        write(6,*) 
        call search_vis_range(calc_k(1),amin,amax,bmin,bmax,cmin,cmax) 
       endif  
       !
       write(6,'(a20,3i10)')'supercell:',amax,bmax,cmax 
       allocate(BF_REALSPACE((amin)*nwx2:(amax)*nwx2-1,
     +                       (bmin)*nwy2:(bmax)*nwy2-1, 
     +                       (cmin)*nwz2:(cmax)*nwz2-1))!<--complex(8)
       !-- 
       do ib=1,N_BAND(ik)
        BF_REALSPACE=0.0d0           
        call make_pw(C0(1,N_BAND_BTM(ik)+ib,ik),wfunc(1),fftwk(1),
     +   NG0(ik),KG0(1,1,ik),NTG,nwx2,nwy2,nwz2,nfft1,nfft2,Nl123,fs,
     +   amin,amax,bmin,bmax,cmin,cmax,SK0(1,ik),
     +   BF_REALSPACE(amin*nwx2,bmin*nwy2,cmin*nwz2)) 
        call phase_fix_bf(nwx2,nwy2,nwz2,amin,amax,bmin,bmax,
     +   cmin,cmax,BF_REALSPACE(amin*nwx2,bmin*nwy2,cmin*nwz2)) 
        !--
        !normalize 
        !BF_REALSPACE=BF_REALSPACE/DBLE(NTK)/DSQRT(VOLUME)
        !--
        !OUTPUT by VESTA format 
        !OPEN(116,W,FILE='dat.bf-realspace-xxx.grd')
        write(filename,"('dat.bf-realspace-',i3.3,'.grd')")ib 
        OPEN(116,FILE=filename) 
        REWIND(116)
        na_grid=nwx2 
        nb_grid=nwy2
        nc_grid=nwz2 
        aa1=a1*dble(amax-amin)
        aa2=a2*dble(bmax-bmin)
        aa3=a3*dble(cmax-cmin)
        call est_latparam(aa1(1),aa2(1),aa3(1),a,b,c,alp,bet,gmm) 
        write(116,'(a40,f10.5)')'Electron density and energy (eV)',
     +  E_EIG(N_BAND_BTM(ik)+ib,ik)*au 
        write(116,'(6f15.10)') a*bohr,b*bohr,c*bohr,alp,bet,gmm  
        write(116,'(3I5)')
     +  (amax-amin)*na_grid,(bmax-bmin)*nb_grid,(cmax-cmin)*nc_grid 
        write(116,'(6g25.16)')(((dble(BF_REALSPACE(i1,i2,i3)),
     +  i3=(cmin)*nc_grid,(cmax)*nc_grid-1),
     +  i2=(bmin)*nb_grid,(bmax)*nb_grid-1),
     +  i1=(amin)*na_grid,(amax)*na_grid-1)
        CLOSE(116)
       enddo!ib 
       deallocate(BF_REALSPACE,fftwk,wfunc)
       !--
       !gencif 20180906  
       call gencif(a1(1),a2(1),a3(1),amin,amax,bmin,bmax,cmin,cmax,
     +      nsi,chemical_species(1),asi(1,1))
       !--
       write(6,*) 
       write(6,*)'==================================='
       write(6,*)'FINIDH VISUALIZE BLOCH FUNCTION ==='
       write(6,*)'!!!!!! We stop calculation !!!!!!!!'
       write(6,*)'==================================='
       write(6,*) 
       stop
      endif!CALC_REAL_SPACE_BLOCH 
      !--
      WRITE(6,*) 
      WRITE(6,*)'========================================='
      WRITE(6,*)'=== WE NOW CALCULATE WANNIER FUNCTION ==='
      WRITE(6,*)'========================================='
      WRITE(6,*) 

!================================
!=== SET B-VECTORS AND WEIGHT ===
!================================

      select case(icell) 
      case(0)!=== CASE(0) ===      
!
!== Automatic generation of b-vectors by TERUMASA TADANO ==
!
       allocate(b_LAT(3,NBMAX));b_LAT(:,:)=0.0D0 
       allocate(VEC_b(3,NBMAX));VEC_b(:,:)=0.0d0                 
       allocate(WEIGHT_b(NBMAX));WEIGHT_b(:)=0.0D0             
       call generate_bvectors(b1(1),b2(1),b3(1),nkb1,nkb2,nkb3, 
     + NBMAX,maxshell,verbosity,NB,b_LAT(1,1),VEC_b(1,1),WEIGHT_b(1))
      !--
      case(1)!=== CASE(1) SIMPLE CUBIC CELL ===      
       NB=6
       allocate(b_LAT(3,NB));b_LAT(:,:)=0.0D0 
       allocate(VEC_b(3,NB));VEC_b(:,:)=0.0d0                 
       allocate(WEIGHT_b(NB));WEIGHT_b(:)=0.0D0             
       b_LAT(:,:) = 0.0D0 
       b_LAT(1,1) =  (1.0D0/DBLE(nkb1))
       b_LAT(2,2) =  (1.0D0/DBLE(nkb1))
       b_LAT(3,3) =  (1.0D0/DBLE(nkb1))
       b_LAT(1,4) = -(1.0D0/DBLE(nkb1))
       b_LAT(2,5) = -(1.0D0/DBLE(nkb1))
       b_LAT(3,6) = -(1.0D0/DBLE(nkb1))
       VEC_b(:,:) = 0.0D0                            
       VEC_b(:,1) =  (1.0D0/DBLE(nkb1)) * b1(:) 
       VEC_b(:,2) =  (1.0D0/DBLE(nkb1)) * b2(:) 
       VEC_b(:,3) =  (1.0D0/DBLE(nkb1)) * b3(:) 
       VEC_b(:,4) = -(1.0D0/DBLE(nkb1)) * b1(:)
       VEC_b(:,5) = -(1.0D0/DBLE(nkb1)) * b2(:) 
       VEC_b(:,6) = -(1.0D0/DBLE(nkb1)) * b3(:) 
       VEC_b2=VEC_b(1,1)**2+VEC_b(2,1)**2+VEC_b(3,1)**2 
       WEIGHT_b(:) = 3.0d0/dble(NB)/VEC_b2              
       !WEIGHT_b(:) = (3.0D0/6.0D0)*(alat*nkb1/tpi)**2 
      !--
      case(2)!=== CASE(2) FCC PRIMITIVE CELL ===      
       NB=8
       allocate(b_LAT(3,NB));b_LAT(:,:)=0.0D0 
       allocate(VEC_b(3,NB));VEC_b(:,:)=0.0d0                 
       allocate(WEIGHT_b(NB));WEIGHT_b(:)=0.0D0             
       b_LAT(:,:) = 0.0D0 
       b_LAT(1,1) =  (1.0D0/DBLE(nkb1))
       b_LAT(2,2) =  (1.0D0/DBLE(nkb1))
       b_LAT(3,3) =  (1.0D0/DBLE(nkb1))
       b_LAT(1,4) =  (1.0D0/DBLE(nkb1))
       b_LAT(2,4) =  (1.0D0/DBLE(nkb1))
       b_LAT(3,4) =  (1.0D0/DBLE(nkb1))
       b_LAT(1,5) = -(1.0D0/DBLE(nkb1))
       b_LAT(2,6) = -(1.0D0/DBLE(nkb1))
       b_LAT(3,7) = -(1.0D0/DBLE(nkb1))
       b_LAT(1,8) = -(1.0D0/DBLE(nkb1))
       b_LAT(2,8) = -(1.0D0/DBLE(nkb1))
       b_LAT(3,8) = -(1.0D0/DBLE(nkb1))
       VEC_b(:,:) = 0.0D0                            
       VEC_b(:,1) =  (1.0D0/DBLE(nkb1)) * b1(:) 
       VEC_b(:,2) =  (1.0D0/DBLE(nkb1)) * b2(:) 
       VEC_b(:,3) =  (1.0D0/DBLE(nkb1)) * b3(:) 
       VEC_b(:,4) =  (1.0D0/DBLE(nkb1)) *(b1(:)+b2(:)+b3(:))  
       VEC_b(:,5) = -(1.0D0/DBLE(nkb1)) * b1(:) 
       VEC_b(:,6) = -(1.0D0/DBLE(nkb1)) * b2(:) 
       VEC_b(:,7) = -(1.0D0/DBLE(nkb1)) * b3(:) 
       VEC_b(:,8) = -(1.0D0/DBLE(nkb1)) *(b1(:)+b2(:)+b3(:))  
       VEC_b2=VEC_b(1,1)**2+VEC_b(2,1)**2+VEC_b(3,1)**2 
       WEIGHT_b(:) = 3.0d0/dble(NB)/VEC_b2              
       !WEIGHT_b(:) = (1.0D0/8.0D0)*(alat*nkb1/tpi)**2 
      !--
      case(3)!=== CASE(3) BCC PRIMITIVE CELL ===      
       NB=12 
       allocate(b_LAT(3,NB));b_LAT(:,:)=0.0D0 
       allocate(VEC_b(3,NB));VEC_b(:,:)=0.0d0                 
       allocate(WEIGHT_b(NB));WEIGHT_b(:)=0.0D0             
       b_LAT(:,:) = 0.0D0 
       b_LAT(1, 1) =  (1.0D0/DBLE(nkb1))
       b_LAT(2, 2) =  (1.0D0/DBLE(nkb1))
       b_LAT(3, 3) =  (1.0D0/DBLE(nkb1))
       b_LAT(1, 4) =  (1.0D0/DBLE(nkb1))
       b_LAT(2, 4) = -(1.0D0/DBLE(nkb1))
       b_LAT(2, 5) =  (1.0D0/DBLE(nkb1))
       b_LAT(3, 5) = -(1.0D0/DBLE(nkb1))
       b_LAT(3, 6) =  (1.0D0/DBLE(nkb1))
       b_LAT(1, 6) = -(1.0D0/DBLE(nkb1))
       b_LAT(1, 7) = -(1.0D0/DBLE(nkb1))
       b_LAT(2, 8) = -(1.0D0/DBLE(nkb1))
       b_LAT(3, 9) = -(1.0D0/DBLE(nkb1))
       b_LAT(1,10) = -(1.0D0/DBLE(nkb1))
       b_LAT(2,10) =  (1.0D0/DBLE(nkb1))
       b_LAT(2,11) = -(1.0D0/DBLE(nkb1))
       b_LAT(3,11) =  (1.0D0/DBLE(nkb1))
       b_LAT(3,12) = -(1.0D0/DBLE(nkb1))
       b_LAT(1,12) =  (1.0D0/DBLE(nkb1))
       VEC_b(:,:) = 0.0D0                            
       VEC_b(:, 1) =  (1.0D0/DBLE(nkb1)) * b1(:) 
       VEC_b(:, 2) =  (1.0D0/DBLE(nkb1)) * b2(:) 
       VEC_b(:, 3) =  (1.0D0/DBLE(nkb1)) * b3(:) 
       VEC_b(:, 4) =  (1.0D0/DBLE(nkb1)) *(b1(:)-b2(:))      
       VEC_b(:, 5) =  (1.0D0/DBLE(nkb1)) *(b2(:)-b3(:))      
       VEC_b(:, 6) =  (1.0D0/DBLE(nkb1)) *(b3(:)-b1(:))      
       VEC_b(:, 7) = -(1.0D0/DBLE(nkb1)) * b1(:) 
       VEC_b(:, 8) = -(1.0D0/DBLE(nkb1)) * b2(:) 
       VEC_b(:, 9) = -(1.0D0/DBLE(nkb1)) * b3(:) 
       VEC_b(:,10) = -(1.0D0/DBLE(nkb1)) *(b1(:)-b2(:))   
       VEC_b(:,11) = -(1.0D0/DBLE(nkb1)) *(b2(:)-b3(:))   
       VEC_b(:,12) = -(1.0D0/DBLE(nkb1)) *(b3(:)-b1(:))   
       VEC_b2=VEC_b(1,1)**2+VEC_b(2,1)**2+VEC_b(3,1)**2 
       WEIGHT_b(:) = 3.0d0/dble(NB)/VEC_b2              
       !WEIGHT_b(:) = (1.0D0/8.0D0)*(alat*nkb1/tpi)**2 
      !--
      case(4)!=== CASE(4) HCP PRIMITIVE CELL ===      
       NB=8 
       allocate(b_LAT(3,NB));b_LAT(:,:)=0.0D0 
       allocate(VEC_b(3,NB));VEC_b(:,:)=0.0d0                 
       allocate(WEIGHT_b(NB));WEIGHT_b(:)=0.0D0             
       b_LAT(:,:) = 0.0D0 
       b_LAT(1,1) =  (1.0D0/DBLE(nkb1))
       b_LAT(2,2) =  (1.0D0/DBLE(nkb2))
       b_LAT(1,3) =  (1.0D0/DBLE(nkb1))
       b_LAT(2,3) = -(1.0D0/DBLE(nkb2))!SiO2(hcp)
       !b_LAT(2,3) =  (1.0D0/DBLE(nkb2))!Re(hcp)
       b_LAT(1,4) = -(1.0D0/DBLE(nkb1))
       b_LAT(2,5) = -(1.0D0/DBLE(nkb2))
       b_LAT(1,6) = -(1.0D0/DBLE(nkb1))
       b_LAT(2,6) =  (1.0D0/DBLE(nkb2))!SiO2(hcp) 
       !b_LAT(2,6) = -(1.0D0/DBLE(nkb2))!Re(hcp)
       b_LAT(3,7) =  (1.0D0/DBLE(nkb3))
       b_LAT(3,8) = -(1.0D0/DBLE(nkb3))
       VEC_b(:,:) = 0.0D0                            
       VEC_b(:,1) =  (1.0D0/DBLE(nkb1)) * b1(:) 
       VEC_b(:,2) =  (1.0D0/DBLE(nkb2)) * b2(:) 
       VEC_b(:,3) =  (1.0D0/DBLE(nkb1)) * b1(:) 
     +            -  (1.0D0/DBLE(nkb2)) * b2(:)!SiO2(hcp) 
     !+            +  (1.0D0/DBLE(nkb2)) * b2(:)!Re(hcp) 
       VEC_b(:,4) = -(1.0D0/DBLE(nkb1)) * b1(:) 
       VEC_b(:,5) = -(1.0D0/DBLE(nkb2)) * b2(:) 
       VEC_b(:,6) = -(1.0D0/DBLE(nkb1)) * b1(:) 
     +            +  (1.0D0/DBLE(nkb2)) * b2(:)!SiO2(hcp)  
     !+            -  (1.0D0/DBLE(nkb2)) * b2(:)!Re(hcp)
       VEC_b(:,7) =  (1.0D0/DBLE(nkb3)) * b3(:) 
       VEC_b(:,8) = -(1.0D0/DBLE(nkb3)) * b3(:)            
       La = DSQRT( a1(1)**2 + a1(2)**2 + a1(3)**2 )     
       Lc = DSQRT( a3(1)**2 + a3(2)**2 + a3(3)**2 )     
       WEIGHT_b(:) = 0.0D0             
       WEIGHT_b(1) = (1.0D0/4.0D0)*(La*nkb1/tpi)**2 
       WEIGHT_b(2) = (1.0D0/4.0D0)*(La*nkb1/tpi)**2 
       WEIGHT_b(3) = (1.0D0/4.0D0)*(La*nkb1/tpi)**2 
       WEIGHT_b(4) = (1.0D0/4.0D0)*(La*nkb1/tpi)**2 
       WEIGHT_b(5) = (1.0D0/4.0D0)*(La*nkb1/tpi)**2 
       WEIGHT_b(6) = (1.0D0/4.0D0)*(La*nkb1/tpi)**2 
       WEIGHT_b(7) = (1.0D0/2.0D0)*(Lc*nkb3/tpi)**2 
       WEIGHT_b(8) = (1.0D0/2.0D0)*(Lc*nkb3/tpi)**2 
      !--
      case(5)!=== CASE(5) ORTHORHOMBIC OR TETRAGONAL CELL ===
       NB=6 
       allocate(b_LAT(3,NB));b_LAT(:,:)=0.0D0 
       allocate(VEC_b(3,NB));VEC_b(:,:)=0.0d0                 
       allocate(WEIGHT_b(NB));WEIGHT_b(:)=0.0D0             
       b_LAT(:,:) = 0.0D0 
       b_LAT(1,1) =  (1.0D0/DBLE(nkb1))
       b_LAT(2,2) =  (1.0D0/DBLE(nkb2))
       b_LAT(3,3) =  (1.0D0/DBLE(nkb3))
       b_LAT(1,4) = -(1.0D0/DBLE(nkb1))
       b_LAT(2,5) = -(1.0D0/DBLE(nkb2))
       b_LAT(3,6) = -(1.0D0/DBLE(nkb3))
       VEC_b(:,:) = 0.0D0                            
       VEC_b(:,1) =  (1.0D0/DBLE(nkb1)) * b1(:) 
       VEC_b(:,2) =  (1.0D0/DBLE(nkb2)) * b2(:) 
       VEC_b(:,3) =  (1.0D0/DBLE(nkb3)) * b3(:) 
       VEC_b(:,4) = -(1.0D0/DBLE(nkb1)) * b1(:)
       VEC_b(:,5) = -(1.0D0/DBLE(nkb2)) * b2(:) 
       VEC_b(:,6) = -(1.0D0/DBLE(nkb3)) * b3(:) 
       La = DSQRT( a1(1)**2 + a1(2)**2 + a1(3)**2 )     
       Lb = DSQRT( a2(1)**2 + a2(2)**2 + a2(3)**2 )     
       Lc = DSQRT( a3(1)**2 + a3(2)**2 + a3(3)**2 )     
       WEIGHT_b(:) = 0.0D0             
       WEIGHT_b(1) = (1.0D0/2.0D0)*(La*nkb1/tpi)**2 
       WEIGHT_b(2) = (1.0D0/2.0D0)*(Lb*nkb2/tpi)**2 
       WEIGHT_b(3) = (1.0D0/2.0D0)*(Lc*nkb3/tpi)**2 
       WEIGHT_b(4) = (1.0D0/2.0D0)*(La*nkb1/tpi)**2 
       WEIGHT_b(5) = (1.0D0/2.0D0)*(Lb*nkb2/tpi)**2 
       WEIGHT_b(6) = (1.0D0/2.0D0)*(Lc*nkb3/tpi)**2 
      !--
      case(6)!=== CASE(6) MONOCLINIC CELL === 
       NB=8
       allocate(b_LAT(3,NB));b_LAT(:,:)=0.0D0 
       allocate(VEC_b(3,NB));VEC_b(:,:)=0.0d0                 
       allocate(WEIGHT_b(NB));WEIGHT_b(:)=0.0D0             
       b_LAT(:,:) = 0.0D0 
       b_LAT(1,1) =  (1.0D0/DBLE(nkb1))
       b_LAT(3,2) =  (1.0D0/DBLE(nkb3))
       b_LAT(1,3) =  (1.0D0/DBLE(nkb1))
       b_LAT(3,3) = -(1.0D0/DBLE(nkb3))
       b_LAT(1,4) = -(1.0D0/DBLE(nkb1))
       b_LAT(3,5) = -(1.0D0/DBLE(nkb3))
       b_LAT(1,6) = -(1.0D0/DBLE(nkb1))
       b_LAT(3,6) =  (1.0D0/DBLE(nkb3))
       b_LAT(2,7) =  (1.0D0/DBLE(nkb2))
       b_LAT(2,8) = -(1.0D0/DBLE(nkb2))
       VEC_b(:,:) = 0.0D0                            
       VEC_b(:,1) = (1.0D0/DBLE(nkb1))*b1(:) 
       VEC_b(:,2) = (1.0D0/DBLE(nkb3))*b3(:) 
       VEC_b(:,3) = VEC_b(:,1)-VEC_b(:,2)
       VEC_b(:,4) =-VEC_b(:,1)
       VEC_b(:,5) =-VEC_b(:,2)
       VEC_b(:,6) =-VEC_b(:,3)
       VEC_b(:,7) = (1.0D0/DBLE(nkb2))*b2(:) 
       VEC_b(:,8) =-VEC_b(:,7)
       s = VEC_b(1,1)
       t = VEC_b(3,1)
       u = VEC_b(3,2)
       v = VEC_b(2,7) 
       WEIGHT_b(:) = 0.0D0             
       WEIGHT_b(1) = (u-t)/(2.0d0*(s**2)*u) 
       WEIGHT_b(2) = (s**2-t*(u-t))/(2.0d0*(s**2)*(u**2))
       WEIGHT_b(3) = t/(2.0d0*(s**2)*u) 
       WEIGHT_b(4) = WEIGHT_b(1)
       WEIGHT_b(5) = WEIGHT_b(2)
       WEIGHT_b(6) = WEIGHT_b(3)
       WEIGHT_b(7) = 1.0d0/(2.0d0*(v**2)) 
       WEIGHT_b(8) = WEIGHT_b(7) 
      !--
      case(8)!=== CASE(8) BCT CELL === 2011/1/4 by Y.Nomura ===
       write(6,*) 
       write(6,*)'==============================='
       write(6,*)' body centered tetragonal cell '
       write(6,*)'==============================='
       write(6,*) 
       NB=10
       allocate(b_LAT(3,NB));b_LAT(:,:)=0.0D0 
       allocate(VEC_b(3,NB));VEC_b(:,:)=0.0d0                 
       allocate(WEIGHT_b(NB));WEIGHT_b(:)=0.0D0             
       if(nkb1/=nkb2)stop'nkb1/=nkb2'
       if(nkb1/=nkb3)stop'nkb1/=nkb3'
       if(nkb2/=nkb3)stop'nkb2/=nkb3'
       b_LAT(:,:)=0.0D0 
       b_LAT(1,1)=(1.0D0/DBLE(nkb1))
       b_LAT(2,1)=(1.0D0/DBLE(nkb2))
       b_LAT(3,1)=-(1.0D0/DBLE(nkb3))
       b_LAT(1,2)=(1.0D0/DBLE(nkb1))
       b_LAT(2,3)=(1.0D0/DBLE(nkb2))
       b_LAT(1,4)=(1.0D0/DBLE(nkb1))
       b_LAT(3,4)=-(1.0D0/DBLE(nkb3))
       b_LAT(2,5)=(1.0D0/DBLE(nkb2))
       b_LAT(3,5)=-(1.0D0/DBLE(nkb3))
       b_LAT(:,6)=-b_LAT(:,1)
       b_LAT(:,7)=-b_LAT(:,2)
       b_LAT(:,8)=-b_LAT(:,3)
       b_LAT(:,9)=-b_LAT(:,4)
       b_LAT(:,10)=-b_LAT(:,5)
       VEC_b(:,:)=0.0D0                            
       VEC_b(:,1)=(1.0D0/DBLE(nkb1))*(b1(:)+b2(:)-b3(:))
       VEC_b(:,2)=(1.0D0/DBLE(nkb1))*b1(:) 
       VEC_b(:,3)=(1.0D0/DBLE(nkb2))*b2(:) 
       VEC_b(:,4)=(1.0D0/DBLE(nkb1))*b1(:)-(1.0D0/DBLE(nkb3))*b3(:)
       VEC_b(:,5)=(1.0D0/DBLE(nkb2))*b2(:)-(1.0D0/DBLE(nkb3))*b3(:)
       VEC_b(:,6)=-VEC_b(:,1)
       VEC_b(:,7)=-VEC_b(:,2)
       VEC_b(:,8)=-VEC_b(:,3)
       VEC_b(:,9)=-VEC_b(:,4)
       VEC_b(:,10)=-VEC_b(:,5)
!--
!NUMPACK
!       NBh=NB/2 
!       allocate(VEC_d(6));VEC_d(:)=0.0d0                  
!       allocate(bb(NBh,6));bb(:,:)=0.0d0 
!       allocate(aa(NBh,6));aa(:,:)=0.0d0 
!       do ib=1,Nbh 
!        ab=0
!        do ix=1,3
!         do jx=ix,3
!          ab=ab+1
!          bb(ib,ab)=VEC_b(ix,ib)*VEC_b(jx,ib)
!         enddo 
!        enddo 
!       enddo 
!       aa=bb 
!!do i=1,NBh 
!! write(6,'(10f20.5)')(bb(i,j),j=1,6)!org 
!!enddo 
!       call inv_ge(NBh,6,aa(1,1)) 
!!do i=1,NBh 
!! write(6,'(10f20.5)')(aa(i,j),j=1,6)!inv
!!enddo 
!!write(6,*) 
!!write(6,*)'CHECK BB*AA=1'       
!!write(6,*) 
!!do i=1,NBh
!! do j=1,NBh 
!!  SUM_REAL=0.0d0
!!  do k=1,6 
!!   SUM_REAL=SUM_REAL+bb(i,k)*aa(j,k)
!!  enddo 
!!  write(6,'(i5,i5,f15.10)') i,j,SUM_REAL 
!! enddo
!!enddo 
!--
!LAPACK
       NBh=NB/2 
       allocate(VEC_d(6));VEC_d(:)=0.0d0!tmp                   
       allocate(aa(6,6));aa(:,:)=0.0d0!matrix A tmp
       allocate(bb(6,6));bb(:,:)=0.0d0!matrix A^-1 tmp
       do ib=1,Nbh 
        ab=0
        do ix=1,3
         do jx=ix,3
          ab=ab+1
          aa(ib,ab)=VEC_b(ix,ib)*VEC_b(jx,ib)
         enddo 
        enddo 
       enddo 
       !--
       !bb<-aa 
       do i=1,6
        do j=1,6
         bb(j,i)=aa(i,j) 
        enddo
       enddo 
!--
!do i=1,n 
! write(6,'(10f15.10)')(bb(i,j),j=1,n)!org 
!enddo 
!--
       call inv_ge_lapack(6,NBh,bb(1,1)) 
!--
!write(6,*) 
!do i=1,n 
! write(6,'(10f15.10)')(bb(i,j),j=1,n)!inv
!enddo 
!--
!      write(6,*) 
!      write(6,*)'CHECK BB*AA=1'       
!      write(6,*) 
!      do i=1,NBh
!       do j=1,NBh 
!        s=0.0d0
!        do k=1,n 
!         s=s+aa(i,k)*bb(j,k)
!        enddo 
!        write(6,'(i5,i5,f15.10)') i,j,s
!       enddo
!      enddo 
!--
       VEC_d(1)=0.5d0!xx
       VEC_d(2)=0.0d0!xy
       VEC_d(3)=0.0d0!xz
       VEC_d(4)=0.5d0!yy 
       VEC_d(5)=0.0d0!yz 
       VEC_d(6)=0.5d0!zz
       !--
       WEIGHT_b(:)=0.0D0             
       do ib=1,Nbh
        ab=0
        SUM_REAL=0.0d0 
        do ix=1,3
         do jx=ix,3
          ab=ab+1
          !SUM_REAL=SUM_REAL+VEC_d(ab)*aa(ib,ab)!NUMPACK
          SUM_REAL=SUM_REAL+VEC_d(ab)*bb(ib,ab)!LAPACK 
         enddo 
        enddo 
        WEIGHT_b(ib)=SUM_REAL 
       enddo 
       WEIGHT_b(6)=WEIGHT_b(1)
       WEIGHT_b(7)=WEIGHT_b(2)
       WEIGHT_b(8)=WEIGHT_b(3) 
       WEIGHT_b(9)=WEIGHT_b(4) 
       WEIGHT_b(10)=WEIGHT_b(5) 
      !--
      case(7)!=== CASE(7) TRICLINIC CELL ===* 
       NB=12 
       allocate(b_LAT(3,NB));b_LAT(:,:)=0.0D0 
       allocate(VEC_b(3,NB));VEC_b(:,:)=0.0d0                 
       allocate(WEIGHT_b(NB));WEIGHT_b(:)=0.0D0             
       NBh=NB/2 
       b_LAT(:,:) = 0.0D0 
       b_LAT(1,1) =  (1.0D0/DBLE(nkb1))
       b_LAT(2,2) =  (1.0D0/DBLE(nkb2))
       b_LAT(3,3) =  (1.0D0/DBLE(nkb3))
       b_LAT(1,4) =  (1.0D0/DBLE(nkb1))
       b_LAT(2,4) = -(1.0D0/DBLE(nkb2))
       b_LAT(2,5) =  (1.0D0/DBLE(nkb2))
       b_LAT(3,5) = -(1.0D0/DBLE(nkb3))
       b_LAT(3,6) =  (1.0D0/DBLE(nkb3))
       b_LAT(1,6) = -(1.0D0/DBLE(nkb1))
       b_LAT(:,7) = -b_LAT(:,1) 
       b_LAT(:,8) = -b_LAT(:,2) 
       b_LAT(:,9) = -b_LAT(:,3) 
       b_LAT(:,10)= -b_LAT(:,4) 
       b_LAT(:,11)= -b_LAT(:,5) 
       b_LAT(:,12)= -b_LAT(:,6) 
       VEC_b(:,:) = 0.0D0                            
       VEC_b(:,1) = (1.0D0/DBLE(nkb1))*b1(:) 
       VEC_b(:,2) = (1.0D0/DBLE(nkb2))*b2(:) 
       VEC_b(:,3) = (1.0D0/DBLE(nkb3))*b3(:) 
       VEC_b(:,4) = VEC_b(:,1)-VEC_b(:,2)
       VEC_b(:,5) = VEC_b(:,2)-VEC_b(:,3)
       VEC_b(:,6) = VEC_b(:,3)-VEC_b(:,1)
       VEC_b(:,7) =-VEC_b(:,1)
       VEC_b(:,8) =-VEC_b(:,2)
       VEC_b(:,9) =-VEC_b(:,3)
       VEC_b(:,10)=-VEC_b(:,4)
       VEC_b(:,11)=-VEC_b(:,5)
       VEC_b(:,12)=-VEC_b(:,6)
       allocate(VEC_d(NBh));VEC_d(:)=0.0d0                  
       allocate(bb(NBh,NBh));bb(:,:)=0.0d0 
       do ib=1,Nbh 
        ab=0
        do ix=1,3
         do jx=ix,3
          ab=ab+1
          bb(ib,ab)=VEC_b(ix,ib)*VEC_b(jx,ib)
         enddo 
        enddo 
       enddo 
       call invmat(Nbh,bb(1,1)) 
       VEC_d(1)=0.5d0 !xx
       VEC_d(2)=0.0d0 !xy
       VEC_d(3)=0.0d0 !xz
       VEC_d(4)=0.5d0 !yy 
       VEC_d(5)=0.0d0 !yz 
       VEC_d(6)=0.5d0 !zz
       WEIGHT_b(:)=0.0D0             
       do ib=1,Nbh
        ab=0
        SUM_REAL=0.0d0 
        do ix=1,3
         do jx=ix,3
          ab=ab+1
          SUM_REAL=SUM_REAL+VEC_d(ab)*bb(ab,ib)
         enddo 
        enddo 
        WEIGHT_b(ib)=SUM_REAL
       enddo 
       WEIGHT_b(7) = WEIGHT_b(1)
       WEIGHT_b(8) = WEIGHT_b(2) 
       WEIGHT_b(9) = WEIGHT_b(3) 
       WEIGHT_b(10)= WEIGHT_b(4) 
       WEIGHT_b(11)= WEIGHT_b(5) 
       WEIGHT_b(12)= WEIGHT_b(6) 
      end select 
      !--
      do ibvec=1,NB 
       write(6,*) b_LAT(:,ibvec)       
      enddo 
      write(6,*) 
      SUM_REAL=0.0D0          
      Do ibvec=1,NB 
       SUM_REAL=SUM_REAL+WEIGHT_b(ibvec) 
      ENDDO           
      WEIGHT=SUM_REAL          
      !--
      write(6,*)'WEIGHT_b=',WEIGHT_b(:) 
      write(6,*)'WEIGHT  =',WEIGHT              
 
!======================
!=== CHECK SUM RULE ===
!======================

      write(6,*) 
      write(6,*)'SUM RULE CHECK OF B-VEC'       
      write(6,*) 
      Do IC=1,3   
       Do JC=1,3   
        SUM_REAL=0.0D0 
        Do ibvec=1,NB      
         SUM_REAL
     +  =SUM_REAL+WEIGHT_b(ibvec)*VEC_b(IC,ibvec)*VEC_b(JC,ibvec)
        ENDDO!ibvec            
        write(6,*) IC,JC,SUM_REAL          
       ENDDO!JC    
      ENDDO!IC    

!=========================
!=== INTERSTATE MATRIX ===
!=========================

      WRITE(6,*) 
      WRITE(6,*)'==================='
      WRITE(6,*)' INTERSTATE MATRIX '
      WRITE(6,*)'==================='
      WRITE(6,*) 
      allocate(OVERLAP(NTB,NTB,NTK,NB));OVERLAP=0.0D0   
      allocate(KPT(NTK,NB));KPT=0  
      nvx=nwx2+1;nvy=nwy2+1;nvz=nwz2+1 
!--
!$OMP PARALLEL PRIVATE(ik,C0_BRA,ig,ibvec,SHIFT_b,ik_ib,C0_KET,
!$OMP&  OVERLAP_TMP)
      allocate(C0_BRA(0:NTG,NTB));C0_BRA(0,:)=0.0D0       
      allocate(C0_KET(0:NTG,NTB));C0_KET(0,:)=0.0D0       
      allocate(SHIFT_b(3));SHIFT_b=0 
      allocate(OVERLAP_TMP(NTB,NTB));OVERLAP_TMP=0.0d0 
!$OMP DO 
      do ik=1,NTK       
!C0_BRA
       C0_BRA(0,:)=0.0D0       
       do ig=1,NG0(ik) 
        C0_BRA(ig,:)=C0(ig,:,ik)  
       enddo 
!SHIFT VECTOR
       do ibvec=1,NB               
        call make_shift_vector(NTK,NB,ik,ibvec,SK0(1,1),b_LAT(1,1),
     +  SHIFT_b(1),ik_ib)
        KPT(ik,ibvec)=ik_ib                
!C0_KET
        C0_KET(0,:)=0.0D0       
        do ig=1,NG0(ik_ib) 
         C0_KET(ig,:)=C0(ig,:,ik_ib)  
        enddo 
!CALC OVERLAP 
        CALL CALC_OVERLAP(NTK,NTB,NTG,nvx,nvy,nvz,ik,ik_ib,NG0(1),
     +  KG0(1,1,1),N_BAND(1),SHIFT_b(1),C0_BRA(0,1),C0_KET(0,1),
     +  OVERLAP_TMP(1,1),N_BAND_BTM(1))  
        OVERLAP(:,:,ik,ibvec)=OVERLAP_TMP(:,:)          
       !--
       enddo!ibvec                
      enddo!ik                           
!$OMP END DO 
      deallocate(C0_BRA,C0_KET,SHIFT_b,OVERLAP_TMP) 
!$OMP END PARALLEL 
      write(6,*)'finish calc ISM' 
      !--
      WRITE(6,*) 
      write(6,*)'================================================='
      write(6,*)' LOG{det[M(k,k+b)]} only when M is square matrix '
      write(6,*)'=================================================' 
      WRITE(6,*) 
      do ik=1,NTK       
       do ibvec=1,NB               
        ik_ib=KPT(ik,ibvec) 
        if(N_BAND(ik)==N_BAND(ik_ib))then 
         nm=N_BAND(ik) 
         allocate(mat_tmp(nm,nm));mat_tmp=0.0d0 
         mat_tmp(1:nm,1:nm)=OVERLAP(1:nm,1:nm,ik,ibvec) 
         call calcdet(nm,mat_tmp(1,1),DET) 
         write(6,'(a8,i6,a8,i6,a8,2f15.10)')
     +   'ik=',ik,'ik_ib=',ik_ib,'log=',log(DET) 
         deallocate(mat_tmp) 
        endif 
       enddo!ibvec                
      enddo!ik                           

!==============================================
!=== INITIAL GUESS FOR OMEGA_I MINIMIZATION ===
!==============================================

      allocate(orbtype(nigs)) 
      allocate(ALPHA_GAUSS(nigs));ALPHA_GAUSS(:)=0.0d0
      allocate(TAU_GAUSS(3,nigs));TAU_GAUSS(:,:)=0.0d0
      allocate(NORM_GAUSS(nigs));NORM_GAUSS(:)=0.0d0
      allocate(loc_x(3,nigs));loc_x=0.0d0
      allocate(loc_y(3,nigs));loc_y=0.0d0
      allocate(loc_z(3,nigs));loc_z=0.0d0
      !--
      do ig=1,nigs 
       orbtype(ig)=vec_ini(ig)%orb
       ALPHA_GAUSS(ig)=vec_ini(ig)%a 
       TAU_GAUSS(1,ig)=vec_ini(ig)%x
       TAU_GAUSS(2,ig)=vec_ini(ig)%y 
       TAU_GAUSS(3,ig)=vec_ini(ig)%z 
       loc_x(:,ig)=vec_ini(ig)%lx(:)
       loc_y(:,ig)=vec_ini(ig)%ly(:)
       loc_z(:,ig)=vec_ini(ig)%lz(:)
      enddo 
!--
!20170406 
      allocate(LGAUSS(nigs));LGAUSS(:)=0 
      allocate(MGAUSS(nigs));MGAUSS(:)=0 
      do ig=1,nigs 
       if(orbtype(ig)=='s')then 
        LGAUSS(ig)=1;MGAUSS(ig)=1   
       endif 
       if(orbtype(ig)=='px')then 
        LGAUSS(ig)=2;MGAUSS(ig)=1   
       endif 
       if(orbtype(ig)=='py')then 
        LGAUSS(ig)=2;MGAUSS(ig)=2   
       endif 
       if(orbtype(ig)=='pz')then 
        LGAUSS(ig)=2;MGAUSS(ig)=3   
       endif 
       if(orbtype(ig)=='dxy')then 
        LGAUSS(ig)=3;MGAUSS(ig)=1   
       endif 
       if(orbtype(ig)=='dyz')then 
        LGAUSS(ig)=3;MGAUSS(ig)=2   
       endif 
       if(orbtype(ig)=='dz2')then 
        LGAUSS(ig)=3;MGAUSS(ig)=3   
       endif 
       if(orbtype(ig)=='dzx')then 
        LGAUSS(ig)=3;MGAUSS(ig)=4   
       endif 
       if(orbtype(ig)=='dx2')then 
        LGAUSS(ig)=3;MGAUSS(ig)=5   
       endif 
      enddo!ig 
!--
!NORMALIZE
      Do ig=1,nigs
       iL=LGAUSS(ig) 
       IF(iL==1)THEN 
        NORM_GAUSS(ig)
     +  =(8.0D0*ALPHA_GAUSS(ig)/pi)**0.25D0
     +  *(4.0D0*ALPHA_GAUSS(ig)/1.0D0)**0.50D0    
       ENDIF 
       IF(iL==2)THEN 
        NORM_GAUSS(ig)
     +  =(8.0D0*ALPHA_GAUSS(ig)/pi)**0.25D0    
     +  *(4.0D0*ALPHA_GAUSS(ig)/1.0D0)**0.50D0    
     +  *(4.0D0*ALPHA_GAUSS(ig)/3.0D0)**0.50D0    
       ENDIF 
       IF(iL==3)THEN 
        NORM_GAUSS(ig)
     +  =(8.0D0*ALPHA_GAUSS(ig)/pi)**0.25D0    
     +  *(4.0D0*ALPHA_GAUSS(ig)/1.0D0)**0.50D0    
     +  *(4.0D0*ALPHA_GAUSS(ig)/3.0D0)**0.50D0    
     +  *(4.0D0*ALPHA_GAUSS(ig)/5.0D0)**0.50D0    
       ENDIF 
      ENDDO!ig 
      !--
      write(6,*) 
      write(6,*)'=================='
      write(6,*)'INITIAL GUESS DATA'
      write(6,*)'=================='
      write(6,*) 
      do ig=1,nigs 
       write(6,'(2i5,4f15.10)')
     + LGAUSS(ig),MGAUSS(ig),ALPHA_GAUSS(ig),TAU_GAUSS(:,ig)
      enddo 
      write(6,*) 
      write(6,*)'=================='
      write(6,*)'COORD IN CARTESIAN'
      write(6,*)'=================='
      write(6,*) 
      do ig=1,nigs
       tmp(:)=TAU_GAUSS(:,ig)  
       TAU_GAUSS(:,ig)=tmp(1)*a1(:)+tmp(2)*a2(:)+tmp(3)*a3(:) 
       write(6,'(3f20.12)') TAU_GAUSS(:,ig) 
      enddo 
      write(6,*) 
!--
!A_MAT
      write(6,*) 
      write(6,*)'=========='
      write(6,*)'CALC A_MAT' 
      write(6,*)'=========='
      write(6,*) 
      allocate(A_MAT(NTB,nigs,NTK));A_MAT=0.0d0          
!--
!$OMP PARALLEL PRIVATE(ik,A_TMP)
      allocate(A_TMP(NTB,nigs));A_TMP=0.0d0 
!$OMP DO 
      do ik=1,NTK         
       call make_amat(nigs,NTB,NTG,NG0(ik),SK0(1,ik),KG0(1,1,ik),
     + b1(1),b2(1),b3(1),TAU_GAUSS(1,1),ALPHA_GAUSS(1),
     + loc_x(1,1),loc_y(1,1),loc_z(1,1),LGAUSS(1),MGAUSS(1),
     + N_BAND(ik),C0(1,1,ik),N_BAND_BTM(ik),NORM_GAUSS(1),VOLUME,
     + A_TMP(1,1)) 
       A_MAT(:,:,ik)=A_TMP(:,:) 
      enddo!ik           
!$OMP END DO 
      deallocate(A_TMP) 
!$OMP END PARALLEL 
!--
!     do ik=1,NTK 
!      WRITE(6,*)'ik=',ik 
!      do i_band=1,N_BAND(ik)       
!       do j_gauss=1,nigs 
!        write(6,'(2i5,2f15.10)')i_band,j_gauss,A_MAT(i_band,j_gauss,ik)
!       enddo!j_gauss
!      enddo!i_band 
!      WRITE(6,*) 
!     enddo!ik  
!     WRITE(6,*) 
!--
!BMAT 
!SP2 WF 
!      B_MAT(:,:) = 0.0D0 
!      B_MAT(1,1) = DSQRT(1.0D0/3.0D0)  
!      B_MAT(1,2) = DSQRT(1.0D0/3.0D0)  
!      B_MAT(1,3) = DSQRT(1.0D0/3.0D0)  
!      B_MAT(2,1) =-DSQRT(2.0D0/3.0D0)  
!      B_MAT(2,2) = DSQRT(1.0D0/6.0D0)  
!      B_MAT(2,3) = DSQRT(1.0D0/6.0D0)  
!      B_MAT(3,4) = 1.0D0 
!      B_MAT(4,2) = DSQRT(1.0D0/2.0D0)  
!      B_MAT(4,3) =-DSQRT(1.0D0/2.0D0)  
!      B_MAT(5,5) = DSQRT(1.0D0/3.0D0)  
!      B_MAT(5,6) = DSQRT(1.0D0/3.0D0)  
!      B_MAT(5,7) = DSQRT(1.0D0/3.0D0)  
!      B_MAT(6,5) = DSQRT(2.0D0/3.0D0)  
!      B_MAT(6,6) =-DSQRT(1.0D0/6.0D0)  
!      B_MAT(6,7) =-DSQRT(1.0D0/6.0D0)  
!      B_MAT(7,8) = 1.0D0 
!      B_MAT(8,6) =-DSQRT(1.0D0/2.0D0)  
!      B_MAT(8,7) = DSQRT(1.0D0/2.0D0)  
!SP3 WF 
!      B_MAT(:,:) = 0.0D0 
!      B_MAT(1,1) =  0.5D0  
!      B_MAT(1,2) =  0.5D0  
!      B_MAT(1,3) =  0.5D0  
!      B_MAT(1,4) =  0.5D0  
!      B_MAT(2,1) = -0.5D0  
!      B_MAT(2,2) =  0.5D0  
!      B_MAT(2,3) = -0.5D0  
!      B_MAT(2,4) =  0.5D0  
!      B_MAT(3,1) =  0.5D0  
!      B_MAT(3,2) = -0.5D0  
!      B_MAT(3,3) = -0.5D0  
!      B_MAT(3,4) =  0.5D0  
!      B_MAT(4,1) = -0.5D0  
!      B_MAT(4,2) = -0.5D0  
!      B_MAT(4,3) =  0.5D0  
!      B_MAT(4,4) =  0.5D0  
!      B_MAT(5,5) =  0.5D0  
!      B_MAT(5,6) =  0.5D0  
!      B_MAT(5,7) =  0.5D0  
!      B_MAT(5,8) =  0.5D0  
!      B_MAT(6,5) =  0.5D0  
!      B_MAT(6,6) = -0.5D0  
!      B_MAT(6,7) =  0.5D0  
!      B_MAT(6,8) = -0.5D0  
!      B_MAT(7,5) = -0.5D0  
!      B_MAT(7,6) =  0.5D0  
!      B_MAT(7,7) =  0.5D0  
!      B_MAT(7,8) = -0.5D0  
!      B_MAT(8,5) =  0.5D0  
!      B_MAT(8,6) =  0.5D0  
!      B_MAT(8,7) = -0.5D0  
!      B_MAT(8,8) = -0.5D0  
!--
!A_MAT=A_MAT*B_MAT
      allocate(A_TMP(NTB,nigs));A_TMP=0.0d0 
      do ik=1,NTK 
       A_TMP(:,:)=A_MAT(:,:,ik) 
       do ib=1,N_BAND(ik)
        do ig=1,NGAUSS
         SUM_CMPX=0.0D0 
         do jg=1,nigs
          SUM_CMPX=SUM_CMPX+A_TMP(ib,jg)*B_MAT(jg,ig) 
         ENDDO 
         A_MAT(ib,ig,ik)=SUM_CMPX 
        ENDDO!ig 
       ENDDO!ib
      ENDDO!ik 
      write(6,*)'finish calc A_MAT' 
      !--
      allocate(S_MAT(NGAUSS,NGAUSS));S_MAT(:,:)=0.0D0
      allocate(S_TMP(NGAUSS,NGAUSS));S_TMP(:,:)=0.0d0
      allocate(S_EIG(NGAUSS));S_EIG(:)=0.0D0
      allocate(X_MAT(NGAUSS,NGAUSS));X_MAT(:,:)=0.0D0
      allocate(S_HALF_MAT(NGAUSS,NGAUSS));S_HALF_MAT(:,:)=0.0D0
      allocate(C_MAT(NTB,NGAUSS,NTK));C_MAT(:,:,:)=0.0D0          
      !--
      do ik=1,NTK        
!S_MAT
       S_MAT(:,:)=0.0d0 
       do ig=1,NGAUSS
        do jg=1,NGAUSS
         SUM_CMPX=0.0D0         
         do ib=1,N_BAND(ik)        
          SUM_CMPX=SUM_CMPX+CONJG(A_MAT(ib,ig,ik))*A_MAT(ib,jg,ik)   
         enddo!ib
         S_MAT(ig,jg)=SUM_CMPX 
        enddo!jg
       enddo!ig       
!--
!      write(6,*)'IK=',ik     
!       do ig=1,NGAUSS
!        do jg=1,NGAUSS
!         write(6,'(2i5,2f15.10)') ig,jg,S_MAT(ig,jg) 
!        enddo!ig       
!       enddo!jg       
!       write(6,*) 
!diag S_MAT
       S_EIG(:)=0.0D0                
       S_TMP(:,:)=S_MAT(:,:) 
       call diagV(NGAUSS,S_TMP(1,1),S_EIG(1)) 
       X_MAT(:,:)=S_TMP(:,:)
       do ig=1,NGAUSS 
        if(S_EIG(ig)<=0.0D0)then 
         write(6,*)'WARNING diag S_MAT' 
         write(6,*) ig,S_EIG(ig) 
         S_EIG(ig)=1.0D-10 
        endif 
       enddo 
!S_HALF_MAT
       S_HALF_MAT(:,:)=0.0D0 
       S_EIG(:)=1.0D0/DSQRT(S_EIG(:)) 
       do ig=1,NGAUSS
        do jg=1,NGAUSS
         SUM_CMPX=0.0D0         
         do kg=1,NGAUSS      
          SUM_CMPX=SUM_CMPX+X_MAT(ig,kg)*S_EIG(kg)*CONJG(X_MAT(jg,kg)) 
         enddo!kg 
         S_HALF_MAT(ig,jg)=SUM_CMPX  
        enddo!jg 
       enddo!ig 
!C_MAT
       do ib=1,N_BAND(ik)         
        do jg=1,NGAUSS
         SUM_CMPX=0.0D0             
         do kg=1,NGAUSS             
          SUM_CMPX=SUM_CMPX+A_MAT(ib,kg,ik)*S_HALF_MAT(kg,jg)        
         enddo 
         C_MAT(ib,jg,ik)=SUM_CMPX         
        enddo!jg       
       enddo!ib    
      enddo!ik ===CONSTRUCT C_MAT ROOP===      
      deallocate(S_MAT,S_TMP,S_EIG,X_MAT,S_HALF_MAT)
      ! 
      !write(6,*)'C_MAT'
      !write(6,*) C_MAT 
      ! 
!--
!      WRITE(6,*) 
!      WRITE(6,*)'==========================='
!      WRITE(6,*)'CHECK OF UNITALITY OF C_MAT'
!      WRITE(6,*)'==========================='
!      WRITE(6,*) 
!      do ik=1,NTK      
!       do i_gauss=1,NGAUSS       
!        do j_gauss=1,NGAUSS
!         SUM_CMPX=0.0D0             
!         do i_band=1,N_BAND(ik)        
!          SUM_CMPX 
!     +   =SUM_CMPX 
!     +   +CONJG(C_MAT(i_band,i_gauss,ik))  
!     +   *      C_MAT(i_band,j_gauss,ik)        
!         enddo 
!         write(6,'(2I5,2F20.15)') i_gauss,j_gauss,SUM_CMPX
!        enddo!j_gauss     
!       enddo!i_gauss     
!       write(6,*) 
!      enddo!ik            
!
!=========================================================
!=== CALCULATE THE NUMBER OF BANDS INSIDE INNER WINDOW ===
!=========================================================

      if(set_inner_window)then 
       allocate(N_BAND_inner(NTK));N_BAND_inner=0
       allocate(N_BAND_BTM_inner(NTK));N_BAND_BTM_inner=0
       allocate(inner(NTB,NTK));inner=0
       write(6,*) 
       write(6,*)'===================================='
       write(6,*)'BAND INFO INSIDE INNER ENERGY WINDOW'
       write(6,*)'===================================='
       write(6,*) 
       do ik=1,NTK 
        SUM_INT=0 
        do ib=1,N_BAND(ik) 
         if(E_LOWER_inner<=E_EIG(N_BAND_BTM(ik)+ib,ik).and.
     +      E_EIG(N_BAND_BTM(ik)+ib,ik)<=E_UPPER_inner)then 
          inner(ib,ik)=1
          SUM_INT=SUM_INT+1 
          if(SUM_INT==1) N_BAND_BTM_inner(ik)=ib-1  
         endif 
        enddo 
        N_BAND_inner(ik)=SUM_INT
!--
!check N_BAND_inner(ik) <= N_wannier 
        if(N_BAND_inner(ik)>N_wannier)then 
         write(6,'(a)')'inner energy window setting: WRONG'
         write(6,'(a)')'NBAND_inner SHOULD BE SMALLER THAN N_wannier'
         write(6,'(a,3i8)')'ik, N_BAND_inner(ik), N_wannier',
     +                      ik, N_BAND_inner(ik), N_wannier 
         stop 
        endif  
        write(6,*) ik,N_BAND_BTM_inner(ik),N_BAND_inner(ik) 
       enddo!ik  
       write(6,*) 
!--
!initial CMAT
       do ik=1,NTK 
!X=QPQ 
        allocate(P_MAT(NTB,NTB));P_MAT=0.0d0
        do ib=1,N_BAND(ik)
         do jb=1,N_BAND(ik)
          SUM_CMPX=0.0d0 
          do jg=1,NGAUSS 
           SUM_CMPX=SUM_CMPX+C_MAT(ib,jg,ik)*CONJG(C_MAT(jb,jg,ik)) 
          enddo 
          P_MAT(ib,jb)=SUM_CMPX 
         enddo!jb  
        enddo!ib
        allocate(Qin(NTB,NTB));Qin=0.0d0
        do ib=1,N_BAND(ik)
         if(inner(ib,ik)==1)cycle 
         Qin(ib,ib)=1.0d0 
        enddo 
!--
!        do ib=1,N_BAND_inner(ik)
!         do jb=1,N_BAND_inner(ik) 
!          Qin(ib+N_BAND_BTM_inner(ik),jb+N_BAND_BTM_inner(ik)) 
!     +   =Qin(ib+N_BAND_BTM_inner(ik),jb+N_BAND_BTM_inner(ik)) 
!     +   -P_MAT(ib+N_BAND_BTM_inner(ik),jb+N_BAND_BTM_inner(ik)) 
!         enddo 
!        enddo 
!--
        allocate(X_MAT(NTB,NTB));X_MAT=0.0d0  
        do ib=1,N_BAND(ik)
         do jb=1,N_BAND(ik) 
          SUM_CMPX=0.0d0 
          do k=1,N_BAND(ik) 
           do l=1,N_BAND(ik) 
            SUM_CMPX=SUM_CMPX+CONJG(Qin(k,ib))*P_MAT(k,l)*Qin(l,jb) 
           enddo!l 
          enddo!k 
          X_MAT(ib,jb)=SUM_CMPX 
         enddo!jb
        enddo!ib  
!XC=Cx
        nm=N_BAND(ik) 
        allocate(mat_tmp(nm,nm));mat_tmp(:,:)=0.0d0 
        allocate(eig_tmp(nm));eig_tmp(:)=0.0d0        
        mat_tmp(1:nm,1:nm)=X_MAT(1:nm,1:nm) 
        call diagV(nm,mat_tmp(1,1),eig_tmp(1)) 
!descending order
        X_MAT=0.0d0 
        do ib=1,nm 
         X_MAT(1:nm,ib)=mat_tmp(1:nm,nm-ib+1) 
        enddo 
        deallocate(mat_tmp,eig_tmp) 
!CMAT 
        do jg=1,NGAUSS-N_BAND_inner(ik) 
         C_MAT(:,jg,ik)=X_MAT(:,jg) 
        enddo               
        do jg=1,N_BAND_inner(ik) 
         C_MAT(:,NGAUSS-N_BAND_inner(ik)+jg,ik)=0.0d0 
        enddo 
        do ib=1,N_BAND_inner(ik)
         i=N_BAND_BTM_inner(ik)+ib 
         j=NGAUSS-N_BAND_inner(ik)+ib 
         C_MAT(i,j,ik)=1.0d0 
        enddo 
        deallocate(P_MAT,Qin,X_MAT) 
       enddo!ik   
      else
       allocate(N_BAND_inner(NTK));N_BAND_inner=0
       allocate(N_BAND_BTM_inner(NTK));N_BAND_BTM_inner=0
       allocate(inner(NTB,NTK));inner=0
      endif!set_inner_window 
!--
!      write(6,*)'INNER'
!      do ik=1,NTK 
!       write(6,*) ik,N_BAND_BTM_inner(ik),N_BAND_inner(ik) 
!      enddo 
!--
      WRITE(6,*) 
      WRITE(6,*)'==========================================='        
      WRITE(6,*)'=== MINIMIZATION OF SPILLAGE FUNCTIONAL ==='
      WRITE(6,*)'==========================================='        
      WRITE(6,*) 
      !--
      allocate(Z_EIG(NTB));Z_EIG(:)=0.0D0            
      allocate(Z_MAT(NTB,NTB,NTK));Z_MAT(:,:,:)=0.0D0        
      allocate(Z_TMP(NTB,NTB));Z_TMP(:,:)=0.0D0               
      allocate(C_TMP(NTB,NTB));C_TMP(:,:)=0.0D0      
      allocate(M_MAT(NTB,NTB,NTK,NB));M_MAT(:,:,:,:)=0.0D0            
!--
!N    : NGAUSS  
!N_k  : N_BAND(ik) 
!N_k+b: N_BAND(ik_ib) 
!M_k  : N_BAND_inner(ik) 
!--
      I_SCF=0           
      DEL_OMEGA_I=1.0D0 
      OMEGA_I_OLD=0.0D0        
      do while(DEL_OMEGA_I>EPS_SPILLAGE) 
!make M_MAT: size = N_k * N_k+b 
       M_MAT(:,:,:,:)=0.0D0            
       do ik=1,NTK         
        do ibvec=1,NB 
         ik_ib=KPT(ik,ibvec)       
         do ib=1,N_BAND(ik)       
          do jg=1,NGAUSS 
           SUM_CMPX=0.0D0             
           do kb=1,N_BAND(ik_ib)         
            SUM_CMPX=SUM_CMPX+OVERLAP(ib,kb,ik,ibvec)*C_MAT(kb,jg,ik_ib)
           enddo 
           M_MAT(ib,jg,ik,ibvec)=SUM_CMPX       
          enddo!jg
         enddo!ib     
        enddo!ibvec  
       enddo!ik         
       !--
       SUM_REAL=0.0D0           
       OMEGA_I_NEW=0.0D0           
       C_MAT(:,:,:)=0.0D0 
       do ik=1,NTK           
!make Z_TMP: size = N_k * N_k 
        Z_TMP(:,:)=0.0D0               
        do ib=1,N_BAND(ik)    
         do jb=1,N_BAND(ik)    
          SUM_CMPX=0.0D0                
          do ibvec=1,NB          
           do jg=1,NGAUSS 
            SUM_CMPX
     +     =SUM_CMPX
     +     +WEIGHT_b(ibvec)
     +     *M_MAT(ib,jg,ik,ibvec)*CONJG(M_MAT(jb,jg,ik,ibvec))  
           enddo 
          enddo 
          Z_TMP(ib,jb)=SUM_CMPX                  
         enddo!jb         
        enddo!ib         
!make Z_MAT: damped Z_TMP
        IF(I_SCF==0) THEN  
         Z_MAT(:,:,ik)=Z_TMP(:,:)          
        ELSE 
         Z_MAT(:,:,ik)=DAMP*Z_TMP(:,:)+(1.0D0-DAMP)*Z_MAT(:,:,ik)   
        ENDIF        
!make X: size = (N_k-M_k)*(N_k-M_k) 
        nm=N_BAND(ik)-N_BAND_inner(ik)!inner window    
        allocate(mat_tmp(nm,nm));mat_tmp(:,:)=0.0d0 
        allocate(eig_tmp(nm));eig_tmp(:)=0.0d0        
        i=0 
        do ib=1,N_BAND(ik) 
         if(inner(ib,ik)==1)cycle 
         i=i+1 
         j=0 
         do jb=1,N_BAND(ik) 
          if(inner(jb,ik)==1)cycle
          j=j+1 
          mat_tmp(i,j)=Z_MAT(ib,jb,ik)
         enddo 
        enddo   
!diag X: XC'=C'x 
        eig_tmp(:)=0.0d0 
        !--
        !20180412 by Y.Nomura 
        !--
        if(nm>0)then
         call diagV(nm,mat_tmp(1,1),eig_tmp(1)) 
        endif 
        !--
!descending order
        Z_EIG(:)=0.0D0            
        C_TMP(:,:)=0.0d0 
        do ib=1,nm
         Z_EIG(ib)=eig_tmp(nm-ib+1)
         C_TMP(1:nm,ib)=mat_tmp(1:nm,nm-ib+1) 
        enddo 
        deallocate(mat_tmp,eig_tmp) 
!make C_MAT: size = N_k * N 
        i=0 
        do ib=1,N_BAND(ik) 
         if(inner(ib,ik)==0)then 
          i=i+1
          !C_MAT(ib,:,ik)=C_TMP(i,:) ###bag### 20171208
          C_MAT(ib,1:NGAUSS-N_BAND_inner(ik),ik)
     +   =C_TMP(i ,1:NGAUSS-N_BAND_inner(ik)) 
         else!(inner(ib,ik)==1) 
          C_MAT(ib,:,ik)=0.0d0 
         endif 
        enddo!ib  
        do jg=1,N_BAND_inner(ik) 
         C_MAT(:,NGAUSS-N_BAND_inner(ik)+jg,ik)=0.0d0 
        enddo 
        do ib=1,N_BAND_inner(ik)
         i=N_BAND_BTM_inner(ik)+ib 
         j=NGAUSS-N_BAND_inner(ik)+ib 
         C_MAT(i,j,ik)=1.0d0 
        enddo 
!--
!SPILLAGE 
        do ig=1,NGAUSS
         do jb=1,N_BAND(ik)
          do ib=1,N_BAND(ik)
           SUM_REAL
     +    =SUM_REAL 
     +    +dble(conjg(C_MAT(ib,ig,ik))*Z_MAT(ib,jb,ik)*C_MAT(jb,ig,ik)) 
          enddo 
         enddo 
        enddo 
        OMEGA_I_NEW=OMEGA_I_NEW+dble(NGAUSS)*WEIGHT 
!--
!        nm=NGAUSS-N_BAND_inner(ik)!inner window
!        do jg=1,nm
!         SUM_REAL=SUM_REAL+Z_EIG(jg) 
!        enddo 
!        OMEGA_I_NEW=OMEGA_I_NEW+dble(nm)*WEIGHT 
!--
       enddo!ik 
       OMEGA_I_NEW=(OMEGA_I_NEW-SUM_REAL)/DBLE(NTK)      
       DEL_OMEGA_I=DABS(OMEGA_I_NEW-OMEGA_I_OLD)            
       OMEGA_I_OLD=OMEGA_I_NEW          
       I_SCF=I_SCF+1                                       
       write(6,'(A8,I5,A15,F20.10,A10,F20.10)')'I_SCF=',I_SCF, 
     +         'DEL_OMEGA_I=',DEL_OMEGA_I,'OMEGA_I',OMEGA_I_NEW        
      enddo!do-while##### OMEGA_I MINIMIZATION ##### 
!--
!      do ik=1,NTK 
!       if(N_BAND_inner(ik)/=0)then 
!        write(6,*)'CMAT(fns) ik=',ik  
!        do ib=1,N_BAND(ik) 
!         write(6,'(100f10.5)')(C_MAT(ib,jg,ik),jg=1,NGAUSS)  
!        enddo 
!        write(6,*) 
!       endif 
!      enddo!ik
!===============================================
!=== PSEUDO EIGENSTATE AND PSEUDO EIGENVALUE === 
!===============================================

      allocate(PSEUDO_EIG(NGAUSS,NTK));PSEUDO_EIG(:,:)=0.0D0    
      allocate(PSEUDO_MAT(NGAUSS,NGAUSS,NTK));PSEUDO_MAT(:,:,:)=0.0D0 
      allocate(P_EIG(NGAUSS));P_EIG(:)=0.0d0 
      allocate(P_MAT_IN(NGAUSS,NGAUSS));P_MAT_IN(:,:)=0.0D0 
      allocate(P_MAT_OUT(NGAUSS,NGAUSS));P_MAT_OUT(:,:)=0.0D0 
      !--
      do ik=1,NTK               
!PSEUDO-KS 
       P_MAT_IN(:,:)=0.0D0          
       do ig=1,NGAUSS            
        do jg=1,NGAUSS           
         SUM_CMPX=0.0D0           
         do kb=1,N_BAND(ik)                 
          SUM_CMPX 
     +   =SUM_CMPX 
     +   +CONJG(C_MAT(kb,ig,ik))          
     +   *E_EIG(kb+N_BAND_BTM(ik),ik)      
     +   *C_MAT(kb,jg,ik)                     
         enddo 
         P_MAT_IN(ig,jg)=SUM_CMPX              
        enddo!jg        
       enddo!ig        
!diag PSEUDO-KS 
       P_EIG(:)=0.0D0            
       call diagV(NGAUSS,P_MAT_IN(1,1),P_EIG(1)) 
       P_MAT_OUT(:,:)=P_MAT_IN(:,:) 
       PSEUDO_MAT(:,:,ik)=P_MAT_OUT(:,:)          
       PSEUDO_EIG(:,ik)=P_EIG(:)          
      enddo!ik                    
      !deallocate(E_EIG) 
!--
!     do ik=1,NTK                  
!      write(6,'(30F15.10)')PSEUDO_EIG(:,ik)*au 
!     enddo 
!
!====================================================
!=== PREPARE INITIAL GUESS FOR OMEGA MINIMIZATION ===
!====================================================

      allocate(S_MAT(NGAUSS,NGAUSS));S_MAT(:,:)=0.0D0
      allocate(S_TMP(NGAUSS,NGAUSS));S_TMP(:,:)=0.0d0
      allocate(S_EIG(NGAUSS));S_EIG(:)=0.0D0
      allocate(X_MAT(NGAUSS,NGAUSS));X_MAT(:,:)=0.0D0
      allocate(S_HALF_MAT(NGAUSS,NGAUSS));S_HALF_MAT(:,:)=0.0D0
      Z_MAT(:,:,:)=0.0D0!INITIAL GUESS <---> CELL-PERIODIC 
      do ik=1,NTK               
!A_TMP
       A_TMP(:,:)=0.0D0             
       do ib=1,n_occ        
        do jg=1,NGAUSS          
         SUM_CMPX=0.0D0           
         do kb=1,N_BAND(ik)          
          do lg=1,NGAUSS    
           SUM_CMPX
     +    =SUM_CMPX
     +    +CONJG(PSEUDO_MAT(lg,ib,ik))          
     +    *CONJG(C_MAT(kb,lg,ik))          
     +    *A_MAT(kb,jg,ik)               
          ENDDO 
         ENDDO 
         A_TMP(ib,jg)=SUM_CMPX               
        ENDDO!jg      
       ENDDO!ib            
!S_MAT
       S_MAT(:,:)=0.0D0                   
       Do ig=1,NGAUSS
        Do jg=1,NGAUSS
         SUM_CMPX=0.0D0         
         Do kb=1,n_occ                      
          SUM_CMPX 
     +   =SUM_CMPX 
     +   +CONJG(A_TMP(kb,ig))
     +   *A_TMP(kb,jg) 
         ENDDO!kb
          S_MAT(ig,jg)=SUM_CMPX 
        ENDDO!jg
       ENDDO!ig       
!diag S_MAT
       S_EIG(:)=0.0D0            
       call diagV(NGAUSS,S_MAT(1,1),S_EIG(1)) 
       X_MAT(:,:)=S_MAT(:,:) 
       Do ig=1,NGAUSS 
        IF(S_EIG(ig)<=0.0D0) THEN 
         write(6,*) ig,S_EIG(ig) 
         S_EIG(ig)=1.0D-10 
        ENDIF       
       ENDDO 
!S_HALF_MAT
       S_EIG(:)=1.0D0/DSQRT(S_EIG(:)) 
       S_HALF_MAT(:,:)=0.0D0         
       Do ig=1,NGAUSS
        Do jg=1,NGAUSS
         SUM_CMPX=0.0D0         
         Do kg=1,NGAUSS      
          SUM_CMPX 
     +   =SUM_CMPX 
     +   +X_MAT(ig,kg) 
     +   *S_EIG(kg) 
     +   *CONJG(X_MAT(jg,kg))          
         ENDDO 
         S_HALF_MAT(ig,jg)=SUM_CMPX  
        ENDDO!jg 
       ENDDO!ig 
!A_MAT
       A_MAT(:,:,ik)=0.0D0         
       Do ib=1,n_occ              
        Do jg=1,NGAUSS
         SUM_CMPX=0.0D0             
         Do kg=1,NGAUSS             
          SUM_CMPX 
     +   =SUM_CMPX 
     +   +A_TMP(ib,kg) 
     +   *S_HALF_MAT(kg,jg)        
         ENDDO              
         A_MAT(ib,jg,ik)=SUM_CMPX         
        ENDDO!jg       
       ENDDO!ib    
!Z_MAT 
       Do ib=1,N_BAND(ik)          
        Do jg=1,NGAUSS 
         SUM_CMPX=0.0D0                   
         Do kg=1,NGAUSS 
          Do llb=1,n_occ                
           SUM_CMPX 
     +    =SUM_CMPX 
     +    +C_MAT(ib,kg,ik) 
     +    *PSEUDO_MAT(kg,llb,ik) 
     +    *A_MAT(llb,jg,ik)
          ENDDO 
         ENDDO 
         Z_MAT(ib,jg,ik)=SUM_CMPX          
        ENDDO!jg         
       ENDDO!ib          
      ENDDO!ik ##### INITIAL GUESS FOR OMEGA MINIMIZATION #####
      deallocate(S_MAT,S_TMP,S_EIG,X_MAT,S_HALF_MAT)
!--
!      WRITE(6,*) 
!      WRITE(6,*)'==========================='
!      WRITE(6,*)'CHECK OF UNITALITY OF Z_MAT'
!      WRITE(6,*)'==========================='
!      WRITE(6,*) 
!      do ik=1,NTK      
!       do i_gauss=1,NGAUSS       
!        do j_gauss=1,NGAUSS
!         SUM_CMPX=0.0D0             
!         do i_band=1,N_BAND(ik)        
!          SUM_CMPX 
!     +   =SUM_CMPX 
!     +   +CONJG(Z_MAT(i_band,i_gauss,ik))  
!     +   *Z_MAT(i_band,j_gauss,ik)        
!         enddo               
!         write(6,'(2i5,2f20.15)') i_gauss,j_gauss,SUM_CMPX
!        enddo!j_gauss     
!       enddo!i_gauss     
!       write(6,*)  
!      enddo!ik            
!
!=========================================
!=== MINIMIZATION OF SPREAD FUNCTIONAL ===
!=========================================

      WRITE(6,*) 
      WRITE(6,*)'====================================='
      WRITE(6,*)'INTERSTATE MATRIX; LOG{det[M(k,k+b)]}'
      WRITE(6,*)'====================================='
      WRITE(6,*) 

      M_MAT(:,:,:,:)=OVERLAP(:,:,:,:)                      
      OVERLAP(:,:,:,:)=0.0D0                                          
      Do ik=1,NTK 
       Do ibvec=1,NB      
        ik_ib=KPT(ik,ibvec)        
        Do i_band=1,n_occ!NOTICE n_occ=NGAUSS  
         Do j_band=1,n_occ!NOTICE n_occ=NGAUSS  
          SUM_CMPX=0.0D0           
          Do k_band=1,N_BAND(ik)      
           Do l_band=1,N_BAND(ik_ib)      
            SUM_CMPX 
     +     =SUM_CMPX 
     +     +CONJG(Z_MAT(k_band,i_band,ik)) 
     +     *      M_MAT(k_band,l_band,ik,ibvec)            
     +     *      Z_MAT(l_band,j_band,ik_ib) 
           ENDDO 
          ENDDO 
          OVERLAP(i_band,j_band,ik,ibvec)=SUM_CMPX          
         ENDDO!j_band        
        ENDDO!i_band        
        !--
        !DET
         nm=n_occ 
         allocate(mat_tmp(nm,nm));mat_tmp(:,:)=0.0d0 
         mat_tmp(1:nm,1:nm)=OVERLAP(1:nm,1:nm,ik,ibvec) 
         call calcdet(nm,mat_tmp(1,1),DET) 
         deallocate(mat_tmp) 
        !--
        write(6,'(a8,i5,a8,i5,a8,2f15.10)')
     +  'IK=',ik,'IK_IB=',ik_ib,'LOG=',LOG(DET) 
       ENDDO!ibvec                
       WRITE(6,*) 
      ENDDO!ik                           

!=================================
!=== INITIAL SPREAD FUNCTIONAL ===
!=================================

      allocate(WF_CHARGE(n_occ,NTK,NB));WF_CHARGE(:,:,:)=0.0d0       
      allocate(TR_GRAD(NTK));TR_GRAD(:)=0.0d0 
      allocate(TR_GRAD_OLD(NTK));TR_GRAD_OLD(:)=0.0d0   
      allocate(GRADIENT(n_occ,n_occ,NTK));GRADIENT(:,:,:)=0.0d0     
      allocate(DIRECTION(n_occ,n_occ,NTK));DIRECTION(:,:,:)=0.0d0      
      allocate(DIRECTION2(n_occ,n_occ,NTK));DIRECTION2(:,:,:)=0.0d0 
      allocate(UNITARY(n_occ,n_occ,NTK));UNITARY(:,:,:)=0.0d0        
      allocate(U_ORG(n_occ,n_occ,NTK));U_ORG(:,:,:)=0.0d0               
      allocate(U_OLD(n_occ,n_occ,NTK));U_OLD(:,:,:)=0.0d0               
      !--
      ILS=1           
      CALL CALC_OMEGA(NTK,NTB,n_occ,NB,VEC_b,ILS,WEIGHT_b,OVERLAP,
     +     WF_CHARGE,OMEGA_I,OMEGA_OD,OMEGA_D)                   
      !--
      write(6,*)'OMEGA_I =',OMEGA_I 
      write(6,*)'OMEGA_OD=',OMEGA_OD      
      write(6,*)'OMEGA_D =',OMEGA_D 
      write(6,*)'OMEGA_total=',OMEGA_I+OMEGA_D+OMEGA_OD
      WRITE(6,*) 
      WRITE(6,*)'========================================='        
      WRITE(6,*)'=== MINIMIZATION OF SPREAD FUNCTIONAL ==='
      WRITE(6,*)'========================================='        
      WRITE(6,*) 
      !--
      M_MAT(:,:,:,:)=OVERLAP(:,:,:,:)                      
      DEL_SPREAD=1.0d0!EPS_SPREAD=1.0D-6 default
      SPREAD_OLD=0.0D0 
      I_STEP=1        
      SPREAD=OMEGA_I+OMEGA_OD+OMEGA_D               
      do while(DEL_SPREAD>EPS_SPREAD) 
!--
!initial step 
      ILS=0 
      DELTA(1)=0.0D0  
      OMEGA(1)=SPREAD       
      write(6,*)'ILS=',ILS,DELTA(1),OMEGA(1)
!--
!CALC GRADIENT
      CALL CALC_GRADIENT(NTK,NTB,n_occ,NB,WF_CHARGE,WEIGHT_b,OVERLAP,
     +                   tci,GRADIENT,TR_GRAD) 
      DIRECTION(:,:,:)=GRADIENT(:,:,:) 
!--
!CONJUGATE GRADIENT
!      IF(I_STEP==1) THEN 
!       DIRECTION(:,:,:)=GRADIENT(:,:,:) 
!      ELSE                 
!       Do ik=1,NTK 
!         DIRECTION(:,:,ik) 
!     + = GRADIENT(:,:,ik) 
!     + + (TR_GRAD(ik)/TR_GRAD_OLD(ik)) 
!     + * DIRECTION(:,:,ik) 
!       ENDDO          
!      ENDIF              
!      TR_GRAD_OLD(:)=TR_GRAD(:) 
!--
!DIRECTION^2
      DIRECTION2(:,:,:)=0.0D0        
      Do ik=1,NTK        
       Do i_band=1,n_occ       
        Do j_band=1,n_occ       
         SUM_CMPX=0.0D0 
         Do k_band=1,n_occ    
          SUM_CMPX 
     +   =SUM_CMPX 
     +   +DIRECTION(i_band,k_band,ik)     
     +   *DIRECTION(k_band,j_band,ik) 
         ENDDO!k_band     
         DIRECTION2(i_band,j_band,ik)=SUM_CMPX                
        ENDDO!j_band           
       ENDDO!i_band           
      ENDDO!ik       
!===
!=== ONE-DIMENSIONAL SEARCH 
!===
      ILS=1           
      STEP_LENGTH=MAX_STEP_LENGTH/4.0D0/WEIGHT 
      DELTA(2)=STEP_LENGTH 
      CALL UNITARY_GEN(NTK,n_occ,ILS,DIRECTION(1,1,1),DIRECTION2(1,1,1),
     +     STEP_LENGTH,UNITARY(1,1,1),U_ORG(1,1,1),U_OLD(1,1,1))   
      OVERLAP(:,:,:,:)=M_MAT(:,:,:,:) 
      CALL UPDATE_OVERLAP(NTK,NTB,n_occ,NB,KPT,UNITARY,OVERLAP)
      CALL CALC_OMEGA(NTK,NTB,n_occ,NB,VEC_b,ILS,      
     +     WEIGHT_b,OVERLAP,WF_CHARGE,OMEGA_I,OMEGA_OD,OMEGA_D) 
      OMEGA(2)=OMEGA_I+OMEGA_D+OMEGA_OD
      write(6,*)'ILS=',ILS,DELTA(2),OMEGA(2)
!20170917
      if(abs(OMEGA(1)-OMEGA(2))<1.0d-10)then
!update M_MAT
       M_MAT(:,:,:,:)=OVERLAP(:,:,:,:)                      
!update A_MAT
       CALL A_MAT_UPDATE(NTK,NTB,NGAUSS,n_occ,UNITARY,A_MAT)      
!CONVERGENCE CHECK
       SPREAD_OLD=OMEGA(1) 
       SPREAD=OMEGA(2) 
       DEL_SPREAD=ABS(SPREAD-SPREAD_OLD) 
       SPREAD_OLD=OMEGA(2)               
       write(6,*) 
       write(6,'(a47)')'==============================================' 
       write(6,'(a47)')' SPREAD HAS ALREADY CONVERGED IN INITIAL STEP ' 
       write(6,'(a47)')'=============================================='
       write(6,*) 
       write(6,'(a26,i8,2f20.10)')'I_STEP SPREAD DEL_SPREAD:', 
     +                             I_STEP,SPREAD,DEL_SPREAD 
       write(6,*)  
       I_STEP=I_STEP+1          
      else
       do while(OMEGA(1)>OMEGA(2)) 
        ILS=ILS+1           
        DELTA(3) = DELTA(2)+STEP_LENGTH 
        CALL UNITARY_GEN(NTK,n_occ,ILS,DIRECTION(1,1,1),
     +  DIRECTION2(1,1,1),STEP_LENGTH,UNITARY(1,1,1),
     +  U_ORG(1,1,1),U_OLD(1,1,1))          
        OVERLAP(:,:,:,:)=M_MAT(:,:,:,:) 
        CALL UPDATE_OVERLAP(NTK,NTB,n_occ,NB,
     +                      KPT,UNITARY,OVERLAP)                       
        CALL CALC_OMEGA(NTK,NTB,n_occ,NB,VEC_b,ILS,          
     +                  WEIGHT_b,OVERLAP,WF_CHARGE,
     +                  OMEGA_I,OMEGA_OD,OMEGA_D)                   
        OMEGA(3)=OMEGA_I+OMEGA_D+OMEGA_OD
        write(6,*) 'ILS=',ILS,DELTA(3),OMEGA(3) 
        IF(OMEGA(2)<OMEGA(3)) EXIT 
        DELTA(1)=DELTA(2) 
        OMEGA(1)=OMEGA(2) 
        DELTA(2)=DELTA(3) 
        OMEGA(2)=OMEGA(3) 
       enddo!do while 
       !--
       DELTA(4)=0.5D0*(OMEGA(1)*(DELTA(2)**2-DELTA(3)**2) 
     +                +OMEGA(2)*(DELTA(3)**2-DELTA(1)**2)
     +                +OMEGA(3)*(DELTA(1)**2-DELTA(2)**2))
     +               /(OMEGA(1)*(DELTA(2)-DELTA(3))        
     +                +OMEGA(2)*(DELTA(3)-DELTA(1))        
     +                +OMEGA(3)*(DELTA(1)-DELTA(2)))       
       ILS=1           
       STEP_LENGTH=DELTA(4) 
       CALL UNITARY_GEN(NTK,n_occ,ILS,DIRECTION(1,1,1),
     +      DIRECTION2(1,1,1),STEP_LENGTH,UNITARY(1,1,1),
     +      U_ORG(1,1,1),U_OLD(1,1,1))      
       OVERLAP(:,:,:,:) = M_MAT(:,:,:,:) 
       CALL UPDATE_OVERLAP(NTK,NTB,n_occ,NB,KPT,UNITARY,OVERLAP) 
       CALL CALC_OMEGA(NTK,NTB,n_occ,NB,VEC_b,ILS,         
     +      WEIGHT_b,OVERLAP,WF_CHARGE,OMEGA_I,OMEGA_OD,OMEGA_D)
       OMEGA(4) = OMEGA_I+OMEGA_D+OMEGA_OD
       write(6,*)'DELTA_OPT=',DELTA(4),'OMEGA_OPT=',OMEGA(4) 
       write(6,*)'OMEGA_I =',OMEGA_I 
       write(6,*)'OMEGA_OD=',OMEGA_OD
       write(6,*)'OMEGA_D =',OMEGA_D
!update M_MAT
       M_MAT(:,:,:,:)=OVERLAP(:,:,:,:)                      
!update A_MAT
       CALL A_MAT_UPDATE(NTK,NTB,NGAUSS,n_occ,UNITARY,A_MAT)      
!CONVERGENCE CHECK
       SPREAD=OMEGA(4) 
       DEL_SPREAD=ABS(SPREAD-SPREAD_OLD) 
       SPREAD_OLD=OMEGA(4)               
       write(6,'(a26,i8,2f20.10)')'I_STEP SPREAD DEL_SPREAD:', 
     +                             I_STEP,SPREAD,DEL_SPREAD 
       write(6,*)  
       I_STEP=I_STEP+1
      endif
      !--           
      ENDDO!##### SPREAD-MINIMIZATION ROOP #####         
!--
!      WRITE(6,*) 
!      WRITE(6,*)'========================'
!      WRITE(6,*)'UNITALITY CHECK OF A_MAT'
!      WRITE(6,*)'========================'
!      WRITE(6,*)  
!      do ik=1,NTK 
!       do ib=1,n_occ   
!        do jb=1,n_occ   
!         SUM_CMPX=0.0D0            
!         do kb=1,n_occ       
!          SUM_CMPX=SUM_CMPX+CONJG(A_MAT(kb,ib,ik))*A_MAT(kb,jb,ik) 
!         enddo 
!         write(6,*) ib,jb,SUM_CMPX        
!        enddo  
!       enddo  
!       write(6,*)       
!      enddo!ik        
!--
!wannier center 
      allocate(wannier_center(3,n_occ));wannier_center(:,:)=0.0d0  
      wannier_center(:,:)=0.0D0            
      Do i_band=1,n_occ       
       tmp(:)=0.0D0          
       Do ik=1,NTK        
        Do ibvec=1,NB           
         tmp(:) 
     +  =tmp(:) 
     +  +WEIGHT_b(ibvec) 
     +  *VEC_b(:,ibvec) 
     +  *AIMAG(LOG(OVERLAP(i_band,i_band,ik,ibvec)))
        ENDDO!ibvec            
       ENDDO!ik       
       wannier_center(:,i_band)=-tmp(:)/DBLE(NTK)   
      ENDDO!i_band           
!--
!OPEN(134,W,FILE='dat.wan-center') 
      OPEN(134,FILE='./dir-wan/dat.wan-center') 
      rewind(134)
      write(134,'(a)')'#Wannier center'
      write(134,'(a)')'#1:x, 2:y, 3:z (in xyz coord)'
      do i_band=1,n_occ 
       write(134,"(3F20.10)")(wannier_center(i,i_band),i=1,3)  
      enddo 
!-- 
!OPEN(113,W,FILE='dat.wan',FORM='unformatted') 
      OPEN(113,FILE='./dir-wan/dat.wan',FORM='unformatted') 
      allocate(C0_TMP_1(NTG,NTB))    
      allocate(C0_TMP_2(NTG,NTB))    
      REWIND(113)       
      write(113) n_occ!20170331(n_occ=NWF) 
      do ik=1,NTK 
       A_TMP(:,:)=0.0D0 
       do ib=1,N_BAND(ik)        
        do jb=1,n_occ        
         SUM_CMPX=0.0D0             
         do k=1,n_occ       
          do l=1,n_occ         
           SUM_CMPX
     +    =SUM_CMPX
     +    +C_MAT(ib,k,ik)*PSEUDO_MAT(k,l,ik)*A_MAT(l,jb,ik)
          enddo!l 
         enddo!k 
         A_TMP(ib,jb)=SUM_CMPX 
        enddo!jb
       enddo!ib
       C0_TMP_1(:,:)=C0(:,:,ik)        
       C0_TMP_2(:,:)=0.0D0     
!--
!write BF 
!       A_TMP(:,:)=0.d0 
!       do k_band=1,N_BAND(ik)
!        A_TMP(k_band,k_band)=1.0d0
!       enddo 
!--
       do ig=1,NG0(ik)  
        do jb=1,n_occ        
         SUM_CMPX=0.0D0             
         do kb=1,N_BAND(ik)       
          SUM_CMPX 
     +   =SUM_CMPX 
     +   +C0_TMP_1(ig,kb+N_BAND_BTM(ik))*A_TMP(kb,jb) 
         enddo!kb 
         C0_TMP_2(ig,jb)=SUM_CMPX  
        enddo!jb 
       enddo!ig             
       write(113)((C0_TMP_2(ig,ib),ig=1,NG0(ik)),ib=1,n_occ) 
       !--
       !20180519
       !do ib=1,n_occ
       ! write(113) (C0_TMP_2(ig,ib),ig=1,NG0(ik))
       !enddo!ib 
       !--
      enddo!ik              
!-- 
!OPEN(113,R,FILE='dat.wan',FORM='unformatted') 
      ierr=CHDIR("./dir-wan") 
      call system('pwd') 
      REWIND(113)       
      read(113) NWF 
      allocate(C0WN(NTG,NWF,NTK));C0WN(:,:,:)=0.0D0 
      do ik=1,NTK 
       read(113)((C0WN(ig,jb,ik),ig=1,NG0(ik)),jb=1,NWF) 
       !--
       !20180519
       !do jb=1,NWF 
       ! read(113)(C0WN(ig,jb,ik),ig=1,NG0(ik))
       !enddo!jb 
       !--
      enddo!ik
      write(6,*)'FINISH REDING C0_WN'
      ierr=CHDIR("..") 
      call system('pwd') 
!--
!OPEN(150,W,FILE='dat.umat')
      OPEN(150,FILE='./dir-wan/dat.umat') 
      rewind(150) 
      Mb=maxval(N_BAND) 
      write(6,*)"Mb=",Mb 
      write(150,*) n_occ!20170331(n_occ=NWF) 
      allocate(UNT(Mb,n_occ,NTK));UNT(:,:,:)=0.0d0 
      do ik=1,NTK 
       do ib=1,N_BAND(ik)
        do jb=1,n_occ 
         SUM_CMPX=0.0d0 
         do ig=1,NG0(ik)  
          SUM_CMPX 
     +   =SUM_CMPX
     +   +conjg(C0(ig,ib+N_BAND_BTM(ik),ik))*C0WN(ig,jb,ik)
         enddo 
         UNT(ib,jb,ik)=SUM_CMPX 
        enddo 
       enddo 
      enddo 
      deallocate(C0WN) 
      do ik=1,NTK 
       do ib=1,N_BAND(ik)
        write(150,*)(UNT(ib,jb,ik),jb=1,n_occ)
       enddo 
      enddo 
!--
!      WRITE(6,*) 
!      WRITE(6,*)'============================='
!      WRITE(6,*)'UNITALITY CHECK OF UMAT:<n|m>'
!      WRITE(6,*)'============================='
!      WRITE(6,*)  
!      do ik=1,NTK 
!       nm=n_occ 
!       allocate(mat_tmp(nm,nm));mat_tmp=0.0d0 
!       do iw=1,nm 
!        do jw=1,nm 
!         SUM_CMPX=0.0D0            
!         do ib=1,N_BAND(ik) 
!          SUM_CMPX=SUM_CMPX+CONJG(UNT(ib,iw,ik))*UNT(ib,jw,ik) 
!         enddo!ib  
!         mat_tmp(iw,jw)=SUM_CMPX 
!        enddo!jw  
!       enddo!iw   
!       do iw=1,nm 
!        write(6,'(300f10.5)')(mat_tmp(iw,jw),jw=1,nm)  
!       enddo 
!       write(6,*)       
!       deallocate(mat_tmp) 
!      enddo!ik        
!--
      WRITE(6,*) 
      WRITE(6,*)'============================='
      WRITE(6,*)'UNITALITY CHECK OF UMAT:<a|b>'
      WRITE(6,*)'============================='
      WRITE(6,*)
      do ik=1,NTK 
       nm=N_BAND(ik) 
       allocate(mat_tmp(nm,nm));mat_tmp=0.0d0 
       do ib=1,nm 
        do jb=1,nm 
         SUM_CMPX=0.0D0            
         do iw=1,n_occ       
          SUM_CMPX=SUM_CMPX+CONJG(UNT(ib,iw,ik))*UNT(jb,iw,ik) 
         enddo!iw  
         mat_tmp(ib,jb)=SUM_CMPX 
        enddo!jb  
       enddo!ib   
       !--
       SUM_CMPX=0.0d0 
       do ib=1,nm 
        SUM_CMPX=SUM_CMPX+mat_tmp(ib,ib) 
       enddo!ib 
       write(6,'(a8,x,i8,2f15.10)')'ik=',ik,SUM_CMPX 
       write(6,'(300f10.5)')(mat_tmp(ib,ib),ib=1,nm)  
       write(6,*)       
       !--
       if(N_BAND_inner(ik)/=0)then  
        do ib=1,nm 
         write(6,'(300f10.5)')(mat_tmp(ib,jb),jb=1,nm)  
        enddo 
        write(6,*)       
       endif 
       !--
       deallocate(mat_tmp) 
      enddo!ik        

!====================================
!=== INTERPOLATED BAND DISPERSION ===
!====================================

      allocate(H_MAT_K(n_occ,n_occ,NTK));H_MAT_K(:,:,:)=0.0D0      
      allocate(E_TMP(n_occ));E_TMP(:)=0.0d0            
      allocate(H_TMP_IN(n_occ,n_occ));H_TMP_IN(:,:)=0.0d0           
      allocate(H_TMP_OUT(n_occ,n_occ));H_TMP_OUT(:,:)=0.0d0  
      allocate(H_MAT_R(n_occ,n_occ,-Na1:Na1,-Na2:Na2,-Na3:Na3)) 
      allocate(WEIGHT_R(-Na1:Na1,-Na2:Na2,-Na3:Na3)) 
!--
!(1) SAMPLE k-POINTS FOR INTERPOLATED BAND DISPERSION 
      NSK_BAND_DISP=Ndiv*(N_sym_points-1)+1
      allocate(SK_BAND_DISP(3,NSK_BAND_DISP));SK_BAND_DISP(:,:)=0.0d0 
      allocate(E_BAND_DISP(n_occ,NSK_BAND_DISP));E_BAND_DISP(:,:)=0.0d0 
      call makekpts(Ndiv,N_sym_points,NSK_BAND_DISP,SK_sym_pts(1,1),
     +              SK_BAND_DISP(1,1))
      write(6,*) 
      write(6,*)'====================================='
      write(6,*)' CALCULATED kpts FOR BAND DISPERSION '
      write(6,*)'====================================='
      write(6,*) 
      write(6,'(a20,i15)')'NSK_BAND_DISP=',NSK_BAND_DISP
      write(6,*) 
      do ik=1,NSK_BAND_DISP 
       write(6,*) SK_BAND_DISP(:,ik)
      enddo 
!--
!(2) H_MAT(k) IN WANNIER BASIS 
      H_MAT_K(:,:,:)=0.0D0            
      Do ik=1,NTK 
       Do i_band=1,n_occ   
        Do j_band=1,n_occ   
         SUM_CMPX=0.0D0            
         Do k_band=1,n_occ       
          SUM_CMPX 
     +   =SUM_CMPX 
     +   +CONJG(A_MAT(k_band,i_band,ik))        
     +   *PSEUDO_EIG(k_band,ik)             
     +   *A_MAT(k_band,j_band,ik) 
         ENDDO 
         H_MAT_K(i_band,j_band,ik)=SUM_CMPX        
        ENDDO!j_band    
       ENDDO!i_band    
      ENDDO!ik        
!--
!(3) H_MAT(R) IN WANNIER BASIS 
      H_MAT_R(:,:,:,:,:)=0.0D0            
      Do ia1=-Na1,Na1!-1
       Do ia2=-Na2,Na2!-1
        Do ia3=-Na3,Na3!-1
         Do i_band=1,n_occ          
          Do j_band=1,n_occ          
           SUM_CMPX=0.0D0                             
           Do ik=1,NTK           
            PHASE=tpi*(SK0(1,ik)*DBLE(ia1) 
     +                +SK0(2,ik)*DBLE(ia2) 
     +                +SK0(3,ik)*DBLE(ia3))                
            PHASE_FACTOR=EXP(-ci*PHASE)         
            SUM_CMPX 
     +     =SUM_CMPX 
     +     +H_MAT_K(i_band,j_band,ik) 
     +     *PHASE_FACTOR                
           ENDDO               
           H_MAT_R(i_band,j_band,ia1,ia2,ia3)=SUM_CMPX/DBLE(NTK)     
          ENDDO!j_band
         ENDDO!i_band
        ENDDO!ia3      
       ENDDO!ia2       
      ENDDO!ia1         
!-- 
!(4) PRINT H_MAT(R) IN WANNIER BASIS 
      !
      !OPEN(122,W,FILE='dat.h_mat_r') 
      !
      OPEN(122,FILE='./dir-wan/dat.h_mat_r') 
      REWIND(122)
      write(122,'(a)')'#transfer'
      write(122,'(a)')'#1:R1, 2:R2, 3:R3 (lattice vector)'
      write(122,'(a)')'#1:i, 2:j, 3:Re(t_ij) [eV], 4:Im(t_ij) [eV]' 
      do ia1=-Na1,Na1!-1
       do ia2=-Na2,Na2!-1
        do ia3=-Na3,Na3!-1
         write(122,*) ia1,ia2,ia3           
         do ib=1,n_occ          
          do jb=1,n_occ  
           write(122,'(i5,i5,2f20.10)') 
     +     ib,jb,H_MAT_R(ib,jb,ia1,ia2,ia3)*au 
          enddo!jb 
         enddo!ib 
         write(122,*) 
        enddo 
       enddo 
      enddo 
      CLOSE(122) 
      !--
      write(6,*) 
      write(6,*)'============='
      write(6,*)'H_MAT_R in eV'
      write(6,*)'============='
      write(6,*) 
      do ia1=-Na1,Na1!-1
       do ia2=-Na2,Na2!-1
        do ia3=-Na3,Na3!-1
         write(6,*) ia1,ia2,ia3           
         do ib=1,n_occ          
          write(6,'(300f15.8)')
     +    (DBLE(H_MAT_R(ib,jb,ia1,ia2,ia3)*au),jb=1,n_occ)
         enddo 
         write(6,*) 
        enddo      
       enddo      
      enddo      
!--
!(5) WEIGHT_R 20170317 by Y.Nomura 
      write(6,'(a40)')'WEIGHT CALCULATED' 
      WEIGHT_R=1.0d0; SUM_REAL=0.0d0 
      do ia1=-Na1,Na1!-1         
       do ia2=-Na2,Na2!-1         
        do ia3=-Na3,Na3!-1         
         if(abs(ia1)==Na1.and.mod(nkb1,2)==0.and.Na1/=0) 
     +    WEIGHT_R(ia1,ia2,ia3)=WEIGHT_R(ia1,ia2,ia3)*0.5d0 
         if(abs(ia2)==Na2.and.mod(nkb2,2)==0.and.Na2/=0) 
     +    WEIGHT_R(ia1,ia2,ia3)=WEIGHT_R(ia1,ia2,ia3)*0.5d0 
         if(abs(ia3)==Na3.and.mod(nkb3,2)==0.and.Na3/=0) 
     +    WEIGHT_R(ia1,ia2,ia3)=WEIGHT_R(ia1,ia2,ia3)*0.5d0 
         SUM_REAL=SUM_REAL+WEIGHT_R(ia1,ia2,ia3)
        enddo
       enddo
      enddo 
      write(6,'(a40,f15.8,i8)')'SUM_WEIGHT,NTK',SUM_REAL,NTK  
      if(abs(SUM_REAL-dble(NTK))>1.0d-6)then 
       write(6,'(a40)')'SUM_WEIGHT/=NTK'
       stop 
      endif 
!--
!(6) OUTPUT Hphi/mVMC DATA 
      write(6,*) 
      write(6,'(a40)')'+++ write transfers in dir-mvmc +++'
      write(6,*) 
      call system('rm -rf dir-mvmc') 
      call system('mkdir dir-mvmc') 
      !
      call wrt_mvmc_hr(Na1,Na2,Na3,n_occ,H_MAT_R(1,1,-Na1,-Na2,-Na3),
     +     WEIGHT_R(-Na1,-Na2,-Na3)) 
      call wrt_mvmc_wcenter(n_occ,a1(1),a2(1),a3(1),wannier_center(1,1))
      call wrt_mvmc_SK_BAND_DISP(Ndiv,N_sym_points,NSK_BAND_DISP,
     +     SK_BAND_DISP(1,1)) 
      call wrt_mvmc_SK0(NTK,SK0(1,1)) 
      call wrt_mvmc_ef(FermiEnergy) 
      !
      !ELECTRON NUMBER == 0.0 -> DMX IS EVALUATED WITH BAND-CALC EF
      !
      if(elec_num==0.0d0)then 
       call calculate_dmx(NTB,NTK,NWF,nkb1,nkb2,nkb3,Na1,Na2,Na3,!NWF=n_occ 
     +      FermiEnergy,a1(1),a2(1),a3(1),b1(1),b2(1),b3(1),SK0(1,1),
     +      E_EIG(1,1),UNT(1,1,1)) 
      endif 
      deallocate(E_EIG) 
!--
!(7) H_MAT(k') IN WANNIER BASIS. AND DIAGONALIZE 
      do ik=1,NSK_BAND_DISP          
       H_TMP_IN(:,:)=0.0D0               
       do ib=1,n_occ      
        do jb=1,n_occ      
         SUM_CMPX=0.0D0                    
         do ia1=-Na1,Na1!-1         
          do ia2=-Na2,Na2!-1         
           do ia3=-Na3,Na3!-1         
            !--
            !-- NEAREST R SEARCH BY YOSHIHIDE YOSHIMOTO 
            !--
            call search_Rmin(ia1,ia2,ia3,nkb1,nkb2,nkb3,
     +      a1(1),a2(1),a3(1),ia1min,ia2min,ia3min)
            PHASE=tpi*(SK_BAND_DISP(1,ik)*DBLE(ia1min) 
     +                +SK_BAND_DISP(2,ik)*DBLE(ia2min) 
     +                +SK_BAND_DISP(3,ik)*DBLE(ia3min))                
            !--
            PHASE_FACTOR=EXP(ci*PHASE)*WEIGHT_R(ia1,ia2,ia3)!20170317 
            SUM_CMPX=SUM_CMPX+H_MAT_R(ib,jb,ia1,ia2,ia3)*PHASE_FACTOR 
           enddo!ia3            
          enddo!ia2                       
         enddo!ia1
         H_TMP_IN(ib,jb)=SUM_CMPX           
        enddo!jb
       enddo!ib
       !--
       E_TMP(:)=0.0D0                
       call diagV(n_occ,H_TMP_IN(1,1),E_TMP(1)) 
       H_TMP_OUT(:,:)=H_TMP_IN(:,:)
       E_BAND_DISP(:,ik)=E_TMP(:)            
      ENDDO!ik           
!--
!      write(6,*) 
!      do ik=1,NSK_BAND_DISP          
!       do i_band=1,n_occ 
!        write(6,*) E_BAND_DISP(i_band,ik)*27.21151D0  
!       enddo 
!      enddo 
!--
!(8) WRITE BAND DISPERSION 
      allocate(DIST_K(NSK_BAND_DISP));DIST_K(:)=0.0d0 
      allocate(dist(0:N_sym_points-1));dist(:)=0.0d0 
      !--
      dist(0)=0.0d0 
      do ix=1,N_sym_points-1 
       ks=(ix-1)*Ndiv+1
       ke=(ix)*Ndiv+1
       do ik=ks,ke 
        DIST_B(:)
     + =(SK_BAND_DISP(1,ik)-SK_BAND_DISP(1,ks))*b1(:)
     + +(SK_BAND_DISP(2,ik)-SK_BAND_DISP(2,ks))*b2(:)
     + +(SK_BAND_DISP(3,ik)-SK_BAND_DISP(3,ks))*b3(:)
        DIST_KSPACE
     + =DSQRT(DIST_B(1)**2+DIST_B(2)**2+DIST_B(3)**2)
        DIST_K(ik)=DIST_KSPACE+dist(ix-1) 
       enddo!ik 
       dist(ix)=dist(ix-1)+DIST_KSPACE
      enddo!ix 
      !
      !OPEN(114,W,FILE='dat.iband') 
      !
      OPEN(114,FILE='./dir-wan/dat.iband') 
      write(114,'(a)')'#Wannier interpolaed band'
      write(114,'(a)')'#1:k, 2:Energy [eV]' 
      REVERSE=.TRUE.        
      do ib=1,n_occ             
       if(REVERSE) then 
        do ik=1,NSK_BAND_DISP                     
         write(114,*) DIST_K(ik)/DIST_K(NSK_BAND_DISP),
     +                E_BAND_DISP(ib,ik)*au
        enddo!ik        
        REVERSE=.FALSE.        
       else         
        do ik=NSK_BAND_DISP,1,-1          
         write(114,*) DIST_K(ik)/DIST_K(NSK_BAND_DISP),
     +                E_BAND_DISP(ib,ik)*au
        enddo!ik        
        REVERSE=.TRUE.        
       endif!REVERSE                   
      enddo!ib                      
!--
!(9) WRITE FERMISURFACE WITH WANNIER INTERPOLATION 
      IF(ALL(dense(1:3)/=0)) THEN
         CALL wrt_frmsf_wan(Na1,Na2,Na3,nkb1,nkb2,nkb3,
     +                      a1,a2,a3,b1,b2,b3,FermiEnergy,
     +                      WEIGHT_R,H_MAT_R)
      ENDIF

!============================================================
!=== DENSITY OF STATE and DENSITY MATRIX IN WANNIER SPACE ===
!============================================================
   
      write(6,*) 
      write(6,'(a40)')'==============================='
      write(6,'(a40)')' CALCULATE WANNIER DOS AND DMX '
      write(6,'(a40)')'==============================='
      write(6,*) 
      allocate(E_EIG(NWF,NTK)); E_EIG(:,:)=0.0d0
      allocate(V_EIG(NWF,NWF,NTK)); V_EIG(:,:,:)=0.0d0
      call calculate_eigenstate(NWF,Na1,Na2,Na3,
     +     H_MAT_R(1,1,-Na1,-Na2,-Na3),NTK,nkb1,nkb2,nkb3,
     +     a1(1),a2(1),a3(1),1,NTK,SK0(1,1),E_EIG(1,1),V_EIG(1,1,1))   
      filename='./dir-wan/dat.dos.wannier'
      !
      !ELECTRON NUMBER == 0.0 -> DOS IS EVALUATED WITH BAND-CALC EF
      !
      if(elec_num==0.0d0)then 
       call calculate_dos(NWF,NTK,nkb1,nkb2,nkb3,elec_num,FermiEnergy,
     +      filename,b1(1),b2(1),b3(1),SK0(1,1),E_EIG(1,1))
      endif 
      !
      !ELECTRON NUMBER /= 0.0 -> DOS & DMX ARE EVALUATED WITH RECALC EF
      !
      if(elec_num/=0.0d0)then  
       call calculate_dos(NWF,NTK,nkb1,nkb2,nkb3,elec_num,FermiEnergy,
     +      filename,b1(1),b2(1),b3(1),SK0(1,1),E_EIG(1,1))
       call calculate_dmx(NWF,NTK,NWF,nkb1,nkb2,nkb3,Na1,Na2,Na3,!NWF=n_occ 
     +      FermiEnergy,a1(1),a2(1),a3(1),b1(1),b2(1),b3(1),SK0(1,1),
     +      E_EIG(1,1),V_EIG(1,1,1)) 
      endif 
      deallocate(E_EIG,V_EIG)

!==================================
!=== VISUALIZE WANNIER FUNCTION ===
!==================================

      if(CALC_REAL_SPACE_WANNIER)then 
       !--
       write(6,*) 
      write(6,'(a40)')'============================='
      write(6,'(a40)')' CALCULATE REALSPACE WANNIER '
      write(6,'(a40)')'============================='
       write(6,*) 
       nfft1=nwx2+1;nfft2=nwy2+1;nfft3=nwz2+1
       Nl123=nfft1*nfft2*nfft3 
       write(6,'(a40,i5,i5,i5)')'nwx2,nwy2,nwz2',nwx2,nwy2,nwz2 
       write(6,'(a40,i5,i5,i5)')'nfft1,nfft2,nfft2',nfft1,nfft2,nfft3 
       write(6,'(a40,i10)')'Nl123',Nl123 
       call fft3_init(nwx2,nwy2,nwz2,nfft1,nfft2,nfft3,fs) 
       !--
       mem_size=dble((xmax-xmin)*nwx2+(ymax-ymin)*nwy2+(zmax-zmin)*nwz2)
       mem_size=mem_size*16.0d0/1024.0d0/1024.0d0/1024.0d0 
       write(6,'(a40,f20.15)')'mem size realspace wannier (GB)',mem_size
       if(mem_size>1.0d0)then 
        write(6,'(a40,f20.15)')'mem size > 1 GB; skipp' 
        go to 9999
       endif  
       !--
       allocate(WANNIER_REALSPACE((xmin)*nwx2:(xmax)*nwx2-1,
     +                            (ymin)*nwy2:(ymax)*nwy2-1, 
     +                            (zmin)*nwz2:(zmax)*nwz2-1))
       !--
       !OPEN(113,R,FILE='dat.wan',FORM='unformatted') 
       ierr=CHDIR("./dir-wan") 
       call system('pwd') 
       REWIND(113)       
       read(113) NWF!20170331(n_occ=NWF) 
       allocate(C0WN(NTG,NWF,NTK));C0WN(:,:,:)=0.0D0 
       do ik=1,NTK 
        read(113)((C0WN(ig,jb,ik),ig=1,NG0(ik)),jb=1,NWF) 
        !
        !20180519
        !do jb=1,NWF 
        ! read(113)(C0WN(ig,jb,ik),ig=1,NG0(ik))
        !enddo!jb 
        !
       enddo!ik 
       !write(6,*)'FINISH REDING C0_WN'
       !--
       do ix=1,N_write_wannier!NWF 
        iw=wrt_list(ix)!20170406  
        WANNIER_REALSPACE=0.0d0           
!$OMP PARALLEL PRIVATE(ik,C0_I,wfunc,fftwk,pw) 
        allocate(C0_I(NTG));C0_I=0.0D0 
        allocate(fftwk(Nl123*2),stat=err) 
        allocate(wfunc(Nl123*2),stat=err) 
        allocate(pw((xmin)*nwx2:(xmax)*nwx2-1,
     +              (ymin)*nwy2:(ymax)*nwy2-1, 
     +              (zmin)*nwz2:(zmax)*nwz2-1));pw=0.0d0           
!$OMP DO 
        do ik=1,NTK
         C0_I(:)=C0WN(:,iw,ik) 
         call make_pw(C0_I(1),wfunc(1),fftwk(1),NG0(ik),KG0(1,1,ik),NTG,
     +    nwx2,nwy2,nwz2,nfft1,nfft2,Nl123,fs,
     +    xmin,xmax,ymin,ymax,zmin,zmax,SK0(1,ik),
     +    pw(xmin*nwx2,ymin*nwy2,zmin*nwz2)) 
        enddo!ik         
!$OMP END DO
!$OMP CRITICAL
        WANNIER_REALSPACE=WANNIER_REALSPACE+pw 
!$OMP END CRITICAL
        deallocate(C0_I,fftwk,wfunc,pw)
!$OMP END PARALLEL
        !--
        !normalize 
        WANNIER_REALSPACE=WANNIER_REALSPACE/DBLE(NTK)/DSQRT(VOLUME)
        !--
        !OUTPUT by VESTA format 
        !OPEN(116,W,FILE='dat.wan-realspace-xxx.grd')
        write(filename,"('dat.wan-realspace-',i3.3,'.grd')") iw 
        OPEN(116,FILE=filename) 
        REWIND(116)
        na_grid=nwx2 
        nb_grid=nwy2
        nc_grid=nwz2 
        aa1=a1*dble(xmax-xmin)
        aa2=a2*dble(ymax-ymin)
        aa3=a3*dble(zmax-zmin)
        call est_latparam(aa1(1),aa2(1),aa3(1),a,b,c,alp,bet,gmm) 
        write(116,*)'Electron density'
        write(116,'(6f15.10)')
     +  a*bohr,b*bohr,c*bohr,alp,bet,gmm!bohr->angstrom
        write(116,'(3I5)') 
     +  (xmax-xmin)*na_grid,(ymax-ymin)*nb_grid,(zmax-zmin)*nc_grid 
        write(116,"(6g25.16)") (((dble(WANNIER_REALSPACE(i1,i2,i3)),
     +  i3=(zmin)*nc_grid,(zmax)*nc_grid-1),
     +  i2=(ymin)*nb_grid,(ymax)*nb_grid-1),
     +  i1=(xmin)*na_grid,(xmax)*na_grid-1)
        CLOSE(116)
       enddo!ix 
       !--
       !gencif 20180906  
       call gencif(a1(1),a2(1),a3(1),xmin,xmax,ymin,ymax,zmin,zmax,
     +  nsi,chemical_species(1),asi(1,1))
       !--
       ierr=CHDIR("..") 
       call system('pwd') 
9999  endif!CALC_REAL_SPACE_WANNIER 
      !--
      STOP 
      END        
