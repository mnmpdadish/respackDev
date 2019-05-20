      PROGRAM MAIN
!--
      use m_rdinput       
      use fft_3d 
      include "config_w3d.h" 
      call read_input 
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
      write(6,*)'FINISH REDING a1,a2,a3'
      CLOSE(105)
      call OUTER_PRODUCT(a2(1),a3(1),b1(1))
      VOLUME=a1(1)*b1(1)+a1(2)*b1(2)+a1(3)*b1(3)
      b1(:)=b1(:)*tpi/VOLUME 
      call OUTER_PRODUCT(a3(1),a1(1),b2(1))
      b2(:)=b2(:)*tpi/VOLUME 
      call OUTER_PRODUCT(a1(1),a2(1),b3(1))
      b3(:)=b3(:)*tpi/VOLUME 
      write(6,*) 
      write(6,*)'======================'
      write(6,*)'CHECK OF (A_VEC*B_VEC)'
      write(6,*)'======================'
      write(6,*) 
      write(6,*)'a1*b1=',b1(1)*a1(1)+b1(2)*a1(2)+b1(3)*a1(3)
      write(6,*)'a2*b2=',b2(1)*a2(1)+b2(2)*a2(2)+b2(3)*a2(3)
      write(6,*)'a3*b3=',b3(1)*a3(1)+b3(2)*a3(2)+b3(3)*a3(3)
      write(6,*)'VOLUME OF UNIT CELL=',VOLUME  
!--
!OPEN(100,R,FILE='dat.symmetry') 
      OPEN(100,FILE='./dir-wfn/dat.symmetry') 
      rewind(100) 
      read(100,*) nsymq 
      read(100,*) nnp
      allocate(rg(3,3,nsymq));rg=0
      allocate(pg(3,nsymq));pg=0
      allocate(rginv(3,3,nsymq));rginv(:,:,:)=0.0d0 
      do iop=1,nsymq
       read(100,*)((rg(i,j,iop),i=1,3),j=1,3) 
       read(100,*)(pg(i,iop),i=1,3)   
      enddo 
      CLOSE(100) 
      rginv=rg 
      do iop=1,nsymq
      call invmat(3,rginv(1,1,iop)) 
      enddo 
      write(6,*)'finish rg'
!check 
      do iop=1,nsymq
       write(6,*) iop
       do i=1,3
        write(6,'(3I5,1x,3F15.10)')
     +  (rg(i,j,iop),j=1,3),(rginv(i,j,iop),j=1,3)
       enddo 
      enddo 
!do iop=1,nsymq
!do i=1,3
!do j=1,3
!s=0.0d0 
!do k=1,3
!s=s+rg(i,k,iop)*rginv(k,j,iop)
!enddo 
!write(6,*) i,j,s
!enddo 
!enddo 
!enddo 
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
      write(6,*)'finish ski'
!20170327 
      !call est_NTK(NK_irr,SKI(1,1),NTK,nkb1,nkb2,nkb3)  
      !write(6,'(a24,4i10)')'nkb1,nkb2,nkb3,NTK=',nkb1,nkb2,nkb3,NTK  
      !Na1=nkb1/2;Na2=nkb2/2;Na3=nkb3/2
      !write(6,'(a24,3i10)')'Na1,Na2,Na3=',Na1,Na2,Na3 
!--
!OPEN(132,R,FILE='dat.nkm')20170419  
      OPEN(132,FILE='./dir-wfn/dat.nkm') 
      allocate(NGI(Nk_irr));NGI(:)=0
!     allocate(NBI(Nk_irr));NBI(:)=0!20170419
      rewind(132)
      do ik=1,Nk_irr 
       read(132,*) NGI(ik)!,NBI(ik)!20170419
      enddo 
      close(132) 
      NTG=maxval(abs(NGI(:))) 
!     NTB=maxval(abs(NBI(:)))!20170419  
      write(6,'(a10,i10)')'NTG=',NTG  
!     write(6,'(a10,i10)')'NTB=',NTB!20170419
!--
!OPEN(104,R,FILE='dat.kg') 
      OPEN(104,FILE='./dir-wfn/dat.kg') 
      rewind(104) 
      allocate(KGI(3,NTG,Nk_irr));KGI(:,:,:)=0 
      do ik=1,Nk_irr 
       read(104,*) NG_for_psi 
       NGI(IK)=NG_for_psi 
       do ig=1,NG_for_psi 
        read(104,*)(KGI(i,ig,ik),i=1,3) 
       enddo 
      enddo  
      close(104)
      write(6,*)'finish kgi'
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
      L1=maxval(abs(KGI(1,:,:)))+1;write(6,*)'L1=',L1 
      L2=maxval(abs(KGI(2,:,:)))+1;write(6,*)'L2=',L2 
      L3=maxval(abs(KGI(3,:,:)))+1;write(6,*)'L3=',L3 
      allocate(packing(-L1:L1,-L2:L2,-L3:L3,Nk_irr))   
      packing(:,:,:,:)=0 
      do ik=1,Nk_irr 
       do ig=1,NGI(ik) 
        i1=KGI(1,ig,ik);j1=KGI(2,ig,ik);k1=KGI(3,ig,ik) 
        packing(i1,j1,k1,ik)=ig 
       enddo 
      enddo 
!--
!fft grid
      !20180309 
      !m1=maxval(abs(KGI(1,:,:))) 
      !m2=maxval(abs(KGI(2,:,:))) 
      !m3=maxval(abs(KGI(3,:,:))) 
      !nwx2=algn235(2*m1) 
      !nwy2=algn235(2*m2) 
      !nwz2=algn235(2*m3) 
      !
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
      !
      !nwx2=algn235(int(qwf/d1)+1,1) 
      !nwy2=algn235(int(qwf/d2)+1,1) 
      !nwz2=algn235(int(qwf/d3)+1,1) 
      !
      nwx2=algn235(int(qwf/d1)+1)!20180412 by Y.Nomura 
      nwy2=algn235(int(qwf/d2)+1)!20180412 by Y.Nomura 
      nwz2=algn235(int(qwf/d3)+1)!20180412 by Y.Nomura 
      !
      nfft1=nwx2+1
      nfft2=nwy2+1
      nfft3=nwz2+1
      Nl123=nfft1*nfft2*nfft3 
      write(6,*)'nwx2=',nwx2 
      write(6,*)'nwy2=',nwy2 
      write(6,*)'nwz2=',nwz2 
      write(6,*)'nfft1=',nfft1 
      write(6,*)'nfft2=',nfft2 
      write(6,*)'nfft3=',nfft3 
      write(6,*)'NL123=',Nl123  
      call fft3_init(nwx2,nwy2,nwz2,nfft1,nfft2,nfft3,fs) 
!--
!20180302  
      call est_NTK(Nk_irr,Nsymq,SKI(1,1),rg(1,1,1),NTK)
      write(6,*)'Estimated NTK=',NTK 
!--
!gen(SK0,numirr,numrot,trs,RW) 
      allocate(SK0(3,NTK));SK0(:,:)=0.0d0
      allocate(numirr(NTK));numirr(:)=0
      allocate(numrot(NTK));numrot(:)=0
      allocate(trs(NTK));trs(:)=0
      allocate(RW(3,NTK));RW(:,:)=0
!20161207 
!      do ik=1,Nk_irr 
!       SK0(:,ik)=SKI(:,ik) 
!       numirr(ik)=ik;numrot(ik)=1;trs(ik)=1;RW(1:3,ik)=0
!      enddo 
!      jk=Nk_irr 
!--
      jk=0
      do ik=1,Nk_irr
      do iop=1,Nsymq
!sym
       ktmp(:)=0.0d0;RWtmp(:)=0  
       ktmp(1)=rg(1,1,iop)*SKI(1,ik)
     +        +rg(1,2,iop)*SKI(2,ik)
     +        +rg(1,3,iop)*SKI(3,ik)
       ktmp(2)=rg(2,1,iop)*SKI(1,ik)
     +        +rg(2,2,iop)*SKI(2,ik)
     +        +rg(2,3,iop)*SKI(3,ik)
       ktmp(3)=rg(3,1,iop)*SKI(1,ik)
     +        +rg(3,2,iop)*SKI(2,ik)
     +        +rg(3,3,iop)*SKI(3,ik)
       call kcheck(ktmp(1),RWtmp(1))!rewind check 
       do iik=1,jk
        if(abs(SK0(1,iik)-ktmp(1))<1.0d-4.and. 
     +     abs(SK0(2,iik)-ktmp(2))<1.0d-4.and. 
     +     abs(SK0(3,iik)-ktmp(3))<1.0d-4)goto 1000
       enddo!iik
       jk=jk+1
       SK0(:,jk)=ktmp(:)
       numirr(jk)=ik;numrot(jk)=iop;trs(jk)=1;RW(:,jk)=RWtmp(:)
!time-reversal
1000   ktmp(:)=0.0d0;RWtmp(:)=0  
       ktmp(1)=rg(1,1,iop)*SKI(1,ik)
     +        +rg(1,2,iop)*SKI(2,ik)
     +        +rg(1,3,iop)*SKI(3,ik)
       ktmp(2)=rg(2,1,iop)*SKI(1,ik)
     +        +rg(2,2,iop)*SKI(2,ik)
     +        +rg(2,3,iop)*SKI(3,ik) 
       ktmp(3)=rg(3,1,iop)*SKI(1,ik)
     +        +rg(3,2,iop)*SKI(2,ik)
     +        +rg(3,3,iop)*SKI(3,ik) 
       call kcheck_trs(ktmp(1),RWtmp(1))!rewind check 20170321 
       do iik=1,jk
        if(abs(SK0(1,iik)-(-ktmp(1)))<1.0d-4.and. 
     +     abs(SK0(2,iik)-(-ktmp(2)))<1.0d-4.and. 
     +     abs(SK0(3,iik)-(-ktmp(3)))<1.0d-4)goto 2000
       enddo!iik
       jk=jk+1
       SK0(:,jk)=-ktmp(:)
       numirr(jk)=ik;numrot(jk)=iop;trs(jk)=-1;RW(:,jk)=RWtmp(:)
2000  enddo 
      enddo 
!--
      if(NTK/=jk)then 
       write(6,*)'ERROR;STOP;NTK should be jk'   
       write(6,*)'NTK=',NTK,'jk=',jk;STOP
      endif 
!--   
       write(6,*) 
       write(6,*)'====='
       write(6,*)' SK0 '
       write(6,*)'====='
       write(6,*) 
       do ik=1,NTK
        write(6,'(I5,3F12.5)') ik,SK0(:,ik) 
       enddo 
!--
      call est_nkbi(NTK,SK0(1,1),nkb1,nkb2,nkb3)  
      write(6,*) 
      write(6,'(a24,4i10)')'nkb1,nkb2,nkb3,NTK=',nkb1,nkb2,nkb3,NTK  
      Na1=nkb1/2;Na2=nkb2/2;Na3=nkb3/2
      write(6,'(a24,3i10)')'Na1,Na2,Na3=',Na1,Na2,Na3 
!--
!gen(NG0,KG0) 
      allocate(KGtmp(3,NTG));KGtmp(:,:)=0 
      allocate(KG0(3,NTG,NTK));KG0(:,:,:)=0 
      allocate(NG0(NTK));NG0(:)=0
!20161207 
!      do jk=1,Nk_irr 
!       NG0(jk)=NGI(jk);KG0(:,:,jk)=KGI(:,:,jk)
!      enddo 
!      do jk=Nk_irr+1,NTK 
!--
      do jk=1,NTK 
       if(trs(jk)==1)then 
        ik=numirr(jk);iop=numrot(jk) 
        ktmp(1)=rg(1,1,iop)*SKI(1,ik)
     +         +rg(1,2,iop)*SKI(2,ik)
     +         +rg(1,3,iop)*SKI(3,ik)+dble(RW(1,jk)) 
        ktmp(2)=rg(2,1,iop)*SKI(1,ik)
     +         +rg(2,2,iop)*SKI(2,ik)
     +         +rg(2,3,iop)*SKI(3,ik)+dble(RW(2,jk))  
        ktmp(3)=rg(3,1,iop)*SKI(1,ik)
     +         +rg(3,2,iop)*SKI(2,ik)
     +         +rg(3,3,iop)*SKI(3,ik)+dble(RW(3,jk))  
        call make_KG0(NTG,b1(1),b2(1),b3(1),Ecut_for_psi,
     +       ktmp(1),ktmp(2),ktmp(3),KG0(1,1,jk),NG_for_psi)
        if(NG_for_psi/=NGI(ik))then 
         write(6,*)'ERROR;STOP;NG_for_psi should be NGI(ik)'   
         write(6,*)'NG_for_psi=',NG_for_psi,'NGI(ik)=',NGI(ik);STOP
        endif 
        NG0(jk)=NG_for_psi  
       elseif(trs(jk)==-1)then  
        ik=numirr(jk); iop=numrot(jk) 
        ktmp(1)=rg(1,1,iop)*SKI(1,ik)
     +         +rg(1,2,iop)*SKI(2,ik)
     +         +rg(1,3,iop)*SKI(3,ik)+dble(RW(1,jk)) 
        ktmp(2)=rg(2,1,iop)*SKI(1,ik)
     +         +rg(2,2,iop)*SKI(2,ik)
     +         +rg(2,3,iop)*SKI(3,ik)+dble(RW(2,jk))  
        ktmp(3)=rg(3,1,iop)*SKI(1,ik)
     +         +rg(3,2,iop)*SKI(2,ik)
     +         +rg(3,3,iop)*SKI(3,ik)+dble(RW(3,jk))  
        KGtmp(:,:)=0 
        call make_KG0(NTG,b1(1),b2(1),b3(1),Ecut_for_psi,
     +       ktmp(1),ktmp(2),ktmp(3),KGtmp(1,1),NG_for_psi)
        if(NG_for_psi/=NGI(ik))then 
         write(6,*)'ERROR;STOP;NG_for_psi should be NGI(ik)'   
         write(6,*)'NG_for_psi=',NG_for_psi,'NGI(ik)=',NGI(ik);STOP
        endif 
        NG0(jk)=NG_for_psi  
        KG0(:,:,jk)=-KGtmp(:,:) 
       endif 
      enddo 
!--
      write(6,*) 
      write(6,*)'====='
      write(6,*)' NG0 '
      write(6,*)'====='
      write(6,*) 
      do ik=1,NTK 
       write(6,'(1I10)') NG0(ik) 
      enddo 
!--
!OPEN(113,R,FILE='dat.wan',FORM='unformatted') 
      OPEN(113,FILE='./dir-wan/dat.wan',FORM='unformatted') 
      REWIND(113)       
      read(113) NWF 
      allocate(C0(NTG,NWF,NTK));C0(:,:,:)=0.0D0 
      do ik=1,NTK
       read(113)((C0(ig,iw,ik),ig=1,NG0(ik)),iw=1,NWF)           
       !
       !20180519 
       !do iw=1,NWF
       ! read(113)(C0(ig,iw,ik),ig=1,NG0(ik))
       !enddo!iw
       !
      enddo!ik
      write(6,*)'FINISH REDING C0'
      CLOSE(113)  
!--
!      do ik=1,Nk_irr 
!       do iw=1,NWF
!        do jw=1,NWF 
!         SUM_CMPX=0.0d0 
!         do ig=1,NG0(ik) 
!          SUM_CMPX=SUM_CMPX+CONJG(C0(ig,iw,ik))*C0(ig,jw,ik) 
!         enddo 
!         write(6,'(3i5,x,2f15.10)') ik,iw,jw,SUM_CMPX 
!        enddo 
!       enddo 
!      enddo 
!--
!OPEN(134,R,FILE='dat.wan-center') 
      OPEN(134,FILE='./dir-wan/dat.wan-center') 
      allocate(coord(3,NWF));coord(:,:)=0.0d0    
      rewind(134)
      read(134,'(a)') dum_char
      read(134,'(a)') dum_char
      do iw=1,NWF
       read(134,*)(coord(i,iw),i=1,3) 
      enddo 
!--
!OPEN(127,R,FILE='dat.sample-q')
!      OPEN(127,FILE='./dir-wfn/dat.sample-q')
!      Nq_irr=Nk_irr 
!      allocate(SQR(3,Nq_irr));SQR(:,:)=0.0D0 
!      rewind(127) 
!      read(127,*) idum 
!      do iq=1,Nq_irr 
!       read(127,*)(SQR(i,iq),i=1,3) 
!      enddo 
!      write(6,*)'FINISH REDING SQR'
!      CLOSE(127) 
!20161207 
!shift[-1/2:1/2]  
!      allocate(SQI(3,Nq_irr));SQI(:,:)=0.0D0 
!      do iq=1,Nq_irr
!       call search_q(SQR(1,iq),SQI(1,iq)) 
!      enddo  
!      write(6,*) 
!      write(6,*)'====================='
!      write(6,*)'shifted SQ (SQR->SQI)'
!      write(6,*)'====================='
!      write(6,*) 
!      do iq=1,Nq_irr 
!       write(6,'(a10,3f10.5,2x,3f10.5)')'SQR->SQI',SQR(:,iq),SQI(:,iq) 
!      enddo 
!--
!OPEN(300,R,FILE='dat.chi_cutoff')!20170402  
       OPEN(300,FILE='./dir-eps/dat.chi_cutoff')!20170402 
       read(300,*) Ecut_for_eps!20170402  
!--
!OPEN(135,R,FILE='W-grid') 
      OPEN(135,FILE='./dir-eps/dat.wgrid') 
      rewind(135) 
      read(135,*) Num_freq_grid!20170402
      ne=Num_freq_grid!20170402
      allocate(em(ne));em(:)=0.0d0 
      do ie=1,ne 
       !read(135,*) em(ie) 
        read(135,'(2f15.10)') em(ie)!20171212
      enddo 
!--
!OPEN(400,R,FILE='dat.log.400')!20170420  
!OPEN(301,R,FILE='dat.sq')!20170420   
      Nq_irr=Nk_irr 
      allocate(SQI(3,Nq_irr));SQI(:,:)=0.0D0 
      ierr=CHDIR("./dir-eps") 
      do iq=1,Nq_irr 
       write(dirname,"('q',i3.3)")iq 
       ierr=CHDIR(dirname) 
       inquire(file='dat.log.400',exist=file_e) 
       if(file_e)then 
        write(6,*)'dat.log.400 exists in ',trim(dirname) 
       else
        write(6,*)'error: no dat.log.400 in ',trim(dirname) 
       endif 
       OPEN(301,FILE='dat.sq') 
       rewind(301)
       inquire(file='dat.sq',exist=file_e) 
       if(file_e)then 
        read(301,*)(SQI(i,iq),i=1,3) 
       else
        write(6,*)'error: no dat.sq in',trim(dirname) 
       endif 
       close(301) 
       ierr=CHDIR("..") 
      enddo!iq 
      ierr=CHDIR("..") 
      call system('pwd') 
!--
!gen(SQ,numirrq,numrotq,trsq,RWq) 
      NTQ=NTK!TOTAL NUMBER OF Q POINTS 
      allocate(SQ(3,NTQ));SQ(:,:)=0.0d0 
      allocate(numirrq(NTQ));numirrq(:)=0 
      allocate(numrotq(NTQ));numrotq(:)=0 
      allocate(trsq(NTQ));trsq(:)=0
      allocate(RWq(3,NTQ));RWq(:,:)=0      
!--
      do iq=1,Nq_irr 
       SQ(:,iq)=SQI(:,iq) 
       numirrq(iq)=iq;numrotq(iq)=1;trsq(iq)=1;RWq(1:3,iq)=0
      enddo 
      jk=Nq_irr 
      do iq=1,Nq_irr
      do iop=1,Nsymq
!sym
       ktmp(:)=0.0d0;RWtmp(:)=0  
       ktmp(1)=rg(1,1,iop)*SQI(1,iq)
     +        +rg(1,2,iop)*SQI(2,iq)
     +        +rg(1,3,iop)*SQI(3,iq)
       ktmp(2)=rg(2,1,iop)*SQI(1,iq)
     +        +rg(2,2,iop)*SQI(2,iq)
     +        +rg(2,3,iop)*SQI(3,iq)
       ktmp(3)=rg(3,1,iop)*SQI(1,iq)
     +        +rg(3,2,iop)*SQI(2,iq)
     +        +rg(3,3,iop)*SQI(3,iq)
       call kcheck(ktmp(1),RWtmp(1))!rewind check 
       do iik=1,jk 
        if(abs(SQ(1,iik)-ktmp(1))<1.0d-4.and. 
     +     abs(SQ(2,iik)-ktmp(2))<1.0d-4.and. 
     +     abs(SQ(3,iik)-ktmp(3))<1.0d-4)goto 1100
       enddo!iik
       jk=jk+1
       SQ(:,jk)=ktmp(:)
       numirrq(jk)=iq;numrotq(jk)=iop;trsq(jk)=1;RWq(:,jk)=RWtmp(:)
!time-reversal
1100   ktmp(:)=0.0d0;RWtmp(:)=0  
       ktmp(1)=rg(1,1,iop)*SQI(1,iq)
     +        +rg(1,2,iop)*SQI(2,iq)
     +        +rg(1,3,iop)*SQI(3,iq)
       ktmp(2)=rg(2,1,iop)*SQI(1,iq)
     +        +rg(2,2,iop)*SQI(2,iq)
     +        +rg(2,3,iop)*SQI(3,iq) 
       ktmp(3)=rg(3,1,iop)*SQI(1,iq)
     +        +rg(3,2,iop)*SQI(2,iq)
     +        +rg(3,3,iop)*SQI(3,iq) 
       call kcheck_trs(ktmp(1),RWtmp(1))!rewind check 20170321 
       do iik=1,jk
        if(abs(SQ(1,iik)-(-ktmp(1)))<1.0d-4.and. 
     +     abs(SQ(2,iik)-(-ktmp(2)))<1.0d-4.and. 
     +     abs(SQ(3,iik)-(-ktmp(3)))<1.0d-4)goto 2100
       enddo!iik
       jk=jk+1
       SQ(:,jk)=-ktmp(:)
       numirrq(jk)=iq;numrotq(jk)=iop;trsq(jk)=-1;RWq(:,jk)=RWtmp(:)
2100  enddo 
      enddo 
!--
      if(NTQ/=jk)then 
       write(6,*)'ERROR;STOP;NTQ should be jk'   
       write(6,*)'NTQ=',NTQ,'jk=',jk;STOP
      endif 
!--
      write(6,*) 
      write(6,*)'============'
      write(6,*)' SQ and SK0 '
      write(6,*)'============'
      write(6,*) 
      do iq=1,NTQ 
       ik=iq
       write(6,'(I5,3F10.5,2x,3F10.5)') iq,SQ(:,iq),SK0(:,ik) 
      enddo 
      write(6,*) 
!--
!      do iq=1,NTQ 
!       write(6,'(i7,x,i7,x,i7,x,i7,x,3i5)')
!     + iq,numirrq(iq),numrotq(iq),trsq(iq),RWq(:,iq) 
!      enddo 
!      write(6,*) 
!--
!gen(NGQ_eps,NGQ_psi,LG0)
      allocate(LG0(3,NTG,NTQ));LG0(:,:,:)=0
      allocate(NGQ_eps(NTQ));NGQ_eps(:)=0
      allocate(NGQ_psi(NTQ));NGQ_psi(:)=0
!--
      do iq=1,Nq_irr 
       q1=SQI(1,iq);q2=SQI(2,iq);q3=SQI(3,iq)
       call make_LG0(NTG,b1(1),b2(1),b3(1),Ecut_for_eps,
     +      Ecut_for_psi,q1,q2,q3,LG0(1,1,iq),NG_for_eps,NG_for_psi) 
       NGQ_eps(iq)=NG_for_eps;NGQ_psi(iq)=NG_for_psi  
       write(6,'(i8,3f10.5,a8,i8,a8,i10)') 
     + iq,q1,q2,q3,'NGeps',NG_for_eps,'NGpsi',NG_for_psi  
      enddo 
      do iq=Nq_irr+1,NTQ 
       if(trsq(iq)==1)then 
        iqir=numirrq(iq);iop=numrotq(iq) 
        q1=rg(1,1,iop)*SQI(1,iqir)
     +    +rg(1,2,iop)*SQI(2,iqir)
     +    +rg(1,3,iop)*SQI(3,iqir)+dble(RWq(1,iq)) 
        q2=rg(2,1,iop)*SQI(1,iqir)
     +    +rg(2,2,iop)*SQI(2,iqir)
     +    +rg(2,3,iop)*SQI(3,iqir)+dble(RWq(2,iq))  
        q3=rg(3,1,iop)*SQI(1,iqir)
     +    +rg(3,2,iop)*SQI(2,iqir)
     +    +rg(3,3,iop)*SQI(3,iqir)+dble(RWq(3,iq))  
        call make_LG0(NTG,b1(1),b2(1),b3(1),Ecut_for_eps,
     +       Ecut_for_psi,q1,q2,q3,LG0(1,1,iq),NG_for_eps,NG_for_psi) 
        NGQ_eps(iq)=NG_for_eps;NGQ_psi(iq)=NG_for_psi  
        write(6,'(i8,3f10.5,a8,i8,a8,i10)') 
     +  iq,q1,q2,q3,'NGeps',NG_for_eps,'NGpsi',NG_for_psi  
       elseif(trsq(iq)==-1)then  
        iqir=numirrq(iq);iop=numrotq(iq) 
        q1=rg(1,1,iop)*SQI(1,iqir)
     +    +rg(1,2,iop)*SQI(2,iqir)
     +    +rg(1,3,iop)*SQI(3,iqir)+dble(RWq(1,iq)) 
        q2=rg(2,1,iop)*SQI(1,iqir)
     +    +rg(2,2,iop)*SQI(2,iqir)
     +    +rg(2,3,iop)*SQI(3,iqir)+dble(RWq(2,iq))  
        q3=rg(3,1,iop)*SQI(1,iqir)
     +    +rg(3,2,iop)*SQI(2,iqir)
     +    +rg(3,3,iop)*SQI(3,iqir)+dble(RWq(3,iq))  
        KGtmp(:,:)=0 
        call make_LG0(NTG,b1(1),b2(1),b3(1),Ecut_for_eps,
     +       Ecut_for_psi,q1,q2,q3,KGtmp(1,1),NG_for_eps,NG_for_psi) 
        NGQ_eps(iq)=NG_for_eps;NGQ_psi(iq)=NG_for_psi  
        write(6,'(i8,3f10.5,a8,i8,a8,i10)') 
     +  iq,q1,q2,q3,'NGeps',NG_for_eps,'NGpsi',NG_for_psi  
        LG0(:,:,iq)=-KGtmp(:,:) 
       endif 
      enddo 
!--
      Lq1=maxval(abs(LG0(1,:,:)))+1;write(6,*)'Lq1=',Lq1 
      Lq2=maxval(abs(LG0(2,:,:)))+1;write(6,*)'Lq2=',Lq2 
      Lq3=maxval(abs(LG0(3,:,:)))+1;write(6,*)'Lq3=',Lq3 
      allocate(packingq(-Lq1:Lq1,-Lq2:Lq2,-Lq3:Lq3,Nq_irr))   
      packingq(:,:,:,:)=0 
      do iq=1,Nq_irr 
       do ig=1,NGQ_eps(iq) 
        i1=LG0(1,ig,iq);j1=LG0(2,ig,iq);k1=LG0(3,ig,iq) 
        packingq(i1,j1,k1,iq)=ig 
        write(6,'(a,5i10)')'iq,ig,i1,j1,k1',iq,ig,i1,j1,k1
       enddo 
      enddo 
!
      write(6,*) 
      do iq=1,NTQ 
       write(6,'(i5,3f15.10,6i5)')iq,(SQ(i1,iq),i1=1,3),
     + numirrq(iq),numrotq(iq),trsq(iq),(RWq(j1,iq),j1=1,3)
      enddo  
!-- 
      NTGQ=maxval(NGQ_eps(:))!TOTAL NUMBER OF G VECTORS FOR EPSILON
      write(6,*)'NTGQ=',NTGQ  
!--
!OPEN(600-,R,FILE='DATA_EPSILON_fRPA',FORM='unformatted') 
      allocate(epstmp(NTGQ,NTGQ,ne));epstmp(:,:,:)=0.0D0!real8
      allocate(epstmpgm(NTGQ,NTGQ,ne,3));epstmpgm(:,:,:,:)=0.0D0!real8
      allocate(epsirr(NTGQ,NTGQ,ne,Nq_irr));epsirr(:,:,:,:)=0.0D0!real4 
!--
      ierr=CHDIR("./dir-eps") 
      call system('pwd') 
      do iq=1,Nq_irr 
       if(abs(SKI(1,iq))<1.0d-5.and.abs(SKI(2,iq))<1.0d-5.and.
     +    abs(SKI(3,iq))<1.0d-5)then 
        write(dirname,"('q',i3.3)")iq
        iqgm=iq 
        ierr=CHDIR(dirname) 
        do ix=1,3 
         file_num=600+(ix-1)  
         write(filename,'("dat.epsqw.",i3.3)')file_num 
         open(file_num,file=filename,form='unformatted') 
         write(6,'(a10,a15,a5,a10)') 
     +   'read: ',trim(filename),'in ',trim(dirname)  
         rewind(file_num) 
         NG_for_eps=NGQ_eps(iq)
         read(file_num)(((epstmpgm(ig,jg,ie,ix),ig=1,NG_for_eps)
     +                          ,jg=1,NG_for_eps),ie=1,ne)
        enddo!ix 
        epsirr(:,:,:,iq)!epsirr:real4,epstmp:real8 
     + =(epstmpgm(:,:,:,1)+epstmpgm(:,:,:,2)+epstmpgm(:,:,:,3))/3.0d0
       else
        write(dirname,"('q',i3.3)")iq
        ierr=CHDIR(dirname) 
        file_num=600 
        write(filename,'("dat.epsqw.",i3.3)')file_num 
        open(file_num,file=filename,form='unformatted') 
        write(6,'(a10,a15,a5,a10)') 
     +  'read: ',trim(filename),'in ',trim(dirname)  
        rewind(file_num) 
        NG_for_eps=NGQ_eps(iq)
        read(file_num)(((epstmp(ig,jg,ie),ig=1,NG_for_eps)
     +                         ,jg=1,NG_for_eps),ie=1,ne)
        epsirr(:,:,:,iq)=epstmp(:,:,:)!epsirr:real4,epstmp:real8 
       endif!q=0 or not 
       ierr=CHDIR("..") 
      enddo!iq   
      ierr=CHDIR("..") 
      call system('pwd') 
!--
!      do iq=1,Nq_irr 
!       NG_for_eps=NGQ_eps(iq)
!       write(6,*) NG_for_eps 
!       do ig=1,NG_for_eps 
!        write(6,*) epsirr(ig,ig,100,iq) 
!       enddo  
!      enddo  
!--
      WRITE(6,*) 
      WRITE(6,*)'====================================='
      WRITE(6,*)'=== INTERACTION CALCULATION START ==='
      WRITE(6,*)'====================================='
      WRITE(6,*) 
!--
!make rho(G,iw,q/=0)  
      allocate(rho(NTG,NWF,NTQ));rho(:,:,:)=0.0D0 
      do iq=1,NTQ 
       q1=SQ(1,iq);q2=SQ(2,iq);q3=SQ(3,iq)
       if(q1.ne.0.0d0.or.q2.ne.0.0d0.or.q3.ne.0.0d0) then 
        NG_for_eps=NGQ_eps(iq)
        NG_for_psi=NGQ_psi(iq)
        allocate(length_qg(NG_for_psi));length_qg(:)=0.0D0 
        do igL=1,NG_for_psi   
         igL1=LG0(1,igL,iq)
         igL2=LG0(2,igL,iq)
         igL3=LG0(3,igL,iq)
         qgL(:)=(q1+dble(igL1))*b1(:)
     +         +(q2+dble(igL2))*b2(:)
     +         +(q3+dble(igL3))*b3(:)
         qgL2=qgL(1)**2+qgL(2)**2+qgL(3)**2
         qgL1=dsqrt(qgL2) 
         length_qg(igL)=qgL1
        enddo 
!$OMP PARALLEL PRIVATE(ik,iw,igL,shift_G,ikq,C0_K,C0_KQ,
!$OMP&  fftwk,wfunc,rho_tmp,prho)
        allocate(prho(NTG,NWF));prho(:,:)=0.0D0 
        allocate(fftwk(Nl123*2));fftwk(:)=0.0d0  
        allocate(wfunc(Nl123*2));wfunc(:)=0.0d0 
        allocate(rho_tmp(NG_for_psi));rho_tmp(:)=0.0D0 
        allocate(C0_K(NTG));C0_K(:)=0.0d0     
        allocate(C0_KQ(NTG));C0_KQ(:)=0.0d0  
!$OMP DO 
        do ik=1,NTK 
         do iw=1,NWF 
          shift_G(:)=0
          call search_kq(NTK,SK0(1,1),q1,q2,q3,ik,ikq,shift_G(1))
          C0_K(:)=C0(:,iw,ik) 
          C0_KQ(:)=C0(:,iw,ikq) 
          rho_tmp(:)=0.0D0 
          call calc_InterStateMatrix(NTK,NTG,NG0(1),KG0(1,1,1),
     +     C0_K(1),C0_KQ(1),ik,ikq,nwx2,nwy2,nwz2,
     +     nfft1,nfft2,Nl123,wfunc(1),fftwk(1),fs,
     +     LG0(1,1,iq),NG_for_psi,shift_G(1),rho_tmp(1))
          do igL=1,NG_for_psi  
           prho(igL,iw)=prho(igL,iw)+rho_tmp(igL)
          enddo!igL
         enddo!iw  
        enddo!ik 
!$OMP END DO 
!$OMP CRITICAL 
        rho(:,:,iq)=rho(:,:,iq)+prho(:,:) 
!$OMP END CRITICAL 
        deallocate(prho,fftwk,wfunc,rho_tmp,C0_K,C0_KQ) 
!$OMP END PARALLEL 
        do iw=1,NWF
         do igL=1,NG_for_psi 
          rho(igL,iw,iq)=rho(igL,iw,iq)/length_qg(igL)
         enddo!igL
        enddo!iw 
        deallocate(length_qg)
        write(6,'(a10,x,2i5)')'FINISH iq',iq
!--
!gen rho(G,iw,q=0); iq=NTQ 
       elseif(q1==0.0d0.and.q2==0.0d0.and.q3==0.0d0)then 
        write(6,'(a,x,3f10.5)')'q1,q2,q3',q1,q2,q3 
        NG_for_eps=NGQ_eps(iq);NG_for_psi=NGQ_psi(iq)
        allocate(length_qg(NG_for_psi));length_qg(:)=0.0d0 
        do igL=1,NG_for_psi  
         igL1=LG0(1,igL,iq)
         igL2=LG0(2,igL,iq)
         igL3=LG0(3,igL,iq)
         if(igL1==0.and.igL2==0.and.igL3==0)then 
          NoG0=igL  
         else
          qgL(:)=(q1+dble(igL1))*b1(:)
     +          +(q2+dble(igL2))*b2(:)
     +          +(q3+dble(igL3))*b3(:)
          qgL2=qgL(1)**2+qgL(2)**2+qgL(3)**2
          qgL1=dsqrt(qgL2) 
          length_qg(igL)=qgL1
         endif  
        enddo 
!$OMP PARALLEL PRIVATE(ik,iw,igL,shift_G,ikq,C0_K,C0_KQ,
!$OMP&  fftwk,wfunc,rho_tmp,prho)
        allocate(prho(NTG,NWF));prho(:,:)=0.0D0 
        allocate(fftwk(Nl123*2));fftwk(:)=0.0d0  
        allocate(wfunc(Nl123*2));wfunc(:)=0.0d0 
        allocate(rho_tmp(NG_for_psi));rho_tmp(:)=0.0D0 
        allocate(C0_K(NTG));C0_K(:)=0.0d0     
        allocate(C0_KQ(NTG));C0_KQ(:)=0.0d0  
!$OMP DO 
        do ik=1,NTK 
         do iw=1,NWF 
          ikq=ik 
          shift_G(:)=0 
          C0_K(:)=C0(:,iw,ik) 
          C0_KQ(:)=C0(:,iw,ikq) 
          rho_tmp(:)=0.0D0 
          call calc_InterStateMatrix(NTK,NTG,NG0(1),KG0(1,1,1),
     +     C0_K(1),C0_KQ(1),ik,ikq,nwx2,nwy2,nwz2,
     +     nfft1,nfft2,Nl123,wfunc(1),fftwk(1),fs,
     +     LG0(1,1,iq),NG_for_psi,shift_G(1),rho_tmp(1))
          rho_tmp(NoG0)=0.0d0!rho=0.0 is correct 
          do igL=1,NG_for_psi  
           prho(igL,iw)=prho(igL,iw)+rho_tmp(igL)
          enddo!igL
         enddo!iw 
        enddo!ik  
!$OMP END DO 
!$OMP CRITICAL 
        rho(:,:,iq)=rho(:,:,iq)+prho(:,:) 
!$OMP END CRITICAL 
        deallocate(prho,fftwk,wfunc,rho_tmp,C0_K,C0_KQ) 
!$OMP END PARALLEL 
        do iw=1,NWF 
         do igL=1,NG_for_psi  
          if(igL/=NoG0)then 
           rho(igL,iw,iq)=rho(igL,iw,iq)/length_qg(igL)
          endif 
         enddo!igL 
        enddo!iw
        deallocate(length_qg)
        write(6,'(a10,x,2i5)')'FINISH iq',iq
       endif 
      enddo!iq 
      rho(:,:,:)=rho(:,:,:)/dble(NTK)
!--
!     do iq=1,NTQ
!      do iw=1,NWF 
!       do ig=1,NGQ_psi(iq) 
!        write(6,'(i5,2f15.10)') ig,rho(ig,iw,iq) 
!       enddo 
!      enddo 
!     enddo 
!     stop 
!--
!CORRECTION TERM TO DIVERGENCE 
!      ie=calc_ifreq          
!      write(6,*)'ie=',ie 
!--
       write(6,*)'HybertsenLouie'!Louie 
       write(6,*)'NoG0=',NoG0
       write(6,*)'iqgm=',iqgm 
       qsz=(6.0D0*(pi**2)/dble(NTQ)/dble(VOLUME))**(1.0D0/3.0D0)  
       write(6,*)'qsz=',qsz
       write(6,*)'2*qsz/pi=',(2.0d0/pi)*qsz
       chead=(2.0D0/pi)*qsz 
       write(6,*)'chead(HL) in eV=',chead*au    
       allocate(cheadw(ne));cheadw(:)=0.0d0 
       do iw=1,ne 
        cheadw(iw)=(2.0D0/pi)*epsirr(NoG0,NoG0,iw,iqgm)*qsz 
        write(6,*)'cheadw(HL) in eV=',cheadw(iw)*au   
        write(6,*)'1/eps=',epsirr(NoG0,NoG0,iw,Nq_irr)
       enddo 
!      endif 
!--
!wing 
!      allocate(cwing(NWF,NWF));cwing(:,:)=0.0D0 
!      allocate(wL(NWF));wL(:)=0.0D0 
!      allocate(wR(NWF));wR(:)=0.0D0 
!      do iw=1,NWF 
!       SUM_CMPX_1=0.0D0 
!       SUM_CMPX_2=0.0D0 
!       do igL=1,NGQ_eps(Nq_irr) 
!        if(igL/=NoG0) then 
!         SUM_CMPX_1
!     +  =SUM_CMPX_1
!     +  +rho(igL,iw,Nq_irr) 
!     +  *epsirr(igL,NoG0,ie,Nq_irr) 
!         SUM_CMPX_2
!     +  =SUM_CMPX_2
!     +  +CONJG(rho(igL,iw,Nq_irr))
!     +  *epsirr(NoG0,igL,ie,Nq_irr) 
!        Endif 
!       enddo 
!       wL(iw)=SUM_CMPX_1
!       wR(iw)=SUM_CMPX_2 
!      enddo 
!      do iw=1,NWF
!       write(6,'(4f15.10)') wL(iw),wR(iw)
!      enddo 
!      do iw=1,NWF
!       do jw=1,NWF 
!        cwing(iw,jw)=(qsz**2)*(wL(iw)+wR(jw))/pi 
!        write(6,'(a20,2f15.10)') 'cwing(eV)=',cwing(iw,jw)*au     
!       enddo 
!      enddo 
!      deallocate(cwing,wL,wR) 
!-- 
!make func(\mu,\nu,q)  
!--
       allocate(func(NWF,NWF,NTQ));func(:,:,:)=0.0D0
       do iq=1,NTQ
        do iw=1,NWF
         do jw=1,NWF
          SUM_CMPX=0.0D0
          do igL=1,NGQ_psi(iq)
           SUM_CMPX
     +    =SUM_CMPX
     +    +rho(igL,iw,iq)
     +    *CONJG(rho(igL,jw,iq))
          enddo 
          func(iw,jw,iq)=SUM_CMPX 
         enddo 
        enddo 
        write(6,'(a,x,i5)')'FINISH iq',iq
       enddo!iq 
       allocate(funcw(NWF,NWF,ne,NTQ));funcw(:,:,:,:)=0.0D0
       do iq=1,NTQ
        allocate(epsmk(NTGQ,NTGQ,ne))
        allocate(packtmp(-Lq1:Lq1,-Lq2:Lq2,-Lq3:Lq3)) 
        epsmk(:,:,:)=0.0d0 
        iqir=numirrq(iq)
        iop=numrotq(iq)
        packtmp(:,:,:)=packingq(:,:,:,iqir) 
        call make_eps(NTG,NTGQ,ne,trsq(iq),NGQ_eps(iq),LG0(1,1,iq), 
     +  RWq(1,iq),rginv(1,1,iop),pg(1,iop),nnp,Lq1,Lq2,Lq3,
     +  packtmp(-Lq1,-Lq2,-Lq3),epsirr(1,1,1,iqir),epsmk(1,1,1)) 
!--
!        if(iq.eq.xxx)then 
!         write(6,'(a5,3f15.10)')'SQ',SQ(1,iq),SQ(2,iq),SQ(3,iq)
!         write(6,'(a5,i10)')'iqir',iqir 
!         write(6,'(a5,i10)')'trsq',trsq(iq) 
!         write(6,'(a5,3i10)')'RWq',RWq(1,iq),RWq(2,iq),RWq(3,iq)
!         write(6,'(a5)')'epsirr' 
!         do ie=1,ne 
!          do igL=1,NGQ_eps(iqir)
!           do jgL=1,NGQ_eps(iqir) 
!            write(6,'(i5,i5,2f15.8)')igL,jgL,epsirr(igL,jgL,ie,iqir)
!           enddo
!          enddo 
!         enddo 
!         write(6,'(a5)')'epsmk' 
!         do ie=1,ne 
!          do igL=1,NGQ_eps(iq)
!           do jgL=1,NGQ_eps(iq) 
!            write(6,'(i5,i5,2f15.8)')igL,jgL,epsmk(igL,jgL,ie)
!           enddo
!          enddo 
!         enddo 
!        endif 
!$OMP PARALLEL PRIVATE(ie,iw,jw,SUM_CMPX,igL,jgL)
!$OMP DO 
        do ie=1,ne 
         do iw=1,NWF
          do jw=1,NWF
           SUM_CMPX=0.0D0
           do igL=1,NGQ_eps(iq)
            do jgL=1,NGQ_eps(iq) 
              SUM_CMPX
     +       =SUM_CMPX
     +       +rho(igL,iw,iq)
     +       *epsmk(igL,jgL,ie)
     +       *CONJG(rho(jgL,jw,iq))
            enddo 
           enddo 
           do igL=NGQ_eps(iq)+1,NGQ_psi(iq)
             SUM_CMPX
     +      =SUM_CMPX
     +      +rho(igL,iw,iq)
     +      *CONJG(rho(igL,jw,iq))
           enddo 
           funcw(iw,jw,ie,iq)=SUM_CMPX 
          enddo!jw 
         enddo!iw
        enddo!ie 
!$OMP END DO 
!$OMP END PARALLEL 
        deallocate(epsmk,packtmp)  
        write(6,'(a,x,i5)')'FINISH iq',iq
       enddo!iq 
!--
!write func  
!      write(6,*) 
!      write(6,*)'----'
!      write(6,*)'func'
!      write(6,*)'----'
!      write(6,*) 
!      Do iq=1,NTQ
!       write(6,*) SQ(1,iq),SQ(2,iq),SQ(3,iq) 
!       Do iw=1,NWF
!        write(6,'(100F12.5)')(func(iw,jw,iq),jw=1,NWF)
!       ENDDO      
!       write(6,*) 
!      ENDDO      
!write funcw  
!      write(6,*) 
!      write(6,*)'-----'
!      write(6,*)'funcw'
!      write(6,*)'-----'
!      write(6,*) 
!      do iq=1,NTQ
!       write(6,*) SQ(1,iq),SQ(2,iq),SQ(3,iq) 
!       do iw=1,NWF
!        write(6,'(100F12.5)')(funcw(iw,jw,1,iq),jw=1,NWF)
!       enddo       
!       write(6,*) 
!      enddo 
!-- 
!make MATRIX_W  
       allocate(V_MAT(NWF,NWF,-Na1:Na1,-Na2:Na2,-Na3:Na3))  
       V_MAT(:,:,:,:,:)=0.0D0 
       do ia1=-Na1,Na1
        do ia2=-Na2,Na2
         do ia3=-Na3,Na3
          do iw=1,NWF
           do jw=1,NWF 
            SUM_CMPX=0.0D0 
            do iq=1,NTQ 
             q1=SQ(1,iq);q2=SQ(2,iq);q3=SQ(3,iq)
             phase=tpi*(q1*dble(ia1)+q2*dble(ia2)+q3*dble(ia3)) 
             pf=exp(-ci*phase) 
             SUM_CMPX=SUM_CMPX+func(iw,jw,iq)*pf
            enddo 
            V_MAT(iw,jw,ia1,ia2,ia3)
     +     =2.0D0*tpi*SUM_CMPX/VOLUME/DBLE(NTQ)+chead 
           enddo 
          enddo 
         enddo 
        enddo 
       enddo 
       allocate(W_MAT(NWF,NWF,-Na1:Na1,-Na2:Na2,-Na3:Na3,ne))  
       W_MAT(:,:,:,:,:,:)=0.0D0 
       do ie=1,ne
        do ia1=-Na1,Na1
         do ia2=-Na2,Na2
          do ia3=-Na3,Na3
           do iw=1,NWF
            do jw=1,NWF 
             SUM_CMPX=0.0D0 
             do iq=1,NTQ 
              q1=SQ(1,iq);q2=SQ(2,iq);q3=SQ(3,iq)
              phase=tpi*(q1*dble(ia1)+q2*dble(ia2)+q3*dble(ia3)) 
              pf=exp(-ci*phase) 
              SUM_CMPX=SUM_CMPX+funcw(iw,jw,ie,iq)*pf
             enddo 
             W_MAT(iw,jw,ia1,ia2,ia3,ie)
     +      =2.0D0*tpi*SUM_CMPX/VOLUME/DBLE(NTQ)+cheadw(ie) 
            enddo 
           enddo 
          enddo 
         enddo 
        enddo 
       enddo 
!--
!WRITE V_MAT  
       write(6,*)'# '
       write(6,*)'#-----'
       write(6,*)'#V_MAT'
       write(6,*)'#-----'
       write(6,*)'# '
       do ia1=-Na1,Na1 
        do ia2=-Na2,Na2 
         do ia3=-Na3,Na3 
          write(6,*) ia1,ia2,ia3           
          do iw=1,NWF
           write(6,'(100F12.5)')(DBLE(V_MAT(iw,jw,ia1,ia2,ia3)*au), 
     +                                         jw=1,NWF)
          enddo       
          write(6,*) 
         enddo 
        enddo 
       enddo 
!--
!OPEN(200,W,FILE='dat.Vmat')
       call system('rm -rf dir-intW') 
       call system('mkdir dir-intW') 
       ierr=CHDIR("./dir-intW") 
       call system("pwd") 
       OPEN(200,FILE='dat.Vmat') 
       rewind(200)
       write(200,'(a)')'#Bare Coulomb: V' 
       write(200,'(a)')'#1:R1, 2:R2, 3:R3 (lattice vector)'
       write(200,'(a)')'#1:i, 2:j, 3:Re(V_ij) [eV], 4:Im(V_ij) [eV]' 
       do ia1=-Na1,Na1
        do ia2=-Na2,Na2
         do ia3=-Na3,Na3
          write(200,*) ia1,ia2,ia3           
          do iw=1,NWF 
           do jw=1,NWF 
            write(200,'(i5,i5,2f20.10)')
     +      iw,jw,V_MAT(iw,jw,ia1,ia2,ia3)*au
           enddo!jw 
          enddo!iw 
          write(200,*)
         enddo 
        enddo 
       enddo 
       close(200) 
!--
!WRITE W_MAT  
!--
       ie=calc_ifreq          
       write(6,*)'# '
       write(6,*)'#-------------------------'
       write(6,*)'# W_MAT(omega=calc_ifreq) '
       write(6,*)'#-------------------------'
       write(6,*)'# '
       do ia1=-Na1,Na1 
        do ia2=-Na2,Na2 
         do ia3=-Na3,Na3 
          write(6,*) ia1,ia2,ia3           
          do iw=1,NWF
           write(6,'(100F12.5)')(DBLE(W_MAT(iw,jw,ia1,ia2,ia3,ie)*au), 
     +                                         jw=1,NWF)
          enddo       
          write(6,*) 
         enddo 
        enddo 
       enddo 
!--
!OPEN(201,W,FILE='./dir-intW/dat.Wmat')
       ie=calc_ifreq          
       OPEN(201,FILE='dat.Wmat') 
       rewind(201)
       write(201,'(a)')'#Screened Coulomb: W at omega=calc_ifreq' 
       write(201,'(a)')'#1:R1, 2:R2, 3:R3 (lattice vector)'
       write(201,'(a)')'#1:i, 2:j, 3:Re(W_ij) [eV], 4:Im(W_ij) [eV]' 
       do ia1=-Na1,Na1
        do ia2=-Na2,Na2
         do ia3=-Na3,Na3
          write(201,*) ia1,ia2,ia3           
          do iw=1,NWF 
           do jw=1,NWF
            write(201,'(i5,i5,2f20.10)')
     +      iw,jw,W_MAT(iw,jw,ia1,ia2,ia3,ie)*au 
           enddo!jw 
          enddo!iw 
          write(201,*)
         enddo 
        enddo 
       enddo 
       close(201)  
!--
!OPEN(1000,W,FILE='./dir-intW/dat.VvsR.xxx')
       do iw=1,NWF 
        write(filename,'("dat.VvsR.",i3.3)') iw 
        file_num=1000+(iw-1) 
        OPEN(file_num,FILE=filename) 
        rewind(file_num) 
        write(file_num,'(a)')'#R dependence of bare Coulomb' 
        write(file_num,'(a)')
     +  '#1:R [AA], 2:Re[V(R)] [eV], 3:Im[V(R)] [eV]'
        do ia1=-Na1,Na1
         do ia2=-Na2,Na2
          do ia3=-Na3,Na3
           do jw=1,NWF 
            rij(:)=dble(ia1)*a1(:)+dble(ia2)*a2(:)+dble(ia3)*a3(:)
     +            +(coord(:,jw)-coord(:,iw)) 
            rij2=rij(1)**2+rij(2)**2+rij(3)**2
            rij1=dsqrt(rij2)
            write(file_num,'(3F12.5)')
     +      rij1*bohr,V_MAT(iw,jw,ia1,ia2,ia3)*au  
           enddo!jw 
          enddo!ia3 
         enddo!ia2 
        enddo!ia1 
        close(file_num) 
       enddo!iw 
!--
!OPEN(2000,W,FILE='./dir-intW/dat.WvsR')
       ie=calc_ifreq          
       do iw=1,NWF 
        write(filename,'("dat.WvsR.",i3.3)') iw 
        file_num=2000+(iw-1) 
        OPEN(file_num,FILE=filename) 
        rewind(file_num) 
        write(file_num,'(a)')
     +  '#R dependence of screened Coulomb at omega=calc_ifreq' 
        write(file_num,'(a)') 
     +  '#1:R [AA], 2:Re[W(R)] [eV], 3:Im[W(R)] [eV]'  
        do ia1=-Na1,Na1
         do ia2=-Na2,Na2
          do ia3=-Na3,Na3
            do jw=1,NWF 
             rij(:)=dble(ia1)*a1(:)+dble(ia2)*a2(:)+dble(ia3)*a3(:)
     +             +(coord(:,jw)-coord(:,iw)) 
             rij2=rij(1)**2+rij(2)**2+rij(3)**2
             rij1=dsqrt(rij2)
             write(file_num,'(3F12.5)')
     +       rij1*bohr,W_MAT(iw,jw,ia1,ia2,ia3,ie)*au    
            enddo!jw 
          enddo!ia3 
         enddo!ia2 
        enddo!ia1  
        close(file_num) 
       enddo!iw  
!--
!OPEN(3000,W,FILE='./dir-intW/dat.UvsE')
      do iw=1,NWF
       do jw=iw,NWF 
        write(filename,'("dat.UvsE.",i3.3,"-",i3.3)') iw,jw  
        file_num=3000+(iw-1)*NWF+(jw-1)  
        OPEN(file_num,FILE=filename) 
        rewind(file_num) 
        write(file_num,'(a)')
     +  '#Omega dependence of onsite Coulomb' 
        write(file_num,'(a)') 
     +  '#1:Re(w) [eV], 2:Im(w) [eV], 3:Re[U(w)] [eV], 4:Im[U(w)] [eV]'
        do ie=1,ne
         write(file_num,'(4F12.5)') em(ie)*au,W_MAT(iw,jw,0,0,0,ie)*au 
        enddo 
        close(file_num) 
       enddo
      enddo
      ierr=CHDIR("..") 
      call system("pwd") 
!--
!write model parameters 
      ierr=CHDIR("./dir-model") 
      call system("pwd") 
      ie=calc_ifreq          
      !wcut_model=wcut_model/au 
      !
      !old: 20180604
      !
      ! call wrt_model(nkb1,nkb2,nkb3,NTK,Na1,Na2,Na3,NWF,wcut_model,
      !+     W_MAT(1,1,-Na1,-Na2,-Na3,ie)) 
      !
      !new: 20181127 
      !
      !WEIGHT_R 20170317 by Y.Nomura 
      !
      allocate(WEIGHT_R(-Na1:Na1,-Na2:Na2,-Na3:Na3)); WEIGHT_R=1.0d0
      SUM_REAL=0.0d0 
      do ia1=-Na1,Na1
       do ia2=-Na2,Na2
        do ia3=-Na3,Na3
         if(abs(ia1)==Na1.and.mod(nkb1,2)==0.and.Na1/=0) 
     +    WEIGHT_R(ia1,ia2,ia3)=WEIGHT_R(ia1,ia2,ia3)*0.5d0 
         if(abs(ia2)==Na2.and.mod(nkb2,2)==0.and.Na2/=0) 
     +    WEIGHT_R(ia1,ia2,ia3)=WEIGHT_R(ia1,ia2,ia3)*0.5d0 
         if(abs(ia3)==Na3.and.mod(nkb3,2)==0.and.Na3/=0) 
     +    WEIGHT_R(ia1,ia2,ia3)=WEIGHT_R(ia1,ia2,ia3)*0.5d0 
         SUM_REAL=SUM_REAL+WEIGHT_R(ia1,ia2,ia3)
        enddo!ia3
       enddo!ia2
      enddo!ia1 
      write(6,'(a20,f15.8,i8)')'SUM_WEIGHT,NTK',SUM_REAL,NTK  
      if(abs(SUM_REAL-dble(NTK))>1.0d-6)then 
       stop 'SUM_WEIGHT/=NTK'
      endif 
      !
      call wrt_model(Na1,Na2,Na3,NWF,W_MAT(1,1,-Na1,-Na2,-Na3,ie),
     +     WEIGHT_R(-Na1,-Na2,-Na3)) 
      !
      ierr=CHDIR("..") 
      call system("pwd") 
      !
      deallocate(C0,epsirr) 
      !
      STOP
      END           
