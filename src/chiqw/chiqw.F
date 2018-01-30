      PROGRAM DIELECTRIC 
      use mpi
      use fft_3d 
      use m_tetrahedron
      use m_rdinput       
      include "config.h"
!mpi global 
      comm_glb=MPI_COMM_WORLD
      call MPI_INIT(ierr)
      call MPI_COMM_RANK(comm_glb,myrank_glb,ierr)
      call MPI_COMM_SIZE(comm_glb,nproc_glb,ierr)
      write(6,*)'myrank_glb=',myrank_glb,'nproc_glb=',nproc_glb
!read input from master 
      if(myrank_glb.eq.0)then 
       call read_input(nproc_glb)!20170926  
      endif 
!--
      call MPI_Bcast(N_CALC_BAND,1,MPI_INTEGER,0,comm_glb,ierr) 
      call MPI_Bcast(Ecut_for_eps,1,MPI_DOUBLE_PRECISION,0,
     +               comm_glb,ierr)
      call MPI_Bcast(shift_ef,1,MPI_DOUBLE_PRECISION,0,comm_glb,ierr)
      call MPI_Bcast(Max_excitation_energy,1,MPI_DOUBLE_PRECISION,0,
     +               comm_glb,ierr) 
      call MPI_Bcast(delta_ex,1,MPI_DOUBLE_PRECISION,0,comm_glb,ierr)
      call MPI_Bcast(Green_func_delt,1,MPI_DOUBLE_PRECISION,0,
     +               comm_glb,ierr) 
      call MPI_Bcast(delt,1,MPI_DOUBLE_PRECISION,0,comm_glb,ierr) 
      call MPI_Bcast(ttrhdrn_dmna,1,MPI_DOUBLE_PRECISION,0,
     +               comm_glb,ierr) 
      call MPI_Bcast(dmna,1,MPI_DOUBLE_PRECISION,0,comm_glb,ierr) 
      call MPI_Bcast(ttrhdrn_dmnr,1,MPI_DOUBLE_PRECISION,0,
     +               comm_glb,ierr) 
      call MPI_Bcast(dmnr,1,MPI_DOUBLE_PRECISION,0,comm_glb,ierr) 
      call MPI_Bcast(MPI_io_rank,1,MPI_INTEGER,0,comm_glb,ierr) 
      call MPI_Bcast(io_rank,1,MPI_INTEGER,0,comm_glb,ierr) 
      call MPI_Bcast(MPI_num_proc_per_qcomm,1,MPI_INTEGER,0,
     +               comm_glb,ierr) 
      call MPI_Bcast(num_of_mpi_per_comm,1,MPI_INTEGER,0,comm_glb,ierr)
      call MPI_Bcast(MPI_num_qcomm,1,MPI_INTEGER,0,comm_glb,ierr)
      call MPI_Bcast(num_of_comm,1,MPI_INTEGER,0,comm_glb,ierr)
      call MPI_Bcast(Num_freq_grid,1,MPI_INTEGER,0,comm_glb,ierr)
      call MPI_Bcast(nen,1,MPI_INTEGER,0,comm_glb,ierr) 
      call MPI_Bcast(Lower_bound_energy_window,1,MPI_DOUBLE_PRECISION,0,
     +               comm_glb,ierr) 
      call MPI_Bcast(E_LOWER,1,MPI_DOUBLE_PRECISION,0,comm_glb,ierr) 
      call MPI_Bcast(Upper_bound_energy_window,1,MPI_DOUBLE_PRECISION,0,
     +               comm_glb,ierr) 
      call MPI_Bcast(E_UPPER,1,MPI_DOUBLE_PRECISION,0,comm_glb,ierr) 
      call MPI_Bcast(flg_cRPA,1,MPI_INTEGER,0,comm_glb,ierr) 
      call MPI_Bcast(flg_cRPA_band,1,MPI_INTEGER,0,comm_glb,ierr) 
      call MPI_Bcast(flg_cRPA_ewin,1,MPI_INTEGER,0,comm_glb,ierr) 
      call MPI_Bcast(flg_calc_type,1,MPI_INTEGER,0,comm_glb,ierr) 
      call MPI_Bcast(file_num_log_start,1,MPI_INTEGER,0,comm_glb,ierr)
      call MPI_Bcast(file_num_chi_start,1,MPI_INTEGER,0,comm_glb,ierr)
      call MPI_Bcast(file_num_eps_start,1,MPI_INTEGER,0,comm_glb,ierr) 
      call MPI_Bcast(file_num_chi_base_start,1,MPI_INTEGER,0,
     +               comm_glb,ierr) 
      call MPI_Bcast(file_num_eps_base_start,1,MPI_INTEGER,0,
     +               comm_glb,ierr) 
      call MPI_Bcast(n_calc_q,1,MPI_INTEGER,0,comm_glb,ierr) 
      if(myrank_glb /= 0) allocate(calc_num_k(n_calc_q))
      call MPI_Bcast(calc_num_k,n_calc_q,MPI_INTEGER,0,comm_glb,ierr)
!--
!      call MPI_FINALIZE(ierr)
!      stop
!--
      shift_ef=shift_ef/au 
      delta_ex=delta_ex/au 
      delt=delt/au 
      id_sub=myrank_glb/num_of_mpi_per_comm 
!mpi split 
      call MPI_BARRIER(comm_glb,ierr)
      call MPI_COMM_SPLIT(comm_glb,id_sub,myrank_glb,comm,ierr)
      call MPI_COMM_RANK(comm,myrank,ierr)
      call MPI_COMM_SIZE(comm,nproc,ierr)
!--
      call MPI_BARRIER(comm_glb,ierr)
      start_time=MPI_Wtime() 
      if(myrank_glb.eq.MPI_io_rank)then 
       WRITE(6,*) 
       WRITE(6,*)'============='
       WRITE(6,*)' setup start '
       WRITE(6,*)'============='
       WRITE(6,*) 
      endif 
!--
!OPEN(117,R,FILE='dat.bandcalc') 
      OPEN(117,FILE='./dir-wfn/dat.bandcalc') 
      rewind(117) 
      read(117,*) Ecut_for_psi 
      read(117,*) FermiEnergy  
      read(117,*) Etot
      if (myrank_glb.eq.0) then 
       write(6,'(a20,f15.10)')'Ecut_for_psi=',Ecut_for_psi 
       write(6,'(a20,f15.10)')'FermiEnergy (eV)=',FermiEnergy*au  
       write(6,'(a20,f15.10)')'Etot (au)=',Etot
      endif 
!--
!artificially fermi level shifted <--- Option, basically shift=0.00 
      FermiEnergy=FermiEnergy+shift_ef 
      if(myrank_glb.eq.MPI_io_rank)then 
       write(6,*) 
       write(6,*)'========================='
       write(6,*)'FermiEnergy=',FermiEnergy
       write(6,*)'========================='
       write(6,*) 
      endif 
!--
!cRPA 
      flg_cRPA_band=0
      flg_cRPA_ewin=0
      if(flg_cRPA==1)then 
       flg_cRPA_band=1!band-cRPA
!      flg_cRPA_ewin=1!window-cRPA 
      endif 
!--
      if(flg_cRPA_ewin==1)then 
       if(E_LOWER==E_UPPER)then 
        write(6,*)'ERROR: wrong setting energy window; stop' 
        stop
       endif 
       if(E_LOWER==0.0d0)then 
        write(6,*)'WARNING: E_LOWER=0.0eV'
       endif 
       if(E_UPPER==0.0d0)then 
        write(6,*)'WARNING: E_UPPER=0.0eV'
       endif 
       E_LOWER=E_LOWER/au 
       E_UPPER=E_UPPER/au 
!20180110 
!      else 
!       E_LOWER=FermiEnergy-0.5d0/au 
!       E_UPPER=FermiEnergy+0.5d0/au 
!       E_LOWER=FermiEnergy-1.5d0/au 
!       E_UPPER=FermiEnergy+1.5d0/au 
!       E_LOWER=E_LOWER/au!20171220 
!       E_UPPER=E_UPPER/au!20171220 
!--
      endif 
      if(myrank_glb.eq.MPI_io_rank)then 
       write(6,'(a20,i5)')'flg_cRPA     =',flg_cRPA
       write(6,'(a20,i5)')'flg_cRPA_band=',flg_cRPA_band
       write(6,'(a20,i5)')'flg_cRPA_ewin=',flg_cRPA_ewin 
       write(6,'(a20,f15.10)')'EPS_UPPER (eV)=',E_UPPER*au 
       write(6,'(a20,f15.10)')'EPS_LOWER (eV)=',E_LOWER*au 
      endif 
!--
!OPEN(105,R,FILE='dat.lattice') 
      OPEN(105,FILE='./dir-wfn/dat.lattice') 
      REWIND(105)
      READ(105,*) a1(1),a1(2),a1(3)!a1 vector
      READ(105,*) a2(1),a2(2),a2(3)!a2 vector
      READ(105,*) a3(1),a3(2),a3(3)!a3 vector
      CLOSE(105)
      call OUTER_PRODUCT(a2(1),a3(1),b1(1))
      VOLUME=a1(1)*b1(1)+a1(2)*b1(2)+a1(3)*b1(3)
      b1(:)=b1(:)*tpi/VOLUME 
      call OUTER_PRODUCT(a3(1),a1(1),b2(1))
      b2(:)=b2(:)*tpi/VOLUME 
      call OUTER_PRODUCT(a1(1),a2(1),b3(1))
      b3(:)=b3(:)*tpi/VOLUME 
      if(myrank_glb.eq.MPI_io_rank)then 
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
      endif 
!---
      if(myrank_glb.eq.MPI_io_rank)then 
       write(6,*) 
       write(6,*)'======================'
       write(6,*)'CHECK OF (A_VEC*B_VEC)'
       write(6,*)'======================'
       write(6,*) 
       write(6,*)'a1*b1=',b1(1)*a1(1)+b1(2)*a1(2)+b1(3)*a1(3)
       write(6,*)'a2*b2=',b2(1)*a2(1)+b2(2)*a2(2)+b2(3)*a2(3)
       write(6,*)'a3*b3=',b3(1)*a3(1)+b3(2)*a3(2)+b3(3)*a3(3)
      endif 
!--
!OPEN(100,R,FILE='dat.symmetry') 
      OPEN(100,FILE='./dir-wfn/dat.symmetry') 
      rewind(100) 
      read(100,*) nsymq 
      read(100,*) nnp 
      allocate(rg(3,3,nsymq));rg(:,:,:)=0 
      allocate(pg(3,nsymq));pg(:,:)=0
      allocate(rginv(3,3,nsymq));rginv(:,:,:)=0.0d0 
      do iop=1,nsymq
       read(100,*)((rg(i,j,iop),i=1,3),j=1,3) 
       read(100,*)(pg(i,iop),i=1,3)   
      enddo 
      close(100) 
      rginv=rg 
      do iop=1,nsymq
       call invmat(3,rginv(1,1,iop)) 
      enddo 
      if(myrank_glb.eq.0)then!master
       write(6,*) 
       write(6,*)'========'
       write(6,*)'SYMMETRY'
       write(6,*)'========'
       write(6,*) 
       write(6,*)'nsymq=',nsymq 
       write(6,*)'nnp=',nnp 
       do iop=1,nsymq  
        write(6,'(a5,i4)')'sym=',iop
        write(6,'(a10)')'rg and pg' 
        do i=1,3
         write(6,'(3i5,5x,i5)')(rg(i,j,iop),j=1,3),pg(i,iop) 
        enddo
        write(6,*) 
       enddo
       rginv=rg 
       do iop=1,nsymq
        call invmat(3,rginv(1,1,iop)) 
       enddo 
       do iop=1,nsymq  
        write(6,'(a5,i4)')'sym=',iop
        write(6,'(a10)')'rginv' 
        do i=1,3
         write(6,'(3f10.5)')(rginv(i,j,iop),j=1,3) 
        enddo
        write(6,*) 
       enddo
       write(6,*)'finish read dat.symmetry'
      endif 
!check 
!do iop=1,nsymq
!write(6,*) iop
!do i=1,3
!write(6,'(3I5,1x,3F15.10)') (rg(i,j,iop),j=1,3),(rginv(i,j,iop),j=1,3)
!enddo 
!enddo 
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
       read(101,*)(SKI(i,ik),i=1,3) 
      enddo 
      close(101) 
      if(myrank_glb.eq.MPI_io_rank)then 
       write(6,*)'finish fort.101' 
      endif 
      call est_NTK(Nk_irr,SKI(1,1),NTK,nkb1,nkb2,nkb3)  
!--
!20170415 
      allocate(SKI_list(3,Nk_irr));SKI_list(:,:)=0.0D0 
      do ik=1,Nk_irr
       call search_list(SKI(1,ik),SKI_list(1,ik)) 
      enddo  
!--
      if(myrank_glb.eq.MPI_io_rank)then 
       write(6,'(a24,4i10)')'nkb1,nkb2,nkb3,NTK=',nkb1,nkb2,nkb3,NTK  
      endif 
!---
!OPEN(132,R,FILE='dat.nkm')20170406 
      OPEN(132,FILE='./dir-wfn/dat.nkm')!20170406 
      allocate(NGI(Nk_irr));NGI(:)=0
      rewind(132)
      do ik=1,Nk_irr 
       read(132,*) NGI(ik) 
      enddo 
      close(132) 
      NTG=maxval(abs(NGI(:))) 
      if(myrank_glb.eq.MPI_io_rank)then 
       write(6,'(a10,i10)')'NTG=',NTG  
      endif 
!---
!OPEN(111,R,FILE='dat.eigenvalue') 
      OPEN(111,FILE='./dir-wfn/dat.eigenvalue') 
      rewind(111)
      read(111,*) NTB!20170406  
      allocate(E_EIGI(NTB,Nk_irr));E_EIGI(:,:)=0.0d0
      do ik=1,Nk_irr 
       do ib=1,NTB 
        read(111,*) E_EIGI(ib,ik) 
       enddo!ib 
      enddo!ik          
      close(111) 
!--
      if(N_CALC_BAND==0) N_CALC_BAND=NTB  
      if(myrank_glb.eq.MPI_io_rank)then 
       write(6,'(a15,i10)')'NTB=',NTB  
       write(6,'(a15,i10)')'N_CALC_BAND=',N_CALC_BAND 
      endif 
!      call MPI_FINALIZE(ierr)
!      stop
!--
!OPEN(131,R,FILE='dat.MAT_NLr_rNL')---not active---
      allocate(MAT_NLr_rNL(3,NTB,NTB,Nk_irr));MAT_NLr_rNL=0.0d0   
!rewind(131) 
!close(131) 
!---
!OPEN(104,R,FILE='dat.kg') 
      OPEN(104,FILE='./dir-wfn/dat.kg') 
      rewind(104) 
      allocate(KGI(3,NTG,Nk_irr));KGI=0 
      do ik=1,Nk_irr 
       read(104,*) NG_for_psi 
       if(NG_for_psi/=NGI(ik))THEN 
        write(6,*)'ERROR;STOP;NG_for_psi should be NGI(ik)'   
        write(6,*)'NG_for_pis=',NG_for_psi,'NGI(ik)=',NGI(ik)
        stop 
       endif 
       do ig=1,NG_for_psi 
        read(104,*)(KGI(i,ig,ik),i=1,3) 
       enddo 
      enddo!ik 
      close(104)
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
      enddo!ik 
!--
      Ecut_for_psi=maxval(LKGI(:,:))+1.0d-8!20170512 
      if(Ecut_for_eps==0.0d0)then!default
       Ecut_for_eps=Ecut_for_psi/10.0d0!20171010  
      endif 
!--
      if(myrank_glb.eq.MPI_io_rank)then 
       do ik=1,Nk_irr 
        write(6,*) maxval(LKGI(:,ik))
       enddo 
       write(6,*) 
       write(6,*) maxval(LKGI(:,:))
       write(6,*)'reset Ecut_for_psi=',Ecut_for_psi 
       write(6,*)'set Ecut_for_eps=',Ecut_for_eps 
      endif 
!--
      L1=maxval(abs(KGI(1,:,:)))+1
      L2=maxval(abs(KGI(2,:,:)))+1
      L3=maxval(abs(KGI(3,:,:)))+1
      if(myrank_glb.eq.MPI_io_rank)then 
       write(6,*)'maxKGI(1)=',L1 
       write(6,*)'maxKGI(2)=',L2 
       write(6,*)'maxKGI(3)=',L3 
      endif 
      allocate(packing(-L1:L1,-L2:L2,-L3:L3,Nk_irr));packing(:,:,:,:)=0 
      do ik=1,Nk_irr 
       do ig=1,NGI(ik) 
        i1=KGI(1,ig,ik);j1=KGI(2,ig,ik);k1=KGI(3,ig,ik) 
        packing(i1,j1,k1,ik)=ig 
       enddo 
      enddo 
!--
!fft grid
      m1=maxval(abs(KGI(1,:,:))) 
      m2=maxval(abs(KGI(2,:,:))) 
      m3=maxval(abs(KGI(3,:,:))) 
      nwx2=algn235(2*m1) 
      nwy2=algn235(2*m2) 
      nwz2=algn235(2*m3) 
      nfft1=nwx2+1
      nfft2=nwy2+1
      nfft3=nwz2+1
      Nl123=nfft1*nfft2*nfft3 
      if(myrank_glb.eq.MPI_io_rank)then 
       write(6,*)'nwx2=',nwx2 
       write(6,*)'nwy2=',nwy2 
       write(6,*)'nwz2=',nwz2 
       write(6,*)'nfft1=',nfft1 
       write(6,*)'nfft2=',nfft2 
       write(6,*)'nfft3=',nfft3 
       write(6,*)'NL123=',Nl123  
      endif 
      call fft3_init(nwx2,nwy2,nwz2,nfft1,nfft2,nfft3,fs) 
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
!       numirr(ik)=ik; numrot(ik)=1; trs(ik)=1; RW(1:3,ik)=0
!      enddo 
!      jk=Nk_irr 
!--
      jk=0 
      do ik=1,Nk_irr
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
       SK0(:,jk)=ktmp(:)
       numirr(jk)=ik;numrot(jk)=iop;trs(jk)=1;RW(:,jk)=RWtmp(:)
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
      if (myrank_glb.eq.0) then 
       write(6,*) 
       write(6,*)'====='
       write(6,*)' SK0 '
       write(6,*)'====='
       write(6,*) 
       DO ik=1,NTK
        WRITE(6,'(I5,3F15.10)') ik,SK0(:,ik) 
       ENDDO 
      endif 
!--
!gen(NB0,NG0,KG0,E_EIG)
      allocate(NG0(NTK));NG0(:)=0
      allocate(NB0(NTK));NB0(:)=0
      allocate(KG0(3,NTG,NTK));KG0(:,:,:)=0 
      allocate(KGtmp(3,NTG));KGtmp(:,:)=0 
!20161207 
!      do jk=1,Nk_irr 
!       NG0(jk)=NGI(jk);NB0(jk)=NTB;KG0(:,:,jk)=KGI(:,:,jk) 
!      enddo 
!      do jk=Nk_irr+1,NTK 
!--
      do jk=1,NTK 
       if(trs(jk)==1)then 
        ik=numirr(jk);iop=numrot(jk) 
        ktmp(1)=rg(1,1,iop)*SKI(1,ik)+rg(1,2,iop)*SKI(2,ik)
     +         +rg(1,3,iop)*SKI(3,ik)+dble(RW(1,jk)) 
        ktmp(2)=rg(2,1,iop)*SKI(1,ik)+rg(2,2,iop)*SKI(2,ik)
     +         +rg(2,3,iop)*SKI(3,ik)+dble(RW(2,jk))  
        ktmp(3)=rg(3,1,iop)*SKI(1,ik)+rg(3,2,iop)*SKI(2,ik)
     +         +rg(3,3,iop)*SKI(3,ik)+dble(RW(3,jk))  
        NB0(jk)=NTB 
        call make_KG0(NTG,b1(1),b2(1),b3(1),Ecut_for_psi,
     +       ktmp(1),ktmp(2),ktmp(3),KG0(1,1,jk),NG_for_psi)
        if(NG_for_psi/=NGI(ik))then 
         write(6,*)'ERROR;STOP;NG_for_psi should be NGI(ik)'   
         write(6,*)'NG_for_psi=',NG_for_psi,'NGI(ik)=',NGI(ik);STOP
        endif 
        NG0(jk)=NG_for_psi  
       elseif(trs(jk)==-1)then  
        ik=numirr(jk);iop=numrot(jk) 
        ktmp(1)=rg(1,1,iop)*SKI(1,ik)+rg(1,2,iop)*SKI(2,ik)
     +         +rg(1,3,iop)*SKI(3,ik)+dble(RW(1,jk)) 
        ktmp(2)=rg(2,1,iop)*SKI(1,ik)+rg(2,2,iop)*SKI(2,ik)
     +         +rg(2,3,iop)*SKI(3,ik)+dble(RW(2,jk))  
        ktmp(3)=rg(3,1,iop)*SKI(1,ik)+rg(3,2,iop)*SKI(2,ik)
     +         +rg(3,3,iop)*SKI(3,ik)+dble(RW(3,jk))  
        NB0(jk)=NTB 
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
      deallocate(KGtmp) 
!--
      if(myrank_glb.eq.MPI_io_rank)then 
       write(6,*) 
       write(6,*)'==========='
       write(6,*)'NB0 and NG0'
       write(6,*)'==========='
       write(6,*) 
       do ik=1,NTK 
        write(6,'(2I20)') NB0(ik),NG0(ik) 
       enddo 
      endif 
!20170418 
!mkdir dir-eps/qxx
      if(myrank_glb.eq.0)then!master 
       !inquire(directory='./dir-eps',exist=dir_e) 
       !if(dir_e)then
       ! write(6,*) 'dir exist'
       !else
       ! call system('rm -rf dir-eps') 
       ! call system('mkdir dir-eps') 
       !endif 
       call system('mkdir -p dir-eps')!20171208
       ierr=CHDIR("./dir-eps") 
       call system('pwd') 
       do ik=1,Nk_irr 
        !write(dirname,"('q',i3.3)")ik 
        !inquire(directory=dirname,exist=dir_e) 
        !if(dir_e)then
        ! write(6,*) 'dir exist'
        !else
        ! write(command,"('mkdir q',i3.3)")ik 
        ! call system(command) 
        !endif 
        write(command,"('mkdir -p q',i3.3)")ik 
        call system(command) 
       enddo  
       ierr=CHDIR("..") 
       call system('pwd') 
      endif!myrank_glb==0 
      flush(6) 
!Nq_irr, calc_qlist
      allocate(calc_qlist(3,Nk_irr));calc_qlist(:,:)=0.0D0 
      if(myrank_glb.eq.0)then!master 
       ierr=CHDIR("./dir-eps") 
       call system('pwd') 
!all-q-calc 
       if(flg_calc_type==0)then!all-q-calc
        iq=0 
        do ik=1,Nk_irr 
         write(dirname,"('q',i3.3)")ik 
         ierr=CHDIR(dirname) 
         inquire(file='dat.log.400',exist=file_e) 
         if(file_e)then
          write(6,*)'dat.log.400 exists in',trim(dirname)  
         else
          iq=iq+1
          calc_qlist(:,iq)=SKI(:,ik) 
         endif!file_e 
         ierr=CHDIR("..") 
        enddo  
        Nq_irr=iq  
       endif!flg_calc_type==0
!gamma-only 
       if(flg_calc_type==1)then!gamma-only 
        iq=0;ix=0  
        do ik=1,Nk_irr 
         if(abs(SKI(1,ik))<1.0d-5.and.abs(SKI(2,ik))<1.0d-5.and.
     +      abs(SKI(3,ik))<1.0d-5)then 
          ix=ik 
         endif 
        enddo 
        write(dirname,"('q',i3.3)")ix  
        ierr=CHDIR(dirname) 
        inquire(file='dat.log.400',exist=file_e) 
        if(file_e)then
         write(6,*)'dat.log.400 exists in',trim(dirname)  
        else
         iq=iq+1
         calc_qlist(:,iq)=SKI(:,ix) 
        endif!file_e 
        ierr=CHDIR("..") 
        Nq_irr=iq  
       endif!flg_calc_type==1
!manual mode 
       if(flg_calc_type==2)then!manual mode 
        iq=0 
        do ix=1,N_calc_q 
         ik=calc_num_k(ix) 
         write(dirname,"('q',i3.3)")ik 
         ierr=CHDIR(dirname) 
         inquire(file='dat.log.400',exist=file_e) 
         if(file_e)then
          write(6,*)'dat.log.400 exists in',trim(dirname)  
         else
          iq=iq+1
          calc_qlist(:,iq)=SKI(:,ik) 
         endif!file_e 
         ierr=CHDIR("..") 
        enddo!ix 
        Nq_irr=iq  
       endif!flg_calc_type==2 
!--
       ierr=CHDIR("..") 
       call system('pwd') 
      endif!myrank_glb==0 
      flush(6) 
      call MPI_Bcast(Nq_irr,1,MPI_INTEGER,0,comm_glb,ierr) 
      call MPI_Bcast(calc_qlist,3*Nk_irr,MPI_DOUBLE_PRECISION,
     +               0,comm_glb,ierr) 
!--
!      call MPI_FINALIZE(ierr)
!      stop
!--
!sample-q
!--
      pnq=Nq_irr/num_of_comm
      if(mod(Nq_irr,num_of_comm).ne.0)then
       nbufq=pnq+1
      else
       nbufq=pnq
      endif
      allocate(SQI(3,nbufq)) 
      allocate(SQ(3,nbufq)) 
      allocate(LG0(3,NTG,nbufq)) 
      allocate(NGQ(nbufq)) 
!--
      if(id_sub.lt.mod(Nq_irr,num_of_comm))then
       pnq=pnq+1
       bnq=pnq*id_sub+1 
       enq=bnq+pnq-1
      else
       bnq=(pnq+1)*mod(Nq_irr,num_of_comm) 
     +    +(pnq)*(id_sub-mod(Nq_irr,num_of_comm))+1 
       enq=bnq+pnq-1
      endif
!--
      do iq=1,Nq_irr
       if((bnq.le.iq).and.(iq.le.enq))then 
        SQI(:,iq-bnq+1)=calc_qlist(:,iq) 
       endif
      enddo 
!--
!OPEN(127,R,FILE='SAMPLE_Qirr')
!      OPEN(127,FILE='./dir-wfn/dat.sample-q')
!      SQI(:,:)=0.0d0 
!      read(127,*) idum 
!      do iq=1,Nq_irr
!       if((bnq.le.iq).and.(iq.le.enq))then 
!        read(127,*)(SQI(i,iq-bnq+1),i=1,3) 
!       else 
!        read(127,*) 
!       endif
!      enddo 
!--
!OPEN(101,R,FILE='dat.sample-k'='dat.sample-q')
!       OPEN(101,FILE='./dir-wfn/dat.sample-k')
!       SQI(:,:)=0.0d0 
!       rewind(101) 
!       read(101,*) idum 
!       do iq=1,Nq_irr
!        if((bnq.le.iq).and.(iq.le.enq))then 
!         read(101,*)(SQI(i,iq-bnq+1),i=1,3) 
!        else 
!         read(101,*) 
!        endif
!       enddo 
!--
!      elseif(flg_calc_type==1)then!gamma only
!       Nq_irr=1
!       pnq=Nq_irr/num_of_comm
!       if(mod(Nq_irr,num_of_comm).ne.0)then
!        nbufq=pnq+1
!       else
!        nbufq=pnq
!       endif
!       allocate(SQI(3,nbufq)) 
!       allocate(SQ(3,nbufq)) 
!       allocate(LG0(3,NTG,nbufq)) 
!       allocate(NGQ(nbufq)) 
!!--
!       if(id_sub.lt.mod(Nq_irr,num_of_comm))then
!        pnq=pnq+1
!        bnq=pnq*id_sub+1 
!        enq=bnq+pnq-1
!       else
!        bnq=(pnq+1)*mod(Nq_irr,num_of_comm) 
!     +     +(pnq)*(id_sub-mod(Nq_irr,num_of_comm))+1 
!        enq=bnq+pnq-1
!       endif
!!      SQI(:,:)=0.0d0 
!       do iq=1,Nq_irr
!       if((bnq.le.iq).and.(iq.le.enq))then 
!        SQI(:,iq-bnq+1)=0.0d0!gamma only 
!       endif
!      endif 
!--
!20161207 
!shift[-1/2:1/2]  
      do iq=1,pnq 
       call search_list(SQI(1,iq),SQ(1,iq)) 
      enddo  
      if(myrank.eq.0)then 
       write(6,*) 
       write(6,*)'======================'
       write(6,*)'calclated SQ (SQI->SQ)'
       write(6,*)'======================'
       write(6,*) 
       write(6,*)'Nq_irr=',Nq_irr 
       write(6,*) 
       do iq=1,pnq 
        write(6,'(a10,3f10.5,2x,3f10.5)')'SQI->SQ',SQI(:,iq),SQ(:,iq) 
       enddo 
      endif  
      flush(6) 
!--
!      call MPI_BARRIER(comm_glb,ierr)
!      call MPI_FINALIZE(ierr)
!      stop
!--
!      if(myrank.eq.0)then!master  
!       ierr=CHDIR("./dir-eps") 
!       call system('pwd') 
!       do iq=1,pnq
!        write(file_num_log,'(a8,3(1x,f10.5),1x,a,1x,i3)') 
!     +  '#SQ=',SQ(1,iq),SQ(2,iq),SQ(3,iq),'myrank_glb=',myrank_glb 
!       enddo 
!       ierr=CHDIR("..") 
!       call system('pwd') 
!      endif 
!      call MPI_FINALIZE(ierr)
!      stop
!--
      if(myrank.eq.0)then 
       write(6,*) 
       write(6,*)'========'
       write(6,*)'make LG0'
       write(6,*)'========'
       write(6,*) 
      endif 
      do iq=1,pnq 
       q1=SQ(1,iq)
       q2=SQ(2,iq)
       q3=SQ(3,iq)
       if(myrank.eq.0)then 
        write(6,'(i5,3f15.8)') iq,q1,q2,q3 
       endif 
       call make_KG0(NTG,b1(1),b2(1),b3(1),Ecut_for_eps,
     +               q1,q2,q3,LG0(1,1,iq),NG_for_eps) 
       if(myrank.eq.0)then 
        write(6,'(a,i10)')'NG_for_eps',NG_for_eps 
       endif 
       NGQ(iq)=NG_for_eps  
      enddo  
!--
      call MPI_BARRIER(comm_glb,ierr)
!--
      allocate(index_kpt(nkb1,nkb2,nkb3));index_kpt(:,:,:)=0     
      call make_index_kpt(NTK,nkb1,nkb2,nkb3,SK0(1,1),index_kpt(1,1,1)) 
!--
!      if(myrank_glb.eq.MPI_io_rank)then 
!       do iz=1,nkb3 
!        do iy=1,nkb2
!         do ix=1,nkb1
!          ik=index_kpt(ix,iy,iz)
!          write(6,'(i8,3f15.10)')ik,SK0(:,ik)
!         enddo
!        enddo
!       enddo
!      endif 
!--
      allocate(WindowInside(NTB));WindowInside(:)=0
      allocate(Tindx(NTB));Tindx(:)=0 
      allocate(E_AVE(NTB));E_AVE(:)=0.0d0 
      allocate(band_max(NTB));band_max(:)=0.0d0
      allocate(band_min(NTB));band_min(:)=0.0d0
      allocate(prob(NTB,NTK));prob(:,:)=0.0d0 
!--
      if(flg_cRPA==0)then 
       call judge_FermiInside(Nk_irr,NTB,E_EIGI(1,1),FermiEnergy,
     + WindowInside(1),Tindx(1),E_AVE(1),band_max(1),band_min(1))
      endif 
!--
      if(flg_cRPA==1.and.flg_cRPA_band==1)then 
       call judge_FermiInside(Nk_irr,NTB,E_EIGI(1,1),FermiEnergy,
     + WindowInside(1),Tindx(1),E_AVE(1),band_max(1),band_min(1))
      endif 
!--
      if(flg_cRPA==1.and.flg_cRPA_ewin==1)then 
       call judge_WindowInside(Nk_irr,NTB,E_EIGI(1,1),E_LOWER,E_UPPER,
     + WindowInside(1),Tindx(1),E_AVE(1),band_max(1),band_min(1)) 
      endif 
!--
      if (myrank_glb.eq.0) then 
       write(6,*) 
       write(6,*)'Winside'
       write(6,*) 
       do ib=1,NTB
        write(6,*) ib,WindowInside(ib)
       enddo 
!
       write(6,*) 
       write(6,*)'Tindx'
       write(6,*) 
       do ib=1,NTB
        write(6,*) ib,Tindx(ib)
       enddo 
!
       write(6,*) 
       write(6,*)'bandMax and bandMin in eV'
       write(6,*) 
       do ib=1,NTB
        write(6,'(i5,2f15.10)') ib,band_max(ib)*au,band_min(ib)*au 
       enddo 
!
       write(6,*) 
       write(6,*)'E_AVE in eV'
       write(6,*) 
       do ib=1,NTB
        write(6,*) ib,E_AVE(ib)*au 
       enddo 
      endif 
!--
!band-cRPA
      if(flg_cRPA_band==1)then 
!--
!OPEN(149,R,FILE='dat.ns-nb')
       OPEN(149,FILE='./dir-wan/dat.ns-nb') 
       rewind(149) 
       allocate(Ns(NTK));Ns=0
       allocate(Nb(NTK));Nb=0
       do ik=1,NTK 
        read(149,*) Ns(ik),Nb(ik) 
       enddo  
       Mb=maxval(Nb) 
       if(myrank_glb.eq.MPI_io_rank)then 
        write(6,*)'Mb=',Mb 
       endif 
!--
!OPEN(150,R,FILE='dat.umat')
       OPEN(150,FILE='./dir-wan/dat.umat') 
       rewind(150) 
       read(150,*) NWF 
       allocate(UNT(Mb,NWF,NTK));UNT(:,:,:)=0.0d0 
       do ik=1,NTK
        do jb=1,Nb(ik)
         read(150,*)(UNT(jb,jw,ik),jw=1,NWF) 
        enddo 
       enddo 
!--
!probability
       prob(:,:)=0.0d0 
       do ik=1,NTK 
        do ib=1,Nb(ik) 
         SUM_CMPX=0.0D0            
         do iw=1,NWF 
          SUM_CMPX=SUM_CMPX+CONJG(UNT(ib,iw,ik))*UNT(ib,iw,ik) 
         enddo!iw  
         prob(Ns(ik)+ib,ik)=SUM_CMPX 
        enddo!ib   
       enddo!ik        
!-- 
       if(myrank_glb.eq.0)then!master  
        WRITE(6,*) 
        WRITE(6,*)'========================'
        WRITE(6,*)'PROBABILITY MATRIX:<a|b>'
        WRITE(6,*)'========================'
        WRITE(6,*) 
        do ik=1,NTK 
         SUM_CMPX=0.0d0 
         do ib=1,NTB 
          if(abs(prob(ib,ik))>1.0d-4)then 
           write(6,'(i5,2f10.5)') ib,prob(ib,ik) 
          endif 
          SUM_CMPX=SUM_CMPX+prob(ib,ik) 
         enddo!ib 
         write(6,'(a10,2f15.10)')'total=',SUM_CMPX 
         write(6,*)       
        enddo!ik   
       endif  
      else 
       if(myrank_glb.eq.0)then!master  
        WRITE(6,*) 
        WRITE(6,*)'======================================='
        WRITE(6,*)'PROBABILITY MATRIX:<a|b> is set to zero'
        WRITE(6,*)'======================================='
        WRITE(6,*) 
       endif  
      endif!band-cRPA 
!--
!SETTING MPI-process PARAMETERS
      call MPI_BARRIER(comm_glb,ierr)
      if(myrank_glb.eq.MPI_io_rank)then 
       write(6,*) 
       write(6,*)'MPI started'
       write(6,*) 
      endif 
      id=0
      ip=0
      iu=0
      do i_band=1,N_CALC_BAND
       fs1=WindowInside(i_band)
       if(fs1==0) id=id+1
       if(fs1==1) ip=ip+1
       if(fs1==2) iu=iu+1
      enddo 
!
      N_docc=id
      N_pocc=ip
      N_uocc=iu 
      if(myrank_glb.eq.MPI_io_rank)then 
       write(6,*)'N_docc=',N_docc 
       write(6,*)'N_pocc=',N_pocc 
       write(6,*)'N_uocc=',N_uocc 
      endif 
!
      N_occ=N_docc+N_pocc
      N_vir=N_uocc+N_pocc  
      if(myrank_glb.eq.MPI_io_rank)then 
       write(6,*)'N_occ=',N_occ 
       write(6,*)'N_vir=',N_vir  
      endif 
!
      pnv=N_vir/nproc
      if(mod(N_vir,nproc).ne.0)then
       nbufv=pnv+1
      else
       nbufv=pnv
      endif
      if(myrank.lt.mod(N_vir,nproc))then
       pnv=pnv+1
       bnv=pnv*myrank+1 
       env=bnv+pnv-1
      else
       bnv=(pnv+1)*mod(N_vir,nproc) 
     +    +pnv*(myrank-mod(N_vir,nproc))+1 
       env=bnv+pnv-1
      endif
!--
!      allocate(datav(nbufv),newdatav(nbufv))
!      do i=1,pnv
!         datav(i)=bnv+i-1
!      end do
!--
      pno=N_occ/nproc
      if(mod(N_occ,nproc).ne.0)then
       nbufo=pno+1
      else
       nbufo=pno
      endif
      if(myrank.lt.mod(N_occ,nproc)) then
       pno=pno+1
       bno=pno*myrank+1 
       eno=bno+pno-1
      else
       bno=(pno+1)*mod(N_occ,nproc) 
     +    +pno*(myrank-mod(N_occ,nproc))+1 
       eno=bno+pno-1
      endif
!--
!      allocate(datao(nbufo),newdatao(nbufo))
!      do i=1,pno
!         datao(i)=bno+i-1
!      end do
!--
!      call MPI_BARRIER(comm,ierr)
!      if (myrank_glb.eq.0) then 
!       write(6,'(a,4(1x,i5),1x,a,100(1x,f10.5))')'myrank=',myrank, 
!     +      bnv,env,pnv,'datav',datav(1:pnv)
!      endif 
!--
!      call MPI_BARRIER(comm,ierr)
!      if (myrank_glb.eq.0) then 
!       write(6,'(a,4(1x,i5),1x,a,100(1x,f10.5))')'myrank=',myrank, 
!     +       bno,eno,pno,'datao',datao(1:pno)
!      endif 
!      write(6,*)'myrank=',myrank,bno,eno,pno,nbufo 
!      write(6,*)'myrank=',myrank,bnv,env,pnv,nbufv
!--
      call MPI_BARRIER(comm_glb,ierr)
      allocate(E_OCC(1:nbufo,Nk_irr));E_OCC(:,:)=0.0d0 
      allocate(E_VIR(1:nbufv,Nk_irr));E_VIR(:,:)=0.0d0 
      allocate(E_VIR_new(1:nbufv,Nk_irr));E_VIR_new(:,:)=0.0d0 
      allocate(W_OCC(1:nbufo));W_OCC(:)=0
      allocate(W_VIR(1:nbufv));W_VIR(:)=0
      allocate(W_VIR_new(1:nbufv));W_VIR_new(:)=0
      allocate(T_OCC(1:nbufo));T_OCC(:)=0
      allocate(T_VIR(1:nbufv));T_VIR(:)=0
      allocate(T_VIR_new(1:nbufv));T_VIR_new(:)=0
      allocate(EO_AVE(1:nbufo));EO_AVE(:)=0.0d0 
      allocate(EV_AVE(1:nbufv));EV_AVE(:)=0.0d0 
      allocate(EV_AVE_new(1:nbufv));EV_AVE_new(:)=0.0d0 
!--
      allocate(P_OCC(1:nbufo,NTK));P_OCC(:,:)=0.0d0 
      allocate(P_VIR(1:nbufv,NTK));P_VIR(:,:)=0.0d0 
      allocate(P_VIR_new(1:nbufv,NTK));P_VIR_new(:,:)=0.0d0 
!--
      do i_band=1,pno!N_occ
       W_OCC(i_band)=WindowInside(bno-1+i_band)
       E_OCC(i_band,:)=E_EIGI(bno-1+i_band,:)
       T_OCC(i_band)=Tindx(bno-1+i_band) 
       EO_AVE(i_band)=E_AVE(bno-1+i_band) 
       P_OCC(i_band,:)=prob(bno-1+i_band,:)
      enddo 
      do i_band=1,pnv!N_vir
       W_VIR(i_band)=WindowInside(N_docc+bnv-1+i_band)
       E_VIR(i_band,:)=E_EIGI(N_docc+bnv-1+i_band,:)
       T_VIR(i_band)=Tindx(N_docc+bnv-1+i_band) 
       EV_AVE(i_band)=E_AVE(N_docc+bnv-1+i_band) 
       P_VIR(i_band,:)=prob(N_docc+bnv-1+i_band,:)
      enddo 
!--
      call MPI_BARRIER(comm_glb,ierr)
      if(myrank_glb.eq.MPI_io_rank)then 
       write(6,*) 
       write(6,*)'W_OCC'
       write(6,*) 
       do i_band=1,pno!N_occ
        write(6,*) i_band,W_OCC(i_band)
       enddo 
       write(6,*) 
       write(6,*)'W_VIR'
       write(6,*) 
       do i_band=1,pnv!N_vir
        write(6,*) i_band,W_VIR(i_band)
       enddo 
       write(6,*) 
       write(6,*)'T_OCC'
       write(6,*) 
       do i_band=1,pno!N_occ
        write(6,*) i_band,T_OCC(i_band)
       enddo 
       write(6,*) 
       write(6,*)'T_VIR'
       write(6,*) 
       do i_band=1,pnv!N_vir
        write(6,*) i_band,T_VIR(i_band)
       enddo 
       write(6,*) 
       write(6,*)'EO_AVE'
       write(6,*) 
       do i_band=1,pno!N_occ
        write(6,*) i_band,EO_AVE(i_band)
       enddo 
       write(6,*) 
       write(6,*)'EV_AVE'
       write(6,*) 
       do i_band=1,pnv!N_vir
        write(6,*) i_band,EV_AVE(i_band)
       enddo 
      endif 
      flush(6) 
!-- 
      call MPI_BARRIER(comm_glb,ierr)
      allocate(OCC(NTG,1:nbufo,Nk_irr));OCC(:,:,:)=0.0d0 
      allocate(VIR(NTG,1:nbufv,Nk_irr));VIR(:,:,:)=0.0d0 
      allocate(VIR_new(NTG,1:nbufv,Nk_irr));VIR_new(:,:,:)=0.0d0 
      call MPI_BARRIER(comm_glb,ierr)
!--
!OPEN(102,FILE='./dir-wfn/dat.wfn',FORM='unformatted') 
      OPEN(102,FILE='./dir-wfn/dat.wfn',FORM='unformatted') 
      rewind(102) 
      read(102) ncomp!20170512  
      if(ncomp/=1)then 
       write(6,*)'This program not suport ncomp/=1; then stop'
       stop
      endif 
      call system('pwd') 
      do ik=1,Nk_irr  
       if(myrank/=0)then!not master  
        dest=0
        pnocld=pno 
        call MPI_Send(pnocld,1,MPI_INTEGER,dest,0,comm,ierr) 
        call MPI_Recv(OCC(1,1,ik),nbufo*NTG,MPI_DOUBLE_COMPLEX,
     +                dest,0,comm,status,ierr) 
        pnvcld=pnv 
        call MPI_Send(pnvcld,1,MPI_INTEGER,dest,0,comm,ierr) 
        call MPI_Recv(VIR(1,1,ik),nbufv*NTG,MPI_DOUBLE_COMPLEX,
     +                dest,0,comm,status,ierr) 
       else 
        do ib=1,pno
         read(102)(OCC(ig,ib,ik),ig=1,NGI(ik))
        enddo 
        allocate(Otmp(NTG,nbufo));Otmp(:,:)=0.0d0  
        do source=1,nproc-1
         write(6,*)'source=',source 
         call MPI_Recv(pnocld,1,MPI_INTEGER,source,0,comm,status,ierr) 
         do ib=1,pnocld 
          read(102)(Otmp(ig,ib),ig=1,NGI(ik))
         enddo 
         call MPI_Send(Otmp(1,1),nbufo*NTG,MPI_DOUBLE_COMPLEX,
     +                 source,0,comm,ierr) 
        enddo!source 
        deallocate(Otmp) 
        do ib=1,N_pocc
         backspace(102)
        enddo 
        do ib=1,pnv 
         read(102)(VIR(ig,ib,ik),ig=1,NGI(ik))
        enddo 
        allocate(Vtmp(NTG,nbufv));Vtmp(:,:)=0.0d0  
        do source=1,nproc-1
         call MPI_Recv(pnvcld,1,MPI_INTEGER,source,0,comm,status,ierr) 
         do ib=1,pnvcld 
          read(102)(Vtmp(ig,ib),ig=1,NGI(ik))
         enddo 
         call MPI_Send(Vtmp(1,1),nbufv*NTG,MPI_DOUBLE_COMPLEX,
     +                 source,0,comm,ierr) 
        enddo!source 
        deallocate(Vtmp) 
        do ib=1,NTB-N_CALC_BAND 
        read(102)
        enddo 
       endif  
      enddo 
!MPI_Broadcast
      IF(num_of_comm.gt.1)then 
       id_sub_perp=MOD(myrank_glb,num_of_mpi_per_comm)  
       write(6,*)'glb=',myrank_glb,'id_sub_perp=',id_sub_perp 
       call MPI_COMM_SPLIT(comm_glb,id_sub_perp,myrank_glb
     +                    ,comm_perp,ierr)
       call MPI_COMM_RANK(comm_perp,myrank_perp,ierr)
       call MPI_COMM_SIZE(comm_perp,nproc_perp,ierr)
       call MPI_Bcast(OCC,nbufo*Nk_irr*NTG,MPI_DOUBLE_COMPLEX,  
     +                0,comm_perp,ierr) 
       call MPI_Bcast(VIR,nbufv*Nk_irr*NTG,MPI_DOUBLE_COMPLEX,  
     +                0,comm_perp,ierr) 
      ENDIF 
!--
!      call MPI_BARRIER(comm_glb,ierr)
!      end_time=MPI_Wtime()
!      diff_time=end_time-start_time 
!      if(myrank.eq.0) then 
!       write(file_num_log,*)'#TOTAL TIME=',diff_time 
!      endif 
!      call MPI_FINALIZE(ierr)
!      stop
!--
      emax=dabs(maxval(band_max(1:N_CALC_BAND))
     +         -minval(band_min(1:N_CALC_BAND))) 
      if(myrank_glb.eq.MPI_io_rank)then 
       write(6,*)'emax(eV)=',emax*au 
       write(6,*)'bandmax(eV)=',maxval(band_max(1:N_CALC_BAND))*au 
       write(6,*)'bandmin(eV)=',minval(band_min(1:N_CALC_BAND))*au 
      endif 
      allocate(em(nen));em(:)=0.0d0 
      call gen_grid(nen,emax,em(1)) 
      if(myrank_glb.eq.0)then!master  
       write(6,*) 
       write(6,*)'==============='
       write(6,*)' freqency grid ' 
       write(6,*)'==============='
       write(6,*) 
       do ie=1,nen
        write(6,'(A5,2F10.5)')'em=',em(ie) 
       enddo 
       write(6,*) 
       ierr=CHDIR("./dir-eps") 
       call system('pwd') 
!--
!OPEN(300,W,FILE='dat.chi_cutoff')!20170402  
       OPEN(300,FILE='./dat.chi_cutoff')!20170402 
       write(300,*) Ecut_for_eps!20170402  
!--
!OPEN(135,R,FILE='dat.wgrid') 
       OPEN(135,FILE='./dat.wgrid') 
       write(135,*) nen!20170331 nen=Num_freq_grid 
       do ie=1,nen 
        write(135,'(2F15.10)') em(ie)
       enddo 
       ierr=CHDIR("..") 
       call system('pwd') 
      endif 
!     call MPI_FINALIZE(ierr)
!     stop
!--
      call MPI_BARRIER(comm,ierr)
      if(myrank_glb.eq.MPI_io_rank)then 
       WRITE(6,*) 
       WRITE(6,*)'====================='
       WRITE(6,*)'=== dielmat start ==='
       WRITE(6,*)'====================='
       WRITE(6,*) 
      endif 
!--
      allocate(imt1(4*nkb1*nkb2*nkb3*6));imt1(:)=0    
!--
!     call ttrhdrn_mkidx(nkb1,nkb2,nkb3,imt1(1)) 
!--
!20170325 for new ttrhdrn
      call ttrhdrn_mkidx(nkb1,nkb2,nkb3,imt1(1),b1(1),b2(1),b3(1))
!--
!this mpi calculates eps(q,w) for bnq to enq points 
!      if(myrank.eq.0)then!master  
!       ierr=CHDIR("./dir-eps") 
!       call system('pwd') 
!       write(file_num_eps) 
!     + bnq+num_start_irrq-1,enq+num_start_irrq-1,pnq  
!       write(6,*)'this mpi calculates eps(q,w) for bnq to enq points',
!     +            bnq+num_start_irrq-1,enq+num_start_irrq-1,pnq  
!       ierr=CHDIR("..") 
!       call system('pwd') 
!      endif 
!--
      do iq=1,pnq 
       q1=SQ(1,iq);q2=SQ(2,iq);q3=SQ(3,iq)
       NG_for_eps=NGQ(iq) 
       if(myrank.eq.0)then!master 
        write(6,'(i5,3f15.8)') iq,q1,q2,q3 
       endif 
!q.ne.0 
      if((q1/=0.0d0).or.(q2/=0.0d0).or.(q3/=0.0d0))then 
       allocate(chiqw(NG_for_eps,NG_for_eps,nen))
       allocate(pchiqw(NG_for_eps,NG_for_eps,nen))
       allocate(i_for_pair(N_CALC_BAND*N_CALC_BAND))
       allocate(j_for_pair(N_CALC_BAND*N_CALC_BAND))
       chiqw(:,:,:)=0.0D0 
       pchiqw(:,:,:)=0.0D0 
!wavefunction rotation loop  
       do irot=0,nproc-1
!double loop-->single loop 
       i_for_pair(:)=0
       j_for_pair(:)=0
       ij=0
       do jb=1,pno 
        do ib=1,pnv 
         fs1=W_VIR(ib)
         fs2=W_OCC(jb)
!constrained excitation 
         bci=EV_AVE(ib) 
         bcj=EO_AVE(jb) 
         if(abs(bci-bcj)>=delta_ex)cycle 
!--
!         if(myrank.eq.0)then 
!          write(6,*) abs(bci-bcj)*au,delta_ex*au  
!         endif 
!--
!1:window-cRPA not default 
        if(fs1==1.and.fs2==1.and.flg_cRPA_ewin==1)cycle 
!--
         ij=ij+1
         i_for_pair(ij)=ib  
         j_for_pair(ij)=jb  
         if(myrank.eq.0)then 
          write(6,*)'ij=',ij,ib,jb 
         endif 
        enddo 
       enddo 
       N_PAIR=ij 
       if(myrank.eq.0)write(6,*)'FOR q/=0,N_PAIR=',N_PAIR 
!$OMP PARALLEL PRIVATE(ij,i_band,j_band
!$OMP&                ,e1_1D,e2_1D,ik,shift_G,ikq
!$OMP&                ,C0_K,C0_KQ,dsgn,kpq
!$OMP&                ,wfunc,fftwk,m_tmp,ism
!$OMP&                ,e1_3D,e2_3D,ikb3,ikb2,ikb1
!$OMP&                ,xow,xow_1D,ca1,ie,ig,jg
!$OMP&                ,pc0d,iomp,ikir,iop)  
       allocate(pc0d(NG_for_eps,NG_for_eps,nen));pc0d=(0.d0,0.d0)
       allocate(fftwk(Nl123*2)) 
       allocate(wfunc(Nl123*2)) 
       allocate(m_tmp(NG_for_eps))
       allocate(ism(NTK,NG_for_eps))
       allocate(C0_K(NTG))    
       allocate(C0_KQ(NTG))    
       allocate(kpq(NTK)) 
       allocate(e1_1D(NTK))
       allocate(e2_1D(NTK))
       allocate(e1_3D(nkb1,nkb2,nkb3)) 
       allocate(e2_3D(nkb1,nkb2,nkb3))
       allocate(xow(nen,nkb1,nkb2,nkb3)) 
       allocate(xow_1D(NTK,nen)) 
       allocate(ca1(4*nen))  
!$OMP DO
       do ij=1,N_PAIR 
        i_band=i_for_pair(ij) 
        j_band=j_for_pair(ij) 
        do ik=1,NTK 
         shift_G(:)=0 
         call search_kq(NTK,SK0(1,1),q1,q2,q3,ik,ikq,shift_G(1))  
         kpq(ik)=ikq
!C0_K=OCC(jb,ik)
         ikir=numirr(ik)
         iop=numrot(ik) 
         C0_K(:)=0.0d0 
         call make_C0(NTG,trs(ik),NG0(ik),KG0(1,1,ik),RW(1,ik),
     +   rginv(1,1,iop),pg(1,iop),nnp,L1,L2,L3,
     +   packing(-L1,-L2,-L3,ikir),OCC(1,j_band,ikir),C0_K(1)) 
!C0_KQ=VIR(ib,ikq)
         ikir=numirr(ikq)
         iop=numrot(ikq) 
         C0_KQ(:)=0.0d0 
         call make_C0(NTG,trs(ikq),NG0(ikq),KG0(1,1,ikq),RW(1,ikq),     
     +   rginv(1,1,iop),pg(1,iop),nnp,L1,L2,L3,
     +   packing(-L1,-L2,-L3,ikir),VIR(1,i_band,ikir),C0_KQ(1)) 
!calc ISM 
         call calc_InterStateMatrix(NTK,NTG,NG0(1),KG0(1,1,1),
     +                              C0_K(1),C0_KQ(1),
     +                              ik,ikq,nwx2,nwy2,nwz2,
     +                              nfft1,nfft2,Nl123,
     +                              wfunc(1),fftwk(1),fs,LG0(1,1,iq),
     +                              NG_for_eps,shift_G(1),m_tmp(1))
         ism(ik,:)=m_tmp(:) 
        enddo!ik 
!Using ttrhdrn_chi0 
        e1_1D(:)=0.0d0
        e2_1D(:)=0.0d0 
        do ik=1,NTK
         ikir=numirr(ik) 
         e1_1D(ik)=cmplx(E_OCC(j_band,ikir),0.0d0)  
         ikq=kpq(ik)
         ikir=numirr(ikq) 
         e2_1D(ik)=cmplx(E_VIR(i_band,ikir),0.0d0)  
        enddo 
        e1_3D(:,:,:)=0.0d0 
        e2_3D(:,:,:)=0.0d0
        do ikb3=1,nkb3
         do ikb2=1,nkb2
          do ikb1=1,nkb1
           ik=index_kpt(ikb1,ikb2,ikb3)
           e1_3D(ikb1,ikb2,ikb3)=e1_1D(ik)
           e2_3D(ikb1,ikb2,ikb3)=e2_1D(ik)
          enddo 
         enddo 
        enddo 
!freq dep 
        xow(:,:,:,:)=0.0d0 
        ca1(:)=0.0d0
        dsgn=1.0d0
        call ttrhdrn_chi0(dmna,dmnr,nkb1,nkb2,nkb3,imt1(1),
     +                   e1_3D(1,1,1),e2_3D(1,1,1),
     +                   FermiEnergy,FermiEnergy,delt,dsgn,
     +                   nen,em(1),ca1(1),xow(1,1,1,1)) 
        dsgn=-1.0d0
        ca1(:)=0.0d0
        call ttrhdrn_chi0(dmna,dmnr,nkb1,nkb2,nkb3,imt1(1),
     +                   e1_3D(1,1,1),e2_3D(1,1,1),
     +                   FermiEnergy,FermiEnergy,delt,dsgn,
     +                   nen,em(1),ca1(1),xow(1,1,1,1)) 
!--
        xow_1D(:,:)=0.0d0 
        do ikb3=1,nkb3
         do ikb2=1,nkb2
          do ikb1=1,nkb1
           ik=index_kpt(ikb1,ikb2,ikb3)
           do ie=1,nen
            xow_1D(ik,ie)=xow(ie,ikb1,ikb2,ikb3)
           enddo 
          enddo 
         enddo 
        enddo 
!--
        do ie=1,nen
         do ig=1,NG_for_eps
          do jg=1,NG_for_eps
           do ik=1,NTK 
            ikq=kpq(ik)
            pc0d(ig,jg,ie)
     +     =pc0d(ig,jg,ie)
     +     +xow_1D(ik,ie)
     +     *CONJG(ism(ik,ig))*ism(ik,jg)/dble(NTK)
     +     *(1.0d0-P_OCC(j_band,ik)*P_VIR(i_band,ikq)) 
           enddo 
          enddo 
         enddo 
        enddo 
!--
        iomp=omp_get_thread_num() 
        if(myrank.eq.0.and.iomp.eq.0)then 
         write(6,*)'#',i_band,j_band 
        endif 
!--
       enddo!ij
!$OMP END DO
!$OMP CRITICAL
       pchiqw=pchiqw+pc0d
!$OMP END CRITICAL
       deallocate(pc0d,fftwk,wfunc,m_tmp,ism)
       deallocate(C0_K,C0_KQ,kpq,e1_1D,e2_1D,e1_3D,e2_3D,xow,xow_1D,ca1) 
!$OMP END PARALLEL
!---
       if(myrank.eq.0)write(6,*)'#MPI start' 
       call MPI_BARRIER(comm,ierr)
       dest=modulo(myrank+1,nproc)
       source=modulo(myrank-1,nproc)
       tmpnv(1)=pnv
       tmpnv(2)=bnv
       tmpnv(3)=env
       call MPI_SENDRECV(   tmpnv,      3, MPI_INTEGER, 
     +                       dest,      0, 
     +                   newtmpnv,      3, MPI_INTEGER,
     +                     source,      0, 
     +                       comm, status,        ierr)
       newpnv=newtmpnv(1)
       newbnv=newtmpnv(2)
       newenv=newtmpnv(3)
       call MPI_SENDRECV(    W_VIR,  nbufv, MPI_INTEGER,
     +                        dest,      0, 
     +                   W_VIR_new,  nbufv, MPI_INTEGER,
     +                      source,      0, 
     +                        comm, status,        ierr)
       call MPI_SENDRECV(    E_VIR, nbufv*Nk_irr, MPI_DOUBLE_PRECISION,
     +                        dest,         0, 
     +                   E_VIR_new, nbufv*Nk_irr, MPI_DOUBLE_PRECISION,
     +                      source,         0, 
     +                        comm,    status,                  ierr)
       call MPI_SENDRECV(    P_VIR, nbufv*NTK, MPI_DOUBLE_COMPLEX,
     +                        dest,         0, 
     +                   P_VIR_new, nbufv*NTK, MPI_DOUBLE_COMPLEX, 
     +                      source,         0, 
     +                        comm,    status,                  ierr)
       call MPI_SENDRECV(    VIR, nbufv*Nk_irr*NTG, MPI_DOUBLE_COMPLEX,
     +                      dest,             0, 
     +                    VIR_new, nbufv*Nk_irr*NTG, MPI_DOUBLE_COMPLEX,
     +                     source,             0,
     +                      comm,        status,             ierr)
       pnv=newpnv
       bnv=newbnv
       env=newenv
       W_VIR(1:nbufv)=W_VIR_new(1:nbufv)
       E_VIR(1:nbufv,1:Nk_irr)=E_VIR_new(1:nbufv,1:Nk_irr)
       P_VIR(1:nbufv,1:NTK)=P_VIR_new(1:nbufv,1:NTK)
       VIR(1:NTG,1:nbufv,1:Nk_irr)=VIR_new(1:NTG,1:nbufv,1:Nk_irr)
       if(myrank.eq.0)write(6,*)'#MPI end' 
       enddo!irot 
!--
       call MPI_REDUCE(pchiqw,chiqw,NG_for_eps*NG_for_eps*nen,
     +      MPI_DOUBLE_COMPLEX,MPI_SUM,0,comm,ierr)
!--
!output 
       if(myrank.eq.0)then 
        chiqw(:,:,:)=chiqw(:,:,:)/VOLUME 
        ierr=CHDIR("./dir-eps") 
        call system('pwd') 
!--
!2017413 
        do ik=1,Nk_irr 
         if(ABS(SKI_list(1,ik)-q1)<1.D-6.and.
     +      ABS(SKI_list(2,ik)-q2)<1.D-6.and. 
     +      ABS(SKI_list(3,ik)-q3)<1.D-6)then 
          qnum=ik 
         endif 
        enddo 
        write(6,*) 'qnum=',qnum
        write(dirname,"('q',i3.3)")qnum 
        ierr=CHDIR(dirname) 
        call system('pwd') 
!--
!OPEN(301,W,FILE='dat.sq') 
        OPEN(301,FILE='dat.sq') 
        rewind(301)
        write(301,*) q1,q2,q3 
        close(301) 
!--
        allocate(epsqw(NG_for_eps,NG_for_eps,nen));epsqw(:,:,:)=0.0d0 
        allocate(eps_rpa(NG_for_eps,NG_for_eps));eps_rpa(:,:)=0.0d0 
        allocate(chi0(NG_for_eps,NG_for_eps));chi0(:,:)=0.0d0 
!--
!OPEN(60000-,W,FILE='epsqw(G,G)) 
        file_num_eps_base=file_num_eps_base_start!20171011  
        do ie=1,nen 
         file_num_eps_ie=file_num_eps_base+ie!-1 
         write(filename,'("dat.epsqw.",i6.6)')file_num_eps_ie   
         OPEN(file_num_eps_ie,FILE=filename) 
         rewind(file_num_eps_ie)
         write(file_num_eps_ie,'(a)')'#epsqw diagonal term' 
         write(file_num_eps_ie,'(a)')
     +   '#1:q+G [au], 2:Re(eps-1), 3:Im(eps-1)'
        enddo!ie 
!--
        do ie=1,nen 
         file_num_chi_ie=0!file_num_chi_base+ie-1
         file_num_eps_ie=file_num_eps_base+ie!-1 
         chi0(:,:)=chiqw(:,:,ie)
         call calc_eps_rpa(NTG,NG_for_eps,LG0(1,1,iq),
     +                     q1,q2,q3,b1(1),b2(1),b3(1),
     +                     chi0(1,1),tpi,eps_rpa(1,1), 
     +                     file_num_chi_ie,file_num_eps_ie)  
         epsqw(:,:,ie)=eps_rpa(:,:)
         close(file_num_eps_ie)
        enddo!ie 
!--
!chiqw 
!       write(file_num_chi)(((chiqw(ig,jg,ie),ig=1,NG_for_eps)
!     +                            ,jg=1,NG_for_eps),ie=1,nen)   
!--
!OPEN(600,W,FILE='epsqw',FORM='unformatted') 
        file_num_eps=file_num_eps_start 
        write(filename,'("dat.epsqw.",i3.3)')file_num_eps  
        OPEN(file_num_eps,FILE=filename,FORM='unformatted') 
        rewind(file_num_eps)
        write(file_num_eps)(((epsqw(ig,jg,ie),ig=1,NG_for_eps)
     +                             ,jg=1,NG_for_eps),ie=1,nen) 
        close(file_num_eps) 
!--
!OPEN(400,W,FILE='dat.log') 
        file_num_log=file_num_log_start 
        write(filename,'("dat.log.",i3.3)')file_num_log 
        OPEN(file_num_log,FILE=filename) 
        rewind(file_num_log)
        write(file_num_log,*)'#'
        write(file_num_log,*)'#finish chi0 evaluation for',qnum!iq 
        write(file_num_log,*)'#' 
        close(file_num_log) 
!--
        deallocate(epsqw,eps_rpa,chi0)
        ierr=CHDIR("../..")!2017413  
        call system('pwd') 
       endif 
       deallocate(chiqw,pchiqw,i_for_pair,j_for_pair) 
!--
!q.eq.0 
      elseif(q1==0.0d0.and.q2==0.0d0.and.q3==0.0d0)then 
       if(myrank.eq.0)then 
        write(6,'(i5,3f15.8)') iq,q1,q2,q3 
       endif 
       do igL=1,NG_for_eps
        igL1=LG0(1,igL,iq) 
        igL2=LG0(2,igL,iq) 
        igL3=LG0(3,igL,iq) 
        if(igL1==0.and.igL2==0.and.igL3==0)then
         No_G_0=igL 
         write(6,*)'No. for G=0',No_G_0
        endif
       enddo!igL
       allocate(chiqwgm(NG_for_eps,NG_for_eps,nen,3))
       chiqwgm(:,:,:,:)=0.0D0        
       allocate(pchiqwgm(NG_for_eps,NG_for_eps,nen,3))
       pchiqwgm(:,:,:,:)=0.0D0 
       allocate(i_for_pair(N_CALC_BAND*N_CALC_BAND))
       allocate(j_for_pair(N_CALC_BAND*N_CALC_BAND))
       allocate(i_for_intra(N_CALC_BAND*N_CALC_BAND))
       allocate(j_for_intra(N_CALC_BAND*N_CALC_BAND))
!wavefunction-rotation loop  
!wpl calc 20130308 
       wpl(:)=0.0d0;pwpl(:)=0.0d0
      do irot=0,nproc-1
!double loop-->single loop 
       i_for_pair(:)=0
       j_for_pair(:)=0
       ij=0
       ij_intra=0
       i_for_intra(:)=0
       j_for_intra(:)=0
       do jb=1,pno 
       do ib=1,pnv
        fs1=W_VIR(ib)
        fs2=W_OCC(jb)
!--
!constrained excitation 
        bci=EV_AVE(ib) 
        bcj=EO_AVE(jb) 
        if(abs(bci-bcj)>=delta_ex)cycle 
!--
!        if(myrank.eq.0)then 
!         write(6,*) abs(bci-bcj)*au,delta_ex*au  
!        endif 
!--
!1:window-cRPA not default
       if(fs1==1.and.fs2==1.and.flg_cRPA_ewin==1)cycle
!--
        fs1=T_VIR(ib)
        fs2=T_OCC(jb)
        if((fs1==fs2).and.fs1.ne.0)then 
!wpl calc 20130308
         ij_intra=ij_intra+1
         write(6,*)'iandj',fs1,fs2 
         i_for_intra(ij_intra)=ib  
         j_for_intra(ij_intra)=jb  
        else 
         ij=ij+1
         i_for_pair(ij)=ib  
         j_for_pair(ij)=jb  
         write(6,*)'ij=',ij,ib,jb 
        endif 
       enddo 
       enddo 
       N_PAIR=ij 
       if(myrank.eq.0)write(6,*)'FOR q=0,N_PAIR=',N_PAIR 
       N_INTRA=ij_intra 
!cRPA
       if(flg_cRPA==1)then 
        if(myrank.eq.0)write(6,*)'In cRPA,N_INTRA=0 forcidly'
        N_INTRA=0
       endif 
       if(myrank.eq.0)write(6,*)'FOR q=0, N_INTRA=',N_INTRA 
!wpl calc 20130308
       allocate(pk_1D(NTK,3))
       allocate(pk_3D(nkb1,nkb2,nkb3,3))  
       allocate(fk_1D(NTK))
       allocate(fk_3D(nkb1,nkb2,nkb3))  
       allocate(gk_1D(NTK)) 
       allocate(gk_3D(nkb1,nkb2,nkb3)) 
       allocate(xo(nkb1,nkb2,nkb3))  
       do ij_intra=1,N_INTRA 
        ib=i_for_intra(ij_intra)
        jb=j_for_intra(ij_intra)
        gk_1D(:)=0.0d0 
        fk_1D(:)=0.0d0
        pk_1D(:,:)=0.0d0
        do ik=1,NTK 
         ikq=ik 
         SKT(:)=SK0(:,ik) 
!C0_K(:)=OCC(:,jb,ik) 
         ikir=numirr(ik)
         iop=numrot(ik) 
         allocate(C0_K(NTG));C0_K(:)=0.0d0 
         call make_C0(NTG,trs(ik),NG0(ik),KG0(1,1,ik),RW(1,ik),
     +   rginv(1,1,iop),pg(1,iop),nnp,L1,L2,L3,
     +   packing(-L1,-L2,-L3,ikir),OCC(1,jb,ikir),C0_K(1)) 
!C0_KQ(:)=VIR(:,ib,ikq) 
         ikir=numirr(ikq)
         iop=numrot(ikq) 
         allocate(C0_KQ(NTG));C0_KQ(:)=0.0d0 
         call make_C0(NTG,trs(ikq),NG0(ikq),KG0(1,1,ikq),RW(1,ikq),     
     +   rginv(1,1,iop),pg(1,iop),nnp,L1,L2,L3,
     +   packing(-L1,-L2,-L3,ikir),VIR(1,ib,ikir),C0_KQ(1)) 
         call calc_VMaa(NTK,NTG,NG0(1),KG0(1,1,1),C0_K(1),C0_KQ(1),ik, 
     +                  b1(1),b2(1),b3(1),SKT(1),vm(1))  
         deallocate(C0_K,C0_KQ) 
!20170424 
        NLvm(:)=0.0d0 
        NLvm(:)=MAT_NLr_rNL(:,ib+(N_docc+bnv-1),jb+(bno-1),numirr(ik))
        pk_1D(ik,:)=-ci*(vm(:)-NLvm(:)) 
!--
        fk_1D(ik)=1.0d0
        ikir=numirr(ik) 
        gk_1D(ik)=cmplx(E_OCC(jb,ikir)-FermiEnergy,-delt)
        enddo!ik  
        pk_3D(:,:,:,:)=0.0d0
        fk_3D(:,:,:)=0.0d0
        gk_3D(:,:,:)=0.0d0 
        do ikb3=1,nkb3
         do ikb2=1,nkb2
          do ikb1=1,nkb1
           ik=index_kpt(ikb1,ikb2,ikb3)
           pk_3D(ikb1,ikb2,ikb3,:)=pk_1D(ik,:)
           fk_3D(ikb1,ikb2,ikb3)=fk_1D(ik)
           gk_3D(ikb1,ikb2,ikb3)=gk_1D(ik)
          enddo 
         enddo 
        enddo 
!--
        xo(:,:,:)=0.0d0 
        call ttrhdrn_simple(dmna,dmnr,nkb1,nkb2,nkb3,imt1(1),
     +                      fk_3D(1,1,1),gk_3D(1,1,1),xo(1,1,1))
!--
        do ikb3=1,nkb3
         do ikb2=1,nkb2
          do ikb1=1,nkb1
           ik=index_kpt(ikb1,ikb2,ikb3)
           pwpl(:)
     +    =pwpl(:)+dimag(xo(ikb1,ikb2,ikb3))
     +    *pk_3D(ikb1,ikb2,ikb3,:)*pk_3D(ikb1,ikb2,ikb3,:) 
          enddo 
         enddo 
        enddo 
!--
!20170324 
!        xo_1D(:)=0.0d0 
!        do ikb3=1,nkb3
!         do ikb2=1,nkb2
!          do ikb1=1,nkb1
!           ik=index_kpt(ikb1,ikb2,ikb3)
!           xo_1D(ik)=xo(ikb1,ikb2,ikb3)
!          enddo 
!         enddo 
!        enddo 
!--
!       xo_1D_av(:)=0.0d0 
!       do ikir=1,Nk_irr 
!        nk=0
!        SUM_CMPX=0.0d0 
!        do ik=1,NTK 
!         if(numirr(ik)==ikir)then 
!          SUM_CMPX=SUM_CMPX+xo_1D(ik)
!          nk=nk+1
!         endif 
!        enddo!ik
!        xo_1D_av(ikir)=SUM_CMPX/dble(nk)
!       enddo!ikir
!--
!        do ik=1,NTK 
!         ikir=numirr(ik) 
!         pwpl=pwpl+dimag(xo_1D(ik))*pk_1D(ik)*pk_1D(ik) 
!!        pwpl=pwpl+dimag(xo_1D_av(ikir))*pk_1D(ik)*pk_1D(ik)!20170324  
!        enddo 
!--
       enddo!ij_intra 
       deallocate(pk_1D,pk_3D,fk_1D,fk_3D,gk_1D,gk_3D,xo)
!$OMP PARALLEL PRIVATE(ij,ib,jb,ix,
!$OMP&                 e1_1D,e2_1D,ik,shift_G,ikq, 
!$OMP&                 C0_K,C0_KQ,dsgn,kpq, 
!$OMP&                 wfunc,fftwk,m_tmp,ismgm,vm,NLvm,del_eps, 
!$OMP&                 e1_3D,e2_3D,ikb3,ikb2,ikb1, 
!$OMP&                 xow,xow_1D,ca1,ie,ig,jg, 
!$OMP&                 pc0dgm,iomp,ikir,iop) 
       allocate(pc0dgm(NG_for_eps,NG_for_eps,nen,3));pc0dgm=(0.d0,0.d0) 
       allocate(ismgm(NTK,NG_for_eps,3))
!--
       allocate(m_tmp(NG_for_eps))
       allocate(fftwk(Nl123*2)) 
       allocate(wfunc(Nl123*2)) 
       allocate(C0_K(NTG))    
       allocate(C0_KQ(NTG))    
       allocate(kpq(NTK)) 
       allocate(e1_1D(NTK))
       allocate(e2_1D(NTK))
       allocate(e1_3D(nkb1,nkb2,nkb3)) 
       allocate(e2_3D(nkb1,nkb2,nkb3))
       allocate(xow(nen,nkb1,nkb2,nkb3)) 
       allocate(xow_1D(NTK,nen)) 
       allocate(ca1(4*nen))  
!$OMP DO
       do ij=1,N_PAIR 
        ib=i_for_pair(ij) 
        jb=j_for_pair(ij) 
        do ik=1,NTK 
         ikq=ik 
         kpq(ik)=ikq
         shift_G(:)=0 
         ikir=numirr(ik)
         iop=numrot(ik) 
!C0_K(:)=OCC(:,jb,ik) 
         C0_K(:)=0.0d0 
         call make_C0(NTG,trs(ik),NG0(ik),KG0(1,1,ik),RW(1,ik),
     +   rginv(1,1,iop),pg(1,iop),nnp,L1,L2,L3,
     +   packing(-L1,-L2,-L3,ikir),OCC(1,jb,ikir),C0_K(1)) 
         ikir=numirr(ikq)
         iop=numrot(ikq) 
!C0_KQ(:)=VIR(:,ib,ikq) 
         C0_KQ(:)=0.0d0 
         call make_C0(NTG,trs(ikq),NG0(ikq),KG0(1,1,ikq),RW(1,ikq),     
     +   rginv(1,1,iop),pg(1,iop),nnp,L1,L2,L3,
     +   packing(-L1,-L2,-L3,ikir),VIR(1,ib,ikir),C0_KQ(1)) 
         call calc_InterStateMatrix(NTK,NTG,NG0(1),KG0(1,1,1),
     +                              C0_K(1),C0_KQ(1),
     +                              ik,ikq,nwx2,nwy2,nwz2,
     +                              nfft1,nfft2,Nl123,
     +                              wfunc(1),fftwk(1),fs,LG0(1,1,iq),
     +                              NG_for_eps,shift_G(1),m_tmp(1))
         call calc_VMab(NTK,NTG,NG0(1),KG0(1,1,1),C0_K(1),C0_KQ(1),ik, 
     +                  b1(1),b2(1),b3(1),vm(1))  
!--
         del_eps=E_VIR(ib,numirr(ikq))-E_OCC(jb,numirr(ik))
!20140915 
         if(abs(del_eps)<delt)then 
          write(6,*)'del_eps=',del_eps 
          del_eps=delt 
         endif 
!2017022 
         NLvm(:)=0.0d0 
         NLvm(:)=MAT_NLr_rNL(:,ib+(N_docc+bnv-1),jb+(bno-1),numirr(ik))
         vm(:)=-ci*(vm(:)-NLvm(:))/del_eps!correct  
!20170424 
         do ix=1,3
          m_tmp(No_G_0)=vm(ix)
          if(abs(del_eps)<=delt)then 
           m_tmp(No_G_0)=0.0d0 
          endif 
          ismgm(ik,:,ix)=m_tmp(:) 
         enddo!ix  
        enddo!ik 
!Using ttrhdrn_chi0 
        e1_1D(:)=0.0d0
        e2_1D(:)=0.0d0 
        do ik=1,NTK
         ikq=kpq(ik)
         e1_1D(ik)=cmplx(E_OCC(jb,numirr(ik)),0.0d0) 
         e2_1D(ik)=cmplx(E_VIR(ib,numirr(ikq)),0.0d0) 
        enddo 
!avoid degeneracy 20120916 
!       do ik=1,NTK 
!        if(abs(e1_1D(ik)-e2_1D(ik))<delt) then 
!         e1_1D(ik)=e1_1D(ik)-0.5d0*delt !*100.0d0 
!         e2_1D(ik)=e2_1D(ik)+0.5d0*delt !*100.0d0
!        endif 
!       enddo 
!--
        e1_3D(:,:,:)=0.0d0 
        e2_3D(:,:,:)=0.0d0
        do ikb3=1,nkb3
         do ikb2=1,nkb2
          do ikb1=1,nkb1
           ik=index_kpt(ikb1,ikb2,ikb3)
           e1_3D(ikb1,ikb2,ikb3)=e1_1D(ik)
           e2_3D(ikb1,ikb2,ikb3)=e2_1D(ik)
          enddo 
         enddo 
        enddo 
!freq dep 
        xow(:,:,:,:)=0.0d0 
        dsgn=1.0d0
        ca1(:)=0.0d0
        call ttrhdrn_chi0(dmna,dmnr,nkb1,nkb2,nkb3,imt1(1),
     +                   e1_3D(1,1,1),e2_3D(1,1,1),
     +                   FermiEnergy,FermiEnergy,delt,dsgn,
     +                   nen,em(1),ca1(1),xow(1,1,1,1)) 
        dsgn=-1.0d0
        ca1(:)=0.0d0
        call ttrhdrn_chi0(dmna,dmnr,nkb1,nkb2,nkb3,imt1(1),
     +                   e1_3D(1,1,1),e2_3D(1,1,1),
     +                   FermiEnergy,FermiEnergy,delt,dsgn,
     +                   nen,em(1),ca1(1),xow(1,1,1,1)) 
!--
        xow_1D(:,:)=0.0d0 
        do ikb3=1,nkb3
         do ikb2=1,nkb2
          do ikb1=1,nkb1
           ik=index_kpt(ikb1,ikb2,ikb3)
           do ie=1,nen 
            xow_1D(ik,ie)=xow(ie,ikb1,ikb2,ikb3)
           enddo 
          enddo 
         enddo 
        enddo 
!org freq dep 
        do ie=1,nen 
         do ig=1,NG_for_eps
          do jg=1,NG_for_eps
           do ik=1,NTK 
            ikq=kpq(ik)
            pc0dgm(ig,jg,ie,:)
     +     =pc0dgm(ig,jg,ie,:)
     +     +xow_1D(ik,ie)
     +     *CONJG(ismgm(ik,ig,:))*ismgm(ik,jg,:)/dble(NTK)
     +     *(1.0d0-P_OCC(jb,ik)*P_VIR(ib,ikq)) 
           enddo!ik 
          enddo!jg 
         enddo!ig 
        enddo!ie
!--
!20170321 
!       xow_1D_av(:,:)=0.0d0 
!       do ie=1,nen 
!        do ikir=1,Nk_irr 
!         nk=0
!         SUM_CMPX=0.0d0 
!         do ik=1,NTK 
!          if(numirr(ik)==ikir)then 
!           SUM_CMPX=SUM_CMPX+xow_1D(ik,ie)
!           nk=nk+1
!          endif 
!         enddo!ik
!         xow_1D_av(ikir,ie)=SUM_CMPX/dble(nk)
!        enddo!ikir
!       enddo!ie 
!--
!20170321 averaging: In new ttrhdrn, averaging is no need.
!       
!       do ie=1,nen 
!        do ig=1,NG_for_eps
!         do jg=1,NG_for_eps
!          do ik=1,NTK 
!           ikir=numirr(ik) 
!           ikq=kpq(ik)
!           pc0d(ig,jg,ie)
!     +    =pc0d(ig,jg,ie)
!     +    +xow_1D_av(ikir,ie)
!!    +    *CONJG(ism(ik,ig))*ism(ik,jg)/dble(NTK)
!     +    *CONJG(ism_st(ik,ig))*ism(ik,jg)/dble(NTK)
!     +    *(1.0d0-P_OCC(jb,ik)*P_VIR(ib,ikq)) 
!          enddo!ik 
!         enddo!jg 
!        enddo!ig 
!       enddo!ie 
!---
       iomp=omp_get_thread_num() 
       if(myrank.eq.0.and.iomp.eq.0)write(6,*)'#',ib,jb
!      if(myrank.eq.0.and.iomp.eq.0)write(6,*)'av'!20170321  
!---
      enddo!ij 
!$OMP END DO
!$OMP CRITICAL
      pchiqwgm=pchiqwgm+pc0dgm 
!$OMP END CRITICAL
      deallocate(pc0dgm,fftwk,wfunc,m_tmp,ismgm)
      deallocate(C0_K,C0_KQ,kpq,e1_1D,e2_1D,e1_3D,e2_3D,xow,xow_1D,ca1) 
!     deallocate(xow_1D_av,m_tmp_st,ism_st)!20170322
!$OMP END PARALLEL
!---
       if(myrank.eq.0)write(6,*)'#MPI start' 
       call MPI_BARRIER(comm,ierr)
       dest=modulo(myrank+1,nproc)
       source=modulo(myrank-1,nproc)
       tmpnv(1)=pnv
       tmpnv(2)=bnv
       tmpnv(3)=env
       call MPI_SENDRECV(   tmpnv,      3, MPI_INTEGER, 
     +                       dest,      0, 
     +                   newtmpnv,      3, MPI_INTEGER,
     +                     source,      0, 
     +                       comm, status,        ierr)
       newpnv = newtmpnv(1)
       newbnv = newtmpnv(2)
       newenv = newtmpnv(3)
       call MPI_SENDRECV(    W_VIR,  nbufv, MPI_INTEGER,
     +                        dest,      0, 
     +                   W_VIR_new,  nbufv, MPI_INTEGER,
     +                      source,      0, 
     +                        comm, status,        ierr)
       call MPI_SENDRECV(    T_VIR,  nbufv, MPI_INTEGER,
     +                        dest,      0, 
     +                   T_VIR_new,  nbufv, MPI_INTEGER,
     +                      source,      0, 
     +                        comm, status,        ierr)
       call MPI_SENDRECV(   EV_AVE,  nbufv, MPI_DOUBLE_PRECISION,
     +                        dest,      0, 
     +                  EV_AVE_new,  nbufv, MPI_DOUBLE_PRECISION,
     +                      source,      0, 
     +                        comm,    status,                  ierr)
       call MPI_SENDRECV(    E_VIR, nbufv*Nk_irr, MPI_DOUBLE_PRECISION,
     +                        dest,         0, 
     +                   E_VIR_new, nbufv*Nk_irr, MPI_DOUBLE_PRECISION,
     +                      source,         0, 
     +                        comm,    status,                  ierr)
       call MPI_SENDRECV(    P_VIR, nbufv*NTK, MPI_DOUBLE_COMPLEX, 
     +                        dest,         0, 
     +                   P_VIR_new, nbufv*NTK, MPI_DOUBLE_COMPLEX, 
     +                      source,         0, 
     +                        comm,    status,                  ierr)
       call MPI_SENDRECV(    VIR, nbufv*Nk_irr*NTG, MPI_DOUBLE_COMPLEX,
     +                      dest,             0, 
     +                   VIR_new, nbufv*Nk_irr*NTG, MPI_DOUBLE_COMPLEX,
     +                    source,             0,
     +                      comm,        status,             ierr)
       pnv=newpnv
       bnv=newbnv
       env=newenv
       W_VIR(1:nbufv)=W_VIR_new(1:nbufv)
       T_VIR(1:nbufv)=T_VIR_new(1:nbufv)
       EV_AVE(1:nbufv)=EV_AVE_new(1:nbufv)
       E_VIR(1:nbufv,1:Nk_irr)=E_VIR_new(1:nbufv,1:Nk_irr)
       P_VIR(1:nbufv,1:NTK)=P_VIR_new(1:nbufv,1:NTK)
       VIR(1:NTG,1:nbufv,1:Nk_irr)=VIR_new(1:NTG,1:nbufv,1:Nk_irr)
       if(myrank.eq.0) write(6,*)'#MPI end' 
      enddo!irot 
!--
      call MPI_REDUCE(pchiqwgm,chiqwgm,NG_for_eps*NG_for_eps*nen*3,
     +     MPI_DOUBLE_COMPLEX,MPI_SUM,0,comm,ierr)
      call MPI_REDUCE(pwpl,wpl,3,MPI_DOUBLE_COMPLEX,MPI_SUM,0,comm,ierr)
!--
!output 
      if(myrank.eq.0)then 
       chiqwgm(:,:,:,:)=chiqwgm(:,:,:,:)/VOLUME 
!--
!wpl calc 20130308 
       wpl(:)=8.0d0*wpl(:)/VOLUME/dble(NTK) 
       wd(:)=dsqrt(dble(wpl(:))) 
       write(6,'(a10,3f15.10)')'wd in eV',wd*au 
!--
       ierr=CHDIR("./dir-eps") 
       call system('pwd') 
!--
!2017413 
       do ik=1,Nk_irr 
        if(ABS(SKI_list(1,ik)-q1)<1.D-6.and.
     +     ABS(SKI_list(2,ik)-q2)<1.D-6.and. 
     +     ABS(SKI_list(3,ik)-q3)<1.D-6)then 
         qnum=ik 
        endif 
       enddo 
       write(6,*)'qnum=',qnum
       write(dirname,"('q',i3.3)")qnum 
       ierr=CHDIR(dirname) 
       call system('pwd') 
!--
!OPEN(301,W,FILE='dat.sq') 
        OPEN(301,FILE='dat.sq') 
        rewind(301)
        write(301,*) q1,q2,q3 
        close(301) 
!--
       allocate(epsqwgm(NG_for_eps,NG_for_eps,nen,3))
       epsqwgm(:,:,:,:)=0.0d0 
       allocate(eps_rpa(NG_for_eps,NG_for_eps));eps_rpa(:,:)=0.0d0 
       allocate(chi0(NG_for_eps,NG_for_eps));chi0(:,:)=0.0d0 
!--
!OPEN(60000-,W,FILE='epsqw(G,G)) 
       do ix=1,3
        file_num_eps_base=file_num_eps_base_start+(ix-1)*1000   
        do ie=1,nen 
         file_num_eps_ie=file_num_eps_base+ie!-1 
         write(filename,'("dat.epsqw.",i6.6)')file_num_eps_ie   
         OPEN(file_num_eps_ie,FILE=filename) 
         rewind(file_num_eps_ie)
         write(file_num_eps_ie,'(a)')'#epsqw diagonal term' 
         write(file_num_eps_ie,'(a)')
     +   '#1:q+G [au], 2:Re(eps-1), 3:Im(eps-1)' 
        enddo!ie 
       enddo!ix
!--
       do ix=1,3
        file_num_eps_base=file_num_eps_base_start+(ix-1)*1000   
        do ie=1,nen 
         file_num_chi_ie=0 
         file_num_eps_ie=file_num_eps_base+ie!-1 
         chi0(:,:)=chiqwgm(:,:,ie,ix)
         call calc_eps_rpa_q_0(NTG,NG_for_eps,LG0(1,1,iq),
     +                     q1,q2,q3,b1(1),b2(1),b3(1),
     +                     chi0(1,1),tpi,No_G_0,eps_rpa(1,1),
     +                     wd(ix),delt,em(ie),     
     +                     file_num_chi_ie,file_num_eps_ie)  
         epsqwgm(:,:,ie,ix)=eps_rpa(:,:)
         close(file_num_eps_ie)
        enddo!ie
       enddo!ix 
!--
!chiqw 
!      write(file_num_chi)(((chiqw(ig,jg,ie),ig=1,NG_for_eps)
!     +                           ,jg=1,NG_for_eps),ie=1,nen) 
!--
!OPEN(600,W,FILE='epsqw',FORM='unformatted') 
       do ix=1,3 
        file_num_eps=file_num_eps_start+(ix-1)  
        write(filename,'("dat.epsqw.",i3.3)')file_num_eps  
        OPEN(file_num_eps,FILE=filename,FORM='unformatted') 
        rewind(file_num_eps)
        write(file_num_eps)(((epsqwgm(ig,jg,ie,ix),ig=1,NG_for_eps)
     +                               ,jg=1,NG_for_eps),ie=1,nen) 
        close(file_num_eps)
       enddo 
!--
!OPEN(400,W,FILE='dat.log') 
       file_num_log=file_num_log_start 
       write(filename,'("dat.log.",i3.3)')file_num_log 
       OPEN(file_num_log,FILE=filename) 
       rewind(file_num_log)
       write(file_num_log,*)'#'
       write(file_num_log,*)'#finish chi0 evaluation for',qnum 
       write(file_num_log,*)'#' 
       close(file_num_log)
!--
       deallocate(epsqwgm,eps_rpa,chi0)
       ierr=CHDIR("../..")!20170413 
       call system('pwd') 
      endif!myrank==0  
      deallocate(chiqwgm,pchiqwgm,i_for_pair,j_for_pair)
      deallocate(i_for_intra,j_for_intra) 
!--
      endif!q=0 or q/=0
      enddo!iq 
!     CLOSE(file_num_chi)
!     CLOSE(file_num_eps)
!--
      call MPI_BARRIER(comm_glb,ierr)
      end_time=MPI_Wtime()
      diff_time=end_time-start_time 
      if(myrank_glb.eq.0) then 
       call system('pwd') 
       write(6,*)'#TOTAL TIME=',diff_time 
      endif 
!--
      call MPI_FINALIZE(ierr)
      STOP 
      END              
