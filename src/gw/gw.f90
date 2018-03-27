PROGRAM GW
!
use mpi
use fft_3d 
use m_rdinput       
use m_rd_dat_wfn
use m_rd_dat_wan 
use m_rd_dat_eps 
include "config.h" 
!
write(6,*)'calc_type is off-diag' 
!
!mpi 
!
comm=MPI_COMM_WORLD
call MPI_INIT(ierr)
call MPI_COMM_RANK(comm,myrank,ierr)
call MPI_COMM_SIZE(comm,nproc,ierr)
write(6,'(a10,i10,a10,i10)')'myrank=',myrank,'nproc=',nproc
file_id=9000+myrank 
call MPI_BARRIER(comm,ierr)
start_time=MPI_Wtime() 
!
!read input
!
if(myrank.eq.0)then 
 call read_input 
 Green_func_delt=Green_func_delt/au  
 idlt=idlt/au
endif 
call MPI_BARRIER(comm,ierr)
!
call MPI_Bcast(CALC_IFREQ,1,MPI_INTEGER,0,comm,ierr) 
call MPI_Bcast(IX_INTJ_MIN,1,MPI_INTEGER,0,comm,ierr) 
call MPI_Bcast(IX_INTJ_MIN,1,MPI_INTEGER,0,comm,ierr) 
call MPI_Bcast(IY_INTJ_MIN,1,MPI_INTEGER,0,comm,ierr) 
call MPI_Bcast(IY_INTJ_MAX,1,MPI_INTEGER,0,comm,ierr) 
call MPI_Bcast(IZ_INTJ_MIN,1,MPI_INTEGER,0,comm,ierr) 
call MPI_Bcast(IZ_INTJ_MAX,1,MPI_INTEGER,0,comm,ierr) 
call MPI_Bcast(WCUT_MVMC,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
call MPI_Bcast(JCUT_MVMC,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
call MPI_Bcast(Green_func_delt,1,MPI_DOUBLE_PRECISION,0,comm,ierr) 
call MPI_Bcast(idlt,1,MPI_DOUBLE_PRECISION,0,comm,ierr) 
call MPI_Bcast(Rc_range_spacing,1,MPI_INTEGER,0,comm,ierr) 
call MPI_Bcast(N_sym_points,1,MPI_INTEGER,0,comm,ierr) 
call MPI_Bcast(Ndiv,1,MPI_INTEGER,0,comm,ierr) 
call MPI_Bcast(reading_sk_format,1,MPI_INTEGER,0,comm,ierr) 
call MPI_BARRIER(comm,ierr)
!
if(myrank/=0) allocate(dense(3))
if(myrank/=0) allocate(SK_sym_pts(3,N_sym_points))  
call MPI_Bcast(dense,3,MPI_INTEGER,0,comm,ierr)
call MPI_Bcast(SK_sym_pts,3*N_sym_points,MPI_DOUBLE_PRECISION,0,comm,ierr)
call MPI_BARRIER(comm,ierr)
!
!read wfn 
!
if(myrank.eq.0)then 
 call rd_dat_symmetry 
 call rd_dat_bandcalc 
 call rd_dat_lattice 
 call rd_dat_sample_k 
 call rd_dat_nkm 
 call rd_dat_kg 
 call rd_dat_eigenvalue  
 call rd_dat_wavefunction 
endif 
call MPI_BARRIER(comm,ierr)
!
!sym(100)
!
call MPI_Bcast(nsymq,1,MPI_INTEGER,0,comm,ierr) 
call MPI_Bcast(nnp,1,MPI_INTEGER,0,comm,ierr) 
!
!bandcalc(117) 
!
call MPI_Bcast(Ecut_for_psi,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
call MPI_Bcast(FermiEnergy,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
call MPI_Bcast(Etot,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
!
!avec(105)
!
call MPI_Bcast(a1,3,MPI_DOUBLE_PRECISION,0,comm,ierr) 
call MPI_Bcast(a2,3,MPI_DOUBLE_PRECISION,0,comm,ierr) 
call MPI_Bcast(a3,3,MPI_DOUBLE_PRECISION,0,comm,ierr) 
call MPI_Bcast(b1,3,MPI_DOUBLE_PRECISION,0,comm,ierr) 
call MPI_Bcast(b2,3,MPI_DOUBLE_PRECISION,0,comm,ierr) 
call MPI_Bcast(b3,3,MPI_DOUBLE_PRECISION,0,comm,ierr) 
call MPI_Bcast(VOLUME,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
call MPI_Bcast(a,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
call MPI_Bcast(b,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
call MPI_Bcast(c,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
call MPI_Bcast(alp,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
call MPI_Bcast(bet,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
call MPI_Bcast(gmm,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
!
!sample-k(101)
!
call MPI_Bcast(Nk_irr,1,MPI_INTEGER,0,comm,ierr) 
call MPI_Bcast(NTK,1,MPI_INTEGER,0,comm,ierr) 
call MPI_Bcast(nkb1,1,MPI_INTEGER,0,comm,ierr) 
call MPI_Bcast(nkb2,1,MPI_INTEGER,0,comm,ierr) 
call MPI_Bcast(nkb3,1,MPI_INTEGER,0,comm,ierr) 
call MPI_Bcast(Na1,1,MPI_INTEGER,0,comm,ierr) 
call MPI_Bcast(Na2,1,MPI_INTEGER,0,comm,ierr) 
call MPI_Bcast(Na3,1,MPI_INTEGER,0,comm,ierr) 
!
!kg(104)
!
call MPI_Bcast(L1,1,MPI_INTEGER,0,comm,ierr) 
call MPI_Bcast(L2,1,MPI_INTEGER,0,comm,ierr) 
call MPI_Bcast(L3,1,MPI_INTEGER,0,comm,ierr) 
!
!fft 
!
call MPI_Bcast(nwx2,1,MPI_INTEGER,0,comm,ierr) 
call MPI_Bcast(nwy2,1,MPI_INTEGER,0,comm,ierr) 
call MPI_Bcast(nwz2,1,MPI_INTEGER,0,comm,ierr) 
!
!nkm(132)  
!
call MPI_Bcast(NTG,1,MPI_INTEGER,0,comm,ierr) 
!
!eigenvalue(111) 
!
call MPI_Bcast(NTB,1,MPI_INTEGER,0,comm,ierr) 
!
!cir(102) 
!
call MPI_Bcast(ncomp,1,MPI_INTEGER,0,comm,ierr) 
call MPI_BARRIER(comm,ierr)
!
if(myrank/=0) allocate(rg(3,3,nsymq))
if(myrank/=0) allocate(pg(3,nsymq))
if(myrank/=0) allocate(rginv(3,3,nsymq))  
if(myrank/=0) allocate(SKI(3,Nk_irr))  
if(myrank/=0) allocate(SK0(3,NTK))  
if(myrank/=0) allocate(numirr(NTK)) 
if(myrank/=0) allocate(numrot(NTK)) 
if(myrank/=0) allocate(trs(NTK)) 
if(myrank/=0) allocate(RW(3,NTK))
if(myrank/=0) allocate(numMK(Nk_irr))!20180316  
if(myrank/=0) allocate(NGI(Nk_irr)) 
if(myrank/=0) allocate(NG0(NTK))
if(myrank/=0) allocate(KGI(3,NTG,Nk_irr)) 
if(myrank/=0) allocate(KG0(3,NTG,NTK))
if(myrank/=0) allocate(packing(-L1:L1,-L2:L2,-L3:L3,Nk_irr))
if(myrank/=0) allocate(E_EIGI(NTB,Nk_irr))  
if(myrank/=0) allocate(CIR(NTG,NTB,Nk_irr)) 
call MPI_BARRIER(comm,ierr)
!
call MPI_Bcast(rg,3*3*nsymq,MPI_INTEGER,0,comm,ierr) 
call MPI_Bcast(pg,3*nsymq,MPI_INTEGER,0,comm,ierr) 
call MPI_Bcast(rginv,3*3*nsymq,MPI_DOUBLE_PRECISION,0,comm,ierr) 
call MPI_Bcast(SKI,3*Nk_irr,MPI_DOUBLE_PRECISION,0,comm,ierr) 
call MPI_Bcast(SK0,3*NTK,MPI_DOUBLE_PRECISION,0,comm,ierr) 
call MPI_Bcast(numirr,NTK,MPI_INTEGER,0,comm,ierr) 
call MPI_Bcast(numrot,NTK,MPI_INTEGER,0,comm,ierr) 
call MPI_Bcast(trs,NTK,MPI_INTEGER,0,comm,ierr) 
call MPI_Bcast(RW,3*NTK,MPI_INTEGER,0,comm,ierr) 
call MPI_Bcast(numMK,Nk_irr,MPI_INTEGER,0,comm,ierr)!20180316  
call MPI_Bcast(NGI,Nk_irr,MPI_INTEGER,0,comm,ierr) 
call MPI_Bcast(NG0,NTK,MPI_INTEGER,0,comm,ierr) 
call MPI_Bcast(KGI,3*NTG*Nk_irr,MPI_INTEGER,0,comm,ierr) 
call MPI_Bcast(KG0,3*NTG*NTK,MPI_INTEGER,0,comm,ierr) 
call MPI_Bcast(packing(-L1,-L2,-L3,1),(2*L1+1)*(2*L2+1)*(2*L3+1)*Nk_irr,MPI_INTEGER,0,comm,ierr) 
call MPI_Bcast(E_EIGI,NTB*Nk_irr,MPI_DOUBLE_PRECISION,0,comm,ierr) 
call MPI_Bcast(CIR,NTG*NTB*Nk_irr,MPI_DOUBLE_COMPLEX,0,comm,ierr) 
call MPI_BARRIER(comm,ierr)
!
!read wan 
!
if(myrank.eq.0)then 
 call rd_dat_ns_nb 
 call rd_dat_umat 
 call rd_dat_wan 
 call rd_dat_hmatr 
endif 
call MPI_BARRIER(comm,ierr)
!
!ns-nb(149)  
!
call MPI_Bcast(Mb,1,MPI_INTEGER,0,comm,ierr) 
call MPI_Bcast(Mt,1,MPI_INTEGER,0,comm,ierr) 
!
!umat(150)  
!
call MPI_Bcast(NWF,1,MPI_INTEGER,0,comm,ierr) 
call MPI_BARRIER(comm,ierr)
!
if(myrank/=0) allocate(Ns(NTK))
if(myrank/=0) allocate(Nb(NTK))
if(myrank/=0) allocate(Nt(NTK))
if(myrank/=0) allocate(UNT(Mb,NWF,NTK))
if(myrank/=0) allocate(C0_WN(NTG,NWF,NTK))
if(myrank/=0) allocate(H_MAT_R(NWF,NWF,-Na1:Na1,-Na2:Na2,-Na3:Na3))
call MPI_BARRIER(comm,ierr)
!
call MPI_Bcast(Ns,NTK,MPI_INTEGER,0,comm,ierr) 
call MPI_Bcast(Nb,NTK,MPI_INTEGER,0,comm,ierr) 
call MPI_Bcast(Nt,NTK,MPI_INTEGER,0,comm,ierr) 
call MPI_Bcast(UNT,Mb*NWF*NTK,MPI_DOUBLE_COMPLEX,0,comm,ierr) 
call MPI_Bcast(C0_WN,NTG*NWF*NTK,MPI_DOUBLE_COMPLEX,0,comm,ierr) 
call MPI_Bcast(H_MAT_R(1,1,-Na1,-Na2,-Na3),NWF*NWF*(2*Na1+1)*(2*Na2+1)*(2*Na3+1),MPI_DOUBLE_COMPLEX,0,comm,ierr) 
call MPI_BARRIER(comm,ierr)
!
!read eps 
!
if(myrank.eq.0)then 
 call rd_dat_chi_cutoff
 call rd_dat_wgrid 
 call rd_dat_sq 
 call rd_dat_eps!20180317  
 ! 
 !20180317 
 !
 do iq=1,Nq_irr 
  NG_for_eps=NGQ_eps(iq) 
  do ie=1,ne 
   do igL=1,NG_for_eps 
    epsirr(igL,igL,ie,iq)=epsirr(igL,igL,ie,iq)-1.0d0    
   enddo 
  enddo 
 enddo 
 !
 !
endif 
!
!chi_cutoff(300)  
!
call MPI_Bcast(Ecut_for_eps,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
!
!wgrid(135)  
!
call MPI_Bcast(Num_freq_grid,1,MPI_INTEGER,0,comm,ierr) 
call MPI_Bcast(ne,1,MPI_INTEGER,0,comm,ierr) 
!
!sq(301)  
!
call MPI_Bcast(Nq_irr,1,MPI_INTEGER,0,comm,ierr) 
call MPI_Bcast(NTQ,1,MPI_INTEGER,0,comm,ierr) 
call MPI_Bcast(NTGQ,1,MPI_INTEGER,0,comm,ierr) 
call MPI_Bcast(Lq1,1,MPI_INTEGER,0,comm,ierr) 
call MPI_Bcast(Lq2,1,MPI_INTEGER,0,comm,ierr) 
call MPI_Bcast(Lq3,1,MPI_INTEGER,0,comm,ierr) 
call MPI_BARRIER(comm,ierr)
!
if(myrank/=0) allocate(em(ne))
if(myrank/=0) allocate(pole_of_chi(ne))!20180317
if(myrank/=0) allocate(mat_b(ne,ne))!20180317
!
if(myrank/=0) allocate(SQI(3,Nq_irr)) 
if(myrank/=0) allocate(SQ(3,NTQ)) 
if(myrank/=0) allocate(numirrq(NTQ)) 
if(myrank/=0) allocate(numrotq(NTQ))  
if(myrank/=0) allocate(trsq(NTQ))  
if(myrank/=0) allocate(RWq(3,NTQ))
if(myrank/=0) allocate(LG0(3,NTG,NTQ))    
if(myrank/=0) allocate(NGQ_eps(NTQ))
if(myrank/=0) allocate(NGQ_psi(NTQ)) 
if(myrank/=0) allocate(packingq(-Lq1:Lq1,-Lq2:Lq2,-Lq3:Lq3,Nq_irr))
if(myrank/=0) allocate(epsirr(NTGQ,NTGQ,ne,Nq_irr))!real4 20180317 
call MPI_BARRIER(comm,ierr)
!
call MPI_Bcast(em,ne,MPI_DOUBLE_COMPLEX,0,comm,ierr) 
call MPI_Bcast(pole_of_chi,ne,MPI_DOUBLE_COMPLEX,0,comm,ierr)!20180317
call MPI_Bcast(mat_b,ne*ne,MPI_DOUBLE_COMPLEX,0,comm,ierr)!20180317
call MPI_Bcast(SQI,3*Nq_irr,MPI_DOUBLE_PRECISION,0,comm,ierr) 
call MPI_Bcast(SQ,3*NTQ,MPI_DOUBLE_PRECISION,0,comm,ierr) 
call MPI_Bcast(numirrq,NTQ,MPI_INTEGER,0,comm,ierr) 
call MPI_Bcast(numrotq,NTQ,MPI_INTEGER,0,comm,ierr) 
call MPI_Bcast(trsq,NTQ,MPI_INTEGER,0,comm,ierr) 
call MPI_Bcast(RWq,3*NTQ,MPI_INTEGER,0,comm,ierr)!bugfix: NTQ -> 3*NTQ !20180324  
call MPI_Bcast(LG0,3*NTG*NTQ,MPI_INTEGER,0,comm,ierr) 
call MPI_Bcast(NGQ_eps,NTQ,MPI_INTEGER,0,comm,ierr) 
call MPI_Bcast(NGQ_psi,NTQ,MPI_INTEGER,0,comm,ierr) 
call MPI_Bcast(packingq(-Lq1,-Lq2,-Lq3,1),(2*Lq1+1)*(2*Lq2+1)*(2*Lq3+1)*Nq_irr,MPI_INTEGER,0,comm,ierr) 
call MPI_Bcast(epsirr,NTGQ*NTGQ*ne*Nq_irr,MPI_COMPLEX,0,comm,ierr)!20180317 
call MPI_BARRIER(comm,ierr)
!
if(myrank.eq.0)then 
 write(6,*) 
 write(6,*)'================================'
 write(6,*)'dir-eps: FREQUENCY GRID IN EPSQW'
 write(6,*)'================================'
 write(6,*) 
 do ie=1,ne 
  write(6,'(2f15.10)') em(ie)
 enddo 
 write(6,*) 
 write(6,*)'==========================='
 write(6,*)'dir-eps: POLE OF MODEL GRID'
 write(6,*)'==========================='
 write(6,*) 
 do ie=1,ne 
  write(6,'(2f15.10)') pole_of_chi(ie)
 enddo 
endif  
!
!Rc, idlt, and Ncalc
!
avec_length(1)=dsqrt(a1(1)**2+a1(2)**2+a1(3)**2)
avec_length(2)=dsqrt(a2(1)**2+a2(2)**2+a2(3)**2)
avec_length(3)=dsqrt(a3(1)**2+a3(2)**2+a3(3)**2)
Rc=Rc_range_spacing*maxval(avec_length) 
!if(Mt<Ncalc) Mt=Ncalc!REMARK
Ncalc=20 
write(6,'(a25,i10)')'Considered bands for self-energy',Ncalc 
write(6,'(a25,f15.10)')'idlt (eV)',idlt*au 
write(6,'(a25,f15.10)')'Rc (AA)=',Rc*bohr   
!
!fbk 
!
allocate(fbk(NTB,NTK));fbk(:,:)=0.0d0 
do ik=1,NTK
 ikir=numirr(ik) 
 do ib=1,NTB
  FACTOR=(FermiEnergy-E_EIGI(ib,ikir))/idlt    
  fbk(ib,ik)=atan(FACTOR)/pi+0.5d0 
 enddo 
enddo 
SUM_REAL=0.0d0 
do ik=1,NTK
 ikir=numirr(ik) 
 do ib=1,NTB
  SUM_REAL=SUM_REAL+fbk(ib,ik) 
  en=(E_EIGI(ib,ikir)-FermiEnergy)*au  
 enddo
enddo
write(6,'(a10,f15.10)')'Nele=',2.0d0*SUM_REAL/dble(NTK)  
!
!OPEN(136,R,FILE='./dat.gwgrid') 
!
OPEN(136,FILE='./dat.gwgrid') 
rewind(136) 
read(136,*) nsgm 
allocate(sgmw(nsgm));sgmw(:)=0.0d0 
do ie=1,nsgm 
 read(136,'(2f15.10)') sgmw(ie) 
enddo 
write(6,*)'finish reading dat.gwgrid'
!
!MPI 
!
pnq=NTQ/nproc 
if(mod(NTQ,nproc).ne.0) then
 nbufq=pnq+1
else
 nbufq=pnq
end if
if(myrank.lt.mod(NTQ,nproc)) then
 pnq=pnq+1
 bnq=pnq*myrank+1 
 enq=bnq+pnq-1
else
 bnq=(pnq+1)*mod(NTQ,nproc)+pnq*(myrank-mod(NTQ,nproc))+1 
 enq=bnq+pnq-1
end if
write(file_id,'(a,3i7)')'bnq,enq,pnq',bnq,enq,pnq 
!
!fft
!
nfft1=nwx2+1; nfft2=nwy2+1; nfft3=nwz2+1; Nl123=nfft1*nfft2*nfft3 
call fft3_init(nwx2,nwy2,nwz2,nfft1,nfft2,nfft3,fs) 
!allocate(fftwk(Nl123*2),stat=err) 
!allocate(wfunc(Nl123*2),stat=err) 
!
!call MPI_FINALIZE(ierr)
!STOP 
!
!if(.true.) goto 9998 
!if(.true.) goto 9999 
!
allocate(SCirr(nsgm,Mb,Mb,Nk_irr)); SCirr(:,:,:,:)=0.0d0 
allocate(pSC(nsgm,Mb,Mb,Nk_irr)); pSC(:,:,:,:)=0.0d0 
!
do iq=1,pnq
 !
 q1=SQ(1,bnq+iq-1);q2=SQ(2,bnq+iq-1);q3=SQ(3,bnq+iq-1)
 !
 !q.ne.0 
 !
 if(q1/=0.0d0.or.q2/=0.0d0.or.q3/=0.0d0)then 
  NG_for_eps=NGQ_eps(bnq+iq-1)
  allocate(length_qg(NG_for_eps)); length_qg=0.0D0 
  allocate(atten_factor(NG_for_eps)); atten_factor=0.0D0 
  allocate(epsmk(NTGQ,NTGQ,ne)); epsmk=0.0d0 
  iqir=numirrq(bnq+iq-1)
  iop=numrotq(bnq+iq-1)
  !
  call make_eps(NTG,NTGQ,ne,trsq(bnq+iq-1),NGQ_eps(bnq+iq-1),LG0(1,1,bnq+iq-1),RWq(1,bnq+iq-1),&
  rginv(1,1,iop),pg(1,iop),nnp,Lq1,Lq2,Lq3,packingq(-Lq1,-Lq2,-Lq3,iqir),epsirr(1,1,1,iqir),&
  epsmk(1,1,1)) 
  !
  length_qg(:)=0.0d0 
  atten_factor(:)=0.0d0 
  do igL=1,NG_for_eps   
   igL1=LG0(1,igL,bnq+iq-1)
   igL2=LG0(2,igL,bnq+iq-1)
   igL3=LG0(3,igL,bnq+iq-1)
   qgL(:)=(q1+dble(igL1))*b1(:)+(q2+dble(igL2))*b2(:)+(q3+dble(igL3))*b3(:)
   qgL2=qgL(1)**2+qgL(2)**2+qgL(3)**2
   qgL1=dsqrt(qgL2) 
   length_qg(igL)=qgL1
   atten_factor(igL)=dsqrt(1.0d0-dcos(qgL1*Rc))   
  enddo!igL 
!---
!$OMP PARALLEL PRIVATE(ib,ik,jb,shift_G,ikq,C0_K,C0_KmQ,rho,wfunc,fftwk,rho_tmp,ikir,iop,kb,ie,SUM_CMPX,&
!$OMP&         ikqir,ikqop,igL,jgL,vecf,je,veca,delta,de,sgn,dnm,pSComp) 
  allocate(rho(NG_for_eps,Mb,Nk_irr));rho=0.0d0
  allocate(fftwk(Nl123*2));fftwk=0.0d0 
  allocate(wfunc(Nl123*2));wfunc=0.0d0 
  allocate(rho_tmp(NG_for_eps));rho_tmp=0.0d0
  allocate(pSComp(nsgm,Mb,Mb,Nk_irr)); pSComp=0.0d0  
  allocate(C0_K(NTG));C0_K=0.0d0
  allocate(C0_KmQ(NTG));C0_KmQ=0.0d0
  allocate(vecf(ne));vecf=0.0d0
  allocate(veca(ne));veca=0.0d0
!$OMP DO  
  do ib=1,Ncalc
   if(myrank.eq.0) write(6,*)'ib=',ib
   !
   rho(:,:,:)=0.0D0 
   !
   do ikir=1,Nk_irr
    !
    ik=numMK(ikir) 
    iop=numrot(ik) 
    !
    do jb=1,Mb  
     !
     !C0_K(:)=CIR(:,jb+Ns(ik),ik)  
     !
     call make_C0_for_given_band(NTG,trs(ik),NG0(ik),KG0(1,1,ik),RW(1,ik),rginv(1,1,iop),pg(1,iop),&
     L1,L2,L3,packing(-L1,-L2,-L3,ikir),CIR(1,jb+Ns(ik),ikir),C0_K(1)) 
     ! 
     shift_G(:)=0
     call search_kq(NTK,SK0(1,1),-q1,-q2,-q3,ik,ikq,shift_G(1))
     shift_G(:)=-shift_G(:)
     !
     ikqir=numirr(ikq)
     ikqop=numrot(ikq) 
     !
     !call make_C0(NTG,trs(ikq),NG0(ikq),KG0(1,1,ikq),RW(1,ikq),rginv(1,1,iop),pg(1,iop),L1,L2,L3,&
     !packing(-L1,-L2,-L3,ikir),CIR(1,ib,ikir),C0_KmQ(1)) 
     !
     call make_C0_for_given_band(NTG,trs(ikq),NG0(ikq),KG0(1,1,ikq),RW(1,ikq),rginv(1,1,ikqop),pg(1,ikqop),&
     L1,L2,L3,packing(-L1,-L2,-L3,ikqir),CIR(1,ib,ikqir),C0_KmQ(1)) 
     !
     call calc_InterStateMatrix(NTK,NTG,NG0(1),KG0(1,1,1),C0_KmQ(1),C0_K(1),ikq,ik,nwx2,nwy2,nwz2,&
     nfft1,nfft2,Nl123,wfunc(1),fftwk(1),fs,LG0(1,1,bnq+iq-1),NG_for_eps,shift_G(1),rho_tmp(1))
     !
     rho(:,jb,ikir)=rho_tmp(:)
     !
    enddo!jb  
   enddo!ikir 
   !
   do ikir=1,Nk_irr 
    !
    ik=numMK(ikir) 
    !
    shift_G(:)=0
    call search_kq(NTK,SK0(1,1),-q1,-q2,-q3,ik,ikq,shift_G(1))
    !
    do jb=1,Mb 
     do kb=1,Mb 
      vecf=0.0d0
      do ie=1,ne 
       SUM_CMPX=0.0d0
       do igL=1,NG_for_eps  
        do jgL=1,NG_for_eps 
         !
         !SUM_CMPX=SUM_CMPX+rho(igL,jb,ik)/length_qg(igL)*atten_factor(igL)*epsmk(igL,jgL,ie)*CONJG(rho(jgL,kb,ik))/length_qg(jgL)*atten_factor(jgL) 
         !
         SUM_CMPX=SUM_CMPX+rho(igL,jb,ikir)/length_qg(igL)*atten_factor(igL)*epsmk(igL,jgL,ie)&
                 *CONJG(rho(jgL,kb,ikir))/length_qg(jgL)*atten_factor(jgL) 
         !
        enddo!jgL 
       enddo!igL 
       vecf(ie)=2.0d0*tpi*SUM_CMPX/dble(NTQ)/VOLUME 
      enddo!ie 
      !
      veca=0.0d0 
      do je=1,ne 
       SUM_CMPX=0.0d0 
       do ie=1,ne 
        SUM_CMPX=SUM_CMPX+mat_b(je,ie)*vecf(ie) 
       enddo!ie 
       veca(je)=SUM_CMPX
      enddo!je 
      !
      do ie=1,nsgm 
       delta=sgmw(ie) 
       !
       !ikir=numirr(ikq) 
       !de=delta-E_EIGI(ib,ikir)
       !
       ikqir=numirr(ikq) 
       de=delta-E_EIGI(ib,ikqir)
       !
       !if(E_EIGI(ib,ikir)>FermiEnergy)then 
       !
       if(E_EIGI(ib,ikqir)>FermiEnergy)then 
        sgn=1.0d0
       else
        sgn=-1.0d0
       endif 
       !
       do je=1,ne 
        !
        dnm=de-(pole_of_chi(je)-ci*idlt)*sgn
        !
        !pSComp(ie,jb,kb,ik)=pSComp(ie,jb,kb,ik)+veca(je)/dnm 
        !
        pSComp(ie,jb,kb,ikir)=pSComp(ie,jb,kb,ikir)+veca(je)/dnm 
        !
       enddo!je 
      enddo!ie  
     enddo!kb 
    enddo!jb 
   enddo!ikir  
  enddo!ib 
!$OMP END DO
!$OMP CRITICAL
  pSC=pSC+pSComp
!$OMP END CRITICAL
  deallocate(fftwk,wfunc,rho_tmp,rho,pSComp,C0_K,C0_KmQ,vecf,veca)
!$OMP END PARALLEL
!
 if(myrank.eq.0) write(6,*)'FINISH iq',iq 
 deallocate(length_qg,atten_factor,epsmk) 
 !
 !q.eq.0 
 !
 elseif(q1==0.0d0.and.q2==0.0d0.and.q3==0.0d0) then 
  NG_for_eps=NGQ_eps(bnq+iq-1)
  allocate(length_qg(NG_for_eps));length_qg=0.0d0 
  allocate(atten_factor(NG_for_eps));atten_factor=0.0d0  
  allocate(rho(NG_for_eps,Mb,Nk_irr));rho=0.0d0 
  allocate(epsmk(NTGQ,NTGQ,ne));epsmk(:,:,:)=0.0d0 
  iqir=numirrq(bnq+iq-1)
  iop=numrotq(bnq+iq-1)
  !
  call make_eps(NTG,NTGQ,ne,trsq(bnq+iq-1),NGQ_eps(bnq+iq-1),LG0(1,1,bnq+iq-1),RWq(1,bnq+iq-1),&
  rginv(1,1,iop),pg(1,iop),nnp,Lq1,Lq2,Lq3,packingq(-Lq1,-Lq2,-Lq3,iqir),epsirr(1,1,1,iqir),&
  epsmk(1,1,1)) 
  !
  length_qg(:)=0.0D0 
  do igL=1,NG_for_eps   
   igL1=LG0(1,igL,bnq+iq-1)
   igL2=LG0(2,igL,bnq+iq-1)
   igL3=LG0(3,igL,bnq+iq-1)
   if(igL1==0.and.igL2==0.and.igL3==0)then 
    No_G_0=igL  
   else
    qgL(:)=(q1+dble(igL1))*b1(:)+(q2+dble(igL2))*b2(:)+(q3+dble(igL3))*b3(:)
    qgL2=qgL(1)**2+qgL(2)**2+qgL(3)**2
    qgL1=dsqrt(qgL2) 
    length_qg(igL)=qgL1
    atten_factor(igL)=dsqrt(1.0d0-dcos(qgL1*Rc))   
   endif 
  enddo!igL 
  !
  do ib=1,Ncalc
   if(myrank.eq.0) write(6,*)'ib=',ib
   rho(:,:,:)=0.0D0 
!$OMP PARALLEL PRIVATE(ik,jb,ikqir,ikqop,shift_G,ikq,C0_K,C0_KmQ,wfunc,fftwk,rho_tmp,ikir,iop) 
   allocate(fftwk(Nl123*2));fftwk=0.0d0 
   allocate(wfunc(Nl123*2));wfunc=0.0d0  
   allocate(rho_tmp(NG_for_eps));rho_tmp=0.0d0 
   allocate(C0_K(NTG));C0_K=0.0d0
   allocate(C0_KmQ(NTG));C0_KmQ=0.0d0
!$OMP DO
   do ikir=1,Nk_irr 
    !
    ik=numMK(ikir) 
    iop=numrot(ik) 
    !
    do jb=1,Mb 
     !
     !C0_K(:)=CIR(:,jb+Ns(ik),ik)  
     !
     call make_C0_for_given_band(NTG,trs(ik),NG0(ik),KG0(1,1,ik),RW(1,ik),rginv(1,1,iop),pg(1,iop),&
     L1,L2,L3,packing(-L1,-L2,-L3,ikir),CIR(1,jb+Ns(ik),ikir),C0_K(1)) 
     !
     shift_G(:)=0
     ikq=ik 
     !
     !ikir=numirr(ikq)
     !iop=numrot(ikq) 
     !
     ikqir=numirr(ikq)
     ikqop=numrot(ikq) 
     !
     !call make_C0(NTG,trs(ikq),NG0(ikq),KG0(1,1,ikq),RW(1,ikq),rginv(1,1,iop),pg(1,iop),L1,L2,L3,packing(-L1,-L2,-L3,ikir),&
     !CIR(1,ib,ikir),C0_KmQ(1)) 
     !
     call make_C0_for_given_band(NTG,trs(ikq),NG0(ikq),KG0(1,1,ikq),RW(1,ikq),rginv(1,1,ikqop),pg(1,ikqop),&
     L1,L2,L3,packing(-L1,-L2,-L3,ikqir),CIR(1,ib,ikqir),C0_KmQ(1)) 
     !
     call calc_InterStateMatrix(NTK,NTG,NG0(1),KG0(1,1,1),C0_KmQ(1),C0_K(1),ikq,ik,nwx2,nwy2,nwz2,&
     nfft1,nfft2,Nl123,wfunc(1),fftwk(1),fs,LG0(1,1,bnq+iq-1),NG_for_eps,shift_G(1),rho_tmp(1))
     !
     !rho(:,jb,ik)=rho_tmp(:)
     !
     rho(:,jb,ikir)=rho_tmp(:)
     !
    enddo!jb
   enddo!ikir 
!$OMP END DO
   deallocate(fftwk,wfunc,rho_tmp,C0_K,C0_KmQ)
!$OMP END PARALLEL
!---
!$OMP PARALLEL PRIVATE(ik,jb,kb,ie,SUM_CMPX,igL,jgL,vecf,je,veca,ikq,delta,ikir,ikqir,de,sgn,dnm) 
  allocate(vecf(ne));vecf=0.0d0
  allocate(veca(ne));veca=0.0d0
!$OMP DO
   do ikir=1,Nk_irr 
    !
    ik=numMK(ikir)
    !
    do jb=1,Mb 
     do kb=1,Mb 
      vecf=0.0d0 
      do ie=1,ne 
       SUM_CMPX=0.0d0
       do igL=1,NG_for_eps 
        if(igL==No_G_0) cycle 
        do jgL=1,NG_for_eps 
         if(jgL==No_G_0) cycle 
         !
         !SUM_CMPX=SUM_CMPX+rho(igL,jb,ik)/length_qg(igL)*atten_factor(igL)*epsmk(igL,jgL,ie)*CONJG(rho(jgL,kb,ik))/length_qg(jgL)*atten_factor(jgL) 
         !
         SUM_CMPX=SUM_CMPX+rho(igL,jb,ikir)/length_qg(igL)*atten_factor(igL)*epsmk(igL,jgL,ie)&
                 *CONJG(rho(jgL,kb,ikir))/length_qg(jgL)*atten_factor(jgL) 
         !
        enddo!jgL 
       enddo!igL  
       vecf(ie)=2.0d0*tpi*SUM_CMPX/VOLUME/dble(NTQ)  
      enddo!ie 
      veca=0.0d0
      do je=1,ne 
       SUM_CMPX=0.0d0 
       do ie=1,ne 
        SUM_CMPX=SUM_CMPX+mat_b(je,ie)*vecf(ie) 
       enddo!ie  
       veca(je)=SUM_CMPX 
      enddo!je  
      do ie=1,nsgm 
       ikq=ik 
       delta=sgmw(ie) 
       !
       !ikir=numirr(ikq) 
       !de=delta-E_EIGI(ib,ikir)
       !
       ikqir=numirr(ikq) 
       de=delta-E_EIGI(ib,ikqir)
       !
       !if(E_EIGI(ib,ikir)>FermiEnergy)then 
       !
       if(E_EIGI(ib,ikqir)>FermiEnergy)then 
        sgn=1.0d0
       else
        sgn=-1.0d0
       endif 
       !
       do je=1,ne 
        dnm=de-(pole_of_chi(je)-ci*idlt)*sgn
        !
        !pSC(ie,jb,kb,ik)=pSC(ie,jb,kb,ik)+veca(je)/dnm 
        !
        pSC(ie,jb,kb,ikir)=pSC(ie,jb,kb,ikir)+veca(je)/dnm 
        !
       enddo!je  
      enddo!ie  
     enddo!kb 
    enddo!jb 
   enddo!ikir 
!$OMP END DO
  deallocate(vecf,veca)
!$OMP END PARALLEL
  enddo!ib 
  deallocate(length_qg,atten_factor,rho,epsmk) 
  !
  !<head contribution> 
  !
  !(i) atten(Spencer Alabi)
  !
  chead=(tpi/dble(NTQ)/VOLUME)*Rc*Rc   
  !
  !(ii) Louie
  !
  !qsz=(6.0D0*(pi**2)/dble(NTQ)/dble(VOLUME))**(1.0D0/3.0D0)  
  !chead=(2.0D0/pi)*qsz 
  !
  write(file_id,*)'bnq+iq-1=',bnq+iq-1 
  write(file_id,*)'correction_head=',chead   
  !
  !cGW
  !
  !chead=0.0d0    
  !write(file_id,*)'cGW, then correction_head=0'
  !
  write(file_id,*)'No_G_0=',No_G_0 
  allocate(vecf(ne));vecf=0.0d0
  allocate(veca(ne));veca=0.0d0
  do ikir=1,Nk_irr 
   !
   !ikir=numirr(ik) 
   !
   ik=numMK(ikir) 
   !
   do jb=1,Mb 
    !
    vecf=0.0d0 
    do ie=1,ne 
     vecf(ie)=epsirr(No_G_0,No_G_0,ie,bnq+iq-1)
    enddo!ie   
    !
    veca=0.0d0 
    do je=1,ne 
     SUM_CMPX=0.0d0 
     do ie=1,ne 
      SUM_CMPX=SUM_CMPX+mat_b(je,ie)*vecf(ie) 
     enddo!ie 
     veca(je)=SUM_CMPX 
    enddo!je  
    !
    do ie=1,nsgm  
     delta=sgmw(ie) 
     !
     !de=delta-E_EIGI(jb+Ns(ikir),ikir) 
     !
     de=delta-E_EIGI(jb+Ns(ik),ikir) 
     !
     !if(E_EIGI(jb+Ns(ikir),ikir)>FermiEnergy)then 
     !
     if(E_EIGI(jb+Ns(ik),ikir)>FermiEnergy)then 
      sgn=1.0d0
     else
      sgn=-1.0d0
     endif 
     do je=1,ne 
      dnm=de-(pole_of_chi(je)-ci*idlt)*sgn
      !
      !pSC(ie,jb,jb,ik)=pSC(ie,jb,jb,ik)+veca(je)/dnm*chead  
      !
      pSC(ie,jb,jb,ikir)=pSC(ie,jb,jb,ikir)+veca(je)/dnm*chead  
      !
     enddo!je
    enddo!ie  
   enddo!jb 
  enddo!ikir 
  deallocate(vecf,veca) 
  if(myrank.eq.0) write(6,*)'FINISH iq',iq 
 endif 
enddo!iq 
!
call MPI_BARRIER(comm,ierr)
write(file_id,*) 'I finished pSC calc' 
call MPI_REDUCE(pSC,SCirr,nsgm*Mb*Mb*Nk_irr,MPI_DOUBLE_COMPLEX,MPI_SUM,0,comm,ierr)
!
if(myrank.eq.0)then 
 write(1500) SCirr  
 end_time=MPI_Wtime()
 diff_time=end_time-start_time 
 write(6,*)'#TOTAL TIME=',diff_time 
endif 
!
call MPI_FINALIZE(ierr)
STOP 
!
9998 write(6,*)'calc SCR' 
if(myrank.eq.0)then 
 allocate(SCirr(nsgm,Mb,Mb,Nk_irr)) 
 read(1500) SCirr 
 !
 allocate(SC(Mb,Mb,NTK,nsgm)); SC(:,:,:,:)=0.0d0
 do ik=1,NTK 
  ikir=numirr(ik) 
  if(trs(ik)==1) then 
   do ie=1,nsgm
    do ib=1,Mb!Nb(ik)
     do jb=1,Mb!Nb(ik) 
      SC(ib,jb,ik,ie)=SCirr(ie,ib,jb,ikir) 
     enddo 
    enddo 
   enddo 
  elseif(trs(ik)==-1) then 
   do ie=1,nsgm  
    do ib=1,Mb!Nb(ik)
     do jb=1,Mb!Nb(ik) 
      SC(ib,jb,ik,ie)=SCirr(ie,jb,ib,ikir) 
     enddo 
    enddo 
   enddo 
  endif 
 enddo 
 deallocate(SCirr) 
 !---
 allocate(pf(-Na1:Na1,-Na2:Na2,-Na3:Na3,NTK));pf=0.0d0     
 do ik=1,NTK 
  do ia3=-Na3,Na3 
   do ia2=-Na2,Na2 
    do ia1=-Na1,Na1 
     phase=tpi*(SK0(1,ik)*dble(ia1)+SK0(2,ik)*dble(ia2)+SK0(3,ik)*dble(ia3)) 
     pf(ia1,ia2,ia3,ik)=exp(-ci*phase) 
    enddo!ia1  
   enddo!ia2  
  enddo!ia3  
 enddo!ik  
 write(6,*)'#finish make pf'
 !
 !
 !OPEN(157,W,FILE='SCR') 
 !
 rewind(157)
 allocate(SCR(nsgm,NWF,NWF,-Na1:Na1,-Na2:Na2,-Na3:Na3));SCR(:,:,:,:,:,:)=0.0d0 
 do ia3=-Na3,Na3
  do ia2=-Na2,Na2
   do ia1=-Na1,Na1
    do iw=1,NWF
     do jw=1,NWF 
      do ie=1,nsgm  
       SUM_CMPX=0.0D0 
!$OMP PARALLEL reduction(+:SUM_CMPX) private(psum,ik,jb,kb)  
       psum=0.0d0 
!$OMP DO 
       do ik=1,NTK 
        do jb=1,Nb(ik) 
         do kb=1,Nb(ik) 
          psum=psum+CONJG(UNT(jb,iw,ik))*SC(jb,kb,ik,ie)*UNT(kb,jw,ik)*pf(ia1,ia2,ia3,ik) 
         enddo!kb 
        enddo!jb 
       enddo!ik 
!$OMP END DO 
       SUM_CMPX=SUM_CMPX+psum 
!$OMP END PARALLEL 
       SCR(ie,iw,jw,ia1,ia2,ia3)=SUM_CMPX/DBLE(NTK)
      enddo!ie
     enddo!jw
    enddo!iw
    write(6,'(a20,x,3i7)')"finish ia1,ia2,ia3",ia1,ia2,ia3 
   enddo!ia1
  enddo!ia2 
 enddo!ia3  
 !
 !do ia1=-Na1,Na1
 ! do ia2=-Na2,Na2 
 !  do ia3=-Na3,Na3 
 !   do jw=1,NWF
 !    do iw=1,NWF 
 !     write(157)(SCR(ie,iw,jw,ia1,ia2,ia3),ie=1,nsgm) 
 !    enddo!iw  
 !   enddo!jw  
 !  enddo!ia3 
 ! enddo!ia2 
 !enddo!ia1 
 !close(157) 
 !
 write(5004) SCR 
 write(6,*)'finish SCR' 
 !
endif!myrank.eq.0 
call MPI_FINALIZE(ierr)
STOP
!
9999 write(6,*)'calc SX' 
allocate(SXirr(Mb,Mb,Nk_irr));SXirr(:,:,:)=0.0d0 
allocate(pSX(Mb,Mb,Nk_irr));pSX(:,:,:)=0.0d0 
do iq=1,pnq 
 !
 q1=SQ(1,bnq+iq-1); q2=SQ(2,bnq+iq-1); q3=SQ(3,bnq+iq-1)
 !
 !q.ne.0 
 !
 if(q1/=0.0d0.or.q2/=0.0d0.or.q3/=0.0d0)then 
  NG_for_psi=NGQ_psi(bnq+iq-1)
  allocate(atten_factor(NG_for_psi));atten_factor(:)=0.0D0 
  allocate(length_qg(NG_for_psi));length_qg(:)=0.0D0 
  allocate(rho(NG_for_psi,Mb,Nk_irr))
  do igL=1,NG_for_psi   
   igL1=LG0(1,igL,bnq+iq-1)
   igL2=LG0(2,igL,bnq+iq-1)
   igL3=LG0(3,igL,bnq+iq-1)
   qgL(:)=(q1+dble(igL1))*b1(:)+(q2+dble(igL2))*b2(:)+(q3+dble(igL3))*b3(:)
   qgL2=qgL(1)**2+qgL(2)**2+qgL(3)**2
   qgL1=dsqrt(qgL2) 
   length_qg(igL)=qgL1
   atten_factor(igL)=1.0d0-dcos(qgL1*Rc)  
  enddo 
  do ib=1,Ncalc!Mt!
   rho(:,:,:)=0.0D0 
!$OMP PARALLEL PRIVATE(ik,shift_G,ikq,C0_K,C0_KmQ,jb,wfunc,fftwk,rho_tmp,ikir,iop,ikqir,ikqop) 
   allocate(fftwk(Nl123*2));fftwk=0.0d0 
   allocate(wfunc(Nl123*2));wfunc=0.0d0 
   allocate(rho_tmp(NG_for_psi));rho_tmp=0.0d0 
   allocate(C0_K(NTG));C0_K=0.0d0
   allocate(C0_KmQ(NTG));C0_KmQ=0.0d0
!$OMP DO
   do ikir=1,Nk_irr 
    !
    ik=numMK(ikir) 
    iop=numrot(ik) 
    !
    do jb=1,Mb 
     !
     !C0_K(:)=CIR(:,jb+Ns(ik),ik)  
     !
     call make_C0_for_given_band(NTG,trs(ik),NG0(ik),KG0(1,1,ik),RW(1,ik),rginv(1,1,iop),pg(1,iop),&
     L1,L2,L3,packing(-L1,-L2,-L3,ikir),CIR(1,jb+Ns(ik),ikir),C0_K(1)) 
     ! 
     shift_G(:)=0
     call search_kq(NTK,SK0(1,1),-q1,-q2,-q3,ik,ikq,shift_G(1))
     shift_G(:)=-shift_G(:)  
     !
     ikqir=numirr(ikq)
     ikqop=numrot(ikq) 
     !
     !call make_C0_for_given_band(NTG,trs(ikq),NG0(ikq),KG0(1,1,ikq),RW(1,ikq),rginv(1,1,iop),pg(1,iop),&
     !L1,L2,L3,packing(-L1,-L2,-L3,ikir),CIR(1,ib,ikir),C0_KmQ(1)) 
     !
     call make_C0_for_given_band(NTG,trs(ikq),NG0(ikq),KG0(1,1,ikq),RW(1,ikq),rginv(1,1,ikqop),pg(1,ikqop),&
     L1,L2,L3,packing(-L1,-L2,-L3,ikqir),CIR(1,ib,ikqir),C0_KmQ(1)) 
     !
     !
     call calc_InterStateMatrix(NTK,NTG,NG0(1),KG0(1,1,1),C0_KmQ(1),C0_K(1),ikq,ik,nwx2,nwy2,nwz2,&
     nfft1,nfft2,Nl123,wfunc(1),fftwk(1),fs,LG0(1,1,bnq+iq-1),NG_for_psi,shift_G(1),rho_tmp(1))
     !
     !
     !rho(:,jb,ik)=rho_tmp(:)
     !
     rho(:,jb,ikir)=rho_tmp(:)
     !
    enddo!jb  
   enddo!ik 
!$OMP END DO 
   deallocate(fftwk,wfunc,rho_tmp,C0_K,C0_KmQ)
!$OMP END PARALLEL 
!$OMP PARALLEL PRIVATE(ik,ikir,shift_G,ikq,jb,kb,SUM_CMPX,igL)  
!$OMP DO
   do ikir=1,Nk_irr 
    !
    ik=numMK(ikir) 
    ! 
    shift_G(:)=0
    call search_kq(NTK,SK0(1,1),-q1,-q2,-q3,ik,ikq,shift_G(1))
    !
    do jb=1,Mb 
     do kb=1,Mb 
      SUM_CMPX=0.0d0
      do igL=1,NG_for_psi 
       !
       !SUM_CMPX=SUM_CMPX+fbk(ib,ikq)*rho(igL,jb,ik)*CONJG(rho(igL,kb,ik))/((length_qg(igL))**2)*atten_factor(igL)
       !
       SUM_CMPX=SUM_CMPX+fbk(ib,ikq)*rho(igL,jb,ikir)&
               *CONJG(rho(igL,kb,ikir))/((length_qg(igL))**2)*atten_factor(igL)
       !
      enddo 
      !
      !pSX(jb,kb,ik)=pSX(jb,kb,ik)+SUM_CMPX 
      !
      pSX(jb,kb,ikir)=pSX(jb,kb,ikir)+SUM_CMPX 
      !
     enddo!kb 
    enddo!jb 
   enddo!ik 
!$OMP END DO 
!$OMP END PARALLEL 
  enddo!ib 
  write(file_id,*)'FINISH iq',iq
  if(myrank.eq.0) write(6,*)'FINISH iq',iq 
  deallocate(length_qg,atten_factor,rho) 
 elseif(q1==0.0d0.and.q2==0.0d0.and.q3==0.0d0)then 
  !q.eq.0 
  write(file_id,'(a,x,3f10.5)')'q1,q2,q3=',q1,q2,q3  
  NG_for_psi=NGQ_psi(bnq+iq-1)
  allocate(length_qg(NG_for_psi));length_qg(:)=0.0D0 
  allocate(atten_factor(NG_for_psi));atten_factor(:)=0.0D0 
  allocate(rho(NG_for_psi,Mb,Nk_irr))
  do igL=1,NG_for_psi 
   igL1=LG0(1,igL,bnq+iq-1)
   igL2=LG0(2,igL,bnq+iq-1)
   igL3=LG0(3,igL,bnq+iq-1)
   if(igL1==0.and.igL2==0.and.igL3==0)then 
    No_G_0=igL  
   else
    qgL(:)=(q1+dble(igL1))*b1(:)+(q2+dble(igL2))*b2(:)+(q3+dble(igL3))*b3(:)
    qgL2=qgL(1)**2+qgL(2)**2+qgL(3)**2
    qgL1=dsqrt(qgL2) 
    length_qg(igL)=qgL1
    atten_factor(igL)=1.0d0-dcos(qgL1*Rc)  
   endif 
  enddo!igL 
  do ib=1,Ncalc!Mt 
   rho(:,:,:)=0.0D0 
!$OMP PARALLEL PRIVATE(ik,shift_G,ikq,C0_K,C0_KmQ,jb,wfunc,fftwk,rho_tmp,ikir,iop,ikqir,ikqop) 
   allocate(fftwk(Nl123*2));fftwk=0.0d0 
   allocate(wfunc(Nl123*2));wfunc=0.0d0 
   allocate(rho_tmp(NG_for_psi));rho_tmp=0.0d0
   allocate(C0_K(NTG));C0_K=0.0d0
   allocate(C0_KmQ(NTG));C0_KmQ=0.0d0
!$OMP DO
   do ikir=1,Nk_irr 
    !
    ik=numMK(ikir) 
    iop=numrot(ik) 
    !
    do jb=1,Mb 
     !
     !C0_K(:)=CIR(:,jb+Ns(ik),ik)  
     !
     call make_C0_for_given_band(NTG,trs(ik),NG0(ik),KG0(1,1,ik),RW(1,ik),rginv(1,1,iop),pg(1,iop),&
     L1,L2,L3,packing(-L1,-L2,-L3,ikir),CIR(1,jb+Ns(ik),ikir),C0_K(1)) 
     !
     shift_G(:)=0
     ikq=ik 
     !
     ikqir=numirr(ikq)
     ikqop=numrot(ikq) 
     !
     !call make_C0_for_given_band(NTG,trs(ikq),NG0(ikq),KG0(1,1,ikq),RW(1,ikq),rginv(1,1,iop),pg(1,iop),&
     !L1,L2,L3,packing(-L1,-L2,-L3,ikir),CIR(1,ib,ikir),C0_KmQ(1)) 
     !
     call make_C0_for_given_band(NTG,trs(ikq),NG0(ikq),KG0(1,1,ikq),RW(1,ikq),rginv(1,1,ikqop),pg(1,ikqop),&
     L1,L2,L3,packing(-L1,-L2,-L3,ikqir),CIR(1,ib,ikqir),C0_KmQ(1)) 
     !
     call calc_InterStateMatrix(NTK,NTG,NG0(1),KG0(1,1,1),C0_KmQ(1),C0_K(1),ikq,ik,nwx2,nwy2,nwz2,&
     nfft1,nfft2,Nl123,wfunc(1),fftwk(1),fs,LG0(1,1,bnq+iq-1),NG_for_psi,shift_G(1),rho_tmp(1))
     !
     !rho(:,jb,ik)=rho_tmp(:)
     !
     rho(:,jb,ikir)=rho_tmp(:)
     !
    enddo!jb 
   enddo!ik  
!$OMP END DO
   deallocate(fftwk,wfunc,rho_tmp,C0_K,C0_KmQ)
!$OMP END PARALLEL 
!$OMP PARALLEL PRIVATE(ikir,ik,ikq,jb,kb,SUM_CMPX,igL) 
!$OMP DO
   do ikir=1,Nk_irr 
    !
    ik=numMK(ikir)
    !
    ikq=ik 
    !
    do jb=1,Mb 
     do kb=1,Mb 
      SUM_CMPX=0.0d0
      do igL=1,NG_for_psi 
       if(igL.ne.No_G_0) then 
        !
        !SUM_CMPX=SUM_CMPX+fbk(ib,ikq)*rho(igL,jb,ik)*CONJG(rho(igL,kb,ik))/((length_qg(igL))**2)*atten_factor(igL)
        !
        SUM_CMPX=SUM_CMPX+fbk(ib,ikq)*rho(igL,jb,ikir)*CONJG(rho(igL,kb,ikir))/((length_qg(igL))**2)*atten_factor(igL)
        !
       endif 
      enddo!igL 
      !
      !pSX(jb,kb,ik)=pSX(jb,kb,ik)+SUM_CMPX 
      !
      pSX(jb,kb,ikir)=pSX(jb,kb,ikir)+SUM_CMPX 
      !
     enddo!kb 
    enddo!jb 
   enddo!ik 
!$OMP END DO
!$OMP END PARALLEL 
  enddo ! ib 
  write(file_id,*)'FINISH iq',iq
  if(myrank.eq.0) write(6,*)'FINISH iq',iq 
  deallocate(length_qg,atten_factor,rho) 
 endif 
enddo!iq 
!---
call MPI_BARRIER(comm,ierr)
write(file_id,*)'I finished pSX calc' 
call MPI_REDUCE(pSX,SXirr,Mb*Mb*Nk_irr,MPI_DOUBLE_COMPLEX,MPI_SUM,0,comm,ierr)
if(myrank.eq.0)then 
 SXirr(:,:,:)=2.0d0*tpi*SXirr(:,:,:)/dble(NTQ)/VOLUME 
 !--
 !head contribution iq=NTQ
 !--
 write(6,*)'Spencer-Alabi'
 chead=(tpi/dble(NTQ)/VOLUME)*Rc*Rc   
 write(6,*)'correction_head=',chead   
 !--
 !write(6,*)'Hybertsen-Louie'
 !qsz=(6.0D0*(pi**2)/dble(NTQ)/dble(VOLUME))**(1.0D0/3.0D0)  
 !chead=(2.0D0/pi)*qsz 
 !write(6,*)'correction_head=',chead   
 !---
 allocate(SX_DIV(Mb,Mb,Nk_irr));SX_DIV(:,:,:)=0.0d0 
 do ikir=1,Nk_irr 
  !
  ik=numMK(ikir)
  ! 
  do jb=1,Mb 
   !
   !SX_DIV(jb,jb,ik)=fbk(jb+Ns(ik),ik)*chead 
   !
   SX_DIV(jb,jb,ikir)=fbk(jb+Ns(ik),ik)*chead 
   !
  enddo!jb 
 enddo!ik 
 !--
 SXirr=SXirr+SX_DIV!fGW
 !--
 !SXirr=SXirr      !cGW
 !write(1600) SXirr
 !deallocate(SXirr) 
 !end_time=MPI_Wtime()
 !diff_time=end_time-start_time 
 !write(6,*)'#TOTAL TIME=',diff_time 
 !endif!myrank.eq.0
 !call MPI_FINALIZE(ierr)
 !STOP 
 !9997 write(6,*)'skip calc sx'
 !if(myrank.eq.0)then 
 !allocate(SXirr(Mb,Mb,Nk_irr))
 !read(1600) SXirr 
 !--
 allocate(SX(Mb,Mb,NTK));SX(:,:,:)=0.0d0 
 do ik=1,NTK 
  ikir=numirr(ik) 
  if(trs(ik)==1) then 
   SX(:,:,ik)=SXirr(:,:,ikir) 
  elseif(trs(ik)==-1) then 
   do ib=1,Nb(ik)
    do jb=1,Nb(ik) 
     SX(ib,jb,ik)=SXirr(jb,ib,ikir) 
    enddo!jb
   enddo!ib 
  endif 
 enddo!ik 
 !
 !OPEN(154,FILE='./dat.energy_vs_sx') 
 !
 OPEN(154,FILE='./dat.energy_vs_sx') 
 rewind(154) 
 do ik=1,Nk_irr 
  do jb=1,Mb 
   en=(E_EIGI(jb+Ns(ik),ik)-FermiEnergy)*au 
   write(154,'(3F20.10)') en,SXirr(jb,jb,ik)*au 
  enddo!jb 
 enddo!ik 
 deallocate(SXirr) 
 !---
 allocate(pf(-Na1:Na1,-Na2:Na2,-Na3:Na3,NTK));pf=0.0d0     
 do ik=1,NTK 
  do ia3=-Na3,Na3 
   do ia2=-Na2,Na2 
    do ia1=-Na1,Na1 
     phase=tpi*(SK0(1,ik)*dble(ia1)+SK0(2,ik)*dble(ia2)+SK0(3,ik)*dble(ia3)) 
     pf(ia1,ia2,ia3,ik)=exp(-ci*phase) 
    enddo  
   enddo  
  enddo  
 enddo  
 write(6,*)'#finish make pf'
!---
 allocate(MAT_SX_R(NWF,NWF,-Na1:Na1,-Na2:Na2,-Na3:Na3));MAT_SX_R(:,:,:,:,:)=0.0d0 
 do ia3=-Na3,Na3
  do ia2=-Na2,Na2
   do ia1=-Na1,Na1
    do iw=1,NWF
     do jw=1,NWF 
      SUM_CMPX=0.0D0 
      if(.true.)then!.true.=off-diag, .false.=diag
       do ik=1,NTK 
        do ib=1,Nb(ik) 
         do jb=1,Nb(ik) 
          SUM_CMPX=SUM_CMPX+CONJG(UNT(ib,iw,ik))*SX(ib,jb,ik)*UNT(jb,jw,ik)*pf(ia1,ia2,ia3,ik) 
         enddo 
        enddo 
       enddo 
      else 
       do ik=1,NTK 
        do ib=1,Nb(ik) 
         SUM_CMPX=SUM_CMPX+CONJG(UNT(ib,iw,ik))*SX(ib,ib,ik)*UNT(ib,jw,ik)*pf(ia1,ia2,ia3,ik) 
        enddo 
       enddo 
      endif 
      MAT_SX_R(iw,jw,ia1,ia2,ia3)=SUM_CMPX/DBLE(NTK)
     enddo 
    enddo 
   enddo 
  enddo 
 enddo 
 deallocate(SX,pf) 
 !
 !OPEN(155,FILE='./dat.sx_mat_r') 
 !
 OPEN(155,FILE='./dat.sx_mat_r') 
 rewind(155)
 write(155,'(a)')'#sx_mat_r'
 write(155,'(a)')'#1:R1, 2:R2, 3:R3 (lattice vector)'
 write(155,'(a)')'#1:i, 2:j, 3:Re(Sx_ij) [eV], 4:Im(Sx_ij) [eV]' 
 do ia1=-Na1,Na1 
  do ia2=-Na2,Na2 
   do ia3=-Na3,Na3
    write(155,*) ia1,ia2,ia3           
    do iw=1,NWF
     do jw=1,NWF 
      write(155,'(i5,i5,2f20.10)') iw,jw,MAT_SX_R(iw,jw,ia1,ia2,ia3)*au 
     enddo!jw       
    enddo!iw       
    write(155,*) 
   enddo!ia3
  enddo!ia2
 enddo!ia1
 close(155) 
 !
 !WRITE MAT_SX_R
 !
 write(6,*) 
 write(6,*)'================'
 write(6,*)' MAT_SX_R in eV '
 write(6,*)'================'
 write(6,*) 
 do ia1=-Na1,Na1 
  do ia2=-Na2,Na2 
   do ia3=-Na3,Na3 
    write(6,*) ia1,ia2,ia3           
    do iw=1,NWF
     write(6,'(300F15.8)')(dble(MAT_SX_R(iw,jw,ia1,ia2,ia3)*au),jw=1,NWF)
    enddo!iw  
    write(6,*)  
   enddo!ia3       
  enddo!ia2 
 enddo!ia1 
 !
 !NSK_BAND_DISP=Ndiv*(N_sym_points-1)+1
 !allocate(SK_BAND_DISP(3,NSK_BAND_DISP)); SK_BAND_DISP(:,:)=0.0d0 
 !call makekpts(Ndiv,N_sym_points,NSK_BAND_DISP,SK_sym_pts(1,1),SK_BAND_DISP(1,1))
 !
 !H_MAT_R=H_MAT_R-MAT_SX_R 
 !call calc_band_disp(Ndiv,N_sym_points,NTK,NSK_BAND_DISP,Na1,Na2,Na3,NWF,SK_BAND_DISP(1,1),&
 !H_MAT_R(1,1,-Na1,-Na2,-Na3),nkb1,nkb2,nkb3,a1(1),a2(1),a3(1),b1(1),b2(1),b3(1))  
 !
endif!myrank.eq.0
!---
!end_time=MPI_Wtime()
!diff_time=end_time-start_time 
!write(6,*) '# TOTAL TIME=',diff_time 
!call MPI_FINALIZE(ierr)
!!---
!STOP
!END           
!
!--
!VXC 
!--
if(myrank.eq.0)then
 !
 !OPEN(151,FILE='./dat.vxc',FORM='unformatted') 
 !
 OPEN(151,FILE='./dat.vxc',FORM='unformatted') 
 REWIND(151)       
 read(151) header 
 read(151) aa,nsnd 
 if(nsnd>1)then
  write(6,*)'nspin*ndens>1; the program not suported; then stop'
  stop
 endif 
 read(151) nrx2,nry2,nrz2; write(6,*)'nrx2,nry2,nrz2=',nrx2,nry2,nrz2 
 allocate(vxcr(nrx2,nry2,nrz2));vxcr(:,:,:)=0.0d0 
 read(151) vxcr(:,:,:) 
 close(151)  
 write(6,*)'FINISH REDING VXCR'
 !fft
 nfft1=nrx2+1; nfft2=nry2+1; nfft3=nrz2+1; Nl123=nfft1*nfft2*nfft3 
 call fft3_init(nrx2,nry2,nrz2,nfft1,nfft2,nfft3,fs) 
 allocate(fftwk(Nl123*2),stat=err) 
 allocate(wfunc(Nl123*2),stat=err) 
 !
 allocate(VXCirr(Mb,Mb,Nk_irr));VXCirr(:,:,:)=0.0d0 
 allocate(MAT_VXC(Mb,Mb));MAT_VXC=0.0d0             
 do ik=1,Nk_irr 
  write(6,*)'ik=',ik 
  MAT_VXC(:,:)=0.0d0 
  call calc_MAT_VXC(Nk_irr,Mb,NTG,NGI(1),KGI(1,1,1),vxcr(1,1,1),CIR(1,1,ik),ik,&
  nrx2,nry2,nrz2,nfft1,nfft2,Nl123,wfunc(1),fftwk(1),fs,MAT_VXC(1,1))
  VXCirr(:,:,ik)=MAT_VXC(:,:)
 enddo!ik
 deallocate(MAT_VXC,vxcr) 
!--
 allocate(VXC(Mb,Mb,NTK));VXC(:,:,:)=0.0d0 
 do ik=1,NTK 
  ikir=numirr(ik) 
  if(trs(ik)==1)then 
   VXC(:,:,ik)=VXCirr(:,:,ikir) 
  elseif(trs(ik)==-1)then 
   do ib=1,Nb(ik)
    do jb=1,Nb(ik) 
     VXC(ib,jb,ik)=VXCirr(jb,ib,ikir) 
    enddo 
   enddo 
  endif 
 enddo 
 write(6,*)'#finish make VXC'
 !
 !OPEN(152,W,FILE='energy_vs_vxc') 
 !
 OPEN(152,FILE='./dat.energy_vs_vxc') 
 rewind(152) 
 do ik=1,Nk_irr 
  do jb=1,Mb 
   en=(E_EIGI(jb+Ns(ik),ik)-FermiEnergy)*au 
   write(152,'(3f20.10)') en,VXCirr(jb,jb,ik)*au 
  enddo 
 enddo 
 deallocate(VXCirr) 
 !
 !--
 allocate(pf(-Na1:Na1,-Na2:Na2,-Na3:Na3,NTK));pf=0.0d0     
 do ik=1,NTK 
  do ia3=-Na3,Na3 
   do ia2=-Na2,Na2 
    do ia1=-Na1,Na1 
     phase=tpi*(SK0(1,ik)*dble(ia1)+SK0(2,ik)*dble(ia2)+SK0(3,ik)*dble(ia3)) 
     pf(ia1,ia2,ia3,ik)=exp(-ci*phase) 
    enddo  
   enddo  
  enddo  
 enddo  
 write(6,*)'#finish make pf'
 !--
 allocate(MAT_VXC_R(NWF,NWF,-Na1:Na1,-Na2:Na2,-Na3:Na3));MAT_VXC_R(:,:,:,:,:)=0.0d0 
 do ia3=-Na3,Na3
  do ia2=-Na2,Na2
   do ia1=-Na1,Na1
    do iw=1,NWF
     do jw=1,NWF 
     SUM_CMPX=0.0D0 
     if(.true.)then!.true.=off-diag, .false.=diag
      do ik=1,NTK 
       do ib=1,Nb(ik) 
        do jb=1,Nb(ik) 
         SUM_CMPX=SUM_CMPX+CONJG(UNT(ib,iw,ik))*VXC(ib,jb,ik)*UNT(jb,jw,ik)*pf(ia1,ia2,ia3,ik) 
        enddo 
       enddo 
      enddo 
     else 
      do ik=1,NTK 
       do ib=1,Nb(ik) 
        SUM_CMPX=SUM_CMPX+CONJG(UNT(ib,iw,ik))*VXC(ib,ib,ik)*UNT(ib,jw,ik)*pf(ia1,ia2,ia3,ik) 
       enddo 
      enddo 
     endif 
     MAT_VXC_R(iw,jw,ia1,ia2,ia3)=SUM_CMPX/DBLE(NTK)
     enddo 
    enddo 
   enddo 
  enddo 
 enddo 
 deallocate(VXC,pf) 
 !
 !OPEN(153,W,FILE='vxc_mat_r') 
 !
 OPEN(153,FILE='./dat.vxc_mat_r') 
 rewind(153)
 write(153,'(a)')'#vxc_mat_r'
 write(153,'(a)')'#1:R1, 2:R2, 3:R3 (lattice vector)'
 write(153,'(a)')'#1:i, 2:j, 3:Re(vxc_ij) [eV], 4:Im(vxc_ij) [eV]' 
 do ia1=-Na1,Na1 
  do ia2=-Na2,Na2 
   do ia3=-Na3,Na3 
    write(153,*) ia1,ia2,ia3           
    do iw=1,NWF
     do jw=1,NWF 
      write(153,'(i5,i5,2f20.10)') iw,jw,MAT_VXC_R(iw,jw,ia1,ia2,ia3)*au 
     enddo!jw  
    enddo!iw  
    write(153,*)  
   enddo!ia3       
  enddo!ia2 
 enddo!ia1  
 close(153) 
!--
 write(6,*) 
 write(6,*)'================='
 write(6,*)' MAT_VXC_R in eV '
 write(6,*)'================='
 write(6,*) 
 do ia1=-Na1,Na1 
  do ia2=-Na2,Na2 
   do ia3=-Na3,Na3 
    write(6,*) ia1,ia2,ia3           
    do iw=1,NWF
     write(6,'(300F15.8)')(dble(MAT_VXC_R(iw,jw,ia1,ia2,ia3)*au),jw=1,NWF)
    enddo!iw  
    write(6,*)  
   enddo!ia3       
  enddo!ia2 
 enddo!ia1 
 !
 NSK_BAND_DISP=Ndiv*(N_sym_points-1)+1
 allocate(SK_BAND_DISP(3,NSK_BAND_DISP)); SK_BAND_DISP(:,:)=0.0d0 
 call makekpts(Ndiv,N_sym_points,NSK_BAND_DISP,SK_sym_pts(1,1),SK_BAND_DISP(1,1))
 !
 write(5001) H_MAT_R 
 write(5002) MAT_VXC_R 
 write(5003) MAT_SX_R 
 ! 
 H_MAT_R=H_MAT_R - MAT_VXC_R - MAT_SX_R  
 !
 call calc_band_disp(Ndiv,N_sym_points,NTK,NSK_BAND_DISP,Na1,Na2,Na3,NWF,SK_BAND_DISP(1,1),&
 H_MAT_R(1,1,-Na1,-Na2,-Na3),nkb1,nkb2,nkb3,a1(1),a2(1),a3(1),b1(1),b2(1),b3(1))  
endif!myrank.eq.0
!
call MPI_BARRIER(comm,ierr)
end_time=MPI_Wtime()
diff_time=end_time-start_time 
if(myrank.eq.0) then 
 call system('pwd') 
 write(6,*)'#TOTAL TIME=',diff_time 
endif 
!
call MPI_FINALIZE(ierr)
STOP 
END              
