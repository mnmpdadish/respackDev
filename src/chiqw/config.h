!CONFIG.H 
!This program uses the generalized tetrahedron method and implements
!an efficient tool in the GW-LMTO package developed by Fujiwara group
!in University of Tokyo. The algorithm of the generalized tetrahedron
!method is written in the following published papers:
!
![Generalized tetrahedron method] Takeo Fujiwara, Susumu Yamamoto, and
!Yasushi Ishii, J. Phys. Soc. Jpn. 72, No.4, 777-780 (2003).
!(cond.mat/0211159)
!
![GW-LMTO package] Yoshiro Nohara, Susumu Yamamoto, and Takeo Fujiwara,
!Phys. Rev. B 79, 195110 (2009)
!
!One can use this part of the program package in one's developed code
!on the condition of citing these two papers and mentioning the name of
!`the generalized tetrahedron method' in the published paper.
!---
      implicit none 
      real(8),parameter::au=27.21151D0!hartree 
      real(8),parameter::pi=dacos(-1.0d0)
      real(8),parameter::tpi=2.0d0*pi 
      complex(8),parameter::ci=(0.0D0,1.0D0) 
      complex(8),parameter::tci=(0.0D0,2.0D0) 
!param-bandcalc 
      integer::nsymq!TOTAL NUMBER OF SYMMETRY OPERATORS  
      integer::NTK!TOTAL NUMBER OF K POINTS IN MK MESHES 
      integer::NTG!TOTAL NUMBER OF G VECTORS
      integer::NTB!TOTAL NUMBER OF CALCULATED BAND IN BAND CALC. 
      integer::nkb1!SAMPLING K POINTS ALONG b1 VECTOR
      integer::nkb2!SAMPLING K POINTS ALONG b2 VECTOR
      integer::nkb3!SAMPLING K POINTS ALONG b3 VECTOR
      real(8)::Ecut_for_psi!cutoff energy (Ry) for wfn 
      real(8)::Etot!Total Energy  
      real(8)::FermiEnergy  
!sample-k 
      integer::Nk_irr 
      real(8),allocatable::SKI(:,:)!SKI(3,Nk_irr)
      real(8),allocatable::SKI_list(:,:)!SKI_list(3,Nk_irr)
      real(8),allocatable::SK0(:,:)!SK0(3,NTK) 
      real(8)::ktmp(3) 
!wfn 
      integer,allocatable::NGI(:)!NGI(Nk_irr)
      integer,allocatable::NG0(:)!NG0(NTK)
      integer,allocatable::NB0(:)!NB0(NTK)
      real(8),allocatable::E_EIGI(:,:)!E_EIGI(NTB,Nk_irr)
!kg 
      integer,allocatable::KGI(:,:,:)!KGI(3,NTG,Nk_irr)            
      integer,allocatable::KGtmp(:,:)!KGtmp(3,NTG)
      integer,allocatable::KG0(:,:,:)!KG0(3,NTG,NTK)
      integer,allocatable::packing(:,:,:,:)
!packing(-L1:L1,-L2:L2,-L3:L3,Nk_irr)
      real(8),allocatable::LKGI(:,:)           
      integer::NG_for_psi  
!sym
      integer,allocatable::rg(:,:,:)!rg(3,3,nsymq)
      integer,allocatable::pg(:,:)!pg(3,nsymq)
      integer,allocatable::numirr(:)!numirr(NTK) 
      integer,allocatable::numrot(:)!numrot(NTK) 
      integer,allocatable::trs(:)!trs(NTK) 
      integer,allocatable::RW(:,:)!RW(3,NTK)
      real(8),allocatable::rginv(:,:,:)!rginv(3,3,nsymq) 
      integer::RWtmp(3) 
      integer::ik,jk,iik,ikir,iop,i,j,ib,ig,jg 
      integer::i1,j1,k1
      integer::L1,L2,L3,nnp  
!fft 
      type(fft3_struct)::fs 
      integer::nwx2,nwy2,nwz2,nfft1,nfft2,nfft3,Nl123 !,m1,m2,m3  
      real(8)::htmp,d1,d2,d3,qwf,h1(3),h2(3),h3(3)   
      integer::algn235
      real(8),allocatable::fftwk(:)!fftwk(Nl123*2) 
      real(8),allocatable::wfunc(:)!wfunc(Nl123*2) 
!dielmat 
      integer,allocatable::index_kpt(:,:,:)!index_kpt(nkb1,nkb2,nkb3)
      integer,allocatable::WindowInside(:)!WindowInside(NTB)
      integer,allocatable::Tindx(:)!Tindx(NTB)
      real(8),allocatable::E_AVE(:)!E_AVE(NTB)
      real(8),allocatable::band_max(:)!band_max(NTB)
      real(8),allocatable::band_min(:)!band_min(NTB)
      real(8)::bci,bcj 
!polarization
      integer::NG_for_eps,fs1,fs2,No_G_0,N_PAIR,shift_G(3)  
      integer::i_band,j_band,igL,igL1,igL2,igL3,iq,ij,ikq,ikb1,ikb2,ikb3  
      real(8)::q1,q2,q3,del_eps,dsgn 
      integer,allocatable::kpq(:)!kpq(NTK) 
      integer,allocatable::i_for_pair(:)
      integer,allocatable::j_for_pair(:)
!i_for_pair(N_CALC_BAND*N_CALC_BAND)
!j_for_pair(N_CALC_BAND*N_CALC_BAND)
      complex(8)::vm(3) 
      complex(8),allocatable::C0_K(:)!C0_K(NTG)    
      complex(8),allocatable::C0_KQ(:)!C0_KQ(NTG)    
      complex(8),allocatable::e1_1D(:)!e1_1D(NTK)
      complex(8),allocatable::e2_1D(:)!e2_1D(NTK)
      complex(8),allocatable::e1_3D(:,:,:)!e1_3D(nkb1,nkb2,nkb3) 
      complex(8),allocatable::e2_3D(:,:,:)!e2_3D(nkb1,nkb2,nkb3)
      complex(8),allocatable::chiqw(:,:,:)!chiqw(NGeps,NGeps,nen)
      complex(8),allocatable::pchiqw(:,:,:)!pchiqw(NGeps,NGeps,nen)
      complex(8),allocatable::pc0d(:,:,:)!pc0d(NGeps,NGeps,nen) 
      complex(8),allocatable::ism(:,:)!ism(NTK,NGeps)
      complex(8),allocatable::m_tmp(:)!m_tmp(NGeps)
      complex(8),allocatable::epsqw(:,:,:)
      complex(8),allocatable::chi0(:,:)!chi0(NGeps,NGeps)
      complex(8),allocatable::eps_rpa(:,:)!eps_rpa(NGeps,NGeps)
      complex(8),allocatable::chiqwgm(:,:,:,:)
      complex(8),allocatable::pchiqwgm(:,:,:,:)
      complex(8),allocatable::pc0dgm(:,:,:,:)
      complex(8),allocatable::ismgm(:,:,:)
      complex(8),allocatable::epsqwgm(:,:,:,:)
      complex(8),allocatable::eps00w(:,:)
!     complex(8),allocatable::m_tmp_st(:)!m_tmp_st(NGeps) 
!     complex(8),allocatable::ism_st(:,:)!ism_st(NTK,NGeps) 
!epsqw(NGeps,NGeps,nen)
!chiqwgm(NGeps,NGeps,nen,3) 
!pchiqwgm(NGeps,NGeps,nen,3)
!pc0dgm(NGeps,NGeps,nen,3) 
!ismgm(NTK,NGeps,3)
!epsqwgm(NGeps,NGeps,nen,3)
!     integer::jdirec!20170322 
!NonLocalCorrection 
      integer::ix,JB 
      complex(8)::NLvm(3) 
      complex(8),allocatable::MAT_NLr_rNL(:,:,:,:)
!MAT_NLr_rNL(3,NTB,NTB,Nk_irr)
!avec 
      real(8)::VOLUME 
      real(8)::a1(3),a2(3),a3(3)
      real(8)::b1(3),b2(3),b3(3)
!mpi 
      integer::id,ip,iu,N_docc,N_pocc,N_uocc,N_occ,N_vir 
      complex(8),allocatable::OCC(:,:,:)!OCC(NTG,1:nbufo,Nk_irr)
      complex(8),allocatable::VIR(:,:,:)!VIR(NTG,1:nbufv,Nk_irr)
      complex(8),allocatable::VIR_new(:,:,:)
!VIR_new(NTG,1:nbufv,Nk_irr) 
      complex(8),allocatable::Otmp(:,:)!Otmp(NTG,1:nbufo) 
      complex(8),allocatable::Vtmp(:,:)!Vtmp(NTG,1:nbufv) 
      real(8),allocatable::E_OCC(:,:)!E_OCC(1:nbufo,Nk_irr) 
      real(8),allocatable::E_VIR(:,:),E_VIR_new(:,:)
      integer,allocatable::W_OCC(:),W_VIR(:),W_VIR_new(:)
      integer,allocatable::T_OCC(:),T_VIR(:),T_VIR_new(:)
      real(8),allocatable::EO_AVE(:),EV_AVE(:),EV_AVE_new(:)
!E_VIR(1:nbufv,Nk_irr),E_VIR_new(1:nbufv,Nk_irr)
!W_OCC(1:nbufo),W_VIR(1:nbufv),W_VIR_new(1:nbufv)
!T_OCC(1:nbufo),T_VIR(1:nbufv),T_VIR_new(1:nbufv)
!EO_AVE(1:nbufo),EV_AVE(1:nbufv),EV_AVE_new(1:nbufv)
      integer::comm_glb,myrank_glb,nproc_glb
      integer::status(MPI_STATUS_SIZE)
      integer::ie,nbufv,nbufo 
      integer::pnv,bnv,env,newpnv,newbnv,newenv 
      integer::pno,bno,eno!,newpno,newbno,neweno
      integer::pnocld,pnvcld
      integer::tmpnv(3),newtmpnv(3)
      integer::irot!,rcnt,scnt,rotrank
      integer::ncomp 
      real(8)::start_time,end_time,diff_time 
!split comm
      integer::comm,myrank,nproc,dest,source
      integer::id_sub
      integer::comm_perp,myrank_perp,nproc_perp 
      integer::id_sub_perp
      integer::nbufq,pnq,bnq,enq
      real(8),allocatable::SQI(:,:)!SQI(3,nbufq) 
      real(8),allocatable::SQ(:,:)!SQ(3,nbufq) 
      integer,allocatable::LG0(:,:,:)!LG0(3,NTG,nbufq) 
      integer,allocatable::NGQ(:)!NGQ(nbufq) 
      integer::iomp,omp_get_thread_num  
!freq dep  
      integer,allocatable::imt1(:)!imt1(4*nkb1*nkb2*nkb3*6)  
      complex(8),allocatable::em(:)!em(nen)
      complex(8),allocatable::ca1(:)!ca1(4*nen) 
      complex(8),allocatable::xow(:,:,:,:)!xow(nen,nkb1,nkb2,nkb3)
      complex(8),allocatable::xow_1D(:,:)!xow_1D(NTK,nen)
!     complex(8),allocatable::xow_1D_av(:,:)!xow_1D_av(Nk_irr,nen)
      real(8)::emax
!plasma frequency  
      integer::ij_intra 
      integer::N_INTRA
      integer,allocatable::i_for_intra(:)
      integer,allocatable::j_for_intra(:)
      complex(8),allocatable::xo(:,:,:)!xo(nkb1,nkb2,nkb3) 
!     complex(8),allocatable::xo_1D(:)!xo_1D(NTK) 
!     complex(8),allocatable::xo_1D_av(:)!xo_1D_av(Nk_irr) 
      complex(8),allocatable::pk_1D(:,:),pk_3D(:,:,:,:)
      complex(8),allocatable::fk_1D(:),fk_3D(:,:,:)
      complex(8),allocatable::gk_1D(:),gk_3D(:,:,:)
!i_for_intra(N_CALC_BAND*N_CALC_BAND)
!j_for_intra(N_CALC_BAND*N_CALC_BAND)
!pk_1D(NTK,3),pk_3D(nkb1,nkb2,nkb3,3) 
!fk_1D(NTK),fk_3D(nkb1,nkb2,nkb3) 
!gk_1D(NTK),gk_3D(nkb1,nkb2,nkb3)
      real(8)::wd(3),SKT(3) 
      complex(8)::wpl(3),pwpl(3) 
!wannier 
      integer,allocatable::Nb(:)!Nb(NTK)         
      integer,allocatable::Ns(:)!Ns(NTK)         
      complex(8),allocatable::UNT(:,:,:) 
      integer::Mb,NWF  
      integer::iw,jw  
      complex(8),allocatable::prob(:,:)!prob(NTB,NTK)           
      complex(8)::SUM_CMPX        
      complex(8),allocatable::P_OCC(:,:)!P_OCC(1:nbufo,NTK) 
      complex(8),allocatable::P_VIR(:,:)!P_VIR(1:nbufv,NTK) 
      complex(8),allocatable::P_VIR_new(:,:)!P_VIR_new(1:nbufv,NTK) 
!file output 
      integer::ierr,chdir
      integer::file_num_log,file_num_eps!,file_num_chi
      integer::file_num_eps_base!,file_num_chi_base 
      integer::file_num_chi_ie,file_num_eps_ie
      character(99)::filename,command,dirname 
      logical::file_e!,dir_e
      real(8),allocatable::calc_qlist(:,:)!calc_qlist(3,Nk_irr)
      integer::qnum 
      integer::Nq_irr!TOTAL NUMBER OF Q IRR
