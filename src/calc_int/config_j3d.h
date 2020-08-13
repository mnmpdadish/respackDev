!config.h calc_matJ 
      implicit none        
      real(8),parameter::au=27.21151d0
      real(8),parameter::bohr=0.529177d0 
      real(8),parameter::pi=dacos(-1.0d0)
      real(8),parameter::tpi=2.0d0*pi 
      complex(8),parameter::ci=(0.0D0,1.0D0) 
      complex(8),parameter::tci=(0.0D0,2.0D0) 
!param-bandcalc 
      integer::nsymq!TOTAL NUMBER OF SYMMETRY OPERATORS  
      integer::NTK!TOTAL NUMBER OF K POINTS IN MK MESHES 
      integer::NTG!TOTAL NUMBER OF G VECTORS
      integer::nkb1!SAMPLING K POINTS ALONG b1 VECTOR
      integer::nkb2!SAMPLING K POINTS ALONG b2 VECTOR
      integer::nkb3!SAMPLING K POINTS ALONG b3 VECTOR
      integer::Na1!LATTICE TRANSLATIONL IN a1 20170327
      integer::Na2!LATTICE TRANSLATIONL IN a2 20170327
      integer::Na3!LATTICE TRANSLATIONL IN a3 20170327
      real(8)::Ecut_for_psi!cutoff energy (Ry) for wfn 
      real(8)::Etot!Total Energy  
      real(8)::FermiEnergy  
!fft 
      type(fft3_struct)::fs 
      integer::nwx2,nwy2,nwz2,nfft1,nfft2,nfft3,Nl123!,m1,m2,m3
      real(8)::htmp,d1,d2,d3,qwf,h1(3),h2(3),h3(3)   
      integer::algn235
      real(8),ALLOCATABLE::fftwk(:)
      real(8),ALLOCATABLE::wfunc(:)
!avec
      real(8)::VOLUME 
      real(8)::a1(3),a2(3),a3(3)
      real(8)::b1(3),b2(3),b3(3)
!sym 
      integer,allocatable::rg(:,:,:)!rg(3,3,nsymq)
      integer,allocatable::pg(:,:)!pg(3,nsymq)
      real(8),allocatable::rginv(:,:,:)!rginv(3,3,nsymq) 
      integer::nnp 
!skirr 
      real(8),allocatable::SKI(:,:)
      integer::Nk_irr
      real(8)::ktmp(3) 
!kg 
      integer,ALLOCATABLE::NGI(:)!NGI(Nk_irr) 
      integer,allocatable::KGI(:,:,:) 
      integer,allocatable::KGtmp(:,:)!KGtmp(3,NTG)
      integer,allocatable::packing(:,:,:,:) 
      integer,allocatable::packtmp(:,:,:) 
      real(8),allocatable::LKGI(:,:)           
!sk0
      real(8),allocatable::SK0(:,:)!SK0(3,NTK)
      integer,allocatable::numirr(:)!numirr(NTK)
      integer,allocatable::numrot(:)!numrot(NTK)
      integer,allocatable::trs(:)!trs(NTK)
      integer,allocatable::RW(:,:)!RW(3,NTK)
      integer::RWtmp(3) 
!kg0 
      integer,allocatable::NG0(:)!NG0(NTK) 
      integer,allocatable::KG0(:,:,:)!KG0(3,NTG,NTK)
!wannier 
      real(8),allocatable::coord(:,:)!coord(3,NWF)  
      complex(8),allocatable::C0(:,:,:,:)!C0(NTG,ncomp,NWF,NTK)
      integer::NWF 
      integer::ncomp,ic,jc
!sqirr 
      real(8),allocatable::SQI(:,:)  
      integer::Nq_irr 
!sq 
      real(8),allocatable::SQ(:,:)!SQ(3,NTQ)  
      integer,allocatable::numirrq(:)!numirrq(NTQ)  
      integer,allocatable::numrotq(:)!numrotq(NTQ)  
      integer,allocatable::trsq(:)!trsq(NTQ)  
      integer,allocatable::RWq(:,:)!RWq(3,NTQ)
      integer::NTQ !=NTK!TOTAL NUMBER OF Q POINTS 
!lg 
      integer,allocatable::LG0(:,:,:)!LG0(3,NTG,NTQ)    
      integer,allocatable::NGQ_eps(:)!NGQ_eps(NTQ) 
      integer,allocatable::NGQ_psi(:)!NGQ_psi(NTQ) 
      integer,allocatable::packingq(:,:,:,:) 
      integer::NG_for_eps,NG_for_psi            
      real(8)::q1,q2,q3,qgL1,qgL2,qgL(3)
      integer::NTGQ !TOTAL NUMBER OF G VECTORS FOR EPSILON
!epsqw 
      complex(8),allocatable::em(:)!em(ne)
      complex(8),allocatable::epstmp(:,:,:) 
      complex(4),allocatable::epsirr(:,:,:,:) 
      complex(8),allocatable::epstmpgm(:,:,:,:) 
      integer::Num_freq_grid!frequency grid!20170420 
      integer::ne!frequency grid!20170420 
      real(8)::Ecut_for_eps!Cutoff energy (Ry) for dielectric matrix
      integer::ierr,chdir
      character(99)::filename,dirname,dum_char 
      logical::file_e 
!matJ
      real(8),ALLOCATABLE::length_qg(:) 
      complex(8),allocatable::C0_K(:,:)!C0_K(NTG,ncomp)    
      complex(8),allocatable::C0_KQ(:,:)!C0_KQ(NTG,ncomp)    
      complex(8),ALLOCATABLE::rho_tmp(:) 
      !
      !20200812 Kazuma Nakamura
      !
      !complex(8),ALLOCATABLE::rho(:,:,:)!prho(NG_for_psi,NWF,NWF) 
      !complex(8),allocatable::prho(:,:,:)!prho(NG_for_psi,NWF,NWF) 
      complex(8),ALLOCATABLE::rho(:,:,:,:)!prho(NG_for_psi,NWF/ncomp,NWF/ncomp,ncomp) 
      complex(8),allocatable::prho(:,:,:,:)!prho(NG_for_psi,NWF/ncomp,NWF/ncomp,ncomp) 
      complex(4),allocatable::epsmk(:,:,:) 
      complex(8),allocatable::func(:,:,:)!func(NWF,NWF,NTQ)
      complex(8),allocatable::funcw(:,:,:,:)!funcw(NWF,NWF,NTQ,ne)
      complex(8),allocatable::X_MAT(:,:,:,:,:)
      complex(8),allocatable::J_MAT(:,:,:,:,:,:)
!X_MAT(NWF,NWF,-Na1:Na1,-Na2:Na2,-Na3:Na3) 
!J_MAT(NWF,NWF,-Na1:Na1,-Na2:Na2,-Na3:Na3,ne) 
      !
      !20200314 Kazuma Nakamura
      !
      complex(8),allocatable::gfunc(:,:,:)!gfunc(NWF,NWF,NTQ)
      complex(8),allocatable::gfuncw(:,:,:,:)!gfuncw(NWF,NWF,NTQ,ne)
      complex(8),allocatable::Y_MAT(:,:,:,:,:)
      complex(8),allocatable::K_MAT(:,:,:,:,:,:)
!Y_MAT(NWF,NWF,-Na1:Na1,-Na2:Na2,-Na3:Na3) 
!K_MAT(NWF,NWF,-Na1:Na1,-Na2:Na2,-Na3:Na3,ne) 
      !
      integer::iqgm,NoG0,shift_G(3)
      integer::file_num 
      real(8)::qsz 
      real(8)::chead 
      real(8),allocatable::cheadw(:)!cheadw(ne)  
      complex(8)::pf,SUM_CMPX 
!model 
      real(8)::SUM_REAL 
      real(8),allocatable::WEIGHT_R(:,:,:) 
!local
      integer::ik,ikq,jk,iik,iqir,iop,i,j
      integer::i1,j1,k1,iw,jw,ix 
      integer::L1,L2,L3,Lq1,Lq2,Lq3  
      integer::ia1,ia2,ia3,igL,jgL,igL1,igL2,igL3 
      integer::iq,ig,jg,ie 
!END CONFIG.H 
