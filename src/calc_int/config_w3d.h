!CONFIG.H 
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
!avec 
      real(8)::VOLUME 
      real(8)::a1(3),a2(3),a3(3)
      real(8)::b1(3),b2(3),b3(3)
!fft 
      type(fft3_struct)::fs 
      integer::nwx2,nwy2,nwz2,nfft1,nfft2,nfft3,Nl123!,m1,m2,m3
      real(8)::htmp,d1,d2,d3,qwf,h1(3),h2(3),h3(3)   
      integer::algn235
      real(8),allocatable::fftwk(:)
      real(8),allocatable::wfunc(:)
!sym 
      integer,allocatable::rg(:,:,:)!rg(3,3,nsymq)
      integer,allocatable::pg(:,:)!pg(3,nsymq)
      real(8),allocatable::rginv(:,:,:)!rginv(3,3,nsymq) 
      integer::nnp 
!sk
      integer,allocatable::numirr(:)!numirr(NTK)
      integer,allocatable::numrot(:)!numrot(NTK)
      integer,allocatable::trs(:)!trs(NTK)
      integer,allocatable::RW(:,:)!RW(3,NTK)
      real(8),allocatable::SKI(:,:)
      real(8),allocatable::SK0(:,:)!SK0(3,NTK)
      integer::Nk_irr,RWtmp(3) 
      real(8)::ktmp(3) 
!kg
      integer,ALLOCATABLE::NGI(:)!NGI(Nk_irr) 
      integer,allocatable::NG0(:)!NG0(NTK)
      integer,allocatable::KGI(:,:,:) 
      integer,allocatable::KGtmp(:,:)!KGtmp(3,NTG)
      integer,allocatable::KG0(:,:,:)!KG0(3,NTG,NTK)
      integer,allocatable::packing(:,:,:,:) 
      integer,allocatable::packingq(:,:,:,:) 
      real(8),allocatable::LKGI(:,:)           
!wannier  
      real(8),allocatable::coord(:,:)!coord(3,NWF)  
      complex(8),allocatable::C0(:,:,:)
      integer::NG_for_psi,shift_G(3)  
      real(8)::rij(3),rij2,rij1 
!sq 
      real(8),allocatable::SQI(:,:)  
      real(8),allocatable::SQ(:,:)!SQ(3,NTQ)  
      integer,allocatable::numirrq(:)!numirrq(NTQ)  
      integer,allocatable::numrotq(:)!numrotq(NTQ)  
      integer,allocatable::trsq(:)!trsq(NTQ)  
      integer,allocatable::RWq(:,:)!RWq(3,NTQ)
      integer::NTQ!=NTK!TOTAL NUMBER OF Q POINTS 
      integer::Nq_irr 
!lg 
      integer,allocatable::LG0(:,:,:)!LG0(3,NTG,NTQ)    
      integer,allocatable::NGQ_eps(:)!NGQ_eps(NTQ)
      integer,allocatable::NGQ_psi(:)!NGQ_psi(NTQ)  
      integer::NTGQ !TOTAL NUMBER OF G VECTORS FOR EPSILON
!epsqw 
      complex(8),allocatable::em(:)!em(ne)
      complex(8),allocatable::epstmp(:,:,:) 
      complex(8),allocatable::epstmpgm(:,:,:,:) 
      complex(4),allocatable::epsirr(:,:,:,:) 
      integer::Num_freq_grid!frequency grid 20170402
      integer::ne!frequency grid 20170402
      real(8)::Ecut_for_eps!Cutoff energy (Ry) for dielectric matrix
      integer::ierr,chdir
      character(99)::filename,dirname,dum_char 
      logical::file_e 
!dielmat
      integer,allocatable::packtmp(:,:,:) 
      real(8),allocatable::length_qg(:) 
      complex(8),allocatable::rho(:,:,:)!rho(NTG,NWF,NTQ)            
      complex(8),allocatable::prho(:,:)!prho(NTG,NWF)            
      complex(8),ALLOCATABLE::rho_tmp(:) 
      complex(8),allocatable::C0_K(:)!C0_K(NTG)    
      complex(8),allocatable::C0_KQ(:)!C0_KQ(NTG)    
      complex(4),allocatable::epsmk(:,:,:) 
      complex(8),allocatable::func(:,:,:)!func(NWF,NWF,NTQ)
      complex(8),allocatable::funcw(:,:,:,:)!func(NWF,NWF,ne,NTQ)
!     complex(8),allocatable::cwing(:,:)!cwing(NWF,NWF)  
!     complex(8),allocatable::wL(:)!wL(NWF)
!     complex(8),allocatable::wR(:)!wR(NWF)  
      integer::NWF 
      integer::NG_for_eps,iq,ig,jg,ie 
      integer::file_num,igL,jgL,igL1,igL2,igL3 
      real(8)::q1,q2,q3,qsz,qgL1,qgL2,qgL(3)
      real(8)::chead
      integer::iqgm,NoG0
      real(8),allocatable::cheadw(:)!cheadw(ne)
!matW
      complex(8),allocatable::V_MAT(:,:,:,:,:)
      complex(8),allocatable::W_MAT(:,:,:,:,:,:)
!V_MAT(NWF,NWF,-Na1:Na1,-Na2:Na2,-Na3:Na3) 
!W_MAT(NWF,NWF,-Na1:Na1,-Na2:Na2,-Na3:Na3,ne) 
      real(8)::phase 
      complex(8)::pf,SUM_CMPX
!model 
      real(8)::SUM_REAL 
      real(8),allocatable::WEIGHT_R(:,:,:) 
!local 
      integer::ik,ikq,jk,iik,iqir,iop,i,j,iw,jw,ix  
      integer::i1,j1,k1,ia1,ia2,ia3
      integer::L1,L2,L3,Lq1,Lq2,Lq3  
!END CONFIG.H 
