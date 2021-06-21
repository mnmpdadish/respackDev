implicit none        
!
!parameters 
!
real(8),parameter::au=27.21151d0
real(8),parameter::bohr=0.529177d0 
real(8),parameter::pi=dacos(-1.0d0)
real(8),parameter::tpi=2.0d0*pi 
complex(8),parameter::ci=(0.0D0,1.0D0) 
complex(8),parameter::tci=(0.0D0,2.0D0) 
!
!MPI 
!
integer::comm,myrank,nproc,ierr 
real(8)::start_time,end_time,diff_time 
!integer::status(MPI_STATUS_SIZE)
integer::nbufq,pnq,bnq,enq
integer::nbufw,pnw,bnw,enw 
!
!integer::file_id 
!
!fft 
!
type(fft3_struct)::fs 
integer::nfft1,nfft2,nfft3,Nl123,err  
real(8),allocatable::fftwk(:)!fftwk(Nl123*2) 
real(8),allocatable::wfunc(:)!wfunc(Nl123*2) 
!
!SC 
!
complex(8),allocatable::vecf(:)!vecf(ne)
complex(8),allocatable::veca(:)!veca(ne) 
complex(8),allocatable::SC(:,:,:,:) 
!
!20180924
!
!complex(8),allocatable::SCirr(:,:,:,:) 
!complex(8),allocatable::pSC(:,:,:,:) 
complex(4),allocatable::SCirr(:,:,:,:) 
complex(4),allocatable::pSC(:,:,:,:) 
!
complex(8),allocatable::pSComp(:,:,:,:) 
complex(4),allocatable::epsmk(:,:,:) 
complex(4),allocatable::pSCR(:,:,:,:,:,:) 
complex(4),allocatable::MAT_SC_R(:,:,:,:,:,:) 
real(8)::sgn,de,delta 
complex(8)::dnm 
complex(8)::psum 
integer::min_ie 
real(8)::min_diff 
!
!SX 
!
!integer::Ncalc 
integer::NG_for_eps,NG_for_psi            
integer::No_G_0 
integer::shift_G(3)  
real(8)::Rc,FACTOR 
real(8)::avec_length(3) 
real(8),allocatable::fbk(:,:)!fbk(NTB,NTK)           
!integer,allocatable::kmq(:)!kmq(NTK) 
real(8),ALLOCATABLE::length_qg(:) 
real(8),ALLOCATABLE::atten_factor(:) 
complex(8),allocatable::SX(:,:,:) 
complex(8),allocatable::SXirr(:,:,:) 
complex(8),allocatable::pSX(:,:,:) 
complex(8),allocatable::SX_DIV(:,:,:) 
complex(8),allocatable::rho(:,:,:) 
complex(8),ALLOCATABLE::rho_tmp(:) 
!complex(8),allocatable::CIRtmp(:,:) 
!complex(8),allocatable::C0(:,:) 
complex(8),allocatable::C0_K(:)!C0_K(NTG)    
complex(8),allocatable::C0_KmQ(:)!C0_KmQ(NTG)    
complex(8),allocatable::MAT_SX_R(:,:,:,:,:)!MAT_SX_R(NWF,NWF,-Na1:Na1,-Na2:Na2,-Na3:Na3) 
real(8)::q1,q2,q3 
real(8)::qgL(3),qgL1,qgL2  
real(8)::chead 
!real(8)::qsz
!complex(8)::psum 
!
!vxcr(151) 
!
character(8)::header 
real(8)::aa(3,3)
integer::nsnd,nrx2,nry2,nrz2 
real(8),allocatable::vxcr(:,:,:)!vxcr(nrx2,nry2,nrz2) 
!
!grid for GW
!
integer::nsgm
integer::nsgmqp
real(8),allocatable::sgmw(:)!sgmw(nsgm)  
real(8),allocatable::sgmwqp(:)!sgmw(nsgmqp)  
real(8)::bandmax,bandmin,diff_band_energy,chiqw_grd_size,emax,emin,ecmax,ecmin 
!
!XC 
!
complex(8),allocatable::MAT_VXC(:,:)
complex(8),allocatable::VXCirr(:,:,:) 
complex(8),allocatable::VXC(:,:,:) 
complex(8),allocatable::MAT_VXC_R(:,:,:,:,:)!MAT_VXC_R(NWF,NWF,-Na1:Na1,-Na2:Na2,-Na3:Na3) 
complex(8),allocatable::pf(:,:,:,:)!pf(-Na1:Na1,-Na2:Na2,-Na3:Na3,NTK)   
real(8)::phase,en 
!
!banddisp
!
real(8),allocatable::SK_BAND_DISP(:,:)!SK_BAND_DISP(3,NSK_BAND_DISP)
integer::NSK_BAND_DISP 
!
!DOS 
!
real(8)::shift_value
real(8),allocatable::ksdos(:)!ksdos(nsgm) 
real(8),allocatable::gwdos(:)!gwdos(nsgm) 
complex(8),allocatable::gw_sigma_dos(:)!gw_sigma_dos(nsgm) 
!
!AKW
!
LOGICAL::REVERSE 
real(8),allocatable::kdata(:)!kdata(NSK_BAND_DISP) 
real(8),allocatable::EKS(:,:)!EKS(NWF,NSK_BAND_DISP) 
real(8),allocatable::gwakw(:,:)!gwakw(NSK_BAND_DISP,nsgm) 
real(8),allocatable::gw_sigma_kw(:,:)!gw_sigma_kw(NSK_BAND_DISP,nsgm) 
complex(8),allocatable::G_MAT_R(:,:,:,:,:)!G_MAT_R(NWF,NWF,-Na1:Na1,-Na2:Na2,-Na3:Na3) 
!
!index 
!
integer::ia1,ia2,ia3,iw,jw,ik,ikir,iop,ib,jb,iq,ikq,kb,ikqir,ikqop,ie,je,iqir   
integer::igL,jgL,igL1,igL2,igL3 
real(8)::SUM_REAL 
real(8)::mem_size          
complex(8)::SUM_CMPX
!
!fileIO
!
integer::chdir
character(99)::filename,command 
complex(4),allocatable::xowtjk(:,:,:)!xowtjk(ne,nsgm,NK_irr) 
complex(4),allocatable::xowtjkq0(:,:,:,:)!xowtjk(ne,nsgm,NK_irr,maxval(Nb)) 
integer::file_num,rec_len,rec_num 
!
!end-of-config.h 
