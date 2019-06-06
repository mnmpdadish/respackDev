!config.h 
      implicit none 
      real(8),parameter::au=27.21151d0
      real(8),parameter::bohr=0.529177249d0!20180412
      real(8),parameter::pi=dacos(-1.0d0)
      real(8),parameter::tpi=2.0d0*pi 
      integer,parameter::NBMAX=12
      integer,parameter::maxshell=20 
      integer,parameter::verbosity=0
      complex(8),parameter::ci=(0.0D0,1.0D0) 
      complex(8),parameter::tci=(0.0D0,2.0D0) 
!param-bandcalc 
      integer::nsymq !TOTAL NUMBER OF SYMMETRY OPERATORS  
      integer::NTK!TOTAL NUMBER OF K POINTS IN MK MESHES 
      integer::NTG!TOTAL NUMBER OF G VECTORS
      integer::NTB_glb!TOTAL NUMBER OF CALCULATED BAND IN BAND CALC. 
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
      real(8)::aa1(3),aa2(3),aa3(3)
      real(8)::b1(3),b2(3),b3(3)
      real(8)::a,b,c,alp,bet,gmm  
!sym
      integer,allocatable::rg(:,:,:)!rg(3,3,nsymq)
      integer,allocatable::pg(:,:)!pg(3,nsymq)
      integer,allocatable::numirr(:)!numirr(NTK) 
!      integer,allocatable::reorder(:)!reorder(NTK) 
!      real(8)::reorder_value=0.0d0,dummy1=0.0d0,dummy2=0.0d0
!      real(8)::reorder_value_max=0.0d0
!      real(8)::reorder_value_min=0.0d0
!      integer::ivalue1 = 0
!      integer,allocatable::is_irr(:)!numirr(NTK) 
      integer,allocatable::numrot(:)!numrot(NTK) 
      integer,allocatable::trs(:)!trs(NTK) 
      integer,allocatable::RW(:,:)!RW(3,NTK)
      real(8),allocatable::rginv(:,:,:)!rginv(3,3,nsymq) 
      complex(8),allocatable::rinv_SO(:,:,:)!rinv_SO(2,2,nsymq) 
      
      integer::ik,jk,iik,iop,i,j,k,l,ib,jb,kb,llb 
      integer::ig,jg,kg,lg,j1,k1,iw 
      integer::L1,L2,L3,nnp  
      integer::RWtmp(3) 
!sample-k
      real(8),allocatable::SKI(:,:) 
      real(8),allocatable::SK0(:,:)!SK0(3,NTK) 
      real(8)::ktmp(3) 
      integer::Nk_irr 
!wfn
      integer,allocatable::NGI(:),KGI(:,:,:) 
      real(8),allocatable::LKGI(:,:)           
      real(8),allocatable::E_EIGI(:,:)           
      real(8),allocatable::E_EIG(:,:)!E_EIG(NTB,NTK)           
      complex(8),allocatable::V_EIG(:,:,:)!V_EIG(NTB,NTB,NTK)           
      complex(8),allocatable::CIR(:,:,:,:)!CIR(NTG,ncomp,NTB,Nk_irr)
      complex(8),allocatable::C0(:,:,:,:)!C0(NTG,ncomp,NTB,NTK)
      complex(8),allocatable::C0_BRA(:,:,:)!C0_BRA(0:NTG,ncomp,NTB)    
      complex(8),allocatable::C0_KET(:,:,:)!C0_KET(0:NTG,ncomp,NTB)    
      complex(8),allocatable::C0_TMP_1(:,:,:)!C0_TMP_1(NTG,ncomp,NTB)    
      complex(8),allocatable::C0_TMP_2(:,:,:)!C0_TMP_2(NTG,ncomp,NTB)    
      complex(8),allocatable::C0WN(:,:,:,:)
      integer::ncomp 
      integer::NWF
      integer::NBAND 
      integer::NB_start!Initial band number for wannier calc
      integer::NB_end!The end band number for wannier calc
!kg 
      integer,allocatable::KG0(:,:,:)!KG0(3,NTG,NTK)
      integer,allocatable::KGtmp(:,:)!KGtmp(3,NTG) 
      integer,allocatable::NG0(:)!NG0(NTK)
      integer,allocatable::NB0(:)!NB0(NTK) 
      integer,allocatable::packing(:,:,:,:) 
      integer::NKG_tmp,NG_for_psi  
!interstate matrix 
      integer,allocatable::KPT(:,:)!KPT(NTK,NB)         
      integer,allocatable::N_BAND(:)!N_BAND(NTK)         
      integer,allocatable::N_BAND_BTM(:)!N_BAND_BTM(NTK) 
      real(8),allocatable::b_LAT(:,:)!b_LAT(3,NB)  
      real(8),allocatable::VEC_b(:,:)!VEC_b(3,NB)  
      real(8),allocatable::WEIGHT_b(:)!WEIGHT_b(NB) 
      complex(8),allocatable::OVERLAP(:,:,:,:)
      complex(8),allocatable::OVERLAP_TMP(:,:)!OVERLAP_TMP(NTB,NTB)
      integer::NTB,NGAUSS,NB
      integer,allocatable::SHIFT_b(:)!SHIFT_b(3)
      integer::ik_ib,ibvec,i_band,j_band,k_band,l_band,Mb
      real(8)::WEIGHT,La,Lb,Lc,s,t,u,v 
      real(8)::VEC_b2
      complex(8)::DET
!inner window 
      integer,allocatable::N_BAND_inner(:)!N_BAND_inner(NTK)
      integer,allocatable::N_BAND_BTM_inner(:)!N_BAND_BTM_inner(NTK)
      integer,allocatable::inner(:,:)!inner(NTB,NTK)
      complex(8),allocatable::P_MAT(:,:)!P_MAT(NTB,NTB)
      complex(8),allocatable::Qin(:,:)!Qin(NTB,NTB)
!triclinic Wannier 
      real(8),allocatable::VEC_d(:)!VEC_d(6)
      real(8),allocatable::bb(:,:)!bb(NBh,6)
      real(8),allocatable::aa(:,:)!aa(NBh,6)
      integer::ab,ix,jx,NBh
!initial guess
      character(10),allocatable::orbtype(:)!orbtype(nigs)20170406 
      integer,allocatable::LGAUSS(:)!LGAUSS(nigs)
      integer,allocatable::MGAUSS(:)!MGAUSS(nigs)
      real(8),allocatable::TAU_GAUSS(:,:)!TAU_GAUSS(3,nigs)
      real(8),allocatable::ALPHA_GAUSS(:)!ALPHA_GAUSS(nigs)
      real(8),allocatable::NORM_GAUSS(:)!NORM_GAUSS(nigs)
      real(8),allocatable::loc_x(:,:)!loc_x(3,nigs)20170914 
      real(8),allocatable::loc_y(:,:)!loc_y(3,nigs)20170914 
      real(8),allocatable::loc_z(:,:)!loc_z(3,nigs)20170914 
      complex(8),allocatable::loc_s_up(:,:)!loc_s_dn(2,nigs)20190210 
      complex(8),allocatable::loc_s_dn(:,:)!loc_s_dn(2,nigs)20190210 
      complex(8),allocatable::A_MAT(:,:,:)!A_MAT(NTB,nigs,NTK)
      complex(8),allocatable::A_TMP(:,:)!A_TMP(NTB,nigs)
!spillage 
      real(8),allocatable::S_EIG(:)!S_EIG(NGS) 
      real(8),allocatable::Z_EIG(:)!Z_EIG(NTB)     
      real(8),allocatable::eig_tmp(:)!eig_tmp(nm) 
      complex(8),allocatable::M_MAT(:,:,:,:) 
      complex(8),allocatable::S_MAT(:,:)!S_MAT(NGS,NGS) 
      complex(8),allocatable::S_TMP(:,:)!S_TMP(NGS,NGS) 
      complex(8),allocatable::X_MAT(:,:)!X_MAT(NGS,NGS) 
      complex(8),allocatable::S_HALF_MAT(:,:)!S_HALF_MAT(NGS,NGS)
      complex(8),allocatable::Z_MAT(:,:,:)!Z_MAT(NTB,NTB,NTK)
      complex(8),allocatable::Z_TMP(:,:)!Z_TMP(NTB,NTB) 
      complex(8),allocatable::C_MAT(:,:,:)!C_MAT(NTB,NGS,NTK)
      complex(8),allocatable::C_TMP(:,:)!C_TMP(NTB,NTB) 
      complex(8),allocatable::mat_tmp(:,:)!mat_tmp(nm,nm) 
      integer::I_SCF,nm  
      real(8)::OMEGA_I_NEW,OMEGA_I_OLD,DEL_OMEGA_I 
!pseudo eigenvalue
      real(8),allocatable::P_EIG(:)!P_EIG(NGS)    
      real(8),allocatable::PSEUDO_EIG(:,:)!PSEUDO_EIG(NGS,NTK) 
      complex(8),allocatable::P_MAT_IN(:,:)!P_MAT_IN(NGS,NGS)
      complex(8),allocatable::P_MAT_OUT(:,:)!P_MAT_OUT(NGS,NGS)
      complex(8),allocatable::PSEUDO_MAT(:,:,:)!PSEUDO(NGS,NGS,NTK)
!spread 
      real(8),allocatable::wannier_center(:,:)!wannier_center(3,n_occ)
      real(8),allocatable::WF_CHARGE(:,:,:)!WF_CHARGE(n_occ,NTK,NB)
      real(8),allocatable::TR_GRAD(:)!TR_GRAD(NTK)
      real(8),allocatable::TR_GRAD_OLD(:)!TR_GRAD_OLD(NTK)
      complex(8),allocatable::GRADIENT(:,:,:)!GRADIENT(n_occ,n_occ,NTK)
      complex(8),allocatable::DIRECTION(:,:,:)!DIREC(n_occ,n_occ,NTK)
      complex(8),allocatable::DIRECTION2(:,:,:)!DIREC2(n_occ,n_occ,NTK)
      complex(8),allocatable::UNITARY(:,:,:)!UNITARY(n_occ,n_occ,NTK)
      complex(8),allocatable::U_ORG(:,:,:)!U_ORG(n_occ,n_occ,NTK)
      complex(8),allocatable::U_OLD(:,:,:)!U_OLD(n_occ,n_occ,NTK)
      complex(8),allocatable::UNT(:,:,:) 
      real(8)::OMEGA_I,OMEGA_OD,OMEGA_D,SPREAD,SPREAD_OLD,DEL_SPREAD
      real(8)::STEP_LENGTH,DELTA(4),OMEGA(4) 
      integer::ILS,I_STEP      
!iband
      real(8),allocatable::SK_BAND_DISP(:,:)!SK_BAND(3,NSK_BAND_DISP)
      real(8),allocatable::E_TMP(:)!E_TMP(n_occ)  
      real(8),allocatable::E_BAND_DISP(:,:)!E_BAND_DISP(n_occ,NTK)
      real(8),allocatable::DIST_K(:)!DIST_K(NTK)
      real(8),allocatable::dist(:)!dist(0:Nblk-1)
      complex(8),allocatable::H_MAT_K(:,:,:)!H_MAT_K(n_occ,n_occ,NTK) 
      complex(8),allocatable::H_TMP_IN(:,:)!H_TMP_IN(n_occ,n_occ)
      complex(8),allocatable::H_TMP_OUT(:,:,:)!H_TMP_OUT(n_occ,n_occ,NSK_BAND_DISP)
      integer::NSK_BAND_DISP,ia1,ia2,ia3,ks,ke,ia1min,ia2min,ia3min 
      real(8)::DIST_B(3),DIST_KSPACE
      real(8)::PHASE 
      complex(8)::PHASE_FACTOR     
      LOGICAL::REVERSE 
!
!vis-WAN(yoshimoto-fft)  
!
      type(fft3_struct)::fs 
      character(99)::filename 
      integer::nwx2,nwy2,nwz2,nfft1,nfft2,nfft3,Nl123,err
      real(8)::htmp,d1,d2,d3,qwf,h1(3),h2(3),h3(3)   
      integer::ndx2,ndy2,ndz2 
      integer::algn235
      integer::nvx,nvy,nvz   
      real(8),allocatable::fftwk(:)
      real(8),allocatable::wfunc(:)
      real(8),allocatable::WEIGHT_R(:,:,:) 
      complex(8),allocatable::C0_I(:) 
      complex(8),ALLOCATABLE::WANNIER_REALSPACE(:,:,:) 
      complex(8),ALLOCATABLE::pw(:,:,:) 
      complex(8),allocatable::H_MAT_R(:,:,:,:,:) 
      integer::i1,i2,i3,na_grid,nb_grid,nc_grid
      real(8)::mem_size          
!vis-Bloch 
      complex(8),ALLOCATABLE::BF_REALSPACE(:,:,:) 
      integer::amin,amax,bmin,bmax,cmin,cmax 
!basic variables 
      integer::iL,IC,JC,SUM_INT,ierr,chdir            
      real(8)::SUM_REAL,tmp(3)
      complex(8)::SUM_CMPX        
!
!20180921 
!
!atom_position
      integer::nsi 
      character(len=2),allocatable::chemical_species(:)!chemical_species(nsi) 
      real(8),allocatable::asi(:,:)!asi(3,nsi) 
!
!end config.h 
