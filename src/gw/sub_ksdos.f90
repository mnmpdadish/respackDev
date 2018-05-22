subroutine calc_dos_KS(NTB,NTK,nkb1,nkb2,nkb3,ndosgrd,dosgrd,E_EIG,SK0,delt,dmnr,dmna,b1,b2,b3,dos) 
  use m_tetrahedron
  implicit none
  integer::NTB,NTK,nkb1,nkb2,nkb3,ndosgrd
  real(8)::dosgrd(ndosgrd) 
  real(8)::E_EIG(NTB,NTK)           
  real(8)::SK0(3,NTK)           
  real(8)::delt,dmnr,dmna 
  real(8)::b1(3),b2(3),b3(3)
  real(8)::dos(ndosgrd) 
  integer::ie,jb,ik,ikb1,ikb2,ikb3
  integer::iomp,omp_get_thread_num  
  real(8)::SUM_REAL 
  !
  !ttrhdrn
  !
  integer,allocatable::imt1(:)!imt1(4*nkb1*nkb2*nkb3*6)  
  integer,allocatable::index_kpt(:,:,:)!index_kpt(nkb1,nkb2,nkb3)    
  complex(8),allocatable::fk_1D(:)!fk_1D(NTK)
  complex(8),allocatable::gk_1D(:)!gk_1D(NTK)
  complex(8),allocatable::fk_3D(:,:,:)!fk_3D(nkb1,nkb2,nkb3) 
  complex(8),allocatable::gk_3D(:,:,:)!gk_3D(nkb1,nkb2,nkb3)
  complex(8),allocatable::xo(:,:,:)!xo(nkb1,nkb2,nkb3) 
  complex(8),allocatable::xow(:,:,:,:)!xow(nsgm,nkb1,nkb2,nkb3)
  !
  real(8),parameter::pi=dacos(-1.0d0)
  real(8),parameter::au=27.21151d0 
  !
  allocate(imt1(4*nkb1*nkb2*nkb3*6)); imt1=0   
  allocate(index_kpt(nkb1,nkb2,nkb3)); index_kpt=0  
  call make_index_kpt(NTK,nkb1,nkb2,nkb3,SK0(1,1),index_kpt(1,1,1)) 
  call ttrhdrn_mkidx(nkb1,nkb2,nkb3,imt1(1),b1(1),b2(1),b3(1))
  !  
  allocate(xow(ndosgrd,nkb1,nkb2,nkb3)); xow=0.0d0 
!$OMP PARALLEL PRIVATE(ie,jb,fk_1D,gk_1D,ik,fk_3D,gk_3D,ikb3,ikb2,ikb1,xo,iomp) 
  allocate(fk_1D(NTK)); fk_1D=0.0d0 
  allocate(gk_1D(NTK)); gk_1D=0.0d0 
  allocate(fk_3D(nkb1,nkb2,nkb3)); fk_3D=0.0d0 
  allocate(gk_3D(nkb1,nkb2,nkb3)); gk_3D=0.0d0 
  allocate(xo(nkb1,nkb2,nkb3)); xo=0.0d0  
!$OMP DO 
  do ie=1,ndosgrd 
   do jb=1,NTB
    fk_1D=0.0d0
    gk_1D=0.0d0
    do ik=1,NTK 
     fk_1D(ik)=1.0d0 
     gk_1D(ik)=cmplx(dosgrd(ie)-E_EIG(jb,ik),-delt) 
    enddo 
    fk_3D=0.0d0 
    gk_3D=0.0d0 
    do ikb3=1,nkb3
     do ikb2=1,nkb2
      do ikb1=1,nkb1 
       ik=index_kpt(ikb1,ikb2,ikb3) 
       fk_3D(ikb1,ikb2,ikb3)=fk_1D(ik)
       gk_3D(ikb1,ikb2,ikb3)=gk_1D(ik)
      enddo 
     enddo 
    enddo 
    xo=0.0d0 
    call ttrhdrn_simple(dmna,dmnr,nkb1,nkb2,nkb3,imt1(1),fk_3D(1,1,1),gk_3D(1,1,1),xo(1,1,1))
    xow(ie,:,:,:)=xow(ie,:,:,:)+xo(:,:,:) 
   enddo!jb 
   iomp=omp_get_thread_num() 
   if(iomp.eq.0) then 
    write(6,*)'#',ie  
   endif 
  enddo!ie
!$OMP END DO 
  deallocate(fk_1D,gk_1D,fk_3D,gk_3D,xo) 
!$OMP END PARALLEL 
  !
  dos=0.0d0 
  do ie=1,ndosgrd 
   SUM_REAL=0.0d0 
   do ikb3=1,nkb3
    do ikb2=1,nkb2
     do ikb1=1,nkb1 
      SUM_REAL=SUM_REAL+dabs(dimag(xow(ie,ikb1,ikb2,ikb3)))/pi  
     enddo 
    enddo 
   enddo 
   dos(ie)=2.0d0*SUM_REAL/dble(NTK)!2 is spin 
  enddo 
  !
  deallocate(imt1,index_kpt,xow) 
return
end
