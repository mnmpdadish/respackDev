subroutine calc_dos_GW(nkb1,nkb2,nkb3,NTK,NWF,nsgm,idlt,dmna,dmnr,shift_value,b1,b2,b3,sgmw,SK0,EMK,AW) 
 !
 use m_tetrahedron
 !
 implicit none 
 integer::nkb1,nkb2,nkb3,NTK,NWF,nsgm 
 real(8)::idlt,dmna,dmnr,shift_value 
 real(8)::b1(3),b2(3),b3(3),SK0(3,NTK),sgmw(nsgm)  
 complex(8)::EMK(NWF,NTK,nsgm)           
 integer::jb,ie,ikb1,ikb2,ikb3,ik   
 real(8)::delta,ReZ,ImZ,SUM_REAL 
 complex(8)::en  
 real(8)::AW(nsgm) 
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
 allocate(xow(nsgm,nkb1,nkb2,nkb3)); xow=0.0d0 
!$OMP PARALLEL PRIVATE(ie,jb,delta,fk_1D,gk_1D,ik,en,ReZ,ImZ,fk_3D,gk_3D,ikb3,ikb2,ikb1,xo) 
 allocate(fk_1D(NTK)); fk_1D=0.0d0 
 allocate(gk_1D(NTK)); gk_1D=0.0d0 
 allocate(fk_3D(nkb1,nkb2,nkb3)); fk_3D=0.0d0 
 allocate(gk_3D(nkb1,nkb2,nkb3)); gk_3D=0.0d0 
 allocate(xo(nkb1,nkb2,nkb3)); xo=0.0d0  
!$OMP DO 
 do ie=1,nsgm 
  do jb=1,NWF 
   !
   delta=sgmw(ie)+shift_value 
   !
   fk_1D=0.0d0
   gk_1D=0.0d0 
   do ik=1,NTK 
    en=EMK(jb,ik,ie)
    ReZ=delta-dble(en) 
    ImZ=-(idlt+dabs(dimag(en)))!retarded 
    gk_1D(ik)=cmplx(ReZ,ImZ)
    fk_1D(ik)=1.0d0 
   enddo!ik 
   !
   fk_3D=0.0d0
   gk_3D=0.0d0 
   do ikb3=1,nkb3
    do ikb2=1,nkb2
     do ikb1=1,nkb1
      ik=index_kpt(ikb1,ikb2,ikb3)
      fk_3D(ikb1,ikb2,ikb3)=fk_1D(ik)
      gk_3D(ikb1,ikb2,ikb3)=gk_1D(ik)
     enddo!ikb1 
    enddo!ikb2 
   enddo!ikb3 
   !
   xo=0.0d0 
   call ttrhdrn_simple(dmna,dmnr,nkb1,nkb2,nkb3,imt1(1),fk_3D(1,1,1),gk_3D(1,1,1),xo(1,1,1))
   xow(ie,:,:,:)=xow(ie,:,:,:)+xo(:,:,:) 
   !
  enddo!jb 
 enddo!ie
!$OMP END DO
 deallocate(fk_1D,gk_1D,fk_3D,gk_3D,xo) 
!$OMP END PARALLEL
 !
 AW=0.0d0 
 do ie=1,nsgm
  SUM_REAL=0.0d0 
  do ikb3=1,nkb3
   do ikb2=1,nkb2
    do ikb1=1,nkb1
     SUM_REAL=SUM_REAL+dabs(dimag(xow(ie,ikb1,ikb2,ikb3))) 
    enddo!ikb1 
   enddo!ikb2 
  enddo!ikb3 
 AW(ie)=2.0d0*SUM_REAL/dble(NTK)/pi 
 enddo!ie  
 !
 deallocate(imt1,index_kpt,xow) 
 !
RETURN  
END 
