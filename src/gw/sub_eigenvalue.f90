subroutine make_eks(NTK,NWF,Na1,Na2,Na3,HmatR,pf,EKS,VKS)  
  implicit none 
  integer::NTK,NWF,Na1,Na2,Na3 
  complex(8)::pf(-Na1:Na1,-Na2:Na2,-Na3:Na3,NTK)   
  complex(8)::HmatR(NWF,NWF,-Na1:Na1,-Na2:Na2,-Na3:Na3)   
  complex(8),allocatable::Hin(:,:)!Hin(NWF,NWF)          
  real(8),allocatable::E_TMP_R(:)!E_TMP_R(NWF)           
  integer::ik,ib,jb,ia1,ia2,ia3
  complex(8)::SUM_CMPX 
  real(8)::EKS(NWF,NTK)           
  complex(8)::VKS(NWF,NWF,NTK)  
  !
  !HKS IN WANNIER BASIS. AND DIAGONALIZE
  !
  allocate(Hin(NWF,NWF)); Hin=0.0d0           
  allocate(E_TMP_R(NWF)); E_TMP_R=0.0d0            
  !
  EKS=0.0d0 
  VKS=0.0d0 
  do ik=1,NTK 
   Hin(:,:)=0.0D0               
   do ib=1,NWF      
    do jb=1,NWF      
     SUM_CMPX=0.0D0                    
     do ia1=-Na1,Na1 
      do ia2=-Na2,Na2 
       do ia3=-Na3,Na3 
        SUM_CMPX=SUM_CMPX+HmatR(ib,jb,ia1,ia2,ia3)*pf(ia1,ia2,ia3,ik)  
       enddo!ia3            
      enddo!ia2                       
     enddo!ia1
     Hin(ib,jb)=SUM_CMPX           
    enddo!jb 
   enddo!ib 
   E_TMP_R(:)=0.0D0 
   call diagV(NWF,Hin(1,1),E_TMP_R(1)) 
   EKS(:,ik)=E_TMP_R(:) 
   VKS(:,:,ik)=Hin(:,:) 
  enddo!ik 
  deallocate(Hin,E_TMP_R) 
  !
return
end

subroutine make_emk(NTK,NWF,nsgm,Na1,Na2,Na3,HmatT,pf,EMK)  
  implicit none 
  integer::NTK,NWF,nsgm,Na1,Na2,Na3 
  complex(8)::pf(-Na1:Na1,-Na2:Na2,-Na3:Na3,NTK)   
  !complex(8)::HmatT(NWF,NWF,-Na1:Na1,-Na2:Na2,-Na3:Na3,nsgm)   
  complex(4)::HmatT(NWF,NWF,-Na1:Na1,-Na2:Na2,-Na3:Na3,nsgm)   
  complex(8)::larger_en
  complex(8),allocatable::Hin(:,:)!Hin(NWF,NWF)          
  complex(8),allocatable::E_TMP_C(:)!E_TMP_C(NWF)           
  integer::LDVL,LDVR,LWORK,ind 
  complex(8),ALLOCATABLE::VL(:,:)
  complex(8),ALLOCATABLE::VR(:,:)
  complex(8),ALLOCATABLE::work_zgeev(:)
  real(8),ALLOCATABLE::rwork_zgeev(:)
  integer::ie,ik,ib,jb,ia1,ia2,ia3,i,j
  complex(8)::SUM_CMPX 
  complex(8)::EMK(NWF,NTK,nsgm)           
  !
  do ie=1,nsgm
   LWORK=5*NWF;LDVL=NWF;LDVR=NWF;ind=0                 
!$OMP PARALLEL PRIVATE(ik,ib,jb,ia1,ia2,ia3,SUM_CMPX,Hin,E_TMP_C,VL,VR,work_zgeev,rwork_zgeev,i,j,larger_en)
   allocate(VL(LDVL,NWF))
   allocate(VR(LDVR,NWF))
   allocate(work_zgeev(LWORK))
   allocate(rwork_zgeev(2*NWF))
   allocate(Hin(NWF,NWF)); Hin=0.0d0           
   allocate(E_TMP_C(NWF)); E_TMP_C=0.0d0            
!$OMP DO
   do ik=1,NTK 
    Hin(:,:)=0.0D0               
    do ib=1,NWF      
     do jb=1,NWF      
      SUM_CMPX=0.0D0                    
      do ia1=-Na1,Na1 
       do ia2=-Na2,Na2 
        do ia3=-Na3,Na3 
         SUM_CMPX=SUM_CMPX+HmatT(ib,jb,ia1,ia2,ia3,ie)*pf(ia1,ia2,ia3,ik)  
        enddo!ia3            
       enddo!ia2                       
      enddo!ia1
      Hin(ib,jb)=SUM_CMPX           
     ENDDO!jb 
    ENDDO!ib 
    !write(6,*)Hin
    E_TMP_C(:)=0.0d0  
    call zgeev("N","N",NWF,Hin,NWF,E_TMP_C,VL,NWF,VR,NWF,work_zgeev,LWORK,rwork_zgeev,ind)
    if(ind/=0) write(6,*) ind,ik,ib 
    !
    do i=1,NWF-1
     do j=i+1,NWF 
      if(dble(E_TMP_C(i))>dble(E_TMP_C(j)))then 
       larger_en=E_TMP_C(i)
       E_TMP_C(i)=E_TMP_C(j) 
       E_TMP_C(j)=larger_en 
      endif 
     enddo!j 
    enddo!i 
    do ib=1,NWF  
     EMK(ib,ik,ie)=E_TMP_C(ib) 
    enddo!ib  
    !
   enddo!ik 
!$OMP END DO
    deallocate(VL,VR,work_zgeev,rwork_zgeev,Hin,E_TMP_C) 
!$OMP END PARALLEL
   write(6,*)'ie=',ie 
  enddo !ie 
  !
  !
return
end
