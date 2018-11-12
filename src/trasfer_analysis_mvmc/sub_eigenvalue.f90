!
!HKS IN WANNIER BASIS. AND DIAGONALIZE
!
subroutine eigenvalue(NTK,NWF,Na1,Na2,Na3,HmatR,pf,EKS,VKS)  
  implicit none 
  integer,intent(in)::NTK,NWF,Na1,Na2,Na3 
  complex(8),intent(in)::pf(-Na1:Na1,-Na2:Na2,-Na3:Na3,NTK)   
  complex(8),intent(in)::HmatR(NWF,NWF,-Na1:Na1,-Na2:Na2,-Na3:Na3)   
  real(8),intent(inout)::EKS(NWF,NTK)           
  complex(8),intent(inout)::VKS(NWF,NWF,NTK)  
  complex(8),allocatable::Hin(:,:)!Hin(NWF,NWF)          
  real(8),allocatable::E_TMP_R(:)!E_TMP_R(NWF)           
  integer::ik,ib,jb,ia1,ia2,ia3
  complex(8)::SUM_CMPX 
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
