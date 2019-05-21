module m_eigenstate  
  implicit none
contains
  subroutine calculate_eigenstate(NWF,ncalck,Na1,Na2,Na3,nkb1,nkb2,nkb3,flg_weight,a1,a2,a3,SK0,KS_R,EKS,VKS) 
    implicit none 
    integer,intent(in)::NWF,ncalck,Na1,Na2,Na3,nkb1,nkb2,nkb3 
    integer,intent(in)::flg_weight 
    real(8),intent(in)::a1(3),a2(3),a3(3)
    real(8),intent(in)::SK0(3,ncalck) 
    complex(8),intent(in)::KS_R(NWF,NWF,-Na1:Na1,-Na2:Na2,-Na3:Na3)   
    real(8),intent(out)::EKS(NWF,ncalck)           
    complex(8),intent(out)::VKS(NWF,NWF,ncalck)   
    real(8),allocatable::WEIGHT_R(:,:,:)!WEIGHT_R(-Na1:Na1,-Na2:Na2,-Na3:Na3)
    complex(8),allocatable::pf(:,:,:,:)!pf(-Na1:Na1,-Na2:Na2,-Na3:Na3,ncalck)   
    integer::ia1min,ia2min,ia3min 
    integer::ia1,ia2,ia3,ik 
    real(8)::PHASE 
    real(8),parameter::au=27.21151d0
    real(8),parameter::tpi=2.0d0*dacos(-1.0d0)
    complex(8),parameter::ci=(0.0D0,1.0D0) 
    !
    !WEIGHT_R BY Y.Nomura NOMURA 
    !
    allocate(WEIGHT_R(-Na1:Na1,-Na2:Na2,-Na3:Na3)); WEIGHT_R=1.0d0
    if(flg_weight.eq.1)then
     write(6,*) 
     write(6,'(a50)')'WEIGHT CALCULATED:' 
     write(6,*) 
     call get_weightR(nkb1,nkb2,nkb3,Na1,Na2,Na3,WEIGHT_R(-Na1,-Na2,-Na3))      
    else
     write(6,*) 
     write(6,'(a50)')'WEIGHT NOT CALCULATED:' 
     write(6,*) 
    endif 
    !
    !phase factor 
    !
    allocate(pf(-Na1:Na1,-Na2:Na2,-Na3:Na3,ncalck)); pf=0.0d0 
    do ik=1,ncalck 
     do ia3=-Na3,Na3 
      do ia2=-Na2,Na2 
       do ia1=-Na1,Na1 
        !
        !NEAREST R SEARCH BY Y.Yoshimoto 
        !
        call search_Rmin(ia1,ia2,ia3,nkb1,nkb2,nkb3,a1(1),a2(1),a3(1),ia1min,ia2min,ia3min)
        PHASE=tpi*(SK0(1,ik)*DBLE(ia1min)+SK0(2,ik)*DBLE(ia2min)+SK0(3,ik)*DBLE(ia3min))           
        pf(ia1,ia2,ia3,ik)=EXP(ci*PHASE)*WEIGHT_R(ia1,ia2,ia3)          
       enddo!ia1 
      enddo!ia2 
     enddo!ia3 
    enddo!ik  
    !
    !HKS IN WANNIER BASIS. AND DIAGONALIZE
    !
    EKS=0.0d0 
    VKS=0.0d0 
    call eigenvalue(ncalck,NWF,Na1,Na2,Na3,KS_R(1,1,-Na1,-Na2,-Na3),pf(-Na1,-Na2,-Na3,1),EKS(1,1),VKS(1,1,1))  
    !
    deallocate(WEIGHT_R,pf) 
    !
    return 
  end subroutine 

  subroutine search_Rmin(i,j,k,nkb1,nkb2,nkb3,a1,a2,a3,imin,jmin,kmin)
    implicit none
    integer::i,j,k,nkb1,nkb2,nkb3
    integer::imin,jmin,kmin
    real(8)::a1(3),a2(3),a3(3) 
    integer::nmin,mmin,lmin
    integer::n,m,l 
    real(8)::R_pos(3),R_abs,R_min,R_bfr
    R_pos(:)=dble(i)*a1(:)+dble(j)*a2(:)+dble(k)*a3(:)
    nmin=0;mmin=0;lmin=0
    R_bfr=dsqrt(R_pos(1)**2+R_pos(2)**2+R_pos(3)**2) 
    R_min=R_bfr 
    do n=-3,3
     do m=-3,3
      do l=-3,3
       R_pos(:)=dble(i+n*nkb1)*a1(:)+dble(j+m*nkb2)*a2(:)+dble(k+l*nkb3)*a3(:)
       R_abs=dsqrt(R_pos(1)**2+R_pos(2)**2+R_pos(3)**2) 
       if(R_min>R_abs)then 
        R_min=R_abs
        nmin=n
        mmin=m
        lmin=l
       endif 
      enddo
     enddo
    enddo
    imin=i+nmin*nkb1
    jmin=j+mmin*nkb2
    kmin=k+lmin*nkb3
    return
  end subroutine 

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
     !
     !H(ib,jb):: ib: basis, jb: eigenvector
     !
     ![NOTE] trnspose 
     !
     !VKS(jb,ib):: jb: eigenvector, ib: basis 
     !
     do ib=1,NWF
      do jb=1,NWF
       VKS(jb,ib,ik)=Hin(ib,jb) 
      enddo
     enddo
     !
    enddo!ik 
    deallocate(Hin,E_TMP_R) 
    return
  end subroutine 

  subroutine diagV(nm,mat,eig)
    implicit none 
    integer,intent(in)::nm
    complex(8),intent(inout)::mat(nm,nm)
    real(8),intent(out)::eig(nm)
    integer::LWORK,LRWORK,LIWORK  
    integer,allocatable::iwork_zheevd(:)
    real(8),allocatable::rwork_zheevd(:)
    complex(8),allocatable::work_zheevd(:)
    integer::ind
    real(8)::eps 
    !
    LWORK= 2*nm+nm**2
    LRWORK=1+12*nm+3*nm**2
    LIWORK=3+10*nm 
    allocate(work_zheevd(LWORK));work_zheevd(:)=0.0d0
    allocate(rwork_zheevd(LRWORK));rwork_zheevd(:)=0.0d0
    allocate(iwork_zheevd(LIWORK));iwork_zheevd(:)=0
    eps=1.0d-18
    ind=0                 
    !
    call zheevd("V","U",nm,mat,nm,eig,work_zheevd,LWORK,rwork_zheevd,LRWORK,iwork_zheevd,LIWORK,ind)
    !
    if(ind/=0)then 
     write(6,'(a50,i10)')'ind:',ind 
     stop
    endif 
    !
    deallocate(work_zheevd,rwork_zheevd,iwork_zheevd) 
    return 
  end subroutine

  subroutine diagN(nm,mat,eig)
    implicit none 
    integer,intent(in)::nm
    complex(8),intent(inout)::mat(nm,nm)
    real(8),intent(out)::eig(nm)
    integer::LWORK,LRWORK,LIWORK  
    integer,allocatable::iwork_zheevd(:)
    real(8),allocatable::rwork_zheevd(:)
    complex(8),allocatable::work_zheevd(:)
    integer::ind
    real(8)::eps 
    !
    LWORK= 2*nm+nm**2
    LRWORK=1+12*nm+3*nm**2
    LIWORK=3+10*nm 
    allocate(work_zheevd(LWORK));work_zheevd(:)=0.0d0
    allocate(rwork_zheevd(LRWORK));rwork_zheevd(:)=0.0d0
    allocate(iwork_zheevd(LIWORK));iwork_zheevd(:)=0
    eps=1.0d-18
    ind=0                 
    !
    call zheevd("N","U",nm,mat,nm,eig,work_zheevd,LWORK,rwork_zheevd,LRWORK,iwork_zheevd,LIWORK,ind)
    !
    if(ind/=0)then 
     write(6,'(a50,i10)')'ind:',ind 
     stop
    endif 
    !
    deallocate(work_zheevd,rwork_zheevd,iwork_zheevd) 
    return 
  end subroutine

  subroutine get_weightR(kgd1,kgd2,kgd3,Na1,Na2,Na3,WEIGHT_R) 
    implicit none 
    integer,intent(in)::kgd1,kgd2,kgd3,Na1,Na2,Na3 
    real(8),intent(out)::WEIGHT_R(-Na1:Na1,-Na2:Na2,-Na3:Na3) 
    integer::ia1,ia2,ia3 
    integer::ncalck 
    real(8)::SUM_REAL
    !
    SUM_REAL=0.0d0 
    WEIGHT_R=1.0d0 
    ncalck=kgd1*kgd2*kgd3 
    do ia1=-Na1,Na1
     do ia2=-Na2,Na2
      do ia3=-Na3,Na3
       if(abs(ia1)==Na1.and.mod(kgd1,2)==0.and.Na1/=0) WEIGHT_R(ia1,ia2,ia3)=WEIGHT_R(ia1,ia2,ia3)*0.5d0 
       if(abs(ia2)==Na2.and.mod(kgd2,2)==0.and.Na2/=0) WEIGHT_R(ia1,ia2,ia3)=WEIGHT_R(ia1,ia2,ia3)*0.5d0 
       if(abs(ia3)==Na3.and.mod(kgd3,2)==0.and.Na3/=0) WEIGHT_R(ia1,ia2,ia3)=WEIGHT_R(ia1,ia2,ia3)*0.5d0 
       SUM_REAL=SUM_REAL+WEIGHT_R(ia1,ia2,ia3)
      enddo!ia3
     enddo!ia2
    enddo!ia1
    write(6,'(a50,f15.8,i8)')'SUM_WEIGHT, ncalck:',SUM_REAL, ncalck 
    if(abs(SUM_REAL-dble(ncalck))>1.0d-6)then 
     stop 'SUM_WEIGHT/=ncalck'
    endif
    return
  end subroutine get_weightR    

end module m_eigenstate 
