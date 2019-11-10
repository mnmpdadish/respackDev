module m_FR   
  implicit none
contains
  
  subroutine calculate_FR_from_FKdiag(Mb,NWF,NTK,Na1,Na2,Na3,Nb,SK0,UNT,FK,FR) 
    implicit none
    !
    !in
    !
    integer,intent(in)::Mb,NWF,NTK,Na1,Na2,Na3
    integer,intent(in)::Nb(NTK)
    real(8),intent(in)::SK0(3,NTK)
    complex(8),intent(in)::FK(Mb,NTK)  
    complex(8),intent(in)::UNT(Mb,NWF,NTK)
    !
    !out  
    !
    complex(8),intent(out)::FR(NWF,NWF,-Na1:Na1,-Na2:Na2,-Na3:Na3) 
    !
    !local 
    !
    integer::ik,ia1,ia2,ia3,iw,jw,ib 
    real(8)::phase 
    complex(8),allocatable::pf(:,:,:,:)!pf(-Na1:Na1,-Na2:Na2,-Na3:Na3,NTK) 
    complex(8)::SUM_CMPX 
    !
    !parameter 
    !
    real(8),parameter::au=27.21151d0
    real(8),parameter::pi=dacos(-1.0d0)
    real(8),parameter::tpi=2.0d0*pi 
    complex(8),parameter::ci=(0.0D0,1.0D0) 
    !
    allocate(pf(-Na1:Na1,-Na2:Na2,-Na3:Na3,NTK));pf=0.0d0     
    !
    do ik=1,NTK 
     do ia3=-Na3,Na3 
      do ia2=-Na2,Na2 
       do ia1=-Na1,Na1 
        phase=tpi*(SK0(1,ik)*dble(ia1)+SK0(2,ik)*dble(ia2)+SK0(3,ik)*dble(ia3)) 
        pf(ia1,ia2,ia3,ik)=exp(-ci*phase) 
       enddo  
      enddo  
     enddo  
    enddo  
    write(6,'(a30)')'#finish make pf'
    !
    FR=0.0d0 
    do ia3=-Na3,Na3
     do ia2=-Na2,Na2
      do ia1=-Na1,Na1
       do iw=1,NWF
        do jw=1,NWF 
         SUM_CMPX=0.0D0 
         do ik=1,NTK 
          do ib=1,Nb(ik) 
           SUM_CMPX=SUM_CMPX+CONJG(UNT(ib,iw,ik))*FK(ib,ik)*UNT(ib,jw,ik)*pf(ia1,ia2,ia3,ik) 
          enddo 
         enddo 
         FR(iw,jw,ia1,ia2,ia3)=SUM_CMPX/DBLE(NTK)
        enddo 
       enddo 
      enddo 
     enddo 
    enddo 
    deallocate(pf) 
    write(6,'(a30)')'#finish make FR' 
    return
  end subroutine calculate_FR_from_FKdiag 
  
  subroutine calculate_FR_from_FKoffdiag(Mb,NWF,NTK,Na1,Na2,Na3,Nb,SK0,UNT,FK,FR) 
    implicit none
    !
    !in
    !
    integer,intent(in)::Mb,NWF,NTK,Na1,Na2,Na3
    integer,intent(in)::Nb(NTK)
    real(8),intent(in)::SK0(3,NTK)
    complex(8),intent(in)::FK(Mb,Mb,NTK)  
    complex(8),intent(in)::UNT(Mb,NWF,NTK)
    !
    !out  
    !
    complex(8),intent(out)::FR(NWF,NWF,-Na1:Na1,-Na2:Na2,-Na3:Na3) 
    !
    !local 
    !
    integer::ik,ia1,ia2,ia3,iw,jw,ib,jb  
    real(8)::phase 
    complex(8),allocatable::pf(:,:,:,:)!pf(-Na1:Na1,-Na2:Na2,-Na3:Na3,NTK) 
    complex(8)::SUM_CMPX 
    !
    !parameter 
    !
    real(8),parameter::au=27.21151d0
    real(8),parameter::pi=dacos(-1.0d0)
    real(8),parameter::tpi=2.0d0*pi 
    complex(8),parameter::ci=(0.0D0,1.0D0) 
    !
    allocate(pf(-Na1:Na1,-Na2:Na2,-Na3:Na3,NTK));pf=0.0d0     
    !
    do ik=1,NTK 
     do ia3=-Na3,Na3 
      do ia2=-Na2,Na2 
       do ia1=-Na1,Na1 
        phase=tpi*(SK0(1,ik)*dble(ia1)+SK0(2,ik)*dble(ia2)+SK0(3,ik)*dble(ia3)) 
        pf(ia1,ia2,ia3,ik)=exp(-ci*phase) 
       enddo  
      enddo  
     enddo  
    enddo  
    write(6,'(a30)')'#finish make pf'
    !
    do ia3=-Na3,Na3
     do ia2=-Na2,Na2
      do ia1=-Na1,Na1
       do iw=1,NWF
        do jw=1,NWF 
         SUM_CMPX=0.0D0 
         do ik=1,NTK 
          do ib=1,Nb(ik) 
           do jb=1,Nb(ik) 
            SUM_CMPX=SUM_CMPX+CONJG(UNT(ib,iw,ik))*FK(ib,jb,ik)*UNT(jb,jw,ik)*pf(ia1,ia2,ia3,ik) 
           enddo 
          enddo 
         enddo 
         FR(iw,jw,ia1,ia2,ia3)=SUM_CMPX/DBLE(NTK)
        enddo 
       enddo 
      enddo 
     enddo 
    enddo 
    deallocate(pf) 
    write(6,'(a30)')'#finish make FR' 
    return
  end subroutine calculate_FR_from_FKoffdiag 
  
end module m_FR  
