module m_truncation 
  implicit none
contains
  subroutine truncation(NWF,Na1,Na2,Na3,threshold,KS_R)
    implicit none 
    integer,intent(in)::NWF,Na1,Na2,Na3 
    real(8),intent(in)::threshold 
    complex(8),intent(inout)::KS_R(NWF,NWF,-Na1:Na1,-Na2:Na2,-Na3:Na3)   
    complex(8)::tij
    integer::ia1,ia2,ia3,ib,jb
    real(8),parameter::au=27.21151d0 
    !
    !OPEN(122,W,FILE='./dat.tr') 
    !
    OPEN(122,FILE='./dat.tr') 
    rewind(122)
    write(122,'(a,x,f10.5)')'#Threshold for transfer (eV)=',threshold*au
    do ia1=-Na1,Na1
     do ia2=-Na2,Na2
      do ia3=-Na3,Na3
       do ib=1,NWF
        do jb=1,NWF 
        tij=KS_R(ib,jb,ia1,ia2,ia3)
         if(threshold>abs(tij))then 
          tij=0.0d0 
         else 
          write(122,'(5i5,2f15.7)') ib,jb,ia1,ia2,ia3,dble(tij)*au,abs(tij)*au
         endif 
         KS_R(ib,jb,ia1,ia2,ia3)=tij 
        enddo!jb 
       enddo!ib 
      enddo 
     enddo 
    enddo 
    return
  end subroutine truncation 
  !
  subroutine set_kgrid(dense,kvec)
    IMPLICIT NONE
    integer,intent(in)::dense(3)!Calculated k-grid 
    REAL(8),intent(out)::kvec(3,PRODUCT(dense(1:3))) 
    INTEGER::ik,i1,i2,i3,s1,s2,s3 
    !
    WRITE(*,*) " Calculated k-grid : ", dense(1:3)
    !
    s1=dense(1)/2; s2=dense(2)/2; s3=dense(3)/2 
    ik = 0
    DO i1=-(s1-mod(dense(1)+1,2)),s1  
       DO i2=-(s2-mod(dense(2)+1,2)),s2  
          DO i3=-(s3-mod(dense(3)+1,2)),s3  
             ik=ik + 1
             kvec(1:3,ik)=DBLE((/i1, i2, i3/))/DBLE(dense(1:3))
          END DO
       END DO
    END DO
    !
    WRITE(*,*) " Total number of k-grid : ", ik 
    return
  end subroutine set_kgrid 
  !
  subroutine set_FR(NWF,nkb1,nkb2,nkb3,kgd1,kgd2,kgd3,Na1,Na2,Na3,La1,La2,La3,HR,FR) 
    implicit none 
    integer,intent(in)::NWF,nkb1,nkb2,nkb3,kgd1,kgd2,kgd3,Na1,Na2,Na3,La1,La2,La3 
    complex(8),intent(inout)::HR(NWF,NWF,-Na1:Na1,-Na2:Na2,-Na3:Na3) 
    complex(8),intent(out)::FR(NWF,NWF,-La1:La1,-La2:La2,-La3:La3) 
    real(8),allocatable::WR(:,:,:) 
    integer::i,j,k,ib,jb  
    complex(8)::tij
    real(8),parameter::au=27.21151d0 
    !
    allocate(WR(-Na1:Na1,-Na2:Na2,-Na3:Na3)); WR=1.0 
    call get_weightR(nkb1,nkb2,nkb3,Na1,Na2,Na3,WR(-Na1,-Na2,-Na3))      
    !
    !rm weight factor
    !
    do i=-Na1,Na1
     do j=-Na2,Na2
      do k=-Na3,Na3 
       HR(:,:,i,j,k)=HR(:,:,i,j,k)/WR(i,j,k) 
      enddo
     enddo
    enddo
    deallocate(WR) 
    !
    !F(R)<--H(R) 
    !
    do i=-min(Na1,La1),min(Na1,La1)  
     do j=-min(Na2,La2),min(Na2,La2)  
      do k=-min(Na3,La3),min(Na3,La3)  
       FR(:,:,i,j,k)=HR(:,:,i,j,k) 
      enddo
     enddo
    enddo 
    !
    allocate(WR(-La1:La1,-La2:La2,-La3:La3)); WR=1.0 
    call get_weightR(kgd1,kgd2,kgd3,La1,La2,La3,WR(-La1,-La2,-La3))      
    ! 
    !add weight factor
    !
    do i=-La1,La1
     do j=-La2,La2
      do k=-La3,La3 
       FR(:,:,i,j,k)=FR(:,:,i,j,k)*WR(i,j,k) 
      enddo
     enddo
    enddo
    deallocate(WR) 
    !
    !output check 
    !
    !do i=-La1,La1
    ! do j=-La2,La2
    !  do k=-La3,La3 
    !   do ib=1,NWF
    !    do jb=1,NWF 
    !     tij=FR(ib,jb,i,j,k)
    !     write(6,'(5i5,2f15.7)') ib,jb,i,j,k,dble(tij)*au,abs(tij)*au
    !    enddo
    !   enddo 
    !  enddo
    ! enddo
    !enddo  
    ! 
    return
  end subroutine set_FR 
  !
  subroutine get_weightR(kgd1,kgd2,kgd3,Na1,Na2,Na3,WEIGHT_R) 
    implicit none 
    integer,intent(in)::kgd1,kgd2,kgd3,Na1,Na2,Na3 
    real(8),intent(out)::WEIGHT_R(-Na1:Na1,-Na2:Na2,-Na3:Na3) 
    integer::ia1,ia2,ia3 
    integer::ncalck 
    real(8)::SUM_REAL
    !
    !WEIGHT_R BY Y.Nomura NOMURA 
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
    write(6,'(a20,f15.8,i8)')'SUM_WEIGHT,ncalck',SUM_REAL,ncalck 
    if(abs(SUM_REAL-dble(ncalck))>1.0d-6)then 
     stop 'SUM_WEIGHT/=ncalck'
    endif 
    return
  end subroutine get_weightR    
  !
end module m_truncation 
