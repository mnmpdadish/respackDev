module m_truncation 
  implicit none
contains
  subroutine truncation(NWF,Na1,Na2,Na3,threshold_e,threshold_r,diff_transfers,a1,a2,a3,wcenter_lat,KS_R)
    implicit none 
    integer,intent(in)::NWF,Na1,Na2,Na3 
    real(8),intent(in)::threshold_e,threshold_r,diff_transfers  
    real(8),intent(in)::a1(3),a2(3),a3(3) 
    real(8),intent(in)::wcenter_lat(3,NWF) 
    complex(8),intent(inout)::KS_R(NWF,NWF,-Na1:Na1,-Na2:Na2,-Na3:Na3)   
    complex(8)::tij
    complex(8),allocatable::tr(:)
    real(8),allocatable::dist(:)
    integer,allocatable::ib_map(:)
    integer,allocatable::jb_map(:)
    integer,allocatable::ia1_map(:)
    integer,allocatable::ia2_map(:)
    integer,allocatable::ia3_map(:)
    integer::ia1,ia2,ia3,ib,jb
    integer::Ntr,itr 
    real(8)::vec_rij(3),rij 
    !
    Ntr=NWF*NWF*(2*Na3+1)*(2*Na2+1)*(2*Na1+1) 
    allocate(tr(Ntr)); tr=0.0d0  
    allocate(dist(Ntr)); dist=0.0d0  
    allocate(ib_map(Ntr)); ib_map=0 
    allocate(jb_map(Ntr)); jb_map=0 
    allocate(ia1_map(Ntr)); ia1_map=0
    allocate(ia2_map(Ntr)); ia2_map=0
    allocate(ia3_map(Ntr)); ia3_map=0
    !
    do ia1=-Na1,Na1
     do ia2=-Na2,Na2
      do ia3=-Na3,Na3
       do ib=1,NWF
        do jb=1,NWF 
         tij=KS_R(ib,jb,ia1,ia2,ia3)
         !
         if(threshold_e>abs(tij))then 
          tij=0.0d0 
         endif 
         !
         vec_rij(:)=(dble(ia1)+wcenter_lat(1,jb)-wcenter_lat(1,ib))*a1(:)&
                   +(dble(ia2)+wcenter_lat(2,jb)-wcenter_lat(2,ib))*a2(:)&
                   +(dble(ia3)+wcenter_lat(3,jb)-wcenter_lat(3,ib))*a3(:)
         rij=dsqrt(vec_rij(1)**2+vec_rij(2)**2+vec_rij(3)**2) 
         if(threshold_r<abs(rij))then 
          !write(6,*)'rij=',rij 
          tij=0.0d0 
         endif 
         !
         KS_R(ib,jb,ia1,ia2,ia3)=tij 
         ! 
         itr=jb +(ib-1)*NWF +(ia3+Na3)*NWF*NWF +(ia2+Na2)*NWF*NWF*(2*Na3+1) +(ia1+Na1)*NWF*NWF*(2*Na3+1)*(2*Na2+1)  
         ib_map(itr)=ib
         jb_map(itr)=jb
         ia3_map(itr)=ia3
         ia2_map(itr)=ia2
         ia1_map(itr)=ia1
         tr(itr)=tij
         dist(itr)=rij
         !
        enddo!jb 
       enddo!ib 
      enddo!ia3 
     enddo!ia2 
    enddo!ia1 
    !
    call wrt_transfer_analysis(threshold_e,threshold_r,diff_transfers,Ntr,ib_map(1),jb_map(1),ia1_map(1),ia2_map(1),ia3_map(1),dist(1),tr(1)) 
    !
    deallocate(tr,dist,ib_map,jb_map,ia1_map,ia2_map,ia3_map)
    !
    return
  end subroutine truncation 

  subroutine wrt_transfer_analysis(threshold_e,threshold_r,diff_transfers,Ntr,ib_map,jb_map,ia1_map,ia2_map,ia3_map,dist,tr) 
    implicit none 
    integer,intent(in)::Ntr
    real(8),intent(in)::threshold_e,threshold_r,diff_transfers   
    integer,intent(in)::ib_map(Ntr)
    integer,intent(in)::jb_map(Ntr)
    integer,intent(in)::ia1_map(Ntr)
    integer,intent(in)::ia2_map(Ntr)
    integer,intent(in)::ia3_map(Ntr)
    real(8),intent(in)::dist(Ntr)
    complex(8),intent(in)::tr(Ntr)
    real(8)::tri,trj 
    integer::itr,jtr
    integer,allocatable::skip(:)
    real(8),parameter::au=27.21151d0 
    real(8),parameter::del_zero=(1.d-5)/au  
    !
    !OPEN(122,W,FILE='./dat.tr') 
    !
    OPEN(122,FILE='./dat.tr') 
    rewind(122)
    write(122,'(a30,x,f10.5)')'#Threshold for transfer (eV):',threshold_e*au
    write(122,'(a30,x,f10.5)')'#Threshold for distance (AA):',threshold_r 
    write(122,'(a30,x,f10.5)')'#difference in transfers (eV):',diff_transfers*au  
    write(122,'(a30)')'#i,j,ia1,ia2,ia3,tr,dist:'
    allocate(skip(Ntr)); skip=0
    do itr=1,Ntr 
     tri=abs(dble(tr(itr)))  
     if(tri<del_zero)cycle 
     if(skip(itr)==1)cycle 
     write(122,*) 
     write(122,'(a30)')'irreducible transfer:' 
     !write(122,'(5i5,2f15.7)') ib_map(itr),jb_map(itr),ia1_map(itr),ia2_map(itr),ia3_map(itr),dble(tr(itr))*au,dist(itr)
     write(122,'(5i3,2f10.5)') ib_map(itr),jb_map(itr),ia1_map(itr),ia2_map(itr),ia3_map(itr),dble(tr(itr))*au,dist(itr)
     write(122,*) 
     do jtr=1,Ntr 
      trj=abs(dble(tr(jtr)))  
      if(itr==jtr)cycle 
      if(skip(jtr)==1)cycle 
      if(trj<del_zero)cycle 
      if(abs(tri-trj)<diff_transfers)then 
       !write(122,'(5i5,2f15.7)') ib_map(jtr),jb_map(jtr),ia1_map(jtr),ia2_map(jtr),ia3_map(jtr),dble(tr(jtr))*au,dist(jtr)
       write(122,'(5i3,2f10.5)') ib_map(jtr),jb_map(jtr),ia1_map(jtr),ia2_map(jtr),ia3_map(jtr),dble(tr(jtr))*au,dist(jtr)
       !write(122,'(2i10,2f15.7)') itr,jtr,dble(tr(itr))*au,dble(tr(jtr))*au  
       skip(jtr)=1
      endif 
     enddo 
     skip(itr)=1 
    enddo 
    deallocate(skip)
    !
    return
  end subroutine wrt_transfer_analysis 
  !
  subroutine set_kgrid(dense,kvec)
    IMPLICIT NONE
    integer,intent(in)::dense(3)!Calculated k-grid 
    REAL(8),intent(out)::kvec(3,PRODUCT(dense(1:3))) 
    INTEGER::ik,i1,i2,i3,s1,s2,s3 
    !
    WRITE(6,'(a50,x,3i10)')'Calculated k-grid:',dense(1:3)
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
    WRITE(6,'(a50,x,i10)')'Total number of k-grid:',ik 
    return
  end subroutine set_kgrid 
  !
  subroutine set_FR(NWF,nkb1,nkb2,nkb3,kgd1,kgd2,kgd3,Na1,Na2,Na3,La1,La2,La3,HR,FR) 
    implicit none 
    integer,intent(in)::NWF,nkb1,nkb2,nkb3,kgd1,kgd2,kgd3,Na1,Na2,Na3,La1,La2,La3 
    complex(8),intent(inout)::HR(NWF,NWF,-Na1:Na1,-Na2:Na2,-Na3:Na3) 
    complex(8),intent(out)::FR(NWF,NWF,-La1:La1,-La2:La2,-La3:La3) 
    real(8),allocatable::WR(:,:,:) 
    integer::i,j,k 
    !integer::ib,jb  
    !complex(8)::tij
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
    write(6,'(a50,f15.8,i8)')'SUM_WEIGHT, ncalck:',SUM_REAL,ncalck 
    if(abs(SUM_REAL-dble(ncalck))>1.0d-6)then 
     stop 'SUM_WEIGHT/=ncalck'
    endif 
    return
  end subroutine get_weightR    
  !
end module m_truncation 
