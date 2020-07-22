module m_SO3_to_SU2_Local  
  implicit none
contains
  subroutine make_SU2_from_local_axis(m,n,axis,SU2) 
    implicit none 
    integer,intent(in)::m 
    integer,intent(in)::n 
    real(8),intent(in)::axis(m,n) 
    complex(8),intent(out)::SU2(2,2) 
    !
    real(8)::S(m),U(m,m),VT(n,n),P(m,m),Pinv(m,m)   
    integer::i,j,k 
    integer::L,lwork,info
    real(8),allocatable::work(:)
    real(8)::SUM_REAL,alp,bet,gmm 
    !
    real(8)::SO3(3,3)
    integer::iloc,jloc,ixyz
    integer::ialp,ibet,igmm
    !
    real(8)::DET
    real(8)::tmp(m,n),tmp_det(m,n) 
    !
    real(8),parameter::pi=dacos(-1.0d0)
    complex(8),parameter::ci=(0.0D0,1.0D0) 
    !
    L=min(m,n)!write(6,*)'L=',L 
    lwork=max(3*L+max(m,n),5*L) 
    lwork=2*lwork
    allocate(work(lwork))
    !
    !do i=1,m
    ! write(6,'(3f15.10)')(axis(i,j),j=1,n)
    !enddo
    !
    !0. Evaluate det(axis) 
    !
    tmp_det=axis
    call calcdet_real(3,tmp_det(1,1),DET) 
    ! 
    !write(6,'(a30,f15.10)')'Determinant of axis=',DET
    !
    tmp=axis
    !
    !write(6,*) 
    !write(6,'(a20)')'axis(before):'
    !write(6,*) 
    !do i=1,m
    ! write(6,'(3f15.10)')(tmp(i,j),j=1,m)
    !enddo 
    !
    if(DET<0.0d0)then
     tmp=-axis
    endif 
    !
    !write(6,*) 
    !write(6,'(a20)')'axis(after):'
    !write(6,*) 
    !do i=1,m
    ! write(6,'(3f15.10)')(tmp(i,j),j=1,m)
    !enddo 
    !
    !1. SVD: axis=U*S*VT
    !
    call dgesvd('A', 'A', m, n, tmp(1,1), m, S(1), U(1,1), m, VT(1,1), n, work(1), lwork, info)
    !
    !write(6,*) 
    !write(6,'(a20)')'axis(after):'
    !write(6,*) 
    !do i=1,m
    ! write(6,'(3f15.10)')(axis(i,j),j=1,m)
    !enddo 
    !write(6,*) 
    !write(6,'(a20)')'U:'
    !write(6,*) 
    !do i=1,m
    ! write(6,'(3f15.10)')(U(i,j),j=1,m)
    !enddo
    !write(6,*) 
    !write(6,'(a20)')'VT:'
    !write(6,*) 
    !do i=1,n
    ! write(6,'(3f15.10)')(VT(i,j),j=1,n)
    !enddo
    !
    !2. P=U*VT
    !
    P=0.0d0 
    do i=1,m
     do j=1,m
      SUM_REAL=0.0d0 
      do k=1,n
       SUM_REAL=SUM_REAL+U(i,k)*VT(k,j)
      enddo 
      P(i,j)=SUM_REAL
     enddo
    enddo
    !
    !Pinv=0.0d0
    !do i=1,m
    ! do j=1,m
    !  Pinv(j,i)=P(i,j)
    ! enddo
    !enddo
    !
    !write(6,*) 
    !write(6,'(a20)')'Orthogonality check:'
    !write(6,*) 
    !do i=1,m
    ! do j=1,m
    !  SUM_REAL=0.0d0
    !  do k=1,m
    !   SUM_REAL=SUM_REAL+P(i,k)*Pinv(k,j)
    !  enddo
    !  write(6,'(a20,2i5,f15.10)')'i,j:',i,j,SUM_REAL 
    ! enddo
    !enddo
    !
    !3. Find Euler angles
    !
    !if(P(1,3).ne.0.0d0.and.P(3,1).ne.0.0d0)then 
    !  alp=atan(P(2,3)/P(1,3)) 
    !  bet=acos(P(3,3)) 
    !  gmm=atan(-P(3,2)/P(3,1)) 
    !endif 
    !
    !if(P(1,3).eq.0.0d0.and.P(3,1).eq.0.0d0)then 
    !  alp=atan(P(2,1)/P(1,1)) 
    !  bet=0.0d0 
    !  gmm=0.0d0
    !endif 
    !--
    !
    !20200719 Kazuma Nakamura; Gimbal Lock
    !
    if(abs(P(1,3))<1.0d-5.and.abs(P(3,1))<1.0d-5)then
     do ialp=1,2
      do ibet=1,2
       do igmm=1,1
        alp=atan(P(2,1)/P(1,1))+dble(ialp-1)*pi 
        bet=0.0d0+dble(ibet-1)*pi  
        gmm=0.0d0 
        !
        SO3(:,:)=0.0d0 
        call make_SO3_matrix(alp,bet,gmm,SO3) 
        !
        sum_real=0.0d0 
        do ixyz=1,3
         do iloc=1,3
          sum_real=sum_real+abs(SO3(ixyz,iloc)-P(ixyz,iloc))
         enddo
        enddo
        !write(6,*)sum_real 
        if(sum_real<1.0d-5)goto 999 
       enddo
      enddo
     enddo
    endif 
    !
    !20200123 Kazuma Nakamura 
    !
    do ialp=1,2
     do ibet=1,2
      do igmm=1,2
       !
       alp=atan(P(2,3)/P(1,3))+dble(ialp-1)*pi 
       bet=acos(P(3,3))+dble(ibet-1)*pi  
       gmm=atan(-P(3,2)/P(3,1))+dble(igmm-1)*pi 
       ! 
       SO3(:,:)=0.0d0 
       call make_SO3_matrix(alp,bet,gmm,SO3) 
       !
       sum_real=0.0d0 
       do ixyz=1,3
        do iloc=1,3
         sum_real=sum_real+abs(SO3(ixyz,iloc)-P(ixyz,iloc))
        enddo
       enddo
       !write(6,*)sum_real 
       if(sum_real<1.0d-5)goto 999 
      enddo
     enddo
    enddo
    !
999 write(6,*) 
    write(6,'(a30,3i5,f10.5)')'ialp, ibet, igmm, sum_real',ialp,ibet,igmm,sum_real 
    write(6,*) 
    write(6,'(a30)')'Euler angle:'
    write(6,*) 
    write(6,'(a30,f15.10)')'alp[deg]:', alp*180.0d0/pi 
    write(6,'(a30,f15.10)')'bet[deg]:', bet*180.0d0/pi 
    write(6,'(a30,f15.10)')'gmm[deg]:', gmm*180.0d0/pi
    write(6,*) 
    !
    !4. Make SU2
    !
    SU2(1,1)= exp(-ci*(alp+gmm)/2.0d0)*cos(bet/2.0d0) 
    SU2(1,2)=-exp(-ci*(alp-gmm)/2.0d0)*sin(bet/2.0d0) 
    SU2(2,1)= exp( ci*(alp-gmm)/2.0d0)*sin(bet/2.0d0) 
    SU2(2,2)= exp( ci*(alp+gmm)/2.0d0)*cos(bet/2.0d0) 
    !
    !(comment out!) SU2=CONJG(SU2)  
    !
    write(6,*) 
    write(6,'(a30)')'SU(2):'
    write(6,*) 
    do i=1,2
     write(6,'(4f15.10)')(SU2(i,j),j=1,2)
    enddo
    !
    return 
  end subroutine make_SU2_from_local_axis 
  !
  !  make SO3 matrix
  !
  subroutine make_SO3_matrix(alp,bet,gmm,SO3) 
    implicit none 
    real(8), intent(in)  :: alp, bet, gmm   
    real(8), intent(out) :: SO3(3,3)
  
    real(8) :: c1, c2, c3, s1, s2, s3
  
    c1 = cos(alp)
    c2 = cos(bet)
    c3 = cos(gmm)
  
    s1 = sin(alp)
    s2 = sin(bet)
    s3 = sin(gmm)
  
    SO3(1,1) = c1*c2*c3 - s1*s3
    SO3(2,1) = s1*c2*c3 + c1*s3
    SO3(3,1) = -s2*c3
  
    SO3(1,2) = -c1*c2*s3 - s1*c3
    SO3(2,2) = -s1*c2*s3 + c1*c3
    SO3(3,2) = s2*s3
  
    SO3(1,3) = c1*s2
    SO3(2,3) = s1*s2
    SO3(3,3) = c2
  
    return
  end subroutine

end module m_SO3_to_SU2_Local  
