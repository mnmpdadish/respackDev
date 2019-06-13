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
    !1. SVD: axis=U*S*VT
    !
    call dgesvd('A', 'A', m, n, axis(1,1), m, S(1), U(1,1), m, VT(1,1), n, work(1), lwork, info)
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
    Pinv=0.0d0
    do i=1,m
     do j=1,m
      Pinv(j,i)=P(i,j)
     enddo
    enddo
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
    if(P(1,3).ne.0.0d0.and.P(3,1).ne.0.0d0)then 
      alp=atan(P(2,3)/P(1,3)) 
      bet=acos(P(3,3)) 
      gmm=atan(-P(3,2)/P(3,1)) 
    endif 
    !
    if(P(1,3).eq.0.0d0.and.P(3,1).eq.0.0d0)then 
      alp=atan(P(2,1)/P(1,1)) 
      bet=0.0d0 
      gmm=0.0d0
    endif 
    !
    write(6,*) 
    write(6,'(a20)')'Euler angle:'
    write(6,*) 
    write(6,'(a20,f15.10)')'alp[deg]:', alp*180.0d0/pi 
    write(6,'(a20,f15.10)')'bet[deg]:', bet*180.0d0/pi 
    write(6,'(a20,f15.10)')'gmm[deg]:', gmm*180.0d0/pi
    write(6,*) 
    !
    !4. Make SU2
    !
    SU2(1,1)= exp(-ci*(alp+gmm)/2.0d0)*cos(bet/2.0d0) 
    SU2(1,2)=-exp(-ci*(alp-gmm)/2.0d0)*sin(bet/2.0d0) 
    SU2(2,1)= exp( ci*(alp-gmm)/2.0d0)*sin(bet/2.0d0) 
    SU2(2,2)= exp( ci*(alp+gmm)/2.0d0)*cos(bet/2.0d0) 
    write(6,*) 
    write(6,'(a20)')'SU(2):'
    write(6,*) 
    do i=1,2
     write(6,'(4f15.10)')(SU2(i,j),j=1,2)
    enddo
    !
    return 
   end subroutine make_SU2_from_local_axis 
  !
end module m_SO3_to_SU2_Local  
