! $Id: bspline.f90,v 1.17 2010/03/29 05:02:49 yosimoto Exp $
! 
! use LAPACK:dgbsv()
!
module bspline
  !private:: bspl, laeq
contains
  subroutine bspl(x,nknot,neknt,gzai,irnk,igz,rn)
    implicit none
    integer,intent(in):: nknot,neknt,irnk
    integer,intent(out):: igz
    real(kind=8),intent(in):: x,gzai(1-neknt:nknot+neknt)
    real(kind=8),intent(out):: rn(0:irnk)

    integer:: j,jn,r

    do j=1,nknot
       if (x .lt. gzai(j)) exit
    end do
    igz = j

    rn = 0.0d0
    rn(1) = 1.0d0
    do r=2,irnk
       do jn=r,1,-1
          j=igz+jn-1
          rn(jn) = &
               ( ( x - gzai(j-r) )*rn(jn-1)/( gzai(j-1) - gzai(j-r) ) &
               + ( gzai(j) - x )*rn(jn)/( gzai(j) - gzai(j-r+1) ) )
       end do
    end do
    return
  end subroutine bspl

  subroutine bsplv(nn,x,nknot,neknt,gzai,irnk,igz,rn)
    implicit none
    integer,intent(in):: nn,nknot,neknt,irnk
    integer,intent(out):: igz(nn)
    real(kind=8),intent(in):: x(nn),gzai(1-neknt:nknot+neknt)
    real(kind=8),intent(out):: rn(nn,0:irnk)

    integer:: i,j,jn,r

    do j=1,nknot
       if (x(1) .lt. gzai(j)) exit
    end do
    igz(1) = j
    do i=2,nn
       do j=igz(i-1),nknot
          if (x(i) .lt. gzai(j)) exit
       end do
       igz(i) = j
    end do

    rn = 0.0d0
    rn(:,1) = 1.0d0
!$omp parallel do private(r,jn,j)
    do i=1,nn
       do r=2,irnk
          do jn=r,1,-1
             j=igz(i)+jn-1
             rn(i,jn) = ( &
                  ( x(i) - gzai(j-r) )*rn(i,jn-1)/( gzai(j-1) - gzai(j-r) ) &
                  + ( gzai(j) - x(i) )*rn(i,jn)/( gzai(j) - gzai(j-r+1) ) )
          end do
       end do
    end do
    return
  end subroutine bsplv

  subroutine laeq(x,nod,nknot,gzai,irnk,la,a,ierr)
    implicit none
    integer,intent(in):: nod,nknot,irnk,la
    integer,intent(out):: ierr
    real(kind=8),intent(in):: x(nod),gzai(1-irnk:nknot+irnk)
    real(kind=8),intent(out):: a(-la-irnk+1:irnk-1,nod)

    integer:: i,j,k,jn,kb,igz
    real(kind=8):: rn(0:irnk)

    ierr = 0
    if (nod.ne.nknot+irnk) then
       ierr = 1
       return
    end if
    do i=1,nod
       if (.not. (gzai(i-irnk).lt.x(i)).and.(x(i).lt.gzai(i)) ) then
          ierr = 2
          return
       end if
    end do
    if (la .lt. 0) then
       ierr = 3
       return
    end if
    a = 0.0d0
    do k=1,nod
       call bspl(x(k),nknot,irnk,gzai,irnk,igz,rn)
       do jn=1,irnk
          j = igz+jn-1
          kb = k-j
          a(kb,j) = rn(jn)
       end do
    end do
    return
  end subroutine laeq

  subroutine evspline(nknot,gzai,irnk,mdrv,c,xx,nn,s)
    implicit none
    integer:: nknot,irnk,mdrv,nn
    real(kind=8),intent(in):: gzai(1-irnk:nknot+irnk)
    real(kind=8),intent(in):: xx(nn),c(nknot+irnk)
    real(kind=8),intent(out):: s(nn)

    integer,allocatable:: igz(:)
    integer:: j,kk,jn,r,inf,sup,m
    real(kind=8):: rn(0:irnk)

    allocate(igz(nn))

!$omp parallel do private(inf,sup,m,j)
    do kk=1,nn
       inf = 1
       sup = nknot
       do while (inf+5 .lt. sup)
          m = (inf+sup)/2
          if (gzai(m) .le. xx(kk)) then
             inf = m
          else
             sup = m
          end if
       end do
       do j=inf,sup
          if (xx(kk) .lt. gzai(j)) exit
       end do
       igz(kk) = j
    end do

!$omp parallel do private(rn,r,jn,j)
!poption tlocal(rn,j)
    do kk=1,nn
       ! begin inline bspl()
       rn = 0.0d0
       rn(1) = 1.0d0
       do r=2,irnk-mdrv
          do jn=r,1,-1
             j=igz(kk)+jn-1
             rn(jn) = ( &
                  ( xx(kk) - gzai(j-r) )*rn(jn-1)/( gzai(j-1) - gzai(j-r) ) &
                  + ( gzai(j) - xx(kk) )*rn(jn)/( gzai(j) - gzai(j-r+1) ) )
          end do
       end do
       ! end inline bspl()

       s(kk) = 0.0d0
       do jn=1,irnk-mdrv
          j = jn+igz(kk)-1
          s(kk) = s(kk)+c(j)*rn(jn)
       end do
    end do

    deallocate(igz)
    return
  end subroutine evspline

  subroutine evsplinev(nknot,gzai,irnk,mdrv,c,xx,nn,s)
    implicit none
    integer:: nknot,irnk,mdrv,nn
    real(kind=8),intent(in):: gzai(1-irnk:nknot+irnk)
    real(kind=8),intent(in):: xx(nn),c(nknot+irnk)
    real(kind=8),intent(out):: s(nn)

    integer:: j,kk,jn
    real(kind=8),allocatable:: rn(:,:)
    integer,allocatable:: igz(:)

    allocate(rn(nn,0:irnk),igz(nn))

    call bsplv(nn,xx,nknot,irnk,gzai,irnk-mdrv,igz,rn)
!$omp parallel do private(j,jn)
    do kk=1,nn
       s(kk) = 0.0d0
       do jn=1,irnk-mdrv
          j = jn+igz(kk)-1
          s(kk) = s(kk)+c(j)*rn(kk,jn)
       end do
    end do

    deallocate(rn,igz)
    return
  end subroutine evsplinev

  subroutine cderiv(nknot,gzai,irnk,mdrv,c)
    implicit none
    integer:: nknot,irnk,mdrv
    real(kind=8),intent(in):: gzai(1-irnk:nknot+irnk)
    real(kind=8),intent(inout):: c(nknot+irnk)

    integer:: i,j

    do j=1,mdrv
       do i=1,nknot+irnk-j
          c(i) = (irnk-j)*(c(i+1)-c(i))/(gzai(i)-gzai(i-irnk+j))
       end do
    end do
    return
  end subroutine cderiv

  subroutine formgzai(nod,irnk,x,nknot,gzai,ierr)
    implicit none
    integer,intent(in):: nod,irnk,nknot
    integer,intent(out):: ierr
    real(kind=8),intent(in):: x(nod)
    real(kind=8),intent(out):: gzai(1-irnk:nknot+irnk)

    real(kind=8):: maxspc
    integer:: i

    ierr = 0

    if (nod.ne.nknot+irnk) then
       ierr = 1
       return
    end if

    maxspc = 0.0d0
    do i=1,irnk
       maxspc = max(maxspc,spacing(x(i)))
    end do
    do i=-irnk+1,0
       gzai(i) = x(1)+maxspc*i
    end do

    maxspc = 0.0d0
    do i=nod-irnk+1,nod
       maxspc = max(maxspc,spacing(x(i)))
    end do
    do i=nknot+1,nknot+irnk
       gzai(i) = x(nod)+maxspc*(i-nknot-1)
    end do

    if (mod(irnk,2).eq.0) then
       do i=1,nknot
          gzai(i) = x(i+irnk/2)
       end do
    else
       do i=1,nknot
          gzai(i) = ( x(i+irnk/2) + x(i+irnk/2+1) )/2
       end do
    end if
    return
  end subroutine formgzai

  subroutine bintrpl(nod,irnk,nknot,gzai,x,m,c,ierr)
    implicit none
    integer,intent(in):: nod,irnk,nknot,m
    integer,intent(out):: ierr
    real(kind=8),intent(in):: x(nod),gzai(1-irnk:nknot+irnk)
    real(kind=8),intent(inout):: c(nod,m)

    external dgbsv

    integer,allocatable:: ipiv(:)
    real(kind=8),allocatable:: a(:,:)

    integer:: la,info

    la = irnk-1
    allocate(a(-la-irnk+1:irnk-1,nod),ipiv(nod))

    call laeq(x,nod,nknot,gzai,irnk,la,a,ierr)
    if (ierr.ne.0) then
       deallocate(a,ipiv)
       return
    end if

    call dgbsv(nod,irnk-1,irnk-1,m,a,2*irnk-1+la,ipiv,c,nod,info)

    deallocate(a,ipiv)
    return
  end subroutine bintrpl

  subroutine bintrpl2(nod1,nod2,irnk,nknot1,nknot2,gzai1,gzai2,x1,x2,c,ierr)
    implicit none
    integer,intent(in):: nod1,nod2,irnk,nknot1,nknot2
    integer,intent(out):: ierr
    real(kind=8),intent(in):: x1(nod1),x2(nod2)
    real(kind=8),intent(in):: gzai1(1-irnk:nknot1+irnk)
    real(kind=8),intent(in):: gzai2(1-irnk:nknot2+irnk)
    real(kind=8),intent(inout):: c(nod1,nod2)

    external dgbsv

    integer,allocatable:: ipiv(:)
    real(kind=8),allocatable:: a(:,:),cc(:,:)

    integer:: la,info
    integer:: i,j

    la = irnk-1
    allocate(a(-la-irnk+1:irnk-1,nod1),ipiv(nod1))

    call laeq(x1,nod1,nknot1,gzai1,irnk,la,a,ierr)
    if (ierr.ne.0) then
       deallocate(a,ipiv)
       return
    end if

    call dgbsv(nod1,irnk-1,irnk-1,nod2,a,2*irnk-1+la,ipiv,c,nod1,info)

    deallocate(a,ipiv)

    allocate(a(-la-irnk+1:irnk-1,nod2),ipiv(nod2),cc(nod2,nod1))

    do j=1,nod2
       do i=1,nod1
          cc(j,i) = c(i,j)
       end do
    end do

    call laeq(x2,nod2,nknot2,gzai2,irnk,la,a,ierr)
    if (ierr.ne.0) then
       deallocate(a,ipiv,cc)
       return
    end if

    call dgbsv(nod2,irnk-1,irnk-1,nod1,a,2*irnk-1+la,ipiv,cc,nod2,info)

    do i=1,nod1
       do j=1,nod2
          c(i,j) = cc(j,i)
       end do
    end do

    deallocate(a,ipiv,cc)
    return
  end subroutine bintrpl2

  subroutine cderiv2(nknot1,nknot2,gzai1,gzai2,irnk,mdrv1,mdrv2,c)
    implicit none
    integer:: nknot1,nknot2,irnk,mdrv1,mdrv2
    real(kind=8),intent(in):: gzai1(1-irnk:nknot1+irnk)
    real(kind=8),intent(in):: gzai2(1-irnk:nknot2+irnk)
    real(kind=8),intent(inout):: c(nknot1+irnk,nknot2+irnk)

    integer:: i,j,p,q

    do q=1,nknot2+irnk
       do j=1,mdrv1
          do i=1,nknot1+irnk-j
             c(i,q) = (irnk-j)*(c(i+1,q)-c(i,q))/(gzai1(i)-gzai1(i-irnk+j))
          end do
       end do
    end do

    do j=1,mdrv2
       do i=1,nknot2+irnk-j
          do p=1,nknot1+irnk
             c(p,i) = (irnk-j)*(c(p,i+1)-c(p,i))/(gzai2(i)-gzai2(i-irnk+j))
          end do
       end do
    end do

    return
  end subroutine cderiv2

  subroutine evspline2(nknot1,nknot2,gzai1,gzai2,irnk,mdrv1,mdrv2,c,xx1,xx2,s)
    implicit none
    integer:: nknot1,nknot2,irnk,mdrv1,mdrv2
    real(kind=8),intent(in):: gzai1(1-irnk:nknot1+irnk)
    real(kind=8),intent(in):: gzai2(1-irnk:nknot2+irnk)
    real(kind=8),intent(in):: xx1,xx2,c(nknot1+irnk,nknot2+irnk)
    real(kind=8),intent(out):: s

    integer:: j,jn,k,kn
    real(kind=8),allocatable:: rn1(:),rn2(:)
    integer:: igz1,igz2

    allocate(rn1(0:irnk),rn2(0:irnk))

    call bspl(xx1,nknot1,irnk,gzai1,irnk-mdrv1,igz1,rn1)
    call bspl(xx2,nknot2,irnk,gzai2,irnk-mdrv2,igz2,rn2)

    s = 0.0d0
!$omp parallel do private(j,k,jn,kn),reduction(+:s),collapse(2)
    do kn=1,irnk-mdrv2
       do jn=1,irnk-mdrv1
          j = jn+igz1-1
          k = kn+igz2-1
          s = s+c(j,k)*rn1(jn)*rn2(kn)
       end do
    end do

    deallocate(rn1,rn2)
    return
  end subroutine evspline2

  subroutine evspline2v(nknot1,nknot2,gzai1,gzai2,irnk,mdrv1,mdrv2, &
       c,xx1,xx2,nn,s)
    implicit none
    integer:: nknot1,nknot2,irnk,mdrv1,mdrv2,nn
    real(kind=8),intent(in):: gzai1(1-irnk:nknot1+irnk)
    real(kind=8),intent(in):: gzai2(1-irnk:nknot2+irnk)
    real(kind=8),intent(in):: xx1(nn),xx2(nn),c(nknot1+irnk,nknot2+irnk)
    real(kind=8),intent(out):: s(nn)

    integer:: j,jn,k,kn,ii
    real(kind=8),allocatable:: rn1(:,:),rn2(:,:)
    integer,allocatable:: igz1(:),igz2(:)

    allocate(rn1(nn,0:irnk),rn2(nn,0:irnk),igz1(nn),igz2(nn))

    call bsplv(nn,xx1,nknot1,irnk,gzai1,irnk-mdrv1,igz1,rn1)
    call bsplv(nn,xx2,nknot2,irnk,gzai2,irnk-mdrv2,igz2,rn2)

!$omp parallel do private(j,k,jn,kn)
    do ii=1,nn
       s(ii) = 0.0d0
       do kn=1,irnk-mdrv2
          do jn=1,irnk-mdrv1
             j = jn+igz1(ii)-1
             k = kn+igz2(ii)-1
             s(ii) = s(ii)+c(j,k)*rn1(ii,jn)*rn2(ii,kn)
          end do
       end do
    end do

    deallocate(rn1,rn2,igz1,igz2)
    return
  end subroutine evspline2v

  subroutine formdelta(nod,irnk,nknot,gzai,x,c,ierr)
    implicit none
    integer,intent(in):: nod,irnk,nknot
    integer,intent(out):: ierr
    real(kind=8),intent(in):: x(nod),gzai(1-irnk:nknot+irnk)
    real(kind=8),intent(out):: c(nod,nod)

    external dgbsv

    integer,allocatable:: ipiv(:)
    real(kind=8),allocatable:: a(:,:)

    integer:: la,info
    integer:: i

    la = irnk-1
    allocate(a(-la-irnk+1:irnk-1,nod),ipiv(nod))

    call laeq(x,nod,nknot,gzai,irnk,la,a,ierr)
    if (ierr.ne.0) then
       deallocate(a,ipiv)
       return
    end if

    c = 0.0d0
    do i = 1,nod
       c(i,i) = 1.0d0
    end do

    call dgbsv(nod,irnk-1,irnk-1,nod,a,2*irnk-1+la,ipiv,c,nod,info)

    deallocate(a,ipiv)
    return
  end subroutine formdelta

  subroutine formconst(nod,irnk,nknot,gzai,x,c,ierr)
    implicit none
    integer,intent(in):: nod,irnk,nknot
    integer,intent(out):: ierr
    real(kind=8),intent(in):: x(nod),gzai(1-irnk:nknot+irnk)
    real(kind=8),intent(out):: c(nod)

    external dgbsv

    integer,allocatable:: ipiv(:)
    real(kind=8),allocatable:: a(:,:)

    integer:: la,info

    la = irnk-1
    allocate(a(-la-irnk+1:irnk-1,nod),ipiv(nod))

    call laeq(x,nod,nknot,gzai,irnk,la,a,ierr)
    if (ierr.ne.0) then
       deallocate(a,ipiv)
       return
    end if

    c = 1.0d0

    call dgbsv(nod,irnk-1,irnk-1,1,a,2*irnk-1+la,ipiv,c,nod,info)

    deallocate(a,ipiv)
    return
  end subroutine formconst

  subroutine bsmooth(nod,x,y,nknot,gzai,irnk,c,ierr)
    implicit none
    integer,intent(in):: nod,nknot,irnk
    integer,intent(out):: ierr
    real(kind=8),intent(in):: x(nod),y(nod),gzai(1-irnk:nknot+irnk)
    real(kind=8),intent(out):: c(nknot+irnk)

    integer:: i,j,k,in,jn,jb,igz
    real(kind=8):: rn(0:irnk)

    real(kind=8),allocatable:: a(:,:)
    integer,allocatable:: ipiv(:)
    integer:: la

    la = irnk-1
    allocate(a(-la-irnk+1:irnk-1,nknot+irnk),ipiv(nknot+irnk))

    a = 0.0d0
    c = 0.0d0
    do k=1,nod
       call bspl(x(k),nknot,irnk,gzai,irnk,igz,rn)
       do in=1,irnk
          i = igz+in-1
          do jn=1,irnk
             j = igz+jn-1
             jb = j-i
             ! a(jb,i) = a(jb,i) + w(k)*rn(jn)*rn(in)
             a(jb,i) = a(jb,i) + rn(jn)*rn(in)
          end do
          ! c(i) = c(i) + w(k)*y(k)*rn(in)
          c(i) = c(i) + y(k)*rn(in)
       end do
    end do

    call dgbsv(nknot+irnk,irnk-1,irnk-1,1,a,2*irnk-1+la,ipiv,c,nknot+irnk,ierr)

    deallocate(a,ipiv)
    return
  end subroutine bsmooth
end module bspline
