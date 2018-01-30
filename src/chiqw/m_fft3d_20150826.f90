! $Id: fft3d.f90,v 1.11 2002/09/18 00:47:12 yosimoto Exp $
module fft_3d
  integer,parameter,private:: nradix=6
  integer,dimension(nradix),parameter,private:: radix=(/6, 8, 3, 4, 2, 5/)
  integer,parameter,private:: nfmax = 30
  logical,parameter,private:: prefer4 = .false.
  private:: fact235,genw,fftv,bftv,fft_tp1,fft_tp,fft_tps1,fft_tps
  private:: fftrd2,bftrd2,fftrd3,bftrd3,fftrd4,bftrd4,fftrd5,bftrd5
  private:: fftrd6,bftrd6,fftrd8,bftrd8
  !
  type fft3_struct
     integer:: nx,ny,nz,nxyz,lx,ly,lz,lxyz
     integer:: nfx,nfy,nfz
     integer,dimension(nfmax):: facx,facy,facz
     real(kind=8),pointer,dimension(:,:):: wx,wy,wz
  end type fft3_struct
contains

  subroutine fact235(n,nf,fac)
    implicit none
    integer,intent(in):: n
    integer,intent(out):: nf,fac(nfmax)
    !
    integer:: p,q,r,j
    integer:: n2,n3,n4,n5,n6,n8
    !
    p = n
    fac(:) = 0
    !
    n2 = 0
    q = p/2
    r = p-2*q
    do while ( (q .gt. 0) .and. (r .eq. 0) )
       n2 = n2 + 1
       p = q
       q = p/2
       r = p-2*q
    end do

    n3 = 0
    q = p/3
    r = p-3*q
    do while ( (q .gt. 0) .and. (r .eq. 0) )
       n3 = n3 + 1
       p = q
       q = p/3
       r = p-3*q
    end do

    n5 = 0
    q = p/5
    r = p-5*q
    do while ( (q .gt. 0) .and. (r .eq. 0) )
       n5 = n5 + 1
       p = q
       q = p/5
       r = p-5*q
    end do

    if (p .ne. 1) then
       write(6,*) n, ' can not be factorized with 2, 3, 5'
       stop
    end if

    n6 = min(n2,n3)
    n2 = n2 - n6
    n3 = n3 - n6
    if (prefer4) then ! prefer radix 4
       n4 = n2/2
       j = mod(n2,2)
       if ( (j .eq. 1) .and. (n4 .ge. 1) ) then
          n8 = 1
          n4 = n4 - 1
       else
          n8 = 0
       end if
    else ! prefer radix 8
       n8 = n2/3
       j = mod(n2,3)
       if (j .eq. 2) then
          n4 = 1
       else if ( (j .eq. 1) .and. (n8 .ge. 1) ) then
          n8 = n8 - 1
          n4 = 2
       else
          n4 = 0
       end if
    end if
       
    n2 = n2 - n4*2 - n8*3
    if (n2 .lt. 0) then
       write(6,*) n2, ' n2 become negative'
       stop
    end if
    if (n2+n3+n4+n5+n6+n8 .gt. nfmax) then
       write(6,*) 'too many factors'
       stop
    end if

    p = 1
    nf = 0
    do j=1,n2
       nf = nf + 1
       fac(nf) = 2
       p = p*2
    end do
    do j=1,n3
       nf = nf + 1
       fac(nf) = 3
       p = p*3
    end do
    do j=1,n4
       nf = nf + 1
       fac(nf) = 4
       p = p*4
    end do
    do j=1,n5
       nf = nf + 1
       fac(nf) = 5
       p = p*5
    end do
    do j=1,n6
       nf = nf + 1
       fac(nf) = 6
       p = p*6
    end do
    do j=1,n8
       nf = nf + 1
       fac(nf) = 8
       p = p*8
    end do

    if (p .ne. n) then
       write(6,*) p, ' is not equal to n = ', n
       stop
    end if

    return
  end subroutine fact235


  subroutine genw(n,nf,fac,w)
    implicit none
    integer,intent(in):: n,nf
    integer:: fac(nfmax)
    real(kind=8),intent(out):: w(2,n)
    !
    real(kind=8),parameter:: tpi=2.0d0*3.14159265358979323846d0
    integer:: iw,darg,ir,rad,narg,k,p
    real(kind=8):: argh,argi,arg
    !
    iw = 1
    darg = n
    do ir = 1,nf
       argh = tpi/darg
       rad = fac(ir)
       narg = darg/rad
       do p = 0,narg-1
          do k = 1,rad-1
             argi = argh*k
             arg = argi*p
             w(1,iw) = cos(arg)
             w(2,iw) = -sin(arg)
             iw = iw + 1
          end do
       end do
       darg = darg/rad
    end do
    return
  end subroutine genw

  subroutine fftrd2(nvk,np,w,xr,xi,yr,yi)
    implicit none
    integer,intent(in):: nvk,np
    real(kind=8),intent(in):: xr(nvk,np,2),xi(nvk,np,2)
    real(kind=8),intent(out):: yr(nvk,2,np),yi(nvk,2,np)
    real(kind=8),intent(in):: w(2,np)
    !
    real(kind=8):: wr1,wi1,tr,ti
    integer:: p,ivk
    !
    do p = 1,np
       wr1 = w(1,p)
       wi1 = w(2,p)
       do ivk = 1,nvk
          yr(ivk,1,p) = xr(ivk,p,1) + xr(ivk,p,2)
          yi(ivk,1,p) = xi(ivk,p,1) + xi(ivk,p,2)
          !
          tr = xr(ivk,p,1) - xr(ivk,p,2)
          ti = xi(ivk,p,1) - xi(ivk,p,2)
          yr(ivk,2,p) = wr1*tr - wi1*ti
          yi(ivk,2,p) = wr1*ti + wi1*tr
       end do
    end do
    return
  end subroutine fftrd2
  !
  subroutine bftrd2(nvk,np,w,xr,xi,yr,yi)
    implicit none
    integer,intent(in):: nvk,np
    real(kind=8),intent(in):: xr(nvk,np,2),xi(nvk,np,2)
    real(kind=8),intent(out):: yr(nvk,2,np),yi(nvk,2,np)
    real(kind=8),intent(in):: w(2,np)
    !
    real(kind=8):: wr1,wi1,tr,ti
    integer:: p,ivk
    !
    do p = 1,np
       wr1 = w(1,p)
       wi1 = w(2,p)
       do ivk = 1,nvk
          yr(ivk,1,p) = xr(ivk,p,1) + xr(ivk,p,2)
          yi(ivk,1,p) = xi(ivk,p,1) + xi(ivk,p,2)
          !
          tr = xr(ivk,p,1) - xr(ivk,p,2)
          ti = xi(ivk,p,1) - xi(ivk,p,2)
          yr(ivk,2,p) = wr1*tr + wi1*ti
          yi(ivk,2,p) = wr1*ti - wi1*tr
       end do
    end do
    return
  end subroutine bftrd2
  !
  subroutine fftrd3(nvk,np,w,xr,xi,yr,yi)
    implicit none
    integer,intent(in):: nvk,np
    real(kind=8),intent(in):: xr(nvk,np,3),xi(nvk,np,3)
    real(kind=8),intent(out):: yr(nvk,3,np),yi(nvk,3,np)
    real(kind=8),intent(in):: w(2,2*np)
    !
    real(kind=8),parameter:: sp3 = 0.866025403784438647d0
    !
    real(kind=8):: wr1,wi1,wr2,wi2,tr1,ti1,tr2,ti2,tr3,ti3,tr,ti
    integer:: p,ivk,iw
    !
    iw = 1
    do p = 1,np
       wr1 = w(1,iw)
       wi1 = w(2,iw)
       iw = iw+1
       wr2 = w(1,iw)
       wi2 = w(2,iw)
       iw = iw+1
       do ivk = 1,nvk
          tr1 = xr(ivk,p,2) + xr(ivk,p,3)
          ti1 = xi(ivk,p,2) + xi(ivk,p,3)
          tr2 = xr(ivk,p,1) - 0.5d0*tr1
          ti2 = xi(ivk,p,1) - 0.5d0*ti1
          tr3 = sp3*(xi(ivk,p,2) - xi(ivk,p,3))
          ti3 = -sp3*(xr(ivk,p,2) - xr(ivk,p,3))
          yr(ivk,1,p) = xr(ivk,p,1) + tr1
          yi(ivk,1,p) = xi(ivk,p,1) + ti1
          tr = tr2 + tr3
          ti = ti2 + ti3
          yr(ivk,2,p) = wr1*tr - wi1*ti
          yi(ivk,2,p) = wr1*ti + wi1*tr
          tr = tr2 - tr3
          ti = ti2 - ti3
          yr(ivk,3,p) = wr2*tr - wi2*ti
          yi(ivk,3,p) = wr2*ti + wi2*tr
       end do
    end do
    return
  end subroutine fftrd3
  !
  subroutine bftrd3(nvk,np,w,xr,xi,yr,yi)
    implicit none
    integer,intent(in):: nvk,np
    real(kind=8),intent(in):: xr(nvk,np,3),xi(nvk,np,3)
    real(kind=8),intent(out):: yr(nvk,3,np),yi(nvk,3,np)
    real(kind=8),intent(in):: w(2,2*np)
    !
    real(kind=8),parameter:: sp3 = 0.866025403784438647d0
    !
    real(kind=8):: wr1,wi1,wr2,wi2,tr1,ti1,tr2,ti2,tr3,ti3,tr,ti
    integer:: p,ivk,iw
    !
    iw = 1
    do p = 1,np
       wr1 = w(1,iw)
       wi1 = w(2,iw)
       iw = iw+1
       wr2 = w(1,iw)
       wi2 = w(2,iw)
       iw = iw+1
       do ivk = 1,nvk
          tr1 = xr(ivk,p,2) + xr(ivk,p,3)
          ti1 = xi(ivk,p,2) + xi(ivk,p,3)
          tr2 = xr(ivk,p,1) - 0.5d0*tr1
          ti2 = xi(ivk,p,1) - 0.5d0*ti1
          tr3 = -sp3*(xi(ivk,p,2) - xi(ivk,p,3))
          ti3 = sp3*(xr(ivk,p,2) - xr(ivk,p,3))
          yr(ivk,1,p) = xr(ivk,p,1) + tr1
          yi(ivk,1,p) = xi(ivk,p,1) + ti1
          tr = tr2 + tr3
          ti = ti2 + ti3
          yr(ivk,2,p) = wr1*tr + wi1*ti
          yi(ivk,2,p) = wr1*ti - wi1*tr
          tr = tr2 - tr3
          ti = ti2 - ti3
          yr(ivk,3,p) = wr2*tr + wi2*ti
          yi(ivk,3,p) = wr2*ti - wi2*tr
       end do
    end do
    return
  end subroutine bftrd3
  !
  subroutine fftrd4(nvk,np,w,xr,xi,yr,yi)
    implicit none
    integer,intent(in):: nvk,np
    real(kind=8),intent(in):: xr(nvk,np,4),xi(nvk,np,4)
    real(kind=8),intent(out):: yr(nvk,4,np),yi(nvk,4,np)
    real(kind=8),intent(in):: w(2,3*np)
    !
    real(kind=8):: wr1,wi1,wr2,wi2,wr3,wi3
    real(kind=8):: tr1,ti1,tr2,ti2,tr3,ti3,tr4,ti4,tr,ti
    integer:: p,ivk,iw
    !
    iw = 1
    do p = 1,np
       wr1 = w(1,iw)
       wi1 = w(2,iw)
       iw = iw+1
       wr2 = w(1,iw)
       wi2 = w(2,iw)
       iw = iw+1
       wr3 = w(1,iw)
       wi3 = w(2,iw)
       iw = iw+1
       do ivk = 1,nvk
          tr1 = xr(ivk,p,1) + xr(ivk,p,3)
          ti1 = xi(ivk,p,1) + xi(ivk,p,3)
          tr2 = xr(ivk,p,1) - xr(ivk,p,3)
          ti2 = xi(ivk,p,1) - xi(ivk,p,3)
          tr3 = xr(ivk,p,2) + xr(ivk,p,4)
          ti3 = xi(ivk,p,2) + xi(ivk,p,4)
          tr4 = xi(ivk,p,2) - xi(ivk,p,4)
          ti4 = xr(ivk,p,4) - xr(ivk,p,2)
          yr(ivk,1,p) = tr1 + tr3
          yi(ivk,1,p) = ti1 + ti3
          tr = tr2 + tr4
          ti = ti2 + ti4
          yr(ivk,2,p) = wr1*tr - wi1*ti
          yi(ivk,2,p) = wr1*ti + wi1*tr
          tr = tr1 - tr3
          ti = ti1 - ti3
          yr(ivk,3,p) = wr2*tr - wi2*ti
          yi(ivk,3,p) = wr2*ti + wi2*tr
          tr = tr2 - tr4
          ti = ti2 - ti4
          yr(ivk,4,p) = wr3*tr - wi3*ti
          yi(ivk,4,p) = wr3*ti + wi3*tr
       end do
    end do
    return
  end subroutine fftrd4
  !
  subroutine bftrd4(nvk,np,w,xr,xi,yr,yi)
    implicit none
    integer,intent(in):: nvk,np
    real(kind=8),intent(in):: xr(nvk,np,4),xi(nvk,np,4)
    real(kind=8),intent(out):: yr(nvk,4,np),yi(nvk,4,np)
    real(kind=8),intent(in):: w(2,3*np)
    !
    real(kind=8):: wr1,wi1,wr2,wi2,wr3,wi3
    real(kind=8):: tr1,ti1,tr2,ti2,tr3,ti3,tr4,ti4,tr,ti
    integer:: p,ivk,iw
    !
    iw = 1
    do p = 1,np
       wr1 = w(1,iw)
       wi1 = w(2,iw)
       iw = iw+1
       wr2 = w(1,iw)
       wi2 = w(2,iw)
       iw = iw+1
       wr3 = w(1,iw)
       wi3 = w(2,iw)
       iw = iw+1
       do ivk = 1,nvk
          tr1 = xr(ivk,p,1) + xr(ivk,p,3)
          ti1 = xi(ivk,p,1) + xi(ivk,p,3)
          tr2 = xr(ivk,p,1) - xr(ivk,p,3)
          ti2 = xi(ivk,p,1) - xi(ivk,p,3)
          tr3 = xr(ivk,p,2) + xr(ivk,p,4)
          ti3 = xi(ivk,p,2) + xi(ivk,p,4)
          tr4 = xi(ivk,p,4) - xi(ivk,p,2)
          ti4 = xr(ivk,p,2) - xr(ivk,p,4)
          yr(ivk,1,p) = tr1 + tr3
          yi(ivk,1,p) = ti1 + ti3
          tr = tr2 + tr4
          ti = ti2 + ti4
          yr(ivk,2,p) = wr1*tr + wi1*ti
          yi(ivk,2,p) = wr1*ti - wi1*tr
          tr = tr1 - tr3
          ti = ti1 - ti3
          yr(ivk,3,p) = wr2*tr + wi2*ti
          yi(ivk,3,p) = wr2*ti - wi2*tr
          tr = tr2 - tr4
          ti = ti2 - ti4
          yr(ivk,4,p) = wr3*tr + wi3*ti
          yi(ivk,4,p) = wr3*ti - wi3*tr
       end do
    end do
    return
  end subroutine bftrd4
  !
  subroutine fftrd5(nvk,np,w,xr,xi,yr,yi)
    implicit none
    integer,intent(in):: nvk,np
    real(kind=8),intent(in):: xr(nvk,np,5),xi(nvk,np,5)
    real(kind=8),intent(out):: yr(nvk,5,np),yi(nvk,5,np)
    real(kind=8),intent(in):: w(2,4*np)
    !
    real(kind=8),parameter:: sp25 = 0.951056516295153572d0
    real(kind=8),parameter:: sq54 = 0.559016994374947424d0
    real(kind=8),parameter:: ss = 0.618033988749894848d0
    !
    real(kind=8):: wr1,wi1,wr2,wi2,wr3,wi3,wr4,wi4
    real(kind=8):: tr1,ti1,tr2,ti2,tr3,ti3,tr4,ti4,tr5,ti5
    real(kind=8):: tr6,ti6,tr7,ti7,tr8,ti8,tr9,ti9,tr10,ti10,tr11,ti11
    real(kind=8):: tr,ti
    integer:: p,ivk,iw
    !
    do ivk = 1,nvk
       tr1 = xr(ivk,1,2) + xr(ivk,1,5)
       ti1 = xi(ivk,1,2) + xi(ivk,1,5)
       tr2 = xr(ivk,1,3) + xr(ivk,1,4)
       ti2 = xi(ivk,1,3) + xi(ivk,1,4)
       tr3 = sp25*(xr(ivk,1,2) - xr(ivk,1,5))
       ti3 = sp25*(xi(ivk,1,2) - xi(ivk,1,5))
       tr4 = sp25*(xr(ivk,1,3) - xr(ivk,1,4))
       ti4 = sp25*(xi(ivk,1,3) - xi(ivk,1,4))
       tr5 = tr1 + tr2
       ti5 = ti1 + ti2
       tr6 = sq54*(tr1 - tr2)
       ti6 = sq54*(ti1 - ti2)
       tr7 = xr(ivk,1,1) - 0.25d0*tr5
       ti7 = xi(ivk,1,1) - 0.25d0*ti5
       tr8 = tr7 + tr6
       ti8 = ti7 + ti6
       tr9 = tr7 - tr6
       ti9 = ti7 - ti6
       tr10 = ti3 + ss*ti4
       ti10 = - tr3 - ss*tr4
       tr11 = ss*ti3 - ti4
       ti11 = - ss*tr3 + tr4
       yr(ivk,1,1) = xr(ivk,1,1) + tr5
       yi(ivk,1,1) = xi(ivk,1,1) + ti5
       yr(ivk,2,1) = tr8 + tr10
       yi(ivk,2,1) = ti8 + ti10
       yr(ivk,3,1) = tr9 + tr11
       yi(ivk,3,1) = ti9 + ti11
       yr(ivk,4,1) = tr9 - tr11
       yi(ivk,4,1) = ti9 - ti11
       yr(ivk,5,1) = tr8 - tr10
       yi(ivk,5,1) = ti8 - ti10
    end do

    iw = 5
    do p = 2,np
       wr1 = w(1,iw)
       wi1 = w(2,iw)
       iw = iw+1
       wr2 = w(1,iw)
       wi2 = w(2,iw)
       iw = iw+1
       wr3 = w(1,iw)
       wi3 = w(2,iw)
       iw = iw+1
       wr4 = w(1,iw)
       wi4 = w(2,iw)
       iw = iw+1
       do ivk = 1,nvk
          tr1 = xr(ivk,p,2) + xr(ivk,p,5)
          ti1 = xi(ivk,p,2) + xi(ivk,p,5)
          tr2 = xr(ivk,p,3) + xr(ivk,p,4)
          ti2 = xi(ivk,p,3) + xi(ivk,p,4)
          tr3 = sp25*(xr(ivk,p,2) - xr(ivk,p,5))
          ti3 = sp25*(xi(ivk,p,2) - xi(ivk,p,5))
          tr4 = sp25*(xr(ivk,p,3) - xr(ivk,p,4))
          ti4 = sp25*(xi(ivk,p,3) - xi(ivk,p,4))
          tr5 = tr1 + tr2
          ti5 = ti1 + ti2
          tr6 = sq54*(tr1 - tr2)
          ti6 = sq54*(ti1 - ti2)
          tr7 = xr(ivk,p,1) - 0.25d0*tr5
          ti7 = xi(ivk,p,1) - 0.25d0*ti5
          tr8 = tr7 + tr6
          ti8 = ti7 + ti6
          tr9 = tr7 - tr6
          ti9 = ti7 - ti6
          tr10 = ti3 + ss*ti4
          ti10 = - tr3 - ss*tr4
          tr11 = ss*ti3 - ti4
          ti11 = - ss*tr3 + tr4
          yr(ivk,1,p) = xr(ivk,p,1) + tr5
          yi(ivk,1,p) = xi(ivk,p,1) + ti5
          tr = tr8 + tr10
          ti = ti8 + ti10
          yr(ivk,2,p) = wr1*tr - wi1*ti
          yi(ivk,2,p) = wr1*ti + wi1*tr
          tr = tr9 + tr11
          ti = ti9 + ti11
          yr(ivk,3,p) = wr2*tr - wi2*ti
          yi(ivk,3,p) = wr2*ti + wi2*tr
          tr = tr9 - tr11
          ti = ti9 - ti11
          yr(ivk,4,p) = wr3*tr - wi3*ti
          yi(ivk,4,p) = wr3*ti + wi3*tr
          tr = tr8 - tr10
          ti = ti8 - ti10
          yr(ivk,5,p) = wr4*tr - wi4*ti
          yi(ivk,5,p) = wr4*ti + wi4*tr
       end do
    end do
    return
  end subroutine fftrd5
  !
  subroutine bftrd5(nvk,np,w,xr,xi,yr,yi)
    implicit none
    integer,intent(in):: nvk,np
    real(kind=8),intent(in):: xr(nvk,np,5),xi(nvk,np,5)
    real(kind=8),intent(out):: yr(nvk,5,np),yi(nvk,5,np)
    real(kind=8),intent(in):: w(2,4*np)
    !
    real(kind=8),parameter:: sp25 = 0.951056516295153572d0
    real(kind=8),parameter:: sq54 = 0.559016994374947424d0
    real(kind=8),parameter:: ss = 0.618033988749894848d0
    !
    real(kind=8):: wr1,wi1,wr2,wi2,wr3,wi3,wr4,wi4
    real(kind=8):: tr1,ti1,tr2,ti2,tr3,ti3,tr4,ti4,tr5,ti5
    real(kind=8):: tr6,ti6,tr7,ti7,tr8,ti8,tr9,ti9,tr10,ti10,tr11,ti11
    real(kind=8):: tr,ti
    integer:: p,ivk,iw
    !
    do ivk = 1,nvk
       tr1 = xr(ivk,1,2) + xr(ivk,1,5)
       ti1 = xi(ivk,1,2) + xi(ivk,1,5)
       tr2 = xr(ivk,1,3) + xr(ivk,1,4)
       ti2 = xi(ivk,1,3) + xi(ivk,1,4)
       tr3 = sp25*(xr(ivk,1,2) - xr(ivk,1,5))
       ti3 = sp25*(xi(ivk,1,2) - xi(ivk,1,5))
       tr4 = sp25*(xr(ivk,1,3) - xr(ivk,1,4))
       ti4 = sp25*(xi(ivk,1,3) - xi(ivk,1,4))
       tr5 = tr1 + tr2
       ti5 = ti1 + ti2
       tr6 = sq54*(tr1 - tr2)
       ti6 = sq54*(ti1 - ti2)
       tr7 = xr(ivk,1,1) - 0.25d0*tr5
       ti7 = xi(ivk,1,1) - 0.25d0*ti5
       tr8 = tr7 + tr6
       ti8 = ti7 + ti6
       tr9 = tr7 - tr6
       ti9 = ti7 - ti6
       tr10 = - ti3 - ss*ti4
       ti10 = tr3 + ss*tr4
       tr11 = - ss*ti3 + ti4
       ti11 = ss*tr3 - tr4
       yr(ivk,1,1) = xr(ivk,1,1) + tr5
       yi(ivk,1,1) = xi(ivk,1,1) + ti5
       yr(ivk,2,1) = tr8 + tr10
       yi(ivk,2,1) = ti8 + ti10
       yr(ivk,3,1) = tr9 + tr11
       yi(ivk,3,1) = ti9 + ti11
       yr(ivk,4,1) = tr9 - tr11
       yi(ivk,4,1) = ti9 - ti11
       yr(ivk,5,1) = tr8 - tr10
       yi(ivk,5,1) = ti8 - ti10
    end do

    iw = 5
    do p = 2,np
       wr1 = w(1,iw)
       wi1 = w(2,iw)
       iw = iw+1
       wr2 = w(1,iw)
       wi2 = w(2,iw)
       iw = iw+1
       wr3 = w(1,iw)
       wi3 = w(2,iw)
       iw = iw+1
       wr4 = w(1,iw)
       wi4 = w(2,iw)
       iw = iw+1
       do ivk = 1,nvk
          tr1 = xr(ivk,p,2) + xr(ivk,p,5)
          ti1 = xi(ivk,p,2) + xi(ivk,p,5)
          tr2 = xr(ivk,p,3) + xr(ivk,p,4)
          ti2 = xi(ivk,p,3) + xi(ivk,p,4)
          tr3 = sp25*(xr(ivk,p,2) - xr(ivk,p,5))
          ti3 = sp25*(xi(ivk,p,2) - xi(ivk,p,5))
          tr4 = sp25*(xr(ivk,p,3) - xr(ivk,p,4))
          ti4 = sp25*(xi(ivk,p,3) - xi(ivk,p,4))
          tr5 = tr1 + tr2
          ti5 = ti1 + ti2
          tr6 = sq54*(tr1 - tr2)
          ti6 = sq54*(ti1 - ti2)
          tr7 = xr(ivk,p,1) - 0.25d0*tr5
          ti7 = xi(ivk,p,1) - 0.25d0*ti5
          tr8 = tr7 + tr6
          ti8 = ti7 + ti6
          tr9 = tr7 - tr6
          ti9 = ti7 - ti6
          tr10 = - ti3 - ss*ti4
          ti10 = tr3 + ss*tr4
          tr11 = - ss*ti3 + ti4
          ti11 = ss*tr3 - tr4
          yr(ivk,1,p) = xr(ivk,p,1) + tr5
          yi(ivk,1,p) = xi(ivk,p,1) + ti5
          tr = tr8 + tr10
          ti = ti8 + ti10
          yr(ivk,2,p) = wr1*tr + wi1*ti
          yi(ivk,2,p) = wr1*ti - wi1*tr
          tr = tr9 + tr11
          ti = ti9 + ti11
          yr(ivk,3,p) = wr2*tr + wi2*ti
          yi(ivk,3,p) = wr2*ti - wi2*tr
          tr = tr9 - tr11
          ti = ti9 - ti11
          yr(ivk,4,p) = wr3*tr + wi3*ti
          yi(ivk,4,p) = wr3*ti - wi3*tr
          tr = tr8 - tr10
          ti = ti8 - ti10
          yr(ivk,5,p) = wr4*tr + wi4*ti
          yi(ivk,5,p) = wr4*ti - wi4*tr
       end do
    end do
    return
  end subroutine bftrd5

  subroutine fftrd6(nvk,np,w,xr,xi,yr,yi)
    implicit none
    integer,intent(in):: nvk,np
    real(kind=8),intent(in):: xr(nvk,np,6),xi(nvk,np,6)
    real(kind=8),intent(out):: yr(nvk,6,np),yi(nvk,6,np)
    real(kind=8),intent(in):: w(2,5*np)
    !
    real(kind=8),parameter:: sp3 = 0.866025403784438647d0
    !
    real(kind=8):: wr1,wi1,wr2,wi2,wr3,wi3,wr4,wi4,wr5,wi5
    real(kind=8):: tr1,ti1,tr2,ti2,tr3,ti3,tr,ti
    real(kind=8):: ur1,ui1,ur2,ui2,ur3,ui3
    real(kind=8):: vr1,vi1,vr2,vi2,vr3,vi3
    integer:: p,ivk,iw
    !
    do ivk = 1,nvk
       tr1 = xr(ivk,1,3) + xr(ivk,1,5)
       ti1 = xi(ivk,1,3) + xi(ivk,1,5)
       tr2 = xr(ivk,1,1) - 0.5d0*tr1
       ti2 = xi(ivk,1,1) - 0.5d0*ti1
       tr3 = -sp3*(xr(ivk,1,3) - xr(ivk,1,5))
       ti3 = -sp3*(xi(ivk,1,3) - xi(ivk,1,5))
       ur1 = xr(ivk,1,1) + tr1
       ui1 = xi(ivk,1,1) + ti1
       ur2 = tr2 - ti3
       ui2 = ti2 + tr3
       ur3 = tr2 + ti3
       ui3 = ti2 - tr3

       tr1 = xr(ivk,1,6) + xr(ivk,1,2)
       ti1 = xi(ivk,1,6) + xi(ivk,1,2)
       tr2 = xr(ivk,1,4) - 0.5d0*tr1
       ti2 = xi(ivk,1,4) - 0.5d0*ti1
       tr3 = -sp3*(xr(ivk,1,6) - xr(ivk,1,2))
       ti3 = -sp3*(xi(ivk,1,6) - xi(ivk,1,2))
       vr1 = xr(ivk,1,4) + tr1
       vi1 = xi(ivk,1,4) + ti1
       vr2 = tr2 - ti3
       vi2 = ti2 + tr3
       vr3 = tr2 + ti3
       vi3 = ti2 - tr3

       yr(ivk,1,1) = ur1 + vr1
       yi(ivk,1,1) = ui1 + vi1
       yr(ivk,5,1) = ur2 + vr2
       yi(ivk,5,1) = ui2 + vi2
       yr(ivk,3,1) = ur3 + vr3
       yi(ivk,3,1) = ui3 + vi3
       yr(ivk,4,1) = ur1 - vr1
       yi(ivk,4,1) = ui1 - vi1
       yr(ivk,2,1) = ur2 - vr2
       yi(ivk,2,1) = ui2 - vi2
       yr(ivk,6,1) = ur3 - vr3
       yi(ivk,6,1) = ui3 - vi3
    end do

    iw = 6
    do p = 2,np
       wr1 = w(1,iw)
       wi1 = w(2,iw)
       iw = iw+1
       wr2 = w(1,iw)
       wi2 = w(2,iw)
       iw = iw+1
       wr3 = w(1,iw)
       wi3 = w(2,iw)
       iw = iw+1
       wr4 = w(1,iw)
       wi4 = w(2,iw)
       iw = iw+1
       wr5 = w(1,iw)
       wi5 = w(2,iw)
       iw = iw+1
       do ivk = 1,nvk
          tr1 = xr(ivk,p,3) + xr(ivk,p,5)
          ti1 = xi(ivk,p,3) + xi(ivk,p,5)
          tr2 = xr(ivk,p,1) - 0.5d0*tr1
          ti2 = xi(ivk,p,1) - 0.5d0*ti1
          tr3 = -sp3*(xr(ivk,p,3) - xr(ivk,p,5))
          ti3 = -sp3*(xi(ivk,p,3) - xi(ivk,p,5))
          ur1 = xr(ivk,p,1) + tr1
          ui1 = xi(ivk,p,1) + ti1
          ur2 = tr2 - ti3
          ui2 = ti2 + tr3
          ur3 = tr2 + ti3
          ui3 = ti2 - tr3

          tr1 = xr(ivk,p,6) + xr(ivk,p,2)
          ti1 = xi(ivk,p,6) + xi(ivk,p,2)
          tr2 = xr(ivk,p,4) - 0.5d0*tr1
          ti2 = xi(ivk,p,4) - 0.5d0*ti1
          tr3 = -sp3*(xr(ivk,p,6) - xr(ivk,p,2))
          ti3 = -sp3*(xi(ivk,p,6) - xi(ivk,p,2))
          vr1 = xr(ivk,p,4) + tr1
          vi1 = xi(ivk,p,4) + ti1
          vr2 = tr2 - ti3
          vi2 = ti2 + tr3
          vr3 = tr2 + ti3
          vi3 = ti2 - tr3

          yr(ivk,1,p) = ur1 + vr1
          yi(ivk,1,p) = ui1 + vi1
          tr = ur2 + vr2
          ti = ui2 + vi2
          yr(ivk,5,p) = wr4*tr - wi4*ti
          yi(ivk,5,p) = wr4*ti + wi4*tr
          tr = ur3 + vr3
          ti = ui3 + vi3
          yr(ivk,3,p) = wr2*tr - wi2*ti
          yi(ivk,3,p) = wr2*ti + wi2*tr
          tr = ur1 - vr1
          ti = ui1 - vi1
          yr(ivk,4,p) = wr3*tr - wi3*ti
          yi(ivk,4,p) = wr3*ti + wi3*tr
          tr = ur2 - vr2
          ti = ui2 - vi2
          yr(ivk,2,p) = wr1*tr - wi1*ti
          yi(ivk,2,p) = wr1*ti + wi1*tr
          tr = ur3 - vr3
          ti = ui3 - vi3
          yr(ivk,6,p) = wr5*tr - wi5*ti
          yi(ivk,6,p) = wr5*ti + wi5*tr
       end do
    end do
    return
  end subroutine fftrd6

  subroutine bftrd6(nvk,np,w,xr,xi,yr,yi)
    implicit none
    integer,intent(in):: nvk,np
    real(kind=8),intent(in):: xr(nvk,np,6),xi(nvk,np,6)
    real(kind=8),intent(out):: yr(nvk,6,np),yi(nvk,6,np)
    real(kind=8),intent(in):: w(2,5*np)
    !
    real(kind=8),parameter:: sp3 = 0.866025403784438647d0
    !
    real(kind=8):: wr1,wi1,wr2,wi2,wr3,wi3,wr4,wi4,wr5,wi5
    real(kind=8):: tr1,ti1,tr2,ti2,tr3,ti3,tr,ti
    real(kind=8):: ur1,ui1,ur2,ui2,ur3,ui3
    real(kind=8):: vr1,vi1,vr2,vi2,vr3,vi3
    integer:: p,ivk,iw
    !
    do ivk = 1,nvk
       tr1 = xr(ivk,1,3) + xr(ivk,1,5)
       ti1 = xi(ivk,1,3) + xi(ivk,1,5)
       tr2 = xr(ivk,1,1) - 0.5d0*tr1
       ti2 = xi(ivk,1,1) - 0.5d0*ti1
       tr3 = sp3*(xr(ivk,1,3) - xr(ivk,1,5))
       ti3 = sp3*(xi(ivk,1,3) - xi(ivk,1,5))
       ur1 = xr(ivk,1,1) + tr1
       ui1 = xi(ivk,1,1) + ti1
       ur2 = tr2 - ti3
       ui2 = ti2 + tr3
       ur3 = tr2 + ti3
       ui3 = ti2 - tr3

       tr1 = xr(ivk,1,6) + xr(ivk,1,2)
       ti1 = xi(ivk,1,6) + xi(ivk,1,2)
       tr2 = xr(ivk,1,4) - 0.5d0*tr1
       ti2 = xi(ivk,1,4) - 0.5d0*ti1
       tr3 = sp3*(xr(ivk,1,6) - xr(ivk,1,2))
       ti3 = sp3*(xi(ivk,1,6) - xi(ivk,1,2))
       vr1 = xr(ivk,1,4) + tr1
       vi1 = xi(ivk,1,4) + ti1
       vr2 = tr2 - ti3
       vi2 = ti2 + tr3
       vr3 = tr2 + ti3
       vi3 = ti2 - tr3

       yr(ivk,1,1) = ur1 + vr1
       yi(ivk,1,1) = ui1 + vi1
       yr(ivk,5,1) = ur2 + vr2
       yi(ivk,5,1) = ui2 + vi2
       yr(ivk,3,1) = ur3 + vr3
       yi(ivk,3,1) = ui3 + vi3
       yr(ivk,4,1) = ur1 - vr1
       yi(ivk,4,1) = ui1 - vi1
       yr(ivk,2,1) = ur2 - vr2
       yi(ivk,2,1) = ui2 - vi2
       yr(ivk,6,1) = ur3 - vr3
       yi(ivk,6,1) = ui3 - vi3
    end do

!   write(6,*) '20150826 np=',np 
    iw = 6
    do p = 2,np
       wr1 = w(1,iw)
       wi1 = w(2,iw)
       iw = iw+1
       wr2 = w(1,iw)
       wi2 = w(2,iw)
       iw = iw+1
       wr3 = w(1,iw)
       wi3 = w(2,iw)
       iw = iw+1
       wr4 = w(1,iw)
       wi4 = w(2,iw)
       iw = iw+1
       wr5 = w(1,iw)
       wi5 = w(2,iw)
       iw = iw+1
       do ivk = 1,nvk
          tr1 = xr(ivk,p,3) + xr(ivk,p,5)
          ti1 = xi(ivk,p,3) + xi(ivk,p,5)
          tr2 = xr(ivk,p,1) - 0.5d0*tr1
          ti2 = xi(ivk,p,1) - 0.5d0*ti1
          tr3 = sp3*(xr(ivk,p,3) - xr(ivk,p,5))
          ti3 = sp3*(xi(ivk,p,3) - xi(ivk,p,5))
          ur1 = xr(ivk,p,1) + tr1
          ui1 = xi(ivk,p,1) + ti1
          ur2 = tr2 - ti3
          ui2 = ti2 + tr3
          ur3 = tr2 + ti3
          ui3 = ti2 - tr3

          tr1 = xr(ivk,p,6) + xr(ivk,p,2)
          ti1 = xi(ivk,p,6) + xi(ivk,p,2)
          tr2 = xr(ivk,p,4) - 0.5d0*tr1
          ti2 = xi(ivk,p,4) - 0.5d0*ti1
          tr3 = sp3*(xr(ivk,p,6) - xr(ivk,p,2))
          ti3 = sp3*(xi(ivk,p,6) - xi(ivk,p,2))
          vr1 = xr(ivk,p,4) + tr1
          vi1 = xi(ivk,p,4) + ti1
          vr2 = tr2 - ti3
          vi2 = ti2 + tr3
          vr3 = tr2 + ti3
          vi3 = ti2 - tr3

          yr(ivk,1,p) = ur1 + vr1
          yi(ivk,1,p) = ui1 + vi1
          tr = ur2 + vr2
          ti = ui2 + vi2
          yr(ivk,5,p) = wr4*tr + wi4*ti
          yi(ivk,5,p) = wr4*ti - wi4*tr
          tr = ur3 + vr3
          ti = ui3 + vi3
          yr(ivk,3,p) = wr2*tr + wi2*ti
          yi(ivk,3,p) = wr2*ti - wi2*tr
          tr = ur1 - vr1
          ti = ui1 - vi1
          yr(ivk,4,p) = wr3*tr + wi3*ti
          yi(ivk,4,p) = wr3*ti - wi3*tr
          tr = ur2 - vr2
          ti = ui2 - vi2
          yr(ivk,2,p) = wr1*tr + wi1*ti
          yi(ivk,2,p) = wr1*ti - wi1*tr
          tr = ur3 - vr3
          ti = ui3 - vi3
          yr(ivk,6,p) = wr5*tr + wi5*ti
          yi(ivk,6,p) = wr5*ti - wi5*tr
       end do
    end do
    return
  end subroutine bftrd6
  !
  subroutine fftrd8(nvk,np,w,xr,xi,yr,yi)
    implicit none
    integer,intent(in):: nvk,np
    real(kind=8),intent(in):: xr(nvk,np,8),xi(nvk,np,8)
    real(kind=8),intent(out):: yr(nvk,8,np),yi(nvk,8,np)
    real(kind=8),intent(in):: w(2,7*np)
    !
    real(kind=8),parameter:: sp2 = 0.70710678118654752440d0
    !
    real(kind=8):: wr1,wi1,wr2,wi2,wr3,wi3,wr4,wi4,wr5,wi5,wr6,wi6,wr7,wi7
    real(kind=8):: tr1,ti1,tr2,ti2,tr3,ti3,tr4,ti4,tr,ti,sr,si
    real(kind=8):: ur1,ui1,ur2,ui2,ur3,ui3,ur4,ui4
    real(kind=8):: vr1,vi1,vr2,vi2,vr3,vi3,vr4,vi4
    integer:: p,ivk,iw
    !
    do ivk = 1,nvk
       tr1 = xr(ivk,1,1) + xr(ivk,1,5)
       ti1 = xi(ivk,1,1) + xi(ivk,1,5)
       tr2 = xr(ivk,1,1) - xr(ivk,1,5)
       ti2 = xi(ivk,1,1) - xi(ivk,1,5)
       tr3 = xr(ivk,1,3) + xr(ivk,1,7)
       ti3 = xi(ivk,1,3) + xi(ivk,1,7)
       tr4 = xi(ivk,1,3) - xi(ivk,1,7)
       ti4 = xr(ivk,1,7) - xr(ivk,1,3)
       ur1 = tr1 + tr3
       ui1 = ti1 + ti3
       ur2 = tr2 + tr4
       ui2 = ti2 + ti4
       ur3 = tr1 - tr3
       ui3 = ti1 - ti3
       ur4 = tr2 - tr4
       ui4 = ti2 - ti4
       tr1 = xr(ivk,1,2) + xr(ivk,1,6)
       ti1 = xi(ivk,1,2) + xi(ivk,1,6)
       tr2 = xr(ivk,1,2) - xr(ivk,1,6)
       ti2 = xi(ivk,1,2) - xi(ivk,1,6)
       tr3 = xr(ivk,1,4) + xr(ivk,1,8)
       ti3 = xi(ivk,1,4) + xi(ivk,1,8)
       tr4 = xi(ivk,1,4) - xi(ivk,1,8)
       ti4 = xr(ivk,1,8) - xr(ivk,1,4)
       vr1 = tr1 + tr3
       vi1 = ti1 + ti3
       vr2 = tr2 + tr4
       vi2 = ti2 + ti4
       vr3 = tr1 - tr3
       vi3 = ti1 - ti3
       vr4 = tr2 - tr4
       vi4 = ti2 - ti4
       yr(ivk,1,1) = ur1 + vr1
       yi(ivk,1,1) = ui1 + vi1
       yr(ivk,5,1) = ur1 - vr1
       yi(ivk,5,1) = ui1 - vi1

       sr = sp2*(vr2 + vi2)
       si = sp2*(vi2 - vr2)
       yr(ivk,2,1) = ur2 + sr
       yi(ivk,2,1) = ui2 + si
       yr(ivk,6,1) = ur2 - sr
       yi(ivk,6,1) = ui2 - si

       yr(ivk,3,1) = ur3 + vi3
       yi(ivk,3,1) = ui3 - vr3
       yr(ivk,7,1) = ur3 - vi3
       yi(ivk,7,1) = ui3 + vr3

       sr = sp2*(vi4 - vr4)
       si = -sp2*(vr4 + vi4)
       yr(ivk,4,1) = ur4 + sr
       yi(ivk,4,1) = ui4 + si
       yr(ivk,8,1) = ur4 - sr
       yi(ivk,8,1) = ui4 - si
    end do

    iw = 8
    do p = 2,np
       wr1 = w(1,iw)
       wi1 = w(2,iw)
       iw = iw+1
       wr2 = w(1,iw)
       wi2 = w(2,iw)
       iw = iw+1
       wr3 = w(1,iw)
       wi3 = w(2,iw)
       iw = iw+1
       wr4 = w(1,iw)
       wi4 = w(2,iw)
       iw = iw+1
       wr5 = w(1,iw)
       wi5 = w(2,iw)
       iw = iw+1
       wr6 = w(1,iw)
       wi6 = w(2,iw)
       iw = iw+1
       wr7 = w(1,iw)
       wi7 = w(2,iw)
       iw = iw+1
       do ivk = 1,nvk
          tr1 = xr(ivk,p,1) + xr(ivk,p,5)
          ti1 = xi(ivk,p,1) + xi(ivk,p,5)
          tr2 = xr(ivk,p,1) - xr(ivk,p,5)
          ti2 = xi(ivk,p,1) - xi(ivk,p,5)
          tr3 = xr(ivk,p,3) + xr(ivk,p,7)
          ti3 = xi(ivk,p,3) + xi(ivk,p,7)
          tr4 = xi(ivk,p,3) - xi(ivk,p,7)
          ti4 = xr(ivk,p,7) - xr(ivk,p,3)
          ur1 = tr1 + tr3
          ui1 = ti1 + ti3
          ur2 = tr2 + tr4
          ui2 = ti2 + ti4
          ur3 = tr1 - tr3
          ui3 = ti1 - ti3
          ur4 = tr2 - tr4
          ui4 = ti2 - ti4
          tr1 = xr(ivk,p,2) + xr(ivk,p,6)
          ti1 = xi(ivk,p,2) + xi(ivk,p,6)
          tr2 = xr(ivk,p,2) - xr(ivk,p,6)
          ti2 = xi(ivk,p,2) - xi(ivk,p,6)
          tr3 = xr(ivk,p,4) + xr(ivk,p,8)
          ti3 = xi(ivk,p,4) + xi(ivk,p,8)
          tr4 = xi(ivk,p,4) - xi(ivk,p,8)
          ti4 = xr(ivk,p,8) - xr(ivk,p,4)
          vr1 = tr1 + tr3
          vi1 = ti1 + ti3
          vr2 = tr2 + tr4
          vi2 = ti2 + ti4
          vr3 = tr1 - tr3
          vi3 = ti1 - ti3
          vr4 = tr2 - tr4
          vi4 = ti2 - ti4
          yr(ivk,1,p) = ur1 + vr1
          yi(ivk,1,p) = ui1 + vi1
          tr = ur1 - vr1
          ti = ui1 - vi1
          yr(ivk,5,p) = wr4*tr - wi4*ti
          yi(ivk,5,p) = wr4*ti + wi4*tr

          sr = sp2*(vr2 + vi2)
          si = sp2*(vi2 - vr2)
          tr = ur2 + sr
          ti = ui2 + si
          yr(ivk,2,p) = wr1*tr - wi1*ti
          yi(ivk,2,p) = wr1*ti + wi1*tr
          tr = ur2 - sr
          ti = ui2 - si
          yr(ivk,6,p) = wr5*tr - wi5*ti
          yi(ivk,6,p) = wr5*ti + wi5*tr

          tr = ur3 + vi3
          ti = ui3 - vr3
          yr(ivk,3,p) = wr2*tr - wi2*ti
          yi(ivk,3,p) = wr2*ti + wi2*tr
          tr = ur3 - vi3
          ti = ui3 + vr3
          yr(ivk,7,p) = wr6*tr - wi6*ti
          yi(ivk,7,p) = wr6*ti + wi6*tr

          sr = sp2*(vi4 - vr4)
          si = -sp2*(vr4 + vi4)
          tr = ur4 + sr
          ti = ui4 + si
          yr(ivk,4,p) = wr3*tr - wi3*ti
          yi(ivk,4,p) = wr3*ti + wi3*tr
          tr = ur4 - sr
          ti = ui4 - si
          yr(ivk,8,p) = wr7*tr - wi7*ti
          yi(ivk,8,p) = wr7*ti + wi7*tr
       end do
    end do
    return
  end subroutine fftrd8
  !
  subroutine bftrd8(nvk,np,w,xr,xi,yr,yi)
    implicit none
    integer,intent(in):: nvk,np
    real(kind=8),intent(in):: xr(nvk,np,8),xi(nvk,np,8)
    real(kind=8),intent(out):: yr(nvk,8,np),yi(nvk,8,np)
    real(kind=8),intent(in):: w(2,7*np)
    !
    real(kind=8),parameter:: sp2 = 0.70710678118654752440d0
    !
    real(kind=8):: wr1,wi1,wr2,wi2,wr3,wi3,wr4,wi4,wr5,wi5,wr6,wi6,wr7,wi7
    real(kind=8):: tr1,ti1,tr2,ti2,tr3,ti3,tr4,ti4,tr,ti,sr,si
    real(kind=8):: ur1,ui1,ur2,ui2,ur3,ui3,ur4,ui4
    real(kind=8):: vr1,vi1,vr2,vi2,vr3,vi3,vr4,vi4
    integer:: p,ivk,iw
    !
    do ivk = 1,nvk
       tr1 = xr(ivk,1,1) + xr(ivk,1,5)
       ti1 = xi(ivk,1,1) + xi(ivk,1,5)
       tr2 = xr(ivk,1,1) - xr(ivk,1,5)
       ti2 = xi(ivk,1,1) - xi(ivk,1,5)
       tr3 = xr(ivk,1,3) + xr(ivk,1,7)
       ti3 = xi(ivk,1,3) + xi(ivk,1,7)
       tr4 = xi(ivk,1,7) - xi(ivk,1,3)
       ti4 = xr(ivk,1,3) - xr(ivk,1,7)
       ur1 = tr1 + tr3
       ui1 = ti1 + ti3
       ur2 = tr2 + tr4
       ui2 = ti2 + ti4
       ur3 = tr1 - tr3
       ui3 = ti1 - ti3
       ur4 = tr2 - tr4
       ui4 = ti2 - ti4
       tr1 = xr(ivk,1,2) + xr(ivk,1,6)
       ti1 = xi(ivk,1,2) + xi(ivk,1,6)
       tr2 = xr(ivk,1,2) - xr(ivk,1,6)
       ti2 = xi(ivk,1,2) - xi(ivk,1,6)
       tr3 = xr(ivk,1,4) + xr(ivk,1,8)
       ti3 = xi(ivk,1,4) + xi(ivk,1,8)
       tr4 = xi(ivk,1,8) - xi(ivk,1,4)
       ti4 = xr(ivk,1,4) - xr(ivk,1,8)
       vr1 = tr1 + tr3
       vi1 = ti1 + ti3
       vr2 = tr2 + tr4
       vi2 = ti2 + ti4
       vr3 = tr1 - tr3
       vi3 = ti1 - ti3
       vr4 = tr2 - tr4
       vi4 = ti2 - ti4
       yr(ivk,1,1) = ur1 + vr1
       yi(ivk,1,1) = ui1 + vi1
       yr(ivk,5,1) = ur1 - vr1
       yi(ivk,5,1) = ui1 - vi1

       sr = sp2*(vr2 - vi2)
       si = sp2*(vi2 + vr2)
       yr(ivk,2,1) = ur2 + sr
       yi(ivk,2,1) = ui2 + si
       yr(ivk,6,1) = ur2 - sr
       yi(ivk,6,1) = ui2 - si

       yr(ivk,3,1) = ur3 - vi3
       yi(ivk,3,1) = ui3 + vr3
       yr(ivk,7,1) = ur3 + vi3
       yi(ivk,7,1) = ui3 - vr3

       sr = -sp2*(vr4 + vi4)
       si = sp2*(vr4 - vi4)
       yr(ivk,4,1) = ur4 + sr
       yi(ivk,4,1) = ui4 + si
       yr(ivk,8,1) = ur4 - sr
       yi(ivk,8,1) = ui4 - si
    end do

    iw = 8
    do p = 2,np
       wr1 = w(1,iw)
       wi1 = w(2,iw)
       iw = iw+1
       wr2 = w(1,iw)
       wi2 = w(2,iw)
       iw = iw+1
       wr3 = w(1,iw)
       wi3 = w(2,iw)
       iw = iw+1
       wr4 = w(1,iw)
       wi4 = w(2,iw)
       iw = iw+1
       wr5 = w(1,iw)
       wi5 = w(2,iw)
       iw = iw+1
       wr6 = w(1,iw)
       wi6 = w(2,iw)
       iw = iw+1
       wr7 = w(1,iw)
       wi7 = w(2,iw)
       iw = iw+1
       do ivk = 1,nvk
          tr1 = xr(ivk,p,1) + xr(ivk,p,5)
          ti1 = xi(ivk,p,1) + xi(ivk,p,5)
          tr2 = xr(ivk,p,1) - xr(ivk,p,5)
          ti2 = xi(ivk,p,1) - xi(ivk,p,5)
          tr3 = xr(ivk,p,3) + xr(ivk,p,7)
          ti3 = xi(ivk,p,3) + xi(ivk,p,7)
          tr4 = xi(ivk,p,7) - xi(ivk,p,3)
          ti4 = xr(ivk,p,3) - xr(ivk,p,7)
          ur1 = tr1 + tr3
          ui1 = ti1 + ti3
          ur2 = tr2 + tr4
          ui2 = ti2 + ti4
          ur3 = tr1 - tr3
          ui3 = ti1 - ti3
          ur4 = tr2 - tr4
          ui4 = ti2 - ti4
          tr1 = xr(ivk,p,2) + xr(ivk,p,6)
          ti1 = xi(ivk,p,2) + xi(ivk,p,6)
          tr2 = xr(ivk,p,2) - xr(ivk,p,6)
          ti2 = xi(ivk,p,2) - xi(ivk,p,6)
          tr3 = xr(ivk,p,4) + xr(ivk,p,8)
          ti3 = xi(ivk,p,4) + xi(ivk,p,8)
          tr4 = xi(ivk,p,8) - xi(ivk,p,4)
          ti4 = xr(ivk,p,4) - xr(ivk,p,8)
          vr1 = tr1 + tr3
          vi1 = ti1 + ti3
          vr2 = tr2 + tr4
          vi2 = ti2 + ti4
          vr3 = tr1 - tr3
          vi3 = ti1 - ti3
          vr4 = tr2 - tr4
          vi4 = ti2 - ti4
          yr(ivk,1,p) = ur1 + vr1
          yi(ivk,1,p) = ui1 + vi1
          tr = ur1 - vr1
          ti = ui1 - vi1
          yr(ivk,5,p) = wr4*tr + wi4*ti
          yi(ivk,5,p) = wr4*ti - wi4*tr

          sr = sp2*(vr2 - vi2)
          si = sp2*(vi2 + vr2)
          tr = ur2 + sr
          ti = ui2 + si
          yr(ivk,2,p) = wr1*tr + wi1*ti
          yi(ivk,2,p) = wr1*ti - wi1*tr
          tr = ur2 - sr
          ti = ui2 - si
          yr(ivk,6,p) = wr5*tr + wi5*ti
          yi(ivk,6,p) = wr5*ti - wi5*tr

          tr = ur3 - vi3
          ti = ui3 + vr3
          yr(ivk,3,p) = wr2*tr + wi2*ti
          yi(ivk,3,p) = wr2*ti - wi2*tr
          tr = ur3 + vi3
          ti = ui3 - vr3
          yr(ivk,7,p) = wr6*tr + wi6*ti
          yi(ivk,7,p) = wr6*ti - wi6*tr

          sr = -sp2*(vr4 + vi4)
          si = sp2*(vr4 - vi4)
          tr = ur4 + sr
          ti = ui4 + si
          yr(ivk,4,p) = wr3*tr + wi3*ti
          yi(ivk,4,p) = wr3*ti - wi3*tr
          tr = ur4 - sr
          ti = ui4 - si
          yr(ivk,8,p) = wr7*tr + wi7*ti
          yi(ivk,8,p) = wr7*ti - wi7*tr
       end do
    end do
    return
  end subroutine bftrd8
  !
  subroutine fftv(nv,n,nf,fac,w,lxyz,dr,di,er,ei,ix)
    implicit none
    integer,intent(in):: nv,n,lxyz
    integer,intent(in):: nf,fac(nfmax)
    real(kind=8),intent(in):: w(2,n)
    real(kind=8),dimension(lxyz):: dr,di,er,ei
    integer:: ix
    !
    integer:: iw,nvk,np,uif,if
    !
    iw = 1
    nvk = nv
    np = n/fac(1)
    uif = mod(nf,2)
    do if = 1,nf-uif,2
       select case (fac(if))
       case (2)
          call fftrd2(nvk,np,w(1,iw),dr,di,er,ei)
       case (3)
          call fftrd3(nvk,np,w(1,iw),dr,di,er,ei)
       case (4)
          call fftrd4(nvk,np,w(1,iw),dr,di,er,ei)
       case (5)
          call fftrd5(nvk,np,w(1,iw),dr,di,er,ei)
       case (6)
          call fftrd6(nvk,np,w(1,iw),dr,di,er,ei)
       case (8)
          call fftrd8(nvk,np,w(1,iw),dr,di,er,ei)
       case default
          write(6,*)'unsupported factor'
          stop
       end select
       iw=iw+np*(fac(if)-1)
       nvk=nvk*fac(if)
       np=np/fac(if+1)
       select case (fac(if+1))
       case (2)
          call fftrd2(nvk,np,w(1,iw),er,ei,dr,di)
       case (3)
          call fftrd3(nvk,np,w(1,iw),er,ei,dr,di)
       case (4)
          call fftrd4(nvk,np,w(1,iw),er,ei,dr,di)
       case (5)
          call fftrd5(nvk,np,w(1,iw),er,ei,dr,di)
       case (6)
          call fftrd6(nvk,np,w(1,iw),er,ei,dr,di)
       case (8)
          call fftrd8(nvk,np,w(1,iw),er,ei,dr,di)
       case default
          write(6,*)'unsupported factor'
          stop
       end select
       iw=iw+np*(fac(if+1)-1)
       nvk=nvk*fac(if+1)
       if (if+2 .le. nf) then
          np=np/fac(if+2)
       end if
    end do
    if (uif .eq. 1) then
       select case (fac(nf))
       case (2)
          call fftrd2(nvk,np,w(1,iw),dr,di,er,ei)
       case (3)
          call fftrd3(nvk,np,w(1,iw),dr,di,er,ei)
       case (4)
          call fftrd4(nvk,np,w(1,iw),dr,di,er,ei)
       case (5)
          call fftrd5(nvk,np,w(1,iw),dr,di,er,ei)
       case (6)
          call fftrd6(nvk,np,w(1,iw),dr,di,er,ei)
       case (8)
          call fftrd8(nvk,np,w(1,iw),dr,di,er,ei)
       case default
          write(6,*)'unsupported factor'
          stop
       end select
       ix = - ix
    end if

    return
  end subroutine fftv
  !
  subroutine bftv(nv,n,nf,fac,w,lxyz,dr,di,er,ei,ix)
    implicit none
    integer,intent(in):: nv,n,lxyz
    integer,intent(in):: nf,fac(nfmax)
    real(kind=8),intent(in):: w(2,n)
    real(kind=8),dimension(lxyz):: dr,di,er,ei
    integer:: ix
    !
    integer:: iw,nvk,np,uif,if
    !
    iw = 1
    nvk = nv
    np = n/fac(1)
    uif = mod(nf,2)
    do if = 1,nf-uif,2
       select case (fac(if))
       case (2)
          call bftrd2(nvk,np,w(1,iw),dr,di,er,ei)
       case (3)
          call bftrd3(nvk,np,w(1,iw),dr,di,er,ei)
       case (4)
          call bftrd4(nvk,np,w(1,iw),dr,di,er,ei)
       case (5)
          call bftrd5(nvk,np,w(1,iw),dr,di,er,ei)
       case (6)
          call bftrd6(nvk,np,w(1,iw),dr,di,er,ei)
       case (8)
          call bftrd8(nvk,np,w(1,iw),dr,di,er,ei)
       case default
          write(6,*)'unsupported factor'
          stop
       end select
       iw=iw+np*(fac(if)-1)
       nvk=nvk*fac(if)
       np=np/fac(if+1)
       select case(fac(if+1))
       case (2)
          call bftrd2(nvk,np,w(1,iw),er,ei,dr,di)
       case (3)
          call bftrd3(nvk,np,w(1,iw),er,ei,dr,di)
       case (4)
          call bftrd4(nvk,np,w(1,iw),er,ei,dr,di)
       case (5)
          call bftrd5(nvk,np,w(1,iw),er,ei,dr,di)
       case (6)
          call bftrd6(nvk,np,w(1,iw),er,ei,dr,di)
       case (8)
          call bftrd8(nvk,np,w(1,iw),er,ei,dr,di)
       case default
          write(6,*)'unsupported factor'
          stop
       end select
       iw=iw+np*(fac(if+1)-1)
       nvk=nvk*fac(if+1)
       if (if+2 .le. nf) then
          np=np/fac(if+2)
       end if
    end do
    if (uif .eq. 1) then
       select case (fac(nf))
       case (2)
          call bftrd2(nvk,np,w(1,iw),dr,di,er,ei)
       case (3)
          call bftrd3(nvk,np,w(1,iw),dr,di,er,ei)
       case (4)
          call bftrd4(nvk,np,w(1,iw),dr,di,er,ei)
       case (5)
          call bftrd5(nvk,np,w(1,iw),dr,di,er,ei)
       case (6)
          call bftrd6(nvk,np,w(1,iw),dr,di,er,ei)
       case (8)
          call bftrd8(nvk,np,w(1,iw),dr,di,er,ei)
       case default
          write(6,*)'unsupported factor'
          stop
       end select
       ix = - ix
    end if

    return
  end subroutine bftv
  !
  subroutine fft_pad_zero(nx,ny,nz,lx,ly,lz,dat)
    implicit none
    integer,intent(in):: nx,ny,nz,lx,ly,lz
    real(kind=8):: dat(lx,ly,lz,2)
    !
    integer:: j,k
    do k=1,2
       do j=nz+1,lz
          dat(:,:,j,k) = 0.0d0
       end do
       do j=nx+1,lx
          dat(j,:,:,k) = 0.0d0
       end do
       do j=ny+1,ly
          dat(:,j,:,k) = 0.0d0
       end do
    end do
    return
  end subroutine fft_pad_zero
  !
  subroutine fft_tp1(lxy,lz,src,dst)
    implicit none
    integer,intent(in):: lxy,lz
    real(kind=8),intent(in):: src(lxy,lz)
    real(kind=8),intent(out):: dst(lz,lxy)
    integer:: i,j,ui,uj
    real(kind=8):: r11,r12,r13,r14
    real(kind=8):: r21,r22,r23,r24
    real(kind=8):: r31,r32,r33,r34
    real(kind=8):: r41,r42,r43,r44

    ui = mod(lz,4)
    uj = mod(lxy,4)
!SOPTION LOOP(MAX .LE. 3)
    do i=1,ui
       do j=1,lxy
          dst(i,j) = src(j,i)
       end do
    end do
!SOPTION LOOP(MAX .LE. 3)
    do j=1,uj
       do i=ui+1,lz
          dst(i,j) = src(j,i)
       end do
    end do
    do i=ui+1,lz,4
       do j=uj+1,lxy,4
          r11 = src(j,i)
          r21 = src(j+1,i)
          r31 = src(j+2,i)
          r41 = src(j+3,i)
          r12 = src(j,i+1)
          r22 = src(j+1,i+1)
          r32 = src(j+2,i+1)
          r42 = src(j+3,i+1)
          r13 = src(j,i+2)
          r23 = src(j+1,i+2)
          r33 = src(j+2,i+2)
          r43 = src(j+3,i+2)
          r14 = src(j,i+3)
          r24 = src(j+1,i+3)
          r34 = src(j+2,i+3)
          r44 = src(j+3,i+3)

          dst(i,j) = r11
          dst(i+1,j) = r12
          dst(i+2,j) = r13
          dst(i+3,j) = r14
          dst(i,j+1) = r21
          dst(i+1,j+1) = r22
          dst(i+2,j+1) = r23
          dst(i+3,j+1) = r24
          dst(i,j+2) = r31
          dst(i+1,j+2) = r32
          dst(i+2,j+2) = r33
          dst(i+3,j+2) = r34
          dst(i,j+3) = r41
          dst(i+1,j+3) = r42
          dst(i+2,j+3) = r43
          dst(i+3,j+3) = r44
       end do
    end do
    return
  end subroutine fft_tp1
  !
  subroutine fft_tp(lxy,lz,src,dst)
    implicit none
    integer,intent(in):: lxy,lz
    real(kind=8),intent(in):: src(lxy,lz,2)
    real(kind=8),intent(out):: dst(lz,lxy,2)
    !
    call fft_tp1(lxy,lz,src(1,1,1),dst(1,1,1))
    call fft_tp1(lxy,lz,src(1,1,2),dst(1,1,2))
    return
  end subroutine fft_tp
  !
  subroutine fft_tps1(lxy,lz,src,dst,scl)
    implicit none
    integer,intent(in):: lxy,lz
    real(kind=8),intent(in):: src(lxy,lz)
    real(kind=8),intent(out):: dst(lz,lxy)
    real(kind=8),intent(in):: scl
    integer:: i,j,ui,uj
    real(kind=8):: r11,r12,r13,r14
    real(kind=8):: r21,r22,r23,r24
    real(kind=8):: r31,r32,r33,r34
    real(kind=8):: r41,r42,r43,r44

    ui = mod(lz,4)
    uj = mod(lxy,4)
!SOPTION LOOP(MAX .LE. 3)
    do i=1,ui
       do j=1,lxy
          dst(i,j) = scl*src(j,i)
       end do
    end do
!SOPTION LOOP(MAX .LE. 3)
    do j=1,uj
       do i=ui+1,lz
          dst(i,j) = scl*src(j,i)
       end do
    end do
    do i=ui+1,lz,4
       do j=uj+1,lxy,4
          r11 = src(j,i)
          r21 = src(j+1,i)
          r31 = src(j+2,i)
          r41 = src(j+3,i)
          r12 = src(j,i+1)
          r22 = src(j+1,i+1)
          r32 = src(j+2,i+1)
          r42 = src(j+3,i+1)
          r13 = src(j,i+2)
          r23 = src(j+1,i+2)
          r33 = src(j+2,i+2)
          r43 = src(j+3,i+2)
          r14 = src(j,i+3)
          r24 = src(j+1,i+3)
          r34 = src(j+2,i+3)
          r44 = src(j+3,i+3)

          dst(i,j) = scl*r11
          dst(i+1,j) = scl*r12
          dst(i+2,j) = scl*r13
          dst(i+3,j) = scl*r14
          dst(i,j+1) = scl*r21
          dst(i+1,j+1) = scl*r22
          dst(i+2,j+1) = scl*r23
          dst(i+3,j+1) = scl*r24
          dst(i,j+2) = scl*r31
          dst(i+1,j+2) = scl*r32
          dst(i+2,j+2) = scl*r33
          dst(i+3,j+2) = scl*r34
          dst(i,j+3) = scl*r41
          dst(i+1,j+3) = scl*r42
          dst(i+2,j+3) = scl*r43
          dst(i+3,j+3) = scl*r44
       end do
    end do
    return
  end subroutine fft_tps1
  !
  subroutine fft_tps(lxy,lz,src,dst,scl)
    implicit none
    integer,intent(in):: lxy,lz
    real(kind=8),intent(in):: src(lxy,lz,2)
    real(kind=8),intent(out):: dst(lz,lxy,2)
    real(kind=8),intent(in):: scl
    !
    call fft_tps1(lxy,lz,src(1,1,1),dst(1,1,1),scl)
    call fft_tps1(lxy,lz,src(1,1,2),dst(1,1,2),scl)
    return
  end subroutine fft_tps
  !
  subroutine fft3_init(mx,my,mz,lmx,lmy,lmz,fs)
    implicit none
    integer,intent(in):: mx,my,mz,lmx,lmy,lmz
    type (fft3_struct):: fs
    !
    fs%nx = mx
    fs%ny = my
    fs%nz = mz
    fs%nxyz = mx*my*mz
    fs%lx = lmx
    fs%ly = lmy
    fs%lz = lmz
    fs%lxyz = lmx*lmy*lmz
    call fact235(mx,fs%nfx,fs%facx)
    call fact235(my,fs%nfy,fs%facy)
    call fact235(mz,fs%nfz,fs%facz)
    allocate(fs%wx(2,mx),fs%wy(2,my),fs%wz(2,mz))
    call genw(mx,fs%nfx,fs%facx,fs%wx)
    call genw(my,fs%nfy,fs%facy,fs%wy)
    call genw(mz,fs%nfz,fs%facz,fs%wz)
    return
  end subroutine fft3_init

  subroutine fft3_done(fs)
    implicit none
    type (fft3_struct):: fs
    !
    deallocate(fs%wx,fs%wy,fs%wz)
    return
  end subroutine fft3_done

  subroutine fft3_fw(fs,dat,tmp)
    implicit none
    type (fft3_struct):: fs
    real(kind=8):: dat(fs%lxyz,2)
    real(kind=8):: tmp(fs%lxyz,2)
    integer:: ix
    real(kind=8):: scl
    !
    call fft_pad_zero(fs%nx,fs%ny,fs%nz,fs%lx,fs%ly,fs%lz,tmp)
    ix = 1
    !
    call fftv(fs%lx*fs%ly, fs%nz, fs%nfz, fs%facz, fs%wz, &
         fs%lxyz, dat(1,1), dat(1,2), tmp(1,1), tmp(1,2), ix)
    !
    if (ix .eq. 1) then
       call fft_tp(fs%lx*fs%ly, fs%lz, dat, tmp)
       ix = - ix
       !
       call fftv(fs%lz*fs%lx, fs%ny, fs%nfy, fs%facy, fs%wy, &
            fs%lxyz, tmp(1,1), tmp(1,2), dat(1,1), dat(1,2), ix)
    else if (ix .eq. -1) then
       call fft_tp(fs%lx*fs%ly, fs%lz, tmp, dat)
       ix = - ix
       !
       call fftv(fs%lz*fs%lx, fs%ny, fs%nfy, fs%facy, fs%wy, &
            fs%lxyz, dat(1,1), dat(1,2), tmp(1,1), tmp(1,2), ix)
    else
       write(6,*) 'fft3_fw: y, abs(ix) != 1'
       stop
    end if
    !
    if (ix .eq. 1) then
       call fft_tp(fs%lz*fs%lx, fs%ly, dat, tmp)
       ix = - ix
       !
       call fftv(fs%ly*fs%lz, fs%nx, fs%nfx, fs%facx, fs%wx, &
            fs%lxyz, tmp(1,1), tmp(1,2), dat(1,1), dat(1,2), ix)
    else if (ix .eq. -1) then
       call fft_tp(fs%lz*fs%lx, fs%ly, tmp, dat)
       ix = - ix
       !
       call fftv(fs%ly*fs%lz, fs%nx, fs%nfx, fs%facx, fs%wx, &
            fs%lxyz, dat(1,1), dat(1,2), tmp(1,1), tmp(1,2), ix)
    else
       write(6,*) 'fft3_fw: x, abs(ix) != 1'
       stop
    end if
    !
    scl = 1.0d0/fs%nxyz
    if (ix .eq. 1) then
       call fft_tp(fs%ly*fs%lz, fs%lx, dat, tmp)
       !
       dat = scl*tmp
    else if (ix .eq. -1) then
       call fft_tp(fs%ly*fs%lz, fs%lx, tmp, dat)
       !
       dat = scl*dat
    else
       write(6,*) 'fft3_fw: f, abs(ix) != 1'
       stop
    end if
    return
  end subroutine fft3_fw

  subroutine fft3_bw(fs,dat,tmp)
    implicit none
    type (fft3_struct):: fs
    real(kind=8):: dat(fs%lxyz,2)
    real(kind=8):: tmp(fs%lxyz,2)
    integer:: ix
    !
    call fft_pad_zero(fs%nx,fs%ny,fs%nz,fs%lx,fs%ly,fs%lz,tmp)
    ix = 1
    !
    call bftv(fs%lx*fs%ly, fs%nz, fs%nfz, fs%facz, fs%wz, &
         fs%lxyz, dat(1,1), dat(1,2), tmp(1,1), tmp(1,2), ix)
    !
    if (ix .eq. 1) then
       call fft_tp(fs%lx*fs%ly, fs%lz, dat, tmp)
       ix = - ix
       !
       call bftv(fs%lz*fs%lx, fs%ny, fs%nfy, fs%facy, fs%wy, &
            fs%lxyz, tmp(1,1), tmp(1,2), dat(1,1), dat(1,2), ix)
    else if (ix .eq. -1) then
       call fft_tp(fs%lx*fs%ly, fs%lz, tmp, dat)
       ix = - ix
       !
       call bftv(fs%lz*fs%lx, fs%ny, fs%nfy, fs%facy, fs%wy, &
            fs%lxyz, dat(1,1), dat(1,2), tmp(1,1), tmp(1,2), ix)
    else
       write(6,*) 'fft3_bw: y, abs(ix) != 1'
       stop
    end if
    !
    if (ix .eq. 1) then
       call fft_tp(fs%lz*fs%lx, fs%ly, dat, tmp)
       ix = - ix
       !
       call bftv(fs%ly*fs%lz, fs%nx, fs%nfx, fs%facx, fs%wx, &
            fs%lxyz, tmp(1,1), tmp(1,2), dat(1,1), dat(1,2), ix)
    else if (ix .eq. -1) then
       call fft_tp(fs%lz*fs%lx, fs%ly, tmp, dat)
       ix = - ix
       !
       call bftv(fs%ly*fs%lz, fs%nx, fs%nfx, fs%facx, fs%wx, &
            fs%lxyz, dat(1,1), dat(1,2), tmp(1,1), tmp(1,2), ix)
    else
       write(6,*) 'fft3_bw: x, abs(ix) != 1'
       stop
    end if
    !
    if (ix .eq. 1) then
       call fft_tp(fs%ly*fs%lz, fs%lx, dat, tmp)
       dat = tmp
    else if (ix .eq. -1) then
       call fft_tp(fs%ly*fs%lz, fs%lx, tmp, dat)
    else
       write(6,*) 'fft3_bw: f, abs(ix) != 1'
       stop
    end if
    return
  end subroutine fft3_bw
  !
  subroutine fft3_fwsx(fs,dat,tmp,ix)
    implicit none
    type (fft3_struct):: fs
    real(kind=8):: dat(fs%lxyz,2)
    real(kind=8):: tmp(fs%lxyz,2)
    integer:: ix
    real(kind=8):: scl
    !
    call fft_pad_zero(fs%nx,fs%ny,fs%nz,fs%lx,fs%ly,fs%lz,tmp)
    !
    if (ix .eq. 1) then
       call fftv(fs%lx*fs%ly, fs%nz, fs%nfz, fs%facz, fs%wz, &
            fs%lxyz, dat(1,1), dat(1,2), tmp(1,1), tmp(1,2), ix)
    else if (ix .eq. -1) then
       call fftv(fs%lx*fs%ly, fs%nz, fs%nfz, fs%facz, fs%wz, &
            fs%lxyz, tmp(1,1), tmp(1,2), dat(1,1), dat(1,2), ix)
    else
       write(6,*) 'fft3_fwsx: z, abs(ix) != 1'
       stop
    end if
    !
    if (ix .eq. 1) then
       call fft_tp(fs%lx*fs%ly, fs%lz, dat, tmp)
       ix = - ix
       !
       call fftv(fs%lz*fs%lx, fs%ny, fs%nfy, fs%facy, fs%wy, &
            fs%lxyz, tmp(1,1), tmp(1,2), dat(1,1), dat(1,2), ix)
    else if (ix .eq. -1) then
       call fft_tp(fs%lx*fs%ly, fs%lz, tmp, dat)
       ix = - ix
       !
       call fftv(fs%lz*fs%lx, fs%ny, fs%nfy, fs%facy, fs%wy, &
            fs%lxyz, dat(1,1), dat(1,2), tmp(1,1), tmp(1,2), ix)
    else
       write(6,*) 'fft3_fwsx: y, abs(ix) != 1'
       stop
    end if
    !
    scl = 1.0d0/fs%nxyz
    !
    if (ix .eq. 1) then
       call fft_tps(fs%lz*fs%lx, fs%ly, dat, tmp, scl)
       ix = - ix
       !
       call fftv(fs%ly*fs%lz, fs%nx, fs%nfx, fs%facx, fs%wx, &
            fs%lxyz, tmp(1,1), tmp(1,2), dat(1,1), dat(1,2), ix)
    else if (ix .eq. -1) then
       call fft_tps(fs%lz*fs%lx, fs%ly, tmp, dat, scl)
       ix = - ix
       !
       call fftv(fs%ly*fs%lz, fs%nx, fs%nfx, fs%facx, fs%wx, &
            fs%lxyz, dat(1,1), dat(1,2), tmp(1,1), tmp(1,2), ix)
    else
       write(6,*) 'fft3_fwsx: x, abs(ix) != 1'
       stop
    end if
    !
    return
  end subroutine fft3_fwsx

  subroutine fft3_bwsx(fs,dat,tmp,ix)
    implicit none
    type (fft3_struct):: fs
    real(kind=8):: dat(fs%lxyz,2)
    real(kind=8):: tmp(fs%lxyz,2)
    integer:: ix
    !
    call fft_pad_zero(fs%ny,fs%nz,fs%nx,fs%ly,fs%lz,fs%lx,tmp)
    !
    if (ix .eq. 1) then
       call bftv(fs%lz*fs%lx, fs%ny, fs%nfy, fs%facy, fs%wy, &
            fs%lxyz, dat(1,1), dat(1,2), tmp(1,1), tmp(1,2), ix)
    else if (ix .eq. -1) then
       call bftv(fs%lz*fs%lx, fs%ny, fs%nfy, fs%facy, fs%wy, &
            fs%lxyz, tmp(1,1), tmp(1,2), dat(1,1), dat(1,2), ix)
    else
       write(6,*) 'fft3_bwsx: x, abs(ix) != 1'
       stop
    end if
    !
    if (ix .eq. 1) then
       call fft_tp(fs%lz*fs%lx, fs%ly, dat, tmp)
       ix = - ix
       !
       call bftv(fs%ly*fs%lz, fs%nx, fs%nfx, fs%facx, fs%wx, &
            fs%lxyz, tmp(1,1), tmp(1,2), dat(1,1), dat(1,2), ix)
    else if (ix .eq. -1) then
       call fft_tp(fs%lz*fs%lx, fs%ly, tmp, dat)
       ix = - ix
       !
       call bftv(fs%ly*fs%lz, fs%nx, fs%nfx, fs%facx, fs%wx, &
            fs%lxyz, dat(1,1), dat(1,2), tmp(1,1), tmp(1,2), ix)
    else
       write(6,*) 'fft3_bwsx: y, abs(ix) != 1'
       stop
    end if
    !
    if (ix .eq. 1) then
       call fft_tp(fs%ly*fs%lz, fs%lx, dat, tmp)
       ix = - ix
       !
       call bftv(fs%lx*fs%ly, fs%nz, fs%nfz, fs%facz, fs%wz, &
         fs%lxyz, tmp(1,1), tmp(1,2), dat(1,1), dat(1,2), ix)
    else if (ix .eq. -1) then
       call fft_tp(fs%ly*fs%lz, fs%lx, tmp, dat)
       ix = - ix
       !
       call bftv(fs%lx*fs%ly, fs%nz, fs%nfz, fs%facz, fs%wz, &
         fs%lxyz, dat(1,1), dat(1,2), tmp(1,1), tmp(1,2), ix)
    else
       write(6,*) 'fft3_bwsx: z, abs(ix) != 1'
       stop
    end if
    !
    return
  end subroutine fft3_bwsx

end module fft_3d
