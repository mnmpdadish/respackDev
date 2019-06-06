! The MIT License (MIT)
! Copyright (c) 2016, Yoshiro Nohara, Takeo Fujiwara, and University of Tokyo.
!
! Permission is hereby granted, free of charge, to any person
! obtaining a copy of this software and associated documentation
! files (the "Software"), to deal in the Software without restriction,
! including without limitation the rights to use, copy, modify, merge,
! publish, distribute, sublicense, and/or sell copies of the Software,
! and to permit persons to whom the Software is furnished to do so,
! subject to the following conditions:
!
! The above copyright notice and this permission notice shall be
! included in all copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
! EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
! OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
! NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
! HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
! WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
! IN THE SOFTWARE.
!
!     This program uses the generalized tetrahedron method and implements
!     an efficient tool in the GW-LMTO package developed by Fujiwara group
!     in University of Tokyo. The algorithm of the generalized tetrahedron
!     method is written in the following published papers:
!
!     [Generalized tetrahedron method] Takeo Fujiwara, Susumu Yamamoto, and
!     Yasushi Ishii, J. Phys. Soc. Jpn. 72, No.4, 777-780 (2003).
!     (cond.mat/0211159)
!
!     [GW-LMTO package] Yoshiro Nohara, Susumu Yamamoto, and Takeo Fujiwara,
!     Phys. Rev. B 79, 195110 (2009)
!
!     One can use this part of the program package in one's developed code
!     on the condition of citing these two papers and mentioning the name of
!     `the generalized tetrahedron method'  in the published paper.
      module m_tetrainteg
      implicit none
      private

      public :: ttritg_mkidx
      public :: ttritg_sum
      public :: ttritg_simple
      public :: ttritg_search
      public :: ttritg_sub1normal
      public :: ttritg_sub1simple
      public :: ttritg_sub1g30g20
      public :: ttritg_sub1g30g21
      public :: ttritg_sub1g30
      public :: ttritg_sub2normal
      public :: ttritg_sub2simple
      public :: ttritg_sub2g30g20
      public :: ttritg_sub2g30g10
      public :: ttritg_sub2g30g21
      public :: ttritg_sub2g30
      public :: ttritg_sub2g10


      integer, parameter, public :: ttritg_vtx=4
      integer, parameter, public :: ttritg_div=6
      integer, parameter :: nvtx=ttritg_vtx
      integer, parameter :: ndiv=ttritg_div
      integer, parameter, public :: ttritg_tbl1(nvtx,ndiv)
     $                       =reshape((/1,2,3,6
     $                                 ,1,3,5,6
     $                                 ,3,5,6,7
     $                                 ,2,3,4,6
     $                                 ,3,4,6,8
     $                                 ,3,6,7,8/)
     $                               ,(/nvtx,ndiv/))

      contains

      subroutine ttritg_search(dmna,dmnr,n1,n2,n3,n4,imt1,fk,gk,xo
     $                        ,dtgt,drg1,drg2,dout)
      real(8), intent(in) :: dmna,dmnr
      integer, intent(in) :: n1,n2,n3,n4
      integer, intent(in) :: imt1(nvtx,n1*n2*n3*ndiv)
      complex(8), intent(in) :: fk(n1*n2*n3,n4)
      complex(8), intent(in) :: gk(n1*n2*n3,n4)
      complex(8), intent(inout) :: xo(n1*n2*n3)
      real(8), intent(in) :: dtgt,drg1,drg2
      real(8), intent(out):: dout
      integer :: i1,i2,i3
      real(8) :: d1,d2,d3,d4,d5,d6
      logical :: b1,b2,b3
      complex(8) :: ca1(nvtx),ca2(1)
      complex(8), allocatable :: ca3(:)
      d3=drg1
      d4=drg2
      xo=0.d0
      ca2=d3
      do i1=1,n4
        call ttritg_sum(dmna,dmnr,n1,n2,n3,imt1
     $                 ,fk(:,i1),gk(:,i1),1,ca2,ca1,xo)
      end do
      d1=abs(dimag(sum(xo)))
      xo=0.d0
      ca2=d4
      do i1=1,n4
        call ttritg_sum(dmna,dmnr,n1,n2,n3,imt1
     $                 ,fk(:,i1),gk(:,i1),1,ca2,ca1,xo)
      end do
      d2=abs(dimag(sum(xo)))
#ifdef MPI
#else
      write(*,"(a,2g25.16)")"#ttritg_search: d1,d2=",d1,d2
#endif
      if((d1-dtgt)*(dtgt-d2).lt.0.d0) 
     $  call errexec("bad drg1,drg2,dtgt")
      b1=.false.
      do while(.not.b1)
        d6=(d3+d4)*.5d0
        b1=abs(d4-d6).lt.1.d-15.or.abs((d4-d6)/(d4+d6))*2.d0.lt.1.d-15
        xo=0.d0
        ca2=d6
!$OMP PARALLEL PRIVATE(ca1,ca3)
        allocate(ca3(n1*n2*n3)); ca3=(0.d0,0.d0)
!$OMP DO
        do i1=1,n4
          call ttritg_sum(dmna,dmnr,n1,n2,n3,imt1
     $                   ,fk(:,i1),gk(:,i1),1,ca2,ca1,ca3)
        end do
!$OMP END DO
!$OMP CRITICAL
        xo=xo+ca3
!$OMP END CRITICAL
        deallocate(ca3)
!$OMP END PARALLEL
        d5=abs(dimag(sum(xo)))
        if((d1-dtgt)*(dtgt-d5).ge.0.d0)then
          d4=d6
        else
          d3=d6
        end if
#ifdef MPI
#else
        write(*,"(a,2g25.16)")"#ttritg_search: d6,d5=",d6,d5
#endif
      end do
      dout=d6
      end subroutine

      subroutine ttritg_simple(dmna,dmnr,n1,n2,n3,imt1,fk,gk,xo)
      real(8), intent(in) :: dmna,dmnr
      integer, intent(in) :: n1
      integer, intent(in) :: n2
      integer, intent(in) :: n3
      integer, intent(in) :: imt1(nvtx,n1*n2*n3*ndiv)
      complex(8), intent(in) :: fk(n1*n2*n3)
      complex(8), intent(in) :: gk(n1*n2*n3)
      complex(8), intent(inout) :: xo(n1*n2*n3)
      complex(8) :: ca1(nvtx),ca2(1)
      ca2=(0.d0,0.d0)
      call ttritg_sum(dmna,dmnr,n1,n2,n3,imt1,fk,gk
     $                ,1,ca2,ca1,xo)
      end subroutine

      subroutine ttritg_sum(dmna,dmnr,n1,n2,n3,imt1,fk,gk,ne,em,ca1,xo)
      real(8), intent(in) :: dmna,dmnr
      integer, intent(in) :: n1
      integer, intent(in) :: n2
      integer, intent(in) :: n3
      integer, intent(in) :: imt1(nvtx,n1*n2*n3*ndiv)
      complex(8), intent(in) :: fk(n1*n2*n3)
      complex(8), intent(in) :: gk(n1*n2*n3)
      integer, intent(in) :: ne
      complex(8), intent(in) :: em(ne)
      complex(8), intent(  out) :: ca1(ne*nvtx)
      complex(8), intent(inout) :: xo(ne,n1*n2*n3)
      integer :: i1,i2,i3,i4,i5,i6,i7
      if(n1.eq.1.and.n2.eq.1.and.n3.eq.1)then
        xo(:,1)=xo(:,1)+fk(1)*log((gk(1)+em)/(-1.d0))
        return
      end if
      do i1=1,n1*n2*n3*ndiv
        call calc(dmna,dmnr,imt1(:,i1),n1*n2*n3
     $           ,fk
     $           ,gk,ne,em,ca1
     $           ,xo)
      end do
      end subroutine

      subroutine calc(dmna,dmnr,ia1,nnn,fk,gk,ne,em,ca1,xo)
      real(8), intent(in) :: dmna,dmnr
      integer, intent(in) :: ia1(nvtx)
      integer, intent(in) :: nnn
      complex(8), intent(in) :: fk(nnn)
      complex(8), intent(in) :: gk(nnn)
      integer, intent(in) :: ne
      complex(8), intent(in) :: em(ne)
      complex(8), intent(  out) :: ca1(ne,nvtx)
      complex(8), intent(inout) :: xo(ne,nnn)
      integer :: i1,i2,i3,i4,i5,i6
      integer :: ia2(4),ia3(4),ia4(4)
      complex(8) :: geq1,geq2
      call chk1(dmna,dmnr,nnn,gk,ia1,ia2,ia3,ia4,i3)
      select case(sum(ia2))
      case(0)
        call ttritg_sub2normal(
     $             gk(ia4(2)),gk(ia4(1)),gk(ia4(3)),gk(ia4(4))
     $            ,ne,em,ca1(:,1))
        call ttritg_sub2normal(
     $             gk(ia4(1)),gk(ia4(2)),gk(ia4(3)),gk(ia4(4))
     $            ,ne,em,ca1(:,2))
        call ttritg_sub2normal(
     $             gk(ia4(1)),gk(ia4(3)),gk(ia4(2)),gk(ia4(4))
     $            ,ne,em,ca1(:,3))
        call ttritg_sub2normal(
     $             gk(ia4(1)),gk(ia4(4)),gk(ia4(2)),gk(ia4(3))
     $            ,ne,em,ca1(:,4))
          xo(:,ia4(1))=xo(:,ia4(1))+fk(ia4(1))*ca1(:,1)
          xo(:,ia4(2))=xo(:,ia4(2))+fk(ia4(2))*ca1(:,2)
          xo(:,ia4(3))=xo(:,ia4(3))+fk(ia4(3))*ca1(:,3)
          xo(:,ia4(4))=xo(:,ia4(4))+fk(ia4(4))*ca1(:,4)
      case(2)
          geq1=sum(gk(ia3(1:2)))/2.d0
          call ttritg_sub2g10   (
     $             geq1      ,geq1      ,gk(ia4(1)),gk(ia4(2))
     $            ,ne,em,ca1(:,1))
          call ttritg_sub2g30   (
     $             geq1      ,gk(ia4(1)),gk(ia4(2)),geq1      
     $            ,ne,em,ca1(:,2))
          call ttritg_sub2g30   (
     $             geq1      ,gk(ia4(2)),gk(ia4(1)),geq1      
     $            ,ne,em,ca1(:,3))
          xo(:,ia3(1))=xo(:,ia3(1))+fk(ia3(1))*ca1(:,1)
          xo(:,ia3(2))=xo(:,ia3(2))+fk(ia3(2))*ca1(:,1)
          xo(:,ia4(1))=xo(:,ia4(1))+fk(ia4(1))*ca1(:,2)
          xo(:,ia4(2))=xo(:,ia4(2))+fk(ia4(2))*ca1(:,3)
      case(3)
          geq1=sum(gk(ia3(1:3)))/3.d0
          call ttritg_sub2g30g10(
     $             geq1      ,geq1      ,gk(ia4(1)),geq1      
     $            ,ne,em,ca1(:,1))
          call ttritg_sub2g30g20(
     $             geq1      ,gk(ia4(1)),geq1      ,geq1      
     $            ,ne,em,ca1(:,2))
          xo(:,ia3(1))=xo(:,ia3(1))+fk(ia3(1))*ca1(:,1)
          xo(:,ia3(2))=xo(:,ia3(2))+fk(ia3(2))*ca1(:,1)
          xo(:,ia3(3))=xo(:,ia3(3))+fk(ia3(3))*ca1(:,1)
          xo(:,ia4(1))=xo(:,ia4(1))+fk(ia4(1))*ca1(:,2)
      case(4)
        select case(i3)
        case(2)
          geq1=sum(gk(ia3(1:2)))/2.d0; geq2=sum(gk(ia3(3:4)))/2.d0
          call ttritg_sub2g30g21(
     $             geq2      ,geq1      ,geq1      ,geq2      
     $            ,ne,em,ca1(:,1))
          call ttritg_sub2g30g21(
     $             geq1      ,geq2      ,geq2      ,geq1      
     $            ,ne,em,ca1(:,2))
          xo(:,ia3(1))=xo(:,ia3(1))+fk(ia3(1))*ca1(:,1)
          xo(:,ia3(2))=xo(:,ia3(2))+fk(ia3(2))*ca1(:,1)
          xo(:,ia3(3))=xo(:,ia3(3))+fk(ia3(3))*ca1(:,2)
          xo(:,ia3(4))=xo(:,ia3(4))+fk(ia3(4))*ca1(:,2)
        case(3:6)
          geq1=sum(gk(ia3(1:4)))/4.d0
          call ttritg_sub2simple(
     $             geq1      ,geq1      ,geq1      ,geq1      
     $            ,ne,em,ca1(:,1))
          xo(:,ia3(1))=xo(:,ia3(1))+ca1(:,1)*fk(ia3(1))
          xo(:,ia3(2))=xo(:,ia3(2))+ca1(:,1)*fk(ia3(2))
          xo(:,ia3(3))=xo(:,ia3(3))+ca1(:,1)*fk(ia3(3))
          xo(:,ia3(4))=xo(:,ia3(4))+ca1(:,1)*fk(ia3(4))
        case default
          call errexec("bad i3")
        end select
      case default
        call errexec("bad sum ia2")
      end select
      end subroutine

      subroutine ttritg_sub1simple(h0,h1,h2,h3,ne,em,x0)
      complex(8), intent(in) :: h0,h1,h2,h3
      integer, intent(in) :: ne
      complex(8), intent(in) :: em(ne)
      complex(8), intent(out) :: x0(ne)
      x0=log((h0+em)/(-1.d0))/6.d0
      end subroutine

      subroutine ttritg_sub1g30g21(h0,h1,h2,h3,ne,em,x0)
      complex(8), intent(in) :: h0,h1,h2,h3
      integer, intent(in) :: ne
      complex(8), intent(in) :: em(ne)
      complex(8), intent(out) :: x0(ne)
      real(8) :: g0,g1,g2,g3,g01,g02,g03,g12,g13,g23
      real(8) :: gmn1,gmx1,gmn2,gmx2,gmn3,gmx3
      g0=dble(h0); g1=dble(h1); g2=dble(h2); g3=dble(h3)
      g01=g0-g1; g02=g0-g2; g03=g0-g3; g12=g1-g2; g13=g1-g3; g23=g2-g3
      gmn1=min(-g0,-g1); gmx1=max(-g0,-g1)
      x0=
     $     +g0**2*(g0-3.d0*g1)*log((h0+em)/(h1+em))
     $     -min(max(dble(em),gmn1),gmx1)*(6.d0*g0*g1
     $     +min(max(dble(em),gmn1),gmx1)*(3.d0*(g0+g1)
     $     +min(max(dble(em),gmn1),gmx1)*2.d0))
     $        *log((-h0-em)/(-h1-em))
      x0=x0/(6.d0*g01**3)
      x0=x0-log((-1.d0)/(h1+em))/6.d0
      end subroutine

      subroutine ttritg_sub1g30g20(h0,h1,h2,h3,ne,em,x0)
      complex(8), intent(in) :: h0,h1,h2,h3
      integer, intent(in) :: ne
      complex(8), intent(in) :: em(ne)
      complex(8), intent(out) :: x0(ne)
      real(8) :: g0,g1,g2,g3,g01,g02,g03,g12,g13,g23
      real(8) :: gmn1,gmx1,gmn2,gmx2,gmn3,gmx3
      g0=dble(h0); g1=dble(h1); g2=dble(h2); g3=dble(h3)
      g01=g0-g1; g02=g0-g2; g03=g0-g3; g12=g1-g2; g13=g1-g3; g23=g2-g3
      gmn1=min(-g0,-g1); gmx1=max(-g0,-g1)
      x0=-(g1+min(max(dble(em),gmn1),gmx1))**3*log((-h1-em)/(-h0-em))
      x0=x0/(6.d0*g01**3)+log(-h0-em)/6.d0
      end subroutine

      subroutine ttritg_sub1g30   (h0,h1,h2,h3,ne,em,x0)
      complex(8), intent(in) :: h0,h1,h2,h3
      integer, intent(in) :: ne
      complex(8), intent(in) :: em(ne)
      complex(8), intent(out) :: x0(ne)
      real(8) :: g0,g1,g2,g3,g01,g02,g03,g12,g13,g23
      real(8) :: gmn1,gmx1,gmn2,gmx2,gmn3,gmx3
      g0=dble(h0); g1=dble(h1); g2=dble(h2); g3=dble(h3)
      g01=g0-g1; g02=g0-g2; g03=g0-g3; g12=g1-g2; g13=g1-g3; g23=g2-g3
      gmn1=min(-g0,-g1); gmx1=max(-g0,-g1)
      gmn2=min(-g1,-g2); gmx2=max(-g1,-g2)
      x0=
     $ -g12*(g0+min(max(dble(em),gmn1),gmx1))**2
     $        *(g0**2-2.d0*g0*(g1+g2)
     $         +3.d0*g1*g2+(-2.d0*g0+g1+g2)
     $         *min(max(dble(em),gmn1),gmx1))*log((h0+em)/(h1+em))
     $ +g01**2*(
     $         +(g2+min(max(dble(em),gmn2),gmx2))**3
     $          *log((-h2-em)/(-h1-em))
     $         +(g0**2*g1-g0**2*g2-2.d0*g0*g1*g2+2.d0*g0*g2**2+g1*g2**2)
     $          *log((-1.d0)/(h1+em)))
     $ +g01**2*g2**3*log(-h1-em)
      x0=-x0/(6.d0*g01**2*g02**2*g12)
      end subroutine

      subroutine ttritg_sub1normal(h0,h1,h2,h3,ne,em,x0)
      complex(8), intent(in) :: h0,h1,h2,h3
      integer, intent(in) :: ne
      complex(8), intent(in) :: em(ne)
      complex(8), intent(out) :: x0(ne)
      real(8) :: g0,g1,g2,g3,g01,g02,g03,g12,g13,g23
      real(8) :: gmn1,gmx1,gmn2,gmx2,gmn3,gmx3
      g0=dble(h0); g1=dble(h1); g2=dble(h2); g3=dble(h3)
      g01=g0-g1; g02=g0-g2; g03=g0-g3; g12=g1-g2; g13=g1-g3; g23=g2-g3
      gmn1=min(-g0,-g1); gmx1=max(-g0,-g1)
      gmn2=min(-g1,-g3); gmx2=max(-g1,-g3)
      gmn3=min(-g1,-g2); gmx3=max(-g1,-g2)
      x0=
     $  +g12*g13*g23*(g0+min(max(dble(em),gmn1),gmx1))**3
     $   *log((h0+em)/(h1+em))
     $  -g01*g02*g12*(g3+min(max(dble(em),gmn2),gmx2))**3
     $   *log((-h3-em)/(-h1-em))
     $  +g01*g03*g13*(g2+min(max(dble(em),gmn3),gmx3))**3
     $   *log((-h2-em)/(-h1-em))
     $  +g01*g23*(g03*g1*g2**2+g12*g3*(g2*g03+g3*g0))*log(-h1-em)
     $  -g01*g23
     $   *(g0**2*(g1*(g1-g2-g3)+g2*g3)
     $    +g0*g1*(g2*g3-g1*(g2+g3))+g1**2*g2*g3)*log((-1.d0)/(h1+em))
      x0=x0/(6.d0*(g01*g12*g23*g13*g02*g03))
      end subroutine

      subroutine ttritg_sub2simple(h0,h1,h2,h3,ne,em,x0)
      complex(8), intent(in) :: h0,h1,h2,h3
      integer, intent(in) :: ne
      complex(8), intent(in) :: em(ne)
      complex(8), intent(out) :: x0(ne)
      x0=log((h0+em)/(-1.d0))/24.d0
      end subroutine

      subroutine ttritg_sub2g30g20(h0,h1,h2,h3,ne,em,x0)
      complex(8), intent(in) :: h0,h1,h2,h3
      integer, intent(in) :: ne
      complex(8), intent(in) :: em(ne)
      complex(8), intent(out) :: x0(ne)
      real(8) :: g0,g1,g2,g3,g01,g02,g03,g12,g13,g23
      real(8) :: gmn1,gmx1,gmn2,gmx2,gmn3,gmx3
      g0=dble(h0); g1=dble(h1); g2=dble(h2); g3=dble(h3)
      g01=g0-g1; g02=g0-g2; g03=g0-g3; g12=g1-g2; g13=g1-g3; g23=g2-g3
      gmn1=min(-g0,-g1); gmx1=max(-g0,-g1)
      x0=
     $  -(g1+min(max(dble(em),gmn1),gmx1))**3
     $   *(4.d0*g0-g1+3.d0*min(max(dble(em),gmn1),gmx1))
     $   *log((-h1-em)/(-h0-em))
      x0=x0/(24.d0*g01**4)+log(-h0-em)/24.d0
      end subroutine

      subroutine ttritg_sub2g30g21(h0,h1,h2,h3,ne,em,x0)
      complex(8), intent(in) :: h0,h1,h2,h3
      integer, intent(in) :: ne
      complex(8), intent(in) :: em(ne)
      complex(8), intent(out) :: x0(ne)
      call ttritg_sub2g30g20(h1,h0,h1,h1,ne,em,x0)
      end subroutine

      subroutine ttritg_sub2g30g10(h0,h1,h2,h3,ne,em,x0)
      complex(8), intent(in) :: h0,h1,h2,h3
      integer, intent(in) :: ne
      complex(8), intent(in) :: em(ne)
      complex(8), intent(out) :: x0(ne)
      real(8) :: g0,g1,g2,g3,g01,g02,g03,g12,g13,g23
      real(8) :: gmn1,gmx1,gmn2,gmx2,gmn3,gmx3
      g0=dble(h0); g1=dble(h1); g2=dble(h2); g3=dble(h3)
      g01=g0-g1; g02=g0-g2; g03=g0-g3; g12=g1-g2; g13=g1-g3; g23=g2-g3
      gmn1=min(-g0,-g2); gmx1=max(-g0,-g2)
      x0=(g2+min(max(dble(em),gmn1),gmx1))**4*log((h2+em)/(h0+em))
     $  -g02**4*log((-1.d0)/(h0+em))
      x0=x0/(24.d0*g02**4)
      end subroutine

      subroutine ttritg_sub2g10   (h0,h1,h2,h3,ne,em,x0)
      complex(8), intent(in) :: h0,h1,h2,h3
      integer, intent(in) :: ne
      complex(8), intent(in) :: em(ne)
      complex(8), intent(out) :: x0(ne)
      real(8) :: g0,g1,g2,g3,g01,g02,g03,g12,g13,g23
      real(8) :: gmn1,gmx1,gmn2,gmx2,gmn3,gmx3
      g0=dble(h0); g1=dble(h1); g2=dble(h2); g3=dble(h3)
      g01=g0-g1; g02=g0-g2; g03=g0-g3; g12=g1-g2; g13=g1-g3; g23=g2-g3
      gmn1=min(-g0,-g2); gmx1=max(-g0,-g2)
      gmn2=min(-g0,-g3); gmx2=max(-g0,-g3)
      x0=-g03**3*(g2+min(max(dble(em),gmn1),gmx1))**4
     $    *log((h2+em)/(h0+em))
     $   +g02**3*(g3+min(max(dble(em),gmn2),gmx2))**4
     $    *log((h3+em)/(h0+em))
      x0=x0/(24.d0*g02**3*g03**3*g23)-log((-1.d0)/(h0+em))/24.d0
      end subroutine

      subroutine ttritg_sub2g30   (h0,h1,h2,h3,ne,em,x0)
      complex(8), intent(in) :: h0,h1,h2,h3
      integer, intent(in) :: ne
      complex(8), intent(in) :: em(ne)
      complex(8), intent(out) :: x0(ne)
      real(8) :: g0,g1,g2,g3,g01,g02,g03,g12,g13,g23
      real(8) :: gmn1,gmx1,gmn2,gmx2,gmn3,gmx3
      g0=dble(h0); g1=dble(h1); g2=dble(h2); g3=dble(h3)
      g01=g0-g1; g02=g0-g2; g03=g0-g3; g12=g1-g2; g13=g1-g3; g23=g2-g3
      gmn1=min(-g0,-g1); gmx1=max(-g0,-g1)
      gmn2=min(-g0,-g2); gmx2=max(-g0,-g2)
      x0=
     $  +g02**2*(g1+min(max(dble(em),gmn1),gmx1))**3
     $   *(g12*(3.d0*g0-g1)-g2*g01
     $    +(2.d0*g02-3.d0*g01)*min(max(dble(em),gmn1),gmx1))
     $                     *log((-h1-em)/(-h0-em))
     $  +g01**3*(g2+min(max(dble(em),gmn2),gmx2))**4
     $   *log((-h2-em)/(-h0-em))
      x0=x0/(24.d0*g01**3*g02**2*g12**2)+log(-h0-em)/24.d0
      end subroutine

      subroutine ttritg_sub2normal(h0,h1,h2,h3,ne,em,x0)
      complex(8), intent(in) :: h0,h1,h2,h3
      integer, intent(in) :: ne
      complex(8), intent(in) :: em(ne)
      complex(8), intent(out) :: x0(ne)
      real(8) :: g0,g1,g2,g3,g01,g02,g03,g12,g13,g23
      real(8) :: gmn1,gmx1,gmn2,gmx2,gmn3,gmx3
      g0=dble(h0); g1=dble(h1); g2=dble(h2); g3=dble(h3)
      gmn1=min(-g1,-g2); gmx1=max(-g1,-g2)
      gmn2=min(-g1,-g3); gmx2=max(-g1,-g3)
      gmn3=min(-g0,-g1); gmx3=max(-g0,-g1)
      g01=g0-g1; g02=g0-g2; g03=g0-g3; g12=g1-g2; g13=g1-g3; g23=g2-g3
      x0=-g01**2*g03*g13**2*(g2+min(max(dble(em),gmn1),gmx1))**4
     $    *log((h2+em)/(h1+em))
     $   +g01**2*g02*g12**2*(g3+min(max(dble(em),gmn2),gmx2))**4
     $    *log((h3+em)/(h1+em))
     $   -g12**2*g13**2*g23*min(max(dble(em),gmn3),gmx3)
     $    *(2.d0*g0+min(max(dble(em),gmn3),gmx3))
     $    *(2.d0*g0**2+min(max(dble(em),gmn3),gmx3)
     $    *(2.d0*g0+min(max(dble(em),gmn3),gmx3)))
     $    *log((-h1-em)/(-h0-em))
     $ +g12**2*g13**2*g23*(
     $   (g01*(2.d0*g0*g1*g03-g0*g2*g13+g1*g2*g3)
     $   +g02*g0**2*(g23+2.d0*g3)+g12*g0*g1*g03
     $   +g0**2*g2*g23+g0*g1*g2*g3
     $   )*log((-1.d0)/(h1+em))+g0**4*log(-h0-em)
     $                    )
      x0=x0/(24.d0*g01**2*g02*g12**2*g03*g13**2*g23)
      end subroutine





      subroutine chk1(dmna,dmnr,nnn,gk,ia1,ia2,ia3,ia4,i3)
      real(8), intent(in) :: dmna,dmnr
      integer, intent(in) :: nnn
      complex(8), intent(in) :: gk(nnn)
      integer, intent(in) :: ia1(4)
      integer, intent(out):: ia2(4)
      integer, intent(out):: ia3(4)
      integer, intent(out):: ia4(4)
      integer, intent(out):: i3
      integer :: i1,i2,i4,i5
      real(8) :: d1,d2
      i3=0
      i4=0
      ia2=0
      ia3=0
      do i1=1,3
        do i2=i1+1,4
          d1=abs(dble(gk(ia1(mod(i1,4)+1))-gk(ia1(mod(i2,4)+1))))
          d2=abs(dble(gk(ia1(mod(i1,4)+1))+gk(ia1(mod(i2,4)+1))))*.5d0
          if(d1.lt.dmna.or.d1/d2.lt.dmnr)then
            i3=i3+1
            if(ia2(mod(i1,4)+1).eq.0)then
              i4=i4+1
              ia3(i4)=ia1(mod(i1,4)+1)
              ia2(mod(i1,4)+1)=1
            end if
            if(ia2(mod(i2,4)+1).eq.0)then
              i4=i4+1
              ia3(i4)=ia1(mod(i2,4)+1)
              ia2(mod(i2,4)+1)=1
            end if
          else
          end if
        end do
      end do
      i5=0
      ia4=0
      do i1=1,4
        if(ia2(mod(i1,4)+1).eq.0) then
          i5=i5+1
          ia4(i5)=ia1(mod(i1,4)+1)
        end if
      end do
      if(sum(ia2).ne.i4) call errexec("bad ia2 and i4")
      if(i4+i5.ne.4)     call errexec("bad i4 and i5")
      end subroutine

      !subroutine ttritg_mkidx(n1,n2,n3,idx1)
      subroutine ttritg_mkidx(n1,n2,n3,idx1,bvec1,bvec2,bvec3)
      integer, intent(in) :: n1
      integer, intent(in) :: n2
      integer, intent(in) :: n3
      integer, intent(out):: idx1(nvtx,n1*n2*n3*ndiv)
      real(8), intent(in), optional :: bvec1(3),bvec2(3),bvec3(3)
      real(8) :: bvec(3,3)
      integer :: i1,i2,i3,i4,i5,i6,i7
      integer :: imt1(8,8)
      integer :: ia1(8)
      if(n1.eq.1.and.n2.eq.1.and.n3.eq.1)return
!--
!      imt1(:,1)=(/0,1,n1,n1+1,n1*n2,n1*n2+1,n1*n2+n1,n1*n2+n1+1/)
!      imt1(:,2)=imt1(:,1)-(/0,n1,0,n1,0,n1,0,n1/)
!      imt1(:,3)=imt1(:,1)
!     $                -(/0,0,n1*n2,n1*n2,0,0,n1*n2,n1*n2/)
!      imt1(:,4)=imt1(:,3)-(/0,n1,0,n1,0,n1,0,n1/)
!      imt1(:,5)=imt1(:,1)
!     $                -(/0,0,0,0,n1*n2*n3,n1*n2*n3,n1*n2*n3,n1*n2*n3/)
!      imt1(:,6)=imt1(:,2)
!     $                -(/0,0,0,0,n1*n2*n3,n1*n2*n3,n1*n2*n3,n1*n2*n3/)
!      imt1(:,7)=imt1(:,3)
!     $                -(/0,0,0,0,n1*n2*n3,n1*n2*n3,n1*n2*n3,n1*n2*n3/)
!      imt1(:,8)=imt1(:,4)
!     $                -(/0,0,0,0,n1*n2*n3,n1*n2*n3,n1*n2*n3,n1*n2*n3/)
!      i6=0
!      i4=0
!      do i3=1,n3
!      do i2=1,n2
!      do i1=1,n1
!        i4=i4+1
!        i5=1
!        if (i1.eq.n1) i5=i5+1
!        if (i2.eq.n2) i5=i5+2
!        if (i3.eq.n3) i5=i5+4
!        ia1=i4+imt1(:,i5)
!        do i7=1,ndiv
!          i6=i6+1
!          idx1(:,i6)=ia1(ttritg_tbl1(:,i7))
!        end do
!      end do
!      end do
!      end do
!      if (i4.ne.n1*n2*n3     ) call errexec("bad i4")
!      if (i6.ne.n1*n2*n3*ndiv) call errexec("bad i6")
!--
      bvec=reshape((/1d0,0d0,0d0,0d0,1d0,0d0,0d0,0d0,1d0/),(/3,3/))
      if(present(bvec1))bvec(:,1)=bvec1
      if(present(bvec2))bvec(:,2)=bvec2
      if(present(bvec3))bvec(:,3)=bvec3
      call libtetrabzinitialize(bvec,(/n1,n2,n3/),idx1)
      end subroutine
        
      subroutine errexec(smsg)
      character(*), intent(in) :: smsg
      integer :: i(1:1)
      integer :: j
      j=0
      write(*,"(a)")"#ERROR in ttritg_: "//smsg
      i(j)=0
      write(*,"('#errexec cannot stop')")
      write(*,"('#please use debug option')")
      stop
      end subroutine

!!!FROM LIBTETRABZ, MIT License, Copyright (c) 2014 Mitsuaki Kawamura!!!
      subroutine libtetrabzinitialize(bvec,ng,indx1)
      real(8), intent(in   ) :: bvec(3,3)
      integer, intent(in   ) :: ng(3)
      integer, intent(  out) :: indx1(4,*)
      integer :: nk,nt, ivvec(3,4,6), ikv(3)
      integer :: itype, i1, i2, i3, it, ii, divvec(4,4), ivvec0(4)
      real(8) :: l(4), bvec2(3,3), bvec3(3,4)
      nk = product(ng(1:3))
      nt = nk * 6
      do i1 = 1, 3
        bvec2(1:3,i1) = bvec(1:3,i1) / dble(ng(i1))
      end do
      bvec3(1:3,1) = -bvec2(1:3,1) + bvec2(1:3,2) + bvec2(1:3,3)
      bvec3(1:3,2) =  bvec2(1:3,1) - bvec2(1:3,2) + bvec2(1:3,3)
      bvec3(1:3,3) =  bvec2(1:3,1) + bvec2(1:3,2) - bvec2(1:3,3)
      bvec3(1:3,4) =  bvec2(1:3,1) + bvec2(1:3,2) + bvec2(1:3,3)
      ! length of delta bvec
      do i1 = 1, 4
        l(i1) = dot_product(bvec3(1:3,i1),bvec3(1:3,i1))
      end do
      itype = minloc(l(1:4),1)
      ! start & last
      ivvec0(1:4) = (/ 0, 0, 0, 0 /)
      divvec(1:4,1) = (/ 1, 0, 0, 0 /)
      divvec(1:4,2) = (/ 0, 1, 0, 0 /)
      divvec(1:4,3) = (/ 0, 0, 1, 0 /)
      divvec(1:4,4) = (/ 0, 0, 0, 1 /)
      ivvec0(itype) = 1
      divvec(itype, itype) = - 1
      ! Corners of tetrahedra
      it = 0
      do i1 = 1, 3
        do i2 = 1, 3
          if(i2 == i1) cycle
          do i3 = 1, 3
            if(i3 == i1 .or. i3 == i2) cycle
            it = it + 1
            ivvec(1:3,1,it) = ivvec0(1:3)
            ivvec(1:3,2,it) = ivvec(1:3,1,it) + divvec(1:3,i1)
            ivvec(1:3,3,it) = ivvec(1:3,2,it) + divvec(1:3,i2)
            ivvec(1:3,4,it) = ivvec(1:3,3,it) + divvec(1:3,i3)
          end do
        end do
      end do
!     allocate(indx1(20, 6 * nk), indx2(20, 6 * nk), indx3(20 * 6 * nk))
      nt = 0
      do i3 = 1, ng(3)
        do i2  = 1, ng(2)
          do i1 = 1, ng(1)
            do it = 1, 6
              nt = nt + 1
              do ii = 1, 4
                ikv(1:3) = (/i1, i2, i3/) + ivvec(1:3,ii,it) - 1
                ikv(1:3) = modulo(ikv(1:3), ng(1:3))
                indx1(ii,nt) = 1 + ikv(1) + ng(1) * ikv(2) 
     $                                    + ng(1) * ng(2) * ikv(3)
              end do
            end do
          end do
        end do
      end do
      end subroutine
!!!FROM LIBTETRABZ, MIT License, Copyright (c) 2014 Mitsuaki Kawamura!!!
        
      end module
