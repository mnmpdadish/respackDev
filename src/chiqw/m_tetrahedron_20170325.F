!
!     Generalized Tetrahedron Method
!
!     (Copyright) 2003 Takeo Fujiwara, Susumu Yamamoto, Yasushi Ishii, and
!     Yoshiro Nohara
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
!
!     For those not following the above, this part of the program package is also
!     released under the MIT License. http://opensource.org/licenses/mit-license.php
!
      module m_tetrahedron
      implicit none
      private

      public :: ttrhdrn_mkidx
      public :: ttrhdrn_chi0
      public :: ttrhdrn_sum
      public :: ttrhdrn_simple
      public :: ttrhdrn_sub1normal
      public :: ttrhdrn_sub1simple
      public :: ttrhdrn_sub1g30g20
      public :: ttrhdrn_sub1g30g21
      public :: ttrhdrn_sub1g30
      public :: ttrhdrn_sub2normal
      public :: ttrhdrn_sub2simple
      public :: ttrhdrn_sub2g30g20
      public :: ttrhdrn_sub2g30g10
      public :: ttrhdrn_sub2g30g21
      public :: ttrhdrn_sub2g30
      public :: ttrhdrn_sub2g10


      integer, parameter, public :: ttrhdrn_vtx=4
      integer, parameter, public :: ttrhdrn_div=6
      integer, parameter :: nvtx=ttrhdrn_vtx
      integer, parameter :: ndiv=ttrhdrn_div
      integer, parameter, public :: ttrhdrn_tbl1(nvtx,ndiv)
     $                       =reshape((/1,2,3,6
     $                                 ,1,3,5,6
     $                                 ,3,5,6,7
     $                                 ,2,3,4,6
     $                                 ,3,4,6,8
     $                                 ,3,6,7,8/)
     $                               ,(/nvtx,ndiv/))

      contains

      subroutine ttrhdrn_simple(dmna,dmnr,n1,n2,n3,imt1,fk,gk,xo)
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
      call ttrhdrn_sum(dmna,dmnr,n1,n2,n3,imt1,fk,gk
     $                ,1,ca2,ca1,xo)
      end subroutine

      subroutine ttrhdrn_sum(dmna,dmnr,n1,n2,n3,imt1,fk,gk,ne,em,ca1,xo)
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
      do i1=1,n1*n2*n3*ndiv
        call calc(dmna,dmnr,imt1(:,i1),n1*n2*n3
     $           ,fk
     $           ,gk,ne,em,ca1
     $           ,xo)
      end do
      end subroutine

      subroutine ttrhdrn_chi0(dmna,dmnr,n1,n2,n3,imt1
     $                       ,e1,e2,f1,f2,delt,dsgn,ne,em,ca1,xo)
      real(8), intent(in) :: dmna,dmnr
      integer, intent(in) :: n1,n2,n3,imt1(nvtx,n1*n2*n3*ndiv)
      complex(8), intent(in) :: e1(n1*n2*n3),e2(n1*n2*n3)
      real(8), intent(in) :: f1,f2,delt,dsgn
      integer, intent(in) :: ne
      complex(8), intent(in) :: em(ne)
      complex(8), intent(  out) :: ca1(ne*nvtx)
      complex(8), intent(inout) :: xo(ne,n1*n2*n3)
      integer :: i1,i2,i3,i4,i5,i6,i7
      real(8) :: ee1(4),ee2(4)
      complex(8) :: ffk(4),ggk(4),xxo(ne,4)
      integer, parameter :: ia1(4)=(/1,2,3,4/),nmax=99
!--- 2013.4.23 
!     real(8), allocatable :: dwgt(:,:,:),dvol(:)
!     allocate(dwgt(4,4,nmax),dvol(nmax))
!---
      real(8) :: dwgt(4,4,nmax),dvol(nmax)
      do i1=1,n1*n2*n3*ndiv
        dwgt(:,:,1)=0.d0; do i2=1,4; dwgt(i2,i2,1)=1.d0; end do
        dvol(1)=1.d0; i2=1
        ee2=dble(e2(imt1(:,i1))); call ttct(f2,ee2,+1,nmax,i2,dwgt,dvol)
        ee1=dble(e1(imt1(:,i1))); call ttct(f1,ee1,-1,nmax,i2,dwgt,dvol)
        do i3=1,i2
          xxo=0.d0; ffk=dvol(i3)
          do i4=1,4
            ggk(i4)=-dcmplx(sum((ee2-ee1)*dwgt(:,i4,i3)),-delt)
          end do
          ffk=ffk*dsgn; ggk=ggk*dsgn
          call calc(dmna,dmnr,ia1,4
     $             ,ffk
     $             ,ggk,ne,em,ca1
     $             ,xxo)
          do i4=1,4; do i5=1,4
            xo(:,imt1(i5,i1))=xo(:,imt1(i5,i1))+xxo(:,i4)*dwgt(i5,i4,i3)
          end do; end do
        end do
      end do
!--- 2013.4.23 
!     deallocate(dwgt,dvol)
!---
      end subroutine

      subroutine ttct(dval,dvtx,imod,nmax,ntet,dwgt,dvol)
      real(8), intent(in) :: dval,dvtx(4)
      integer, intent(in) :: imod,nmax
      integer, intent(inout) :: ntet
      real(8), intent(inout) :: dwgt(4,4,nmax),dvol(nmax)
      logical :: bvtx0(4)
      logical :: b1
      integer :: i1,i2
      real(8) :: d1,d2
      real(8) :: dvtx0(4),dvol0
      integer :: itet
!--- 2013.4.23 
!      logical, allocatable :: bvtx(:,:),btsk(:)
!      allocate(bvtx(4,nmax),btsk(nmax)); bvtx=.false.; btsk=.false.
!---
      logical :: bvtx(4,nmax),btsk(nmax)
      bvtx=.false.
      btsk=.false.
!---
      itet=0
      do while(itet.lt.ntet); itet=itet+1
        do i1=1,4; dvtx0(i1)=sum(dwgt(:,i1,itet)*dvtx(:)); end do
        dvol0=dvol(itet); bvtx0=bvtx(:,itet); b1=.false.
        do i1=   1,3; if(bvtx0(i1))cycle
        do i2=i1+1,4; if(bvtx0(i2))cycle
          b1=(dvtx0(i1)-dval)*(dval-dvtx0(i2)).gt.0.d0; if(b1)exit
        end do; if(b1)exit
        end do
        if(b1)then
!         write(*,*)"divide",dvtx0(i2),dvtx0(i1)
          d1=(dvtx0(i2)-dval)/(dvtx0(i2)-dvtx0(i1))
          d2=(dval-dvtx0(i1))/(dvtx0(i2)-dvtx0(i1))
          ntet=ntet+1; if(ntet.gt.nmax)stop "ntet.gt.nmax"
          bvtx(:,ntet)=bvtx(:,itet); dwgt(:,:,ntet)=dwgt(:,:,itet)
          bvtx(i1,ntet)=.true.
          dwgt(:,i1,ntet)=dwgt(:,i1,itet)*d1+dwgt(:,i2,itet)*d2
          dvol(ntet)=dvol(itet)*d1
          bvtx(i2,itet)=.true.
          dwgt(:,i2,itet)=dwgt(:,i1,itet)*d1+dwgt(:,i2,itet)*d2
          dvol(itet)=dvol(itet)*d2
          itet=itet-1
        else
          do i1=1,4; if(bvtx0(i1))cycle
            btsk(itet)=btsk(itet)
     $        .or.dble(imod)*(dvtx0(i1)-dval).gt.0.d0
!           if(btsk(itet))write(*,"(4g15.7)")dvtx0
!           btsk(itet)=.true.!!!This makes last checking values become 1!!!
          end do
        end if
      end do
!     write(*,*)"ntet=",ntet
!     write(*,*)"btet=",btsk(1:ntet)
      itet=0! dvtx0=0.d0
      do i1=1,ntet; if(.not.btsk(i1))cycle
        itet=itet+1; dwgt(:,:,itet)=dwgt(:,:,i1);dvol(itet)=dvol(i1)
!       write(*,"(4g15.7)")itet,dvol(itet)
!       write(*,"(4g15.7)")dwgt(:,:,itet)
!       dvtx0=dvtx0+dvol(itet)
!    $  *(dwgt(:,1,itet)+dwgt(:,2,itet)+dwgt(:,3,itet)+dwgt(:,4,itet))
      end do; ntet=itet
!     write(*,*)"check",dvtx0
!--- 2013.4.23 
!      deallocate(bvtx,btsk)
!---
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
        call ttrhdrn_sub2normal(
     $             gk(ia4(2)),gk(ia4(1)),gk(ia4(3)),gk(ia4(4))
     $            ,ne,em,ca1(:,1))
        call ttrhdrn_sub2normal(
     $             gk(ia4(1)),gk(ia4(2)),gk(ia4(3)),gk(ia4(4))
     $            ,ne,em,ca1(:,2))
        call ttrhdrn_sub2normal(
     $             gk(ia4(1)),gk(ia4(3)),gk(ia4(2)),gk(ia4(4))
     $            ,ne,em,ca1(:,3))
        call ttrhdrn_sub1normal(
     $             gk(ia4(1)),gk(ia4(4)),gk(ia4(2)),gk(ia4(3))
     $            ,ne,em,ca1(:,4))
          xo(:,ia4(1))=xo(:,ia4(1))+fk(ia4(1))*ca1(:,1)
          xo(:,ia4(2))=xo(:,ia4(2))+fk(ia4(2))*ca1(:,2)
          xo(:,ia4(3))=xo(:,ia4(3))+fk(ia4(3))*ca1(:,3)
          xo(:,ia4(4))=xo(:,ia4(4))+fk(ia4(4))
     $            *(ca1(:,4)-ca1(:,3)-ca1(:,2)-ca1(:,1))
      case(2)
          geq1=sum(gk(ia3(1:2)))/2.d0
          call ttrhdrn_sub2g10   (
     $             geq1      ,geq1      ,gk(ia4(1)),gk(ia4(2))
     $            ,ne,em,ca1(:,1))
          call ttrhdrn_sub2g30   (
     $             geq1      ,gk(ia4(1)),gk(ia4(2)),geq1      
     $            ,ne,em,ca1(:,2))
          call ttrhdrn_sub1g30   (
     $             geq1      ,gk(ia4(2)),gk(ia4(1)),geq1      
     $            ,ne,em,ca1(:,3))
          xo(:,ia3(1))=xo(:,ia3(1))+fk(ia3(1))*ca1(:,1)
          xo(:,ia3(2))=xo(:,ia3(2))+fk(ia3(2))*ca1(:,1)
          xo(:,ia4(1))=xo(:,ia4(1))+fk(ia4(1))*ca1(:,2)
          xo(:,ia4(2))=xo(:,ia4(2))+fk(ia4(2))
     $            *(ca1(:,3)-ca1(:,2)-ca1(:,1)-ca1(:,1))
      case(3)
          geq1=sum(gk(ia3(1:3)))/3.d0
          call ttrhdrn_sub2g30g10(
     $             geq1      ,geq1      ,gk(ia4(1)),geq1      
     $            ,ne,em,ca1(:,1))
          call ttrhdrn_sub1g30g20(
     $             geq1      ,gk(ia4(1)),geq1      ,geq1      
     $            ,ne,em,ca1(:,2))
          xo(:,ia3(1))=xo(:,ia3(1))+fk(ia3(1))*ca1(:,1)
          xo(:,ia3(2))=xo(:,ia3(2))+fk(ia3(2))*ca1(:,1)
          xo(:,ia3(3))=xo(:,ia3(3))+fk(ia3(3))*ca1(:,1)
          xo(:,ia4(1))=xo(:,ia4(1))+fk(ia4(1))
     $            *(ca1(:,2)-ca1(:,1)-ca1(:,1)-ca1(:,1))
      case(4)
        select case(i3)
        case(2)
          geq1=sum(gk(ia3(1:2)))/2.d0; geq2=sum(gk(ia3(3:4)))/2.d0
          call ttrhdrn_sub2g30g21(
     $             geq2      ,geq1      ,geq1      ,geq2      
     $            ,ne,em,ca1(:,1))
          call ttrhdrn_sub1g30g21(
     $             geq1      ,geq2      ,geq2      ,geq1      
     $            ,ne,em,ca1(:,2))
          xo(:,ia3(1))=xo(:,ia3(1))+fk(ia3(1))*ca1(:,1)
          xo(:,ia3(2))=xo(:,ia3(2))+fk(ia3(2))*ca1(:,1)
          xo(:,ia3(3))=xo(:,ia3(3))+fk(ia3(3))
     $                     *(ca1(:,2)-ca1(:,1)-ca1(:,1))*.5d0
          xo(:,ia3(4))=xo(:,ia3(4))+fk(ia3(4))
     $                     *(ca1(:,2)-ca1(:,1)-ca1(:,1))*.5d0
        case(3:6)
          geq1=sum(gk(ia3(1:4)))/4.d0
          call ttrhdrn_sub2simple(
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

      subroutine ttrhdrn_sub1simple(g0,g1,g2,g3,ne,em,x0)
      complex(8), intent(in) :: g0,g1,g2,g3
      integer, intent(in) :: ne
      complex(8), intent(in) :: em(ne)
      complex(8), intent(out) :: x0(ne)
      x0=1.d0/((g0+em)*6.d0)
      end subroutine

      subroutine ttrhdrn_sub1g30g21(g0,g1,g2,g3,ne,em,x0)
      complex(8), intent(in) :: g0,g1,g2,g3
      integer, intent(in) :: ne
      complex(8), intent(in) :: em(ne)
      complex(8), intent(out) :: x0(ne)
      x0=   (.5d0*(g0+g1)+em)/(g0-g1)**2
      x0=x0-   log((g0+em)/(g1+em))*(g0+em)*(g1+em)/(g0-g1)**3
      end subroutine

      subroutine ttrhdrn_sub1g30g20(g0,g1,g2,g3,ne,em,x0)
      complex(8), intent(in) :: g0,g1,g2,g3
      integer, intent(in) :: ne
      complex(8), intent(in) :: em(ne)
      complex(8), intent(out) :: x0(ne)
      x0=   (g0-3.d0*g1-2.d0*em)*.25d0/(g0-g1)**2
      x0=x0+   (g1+em)**2*log((g0+em)/(g1+em))*.5d0/(g0-g1)**3
      end subroutine

      subroutine ttrhdrn_sub1g30   (g0,g1,g2,g3,ne,em,x0)
      complex(8), intent(in) :: g0,g1,g2,g3
      integer, intent(in) :: ne
      complex(8), intent(in) :: em(ne)
      complex(8), intent(out) :: x0(ne)
      x0=   .5d0/(g0-g2)
      x0=x0-   .5d0/((g0-g2)*(g0-g1))**2
     $        *(g0+em)
     $        *(-2.d0*(g1+em)*(g2+em)+(g0+em)*(g1+g2+2.d0*em))
     $        *log((g0+em)/(g1+em))
      x0=x0+   .5d0/(g0-g2)**2
     $        *(g2+em)**2*log((g1+em)/(g2+em))/(g1-g2)
      x0=x0+   .5d0/((g0-g2)*(g0-g1))
     $        *(g1+em)
      end subroutine

      subroutine ttrhdrn_sub1normal(g0,g1,g2,g3,ne,em,x0)
      complex(8), intent(in) :: g0,g1,g2,g3
      integer, intent(in) :: ne
      complex(8), intent(in) :: em(ne)
      complex(8), intent(out) :: x0(ne)
!     x0=   (-g0-g1-2.d0*em)/((g0-g2)*(g0-g3))*.25d0
      x0=  +   (
     $         -(g1+em)**2*log((g1+em)/(g0+em)))
     $        /((g0-g2)*(g0-g3)*(g0-g1))*.5d0
!     x0=x0-   (-g1-g2-2.d0*em)/((g0-g2)*(g2-g3))*.25d0
      x0=x0-   ((g1+em)**2*log((g1+em)/(g0+em))
     $         -(g2+em)**2*log((g2+em)/(g0+em)))
     $        /((g0-g2)*(g1-g2)*(g2-g3))*.5d0
!     x0=x0+   (-g1-g3-2.d0*em)/((g0-g3)*(g2-g3))*.25d0
      x0=x0+   ((g1+em)**2*log((g1+em)/(g0+em))
     $         -(g3+em)**2*log((g3+em)/(g0+em)))
     $        /((g1-g3)*(g0-g3)*(g2-g3))*.5d0
      end subroutine

      subroutine ttrhdrn_sub2simple(g0,g1,g2,g3,ne,em,x0)
      complex(8), intent(in) :: g0,g1,g2,g3
      integer, intent(in) :: ne
      complex(8), intent(in) :: em(ne)
      complex(8), intent(out) :: x0(ne)
      x0=1.d0/((g0+em)*24.d0)
      end subroutine

      subroutine ttrhdrn_sub2g30g20(g0,g1,g2,g3,ne,em,x0)
      complex(8), intent(in) :: g0,g1,g2,g3
      integer, intent(in) :: ne
      complex(8), intent(in) :: em(ne)
      complex(8), intent(out) :: x0(ne)
      x0=   ((g0-g1)*(g0+2.d0*g1+3.d0*em)
     $         -6.d0*(g0+em)*(g1+em))
     $        /(12.d0*(g0-g1)**3)
      x0=x0+   (g0+em)*(g1+em)**2*log((g0+em)/(g1+em))
     $        /(g0-g1)**4*.5d0
      end subroutine

      subroutine ttrhdrn_sub2g30g21(g0,g1,g2,g3,ne,em,x0)
      complex(8), intent(in) :: g0,g1,g2,g3
      integer, intent(in) :: ne
      complex(8), intent(in) :: em(ne)
      complex(8), intent(out) :: x0(ne)
      call ttrhdrn_sub2g30g20(g1,g0,g1,g1,ne,em,x0)
      end subroutine

      subroutine ttrhdrn_sub2g30g10(g0,g1,g2,g3,ne,em,x0)
      complex(8), intent(in) :: g0,g1,g2,g3
      integer, intent(in) :: ne
      complex(8), intent(in) :: em(ne)
      complex(8), intent(out) :: x0(ne)
      x0=   (2.d0*(g0+em)**2-7.d0*(g0+em)*(g2+em)+11.d0*(g2+em)**2)
     $        /(36.d0*(g0-g2)**3)
      x0=x0-   (g2+em)**3*log((g0+em)/(g2+em))
     $        /(6.d0*(g0-g2)**4)
      end subroutine

      subroutine ttrhdrn_sub2g10   (g0,g1,g2,g3,ne,em,x0)
      complex(8), intent(in) :: g0,g1,g2,g3
      integer, intent(in) :: ne
      complex(8), intent(in) :: em(ne)
      complex(8), intent(out) :: x0(ne)
      x0=   (g0+em)*((g0+em)**2+5.d0*(g2+em)*(g3+em)
     $                 -3.d0*(g0+em)*(g2+g3+2.d0*em))
     $        /(12.d0*(g0-g2)**2*(g0-g3)**2)
!     x0=x0+   (g0+em)*(3.d0*(g2+em)*(g3+em)*(g0-g2)*(g0-g3)
!    $                 +((g0+em)*(g2-g3))**2)
!    $        *log(g0+em)/(6.d0*(g0-g2)**3*(g0-g3)**3)
      x0=x0-   (g2+em)**3*log((g2+em)/(g0+em))/(6.d0*(g0-g2)**3*(g2-g3))
      x0=x0+   (g3+em)**3*log((g3+em)/(g0+em))/(6.d0*(g0-g3)**3*(g2-g3))
      end subroutine

      subroutine ttrhdrn_sub2g30   (g0,g1,g2,g3,ne,em,x0)
      complex(8), intent(in) :: g0,g1,g2,g3
      integer, intent(in) :: ne
      complex(8), intent(in) :: em(ne)
      complex(8), intent(out) :: x0(ne)
      x0=   ((g0+em)**2*(g1-g2)+(g1+em)**2*(g0-g2))
     $        /(6.d0*(g0-g2)*(g0-g1)**2*(g1-g2))
      x0=x0-   (g0+em)**2*(2.d0*(g1+em)*(g0-g2)+(g2+em)*(g0-g1))
     $                   *log((g0+em)/(g2+em))
     $        /(6.d0*(g0-g2)**2*(g0-g1)**3)
      x0=x0+   (g1+em)**2*(2.d0*(g0+em)*(g1-g2)-(g2+em)*(g0-g1))
     $                   *log((g1+em)/(g2+em))
     $        /(6.d0*(g0-g1)**3*(g1-g2)**2)
!     x0=x0+   (g2+em)**3*log(g2+em)/(6.d0*(g0-g2)**2*(g1-g2)**2)
      end subroutine

      subroutine ttrhdrn_sub2normal(g0,g1,g2,g3,ne,em,x0)
      complex(8), intent(in) :: g0,g1,g2,g3
      integer, intent(in) :: ne
      complex(8), intent(in) :: em(ne)
      complex(8), intent(out) :: x0(ne)
      x0=-   (g1+em)**2/(6.d0*(g0-g1)*(g1-g2)*(g1-g3))
!     x0=x0+   ((g1+em)*(g2+em)**2*(g1-g3)*(g0-g3)
!    $         +(g1+em)*(g2+em)*(g3+em)*(g1-g2)*(g0-g3)
!    $         +(g0+em)*(g3+em)**2*(g1-g2)**2)*log(g1+em)
!    $        /(6.d0*(g0-g2)*(g0-g3)*(g1-g2)**2*(g1-g3)**2)
      x0=x0-   (g2+em)**3*log((g2+em)/(g1+em))
     $        /(6.d0*(g0-g2)*(g1-g2)**2*(g2-g3))
      x0=x0+   (g3+em)**3*log((g3+em)/(g1+em))
     $        /(6.d0*(g0-g3)*(g1-g3)**2*(g2-g3))
      x0=x0+   (g0+em)**3*log((g0+em)/(g1+em))
     $        /(6.d0*(g0-g1)**2*(g0-g2)*(g0-g3))
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
          d1=abs(gk(ia1(mod(i1,4)+1))-gk(ia1(mod(i2,4)+1)))
          d2=abs(gk(ia1(mod(i1,4)+1))+gk(ia1(mod(i2,4)+1)))*.5d0
           if(d1.lt.dmna.or.d1/d2.lt.dmnr)then
!          if(d1.lt.dmna.and.d1/d2.lt.dmnr)then
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

      subroutine ttrhdrn_mkidx(n1,n2,n3,idx1,bvec1,bvec2,bvec3)
      integer, intent(in) :: n1
      integer, intent(in) :: n2
      integer, intent(in) :: n3
      integer, intent(out):: idx1(nvtx,n1*n2*n3*ndiv)
      real(8), intent(in), optional :: bvec1(3),bvec2(3),bvec3(3)
      real(8) :: bvec(3,3)
      integer :: i1,i2,i3,i4,i5,i6,i7
      integer :: imt1(8,8)
      integer :: ia1(8)
      !
      integer :: nvec(3)!20180413 Y.Nohara 
      !
      bvec=reshape((/1d0,0d0,0d0,0d0,1d0,0d0,0d0,0d0,1d0/),(/3,3/))
      if(present(bvec1))bvec(:,1)=bvec1
      if(present(bvec2))bvec(:,2)=bvec2
      if(present(bvec3))bvec(:,3)=bvec3
      !
      !20180413 Y.Nohara    
      !
      !(/n1,n2,n3/) gives a warning at run time
      !
      !call libtetrabzinitialize(bvec,(/n1,n2,n3/),idx1)
      !
      nvec(1)=n1; nvec(2)=n2; nvec(3)=n3
      call libtetrabzinitialize(bvec,nvec(1),idx1)
      !
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
        
      subroutine errexec(smsg)
      character(*), intent(in) :: smsg
      integer :: i(1:1)
      integer :: j
      j=0
      write(*,"(a)")"#ERROR in ttrhdrn_ : "//smsg
      i(j)=0
      write(*,"('#errexec cannot stop')")
      write(*,"('#please use debug option')")
      stop
      end subroutine

      end module
