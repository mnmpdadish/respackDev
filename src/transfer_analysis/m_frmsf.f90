MODULE m_frmsf !(c) Mitsuaki Kawamura 
  !
  IMPLICIT NONE
  !
CONTAINS
  !
!org SUBROUTINE wrt_frmsf_wan(Na1,Na2,Na3,nkb1,nkb2,nkb3,a1,a2,a3,b1,b2,b3,FermiEnergy,WEIGHT_R,H_MAT_R)
SUBROUTINE wrt_frmsf(NWF,dense,Na1,Na2,Na3,nkb1,nkb2,nkb3,a1,a2,a3,b1,b2,b3,FermiEnergy,H_MAT_R)
  !
  !USE m_rdinput, ONLY : n_occ, dense
  !
  IMPLICIT NONE
  !
  integer,intent(in) :: NWF!=n_occ  
  integer,intent(in) :: dense(3)!Dense k-grid for the Wnnier-interpolated FS
  !
  INTEGER,INTENT(IN) :: Na1, Na2, Na3, nkb1, nkb2, nkb3
  REAL(8),INTENT(IN) :: a1(3), a2(3), a3(3), b1(3), b2(3), b3(3), FermiEnergy!, WEIGHT_R(-Na1:Na1,-Na2:Na2,-Na3:Na3) 
  !
  !COMPLEX(8),INTENT(IN) :: H_MAT_R(n_occ,n_occ,-Na1:Na1,-Na2:Na2,-Na3:Na3) 
  COMPLEX(8),INTENT(IN) :: H_MAT_R(NWF,NWF,-Na1:Na1,-Na2:Na2,-Na3:Na3) 
  !
  !INTEGER(8) :: ik, i1, i2, i3, ib, jb, fo = 21, nk, i1min, i2min, i3min
  INTEGER :: ik, i1, i2, i3, ib, jb, fo = 21, nk, i1min, i2min, i3min
  REAL(8) :: kvec(3,PRODUCT(dense(1:3))), phase, tpi=2.0d0*acos(-1.0d0), &
  !&          Ek(n_occ,PRODUCT(dense(1:3))), proj(n_occ,n_occ,PRODUCT(dense(1:3)))
  &          Ek(NWF,PRODUCT(dense(1:3))), proj(NWF,NWF,PRODUCT(dense(1:3)))
  !
  real(8),allocatable :: WEIGHT_R(:,:,:)!WEIGHT_R(-Na1:Na1,-Na2:Na2,-Na3:Na3) 
  !
  !COMPLEX(8) :: Hk(n_occ,n_occ), phase_factor, ci=CMPLX(0.0d0,1.0d0)
  COMPLEX(8) :: Hk(NWF,NWF), phase_factor, ci=CMPLX(0.0d0,1.0d0)
  CHARACTER(256) :: fname
  !
  WRITE(*,*) "Wannnier-interpolated Fermi surface"
  !
  WRITE(*,*) "  Dense grid : ", dense(1:3)
  !
  nk = PRODUCT(dense(1:3))
  !
  ik = 0
  DO i1 = 1, dense(1)
     DO i2 = 1, dense(2)
        DO i3 = 1, dense(3)
           ik = ik + 1
           kvec(1:3,ik) = DBLE((/i1, i2, i3/) - 1) / DBLE(dense(1:3))
        END DO
     END DO
  END DO
  !
  allocate(WEIGHT_R(-Na1:Na1,-Na2:Na2,-Na3:Na3)); WEIGHT_R=1.0d0
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(nk,NWF,Na1,Na2,Na3,nkb1,nkb2,nkb3,a1,a2,a3,tpi,ci, &
  !$OMP &        WEIGHT_R,H_MAT_R,kvec,Ek,proj) &
  !$OMP & PRIVATE(ik,ib,jb,i1,i2,i3,i1min,i2min,i3min,phase,Hk,phase_factor)
  !
  !Ek(1:n_occ,        1:nk) = 0.0d0
  Ek(1:NWF,        1:nk) = 0.0d0
  !
  !$OMP DO
  DO ik = 1, nk
     Hk(1:NWF,1:NWF) = 0.0d0
     DO ib = 1, NWF!n_occ
        DO jb = 1, NWF!n_occ
           DO i1 = -Na1, Na1!-1
              DO i2 = -Na2, Na2!-1
                 DO i3 = -Na3, Na3!-1
                    !--
                    !-- NEAREST R SEARCH
                    CALL search_Rmin(i1,i2,i3,nkb1,nkb2,nkb3, &
                    &                a1,a2,a3,i1min,i2min,i3min)
                    phase = tpi * DOT_PRODUCT(kvec(1:3,ik), DBLE((/i1min, i2min, i3min/)))
                    !--
                    phase_factor = EXP(ci * phase) * WEIGHT_R(i1,i2,i3)
                    Hk(jb,ib) = Hk(jb,ib) &
                    &         + H_MAT_R(jb,ib,i1,i2,i3) * phase_factor
                 END DO!i3
              END DO!i2
           END DO!i1
        END DO!jb
     END DO!ib
     !
     !CALL diagV(n_occ,Hk(1:n_occ,1:n_occ),Ek(1:n_occ,ik))
     CALL diagV(NWF,Hk(1:NWF,1:NWF),Ek(1:NWF,ik))
     !
     !proj(1:n_occ,1:n_occ,ik) = DBLE(CONJG(Hk(1:n_occ,1:n_occ)) * Hk(1:n_occ,1:n_occ))
     proj(1:NWF,1:NWF,ik) = DBLE(CONJG(Hk(1:NWF,1:NWF)) * Hk(1:NWF,1:NWF))
     !
  END DO !ik
  !$OMP END DO
  !$OMP END PARALLEL
  ! 
  !20190421 Kazuma Nakamura 
  !
  WRITE(fname,'(a)')"./dat.frmsf"
  OPEN(fo,FILE=TRIM(fname)) 
  rewind(fo) 
  WRITE(fo,*) dense(1:3)
  WRITE(fo,*) 1
  WRITE(fo,*) NWF 
  WRITE(fo,*) REAL(b1(1:3))
  WRITE(fo,*) REAL(b2(1:3))
  WRITE(fo,*) REAL(b3(1:3))
  DO ib=1,NWF 
     DO ik=1,nk
        WRITE(fo,*) REAL(Ek(ib,ik)-FermiEnergy)
     ENDDO!ik
  ENDDO!ib
  DO ib=1,NWF 
     DO ik=1,nk
        WRITE(fo,*) REAL(ib) 
     ENDDO!ik
  ENDDO!ib
  CLOSE(fo)
  !
  WRITE(*,*) TRIM(fname)," finish" 
  !
  ! Write to file
  !
  DO ib = 1, NWF!n_occ
     !
     !WRITE(fname,'(a,i0,a)')"./dir-wan/orb",ib,".frmsf"
     WRITE(fname,'(a,i3.3,a)')"./dat.orb-",ib,".frmsf"       
     OPEN(fo,FILE=TRIM(fname)) 
     WRITE(fo,*) dense(1:3)
     WRITE(fo,*) 1
     WRITE(fo,*) NWF!n_occ 
     WRITE(fo,*) REAL(b1(1:3))
     WRITE(fo,*) REAL(b2(1:3))
     WRITE(fo,*) REAL(b3(1:3))
     DO jb = 1, NWF!n_occ
        DO ik = 1, nk
           WRITE(fo,*) REAL(Ek(jb,ik) - FermiEnergy)
        END DO
     END DO
     DO jb = 1, NWF!n_occ
        DO ik = 1, nk
           WRITE(fo,*) REAL(proj(ib,jb,ik))
        END DO
     END DO
     CLOSE(fo)
     !
     WRITE(*,*) TRIM(fname), " for orbital ", ib
     !
  END DO
  !
!END SUBROUTINE wrt_frmsf_wan
END SUBROUTINE wrt_frmsf 
!
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
    !
end subroutine search_Rmin 

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
     write(6,*)'ind=',ind 
     stop
    endif 
    !
    deallocate(work_zheevd,rwork_zheevd,iwork_zheevd) 
    return 
end subroutine diagV
!
END MODULE m_frmsf 
