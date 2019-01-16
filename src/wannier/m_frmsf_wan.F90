MODULE m_frmsf_wan
  !
  IMPLICIT NONE
  !
CONTAINS
  !
SUBROUTINE wrt_frmsf_wan(Na1,Na2,Na3,nkb1,nkb2,nkb3,a1,a2,a3,b1,b2,b3,FermiEnergy,WEIGHT_R,H_MAT_R)
  !
  USE m_rdinput, ONLY : n_occ, dense
  !
  IMPLICIT NONE
  !
  INTEGER,INTENT(IN) :: Na1, Na2, Na3, nkb1, nkb2, nkb3
  REAL(8),INTENT(IN) :: a1(3), a2(3), a3(3), b1(3), b2(3), b3(3), FermiEnergy, &
  &                     WEIGHT_R(-Na1:Na1,-Na2:Na2,-Na3:Na3) 
  COMPLEX(8),INTENT(IN) :: H_MAT_R(n_occ,n_occ,-Na1:Na1,-Na2:Na2,-Na3:Na3) 
  !
  !INTEGER(8) :: ik, i1, i2, i3, ib, jb, fo = 21, nk, i1min, i2min, i3min
  INTEGER :: ik, i1, i2, i3, ib, jb, fo = 21, nk, i1min, i2min, i3min
  REAL(8) :: kvec(3,PRODUCT(dense(1:3))), phase, tpi=2.0d0*acos(-1.0d0), &
  &          Ek(n_occ,PRODUCT(dense(1:3))), proj(n_occ,n_occ,PRODUCT(dense(1:3)))
  COMPLEX(8) :: Hk(n_occ,n_occ), phase_factor, ci=CMPLX(0.0d0,1.0d0)
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
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(nk,n_occ,Na1,Na2,Na3,nkb1,nkb2,nkb3,a1,a2,a3,tpi,ci, &
  !$OMP &        WEIGHT_R,H_MAT_R,kvec,Ek,proj) &
  !$OMP & PRIVATE(ik,ib,jb,i1,i2,i3,i1min,i2min,i3min,phase,Hk,phase_factor)
  !
  Ek(1:n_occ,        1:nk) = 0.0d0
  !
  !$OMP DO
  DO ik = 1, nk
     Hk(1:n_occ,1:n_occ) = 0.0d0
     DO ib = 1, n_occ
        DO jb = 1, n_occ
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
     CALL diagV(n_occ,Hk(1:n_occ,1:n_occ),Ek(1:n_occ,ik))
     !
     proj(1:n_occ,1:n_occ,ik) = DBLE(CONJG(Hk(1:n_occ,1:n_occ)) * Hk(1:n_occ,1:n_occ))
     !
  END DO !ik
  !$OMP END DO
  !$OMP END PARALLEL
  !
  ! Write to file
  !
  DO ib = 1, n_occ
     !
     WRITE(fname,'(a,i0,a)') "./dir-wan/orb", ib, ".frmsf"
     OPEN(fo,FILE=TRIM(fname)) 
     WRITE(fo,*) dense(1:3)
     WRITE(fo,*) 1
     WRITE(fo,*) n_occ 
     WRITE(fo,*) REAL(b1(1:3))
     WRITE(fo,*) REAL(b2(1:3))
     WRITE(fo,*) REAL(b3(1:3))
     DO jb = 1, n_occ
        DO ik = 1, nk
           WRITE(fo,*) REAL(Ek(jb,ik) - FermiEnergy)
        END DO
     END DO
     DO jb = 1, n_occ
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
END SUBROUTINE wrt_frmsf_wan
!
END MODULE m_frmsf_wan
