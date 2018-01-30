! m_bvector.f90 -- a simple generator of b-vectors
!
! Copyright(C) 2017 Terumasa Tadano
!
! This file is distributed under the terms of the MIT license.
! Please see http://opensource.org/licenses/mit-license.php for information.
!
module m_bvector
  implicit none
contains
  !
  subroutine generate_bvectors(b1, b2, b3, nkb1, nkb2, nkb3, &
       nbmax, maxshell, verbosity, &
       nb, bvec_frac, bvec_cart, weight) 
    !
    ! Automatic generation of b-vectors
    !
    real(8), intent(in) :: b1(3), b2(3), b3(3)
    integer, intent(in) :: nkb1, nkb2, nkb3
    integer, intent(in) :: nbmax, maxshell
    integer, intent(in) :: verbosity ! expect 0 or 1
    integer, intent(out) :: nb
    real(8), intent(out) :: bvec_frac(3,nbmax), bvec_cart(3,nbmax)
    real(8), intent(out) :: weight(nbmax)
    ! local variables
    integer :: i, j, k, m, l
    integer :: nk, nbvec
    integer :: nmulti, nshell
    integer :: nbvec_now, nbvec_shell
    integer :: nrank, nmin, shift
    integer, allocatable :: multiplicity(:)
    integer, allocatable :: index_sort(:), index_k(:)
    integer, allocatable :: index_k_shell(:), index_k_tmp(:)
    logical :: is_new, is_fullrank, satisfy_orthogonal, success
    real(8), parameter :: tolerance = 1.0d-8
    real(8) :: kval(3), dist_tmp, sigma_tmp(6), bvec(3)
    real(8), parameter :: sigma_vec(6) = (/1.d0, 0.d0, 0.d0, 1.d0, 0.d0, 1.d0/)
    real(8), allocatable :: xk(:,:), xk_cart(:,:), dist_k(:)
    real(8), allocatable :: bvec_merged(:,:), bvec_shell(:,:)
    real(8), allocatable :: bvec_tmp(:,:), bvec_now(:,:)
    real(8), allocatable :: bbmat(:,:), weight_tmp(:)
    real(8), allocatable :: Umat(:,:), Sval(:), VTmat(:,:), Sinv(:,:)
    !
    if (verbosity > 0) then
       write(6,*) ' '
       write(6,*) '============================================='
       write(6,*) '     Auatomatic generation of b-vectors      '
       write(6,*) '============================================='
       write(6,*) ' '
    endif
    !
    bvec_frac = 0.0d0
    bvec_cart = 0.0d0
    weight = 0.0d0
    nb = 0
    !
    ! Generate gamma-centered mesh points in fractional coordinate 
    ! in the range of [-0.5:0.5)
    !
    nk = nkb1 * nkb2 * nkb3 + 26
    allocate(xk(3, nk))
    m = 0
    do i = 1, nkb1
       kval(1) = dble(i-1) / dble(nkb1)
       if (kval(1) >= 0.5d0) then
          kval(1) = kval(1) - 1.0d0
       endif
       do j = 1, nkb2
          kval(2) = dble(j-1) / dble(nkb2)
          if (kval(2) >= 0.5d0) then
             kval(2) = kval(2) - 1.0d0
          endif
          do k = 1, nkb3
             kval(3) = dble(k-1) / dble(nkb3)
             if (kval(3) >= 0.5d0) then
                kval(3) = kval(3) - 1.0d0
             endif
             m = m + 1
             xk(1:3, m) = kval(1:3)
          enddo
       enddo
    enddo
    !
    ! Also add nearest G-vectors in each direction
    !
    do i = -1, 1
       do j = -1, 1
          do k = -1, 1
             if (i == 0 .and. j == 0 .and. k == 0) then 
                cycle
             endif
             m = m + 1
             xk(1:3, m) = (/dble(i), dble(j), dble(k)/) 
          enddo
       enddo
    enddo
    !
    ! Convert mesh points in Cartesian coordinate and
    ! calculate distance from gamma point
    !
    allocate(xk_cart(3, nk), dist_k(nk))
    allocate(index_sort(nk))
    do i = 1, nk
       xk_cart(:,i) = xk(1,i)*b1(:) + xk(2,i)*b2(:) + xk(3,i)*b3(:)
       dist_k(i) = sqrt(dot_product(xk_cart(:,i), xk_cart(:,i)))
       index_sort(i) = i
    enddo
    !
    ! Sort distance in ascending order (dist_k(1) <= dist_k(2) <= ...)
    ! 
    call quicksort_double_arg(nk, dist_k, index_sort, 1, nk)
    !
    ! Construct information about nearest-neighbor shells.
    ! multiplicity(i) represents the number k points in the (i-1)th shelll.
    !
    allocate(multiplicity(maxshell+1))
    nmulti = 0
    dist_tmp = 0.0d0
    multiplicity = 0
    m = 1 
    do i = 1, nk
       if (abs(dist_k(i) - dist_tmp) < tolerance) then
          nmulti = nmulti + 1
       else
          dist_tmp = dist_k(i)
          multiplicity(m) = nmulti
          ! Reset value for the next shell
          nmulti = 1
          m = m + 1
          if (m > maxshell + 1) then
             exit
          endif
       endif
    enddo
    !
    nshell = 0
    do i = 1, maxshell + 1
       if (multiplicity(i) > 0) then
          nshell = nshell + 1
       endif
    enddo
    !
    if (verbosity > 0) then
       write(6,*) "  Shell        Distance        Multiplicity"
       write(6,*) " -------      ----------       ------------"
       m = 1
       do i = 1, nshell
          dist_tmp = dist_k(m)
          ! skip 0th shell (distance = 0.0)
          if (dist_tmp > tolerance) then
             write(6,'(i5,F20.8,i18)') i - 1, dist_tmp, multiplicity(i) 
          endif
          m = m + multiplicity(i)
       enddo
       write(6,*)
    endif
    deallocate(dist_k)
    !
    ! index_k: contains k point index of b-vectors 
    ! bvec_now: contains Cartesian coordinates of b-vectors
    !
    allocate(index_k(nk))
    allocate(bvec_now(3,nk))
    nbvec_now = 0
    index_k = 0
    bvec_now = 0.0d0
    !
    ! Main loop over the shell index
    !
    do i = 2, nshell 
       !
       nmulti = multiplicity(i)
       !
       allocate(bvec_shell(3, nmulti)) 
       allocate(bvec_tmp(3, nmulti))
       allocate(index_k_shell(nmulti))
       allocate(index_k_tmp(nmulti))
       bvec_shell = 0.0d0
       index_k_shell = 0
       !
       shift = sum(multiplicity(1:i-1))
       index_k_tmp(:) = index_sort(shift+1:shift+nmulti)
       bvec_tmp(:,1:nmulti) = xk_cart(:,index_k_tmp(1:nmulti))
       !
       ! Search unique set of b-vectors in this shell and
       ! save them in bvec_shell.
       !
       m = 0
       do j = 1, nmulti
          !
          is_new = .true.
          !
          bvec(:) = -bvec_tmp(:,j)
          do k = 1, m
             if (sqrt(dot_product(bvec_shell(:,k)-bvec(:), &
                  bvec_shell(:,k)-bvec(:))) < tolerance) then
                is_new = .false.
                exit
             endif
          enddo
          !
          if (is_new) then
             m = m + 1
             bvec_shell(:,m) = bvec_tmp(:,j)
             index_k_shell(m) = index_k_tmp(j)
          endif
       enddo
       nbvec_shell = m
       nbvec = nbvec_now + nbvec_shell
       deallocate(index_k_tmp, bvec_tmp)
       !
       ! Merge b-vectors for making the b*b matrix
       !
       allocate(bvec_merged(3, nbvec))
       bvec_merged = 0.0d0
       bvec_merged(:,1:nbvec_now) = bvec_now(:,1:nbvec_now)
       bvec_merged(:,nbvec_now+1:nbvec) = bvec_shell(:,1:nbvec_shell)
       !
       if (verbosity > 0) then
          write(6,*)
          write(6,*) " ------------------------------------------------ "
          write(6,'(2x,a,i4)') " New shell index : ", i - 1
          write(6,*)
          write(6,'(2x,a,i4)') " Number of candidate b-vectors : ", nbvec
          write(6,*) "  List of b-vectors (Cartesian) :"
          do j = 1, nbvec
             write(6,'(3F15.8)') bvec_merged(:,j) 
          enddo
       endif
       !
       ! Generate b*b matrix
       !
       allocate(bbmat(6, nbvec))
       allocate(weight_tmp(nbvec))
       !
       m = 0
       do j = 1, 3
          do k = j, 3
             m = m + 1
             do l = 1, nbvec
                bbmat(m, l) = bvec_merged(j,l) * bvec_merged(k,l)
             enddo
          enddo
       enddo
       !
       nmin = min(6, nbvec)
       allocate(Sval(nmin), Umat(6, 6), VTmat(nbvec, nbvec))
       !
       call SVD(6, nbvec, nmin, bbmat, nrank, Sval, Umat, VTmat)
       !
       is_fullrank = (nrank == nmin)
       !
       if (is_fullrank) then
          !
          ! Accept the current shell and update bvec_now and index_k.
          !
          bvec_now(:,nbvec_now+1:nbvec) = bvec_shell(:,1:nbvec_shell)
          index_k(nbvec_now+1:nbvec) = index_k_shell(1:nbvec_shell)
          nbvec_now = nbvec
          !
          allocate(Sinv(6, nbvec))
          Sinv = 0.0d0
          do j = 1, nmin
             Sinv(j,j) = 1.0d0 / Sval(j)
          enddo
          !
          weight_tmp = matmul(transpose(matmul(matmul(Umat, Sinv), VTmat)), sigma_vec)
          !
          deallocate(Sinv) 
          !
          if (verbosity > 0) then
             write(6,*)
             write(6,*) "  Matrix is full-rank. This shell is"
             write(6,*) "  linearly independent on the existing shells."
             write(*,*) "  Weights of each b-vectors below:"
             do j = 1, nbvec
                write(6,'(3x, F15.8)') weight_tmp(j)*0.5d0
             enddo
          endif
          ! 
          ! Check if the b-vectors satisfy the orthogonality relation
          !
          m = 0
          do j = 1, 3
             do k = j, 3
                m = m + 1
                sigma_tmp(m) = 0.0d0
                do l = 1, nbvec
                   sigma_tmp(m) = sigma_tmp(m) &
                        + weight_tmp(l)*bvec_merged(j,l)*bvec_merged(k,l)
                enddo
             enddo
          enddo
          !
          sigma_tmp = sigma_tmp - sigma_vec
          if (sqrt(dot_product(sigma_tmp,sigma_tmp)) < tolerance) then
             satisfy_orthogonal = .true.
          else
             satisfy_orthogonal = .false.
          endif
       else
          if (verbosity > 0) then
             write(6,*)
             write(6,*) "  Matrix is rank-deficient. This shell is"
             write(6,*) "  linearly dependent on the existing shells."
             write(6,*) "  Skipping this shell..."
          endif
       endif
       !
       deallocate(Sval, Umat, VTmat)
       deallocate(bbmat)
       deallocate(bvec_shell,bvec_merged)
       deallocate(index_k_shell)
       !  
       success = .false.
       if (is_fullrank) then
          if (satisfy_orthogonal) then
             !
             success = .true.
             !
             if (verbosity > 0) then
                write(6,*)
                write(6,*) "  Weighted orthogonality is satisfied."
             endif
             !
             ! Store data in the return arrays
             !
             nb = 2 * nbvec
             do j = 1, nbvec
                bvec_frac(:, 2*(j-1) + 1) =  xk(:, index_k(j))
                bvec_frac(:, 2*j        ) = -xk(:, index_k(j))
                bvec_cart(:, 2*(j-1) + 1) =  xk_cart(:, index_k(j))
                bvec_cart(:, 2*j        ) = -xk_cart(:, index_k(j))
                weight(2*(j-1) + 1) = weight_tmp(j)*0.5d0
                weight(2*j        ) = weight_tmp(j)*0.5d0
             enddo
          else
             if (verbosity > 0) then
                write(6,*)
                write(6,*) "  Weighted orthogonality is NOT satisfied."
             endif
          endif
       endif
       !
       deallocate(weight_tmp)
       !
       if (success) exit
    enddo
    !
    deallocate(xk, xk_cart, index_k, bvec_now)
    !
    if (success) then
       if (verbosity > 0) then
          write(6,*)
          write(6,*) " b-vectors are created successfully."
          write(6,*) ' '
          write(6,*) '============================================='
          write(6,*) '        Finish automatic generation          '
          write(6,*) '============================================='
          write(6,*) ' '
       endif
    else
       stop ' Could not find b-vectors'
    endif
    !
    return
  end subroutine generate_bvectors
  !
  !
  subroutine SVD(m, n, k, A, nrank, S, U, VT)
    !
    ! Singular value decomposition of matrix A.
    ! 
    implicit none
    integer, intent(in) :: m, n, k
    real(8), intent(in) :: A(m, n)
    integer, intent(out) :: nrank
    real(8), intent(out) :: S(k), U(m, m), VT(n, n)
    ! local variable
    integer :: i
    integer :: lwork, info
    real(8), allocatable :: work(:)
    real(8), allocatable :: Acopy(:,:)
    !
    lwork = max(3*k+max(m,n), 5*k) 
    lwork = 2 * lwork
    !
    allocate(Acopy(m, n), work(lwork))
    Acopy = A
    !
    call dgesvd('A', 'A', m, n, Acopy, m, S, U, m, &
         VT, n, work, lwork, info)
    !
    nrank = 0
    do i = 1, k
       if (S(1) > 1.0d-12) then
          if (S(i) > S(1)*1.0d-8) then
             nrank = nrank + 1
          endif
       endif
    enddo
    !
    deallocate(Acopy, work)
    !
    return
  end subroutine SVD
  ! 
  recursive subroutine quicksort_double_arg(n, array, ind, first, last)
    !
    ! Quick sort with index (similar to argsort in python)
    !
    implicit none
    integer, intent(in) :: n, first, last
    integer, intent(inout) :: ind(n)
    real(8), intent(inout) :: array(n)
    integer :: i, j, ind_tmp
    real(8) :: tmp, pivot
    pivot= array((first+last)/2)
    i= first
    j= last
    do
       do while (array(i) < pivot)
          i = i + 1
       end do
       do while (pivot < array(j))
          j = j - 1
       end do
       if (i >= j) exit
       ! 
       tmp = array(i)
       array(i) = array(j)
       array(j) = tmp
       ind_tmp = ind(i)
       ind(i) = ind(j)
       ind(j) = ind_tmp
       ! 
       i= i + 1
       j= j - 1
    end do
    if (first < i - 1) call quicksort_double_arg(n, array, ind, first, i - 1)
    if (j + 1 < last)  call quicksort_double_arg(n, array, ind, j + 1, last)
  end subroutine quicksort_double_arg
  !
end module m_bvector
