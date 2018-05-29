subroutine SO3_to_SU2_rotation(SO3_matrix, SU2_matrix)
  implicit none
  real(8),   intent(in   ):: SO3_matrix(3,3)
  complex(8),intent(inout):: SU2_matrix(2,2)
  
  logical:: found=.true., inv_found=.true.
  integer:: ii,jj,kk
  real(8), parameter :: sin3 = 0.866025403784438597d0, &
                        cos3 = 0.5d0
                        
  real(8), parameter :: s2m1 = 0.707106781186547524d0, &
                        s3m1 = 0.57735026918962573d0, &
                        pi   = 3.141592653589793d0, &
                        pi_2 = 1.5707963267948966d0, &
                        pi2_3 = 2.0943951023931953, &
                        pi_3 = 1.0471975511965976d0, &
                        pi_6 = 0.5235987755982988d0
                        
  real(8)::s0(3, 3, 32)
  real(8)::n_vector(3, 32),nx,ny,nz
  real(8)::angles(32),theta_2
  real(8)::tol=1e-8
  complex(8)::ci=(0.0d0,1.0d0)
  
  !from quantum espresso:
  data s0/ 1.d0,  0.d0,  0.d0,  0.d0,  1.d0,  0.d0,  0.d0,  0.d0,  1.d0, &
          -1.d0,  0.d0,  0.d0,  0.d0, -1.d0,  0.d0,  0.d0,  0.d0,  1.d0, &
          -1.d0,  0.d0,  0.d0,  0.d0,  1.d0,  0.d0,  0.d0,  0.d0, -1.d0, &
           1.d0,  0.d0,  0.d0,  0.d0, -1.d0,  0.d0,  0.d0,  0.d0, -1.d0, &
           0.d0,  1.d0,  0.d0,  1.d0,  0.d0,  0.d0,  0.d0,  0.d0, -1.d0, &
           0.d0, -1.d0,  0.d0, -1.d0,  0.d0,  0.d0,  0.d0,  0.d0, -1.d0, &
           0.d0, -1.d0,  0.d0,  1.d0,  0.d0,  0.d0,  0.d0,  0.d0,  1.d0, &
           0.d0,  1.d0,  0.d0, -1.d0,  0.d0,  0.d0,  0.d0,  0.d0,  1.d0, &
           0.d0,  0.d0,  1.d0,  0.d0, -1.d0,  0.d0,  1.d0,  0.d0,  0.d0, &
           0.d0,  0.d0, -1.d0,  0.d0, -1.d0,  0.d0, -1.d0,  0.d0,  0.d0, &
           0.d0,  0.d0, -1.d0,  0.d0,  1.d0,  0.d0,  1.d0,  0.d0,  0.d0, &
           0.d0,  0.d0,  1.d0,  0.d0,  1.d0,  0.d0, -1.d0,  0.d0,  0.d0, &
          -1.d0,  0.d0,  0.d0,  0.d0,  0.d0,  1.d0,  0.d0,  1.d0,  0.d0, &
          -1.d0,  0.d0,  0.d0,  0.d0,  0.d0, -1.d0,  0.d0, -1.d0,  0.d0, &
           1.d0,  0.d0,  0.d0,  0.d0,  0.d0, -1.d0,  0.d0,  1.d0,  0.d0, &
           1.d0,  0.d0,  0.d0,  0.d0,  0.d0,  1.d0,  0.d0, -1.d0,  0.d0, &
           0.d0,  0.d0,  1.d0,  1.d0,  0.d0,  0.d0,  0.d0,  1.d0,  0.d0, &
           0.d0,  0.d0, -1.d0, -1.d0,  0.d0,  0.d0,  0.d0,  1.d0,  0.d0, &
           0.d0,  0.d0, -1.d0,  1.d0,  0.d0,  0.d0,  0.d0, -1.d0,  0.d0, &
           0.d0,  0.d0,  1.d0, -1.d0,  0.d0,  0.d0,  0.d0, -1.d0,  0.d0, &
           0.d0,  1.d0,  0.d0,  0.d0,  0.d0,  1.d0,  1.d0,  0.d0,  0.d0, &
           0.d0, -1.d0,  0.d0,  0.d0,  0.d0, -1.d0,  1.d0,  0.d0,  0.d0, &
           0.d0, -1.d0,  0.d0,  0.d0,  0.d0,  1.d0, -1.d0,  0.d0,  0.d0, &
           0.d0,  1.d0,  0.d0,  0.d0,  0.d0, -1.d0, -1.d0,  0.d0,  0.d0, &
           cos3,  sin3,  0.d0, -sin3,  cos3,  0.d0,  0.d0,  0.d0,  1.d0, &
           cos3, -sin3,  0.d0,  sin3,  cos3,  0.d0,  0.d0,  0.d0,  1.d0, &
          -cos3,  sin3,  0.d0, -sin3, -cos3,  0.d0,  0.d0,  0.d0,  1.d0, &
          -cos3, -sin3,  0.d0,  sin3, -cos3,  0.d0,  0.d0,  0.d0,  1.d0, &
           cos3, -sin3,  0.d0, -sin3, -cos3,  0.d0,  0.d0,  0.d0, -1.d0, &
           cos3,  sin3,  0.d0,  sin3, -cos3,  0.d0,  0.d0,  0.d0, -1.d0, &
          -cos3, -sin3,  0.d0, -sin3,  cos3,  0.d0,  0.d0,  0.d0, -1.d0, &
          -cos3,  sin3,  0.d0,  sin3,  cos3,  0.d0,  0.d0,  0.d0, -1.d0 /
  
  data angles/   0.,    pi,    pi,    pi, &
                 pi,    pi,  pi_2,  pi_2, & 
                 pi,    pi,  pi_2,  pi_2, &
                 pi,    pi,  pi_2,  pi_2, &
              pi2_3, pi2_3, pi2_3, pi2_3, &
              pi2_3, pi2_3, pi2_3, pi2_3, &
               pi_3,  pi_3, pi2_3, pi2_3, &
                 pi,    pi,    pi,    pi /

  
  data n_vector/ 1.d0,  0.d0,  0.d0, &
           0.d0,  0.d0,  1.d0, &
           0.d0,  1.d0,  0.d0, &
           1.d0,  0.d0,  0.d0, &
           s2m1,  s2m1,  0.d0, &
           s2m1, -s2m1,  0.d0, &
           0.d0,  0.d0, -1.d0, &
           0.d0,  0.d0,  1.d0, &
           s2m1,  0.d0,  s2m1, &
          -s2m1,  0.d0,  s2m1, &
           0.d0,  1.d0,  0.d0, &
           0.d0, -1.d0,  0.d0, &
           0.d0,  s2m1,  s2m1, &
           0.d0,  s2m1, -s2m1, &
          -1.d0,  0.d0,  0.d0, &
           1.d0,  0.d0,  0.d0, &
          -s3m1, -s3m1, -s3m1, &
          -s3m1,  s3m1,  s3m1, &
           s3m1,  s3m1, -s3m1, &
           s3m1, -s3m1,  s3m1, &
           s3m1,  s3m1,  s3m1, &
          -s3m1,  s3m1, -s3m1, &
           s3m1, -s3m1, -s3m1, &
          -s3m1, -s3m1,  s3m1, &
           0.d0,  0.d0,  1.d0, &
           0.d0,  0.d0, -1.d0, &
           0.d0,  0.d0,  1.d0, &
           0.d0,  0.d0, -1.d0, &
           sin3, -cos3,  0.d0, &
           sin3,  cos3,  0.d0, &
           cos3, -sin3,  0.d0, &
           cos3,  sin3,  0.d0 /

403 format(3(1x,F10.5,"        ",1x))
404 format(2(1x,F20.5,SP,F20.5,"i        ",1x))

  
  !do kk=1,32
  !  call SU2_rotation(n_vector(1,kk),n_vector(2,kk),n_vector(3,kk),angles(kk)/2.0d0, SU2_matrix)
  !  call SO3_rotation(n_vector(1,kk),n_vector(2,kk),n_vector(3,kk),angles(kk), SO3_matrix)
  !  write(6,*)
  !  write(6,*)
  !  write(6,*)
  !  write(6,404) ((SU2_matrix(ii,jj), jj = 1,2), ii = 1,2)
  !  
  !  call SU2_rotation(n_vector(1,kk),n_vector(2,kk),n_vector(3,kk),angles(kk)/(-2.0d0), SU2_matrix)
  !  write(6,*)
  !  write(6,404) ((SU2_matrix(ii,jj), jj = 1,2), ii = 1,2)
  !  write(6,*)
  !  write(6,403) ((SO3_matrix(ii,jj), jj = 1,3), ii = 1,3)
  !  write(6,*)
  !  write(6,403) ((s0(ii,jj,kk), jj = 1,3), ii = 1,3)
  !enddo
  !stop
  
  do kk=1,32
    !write(6,*) kk
    found=.true.
    inv_found=.true.
    do ii=1,3
      do jj=1,3
        !write(6,*) SO3_matrix(ii,jj), s0(ii,jj,kk)
        if(abs(SO3_matrix(ii,jj)-s0(ii,jj,kk))>tol) found=.false.
        if(abs(SO3_matrix(ii,jj)+s0(ii,jj,kk))>tol) inv_found=.false.
      enddo
    enddo
    if((found==.true.).or.(inv_found==.true.)) then
      exit
      !write(6,*) 'symmetry found.'
    elseif(kk==32) then
      write(6,*) 'ERROR. rotation symmetry not found.'
      stop 
    endif
  enddo
  
  if(found==.true.) then
    call SU2_rotation(n_vector(1,kk),n_vector(2,kk),n_vector(3,kk),angles(kk)/2.0d0, SU2_matrix)
  elseif(inv_found==.true.) then
    call SU2_rotation(n_vector(1,kk),n_vector(2,kk),n_vector(3,kk),-angles(kk)/2.0d0, SU2_matrix)
  else
    write(6,*) 'ERROR. rotation symmetry not found.'
    stop 
  endif
  !do ii=1,3
  !  do jj=1,3
  !    if(abs(SO3_matrix(ii,jj)-s0(ii,jj,kk))<tol) found=.false.
  !    if(abs(SO3_matrix(ii,jj)+s0(ii,jj,kk))<tol) inv_found=.false.
  !  enddo
  !enddo
  
  write(6,*)
  write(6,*)
  write(6,403) ((SO3_matrix(ii,jj), jj = 1,3), ii = 1,3)
  write(6,*)
  write(6,404) ((SU2_matrix(ii,jj), jj = 1,2), ii = 1,2)
  
  return
end subroutine



subroutine SU2_rotation(nx,ny,nz,theta_2, SU2_matrix)
  implicit none
                        
  real(8),intent(in)::nx,ny,nz
  real(8),intent(in)::theta_2
  complex(8),intent(inout):: SU2_matrix(2,2)
  
  real(8)::tol=1e-8
  complex(8)::ci=(0.0d0,1.0d0)

  !theta_2 = angles(kk)/2.0d0
  !nx = n_vector(1,kk)
  !ny = n_vector(2,kk)
  !nz = n_vector(3,kk)
  SU2_matrix(1,1) = cos(theta_2)-ci*nz*sin(theta_2)
  SU2_matrix(1,2) = (-ci*nx-ny)*sin(theta_2) 
  SU2_matrix(2,1) = (-ci*nx+ny)*sin(theta_2) 
  SU2_matrix(2,2) = cos(theta_2)+ci*nz*sin(theta_2)
  
  return
end subroutine



subroutine SO3_rotation(nx,ny,nz,theta, SO3_matrix)
  implicit none
                        
  real(8),intent(in)::nx,ny,nz
  real(8),intent(in)::theta
  real(8),intent(inout):: SO3_matrix(3,3)
  
  real(8)::tol=1e-8
  
  SO3_matrix(1,1) = (1.0d0 - cos(theta))*nx*nx + cos(theta)
  SO3_matrix(1,2) = (1.0d0 - cos(theta))*nx*ny - nz*sin(theta)
  SO3_matrix(1,3) = (1.0d0 - cos(theta))*nx*nz + ny*sin(theta)
  
  SO3_matrix(2,1) = (1.0d0 - cos(theta))*ny*nx + nz*sin(theta)
  SO3_matrix(2,2) = (1.0d0 - cos(theta))*ny*ny + cos(theta)
  SO3_matrix(2,3) = (1.0d0 - cos(theta))*ny*nz - nx*sin(theta)
  
  SO3_matrix(3,1) = (1.0d0 - cos(theta))*nz*nx - ny*sin(theta)
  SO3_matrix(3,2) = (1.0d0 - cos(theta))*nz*ny + nx*sin(theta)
  SO3_matrix(3,3) = (1.0d0 - cos(theta))*nz*nz + cos(theta)
   
  
  return
end subroutine


