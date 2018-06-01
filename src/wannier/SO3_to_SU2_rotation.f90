subroutine SO3_to_SU2_rotation(b1,b2,b3,nsymq,rg,rinv_SO)
  implicit none
  real(8), intent(in)::b1(3),b2(3),b3(3)
  integer, intent(in)::nsymq
  integer, intent(in)::rg(3,3,nsymq)
  complex(8),intent(inout)::rinv_SO(2,2,nsymq) 

  integer:: ii,jj,kk,ll, iop
  real(8)::mat_b(3,3)
  real(8)::mat_b_inv(3,3)
  real(8)::SO3_matrix(3,3)
  complex(8)::SU2_matrix(2,2)
  
  do ii = 1, 3 
    mat_b(ii,1)=b1(ii)
    mat_b(ii,2)=b2(ii)
    mat_b(ii,3)=b3(ii)
  end do
      
  mat_b_inv=mat_b
  call invmat(3,mat_b_inv(1,1))
    
  do iop=1,nsymq
   do ii = 1, 3 
    do jj = 1, 3 
     SO3_matrix(ii,jj) = 0.0d0
     do kk = 1, 3 
      do ll = 1, 3   ! double matrix multiplication
       SO3_matrix(ii,jj) = SO3_matrix(ii,jj) + mat_b(ii,kk)*rg(kk,ll,iop)*mat_b_inv(ll,jj)
      end do
     end do
    end do
   end do
   
   !call SO3_to_SU2_rotation(SO3_matrix(1,1), SU2_matrix(1,1))
   call from_SO3_matrix_to_SU2_matrix(SO3_matrix,SU2_matrix)
   
   write(6,*)
   write(6,*)
   write(6,fmt='(3(1x,F10.5,"      ",1x))') ((SO3_matrix(ii,jj), jj = 1,3), ii = 1,3)
   write(6,*)
   write(6,fmt='(2(1x,F20.2,SP,F20.2,"i   ",1x))') ((SU2_matrix(ii,jj), jj = 1,2), ii = 1,2)
   rinv_SO(:,:,iop)=SU2_matrix(:,:)
  enddo
  
  
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


subroutine from_SO3_matrix_to_SU2_matrix(SO3_matrix,SU2_matrix)
  implicit none
                        
  real(8),intent(in):: SO3_matrix(3,3)
  complex(8),intent(out)::SU2_matrix(2,2)
  
  real(8):: s(3,3),det
  
  real(8)::nx,ny,nz
  real(8)::theta,sintheta,costheta
  real(8)::tol=1e-3,diff
  
  real(8)::matrix1(3,3),matrix2(3,3), Imatrix(3,3)
  integer::ii,jj,kk,nn
  real(8),parameter::pi=dacos(-1.0d0)
  logical::isInv=.false.
      
  
  s=SO3_matrix
  det = s(1,1) * ( s(2,2) * s(3,3) - s(3,2) * s(2,3) )-   &
        s(1,2) * ( s(2,1) * s(3,3) - s(3,1) * s(2,3) )+   &
        s(1,3) * ( s(2,1) * s(3,2) - s(3,1) * s(2,2) )
  
  if(abs(det+1.d0) < tol) then
    s=-SO3_matrix
    isInv=.true.
  endif
  
  
  Imatrix(:,:) = 0.0d0 
  do ii=1,3
    Imatrix(ii,ii) = 1.0d0
  enddo
  
  matrix1 = s
  matrix2 = s
  
  do nn=1,10
    
    diff=0.0d0
    do ii=1,3
      do jj=1,3
        diff = diff + abs(matrix2(ii,jj) - Imatrix(ii,jj))
      enddo
    enddo
    
    if(diff < tol) exit
    matrix2(:,:) = 0.0d0 
    
    do ii=1,3
      do jj=1,3
        do kk=1,3
          matrix2(ii,jj) = matrix2(ii,jj) + matrix1(ii,kk)*s(kk,jj)
        enddo
      enddo
    enddo
    
    !write(6,*) 'nn=', nn, diff
    !write(6,*)
    !write(6,fmt='(3(1x,F10.5,"        ",1x))') ((SO3_matrix(ii,jj), jj = 1,3), ii = 1,3)
    !write(6,*)
    !write(6,fmt='(3(1x,F10.5,"        ",1x))') ((matrix2(ii,jj), jj = 1,3), ii = 1,3)
    !write(6,*) 
    !write(6,fmt='(3(1x,F10.5,"        ",1x))') ((matrix1(ii,jj), jj = 1,3), ii = 1,3)
    !write(6,*) 
    
    matrix1(:,:)=matrix2(:,:)
  enddo
  
  theta = pi*2.0d0/(nn)
  
  !write(6,*)
  !write(6,*) 'N=', nn, theta
  !SO3_matrix(1,1)-SO3_matrix(1,1)
  
  sintheta = sin(theta)
  costheta = cos(theta)
  
  if(nn==1) then
    nx = 0
    ny = 0 
    nz = 0  !does not matter
    theta=0.0d0 !matter

  elseif(nn==2) then !more complicated
    matrix2(:,:) = s(:,:)
    do ii=1,3
      do jj=1,3
        matrix2(ii,jj) = 0.5d0*(matrix2(ii,jj) + Imatrix(ii,jj))
      enddo
    enddo
    if(matrix2(1,1)>tol) then
      nx = sqrt(matrix2(1,1)) ! choose positive (can for 180 angle)
      ny = matrix2(1,2) / nx
      nz = matrix2(1,3) / nx
    elseif(matrix2(2,2)>tol) then
      ny = sqrt(matrix2(2,2)) ! same
      nx = matrix2(2,1) / ny
      nz = matrix2(2,3) / ny    
    elseif(matrix2(3,3)>tol) then
      nz = sqrt(matrix2(3,3)) ! same
      nx = matrix2(3,1) / nz
      ny = matrix2(3,2) / nz        
    else
      write(6,*) 'error with SO3 matrix'
      write(6,fmt='(3(1x,F10.5,"        ",1x))') ((s(ii,jj), jj = 1,3), ii = 1,3)
      write(6,*) 'error with SO3 matrix'
      write(6,fmt='(3(1x,F10.5,"        ",1x))') ((matrix2(ii,jj), jj = 1,3), ii = 1,3)
      stop
    endif
        
  else !nn>2
  
    nx = 0.5d0*(s(3,2) - s(2,3)) / sintheta
    ny = 0.5d0*(s(1,3) - s(3,1)) / sintheta
    nz = 0.5d0*(s(2,1) - s(1,2)) / sintheta
    
  endif
  
  !write(6,*) 'n=', nx,ny,nz
  
  if(isInv) then
    call SU2_rotation(nx,ny,nz,theta/2.0d0,SU2_matrix)
  else
    call SU2_rotation(nx,ny,nz,theta/2.0d0,SU2_matrix)
  endif
  
  return
end subroutine


