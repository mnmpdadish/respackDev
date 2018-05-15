!Copyright (c) 2017 Yusuke Nomura, Ryosuke Akashi, Terumasa Tadano 
program convert101
  use iotk_module
  implicit none
  character(150) :: data_dir 
  character(150) :: data_xml 
  real(8), parameter :: pi = dacos(-1.0d0)

  character(iotk_attlenx) :: attr
  real(8), allocatable :: k_vec(:,:), k_tmp(:,:)
  integer, allocatable :: num_Gk(:)
  integer, allocatable :: Gk_grid(:,:)
  complex(8), allocatable :: evc(:,:)
  complex(8):: c1
  real(8), allocatable :: eig(:)
  character(100), allocatable :: data_wfc1(:), data_wfc2(:), data_Gk(:)
  character(100), allocatable :: data_eig(:)

  real(8) :: a1(3), a2(3), a3(3)
  real(8) :: b1(3), b2(3), b3(3)
  real(8) :: mat_b(3,3)
  real(8) :: mat_b_inv(3,3)
  real(8) :: mat_tmp1(3,3)
  real(8) :: mat_tmp2(3,3)
  
  real(8) :: eFermi, celldm
  real(8) :: Ecut_for_psi, Etot
  integer :: i, j, jj, num_k, num_b
  integer :: ncomp = 1 
  character(5) :: cwrk
  
  logical :: is_SpinOrbit = .false.

  i = command_argument_count()
  if( i /= 1 ) then 
    stop "please specify single directory"
  else 
    call get_command_argument(1,data_dir)
  end if 

  data_xml = trim(data_dir)//"/data-file.xml"
  write(*,'(a)') "Read data from "
  write(*,*) data_xml

  call iotk_open_read(11,data_xml)
  call iotk_scan_begin(11,"BAND_STRUCTURE_INFO",attr=attr)
    call iotk_scan_dat(11,"NUMBER_OF_K-POINTS",num_k)
    call iotk_scan_dat(11,"NUMBER_OF_BANDS",num_b)
    call iotk_scan_dat(11,"FERMI_ENERGY",eFermi)
  call iotk_scan_end(11,"BAND_STRUCTURE_INFO")
  
  call iotk_scan_begin(11,"SPIN",attr=attr)
    call iotk_scan_dat(11,"SPIN-ORBIT_CALCULATION",is_SpinOrbit)
  call iotk_scan_end(11,"SPIN")
  
  if(is_SpinOrbit) ncomp=2 ! else default ncomp=1
  write(6,*) 'is_SpinOrbit=', is_SpinOrbit
  
      
  allocate(k_vec(3,num_k),k_tmp(3,num_k))
  allocate(data_wfc1(num_k),data_wfc2(num_k),data_Gk(num_k),data_eig(num_k))
  call iotk_scan_begin(11,"EIGENVALUES",attr=attr)
  do i = 1, num_k
    !
    ! generate filenames 
    !
    cwrk = "     "
    write(cwrk,'(i5)') i
    do j = 1, 5
      if(cwrk(j:j).eq." ") cwrk(j:j) = "0"
    end do ! j 
    !
    ! generate wfc filename
    !
    if(is_SpinOrbit) then
      data_wfc1(i) = "/K" // cwrk
      data_wfc1(i) = trim(data_wfc1(i)) // "/evc1.dat"
      data_wfc1(i) = trim(data_dir) // trim(data_wfc1(i))
      data_wfc2(i) = "/K" // cwrk
      data_wfc2(i) = trim(data_wfc2(i)) // "/evc2.dat"
      data_wfc2(i) = trim(data_dir) // trim(data_wfc2(i))
    else
      data_wfc1(i) = "/K" // cwrk
      data_wfc1(i) = trim(data_wfc1(i)) // "/evc.dat"
      data_wfc1(i) = trim(data_dir) // trim(data_wfc1(i))
    endif
    
    write(*,*) data_wfc1(i)
    write(*,*) data_wfc2(i)
    !
    ! generate Gk filename
    !
    data_Gk(i) = "/K" // cwrk
    data_Gk(i) = trim(data_Gk(i)) // "/gkvectors.dat"
    data_Gk(i) = trim(data_dir) // trim(data_Gk(i))
    write(*,*) data_Gk(i)
    ! 
    ! generate eigen filename
    !
    data_eig(i) = "/K" // cwrk
    data_eig(i) = trim(data_eig(i)) // "/eigenval.xml"
    data_eig(i) = trim(data_dir) // trim(data_eig(i))
    write(*,*) data_eig(i)

    call iotk_scan_begin(11,"K-POINT"//iotk_index(i),attr=attr)
      call iotk_scan_dat(11,"K-POINT_COORDS",k_vec(:,i))
    call iotk_scan_end(11,"K-POINT"//iotk_index(i))
  end do ! i  
  call iotk_scan_end(11,"EIGENVALUES")

  allocate(num_Gk(num_k))
  call iotk_scan_begin(11,"EIGENVECTORS",attr=attr)
  do i= 1, num_k
    call iotk_scan_begin(11,"K-POINT"//iotk_index(i),attr=attr)
      call iotk_scan_dat(11,"NUMBER_OF_GK-VECTORS",num_Gk(i))
    call iotk_scan_end(11,"K-POINT"//iotk_index(i))
  end do ! i 
  call iotk_scan_end(11,"EIGENVECTORS")
  ! 
  ! Lattice vectors 
  ! 
  call iotk_scan_begin(11,"CELL",attr=attr)

    call iotk_scan_dat(11,"LATTICE_PARAMETER",celldm)

    call iotk_scan_begin(11,"DIRECT_LATTICE_VECTORS",attr=attr)
      call iotk_scan_dat(11,"a1",a1)
      call iotk_scan_dat(11,"a2",a2)
      call iotk_scan_dat(11,"a3",a3)
    call iotk_scan_end(11,"DIRECT_LATTICE_VECTORS")

  call iotk_scan_end(11,"CELL")

  call iotk_scan_begin(11,"PLANE_WAVES",attr=attr)
    call iotk_scan_dat(11,"WFC_CUTOFF",Ecut_for_psi)
  call iotk_scan_end(11,"PLANE_WAVES")
  
  
  

  call iotk_close_read(11)
 
  k_tmp(1,:) = a1(1)*k_vec(1,:) + a1(2)*k_vec(2,:) + a1(3)*k_vec(3,:) 
  k_tmp(2,:) = a2(1)*k_vec(1,:) + a2(2)*k_vec(2,:) + a2(3)*k_vec(3,:) 
  k_tmp(3,:) = a3(1)*k_vec(1,:) + a3(2)*k_vec(2,:) + a3(3)*k_vec(3,:) 
  k_tmp(:,:) = k_tmp(:,:)/celldm
  k_vec(:,:) = k_tmp(:,:)

  open(101,file='./dir-wfn/dat.sample-k') 
    rewind(101) 
    write(101,*) num_k
    do i = 1, num_k
      write(101,*) k_vec(:,i)
    end do 
  close(101)
  deallocate(k_vec,k_tmp)

  open(105,file='./dir-wfn/dat.lattice') 
    rewind(105)
    write(105,*) a1(:)
    write(105,*) a2(:)
    write(105,*) a3(:)
  close(105)
  
  open(117,file='./dir-wfn/dat.bandcalc') 
    rewind(117) 
    Etot = 0d0
    write(117,*) Ecut_for_psi*2d0 ! Hartree => Rydberg  
    write(117,*) eFermi       
    write(117,*) Etot         
  close(117)
 
  open(132,file='./dir-wfn/dat.nkm') 
    rewind(132)
    do i = 1, num_k
      write(132,*) num_Gk(i)
    end do ! i  
  close(132)

  open(104,file='./dir-wfn/dat.kg') 
  rewind(104) 
  do i = 1, num_k
    allocate(GK_grid(3,num_Gk(i)))
    call iotk_open_read(12,data_Gk(i))
      call iotk_scan_dat(12,"GRID",GK_grid)
      write(104,*) num_Gk(i)
      do j = 1, num_Gk(i)
        write(104,*) GK_grid(:,j)
      end do ! j 
    call iotk_close_read(12)
    deallocate(Gk_grid)
  end do ! i 

  open(102,file='./dir-wfn/dat.wfn',FORM='unformatted') 
    rewind(102)
    write(102) ncomp
    do i = 1, num_k
      allocate(evc(num_Gk(i),num_b))
      call iotk_open_read(13, data_wfc1(i))
        do j = 1, num_b
          call iotk_scan_dat(13,"evc"//iotk_index(j),evc(:,j))
        end do ! j  
        do j = 1, num_b
          write(102) evc(:,j)
        end do ! j 
      call iotk_close_read(13)
      
      !Maxime Charlebois
      if(is_SpinOrbit) then
        call iotk_open_read(13, data_wfc2(i))
          do j = 1, num_b
            call iotk_scan_dat(13,"evc"//iotk_index(j),evc(:,j))
          end do ! j  
          do j = 1, num_b
            write(102) evc(:,j)
          end do ! j 
        call iotk_close_read(13)
      endif
      !
      deallocate(evc)
    end do ! i 
  close(102)

!  open(103,file='./dir-wfn/dat.wfn_ascii')
!    rewind(103)
!    write(103,*) ncomp
!    do i = 1, num_k
!      write(103,*) 'i= ', i
!      allocate(evc(num_Gk(i),num_b))
!      call iotk_open_read(13, data_wfc1(i))
!        write(103,*) ' '
!        write(103,*) ' '
!        write(103,*) 'evc1 '
!        do j = 1, 10
!          call iotk_scan_dat(13,"evc"//iotk_index(1),evc(:,1))
!        end do ! j  
!        do j = 1, 10
!          c1 = evc(j,1)/evc(1,1)*abs(evc(1,1))
!          write(103,*) i, j, dsqrt( real(c1)*real(c1)+ imag(c1)*imag(c1) ), real(c1), imag(c1)
!        end do ! j 
!      call iotk_close_read(13)
!      deallocate(evc)
!    end do ! i 
!  close(103)
  deallocate(num_Gk)


  open(111,file='./dir-wfn/dat.eigenvalue') 
    rewind(111)
    write(111,*) num_b 
    do i = 1, num_k
      allocate(eig(num_b))
      call iotk_open_read(14,data_eig(i))
        call iotk_scan_dat(14,"EIGENVALUES",eig)
        do j = 1, num_b
          write(111,*) eig(j)
        end do ! j 
      call iotk_close_read(14)
      deallocate(eig)
    end do ! i  
  close(111)

  call convert_sym(data_xml)

  deallocate(data_wfc1,data_wfc2,data_Gk,data_eig)

end program



subroutine convert_sym(data_xml)
  use iotk_module
  implicit none

  character(150), intent(in) :: data_xml         
  real(8), parameter :: pi = dacos(-1.0d0)
  character(iotk_attlenx) :: attr
  real(8) , allocatable :: ftau(:,:)
  integer , allocatable :: tau(:,:)
  real(8) :: tmp, mat(3,3) 
  integer i, j, k, l, n_sym
  integer :: torf
  integer , allocatable :: mat_sym(:,:)
  
  real(8) :: b1(3), b2(3), b3(3)
  real(8) :: mat_b(3,3)
  real(8) :: mat_b_inv(3,3)
  real(8) :: mat_tmp1(3,3)
  real(8) :: mat_tmp2(3,3)
  integer:: ii,jj,kk,ll

404 format(3(1x,F10.5,"        ",1x))


  open(100,file='./dir-wfn/dat.symmetry') 
  rewind(100) 
  
  
  call iotk_open_read(11,data_xml)
    call iotk_scan_begin(11,"CELL",attr=attr)
    call iotk_scan_begin(11,"RECIPROCAL_LATTICE_VECTORS",attr=attr)
      call iotk_scan_dat(11,"b1",b1)
      call iotk_scan_dat(11,"b2",b2)
      call iotk_scan_dat(11,"b3",b3)
    call iotk_scan_end(11,"RECIPROCAL_LATTICE_VECTORS")
    call iotk_scan_end(11,"CELL")
    
  call iotk_scan_begin(11,"SYMMETRIES",attr=attr)
    call iotk_scan_dat(11,"NUMBER_OF_SYMMETRIES",n_sym)
    allocate(ftau(3,n_sym))
    allocate(tau(3,n_sym))
    allocate(mat_sym(9,n_sym))
    do i = 1, n_sym
      call iotk_scan_begin(11,"SYMM"//iotk_index(i),attr=attr)
        call iotk_scan_dat(11,"ROTATION",mat_sym(:,i))
        call iotk_scan_dat(11,"FRACTIONAL_TRANSLATION",ftau(:,i))
      call iotk_scan_end(11,"SYMM"//iotk_index(i))
    end do ! i 
    

    do i = 1, 1000
      torf = 1
      do j = 1, n_sym
        do k = 1, 3 
          tmp = ftau(k,j)*dble(i) - idnint(ftau(k,j)*dble(i))
          if(abs(tmp)>1.0d-7) torf = 0
        end do ! k 
      end do ! j  
      if( torf == 1 ) then 
        write(100,'(I5)') n_sym
        write(100,'(I5)') i
        do j = 1, n_sym
          do k = 1, 3 
            tau(k,j) = -idnint(ftau(k,j)*dble(i))
          end do ! k  
        end do ! j  
        exit 
      end if 
    end do ! i

    do i = 1, n_sym 
      write(*,*)
      do j = 1, 3 
        write(*,'(9I3)') (mat_sym(k+3*(j-1),i),k=1,3)
      end do 
    end do 

    do i = 1, 3 
      mat_b(i,1)=b1(i)
      mat_b(i,2)=b2(i)
      mat_b(i,3)=b3(i)
      
      mat_b_inv(i,1)=b1(i)
      mat_b_inv(i,2)=b2(i)
      mat_b_inv(i,3)=b3(i)
    end do
    
    write(*,*) 
    write(*,*) 'mat_B'
    write(6,404), ((mat_b_inv(i,k),i=1,3), k=1,3) 
    call inv(3,mat_b_inv)
    write(*,*) 
    write(*,*) 'mat_B_inv'
    write(6,404), ((mat_b_inv(i,k),i=1,3), k=1,3) 
    write(*,*)
    write(*,*)
    write(*,*)
      
    do i = 1, n_sym 
      
      
      write(*,*)
      write(*,*) "S"//iotk_index(i)//"="
      write(6,'(3I3)'), ((mat_sym(ii+3*(kk-1),i),ii=1,3), kk=1,3) 
    
      do ii = 1, 3 
        do jj = 1, 3
          mat_tmp1(ii,jj) = 0.0d0 
          do kk = 1, 3 
            mat_tmp1(ii,jj) = mat_tmp1(ii,jj) + mat_sym(ii+3*(kk-1),i) * mat_b_inv(kk,jj)
          end do
        end do
      end do
      
      do ii = 1, 3 
        do jj = 1, 3 
          mat_tmp2(ii,jj) = 0.0d0
          do kk = 1, 3 
            mat_tmp2(ii,jj) = mat_tmp2(ii,jj) + mat_b(ii,kk) *  mat_tmp1(kk,jj) !mat_b_inv(kk,jj)
          end do
        end do  
      end do
      
      
      do ii = 1, 3 
        do jj = 1, 3 
          mat_tmp1(ii,jj) = 0.0d0
          do kk = 1, 3 
            do ll = 1, 3 
              mat_tmp1(ii,jj) = mat_tmp1(ii,jj) + mat_b(ii,kk) * mat_sym(kk+3*(ll-1),i) * mat_b_inv(ll,jj)
            end do
          end do
        end do  
      end do
      
      
      write(*,*)
      write(*,*) "R"//iotk_index(i)//"="
      do j = 1, 3 
        write(*,404) (mat_tmp2(j,k),k=1,3) 
      end do 
      
      write(*,*)
      write(*,*) "autre facton"//iotk_index(i)//"="
      do j = 1, 3 
        write(*,404) (mat_tmp1(j,k),k=1,3) 
      end do 
      
      
      
      mat = 0.0d0
      l = 0
      !write(*,*)
      do j = 1, 3 
        !write(*,'(9I3)') (mat_sym(k+3*(j-1),i),k=1,3)
        do k = 1, 3 
          l = l+1
          mat(k,j) = dble(mat_sym(l,i))
        end do 
      end do 
      call inv(3,mat) 
      do j = 1, 3 
      do k = 1, 3 
        tmp = mat(k,j) - dble(idnint(mat(k,j)))
        if(tmp>1.0d-10) stop 'Error: something wrong in symmetry matrix'
      end do ! k 
      end do ! j 
      write(100,'(9I3)') ((idnint(mat(k,j)),k=1,3),j=1,3)
      write(100,'(3I3)') tau(:,i)
    end do ! i 
    call iotk_scan_end(11,"SYMMETRIES")
  call iotk_close_read(11)
  close(100)

  deallocate(ftau,tau,mat_sym)

  return

end subroutine


subroutine inv(nm,mat)
  implicit none 
  integer , intent(in) :: nm
  real(8) , intent(inout) :: mat(nm,nm)
  integer :: ipiv(nm)
  integer :: Lwork 
  real(8) , allocatable :: work(:)
  integer :: info 

  Lwork = 10*nm
  allocate (work(Lwork))
  info = 0
  call dgetrf(nm,nm,mat,nm,ipiv,info)
  call dgetri(nm,mat,nm,ipiv,work,Lwork,info)

  if(info /= 0) then
    write(6,*) 'Error: inverse matrix calculation fails, & 
                info (subrouitine inv):' , info
    stop
  end if 
  deallocate(work)

  return 
end subroutine



