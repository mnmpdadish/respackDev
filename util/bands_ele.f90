      program branch_phonon
      implicit none
      character(5)  gomi1,gomi2,gomi3
      !input
      character(100) :: inputfile
      !/input
      integer  nbnd, nks   
      integer  nk_dummy
      character(10) char_dummy
      integer  nk,i,ibnd,i_ene
      logical fileExists
 
      real(8),allocatable :: branch_q(:,:),abs_q(:),branch_w(:,:)
      real(8),allocatable :: gam(:,:,:)
      real(8),allocatable :: broad(:)
      real(8) ene
      real(8)  del_ene 
      !
      !read phonon branches
      !

      if(iargc()/=1) then
        write(6,*) "hey, must have one argument."
        write(6,*) 
        write(6,*) "example of use:"
        write(6,*) "/a.out Cu.bands"
        stop
      endif

      call getarg(1, inputfile)
      inquire( file=inputfile, EXIST=fileExists ) 
      if(.not.fileExists) then
        write(6,*) "error: file does not exist."
        write(6,*) "terminated."
        stop 
      endif 

      open(unit=100,file=trim(inputfile),status="unknown")
      read(100,*)gomi1,gomi2,nbnd,gomi3,nks

      write(*,*)"gomi1=",gomi1,"  gomi2=",gomi2,"  gomi3=",gomi3
      write(*,*)"nbnd=",nbnd,"nks=",nks

      allocate(branch_q(nks,3))
      allocate(branch_w(nks,nbnd))


      DO nk = 1,nks 
        read(100,*) (branch_q(nk,i),i=1,3)
        read(100,*) (branch_w(nk,ibnd),ibnd=1,nbnd)
      ENDDO

      !write(*,*)"branch_q(2,*)=",(branch_q(2,i),i=1,3)
      !write(*,*)"branch_w(2,*)=",(branch_w(2,ibnd),ibnd=1,nbnd)

      close(100)
      !
      !read gammas
      !
    
      allocate(abs_q(nks))

      abs_q(1) = 0d0
      DO nk=2 , nks
        abs_q(nk)=abs_q(nk-1) &
                   + sqrt( ( branch_q(nk,1)-branch_q(nk-1,1) )**2d0 &
           &              +( branch_q(nk,2)-branch_q(nk-1,2) )**2d0 &
           &              +( branch_q(nk,3)-branch_q(nk-1,3) )**2d0 )   
      ENDDO


      open(unit=300,file="out.dat",status="unknown") 
      DO nk=1,nks
        write(300,'(f16.8,100f16.8)')abs_q(nk),                         &
     &                        (branch_w(nk,ibnd),ibnd=1,nbnd)  
      ENDDO

      close(300)

      deallocate(branch_q)
      deallocate(branch_w)
      deallocate(abs_q)

      call bands_sort(nbnd,nks)

      end program



       subroutine bands_sort(N_band,N_data)
       implicit none
       integer, intent(in) :: N_band 
       integer, intent(in) :: N_data 
       real(8), allocatable :: k_contour(:), ene(:,:)

       real(8) swap
       integer i,j,k

      
       allocate(k_contour(N_data), ene(N_band,N_data))
       open(unit=11, file="out.dat", status="unknown")
       DO i=1, N_data
         read(11,*)k_contour(i), ene(:,i)
       ENDDO
       close(11)

       ! swap
       DO i=1, N_data
         DO j=1,N_band-1
           DO k=1, N_band -j
             IF(ene(k,i).gt.ene(k+1,i))THEN
               swap = ene(k,i)
               ene(k,i)= ene(k+1,i)
               ene(k+1,i) = swap
             ENDIF
           ENDDO
         ENDDO
       ENDDO

       !/swap

       open(unit=12,file="out_sort.dat",status="unknown")
       do j = 1, N_band
       DO i=1, N_data
         write(12,'(f16.8,f16.8)')                                      &
     &       k_contour(i)/k_contour(N_data), ene(j,i)
       ENDDO
       write(12,*)
       end do 
       close(12)
       return
       end subroutine

