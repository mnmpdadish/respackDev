module m_rdinput 
  implicit none
  private
  public::read_input
  !--
  !&param-wannier 
  integer,public::icell!0:gen,1:sc,2:fcc,3:bcc,4:hcp,5:ttr&ort,6:mon,7:tri,8:bct  
  integer,public::n_occ!Total number of considerd band in wannier calc 
  integer,public::N_wannier!Total number of considerd band in wannier calc 
  integer,public::nigs!TOTAL NUMBER OF INITIAL GUESS GAUSSIAN       
  integer,public::N_initial_guess!TOTAL NUMBER OF INITIAL GUESS GAUSSIAN       
  real(8),public::EPS_SPILLAGE!Threshold for spilage minimization 
  real(8),public::DAMP!damping in spilage minimization 
  real(8),public::EPS_SPREAD!Threshold for spread minimization 
  real(8),public::MAX_STEP_LENGTH!Step length in line search 
  real(8),public::Lower_energy_window!E_LOWER!LOWER BOUND OF ENERGY WINDOW (eV)
  real(8),public::E_LOWER!LOWER BOUND OF ENERGY WINDOW (eV)
  real(8),public::Upper_energy_window!E_UPPER!UPPER BOUND OF ENERGY WINDOW (eV)
  real(8),public::E_UPPER!UPPER BOUND OF ENERGY WINDOW (eV)
  logical,public::set_inner_window!flag for inner window 
  real(8),public::Upper_inner_window!Upper inner energy window for wannier calc (eV)
  real(8),public::E_UPPER_inner!Upper inner energy window for wannier calc (eV)
  real(8),public::Lower_inner_window!Lower inner energy window for wannier calc (eV)
  real(8),public::E_LOWER_inner!Lower inner energy window for wannier calc (eV)
  integer,public::flg_BMAT!0:BMAT=unit matrix, 1:BMAT from input 
  integer,public::flg_initial_guess_direc!0:use global coordinate, 1:use local coordinate  
  real(8),public::electron_number_wannier_space!electron number in wannier_space for density matrix 
  real(8),public::elec_num!electron number in wannier_space for density matrix 
  integer,public::flg_fermisurface!flag for fermisurface  
  logical,public::CALC_FERMISURFACE!flag for fermisurface  
  integer,public::flg_global_dos!flag for global dos 
  logical,public::CALC_GLOBAL_DOS!flag for global dos 
  !hidden parameters 
  integer,public::flg_vis_bloch!flag for visualization  of Bloch function
  logical,public::CALC_REAL_SPACE_BLOCH!flag for visualization of Bloch function 
  logical,public::CALC_BLOCH_VIS_RANGE!flag to calculate visualization range of Bloch function 
  real(8),public::calc_k(3)!k points to be calculated 
  integer,public::reading_bmat_format!To specify bmat reading format  
  !&param-interpolation   
  integer,public::N_sym_points!The number of k-point points in symmetry line
  integer,public::Ndiv!Separation between symmetry points  
  integer,public::reading_sk_format!20170709 
  integer,public::flg_3d_fatband!flag for color plot of fat band 
  logical,public::CALC_3D_FATBAND!flag for color plot of fat band 
  !&param-visualization   
  integer,public::flg_vis_wannier!flag for visualization 
  logical,public::CALC_REAL_SPACE_WANNIER!flag for visualization 
  integer,public::N_write_wannier!Total number of visalization wannier calculation (manual mode)
  integer,public::ix_vis_min!x-range for visualization 
  integer,public::ix_vis_max!x-range for visualization 
  integer,public::iy_vis_min!y-range for visualization 
  integer,public::iy_vis_max!y-range for visualization 
  integer,public::iz_vis_min!z-range for visualization 
  integer,public::iz_vis_max!z-range for visualization 
  integer,public::xmin!x-range for visualization 
  integer,public::xmax!x-range for visualization 
  integer,public::ymin!y-range for visualization 
  integer,public::ymax!y-range for visualization 
  integer,public::zmin!z-range for visualization 
  integer,public::zmax!z-range for visualization
  integer,public::dense(3)!Dense k-grid for the Wnnier-interpolated FS
  !initial_guess 
  type initial_guess 
  !integer::i,j!20170406
  character(10)::orb!20170406
  real(8)::a,x,y,z 
  real(8)::lx(3)!20170914
  real(8)::ly(3)!20170914
  real(8)::lz(3)!20170914
  complex(8)::s_up(2)!20190210 Maxime Charlebois
  complex(8)::s_dn(2)!20190210 Maxime Charlebois
  end type initial_guess  
  type(initial_guess),public,allocatable::vec_ini(:) 
  !BMAT  
  real(8),public,allocatable::B_MAT(:,:)!B_MAT(nigs,n_occ)
  !symmetry line in BZ 
  real(8),public,allocatable::SK_sym_pts(:,:) 
  !wrt_list 
  integer,public,allocatable::wrt_list(:)!wrt_list(N_write_wannier) 
  !--
  namelist/param_wannier/icell,N_wannier,N_initial_guess,EPS_SPILLAGE,DAMP,EPS_SPREAD,MAX_STEP_LENGTH,&
  Lower_energy_window,Upper_energy_window,set_inner_window,Upper_inner_window,Lower_inner_window,flg_BMAT,&
  reading_bmat_format,flg_initial_guess_direc,flg_vis_bloch,calc_k,electron_number_wannier_space,&
  flg_fermisurface,flg_global_dos
  namelist/param_interpolation/N_sym_points,Ndiv,reading_sk_format,dense,&!20170709  
  flg_3d_fatband
  namelist/param_visualization/flg_vis_wannier,N_write_wannier,& 
  ix_vis_min,ix_vis_max,iy_vis_min,iy_vis_max,iz_vis_min,iz_vis_max
contains
  !
  subroutine read_input
    use m_SO3_to_SU2_Local, only: make_SU2_from_local_axis 
    implicit none
    integer::igs,ix,ik  
    !
    !20190611 Kazuma Nakamura 
    !
    !SO3_to_SU2_Local 
    real(8),allocatable::local(:,:)!local(3,3)
    complex(8),allocatable::SU2(:,:)!SU2(2,2)
    !
    !--
    !&param_wannier
    !--
    !default
    ICELL=0
    EPS_SPILLAGE=1.0d-4
    DAMP=0.1d0 
    EPS_SPREAD=1.0d-4
    MAX_STEP_LENGTH=4.0d0 
    SET_INNER_WINDOW=.FALSE.
    UPPER_INNER_WINDOW=0.0d0 
    LOWER_INNER_WINDOW=0.0d0 
    FLG_BMAT=0
    READING_BMAT_FORMAT=0
    FLG_INITIAL_GUESS_DIREC=0
    electron_number_wannier_space=0.0d0 
    flg_fermisurface=0 
    flg_global_dos=0 
    flg_vis_bloch=0
    CALC_BLOCH_VIS_RANGE=.FALSE.
    calc_k(:)=0.0d0
    !---
    !open(999,file='input.in')
    !read(999,nml=param_wannier)
    read(5,nml=param_wannier)
    !--
    n_occ=N_wannier 
    nigs=N_initial_guess 
    E_LOWER=Lower_energy_window
    E_UPPER=Upper_energy_window
    E_LOWER_inner=Lower_inner_window
    E_UPPER_inner=Upper_inner_window
    elec_num=electron_number_wannier_space 
    !--
    !check ENERGY WINDOW
    if(E_LOWER>=E_UPPER)then 
     write(6,'(a)')'wrong input: energy windows; stop'
     write(6,'(a)')'Lower_energy_window >= Upper_energy_window: WRONG'
     write(6,'(a,2f20.9)')'Lower_energy_window, Upper_energy_window',Lower_energy_window,Upper_energy_window
     stop
    endif 
    !--
    !default set_inner_window==F
    !if(set_inner_window.eqv..true..and.E_UPPER_inner==0.0d0.and.E_LOWER_inner==0.0d0)then!20171212 
    !if(set_inner_window==.true..and.E_UPPER_inner==0.0d0.and.E_LOWER_inner==0.0d0)then!20180822 
    if(set_inner_window.and.E_UPPER_inner==0.0d0.and.E_LOWER_inner==0.0d0)then!20181009
     write(6,'(a)')'wrong input: inner energy windows; stop'
     write(6,'(a)')'Lower_inner_window = Upper_inner_window = 0.0d0: WRONG'
     write(6,'(a,2f20.9)')'Upper_inner_window,Lower_inner_window',Upper_inner_window,Lower_inner_window
     stop
    endif 
    !if(set_inner_window.eqv..true.)then!20171212 
    !if(set_inner_window==.true.)then!20180822 
    if(set_inner_window)then!20181009 
     if(E_LOWER_inner>=E_UPPER_inner)then 
      write(6,'(a)')'wrong input: inner energy windows; stop'
      write(6,'(a)')'Lower_inner_window >= Upper_inner_window: WRONG'
      write(6,'(a,2f20.9)')'Lower_inner_window,Upper_inner_window',Lower_inner_window,Upper_inner_window
      stop
     endif 
     if(E_UPPER<=E_UPPER_inner)then 
      write(6,'(a)')'wrong input: inner energy windows; stop'
      write(6,'(a)')'Upper_energy_window <= Upper_inner_window: WRONG'
      write(6,'(a,2f20.9)')'Upper_energy_window,Upper_inner_window',Upper_energy_window,Upper_inner_window
      stop
     endif  
     if(E_LOWER>=E_LOWER_inner)then 
      write(6,'(a)')'wrong input: inner energy windows; stop'
      write(6,'(a)')'Lower_energy_window >= Lower_inner_window: WRONG'
      write(6,'(a,2f20.9)')'Lower_energy_window,Lower_inner_window',Lower_energy_window,Lower_inner_window
      stop
     endif 
    endif 
    !--
    !INITIAL GUESS
    allocate(vec_ini(nigs)) 
    !--
    if(flg_initial_guess_direc==0)then!default 20170914 
     write(6,'(a)')'DIRECTION OF INITIAL GUESS IS REPRESENTED IN GLOBAL COORDINATE' 
     do igs=1,nigs
     !read(999,*)vec_ini(igs)
     !read(999,*)vec_ini(igs)%orb,vec_ini(igs)%a,vec_ini(igs)%x,vec_ini(igs)%y,vec_ini(igs)%z
      read(5,*)vec_ini(igs)%orb,vec_ini(igs)%a,vec_ini(igs)%x,vec_ini(igs)%y,vec_ini(igs)%z
     enddo 
     do igs=1,nigs 
      vec_ini(igs)%lx(:)=0.0d0  
      vec_ini(igs)%ly(:)=0.0d0  
      vec_ini(igs)%lz(:)=0.0d0  
      vec_ini(igs)%lx(1)=1.0d0  
      vec_ini(igs)%ly(2)=1.0d0  
      vec_ini(igs)%lz(3)=1.0d0  
      vec_ini(igs)%s_up(:)=0.0d0  
      vec_ini(igs)%s_dn(:)=0.0d0  
      vec_ini(igs)%s_up(1)=1.0d0  
      vec_ini(igs)%s_dn(2)=1.0d0  
     enddo 
    endif 
    !--
    if(flg_initial_guess_direc==1)then! 20170914 
     write(6,'(a)')'DIRECTION OF INITIAL GUESS IS REPRESENTED IN LOCAL COORDINATE' 
     do igs=1,nigs
      !read(999,*)vec_ini(igs)
      !read(5,*)vec_ini(igs)
       read(5,*)vec_ini(igs)%orb,vec_ini(igs)%a,vec_ini(igs)%x,vec_ini(igs)%y,vec_ini(igs)%z,&
               &vec_ini(igs)%lx,vec_ini(igs)%ly,vec_ini(igs)%lz
     enddo 
     !
     !20190613 Kazuma Nakamura
     !
     !do igs=1,nigs 
     ! vec_ini(igs)%s_up(:)=0.0d0  
     ! vec_ini(igs)%s_dn(:)=0.0d0  
     ! vec_ini(igs)%s_up(1)=1.0d0  
     ! vec_ini(igs)%s_dn(2)=1.0d0  
     !enddo 
     !
     allocate(local(3,3));local=0.0d0 
     allocate(SU2(2,2));SU2=0.0d0 
     do igs=1,nigs 
       local(1,:)=vec_ini(igs)%lx(:)
       local(2,:)=vec_ini(igs)%ly(:)
       local(3,:)=vec_ini(igs)%lz(:)
       call make_SU2_from_local_axis(3,3,local(1,1),SU2(1,1)) 
       vec_ini(igs)%s_up(:)=SU2(:,1) 
       vec_ini(igs)%s_dn(:)=SU2(:,2) 
     enddo 
     deallocate(local,SU2) 
     !
    endif 
    !--
    if(flg_initial_guess_direc==2)then 
     write(6,'(a)')'DIRECTION OF INITIAL GUESS IS REPRESENTED IN LOCAL COORDINATE'
     write(6,'(a)')'ROTATION OF THE SPIN QUANTIZATION AXIS' 
     do igs=1,nigs
      !read(999,*)vec_ini(igs)
      read(5,*)vec_ini(igs)
     enddo 
    endif 
    !--
    !B_MAT
    allocate(B_MAT(nigs,n_occ));B_MAT(:,:)=0.0D0 
    !default flg_BMAT==0
    if(flg_BMAT==0)then  
     write(6,'(a)')'BMAT is 0'
     if(nigs/=n_occ)then 
      write(6,'(a)')'wrong input: nigs/=n_occ'
      stop
     endif 
     do igs=1,nigs!=n_occ 
      B_MAT(igs,igs)=1.0d0 
     enddo 
    elseif(flg_BMAT==1)then 
     write(6,'(a)')'BMAT from input.in'
     if(reading_bmat_format==0)then 
      write(6,'(a)')'READING BMAT format=0: respack'
      do igs=1,nigs 
       !read(999,*)(B_MAT(igs,ix),ix=1,n_occ)
        read(5,*)(B_MAT(igs,ix),ix=1,n_occ)
      enddo
     endif 
     if(reading_bmat_format==1)then 
      write(6,'(a)')'READING BMAT format=1: xtapp'
      do ix=1,n_occ 
       !read(999,*)(B_MAT(igs,ix),igs=1,nigs)
        read(5,*)(B_MAT(igs,ix),igs=1,nigs)
      enddo
     endif 
    endif 
    !--
    if(flg_vis_bloch==0) CALC_REAL_SPACE_BLOCH=.false.!default 
    if(flg_vis_bloch==1) CALC_REAL_SPACE_BLOCH=.true. 
    !--
    if(flg_fermisurface==0) CALC_FERMISURFACE=.false.!default 
    if(flg_fermisurface==1) CALC_FERMISURFACE=.true. 
    !--
    if(flg_global_dos==0) CALC_GLOBAL_DOS=.false.!default 
    if(flg_global_dos==1) CALC_GLOBAL_DOS=.true. 
    !--
    write(6,param_wannier) 
    do igs=1,nigs 
     !write(6,'(a3,x,4f8.4,9f8.4)') vec_ini(igs) 
     !
     !20190609 Kazuma Nakamura
     !
     write(6,'(a3,x,4f8.4,x,9f8.4,x,8f10.5)') vec_ini(igs) 
     !
    enddo 
    do igs=1,nigs 
     write(6,'(100f10.5)')(B_MAT(igs,ix),ix=1,n_occ)
    enddo
    write(6,*) 
    !--
    !&param_interpolation
    !--
    !default
    Ndiv=40!Separation between symmetry points  
    dense(1:3)=0
    !
    !20190616 Kazuma Nakamura
    !
    flg_3d_fatband=0 
    !--
    !read(999,nml=param_interpolation)
    read(5,nml=param_interpolation)
    !k for band dispersion
    allocate(SK_sym_pts(3,N_sym_points))  
    if(reading_sk_format==0)then 
     write(6,'(a)')'READING SK_sym_pts format=0: respack'
     do ik=1,N_sym_points 
      !read(999,*)(SK_sym_pts(ix,ik),ix=1,3) 
      read(5,*)(SK_sym_pts(ix,ik),ix=1,3) 
     enddo 
    endif 
    if(reading_sk_format==1)then 
     write(6,'(a)')'READING SK_sym_pts format=1: xtapp'
     do ix=1,3 
      !read(999,*)(SK_sym_pts(ix,ik),ik=1,N_sym_points) 
      read(5,*)(SK_sym_pts(ix,ik),ik=1,N_sym_points) 
     enddo 
    endif 
    !
    !20190616 Kazuma Nakamura
    !
    if(flg_3d_fatband==0) CALC_3D_FATBAND=.false.!default 
    if(flg_3d_fatband==1) CALC_3D_FATBAND=.true. 
    !--
    write(6,param_interpolation) 
    do ik=1,N_sym_points 
     write(6,'(3f10.5)')(SK_sym_pts(ix,ik),ix=1,3) 
    enddo 
    write(6,*) 
    !--
    !&param_visualization)
    !--
    !default
    ix_vis_min=-1
    ix_vis_max=1
    iy_vis_min=-1
    iy_vis_max=1
    iz_vis_min=-1
    iz_vis_max=1
    !--
    !read(999,nml=param_visualization)
    read(5,nml=param_visualization)
    !default
    if(flg_vis_wannier==0) CALC_REAL_SPACE_WANNIER=.false.!default 
    if(flg_vis_wannier==1) CALC_REAL_SPACE_WANNIER=.true. 
    !write vis wannier 
    if(N_write_wannier==0)then
     N_write_wannier=N_wannier 
     allocate(wrt_list(N_write_wannier)) 
     do ix=1,N_write_wannier 
      wrt_list(ix)=ix
     enddo 
    elseif(N_write_wannier/=0)then 
     allocate(wrt_list(N_write_wannier)) 
     !read(999,*)(wrt_list(ix),ix=1,N_write_wannier)
     read(5,*)(wrt_list(ix),ix=1,N_write_wannier)
    endif  
    xmin=ix_vis_min
    xmax=ix_vis_max
    ymin=iy_vis_min
    ymax=iy_vis_max
    zmin=iz_vis_min
    zmax=iz_vis_max
    !--
    write(6,param_visualization) 
    write(6,'(100i5)')(wrt_list(ix),ix=1,N_write_wannier) 
    write(6,*) 
    !--
    end subroutine read_input
    !
end module m_rdinput 
