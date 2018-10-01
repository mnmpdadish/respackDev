!
! Copyright (c) 2013-2015,2018 Yoshihide Yoshimoto
!
module subr_fmtconv
  character(len=2):: atom_name(0:102) = (/ 'Un', 'H ','He', &
       'Li','Be','B ','C ','N ','O ','F ','Ne', &
       'Na','Mg','Al','Si','P ','S ','Cl','Ar', &
       'K ','Ca', &
       'Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn', &
       'Ga','Ge','As','Se','Br','Kr', &
       'Rb','Sr', &
       'Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd', &
       'In','Sn','Sb','Te','I ','Xe', &
       'Cs','Ba', &
       'La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb', &
       'Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg', &
       'Tl','Pb','Bi','Po','At','Rn', &
       'Fr','Ra', &
       'Ac','Th','Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No' /)
  real(kind=8),parameter:: oneau = 0.52917721067d0 ! by 2014 CODATA
  real(kind=8),parameter:: pi = 3.14159265358979323846d0

contains
  subroutine printcif(cbuf,aa,nkd,zo,zn,nsi,kd,asi)
    implicit none
    character(len=*):: cbuf
    real(kind=8),intent(in):: aa(3,3)
    integer,intent(in):: nkd,nsi,kd(nsi)
    real(kind=8),pointer:: zo(:),zn(:)
    real(kind=8),pointer:: asi(:,:)

    integer:: nsikd(nkd)
    character(len=10240):: chem_abs_conf
    character(len=6):: cnsikd
    integer:: it,ikd,izn
    real(kind=8):: cla,clb,clc,tmp

    nsikd = 0
    do it=1,nsi
       nsikd(kd(it)) = nsikd(kd(it))+1
    end do
    chem_abs_conf = ''''
    do ikd=1,nkd
       izn = nint(zn(ikd))
       write(cnsikd,'(i0)') nsikd(ikd)
       chem_abs_conf = trim(chem_abs_conf) // ' ' &
            // trim(atom_name(izn)) // trim(cnsikd)
    end do
    chem_abs_conf = trim(chem_abs_conf)//''''

    write(6,'(a)') 'data_global'
    write(6,'(a)') '_chemical_name ''' // trim(cbuf) // ''''
    write(6,'(a)') '_chemical_absolute_configuration ' // trim(chem_abs_conf)
    write(6,'(a)') '_symmetry_space_group_name ''P1'''
    
    cla = aa(1,1)*aa(1,1)+aa(2,1)*aa(2,1)+aa(3,1)*aa(3,1)
    cla = sqrt(cla)
    write(6,'(a,f15.7)') '_cell_length_a ', oneau*cla

    clb = aa(1,2)*aa(1,2)+aa(2,2)*aa(2,2)+aa(3,2)*aa(3,2)
    clb = sqrt(clb)
    write(6,'(a,f15.7)') '_cell_length_b ', oneau*clb

    clc = aa(1,3)*aa(1,3)+aa(2,3)*aa(2,3)+aa(3,3)*aa(3,3)
    clc = sqrt(clc)
    write(6,'(a,f15.7)') '_cell_length_c ', oneau*clc

    tmp = aa(1,2)*aa(1,3)+aa(2,2)*aa(2,3)+aa(3,2)*aa(3,3)
    tmp = tmp/(clb*clc)
    tmp = acos(tmp)/pi*180.0
    write(6,'(a,f11.7)') '_cell_angle_alpha ', tmp

    tmp = aa(1,3)*aa(1,1)+aa(2,3)*aa(2,1)+aa(3,3)*aa(3,1)
    tmp = tmp/(clc*cla)
    tmp = acos(tmp)/pi*180.0
    write(6,'(a,f11.7)') '_cell_angle_beta ', tmp

    tmp = aa(1,1)*aa(1,2)+aa(2,1)*aa(2,2)+aa(3,1)*aa(3,2)
    tmp = tmp/(cla*clb)
    tmp = acos(tmp)/pi*180.0
    write(6,'(a,f11.7)') '_cell_angle_gamma ', tmp

    write(6,'(a,f15.7)') '_atom_sites_Cartn_tran_matrix_11 ', oneau*aa(1,1)
    write(6,'(a,f15.7)') '_atom_sites_Cartn_tran_matrix_21 ', oneau*aa(2,1)
    write(6,'(a,f15.7)') '_atom_sites_Cartn_tran_matrix_31 ', oneau*aa(3,1)
    write(6,'(a,f15.7)') '_atom_sites_Cartn_tran_matrix_12 ', oneau*aa(1,2)
    write(6,'(a,f15.7)') '_atom_sites_Cartn_tran_matrix_22 ', oneau*aa(2,2)
    write(6,'(a,f15.7)') '_atom_sites_Cartn_tran_matrix_32 ', oneau*aa(3,2)
    write(6,'(a,f15.7)') '_atom_sites_Cartn_tran_matrix_13 ', oneau*aa(1,3)
    write(6,'(a,f15.7)') '_atom_sites_Cartn_tran_matrix_23 ', oneau*aa(2,3)
    write(6,'(a,f15.7)') '_atom_sites_Cartn_tran_matrix_33 ', oneau*aa(3,3)
    write(6,*)

    write(6,'(a)') 'loop_'
    write(6,'(a)') '_symmetry_equiv_pos_as_xyz'
    write(6,'(a)') '''x, y, z'''
    write(6,*)

    write(6,'(a)') 'loop_'
    write(6,'(a)') '_atom_site_type_symbol'
    write(6,'(a)') '_atom_site_label'
    write(6,'(a)') '_atom_site_fract_x'
    write(6,'(a)') '_atom_site_fract_y'
    write(6,'(a)') '_atom_site_fract_z'

    nsikd = 0
    do it=1,nsi
       izn = nint(zn(kd(it)))
       nsikd(kd(it)) = nsikd(kd(it))+1
       write(cnsikd,'(i0)') nsikd(kd(it))
       write(6,'(a,1x,a,3(1x,f13.9))') trim(atom_name(izn)), &
            trim(atom_name(izn))//trim(cnsikd), &
            modulo(asi(1,it),1.0d0), &
            modulo(asi(2,it),1.0d0), &
            modulo(asi(3,it),1.0d0)
    end do
    return
  end subroutine printcif

  subroutine printxyz(cbuf,aa,nkd,zo,zn,nsi,kd,asi)
    implicit none
    character(len=*):: cbuf
    real(kind=8),intent(in):: aa(3,3)
    integer,intent(in):: nkd,nsi,kd(nsi)
    real(kind=8),pointer:: zo(:),zn(:)
    real(kind=8),pointer:: asi(:,:)

    integer:: it,izn
    real(kind=8):: x,y,z

    write(6,'(i0)') nsi
    write(6,'(a)') trim(cbuf)

    do it=1,nsi
       izn = nint(zn(kd(it)))
       x = aa(1,1)*asi(1,it)+aa(1,2)*asi(2,it)+aa(1,3)*asi(3,it)
       y = aa(2,1)*asi(1,it)+aa(2,2)*asi(2,it)+aa(2,3)*asi(3,it)
       z = aa(3,1)*asi(1,it)+aa(3,2)*asi(2,it)+aa(3,3)*asi(3,it)
       write(6,'(a,3(1x,e15.7))') trim(atom_name(izn)), &
            oneau*x, oneau*y, oneau*z
    end do
    return
  end subroutine printxyz

  subroutine printxyzpr(cbuf,aa,nkd,zo,zn,nsi,kd,asi)
    implicit none
    character(len=*):: cbuf
    real(kind=8),intent(in):: aa(3,3)
    integer,intent(in):: nkd,nsi,kd(nsi)
    real(kind=8),pointer:: zo(:),zn(:)
    real(kind=8),pointer:: asi(:,:)

    integer:: it,izn,i,i1,i2,i3,nsipr
    real(kind=8):: x,y,z,a,b,c,delta(3)
    real(kind=8),parameter:: eps = 0.1d0

    do i=1,3
       delta(i) = eps/sqrt( aa(1,i)*aa(1,i)+aa(2,i)*aa(2,i)+aa(3,i)*aa(3,i) )
    end do

    nsipr = 0
    do it=1,nsi
       izn = nint(zn(kd(it)))
       do i3=-1,1
          do i2=-1,1
             do i1=-1,1
                a = modulo(asi(1,it),1.0d0) + i1
                b = modulo(asi(2,it),1.0d0) + i2
                c = modulo(asi(3,it),1.0d0) + i3
                if ( -delta(1).lt.a .and. a.lt.1.0d0+delta(1) .and. &
                     -delta(2).lt.b .and. b.lt.1.0d0+delta(2) .and. &
                     -delta(3).lt.c .and. c.lt.1.0d0+delta(3) ) then
                   nsipr = nsipr + 1
                end if
             end do
          end do
       end do
    end do

    write(6,'(i0)') nsipr
    write(6,'(a)') trim(cbuf)

    i = 0
    do it=1,nsi
       izn = nint(zn(kd(it)))
       do i3=-1,1
          do i2=-1,1
             do i1=-1,1
                a = modulo(asi(1,it),1.0d0) + i1
                b = modulo(asi(2,it),1.0d0) + i2
                c = modulo(asi(3,it),1.0d0) + i3
                if ( -delta(1).lt.a .and. a.lt.1.0d0+delta(1) .and. &
                     -delta(2).lt.b .and. b.lt.1.0d0+delta(2) .and. &
                     -delta(3).lt.c .and. c.lt.1.0d0+delta(3) ) then
                   x = aa(1,1)*a + aa(1,2)*b + aa(1,3)*c
                   y = aa(2,1)*a + aa(2,2)*b + aa(2,3)*c
                   z = aa(3,1)*a + aa(3,2)*b + aa(3,3)*c
                   write(6,'(a,3(1x,e15.7))') trim(atom_name(izn)), &
                        oneau*x, oneau*y, oneau*z
                   i = i + 1 
                end if
             end do
          end do
       end do
    end do

    if (i.ne.nsipr) then
       write(6,*) 'error in printxyzpr()'
       stop
    end if
    return
  end subroutine printxyzpr

  subroutine printxyzmol(cbuf,aa,nkd,zo,zn,nsi,kd,asi)
    implicit none
    character(len=*):: cbuf
    real(kind=8),intent(in):: aa(3,3)
    integer,intent(in):: nkd,nsi,kd(nsi)
    real(kind=8),pointer:: zo(:),zn(:)
    real(kind=8),pointer:: asi(:,:)

    integer:: it,izn,i1,i2,i3
    real(kind=8):: x,y,z,c(3),as(3),am(3),sc(3),d,td

    write(6,'(i0)') nsi
    write(6,'(a)') trim(cbuf)

    it = 1
    sc = asi(:,it)
    izn = nint(zn(kd(it)))
    x = aa(1,1)*asi(1,it)+aa(1,2)*asi(2,it)+aa(1,3)*asi(3,it)
    y = aa(2,1)*asi(1,it)+aa(2,2)*asi(2,it)+aa(2,3)*asi(3,it)
    z = aa(3,1)*asi(1,it)+aa(3,2)*asi(2,it)+aa(3,3)*asi(3,it)
    write(6,'(a,3(1x,e15.7))') trim(atom_name(izn)), &
         oneau*x, oneau*y, oneau*z

    do it=2,nsi
       izn = nint(zn(kd(it)))
       c = sc/(it-1)
       d = huge(1.0d0)
       am = huge(1.0d0)
       do i3=-1,1
          do i2=-1,1
             do i1=-1,1
                as(1) = modulo(asi(1,it),1.0d0) + i1
                as(2) = modulo(asi(2,it),1.0d0) + i2
                as(3) = modulo(asi(3,it),1.0d0) + i3
                td = abs(as(1)-c(1))+abs(as(2)-c(2))+abs(as(3)-c(3))
                if (td.lt.d) then
                   d = td
                   am = as
                end if
             end do
          end do
       end do
       sc = sc + am

       x = aa(1,1)*am(1)+aa(1,2)*am(2)+aa(1,3)*am(3)
       y = aa(2,1)*am(1)+aa(2,2)*am(2)+aa(2,3)*am(3)
       z = aa(3,1)*am(1)+aa(3,2)*am(2)+aa(3,3)*am(3)
       write(6,'(a,3(1x,e15.7))') trim(atom_name(izn)), &
            oneau*x, oneau*y, oneau*z
    end do
    return
  end subroutine printxyzmol

  subroutine printposcar(nf,ax,aa,nkd,zo,zn,nsi,kd,asi)
    implicit none
    integer,intent(in):: nf
    real(kind=8),intent(in):: ax,aa(3,3)
    integer,intent(in):: nkd,nsi,kd(nsi)
    real(kind=8),pointer:: zo(:),zn(:)
    real(kind=8),pointer:: asi(:,:)

    integer:: it,ikd
    integer:: nsikd(nkd)

    do ikd=1,nkd
       write(nf,'(1x,a)',advance='no')  atom_name(nint(zn(ikd)))
    end do
    write(nf,*)
    write(nf,'(1x,e18.10)') ax*oneau
    write(nf,'(3(1x,e18.10))') aa(:,1)/ax
    write(nf,'(3(1x,e18.10))') aa(:,2)/ax
    write(nf,'(3(1x,e18.10))') aa(:,3)/ax

    nsikd = 0
    do it=1,nsi
       nsikd(kd(it)) = nsikd(kd(it)) + 1
    end do
    do ikd=1,nkd
       write(nf,'(1x,i6)',advance='no') nsikd(ikd)
    end do
    write(nf,*)

    write(nf,'(a)') 'Direct'
    do ikd=1,nkd
       do it=1,nsi
          if (kd(it).eq.ikd) then
             write(nf,'(1x,3f16.12)') &
                  asi(1,it),asi(2,it),asi(3,it)
          end if
       end do
    end do
    return
  end subroutine printposcar

  subroutine printrespack(nf,cbuf,aa,nkd,zo,zn,nsi,kd,asi,etotal,efermi)
    implicit none
    integer,intent(in):: nf
    character(len=*):: cbuf
    real(kind=8),intent(in):: aa(3,3)
    integer,intent(in):: nkd,nsi,kd(nsi)
    real(kind=8),pointer:: zo(:),zn(:)
    real(kind=8),pointer:: asi(:,:)
    real(kind=8),intent(in):: etotal,efermi(2)

    integer:: it,izn,ios

    open(nf, file='dat.atom_position', status='replace', &
         action='write', iostat=ios)
    if (ios.ne.0) then
       write(6,*) 'error in opening '//trim(cbuf)
       stop
    end if
    
    write(nf,'(1x,i0)') nsi
    do it=1,nsi
       izn = nint(zn(kd(it)))
       write(nf,'(1x,a,1x,3(1x,f13.9))') trim(atom_name(izn)), &
            modulo(asi(1,it),1.0d0), &
            modulo(asi(2,it),1.0d0), &
            modulo(asi(3,it),1.0d0)
    end do
    close(nf)

    open(nf, file='dat.bandcalc', status='old', &
         action='write', position='append', iostat=ios)
    if (ios.ne.0) then
       write(6,*) 'error in opening '//trim(cbuf)
       stop
    end if
    write(nf,'(1x,e19.10,1x,e19.10)') efermi(1),efermi(2)
    write(nf,'(1x,e21.12)') etotal
    close(nf)
    
    return
  end subroutine printrespack

  subroutine printcforce(nf,nsi,cforce)
    implicit none
    integer,intent(in):: nf
    integer,intent(in):: nsi
    real(kind=8),pointer:: cforce(:,:)

    integer:: it

    do it=1,nsi
       write(nf,'(3(1x,e18.10))') cforce(1,it),cforce(2,it),cforce(3,it)
    end do
    return
  end subroutine printcforce
end module subr_fmtconv
