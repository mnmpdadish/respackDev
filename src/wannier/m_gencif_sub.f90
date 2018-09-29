!
! Copyright (c) 2013-2018 Yoshihide Yoshimoto and Kazuma Nakamura 
!
module m_gencif_sub!subr_fmtconv: original name in xtapp
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
  real(8),parameter::oneau=0.52917721067d0!by 2014 CODATA
  real(8),parameter::pi=3.14159265358979323846d0
contains
  subroutine printcif(cbuf,aa,nkd,zn,nsi,kd,asi,na,nb,nc) 
    implicit none
    character(len=*)::cbuf
    real(8),intent(in)::aa(3,3)
    integer,intent(in)::nkd,nsi,kd(nsi)
    real(8)::zn(nkd)
    real(8)::asi(3,nsi)
    integer::nsikd(nkd)
    character(len=10240)::chem_abs_conf
    character(len=6)::cnsikd
    integer::it,ikd,izn
    real(8)::cla,clb,clc,tmp
    !
    character(99)::filename 
    integer,parameter::file_num=304
    integer::na,nb,nc 
    !
    !OUTPUT by cif format 
    !OPEN(file_num,W,FILE='dat.supercell-xxx.cif') 
    write(filename,"('dat.supercell-',i3.3,'x',i3.3,'x',i3.3,'.cif')")na,nb,nc
    write(6,*)'filename',filename 
    OPEN(file_num,FILE=filename) 
    REWIND(file_num)
    ! 
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
    !
    write(file_num,'(a)') 'data_global'
    write(file_num,'(a)') '_chemical_name ''' // trim(cbuf) // ''''
    write(file_num,'(a)') '_chemical_absolute_configuration ' // trim(chem_abs_conf)
    write(file_num,'(a)') '_symmetry_space_group_name ''P1'''
    ! 
    cla = aa(1,1)*aa(1,1)+aa(2,1)*aa(2,1)+aa(3,1)*aa(3,1)
    cla = sqrt(cla)
    write(file_num,'(a,f15.7)') '_cell_length_a ', oneau*cla
    !
    clb = aa(1,2)*aa(1,2)+aa(2,2)*aa(2,2)+aa(3,2)*aa(3,2)
    clb = sqrt(clb)
    write(file_num,'(a,f15.7)') '_cell_length_b ', oneau*clb
    !
    clc = aa(1,3)*aa(1,3)+aa(2,3)*aa(2,3)+aa(3,3)*aa(3,3)
    clc = sqrt(clc)
    write(file_num,'(a,f15.7)') '_cell_length_c ', oneau*clc
    !
    tmp = aa(1,2)*aa(1,3)+aa(2,2)*aa(2,3)+aa(3,2)*aa(3,3)
    tmp = tmp/(clb*clc)
    tmp = acos(tmp)/pi*180.0
    write(file_num,'(a,f11.7)') '_cell_angle_alpha ', tmp
    !
    tmp = aa(1,3)*aa(1,1)+aa(2,3)*aa(2,1)+aa(3,3)*aa(3,1)
    tmp = tmp/(clc*cla)
    tmp = acos(tmp)/pi*180.0
    write(file_num,'(a,f11.7)') '_cell_angle_beta ', tmp
    !
    tmp = aa(1,1)*aa(1,2)+aa(2,1)*aa(2,2)+aa(3,1)*aa(3,2)
    tmp = tmp/(cla*clb)
    tmp = acos(tmp)/pi*180.0
    write(file_num,'(a,f11.7)') '_cell_angle_gamma ', tmp
    !
    write(file_num,'(a,f15.7)') '_atom_sites_Cartn_tran_matrix_11 ', oneau*aa(1,1)
    write(file_num,'(a,f15.7)') '_atom_sites_Cartn_tran_matrix_21 ', oneau*aa(2,1)
    write(file_num,'(a,f15.7)') '_atom_sites_Cartn_tran_matrix_31 ', oneau*aa(3,1)
    write(file_num,'(a,f15.7)') '_atom_sites_Cartn_tran_matrix_12 ', oneau*aa(1,2)
    write(file_num,'(a,f15.7)') '_atom_sites_Cartn_tran_matrix_22 ', oneau*aa(2,2)
    write(file_num,'(a,f15.7)') '_atom_sites_Cartn_tran_matrix_32 ', oneau*aa(3,2)
    write(file_num,'(a,f15.7)') '_atom_sites_Cartn_tran_matrix_13 ', oneau*aa(1,3)
    write(file_num,'(a,f15.7)') '_atom_sites_Cartn_tran_matrix_23 ', oneau*aa(2,3)
    write(file_num,'(a,f15.7)') '_atom_sites_Cartn_tran_matrix_33 ', oneau*aa(3,3)
    write(file_num,*)
    !
    write(file_num,'(a)') 'loop_'
    write(file_num,'(a)') '_symmetry_equiv_pos_as_xyz'
    write(file_num,'(a)') '''x, y, z'''
    write(file_num,*)
    !
    write(file_num,'(a)') 'loop_'
    write(file_num,'(a)') '_atom_site_type_symbol'
    write(file_num,'(a)') '_atom_site_label'
    write(file_num,'(a)') '_atom_site_fract_x'
    write(file_num,'(a)') '_atom_site_fract_y'
    write(file_num,'(a)') '_atom_site_fract_z'
    !
    nsikd = 0
    do it=1,nsi
       izn = nint(zn(kd(it)))
       nsikd(kd(it)) = nsikd(kd(it))+1
       write(cnsikd,'(i0)') nsikd(kd(it))
       write(file_num,'(a,1x,a,3(1x,f13.9))') trim(atom_name(izn)), &
            trim(atom_name(izn))//trim(cnsikd), &
            modulo(asi(1,it),1.0d0), &
            modulo(asi(2,it),1.0d0), &
            modulo(asi(3,it),1.0d0)
    end do
    !
    close(file_num)
    !
    return
  end subroutine printcif
end module m_gencif_sub
