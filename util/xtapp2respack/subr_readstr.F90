!
! Copyright (c) 2013,2016 Yoshihide Yoshimoto
!
module subr_readstr
contains
  subroutine read_strfile(nf,nkd,nsi,ax,aa,kd,asi,cforce,zo,zn, &
    o_total_energy,o_fermi_energy)
    implicit none
    integer,intent(in):: nf
    real(kind=8):: ax
    integer:: nkd,nsi
    real(kind=8):: aa(3,3)
    integer,pointer:: kd(:)
    real(kind=8),pointer:: asi(:,:),cforce(:,:)
    real(kind=8),pointer:: zo(:),zn(:)
    real(kind=8),optional:: o_total_energy,o_fermi_energy(2)

    real(kind=8):: lattice_factor, lattice_list(9), total_energy
    real(kind=8):: stress_tensor(6), fermi_energy(2)
    integer:: number_element, number_atom
    real(kind=8):: spin_polarization, abs_spin_polarization
    real(kind=8):: spin_polarization_vector(3)
    real(kind=8):: abs_spin_polarization_vector(3)
    integer:: ikd,it

    namelist /struct_data/ lattice_factor, lattice_list, total_energy, &
         stress_tensor, fermi_energy, number_element, number_atom, &
         spin_polarization, spin_polarization_vector, abs_spin_polarization, &
         abs_spin_polarization_vector

    read(nf, struct_data, err=100, end=100)

    if (present(o_total_energy)) then
       o_total_energy = total_energy
    end if
    if (present(o_fermi_energy)) then
       o_fermi_energy = fermi_energy
    end if

    ax = lattice_factor
    aa(:,1) = lattice_list(1:3)
    aa(:,2) = lattice_list(4:6)
    aa(:,3) = lattice_list(7:9)
    nkd = number_element
    nsi = number_atom
    allocate(zo(nkd),zn(nkd))
    allocate(kd(nsi),asi(3,nsi),cforce(3,nsi))

    read(nf,*,err=200)
    read(nf,*,err=200)
    do ikd=1,nkd
       read(nf,*,err=200) zo(ikd),zn(ikd)
    end do

    read(nf,*,err=300)
    read(nf,*,err=300)
    do it=1,nsi
       read(nf,*,err=300) kd(it),asi(1,it),asi(2,it),asi(3,it)
    end do

    read(nf,*,err=400)
    read(nf,*,err=400)
    do it=1,nsi
       read(nf,*,err=400) cforce(1,it),cforce(2,it),cforce(3,it)
    end do

    return

100 continue
    write(6,'(a)')'error in reading struct_data'
    stop
200 continue
    write(6,'(a)')'error in reading valence and nucleus charge'
    stop
300 continue
    write(6,'(a)')'error in reading atom kind and position'
    stop
400 continue
    write(6,'(a)')'error in reading force'
    stop
  end subroutine read_strfile

  subroutine xread_strfile(nf,nkd,nsi,ax,aa,kd,asi,cforce,zo,zn, &
       etotal,stress,ef,spinpol,aspinpol)
    implicit none
    integer,intent(in):: nf
    real(kind=8):: ax
    integer:: nkd,nsi
    real(kind=8):: aa(3,3)
    integer,pointer:: kd(:)
    real(kind=8),pointer:: asi(:,:),cforce(:,:)
    real(kind=8),pointer:: zo(:),zn(:)
    real(kind=8):: etotal,stress(6),ef(2),spinpol,aspinpol

    real(kind=8):: lattice_factor, lattice_list(9), total_energy
    real(kind=8):: stress_tensor(6), fermi_energy(2)
    integer:: number_element, number_atom
    real(kind=8):: spin_polarization, abs_spin_polarization
    real(kind=8):: spin_polarization_vector(3)
    real(kind=8):: abs_spin_polarization_vector(3)
    integer:: ikd,it

    namelist /struct_data/ lattice_factor, lattice_list, total_energy, &
         stress_tensor, fermi_energy, number_element, number_atom, &
         spin_polarization, abs_spin_polarization, &
         spin_polarization_vector, abs_spin_polarization_vector

    read(nf, struct_data, err=100, end=100)

    ax = lattice_factor
    aa(:,1) = lattice_list(1:3)
    aa(:,2) = lattice_list(4:6)
    aa(:,3) = lattice_list(7:9)
    nkd = number_element
    nsi = number_atom
    etotal = total_energy
    stress = stress_tensor
    ef = fermi_energy
    spinpol = spin_polarization
    aspinpol = abs_spin_polarization

    allocate(zo(nkd),zn(nkd))
    allocate(kd(nsi),asi(3,nsi),cforce(3,nsi))

    read(nf,*,err=200)
    read(nf,*,err=200)
    do ikd=1,nkd
       read(nf,*,err=200) zo(ikd),zn(ikd)
    end do

    read(nf,*,err=300)
    read(nf,*,err=300)
    do it=1,nsi
       read(nf,*,err=300) kd(it),asi(1,it),asi(2,it),asi(3,it)
    end do

    read(nf,*,err=400)
    read(nf,*,err=400)
    do it=1,nsi
       read(nf,*,err=400) cforce(1,it),cforce(2,it),cforce(3,it)
    end do

    return

100 continue
    write(6,'(a)')'error in reading struct_data'
    stop
200 continue
    write(6,'(a)')'error in reading valence and nucleus charge'
    stop
300 continue
    write(6,'(a)')'error in reading atom kind and position'
    stop
400 continue
    write(6,'(a)')'error in reading force'
    stop
  end subroutine xread_strfile

end module subr_readstr
