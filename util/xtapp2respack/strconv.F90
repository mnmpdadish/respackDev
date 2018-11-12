!
! Copyright (c) 2013-2015,2018 Yoshihide Yoshimoto
!
program strconv
  use subr_fmtconv
  use subr_readstr
  implicit none
  real(kind=8):: ax
  integer:: nkd,nsi
  integer,parameter:: nfin = 10, nfout = 11
  !
  real(kind=8):: aa(3,3)
  integer,pointer:: kd(:)
  real(kind=8),pointer:: asi(:,:),cforce(:,:)

  real(kind=8),pointer:: zo(:),zn(:)
  real(kind=8):: etotal,efermi(2)

  character(len=256):: cbuf
  integer:: ftype

  integer:: tatm
  real(kind=8):: dspl(3)

  integer:: ios

  call getarg(1, cbuf)

  if (trim(cbuf).eq.'cif') then
     ftype = 1
  else if (trim(cbuf).eq.'xyz') then
     ftype = 2
  else if (trim(cbuf).eq.'xyzpr') then
     ftype = 3
  else if (trim(cbuf).eq.'xyzmol') then
     ftype = 4
  else if (trim(cbuf).eq.'poscar') then
     ftype = 5
  else if (trim(cbuf).eq.'cforce') then
     ftype = 6
  else if (trim(cbuf).eq.'respack') then
     ftype = 7
  else
     write(6,*) 'usage: strconv (cif|xyz|xyzpr|xyzmol|poscar|cforce|respack) fname'
     stop
  end if

  call getarg(2, cbuf)
  open(nfin, file=cbuf, status='old', action='read', iostat=ios)
  if (ios.ne.0) then
     write(6,*) 'error in opening '//trim(cbuf)
     stop
  end if

  call read_strfile(nfin,nkd,nsi,ax,aa,kd,asi,cforce,zo,zn,etotal,efermi)
  close(nfin)

  aa = ax*aa

  cbuf = 'xTAPP_str: ' // cbuf

  if (ftype.eq.1) then
     call printcif(cbuf,aa,nkd,zo,zn,nsi,kd,asi)
  else if (ftype.eq.2) then
     call printxyz(cbuf,aa,nkd,zo,zn,nsi,kd,asi)
  else if (ftype.eq.3) then
     call printxyzpr(cbuf,aa,nkd,zo,zn,nsi,kd,asi)
  else if (ftype.eq.4) then
     call printxyzmol(cbuf,aa,nkd,zo,zn,nsi,kd,asi)
  else if (ftype.eq.5) then
     call printposcar(6,ax,aa,nkd,zo,zn,nsi,kd,asi)
  else if (ftype.eq.6) then
     call printcforce(6,nsi,cforce)
  else if (ftype.eq.7) then
     call printrespack(nfout,cbuf,aa,nkd,zo,zn,nsi,kd,asi,etotal,efermi)
  end if

  stop

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
end program strconv
