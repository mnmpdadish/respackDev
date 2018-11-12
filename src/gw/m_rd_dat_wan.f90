module m_rd_dat_wan
use m_rd_dat_wfn
implicit none
private 
public::rd_dat_ns_nb 
public::rd_dat_umat 
public::rd_dat_wan 
public::rd_dat_hmatr 
!ns-nb(149)  
integer,public::Mb,Mt  
integer,public,allocatable::Ns(:)!Ns(NTK) 
integer,public,allocatable::Nb(:)!Nb(NTK) 
integer,public,allocatable::Nt(:)!Nt(NTK) 
!umat(150)  
integer,public::NWF 
complex(8),public,allocatable::UNT(:,:,:)!UNT(Mb,NWF,NTK) 
!c0_wn(113)  
complex(8),public,allocatable::C0_WN(:,:,:)!C0_WN(NTG,NWF,NTK) 
!h_mat_r(122) 
complex(8),public,allocatable::H_MAT_R(:,:,:,:,:)!H_MAT_R(NWF,NWF,-Na1:Na1,-Na2:Na2,-Na3:Na3) 
contains
subroutine rd_dat_ns_nb 
implicit none 
integer::ik 
!--
OPEN(149,FILE='./dir-wan/dat.ns-nb') 
rewind(149) 
allocate(Ns(NTK));Ns=0
allocate(Nb(NTK));Nb=0
allocate(Nt(NTK));Nt=0
do ik=1,NTK 
 read(149,*) Ns(ik),Nb(ik) 
enddo  
Mb=maxval(Nb) 
!write(6,*)'Mb=',Mb 
Nt(:)=0
do ik=1,NTK 
 Nt(ik)=Ns(ik)+Nb(ik) 
enddo 
Mt=maxval(Nt) 
!write(6,*) 'Mt=',Mt 
!---
end subroutine
!--
subroutine rd_dat_umat  
implicit none 
integer::ik,jb,jw
OPEN(150,FILE='./dir-wan/dat.umat') 
rewind(150) 
read(150,*) NWF 
allocate(UNT(Mb,NWF,NTK));UNT(:,:,:)=0.0d0 
do ik=1,NTK
 do jb=1,Nb(ik)
  read(150,*)(UNT(jb,jw,ik),jw=1,NWF) 
 enddo 
enddo 
!write(6,*)'NWF=',NWF 
end subroutine
!--
subroutine rd_dat_wan 
implicit none
integer::ik,iw,ig,nwftmp 
OPEN(113,FILE='./dir-wan/dat.wan',FORM='unformatted') 
REWIND(113)       
read(113) nwftmp 
if(nwftmp/=NWF)then
 write(6,*)'ERROR; STOP; nwftmp should be NWF'   
 write(6,*)'nwftmp=',nwftmp,'NWF=',NWF;STOP
endif 
allocate(C0_WN(NTG,NWF,NTK));C0_WN(:,:,:)=0.0D0 
do ik=1,NTK
 read(113)((C0_WN(ig,iw,ik),ig=1,NG0(ik)),iw=1,NWF)           
 !
 !20180519 
 !do iw=1,NWF
 ! read(113)(C0_WN(ig,iw,ik),ig=1,NG0(ik))
 !enddo!iw
 !
enddo!ik 
CLOSE(113)  
!write(6,*)'FINISH REDING C0_WN'
!--
!do ik=1,Nk_irr 
! do iw=1,NWF
!  do jw=1,NWF 
!   SUM_CMPX=0.0d0 
!   do ig=1,NG0(ik) 
!    SUM_CMPX=SUM_CMPX+CONJG(C0_WN(ig,iw,ik))*C0_WN(ig,jw,ik) 
!   enddo 
!   write(6,'(3i5,x,2f15.10)') ik,iw,jw,SUM_CMPX 
!  enddo 
! enddo 
!enddo 
end subroutine
!--
subroutine rd_dat_hmatr 
implicit none
character(99)::header1,header2,header3 
integer::ia1,ia2,ia3,idum1,idum2,idum3,ib,jb 
real(8),parameter::au=27.21151d0
OPEN(122,FILE='./dir-wan/dat.h_mat_r') 
allocate(H_MAT_R(NWF,NWF,-Na1:Na1,-Na2:Na2,-Na3:Na3));H_MAT_R(:,:,:,:,:)=0.0d0  
REWIND(122)
read(122,'(a)') header1 
read(122,'(a)') header2 
read(122,'(a)') header3 
do ia1=-Na1,Na1!-1
 do ia2=-Na2,Na2!-1
  do ia3=-Na3,Na3!-1
   read(122,*) idum1,idum2,idum3 
   do ib=1,NWF 
    do jb=1,NWF 
     read(122,'(i5,i5,2f20.10)') idum1,idum2,H_MAT_R(ib,jb,ia1,ia2,ia3) 
    enddo!jb 
   enddo!ib 
   read(122,*) 
  enddo 
 enddo 
enddo 
CLOSE(122) 
H_MAT_R=H_MAT_R/au 
end subroutine
!--
end module 
