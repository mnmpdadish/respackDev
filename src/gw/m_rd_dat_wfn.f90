module m_rd_dat_wfn
implicit none
public::rd_dat_symmetry
public::rd_dat_bandcalc 
public::rd_dat_lattice 
public::rd_dat_sample_k 
public::rd_dat_nkm 
public::rd_dat_kg 
public::rd_dat_eigenvalue 
public::rd_dat_wavefunction 
!sym(100)
integer,public::nsymq,nnp
integer,public,allocatable::rg(:,:,:)!3,3,nsymq)
integer,public,allocatable::pg(:,:)!pg(3,nsymq)
real(8),public,allocatable::rginv(:,:,:)!rginv(3,3,nsymq) 
!bandcalc(117) 
real(8),public::Ecut_for_psi 
real(8),public::FermiEnergy  
real(8),public::Etot
!avec(105)
real(8),public::a1(3),a2(3),a3(3)
real(8),public::b1(3),b2(3),b3(3)
real(8),public::VOLUME,a,b,c,alp,bet,gmm  
!sample-k(101)
integer,public::Nk_irr,NTK 
integer,public::nkb1!SAMPLING K POINTS ALONG b1 VECTOR
integer,public::nkb2!SAMPLING K POINTS ALONG b2 VECTOR
integer,public::nkb3!SAMPLING K POINTS ALONG b3 VECTOR
integer,public::Na1!LATTICE TRANSLATIONL IN a1 20170327 
integer,public::Na2!LATTICE TRANSLATIONL IN a2 20170327 
integer,public::Na3!LATTICE TRANSLATIONL IN a3 20170327 
real(8),public,allocatable::SKI(:,:)!SKI(3,Nk_irr)  
real(8),public,allocatable::SK0(:,:)!SK0(3,NTK)  
integer,public,allocatable::numirr(:)!numirr(NTK) 
integer,public,allocatable::numrot(:)!numrot(NTK) 
integer,public,allocatable::trs(:)!trs(NTK) 
integer,public,allocatable::RW(:,:)!RW(3,NTK)
integer,public,allocatable::numMK(:)!numMK(Nk_irr)!20180316  
!eigenvalue(111) 
integer,public::NTB 
real(8),public,allocatable::E_EIGI(:,:)!E_EIGI(NTB,Nk_irr) 
!nkm(132)  
integer,public::NTG 
integer,public,allocatable::NGI(:)!NG0(Nk_irr)
!kg(104)
integer,public::L1,L2,L3 
integer,public,allocatable::KGI(:,:,:)!KGI(3,NTG,Nk_irr) 
integer,public,allocatable::KG0(:,:,:)!KG0(3,NTG,NTK)
integer,public,allocatable::NG0(:)!NG0(NTK)
integer,public,allocatable::packing(:,:,:,:)!packing(-L1:L1,-L2:L2,-L3:L3,Nk_irr)
!fft 
integer,public::nwx2,nwy2,nwz2!,nfft1,nfft2,nfft3,Nl123 
!wfn(102)
integer,public::ncomp 
!complex(8),public,allocatable::CIR(:,:,:)!CIR(NTG,NTB,Nk_irr) 
complex(4),public,allocatable::CIR(:,:,:)!CIR(NTG,NTB,Nk_irr) 
complex(8),allocatable::CIRtmp(:)!CIRtmp(NTG)  
contains
subroutine rd_dat_symmetry 
implicit none 
integer::i,j,iop
OPEN(100,FILE='./dir-wfn/dat.symmetry') 
rewind(100) 
read(100,*)nsymq 
read(100,*)nnp 
allocate(rg(3,3,nsymq));rg=0
allocate(pg(3,nsymq));pg=0
allocate(rginv(3,3,nsymq));rginv=0.0d0 
do iop=1,nsymq
 read(100,*)((rg(i,j,iop),i=1,3),j=1,3) 
 read(100,*)(pg(i,iop),i=1,3)   
enddo 
close(100) 
!
rginv=rg 
do iop=1,nsymq
 call invmat(3,rginv(1,1,iop)) 
enddo 
!
do iop=1,nsymq
 write(6,*) iop
 do i=1,3
  write(6,'(3I5,1x,3F15.10)')(rg(i,j,iop),j=1,3),(rginv(i,j,iop),j=1,3)
 enddo 
enddo 
!
end subroutine
!
subroutine rd_dat_bandcalc 
implicit none 
OPEN(117,FILE='./dir-wfn/dat.bandcalc') 
rewind(117) 
read(117,*)Ecut_for_psi 
read(117,*)FermiEnergy  
read(117,*)Etot
!write(6,*)'Ecut_for_psi=',Ecut_for_psi 
!write(6,*)'FermiEnergy=',FermiEnergy  
!write(6,*)'Etot=',Etot
end subroutine
! 
subroutine rd_dat_lattice 
implicit none 
real(8),parameter::tpi=2.0d0*dacos(-1.0d0)
OPEN(105,FILE='./dir-wfn/dat.lattice') 
REWIND(105)
READ(105,*) a1(1),a1(2),a1(3)!a1 vector
READ(105,*) a2(1),a2(2),a2(3)!a2 vector
READ(105,*) a3(1),a3(2),a3(3)!a3 vector
CLOSE(105)
!--
call OUTER_PRODUCT(a2(1),a3(1),b1(1))
VOLUME=a1(1)*b1(1)+a1(2)*b1(2)+a1(3)*b1(3)
b1(:)=b1(:)*tpi/VOLUME 
call OUTER_PRODUCT(a3(1),a1(1),b2(1))
b2(:)=b2(:)*tpi/VOLUME 
call OUTER_PRODUCT(a1(1),a2(1),b3(1))
b3(:)=b3(:)*tpi/VOLUME 
call est_latparam(a1(1),a2(1),a3(1),a,b,c,alp,bet,gmm) 
end subroutine
!--
subroutine rd_dat_sample_k 
implicit none 
integer::i,ik,jk,iop,iik   
real(8)::ktmp(3) 
integer::RWtmp(3) 
integer::initial_flg!20180316 
OPEN(101,FILE='./dir-wfn/dat.sample-k') 
rewind(101) 
read(101,*)Nk_irr 
allocate(SKI(3,Nk_irr));SKI(:,:)=0.0D0 
do ik=1,Nk_irr 
 read(101,*)(SKI(i,ik),i=1,3) 
enddo 
close(101) 
!--
call est_NTK(Nk_irr,nsymq,SKI(1,1),rg(1,1,1),NTK)
write(6,*)'Estimated NTK=',NTK 
!--
!SK0,numirr,numrot,trs,RW
allocate(SK0(3,NTK));SK0(:,:)=0.0d0
allocate(numirr(NTK));numirr(:)=0
allocate(numrot(NTK));numrot(:)=0
allocate(trs(NTK));trs(:)=0
allocate(RW(3,NTK));RW(:,:)=0
allocate(numMK(Nk_irr));numMK=0!20180316 
jk=0
do ik=1,Nk_irr
 initial_flg=0!20180316 
 do iop=1,nsymq
!sym
  ktmp(:)=0.0d0; RWtmp(:)=0  
  ktmp(1)=dble(rg(1,1,iop))*SKI(1,ik)+dble(rg(1,2,iop))*SKI(2,ik)+dble(rg(1,3,iop))*SKI(3,ik)
  ktmp(2)=dble(rg(2,1,iop))*SKI(1,ik)+dble(rg(2,2,iop))*SKI(2,ik)+dble(rg(2,3,iop))*SKI(3,ik)
  ktmp(3)=dble(rg(3,1,iop))*SKI(1,ik)+dble(rg(3,2,iop))*SKI(2,ik)+dble(rg(3,3,iop))*SKI(3,ik)
  call kcheck(ktmp(1),RWtmp(1))!rewind check 
  do iik=1,jk
   if(abs(SK0(1,iik)-ktmp(1))<1.0d-4.and.abs(SK0(2,iik)-ktmp(2))<1.0d-4.and.abs(SK0(3,iik)-ktmp(3))<1.0d-4) goto 1000
  enddo!iik
  jk=jk+1
  SK0(:,jk)=ktmp(:)
  numirr(jk)=ik
  numrot(jk)=iop
  trs(jk)=1
  RW(:,jk)=RWtmp(:)
  !
  !20180316
  !
  if(initial_flg.eq.0)then
   numMK(ik)=jk
   initial_flg=1 
  endif 
  !
!time-reversal
1000 ktmp(:)=0.0d0;RWtmp(:)=0  
  ktmp(1)=dble(rg(1,1,iop))*SKI(1,ik)+dble(rg(1,2,iop))*SKI(2,ik)+dble(rg(1,3,iop))*SKI(3,ik)
  ktmp(2)=dble(rg(2,1,iop))*SKI(1,ik)+dble(rg(2,2,iop))*SKI(2,ik)+dble(rg(2,3,iop))*SKI(3,ik)
  ktmp(3)=dble(rg(3,1,iop))*SKI(1,ik)+dble(rg(3,2,iop))*SKI(2,ik)+dble(rg(3,3,iop))*SKI(3,ik)
  call kcheck_trs(ktmp(1),RWtmp(1))!rewind check modified 20170316  
  do iik=1,jk
   if(abs(SK0(1,iik)-(-ktmp(1)))<1.0d-4.and.abs(SK0(2,iik)-(-ktmp(2)))<1.0d-4.and.abs(SK0(3,iik)-(-ktmp(3)))<1.0d-4) goto 2000
  enddo!iik
  jk=jk+1
  SK0(:,jk)=-ktmp(:) 
  numirr(jk)=ik
  numrot(jk)=iop
  trs(jk)=-1
  RW(:,jk)=RWtmp(:) 
2000 enddo!iop  
enddo!ik
call est_nkbi(NTK,SK0(1,1),nkb1,nkb2,nkb3)  
Na1=nkb1/2; Na2=nkb2/2; Na3=nkb3/2
!--
if(NTK/=jk) then 
 write(6,*)'ERROR;STOP;NTK should be jk'   
 write(6,*)'NTK=',NTK,'jk=',jk;STOP
endif 
end subroutine
!--
subroutine rd_dat_nkm 
implicit none 
integer::ik 
OPEN(132,FILE='./dir-wfn/dat.nkm') 
allocate(NGI(Nk_irr));NGI(:)=0
rewind(132)
do ik=1,Nk_irr 
 read(132,*)NGI(ik) 
enddo 
close(132) 
NTG=maxval(abs(NGI(:))) 
end subroutine
!--
subroutine rd_dat_kg 
implicit none 
integer::i,ig,ik,jk,iop,i1,j1,k1
integer::algn235
real(8)::ktmp(3)
real(8)::tmp,d1,d2,d3,qwf   
real(8)::h1(3),h2(3),h3(3)   
integer::NG_for_psi 
integer,allocatable::KGtmp(:,:)!KGtmp(3,NTG) 
real(8),allocatable::LKGI(:,:)!LKGI(NTG,Nk_irr) 
!
OPEN(104,FILE='./dir-wfn/dat.kg') 
rewind(104) 
allocate(KGI(3,NTG,Nk_irr));KGI(:,:,:)=0 
do ik=1,Nk_irr 
 read(104,*)NG_for_psi 
 !NGI(ik)=NG_for_psi 
 do ig=1,NGI(ik)!NG_for_psi 
  read(104,*)(KGI(i,ig,ik),i=1,3) 
 enddo 
enddo  
close(104)
!--
allocate(LKGI(NTG,Nk_irr));LKGI=0.0d0 
do ik=1,Nk_irr 
 do ig=1,NGI(ik) 
  ktmp(1)=(SKI(1,ik)+dble(KGI(1,ig,ik)))*b1(1)+(SKI(2,ik)+dble(KGI(2,ig,ik)))*b2(1)+(SKI(3,ik)+dble(KGI(3,ig,ik)))*b3(1) 
  ktmp(2)=(SKI(1,ik)+dble(KGI(1,ig,ik)))*b1(2)+(SKI(2,ik)+dble(KGI(2,ig,ik)))*b2(2)+(SKI(3,ik)+dble(KGI(3,ig,ik)))*b3(2) 
  ktmp(3)=(SKI(1,ik)+dble(KGI(1,ig,ik)))*b1(3)+(SKI(2,ik)+dble(KGI(2,ig,ik)))*b2(3)+(SKI(3,ik)+dble(KGI(3,ig,ik)))*b3(3) 
  LKGI(ig,ik)=ktmp(1)**2+ktmp(2)**2+ktmp(3)**2
 enddo!ig 
 !write(6,*)maxval(LKGI(:,ik))
enddo!ik 
!write(6,*) 
!write(6,*)maxval(LKGI(:,:))
Ecut_for_psi=maxval(LKGI(:,:))+1.0d-8
deallocate(LKGI) 
!--
L1=maxval(abs(KGI(1,:,:)))+1;write(6,*)'L1=',L1 
L2=maxval(abs(KGI(2,:,:)))+1;write(6,*)'L2=',L2 
L3=maxval(abs(KGI(3,:,:)))+1;write(6,*)'L3=',L3 
allocate(packing(-L1:L1,-L2:L2,-L3:L3,Nk_irr)); packing(:,:,:,:)=0 
do ik=1,Nk_irr 
 do ig=1,NGI(ik) 
  i1=KGI(1,ig,ik);j1=KGI(2,ig,ik);k1=KGI(3,ig,ik) 
  packing(i1,j1,k1,ik)=ig 
 enddo 
enddo 
!--
!fft grid
tmp=dsqrt(dot_product(a1,a1))
h1(:)=a1(:)/tmp 
tmp=dsqrt(dot_product(a2,a2))
h2(:)=a2(:)/tmp 
tmp=dsqrt(dot_product(a3,a3)) 
h3(:)=a3(:)/tmp 
d1=abs(dot_product(b1,h1)) 
d2=abs(dot_product(b2,h2)) 
d3=abs(dot_product(b3,h3)) 
!write(6,*)'d1=',d1
!write(6,*)'d2=',d2
!write(6,*)'d3=',d3 
qwf=2.0d0*dsqrt(Ecut_for_psi) 
!nwx2=algn235(int(qwf/d1)+1,1)!error 
!nwy2=algn235(int(qwf/d2)+1,1)!error
!nwz2=algn235(int(qwf/d3)+1,1)!error
nwx2=algn235(int(qwf/d1)+1) 
nwy2=algn235(int(qwf/d2)+1) 
nwz2=algn235(int(qwf/d3)+1) 
!nfft1=nwx2+1
!nfft2=nwy2+1
!nfft3=nwz2+1
!Nl123=nfft1*nfft2*nfft3 
!--
!NG0,KG0 
allocate(NG0(NTK));NG0(:)=0
allocate(KG0(3,NTG,NTK));KG0(:,:,:)=0 
allocate(KGtmp(3,NTG));KGtmp(:,:)=0 
!---
do jk=1,NTK 
 if(trs(jk)==1) then 
  ik=numirr(jk)
  iop=numrot(jk) 
  ktmp(1)=dble(rg(1,1,iop))*SKI(1,ik)+dble(rg(1,2,iop))*SKI(2,ik)+dble(rg(1,3,iop))*SKI(3,ik)+dble(RW(1,jk)) 
  ktmp(2)=dble(rg(2,1,iop))*SKI(1,ik)+dble(rg(2,2,iop))*SKI(2,ik)+dble(rg(2,3,iop))*SKI(3,ik)+dble(RW(2,jk))  
  ktmp(3)=dble(rg(3,1,iop))*SKI(1,ik)+dble(rg(3,2,iop))*SKI(2,ik)+dble(rg(3,3,iop))*SKI(3,ik)+dble(RW(3,jk))  
!--
!20170316 KG0s are generated by original KGI
!  NG0(jk)=NGI(ik) 
!  do ig=1,NG0(jk)  
!   i1=rg(1,1,iop)*KGI(1,ig,ik)+rg(1,2,iop)*KGI(2,ig,ik)+rg(1,3,iop)*KGI(3,ig,ik)-RW(1,jk) 
!   i2=rg(2,1,iop)*KGI(1,ig,ik)+rg(2,2,iop)*KGI(2,ig,ik)+rg(2,3,iop)*KGI(3,ig,ik)-RW(2,jk) 
!   i3=rg(3,1,iop)*KGI(1,ig,ik)+rg(3,2,iop)*KGI(2,ig,ik)+rg(3,3,iop)*KGI(3,ig,ik)-RW(3,jk) 
!   KG0(1,ig,jk)=i1 
!   KG0(2,ig,jk)=i2 
!   KG0(3,ig,jk)=i3 
!   phase=tpi*(dble(i1)*dble(pg(1,iop))+dble(i2)*dble(pg(2,iop))+dble(i3)*dble(pg(3,iop))) 
!   pf=exp(-ci*phase/dble(nnp)) 
!   C0(ig,:,jk)=CIR(ig,:,ik)*pf 
!  enddo!ig  
! write(6,'(i5,3f15.10,i5)') jk,ktmp(:),trs(jk)  
!--
  call make_KG0(NTG,b1(1),b2(1),b3(1),Ecut_for_psi,ktmp(1),ktmp(2),ktmp(3),KG0(1,1,jk),NG_for_psi)
  if(NG_for_psi/=NGI(ik)) then 
   write(6,*)'ERROR; STOP; NG_for_psi should be NGI(ik)'   
   write(6,*)'NG_for_psi=',NG_for_psi,'NG0(ik)=',NGI(ik)
   write(6,*)'ik,jk',ik,jk;STOP
  endif 
  NG0(jk)=NG_for_psi  
!--
 elseif(trs(jk)==-1)then  
  ik=numirr(jk);iop=numrot(jk) 
  ktmp(1)=rg(1,1,iop)*SKI(1,ik)+rg(1,2,iop)*SKI(2,ik)+rg(1,3,iop)*SKI(3,ik)+dble(RW(1,jk)) 
  ktmp(2)=rg(2,1,iop)*SKI(1,ik)+rg(2,2,iop)*SKI(2,ik)+rg(2,3,iop)*SKI(3,ik)+dble(RW(2,jk))  
  ktmp(3)=rg(3,1,iop)*SKI(1,ik)+rg(3,2,iop)*SKI(2,ik)+rg(3,3,iop)*SKI(3,ik)+dble(RW(3,jk))  
  KGtmp(:,:)=0 
  call make_KG0(NTG,b1(1),b2(1),b3(1),Ecut_for_psi,ktmp(1),ktmp(2),ktmp(3),KGtmp(1,1),NG_for_psi)
  if(NG_for_psi/=NGI(ik))then 
   write(6,*)'ERROR; STOP; NG_for_psi should be NGI(ik)'   
   write(6,*)'NG_for_psi=',NG_for_psi,'NGI(ik)=',NGI(ik);STOP
  endif 
  NG0(jk)=NG_for_psi  
  KG0(:,:,jk)=-KGtmp(:,:)!notice on '-' sign 
 endif 
enddo!jk 
deallocate(KGtmp) 
end subroutine
!--
subroutine rd_dat_eigenvalue 
implicit none 
integer::ik,ib 
OPEN(111,FILE='./dir-wfn/dat.eigenvalue') 
rewind(111)
read(111,*)NTB 
allocate(E_EIGI(NTB,Nk_irr));E_EIGI=0.0d0
do ik=1,Nk_irr 
 do ib=1,NTB  
  read(111,*)E_EIGI(ib,ik)
 enddo!ib
enddo!ik          
end subroutine
!--
subroutine rd_dat_wavefunction 
implicit none 
integer::ik,ib,ig 
OPEN(102,FILE='./dir-wfn/dat.wfn',FORM='unformatted') 
rewind(102)
read(102)ncomp 
if(ncomp/=1)then 
 write(6,*)'This program not suport ncomp/=1; then stop'
 stop
endif 
allocate(CIR(NTG,NTB,Nk_irr));CIR=0.0d0  
do ik=1,Nk_irr 
 do ib=1,NTB 
  !
  !20180922 
  !
  !read(102)(CIR(ig,ib,ik),ig=1,NGI(ik))
  !
  allocate(CIRtmp(NTG));CIRtmp=0.0d0  
  read(102)(CIRtmp(ig),ig=1,NGI(ik))
  CIR(:,ib,ik)=CIRtmp(:) 
  deallocate(CIRtmp) 
  !
  !
 enddo!ib 
enddo!ik          
close(102) 
end subroutine
!--
end module 
