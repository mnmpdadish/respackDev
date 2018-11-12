module m_rd_dat_eps 
use m_rd_dat_wfn
implicit none
public::rd_dat_chi_cutoff
public::rd_dat_ttrhdrn 
public::rd_dat_wgrid 
public::rd_dat_sq 
public::rd_dat_eps 
!chi_cutoff(300)  
real(8),public::Ecut_for_eps
!ttrhdrn(302)  
real(8),public::idlt!Green's function delt (au)!ttrhdrn
real(8),public::dmna!dmna (au)!ttrhdrn
real(8),public::dmnr!dmnr (au)!ttrhdrn
!wgrid(135)  
integer,public::Num_freq_grid
integer,public::ne 
complex(8),public,allocatable::em(:)!em(ne)
complex(8),allocatable::pole_of_chi(:)!pole_of_chi(ne)  
complex(8),allocatable::mat_b(:,:)!mat_b(ne,ne) 
!sq(301)  
integer,public::Nq_irr!TOTAL NUMBER OF irreducible q POINTS 
integer,public::NTQ!TOTAL NUMBER OF Q POINTS IN MK MESHES 
integer,public::NTGQ!TOTAL NUMBER OF G VECTORS FOR EPSILON
integer,public::Lq1,Lq2,Lq3  
real(8),public,allocatable::SQI(:,:)!SQI(3,Nq_irr) 
real(8),public,allocatable::SQ(:,:)!SQ(3,NTQ)  
integer,public,allocatable::numirrq(:)!numirrq(NTQ)  
integer,public,allocatable::numrotq(:)!numrotq(NTQ)  
integer,public,allocatable::trsq(:)!trsq(NTQ)  
integer,public,allocatable::RWq(:,:)!RWq(3,NTQ)
integer,public,allocatable::LG0(:,:,:)!LG0(3,NTG,NTQ)    
integer,public,allocatable::NGQ_eps(:)!NGQ_eps(NTQ)
integer,public,allocatable::NGQ_psi(:)!NGQ_psi(NTQ)  
integer,public,allocatable::packingq(:,:,:,:)!packingq(-Lq1:Lq1,-Lq2:Lq2,-Lq3:Lq3,Nq_irr) 
!epsqw(600-) 
!complex(8),public,allocatable::epstmp(:,:,:)!epstmp(NTGQ,NTGQ,ne)
!complex(8),public,allocatable::epstmpgm(:,:,:,:)!epstmpgm(NTGQ,NTGQ,ne,3)
complex(4),public,allocatable::epsirr(:,:,:,:)!epsirr(NTGQ,NTGQ,ne,Nq_irr) 
contains
!--
subroutine rd_dat_chi_cutoff
implicit none 
OPEN(300,FILE='./dir-eps/dat.chi_cutoff')
read(300,*) Ecut_for_eps
end subroutine
!--
!20180423 
!--
subroutine rd_dat_ttrhdrn
implicit none 
OPEN(302,FILE='./dir-eps/dat.ttrhdrn')
read(302,*) idlt!ttrhdrn Green's function delt (au)
read(302,*) dmna!(au) 
read(302,*) dmnr!(au) 
end subroutine
!--
subroutine rd_dat_wgrid 
implicit none 
integer::ie,je,ke 
real(8)::x,y
complex(8),allocatable::emp(:)!emp(ne+1)
complex(8),allocatable::mat_c(:,:)!mat_c(ne,ne) 
complex(8)::sum_cmpx 
!
OPEN(135,FILE='./dir-eps/dat.wgrid') 
rewind(135) 
read(135,*) Num_freq_grid
ne=Num_freq_grid
allocate(em(ne)); em(:)=0.0d0 
do ie=1,ne 
 read(135,'(2f15.10)') em(ie)
enddo 
!
allocate(emp(ne+1));emp(:)=0.0d0 
do ie=1,ne
 emp(ie)=em(ie)
enddo
emp(ne+1)=em(ne)+1.3d0*(em(ne)-em(ne-1)) 
!
!pole of chi
!
allocate(pole_of_chi(ne)); pole_of_chi(:)=0.0d0 
do ie=1,ne 
 x=dble((emp(ie+1)+emp(ie))/2.0d0) 
 y=-dble(1.5d0*(emp(ie+1)-emp(ie))) 
 pole_of_chi(ie)=cmplx(x,y) 
enddo 
!
!mat_b
!
allocate(mat_b(ne,ne)); mat_b=0.0d0 
do ie=1,ne
 do je=1,ne 
  mat_b(ie,je)=1.0d0/(em(ie)-pole_of_chi(je))-1.0d0/(em(ie)+pole_of_chi(je)) 
 enddo 
enddo 
!
!check: mat_b * mat_c = 1
!
allocate(mat_c(ne,ne)); mat_c=0.0d0 
mat_c(:,:)=mat_b(:,:)
!
call invmat_complex(ne,mat_b) 
!
do ie=1,ne
 do je=1,ne
  sum_cmpx=0.0d0
  do ke=1,ne
   sum_cmpx=sum_cmpx+mat_b(ie,ke)*mat_c(ke,je)
  enddo 
  !write(6,'(2I5,2F15.10)') ie,je,sum_cmpx 
 enddo 
enddo 
end subroutine
!--
subroutine rd_dat_sq 
implicit none 
integer::iq,i,iik,ierr,chdir,jk,iop,iqir,NG_for_eps,NG_for_psi,ig,i1,j1,k1  
character(99)::dirname 
logical::file_e 
real(8)::q1,q2,q3,ktmp(3)  
integer::RWtmp(3) 
integer,allocatable::KGtmp(:,:)!KGtmp(3,NTG)
!
Nq_irr=Nk_irr 
allocate(SQI(3,Nq_irr));SQI(:,:)=0.0D0 
ierr=CHDIR("./dir-eps") 
do iq=1,Nq_irr 
 write(dirname,"('q',i3.3)")iq 
 ierr=CHDIR(dirname) 
 inquire(file='dat.log.400',exist=file_e) 
 if(file_e)then 
  write(6,*)'dat.log.400 exists in ',trim(dirname) 
 else
  write(6,*)'error: no dat.log.400 in ',trim(dirname) 
 endif 
 OPEN(301,FILE='dat.sq') 
 rewind(301)
 inquire(file='dat.sq',exist=file_e) 
 if(file_e)then 
  read(301,*)(SQI(i,iq),i=1,3) 
 else
  write(6,*)'error: no dat.sq in',trim(dirname) 
 endif 
 close(301) 
 ierr=CHDIR("..") 
enddo!iq 
ierr=CHDIR("..") 
!call system('pwd') 
!--
!gen(SQ,numirrq,numrotq,trsq,RWq) 
NTQ=NTK
allocate(SQ(3,NTQ));SQ(:,:)=0.0d0 
allocate(numirrq(NTQ));numirrq(:)=0 
allocate(numrotq(NTQ));numrotq(:)=0 
allocate(trsq(NTQ));trsq(:)=0
allocate(RWq(3,NTQ));RWq(:,:)=0      
!
do iq=1,Nq_irr 
 SQ(:,iq)=SQI(:,iq) 
 numirrq(iq)=iq
 numrotq(iq)=1
 trsq(iq)=1
 RWq(1:3,iq)=0
enddo 
jk=Nq_irr 
do iq=1,Nq_irr
 do iop=1,Nsymq
  !sym
  ktmp(:)=0.0d0;RWtmp(:)=0  
  ktmp(1)=dble(rg(1,1,iop))*SQI(1,iq)+dble(rg(1,2,iop))*SQI(2,iq)+dble(rg(1,3,iop))*SQI(3,iq)
  ktmp(2)=dble(rg(2,1,iop))*SQI(1,iq)+dble(rg(2,2,iop))*SQI(2,iq)+dble(rg(2,3,iop))*SQI(3,iq)
  ktmp(3)=dble(rg(3,1,iop))*SQI(1,iq)+dble(rg(3,2,iop))*SQI(2,iq)+dble(rg(3,3,iop))*SQI(3,iq)
  !
  call kcheck(ktmp(1),RWtmp(1))!rewind check 
  !
  do iik=1,jk 
   if(abs(SQ(1,iik)-ktmp(1))<1.0d-4.and.abs(SQ(2,iik)-ktmp(2))<1.0d-4.and.abs(SQ(3,iik)-ktmp(3))<1.0d-4) goto 1100
  enddo!iik
  jk=jk+1
  SQ(:,jk)=ktmp(:)
  numirrq(jk)=iq
  numrotq(jk)=iop
  trsq(jk)=1
  RWq(:,jk)=RWtmp(:)
!time-reversal
1100 ktmp(:)=0.0d0;RWtmp(:)=0  
  ktmp(1)=dble(rg(1,1,iop))*SQI(1,iq)+dble(rg(1,2,iop))*SQI(2,iq)+dble(rg(1,3,iop))*SQI(3,iq)
  ktmp(2)=dble(rg(2,1,iop))*SQI(1,iq)+dble(rg(2,2,iop))*SQI(2,iq)+dble(rg(2,3,iop))*SQI(3,iq) 
  ktmp(3)=dble(rg(3,1,iop))*SQI(1,iq)+dble(rg(3,2,iop))*SQI(2,iq)+dble(rg(3,3,iop))*SQI(3,iq) 
  !
  call kcheck_trs(ktmp(1),RWtmp(1))!rewind check 20170321 
  !
  do iik=1,jk
   if(abs(SQ(1,iik)-(-ktmp(1)))<1.0d-4.and.abs(SQ(2,iik)-(-ktmp(2)))<1.0d-4.and.abs(SQ(3,iik)-(-ktmp(3)))<1.0d-4) goto 2100
  enddo!iik
  jk=jk+1
  SQ(:,jk)=-ktmp(:)
  numirrq(jk)=iq
  numrotq(jk)=iop
  trsq(jk)=-1
  RWq(:,jk)=RWtmp(:)
2100 enddo!iop 
enddo!iq  
!--
if(NTQ/=jk)then 
 write(6,*)'ERROR;STOP;NTQ should be jk'   
 write(6,*)'NTQ=',NTQ,'jk=',jk
 stop 
endif 
!--
!gen(NGQ_eps,NGQ_psi,LG0)
allocate(KGtmp(3,NTG));KGtmp(:,:)=0 
allocate(LG0(3,NTG,NTQ));LG0(:,:,:)=0
allocate(NGQ_eps(NTQ));NGQ_eps(:)=0
allocate(NGQ_psi(NTQ));NGQ_psi(:)=0
!
do iq=1,Nq_irr 
 q1=SQI(1,iq); q2=SQI(2,iq); q3=SQI(3,iq)
 call make_LG0(NTG,b1(1),b2(1),b3(1),Ecut_for_eps,Ecut_for_psi,q1,q2,q3,LG0(1,1,iq),NG_for_eps,NG_for_psi) 
 NGQ_eps(iq)=NG_for_eps
 NGQ_psi(iq)=NG_for_psi  
 write(6,'(i8,3f10.5,a8,i8,a8,i10)') iq,q1,q2,q3,'NGeps',NG_for_eps,'NGpsi',NG_for_psi  
enddo 
do iq=Nq_irr+1,NTQ 
 if(trsq(iq)==1)then 
  iqir=numirrq(iq)
  iop=numrotq(iq) 
  !
  !q1=dble(rg(1,1,iop))*SQI(1,iqir)+dble(rg(1,2,iop))*SQI(2,iqir)+dble(rg(1,3,iop))*SQI(3,iqir)+dble(RWq(1,iq)) 
  !q2=dble(rg(2,1,iop))*SQI(1,iqir)+dble(rg(2,2,iop))*SQI(2,iqir)+dble(rg(2,3,iop))*SQI(3,iqir)+dble(RWq(2,iq))  
  !q3=dble(rg(3,1,iop))*SQI(1,iqir)+dble(rg(3,2,iop))*SQI(2,iqir)+dble(rg(3,3,iop))*SQI(3,iqir)+dble(RWq(3,iq))  
  !
  q1=rg(1,1,iop)*SQI(1,iqir)+rg(1,2,iop)*SQI(2,iqir)+rg(1,3,iop)*SQI(3,iqir)+dble(RWq(1,iq)) 
  q2=rg(2,1,iop)*SQI(1,iqir)+rg(2,2,iop)*SQI(2,iqir)+rg(2,3,iop)*SQI(3,iqir)+dble(RWq(2,iq))  
  q3=rg(3,1,iop)*SQI(1,iqir)+rg(3,2,iop)*SQI(2,iqir)+rg(3,3,iop)*SQI(3,iqir)+dble(RWq(3,iq))  
  !
  call make_LG0(NTG,b1(1),b2(1),b3(1),Ecut_for_eps,Ecut_for_psi,q1,q2,q3,LG0(1,1,iq),NG_for_eps,NG_for_psi) 
  NGQ_eps(iq)=NG_for_eps
  NGQ_psi(iq)=NG_for_psi  
  write(6,'(i8,3f10.5,a8,i8,a8,i10)') iq,q1,q2,q3,'NGeps',NG_for_eps,'NGpsi',NG_for_psi  
 elseif(trsq(iq)==-1)then  
  iqir=numirrq(iq)
  iop=numrotq(iq) 
  !
  !q1=dble(rg(1,1,iop))*SQI(1,iqir)+dble(rg(1,2,iop))*SQI(2,iqir)+dble(rg(1,3,iop))*SQI(3,iqir)+dble(RWq(1,iq)) 
  !q2=dble(rg(2,1,iop))*SQI(1,iqir)+dble(rg(2,2,iop))*SQI(2,iqir)+dble(rg(2,3,iop))*SQI(3,iqir)+dble(RWq(2,iq))  
  !q3=dble(rg(3,1,iop))*SQI(1,iqir)+dble(rg(3,2,iop))*SQI(2,iqir)+dble(rg(3,3,iop))*SQI(3,iqir)+dble(RWq(3,iq))  
  !
  q1=rg(1,1,iop)*SQI(1,iqir)+rg(1,2,iop)*SQI(2,iqir)+rg(1,3,iop)*SQI(3,iqir)+dble(RWq(1,iq)) 
  q2=rg(2,1,iop)*SQI(1,iqir)+rg(2,2,iop)*SQI(2,iqir)+rg(2,3,iop)*SQI(3,iqir)+dble(RWq(2,iq))  
  q3=rg(3,1,iop)*SQI(1,iqir)+rg(3,2,iop)*SQI(2,iqir)+rg(3,3,iop)*SQI(3,iqir)+dble(RWq(3,iq))  
  !
  KGtmp(:,:)=0 
  call make_LG0(NTG,b1(1),b2(1),b3(1),Ecut_for_eps,Ecut_for_psi,q1,q2,q3,KGtmp(1,1),NG_for_eps,NG_for_psi) 
  NGQ_eps(iq)=NG_for_eps 
  NGQ_psi(iq)=NG_for_psi  
  write(6,'(i8,3f10.5,a8,i8,a8,i10)') iq,q1,q2,q3,'NGeps',NG_for_eps,'NGpsi',NG_for_psi  
  LG0(:,:,iq)=-KGtmp(:,:) 
 endif 
enddo 
!--
Lq1=maxval(abs(LG0(1,:,:)))+1
Lq2=maxval(abs(LG0(2,:,:)))+1
Lq3=maxval(abs(LG0(3,:,:)))+1
allocate(packingq(-Lq1:Lq1,-Lq2:Lq2,-Lq3:Lq3,Nq_irr));packingq(:,:,:,:)=0 
do iq=1,Nq_irr 
 do ig=1,NGQ_eps(iq) 
  i1=LG0(1,ig,iq)
  j1=LG0(2,ig,iq)
  k1=LG0(3,ig,iq) 
  packingq(i1,j1,k1,iq)=ig 
  !write(6,'(a,5i10)')'iq,ig,i1,j1,k1',iq,ig,i1,j1,k1
 enddo 
enddo 
!-- 
NTGQ=maxval(NGQ_eps(:))
write(6,*)'Lq1=',Lq1 
write(6,*)'Lq2=',Lq2 
write(6,*)'Lq3=',Lq3 
write(6,*)'NTGQ=',NTGQ  
!--
do iq=1,NTQ 
 write(6,'(i5,3f15.10,6i5)')iq,(SQ(i1,iq),i1=1,3),numirrq(iq),numrotq(iq),trsq(iq),(RWq(j1,iq),j1=1,3) 
enddo  
!--
end subroutine
!--
subroutine rd_dat_eps 
implicit none 
integer::ierr,chdir,iq,iqgm,ix,file_num,NG_for_eps,ig,jg,ie
character(99)::filename,dirname 
complex(8),allocatable::epstmp(:,:,:)!epstmp(NTGQ,NTGQ,ne)
complex(8),allocatable::epstmpgm(:,:,:,:)!epstmpgm(NTGQ,NTGQ,ne,3)
!--
!OPEN(600-,R,FILE='dat.epsqw',FORM='unformatted') 
allocate(epstmp(NTGQ,NTGQ,ne));epstmp(:,:,:)=0.0d0!real8
allocate(epstmpgm(NTGQ,NTGQ,ne,3));epstmpgm(:,:,:,:)=0.0d0!real8
allocate(epsirr(NTGQ,NTGQ,ne,Nq_irr));epsirr(:,:,:,:)=0.0d0!real4 
!
ierr=CHDIR("./dir-eps") 
call system('pwd') 
do iq=1,Nq_irr 
 if(abs(SKI(1,iq))<1.0d-5.and.abs(SKI(2,iq))<1.0d-5.and.abs(SKI(3,iq))<1.0d-5)then 
  write(dirname,"('q',i3.3)")iq
  iqgm=iq 
  ierr=CHDIR(dirname) 
  do ix=1,3 
   file_num=600+(ix-1)  
   write(filename,'("dat.epsqw.",i3.3)')file_num 
   open(file_num,file=filename,form='unformatted') 
   write(6,'(a10,a15,a5,a10)')'read: ',trim(filename),'in ',trim(dirname)  
   rewind(file_num) 
   NG_for_eps=NGQ_eps(iq)
   read(file_num)(((epstmpgm(ig,jg,ie,ix),ig=1,NG_for_eps),jg=1,NG_for_eps),ie=1,ne)
  enddo!ix 
  !epsirr: real4
  !epstmp: real8 
  epsirr(:,:,:,iq)=(epstmpgm(:,:,:,1)+epstmpgm(:,:,:,2)+epstmpgm(:,:,:,3))/3.0d0
 else
  write(dirname,"('q',i3.3)")iq
  ierr=CHDIR(dirname) 
  file_num=600 
  write(filename,'("dat.epsqw.",i3.3)')file_num 
  open(file_num,file=filename,form='unformatted') 
  write(6,'(a10,a15,a5,a10)')'read: ',trim(filename),'in ',trim(dirname)  
  rewind(file_num) 
  NG_for_eps=NGQ_eps(iq)
  read(file_num)(((epstmp(ig,jg,ie),ig=1,NG_for_eps),jg=1,NG_for_eps),ie=1,ne)
  !epsirr: real4 
  !epstmp: real8 
  epsirr(:,:,:,iq)=epstmp(:,:,:)
 endif!q=0 or not 
 ierr=CHDIR("..") 
enddo!iq   
ierr=CHDIR("..") 
call system('pwd') 
!
!do iq=1,Nq_irr 
! NG_for_eps=NGQ_eps(iq)
! write(6,*) NG_for_eps 
! do ig=1,NG_for_eps 
!  write(6,*) epsirr(ig,ig,100,iq) 
! enddo  
!enddo  
!
!do iq=1,Nq_irr 
! NG_for_eps=NGQ_eps(iq) 
! do ie=1,nen 
!  do igL=1,NG_for_eps 
!   epsirr(igL,igL,ie,iq)=epsirr(igL,igL,ie,iq)-1.0d0    
!  enddo 
! enddo 
!enddo 
!write(6,*) 'FINISH READING EPSILON'
!
end subroutine
!--
end module 
