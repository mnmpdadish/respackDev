module m_rd_dat_eps 
use m_rd_dat_wfn
implicit none
public::rd_dat_chi_cutoff
public::rd_dat_wgrid 
public::rd_dat_sq 
!chi_cutoff(300)  
real(8),public::Ecut_for_eps
!wgrid(135)  
integer,public::Num_freq_grid
integer,public::ne 
complex(8),public,allocatable::em(:)!em(ne)
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
!epsqw 
contains
!--
subroutine rd_dat_chi_cutoff
implicit none 
OPEN(300,FILE='./dir-eps/dat.chi_cutoff')
read(300,*) Ecut_for_eps
end subroutine
!--
subroutine rd_dat_wgrid 
implicit none 
integer::ie  
OPEN(135,FILE='./dir-eps/dat.wgrid') 
rewind(135) 
read(135,*) Num_freq_grid
ne=Num_freq_grid
allocate(em(ne)); em(:)=0.0d0 
do ie=1,ne 
 read(135,'(2f15.10)') em(ie)
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
  call kcheck(ktmp(1),RWtmp(1))!rewind check 
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
  call kcheck_trs(ktmp(1),RWtmp(1))!rewind check 20170321 
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
 !write(6,'(i8,3f10.5,a8,i8,a8,i10)') iq,q1,q2,q3,'NGeps',NG_for_eps,'NGpsi',NG_for_psi  
enddo 
do iq=Nq_irr+1,NTQ 
 if(trsq(iq)==1)then 
  iqir=numirrq(iq)
  iop=numrotq(iq) 
  q1=dble(rg(1,1,iop))*SQI(1,iqir)+dble(rg(1,2,iop))*SQI(2,iqir)+dble(rg(1,3,iop))*SQI(3,iqir)+dble(RWq(1,iq)) 
  q2=dble(rg(2,1,iop))*SQI(1,iqir)+dble(rg(2,2,iop))*SQI(2,iqir)+dble(rg(2,3,iop))*SQI(3,iqir)+dble(RWq(2,iq))  
  q3=dble(rg(3,1,iop))*SQI(1,iqir)+dble(rg(3,2,iop))*SQI(2,iqir)+dble(rg(3,3,iop))*SQI(3,iqir)+dble(RWq(3,iq))  
  call make_LG0(NTG,b1(1),b2(1),b3(1),Ecut_for_eps,Ecut_for_psi,q1,q2,q3,LG0(1,1,iq),NG_for_eps,NG_for_psi) 
  NGQ_eps(iq)=NG_for_eps
  NGQ_psi(iq)=NG_for_psi  
  !write(6,'(i8,3f10.5,a8,i8,a8,i10)') iq,q1,q2,q3,'NGeps',NG_for_eps,'NGpsi',NG_for_psi  
 elseif(trsq(iq)==-1)then  
  iqir=numirrq(iq)
  iop=numrotq(iq) 
  q1=dble(rg(1,1,iop))*SQI(1,iqir)+dble(rg(1,2,iop))*SQI(2,iqir)+dble(rg(1,3,iop))*SQI(3,iqir)+dble(RWq(1,iq)) 
  q2=dble(rg(2,1,iop))*SQI(1,iqir)+dble(rg(2,2,iop))*SQI(2,iqir)+dble(rg(2,3,iop))*SQI(3,iqir)+dble(RWq(2,iq))  
  q3=dble(rg(3,1,iop))*SQI(1,iqir)+dble(rg(3,2,iop))*SQI(2,iqir)+dble(rg(3,3,iop))*SQI(3,iqir)+dble(RWq(3,iq))  
  KGtmp(:,:)=0 
  call make_LG0(NTG,b1(1),b2(1),b3(1),Ecut_for_eps,Ecut_for_psi,q1,q2,q3,KGtmp(1,1),NG_for_eps,NG_for_psi) 
  NGQ_eps(iq)=NG_for_eps 
  NGQ_psi(iq)=NG_for_psi  
  !write(6,'(i8,3f10.5,a8,i8,a8,i10)') iq,q1,q2,q3,'NGeps',NG_for_eps,'NGpsi',NG_for_psi  
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
 enddo 
enddo 
!-- 
NTGQ=maxval(NGQ_eps(:))
!write(6,*)'Lq1=',Lq1 
!write(6,*)'Lq2=',Lq2 
!write(6,*)'Lq3=',Lq3 
!write(6,*)'NTGQ=',NTGQ  
!--
end subroutine
!--
end module 
