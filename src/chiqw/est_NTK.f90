subroutine est_nkbi(N,SK,nkb1,nkb2,nkb3)  
implicit none 
integer::N,nkb1,nkb2,nkb3,NTK  
real(8)::SK(3,N) 
integer::i 
real(8)::x 
!---
x=1.0d0 
do i=1,N
 if(abs(SK(1,i))<1.0d-7) cycle 
 if(abs(SK(1,i))<x) then 
  x=abs(SK(1,i))  
 endif 
enddo    
nkb1=nint(1.0d0/x)  
!---
x=1.0d0 
do i=1,N
 if(abs(SK(2,i))<1.0d-7) cycle 
 if(abs(SK(2,i))<x) then 
  x=abs(SK(2,i))  
 endif 
enddo    
nkb2=nint(1.0d0/x)  
!---
x=1.0d0 
do i=1,N
 if(abs(SK(3,i))<1.0d-7) cycle 
 if(abs(SK(3,i))<x) then 
  x=abs(SK(3,i))  
 endif 
enddo    
nkb3=nint(1.0d0/x)  
!---
NTK=nkb1*nkb2*nkb3 
!write(6,'(a24,4i10)') 'nkb1,nkb2,nkb3,NTK=',nkb1,nkb2,nkb3,NTK  
!---
return 
end 
!--
!20180301 
subroutine est_NTK(Nk_irr,Nsymq,SKI,rg,NTK)
implicit none 
integer::Nk_irr,Nsymq,N  
real(8)::SKI(3,Nk_irr) 
integer::rg(3,3,Nsymq) 
real(8),allocatable::SK0(:,:)!SK0(3,N) 
real(8)::ktmp(3)
integer::RWtmp(3)
integer::jk,ik,iop,iik 
integer::NTK 
!---
N=Nk_irr*Nsymq*2
!write(6,*)'N=',N
!--
!SK0,numirr,numrot,trs,RW
allocate(SK0(3,N));SK0(:,:)=0.0d0
!allocate(numirr(NTK));numirr(:)=0
!allocate(numrot(NTK));numrot(:)=0
!allocate(trs(NTK));trs(:)=0
!allocate(RW(3,NTK));RW(:,:)=0
jk=0
do ik=1,Nk_irr
 do iop=1,Nsymq
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
  !numirr(jk)=ik;numrot(jk)=iop;trs(jk)=1;RW(:,jk)=RWtmp(:)
1000 ktmp(:)=0.0d0; RWtmp(:)=0  
  ktmp(1)=dble(rg(1,1,iop))*SKI(1,ik)+dble(rg(1,2,iop))*SKI(2,ik)+dble(rg(1,3,iop))*SKI(3,ik)
  ktmp(2)=dble(rg(2,1,iop))*SKI(1,ik)+dble(rg(2,2,iop))*SKI(2,ik)+dble(rg(2,3,iop))*SKI(3,ik)
  ktmp(3)=dble(rg(3,1,iop))*SKI(1,ik)+dble(rg(3,2,iop))*SKI(2,ik)+dble(rg(3,3,iop))*SKI(3,ik)
  call kcheck_trs(ktmp(1),RWtmp(1))!rewind check modified 20170316  
  do iik=1,jk
   if(abs(SK0(1,iik)-(-ktmp(1)))<1.0d-4.and.abs(SK0(2,iik)-(-ktmp(2)))<1.0d-4.and.abs(SK0(3,iik)-(-ktmp(3)))<1.0d-4) goto 2000
  enddo!iik
  jk=jk+1
  SK0(:,jk)=-ktmp(:) 
  !numirr(jk)=ik;numrot(jk)=iop;trs(jk)=-1;RW(:,jk)=RWtmp(:) 
2000 enddo!iop 
enddo!ik 
!--
NTK=jk 
if(NTK>N)then 
 write(6,*)'Estimated NTK is too large; stop' 
 write(6,*)'NTK, N=',NTK, N 
 stop
endif 
!write(6,*)'Estimated NTK=',NTK 
return
end
