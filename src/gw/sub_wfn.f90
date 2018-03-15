subroutine kcheck(ktmp,RWtmp) 
implicit none 
real(8),intent(inout)::ktmp(3)
integer,intent(out)::RWtmp(3) 
real(8),parameter::dlt_BZ=1.0d-6 
!--
if(ktmp(1)>1.50d0+dlt_BZ)then 
 ktmp(1)=ktmp(1)-2.0d0;RWtmp(1)=-2
endif 
if(ktmp(1)>0.50d0+dlt_BZ)then 
 ktmp(1)=ktmp(1)-1.0d0;RWtmp(1)=-1
endif 
if(ktmp(1)<=-1.50d0+dlt_BZ)then 
 ktmp(1)=ktmp(1)+2.0d0;RWtmp(1)=2 
endif 
if(ktmp(1)<=-0.50d0+dlt_BZ)then 
 ktmp(1)=ktmp(1)+1.0d0;RWtmp(1)=1
endif 
!--
if(ktmp(2)>1.50d0+dlt_BZ)then 
 ktmp(2)=ktmp(2)-2.0d0;RWtmp(2)=-2 
endif 
if(ktmp(2)>0.50d0+dlt_BZ)then 
 ktmp(2)=ktmp(2)-1.0d0;RWtmp(2)=-1
endif 
if(ktmp(2)<=-1.50d0+dlt_BZ)then 
 ktmp(2)=ktmp(2)+2.0d0;RWtmp(2)=2 
endif 
if(ktmp(2)<=-0.50d0+dlt_BZ)then 
 ktmp(2)=ktmp(2)+1.0d0;RWtmp(2)=1
endif 
!--
if(ktmp(3)>1.50d0+dlt_BZ)then 
 ktmp(3)=ktmp(3)-2.0d0;RWtmp(3)=-2 
endif 
if(ktmp(3)>0.50d0+dlt_BZ)then 
 ktmp(3)=ktmp(3)-1.0d0;RWtmp(3)=-1
endif 
if(ktmp(3)<=-1.50d0+dlt_BZ)then 
 ktmp(3)=ktmp(3)+2.0d0;RWtmp(3)=2 
endif 
if(ktmp(3)<=-0.50d0+dlt_BZ)then 
 ktmp(3)=ktmp(3)+1.0d0;RWtmp(3)=1
endif 
return
end 
!--
subroutine kcheck_trs(ktmp,RWtmp) 
implicit none 
real(8),intent(inout)::ktmp(3)
integer,intent(out)::RWtmp(3) 
real(8),parameter::dlt_BZ=-1.0d-6 
!--
if(ktmp(1)>=1.50d0+dlt_BZ)then 
 ktmp(1)=ktmp(1)-2.0d0;RWtmp(1)=-2
endif 
if(ktmp(1)>=0.50d0+dlt_BZ)then 
 ktmp(1)=ktmp(1)-1.0d0;RWtmp(1)=-1
endif 
if(ktmp(1)<-1.50d0+dlt_BZ)then 
 ktmp(1)=ktmp(1)+2.0d0;RWtmp(1)=2 
endif 
if(ktmp(1)<-0.50d0+dlt_BZ)then 
 ktmp(1)=ktmp(1)+1.0d0;RWtmp(1)=1
endif 
!--
if(ktmp(2)>=1.50d0+dlt_BZ)then 
 ktmp(2)=ktmp(2)-2.0d0;RWtmp(2)=-2 
endif 
if(ktmp(2)>=0.50d0+dlt_BZ)then 
 ktmp(2)=ktmp(2)-1.0d0;RWtmp(2)=-1
endif 
if(ktmp(2)<-1.50d0+dlt_BZ)then 
 ktmp(2)=ktmp(2)+2.0d0;RWtmp(2)=2 
endif 
if(ktmp(2)<-0.50d0+dlt_BZ)then 
 ktmp(2)=ktmp(2)+1.0d0;RWtmp(2)=1
endif 
!--
if(ktmp(3)>=1.50d0+dlt_BZ)then 
 ktmp(3)=ktmp(3)-2.0d0;RWtmp(3)=-2 
endif 
if(ktmp(3)>=0.50d0+dlt_BZ)then 
 ktmp(3)=ktmp(3)-1.0d0;RWtmp(3)=-1
endif 
if(ktmp(3)<-1.50d0+dlt_BZ)then 
 ktmp(3)=ktmp(3)+2.0d0;RWtmp(3)=2 
endif 
if(ktmp(3)<-0.50d0+dlt_BZ)then 
 ktmp(3)=ktmp(3)+1.0d0;RWtmp(3)=1
endif 
return
end 
!--
subroutine make_KG0(NTG,b1,b2,b3,Gcut,q1,q2,q3,KG0,NG) 
implicit none 
integer,intent(in)::NTG
real(8),intent(in)::b1(3),b2(3),b3(3) 
real(8),intent(in)::Gcut,q1,q2,q3        
integer,intent(out)::KG0(3,NTG),NG 
integer::igL,igL1,igL2,igL3
real(8)::qgL(3),qgL2  
integer,parameter::NGL1=100
integer,parameter::NGL2=100 
integer,parameter::NGL3=100
!
igL=0
do igL1=-NGL1,NGL1 
 do igL2=-NGL2,NGL2 
  do igL3=-NGL3,NGL3 
   qgL(:)=(q1+dble(igL1))*b1(:)+(q2+dble(igL2))*b2(:)+(q3+dble(igL3))*b3(:)    
   qgL2=qgL(1)**2+qgL(2)**2+qgL(3)**2
   if(qgL2<=Gcut) then 
    igL=igL+1 
    KG0(1,igL)=igL1
    KG0(2,igL)=igL2
    KG0(3,igL)=igL3 
   endif  
  enddo 
 enddo 
enddo 
NG=igL 
!--
RETURN 
END 
!--

subroutine make_C0(NTG,NTB,itrs,NG,KGtmp,RWtmp,rginvtmp,pgtmp,nnp,L1,L2,L3,packtmp,C0irr,C0_K) 
implicit none 
integer::NTG,NTB,itrs,NG,L1,L2,L3,nnp 
integer::KGtmp(3,NTG),RWtmp(3),pgtmp(3) 
integer::packtmp(-L1:L1,-L2:L2,-L3:L3) 
integer::ig,jg,i1,i2,i3,j1,j2,j3,k1,k2,k3,ib  
real(8)::rginvtmp(3,3),phase 
complex(8)::C0irr(NTG,NTB)
complex(8)::pf
complex(8),intent(out)::C0_K(NTG,NTB) 
real(8),parameter::tpi=2.0d0*dacos(-1.0d0)
complex(8),parameter::ci=(0.0D0,1.0D0) 
!--
C0_K(:,:)=0.0d0 
write(6,*)'nnp=',nnp 
select case(itrs) 
case(1)!=== not time-reversal ===*      
 do ig=1,NG 
  i1=KGtmp(1,ig); j1=KGtmp(2,ig); k1=KGtmp(3,ig) 
  i2=i1+RWtmp(1); j2=j1+RWtmp(2); k2=k1+RWtmp(3) 
  i3=int(rginvtmp(1,1))*i2+int(rginvtmp(1,2))*j2+int(rginvtmp(1,3))*k2
  j3=int(rginvtmp(2,1))*i2+int(rginvtmp(2,2))*j2+int(rginvtmp(2,3))*k2 
  k3=int(rginvtmp(3,1))*i2+int(rginvtmp(3,2))*j2+int(rginvtmp(3,3))*k2 
  jg=packtmp(i3,j3,k3) 
  phase=tpi*(dble(i1)*dble(pgtmp(1))+dble(j1)*dble(pgtmp(2))+dble(k1)*dble(pgtmp(3))) 
  pf=exp(-ci*phase/dble(nnp)) 
  write(6,*)'pf=',pf 
  !C0_K(ig,:)=C0irr(jg,:)*pf 
  do ib=1,NTB 
   write(6,*) ig,jg 
   C0_K(ig,ib)=C0irr(jg,ib)*pf 
  enddo!ib
 enddo!ig 
case(-1)!=== time-reversal ===*      
 do ig=1,NG 
  i1=KGtmp(1,ig); j1=KGtmp(2,ig); k1=KGtmp(3,ig) 
  i2=i1+RWtmp(1); j2=j1+RWtmp(2); k2=k1+RWtmp(3) 
  i3=int(rginvtmp(1,1))*i2+int(rginvtmp(1,2))*j2+int(rginvtmp(1,3))*k2
  j3=int(rginvtmp(2,1))*i2+int(rginvtmp(2,2))*j2+int(rginvtmp(2,3))*k2 
  k3=int(rginvtmp(3,1))*i2+int(rginvtmp(3,2))*j2+int(rginvtmp(3,3))*k2 
  jg=packtmp(i3,j3,k3) 
  phase=tpi*(dble(i1)*dble(pgtmp(1))+dble(j1)*dble(pgtmp(2))+dble(k1)*dble(pgtmp(3))) 
  pf=exp(-ci*phase/dble(nnp)) 
  !C0_K(ig,:)=C0irr(jg,:)*pf 
  do ib=1,NTB 
   write(6,*) ig,jg 
   C0_K(ig,ib)=C0irr(jg,ib)*pf 
  enddo!ib
 enddo!ig 
 C0_K(:,:)=conjg(C0_K(:,:))  
end select 
!--
return 
end 
!--
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
