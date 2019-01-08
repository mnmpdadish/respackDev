subroutine kcheck(ktmp,RWtmp) 
implicit none 
real(8),intent(inout)::ktmp(3)
integer,intent(out)::RWtmp(3) 
real(8),parameter::dlt_BZ=1.0d-6 
!
if(ktmp(1)>1.50d0+dlt_BZ)then 
 ktmp(1)=ktmp(1)-2.0d0
 RWtmp(1)=-2
endif 
if(ktmp(1)>0.50d0+dlt_BZ)then 
 ktmp(1)=ktmp(1)-1.0d0
 RWtmp(1)=-1
endif 
if(ktmp(1)<=-1.50d0+dlt_BZ)then 
 ktmp(1)=ktmp(1)+2.0d0
 RWtmp(1)=2 
endif 
if(ktmp(1)<=-0.50d0+dlt_BZ)then 
 ktmp(1)=ktmp(1)+1.0d0
 RWtmp(1)=1
endif 
!
if(ktmp(2)>1.50d0+dlt_BZ)then 
 ktmp(2)=ktmp(2)-2.0d0
 RWtmp(2)=-2 
endif 
if(ktmp(2)>0.50d0+dlt_BZ)then 
 ktmp(2)=ktmp(2)-1.0d0
 RWtmp(2)=-1
endif 
if(ktmp(2)<=-1.50d0+dlt_BZ)then 
 ktmp(2)=ktmp(2)+2.0d0
 RWtmp(2)=2 
endif 
if(ktmp(2)<=-0.50d0+dlt_BZ)then 
 ktmp(2)=ktmp(2)+1.0d0
 RWtmp(2)=1
endif 
!
if(ktmp(3)>1.50d0+dlt_BZ)then 
 ktmp(3)=ktmp(3)-2.0d0
 RWtmp(3)=-2 
endif 
if(ktmp(3)>0.50d0+dlt_BZ)then 
 ktmp(3)=ktmp(3)-1.0d0
 RWtmp(3)=-1
endif 
if(ktmp(3)<=-1.50d0+dlt_BZ)then 
 ktmp(3)=ktmp(3)+2.0d0
 RWtmp(3)=2 
endif 
if(ktmp(3)<=-0.50d0+dlt_BZ)then 
 ktmp(3)=ktmp(3)+1.0d0
 RWtmp(3)=1
endif 
return
end 
!
subroutine kcheck_trs(ktmp,RWtmp) 
implicit none 
real(8),intent(inout)::ktmp(3)
integer,intent(out)::RWtmp(3) 
real(8),parameter::dlt_BZ=-1.0d-6 
!
if(ktmp(1)>=1.50d0+dlt_BZ)then 
 ktmp(1)=ktmp(1)-2.0d0
 RWtmp(1)=-2
endif 
if(ktmp(1)>=0.50d0+dlt_BZ)then 
 ktmp(1)=ktmp(1)-1.0d0
 RWtmp(1)=-1
endif 
if(ktmp(1)<-1.50d0+dlt_BZ)then 
 ktmp(1)=ktmp(1)+2.0d0
 RWtmp(1)=2 
endif 
if(ktmp(1)<-0.50d0+dlt_BZ)then 
 ktmp(1)=ktmp(1)+1.0d0
 RWtmp(1)=1
endif 
!
if(ktmp(2)>=1.50d0+dlt_BZ)then 
 ktmp(2)=ktmp(2)-2.0d0
 RWtmp(2)=-2 
endif 
if(ktmp(2)>=0.50d0+dlt_BZ)then 
 ktmp(2)=ktmp(2)-1.0d0
 RWtmp(2)=-1
endif 
if(ktmp(2)<-1.50d0+dlt_BZ)then 
 ktmp(2)=ktmp(2)+2.0d0
 RWtmp(2)=2 
endif 
if(ktmp(2)<-0.50d0+dlt_BZ)then 
 ktmp(2)=ktmp(2)+1.0d0
 RWtmp(2)=1
endif 
!
if(ktmp(3)>=1.50d0+dlt_BZ)then 
 ktmp(3)=ktmp(3)-2.0d0
 RWtmp(3)=-2 
endif 
if(ktmp(3)>=0.50d0+dlt_BZ)then 
 ktmp(3)=ktmp(3)-1.0d0
 RWtmp(3)=-1
endif 
if(ktmp(3)<-1.50d0+dlt_BZ)then 
 ktmp(3)=ktmp(3)+2.0d0
 RWtmp(3)=2 
endif 
if(ktmp(3)<-0.50d0+dlt_BZ)then 
 ktmp(3)=ktmp(3)+1.0d0
 RWtmp(3)=1
endif 
return
end 
!
subroutine make_KG0(NTG,b1,b2,b3,Gcut,q1,q2,q3,KG0,NG) 
implicit none 
integer,intent(in)::NTG
real(8),intent(in)::b1(3),b2(3),b3(3) 
real(8),intent(in)::Gcut,q1,q2,q3        
integer,intent(out)::KG0(3,NTG),NG 
integer::igL,igL1,igL2,igL3
real(8)::qgL(3),qgL2  
integer,parameter::NGL1=150!100
integer,parameter::NGL2=150!100 
integer,parameter::NGL3=150!100
!
igL=0
do igL1=-NGL1,NGL1 
 do igL2=-NGL2,NGL2 
  do igL3=-NGL3,NGL3 
   qgL(:)=(q1+dble(igL1))*b1(:)+(q2+dble(igL2))*b2(:)+(q3+dble(igL3))*b3(:)    
   qgL2=qgL(1)**2+qgL(2)**2+qgL(3)**2
   if(qgL2<=Gcut)then 
    igL=igL+1 
    KG0(1,igL)=igL1
    KG0(2,igL)=igL2
    KG0(3,igL)=igL3 
   endif  
  enddo 
 enddo 
enddo 
NG=igL 
!
RETURN 
END 
!
subroutine est_nkbi(N,SK,nkb1,nkb2,nkb3)  
implicit none 
integer::N,nkb1,nkb2,nkb3,NTK  
real(8)::SK(3,N) 
integer::i 
real(8)::x 
!
x=1.0d0 
do i=1,N
 if(abs(SK(1,i))<1.0d-7)cycle 
 if(abs(SK(1,i))<x)then 
  x=abs(SK(1,i))  
 endif 
enddo    
nkb1=nint(1.0d0/x)  
!
x=1.0d0 
do i=1,N
 if(abs(SK(2,i))<1.0d-7)cycle 
 if(abs(SK(2,i))<x)then 
  x=abs(SK(2,i))  
 endif 
enddo    
nkb2=nint(1.0d0/x)  
!
x=1.0d0 
do i=1,N
 if(abs(SK(3,i))<1.0d-7)cycle 
 if(abs(SK(3,i))<x)then 
  x=abs(SK(3,i))  
 endif 
enddo    
nkb3=nint(1.0d0/x)  
!
NTK=nkb1*nkb2*nkb3 
!
!write(6,'(a24,4i10)')'nkb1,nkb2,nkb3,NTK=',nkb1,nkb2,nkb3,NTK  
!
return 
end 
!
subroutine est_NTK(Nk_irr,Nsymq,SKI,rg,NTK)
!
!20180301 
!
implicit none 
integer::Nk_irr,Nsymq,N  
real(8)::SKI(3,Nk_irr) 
integer::rg(3,3,Nsymq) 
real(8),allocatable::SK0(:,:)!SK0(3,N) 
real(8)::ktmp(3)
integer::RWtmp(3)
integer::jk,ik,iop,iik 
integer::NTK 
!
N=Nk_irr*Nsymq*2
!
!write(6,*)'N=',N
!
!SK0
!
allocate(SK0(3,N));SK0(:,:)=0.0d0
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
2000 enddo!iop 
enddo!ik 
!
NTK=jk 
if(NTK>N)then 
 write(6,*)'Estimated NTK is too large; stop' 
 write(6,*)'NTK, N=',NTK, N 
 stop
endif 
!
!write(6,*)'Estimated NTK=',NTK 
!
return
end
!
subroutine est_latparam(a1,a2,a3,a,b,c,alp,bet,gmm)   
implicit none 
real(8)::a1(3),a2(3),a3(3) 
real(8)::a,b,c,alp,bet,gmm  
real(8),parameter::pi=DACOS(-1.0d0)
!
a=0.0d0 
b=0.0d0 
c=0.0d0 
a=dsqrt(a1(1)**2+a1(2)**2+a1(3)**2) 
b=dsqrt(a2(1)**2+a2(2)**2+a2(3)**2) 
c=dsqrt(a3(1)**2+a3(2)**2+a3(3)**2) 
alp=(a2(1)*a3(1)+a2(2)*a3(2)+a2(3)*a3(3))/b/c 
bet=(a3(1)*a1(1)+a3(2)*a1(2)+a3(3)*a1(3))/c/a
gmm=(a1(1)*a2(1)+a1(2)*a2(2)+a1(3)*a2(3))/a/b 
alp=dacos(alp)*180.0d0/pi  
bet=dacos(bet)*180.0d0/pi  
gmm=dacos(gmm)*180.0d0/pi  
!
!write(6,'(6f15.10)') a,b,c,alp,bet,gmm 
!
return 
end 
!
integer function algn235(inr)
implicit none
integer,intent(in)::inr
integer:: nr,m2,m3,m5,info
nr=inr
call fctck(nr,m2,m3,m5,info)
do while (info .eq. 1)
   nr = nr + 1
   call fctck(nr,m2,m3,m5,info)
end do
algn235 = nr 
return
end function algn235
!
subroutine fctck(n,m2,m3,m5,info)
implicit none
integer,intent(in):: n
integer,intent(out):: m2, m3, m5, info
integer:: i
i=n
m2 = 0
m3 = 0
m5 = 0
info = 0
do while (i .ne. 1)
   if (mod(i,2) .eq. 0) then
      m2 = m2 + 1
      i = i / 2
   else if (mod(i,3) .eq. 0) then
      m3 = m3 + 1
      i = i / 3
   else if (mod(i,5) .eq. 0) then
      m5 = m5 + 1
      i = i / 5
   else
      info = 1
      exit
   end if
end do
return
end subroutine fctck
!
!subroutine make_LG0(NTG,b1,b2,b3,Gcut_for_eps,Gcut_for_psi,q1,q2,q3,LG0,NG_for_eps,NG_for_psi)
!implicit none 
!integer::NTG,igL,igL1,igL2,igL3
!integer::NG_for_eps,NG_for_psi            
!integer::LG0(3,NTG)    
!real(8)::b1(3),b2(3),b3(3) 
!real(8)::Gcut_for_eps,Gcut_for_psi  
!real(8)::q1,q2,q3,qgL2,qgL(3)
!integer,parameter::NGL1=100 
!integer,parameter::NGL2=100 
!integer,parameter::NGL3=100 
!!
!igL=0
!do igL1=-NGL1,NGL1 
! do igL2=-NGL2,NGL2 
!  do igL3=-NGL3,NGL3 
!   qgL(:)=(q1+dble(igL1))*b1(:)+(q2+dble(igL2))*b2(:)+(q3+dble(igL3))*b3(:)    
!   qgL2=qgL(1)**2+qgL(2)**2+qgL(3)**2
!   if(qgL2<=Gcut_for_eps)then 
!    igL=igL+1 
!    LG0(1,igL)=igL1
!    LG0(2,igL)=igL2
!    LG0(3,igL)=igL3 
!    !write(6,*) igL
!   endif  
!  enddo 
! enddo 
!enddo 
!NG_for_eps=igL 
!!
!do igL1=-NGL1,NGL1 
! do igL2=-NGL2,NGL2 
!  do igL3=-NGL3,NGL3 
!   qgL(:)=(q1+dble(igL1))*b1(:)+(q2+dble(igL2))*b2(:)+(q3+dble(igL3))*b3(:)    
!   qgL2=qgL(1)**2+qgL(2)**2+qgL(3)**2
!   if(qgL2>Gcut_for_eps.and.qgL2<=Gcut_for_psi)then 
!    igL=igL+1 
!    LG0(1,igL)=igL1
!    LG0(2,igL)=igL2
!    LG0(3,igL)=igL3 
!    !write(6,*) igL
!   endif  
!  enddo 
! enddo 
!enddo 
!NG_for_psi=igL 
!!
!!do igL=1,NG_for_eps 
!! write(6,*) igL,LG0(:,igL)
!!enddo 
!!
!RETURN 
!END 
!
!
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
case(1)!=== not time-reversal ===   
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
case(-1)!=== time-reversal ===      
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
!
subroutine make_C0_for_given_band(NTG,itrs,NG,KGtmp,RWtmp,rginvtmp,pgtmp,nnp,L1,L2,L3,packtmp,OCCtmp,C0_K) 
implicit none 
integer::NTG,L1,L2,L3,nnp 
integer::itrs 
integer::NG
integer::KGtmp(3,NTG) 
integer::RWtmp(3) 
real(8)::rginvtmp(3,3) 
integer::pgtmp(3) 
integer::packtmp(-L1:L1,-L2:L2,-L3:L3) 
!
!20180922
!
!complex(8)::OCCtmp(NTG) 
complex(4)::OCCtmp(NTG) 
!
integer::ig,jg,i1,i2,i3,j1,j2,j3,k1,k2,k3 
real(8),parameter::pi=dacos(-1.0d0)
real(8),parameter::tpi=2.0d0*pi 
complex(8),parameter::ci=(0.0D0,1.0D0) 
real(8)::phase 
complex(8)::pf 
complex(8)::C0_K(NTG) 
!
C0_K(:)=0.0d0 
select case(itrs) 
case(1)!=== not time-reversal ===      
 do ig=1,NG 
  i1=KGtmp(1,ig); j1=KGtmp(2,ig); k1=KGtmp(3,ig) 
  i2=i1+RWtmp(1); j2=j1+RWtmp(2); k2=k1+RWtmp(3) 
  i3=int(rginvtmp(1,1))*i2+int(rginvtmp(1,2))*j2+int(rginvtmp(1,3))*k2 
  j3=int(rginvtmp(2,1))*i2+int(rginvtmp(2,2))*j2+int(rginvtmp(2,3))*k2 
  k3=int(rginvtmp(3,1))*i2+int(rginvtmp(3,2))*j2+int(rginvtmp(3,3))*k2 
  jg=packtmp(i3,j3,k3) 
  phase=tpi*(dble(i1)*dble(pgtmp(1))+dble(j1)*dble(pgtmp(2))+dble(k1)*dble(pgtmp(3))) 
  pf=exp(-ci*phase/dble(nnp)) 
  C0_K(ig)=OCCtmp(jg)*pf 
 enddo!ig 
case(-1)!=== time-reversal ===      
 do ig=1,NG 
  i1=-KGtmp(1,ig); j1=-KGtmp(2,ig); k1=-KGtmp(3,ig) 
  i2=i1+RWtmp(1); j2=j1+RWtmp(2); k2=k1+RWtmp(3) 
  i3=int(rginvtmp(1,1))*i2+int(rginvtmp(1,2))*j2+int(rginvtmp(1,3))*k2 
  j3=int(rginvtmp(2,1))*i2+int(rginvtmp(2,2))*j2+int(rginvtmp(2,3))*k2 
  k3=int(rginvtmp(3,1))*i2+int(rginvtmp(3,2))*j2+int(rginvtmp(3,3))*k2 
  jg=packtmp(i3,j3,k3) 
  phase=tpi*(dble(i1)*dble(pgtmp(1))+dble(j1)*dble(pgtmp(2))+dble(k1)*dble(pgtmp(3))) 
  pf=exp(-ci*phase/dble(nnp)) 
  C0_K(ig)=OCCtmp(jg)*pf 
 enddo !ig 
 C0_K(:)=conjg(C0_K(:))  
end select 
return 
end 
