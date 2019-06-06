subroutine make_LG0(NTG,b1,b2,b3,Gcut_for_eps,Gcut_for_psi,q1,q2,q3,LG0,NG_for_eps,NG_for_psi)
  implicit none 
  integer::NTG,igL,igL1,igL2,igL3
  integer::NG_for_eps,NG_for_psi            
  integer::LG0(3,NTG)    
  real(8)::b1(3),b2(3),b3(3) 
  real(8)::Gcut_for_eps,Gcut_for_psi  
  real(8)::q1,q2,q3,qgL2,qgL(3)
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
     if(qgL2<=Gcut_for_eps)then 
      igL=igL+1 
      LG0(1,igL)=igL1
      LG0(2,igL)=igL2
      LG0(3,igL)=igL3 
      !write(6,*) igL
     endif  
    enddo 
   enddo 
  enddo 
  NG_for_eps=igL 
  !
  do igL1=-NGL1,NGL1 
   do igL2=-NGL2,NGL2 
    do igL3=-NGL3,NGL3 
     qgL(:)=(q1+dble(igL1))*b1(:)+(q2+dble(igL2))*b2(:)+(q3+dble(igL3))*b3(:)    
     qgL2=qgL(1)**2+qgL(2)**2+qgL(3)**2
     if(qgL2>Gcut_for_eps.and.qgL2<=Gcut_for_psi)then 
      igL=igL+1 
      LG0(1,igL)=igL1
      LG0(2,igL)=igL2
      LG0(3,igL)=igL3 
      !write(6,*) igL
     endif  
    enddo 
   enddo 
  enddo 
  NG_for_psi=igL 
  !
  !do igL=1,NG_for_eps 
  ! write(6,*) igL,LG0(:,igL)
  !enddo 
  !
  RETURN 
END subroutine  
!
subroutine make_eps(NTG,NTGQ,ne,itrs,NG,LGtmp,RWtmp,rginvtmp,pgtmp,nnp,L1,L2,L3,packtmp,epsirr,epsmk) 
  implicit none 
  integer::NTG,NTGQ,itrs,NG,ne,L1,L2,L3,nnp 
  integer::LGtmp(3,NTG) 
  integer::RWtmp(3) 
  real(8)::rginvtmp(3,3) 
  integer::pgtmp(3) 
  integer::packtmp(-L1:L1,-L2:L2,-L3:L3) 
  complex(4)::epsirr(NTGQ,NTGQ,ne) 
  integer::ig,jg,i1,i2,i3,j1,j2,j3,k1,k2,k3,iig,jjg 
  real(8)::phase 
  complex(8)::pf1,pf2 
  complex(4)::epsmk(NTGQ,NTGQ,ne) 
  !
  real(8),parameter::pi=dacos(-1.0d0)
  real(8),parameter::tpi=2.0d0*pi 
  complex(8),parameter::ci=(0.0D0,1.0D0) 
  !
  epsmk(:,:,:)=0.0d0 
  select case(itrs) 
  case(1)!=== not time-reversal ===    
  do ig=1,NG 
   i1=LGtmp(1,ig); j1=LGtmp(2,ig); k1=LGtmp(3,ig) 
   i2=i1+RWtmp(1); j2=j1+RWtmp(2); k2=k1+RWtmp(3) 
   i3=int(rginvtmp(1,1))*i2+int(rginvtmp(1,2))*j2+int(rginvtmp(1,3))*k2 
   j3=int(rginvtmp(2,1))*i2+int(rginvtmp(2,2))*j2+int(rginvtmp(2,3))*k2 
   k3=int(rginvtmp(3,1))*i2+int(rginvtmp(3,2))*j2+int(rginvtmp(3,3))*k2 
   iig=packtmp(i3,j3,k3) 
   phase=tpi*(dble(i1)*dble(pgtmp(1))+dble(j1)*dble(pgtmp(2))+dble(k1)*dble(pgtmp(3))) 
   pf1=exp(-ci*phase/dble(nnp)) 
   do jg=1,NG 
    i1=LGtmp(1,jg); j1=LGtmp(2,jg); k1=LGtmp(3,jg) 
    i2=i1+RWtmp(1); j2=j1+RWtmp(2); k2=k1+RWtmp(3) 
    i3=int(rginvtmp(1,1))*i2+int(rginvtmp(1,2))*j2+int(rginvtmp(1,3))*k2 
    j3=int(rginvtmp(2,1))*i2+int(rginvtmp(2,2))*j2+int(rginvtmp(2,3))*k2 
    k3=int(rginvtmp(3,1))*i2+int(rginvtmp(3,2))*j2+int(rginvtmp(3,3))*k2 
    jjg=packtmp(i3,j3,k3) 
    phase=tpi*(dble(i1)*dble(pgtmp(1))+dble(j1)*dble(pgtmp(2))+dble(k1)*dble(pgtmp(3))) 
    pf2=exp(ci*phase/dble(nnp)) 
    !write(6,'(a,x,2i5)') 'iig jjg',iig,jjg  
    epsmk(ig,jg,:)=epsirr(iig,jjg,:)*pf1*pf2 
   enddo!jg 
  enddo!ig 
  case(-1)!=== time-reversal ===*      
  do ig=1,NG 
   i1=-LGtmp(1,ig); j1=-LGtmp(2,ig); k1=-LGtmp(3,ig) 
   i2=i1+RWtmp(1); j2=j1+RWtmp(2); k2=k1+RWtmp(3) 
   i3=int(rginvtmp(1,1))*i2+int(rginvtmp(1,2))*j2+int(rginvtmp(1,3))*k2 
   j3=int(rginvtmp(2,1))*i2+int(rginvtmp(2,2))*j2+int(rginvtmp(2,3))*k2 
   k3=int(rginvtmp(3,1))*i2+int(rginvtmp(3,2))*j2+int(rginvtmp(3,3))*k2 
   iig=packtmp(i3,j3,k3) 
   phase=tpi*(dble(i1)*dble(pgtmp(1))+dble(j1)*dble(pgtmp(2))+dble(k1)*dble(pgtmp(3))) 
   pf1=exp(-ci*phase/dble(nnp)) 
   do jg=1,NG 
    i1=-LGtmp(1,jg); j1=-LGtmp(2,jg); k1=-LGtmp(3,jg) 
    i2=i1+RWtmp(1); j2=j1+RWtmp(2); k2=k1+RWtmp(3) 
    i3=int(rginvtmp(1,1))*i2+int(rginvtmp(1,2))*j2+int(rginvtmp(1,3))*k2 
    j3=int(rginvtmp(2,1))*i2+int(rginvtmp(2,2))*j2+int(rginvtmp(2,3))*k2 
    k3=int(rginvtmp(3,1))*i2+int(rginvtmp(3,2))*j2+int(rginvtmp(3,3))*k2 
    jjg=packtmp(i3,j3,k3) 
    phase=tpi*(dble(i1)*dble(pgtmp(1))+dble(j1)*dble(pgtmp(2))+dble(k1)*dble(pgtmp(3))) 
    pf2=exp(ci*phase/dble(nnp)) 
    epsmk(ig,jg,:)=epsirr(jjg,iig,:)*pf1*pf2 
   enddo!jg 
  enddo!ig 
  end select 
  !
  !write(6,*)'finish make_eps'
  !
  return 
end subroutine 
