subroutine wrt_mvmc(nkb1,nkb2,nkb3,NTK,Na1,Na2,Na3,n_occ,cut,HR) 
implicit none 
integer::nkb1,nkb2,nkb3,NTK,Na1,Na2,Na3,n_occ 
real(8)::cut 
integer::L1,L2,L3,Nsite 
integer::ia1,ia2,ia3,ja1,ja2,ja3,t1,t2,t3,isite,jsite,ib,jb,i,j,k,is     
complex(8)::HR(n_occ,n_occ,-Na1:Na1,-Na2:Na2,-Na3:Na3)
complex(8),allocatable::H2D(:,:,:,:) 
complex(8)::en!20171107 check misawa 
real(8)::en_abs!20171107 check misawa 
real(8),parameter::au=27.21151d0
!--
L1=nkb1
L2=nkb2
L3=nkb3
Nsite=NTK 
allocate(H2d(n_occ,n_occ,Nsite,Nsite)) 
!--
if(mod(Nsite,2)==0)then!even 
 write(6,*)'Nsite even'
 do ia1=-Na1+1,Na1
  do ia2=-Na2+1,Na2
   do ia3=-Na3+1,Na3 
    isite=(ia3+Na3-1)+(ia2+Na2-1)*L3+(ia1+Na1-1)*L3*L2+1 
    do ja1=-Na1+1,Na1
     do ja2=-Na2+1,Na2
      do ja3=-Na3+1,Na3 
       jsite=(ja3+Na3-1)+(ja2+Na2-1)*L3+(ja1+Na1-1)*L3*L2+1  
       t1=ja1-ia1
       t2=ja2-ia2 
       t3=ja3-ia3 
       call rwindvec_even(t1,Na1,L1)
       call rwindvec_even(t2,Na2,L2)
       call rwindvec_even(t3,Na3,L3)
       H2d(:,:,isite,jsite)=HR(:,:,t1,t2,t3)   
      enddo 
     enddo 
    enddo 
   enddo 
  enddo 
 enddo 
else!odd
 write(6,*)'Nsite odd'
 do ia1=-Na1,Na1
  do ia2=-Na2,Na2
   do ia3=-Na3,Na3 
    isite=(ia3+Na3)+(ia2+Na2)*L3+(ia1+Na1)*L3*L2+1 
    do ja1=-Na1,Na1
     do ja2=-Na2,Na2
      do ja3=-Na3,Na3 
       jsite=(ja3+Na3)+(ja2+Na2)*L3+(ja1+Na1)*L3*L2+1  
       t1=ja1-ia1
       t2=ja2-ia2 
       t3=ja3-ia3 
       call rwindvec_odd(t1,Na1,L1)
       call rwindvec_odd(t2,Na2,L2)
       call rwindvec_odd(t3,Na3,L3)
       H2d(:,:,isite,jsite)=HR(:,:,t1,t2,t3)   
      enddo 
     enddo 
    enddo 
   enddo 
  enddo 
 enddo 
endif 
!--
k=0 
do isite=1,Nsite
 do ib=1,n_occ 
  i=ib+(isite-1)*n_occ 
  do jsite=1,Nsite 
   do jb=1,n_occ 
    j=jb+(jsite-1)*n_occ 
    en_abs=abs(H2d(ib,jb,isite,jsite))!20171107 check misawa  
    if(en_abs<cut)cycle!20171107 check misawa 
     do is=1,2
      k=k+1
     enddo!is 
   enddo!jb 
  enddo!jsite 
 enddo!ib
enddo!isite  
!--
!OPEN(300,W,FILE='ztransfer.def') 
OPEN(300,FILE='./dir-mvmc/ztransfer.def') 
write(300,'(a30)')'=============================='
write(300,'(a10,i20)')'NTransfer',k
write(300,'(a30)')'=============================='
write(300,'(a10,f10.5,a10)')'=== tcut >',cut*au,'eV ======='
write(300,'(a30)')'=============================='
do isite=1,Nsite
 do ib=1,n_occ 
  i=ib+(isite-1)*n_occ 
  do jsite=1,Nsite 
   do jb=1,n_occ 
    j=jb+(jsite-1)*n_occ 
    en=H2d(ib,jb,isite,jsite)!20171107 check misawa  
    en_abs=abs(H2d(ib,jb,isite,jsite))!20171107 check misawa  
    if(en_abs<cut)cycle!20171107 check misawa 
     do is=1,2
     write(300,'(4i5,2f15.10)') i-1,is-1,j-1,is-1,-en*au!20171107 check misawa 
     enddo!is 
   enddo!jb 
  enddo!jsite 
 enddo!ib
enddo!isite  
close(300) 
deallocate(H2d) 
!--
!OPEN(301,W,FILE='mvmc.info') 
OPEN(301,FILE='./dir-mvmc/mvmc.info') 
write(301,'(i5,a45)') L1,':The total number of lattices in a1 direc'
write(301,'(i5,a45)') L2,':The total number of lattices in a2 direc'
write(301,'(i5,a45)') L3,':The total number of lattices in a3 direc'
write(301,'(i5,a45)') n_occ,':The total number of wanniers in unitcell'
write(301,'(i5,a45)') n_occ*Nsite,':The total degree of freedom  in unitcell'
write(301,'(a5,a45)') 'cubic',':Lattice type (default) check yourself   '
close(301)
!--
return 
end 
!
subroutine rwindvec_odd(t,N,L)
implicit none 
integer::t,N,L
if(t<-N) t=t+L 
if(t> N) t=t-L 
return
end 
!
subroutine rwindvec_even(t,N,L)
implicit none 
integer::t,N,L
if(t<-N+1) t=t+L 
if(t> N) t=t-L 
return
end 
