subroutine calc_MAT_VXC(NTK,NWF,NTG,NG0,KG0,vhxc,C0_K,ik,nrx2,nry2,nrz2,nfft1,nfft2,Nl123,&
wfunc,fftwk,fs,MAT_VHXC)
use fft_3d 
implicit none 
type(fft3_struct)::fs 
integer::ik,NTK,NWF,NTG 
integer::NG0(NTK)
integer::KG0(3,NTG,NTK) 
complex(8)::C0_K(NTG,NWF)
integer::nrx2,nry2,nrz2,nfft1,nfft2,Nl123
real(8)::wfunc(Nl123*2)
real(8)::fftwk(Nl123*2) 
real(8)::vhxc(nrx2,nry2,nrz2) 
integer::ig,igb1,igb2,igb3,ind
integer::ir,ir1,ir2,ir3,iw,jw
complex(8)::SUM_CMPX
complex(8),allocatable::u_1D(:)!u_1D(Nl123)  
complex(8),allocatable::u_3D(:,:,:,:)!u_3D(nrx2,nry2,nrz2,NWF)  
complex(8),intent(out)::MAT_VHXC(NWF,NWF) 
!
allocate(u_1D(Nl123));u_1D=0.0d0 
allocate(u_3D(nrx2,nry2,nrz2,NWF));u_3D=0.0d0 
do iw=1,NWF 
 wfunc=0.0d0
 fftwk=0.0d0
 do ig=1,NG0(ik) 
  igb1=KG0(1,ig,ik) 
  igb2=KG0(2,ig,ik) 
  igb3=KG0(3,ig,ik) 
  igb1=MOD(nrx2+igb1,nrx2)+1
  igb2=MOD(nry2+igb2,nry2)+1
  igb3=MOD(nrz2+igb3,nrz2)+1
  ind=igb1+(igb2-1)*nfft1+(igb3-1)*nfft1*nfft2 
  wfunc(ind)=dble(C0_K(ig,iw))
  wfunc(ind+Nl123)=dimag(C0_K(ig,iw))
 enddo!ig 
 call fft3_bw(fs,wfunc(1),fftwk(1)) 
 u_1D(:)=0.0d0 
 do ir=1,Nl123 
  u_1D(ir)=cmplx(wfunc(ir),wfunc(ir+Nl123)) 
 enddo!ir 
 do ir3=1,nrz2
  do ir2=1,nry2
   do ir1=1,nrx2
    ir=ir1+(ir2-1)*nfft1+(ir3-1)*nfft1*nfft2 
    u_3D(ir1,ir2,ir3,iw)=u_1D(ir) 
   enddo!ir1 
  enddo!ir2 
 enddo!ir3 
enddo!iw  
MAT_VHXC(:,:)=0.0d0 
do iw=1,NWF
 do jw=1,NWF
  SUM_CMPX=0.0D0 
  do ir3=1,nrz2
   do ir2=1,nry2
    do ir1=1,nrx2
     SUM_CMPX=SUM_CMPX+CONJG(u_3D(ir1,ir2,ir3,iw))*vhxc(ir1,ir2,ir3)*u_3D(ir1,ir2,ir3,jw)
    enddo!ir1 
   enddo!ir2 
  enddo!ir3 
  MAT_VHXC(iw,jw)=SUM_CMPX/dble(nrx2)/dble(nry2)/dble(nrz2)
 enddo!jw 
enddo!iw  
!write(6,*)'CHECK MAT_VHXC'
!do iw=1,NWF 
! write(6,*)(MAT_VHXC(iw,jw),jw=1,NWF)
!enddo 
!write(6,*)
deallocate(u_1D,u_3D)
RETURN  
END 
