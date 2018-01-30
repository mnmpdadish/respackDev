subroutine makekpts(Ndiv,Nblk,NSK_BAND_DISP,SKI,SK)  
implicit none 
integer::Nblk,Ndiv,NSK_BAND_DISP 
real(8)::SKI(3,Nblk) 
real(8)::SK(3,NSK_BAND_DISP)
integer::i,k  
do k=1,Nblk-1 
do i=1,Ndiv 
 SK(:,i+Ndiv*(k-1))& 
=SKI(:,k)+(SKI(:,k+1)-SKI(:,k))/dble(Ndiv)*dble(i-1)  
enddo 
enddo  
SK(:,Ndiv*(Nblk-1)+1)=SKI(:,Nblk) 
!write(6,*) NSK_BAND_DISP
!do i=1,Ndiv*(Nblk-1)+1
!write(6,*) SK(:,i)
!enddo 
return 
end 
