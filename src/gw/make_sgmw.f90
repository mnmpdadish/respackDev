subroutine estimate_nsgm(ecmin,emin,emax,ecmax,gw_grid_separation,nproc,nsgm)
 ! 
 implicit none
 real(8)::ecmin,emin,emax,ecmax,gw_grid_separation 
 real(8)::omega,grd_separation 
 integer::ie,je,nproc
 integer::nsgm 
 integer,parameter::nsgm_max=100000 
 real(8),parameter::expansion=10.0d0!50.0d0  
 !
 je=0
 grd_separation=expansion*gw_grid_separation 
 omega=ecmin 
 do ie=1,nsgm_max
  !
  if(ecmin<=omega.and.omega<emin)then
   grd_separation=expansion*gw_grid_separation 
   je=je+1
  endif 
  !
  if(emin<=omega.and.omega<=emax)then
   grd_separation=gw_grid_separation 
   je=je+1
  endif 
  !
  if(emax<omega)then
   grd_separation=expansion*gw_grid_separation 
   je=je+1
  endif 
  !
  if(omega>ecmax.and.(mod(je,nproc).eq.0))exit 
  !
  omega=omega+grd_separation 
  !
 enddo!ie
 ! 
 !write(6000,*)'#je=',je  
 ! 
 nsgm=je 
 !
 return
end 
!
subroutine make_sgmw(ecmin,emin,emax,gw_grid_separation,nsgm,sgmw)
 !
 implicit none
 real(8)::ecmin,emin,emax,gw_grid_separation
 real(8)::omega,grd_separation 
 integer::nsgm 
 real(8)::sgmw(nsgm) 
 integer::ie 
 real(8),parameter::expansion=10.0d0!50.0d0  
 !
 sgmw=0.0d0 
 !
 grd_separation=expansion*gw_grid_separation 
 omega=ecmin 
 do ie=1,nsgm 
  !
  if(ecmin<=omega.and.omega<emin)then
   sgmw(ie)=omega
   grd_separation=expansion*gw_grid_separation 
  endif 
  !
  if(emin<=omega.and.omega<=emax)then
   sgmw(ie)=omega
   grd_separation=gw_grid_separation
  endif 
  !
  if(emax<omega)then
   sgmw(ie)=omega
   grd_separation=expansion*gw_grid_separation 
  endif 
  !
  omega=omega+grd_separation 
  !
 enddo!ie
 !
 !do ie=1,nsgm 
 ! write(6000,*) ie,sgmw(ie)!*au 
 !enddo 
 !
return
end 
