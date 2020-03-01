module m_fat_band
  implicit none
contains
  !
  subroutine calc_fat_band(NWF,NSK_BAND_DISP,DIST,EKS,VKS)
    implicit none 
    integer,intent(in)::NWF,NSK_BAND_DISP  
    real(8),intent(in)::DIST(NSK_BAND_DISP) 
    real(8),intent(in)::EKS(NWF,NSK_BAND_DISP) 
    complex(8),intent(in)::VKS(NWF,NWF,NSK_BAND_DISP)
    !
    character(99)::filename
    if(.true.)then 
      filename='./dir-wan/dat.iband.fat-'
      call wrt_fat_band(NWF,NSK_BAND_DISP,filename,DIST,EKS,VKS)
    endif 
    if(.true.)then 
      filename='./dir-wan/dat.iband.fat.pm3d-'
      call calc_fat_band_pm3d(NWF,NSK_BAND_DISP,filename,DIST,EKS,VKS)
    endif 
    return
  end subroutine calc_fat_band

  subroutine wrt_fat_band(NWF,NSK_BAND_DISP,filename,DIST,EKS,VKS)
    implicit none 
    integer,intent(in)::NWF,NSK_BAND_DISP  
    character(99),intent(in)::filename
    real(8),intent(in)::DIST(NSK_BAND_DISP) 
    real(8),intent(in)::EKS(NWF,NSK_BAND_DISP) 
    complex(8),intent(in)::VKS(NWF,NWF,NSK_BAND_DISP)
    !
    integer::iw,ib,ik 
    LOGICAL::REVERSE 
    character(99)::fname 
    real(8),parameter::au=27.21151d0 
    !
    !OPEN(400,W,file='./dir-wan/dat.iband.fat-ib')
    !
    do iw=1,NWF 
     write(fname,'(a,i3.3)') trim(filename),iw 
     OPEN(400,FILE=TRIM(fname)) 
     write(400,'(a)')'#Wannier-interpolaed fat band'
     write(400,'(a)')'#1:k, 2:Energy [eV], 3:|VKS(ib,iw,ik))|^2 ' 
     REVERSE=.TRUE.        
     do ib=1,NWF 
      if(REVERSE)then 
       do ik=1,NSK_BAND_DISP                     
        write(400,'(3f20.10)') DIST(ik)/DIST(NSK_BAND_DISP),EKS(ib,ik)*au,abs(VKS(ib,iw,ik))**2
       enddo!ik        
       REVERSE=.FALSE.        
      else         
       do ik=NSK_BAND_DISP,1,-1          
        write(400,'(3f20.10)') DIST(ik)/DIST(NSK_BAND_DISP),EKS(ib,ik)*au,abs(VKS(ib,iw,ik))**2
       enddo!ik        
       REVERSE=.TRUE.        
      endif!REVERSE                   
     enddo!ib 
     close(400) 
    enddo!iw 
    !
    return 
  end subroutine   
  
  subroutine calc_fat_band_pm3d(NWF,NSK_BAND_DISP,filename,DIST,EKS,VKS)
    implicit none 
    integer,intent(in)::NWF,NSK_BAND_DISP  
    character(99),intent(in)::filename
    real(8),intent(in)::DIST(NSK_BAND_DISP) 
    real(8),intent(in)::EKS(NWF,NSK_BAND_DISP) 
    complex(8),intent(in)::VKS(NWF,NWF,NSK_BAND_DISP)
    !
    real(8),parameter::au=27.21151d0
    real(8),parameter::delt=0.01d0/au!Greens function delt in au 
    real(8),parameter::delw=2.0d0*delt!Grid width in au 
    complex(8)::w,en,d 
    integer::iw,ik,jb,ie,fo 
    real(8),allocatable::AKW(:,:,:)!AKW(NSK_BAND_DISP,ndosgrd,NWF)           
    !
    real(8)::emax!=maxval(EIG)
    real(8)::emin!=minval(EIG)
    integer::ndosgrd!=int(1.2d0*diff/dlt)+1
    real(8),allocatable::dosgrd(:)!dosgrd(ndosgrd) 
    character(99)::fname 
    ! 
    !dos-grid
    !
    call est_ndosgrd(NWF,NSK_BAND_DISP,EKS(1,1),delw,emin,emax,ndosgrd)  
    allocate(dosgrd(ndosgrd));dosgrd=0.0d0 
    call make_dosgrd(emin,delw,ndosgrd,dosgrd(1))  
    !
    !calc spectranl function 
    !
    allocate(AKW(NSK_BAND_DISP,ndosgrd,NWF)); AKW=0.0d0
    do iw=1,NWF
     do ik=1,NSK_BAND_DISP 
      do jb=1,NWF 
       do ie=1,ndosgrd 
        w=cmplx(dosgrd(ie),-delt) 
        en=cmplx(EKS(jb,ik),0.0d0) 
        d=1.0d0/(w-en) 
        AKW(ik,ie,iw)=AKW(ik,ie,iw)+abs(imag(d))*abs(VKS(jb,iw,ik))**2  
       enddo!ie 
      enddo!jb 
     enddo!ik 
    enddo!iw  
    !
    !plt
    !
    !OPEN(400,W,file='./dir-wan/dat.iband.fat.pm3d-ib')
    !
    do iw=1,NWF 
     write(fname,'(a,i3.3)') trim(filename),iw 
     OPEN(400,FILE=TRIM(fname)) 
     write(400,'(a)')'#Wannier-interpolaed fat band for pm3d'
     write(400,'(a)')'#1:k, 2:Energy [eV], 3:AKW(ik,ie,iw)' 
     do ie=1,ndosgrd 
      do ik=1,NSK_BAND_DISP 
       if(AKW(ik,ie,iw)>=500.0d0) AKW(ik,ie,iw)=500.0d0 
       write(400,'(3f20.10)') DIST(ik)/DIST(NSK_BAND_DISP),dosgrd(ie)*au,AKW(ik,ie,iw) 
      enddo!ik 
      write(400,*) 
     enddo!ie 
     close(400)
    enddo!iw  
    !
    deallocate(dosgrd,AKW)
    !
    return
  end subroutine 
  ! 
  subroutine est_ndosgrd(NTB,NTK,E_EIG,deltw,emin,emax,ndosgrd)  
    implicit none
    integer,intent(in)::NTB,NTK
    real(8),intent(in)::E_EIG(NTB,NTK)
    real(8),intent(in)::deltw 
    real(8),intent(out)::emax,emin
    integer,intent(out)::ndosgrd
    real(8)::diff 
    real(8),parameter::au=27.21151d0
    !
    emax=maxval(E_EIG)
    emin=minval(E_EIG)
    diff=emax-emin 
    !
    !define grid range
    !
    ndosgrd=int(1.2d0*diff/deltw)+1
    emax=emax+0.10d0*diff
    emin=emin-0.10d0*diff
    !
    write(6,*)
    write(6,'(a40)')'+++ m_fat_band: est_ndosgrd +++'
    write(6,'(a40)')'GRID DATA FOR FAT BAND'
    write(6,'(a40,f15.8)')'emin(eV)',emin*au 
    write(6,'(a40,f15.8)')'emax(eV)',emax*au 
    write(6,'(a40,f15.8)')'deltw(eV)',deltw*au 
    write(6,'(a40,i15)')'ndosgrd',ndosgrd 
    write(6,*)
    return
  end subroutine
  !
  subroutine make_dosgrd(emin,deltw,ndosgrd,dosgrd)  
    implicit none
    real(8),intent(in)::emin,deltw 
    integer,intent(in)::ndosgrd
    real(8),intent(out)::dosgrd(ndosgrd) 
    integer::ie 
    !
    dosgrd=0.0d0 
    do ie=1,ndosgrd
     dosgrd(ie)=emin+dble(ie)*deltw 
    enddo 
    return
  end subroutine 
  !
end module m_fat_band 
