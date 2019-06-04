      PROGRAM WANNIER          
      implicit none

      integer::ik,ib,ic,ig, ncomp, NBAND, Nk_irr, NMAX,NN,NTG
      complex(8),allocatable::CIR(:,:,:,:),value(:)
      integer,allocatable::NGI(:)

      !input
      character(100) :: inputfile
      logical fileExists
      !/input

      if(iargc()/=1) then
        write(6,*) "hey, must have one argument."
        write(6,*) 
        write(6,*) "example of use:"
        write(6,*) "/a.out work/prefix.save/"
        stop
      endif

      call getarg(1, inputfile)
      inquire( file = (TRIM(inputfile)//"/dat.sample-k"), 
     +                                     EXIST=fileExists ) 
      if(.not.fileExists) then
        write(6,*) "error: file does not exist."
        write(6,*) "file: "//(TRIM(inputfile)//"/dat.sample-k")
        write(6,*) "terminated."
        stop 
      endif 
      
      OPEN(101,FILE=(TRIM(inputfile)//"/dat.sample-k")) 
      rewind(101) 
      read(101,*) Nk_irr 
      close(101) 


      OPEN(132,FILE=TRIM(inputfile)//'/dat.nkm')
      allocate(NGI(Nk_irr));NGI(:)=0
      rewind(132)
      do ik=1,Nk_irr 
       read(132,*) NGI(ik) 
      enddo 
      close(132) 
      NTG=maxval(abs(NGI(:))) 
      write(6,'(a10,i10)')'NTG=',NTG  
            
      OPEN(111,FILE=TRIM(inputfile)//'/dat.eigenvalue')
      rewind(111)
      read(111,*) NBAND 
      close(111)
            
      NMAX=20
      NN=0
      OPEN(102,FILE=(TRIM(inputfile)//'/dat.wfn'),FORM='unformatted') 
      rewind(102)
      read(102) ncomp 
      write(6,*) 'NTG,ncomp,NTB,Nk_irr ', NTG, ncomp, NBAND, Nk_irr
      allocate(CIR(NTG,ncomp,NBAND,Nk_irr));CIR=0.0d0  
      allocate(value(1));value=0.0d0  
      do ik=1,Nk_irr
       write(6,*) 'NGI(ik) ', NGI(ik)
       do ib=1,NBAND
        !do ic=1,ncomp
         
      !do ik=1,2
       !write(6,*) 'NGI(ik) ', NGI(ik)
       !do ib=1,2
        !do ic=1,ncomp
        read(102)((CIR(IG,ic,ib,ik),IG=1,NGI(ik)),ic=1,ncomp)
        !enddo!ic 
        do ic=1,ncomp

         !write(6,*) 'CIR ', ik,ib,ncomp, NGI(ik)
         !do ig=1,NGI(ik)
          !read(102) value(1)
          !write(6,'(2F10.5)') value(1)
          
          
          
          !read(102)(CIR(IG,ic,ib,ik),IG=1,NGI(ik))
          write(6,*) 'CIR,   ik,  ib,  ic,  NGI(ik)'
          write(6,*) 'CIR ', ik,ib,ic, NGI(ik)
          write(6,'(2F10.5)') (CIR(IG,ic,ib,ik),IG=1,NGI(ik))
         !enddo!ig 
        enddo!ic 
       enddo!ib
      enddo!ik
      close(102)
      
      END
