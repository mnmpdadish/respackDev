      SUBROUTINE CALC_GRADIENT(NTK,NTB,n_occ,NB,WF_CHARGE,WEIGHT_b,
     + OVERLAP,tci,GRADIENT,TR_GRAD)                  
      implicit none           
      integer::NTK,NTB,n_occ,NB         
      integer::ik,ibvec,i_band,j_band
      real(8)::WF_CHARGE(n_occ,NTK,NB),WEIGHT_b(NB)                 
      real(8)::TR_GRAD(NTK)              
      real(8)::SUM_REAL                 
      complex(8)::OVERLAP(NTB,NTB,NTK,NB)         
      complex(8)::GRADIENT(n_occ,n_occ,NTK)              
      complex(8)::WF_R(n_occ,n_occ,NB)         
      complex(8)::WF_T(n_occ,n_occ,NB)         
      complex(8)::SUM_CMPX,tci       
!--
      GRADIENT(:,:,:)=0.0D0        
      Do ik=1,NTK        
       Do i_band=1,n_occ       
        Do j_band=1,n_occ       
         Do ibvec=1,NB           
          WF_R(i_band,j_band,ibvec)          
     +   =OVERLAP(i_band,j_band,ik,ibvec)
     +   *CONJG(OVERLAP(j_band,j_band,ik,ibvec)) !modified 2007 6 26 
          WF_T(i_band,j_band,ibvec)          
     +   =OVERLAP(i_band,j_band,ik,ibvec)
     +   /OVERLAP(j_band,j_band,ik,ibvec) !modified 2007 6 26 
     +   *WF_CHARGE(j_band,ik,ibvec) !modified 2007 6 26 
         ENDDO !ibvec            
        ENDDO !j_band           
       ENDDO !i_band           
       Do i_band=1,n_occ       
        Do j_band=1,n_occ       
         SUM_CMPX=0.0D0 
         Do ibvec=1,NB           
          SUM_CMPX 
     +   =SUM_CMPX 
     +   +4.0D0 
     +   *WEIGHT_b(ibvec) 
     +   *(((WF_R(i_band,j_band,ibvec)          
     +    -CONJG(WF_R(j_band,i_band,ibvec)))/2.0D0)            
     +    -((WF_T(i_band,j_band,ibvec)          
     +    +CONJG(WF_T(j_band,i_band,ibvec)))/tci))           
         ENDDO !ibvec            
         GRADIENT(i_band,j_band,ik)=SUM_CMPX/CMPLX(NTK) 
        ENDDO !j_band           
       ENDDO !i_band           
      ENDDO !ik       
!-- block diagonal constraint --
!      G(:,:,:)=GRADIENT(:,:,:) 
!      GRADIENT=0.0d0 
!      GRADIENT(1:2,1:2,:)=G(1:2,1:2,:) 
!      GRADIENT(3:4,3:4,:)=G(3:4,3:4,:) 
!      GRADIENT(5:6,5:6,:)=G(5:6,5:6,:) 
!-- Tr(G^2) --                 
      TR_GRAD(:)=0.0D0      
      SUM_REAL=0.0D0                          
      Do ik=1,NTK    
       SUM_CMPX=0.0D0        
       Do i_band=1,n_occ 
        Do j_band=1,n_occ 
         SUM_CMPX 
     +  =SUM_CMPX 
     +  +GRADIENT(i_band,j_band,ik)     
     +  *GRADIENT(j_band,i_band,ik)     
        ENDDO!j_band   
       ENDDO!i_band   
       TR_GRAD(ik)=DBLE(SUM_CMPX)          
       SUM_REAL=SUM_REAL+TR_GRAD(ik)          
      ENDDO !ik         
      write(6,*) 'Tr(G^2)=',SUM_REAL          
      RETURN
      END                
!--
      SUBROUTINE CALC_OMEGA(NTK,NTB,n_occ,NB,VEC_b,ILS,       
     + WEIGHT_b,OVERLAP,WF_CHARGE,OMEGA_I,OMEGA_OD,OMEGA_D)        
      implicit none 
      integer::NTK,NTB,n_occ,NB,ILS                 
      integer::ik,ibvec,i_band,j_band
      real(8)::VEC_b(3,NB),WEIGHT_b(NB)                 
      real(8)::WF_CHARGE(n_occ,NTK,NB)          
      real(8)::OMEGA_I,OMEGA_OD,OMEGA_D                      
      real(8)::WF_POSITION(3,n_occ)         
      real(8)::VEC_SUM_REAL(3),SUM_REAL,SUM_REAL1,SUM_REAL2 
      complex(8)::LOG_OVERLAP      
      complex(8)::OVERLAP(NTB,NTB,NTK,NB)         

!(1) WANNIER POSITION

      WF_POSITION(:,:) = 0.0D0            
      Do i_band = 1 , n_occ       
      VEC_SUM_REAL(:) = 0.0D0          
      Do ik     = 1 , NTK        
      Do ibvec  = 1 , NB           
       LOG_OVERLAP=LOG(OVERLAP(i_band,i_band,ik,ibvec)) 
         VEC_SUM_REAL(:) 
     + = VEC_SUM_REAL(:) 
     + + WEIGHT_b(ibvec) 
     + * VEC_b (:,ibvec) 
     + * AIMAG( LOG_OVERLAP ) 
      ENDDO!ibvec            
      ENDDO!ik       
      WF_POSITION(:,i_band)=-VEC_SUM_REAL(:)/DBLE(NTK)   
      ENDDO!i_band           

!(2) WANNIER CHARGE

      WF_CHARGE(:,:,:) = 0.0D0            
      Do ik     = 1 , NTK        
      Do ibvec  = 1 , NB           
      Do i_band = 1 , n_occ       
       SUM_REAL = VEC_b(1,ibvec)*WF_POSITION(1,i_band)         
     +          + VEC_b(2,ibvec)*WF_POSITION(2,i_band)         
     +          + VEC_b(3,ibvec)*WF_POSITION(3,i_band)         
       LOG_OVERLAP = LOG ( OVERLAP(i_band,i_band,ik,ibvec) ) 
         WF_CHARGE(i_band,ik,ibvec) 
     + = AIMAG( LOG_OVERLAP )
     + + SUM_REAL       
      ENDDO!i_band           
      ENDDO!ibvec            
      ENDDO!ik       

!(3-a) OMEGA_I 

      SUM_REAL1 = 0.0D0            
      Do ik     = 1 , NTK        
      Do ibvec  = 1 , NB           
      SUM_REAL2 = 0.0D0          
      Do i_band = 1 , n_occ       
      Do j_band = 1 , n_occ       
         SUM_REAL2 
     + = SUM_REAL2 
     + + ( ABS( OVERLAP(i_band,j_band,ik,ibvec) ) )**2     
      ENDDO!j_band           
      ENDDO!i_band           
         SUM_REAL1 
     + = SUM_REAL1 
     + + WEIGHT_b(ibvec) 
     + * ( DBLE(n_occ) - SUM_REAL2 )       
      ENDDO!ibvec            
      ENDDO!ik       
      OMEGA_I = SUM_REAL1 / DBLE(NTK)        

!(3-b) OMEGA_OD 

      SUM_REAL1 = 0.0D0            
      Do ik     = 1 , NTK        
      Do ibvec  = 1 , NB           
      SUM_REAL2 = 0.0D0          
      Do i_band = 1 , n_occ       
      Do j_band = 1 , n_occ       
      IF(i_band /= j_band) THEN          
         SUM_REAL2 
     + = SUM_REAL2 
     + + ( ABS( OVERLAP(i_band,j_band,ik,ibvec) ) )**2     
      ENDIF            
      ENDDO!j_band           
      ENDDO!i_band           
         SUM_REAL1 
     + = SUM_REAL1 
     + + WEIGHT_b(ibvec) 
     + * SUM_REAL2                  
      ENDDO!ibvec            
      ENDDO!ik       
      OMEGA_OD = SUM_REAL1 / DBLE(NTK)        

!(3-c) OMEGA_D 

      SUM_REAL1 = 0.0D0            
      Do ik     = 1 , NTK        
      Do ibvec  = 1 , NB           
      SUM_REAL2 = 0.0D0          
      Do i_band = 1 , n_occ       
         SUM_REAL2 
     + = SUM_REAL2 
     + + ( WF_CHARGE(i_band,ik,ibvec) )**2     
      ENDDO!i_band           
         SUM_REAL1 
     + = SUM_REAL1 
     + + WEIGHT_b(ibvec) 
     + * SUM_REAL2                  
      ENDDO!ibvec            
      ENDDO!ik       
      OMEGA_D = SUM_REAL1 / DBLE(NTK)        

!(4) OMEGA_I,OD AND OMEGA_D OF EACH BAND 

      IF(ILS == 1) THEN 
       Do i_band = 1 , n_occ       
        SUM_REAL1 = 0.0D0            
        SUM_REAL2 = 0.0D0            
        Do ik     = 1 , NTK        
         Do ibvec  = 1 , NB           
          SUM_REAL1 
     +   =SUM_REAL1 
     +   +WEIGHT_b(ibvec) 
     +   *(1.0D0-(ABS(OVERLAP(i_band,i_band,ik,ibvec)))**2)       
          SUM_REAL2 
     +   =SUM_REAL2 
     +   +WEIGHT_b(ibvec)                  
     +   *(WF_CHARGE(i_band,ik,ibvec))**2     
         ENDDO!ibvec            
        ENDDO!ik       
        SUM_REAL1=SUM_REAL1/DBLE(NTK) 
        SUM_REAL2=SUM_REAL2/DBLE(NTK) 
        write(6,*) i_band,SUM_REAL1,SUM_REAL2 
       ENDDO!i_band           
      ENDIF               
!--
      RETURN
      END
!--
      SUBROUTINE UNITARY_GEN(NTK,n_occ,ILS, 
     + DIRECTION,DIRECTION2,STEP_LENGTH,UNITARY,U_ORG,U_OLD)          
      implicit none           
      integer::NTK,n_occ,ILS          
      integer::ik,i_band,j_band,k_band   
      real(8)::STEP_LENGTH         
      real(8)::D_EIG(n_occ),D_INV(n_occ),D_COS(n_occ),D_SIN(n_occ) 
      complex(8),INTENT(IN)::DIRECTION (n_occ,n_occ,NTK)              
      complex(8),INTENT(IN)::DIRECTION2(n_occ,n_occ,NTK)              
      complex(8)::UNITARY(n_occ,n_occ,NTK)              
      complex(8)::U_ORG(n_occ,n_occ,NTK)              
      complex(8)::U_OLD(n_occ,n_occ,NTK)              
      complex(8)::D_MAT(n_occ,n_occ) 
      complex(8)::R_MAT(n_occ,n_occ) 
      complex(8)::RC_MAT(n_occ,n_occ) 
      complex(8)::RS_MAT(n_occ,n_occ) 
      complex(8)::G_MAT(n_occ,n_occ) 
      complex(8)::SUM_CMPX

      UNITARY(:,:,:)=0.0D0 
      IF(ILS==1) THEN 
       Do ik=1,NTK           
!STEP(1) diag D_MAT
        D_MAT(:,:)=DIRECTION2(:,:,ik)              
        D_EIG(:)=0.0D0                
        call diagV(n_occ,D_MAT(1,1),D_EIG(1)) 
        R_MAT(:,:)=D_MAT(:,:)

!STEP(2) CONSTRUCT DIAGONAL MATRICES
        D_INV(:)=0.0D0           
        D_COS(:)=0.0D0 
        D_SIN(:)=0.0D0   
        Do i_band=1,n_occ         
         IF(D_EIG(i_band)==0.0D0) THEN 
          D_EIG(i_band)=1.0D-20 
         ENDIF          
         D_EIG(i_band)=DSQRT(ABS(D_EIG(i_band)))          
         D_INV(i_band)=1.0D0/D_EIG(i_band)        
         D_COS(i_band)=DCOS(D_EIG(i_band)*STEP_LENGTH)       
         D_SIN(i_band)=DSIN(D_EIG(i_band)*STEP_LENGTH)       
        ENDDO 

!STEP(3) CONSTRUCT RC & RS MATRIX 
        RC_MAT(:,:)=0.0D0    
        RS_MAT(:,:)=0.0D0    
        Do i_band=1,n_occ   
         Do j_band=1,n_occ   
           RC_MAT(i_band,j_band) 
     +    =R_MAT(i_band,j_band) 
     +    *D_COS(j_band)       
           RS_MAT(i_band,j_band) 
     +    =R_MAT(i_band,j_band) 
     +    *D_SIN(j_band) 
     +    *D_INV(j_band)           
         ENDDO 
        ENDDO 

!STEP(4) CONSTRUCT RC MATRIX 
        G_MAT(:,:)=DIRECTION(:,:,ik)         
        Do i_band=1,n_occ   
         Do j_band=1,n_occ   
          SUM_CMPX=0.0D0            
          Do k_band=1,n_occ       
           SUM_CMPX 
     +    =SUM_CMPX 
     +    +G_MAT(i_band,k_band) 
     +    *RS_MAT(k_band,j_band)     
          ENDDO 
          RC_MAT(i_band,j_band) 
     +   =RC_MAT(i_band,j_band)+SUM_CMPX        
         ENDDO 
        ENDDO 

!STEP(5) CONSTRUCT U MATRIX 
        Do i_band=1,n_occ   
         Do j_band=1,n_occ   
          SUM_CMPX=0.0D0            
          Do k_band=1,n_occ       
           SUM_CMPX 
     +    =SUM_CMPX 
     +    +RC_MAT(i_band,k_band) 
     +    *CONJG(R_MAT(j_band,k_band))        
          ENDDO 
          UNITARY(i_band,j_band,ik)=SUM_CMPX        
         ENDDO 
        ENDDO 
       ENDDO!ik === UNITARY MATRIX GENERATION ROOP ===     
       U_ORG(:,:,:)=UNITARY(:,:,:)
       U_OLD(:,:,:)=UNITARY(:,:,:)
      ELSE 
       Do ik=1,NTK 
        DO i_band=1,n_occ 
         DO j_band=1,n_occ 
          SUM_CMPX=0.0D0            
          Do k_band=1,n_occ       
           SUM_CMPX 
     +    =SUM_CMPX 
     +    +U_OLD(i_band,k_band,ik) 
     +    *U_ORG(k_band,j_band,ik)  
          ENDDO 
          UNITARY(i_band,j_band,ik)=SUM_CMPX        
         ENDDO!j_band 
        ENDDO!i_band 
       ENDDO!ik 
       U_OLD(:,:,:)=UNITARY(:,:,:)
      ENDIF
      RETURN
      END        
!--
      SUBROUTINE UPDATE_OVERLAP(NTK,NTB,n_occ,NB,KPT,UNITARY,OVERLAP)  
      implicit none           
      integer::NTK,NTB,n_occ,NB
      integer::KPT(NTK,NB)         
      integer::ik,ik_ib,ibvec,i_band,j_band,k_band   
      complex(8)::UNITARY(n_occ,n_occ,NTK)              
      complex(8)::OVERLAP(NTB,NTB,NTK,NB)         
      complex(8)::U_MAT(n_occ,n_occ) 
      complex(8)::SUM_CMPX

      Do ik=1,NTK 
       Do ibvec=1,NB      
        ik_ib=KPT(ik,ibvec)        
        Do i_band=1,n_occ          
         Do j_band=1,n_occ          
          SUM_CMPX=0.0D0          
          Do k_band=1,n_occ 
           SUM_CMPX 
     +    =SUM_CMPX 
     +    +OVERLAP(i_band,k_band,ik,ibvec)
     +    *UNITARY(k_band,j_band,ik_ib)
          ENDDO 
          U_MAT(i_band,j_band)=SUM_CMPX         
         ENDDO 
        ENDDO 
        Do i_band=1,n_occ          
         Do j_band=1,n_occ          
          SUM_CMPX=0.0D0          
          Do k_band=1,n_occ 
           SUM_CMPX 
     +    =SUM_CMPX 
     +    +CONJG(UNITARY(k_band,i_band,ik))        
     +    *U_MAT(k_band,j_band)
          ENDDO 
          OVERLAP(i_band,j_band,ik,ibvec)=SUM_CMPX         
         ENDDO 
        ENDDO 
       ENDDO !ibvec                   
      ENDDO !ik          
      RETURN
      END
!---
      SUBROUTINE A_MAT_UPDATE(NTK,NTB,NGAUSS,n_occ,UNITARY,A_MAT)   
      implicit none           
      integer::NTK,NTB,NGAUSS,n_occ
      integer::ik,i_band,j_band,k_band   
      complex(8)::UNITARY(n_occ,n_occ,NTK)              
      complex(8)::A_MAT(NTB,NGAUSS,NTK)           
      complex(8)::A_TMP(NTB,NGAUSS)        
      complex(8)::SUM_CMPX

      Do ik=1,NTK        
       A_TMP(:,:)=0.0D0             
       A_TMP(:,:)=A_MAT(:,:,ik)            
       Do i_band=1,n_occ       
        Do j_band=1,n_occ       
         SUM_CMPX=0.0D0 
         Do k_band=1,n_occ           
          SUM_CMPX 
     +   =SUM_CMPX 
     +   +A_TMP(i_band,k_band) 
     +   *UNITARY(k_band,j_band,ik)        
         ENDDO              
         A_MAT(i_band,j_band,ik)=SUM_CMPX     
        ENDDO !j_band
       ENDDO !i_band
      ENDDO !ik          
      RETURN
      END    
!---
      SUBROUTINE CALC_OVERLAP(NTK,NTB,NTG,nvx,nvy,nvz,ik,ik_ib,NG0,
     + KG0,N_BAND,SHIFT_b,C0_BRA,C0_KET,OVERLAP_TMP,N_BAND_BOTTOM)
      implicit none           
      integer::NTK,NTB,NTG,nvx,nvy,nvz       
      integer::NG0(NTK),KG0(3,NTG,NTK),N_BAND(NTK)         
      integer::SHIFT_b(3),N_BAND_BOTTOM(NTK)         
      integer::PULLBACK(NTG)           
      integer::MAP(-nvx:nvx,-nvy:nvy,-nvz:nvz)
      integer::ik,ik_ib,ig,ig1,ig2,i_band,j_band
      complex(8)::C0_BRA(0:NTG,NTB)    
      complex(8)::C0_KET(0:NTG,NTB)    
      complex(8)::OVERLAP_TMP(NTB,NTB)            
      complex(8)::ELEMENT_OVERLAP
!MAP FUNCTION
      MAP(:,:,:)=0
      do ig=1,NG0(ik) 
       MAP(KG0(1,ig,ik)+SHIFT_b(1),
     +     KG0(2,ig,ik)+SHIFT_b(2),
     +     KG0(3,ig,ik)+SHIFT_b(3))=ig                  
      enddo       
!PULLBACK 
      PULLBACK(:)=0
      do ig=1,NG0(ik_ib) 
         PULLBACK(ig) 
     + = MAP(KG0(1,ig,ik_ib) 
     +      ,KG0(2,ig,ik_ib) 
     +      ,KG0(3,ig,ik_ib))     
      enddo       
!INTERSTATE MATRIX 
      OVERLAP_TMP(:,:)=0.0D0         
      do i_band=1,N_BAND(ik)        
       do j_band=1,N_BAND(ik_ib)       
        ELEMENT_OVERLAP=0.0D0 
        do ig2=1,NG0(ik_ib)     
         ig1=PULLBACK(ig2)   
          ELEMENT_OVERLAP 
     +  = ELEMENT_OVERLAP 
     +  + CONJG(C0_BRA(ig1,i_band+N_BAND_BOTTOM(ik)))
     +  *       C0_KET(ig2,j_band+N_BAND_BOTTOM(ik_ib))
        enddo!IG2 
        OVERLAP_TMP(i_band,j_band)=ELEMENT_OVERLAP      
       enddo!j_band     
      enddo!i_band     
      RETURN 
      END              
!--
      subroutine kcheck(ktmp,RWtmp) 
      implicit none 
      real(8),intent(inout)::ktmp(3)
      integer,intent(out)::RWtmp(3) 
      real(8),parameter::dlt_BZ=1.0d-6 
!--
!      write(6,'(a8,3f15.10)') 'before',ktmp(1),ktmp(2),ktmp(3)
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
!--
!      write(6,'(a8,3f15.10)') 'after',ktmp(1),ktmp(2),ktmp(3)
!--
      return
      end 
!--
!20170316 
      subroutine kcheck_trs(ktmp,RWtmp) 
      implicit none 
      real(8),intent(inout)::ktmp(3)
      integer,intent(out)::RWtmp(3) 
      real(8),parameter::dlt_BZ=-1.0d-6 
!--
!      write(6,'(a8,3f15.10)') 'before',ktmp(1),ktmp(2),ktmp(3)
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
!--
!      write(6,'(a8,3f15.10)') 'after',ktmp(1),ktmp(2),ktmp(3)
!--
      return
      end 
!--
      subroutine kcheck_old(ktmp,RWtmp) 
      implicit none 
      real(8),intent(inout)::ktmp(3)
      integer,intent(out)::RWtmp(3) 
      if(ktmp(1)>=0.5d0)then 
       ktmp(1)=ktmp(1)-1.0d0;RWtmp(1)=-1
      endif 
      if(ktmp(1)<-0.5d0)then 
       ktmp(1)=ktmp(1)+1.0d0;RWtmp(1)=1
      endif 
      if(ktmp(2)>=0.5d0)then 
       ktmp(2)=ktmp(2)-1.0d0;RWtmp(2)=-1
      endif 
      if(ktmp(2)<-0.5d0)then 
       ktmp(2)=ktmp(2)+1.0d0;RWtmp(2)=1
      endif 
      if(ktmp(3)>=0.5d0)then 
       ktmp(3)=ktmp(3)-1.0d0;RWtmp(3)=-1
      endif 
      if(ktmp(3)<-0.5d0)then 
       ktmp(3)=ktmp(3)+1.0d0;RWtmp(3)=1
      endif 
      return
      end 
!---
      subroutine make_KG0(NTG,b1,b2,b3,Gcut,q1,q2,q3,KG0,NG) 
      implicit none 
      integer,intent(in)::NTG
      real(8),intent(in)::b1(3),b2(3),b3(3) 
      real(8),intent(in)::Gcut,q1,q2,q3        
      integer,intent(out)::KG0(3,NTG),NG 
      integer::igL,igL1,igL2,igL3
      real(8)::qgL(3),qgL2  
      integer,parameter::NGL1=150!100!300
      integer,parameter::NGL2=150!100!300 
      integer,parameter::NGL3=150!100!300

      igL=0
      do igL1=-NGL1,NGL1 
      do igL2=-NGL2,NGL2 
      do igL3=-NGL3,NGL3 
      qgL(:)=(q1+dble(igL1))*b1(:)+(q2+dble(igL2))*b2(:)
     +      +(q3+dble(igL3))*b3(:)    
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
!write(6,*) "maxabs KG0_1=",maxval(abs(KG0(1,:))) 
!write(6,*) "maxabs KG0_2=",maxval(abs(KG0(2,:))) 
!write(6,*) "maxabs KG0_3=",maxval(abs(KG0(3,:))) 
!--
!      do ig=1,NG 
!       write(6,'(4i8)') ig,KG0(:,ig)
!      enddo 
!--
      RETURN 
      END 
!--
      subroutine make_C0(NTG,NTB,itrs,NG,KGtmp,RWtmp,rginvtmp,pgtmp,
     +           nnp,L1,L2,L3,packtmp,C0irr,C0_K) 
      implicit none 

      integer::NTG,NTB,itrs,NG,L1,L2,L3,nnp 
      integer::KGtmp(3,NTG),RWtmp(3),pgtmp(3) 
      integer::packtmp(-L1:L1,-L2:L2,-L3:L3) 
      integer::ig,jg,i1,i2,i3,j1,j2,j3,k1,k2,k3 
      real(8)::rginvtmp(3,3),phase 
      complex(8)::C0irr(NTG,NTB),C0_K(NTG,NTB) 
      complex(8)::pf
      real(8),parameter::pi=dacos(-1.0d0)
      real(8),parameter::tpi=2.0d0*pi 
      complex(8),parameter::ci=(0.0D0,1.0D0) 
!--
      C0_K(:,:)=0.0d0 
      select case(itrs) 
      case(1)!=== not time-reversal ===*      
      do ig=1,NG 
       i1=KGtmp(1,ig); j1=KGtmp(2,ig); k1=KGtmp(3,ig) 
       i2=i1+RWtmp(1); j2=j1+RWtmp(2); k2=k1+RWtmp(3) 
       i3=int(rginvtmp(1,1))*i2+int(rginvtmp(1,2))*j2 
     +   +int(rginvtmp(1,3))*k2 
       j3=int(rginvtmp(2,1))*i2+int(rginvtmp(2,2))*j2 
     +   +int(rginvtmp(2,3))*k2 
       k3=int(rginvtmp(3,1))*i2+int(rginvtmp(3,2))*j2 
     +   +int(rginvtmp(3,3))*k2 
       jg=packtmp(i3,j3,k3) 
       phase=tpi*(dble(i1)*dble(pgtmp(1))
     +           +dble(j1)*dble(pgtmp(2))
     +           +dble(k1)*dble(pgtmp(3))) 
       pf=exp(-ci*phase/dble(nnp)) 
       C0_K(ig,:)=C0irr(jg,:)*pf 
      enddo!ig 
      case(-1)!=== time-reversal ===*      
      do ig=1,NG 
       i1=KGtmp(1,ig); j1=KGtmp(2,ig); k1=KGtmp(3,ig) 
       i2=i1+RWtmp(1); j2=j1+RWtmp(2); k2=k1+RWtmp(3) 
       i3=int(rginvtmp(1,1))*i2+int(rginvtmp(1,2))*j2 
     +   +int(rginvtmp(1,3))*k2 
       j3=int(rginvtmp(2,1))*i2+int(rginvtmp(2,2))*j2 
     +   +int(rginvtmp(2,3))*k2 
       k3=int(rginvtmp(3,1))*i2+int(rginvtmp(3,2))*j2 
     +   +int(rginvtmp(3,3))*k2 
       jg=packtmp(i3,j3,k3) 
       phase=tpi*(dble(i1)*dble(pgtmp(1))
     +           +dble(j1)*dble(pgtmp(2))
     +           +dble(k1)*dble(pgtmp(3))) 
       pf=exp(-ci*phase/dble(nnp)) 
       C0_K(ig,:)=C0irr(jg,:)*pf 
      enddo!ig 
      C0_K(:,:)=conjg(C0_K(:,:))  
      end select 
!--
!      SUM_CMPX=0.0d0 
!      do ig=1,NG 
!       SUM_CMPX=SUM_CMPX+CONJG(C0_K(ig))*C0_K(ig) 
!      enddo 
!      write(6,'(a,x,2f15.10)') 'sumC0_K',SUM_CMPX 
!--
      return 
      end 
!--
      SUBROUTINE OUTER_PRODUCT(vec_x,vec_y,vec_z)
      implicit none 
      real(8)::vec_x(3),vec_y(3),vec_z(3) 
      vec_z(1)=vec_x(2)*vec_y(3)-vec_x(3)*vec_y(2)
      vec_z(2)=vec_x(3)*vec_y(1)-vec_x(1)*vec_y(3) 
      vec_z(3)=vec_x(1)*vec_y(2)-vec_x(2)*vec_y(1)
      RETURN
      END
!--
      subroutine BZcheck(ktmp,kflg) 
      implicit none 
      real(8),intent(in)::ktmp(3)
      real(8)::kx,ky,kz 
      integer,intent(out)::kflg 
      kx=-ktmp(1) 
      ky=-ktmp(2) 
      kz=-ktmp(3) 
      kflg=0 
      if(kx<=-0.5d0) kflg=1
      if(kx>0.5d0) kflg=1
      if(ky<=-0.5d0) kflg=1
      if(ky>0.5d0) kflg=1
      if(kz<=-0.5d0) kflg=1
      if(kz>0.5d0) kflg=1
      return
      end
!--
!set NB_start,NB_end 
      subroutine calc_NB_start_end(NK,NB,EIG,EL,EU,NB_start,NB_end) 
      implicit none
      integer::NK,NB
      real(8)::EIG(NB,NK) 
      real(8)::EL,EU  
      integer::NB_start,NB_end
      integer::Ein(NB,NK) 
      integer::ik,ib
!
      Ein=0 
      do ik=1,NK 
       do ib=1,NB 
        if(EL<=EIG(ib,ik).and.EIG(ib,ik)<=EU)then 
         Ein(ib,ik)=ib
        endif 
       enddo
      enddo   
      NB_end=maxval(Ein)
      do ik=1,NK 
       do ib=1,NB 
        if(Ein(ib,ik)==0) Ein(ib,ik)=10000
       enddo
      enddo 
      NB_start=minval(Ein) 
      return
      end 
!--
!SHIFT VECTOR
      subroutine make_shift_vector(NTK,NB,ik,ibvec,SK0,b_LAT,
     +SHIFT_b,ik_ib)
      implicit none
      integer::NTK,NB,ik,ibvec 
      real(8)::SK0(3,NTK) 
      real(8)::b_LAT(3,NB)  
      integer::SHIFT_b(3)
      integer::ik_ib 
      real(8)::SL(3) 
      integer::jk 
      real(8),parameter::dlt_BZ=1.0d-6 
!--
      SL(:)=SK0(:,ik)+b_LAT(:,ibvec)         
!     write(6,'(a8,3f15.10)')'SLbfr=',SL(1),SL(2),SL(3) 
!     dlt_BZ=1.d-6
      SHIFT_b(:)=0           
      IF(SL(1)<=-0.5D0+dlt_BZ)THEN!20170210 
       SL(1)=SL(1)+1.0D0 
       SHIFT_b(1)=-1                
      ENDIF           
      IF(SL(2)<=-0.5D0+dlt_BZ)THEN!20170210
       SL(2)=SL(2)+1.0D0 
       SHIFT_b(2)=-1                
      ENDIF           
      IF(SL(3)<=-0.5D0+dlt_BZ)THEN!20170210 
       SL(3)=SL(3)+1.0D0 
       SHIFT_b(3)=-1                
      ENDIF           
      IF(SL(1)>0.5D0+dlt_BZ)THEN
       SL(1)=SL(1)-1.0D0 
       SHIFT_b(1)=1                
      ENDIF           
      IF(SL(2)>0.5D0+dlt_BZ)THEN
       SL(2)=SL(2)-1.0D0 
       SHIFT_b(2)=1                
      ENDIF           
      IF(SL(3)>0.5D0+dlt_BZ)THEN
       SL(3)=SL(3)-1.0D0 
       SHIFT_b(3)=1                
      ENDIF           
!     write(6,'(a8,3f15.10)')'SLaft=',SL(1),SL(2),SL(3) 
!IK_IB
      do jk=1,NTK          
       IF(ABS(SK0(1,jk)-SL(1))<=1.d-4.AND.
     +    ABS(SK0(2,jk)-SL(2))<=1.d-4.AND.
     +    ABS(SK0(3,jk)-SL(3))<=1.d-4)THEN           
        ik_ib=jk           
       ENDIF        
      enddo!jk               
      return
      end 
!--
!       call make_amat(nigs,NTB,NTG,NG0(ik),SK0(1,ik),KG0(1,1,ik),
!     + b1(1),b2(1),b3(1),TAU_GAUSS(1,1),ALPHA_GAUSS(1),
!     + loc_x(1,1),loc_y(1,1),loc_z(1,1),LGAUSS(1),MGAUSS(1),
!     + N_BAND(ik),C0(1,1,ik),N_BAND_BTM(ik),NORM_GAUSS(1),VOLUME,
!     + A_TMP(1,1)) 
!--
      subroutine make_amat(nigs,NTB,NTG,NG0,SK0,KG0,b1,b2,b3,    
     + TAU_GAUSS,ALPHA_GAUSS,loc_x,loc_y,loc_z,LGAUSS,MGAUSS,
     + N_BAND,C0,N_BAND_BTM,NORM_GAUSS,VOLUME,A_TMP) 
      implicit none
      integer::nigs,NTB,NTG
      integer::NG0
      real(8)::SK0(3)
      integer::KG0(3,NTG)
      real(8)::b1(3),b2(3),b3(3) 
      real(8)::VOLUME 
      real(8)::TAU_GAUSS(3,nigs)
      real(8)::ALPHA_GAUSS(nigs)
      real(8)::NORM_GAUSS(nigs)
      real(8)::loc_x(3,nigs),loc_y(3,nigs),loc_z(3,nigs) 
      integer::LGAUSS(nigs),MGAUSS(nigs)
      integer::N_BAND,N_BAND_BTM 
      complex(8)::C0(NTG,NTB)
      complex(8)::A_TMP(NTB,nigs)  
      !--
      real(8),parameter::pi=dacos(-1.0d0)
      complex(8),parameter::ci=(0.0D0,1.0D0) 
      integer::j_gauss,ig,i_band,jL,jM 
      real(8)::PHASE,ZETA1,ZETA2,YLM
      real(8)::VEC_KG(3,NTG)  
      real(8)::KG2(NTG) 
      complex(8)::PF1(NTG)
      real(8)::PF2(NTG) 
      real(8)::ktmp(3),tmp(3)
      complex(8)::SUM_CMPX        
      !
      !write(6,*) C0 
      !write(6,*) VOLUME  
      !write(6,*) ALPHA_GAUSS 
      !write(6,*) TAU_GAUSS 
      !write(6,*) NORM_GAUSS 
      !write(6,*) loc_x 
      !write(6,*) loc_y 
      !write(6,*) loc_z 
      !
      A_TMP=0.0d0 
      do j_gauss=1,nigs
       !--
       VEC_KG(:,:)=0.0D0             
       PF1(:)=0.0D0 
       PF2(:)=0.0D0 
       KG2(:)=0.0D0
       do ig=1,NG0 
        VEC_KG(:,ig)=(SK0(1)+DBLE(KG0(1,ig)))*b1(:)   
     +              +(SK0(2)+DBLE(KG0(2,ig)))*b2(:)   
     +              +(SK0(3)+DBLE(KG0(3,ig)))*b3(:)   
        KG2(ig)=VEC_KG(1,ig)**2+VEC_KG(2,ig)**2+VEC_KG(3,ig)**2
        PHASE=VEC_KG(1,ig)*TAU_GAUSS(1,j_gauss) 
     +       +VEC_KG(2,ig)*TAU_GAUSS(2,j_gauss) 
     +       +VEC_KG(3,ig)*TAU_GAUSS(3,j_gauss)                 
        PF1(ig)=EXP(-ci*PHASE )         
        PF2(ig)=EXP(-KG2(ig)/4.0D0/ALPHA_GAUSS(j_gauss))
       enddo!ig  
       !
       !write(6,*)'PF1'
       !write(6,*)PF1
       !write(6,*)'PF2'
       !write(6,*)PF2
       !
       do ig=1,NG0
        ktmp(:)=VEC_KG(:,ig) 
        tmp(1)=ktmp(1)*loc_x(1,j_gauss)
     +        +ktmp(2)*loc_x(2,j_gauss)
     +        +ktmp(3)*loc_x(3,j_gauss)   
        tmp(2)=ktmp(1)*loc_y(1,j_gauss)
     +        +ktmp(2)*loc_y(2,j_gauss)
     +        +ktmp(3)*loc_y(3,j_gauss) 
        tmp(3)=ktmp(1)*loc_z(1,j_gauss)
     +        +ktmp(2)*loc_z(2,j_gauss)
     +        +ktmp(3)*loc_z(3,j_gauss) 
        VEC_KG(:,ig)=tmp(:) 
       enddo!ig 
       !
       !write(6,*)'VEC_KG'
       !write(6,*)VEC_KG 
       !
       jL=LGAUSS(j_gauss) 
       jM=MGAUSS(j_gauss)     
       !--
       do i_band=1,N_BAND 
        SUM_CMPX=0.0D0                   
!<Psi|s>
        IF(jL==1)THEN                          
         do ig=1,NG0 
          ZETA1=(pi/ALPHA_GAUSS(j_gauss))**1.5D0          
          YLM=DSQRT(1.0D0/4.0D0/pi)        
          SUM_CMPX 
     +   =SUM_CMPX 
     +   +CONJG(C0(ig,i_band+N_BAND_BTM)) 
     +   *PF1(ig) 
     +   *YLM                                      
     +   *NORM_GAUSS(j_gauss) 
     +   *ZETA1*PF2(ig)       
         enddo!ig          
        ENDIF 
!<Psi|py>
        IF(jL==2.AND.jM==2)THEN
         do ig=1,NG0 
          ZETA1=(pi/ALPHA_GAUSS(j_gauss))**1.5D0          
          ZETA2=1.0D0/2.0D0/ALPHA_GAUSS(j_gauss) 
          YLM=DSQRT(3.0D0/4.0D0/pi)*VEC_KG(2,ig)
          SUM_CMPX 
     +   =SUM_CMPX 
     +   +CONJG(C0(ig,i_band+N_BAND_BTM)) 
     +   *PF1(ig)  
     +   *(-ci)*YLM                                      
     +   *NORM_GAUSS(j_gauss) 
     +   *ZETA1*ZETA2*PF2(ig)       
         enddo!ig          
        ENDIF 
!<Psi|pz>
        IF(jL==2.AND.jM==3)THEN
         do ig=1,NG0 
          ZETA1=(pi/ALPHA_GAUSS(j_gauss))**1.5D0          
          ZETA2=1.0D0/2.0D0/ALPHA_GAUSS(j_gauss) 
          YLM=DSQRT(3.0D0/4.0D0/pi)*VEC_KG(3,ig) 
          SUM_CMPX 
     +   =SUM_CMPX 
     +   +CONJG(C0(ig,i_band+N_BAND_BTM)) 
     +   *PF1(ig)            
     +   *(-ci)*YLM                                      
     +   *NORM_GAUSS(j_gauss) 
     +   *ZETA1*ZETA2*PF2(ig)  
         enddo!ig          
        ENDIF 
!<Psi|px>
        IF(jL==2.AND.jM==1)THEN
         do ig=1,NG0
          ZETA1=(pi/ALPHA_GAUSS(j_gauss))**1.5D0          
          ZETA2=1.0D0/2.0D0/ALPHA_GAUSS(j_gauss) 
          YLM=DSQRT(3.d0/4.d0/pi)*VEC_KG(1,ig)
          SUM_CMPX 
     +   =SUM_CMPX 
     +   +CONJG(C0(ig,i_band+N_BAND_BTM)) 
     +   *PF1(ig)  
     +   *(-ci)*YLM                                      
     +   *NORM_GAUSS(j_gauss) 
     +   *ZETA1*ZETA2*PF2(ig)     
         enddo!ig          
        ENDIF 
!<Psi|dxy>
        IF(jL==3.AND.jM==1)THEN                          
         do ig=1,NG0          
          ZETA1=(pi/ALPHA_GAUSS(j_gauss))**1.5D0          
          ZETA2=(1.0D0/2.0D0/ALPHA_GAUSS(j_gauss))**2.0D0   
          YLM=0.5d0*DSQRT(15.d0/pi)*VEC_KG(1,ig)*VEC_KG(2,ig)       
          SUM_CMPX 
     +   =SUM_CMPX 
     +   +CONJG(C0(ig,i_band+N_BAND_BTM)) 
     +   *PF1(ig)     
     +   *(-YLM)                                             
     +   *NORM_GAUSS(j_gauss) 
     +   *ZETA1*ZETA2*PF2(ig) 
         enddo!ig          
        ENDIF 
!<Psi|dyz>
        IF(jL==3.AND.jM==2)THEN                          
         do ig=1,NG0
          ZETA1=(pi/ALPHA_GAUSS(j_gauss))**1.5D0          
          ZETA2=(1.0D0/2.0D0/ALPHA_GAUSS(j_gauss))**2.0D0   
          YLM=(0.5d0)*DSQRT(15.d0/pi)*VEC_KG(2,ig)*VEC_KG(3,ig)
          SUM_CMPX 
     +   =SUM_CMPX 
     +   +CONJG(C0(ig,i_band+N_BAND_BTM)) 
     +   *PF1(ig)  
     +   *(-YLM)                                             
     +   *NORM_GAUSS(j_gauss) 
     +   *ZETA1*ZETA2*PF2(ig) 
         enddo!ig          
        ENDIF 
!<Psi|dz2>
        IF(jL==3.AND.jM==3)THEN                          
         do ig=1,NG0
          ZETA1=(pi/ALPHA_GAUSS(j_gauss))**1.5D0          
          ZETA2=(1.0D0/2.0D0/ ALPHA_GAUSS(j_gauss))**2.0D0     
          YLM=(0.25d0)*DSQRT(5.d0/pi)*(3.d0*VEC_KG(3,ig)**2-KG2(ig))
          SUM_CMPX 
     +   =SUM_CMPX 
     +   +CONJG(C0(ig,i_band+N_BAND_BTM)) 
     +   *PF1(ig)  
     +   *(-YLM)                                             
     +   *NORM_GAUSS(j_gauss) 
     +   *ZETA1*ZETA2*PF2(ig) 
         enddo!ig          
        ENDIF 
!<Psi|dzx>
        IF(jL==3.AND.jM==4)THEN                          
         do ig=1,NG0 
          ZETA1=(pi/ALPHA_GAUSS(j_gauss))**1.5D0          
          ZETA2=(1.0D0/2.0D0/ALPHA_GAUSS(j_gauss))**2.0D0    
          YLM=(0.5d0)*DSQRT(15.d0/pi)*VEC_KG(3,ig)*VEC_KG(1,ig)
          SUM_CMPX 
     +   =SUM_CMPX 
     +   +CONJG(C0(ig,i_band+N_BAND_BTM)) 
     +   *PF1(ig)           
     +   *(-YLM)                                             
     +   *NORM_GAUSS(j_gauss) 
     +   *ZETA1*ZETA2*PF2(ig) 
         enddo!ig          
        ENDIF 
!<Psi|dx2>
        IF(jL==3.AND.jM==5)THEN                          
         do ig=1,NG0 
          ZETA1=(pi/ALPHA_GAUSS(j_gauss))**1.5D0          
          ZETA2=(1.0D0/2.0D0/ALPHA_GAUSS(j_gauss))**2.0D0    
          YLM=(0.25d0)*DSQRT(15.d0/pi)
     +       *(VEC_KG(1,ig)**2-VEC_KG(2,ig)**2) 
          SUM_CMPX 
     +   =SUM_CMPX 
     +   +CONJG(C0(ig,i_band+N_BAND_BTM)) 
     +   *PF1(ig)              
     +   *(-YLM)                                             
     +   *NORM_GAUSS(j_gauss) 
     +   *ZETA1*ZETA2*PF2(ig)     
         enddo!ig          
        ENDIF 
        !--
        A_TMP(i_band,j_gauss)=SUM_CMPX/DSQRT(VOLUME)             
       enddo!i_band         
      enddo!j_gauss 
      ! 
      !do i_band=1,N_BAND
      ! write(6,'(100f12.5)')(A_TMP(i_band,j_gauss),j_gauss=1,nigs)
      !enddo 
      !write(6,*) 
      ! 
      return
      end 
      !--
      subroutine make_pw(C0_I,wfunc,fftwk,NG0,KG0,NTG,nwx2,nwy2,nwz2,
     + nfft1,nfft2,Nl123,fs,xmin,xmax,ymin,ymax,zmin,zmax,SK0,pw) 
      use fft_3d 
      implicit none
      type(fft3_struct)::fs 
      integer::NTG,nwx2,nwy2,nwz2,nfft1,nfft2,Nl123
      integer::xmin,xmax,ymin,ymax,zmin,zmax
      integer::KG0(3,NTG)
      integer::NG0 
      real(8)::fftwk(Nl123*2)
      real(8)::wfunc(Nl123*2) 
      real(8)::SK0(3) 
      complex(8)::C0_I(NTG) 
      complex(8)::pw((xmin)*nwx2:(xmax)*nwx2-1,
     +               (ymin)*nwy2:(ymax)*nwy2-1, 
     +               (zmin)*nwz2:(zmax)*nwz2-1)
      integer::ig,igb1,igb2,igb3,ind,ir,ir1,ir2,ir3 
      integer::ir1_SHIFT,ir2_SHIFT,ir3_SHIFT,jr1,jr2,jr3  
      real(8)::PHASE
      complex(8)::pf 
      real(8),parameter::tpi=2.0d0*dacos(-1.0d0)
      complex(8),parameter::ci=(0.0D0,1.0D0) 
      complex(8)::u_1D(Nl123) 
      complex(8)::u_3D(0:nwx2-1,0:nwy2-1,0:nwz2-1) 
      !--
      u_3D=0.0D0 
      u_1D=0.0D0 
      wfunc=0.0d0 
      fftwk=0.0d0 
      do ig=1,NG0 
       igb1=KG0(1,ig) 
       igb2=KG0(2,ig) 
       igb3=KG0(3,ig) 
       igb1=MOD(nwx2+igb1,nwx2)+1
       igb2=MOD(nwy2+igb2,nwy2)+1
       igb3=MOD(nwz2+igb3,nwz2)+1
       ind=igb1+(igb2-1)*nfft1+(igb3-1)*nfft1*nfft2 
       wfunc(ind)=dble(C0_I(ig)) 
       wfunc(ind+Nl123)=dimag(C0_I(ig)) 
      enddo!ig 
      call fft3_bw(fs,wfunc(1),fftwk(1)) 
      do ir=1,Nl123
       u_1D(ir)=cmplx(wfunc(ir),wfunc(ir+Nl123)) 
      enddo!ir 
      !--
      do ir3=1,nwz2
       do ir2=1,nwy2
        do ir1=1,nwx2 
         ir=ir1+(ir2-1)*nfft1+(ir3-1)*nfft1*nfft2  
         u_3D(ir1-1,ir2-1,ir3-1)=u_1D(ir) 
        enddo 
       enddo 
      enddo 
      !--
      do ir1=(xmin)*nwx2,(xmax)*nwx2-1         
       do ir2=(ymin)*nwy2,(ymax)*nwy2-1        
        do ir3=(zmin)*nwz2,(zmax)*nwz2-1         
         PHASE=tpi*(DBLE(ir1)*SK0(1)/DBLE(nwx2) 
     +             +DBLE(ir2)*SK0(2)/DBLE(nwy2) 
     +             +DBLE(ir3)*SK0(3)/DBLE(nwz2)) 
         pf=EXP(ci*PHASE)          
         !(ir1,ir2,ir3)->(jr1,jr2,jr3)
         ir1_SHIFT=ir1       
         ir2_SHIFT=ir2       
         ir3_SHIFT=ir3       
         DO WHILE(ir1_SHIFT<0) 
          ir1_SHIFT=ir1_SHIFT+nwx2               
         ENDDO        
         DO WHILE(ir2_SHIFT<0) 
          ir2_SHIFT=ir2_SHIFT+nwy2                              
         ENDDO        
         DO WHILE(ir3_SHIFT<0) 
          ir3_SHIFT=ir3_SHIFT+nwz2            
         ENDDO        
         jr1=MOD(ir1_SHIFT,nwx2)       
         jr2=MOD(ir2_SHIFT,nwy2)       
         jr3=MOD(ir3_SHIFT,nwz2)       
         !CALC wannier(ir1,ir2,ir3)
         pw(ir1,ir2,ir3)=pw(ir1,ir2,ir3)+u_3D(jr1,jr2,jr3)*pf 
        enddo!ir3         
       enddo!ir2 
      enddo!ir1 
      !-- 
      return
      end
      !--
      subroutine make_bf(C0_I,wfunc,fftwk,NG0,KG0,NTG,nwx2,nwy2,nwz2,
     + nfft1,nfft2,Nl123,fs,pw) 
      use fft_3d 
      implicit none
      type(fft3_struct)::fs 
      integer::NTG,nwx2,nwy2,nwz2,nfft1,nfft2,Nl123
      integer::KG0(3,NTG)
      integer::NG0 
      real(8)::fftwk(Nl123*2)
      real(8)::wfunc(Nl123*2) 
      complex(8)::C0_I(NTG) 
      integer::ig,igb1,igb2,igb3,ind,ir,ir1,ir2,ir3 
      complex(8)::u_1D(Nl123) 
      real(8)::pw(0:nwx2-1,0:nwy2-1,0:nwz2-1)
      !--
      u_1D=0.0d0 
      wfunc=0.0d0 
      fftwk=0.0d0 
      do ig=1,NG0 
       igb1=KG0(1,ig) 
       igb2=KG0(2,ig) 
       igb3=KG0(3,ig) 
       igb1=MOD(nwx2+igb1,nwx2)+1
       igb2=MOD(nwy2+igb2,nwy2)+1
       igb3=MOD(nwz2+igb3,nwz2)+1
       ind=igb1+(igb2-1)*nfft1+(igb3-1)*nfft1*nfft2 
       wfunc(ind)=dble(C0_I(ig)) 
       wfunc(ind+Nl123)=dimag(C0_I(ig)) 
      enddo!ig 
      call fft3_bw(fs,wfunc(1),fftwk(1)) 
      do ir=1,Nl123
       u_1D(ir)=cmplx(wfunc(ir),wfunc(ir+Nl123)) 
      enddo!ir 
      !--
      do ir3=1,nwz2
       do ir2=1,nwy2
        do ir1=1,nwx2 
         ir=ir1+(ir2-1)*nfft1+(ir3-1)*nfft1*nfft2  
         pw(ir1-1,ir2-1,ir3-1)=abs(u_1D(ir))**2 
        enddo 
       enddo 
      enddo 
      !--
      return
      end
      !--
      subroutine search_ik(NTK,SK0,ktmp,ik) 
      implicit none 
      integer::NTK
      real(8)::SK0(3,NTK) 
      real(8)::ktmp(3)
      integer::iik,ik
      real(8),parameter::dlt_BZ=1.0d-5 
      !--
      ik=0 
      do iik=1,NTK 
       !write(6,*)SK0(:,iik) 
       if( abs(SK0(1,iik)-ktmp(1))<dlt_BZ.and. 
     +     abs(SK0(2,iik)-ktmp(2))<dlt_BZ.and. 
     +     abs(SK0(3,iik)-ktmp(3))<dlt_BZ )then 
        ik=iik 
       endif 
      enddo!iik
      if(ik==0)then 
       write(6,*)'Failed search ik:'
       write(6,*)'There exits no k pts to be calculated:'
       write(6,*)'Then, stop'
       stop 
      endif 
      !--
      return
      end 
      !
      !20180817
      !
      subroutine search_vis_range(ktmp,amin,amax,bmin,bmax,cmin,cmax) 
      implicit none 
      real(8),intent(in)::ktmp(3)
      integer,intent(out)::amin,amax,bmin,bmax,cmin,cmax 
      real(8)::denom
      real(8),parameter::dlt_BZ=1.0d-5 
      amin=0; amax=0; bmin=0; bmax=0; cmin=0; cmax=0
      !
      if(dabs(ktmp(1))<dlt_BZ)then 
       denom=1.0d0 
      else
       denom=ktmp(1) 
      endif 
      amax=nint(1.0d0/denom)
      !
      if(dabs(ktmp(2))<dlt_BZ)then 
       denom=1.0d0 
      else
       denom=ktmp(2) 
      endif 
      bmax=nint(1.0d0/denom)
      !
      if(dabs(ktmp(3))<dlt_BZ)then 
       denom=1.0d0 
      else
       denom=ktmp(3) 
      endif 
      cmax=nint(1.0d0/denom)
      !
      return 
      end
      !
      !20180820
      !
      subroutine phase_fix_bf(nwx2,nwy2,nwz2,xmin,xmax,ymin,ymax,
     + zmin,zmax,pw) 
      implicit none
      integer,intent(in)::nwx2,nwy2,nwz2,xmin,xmax,ymin,ymax,zmin,zmax
      complex(8),intent(inout)::pw((xmin)*nwx2:(xmax)*nwx2-1,
     +                             (ymin)*nwy2:(ymax)*nwy2-1, 
     +                             (zmin)*nwz2:(zmax)*nwz2-1)
      integer::ir1,ir2,ir3 
      real(8)::x,y,r,rmax 
      complex(8)::glb_phase 
      complex(8),parameter::ci=(0.0d0,1.0d0) 
      ! 
      rmax=1.0d-6 
      do ir1=(xmin)*nwx2,(xmax)*nwx2-1         
       do ir2=(ymin)*nwy2,(ymax)*nwy2-1        
        do ir3=(zmin)*nwz2,(zmax)*nwz2-1         
         x=dble(pw(ir1,ir2,ir3))  
         y=dimag(pw(ir1,ir2,ir3))  
         r=dsqrt(x**2+y**2)
         if(r>rmax)then
          glb_phase=(x-ci*y)/r 
          rmax=r 
         endif 
        enddo!ir3         
       enddo!ir2 
      enddo!ir1 
      !write(6,*)'rmax=',rmax 
      pw=pw*glb_phase 
      !
      return
      end
