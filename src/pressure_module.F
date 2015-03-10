!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!      This module is for adding non-hydrostatic pressure to the calculation using 
!        Green-Nhagdhi Pressure Poisson Finite Difference formulations
!
!       
!      pressure_module.F
!
!      Jessica Meixner and Aaron Donahue  Feb 2015 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      MODULE pressure_module

      USE SIZES,    ONLY: SZ

      INTEGER LengthX, LengthY 
      INTEGER NumCols, NumRows 


      !Variables that will be passed from DG-SWEM to Pressure Poisson 
      REAL(SZ), ALLOCATABLE :: totaldepth_pp(:,:) !total depth	
      REAL(SZ), ALLOCATABLE :: eta_pp(:,:)        !water surface elevation
      REAL(SZ), ALLOCATABLE :: u_pp(:,:)          !u- current in x direction
      REAL(SZ), ALLOCATABLE :: v_pp(:,:)          !v- current in y direction
      REAL(SZ), ALLOCATABLE :: depth_pp(:,:)       !bathymetry 
      INTEGER,  ALLOCATABLE :: wetdry_pp(:,:)     !wet (1) dry (0) flag 

      REAL(SZ), ALLOCATABLE :: Etax_pp(:,:),Etay_pp(:,:),Etaxx_pp(:,:),Etayy_pp(:,:) !derivatives of surface elevation 
      REAL(SZ), ALLOCATABLE :: Ux_pp(:,:),Uy_pp(:,:),Vx_pp(:,:),Vy_pp(:,:)        !derivatives of current 
      REAL(SZ), ALLOCATABLE :: Depx_pp(:,:),Depy_pp(:,:),Depxx_pp(:,:),Depyy_pp(:,:),Depxy_pp(:,:)   !derivatevs of bathymetry 
      REAL(SZ), ALLOCATABLE :: p1_pp(:,:),p1x_pp(:,:),p1y_pp(:,:)

      REAL(SZ)              :: DX,DY
      REAL(SZ), ALLOCATABLE :: pd_pp(:,:),pb_pp(:,:)
      REAL(SZ), ALLOCATABLE :: p1old(:),p1oldest(:)  
     
      CONTAINS 

      SUBROUTINE  read_and_initialize_PPvalues 

      USE global,   ONLY: DP 

      IMPLICIT NONE 

      INTEGER i, j, inode


      !read in what is LengthInX, LengthInY 
      !read in what is NumCols, NumRows
!      NumCols = 
!      NumRows = 
!      DX = 
!      DY = 
  
      !read in what are the different boundary types: 
         !only 1 (reflective) boundary type for now. 

      ALLOCATE(totaldepth_pp(NumRows,NumCols)) 
      ALLOCATE(eta_pp(NumRows,NumCols))
      ALLOCATE(u_pp(NumRows,NumCols))
      ALLOCATE(v_pp(NumRows,NumCols))
      ALLOCATE(depth_pp(NumRows,NumCols))
      ALLOCATE(wetdry_pp(NumRows,NumCols))
      ALLOCATE(Etax_pp(NumRows,NumCols),Etay_pp(NumRows,NumCols))
      ALLOCATE(Etaxx_pp(NumRows,NumCols),Etayy_pp(NumRows,NumCols))
      ALLOCATE(Ux_pp(NumRows,NumCols),Uy_pp(NumRows,NumCols))
      ALLOCATE(Vx_pp(NumRows,NumCols),Vy_pp(NumRows,NumCols))      
      ALLOCATE(Depx_pp(NumRows,NumCols),Depy_pp(NumRows,NumCols))
      ALLOCATE(Depxx_pp(NumRows,NumCols),Depyy_pp(NumRows,NumCols))
      ALLOCATE(Depxy_pp(NumRows,NumCols))   
      ALLOCATE(p1_pp(NumRows,NumCols),p1x_pp(NumRows,NumCols),p1y_pp(NumRows,NumCols))
      ALLOCATE(pd_pp(NumRows,NumCols),pb_pp(NumRows,NumCols))
      ALLOCATE(p1old(NumRows*NumCols),p1oldest(NumRows*NumCols))  

 
      !Inititalize values:     
      totaldepth_pp = 0.d0 
      eta_pp = 0.d0 
      u_pp = 0.d0 
      v_pp = 0.d0 
      wetdry_pp = 1
      Etax_pp=0.d0; Etay_pp=0.d0 
      Etaxx_pp=0.d0; Etayy_pp=0.d0
      Ux_pp=0.d0;Uy_pp=0.d0
      Vx_pp=0.d0;Vy_pp=0.d0 
      Depx_pp=0.d0;Depy_pp=0.d0
      Depxx_pp=0.d0; Depyy_pp=0.d0
      Depxy_pp=0.d0
      p1_pp=0.d0;p1x_pp=0.d0;p1y_pp=0.d0
      pd_pp=0.d0;pb_pp=0.d0
      p1ol=0.d0;p1oldest=0.d0 
 
            


 
      inode = 1
      DO i = 1, NumRows
        DO j = 1, NumCols
           depth_pp(i,j) = DP(inode)
           inode = inode + 1
        END DO !i 
      END DO !j 

      END SUBROUTINE  read_and_initialize_PPvalues 

      SUBROUTINE  from_DG_to_Pressure 

      !Need to obtain u,v,eta,depth and wet/dry at each node.  then need to pass to (i,j) type format 

      USE GLOBAL, ONLY: MNP, EL_COUNT, ELETAB, NM,
     &                  ZE_DG, QX_DG, QY_DG, HB_DG, 
     &                  AREAS, ETA2, UU2, VV2, 
     &                  NLEQ, LEQ 
      USE DG,     ONLY: WDFLG, ZE, QX, QY, HB, DOFH, 
     &                  PHI_CORNER, ph   

      IMPLICIT NONE 

C.....Declare local variables

      INTEGER I, NO_NBORS, J,K,N1,KK, NBOR_EL, inode 
      REAL(SZ) AREA_SUM, AREA, DEPTH, FH_NL
      REAL(SZ), ALLOCATABLE :: TOTALDEPTH(:)
      INTEGER, ALLOCATABLE :: WDFLG_VERT(:) 

      ALLOCATE(TOTALDEPTH(MNP))
      ALLOCATE(WDFLG_VERT(MNP))
      WDFLG_VERT = 1 !Start as all being wet, will set to dry if any element connected to vertex is dry 

      !Obtain u,v,eta, depth and wet/dry flag at each node 

      DO I = 1,MNP

        NO_NBORS = EL_COUNT(I)

        AREA_SUM = 0
        DO 333 J = 1,NO_NBORS
          NBOR_EL = ELETAB(I,1+J)

          IF(WDFLG(NBOR_EL).EQ.0) THEN 
            WDFLG_VERT(I) = 0 !if any neighboring elem is dry, consider vertex dry 
            CYCLE ! DON'T COUNT DRY ELEMENTS  sb 02/26/07
          END IF 

          DO K = 1,3
            N1 = NM(NBOR_EL,K)
            IF (N1.EQ.I) THEN
              ZE_DG(J) = ZE(1,NBOR_EL,1)
              QX_DG(J) = QX(1,NBOR_EL,1)
              QY_DG(J) = QY(1,NBOR_EL,1)
              HB_DG(J) = HB(1,NBOR_EL,1)
              DO KK = 2,DOFH
                ZE_DG(J) = ZE_DG(J) + PHI_CORNER(KK,K,ph)*ZE(KK,NBOR_EL,1)
                QX_DG(J) = QX_DG(J) + PHI_CORNER(KK,K,ph)*QX(KK,NBOR_EL,1)
                QY_DG(J) = QY_DG(J) + PHI_CORNER(KK,K,ph)*QY(KK,NBOR_EL,1)
                HB_DG(J) = HB_DG(J) + PHI_CORNER(KK,K,ph)*HB(KK,NBOR_EL,1)
              ENDDO
              AREA = 0.5*AREAS(NBOR_EL)
              AREA_SUM = AREA_SUM + AREA
              GOTO 333
            ENDIF
          ENDDO
 333    CONTINUE

        ETA2(I) = 0.D0
        UU2(I)  = 0.D0
        VV2(I)  = 0.D0
        TOTALDEPTH(I) = 0.d0 


        DO J = 1,NO_NBORS
          NBOR_EL = ELETAB(I,1+J)

          IF(WDFLG(NBOR_EL).EQ.0) CYCLE ! DON'T COUNT DRY ELEMENTS  sb 02/26/07

          AREA = 0.5*AREAS(NBOR_EL)/AREA_SUM
          DEPTH = ZE_DG(J) + HB_DG(J)
          FH_NL = 1.D0/(NLEQ*DEPTH + LEQ)
          ETA2(I) = ETA2(I) + AREA*ZE_DG(J)
          UU2(I)  = UU2(I)  + AREA*QX_DG(J)*FH_NL
          VV2(I)  = VV2(I)  + AREA*QY_DG(J)*FH_NL
          TOTALDEPTH(I) = TOTALDEPTH(I) + AREA*DEPTH  
        ENDDO
      ENDDO



      !now that we have ETA2,UU2,VV2,TOTALDEPTH, and WDFLG_VERT 
      !need to go from a node numbering to a (i,j) vertex numbering PP FD can use.  
      inode = 1 
      DO i = 1, NumRows 
        DO j = 1, NumCols 
           totaldepth_pp(i,j) = TOTALDEPTH(inode) 
           eta_pp(i,j) = ETA2(inode) 
           u_pp(i,j) = UU2(inode) 
           v_pp(i,j) = VV2(inode) 
           wetdry_pp(i,j) = WDFLG_VERT(inode)
           inode = inode + 1 
        END DO !i 
      END DO !j 


      CALL pressure_WASUPP 
      
      END SUBROUTINE from_DG_to_Pressure      


      END MODULE

