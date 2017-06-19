C***********************************************************************
C     
C     SUBROUTINE INTERNAL_EDGE_HYDRO( )
C     
C     This subroutine does the following:
C     
C     1.  Calculates the values of the necessary variables at the edge
C     gauss points for INTERNAL edges
C     2.  Calls the appropriate subroutine to compute the flux at
C     these points.
C     3.  Calls the appropriate subroutine to compute the boundary
C     integrals.
C     
C     Written by Ethan Kubatko (06-11-2004)
C     
C     01-10-2011 - cem - adapted for p_enrichment and multicomponent
C     11-11-2011 - cem - adapted for layered sediment
C     
C-----------------------------------------------------------------------
C     
C     01-02-2007, sb, Modified for LDG
C     C***********************************************************************

      SUBROUTINE INTERNAL_EDGE_HYDRO(IT)

C.....Use appropriate modules

      USE GLOBAL
      USE DG
      USE NodalAttributes, ONLY :  ESLM
      use sizes, only:  sz, myproc, layers

      IMPLICIT NONE

C.....Declare local variables

      INTEGER L, LED_IN, LED_EX, GED, GP_IN, GP_EX,k,i,ll,IT
      REAL(SZ), PARAMETER :: ZERO = 1.D-12
      REAL(SZ) TX, TY, W_IN, W_EX
      REAL(SZ) LZ_XX_IN, LZ_XY_IN, LZ_YX_IN, LZ_YY_IN
      REAL(SZ) LZ_XX_EX, LZ_XY_EX, LZ_YX_EX, LZ_YY_EX
      REAL(SZ) HZ_X_EX, HZ_Y_EX, HZ_X_IN, HZ_Y_IN
      REAL(SZ) TZ_X_EX, TZ_Y_EX, TZ_X_IN, TZ_Y_IN
      REAL(SZ) EDFAC_IN, EDFAC_EX, DEN
      REAL(SZ) XLEN_EL_IN, XLEN_EL_EX,MZ_X_EX(layers),MZ_Y_EX(layers)
      REAL(SZ) MASS_EL_IN, MASS_EL_EX,MZ_X_IN(layers),MZ_Y_IN(layers)
      REAL(SZ), SAVE, ALLOCATABLE :: 
     &     RHS_ZE_IN(:), RHS_QX_IN(:), RHS_QY_IN(:),
     &     RHS_ZE_EX(:), RHS_QX_EX(:), RHS_QY_EX(:)
#ifdef TRACE     
     &     ,RHS_iota_IN(:), RHS_iota_EX(:)
#endif
#ifdef CHEM
     &     ,RHS_iota_IN(:), RHS_iota_EX(:)
     &     ,RHS_iota2_IN(:), RHS_iota2_EX(:)
#endif
#ifdef DYNP     
     &     ,RHS_dynP_IN(:), RHS_dynP_EX(:)
#endif

      REAL(SZ) ARK, BRK
      REAL(SZ) MAX_BOA          ! Maximum of beta_il/alpha_il for all l
      REAL(SZ) NLEQG_TMP, G_TMP
      REAL(SZ) F_HAT_O, G_HAT_O, H_HAT_O,i_hat_o,j_hat_o
      REAL(SZ) G_HAT_IN, H_HAT_IN
      REAL(SZ) G_HAT_EX, H_HAT_EX
      REAL(SZ) K_HAT_O
      Real(SZ) bed_HAT_IN(layers),bed_HAT_EX(layers)
c     

      IF(.NOT.ALLOCATED(RHS_ZE_IN)) THEN
         ALLOCATE ( RHS_ZE_IN(DOFH),RHS_QX_IN(DOFH),RHS_QY_IN(DOFH) )
         ALLOCATE ( RHS_ZE_EX(DOFH),RHS_QX_EX(DOFH),RHS_QY_EX(DOFH) )
#ifdef TRACE
         Allocate ( rhs_iota_in(dofh),rhs_iota_ex(dofh) )
#endif
#ifdef CHEM
         Allocate ( rhs_iota_in(dofh),rhs_iota2_in(dofh) )
         Allocate ( rhs_iota_ex(dofh),rhs_iota2_ex(dofh) )
#endif
#ifdef DYNP
         Allocate ( rhs_dynP_IN(dofh),rhs_dynP_EX(dofh) )
#endif
      ENDIF

      test_el = 0
      DO 1000 L = 1,NIEDS

C.......Retrieve the global and local edge number

         GED = NIEDN(L)
         LED_IN = NEDSD(1,GED)
         LED_EX = NEDSD(2,GED)

C.......Retrieve the elements which share the edge

         EL_IN = NEDEL(1,GED)
         EL_EX = NEDEL(2,GED)

         EL = EL_EX

         if (DOFS(EL_EX).LT.DOFS(EL_IN)) then
            EL = EL_IN
         endif
         
         pa = PDG_EL(EL)

#ifdef P0
         if (pa.eq.0) then
            pa = 1
         endif
#endif

C.......If both elements on either side of edge are dry then skip

         IF((WDFLG(EL_IN).EQ.0).AND.(WDFLG(EL_EX).EQ.0)) GOTO 1000

         !test_el = test_el+1

C.....Save the current RHS

         DO K = 1,DOFS(EL)

            RHS_ZE_IN(K) = RHS_ZE(K,EL_IN,IRK)
            RHS_QX_IN(K) = RHS_QX(K,EL_IN,IRK)
            RHS_QY_IN(K) = RHS_QY(K,EL_IN,IRK)

            RHS_ZE_EX(K) = RHS_ZE(K,EL_EX,IRK)
            RHS_QX_EX(K) = RHS_QX(K,EL_EX,IRK)
            RHS_QY_EX(K) = RHS_QY(K,EL_EX,IRK)

#ifdef TRACE
            RHS_iota_IN(K) = RHS_iota(K,EL_IN,IRK)
            RHS_iota_EX(K) = RHS_iota(K,EL_EX,IRK)
#endif

#ifdef CHEM
            RHS_iota_IN(K) = RHS_iota(K,EL_IN,IRK)
            RHS_iota2_IN(K) = RHS_iota2(K,EL_IN,IRK)
            RHS_iota_EX(K) = RHS_iota(K,EL_EX,IRK)
            RHS_iota2_EX(K) = RHS_iota2(K,EL_EX,IRK)
#endif

#ifdef DYNP
            RHS_dynP_IN(K) = RHS_dynP(K,EL_IN,IRK)
            RHS_dynP_EX(K) = RHS_dynP(K,EL_EX,IRK)
#endif

#ifdef SED_LAY
            do ll=1,layers
               RHS_bed_IN(K,ll) = RHS_bed(k,EL_IN,IRK,ll)
               RHS_bed_EX(K,ll) = RHS_bed(k,EL_EX,IRK,ll)
            enddo
#endif

         ENDDO

C.....Compute the sum of the lengths of three edges

         XLEN_EL_IN = XLEN(NELED(1,EL_IN))
         XLEN_EL_IN = XLEN_EL_IN + XLEN(NELED(2,EL_IN))
         XLEN_EL_IN = XLEN_EL_IN + XLEN(NELED(3,EL_IN))

         XLEN_EL_EX = XLEN(NELED(1,EL_EX))
         XLEN_EL_EX = XLEN_EL_EX + XLEN(NELED(2,EL_EX))
         XLEN_EL_EX = XLEN_EL_EX + XLEN(NELED(3,EL_EX))

C.....Compute the total mass in the elements

         MASS_EL_IN = (ZE(1,EL_IN,IRK)+HB(1,EL_IN,1))*AREAS(EL_IN)*0.5D0
         MASS_EL_EX = (ZE(1,EL_EX,IRK)+HB(1,EL_EX,1))*AREAS(EL_EX)*0.5D0

#ifdef SED_LAY
         HB(1,EL_IN,irk) = 0.D0
         HB(1,EL_EX,irk) = 0.D0
         do ll=1,layers
            HB(1,EL_IN,irk) = HB(1,EL_IN,irk) + bed(1,EL_IN,irk,ll)
            HB(1,EL_EX,irk) = HB(1,EL_EX,irk) + bed(1,EL_EX,irk,ll)
         enddo
         MASS_EL_IN = (ZE(1,EL_IN,IRK)+HB(1,EL_IN,1))*AREAS(EL_IN)*0.5D0
         MASS_EL_EX = (ZE(1,EL_EX,IRK)+HB(1,EL_EX,1))*AREAS(EL_EX)*0.5D0
#endif
         
C.....Retrieve the components of the normal vector to the edge

         NX = COSNX(GED)
         NY = SINNX(GED)
         
C.....Set the components for the tangential vector to the edge

         TX = -NY
         TY =  NX
         
         EDFAC_IN = XLEN(GED)/AREAS(EL_IN)
         EDFAC_EX = XLEN(GED)/AREAS(EL_EX)

C.....Compute ZE, QX, QY, and HB at each edge Gauss quadrature point

         DO I = 1,NEGP(pa)

            GP_IN = I
            GP_EX = NEGP(pa) - I + 1

            ZE_IN = ZE(1,EL_IN,IRK)
            QX_IN = QX(1,EL_IN,IRK)
            QY_IN = QY(1,EL_IN,IRK)
            HB_IN = BATHED(GP_IN,LED_IN,EL_IN,pa)

#ifdef SED_LAY
            do ll = 1,layers
               bed_IN(ll) = bed(1,EL_IN,irk,ll)
            enddo
            HB_IN = HB(1,EL_IN,irk)
#endif

            SFAC_IN = SFACED(GP_IN,LED_IN,EL_IN,pa)
            
            ZE_EX = ZE(1,EL_EX,IRK)
            QX_EX = QX(1,EL_EX,IRK)
            QY_EX = QY(1,EL_EX,IRK)
            HB_EX = BATHED(GP_EX,LED_EX,EL_EX,pa)

#ifdef SED_LAY 
            do ll = 1,layers
               bed_EX(ll) = bed(1,EL_EX,irk,ll)
            enddo
            HB_EX = HB(1,EL_EX,irk)
#endif

            SFAC_EX = SFACED(GP_EX,LED_EX,EL_EX,pa)


#ifdef TRACE
            iota_IN = iota(1,EL_IN,IRK)
            iota_EX = iota(1,EL_EX,IRK)

            TZ_X_IN = TZ(1,1,1,EL_IN)
            TZ_Y_IN = TZ(1,2,2,EL_IN)
            TZ_X_EX = TZ(1,1,1,EL_EX)
            TZ_Y_EX = TZ(1,2,2,EL_EX)
#endif

#ifdef CHEM
            iota_IN = iota(1,EL_IN,IRK)
            iota2_IN = iota2(1,EL_IN,IRK)
            iota_EX = iota(1,EL_EX,IRK)
            iota2_EX = iota2(1,EL_EX,IRK)
#endif

#ifdef DYNP
            dynP_IN = dynP(1,EL_IN,IRK)
            dynP_EX = dynP(1,EL_EX,IRK)
#endif


#ifdef SED_LAY
            do ll = 1,layers
               bed_IN(ll) = bed(1,EL_IN,IRK,ll)
               bed_EX(ll) = bed(1,EL_EX,IRK,ll)
               
               MZ_X_IN(ll) = MZ(1,1,ll,EL_IN)
               MZ_Y_IN(ll) = MZ(1,2,ll,EL_IN)  

               MZ_X_EX(ll) = MZ(1,1,ll,EL_EX)
               MZ_Y_EX(ll) = MZ(1,2,ll,EL_EX)  
            enddo
#endif

                                !LDG terms
#ifdef WAVE_DIF
            HZ_X_IN = HZ(1,1,1,EL_IN)
            HZ_Y_IN = HZ(1,2,2,EL_IN)
            HZ_X_EX = HZ(1,1,1,EL_EX)
            HZ_Y_EX = HZ(1,2,2,EL_EX)
#endif
            LZ_XX_IN = LZ(1,1,1,EL_IN)
            LZ_XY_IN = LZ(1,1,2,EL_IN)
            LZ_YX_IN = LZ(1,2,1,EL_IN)
            LZ_YY_IN = LZ(1,2,2,EL_IN)

            LZ_XX_EX = LZ(1,1,1,EL_EX)
            LZ_XY_EX = LZ(1,1,2,EL_EX)
            LZ_YX_EX = LZ(1,2,1,EL_EX)
            LZ_YY_EX = LZ(1,2,2,EL_EX)

            DO K = 2,DOFS(EL)

             ZE_IN = ZE_IN + ZE(K,EL_IN,IRK)*PHI_EDGE(K,GP_IN,LED_IN,pa)
             QX_IN = QX_IN + QX(K,EL_IN,IRK)*PHI_EDGE(K,GP_IN,LED_IN,pa)
             QY_IN = QY_IN + QY(K,EL_IN,IRK)*PHI_EDGE(K,GP_IN,LED_IN,pa)

             ZE_EX = ZE_EX + ZE(K,EL_EX,IRK)*PHI_EDGE(K,GP_EX,LED_EX,pa)
             QX_EX = QX_EX + QX(K,EL_EX,IRK)*PHI_EDGE(K,GP_EX,LED_EX,pa)
             QY_EX = QY_EX + QY(K,EL_EX,IRK)*PHI_EDGE(K,GP_EX,LED_EX,pa)

#ifdef TRACE
       iota_IN = iota_IN + iota(K,EL_IN,IRK)*PHI_EDGE(K,GP_IN,LED_IN,pa)
       iota_EX = iota_EX + iota(K,EL_EX,IRK)*PHI_EDGE(K,GP_EX,LED_EX,pa)
#endif

#ifdef CHEM
      iota_IN = iota_IN  + iota(K,EL_IN,IRK)*PHI_EDGE(K,GP_IN,LED_IN,pa)
      iota_EX = iota_EX  + iota(K,EL_EX,IRK)*PHI_EDGE(K,GP_EX,LED_EX,pa)

      iota2_IN = iota2_IN+iota2(K,EL_IN,IRK)*PHI_EDGE(K,GP_IN,LED_IN,pa)
      iota2_EX=iota2_EX + iota2(K,EL_EX,IRK)*PHI_EDGE(K,GP_EX,LED_EX,pa)
#endif

#ifdef DYNP
       dynP_IN = dynP_IN + dynP(K,EL_IN,IRK)*PHI_EDGE(K,GP_IN,LED_IN,pa)
       dynP_EX = dynP_EX + dynP(K,EL_EX,IRK)*PHI_EDGE(K,GP_EX,LED_EX,pa)
#endif

#ifdef SED_LAY
               do ll=1,layers
                  bed_IN(ll) = bed_IN(ll) +
     & bed(K,EL_IN,IRK,ll)*PHI_EDGE(K,GP_IN,LED_IN,pa)
                  bed_EX(ll) = bed_EX(ll) +
     & bed(K,EL_EX,IRK,ll)*PHI_EDGE(K,GP_EX,LED_EX,pa)
                  HB_IN = HB_IN +
     & bed(k,EL_IN,irk,ll)*PHI_EDGE(K,GP_IN,LED_IN,pa)
                  HB_EX = HB_EX +
     & bed(k,EL_EX,irk,ll)*PHI_EDGE(K,GP_EX,LED_EX,pa)

                  MZ_X_IN(ll) = MZ_X_IN(ll) +
     & MZ(K,1,ll,EL_IN)*PHI_EDGE(K,GP_IN,LED_IN,pa)
                  MZ_Y_IN(ll) = MZ_Y_IN(ll) +
     & MZ(K,2,ll,EL_IN)*PHI_EDGE(K,GP_IN,LED_IN,pa)
                  
                  MZ_X_EX(ll) = MZ_X_EX(ll) +
     & MZ(K,1,ll,EL_EX)*PHI_EDGE(K,GP_EX,LED_EX,pa)
                  MZ_Y_EX(ll) = MZ_Y_EX(ll) +
     & MZ(K,2,ll,EL_EX)*PHI_EDGE(K,GP_EX,LED_EX,pa)
               enddo
#endif

                                ! LDG terms
#ifdef WAVE_DIF
               HZ_X_IN = HZ_X_IN 
     &              + HZ(K,1,1,EL_IN)*PHI_EDGE(K,GP_IN,LED_IN,pa)
               HZ_Y_IN = HZ_Y_IN 
     &              + HZ(K,2,2,EL_IN)*PHI_EDGE(K,GP_IN,LED_IN,pa)
               HZ_X_EX = HZ_X_EX 
     &              + HZ(K,1,1,EL_EX)*PHI_EDGE(K,GP_EX,LED_EX,pa)
               HZ_Y_EX = HZ_Y_EX 
     &              + HZ(K,2,2,EL_EX)*PHI_EDGE(K,GP_EX,LED_EX,pa)
#endif
               LZ_XX_IN = LZ_XX_IN 
     &              + LZ(K,1,1,EL_IN)*PHI_EDGE(K,GP_IN,LED_IN,pa)
               LZ_XY_IN = LZ_XY_IN 
     &              + LZ(K,1,2,EL_IN)*PHI_EDGE(K,GP_IN,LED_IN,pa)
               LZ_YX_IN = LZ_YX_IN 
     &              + LZ(K,2,1,EL_IN)*PHI_EDGE(K,GP_IN,LED_IN,pa)
               LZ_YY_IN = LZ_YY_IN 
     &              + LZ(K,2,2,EL_IN)*PHI_EDGE(K,GP_IN,LED_IN,pa)

               LZ_XX_EX = LZ_XX_EX 
     &              + LZ(K,1,1,EL_EX)*PHI_EDGE(K,GP_EX,LED_EX,pa)
               LZ_XY_EX = LZ_XY_EX 
     &              + LZ(K,1,2,EL_EX)*PHI_EDGE(K,GP_EX,LED_EX,pa)
               LZ_YX_EX = LZ_YX_EX 
     &              + LZ(K,2,1,EL_EX)*PHI_EDGE(K,GP_EX,LED_EX,pa)
               LZ_YY_EX = LZ_YY_EX 
     &              + LZ(K,2,2,EL_EX)*PHI_EDGE(K,GP_EX,LED_EX,pa)

#ifdef TRACE
               TZ_X_IN = TZ_X_IN 
     &              + TZ(K,1,1,EL_IN)*PHI_EDGE(K,GP_IN,LED_IN,pa)
               TZ_Y_IN = TZ_Y_IN 
     &              + TZ(K,2,2,EL_IN)*PHI_EDGE(K,GP_IN,LED_IN,pa)
               TZ_X_EX = TZ_X_EX 
     &              + TZ(K,1,1,EL_EX)*PHI_EDGE(K,GP_EX,LED_EX,pa)
               TZ_Y_EX = TZ_Y_EX 
     &              + TZ(K,2,2,EL_EX)*PHI_EDGE(K,GP_EX,LED_EX,pa)

#endif

            ENDDO

C.....Compute the numerical flux
            
            CALL NUMERICAL_FLUX(IT,test_el)

C......... dummy
            F_HAT_O  = F_HAT

#ifdef TRACE
            I_HAT_O  = I_HAT
#endif

#ifdef CHEM
            I_HAT_O  = I_HAT
            J_HAT_O  = J_HAT
#endif

#ifdef DYNP
            K_HAT_O  = K_HAT
#endif

#ifdef SED_LAY
            bed_HAT_IN(:) = bed_HAT(:)
            bed_HAT_EX(:) = bed_HAT(:)
#endif

            G_HAT_IN = G_HAT
            H_HAT_IN = H_HAT
            G_HAT_EX = G_HAT
            H_HAT_EX = H_HAT

C.....Check if the flux is large enough to dry up the elements
C.....1.01D0 is a safty factor.

          IF ( (1.01*F_HAT*XLEN_EL_IN*MAX_BOA_DT(IRK).GE.MASS_EL_IN).OR.
     &           (1.01*F_HAT*XLEN_EL_EX*MAX_BOA_DT(IRK)*(-1.D0).GE.
     &           MASS_EL_EX)) THEN

C...........Put back the saved RHS

               DO K = 1,DOFS(EL)
                  RHS_ZE(K,EL_IN,IRK) = RHS_ZE_IN(K)
                  RHS_QX(K,EL_IN,IRK) = RHS_QX_IN(K)
                  RHS_QY(K,EL_IN,IRK) = RHS_QY_IN(K)

                  RHS_ZE(K,EL_EX,IRK) = RHS_ZE_EX(K)
                  RHS_QX(K,EL_EX,IRK) = RHS_QX_EX(K)
                  RHS_QY(K,EL_EX,IRK) = RHS_QY_EX(K)

#ifdef TRACE
                  RHS_iota(K,EL_IN,IRK) = RHS_iota_IN(K)
                  RHS_iota(K,EL_EX,IRK) = RHS_iota_EX(K)
#endif

#ifdef CHEM
                  RHS_iota(K,EL_IN,IRK) = RHS_iota_IN(K)
                  RHS_iota2(K,EL_IN,IRK) = RHS_iota2_IN(K)
                  RHS_iota(K,EL_EX,IRK) = RHS_iota_EX(K)
                  RHS_iota2(K,EL_EX,IRK) = RHS_iota2_EX(K)
#endif

#ifdef DYNP
                  RHS_dynP(K,EL_IN,IRK) = RHS_dynP_IN(K)
                  RHS_dynP(K,EL_EX,IRK) = RHS_dynP_EX(K)
#endif

#ifdef SED_LAY
                  do ll = 1,layers
                     RHS_bed(K,EL_IN,IRK,ll) = RHS_bed_IN(K,ll)
                     RHS_bed(K,EL_EX,IRK,ll) = RHS_bed_EX(K,ll)
                  enddo
#endif

               ENDDO

               GOTO 100
            ENDIF
            
C........Check to make sure mass flux is not coming from a dry element

            IF (abs(F_HAT).gt.1.d-12) THEN
               !
               IF (WDFLG(EL_IN).EQ.0) THEN
                  ! EL_IN is dry !
                  IF (F_HAT.GT.0) THEN
                     ! Flux going from the dry element (in)
                     ! On the wet side (ex): reflect boundary
                     Q_N_EXT = QX_EX*NX + QY_EX*NY
                     Q_T_EXT = QX_EX*TX + QY_EX*TY
                     Q_N_INT = -Q_N_EXT
                     Q_T_INT =  Q_T_EXT
                     DEN = 1.D0/(NX*TY - NY*TX)
                     QX_IN = ( TY*Q_N_INT - NY*Q_T_INT)*DEN
                     QY_IN = (-TX*Q_N_INT + NX*Q_T_INT)*DEN
                     ZE_IN = ZE_EX

#ifdef TRACE
                     iota_IN = iota_EX
#endif

#ifdef CHEM
                     iota_IN = iota_EX
                     iota2_IN = iota2_EX
#endif

#ifdef DYNP
                     dynP_IN = dynP_EX
#endif

#ifdef SED_LAY
                     bed_IN(:) = bed_EX(:)
#endif

                     HB_IN = HB_EX
                     SFAC_IN = SFAC_EX
                     CALL NUMERICAL_FLUX(IT,test_el)

                     F_HAT_O  = F_HAT
                     G_HAT_EX = G_HAT
                     H_HAT_EX = H_HAT

                     ! on the dry side (in): do nothing
                     G_HAT_IN = G_HAT
                     H_HAT_IN = H_HAT
                     G_HAT_IN = 0.D0 
                     H_HAT_IN = 0.D0 

#ifdef TRACE
                     I_HAT_O  = I_HAT
#endif

#ifdef CHEM
                     I_HAT_O  = I_HAT
                     J_HAT_O  = J_HAT
#endif

#ifdef DYNP
                     K_HAT_O  = K_HAT
#endif

#ifdef SED_LAY
                     bed_HAT_IN(:) = bed_HAT(:)
                     bed_HAT_EX(:) = bed_HAT(:)
#endif

                  ELSE
                     ! Flux coming from wet side (ex)
                     ! Leave fluxes on ex side intact and use zero gravtity flux for dry side      
                     NLEQG_TMP = NLEQG
                     NLEQG = 0.D0
                     G_TMP = G
                     G = 0.D0
                     CALL NUMERICAL_FLUX(IT,test_el)
                     NLEQG = NLEQG_TMP
                     G = G_TMP
                     G_HAT_IN = G_HAT
                     H_HAT_IN = H_HAT
#ifdef SED_LAY
                     bed_HAT_IN(:) = bed_HAT(:)
#endif
                  ENDIF
                  
               ELSEIF (WDFLG(EL_EX).EQ.0) THEN
                  
                  ! EL_EX is dry
                  IF (F_HAT.LT.0) THEN
                     ! Flux comming from dry size (ex)
                     ! On the wet side (in): reflect boundary
                     Q_N_INT = QX_IN*NX + QY_IN*NY
                     Q_T_INT = QX_IN*TX + QY_IN*TY
                     Q_N_EXT = -Q_N_INT
                     Q_T_EXT =  Q_T_INT
                     DEN = 1.D0/(NX*TY - NY*TX)
                     QX_EX = ( TY*Q_N_EXT - NY*Q_T_EXT)*DEN
                     QY_EX = (-TX*Q_N_EXT + NX*Q_T_EXT)*DEN
                     ZE_EX = ZE_IN
                     HB_EX = HB_IN
                     SFAC_EX = SFAC_IN

#ifdef TRACE
                     iota_EX = iota_IN
#endif

#ifdef CHEM
                     iota_EX = iota_IN
                     iota2_EX = iota2_IN
#endif

#ifdef DYNP
                     dynP_EX = dynP_IN
#endif

#ifdef SED_LAY
                     bed_EX(:) = bed_IN(:)
#endif

                     CALL NUMERICAL_FLUX(IT,test_el)
                     F_HAT_O  = F_HAT

#ifdef TRACE
                     I_HAT_O  = I_HAT
#endif

#ifdef CHEM
                     I_HAT_O  = I_HAT
                     J_HAT_O  = J_HAT
#endif

#ifdef DYNP
                     K_HAT_O  = K_HAT
#endif

#ifdef SED_LAY
                     bed_HAT_IN(:) = bed_HAT(:)
                     bed_HAT_EX(:) = bed_HAT(:)
#endif
                     
                     G_HAT_IN = G_HAT
                     H_HAT_IN = H_HAT

                     ! zero out mentum flux on the dry side (ex)
                     G_HAT_EX = G_HAT
                     H_HAT_EX = H_HAT
                     G_HAT_EX = 0.D0 
                     H_HAT_EX = 0.D0 
                  ELSE
                     ! Flux comming from the wet side (in)
                     ! Leave fluxes on in side intact and use zero gravtity for flux for dry (ex) side     
                     NLEQG_TMP = NLEQG
                     NLEQG = 0.D0
                     G_TMP = G
                     G = 0.D0
                     CALL NUMERICAL_FLUX(IT,test_el)
                     NLEQG = NLEQG_TMP
                     G = G_TMP
                     G_HAT_EX = G_HAT
                     H_HAT_EX = H_HAT
#ifdef TRACE
                     I_HAT_O  = I_HAT
#endif
#ifdef CHEM
                     I_HAT_O  = I_HAT
                     J_HAT_O  = J_HAT
#endif

#ifdef SED_LAY
                     bed_HAT_EX(:) =  bed_HAT(:)
#endif
                  ENDIF
               ENDIF
            ENDIF

                                !Let us add some diffusive penalties
#ifdef WAVE_DIF
        F_HAT_O = F_HAT_O + 0.5D0*(HZ_X_IN*SFAC_IN+HZ_X_EX*SFAC_EX)*NX +
     &           0.5D0*(HZ_Y_IN+HZ_Y_EX)*NY
            !F_HAT_EX = F_HAT_EX + 0.5D0*(HZ_X_IN+HZ_X_EX)*SFAC_IN*NX + 
            !&           0.5D0*(HZ_Y_EX+HZ_Y_EX)*NY
#endif
            G_HAT_IN = G_HAT_IN + 0.5D0*(LZ_XX_IN*SFAC_IN + 
     &           LZ_XX_EX*SFAC_EX)*NX
     &           + 0.5D0*(LZ_XY_IN + LZ_XY_EX)*NY
cnd try adding the penalty term in the LDG method
      !$           +ESLM*.01D0*(QX_IN-QX_EX)
            H_HAT_IN = H_HAT_IN + 0.5D0*(LZ_YX_IN*SFAC_IN + 
     &           LZ_YX_EX*SFAC_EX)*NX
     &           + 0.5D0*(LZ_YY_IN + LZ_YY_EX)*NY
cnd try adding the penalty term in the LDG method
      !$           +ESLM*.01D0*(QY_IN-QY_EX)
            G_HAT_EX = G_HAT_EX + 0.5D0*(LZ_XX_IN*SFAC_IN
     &           + LZ_XX_EX*SFAC_EX)*NX
     &           + 0.5D0*(LZ_XY_IN + LZ_XY_EX)*NY
cnd try adding the penalty term in the LDG method
      !$           +ESLM*.01D0*(QX_IN-QX_EX)
            H_HAT_EX = H_HAT_EX + 0.5D0*(LZ_YX_IN*SFAC_IN
     &           + LZ_YX_EX*SFAC_EX)*NX
     &           + 0.5D0*(LZ_YY_IN + LZ_YY_EX)*NY
cnd try adding the penalty term in the LDG method
      !$           +ESLM*.01D0*(QY_IN-QY_EX) 

#ifdef TRACE

        I_HAT_O = I_HAT_O + 0.5D0*(TZ_X_IN*SFAC_IN+TZ_X_EX*SFAC_EX)*NX +
     &           0.5D0*(TZ_Y_IN+TZ_Y_EX)*NY
      !&           +ESLM*.01D0*(iota_IN-iota_EX)

#endif
     

                                !Let us add some sediment diffusion
#ifdef SED_LAY
            do ll=1,layers
               bed_HAT_IN(ll) = bed_HAT_IN(ll) + 0.5D0 * ( (MZ_X_IN(ll)+
     &             MZ_X_EX(ll))*SFAC_IN*NX+(MZ_Y_IN(ll)+MZ_Y_EX(ll))*NY)

               bed_HAT_EX(ll) = bed_HAT_EX(ll) + 0.5D0 * ( (MZ_X_IN(ll)+
     &             MZ_X_EX(ll))*SFAC_IN*NX+(MZ_Y_IN(ll)+MZ_Y_EX(ll))*NY)
            enddo
#endif

C.....Compute the edge integral

            DO K = 1,DOFS(EL)
               
               W_IN = EDFAC_IN*EDGEQ(K,GP_IN,LED_IN,pa)
               W_EX = EDFAC_EX*EDGEQ(K,GP_EX,LED_EX,pa)

               RHS_ZE(K,EL_IN,IRK) = RHS_ZE(K,EL_IN,IRK) - W_IN*F_HAT_O
               RHS_QX(K,EL_IN,IRK) = RHS_QX(K,EL_IN,IRK) - W_IN*G_HAT_IN
               RHS_QY(K,EL_IN,IRK) = RHS_QY(K,EL_IN,IRK) - W_IN*H_HAT_IN

               RHS_ZE(K,EL_EX,IRK) = RHS_ZE(K,EL_EX,IRK) + W_EX*F_HAT_O
               RHS_QX(K,EL_EX,IRK) = RHS_QX(K,EL_EX,IRK) + W_EX*G_HAT_EX
               RHS_QY(K,EL_EX,IRK) = RHS_QY(K,EL_EX,IRK) + W_EX*H_HAT_EX

#ifdef TRACE

            RHS_iota(K,EL_IN,IRK) = RHS_iota(K,EL_IN,IRK) - W_IN*I_HAT_O
            RHS_iota(K,EL_EX,IRK) = RHS_iota(K,EL_EX,IRK) + W_EX*I_HAT_O
#endif

#ifdef CHEM
          RHS_iota(K,EL_IN,IRK)  = RHS_iota(K,EL_IN,IRK)  - W_IN*I_HAT_O
          RHS_iota2(K,EL_IN,IRK) = RHS_iota2(K,EL_IN,IRK) - W_IN*J_HAT_O
          RHS_iota(K,EL_EX,IRK)  = RHS_iota(K,EL_EX,IRK)  + W_EX*I_HAT_O
          RHS_iota2(K,EL_EX,IRK) = RHS_iota2(K,EL_EX,IRK) + W_EX*J_HAT_O
#endif

#ifdef DYNP
           RHS_dynP(K,EL_IN,IRK) = RHS_dynP(K,EL_IN,IRK) - W_IN*K_HAT_O
           RHS_dynP(K,EL_EX,IRK) = RHS_dynP(K,EL_EX,IRK) + W_EX*K_HAT_O
#endif


#ifdef SED_LAY
               do ll = 1,layers
                  RHS_bed(K,EL_IN,IRK,ll) = RHS_bed(K,EL_IN,IRK,ll) -
     & W_IN*bed_HAT_IN(ll)
                  RHS_bed(K,EL_EX,IRK,ll) = RHS_bed(K,EL_EX,IRK,ll) +
     & W_EX*bed_HAT_EX(ll)
               enddo
#endif

            ENDDO

         ENDDO

         GOTO 1000

 100     CONTINUE

C--------------------------------------------------------------------
C     Compute reflection flux
C--------------------------------------------------------------------

C.....Set the components for the tangential vector to the edge

         TX = -NY
         TY =  NX

C.....Compute ZE, QX, QY, bed, and HB at each edge Gauss quadrature point

         DO I = 1,NEGP(pa)

            GP_IN = I
            GP_EX = NEGP(pa) - I + 1

C     ----------------------------- EL_IN ------------------------------

            ZE_IN = ZE(1,EL_IN,IRK)
            QX_IN = QX(1,EL_IN,IRK)
            QY_IN = QY(1,EL_IN,IRK)
            HB_IN = BATHED(GP_IN,LED_IN,EL_IN,pa)

            SFAC_IN = SFACED(GP_IN,LED_IN,EL_IN,pa)

#ifdef TRACE
            iota_IN = iota(1,EL_IN,IRK)
#endif

#ifdef CHEM
            iota_IN = iota(1,EL_IN,IRK)
            iota2_IN = iota2(1,EL_IN,IRK)
#endif

#ifdef DYNP
            dynP_IN = dynP(1,EL_IN,IRK)
#endif

#ifdef SED_LAY
            HB(1,EL_IN,irk) = 0.D0
            do ll=1,layers
               bed_IN(ll) = bed(1,EL_IN,irk,ll)
               HB(1,EL_IN,irk) = HB(1,EL_IN,irk) + bed(1,EL_IN,irk,ll)

               !sediment diffusion
               MZ_X_IN(ll) = MZ(1,1,ll,EL_IN)
               MZ_Y_IN(ll) = MZ(1,2,ll,EL_IN)  
            enddo
            HB_IN = HB(1,EL_IN,irk)
#endif

            !diffusion
#ifdef WAVE_DIF
            HZ_X_IN = HZ(1,1,1,EL_IN)
            HZ_Y_IN = HZ(1,2,2,EL_IN)
#endif
            LZ_XX_IN = LZ(1,1,1,EL_IN)
            LZ_XY_IN = LZ(1,1,2,EL_IN)
            LZ_YX_IN = LZ(1,2,1,EL_IN)
            LZ_YY_IN = LZ(1,2,2,EL_IN)

#ifdef TRACE

            TZ_X_IN = TZ(1,1,1,EL_IN)
            TZ_Y_IN = TZ(1,2,2,EL_IN)

#endif

            DO K = 2,DOFS(EL)

             ZE_IN = ZE_IN + ZE(K,EL_IN,IRK)*PHI_EDGE(K,GP_IN,LED_IN,pa)
             QX_IN = QX_IN + QX(K,EL_IN,IRK)*PHI_EDGE(K,GP_IN,LED_IN,pa)
             QY_IN = QY_IN + QY(K,EL_IN,IRK)*PHI_EDGE(K,GP_IN,LED_IN,pa)

#ifdef TRACE
       iota_IN = iota_IN + iota(K,EL_IN,IRK)*PHI_EDGE(K,GP_IN,LED_IN,PA)
#endif

#ifdef CHEM
      iota_IN = iota_IN  + iota(K,EL_IN,IRK)*PHI_EDGE(K,GP_IN,LED_IN,PA)
      iota2_IN =iota2_IN+ iota2(K,EL_IN,IRK)*PHI_EDGE(K,GP_IN,LED_IN,PA)
#endif

#ifdef DYNP
      dynP_IN = dynP_IN + dynP(K,EL_IN,IRK)*PHI_EDGE(K,GP_IN,LED_IN,PA)
#endif

#ifdef SED_LAY
               do ll=1,layers
                  bed_IN(ll) = bed_IN(ll) +
     & bed(k,EL_IN,irk,ll)*PHI_EDGE(K,GP_IN,LED_IN,pa)
                  HB_IN = HB_IN +
     & bed(k,EL_IN,irk,ll)*PHI_EDGE(K,GP_IN,LED_IN,pa)
                  
                  MZ_X_IN(ll) = MZ_X_IN(ll) +
     & MZ(K,1,ll,EL_IN)*PHI_EDGE(K,GP_IN,LED_IN,pa)
                  MZ_Y_IN(ll) = MZ_Y_IN(ll) +
     & MZ(K,2,ll,EL_IN)*PHI_EDGE(K,GP_IN,LED_IN,pa)
               enddo
#endif
                                ! LDG terms
#ifdef WAVE_DIF
               HZ_X_IN = HZ_X_IN
     &              + HZ(K,1,1,EL_IN)*PHI_EDGE(K,GP_IN,LED_IN,pa)
               HZ_Y_IN = HZ_Y_IN
     &              + HZ(K,2,2,EL_IN)*PHI_EDGE(K,GP_IN,LED_IN,pa)
#endif
               LZ_XX_IN = LZ_XX_IN
     &              + LZ(K,1,1,EL_IN)*PHI_EDGE(K,GP_IN,LED_IN,pa)
               LZ_XY_IN = LZ_XY_IN
     &              + LZ(K,1,2,EL_IN)*PHI_EDGE(K,GP_IN,LED_IN,pa)
               LZ_YX_IN = LZ_YX_IN
     &              + LZ(K,2,1,EL_IN)*PHI_EDGE(K,GP_IN,LED_IN,pa)
               LZ_YY_IN = LZ_YY_IN
     &              + LZ(K,2,2,EL_IN)*PHI_EDGE(K,GP_IN,LED_IN,pa)

#ifdef TRACE
               TZ_X_IN = TZ_X_IN
     &              + TZ(K,1,1,EL_IN)*PHI_EDGE(K,GP_IN,LED_IN,pa)
               TZ_Y_IN = TZ_Y_IN
     &              + TZ(K,2,2,EL_IN)*PHI_EDGE(K,GP_IN,LED_IN,pa)
#endif

            ENDDO

C.....Compute the velocity in the normal and tangental direction

            Q_N_INT = QX_IN*NX + QY_IN*NY
            Q_T_INT = QX_IN*TX + QY_IN*TY

C.....Reflect the velocity in the normal direction

            Q_N_EXT = -Q_N_INT
            Q_T_EXT =  Q_T_INT

C.....Compute the x and y components of the external state flow

            DEN = 1.D0/(NX*TY - NY*TX)
            QX_EX = ( TY*Q_N_EXT - NY*Q_T_EXT)*DEN
            QY_EX = (-TX*Q_N_EXT + NX*Q_T_EXT)*DEN

            ZE_EX = ZE_IN
            HB_EX = HB_IN
            SFAC_EX = SFAC_IN

#ifdef TRACE
            iota_EX = iota_IN
#endif

#ifdef CHEM
            iota_EX = iota_IN
            iota2_EX = iota2_IN
#endif

#ifdef DYNP
            dynP_EX = dynP_IN
#endif

#ifdef SED_LAY
            bed_EX(:) = bed_IN(:)
#endif

C.....Compute the numerical flux

C     CALL NUMERICAL_FLUX(IT)
            IF(WDFLG(EL_IN).EQ.0) THEN
               NLEQG_TMP = NLEQG
               NLEQG = 0.D0
               G_TMP = G
               G = 0.D0
               CALL NUMERICAL_FLUX(IT,test_el)
               NLEQG = NLEQG_TMP
               G = G_TMP
            ELSE
               CALL NUMERICAL_FLUX(IT,test_el)
            ENDIF
            F_HAT = 0.D0        ! Ensure zero mass flux

C.....Add LDG terms for eddy viscosity

#ifdef WAVE_DIF
            F_HAT = F_HAT + HZ_X_IN*NX*SFAC_IN
     &           + HZ_Y_IN*NY
#endif
            G_HAT = G_HAT + LZ_XX_IN*NX*SFAC_IN
     &           + LZ_XY_IN*NY
            H_HAT = H_HAT + LZ_YX_IN*NX*SFAC_IN
     &           + LZ_YY_IN*NY

#ifdef TRACE

            I_HAT = I_HAT + TZ_X_IN*NX*SFAC_IN
     &           + TZ_Y_IN*NY

#endif

C.....Add LDG terms for sediment diffusion

#ifdef SED_LAY
            do ll=1,layers
               bed_HAT(ll) = bed_HAT(ll) + MZ_X_IN(ll)*SFAC_IN*NX 
     &              + MZ_Y_IN(ll)*NY
            enddo
#endif

C.....Compute the edge integral
            DO K = 1,DOFS(EL)

               W_IN = EDFAC_IN*EDGEQ(K,GP_IN,LED_IN,pa)

               RHS_ZE(K,EL_IN,IRK) = RHS_ZE(K,EL_IN,IRK) - W_IN*F_HAT
               RHS_QX(K,EL_IN,IRK) = RHS_QX(K,EL_IN,IRK) - W_IN*G_HAT
               RHS_QY(K,EL_IN,IRK) = RHS_QY(K,EL_IN,IRK) - W_IN*H_HAT

#ifdef TRACE
             RHS_iota(K,EL_IN,IRK) = RHS_iota(K,EL_IN,IRK) - W_IN*I_HAT
#endif

#ifdef CHEM
            RHS_iota(K,EL_IN,IRK)  = RHS_iota(K,EL_IN,IRK)  - W_IN*I_HAT
            RHS_iota2(K,EL_IN,IRK) = RHS_iota2(K,EL_IN,IRK) - W_IN*J_HAT
#endif

#ifdef DYNP
              RHS_dynP(K,EL_IN,IRK) = RHS_dynP(K,EL_IN,IRK) - W_IN*K_HAT
#endif

#ifdef SED_LAY
               do ll = 1,layers
                  RHS_bed(K,EL_IN,IRK,ll) = RHS_bed(K,EL_IN,IRK,ll) -
     & W_IN*bed_HAT(ll)
               enddo
#endif

            ENDDO

C     ----------------------------- EL_EX ------------------------------

            ZE_EX = ZE(1,EL_EX,IRK)
            QX_EX = QX(1,EL_EX,IRK)
            QY_EX = QY(1,EL_EX,IRK)
            HB_EX = BATHED(GP_EX,LED_EX,EL_EX,pa)

            SFAC_EX = SFACED(GP_EX,LED_EX,EL_EX,pa)

#ifdef TRACE
            iota_EX = iota(1,EL_EX,IRK)
#endif

#ifdef CHEM
            iota_EX = iota(1,EL_EX,IRK)
            iota2_EX = iota2(1,EL_EX,IRK)
#endif

#ifdef DYNP
            dynP_EX = dynP(1,EL_EX,IRK)
#endif

#ifdef SED_LAY
            HB(1,EL_EX,irk) = 0.D0
            do ll=1,layers
               bed_EX(ll) = bed(1,EL_EX,irk,ll)
               HB(1,EL_EX,irk) = HB(1,EL_EX,irk) + bed(1,EL_EX,irk,ll)
               !sediment diffusion
               MZ_X_EX(ll) = MZ(1,1,ll,EL_EX)
               MZ_Y_EX(ll) = MZ(1,2,ll,EL_EX)  
            enddo
            HB_EX = HB(1,EL_EX,irk)
#endif

            !LDG
#ifdef WAVE_DIF
            HZ_X_EX = HZ(1,1,1,EL_EX)
            HZ_Y_EX = HZ(1,2,2,EL_EX)
#endif
            LZ_XX_EX = LZ(1,1,1,EL_EX)
            LZ_XY_EX = LZ(1,1,2,EL_EX)
            LZ_YX_EX = LZ(1,2,1,EL_EX)
            LZ_YY_EX = LZ(1,2,2,EL_EX)

#ifdef TRACE

            TZ_X_EX = TZ(1,1,1,EL_EX)
            TZ_Y_EX = TZ(1,2,2,EL_EX)

#endif

            DO K = 2,DOFS(EL)

             ZE_EX = ZE_EX + ZE(K,EL_EX,IRK)*PHI_EDGE(K,GP_EX,LED_EX,pa)
             QX_EX = QX_EX + QX(K,EL_EX,IRK)*PHI_EDGE(K,GP_EX,LED_EX,pa)
             QY_EX = QY_EX + QY(K,EL_EX,IRK)*PHI_EDGE(K,GP_EX,LED_EX,pa)

#ifdef TRACE
       iota_EX = iota_EX + iota(K,EL_EX,IRK)*PHI_EDGE(K,GP_EX,LED_EX,pa)
#endif

#ifdef CHEM
      iota_EX = iota_EX  + iota(K,EL_EX,IRK)*PHI_EDGE(K,GP_EX,LED_EX,pa)
      iota2_EX = iota2_EX+iota2(K,EL_EX,IRK)*PHI_EDGE(K,GP_EX,LED_EX,pa)
#endif

#ifdef DYNP
      dynP_EX = dynP_EX + dynP(K,EL_EX,IRK)*PHI_EDGE(K,GP_EX,LED_EX,pa)
#endif

#ifdef SED_LAY
               do ll=1,layers
                  bed_EX(ll) = bed_EX(ll) +
     & bed(k,EL_EX,irk,ll)*PHI_EDGE(K,GP_EX,LED_EX,pa)
                  HB_EX = HB_EX +
     & bed(k,EL_EX,irk,ll)*PHI_EDGE(K,GP_EX,LED_EX,pa)
                  
                  MZ_X_EX(ll) = MZ_X_EX(ll) +
     & MZ(K,1,ll,EL_EX)*PHI_EDGE(K,GP_EX,LED_EX,pa)
                  MZ_Y_EX(ll) = MZ_Y_EX(ll) +
     & MZ(K,2,ll,EL_EX)*PHI_EDGE(K,GP_EX,LED_EX,pa)
               enddo
#endif

                                ! LDG terms
#ifdef WAVE_DIF
               HZ_X_EX = HZ_X_EX
     &              + HZ(K,1,1,EL_EX)*PHI_EDGE(K,GP_EX,LED_EX,pa)
               HZ_Y_EX = HZ_Y_EX
     &              + HZ(K,2,2,EL_EX)*PHI_EDGE(K,GP_EX,LED_EX,pa)
#endif
               LZ_XX_EX = LZ_XX_EX
     &              + LZ(K,1,1,EL_EX)*PHI_EDGE(K,GP_EX,LED_EX,pa)
               LZ_XY_EX = LZ_XY_EX
     &              + LZ(K,1,2,EL_EX)*PHI_EDGE(K,GP_EX,LED_EX,pa)
               LZ_YX_EX = LZ_YX_EX
     &              + LZ(K,2,1,EL_EX)*PHI_EDGE(K,GP_EX,LED_EX,pa)
               LZ_YY_EX = LZ_YY_EX
     &              + LZ(K,2,2,EL_EX)*PHI_EDGE(K,GP_EX,LED_EX,pa)

#ifdef TRACE

               TZ_X_EX = TZ_X_EX
     &              + TZ(K,1,1,EL_EX)*PHI_EDGE(K,GP_EX,LED_EX,pa)
               TZ_Y_EX = TZ_Y_EX
     &              + TZ(K,2,2,EL_EX)*PHI_EDGE(K,GP_EX,LED_EX,pa)

#endif

                                ! LDG terms for sediment diffusion               

            ENDDO

C.....Compute the velocity in the normal and tangental direction

            Q_N_EXT = QX_EX*NX + QY_EX*NY
            Q_T_EXT = QX_EX*TX + QY_EX*TY

C.....Reflect the velocity in the normal direction

            Q_N_INT = -Q_N_EXT
            Q_T_INT =  Q_T_EXT

C.....Compute the x and y components of the external state flow

            DEN = 1.D0/(NX*TY - NY*TX)
            QX_IN = ( TY*Q_N_INT - NY*Q_T_INT)*DEN
            QY_IN = (-TX*Q_N_INT + NX*Q_T_INT)*DEN

            ZE_IN = ZE_EX
            HB_IN = HB_EX

            SFAC_IN = SFAC_EX

#ifdef TRACE
            iota_IN = iota_EX
#endif

#ifdef CHEM
            iota_IN = iota_EX
            iota2_IN = iota2_EX
#endif

#ifdef DYNP
            dynP_IN = dynP_EX
#endif

#ifdef SED_LAY
            bed_IN(:) = bed_EX(:)
#endif

C.....Compute the numerical flux

C     CALL NUMERICAL_FLUX(IT,L)
            IF(WDFLG(EL_EX).EQ.0) THEN
               NLEQG_TMP = NLEQG
               NLEQG = 0.D0
               G_TMP = G
               G = 0.D0
               CALL NUMERICAL_FLUX(IT,test_el)
               NLEQG = NLEQG_TMP
               G = G_TMP
            ELSE
               CALL NUMERICAL_FLUX(IT,test_el)
            ENDIF
            F_HAT = 0.D0        ! Ensure zero mass flux

C.....Add LDG terms for eddy viscosity
#ifdef WAVE_DIF
            F_HAT = F_HAT + HZ_X_EX*NX*SFAC_EX
     &           + HZ_Y_EX*NY
#endif
            G_HAT = G_HAT + LZ_XX_EX*NX*SFAC_EX
     &           + LZ_XY_EX*NY
            H_HAT = H_HAT + LZ_YX_EX*NX*SFAC_EX
     &           + LZ_YY_EX*NY

#ifdef TRACE

            I_HAT = I_HAT + TZ_X_EX*NX*SFAC_EX
     &           + TZ_Y_EX*NY

#endif

C.....Add LDG terms for sediment diffusion

#ifdef SED_LAY
            do ll=1,layers
               bed_HAT(ll) = bed_HAT(ll) + MZ_X_EX(ll)*SFAC_EX*NX 
     &              + MZ_Y_EX(ll)*NY
            enddo
#endif

C.....Compute the edge integral
            DO K = 1,DOFS(EL)

               W_EX = EDFAC_EX*EDGEQ(K,GP_EX,LED_EX,pa)

               RHS_ZE(K,EL_EX,IRK) = RHS_ZE(K,EL_EX,IRK) + W_EX*F_HAT
               RHS_QX(K,EL_EX,IRK) = RHS_QX(K,EL_EX,IRK) + W_EX*G_HAT
               RHS_QY(K,EL_EX,IRK) = RHS_QY(K,EL_EX,IRK) + W_EX*H_HAT

#ifdef TRACE
              RHS_iota(K,EL_EX,IRK) = RHS_iota(K,EL_EX,IRK) + W_EX*I_HAT
#endif

#ifdef CHEM
            RHS_iota(K,EL_EX,IRK)  = RHS_iota(K,EL_EX,IRK)  + W_EX*I_HAT
            RHS_iota2(K,EL_EX,IRK) = RHS_iota2(K,EL_EX,IRK) + W_EX*J_HAT
#endif

#ifdef DYNP
             RHS_dynP(K,EL_EX,IRK) = RHS_dynP(K,EL_EX,IRK) + W_EX*K_HAT
#endif

#ifdef SED_LAY
               do ll = 1,layers
                  RHS_bed(K,EL_EX,IRK,ll) = RHS_bed(K,EL_EX,IRK,ll) +
     & W_EX*bed_HAT(ll)
               enddo
#endif

            ENDDO

         ENDDO
         
 1000 CONTINUE

      RETURN
      END SUBROUTINE
