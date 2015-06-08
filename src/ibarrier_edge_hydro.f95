C***********************************************************************
C     
C     SUBROUTINE IBARRIER_EDGE_HYDRO()
C     
C     This subroutine does the following:
C     
C     1.  Calculates the values of the necessary variables at the edge
C     gauss points for INTERNAL BARRIER edges
C     2.  Calls the appropriate subroutine to compute the flux at
C     these points.
C     3.  Computes the boundary integrals.
C     
C     Written by Ethan Kubatko (08-11-2008)
C     
C     01-10-2011 - cem - adapted for p_enrichment and multicomponent
C     01-20-2012 - cem - multiphase added, we do not (natively)
C                        transport sediment across weirs
C     
C***********************************************************************

      SUBROUTINE IBARRIER_EDGE_HYDRO(IT)
      
C.....Use appropriate modules

      USE GLOBAL
      USE DG
      USE SIZES, ONLY:  MYPROC, layers

      IMPLICIT NONE

C.....Declare local variables

      INTEGER GEDB, GEDF, LEDB, LEDF, ELB, ELF, GPB, GPF, L, WEIR_FLOW
      INTEGER IT, NB1, NB2, NF1, NF2,i,k,ll
      REAL(SZ) ZEB, QXB, QYB, HBB, ZEF, QXF, QYF, HBF
      REAL(SZ) SFACB,SFACF
      Real(SZ) iotaB, iota2B, iotaF, iota2F, dynPB, dynPF
      REAL(SZ) FB_HAT, GB_HAT, HB_HAT, FF_HAT, GF_HAT, HF_HAT
      REAL(SZ) ib_hat,if_hat,jb_hat,jf_hat, kf_hat, kb_hat
      REAL(SZ) ABOVEB, ABOVEF
      REAL(SZ) QNIB(2)
      REAL(SZ) NLEQG_TMP, G_TMP
      REAL(SZ) SUBSUPB, SUBSUPF, WEGPB, WEGPF
      REAL(SZ) NXB, NYB, NXF, NYF, TXB, TYB, TXF, TYF
      REAL(SZ) QB_N_INT, QF_N_INT, QB_N_EXT, QF_N_EXT
      REAL(SZ) QB_T_INT, QF_T_INT, QB_T_EXT, QF_T_EXT

C.....Loop over the internal barrier segments (note: an internal barrier
C.....segment consists of two internal barrier edges -- a "front" edge &
C.....a "back" edge.
      test_el = 0
      DO 1000 L = 1,NIBSEG
         
C.......Obtain the global and local edges of the back and front sides

         GEDB = NIBSEGN(1,L)
         GEDF = NIBSEGN(2,L)

         if (gedf.eq.0) go to 1000
         
         LEDB = NEDSD(1,GEDB)
         LEDF = NEDSD(1,GEDF)
         
C.......Obtain the elements of the back and front sides

         ELB = NEDEL(1,GEDB)
         ELF = NEDEL(1,GEDF)

         if (DOFS(ELB).LT.DOFS(ELF)) then
            EL = ELF
         endif

         pa = PDG_EL(EL)

#ifdef P0         
         if (pa.eq.0) then
            pa = 1
         endif
#endif
         
         NB1 = NEDNO(1,GEDB)
         NB2 = NEDNO(2,GEDB)
         
         NF1 = NEDNO(1,GEDF)
         NF2 = NEDNO(2,GEDF)
         IF((WDFLG(ELB).EQ.0).AND.(WDFLG(ELF).EQ.0)) GOTO 1000

         test_el = test_el+1

C.......Retrieve the components of the normal vector to the edge

         NXB = COSNX(GEDB)
         NYB = SINNX(GEDB)
         
         NXF = COSNX(GEDF)
         NYF = SINNX(GEDF)
         
C.....Set the components for the tangential vector to the edge

         TXB = -NYB
         TYB =  NXB

         TXF = -NYF
         TYF =  NXF
         
C.......Compute the variables at the quadrature points

         DO I = 1,NEGP(pa)
            
            GPB = I
            GPF = NEGP(pa) - I + 1
            
C.........Obtain the height of the barrier at the quadrature point

            ZEB = ZE(1,ELB,IRK)
            QXB = QX(1,ELB,IRK)
            QYB = QY(1,ELB,IRK)
            HBB = BATHED(GPB,LEDB,ELB,pa)
            
            ZEF = ZE(1,ELF,IRK)
            QXF = QX(1,ELF,IRK)
            QYF = QY(1,ELF,IRK)
            HBF = BATHED(GPF,LEDF,ELF,pa)

#ifdef TRACE
            iotaB = iota(1,ELB,IRK)
            iotaF = iota(1,ELF,IRK)
#endif

#ifdef CHEM
            iotaB = iota(1,ELB,IRK)
            iota2B = iota2(1,ELB,IRK)
            iotaF = iota(1,ELF,IRK)
            iota2F = iota2(1,ELF,IRK)
#endif

#ifdef DYNP
            dynPB = dynP(1,ELB,IRK)
            dynPF = dynP(1,ELF,IRK)
#endif

#ifdef SED_LAY !This does not affect much here, hb and bed do not change at weirs
            HB(1,ELB,irk) = 0.D0
            do ll=1,layers
               HB(1,ELB,irk) = HB(1,ELB,irk) + bed(1,ELB,irk,ll)
            enddo
            HBB = HB(1,ELB,irk)
            HB(1,ELF,irk) = 0.D0
            do ll=1,layers
               HB(1,ELF,irk) = HB(1,ELF,irk) + bed(1,ELF,irk,ll)
            enddo
            HBF = HB(1,ELF,irk)
#endif



            DO K = 2,dofs(EL)
               ZEB = ZEB + ZE(K,ELB,IRK)*PHI_EDGE(K,GPB,LEDB,pa)
               QXB = QXB + QX(K,ELB,IRK)*PHI_EDGE(K,GPB,LEDB,pa)
               QYB = QYB + QY(K,ELB,IRK)*PHI_EDGE(K,GPB,LEDB,pa)
                                !HBB = HBB + HB(K,ELB,1  )*PHI_EDGE(K,GPB,LEDB,pa)

               ZEF = ZEF + ZE(K,ELF,IRK)*PHI_EDGE(K,GPF,LEDF,pa)
               QXF = QXF + QX(K,ELF,IRK)*PHI_EDGE(K,GPF,LEDF,pa)
               QYF = QYF + QY(K,ELF,IRK)*PHI_EDGE(K,GPF,LEDF,pa)
                                !HBF = HBF + HB(K,ELF,1  )*PHI_EDGE(K,GPF,LEDF,pa)
#ifdef TRACE
               iotaB = iotaB + iota(K,ELB,IRK)*PHI_EDGE(K,GPB,LEDB,pa)
               iotaF = iotaF + iota(K,ELF,IRK)*PHI_EDGE(K,GPF,LEDF,pa)
#endif

#ifdef CHEM
               iotaB = iotaB + iota(K,ELB,IRK)*PHI_EDGE(K,GPB,LEDB,pa)
               iotaF = iotaF + iota(K,ELF,IRK)*PHI_EDGE(K,GPF,LEDF,pa)
               iota2B = iota2B + iota2(K,ELB,IRK)*PHI_EDGE(K,GPB,LEDB,pa)
               iota2B = iota2F + iota2(K,ELF,IRK)*PHI_EDGE(K,GPF,LEDF,pa)
#endif

#ifdef DYNP
               dynPB = dynPB + dynP(K,ELB,IRK)*PHI_EDGE(K,GPB,LEDB,pa)
               dynPF = dynPF + dynP(K,ELF,IRK)*PHI_EDGE(K,GPF,LEDF,pa)
#endif

            ENDDO


            SFACB = SFACED(I,LEDB,ELB,pa)
            SFACF = SFACED(I,LEDF,ELF,pa)
            
            ABOVEB = 0
            ABOVEF = 0 
            IF (WDFLG(ELB).EQ.1) ABOVEB = ZEB - IBHT(L)
            IF (WDFLG(ELF).EQ.1) ABOVEF = ZEF - IBHT(L)
            
C.........Case 1:  Water is below barrier on both sides
C     ---------------------------------------------
            
            IF ((ABOVEF.LE.BARMIN).AND.(ABOVEB.LE.BARMIN)) THEN
               QB_N_INT = -(QXB*NXB + QYB*NYB)
               QF_N_INT = -(QXF*NXF + QYF*NYF)
               QB_T_INT = QXB*TXB + QYB*TYB
               QF_T_INT = QXF*TXF + QYF*TYF
               WEIR_FLOW = 0
               GOTO 100
               
C.........Case 2:  Water is above on both sides and equal (within tol)
C     ------------------------------------------------------------
               
            ELSEIF (ABS(ABOVEF-ABOVEB).LT.BARMIN) THEN
               QB_N_INT = -(QXB*NXB + QYB*NYB)
               QF_N_INT = -(QXF*NXF + QYF*NYF)
               QB_T_INT = QXB*TXB + QYB*TYB
               QF_T_INT = QXF*TXF + QYF*TYF
               WEIR_FLOW = 0
               GOTO 100
            ENDIF

            SUBSUPF = 2.D0*ABOVEF/3.D0
            SUBSUPB = 2.D0*ABOVEB/3.D0

C.........Case 3: Overtopping of barrier with water higher on front side
C     Flow from front to back side
C     -------------------------------------------------------------

            IF ((ABOVEF.GT.ABOVEB).AND.(ABOVEF.GT.BARMIN)) THEN
               
               WEIR_FLOW = 1
               
C...........Case 3a) Subcritical flow
C     -------------------------

               IF (ABOVEB.GT.SUBSUPF) THEN
                  QF_N_INT = RAMPDG*IBCFSB(L)*ABOVEB
     &                 *SQRT((2.D0*G*(ABOVEF-ABOVEB)))
                  QF_T_INT = 0.D0
                  
C...........Case 3b) Supercritical flow
C     ---------------------------

               ELSE
                  QF_N_INT = RAMPDG*IBCFSP(L)*SUBSUPF*SQRT(SUBSUPF*G)
                  QF_T_INT = 0.D0
               ENDIF
               GOTO 100

            ENDIF
            
C.........Case 4: Overtopping of barrier with water higher on back side
C     Flow from back to front side
C     --------------------------------------------------------------

            IF ((ABOVEB.GT.ABOVEF).AND.(ABOVEB.GT.BARMIN)) THEN
               
               WEIR_FLOW = -1
               
C...........Case 4a) Subcritical flow
C     -------------------------

               IF (ABOVEF.GT.SUBSUPB) THEN
                  QB_N_INT = RAMPDG*IBCFSB(L)*ABOVEF
     &                 *SQRT((2.D0*G*(ABOVEB-ABOVEF)))
                  QB_T_INT = 0.D0
                  
C...........Case 4b) Supercritical flow
C     ---------------------------
                  
               ELSE
                  QB_N_INT = RAMPDG*IBCFSP(L)*SUBSUPB*SQRT(SUBSUPB*G)
                  QB_T_INT = 0.D0
               ENDIF
               GOTO 100
            ENDIF
            
 100        CONTINUE

            IF (WEIR_FLOW.LE.0) THEN

C...........Compute the numerical flux for the back side edge

               ZE_IN = ZEB
               QX_IN = QXB
               QY_IN = QYB
               HB_IN = HBB

               SFAC_IN = SFACB

#ifdef TRACE
               iota_IN = iotaB
#endif

#ifdef CHEM
               iota_IN = iotaB
               iota2_IN = iota2B
#endif

#ifdef DYNP
               dynP_IN = dynPB
#endif
               
               ZE_EX = ZEB
               HB_EX = HBB

               SFAC_EX = SFACB

#ifdef TRACE
               iota_EX = iotaB
#endif

#ifdef CHEM
               iota_EX = iotaB
               iota2_EX = iota2B
#endif

#ifdef DYNP
               dynP_EX = dynPB
#endif
               
               NX = NXB
               NY = NYB

C...........Reflect the velocity in the normal direction

               Q_N_EXT = QB_N_INT
               Q_T_EXT = QB_T_INT

C...........Compute the x and y components of the external state flow

               QX_EX = ( TYB*Q_N_EXT - NYB*Q_T_EXT)/(NXB*TYB - NYB*TXB)
               QY_EX = (-TXB*Q_N_EXT + NXB*Q_T_EXT)/(NXB*TYB - NYB*TXB)

C...........Compute the numerical flux
               
               CALL NUMERICAL_FLUX(IT)
               FB_HAT = F_HAT
               GB_HAT = G_HAT
               HB_HAT = H_HAT

#ifdef TRACE
               IB_HAT = I_HAT
#endif

#ifdef CHEM
               IB_HAT = I_HAT
               JB_HAT = J_HAT
#endif

#ifdef DYNP
               KB_HAT = K_HAT
#endif

               IF (WEIR_FLOW.LT.0) THEN
                  ZE_IN = ZEF
                  QX_IN = QXF
                  QY_IN = QYF
                  HB_IN = HBF
                  ZE_EX = ZEF
                  HB_EX = HBF

                  SFAC_IN = SFACF
                  SFAC_EX = SFACF


                  NX = NXF
                  NY = NYF
                  
#ifdef TRACE
                  iota_IN = iotaF
                  iota_EX = iotaF
#endif

#ifdef CHEM
                  iota_IN = iotaF
                  iota2_IN = iota2F
                  iota_EX = iotaF
                  iota2_EX = iota2F
#endif

#ifdef DYNP
                  dynP_IN = iotaF
                  dynP_EX = iotaF
#endif

                  IF (WDFLG(ELF).EQ.0) THEN
                     NLEQG_TMP = NLEQG
                     NLEQG = 0.D0
                     G_TMP = G
                     G = 0.D0
                  ENDIF
                  CALL NUMERICAL_FLUX(IT)
                  FF_HAT = F_HAT
                  GF_HAT = G_HAT
                  HF_HAT = H_HAT

#ifdef TRACE
                  IF_HAT = I_HAT
#endif

#ifdef CHEM
                  IF_HAT = I_HAT
                  JF_HAT = J_HAT
#endif

#ifdef DYNP
                  KF_HAT = K_HAT
#endif

                  IF (WDFLG(ELF).EQ.0) THEN
                     NLEQG = NLEQG_TMP
                     G = G_TMP
                  ENDIF
                  GOTO 200
               ENDIF
            ENDIF

            IF (WEIR_FLOW.GE.0) THEN
               
C...........Compute the numerical flux for the front side edge

               ZE_IN = ZEF
               QX_IN = QXF
               QY_IN = QYF
               HB_IN = HBF

               ZE_EX = ZEF
               HB_EX = HBF

               SFAC_IN = SFACF
               SFAC_EX = SFACF

#ifdef TRACE
               iota_IN = iotaF
               iota_EX = iotaF
#endif

#ifdef CHEM
               iota_IN = iotaF
               iota_EX = iotaF
               iota2_IN = iota2F
               iota2_EX = iota2F
#endif

#ifdef DYNP
               dynP_IN = dynPF
               dynP_EX = dynPF
#endif

               NX = NXF
               NY = NYF
               
C...........Reflect the velocity in the normal direction

               Q_N_EXT = QF_N_INT
               Q_T_EXT = QF_T_INT

C...........Compute the x and y components of the external state flow

               QX_EX = ( TYF*Q_N_EXT - NYF*Q_T_EXT)/(NXF*TYF - NYF*TXF)
               QY_EX = (-TXF*Q_N_EXT + NXF*Q_T_EXT)/(NXF*TYF - NYF*TXF)

               CALL NUMERICAL_FLUX(IT)
               FF_HAT = F_HAT
               GF_HAT = G_HAT
               HF_HAT = H_HAT

#ifdef TRACE
               IF_HAT = I_HAT
#endif

#ifdef CHEM
               IF_HAT = I_HAT
               JF_HAT = J_HAT
#endif

#ifdef DYNP
               KF_HAT = K_HAT
#endif
               
               IF (WEIR_FLOW.GT.0) THEN
                  ZE_IN = ZEB
                  QX_IN = QXB
                  QY_IN = QYB
                  HB_IN = HBB
                  ZE_EX = ZEB
                  HB_EX = HBB

                  SFAC_IN = SFACB
                  SFAC_EX = SFACB

#ifdef TRACE
                  iota_IN = iotaB
                  iota_EX = iotaB
#endif

#ifdef CHEM
                  iota_IN = iotaB
                  iota2_IN = iota2B
                  iota_EX = iotaB
                  iota2_EX = iota2B
#endif

#ifdef DYNP
                  dynP_IN = dynPB
                  dynP_EX = dynPB
#endif

                  NX = NXB
                  NY = NYB
                  IF (WDFLG(ELB).EQ.0) THEN
                     NLEQG_TMP = NLEQG
                     NLEQG = 0.D0
                     G_TMP = G
                     G = 0.D0
                  ENDIF
                  CALL NUMERICAL_FLUX(IT)
                  FB_HAT = F_HAT
                  GB_HAT = G_HAT

#ifdef TRACE
                  IB_HAT = I_HAT
#endif

#ifdef CHEM
                  IB_HAT = I_HAT
                  JB_HAT = J_HAT
#endif

#ifdef DYNP
                  KB_HAT = K_HAT
#endif

                  IF (WDFLG(ELB).EQ.0) THEN
                     NLEQG = NLEQG_TMP
                     G = G_TMP
                  ENDIF
                  HB_HAT = H_HAT
                  GOTO 200
               ENDIF
            ENDIF
            
 200        CONTINUE
c     
            DO K = 1,DOFS(el)

               WEGPB = 2.0*M_INV(K,pa)/AREAS(ELB)*XLEN(GEDB)
     &              *PHI_EDGE(K,GPB,LEDB,pa)*WEGP(GPB,pa)
               WEGPF = 2.0*M_INV(K,pa)/AREAS(ELF)*XLEN(GEDF)
     &              *PHI_EDGE(K,GPF,LEDF,pa)*WEGP(GPF,pa)

               RHS_ZE(K,ELB,IRK) = RHS_ZE(K,ELB,IRK) - WEGPB*FB_HAT
               RHS_QX(K,ELB,IRK) = RHS_QX(K,ELB,IRK) - WEGPB*GB_HAT
               RHS_QY(K,ELB,IRK) = RHS_QY(K,ELB,IRK) - WEGPB*HB_HAT

               RHS_ZE(K,ELF,IRK) = RHS_ZE(K,ELF,IRK) - WEGPF*FF_HAT
               RHS_QX(K,ELF,IRK) = RHS_QX(K,ELF,IRK) - WEGPF*GF_HAT
               RHS_QY(K,ELF,IRK) = RHS_QY(K,ELF,IRK) - WEGPF*HF_HAT

#ifdef TRACE
               RHS_iota(K,ELB,IRK) = RHS_iota(K,ELB,IRK) - WEGPB*IB_HAT
               RHS_iota(K,ELF,IRK) = RHS_iota(K,ELF,IRK) - WEGPB*IF_HAT
#endif

#ifdef CHEM
               RHS_iota(K,ELB,IRK) = RHS_iota(K,ELB,IRK) - WEGPB*IB_HAT
               RHS_iota(K,ELF,IRK) = RHS_iota(K,ELF,IRK) - WEGPB*IF_HAT
               
               RHS_iota2(K,ELB,IRK) = RHS_iota2(K,ELB,IRK) - WEGPB*JB_HAT
               RHS_iota2(K,ELF,IRK) = RHS_iota2(K,ELF,IRK) - WEGPB*JF_HAT
#endif

#ifdef DYNP
               RHS_dynP(K,ELB,IRK) = RHS_dynP(K,ELB,IRK) - WEGPB*KB_HAT
               RHS_dynP(K,ELF,IRK) = RHS_dynP(K,ELF,IRK) - WEGPB*KF_HAT
#endif

            ENDDO
         ENDDO
 1000 CONTINUE
      RETURN
      END SUBROUTINE IBARRIER_EDGE_HYDRO
