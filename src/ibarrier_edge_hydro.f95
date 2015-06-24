!***********************************************************************
!     
!     SUBROUTINE IBARRIER_EDGE_HYDRO()
!     
!     This subroutine does the following:
!     
!     1.  Calculates the values of the necessary variables at the edge
!     gauss points for INTERNAL BARRIER edges
!     2.  Calls the appropriate subroutine to compute the flux at
!     these points.
!     3.  Computes the boundary integrals.
!     
!     Written by Ethan Kubatko (08-11-2008)
!     
!     01-10-2011 - cem - adapted for p_enrichment and multicomponent
!     01-20-2012 - cem - multiphase added, we do not (natively)
!                        transport sediment across weirs
!     
!***********************************************************************

      SUBROUTINE IBARRIER_EDGE_HYDRO(s,IT)
      
!.....Use appropriate modules

      USE GLOBAL
      USE DG
      USE SIZES

      IMPLICIT NONE

      type (sizes_type) :: s

!.....Declare local variables

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

!.....Loop over the internal barrier segments (note: an internal barrier
!.....segment consists of two internal barrier edges -- a "front" edge &
!.....a "back" edge.
      dg%test_el = 0
      DO 1000 L = 1,dg%NIBSEG
         
!.......Obtain the global and local edges of the back and front sides

         GEDB = dg%NIBSEGN(1,L)
         GEDF = dg%NIBSEGN(2,L)

         if (gedf.eq.0) go to 1000
         
         LEDB = dg%NEDSD(1,GEDB)
         LEDF = dg%NEDSD(1,GEDF)
         
!.......Obtain the elements of the back and front sides

         ELB = dg%NEDEL(1,GEDB)
         ELF = dg%NEDEL(1,GEDF)

         if (dg%DOFS(ELB).LT.dg%DOFS(ELF)) then
            dg%EL = ELF
         endif

         dg%pa = PDG_EL(dg%EL)

#ifdef P0         
         if (dg%pa.eq.0) then
            dg%pa = 1
         endif
#endif
         
         NB1 = dg%NEDNO(1,GEDB)
         NB2 = dg%NEDNO(2,GEDB)
         
         NF1 = dg%NEDNO(1,GEDF)
         NF2 = dg%NEDNO(2,GEDF)
         IF((dg%WDFLG(ELB).EQ.0).AND.(dg%WDFLG(ELF).EQ.0)) GOTO 1000

         dg%test_el = dg%test_el+1

!.......Retrieve the components of the normal vector to the edge

         NXB = dg%COSNX(GEDB)
         NYB = dg%SINNX(GEDB)
         
         NXF = dg%COSNX(GEDF)
         NYF = dg%SINNX(GEDF)
         
!.....Set the components for the tangential vector to the edge

         TXB = -NYB
         TYB =  NXB

         TXF = -NYF
         TYF =  NXF
         
!.......Compute the variables at the quadrature points

         DO I = 1,dg%NEGP(dg%pa)
            
            GPB = I
            GPF = dg%NEGP(dg%pa) - I + 1
            
!.........Obtain the height of the barrier at the quadrature point

            ZEB = dg%ZE(1,ELB,dg%IRK)
            QXB = dg%QX(1,ELB,dg%IRK)
            QYB = dg%QY(1,ELB,dg%IRK)
            HBB = dg%BATHED(GPB,LEDB,ELB,dg%pa)
            
            ZEF = dg%ZE(1,ELF,dg%IRK)
            QXF = dg%QX(1,ELF,dg%IRK)
            QYF = dg%QY(1,ELF,dg%IRK)
            HBF = dg%BATHED(GPF,LEDF,ELF,dg%pa)

#ifdef TRACE
            iotaB = dg%iota(1,ELB,dg%IRK)
            iotaF = dg%iota(1,ELF,dg%IRK)
#endif

#ifdef CHEM
            iotaB = dg%iota(1,ELB,dg%IRK)
            iota2B = dg%iota2(1,ELB,dg%IRK)
            iotaF = dg%iota(1,ELF,dg%IRK)
            iota2F = dg%iota2(1,ELF,dg%IRK)
#endif

#ifdef DYNP
            dynPB = dg%dynP(1,ELB,dg%IRK)
            dynPF = dg%dynP(1,ELF,dg%IRK)
#endif

!This does not affect much here, dg%hb and dg%bed do not change at weirs
#ifdef SED_LAY 
            dg%HB(1,ELB,dg%irk) = 0.D0
            do ll=1,s%layers
               dg%HB(1,ELB,dg%irk) = dg%HB(1,ELB,dg%irk) + dg%bed(1,ELB,dg%irk,ll)
            enddo
            HBB = dg%HB(1,ELB,dg%irk)
            dg%HB(1,ELF,dg%irk) = 0.D0
            do ll=1,s%layers 
               dg%HB(1,ELF,dg%irk) = dg%HB(1,ELF,dg%irk) + dg%bed(1,ELF,dg%irk,ll)
            enddo
            HBF = dg%HB(1,ELF,dg%irk)
#endif



            DO K = 2,dg%dofs(dg%EL)
               ZEB = ZEB + dg%ZE(K,ELB,dg%IRK)*dg%PHI_EDGE(K,GPB,LEDB,dg%pa)
               QXB = QXB + dg%QX(K,ELB,dg%IRK)*dg%PHI_EDGE(K,GPB,LEDB,dg%pa)
               QYB = QYB + dg%QY(K,ELB,dg%IRK)*dg%PHI_EDGE(K,GPB,LEDB,dg%pa)
                                !HBB = HBB + dg%HB(K,ELB,1  )*dg%PHI_EDGE(K,GPB,LEDB,dg%pa)

               ZEF = ZEF + dg%ZE(K,ELF,dg%IRK)*dg%PHI_EDGE(K,GPF,LEDF,dg%pa)
               QXF = QXF + dg%QX(K,ELF,dg%IRK)*dg%PHI_EDGE(K,GPF,LEDF,dg%pa)
               QYF = QYF + dg%QY(K,ELF,dg%IRK)*dg%PHI_EDGE(K,GPF,LEDF,dg%pa)
                                !HBF = HBF + dg%HB(K,ELF,1  )*dg%PHI_EDGE(K,GPF,LEDF,dg%pa)
#ifdef TRACE
               iotaB = iotaB + dg%iota(K,ELB,dg%IRK)*dg%PHI_EDGE(K,GPB,LEDB,dg%pa)
               iotaF = iotaF + dg%iota(K,ELF,dg%IRK)*dg%PHI_EDGE(K,GPF,LEDF,dg%pa)
#endif

#ifdef CHEM
               iotaB = iotaB + dg%iota(K,ELB,dg%IRK)*dg%PHI_EDGE(K,GPB,LEDB,dg%pa)
               iotaF = iotaF + dg%iota(K,ELF,dg%IRK)*dg%PHI_EDGE(K,GPF,LEDF,dg%pa)
               iota2B = iota2B + dg%iota2(K,ELB,dg%IRK)*dg%PHI_EDGE(K,GPB,LEDB,dg%pa)
               iota2B = iota2F + dg%iota2(K,ELF,dg%IRK)*dg%PHI_EDGE(K,GPF,LEDF,dg%pa)
#endif

#ifdef DYNP
               dynPB = dynPB + dg%dynP(K,ELB,dg%IRK)*dg%PHI_EDGE(K,GPB,LEDB,dg%pa)
               dynPF = dynPF + dg%dynP(K,ELF,dg%IRK)*dg%PHI_EDGE(K,GPF,LEDF,dg%pa)
#endif

            ENDDO


            SFACB = dg%SFACED(I,LEDB,ELB,dg%pa)
            SFACF = dg%SFACED(I,LEDF,ELF,dg%pa)
            
            ABOVEB = 0
            ABOVEF = 0 
            IF (dg%WDFLG(ELB).EQ.1) ABOVEB = ZEB - dg%IBHT(L)
            IF (dg%WDFLG(ELF).EQ.1) ABOVEF = ZEF - dg%IBHT(L)
            
!.........Case 1:  Water is below barrier on both sides
!     ---------------------------------------------
            
            IF ((ABOVEF.LE.BARMIN).AND.(ABOVEB.LE.BARMIN)) THEN
               QB_N_INT = -(QXB*NXB + QYB*NYB)
               QF_N_INT = -(QXF*NXF + QYF*NYF)
               QB_T_INT = QXB*TXB + QYB*TYB
               QF_T_INT = QXF*TXF + QYF*TYF
               WEIR_FLOW = 0
               GOTO 100
               
!.........Case 2:  Water is above on both sides and equal (within tol)
!     ------------------------------------------------------------
               
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

!.........Case 3: Overtopping of barrier with water higher on front side
!     Flow from front to back side
!     -------------------------------------------------------------

            IF ((ABOVEF.GT.ABOVEB).AND.(ABOVEF.GT.BARMIN)) THEN
               
               WEIR_FLOW = 1
               
!...........Case 3a) Subcritical flow
!     -------------------------

               IF (ABOVEB.GT.SUBSUPF) THEN
                  QF_N_INT = dg%RAMPDG*dg%IBCFSB(L)*ABOVEB&
                      *SQRT((2.D0*G*(ABOVEF-ABOVEB)))
                  QF_T_INT = 0.D0
                  
!...........Case 3b) Supercritical flow
!     ---------------------------

               ELSE
                  QF_N_INT = dg%RAMPDG*dg%IBCFSP(L)*SUBSUPF*SQRT(SUBSUPF*G)
                  QF_T_INT = 0.D0
               ENDIF
               GOTO 100

            ENDIF
            
!.........Case 4: Overtopping of barrier with water higher on back side
!     Flow from back to front side
!     --------------------------------------------------------------

            IF ((ABOVEB.GT.ABOVEF).AND.(ABOVEB.GT.BARMIN)) THEN
               
               WEIR_FLOW = -1
               
!...........Case 4a) Subcritical flow
!     -------------------------

               IF (ABOVEF.GT.SUBSUPB) THEN
                  QB_N_INT = dg%RAMPDG*dg%IBCFSB(L)*ABOVEF&
                      *SQRT((2.D0*G*(ABOVEB-ABOVEF)))
                  QB_T_INT = 0.D0
                  
!...........Case 4b) Supercritical flow
!     ---------------------------
                  
               ELSE
                  QB_N_INT = dg%RAMPDG*dg%IBCFSP(L)*SUBSUPB*SQRT(SUBSUPB*G)
                  QB_T_INT = 0.D0
               ENDIF
               GOTO 100
            ENDIF
            
 100        CONTINUE

            IF (WEIR_FLOW.LE.0) THEN

!...........Compute the numerical flux for the back side edge

               dg%ZE_IN = ZEB
               dg%QX_IN = QXB
               dg%QY_IN = QYB
               dg%HB_IN = HBB

               dg%SFAC_IN = SFACB

#ifdef TRACE
               dg%iota_IN = iotaB
#endif

#ifdef CHEM
               dg%iota_IN = iotaB
               dg%iota2_IN = iota2B
#endif

#ifdef DYNP
               dynP_IN = dynPB
#endif
               
               dg%ZE_EX = ZEB
               dg%HB_EX = HBB

               dg%SFAC_EX = SFACB

#ifdef TRACE
               dg%iota_EX = iotaB
#endif

#ifdef CHEM
               dg%iota_EX = iotaB
               dg%iota2_EX = iota2B
#endif

#ifdef DYNP
               dynP_EX = dynPB
#endif
               
               dg%NX = NXB
               dg%NY = NYB

!...........Reflect the velocity in the normal direction

               Q_N_EXT = QB_N_INT
               Q_T_EXT = QB_T_INT

!...........Compute the x and y components of the external state flow

               dg%QX_EX = ( TYB*Q_N_EXT - NYB*Q_T_EXT)/(NXB*TYB - NYB*TXB)
               dg%QY_EX = (-TXB*Q_N_EXT + NXB*Q_T_EXT)/(NXB*TYB - NYB*TXB)

!...........Compute the numerical flux
               
               CALL NUMERICAL_FLUX(s,IT)
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
                  dg%ZE_IN = ZEF
                  dg%QX_IN = QXF
                  dg%QY_IN = QYF
                  dg%HB_IN = HBF
                  dg%ZE_EX = ZEF
                  dg%HB_EX = HBF

                  dg%SFAC_IN = SFACF
                  dg%SFAC_EX = SFACF


                  dg%NX = NXF
                  dg%NY = NYF
                  
#ifdef TRACE
                  dg%iota_IN = iotaF
                  dg%iota_EX = iotaF
#endif

#ifdef CHEM
                  dg%iota_IN = iotaF
                  dg%iota2_IN = iota2F
                  dg%iota_EX = iotaF
                  dg%iota2_EX = iota2F
#endif

#ifdef DYNP
                  dynP_IN = iotaF
                  dynP_EX = iotaF
#endif

                  IF (dg%WDFLG(ELF).EQ.0) THEN
                     NLEQG_TMP = NLEQG
                     NLEQG = 0.D0
                     G_TMP = G
                     G = 0.D0
                  ENDIF
                  CALL NUMERICAL_FLUX(s,IT)
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

                  IF (dg%WDFLG(ELF).EQ.0) THEN
                     NLEQG = NLEQG_TMP
                     G = G_TMP
                  ENDIF
                  GOTO 200
               ENDIF
            ENDIF

            IF (WEIR_FLOW.GE.0) THEN
               
!...........Compute the numerical flux for the front side edge

               dg%ZE_IN = ZEF
               dg%QX_IN = QXF
               dg%QY_IN = QYF
               dg%HB_IN = HBF

               dg%ZE_EX = ZEF
               dg%HB_EX = HBF

               dg%SFAC_IN = SFACF
               dg%SFAC_EX = SFACF

#ifdef TRACE
               dg%iota_IN = iotaF
               dg%iota_EX = iotaF
#endif

#ifdef CHEM
               dg%iota_IN = iotaF
               dg%iota_EX = iotaF
               dg%iota2_IN = iota2F
               dg%iota2_EX = iota2F
#endif

#ifdef DYNP
               dynP_IN = dynPF
               dynP_EX = dynPF
#endif

               dg%NX = NXF
               dg%NY = NYF
               
!...........Reflect the velocity in the normal direction

               Q_N_EXT = QF_N_INT
               Q_T_EXT = QF_T_INT

!...........Compute the x and y components of the external state flow

               dg%QX_EX = ( TYF*Q_N_EXT - NYF*Q_T_EXT)/(NXF*TYF - NYF*TXF)
               dg%QY_EX = (-TXF*Q_N_EXT + NXF*Q_T_EXT)/(NXF*TYF - NYF*TXF)

               CALL NUMERICAL_FLUX(s,IT)
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
                  dg%ZE_IN = ZEB
                  dg%QX_IN = QXB
                  dg%QY_IN = QYB
                  dg%HB_IN = HBB
                  dg%ZE_EX = ZEB
                  dg%HB_EX = HBB

                  dg%SFAC_IN = SFACB
                  dg%SFAC_EX = SFACB

#ifdef TRACE
                  dg%iota_IN = iotaB
                  dg%iota_EX = iotaB
#endif

#ifdef CHEM
                  dg%iota_IN = iotaB
                  dg%iota2_IN = iota2B
                  dg%iota_EX = iotaB
                  dg%iota2_EX = iota2B
#endif

#ifdef DYNP
                  dynP_IN = dynPB
                  dynP_EX = dynPB
#endif

                  dg%NX = NXB
                  dg%NY = NYB
                  IF (dg%WDFLG(ELB).EQ.0) THEN
                     NLEQG_TMP = NLEQG
                     NLEQG = 0.D0
                     G_TMP = G
                     G = 0.D0
                  ENDIF
                  CALL NUMERICAL_FLUX(s,IT)
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

                  IF (dg%WDFLG(ELB).EQ.0) THEN
                     NLEQG = NLEQG_TMP
                     G = G_TMP
                  ENDIF
                  HB_HAT = H_HAT
                  GOTO 200
               ENDIF
            ENDIF
            
 200        CONTINUE
!     
            DO K = 1,dg%DOFS(dg%el)

               WEGPB = 2.0*dg%M_INV(K,dg%pa)/AREAS(ELB)*dg%XLEN(GEDB)&
                   *dg%PHI_EDGE(K,GPB,LEDB,dg%pa)*dg%WEGP(GPB,dg%pa)
               WEGPF = 2.0*dg%M_INV(K,dg%pa)/AREAS(ELF)*dg%XLEN(GEDF)&
                   *dg%PHI_EDGE(K,GPF,LEDF,dg%pa)*dg%WEGP(GPF,dg%pa)

               dg%RHS_ZE(K,ELB,dg%IRK) = dg%RHS_ZE(K,ELB,dg%IRK) - WEGPB*FB_HAT
               dg%RHS_QX(K,ELB,dg%IRK) = dg%RHS_QX(K,ELB,dg%IRK) - WEGPB*GB_HAT
               dg%RHS_QY(K,ELB,dg%IRK) = dg%RHS_QY(K,ELB,dg%IRK) - WEGPB*HB_HAT

               dg%RHS_ZE(K,ELF,dg%IRK) = dg%RHS_ZE(K,ELF,dg%IRK) - WEGPF*FF_HAT
               dg%RHS_QX(K,ELF,dg%IRK) = dg%RHS_QX(K,ELF,dg%IRK) - WEGPF*GF_HAT
               dg%RHS_QY(K,ELF,dg%IRK) = dg%RHS_QY(K,ELF,dg%IRK) - WEGPF*HF_HAT

#ifdef TRACE
               dg%RHS_iota(K,ELB,dg%IRK) = dg%RHS_iota(K,ELB,dg%IRK) - WEGPB*IB_HAT
               dg%RHS_iota(K,ELF,dg%IRK) = dg%RHS_iota(K,ELF,dg%IRK) - WEGPB*IF_HAT
#endif

#ifdef CHEM
               dg%RHS_iota(K,ELB,dg%IRK) = dg%RHS_iota(K,ELB,dg%IRK) - WEGPB*IB_HAT
               dg%RHS_iota(K,ELF,dg%IRK) = dg%RHS_iota(K,ELF,dg%IRK) - WEGPB*IF_HAT
               
               dg%RHS_iota2(K,ELB,dg%IRK) = dg%RHS_iota2(K,ELB,dg%IRK) - WEGPB*JB_HAT
               dg%RHS_iota2(K,ELF,dg%IRK) = dg%RHS_iota2(K,ELF,dg%IRK) - WEGPB*JF_HAT
#endif

#ifdef DYNP
               dg%RHS_dynP(K,ELB,dg%IRK) = dg%RHS_dynP(K,ELB,dg%IRK) - WEGPB*KB_HAT
               dg%RHS_dynP(K,ELF,dg%IRK) = dg%RHS_dynP(K,ELF,dg%IRK) - WEGPB*KF_HAT
#endif

            ENDDO
         ENDDO
 1000 CONTINUE
      RETURN
      END SUBROUTINE IBARRIER_EDGE_HYDRO
