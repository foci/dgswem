!***********************************************************************
!     
!     SUBROUTINE IBARRIER_EDGE_HYDRO()
!     
!     This subroutine does the following:
!     
!     1.  Calculates the values of the necessary variables at the global_here%edge
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

      SUBROUTINE IBARRIER_EDGE_HYDRO(s,dg_here,global_here,IT)
      
!.....Use appropriate modules

      USE GLOBAL
      USE DG
      USE SIZES

      IMPLICIT NONE

      type (sizes_type) :: s
      type (dg_type) :: dg_here
      type (global_type) :: global_here

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
!.....segment consists of two internal barrier edges -- a "front" global_here%edge &
!.....a "back" global_here%edge.
      dg_here%test_el = 0
      DO 1000 L = 1,dg_here%NIBSEG
         
!.......Obtain the global and local edges of the back and front sides

         GEDB = dg_here%NIBSEGN(1,L)
         GEDF = dg_here%NIBSEGN(2,L)

         if (gedf.eq.0) go to 1000
         
         LEDB = dg_here%NEDSD(1,GEDB)
         LEDF = dg_here%NEDSD(1,GEDF)
         
!.......Obtain the elements of the back and front sides

         ELB = dg_here%NEDEL(1,GEDB)
         ELF = dg_here%NEDEL(1,GEDF)

         if (dg_here%DOFS(ELB).LT.dg_here%DOFS(ELF)) then
            dg_here%EL = ELF
         endif

         dg_here%pa = global_here%PDG_EL(dg_here%EL)

#ifdef P0         
         if (dg_here%pa.eq.0) then
            dg_here%pa = 1
         endif
#endif
         
         NB1 = dg_here%NEDNO(1,GEDB)
         NB2 = dg_here%NEDNO(2,GEDB)
         
         NF1 = dg_here%NEDNO(1,GEDF)
         NF2 = dg_here%NEDNO(2,GEDF)
         IF((dg_here%WDFLG(ELB).EQ.0).AND.(dg_here%WDFLG(ELF).EQ.0)) GOTO 1000

         dg_here%test_el = dg_here%test_el+1

!.......Retrieve the components of the normal vector to the global_here%edge

         NXB = dg_here%COSNX(GEDB)
         NYB = dg_here%SINNX(GEDB)
         
         NXF = dg_here%COSNX(GEDF)
         NYF = dg_here%SINNX(GEDF)
         
!.....Set the components for the tangential vector to the global_here%edge

         TXB = -NYB
         TYB =  NXB

         TXF = -NYF
         TYF =  NXF
         
!.......Compute the variables at the quadrature points

         DO I = 1,dg_here%NEGP(dg_here%pa)
            
            GPB = I
            GPF = dg_here%NEGP(dg_here%pa) - I + 1
            
!.........Obtain the height of the barrier at the quadrature point

            ZEB = dg_here%ZE(1,ELB,dg_here%IRK)
            QXB = dg_here%QX(1,ELB,dg_here%IRK)
            QYB = dg_here%QY(1,ELB,dg_here%IRK)
            HBB = dg_here%BATHED(GPB,LEDB,ELB,dg_here%pa)
            
            ZEF = dg_here%ZE(1,ELF,dg_here%IRK)
            QXF = dg_here%QX(1,ELF,dg_here%IRK)
            QYF = dg_here%QY(1,ELF,dg_here%IRK)
            HBF = dg_here%BATHED(GPF,LEDF,ELF,dg_here%pa)

#ifdef TRACE
            iotaB = dg_here%iota(1,ELB,dg_here%IRK)
            iotaF = dg_here%iota(1,ELF,dg_here%IRK)
#endif

#ifdef CHEM
            iotaB = dg_here%iota(1,ELB,dg_here%IRK)
            iota2B = dg_here%iota2(1,ELB,dg_here%IRK)
            iotaF = dg_here%iota(1,ELF,dg_here%IRK)
            iota2F = dg_here%iota2(1,ELF,dg_here%IRK)
#endif

#ifdef DYNP
            dynPB = dg_here%dynP(1,ELB,dg_here%IRK)
            dynPF = dg_here%dynP(1,ELF,dg_here%IRK)
#endif

!This does not affect much here, dg_here%hb and dg_here%bed do not change at weirs
#ifdef SED_LAY 
            dg_here%HB(1,ELB,dg_here%irk) = 0.D0
            do ll=1,s%layers
               dg_here%HB(1,ELB,dg_here%irk) = dg_here%HB(1,ELB,dg_here%irk) + dg_here%bed(1,ELB,dg_here%irk,ll)
            enddo
            HBB = dg_here%HB(1,ELB,dg_here%irk)
            dg_here%HB(1,ELF,dg_here%irk) = 0.D0
            do ll=1,s%layers 
               dg_here%HB(1,ELF,dg_here%irk) = dg_here%HB(1,ELF,dg_here%irk) + dg_here%bed(1,ELF,dg_here%irk,ll)
            enddo
            HBF = dg_here%HB(1,ELF,dg_here%irk)
#endif



            DO K = 2,dg_here%dofs(dg_here%EL)
               ZEB = ZEB + dg_here%ZE(K,ELB,dg_here%IRK)*dg_here%PHI_EDGE(K,GPB,LEDB,dg_here%pa)
               QXB = QXB + dg_here%QX(K,ELB,dg_here%IRK)*dg_here%PHI_EDGE(K,GPB,LEDB,dg_here%pa)
               QYB = QYB + dg_here%QY(K,ELB,dg_here%IRK)*dg_here%PHI_EDGE(K,GPB,LEDB,dg_here%pa)
                                !HBB = HBB + dg_here%HB(K,ELB,1  )*dg_here%PHI_EDGE(K,GPB,LEDB,dg_here%pa)

               ZEF = ZEF + dg_here%ZE(K,ELF,dg_here%IRK)*dg_here%PHI_EDGE(K,GPF,LEDF,dg_here%pa)
               QXF = QXF + dg_here%QX(K,ELF,dg_here%IRK)*dg_here%PHI_EDGE(K,GPF,LEDF,dg_here%pa)
               QYF = QYF + dg_here%QY(K,ELF,dg_here%IRK)*dg_here%PHI_EDGE(K,GPF,LEDF,dg_here%pa)
                                !HBF = HBF + dg_here%HB(K,ELF,1  )*dg_here%PHI_EDGE(K,GPF,LEDF,dg_here%pa)
#ifdef TRACE
               iotaB = iotaB + dg_here%iota(K,ELB,dg_here%IRK)*dg_here%PHI_EDGE(K,GPB,LEDB,dg_here%pa)
               iotaF = iotaF + dg_here%iota(K,ELF,dg_here%IRK)*dg_here%PHI_EDGE(K,GPF,LEDF,dg_here%pa)
#endif

#ifdef CHEM
               iotaB = iotaB + dg_here%iota(K,ELB,dg_here%IRK)*dg_here%PHI_EDGE(K,GPB,LEDB,dg_here%pa)
               iotaF = iotaF + dg_here%iota(K,ELF,dg_here%IRK)*dg_here%PHI_EDGE(K,GPF,LEDF,dg_here%pa)
               iota2B = iota2B + dg_here%iota2(K,ELB,dg_here%IRK)*dg_here%PHI_EDGE(K,GPB,LEDB,dg_here%pa)
               iota2B = iota2F + dg_here%iota2(K,ELF,dg_here%IRK)*dg_here%PHI_EDGE(K,GPF,LEDF,dg_here%pa)
#endif

#ifdef DYNP
               dynPB = dynPB + dg_here%dynP(K,ELB,dg_here%IRK)*dg_here%PHI_EDGE(K,GPB,LEDB,dg_here%pa)
               dynPF = dynPF + dg_here%dynP(K,ELF,dg_here%IRK)*dg_here%PHI_EDGE(K,GPF,LEDF,dg_here%pa)
#endif

            ENDDO


            SFACB = dg_here%SFACED(I,LEDB,ELB,dg_here%pa)
            SFACF = dg_here%SFACED(I,LEDF,ELF,dg_here%pa)
            
            ABOVEB = 0
            ABOVEF = 0 
            IF (dg_here%WDFLG(ELB).EQ.1) ABOVEB = ZEB - dg_here%IBHT(L)
            IF (dg_here%WDFLG(ELF).EQ.1) ABOVEF = ZEF - dg_here%IBHT(L)
            
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

!.........Case 3: Overtopping of barrier with water higher on front global_here%side
!     Flow from front to back global_here%side
!     -------------------------------------------------------------

            IF ((ABOVEF.GT.ABOVEB).AND.(ABOVEF.GT.BARMIN)) THEN
               
               WEIR_FLOW = 1
               
!...........Case 3a) Subcritical flow
!     -------------------------

               IF (ABOVEB.GT.SUBSUPF) THEN
                  QF_N_INT = dg_here%RAMPDG*dg_here%IBCFSB(L)*ABOVEB&
                      *SQRT((2.D0*global_here%G*(ABOVEF-ABOVEB)))
                  QF_T_INT = 0.D0
                  
!...........Case 3b) Supercritical flow
!     ---------------------------

               ELSE
                  QF_N_INT = dg_here%RAMPDG*dg_here%IBCFSP(L)*SUBSUPF*SQRT(SUBSUPF*global_here%G)
                  QF_T_INT = 0.D0
               ENDIF
               GOTO 100

            ENDIF
            
!.........Case 4: Overtopping of barrier with water higher on back global_here%side
!     Flow from back to front global_here%side
!     --------------------------------------------------------------

            IF ((ABOVEB.GT.ABOVEF).AND.(ABOVEB.GT.BARMIN)) THEN
               
               WEIR_FLOW = -1
               
!...........Case 4a) Subcritical flow
!     -------------------------

               IF (ABOVEF.GT.SUBSUPB) THEN
                  QB_N_INT = dg_here%RAMPDG*dg_here%IBCFSB(L)*ABOVEF&
                      *SQRT((2.D0*global_here%G*(ABOVEB-ABOVEF)))
                  QB_T_INT = 0.D0
                  
!...........Case 4b) Supercritical flow
!     ---------------------------
                  
               ELSE
                  QB_N_INT = dg_here%RAMPDG*dg_here%IBCFSP(L)*SUBSUPB*SQRT(SUBSUPB*global_here%G)
                  QB_T_INT = 0.D0
               ENDIF
               GOTO 100
            ENDIF
            
 100        CONTINUE

            IF (WEIR_FLOW.LE.0) THEN

!...........Compute the numerical flux for the back global_here%side global_here%edge

               dg_here%ZE_IN = ZEB
               dg_here%QX_IN = QXB
               dg_here%QY_IN = QYB
               dg_here%HB_IN = HBB

               dg_here%SFAC_IN = SFACB

#ifdef TRACE
               dg_here%iota_IN = iotaB
#endif

#ifdef CHEM
               dg_here%iota_IN = iotaB
               dg_here%iota2_IN = iota2B
#endif

#ifdef DYNP
               dynP_IN = dynPB
#endif
               
               dg_here%ZE_EX = ZEB
               dg_here%HB_EX = HBB

               dg_here%SFAC_EX = SFACB

#ifdef TRACE
               dg_here%iota_EX = iotaB
#endif

#ifdef CHEM
               dg_here%iota_EX = iotaB
               dg_here%iota2_EX = iota2B
#endif

#ifdef DYNP
               dynP_EX = dynPB
#endif
               
               dg_here%NX = NXB
               dg_here%NY = NYB

!...........Reflect the velocity in the normal direction

               global_here%Q_N_EXT = QB_N_INT
               global_here%Q_T_EXT = QB_T_INT

!...........Compute the global_here%x and global_here%y components of the external state flow

               dg_here%QX_EX = ( TYB*global_here%Q_N_EXT - NYB*global_here%Q_T_EXT)/(NXB*TYB - NYB*TXB)
               dg_here%QY_EX = (-TXB*global_here%Q_N_EXT + NXB*global_here%Q_T_EXT)/(NXB*TYB - NYB*TXB)

!...........Compute the numerical flux
               
               CALL NUMERICAL_FLUX(s,dg_here,global_here,IT)
               FB_HAT = global_here%F_HAT
               GB_HAT = global_here%G_HAT
               HB_HAT = global_here%H_HAT

#ifdef TRACE
               IB_HAT = global_here%I_HAT
#endif

#ifdef CHEM
               IB_HAT = global_here%I_HAT
               JB_HAT = global_here%J_HAT
#endif

#ifdef DYNP
               KB_HAT = global_here%K_HAT
#endif

               IF (WEIR_FLOW.LT.0) THEN
                  dg_here%ZE_IN = ZEF
                  dg_here%QX_IN = QXF
                  dg_here%QY_IN = QYF
                  dg_here%HB_IN = HBF
                  dg_here%ZE_EX = ZEF
                  dg_here%HB_EX = HBF

                  dg_here%SFAC_IN = SFACF
                  dg_here%SFAC_EX = SFACF


                  dg_here%NX = NXF
                  dg_here%NY = NYF
                  
#ifdef TRACE
                  dg_here%iota_IN = iotaF
                  dg_here%iota_EX = iotaF
#endif

#ifdef CHEM
                  dg_here%iota_IN = iotaF
                  dg_here%iota2_IN = iota2F
                  dg_here%iota_EX = iotaF
                  dg_here%iota2_EX = iota2F
#endif

#ifdef DYNP
                  dynP_IN = iotaF
                  dynP_EX = iotaF
#endif

                  IF (dg_here%WDFLG(ELF).EQ.0) THEN
                     NLEQG_TMP = global_here%NLEQG
                     global_here%NLEQG = 0.D0
                     G_TMP = global_here%G
                     global_here%G = 0.D0
                  ENDIF
                  call numerical_flux(s,dg_here,global_here,IT)
                  FF_HAT = global_here%F_HAT
                  GF_HAT = global_here%G_HAT
                  HF_HAT = global_here%H_HAT

#ifdef TRACE
                  IF_HAT = global_here%I_HAT
#endif

#ifdef CHEM
                  IF_HAT = global_here%I_HAT
                  JF_HAT = global_here%J_HAT
#endif

#ifdef DYNP
                  KF_HAT = global_here%K_HAT
#endif

                  IF (dg_here%WDFLG(ELF).EQ.0) THEN
                     global_here%NLEQG = NLEQG_TMP
                     global_here%G = G_TMP
                  ENDIF
                  GOTO 200
               ENDIF
            ENDIF

            IF (WEIR_FLOW.GE.0) THEN
               
!...........Compute the numerical flux for the front global_here%side global_here%edge

               dg_here%ZE_IN = ZEF
               dg_here%QX_IN = QXF
               dg_here%QY_IN = QYF
               dg_here%HB_IN = HBF

               dg_here%ZE_EX = ZEF
               dg_here%HB_EX = HBF

               dg_here%SFAC_IN = SFACF
               dg_here%SFAC_EX = SFACF

#ifdef TRACE
               dg_here%iota_IN = iotaF
               dg_here%iota_EX = iotaF
#endif

#ifdef CHEM
               dg_here%iota_IN = iotaF
               dg_here%iota_EX = iotaF
               dg_here%iota2_IN = iota2F
               dg_here%iota2_EX = iota2F
#endif

#ifdef DYNP
               dynP_IN = dynPF
               dynP_EX = dynPF
#endif

               dg_here%NX = NXF
               dg_here%NY = NYF
               
!...........Reflect the velocity in the normal direction

               global_here%Q_N_EXT = QF_N_INT
               global_here%Q_T_EXT = QF_T_INT

!...........Compute the global_here%x and global_here%y components of the external state flow

               dg_here%QX_EX = ( TYF*global_here%Q_N_EXT - NYF*global_here%Q_T_EXT)/(NXF*TYF - NYF*TXF)
               dg_here%QY_EX = (-TXF*global_here%Q_N_EXT + NXF*global_here%Q_T_EXT)/(NXF*TYF - NYF*TXF)

               call numerical_flux(s,dg_here,global_here,IT)
               FF_HAT = global_here%F_HAT
               GF_HAT = global_here%G_HAT
               HF_HAT = global_here%H_HAT

#ifdef TRACE
               IF_HAT = global_here%I_HAT
#endif

#ifdef CHEM
               IF_HAT = global_here%I_HAT
               JF_HAT = global_here%J_HAT
#endif

#ifdef DYNP
               KF_HAT = global_here%K_HAT
#endif
               
               IF (WEIR_FLOW.GT.0) THEN
                  dg_here%ZE_IN = ZEB
                  dg_here%QX_IN = QXB
                  dg_here%QY_IN = QYB
                  dg_here%HB_IN = HBB
                  dg_here%ZE_EX = ZEB
                  dg_here%HB_EX = HBB

                  dg_here%SFAC_IN = SFACB
                  dg_here%SFAC_EX = SFACB

#ifdef TRACE
                  dg_here%iota_IN = iotaB
                  dg_here%iota_EX = iotaB
#endif

#ifdef CHEM
                  dg_here%iota_IN = iotaB
                  dg_here%iota2_IN = iota2B
                  dg_here%iota_EX = iotaB
                  dg_here%iota2_EX = iota2B
#endif

#ifdef DYNP
                  dynP_IN = dynPB
                  dynP_EX = dynPB
#endif

                  dg_here%NX = NXB
                  dg_here%NY = NYB
                  IF (dg_here%WDFLG(ELB).EQ.0) THEN
                     NLEQG_TMP = global_here%NLEQG
                     global_here%NLEQG = 0.D0
                     G_TMP = global_here%G
                     global_here%G = 0.D0
                  ENDIF
                  call numerical_flux(s,dg_here,global_here,IT)
                  FB_HAT = global_here%F_HAT
                  GB_HAT = global_here%G_HAT

#ifdef TRACE
                  IB_HAT = global_here%I_HAT
#endif

#ifdef CHEM
                  IB_HAT = global_here%I_HAT
                  JB_HAT = global_here%J_HAT
#endif

#ifdef DYNP
                  KB_HAT = global_here%K_HAT
#endif

                  IF (dg_here%WDFLG(ELB).EQ.0) THEN
                     global_here%NLEQG = NLEQG_TMP
                     global_here%G = G_TMP
                  ENDIF
                  HB_HAT = global_here%H_HAT
                  GOTO 200
               ENDIF
            ENDIF
            
 200        CONTINUE
!     
            DO K = 1,dg_here%DOFS(dg_here%el)

               WEGPB = 2.0*dg_here%M_INV(K,dg_here%pa)/global_here%AREAS(ELB)*dg_here%XLEN(GEDB)&
                   *dg_here%PHI_EDGE(K,GPB,LEDB,dg_here%pa)*dg_here%WEGP(GPB,dg_here%pa)
               WEGPF = 2.0*dg_here%M_INV(K,dg_here%pa)/global_here%AREAS(ELF)*dg_here%XLEN(GEDF)&
                   *dg_here%PHI_EDGE(K,GPF,LEDF,dg_here%pa)*dg_here%WEGP(GPF,dg_here%pa)

               dg_here%RHS_ZE(K,ELB,dg_here%IRK) = dg_here%RHS_ZE(K,ELB,dg_here%IRK) - WEGPB*FB_HAT
               dg_here%RHS_QX(K,ELB,dg_here%IRK) = dg_here%RHS_QX(K,ELB,dg_here%IRK) - WEGPB*GB_HAT
               dg_here%RHS_QY(K,ELB,dg_here%IRK) = dg_here%RHS_QY(K,ELB,dg_here%IRK) - WEGPB*HB_HAT

               dg_here%RHS_ZE(K,ELF,dg_here%IRK) = dg_here%RHS_ZE(K,ELF,dg_here%IRK) - WEGPF*FF_HAT
               dg_here%RHS_QX(K,ELF,dg_here%IRK) = dg_here%RHS_QX(K,ELF,dg_here%IRK) - WEGPF*GF_HAT
               dg_here%RHS_QY(K,ELF,dg_here%IRK) = dg_here%RHS_QY(K,ELF,dg_here%IRK) - WEGPF*HF_HAT

#ifdef TRACE
               dg_here%RHS_iota(K,ELB,dg_here%IRK) = dg_here%RHS_iota(K,ELB,dg_here%IRK) - WEGPB*IB_HAT
               dg_here%RHS_iota(K,ELF,dg_here%IRK) = dg_here%RHS_iota(K,ELF,dg_here%IRK) - WEGPB*IF_HAT
#endif

#ifdef CHEM
               dg_here%RHS_iota(K,ELB,dg_here%IRK) = dg_here%RHS_iota(K,ELB,dg_here%IRK) - WEGPB*IB_HAT
               dg_here%RHS_iota(K,ELF,dg_here%IRK) = dg_here%RHS_iota(K,ELF,dg_here%IRK) - WEGPB*IF_HAT
               
               dg_here%RHS_iota2(K,ELB,dg_here%IRK) = dg_here%RHS_iota2(K,ELB,dg_here%IRK) - WEGPB*JB_HAT
               dg_here%RHS_iota2(K,ELF,dg_here%IRK) = dg_here%RHS_iota2(K,ELF,dg_here%IRK) - WEGPB*JF_HAT
#endif

#ifdef DYNP
               dg_here%RHS_dynP(K,ELB,dg_here%IRK) = dg_here%RHS_dynP(K,ELB,dg_here%IRK) - WEGPB*KB_HAT
               dg_here%RHS_dynP(K,ELF,dg_here%IRK) = dg_here%RHS_dynP(K,ELF,dg_here%IRK) - WEGPB*KF_HAT
#endif

            ENDDO
         ENDDO
 1000 CONTINUE
      RETURN
      END SUBROUTINE IBARRIER_EDGE_HYDRO
