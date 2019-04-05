C***********************************************************************
C     
C     SUBROUTINE LAND_EDGE_HYDRO( )
C     
C     This subroutine does the following:
C     
C     1.  Calculates the values of the necessary variables at the edge
C     gauss points for NO-NORMAL FLOW edges
C     2.  Calls the appropriate subroutine to compute the flux at
C     these points.
C     3.  Calls the appropriate subroutine to compute the boundary
C     integrals.
C     
C     Written by Ethan Kubatko (06-11-2004)
C     
C-----------------------------------------------------------------------
C     
C     01-02-2007, sb, Modified for LDG
C     08-xx-2005, sb, Modifications for wetting/drying
C     
C     01-10-2011 - cem - adapted for p_enrichment and multicomponent
C     06-02-2012 - cem - adapted for sediment
C     
C***********************************************************************

      SUBROUTINE LAND_EDGE_HYDRO(IT)

C.....Use appropriate modules

      USE GLOBAL
      USE DG
      Use SIZES, only: layers

      IMPLICIT NONE

C.....Declare local variables

      INTEGER L, LED, GED, i,k,kk,j,jj,II,ll,IT,mm
      Real(SZ) DEN2,U_AVG,V_AVG,VEL_NORMAL,q_RoeX, q_RoeY, q_Roe
      REAL(SZ) TX, TY, DEN,ell_1,ell_2,ell_3,HZ_X_IN,HZ_Y_IN
      REAL(SZ) LZ_XX_IN, LZ_XY_IN, LZ_YX_IN, LZ_YY_IN,W_IN
      Real(SZ) MZ_X_IN(layers),MZ_Y_IN(layers),TZ_X_IN,TZ_Y_IN
      Real(SZ), allocatable :: vtemp1(:), vtemp2(:)

      test_el = 0
      DO 1000 L = 1,NLEDS
         
C.....Retrieve the global and local edge number

         GED = NLEDN(L)
         LED = NEDSD(1,GED)

C.....Retrieve the elements which share the edge

         EL_IN = NEDEL(1,GED)
         J = EDN(EL_IN) - 2
         pa = PDG_EL(EL_IN)
         
         if (J .eq. 2) then
            allocate( vtemp1(DOFS(EL_IN)),vtemp2(DOFS(EL_IN)) )
         endif

#ifdef P0
         if (pa.eq.0) then
            pa = 1
         endif
#endif
         
C.....If the element is dry then skip the edge calculation

         IF (WDFLG(EL_IN).EQ.0) GOTO 1000

         test_el = test_el+1

C.....Retrieve the components of the normal vector to the edge
         
         NX = COSNX(GED)
         NY = SINNX(GED)
         
C.....Set the components for the tangential vector to the edge
         
         TX = -NY
         TY =  NX

C.....Compute ZE, QX, QY, and HB at each Gauss point

         DO I = 1,NEGP(pa)
            
            ZE_IN = ZE(1,EL_IN,IRK)
            QX_IN = QX(1,EL_IN,IRK)
            QY_IN = QY(1,EL_IN,IRK)
            HB_IN = BATHED(I,LED,EL_IN,pa)

            SFAC_IN = SFACED(I,LED,EL_IN,pa)

#ifdef TRACE
            iota_IN = iota(1,EL_IN,IRK)
#endif

#ifdef CHEM
            iota_IN = 0.D0 !iota(1,EL_IN,IRK)
            iota2_IN = 0.D0 !iota2(1,EL_IN,IRK)
#endif

#ifdef DYNP
            dynP_IN = dynP(1,EL_IN,IRK)
#endif

#ifdef SED_LAY
            HB(1,EL_IN,irk) = 0.D0
            do ll=1,layers
               HB(1,EL_IN,irk) = HB(1,EL_IN,irk) + bed(1,EL_IN,irk,ll)

               MZ_X_IN(ll) =  MZ(1,1,ll,EL_IN)
               MZ_Y_IN(ll) =  MZ(1,2,ll,EL_IN)
            enddo
            bed_IN(:) = bed(1,EL_IN,irk,:)
            HB_IN = HB(1,EL_IN,irk)

#endif

            !...do it and do it again ... 
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

            !Compute sediment diffusion contribution

C.....Compute the solution at the interior state
C.....(modified for wetting and drying)

            DO K = 2,DOFS(EL_IN)
               
               ZE_IN = ZE_IN + ZE(K,EL_IN,IRK)*PHI_EDGE(K,I,LED,pa,J)
               QX_IN = QX_IN + QX(K,EL_IN,IRK)*PHI_EDGE(K,I,LED,pa,J)
               QY_IN = QY_IN + QY(K,EL_IN,IRK)*PHI_EDGE(K,I,LED,pa,J)

                                !LDG terms
#ifdef WAVE_DIF
              HZ_X_IN = HZ_X_IN + HZ(K,1,1,EL_IN)*PHI_EDGE(K,I,LED,pa,J)
              HZ_Y_IN = HZ_Y_IN + HZ(K,2,2,EL_IN)*PHI_EDGE(K,I,LED,pa,J)
#endif
            LZ_XX_IN = LZ_XX_IN + LZ(K,1,1,EL_IN)*PHI_EDGE(K,I,LED,pa,J)
            LZ_XY_IN = LZ_XY_IN + LZ(K,1,2,EL_IN)*PHI_EDGE(K,I,LED,pa,J)
            LZ_YX_IN = LZ_YX_IN + LZ(K,2,1,EL_IN)*PHI_EDGE(K,I,LED,pa,J)
            LZ_YY_IN = LZ_YY_IN + LZ(K,2,2,EL_IN)*PHI_EDGE(K,I,LED,pa,J)

#ifdef TRACE

              TZ_X_IN = TZ_X_IN + TZ(K,1,1,EL_IN)*PHI_EDGE(K,I,LED,pa,J)
              TZ_Y_IN = TZ_Y_IN + TZ(K,2,2,EL_IN)*PHI_EDGE(K,I,LED,pa,J)

#endif

                                !LDG terms for sediment diffusion

#ifdef TRACE
            iota_IN = iota_IN + iota(K,EL_IN,IRK)*PHI_EDGE(K,I,LED,pa,J)
#endif

#ifdef CHEM
            iota_IN = iota_IN + iota(K,EL_IN,IRK)*PHI_EDGE(K,I,LED,pa,J)
         iota2_IN = iota2_IN + iota2(K,EL_IN,IRK)*PHI_EDGE(K,I,LED,pa,J)
#endif

#ifdef DYNP
            dynP_IN = dynP_IN + dynP(K,EL_IN,IRK)*PHI_EDGE(K,I,LED,pa,J)
#endif

#ifdef SED_LAY
               do ll = 1,layers
                  bed_IN(ll) = bed_IN(ll) +
     &                  bed(K,EL_IN,IRK,ll)*PHI_EDGE(K,I,LED,pa,J)
                HB_IN=HB_IN + bed(k,EL_IN,irk,ll)*PHI_EDGE(K,I,LED,pa,J)

                  MZ_X_IN(ll) = MZ_X_IN(ll) +
     & MZ(K,1,ll,EL_IN)*PHI_EDGE(K,I,LED,pa,J)
                  MZ_Y_IN(ll) = MZ_Y_IN(ll) +
     & MZ(K,2,ll,EL_IN)*PHI_EDGE(K,I,LED,pa,J)

               enddo
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
            iota_EX = 0.D0 !iota_IN
            !print*, 'test'
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

C.....Compute the Roe flux

            CALL NUMERICAL_FLUX(IT,test_el)

C.....Add LDG terms for viscosity

#ifdef WAVE_DIF
            F_HAT = F_HAT + HZ_X_IN*NX*SFAC_IN + HZ_Y_IN*NY
#endif
            G_HAT = G_HAT + LZ_XX_IN*NX*SFAC_IN + LZ_XY_IN*NY
            H_HAT = H_HAT + LZ_YX_IN*NX*SFAC_IN + LZ_YY_IN*NY

#ifdef TRACE
            I_HAT = I_HAT + TZ_X_IN*NX*SFAC_IN + TZ_Y_IN*NY
#endif

C.....Add LDG terms for sediment

#ifdef SED_LAY
            do ll=1,layers
               bed_HAT(ll) = bed_HAT(ll) +
     & MZ_X_IN(ll)*NX*SFAC_IN + MZ_Y_IN(ll)*NY
            enddo
#endif
            
C.....Compute the edge integral
C.....(modified for wetting and drying) 

C.....Compute the edge integral
            DO K = 1,DOFS(EL_IN)

               if (EDN(EL_IN) .eq. 3) then
                  W_IN = 2.0d0*M_INV(K,pa)/AREAS(EL_IN)*XLEN(GED)*
     &               PHI_EDGE(K,I,LED,pa,J)*WEGP(I,pa)
               elseif (EDN(EL_IN) .eq. 4) then
                  do kk = 1,DOFS(EL_IN)
                     vtemp1(kk) = M_INVM(EL_IN,K,kk,pa)
                     vtemp2(kk) = PHI_EDGE(kk,I,LED,pa,J)
                  enddo
                  W_IN = dot_product(vtemp1,vtemp2)
                  W_IN = W_IN*0.5d0*XLEN(GED)*WEGP(I,pa)
               endif

               RHS_ZE(K,EL_IN,IRK) = RHS_ZE(K,EL_IN,IRK) - W_IN*F_HAT
               RHS_QX(K,EL_IN,IRK) = RHS_QX(K,EL_IN,IRK) - W_IN*G_HAT
               RHS_QY(K,EL_IN,IRK) = RHS_QY(K,EL_IN,IRK) - W_IN*H_HAT

#ifdef TRACE
               RHS_iota(K,EL_IN,IRK) = RHS_iota(K,EL_IN,IRK) - W_IN*I_HAT
#endif

#ifdef CHEM
               RHS_iota(K,EL_IN,IRK) = RHS_iota(K,EL_IN,IRK) - W_IN*I_HAT
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

         ENDDO
         if (J .eq. 2) deallocate( vtemp1,vtemp2 )
 1000 CONTINUE
      
      RETURN
      END SUBROUTINE


C***********************************************************************
C     
C     SUBROUTINE LAND_EDGE_HYDRO_POST( )
C     
C     This subroutine does the following:
C     
C     1.  Direct velocity at each node on land edges toward
C     the tangential direction of the edge.
C     (This procedure needs to be called after the wetting and 
C     drying post-process.)
C     
C     Written by Shintaro Bunya (02-13-2006)
C     
C***********************************************************************

C      SUBROUTINE LAND_EDGE_HYDRO_POST()
C
C.....Use appropriate modules
C
C      USE SIZES,ONLY : SZ,MYPROC,layers
C      USE GLOBAL,ONLY : pdg_el 
C      USE DG,ONLY : NLEDS,NLEDN,NEDSD,NEDEL,COSNX,SINNX,QX,QY,
C     &     PHI_CORNER, IRK, EL_UPDATED,pa,DOFS
C
C      IMPLICIT NONE
C
C.....Declare local variables
C
C      INTEGER L, LED, GED, EL_IN, K, KK, NOD1, NOD2, NOD3, NEDGES
C      REAL(SZ) NX, NY
C      REAL(SZ) TX, TY
C      REAL(SZ) DIR
C      REAL(SZ) QXP(3),QXPN(3)
C      REAL(SZ) QYP(3),QYPN(3)
C      REAL(SZ) QP(3)
C      
C      DO 1000 L=1, NLEDS
C
C.....Retrieve the global and local edge number
C
C         GED = NLEDN(L)
C         LED = NEDSD(1,GED)
C
C.....Retrieve the elements which share the edge
C
C         EL_IN = NEDEL(1,GED)
C
C         PA = PDG_EL(GED)
C
C#ifdef P0
C         if (pa.eq.0) then
C            pa = 1
C         endif
C#endif
C
C.....Apply the following sequences only on dry, drying, or wetting elements
C
C         IF(EL_UPDATED(EL_IN).EQ.0) CYCLE
C
C.....Retrieve the components of the normal vector to the edge
C         
C         NX = COSNX(GED)
C         NY = SINNX(GED)
C         
C.....Set the components for the tangential vector to the edge
C         
C         TX = -NY
C         TY =  NX
C
C.....Compute nodal QX and QY
C
C         QXP(:) = 0.D0
C         QYP(:) = 0.D0
C
C         DO K = 1,3
C            DO KK = 1,DOFS(EL_IN)
C               QXP(K) = QXP(K) + PHI_CORNER(KK,K,pa)*QX(KK,EL_IN,IRK+1)
C               QYP(K) = QYP(K) + PHI_CORNER(KK,K,pa)*QY(KK,EL_IN,IRK+1)
C            ENDDO
C            QP(K) = SQRT(QXP(K)*QXP(K)+QYP(K)*QYP(K))
C            
C         ENDDO
C
C.....Retrieve the end nodes of the edge and the other node
C
C         NOD1 = MOD(LED+0,3)+1
C         NOD2 = MOD(LED+1,3)+1
C         NOD3 = MOD(LED+2,3)+1
C         
C.....Modify QX and QY at the end nodes

C.......Node 1

C         IF(ABS( QXP(NOD1)*TX + QYP(NOD1)*TY ).GT.(10.D0**(-20.D0)))THEN
C            DIR = ( QXP(NOD1)*TX + QYP(NOD1)*TY )
C     &           / ABS( QXP(NOD1)*TX + QYP(NOD1)*TY )
C
C#if 1
C            QXPN(NOD1) = DIR*QP(NOD1)*TX
C            QYPN(NOD1) = DIR*QP(NOD1)*TY
C#else
C            QXPN(NOD1) = (QXP(NOD1)*TX + QYP(NOD1)*TY)*TX
C            QYPN(NOD1) = (QXP(NOD1)*TX + QYP(NOD1)*TY)*TY
C#endif
C         ELSE
C            QXPN(NOD1) = 0.D0
C            QYPN(NOD1) = 0.D0
C         ENDIF
C
C.......Node 2
C
C         IF(ABS( QXP(NOD2)*TX + QYP(NOD2)*TY ).GT.(10.D0**(-20.D0)))THEN
C            DIR = ( QXP(NOD2)*TX + QYP(NOD2)*TY )
C     &           / ABS( QXP(NOD2)*TX + QYP(NOD2)*TY )
C
C#if 1
C            QXPN(NOD2) = DIR*QP(NOD2)*TX
C            QYPN(NOD2) = DIR*QP(NOD2)*TY
C#else
C            QXPN(NOD2) = (QXP(NOD2)*TX + QYP(NOD2)*TY)*TX
C            QYPN(NOD2) = (QXP(NOD2)*TX + QYP(NOD2)*TY)*TY
C#endif
C         ELSE
C            QXPN(NOD2) = 0.D0
C            QYPN(NOD2) = 0.D0
C         ENDIF
C
C.......Node 3

C         QXPN(NOD3) = QXP(NOD3)
C         QYPN(NOD3) = QYP(NOD3)

C.....Inverse the modified QX and QY

C         QX(1,EL_IN,IRK+1) = 1.D0/3.D0*(QXPN(1)+QXPN(2)+QXPN(3))
C       QX(2,EL_IN,IRK+1) =-1.D0/6.D0*(QXPN(1)+QXPN(2))+1.D0/3.D0*QXPN(3)
C         QX(3,EL_IN,IRK+1) =-0.5D0*QXPN(1)+0.5D0*QXPN(2)
C
C         QY(1,EL_IN,IRK+1) = 1.D0/3.D0*(QYPN(1)+QYPN(2)+QYPN(3))
C       QY(2,EL_IN,IRK+1) =-1.D0/6.D0*(QYPN(1)+QYPN(2))+1.D0/3.D0*QYPN(3)
C         QY(3,EL_IN,IRK+1) =-0.5D0*QYPN(1)+0.5D0*QYPN(2)
C
C 1000 CONTINUE
C      RETURN
C      END SUBROUTINE
