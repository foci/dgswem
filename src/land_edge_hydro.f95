!***********************************************************************
!     
!     SUBROUTINE LAND_EDGE_HYDRO( )
!     
!     This subroutine does the following:
!     
!     1.  Calculates the values of the necessary variables at the edge
!     gauss points for NO-NORMAL FLOW edges
!     2.  Calls the appropriate subroutine to compute the flux at
!     these points.
!     3.  Calls the appropriate subroutine to compute the boundary
!     integrals.
!     
!     Written by Ethan Kubatko (06-11-2004)
!     
!-----------------------------------------------------------------------
!     
!     01-02-2007, sb, Modified for LDG
!     08-xx-2005, sb, Modifications for wetting/drying
!     
!     01-10-2011 - cem - adapted for p_enrichment and multicomponent
!     06-02-2012 - cem - adapted for sediment
!     
!***********************************************************************

      SUBROUTINE LAND_EDGE_HYDRO(s,IT)

!.....Use appropriate modules

      USE GLOBAL
      USE DG
      Use SIZES

      IMPLICIT NONE

      type (sizes_type) :: s

!.....Declare local variables

      INTEGER L, LED, GED, i,k,jj,II,ll,IT,mm
      Real(SZ) DEN2,U_AVG,V_AVG,VEL_NORMAL,q_RoeX, q_RoeY, q_Roe
      REAL(SZ) TX, TY, DEN,ell_1,ell_2,ell_3,HZ_X_IN,HZ_Y_IN
      REAL(SZ) LZ_XX_IN, LZ_XY_IN, LZ_YX_IN, LZ_YY_IN,W_IN
      Real(SZ) MZ_X_IN(s%layers),MZ_Y_IN(s%layers),TZ_X_IN,TZ_Y_IN

      dg%test_el = 0
      DO 1000 L = 1,dg%NLEDS
         
!.....Retrieve the global and local edge number

         GED = dg%NLEDN(L)
         LED = dg%NEDSD(1,GED)

!.....Retrieve the elements which share the edge

         EL_IN = dg%NEDEL(1,GED)

         dg%pa = PDG_EL(EL_IN)

#ifdef P0
         if (dg%pa.eq.0) then
            dg%pa = 1
         endif
#endif
         
!.....If the element is dry then skip the edge calculation

         IF (dg%WDFLG(EL_IN).EQ.0) GOTO 1000

         dg%test_el = dg%test_el+1

!.....Retrieve the components of the normal vector to the edge
         
         dg%NX = dg%COSNX(GED)
         dg%NY = dg%SINNX(GED)
         
!.....Set the components for the tangential vector to the edge
         
         TX = -dg%NY
         TY =  dg%NX

!.....Compute dg%ZE, `dg%QX, dg%QY, and dg%HB at each Gauss point

         DO I = 1,dg%NEGP(dg%pa)
            
            dg%ZE_IN = dg%ZE(1,EL_IN,dg%IRK)
            dg%QX_IN = dg%QX(1,EL_IN,dg%IRK)
            dg%QY_IN = dg%QY(1,EL_IN,dg%IRK)
            dg%HB_IN = dg%BATHED(I,LED,EL_IN,dg%pa)

            dg%SFAC_IN = dg%SFACED(I,LED,EL_IN,dg%pa)

#ifdef TRACE
            dg%iota_IN = dg%iota(1,EL_IN,dg%IRK)
#endif

#ifdef CHEM
            dg%iota_IN = 0.D0 !dg%iota(1,EL_IN,dg%IRK)
            dg%iota2_IN = 0.D0 !dg%iota2(1,EL_IN,dg%IRK)
#endif

#ifdef dynp
            dynP_IN = dg%dynP(1,EL_IN,dg%IRK)
#endif

            !When layered, these change
#ifdef SED_LAY
            dg%HB(1,EL_IN,dg%irk) = 0.D0
            do ll=1,s%layers
               dg%HB(1,EL_IN,dg%irk) = dg%HB(1,EL_IN,dg%irk) + dg%bed(1,EL_IN,dg%irk,ll)

               MZ_X_IN(ll) =  dg%MZ(1,1,ll,EL_IN)
               MZ_Y_IN(ll) =  dg%MZ(1,2,ll,EL_IN)
            enddo
            dg%bed_IN(:) = dg%bed(1,EL_IN,dg%irk,:)
            dg%HB_IN = dg%HB(1,EL_IN,dg%irk)

#endif

            !...do it and do it again ... 
#ifdef WAVE_DIF
            HZ_X_IN = dg%HZ(1,1,1,EL_IN)
            HZ_Y_IN = dg%HZ(1,2,2,EL_IN)
#endif
            LZ_XX_IN = dg%LZ(1,1,1,EL_IN)
            LZ_XY_IN = dg%LZ(1,1,2,EL_IN)
            LZ_YX_IN = dg%LZ(1,2,1,EL_IN)
            LZ_YY_IN = dg%LZ(1,2,2,EL_IN)

#ifdef TRACE

            TZ_X_IN = dg%TZ(1,1,1,EL_IN)
            TZ_Y_IN = dg%TZ(1,2,2,EL_IN)

#endif

            !Compute sediment diffusion contribution

!.....Compute the solution at the interior state
!.....(modified for wetting and drying)

            DO K = 2,dg%DOFS(EL_IN)
               
               dg%ZE_IN = dg%ZE_IN + dg%ZE(K,EL_IN,dg%IRK)*dg%PHI_EDGE(K,I,LED,dg%pa)
               dg%QX_IN = dg%QX_IN + dg%QX(K,EL_IN,dg%IRK)*dg%PHI_EDGE(K,I,LED,dg%pa)
               dg%QY_IN = dg%QY_IN + dg%QY(K,EL_IN,dg%IRK)*dg%PHI_EDGE(K,I,LED,dg%pa)

                                !LDG terms
#ifdef WAVE_DIF
               HZ_X_IN = HZ_X_IN + dg%HZ(K,1,1,EL_IN)*dg%PHI_EDGE(K,I,LED,dg%pa)
               HZ_Y_IN = HZ_Y_IN + dg%HZ(K,2,2,EL_IN)*dg%PHI_EDGE(K,I,LED,dg%pa)
#endif
               LZ_XX_IN = LZ_XX_IN + dg%LZ(K,1,1,EL_IN)*dg%PHI_EDGE(K,I,LED,dg%pa)
               LZ_XY_IN = LZ_XY_IN + dg%LZ(K,1,2,EL_IN)*dg%PHI_EDGE(K,I,LED,dg%pa)
               LZ_YX_IN = LZ_YX_IN + dg%LZ(K,2,1,EL_IN)*dg%PHI_EDGE(K,I,LED,dg%pa)
               LZ_YY_IN = LZ_YY_IN + dg%LZ(K,2,2,EL_IN)*dg%PHI_EDGE(K,I,LED,dg%pa)

#ifdef TRACE

               TZ_X_IN = TZ_X_IN + dg%TZ(K,1,1,EL_IN)*dg%PHI_EDGE(K,I,LED,dg%pa)
               TZ_Y_IN = TZ_Y_IN + dg%TZ(K,2,2,EL_IN)*dg%PHI_EDGE(K,I,LED,dg%pa)

#endif

                                !LDG terms for sediment diffusion

#ifdef TRACE
               dg%iota_IN = dg%iota_IN + dg%iota(K,EL_IN,dg%IRK)*dg%PHI_EDGE(K,I,LED,dg%pa)
#endif

#ifdef CHEM
               dg%iota_IN = dg%iota_IN + dg%iota(K,EL_IN,dg%IRK)*dg%PHI_EDGE(K,I,LED,dg%pa)
               dg%iota2_IN = dg%iota2_IN + dg%iota2(K,EL_IN,dg%IRK)*dg%PHI_EDGE(K,I,LED,dg%pa)
#endif

#ifdef dynp
               dynP_IN = dynP_IN + dg%dynP(K,EL_IN,dg%IRK)*dg%PHI_EDGE(K,I,LED,dg%pa)
#endif

#ifdef SED_LAY
               do ll = 1,s%layers
                  dg%bed_IN(ll) = dg%bed_IN(ll) + dg%bed(K,EL_IN,dg%IRK,ll)*dg%PHI_EDGE(K,I,LED,dg%pa)
                  dg%HB_IN = dg%HB_IN + dg%bed(k,EL_IN,dg%irk,ll)*dg%PHI_EDGE(K,I,LED,dg%pa)

                  MZ_X_IN(ll) = MZ_X_IN(ll) + dg%MZ(K,1,ll,EL_IN)*dg%PHI_EDGE(K,I,LED,dg%pa)
                  MZ_Y_IN(ll) = MZ_Y_IN(ll) + dg%MZ(K,2,ll,EL_IN)*dg%PHI_EDGE(K,I,LED,dg%pa)

               enddo
#endif

            ENDDO


!.....Compute the velocity in the normal and tangental direction
      
            Q_N_INT = dg%QX_IN*dg%NX + dg%QY_IN*dg%NY
            Q_T_INT = dg%QX_IN*TX + dg%QY_IN*TY

!.....Reflect the velocity in the normal direction

            Q_N_EXT = -Q_N_INT
            Q_T_EXT =  Q_T_INT
            
!.....Compute the x and y components of the external state flow

 
            DEN = 1.D0/(dg%NX*TY - dg%NY*TX)
            dg%QX_EX = ( TY*Q_N_EXT - dg%NY*Q_T_EXT)*DEN
            dg%QY_EX = (-TX*Q_N_EXT + dg%NX*Q_T_EXT)*DEN

            dg%ZE_EX = dg%ZE_IN
            dg%HB_EX = dg%HB_IN
            dg%SFAC_EX = dg%SFAC_IN

#ifdef TRACE
            dg%iota_EX = 0.D0 !dg%iota_IN
            !print*, 'test'
#endif

#ifdef CHEM
            dg%iota_EX = dg%iota_IN
            dg%iota2_EX = dg%iota2_IN
#endif

#ifdef dynp
            dynP_EX = dynP_IN
#endif

#ifdef SED_LAY
            dg%bed_EX(:) = dg%bed_IN(:)
#endif

!.....Compute the Roe flux

            CALL NUMERICAL_FLUX(s,IT,dg%test_el)

!.....Add LDG terms for viscosity

#ifdef WAVE_DIF
            F_HAT = F_HAT + HZ_X_IN*dg%NX*dg%SFAC_IN + HZ_Y_IN*dg%NY
#endif
            G_HAT = G_HAT + LZ_XX_IN*dg%NX*dg%SFAC_IN + LZ_XY_IN*dg%NY
            H_HAT = H_HAT + LZ_YX_IN*dg%NX*dg%SFAC_IN + LZ_YY_IN*dg%NY

#ifdef TRACE
            I_HAT = I_HAT + TZ_X_IN*dg%NX*dg%SFAC_IN + TZ_Y_IN*dg%NY
#endif

!.....Add LDG terms for sediment

#ifdef SED_LAY
            do ll=1,s%layers
               dg%bed_HAT(ll) = dg%bed_HAT(ll) + MZ_X_IN(ll)*dg%NX*dg%SFAC_IN + MZ_Y_IN(ll)*dg%NY
            enddo
#endif
            
!.....Compute the edge integral
!.....(modified for wetting and drying) 

!.....Compute the edge integral
            DO K = 1,dg%DOFS(EL_IN)

               W_IN = 2.0*dg%M_INV(K,dg%pa)/AREAS(EL_IN)*dg%XLEN(GED)*&
              dg%PHI_EDGE(K,I,LED,dg%pa)*dg%WEGP(I,dg%pa)

               dg%RHS_ZE(K,EL_IN,dg%IRK) = dg%RHS_ZE(K,EL_IN,dg%IRK) - W_IN*F_HAT
               dg%RHS_QX(K,EL_IN,dg%IRK) = dg%RHS_QX(K,EL_IN,dg%IRK) - W_IN*G_HAT
               dg%RHS_QY(K,EL_IN,dg%IRK) = dg%RHS_QY(K,EL_IN,dg%IRK) - W_IN*H_HAT

#ifdef TRACE
               dg%RHS_iota(K,EL_IN,dg%IRK) = dg%RHS_iota(K,EL_IN,dg%IRK) - W_IN*I_HAT
#endif

#ifdef CHEM
               dg%RHS_iota(K,EL_IN,dg%IRK) = dg%RHS_iota(K,EL_IN,dg%IRK) - W_IN*I_HAT
               dg%RHS_iota2(K,EL_IN,dg%IRK) = dg%RHS_iota2(K,EL_IN,dg%IRK) - W_IN*J_HAT
#endif

#ifdef dynp
               dg%RHS_dynP(K,EL_IN,dg%IRK) = dg%RHS_dynP(K,EL_IN,dg%IRK) - W_IN*K_HAT
#endif

#ifdef SED_LAY
               do ll = 1,s%layers
                  dg%RHS_bed(K,EL_IN,dg%IRK,ll) = dg%RHS_bed(K,EL_IN,dg%IRK,ll) - W_IN*dg%bed_HAT(ll)
               enddo
#endif

            ENDDO

         ENDDO

 1000 CONTINUE
      
      RETURN
      END SUBROUTINE


!***********************************************************************
!     
!     SUBROUTINE LAND_EDGE_HYDRO_POST( )
!     
!     This subroutine does the following:
!     
!     1.  Direct velocity at each node on land edges toward
!     the tangential direction of the edge.
!     (This procedure dg%needs to be called after the wetting and 
!     drying post-process.)
!     
!     Written by Shintaro Bunya (02-13-2006)
!     
!***********************************************************************

      SUBROUTINE LAND_EDGE_HYDRO_POST(s)

!.....Use appropriate modules

      USE SIZES
      USE GLOBAL,ONLY : pdg_el 
      USE DG,ONLY : NLEDS,NLEDN,NEDSD,NEDEL,COSNX,SINNX,QX,QY,&
     PHI_CORNER, IRK, EL_UPDATED,pa,DOFS

      IMPLICIT NONE

      type (sizes_type) :: s
!.....Declare local variables

      INTEGER L, LED, GED, EL_IN, K, KK, NOD1, NOD2, NOD3, NEDGES
      REAL(SZ) NX, NY
      REAL(SZ) TX, TY
      REAL(SZ) DIR
      REAL(SZ) QXP(3),QXPN(3)
      REAL(SZ) QYP(3),QYPN(3)
      REAL(SZ) QP(3)
      
      DO 1000 L=1, NLEDS

!.....Retrieve the global and local edge number

         GED = NLEDN(L)
         LED = NEDSD(1,GED)

!.....Retrieve the elements which share the edge

         EL_IN = NEDEL(1,GED)

         PA = PDG_EL(GED)

#ifdef P0
         if (pa.eq.0) then
            pa = 1
         endif
#endif

!.....Apply the following sequences only on dry, drying, or wetting elements

         IF(EL_UPDATED(EL_IN).EQ.0) CYCLE

!.....Retrieve the components of the normal vector to the edge
         
         NX = COSNX(GED)
         NY = SINNX(GED)
         
!.....Set the components for the tangential vector to the edge
         
         TX = -NY
         TY =  NX

!.....Compute nodal QX and QY

         QXP(:) = 0.D0
         QYP(:) = 0.D0

         DO K = 1,3
            DO KK = 1,DOFS(EL_IN)
               QXP(K) = QXP(K) + PHI_CORNER(KK,K,pa)*QX(KK,EL_IN,IRK+1)
               QYP(K) = QYP(K) + PHI_CORNER(KK,K,pa)*QY(KK,EL_IN,IRK+1)
            ENDDO
            QP(K) = SQRT(QXP(K)*QXP(K)+QYP(K)*QYP(K))
            
         ENDDO

!.....Retrieve the end nodes of the edge and the other node

         NOD1 = MOD(LED+0,3)+1
         NOD2 = MOD(LED+1,3)+1
         NOD3 = MOD(LED+2,3)+1
         
!.....Modify QX and QY at the end nodes

!.......Node 1

         IF(ABS( QXP(NOD1)*TX + QYP(NOD1)*TY ).GT.(10.D0**(-20.D0))) THEN
            DIR = ( QXP(NOD1)*TX + QYP(NOD1)*TY )&
           / ABS( QXP(NOD1)*TX + QYP(NOD1)*TY )

#if 1
            QXPN(NOD1) = DIR*QP(NOD1)*TX
            QYPN(NOD1) = DIR*QP(NOD1)*TY
#else
            QXPN(NOD1) = (QXP(NOD1)*TX + QYP(NOD1)*TY)*TX
            QYPN(NOD1) = (QXP(NOD1)*TX + QYP(NOD1)*TY)*TY
#endif
         ELSE
            QXPN(NOD1) = 0.D0
            QYPN(NOD1) = 0.D0
         ENDIF

!.......Node 2

         IF(ABS( QXP(NOD2)*TX + QYP(NOD2)*TY ).GT.(10.D0**(-20.D0))) THEN
            DIR = ( QXP(NOD2)*TX + QYP(NOD2)*TY )&
           / ABS( QXP(NOD2)*TX + QYP(NOD2)*TY )

#if 1
            QXPN(NOD2) = DIR*QP(NOD2)*TX
            QYPN(NOD2) = DIR*QP(NOD2)*TY
#else
            QXPN(NOD2) = (QXP(NOD2)*TX + QYP(NOD2)*TY)*TX
            QYPN(NOD2) = (QXP(NOD2)*TX + QYP(NOD2)*TY)*TY
#endif
         ELSE
            QXPN(NOD2) = 0.D0
            QYPN(NOD2) = 0.D0
         ENDIF

!.......Node 3

         QXPN(NOD3) = QXP(NOD3)
         QYPN(NOD3) = QYP(NOD3)

!.....Inverse the modified QX and QY

         QX(1,EL_IN,IRK+1) = 1.D0/3.D0*(QXPN(1)+QXPN(2)+QXPN(3))
         QX(2,EL_IN,IRK+1) =-1.D0/6.D0*(QXPN(1)+QXPN(2))+1.D0/3.D0*QXPN(3)
         QX(3,EL_IN,IRK+1) =-0.5D0*QXPN(1)+0.5D0*QXPN(2)

         QY(1,EL_IN,IRK+1) = 1.D0/3.D0*(QYPN(1)+QYPN(2)+QYPN(3))
         QY(2,EL_IN,IRK+1) =-1.D0/6.D0*(QYPN(1)+QYPN(2))+1.D0/3.D0*QYPN(3)
         QY(3,EL_IN,IRK+1) =-0.5D0*QYPN(1)+0.5D0*QYPN(2)

 1000 CONTINUE
      RETURN
      END SUBROUTINE
