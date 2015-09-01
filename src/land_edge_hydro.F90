!***********************************************************************
!     
!     SUBROUTINE LAND_EDGE_HYDRO( )
!     
!     This subroutine does the following:
!     
!     1.  Calculates the values of the necessary variables at the global_here%edge
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

      SUBROUTINE LAND_EDGE_HYDRO(s,dg_here,global_here,IT)

!.....Use appropriate modules

      USE GLOBAL
      USE DG
      Use SIZES

      IMPLICIT NONE

      type (sizes_type) :: s
      type (dg_type) :: dg_here
      type (global_type) :: global_here

!.....Declare local variables

      INTEGER L, LED, GED, i,k,jj,II,ll,IT,mm
      Real(SZ) DEN2,U_AVG,V_AVG,VEL_NORMAL,q_RoeX, q_RoeY, q_Roe
      REAL(SZ) TX, TY, DEN,ell_1,ell_2,ell_3,HZ_X_IN,HZ_Y_IN
      REAL(SZ) LZ_XX_IN, LZ_XY_IN, LZ_YX_IN, LZ_YY_IN,W_IN
      Real(SZ) MZ_X_IN(s%layers),MZ_Y_IN(s%layers),TZ_X_IN,TZ_Y_IN

      dg_here%test_el = 0
      DO 1000 L = 1,dg_here%NLEDS
         
!.....Retrieve the global and local global_here%edge number

         GED = dg_here%NLEDN(L)
         LED = dg_here%NEDSD(1,GED)

!.....Retrieve the elements which share the global_here%edge

         global_here%EL_IN = dg_here%NEDEL(1,GED)

         dg_here%pa = global_here%PDG_EL(global_here%EL_IN)

#ifdef P0
         if (dg_here%pa.eq.0) then
            dg_here%pa = 1
         endif
#endif
         
!.....If the element is dry then skip the global_here%edge calculation

         IF (dg_here%WDFLG(global_here%EL_IN).EQ.0) GOTO 1000

         dg_here%test_el = dg_here%test_el+1

!.....Retrieve the components of the normal vector to the global_here%edge
         
         dg_here%NX = dg_here%COSNX(GED)
         dg_here%NY = dg_here%SINNX(GED)
         
!.....Set the components for the tangential vector to the global_here%edge
         
         TX = -dg_here%NY
         TY =  dg_here%NX

!.....Compute dg_here%ZE, `dg_here%QX, dg_here%QY, and dg_here%HB at each Gauss point

         DO I = 1,dg_here%NEGP(dg_here%pa)
            
            dg_here%ZE_IN = dg_here%ZE(1,global_here%EL_IN,dg_here%IRK)
            dg_here%QX_IN = dg_here%QX(1,global_here%EL_IN,dg_here%IRK)
            dg_here%QY_IN = dg_here%QY(1,global_here%EL_IN,dg_here%IRK)
            dg_here%HB_IN = dg_here%BATHED(I,LED,global_here%EL_IN,dg_here%pa)

            dg_here%SFAC_IN = dg_here%SFACED(I,LED,global_here%EL_IN,dg_here%pa)

#ifdef TRACE
            dg_here%iota_IN = dg_here%iota(1,global_here%EL_IN,dg_here%IRK)
#endif

#ifdef CHEM
            dg_here%iota_IN = 0.D0 !dg_here%iota(1,global_here%EL_IN,dg_here%IRK)
            dg_here%iota2_IN = 0.D0 !dg_here%iota2(1,global_here%EL_IN,dg_here%IRK)
#endif

#ifdef dynp
            dynP_IN = dg_here%dynP(1,global_here%EL_IN,dg_here%IRK)
#endif

            !When layered, these change
#ifdef SED_LAY
            dg_here%HB(1,global_here%EL_IN,dg_here%irk) = 0.D0
            do ll=1,s%layers
               dg_here%HB(1,global_here%EL_IN,dg_here%irk) = dg_here%HB(1,global_here%EL_IN,dg_here%irk) + dg_here%bed(1,global_here%EL_IN,dg_here%irk,ll)

               MZ_X_IN(ll) =  dg_here%MZ(1,1,ll,global_here%EL_IN)
               MZ_Y_IN(ll) =  dg_here%MZ(1,2,ll,global_here%EL_IN)
            enddo
            dg_here%bed_IN(:) = dg_here%bed(1,global_here%EL_IN,dg_here%irk,:)
            dg_here%HB_IN = dg_here%HB(1,global_here%EL_IN,dg_here%irk)

#endif

            !...do it and do it again ... 
#ifdef WAVE_DIF
            HZ_X_IN = dg_here%HZ(1,1,1,global_here%EL_IN)
            HZ_Y_IN = dg_here%HZ(1,2,2,global_here%EL_IN)
#endif
            LZ_XX_IN = dg_here%LZ(1,1,1,global_here%EL_IN)
            LZ_XY_IN = dg_here%LZ(1,1,2,global_here%EL_IN)
            LZ_YX_IN = dg_here%LZ(1,2,1,global_here%EL_IN)
            LZ_YY_IN = dg_here%LZ(1,2,2,global_here%EL_IN)

#ifdef TRACE

            TZ_X_IN = dg_here%TZ(1,1,1,global_here%EL_IN)
            TZ_Y_IN = dg_here%TZ(1,2,2,global_here%EL_IN)

#endif

            !Compute sediment diffusion contribution

!.....Compute the solution at the interior state
!.....(modified for wetting and drying)

            DO K = 2,dg_here%DOFS(global_here%EL_IN)
               
               dg_here%ZE_IN = dg_here%ZE_IN + dg_here%ZE(K,global_here%EL_IN,dg_here%IRK)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
               dg_here%QX_IN = dg_here%QX_IN + dg_here%QX(K,global_here%EL_IN,dg_here%IRK)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
               dg_here%QY_IN = dg_here%QY_IN + dg_here%QY(K,global_here%EL_IN,dg_here%IRK)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)

                                !LDG terms
#ifdef WAVE_DIF
               HZ_X_IN = HZ_X_IN + dg_here%HZ(K,1,1,global_here%EL_IN)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
               HZ_Y_IN = HZ_Y_IN + dg_here%HZ(K,2,2,global_here%EL_IN)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
#endif
               LZ_XX_IN = LZ_XX_IN + dg_here%LZ(K,1,1,global_here%EL_IN)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
               LZ_XY_IN = LZ_XY_IN + dg_here%LZ(K,1,2,global_here%EL_IN)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
               LZ_YX_IN = LZ_YX_IN + dg_here%LZ(K,2,1,global_here%EL_IN)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
               LZ_YY_IN = LZ_YY_IN + dg_here%LZ(K,2,2,global_here%EL_IN)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)

#ifdef TRACE

               TZ_X_IN = TZ_X_IN + dg_here%TZ(K,1,1,global_here%EL_IN)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
               TZ_Y_IN = TZ_Y_IN + dg_here%TZ(K,2,2,global_here%EL_IN)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)

#endif

                                !LDG terms for sediment diffusion

#ifdef TRACE
               dg_here%iota_IN = dg_here%iota_IN + dg_here%iota(K,global_here%EL_IN,dg_here%IRK)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
#endif

#ifdef CHEM
               dg_here%iota_IN = dg_here%iota_IN + dg_here%iota(K,global_here%EL_IN,dg_here%IRK)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
               dg_here%iota2_IN = dg_here%iota2_IN + dg_here%iota2(K,global_here%EL_IN,dg_here%IRK)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
#endif

#ifdef dynp
               dynP_IN = dynP_IN + dg_here%dynP(K,global_here%EL_IN,dg_here%IRK)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
#endif

#ifdef SED_LAY
               do ll = 1,s%layers
                  dg_here%bed_IN(ll) = dg_here%bed_IN(ll) + dg_here%bed(K,global_here%EL_IN,dg_here%IRK,ll)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
                  dg_here%HB_IN = dg_here%HB_IN + dg_here%bed(k,global_here%EL_IN,dg_here%irk,ll)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)

                  MZ_X_IN(ll) = MZ_X_IN(ll) + dg_here%MZ(K,1,ll,global_here%EL_IN)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
                  MZ_Y_IN(ll) = MZ_Y_IN(ll) + dg_here%MZ(K,2,ll,global_here%EL_IN)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)

               enddo
#endif

            ENDDO


!.....Compute the velocity in the normal and tangental direction
      
            global_here%Q_N_INT = dg_here%QX_IN*dg_here%NX + dg_here%QY_IN*dg_here%NY
            global_here%Q_T_INT = dg_here%QX_IN*TX + dg_here%QY_IN*TY

!.....Reflect the velocity in the normal direction

            global_here%Q_N_EXT = -global_here%Q_N_INT
            global_here%Q_T_EXT =  global_here%Q_T_INT
            
!.....Compute the global_here%x and global_here%y components of the external state flow

 
            DEN = 1.D0/(dg_here%NX*TY - dg_here%NY*TX)
            dg_here%QX_EX = ( TY*global_here%Q_N_EXT - dg_here%NY*global_here%Q_T_EXT)*DEN
            dg_here%QY_EX = (-TX*global_here%Q_N_EXT + dg_here%NX*global_here%Q_T_EXT)*DEN

            dg_here%ZE_EX = dg_here%ZE_IN
            dg_here%HB_EX = dg_here%HB_IN
            dg_here%SFAC_EX = dg_here%SFAC_IN

#ifdef TRACE
            dg_here%iota_EX = 0.D0 !dg_here%iota_IN
            !print*, 'test'
#endif

#ifdef CHEM
            dg_here%iota_EX = dg_here%iota_IN
            dg_here%iota2_EX = dg_here%iota2_IN
#endif

#ifdef dynp
            dynP_EX = dynP_IN
#endif

#ifdef SED_LAY
            dg_here%bed_EX(:) = dg_here%bed_IN(:)
#endif

!.....Compute the Roe flux

            CALL NUMERICAL_FLUX(s,dg_here,global_here,IT,dg_here%test_el)

!.....Add LDG terms for viscosity

#ifdef WAVE_DIF
            global_here%F_HAT = global_here%F_HAT + HZ_X_IN*dg_here%NX*dg_here%SFAC_IN + HZ_Y_IN*dg_here%NY
#endif
            global_here%G_HAT = global_here%G_HAT + LZ_XX_IN*dg_here%NX*dg_here%SFAC_IN + LZ_XY_IN*dg_here%NY
            global_here%H_HAT = global_here%H_HAT + LZ_YX_IN*dg_here%NX*dg_here%SFAC_IN + LZ_YY_IN*dg_here%NY

#ifdef TRACE
            global_here%I_HAT = global_here%I_HAT + TZ_X_IN*dg_here%NX*dg_here%SFAC_IN + TZ_Y_IN*dg_here%NY
#endif

!.....Add LDG terms for sediment

#ifdef SED_LAY
            do ll=1,s%layers
               dg_here%bed_HAT(ll) = dg_here%bed_HAT(ll) + MZ_X_IN(ll)*dg_here%NX*dg_here%SFAC_IN + MZ_Y_IN(ll)*dg_here%NY
            enddo
#endif
            
!.....Compute the global_here%edge integral
!.....(modified for wetting and drying) 

!.....Compute the global_here%edge integral
            DO K = 1,dg_here%DOFS(global_here%EL_IN)

               W_IN = 2.0*dg_here%M_INV(K,dg_here%pa)/global_here%AREAS(global_here%EL_IN)*dg_here%XLEN(GED)*&
              dg_here%PHI_EDGE(K,I,LED,dg_here%pa)*dg_here%WEGP(I,dg_here%pa)

               dg_here%RHS_ZE(K,global_here%EL_IN,dg_here%IRK) = dg_here%RHS_ZE(K,global_here%EL_IN,dg_here%IRK) - W_IN*global_here%F_HAT
               dg_here%RHS_QX(K,global_here%EL_IN,dg_here%IRK) = dg_here%RHS_QX(K,global_here%EL_IN,dg_here%IRK) - W_IN*global_here%G_HAT
               dg_here%RHS_QY(K,global_here%EL_IN,dg_here%IRK) = dg_here%RHS_QY(K,global_here%EL_IN,dg_here%IRK) - W_IN*global_here%H_HAT

#ifdef TRACE
               dg_here%RHS_iota(K,global_here%EL_IN,dg_here%IRK) = dg_here%RHS_iota(K,global_here%EL_IN,dg_here%IRK) - W_IN*global_here%I_HAT
#endif

#ifdef CHEM
               dg_here%RHS_iota(K,global_here%EL_IN,dg_here%IRK) = dg_here%RHS_iota(K,global_here%EL_IN,dg_here%IRK) - W_IN*global_here%I_HAT
               dg_here%RHS_iota2(K,global_here%EL_IN,dg_here%IRK) = dg_here%RHS_iota2(K,global_here%EL_IN,dg_here%IRK) - W_IN*global_here%J_HAT
#endif

#ifdef dynp
               dg_here%RHS_dynP(K,global_here%EL_IN,dg_here%IRK) = dg_here%RHS_dynP(K,global_here%EL_IN,dg_here%IRK) - W_IN*global_here%K_HAT
#endif

#ifdef SED_LAY
               do ll = 1,s%layers
                  dg_here%RHS_bed(K,global_here%EL_IN,dg_here%IRK,ll) = dg_here%RHS_bed(K,global_here%EL_IN,dg_here%IRK,ll) - W_IN*dg_here%bed_HAT(ll)
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
!     (This procedure dg_here%needs to be called after the wetting and 
!     drying post-process.)
!     
!     Written by Shintaro Bunya (02-13-2006)
!     
!***********************************************************************

#if 0
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
#endif
