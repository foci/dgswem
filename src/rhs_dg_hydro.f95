!***********************************************************************
!     
!     SUBROUTINE RHS_DG_HYDRO()
!     
!     This subroutine computes the area integrals for the DG hydro and
!     adds them into the RHS.
!     
!     Written by Ethan Kubatko (06-11-2004)
!     
!-----------------------------------------------------------------------
!     
!     Feb 23, 2007, sb, Modified for better performance
!     Jan 02, 2007, sb, Modified for LDG
!     Aug xx, 2005, sb, Modified for wetting/drying
!     01-10-2011 - cem - adapted for p_enrichment and multicomponent
!     
!***********************************************************************
      SUBROUTINE RHS_DG_HYDRO(s,dg_here)
      
!.....Use appropriate modules

      USE GLOBAL
      USE DG
      USE NodalAttributes, ONLY : TAU, IFLINBF, IFHYBF, HBREAK, FTHETA,&
     FGAMMA,LoadManningsN,ManningsN,CF

      USE sizes

      IMPLICIT NONE

      type (sizes_type) :: s
      type (dg_type) :: dg_here
      
!.....Declare local variables

      INTEGER L,k,i,ll
      REAL(SZ) DPSIDX(3), DPSIDY(3)
      REAL(SZ) AREA, IMASS, TKX, TKY, Xpart, Ypart,H_0,C_1
      REAL(SZ) PHI_AREA_KI,MN_IN, MassAction1st,MassAction2nd,fx,fy
      REAL(SZ) LZ_XX, LZ_XY, LZ_YX, LZ_YY, rate, s_mass, s_sed,b_0
      REAL(SZ) DEPTH, F1_NL, FU_NL, FV_NL, FG_NL, FH_NL, FW_NL
      REAL(SZ) HUU, HVV, HUV, GH2,MZ_X(s%layers),MZ_Y(s%layers), fgauss, sig
      REAL(SZ) DEPTH_C, FH_NL_C, UX_C, UY_C, UMAG_C, DTDPH,SFACQUAD
      Real(SZ) discharge_modelX_IN,discharge_modelY_IN
      Real(SZ) DH_X,DH_Y,phi_tot,C_0,HZ_X,HZ_Y,TZ_X,TZ_Y

      Real(SZ),allocatable :: XBCbt(:),YBCbt(:)

      Allocate ( XBCbt(s%MNE),YBCbt(s%MNE) )

      DTDPH = 1.D0/DTDP
      DO 1000 L = 1, NE
!     nd
         dg_here%advectqx(l)=0.0
         dg_here%advectqx(l)=0.0
         dg_here%sourceqx(l)=0.0
         dg_here%sourceqy(l)=0.0
!     nd
         
!.......Adjust the p values for constants
         
         dg_here%pa = PDG_EL(L)

#ifdef P0
         if (dg_here%pa.eq.0) then
            dg_here%pa = 1
         endif
#endif
         
!.......If element is dry then skip calculations
         
         IF (dg_here%WDFLG(L).EQ.0) then
            GOTO 1000
         endif
         
!.......Retrieve the global node numbers for the element

         N1 = NM(L,1)
         N2 = NM(L,2)
         N3 = NM(L,3)

!.......Compute avaraged values
!.......These will be used later when bottom friction is computed

         DEPTH_C = dg_here%HB(1,L,1) + dg_here%ZE(1,L,dg_here%IRK)
         FH_NL_C = 1.D0/(NLEQ*DEPTH_C + LEQ)
         UX_C = dg_here%QX(1,L,dg_here%IRK)*FH_NL_C
         UY_C = dg_here%QY(1,L,dg_here%IRK)*FH_NL_C
         UMAG_C = SQRT(UX_C*UX_C + UY_C*UY_C)
         
!.......Compute derivatives of Lagrange basis functions at nodes
         
         IF ((NWS.NE.0).OR.(NTIP.NE.0)) THEN
            DPSIDX(1) = dg_here%DRPSI(1)*dg_here%DRDX(L) + dg_here%DSPSI(1)*dg_here%DSDX(L)
            DPSIDX(2) = dg_here%DRPSI(2)*dg_here%DRDX(L) + dg_here%DSPSI(2)*dg_here%DSDX(L)
            DPSIDX(3) = dg_here%DRPSI(3)*dg_here%DRDX(L) + dg_here%DSPSI(3)*dg_here%DSDX(L)
            DPSIDY(1) = dg_here%DRPSI(1)*dg_here%DRDY(L) + dg_here%DSPSI(1)*dg_here%DSDY(L)
            DPSIDY(2) = dg_here%DRPSI(2)*dg_here%DRDY(L) + dg_here%DSPSI(2)*dg_here%DSDY(L)
            DPSIDY(3) = dg_here%DRPSI(3)*dg_here%DRDY(L) + dg_here%DSPSI(3)*dg_here%DSDY(L)
         ENDIF

!.......Compute dg_here%ZE, dg_here%QX, dg_here%QY, and dg_here%HB at each area Gauss quadrature point
         
         DO I = 1,dg_here%NAGP(dg_here%pa)
            
            dg_here%ZE_IN = dg_here%ZE(1,L,dg_here%IRK)
            dg_here%QX_IN = dg_here%QX(1,L,dg_here%IRK)
            dg_here%QY_IN = dg_here%QY(1,L,dg_here%IRK)

#ifdef TRACE
            dg_here%iota_IN = dg_here%iota(1,L,dg_here%IRK)
#endif

#ifdef CHEM
            dg_here%iota_IN = dg_here%iota(1,L,dg_here%IRK)
            dg_here%iota2_IN = dg_here%iota2(1,L,dg_here%IRK)
#endif

#ifdef DYNP
            dynP_IN = dg_here%dynP(1,L,dg_here%IRK)
#endif
            
            dg_here%HB_IN = dg_here%BATH(I,L,dg_here%pa)
            dg_here%DHB_X = dg_here%DBATHDX(I,L,dg_here%pa)
            dg_here%DHB_Y = dg_here%DBATHDY(I,L,dg_here%pa)

            !When layered, these change
#ifdef SED_LAY
            dg_here%HB(:,L,dg_here%irk) = 0.D0
            do ll = 1,s%layers
               dg_here%HB(1,L,dg_here%irk) = dg_here%HB(1,L,dg_here%irk) + dg_here%bed(1,L,dg_here%irk,ll)

               MZ_X(ll) =  dg_here%MZ(1,1,ll,L)
               MZ_Y(ll) =  dg_here%MZ(1,2,ll,L)
            enddo
            dg_here%HB_IN = dg_here%HB(1,L,dg_here%irk)
            dg_here%DHB_X = 0.D0
            dg_here%DHB_Y = 0.D0
            DH_Y = 0.D0
            DH_X = 0.D0
            dg_here%DPHIDX = 0.D0
            dg_here%DPHIDY = 0.D0
            dg_here%HB(1,L,dg_here%irk) = 0.D0
            do K = 1,dg_here%DOFS(L)
               do ll = 1,s%layers
                  dg_here%HB(k,L,dg_here%irk) = dg_here%HB(k,L,dg_here%irk) + dg_here%bed(k,L,dg_here%irk,ll)
               enddo
               dg_here%DPHIDX = dg_here%DRPHI(K,I,dg_here%pa)*dg_here%DRDX(L) + dg_here%DSPHI(K,I,dg_here%pa)*dg_here%DSDX(L)
               dg_here%DPHIDY = dg_here%DRPHI(K,I,dg_here%pa)*dg_here%DRDY(L) + dg_here%DSPHI(K,I,dg_here%pa)*dg_here%DSDY(L)
               dg_here%DHB_X = dg_here%DHB_X + dg_here%HB(K,L,dg_here%irk)*dg_here%DPHIDX
               dg_here%DHB_Y = dg_here%DHB_Y + dg_here%HB(K,L,dg_here%irk)*dg_here%DPHIDY
               DH_Y = DH_Y + (dg_here%HB(K,L,dg_here%irk)+dg_here%ZE(K,L,dg_here%irk))*dg_here%DPHIDY
               DH_X = DH_X + (dg_here%HB(K,L,dg_here%irk)+dg_here%ZE(K,L,dg_here%irk))*dg_here%DPHIDX
            enddo
#endif

#ifdef WAVE_DIF 
            HZ_X = dg_here%HZ(1,1,1,L)
            HZ_Y = dg_here%HZ(1,2,2,L)
#endif
            
            LZ_XX = dg_here%LZ(1,1,1,L)
            LZ_XY = dg_here%LZ(1,1,2,L)
            LZ_YX = dg_here%LZ(1,2,1,L)
            LZ_YY = dg_here%LZ(1,2,2,L)

#ifdef TRACE

            TZ_X = dg_here%TZ(1,1,1,L)
            TZ_Y = dg_here%TZ(1,2,2,L)

#endif
            
            SFACQUAD = dg_here%SFAC_ELEM(I,L,dg_here%pa)
            
            DO K = 2,dg_here%DOFS(L)
               
               dg_here%ZE_IN = dg_here%ZE_IN + dg_here%ZE(K,L,dg_here%IRK)*dg_here%PHI_AREA(K,I,dg_here%pa)
               dg_here%QX_IN = dg_here%QX_IN + dg_here%QX(K,L,dg_here%IRK)*dg_here%PHI_AREA(K,I,dg_here%pa)
               dg_here%QY_IN = dg_here%QY_IN + dg_here%QY(K,L,dg_here%IRK)*dg_here%PHI_AREA(K,I,dg_here%pa)

#ifdef WAVE_DIF 
               HZ_X = HZ_X + dg_here%HZ(K,1,1,L)*dg_here%PHI_AREA(K,I,dg_here%pa)
               HZ_Y = HZ_Y + dg_here%HZ(K,2,2,L)*dg_here%PHI_AREA(K,I,dg_here%pa)
#endif

               LZ_XX = LZ_XX + dg_here%LZ(K,1,1,L)*dg_here%PHI_AREA(K,I,dg_here%pa)
               LZ_XY = LZ_XY + dg_here%LZ(K,1,2,L)*dg_here%PHI_AREA(K,I,dg_here%pa)
               LZ_YX = LZ_YX + dg_here%LZ(K,2,1,L)*dg_here%PHI_AREA(K,I,dg_here%pa)
               LZ_YY = LZ_YY + dg_here%LZ(K,2,2,L)*dg_here%PHI_AREA(K,I,dg_here%pa)

#ifdef TRACE
               TZ_X = TZ_X + dg_here%TZ(K,1,1,L)*dg_here%PHI_AREA(K,I,dg_here%pa)
               TZ_Y = TZ_Y + dg_here%TZ(K,2,2,L)*dg_here%PHI_AREA(K,I,dg_here%pa)

               dg_here%iota_IN = dg_here%iota_IN + dg_here%iota(K,L,dg_here%IRK)*dg_here%PHI_AREA(K,I,dg_here%pa)
#endif

#ifdef CHEM
               dg_here%iota_IN = dg_here%iota_IN + dg_here%iota(K,L,dg_here%IRK)*dg_here%PHI_AREA(K,I,dg_here%pa)
               dg_here%iota2_IN = dg_here%iota2_IN + dg_here%iota2(K,L,dg_here%IRK)*dg_here%PHI_AREA(K,I,dg_here%pa)
#endif

#ifdef DYNP
               dynP_IN = dynP_IN + dg_here%dynP(K,L,dg_here%IRK)*dg_here%PHI_AREA(K,I,dg_here%pa)
#endif

               DEPTH = dg_here%ZE_IN + dg_here%HB_IN

#ifdef SED_LAY
               do ll = 1,s%layers
                  dg_here%bed_IN(ll) = dg_here%bed_IN(ll) + dg_here%bed(K,L,dg_here%IRK,ll)*dg_here%PHI_AREA(K,I,dg_here%pa)
                  dg_here%HB_IN = dg_here%HB_IN + dg_here%bed(K,L,dg_here%irk,ll)*dg_here%PHI_AREA(K,I,dg_here%pa)

                  MZ_X(ll) =  MZ_X(ll) + dg_here%MZ(K,1,ll,L)*dg_here%PHI_AREA(K,I,dg_here%pa)
                  MZ_Y(ll) =  MZ_Y(ll) + dg_here%MZ(K,2,ll,L)*dg_here%PHI_AREA(K,I,dg_here%pa)
               enddo
               DEPTH = 0.D0
               DEPTH = dg_here%ZE_IN + dg_here%HB_IN
!.........Compute sediment discharge model                                                                           

               !Note that the choice of linearization can require this to be changed
            dg_here%QMag_IN = (dg_here%QX_IN*dg_here%QX_IN/(DEPTH**2) + dg_here%QY_IN*dg_here%QY_IN/(DEPTH)**2)**(1/2)
            discharge_modelX_IN = dg_here%porosity * DEPTH**(-1) * dg_here%QMag_IN**(2) * dg_here%QX_IN*SFACQUAD                             
            discharge_modelY_IN = dg_here%porosity * DEPTH**(-1) * dg_here%QMag_IN**(2) *dg_here%QY_IN  

#endif
               
            ENDDO
 
!.........Compute continuity fluxes

            F1_NL = NLEQ + LEQ*dg_here%HB_IN

#ifdef WAVE_DIF 
            FX_IN = (dg_here%QX_IN+HZ_X)*F1_NL*SFACQUAD
#else
            FX_IN = dg_here%QX_IN*F1_NL*SFACQUAD
#endif

#ifdef WAVE_DIF 
            FY_IN = (dg_here%QY_IN+HZ_Y)*F1_NL
#else
            FY_IN = dg_here%QY_IN*F1_NL
#endif
!.........Compute momentum flux terms

            FU_NL = NLEQ*dg_here%QX_IN
            FV_NL = NLEQ*dg_here%QY_IN
            FG_NL = NLEQG*dg_here%ZE_IN*dg_here%WDFLG(L)
            FH_NL = 1.D0/(NLEQ*DEPTH + LEQ)
            U_IN  = dg_here%QX_IN*FH_NL
            V_IN  = dg_here%QY_IN*FH_NL

            HUU = FU_NL*U_IN
            HVV = FV_NL*V_IN
            HUV = FU_NL*V_IN
            GH2 = FG_NL*(0.5D0*dg_here%ZE_IN + dg_here%HB_IN) + dg_here%FG_L*dg_here%ZE_IN
#ifdef SED_LAY
            !Not well-balanced, be careful here!
            !GH2 =  0.D0
            !GH2 =  0.5D0*G*(DEPTH**2)
            !GH2 =  dg_here%WDFLG(L)*0.5D0*G*(DEPTH**2)
            !GH2 = FG_NL*(0.5D0*dg_here%ZE_IN + 2.D0*dg_here%HB_IN) + 0.5D0*dg_here%FG_L* (dg_here%ZE_IN**2 - DEPTH**2
#endif
           
!.........Compute x momentum fluxes

            GX_IN = (HUU + GH2 + LZ_XX)*SFACQUAD
            GY_IN = HUV + LZ_XY

            dg_here%advectqx(l)=dg_here%advectqx(l)+gx_in+gy_in

!.........Compute y momentum fluxes

            HX_IN = (HUV + LZ_YX)*SFACQUAD
            HY_IN = HVV + GH2 + LZ_YY

            dg_here%advectqy(l)=dg_here%advectqy(l)+hx_in+hy_in

!.........Compute the friction factor
            if (LoadManningsN) then
!     MN_IN=dg_here%MANN(1,L)
!     do k=2,dg_here%dof
!     MN_IN=MN_IN + dg_here%MANN(K,L)*dg_here%PHI_AREA(K,I)
!     enddo
               dg_here%fric_el(L)=G*&
                   ((ManningsN(n1)+ManningsN(n2)+ManningsN(n3))/3.)**2&
!     $            MN_IN**2
                   /(DEPTH**(1.d0/3.d0))
               if (dg_here%fric_el(L).lt.CF) dg_here%fric_el(L)=CF
            endif
            TAU = dg_here%FRIC_EL(L)
!     IF (IFLINBF.EQ.0) THEN
!     dg_here%UMAG = SQRT( U_IN*U_IN + V_IN*V_IN )
!     TAU  = TAU*dg_here%UMAG*FH_NL
!     IF (IFHYBF.EQ.1) TAU = TAU*
!     &             (1.D0  + (HBREAK*FH_NL)**FTHETA)**(FGAMMA/FTHETA)
!     ENDIF
!     Modified to compute TAU using elemental averages.
!     This seems necessary to avoid exessive bottom friction
!     at wetting-drying fronts where the total column height is very
!     small. S.B. 9-Feb-2008

            IF (IFLINBF.EQ.0) THEN
               dg_here%UMAG = SQRT( U_IN*U_IN + V_IN*V_IN )
!     cnd modified 4/23/10 to test friction 
!     TAU  = TAU*UMAG_C*FH_NL_C
               TAU  = TAU*dg_here%UMAG*FH_NL
               IF (IFHYBF.EQ.1) TAU = TAU*&
              (1.D0  + (HBREAK*FH_NL_C)**FTHETA)**(FGAMMA/FTHETA)
!     &            (1.D0  + (HBREAK*FH_NL)**FTHETA)**(FGAMMA/FTHETA)
!     It is numerically probable that the bottom friction is large enoght
!     to reverse the direction of currents backward due to a too small column
!     height even though it does not happen in reality. To avoid this, the MIN
!     function bellow is added. It is expected that this MIN function upper-limits
!     TAU so the bottom friction force does not reverse the currents within
!     half a time step.  S.B. 9-Feb-2008
!     IF(TAU.GT.2.D0*DTDPH) PRINT *, "TAU = ", TAU, 2.D0*DTDPH,
!     $           u_in,v_in,myproc,L
!     TAU = MIN(TAU, 2.D0*DTDPH)
               TAU = MIN(TAU, .9D0*DTDPH)
            ENDIF
!     IF (dg_here%RAMPDG.LT.1.D0) TAU = MAX(TAU,0.001)

!.........Compute the x momentum source/sink terms

            dg_here%SOURCE_X = &

!.........1.) Friction term

               - TAU*dg_here%QX_IN&
      
!.........2.) Bathymetric slope term

                + FG_NL*dg_here%DHB_X*SFACQUAD&

!.........3.) Coriolis force

                + dg_here%CORI_EL(L)*dg_here%QY_IN

!.....Compute the y momentum source/sink terms

            dg_here%SOURCE_Y =&

!.........1.) Friction term

                - TAU*dg_here%QY_IN&

!.........2) Bathymetric slope term

      
                + FG_NL*dg_here%DHB_Y&
     
!.........3.) Coriolis force

                - dg_here%CORI_EL(L)*dg_here%QX_IN
            
!.........4.) Wind and pressure forcing (in x and y)

            IF (NWS.NE.0) THEN
               FW_NL = 1.D0/F1_NL
               dg_here%SOURCE_X = dg_here%SOURCE_X + FW_NL*( WSX2(N1)*dg_here%PSI1(I,dg_here%pa)&
                    + WSX2(N2)*dg_here%PSI2(I,dg_here%pa)  + WSX2(N3)*dg_here%PSI3(I,dg_here%pa) )&
                    - G*SFACQUAD*DEPTH&
                    *( PR2(N1)*DPSIDX(1)&
                    + PR2(N2)*DPSIDX(2) + PR2(N3)*DPSIDX(3))
               dg_here%SOURCE_Y = dg_here%SOURCE_Y + FW_NL*( WSY2(N1)*dg_here%PSI1(I,dg_here%pa)&
                    + WSY2(N2)*dg_here%PSI2(I,dg_here%pa)  + WSY2(N3)*dg_here%PSI3(I,dg_here%pa) )&
                    - G*DEPTH*( PR2(N1)*DPSIDY(1)&
                    + PR2(N2)*DPSIDY(2) + PR2(N3)*DPSIDY(3))
            ENDIF

!     if (myproc.eq.1.and.l.eq.440.and.i.eq.1) then
!     write(440,*) 'tau ',tau,dg_here%umag,dg_here%fric_el(l)
!     endif

!.........5) Tidal potential forcing (in x and y)

            IF (NTIP.NE.0) THEN
!$$$               dg_here%SOURCE_X = dg_here%SOURCE_X + dg_here%RAMPDG*G*DEPTH*dg_here%SFAC_ELEM(I,L,dg_here%pa)*
!$$$     &              ( DPSIDX(1)*TIP2(N1)
!$$$     &              + DPSIDX(2)*TIP2(N2) + DPSIDX(3)*TIP2(N3) )

               dg_here%SOURCE_X = dg_here%SOURCE_X + &
                    dg_here%RAMPDG*G*DEPTH*SFACQUAD*( DPSIDX(1)*TIP2(N1)&
                    + DPSIDX(2)*TIP2(N2) + DPSIDX(3)*TIP2(N3) )
               
               dg_here%SOURCE_Y = dg_here%SOURCE_Y + dg_here%RAMPDG*G*DEPTH*( DPSIDY(1)*TIP2(N1)&
                    + DPSIDY(2)*TIP2(N2) + DPSIDY(3)*TIP2(N3) )
            ENDIF

            dg_here%sourceqx(l)=dg_here%sourceqx(l)+dg_here%source_x
            dg_here%sourceqy(l)=dg_here%sourceqy(l)+dg_here%source_y

!.........6) Chemical mass action  

#ifdef CHEM     
            MassAction1st = 0.D0
            MassAction2nd = 0.D0

            MassMax(L) = 0.D0
            
            rate = reaction_rate*1.D0/86400.D0

            MassAction1st =  rate * (max(dg_here%iota2_IN,0.D0)*max(dg_here%iota_IN,0.D0))
            MassAction2nd =  rate * (max(dg_here%iota2_IN,0.D0)*max(dg_here%iota_IN,0.D0))

            MassMax(L) = min(MassAction1st,MassAction2nd)

!.........Build the rhs

            dg_here%RHS_iota(1,L,dg_here%IRK) = dg_here%RHS_iota(1,L,dg_here%IRK) &
           + MassAction1st*dg_here%SRFAC(1,I,L,dg_here%pa)*FH_NL
            dg_here%RHS_iota2(1,L,dg_here%IRK) = dg_here%RHS_iota2(1,L,dg_here%IRK) &
           + MassAction2nd*dg_here%SRFAC(1,I,L,dg_here%pa)*FH_NL  
#endif


            dg_here%RHS_QX(1,L,dg_here%IRK) = dg_here%RHS_QX(1,L,dg_here%IRK) + dg_here%SRFAC(1,I,L,dg_here%pa)*dg_here%SOURCE_X
            dg_here%RHS_QY(1,L,dg_here%IRK) = dg_here%RHS_QY(1,L,dg_here%IRK) + dg_here%SRFAC(1,I,L,dg_here%pa)*dg_here%SOURCE_Y
            
            DO K = 2,dg_here%DOFS(L)

               dg_here%RHS_ZE(K,L,dg_here%IRK) = dg_here%RHS_ZE(K,L,dg_here%IRK) + dg_here%XFAC(K,I,L,dg_here%pa)*FX_IN&
              + dg_here%YFAC(K,I,L,dg_here%pa)*FY_IN 
               dg_here%RHS_QX(K,L,dg_here%IRK) = dg_here%RHS_QX(K,L,dg_here%IRK) + dg_here%XFAC(K,I,L,dg_here%pa)*GX_IN&
              + dg_here%YFAC(K,I,L,dg_here%pa)*GY_IN + dg_here%SRFAC(K,I,L,dg_here%pa)*dg_here%SOURCE_X
               dg_here%RHS_QY(K,L,dg_here%IRK) = dg_here%RHS_QY(K,L,dg_here%IRK) + dg_here%XFAC(K,I,L,dg_here%pa)*HX_IN&
              + dg_here%YFAC(K,I,L,dg_here%pa)*HY_IN + dg_here%SRFAC(K,I,L,dg_here%pa)*dg_here%SOURCE_Y

#ifdef SED_LAY

               do ll = 1,s%layers !only really makes sense for single layer

                  dg_here%RHS_bed(K,L,dg_here%IRK,ll) = dg_here%RHS_bed(K,L,dg_here%IRK,ll) &
                 + dg_here%XFAC(K,I,L,dg_here%pa)*(discharge_modelX_IN+MZ_X(ll)*SFACQUAD)&
                 + dg_here%YFAC(K,I,L,dg_here%pa)*(discharge_modelY_IN+MZ_Y(ll))

               enddo
#endif

#ifdef TRACE
               dg_here%RHS_iota(K,L,dg_here%IRK) = dg_here%RHS_iota(K,L,dg_here%IRK) &
              + dg_here%XFAC(K,I,L,dg_here%pa)*(dg_here%iota_IN*dg_here%QX_IN*FH_NL+TZ_X*SFACQUAD) &
              + dg_here%YFAC(K,I,L,dg_here%pa)*(dg_here%iota_IN*dg_here%QY_IN*FH_NL+TZ_Y) 
#endif

#ifdef CHEM
               dg_here%RHS_iota(K,L,dg_here%IRK) = dg_here%RHS_iota(K,L,dg_here%IRK) &
              + dg_here%XFAC(K,I,L,dg_here%pa)*dg_here%iota_IN*dg_here%QX_IN*FH_NL &
              + dg_here%YFAC(K,I,L,dg_here%pa)*dg_here%iota_IN*dg_here%QY_IN*FH_NL &
              + MassAction1st*dg_here%SRFAC(K,I,L,dg_here%pa)*FH_NL 
               
               dg_here%RHS_iota2(K,L,dg_here%IRK) = dg_here%RHS_iota2(K,L,dg_here%IRK) &
              + dg_here%XFAC(K,I,L,dg_here%pa)*dg_here%iota2_IN*dg_here%QX_IN*FH_NL      &
              + dg_here%YFAC(K,I,L,dg_here%pa)*dg_here%iota2_IN*dg_here%QY_IN*FH_NL &
              + MassAction2nd*dg_here%SRFAC(K,I,L,dg_here%pa)*FH_NL 
#endif

#ifdef DYNP
               dg_here%RHS_dynP(K,L,dg_here%IRK) = dg_here%RHS_dynP(K,L,dg_here%IRK) &
              + dg_here%XFAC(K,I,L,dg_here%pa)*dynP_IN*dg_here%QX_IN*FH_NL &
              + dg_here%YFAC(K,I,L,dg_here%pa)*dynP_IN*dg_here%QY_IN*FH_NL&
              + dg_here%subphi_IN*dg_here%SRFAC(K,I,L,dg_here%pa)*FH_NL 
#endif
               
            ENDDO
            
         ENDDO  

         
 1000 CONTINUE      
      RETURN
      END SUBROUTINE
