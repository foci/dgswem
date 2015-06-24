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
      SUBROUTINE RHS_DG_HYDRO(s,dg)
      
!.....Use appropriate modules

      USE GLOBAL
      USE DG
      USE NodalAttributes, ONLY : TAU, IFLINBF, IFHYBF, HBREAK, FTHETA,&
     FGAMMA,LoadManningsN,ManningsN,CF

      USE sizes

      IMPLICIT NONE

      type (sizes_type) :: s
      type (dg_type) :: dg
      
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
         dg%advectqx(l)=0.0
         dg%advectqx(l)=0.0
         dg%sourceqx(l)=0.0
         dg%sourceqy(l)=0.0
!     nd
         
!.......Adjust the p values for constants
         
         dg%pa = PDG_EL(L)

#ifdef P0
         if (dg%pa.eq.0) then
            dg%pa = 1
         endif
#endif
         
!.......If element is dry then skip calculations
         
         IF (dg%WDFLG(L).EQ.0) then
            GOTO 1000
         endif
         
!.......Retrieve the global node numbers for the element

         N1 = NM(L,1)
         N2 = NM(L,2)
         N3 = NM(L,3)

!.......Compute avaraged values
!.......These will be used later when bottom friction is computed

         DEPTH_C = dg%HB(1,L,1) + dg%ZE(1,L,dg%IRK)
         FH_NL_C = 1.D0/(NLEQ*DEPTH_C + LEQ)
         UX_C = dg%QX(1,L,dg%IRK)*FH_NL_C
         UY_C = dg%QY(1,L,dg%IRK)*FH_NL_C
         UMAG_C = SQRT(UX_C*UX_C + UY_C*UY_C)
         
!.......Compute derivatives of Lagrange basis functions at nodes
         
         IF ((NWS.NE.0).OR.(NTIP.NE.0)) THEN
            DPSIDX(1) = dg%DRPSI(1)*dg%DRDX(L) + dg%DSPSI(1)*dg%DSDX(L)
            DPSIDX(2) = dg%DRPSI(2)*dg%DRDX(L) + dg%DSPSI(2)*dg%DSDX(L)
            DPSIDX(3) = dg%DRPSI(3)*dg%DRDX(L) + dg%DSPSI(3)*dg%DSDX(L)
            DPSIDY(1) = dg%DRPSI(1)*dg%DRDY(L) + dg%DSPSI(1)*dg%DSDY(L)
            DPSIDY(2) = dg%DRPSI(2)*dg%DRDY(L) + dg%DSPSI(2)*dg%DSDY(L)
            DPSIDY(3) = dg%DRPSI(3)*dg%DRDY(L) + dg%DSPSI(3)*dg%DSDY(L)
         ENDIF

!.......Compute dg%ZE, dg%QX, dg%QY, and dg%HB at each area Gauss quadrature point
         
         DO I = 1,dg%NAGP(dg%pa)
            
            dg%ZE_IN = dg%ZE(1,L,dg%IRK)
            dg%QX_IN = dg%QX(1,L,dg%IRK)
            dg%QY_IN = dg%QY(1,L,dg%IRK)

#ifdef TRACE
            dg%iota_IN = dg%iota(1,L,dg%IRK)
#endif

#ifdef CHEM
            dg%iota_IN = dg%iota(1,L,dg%IRK)
            dg%iota2_IN = dg%iota2(1,L,dg%IRK)
#endif

#ifdef DYNP
            dynP_IN = dg%dynP(1,L,dg%IRK)
#endif
            
            dg%HB_IN = dg%BATH(I,L,dg%pa)
            dg%DHB_X = dg%DBATHDX(I,L,dg%pa)
            dg%DHB_Y = dg%DBATHDY(I,L,dg%pa)

            !When layered, these change
#ifdef SED_LAY
            dg%HB(:,L,dg%irk) = 0.D0
            do ll = 1,s%layers
               dg%HB(1,L,dg%irk) = dg%HB(1,L,dg%irk) + dg%bed(1,L,dg%irk,ll)

               MZ_X(ll) =  dg%MZ(1,1,ll,L)
               MZ_Y(ll) =  dg%MZ(1,2,ll,L)
            enddo
            dg%HB_IN = dg%HB(1,L,dg%irk)
            dg%DHB_X = 0.D0
            dg%DHB_Y = 0.D0
            DH_Y = 0.D0
            DH_X = 0.D0
            dg%DPHIDX = 0.D0
            dg%DPHIDY = 0.D0
            dg%HB(1,L,dg%irk) = 0.D0
            do K = 1,dg%DOFS(L)
               do ll = 1,s%layers
                  dg%HB(k,L,dg%irk) = dg%HB(k,L,dg%irk) + dg%bed(k,L,dg%irk,ll)
               enddo
               dg%DPHIDX = dg%DRPHI(K,I,dg%pa)*dg%DRDX(L) + dg%DSPHI(K,I,dg%pa)*dg%DSDX(L)
               dg%DPHIDY = dg%DRPHI(K,I,dg%pa)*dg%DRDY(L) + dg%DSPHI(K,I,dg%pa)*dg%DSDY(L)
               dg%DHB_X = dg%DHB_X + dg%HB(K,L,dg%irk)*dg%DPHIDX
               dg%DHB_Y = dg%DHB_Y + dg%HB(K,L,dg%irk)*dg%DPHIDY
               DH_Y = DH_Y + (dg%HB(K,L,dg%irk)+dg%ZE(K,L,dg%irk))*dg%DPHIDY
               DH_X = DH_X + (dg%HB(K,L,dg%irk)+dg%ZE(K,L,dg%irk))*dg%DPHIDX
            enddo
#endif

#ifdef WAVE_DIF 
            HZ_X = dg%HZ(1,1,1,L)
            HZ_Y = dg%HZ(1,2,2,L)
#endif
            
            LZ_XX = dg%LZ(1,1,1,L)
            LZ_XY = dg%LZ(1,1,2,L)
            LZ_YX = dg%LZ(1,2,1,L)
            LZ_YY = dg%LZ(1,2,2,L)

#ifdef TRACE

            TZ_X = dg%TZ(1,1,1,L)
            TZ_Y = dg%TZ(1,2,2,L)

#endif
            
            SFACQUAD = dg%SFAC_ELEM(I,L,dg%pa)
            
            DO K = 2,dg%DOFS(L)
               
               dg%ZE_IN = dg%ZE_IN + dg%ZE(K,L,dg%IRK)*dg%PHI_AREA(K,I,dg%pa)
               dg%QX_IN = dg%QX_IN + dg%QX(K,L,dg%IRK)*dg%PHI_AREA(K,I,dg%pa)
               dg%QY_IN = dg%QY_IN + dg%QY(K,L,dg%IRK)*dg%PHI_AREA(K,I,dg%pa)

#ifdef WAVE_DIF 
               HZ_X = HZ_X + dg%HZ(K,1,1,L)*dg%PHI_AREA(K,I,dg%pa)
               HZ_Y = HZ_Y + dg%HZ(K,2,2,L)*dg%PHI_AREA(K,I,dg%pa)
#endif

               LZ_XX = LZ_XX + dg%LZ(K,1,1,L)*dg%PHI_AREA(K,I,dg%pa)
               LZ_XY = LZ_XY + dg%LZ(K,1,2,L)*dg%PHI_AREA(K,I,dg%pa)
               LZ_YX = LZ_YX + dg%LZ(K,2,1,L)*dg%PHI_AREA(K,I,dg%pa)
               LZ_YY = LZ_YY + dg%LZ(K,2,2,L)*dg%PHI_AREA(K,I,dg%pa)

#ifdef TRACE
               TZ_X = TZ_X + dg%TZ(K,1,1,L)*dg%PHI_AREA(K,I,dg%pa)
               TZ_Y = TZ_Y + dg%TZ(K,2,2,L)*dg%PHI_AREA(K,I,dg%pa)

               dg%iota_IN = dg%iota_IN + dg%iota(K,L,dg%IRK)*dg%PHI_AREA(K,I,dg%pa)
#endif

#ifdef CHEM
               dg%iota_IN = dg%iota_IN + dg%iota(K,L,dg%IRK)*dg%PHI_AREA(K,I,dg%pa)
               dg%iota2_IN = dg%iota2_IN + dg%iota2(K,L,dg%IRK)*dg%PHI_AREA(K,I,dg%pa)
#endif

#ifdef DYNP
               dynP_IN = dynP_IN + dg%dynP(K,L,dg%IRK)*dg%PHI_AREA(K,I,dg%pa)
#endif

               DEPTH = dg%ZE_IN + dg%HB_IN

#ifdef SED_LAY
               do ll = 1,s%layers
                  dg%bed_IN(ll) = dg%bed_IN(ll) + dg%bed(K,L,dg%IRK,ll)*dg%PHI_AREA(K,I,dg%pa)
                  dg%HB_IN = dg%HB_IN + dg%bed(K,L,dg%irk,ll)*dg%PHI_AREA(K,I,dg%pa)

                  MZ_X(ll) =  MZ_X(ll) + dg%MZ(K,1,ll,L)*dg%PHI_AREA(K,I,dg%pa)
                  MZ_Y(ll) =  MZ_Y(ll) + dg%MZ(K,2,ll,L)*dg%PHI_AREA(K,I,dg%pa)
               enddo
               DEPTH = 0.D0
               DEPTH = dg%ZE_IN + dg%HB_IN
!.........Compute sediment discharge model                                                                           

               !Note that the choice of linearization can require this to be changed
            dg%QMag_IN = (dg%QX_IN*dg%QX_IN/(DEPTH**2) + dg%QY_IN*dg%QY_IN/(DEPTH)**2)**(1/2)
            discharge_modelX_IN = dg%porosity * DEPTH**(-1) * dg%QMag_IN**(2) * dg%QX_IN*SFACQUAD                             
            discharge_modelY_IN = dg%porosity * DEPTH**(-1) * dg%QMag_IN**(2) *dg%QY_IN  

#endif
               
            ENDDO
 
!.........Compute continuity fluxes

            F1_NL = NLEQ + LEQ*dg%HB_IN

#ifdef WAVE_DIF 
            FX_IN = (dg%QX_IN+HZ_X)*F1_NL*SFACQUAD
#else
            FX_IN = dg%QX_IN*F1_NL*SFACQUAD
#endif

#ifdef WAVE_DIF 
            FY_IN = (dg%QY_IN+HZ_Y)*F1_NL
#else
            FY_IN = dg%QY_IN*F1_NL
#endif
!.........Compute momentum flux terms

            FU_NL = NLEQ*dg%QX_IN
            FV_NL = NLEQ*dg%QY_IN
            FG_NL = NLEQG*dg%ZE_IN*dg%WDFLG(L)
            FH_NL = 1.D0/(NLEQ*DEPTH + LEQ)
            U_IN  = dg%QX_IN*FH_NL
            V_IN  = dg%QY_IN*FH_NL

            HUU = FU_NL*U_IN
            HVV = FV_NL*V_IN
            HUV = FU_NL*V_IN
            GH2 = FG_NL*(0.5D0*dg%ZE_IN + dg%HB_IN) + dg%FG_L*dg%ZE_IN
#ifdef SED_LAY
            !Not well-balanced, be careful here!
            !GH2 =  0.D0
            !GH2 =  0.5D0*G*(DEPTH**2)
            !GH2 =  dg%WDFLG(L)*0.5D0*G*(DEPTH**2)
            !GH2 = FG_NL*(0.5D0*dg%ZE_IN + 2.D0*dg%HB_IN) + 0.5D0*dg%FG_L* (dg%ZE_IN**2 - DEPTH**2
#endif
           
!.........Compute x momentum fluxes

            GX_IN = (HUU + GH2 + LZ_XX)*SFACQUAD
            GY_IN = HUV + LZ_XY

            dg%advectqx(l)=dg%advectqx(l)+gx_in+gy_in

!.........Compute y momentum fluxes

            HX_IN = (HUV + LZ_YX)*SFACQUAD
            HY_IN = HVV + GH2 + LZ_YY

            dg%advectqy(l)=dg%advectqy(l)+hx_in+hy_in

!.........Compute the friction factor
            if (LoadManningsN) then
!     MN_IN=dg%MANN(1,L)
!     do k=2,dg%dof
!     MN_IN=MN_IN + dg%MANN(K,L)*dg%PHI_AREA(K,I)
!     enddo
               dg%fric_el(L)=G*&
                   ((ManningsN(n1)+ManningsN(n2)+ManningsN(n3))/3.)**2&
!     $            MN_IN**2
                   /(DEPTH**(1.d0/3.d0))
               if (dg%fric_el(L).lt.CF) dg%fric_el(L)=CF
            endif
            TAU = dg%FRIC_EL(L)
!     IF (IFLINBF.EQ.0) THEN
!     dg%UMAG = SQRT( U_IN*U_IN + V_IN*V_IN )
!     TAU  = TAU*dg%UMAG*FH_NL
!     IF (IFHYBF.EQ.1) TAU = TAU*
!     &             (1.D0  + (HBREAK*FH_NL)**FTHETA)**(FGAMMA/FTHETA)
!     ENDIF
!     Modified to compute TAU using elemental averages.
!     This seems necessary to avoid exessive bottom friction
!     at wetting-drying fronts where the total column height is very
!     small. S.B. 9-Feb-2008

            IF (IFLINBF.EQ.0) THEN
               dg%UMAG = SQRT( U_IN*U_IN + V_IN*V_IN )
!     cnd modified 4/23/10 to test friction 
!     TAU  = TAU*UMAG_C*FH_NL_C
               TAU  = TAU*dg%UMAG*FH_NL
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
!     IF (dg%RAMPDG.LT.1.D0) TAU = MAX(TAU,0.001)

!.........Compute the x momentum source/sink terms

            dg%SOURCE_X = &

!.........1.) Friction term

               - TAU*dg%QX_IN&
      
!.........2.) Bathymetric slope term

                + FG_NL*dg%DHB_X*SFACQUAD&

!.........3.) Coriolis force

                + dg%CORI_EL(L)*dg%QY_IN

!.....Compute the y momentum source/sink terms

            dg%SOURCE_Y =&

!.........1.) Friction term

                - TAU*dg%QY_IN&

!.........2) Bathymetric slope term

      
                + FG_NL*dg%DHB_Y&
     
!.........3.) Coriolis force

                - dg%CORI_EL(L)*dg%QX_IN
            
!.........4.) Wind and pressure forcing (in x and y)

            IF (NWS.NE.0) THEN
               FW_NL = 1.D0/F1_NL
               dg%SOURCE_X = dg%SOURCE_X + FW_NL*( WSX2(N1)*dg%PSI1(I,dg%pa)&
                    + WSX2(N2)*dg%PSI2(I,dg%pa)  + WSX2(N3)*dg%PSI3(I,dg%pa) )&
                    - G*SFACQUAD*DEPTH&
                    *( PR2(N1)*DPSIDX(1)&
                    + PR2(N2)*DPSIDX(2) + PR2(N3)*DPSIDX(3))
               dg%SOURCE_Y = dg%SOURCE_Y + FW_NL*( WSY2(N1)*dg%PSI1(I,dg%pa)&
                    + WSY2(N2)*dg%PSI2(I,dg%pa)  + WSY2(N3)*dg%PSI3(I,dg%pa) )&
                    - G*DEPTH*( PR2(N1)*DPSIDY(1)&
                    + PR2(N2)*DPSIDY(2) + PR2(N3)*DPSIDY(3))
            ENDIF

!     if (myproc.eq.1.and.l.eq.440.and.i.eq.1) then
!     write(440,*) 'tau ',tau,dg%umag,dg%fric_el(l)
!     endif

!.........5) Tidal potential forcing (in x and y)

            IF (NTIP.NE.0) THEN
!$$$               dg%SOURCE_X = dg%SOURCE_X + dg%RAMPDG*G*DEPTH*dg%SFAC_ELEM(I,L,dg%pa)*
!$$$     &              ( DPSIDX(1)*TIP2(N1)
!$$$     &              + DPSIDX(2)*TIP2(N2) + DPSIDX(3)*TIP2(N3) )

               dg%SOURCE_X = dg%SOURCE_X + &
                    dg%RAMPDG*G*DEPTH*SFACQUAD*( DPSIDX(1)*TIP2(N1)&
                    + DPSIDX(2)*TIP2(N2) + DPSIDX(3)*TIP2(N3) )
               
               dg%SOURCE_Y = dg%SOURCE_Y + dg%RAMPDG*G*DEPTH*( DPSIDY(1)*TIP2(N1)&
                    + DPSIDY(2)*TIP2(N2) + DPSIDY(3)*TIP2(N3) )
            ENDIF

            dg%sourceqx(l)=dg%sourceqx(l)+dg%source_x
            dg%sourceqy(l)=dg%sourceqy(l)+dg%source_y

!.........6) Chemical mass action  

#ifdef CHEM     
            MassAction1st = 0.D0
            MassAction2nd = 0.D0

            MassMax(L) = 0.D0
            
            rate = reaction_rate*1.D0/86400.D0

            MassAction1st =  rate * (max(dg%iota2_IN,0.D0)*max(dg%iota_IN,0.D0))
            MassAction2nd =  rate * (max(dg%iota2_IN,0.D0)*max(dg%iota_IN,0.D0))

            MassMax(L) = min(MassAction1st,MassAction2nd)

!.........Build the rhs

            dg%RHS_iota(1,L,dg%IRK) = dg%RHS_iota(1,L,dg%IRK) &
           + MassAction1st*dg%SRFAC(1,I,L,dg%pa)*FH_NL
            dg%RHS_iota2(1,L,dg%IRK) = dg%RHS_iota2(1,L,dg%IRK) &
           + MassAction2nd*dg%SRFAC(1,I,L,dg%pa)*FH_NL  
#endif


            dg%RHS_QX(1,L,dg%IRK) = dg%RHS_QX(1,L,dg%IRK) + dg%SRFAC(1,I,L,dg%pa)*dg%SOURCE_X
            dg%RHS_QY(1,L,dg%IRK) = dg%RHS_QY(1,L,dg%IRK) + dg%SRFAC(1,I,L,dg%pa)*dg%SOURCE_Y
            
            DO K = 2,dg%DOFS(L)

               dg%RHS_ZE(K,L,dg%IRK) = dg%RHS_ZE(K,L,dg%IRK) + dg%XFAC(K,I,L,dg%pa)*FX_IN&
              + dg%YFAC(K,I,L,dg%pa)*FY_IN 
               dg%RHS_QX(K,L,dg%IRK) = dg%RHS_QX(K,L,dg%IRK) + dg%XFAC(K,I,L,dg%pa)*GX_IN&
              + dg%YFAC(K,I,L,dg%pa)*GY_IN + dg%SRFAC(K,I,L,dg%pa)*dg%SOURCE_X
               dg%RHS_QY(K,L,dg%IRK) = dg%RHS_QY(K,L,dg%IRK) + dg%XFAC(K,I,L,dg%pa)*HX_IN&
              + dg%YFAC(K,I,L,dg%pa)*HY_IN + dg%SRFAC(K,I,L,dg%pa)*dg%SOURCE_Y

#ifdef SED_LAY

               do ll = 1,s%layers !only really makes sense for single layer

                  dg%RHS_bed(K,L,dg%IRK,ll) = dg%RHS_bed(K,L,dg%IRK,ll) &
                 + dg%XFAC(K,I,L,dg%pa)*(discharge_modelX_IN+MZ_X(ll)*SFACQUAD)&
                 + dg%YFAC(K,I,L,dg%pa)*(discharge_modelY_IN+MZ_Y(ll))

               enddo
#endif

#ifdef TRACE
               dg%RHS_iota(K,L,dg%IRK) = dg%RHS_iota(K,L,dg%IRK) &
              + dg%XFAC(K,I,L,dg%pa)*(dg%iota_IN*dg%QX_IN*FH_NL+TZ_X*SFACQUAD) &
              + dg%YFAC(K,I,L,dg%pa)*(dg%iota_IN*dg%QY_IN*FH_NL+TZ_Y) 
#endif

#ifdef CHEM
               dg%RHS_iota(K,L,dg%IRK) = dg%RHS_iota(K,L,dg%IRK) &
              + dg%XFAC(K,I,L,dg%pa)*dg%iota_IN*dg%QX_IN*FH_NL &
              + dg%YFAC(K,I,L,dg%pa)*dg%iota_IN*dg%QY_IN*FH_NL &
              + MassAction1st*dg%SRFAC(K,I,L,dg%pa)*FH_NL 
               
               dg%RHS_iota2(K,L,dg%IRK) = dg%RHS_iota2(K,L,dg%IRK) &
              + dg%XFAC(K,I,L,dg%pa)*dg%iota2_IN*dg%QX_IN*FH_NL      &
              + dg%YFAC(K,I,L,dg%pa)*dg%iota2_IN*dg%QY_IN*FH_NL &
              + MassAction2nd*dg%SRFAC(K,I,L,dg%pa)*FH_NL 
#endif

#ifdef DYNP
               dg%RHS_dynP(K,L,dg%IRK) = dg%RHS_dynP(K,L,dg%IRK) &
              + dg%XFAC(K,I,L,dg%pa)*dynP_IN*dg%QX_IN*FH_NL &
              + dg%YFAC(K,I,L,dg%pa)*dynP_IN*dg%QY_IN*FH_NL&
              + dg%subphi_IN*dg%SRFAC(K,I,L,dg%pa)*FH_NL 
#endif
               
            ENDDO
            
         ENDDO  

         
 1000 CONTINUE      
      RETURN
      END SUBROUTINE
