C***********************************************************************
C     
C     SUBROUTINE RHS_DG_HYDRO()
C     
C     This subroutine computes the area integrals for the DG hydro and
C     adds them into the RHS.
C     
C     Written by Ethan Kubatko (06-11-2004)
C     
C-----------------------------------------------------------------------
C     
C     Feb 23, 2007, sb, Modified for better performance
C     Jan 02, 2007, sb, Modified for LDG
C     Aug xx, 2005, sb, Modified for wetting/drying
C     01-10-2011 - cem - adapted for p_enrichment and multicomponent
C     
C***********************************************************************
      SUBROUTINE RHS_DG_HYDRO()
      
C.....Use appropriate modules

      USE GLOBAL
      USE DG
      USE NodalAttributes, ONLY : TAU, IFLINBF, IFHYBF, HBREAK, FTHETA,
     &     FGAMMA,LoadManningsN,ManningsN,CF

      USE sizes, ONLY: myproc,layers

      IMPLICIT NONE
      
C.....Declare local variables

      INTEGER L,k,i,ll
      REAL(SZ) DPSIDX(3), DPSIDY(3)
      REAL(SZ) AREA, IMASS, TKX, TKY, Xpart, Ypart,H_0,C_1
      REAL(SZ) PHI_AREA_KI,MN_IN, MassAction1st,MassAction2nd,fx,fy
      REAL(SZ) LZ_XX, LZ_XY, LZ_YX, LZ_YY, rate, s_mass, s_sed,b_0
      REAL(SZ) DEPTH, F1_NL, FU_NL, FV_NL, FG_NL, FH_NL, FW_NL
      REAL(SZ) HUU, HVV, HUV, GH2,MZ_X(layers),MZ_Y(layers), fgauss, sig
      REAL(SZ) DEPTH_C, FH_NL_C, UX_C, UY_C, UMAG_C, DTDPH,SFACQUAD
      Real(SZ) discharge_modelX_IN,discharge_modelY_IN
      Real(SZ) DH_X,DH_Y,phi_tot,C_0,HZ_X,HZ_Y,TZ_X,TZ_Y

      Real(SZ),allocatable :: XBCbt(:),YBCbt(:)

      Allocate ( XBCbt(MNE),YBCbt(MNE) )

      DTDPH = 1.D0/DTDP
      DO 1000 L = 1, NE
c     nd
         advectqx(l)=0.0
         advectqx(l)=0.0
         sourceqx(l)=0.0
         sourceqy(l)=0.0
c     nd
         
C.......Adjust the p values for constants
         
         pa = PDG_EL(L)

#ifdef P0
         if (pa.eq.0) then
            pa = 1
         endif
#endif
         
C.......If element is dry then skip calculations
         
         IF (WDFLG(L).EQ.0) then
            GOTO 1000
         endif
         
C.......Retrieve the global node numbers for the element

         N1 = NM(L,1)
         N2 = NM(L,2)
         N3 = NM(L,3)

C.......Compute avaraged values
C.......These will be used later when bottom friction is computed

         DEPTH_C = HB(1,L,1) + ZE(1,L,IRK)
         FH_NL_C = 1.D0/(NLEQ*DEPTH_C + LEQ)
         UX_C = QX(1,L,IRK)*FH_NL_C
         UY_C = QY(1,L,IRK)*FH_NL_C
         UMAG_C = SQRT(UX_C*UX_C + UY_C*UY_C)
         
C.......Compute derivatives of Lagrange basis functions at nodes
         
         IF ((NWS.NE.0).OR.(NTIP.NE.0)) THEN
            DPSIDX(1) = DRPSI(1)*DRDX(L) + DSPSI(1)*DSDX(L)
            DPSIDX(2) = DRPSI(2)*DRDX(L) + DSPSI(2)*DSDX(L)
            DPSIDX(3) = DRPSI(3)*DRDX(L) + DSPSI(3)*DSDX(L)
            DPSIDY(1) = DRPSI(1)*DRDY(L) + DSPSI(1)*DSDY(L)
            DPSIDY(2) = DRPSI(2)*DRDY(L) + DSPSI(2)*DSDY(L)
            DPSIDY(3) = DRPSI(3)*DRDY(L) + DSPSI(3)*DSDY(L)
         ENDIF

C.......Compute ZE, QX, QY, and HB at each area Gauss quadrature point
         
         DO I = 1,NAGP(pa)
            
            ZE_IN = ZE(1,L,IRK)
            QX_IN = QX(1,L,IRK)
            QY_IN = QY(1,L,IRK)

#ifdef TRACE
            iota_IN = iota(1,L,IRK)
#endif

#ifdef CHEM
            iota_IN = iota(1,L,IRK)
            iota2_IN = iota2(1,L,IRK)
#endif

#ifdef DYNP
            dynP_IN = dynP(1,L,IRK)
#endif
            
            HB_IN = BATH(I,L,pa)
            DHB_X = DBATHDX(I,L,pa)
            DHB_Y = DBATHDY(I,L,pa)

#ifdef SED_LAY                  !When layered, these change
            HB(:,L,irk) = 0.D0
            do ll = 1,layers
               HB(1,L,irk) = HB(1,L,irk) + bed(1,L,irk,ll)

               MZ_X(ll) =  MZ(1,1,ll,L)
               MZ_Y(ll) =  MZ(1,2,ll,L)
            enddo
            HB_IN = HB(1,L,irk)
            DHB_X = 0.D0
            DHB_Y = 0.D0
            DH_Y = 0.D0
            DH_X = 0.D0
            DPHIDX = 0.D0
            DPHIDY = 0.D0
            HB(1,L,irk) = 0.D0
            do K = 1,DOFS(L)
               do ll = 1,layers
                  HB(k,L,irk) = HB(k,L,irk) + bed(k,L,irk,ll)
               enddo
               DPHIDX = DRPHI(K,I,pa)*DRDX(L) + DSPHI(K,I,pa)*DSDX(L)
               DPHIDY = DRPHI(K,I,pa)*DRDY(L) + DSPHI(K,I,pa)*DSDY(L)
               DHB_X = DHB_X + HB(K,L,irk)*DPHIDX
               DHB_Y = DHB_Y + HB(K,L,irk)*DPHIDY
               DH_Y = DH_Y + (HB(K,L,irk)+ZE(K,L,irk))*DPHIDY
               DH_X = DH_X + (HB(K,L,irk)+ZE(K,L,irk))*DPHIDX
            enddo
#endif

#ifdef WAVE_DIF 
            HZ_X = HZ(1,1,1,L)
            HZ_Y = HZ(1,2,2,L)
#endif
            
            LZ_XX = LZ(1,1,1,L)
            LZ_XY = LZ(1,1,2,L)
            LZ_YX = LZ(1,2,1,L)
            LZ_YY = LZ(1,2,2,L)

#ifdef TRACE

            TZ_X = TZ(1,1,1,L)
            TZ_Y = TZ(1,2,2,L)

#endif
            
            SFACQUAD = SFAC_ELEM(I,L,pa)
            
            DO K = 2,DOFS(L)
               
               ZE_IN = ZE_IN + ZE(K,L,IRK)*PHI_AREA(K,I,pa)
               QX_IN = QX_IN + QX(K,L,IRK)*PHI_AREA(K,I,pa)
               QY_IN = QY_IN + QY(K,L,IRK)*PHI_AREA(K,I,pa)

#ifdef WAVE_DIF 
               HZ_X = HZ_X + HZ(K,1,1,L)*PHI_AREA(K,I,pa)
               HZ_Y = HZ_Y + HZ(K,2,2,L)*PHI_AREA(K,I,pa)
#endif

               LZ_XX = LZ_XX + LZ(K,1,1,L)*PHI_AREA(K,I,pa)
               LZ_XY = LZ_XY + LZ(K,1,2,L)*PHI_AREA(K,I,pa)
               LZ_YX = LZ_YX + LZ(K,2,1,L)*PHI_AREA(K,I,pa)
               LZ_YY = LZ_YY + LZ(K,2,2,L)*PHI_AREA(K,I,pa)

#ifdef TRACE
               TZ_X = TZ_X + TZ(K,1,1,L)*PHI_AREA(K,I,pa)
               TZ_Y = TZ_Y + TZ(K,2,2,L)*PHI_AREA(K,I,pa)

               iota_IN = iota_IN + iota(K,L,IRK)*PHI_AREA(K,I,pa)
#endif

#ifdef CHEM
               iota_IN = iota_IN + iota(K,L,IRK)*PHI_AREA(K,I,pa)
               iota2_IN = iota2_IN + iota2(K,L,IRK)*PHI_AREA(K,I,pa)
#endif

#ifdef DYNP
               dynP_IN = dynP_IN + dynP(K,L,IRK)*PHI_AREA(K,I,pa)
#endif

               DEPTH = ZE_IN + HB_IN

#ifdef SED_LAY
               do ll = 1,layers
                  bed_IN(ll) = bed_IN(ll) + bed(K,L,IRK,ll)*PHI_AREA(K,I,pa)
                  HB_IN = HB_IN + bed(K,L,irk,ll)*PHI_AREA(K,I,pa)

                  MZ_X(ll) =  MZ_X(ll) + MZ(K,1,ll,L)*PHI_AREA(K,I,pa)
                  MZ_Y(ll) =  MZ_Y(ll) + MZ(K,2,ll,L)*PHI_AREA(K,I,pa)
               enddo
               DEPTH = 0.D0
               DEPTH = ZE_IN + HB_IN
C.........Compute sediment discharge model                                                                           

               !Note that the choice of linearization can require this to be changed
            QMag_IN = (QX_IN*QX_IN/(DEPTH**2) + QY_IN*QY_IN/(DEPTH)**2)**(1/2)
            discharge_modelX_IN = porosity * DEPTH**(-1) * QMag_IN**(2) * QX_IN*SFACQUAD                             
            discharge_modelY_IN = porosity * DEPTH**(-1) * QMag_IN**(2) *QY_IN  

#endif
               
            ENDDO
 
C.........Compute continuity fluxes

            F1_NL = NLEQ + LEQ*HB_IN

#ifdef WAVE_DIF 
            FX_IN = (QX_IN+HZ_X)*F1_NL*SFACQUAD
#else
            FX_IN = QX_IN*F1_NL*SFACQUAD
#endif

#ifdef WAVE_DIF 
            FY_IN = (QY_IN+HZ_Y)*F1_NL
#else
            FY_IN = QY_IN*F1_NL
#endif
C.........Compute momentum flux terms

            FU_NL = NLEQ*QX_IN
            FV_NL = NLEQ*QY_IN
            FG_NL = NLEQG*ZE_IN*WDFLG(L)
            FH_NL = 1.D0/(NLEQ*DEPTH + LEQ)
            U_IN  = QX_IN*FH_NL
            V_IN  = QY_IN*FH_NL

            HUU = FU_NL*U_IN
            HVV = FV_NL*V_IN
            HUV = FU_NL*V_IN
            GH2 = FG_NL*(0.5D0*ZE_IN + HB_IN) + FG_L*ZE_IN
#ifdef SED_LAY
            !Not well-balanced, be careful here!
            !GH2 =  0.D0
            !GH2 =  0.5D0*G*(DEPTH**2)
            !GH2 =  WDFLG(L)*0.5D0*G*(DEPTH**2)
            !GH2 = FG_NL*(0.5D0*ZE_IN + 2.D0*HB_IN) + 0.5D0*FG_L* (ZE_IN**2 - DEPTH**2
#endif
           
C.........Compute x momentum fluxes

            GX_IN = (HUU + GH2 + LZ_XX)*SFACQUAD
            GY_IN = HUV + LZ_XY

            advectqx(l)=advectqx(l)+gx_in+gy_in

C.........Compute y momentum fluxes

            HX_IN = (HUV + LZ_YX)*SFACQUAD
            HY_IN = HVV + GH2 + LZ_YY

            advectqy(l)=advectqy(l)+hx_in+hy_in

C.........Compute the friction factor
            if (LoadManningsN) then
c     MN_IN=MANN(1,L)
c     do k=2,dof
c     MN_IN=MN_IN + MANN(K,L)*PHI_AREA(K,I)
c     enddo
               fric_el(L)=G*
     $              ((ManningsN(n1)+ManningsN(n2)+ManningsN(n3))/3.)**2
c     $            MN_IN**2
     $              /(DEPTH**(1.d0/3.d0))
               if (fric_el(L).lt.CF) fric_el(L)=CF
            endif
            TAU = FRIC_EL(L)
C     IF (IFLINBF.EQ.0) THEN
C     UMAG = SQRT( U_IN*U_IN + V_IN*V_IN )
C     TAU  = TAU*UMAG*FH_NL
C     IF (IFHYBF.EQ.1) TAU = TAU*
C     &             (1.D0  + (HBREAK*FH_NL)**FTHETA)**(FGAMMA/FTHETA)
C     ENDIF
C     Modified to compute TAU using elemental averages.
C     This seems necessary to avoid exessive bottom friction
C     at wetting-drying fronts where the total column height is very
C     small. S.B. 9-Feb-2008

            IF (IFLINBF.EQ.0) THEN
               UMAG = SQRT( U_IN*U_IN + V_IN*V_IN )
c     cnd modified 4/23/10 to test friction 
c     TAU  = TAU*UMAG_C*FH_NL_C
               TAU  = TAU*UMAG*FH_NL
               IF (IFHYBF.EQ.1) TAU = TAU*
     &              (1.D0  + (HBREAK*FH_NL_C)**FTHETA)**(FGAMMA/FTHETA)
c     &            (1.D0  + (HBREAK*FH_NL)**FTHETA)**(FGAMMA/FTHETA)
C     It is numerically probable that the bottom friction is large enoght
C     to reverse the direction of currents backward due to a too small column
C     height even though it does not happen in reality. To avoid this, the MIN
C     function bellow is added. It is expected that this MIN function upper-limits
C     TAU so the bottom friction force does not reverse the currents within
C     half a time step.  S.B. 9-Feb-2008
c     IF(TAU.GT.2.D0*DTDPH) PRINT *, "TAU = ", TAU, 2.D0*DTDPH,
c     $           u_in,v_in,myproc,L
c     TAU = MIN(TAU, 2.D0*DTDPH)
               TAU = MIN(TAU, .9D0*DTDPH)
            ENDIF
c     IF (RAMPDG.LT.1.D0) TAU = MAX(TAU,0.001)

C.........Compute the x momentum source/sink terms

            SOURCE_X =

C.........1.) Friction term

     &           - TAU*QX_IN
      
C.........2.) Bathymetric slope term

     &           + FG_NL*DHB_X*SFACQUAD

C.........3.) Coriolis force

     &           + CORI_EL(L)*QY_IN

C.....Compute the y momentum source/sink terms

            SOURCE_Y =

C.........1.) Friction term

     &           - TAU*QY_IN

C.........2) Bathymetric slope term

      
     &           + FG_NL*DHB_Y
     
C.........3.) Coriolis force

     &           - CORI_EL(L)*QX_IN
            
C.........4.) Wind and pressure forcing (in x and y)

            IF (NWS.NE.0) THEN
               FW_NL = 1.D0/F1_NL
               SOURCE_X = SOURCE_X + FW_NL*( WSX2(N1)*PSI1(I,pa)
     &              + WSX2(N2)*PSI2(I,pa)  + WSX2(N3)*PSI3(I,pa) )
c     &                          - G*DEPTH*( PR2(N1)*DPSIDX(1)
c     nd
     &              - G*SFACQUAD*DEPTH
     $              *( PR2(N1)*DPSIDX(1)
     &              + PR2(N2)*DPSIDX(2) + PR2(N3)*DPSIDX(3))
               SOURCE_Y = SOURCE_Y + FW_NL*( WSY2(N1)*PSI1(I,pa)
     &              + WSY2(N2)*PSI2(I,pa)  + WSY2(N3)*PSI3(I,pa) )
     &              - G*DEPTH*( PR2(N1)*DPSIDY(1)
     &              + PR2(N2)*DPSIDY(2) + PR2(N3)*DPSIDY(3))
            ENDIF

c     if (myproc.eq.1.and.l.eq.440.and.i.eq.1) then
c     write(440,*) 'tau ',tau,umag,fric_el(l)
c     endif

C.........5) Tidal potential forcing (in x and y)

            IF (NTIP.NE.0) THEN
c$$$               SOURCE_X = SOURCE_X + RAMPDG*G*DEPTH*SFAC_ELEM(I,L,pa)*
c$$$     &              ( DPSIDX(1)*TIP2(N1)
c$$$     &              + DPSIDX(2)*TIP2(N2) + DPSIDX(3)*TIP2(N3) )

               SOURCE_X = SOURCE_X + 
c     $            RAMPDG*G*DEPTH*( DPSIDX(1)*TIP2(N1)
c     nd
     $              RAMPDG*G*DEPTH*SFACQUAD*( DPSIDX(1)*TIP2(N1)
     &              + DPSIDX(2)*TIP2(N2) + DPSIDX(3)*TIP2(N3) )
               
               SOURCE_Y = SOURCE_Y + RAMPDG*G*DEPTH*( DPSIDY(1)*TIP2(N1)
     &              + DPSIDY(2)*TIP2(N2) + DPSIDY(3)*TIP2(N3) )
            ENDIF

            sourceqx(l)=sourceqx(l)+source_x
            sourceqy(l)=sourceqy(l)+source_y

C.........6) Chemical mass action  

#ifdef CHEM     
            MassAction1st = 0.D0
            MassAction2nd = 0.D0

            MassMax(L) = 0.D0
            
            rate = reaction_rate*1.D0/86400.D0

            MassAction1st =  rate * (max(iota2_IN,0.D0)*max(iota_IN,0.D0))
            MassAction2nd =  rate * (max(iota2_IN,0.D0)*max(iota_IN,0.D0))

            MassMax(L) = min(MassAction1st,MassAction2nd)

C.........Build the rhs

            RHS_iota(1,L,IRK) = RHS_iota(1,L,IRK) 
     &           + MassAction1st*SRFAC(1,I,L,pa)*FH_NL
            RHS_iota2(1,L,IRK) = RHS_iota2(1,L,IRK) 
     &           + MassAction2nd*SRFAC(1,I,L,pa)*FH_NL  
#endif


            RHS_QX(1,L,IRK) = RHS_QX(1,L,IRK) + SRFAC(1,I,L,pa)*SOURCE_X
            RHS_QY(1,L,IRK) = RHS_QY(1,L,IRK) + SRFAC(1,I,L,pa)*SOURCE_Y
            
            DO K = 2,DOFS(L)

               RHS_ZE(K,L,IRK) = RHS_ZE(K,L,IRK) + XFAC(K,I,L,pa)*FX_IN
     &              + YFAC(K,I,L,pa)*FY_IN 
               RHS_QX(K,L,IRK) = RHS_QX(K,L,IRK) + XFAC(K,I,L,pa)*GX_IN
     &              + YFAC(K,I,L,pa)*GY_IN + SRFAC(K,I,L,pa)*SOURCE_X
               RHS_QY(K,L,IRK) = RHS_QY(K,L,IRK) + XFAC(K,I,L,pa)*HX_IN
     &              + YFAC(K,I,L,pa)*HY_IN + SRFAC(K,I,L,pa)*SOURCE_Y

#ifdef SED_LAY

               do ll = 1,layers !only really makes sense for single layer

                  RHS_bed(K,L,IRK,ll) = RHS_bed(K,L,IRK,ll) 
     &                 + XFAC(K,I,L,pa)*(discharge_modelX_IN+MZ_X(ll)*SFACQUAD)
     &                 + YFAC(K,I,L,pa)*(discharge_modelY_IN+MZ_Y(ll))

               enddo
#endif

#ifdef TRACE
               RHS_iota(K,L,IRK) = RHS_iota(K,L,IRK) 
     &              + XFAC(K,I,L,pa)*(iota_IN*QX_IN*FH_NL+TZ_X*SFACQUAD) 
     &              + YFAC(K,I,L,pa)*(iota_IN*QY_IN*FH_NL+TZ_Y) 
#endif

#ifdef CHEM
               RHS_iota(K,L,IRK) = RHS_iota(K,L,IRK) 
     &              + XFAC(K,I,L,pa)*iota_IN*QX_IN*FH_NL 
     &              + YFAC(K,I,L,pa)*iota_IN*QY_IN*FH_NL 
     &              + MassAction1st*SRFAC(K,I,L,pa)*FH_NL 
               
               RHS_iota2(K,L,IRK) = RHS_iota2(K,L,IRK) 
     &              + XFAC(K,I,L,pa)*iota2_IN*QX_IN*FH_NL      
     &              + YFAC(K,I,L,pa)*iota2_IN*QY_IN*FH_NL 
     &              + MassAction2nd*SRFAC(K,I,L,pa)*FH_NL 
#endif

#ifdef DYNP
               RHS_dynP(K,L,IRK) = RHS_dynP(K,L,IRK) 
     &              + XFAC(K,I,L,pa)*dynP_IN*QX_IN*FH_NL 
     &              + YFAC(K,I,L,pa)*dynP_IN*QY_IN*FH_NL
     &              + subphi_IN*SRFAC(K,I,L,pa)*FH_NL 
#endif
               
            ENDDO
            
         ENDDO  

         
 1000 CONTINUE      
      RETURN
      END SUBROUTINE
