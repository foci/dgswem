C***********************************************************************
C
C     SUBROUTINE SLOPELIMITER()
C
C     This subroutine selects either SLOPELIMITER1() or SLOPELIMITER2()
C     according to SLOPEFLAG. SLOPEFLAG is specified in fort.dg.
C   
C     All routines rewritten for p_adaptive multicomponent version
C     Slopelimiters 2 and 3 are not compatible. 1,4,5,6,7,8,9,10 are.
C     -- cem, 2011
C
C***********************************************************************

      SUBROUTINE SLOPELIMITER()

      USE DG, ONLY : SLOPEFLAG
      
      IMPLICIT NONE

#ifdef SLOPEALL
      IF (SLOPEFLAG.EQ.1) THEN
        CALL SLOPELIMITER1()
      ELSE IF (SLOPEFLAG.EQ.2) THEN
        CALL SLOPELIMITER2()
      ELSE IF (SLOPEFLAG.EQ.3) THEN
        CALL SLOPELIMITER3()
      ELSE IF (SLOPEFLAG.EQ.4) THEN
        CALL SLOPELIMITER4()
      ELSE IF (SLOPEFLAG.EQ.5) THEN
        CALL SLOPELIMITER5()
      ELSE IF (SLOPEFLAG.EQ.6) THEN
        CALL SLOPELIMITER6()
      ELSE IF (SLOPEFLAG.EQ.7) THEN
        CALL SLOPELIMITER7()
      ELSE IF (SLOPEFLAG.EQ.8) THEN
        CALL SLOPELIMITER8()
      ELSE IF (SLOPEFLAG.EQ.9) THEN
        CALL SLOPELIMITER9()
      ELSE IF (SLOPEFLAG.EQ.10) THEN
        CALL SLOPELIMITER10()
      ENDIF
#endif

#ifdef SLOPE5 
      IF (SLOPEFLAG .NE. 0) THEN
        CALL SLOPELIMITER5()
      ENDIF
#endif

#ifdef STBLZR
        CALL SLOPELIMITER5()
#endif
    

      RETURN
      END SUBROUTINE

#ifdef SLOPEALL

C***********************************************************************
C     
C     SUBROUTINE SLOPELIMITER1()
C     
C     Written by Ethan Kubatko
C     
C     08-Feb-2008
C     - Modified to apply this slopelimiter on domain 
C     boundaries. S.B.
C     28-Jun-2010 Modified for transport and chemistry and p enrichment
C     -cem
C     
C***********************************************************************

      SUBROUTINE SLOPELIMITER1()

C.....Use appropriate modules

      USE GLOBAL
      USE DG

      IMPLICIT NONE

C.....Declare local variables

      INTEGER EL2, EL3, L,i,j,kk,k
      
      REAL(SZ) DEN, DZEDR, DZEDS, DQXDR, DQXDS, DQYDR, DQYDS
      Real(SZ) DiotaDR,DiotaDS,Diota2DR,Diota2DS
      REAL(SZ) DZEDX(4), DZEDY(4)
      REAL(SZ) DQXDX(4), DQXDY(4)
      REAL(SZ) DQYDX(4), DQYDY(4)
      Real(SZ) DiotaDX(4), DiotaDY(4)
      Real(SZ) Diota2DX(4), Diota2DY(4)
      REAL(SZ) GRADZE(4), GRADQX(4), GRADQY(4)
      Real(SZ) GRADiota(4), GRADiota2(4)
      REAL(SZ) SL1, SL2
      REAL(SZ) ZELIM, ZELIMX, ZELIMY
      REAL(SZ) QXLIM, QXLIMX, QXLIMY
      REAL(SZ) QYLIM, QYLIMX, QYLIMY
      Real(SZ) iotaLIM, iotaLIMX, iotaLIMY
      Real(SZ) iota2LIM, iota2LIMX, iota2LIMY

C.....Save the original values  02/28/2007 sb

      DO L=1,NE
         if (dofs(l).eq.3) then
            DO K = 1,DOFS(L)
               ZE(K,L,NRK+2) = ZE(K,L,IRK+1)
               QX(K,L,NRK+2) = QX(K,L,IRK+1)
               QY(K,L,NRK+2) = QY(K,L,IRK+1)
               
#ifdef TRACE           
               iota(K,L,NRK+2) = iota(K,L,IRK+1)
#endif

#ifdef CHEM
               iota(K,L,NRK+2) = iota(K,L,IRK+1)
               iota2(K,L,NRK+2) = iota2(K,L,IRK+1)
#endif

            ENDDO
         endif
      ENDDO      

      DO 1000 L=1,ne

         if (dofs(l).eq.3) then

C.....Retrieve the barycenter coordinates of the element
            
            X1 = XBC(L)
            Y1 = YBC(L)

C.....Compute the x and y derivatives of the P1 part of ZE

            DZEDR = ZE(3,L,NRK+2)
            DZEDS = 3.D0/2.D0*ZE(2,L,NRK+2) + 1.D0/2.D0*ZE(3,L,NRK+2)

            DZEDX(4) = DZEDR*DRDX(L) + DZEDS*DSDX(L)
            DZEDY(4) = DZEDR*DRDY(L) + DZEDS*DSDY(L)
            
            GRADZE(4) = SQRT(DZEDX(4)**2.D0 + DZEDY(4)**2.D0)
            
C.....Compute the x and y derivatives of the P1 part of QX
            
            DQXDR = QX(3,L,NRK+2)
            DQXDS = 3.D0/2.D0*QX(2,L,NRK+2) + 1.D0/2.D0*QX(3,L,NRK+2)

            DQXDX(4) = DQXDR*DRDX(L) + DQXDS*DSDX(L)
            DQXDY(4) = DQXDR*DRDY(L) + DQXDS*DSDY(L)
            
            GRADQX(4) = SQRT(DQXDX(4)**2.D0 + DQXDY(4)**2.D0)
            
C.....Compute the x and y derivatives of the P1 part of QY
            
            DQYDR = QY(3,L,NRK+2)
            DQYDS = 3.D0/2.D0*QY(2,L,NRK+2) + 1.D0/2.D0*QY(3,L,NRK+2)

            DQYDX(4) = DQYDR*DRDX(L) + DQYDS*DSDX(L)
            DQYDY(4) = DQYDR*DRDY(L) + DQYDS*DSDY(L)
            
            GRADQY(4) = SQRT(DQYDX(4)**2.D0 + DQYDY(4)**2.D0)


C.....Compute the x and y derivatives of the P1 part of iota & iota2

#ifdef TRACE     
            DiotaDR = iota(3,L,NRK+2)
            DiotaDS =3.D0/2.D0*iota(2,L,NRK+2)+1.D0/2.D0*iota(3,L,NRK+2)

            DiotaDX(4) = DiotaDR*DRDX(L) + DiotaDS*DSDX(L)
            DiotaDY(4) = DiotaDR*DRDY(L) + DiotaDS*DSDY(L)
            
            GRADiota(4) = SQRT(DiotaDX(4)**2.D0 + DiotaDY(4)**2.D0)
#endif

#ifdef CHEM
            DiotaDR = iota(3,L,NRK+2)
            DiotaDS =3.D0/2.D0*iota(2,L,NRK+2)+1.D0/2.D0*iota(3,L,NRK+2)

            Diota2DR = iota2(3,L,NRK+2)
            Diota2DS = 3.D0/2.D0*iota2(2,L,NRK+2) + 1.D0/2.D0*
     &         iota2(3,L,NRK+2)

            DiotaDX(4) = DiotaDR*DRDX(L) + DiotaDS*DSDX(L)
            DiotaDY(4) = DiotaDR*DRDY(L) + DiotaDS*DSDY(L)

            Diota2DX(4) = Diota2DR*DRDX(L) + Diota2DS*DSDX(L)
            Diota2DY(4) = Diota2DR*DRDY(L) + Diota2DS*DSDY(L)
            
            GRADiota(4) = SQRT(DiotaDX(4)**2.D0 + DiotaDY(4)**2.D0)
            GRADiota2(4) = SQRT(Diota2DX(4)**2.D0 + Diota2DY(4)**2.D0)
#endif
            
            DO I = 1,3
               
               EL2 = EL_NBORS(I,L)
               EL3 = EL_NBORS(I+1,L)

C.......Special treatment for elements that do not have three neighboring elements
               IF((EL2.EQ.0).OR.(EL3.EQ.0).OR.(SL3(I,L).LE.0.D0)) THEN
                  GRADZE(I) = 1.D10
                  GRADQX(I) = 1.D10
                  GRADQY(I) = 1.D10

#ifdef TRACE
                  GRADiota(I) = 1.D10
#endif

#ifdef CHEM
                  GRADiota(I) = 1.D10
                  GRADiota2(I) = 1.D10
#endif

                  CYCLE
               ENDIF
               
               X2 = XBC(EL2)
               Y2 = YBC(EL2)
               
               X3 = XBC(EL3)
               Y3 = YBC(EL3)
               
C.....Compute limiting slopes for the surface elevation

               SL1 = ZE(1,L,NRK+2)*(Y3 - Y2) + ZE(1,EL2,NRK+2)*(Y1 - Y3)
     &              + ZE(1,EL3,NRK+2)*(Y2 - Y1)
               
               SL2 = ZE(1,L,NRK+2)*(X2 - X3) + ZE(1,EL2,NRK+2)*(X3 - X1)
     &              + ZE(1,EL3,NRK+2)*(X1 - X2)
               
               DZEDX(I) = -SL1/SL3(I,L)
               
               DZEDY(I) = -SL2/SL3(I,L)
               
               GRADZE(I) = SQRT(DZEDX(I)**2.D0 + DZEDY(I)**2.D0)
               
C.....Compute limiting slopes for the X flow

               SL1 = QX(1,L,NRK+2)*(Y3 - Y2) + QX(1,EL2,NRK+2)*(Y1 - Y3)
     &              + QX(1,EL3,NRK+2)*(Y2 - Y1)

               SL2 = QX(1,L,NRK+2)*(X2 - X3) + QX(1,EL2,NRK+2)*(X3 - X1)
     &              + QX(1,EL3,NRK+2)*(X1 - X2)

               DQXDX(I) = -SL1/SL3(I,L)

               DQXDY(I) = -SL2/SL3(I,L)
               
               GRADQX(I) = SQRT(DQXDX(I)**2.D0 + DQXDY(I)**2.D0)
               
C.....Compute limiting slopes for the Y flow

               SL1 = QY(1,L,NRK+2)*(Y3 - Y2) + QY(1,EL2,NRK+2)*(Y1 - Y3)
     &              + QY(1,EL3,NRK+2)*(Y2 - Y1)

               SL2 = QY(1,L,NRK+2)*(X2 - X3) + QY(1,EL2,NRK+2)*(X3 - X1)
     &              + QY(1,EL3,NRK+2)*(X1 - X2)

               DQYDX(I) = -SL1/SL3(I,L)

               DQYDY(I) = -SL2/SL3(I,L)
               
               GRADQY(I) = SQRT(DQYDX(I)**2.D0 + DQYDY(I)**2.D0)
               

C.....Compute limiting slopes for the concentrations

#ifdef TRACE
               SL1 = iota(1,L,NRK+2)*(Y3 - Y2)+iota(1,EL2,NRK+2)*(Y1-Y3)
     &              + iota(1,EL3,NRK+2)*(Y2 - Y1)
               
               SL2 = iota(1,L,NRK+2)*(X2 - X3)+iota(1,EL2,NRK+2)*(X3-X1)
     &              + iota(1,EL3,NRK+2)*(X1 - X2)
               
               DiotaDX(I) = -SL1/SL3(I,L)
               
               DiotaDY(I) = -SL2/SL3(I,L)
               
               GRADiota(I) = SQRT(DiotaDX(I)**2.D0 + DiotaDY(I)**2.D0)
#endif

#ifdef CHEM
               SL1 = iota(1,L,NRK+2)*(Y3 - Y2)+iota(1,EL2,NRK+2)*(Y1-Y3)
     &              + iota(1,EL3,NRK+2)*(Y2 - Y1)
               
               SL2 = iota(1,L,NRK+2)*(X2 - X3)+iota(1,EL2,NRK+2)*(X3-X1)
     &              + iota(1,EL3,NRK+2)*(X1 - X2)
               
               DiotaDX(I) = -SL1/SL3(I,L)
               
               DiotaDY(I) = -SL2/SL3(I,L)
               
               GRADiota(I) = SQRT(DiotaDX(I)**2.D0 + DiotaDY(I)**2.D0)

               SL1=iota2(1,L,NRK+2)*(Y3 - Y2)+iota2(1,EL2,NRK+2)*(Y1-Y3)
     &              + iota2(1,EL3,NRK+2)*(Y2 - Y1)
               
               SL2=iota2(1,L,NRK+2)*(X2 - X3)+iota2(1,EL2,NRK+2)*(X3-X1)
     &              + iota2(1,EL3,NRK+2)*(X1 - X2)
               
               Diota2DX(I) = -SL1/SL3(I,L)
               
               Diota2DY(I) = -SL2/SL3(I,L)
               
               GRADiota2(I) = SQRT(Diota2DX(I)**2.D0+Diota2DY(I)**2.D0)
#endif

            ENDDO
            
C.....Compute the (possibly) limited gradient

            ZELIM = MIN( GRADZE(1), GRADZE(2), GRADZE(3), GRADZE(4) )
            QXLIM = MIN( GRADQX(1), GRADQX(2), GRADQX(3), GRADQX(4) )
            QYLIM = MIN( GRADQY(1), GRADQY(2), GRADQY(3), GRADQY(4) )

#ifdef TRACE
           iotaLIM = MIN(GRADiota(1), GRADiota(2), GRADiota(3),
     &         GRADiota(4))
#endif

#ifdef CHEM
            iotaLIM = MIN( GRADiota(1), GRADiota(2), GRADiota(3), 
     &         GRADiota(4) )
            iota2LIM = MIN( GRADiota2(1), GRADiota2(2), GRADiota2(3), 
     &         GRADiota2(4) )
#endif
            
            DO I = 1,4
               
               IF (ZELIM.EQ.GRADZE(I)) THEN
                  ZELIMX = DZEDX(I)
                  ZELIMY = DZEDY(I)
               ENDIF
               
               IF (QXLIM.EQ.GRADQX(I)) THEN
                  QXLIMX = DQXDX(I)
                  QXLIMY = DQXDY(I)
               ENDIF
               
               IF (QYLIM.EQ.GRADQY(I)) THEN
                  QYLIMX = DQYDX(I)
                  QYLIMY = DQYDY(I)
               ENDIF

#ifdef TRACE
               IF (iotaLIM.EQ.GRADiota(I)) THEN
                  iotaLIMX = DiotaDX(I)
                  iotaLIMY = DiotaDY(I)
               ENDIF
#endif

#ifdef CHEM
               IF (iotaLIM.EQ.GRADiota(I)) THEN
                  iotaLIMX = DiotaDX(I)
                  iotaLIMY = DiotaDY(I)
               ENDIF

               IF (iota2LIM.EQ.GRADiota2(I)) THEN
                  iota2LIMX = Diota2DX(I)
                  iota2LIMY = Diota2DY(I)
               ENDIF
#endif
            ENDDO

            DEN = DRDX(L)*DSDY(L) - DRDY(L)*DSDX(L)

            ZE(2,L,IRK+1) = 1.D0/3.D0*(2.D0*ZELIMY*DRDX(L)
     &           -2.D0*DRDY(L)*ZELIMX-ZELIMX*DSDY(L)+DSDX(L)*ZELIMY)/DEN
            QX(2,L,IRK+1) = 1.D0/3.D0*(2.D0*QXLIMY*DRDX(L)
     &           -2.D0*DRDY(L)*QXLIMX-QXLIMX*DSDY(L)+DSDX(L)*QXLIMY)/DEN
            QY(2,L,IRK+1) = 1.D0/3.D0*(2.D0*QYLIMY*DRDX(L)
     &           -2.D0*DRDY(L)*QYLIMX-QYLIMX*DSDY(L)+DSDX(L)*QYLIMY)/DEN

            ZE(3,L,IRK+1) = (ZELIMX*DSDY(L) - DSDX(L)*ZELIMY)/DEN
            QX(3,L,IRK+1) = (QXLIMX*DSDY(L) - DSDX(L)*QXLIMY)/DEN
            QY(3,L,IRK+1) = (QYLIMX*DSDY(L) - DSDX(L)*QYLIMY)/DEN

#ifdef TRACE
            iota(2,L,IRK+1) = 1.D0/3.D0*(2.D0*iotaLIMY*DRDX(L)
     &           - 2.D0*DRDY(L)*iotaLIMX - iotaLIMX*DSDY(L) + DSDX(L)*
     &           iotaLIMY)/DEN

            iota(3,L,IRK+1) = (iotaLIMX*DSDY(L) - DSDX(L)*iotaLIMY)/DEN
#endif

#ifdef CHEM
            iota(2,L,IRK+1) = 1.D0/3.D0*(2.D0*iotaLIMY*DRDX(L)
     &           - 2.D0*DRDY(L)*iotaLIMX - iotaLIMX*DSDY(L) + DSDX(L)*
     &           iotaLIMY)/DEN
            iota2(2,L,IRK+1) = 1.D0/3.D0*(2.D0*iota2LIMY*DRDX(L)
     &           - 2.D0*DRDY(L)*iota2LIMX - iota2LIMX*DSDY(L) + DSDX(L)*
     &           iota2LIMY)/DEN

            iota(3,L,IRK+1) = (iotaLIMX*DSDY(L) - DSDX(L)*iotaLIMY)/DEN
            iota2(3,L,IRK+1) = (iota2LIMX*DSDY(L)-DSDX(L)*iota2LIMY)/DEN
#endif
            
 111        continue

         endif

 1000 CONTINUE
      RETURN
      END SUBROUTINE

      
C***********************************************************************
C     
C     SUBROUTINE SLOPELIMITER2()
C     
C     Written by Shintaro Bunya - 08 Feb 2008
C     
C     This subroutine is an implementation of the slope limiter used in
C     B. Cockburn and C-W. Shu, "The Runge-Kutta Discontinuous Galerkin
C     Method for Conservation Laws V, Multidimensional Systems,"
C     Journal of Computational Physics 141, 199-224 (1998).
C     
C***********************************************************************

      SUBROUTINE SLOPELIMITER2()

C.....Use appropriate modules

      USE SIZES, ONLY : SZ
      USE GLOBAL
      USE DG

      IMPLICIT NONE

C.....Declare local variables

      INTEGER L, GED, II,i,j,kk,k
      INTEGER EL2, EL3, EL_N, EL_N1, EL_N2
      
      REAL(SZ) DTM, DZEDR, DZEDS, DQXDR, DQXDS, DQYDR, DQYDS
      REAL(SZ) DZEDX(4), DZEDY(4)
      REAL(SZ) DQXDX(4), DQXDY(4)
      REAL(SZ) DQYDX(4), DQYDY(4)
      REAL(SZ) DiotaDX(4), DiotaDY(4)
      REAL(SZ) GRADZE(4), GRADQX(4), GRADQY(4)
      REAL(SZ) SL1, SL2
      REAL(SZ) ZELIM, ZELIMX, ZELIMY
      REAL(SZ) QXLIM, QXLIMX, QXLIMY
      REAL(SZ) QYLIM, QYLIMX, QYLIMY
      REAL(SZ) XB(0:3), YB(0:3)
      REAL(SZ) XM(1:3), YM(1:3)
      REAL(SZ) ZB(0:3),DB(0:3),QXB(0:3),QYB(0:3) ! Values at the barycenters
      REAL(SZ) ZC(1:3),DC(1:3),QXC(1:3),QYC(1:3) ! Values at the corner nodes
      REAL(SZ) ZM(1:3),DM(1:3),QXM(1:3),QYM(1:3) ! Values at the edge midpoints
      REAL(SZ) HTB              ! Total height at the barycenter
      REAL(SZ) NXX(3), NYY(3)
      REAL(SZ) ALPHA1(3), ALPHA2(3)
      REAL(SZ) A(2,2), B(2)
      REAL(SZ) TE(3,3)
      REAL(SZ) WB_O(3,0:2), WM_O(3) ! Variable vectors in original space
      REAL(SZ) WB_C(3,0:2), WM_C(3) ! Variable vectors in characteristic space
      REAL(SZ) WC_O(3,3)        ! Variable vectors in original space at the corner nodes
      REAL(SZ) W_TILDA(3), DELTA_W(3)
      REAL(SZ) R_SQ             ! Squared r
      REAL(SZ) S
      REAL(SZ) DELTA_O(3,3), DELTA_C(3,3), DELTA_HAT(3,3)
      INTEGER N(3)
      REAL(SZ) POS, NEG, THETA_P, THETA_N
      LOGICAL LIMITIT(3)
      REAL(SZ), PARAMETER :: ZERO = 1.D-15

C.....Save the original values  02/28/2007 sb
      DO L=1,NE
         DO K = 1,DOFS(L)
            ZE(K,L,NRK+2) = ZE(K,L,IRK+1)
            QX(K,L,NRK+2) = QX(K,L,IRK+1)
            QY(K,L,NRK+2) = QY(K,L,IRK+1)
         ENDDO
      ENDDO      

      DO 1000 L=1, NE

         pa = PDG_EL(L)

         if (pa.eq.0) then
            pa = 1
         endif

C.....Retrieve the barycenter info of this element
         XB(0) = XBC(L)
         YB(0) = YBC(L)
         ZB(0) = ZE(1,L,IRK+1)
         DB(0) = HB(1,L,1)
         QXB(0) = QX(1,L,IRK+1)
         QYB(0) = QY(1,L,IRK+1)
         
C.....Retrieve the nodal info

         N(1) = NM(L,1)
         N(2) = NM(L,2)
         N(3) = NM(L,3)

         DO K = 1,3
            ZC(K) = ZE(1,L,IRK+1)
            DC(K) = HB(1,L,1)
            QXC(K) = QX(1,L,IRK+1)
            QYC(K) = QY(1,L,IRK+1)
            DO KK = 2,3         ! This slope limiter assumes the functions to limit as a linear function.
               ZC(K) = ZC(K) + PHI_CORNER(KK,K,pa)*ZE(KK,L,IRK+1)
               DC(K) = DC(K) + PHI_CORNER(KK,K,pa)*HB(KK,L,1)
               QXC(K) = QXC(K) + PHI_CORNER(KK,K,pa)*QX(KK,L,IRK+1)
               QYC(K) = QYC(K) + PHI_CORNER(KK,K,pa)*QY(KK,L,IRK+1)
            ENDDO
         ENDDO

C.....Retrieve info at barycenters and edge midpoints
         DO I = 1,3
            EL_N = EL_NBORS(I,L)

C.......Edge midpoints
            N1 = N(MOD(I+0,3)+1)
            N2 = N(MOD(I+1,3)+1)
            XM(I) = 0.5D0*(X(N1) + X(N2))
            YM(I) = 0.5D0*(Y(N1) + Y(N2))

            N1 = MOD(I+0,3)+1
            N2 = MOD(I+1,3)+1
            ZM(I) = 0.5D0*(ZC(N1) + ZC(N2))
            DM(I) = 0.5D0*(DC(N1) + DC(N2))
            QXM(I) = 0.5D0*(QXC(N1) + QXC(N2))
            QYM(I) = 0.5D0*(QYC(N1) + QYC(N2))
            
C.......Barycenters
            IF(EL_N == 0) THEN
C     If there is no neighboring element on this side (e.g. domain boundary),
C     locate a barycenter at the mirroring point and use the values at the 
C     baricenter of this element.
               XB(I) = 2.D0*XM(I) - XB(0)
               YB(I) = 2.D0*YM(I) - YB(0)
               ZB(I) = ZB(0)
               DB(I) = DB(0)
               QXB(I) = QXB(0)
               QYB(I) = QYB(0)
            ELSE
               XB(I) = XBC(EL_N)
               YB(I) = YBC(EL_N)
               ZB(I) = ZE(1,EL_N,IRK+1)
               DB(I) = HB(1,EL_N,1)
               QXB(I) = QX(1,EL_N,IRK+1)
               QYB(I) = QY(1,EL_N,IRK+1)
            ENDIF

C.......Edge normal vector
            GED = NELED(I,L)
            NXX(I) = COSNX(GED)
            NYY(I) = SINNX(GED)
         ENDDO

C.....Compute alpha1 and alpha2      
         DO I = 1,3
            EL_N1 = I
            EL_N2 = MOD(I,3)+1
            A(1,1) = XB(EL_N1) - XB(0)
            A(1,2) = XB(EL_N2) - XB(0)
            A(2,1) = YB(EL_N1) - YB(0)
            A(2,2) = YB(EL_N2) - YB(0)
            B(1) = XM(I) - XB(0)
            B(2) = YM(I) - YB(0)
            DTM = A(1,1)*A(2,2) - A(1,2)*A(2,1) ! Compute the determinant
            ALPHA1(I) = ( A(2,2)*B(1) - A(1,2)*B(2))/DTM
            ALPHA2(I) = (-A(2,1)*B(1) + A(1,1)*B(2))/DTM
         ENDDO


C.....Set the variable vectors at the baricenter of element 0
         WB_O(1,0) = ZB(0)
         WB_O(2,0) = QXB(0)
         WB_O(3,0) = QYB(0)

C.....Start computing deltas for each edge
         DO I = 1,3

C.......Compute eigen values and vectors at the midpoint using the variables at the barycenter
            HTB = ZB(0)+DB(0)
            CALL HYDRO_EIGEN_VALUES(HTB, QXB(0)/HTB, QYB(0)/HTB,
     &           NXX(I), NYY(I), EIGVAL, RI, LE)

C.......Set variable vectors
            WM_O(1) = ZM(I)
            WM_O(2) = QXM(I)
            WM_O(3) = QYM(I)

            EL_N1 = I
            EL_N2 = MOD(I,3)+1

            WB_O(1,1) = ZB(EL_N1)
            WB_O(2,1) = QXB(EL_N1)
            WB_O(3,1) = QYB(EL_N1)
            
            WB_O(1,2) = ZB(EL_N2)
            WB_O(2,2) = QXB(EL_N2)
            WB_O(3,2) = QYB(EL_N2)

C.......Transform original variables into the characteristic space
            DO J = 0,2
               WB_C(1,J)
     &            =LE(1,1)*WB_O(1,J)+LE(1,2)*WB_O(2,J)+LE(1,3)*WB_O(3,J)
               WB_C(2,J)
     &            =LE(2,1)*WB_O(1,J)+LE(2,2)*WB_O(2,J)+LE(2,3)*WB_O(3,J)
               WB_C(3,J)
     &            =LE(3,1)*WB_O(1,J)+LE(3,2)*WB_O(2,J)+LE(3,3)*WB_O(3,J)
            ENDDO
            WM_C(1)= LE(1,1)*WM_O(1) + LE(1,2)*WM_O(2) + LE(1,3)*WM_O(3)
            WM_C(2)= LE(2,1)*WM_O(1) + LE(2,2)*WM_O(2) + LE(2,3)*WM_O(3)
            WM_C(3)= LE(3,1)*WM_O(1) + LE(3,2)*WM_O(2) + LE(3,3)*WM_O(3)

C.......Compute W_TILDA
            W_TILDA(1) = WM_C(1) - WB_C(1,0)
            W_TILDA(2) = WM_C(2) - WB_C(2,0)
            W_TILDA(3) = WM_C(3) - WB_C(3,0)
            
C.......Compute DELTA_W
            DELTA_W(1)
     &           = ALPHA1(I)*(WB_C(1,1) - WB_C(1,0))
     &           + ALPHA2(I)*(WB_C(1,2) - WB_C(1,0))
            DELTA_W(2)
     &           = ALPHA1(I)*(WB_C(2,1) - WB_C(2,0))
     &           + ALPHA2(I)*(WB_C(2,2) - WB_C(2,0))
            DELTA_W(3)
     &           = ALPHA1(I)*(WB_C(3,1) - WB_C(3,0))
     &           + ALPHA2(I)*(WB_C(3,2) - WB_C(3,0))

C.......Apply the TVB modified minmod function
            R_SQ = (XM(I) - XB(0))*(XM(I) - XB(0))
     &           + (YM(I) - YB(0))*(YM(I) - YB(0))

            DO K = 1,3          ! Loop for components of variable vectors
               IF(ABS(W_TILDA(K)) < SL2_M*R_SQ) THEN
                  DELTA_C(K,I) = W_TILDA(K)
               ELSE
                  IF((W_TILDA(K)*DELTA_W(K) >= 0.D0)
     &                 .AND.(ABS(W_TILDA(K)) > ZERO)) THEN
                     S = W_TILDA(K)/ABS(W_TILDA(K))
                     DELTA_C(K,I) = S*MIN(
     &                    ABS(W_TILDA(K)),ABS(SL2_NYU*DELTA_W(K)))
                  ELSE
                     DELTA_C(K,I) = 0.D0
                  ENDIF
               ENDIF
            ENDDO

C.......Transform back to the original space
            DELTA_O(1,I) = RI(1,1)*DELTA_C(1,I) + RI(1,2)*DELTA_C(2,I)
     &           + RI(1,3)*DELTA_C(3,I)
            DELTA_O(2,I) = RI(2,1)*DELTA_C(1,I) + RI(2,2)*DELTA_C(2,I)
     &           + RI(2,3)*DELTA_C(3,I)
            DELTA_O(3,I) = RI(3,1)*DELTA_C(1,I) + RI(3,2)*DELTA_C(2,I)
     &           + RI(3,3)*DELTA_C(3,I)

         ENDDO                  ! Loop for edges ends here

C     Set the variable vector at the midpoint 1. This is used later to judge
C     whether the variables should be limited.
         WM_O(1) = ZM(1)
         WM_O(2) = QXM(1)
         WM_O(3) = QYM(1)

C     Find the possibly limited DELTA and variable vectors at the corner nodes
         DO I = 1,3             ! Loop for variable components
C     Find the possibly limited DELTA
            POS = 0.D0
            NEG = 0.D0
            DO J = 1,3          ! Loop for edges
               POS = POS + MAX(0.D0, DELTA_O(I,J))
               NEG = NEG + MAX(0.D0,-DELTA_O(I,J))
            ENDDO
            IF((POS < ZERO).OR.(NEG < ZERO)) THEN 
               DELTA_HAT(I,1) = 0.D0
               DELTA_HAT(I,2) = 0.D0
               DELTA_HAT(I,3) = 0.D0
            ELSE
               THETA_P = MIN(1.D0,NEG/POS)
               THETA_N = MIN(1.D0,POS/NEG)
               
               DELTA_HAT(I,1)
     &              = THETA_P*MAX(0.D0, DELTA_O(I,1))
     &              - THETA_N*MAX(0.D0,-DELTA_O(I,1))
               DELTA_HAT(I,2)
     &              = THETA_P*MAX(0.D0, DELTA_O(I,2))
     &              - THETA_N*MAX(0.D0,-DELTA_O(I,2))
               DELTA_HAT(I,3)
     &              = THETA_P*MAX(0.D0, DELTA_O(I,3))
     &              - THETA_N*MAX(0.D0,-DELTA_O(I,3))
            ENDIF

C     Compute the possibly limited variable vectors at the corner nodes
            IF(( DELTA_HAT(I,1) < (WM_O(I) - WB_O(I,0) - ZERO) ).OR.
     &           ( DELTA_HAT(I,1) > (WM_O(I) - WB_O(I,0) + ZERO) )) THEN
C     Update the variables only if they need to be limited.
C     If the condition in the if statement is false, the following 
C     isn't applied and therefore the quradratic and higher order 
C     coefficients are preserved.

C     Computing limited variables at the corner nodes
               WC_O(I,1) = WB_O(I,0)
     &              - DELTA_HAT(I,1)
     &              + DELTA_HAT(I,2)
     &              + DELTA_HAT(I,3)
               WC_O(I,2) = WB_O(I,0)
     &              + DELTA_HAT(I,1)
     &              - DELTA_HAT(I,2)
     &              + DELTA_HAT(I,3)
               WC_O(I,3) = WB_O(I,0)
     &              + DELTA_HAT(I,1)
     &              + DELTA_HAT(I,2)
     &              - DELTA_HAT(I,3)
               LIMITIT(I) = .TRUE.
            ELSE
               LIMITIT(I) = .FALSE.
            ENDIF
         ENDDO

         IF(LIMITIT(1)) THEN
            ZE(2,L,IRK+1) = -1.D0/6.D0*(WC_O(1,1) + WC_O(1,2))
     &           + 1.D0/3.D0*WC_O(1,3)
            ZE(3,L,IRK+1) = -0.5D0*WC_O(1,1) + 0.5D0*WC_O(1,2)
         ENDIF
         IF(LIMITIT(2)) THEN
            QX(2,L,IRK+1) = -1.D0/6.D0*(WC_O(2,1) + WC_O(2,2))
     &           + 1.D0/3.D0*WC_O(2,3)
            QX(3,L,IRK+1) = -0.5D0*WC_O(2,1) + 0.5D0*WC_O(2,2)
         ENDIF
         IF(LIMITIT(3)) THEN
            QY(2,L,IRK+1) = -1.D0/6.D0*(WC_O(3,1) + WC_O(3,2))
     &           + 1.D0/3.D0*WC_O(3,3)
            QY(3,L,IRK+1) = -0.5D0*WC_O(3,1) + 0.5D0*WC_O(3,2)
         ENDIF

 1000 CONTINUE

      RETURN
      END SUBROUTINE

C***********************************************************************
C     
C     SUBROUTINE SLOPELIMITER3()
C     
C     Written by Shintaro Bunya - 27 Feb 2008
C     
C     This subroutine is an extentioin of the slope limiter used in
C     B. Cockburn and C-W. Shu, "The Runge-Kutta Discontinuous Galerkin
C     Method for Conservation Laws V, Multidimensional Systems,"
C     Journal of Computational Physics 141, 199-224 (1998).
C     
C     An additional parameter SL3_MD is used in this slope limiter to
C     make sure that slopes are limited in elements whose water depth
C     is small. (A special treatment for wetting and drying.
C     
C***********************************************************************

      SUBROUTINE SLOPELIMITER3()

C.....Use appropriate modules

      USE SIZES, ONLY : SZ
      USE GLOBAL
      USE DG

      IMPLICIT NONE

C.....Declare local variables

      INTEGER L, NNBORS,i,j,kk

      DO L= 1,NE
c     
         NNBORS = 0
         DO I=1,3
            IF(EL_NBORS(I,L).NE.0) THEN
               NNBORS = NNBORS + 1
            ENDIF
         ENDDO
         IF(NNBORS.EQ.3) THEN
            CALL SLOPELIMITER3_3NBORS(L)
         ELSE IF(NNBORS.EQ.2) THEN
            CALL SLOPELIMITER3_2NBORS(L)
         ELSE IF(NNBORS.EQ.1) THEN
            CALL SLOPELIMITER3_1NBORS(L)
#ifdef CMPI
C     CYCLE
#else
C     STOP 'STILL UNDER CONSTRUCTION!!!'
#endif
         ENDIF
      ENDDO

      RETURN
      END SUBROUTINE

      SUBROUTINE SLOPELIMITER3_3NBORS(L)

C.....Use appropriate modules

      USE SIZES, ONLY : SZ,myproc
      USE GLOBAL
      USE DG


      IMPLICIT NONE

C.....Declare local variables

      INTEGER,INTENT(IN) ::  L
      INTEGER GED, II,i,j,kk,k
      INTEGER EL2, EL3, EL_N, EL_N1, EL_N2
      
      REAL(SZ) DTM, DZEDR, DZEDS, DQXDR, DQXDS, DQYDR, DQYDS
      REAL(SZ) DZEDX(4), DZEDY(4)
      REAL(SZ) DQXDX(4), DQXDY(4)
      REAL(SZ) DQYDX(4), DQYDY(4)
      REAL(SZ) GRADZE(4), GRADQX(4), GRADQY(4)
      REAL(SZ) SL1, SL2
      REAL(SZ) ZELIM, ZELIMX, ZELIMY
      REAL(SZ) QXLIM, QXLIMX, QXLIMY
      REAL(SZ) QYLIM, QYLIMX, QYLIMY
      REAL(SZ) XB(0:3), YB(0:3)
      REAL(SZ) XM(1:3), YM(1:3)
      REAL(SZ) ZB(0:3),DB(0:3),QXB(0:3),QYB(0:3) ! Values at the barycenters
      REAL(SZ) ZC(1:3),DC(1:3),QXC(1:3),QYC(1:3) ! Values at the corner nodes
      REAL(SZ) ZM(1:3),DM(1:3),QXM(1:3),QYM(1:3) ! Values at the edge midpoints
      REAL(SZ) HTB              ! Total height at the barycenter
      REAL(SZ) NXX(3), NYY(3)
      REAL(SZ) ALPHA1(3), ALPHA2(3)
      REAL(SZ) A(2,2), B(2)
      REAL(SZ) TE(3,3)
      REAL(SZ) WB_O(3,0:2), WM_O(3) ! Variable vectors in original space
      REAL(SZ) WB_C(3,0:2), WM_C(3) ! Variable vectors in characteristic space
      REAL(SZ) WC_O(3,3), WC_C(3,3) ! Variable vectors in original space at the corner nodes
      REAL(SZ) W_TILDA(3), DELTA_W(3)
      REAL(SZ) R_SQ             ! Squared r
      REAL(SZ) HMIN, HMAX, HAVG, HLIMIT
      REAL(SZ) S
      REAL(SZ) DELTA_O(3,3), DELTA_C(3,3), DELTA_HAT(3,3)
      INTEGER N(3), NOT_LIMITED(3)
      REAL(SZ) POS, NEG, THETA_P, THETA_N
      LOGICAL LIMITIT(3)
      REAL(SZ), PARAMETER :: ZERO = 1.D-15

C.....Do nothing if this or a neighboring element is dry.
C     IF(WDFLG(L).EQ.0) RETURN
C     DO I = 1,3
C     EL_N = EL_NBORS(I,L)
C     IF(WDFLG(EL_N).EQ.0) RETURN
C     ENDDO
      
C.....Retrieve the barycenter info of this element
      XB(0) = XBC(L)
      YB(0) = YBC(L)
      ZB(0) = ZE(1,L,IRK+1)
      DB(0) = HB(1,L,1)
      QXB(0) = QX(1,L,IRK+1)
      QYB(0) = QY(1,L,IRK+1)

      pa = PDG_EL(L)

      if (pa.eq.0) then
         pa = 1
      endif

C.....Retrieve the nodal info

      N(1) = NM(L,1)
      N(2) = NM(L,2)
      N(3) = NM(L,3)

      DO K = 1,3
         ZC(K) = ZE(1,L,IRK+1)
         DC(K) = HB(1,L,1)
         QXC(K) = QX(1,L,IRK+1)
         QYC(K) = QY(1,L,IRK+1)         

         DO KK = 2,3            ! This slope limiter assumes the functions to limit as a linear function.
            ZC(K) = ZC(K) + PHI_CORNER(KK,K,pa)*ZE(KK,L,IRK+1)
            DC(K) = DC(K) + PHI_CORNER(KK,K,pa)*HB(KK,L,1)
            QXC(K) = QXC(K) + PHI_CORNER(KK,K,pa)*QX(KK,L,IRK+1)
            QYC(K) = QYC(K) + PHI_CORNER(KK,K,pa)*QY(KK,L,IRK+1)

         ENDDO
      ENDDO

      HMIN = MAX(0.D0,MIN(ZC(1)+DC(1),ZC(2)+DC(2),ZC(3)+DC(3)))
      HMAX = MAX(ZC(1)+DC(1),ZC(2)+DC(2),ZC(3)+DC(3))
      HAVG = ZE(1,L,IRK+1) + HB(1,L,1)
      HLIMIT = MIN(1.D0,(HMIN*HMIN)/(SL3_MD*SL3_MD))
      HLIMIT = 1.D0
      IF(nolifa.eq.2.and.HMIN.LE.H0) RETURN
c     IF(HMIN.LE.0.01D0) RETURN

C.....Retrieve info at barycenters and edge midpoints
      DO I = 1,3
         EL_N = EL_NBORS(I,L)

C.......Edge midpoints
         N1 = N(MOD(I+0,3)+1)
         N2 = N(MOD(I+1,3)+1)
         XM(I) = 0.5D0*(X(N1) + X(N2))
         YM(I) = 0.5D0*(Y(N1) + Y(N2))

         N1 = MOD(I+0,3)+1
         N2 = MOD(I+1,3)+1
         ZM(I) = 0.5D0*(ZC(N1) + ZC(N2))
         DM(I) = 0.5D0*(DC(N1) + DC(N2))
         QXM(I) = 0.5D0*(QXC(N1) + QXC(N2))
         QYM(I) = 0.5D0*(QYC(N1) + QYC(N2))
         
C.......Barycenters
C     IF(EL_N == 0) THEN
C     If there is no neighboring element on this side (e.g. domain boundary),
C     locate a barycenter at the mirroring point and use extrapolated values.
C     but all elements in this subroutine will have 3 neighbors
C     XB(I) = 2.D0*XM(I) - XB(0)
C     YB(I) = 2.D0*YM(I) - YB(0)
C     ZB(I) = 2.D0*ZM(I) - ZB(0)
C     DB(I) = 2.D0*DM(I) - DB(0)
C     QXB(I) = 2.D0*QXM(I) - QXB(0)
C     QYB(I) = 2.D0*QYM(I) - QYB(0)
C     ELSE
         XB(I) = XBC(EL_N)
         YB(I) = YBC(EL_N)
         ZB(I) = ZE(1,EL_N,IRK+1)
         DB(I) = HB(1,EL_N,1)
         QXB(I) = QX(1,EL_N,IRK+1)
         QYB(I) = QY(1,EL_N,IRK+1)
C     ENDIF

C.......Edge normal vector
         GED = NELED(I,L)
C.......Shouldn't be normal vector (EJK) see Cockburn and Shu reference
C     NXX(I) = COSNX(GED)
C     NYY(I) = SINNX(GED)
         NXX(I) = XM(I)-XB(I)
         IF (NXX(I).NE.0) NXX(I) = NXX(I)/ABS(NXX(I))
         NYY(I) = YM(I)-YB(I)
         IF (NYY(I).NE.0) NYY(I) = NYY(I)/ABS(NYY(I))
c     if (myproc.eq.42.and.l.eq.4229) then
c     write(42,*) 'slope limiter ',i,el_n,qxm(i),qxb(i)
c     endif
      ENDDO


C.....Compute alpha1 and alpha2      

      DO I = 1,3
         EL_N1 = I
         EL_N2 = MOD(I,3)+1
         A(1,1) = XB(EL_N1) - XB(0)
         A(1,2) = XB(EL_N2) - XB(0)
         A(2,1) = YB(EL_N1) - YB(0)
         A(2,2) = YB(EL_N2) - YB(0)
         B(1) = XM(I) - XB(0)
         B(2) = YM(I) - YB(0)
         DTM = A(1,1)*A(2,2) - A(1,2)*A(2,1) ! Compute the determinant
         
C.......Compute the geometric factors alpha1 and alpha2. If either one
C.......is negative that element must be set to p = 0.

         ALPHA1(I) = ( A(2,2)*B(1) - A(1,2)*B(2))/DTM
         IF (ALPHA1(I).LT.0) THEN
            ZE(2:DOF,L,IRK+1) = 0.D0
            QX(2:DOF,L,IRK+1) = 0.D0
            QY(2:DOF,L,IRK+1) = 0.D0
            RETURN
         ENDIF
         ALPHA2(I) = (-A(2,1)*B(1) + A(1,1)*B(2))/DTM
         IF (ALPHA2(I).LT.0) THEN
            ZE(2:DOF,L,IRK+1) = 0.D0
            QX(2:DOF,L,IRK+1) = 0.D0
            QY(2:DOF,L,IRK+1) = 0.D0
            RETURN
         ENDIF
      ENDDO
      
C.....Set the variable vectors at the baricenter of element 0

      WB_O(1,0) = ZB(0)
      WB_O(2,0) = QXB(0)
      WB_O(3,0) = QYB(0)

C.....Start computing deltas for each edge

      DO I = 1,3

C.......Compute eigen values and vectors at the midpoint using the variables at the barycenter

         HTB = ZB(0)+DB(0)
         CALL HYDRO_EIGEN_VALUES(HTB, QXB(0)/HTB, QYB(0)/HTB,
     &        NXX(I), NYY(I), EIGVAL, RI, LE)

C.......Set variable vectors

         WM_O(1) = ZM(I)
         WM_O(2) = QXM(I)
         WM_O(3) = QYM(I)
         EL_N1 = I
         EL_N2 = MOD(I,3)+1
         WB_O(1,1) = ZB(EL_N1)
         WB_O(2,1) = QXB(EL_N1)
         WB_O(3,1) = QYB(EL_N1)
         
         WB_O(1,2) = ZB(EL_N2)
         WB_O(2,2) = QXB(EL_N2)
         WB_O(3,2) = QYB(EL_N2)

C.......Transform original variables into the characteristic space

         DO J = 0,2
            WB_C(1,J)
     &       = LE(1,1)*WB_O(1,J) + LE(1,2)*WB_O(2,J) + LE(1,3)*WB_O(3,J)
            WB_C(2,J)
     &       = LE(2,1)*WB_O(1,J) + LE(2,2)*WB_O(2,J) + LE(2,3)*WB_O(3,J)
            WB_C(3,J)
     &       = LE(3,1)*WB_O(1,J) + LE(3,2)*WB_O(2,J) + LE(3,3)*WB_O(3,J)
         ENDDO
         WM_C(1)= LE(1,1)*WM_O(1) + LE(1,2)*WM_O(2) + LE(1,3)*WM_O(3)
         WM_C(2)= LE(2,1)*WM_O(1) + LE(2,2)*WM_O(2) + LE(2,3)*WM_O(3)
         WM_C(3)= LE(3,1)*WM_O(1) + LE(3,2)*WM_O(2) + LE(3,3)*WM_O(3)

C.......Compute W_TILDA
         W_TILDA(1) = WM_C(1) - WB_C(1,0)
         W_TILDA(2) = WM_C(2) - WB_C(2,0)
         W_TILDA(3) = WM_C(3) - WB_C(3,0)
         
C.......Compute DELTA_W
         DELTA_W(1)
     &        = ALPHA1(I)*(WB_C(1,1) - WB_C(1,0))
     &        + ALPHA2(I)*(WB_C(1,2) - WB_C(1,0))
         DELTA_W(2)
     &        = ALPHA1(I)*(WB_C(2,1) - WB_C(2,0))
     &        + ALPHA2(I)*(WB_C(2,2) - WB_C(2,0))
         DELTA_W(3)
     &        = ALPHA1(I)*(WB_C(3,1) - WB_C(3,0))
     &        + ALPHA2(I)*(WB_C(3,2) - WB_C(3,0))

C.......Apply the TVB modified minmod function
C.......Switched to minmod function (EJK)

C     R_SQ = (XM(I) - XB(0))*(XM(I) - XB(0))
C     &       + (YM(I) - YB(0))*(YM(I) - YB(0))

         NOT_LIMITED = 0
c     if (myproc.eq.42.and.l.eq.4229) then
c     write(42,*) 'slope limiter ',delta_w(1),delta_w(2),
c     $          delta_w(3),w_tilda(1),w_tilda(2),w_tilda(3)
c     endif
         DO K = 1,3             ! Loop for components of variable vectors
C     IF(ABS(W_TILDA(K)) < SL2_M*HLIMIT*R_SQ) THEN
C     DELTA_C(K,I) = W_TILDA(K)
C     ELSE
            IF((W_TILDA(K)*DELTA_W(K) >= 0.D0)
     &           .AND.(ABS(W_TILDA(K)) > ZERO)) THEN
               S = W_TILDA(K)/ABS(W_TILDA(K))
C.......Hard-wired SL2_NYU = 1.5 as taken in Cockburn and Shu reference
               DELTA_C(K,I) = S*MIN(
     &              ABS(W_TILDA(K)),ABS(1.5*DELTA_W(K)))
               IF (DELTA_C(K,I).EQ.(S*ABS(W_TILDA(K))))
     &              NOT_LIMITED(K) = NOT_LIMITED(K) + 1
            ELSE
               DELTA_C(K,I) = 0.D0
            ENDIF
C     ENDIF
         ENDDO
      ENDDO                     ! Loop for edges ends here
      
C.....If none of the variables were limited then end
      
      IF ((NOT_LIMITED(1).EQ.3).AND.(NOT_LIMITED(2).EQ.3)
     &     .AND.(NOT_LIMITED(3).EQ.3)) THEN
         RETURN
         
C.....Else perform limiting on characteristic variables

      ELSE
         DO I = 1,3             ! Loop for variable components
            POS = 0.D0
            NEG = 0.D0
            DO J = 1,3          ! Loop for edges
               POS = POS + MAX(0.D0, DELTA_C(I,J))
               NEG = NEG + MAX(0.D0,-DELTA_C(I,J))
            ENDDO
            IF((POS < ZERO).OR.(NEG < ZERO)) THEN
               DELTA_HAT(I,1) = 0.D0
               DELTA_HAT(I,2) = 0.D0
               DELTA_HAT(I,3) = 0.D0
            ELSE
               THETA_P = MIN(1.D0,NEG/POS)
               THETA_N = MIN(1.D0,POS/NEG)
               DELTA_HAT(I,1) = THETA_P*MAX(0.D0, DELTA_C(I,1))
     &              - THETA_N*MAX(0.D0,-DELTA_C(I,1))
               DELTA_HAT(I,2) = THETA_P*MAX(0.D0, DELTA_C(I,2))
     &              - THETA_N*MAX(0.D0,-DELTA_C(I,2))
               DELTA_HAT(I,3) = THETA_P*MAX(0.D0, DELTA_C(I,3))
     &              - THETA_N*MAX(0.D0,-DELTA_C(I,3))
            ENDIF

C.........Computed limited characteristic variables at the corner nodes

            WC_C(I,1) = WB_C(I,0) - DELTA_HAT(I,1) + DELTA_HAT(I,2)
     &           + DELTA_HAT(I,3)
            WC_C(I,2) = WB_C(I,0) + DELTA_HAT(I,1) - DELTA_HAT(I,2)
     &           + DELTA_HAT(I,3)
            WC_C(I,3) = WB_C(I,0) + DELTA_HAT(I,1) + DELTA_HAT(I,2)
     &           - DELTA_HAT(I,3)
         enddo
C.........Transform back to the original space
         do i=1,3               !loop over corners
            WC_O(1,I) = RI(1,1)*WC_C(1,I) + RI(1,2)*WC_C(2,I)
     &           + RI(1,3)*WC_C(3,I)
            WC_O(2,I) = RI(2,1)*WC_C(1,I) + RI(2,2)*WC_C(2,I)
     &           + RI(2,3)*WC_C(3,I)
            WC_O(3,I) = RI(3,1)*WC_C(1,I) + RI(3,2)*WC_C(2,I)
     &           + RI(3,3)*WC_C(3,I)

c     if (myproc.eq.42.and.l.eq.4229) then
c     write(42,*) 'slope limiter ',wc_c(i,1),wc_c(i,2),wc_c(i,3)
c     endif
         ENDDO                 

         
C.......If the limited water depth < 0, don't apply this slope limiting

         IF((WC_O(1,1)+DC(1)).LE.0.D0.OR.(WC_O(1,2)+DC(2)).LE.0.D0.OR.
     &        (WC_O(1,2)+DC(2)).LE.0.D0) RETURN
         
C.......Else compute new modal dofs

         ZE(2,L,IRK+1) = -1.D0/6.D0*(WC_O(1,1) + WC_O(1,2))
     &        +  1.D0/3.D0*WC_O(1,3)
         ZE(3,L,IRK+1) = -0.5D0*WC_O(1,1) + 0.5D0*WC_O(1,2)
         QX(2,L,IRK+1) = -1.D0/6.D0*(WC_O(2,1) + WC_O(2,2))
     &        +  1.D0/3.D0*WC_O(2,3)
         QX(3,L,IRK+1) = -0.5D0*WC_O(2,1) + 0.5D0*WC_O(2,2)
         QY(2,L,IRK+1) = -1.D0/6.D0*(WC_O(3,1) + WC_O(3,2))
     &        +  1.D0/3.D0*WC_O(3,3)
         QY(3,L,IRK+1) = -0.5D0*WC_O(3,1) + 0.5D0*WC_O(3,2)

      ENDIF

      RETURN
      END SUBROUTINE

      SUBROUTINE SLOPELIMITER3_2NBORS(L)

C.....Use appropriate modules

      USE SIZES, ONLY : SZ
      USE GLOBAL
      USE DG

      IMPLICIT NONE

C.....Declare local variables

      INTEGER,INTENT(IN) ::  L
      INTEGER GED, II,i,j,kk,k
      INTEGER EL2, EL3, EL_N, EL_N1, EL_N2, EL_NN(2), NNBORS
      
      REAL(SZ) DTM, DZEDR, DZEDS, DQXDR, DQXDS, DQYDR, DQYDS
      REAL(SZ) DZEDX(4), DZEDY(4)
      REAL(SZ) DQXDX(4), DQXDY(4)
      REAL(SZ) DQYDX(4), DQYDY(4)
      REAL(SZ) GRADZE(4), GRADQX(4), GRADQY(4)
      REAL(SZ) SL1, SL2
      REAL(SZ) ZELIM, ZELIMX, ZELIMY
      REAL(SZ) QXLIM, QXLIMX, QXLIMY
      REAL(SZ) QYLIM, QYLIMX, QYLIMY
      REAL(SZ) XB(0:2), YB(0:2)
      REAL(SZ) XM(1:3), YM(1:3)
      REAL(SZ) ZB(0:2), DB(0:2), QXB(0:2), QYB(0:2) ! Values at the barycenters
      REAL(SZ) ZC(1:3), DC(1:3), QXC(1:3), QYC(1:3) ! Values at the corner nodes
      REAL(SZ) ZM(1:3), DM(1:3), QXM(1:3), QYM(1:3) ! Values at the edge midpoints
      REAL(SZ) HTB              ! Total height at the barycenter
      REAL(SZ) NXX(3), NYY(3)
      REAL(SZ) ALPHA1(3), ALPHA2(3)
      REAL(SZ) A(2,2), B(2)
      REAL(SZ) TE(3,3)
      REAL(SZ) WB_O(3,0:2), WM_O(3) ! Variable vectors in original space
      REAL(SZ) WB_C(3,0:2), WM_C(3) ! Variable vectors in characteristic space
      REAL(SZ) WC_O(3,3), WC_C(3,3) ! Variable vectors in original space at the corner nodes
      REAL(SZ) W_TILDA(3), DELTA_W(3)
      REAL(SZ) R_SQ             ! Squared r
      REAL(SZ) HMIN, HMAX, HAVG, HLIMIT
      REAL(SZ) S
      REAL(SZ) DELTA_O(3,3), DELTA_C(3,3), DELTA_HAT(3,3)
      INTEGER N(3), NOT_LIMITED(3)
      REAL(SZ) POS, NEG, THETA_P, THETA_N
      LOGICAL LIMITIT(3)
      REAL(SZ), PARAMETER :: ZERO = 1.D-15


C.....Find two neighboring elements
      NNBORS = 0
      DO I = 1,3
         EL_N = EL_NBORS(I,L)
         IF(EL_N.NE.0) THEN
            NNBORS = NNBORS + 1
            EL_NN(NNBORS) = EL_N
         ENDIF
      ENDDO

C.....Do nothing if a neighboring element is dry.
C     IF(WDFLG(L).EQ.0) RETURN
C     DO I = 1,NNBORS
C     EL_N = EL_NN(I)
C     IF(WDFLG(EL_N).EQ.0) RETURN
C     ENDDO
      
C.....Retrieve the barycenter info of this element
      XB(0) = XBC(L)
      YB(0) = YBC(L)
      ZB(0) = ZE(1,L,IRK+1)
      DB(0) = HB(1,L,1)
      QXB(0) = QX(1,L,IRK+1)
      QYB(0) = QY(1,L,IRK+1)
      
C.....Retrieve the nodal info

      N(1) = NM(L,1)
      N(2) = NM(L,2)
      N(3) = NM(L,3)

      pa = PDG_EL(L)

      if (pa.eq.0) then
         pa = 1
      endif

      DO K = 1,3
         ZC(K) = ZE(1,L,IRK+1)
         DC(K) = HB(1,L,1)
         QXC(K) = QX(1,L,IRK+1)
         QYC(K) = QY(1,L,IRK+1)
         DO KK = 2,3            ! This slope limiter assumes the functions to limit as a linear function.
            ZC(K) = ZC(K) + PHI_CORNER(KK,K,pa)*ZE(KK,L,IRK+1)
            DC(K) = DC(K) + PHI_CORNER(KK,K,pa)*HB(KK,L,1)
            QXC(K) = QXC(K) + PHI_CORNER(KK,K,pa)*QX(KK,L,IRK+1)
            QYC(K) = QYC(K) + PHI_CORNER(KK,K,pa)*QY(KK,L,IRK+1)
         ENDDO
      ENDDO

      HMIN = MAX(0.D0,MIN(ZC(1)+DC(1),ZC(2)+DC(2),ZC(3)+DC(3)))
      HMAX = MAX(ZC(1)+DC(1),ZC(2)+DC(2),ZC(3)+DC(3))
      HAVG = ZE(1,L,IRK+1) + HB(1,L,1)
      HLIMIT = MIN(1.D0,(HMIN*HMIN)/(SL3_MD*SL3_MD))
      HLIMIT = 1.D0
      IF(nolifa.eq.2.and.HMIN.LE.H0) RETURN
c     IF(HMIN.LE.0.01D0) RETURN

C.....Retrieve info at edge midpoints
      DO I = 1,3
C.......Edge midpoints
         N1 = N(MOD(I+0,3)+1)
         N2 = N(MOD(I+1,3)+1)
         XM(I) = 0.5D0*(X(N1) + X(N2))
         YM(I) = 0.5D0*(Y(N1) + Y(N2))

         N1 = MOD(I+0,3)+1
         N2 = MOD(I+1,3)+1
         ZM(I) = 0.5D0*(ZC(N1) + ZC(N2))
         DM(I) = 0.5D0*(DC(N1) + DC(N2))
         QXM(I) = 0.5D0*(QXC(N1) + QXC(N2))
         QYM(I) = 0.5D0*(QYC(N1) + QYC(N2))
         
C.......Edge normal vector
         GED = NELED(I,L)
c     NXX(I) = -COSNX(GED)
c     NYY(I) = -SINNX(GED)
         NXX(I) = XM(I)-XB(I)
         IF (NXX(I).NE.0) NXX(I) = NXX(I)/ABS(NXX(I))
         NYY(I) = YM(I)-YB(I)
         IF (NYY(I).NE.0) NYY(I) = NYY(I)/ABS(NYY(I))
      ENDDO

C.....Retreive barycenter info
      DO I = 1,NNBORS
         EL_N = EL_NN(I)
         XB(I) = XBC(EL_N)
         YB(I) = YBC(EL_N)
         ZB(I) = ZE(1,EL_N,IRK+1)
         DB(I) = HB(1,EL_N,1)
         QXB(I) = QX(1,EL_N,IRK+1)
         QYB(I) = QY(1,EL_N,IRK+1)
      ENDDO

C.....Compute alpha1 and alpha2      
      DO I = 1,3
         EL_N1 = 1
         EL_N2 = 2
         A(1,1) = XB(EL_N1) - XB(0)
         A(1,2) = XB(EL_N2) - XB(0)
         A(2,1) = YB(EL_N1) - YB(0)
         A(2,2) = YB(EL_N2) - YB(0)
         B(1) = XM(I) - XB(0)
         B(2) = YM(I) - YB(0)
         DTM = A(1,1)*A(2,2) - A(1,2)*A(2,1) ! Compute the determinant
         ALPHA1(I) = ( A(2,2)*B(1) - A(1,2)*B(2))/DTM
         IF (ALPHA1(I).LT.0) THEN
            ZE(2:DOF,L,IRK+1) = 0.D0
            QX(2:DOF,L,IRK+1) = 0.D0
            QY(2:DOF,L,IRK+1) = 0.D0
            RETURN
         ENDIF
         ALPHA2(I) = (-A(2,1)*B(1) + A(1,1)*B(2))/DTM
         IF (ALPHA2(I).LT.0) THEN
            ZE(2:DOF,L,IRK+1) = 0.D0
            QX(2:DOF,L,IRK+1) = 0.D0
            QY(2:DOF,L,IRK+1) = 0.D0
            RETURN
         ENDIF
      ENDDO

C.....Set the variable vectors at the baricenter of element 0
      WB_O(1,0) = ZB(0)
      WB_O(2,0) = QXB(0)
      WB_O(3,0) = QYB(0)

C.....Start computing deltas for each edge
      DO I = 1,3

C.......Compute eigen values and vectors at the midpoint using the variables at the barycenter
         HTB = ZB(0)+DB(0)
         CALL HYDRO_EIGEN_VALUES(HTB, QXB(0)/HTB, QYB(0)/HTB,
     &        NXX(I), NYY(I), EIGVAL, RI, LE)

C.......Set variable vectors
         WM_O(1) = ZM(I)
         WM_O(2) = QXM(I)
         WM_O(3) = QYM(I)

         EL_N1 = 1
         EL_N2 = 2

         WB_O(1,1) = ZB(EL_N1)
         WB_O(2,1) = QXB(EL_N1)
         WB_O(3,1) = QYB(EL_N1)
         
         WB_O(1,2) = ZB(EL_N2)
         WB_O(2,2) = QXB(EL_N2)
         WB_O(3,2) = QYB(EL_N2)

C.......Transform original variables into the characteristic space
         DO J = 0,2
            WB_C(1,J)
     &       = LE(1,1)*WB_O(1,J) + LE(1,2)*WB_O(2,J) + LE(1,3)*WB_O(3,J)
            WB_C(2,J)
     &       = LE(2,1)*WB_O(1,J) + LE(2,2)*WB_O(2,J) + LE(2,3)*WB_O(3,J)
            WB_C(3,J)
     &       = LE(3,1)*WB_O(1,J) + LE(3,2)*WB_O(2,J) + LE(3,3)*WB_O(3,J)
         ENDDO
         WM_C(1)= LE(1,1)*WM_O(1) + LE(1,2)*WM_O(2) + LE(1,3)*WM_O(3)
         WM_C(2)= LE(2,1)*WM_O(1) + LE(2,2)*WM_O(2) + LE(2,3)*WM_O(3)
         WM_C(3)= LE(3,1)*WM_O(1) + LE(3,2)*WM_O(2) + LE(3,3)*WM_O(3)

C.......Compute W_TILDA
         W_TILDA(1) = WM_C(1) - WB_C(1,0)
         W_TILDA(2) = WM_C(2) - WB_C(2,0)
         W_TILDA(3) = WM_C(3) - WB_C(3,0)
         
C.......Compute DELTA_W
         DELTA_W(1)
     &        = ALPHA1(I)*(WB_C(1,1) - WB_C(1,0))
     &        + ALPHA2(I)*(WB_C(1,2) - WB_C(1,0))
         DELTA_W(2)
     &        = ALPHA1(I)*(WB_C(2,1) - WB_C(2,0))
     &        + ALPHA2(I)*(WB_C(2,2) - WB_C(2,0))
         DELTA_W(3)
     &        = ALPHA1(I)*(WB_C(3,1) - WB_C(3,0))
     &        + ALPHA2(I)*(WB_C(3,2) - WB_C(3,0))

C.......Apply the TVB modified minmod function
c     R_SQ = (XM(I) - XB(0))*(XM(I) - XB(0))
c     &       + (YM(I) - YB(0))*(YM(I) - YB(0))
         NOT_LIMITED=0

         DO K = 1,3             ! Loop for components of variable vectors
c     IF(ABS(W_TILDA(K)) < SL2_M*HLIMIT*R_SQ) THEN
c     DELTA_C(K,I) = W_TILDA(K)
c     ELSE
            IF((W_TILDA(K)*DELTA_W(K) >= 0.D0)
     &           .AND.(ABS(W_TILDA(K)) > ZERO)) THEN
               S = W_TILDA(K)/ABS(W_TILDA(K))
               DELTA_C(K,I) = S*MIN(
     &              ABS(W_TILDA(K)),ABS(SL2_NYU*DELTA_W(K)))
               IF (DELTA_C(K,I).EQ.(S*ABS(W_TILDA(K))))
     &              NOT_LIMITED(K) = NOT_LIMITED(K) + 1
            ELSE
               DELTA_C(K,I) = 0.D0
            ENDIF
c     ENDIF
         ENDDO
      ENDDO
C.....If none of the variables were limited then end
      
      IF ((NOT_LIMITED(1).EQ.3).AND.(NOT_LIMITED(2).EQ.3)
     &     .AND.(NOT_LIMITED(3).EQ.3)) THEN
         RETURN
         
C.....Else perform limiting on characteristic variables

      ELSE
         DO I = 1,3             ! Loop for variable components
            POS = 0.D0
            NEG = 0.D0
            DO J = 1,3          ! Loop for edges
               POS = POS + MAX(0.D0, DELTA_C(I,J))
               NEG = NEG + MAX(0.D0,-DELTA_C(I,J))
            ENDDO
            IF((POS < ZERO).OR.(NEG < ZERO)) THEN
               DELTA_HAT(I,1) = 0.D0
               DELTA_HAT(I,2) = 0.D0
               DELTA_HAT(I,3) = 0.D0
            ELSE
               THETA_P = MIN(1.D0,NEG/POS)
               THETA_N = MIN(1.D0,POS/NEG)
               DELTA_HAT(I,1) = THETA_P*MAX(0.D0, DELTA_C(I,1))
     &              - THETA_N*MAX(0.D0,-DELTA_C(I,1))
               DELTA_HAT(I,2) = THETA_P*MAX(0.D0, DELTA_C(I,2))
     &              - THETA_N*MAX(0.D0,-DELTA_C(I,2))
               DELTA_HAT(I,3) = THETA_P*MAX(0.D0, DELTA_C(I,3))
     &              - THETA_N*MAX(0.D0,-DELTA_C(I,3))
            ENDIF

C.........Computed limited characteristic variables at the corner nodes

            WC_C(I,1) = WB_C(I,0) - DELTA_HAT(I,1) + DELTA_HAT(I,2)
     &           + DELTA_HAT(I,3)
            WC_C(I,2) = WB_C(I,0) + DELTA_HAT(I,1) - DELTA_HAT(I,2)
     &           + DELTA_HAT(I,3)
            WC_C(I,3) = WB_C(I,0) + DELTA_HAT(I,1) + DELTA_HAT(I,2)
     &           - DELTA_HAT(I,3)
         enddo                  ! End variable component loop
C.........Transform back to the original space
         do i=1,3

            WC_O(1,I) = RI(1,1)*WC_C(1,I) + RI(1,2)*WC_C(2,I)
     &           + RI(1,3)*WC_C(3,I)
            WC_O(2,I) = RI(2,1)*WC_C(1,I) + RI(2,2)*WC_C(2,I)
     &           + RI(2,3)*WC_C(3,I)
            WC_O(3,I) = RI(3,1)*WC_C(1,I) + RI(3,2)*WC_C(2,I)
     &           + RI(3,3)*WC_C(3,I)

         ENDDO 
         
C.......If the limited water depth < 0, don't apply this slope limiting

         IF((WC_O(1,1)+DC(1)).LE.0.D0.OR.(WC_O(1,2)+DC(2)).LE.0.D0.OR.
     &        (WC_O(1,2)+DC(2)).LE.0.D0) RETURN
         
C.......Else compute new modal dofs

         ZE(2,L,IRK+1) = -1.D0/6.D0*(WC_O(1,1) + WC_O(1,2))
     &        +  1.D0/3.D0*WC_O(1,3)
         ZE(3,L,IRK+1) = -0.5D0*WC_O(1,1) + 0.5D0*WC_O(1,2)
         QX(2,L,IRK+1) = -1.D0/6.D0*(WC_O(2,1) + WC_O(2,2))
     &        +  1.D0/3.D0*WC_O(2,3)
         QX(3,L,IRK+1) = -0.5D0*WC_O(2,1) + 0.5D0*WC_O(2,2)
         QY(2,L,IRK+1) = -1.D0/6.D0*(WC_O(3,1) + WC_O(3,2))
     &        +  1.D0/3.D0*WC_O(3,3)
         QY(3,L,IRK+1) = -0.5D0*WC_O(3,1) + 0.5D0*WC_O(3,2)

      ENDIF



      RETURN
      END SUBROUTINE

      SUBROUTINE SLOPELIMITER3_1NBORS(L)

C.....Use appropriate modules

      USE SIZES, ONLY : SZ
      USE GLOBAL
      USE DG

      IMPLICIT NONE

C.....Declare local variables

      INTEGER,INTENT(IN) ::  L
      INTEGER GED, II,i,j,kk,k
      INTEGER EL2, EL3, EL_N, EL_N1, EL_N2
      
      REAL(SZ) DTM, DZEDR, DZEDS, DQXDR, DQXDS, DQYDR, DQYDS
      REAL(SZ) DZEDX(4), DZEDY(4)
      REAL(SZ) DQXDX(4), DQXDY(4)
      REAL(SZ) DQYDX(4), DQYDY(4)
      REAL(SZ) GRADZE(4), GRADQX(4), GRADQY(4)
      REAL(SZ) SL1, SL2
      REAL(SZ) ZELIM, ZELIMX, ZELIMY
      REAL(SZ) QXLIM, QXLIMX, QXLIMY
      REAL(SZ) QYLIM, QYLIMX, QYLIMY
      REAL(SZ) XB(0:3), YB(0:3)
      REAL(SZ) XM(1:3), YM(1:3)
      REAL(SZ) ZB(0:3), DB(0:3), QXB(0:3), QYB(0:3) ! Values at the barycenters
      REAL(SZ) ZC(1:3), DC(1:3), QXC(1:3), QYC(1:3) ! Values at the corner nodes
      REAL(SZ) ZM(1:3), DM(1:3), QXM(1:3), QYM(1:3) ! Values at the edge midpoints
      REAL(SZ) HTB              ! Total height at the barycenter
      REAL(SZ) NXX(3), NYY(3)
      REAL(SZ) ALPHA1(3), ALPHA2(3)
      REAL(SZ) A(2,2), B(2)
      REAL(SZ) TE(3,3)
      REAL(SZ) WB_O(3,0:2), WM_O(3) ! Variable vectors in original space
      REAL(SZ) WB_C(3,0:2), WM_C(3) ! Variable vectors in characteristic space
      REAL(SZ) WC_O(3,3)        ! Variable vectors in original space at the corner nodes
      REAL(SZ) W_TILDA(3), DELTA_W(3)
      REAL(SZ) R_SQ             ! Squared r
      REAL(SZ) HMIN, HMAX, HAVG, HLIMIT
      REAL(SZ) S
      REAL(SZ) DELTA_O(3,3), DELTA_C(3,3), DELTA_HAT(3,3)
      INTEGER N(3)
      REAL(SZ) POS, NEG, THETA_P, THETA_N
      LOGICAL LIMITIT(3)
      REAL(SZ), PARAMETER :: ZERO = 1.D-15

C.....Retrieve the barycenter info of this element
      XB(0) = XBC(L)
      YB(0) = YBC(L)
      ZB(0) = ZE(1,L,IRK+1)
      DB(0) = HB(1,L,1)
      QXB(0) = QX(1,L,IRK+1)
      QYB(0) = QY(1,L,IRK+1)
      
C.....Retrieve the nodal info
      limitit(1)=.true.
      limitit(2)=.true.
      limitit(3)=.true.
      IF(LIMITIT(1)) THEN
         ZE(2,L,IRK+1) = 0.0
         ZE(3,L,IRK+1) = 0.0
      ENDIF
      IF(LIMITIT(2)) THEN
         QX(2,L,IRK+1) = 0.0
         QX(3,L,IRK+1) = 0.0
      ENDIF
      IF(LIMITIT(3)) THEN
         QY(2,L,IRK+1) = 0.0
         QY(3,L,IRK+1) = 0.0
         return
      ENDIF

      N(1) = NM(L,1)
      N(2) = NM(L,2)
      N(3) = NM(L,3)

      pa = PDG_EL(L)

      if (pa.eq.0) then
         pa = 1
      endif

      DO K = 1,3
         ZC(K) = ZE(1,L,IRK+1)
         DC(K) = HB(1,L,1)
         QXC(K) = QX(1,L,IRK+1)
         QYC(K) = QY(1,L,IRK+1)
         DO KK = 2,3            ! This slope limiter assumes the functions to limit as a linear function.
            ZC(K) = ZC(K) + PHI_CORNER(KK,K,pa)*ZE(KK,L,IRK+1)
            DC(K) = DC(K) + PHI_CORNER(KK,K,pa)*HB(KK,L,1)
            QXC(K) = QXC(K) + PHI_CORNER(KK,K,pa)*QX(KK,L,IRK+1)
            QYC(K) = QYC(K) + PHI_CORNER(KK,K,pa)*QY(KK,L,IRK+1)
         ENDDO
      ENDDO

      HMIN = MAX(0.D0,MIN(ZC(1)+DC(1),ZC(2)+DC(2),ZC(3)+DC(3)))
      HMAX = MAX(ZC(1)+DC(1),ZC(2)+DC(2),ZC(3)+DC(3))
      HAVG = ZE(1,L,IRK+1) + HB(1,L,1)
      HLIMIT = MIN(1.D0,(HMIN*HMIN)/(SL3_MD*SL3_MD))
      HLIMIT = 1.D0
      IF(NOLIFA.eq.2.and.HMIN.LE.H0) RETURN

C.....Retrieve info at barycenters and edge midpoints
      DO I = 1,3
         EL_N = EL_NBORS(I,L)

C.......Edge midpoints
         N1 = N(MOD(I+0,3)+1)
         N2 = N(MOD(I+1,3)+1)
         XM(I) = 0.5D0*(X(N1) + X(N2))
         YM(I) = 0.5D0*(Y(N1) + Y(N2))

         N1 = MOD(I+0,3)+1
         N2 = MOD(I+1,3)+1
         ZM(I) = 0.5D0*(ZC(N1) + ZC(N2))
         DM(I) = 0.5D0*(DC(N1) + DC(N2))
         QXM(I) = 0.5D0*(QXC(N1) + QXC(N2))
         QYM(I) = 0.5D0*(QYC(N1) + QYC(N2))
         
C.......Barycenters
         IF(EL_N == 0) THEN
C     If there is no neighboring element on this side (e.g. domain boundary),
C     locate a barycenter at the mirroring point and use extrapolated values.
            XB(I) = 2.D0*XM(I) - XB(0)
            YB(I) = 2.D0*YM(I) - YB(0)
            ZB(I) = 2.D0*ZM(I) - ZB(0)
            DB(I) = 2.D0*DM(I) - DB(0)
            QXB(I) = 2.D0*QXM(I) - QXB(0)
            QYB(I) = 2.D0*QYM(I) - QYB(0)
         ELSE
            XB(I) = XBC(EL_N)
            YB(I) = YBC(EL_N)
            ZB(I) = ZE(1,EL_N,IRK+1)
            DB(I) = HB(1,EL_N,1)
            QXB(I) = QX(1,EL_N,IRK+1)
            QYB(I) = QY(1,EL_N,IRK+1)
         ENDIF

C.......Edge normal vector
         GED = NELED(I,L)
         NXX(I) = COSNX(GED)
         NYY(I) = SINNX(GED)
      ENDDO

C.....Compute alpha1 and alpha2      
      DO I = 1,3
         EL_N1 = I
         EL_N2 = MOD(I,3)+1
         A(1,1) = XB(EL_N1) - XB(0)
         A(1,2) = XB(EL_N2) - XB(0)
         A(2,1) = YB(EL_N1) - YB(0)
         A(2,2) = YB(EL_N2) - YB(0)
         B(1) = XM(I) - XB(0)
         B(2) = YM(I) - YB(0)
         DTM = A(1,1)*A(2,2) - A(1,2)*A(2,1) ! Compute the determinant
         ALPHA1(I) = ( A(2,2)*B(1) - A(1,2)*B(2))/DTM
         ALPHA2(I) = (-A(2,1)*B(1) + A(1,1)*B(2))/DTM
      ENDDO


C.....Set the variable vectors at the baricenter of element 0
      WB_O(1,0) = ZB(0)
      WB_O(2,0) = QXB(0)
      WB_O(3,0) = QYB(0)

C.....Start computing deltas for each edge
      DO I = 1,3

C.......Compute eigen values and vectors at the midpoint using the variables at the barycenter
         HTB = ZB(0)+DB(0)
         CALL HYDRO_EIGEN_VALUES(HTB, QXB(0)/HTB, QYB(0)/HTB,
     &        NXX(I), NYY(I), EIGVAL, RI, LE)

C.......Set variable vectors
         WM_O(1) = ZM(I)
         WM_O(2) = QXM(I)
         WM_O(3) = QYM(I)

         EL_N1 = I
         EL_N2 = MOD(I,3)+1

         WB_O(1,1) = ZB(EL_N1)
         WB_O(2,1) = QXB(EL_N1)
         WB_O(3,1) = QYB(EL_N1)
         
         WB_O(1,2) = ZB(EL_N2)
         WB_O(2,2) = QXB(EL_N2)
         WB_O(3,2) = QYB(EL_N2)

C.......Transform original variables into the characteristic space
         DO J = 0,2
            WB_C(1,J)
     &       = LE(1,1)*WB_O(1,J) + LE(1,2)*WB_O(2,J) + LE(1,3)*WB_O(3,J)
            WB_C(2,J)
     &       = LE(2,1)*WB_O(1,J) + LE(2,2)*WB_O(2,J) + LE(2,3)*WB_O(3,J)
            WB_C(3,J)
     &       = LE(3,1)*WB_O(1,J) + LE(3,2)*WB_O(2,J) + LE(3,3)*WB_O(3,J)
         ENDDO
         WM_C(1)= LE(1,1)*WM_O(1) + LE(1,2)*WM_O(2) + LE(1,3)*WM_O(3)
         WM_C(2)= LE(2,1)*WM_O(1) + LE(2,2)*WM_O(2) + LE(2,3)*WM_O(3)
         WM_C(3)= LE(3,1)*WM_O(1) + LE(3,2)*WM_O(2) + LE(3,3)*WM_O(3)

C.......Compute W_TILDA
         W_TILDA(1) = WM_C(1) - WB_C(1,0)
         W_TILDA(2) = WM_C(2) - WB_C(2,0)
         W_TILDA(3) = WM_C(3) - WB_C(3,0)
         
C.......Compute DELTA_W
         DELTA_W(1)
     &        = ALPHA1(I)*(WB_C(1,1) - WB_C(1,0))
     &        + ALPHA2(I)*(WB_C(1,2) - WB_C(1,0))
         DELTA_W(2)
     &        = ALPHA1(I)*(WB_C(2,1) - WB_C(2,0))
     &        + ALPHA2(I)*(WB_C(2,2) - WB_C(2,0))
         DELTA_W(3)
     &        = ALPHA1(I)*(WB_C(3,1) - WB_C(3,0))
     &        + ALPHA2(I)*(WB_C(3,2) - WB_C(3,0))

C.......Apply the TVB modified minmod function
         R_SQ = (XM(I) - XB(0))*(XM(I) - XB(0))
     &        + (YM(I) - YB(0))*(YM(I) - YB(0))

         DO K = 1,3             ! Loop for components of variable vectors
            IF(ABS(W_TILDA(K)) < SL2_M*HLIMIT*R_SQ) THEN
               DELTA_C(K,I) = W_TILDA(K)
            ELSE
               IF((W_TILDA(K)*DELTA_W(K) >= 0.D0)
     &              .AND.(ABS(W_TILDA(K)) > ZERO)) THEN
                  S = W_TILDA(K)/ABS(W_TILDA(K))
                  DELTA_C(K,I) = S*MIN(
     &                 ABS(W_TILDA(K)),ABS(SL2_NYU*DELTA_W(K)))
               ELSE
                  DELTA_C(K,I) = 0.D0
               ENDIF
            ENDIF
         ENDDO

C.......Transform back to the original space
         DELTA_O(1,I)
     &        = RI(1,1)*DELTA_C(1,I)
     &        + RI(1,2)*DELTA_C(2,I)
     &        + RI(1,3)*DELTA_C(3,I)
         DELTA_O(2,I)
     &        = RI(2,1)*DELTA_C(1,I)
     &        + RI(2,2)*DELTA_C(2,I)
     &        + RI(2,3)*DELTA_C(3,I)
         DELTA_O(3,I)
     &        = RI(3,1)*DELTA_C(1,I)
     &        + RI(3,2)*DELTA_C(2,I)
     &        + RI(3,3)*DELTA_C(3,I)
         

      ENDDO                     ! Loop for edges ends here

C     Set the variable vector at the midpoint 1. This is used later to judge
C     whether the variables should be limited.
      WM_O(1) = ZM(1)
      WM_O(2) = QXM(1)
      WM_O(3) = QYM(1)

C     Find the possibly limited DELTA and variable vectors at the corner nodes
      DO I = 1,3                ! Loop for variable components
C     Find the possibly limited DELTA
         POS = 0.D0
         NEG = 0.D0
         DO J = 1,3             ! Loop for edges
            POS = POS + MAX(0.D0, DELTA_O(I,J))
            NEG = NEG + MAX(0.D0,-DELTA_O(I,J))
         ENDDO
         IF((POS < ZERO).OR.(NEG < ZERO)) THEN 
            DELTA_HAT(I,1) = 0.D0
            DELTA_HAT(I,2) = 0.D0
            DELTA_HAT(I,3) = 0.D0
         ELSE
            THETA_P = MIN(1.D0,NEG/POS)
            THETA_N = MIN(1.D0,POS/NEG)
            
            DELTA_HAT(I,1)
     &           = THETA_P*MAX(0.D0, DELTA_O(I,1))
     &           - THETA_N*MAX(0.D0,-DELTA_O(I,1))
            DELTA_HAT(I,2)
     &           = THETA_P*MAX(0.D0, DELTA_O(I,2))
     &           - THETA_N*MAX(0.D0,-DELTA_O(I,2))
            DELTA_HAT(I,3)
     &           = THETA_P*MAX(0.D0, DELTA_O(I,3))
     &           - THETA_N*MAX(0.D0,-DELTA_O(I,3))
         ENDIF

C     Compute the possibly limited variable vectors at the corner nodes
         IF(( DELTA_HAT(I,1) < (WM_O(I) - WB_O(I,0) - ZERO) ).OR.
     &        ( DELTA_HAT(I,1) > (WM_O(I) - WB_O(I,0) + ZERO) )) THEN
C     Update the variables only if they need to be limited.
C     If the condition in the if statement is false, the following 
C     isn't applied and therefore the quradratic and higher order 
C     coefficients are preserved.

C     Computing limited variables at the corner nodes
            WC_O(I,1) = WB_O(I,0)
     &           - DELTA_HAT(I,1)
     &           + DELTA_HAT(I,2)
     &           + DELTA_HAT(I,3)
            WC_O(I,2) = WB_O(I,0)
     &           + DELTA_HAT(I,1)
     &           - DELTA_HAT(I,2)
     &           + DELTA_HAT(I,3)
            WC_O(I,3) = WB_O(I,0)
     &           + DELTA_HAT(I,1)
     &           + DELTA_HAT(I,2)
     &           - DELTA_HAT(I,3)
            LIMITIT(I) = .TRUE.
         ELSE
            LIMITIT(I) = .FALSE.
         ENDIF
      ENDDO

C     If the limited water depth is less 0, don't apply this slope limiting.
      IF((WC_O(1,1)+DC(1)).LE.0.D0.OR.
     &     (WC_O(1,2)+DC(2)).LE.0.D0.OR.
     &     (WC_O(1,2)+DC(2)).LE.0.D0) THEN
         GOTO 10000
      ENDIF

      IF(LIMITIT(1)) THEN
         ZE(2,L,IRK+1) = -1.D0/6.D0*(WC_O(1,1) + WC_O(1,2))
     &        + 1.D0/3.D0*WC_O(1,3)
         ZE(3,L,IRK+1) = -0.5D0*WC_O(1,1) + 0.5D0*WC_O(1,2)
      ENDIF
      IF(LIMITIT(2)) THEN
         QX(2,L,IRK+1) = -1.D0/6.D0*(WC_O(2,1) + WC_O(2,2))
     &        + 1.D0/3.D0*WC_O(2,3)
         QX(3,L,IRK+1) = -0.5D0*WC_O(2,1) + 0.5D0*WC_O(2,2)
      ENDIF
      IF(LIMITIT(3)) THEN
         QY(2,L,IRK+1) = -1.D0/6.D0*(WC_O(3,1) + WC_O(3,2))
     &        + 1.D0/3.D0*WC_O(3,3)
         QY(3,L,IRK+1) = -0.5D0*WC_O(3,1) + 0.5D0*WC_O(3,2)
      ENDIF

10000 CONTINUE

      RETURN
      END SUBROUTINE

c$$$C***********************************************************************
c$$$C     
c$$$C     SUBROUTINE SLOPELIMITER3_OLD()
c$$$C     
c$$$C     Written by Shintaro Bunya - 27 Feb 2008
c$$$C     
c$$$C     This subroutine is an extentioin of the slope limiter used in
c$$$C     B. Cockburn and C-W. Shu, "The Runge-Kutta Discontinuous Galerkin
c$$$C     Method for Conservation Laws V, Multidimensional Systems,"
c$$$C     Journal of Computational Physics 141, 199-224 (1998).
c$$$C     
c$$$C     An additional parameter SL3_MD is used in this slope limiter to
c$$$C     make sure that slopes are limited in elements whose water depth
c$$$C     is small. (A special treatment for wetting and drying.
c$$$C     
c$$$C***********************************************************************
c$$$
c$$$      SUBROUTINE SLOPELIMITER3_OLD()
c$$$
c$$$C.....Use appropriate modules
c$$$
c$$$      USE SIZES, ONLY : SZ
c$$$      USE GLOBAL
c$$$      USE DG
c$$$
c$$$      IMPLICIT NONE
c$$$
c$$$C.....Declare local variables
c$$$
c$$$      INTEGER L, GED, II,i,j,kk,k
c$$$      INTEGER EL2, EL3, EL_N, EL_N1, EL_N2
c$$$      
c$$$      REAL(SZ) DTM, DZEDR, DZEDS, DQXDR, DQXDS, DQYDR, DQYDS
c$$$      REAL(SZ) DZEDX(4), DZEDY(4)
c$$$      REAL(SZ) DQXDX(4), DQXDY(4)
c$$$      REAL(SZ) DQYDX(4), DQYDY(4)
c$$$      REAL(SZ) GRADZE(4), GRADQX(4), GRADQY(4)
c$$$      REAL(SZ) SL1, SL2
c$$$      REAL(SZ) ZELIM, ZELIMX, ZELIMY
c$$$      REAL(SZ) QXLIM, QXLIMX, QXLIMY
c$$$      REAL(SZ) QYLIM, QYLIMX, QYLIMY
c$$$      REAL(SZ) XB(0:3), YB(0:3)
c$$$      REAL(SZ) XM(1:3), YM(1:3)
c$$$      REAL(SZ) ZB(0:3), DB(0:3), QXB(0:3), QYB(0:3) ! Values at the barycenters
c$$$      REAL(SZ) ZC(1:3), DC(1:3), QXC(1:3), QYC(1:3) ! Values at the corner nodes
c$$$      REAL(SZ) ZM(1:3), DM(1:3), QXM(1:3), QYM(1:3) ! Values at the edge midpoints
c$$$      REAL(SZ) HTB              ! Total height at the barycenter
c$$$      REAL(SZ) NXX(3), NYY(3)
c$$$      REAL(SZ) ALPHA1(3), ALPHA2(3)
c$$$      REAL(SZ) A(2,2), B(2)
c$$$      REAL(SZ) TE(3,3)
c$$$      REAL(SZ) WB_O(3,0:2), WM_O(3) ! Variable vectors in original space
c$$$      REAL(SZ) WB_C(3,0:2), WM_C(3) ! Variable vectors in characteristic space
c$$$      REAL(SZ) WC_O(3,3)        ! Variable vectors in original space at the corner nodes
c$$$      REAL(SZ) W_TILDA(3), DELTA_W(3)
c$$$      REAL(SZ) R_SQ             ! Squared r
c$$$      REAL(SZ) HMIN, HMAX, HAVG, HLIMIT
c$$$      REAL(SZ) S
c$$$      REAL(SZ) DELTA_O(3,3), DELTA_C(3,3), DELTA_HAT(3,3)
c$$$      INTEGER N(3)
c$$$      REAL(SZ) POS, NEG, THETA_P, THETA_N
c$$$      LOGICAL LIMITIT(3)
c$$$      REAL(SZ), PARAMETER :: ZERO = 1.D-15
c$$$
c$$$C.....Save the original values  02/28/2007 sb
c$$$      DO L=1,NE
c$$$         DO K = 1,DOF
c$$$            ZE(K,L,NRK+2) = ZE(K,L,IRK+1)
c$$$            QX(K,L,NRK+2) = QX(K,L,IRK+1)
c$$$            QY(K,L,NRK+2) = QY(K,L,IRK+1)
c$$$         ENDDO
c$$$      ENDDO      
c$$$
c$$$      DO 1000 L=1, NE
c$$$
c$$$C.....Retrieve the barycenter info of this element
c$$$         XB(0) = XBC(L)
c$$$         YB(0) = YBC(L)
c$$$         ZB(0) = ZE(1,L,IRK+1)
c$$$         DB(0) = HB(1,L,1)
c$$$         QXB(0) = QX(1,L,IRK+1)
c$$$         QYB(0) = QY(1,L,IRK+1)
c$$$         
c$$$C.....Retrieve the nodal info
c$$$
c$$$         N(1) = NM(L,1)
c$$$         N(2) = NM(L,2)
c$$$         N(3) = NM(L,3)
c$$$
c$$$         pa = PDG_EL(L)
c$$$
c$$$         if (pa.eq.0) then
c$$$            pa = 1
c$$$         endif
c$$$
c$$$         DO K = 1,3
c$$$            ZC(K) = ZE(1,L,IRK+1)
c$$$            DC(K) = HB(1,L,1)
c$$$            QXC(K) = QX(1,L,IRK+1)
c$$$            QYC(K) = QY(1,L,IRK+1)
c$$$            DO KK = 2,3         ! This slope limiter assumes the functions to limit as a linear function.
c$$$               ZC(K) = ZC(K) + PHI_CORNER(KK,K,pa)*ZE(KK,L,IRK+1)
c$$$               DC(K) = DC(K) + PHI_CORNER(KK,K,pa)*HB(KK,L,1)
c$$$               QXC(K) = QXC(K) + PHI_CORNER(KK,K,pa)*QX(KK,L,IRK+1)
c$$$               QYC(K) = QYC(K) + PHI_CORNER(KK,K,pa)*QY(KK,L,IRK+1)
c$$$            ENDDO
c$$$         ENDDO
c$$$
c$$$         HMIN = MAX(0.D0,MIN(ZC(1)+DC(1),ZC(2)+DC(2),ZC(3)+DC(3)))
c$$$         HMAX = MAX(ZC(1)+DC(1),ZC(2)+DC(2),ZC(3)+DC(3))
c$$$         HAVG = ZE(1,L,IRK+1) + HB(1,L,1)
c$$$         HLIMIT = MIN(1.D0,(HMIN*HMIN)/(SL3_MD*SL3_MD))
c$$$         HLIMIT = 1.D0
c$$$         IF(HMIN.LE.H0) RETURN
c$$$
c$$$C.....Retrieve info at barycenters and edge midpoints
c$$$         DO I = 1,3
c$$$            EL_N = EL_NBORS(I,L)
c$$$
c$$$C.......Edge midpoints
c$$$            N1 = N(MOD(I+0,3)+1)
c$$$            N2 = N(MOD(I+1,3)+1)
c$$$            XM(I) = 0.5D0*(X(N1) + X(N2))
c$$$            YM(I) = 0.5D0*(Y(N1) + Y(N2))
c$$$
c$$$            N1 = MOD(I+0,3)+1
c$$$            N2 = MOD(I+1,3)+1
c$$$            ZM(I) = 0.5D0*(ZC(N1) + ZC(N2))
c$$$            DM(I) = 0.5D0*(DC(N1) + DC(N2))
c$$$            QXM(I) = 0.5D0*(QXC(N1) + QXC(N2))
c$$$            QYM(I) = 0.5D0*(QYC(N1) + QYC(N2))
c$$$            
c$$$C.......Barycenters
c$$$            IF(EL_N == 0) THEN
c$$$C     If there is no neighboring element on this side (e.g. domain boundary),
c$$$C     locate a barycenter at the mirroring point and use extrapolated values.
c$$$               XB(I) = 2.D0*XM(I) - XB(0)
c$$$               YB(I) = 2.D0*YM(I) - YB(0)
c$$$               ZB(I) = 2.D0*ZM(I) - ZB(0)
c$$$               DB(I) = 2.D0*DM(I) - DB(0)
c$$$               QXB(I) = 2.D0*QXM(I) - QXB(0)
c$$$               QYB(I) = 2.D0*QYM(I) - QYB(0)
c$$$            ELSE
c$$$               XB(I) = XBC(EL_N)
c$$$               YB(I) = YBC(EL_N)
c$$$               ZB(I) = ZE(1,EL_N,IRK+1)
c$$$               DB(I) = HB(1,EL_N,1)
c$$$               QXB(I) = QX(1,EL_N,IRK+1)
c$$$               QYB(I) = QY(1,EL_N,IRK+1)
c$$$            ENDIF
c$$$
c$$$C.......Edge normal vector
c$$$            GED = NELED(I,L)
c$$$            NXX(I) = COSNX(GED)
c$$$            NYY(I) = SINNX(GED)
c$$$         ENDDO
c$$$
c$$$C.....Compute alpha1 and alpha2      
c$$$         DO I = 1,3
c$$$            EL_N1 = I
c$$$            EL_N2 = MOD(I,3)+1
c$$$            A(1,1) = XB(EL_N1) - XB(0)
c$$$            A(1,2) = XB(EL_N2) - XB(0)
c$$$            A(2,1) = YB(EL_N1) - YB(0)
c$$$            A(2,2) = YB(EL_N2) - YB(0)
c$$$            B(1) = XM(I) - XB(0)
c$$$            B(2) = YM(I) - YB(0)
c$$$            DTM = A(1,1)*A(2,2) - A(1,2)*A(2,1) ! Compute the determinant
c$$$            ALPHA1(I) = ( A(2,2)*B(1) - A(1,2)*B(2))/DTM
c$$$            ALPHA2(I) = (-A(2,1)*B(1) + A(1,1)*B(2))/DTM
c$$$         ENDDO
c$$$
c$$$
c$$$C.....Set the variable vectors at the baricenter of element 0
c$$$         WB_O(1,0) = ZB(0)
c$$$         WB_O(2,0) = QXB(0)
c$$$         WB_O(3,0) = QYB(0)
c$$$
c$$$C.....Start computing deltas for each edge
c$$$         DO I = 1,3
c$$$
c$$$C.......Compute eigen values and vectors at the midpoint using the variables at the barycenter
c$$$            HTB = ZB(0)+DB(0)
c$$$            CALL HYDRO_EIGEN_VALUES(HTB, QXB(0)/HTB, QYB(0)/HTB,
c$$$     &           NXX(I), NYY(I), EIGVAL, RI, LE)
c$$$
c$$$C.......Set variable vectors
c$$$            WM_O(1) = ZM(I)
c$$$            WM_O(2) = QXM(I)
c$$$            WM_O(3) = QYM(I)
c$$$
c$$$            EL_N1 = I
c$$$            EL_N2 = MOD(I,3)+1
c$$$
c$$$            WB_O(1,1) = ZB(EL_N1)
c$$$            WB_O(2,1) = QXB(EL_N1)
c$$$            WB_O(3,1) = QYB(EL_N1)
c$$$            
c$$$            WB_O(1,2) = ZB(EL_N2)
c$$$            WB_O(2,2) = QXB(EL_N2)
c$$$            WB_O(3,2) = QYB(EL_N2)
c$$$
c$$$C.......Transform original variables into the characteristic space
c$$$            DO J = 0,2
c$$$               WB_C(1,J)
c$$$     &              = LE(1,1)*WB_O(1,J) + LE(1,2)*WB_O(2,J) + LE(1,3)*WB_O(3,J)
c$$$               WB_C(2,J)
c$$$     &              = LE(2,1)*WB_O(1,J) + LE(2,2)*WB_O(2,J) + LE(2,3)*WB_O(3,J)
c$$$               WB_C(3,J)
c$$$     &              = LE(3,1)*WB_O(1,J) + LE(3,2)*WB_O(2,J) + LE(3,3)*WB_O(3,J)
c$$$            ENDDO
c$$$            WM_C(1)= LE(1,1)*WM_O(1) + LE(1,2)*WM_O(2) + LE(1,3)*WM_O(3)
c$$$            WM_C(2)= LE(2,1)*WM_O(1) + LE(2,2)*WM_O(2) + LE(2,3)*WM_O(3)
c$$$            WM_C(3)= LE(3,1)*WM_O(1) + LE(3,2)*WM_O(2) + LE(3,3)*WM_O(3)
c$$$
c$$$C.......Compute W_TILDA
c$$$            W_TILDA(1) = WM_C(1) - WB_C(1,0)
c$$$            W_TILDA(2) = WM_C(2) - WB_C(2,0)
c$$$            W_TILDA(3) = WM_C(3) - WB_C(3,0)
c$$$            
c$$$C.......Compute DELTA_W
c$$$            DELTA_W(1)
c$$$     &           = ALPHA1(I)*(WB_C(1,1) - WB_C(1,0))
c$$$     &           + ALPHA2(I)*(WB_C(1,2) - WB_C(1,0))
c$$$            DELTA_W(2)
c$$$     &           = ALPHA1(I)*(WB_C(2,1) - WB_C(2,0))
c$$$     &           + ALPHA2(I)*(WB_C(2,2) - WB_C(2,0))
c$$$            DELTA_W(3)
c$$$     &           = ALPHA1(I)*(WB_C(3,1) - WB_C(3,0))
c$$$     &           + ALPHA2(I)*(WB_C(3,2) - WB_C(3,0))
c$$$
c$$$C.......Apply the TVB modified minmod function
c$$$            R_SQ = (XM(I) - XB(0))*(XM(I) - XB(0))
c$$$     &           + (YM(I) - YB(0))*(YM(I) - YB(0))
c$$$
c$$$            DO K = 1,3          ! Loop for components of variable vectors
c$$$               IF(ABS(W_TILDA(K)) < SL2_M*HLIMIT*R_SQ) THEN
c$$$                  DELTA_C(K,I) = W_TILDA(K)
c$$$               ELSE
c$$$                  IF((W_TILDA(K)*DELTA_W(K) >= 0.D0)
c$$$     &                 .AND.(ABS(W_TILDA(K)) > ZERO)) THEN
c$$$                     S = W_TILDA(K)/ABS(W_TILDA(K))
c$$$                     DELTA_C(K,I) = S*MIN(
c$$$     &                    ABS(W_TILDA(K)),ABS(SL2_NYU*DELTA_W(K)))
c$$$                  ELSE
c$$$                     DELTA_C(K,I) = 0.D0
c$$$                  ENDIF
c$$$               ENDIF
c$$$            ENDDO
c$$$
c$$$C.......Transform back to the original space
c$$$            DELTA_O(1,I)
c$$$     &           = RI(1,1)*DELTA_C(1,I)
c$$$     &           + RI(1,2)*DELTA_C(2,I)
c$$$     &           + RI(1,3)*DELTA_C(3,I)
c$$$            DELTA_O(2,I)
c$$$     &           = RI(2,1)*DELTA_C(1,I)
c$$$     &           + RI(2,2)*DELTA_C(2,I)
c$$$     &           + RI(2,3)*DELTA_C(3,I)
c$$$            DELTA_O(3,I)
c$$$     &           = RI(3,1)*DELTA_C(1,I)
c$$$     &           + RI(3,2)*DELTA_C(2,I)
c$$$     &           + RI(3,3)*DELTA_C(3,I)
c$$$            
c$$$
c$$$         ENDDO                  ! Loop for edges ends here
c$$$
c$$$C     Set the variable vector at the midpoint 1. This is used later to judge
c$$$C     whether the variables should be limited.
c$$$         WM_O(1) = ZM(1)
c$$$         WM_O(2) = QXM(1)
c$$$         WM_O(3) = QYM(1)
c$$$
c$$$C     Find the possibly limited DELTA and variable vectors at the corner nodes
c$$$         DO I = 1,3             ! Loop for variable components
c$$$C     Find the possibly limited DELTA
c$$$            POS = 0.D0
c$$$            NEG = 0.D0
c$$$            DO J = 1,3          ! Loop for edges
c$$$               POS = POS + MAX(0.D0, DELTA_O(I,J))
c$$$               NEG = NEG + MAX(0.D0,-DELTA_O(I,J))
c$$$            ENDDO
c$$$            IF((POS < ZERO).OR.(NEG < ZERO)) THEN 
c$$$               DELTA_HAT(I,1) = 0.D0
c$$$               DELTA_HAT(I,2) = 0.D0
c$$$               DELTA_HAT(I,3) = 0.D0
c$$$            ELSE
c$$$               THETA_P = MIN(1.D0,NEG/POS)
c$$$               THETA_N = MIN(1.D0,POS/NEG)
c$$$               
c$$$               DELTA_HAT(I,1)
c$$$     &              = THETA_P*MAX(0.D0, DELTA_O(I,1))
c$$$     &              - THETA_N*MAX(0.D0,-DELTA_O(I,1))
c$$$               DELTA_HAT(I,2)
c$$$     &              = THETA_P*MAX(0.D0, DELTA_O(I,2))
c$$$     &              - THETA_N*MAX(0.D0,-DELTA_O(I,2))
c$$$               DELTA_HAT(I,3)
c$$$     &              = THETA_P*MAX(0.D0, DELTA_O(I,3))
c$$$     &              - THETA_N*MAX(0.D0,-DELTA_O(I,3))
c$$$            ENDIF
c$$$
c$$$C     Compute the possibly limited variable vectors at the corner nodes
c$$$            IF(( DELTA_HAT(I,1) < (WM_O(I) - WB_O(I,0) - ZERO) ).OR.
c$$$     &           ( DELTA_HAT(I,1) > (WM_O(I) - WB_O(I,0) + ZERO) )) THEN
c$$$C     Update the variables only if they need to be limited.
c$$$C     If the condition in the if statement is false, the following 
c$$$C     isn't applied and therefore the quradratic and higher order 
c$$$C     coefficients are preserved.
c$$$
c$$$C     Computing limited variables at the corner nodes
c$$$               WC_O(I,1) = WB_O(I,0)
c$$$     &              - DELTA_HAT(I,1)
c$$$     &              + DELTA_HAT(I,2)
c$$$     &              + DELTA_HAT(I,3)
c$$$               WC_O(I,2) = WB_O(I,0)
c$$$     &              + DELTA_HAT(I,1)
c$$$     &              - DELTA_HAT(I,2)
c$$$     &              + DELTA_HAT(I,3)
c$$$               WC_O(I,3) = WB_O(I,0)
c$$$     &              + DELTA_HAT(I,1)
c$$$     &              + DELTA_HAT(I,2)
c$$$     &              - DELTA_HAT(I,3)
c$$$               LIMITIT(I) = .TRUE.
c$$$            ELSE
c$$$               LIMITIT(I) = .FALSE.
c$$$            ENDIF
c$$$         ENDDO
c$$$
c$$$         IF(LIMITIT(1)) THEN
c$$$            ZE(2,L,IRK+1) = -1.D0/6.D0*(WC_O(1,1) + WC_O(1,2))
c$$$     &           + 1.D0/3.D0*WC_O(1,3)
c$$$            ZE(3,L,IRK+1) = -0.5D0*WC_O(1,1) + 0.5D0*WC_O(1,2)
c$$$         ENDIF
c$$$         IF(LIMITIT(2)) THEN
c$$$            QX(2,L,IRK+1) = -1.D0/6.D0*(WC_O(2,1) + WC_O(2,2))
c$$$     &           + 1.D0/3.D0*WC_O(2,3)
c$$$            QX(3,L,IRK+1) = -0.5D0*WC_O(2,1) + 0.5D0*WC_O(2,2)
c$$$         ENDIF
c$$$         IF(LIMITIT(3)) THEN
c$$$            QY(2,L,IRK+1) = -1.D0/6.D0*(WC_O(3,1) + WC_O(3,2))
c$$$     &           + 1.D0/3.D0*WC_O(3,3)
c$$$            QY(3,L,IRK+1) = -0.5D0*WC_O(3,1) + 0.5D0*WC_O(3,2)
c$$$         ENDIF
c$$$
c$$$ 1000 CONTINUE
c$$$
c$$$      RETURN
c$$$      END SUBROUTINE

***********************************************************************
C     
C     SUBROUTINE HYDRO_EIGEN_VALUES
C     
C     This subroutine computes the eigen values and vectors of hydro
C     dynamics equations.
C     
C     Borrowed from ROE_FLUX() - 06 Feb 2008, S.B.
C     
C***********************************************************************
      SUBROUTINE HYDRO_EIGEN_VALUES(H, U, V, NX, NY, EIGVAL, RI, LE)
      USE SIZES, ONLY : SZ
      USE GLOBAL, ONLY : IFNLCT, G

      IMPLICIT NONE

      REAL(SZ) H,U,V,NX,NY,C,DTM
      REAL(SZ) EIGVAL(3)        ! eigen values
      REAL(SZ) RI(3,3)          ! columns are the right eigen vectors
      REAL(SZ) LE(3,3)          ! columns are the left eigen vectors

      C = SQRT(G*H)

C.....Evaluate the eigenvalues at the Roe averaged variables

      EIGVAL(2) = (U*NX + V*NY)*IFNLCT
      EIGVAL(1) = EIGVAL(2) + C
      EIGVAL(3) = EIGVAL(2) - C
      
C.....Evaluate right eigenvectors at Roe averaged variables

      RI(1,1) = 1.D0
      RI(2,1) = U*IFNLCT + (IFNLCT*C + (1-IFNLCT)*G/C)*NX
      RI(3,1) = V*IFNLCT + (IFNLCT*C + (1-IFNLCT)*G/C)*NY

      RI(1,2) = 0.D0
      RI(2,2) = -NY
      RI(3,2) =  NX

      RI(1,3) = 1.D0
      RI(2,3) = U*IFNLCT - (IFNLCT*C + (1-IFNLCT)*G/C)*NX
      RI(3,3) = V*IFNLCT - (IFNLCT*C + (1-IFNLCT)*G/C)*NY

C.....Evaluate left eigenvectors at Roe averaged variables

C.....Determinant
      DTM = 1.0d0/(RI(2,2)*RI(3,3)-RI(2,3)*RI(3,2)+RI(2,1)*RI(3,2)
     &     -RI(3,1)*RI(2,2))

      LE(1,1) =  DTM*( RI(2,2)*RI(3,3)-RI(2,3)*RI(3,2) )
      LE(1,2) =  DTM*RI(3,2)
      LE(1,3) = -DTM*RI(2,2)

      LE(2,1) = -DTM*( RI(2,1)*RI(3,3)-RI(2,3)*RI(3,1) )
      LE(2,2) =  DTM*( RI(3,3)-RI(3,1) )
      LE(2,3) = -DTM*( RI(2,3)-RI(2,1) )

      LE(3,1) =  DTM*( RI(2,1)*RI(3,2) - RI(2,2)*RI(3,1) )
      LE(3,2) = -DTM*RI(3,2)
      LE(3,3) =  DTM*RI(2,2)

      RETURN
      END SUBROUTINE

c$$$  C invert a 3x3 matrix
c$$$  subroutine invert(a,ainv,succeed)
c$$$  sizes,only : sz
c$$$  implicit none
c$$$  real(sz),intent(in) ::  a(3,3)
c$$$  real(sz),intent(out) ::  ainv(3,3)
c$$$  logical,intent(out) :: succeed
c$$$  integer :: i, j, k
c$$$  real(sz) :: det
c$$$  real(sz) :: e(3,3)
c$$$  
c$$$  det = 
c$$$  &     a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2)) +
c$$$  &     a(1,2)*(a(2,3)*a(3,1)-a(2,1)*a(3,3)) +
c$$$  &     a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))
c$$$  
c$$$  if(det.eq.0.D0) then
c$$$  succeed = .false.
c$$$  return
c$$$  endif
c$$$  
c$$$  Ainv(1,1) = (A(1,1)*A(2,2) - A(1,2)*A(2,1))/det
c$$$  Ainv(1,2) = (A(0,2)*A(2,1) - A(0,1)*A(2,2))/det
c$$$  Ainv(1,3) = (A(0,1)*A(1,2) - A(0,2)*A(1,1))/det
c$$$  
c$$$  Ainv(1,0) = (A(1,2)*A(2,0) - A(1,0)*A(2,2))/det
c$$$  Ainv(1,1) = (A(0,0)*A(2,2) - A(0,2)*A(2,0))/det
c$$$  Ainv(1,2) = (A(0,2)*A(1,0) - A(0,0)*A(1,2))/det
c$$$  
c$$$  Ainv(2,0) = (A(1,0)*A(2,1) - A(1,1)*A(2,0))/det
c$$$  Ainv(2,1) = (A(0,1)*A(2,0) - A(0,0)*A(2,1))/det
c$$$  Ainv(2,2) = (A(0,0)*A(1,1) - A(0,1)*A(1,0))/det
c$$$  
c$$$  do i = 1,3
c$$$  do j = 1,3
c$$$  E(i,j) = 0.0
c$$$  do k = 1,3
c$$$  E(i,j) += A(i,k)*Ainv(k,j)
c$$$  enddo
c$$$  enddo
c$$$  enddo
c$$$  end subroutine
c$$$  
c$$$  /* compute projected velocity u */
c$$$  void compute_projection_u
c$$$  (double *H, double *q, 
c$$$  int num_pts, double *area_coords[3], double *weights,
c$$$  double *u)
c$$$  {
c$$$  double A[3][3], Ainv[3][3];
c$$$  double b[3];
c$$$  
c$$$  int i, j, k, p;
c$$$  
c$$$  /* initialize */
c$$$  for(j = 0; j < 3; j++) {
c$$$  for(i = 0; i < 3; i++) {
c$$$  A[j][i] = 0.0;
c$$$  }
c$$$  b[j] = 0.0;
c$$$  }
c$$$  
c$$$  /* determine matrix A */
c$$$  for(j = 0; j < 3; j++) {
c$$$  for(i = 0; i < 3; i++) {
c$$$  for(k = 0; k < 3; k++) {
c$$$  for(p = 0; p < num_pts; p++) {
c$$$  A[j][i] += 
c$$$  H[k]*
c$$$  area_coords[i][p]*
c$$$  area_coords[j][p]*
c$$$  area_coords[k][p]*
c$$$  weights[p];
c$$$  }
c$$$  }
c$$$  }
c$$$  }
c$$$  
c$$$  /* determine vector b */
c$$$  for(j = 0; j < 3; j++) {
c$$$  for(i = 0; i < 3; i++) {
c$$$  for(p = 0; p < num_pts; p++) {
c$$$  b[j] += 
c$$$  q[i]*
c$$$  area_coords[i][p]*
c$$$  area_coords[j][p]*
c$$$  weights[p];
c$$$  }
c$$$  }
c$$$  }
c$$$  
c$$$  /* obtain the inverse of A */
c$$$  invert(A, Ainv);
c$$$  
c$$$  /* obtain projection u */
c$$$  for(j = 0; j < 3; j++) {
c$$$  u[j] = 0;
c$$$  for(i = 0; i < 3; i++) {
c$$$  u[j] += Ainv[j][i]*b[i];
c$$$  }
c$$$  }
c$$$  }



C***********************************************************************
C     
C     SUBROUTINE SLOPELIMITER4()
C     
C     Written 2011
C     
C     This subroutine selects the vertex based slope limiter based on
C     a Taylor Polynomial basis, and is consistent with p_adaptation.F
C     
C     cem
C     
C***********************************************************************

      SUBROUTINE SLOPELIMITER4() 

C.....Use appropriate modules

      USE GLOBAL
      USE DG

#ifdef CMPI
      USE MESSENGER
      USE MESSENGER_ELEM
#endif

      IMPLICIT NONE

      Integer k,ll,ss,lll,ell,bb,maxneigh

      REAL(SZ) fd,marea,fde

C.....We work over the master element
C.....Set initial values

      fd = slope_weight         ! reduces diffusion fd = 1 => full diffusion
      fde = fd                  ! add weight for lower order pieces (fd<1 => stronger limiting)     

      DO k=1,mne

         if (dofs(k).gt.1) then

            DO ll = 1,dofs(k)

               ZE(ll,k,NRK+2) = ZE(ll,k,IRK+1)
               QX(ll,k,NRK+2) = QX(ll,k,IRK+1)
               QY(ll,k,NRK+2) = QY(ll,k,IRK+1)

#ifdef TRACE
               iota(ll,k,NRK+2) = iota(ll,k,IRK+1)
#endif

#ifdef CHEM
               iota(ll,k,NRK+2) = iota(ll,k,IRK+1)
               iota2(ll,k,NRK+2) = iota2(ll,k,IRK+1)
#endif

            ENDDO

         elseif (dofs(k).eq.1) then

            cycle

         endif

      ENDDO   

      marea = 2.D0              !master elements area

C.....Convert initial values to the Taylor basis (multiply by Nmatrix) on base element


      ZEtaylor = 0.D0 
      QXtaylor = 0.D0
      QYtaylor = 0.D0

#ifdef TRACE
      iotataylor = 0.D0
#endif

#ifdef CHEM
      iotataylor = 0.D0
      iota2taylor = 0.D0
#endif

      do k=1,mne

         if (dofs(k).gt.1) then

            do ll=1,dofs(k)
               
               do ss=1,dofs(k)

                  ZEtaylor(k,ll,1) = ZEtaylor(k,ll,1) + Nmatrix(k,ll,ss,
     &               dofs(k)) * ZE(ss,k,nrk+2)
                  QXtaylor(k,ll,1) = QXtaylor(k,ll,1) + Nmatrix(k,ll,ss,
     &               dofs(k))* QX(ss,k,nrk+2)
                  QYtaylor(k,ll,1) = QYtaylor(k,ll,1) + Nmatrix(k,ll,ss,
     &               dofs(k)) * QY(ss,k,nrk+2)

#ifdef TRACE
                  iotataylor(k,ll,1) = iotataylor(k,ll,1)+ Nmatrix(k,ll,
     &               ss,dofs(k)) * iota(ss,k,nrk+2)
#endif                     
                  
#ifdef CHEM
                  iotataylor(k,ll,1) = iotataylor(k,ll,1)+ Nmatrix(k,ll,
     &               ss,dofs(k)) * iota(ss,k,nrk+2)
                  iota2taylor(k,ll,1) = iota2taylor(k,ll,1) + Nmatrix(k,
     &               ll,ss,dofs(k)) * iota2(ss,k,nrk+2)
#endif
                  
               enddo

            enddo

         elseif (dofs(k).eq.1) then

            cycle

         endif

      enddo

C.....Find values at vertices of base elements and neighbors


      ZEmax = -100.D0
      QXmax = -100.D0
      QYmax = -100.D0
      ZEmin = 100.D0
      QXmin = 100.D0
      QYmin = 100.D0

#ifdef TRACE
      iotamax = -100.D0
      iotamin = 100.D0
#endif

#ifdef CHEM
      iotamax = -100.D0
      iota2max = -100.D0
      iotamin = 100.D0
      iota2min = 100.D0
#endif


#ifdef CMPI

      CALL UPDATER_ELEM_MOD2(ZEtaylor,QXtaylor,QYtaylor,1,3)

#ifdef TRACE
      CALL UPDATER_ELEM_MOD2(iotataylor,QXtaylor,QYtaylor,1,2 )
#endif

#ifdef CHEM
      CALL UPDATER_ELEM_MOD2(iotataylor,iota2taylor,QYtaylor,1,2 )
#endif

      CALL UPDATERV(ZEmin)
      CALL UPDATERV(ZEmax)
      CALL UPDATERV(QXmin)
      CALL UPDATERV(QXmax)
      CALL UPDATERV(QYmin)
      CALL UPDATERV(QYmax)

#ifdef TRACE
      CALL UPDATERV(iotamax)
      CALL UPDATERV(iotamin)
#endif

#ifdef CHEM
      CALL UPDATERV(iotamax)
      CALL UPDATERV(iotamin)
      CALL UPDATERV(iota2max)
      CALL UPDATERV(iota2min)
#endif

#endif

      do ell=1,mnp

         do ll=1,minval(dofs(neigh_elem(ell,1:nneigh_elem(ell))))

C.....Find max and min values over polynomial coefficients

            ZEmax(ell,ll) = max(maxval( ZEtaylor(NEIGH_ELEM(ell,
     &         1:nneigh_elem(ell)),ll,1) ), ZEmax(ell,ll))
            QXmax(ell,ll) = max(maxval( QXtaylor(NEIGH_ELEM(ell,
     &         1:nneigh_elem(ell)),ll,1) ) , QXmax(ell,ll))
            QYmax(ell,ll) = max(maxval( QYtaylor(NEIGH_ELEM(ell,1:nneigh_elem(ell)),ll,1) )
     &           , QYmax(ell,ll))
            ZEmin(ell,ll) = min(minval( ZEtaylor(NEIGH_ELEM(ell,1:nneigh_elem(ell)),ll,1) )
     &           , ZEmin(ell,ll))
            QXmin(ell,ll) = min(minval( QXtaylor(NEIGH_ELEM(ell,1:nneigh_elem(ell)),ll,1) )
     &           , QXmin(ell,ll))
            QYmin(ell,ll) = min(minval( QYtaylor(NEIGH_ELEM(ell,1:nneigh_elem(ell)),ll,1) )
     &           , QYmin(ell,ll))

#ifdef TRACE
            iotamax(ell,ll) = max(maxval( iotataylor(NEIGH_ELEM(ell,1:nneigh_elem(ell)),ll,1) )
     &           , iotamax(ell,ll))
            iotamin(ell,ll) = min(minval( iotataylor(NEIGH_ELEM(ell,1:nneigh_elem(ell)),ll,1) )
     &           , iotamin(ell,ll))
#endif

#ifdef CHEM
            iotamax(ell,ll) = max(maxval( iotataylor(NEIGH_ELEM(ell,1:nneigh_elem(ell)),ll,1) )
     &           , iotamax(ell,ll))
            iota2max(ell,ll) = max(maxval( iota2taylor(NEIGH_ELEM(ell,1:nneigh_elem(ell)),ll,1) )
     &           , iota2max(ell,ll))
            iotamin(ell,ll) = min(minval( iotataylor(NEIGH_ELEM(ell,1:nneigh_elem(ell)),ll,1) )
     &           , iotamin(ell,ll))
            iota2min(ell,ll) = min(minval( iota2taylor(NEIGH_ELEM(ell,1:nneigh_elem(ell)),ll,1) )
     &           , iota2min(ell,ll))
#endif
            
         enddo

      enddo


#ifdef CMPI


      CALL UPDATERV(ZEmin)
      CALL UPDATERV(ZEmax)
      CALL UPDATERV(QXmin)
      CALL UPDATERV(QXmax)
      CALL UPDATERV(QYmin)
      CALL UPDATERV(QYmax)

#ifdef TRACE
      CALL UPDATERV(iotamax)
      CALL UPDATERV(iotamin)
#endif

#ifdef CHEM
      CALL UPDATERV(iotamax)
      CALL UPDATERV(iotamin)
      CALL UPDATERV(iota2max)
      CALL UPDATERV(iota2min)
#endif

#endif


C.....Must generate linear recostructions at vertices

      ZEtaylorvert = 0.D0
      QXtaylorvert = 0.D0
      Qytaylorvert = 0.D0

#ifdef TRACE
      iotataylorvert = 0.D0
#endif

#ifdef CHEM
      iotataylorvert = 0.D0
      iota2taylorvert = 0.D0
#endif

      do k=1,mne

         if (dofs(k).gt.1) then

            do ll=1,dofs(k)

               do lll=1,3

                  if (ll.eq.1) then

                     ZEtaylorvert(k,ll,lll) = ZEtaylor(k,ll,1) + ZEtaylor(k,ll+1,1)*( xi2vert(k,lll) -xi2BCb(k) )
     &                    + ZEtaylor(k,ll+2,1)*( xi1vert(k,lll) - xi1BCb(k) )
                     QXtaylorvert(k,ll,lll) = QXtaylor(k,ll,1) + QXtaylor(k,ll+1,1)*( xi2vert(k,lll) -xi2BCb(k) )
     &                    + QXtaylor(k,ll+2,1)*( xi1vert(k,lll) - xi1BCb(k) )
                     Qytaylorvert(k,ll,lll) = QYtaylor(k,ll,1) + QYtaylor(k,ll+1,1)*( xi2vert(k,lll) -xi2BCb(k) ) 
     &                    + QYtaylor(k,ll+2,1)*( xi1vert(k,lll) - xi1BCb(k) )

#ifdef TRACE
                     iotataylorvert(k,ll,lll) = iotataylor(k,ll,1) + 
     &                    iotataylor(k,ll+1,1)*( xi2vert(k,lll) - xi2BCb(k) )
     &                    + iotataylor(k,ll+2,1)*( xi1vert(k,lll) - xi1BCb(k) )
#endif

#ifdef CHEM
                     iotataylorvert(k,ll,lll) = iotataylor(k,ll,1) + 
     &                    iotataylor(k,ll+1,1)*( xi2vert(k,lll) - xi2BCb(k) )
     &                    + iotataylor(k,ll+2,1)*( xi1vert(k,lll) - xi1BCb(k) )
                     iota2taylorvert(k,ll,lll) = iota2taylor(k,ll,1) + 
     &                    iota2taylor(k,ll+1,1)*( xi2vert(k,lll) - xi2BCb(k) )
     &                    + iota2taylor(k,ll+2,1)*( xi1vert(k,lll) - xi1BCb(k) )
#endif

                  elseif (ll.gt.1) then

                     ZEtaylorvert(k,ll,lll) = ZEtaylor(k,ll,1) + 
     &                    ZEtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) ),1)*( xi2vert(k,lll) - xi2BCb(k) ) 
     &                    + ZEtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) )+1,1)*( xi1vert(k,lll) - xi1BCb(k) )
                     QXtaylorvert(k,ll,lll) = QXtaylor(k,ll,1) + 
     &                    QXtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) ),1)*( xi2vert(k,lll) - xi2BCb(k) ) 
     &                    + QXtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) )+1,1)*( xi1vert(k,lll) - xi1BCb(k) )
                     Qytaylorvert(k,ll,lll) = QYtaylor(k,ll,1) + 
     &                    QYtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) ),1)*( xi2vert(k,lll) - xi2BCb(k) ) 
     &                    + QYtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) )+1,1)*( xi1vert(k,lll) - xi1BCb(k) )

#ifdef TRACE
                     iotataylorvert(k,ll,lll) = iotataylor(k,ll,1)+
     &                    iotataylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) ),1)*( xi2vert(k,lll) - xi2BCb(k))
     &                    + iotataylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) )+1,1)*( xi1vert(k,lll) - xi1BCb(k) )
#endif

#ifdef CHEM
                     iotataylorvert(k,ll,lll) = iotataylor(k,ll,1)+
     &                    iotataylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) ),1)*( xi2vert(k,lll) - xi2BCb(k))
     &                    + iotataylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) )+1,1)*( xi1vert(k,lll) - xi1BCb(k) )
                     iota2taylorvert(k,ll,lll) = iota2taylor(k,ll,1)+
     &                    iotataylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) ),1)*( xi2vert(k,lll) - xi2BCb(k))
     &                    + iota2taylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) )+1,1)*( xi1vert(k,lll) - xi1BCb(k) )
#endif

                  else

                  endif

               enddo
               
            enddo
            
         elseif (dofs(k).eq.1) then

            cycle

         endif

      enddo

C.....Compute alphas for each variable in each order derivitive


      alphaZE0 = 0.D0
      alphaQX0 = 0.D0
      alphaQY0 = 0.D0

#ifdef TRACE
      alphaiota0 = 0.D0
#endif

#ifdef CHEM
      alphaiota0 = 0.D0
      alphaiota20 = 0.D0
#endif

      do k=1,mne

         if (dofs(k).gt.1) then

            do lll=1,3
               
               do ll=1,dofs(k)

                  if (ZEmin(nm(k,lll),ll).ne.ZEmax(nm(k,lll),ll)) then

                     if ( ( ZEtaylorvert(k,ll,lll).gt.ZEtaylor(k,ll,1) ).and.
     &                    ( abs(ZEtaylorvert(k,ll,lll)-ZEtaylor(k,ll,1)).gt.1.0E-15 ).and.
     &                    ( ZEmax(nm(k,lll),ll).ne.ZEtaylor(k,ll,1) ) ) then  

                        alphaZE0(k,ll,lll) = min(1.D0,  ( ZEmax(nm(k,lll),ll)
     &                       - ZEtaylor(k,ll,1) )/ (ZEtaylorvert(k,ll,lll) - ZEtaylor(k,ll,1)))

                     elseif ( (ZEtaylorvert(k,ll,lll).lt.ZEtaylor(k,ll,1) )
     &                       .and.( abs(ZEtaylorvert(k,ll,lll)-ZEtaylor(k,ll,1)).gt.1.0E-15 ).and.
     &                       ( ZEmin(nm(k,lll),ll).ne.ZEtaylor(k,ll,1) ) ) then 

                        alphaZE0(k,ll,lll) = min( 1.D0,( ZEmin(nm(k,lll),ll)
     &                       - ZEtaylor(k,ll,1) )/( ZEtaylorvert(k,ll,lll)-ZEtaylor(k,ll,1)))
                        
                     elseif ( ( ZEtaylorvert(k,ll,lll).eq.ZEtaylor(k,ll,1) ).or.
     &                       ( abs(ZEtaylorvert(k,ll,lll)-ZEtaylor(k,ll,1)).le.1.0E-15 ) ) then

                        alphaZE0(k,ll,lll) = 1.D0

                     endif

                  else

                     alphaZE0(k,ll,lll) = 1.D0

                  endif

#ifdef TRACE
                  if (iotamin(nm(k,lll),ll).ne.iotamax(nm(k,lll),ll)) then

                     if ( ( iotataylorvert(k,ll,lll).gt.iotataylor(k,ll,1) ).and.
     &                    ( abs(iotataylorvert(k,ll,lll)-iotataylor(k,ll,1)).gt.1.0E-15 ).and.
     &                    ( iotamax(nm(k,lll),ll).ne.iotataylor(k,ll,1) ) ) then  

                        alphaiota0(k,ll,lll) = min(1.D0,  ( iotamax(nm(k,lll),ll)
     &                       - iotataylor(k,ll,1) )/ (iotataylorvert(k,ll,lll) - iotataylor(k,ll,1))) 
                        
                     elseif ( (iotataylorvert(k,ll,lll).lt.iotataylor(k,ll,1) )
     &                       .and.( abs(iotataylorvert(k,ll,lll)-iotataylor(k,ll,1)).gt.1.0E-15 ).and.
     &                       ( iotamin(nm(k,lll),ll).ne.iotataylor(k,ll,1) ) ) then 

                        alphaiota0(k,ll,lll) = min( 1.D0,( iotamin(nm(k,lll),ll)
     &                       - iotataylor(k,ll,1) )/( iotataylorvert(k,ll,lll)-iotataylor(k,ll,1)))

                     elseif ( ( iotataylorvert(k,ll,lll).eq.iotataylor(k,ll,1) ).or.
     &                       ( abs(iotataylorvert(k,ll,lll)-iotataylor(k,ll,1)).le.1.0E-15 ) ) then

                        alphaiota0(k,ll,lll) = 1.D0

                     endif

                  else

                     alphaiota0(k,ll,lll) = 1.D0

                  endif
#endif
                  
#ifdef CHEM                 
                  if (iotamin(nm(k,lll),ll).ne.iotamax(nm(k,lll),ll)) then

                     if ( ( iotataylorvert(k,ll,lll).gt.iotataylor(k,ll,1) ).and.
     &                    ( abs(iotataylorvert(k,ll,lll)-iotataylor(k,ll,1)).gt.1.0E-15 ).and.
     &                    ( iotamax(nm(k,lll),ll).ne.iotataylor(k,ll,1) ) ) then  

                        alphaiota0(k,ll,lll) = min(1.D0,  ( iotamax(nm(k,lll),ll)
     &                       - iotataylor(k,ll,1) )/ (iotataylorvert(k,ll,lll) - iotataylor(k,ll,1)))

                     elseif ( (iotataylorvert(k,ll,lll).lt.iotataylor(k,ll,1) )
     &                       .and.( abs(iotataylorvert(k,ll,lll)-iotataylor(k,ll,1)).gt.1.0E-15 ).and.
     &                       ( iotamin(nm(k,lll),ll).ne.iotataylor(k,ll,1) ) ) then 

                        alphaiota0(k,ll,lll) = min( 1.D0,( iotamin(nm(k,lll),ll)
     &                       - iotataylor(k,ll,1) )/( iotataylorvert(k,ll,lll)-iotataylor(k,ll,1)))
                        
                     elseif ( ( iotataylorvert(k,ll,lll).eq.iotataylor(k,ll,1) ).or.
     &                       ( abs(iotataylorvert(k,ll,lll)-iotataylor(k,ll,1)).le.1.0E-15 ) ) then

                        alphaiota0(k,ll,lll) = 1.D0

                     endif

                  else

                     alphaiota0(k,ll,lll) = 1.D0

                  endif

                  if (iota2min(nm(k,lll),ll).ne.iota2max(nm(k,lll),ll)) then

                     if ( ( iota2taylorvert(k,ll,lll).gt.iota2taylor(k,ll,1) ).and.
     &                    ( abs(iota2taylorvert(k,ll,lll)-iota2taylor(k,ll,1)).gt.1.0E-15 ).and.
     &                    ( iota2max(nm(k,lll),ll).ne.iota2taylor(k,ll,1) ) ) then  

                        alphaiota20(k,ll,lll) = min(1.D0,  ( iota2max(nm(k,lll),ll)
     &                       - iota2taylor(k,ll,1) )/ (iota2taylorvert(k,ll,lll) - iota2taylor(k,ll,1)))

                     elseif ( (iota2taylorvert(k,ll,lll).lt.iota2taylor(k,ll,1) )
     &                       .and.( abs(iota2taylorvert(k,ll,lll)-iota2taylor(k,ll,1)).gt.1.0E-15 ).and.
     &                       ( iota2min(nm(k,lll),ll).ne.iota2taylor(k,ll,1) ) ) then 

                        alphaiota20(k,ll,lll) = min( 1.D0,( iota2min(nm(k,lll),ll)
     &                       - iota2taylor(k,ll,1) )/( iota2taylorvert(k,ll,lll)-iota2taylor(k,ll,1)))
                        
                     elseif ( ( iota2taylorvert(k,ll,lll).eq.iota2taylor(k,ll,1) ).or.
     &                       ( abs(iota2taylorvert(k,ll,lll)-iota2taylor(k,ll,1)).le.1.0E-15 ) ) then

                        alphaiota20(k,ll,lll) = 1.D0

                     endif

                  else

                     alphaiota20(k,ll,lll) = 1.D0

                  endif
#endif                 

                  if (QXmin(nm(k,lll),ll).ne.QXmax(nm(k,lll),ll)) then

                     if ( ( QXtaylorvert(k,ll,lll).gt.QXtaylor(k,ll,1) ).and.
     &                    ( abs(QXtaylorvert(k,ll,lll)-QXtaylor(k,ll,1)).gt.1.0E-15 ).and.
     &                    ( QXmax(nm(k,lll),ll).ne.QXtaylor(k,ll,1) ) ) then  

                        alphaQX0(k,ll,lll) = min(1.D0,  ( QXmax(nm(k,lll),ll)
     &                       - QXtaylor(k,ll,1) )/ (QXtaylorvert(k,ll,lll) - QXtaylor(k,ll,1)))

                     elseif ( (QXtaylorvert(k,ll,lll).lt.QXtaylor(k,ll,1) )
     &                       .and.( abs(QXtaylorvert(k,ll,lll)-QXtaylor(k,ll,1)).gt.1.0E-15 ).and.
     &                       ( QXmin(nm(k,lll),ll).ne.QXtaylor(k,ll,1) ) ) then 

                        alphaQX0(k,ll,lll) = min( 1.D0,( QXmin(nm(k,lll),ll)
     &                       - QXtaylor(k,ll,1) )/( QXtaylorvert(k,ll,lll)-QXtaylor(k,ll,1)))
                        
                     elseif ( ( QXtaylorvert(k,ll,lll).eq.QXtaylor(k,ll,1) ).or.
     &                       ( abs(QXtaylorvert(k,ll,lll)-QXtaylor(k,ll,1)).le.1.0E-15 ) ) then

                        alphaQX0(k,ll,lll) = 1.D0

                     endif

                  else

                     alphaQX0(k,ll,lll) = 1.D0

                  endif


                  if (QYmin(nm(k,lll),ll).ne.QYmax(nm(k,lll),ll)) then

                     if ( ( QYtaylorvert(k,ll,lll).gt.QYtaylor(k,ll,1) ).and.
     &                    ( abs(QYtaylorvert(k,ll,lll)-QYtaylor(k,ll,1)).gt.1.0E-15 ).and.
     &                    ( QYmax(nm(k,lll),ll).ne.QYtaylor(k,ll,1) ) ) then  

                        alphaQY0(k,ll,lll) = min(1.D0,  ( QYmax(nm(k,lll),ll)
     &                       - QYtaylor(k,ll,1) )/ (QYtaylorvert(k,ll,lll) - QYtaylor(k,ll,1)))

                     elseif ( (QYtaylorvert(k,ll,lll).lt.QYtaylor(k,ll,1) )
     &                       .and.( abs(QYtaylorvert(k,ll,lll)-QYtaylor(k,ll,1)).gt.1.0E-15 ).and.
     &                       ( QYmin(nm(k,lll),ll).ne.QYtaylor(k,ll,1) ) ) then 

                        alphaQY0(k,ll,lll) = min( 1.D0,( QYmin(nm(k,lll),ll)
     &                       - QYtaylor(k,ll,1) )/( QYtaylorvert(k,ll,lll)-QYtaylor(k,ll,1)))                        
                        
                     elseif ( ( QYtaylorvert(k,ll,lll).eq.QYtaylor(k,ll,1) ).or.
     &                       ( abs(QYtaylorvert(k,ll,lll)-QYtaylor(k,ll,1)).le.1.0E-15 ) ) then

                        alphaQY0(k,ll,lll) = 1.D0

                     endif

                  else

                     alphaQY0(k,ll,lll) = 1.D0

                  endif

               enddo            !lll

            enddo               !ll

         elseif (dofs(k).eq.1) then

            cycle

         endif

      enddo

C.....Find the prescribed higher limiters by finding smallest local value

      alphaZE = 0.D0
      alphaQX = 0.D0
      alphaQY = 0.D0

#ifdef TRACE
      alphaiota = 0.D0
#endif

#ifdef CHEM
      alphaiota = 0.D0
      alphaiota2 = 0.D0
#endif

      do k=1,mne

         if (dofs(k).gt.1) then

            do ll=1,dofs(k)		

               alphaZE(k,ll) = minval( alphaZE0(k,ll,:) )
               alphaQX(k,ll) = minval( alphaQX0(k,ll,:) )
               alphaQY(k,ll) = minval( alphaQY0(k,ll,:) )

#ifdef TRACE
               alphaiota(k,ll) = minval( alphaiota0(k,ll,:) )
#endif

#ifdef CHEM
               alphaiota(k,ll) = minval( alphaiota0(k,ll,:) )
               alphaiota2(k,ll) = minval( alphaiota20(k,ll,:) )
#endif

            enddo

         elseif (dofs(k).eq.1) then

            cycle

         endif

      enddo

C.... Choose smallest (minimum) alpha for derivative in x or y

      alphaZEm = 0.D0
      alphaQXm = 0.D0
      alphaQYm = 0.D0

#ifdef TRACE
      alphaiotam = 0.D0
#endif

#ifdef CHEM
      alphaiotam = 0.D0
      alphaiota2m = 0.D0
#endif

      do k=1,mne

         if (dofs(k).gt.1) then
            
            do bb = 1,pdg_el(k)

               if( (bb+1)*(bb+2)/2.le.dofs(k) ) then

                  alphaZEm(k,(bb*(bb+1))/2 + 1) = 
     &                 minval( alphaZE(k,(bb*(bb-1))/2 +1:(bb*(bb+1)) / 2  ) )
                  alphaQXm(k,(bb*(bb+1))/2 + 1) = 
     &                 minval( alphaQX(k,(bb*(bb-1))/2 +1:(bb*(bb+1)) / 2  ) )
                  alphaQYm(k,(bb*(bb+1))/2 + 1) = 
     &                 minval( alphaQY(k,(bb*(bb-1))/2 +1:(bb*(bb+1)) / 2  ) )

#ifdef TRACE
                  alphaiotam(k,(bb*(bb+1))/2 + 1) = 
     &                 minval( alphaiota(k,(bb*(bb-1))/2 +1:(bb*(bb+1)) / 2 ) )
#endif

#ifdef CHEM
                  alphaiotam(k,(bb*(bb+1))/2 + 1) = 
     &                 minval( alphaiota(k,(bb*(bb-1))/2 +1:(bb*(bb+1)) / 2  ) )
                  alphaiota2m(k,(bb*(bb+1))/2 + 1) = 
     &                 minval( alphaiota2(k,(bb*(bb-1))/2 +1:(bb*(bb+1)) / 2 ) )
#endif

               endif
               
            enddo

         elseif (dofs(k).eq.1) then

            cycle

         endif

      enddo

C.....Use max higher derivative values for final limiter value

      alphaZE_max = 0.D0
      alphaQX_max = 0.D0
      alphaQY_max = 0.D0

#ifdef TRACE
      alphaiota_max = 0.D0
#endif

#ifdef CHEM
      alphaiota_max = 0.D0
      alphaiota2_max = 0.D0
#endif

      do k=1,mne

         if (dofs(k).gt.1) then

            do bb =1,pdg_el(k)

               if( (bb+1)*(bb+2)/2.le.dofs(k)) then

                  alphaZE_max(k,(bb*(bb+1))/2 + 1) = fde*maxval( 
     &                 alphaZEm(k,(bb*(bb+1))/2 + 1:dofs(k)) )
                  alphaQX_max(k,(bb*(bb+1))/2 + 1) = fde*maxval( 
     &                 alphaQXm(k,(bb*(bb+1))/2 + 1:dofs(k)) )
                  alphaQY_max(k,(bb*(bb+1))/2 + 1) = fde*maxval( 
     &                 alphaQYm(k,(bb*(bb+1))/2 + 1:dofs(k)) )

#ifdef TRACE
                  alphaiota_max(k,(bb*(bb+1))/2 + 1) = fde*maxval( 
     &                 alphaiotam(k,(bb*(bb+1)/2) +1:dofs(k) ) )
#endif

#ifdef CHEM
                  alphaiota_max(k,(bb*(bb+1))/2 + 1) = fde*maxval( 
     &                 alphaiotam(k,(bb*(bb+1))/2 + 1:dofs(k)) )
                  alphaiota2_max(k,(bb*(bb+1))/2 + 1) = fde*maxval( 
     &                 alphaiota2m(k,(bb*(bb+1))/2 + 1:dofs(k)) )
#endif

               endif

            enddo

         elseif (dofs(k).eq.1) then

            cycle

         endif

      enddo

C.....Limit on the Master element in the Taylor basis, via reconstruction 
C.....of unconstrained solutions with alpha constraints


      limitZE = 0.D0
      limitQX = 0.D0
      limitQY = 0.D0

      lim_count_roll = 0

#ifdef TRACE
      limitiota = 0.D0
#endif

#ifdef CHEM
      limitiota = 0.D0
      limitiota2 = 0.D0
#endif

      do k=1,mne

         lim_count = 0

         if (dofs(k).gt.1) then

            do ll=1,dofs(k)

               if ( ll.eq.1 ) then

                  limitZE(k,ll) = ZEtaylor(k,ll,1)
                  limitQX(k,ll) = QXtaylor(k,ll,1)
                  limitQY(k,ll) = QYtaylor(k,ll,1) 

#ifdef TRACE
                  limitiota(k,ll) = iotataylor(k,ll,1)
#endif

#ifdef CHEM
                  limitiota(k,ll) = iotataylor(k,ll,1)
                  limitiota2(k,ll) = iota2taylor(k,ll,1)
#endif

               elseif ( ll.ge.2 ) then
                  
                  do bb=1,pdg_el(k)

                     if ( ll.le.( ( (bb+1)*(bb+2)) / 2.D0 ).and.(ll.gt.
     &                    (bb*(bb+1)/2.D0) ) ) then

                        limitZE(k,ll) = alphaZE_max(k,(bb*(bb+1))/2 + 1) 
     &                       * ZEtaylor(k,ll,1)
                        limitQX(k,ll) = alphaQX_max(k,(bb*(bb+1))/2 + 1) 
     &                       * ZEtaylor(k,ll,1)
                        limitQY(k,ll) = alphaQY_max(k,(bb*(bb+1))/2 + 1) 
     &                       * ZEtaylor(k,ll,1)


#ifdef TRACE
                        limitiota(k,ll) = alphaiota_max(k,(bb*(bb+1))/2 + 1) 
     &                       * iotataylor(k,ll,1)
#endif

#ifdef CHEM
                        limitiota(k,ll) = alphaiota_max(k,(bb*(bb+1))/2 + 1) 
     &                       * iotataylor(k,ll,1)
                        limitiota2(k,ll) = alphaiota2_max(k,(bb*(bb+1))/2 + 1) 
     &                       * iota2taylor(k,ll,1)
#endif


c$$$  ! Make a counter to track limiting
c$$$  
c$$$  if ( ( alphaZE_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
c$$$  &                       alphaQX_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
c$$$  &                       alphaQY_max(k,(bb*(bb+1))/2 + 1).lt.1.D0 ).and.
c$$$  &                       chem_flag.ne.1.and.tracer_flag.ne.1 ) then
c$$$  
c$$$  lim_count = 1  
c$$$  
c$$$  elseif ( (alphaZE_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
c$$$  &                          alphaQX_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
c$$$  &                          alphaQY_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
c$$$  &                          alphaiota_max(k,(bb*(bb+1))/2 + 1).lt.1.D0 ).and.
c$$$  &                          chem_flag.ne.1.and.tracer_flag.eq.1 ) then
c$$$  
c$$$  lim_count = 1   
c$$$  
c$$$  elseif ( (alphaZE_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
c$$$  &                          alphaQX_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
c$$$  &                          alphaQY_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
c$$$  &                          alphaiota_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
c$$$  &                          alphaiota2_max(k,(bb*(bb+1))/2 + 1).lt.1.D0 ).and.
c$$$  &                          chem_flag.eq.1 ) then
c$$$  
c$$$  lim_count = 1 
c$$$  
c$$$  endif
                        
                     endif
                     
                  enddo
                  
               endif

            enddo

         elseif (dofs(k).eq.1) then

            cycle

         endif

                                !lim_count_roll = lim_count_roll + lim_count

      enddo

C.....Transform back to the Dubiner basis (multiply by NmatrixInv),

      ZEconst =  0.D0
      QXconst =  0.D0
      QYconst =  0.D0

#ifdef TRACE
      iotaconst =  0.D0
#endif

#ifdef CHEM
      iotaconst =  0.D0
      iota2const =  0.D0
#endif

      do k=1,mne

         if (dofs(k).gt.1) then
                                !do lll=1,3

            do ll=1,dofs(k)

               do ss=1,dofs(k)

                  ZEconst(k,ll) = ZEconst(k,ll) + NmatrixInv(k,ll,ss,dofs(k)) 
     &                 * limitZE(k,ss)
                  QXconst(k,ll) = QXconst(k,ll) + NmatrixInv(k,ll,ss,dofs(k)) 
     &                 * limitQX(k,ss)
                  QYconst(k,ll) = QYconst(k,ll) + NmatrixInv(k,ll,ss,dofs(k)) 
     &                 * limitQY(k,ss)

#ifdef TRACE
                  iotaconst(k,ll) = iotaconst(k,ll) + 
     &                 NmatrixInv(k,ll,ss,dofs(k)) * limitiota(k,ss)
#endif

#ifdef CHEM
                  iotaconst(k,ll) = iotaconst(k,ll) + 
     &                 NmatrixInv(k,ll,ss,dofs(k)) * limitiota(k,ss)
                  iota2const(k,ll) = iota2const(k,ll) + 
     &                 NmatrixInv(k,ll,ss,dofs(k)) * limitiota2(k,ss)
#endif


               enddo

            enddo

         elseif (dofs(k).eq.1) then

            cycle

         endif

      enddo

C.....Set limit values

      do k=1,mne

         if (dofs(k).gt.1) then

            do ll = 1,dofs(k)

               ZE(ll,k,irk+1) = ZEconst(k,ll)
               QX(ll,k,irk+1) = QXconst(k,ll)
               QY(ll,k,irk+1) = QYconst(k,ll)

#ifdef TRACE
               iota(ll,k,irk+1) = iotaconst(k,ll)
#endif

#ifdef CHEM
               iota(ll,k,irk+1) = iotaconst(k,ll)
               iota2(ll,k,irk+1) = iota2const(k,ll)
#endif

            enddo

         elseif (dofs(k).eq.1) then

            cycle

         endif

      enddo
      
      return
      end subroutine

#endif

#ifdef SLOPE5

C***********************************************************************
C     
C     SUBROUTINE SLOPELIMITER5()
C     
C     Written by Clint Dawson - 30 June 2010
C     01-10-2011 - cem - adapted for p_enrichment and multicomponent
C     
C     
C***********************************************************************

      SUBROUTINE SLOPELIMITER5()

C.....Use appropriate modules

      USE SIZES, ONLY : SZ,layers
      USE GLOBAL
      USE DG

#ifdef CMPI
      USE MESSENGER
#endif

      IMPLICIT NONE

C.....Declare local variables

      INTEGER L, LL, INC1,INC2,INC3,KDP,NN,IVAR,I,J,kk,k,varnum,bb,varnum_prev
      REAL(SZ) ZEC(3),ZEVERTEX(3),DIF(3),SUMLOC,SUMDIF,SIGNDIF,
     $     DIV,REDFAC,REDMAX,tmp1,tmp2,tmp3,bound
      Real(SZ) ZEMIN1(3),ZEMAX1(3),QXMIN1(3),QXMAX1(3)
      Real(SZ) QYMIN1(3),QYMAX1(3)
      Real(SZ) iotaMIN1(3),iotaMAX1(3)
      Real(SZ) iota2MIN1(3),iota2MAX1(3)
      Real(SZ), Allocatable :: ZE_MIN1(:),ZE_MAX1(:),QX_MIN1(:),QX_MAX1(:)
      Real(SZ), Allocatable :: QY_MIN1(:),QY_MAX1(:)
      Real(SZ), Allocatable :: iota_MIN1(:),iota_MAX1(:)
      Real(SZ), Allocatable :: iota2_MIN1(:),iota2_MAX1(:)
      Real(SZ), Allocatable, target :: bed_min1(:,:), bed_max1(:,:)
      Real(SZ), pointer:: arraymin(:),arraymax(:)

      Allocate ( ZE_MIN1(NP),ZE_MAX1(NP),QX_MIN1(NP) )
      Allocate ( QY_MIN1(NP),QY_MAX1(NP),QX_MAX1(NP) )
      Allocate ( iota_MIN1(NP),iota_MAX1(NP) )
      Allocate ( iota2_MIN1(NP),iota2_MAX1(NP) )
      Allocate ( bed_MIN1(NP,layers),bed_MAX1(NP,layers) )

C     FIND THE MAXIMUM AND MINIMUM OF EACH VARIABLE OVER ALL ELEMENTS 
C     SHARING A NODE

      bound = 0.0D0

      DO I = 1,NP
         ZE_MIN1(I)=99999.
         ZE_MAX1(I)=-99999.
         QX_MIN1(I)=99999.
         QX_MAX1(I)=-99999.
         QY_MIN1(I)=99999.
         QY_MAX1(I)=-99999.

#ifdef TRACE
         iota_MIN1(I)=99999.
         iota_MAX1(I)=-99999.
#endif

#ifdef CHEM
         iota_MIN1(I)=99999.
         iota_MAX1(I)=-99999.
         iota2_MIN1(I)=99999.
         iota2_MAX1(I)=-99999.
#endif

#ifdef SED_LAY
         do l=1,layers
            bed_MIN1(I,l) = 99999.
            bed_MAX1(I,l) =-99999.
         enddo
#endif

         NO_NBORS = EL_COUNT(I)

         DO J = 1,NO_NBORS
            NBOR_EL = ELETAB(I,1+J)

c     IF(WDFLG(NBOR_EL).EQ.0) CYCLE ! DON'T COUNT DRY ELEMENTS  sb 02/26/07

            ZE_DG(J) = ZE(1,NBOR_EL,IRK+1)
            QX_DG(J) = QX(1,NBOR_EL,IRK+1)
            QY_DG(J) = QY(1,NBOR_EL,IRK+1)

#ifdef TRACE
            iota_DG(J) = iota(1,NBOR_EL,IRK+1)
#endif

#ifdef CHEM
            iota_DG(J) = iota(1,NBOR_EL,IRK+1)
            iota2_DG(J) = iota2(1,NBOR_EL,IRK+1)
#endif


#ifdef SED_LAY
            do l=1,layers
               bed_DG(J,l) = bed(1,NBOR_EL,IRK+1,l)
            enddo
#endif

C     
            IF (ZE_DG(J).LT.ZE_MIN1(I))THEN
               ZE_MIN1(I)=ZE_DG(J)
            ENDIF
            IF (ZE_DG(J).GT.ZE_MAX1(I)) THEN
               ZE_MAX1(I)=ZE_DG(J)
            ENDIF
            IF (QX_DG(J).LT.QX_MIN1(I))THEN
               QX_MIN1(I)=QX_DG(J)
            ENDIF
            IF (QX_DG(J).GT.QX_MAX1(I)) THEN
               QX_MAX1(I)=QX_DG(J)
            ENDIF
            IF (QY_DG(J).LT.QY_MIN1(I))THEN
               QY_MIN1(I)=QY_DG(J)
            ENDIF
            IF (QY_DG(J).GT.QY_MAX1(I)) THEN
               QY_MAX1(I)=QY_DG(J)
            ENDIF

#ifdef TRACE
            IF (iota_DG(J).LT.iota_MIN1(I))THEN
               iota_MIN1(I)=iota_DG(J)
            ENDIF

            IF (iota_DG(J).GT.iota_MAX1(I)) THEN
               iota_MAX1(I)=iota_DG(J)
            ENDIF
#endif

#ifdef CHEM
            IF (iota_DG(J).LT.iota_MIN1(I))THEN
               iota_MIN1(I)=iota_DG(J)
            ENDIF

            IF (iota_DG(J).GT.iota_MAX1(I)) THEN
               iota_MAX1(I)=iota_DG(J)
            ENDIF

            IF (iota2_DG(J).LT.iota2_MIN1(I))THEN
               iota2_MIN1(I)=iota2_DG(J)
            ENDIF

            IF (iota2_DG(J).GT.iota2_MAX1(I)) THEN
               iota2_MAX1(I)=iota2_DG(J)
            ENDIF
#endif

#ifdef SED_LAY
            do l=1,layers
               IF (bed_DG(J,l).LT.bed_MIN1(I,l))THEN
                  bed_MIN1(I,l)=bed_DG(J,l)
               ENDIF
               
               IF (bed_DG(J,l).GT.bed_MAX1(I,l)) THEN
                  bed_MAX1(I,l)=bed_DG(J,l)
               ENDIF
            enddo
#endif

         ENDDO
      ENDDO
#ifdef CMPI

      CALL UPDATER(ZE_MIN1,ZE_MAX1,QX_MIN1,3)
      CALL UPDATER(QX_MAX1,QY_MIN1,QY_MAX1,3)

#ifdef TRACE
      CALL UPDATER(iota_MAX1,iota_MIN1,QY_MAX1,2)
#endif

#ifdef CHEM
      CALL UPDATER(iota_MAX1,iota_MIN1,QY_MAX1,2)
      CALL UPDATER(iota2_MAX1,iota2_MIN1,QY_MAX1,2)
#endif

#ifdef SED_LAY
         do l = 1,layers
            arraymax => bed_max1(:,l)
            arraymin => bed_min1(:,l)
            CALL UPDATER(arraymax,arraymin,QY_MAX1,2 )
         enddo
#endif

#endif
C     
C     LOOP OVER ELEMENTS TO CALCULATE NEW VERTEX VALUES
C     

      bb = 1

      DO I=1,NE 
         !IF(WDFLG(I).EQ.0) CYCLE ! DON'T COUNT DRY ELEMENTS  sb 02/26/07
         N1=NM(I,1)
         N2=NM(I,2)
         N3=NM(I,3) 
         
         varnum = 3

#ifdef TRACE
         varnum = 4
#endif

#ifdef CHEM
         varnum = 5
#endif

#ifdef SED_LAY
         varnum_prev = varnum
         varnum = varnum_prev + layers
#endif

         DO IVAR=1,varnum

            IF (IVAR.EQ.1) THEN
               ZEC(1)=ZE(1,I,IRK+1)
               ZEC(2)=ZE(2,I,IRK+1)
               ZEC(3)=ZE(3,I,IRK+1)
               ZEMAX1(1)=ZE_MAX1(N1)
               ZEMIN1(1)=ZE_MIN1(N1)
               ZEMAX1(2)=ZE_MAX1(N2)
               ZEMIN1(2)=ZE_MIN1(N2)
               ZEMAX1(3)=ZE_MAX1(N3)
               ZEMIN1(3)=ZE_MIN1(N3)
            ENDIF

            IF (IVAR.EQ.2) THEN
               ZEC(1)=QX(1,I,IRK+1)
               ZEC(2)=QX(2,I,IRK+1)
               ZEC(3)=QX(3,I,IRK+1)
               ZEMAX1(1)=QX_MAX1(N1)
               ZEMIN1(1)=QX_MIN1(N1)
               ZEMAX1(2)=QX_MAX1(N2)
               ZEMIN1(2)=QX_MIN1(N2)
               ZEMAX1(3)=QX_MAX1(N3)
               ZEMIN1(3)=QX_MIN1(N3)
            ENDIF

            IF (IVAR.EQ.3) THEN
               ZEC(1)=QY(1,I,IRK+1)
               ZEC(2)=QY(2,I,IRK+1)
               ZEC(3)=QY(3,I,IRK+1)
               ZEMAX1(1)=QY_MAX1(N1)
               ZEMIN1(1)=QY_MIN1(N1)
               ZEMAX1(2)=QY_MAX1(N2)
               ZEMIN1(2)=QY_MIN1(N2)
               ZEMAX1(3)=QY_MAX1(N3)
               ZEMIN1(3)=QY_MIN1(N3)
            ENDIF

#ifdef TRACE
            IF (IVAR.EQ.4) THEN
               ZEC(1)=iota(1,I,IRK+1)
               ZEC(2)=iota(2,I,IRK+1)
               ZEC(3)=iota(3,I,IRK+1)
               ZEMAX1(1)=iota_MAX1(N1)
               ZEMIN1(1)=iota_MIN1(N1)
               ZEMAX1(2)=iota_MAX1(N2)
               ZEMIN1(2)=iota_MIN1(N2)
               ZEMAX1(3)=iota_MAX1(N3)
               ZEMIN1(3)=iota_MIN1(N3)
            ENDIF
#endif

#ifdef CHEM
            IF (IVAR.EQ.4) THEN
               ZEC(1)=iota(1,I,IRK+1)
               ZEC(2)=iota(2,I,IRK+1)
               ZEC(3)=iota(3,I,IRK+1)
               ZEMAX1(1)=iota_MAX1(N1)
               ZEMIN1(1)=iota_MIN1(N1)
               ZEMAX1(2)=iota_MAX1(N2)
               ZEMIN1(2)=iota_MIN1(N2)
               ZEMAX1(3)=iota_MAX1(N3)
               ZEMIN1(3)=iota_MIN1(N3)
            ENDIF

            IF (IVAR.EQ.5) THEN
               ZEC(1)=iota2(1,I,IRK+1)
               ZEC(2)=iota2(2,I,IRK+1)
               ZEC(3)=iota2(3,I,IRK+1)
               ZEMAX1(1)=iota2_MAX1(N1)
               ZEMIN1(1)=iota2_MIN1(N1)
               ZEMAX1(2)=iota2_MAX1(N2)
               ZEMIN1(2)=iota2_MIN1(N2)
               ZEMAX1(3)=iota2_MAX1(N3)
               ZEMIN1(3)=iota2_MIN1(N3)
            ENDIF
#endif

#ifdef SED_LAY
            do l=1,layers
               if (IVAR.eq.varnum_prev+l) then
                  ZEC(1)=bed(1,I,IRK+1,l)
                  ZEC(2)=bed(2,I,IRK+1,l)
                  ZEC(3)=bed(3,I,IRK+1,l)
                  ZEMAX1(1)=bed_MAX1(N1,l)
                  ZEMIN1(1)=bed_MIN1(N1,l)
                  ZEMAX1(2)=bed_MAX1(N2,l)
                  ZEMIN1(2)=bed_MIN1(N2,l)
                  ZEMAX1(3)=bed_MAX1(N3,l)
                  ZEMIN1(3)=bed_MIN1(N3,l)
               endif
            enddo
               
            
#endif
            

C     COMPUTE THE VERTEX VALUES

            ZEVERTEX(1)=ZEC(1)
            ZEVERTEX(2)=ZEC(1)
            ZEVERTEX(3)=ZEC(1)
            DO KK=2,3
               ZEVERTEX(1)=ZEVERTEX(1)+ PHI_CORNER(KK,1,1)*ZEC(KK)
               ZEVERTEX(2)=ZEVERTEX(2)+ PHI_CORNER(KK,2,1)*ZEC(KK)
               ZEVERTEX(3)=ZEVERTEX(3)+ PHI_CORNER(KK,3,1)*ZEC(KK)
            ENDDO

            
C     RESET THE VERTEX VALUE TO BE LESS THAN OR EQUAL TO THE MAX AND
C     GREATER THAN OR EQUAL TO THE MIN AT THAT VERTEX
C     
            ZEVERTEX(1)=DMAX1(DMIN1(ZEVERTEX(1),ZEMAX1(1)),ZEMIN1(1))
            ZEVERTEX(2)=DMAX1(DMIN1(ZEVERTEX(2),ZEMAX1(2)),ZEMIN1(2))
            ZEVERTEX(3)=DMAX1(DMIN1(ZEVERTEX(3),ZEMAX1(3)),ZEMIN1(3))

            tmp1 = ZEVERTEX(1)
            tmp2 = ZEVERTEX(2)
            tmp3 = ZEVERTEX(3)


C     LOOP OVER THE VERTICES 3 TIMES
C     IF THE VALUE AT THE VERTEX IS ABOVE (BELOW) THE MAX (MIN) AT THAT
C     VERTEX THEN SUBTRACT OFF THE DIFFERENCE AND ADD IT TO THE OTHER
C     VERTICES
C     
            DO LL=1,3
               SUMLOC=(ZEVERTEX(1)+ZEVERTEX(2)+ZEVERTEX(3))/3.0D0
               SUMDIF=(SUMLOC-ZEC(1))*3.0D0
               SIGNDIF=DSIGN(1.D0,SUMDIF)
               DIF(1)=(ZEVERTEX(1)-ZEC(1))*SIGNDIF
               DIF(2)=(ZEVERTEX(2)-ZEC(1))*SIGNDIF
               DIF(3)=(ZEVERTEX(3)-ZEC(1))*SIGNDIF
               INC1=0
               IF (DIF(1).GT.0) INC1=1
               INC2=0
               IF (DIF(2).GT.0) INC2=1
               INC3=0
               IF (DIF(3).GT.0) INC3=1
               KDP=INC1+INC2+INC3
C     
               DO K=1,3
                  DIV=DMAX1(1.D0,DFLOAT(KDP))
                  IF (DIF(K).GT.0) THEN
                     REDFAC=SUMDIF*SIGNDIF/DIV
                     KDP=KDP-1
                  ELSE
                     REDFAC=0
                  ENDIF
                  IF (SIGNDIF.GT.0) THEN
                     REDMAX=ZEVERTEX(K)-ZEMIN1(K)
                  ELSE
                     REDMAX=ZEMAX1(K)-ZEVERTEX(K)
                  ENDIF
                  REDFAC=DMIN1(REDFAC,REDMAX)
                  SUMDIF=SUMDIF-REDFAC*SIGNDIF
                  ZEVERTEX(K)=ZEVERTEX(K)-REDFAC*SIGNDIF
               ENDDO
            ENDDO
            IF (IVAR.EQ.1) THEN

c$$$               if( (abs(tmp1-ZEVERTEX(1)).ge.bound).or.
c$$$     &              (abs(tmp2-ZEVERTEX(2)).ge.bound).or.
c$$$     &              (abs(tmp3-ZEVERTEX(3)).ge.bound).and.
c$$$     &              (slopeflag.eq.5) ) then
c$$$
c$$$                  ZE(4:dofs(i),i,irk+1)=0.D0
c$$$
c$$$                  if (padapt.eq.1) then
c$$$                     pdg_el(i) = 1
c$$$                  endif 
c$$$
c$$$               endif

               ZE(2,I,IRK+1)=-1.d0/6.d0*(ZEVERTEX(1)+ZEVERTEX(2))
     $              +1.d0/3.d0*ZEVERTEX(3)
               ZE(3,I,IRK+1)=-.5d0*ZEVERTEX(1)+.5d0*ZEVERTEX(2)
            ENDIF
            IF (IVAR.EQ.2) THEN

c$$$               if( (abs(tmp1-ZEVERTEX(1)).ge.bound).or.
c$$$     &              (abs(tmp2-ZEVERTEX(2)).ge.bound).or.
c$$$     &              (abs(tmp3-ZEVERTEX(3)).ge.bound).and.
c$$$     &              (slopeflag.eq.5) ) then
c$$$
c$$$
c$$$                  QX(4:dofs(i),i,irk+1)=0.D0
c$$$
c$$$                  if (padapt.eq.1) then
c$$$                     pdg_el(i) = 1
c$$$                  endif 
c$$$
c$$$               endif

               QX(2,I,IRK+1)=-1.d0/6.d0*(ZEVERTEX(1)+ZEVERTEX(2))
     $              +1.d0/3.d0*ZEVERTEX(3)
               QX(3,I,IRK+1)=-.5d0*ZEVERTEX(1)+.5d0*ZEVERTEX(2)
            ENDIF
            IF (IVAR.EQ.3) THEN

c$$$               if( (abs(tmp1-ZEVERTEX(1)).ge.bound).or.
c$$$     &              (abs(tmp2-ZEVERTEX(2)).ge.bound).or.
c$$$     &              (abs(tmp3-ZEVERTEX(3)).ge.bound).and.
c$$$     &              (slopeflag.eq.5) ) then
c$$$
c$$$                  QY(4:dofs(i),i,irk+1)=0.D0
c$$$
c$$$                  if (padapt.eq.1) then
c$$$                     pdg_el(i) = 1
c$$$                  endif 
c$$$
c$$$               endif

               QY(2,I,IRK+1)=-1.d0/6.d0*(ZEVERTEX(1)+ZEVERTEX(2))
     $              +1.d0/3.d0*ZEVERTEX(3)
               QY(3,I,IRK+1)=-.5d0*ZEVERTEX(1)+.5d0*ZEVERTEX(2)
            ENDIF

#ifdef TRACE
            IF (IVAR.EQ.4) THEN
               
c$$$               if( (abs(tmp1-ZEVERTEX(1)).ge.bound).or.
c$$$     &              (abs(tmp2-ZEVERTEX(2)).ge.bound).or.
c$$$     &              (abs(tmp3-ZEVERTEX(3)).ge.bound).and.
c$$$     &              (slopeflag.eq.5) ) then
c$$$
c$$$
c$$$                  iota(4:dofs(i),i,irk+1)=0.D0
c$$$
c$$$                  if (padapt.eq.1) then
c$$$                     pdg_el(i) = 1
c$$$                  endif 
c$$$
c$$$               endif
               
               iota(2,I,IRK+1)=-1.d0/6.d0*(ZEVERTEX(1)+ZEVERTEX(2))
     $              +1.d0/3.d0*ZEVERTEX(3)
               iota(3,I,IRK+1)=-.5d0*ZEVERTEX(1)+.5d0*ZEVERTEX(2)

            ENDIF
#endif

#ifdef CHEM
            IF (IVAR.EQ.4) THEN
               
c$$$               if( (abs(tmp1-ZEVERTEX(1)).ge.bound).or.
c$$$     &              (abs(tmp2-ZEVERTEX(2)).ge.bound).or.
c$$$     &              (abs(tmp3-ZEVERTEX(3)).ge.bound).and.
c$$$     &              (slopeflag.eq.5) ) then
c$$$
c$$$
c$$$                  iota(4:dofs(i),i,irk+1)=0.D0
c$$$
c$$$                  if (padapt.eq.1) then
c$$$                     pdg_el(i) = 1
c$$$                  endif 
c$$$
c$$$               endif
               
               iota(2,I,IRK+1)=-1.d0/6.d0*(ZEVERTEX(1)+ZEVERTEX(2))
     $              +1.d0/3.d0*ZEVERTEX(3)
               iota(3,I,IRK+1)=-.5d0*ZEVERTEX(1)+.5d0*ZEVERTEX(2)

            ENDIF
            IF (IVAR.EQ.5) THEN

c$$$               if( (abs(tmp1-ZEVERTEX(1)).ge.bound).or.
c$$$     &              (abs(tmp2-ZEVERTEX(2)).ge.bound).or.
c$$$     &              (abs(tmp3-ZEVERTEX(3)).ge.bound).and.
c$$$     &              (slopeflag.eq.5) ) then
c$$$
c$$$
c$$$                  iota2(4:dofs(i),i,irk+1)=0.D0
c$$$
c$$$                  if (padapt.eq.1) then
c$$$                     pdg_el(i) = 1
c$$$                  endif 
c$$$
c$$$               endif

               iota2(2,I,IRK+1)=-1.d0/6.d0*(ZEVERTEX(1)+ZEVERTEX(2))
     $              +1.d0/3.d0*ZEVERTEX(3)
               iota2(3,I,IRK+1)=-.5d0*ZEVERTEX(1)+.5d0*ZEVERTEX(2)
            ENDIF
#endif
         ENDDO

#ifdef SED_LAY

         do l=1,layers
            if (IVAR.eq.varnum_prev+l) then

c$$$              if( (abs(tmp1-ZEVERTEX(1)).ge.bound).or.
c$$$     &              (abs(tmp2-ZEVERTEX(2)).ge.bound).or.
c$$$     &              (abs(tmp3-ZEVERTEX(3)).ge.bound).and.
c$$$     &              (slopeflag.eq.5) ) then
c$$$
c$$$
c$$$                  bed(4:dofs(i),i,irk+1,l)=0.D0
c$$$
c$$$                  if (padapt.eq.1) then
c$$$                     pdg_el(i) = 1
c$$$                  endif 
c$$$
c$$$               endif
               
               bed(2,I,IRK+1,l)=-1.d0/6.d0*(ZEVERTEX(1)+ZEVERTEX(2))
     $              +1.d0/3.d0*ZEVERTEX(3)
               bed(3,I,IRK+1,l)=-.5d0*ZEVERTEX(1)+.5d0*ZEVERTEX(2)
            
            
            endif
         enddo

#endif

                                !endif

      ENDDO
      RETURN
      END SUBROUTINE 

#endif

#ifdef STBLZR

C***********************************************************************
C     
C     SUBROUTINE SLOPELIMITER5()
C     
C     Written by Clint Dawson - 30 June 2010
C     01-10-2011 - cem - adapted for p_enrichment and multicomponent
C     
C     
C***********************************************************************

      SUBROUTINE SLOPELIMITER5()

C.....Use appropriate modules

      USE SIZES, ONLY : SZ,layers
      USE GLOBAL
      USE DG

#ifdef CMPI
      USE MESSENGER
#endif

      IMPLICIT NONE

C.....Declare local variables

      INTEGER L, LL, INC1,INC2,INC3,KDP,NN,IVAR,I,J,kk,k,varnum,bb,varnum_prev
      REAL(SZ) ZEC(3),ZEVERTEX(3),DIF(3),SUMLOC,SUMDIF,SIGNDIF,
     $     DIV,REDFAC,REDMAX,tmp1,tmp2,tmp3,bound
      Real(SZ) ZEMIN1(3),ZEMAX1(3),QXMIN1(3),QXMAX1(3)
      Real(SZ) QYMIN1(3),QYMAX1(3)
      Real(SZ) iotaMIN1(3),iotaMAX1(3)
      Real(SZ) iota2MIN1(3),iota2MAX1(3)
      Real(SZ), Allocatable :: ZE_MIN1(:),ZE_MAX1(:),QX_MIN1(:),QX_MAX1(:)
      Real(SZ), Allocatable :: QY_MIN1(:),QY_MAX1(:)
      Real(SZ), Allocatable :: iota_MIN1(:),iota_MAX1(:)
      Real(SZ), Allocatable :: iota2_MIN1(:),iota2_MAX1(:)
      Real(SZ), Allocatable, target :: bed_min1(:,:), bed_max1(:,:)
      Real(SZ), pointer:: arraymin(:),arraymax(:)

      Allocate ( ZE_MIN1(NP),ZE_MAX1(NP),QX_MIN1(NP) )
      Allocate ( QY_MIN1(NP),QY_MAX1(NP),QX_MAX1(NP) )
      Allocate ( iota_MIN1(NP),iota_MAX1(NP) )
      Allocate ( iota2_MIN1(NP),iota2_MAX1(NP) )
      Allocate ( bed_MIN1(NP,layers),bed_MAX1(NP,layers) )

C     FIND THE MAXIMUM AND MINIMUM OF EACH VARIABLE OVER ALL ELEMENTS 
C     SHARING A NODE

      bound = 0.0D0

      DO I = 1,NP
         ZE_MIN1(I)=99999.
         ZE_MAX1(I)=-99999.
         QX_MIN1(I)=99999.
         QX_MAX1(I)=-99999.
         QY_MIN1(I)=99999.
         QY_MAX1(I)=-99999.

#ifdef TRACE
         iota_MIN1(I)=99999.
         iota_MAX1(I)=-99999.
#endif

#ifdef CHEM
         iota_MIN1(I)=99999.
         iota_MAX1(I)=-99999.
         iota2_MIN1(I)=99999.
         iota2_MAX1(I)=-99999.
#endif

#ifdef SED_LAY
         do l=1,layers
            bed_MIN1(I,l) = 99999.
            bed_MAX1(I,l) =-99999.
         enddo
#endif

         NO_NBORS = EL_COUNT(I)

         DO J = 1,NO_NBORS
            NBOR_EL = ELETAB(I,1+J)

c     IF(WDFLG(NBOR_EL).EQ.0) CYCLE ! DON'T COUNT DRY ELEMENTS  sb 02/26/07

            ZE_DG(J) = ZE(1,NBOR_EL,IRK+1)
            QX_DG(J) = QX(1,NBOR_EL,IRK+1)
            QY_DG(J) = QY(1,NBOR_EL,IRK+1)

#ifdef TRACE
            iota_DG(J) = iota(1,NBOR_EL,IRK+1)
#endif

#ifdef CHEM
            iota_DG(J) = iota(1,NBOR_EL,IRK+1)
            iota2_DG(J) = iota2(1,NBOR_EL,IRK+1)
#endif


#ifdef SED_LAY
            do l=1,layers
               bed_DG(J,l) = bed(1,NBOR_EL,IRK+1,l)
            enddo
#endif

C     
            IF (ZE_DG(J).LT.ZE_MIN1(I))THEN
               ZE_MIN1(I)=ZE_DG(J)
            ENDIF
            IF (ZE_DG(J).GT.ZE_MAX1(I)) THEN
               ZE_MAX1(I)=ZE_DG(J)
            ENDIF
            IF (QX_DG(J).LT.QX_MIN1(I))THEN
               QX_MIN1(I)=QX_DG(J)
            ENDIF
            IF (QX_DG(J).GT.QX_MAX1(I)) THEN
               QX_MAX1(I)=QX_DG(J)
            ENDIF
            IF (QY_DG(J).LT.QY_MIN1(I))THEN
               QY_MIN1(I)=QY_DG(J)
            ENDIF
            IF (QY_DG(J).GT.QY_MAX1(I)) THEN
               QY_MAX1(I)=QY_DG(J)
            ENDIF

#ifdef TRACE
            IF (iota_DG(J).LT.iota_MIN1(I))THEN
               iota_MIN1(I)=iota_DG(J)
            ENDIF

            IF (iota_DG(J).GT.iota_MAX1(I)) THEN
               iota_MAX1(I)=iota_DG(J)
            ENDIF
#endif

#ifdef CHEM
            IF (iota_DG(J).LT.iota_MIN1(I))THEN
               iota_MIN1(I)=iota_DG(J)
            ENDIF

            IF (iota_DG(J).GT.iota_MAX1(I)) THEN
               iota_MAX1(I)=iota_DG(J)
            ENDIF

            IF (iota2_DG(J).LT.iota2_MIN1(I))THEN
               iota2_MIN1(I)=iota2_DG(J)
            ENDIF

            IF (iota2_DG(J).GT.iota2_MAX1(I)) THEN
               iota2_MAX1(I)=iota2_DG(J)
            ENDIF
#endif

#ifdef SED_LAY
            do l=1,layers
               IF (bed_DG(J,l).LT.bed_MIN1(I,l))THEN
                  bed_MIN1(I,l)=bed_DG(J,l)
               ENDIF
               
               IF (bed_DG(J,l).GT.bed_MAX1(I,l)) THEN
                  bed_MAX1(I,l)=bed_DG(J,l)
               ENDIF
            enddo
#endif

         ENDDO
      ENDDO
#ifdef CMPI

      CALL UPDATER(ZE_MIN1,ZE_MAX1,QX_MIN1,3)
      CALL UPDATER(QX_MAX1,QY_MIN1,QY_MAX1,3)

#ifdef TRACE
      CALL UPDATER(iota_MAX1,iota_MIN1,QY_MAX1,2)
#endif

#ifdef CHEM
      CALL UPDATER(iota_MAX1,iota_MIN1,QY_MAX1,2)
      CALL UPDATER(iota2_MAX1,iota2_MIN1,QY_MAX1,2)
#endif

#ifdef SED_LAY
         do l = 1,layers
            arraymax => bed_max1(:,l)
            arraymin => bed_min1(:,l)
            CALL UPDATER(arraymax,arraymin,QY_MAX1,2 )
         enddo
#endif

#endif
C     
C     LOOP OVER ELEMENTS TO CALCULATE NEW VERTEX VALUES
C     

      bb = 1

      DO I=1,NE 
         !IF(WDFLG(I).EQ.0) CYCLE ! DON'T COUNT DRY ELEMENTS  sb 02/26/07
         N1=NM(I,1)
         N2=NM(I,2)
         N3=NM(I,3) 
         
         varnum = 3

#ifdef TRACE
         varnum = 4
#endif

#ifdef CHEM
         varnum = 5
#endif

#ifdef SED_LAY
         varnum_prev = varnum
         varnum = varnum_prev + layers
#endif

         DO IVAR=1,varnum

            IF (IVAR.EQ.1) THEN
               ZEC(1)=ZE(1,I,IRK+1)
               ZEC(2)=ZE(2,I,IRK+1)
               ZEC(3)=ZE(3,I,IRK+1)
               ZEMAX1(1)=ZE_MAX1(N1)
               ZEMIN1(1)=ZE_MIN1(N1)
               ZEMAX1(2)=ZE_MAX1(N2)
               ZEMIN1(2)=ZE_MIN1(N2)
               ZEMAX1(3)=ZE_MAX1(N3)
               ZEMIN1(3)=ZE_MIN1(N3)
            ENDIF

            IF (IVAR.EQ.2) THEN
               ZEC(1)=QX(1,I,IRK+1)
               ZEC(2)=QX(2,I,IRK+1)
               ZEC(3)=QX(3,I,IRK+1)
               ZEMAX1(1)=QX_MAX1(N1)
               ZEMIN1(1)=QX_MIN1(N1)
               ZEMAX1(2)=QX_MAX1(N2)
               ZEMIN1(2)=QX_MIN1(N2)
               ZEMAX1(3)=QX_MAX1(N3)
               ZEMIN1(3)=QX_MIN1(N3)
            ENDIF

            IF (IVAR.EQ.3) THEN
               ZEC(1)=QY(1,I,IRK+1)
               ZEC(2)=QY(2,I,IRK+1)
               ZEC(3)=QY(3,I,IRK+1)
               ZEMAX1(1)=QY_MAX1(N1)
               ZEMIN1(1)=QY_MIN1(N1)
               ZEMAX1(2)=QY_MAX1(N2)
               ZEMIN1(2)=QY_MIN1(N2)
               ZEMAX1(3)=QY_MAX1(N3)
               ZEMIN1(3)=QY_MIN1(N3)
            ENDIF

#ifdef TRACE
            IF (IVAR.EQ.4) THEN
               ZEC(1)=iota(1,I,IRK+1)
               ZEC(2)=iota(2,I,IRK+1)
               ZEC(3)=iota(3,I,IRK+1)
               ZEMAX1(1)=iota_MAX1(N1)
               ZEMIN1(1)=iota_MIN1(N1)
               ZEMAX1(2)=iota_MAX1(N2)
               ZEMIN1(2)=iota_MIN1(N2)
               ZEMAX1(3)=iota_MAX1(N3)
               ZEMIN1(3)=iota_MIN1(N3)
            ENDIF
#endif

#ifdef CHEM
            IF (IVAR.EQ.4) THEN
               ZEC(1)=iota(1,I,IRK+1)
               ZEC(2)=iota(2,I,IRK+1)
               ZEC(3)=iota(3,I,IRK+1)
               ZEMAX1(1)=iota_MAX1(N1)
               ZEMIN1(1)=iota_MIN1(N1)
               ZEMAX1(2)=iota_MAX1(N2)
               ZEMIN1(2)=iota_MIN1(N2)
               ZEMAX1(3)=iota_MAX1(N3)
               ZEMIN1(3)=iota_MIN1(N3)
            ENDIF

            IF (IVAR.EQ.5) THEN
               ZEC(1)=iota2(1,I,IRK+1)
               ZEC(2)=iota2(2,I,IRK+1)
               ZEC(3)=iota2(3,I,IRK+1)
               ZEMAX1(1)=iota2_MAX1(N1)
               ZEMIN1(1)=iota2_MIN1(N1)
               ZEMAX1(2)=iota2_MAX1(N2)
               ZEMIN1(2)=iota2_MIN1(N2)
               ZEMAX1(3)=iota2_MAX1(N3)
               ZEMIN1(3)=iota2_MIN1(N3)
            ENDIF
#endif

#ifdef SED_LAY
            do l=1,layers
               if (IVAR.eq.varnum_prev+l) then
                  ZEC(1)=bed(1,I,IRK+1,l)
                  ZEC(2)=bed(2,I,IRK+1,l)
                  ZEC(3)=bed(3,I,IRK+1,l)
                  ZEMAX1(1)=bed_MAX1(N1,l)
                  ZEMIN1(1)=bed_MIN1(N1,l)
                  ZEMAX1(2)=bed_MAX1(N2,l)
                  ZEMIN1(2)=bed_MIN1(N2,l)
                  ZEMAX1(3)=bed_MAX1(N3,l)
                  ZEMIN1(3)=bed_MIN1(N3,l)
               endif
            enddo
               
            
#endif
            

C     COMPUTE THE VERTEX VALUES

            ZEVERTEX(1)=ZEC(1)
            ZEVERTEX(2)=ZEC(1)
            ZEVERTEX(3)=ZEC(1)
            DO KK=2,3
               ZEVERTEX(1)=ZEVERTEX(1)+ PHI_CORNER(KK,1,1)*ZEC(KK)
               ZEVERTEX(2)=ZEVERTEX(2)+ PHI_CORNER(KK,2,1)*ZEC(KK)
               ZEVERTEX(3)=ZEVERTEX(3)+ PHI_CORNER(KK,3,1)*ZEC(KK)
            ENDDO

            
C     RESET THE VERTEX VALUE TO BE LESS THAN OR EQUAL TO THE MAX AND
C     GREATER THAN OR EQUAL TO THE MIN AT THAT VERTEX
C     
            ZEVERTEX(1)=DMAX1(DMIN1(ZEVERTEX(1),ZEMAX1(1)),ZEMIN1(1))
            ZEVERTEX(2)=DMAX1(DMIN1(ZEVERTEX(2),ZEMAX1(2)),ZEMIN1(2))
            ZEVERTEX(3)=DMAX1(DMIN1(ZEVERTEX(3),ZEMAX1(3)),ZEMIN1(3))

            tmp1 = ZEVERTEX(1)
            tmp2 = ZEVERTEX(2)
            tmp3 = ZEVERTEX(3)


C     LOOP OVER THE VERTICES 3 TIMES
C     IF THE VALUE AT THE VERTEX IS ABOVE (BELOW) THE MAX (MIN) AT THAT
C     VERTEX THEN SUBTRACT OFF THE DIFFERENCE AND ADD IT TO THE OTHER
C     VERTICES
C     
            DO LL=1,3
               SUMLOC=(ZEVERTEX(1)+ZEVERTEX(2)+ZEVERTEX(3))/3.0D0
               SUMDIF=(SUMLOC-ZEC(1))*3.0D0
               SIGNDIF=DSIGN(1.D0,SUMDIF)
               DIF(1)=(ZEVERTEX(1)-ZEC(1))*SIGNDIF
               DIF(2)=(ZEVERTEX(2)-ZEC(1))*SIGNDIF
               DIF(3)=(ZEVERTEX(3)-ZEC(1))*SIGNDIF
               INC1=0
               IF (DIF(1).GT.0) INC1=1
               INC2=0
               IF (DIF(2).GT.0) INC2=1
               INC3=0
               IF (DIF(3).GT.0) INC3=1
               KDP=INC1+INC2+INC3
C     
               DO K=1,3
                  DIV=DMAX1(1.D0,DFLOAT(KDP))
                  IF (DIF(K).GT.0) THEN
                     REDFAC=SUMDIF*SIGNDIF/DIV
                     KDP=KDP-1
                  ELSE
                     REDFAC=0
                  ENDIF
                  IF (SIGNDIF.GT.0) THEN
                     REDMAX=ZEVERTEX(K)-ZEMIN1(K)
                  ELSE
                     REDMAX=ZEMAX1(K)-ZEVERTEX(K)
                  ENDIF
                  REDFAC=DMIN1(REDFAC,REDMAX)
                  SUMDIF=SUMDIF-REDFAC*SIGNDIF
                  ZEVERTEX(K)=ZEVERTEX(K)-REDFAC*SIGNDIF
               ENDDO
            ENDDO
            IF (IVAR.EQ.1) THEN

c$$$               if( (abs(tmp1-ZEVERTEX(1)).ge.bound).or.
c$$$     &              (abs(tmp2-ZEVERTEX(2)).ge.bound).or.
c$$$     &              (abs(tmp3-ZEVERTEX(3)).ge.bound).and.
c$$$     &              (slopeflag.eq.5) ) then
c$$$
c$$$                  ZE(4:dofs(i),i,irk+1)=0.D0
c$$$
c$$$                  if (padapt.eq.1) then
c$$$                     pdg_el(i) = 1
c$$$                  endif 
c$$$
c$$$               endif

               ZE(2,I,IRK+1)=-1.d0/6.d0*(ZEVERTEX(1)+ZEVERTEX(2))
     $              +1.d0/3.d0*ZEVERTEX(3)
               ZE(3,I,IRK+1)=-.5d0*ZEVERTEX(1)+.5d0*ZEVERTEX(2)
            ENDIF
            IF (IVAR.EQ.2) THEN

c$$$               if( (abs(tmp1-ZEVERTEX(1)).ge.bound).or.
c$$$     &              (abs(tmp2-ZEVERTEX(2)).ge.bound).or.
c$$$     &              (abs(tmp3-ZEVERTEX(3)).ge.bound).and.
c$$$     &              (slopeflag.eq.5) ) then
c$$$
c$$$
c$$$                  QX(4:dofs(i),i,irk+1)=0.D0
c$$$
c$$$                  if (padapt.eq.1) then
c$$$                     pdg_el(i) = 1
c$$$                  endif 
c$$$
c$$$               endif

               QX(2,I,IRK+1)=-1.d0/6.d0*(ZEVERTEX(1)+ZEVERTEX(2))
     $              +1.d0/3.d0*ZEVERTEX(3)
               QX(3,I,IRK+1)=-.5d0*ZEVERTEX(1)+.5d0*ZEVERTEX(2)
            ENDIF
            IF (IVAR.EQ.3) THEN

c$$$               if( (abs(tmp1-ZEVERTEX(1)).ge.bound).or.
c$$$     &              (abs(tmp2-ZEVERTEX(2)).ge.bound).or.
c$$$     &              (abs(tmp3-ZEVERTEX(3)).ge.bound).and.
c$$$     &              (slopeflag.eq.5) ) then
c$$$
c$$$                  QY(4:dofs(i),i,irk+1)=0.D0
c$$$
c$$$                  if (padapt.eq.1) then
c$$$                     pdg_el(i) = 1
c$$$                  endif 
c$$$
c$$$               endif

               QY(2,I,IRK+1)=-1.d0/6.d0*(ZEVERTEX(1)+ZEVERTEX(2))
     $              +1.d0/3.d0*ZEVERTEX(3)
               QY(3,I,IRK+1)=-.5d0*ZEVERTEX(1)+.5d0*ZEVERTEX(2)
            ENDIF

#ifdef TRACE
            IF (IVAR.EQ.4) THEN
               
c$$$               if( (abs(tmp1-ZEVERTEX(1)).ge.bound).or.
c$$$     &              (abs(tmp2-ZEVERTEX(2)).ge.bound).or.
c$$$     &              (abs(tmp3-ZEVERTEX(3)).ge.bound).and.
c$$$     &              (slopeflag.eq.5) ) then
c$$$
c$$$
c$$$                  iota(4:dofs(i),i,irk+1)=0.D0
c$$$
c$$$                  if (padapt.eq.1) then
c$$$                     pdg_el(i) = 1
c$$$                  endif 
c$$$
c$$$               endif
               
               iota(2,I,IRK+1)=-1.d0/6.d0*(ZEVERTEX(1)+ZEVERTEX(2))
     $              +1.d0/3.d0*ZEVERTEX(3)
               iota(3,I,IRK+1)=-.5d0*ZEVERTEX(1)+.5d0*ZEVERTEX(2)

            ENDIF
#endif

#ifdef CHEM
            IF (IVAR.EQ.4) THEN
               
c$$$               if( (abs(tmp1-ZEVERTEX(1)).ge.bound).or.
c$$$     &              (abs(tmp2-ZEVERTEX(2)).ge.bound).or.
c$$$     &              (abs(tmp3-ZEVERTEX(3)).ge.bound).and.
c$$$     &              (slopeflag.eq.5) ) then
c$$$
c$$$
c$$$                  iota(4:dofs(i),i,irk+1)=0.D0
c$$$
c$$$                  if (padapt.eq.1) then
c$$$                     pdg_el(i) = 1
c$$$                  endif 
c$$$
c$$$               endif
               
               iota(2,I,IRK+1)=-1.d0/6.d0*(ZEVERTEX(1)+ZEVERTEX(2))
     $              +1.d0/3.d0*ZEVERTEX(3)
               iota(3,I,IRK+1)=-.5d0*ZEVERTEX(1)+.5d0*ZEVERTEX(2)

            ENDIF
            IF (IVAR.EQ.5) THEN

c$$$               if( (abs(tmp1-ZEVERTEX(1)).ge.bound).or.
c$$$     &              (abs(tmp2-ZEVERTEX(2)).ge.bound).or.
c$$$     &              (abs(tmp3-ZEVERTEX(3)).ge.bound).and.
c$$$     &              (slopeflag.eq.5) ) then
c$$$
c$$$
c$$$                  iota2(4:dofs(i),i,irk+1)=0.D0
c$$$
c$$$                  if (padapt.eq.1) then
c$$$                     pdg_el(i) = 1
c$$$                  endif 
c$$$
c$$$               endif

               iota2(2,I,IRK+1)=-1.d0/6.d0*(ZEVERTEX(1)+ZEVERTEX(2))
     $              +1.d0/3.d0*ZEVERTEX(3)
               iota2(3,I,IRK+1)=-.5d0*ZEVERTEX(1)+.5d0*ZEVERTEX(2)
            ENDIF
#endif
         ENDDO

#ifdef SED_LAY

         do l=1,layers
            if (IVAR.eq.varnum_prev+l) then

c$$$              if( (abs(tmp1-ZEVERTEX(1)).ge.bound).or.
c$$$     &              (abs(tmp2-ZEVERTEX(2)).ge.bound).or.
c$$$     &              (abs(tmp3-ZEVERTEX(3)).ge.bound).and.
c$$$     &              (slopeflag.eq.5) ) then
c$$$
c$$$
c$$$                  bed(4:dofs(i),i,irk+1,l)=0.D0
c$$$
c$$$                  if (padapt.eq.1) then
c$$$                     pdg_el(i) = 1
c$$$                  endif 
c$$$
c$$$               endif
               
               bed(2,I,IRK+1,l)=-1.d0/6.d0*(ZEVERTEX(1)+ZEVERTEX(2))
     $              +1.d0/3.d0*ZEVERTEX(3)
               bed(3,I,IRK+1,l)=-.5d0*ZEVERTEX(1)+.5d0*ZEVERTEX(2)
            
            
            endif
         enddo

#endif

                                !endif

      ENDDO
      RETURN
      END SUBROUTINE 

#endif


#ifdef SLOPEALL

C***********************************************************************
C     
C     SUBROUTINE SLOPELIMITER5()
C     
C     Written by Clint Dawson - 30 June 2010
C     01-10-2011 - cem - adapted for p_enrichment and multicomponent
C     
C     
C***********************************************************************

      SUBROUTINE SLOPELIMITER5()

C.....Use appropriate modules

      USE SIZES, ONLY : SZ,layers
      USE GLOBAL
      USE DG

#ifdef CMPI
      USE MESSENGER
#endif

      IMPLICIT NONE

C.....Declare local variables

      INTEGER L, LL, INC1,INC2,INC3,KDP,NN,IVAR,I,J,kk,k,varnum,bb,varnum_prev
      REAL(SZ) ZEC(3),ZEVERTEX(3),DIF(3),SUMLOC,SUMDIF,SIGNDIF,
     $     DIV,REDFAC,REDMAX,tmp1,tmp2,tmp3,bound
      Real(SZ) ZEMIN1(3),ZEMAX1(3),QXMIN1(3),QXMAX1(3)
      Real(SZ) QYMIN1(3),QYMAX1(3)
      Real(SZ) iotaMIN1(3),iotaMAX1(3)
      Real(SZ) iota2MIN1(3),iota2MAX1(3)
      Real(SZ), Allocatable :: ZE_MIN1(:),ZE_MAX1(:),QX_MIN1(:),QX_MAX1(:)
      Real(SZ), Allocatable :: QY_MIN1(:),QY_MAX1(:)
      Real(SZ), Allocatable :: iota_MIN1(:),iota_MAX1(:)
      Real(SZ), Allocatable :: iota2_MIN1(:),iota2_MAX1(:)
      Real(SZ), Allocatable, target :: bed_min1(:,:), bed_max1(:,:)
      Real(SZ), pointer:: arraymin(:),arraymax(:)

      Allocate ( ZE_MIN1(NP),ZE_MAX1(NP),QX_MIN1(NP) )
      Allocate ( QY_MIN1(NP),QY_MAX1(NP),QX_MAX1(NP) )
      Allocate ( iota_MIN1(NP),iota_MAX1(NP) )
      Allocate ( iota2_MIN1(NP),iota2_MAX1(NP) )
      Allocate ( bed_MIN1(NP,layers),bed_MAX1(NP,layers) )

C     FIND THE MAXIMUM AND MINIMUM OF EACH VARIABLE OVER ALL ELEMENTS 
C     SHARING A NODE

      bound = 1.0E-5 

      DO I = 1,NP
         ZE_MIN1(I)=99999.
         ZE_MAX1(I)=-99999.
         QX_MIN1(I)=99999.
         QX_MAX1(I)=-99999.
         QY_MIN1(I)=99999.
         QY_MAX1(I)=-99999.

#ifdef TRACE
         iota_MIN1(I)=99999.
         iota_MAX1(I)=-99999.
#endif

#ifdef CHEM
         iota_MIN1(I)=99999.
         iota_MAX1(I)=-99999.
         iota2_MIN1(I)=99999.
         iota2_MAX1(I)=-99999.
#endif

#ifdef SED_LAY
         do l=1,layers
            bed_MIN1(I,l) = 99999.
            bed_MAX1(I,l) =-99999.
         enddo
#endif

         NO_NBORS = EL_COUNT(I)

         DO J = 1,NO_NBORS
            NBOR_EL = ELETAB(I,1+J)

c     IF(WDFLG(NBOR_EL).EQ.0) CYCLE ! DON'T COUNT DRY ELEMENTS  sb 02/26/07

            ZE_DG(J) = ZE(1,NBOR_EL,IRK+1)
            QX_DG(J) = QX(1,NBOR_EL,IRK+1)
            QY_DG(J) = QY(1,NBOR_EL,IRK+1)

#ifdef TRACE
            iota_DG(J) = iota(1,NBOR_EL,IRK+1)
#endif

#ifdef CHEM
            iota_DG(J) = iota(1,NBOR_EL,IRK+1)
            iota2_DG(J) = iota2(1,NBOR_EL,IRK+1)
#endif


#ifdef SED_LAY
            do l=1,layers
               bed_DG(J,l) = bed(1,NBOR_EL,IRK+1,l)
            enddo
#endif

C     
            IF (ZE_DG(J).LT.ZE_MIN1(I))THEN
               ZE_MIN1(I)=ZE_DG(J)
            ENDIF
            IF (ZE_DG(J).GT.ZE_MAX1(I)) THEN
               ZE_MAX1(I)=ZE_DG(J)
            ENDIF
            IF (QX_DG(J).LT.QX_MIN1(I))THEN
               QX_MIN1(I)=QX_DG(J)
            ENDIF
            IF (QX_DG(J).GT.QX_MAX1(I)) THEN
               QX_MAX1(I)=QX_DG(J)
            ENDIF
            IF (QY_DG(J).LT.QY_MIN1(I))THEN
               QY_MIN1(I)=QY_DG(J)
            ENDIF
            IF (QY_DG(J).GT.QY_MAX1(I)) THEN
               QY_MAX1(I)=QY_DG(J)
            ENDIF

#ifdef TRACE
            IF (iota_DG(J).LT.iota_MIN1(I))THEN
               iota_MIN1(I)=iota_DG(J)
            ENDIF

            IF (iota_DG(J).GT.iota_MAX1(I)) THEN
               iota_MAX1(I)=iota_DG(J)
            ENDIF
#endif

#ifdef CHEM
            IF (iota_DG(J).LT.iota_MIN1(I))THEN
               iota_MIN1(I)=iota_DG(J)
            ENDIF

            IF (iota_DG(J).GT.iota_MAX1(I)) THEN
               iota_MAX1(I)=iota_DG(J)
            ENDIF

            IF (iota2_DG(J).LT.iota2_MIN1(I))THEN
               iota2_MIN1(I)=iota2_DG(J)
            ENDIF

            IF (iota2_DG(J).GT.iota2_MAX1(I)) THEN
               iota2_MAX1(I)=iota2_DG(J)
            ENDIF
#endif

#ifdef SED_LAY
            do l=1,layers
               IF (bed_DG(J,l).LT.bed_MIN1(I,l))THEN
                  bed_MIN1(I,l)=bed_DG(J,l)
               ENDIF
               
               IF (bed_DG(J,l).GT.bed_MAX1(I,l)) THEN
                  bed_MAX1(I,l)=bed_DG(J,l)
               ENDIF
            enddo
#endif

         ENDDO
      ENDDO
#ifdef CMPI

      CALL UPDATER(ZE_MIN1,ZE_MAX1,QX_MIN1,3)
      CALL UPDATER(QX_MAX1,QY_MIN1,QY_MAX1,3)

#ifdef TRACE
      CALL UPDATER(iota_MAX1,iota_MIN1,QY_MAX1,2)
#endif

#ifdef CHEM
      CALL UPDATER(iota_MAX1,iota_MIN1,QY_MAX1,2)
      CALL UPDATER(iota2_MAX1,iota2_MIN1,QY_MAX1,2)
#endif

#ifdef SED_LAY
         do l = 1,layers
            arraymax => bed_max1(:,l)
            arraymin => bed_min1(:,l)
            CALL UPDATER(arraymax,arraymin,QY_MAX1,2 )
         enddo
#endif

#endif
C     
C     LOOP OVER ELEMENTS TO CALCULATE NEW VERTEX VALUES
C     

      bb = 1

      DO I=1,NE 

                                !if (dofs(i).eq.3) then

         !IF(WDFLG(I).EQ.0) CYCLE ! DON'T COUNT DRY ELEMENTS  sb 02/26/07
         N1=NM(I,1)
         N2=NM(I,2)
         N3=NM(I,3) 

         varnum = 3

#ifdef TRACE
         varnum = 4
#endif

#ifdef CHEM
         varnum = 5
#endif

#ifdef SED_LAY
         varnum_prev = varnum
         varnum = varnum_prev + layers
#endif

         DO IVAR=1,varnum

            IF (IVAR.EQ.1) THEN
               ZEC(1)=ZE(1,I,IRK+1)
               ZEC(2)=ZE(2,I,IRK+1)
               ZEC(3)=ZE(3,I,IRK+1)
               ZEMAX1(1)=ZE_MAX1(N1)
               ZEMIN1(1)=ZE_MIN1(N1)
               ZEMAX1(2)=ZE_MAX1(N2)
               ZEMIN1(2)=ZE_MIN1(N2)
               ZEMAX1(3)=ZE_MAX1(N3)
               ZEMIN1(3)=ZE_MIN1(N3)
            ENDIF

            IF (IVAR.EQ.2) THEN
               ZEC(1)=QX(1,I,IRK+1)
               ZEC(2)=QX(2,I,IRK+1)
               ZEC(3)=QX(3,I,IRK+1)
               ZEMAX1(1)=QX_MAX1(N1)
               ZEMIN1(1)=QX_MIN1(N1)
               ZEMAX1(2)=QX_MAX1(N2)
               ZEMIN1(2)=QX_MIN1(N2)
               ZEMAX1(3)=QX_MAX1(N3)
               ZEMIN1(3)=QX_MIN1(N3)
            ENDIF

            IF (IVAR.EQ.3) THEN
               ZEC(1)=QY(1,I,IRK+1)
               ZEC(2)=QY(2,I,IRK+1)
               ZEC(3)=QY(3,I,IRK+1)
               ZEMAX1(1)=QY_MAX1(N1)
               ZEMIN1(1)=QY_MIN1(N1)
               ZEMAX1(2)=QY_MAX1(N2)
               ZEMIN1(2)=QY_MIN1(N2)
               ZEMAX1(3)=QY_MAX1(N3)
               ZEMIN1(3)=QY_MIN1(N3)
            ENDIF

#ifdef TRACE
            IF (IVAR.EQ.4) THEN
               ZEC(1)=iota(1,I,IRK+1)
               ZEC(2)=iota(2,I,IRK+1)
               ZEC(3)=iota(3,I,IRK+1)
               ZEMAX1(1)=iota_MAX1(N1)
               ZEMIN1(1)=iota_MIN1(N1)
               ZEMAX1(2)=iota_MAX1(N2)
               ZEMIN1(2)=iota_MIN1(N2)
               ZEMAX1(3)=iota_MAX1(N3)
               ZEMIN1(3)=iota_MIN1(N3)
            ENDIF
#endif

#ifdef CHEM
            IF (IVAR.EQ.4) THEN
               ZEC(1)=iota(1,I,IRK+1)
               ZEC(2)=iota(2,I,IRK+1)
               ZEC(3)=iota(3,I,IRK+1)
               ZEMAX1(1)=iota_MAX1(N1)
               ZEMIN1(1)=iota_MIN1(N1)
               ZEMAX1(2)=iota_MAX1(N2)
               ZEMIN1(2)=iota_MIN1(N2)
               ZEMAX1(3)=iota_MAX1(N3)
               ZEMIN1(3)=iota_MIN1(N3)
            ENDIF

            IF (IVAR.EQ.5) THEN
               ZEC(1)=iota2(1,I,IRK+1)
               ZEC(2)=iota2(2,I,IRK+1)
               ZEC(3)=iota2(3,I,IRK+1)
               ZEMAX1(1)=iota2_MAX1(N1)
               ZEMIN1(1)=iota2_MIN1(N1)
               ZEMAX1(2)=iota2_MAX1(N2)
               ZEMIN1(2)=iota2_MIN1(N2)
               ZEMAX1(3)=iota2_MAX1(N3)
               ZEMIN1(3)=iota2_MIN1(N3)
            ENDIF
#endif

#ifdef SED_LAY
            do l=1,layers
               if (IVAR.eq.varnum_prev+l) then
                  ZEC(1)=bed(1,I,IRK+1,l)
                  ZEC(2)=bed(2,I,IRK+1,l)
                  ZEC(3)=bed(3,I,IRK+1,l)
                  ZEMAX1(1)=bed_MAX1(N1,l)
                  ZEMIN1(1)=bed_MIN1(N1,l)
                  ZEMAX1(2)=bed_MAX1(N2,l)
                  ZEMIN1(2)=bed_MIN1(N2,l)
                  ZEMAX1(3)=bed_MAX1(N3,l)
                  ZEMIN1(3)=bed_MIN1(N3,l)
               endif
            enddo
               
            
#endif
            

C     COMPUTE THE VERTEX VALUES

            ZEVERTEX(1)=ZEC(1)
            ZEVERTEX(2)=ZEC(1)
            ZEVERTEX(3)=ZEC(1)
            DO KK=2,3
               ZEVERTEX(1)=ZEVERTEX(1)+ PHI_CORNER(KK,1,1)*ZEC(KK)
               ZEVERTEX(2)=ZEVERTEX(2)+ PHI_CORNER(KK,2,1)*ZEC(KK)
               ZEVERTEX(3)=ZEVERTEX(3)+ PHI_CORNER(KK,3,1)*ZEC(KK)
            ENDDO

            
C     RESET THE VERTEX VALUE TO BE LESS THAN OR EQUAL TO THE MAX AND
C     GREATER THAN OR EQUAL TO THE MIN AT THAT VERTEX
C     
            ZEVERTEX(1)=DMAX1(DMIN1(ZEVERTEX(1),ZEMAX1(1)),ZEMIN1(1))
            ZEVERTEX(2)=DMAX1(DMIN1(ZEVERTEX(2),ZEMAX1(2)),ZEMIN1(2))
            ZEVERTEX(3)=DMAX1(DMIN1(ZEVERTEX(3),ZEMAX1(3)),ZEMIN1(3))

            tmp1 = ZEVERTEX(1)
            tmp2 = ZEVERTEX(2)
            tmp3 = ZEVERTEX(3)


C     LOOP OVER THE VERTICES 3 TIMES
C     IF THE VALUE AT THE VERTEX IS ABOVE (BELOW) THE MAX (MIN) AT THAT
C     VERTEX THEN SUBTRACT OFF THE DIFFERENCE AND ADD IT TO THE OTHER
C     VERTICES
C     
            DO LL=1,3
               SUMLOC=(ZEVERTEX(1)+ZEVERTEX(2)+ZEVERTEX(3))/3.0D0
               SUMDIF=(SUMLOC-ZEC(1))*3.0D0
               SIGNDIF=DSIGN(1.D0,SUMDIF)
               DIF(1)=(ZEVERTEX(1)-ZEC(1))*SIGNDIF
               DIF(2)=(ZEVERTEX(2)-ZEC(1))*SIGNDIF
               DIF(3)=(ZEVERTEX(3)-ZEC(1))*SIGNDIF
               INC1=0
               IF (DIF(1).GT.0) INC1=1
               INC2=0
               IF (DIF(2).GT.0) INC2=1
               INC3=0
               IF (DIF(3).GT.0) INC3=1
               KDP=INC1+INC2+INC3
C     
               DO K=1,3
                  DIV=DMAX1(1.D0,DFLOAT(KDP))
                  IF (DIF(K).GT.0) THEN
                     REDFAC=SUMDIF*SIGNDIF/DIV
                     KDP=KDP-1
                  ELSE
                     REDFAC=0
                  ENDIF
                  IF (SIGNDIF.GT.0) THEN
                     REDMAX=ZEVERTEX(K)-ZEMIN1(K)
                  ELSE
                     REDMAX=ZEMAX1(K)-ZEVERTEX(K)
                  ENDIF
                  REDFAC=DMIN1(REDFAC,REDMAX)
                  SUMDIF=SUMDIF-REDFAC*SIGNDIF
                  ZEVERTEX(K)=ZEVERTEX(K)-REDFAC*SIGNDIF
               ENDDO
            ENDDO
            IF (IVAR.EQ.1) THEN

c$$$               if( (abs(tmp1-ZEVERTEX(1)).gt.bound).or.
c$$$     &              (abs(tmp2-ZEVERTEX(2)).gt.bound).or.
c$$$     &              (abs(tmp3-ZEVERTEX(3)).gt.bound).and.
c$$$     &              (slopeflag.eq.5) ) then
c$$$
c$$$                  ZE(4:dofs(i),i,irk+1)=0.D0
c$$$
c$$$                  if (padapt.eq.1) then
c$$$                     pdg_el(i) = 1
c$$$                  endif 
c$$$
c$$$               endif

               ZE(2,I,IRK+1)=-1.d0/6.d0*(ZEVERTEX(1)+ZEVERTEX(2))
     $              +1.d0/3.d0*ZEVERTEX(3)
               ZE(3,I,IRK+1)=-.5d0*ZEVERTEX(1)+.5d0*ZEVERTEX(2)
            ENDIF
            IF (IVAR.EQ.2) THEN

c$$$               if( (abs(tmp1-ZEVERTEX(1)).gt.bound).or.
c$$$     &              (abs(tmp2-ZEVERTEX(2)).gt.bound).or.
c$$$     &              (abs(tmp3-ZEVERTEX(3)).gt.bound).and.
c$$$     &              (slopeflag.eq.5) ) then
c$$$
c$$$
c$$$                  QX(4:dofs(i),i,irk+1)=0.D0
c$$$
c$$$                  if (padapt.eq.1) then
c$$$                     pdg_el(i) = 1
c$$$                  endif 
c$$$
c$$$               endif

               QX(2,I,IRK+1)=-1.d0/6.d0*(ZEVERTEX(1)+ZEVERTEX(2))
     $              +1.d0/3.d0*ZEVERTEX(3)
               QX(3,I,IRK+1)=-.5d0*ZEVERTEX(1)+.5d0*ZEVERTEX(2)
            ENDIF
            IF (IVAR.EQ.3) THEN

c$$$               if( (abs(tmp1-ZEVERTEX(1)).gt.bound).or.
c$$$     &              (abs(tmp2-ZEVERTEX(2)).gt.bound).or.
c$$$     &              (abs(tmp3-ZEVERTEX(3)).gt.bound).and.
c$$$     &              (slopeflag.eq.5) ) then
c$$$
c$$$                  QY(4:dofs(i),i,irk+1)=0.D0
c$$$
c$$$                  if (padapt.eq.1) then
c$$$                     pdg_el(i) = 1
c$$$                  endif 
c$$$
c$$$               endif

               QY(2,I,IRK+1)=-1.d0/6.d0*(ZEVERTEX(1)+ZEVERTEX(2))
     $              +1.d0/3.d0*ZEVERTEX(3)
               QY(3,I,IRK+1)=-.5d0*ZEVERTEX(1)+.5d0*ZEVERTEX(2)
            ENDIF

#ifdef TRACE
            IF (IVAR.EQ.4) THEN
               
c$$$               if( (abs(tmp1-ZEVERTEX(1)).gt.bound).or.
c$$$     &              (abs(tmp2-ZEVERTEX(2)).gt.bound).or.
c$$$     &              (abs(tmp3-ZEVERTEX(3)).gt.bound).and.
c$$$     &              (slopeflag.eq.5) ) then
c$$$
c$$$
c$$$                  iota(4:dofs(i),i,irk+1)=0.D0
c$$$
c$$$                  if (padapt.eq.1) then
c$$$                     pdg_el(i) = 1
c$$$                  endif 
c$$$
c$$$               endif
               
               iota(2,I,IRK+1)=-1.d0/6.d0*(ZEVERTEX(1)+ZEVERTEX(2))
     $              +1.d0/3.d0*ZEVERTEX(3)
               iota(3,I,IRK+1)=-.5d0*ZEVERTEX(1)+.5d0*ZEVERTEX(2)

            ENDIF
#endif

#ifdef CHEM
            IF (IVAR.EQ.4) THEN
               
c$$$               if( (abs(tmp1-ZEVERTEX(1)).gt.bound).or.
c$$$     &              (abs(tmp2-ZEVERTEX(2)).gt.bound).or.
c$$$     &              (abs(tmp3-ZEVERTEX(3)).gt.bound).and.
c$$$     &              (slopeflag.eq.5) ) then
c$$$
c$$$
c$$$                  iota(4:dofs(i),i,irk+1)=0.D0
c$$$
c$$$                  if (padapt.eq.1) then
c$$$                     pdg_el(i) = 1
c$$$                  endif 
c$$$
c$$$               endif
               
               iota(2,I,IRK+1)=-1.d0/6.d0*(ZEVERTEX(1)+ZEVERTEX(2))
     $              +1.d0/3.d0*ZEVERTEX(3)
               iota(3,I,IRK+1)=-.5d0*ZEVERTEX(1)+.5d0*ZEVERTEX(2)

            ENDIF
            IF (IVAR.EQ.5) THEN

c$$$               if( (abs(tmp1-ZEVERTEX(1)).gt.bound).or.
c$$$     &              (abs(tmp2-ZEVERTEX(2)).gt.bound).or.
c$$$     &              (abs(tmp3-ZEVERTEX(3)).gt.bound).and.
c$$$     &              (slopeflag.eq.5) ) then
c$$$
c$$$
c$$$                  iota2(4:dofs(i),i,irk+1)=0.D0
c$$$
c$$$                  if (padapt.eq.1) then
c$$$                     pdg_el(i) = 1
c$$$                  endif 
c$$$
c$$$               endif

               iota2(2,I,IRK+1)=-1.d0/6.d0*(ZEVERTEX(1)+ZEVERTEX(2))
     $              +1.d0/3.d0*ZEVERTEX(3)
               iota2(3,I,IRK+1)=-.5d0*ZEVERTEX(1)+.5d0*ZEVERTEX(2)
            ENDIF
#endif
         ENDDO

#ifdef SED_LAY

         do l=1,layers
            if (IVAR.eq.varnum_prev+l) then

c$$$              if( (abs(tmp1-ZEVERTEX(1)).gt.bound).or.
c$$$     &              (abs(tmp2-ZEVERTEX(2)).gt.bound).or.
c$$$     &              (abs(tmp3-ZEVERTEX(3)).gt.bound).and.
c$$$     &              (slopeflag.eq.5) ) then
c$$$
c$$$
c$$$                  bed(4:dofs(i),i,irk+1,l)=0.D0
c$$$
c$$$                  if (padapt.eq.1) then
c$$$                     pdg_el(i) = 1
c$$$                  endif 
c$$$
c$$$               endif
               
               bed(2,I,IRK+1,l)=-1.d0/6.d0*(ZEVERTEX(1)+ZEVERTEX(2))
     $              +1.d0/3.d0*ZEVERTEX(3)
               bed(3,I,IRK+1,l)=-.5d0*ZEVERTEX(1)+.5d0*ZEVERTEX(2)
            
            
            endif
         enddo

#endif

                                !endif

      ENDDO
      RETURN
      END SUBROUTINE 

C***********************************************************************
C     
C     SUBROUTINE SLOPELIMITER6()
C     
C     Written 2011
C     
C     This subroutine selects the Barth--Jespersen slope limiter based on
C     a Taylor Polynomial basis, and is consistent with p-adaptation
C     to arbitrary order p
C     
C     -- cem
C     
C***********************************************************************

      SUBROUTINE SLOPELIMITER6()

C.....Use appropriate modules

      USE GLOBAL
      USE DG

#ifdef CMPI
      USE MESSENGER_ELEM
#endif

      IMPLICIT NONE

      Integer k,ll,ss,lll,ell,bb

      REAL(SZ) fd,marea
      Real(SZ), Allocatable :: ZEminel(:,:),ZEmaxel(:,:),QXminel(:,:),QXmaxel(:,:)
      Real(SZ), Allocatable :: QYminel(:,:),QYmaxel(:,:)
      Real(SZ), Allocatable :: iotaminel(:,:),iotamaxel(:,:)
      Real(SZ), Allocatable :: iota2minel(:,:),iota2maxel(:,:)

      Allocate ( ZEminel(mne,dofh),ZEmaxel(mne,dofh),QXminel(mne,dofh) )
      Allocate ( QYminel(mne,dofh),QYmaxel(mne,dofh),QXmaxel(mne,dofh) )
      Allocate ( iotaminel(mne,dofh),iotamaxel(mne,dofh) )
      Allocate ( iota2minel(mne,dofh),iota2maxel(mne,dofh) )


C.....We work over the master element
C.....Set initial values

      fd = slope_weight         ! add weight for lower order pieces (fd<1 => stronger limiting)
      

      DO k=1,NE

         if (dofs(k).gt.1) then

            DO ll = 1,dofs(k)

               ZE(ll,k,NRK+2) = ZE(ll,k,IRK+1)
               QX(ll,k,NRK+2) = QX(ll,k,IRK+1)
               QY(ll,k,NRK+2) = QY(ll,k,IRK+1)

#ifdef TRACE
               iota(ll,k,NRK+2) = iota(ll,k,IRK+1)
#endif

#ifdef CHEM
               iota(ll,k,NRK+2) = iota(ll,k,IRK+1)
               iota2(ll,k,NRK+2) = iota2(ll,k,IRK+1)
#endif


            ENDDO

         elseif (dofs(k).eq.1) then

            cycle

         endif

      ENDDO   

      marea = 2.D0              !master elements area

C.....Convert initial values to the Taylor basis (multiply by Nmatrix) on base element


      ZEtaylor = 0.D0 
      QXtaylor = 0.D0
      QYtaylor = 0.D0

#ifdef TRACE
      iotataylor = 0.D0
#endif

#ifdef CHEM
      iotataylor = 0.D0
      iota2taylor = 0.D0
#endif

      do k=1,ne

         if (dofs(k).gt.1) then

            do ll=1,dofs(k)
               
               do ss=1,dofs(k)

                  ZEtaylor(k,ll,1) = ZEtaylor(k,ll,1) + Nmatrix(k,ll,ss,dofs(k)) * ZE(ss,k,nrk+2)
                  QXtaylor(k,ll,1) = QXtaylor(k,ll,1) + Nmatrix(k,ll,ss,dofs(k)) * QX(ss,k,nrk+2)
                  QYtaylor(k,ll,1) = QYtaylor(k,ll,1) + Nmatrix(k,ll,ss,dofs(k)) * QY(ss,k,nrk+2)

#ifdef TRACE
                  iotataylor(k,ll,1) = iotataylor(k,ll,1) + Nmatrix(k,ll,ss,dofs(k))*iota(ss,k,nrk+2)
#endif                     
                  
#ifdef CHEM
                  iotataylor(k,ll,1) = iotataylor(k,ll,1) + Nmatrix(k,ll,ss,dofs(k))*iota(ss,k,nrk+2)
                  iota2taylor(k,ll,1) = iota2taylor(k,ll,1) + Nmatrix(k,ll,ss,dofs(k))*iota2(ss,k,nrk+2)
#endif
                  
               enddo

            enddo

         elseif (dofs(k).eq.1) then

            cycle

         endif

      enddo

C.....Find values at vertices of base elements and neighbors


      ZEmaxel = -100.D0
      QXmaxel = -100.D0
      QYmaxel = -100.D0
      ZEminel = 100.D0
      QXminel = 100.D0
      QYminel = 100.D0

#ifdef TRACE
      iotamaxel = -100.D0
      iotaminel = 100.D0
#endif

#ifdef CHEM
      iotamaxel = -100.D0
      iota2maxel = -100.D0
      iotaminel = 100.D0
      iota2minel = 100.D0
#endif

#ifdef CMPI

      CALL UPDATER_ELEM_MOD2(ZEtaylor,QXtaylor,QYtaylor,1,3)

#ifdef TRACE
      CALL UPDATER_ELEM_MOD2(iotataylor,QXtaylor,QYtaylor,1,2 )  
#endif

#ifdef CHEM
      CALL UPDATER_ELEM_MOD2(iotataylor,iota2taylor,QYtaylor,1,2 )  
#endif

#endif

      do k = 1,ne

         do ll = 1,dofs(k) 

            do ell = 1,3        ! Number of edge neighbors for a triangle

C.....Find max and min values over polynomial coefficients
               
               ZEmaxel(k,ll) = max( ZEtaylor(k,ll,1),ZEtaylor(EL_NBORS(ell,k),ll,1), ZEmaxel(k,ll) )
               QXmaxel(k,ll) = max( QXtaylor(k,ll,1),QXtaylor(EL_NBORS(ell,k),ll,1), QXmaxel(k,ll) )
               QYmaxel(k,ll) = max( QYtaylor(k,ll,1),QYtaylor(EL_NBORS(ell,k),ll,1), QYmaxel(k,ll) )

               ZEminel(k,ll) = min( ZEtaylor(k,ll,1),ZEtaylor(EL_NBORS(ell,k),ll,1), ZEminel(k,ll) )
               QXminel(k,ll) = min( QXtaylor(k,ll,1),QXtaylor(EL_NBORS(ell,k),ll,1), QXminel(k,ll) )
               QYminel(k,ll) = min( QYtaylor(k,ll,1),QYtaylor(EL_NBORS(ell,k),ll,1), QYminel(k,ll) )

#ifdef TRACE
               iotamaxel(k,ll) = max( iotataylor(k,ll,1),iotataylor(EL_NBORS(ell,k),ll,1), 
     &              iotamaxel(k,ll) )
               iotaminel(k,ll) = min( iotataylor(k,ll,1),iotataylor(EL_NBORS(ell,k),ll,1), 
     &              iotaminel(k,ll) )
#endif

#ifdef CHEM
               iotamaxel(k,ll) = max( iotataylor(k,ll,1),iotataylor(EL_NBORS(ell,k),ll,1), 
     &              iotamaxel(k,ll) )
               iota2maxel(k,ll) = max( iota2taylor(k,ll,1),iota2taylor(EL_NBORS(ell,k),ll,1),
     &              iota2maxel(k,ll) )
               iotaminel(k,ll) = min( iotataylor(k,ll,1),iotataylor(EL_NBORS(ell,k),ll,1),
     &              iotaminel(k,ll) )
               iota2minel(k,ll) = min( iota2taylor(k,ll,1),iota2taylor(EL_NBORS(ell,k),ll,1),
     &              iota2minel(k,ll) )
#endif
               
            enddo
            
         enddo

      enddo



C.....Must generate linear recostructions at vertices

      ZEtaylorvert = 0.D0
      QXtaylorvert = 0.D0
      Qytaylorvert = 0.D0

#ifdef TRACE
      iotataylorvert = 0.D0
#endif

#ifdef CHEM
      iotataylorvert = 0.D0
      iota2taylorvert = 0.D0
#endif

      do k=1,ne

         if (dofs(k).gt.1) then

            do ll=1,dofs(k)

               do lll=1,3

                  if (ll.eq.1) then

                     ZEtaylorvert(k,ll,lll) = ZEtaylor(k,ll,1) + ZEtaylor(k,ll+1,1)*( xi2vert(k,lll) -xi2BCb(k) )
     &                    + ZEtaylor(k,ll+2,1)*( xi1vert(k,lll) - xi1BCb(k) )
                     QXtaylorvert(k,ll,lll) = QXtaylor(k,ll,1) + QXtaylor(k,ll+1,1)*( xi2vert(k,lll) -xi2BCb(k) )
     &                    + QXtaylor(k,ll+2,1)*( xi1vert(k,lll) - xi1BCb(k) )
                     Qytaylorvert(k,ll,lll) = QYtaylor(k,ll,1) + QYtaylor(k,ll+1,1)*( xi2vert(k,lll) -xi2BCb(k) ) 
     &                    + QYtaylor(k,ll+2,1)*( xi1vert(k,lll) - xi1BCb(k) )

#ifdef TRACE
                     iotataylorvert(k,ll,lll) = iotataylor(k,ll,1) + 
     &                    iotataylor(k,ll+1,1)*( xi2vert(k,lll) - xi2BCb(k) )
     &                    + iotataylor(k,ll+2,1)*( xi1vert(k,lll) - xi1BCb(k) )
#endif

#ifdef CHEM
                     iotataylorvert(k,ll,lll) = iotataylor(k,ll,1) + 
     &                    iotataylor(k,ll+1,1)*( xi2vert(k,lll) - xi2BCb(k) )
     &                    + iotataylor(k,ll+2,1)*( xi1vert(k,lll) - xi1BCb(k) )
                     iota2taylorvert(k,ll,lll) = iota2taylor(k,ll,1) + 
     &                    iota2taylor(k,ll+1,1)*( xi2vert(k,lll) - xi2BCb(k) )
     &                    + iota2taylor(k,ll+2,1)*( xi1vert(k,lll) - xi1BCb(k) )
#endif

                  elseif (ll.gt.1) then

                     ZEtaylorvert(k,ll,lll) = ZEtaylor(k,ll,1) + 
     &                    ZEtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) ),1)*( xi2vert(k,lll) - xi2BCb(k) ) 
     &                    + ZEtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) )+1,1)*( xi1vert(k,lll) - xi1BCb(k) )
                     QXtaylorvert(k,ll,lll) = QXtaylor(k,ll,1) + 
     &                    QXtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) ),1)*( xi2vert(k,lll) - xi2BCb(k) ) 
     &                    + QXtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) )+1,1)*( xi1vert(k,lll) - xi1BCb(k) )
                     Qytaylorvert(k,ll,lll) = QYtaylor(k,ll,1) + 
     &                    QYtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) ),1)*( xi2vert(k,lll) - xi2BCb(k) ) 
     &                    + QYtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) )+1,1)*( xi1vert(k,lll) - xi1BCb(k) )

#ifdef TRACE
                     iotataylorvert(k,ll,lll) = iotataylor(k,ll,1)+
     &                    iotataylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) ),1)*( xi2vert(k,lll) - xi2BCb(k))
     &                    + iotataylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) )+1,1)*( xi1vert(k,lll) - xi1BCb(k) )
#endif

#ifdef CHEM
                     iotataylorvert(k,ll,lll) = iotataylor(k,ll,1)+
     &                    iotataylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) ),1)*( xi2vert(k,lll) - xi2BCb(k))
     &                    + iotataylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) )+1,1)*( xi1vert(k,lll) - xi1BCb(k) )
                     iota2taylorvert(k,ll,lll) = iota2taylor(k,ll,1)+
     &                    iotataylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) ),1)*( xi2vert(k,lll) - xi2BCb(k))
     &                    + iota2taylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) )+1,1)*( xi1vert(k,lll) - xi1BCb(k) )
#endif

                  else

                  endif

               enddo
               
            enddo
            
         elseif (dofs(k).eq.1) then

            cycle

         endif

      enddo

C.....Compute alphas for each variable in each order derivitive

      alphaZE0 = 0.D0
      alphaQX0 = 0.D0
      alphaQY0 = 0.D0

#ifdef TRACE
      alphaiota0 = 0.D0
#endif

#ifdef CHEM
      alphaiota0 = 0.D0
      alphaiota20 = 0.D0
#endif

      do k = 1,ne

         if (dofs(k).gt.1) then

            do ll = 1,dofs(k)

               do lll = 1,3

                  if ( ZEtaylorvert(k,ll,lll).gt.ZEtaylor(k,ll,1).and.
     &                 abs((ZEtaylorvert(k,ll,lll)-ZEtaylor(k,ll,1))).gt.1.0E-15 ) then 

                     alphaZE0(k,ll,lll) = min( 1.D0, ( ZEmaxel(k,ll) - 
     &                    ZEtaylor(k,ll,1) ) / (  ZEtaylorvert(k,ll,lll) - ZEtaylor(k,ll,1) ) )
                     
                  elseif ( ZEtaylorvert(k,ll,lll).eq.ZEtaylor(k,ll,1).or.
     &                    abs((ZEtaylorvert(k,ll,lll)-ZEtaylor(k,ll,1))).le.1.0E-15 ) then

                     alphaZE0(k,ll,lll) = 1.D0

                  elseif (  ZEtaylorvert(k,ll,lll).lt.ZEtaylor(k,ll,1).and.
     &                    abs((ZEtaylorvert(k,ll,lll)-ZEtaylor(k,ll,1))).gt.1.0E-15 ) then

                     alphaZE0(k,ll,lll) = min( 1.D0, ( ZEminel(k,ll)
     &                    - ZEtaylor(k,ll,1) ) / ( ZEtaylorvert(k,ll,lll) - ZEtaylor(k,ll,1) ) )
                     
                  endif

#ifdef TRACE
                  if ( iotataylorvert(k,ll,lll).gt.iotataylor(k,ll,1).and.
     &                 abs((iotataylorvert(k,ll,lll)-iotataylor(k,ll,1))).gt.1.0E-15 ) then 

                     alphaiota0(k,ll,lll) = min(1.D0,( iotamaxel(k,ll)
     &                    -iotataylor(k,ll,1) )/ (iotataylorvert(k,ll,lll) - iotataylor(k,ll,1)))
                     
                  elseif (iotataylorvert(k,ll,lll).eq.iotataylor(k,ll,1).or.
     &                    abs((iotataylorvert(k,ll,lll)-iotataylor(k,ll,1))).le.1.0E-15 ) then

                     alphaiota0(k,ll,lll) = 1.D0

                  elseif (iotataylorvert(k,ll,lll).lt.iotataylor(k,ll,1).
     &                    and.abs((iotataylorvert(k,ll,lll)-iotataylor(k,ll,1))).gt.1.0E-15) then

                     alphaiota0(k,ll,lll) = min( 1.D0,( iotaminel(k,ll)
     &                    -iotataylor(k,ll,1) )/( iotataylorvert(k,ll,lll)-iotataylor(k,ll,1)))
                     
                  endif
#endif
                  
#ifdef CHEM        
                  if ( iotataylorvert(k,ll,lll).gt.iotataylor(k,ll,1).and.
     &                 abs((iotataylorvert(k,ll,lll)-iotataylor(k,ll,1))).gt.1.0E-15 ) then 

                     alphaiota0(k,ll,lll) = min(1.D0,( iotamaxel(k,ll)-
     &                    iotataylor(k,ll,1) )/ (iotataylorvert(k,ll,lll) - iotataylor(k,ll,1)))
                     
                  elseif (iotataylorvert(k,ll,lll).eq.iotataylor(k,ll,1).or.
     &                    abs((iotataylorvert(k,ll,lll)-iotataylor(k,ll,1))).le.1.0E-15 ) then

                     alphaiota0(k,ll,lll) = 1.D0

                  elseif (iotataylorvert(k,ll,lll).lt.iotataylor(k,ll,1).and.
     &                    abs((iotataylorvert(k,ll,lll)-iotataylor(k,ll,1))).gt.1.0E-15) then

                     alphaiota0(k,ll,lll) = min( 1.D0,( iotaminel(k,ll)
     &                    -iotataylor(k,ll,1) )/( iotataylorvert(k,ll,lll)-iotataylor(k,ll,1)))
                     
                  endif

                  if ( iota2taylorvert(k,ll,lll).gt.iota2taylor(k,ll,1).and.
     &                 abs((iota2taylorvert(k,ll,lll)-iota2taylor(k,ll,1))).gt.1.0E-15 ) then 

                     alphaiota20(k,ll,lll) = min(1.D0,( iota2maxel(k,ll)
     &                    -iota2taylor(k,ll,1) )/ (iota2taylorvert(k,ll,lll) - iota2taylor(k,ll,1)))
                     
                  elseif (iota2taylorvert(k,ll,lll).eq.iota2taylor(k,ll,1).or.
     &                    abs((iota2taylorvert(k,ll,lll)-iota2taylor(k,ll,1))).le.1.0E-15 ) then

                     alphaiota20(k,ll,lll) = 1.D0

                  elseif (iota2taylorvert(k,ll,lll).lt.iota2taylor(k,ll,1).and.
     &                    abs((iota2taylorvert(k,ll,lll)-iota2taylor(k,ll,1))).gt.1.0E-15) then

                     alphaiota20(k,ll,lll) = min( 1.D0,( iota2minel(k,ll)
     &                    -iota2taylor(k,ll,1) )/( iota2taylorvert(k,ll,lll)-iota2taylor(k,ll,1)))
                     
                  endif
#endif

                  if ( QXtaylorvert(k,ll,lll).gt.QXtaylor(k,ll,1).and.
     &                 (QXtaylorvert(k,ll,lll)-QXtaylor(k,ll,1)).gt.1.0E-15 ) then !in xi1

                     alphaQX0(k,ll,lll) = min( 1.D0, ( QXmaxel(k,ll) 
     &                    - QXtaylor(k,ll,1) ) / ( QXtaylorvert(k,ll,lll) - QXtaylor(k,ll,1) ) )
                     
                  elseif ( QXtaylorvert(k,ll,lll).eq.QXtaylor(k,ll,1).or.
     &                    (QXtaylorvert(k,ll,lll)-QXtaylor(k,ll,1)).le.1.0E-15  ) then

                     alphaQX0(k,ll,lll) = 1.D0

                  elseif ( QXtaylorvert(k,ll,lll).lt.QXtaylor(k,ll,1).and.
     &                    (QXtaylorvert(k,ll,lll)-QXtaylor(k,ll,1)).gt.1.0E-15 ) then

                     alphaQX0(k,ll,lll) = min( 1.D0, ( QXminel(k,ll)
     &                    - QXtaylor(k,ll,1) ) / ( QXtaylorvert(k,ll,lll) - QXtaylor(k,ll,1) ) )

                  endif

                  if ( QYtaylorvert(k,ll,lll).gt.QYtaylor(k,ll,1).and.
     &                 (QYtaylorvert(k,ll,lll)-QYtaylor(k,ll,1)).gt.1.0E-15 ) then !in xi1

                     alphaQY0(k,ll,lll) = min( 1.D0, ( QYmaxel(k,ll) 
     &                    - QYtaylor(k,ll,1) ) / ( QYtaylorvert(k,ll,lll) - QYtaylor(k,ll,1) ) )
                     
                  elseif ( QYtaylorvert(k,ll,lll).eq.QYtaylor(k,ll,1).or.
     &                    (QYtaylorvert(k,ll,lll)-QYtaylor(k,ll,1)).le.1.0E-15  ) then

                     alphaQY0(k,ll,lll) = 1.D0

                  elseif ( QYtaylorvert(k,ll,lll).lt.QYtaylor(k,ll,1).and.
     &                    (QYtaylorvert(k,ll,lll)-QYtaylor(k,ll,1)).gt.1.0E-15 ) then

                     alphaQY0(k,ll,lll) = min( 1.D0, ( QYminel(k,ll) 
     &                    - QYtaylor(k,ll,1) ) / ( QYtaylorvert(k,ll,lll)  - QYtaylor(k,ll,1) ) )

                  endif

               enddo            !lll

            enddo               !ll

         elseif (dofs(k).eq.1) then

            cycle

         endif

      enddo

C.....Find the prescribed higher limiters by finding smallest local value

      alphaZE = 0.D0
      alphaQX = 0.D0
      alphaQY = 0.D0

#ifdef TRACE
      alphaiota = 0.D0
#endif

#ifdef CHEM
      alphaiota = 0.D0
      alphaiota2 = 0.D0
#endif

      do k=1,ne

         if (dofs(k).gt.1) then

            do ll=1,dofs(k)		

               alphaZE(k,ll) = minval( alphaZE0(k,ll,:) )
               alphaQX(k,ll) = minval( alphaQX0(k,ll,:) )
               alphaQY(k,ll) = minval( alphaQY0(k,ll,:) )

#ifdef TRACE
               alphaiota(k,ll) = minval( alphaiota0(k,ll,:) )
#endif

#ifdef CHEM
               alphaiota(k,ll) = minval( alphaiota0(k,ll,:) )
               alphaiota2(k,ll) = minval( alphaiota20(k,ll,:) )
#endif

            enddo

         elseif (dofs(k).eq.1) then

            cycle

         endif

      enddo

C.... Choose smallest (minimum) alpha for derivative in x or y

      alphaZEm = 0.D0
      alphaQXm = 0.D0
      alphaQYm = 0.D0

#ifdef TRACE
      alphaiotam = 0.D0
#endif

#ifdef CHEM
      alphaiotam = 0.D0
      alphaiota2m = 0.D0
#endif

      do k=1,ne

         if (dofs(k).gt.1) then
            
            do bb = 1,pdg_el(k)

               if( (bb+1)*(bb+2)/2.le.dofs(k) ) then

                  alphaZEm(k,(bb*(bb+1))/2 + 1) = 
     &                 minval( alphaZE(k,(bb*(bb-1))/2 +1:(bb*(bb+1)) / 2  ) )
                  alphaQXm(k,(bb*(bb+1))/2 + 1) = 
     &                 minval( alphaQX(k,(bb*(bb-1))/2 +1:(bb*(bb+1)) / 2  ) )
                  alphaQYm(k,(bb*(bb+1))/2 + 1) = 
     &                 minval( alphaQY(k,(bb*(bb-1))/2 +1:(bb*(bb+1)) / 2  ) )

#ifdef TRACE
                  alphaiotam(k,(bb*(bb+1))/2 + 1) = 
     &                 minval( alphaiota(k,(bb*(bb-1))/2 +1:(bb*(bb+1)) / 2 ) )
#endif

#ifdef CHEM
                  alphaiotam(k,(bb*(bb+1))/2 + 1) = 
     &                 minval( alphaiota(k,(bb*(bb-1))/2 +1:(bb*(bb+1)) / 2  ) )
                  alphaiota2m(k,(bb*(bb+1))/2 + 1) = 
     &                 minval( alphaiota2(k,(bb*(bb-1))/2 +1:(bb*(bb+1)) / 2 ) )
#endif

               endif
               
            enddo

         elseif (dofs(k).eq.1) then

            cycle

         endif

      enddo

C.....Use max higher derivative values for final limiter value

      alphaZE_max = 0.D0
      alphaQX_max = 0.D0
      alphaQY_max = 0.D0

#ifdef TRACE
      alphaiota_max = 0.D0
#endif

#ifdef CHEM
      alphaiota_max = 0.D0
      alphaiota2_max = 0.D0
#endif

      do k=1,ne

         if (dofs(k).gt.1) then

            do bb =1,pdg_el(k)

               if( (bb+1)*(bb+2)/2.le.dofs(k)) then

                  alphaZE_max(k,(bb*(bb+1))/2 + 1) = fd*maxval( 
     &                 alphaZEm(k,(bb*(bb+1))/2 + 1:dofs(k)) )
                  alphaQX_max(k,(bb*(bb+1))/2 + 1) = fd*maxval( 
     &                 alphaQXm(k,(bb*(bb+1))/2 + 1:dofs(k)) )
                  alphaQY_max(k,(bb*(bb+1))/2 + 1) = fd*maxval( 
     &                 alphaQYm(k,(bb*(bb+1))/2 + 1:dofs(k)) )

#ifdef TRACE
                  alphaiota_max(k,(bb*(bb+1))/2 + 1) = fd*maxval( 
     &                 alphaiotam(k,(bb*(bb+1))/2 + 1:dofs(k)) )
#endif

#ifdef CHEM
                  alphaiota_max(k,(bb*(bb+1))/2 + 1) = fd*maxval( 
     &                 alphaiotam(k,(bb*(bb+1))/2 + 1:dofs(k)) )
                  alphaiota2_max(k,(bb*(bb+1))/2 + 1) = fd*maxval( 
     &                 alphaiota2m(k,(bb*(bb+1))/2 + 1:dofs(k)) )
#endif

               endif

            enddo

         elseif (dofs(k).eq.1) then

            cycle

         endif

      enddo

C.....Limit on the Master element in the Taylor basis, via reconstruction 
C.....of unconstrained solutions with alpha constraints

      limitZE = 0.D0
      limitQX = 0.D0
      limitQY = 0.D0

      lim_count_roll = 0

#ifdef TRACE
      limitiota = 0.D0
#endif

#ifdef CHEM
      limitiota = 0.D0
      limitiota2 = 0.D0
#endif

      do k=1,ne

         lim_count = 0

         if (dofs(k).gt.1) then

            do ll=1,dofs(k)

               if ( ll.eq.1 ) then

                  limitZE(k,ll) = ZEtaylor(k,ll,1)
                  limitQX(k,ll) = QXtaylor(k,ll,1)
                  limitQY(k,ll) = QYtaylor(k,ll,1) 

#ifdef TRACE
                  limitiota(k,ll) = iotataylor(k,ll,1)
#endif

#ifdef CHEM
                  limitiota(k,ll) = iotataylor(k,ll,1)
                  limitiota2(k,ll) = iota2taylor(k,ll,1)
#endif

               elseif ( ll.ge.2 ) then
                  
                  do bb=1,pdg_el(k)

                     if ( ll.le.( ( (bb+1)*(bb+2)) / 2.D0 ).and.(ll.gt.
     &                    (bb*(bb+1)/2.D0) ) ) then

                        limitZE(k,ll) = alphaZE_max(k,(bb*(bb+1))/2 + 1) 
     &                       * ZEtaylor(k,ll,1)
                        limitQX(k,ll) = alphaQX_max(k,(bb*(bb+1))/2 + 1) 
     &                       * ZEtaylor(k,ll,1)
                        limitQY(k,ll) = alphaQY_max(k,(bb*(bb+1))/2 + 1) 
     &                       * ZEtaylor(k,ll,1)

#ifdef TRACE
                        limitiota(k,ll) = alphaiota_max(k,(bb*(bb+1))/2 + 1) 
     &                       * iotataylor(k,ll,1)
#endif

#ifdef CHEM
                        limitiota(k,ll) = alphaiota_max(k,(bb*(bb+1))/2 + 1) 
     &                       * iotataylor(k,ll,1)
                        limitiota2(k,ll) = alphaiota2_max(k,(bb*(bb+1))/2 + 1) 
     &                       * iota2taylor(k,ll,1)
#endif


c$$$  ! Make a counter to track limiting
c$$$  
c$$$  if ( ( alphaZE_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
c$$$  &                       alphaQX_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
c$$$  &                       alphaQY_max(k,(bb*(bb+1))/2 + 1).lt.1.D0 ).and.
c$$$  &                       chem_flag.ne.1.and.tracer_flag.ne.1 ) then
c$$$  
c$$$  lim_count = 1  
c$$$  
c$$$  elseif ( (alphaZE_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
c$$$  &                          alphaQX_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
c$$$  &                          alphaQY_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
c$$$  &                          alphaiota_max(k,(bb*(bb+1))/2 + 1).lt.1.D0 ).and.
c$$$  &                          chem_flag.ne.1.and.tracer_flag.eq.1 ) then
c$$$  
c$$$  lim_count = 1   
c$$$  
c$$$  elseif ( (alphaZE_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
c$$$  &                          alphaQX_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
c$$$  &                          alphaQY_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
c$$$  &                          alphaiota_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
c$$$  &                          alphaiota2_max(k,(bb*(bb+1))/2 + 1).lt.1.D0 ).and.
c$$$  &                          chem_flag.eq.1 ) then
c$$$  
c$$$  lim_count = 1 
c$$$  
c$$$  endif
                        
                     endif
                     
                  enddo
                  
               endif

            enddo

         elseif (dofs(k).eq.1) then

            cycle

         endif

                                !lim_count_roll = lim_count_roll + lim_count

      enddo

C.....Transform back to the Dubiner basis (multiply by NmatrixInv),

      ZEconst =  0.D0
      QXconst =  0.D0
      QYconst =  0.D0

#ifdef TRACE
      iotaconst =  0.D0
#endif

#ifdef CHEM
      iotaconst =  0.D0
      iota2const =  0.D0
#endif

      do k=1,ne

         if (dofs(k).gt.1) then
                                !do lll=1,3

            do ll=1,dofs(k)

               do ss=1,dofs(k)

                  ZEconst(k,ll) = ZEconst(k,ll) + NmatrixInv(k,ll,ss,dofs(k)) 
     &                 * limitZE(k,ss)
                  QXconst(k,ll) = QXconst(k,ll) + NmatrixInv(k,ll,ss,dofs(k)) 
     &                 * limitQX(k,ss)
                  QYconst(k,ll) = QYconst(k,ll) + NmatrixInv(k,ll,ss,dofs(k)) 
     &                 * limitQY(k,ss)

#ifdef TRACE
                  iotaconst(k,ll) = iotaconst(k,ll) + 
     &                 NmatrixInv(k,ll,ss,dofs(k)) * limitiota(k,ss)
#endif

#ifdef CHEM
                  iotaconst(k,ll) = iotaconst(k,ll) + 
     &                 NmatrixInv(k,ll,ss,dofs(k)) * limitiota(k,ss)
                  iota2const(k,ll) = iota2const(k,ll) + 
     &                 NmatrixInv(k,ll,ss,dofs(k)) * limitiota2(k,ss)
#endif


               enddo

            enddo

         elseif (dofs(k).eq.1) then

            cycle

         endif

      enddo

C.....Set limit values

      do k=1,ne

         if (dofs(k).gt.1) then

            do ll = 1,dofs(k)

               ZE(ll,k,irk+1) = ZEconst(k,ll)
               QX(ll,k,irk+1) = QXconst(k,ll)
               QY(ll,k,irk+1) = QYconst(k,ll)

#ifdef TRACE
               iota(ll,k,irk+1) = iotaconst(k,ll)
#endif

#ifdef CHEM
               iota(ll,k,irk+1) = iotaconst(k,ll)
               iota2(ll,k,irk+1) = iota2const(k,ll)
#endif

            enddo

         elseif (dofs(k).eq.1) then

            cycle

         endif

      enddo
      
      return
      end subroutine

C***********************************************************************
C     
C     SUBROUTINE SLOPELIMITER10()
C     
C     Written 2010
C     
C     This subroutine selects the first *adapted* vertex limiter based on
C     a Taylor Polynomial basis, and is consistent with p_adaptation
C     and works for any p
C     
C     -- cem
C     
C***********************************************************************

      SUBROUTINE SLOPELIMITER10()

C.....Use appropriate modules

      USE GLOBAL
      USE DG

#ifdef CMPI
      USE MESSENGER
      USE MESSENGER_ELEM
#endif

      IMPLICIT NONE

      Integer k,ll,ss,lll,ell,bb,maxneigh

      REAL(SZ) fd,marea,fde

C.....We work over the master element
C.....Set initial values

      fd = slope_weight         ! reduces diffusion fd = 1 => full diffusion
      fde = fd                  ! add weight for lower order pieces (fd<1 => stronger limiting)     

      DO k=1,mne

         if (dofs(k).gt.1) then

            DO ll = 1,dofs(k)

               ZE(ll,k,NRK+2) = ZE(ll,k,IRK+1)
               QX(ll,k,NRK+2) = QX(ll,k,IRK+1)
               QY(ll,k,NRK+2) = QY(ll,k,IRK+1)

#ifdef TRACE
               iota(ll,k,NRK+2) = iota(ll,k,IRK+1)
#endif

#ifdef CHEM
               iota(ll,k,NRK+2) = iota(ll,k,IRK+1)
               iota2(ll,k,NRK+2) = iota2(ll,k,IRK+1)
#endif


            ENDDO

         elseif (dofs(k).eq.1) then

            cycle

         endif

      ENDDO   

      marea = 2.D0              !master elements area

C.....Convert initial values to the Taylor basis (multiply by Nmatrix) on base element


      ZEtaylor = 0.D0 
      QXtaylor = 0.D0
      QYtaylor = 0.D0

#ifdef TRACE
      iotataylor = 0.D0
#endif

#ifdef CHEM
      iotataylor = 0.D0
      iota2taylor = 0.D0
#endif

      do k=1,mne

         if (dofs(k).gt.1) then

            do ll=1,dofs(k)
               
               do ss=1,dofs(k)

                  ZEtaylor(k,ll,1) = ZEtaylor(k,ll,1) + Nmatrix(k,ll,ss,dofs(k)) * ZE(ss,k,nrk+2)
                  QXtaylor(k,ll,1) = QXtaylor(k,ll,1) + Nmatrix(k,ll,ss,dofs(k))* QX(ss,k,nrk+2)
                  QYtaylor(k,ll,1) = QYtaylor(k,ll,1) + Nmatrix(k,ll,ss,dofs(k)) * QY(ss,k,nrk+2)

#ifdef TRACE
                  iotataylor(k,ll,1) = iotataylor(k,ll,1) + Nmatrix(k,ll,ss,dofs(k)) * iota(ss,k,nrk+2)
#endif

#ifdef CHEM
                  iotataylor(k,ll,1) = iotataylor(k,ll,1) + Nmatrix(k,ll,ss,dofs(k)) * iota(ss,k,nrk+2)
                  iota2taylor(k,ll,1) = iota2taylor(k,ll,1) + Nmatrix(k,ll,ss,dofs(k)) * iota2(ss,k,nrk+2)
#endif
                  
               enddo

            enddo

         elseif (dofs(k).eq.1) then

            cycle

         endif

      enddo

C.....Find values at vertices of base elements and neighbors


      ZEmax = -100.D0
      QXmax = -100.D0
      QYmax = -100.D0
      ZEmin = 100.D0
      QXmin = 100.D0
      QYmin = 100.D0

#ifdef TRACE
      iotamax = -100.D0
      iotamin = 100.D0
#endif

#ifdef CHEM
      iotamax = -100.D0
      iota2max = -100.D0
      iotamin = 100.D0
      iota2min = 100.D0
#endif


#ifdef CMPI

      CALL UPDATER_ELEM_MOD2(ZEtaylor,QXtaylor,QYtaylor,1,3)

#ifdef TRACE
      CALL UPDATER_ELEM_MOD2(iotataylor,QXtaylor,QYtaylor,1,2 )  
#endif

#ifdef CHEM
      CALL UPDATER_ELEM_MOD2(iotataylor,iota2taylor,QYtaylor,1,2 )  
#endif

      CALL UPDATERV(ZEmin)
      CALL UPDATERV(ZEmax)
      CALL UPDATERV(QXmin)
      CALL UPDATERV(QXmax)
      CALL UPDATERV(QYmin)
      CALL UPDATERV(QYmax)

#ifdef TRACE
         CALL UPDATERV(iotamax)
         CALL UPDATERV(iotamin)
#endif

#ifdef CHEM
         CALL UPDATERV(iotamax)
         CALL UPDATERV(iotamin)
         CALL UPDATERV(iota2max)
         CALL UPDATERV(iota2min)
#endif

#endif

      do ell=1,mnp

         do ll=1,minval(dofs(neigh_elem(ell,1:nneigh_elem(ell))))

C.....Find max and min values over polynomial coefficients

            ZEmax(ell,ll) = max(maxval( ZEtaylor(NEIGH_ELEM(ell,1:nneigh_elem(ell)),ll,1) )
     &           , ZEmax(ell,ll))
            QXmax(ell,ll) = max(maxval( QXtaylor(NEIGH_ELEM(ell,1:nneigh_elem(ell)),ll,1) )
     &           , QXmax(ell,ll))
            QYmax(ell,ll) = max(maxval( QYtaylor(NEIGH_ELEM(ell,1:nneigh_elem(ell)),ll,1) )
     &           , QYmax(ell,ll))
            ZEmin(ell,ll) = min(minval( ZEtaylor(NEIGH_ELEM(ell,1:nneigh_elem(ell)),ll,1) )
     &           , ZEmin(ell,ll))
            QXmin(ell,ll) = min(minval( QXtaylor(NEIGH_ELEM(ell,1:nneigh_elem(ell)),ll,1) )
     &           , QXmin(ell,ll))
            QYmin(ell,ll) = min(minval( QYtaylor(NEIGH_ELEM(ell,1:nneigh_elem(ell)),ll,1) )
     &           , QYmin(ell,ll))

#ifdef TRACE
            iotamax(ell,ll) = max(maxval( iotataylor(NEIGH_ELEM(ell,1:nneigh_elem(ell)),ll,1) )
     &           , iotamax(ell,ll))
            iotamin(ell,ll) = min(minval( iotataylor(NEIGH_ELEM(ell,1:nneigh_elem(ell)),ll,1) )
     &           , iotamin(ell,ll))
#endif

#ifdef CHEM
            iotamax(ell,ll) = max(maxval( iotataylor(NEIGH_ELEM(ell,1:nneigh_elem(ell)),ll,1) )
     &           , iotamax(ell,ll))
            iota2max(ell,ll) = max(maxval( iota2taylor(NEIGH_ELEM(ell,1:nneigh_elem(ell)),ll,1) )
     &           , iota2max(ell,ll))
            iotamin(ell,ll) = min(minval( iotataylor(NEIGH_ELEM(ell,1:nneigh_elem(ell)),ll,1) )
     &           , iotamin(ell,ll))
            iota2min(ell,ll) = min(minval( iota2taylor(NEIGH_ELEM(ell,1:nneigh_elem(ell)),ll,1) )
     &           , iota2min(ell,ll))
#endif
            
         enddo

      enddo


#ifdef CMPI

      CALL UPDATERV(ZEmin)
      CALL UPDATERV(ZEmax)
      CALL UPDATERV(QXmin)
      CALL UPDATERV(QXmax)
      CALL UPDATERV(QYmin)
      CALL UPDATERV(QYmax)

#ifdef TRACE
      CALL UPDATERV(iotamax)
      CALL UPDATERV(iotamin)
#endif

#ifdef CHEM
      CALL UPDATERV(iotamax)
      CALL UPDATERV(iotamin)
      CALL UPDATERV(iota2max)
      CALL UPDATERV(iota2min)
#endif

#endif


C.....Must generate linear recostructions at vertices

      ZEtaylorvert = 0.D0
      QXtaylorvert = 0.D0
      Qytaylorvert = 0.D0

#ifdef TRACE
      iotataylorvert = 0.D0
#endif

#ifdef CHEM
      iotataylorvert = 0.D0
      iota2taylorvert = 0.D0
#endif

      do k=1,mne

         if (dofs(k).gt.1) then

            do ll=1,dofs(k)

               do lll=1,3

                  if (ll.eq.1) then

                     ZEtaylorvert(k,ll,lll) = ZEtaylor(k,ll,1) + ZEtaylor(k,ll+1,1)*( xi2vert(k,lll) -xi2BCb(k) )
     &                    + ZEtaylor(k,ll+2,1)*( xi1vert(k,lll) - xi1BCb(k) )
                     QXtaylorvert(k,ll,lll) = QXtaylor(k,ll,1) + QXtaylor(k,ll+1,1)*( xi2vert(k,lll) -xi2BCb(k) )
     &                    + QXtaylor(k,ll+2,1)*( xi1vert(k,lll) - xi1BCb(k) )
                     Qytaylorvert(k,ll,lll) = QYtaylor(k,ll,1) + QYtaylor(k,ll+1,1)*( xi2vert(k,lll) -xi2BCb(k) ) 
     &                    + QYtaylor(k,ll+2,1)*( xi1vert(k,lll) - xi1BCb(k) )

#ifdef TRACE
                     iotataylorvert(k,ll,lll) = iotataylor(k,ll,1) + 
     &                    iotataylor(k,ll+1,1)*( xi2vert(k,lll) - xi2BCb(k) )
     &                    + iotataylor(k,ll+2,1)*( xi1vert(k,lll) - xi1BCb(k) )
#endif

#ifdef CHEM
                     iotataylorvert(k,ll,lll) = iotataylor(k,ll,1) + 
     &                    iotataylor(k,ll+1,1)*( xi2vert(k,lll) - xi2BCb(k) )
     &                    + iotataylor(k,ll+2,1)*( xi1vert(k,lll) - xi1BCb(k) )
                     iota2taylorvert(k,ll,lll) = iota2taylor(k,ll,1) + 
     &                    iota2taylor(k,ll+1,1)*( xi2vert(k,lll) - xi2BCb(k) )
     &                    + iota2taylor(k,ll+2,1)*( xi1vert(k,lll) - xi1BCb(k) )
#endif

                  elseif (ll.gt.1) then

                     ZEtaylorvert(k,ll,lll) = ZEtaylor(k,ll,1) + 
     &                    ZEtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) ),1)*( xi2vert(k,lll) - xi2BCb(k) ) 
     &                    + ZEtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) )+1,1)*( xi1vert(k,lll) - xi1BCb(k) )
                     QXtaylorvert(k,ll,lll) = QXtaylor(k,ll,1) + 
     &                    QXtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) ),1)*( xi2vert(k,lll) - xi2BCb(k) ) 
     &                    + QXtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) )+1,1)*( xi1vert(k,lll) - xi1BCb(k) )
                     Qytaylorvert(k,ll,lll) = QYtaylor(k,ll,1) + 
     &                    QYtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) ),1)*( xi2vert(k,lll) - xi2BCb(k) ) 
     &                    + QYtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) )+1,1)*( xi1vert(k,lll) - xi1BCb(k) )

#ifdef TRACE
                     iotataylorvert(k,ll,lll) = iotataylor(k,ll,1)+
     &                    iotataylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) ),1)*( xi2vert(k,lll) - xi2BCb(k))
     &                    + iotataylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) )+1,1)*( xi1vert(k,lll) - xi1BCb(k) )
#endif

#ifdef CHEM
                     iotataylorvert(k,ll,lll) = iotataylor(k,ll,1)+
     &                    iotataylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) ),1)*( xi2vert(k,lll) - xi2BCb(k))
     &                    + iotataylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) )+1,1)*( xi1vert(k,lll) - xi1BCb(k) )
                     iota2taylorvert(k,ll,lll) = iota2taylor(k,ll,1)+
     &                    iotataylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) ),1)*( xi2vert(k,lll) - xi2BCb(k))
     &                    + iota2taylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) )+1,1)*( xi1vert(k,lll) - xi1BCb(k) )
#endif

                  else

                  endif

               enddo
               
            enddo
            
         elseif (dofs(k).eq.1) then

            cycle

         endif

      enddo

C.....Compute alphas for each variable in each order derivitive


      alphaZE0 = 0.D0
      alphaQX0 = 0.D0
      alphaQY0 = 0.D0

#ifdef TRACE
      alphaiota0 = 0.D0
#endif

#ifdef CHEM
      alphaiota0 = 0.D0
      alphaiota20 = 0.D0
#endif

      do k=1,mne

         if (dofs(k).gt.1) then

            do lll=1,3
               
               do ll=1,dofs(k)

                  if (ZEmin(nm(k,lll),ll).ne.ZEmax(nm(k,lll),ll)) then

                     if ( ( ZEtaylorvert(k,ll,lll).gt.ZEtaylor(k,ll,1) ).and.
     &                    ( abs(ZEtaylorvert(k,ll,lll)-ZEtaylor(k,ll,1)).gt.1.0E-15 ).and.
     &                    ( ZEmax(nm(k,lll),ll).ne.ZEtaylor(k,ll,1) ) ) then  

                        alphaZE0(k,ll,lll) = min(1.D0,  ( ZEmax(nm(k,lll),ll)
     &                       - ZEtaylor(k,ll,1) )/ (ZEtaylorvert(k,ll,lll) - ZEtaylor(k,ll,1)))

                                !adapted part

                     elseif ( ( ZEtaylorvert(k,ll,lll).gt.ZEtaylor(k,ll,1) ).and.
     &                       ( abs(ZEtaylorvert(k,ll,lll)-ZEtaylor(k,ll,1)).gt.1.0E-15 ).and.
     &                       ( ZEmax(nm(k,lll),ll).eq.ZEtaylor(k,ll,1) ) ) then 

                        alphaZE0(k,ll,lll) = min(fd, abs( ( ZEmax(nm(k,lll),ll)
     &                       - ZEmin(nm(k,lll),ll) )/(ZEtaylorvert(k,ll,lll) - ZEmax(nm(k,lll),ll)) ) ) 
                        

                     elseif ( (ZEtaylorvert(k,ll,lll).lt.ZEtaylor(k,ll,1) )
     &                       .and.( abs(ZEtaylorvert(k,ll,lll)-ZEtaylor(k,ll,1)).gt.1.0E-15 ).and.
     &                       ( ZEmin(nm(k,lll),ll).ne.ZEtaylor(k,ll,1) ) ) then 

                        alphaZE0(k,ll,lll) = min( 1.D0,( ZEmin(nm(k,lll),ll)
     &                       - ZEtaylor(k,ll,1) )/( ZEtaylorvert(k,ll,lll)-ZEtaylor(k,ll,1)))

                                !adapted part

                     elseif ( (ZEtaylorvert(k,ll,lll).lt.ZEtaylor(k,ll,1) )
     &                       .and.( abs(ZEtaylorvert(k,ll,lll)-ZEtaylor(k,ll,1)).gt.1.0E-15 ).and.
     &                       ( ZEmin(nm(k,lll),ll).eq.ZEtaylor(k,ll,1) ) ) then 

                        alphaZE0(k,ll,lll) = min( fd,abs( ( ZEmin(nm(k,lll),ll)
     &                       - ZEmax(nm(k,lll),ll) )/( ZEtaylorvert(k,ll,lll)- ZEmin(nm(k,lll),ll))) )
                        
                        
                     elseif ( ( ZEtaylorvert(k,ll,lll).eq.ZEtaylor(k,ll,1) ).or.
     &                       ( abs(ZEtaylorvert(k,ll,lll)-ZEtaylor(k,ll,1)).le.1.0E-15 ) ) then

                        alphaZE0(k,ll,lll) = 1.D0

                     endif

                  else

                     alphaZE0(k,ll,lll) = 1.D0

                  endif
#ifdef TRACE
                  if (iotamin(nm(k,lll),ll).ne.iotamax(nm(k,lll),ll)) then

                     if ( ( iotataylorvert(k,ll,lll).gt.iotataylor(k,ll,1) ).and.
     &                    ( abs(iotataylorvert(k,ll,lll)-iotataylor(k,ll,1)).gt.1.0E-15 ).and.
     &                    ( iotamax(nm(k,lll),ll).ne.iotataylor(k,ll,1) ) ) then  

                        alphaiota0(k,ll,lll) = min(1.D0,  ( iotamax(nm(k,lll),ll)
     &                       - iotataylor(k,ll,1) )/ (iotataylorvert(k,ll,lll) - iotataylor(k,ll,1)))

                                !adapted part

                     elseif ( ( iotataylorvert(k,ll,lll).gt.iotataylor(k,ll,1) ).and.
     &                       ( abs(iotataylorvert(k,ll,lll)-iotataylor(k,ll,1)).gt.1.0E-15 ).and.
     &                       ( iotamax(nm(k,lll),ll).eq.iotataylor(k,ll,1) ) ) then 

                        alphaiota0(k,ll,lll) = min(fd, abs( ( iotamax(nm(k,lll),ll)
     &                       - iotamin(nm(k,lll),ll) )/(iotataylorvert(k,ll,lll) - iotamax(nm(k,lll),ll)) ) ) 
                        

                     elseif ( (iotataylorvert(k,ll,lll).lt.iotataylor(k,ll,1) )
     &                       .and.( abs(iotataylorvert(k,ll,lll)-iotataylor(k,ll,1)).gt.1.0E-15 ).and.
     &                       ( iotamin(nm(k,lll),ll).ne.iotataylor(k,ll,1) ) ) then 

                        alphaiota0(k,ll,lll) = min( 1.D0,( iotamin(nm(k,lll),ll)
     &                       - iotataylor(k,ll,1) )/( iotataylorvert(k,ll,lll)-iotataylor(k,ll,1)))

                                !adapted part

                     elseif ( (iotataylorvert(k,ll,lll).lt.iotataylor(k,ll,1) )
     &                       .and.( abs(iotataylorvert(k,ll,lll)-iotataylor(k,ll,1)).gt.1.0E-15 ).and.
     &                       ( iotamin(nm(k,lll),ll).eq.iotataylor(k,ll,1) ) ) then 

                        alphaiota0(k,ll,lll) = min( fd,abs( ( iotamin(nm(k,lll),ll)
     &                       - iotamax(nm(k,lll),ll) )/( iotataylorvert(k,ll,lll)- iotamin(nm(k,lll),ll))) )
                        
                        
                     elseif ( ( iotataylorvert(k,ll,lll).eq.iotataylor(k,ll,1) ).or.
     &                       ( abs(iotataylorvert(k,ll,lll)-iotataylor(k,ll,1)).le.1.0E-15 ) ) then

                        alphaiota0(k,ll,lll) = 1.D0

                     endif

                  else

                     alphaiota0(k,ll,lll) = 1.D0

                  endif
#endif
                  
#ifdef CHEM
                  if (iotamin(nm(k,lll),ll).ne.iotamax(nm(k,lll),ll)) then

                     if ( ( iotataylorvert(k,ll,lll).gt.iotataylor(k,ll,1) ).and.
     &                    ( abs(iotataylorvert(k,ll,lll)-iotataylor(k,ll,1)).gt.1.0E-15 ).and.
     &                    ( iotamax(nm(k,lll),ll).ne.iotataylor(k,ll,1) ) ) then  

                        alphaiota0(k,ll,lll) = min(1.D0,  ( iotamax(nm(k,lll),ll)
     &                       - iotataylor(k,ll,1) )/ (iotataylorvert(k,ll,lll) - iotataylor(k,ll,1)))

                                !adapted part

                     elseif ( ( iotataylorvert(k,ll,lll).gt.iotataylor(k,ll,1) ).and.
     &                       ( abs(iotataylorvert(k,ll,lll)-iotataylor(k,ll,1)).gt.1.0E-15 ).and.
     &                       ( iotamax(nm(k,lll),ll).eq.iotataylor(k,ll,1) ) ) then 

                        alphaiota0(k,ll,lll) = min(fd, abs( ( iotamax(nm(k,lll),ll)
     &                       - iotamin(nm(k,lll),ll) )/(iotataylorvert(k,ll,lll) - iotamax(nm(k,lll),ll)) ) ) 
                        

                     elseif ( (iotataylorvert(k,ll,lll).lt.iotataylor(k,ll,1) )
     &                       .and.( abs(iotataylorvert(k,ll,lll)-iotataylor(k,ll,1)).gt.1.0E-15 ).and.
     &                       ( iotamin(nm(k,lll),ll).ne.iotataylor(k,ll,1) ) ) then 

                        alphaiota0(k,ll,lll) = min( 1.D0,( iotamin(nm(k,lll),ll)
     &                       - iotataylor(k,ll,1) )/( iotataylorvert(k,ll,lll)-iotataylor(k,ll,1)))

                                !adapted part

                     elseif ( (iotataylorvert(k,ll,lll).lt.iotataylor(k,ll,1) )
     &                       .and.( abs(iotataylorvert(k,ll,lll)-iotataylor(k,ll,1)).gt.1.0E-15 ).and.
     &                       ( iotamin(nm(k,lll),ll).eq.iotataylor(k,ll,1) ) ) then 

                        alphaiota0(k,ll,lll) = min( fd,abs( ( iotamin(nm(k,lll),ll)
     &                       - iotamax(nm(k,lll),ll) )/( iotataylorvert(k,ll,lll)- iotamin(nm(k,lll),ll))) )
                        
                        
                     elseif ( ( iotataylorvert(k,ll,lll).eq.iotataylor(k,ll,1) ).or.
     &                       ( abs(iotataylorvert(k,ll,lll)-iotataylor(k,ll,1)).le.1.0E-15 ) ) then

                        alphaiota0(k,ll,lll) = 1.D0

                     endif

                  else

                     alphaiota0(k,ll,lll) = 1.D0

                  endif

                  if (iota2min(nm(k,lll),ll).ne.iota2max(nm(k,lll),ll)) then

                     if ( ( iota2taylorvert(k,ll,lll).gt.iota2taylor(k,ll,1) ).and.
     &                    ( abs(iota2taylorvert(k,ll,lll)-iota2taylor(k,ll,1)).gt.1.0E-15 ).and.
     &                    ( iota2max(nm(k,lll),ll).ne.iota2taylor(k,ll,1) ) ) then  

                        alphaiota20(k,ll,lll) = min(1.D0,  ( iota2max(nm(k,lll),ll)
     &                       - iota2taylor(k,ll,1) )/ (iota2taylorvert(k,ll,lll) - iota2taylor(k,ll,1)))

                                !adapted part

                     elseif ( ( iota2taylorvert(k,ll,lll).gt.iota2taylor(k,ll,1) ).and.
     &                       ( abs(iota2taylorvert(k,ll,lll)-iota2taylor(k,ll,1)).gt.1.0E-15 ).and.
     &                       ( iota2max(nm(k,lll),ll).eq.iota2taylor(k,ll,1) ) ) then 

                        alphaiota20(k,ll,lll) = min(fd, abs( ( iota2max(nm(k,lll),ll)
     &                       - iota2min(nm(k,lll),ll) )/(iota2taylorvert(k,ll,lll) - iota2max(nm(k,lll),ll)) ) ) 
                        

                     elseif ( (iota2taylorvert(k,ll,lll).lt.iota2taylor(k,ll,1) )
     &                       .and.( abs(iota2taylorvert(k,ll,lll)-iota2taylor(k,ll,1)).gt.1.0E-15 ).and.
     &                       ( iota2min(nm(k,lll),ll).ne.iota2taylor(k,ll,1) ) ) then 

                        alphaiota20(k,ll,lll) = min( 1.D0,( iota2min(nm(k,lll),ll)
     &                       - iota2taylor(k,ll,1) )/( iota2taylorvert(k,ll,lll)-iota2taylor(k,ll,1)))

                                !adapted part

                     elseif ( (iota2taylorvert(k,ll,lll).lt.iota2taylor(k,ll,1) )
     &                       .and.( abs(iota2taylorvert(k,ll,lll)-iota2taylor(k,ll,1)).gt.1.0E-15 ).and.
     &                       ( iota2min(nm(k,lll),ll).eq.iota2taylor(k,ll,1) ) ) then 

                        alphaiota20(k,ll,lll) = min( fd,abs( ( iota2min(nm(k,lll),ll)
     &                       - iota2max(nm(k,lll),ll) )/( iota2taylorvert(k,ll,lll)- iota2min(nm(k,lll),ll))) )
                        
                        
                     elseif ( ( iota2taylorvert(k,ll,lll).eq.iota2taylor(k,ll,1) ).or.
     &                       ( abs(iota2taylorvert(k,ll,lll)-iota2taylor(k,ll,1)).le.1.0E-15 ) ) then

                        alphaiota20(k,ll,lll) = 1.D0

                     endif

                  else

                     alphaiota20(k,ll,lll) = 1.D0

                  endif
#endif

                  if (QXmin(nm(k,lll),ll).ne.QXmax(nm(k,lll),ll)) then

                     if ( ( QXtaylorvert(k,ll,lll).gt.QXtaylor(k,ll,1) ).and.
     &                    ( abs(QXtaylorvert(k,ll,lll)-QXtaylor(k,ll,1)).gt.1.0E-15 ).and.
     &                    ( QXmax(nm(k,lll),ll).ne.QXtaylor(k,ll,1) ) ) then  

                        alphaQX0(k,ll,lll) = min(1.D0,  ( QXmax(nm(k,lll),ll)
     &                       - QXtaylor(k,ll,1) )/ (QXtaylorvert(k,ll,lll) - QXtaylor(k,ll,1)))

                                !adapted part

                     elseif ( ( QXtaylorvert(k,ll,lll).gt.QXtaylor(k,ll,1) ).and.
     &                       ( abs(QXtaylorvert(k,ll,lll)-QXtaylor(k,ll,1)).gt.1.0E-15 ).and.
     &                       ( QXmax(nm(k,lll),ll).eq.QXtaylor(k,ll,1) ) ) then 

                        alphaQX0(k,ll,lll) = min(fd, abs( ( QXmax(nm(k,lll),ll)
     &                       - QXmin(nm(k,lll),ll) )/(QXtaylorvert(k,ll,lll) - QXmax(nm(k,lll),ll)) ) ) 
                        

                     elseif ( (QXtaylorvert(k,ll,lll).lt.QXtaylor(k,ll,1) )
     &                       .and.( abs(QXtaylorvert(k,ll,lll)-QXtaylor(k,ll,1)).gt.1.0E-15 ).and.
     &                       ( QXmin(nm(k,lll),ll).ne.QXtaylor(k,ll,1) ) ) then 

                        alphaQX0(k,ll,lll) = min( 1.D0,( QXmin(nm(k,lll),ll)
     &                       - QXtaylor(k,ll,1) )/( QXtaylorvert(k,ll,lll)-QXtaylor(k,ll,1)))

                                !adapted part

                     elseif ( (QXtaylorvert(k,ll,lll).lt.QXtaylor(k,ll,1) )
     &                       .and.( abs(QXtaylorvert(k,ll,lll)-QXtaylor(k,ll,1)).gt.1.0E-15 ).and.
     &                       ( QXmin(nm(k,lll),ll).eq.QXtaylor(k,ll,1) ) ) then 

                        alphaQX0(k,ll,lll) = min( fd,abs( ( QXmin(nm(k,lll),ll)
     &                       - QXmax(nm(k,lll),ll) )/( QXtaylorvert(k,ll,lll)- QXmin(nm(k,lll),ll))) )
                        
                        
                     elseif ( ( QXtaylorvert(k,ll,lll).eq.QXtaylor(k,ll,1) ).or.
     &                       ( abs(QXtaylorvert(k,ll,lll)-QXtaylor(k,ll,1)).le.1.0E-15 ) ) then

                        alphaQX0(k,ll,lll) = 1.D0

                     endif

                  else

                     alphaQX0(k,ll,lll) = 1.D0

                  endif


                  if (QYmin(nm(k,lll),ll).ne.QYmax(nm(k,lll),ll)) then

                     if ( ( QYtaylorvert(k,ll,lll).gt.QYtaylor(k,ll,1) ).and.
     &                    ( abs(QYtaylorvert(k,ll,lll)-QYtaylor(k,ll,1)).gt.1.0E-15 ).and.
     &                    ( QYmax(nm(k,lll),ll).ne.QYtaylor(k,ll,1) ) ) then  

                        alphaQY0(k,ll,lll) = min(1.D0,  ( QYmax(nm(k,lll),ll)
     &                       - QYtaylor(k,ll,1) )/ (QYtaylorvert(k,ll,lll) - QYtaylor(k,ll,1)))

                                !adapted part

                     elseif ( ( QYtaylorvert(k,ll,lll).gt.QYtaylor(k,ll,1) ).and.
     &                       ( abs(QYtaylorvert(k,ll,lll)-QYtaylor(k,ll,1)).gt.1.0E-15 ).and.
     &                       ( QYmax(nm(k,lll),ll).eq.QYtaylor(k,ll,1) ) ) then 

                        alphaQY0(k,ll,lll) = min(fd, abs( ( QYmax(nm(k,lll),ll)
     &                       - QYmin(nm(k,lll),ll) )/(QYtaylorvert(k,ll,lll) - QYmax(nm(k,lll),ll)) ) ) 
                        

                     elseif ( (QYtaylorvert(k,ll,lll).lt.QYtaylor(k,ll,1) )
     &                       .and.( abs(QYtaylorvert(k,ll,lll)-QYtaylor(k,ll,1)).gt.1.0E-15 ).and.
     &                       ( QYmin(nm(k,lll),ll).ne.QYtaylor(k,ll,1) ) ) then 

                        alphaQY0(k,ll,lll) = min( 1.D0,( QYmin(nm(k,lll),ll)
     &                       - QYtaylor(k,ll,1) )/( QYtaylorvert(k,ll,lll)-QYtaylor(k,ll,1)))

                                !adapted part

                     elseif ( (QYtaylorvert(k,ll,lll).lt.QYtaylor(k,ll,1) )
     &                       .and.( abs(QYtaylorvert(k,ll,lll)-QYtaylor(k,ll,1)).gt.1.0E-15 ).and.
     &                       ( QYmin(nm(k,lll),ll).eq.QYtaylor(k,ll,1) ) ) then 

                        alphaQY0(k,ll,lll) = min( fd,abs( ( QYmin(nm(k,lll),ll)
     &                       - QYmax(nm(k,lll),ll) )/( QYtaylorvert(k,ll,lll)- QYmin(nm(k,lll),ll))) )
                        
                        
                     elseif ( ( QYtaylorvert(k,ll,lll).eq.QYtaylor(k,ll,1) ).or.
     &                       ( abs(QYtaylorvert(k,ll,lll)-QYtaylor(k,ll,1)).le.1.0E-15 ) ) then

                        alphaQY0(k,ll,lll) = 1.D0

                     endif

                  else

                     alphaQY0(k,ll,lll) = 1.D0

                  endif

               enddo            !lll

            enddo               !ll

         elseif (dofs(k).eq.1) then

            cycle

         endif

      enddo

C.....Find the prescribed higher limiters by finding smallest local value

      alphaZE = 0.D0
      alphaQX = 0.D0
      alphaQY = 0.D0

#ifdef TRACE
      alphaiota = 0.D0
#endif

#ifdef CHEM
      alphaiota = 0.D0
      alphaiota2 = 0.D0
#endif

      do k=1,mne

         if (dofs(k).gt.1) then

            do ll=1,dofs(k)		

               alphaZE(k,ll) = minval( alphaZE0(k,ll,:) )
               alphaQX(k,ll) = minval( alphaQX0(k,ll,:) )
               alphaQY(k,ll) = minval( alphaQY0(k,ll,:) )

#ifdef TRACE
               alphaiota(k,ll) = minval( alphaiota0(k,ll,:) )
#endif

#ifdef CHEM
               alphaiota(k,ll) = minval( alphaiota0(k,ll,:) )
               alphaiota2(k,ll) = minval( alphaiota20(k,ll,:) )
#endif

            enddo

         elseif (dofs(k).eq.1) then

            cycle

         endif

      enddo

C.... Choose smallest (minimum) alpha for derivative in x or y

      alphaZEm = 0.D0
      alphaQXm = 0.D0
      alphaQYm = 0.D0

#ifdef TRACE
      alphaiotam = 0.D0
#endif

#ifdef CHEM
      alphaiotam = 0.D0
      alphaiota2m = 0.D0
#endif

      do k=1,mne

         if (dofs(k).gt.1) then
            
            do bb = 1,pdg_el(k)

               if( (bb+1)*(bb+2)/2.le.dofs(k) ) then

                  alphaZEm(k,(bb*(bb+1))/2 + 1) = 
     &                 minval( alphaZE(k,(bb*(bb-1))/2 +1:(bb*(bb+1)) / 2  ) )
                  alphaQXm(k,(bb*(bb+1))/2 + 1) = 
     &                 minval( alphaQX(k,(bb*(bb-1))/2 +1:(bb*(bb+1)) / 2  ) )
                  alphaQYm(k,(bb*(bb+1))/2 + 1) = 
     &                 minval( alphaQY(k,(bb*(bb-1))/2 +1:(bb*(bb+1)) / 2  ) )

#ifdef TRACE
                  alphaiotam(k,(bb*(bb+1))/2 + 1) = 
     &                 minval( alphaiota(k,(bb*(bb-1))/2 +1:(bb*(bb+1)) / 2 ) )
#endif

#ifdef CHEM
                  alphaiotam(k,(bb*(bb+1))/2 + 1) = 
     &                 minval( alphaiota(k,(bb*(bb-1))/2 +1:(bb*(bb+1)) / 2  ) )
                  alphaiota2m(k,(bb*(bb+1))/2 + 1) = 
     &                 minval( alphaiota2(k,(bb*(bb-1))/2 +1:(bb*(bb+1)) / 2 ) )
#endif

               endif
               
            enddo

         elseif (dofs(k).eq.1) then

            cycle

         endif

      enddo

C.....Use max higher derivative values for final limiter value

      alphaZE_max = 0.D0
      alphaQX_max = 0.D0
      alphaQY_max = 0.D0

#ifdef TRACE
      alphaiota_max = 0.D0
#endif

#ifdef CHEM
      alphaiota_max = 0.D0
      alphaiota2_max = 0.D0
#endif

      do k=1,mne

         if (dofs(k).gt.1) then

            do bb =1,pdg_el(k)

               if( (bb+1)*(bb+2)/2.le.dofs(k)) then

                  alphaZE_max(k,(bb*(bb+1))/2 + 1) = fde*maxval( 
     &                 alphaZEm(k,(bb*(bb+1))/2 + 1:dofs(k)) )
                  alphaQX_max(k,(bb*(bb+1))/2 + 1) = fde*maxval( 
     &                 alphaQXm(k,(bb*(bb+1))/2 + 1:dofs(k)) )
                  alphaQY_max(k,(bb*(bb+1))/2 + 1) = fde*maxval( 
     &                 alphaQYm(k,(bb*(bb+1))/2 + 1:dofs(k)) )

#ifdef TRACE
                  alphaiota_max(k,(bb*(bb+1))/2 + 1) = fde*maxval( 
     &                 alphaiotam(k,(bb*(bb+1)/2) +1 : dofs(k)))
#endif

#ifdef CHEM
                  alphaiota_max(k,(bb*(bb+1))/2 + 1) = fde*maxval( 
     &                 alphaiotam(k,(bb*(bb+1))/2 + 1:dofs(k)) )
                  alphaiota2_max(k,(bb*(bb+1))/2 + 1) = fde*maxval( 
     &                 alphaiota2m(k,(bb*(bb+1))/2 + 1:dofs(k)) )
#endif

               endif

            enddo

         elseif (dofs(k).eq.1) then

            cycle

         endif

      enddo

C.....Limit on the Master element in the Taylor basis, via reconstruction 
C.....of unconstrained solutions with alpha constraints


      limitZE = 0.D0
      limitQX = 0.D0
      limitQY = 0.D0

      lim_count_roll = 0

#ifdef TRACE
      limitiota = 0.D0
#endif

#ifdef CHEM
      limitiota = 0.D0
      limitiota2 = 0.D0
#endif

      do k=1,mne

         lim_count = 0

         if (dofs(k).gt.1) then

            do ll=1,dofs(k)

               if ( ll.eq.1 ) then

                  limitZE(k,ll) = ZEtaylor(k,ll,1)
                  limitQX(k,ll) = QXtaylor(k,ll,1)
                  limitQY(k,ll) = QYtaylor(k,ll,1) 

#ifdef TRACE
                  limitiota(k,ll) = iotataylor(k,ll,1)
#endif

#ifdef CHEM
                  limitiota(k,ll) = iotataylor(k,ll,1)
                  limitiota2(k,ll) = iota2taylor(k,ll,1)
#endif

               elseif ( ll.ge.2 ) then
                  
                  do bb=1,pdg_el(k)

                     if ( ll.le.( ( (bb+1)*(bb+2)) / 2.D0 ).and.(ll.gt.
     &                    (bb*(bb+1)/2.D0) ) ) then

                        limitZE(k,ll) = alphaZE_max(k,(bb*(bb+1))/2 + 1) 
     &                       * ZEtaylor(k,ll,1)
                        limitQX(k,ll) = alphaQX_max(k,(bb*(bb+1))/2 + 1) 
     &                       * ZEtaylor(k,ll,1)
                        limitQY(k,ll) = alphaQY_max(k,(bb*(bb+1))/2 + 1) 
     &                       * ZEtaylor(k,ll,1)


#ifdef TRACE
                        limitiota(k,ll) = alphaiota_max(k,(bb*(bb+1))/2 + 1) 
     &                       * iotataylor(k,ll,1)
#endif

#ifdef CHEM
                        limitiota(k,ll) = alphaiota_max(k,(bb*(bb+1))/2 + 1) 
     &                       * iotataylor(k,ll,1)
                        limitiota2(k,ll) = alphaiota2_max(k,(bb*(bb+1))/2 + 1) 
     &                       * iota2taylor(k,ll,1)
#endif


c$$$  ! Make a counter to track limiting
c$$$  
c$$$  if ( ( alphaZE_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
c$$$  &                       alphaQX_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
c$$$  &                       alphaQY_max(k,(bb*(bb+1))/2 + 1).lt.1.D0 ).and.
c$$$  &                       chem_flag.ne.1.and.tracer_flag.ne.1 ) then
c$$$  
c$$$  lim_count = 1  
c$$$  
c$$$  elseif ( (alphaZE_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
c$$$  &                          alphaQX_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
c$$$  &                          alphaQY_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
c$$$  &                          alphaiota_max(k,(bb*(bb+1))/2 + 1).lt.1.D0 ).and.
c$$$  &                          chem_flag.ne.1.and.tracer_flag.eq.1 ) then
c$$$  
c$$$  lim_count = 1   
c$$$  
c$$$  elseif ( (alphaZE_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
c$$$  &                          alphaQX_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
c$$$  &                          alphaQY_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
c$$$  &                          alphaiota_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
c$$$  &                          alphaiota2_max(k,(bb*(bb+1))/2 + 1).lt.1.D0 ).and.
c$$$  &                          chem_flag.eq.1 ) then
c$$$  
c$$$  lim_count = 1 
c$$$  
c$$$  endif
                        
                     endif
                     
                  enddo
                  
               endif

            enddo

         elseif (dofs(k).eq.1) then

            cycle

         endif

                                !lim_count_roll = lim_count_roll + lim_count

      enddo

C.....Transform back to the Dubiner basis (multiply by NmatrixInv),

      ZEconst =  0.D0
      QXconst =  0.D0
      QYconst =  0.D0

#ifdef TRACE
      iotaconst =  0.D0
#endif

#ifdef CHEM
      iotaconst =  0.D0
      iota2const =  0.D0
#endif

      do k=1,mne

         if (dofs(k).gt.1) then
                                !do lll=1,3

            do ll=1,dofs(k)

               do ss=1,dofs(k)

                  ZEconst(k,ll) = ZEconst(k,ll) + NmatrixInv(k,ll,ss,dofs(k)) 
     &                 * limitZE(k,ss)
                  QXconst(k,ll) = QXconst(k,ll) + NmatrixInv(k,ll,ss,dofs(k)) 
     &                 * limitQX(k,ss)
                  QYconst(k,ll) = QYconst(k,ll) + NmatrixInv(k,ll,ss,dofs(k)) 
     &                 * limitQY(k,ss)

#ifdef TRACE
                  iotaconst(k,ll) = iotaconst(k,ll) + 
     &                 NmatrixInv(k,ll,ss,dofs(k)) * limitiota(k,ss)
#endif

#ifdef CHEM
                  iotaconst(k,ll) = iotaconst(k,ll) + 
     &                 NmatrixInv(k,ll,ss,dofs(k)) * limitiota(k,ss)
                  iota2const(k,ll) = iota2const(k,ll) + 
     &                 NmatrixInv(k,ll,ss,dofs(k)) * limitiota2(k,ss)
#endif

               enddo

            enddo

         elseif (dofs(k).eq.1) then

            cycle

         endif

      enddo

C.....Set limit values

      do k=1,mne

         if (dofs(k).gt.1) then

            do ll = 1,dofs(k)

               ZE(ll,k,irk+1) = ZEconst(k,ll)
               QX(ll,k,irk+1) = QXconst(k,ll)
               QY(ll,k,irk+1) = QYconst(k,ll)

#ifdef TRACE
               iota(ll,k,irk+1) = iotaconst(k,ll)
#endif

#ifdef CHEM
               iota(ll,k,irk+1) = iotaconst(k,ll)
               iota2(ll,k,irk+1) = iota2const(k,ll)
#endif

            enddo

         elseif (dofs(k).eq.1) then

            cycle

         endif

      enddo
      
      return
      end subroutine

C******************************************************************************
C     
C     SUBROUTINE SLOPELIMITER9() -- Needs to be updated for ZE,QX,QY and iota2
C     
C     Written 2011
C     
C     This subroutine selects the first adapted Barth--Jespersen limiter over
C     a Taylor Polynomial basis, and is consistent with p-adaptation
C     to arbitrary order p
C     
C     - cem
C     
C******************************************************************************

      SUBROUTINE SLOPELIMITER9()

C.....Use appropriate modules

      USE GLOBAL
      USE DG

#ifdef CMPI
      USE MESSENGER_ELEM
#endif

      IMPLICIT NONE

      Integer k,ll,ss,lll,ell,bb

      REAL(SZ) fd,marea,fde
      Real(SZ), Allocatable :: ZEminel(:,:),ZEmaxel(:,:),QXminel(:,:),QXmaxel(:,:)
      Real(SZ), Allocatable :: QYminel(:,:),QYmaxel(:,:)
      Real(SZ), Allocatable :: iotaminel(:,:),iotamaxel(:,:)
      Real(SZ), Allocatable :: iota2minel(:,:),iota2maxel(:,:)

      Allocate ( ZEminel(mne,dofh),ZEmaxel(mne,dofh),QXminel(mne,dofh) )
      Allocate ( QYminel(mne,dofh),QYmaxel(mne,dofh),QXmaxel(mne,dofh) )
      Allocate ( iotaminel(mne,dofh),iotamaxel(mne,dofh) )
      Allocate ( iota2minel(mne,dofh),iota2maxel(mne,dofh) )


C.....We work over the master element
C.....Set initial values

      fd = slope_weight         ! add weight for lower order pieces (fd<1 => stronger limiting)
      

      DO k=1,NE

         if (dofs(k).gt.1) then

            DO ll = 1,dofs(k)

               ZE(ll,k,NRK+2) = ZE(ll,k,IRK+1)
               QX(ll,k,NRK+2) = QX(ll,k,IRK+1)
               QY(ll,k,NRK+2) = QY(ll,k,IRK+1)

#ifdef TRACE
               iota(ll,k,NRK+2) = iota(ll,k,IRK+1)
#endif

#ifdef CHEM
               iota(ll,k,NRK+2) = iota(ll,k,IRK+1)
               iota2(ll,k,NRK+2) = iota2(ll,k,IRK+1)
#endif


            ENDDO

         elseif (dofs(k).eq.1) then

            cycle

         endif

      ENDDO   

      marea = 2.D0              !master elements area

C.....Convert initial values to the Taylor basis (multiply by Nmatrix) on base element


      ZEtaylor = 0.D0 
      QXtaylor = 0.D0
      QYtaylor = 0.D0

#ifdef TRACE
      iotataylor = 0.D0
#endif

#ifdef CHEM
      iotataylor = 0.D0
      iota2taylor = 0.D0
#endif

      do k=1,ne

         if (dofs(k).gt.1) then

            do ll=1,dofs(k)
               
               do ss=1,dofs(k)

                  ZEtaylor(k,ll,1) = ZEtaylor(k,ll,1) + Nmatrix(k,ll,ss,dofs(k)) * ZE(ss,k,nrk+2)
                  QXtaylor(k,ll,1) = QXtaylor(k,ll,1) + Nmatrix(k,ll,ss,dofs(k)) * QX(ss,k,nrk+2)
                  QYtaylor(k,ll,1) = QYtaylor(k,ll,1) + Nmatrix(k,ll,ss,dofs(k)) * QY(ss,k,nrk+2)

#ifdef TRACE
                  iotataylor(k,ll,1) = iotataylor(k,ll,1) + Nmatrix(k,ll,ss,dofs(k))*iota(ss,k,nrk+2)
#endif

#ifdef CHEM
                  iotataylor(k,ll,1) = iotataylor(k,ll,1) + Nmatrix(k,ll,ss,dofs(k))*iota(ss,k,nrk+2)
                  iota2taylor(k,ll,1) = iota2taylor(k,ll,1) + Nmatrix(k,ll,ss,dofs(k))*iota2(ss,k,nrk+2)
#endif
                  
               enddo

            enddo

         elseif (dofs(k).eq.1) then

            cycle

         endif

      enddo

C.....Find values at vertices of base elements and neighbors


      ZEmaxel = -100.D0
      QXmaxel = -100.D0
      QYmaxel = -100.D0
      ZEminel = 100.D0
      QXminel = 100.D0
      QYminel = 100.D0

#ifdef TRACE
      iotamaxel = -100.D0
      iotaminel = 100.D0
#endif

#ifdef CHEM
      iotamaxel = -100.D0
      iota2maxel = -100.D0
      iotaminel = 100.D0
      iota2minel = 100.D0
#endif

#ifdef CMPI

      CALL UPDATER_ELEM_MOD2(ZEtaylor,QXtaylor,QYtaylor,1,3)

#ifdef TRACE
      CALL UPDATER_ELEM_MOD2(iotataylor,QXtaylor,QYtaylor,1,2 )  
#endif

#ifdef CHEM
      CALL UPDATER_ELEM_MOD2(iotataylor,iota2taylor,QYtaylor,1,2 )  
#endif

#endif

      do k = 1,ne

         do ll = 1,dofs(k) 

            do ell = 1,3        ! Number of edge neighbors for a triangle

C.....Find max and min values over polynomial coefficients
               
               ZEmaxel(k,ll) = max( ZEtaylor(k,ll,1),ZEtaylor(EL_NBORS(ell,k),ll,1), ZEmaxel(k,ll) )
               QXmaxel(k,ll) = max( QXtaylor(k,ll,1),QXtaylor(EL_NBORS(ell,k),ll,1), QXmaxel(k,ll) )
               QYmaxel(k,ll) = max( QYtaylor(k,ll,1),QYtaylor(EL_NBORS(ell,k),ll,1), QYmaxel(k,ll) )

               ZEminel(k,ll) = min( ZEtaylor(k,ll,1),ZEtaylor(EL_NBORS(ell,k),ll,1), ZEminel(k,ll) )
               QXminel(k,ll) = min( QXtaylor(k,ll,1),QXtaylor(EL_NBORS(ell,k),ll,1), QXminel(k,ll) )
               QYminel(k,ll) = min( QYtaylor(k,ll,1),QYtaylor(EL_NBORS(ell,k),ll,1), QYminel(k,ll) )

#ifdef TRACE
               iotamaxel(k,ll) = max( iotataylor(k,ll,1),iotataylor(EL_NBORS(ell,k),ll,1), 
     &              iotamaxel(k,ll) )
               iotaminel(k,ll) = min( iotataylor(k,ll,1),iotataylor(EL_NBORS(ell,k),ll,1), 
     &              iotaminel(k,ll) )
#endif

#ifdef CHEM
               iotamaxel(k,ll) = max( iotataylor(k,ll,1),iotataylor(EL_NBORS(ell,k),ll,1), 
     &              iotamaxel(k,ll) )
               iota2maxel(k,ll) = max( iota2taylor(k,ll,1),iota2taylor(EL_NBORS(ell,k),ll,1),
     &              iota2maxel(k,ll) )
               iotaminel(k,ll) = min( iotataylor(k,ll,1),iotataylor(EL_NBORS(ell,k),ll,1),
     &              iotaminel(k,ll) )
               iota2minel(k,ll) = min( iota2taylor(k,ll,1),iota2taylor(EL_NBORS(ell,k),ll,1),
     &              iota2minel(k,ll) )
#endif
               
            enddo
            
         enddo

      enddo



C.....Must generate linear recostructions at vertices



      ZEtaylorvert = 0.D0
      QXtaylorvert = 0.D0
      Qytaylorvert = 0.D0

#ifdef TRACE
      iotataylorvert = 0.D0
#endif

#ifdef CHEM
      iotataylorvert = 0.D0
      iota2taylorvert = 0.D0
#endif

      do k=1,ne

         if (dofs(k).gt.1) then

            do ll=1,dofs(k)

               do lll=1,3

                  if (ll.eq.1) then

                     ZEtaylorvert(k,ll,lll) = ZEtaylor(k,ll,1) + ZEtaylor(k,ll+1,1)*( xi2vert(k,lll) -xi2BCb(k) )
     &                    + ZEtaylor(k,ll+2,1)*( xi1vert(k,lll) - xi1BCb(k) )
                     QXtaylorvert(k,ll,lll) = QXtaylor(k,ll,1) + QXtaylor(k,ll+1,1)*( xi2vert(k,lll) -xi2BCb(k) )
     &                    + QXtaylor(k,ll+2,1)*( xi1vert(k,lll) - xi1BCb(k) )
                     Qytaylorvert(k,ll,lll) = QYtaylor(k,ll,1) + QYtaylor(k,ll+1,1)*( xi2vert(k,lll) -xi2BCb(k) ) 
     &                    + QYtaylor(k,ll+2,1)*( xi1vert(k,lll) - xi1BCb(k) )

#ifdef TRACE
                     iotataylorvert(k,ll,lll) = iotataylor(k,ll,1) + 
     &                    iotataylor(k,ll+1,1)*( xi2vert(k,lll) - xi2BCb(k) )
     &                    + iotataylor(k,ll+2,1)*( xi1vert(k,lll) - xi1BCb(k) )
#endif

#ifdef CHEM
                     iotataylorvert(k,ll,lll) = iotataylor(k,ll,1) + 
     &                    iotataylor(k,ll+1,1)*( xi2vert(k,lll) - xi2BCb(k) )
     &                    + iotataylor(k,ll+2,1)*( xi1vert(k,lll) - xi1BCb(k) )
                     iota2taylorvert(k,ll,lll) = iota2taylor(k,ll,1) + 
     &                    iota2taylor(k,ll+1,1)*( xi2vert(k,lll) - xi2BCb(k) )
     &                    + iota2taylor(k,ll+2,1)*( xi1vert(k,lll) - xi1BCb(k) )
#endif

                  elseif (ll.gt.1) then

                     ZEtaylorvert(k,ll,lll) = ZEtaylor(k,ll,1) + 
     &                    ZEtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) ),1)*( xi2vert(k,lll) - xi2BCb(k) ) 
     &                    + ZEtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) )+1,1)*( xi1vert(k,lll) - xi1BCb(k) )
                     QXtaylorvert(k,ll,lll) = QXtaylor(k,ll,1) + 
     &                    QXtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) ),1)*( xi2vert(k,lll) - xi2BCb(k) ) 
     &                    + QXtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) )+1,1)*( xi1vert(k,lll) - xi1BCb(k) )
                     Qytaylorvert(k,ll,lll) = QYtaylor(k,ll,1) + 
     &                    QYtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) ),1)*( xi2vert(k,lll) - xi2BCb(k) ) 
     &                    + QYtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) )+1,1)*( xi1vert(k,lll) - xi1BCb(k) )

#ifdef TRACE
                     iotataylorvert(k,ll,lll) = iotataylor(k,ll,1)+
     &                    iotataylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) ),1)*( xi2vert(k,lll) - xi2BCb(k))
     &                    + iotataylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) )+1,1)*( xi1vert(k,lll) - xi1BCb(k) )
#endif

#ifdef CHEM
                     iotataylorvert(k,ll,lll) = iotataylor(k,ll,1)+
     &                    iotataylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) ),1)*( xi2vert(k,lll) - xi2BCb(k))
     &                    + iotataylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) )+1,1)*( xi1vert(k,lll) - xi1BCb(k) )
                     iota2taylorvert(k,ll,lll) = iota2taylor(k,ll,1)+
     &                    iotataylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) ),1)*( xi2vert(k,lll) - xi2BCb(k))
     &                    + iota2taylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) )+1,1)*( xi1vert(k,lll) - xi1BCb(k) )
#endif

                  else

                  endif

               enddo
               
            enddo
            
         elseif (dofs(k).eq.1) then

            cycle

         endif

      enddo

C.....Compute alphas for each variable in each order derivitive


      alphaZE0 = 0.D0
      alphaQX0 = 0.D0
      alphaQY0 = 0.D0

#ifdef TRACE
      alphaiota0 = 0.D0
#endif

#ifdef CHEM
      alphaiota0 = 0.D0
      alphaiota20 = 0.D0
#endif

      do k = 1,ne

         if (dofs(k).gt.1) then

            do ll = 1,dofs(k)

               do lll = 1,3

                  if ( ZEtaylorvert(k,ll,lll).gt.ZEtaylor(k,ll,1).and.
     &                 abs((ZEtaylorvert(k,ll,lll)-ZEtaylor(k,ll,1))).gt.1.0E-15 ) then 

                     alphaZE0(k,ll,lll) = min( 1.D0, ( ZEmaxel(k,ll) - 
     &                    ZEtaylor(k,ll,1) ) / (  ZEtaylorvert(k,ll,lll) - ZEtaylor(k,ll,1) ) )
                     
                  elseif ( ZEtaylorvert(k,ll,lll).eq.ZEtaylor(k,ll,1).or.
     &                    abs((ZEtaylorvert(k,ll,lll)-ZEtaylor(k,ll,1))).le.1.0E-15 ) then

                     alphaZE0(k,ll,lll) = 1.D0

                  elseif (  ZEtaylorvert(k,ll,lll).lt.ZEtaylor(k,ll,1).and.
     &                    abs((ZEtaylorvert(k,ll,lll)-ZEtaylor(k,ll,1))).gt.1.0E-15 ) then

                     alphaZE0(k,ll,lll) = min( 1.D0, ( ZEminel(k,ll)
     &                    - ZEtaylor(k,ll,1) ) / ( ZEtaylorvert(k,ll,lll) - ZEtaylor(k,ll,1) ) )
                     
                  endif

#ifdef TRACE
                  if (iotaminel(k,ll).ne.iotamaxel(k,ll)) then
                     
                     if ( iotataylorvert(k,ll,lll).gt.iotataylor(k,ll,1).and.
     &                    abs((iotataylorvert(k,ll,lll)-iotataylor(k,ll,1))).gt.1.0E-15 ) then 

                        alphaiota0(k,ll,lll) = min(1.D0,( iotamaxel(k,ll)
     &                       -iotataylor(k,ll,1) )/ (iotataylorvert(k,ll,lll) - iotataylor(k,ll,1)))

                                !adapted part

                     elseif ( ( iotataylorvert(k,ll,lll).gt.iotataylor(k,ll,1) ).and.
     &                       ( abs(iotataylorvert(k,ll,lll)-iotataylor(k,ll,1)).gt.1.0E-15 ).and.
     &                       ( iotamaxel(k,ll).eq.iotataylor(k,ll,1) ) ) then 

                        alphaiota0(k,ll,lll) = min(fd, abs( ( iotamaxel(k,ll)
     &                       - iotaminel(k,ll) )/(iotataylorvert(k,ll,lll) - iotamaxel(k,ll)) ) ) 

                                !adapted part

                     elseif ( ( iotataylorvert(k,ll,lll).gt.iotataylor(k,ll,1) ).and.
     &                       ( abs(iotataylorvert(k,ll,lll)-iotataylor(k,ll,1)).gt.1.0E-15 ).and.
     &                       ( iotaminel(k,ll).eq.iotataylor(k,ll,1) ) ) then 

                        alphaiota0(k,ll,lll) = min(fd, abs( ( iotaminel(k,ll)
     &                       - iotamaxel(k,ll) )/(iotataylorvert(k,ll,lll) - iotaminel(k,ll)) ) )
                        
                     elseif (iotataylorvert(k,ll,lll).eq.iotataylor(k,ll,1).or.
     &                       abs((iotataylorvert(k,ll,lll)-iotataylor(k,ll,1))).le.1.0E-15 ) then

                        alphaiota0(k,ll,lll) = 1.D0

                     elseif (iotataylorvert(k,ll,lll).lt.iotataylor(k,ll,1).
     &                       and.abs((iotataylorvert(k,ll,lll)-iotataylor(k,ll,1))).gt.1.0E-15) then

                        alphaiota0(k,ll,lll) = min( 1.D0,( iotaminel(k,ll)
     &                       -iotataylor(k,ll,1) )/( iotataylorvert(k,ll,lll)-iotataylor(k,ll,1)))
                        
                     endif

                  else
                     
                     alphaiota0(k,ll,lll) = 1.D0

                  endif
#endif

#ifdef CHEM                  
                  if ( iotataylorvert(k,ll,lll).gt.iotataylor(k,ll,1).and.
     &                 abs((iotataylorvert(k,ll,lll)-iotataylor(k,ll,1))).gt.1.0E-15 ) then 

                     alphaiota0(k,ll,lll) = min(1.D0,( iotamaxel(k,ll)-
     &                    iotataylor(k,ll,1) )/ (iotataylorvert(k,ll,lll) - iotataylor(k,ll,1)))
                     
                  elseif (iotataylorvert(k,ll,lll).eq.iotataylor(k,ll,1).or.
     &                    abs((iotataylorvert(k,ll,lll)-iotataylor(k,ll,1))).le.1.0E-15 ) then

                     alphaiota0(k,ll,lll) = 1.D0

                  elseif (iotataylorvert(k,ll,lll).lt.iotataylor(k,ll,1).and.
     &                    abs((iotataylorvert(k,ll,lll)-iotataylor(k,ll,1))).gt.1.0E-15) then

                     alphaiota0(k,ll,lll) = min( 1.D0,( iotaminel(k,ll)
     &                    -iotataylor(k,ll,1) )/( iotataylorvert(k,ll,lll)-iotataylor(k,ll,1)))
                     
                  endif

                  if ( iota2taylorvert(k,ll,lll).gt.iota2taylor(k,ll,1).and.
     &                 abs((iota2taylorvert(k,ll,lll)-iota2taylor(k,ll,1))).gt.1.0E-15 ) then 

                     alphaiota20(k,ll,lll) = min(1.D0,( iota2maxel(k,ll)
     &                    -iota2taylor(k,ll,1) )/ (iota2taylorvert(k,ll,lll) - iota2taylor(k,ll,1)))
                     
                  elseif (iota2taylorvert(k,ll,lll).eq.iota2taylor(k,ll,1).or.
     &                    abs((iota2taylorvert(k,ll,lll)-iota2taylor(k,ll,1))).le.1.0E-15 ) then

                     alphaiota20(k,ll,lll) = 1.D0

                  elseif (iota2taylorvert(k,ll,lll).lt.iota2taylor(k,ll,1).and.
     &                    abs((iota2taylorvert(k,ll,lll)-iota2taylor(k,ll,1))).gt.1.0E-15) then

                     alphaiota20(k,ll,lll) = min( 1.D0,( iota2minel(k,ll)
     &                    -iota2taylor(k,ll,1) )/( iota2taylorvert(k,ll,lll)-iota2taylor(k,ll,1)))
                     
                  endif
#endif

                  if ( QXtaylorvert(k,ll,lll).gt.QXtaylor(k,ll,1).and.
     &                 (QXtaylorvert(k,ll,lll)-QXtaylor(k,ll,1)).gt.1.0E-15 ) then !in xi1

                     alphaQX0(k,ll,lll) = min( 1.D0, ( QXmaxel(k,ll) 
     &                    - QXtaylor(k,ll,1) ) / ( QXtaylorvert(k,ll,lll) - QXtaylor(k,ll,1) ) )
                     
                  elseif ( QXtaylorvert(k,ll,lll).eq.QXtaylor(k,ll,1).or.
     &                    (QXtaylorvert(k,ll,lll)-QXtaylor(k,ll,1)).le.1.0E-15  ) then

                     alphaQX0(k,ll,lll) = 1.D0

                  elseif ( QXtaylorvert(k,ll,lll).lt.QXtaylor(k,ll,1).and.
     &                    (QXtaylorvert(k,ll,lll)-QXtaylor(k,ll,1)).gt.1.0E-15 ) then

                     alphaQX0(k,ll,lll) = min( 1.D0, ( QXminel(k,ll)
     &                    - QXtaylor(k,ll,1) ) / ( QXtaylorvert(k,ll,lll) - QXtaylor(k,ll,1) ) )

                  endif

                  if ( QYtaylorvert(k,ll,lll).gt.QYtaylor(k,ll,1).and.
     &                 (QYtaylorvert(k,ll,lll)-QYtaylor(k,ll,1)).gt.1.0E-15 ) then !in xi1

                     alphaQY0(k,ll,lll) = min( 1.D0, ( QYmaxel(k,ll) 
     &                    - QYtaylor(k,ll,1) ) / ( QYtaylorvert(k,ll,lll) - QYtaylor(k,ll,1) ) )
                     
                  elseif ( QYtaylorvert(k,ll,lll).eq.QYtaylor(k,ll,1).or.
     &                    (QYtaylorvert(k,ll,lll)-QYtaylor(k,ll,1)).le.1.0E-15  ) then

                     alphaQY0(k,ll,lll) = 1.D0

                  elseif ( QYtaylorvert(k,ll,lll).lt.QYtaylor(k,ll,1).and.
     &                    (QYtaylorvert(k,ll,lll)-QYtaylor(k,ll,1)).gt.1.0E-15 ) then

                     alphaQY0(k,ll,lll) = min( 1.D0, ( QYminel(k,ll) 
     &                    - QYtaylor(k,ll,1) ) / ( QYtaylorvert(k,ll,lll)  - QYtaylor(k,ll,1) ) )

                  endif

               enddo            !lll

            enddo               !ll

         elseif (dofs(k).eq.1) then

            cycle

         endif

      enddo

C.....Find the prescribed higher limiters by finding smallest local value

      alphaZE = 0.D0
      alphaQX = 0.D0
      alphaQY = 0.D0

#ifdef TRACE
      alphaiota = 0.D0
#endif

#ifdef CHEM
      alphaiota = 0.D0
      alphaiota2 = 0.D0
#endif

      do k=1,ne

         if (dofs(k).gt.1) then

            do ll=1,dofs(k)		

               alphaZE(k,ll) = minval( alphaZE0(k,ll,:) )
               alphaQX(k,ll) = minval( alphaQX0(k,ll,:) )
               alphaQY(k,ll) = minval( alphaQY0(k,ll,:) )

#ifdef TRACE
               alphaiota(k,ll) = minval( alphaiota0(k,ll,:) )
#endif

#ifdef CHEM
               alphaiota(k,ll) = minval( alphaiota0(k,ll,:) )
               alphaiota2(k,ll) = minval( alphaiota20(k,ll,:) )
#endif

            enddo

         elseif (dofs(k).eq.1) then

            cycle

         endif

      enddo

C.... Choose smallest (minimum) alpha for derivative in x or y

      alphaZEm = 0.D0
      alphaQXm = 0.D0
      alphaQYm = 0.D0

#ifdef TRACE
      alphaiotam = 0.D0
#endif

#ifdef CHEM
      alphaiotam = 0.D0
      alphaiota2m = 0.D0
#endif

      do k=1,ne

         if (dofs(k).gt.1) then
            
            do bb = 1,pdg_el(k)

               if( (bb+1)*(bb+2)/2.le.dofs(k) ) then

                  alphaZEm(k,(bb*(bb+1))/2 + 1) = 
     &                 minval( alphaZE(k,(bb*(bb-1))/2 +1:(bb*(bb+1)) / 2  ) )
                  alphaQXm(k,(bb*(bb+1))/2 + 1) = 
     &                 minval( alphaQX(k,(bb*(bb-1))/2 +1:(bb*(bb+1)) / 2  ) )
                  alphaQYm(k,(bb*(bb+1))/2 + 1) = 
     &                 minval( alphaQY(k,(bb*(bb-1))/2 +1:(bb*(bb+1)) / 2  ) )

#ifdef TRACE
                  alphaiotam(k,(bb*(bb+1))/2 + 1) = 
     &                 minval( alphaiota(k,(bb*(bb-1))/2 +1:(bb*(bb+1)) / 2 ) )
#endif

#ifdef CHEM
                  alphaiotam(k,(bb*(bb+1))/2 + 1) = 
     &                 minval( alphaiota(k,(bb*(bb-1))/2 +1:(bb*(bb+1)) / 2  ) )
                  alphaiota2m(k,(bb*(bb+1))/2 + 1) = 
     &                 minval( alphaiota2(k,(bb*(bb-1))/2 +1:(bb*(bb+1)) / 2 ) )
#endif

               endif
               
            enddo

         elseif (dofs(k).eq.1) then

            cycle

         endif

      enddo

C.....Use max higher derivative values for final limiter value

      alphaZE_max = 0.D0
      alphaQX_max = 0.D0
      alphaQY_max = 0.D0

#ifdef TRACE
      alphaiota_max = 0.D0
#endif

#ifdef CHEM
      alphaiota_max = 0.D0
      alphaiota2_max = 0.D0
#endif

      do k=1,ne

         if (dofs(k).gt.1) then

            do bb =1,pdg_el(k)

               if( (bb+1)*(bb+2)/2.le.dofs(k)) then

                  alphaZE_max(k,(bb*(bb+1))/2 + 1) = fd*maxval( 
     &                 alphaZEm(k,(bb*(bb+1))/2 + 1:dofs(k)) )
                  alphaQX_max(k,(bb*(bb+1))/2 + 1) = fd*maxval( 
     &                 alphaQXm(k,(bb*(bb+1))/2 + 1:dofs(k)) )
                  alphaQY_max(k,(bb*(bb+1))/2 + 1) = fd*maxval( 
     &                 alphaQYm(k,(bb*(bb+1))/2 + 1:dofs(k)) )

#ifdef TRACE
                  alphaiota_max(k,(bb*(bb+1))/2 + 1) = fd*maxval( 
     &                 alphaiotam(k,(bb*(bb+1))/2 + 1:dofs(k)) )
#endif

#ifdef CHEM
                  alphaiota_max(k,(bb*(bb+1))/2 + 1) = fd*maxval( 
     &                 alphaiotam(k,(bb*(bb+1))/2 + 1:dofs(k)) )
                  alphaiota2_max(k,(bb*(bb+1))/2 + 1) = fd*maxval( 
     &                 alphaiota2m(k,(bb*(bb+1))/2 + 1:dofs(k)) )
#endif

               endif

            enddo

         elseif (dofs(k).eq.1) then

            cycle

         endif

      enddo

C.....Limit on the Master element in the Taylor basis, via reconstruction 
C.....of unconstrained solutions with alpha constraints


      limitZE = 0.D0
      limitQX = 0.D0
      limitQY = 0.D0

      lim_count_roll = 0

#ifdef TRACE
      limitiota = 0.D0
#endif

#ifdef CHEM
      limitiota = 0.D0
      limitiota2 = 0.D0
#endif

      do k=1,ne

         lim_count = 0

         if (dofs(k).gt.1) then

            do ll=1,dofs(k)

               if ( ll.eq.1 ) then

                  limitZE(k,ll) = ZEtaylor(k,ll,1)
                  limitQX(k,ll) = QXtaylor(k,ll,1)
                  limitQY(k,ll) = QYtaylor(k,ll,1) 

#ifdef TRACE
                  limitiota(k,ll) = iotataylor(k,ll,1)
#endif

#ifdef CHEM
                  limitiota(k,ll) = iotataylor(k,ll,1)
                  limitiota2(k,ll) = iota2taylor(k,ll,1)
#endif

               elseif ( ll.ge.2 ) then
                  
                  do bb=1,pdg_el(k)

                     if ( ll.le.( ( (bb+1)*(bb+2)) / 2.D0 ).and.(ll.gt.
     &                    (bb*(bb+1)/2.D0) ) ) then

                        limitZE(k,ll) = alphaZE_max(k,(bb*(bb+1))/2 + 1) 
     &                       * ZEtaylor(k,ll,1)
                        limitQX(k,ll) = alphaQX_max(k,(bb*(bb+1))/2 + 1) 
     &                       * ZEtaylor(k,ll,1)
                        limitQY(k,ll) = alphaQY_max(k,(bb*(bb+1))/2 + 1) 
     &                       * ZEtaylor(k,ll,1)

#ifdef TRACE
                        limitiota(k,ll) = alphaiota_max(k,(bb*(bb+1))/2 + 1) 
     &                       * iotataylor(k,ll,1)
#endif

#ifdef CHEM
                        limitiota(k,ll) = alphaiota_max(k,(bb*(bb+1))/2 + 1) 
     &                       * iotataylor(k,ll,1)
                        limitiota2(k,ll) = alphaiota2_max(k,(bb*(bb+1))/2 + 1) 
     &                       * iota2taylor(k,ll,1)
#endif


c$$$  ! Make a counter to track limiting
c$$$  
c$$$  if ( ( alphaZE_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
c$$$  &                       alphaQX_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
c$$$  &                       alphaQY_max(k,(bb*(bb+1))/2 + 1).lt.1.D0 ).and.
c$$$  &                       chem_flag.ne.1.and.tracer_flag.ne.1 ) then
c$$$  
c$$$  lim_count = 1  
c$$$  
c$$$  elseif ( (alphaZE_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
c$$$  &                          alphaQX_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
c$$$  &                          alphaQY_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
c$$$  &                          alphaiota_max(k,(bb*(bb+1))/2 + 1).lt.1.D0 ).and.
c$$$  &                          chem_flag.ne.1.and.tracer_flag.eq.1 ) then
c$$$  
c$$$  lim_count = 1   
c$$$  
c$$$  elseif ( (alphaZE_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
c$$$  &                          alphaQX_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
c$$$  &                          alphaQY_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
c$$$  &                          alphaiota_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
c$$$  &                          alphaiota2_max(k,(bb*(bb+1))/2 + 1).lt.1.D0 ).and.
c$$$  &                          chem_flag.eq.1 ) then
c$$$  
c$$$  lim_count = 1 
c$$$  
c$$$  endif
                        
                     endif
                     
                  enddo
                  
               endif

            enddo

         elseif (dofs(k).eq.1) then

            cycle

         endif

                                !lim_count_roll = lim_count_roll + lim_count

      enddo

C.....Transform back to the Dubiner basis (multiply by NmatrixInv),

      ZEconst =  0.D0
      QXconst =  0.D0
      QYconst =  0.D0

#ifdef TRACE
      iotaconst =  0.D0
#endif

#ifdef CHEM
      iotaconst =  0.D0
      iota2const =  0.D0
#endif

      do k=1,ne

         if (dofs(k).gt.1) then
                                !do lll=1,3

            do ll=1,dofs(k)

               do ss=1,dofs(k)

                  ZEconst(k,ll) = ZEconst(k,ll) + NmatrixInv(k,ll,ss,dofs(k)) 
     &                 * limitZE(k,ss)
                  QXconst(k,ll) = QXconst(k,ll) + NmatrixInv(k,ll,ss,dofs(k)) 
     &                 * limitQX(k,ss)
                  QYconst(k,ll) = QYconst(k,ll) + NmatrixInv(k,ll,ss,dofs(k)) 
     &                 * limitQY(k,ss)

#ifdef TRACE
                  iotaconst(k,ll) = iotaconst(k,ll) + 
     &                 NmatrixInv(k,ll,ss,dofs(k)) * limitiota(k,ss)
#endif

#ifdef CHEM
                  iotaconst(k,ll) = iotaconst(k,ll) + 
     &                 NmatrixInv(k,ll,ss,dofs(k)) * limitiota(k,ss)
                  iota2const(k,ll) = iota2const(k,ll) + 
     &                 NmatrixInv(k,ll,ss,dofs(k)) * limitiota2(k,ss)
#endif

               enddo

            enddo

         elseif (dofs(k).eq.1) then

            cycle

         endif

      enddo

C.....Set limit values

      do k=1,ne

         if (dofs(k).gt.1) then

            do ll = 1,dofs(k)

               ZE(ll,k,irk+1) = ZEconst(k,ll)
               QX(ll,k,irk+1) = QXconst(k,ll)
               QY(ll,k,irk+1) = QYconst(k,ll)

#ifdef TRACE
               iota(ll,k,irk+1) = iotaconst(k,ll)
#endif

#ifdef CHEM
               iota(ll,k,irk+1) = iotaconst(k,ll)
               iota2(ll,k,irk+1) = iota2const(k,ll)
#endif

            enddo

         elseif (dofs(k).eq.1) then

            cycle

         endif

      enddo
      
      return
      end subroutine


C***********************************************************************
C     
C     SUBROUTINE SLOPELIMITER7()
C     
C     Written 2011
C     
C     This subroutine selects the hierarchical reconstruction method over
C     a Taylor Polynomial basis, and is consistent with p-adaptation
C     to arbitrary order p
C     
C     --cem
C     
C***********************************************************************

      SUBROUTINE SLOPELIMITER7()

C.....Use appropriate modules

      USE GLOBAL
      USE DG

#ifdef CMPI
      USE MESSENGER_ELEM
      USE MESSENGER
#endif

      IMPLICIT NONE

      Integer k,ll,ss,lll,ell,bb,MUSCL,ENO,Mixed_MUSCL,Mixed_ENO,b_index,i,j,mm,neigh_flag

      REAL(SZ) fd,marea,fde,new_iota,old_iota
      Real(SZ), Allocatable :: ZE_cand(:,:,:),QX_cand(:,:,:),QY_cand(:,:,:)
      Real(SZ), Allocatable :: iota_cand(:,:,:),iota2_cand(:,:,:)
      Real(SZ), Allocatable :: ZE_minmod(:,:,:),QX_minmod(:,:,:),QY_minmod(:,:,:)
      Real(SZ), Allocatable :: iota_minmod(:,:,:),iota2_minmod(:,:,:)
      Real(SZ), Allocatable :: ZEtaylor0slope(:,:,:),ZEtaylor0(:,:,:)
      Real(SZ), Allocatable :: QXtaylor0slope(:,:,:),QXtaylor0(:,:,:)
      Real(SZ), Allocatable :: QYtaylor0slope(:,:,:),QYtaylor0(:,:,:)
      Real(SZ), Allocatable :: iotataylor0slope(:,:,:),iotataylor0(:,:,:)
      Real(SZ), Allocatable :: iota2taylor0slope(:,:,:),iota2taylor0(:,:,:)


      Allocate ( ZE_cand(mne,dofh,1),QX_cand(mne,dofh,1) )
      Allocate ( QY_cand(mne,dofh,1) ) 
      Allocate ( iota_cand(mne,dofh,1),iota2_cand(mne,dofh,1) )
      Allocate ( ZE_minmod(mne,dofh,1),QX_minmod(mne,dofh,1) )
      Allocate ( QY_minmod(mne,dofh,1) ) 
      Allocate ( iota_minmod(mne,dofh,1),iota2_minmod(mne,dofh,1) )
      Allocate ( ZEtaylor0slope(mne,dofh,1),ZEtaylor0(MNE,dofh,1))
      Allocate ( QXtaylor0slope(mne,dofh,1),QXtaylor0(MNE,dofh,1))
      Allocate ( QYtaylor0slope(mne,dofh,1),QYtaylor0(MNE,dofh,1))
      Allocate ( iotataylor0slope(mne,dofh,1),iotataylor0(MNE,dofh,1))
      Allocate ( iota2taylor0slope(mne,dofh,1),iota2taylor0(MNE,dofh,1))





C.....We work over the master element
C.....Set initial values

      MUSCL = 0                 ! add weight for lower order pieces (fd<1 => stronger limiting)
      ENO = 1
      Mixed_MUSCL = 0
      Mixed_ENO = 0
      neigh_flag = 1
      
      DO k=1,NE

         if (dofs(k).gt.1) then

            DO ll = 1,dofs(k)

               ZE(ll,k,NRK+2) = ZE(ll,k,IRK+1)
               QX(ll,k,NRK+2) = QX(ll,k,IRK+1)
               QY(ll,k,NRK+2) = QY(ll,k,IRK+1)

#ifdef TRACE
               iota(ll,k,NRK+2) = iota(ll,k,IRK+1)
#endif

#ifdef CHEM
               iota(ll,k,NRK+2) = iota(ll,k,IRK+1)
               iota2(ll,k,NRK+2) = iota2(ll,k,IRK+1)
#endif

            ENDDO

         elseif (dofs(k).eq.1) then

            cycle

         endif

      ENDDO   

      marea = 2.D0              !master elements area

C.....Convert initial values to the Taylor basis (multiply by Nmatrix) on base element


      ZEtaylor = 0.D0 
      QXtaylor = 0.D0
      QYtaylor = 0.D0

#ifdef TRACE
      iotataylor = 0.D0      
#endif

#ifdef CHEM
      iotataylor = 0.D0
      iota2taylor = 0.D0
#endif

      do k=1,ne

         if (dofs(k).gt.1) then

            do ll=1,dofs(k)
               
               do ss=1,dofs(k)

                  ZEtaylor(k,ll,1) = ZEtaylor(k,ll,1) + Nmatrix(k,ll,ss,dofs(k)) * ZE(ss,k,nrk+2)
                  QXtaylor(k,ll,1) = QXtaylor(k,ll,1) + Nmatrix(k,ll,ss,dofs(k)) * QX(ss,k,nrk+2)
                  QYtaylor(k,ll,1) = QYtaylor(k,ll,1) + Nmatrix(k,ll,ss,dofs(k)) * QY(ss,k,nrk+2)

#ifdef TRACE
                  iotataylor(k,ll,1) = iotataylor(k,ll,1) + Nmatrix(k,ll,ss,dofs(k))*iota(ss,k,nrk+2)
#endif

#ifdef CHEM
                  iotataylor(k,ll,1) = iotataylor(k,ll,1) + Nmatrix(k,ll,ss,dofs(k))*iota(ss,k,nrk+2)
                  iota2taylor(k,ll,1) = iota2taylor(k,ll,1) + Nmatrix(k,ll,ss,dofs(k))*iota2(ss,k,nrk+2)
#endif
                  
               enddo

            enddo

         elseif (dofs(k).eq.1) then

            cycle

         endif

      enddo


C.....Generate linear recostructions at vertices

      ZEtaylor0 = 0.D0
      QXtaylor0 = 0.D0
      Qytaylor0 = 0.D0

#ifdef TRACE
      iotataylor0 = 0.D0
      iotataylor0slope = 0.D0
#endif

#ifdef CHEM
      iotataylor0 = 0.D0
      iotataylor0slope = 0.D0
      iota2taylor0 = 0.D0
      iota2taylor0slope = 0.D0
#endif

#ifdef CMPI

      CALL UPDATER_ELEM_MOD2(ZEtaylor,QXtaylor,QYtaylor,1,3)


#ifdef TRACE
      CALL UPDATER_ELEM_MOD2(iotataylor,QXtaylor,QYtaylor,1,2 )  
#endif

#ifdef CHEM
      CALL UPDATER_ELEM_MOD2(iotataylor,iota2taylor,QYtaylor,1,2 )  
#endif

#endif


      do k=1,ne

         do ll = dofs(k),3,-1

            ZEtaylor0(k,ll,1) = ZEtaylor(k,ll,1)
            ZEtaylor0slope(k,ll,1) = 0.D0
            QXtaylor0(k,ll,1) = QXtaylor(k,ll,1)
            QXtaylor0slope(k,ll,1) = 0.D0
            QYtaylor0(k,ll,1) = QYtaylor(k,ll,1)
            QYtaylor0slope(k,ll,1) = 0.D0
            
#ifdef TRACE
            iotataylor0(k,ll,1) = iotataylor(k,ll,1)
            iotataylor0slope(k,ll,1) = 0.D0
#endif

#ifdef CHEM            
            iotataylor0(k,ll,1) = iotataylor(k,ll,1)
            iotataylor0slope(k,ll,1) = 0.D0
            iota2taylor0(k,ll,1) = iota2taylor(k,ll,1)
            iota2taylor0slope(k,ll,1) = 0.D0
#endif

            if (ll.gt.3) then

                                ! Construct the coeffs in nonlinear reconstruction

               do i = 0,pdg_el(k)

                  do j = 0,pdg_el(k)

                     b_index = (bj(ll)+j) * (bj(ll) + 1 + j) / 2 + ( bj(ll)+j ) * ( bi(ll)+i ) 
     &                    + ( bi(ll) + i )*(bi(ll) + 3 + i)/2 + 1

                     if ( b_index.le.dofs(k).and.(i+j).gt.0) then


                        Call factorial(i,fact(i))
                        Call factorial(j,fact(j))
                        
                        do mm = 1,nagp(pdg_el(k))

                           ZEtaylor0(k,ll,1) =   ZEtaylor0(k,ll,1) + ZEtaylor(k,b_index,1) * 
     &                          (1.D0/fact(i)*fact(j))* ( ( yagp(mm,pdg_el(k)) - xi2BCb(k) )**i
     &                          *( xagp(mm,pdg_el(k)) - xi1BCb(k) )**j  )  * wagp(mm,pdg_el(k))

                           ZEtaylor0slope(k,ll,1) =  ZEtaylor0slope(k,ll,1) + ZEtaylor(k,b_index,1) * 
     &                          (1.D0/fact(i)*fact(j))* ( ( yagp(mm,pdg_el(k)) - xi2BCb(k) )**i
     &                          *( xagp(mm,pdg_el(k)) - xi1BCb(k) )**j  )  * wagp(mm,pdg_el(k))

                           QXtaylor0(k,ll,1) =   QXtaylor0(k,ll,1) + QXtaylor(k,b_index,1) * 
     &                          (1.D0/fact(i)*fact(j))* ( ( yagp(mm,pdg_el(k)) - xi2BCb(k) )**i
     &                          *( xagp(mm,pdg_el(k)) - xi1BCb(k) )**j  )  * wagp(mm,pdg_el(k))

                           QXtaylor0slope(k,ll,1) =  QXtaylor0slope(k,ll,1) + QXtaylor(k,b_index,1) * 
     &                          (1.D0/fact(i)*fact(j))* ( ( yagp(mm,pdg_el(k)) - xi2BCb(k) )**i
     &                          *( xagp(mm,pdg_el(k)) - xi1BCb(k) )**j  )  * wagp(mm,pdg_el(k))

                           QYtaylor0(k,ll,1) =   QYtaylor0(k,ll,1) + QYtaylor(k,b_index,1) * 
     &                          (1.D0/fact(i)*fact(j))* ( ( yagp(mm,pdg_el(k)) - xi2BCb(k) )**i
     &                          *( xagp(mm,pdg_el(k)) - xi1BCb(k) )**j  )  * wagp(mm,pdg_el(k))

                           QYtaylor0slope(k,ll,1) =  QYtaylor0slope(k,ll,1) + QYtaylor(k,b_index,1) * 
     &                          (1.D0/fact(i)*fact(j))* ( ( yagp(mm,pdg_el(k)) - xi2BCb(k) )**i
     &                          *( xagp(mm,pdg_el(k)) - xi1BCb(k) )**j  )  * wagp(mm,pdg_el(k))

#ifdef TRACE
                           iotataylor0(k,ll,1) =   iotataylor0(k,ll,1) + iotataylor(k,b_index,1) * 
     &                          (1.D0/fact(i)*fact(j))* ( ( yagp(mm,pdg_el(k)) - xi2BCb(k) )**i
     &                          *( xagp(mm,pdg_el(k)) - xi1BCb(k) )**j  )  * wagp(mm,pdg_el(k))

                           iotataylor0slope(k,ll,1) =  iotataylor0slope(k,ll,1) + iotataylor(k,b_index,1) * 
     &                          (1.D0/fact(i)*fact(j))* ( ( yagp(mm,pdg_el(k)) - xi2BCb(k) )**i
     &                          *( xagp(mm,pdg_el(k)) - xi1BCb(k) )**j  )  * wagp(mm,pdg_el(k))
#endif

#ifdef CHEM
                           iotataylor0(k,ll,1) =   iotataylor0(k,ll,1) + iotataylor(k,b_index,1) * 
     &                          (1.D0/fact(i)*fact(j))* ( ( yagp(mm,pdg_el(k)) - xi2BCb(k) )**i
     &                          *( xagp(mm,pdg_el(k)) - xi1BCb(k) )**j  )  * wagp(mm,pdg_el(k))

                           iotataylor0slope(k,ll,1) =  iotataylor0slope(k,ll,1) + iotataylor(k,b_index,1) * 
     &                          (1.D0/fact(i)*fact(j))* ( ( yagp(mm,pdg_el(k)) - xi2BCb(k) )**i
     &                          *( xagp(mm,pdg_el(k)) - xi1BCb(k) )**j  )  * wagp(mm,pdg_el(k))

                           iota2taylor0(k,ll,1) =   iota2taylor0(k,ll,1) + iota2taylor(k,b_index,1) * 
     &                          (1.D0/fact(i)*fact(j))* ( ( yagp(mm,pdg_el(k)) - xi2BCb(k) )**i
     &                          *( xagp(mm,pdg_el(k)) - xi1BCb(k) )**j  )  * wagp(mm,pdg_el(k))

                           iota2taylor0slope(k,ll,1) =  iota2taylor0slope(k,ll,1) + iota2taylor(k,b_index,1) * 
     &                          (1.D0/fact(i)*fact(j))* ( ( yagp(mm,pdg_el(k)) - xi2BCb(k) )**i
     &                          *( xagp(mm,pdg_el(k)) - xi1BCb(k) )**j  )  * wagp(mm,pdg_el(k))
#endif

                        enddo
                        

                     endif

                  enddo
                  
               enddo

            endif


         enddo


      enddo

                                !Find Candidates

      ZE_cand = 0.D0
      QX_cand = 0.D0
      QY_cand = 0.D0

      iota_cand = 0.D0

      iota2_cand = 0.D0

      do k=1,mne

         if (dofs(k).gt.3) then
            
            do ll=dofs(k),3,-1

               ZE_cand(k,ll,1) =  ZEtaylor0(k,ll,1) - ZEtaylor0slope(k,ll,1)

               QX_cand(k,ll,1) =  QXtaylor0(k,ll,1) - QXtaylor0slope(k,ll,1)

               QY_cand(k,ll,1) =  QYtaylor0(k,ll,1) - QYtaylor0slope(k,ll,1)

#ifdef TRACE
               iota_cand(k,ll,1) =  iotataylor0(k,ll,1) - iotataylor0slope(k,ll,1)
#endif

#ifdef CHEM
               iota_cand(k,ll,1) =  iotataylor0(k,ll,1) - iotataylor0slope(k,ll,1)
               iota2_cand(k,ll,1) =  iota2taylor0(k,ll,1) - iota2taylor0slope(k,ll,1)
#endif

            enddo

         endif
         
      enddo


#ifdef CMPI

      CALL UPDATER_ELEM_MOD2(ZE_cand,QX_cand,QY_cand,1,3)


#ifdef TRACE
      CALL UPDATER_ELEM_MOD2(iota_cand,iota_cand,QY_cand,1,2 )
#endif

#ifdef CHEM
      CALL UPDATER_ELEM_MOD2(iota_cand,iota2_cand,QY_cand,1,2 )
#endif

#endif


                                !Find the minmod functions

      ZE_minmod(k,ll,1) = 0.D0
      QX_minmod(k,ll,1) = 0.D0
      QY_minmod(k,ll,1) = 0.D0
      
#ifdef TRACE
      iota_minmod(k,ll,1) = 0.D0
#endif

#ifdef CHEM
      iota_minmod(k,ll,1) = 0.D0
      iota2_minmod(k,ll,1) = 0.D0
#endif

c$$$  #ifdef CMPI
c$$$  
c$$$  CALL UPDATER_ELEM_MOD2(ZE_minmod,QX_minmod,QY_minmod,1,3)
c$$$  
c$$$  #ifdef TRACE
c$$$  CALL UPDATER_ELEM_MOD2(iota_minmod,iota_minmod,QY_minmod,1,2 )
c$$$  #endif
c$$$  
c$$$  #ifdef CHEM
c$$$  CALL UPDATER_ELEM_MOD2(iota_minmod,iota2_minmod,QY_minmod,1,2 )
c$$$  #endif
c$$$  
c$$$  #endif

                                !Apply the friendly minmod regime
      
      do k=1,mne

         do ll=1,dofs(k)

            ZE_minmod(k,ll,1) = ZEtaylor(k,ll,1)
            QX_minmod(k,ll,1) = QXtaylor(k,ll,1)
            QY_minmod(k,ll,1) = QYtaylor(k,ll,1)

#ifdef TRACE
            iota_minmod(k,ll,1) = iotataylor(k,ll,1)
#endif

#ifdef CHEM
            iota_minmod(k,ll,1) = iotataylor(k,ll,1)
            iota2_minmod(k,ll,1) = iota2taylor(k,ll,1)
#endif
            if (ll.gt.1) then

               if (neigh_flag.eq.1) then !focal neighbors

                  if (MUSCL.eq.1) then !construct MUSCL minmod


                     if (minval(ZE_cand(focal_neigh(k,1:focal_up(k)),ll,1)).gt.0.D0 ) then

                        ZE_minmod(k,ll,1) = minval(ZE_cand(focal_neigh(k,1:focal_up(k)),ll,1))

                     elseif (maxval(ZE_cand(focal_neigh(k,1:focal_up(k)),ll,1)).lt.0.D0 ) then

                        ZE_minmod(k,ll,1) = maxval(ZE_cand(focal_neigh(k,1:focal_up(k)),ll,1))

                     else

                        ZE_minmod(k,ll,1) = 0.D0

                     endif

                     if (minval(QX_cand(focal_neigh(k,1:focal_up(k)),ll,1)).gt.0.D0 ) then

                        QX_minmod(k,ll,1) = minval(QX_cand(focal_neigh(k,1:focal_up(k)),ll,1))

                     elseif (maxval(QX_cand(focal_neigh(k,1:focal_up(k)),ll,1)).lt.0.D0 ) then

                        QX_minmod(k,ll,1) = maxval(QX_cand(focal_neigh(k,1:focal_up(k)),ll,1))

                     else

                        QX_minmod(k,ll,1) = 0.D0

                     endif

                     if (minval(QY_cand(focal_neigh(k,1:focal_up(k)),ll,1)).gt.0.D0 ) then

                        QY_minmod(k,ll,1) = minval(QY_cand(focal_neigh(k,1:focal_up(k)),ll,1))

                     elseif (maxval(QY_cand(focal_neigh(k,1:focal_up(k)),ll,1)).lt.0.D0 ) then

                        QY_minmod(k,ll,1) = maxval(QY_cand(focal_neigh(k,1:focal_up(k)),ll,1))

                     else

                        QY_minmod(k,ll,1) = 0.D0

                     endif
                     
                     
                  elseif (ENO.eq.1) then !construct ENO minmod

                     ZE_minmod(k,ll,1) = minval(abs(ZE_cand(focal_neigh(k,1:focal_up(k)),ll,1)))
                     QX_minmod(k,ll,1) = minval(abs(ZE_cand(focal_neigh(k,1:focal_up(k)),ll,1)))
                     QY_minmod(k,ll,1) = minval(abs(ZE_cand(focal_neigh(k,1:focal_up(k)),ll,1)))

                  elseif (Mixed_MUSCL.eq.1) then !construct Mixed minmod

                  elseif (Mixed_ENO.eq.1) then !construct Mixed minmod

                  endif

#ifdef TRACE
                  if (MUSCL.eq.1) then !construct MUSCL minmod


                     if (minval(iota_cand(focal_neigh(k,1:focal_up(k)),ll,1)).gt.0.D0 ) then

                        iota_minmod(k,ll,1) = minval(iota_cand(focal_neigh(k,1:focal_up(k)),ll,1))

                     elseif (maxval(iota_cand(focal_neigh(k,1:focal_up(k)),ll,1)).lt.0.D0 ) then

                        iota_minmod(k,ll,1) = maxval(iota_cand(focal_neigh(k,1:focal_up(k)),ll,1))

                     else

                        iota_minmod(k,ll,1) = 0.D0

                     endif
                     
                     
                  elseif (ENO.eq.1) then !construct ENO minmod

                     iota_minmod(k,ll,1) = minval(abs(iota_cand(focal_neigh(k,1:focal_up(k)),ll,1)))

                  elseif (Mixed_MUSCL.eq.1) then !construct Mixed minmod

                  elseif (Mixed_ENO.eq.1) then !construct Mixed minmod

                  endif
#endif

#ifdef CHEM
                  if (MUSCL.eq.1) then !construct MUSCL minmod


                     if (minval(iota_cand(focal_neigh(k,1:focal_up(k)),ll,1)).gt.0.D0 ) then

                        iota_minmod(k,ll,1) = minval(iota_cand(focal_neigh(k,1:focal_up(k)),ll,1))

                     elseif (maxval(iota_cand(focal_neigh(k,1:focal_up(k)),ll,1)).lt.0.D0 ) then

                        iota_minmod(k,ll,1) = maxval(iota_cand(focal_neigh(k,1:focal_up(k)),ll,1))

                     else

                        iota_minmod(k,ll,1) = 0.D0

                     endif
                     
                     
                  elseif (ENO.eq.1) then !construct ENO minmod

                     iota_minmod(k,ll,1) = minval(abs(iota_cand(focal_neigh(k,1:focal_up(k)),ll,1)))

                  elseif (Mixed_MUSCL.eq.1) then !construct Mixed minmod

                  elseif (Mixed_ENO.eq.1) then !construct Mixed minmod

                  endif

                  if (MUSCL.eq.1) then !construct MUSCL minmod


                     if (minval(iota2_cand(focal_neigh(k,1:focal_up(k)),ll,1)).gt.0.D0 ) then

                        iota2_minmod(k,ll,1) = minval(iota2_cand(focal_neigh(k,1:focal_up(k)),ll,1))

                     elseif (maxval(iota2_cand(focal_neigh(k,1:focal_up(k)),ll,1)).lt.0.D0 ) then

                        iota2_minmod(k,ll,1) = maxval(iota2_cand(focal_neigh(k,1:focal_up(k)),ll,1))

                     else

                        iota2_minmod(k,ll,1) = 0.D0

                     endif
                     
                     
                  elseif (ENO.eq.1) then !construct ENO minmod

                     iota2_minmod(k,ll,1) = minval(abs(iota2_cand(focal_neigh(k,1:focal_up(k)),ll,1)))

                  elseif (Mixed_MUSCL.eq.1) then !construct Mixed minmod

                  elseif (Mixed_ENO.eq.1) then !construct Mixed minmod

                  endif
#endif

               elseif (neigh_flag.eq.2) then !edge neighbors

                  if (MUSCL.eq.1) then !construct MUSCL minmod

                     if (ZE_cand(k,ll,1).gt.0.D0 ) then

                        ZE_minmod(k,ll,1) = minval(ZE_cand(El_nbors(:,k),ll,1))

                     elseif (ZE_cand(k,ll,1).lt.0.D0 ) then

                        ZE_minmod(k,ll,1) = maxval(ZE_cand(El_nbors(:,k),ll,1))

                     else

                        ZE_minmod(k,ll,1) = 0.D0

                     endif

                     if (QX_cand(k,ll,1).gt.0.D0 ) then

                        QX_minmod(k,ll,1) = minval(QX_cand(El_nbors(:,k),ll,1))

                     elseif (QX_cand(k,ll,1).lt.0.D0 ) then

                        QX_minmod(k,ll,1) = maxval(QX_cand(El_nbors(:,k),ll,1))

                     else

                        QX_minmod(k,ll,1) = 0.D0

                     endif

                     if (QY_cand(k,ll,1).gt.0.D0 ) then

                        QY_minmod(k,ll,1) = minval(QY_cand(El_nbors(:,k),ll,1))

                     elseif (QY_cand(k,ll,1).lt.0.D0 ) then

                        QY_minmod(k,ll,1) = maxval(QY_cand(El_nbors(:,k),ll,1))

                     else

                        QY_minmod(k,ll,1) = 0.D0

                     endif

                  elseif (ENO.eq.1) then !construct ENO minmod

                     ZE_minmod(k,ll,1) = minval(abs(ZE_cand(El_nbors(:,k),ll,1)))
                     QX_minmod(k,ll,1) = minval(abs(QX_cand(El_nbors(:,k),ll,1)))
                     QY_minmod(k,ll,1) = minval(abs(QY_cand(El_nbors(:,k),ll,1)))

                  elseif (Mixed_MUSCL.eq.1) then !construct Mixed minmod

                  elseif (Mixed_ENO.eq.1) then !construct Mixed minmod

                  endif

#ifdef TRACE
                  if (MUSCL.eq.1) then !construct MUSCL minmod

                     if (iota_cand(k,ll,1).gt.0.D0 ) then

                        iota_minmod(k,ll,1) = minval(iota_cand(El_nbors(:,k),ll,1))

                     elseif (iota_cand(k,ll,1).lt.0.D0 ) then

                        iota_minmod(k,ll,1) = maxval(iota_cand(El_nbors(:,k),ll,1))

                     else

                        iota_minmod(k,ll,1) = 0.D0

                     endif

                  elseif (ENO.eq.1) then !construct ENO minmod

                     iota_minmod(k,ll,1) = minval(abs(iota_cand(El_nbors(:,k),ll,1)))

                  elseif (Mixed_MUSCL.eq.1) then !construct Mixed minmod

                  elseif (Mixed_ENO.eq.1) then !construct Mixed minmod

                  endif
#endif

#ifdef CHEM
                  if (MUSCL.eq.1) then !construct MUSCL minmod

                     if (iota_cand(k,ll,1).gt.0.D0 ) then

                        iota_minmod(k,ll,1) = minval(iota_cand(El_nbors(:,k),ll,1))

                     elseif (iota_cand(k,ll,1).lt.0.D0 ) then

                        iota_minmod(k,ll,1) = maxval(iota_cand(El_nbors(:,k),ll,1))

                     else

                        iota_minmod(k,ll,1) = 0.D0

                     endif

                  elseif (ENO.eq.1) then !construct ENO minmod

                     iota_minmod(k,ll,1) = minval(abs(iota_cand(El_nbors(:,k),ll,1)))

                  elseif (Mixed_MUSCL.eq.1) then !construct Mixed minmod

                  elseif (Mixed_ENO.eq.1) then !construct Mixed minmod

                  endif

                  if (MUSCL.eq.1) then !construct MUSCL minmod

                     if (iota2_cand(k,ll,1).gt.0.D0 ) then

                        iota2_minmod(k,ll,1) = minval(iota2_cand(El_nbors(:,k),ll,1))

                     elseif (iota2_cand(k,ll,1).lt.0.D0 ) then

                        iota2_minmod(k,ll,1) = maxval(iota2_cand(El_nbors(:,k),ll,1))

                     else

                        iota2_minmod(k,ll,1) = 0.D0

                     endif

                  elseif (ENO.eq.1) then !construct ENO minmod

                     iota2_minmod(k,ll,1) = minval(abs(iota2_cand(El_nbors(:,k),ll,1)))

                  elseif (Mixed_MUSCL.eq.1) then !construct Mixed minmod

                  elseif (Mixed_ENO.eq.1) then !construct Mixed minmod

                  endif
#endif

               endif
               
            endif

         enddo

      enddo

#ifdef CMPI

      CALL UPDATER_ELEM_MOD2(ZE_minmod,QX_minmod,QY_minmod,1,3)

#ifdef TRACE
      CALL UPDATER_ELEM_MOD2(iota_minmod,iota_minmod,QY_minmod,1,2 )
#endif

#ifdef CHEM
      CALL UPDATER_ELEM_MOD2(iota_minmod,iota2_minmod,QY_minmod,1,2 )
#endif

#endif

C.....Transform back to the Dubiner basis (multiply by NmatrixInv),

      ZEconst =  0.D0
      QXconst =  0.D0
      QYconst =  0.D0

#ifdef TRACE
      iotaconst =  0.D0
#endif

#ifdef CHEM
      iotaconst =  0.D0
      iota2const =  0.D0
#endif

      do k=1,ne

         if (dofs(k).gt.1) then
                                !do lll=1,3

            do ll=1,dofs(k)

               do ss=1,dofs(k)

                  ZEconst(k,ll) = ZEconst(k,ll) + NmatrixInv(k,ll,ss,dofs(k)) 
     &                 * ZE_minmod(k,ss,1)
                  QXconst(k,ll) = QXconst(k,ll) + NmatrixInv(k,ll,ss,dofs(k)) 
     &                 * QX_minmod(k,ss,1)
                  QYconst(k,ll) = QYconst(k,ll) + NmatrixInv(k,ll,ss,dofs(k)) 
     &                 * QX_minmod(k,ss,1)

#ifdef TRACE
                  iotaconst(k,ll) = iotaconst(k,ll) + 
     &                 NmatrixInv(k,ll,ss,dofs(k)) * iota_minmod(k,ss,1) 
#endif

#ifdef CHEM
                  iotaconst(k,ll) = iotaconst(k,ll) + 
     &                 NmatrixInv(k,ll,ss,dofs(k)) *  iota_minmod(k,ss,1)
                  iota2const(k,ll) = iota2const(k,ll) + 
     &                 NmatrixInv(k,ll,ss,dofs(k)) * iota2_minmod(k,ss,1)
#endif

               enddo

            enddo

         elseif (dofs(k).eq.1) then

            cycle

         endif

      enddo

C.....Set limit values
      
      CALL SLOPELIMITER5()

      do k=1,ne

         if (dofs(k).gt.3) then

            do ll = 3,dofs(k)

               ZE(ll,k,irk+1) = ZEconst(k,ll)
               QX(ll,k,irk+1) = QXconst(k,ll)
               QY(ll,k,irk+1) = QYconst(k,ll)

#ifdef TRACE
               iota(ll,k,irk+1) = iotaconst(k,ll)
#endif

#ifdef CHEM
               iota(ll,k,irk+1) = iotaconst(k,ll)
               iota2(ll,k,irk+1) = iota2const(k,ll)
#endif

            enddo

         elseif (dofs(k).eq.1) then

            cycle

         endif

      enddo
      
      return
      end subroutine


C***********************************************************************
C     
C     SUBROUTINE SLOPELIMITER8()
C     
C     Written 2011
C     
C     The hierarchic recombination method, works for adapting p
C     
C     -cem
C     
C     
C***********************************************************************

      SUBROUTINE SLOPELIMITER8()

C.....Use appropriate modules

      use sizes, only : sz
      use global
      use dg

#ifdef CMPI
      USE MESSENGER_ELEM
      USE MESSENGER
#endif

      IMPLICIT NONE

C.....Declare local variables

      INTEGER L, LL, INC1,INC2,INC3,KDP,NN,IVAR,I,J,kk,k,varnum,ss
      REAL(SZ) temp2(MNE)
      Real(SZ), Allocatable :: ZE_MIN1(:),ZE_MAX1(:),QX_MIN1(:),QX_MAX1(:)
      Real(SZ), Allocatable :: QY_MIN1(:),QY_MAX1(:)
      Real(SZ), Allocatable :: iota_MIN1(:),iota_MAX1(:)
      Real(SZ), Allocatable :: iota2_MIN1(:),iota2_MAX1(:),temparray1(:,:,:)
      Real(SZ), Allocatable :: temparray2(:,:,:),temparray3(:,:,:)
      Real(SZ), Allocatable :: temparray4(:,:,:),temparray5(:,:,:)

      Allocate ( ZE_MIN1(NP),ZE_MAX1(NP),QX_MIN1(NP) )
      Allocate ( QY_MIN1(NP),QY_MAX1(NP),QX_MAX1(NP) )
      Allocate ( iota_MIN1(NP),iota_MAX1(NP) )
      Allocate ( iota2_MIN1(NP),iota2_MAX1(NP), temparray1(MNE,dofh,1) )
      Allocate (  temparray2(MNE,dofh,1), temparray3(MNE,dofh,1) )
      Allocate (  temparray4(MNE,dofh,1), temparray5(MNE,dofh,1) )

C     FIND THE MAXIMUM AND MINIMUM OF EACH VARIABLE OVER ALL ELEMENTS 
C     SHARING A NODE

      do k=1,ne

         if (dofs(k).gt.1) then

            do ll = 1,dofs(k)

               ze(ll,k,nrk+2) = ze(ll,k,irk+1)
               qx(ll,k,nrk+2) = qx(ll,k,irk+1)
               qy(ll,k,nrk+2) = qy(ll,k,irk+1)

#ifdef TRACE
               iota(ll,k,nrk+2) = iota(ll,k,irk+1)
#endif

#ifdef CHEM
               iota(ll,k,nrk+2) = iota(ll,k,irk+1)
               iota2(ll,k,nrk+2) = iota2(ll,k,irk+1)
#endif

            enddo

         elseif (dofs(k).eq.1) then

            cycle

         endif

      enddo   

c.....convert initial values to the taylor basis (multiply by nmatrix) on base element


      zetaylor = 0.d0 
      qxtaylor = 0.d0
      qytaylor = 0.d0

#ifdef TRACE
      iotataylor = 0.d0
#endif

#ifdef CHEM
      iotataylor = 0.d0
      iota2taylor = 0.d0
#endif

      do k=1,ne

                                !if (dofs(k).gt.1) then

         do ll=1,dofs(k)
            
            do ss=1,dofs(k)

               zetaylor(k,ll,1) = zetaylor(k,ll,1) + nmatrix(k,ll,ss,dofs(k)) * ze(ss,k,irk+1)
               qxtaylor(k,ll,1) = qxtaylor(k,ll,1) + nmatrix(k,ll,ss,dofs(k)) * qx(ss,k,irk+1)
               qytaylor(k,ll,1) = qytaylor(k,ll,1) + nmatrix(k,ll,ss,dofs(k)) * qy(ss,k,irk+1)

#ifdef TRACE
               iotataylor(k,ll,1) = iotataylor(k,ll,1) + 
     &              nmatrix(k,ll,ss,dofs(k))*iota(ss,k,irk+1)
#endif

#ifdef CHEM
               iotataylor(k,ll,1) = iotataylor(k,ll,1) + nmatrix(k,ll,ss,dofs(k))*iota(ss,k,irk+1)
               iota2taylor(k,ll,1) = iota2taylor(k,ll,1) + nmatrix(k,ll,ss,dofs(k))*iota2(ss,k,irk+1)
#endif
               
            enddo

         enddo

      enddo


#ifdef CMPI

      CALL UPDATER_ELEM_MOD2(ZEtaylor,QXtaylor,QYtaylor,1,3)

#ifdef TRACE
      CALL UPDATER_ELEM_MOD2(iotataylor,iotataylor,QYtaylor,1,2 )
#endif

#ifdef CHEM
      CALL UPDATER_ELEM_MOD2(iotataylor,iota2taylor,QYtaylor,1,2 )
#endif

#endif


      temparray1 = 0.D0
      temparray2 = 0.D0
      temparray3 = 0.D0

#ifdef TRACE
      temparray4 = 0.D0
#endif

#ifdef CHEM
      temparray4 = 0.D0
      temparray5 = 0.D0
#endif

      do ll=dofh,1,-1

         if(ll+floor(  0.5D0 + sqrt(real(2*ll)) )+1.le.dofh) then


            Call Slopelim_Aux(ZEtaylor(:,:,1),QXtaylor(:,:,1),QYtaylor(:,:,1),   
     &           iotataylor(:,:,1), iota2taylor(:,:,1),ll,
     &           ll+floor(  0.5D0 + sqrt(real(2*ll)) ),
     &           ll+floor(  0.5D0 + sqrt(real(2*ll)) )+1 )

            if (ll-1+floor(  0.5D0 + sqrt(real(2*(ll-1))) +1).eq.ll+floor(  0.5D0 + sqrt(real(2*(ll))) )) then

               temparray1(:,ll+floor(  0.5D0 + sqrt(real(2*ll))),1) = 
     &              ZEtaylor( :,ll+floor( 0.5D0 + sqrt(real(2*ll)) ),1 )
               temparray2(:,ll+floor(  0.5D0 + sqrt(real(2*ll))),1) = 
     &              QXtaylor( :,ll+floor( 0.5D0 + sqrt(real(2*ll)) ),1 )
               temparray3(:,ll+floor(  0.5D0 + sqrt(real(2*ll))),1) = 
     &              QYtaylor( :,ll+floor( 0.5D0 + sqrt(real(2*ll)) ),1 )

#ifdef TRACE
               temparray4(:,ll+floor(  0.5D0 + sqrt(real(2*ll))),1) = 
     &              iotataylor( :,ll+floor( 0.5D0 + sqrt(real(2*ll)) ),1 )
#endif

#ifdef CHEM
               temparray4(:,ll+floor(  0.5D0 + sqrt(real(2*ll))),1) = 
     &              iotataylor( :,ll+floor( 0.5D0 + sqrt(real(2*ll)) ),1 )
               temparray5(:,ll+floor(  0.5D0 + sqrt(real(2*ll))),1) = 
     &              iota2taylor( :,ll+floor( 0.5D0 + sqrt(real(2*ll)) ),1 )
#endif

            endif

            if (ll+floor(  0.5D0 + sqrt(real(2*(ll))) +1).eq.
     &           ll+1+floor(  0.5D0 + sqrt(real(2*(ll+1))) )) then

               do k = 1,ne

                  if (ZEtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll))+1 ),1 ).gt.0.D0.and.
     &                 temparray1(k,ll+1+floor(  0.5D0 + sqrt(real(2*(ll+1)))),1).gt.0.D0 ) then

                     ZEtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll))+1 ),1 ) = 
     &                    min(ZEtaylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) + 1 ),1),
     &                    temparray1(k,ll+1+floor(  0.5D0 + sqrt(real(2*(ll+1)))),1))

                  elseif (ZEtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) +1 ),1 ).lt.0.D0.and.
     &                    temparray1(k,ll+1+floor(  0.5D0 + sqrt(real(2*(ll+1)))),1).lt.0.D0 ) then
                     
                     ZEtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll))+1 ),1 ) =  
     &                    max(ZEtaylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) + 1),1),
     &                    temparray1(k,ll+1+floor(  0.5D0 + sqrt(real(2*(ll+1)))),1))

                  elseif (temparray1(k,ll+1+floor(  0.5D0 + sqrt(real(2*(ll+1)))),1).ne.0.D0) then

                     ZEtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) +1),1 ) = 0.D0

                  endif

                  if (QXtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll))+1 ),1 ).gt.0.D0.and.
     &                 temparray2(k,ll+1+floor(  0.5D0 + sqrt(real(2*(ll+1)))),1).gt.0.D0 ) then

                     QXtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll))+1 ),1 ) = 
     &                    min(QXtaylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) + 1 ),1),
     &                    temparray2(k,ll+1+floor(  0.5D0 + sqrt(real(2*(ll+1)))),1))

                  elseif (QXtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) +1 ),1 ).lt.0.D0.and.
     &                    temparray2(k,ll+1+floor(  0.5D0 + sqrt(real(2*(ll+1)))),1).lt.0.D0 ) then
                     
                     QXtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll))+1 ),1 ) =  
     &                    max(QXtaylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) + 1),1),
     &                    temparray2(k,ll+1+floor(  0.5D0 + sqrt(real(2*(ll+1)))),1))

                  elseif (temparray2(k,ll+1+floor(  0.5D0 + sqrt(real(2*(ll+1)))),1).ne.0.D0) then

                     QXtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) +1),1 ) = 0.D0

                  endif

                  if (QYtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll))+1 ),1 ).gt.0.D0.and.
     &                 temparray3(k,ll+1+floor(  0.5D0 + sqrt(real(2*(ll+1)))),1).gt.0.D0 ) then

                     QYtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll))+1 ),1 ) = 
     &                    min(QYtaylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) + 1 ),1),
     &                    temparray3(k,ll+1+floor(  0.5D0 + sqrt(real(2*(ll+1)))),1))

                  elseif (QYtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) +1 ),1 ).lt.0.D0.and.
     &                    temparray3(k,ll+1+floor(  0.5D0 + sqrt(real(2*(ll+1)))),1).lt.0.D0 ) then
                     
                     QYtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll))+1 ),1 ) =  
     &                    max(QYtaylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) + 1),1),
     &                    temparray3(k,ll+1+floor(  0.5D0 + sqrt(real(2*(ll+1)))),1))

                  elseif (temparray3(k,ll+1+floor(  0.5D0 + sqrt(real(2*(ll+1)))),1).ne.0.D0) then

                     QYtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) +1),1 ) = 0.D0

                  endif

#ifdef TRACE
                  if (iotataylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll))+1 ),1 ).gt.0.D0.and.
     &                 temparray4(k,ll+1+floor(  0.5D0 + sqrt(real(2*(ll+1)))),1).gt.0.D0 ) then

                     iotataylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll))+1 ),1 ) = 
     &                    min(iotataylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) + 1 ),1),
     &                    temparray4(k,ll+1+floor(  0.5D0 + sqrt(real(2*(ll+1)))),1))

                  elseif (iotataylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) +1 ),1 ).lt.0.D0.and.
     &                    temparray4(k,ll+1+floor(  0.5D0 + sqrt(real(2*(ll+1)))),1).lt.0.D0 ) then
                     
                     iotataylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll))+1 ),1 ) =  
     &                    max(iotataylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) + 1),1),
     &                    temparray4(k,ll+1+floor(  0.5D0 + sqrt(real(2*(ll+1)))),1))

                  elseif (temparray4(k,ll+1+floor(  0.5D0 + sqrt(real(2*(ll+1)))),1).ne.0.D0) then

                     iotataylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) +1),1 ) = 0.D0

                  endif
#endif

#ifdef CHEM
                  if (iotataylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll))+1 ),1 ).gt.0.D0.and.
     &                 temparray4(k,ll+1+floor(  0.5D0 + sqrt(real(2*(ll+1)))),1).gt.0.D0 ) then

                     iotataylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll))+1 ),1 ) = 
     &                    min(iotataylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) + 1 ),1),
     &                    temparray4(k,ll+1+floor(  0.5D0 + sqrt(real(2*(ll+1)))),1))

                  elseif (iotataylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) +1 ),1 ).lt.0.D0.and.
     &                    temparray4(k,ll+1+floor(  0.5D0 + sqrt(real(2*(ll+1)))),1).lt.0.D0 ) then
                     
                     iotataylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll))+1 ),1 ) =  
     &                    max(iotataylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) + 1),1),
     &                    temparray4(k,ll+1+floor(  0.5D0 + sqrt(real(2*(ll+1)))),1))

                  elseif (temparray4(k,ll+1+floor(  0.5D0 + sqrt(real(2*(ll+1)))),1).ne.0.D0) then

                     iotataylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) +1),1 ) = 0.D0

                  endif

                  if (iota2taylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll))+1 ),1 ).gt.0.D0.and.
     &                 temparray5(k,ll+1+floor(  0.5D0 + sqrt(real(2*(ll+1)))),1).gt.0.D0 ) then

                     iota2taylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll))+1 ),1 ) = 
     &                    min(iota2taylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) + 1 ),1),
     &                    temparray5(k,ll+1+floor(  0.5D0 + sqrt(real(2*(ll+1)))),1))

                  elseif (iota2taylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) +1 ),1 ).lt.0.D0.and.
     &                    temparray5(k,ll+1+floor(  0.5D0 + sqrt(real(2*(ll+1)))),1).lt.0.D0 ) then
                     
                     iota2taylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll))+1 ),1 ) =  
     &                    max(iota2taylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) + 1),1),
     &                    temparray5(k,ll+1+floor(  0.5D0 + sqrt(real(2*(ll+1)))),1))

                  elseif (temparray5(k,ll+1+floor(  0.5D0 + sqrt(real(2*(ll+1)))),1).ne.0.D0) then

                     iota2taylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) +1),1 ) = 0.D0

                  endif
#endif

               enddo

            endif

         endif
         
      enddo

#ifdef CMPI

      CALL UPDATER_ELEM_MOD2(ZEtaylor,QXtaylor,QYtaylor,1,3)

#ifdef TRACE
      CALL UPDATER_ELEM_MOD2(iotataylor,iotataylor,QYtaylor,1,2 )
#endif

#ifdef CHEM
      CALL UPDATER_ELEM_MOD2(iotataylor,iota2taylor,QYtaylor,1,2 )
#endif

#endif

C.....Transform back to the Dubiner basis (multiply by NmatrixInv),

      ZEconst =  0.D0
      QXconst =  0.D0
      QYconst =  0.D0

#ifdef TRACE
      iotaconst =  0.D0
#endif

#ifdef CHEM
      iotaconst =  0.D0
      iota2const =  0.D0
#endif

      do k=1,ne

         if (dofs(k).gt.1) then
                                !do lll=1,3

            do ll=1,dofs(k)

               do ss=1,dofs(k)

                  ZEconst(k,ll) = ZEconst(k,ll) + NmatrixInv(k,ll,ss,dofs(k)) 
     &                 * ZEtaylor(k,ss,1)
                  QXconst(k,ll) = QXconst(k,ll) + NmatrixInv(k,ll,ss,dofs(k)) 
     &                 * QXtaylor(k,ss,1)
                  QYconst(k,ll) = QYconst(k,ll) + NmatrixInv(k,ll,ss,dofs(k)) 
     &                 * QXtaylor(k,ss,1)

#ifdef TRACE
                  iotaconst(k,ll) = iotaconst(k,ll) + 
     &                 NmatrixInv(k,ll,ss,dofs(k)) * iotataylor(k,ss,1) 
#endif

#ifdef CHEM
                  iotaconst(k,ll) = iotaconst(k,ll) + 
     &                 NmatrixInv(k,ll,ss,dofs(k)) *  iotataylor(k,ss,1)
                  iota2const(k,ll) = iota2const(k,ll) + 
     &                 NmatrixInv(k,ll,ss,dofs(k)) * iota2taylor(k,ss,1)
#endif

               enddo

            enddo

         elseif (dofs(k).eq.1) then

            cycle

         endif

      enddo

C.....Set limit values

      do k=1,ne

         if (dofs(k).gt.1) then

            do ll = 1,dofs(k)

               ZE(ll,k,irk+1) = ZEconst(k,ll)
               QX(ll,k,irk+1) = QXconst(k,ll)
               QY(ll,k,irk+1) = QYconst(k,ll)

#ifdef TRACE
               iota(ll,k,irk+1) = iotaconst(k,ll)
#endif

#ifdef CHEM
               iota(ll,k,irk+1) = iotaconst(k,ll)
               iota2(ll,k,irk+1) = iota2const(k,ll)
#endif

            enddo

         elseif (dofs(k).eq.1) then

            cycle

         endif

      enddo
      
      return
      end subroutine



C***********************************
C     *
C     SUBROUTINE Slopelim_Aux()    *
C     -cem
C     *
C***********************************

      SUBROUTINE Slopelim_Aux(ZEder,QXder,QYder,
     &     iotader,iota2der,l1,l2,l3)

C.....Use appropriate modules

      USE SIZES, ONLY : SZ
      USE GLOBAL
      USE DG

#ifdef CMPI
      USE MESSENGER
#endif

      IMPLICIT NONE

C.....Declare local variables

      INTEGER L,LL,INC1,INC2,INC3,KDP,NN,IVAR,I,J,kk,k,varnum,ss,lll,l1,l2,l3
      REAL(SZ) ZEC(3),ZEVERTEX(3),DIF(3),SUMLOC,SUMDIF,SIGNDIF,ZEVERTEX2(3),
     $     DIV,REDFAC,REDMAX,zek(MNE,3,1),zecc(MNE,3),zeve(MNE,3)
      Real(SZ), intent(inout) :: ZEder(MNE,dofh,1)
      Real(SZ), intent(inout) :: QXder(MNE,dofh,1)
      Real(SZ), intent(inout) :: QYder(MNE,dofh,1)
      Real(SZ), intent(inout) :: iotader(MNE,dofh,1)
      Real(SZ), intent(inout) :: iota2der(MNE,dofh,1)
      Real(SZ) ZEMIN1(3),ZEMAX1(3),QXMIN1(3),QXMAX1(3)
      Real(SZ) QYMIN1(3),QYMAX1(3)
      Real(SZ) iotaMIN1(3),iotaMAX1(3),zztop
      Real(SZ) iota2MIN1(3),iota2MAX1(3)
      Real(SZ), Allocatable :: ZE_MIN1(:),ZE_MAX1(:),QX_MIN1(:),QX_MAX1(:)
      Real(SZ), Allocatable :: QY_MIN1(:),QY_MAX1(:)
      Real(SZ), Allocatable :: iota_MIN1(:),iota_MAX1(:), temp2(:,:)
      Real(SZ), Allocatable :: iota2_MIN1(:),iota2_MAX1(:)
      Real(SZ), Allocatable :: Zhi(:,:),ZhiInv(:,:), A(:,:),tempa(:,:)

      Allocate ( ZE_MIN1(NP),ZE_MAX1(NP),QX_MIN1(NP) )
      Allocate ( QY_MIN1(NP),QY_MAX1(NP),QX_MAX1(NP) )
      Allocate ( iota_MIN1(NP),iota_MAX1(NP) )
      Allocate ( iota2_MIN1(NP),iota2_MAX1(NP), temp2(3,3) )
      Allocate ( Zhi(3,3),ZhiInv(3,3),A(3,3),tempa(3,3) )



C     FIND THE MAXIMUM AND MINIMUM OF EACH VARIABLE OVER ALL ELEMENTS 
C     SHARING A NODE


      DO I = 1,NP
         ZE_MIN1(I)=99999.
         ZE_MAX1(I)=-99999.
         QX_MIN1(I)=99999.
         QX_MAX1(I)=-99999.
         QY_MIN1(I)=99999.
         QY_MAX1(I)=-99999.

#ifdef TRACE
         iota_MIN1(I)=99999.
         iota_MAX1(I)=-99999.
#endif

#ifdef CHEM
         iota_MIN1(I)=99999.
         iota_MAX1(I)=-99999.
         iota2_MIN1(I)=99999.
         iota2_MAX1(I)=-99999.
#endif

         NO_NBORS = EL_COUNT(I)

         DO J = 1,NO_NBORS
            NBOR_EL = ELETAB(I,1+J)

c     IF(WDFLG(NBOR_EL).EQ.0) CYCLE ! DON'T COUNT DRY ELEMENTS  sb 02/26/07

            ZE_DG(J) = ZEder(NBOR_EL,l1,1)
            QX_DG(J) = QXder(NBOR_EL,l1,1)
            QY_DG(J) = QYder(NBOR_EL,l1,1)

#ifdef TRACE
            iota_DG(J) = iotader(NBOR_EL,l1,1)
#endif

#ifdef CHEM
            iota_DG(J) = iotader(NBOR_EL,l1,1)
            iota2_DG(J) = iota2der(NBOR_EL,l1,1)
#endif

C     
            IF (ZE_DG(J).LT.ZE_MIN1(I))THEN
               ZE_MIN1(I)=ZE_DG(J)
            ENDIF
            IF (ZE_DG(J).GT.ZE_MAX1(I)) THEN
               ZE_MAX1(I)=ZE_DG(J)
            ENDIF
            IF (QX_DG(J).LT.QX_MIN1(I))THEN
               QX_MIN1(I)=QX_DG(J)
            ENDIF
            IF (QX_DG(J).GT.QX_MAX1(I)) THEN
               QX_MAX1(I)=QX_DG(J)
            ENDIF
            IF (QY_DG(J).LT.QY_MIN1(I))THEN
               QY_MIN1(I)=QY_DG(J)
            ENDIF
            IF (QY_DG(J).GT.QY_MAX1(I)) THEN
               QY_MAX1(I)=QY_DG(J)
            ENDIF

#ifdef TRACE
            IF (iota_DG(J).LT.iota_MIN1(I))THEN
               iota_MIN1(I)=iota_DG(J)
            ENDIF

            IF (iota_DG(J).GT.iota_MAX1(I)) THEN
               iota_MAX1(I)=iota_DG(J)
            ENDIF
#endif

#ifdef CHEM
            IF (iota_DG(J).LT.iota_MIN1(I))THEN
               iota_MIN1(I)=iota_DG(J)
            ENDIF


            IF (iota_DG(J).GT.iota_MAX1(I)) THEN
               iota_MAX1(I)=iota_DG(J)
            ENDIF

            IF (iota2_DG(J).LT.iota2_MIN1(I))THEN
               iota2_MIN1(I)=iota2_DG(J)
            ENDIF

            IF (iota2_DG(J).GT.iota2_MAX1(I)) THEN
               iota2_MAX1(I)=iota2_DG(J)
            ENDIF
#endif

         ENDDO
      ENDDO
#ifdef CMPI

      CALL UPDATER(ZE_MIN1,ZE_MAX1,QX_MIN1,3)
      CALL UPDATER(QX_MAX1,QY_MIN1,QY_MAX1,3)

#ifdef TRACE
      CALL UPDATER(iota_MAX1,iota_MIN1,QY_MAX1,2)
#endif

#ifdef CHEM
      CALL UPDATER(iota_MAX1,iota_MIN1,QY_MAX1,2)
      CALL UPDATER(iota2_MAX1,iota2_MIN1,QY_MAX1,2)
#endif

#endif
C     
C     LOOP OVER ELEMENTS TO CALCULATE NEW VERTEX VALUES
C     

      zek = 0.D0
      zecc = 0.D0
      zeve = 0.D0

      DO I=1,NE 

         !IF(WDFLG(I).EQ.0) CYCLE ! DON'T COUNT DRY ELEMENTS  sb 02/26/07
         N1=NM(I,1)
         N2=NM(I,2)
         N3=NM(I,3)

         varnum = 3

#ifdef TRACE
         varnum = 4
#endif

#ifdef CHEM
         varnum = 5
#endif


         DO IVAR=1,varnum

            IF (IVAR.EQ.1) THEN
               ZEC(1)=ZEder(i,l1,1)
               ZEC(2)=ZEder(i,l2,1)
               ZEC(3)=ZEder(i,l3,1)
               ZEMAX1(1)=ZE_MAX1(N1)
               ZEMIN1(1)=ZE_MIN1(N1)
               ZEMAX1(2)=ZE_MAX1(N2)
               ZEMIN1(2)=ZE_MIN1(N2)
               ZEMAX1(3)=ZE_MAX1(N3)
               ZEMIN1(3)=ZE_MIN1(N3)
            ENDIF

            IF (IVAR.EQ.2) THEN
               ZEC(1)=QXder(i,l1,1)
               ZEC(2)=QXder(i,l2,1)
               ZEC(3)=QXder(i,l3,1)
               ZEMAX1(1)=QX_MAX1(N1)
               ZEMIN1(1)=QX_MIN1(N1)
               ZEMAX1(2)=QX_MAX1(N2)
               ZEMIN1(2)=QX_MIN1(N2)
               ZEMAX1(3)=QX_MAX1(N3)
               ZEMIN1(3)=QX_MIN1(N3)
            ENDIF

            IF (IVAR.EQ.3) THEN
               ZEC(1)=QYder(i,l1,1)
               ZEC(2)=QYder(i,l2,1)
               ZEC(3)=QYder(i,l3,1)
               ZEMAX1(1)=QY_MAX1(N1)
               ZEMIN1(1)=QY_MIN1(N1)
               ZEMAX1(2)=QY_MAX1(N2)
               ZEMIN1(2)=QY_MIN1(N2)
               ZEMAX1(3)=QY_MAX1(N3)
               ZEMIN1(3)=QY_MIN1(N3)
            ENDIF

            IF (IVAR.EQ.4) THEN
               ZEC(1)=iotader(i,l1,1)
               ZEC(2)=iotader(i,l2,1)
               ZEC(3)=iotader(i,l3,1)
               ZEMAX1(1)=iota_MAX1(N1)
               ZEMIN1(1)=iota_MIN1(N1)
               ZEMAX1(2)=iota_MAX1(N2)
               ZEMIN1(2)=iota_MIN1(N2)
               ZEMAX1(3)=iota_MAX1(N3)
               ZEMIN1(3)=iota_MIN1(N3)
            ENDIF

            IF (IVAR.EQ.5) THEN
               ZEC(1)=iota2der(i,l1,1)
               ZEC(2)=iota2der(i,l2,1)
               ZEC(3)=iota2der(i,l3,1)
               ZEMAX1(1)=iota2_MAX1(N1)
               ZEMIN1(1)=iota2_MIN1(N1)
               ZEMAX1(2)=iota2_MAX1(N2)
               ZEMIN1(2)=iota2_MIN1(N2)
               ZEMAX1(3)=iota2_MAX1(N3)
               ZEMIN1(3)=iota2_MIN1(N3)
            ENDIF


            ZEVERTEX = 0.D0
            zevertex2 = 0.D0

            Zhi = 0.D0
            ZhiInv = 0.D0

            do lll = 1,3

               Zhi(lll,1) =  var2sigmav(i,lll,1) 
               Zhi(lll,2) =  var2sigmav(i,lll,2) 
               Zhi(lll,3) =  var2sigmav(i,lll,3) 

            enddo

            DO KK=1,3
               ZEVERTEX(1)=ZEVERTEX(1)+ Zhi(1,kk) * ZEC(KK)
               ZEVERTEX(2)=ZEVERTEX(2)+ Zhi(2,kk) * ZEC(KK)
               ZEVERTEX(3)=ZEVERTEX(3)+ Zhi(3,kk) * ZEC(KK)
            ENDDO

            ZEVERTEX2(1)=iota(1,I,IRK+1)
            ZEVERTEX2(2)=iota(1,I,IRK+1)
            ZEVERTEX2(3)=iota(1,I,IRK+1)

            DO KK=2,3
               ZEVERTEX2(1)=ZEVERTEX2(1)+ PHI_CORNER(KK,1,1)*iota(kk,I,IRK+1)
               ZEVERTEX2(2)=ZEVERTEX2(2)+ PHI_CORNER(KK,2,1)*iota(kk,I,IRK+1)
               ZEVERTEX2(3)=ZEVERTEX2(3)+ PHI_CORNER(KK,3,1)*iota(kk,I,IRK+1)
            ENDDO

C     
C     RESET THE VERTEX VALUE TO BE LESS THAN OR EQUAL TO THE MAX AND
C     GREATER THAN OR EQUAL TO THE MIN AT THAT VERTEX
C     
            ZEVERTEX(1)=DMAX1(DMIN1(ZEVERTEX(1),ZEMAX1(1)),ZEMIN1(1))
            ZEVERTEX(2)=DMAX1(DMIN1(ZEVERTEX(2),ZEMAX1(2)),ZEMIN1(2))
            ZEVERTEX(3)=DMAX1(DMIN1(ZEVERTEX(3),ZEMAX1(3)),ZEMIN1(3))

C     LOOP OVER THE VERTICES 3 TIMES
C     IF THE VALUE AT THE VERTEX IS ABOVE (BELOW) THE MAX (MIN) AT THAT
C     VERTEX THEN SUBTRACT OFF THE DIFFERENCE AND ADD IT TO THE OTHER
C     VERTICES
C     
            DO LL=1,3
               SUMLOC=(ZEVERTEX(1)+ZEVERTEX(2)+ZEVERTEX(3))/3.0D0
               SUMDIF=(SUMLOC-ZEC(1))*3.0D0
               SIGNDIF=DSIGN(1.D0,SUMDIF)
               DIF(1)=(ZEVERTEX(1)-ZEC(1))*SIGNDIF
               DIF(2)=(ZEVERTEX(2)-ZEC(1))*SIGNDIF
               DIF(3)=(ZEVERTEX(3)-ZEC(1))*SIGNDIF
               INC1=0
               IF (DIF(1).GT.0) INC1=1
               INC2=0
               IF (DIF(2).GT.0) INC2=1
               INC3=0
               IF (DIF(3).GT.0) INC3=1
               KDP=INC1+INC2+INC3
C     
               DO K=1,3
                  DIV=DMAX1(1.D0,DFLOAT(KDP))
                  IF (DIF(K).GT.0) THEN
                     REDFAC=SUMDIF*SIGNDIF/DIV
                     KDP=KDP-1
                  ELSE
                     REDFAC=0
                  ENDIF
                  IF (SIGNDIF.GT.0) THEN
                     REDMAX=ZEVERTEX(K)-ZEMIN1(K)
                  ELSE
                     REDMAX=ZEMAX1(K)-ZEVERTEX(K)
                  ENDIF
                  REDFAC=DMIN1(REDFAC,REDMAX)
                  SUMDIF=SUMDIF-REDFAC*SIGNDIF
                  ZEVERTEX(K)=ZEVERTEX(K)-REDFAC*SIGNDIF
               ENDDO
            ENDDO
            
            Call Inv2(Zhi, ZhiInv, 3)

            IF (IVAR.EQ.1) THEN

               ZEder(i,l2,1)= 0.D0
               ZEder(i,l3,1)= 0.D0

               do lll = 1,3


                  ZEder(i,l2,1) =  ZEder(i,l2,1) + ZhiInv(2,lll) 
     &                 * ZEVERTEX(lll)

                  ZEder(i,l3,1) =  ZEder(i,l3,1) + ZhiInv(3,lll) 
     &                 * ZEVERTEX(lll)

               enddo

               ZEder(i,l1,1) = ZEder(i,l1,1)

            ENDIF

            IF (IVAR.EQ.2) THEN

               QXder(i,l2,1)= 0.D0
               QXder(i,l3,1)= 0.D0

               do lll = 1,3

                  QXder(i,l2,1) =  QXder(i,l2,1) + ZhiInv(2,lll) 
     &                 * ZEVERTEX(lll)

                  QXder(i,l3,1) =  QXder(i,l3,1) + ZhiInv(3,lll) 
     &                 * ZEVERTEX(lll)

               enddo

               QXder(i,l1,1) = QXder(i,l1,1)

            ENDIF

            IF (IVAR.EQ.3) THEN

               QYder(i,l2,1)= 0.D0
               QYder(i,l3,1)= 0.D0

               do lll = 1,3

                  QYder(i,l2,1) =  QYder(i,l2,1) + ZhiInv(2,lll) 
     &                 * ZEVERTEX(lll)

                  QYder(i,l3,1) =  QYder(i,l3,1) + ZhiInv(3,lll) 
     &                 * ZEVERTEX(lll)

               enddo

               QYder(i,l1,1) = QYder(i,l1,1)

            ENDIF

            IF (IVAR.EQ.4) THEN

                                !Call Inv2(Zhi, ZhiInv, 3) 
               
               iotader(i,l2,1)= 0.D0
               iotader(i,l3,1)= 0.D0

               do lll = 1,3

                  iotader(i,l2,1) =  iotader(i,l2,1) + ZhiInv(2,lll) 
     &                 * ZEVERTEX(lll)

                  iotader(i,l3,1) =  iotader(i,l3,1) + ZhiInv(3,lll) 
     &                 * ZEVERTEX(lll)

               enddo

               iotader(i,l1,1) = iotader(i,l1,1)

            ENDIF

            IF (IVAR.EQ.5) THEN
               
               iota2der(i,l2,1)= 0.D0
               iota2der(i,l3,1)= 0.D0
               
               do lll = 1,3
                  
                  iota2der(i,l2,1) =  iota2der(i,l2,1) + ZhiInv(2,lll) 
     &                 * ZEVERTEX(lll)
                  
                  iota2der(i,l3,1) =  QYder(i,l3,1) + ZhiInv(3,lll) 
     &                 * ZEVERTEX(lll)
                  
               enddo
               
               iota2der(i,l1,1) = iota2der(i,l1,1)

            ENDIF
            
         ENDDO

      ENDDO
      
      RETURN
      END SUBROUTINE 

c$$$C*******************************************************************************
c$$$C     
c$$$C     SUBROUTINE SLOPELIMITER9()old
c$$$C     
c$$$C     Written 27 Oct 2009
c$$$C     
c$$$C     This subroutine selects an adapted (2) vertex based slope limiter based on
c$$$C     a Taylor Polynomial basis, and is consistent with p_adaptation.F
c$$$C
c$$$C     This needs to be fixed
c$$$C     
c$$$C*******************************************************************************
c$$$
c$$$      SUBROUTINE SLOPELIMITER9()
c$$$
c$$$C.....Use appropriate modules
c$$$
c$$$      USE GLOBAL
c$$$      USE DG
c$$$
c$$$#ifdef CMPI
c$$$      USE MESSENGER
c$$$      USE MESSENGER_ELEM
c$$$#endif
c$$$
c$$$      IMPLICIT NONE
c$$$
c$$$      Integer k,ll,ss,lll,ell,bb
c$$$
c$$$      REAL(SZ) fd,marea
c$$$
c$$$C.....We work over the master element
c$$$C.....Set initial values
c$$$
c$$$      fd = slope_weight         ! add weight for lower order pieces (fd<1 => stronger limiting)
c$$$      
c$$$
c$$$      DO k=1,NE
c$$$
c$$$         if (dofs(k).gt.1) then
c$$$
c$$$            DO ll = 1,dofs(k)
c$$$
c$$$               ZE(ll,k,NRK+2) = ZE(ll,k,IRK+1)
c$$$               QX(ll,k,NRK+2) = QX(ll,k,IRK+1)
c$$$               QY(ll,k,NRK+2) = QY(ll,k,IRK+1)
c$$$
c$$$               if (tracer_flag.eq.1) then
c$$$
c$$$                  iota(ll,k,NRK+2) = iota(ll,k,IRK+1)
c$$$
c$$$               endif
c$$$
c$$$               if (chem_flag.eq.1) then
c$$$
c$$$                  iota(ll,k,NRK+2) = iota(ll,k,IRK+1)
c$$$                  iota2(ll,k,NRK+2) = iota2(ll,k,IRK+1)
c$$$
c$$$               endif
c$$$
c$$$
c$$$            ENDDO
c$$$
c$$$         elseif (dofs(k).eq.1) then
c$$$
c$$$            cycle
c$$$
c$$$         endif
c$$$
c$$$      ENDDO   
c$$$
c$$$      marea = 2.D0              !master elements area
c$$$
c$$$C.....Convert initial values to the Taylor basis (multiply by Nmatrix) on base element
c$$$
c$$$
c$$$      ZEtaylor = 0.D0 
c$$$      QXtaylor = 0.D0
c$$$      QYtaylor = 0.D0
c$$$
c$$$      if (tracer_flag.eq.1.or.chem_flag.eq.1) then
c$$$
c$$$         iotataylor = 0.D0
c$$$         iota2taylor = 0.D0
c$$$
c$$$      endif
c$$$
c$$$      do k=1,ne
c$$$
c$$$         if (dofs(k).gt.1) then
c$$$
c$$$            do ll=1,dofs(k)
c$$$               
c$$$               do ss=1,dofs(k)
c$$$
c$$$                  ZEtaylor(k,ll,1) = ZEtaylor(k,ll,1) + Nmatrix(k,ll,ss,dofs(k)) * ZE(ss,k,nrk+2)
c$$$                  QXtaylor(k,ll,1) = QXtaylor(k,ll,1) + Nmatrix(k,ll,ss,dofs(k))* QX(ss,k,nrk+2)
c$$$                  QYtaylor(k,ll,1) = QYtaylor(k,ll,1) + Nmatrix(k,ll,ss,dofs(k)) * QY(ss,k,nrk+2)
c$$$
c$$$                  if (tracer_flag.eq.1) then
c$$$
c$$$                     iotataylor(k,ll,1) = iotataylor(k,ll,1) + Nmatrix(k,ll,ss,dofs(k)) * iota(ss,k,nrk+2)
c$$$                     
c$$$                  endif
c$$$
c$$$                  if (chem_flag.eq.1) then
c$$$
c$$$                     iotataylor(k,ll,1) = iotataylor(k,ll,1) + Nmatrix(k,ll,ss,dofs(k)) * iota(ss,k,nrk+2)
c$$$                     iota2taylor(k,ll,1) = iota2taylor(k,ll,1) + Nmatrix(k,ll,ss,dofs(k)) * iota2(ss,k,nrk+2)
c$$$
c$$$                  endif
c$$$                  
c$$$               enddo
c$$$
c$$$            enddo
c$$$
c$$$         elseif (dofs(k).eq.1) then
c$$$
c$$$            cycle
c$$$
c$$$         endif
c$$$
c$$$      enddo
c$$$
c$$$C.....Find values at vertices of base elements and neighbors
c$$$
c$$$
c$$$      ZEmax = -100.D0
c$$$      QXmax = -100.D0
c$$$      QYmax = -100.D0
c$$$      ZEmin = 100.D0
c$$$      QXmin = 100.D0
c$$$      QYmin = 100.D0
c$$$
c$$$      if (tracer_flag.eq.1) then
c$$$
c$$$         iotamax = -100.D0
c$$$         iotamin = 100.D0
c$$$
c$$$      endif
c$$$
c$$$      if (chem_flag.eq.1) then
c$$$
c$$$         iotamax = -100.D0
c$$$         iota2max = -100.D0
c$$$         iotamin = 100.D0
c$$$         iota2min = 100.D0
c$$$
c$$$      endif
c$$$
c$$$
c$$$#ifdef CMPI
c$$$
c$$$      CALL UPDATER_ELEM_MOD2(ZEtaylor,QXtaylor,QYtaylor,1,3)
c$$$
c$$$      if(tracer_flag.eq.1.or.chem_flag.eq.1) then
c$$$
c$$$         CALL UPDATER_ELEM_MOD2(iotataylor,iota2taylor,QYtaylor,1,2 )
c$$$
c$$$      endif
c$$$
c$$$#endif
c$$$
c$$$
c$$$
c$$$      !Can likely omit this pre-call
c$$$
c$$$#ifdef CMPI
c$$$
c$$$      CALL UPDATERV(ZEmin)
c$$$      CALL UPDATERV(ZEmax)
c$$$      CALL UPDATERV(QXmin)
c$$$      CALL UPDATERV(QXmax)
c$$$      CALL UPDATERV(QYmin)
c$$$      CALL UPDATERV(QYmax)
c$$$
c$$$      if (tracer_flag.eq.1) then
c$$$
c$$$         CALL UPDATERV(iotamax)
c$$$         CALL UPDATERV(iotamin)
c$$$
c$$$      endif
c$$$
c$$$      if (chem_flag.eq.1) then
c$$$
c$$$         CALL UPDATERV(iotamax)
c$$$         CALL UPDATERV(iotamin)
c$$$         CALL UPDATERV(iota2max)
c$$$         CALL UPDATERV(iota2min)
c$$$
c$$$      endif
c$$$
c$$$#endif
c$$$
c$$$      do ell=1,mnp
c$$$
c$$$         do ll=1,minval(dofs(neigh_elem(ell,1:nneigh_elem(ell))))
c$$$
c$$$C.....Find max and min values over polynomial coefficients
c$$$
c$$$            ZEmax(ell,ll) = max(maxval( ZEtaylor(NEIGH_ELEM(ell,1:nneigh_elem(ell)),ll,1) )
c$$$     &           , ZEmax(ell,ll))
c$$$            QXmax(ell,ll) = max(maxval( QXtaylor(NEIGH_ELEM(ell,1:nneigh_elem(ell)),ll,1) )
c$$$     &           , QXmax(ell,ll))
c$$$            QYmax(ell,ll) = max(maxval( QYtaylor(NEIGH_ELEM(ell,1:nneigh_elem(ell)),ll,1) )
c$$$     &           , QYmax(ell,ll))
c$$$            ZEmin(ell,ll) = min(minval( ZEtaylor(NEIGH_ELEM(ell,1:nneigh_elem(ell)),ll,1) )
c$$$     &           , ZEmin(ell,ll))
c$$$            QXmin(ell,ll) = min(minval( QXtaylor(NEIGH_ELEM(ell,1:nneigh_elem(ell)),ll,1) )
c$$$     &           , QXmin(ell,ll))
c$$$            QYmin(ell,ll) = min(minval( QYtaylor(NEIGH_ELEM(ell,1:nneigh_elem(ell)),ll,1) )
c$$$     &           , QYmin(ell,ll))
c$$$
c$$$            if (tracer_flag.eq.1) then
c$$$
c$$$               iotamax(ell,ll) = max(maxval( iotataylor(NEIGH_ELEM(ell,1:nneigh_elem(ell)),ll,1) )
c$$$     &              , iotamax(ell,ll))
c$$$               iotamin(ell,ll) = min(minval( iotataylor(NEIGH_ELEM(ell,1:nneigh_elem(ell)),ll,1) )
c$$$     &              , iotamin(ell,ll))
c$$$
c$$$
c$$$            endif
c$$$
c$$$            if (chem_flag.eq.1) then
c$$$
c$$$               iotamax(ell,ll) = max(maxval( iotataylor(NEIGH_ELEM(ell,1:nneigh_elem(ell)),ll,1) )
c$$$     &              , iotamax(ell,ll))
c$$$               iota2max(ell,ll) = max(maxval( iota2taylor(NEIGH_ELEM(ell,1:nneigh_elem(ell)),ll,1) )
c$$$     &              , iota2max(ell,ll))
c$$$               iotamin(ell,ll) = min(minval( iotataylor(NEIGH_ELEM(ell,1:nneigh_elem(ell)),ll,1) )
c$$$     &              , iotamax(ell,ll))
c$$$               iota2min(ell,ll) = min(minval( iota2taylor(NEIGH_ELEM(ell,1:nneigh_elem(ell)),ll,1) )
c$$$     &              , iota2max(ell,ll))
c$$$               
c$$$            endif
c$$$            
c$$$         enddo
c$$$
c$$$      enddo
c$$$
c$$$
c$$$#ifdef CMPI
c$$$
c$$$
c$$$      CALL UPDATERV(ZEmin)
c$$$      CALL UPDATERV(ZEmax)
c$$$      CALL UPDATERV(QXmin)
c$$$      CALL UPDATERV(QXmax)
c$$$      CALL UPDATERV(QYmin)
c$$$      CALL UPDATERV(QYmax)
c$$$
c$$$      if (tracer_flag.eq.1) then
c$$$
c$$$         CALL UPDATERV(iotamax)
c$$$         CALL UPDATERV(iotamin)
c$$$
c$$$      endif
c$$$
c$$$      if (chem_flag.eq.1) then
c$$$
c$$$         CALL UPDATERV(iotamax)
c$$$         CALL UPDATERV(iotamin)
c$$$         CALL UPDATERV(iota2max)
c$$$         CALL UPDATERV(iota2min)
c$$$
c$$$      endif
c$$$
c$$$#endif
c$$$
c$$$C.....Must generate linear recostructions at vertices
c$$$
c$$$      ZEtaylorvert = 0.D0
c$$$      QXtaylorvert = 0.D0
c$$$      Qytaylorvert = 0.D0
c$$$
c$$$      if (tracer_flag.eq.1) then
c$$$
c$$$         iotataylorvert = 0.D0
c$$$
c$$$      endif
c$$$
c$$$      if (chem_flag.eq.1) then
c$$$
c$$$         iotataylorvert = 0.D0
c$$$         iota2taylorvert = 0.D0
c$$$
c$$$      endif
c$$$
c$$$      do k=1,ne
c$$$
c$$$         if (dofs(k).gt.1) then
c$$$
c$$$            do ll=1,dofs(k)
c$$$
c$$$               do lll=1,3
c$$$
c$$$                  if (ll.eq.1) then
c$$$
c$$$                     ZEtaylorvert(k,ll,lll) = ZEtaylor(k,ll,1) + ZEtaylor(k,ll+1,1)*( xi2vert(k,lll) -xi2BCb(k) )
c$$$     &                    + ZEtaylor(k,ll+2,1)*( xi1vert(k,lll) - xi1BCb(k) )
c$$$                     QXtaylorvert(k,ll,lll) = QXtaylor(k,ll,1) + QXtaylor(k,ll+1,1)*( xi2vert(k,lll) -xi2BCb(k) )
c$$$     &                    + QXtaylor(k,ll+2,1)*( xi1vert(k,lll) - xi1BCb(k) )
c$$$                     Qytaylorvert(k,ll,lll) = QYtaylor(k,ll,1) + QYtaylor(k,ll+1,1)*( xi2vert(k,lll) -xi2BCb(k) ) 
c$$$     &                    + QYtaylor(k,ll+2,1)*( xi1vert(k,lll) - xi1BCb(k) )
c$$$
c$$$                     if (tracer_flag.eq.1) then
c$$$
c$$$                        iotataylorvert(k,ll,lll) = iotataylor(k,ll,1) + 
c$$$     &                       iotataylor(k,ll+1,1)*( xi2vert(k,lll) - xi2BCb(k) )
c$$$     &                       + iotataylor(k,ll+2,1)*( xi1vert(k,lll) - xi1BCb(k) )
c$$$
c$$$                     endif
c$$$
c$$$                     if (chem_flag.eq.1) then
c$$$
c$$$                        iotataylorvert(k,ll,lll) = iotataylor(k,ll,1) + 
c$$$     &                       iotataylor(k,ll+1,1)*( xi2vert(k,lll) - xi2BCb(k) )
c$$$     &                       + iotataylor(k,ll+2,1)*( xi1vert(k,lll) - xi1BCb(k) )
c$$$                        iota2taylorvert(k,ll,lll) = iota2taylor(k,ll,1) + 
c$$$     &                       iota2taylor(k,ll+1,1)*( xi2vert(k,lll) - xi2BCb(k) )
c$$$     &                       + iota2taylor(k,ll+2,1)*( xi1vert(k,lll) - xi1BCb(k) )
c$$$
c$$$                     endif
c$$$
c$$$                  elseif (ll.gt.3) then
c$$$
c$$$                     ZEtaylorvert(k,ll,lll) = ZEtaylor(k,ll,1) + 
c$$$     &                    ZEtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) ),1)*( xi2vert(k,lll) - xi2BCb(k) ) 
c$$$     &                    + ZEtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) )+1,1)*( xi1vert(k,lll) - xi1BCb(k) )
c$$$                     QXtaylorvert(k,ll,lll) = QXtaylor(k,ll,1) + 
c$$$     &                    QXtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) ),1)*( xi2vert(k,lll) - xi2BCb(k) ) 
c$$$     &                    + QXtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) )+1,1)*( xi1vert(k,lll) - xi1BCb(k) )
c$$$                     Qytaylorvert(k,ll,lll) = QYtaylor(k,ll,1) + 
c$$$     &                    QYtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) ),1)*( xi2vert(k,lll) - xi2BCb(k) ) 
c$$$     &                    + QYtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) )+1,1)*( xi1vert(k,lll) - xi1BCb(k) )
c$$$
c$$$                     if (tracer_flag.eq.1) then
c$$$
c$$$                        iotataylorvert(k,ll,lll) = iotataylor(k,ll,1)+
c$$$     &                       iotataylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) ),1)*( xi2vert(k,lll) - xi2BCb(k))
c$$$     &                       + iotataylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) )+1,1)*( xi1vert(k,lll) - xi1BCb(k) )
c$$$
c$$$
c$$$                     endif
c$$$
c$$$                     if (chem_flag.eq.1) then
c$$$
c$$$                        iotataylorvert(k,ll,lll) = iotataylor(k,ll,1)+
c$$$     &                       iotataylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) ),1)*( xi2vert(k,lll) - xi2BCb(k))
c$$$     &                       + iotataylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) )+1,1)*( xi1vert(k,lll) - xi1BCb(k) )
c$$$                        iota2taylorvert(k,ll,lll) = iota2taylor(k,ll,1)+
c$$$     &                       iotataylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) ),1)*( xi2vert(k,lll) - xi2BCb(k))
c$$$     &                       + iota2taylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) )+1,1)*( xi1vert(k,lll) - xi1BCb(k) )
c$$$
c$$$
c$$$                     endif
c$$$
c$$$                  else
c$$$
c$$$                  endif
c$$$
c$$$               enddo
c$$$               
c$$$            enddo
c$$$            
c$$$         elseif (dofs(k).eq.1) then
c$$$
c$$$            cycle
c$$$
c$$$         endif
c$$$
c$$$      enddo
c$$$
c$$$C.....Compute alphas for each variable in each order derivitive
c$$$
c$$$
c$$$      alphaZE0 = 0.D0
c$$$      alphaQX0 = 0.D0
c$$$      alphaQY0 = 0.D0
c$$$
c$$$      if (tracer_flag.eq.1) then
c$$$
c$$$         alphaiota0 = 0.D0
c$$$
c$$$      endif
c$$$
c$$$      if (chem_flag.eq.1) then
c$$$
c$$$         alphaiota0 = 0.D0
c$$$         alphaiota20 = 0.D0
c$$$
c$$$      endif
c$$$
c$$$      do k=1,mne
c$$$
c$$$         if (dofs(k).gt.1) then
c$$$
c$$$            do lll=1,3
c$$$               
c$$$               do ll=1,dofs(k)
c$$$
c$$$                  if (ZEmin(nm(k,lll),ll).ne.ZEmax(nm(k,lll),ll)) then
c$$$
c$$$                     if ( ( ZEtaylorvert(k,ll,lll).gt.ZEtaylor(k,ll,1) ).and.
c$$$     &                    ( abs(ZEtaylorvert(k,ll,lll)-ZEtaylor(k,ll,1)).gt.1.0E-15 ).and.
c$$$     &                    ( ZEmax(nm(k,lll),ll).ne.ZEtaylor(k,ll,1) ) ) then  
c$$$
c$$$                        alphaZE0(k,ll,lll) = min(1.D0,  ( ZEmax(nm(k,lll),ll)
c$$$     &                       - ZEtaylor(k,ll,1) )/ (ZEtaylorvert(k,ll,lll) - ZEtaylor(k,ll,1)))
c$$$
c$$$                                !adapted part
c$$$
c$$$                     elseif ( ( ZEtaylorvert(k,ll,lll).gt.ZEtaylor(k,ll,1) ).and.
c$$$     &                       ( abs(ZEtaylorvert(k,ll,lll)-ZEtaylor(k,ll,1)).gt.1.0E-15 ).and.
c$$$     &                       ( ZEmax(nm(k,lll),ll).eq.ZEtaylor(k,ll,1) ) ) then 
c$$$
c$$$                        alphaZE0(k,ll,lll) = min(fd, abs( ( ZEmax(nm(k,lll),ll)
c$$$     &                       - ZEmin(nm(k,lll),ll) )/(ZEtaylorvert(k,ll,lll) - ZEmax(nm(k,lll),ll)) ) ) 
c$$$                        
c$$$
c$$$                     elseif ( (ZEtaylorvert(k,ll,lll).lt.ZEtaylor(k,ll,1) )
c$$$     &                       .and.( abs(ZEtaylorvert(k,ll,lll)-ZEtaylor(k,ll,1)).gt.1.0E-15 ).and.
c$$$     &                       ( ZEmin(nm(k,lll),ll).ne.ZEtaylor(k,ll,1) ) ) then 
c$$$
c$$$                        alphaZE0(k,ll,lll) = min( 1.D0,( ZEmin(nm(k,lll),ll)
c$$$     &                       - ZEtaylor(k,ll,1) )/( ZEtaylorvert(k,ll,lll)-ZEtaylor(k,ll,1)))
c$$$
c$$$                                !adapted part
c$$$
c$$$                     elseif ( (ZEtaylorvert(k,ll,lll).lt.ZEtaylor(k,ll,1) )
c$$$     &                       .and.( abs(ZEtaylorvert(k,ll,lll)-ZEtaylor(k,ll,1)).gt.1.0E-15 ).and.
c$$$     &                       ( ZEmin(nm(k,lll),ll).eq.ZEtaylor(k,ll,1) ) ) then 
c$$$
c$$$                        alphaZE0(k,ll,lll) = min( fd,abs( ( ZEmin(nm(k,lll),ll)
c$$$     &                       - ZEmax(nm(k,lll),ll) )/( ZEtaylorvert(k,ll,lll)- ZEmin(nm(k,lll),ll))) )
c$$$                        
c$$$                        
c$$$                     elseif ( ( ZEtaylorvert(k,ll,lll).eq.ZEtaylor(k,ll,1) ).or.
c$$$     &                       ( abs(ZEtaylorvert(k,ll,lll)-ZEtaylor(k,ll,1)).le.1.0E-15 ) ) then
c$$$
c$$$                        alphaZE0(k,ll,lll) = 1.D0
c$$$
c$$$                     endif
c$$$
c$$$                  else
c$$$
c$$$                     alphaZE0(k,ll,lll) = 1.D0
c$$$
c$$$                  endif
c$$$
c$$$                  if (tracer_flag.eq.1) then
c$$$
c$$$                     if (iotamin(nm(k,lll),ll).ne.iotamax(nm(k,lll),ll)) then
c$$$
c$$$                        if ( ( iotataylorvert(k,ll,lll).gt.iotataylor(k,ll,1) ).and.
c$$$     &                       ( abs(iotataylorvert(k,ll,lll)-iotataylor(k,ll,1)).gt.1.0E-15 ).and.
c$$$     &                       ( iotamax(nm(k,lll),ll).ne.iotataylor(k,ll,1) ) ) then  
c$$$
c$$$                           alphaiota0(k,ll,lll) = min(1.D0,  ( iotamax(nm(k,lll),ll)
c$$$     &                          - iotataylor(k,ll,1) )/ (iotataylorvert(k,ll,lll) - iotataylor(k,ll,1)))
c$$$  
c$$$                        elseif ( ( iotataylorvert(k,ll,lll).lt.iotataylor(k,ll,1) ).and.
c$$$     &                       ( abs(iotataylorvert(k,ll,lll)-iotataylor(k,ll,1)).gt.1.0E-15 ).and.
c$$$     &                       ( iotamax(nm(k,lll),ll).ne.iotataylor(k,ll,1) ) ) then  
c$$$
c$$$                           alphaiota0(k,ll,lll) = min(1.D0,  ( iotamin(nm(k,lll),ll)
c$$$     &                          - iotataylor(k,ll,1) )/ (iotataylorvert(k,ll,lll) - iotataylor(k,ll,1)))
c$$$
c$$$                                !adapted part
c$$$
c$$$                        elseif ( ( iotataylorvert(k,ll,lll).gt.iotataylor(k,ll,1) ).and.
c$$$     &                          ( abs(iotataylorvert(k,ll,lll)-iotataylor(k,ll,1)).gt.1.0E-15 ).and.
c$$$     &                          ( iotamin(nm(k,lll),ll).ne.iotataylor(k,ll,1) ).or.
c$$$     &                          ( iotamax(nm(k,lll),ll).eq.iotataylor(k,ll,1) ) ) then 
c$$$
c$$$                           alphaiota0(k,ll,lll) = min( fd,( iotamax(nm(k,lll),ll)
c$$$     &                          - iotamin(nm(k,lll),ll) )/(iotataylorvert(k,ll,lll) - iotamax(nm(k,lll),ll)) )
c$$$                           
c$$$                        elseif ( (iotataylorvert(k,ll,lll).lt.iotataylor(k,ll,1) )
c$$$     &                          .and.( abs(iotataylorvert(k,ll,lll)-iotataylor(k,ll,1)).gt.1.0E-15 ).and.
c$$$     &                          ( iotamin(nm(k,lll),ll).ne.iotataylor(k,ll,1) ) ) then 
c$$$
c$$$                           alphaiota0(k,ll,lll) = min( fd,( iotamax(nm(k,lll),ll) - iotamin(nm(k,lll),ll))
c$$$     &                          / (iotataylor(k,ll,1)- iotataylorvert(k,ll,lll) ) )
c$$$
c$$$                                !adapted part
c$$$            
c$$$                        elseif ( (iotataylorvert(k,ll,lll).lt.iotataylor(k,ll,1) )
c$$$     &                          .and.( abs(iotataylorvert(k,ll,lll)-iotataylor(k,ll,1)).gt.1.0E-15 ).and.
c$$$     &                          ( iotamin(nm(k,lll),ll).eq.iotataylor(k,ll,1) ) ) then 
c$$$
c$$$                           alphaiota0(k,ll,lll) =  min( fd, ( iotamin(nm(k,lll),ll)
c$$$     &                          - iotamax(nm(k,lll),ll) )/(iotataylorvert(k,ll,lll) - iotamin(nm(k,lll),ll)) )
c$$$
c$$$                        elseif ( ( iotataylorvert(k,ll,lll).eq.iotataylor(k,ll,1) ).or.
c$$$     &                          ( abs(iotataylorvert(k,ll,lll)-iotataylor(k,ll,1)).le.1.0E-15 ).or.
c$$$     &                           iotamin(nm(k,lll),ll).ne.iotamax(nm(k,lll),ll) ) then
c$$$
c$$$                           alphaiota0(k,ll,lll) = 1.D0
c$$$
c$$$                        endif
c$$$
c$$$                     else
c$$$
c$$$                        alphaiota0(k,ll,lll) = 1.D0
c$$$
c$$$                     endif
c$$$
c$$$                     if (alphaiota0(k,ll,lll).gt.1.d0.or.alphaiota0(k,ll,lll).lt.0.d0) then
c$$$
c$$$                        print*,ll
c$$$
c$$$                     endif
c$$$
c$$$                  endif
c$$$                  
c$$$                  if (chem_flag.eq.1) then
c$$$
c$$$                     if (iotamin(nm(k,lll),ll).ne.iotamax(nm(k,lll),ll)) then
c$$$
c$$$                        if ( ( iotataylorvert(k,ll,lll).gt.iotataylor(k,ll,1) ).and.
c$$$     &                       ( abs(iotataylorvert(k,ll,lll)-iotataylor(k,ll,1)).gt.1.0E-15 ).and.
c$$$     &                       ( iotamax(nm(k,lll),ll).ne.iotataylor(k,ll,1) ) ) then  
c$$$
c$$$                           alphaiota0(k,ll,lll) = min(1.D0,  ( iotamax(nm(k,lll),ll)
c$$$     &                          - iotataylor(k,ll,1) )/ (iotataylorvert(k,ll,lll) - iotataylor(k,ll,1)))
c$$$
c$$$                                !adapted part
c$$$
c$$$                        elseif ( ( iotataylorvert(k,ll,lll).gt.iotataylor(k,ll,1) ).and.
c$$$     &                          ( abs(iotataylorvert(k,ll,lll)-iotataylor(k,ll,1)).gt.1.0E-15 ).and.
c$$$     &                          ( iotamax(nm(k,lll),ll).eq.iotataylor(k,ll,1) ) ) then 
c$$$
c$$$                           alphaiota0(k,ll,lll) = min(fd, abs( ( iotamax(nm(k,lll),ll)
c$$$     &                          - iotamin(nm(k,lll),ll) )/(iotataylorvert(k,ll,lll) - iotamax(nm(k,lll),ll)) ) ) 
c$$$                           
c$$$
c$$$                        elseif ( (iotataylorvert(k,ll,lll).lt.iotataylor(k,ll,1) )
c$$$     &                          .and.( abs(iotataylorvert(k,ll,lll)-iotataylor(k,ll,1)).gt.1.0E-15 ).and.
c$$$     &                          ( iotamin(nm(k,lll),ll).ne.iotataylor(k,ll,1) ) ) then 
c$$$
c$$$                           alphaiota0(k,ll,lll) = min( 1.D0,( iotamin(nm(k,lll),ll)
c$$$     &                          - iotataylor(k,ll,1) )/( iotataylorvert(k,ll,lll)-iotataylor(k,ll,1)))
c$$$
c$$$                                !adapted part
c$$$
c$$$                        elseif ( (iotataylorvert(k,ll,lll).lt.iotataylor(k,ll,1) )
c$$$     &                          .and.( abs(iotataylorvert(k,ll,lll)-iotataylor(k,ll,1)).gt.1.0E-15 ).and.
c$$$     &                          ( iotamin(nm(k,lll),ll).eq.iotataylor(k,ll,1) ) ) then 
c$$$
c$$$                           alphaiota0(k,ll,lll) = min( fd,abs( ( iotamin(nm(k,lll),ll)
c$$$     &                          - iotamax(nm(k,lll),ll) )/( iotataylorvert(k,ll,lll)- iotamin(nm(k,lll),ll))) )
c$$$                           
c$$$                           
c$$$                        elseif ( ( iotataylorvert(k,ll,lll).eq.iotataylor(k,ll,1) ).or.
c$$$     &                          ( abs(iotataylorvert(k,ll,lll)-iotataylor(k,ll,1)).le.1.0E-15 ) ) then
c$$$
c$$$                           alphaiota0(k,ll,lll) = 1.D0
c$$$
c$$$                        endif
c$$$
c$$$                     else
c$$$
c$$$                        alphaiota0(k,ll,lll) = 1.D0
c$$$
c$$$                     endif
c$$$
c$$$                     if (iota2min(nm(k,lll),ll).ne.iota2max(nm(k,lll),ll)) then
c$$$
c$$$                        if ( ( iota2taylorvert(k,ll,lll).gt.iota2taylor(k,ll,1) ).and.
c$$$     &                       ( abs(iota2taylorvert(k,ll,lll)-iota2taylor(k,ll,1)).gt.1.0E-15 ).and.
c$$$     &                       ( iota2max(nm(k,lll),ll).ne.iota2taylor(k,ll,1) ) ) then  
c$$$
c$$$                           alphaiota20(k,ll,lll) = min(1.D0,  ( iota2max(nm(k,lll),ll)
c$$$     &                          - iota2taylor(k,ll,1) )/ (iota2taylorvert(k,ll,lll) - iota2taylor(k,ll,1)))
c$$$
c$$$                                !adapted part
c$$$
c$$$                        elseif ( ( iota2taylorvert(k,ll,lll).gt.iota2taylor(k,ll,1) ).and.
c$$$     &                          ( abs(iota2taylorvert(k,ll,lll)-iota2taylor(k,ll,1)).gt.1.0E-15 ).and.
c$$$     &                          ( iota2max(nm(k,lll),ll).eq.iota2taylor(k,ll,1) ) ) then 
c$$$
c$$$                           alphaiota20(k,ll,lll) = min(fd, abs( ( iota2max(nm(k,lll),ll)
c$$$     &                          - iota2min(nm(k,lll),ll) )/(iota2taylorvert(k,ll,lll) - iota2max(nm(k,lll),ll)) ) ) 
c$$$                           
c$$$
c$$$                        elseif ( (iota2taylorvert(k,ll,lll).lt.iota2taylor(k,ll,1) )
c$$$     &                          .and.( abs(iota2taylorvert(k,ll,lll)-iota2taylor(k,ll,1)).gt.1.0E-15 ).and.
c$$$     &                          ( iota2min(nm(k,lll),ll).ne.iota2taylor(k,ll,1) ) ) then 
c$$$
c$$$                           alphaiota20(k,ll,lll) = min( 1.D0,( iota2min(nm(k,lll),ll)
c$$$     &                          - iota2taylor(k,ll,1) )/( iota2taylorvert(k,ll,lll)-iota2taylor(k,ll,1)))
c$$$
c$$$                                !adapted part
c$$$
c$$$                        elseif ( (iota2taylorvert(k,ll,lll).lt.iota2taylor(k,ll,1) )
c$$$     &                          .and.( abs(iota2taylorvert(k,ll,lll)-iota2taylor(k,ll,1)).gt.1.0E-15 ).and.
c$$$     &                          ( iota2min(nm(k,lll),ll).eq.iota2taylor(k,ll,1) ) ) then 
c$$$
c$$$                           alphaiota20(k,ll,lll) = min( fd,abs( ( iota2min(nm(k,lll),ll)
c$$$     &                          - iota2max(nm(k,lll),ll) )/( iota2taylorvert(k,ll,lll)- iota2min(nm(k,lll),ll))) )
c$$$                           
c$$$                           
c$$$                        elseif ( ( iota2taylorvert(k,ll,lll).eq.iota2taylor(k,ll,1) ).or.
c$$$     &                          ( abs(iota2taylorvert(k,ll,lll)-iota2taylor(k,ll,1)).le.1.0E-15 ) ) then
c$$$
c$$$                           alphaiota20(k,ll,lll) = 1.D0
c$$$
c$$$                        endif
c$$$
c$$$                     else
c$$$
c$$$                        alphaiota20(k,ll,lll) = 1.D0
c$$$
c$$$                     endif
c$$$
c$$$                  endif
c$$$
c$$$                  
c$$$
c$$$                  if (QXmin(nm(k,lll),ll).ne.QXmax(nm(k,lll),ll)) then
c$$$
c$$$                     if ( ( QXtaylorvert(k,ll,lll).gt.QXtaylor(k,ll,1) ).and.
c$$$     &                    ( abs(QXtaylorvert(k,ll,lll)-QXtaylor(k,ll,1)).gt.1.0E-15 ).and.
c$$$     &                    ( QXmax(nm(k,lll),ll).ne.QXtaylor(k,ll,1) ) ) then  
c$$$
c$$$                        alphaQX0(k,ll,lll) = min(1.D0,  ( QXmax(nm(k,lll),ll)
c$$$     &                       - QXtaylor(k,ll,1) )/ (QXtaylorvert(k,ll,lll) - QXtaylor(k,ll,1)))
c$$$
c$$$                                !adapted part
c$$$
c$$$                     elseif ( ( QXtaylorvert(k,ll,lll).gt.QXtaylor(k,ll,1) ).and.
c$$$     &                       ( abs(QXtaylorvert(k,ll,lll)-QXtaylor(k,ll,1)).gt.1.0E-15 ).and.
c$$$     &                       ( QXmax(nm(k,lll),ll).eq.QXtaylor(k,ll,1) ) ) then 
c$$$
c$$$                        alphaQX0(k,ll,lll) = min(fd, abs( ( QXmax(nm(k,lll),ll)
c$$$     &                       - QXmin(nm(k,lll),ll) )/(QXtaylorvert(k,ll,lll) - QXmax(nm(k,lll),ll)) ) ) 
c$$$                        
c$$$
c$$$                     elseif ( (QXtaylorvert(k,ll,lll).lt.QXtaylor(k,ll,1) )
c$$$     &                       .and.( abs(QXtaylorvert(k,ll,lll)-QXtaylor(k,ll,1)).gt.1.0E-15 ).and.
c$$$     &                       ( QXmin(nm(k,lll),ll).ne.QXtaylor(k,ll,1) ) ) then 
c$$$
c$$$                        alphaQX0(k,ll,lll) = min( 1.D0,( QXmin(nm(k,lll),ll)
c$$$     &                       - QXtaylor(k,ll,1) )/( QXtaylorvert(k,ll,lll)-QXtaylor(k,ll,1)))
c$$$
c$$$                                !adapted part
c$$$
c$$$                     elseif ( (QXtaylorvert(k,ll,lll).lt.QXtaylor(k,ll,1) )
c$$$     &                       .and.( abs(QXtaylorvert(k,ll,lll)-QXtaylor(k,ll,1)).gt.1.0E-15 ).and.
c$$$     &                       ( QXmin(nm(k,lll),ll).eq.QXtaylor(k,ll,1) ) ) then 
c$$$
c$$$                        alphaQX0(k,ll,lll) = min( fd,abs( ( QXmin(nm(k,lll),ll)
c$$$     &                       - QXmax(nm(k,lll),ll) )/( QXtaylorvert(k,ll,lll)- QXmin(nm(k,lll),ll))) )
c$$$                        
c$$$                        
c$$$                     elseif ( ( QXtaylorvert(k,ll,lll).eq.QXtaylor(k,ll,1) ).or.
c$$$     &                       ( abs(QXtaylorvert(k,ll,lll)-QXtaylor(k,ll,1)).le.1.0E-15 ) ) then
c$$$
c$$$                        alphaQX0(k,ll,lll) = 1.D0
c$$$
c$$$                     endif
c$$$
c$$$                  else
c$$$
c$$$                     alphaQX0(k,ll,lll) = 1.D0
c$$$
c$$$                  endif
c$$$
c$$$
c$$$                  if (QYmin(nm(k,lll),ll).ne.QYmax(nm(k,lll),ll)) then
c$$$
c$$$                     if ( ( QYtaylorvert(k,ll,lll).gt.QYtaylor(k,ll,1) ).and.
c$$$     &                    ( abs(QYtaylorvert(k,ll,lll)-QYtaylor(k,ll,1)).gt.1.0E-15 ).and.
c$$$     &                    ( QYmax(nm(k,lll),ll).ne.QYtaylor(k,ll,1) ) ) then  
c$$$
c$$$                        alphaQY0(k,ll,lll) = min(1.D0,  ( QYmax(nm(k,lll),ll)
c$$$     &                       - QYtaylor(k,ll,1) )/ (QYtaylorvert(k,ll,lll) - QYtaylor(k,ll,1)))
c$$$
c$$$                                !adapted part
c$$$
c$$$                     elseif ( ( QYtaylorvert(k,ll,lll).gt.QYtaylor(k,ll,1) ).and.
c$$$     &                       ( abs(QYtaylorvert(k,ll,lll)-QYtaylor(k,ll,1)).gt.1.0E-15 ).and.
c$$$     &                       ( QYmax(nm(k,lll),ll).eq.QYtaylor(k,ll,1) ) ) then 
c$$$
c$$$                        alphaQY0(k,ll,lll) = min(fd, abs( ( QYmax(nm(k,lll),ll)
c$$$     &                       - QYmin(nm(k,lll),ll) )/(QYtaylorvert(k,ll,lll) - QYmax(nm(k,lll),ll)) ) ) 
c$$$                        
c$$$
c$$$                     elseif ( (QYtaylorvert(k,ll,lll).lt.QYtaylor(k,ll,1) )
c$$$     &                       .and.( abs(QYtaylorvert(k,ll,lll)-QYtaylor(k,ll,1)).gt.1.0E-15 ).and.
c$$$     &                       ( QYmin(nm(k,lll),ll).ne.QYtaylor(k,ll,1) ) ) then 
c$$$
c$$$                        alphaQY0(k,ll,lll) = min( 1.D0,( QYmin(nm(k,lll),ll)
c$$$     &                       - QYtaylor(k,ll,1) )/( QYtaylorvert(k,ll,lll)-QYtaylor(k,ll,1)))
c$$$
c$$$                                !adapted part
c$$$
c$$$                     elseif ( (QYtaylorvert(k,ll,lll).lt.QYtaylor(k,ll,1) )
c$$$     &                       .and.( abs(QYtaylorvert(k,ll,lll)-QYtaylor(k,ll,1)).gt.1.0E-15 ).and.
c$$$     &                       ( QYmin(nm(k,lll),ll).eq.QYtaylor(k,ll,1) ) ) then 
c$$$
c$$$                        alphaQY0(k,ll,lll) = min( fd,abs( ( QYmin(nm(k,lll),ll)
c$$$     &                       - QYmax(nm(k,lll),ll) )/( QYtaylorvert(k,ll,lll)- QYmin(nm(k,lll),ll))) )
c$$$                        
c$$$                        
c$$$                     elseif ( ( QYtaylorvert(k,ll,lll).eq.QYtaylor(k,ll,1) ).or.
c$$$     &                       ( abs(QYtaylorvert(k,ll,lll)-QYtaylor(k,ll,1)).le.1.0E-15 ) ) then
c$$$
c$$$                        alphaQY0(k,ll,lll) = 1.D0
c$$$
c$$$                     endif
c$$$
c$$$                  else
c$$$
c$$$                     alphaQY0(k,ll,lll) = 1.D0
c$$$
c$$$                  endif
c$$$
c$$$               enddo            !lll
c$$$
c$$$            enddo               !ll
c$$$
c$$$         elseif (dofs(k).eq.1) then
c$$$
c$$$            cycle
c$$$
c$$$         endif
c$$$
c$$$      enddo
c$$$
c$$$C.... Choose smallest (minimum) alpha for derivative in x or y
c$$$
c$$$      alphaZEm = 0.D0
c$$$      alphaQXm = 0.D0
c$$$      alphaQYm = 0.D0
c$$$
c$$$      if (tracer_flag.eq.1) then
c$$$
c$$$         alphaiotam = 0.D0
c$$$
c$$$      endif
c$$$
c$$$      if (chem_flag.eq.1) then
c$$$
c$$$         alphaiotam = 0.D0
c$$$         alphaiota2m = 0.D0
c$$$
c$$$      endif
c$$$
c$$$      do k=1,ne
c$$$
c$$$         if (dofs(k).gt.1) then
c$$$            
c$$$            do bb = 1,pdg_el(k)
c$$$
c$$$               if( (bb+1)*(bb+2)/2.le.dofs(k) ) then
c$$$
c$$$                  alphaZEm(k,(bb*(bb+1))/2 + 1) = 
c$$$     &                 minval( alphaZE(k,(bb*(bb-1))/2 +1:(bb*(bb+1)) / 2  ) )
c$$$                  alphaQXm(k,(bb*(bb+1))/2 + 1) = 
c$$$     &                 minval( alphaQX(k,(bb*(bb-1))/2 +1:(bb*(bb+1)) / 2  ) )
c$$$                  alphaQYm(k,(bb*(bb+1))/2 + 1) = 
c$$$     &                 minval( alphaQY(k,(bb*(bb-1))/2 +1:(bb*(bb+1)) / 2  ) )
c$$$
c$$$                  if (tracer_flag.eq.1) then
c$$$
c$$$                     alphaiotam(k,(bb*(bb+1))/2 + 1) = 
c$$$     &                    minval( alphaiota(k,( (bb*(bb-1))/2 +1) :(bb*(bb+1)) / 2 ) )
c$$$
c$$$                  endif
c$$$
c$$$                  if (chem_flag.eq.1) then
c$$$
c$$$                     alphaiotam(k,(bb*(bb+1))/2 + 1) = 
c$$$     &                    minval( alphaiota(k,(bb*(bb-1))/2 +1:(bb*(bb+1)) / 2  ) )
c$$$                     alphaiota2m(k,(bb*(bb+1))/2 + 1) = 
c$$$     &                    minval( alphaiota2(k,(bb*(bb-1))/2 +1:(bb*(bb+1)) / 2 ) )
c$$$
c$$$                  endif
c$$$
c$$$               endif
c$$$               
c$$$            enddo
c$$$
c$$$         elseif (dofs(k).eq.1) then
c$$$
c$$$            cycle
c$$$
c$$$         endif
c$$$
c$$$      enddo
c$$$
c$$$C.....Use max higher derivative values for final limiter value
c$$$
c$$$      alphaZE_max = 0.D0
c$$$      alphaQX_max = 0.D0
c$$$      alphaQY_max = 0.D0
c$$$
c$$$      if (tracer_flag.eq.1) then
c$$$
c$$$         alphaiota_max = 0.D0
c$$$
c$$$      endif
c$$$
c$$$      if (chem_flag.eq.1) then
c$$$
c$$$         alphaiota_max = 0.D0
c$$$         alphaiota2_max = 0.D0
c$$$
c$$$      endif
c$$$
c$$$      do k=1,ne
c$$$
c$$$         if (dofs(k).gt.1) then
c$$$
c$$$            do bb =1,pdg_el(k)
c$$$
c$$$               if( (bb+1)*(bb+2)/2.le.dofs(k)) then
c$$$
c$$$                  alphaZE_max(k,(bb*(bb+1))/2 + 1) = fd*maxval( 
c$$$     &                 alphaZEm(k,(bb*(bb+1))/2 + 1:dofs(k)) )
c$$$                  alphaQX_max(k,(bb*(bb+1))/2 + 1) = fd*maxval( 
c$$$     &                 alphaQXm(k,(bb*(bb+1))/2 + 1:dofs(k)) )
c$$$                  alphaQY_max(k,(bb*(bb+1))/2 + 1) = fd*maxval( 
c$$$     &                 alphaQYm(k,(bb*(bb+1))/2 + 1:dofs(k)) )
c$$$
c$$$                  if (tracer_flag.eq.1) then
c$$$
c$$$                     alphaiota_max(k,(bb*(bb+1))/2 + 1) = maxval( 
c$$$     &                    alphaiotam(k,( (bb*(bb+1))/2 + 1):dofs(k)) )
c$$$
c$$$                  endif
c$$$
c$$$                  if (chem_flag.eq.1) then
c$$$
c$$$                     alphaiota_max(k,(bb*(bb+1))/2 + 1) = fd*maxval( 
c$$$     &                    alphaiotam(k,(bb*(bb+1))/2 + 1:dofs(k)) )
c$$$                     alphaiota2_max(k,(bb*(bb+1))/2 + 1) = fd*maxval( 
c$$$     &                    alphaiota2m(k,(bb*(bb+1))/2 + 1:dofs(k)) )
c$$$
c$$$                  endif
c$$$
c$$$               endif
c$$$
c$$$            enddo
c$$$
c$$$         elseif (dofs(k).eq.1) then
c$$$
c$$$            cycle
c$$$
c$$$         endif
c$$$
c$$$      enddo
c$$$
c$$$C.....Limit on the Master element in the Taylor basis, via reconstruction 
c$$$C.....of unconstrained solutions with alpha constraints
c$$$
c$$$
c$$$      limitZE = 0.D0
c$$$      limitQX = 0.D0
c$$$      limitQY = 0.D0
c$$$
c$$$      lim_count_roll = 0.D0
c$$$
c$$$      if (tracer_flag.eq.1) then
c$$$
c$$$         limitiota = 0.D0
c$$$
c$$$      endif
c$$$
c$$$      if (chem_flag.eq.1) then
c$$$
c$$$         limitiota = 0.D0
c$$$         limitiota2 = 0.D0
c$$$
c$$$      endif
c$$$
c$$$      do k=1,ne
c$$$
c$$$         lim_count = 0
c$$$
c$$$         if (dofs(k).gt.1) then
c$$$
c$$$            do ll=1,dofs(k)
c$$$
c$$$               if ( ll.eq.1 ) then
c$$$
c$$$                  limitZE(k,ll) = ZEtaylor(k,ll,1)
c$$$                  limitQX(k,ll) = QXtaylor(k,ll,1)
c$$$                  limitQY(k,ll) = QYtaylor(k,ll,1) 
c$$$
c$$$                  if (tracer_flag.eq.1) then
c$$$
c$$$                     limitiota(k,ll) = iotataylor(k,ll,1)
c$$$
c$$$                  endif
c$$$
c$$$                  if (chem_flag.eq.1) then
c$$$
c$$$                     limitiota(k,ll) = iotataylor(k,ll,1)
c$$$                     limitiota2(k,ll) = iota2taylor(k,ll,1)
c$$$
c$$$                  endif
c$$$
c$$$               elseif ( ll.ge.2 ) then
c$$$                  
c$$$                  do bb=1,pdg_el(k)
c$$$
c$$$                     if ( ll.le.( ( (bb+1)*(bb+2)) / 2.D0 ).and.(ll.gt.
c$$$     &                    (bb*(bb+1)/2.D0) ) ) then
c$$$
c$$$                        limitZE(k,ll) = alphaZE_max(k,(bb*(bb+1))/2 + 1) 
c$$$     &                       * ZEtaylor(k,ll,1)
c$$$                        limitQX(k,ll) = alphaQX_max(k,(bb*(bb+1))/2 + 1) 
c$$$     &                       * ZEtaylor(k,ll,1)
c$$$                        limitQY(k,ll) = alphaQY_max(k,(bb*(bb+1))/2 + 1) 
c$$$     &                       * ZEtaylor(k,ll,1)
c$$$
c$$$
c$$$
c$$$                        if (tracer_flag.eq.1) then
c$$$
c$$$                           limitiota(k,ll) = alphaiota_max(k,(bb*(bb+1))/2 + 1) 
c$$$     &                          * iotataylor(k,ll,1)
c$$$
c$$$                        endif
c$$$
c$$$                        if (chem_flag.eq.1) then
c$$$
c$$$                           limitiota(k,ll) = alphaiota_max(k,(bb*(bb+1))/2 + 1) 
c$$$     &                          * iotataylor(k,ll,1)
c$$$                           limitiota2(k,ll) = alphaiota2_max(k,(bb*(bb+1))/2 + 1) 
c$$$     &                          * iota2taylor(k,ll,1)
c$$$
c$$$                        endif
c$$$
c$$$
c$$$                                ! Make a counter to track limiting
c$$$
c$$$                        if ( ( alphaZE_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
c$$$     &                       alphaQX_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
c$$$     &                       alphaQY_max(k,(bb*(bb+1))/2 + 1).lt.1.D0 ).and.
c$$$     &                       chem_flag.ne.1.and.tracer_flag.ne.1 ) then
c$$$                           
c$$$                           lim_count = 0  
c$$$
c$$$                        elseif ( ( alphaiota_max(k,(bb*(bb+1))/2 + 1).lt.1.D0 ).and.
c$$$     &                          tracer_flag.eq.1 ) then
c$$$                           
c$$$                           !print*,alphaiota_max(k,(bb*(bb+1))/2 + 1)
c$$$
c$$$                           lim_count = 1   
c$$$                           
c$$$                        elseif ( (alphaZE_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
c$$$     &                          alphaQX_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
c$$$     &                          alphaQY_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
c$$$     &                          alphaiota_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
c$$$     &                          alphaiota2_max(k,(bb*(bb+1))/2 + 1).lt.1.D0 ).and.
c$$$     &                          chem_flag.eq.1 ) then
c$$$
c$$$                           lim_count = 0 
c$$$                           
c$$$                        endif
c$$$                        
c$$$                     endif
c$$$                     
c$$$                  enddo
c$$$                  
c$$$               endif
c$$$
c$$$            enddo
c$$$
c$$$         elseif (dofs(k).eq.1) then
c$$$
c$$$            cycle
c$$$
c$$$         endif
c$$$
c$$$         lim_count_roll = lim_count_roll + lim_count
c$$$
c$$$      enddo
c$$$
c$$$C.....Transform back to the Dubiner basis (multiply by NmatrixInv),
c$$$
c$$$      ZEconst =  0.D0
c$$$      QXconst =  0.D0
c$$$      QYconst =  0.D0
c$$$
c$$$      if (tracer_flag.eq.1) then
c$$$
c$$$         iotaconst =  0.D0
c$$$
c$$$      endif
c$$$
c$$$      if (chem_flag.eq.1) then
c$$$
c$$$         iotaconst =  0.D0
c$$$         iota2const =  0.D0
c$$$
c$$$      endif
c$$$
c$$$      do k=1,ne
c$$$
c$$$         if (dofs(k).gt.1) then
c$$$                                !do lll=1,3
c$$$
c$$$            do ll=1,dofs(k)
c$$$
c$$$               do ss=1,dofs(k)
c$$$
c$$$                  ZEconst(k,ll) = ZEconst(k,ll) + NmatrixInv(k,ll,ss,dofs(k)) 
c$$$     &                 * limitZE(k,ss)
c$$$                  QXconst(k,ll) = QXconst(k,ll) + NmatrixInv(k,ll,ss,dofs(k)) 
c$$$     &                 * limitQX(k,ss)
c$$$                  QYconst(k,ll) = QYconst(k,ll) + NmatrixInv(k,ll,ss,dofs(k)) 
c$$$     &                 * limitQY(k,ss)
c$$$
c$$$                  if (tracer_flag.eq.1) then
c$$$
c$$$                     iotaconst(k,ll) = iotaconst(k,ll) + 
c$$$     &                    NmatrixInv(k,ll,ss,dofs(k)) * limitiota(k,ss)
c$$$                     
c$$$                  endif
c$$$
c$$$                  if (chem_flag.eq.1) then
c$$$
c$$$                     iotaconst(k,ll) = iotaconst(k,ll) + 
c$$$     &                    NmatrixInv(k,ll,ss,dofs(k)) * limitiota(k,ss)
c$$$                     iota2const(k,ll) = iota2const(k,ll) + 
c$$$     &                    NmatrixInv(k,ll,ss,dofs(k)) * limitiota2(k,ss)
c$$$
c$$$                  endif
c$$$
c$$$
c$$$               enddo
c$$$
c$$$            enddo
c$$$
c$$$         elseif (dofs(k).eq.1) then
c$$$
c$$$            cycle
c$$$
c$$$         endif
c$$$
c$$$      enddo
c$$$
c$$$C.....Set limit values
c$$$
c$$$      do k=1,ne
c$$$
c$$$         if (dofs(k).gt.1) then
c$$$
c$$$            do ll = 1,dofs(k)
c$$$
c$$$               ZE(ll,k,irk+1) = ZEconst(k,ll)
c$$$               QX(ll,k,irk+1) = QXconst(k,ll)
c$$$               QY(ll,k,irk+1) = QYconst(k,ll)
c$$$
c$$$               if (tracer_flag.eq.1) then
c$$$
c$$$                  iota(ll,k,irk+1) = iotaconst(k,ll)
c$$$
c$$$               endif
c$$$
c$$$               if (chem_flag.eq.1) then
c$$$
c$$$                  iota(ll,k,irk+1) = iotaconst(k,ll)
c$$$                  iota2(ll,k,irk+1) = iota2const(k,ll)
c$$$
c$$$               endif
c$$$
c$$$            enddo
c$$$
c$$$         elseif (dofs(k).eq.1) then
c$$$
c$$$            cycle
c$$$
c$$$         endif
c$$$
c$$$      enddo
c$$$      
c$$$      return
c$$$      end subroutine




C.....Subroutine to find the inverse of a square matrix by Guass-Jordan elimination
C.....cem

      subroutine Inv2(matrix, inverse, n)
      Use sizes, only : sz

      implicit none
      integer n
      real(sz), dimension(n,n), intent(in) :: matrix
      real(sz), dimension(n,n), intent(inout) :: inverse
      
      integer :: i, j, k, l
      real(sz) :: m
      real(sz), dimension(n,2*n) :: augmatrix !augmented matrix
      
                                !Augment input matrix with an identity matrix

      do i = 1, n

         do j = 1, 2*n

            if ( j.le.n ) then

               augmatrix(i,j) = matrix(i,j)

            else if ((i+n) == j) then

               augmatrix(i,j) = 1

            else

               augmatrix(i,j) = 0

            endif

         enddo

      enddo
   
                                !Reduce augmented matrix to upper traingular form

      do k =1, n-1

         if (augmatrix(k,k) == 0) then


            do i = k+1, n

               if (augmatrix(i,k) /= 0) then

                  do j = 1,2*n

                     augmatrix(k,j) = augmatrix(k,j)+augmatrix(i,j)

                  enddo

               endif


            enddo

         endif

         do j = k+1, n			

            m = augmatrix(j,k)/augmatrix(k,k)

            do i = k, 2*n

               augmatrix(j,i) = augmatrix(j,i) - m*augmatrix(k,i)

            enddo

         enddo

      enddo 
      
                                !Test for invertibility

      do i = 1, n

         if (augmatrix(i,i) == 0) then

            inverse = 0

            return

         endif

      enddo
      
                                !Make diagonal elements as 1

      do i = 1 , n

         m = augmatrix(i,i)

         do j = i, (2 * n)				

            augmatrix(i,j) = (augmatrix(i,j) / m)

         enddo

      enddo
      
                                !Reduced right side half of augmented matrix to identity matrix

      do k = n-1, 1, -1

         do i =1, k

            m = augmatrix(i,k+1)

            do j = k, (2*n)

               augmatrix(i,j) = augmatrix(i,j) -augmatrix(k+1,j) * m

            enddo

         enddo

      enddo				
      
                                !Compute answer

      do i =1, n

         do j = 1, n

            inverse(i,j) = augmatrix(i,j+n)

         enddo

      enddo

      end subroutine Inv2

#endif
