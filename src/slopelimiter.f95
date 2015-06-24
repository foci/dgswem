!***********************************************************************
!
!     SUBROUTINE SLOPELIMITER()
!
!     This subroutine selects either SLOPELIMITER1() or SLOPELIMITER2()
!     according to SLOPEFLAG. SLOPEFLAG is specified in fort.dg.
!   
!     All routines rewritten for p_adaptive multicomponent version
!     Slopelimiters 2 and 3 are not compatible. 1,4,5,6,7,8,9,10 are.
!     -- cem, 2011
!
!***********************************************************************

      SUBROUTINE SLOPELIMITER(s)

      USE SIZES
      USE DG, ONLY : SLOPEFLAG
      
      IMPLICIT NONE
      
      type (sizes_type) :: s

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
        CALL SLOPELIMITER5(s)
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
        CALL SLOPELIMITER5(s)
      ENDIF
#endif

#ifdef STBLZR
        CALL SLOPELIMITER5(s)
#endif
    

      RETURN
      END SUBROUTINE

#ifdef SLOPEALL

!***********************************************************************
!     
!     SUBROUTINE SLOPELIMITER1()
!     
!     Written by Ethan Kubatko
!     
!     08-Feb-2008
!     - Modified to apply this slopelimiter on domain 
!     boundaries. S.B.
!     28-Jun-2010 Modified for transport and chemistry and p enrichment
!     -cem
!     
!***********************************************************************

      SUBROUTINE SLOPELIMITER1()

!.....Use appropriate modules

      USE GLOBAL
      USE DG

      IMPLICIT NONE

!.....Declare local variables

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

!.....Save the original values  02/28/2007 sb

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

!.....Retrieve the barycenter coordinates of the element
            
            X1 = XBC(L)
            Y1 = YBC(L)

!.....Compute the x and y derivatives of the P1 part of ZE

            DZEDR = ZE(3,L,NRK+2)
            DZEDS = 3.D0/2.D0*ZE(2,L,NRK+2) + 1.D0/2.D0*ZE(3,L,NRK+2)

            DZEDX(4) = DZEDR*DRDX(L) + DZEDS*DSDX(L)
            DZEDY(4) = DZEDR*DRDY(L) + DZEDS*DSDY(L)
            
            GRADZE(4) = SQRT(DZEDX(4)**2.D0 + DZEDY(4)**2.D0)
            
!.....Compute the x and y derivatives of the P1 part of QX
            
            DQXDR = QX(3,L,NRK+2)
            DQXDS = 3.D0/2.D0*QX(2,L,NRK+2) + 1.D0/2.D0*QX(3,L,NRK+2)

            DQXDX(4) = DQXDR*DRDX(L) + DQXDS*DSDX(L)
            DQXDY(4) = DQXDR*DRDY(L) + DQXDS*DSDY(L)
            
            GRADQX(4) = SQRT(DQXDX(4)**2.D0 + DQXDY(4)**2.D0)
            
!.....Compute the x and y derivatives of the P1 part of QY
            
            DQYDR = QY(3,L,NRK+2)
            DQYDS = 3.D0/2.D0*QY(2,L,NRK+2) + 1.D0/2.D0*QY(3,L,NRK+2)

            DQYDX(4) = DQYDR*DRDX(L) + DQYDS*DSDX(L)
            DQYDY(4) = DQYDR*DRDY(L) + DQYDS*DSDY(L)
            
            GRADQY(4) = SQRT(DQYDX(4)**2.D0 + DQYDY(4)**2.D0)


!.....Compute the x and y derivatives of the P1 part of iota & iota2

#ifdef TRACE     
            DiotaDR = iota(3,L,NRK+2)
            DiotaDS = 3.D0/2.D0*iota(2,L,NRK+2) + 1.D0/2.D0*iota(3,L,NRK+2)

            DiotaDX(4) = DiotaDR*DRDX(L) + DiotaDS*DSDX(L)
            DiotaDY(4) = DiotaDR*DRDY(L) + DiotaDS*DSDY(L)
            
            GRADiota(4) = SQRT(DiotaDX(4)**2.D0 + DiotaDY(4)**2.D0)
#endif

#ifdef CHEM
            DiotaDR = iota(3,L,NRK+2)
            DiotaDS = 3.D0/2.D0*iota(2,L,NRK+2) + 1.D0/2.D0*iota(3,L,NRK+2)

            Diota2DR = iota2(3,L,NRK+2)
            Diota2DS = 3.D0/2.D0*iota2(2,L,NRK+2) + 1.D0/2.D0*iota2(3,L,NRK+2)

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

!.......Special treatment for elements that do not have three neighboring elements
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
               
!.....Compute limiting slopes for the surface elevation

               SL1 = ZE(1,L,NRK+2)*(Y3 - Y2) + ZE(1,EL2,NRK+2)*(Y1 - Y3)
     &              + ZE(1,EL3,NRK+2)*(Y2 - Y1)
               
               SL2 = ZE(1,L,NRK+2)*(X2 - X3) + ZE(1,EL2,NRK+2)*(X3 - X1)
     &              + ZE(1,EL3,NRK+2)*(X1 - X2)
               
               DZEDX(I) = -SL1/SL3(I,L)
               
               DZEDY(I) = -SL2/SL3(I,L)
               
               GRADZE(I) = SQRT(DZEDX(I)**2.D0 + DZEDY(I)**2.D0)
               
!.....Compute limiting slopes for the X flow

               SL1 = QX(1,L,NRK+2)*(Y3 - Y2) + QX(1,EL2,NRK+2)*(Y1 - Y3)
     &              + QX(1,EL3,NRK+2)*(Y2 - Y1)

               SL2 = QX(1,L,NRK+2)*(X2 - X3) + QX(1,EL2,NRK+2)*(X3 - X1)
     &              + QX(1,EL3,NRK+2)*(X1 - X2)

               DQXDX(I) = -SL1/SL3(I,L)

               DQXDY(I) = -SL2/SL3(I,L)
               
               GRADQX(I) = SQRT(DQXDX(I)**2.D0 + DQXDY(I)**2.D0)
               
!.....Compute limiting slopes for the Y flow

               SL1 = QY(1,L,NRK+2)*(Y3 - Y2) + QY(1,EL2,NRK+2)*(Y1 - Y3)
     &              + QY(1,EL3,NRK+2)*(Y2 - Y1)

               SL2 = QY(1,L,NRK+2)*(X2 - X3) + QY(1,EL2,NRK+2)*(X3 - X1)
     &              + QY(1,EL3,NRK+2)*(X1 - X2)

               DQYDX(I) = -SL1/SL3(I,L)

               DQYDY(I) = -SL2/SL3(I,L)
               
               GRADQY(I) = SQRT(DQYDX(I)**2.D0 + DQYDY(I)**2.D0)
               

!.....Compute limiting slopes for the concentrations

#ifdef TRACE
               SL1 = iota(1,L,NRK+2)*(Y3 - Y2) + iota(1,EL2,NRK+2)*(Y1 - Y3)
     &              + iota(1,EL3,NRK+2)*(Y2 - Y1)
               
               SL2 = iota(1,L,NRK+2)*(X2 - X3) + iota(1,EL2,NRK+2)*(X3 - X1)
     &              + iota(1,EL3,NRK+2)*(X1 - X2)
               
               DiotaDX(I) = -SL1/SL3(I,L)
               
               DiotaDY(I) = -SL2/SL3(I,L)
               
               GRADiota(I) = SQRT(DiotaDX(I)**2.D0 + DiotaDY(I)**2.D0)
#endif

#ifdef CHEM
               SL1 = iota(1,L,NRK+2)*(Y3 - Y2) + iota(1,EL2,NRK+2)*(Y1 - Y3)
     &              + iota(1,EL3,NRK+2)*(Y2 - Y1)
               
               SL2 = iota(1,L,NRK+2)*(X2 - X3) + iota(1,EL2,NRK+2)*(X3 - X1)
     &              + iota(1,EL3,NRK+2)*(X1 - X2)
               
               DiotaDX(I) = -SL1/SL3(I,L)
               
               DiotaDY(I) = -SL2/SL3(I,L)
               
               GRADiota(I) = SQRT(DiotaDX(I)**2.D0 + DiotaDY(I)**2.D0)

               SL1 = iota2(1,L,NRK+2)*(Y3 - Y2) + iota2(1,EL2,NRK+2)*(Y1 - Y3)
     &              + iota2(1,EL3,NRK+2)*(Y2 - Y1)
               
               SL2 = iota2(1,L,NRK+2)*(X2 - X3) + iota2(1,EL2,NRK+2)*(X3 - X1)
     &              + iota2(1,EL3,NRK+2)*(X1 - X2)
               
               Diota2DX(I) = -SL1/SL3(I,L)
               
               Diota2DY(I) = -SL2/SL3(I,L)
               
               GRADiota2(I) = SQRT(Diota2DX(I)**2.D0 + Diota2DY(I)**2.D0)
#endif

            ENDDO
            
!.....Compute the (possibly) limited gradient

            ZELIM = MIN( GRADZE(1), GRADZE(2), GRADZE(3), GRADZE(4) )
            QXLIM = MIN( GRADQX(1), GRADQX(2), GRADQX(3), GRADQX(4) )
            QYLIM = MIN( GRADQY(1), GRADQY(2), GRADQY(3), GRADQY(4) )

#ifdef TRACE
            iotaLIM = MIN( GRADiota(1), GRADiota(2), GRADiota(3), GRADiota(4) )
#endif

#ifdef CHEM
            iotaLIM = MIN( GRADiota(1), GRADiota(2), GRADiota(3), GRADiota(4) )
            iota2LIM = MIN( GRADiota2(1), GRADiota2(2), GRADiota2(3), GRADiota2(4) )
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
     &           - 2.D0*DRDY(L)*ZELIMX - ZELIMX*DSDY(L) + DSDX(L)*ZELIMY)/DEN
            QX(2,L,IRK+1) = 1.D0/3.D0*(2.D0*QXLIMY*DRDX(L)
     &           - 2.D0*DRDY(L)*QXLIMX - QXLIMX*DSDY(L) + DSDX(L)*QXLIMY)/DEN
            QY(2,L,IRK+1) = 1.D0/3.D0*(2.D0*QYLIMY*DRDX(L)
     &           - 2.D0*DRDY(L)*QYLIMX - QYLIMX*DSDY(L) + DSDX(L)*QYLIMY)/DEN

            ZE(3,L,IRK+1) = (ZELIMX*DSDY(L) - DSDX(L)*ZELIMY)/DEN
            QX(3,L,IRK+1) = (QXLIMX*DSDY(L) - DSDX(L)*QXLIMY)/DEN
            QY(3,L,IRK+1) = (QYLIMX*DSDY(L) - DSDX(L)*QYLIMY)/DEN

#ifdef TRACE
            iota(2,L,IRK+1) = 1.D0/3.D0*(2.D0*iotaLIMY*DRDX(L)
     &           - 2.D0*DRDY(L)*iotaLIMX - iotaLIMX*DSDY(L) + DSDX(L)*iotaLIMY)/DEN

            iota(3,L,IRK+1) = (iotaLIMX*DSDY(L) - DSDX(L)*iotaLIMY)/DEN
#endif

#ifdef CHEM
            iota(2,L,IRK+1) = 1.D0/3.D0*(2.D0*iotaLIMY*DRDX(L)
     &           - 2.D0*DRDY(L)*iotaLIMX - iotaLIMX*DSDY(L) + DSDX(L)*iotaLIMY)/DEN
            iota2(2,L,IRK+1) = 1.D0/3.D0*(2.D0*iota2LIMY*DRDX(L)
     &           - 2.D0*DRDY(L)*iota2LIMX - iota2LIMX*DSDY(L) + DSDX(L)*iota2LIMY)/DEN

            iota(3,L,IRK+1) = (iotaLIMX*DSDY(L) - DSDX(L)*iotaLIMY)/DEN
            iota2(3,L,IRK+1) = (iota2LIMX*DSDY(L) - DSDX(L)*iota2LIMY)/DEN
#endif
            
 111        continue

         endif

 1000 CONTINUE
      RETURN
      END SUBROUTINE

      
!***********************************************************************
!     
!     SUBROUTINE SLOPELIMITER2()
!     
!     Written by Shintaro Bunya - 08 Feb 2008
!     
!     This subroutine is an implementation of the slope limiter used in
!     B. Cockburn and C-W. Shu, "The Runge-Kutta Discontinuous Galerkin
!     Method for Conservation Laws V, Multidimensional Systems,"
!     Journal of Computational Physics 141, 199-224 (1998).
!     
!***********************************************************************

      SUBROUTINE SLOPELIMITER2()

!.....Use appropriate modules

      USE SIZES, ONLY : SZ
      USE GLOBAL
      USE DG

      IMPLICIT NONE

!.....Declare local variables

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

!.....Save the original values  02/28/2007 sb
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

!.....Retrieve the barycenter info of this element
         XB(0) = XBC(L)
         YB(0) = YBC(L)
         ZB(0) = ZE(1,L,IRK+1)
         DB(0) = HB(1,L,1)
         QXB(0) = QX(1,L,IRK+1)
         QYB(0) = QY(1,L,IRK+1)
         
!.....Retrieve the nodal info

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

!.....Retrieve info at barycenters and edge midpoints
         DO I = 1,3
            EL_N = EL_NBORS(I,L)

!.......Edge midpoints
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
            
!.......Barycenters
            IF(EL_N == 0) THEN
!     If there is no neighboring element on this side (e.g. domain boundary),
!     locate a barycenter at the mirroring point and use the values at the 
!     baricenter of this element.
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

!.......Edge normal vector
            GED = NELED(I,L)
            NXX(I) = COSNX(GED)
            NYY(I) = SINNX(GED)
         ENDDO

!.....Compute alpha1 and alpha2      
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


!.....Set the variable vectors at the baricenter of element 0
         WB_O(1,0) = ZB(0)
         WB_O(2,0) = QXB(0)
         WB_O(3,0) = QYB(0)

!.....Start computing deltas for each edge
         DO I = 1,3

!.......Compute eigen values and vectors at the midpoint using the variables at the barycenter
            HTB = ZB(0)+DB(0)
            CALL HYDRO_EIGEN_VALUES(HTB, QXB(0)/HTB, QYB(0)/HTB,
     &           NXX(I), NYY(I), EIGVAL, RI, LE)

!.......Set variable vectors
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

!.......Transform original variables into the characteristic space
            DO J = 0,2
               WB_C(1,J)
     &              = LE(1,1)*WB_O(1,J) + LE(1,2)*WB_O(2,J) + LE(1,3)*WB_O(3,J)
               WB_C(2,J)
     &              = LE(2,1)*WB_O(1,J) + LE(2,2)*WB_O(2,J) + LE(2,3)*WB_O(3,J)
               WB_C(3,J)
     &              = LE(3,1)*WB_O(1,J) + LE(3,2)*WB_O(2,J) + LE(3,3)*WB_O(3,J)
            ENDDO
            WM_C(1)= LE(1,1)*WM_O(1) + LE(1,2)*WM_O(2) + LE(1,3)*WM_O(3)
            WM_C(2)= LE(2,1)*WM_O(1) + LE(2,2)*WM_O(2) + LE(2,3)*WM_O(3)
            WM_C(3)= LE(3,1)*WM_O(1) + LE(3,2)*WM_O(2) + LE(3,3)*WM_O(3)

!.......Compute W_TILDA
            W_TILDA(1) = WM_C(1) - WB_C(1,0)
            W_TILDA(2) = WM_C(2) - WB_C(2,0)
            W_TILDA(3) = WM_C(3) - WB_C(3,0)
            
!.......Compute DELTA_W
            DELTA_W(1)
     &           = ALPHA1(I)*(WB_C(1,1) - WB_C(1,0))
     &           + ALPHA2(I)*(WB_C(1,2) - WB_C(1,0))
            DELTA_W(2)
     &           = ALPHA1(I)*(WB_C(2,1) - WB_C(2,0))
     &           + ALPHA2(I)*(WB_C(2,2) - WB_C(2,0))
            DELTA_W(3)
     &           = ALPHA1(I)*(WB_C(3,1) - WB_C(3,0))
     &           + ALPHA2(I)*(WB_C(3,2) - WB_C(3,0))

!.......Apply the TVB modified minmod function
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

!.......Transform back to the original space
            DELTA_O(1,I) = RI(1,1)*DELTA_C(1,I) + RI(1,2)*DELTA_C(2,I)
     &           + RI(1,3)*DELTA_C(3,I)
            DELTA_O(2,I) = RI(2,1)*DELTA_C(1,I) + RI(2,2)*DELTA_C(2,I)
     &           + RI(2,3)*DELTA_C(3,I)
            DELTA_O(3,I) = RI(3,1)*DELTA_C(1,I) + RI(3,2)*DELTA_C(2,I)
     &           + RI(3,3)*DELTA_C(3,I)

         ENDDO                  ! Loop for edges ends here

!     Set the variable vector at the midpoint 1. This is used later to judge
!     whether the variables should be limited.
         WM_O(1) = ZM(1)
         WM_O(2) = QXM(1)
         WM_O(3) = QYM(1)

!     Find the possibly limited DELTA and variable vectors at the corner nodes
         DO I = 1,3             ! Loop for variable components
!     Find the possibly limited DELTA
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

!     Compute the possibly limited variable vectors at the corner nodes
            IF(( DELTA_HAT(I,1) < (WM_O(I) - WB_O(I,0) - ZERO) ).OR.
     &           ( DELTA_HAT(I,1) > (WM_O(I) - WB_O(I,0) + ZERO) )) THEN
!     Update the variables only if they need to be limited.
!     If the condition in the if statement is false, the following 
!     isn't applied and therefore the quradratic and higher order 
!     coefficients are preserved.

!     Computing limited variables at the corner nodes
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

!***********************************************************************
!     
!     SUBROUTINE SLOPELIMITER3()
!     
!     Written by Shintaro Bunya - 27 Feb 2008
!     
!     This subroutine is an extentioin of the slope limiter used in
!     B. Cockburn and C-W. Shu, "The Runge-Kutta Discontinuous Galerkin
!     Method for Conservation Laws V, Multidimensional Systems,"
!     Journal of Computational Physics 141, 199-224 (1998).
!     
!     An additional parameter SL3_MD is used in this slope limiter to
!     make sure that slopes are limited in elements whose water depth
!     is small. (A special treatment for wetting and drying.
!     
!***********************************************************************

      SUBROUTINE SLOPELIMITER3()

!.....Use appropriate modules

      USE SIZES, ONLY : SZ
      USE GLOBAL
      USE DG

      IMPLICIT NONE

!.....Declare local variables

      INTEGER L, NNBORS,i,j,kk

      DO L= 1,NE
!     
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
!     CYCLE
#else
!     STOP 'STILL UNDER CONSTRUCTION!!!'
#endif
         ENDIF
      ENDDO

      RETURN
      END SUBROUTINE

      SUBROUTINE SLOPELIMITER3_3NBORS(L)

!.....Use appropriate modules

      USE SIZES, ONLY : SZ,myproc
      USE GLOBAL
      USE DG


      IMPLICIT NONE

!.....Declare local variables

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

!.....Do nothing if this or a neighboring element is dry.
!     IF(WDFLG(L).EQ.0) RETURN
!     DO I = 1,3
!     EL_N = EL_NBORS(I,L)
!     IF(WDFLG(EL_N).EQ.0) RETURN
!     ENDDO
      
!.....Retrieve the barycenter info of this element
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

!.....Retrieve the nodal info

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
!     IF(HMIN.LE.0.01D0) RETURN

!.....Retrieve info at barycenters and edge midpoints
      DO I = 1,3
         EL_N = EL_NBORS(I,L)

!.......Edge midpoints
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
         
!.......Barycenters
!     IF(EL_N == 0) THEN
!     If there is no neighboring element on this side (e.g. domain boundary),
!     locate a barycenter at the mirroring point and use extrapolated values.
!     but all elements in this subroutine will have 3 neighbors
!     XB(I) = 2.D0*XM(I) - XB(0)
!     YB(I) = 2.D0*YM(I) - YB(0)
!     ZB(I) = 2.D0*ZM(I) - ZB(0)
!     DB(I) = 2.D0*DM(I) - DB(0)
!     QXB(I) = 2.D0*QXM(I) - QXB(0)
!     QYB(I) = 2.D0*QYM(I) - QYB(0)
!     ELSE
         XB(I) = XBC(EL_N)
         YB(I) = YBC(EL_N)
         ZB(I) = ZE(1,EL_N,IRK+1)
         DB(I) = HB(1,EL_N,1)
         QXB(I) = QX(1,EL_N,IRK+1)
         QYB(I) = QY(1,EL_N,IRK+1)
!     ENDIF

!.......Edge normal vector
         GED = NELED(I,L)
!.......Shouldn't be normal vector (EJK) see Cockburn and Shu reference
!     NXX(I) = COSNX(GED)
!     NYY(I) = SINNX(GED)
         NXX(I) = XM(I)-XB(I)
         IF (NXX(I).NE.0) NXX(I) = NXX(I)/ABS(NXX(I))
         NYY(I) = YM(I)-YB(I)
         IF (NYY(I).NE.0) NYY(I) = NYY(I)/ABS(NYY(I))
!     if (myproc.eq.42.and.l.eq.4229) then
!     write(42,*) 'slope limiter ',i,el_n,qxm(i),qxb(i)
!     endif
      ENDDO


!.....Compute alpha1 and alpha2      

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
         
!.......Compute the geometric factors alpha1 and alpha2. If either one
!.......is negative that element must be set to p = 0.

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
      
!.....Set the variable vectors at the baricenter of element 0

      WB_O(1,0) = ZB(0)
      WB_O(2,0) = QXB(0)
      WB_O(3,0) = QYB(0)

!.....Start computing deltas for each edge

      DO I = 1,3

!.......Compute eigen values and vectors at the midpoint using the variables at the barycenter

         HTB = ZB(0)+DB(0)
         CALL HYDRO_EIGEN_VALUES(HTB, QXB(0)/HTB, QYB(0)/HTB,
     &        NXX(I), NYY(I), EIGVAL, RI, LE)

!.......Set variable vectors

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

!.......Transform original variables into the characteristic space

         DO J = 0,2
            WB_C(1,J)
     &           = LE(1,1)*WB_O(1,J) + LE(1,2)*WB_O(2,J) + LE(1,3)*WB_O(3,J)
            WB_C(2,J)
     &           = LE(2,1)*WB_O(1,J) + LE(2,2)*WB_O(2,J) + LE(2,3)*WB_O(3,J)
            WB_C(3,J)
     &           = LE(3,1)*WB_O(1,J) + LE(3,2)*WB_O(2,J) + LE(3,3)*WB_O(3,J)
         ENDDO
         WM_C(1)= LE(1,1)*WM_O(1) + LE(1,2)*WM_O(2) + LE(1,3)*WM_O(3)
         WM_C(2)= LE(2,1)*WM_O(1) + LE(2,2)*WM_O(2) + LE(2,3)*WM_O(3)
         WM_C(3)= LE(3,1)*WM_O(1) + LE(3,2)*WM_O(2) + LE(3,3)*WM_O(3)

!.......Compute W_TILDA
         W_TILDA(1) = WM_C(1) - WB_C(1,0)
         W_TILDA(2) = WM_C(2) - WB_C(2,0)
         W_TILDA(3) = WM_C(3) - WB_C(3,0)
         
!.......Compute DELTA_W
         DELTA_W(1)
     &        = ALPHA1(I)*(WB_C(1,1) - WB_C(1,0))
     &        + ALPHA2(I)*(WB_C(1,2) - WB_C(1,0))
         DELTA_W(2)
     &        = ALPHA1(I)*(WB_C(2,1) - WB_C(2,0))
     &        + ALPHA2(I)*(WB_C(2,2) - WB_C(2,0))
         DELTA_W(3)
     &        = ALPHA1(I)*(WB_C(3,1) - WB_C(3,0))
     &        + ALPHA2(I)*(WB_C(3,2) - WB_C(3,0))

!.......Apply the TVB modified minmod function
!.......Switched to minmod function (EJK)

!     R_SQ = (XM(I) - XB(0))*(XM(I) - XB(0))
!     &       + (YM(I) - YB(0))*(YM(I) - YB(0))

         NOT_LIMITED = 0
!     if (myproc.eq.42.and.l.eq.4229) then
!     write(42,*) 'slope limiter ',delta_w(1),delta_w(2),
!     $          delta_w(3),w_tilda(1),w_tilda(2),w_tilda(3)
!     endif
         DO K = 1,3             ! Loop for components of variable vectors
!     IF(ABS(W_TILDA(K)) < SL2_M*HLIMIT*R_SQ) THEN
!     DELTA_C(K,I) = W_TILDA(K)
!     ELSE
            IF((W_TILDA(K)*DELTA_W(K) >= 0.D0)
     &           .AND.(ABS(W_TILDA(K)) > ZERO)) THEN
               S = W_TILDA(K)/ABS(W_TILDA(K))
!.......Hard-wired SL2_NYU = 1.5 as taken in Cockburn and Shu reference
               DELTA_C(K,I) = S*MIN(
     &              ABS(W_TILDA(K)),ABS(1.5*DELTA_W(K)))
               IF (DELTA_C(K,I).EQ.(S*ABS(W_TILDA(K))))
     &              NOT_LIMITED(K) = NOT_LIMITED(K) + 1
            ELSE
               DELTA_C(K,I) = 0.D0
            ENDIF
!     ENDIF
         ENDDO
      ENDDO                     ! Loop for edges ends here
      
!.....If none of the variables were limited then end
      
      IF ((NOT_LIMITED(1).EQ.3).AND.(NOT_LIMITED(2).EQ.3)
     &     .AND.(NOT_LIMITED(3).EQ.3)) THEN
         RETURN
         
!.....Else perform limiting on characteristic variables

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

!.........Computed limited characteristic variables at the corner nodes

            WC_C(I,1) = WB_C(I,0) - DELTA_HAT(I,1) + DELTA_HAT(I,2)
     &           + DELTA_HAT(I,3)
            WC_C(I,2) = WB_C(I,0) + DELTA_HAT(I,1) - DELTA_HAT(I,2)
     &           + DELTA_HAT(I,3)
            WC_C(I,3) = WB_C(I,0) + DELTA_HAT(I,1) + DELTA_HAT(I,2)
     &           - DELTA_HAT(I,3)
         enddo
!.........Transform back to the original space
         do i=1,3               !loop over corners
            WC_O(1,I) = RI(1,1)*WC_C(1,I) + RI(1,2)*WC_C(2,I)
     &           + RI(1,3)*WC_C(3,I)
            WC_O(2,I) = RI(2,1)*WC_C(1,I) + RI(2,2)*WC_C(2,I)
     &           + RI(2,3)*WC_C(3,I)
            WC_O(3,I) = RI(3,1)*WC_C(1,I) + RI(3,2)*WC_C(2,I)
     &           + RI(3,3)*WC_C(3,I)

!     if (myproc.eq.42.and.l.eq.4229) then
!     write(42,*) 'slope limiter ',wc_c(i,1),wc_c(i,2),wc_c(i,3)
!     endif
         ENDDO                 

         
!.......If the limited water depth < 0, don't apply this slope limiting

         IF((WC_O(1,1)+DC(1)).LE.0.D0.OR.(WC_O(1,2)+DC(2)).LE.0.D0.OR.
     &        (WC_O(1,2)+DC(2)).LE.0.D0) RETURN
         
!.......Else compute new modal dofs

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

!.....Use appropriate modules

      USE SIZES, ONLY : SZ
      USE GLOBAL
      USE DG

      IMPLICIT NONE

!.....Declare local variables

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


!.....Find two neighboring elements
      NNBORS = 0
      DO I = 1,3
         EL_N = EL_NBORS(I,L)
         IF(EL_N.NE.0) THEN
            NNBORS = NNBORS + 1
            EL_NN(NNBORS) = EL_N
         ENDIF
      ENDDO

!.....Do nothing if a neighboring element is dry.
!     IF(WDFLG(L).EQ.0) RETURN
!     DO I = 1,NNBORS
!     EL_N = EL_NN(I)
!     IF(WDFLG(EL_N).EQ.0) RETURN
!     ENDDO
      
!.....Retrieve the barycenter info of this element
      XB(0) = XBC(L)
      YB(0) = YBC(L)
      ZB(0) = ZE(1,L,IRK+1)
      DB(0) = HB(1,L,1)
      QXB(0) = QX(1,L,IRK+1)
      QYB(0) = QY(1,L,IRK+1)
      
!.....Retrieve the nodal info

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
!     IF(HMIN.LE.0.01D0) RETURN

!.....Retrieve info at edge midpoints
      DO I = 1,3
!.......Edge midpoints
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
         
!.......Edge normal vector
         GED = NELED(I,L)
!     NXX(I) = -COSNX(GED)
!     NYY(I) = -SINNX(GED)
         NXX(I) = XM(I)-XB(I)
         IF (NXX(I).NE.0) NXX(I) = NXX(I)/ABS(NXX(I))
         NYY(I) = YM(I)-YB(I)
         IF (NYY(I).NE.0) NYY(I) = NYY(I)/ABS(NYY(I))
      ENDDO

!.....Retreive barycenter info
      DO I = 1,NNBORS
         EL_N = EL_NN(I)
         XB(I) = XBC(EL_N)
         YB(I) = YBC(EL_N)
         ZB(I) = ZE(1,EL_N,IRK+1)
         DB(I) = HB(1,EL_N,1)
         QXB(I) = QX(1,EL_N,IRK+1)
         QYB(I) = QY(1,EL_N,IRK+1)
      ENDDO

!.....Compute alpha1 and alpha2      
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

!.....Set the variable vectors at the baricenter of element 0
      WB_O(1,0) = ZB(0)
      WB_O(2,0) = QXB(0)
      WB_O(3,0) = QYB(0)

!.....Start computing deltas for each edge
      DO I = 1,3

!.......Compute eigen values and vectors at the midpoint using the variables at the barycenter
         HTB = ZB(0)+DB(0)
         CALL HYDRO_EIGEN_VALUES(HTB, QXB(0)/HTB, QYB(0)/HTB,
     &        NXX(I), NYY(I), EIGVAL, RI, LE)

!.......Set variable vectors
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

!.......Transform original variables into the characteristic space
         DO J = 0,2
            WB_C(1,J)
     &           = LE(1,1)*WB_O(1,J) + LE(1,2)*WB_O(2,J) + LE(1,3)*WB_O(3,J)
            WB_C(2,J)
     &           = LE(2,1)*WB_O(1,J) + LE(2,2)*WB_O(2,J) + LE(2,3)*WB_O(3,J)
            WB_C(3,J)
     &           = LE(3,1)*WB_O(1,J) + LE(3,2)*WB_O(2,J) + LE(3,3)*WB_O(3,J)
         ENDDO
         WM_C(1)= LE(1,1)*WM_O(1) + LE(1,2)*WM_O(2) + LE(1,3)*WM_O(3)
         WM_C(2)= LE(2,1)*WM_O(1) + LE(2,2)*WM_O(2) + LE(2,3)*WM_O(3)
         WM_C(3)= LE(3,1)*WM_O(1) + LE(3,2)*WM_O(2) + LE(3,3)*WM_O(3)

!.......Compute W_TILDA
         W_TILDA(1) = WM_C(1) - WB_C(1,0)
         W_TILDA(2) = WM_C(2) - WB_C(2,0)
         W_TILDA(3) = WM_C(3) - WB_C(3,0)
         
!.......Compute DELTA_W
         DELTA_W(1)
     &        = ALPHA1(I)*(WB_C(1,1) - WB_C(1,0))
     &        + ALPHA2(I)*(WB_C(1,2) - WB_C(1,0))
         DELTA_W(2)
     &        = ALPHA1(I)*(WB_C(2,1) - WB_C(2,0))
     &        + ALPHA2(I)*(WB_C(2,2) - WB_C(2,0))
         DELTA_W(3)
     &        = ALPHA1(I)*(WB_C(3,1) - WB_C(3,0))
     &        + ALPHA2(I)*(WB_C(3,2) - WB_C(3,0))

!.......Apply the TVB modified minmod function
!     R_SQ = (XM(I) - XB(0))*(XM(I) - XB(0))
!     &       + (YM(I) - YB(0))*(YM(I) - YB(0))
         NOT_LIMITED=0

         DO K = 1,3             ! Loop for components of variable vectors
!     IF(ABS(W_TILDA(K)) < SL2_M*HLIMIT*R_SQ) THEN
!     DELTA_C(K,I) = W_TILDA(K)
!     ELSE
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
!     ENDIF
         ENDDO
      ENDDO
!.....If none of the variables were limited then end
      
      IF ((NOT_LIMITED(1).EQ.3).AND.(NOT_LIMITED(2).EQ.3)
     &     .AND.(NOT_LIMITED(3).EQ.3)) THEN
         RETURN
         
!.....Else perform limiting on characteristic variables

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

!.........Computed limited characteristic variables at the corner nodes

            WC_C(I,1) = WB_C(I,0) - DELTA_HAT(I,1) + DELTA_HAT(I,2)
     &           + DELTA_HAT(I,3)
            WC_C(I,2) = WB_C(I,0) + DELTA_HAT(I,1) - DELTA_HAT(I,2)
     &           + DELTA_HAT(I,3)
            WC_C(I,3) = WB_C(I,0) + DELTA_HAT(I,1) + DELTA_HAT(I,2)
     &           - DELTA_HAT(I,3)
         enddo                  ! End variable component loop
!.........Transform back to the original space
         do i=1,3

            WC_O(1,I) = RI(1,1)*WC_C(1,I) + RI(1,2)*WC_C(2,I)
     &           + RI(1,3)*WC_C(3,I)
            WC_O(2,I) = RI(2,1)*WC_C(1,I) + RI(2,2)*WC_C(2,I)
     &           + RI(2,3)*WC_C(3,I)
            WC_O(3,I) = RI(3,1)*WC_C(1,I) + RI(3,2)*WC_C(2,I)
     &           + RI(3,3)*WC_C(3,I)

         ENDDO 
         
!.......If the limited water depth < 0, don't apply this slope limiting

         IF((WC_O(1,1)+DC(1)).LE.0.D0.OR.(WC_O(1,2)+DC(2)).LE.0.D0.OR.
     &        (WC_O(1,2)+DC(2)).LE.0.D0) RETURN
         
!.......Else compute new modal dofs

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

!.....Use appropriate modules

      USE SIZES, ONLY : SZ
      USE GLOBAL
      USE DG

      IMPLICIT NONE

!.....Declare local variables

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

!.....Retrieve the barycenter info of this element
      XB(0) = XBC(L)
      YB(0) = YBC(L)
      ZB(0) = ZE(1,L,IRK+1)
      DB(0) = HB(1,L,1)
      QXB(0) = QX(1,L,IRK+1)
      QYB(0) = QY(1,L,IRK+1)
      
!.....Retrieve the nodal info
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

!.....Retrieve info at barycenters and edge midpoints
      DO I = 1,3
         EL_N = EL_NBORS(I,L)

!.......Edge midpoints
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
         
!.......Barycenters
         IF(EL_N == 0) THEN
!     If there is no neighboring element on this side (e.g. domain boundary),
!     locate a barycenter at the mirroring point and use extrapolated values.
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

!.......Edge normal vector
         GED = NELED(I,L)
         NXX(I) = COSNX(GED)
         NYY(I) = SINNX(GED)
      ENDDO

!.....Compute alpha1 and alpha2      
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


!.....Set the variable vectors at the baricenter of element 0
      WB_O(1,0) = ZB(0)
      WB_O(2,0) = QXB(0)
      WB_O(3,0) = QYB(0)

!.....Start computing deltas for each edge
      DO I = 1,3

!.......Compute eigen values and vectors at the midpoint using the variables at the barycenter
         HTB = ZB(0)+DB(0)
         CALL HYDRO_EIGEN_VALUES(HTB, QXB(0)/HTB, QYB(0)/HTB,
     &        NXX(I), NYY(I), EIGVAL, RI, LE)

!.......Set variable vectors
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

!.......Transform original variables into the characteristic space
         DO J = 0,2
            WB_C(1,J)
     &           = LE(1,1)*WB_O(1,J) + LE(1,2)*WB_O(2,J) + LE(1,3)*WB_O(3,J)
            WB_C(2,J)
     &           = LE(2,1)*WB_O(1,J) + LE(2,2)*WB_O(2,J) + LE(2,3)*WB_O(3,J)
            WB_C(3,J)
     &           = LE(3,1)*WB_O(1,J) + LE(3,2)*WB_O(2,J) + LE(3,3)*WB_O(3,J)
         ENDDO
         WM_C(1)= LE(1,1)*WM_O(1) + LE(1,2)*WM_O(2) + LE(1,3)*WM_O(3)
         WM_C(2)= LE(2,1)*WM_O(1) + LE(2,2)*WM_O(2) + LE(2,3)*WM_O(3)
         WM_C(3)= LE(3,1)*WM_O(1) + LE(3,2)*WM_O(2) + LE(3,3)*WM_O(3)

!.......Compute W_TILDA
         W_TILDA(1) = WM_C(1) - WB_C(1,0)
         W_TILDA(2) = WM_C(2) - WB_C(2,0)
         W_TILDA(3) = WM_C(3) - WB_C(3,0)
         
!.......Compute DELTA_W
         DELTA_W(1)
     &        = ALPHA1(I)*(WB_C(1,1) - WB_C(1,0))
     &        + ALPHA2(I)*(WB_C(1,2) - WB_C(1,0))
         DELTA_W(2)
     &        = ALPHA1(I)*(WB_C(2,1) - WB_C(2,0))
     &        + ALPHA2(I)*(WB_C(2,2) - WB_C(2,0))
         DELTA_W(3)
     &        = ALPHA1(I)*(WB_C(3,1) - WB_C(3,0))
     &        + ALPHA2(I)*(WB_C(3,2) - WB_C(3,0))

!.......Apply the TVB modified minmod function
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

!.......Transform back to the original space
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

!     Set the variable vector at the midpoint 1. This is used later to judge
!     whether the variables should be limited.
      WM_O(1) = ZM(1)
      WM_O(2) = QXM(1)
      WM_O(3) = QYM(1)

!     Find the possibly limited DELTA and variable vectors at the corner nodes
      DO I = 1,3                ! Loop for variable components
!     Find the possibly limited DELTA
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

!     Compute the possibly limited variable vectors at the corner nodes
         IF(( DELTA_HAT(I,1) < (WM_O(I) - WB_O(I,0) - ZERO) ).OR.
     &        ( DELTA_HAT(I,1) > (WM_O(I) - WB_O(I,0) + ZERO) )) THEN
!     Update the variables only if they need to be limited.
!     If the condition in the if statement is false, the following 
!     isn't applied and therefore the quradratic and higher order 
!     coefficients are preserved.

!     Computing limited variables at the corner nodes
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

!     If the limited water depth is less 0, don't apply this slope limiting.
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

!***********************************************************************
!     
!     SUBROUTINE HYDRO_EIGEN_VALUES
!     
!     This subroutine computes the eigen values and vectors of hydro
!     dynamics equations.
!     
!     Borrowed from ROE_FLUX() - 06 Feb 2008, S.B.
!     
!***********************************************************************
      SUBROUTINE HYDRO_EIGEN_VALUES(H, U, V, NX, NY, EIGVAL, RI, LE)
      USE SIZES, ONLY : SZ
      USE GLOBAL, ONLY : IFNLCT, G

      IMPLICIT NONE

      REAL(SZ) H,U,V,NX,NY,C,DTM
      REAL(SZ) EIGVAL(3)        ! eigen values
      REAL(SZ) RI(3,3)          ! columns are the right eigen vectors
      REAL(SZ) LE(3,3)          ! columns are the left eigen vectors

      C = SQRT(G*H)

!.....Evaluate the eigenvalues at the Roe averaged variables

      EIGVAL(2) = (U*NX + V*NY)*IFNLCT
      EIGVAL(1) = EIGVAL(2) + C
      EIGVAL(3) = EIGVAL(2) - C
      
!.....Evaluate right eigenvectors at Roe averaged variables

      RI(1,1) = 1.D0
      RI(2,1) = U*IFNLCT + (IFNLCT*C + (1-IFNLCT)*G/C)*NX
      RI(3,1) = V*IFNLCT + (IFNLCT*C + (1-IFNLCT)*G/C)*NY

      RI(1,2) = 0.D0
      RI(2,2) = -NY
      RI(3,2) =  NX

      RI(1,3) = 1.D0
      RI(2,3) = U*IFNLCT - (IFNLCT*C + (1-IFNLCT)*G/C)*NX
      RI(3,3) = V*IFNLCT - (IFNLCT*C + (1-IFNLCT)*G/C)*NY

!.....Evaluate left eigenvectors at Roe averaged variables

!.....Determinant
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




!***********************************************************************
!     
!     SUBROUTINE SLOPELIMITER4()
!     
!     Written 2011
!     
!     This subroutine selects the vertex based slope limiter based on
!     a Taylor Polynomial basis, and is consistent with p_adaptation.F
!     
!     cem
!     
!***********************************************************************

      SUBROUTINE SLOPELIMITER4() 

!.....Use appropriate modules

      USE GLOBAL
      USE DG

#ifdef CMPI
      USE MESSENGER
      USE MESSENGER_ELEM
#endif

      IMPLICIT NONE

      Integer k,ll,ss,lll,ell,bb,maxneigh

      REAL(SZ) fd,marea,fde

!.....We work over the master element
!.....Set initial values

      fd = dg%slope_weight         ! reduces diffusion fd = 1 => full diffusion
      fde = fd                  ! add weight for lower order pieces (fd<1 => stronger limiting)     

      DO k=1,mne

         if (dg%dofs(k).gt.1) then

            DO ll = 1,dg%dofs(k)

               dg%ZE(ll,k,dg%NRK+2) = dg%ZE(ll,k,dg%IRK+1)
               dg%QX(ll,k,dg%NRK+2) = dg%QX(ll,k,dg%IRK+1)
               dg%QY(ll,k,dg%NRK+2) = dg%QY(ll,k,dg%IRK+1)

#ifdef TRACE
               dg%iota(ll,k,dg%NRK+2) = dg%iota(ll,k,dg%IRK+1)
#endif

#ifdef CHEM
               dg%iota(ll,k,dg%NRK+2) = dg%iota(ll,k,dg%IRK+1)
               dg%iota2(ll,k,dg%NRK+2) = dg%iota2(ll,k,dg%IRK+1)
#endif

            ENDDO

         elseif (dg%dofs(k).eq.1) then

            cycle

         endif

      ENDDO   

      marea = 2.D0              !master elements area

!.....Convert initial values to the Taylor basis (multiply by dg%Nmatrix) on base element


      dg%ZEtaylor = 0.D0 
      dg%QXtaylor = 0.D0
      dg%QYtaylor = 0.D0

#ifdef TRACE
      dg%iotataylor = 0.D0
#endif

#ifdef CHEM
      dg%iotataylor = 0.D0
      dg%iota2taylor = 0.D0
#endif

      do k=1,mne

         if (dg%dofs(k).gt.1) then

            do ll=1,dg%dofs(k)
               
               do ss=1,dg%dofs(k)

                  dg%ZEtaylor(k,ll,1) = dg%ZEtaylor(k,ll,1) + dg%Nmatrix(k,ll,ss,dg%dofs(k)) * dg%ZE(ss,k,dg%nrk+2)
                  dg%QXtaylor(k,ll,1) = dg%QXtaylor(k,ll,1) + dg%Nmatrix(k,ll,ss,dg%dofs(k))* dg%QX(ss,k,dg%nrk+2)
                  dg%QYtaylor(k,ll,1) = dg%QYtaylor(k,ll,1) + dg%Nmatrix(k,ll,ss,dg%dofs(k)) * dg%QY(ss,k,dg%nrk+2)

#ifdef TRACE
                  dg%iotataylor(k,ll,1) = dg%iotataylor(k,ll,1) + dg%Nmatrix(k,ll,ss,dg%dofs(k)) * dg%iota(ss,k,dg%nrk+2)
#endif                     
                  
#ifdef CHEM
                  dg%iotataylor(k,ll,1) = dg%iotataylor(k,ll,1) + dg%Nmatrix(k,ll,ss,dg%dofs(k)) * dg%iota(ss,k,dg%nrk+2)
                  dg%iota2taylor(k,ll,1) = dg%iota2taylor(k,ll,1) + dg%Nmatrix(k,ll,ss,dg%dofs(k)) * dg%iota2(ss,k,dg%nrk+2)
#endif
                  
               enddo

            enddo

         elseif (dg%dofs(k).eq.1) then

            cycle

         endif

      enddo

!.....Find values at vertices of base elements and neighbors


      dg%ZEmax = -100.D0
      dg%QXmax = -100.D0
      dg%QYmax = -100.D0
      dg%ZEmin = 100.D0
      dg%QXmin = 100.D0
      dg%QYmin = 100.D0

#ifdef TRACE
      dg%iotamax = -100.D0
      dg%iotamin = 100.D0
#endif

#ifdef CHEM
      dg%iotamax = -100.D0
      dg%iota2max = -100.D0
      dg%iotamin = 100.D0
      dg%iota2min = 100.D0
#endif


#ifdef CMPI

      CALL UPDATER_ELEM_MOD2(dg%ZEtaylor,dg%QXtaylor,dg%QYtaylor,1,3)

#ifdef TRACE
      CALL UPDATER_ELEM_MOD2(dg%iotataylor,dg%QXtaylor,dg%QYtaylor,1,2 )
#endif

#ifdef CHEM
      CALL UPDATER_ELEM_MOD2(dg%iotataylor,dg%iota2taylor,dg%QYtaylor,1,2 )
#endif

      CALL UPDATERV(dg%ZEmin)
      CALL UPDATERV(dg%ZEmax)
      CALL UPDATERV(dg%QXmin)
      CALL UPDATERV(dg%QXmax)
      CALL UPDATERV(dg%QYmin)
      CALL UPDATERV(dg%QYmax)

#ifdef TRACE
      CALL UPDATERV(dg%iotamax)
      CALL UPDATERV(dg%iotamin)
#endif

#ifdef CHEM
      CALL UPDATERV(dg%iotamax)
      CALL UPDATERV(dg%iotamin)
      CALL UPDATERV(dg%iota2max)
      CALL UPDATERV(dg%iota2min)
#endif

#endif

      do ell=1,mnp

         do ll=1,minval(dg%dofs(neigh_elem(ell,1:nneigh_elem(ell))))

!.....Find max and min values over polynomial coefficients

            dg%ZEmax(ell,ll) = max(maxval( dg%ZEtaylor(NEIGH_ELEM(ell,1:nneigh_elem(ell)),ll,1) )
     &           , dg%ZEmax(ell,ll))
            dg%QXmax(ell,ll) = max(maxval( dg%QXtaylor(NEIGH_ELEM(ell,1:nneigh_elem(ell)),ll,1) )
     &           , dg%QXmax(ell,ll))
            dg%QYmax(ell,ll) = max(maxval( dg%QYtaylor(NEIGH_ELEM(ell,1:nneigh_elem(ell)),ll,1) )
     &           , dg%QYmax(ell,ll))
            dg%ZEmin(ell,ll) = min(minval( dg%ZEtaylor(NEIGH_ELEM(ell,1:nneigh_elem(ell)),ll,1) )
     &           , dg%ZEmin(ell,ll))
            dg%QXmin(ell,ll) = min(minval( dg%QXtaylor(NEIGH_ELEM(ell,1:nneigh_elem(ell)),ll,1) )
     &           , dg%QXmin(ell,ll))
            dg%QYmin(ell,ll) = min(minval( dg%QYtaylor(NEIGH_ELEM(ell,1:nneigh_elem(ell)),ll,1) )
     &           , dg%QYmin(ell,ll))

#ifdef TRACE
            dg%iotamax(ell,ll) = max(maxval( dg%iotataylor(NEIGH_ELEM(ell,1:nneigh_elem(ell)),ll,1) )
     &           , dg%iotamax(ell,ll))
            dg%iotamin(ell,ll) = min(minval( dg%iotataylor(NEIGH_ELEM(ell,1:nneigh_elem(ell)),ll,1) )
     &           , dg%iotamin(ell,ll))
#endif

#ifdef CHEM
            dg%iotamax(ell,ll) = max(maxval( dg%iotataylor(NEIGH_ELEM(ell,1:nneigh_elem(ell)),ll,1) )
     &           , dg%iotamax(ell,ll))
            dg%iota2max(ell,ll) = max(maxval( dg%iota2taylor(NEIGH_ELEM(ell,1:nneigh_elem(ell)),ll,1) )
     &           , dg%iota2max(ell,ll))
            dg%iotamin(ell,ll) = min(minval( dg%iotataylor(NEIGH_ELEM(ell,1:nneigh_elem(ell)),ll,1) )
     &           , dg%iotamin(ell,ll))
            dg%iota2min(ell,ll) = min(minval( dg%iota2taylor(NEIGH_ELEM(ell,1:nneigh_elem(ell)),ll,1) )
     &           , dg%iota2min(ell,ll))
#endif
            
         enddo

      enddo


#ifdef CMPI


      CALL UPDATERV(dg%ZEmin)
      CALL UPDATERV(dg%ZEmax)
      CALL UPDATERV(dg%QXmin)
      CALL UPDATERV(dg%QXmax)
      CALL UPDATERV(dg%QYmin)
      CALL UPDATERV(dg%QYmax)

#ifdef TRACE
      CALL UPDATERV(dg%iotamax)
      CALL UPDATERV(dg%iotamin)
#endif

#ifdef CHEM
      CALL UPDATERV(dg%iotamax)
      CALL UPDATERV(dg%iotamin)
      CALL UPDATERV(dg%iota2max)
      CALL UPDATERV(dg%iota2min)
#endif

#endif


!.....Must generate linear recostructions at vertices

      dg%ZEtaylorvert = 0.D0
      dg%QXtaylorvert = 0.D0
      dg%Qytaylorvert = 0.D0

#ifdef TRACE
      dg%iotataylorvert = 0.D0
#endif

#ifdef CHEM
      dg%iotataylorvert = 0.D0
      dg%iota2taylorvert = 0.D0
#endif

      do k=1,mne

         if (dg%dofs(k).gt.1) then

            do ll=1,dg%dofs(k)

               do lll=1,3

                  if (ll.eq.1) then

                     dg%ZEtaylorvert(k,ll,lll) = dg%ZEtaylor(k,ll,1) + dg%ZEtaylor(k,ll+1,1)*( dg%xi2vert(k,lll) -dg%xi2BCb(k) )
     &                    + dg%ZEtaylor(k,ll+2,1)*( dg%xi1vert(k,lll) - dg%xi1BCb(k) )
                     dg%QXtaylorvert(k,ll,lll) = dg%QXtaylor(k,ll,1) + dg%QXtaylor(k,ll+1,1)*( dg%xi2vert(k,lll) -dg%xi2BCb(k) )
     &                    + dg%QXtaylor(k,ll+2,1)*( dg%xi1vert(k,lll) - dg%xi1BCb(k) )
                     dg%Qytaylorvert(k,ll,lll) = dg%QYtaylor(k,ll,1) + dg%QYtaylor(k,ll+1,1)*( dg%xi2vert(k,lll) -dg%xi2BCb(k) ) 
     &                    + dg%QYtaylor(k,ll+2,1)*( dg%xi1vert(k,lll) - dg%xi1BCb(k) )

#ifdef TRACE
                     dg%iotataylorvert(k,ll,lll) = dg%iotataylor(k,ll,1) + 
     &                    dg%iotataylor(k,ll+1,1)*( dg%xi2vert(k,lll) - dg%xi2BCb(k) )
     &                    + dg%iotataylor(k,ll+2,1)*( dg%xi1vert(k,lll) - dg%xi1BCb(k) )
#endif

#ifdef CHEM
                     dg%iotataylorvert(k,ll,lll) = dg%iotataylor(k,ll,1) + 
     &                    dg%iotataylor(k,ll+1,1)*( dg%xi2vert(k,lll) - dg%xi2BCb(k) )
     &                    + dg%iotataylor(k,ll+2,1)*( dg%xi1vert(k,lll) - dg%xi1BCb(k) )
                     dg%iota2taylorvert(k,ll,lll) = dg%iota2taylor(k,ll,1) + 
     &                    dg%iota2taylor(k,ll+1,1)*( dg%xi2vert(k,lll) - dg%xi2BCb(k) )
     &                    + dg%iota2taylor(k,ll+2,1)*( dg%xi1vert(k,lll) - dg%xi1BCb(k) )
#endif

                  elseif (ll.gt.1) then

                     dg%ZEtaylorvert(k,ll,lll) = dg%ZEtaylor(k,ll,1) + 
     &                    dg%ZEtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) ),1)*( dg%xi2vert(k,lll) - dg%xi2BCb(k) ) 
     &                    + dg%ZEtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) )+1,1)*( dg%xi1vert(k,lll) - dg%xi1BCb(k) )
                     dg%QXtaylorvert(k,ll,lll) = dg%QXtaylor(k,ll,1) + 
     &                    dg%QXtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) ),1)*( dg%xi2vert(k,lll) - dg%xi2BCb(k) ) 
     &                    + dg%QXtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) )+1,1)*( dg%xi1vert(k,lll) - dg%xi1BCb(k) )
                     dg%Qytaylorvert(k,ll,lll) = dg%QYtaylor(k,ll,1) + 
     &                    dg%QYtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) ),1)*( dg%xi2vert(k,lll) - dg%xi2BCb(k) ) 
     &                    + dg%QYtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) )+1,1)*( dg%xi1vert(k,lll) - dg%xi1BCb(k) )

#ifdef TRACE
                     dg%iotataylorvert(k,ll,lll) = dg%iotataylor(k,ll,1)+
     &                    dg%iotataylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) ),1)*( dg%xi2vert(k,lll) - dg%xi2BCb(k))
     &                    + dg%iotataylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) )+1,1)*( dg%xi1vert(k,lll) - dg%xi1BCb(k) )
#endif

#ifdef CHEM
                     dg%iotataylorvert(k,ll,lll) = dg%iotataylor(k,ll,1)+
     &                    dg%iotataylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) ),1)*( dg%xi2vert(k,lll) - dg%xi2BCb(k))
     &                    + dg%iotataylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) )+1,1)*( dg%xi1vert(k,lll) - dg%xi1BCb(k) )
                     dg%iota2taylorvert(k,ll,lll) = dg%iota2taylor(k,ll,1)+
     &                    dg%iotataylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) ),1)*( dg%xi2vert(k,lll) - dg%xi2BCb(k))
     &                    + dg%iota2taylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) )+1,1)*( dg%xi1vert(k,lll) - dg%xi1BCb(k) )
#endif

                  else

                  endif

               enddo
               
            enddo
            
         elseif (dg%dofs(k).eq.1) then

            cycle

         endif

      enddo

!.....Compute alphas for each variable in each order derivitive


      dg%alphaZE0 = 0.D0
      dg%alphaQX0 = 0.D0
      dg%alphaQY0 = 0.D0

#ifdef TRACE
      dg%alphaiota0 = 0.D0
#endif

#ifdef CHEM
      dg%alphaiota0 = 0.D0
      dg%alphaiota20 = 0.D0
#endif

      do k=1,mne

         if (dg%dofs(k).gt.1) then

            do lll=1,3
               
               do ll=1,dg%dofs(k)

                  if (dg%ZEmin(nm(k,lll),ll).ne.dg%ZEmax(nm(k,lll),ll)) then

                     if ( ( dg%ZEtaylorvert(k,ll,lll).gt.dg%ZEtaylor(k,ll,1) ).and.
     &                    ( abs(dg%ZEtaylorvert(k,ll,lll)-dg%ZEtaylor(k,ll,1)).gt.1.0E-15 ).and.
     &                    ( dg%ZEmax(nm(k,lll),ll).ne.dg%ZEtaylor(k,ll,1) ) ) then  

                        dg%alphaZE0(k,ll,lll) = min(1.D0,  ( dg%ZEmax(nm(k,lll),ll)
     &                       - dg%ZEtaylor(k,ll,1) )/ (dg%ZEtaylorvert(k,ll,lll) - dg%ZEtaylor(k,ll,1)))

                     elseif ( (dg%ZEtaylorvert(k,ll,lll).lt.dg%ZEtaylor(k,ll,1) )
     &                       .and.( abs(dg%ZEtaylorvert(k,ll,lll)-dg%ZEtaylor(k,ll,1)).gt.1.0E-15 ).and.
     &                       ( dg%ZEmin(nm(k,lll),ll).ne.dg%ZEtaylor(k,ll,1) ) ) then 

                        dg%alphaZE0(k,ll,lll) = min( 1.D0,( dg%ZEmin(nm(k,lll),ll)
     &                       - dg%ZEtaylor(k,ll,1) )/( dg%ZEtaylorvert(k,ll,lll)-dg%ZEtaylor(k,ll,1)))
                        
                     elseif ( ( dg%ZEtaylorvert(k,ll,lll).eq.dg%ZEtaylor(k,ll,1) ).or.
     &                       ( abs(dg%ZEtaylorvert(k,ll,lll)-dg%ZEtaylor(k,ll,1)).le.1.0E-15 ) ) then

                        dg%alphaZE0(k,ll,lll) = 1.D0

                     endif

                  else

                     dg%alphaZE0(k,ll,lll) = 1.D0

                  endif

#ifdef TRACE
                  if (dg%iotamin(nm(k,lll),ll).ne.dg%iotamax(nm(k,lll),ll)) then

                     if ( ( dg%iotataylorvert(k,ll,lll).gt.dg%iotataylor(k,ll,1) ).and.
     &                    ( abs(dg%iotataylorvert(k,ll,lll)-dg%iotataylor(k,ll,1)).gt.1.0E-15 ).and.
     &                    ( dg%iotamax(nm(k,lll),ll).ne.dg%iotataylor(k,ll,1) ) ) then  

                        dg%alphaiota0(k,ll,lll) = min(1.D0,  ( dg%iotamax(nm(k,lll),ll)
     &                       - dg%iotataylor(k,ll,1) )/ (dg%iotataylorvert(k,ll,lll) - dg%iotataylor(k,ll,1))) 
                        
                     elseif ( (dg%iotataylorvert(k,ll,lll).lt.dg%iotataylor(k,ll,1) )
     &                       .and.( abs(dg%iotataylorvert(k,ll,lll)-dg%iotataylor(k,ll,1)).gt.1.0E-15 ).and.
     &                       ( dg%iotamin(nm(k,lll),ll).ne.dg%iotataylor(k,ll,1) ) ) then 

                        dg%alphaiota0(k,ll,lll) = min( 1.D0,( dg%iotamin(nm(k,lll),ll)
     &                       - dg%iotataylor(k,ll,1) )/( dg%iotataylorvert(k,ll,lll)-dg%iotataylor(k,ll,1)))

                     elseif ( ( dg%iotataylorvert(k,ll,lll).eq.dg%iotataylor(k,ll,1) ).or.
     &                       ( abs(dg%iotataylorvert(k,ll,lll)-dg%iotataylor(k,ll,1)).le.1.0E-15 ) ) then

                        dg%alphaiota0(k,ll,lll) = 1.D0

                     endif

                  else

                     dg%alphaiota0(k,ll,lll) = 1.D0

                  endif
#endif
                  
#ifdef CHEM                 
                  if (dg%iotamin(nm(k,lll),ll).ne.dg%iotamax(nm(k,lll),ll)) then

                     if ( ( dg%iotataylorvert(k,ll,lll).gt.dg%iotataylor(k,ll,1) ).and.
     &                    ( abs(dg%iotataylorvert(k,ll,lll)-dg%iotataylor(k,ll,1)).gt.1.0E-15 ).and.
     &                    ( dg%iotamax(nm(k,lll),ll).ne.dg%iotataylor(k,ll,1) ) ) then  

                        dg%alphaiota0(k,ll,lll) = min(1.D0,  ( dg%iotamax(nm(k,lll),ll)
     &                       - dg%iotataylor(k,ll,1) )/ (dg%iotataylorvert(k,ll,lll) - dg%iotataylor(k,ll,1)))

                     elseif ( (dg%iotataylorvert(k,ll,lll).lt.dg%iotataylor(k,ll,1) )
     &                       .and.( abs(dg%iotataylorvert(k,ll,lll)-dg%iotataylor(k,ll,1)).gt.1.0E-15 ).and.
     &                       ( dg%iotamin(nm(k,lll),ll).ne.dg%iotataylor(k,ll,1) ) ) then 

                        dg%alphaiota0(k,ll,lll) = min( 1.D0,( dg%iotamin(nm(k,lll),ll)
     &                       - dg%iotataylor(k,ll,1) )/( dg%iotataylorvert(k,ll,lll)-dg%iotataylor(k,ll,1)))
                        
                     elseif ( ( dg%iotataylorvert(k,ll,lll).eq.dg%iotataylor(k,ll,1) ).or.
     &                       ( abs(dg%iotataylorvert(k,ll,lll)-dg%iotataylor(k,ll,1)).le.1.0E-15 ) ) then

                        dg%alphaiota0(k,ll,lll) = 1.D0

                     endif

                  else

                     dg%alphaiota0(k,ll,lll) = 1.D0

                  endif

                  if (dg%iota2min(nm(k,lll),ll).ne.dg%iota2max(nm(k,lll),ll)) then

                     if ( ( dg%iota2taylorvert(k,ll,lll).gt.dg%iota2taylor(k,ll,1) ).and.
     &                    ( abs(dg%iota2taylorvert(k,ll,lll)-dg%iota2taylor(k,ll,1)).gt.1.0E-15 ).and.
     &                    ( dg%iota2max(nm(k,lll),ll).ne.dg%iota2taylor(k,ll,1) ) ) then  

                        dg%alphaiota20(k,ll,lll) = min(1.D0,  ( dg%iota2max(nm(k,lll),ll)
     &                       - dg%iota2taylor(k,ll,1) )/ (dg%iota2taylorvert(k,ll,lll) - dg%iota2taylor(k,ll,1)))

                     elseif ( (dg%iota2taylorvert(k,ll,lll).lt.dg%iota2taylor(k,ll,1) )
     &                       .and.( abs(dg%iota2taylorvert(k,ll,lll)-dg%iota2taylor(k,ll,1)).gt.1.0E-15 ).and.
     &                       ( dg%iota2min(nm(k,lll),ll).ne.dg%iota2taylor(k,ll,1) ) ) then 

                        dg%alphaiota20(k,ll,lll) = min( 1.D0,( dg%iota2min(nm(k,lll),ll)
     &                       - dg%iota2taylor(k,ll,1) )/( dg%iota2taylorvert(k,ll,lll)-dg%iota2taylor(k,ll,1)))
                        
                     elseif ( ( dg%iota2taylorvert(k,ll,lll).eq.dg%iota2taylor(k,ll,1) ).or.
     &                       ( abs(dg%iota2taylorvert(k,ll,lll)-dg%iota2taylor(k,ll,1)).le.1.0E-15 ) ) then

                        dg%alphaiota20(k,ll,lll) = 1.D0

                     endif

                  else

                     dg%alphaiota20(k,ll,lll) = 1.D0

                  endif
#endif                 

                  if (dg%QXmin(nm(k,lll),ll).ne.dg%QXmax(nm(k,lll),ll)) then

                     if ( ( dg%QXtaylorvert(k,ll,lll).gt.dg%QXtaylor(k,ll,1) ).and.
     &                    ( abs(dg%QXtaylorvert(k,ll,lll)-dg%QXtaylor(k,ll,1)).gt.1.0E-15 ).and.
     &                    ( dg%QXmax(nm(k,lll),ll).ne.dg%QXtaylor(k,ll,1) ) ) then  

                        dg%alphaQX0(k,ll,lll) = min(1.D0,  ( dg%QXmax(nm(k,lll),ll)
     &                       - dg%QXtaylor(k,ll,1) )/ (dg%QXtaylorvert(k,ll,lll) - dg%QXtaylor(k,ll,1)))

                     elseif ( (dg%QXtaylorvert(k,ll,lll).lt.dg%QXtaylor(k,ll,1) )
     &                       .and.( abs(dg%QXtaylorvert(k,ll,lll)-dg%QXtaylor(k,ll,1)).gt.1.0E-15 ).and.
     &                       ( dg%QXmin(nm(k,lll),ll).ne.dg%QXtaylor(k,ll,1) ) ) then 

                        dg%alphaQX0(k,ll,lll) = min( 1.D0,( dg%QXmin(nm(k,lll),ll)
     &                       - dg%QXtaylor(k,ll,1) )/( dg%QXtaylorvert(k,ll,lll)-dg%QXtaylor(k,ll,1)))
                        
                     elseif ( ( dg%QXtaylorvert(k,ll,lll).eq.dg%QXtaylor(k,ll,1) ).or.
     &                       ( abs(dg%QXtaylorvert(k,ll,lll)-dg%QXtaylor(k,ll,1)).le.1.0E-15 ) ) then

                        dg%alphaQX0(k,ll,lll) = 1.D0

                     endif

                  else

                     dg%alphaQX0(k,ll,lll) = 1.D0

                  endif


                  if (dg%QYmin(nm(k,lll),ll).ne.dg%QYmax(nm(k,lll),ll)) then

                     if ( ( dg%QYtaylorvert(k,ll,lll).gt.dg%QYtaylor(k,ll,1) ).and.
     &                    ( abs(dg%QYtaylorvert(k,ll,lll)-dg%QYtaylor(k,ll,1)).gt.1.0E-15 ).and.
     &                    ( dg%QYmax(nm(k,lll),ll).ne.dg%QYtaylor(k,ll,1) ) ) then  

                        dg%alphaQY0(k,ll,lll) = min(1.D0,  ( dg%QYmax(nm(k,lll),ll)
     &                       - dg%QYtaylor(k,ll,1) )/ (dg%QYtaylorvert(k,ll,lll) - dg%QYtaylor(k,ll,1)))

                     elseif ( (dg%QYtaylorvert(k,ll,lll).lt.dg%QYtaylor(k,ll,1) )
     &                       .and.( abs(dg%QYtaylorvert(k,ll,lll)-dg%QYtaylor(k,ll,1)).gt.1.0E-15 ).and.
     &                       ( dg%QYmin(nm(k,lll),ll).ne.dg%QYtaylor(k,ll,1) ) ) then 

                        dg%alphaQY0(k,ll,lll) = min( 1.D0,( dg%QYmin(nm(k,lll),ll)
     &                       - dg%QYtaylor(k,ll,1) )/( dg%QYtaylorvert(k,ll,lll)-dg%QYtaylor(k,ll,1)))                        
                        
                     elseif ( ( dg%QYtaylorvert(k,ll,lll).eq.dg%QYtaylor(k,ll,1) ).or.
     &                       ( abs(dg%QYtaylorvert(k,ll,lll)-dg%QYtaylor(k,ll,1)).le.1.0E-15 ) ) then

                        dg%alphaQY0(k,ll,lll) = 1.D0

                     endif

                  else

                     dg%alphaQY0(k,ll,lll) = 1.D0

                  endif

               enddo            !lll

            enddo               !ll

         elseif (dg%dofs(k).eq.1) then

            cycle

         endif

      enddo

!.....Find the prescribed higher limiters by finding smallest local value

      dg%alphaZE = 0.D0
      dg%alphaQX = 0.D0
      dg%alphaQY = 0.D0

#ifdef TRACE
      dg%alphaiota = 0.D0
#endif

#ifdef CHEM
      dg%alphaiota = 0.D0
      dg%alphaiota2 = 0.D0
#endif

      do k=1,mne

         if (dg%dofs(k).gt.1) then

            do ll=1,dg%dofs(k)		

               dg%alphaZE(k,ll) = minval( dg%alphaZE0(k,ll,:) )
               dg%alphaQX(k,ll) = minval( dg%alphaQX0(k,ll,:) )
               dg%alphaQY(k,ll) = minval( dg%alphaQY0(k,ll,:) )

#ifdef TRACE
               dg%alphaiota(k,ll) = minval( dg%alphaiota0(k,ll,:) )
#endif

#ifdef CHEM
               dg%alphaiota(k,ll) = minval( dg%alphaiota0(k,ll,:) )
               dg%alphaiota2(k,ll) = minval( dg%alphaiota20(k,ll,:) )
#endif

            enddo

         elseif (dg%dofs(k).eq.1) then

            cycle

         endif

      enddo

!.... Choose smallest (minimum) alpha for derivative in x or y

      dg%alphaZEm = 0.D0
      dg%alphaQXm = 0.D0
      dg%alphaQYm = 0.D0

#ifdef TRACE
      dg%alphaiotam = 0.D0
#endif

#ifdef CHEM
      dg%alphaiotam = 0.D0
      dg%alphaiota2m = 0.D0
#endif

      do k=1,mne

         if (dg%dofs(k).gt.1) then
            
            do bb = 1,pdg_el(k)

               if( (bb+1)*(bb+2)/2.le.dg%dofs(k) ) then

                  dg%alphaZEm(k,(bb*(bb+1))/2 + 1) = 
     &                 minval( dg%alphaZE(k,(bb*(bb-1))/2 +1:(bb*(bb+1)) / 2  ) )
                  dg%alphaQXm(k,(bb*(bb+1))/2 + 1) = 
     &                 minval( dg%alphaQX(k,(bb*(bb-1))/2 +1:(bb*(bb+1)) / 2  ) )
                  dg%alphaQYm(k,(bb*(bb+1))/2 + 1) = 
     &                 minval( dg%alphaQY(k,(bb*(bb-1))/2 +1:(bb*(bb+1)) / 2  ) )

#ifdef TRACE
                  dg%alphaiotam(k,(bb*(bb+1))/2 + 1) = 
     &                 minval( dg%alphaiota(k,(bb*(bb-1))/2 +1:(bb*(bb+1)) / 2 ) )
#endif

#ifdef CHEM
                  dg%alphaiotam(k,(bb*(bb+1))/2 + 1) = 
     &                 minval( dg%alphaiota(k,(bb*(bb-1))/2 +1:(bb*(bb+1)) / 2  ) )
                  dg%alphaiota2m(k,(bb*(bb+1))/2 + 1) = 
     &                 minval( dg%alphaiota2(k,(bb*(bb-1))/2 +1:(bb*(bb+1)) / 2 ) )
#endif

               endif
               
            enddo

         elseif (dg%dofs(k).eq.1) then

            cycle

         endif

      enddo

!.....Use max higher derivative values for final limiter value

      dg%alphaZE_max = 0.D0
      dg%alphaQX_max = 0.D0
      dg%alphaQY_max = 0.D0

#ifdef TRACE
      dg%alphaiota_max = 0.D0
#endif

#ifdef CHEM
      dg%alphaiota_max = 0.D0
      dg%alphaiota2_max = 0.D0
#endif

      do k=1,mne

         if (dg%dofs(k).gt.1) then

            do bb =1,pdg_el(k)

               if( (bb+1)*(bb+2)/2.le.dg%dofs(k)) then

                  dg%alphaZE_max(k,(bb*(bb+1))/2 + 1) = fde*maxval( 
     &                 dg%alphaZEm(k,(bb*(bb+1))/2 + 1:dg%dofs(k)) )
                  dg%alphaQX_max(k,(bb*(bb+1))/2 + 1) = fde*maxval( 
     &                 dg%alphaQXm(k,(bb*(bb+1))/2 + 1:dg%dofs(k)) )
                  dg%alphaQY_max(k,(bb*(bb+1))/2 + 1) = fde*maxval( 
     &                 dg%alphaQYm(k,(bb*(bb+1))/2 + 1:dg%dofs(k)) )

#ifdef TRACE
                  dg%alphaiota_max(k,(bb*(bb+1))/2 + 1) = fde*maxval( 
     &                 dg%alphaiotam(k,(bb*(bb+1)/2) +1:dg%dofs(k) ) )
#endif

#ifdef CHEM
                  dg%alphaiota_max(k,(bb*(bb+1))/2 + 1) = fde*maxval( 
     &                 dg%alphaiotam(k,(bb*(bb+1))/2 + 1:dg%dofs(k)) )
                  dg%alphaiota2_max(k,(bb*(bb+1))/2 + 1) = fde*maxval( 
     &                 dg%alphaiota2m(k,(bb*(bb+1))/2 + 1:dg%dofs(k)) )
#endif

               endif

            enddo

         elseif (dg%dofs(k).eq.1) then

            cycle

         endif

      enddo

!.....Limit on the Master element in the Taylor basis, via reconstruction 
!.....of unconstrained solutions with alpha constraints


      dg%limitZE = 0.D0
      dg%limitQX = 0.D0
      dg%limitQY = 0.D0

      dg%lim_count_roll = 0

#ifdef TRACE
      dg%limitiota = 0.D0
#endif

#ifdef CHEM
      dg%limitiota = 0.D0
      dg%limitiota2 = 0.D0
#endif

      do k=1,mne

         dg%lim_count = 0

         if (dg%dofs(k).gt.1) then

            do ll=1,dg%dofs(k)

               if ( ll.eq.1 ) then

                  dg%limitZE(k,ll) = dg%ZEtaylor(k,ll,1)
                  dg%limitQX(k,ll) = dg%QXtaylor(k,ll,1)
                  dg%limitQY(k,ll) = dg%QYtaylor(k,ll,1) 

#ifdef TRACE
                  dg%limitiota(k,ll) = dg%iotataylor(k,ll,1)
#endif

#ifdef CHEM
                  dg%limitiota(k,ll) = dg%iotataylor(k,ll,1)
                  dg%limitiota2(k,ll) = dg%iota2taylor(k,ll,1)
#endif

               elseif ( ll.ge.2 ) then
                  
                  do bb=1,pdg_el(k)

                     if ( ll.le.( ( (bb+1)*(bb+2)) / 2.D0 ).and.(ll.gt.
     &                    (bb*(bb+1)/2.D0) ) ) then

                        dg%limitZE(k,ll) = dg%alphaZE_max(k,(bb*(bb+1))/2 + 1) 
     &                       * dg%ZEtaylor(k,ll,1)
                        dg%limitQX(k,ll) = dg%alphaQX_max(k,(bb*(bb+1))/2 + 1) 
     &                       * dg%ZEtaylor(k,ll,1)
                        dg%limitQY(k,ll) = dg%alphaQY_max(k,(bb*(bb+1))/2 + 1) 
     &                       * dg%ZEtaylor(k,ll,1)


#ifdef TRACE
                        dg%limitiota(k,ll) = dg%alphaiota_max(k,(bb*(bb+1))/2 + 1) 
     &                       * dg%iotataylor(k,ll,1)
#endif

#ifdef CHEM
                        dg%limitiota(k,ll) = dg%alphaiota_max(k,(bb*(bb+1))/2 + 1) 
     &                       * dg%iotataylor(k,ll,1)
                        dg%limitiota2(k,ll) = dg%alphaiota2_max(k,(bb*(bb+1))/2 + 1) 
     &                       * dg%iota2taylor(k,ll,1)
#endif


!$$$  ! Make a counter to track limiting
!$$$  
!$$$  if ( ( dg%alphaZE_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
!$$$  &                       dg%alphaQX_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
!$$$  &                       dg%alphaQY_max(k,(bb*(bb+1))/2 + 1).lt.1.D0 ).and.
!$$$  &                       chem_flag.ne.1.and.tracer_flag.ne.1 ) then
!$$$  
!$$$  dg%lim_count = 1  
!$$$  
!$$$  elseif ( (dg%alphaZE_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
!$$$  &                          dg%alphaQX_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
!$$$  &                          dg%alphaQY_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
!$$$  &                          dg%alphaiota_max(k,(bb*(bb+1))/2 + 1).lt.1.D0 ).and.
!$$$  &                          chem_flag.ne.1.and.tracer_flag.eq.1 ) then
!$$$  
!$$$  dg%lim_count = 1   
!$$$  
!$$$  elseif ( (dg%alphaZE_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
!$$$  &                          dg%alphaQX_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
!$$$  &                          dg%alphaQY_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
!$$$  &                          dg%alphaiota_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
!$$$  &                          dg%alphaiota2_max(k,(bb*(bb+1))/2 + 1).lt.1.D0 ).and.
!$$$  &                          chem_flag.eq.1 ) then
!$$$  
!$$$  dg%lim_count = 1 
!$$$  
!$$$  endif
                        
                     endif
                     
                  enddo
                  
               endif

            enddo

         elseif (dg%dofs(k).eq.1) then

            cycle

         endif

                                !dg%lim_count_roll = dg%lim_count_roll + dg%lim_count

      enddo

!.....Transform back to the Dubiner basis (multiply by dg%NmatrixInv),

      dg%ZEconst =  0.D0
      dg%QXconst =  0.D0
      dg%QYconst =  0.D0

#ifdef TRACE
      dg%iotaconst =  0.D0
#endif

#ifdef CHEM
      dg%iotaconst =  0.D0
      dg%iota2const =  0.D0
#endif

      do k=1,mne

         if (dg%dofs(k).gt.1) then
                                !do lll=1,3

            do ll=1,dg%dofs(k)

               do ss=1,dg%dofs(k)

                  dg%ZEconst(k,ll) = dg%ZEconst(k,ll) + dg%NmatrixInv(k,ll,ss,dg%dofs(k)) 
     &                 * dg%limitZE(k,ss)
                  dg%QXconst(k,ll) = dg%QXconst(k,ll) + dg%NmatrixInv(k,ll,ss,dg%dofs(k)) 
     &                 * dg%limitQX(k,ss)
                  dg%QYconst(k,ll) = dg%QYconst(k,ll) + dg%NmatrixInv(k,ll,ss,dg%dofs(k)) 
     &                 * dg%limitQY(k,ss)

#ifdef TRACE
                  dg%iotaconst(k,ll) = dg%iotaconst(k,ll) + 
     &                 dg%NmatrixInv(k,ll,ss,dg%dofs(k)) * dg%limitiota(k,ss)
#endif

#ifdef CHEM
                  dg%iotaconst(k,ll) = dg%iotaconst(k,ll) + 
     &                 dg%NmatrixInv(k,ll,ss,dg%dofs(k)) * dg%limitiota(k,ss)
                  dg%iota2const(k,ll) = dg%iota2const(k,ll) + 
     &                 dg%NmatrixInv(k,ll,ss,dg%dofs(k)) * dg%limitiota2(k,ss)
#endif


               enddo

            enddo

         elseif (dg%dofs(k).eq.1) then

            cycle

         endif

      enddo

!.....Set limit values

      do k=1,mne

         if (dg%dofs(k).gt.1) then

            do ll = 1,dg%dofs(k)

               dg%ZE(ll,k,dg%irk+1) = dg%ZEconst(k,ll)
               dg%QX(ll,k,dg%irk+1) = dg%QXconst(k,ll)
               dg%QY(ll,k,dg%irk+1) = dg%QYconst(k,ll)

#ifdef TRACE
               dg%iota(ll,k,dg%irk+1) = dg%iotaconst(k,ll)
#endif

#ifdef CHEM
               dg%iota(ll,k,dg%irk+1) = dg%iotaconst(k,ll)
               dg%iota2(ll,k,dg%irk+1) = dg%iota2const(k,ll)
#endif

            enddo

         elseif (dg%dofs(k).eq.1) then

            cycle

         endif

      enddo
      
      return
      end subroutine

#endif

#ifdef SLOPE5

!***********************************************************************
!     
!     SUBROUTINE SLOPELIMITER5()
!     
!     Written by Clint Dawson - 30 June 2010
!     01-10-2011 - cem - adapted for p_enrichment and multicomponent
!     
!     
!***********************************************************************

      SUBROUTINE SLOPELIMITER5(s)

!.....Use appropriate modules

      USE SIZES
      USE GLOBAL
      USE DG

#ifdef CMPI
      USE MESSENGER
#endif

      IMPLICIT NONE

      type (sizes_type) :: s

!.....Declare local variables

      INTEGER L, LL, INC1,INC2,INC3,KDP,NN,IVAR,I,J,kk,k,varnum,bb,varnum_prev
      REAL(SZ) ZEC(3),ZEVERTEX(3),DIF(3),SUMLOC,SUMDIF,SIGNDIF
      REAL(SZ)    DIV,REDFAC,REDMAX,tmp1,tmp2,tmp3,bound
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
      Allocate ( bed_MIN1(NP,s%layers),bed_MAX1(NP,s%layers) )

!     FIND THE MAXIMUM AND MINIMUM OF EACH VARIABLE OVER ALL ELEMENTS 
!     SHARING A NODE

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
         do l=1,s%layers
            bed_MIN1(I,l) = 99999.
            bed_MAX1(I,l) =-99999.
         enddo
#endif

         NO_NBORS = EL_COUNT(I)

         DO J = 1,NO_NBORS
            NBOR_EL = ELETAB(I,1+J)

!     IF(dg%WDFLG(NBOR_EL).EQ.0) CYCLE ! DON'T COUNT DRY ELEMENTS  sb 02/26/07

            ZE_DG(J) = dg%ZE(1,NBOR_EL,dg%IRK+1)
            QX_DG(J) = dg%QX(1,NBOR_EL,dg%IRK+1)
            QY_DG(J) = dg%QY(1,NBOR_EL,dg%IRK+1)

#ifdef TRACE
            iota_DG(J) = dg%iota(1,NBOR_EL,dg%IRK+1)
#endif

#ifdef CHEM
            iota_DG(J) = dg%iota(1,NBOR_EL,dg%IRK+1)
            iota2_DG(J) = dg%iota2(1,NBOR_EL,dg%IRK+1)
#endif


#ifdef SED_LAY
            do l=1,s%layers
               bed_DG(J,l) = dg%bed(1,NBOR_EL,dg%IRK+1,l)
            enddo
#endif

!     
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
            do l=1,s%layers
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
         do l = 1,s%layers
            arraymax => bed_max1(:,l)
            arraymin => bed_min1(:,l)
            CALL UPDATER(arraymax,arraymin,QY_MAX1,2 )
         enddo
#endif

#endif
!     
!     LOOP OVER ELEMENTS TO CALCULATE NEW VERTEX VALUES
!     

      bb = 1

      DO I=1,NE 
         !IF(dg%WDFLG(I).EQ.0) CYCLE ! DON'T COUNT DRY ELEMENTS  sb 02/26/07
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
         varnum = varnum_prev + s%layers
#endif

         DO IVAR=1,varnum

            IF (IVAR.EQ.1) THEN
               ZEC(1)=dg%ZE(1,I,dg%IRK+1)
               ZEC(2)=dg%ZE(2,I,dg%IRK+1)
               ZEC(3)=dg%ZE(3,I,dg%IRK+1)
               ZEMAX1(1)=ZE_MAX1(N1)
               ZEMIN1(1)=ZE_MIN1(N1)
               ZEMAX1(2)=ZE_MAX1(N2)
               ZEMIN1(2)=ZE_MIN1(N2)
               ZEMAX1(3)=ZE_MAX1(N3)
               ZEMIN1(3)=ZE_MIN1(N3)
            ENDIF

            IF (IVAR.EQ.2) THEN
               ZEC(1)=dg%QX(1,I,dg%IRK+1)
               ZEC(2)=dg%QX(2,I,dg%IRK+1)
               ZEC(3)=dg%QX(3,I,dg%IRK+1)
               ZEMAX1(1)=QX_MAX1(N1)
               ZEMIN1(1)=QX_MIN1(N1)
               ZEMAX1(2)=QX_MAX1(N2)
               ZEMIN1(2)=QX_MIN1(N2)
               ZEMAX1(3)=QX_MAX1(N3)
               ZEMIN1(3)=QX_MIN1(N3)
            ENDIF

            IF (IVAR.EQ.3) THEN
               ZEC(1)=dg%QY(1,I,dg%IRK+1)
               ZEC(2)=dg%QY(2,I,dg%IRK+1)
               ZEC(3)=dg%QY(3,I,dg%IRK+1)
               ZEMAX1(1)=QY_MAX1(N1)
               ZEMIN1(1)=QY_MIN1(N1)
               ZEMAX1(2)=QY_MAX1(N2)
               ZEMIN1(2)=QY_MIN1(N2)
               ZEMAX1(3)=QY_MAX1(N3)
               ZEMIN1(3)=QY_MIN1(N3)
            ENDIF

#ifdef TRACE
            IF (IVAR.EQ.4) THEN
               ZEC(1)=dg%iota(1,I,dg%IRK+1)
               ZEC(2)=dg%iota(2,I,dg%IRK+1)
               ZEC(3)=dg%iota(3,I,dg%IRK+1)
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
               ZEC(1)=dg%iota(1,I,dg%IRK+1)
               ZEC(2)=dg%iota(2,I,dg%IRK+1)
               ZEC(3)=dg%iota(3,I,dg%IRK+1)
               ZEMAX1(1)=iota_MAX1(N1)
               ZEMIN1(1)=iota_MIN1(N1)
               ZEMAX1(2)=iota_MAX1(N2)
               ZEMIN1(2)=iota_MIN1(N2)
               ZEMAX1(3)=iota_MAX1(N3)
               ZEMIN1(3)=iota_MIN1(N3)
            ENDIF

            IF (IVAR.EQ.5) THEN
               ZEC(1)=dg%iota2(1,I,dg%IRK+1)
               ZEC(2)=dg%iota2(2,I,dg%IRK+1)
               ZEC(3)=dg%iota2(3,I,dg%IRK+1)
               ZEMAX1(1)=iota2_MAX1(N1)
               ZEMIN1(1)=iota2_MIN1(N1)
               ZEMAX1(2)=iota2_MAX1(N2)
               ZEMIN1(2)=iota2_MIN1(N2)
               ZEMAX1(3)=iota2_MAX1(N3)
               ZEMIN1(3)=iota2_MIN1(N3)
            ENDIF
#endif

#ifdef SED_LAY
            do l=1,s%layers
               if (IVAR.eq.varnum_prev+l) then
                  ZEC(1)=dg%bed(1,I,dg%IRK+1,l)
                  ZEC(2)=dg%bed(2,I,dg%IRK+1,l)
                  ZEC(3)=dg%bed(3,I,dg%IRK+1,l)
                  ZEMAX1(1)=bed_MAX1(N1,l)
                  ZEMIN1(1)=bed_MIN1(N1,l)
                  ZEMAX1(2)=bed_MAX1(N2,l)
                  ZEMIN1(2)=bed_MIN1(N2,l)
                  ZEMAX1(3)=bed_MAX1(N3,l)
                  ZEMIN1(3)=bed_MIN1(N3,l)
               endif
            enddo
               
            
#endif
            

!     COMPUTE THE VERTEX VALUES

            ZEVERTEX(1)=ZEC(1)
            ZEVERTEX(2)=ZEC(1)
            ZEVERTEX(3)=ZEC(1)
            DO KK=2,3
               ZEVERTEX(1)=ZEVERTEX(1)+ dg%PHI_CORNER(KK,1,1)*ZEC(KK)
               ZEVERTEX(2)=ZEVERTEX(2)+ dg%PHI_CORNER(KK,2,1)*ZEC(KK)
               ZEVERTEX(3)=ZEVERTEX(3)+ dg%PHI_CORNER(KK,3,1)*ZEC(KK)
            ENDDO

            
!     RESET THE VERTEX VALUE TO BE LESS THAN OR EQUAL TO THE MAX AND
!     GREATER THAN OR EQUAL TO THE MIN AT THAT VERTEX
!     
            ZEVERTEX(1)=DMAX1(DMIN1(ZEVERTEX(1),ZEMAX1(1)),ZEMIN1(1))
            ZEVERTEX(2)=DMAX1(DMIN1(ZEVERTEX(2),ZEMAX1(2)),ZEMIN1(2))
            ZEVERTEX(3)=DMAX1(DMIN1(ZEVERTEX(3),ZEMAX1(3)),ZEMIN1(3))

            tmp1 = ZEVERTEX(1)
            tmp2 = ZEVERTEX(2)
            tmp3 = ZEVERTEX(3)


!     LOOP OVER THE VERTICES 3 TIMES
!     IF THE VALUE AT THE VERTEX IS ABOVE (BELOW) THE MAX (MIN) AT THAT
!     VERTEX THEN SUBTRACT OFF THE DIFFERENCE AND ADD IT TO THE OTHER
!     VERTICES
!     
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
!     
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

!$$$               if( (abs(tmp1-ZEVERTEX(1)).ge.bound).or.
!$$$     &              (abs(tmp2-ZEVERTEX(2)).ge.bound).or.
!$$$     &              (abs(tmp3-ZEVERTEX(3)).ge.bound).and.
!$$$     &              (dg%slopeflag.eq.5) ) then
!$$$
!$$$                  dg%ZE(4:dg%dofs(i),i,dg%irk+1)=0.D0
!$$$
!$$$                  if (dg%padapt.eq.1) then
!$$$                     pdg_el(i) = 1
!$$$                  endif 
!$$$
!$$$               endif

               dg%ZE(2,I,dg%IRK+1)=-1.d0/6.d0*(ZEVERTEX(1)+ZEVERTEX(2))&
                   +1.d0/3.d0*ZEVERTEX(3)
               dg%ZE(3,I,dg%IRK+1)=-.5d0*ZEVERTEX(1)+.5d0*ZEVERTEX(2)
            ENDIF
            IF (IVAR.EQ.2) THEN

!$$$               if( (abs(tmp1-ZEVERTEX(1)).ge.bound).or.
!$$$     &              (abs(tmp2-ZEVERTEX(2)).ge.bound).or.
!$$$     &              (abs(tmp3-ZEVERTEX(3)).ge.bound).and.
!$$$     &              (dg%slopeflag.eq.5) ) then
!$$$
!$$$
!$$$                  dg%QX(4:dg%dofs(i),i,dg%irk+1)=0.D0
!$$$
!$$$                  if (dg%padapt.eq.1) then
!$$$                     pdg_el(i) = 1
!$$$                  endif 
!$$$
!$$$               endif

               dg%QX(2,I,dg%IRK+1)=-1.d0/6.d0*(ZEVERTEX(1)+ZEVERTEX(2))&
              +1.d0/3.d0*ZEVERTEX(3)
               dg%QX(3,I,dg%IRK+1)=-.5d0*ZEVERTEX(1)+.5d0*ZEVERTEX(2)
            ENDIF
            IF (IVAR.EQ.3) THEN

!$$$               if( (abs(tmp1-ZEVERTEX(1)).ge.bound).or.
!$$$     &              (abs(tmp2-ZEVERTEX(2)).ge.bound).or.
!$$$     &              (abs(tmp3-ZEVERTEX(3)).ge.bound).and.
!$$$     &              (dg%slopeflag.eq.5) ) then
!$$$
!$$$                  dg%QY(4:dg%dofs(i),i,dg%irk+1)=0.D0
!$$$
!$$$                  if (dg%padapt.eq.1) then
!$$$                     pdg_el(i) = 1
!$$$                  endif 
!$$$
!$$$               endif

               dg%QY(2,I,dg%IRK+1)=-1.d0/6.d0*(ZEVERTEX(1)+ZEVERTEX(2))&
                   +1.d0/3.d0*ZEVERTEX(3)
               dg%QY(3,I,dg%IRK+1)=-.5d0*ZEVERTEX(1)+.5d0*ZEVERTEX(2)
            ENDIF

#ifdef TRACE
            IF (IVAR.EQ.4) THEN
               
!$$$               if( (abs(tmp1-ZEVERTEX(1)).ge.bound).or.
!$$$     &              (abs(tmp2-ZEVERTEX(2)).ge.bound).or.
!$$$     &              (abs(tmp3-ZEVERTEX(3)).ge.bound).and.
!$$$     &              (dg%slopeflag.eq.5) ) then
!$$$
!$$$
!$$$                  dg%iota(4:dg%dofs(i),i,dg%irk+1)=0.D0
!$$$
!$$$                  if (dg%padapt.eq.1) then
!$$$                     pdg_el(i) = 1
!$$$                  endif 
!$$$
!$$$               endif
               
               dg%iota(2,I,dg%IRK+1)=-1.d0/6.d0*(ZEVERTEX(1)+ZEVERTEX(2))&
              +1.d0/3.d0*ZEVERTEX(3)
               dg%iota(3,I,dg%IRK+1)=-.5d0*ZEVERTEX(1)+.5d0*ZEVERTEX(2)

            ENDIF
#endif

#ifdef CHEM
            IF (IVAR.EQ.4) THEN
               
!$$$               if( (abs(tmp1-ZEVERTEX(1)).ge.bound).or.
!$$$     &              (abs(tmp2-ZEVERTEX(2)).ge.bound).or.
!$$$     &              (abs(tmp3-ZEVERTEX(3)).ge.bound).and.
!$$$     &              (dg%slopeflag.eq.5) ) then
!$$$
!$$$
!$$$                  dg%iota(4:dg%dofs(i),i,dg%irk+1)=0.D0
!$$$
!$$$                  if (dg%padapt.eq.1) then
!$$$                     pdg_el(i) = 1
!$$$                  endif 
!$$$
!$$$               endif
               
               dg%iota(2,I,dg%IRK+1)=-1.d0/6.d0*(ZEVERTEX(1)+ZEVERTEX(2))&
                   +1.d0/3.d0*ZEVERTEX(3)
               dg%iota(3,I,dg%IRK+1)=-.5d0*ZEVERTEX(1)+.5d0*ZEVERTEX(2)

            ENDIF
            IF (IVAR.EQ.5) THEN

!$$$               if( (abs(tmp1-ZEVERTEX(1)).ge.bound).or.
!$$$     &              (abs(tmp2-ZEVERTEX(2)).ge.bound).or.
!$$$     &              (abs(tmp3-ZEVERTEX(3)).ge.bound).and.
!$$$     &              (dg%slopeflag.eq.5) ) then
!$$$
!$$$
!$$$                  dg%iota2(4:dg%dofs(i),i,dg%irk+1)=0.D0
!$$$
!$$$                  if (dg%padapt.eq.1) then
!$$$                     pdg_el(i) = 1
!$$$                  endif 
!$$$
!$$$               endif

               dg%iota2(2,I,dg%IRK+1)=-1.d0/6.d0*(ZEVERTEX(1)+ZEVERTEX(2))&
                   +1.d0/3.d0*ZEVERTEX(3)
               dg%iota2(3,I,dg%IRK+1)=-.5d0*ZEVERTEX(1)+.5d0*ZEVERTEX(2)
            ENDIF
#endif
         ENDDO

#ifdef SED_LAY

         do l=1,s%layers
            if (IVAR.eq.varnum_prev+l) then

!$$$              if( (abs(tmp1-ZEVERTEX(1)).ge.bound).or.
!$$$     &              (abs(tmp2-ZEVERTEX(2)).ge.bound).or.
!$$$     &              (abs(tmp3-ZEVERTEX(3)).ge.bound).and.
!$$$     &              (dg%slopeflag.eq.5) ) then
!$$$
!$$$
!$$$                  dg%bed(4:dg%dofs(i),i,dg%irk+1,l)=0.D0
!$$$
!$$$                  if (dg%padapt.eq.1) then
!$$$                     pdg_el(i) = 1
!$$$                  endif 
!$$$
!$$$               endif
               
               dg%bed(2,I,dg%IRK+1,l)=-1.d0/6.d0*(ZEVERTEX(1)+ZEVERTEX(2))&
                   +1.d0/3.d0*ZEVERTEX(3)
               dg%bed(3,I,dg%IRK+1,l)=-.5d0*ZEVERTEX(1)+.5d0*ZEVERTEX(2)
            
            
            endif
         enddo

#endif

                                !endif

      ENDDO
      RETURN
      END SUBROUTINE 

#endif

#ifdef STBLZR

!***********************************************************************
!     
!     SUBROUTINE SLOPELIMITER5()
!     
!     Written by Clint Dawson - 30 June 2010
!     01-10-2011 - cem - adapted for p_enrichment and multicomponent
!     
!     
!***********************************************************************

      SUBROUTINE SLOPELIMITER5(s)

!.....Use appropriate modules

      USE SIZES
      USE GLOBAL
      USE DG

#ifdef CMPI
      USE MESSENGER
#endif

      IMPLICIT NONE
      
      type (sizes_type) :: s

!.....Declare local variables

      INTEGER L, LL, INC1,INC2,INC3,KDP,NN,IVAR,I,J,kk,k,varnum,bb,varnum_prev
      REAL(SZ) ZEC(3),ZEVERTEX(3),DIF(3),SUMLOC,SUMDIF,SIGNDIF
      REAL(SZ)     DIV,REDFAC,REDMAX,tmp1,tmp2,tmp3,bound
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
      Allocate ( bed_MIN1(NP,s%layers),bed_MAX1(NP,s%layers) )

!     FIND THE MAXIMUM AND MINIMUM OF EACH VARIABLE OVER ALL ELEMENTS 
!     SHARING A NODE

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
         do l=1,s%layers
            bed_MIN1(I,l) = 99999.
            bed_MAX1(I,l) =-99999.
         enddo
#endif

         NO_NBORS = EL_COUNT(I)

         DO J = 1,NO_NBORS
            NBOR_EL = ELETAB(I,1+J)

!     IF(dg%WDFLG(NBOR_EL).EQ.0) CYCLE ! DON'T COUNT DRY ELEMENTS  sb 02/26/07

            ZE_DG(J) = dg%ZE(1,NBOR_EL,dg%IRK+1)
            QX_DG(J) = dg%QX(1,NBOR_EL,dg%IRK+1)
            QY_DG(J) = dg%QY(1,NBOR_EL,dg%IRK+1)

#ifdef TRACE
            iota_DG(J) = dg%iota(1,NBOR_EL,dg%IRK+1)
#endif

#ifdef CHEM
            iota_DG(J) = dg%iota(1,NBOR_EL,dg%IRK+1)
            iota2_DG(J) = dg%iota2(1,NBOR_EL,dg%IRK+1)
#endif


#ifdef SED_LAY
            do l=1,s%layers
               bed_DG(J,l) = dg%bed(1,NBOR_EL,dg%IRK+1,l)
            enddo
#endif

!     
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
            do l=1,s%layers
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
         do l = 1,s%layers
            arraymax => bed_max1(:,l)
            arraymin => bed_min1(:,l)
            CALL UPDATER(arraymax,arraymin,QY_MAX1,2 )
         enddo
#endif

#endif
!     
!     LOOP OVER ELEMENTS TO CALCULATE NEW VERTEX VALUES
!     

      bb = 1

      DO I=1,NE 
         !IF(dg%WDFLG(I).EQ.0) CYCLE ! DON'T COUNT DRY ELEMENTS  sb 02/26/07
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
         varnum = varnum_prev + s%layers
#endif

         DO IVAR=1,varnum

            IF (IVAR.EQ.1) THEN
               ZEC(1)=dg%ZE(1,I,dg%IRK+1)
               ZEC(2)=dg%ZE(2,I,dg%IRK+1)
               ZEC(3)=dg%ZE(3,I,dg%IRK+1)
               ZEMAX1(1)=ZE_MAX1(N1)
               ZEMIN1(1)=ZE_MIN1(N1)
               ZEMAX1(2)=ZE_MAX1(N2)
               ZEMIN1(2)=ZE_MIN1(N2)
               ZEMAX1(3)=ZE_MAX1(N3)
               ZEMIN1(3)=ZE_MIN1(N3)
            ENDIF

            IF (IVAR.EQ.2) THEN
               ZEC(1)=dg%QX(1,I,dg%IRK+1)
               ZEC(2)=dg%QX(2,I,dg%IRK+1)
               ZEC(3)=dg%QX(3,I,dg%IRK+1)
               ZEMAX1(1)=QX_MAX1(N1)
               ZEMIN1(1)=QX_MIN1(N1)
               ZEMAX1(2)=QX_MAX1(N2)
               ZEMIN1(2)=QX_MIN1(N2)
               ZEMAX1(3)=QX_MAX1(N3)
               ZEMIN1(3)=QX_MIN1(N3)
            ENDIF

            IF (IVAR.EQ.3) THEN
               ZEC(1)=dg%QY(1,I,dg%IRK+1)
               ZEC(2)=dg%QY(2,I,dg%IRK+1)
               ZEC(3)=dg%QY(3,I,dg%IRK+1)
               ZEMAX1(1)=QY_MAX1(N1)
               ZEMIN1(1)=QY_MIN1(N1)
               ZEMAX1(2)=QY_MAX1(N2)
               ZEMIN1(2)=QY_MIN1(N2)
               ZEMAX1(3)=QY_MAX1(N3)
               ZEMIN1(3)=QY_MIN1(N3)
            ENDIF

#ifdef TRACE
            IF (IVAR.EQ.4) THEN
               ZEC(1)=dg%iota(1,I,dg%IRK+1)
               ZEC(2)=dg%iota(2,I,dg%IRK+1)
               ZEC(3)=dg%iota(3,I,dg%IRK+1)
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
               ZEC(1)=dg%iota(1,I,dg%IRK+1)
               ZEC(2)=dg%iota(2,I,dg%IRK+1)
               ZEC(3)=dg%iota(3,I,dg%IRK+1)
               ZEMAX1(1)=iota_MAX1(N1)
               ZEMIN1(1)=iota_MIN1(N1)
               ZEMAX1(2)=iota_MAX1(N2)
               ZEMIN1(2)=iota_MIN1(N2)
               ZEMAX1(3)=iota_MAX1(N3)
               ZEMIN1(3)=iota_MIN1(N3)
            ENDIF

            IF (IVAR.EQ.5) THEN
               ZEC(1)=dg%iota2(1,I,dg%IRK+1)
               ZEC(2)=dg%iota2(2,I,dg%IRK+1)
               ZEC(3)=dg%iota2(3,I,dg%IRK+1)
               ZEMAX1(1)=iota2_MAX1(N1)
               ZEMIN1(1)=iota2_MIN1(N1)
               ZEMAX1(2)=iota2_MAX1(N2)
               ZEMIN1(2)=iota2_MIN1(N2)
               ZEMAX1(3)=iota2_MAX1(N3)
               ZEMIN1(3)=iota2_MIN1(N3)
            ENDIF
#endif

#ifdef SED_LAY
            do l=1,s%layers
               if (IVAR.eq.varnum_prev+l) then
                  ZEC(1)=dg%bed(1,I,dg%IRK+1,l)
                  ZEC(2)=dg%bed(2,I,dg%IRK+1,l)
                  ZEC(3)=dg%bed(3,I,dg%IRK+1,l)
                  ZEMAX1(1)=bed_MAX1(N1,l)
                  ZEMIN1(1)=bed_MIN1(N1,l)
                  ZEMAX1(2)=bed_MAX1(N2,l)
                  ZEMIN1(2)=bed_MIN1(N2,l)
                  ZEMAX1(3)=bed_MAX1(N3,l)
                  ZEMIN1(3)=bed_MIN1(N3,l)
               endif
            enddo
               
            
#endif
            

!     COMPUTE THE VERTEX VALUES

            ZEVERTEX(1)=ZEC(1)
            ZEVERTEX(2)=ZEC(1)
            ZEVERTEX(3)=ZEC(1)
            DO KK=2,3
               ZEVERTEX(1)=ZEVERTEX(1)+ dg%PHI_CORNER(KK,1,1)*ZEC(KK)
               ZEVERTEX(2)=ZEVERTEX(2)+ dg%PHI_CORNER(KK,2,1)*ZEC(KK)
               ZEVERTEX(3)=ZEVERTEX(3)+ dg%PHI_CORNER(KK,3,1)*ZEC(KK)
            ENDDO

            
!     RESET THE VERTEX VALUE TO BE LESS THAN OR EQUAL TO THE MAX AND
!     GREATER THAN OR EQUAL TO THE MIN AT THAT VERTEX
!     
            ZEVERTEX(1)=DMAX1(DMIN1(ZEVERTEX(1),ZEMAX1(1)),ZEMIN1(1))
            ZEVERTEX(2)=DMAX1(DMIN1(ZEVERTEX(2),ZEMAX1(2)),ZEMIN1(2))
            ZEVERTEX(3)=DMAX1(DMIN1(ZEVERTEX(3),ZEMAX1(3)),ZEMIN1(3))

            tmp1 = ZEVERTEX(1)
            tmp2 = ZEVERTEX(2)
            tmp3 = ZEVERTEX(3)


!     LOOP OVER THE VERTICES 3 TIMES
!     IF THE VALUE AT THE VERTEX IS ABOVE (BELOW) THE MAX (MIN) AT THAT
!     VERTEX THEN SUBTRACT OFF THE DIFFERENCE AND ADD IT TO THE OTHER
!     VERTICES
!     
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
!     
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

!$$$               if( (abs(tmp1-ZEVERTEX(1)).ge.bound).or.
!$$$     &              (abs(tmp2-ZEVERTEX(2)).ge.bound).or.
!$$$     &              (abs(tmp3-ZEVERTEX(3)).ge.bound).and.
!$$$     &              (dg%slopeflag.eq.5) ) then
!$$$
!$$$                  dg%ZE(4:dg%dofs(i),i,dg%irk+1)=0.D0
!$$$
!$$$                  if (dg%padapt.eq.1) then
!$$$                     pdg_el(i) = 1
!$$$                  endif 
!$$$
!$$$               endif

               dg%ZE(2,I,dg%IRK+1)=-1.d0/6.d0*(ZEVERTEX(1)+ZEVERTEX(2))&
                   +1.d0/3.d0*ZEVERTEX(3)
               dg%ZE(3,I,dg%IRK+1)=-.5d0*ZEVERTEX(1)+.5d0*ZEVERTEX(2)
            ENDIF
            IF (IVAR.EQ.2) THEN

!$$$               if( (abs(tmp1-ZEVERTEX(1)).ge.bound).or.
!$$$     &              (abs(tmp2-ZEVERTEX(2)).ge.bound).or.
!$$$     &              (abs(tmp3-ZEVERTEX(3)).ge.bound).and.
!$$$     &              (dg%slopeflag.eq.5) ) then
!$$$
!$$$
!$$$                  dg%QX(4:dg%dofs(i),i,dg%irk+1)=0.D0
!$$$
!$$$                  if (dg%padapt.eq.1) then
!$$$                     pdg_el(i) = 1
!$$$                  endif 
!$$$
!$$$               endif

               dg%QX(2,I,dg%IRK+1)=-1.d0/6.d0*(ZEVERTEX(1)+ZEVERTEX(2))&
                   +1.d0/3.d0*ZEVERTEX(3)
               dg%QX(3,I,dg%IRK+1)=-.5d0*ZEVERTEX(1)+.5d0*ZEVERTEX(2)
            ENDIF
            IF (IVAR.EQ.3) THEN

!$$$               if( (abs(tmp1-ZEVERTEX(1)).ge.bound).or.
!$$$     &              (abs(tmp2-ZEVERTEX(2)).ge.bound).or.
!$$$     &              (abs(tmp3-ZEVERTEX(3)).ge.bound).and.
!$$$     &              (dg%slopeflag.eq.5) ) then
!$$$
!$$$                  dg%QY(4:dg%dofs(i),i,dg%irk+1)=0.D0
!$$$
!$$$                  if (dg%padapt.eq.1) then
!$$$                     pdg_el(i) = 1
!$$$                  endif 
!$$$
!$$$               endif

               dg%QY(2,I,dg%IRK+1)=-1.d0/6.d0*(ZEVERTEX(1)+ZEVERTEX(2))&
              +1.d0/3.d0*ZEVERTEX(3)
               dg%QY(3,I,dg%IRK+1)=-.5d0*ZEVERTEX(1)+.5d0*ZEVERTEX(2)
            ENDIF

#ifdef TRACE
            IF (IVAR.EQ.4) THEN
               
!$$$               if( (abs(tmp1-ZEVERTEX(1)).ge.bound).or.
!$$$     &              (abs(tmp2-ZEVERTEX(2)).ge.bound).or.
!$$$     &              (abs(tmp3-ZEVERTEX(3)).ge.bound).and.
!$$$     &              (dg%slopeflag.eq.5) ) then
!$$$
!$$$
!$$$                  dg%iota(4:dg%dofs(i),i,dg%irk+1)=0.D0
!$$$
!$$$                  if (dg%padapt.eq.1) then
!$$$                     pdg_el(i) = 1
!$$$                  endif 
!$$$
!$$$               endif
               
               dg%iota(2,I,dg%IRK+1)=-1.d0/6.d0*(ZEVERTEX(1)+ZEVERTEX(2))&
              +1.d0/3.d0*ZEVERTEX(3)
               dg%iota(3,I,dg%IRK+1)=-.5d0*ZEVERTEX(1)+.5d0*ZEVERTEX(2)

            ENDIF
#endif

#ifdef CHEM
            IF (IVAR.EQ.4) THEN
               
!$$$               if( (abs(tmp1-ZEVERTEX(1)).ge.bound).or.
!$$$     &              (abs(tmp2-ZEVERTEX(2)).ge.bound).or.
!$$$     &              (abs(tmp3-ZEVERTEX(3)).ge.bound).and.
!$$$     &              (dg%slopeflag.eq.5) ) then
!$$$
!$$$
!$$$                  dg%iota(4:dg%dofs(i),i,dg%irk+1)=0.D0
!$$$
!$$$                  if (dg%padapt.eq.1) then
!$$$                     pdg_el(i) = 1
!$$$                  endif 
!$$$
!$$$               endif
               
               dg%iota(2,I,dg%IRK+1)=-1.d0/6.d0*(ZEVERTEX(1)+ZEVERTEX(2))&
              +1.d0/3.d0*ZEVERTEX(3)
               dg%iota(3,I,dg%IRK+1)=-.5d0*ZEVERTEX(1)+.5d0*ZEVERTEX(2)

            ENDIF
            IF (IVAR.EQ.5) THEN

!$$$               if( (abs(tmp1-ZEVERTEX(1)).ge.bound).or.
!$$$     &              (abs(tmp2-ZEVERTEX(2)).ge.bound).or.
!$$$     &              (abs(tmp3-ZEVERTEX(3)).ge.bound).and.
!$$$     &              (dg%slopeflag.eq.5) ) then
!$$$
!$$$
!$$$                  dg%iota2(4:dg%dofs(i),i,dg%irk+1)=0.D0
!$$$
!$$$                  if (dg%padapt.eq.1) then
!$$$                     pdg_el(i) = 1
!$$$                  endif 
!$$$
!$$$               endif

               dg%iota2(2,I,dg%IRK+1)=-1.d0/6.d0*(ZEVERTEX(1)+ZEVERTEX(2))&
              +1.d0/3.d0*ZEVERTEX(3)
               dg%iota2(3,I,dg%IRK+1)=-.5d0*ZEVERTEX(1)+.5d0*ZEVERTEX(2)
            ENDIF
#endif
         ENDDO

#ifdef SED_LAY

         do l=1,s%layers
            if (IVAR.eq.varnum_prev+l) then

!$$$              if( (abs(tmp1-ZEVERTEX(1)).ge.bound).or.
!$$$     &              (abs(tmp2-ZEVERTEX(2)).ge.bound).or.
!$$$     &              (abs(tmp3-ZEVERTEX(3)).ge.bound).and.
!$$$     &              (dg%slopeflag.eq.5) ) then
!$$$
!$$$
!$$$                  dg%bed(4:dg%dofs(i),i,dg%irk+1,l)=0.D0
!$$$
!$$$                  if (dg%padapt.eq.1) then
!$$$                     pdg_el(i) = 1
!$$$                  endif 
!$$$
!$$$               endif
               
               dg%bed(2,I,dg%IRK+1,l)=-1.d0/6.d0*(ZEVERTEX(1)+ZEVERTEX(2))&
              +1.d0/3.d0*ZEVERTEX(3)
               dg%bed(3,I,dg%IRK+1,l)=-.5d0*ZEVERTEX(1)+.5d0*ZEVERTEX(2)
            
            
            endif
         enddo

#endif

                                !endif

      ENDDO
      RETURN
      END SUBROUTINE 

#endif


#ifdef SLOPEALL

!***********************************************************************
!     
!     SUBROUTINE SLOPELIMITER5()
!     
!     Written by Clint Dawson - 30 June 2010
!     01-10-2011 - cem - adapted for p_enrichment and multicomponent
!     
!     
!***********************************************************************

      SUBROUTINE SLOPELIMITER5()

!.....Use appropriate modules

      USE SIZES
      USE GLOBAL
      USE DG

#ifdef CMPI
      USE MESSENGER
#endif

      IMPLICIT NONE

      type (sizes_type) :: s

!.....Declare local variables

      INTEGER L, LL, INC1,INC2,INC3,KDP,NN,IVAR,I,J,kk,k,varnum,bb,varnum_prev
      REAL(SZ) ZEC(3),ZEVERTEX(3),DIF(3),SUMLOC,SUMDIF,SIGNDIF
      REAL(SZ) DIV,REDFAC,REDMAX,tmp1,tmp2,tmp3,bound
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
      Allocate ( bed_MIN1(NP,s%layers),bed_MAX1(NP,s%layers) )

!     FIND THE MAXIMUM AND MINIMUM OF EACH VARIABLE OVER ALL ELEMENTS 
!     SHARING A NODE

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
         do l=1,s%layers
            bed_MIN1(I,l) = 99999.
            bed_MAX1(I,l) =-99999.
         enddo
#endif

         NO_NBORS = EL_COUNT(I)

         DO J = 1,NO_NBORS
            NBOR_EL = ELETAB(I,1+J)

!     IF(dg%WDFLG(NBOR_EL).EQ.0) CYCLE ! DON'T COUNT DRY ELEMENTS  sb 02/26/07

            ZE_DG(J) = dg%ZE(1,NBOR_EL,dg%IRK+1)
            QX_DG(J) = dg%QX(1,NBOR_EL,dg%IRK+1)
            QY_DG(J) = dg%QY(1,NBOR_EL,dg%IRK+1)

#ifdef TRACE
            iota_DG(J) = dg%iota(1,NBOR_EL,dg%IRK+1)
#endif

#ifdef CHEM
            iota_DG(J) = dg%iota(1,NBOR_EL,dg%IRK+1)
            iota2_DG(J) = dg%iota2(1,NBOR_EL,dg%IRK+1)
#endif


#ifdef SED_LAY
            do l=1,s%layers
               bed_DG(J,l) = dg%bed(1,NBOR_EL,dg%IRK+1,l)
            enddo
#endif

!     
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
            do l=1,s%layers
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
         do l = 1,s%layers
            arraymax => bed_max1(:,l)
            arraymin => bed_min1(:,l)
            CALL UPDATER(arraymax,arraymin,QY_MAX1,2 )
         enddo
#endif

#endif
!     
!     LOOP OVER ELEMENTS TO CALCULATE NEW VERTEX VALUES
!     

      bb = 1

      DO I=1,NE 

                                !if (dg%dofs(i).eq.3) then

         !IF(dg%WDFLG(I).EQ.0) CYCLE ! DON'T COUNT DRY ELEMENTS  sb 02/26/07
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
         varnum = varnum_prev + s%layers
#endif

         DO IVAR=1,varnum

            IF (IVAR.EQ.1) THEN
               ZEC(1)=dg%ZE(1,I,dg%IRK+1)
               ZEC(2)=dg%ZE(2,I,dg%IRK+1)
               ZEC(3)=dg%ZE(3,I,dg%IRK+1)
               ZEMAX1(1)=ZE_MAX1(N1)
               ZEMIN1(1)=ZE_MIN1(N1)
               ZEMAX1(2)=ZE_MAX1(N2)
               ZEMIN1(2)=ZE_MIN1(N2)
               ZEMAX1(3)=ZE_MAX1(N3)
               ZEMIN1(3)=ZE_MIN1(N3)
            ENDIF

            IF (IVAR.EQ.2) THEN
               ZEC(1)=dg%QX(1,I,dg%IRK+1)
               ZEC(2)=dg%QX(2,I,dg%IRK+1)
               ZEC(3)=dg%QX(3,I,dg%IRK+1)
               ZEMAX1(1)=QX_MAX1(N1)
               ZEMIN1(1)=QX_MIN1(N1)
               ZEMAX1(2)=QX_MAX1(N2)
               ZEMIN1(2)=QX_MIN1(N2)
               ZEMAX1(3)=QX_MAX1(N3)
               ZEMIN1(3)=QX_MIN1(N3)
            ENDIF

            IF (IVAR.EQ.3) THEN
               ZEC(1)=dg%QY(1,I,dg%IRK+1)
               ZEC(2)=dg%QY(2,I,dg%IRK+1)
               ZEC(3)=dg%QY(3,I,dg%IRK+1)
               ZEMAX1(1)=QY_MAX1(N1)
               ZEMIN1(1)=QY_MIN1(N1)
               ZEMAX1(2)=QY_MAX1(N2)
               ZEMIN1(2)=QY_MIN1(N2)
               ZEMAX1(3)=QY_MAX1(N3)
               ZEMIN1(3)=QY_MIN1(N3)
            ENDIF

#ifdef TRACE
            IF (IVAR.EQ.4) THEN
               ZEC(1)=dg%iota(1,I,dg%IRK+1)
               ZEC(2)=dg%iota(2,I,dg%IRK+1)
               ZEC(3)=dg%iota(3,I,dg%IRK+1)
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
               ZEC(1)=dg%iota(1,I,dg%IRK+1)
               ZEC(2)=dg%iota(2,I,dg%IRK+1)
               ZEC(3)=dg%iota(3,I,dg%IRK+1)
               ZEMAX1(1)=iota_MAX1(N1)
               ZEMIN1(1)=iota_MIN1(N1)
               ZEMAX1(2)=iota_MAX1(N2)
               ZEMIN1(2)=iota_MIN1(N2)
               ZEMAX1(3)=iota_MAX1(N3)
               ZEMIN1(3)=iota_MIN1(N3)
            ENDIF

            IF (IVAR.EQ.5) THEN
               ZEC(1)=dg%iota2(1,I,dg%IRK+1)
               ZEC(2)=dg%iota2(2,I,dg%IRK+1)
               ZEC(3)=dg%iota2(3,I,dg%IRK+1)
               ZEMAX1(1)=iota2_MAX1(N1)
               ZEMIN1(1)=iota2_MIN1(N1)
               ZEMAX1(2)=iota2_MAX1(N2)
               ZEMIN1(2)=iota2_MIN1(N2)
               ZEMAX1(3)=iota2_MAX1(N3)
               ZEMIN1(3)=iota2_MIN1(N3)
            ENDIF
#endif

#ifdef SED_LAY
            do l=1,s%layers
               if (IVAR.eq.varnum_prev+l) then
                  ZEC(1)=dg%bed(1,I,dg%IRK+1,l)
                  ZEC(2)=dg%bed(2,I,dg%IRK+1,l)
                  ZEC(3)=dg%bed(3,I,dg%IRK+1,l)
                  ZEMAX1(1)=bed_MAX1(N1,l)
                  ZEMIN1(1)=bed_MIN1(N1,l)
                  ZEMAX1(2)=bed_MAX1(N2,l)
                  ZEMIN1(2)=bed_MIN1(N2,l)
                  ZEMAX1(3)=bed_MAX1(N3,l)
                  ZEMIN1(3)=bed_MIN1(N3,l)
               endif
            enddo
               
            
#endif
            

!     COMPUTE THE VERTEX VALUES

            ZEVERTEX(1)=ZEC(1)
            ZEVERTEX(2)=ZEC(1)
            ZEVERTEX(3)=ZEC(1)
            DO KK=2,3
               ZEVERTEX(1)=ZEVERTEX(1)+ dg%PHI_CORNER(KK,1,1)*ZEC(KK)
               ZEVERTEX(2)=ZEVERTEX(2)+ dg%PHI_CORNER(KK,2,1)*ZEC(KK)
               ZEVERTEX(3)=ZEVERTEX(3)+ dg%PHI_CORNER(KK,3,1)*ZEC(KK)
            ENDDO

            
!     RESET THE VERTEX VALUE TO BE LESS THAN OR EQUAL TO THE MAX AND
!     GREATER THAN OR EQUAL TO THE MIN AT THAT VERTEX
!     
            ZEVERTEX(1)=DMAX1(DMIN1(ZEVERTEX(1),ZEMAX1(1)),ZEMIN1(1))
            ZEVERTEX(2)=DMAX1(DMIN1(ZEVERTEX(2),ZEMAX1(2)),ZEMIN1(2))
            ZEVERTEX(3)=DMAX1(DMIN1(ZEVERTEX(3),ZEMAX1(3)),ZEMIN1(3))

            tmp1 = ZEVERTEX(1)
            tmp2 = ZEVERTEX(2)
            tmp3 = ZEVERTEX(3)


!     LOOP OVER THE VERTICES 3 TIMES
!     IF THE VALUE AT THE VERTEX IS ABOVE (BELOW) THE MAX (MIN) AT THAT
!     VERTEX THEN SUBTRACT OFF THE DIFFERENCE AND ADD IT TO THE OTHER
!     VERTICES
!     
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
!     
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

!$$$               if( (abs(tmp1-ZEVERTEX(1)).gt.bound).or.
!$$$     &              (abs(tmp2-ZEVERTEX(2)).gt.bound).or.
!$$$     &              (abs(tmp3-ZEVERTEX(3)).gt.bound).and.
!$$$     &              (dg%slopeflag.eq.5) ) then
!$$$
!$$$                  dg%ZE(4:dg%dofs(i),i,dg%irk+1)=0.D0
!$$$
!$$$                  if (dg%padapt.eq.1) then
!$$$                     pdg_el(i) = 1
!$$$                  endif 
!$$$
!$$$               endif

               dg%ZE(2,I,dg%IRK+1)=-1.d0/6.d0*(ZEVERTEX(1)+ZEVERTEX(2))&
              +1.d0/3.d0*ZEVERTEX(3)
               dg%ZE(3,I,dg%IRK+1)=-.5d0*ZEVERTEX(1)+.5d0*ZEVERTEX(2)
            ENDIF
            IF (IVAR.EQ.2) THEN

!$$$               if( (abs(tmp1-ZEVERTEX(1)).gt.bound).or.
!$$$     &              (abs(tmp2-ZEVERTEX(2)).gt.bound).or.
!$$$     &              (abs(tmp3-ZEVERTEX(3)).gt.bound).and.
!$$$     &              (dg%slopeflag.eq.5) ) then
!$$$
!$$$
!$$$                  dg%QX(4:dg%dofs(i),i,dg%irk+1)=0.D0
!$$$
!$$$                  if (dg%padapt.eq.1) then
!$$$                     pdg_el(i) = 1
!$$$                  endif 
!$$$
!$$$               endif

               dg%QX(2,I,dg%IRK+1)=-1.d0/6.d0*(ZEVERTEX(1)+ZEVERTEX(2))&
              +1.d0/3.d0*ZEVERTEX(3)
               dg%QX(3,I,dg%IRK+1)=-.5d0*ZEVERTEX(1)+.5d0*ZEVERTEX(2)
            ENDIF
            IF (IVAR.EQ.3) THEN

!$$$               if( (abs(tmp1-ZEVERTEX(1)).gt.bound).or.
!$$$     &              (abs(tmp2-ZEVERTEX(2)).gt.bound).or.
!$$$     &              (abs(tmp3-ZEVERTEX(3)).gt.bound).and.
!$$$     &              (dg%slopeflag.eq.5) ) then
!$$$
!$$$                  dg%QY(4:dg%dofs(i),i,dg%irk+1)=0.D0
!$$$
!$$$                  if (dg%padapt.eq.1) then
!$$$                     pdg_el(i) = 1
!$$$                  endif 
!$$$
!$$$               endif

               dg%QY(2,I,dg%IRK+1)=-1.d0/6.d0*(ZEVERTEX(1)+ZEVERTEX(2))&
              +1.d0/3.d0*ZEVERTEX(3)
               dg%QY(3,I,dg%IRK+1)=-.5d0*ZEVERTEX(1)+.5d0*ZEVERTEX(2)
            ENDIF

#ifdef TRACE
            IF (IVAR.EQ.4) THEN
               
!$$$               if( (abs(tmp1-ZEVERTEX(1)).gt.bound).or.
!$$$     &              (abs(tmp2-ZEVERTEX(2)).gt.bound).or.
!$$$     &              (abs(tmp3-ZEVERTEX(3)).gt.bound).and.
!$$$     &              (dg%slopeflag.eq.5) ) then
!$$$
!$$$
!$$$                  dg%iota(4:dg%dofs(i),i,dg%irk+1)=0.D0
!$$$
!$$$                  if (dg%padapt.eq.1) then
!$$$                     pdg_el(i) = 1
!$$$                  endif 
!$$$
!$$$               endif
               
               dg%iota(2,I,dg%IRK+1)=-1.d0/6.d0*(ZEVERTEX(1)+ZEVERTEX(2))&
              +1.d0/3.d0*ZEVERTEX(3)
               dg%iota(3,I,dg%IRK+1)=-.5d0*ZEVERTEX(1)+.5d0*ZEVERTEX(2)

            ENDIF
#endif

#ifdef CHEM
            IF (IVAR.EQ.4) THEN
               
!$$$               if( (abs(tmp1-ZEVERTEX(1)).gt.bound).or.
!$$$     &              (abs(tmp2-ZEVERTEX(2)).gt.bound).or.
!$$$     &              (abs(tmp3-ZEVERTEX(3)).gt.bound).and.
!$$$     &              (dg%slopeflag.eq.5) ) then
!$$$
!$$$
!$$$                  dg%iota(4:dg%dofs(i),i,dg%irk+1)=0.D0
!$$$
!$$$                  if (dg%padapt.eq.1) then
!$$$                     pdg_el(i) = 1
!$$$                  endif 
!$$$
!$$$               endif
               
               dg%iota(2,I,dg%IRK+1)=-1.d0/6.d0*(ZEVERTEX(1)+ZEVERTEX(2))&
              +1.d0/3.d0*ZEVERTEX(3)
               dg%iota(3,I,dg%IRK+1)=-.5d0*ZEVERTEX(1)+.5d0*ZEVERTEX(2)

            ENDIF
            IF (IVAR.EQ.5) THEN

!$$$               if( (abs(tmp1-ZEVERTEX(1)).gt.bound).or.
!$$$     &              (abs(tmp2-ZEVERTEX(2)).gt.bound).or.
!$$$     &              (abs(tmp3-ZEVERTEX(3)).gt.bound).and.
!$$$     &              (dg%slopeflag.eq.5) ) then
!$$$
!$$$
!$$$                  dg%iota2(4:dg%dofs(i),i,dg%irk+1)=0.D0
!$$$
!$$$                  if (dg%padapt.eq.1) then
!$$$                     pdg_el(i) = 1
!$$$                  endif 
!$$$
!$$$               endif

               dg%iota2(2,I,dg%IRK+1)=-1.d0/6.d0*(ZEVERTEX(1)+ZEVERTEX(2))&
              +1.d0/3.d0*ZEVERTEX(3)
               dg%iota2(3,I,dg%IRK+1)=-.5d0*ZEVERTEX(1)+.5d0*ZEVERTEX(2)
            ENDIF
#endif
         ENDDO

#ifdef SED_LAY

         do l=1,s%layers
            if (IVAR.eq.varnum_prev+l) then

!$$$              if( (abs(tmp1-ZEVERTEX(1)).gt.bound).or.
!$$$     &              (abs(tmp2-ZEVERTEX(2)).gt.bound).or.
!$$$     &              (abs(tmp3-ZEVERTEX(3)).gt.bound).and.
!$$$     &              (dg%slopeflag.eq.5) ) then
!$$$
!$$$
!$$$                  dg%bed(4:dg%dofs(i),i,dg%irk+1,l)=0.D0
!$$$
!$$$                  if (dg%padapt.eq.1) then
!$$$                     pdg_el(i) = 1
!$$$                  endif 
!$$$
!$$$               endif
               
               dg%bed(2,I,dg%IRK+1,l)=-1.d0/6.d0*(ZEVERTEX(1)+ZEVERTEX(2))&
              +1.d0/3.d0*ZEVERTEX(3)
               dg%bed(3,I,dg%IRK+1,l)=-.5d0*ZEVERTEX(1)+.5d0*ZEVERTEX(2)
            
            
            endif
         enddo

#endif

                                !endif

      ENDDO
      RETURN
      END SUBROUTINE 

!***********************************************************************
!     
!     SUBROUTINE SLOPELIMITER6()
!     
!     Written 2011
!     
!     This subroutine selects the Barth--Jespersen slope limiter based on
!     a Taylor Polynomial basis, and is consistent with p-adaptation
!     to arbitrary order p
!     
!     -- cem
!     
!***********************************************************************

      SUBROUTINE SLOPELIMITER6()

!.....Use appropriate modules

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

      Allocate ( ZEminel(mne,dg%dofh),ZEmaxel(mne,dg%dofh),QXminel(mne,dg%dofh) )
      Allocate ( QYminel(mne,dg%dofh),QYmaxel(mne,dg%dofh),QXmaxel(mne,dg%dofh) )
      Allocate ( iotaminel(mne,dg%dofh),iotamaxel(mne,dg%dofh) )
      Allocate ( iota2minel(mne,dg%dofh),iota2maxel(mne,dg%dofh) )


!.....We work over the master element
!.....Set initial values

      fd = dg%slope_weight         ! add weight for lower order pieces (fd<1 => stronger limiting)
      

      DO k=1,NE

         if (dg%dofs(k).gt.1) then

            DO ll = 1,dg%dofs(k)

               dg%ZE(ll,k,dg%NRK+2) = dg%ZE(ll,k,dg%IRK+1)
               dg%QX(ll,k,dg%NRK+2) = dg%QX(ll,k,dg%IRK+1)
               dg%QY(ll,k,dg%NRK+2) = dg%QY(ll,k,dg%IRK+1)

#ifdef TRACE
               dg%iota(ll,k,dg%NRK+2) = dg%iota(ll,k,dg%IRK+1)
#endif

#ifdef CHEM
               dg%iota(ll,k,dg%NRK+2) = dg%iota(ll,k,dg%IRK+1)
               dg%iota2(ll,k,dg%NRK+2) = dg%iota2(ll,k,dg%IRK+1)
#endif


            ENDDO

         elseif (dg%dofs(k).eq.1) then

            cycle

         endif

      ENDDO   

      marea = 2.D0              !master elements area

!.....Convert initial values to the Taylor basis (multiply by dg%Nmatrix) on base element


      dg%ZEtaylor = 0.D0 
      dg%QXtaylor = 0.D0
      dg%QYtaylor = 0.D0

#ifdef TRACE
      dg%iotataylor = 0.D0
#endif

#ifdef CHEM
      dg%iotataylor = 0.D0
      dg%iota2taylor = 0.D0
#endif

      do k=1,ne

         if (dg%dofs(k).gt.1) then

            do ll=1,dg%dofs(k)
               
               do ss=1,dg%dofs(k)

                  dg%ZEtaylor(k,ll,1) = dg%ZEtaylor(k,ll,1) + dg%Nmatrix(k,ll,ss,dg%dofs(k)) * dg%ZE(ss,k,dg%nrk+2)
                  dg%QXtaylor(k,ll,1) = dg%QXtaylor(k,ll,1) + dg%Nmatrix(k,ll,ss,dg%dofs(k)) * dg%QX(ss,k,dg%nrk+2)
                  dg%QYtaylor(k,ll,1) = dg%QYtaylor(k,ll,1) + dg%Nmatrix(k,ll,ss,dg%dofs(k)) * dg%QY(ss,k,dg%nrk+2)

#ifdef TRACE
                  dg%iotataylor(k,ll,1) = dg%iotataylor(k,ll,1) + dg%Nmatrix(k,ll,ss,dg%dofs(k))*dg%iota(ss,k,dg%nrk+2)
#endif                     
                  
#ifdef CHEM
                  dg%iotataylor(k,ll,1) = dg%iotataylor(k,ll,1) + dg%Nmatrix(k,ll,ss,dg%dofs(k))*dg%iota(ss,k,dg%nrk+2)
                  dg%iota2taylor(k,ll,1) = dg%iota2taylor(k,ll,1) + dg%Nmatrix(k,ll,ss,dg%dofs(k))*dg%iota2(ss,k,dg%nrk+2)
#endif
                  
               enddo

            enddo

         elseif (dg%dofs(k).eq.1) then

            cycle

         endif

      enddo

!.....Find values at vertices of base elements and neighbors


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

      CALL UPDATER_ELEM_MOD2(dg%ZEtaylor,dg%QXtaylor,dg%QYtaylor,1,3)

#ifdef TRACE
      CALL UPDATER_ELEM_MOD2(dg%iotataylor,dg%QXtaylor,dg%QYtaylor,1,2 )  
#endif

#ifdef CHEM
      CALL UPDATER_ELEM_MOD2(dg%iotataylor,dg%iota2taylor,dg%QYtaylor,1,2 )  
#endif

#endif

      do k = 1,ne

         do ll = 1,dg%dofs(k) 

            do ell = 1,3        ! Number of edge neighbors for a triangle

!.....Find max and min values over polynomial coefficients
               
               ZEmaxel(k,ll) = max( dg%ZEtaylor(k,ll,1),dg%ZEtaylor(dg%EL_NBORS(ell,k),ll,1), ZEmaxel(k,ll) )
               QXmaxel(k,ll) = max( dg%QXtaylor(k,ll,1),dg%QXtaylor(dg%EL_NBORS(ell,k),ll,1), QXmaxel(k,ll) )
               QYmaxel(k,ll) = max( dg%QYtaylor(k,ll,1),dg%QYtaylor(dg%EL_NBORS(ell,k),ll,1), QYmaxel(k,ll) )

               ZEminel(k,ll) = min( dg%ZEtaylor(k,ll,1),dg%ZEtaylor(dg%EL_NBORS(ell,k),ll,1), ZEminel(k,ll) )
               QXminel(k,ll) = min( dg%QXtaylor(k,ll,1),dg%QXtaylor(dg%EL_NBORS(ell,k),ll,1), QXminel(k,ll) )
               QYminel(k,ll) = min( dg%QYtaylor(k,ll,1),dg%QYtaylor(dg%EL_NBORS(ell,k),ll,1), QYminel(k,ll) )

#ifdef TRACE
               iotamaxel(k,ll) = max( dg%iotataylor(k,ll,1),dg%iotataylor(dg%EL_NBORS(ell,k),ll,1), 
     &              iotamaxel(k,ll) )
               iotaminel(k,ll) = min( dg%iotataylor(k,ll,1),dg%iotataylor(dg%EL_NBORS(ell,k),ll,1), 
     &              iotaminel(k,ll) )
#endif

#ifdef CHEM
               iotamaxel(k,ll) = max( dg%iotataylor(k,ll,1),dg%iotataylor(dg%EL_NBORS(ell,k),ll,1), 
     &              iotamaxel(k,ll) )
               iota2maxel(k,ll) = max( dg%iota2taylor(k,ll,1),dg%iota2taylor(dg%EL_NBORS(ell,k),ll,1),
     &              iota2maxel(k,ll) )
               iotaminel(k,ll) = min( dg%iotataylor(k,ll,1),dg%iotataylor(dg%EL_NBORS(ell,k),ll,1),
     &              iotaminel(k,ll) )
               iota2minel(k,ll) = min( dg%iota2taylor(k,ll,1),dg%iota2taylor(dg%EL_NBORS(ell,k),ll,1),
     &              iota2minel(k,ll) )
#endif
               
            enddo
            
         enddo

      enddo



!.....Must generate linear recostructions at vertices

      dg%ZEtaylorvert = 0.D0
      dg%QXtaylorvert = 0.D0
      dg%Qytaylorvert = 0.D0

#ifdef TRACE
      dg%iotataylorvert = 0.D0
#endif

#ifdef CHEM
      dg%iotataylorvert = 0.D0
      dg%iota2taylorvert = 0.D0
#endif

      do k=1,ne

         if (dg%dofs(k).gt.1) then

            do ll=1,dg%dofs(k)

               do lll=1,3

                  if (ll.eq.1) then

                     dg%ZEtaylorvert(k,ll,lll) = dg%ZEtaylor(k,ll,1) + dg%ZEtaylor(k,ll+1,1)*( dg%xi2vert(k,lll) -dg%xi2BCb(k) )
     &                    + dg%ZEtaylor(k,ll+2,1)*( dg%xi1vert(k,lll) - dg%xi1BCb(k) )
                     dg%QXtaylorvert(k,ll,lll) = dg%QXtaylor(k,ll,1) + dg%QXtaylor(k,ll+1,1)*( dg%xi2vert(k,lll) -dg%xi2BCb(k) )
     &                    + dg%QXtaylor(k,ll+2,1)*( dg%xi1vert(k,lll) - dg%xi1BCb(k) )
                     dg%Qytaylorvert(k,ll,lll) = dg%QYtaylor(k,ll,1) + dg%QYtaylor(k,ll+1,1)*( dg%xi2vert(k,lll) -dg%xi2BCb(k) ) 
     &                    + dg%QYtaylor(k,ll+2,1)*( dg%xi1vert(k,lll) - dg%xi1BCb(k) )

#ifdef TRACE
                     dg%iotataylorvert(k,ll,lll) = dg%iotataylor(k,ll,1) + 
     &                    dg%iotataylor(k,ll+1,1)*( dg%xi2vert(k,lll) - dg%xi2BCb(k) )
     &                    + dg%iotataylor(k,ll+2,1)*( dg%xi1vert(k,lll) - dg%xi1BCb(k) )
#endif

#ifdef CHEM
                     dg%iotataylorvert(k,ll,lll) = dg%iotataylor(k,ll,1) + 
     &                    dg%iotataylor(k,ll+1,1)*( dg%xi2vert(k,lll) - dg%xi2BCb(k) )
     &                    + dg%iotataylor(k,ll+2,1)*( dg%xi1vert(k,lll) - dg%xi1BCb(k) )
                     dg%iota2taylorvert(k,ll,lll) = dg%iota2taylor(k,ll,1) + 
     &                    dg%iota2taylor(k,ll+1,1)*( dg%xi2vert(k,lll) - dg%xi2BCb(k) )
     &                    + dg%iota2taylor(k,ll+2,1)*( dg%xi1vert(k,lll) - dg%xi1BCb(k) )
#endif

                  elseif (ll.gt.1) then

                     dg%ZEtaylorvert(k,ll,lll) = dg%ZEtaylor(k,ll,1) + 
     &                    dg%ZEtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) ),1)*( dg%xi2vert(k,lll) - dg%xi2BCb(k) ) 
     &                    + dg%ZEtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) )+1,1)*( dg%xi1vert(k,lll) - dg%xi1BCb(k) )
                     dg%QXtaylorvert(k,ll,lll) = dg%QXtaylor(k,ll,1) + 
     &                    dg%QXtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) ),1)*( dg%xi2vert(k,lll) - dg%xi2BCb(k) ) 
     &                    + dg%QXtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) )+1,1)*( dg%xi1vert(k,lll) - dg%xi1BCb(k) )
                     dg%Qytaylorvert(k,ll,lll) = dg%QYtaylor(k,ll,1) + 
     &                    dg%QYtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) ),1)*( dg%xi2vert(k,lll) - dg%xi2BCb(k) ) 
     &                    + dg%QYtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) )+1,1)*( dg%xi1vert(k,lll) - dg%xi1BCb(k) )

#ifdef TRACE
                     dg%iotataylorvert(k,ll,lll) = dg%iotataylor(k,ll,1)+
     &                    dg%iotataylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) ),1)*( dg%xi2vert(k,lll) - dg%xi2BCb(k))
     &                    + dg%iotataylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) )+1,1)*( dg%xi1vert(k,lll) - dg%xi1BCb(k) )
#endif

#ifdef CHEM
                     dg%iotataylorvert(k,ll,lll) = dg%iotataylor(k,ll,1)+
     &                    dg%iotataylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) ),1)*( dg%xi2vert(k,lll) - dg%xi2BCb(k))
     &                    + dg%iotataylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) )+1,1)*( dg%xi1vert(k,lll) - dg%xi1BCb(k) )
                     dg%iota2taylorvert(k,ll,lll) = dg%iota2taylor(k,ll,1)+
     &                    dg%iotataylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) ),1)*( dg%xi2vert(k,lll) - dg%xi2BCb(k))
     &                    + dg%iota2taylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) )+1,1)*( dg%xi1vert(k,lll) - dg%xi1BCb(k) )
#endif

                  else

                  endif

               enddo
               
            enddo
            
         elseif (dg%dofs(k).eq.1) then

            cycle

         endif

      enddo

!.....Compute alphas for each variable in each order derivitive

      dg%alphaZE0 = 0.D0
      dg%alphaQX0 = 0.D0
      dg%alphaQY0 = 0.D0

#ifdef TRACE
      dg%alphaiota0 = 0.D0
#endif

#ifdef CHEM
      dg%alphaiota0 = 0.D0
      dg%alphaiota20 = 0.D0
#endif

      do k = 1,ne

         if (dg%dofs(k).gt.1) then

            do ll = 1,dg%dofs(k)

               do lll = 1,3

                  if ( dg%ZEtaylorvert(k,ll,lll).gt.dg%ZEtaylor(k,ll,1).and.
     &                 abs((dg%ZEtaylorvert(k,ll,lll)-dg%ZEtaylor(k,ll,1))).gt.1.0E-15 ) then 

                     dg%alphaZE0(k,ll,lll) = min( 1.D0, ( ZEmaxel(k,ll) - 
     &                    dg%ZEtaylor(k,ll,1) ) / (  dg%ZEtaylorvert(k,ll,lll) - dg%ZEtaylor(k,ll,1) ) )
                     
                  elseif ( dg%ZEtaylorvert(k,ll,lll).eq.dg%ZEtaylor(k,ll,1).or.
     &                    abs((dg%ZEtaylorvert(k,ll,lll)-dg%ZEtaylor(k,ll,1))).le.1.0E-15 ) then

                     dg%alphaZE0(k,ll,lll) = 1.D0

                  elseif (  dg%ZEtaylorvert(k,ll,lll).lt.dg%ZEtaylor(k,ll,1).and.
     &                    abs((dg%ZEtaylorvert(k,ll,lll)-dg%ZEtaylor(k,ll,1))).gt.1.0E-15 ) then

                     dg%alphaZE0(k,ll,lll) = min( 1.D0, ( ZEminel(k,ll)
     &                    - dg%ZEtaylor(k,ll,1) ) / ( dg%ZEtaylorvert(k,ll,lll) - dg%ZEtaylor(k,ll,1) ) )
                     
                  endif

#ifdef TRACE
                  if ( dg%iotataylorvert(k,ll,lll).gt.dg%iotataylor(k,ll,1).and.
     &                 abs((dg%iotataylorvert(k,ll,lll)-dg%iotataylor(k,ll,1))).gt.1.0E-15 ) then 

                     dg%alphaiota0(k,ll,lll) = min(1.D0,( iotamaxel(k,ll)
     &                    -dg%iotataylor(k,ll,1) )/ (dg%iotataylorvert(k,ll,lll) - dg%iotataylor(k,ll,1)))
                     
                  elseif (dg%iotataylorvert(k,ll,lll).eq.dg%iotataylor(k,ll,1).or.
     &                    abs((dg%iotataylorvert(k,ll,lll)-dg%iotataylor(k,ll,1))).le.1.0E-15 ) then

                     dg%alphaiota0(k,ll,lll) = 1.D0

                  elseif (dg%iotataylorvert(k,ll,lll).lt.dg%iotataylor(k,ll,1).
     &                    and.abs((dg%iotataylorvert(k,ll,lll)-dg%iotataylor(k,ll,1))).gt.1.0E-15) then

                     dg%alphaiota0(k,ll,lll) = min( 1.D0,( iotaminel(k,ll)
     &                    -dg%iotataylor(k,ll,1) )/( dg%iotataylorvert(k,ll,lll)-dg%iotataylor(k,ll,1)))
                     
                  endif
#endif
                  
#ifdef CHEM        
                  if ( dg%iotataylorvert(k,ll,lll).gt.dg%iotataylor(k,ll,1).and.
     &                 abs((dg%iotataylorvert(k,ll,lll)-dg%iotataylor(k,ll,1))).gt.1.0E-15 ) then 

                     dg%alphaiota0(k,ll,lll) = min(1.D0,( iotamaxel(k,ll)-
     &                    dg%iotataylor(k,ll,1) )/ (dg%iotataylorvert(k,ll,lll) - dg%iotataylor(k,ll,1)))
                     
                  elseif (dg%iotataylorvert(k,ll,lll).eq.dg%iotataylor(k,ll,1).or.
     &                    abs((dg%iotataylorvert(k,ll,lll)-dg%iotataylor(k,ll,1))).le.1.0E-15 ) then

                     dg%alphaiota0(k,ll,lll) = 1.D0

                  elseif (dg%iotataylorvert(k,ll,lll).lt.dg%iotataylor(k,ll,1).and.
     &                    abs((dg%iotataylorvert(k,ll,lll)-dg%iotataylor(k,ll,1))).gt.1.0E-15) then

                     dg%alphaiota0(k,ll,lll) = min( 1.D0,( iotaminel(k,ll)
     &                    -dg%iotataylor(k,ll,1) )/( dg%iotataylorvert(k,ll,lll)-dg%iotataylor(k,ll,1)))
                     
                  endif

                  if ( dg%iota2taylorvert(k,ll,lll).gt.dg%iota2taylor(k,ll,1).and.
     &                 abs((dg%iota2taylorvert(k,ll,lll)-dg%iota2taylor(k,ll,1))).gt.1.0E-15 ) then 

                     dg%alphaiota20(k,ll,lll) = min(1.D0,( iota2maxel(k,ll)
     &                    -dg%iota2taylor(k,ll,1) )/ (dg%iota2taylorvert(k,ll,lll) - dg%iota2taylor(k,ll,1)))
                     
                  elseif (dg%iota2taylorvert(k,ll,lll).eq.dg%iota2taylor(k,ll,1).or.
     &                    abs((dg%iota2taylorvert(k,ll,lll)-dg%iota2taylor(k,ll,1))).le.1.0E-15 ) then

                     dg%alphaiota20(k,ll,lll) = 1.D0

                  elseif (dg%iota2taylorvert(k,ll,lll).lt.dg%iota2taylor(k,ll,1).and.
     &                    abs((dg%iota2taylorvert(k,ll,lll)-dg%iota2taylor(k,ll,1))).gt.1.0E-15) then

                     dg%alphaiota20(k,ll,lll) = min( 1.D0,( iota2minel(k,ll)
     &                    -dg%iota2taylor(k,ll,1) )/( dg%iota2taylorvert(k,ll,lll)-dg%iota2taylor(k,ll,1)))
                     
                  endif
#endif

                  if ( dg%QXtaylorvert(k,ll,lll).gt.dg%QXtaylor(k,ll,1).and.
     &                 (dg%QXtaylorvert(k,ll,lll)-dg%QXtaylor(k,ll,1)).gt.1.0E-15 ) then !in dg%xi1

                     dg%alphaQX0(k,ll,lll) = min( 1.D0, ( QXmaxel(k,ll) 
     &                    - dg%QXtaylor(k,ll,1) ) / ( dg%QXtaylorvert(k,ll,lll) - dg%QXtaylor(k,ll,1) ) )
                     
                  elseif ( dg%QXtaylorvert(k,ll,lll).eq.dg%QXtaylor(k,ll,1).or.
     &                    (dg%QXtaylorvert(k,ll,lll)-dg%QXtaylor(k,ll,1)).le.1.0E-15  ) then

                     dg%alphaQX0(k,ll,lll) = 1.D0

                  elseif ( dg%QXtaylorvert(k,ll,lll).lt.dg%QXtaylor(k,ll,1).and.
     &                    (dg%QXtaylorvert(k,ll,lll)-dg%QXtaylor(k,ll,1)).gt.1.0E-15 ) then

                     dg%alphaQX0(k,ll,lll) = min( 1.D0, ( QXminel(k,ll)
     &                    - dg%QXtaylor(k,ll,1) ) / ( dg%QXtaylorvert(k,ll,lll) - dg%QXtaylor(k,ll,1) ) )

                  endif

                  if ( dg%QYtaylorvert(k,ll,lll).gt.dg%QYtaylor(k,ll,1).and.
     &                 (dg%QYtaylorvert(k,ll,lll)-dg%QYtaylor(k,ll,1)).gt.1.0E-15 ) then !in dg%xi1

                     dg%alphaQY0(k,ll,lll) = min( 1.D0, ( QYmaxel(k,ll) 
     &                    - dg%QYtaylor(k,ll,1) ) / ( dg%QYtaylorvert(k,ll,lll) - dg%QYtaylor(k,ll,1) ) )
                     
                  elseif ( dg%QYtaylorvert(k,ll,lll).eq.dg%QYtaylor(k,ll,1).or.
     &                    (dg%QYtaylorvert(k,ll,lll)-dg%QYtaylor(k,ll,1)).le.1.0E-15  ) then

                     dg%alphaQY0(k,ll,lll) = 1.D0

                  elseif ( dg%QYtaylorvert(k,ll,lll).lt.dg%QYtaylor(k,ll,1).and.
     &                    (dg%QYtaylorvert(k,ll,lll)-dg%QYtaylor(k,ll,1)).gt.1.0E-15 ) then

                     dg%alphaQY0(k,ll,lll) = min( 1.D0, ( QYminel(k,ll) 
     &                    - dg%QYtaylor(k,ll,1) ) / ( dg%QYtaylorvert(k,ll,lll)  - dg%QYtaylor(k,ll,1) ) )

                  endif

               enddo            !lll

            enddo               !ll

         elseif (dg%dofs(k).eq.1) then

            cycle

         endif

      enddo

!.....Find the prescribed higher limiters by finding smallest local value

      dg%alphaZE = 0.D0
      dg%alphaQX = 0.D0
      dg%alphaQY = 0.D0

#ifdef TRACE
      dg%alphaiota = 0.D0
#endif

#ifdef CHEM
      dg%alphaiota = 0.D0
      dg%alphaiota2 = 0.D0
#endif

      do k=1,ne

         if (dg%dofs(k).gt.1) then

            do ll=1,dg%dofs(k)		

               dg%alphaZE(k,ll) = minval( dg%alphaZE0(k,ll,:) )
               dg%alphaQX(k,ll) = minval( dg%alphaQX0(k,ll,:) )
               dg%alphaQY(k,ll) = minval( dg%alphaQY0(k,ll,:) )

#ifdef TRACE
               dg%alphaiota(k,ll) = minval( dg%alphaiota0(k,ll,:) )
#endif

#ifdef CHEM
               dg%alphaiota(k,ll) = minval( dg%alphaiota0(k,ll,:) )
               dg%alphaiota2(k,ll) = minval( dg%alphaiota20(k,ll,:) )
#endif

            enddo

         elseif (dg%dofs(k).eq.1) then

            cycle

         endif

      enddo

!.... Choose smallest (minimum) alpha for derivative in x or y

      dg%alphaZEm = 0.D0
      dg%alphaQXm = 0.D0
      dg%alphaQYm = 0.D0

#ifdef TRACE
      dg%alphaiotam = 0.D0
#endif

#ifdef CHEM
      dg%alphaiotam = 0.D0
      dg%alphaiota2m = 0.D0
#endif

      do k=1,ne

         if (dg%dofs(k).gt.1) then
            
            do bb = 1,pdg_el(k)

               if( (bb+1)*(bb+2)/2.le.dg%dofs(k) ) then

                  dg%alphaZEm(k,(bb*(bb+1))/2 + 1) = 
     &                 minval( dg%alphaZE(k,(bb*(bb-1))/2 +1:(bb*(bb+1)) / 2  ) )
                  dg%alphaQXm(k,(bb*(bb+1))/2 + 1) = 
     &                 minval( dg%alphaQX(k,(bb*(bb-1))/2 +1:(bb*(bb+1)) / 2  ) )
                  dg%alphaQYm(k,(bb*(bb+1))/2 + 1) = 
     &                 minval( dg%alphaQY(k,(bb*(bb-1))/2 +1:(bb*(bb+1)) / 2  ) )

#ifdef TRACE
                  dg%alphaiotam(k,(bb*(bb+1))/2 + 1) = 
     &                 minval( dg%alphaiota(k,(bb*(bb-1))/2 +1:(bb*(bb+1)) / 2 ) )
#endif

#ifdef CHEM
                  dg%alphaiotam(k,(bb*(bb+1))/2 + 1) = 
     &                 minval( dg%alphaiota(k,(bb*(bb-1))/2 +1:(bb*(bb+1)) / 2  ) )
                  dg%alphaiota2m(k,(bb*(bb+1))/2 + 1) = 
     &                 minval( dg%alphaiota2(k,(bb*(bb-1))/2 +1:(bb*(bb+1)) / 2 ) )
#endif

               endif
               
            enddo

         elseif (dg%dofs(k).eq.1) then

            cycle

         endif

      enddo

!.....Use max higher derivative values for final limiter value

      dg%alphaZE_max = 0.D0
      dg%alphaQX_max = 0.D0
      dg%alphaQY_max = 0.D0

#ifdef TRACE
      dg%alphaiota_max = 0.D0
#endif

#ifdef CHEM
      dg%alphaiota_max = 0.D0
      dg%alphaiota2_max = 0.D0
#endif

      do k=1,ne

         if (dg%dofs(k).gt.1) then

            do bb =1,pdg_el(k)

               if( (bb+1)*(bb+2)/2.le.dg%dofs(k)) then

                  dg%alphaZE_max(k,(bb*(bb+1))/2 + 1) = fd*maxval( 
     &                 dg%alphaZEm(k,(bb*(bb+1))/2 + 1:dg%dofs(k)) )
                  dg%alphaQX_max(k,(bb*(bb+1))/2 + 1) = fd*maxval( 
     &                 dg%alphaQXm(k,(bb*(bb+1))/2 + 1:dg%dofs(k)) )
                  dg%alphaQY_max(k,(bb*(bb+1))/2 + 1) = fd*maxval( 
     &                 dg%alphaQYm(k,(bb*(bb+1))/2 + 1:dg%dofs(k)) )

#ifdef TRACE
                  dg%alphaiota_max(k,(bb*(bb+1))/2 + 1) = fd*maxval( 
     &                 dg%alphaiotam(k,(bb*(bb+1))/2 + 1:dg%dofs(k)) )
#endif

#ifdef CHEM
                  dg%alphaiota_max(k,(bb*(bb+1))/2 + 1) = fd*maxval( 
     &                 dg%alphaiotam(k,(bb*(bb+1))/2 + 1:dg%dofs(k)) )
                  dg%alphaiota2_max(k,(bb*(bb+1))/2 + 1) = fd*maxval( 
     &                 dg%alphaiota2m(k,(bb*(bb+1))/2 + 1:dg%dofs(k)) )
#endif

               endif

            enddo

         elseif (dg%dofs(k).eq.1) then

            cycle

         endif

      enddo

!.....Limit on the Master element in the Taylor basis, via reconstruction 
!.....of unconstrained solutions with alpha constraints

      dg%limitZE = 0.D0
      dg%limitQX = 0.D0
      dg%limitQY = 0.D0

      dg%lim_count_roll = 0

#ifdef TRACE
      dg%limitiota = 0.D0
#endif

#ifdef CHEM
      dg%limitiota = 0.D0
      dg%limitiota2 = 0.D0
#endif

      do k=1,ne

         dg%lim_count = 0

         if (dg%dofs(k).gt.1) then

            do ll=1,dg%dofs(k)

               if ( ll.eq.1 ) then

                  dg%limitZE(k,ll) = dg%ZEtaylor(k,ll,1)
                  dg%limitQX(k,ll) = dg%QXtaylor(k,ll,1)
                  dg%limitQY(k,ll) = dg%QYtaylor(k,ll,1) 

#ifdef TRACE
                  dg%limitiota(k,ll) = dg%iotataylor(k,ll,1)
#endif

#ifdef CHEM
                  dg%limitiota(k,ll) = dg%iotataylor(k,ll,1)
                  dg%limitiota2(k,ll) = dg%iota2taylor(k,ll,1)
#endif

               elseif ( ll.ge.2 ) then
                  
                  do bb=1,pdg_el(k)

                     if ( ll.le.( ( (bb+1)*(bb+2)) / 2.D0 ).and.(ll.gt.
     &                    (bb*(bb+1)/2.D0) ) ) then

                        dg%limitZE(k,ll) = dg%alphaZE_max(k,(bb*(bb+1))/2 + 1) 
     &                       * dg%ZEtaylor(k,ll,1)
                        dg%limitQX(k,ll) = dg%alphaQX_max(k,(bb*(bb+1))/2 + 1) 
     &                       * dg%ZEtaylor(k,ll,1)
                        dg%limitQY(k,ll) = dg%alphaQY_max(k,(bb*(bb+1))/2 + 1) 
     &                       * dg%ZEtaylor(k,ll,1)

#ifdef TRACE
                        dg%limitiota(k,ll) = dg%alphaiota_max(k,(bb*(bb+1))/2 + 1) 
     &                       * dg%iotataylor(k,ll,1)
#endif

#ifdef CHEM
                        dg%limitiota(k,ll) = dg%alphaiota_max(k,(bb*(bb+1))/2 + 1) 
     &                       * dg%iotataylor(k,ll,1)
                        dg%limitiota2(k,ll) = dg%alphaiota2_max(k,(bb*(bb+1))/2 + 1) 
     &                       * dg%iota2taylor(k,ll,1)
#endif


!$$$  ! Make a counter to track limiting
!$$$  
!$$$  if ( ( dg%alphaZE_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
!$$$  &                       dg%alphaQX_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
!$$$  &                       dg%alphaQY_max(k,(bb*(bb+1))/2 + 1).lt.1.D0 ).and.
!$$$  &                       chem_flag.ne.1.and.tracer_flag.ne.1 ) then
!$$$  
!$$$  dg%lim_count = 1  
!$$$  
!$$$  elseif ( (dg%alphaZE_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
!$$$  &                          dg%alphaQX_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
!$$$  &                          dg%alphaQY_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
!$$$  &                          dg%alphaiota_max(k,(bb*(bb+1))/2 + 1).lt.1.D0 ).and.
!$$$  &                          chem_flag.ne.1.and.tracer_flag.eq.1 ) then
!$$$  
!$$$  dg%lim_count = 1   
!$$$  
!$$$  elseif ( (dg%alphaZE_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
!$$$  &                          dg%alphaQX_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
!$$$  &                          dg%alphaQY_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
!$$$  &                          dg%alphaiota_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
!$$$  &                          dg%alphaiota2_max(k,(bb*(bb+1))/2 + 1).lt.1.D0 ).and.
!$$$  &                          chem_flag.eq.1 ) then
!$$$  
!$$$  dg%lim_count = 1 
!$$$  
!$$$  endif
                        
                     endif
                     
                  enddo
                  
               endif

            enddo

         elseif (dg%dofs(k).eq.1) then

            cycle

         endif

                                !dg%lim_count_roll = dg%lim_count_roll + dg%lim_count

      enddo

!.....Transform back to the Dubiner basis (multiply by dg%NmatrixInv),

      dg%ZEconst =  0.D0
      dg%QXconst =  0.D0
      dg%QYconst =  0.D0

#ifdef TRACE
      dg%iotaconst =  0.D0
#endif

#ifdef CHEM
      dg%iotaconst =  0.D0
      dg%iota2const =  0.D0
#endif

      do k=1,ne

         if (dg%dofs(k).gt.1) then
                                !do lll=1,3

            do ll=1,dg%dofs(k)

               do ss=1,dg%dofs(k)

                  dg%ZEconst(k,ll) = dg%ZEconst(k,ll) + dg%NmatrixInv(k,ll,ss,dg%dofs(k)) 
     &                 * dg%limitZE(k,ss)
                  dg%QXconst(k,ll) = dg%QXconst(k,ll) + dg%NmatrixInv(k,ll,ss,dg%dofs(k)) 
     &                 * dg%limitQX(k,ss)
                  dg%QYconst(k,ll) = dg%QYconst(k,ll) + dg%NmatrixInv(k,ll,ss,dg%dofs(k)) 
     &                 * dg%limitQY(k,ss)

#ifdef TRACE
                  dg%iotaconst(k,ll) = dg%iotaconst(k,ll) + 
     &                 dg%NmatrixInv(k,ll,ss,dg%dofs(k)) * dg%limitiota(k,ss)
#endif

#ifdef CHEM
                  dg%iotaconst(k,ll) = dg%iotaconst(k,ll) + 
     &                 dg%NmatrixInv(k,ll,ss,dg%dofs(k)) * dg%limitiota(k,ss)
                  dg%iota2const(k,ll) = dg%iota2const(k,ll) + 
     &                 dg%NmatrixInv(k,ll,ss,dg%dofs(k)) * dg%limitiota2(k,ss)
#endif


               enddo

            enddo

         elseif (dg%dofs(k).eq.1) then

            cycle

         endif

      enddo

!.....Set limit values

      do k=1,ne

         if (dg%dofs(k).gt.1) then

            do ll = 1,dg%dofs(k)

               dg%ZE(ll,k,dg%irk+1) = dg%ZEconst(k,ll)
               dg%QX(ll,k,dg%irk+1) = dg%QXconst(k,ll)
               dg%QY(ll,k,dg%irk+1) = dg%QYconst(k,ll)

#ifdef TRACE
               dg%iota(ll,k,dg%irk+1) = dg%iotaconst(k,ll)
#endif

#ifdef CHEM
               dg%iota(ll,k,dg%irk+1) = dg%iotaconst(k,ll)
               dg%iota2(ll,k,dg%irk+1) = dg%iota2const(k,ll)
#endif

            enddo

         elseif (dg%dofs(k).eq.1) then

            cycle

         endif

      enddo
      
      return
      end subroutine

!***********************************************************************
!     
!     SUBROUTINE SLOPELIMITER10()
!     
!     Written 2010
!     
!     This subroutine selects the first *adapted* vertex limiter based on
!     a Taylor Polynomial basis, and is consistent with p_adaptation
!     and works for any p
!     
!     -- cem
!     
!***********************************************************************

      SUBROUTINE SLOPELIMITER10()

!.....Use appropriate modules

      USE GLOBAL
      USE DG

#ifdef CMPI
      USE MESSENGER
      USE MESSENGER_ELEM
#endif

      IMPLICIT NONE

      Integer k,ll,ss,lll,ell,bb,maxneigh

      REAL(SZ) fd,marea,fde

!.....We work over the master element
!.....Set initial values

      fd = dg%slope_weight         ! reduces diffusion fd = 1 => full diffusion
      fde = fd                  ! add weight for lower order pieces (fd<1 => stronger limiting)     

      DO k=1,mne

         if (dg%dofs(k).gt.1) then

            DO ll = 1,dg%dofs(k)

               dg%ZE(ll,k,dg%NRK+2) = dg%ZE(ll,k,dg%IRK+1)
               dg%QX(ll,k,dg%NRK+2) = dg%QX(ll,k,dg%IRK+1)
               dg%QY(ll,k,dg%NRK+2) = dg%QY(ll,k,dg%IRK+1)

#ifdef TRACE
               dg%iota(ll,k,dg%NRK+2) = dg%iota(ll,k,dg%IRK+1)
#endif

#ifdef CHEM
               dg%iota(ll,k,dg%NRK+2) = dg%iota(ll,k,dg%IRK+1)
               dg%iota2(ll,k,dg%NRK+2) = dg%iota2(ll,k,dg%IRK+1)
#endif


            ENDDO

         elseif (dg%dofs(k).eq.1) then

            cycle

         endif

      ENDDO   

      marea = 2.D0              !master elements area

!.....Convert initial values to the Taylor basis (multiply by dg%Nmatrix) on base element


      dg%ZEtaylor = 0.D0 
      dg%QXtaylor = 0.D0
      dg%QYtaylor = 0.D0

#ifdef TRACE
      dg%iotataylor = 0.D0
#endif

#ifdef CHEM
      dg%iotataylor = 0.D0
      dg%iota2taylor = 0.D0
#endif

      do k=1,mne

         if (dg%dofs(k).gt.1) then

            do ll=1,dg%dofs(k)
               
               do ss=1,dg%dofs(k)

                  dg%ZEtaylor(k,ll,1) = dg%ZEtaylor(k,ll,1) + dg%Nmatrix(k,ll,ss,dg%dofs(k)) * dg%ZE(ss,k,dg%nrk+2)
                  dg%QXtaylor(k,ll,1) = dg%QXtaylor(k,ll,1) + dg%Nmatrix(k,ll,ss,dg%dofs(k))* dg%QX(ss,k,dg%nrk+2)
                  dg%QYtaylor(k,ll,1) = dg%QYtaylor(k,ll,1) + dg%Nmatrix(k,ll,ss,dg%dofs(k)) * dg%QY(ss,k,dg%nrk+2)

#ifdef TRACE
                  dg%iotataylor(k,ll,1) = dg%iotataylor(k,ll,1) + dg%Nmatrix(k,ll,ss,dg%dofs(k)) * dg%iota(ss,k,dg%nrk+2)
#endif

#ifdef CHEM
                  dg%iotataylor(k,ll,1) = dg%iotataylor(k,ll,1) + dg%Nmatrix(k,ll,ss,dg%dofs(k)) * dg%iota(ss,k,dg%nrk+2)
                  dg%iota2taylor(k,ll,1) = dg%iota2taylor(k,ll,1) + dg%Nmatrix(k,ll,ss,dg%dofs(k)) * dg%iota2(ss,k,dg%nrk+2)
#endif
                  
               enddo

            enddo

         elseif (dg%dofs(k).eq.1) then

            cycle

         endif

      enddo

!.....Find values at vertices of base elements and neighbors


      dg%ZEmax = -100.D0
      dg%QXmax = -100.D0
      dg%QYmax = -100.D0
      dg%ZEmin = 100.D0
      dg%QXmin = 100.D0
      dg%QYmin = 100.D0

#ifdef TRACE
      dg%iotamax = -100.D0
      dg%iotamin = 100.D0
#endif

#ifdef CHEM
      dg%iotamax = -100.D0
      dg%iota2max = -100.D0
      dg%iotamin = 100.D0
      dg%iota2min = 100.D0
#endif


#ifdef CMPI

      CALL UPDATER_ELEM_MOD2(dg%ZEtaylor,dg%QXtaylor,dg%QYtaylor,1,3)

#ifdef TRACE
      CALL UPDATER_ELEM_MOD2(dg%iotataylor,dg%QXtaylor,dg%QYtaylor,1,2 )  
#endif

#ifdef CHEM
      CALL UPDATER_ELEM_MOD2(dg%iotataylor,dg%iota2taylor,dg%QYtaylor,1,2 )  
#endif

      CALL UPDATERV(dg%ZEmin)
      CALL UPDATERV(dg%ZEmax)
      CALL UPDATERV(dg%QXmin)
      CALL UPDATERV(dg%QXmax)
      CALL UPDATERV(dg%QYmin)
      CALL UPDATERV(dg%QYmax)

#ifdef TRACE
         CALL UPDATERV(dg%iotamax)
         CALL UPDATERV(dg%iotamin)
#endif

#ifdef CHEM
         CALL UPDATERV(dg%iotamax)
         CALL UPDATERV(dg%iotamin)
         CALL UPDATERV(dg%iota2max)
         CALL UPDATERV(dg%iota2min)
#endif

#endif

      do ell=1,mnp

         do ll=1,minval(dg%dofs(neigh_elem(ell,1:nneigh_elem(ell))))

!.....Find max and min values over polynomial coefficients

            dg%ZEmax(ell,ll) = max(maxval( dg%ZEtaylor(NEIGH_ELEM(ell,1:nneigh_elem(ell)),ll,1) )
     &           , dg%ZEmax(ell,ll))
            dg%QXmax(ell,ll) = max(maxval( dg%QXtaylor(NEIGH_ELEM(ell,1:nneigh_elem(ell)),ll,1) )
     &           , dg%QXmax(ell,ll))
            dg%QYmax(ell,ll) = max(maxval( dg%QYtaylor(NEIGH_ELEM(ell,1:nneigh_elem(ell)),ll,1) )
     &           , dg%QYmax(ell,ll))
            dg%ZEmin(ell,ll) = min(minval( dg%ZEtaylor(NEIGH_ELEM(ell,1:nneigh_elem(ell)),ll,1) )
     &           , dg%ZEmin(ell,ll))
            dg%QXmin(ell,ll) = min(minval( dg%QXtaylor(NEIGH_ELEM(ell,1:nneigh_elem(ell)),ll,1) )
     &           , dg%QXmin(ell,ll))
            dg%QYmin(ell,ll) = min(minval( dg%QYtaylor(NEIGH_ELEM(ell,1:nneigh_elem(ell)),ll,1) )
     &           , dg%QYmin(ell,ll))

#ifdef TRACE
            dg%iotamax(ell,ll) = max(maxval( dg%iotataylor(NEIGH_ELEM(ell,1:nneigh_elem(ell)),ll,1) )
     &           , dg%iotamax(ell,ll))
            dg%iotamin(ell,ll) = min(minval( dg%iotataylor(NEIGH_ELEM(ell,1:nneigh_elem(ell)),ll,1) )
     &           , dg%iotamin(ell,ll))
#endif

#ifdef CHEM
            dg%iotamax(ell,ll) = max(maxval( dg%iotataylor(NEIGH_ELEM(ell,1:nneigh_elem(ell)),ll,1) )
     &           , dg%iotamax(ell,ll))
            dg%iota2max(ell,ll) = max(maxval( dg%iota2taylor(NEIGH_ELEM(ell,1:nneigh_elem(ell)),ll,1) )
     &           , dg%iota2max(ell,ll))
            dg%iotamin(ell,ll) = min(minval( dg%iotataylor(NEIGH_ELEM(ell,1:nneigh_elem(ell)),ll,1) )
     &           , dg%iotamin(ell,ll))
            dg%iota2min(ell,ll) = min(minval( dg%iota2taylor(NEIGH_ELEM(ell,1:nneigh_elem(ell)),ll,1) )
     &           , dg%iota2min(ell,ll))
#endif
            
         enddo

      enddo


#ifdef CMPI

      CALL UPDATERV(dg%ZEmin)
      CALL UPDATERV(dg%ZEmax)
      CALL UPDATERV(dg%QXmin)
      CALL UPDATERV(dg%QXmax)
      CALL UPDATERV(dg%QYmin)
      CALL UPDATERV(dg%QYmax)

#ifdef TRACE
      CALL UPDATERV(dg%iotamax)
      CALL UPDATERV(dg%iotamin)
#endif

#ifdef CHEM
      CALL UPDATERV(dg%iotamax)
      CALL UPDATERV(dg%iotamin)
      CALL UPDATERV(dg%iota2max)
      CALL UPDATERV(dg%iota2min)
#endif

#endif


!.....Must generate linear recostructions at vertices

      dg%ZEtaylorvert = 0.D0
      dg%QXtaylorvert = 0.D0
      dg%Qytaylorvert = 0.D0

#ifdef TRACE
      dg%iotataylorvert = 0.D0
#endif

#ifdef CHEM
      dg%iotataylorvert = 0.D0
      dg%iota2taylorvert = 0.D0
#endif

      do k=1,mne

         if (dg%dofs(k).gt.1) then

            do ll=1,dg%dofs(k)

               do lll=1,3

                  if (ll.eq.1) then

                     dg%ZEtaylorvert(k,ll,lll) = dg%ZEtaylor(k,ll,1) + dg%ZEtaylor(k,ll+1,1)*( dg%xi2vert(k,lll) -dg%xi2BCb(k) )
     &                    + dg%ZEtaylor(k,ll+2,1)*( dg%xi1vert(k,lll) - dg%xi1BCb(k) )
                     dg%QXtaylorvert(k,ll,lll) = dg%QXtaylor(k,ll,1) + dg%QXtaylor(k,ll+1,1)*( dg%xi2vert(k,lll) -dg%xi2BCb(k) )
     &                    + dg%QXtaylor(k,ll+2,1)*( dg%xi1vert(k,lll) - dg%xi1BCb(k) )
                     dg%Qytaylorvert(k,ll,lll) = dg%QYtaylor(k,ll,1) + dg%QYtaylor(k,ll+1,1)*( dg%xi2vert(k,lll) -dg%xi2BCb(k) ) 
     &                    + dg%QYtaylor(k,ll+2,1)*( dg%xi1vert(k,lll) - dg%xi1BCb(k) )

#ifdef TRACE
                     dg%iotataylorvert(k,ll,lll) = dg%iotataylor(k,ll,1) + 
     &                    dg%iotataylor(k,ll+1,1)*( dg%xi2vert(k,lll) - dg%xi2BCb(k) )
     &                    + dg%iotataylor(k,ll+2,1)*( dg%xi1vert(k,lll) - dg%xi1BCb(k) )
#endif

#ifdef CHEM
                     dg%iotataylorvert(k,ll,lll) = dg%iotataylor(k,ll,1) + 
     &                    dg%iotataylor(k,ll+1,1)*( dg%xi2vert(k,lll) - dg%xi2BCb(k) )
     &                    + dg%iotataylor(k,ll+2,1)*( dg%xi1vert(k,lll) - dg%xi1BCb(k) )
                     dg%iota2taylorvert(k,ll,lll) = dg%iota2taylor(k,ll,1) + 
     &                    dg%iota2taylor(k,ll+1,1)*( dg%xi2vert(k,lll) - dg%xi2BCb(k) )
     &                    + dg%iota2taylor(k,ll+2,1)*( dg%xi1vert(k,lll) - dg%xi1BCb(k) )
#endif

                  elseif (ll.gt.1) then

                     dg%ZEtaylorvert(k,ll,lll) = dg%ZEtaylor(k,ll,1) + 
     &                    dg%ZEtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) ),1)*( dg%xi2vert(k,lll) - dg%xi2BCb(k) ) 
     &                    + dg%ZEtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) )+1,1)*( dg%xi1vert(k,lll) - dg%xi1BCb(k) )
                     dg%QXtaylorvert(k,ll,lll) = dg%QXtaylor(k,ll,1) + 
     &                    dg%QXtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) ),1)*( dg%xi2vert(k,lll) - dg%xi2BCb(k) ) 
     &                    + dg%QXtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) )+1,1)*( dg%xi1vert(k,lll) - dg%xi1BCb(k) )
                     dg%Qytaylorvert(k,ll,lll) = dg%QYtaylor(k,ll,1) + 
     &                    dg%QYtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) ),1)*( dg%xi2vert(k,lll) - dg%xi2BCb(k) ) 
     &                    + dg%QYtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) )+1,1)*( dg%xi1vert(k,lll) - dg%xi1BCb(k) )

#ifdef TRACE
                     dg%iotataylorvert(k,ll,lll) = dg%iotataylor(k,ll,1)+
     &                    dg%iotataylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) ),1)*( dg%xi2vert(k,lll) - dg%xi2BCb(k))
     &                    + dg%iotataylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) )+1,1)*( dg%xi1vert(k,lll) - dg%xi1BCb(k) )
#endif

#ifdef CHEM
                     dg%iotataylorvert(k,ll,lll) = dg%iotataylor(k,ll,1)+
     &                    dg%iotataylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) ),1)*( dg%xi2vert(k,lll) - dg%xi2BCb(k))
     &                    + dg%iotataylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) )+1,1)*( dg%xi1vert(k,lll) - dg%xi1BCb(k) )
                     dg%iota2taylorvert(k,ll,lll) = dg%iota2taylor(k,ll,1)+
     &                    dg%iotataylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) ),1)*( dg%xi2vert(k,lll) - dg%xi2BCb(k))
     &                    + dg%iota2taylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) )+1,1)*( dg%xi1vert(k,lll) - dg%xi1BCb(k) )
#endif

                  else

                  endif

               enddo
               
            enddo
            
         elseif (dg%dofs(k).eq.1) then

            cycle

         endif

      enddo

!.....Compute alphas for each variable in each order derivitive


      dg%alphaZE0 = 0.D0
      dg%alphaQX0 = 0.D0
      dg%alphaQY0 = 0.D0

#ifdef TRACE
      dg%alphaiota0 = 0.D0
#endif

#ifdef CHEM
      dg%alphaiota0 = 0.D0
      dg%alphaiota20 = 0.D0
#endif

      do k=1,mne

         if (dg%dofs(k).gt.1) then

            do lll=1,3
               
               do ll=1,dg%dofs(k)

                  if (dg%ZEmin(nm(k,lll),ll).ne.dg%ZEmax(nm(k,lll),ll)) then

                     if ( ( dg%ZEtaylorvert(k,ll,lll).gt.dg%ZEtaylor(k,ll,1) ).and.
     &                    ( abs(dg%ZEtaylorvert(k,ll,lll)-dg%ZEtaylor(k,ll,1)).gt.1.0E-15 ).and.
     &                    ( dg%ZEmax(nm(k,lll),ll).ne.dg%ZEtaylor(k,ll,1) ) ) then  

                        dg%alphaZE0(k,ll,lll) = min(1.D0,  ( dg%ZEmax(nm(k,lll),ll)
     &                       - dg%ZEtaylor(k,ll,1) )/ (dg%ZEtaylorvert(k,ll,lll) - dg%ZEtaylor(k,ll,1)))

                                !adapted part

                     elseif ( ( dg%ZEtaylorvert(k,ll,lll).gt.dg%ZEtaylor(k,ll,1) ).and.
     &                       ( abs(dg%ZEtaylorvert(k,ll,lll)-dg%ZEtaylor(k,ll,1)).gt.1.0E-15 ).and.
     &                       ( dg%ZEmax(nm(k,lll),ll).eq.dg%ZEtaylor(k,ll,1) ) ) then 

                        dg%alphaZE0(k,ll,lll) = min(fd, abs( ( dg%ZEmax(nm(k,lll),ll)
     &                       - dg%ZEmin(nm(k,lll),ll) )/(dg%ZEtaylorvert(k,ll,lll) - dg%ZEmax(nm(k,lll),ll)) ) ) 
                        

                     elseif ( (dg%ZEtaylorvert(k,ll,lll).lt.dg%ZEtaylor(k,ll,1) )
     &                       .and.( abs(dg%ZEtaylorvert(k,ll,lll)-dg%ZEtaylor(k,ll,1)).gt.1.0E-15 ).and.
     &                       ( dg%ZEmin(nm(k,lll),ll).ne.dg%ZEtaylor(k,ll,1) ) ) then 

                        dg%alphaZE0(k,ll,lll) = min( 1.D0,( dg%ZEmin(nm(k,lll),ll)
     &                       - dg%ZEtaylor(k,ll,1) )/( dg%ZEtaylorvert(k,ll,lll)-dg%ZEtaylor(k,ll,1)))

                                !adapted part

                     elseif ( (dg%ZEtaylorvert(k,ll,lll).lt.dg%ZEtaylor(k,ll,1) )
     &                       .and.( abs(dg%ZEtaylorvert(k,ll,lll)-dg%ZEtaylor(k,ll,1)).gt.1.0E-15 ).and.
     &                       ( dg%ZEmin(nm(k,lll),ll).eq.dg%ZEtaylor(k,ll,1) ) ) then 

                        dg%alphaZE0(k,ll,lll) = min( fd,abs( ( dg%ZEmin(nm(k,lll),ll)
     &                       - dg%ZEmax(nm(k,lll),ll) )/( dg%ZEtaylorvert(k,ll,lll)- dg%ZEmin(nm(k,lll),ll))) )
                        
                        
                     elseif ( ( dg%ZEtaylorvert(k,ll,lll).eq.dg%ZEtaylor(k,ll,1) ).or.
     &                       ( abs(dg%ZEtaylorvert(k,ll,lll)-dg%ZEtaylor(k,ll,1)).le.1.0E-15 ) ) then

                        dg%alphaZE0(k,ll,lll) = 1.D0

                     endif

                  else

                     dg%alphaZE0(k,ll,lll) = 1.D0

                  endif
#ifdef TRACE
                  if (dg%iotamin(nm(k,lll),ll).ne.dg%iotamax(nm(k,lll),ll)) then

                     if ( ( dg%iotataylorvert(k,ll,lll).gt.dg%iotataylor(k,ll,1) ).and.
     &                    ( abs(dg%iotataylorvert(k,ll,lll)-dg%iotataylor(k,ll,1)).gt.1.0E-15 ).and.
     &                    ( dg%iotamax(nm(k,lll),ll).ne.dg%iotataylor(k,ll,1) ) ) then  

                        dg%alphaiota0(k,ll,lll) = min(1.D0,  ( dg%iotamax(nm(k,lll),ll)
     &                       - dg%iotataylor(k,ll,1) )/ (dg%iotataylorvert(k,ll,lll) - dg%iotataylor(k,ll,1)))

                                !adapted part

                     elseif ( ( dg%iotataylorvert(k,ll,lll).gt.dg%iotataylor(k,ll,1) ).and.
     &                       ( abs(dg%iotataylorvert(k,ll,lll)-dg%iotataylor(k,ll,1)).gt.1.0E-15 ).and.
     &                       ( dg%iotamax(nm(k,lll),ll).eq.dg%iotataylor(k,ll,1) ) ) then 

                        dg%alphaiota0(k,ll,lll) = min(fd, abs( ( dg%iotamax(nm(k,lll),ll)
     &                       - dg%iotamin(nm(k,lll),ll) )/(dg%iotataylorvert(k,ll,lll) - dg%iotamax(nm(k,lll),ll)) ) ) 
                        

                     elseif ( (dg%iotataylorvert(k,ll,lll).lt.dg%iotataylor(k,ll,1) )
     &                       .and.( abs(dg%iotataylorvert(k,ll,lll)-dg%iotataylor(k,ll,1)).gt.1.0E-15 ).and.
     &                       ( dg%iotamin(nm(k,lll),ll).ne.dg%iotataylor(k,ll,1) ) ) then 

                        dg%alphaiota0(k,ll,lll) = min( 1.D0,( dg%iotamin(nm(k,lll),ll)
     &                       - dg%iotataylor(k,ll,1) )/( dg%iotataylorvert(k,ll,lll)-dg%iotataylor(k,ll,1)))

                                !adapted part

                     elseif ( (dg%iotataylorvert(k,ll,lll).lt.dg%iotataylor(k,ll,1) )
     &                       .and.( abs(dg%iotataylorvert(k,ll,lll)-dg%iotataylor(k,ll,1)).gt.1.0E-15 ).and.
     &                       ( dg%iotamin(nm(k,lll),ll).eq.dg%iotataylor(k,ll,1) ) ) then 

                        dg%alphaiota0(k,ll,lll) = min( fd,abs( ( dg%iotamin(nm(k,lll),ll)
     &                       - dg%iotamax(nm(k,lll),ll) )/( dg%iotataylorvert(k,ll,lll)- dg%iotamin(nm(k,lll),ll))) )
                        
                        
                     elseif ( ( dg%iotataylorvert(k,ll,lll).eq.dg%iotataylor(k,ll,1) ).or.
     &                       ( abs(dg%iotataylorvert(k,ll,lll)-dg%iotataylor(k,ll,1)).le.1.0E-15 ) ) then

                        dg%alphaiota0(k,ll,lll) = 1.D0

                     endif

                  else

                     dg%alphaiota0(k,ll,lll) = 1.D0

                  endif
#endif
                  
#ifdef CHEM
                  if (dg%iotamin(nm(k,lll),ll).ne.dg%iotamax(nm(k,lll),ll)) then

                     if ( ( dg%iotataylorvert(k,ll,lll).gt.dg%iotataylor(k,ll,1) ).and.
     &                    ( abs(dg%iotataylorvert(k,ll,lll)-dg%iotataylor(k,ll,1)).gt.1.0E-15 ).and.
     &                    ( dg%iotamax(nm(k,lll),ll).ne.dg%iotataylor(k,ll,1) ) ) then  

                        dg%alphaiota0(k,ll,lll) = min(1.D0,  ( dg%iotamax(nm(k,lll),ll)
     &                       - dg%iotataylor(k,ll,1) )/ (dg%iotataylorvert(k,ll,lll) - dg%iotataylor(k,ll,1)))

                                !adapted part

                     elseif ( ( dg%iotataylorvert(k,ll,lll).gt.dg%iotataylor(k,ll,1) ).and.
     &                       ( abs(dg%iotataylorvert(k,ll,lll)-dg%iotataylor(k,ll,1)).gt.1.0E-15 ).and.
     &                       ( dg%iotamax(nm(k,lll),ll).eq.dg%iotataylor(k,ll,1) ) ) then 

                        dg%alphaiota0(k,ll,lll) = min(fd, abs( ( dg%iotamax(nm(k,lll),ll)
     &                       - dg%iotamin(nm(k,lll),ll) )/(dg%iotataylorvert(k,ll,lll) - dg%iotamax(nm(k,lll),ll)) ) ) 
                        

                     elseif ( (dg%iotataylorvert(k,ll,lll).lt.dg%iotataylor(k,ll,1) )
     &                       .and.( abs(dg%iotataylorvert(k,ll,lll)-dg%iotataylor(k,ll,1)).gt.1.0E-15 ).and.
     &                       ( dg%iotamin(nm(k,lll),ll).ne.dg%iotataylor(k,ll,1) ) ) then 

                        dg%alphaiota0(k,ll,lll) = min( 1.D0,( dg%iotamin(nm(k,lll),ll)
     &                       - dg%iotataylor(k,ll,1) )/( dg%iotataylorvert(k,ll,lll)-dg%iotataylor(k,ll,1)))

                                !adapted part

                     elseif ( (dg%iotataylorvert(k,ll,lll).lt.dg%iotataylor(k,ll,1) )
     &                       .and.( abs(dg%iotataylorvert(k,ll,lll)-dg%iotataylor(k,ll,1)).gt.1.0E-15 ).and.
     &                       ( dg%iotamin(nm(k,lll),ll).eq.dg%iotataylor(k,ll,1) ) ) then 

                        dg%alphaiota0(k,ll,lll) = min( fd,abs( ( dg%iotamin(nm(k,lll),ll)
     &                       - dg%iotamax(nm(k,lll),ll) )/( dg%iotataylorvert(k,ll,lll)- dg%iotamin(nm(k,lll),ll))) )
                        
                        
                     elseif ( ( dg%iotataylorvert(k,ll,lll).eq.dg%iotataylor(k,ll,1) ).or.
     &                       ( abs(dg%iotataylorvert(k,ll,lll)-dg%iotataylor(k,ll,1)).le.1.0E-15 ) ) then

                        dg%alphaiota0(k,ll,lll) = 1.D0

                     endif

                  else

                     dg%alphaiota0(k,ll,lll) = 1.D0

                  endif

                  if (dg%iota2min(nm(k,lll),ll).ne.dg%iota2max(nm(k,lll),ll)) then

                     if ( ( dg%iota2taylorvert(k,ll,lll).gt.dg%iota2taylor(k,ll,1) ).and.
     &                    ( abs(dg%iota2taylorvert(k,ll,lll)-dg%iota2taylor(k,ll,1)).gt.1.0E-15 ).and.
     &                    ( dg%iota2max(nm(k,lll),ll).ne.dg%iota2taylor(k,ll,1) ) ) then  

                        dg%alphaiota20(k,ll,lll) = min(1.D0,  ( dg%iota2max(nm(k,lll),ll)
     &                       - dg%iota2taylor(k,ll,1) )/ (dg%iota2taylorvert(k,ll,lll) - dg%iota2taylor(k,ll,1)))

                                !adapted part

                     elseif ( ( dg%iota2taylorvert(k,ll,lll).gt.dg%iota2taylor(k,ll,1) ).and.
     &                       ( abs(dg%iota2taylorvert(k,ll,lll)-dg%iota2taylor(k,ll,1)).gt.1.0E-15 ).and.
     &                       ( dg%iota2max(nm(k,lll),ll).eq.dg%iota2taylor(k,ll,1) ) ) then 

                        dg%alphaiota20(k,ll,lll) = min(fd, abs( ( dg%iota2max(nm(k,lll),ll)
     &                       - dg%iota2min(nm(k,lll),ll) )/(dg%iota2taylorvert(k,ll,lll) - dg%iota2max(nm(k,lll),ll)) ) ) 
                        

                     elseif ( (dg%iota2taylorvert(k,ll,lll).lt.dg%iota2taylor(k,ll,1) )
     &                       .and.( abs(dg%iota2taylorvert(k,ll,lll)-dg%iota2taylor(k,ll,1)).gt.1.0E-15 ).and.
     &                       ( dg%iota2min(nm(k,lll),ll).ne.dg%iota2taylor(k,ll,1) ) ) then 

                        dg%alphaiota20(k,ll,lll) = min( 1.D0,( dg%iota2min(nm(k,lll),ll)
     &                       - dg%iota2taylor(k,ll,1) )/( dg%iota2taylorvert(k,ll,lll)-dg%iota2taylor(k,ll,1)))

                                !adapted part

                     elseif ( (dg%iota2taylorvert(k,ll,lll).lt.dg%iota2taylor(k,ll,1) )
     &                       .and.( abs(dg%iota2taylorvert(k,ll,lll)-dg%iota2taylor(k,ll,1)).gt.1.0E-15 ).and.
     &                       ( dg%iota2min(nm(k,lll),ll).eq.dg%iota2taylor(k,ll,1) ) ) then 

                        dg%alphaiota20(k,ll,lll) = min( fd,abs( ( dg%iota2min(nm(k,lll),ll)
     &                       - dg%iota2max(nm(k,lll),ll) )/( dg%iota2taylorvert(k,ll,lll)- dg%iota2min(nm(k,lll),ll))) )
                        
                        
                     elseif ( ( dg%iota2taylorvert(k,ll,lll).eq.dg%iota2taylor(k,ll,1) ).or.
     &                       ( abs(dg%iota2taylorvert(k,ll,lll)-dg%iota2taylor(k,ll,1)).le.1.0E-15 ) ) then

                        dg%alphaiota20(k,ll,lll) = 1.D0

                     endif

                  else

                     dg%alphaiota20(k,ll,lll) = 1.D0

                  endif
#endif

                  if (dg%QXmin(nm(k,lll),ll).ne.dg%QXmax(nm(k,lll),ll)) then

                     if ( ( dg%QXtaylorvert(k,ll,lll).gt.dg%QXtaylor(k,ll,1) ).and.
     &                    ( abs(dg%QXtaylorvert(k,ll,lll)-dg%QXtaylor(k,ll,1)).gt.1.0E-15 ).and.
     &                    ( dg%QXmax(nm(k,lll),ll).ne.dg%QXtaylor(k,ll,1) ) ) then  

                        dg%alphaQX0(k,ll,lll) = min(1.D0,  ( dg%QXmax(nm(k,lll),ll)
     &                       - dg%QXtaylor(k,ll,1) )/ (dg%QXtaylorvert(k,ll,lll) - dg%QXtaylor(k,ll,1)))

                                !adapted part

                     elseif ( ( dg%QXtaylorvert(k,ll,lll).gt.dg%QXtaylor(k,ll,1) ).and.
     &                       ( abs(dg%QXtaylorvert(k,ll,lll)-dg%QXtaylor(k,ll,1)).gt.1.0E-15 ).and.
     &                       ( dg%QXmax(nm(k,lll),ll).eq.dg%QXtaylor(k,ll,1) ) ) then 

                        dg%alphaQX0(k,ll,lll) = min(fd, abs( ( dg%QXmax(nm(k,lll),ll)
     &                       - dg%QXmin(nm(k,lll),ll) )/(dg%QXtaylorvert(k,ll,lll) - dg%QXmax(nm(k,lll),ll)) ) ) 
                        

                     elseif ( (dg%QXtaylorvert(k,ll,lll).lt.dg%QXtaylor(k,ll,1) )
     &                       .and.( abs(dg%QXtaylorvert(k,ll,lll)-dg%QXtaylor(k,ll,1)).gt.1.0E-15 ).and.
     &                       ( dg%QXmin(nm(k,lll),ll).ne.dg%QXtaylor(k,ll,1) ) ) then 

                        dg%alphaQX0(k,ll,lll) = min( 1.D0,( dg%QXmin(nm(k,lll),ll)
     &                       - dg%QXtaylor(k,ll,1) )/( dg%QXtaylorvert(k,ll,lll)-dg%QXtaylor(k,ll,1)))

                                !adapted part

                     elseif ( (dg%QXtaylorvert(k,ll,lll).lt.dg%QXtaylor(k,ll,1) )
     &                       .and.( abs(dg%QXtaylorvert(k,ll,lll)-dg%QXtaylor(k,ll,1)).gt.1.0E-15 ).and.
     &                       ( dg%QXmin(nm(k,lll),ll).eq.dg%QXtaylor(k,ll,1) ) ) then 

                        dg%alphaQX0(k,ll,lll) = min( fd,abs( ( dg%QXmin(nm(k,lll),ll)
     &                       - dg%QXmax(nm(k,lll),ll) )/( dg%QXtaylorvert(k,ll,lll)- dg%QXmin(nm(k,lll),ll))) )
                        
                        
                     elseif ( ( dg%QXtaylorvert(k,ll,lll).eq.dg%QXtaylor(k,ll,1) ).or.
     &                       ( abs(dg%QXtaylorvert(k,ll,lll)-dg%QXtaylor(k,ll,1)).le.1.0E-15 ) ) then

                        dg%alphaQX0(k,ll,lll) = 1.D0

                     endif

                  else

                     dg%alphaQX0(k,ll,lll) = 1.D0

                  endif


                  if (dg%QYmin(nm(k,lll),ll).ne.dg%QYmax(nm(k,lll),ll)) then

                     if ( ( dg%QYtaylorvert(k,ll,lll).gt.dg%QYtaylor(k,ll,1) ).and.
     &                    ( abs(dg%QYtaylorvert(k,ll,lll)-dg%QYtaylor(k,ll,1)).gt.1.0E-15 ).and.
     &                    ( dg%QYmax(nm(k,lll),ll).ne.dg%QYtaylor(k,ll,1) ) ) then  

                        dg%alphaQY0(k,ll,lll) = min(1.D0,  ( dg%QYmax(nm(k,lll),ll)
     &                       - dg%QYtaylor(k,ll,1) )/ (dg%QYtaylorvert(k,ll,lll) - dg%QYtaylor(k,ll,1)))

                                !adapted part

                     elseif ( ( dg%QYtaylorvert(k,ll,lll).gt.dg%QYtaylor(k,ll,1) ).and.
     &                       ( abs(dg%QYtaylorvert(k,ll,lll)-dg%QYtaylor(k,ll,1)).gt.1.0E-15 ).and.
     &                       ( dg%QYmax(nm(k,lll),ll).eq.dg%QYtaylor(k,ll,1) ) ) then 

                        dg%alphaQY0(k,ll,lll) = min(fd, abs( ( dg%QYmax(nm(k,lll),ll)
     &                       - dg%QYmin(nm(k,lll),ll) )/(dg%QYtaylorvert(k,ll,lll) - dg%QYmax(nm(k,lll),ll)) ) ) 
                        

                     elseif ( (dg%QYtaylorvert(k,ll,lll).lt.dg%QYtaylor(k,ll,1) )
     &                       .and.( abs(dg%QYtaylorvert(k,ll,lll)-dg%QYtaylor(k,ll,1)).gt.1.0E-15 ).and.
     &                       ( dg%QYmin(nm(k,lll),ll).ne.dg%QYtaylor(k,ll,1) ) ) then 

                        dg%alphaQY0(k,ll,lll) = min( 1.D0,( dg%QYmin(nm(k,lll),ll)
     &                       - dg%QYtaylor(k,ll,1) )/( dg%QYtaylorvert(k,ll,lll)-dg%QYtaylor(k,ll,1)))

                                !adapted part

                     elseif ( (dg%QYtaylorvert(k,ll,lll).lt.dg%QYtaylor(k,ll,1) )
     &                       .and.( abs(dg%QYtaylorvert(k,ll,lll)-dg%QYtaylor(k,ll,1)).gt.1.0E-15 ).and.
     &                       ( dg%QYmin(nm(k,lll),ll).eq.dg%QYtaylor(k,ll,1) ) ) then 

                        dg%alphaQY0(k,ll,lll) = min( fd,abs( ( dg%QYmin(nm(k,lll),ll)
     &                       - dg%QYmax(nm(k,lll),ll) )/( dg%QYtaylorvert(k,ll,lll)- dg%QYmin(nm(k,lll),ll))) )
                        
                        
                     elseif ( ( dg%QYtaylorvert(k,ll,lll).eq.dg%QYtaylor(k,ll,1) ).or.
     &                       ( abs(dg%QYtaylorvert(k,ll,lll)-dg%QYtaylor(k,ll,1)).le.1.0E-15 ) ) then

                        dg%alphaQY0(k,ll,lll) = 1.D0

                     endif

                  else

                     dg%alphaQY0(k,ll,lll) = 1.D0

                  endif

               enddo            !lll

            enddo               !ll

         elseif (dg%dofs(k).eq.1) then

            cycle

         endif

      enddo

!.....Find the prescribed higher limiters by finding smallest local value

      dg%alphaZE = 0.D0
      dg%alphaQX = 0.D0
      dg%alphaQY = 0.D0

#ifdef TRACE
      dg%alphaiota = 0.D0
#endif

#ifdef CHEM
      dg%alphaiota = 0.D0
      dg%alphaiota2 = 0.D0
#endif

      do k=1,mne

         if (dg%dofs(k).gt.1) then

            do ll=1,dg%dofs(k)		

               dg%alphaZE(k,ll) = minval( dg%alphaZE0(k,ll,:) )
               dg%alphaQX(k,ll) = minval( dg%alphaQX0(k,ll,:) )
               dg%alphaQY(k,ll) = minval( dg%alphaQY0(k,ll,:) )

#ifdef TRACE
               dg%alphaiota(k,ll) = minval( dg%alphaiota0(k,ll,:) )
#endif

#ifdef CHEM
               dg%alphaiota(k,ll) = minval( dg%alphaiota0(k,ll,:) )
               dg%alphaiota2(k,ll) = minval( dg%alphaiota20(k,ll,:) )
#endif

            enddo

         elseif (dg%dofs(k).eq.1) then

            cycle

         endif

      enddo

!.... Choose smallest (minimum) alpha for derivative in x or y

      dg%alphaZEm = 0.D0
      dg%alphaQXm = 0.D0
      dg%alphaQYm = 0.D0

#ifdef TRACE
      dg%alphaiotam = 0.D0
#endif

#ifdef CHEM
      dg%alphaiotam = 0.D0
      dg%alphaiota2m = 0.D0
#endif

      do k=1,mne

         if (dg%dofs(k).gt.1) then
            
            do bb = 1,pdg_el(k)

               if( (bb+1)*(bb+2)/2.le.dg%dofs(k) ) then

                  dg%alphaZEm(k,(bb*(bb+1))/2 + 1) = 
     &                 minval( dg%alphaZE(k,(bb*(bb-1))/2 +1:(bb*(bb+1)) / 2  ) )
                  dg%alphaQXm(k,(bb*(bb+1))/2 + 1) = 
     &                 minval( dg%alphaQX(k,(bb*(bb-1))/2 +1:(bb*(bb+1)) / 2  ) )
                  dg%alphaQYm(k,(bb*(bb+1))/2 + 1) = 
     &                 minval( dg%alphaQY(k,(bb*(bb-1))/2 +1:(bb*(bb+1)) / 2  ) )

#ifdef TRACE
                  dg%alphaiotam(k,(bb*(bb+1))/2 + 1) = 
     &                 minval( dg%alphaiota(k,(bb*(bb-1))/2 +1:(bb*(bb+1)) / 2 ) )
#endif

#ifdef CHEM
                  dg%alphaiotam(k,(bb*(bb+1))/2 + 1) = 
     &                 minval( dg%alphaiota(k,(bb*(bb-1))/2 +1:(bb*(bb+1)) / 2  ) )
                  dg%alphaiota2m(k,(bb*(bb+1))/2 + 1) = 
     &                 minval( dg%alphaiota2(k,(bb*(bb-1))/2 +1:(bb*(bb+1)) / 2 ) )
#endif

               endif
               
            enddo

         elseif (dg%dofs(k).eq.1) then

            cycle

         endif

      enddo

!.....Use max higher derivative values for final limiter value

      dg%alphaZE_max = 0.D0
      dg%alphaQX_max = 0.D0
      dg%alphaQY_max = 0.D0

#ifdef TRACE
      dg%alphaiota_max = 0.D0
#endif

#ifdef CHEM
      dg%alphaiota_max = 0.D0
      dg%alphaiota2_max = 0.D0
#endif

      do k=1,mne

         if (dg%dofs(k).gt.1) then

            do bb =1,pdg_el(k)

               if( (bb+1)*(bb+2)/2.le.dg%dofs(k)) then

                  dg%alphaZE_max(k,(bb*(bb+1))/2 + 1) = fde*maxval( 
     &                 dg%alphaZEm(k,(bb*(bb+1))/2 + 1:dg%dofs(k)) )
                  dg%alphaQX_max(k,(bb*(bb+1))/2 + 1) = fde*maxval( 
     &                 dg%alphaQXm(k,(bb*(bb+1))/2 + 1:dg%dofs(k)) )
                  dg%alphaQY_max(k,(bb*(bb+1))/2 + 1) = fde*maxval( 
     &                 dg%alphaQYm(k,(bb*(bb+1))/2 + 1:dg%dofs(k)) )

#ifdef TRACE
                  dg%alphaiota_max(k,(bb*(bb+1))/2 + 1) = fde*maxval( 
     &                 dg%alphaiotam(k,(bb*(bb+1)/2) +1 : dg%dofs(k)))
#endif

#ifdef CHEM
                  dg%alphaiota_max(k,(bb*(bb+1))/2 + 1) = fde*maxval( 
     &                 dg%alphaiotam(k,(bb*(bb+1))/2 + 1:dg%dofs(k)) )
                  dg%alphaiota2_max(k,(bb*(bb+1))/2 + 1) = fde*maxval( 
     &                 dg%alphaiota2m(k,(bb*(bb+1))/2 + 1:dg%dofs(k)) )
#endif

               endif

            enddo

         elseif (dg%dofs(k).eq.1) then

            cycle

         endif

      enddo

!.....Limit on the Master element in the Taylor basis, via reconstruction 
!.....of unconstrained solutions with alpha constraints


      dg%limitZE = 0.D0
      dg%limitQX = 0.D0
      dg%limitQY = 0.D0

      dg%lim_count_roll = 0

#ifdef TRACE
      dg%limitiota = 0.D0
#endif

#ifdef CHEM
      dg%limitiota = 0.D0
      dg%limitiota2 = 0.D0
#endif

      do k=1,mne

         dg%lim_count = 0

         if (dg%dofs(k).gt.1) then

            do ll=1,dg%dofs(k)

               if ( ll.eq.1 ) then

                  dg%limitZE(k,ll) = dg%ZEtaylor(k,ll,1)
                  dg%limitQX(k,ll) = dg%QXtaylor(k,ll,1)
                  dg%limitQY(k,ll) = dg%QYtaylor(k,ll,1) 

#ifdef TRACE
                  dg%limitiota(k,ll) = dg%iotataylor(k,ll,1)
#endif

#ifdef CHEM
                  dg%limitiota(k,ll) = dg%iotataylor(k,ll,1)
                  dg%limitiota2(k,ll) = dg%iota2taylor(k,ll,1)
#endif

               elseif ( ll.ge.2 ) then
                  
                  do bb=1,pdg_el(k)

                     if ( ll.le.( ( (bb+1)*(bb+2)) / 2.D0 ).and.(ll.gt.
     &                    (bb*(bb+1)/2.D0) ) ) then

                        dg%limitZE(k,ll) = dg%alphaZE_max(k,(bb*(bb+1))/2 + 1) 
     &                       * dg%ZEtaylor(k,ll,1)
                        dg%limitQX(k,ll) = dg%alphaQX_max(k,(bb*(bb+1))/2 + 1) 
     &                       * dg%ZEtaylor(k,ll,1)
                        dg%limitQY(k,ll) = dg%alphaQY_max(k,(bb*(bb+1))/2 + 1) 
     &                       * dg%ZEtaylor(k,ll,1)


#ifdef TRACE
                        dg%limitiota(k,ll) = dg%alphaiota_max(k,(bb*(bb+1))/2 + 1) 
     &                       * dg%iotataylor(k,ll,1)
#endif

#ifdef CHEM
                        dg%limitiota(k,ll) = dg%alphaiota_max(k,(bb*(bb+1))/2 + 1) 
     &                       * dg%iotataylor(k,ll,1)
                        dg%limitiota2(k,ll) = dg%alphaiota2_max(k,(bb*(bb+1))/2 + 1) 
     &                       * dg%iota2taylor(k,ll,1)
#endif


!$$$  ! Make a counter to track limiting
!$$$  
!$$$  if ( ( dg%alphaZE_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
!$$$  &                       dg%alphaQX_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
!$$$  &                       dg%alphaQY_max(k,(bb*(bb+1))/2 + 1).lt.1.D0 ).and.
!$$$  &                       chem_flag.ne.1.and.tracer_flag.ne.1 ) then
!$$$  
!$$$  dg%lim_count = 1  
!$$$  
!$$$  elseif ( (dg%alphaZE_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
!$$$  &                          dg%alphaQX_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
!$$$  &                          dg%alphaQY_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
!$$$  &                          dg%alphaiota_max(k,(bb*(bb+1))/2 + 1).lt.1.D0 ).and.
!$$$  &                          chem_flag.ne.1.and.tracer_flag.eq.1 ) then
!$$$  
!$$$  dg%lim_count = 1   
!$$$  
!$$$  elseif ( (dg%alphaZE_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
!$$$  &                          dg%alphaQX_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
!$$$  &                          dg%alphaQY_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
!$$$  &                          dg%alphaiota_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
!$$$  &                          dg%alphaiota2_max(k,(bb*(bb+1))/2 + 1).lt.1.D0 ).and.
!$$$  &                          chem_flag.eq.1 ) then
!$$$  
!$$$  dg%lim_count = 1 
!$$$  
!$$$  endif
                        
                     endif
                     
                  enddo
                  
               endif

            enddo

         elseif (dg%dofs(k).eq.1) then

            cycle

         endif

                                !dg%lim_count_roll = dg%lim_count_roll + dg%lim_count

      enddo

!.....Transform back to the Dubiner basis (multiply by dg%NmatrixInv),

      dg%ZEconst =  0.D0
      dg%QXconst =  0.D0
      dg%QYconst =  0.D0

#ifdef TRACE
      dg%iotaconst =  0.D0
#endif

#ifdef CHEM
      dg%iotaconst =  0.D0
      dg%iota2const =  0.D0
#endif

      do k=1,mne

         if (dg%dofs(k).gt.1) then
                                !do lll=1,3

            do ll=1,dg%dofs(k)

               do ss=1,dg%dofs(k)

                  dg%ZEconst(k,ll) = dg%ZEconst(k,ll) + dg%NmatrixInv(k,ll,ss,dg%dofs(k)) 
     &                 * dg%limitZE(k,ss)
                  dg%QXconst(k,ll) = dg%QXconst(k,ll) + dg%NmatrixInv(k,ll,ss,dg%dofs(k)) 
     &                 * dg%limitQX(k,ss)
                  dg%QYconst(k,ll) = dg%QYconst(k,ll) + dg%NmatrixInv(k,ll,ss,dg%dofs(k)) 
     &                 * dg%limitQY(k,ss)

#ifdef TRACE
                  dg%iotaconst(k,ll) = dg%iotaconst(k,ll) + 
     &                 dg%NmatrixInv(k,ll,ss,dg%dofs(k)) * dg%limitiota(k,ss)
#endif

#ifdef CHEM
                  dg%iotaconst(k,ll) = dg%iotaconst(k,ll) + 
     &                 dg%NmatrixInv(k,ll,ss,dg%dofs(k)) * dg%limitiota(k,ss)
                  dg%iota2const(k,ll) = dg%iota2const(k,ll) + 
     &                 dg%NmatrixInv(k,ll,ss,dg%dofs(k)) * dg%limitiota2(k,ss)
#endif

               enddo

            enddo

         elseif (dg%dofs(k).eq.1) then

            cycle

         endif

      enddo

!.....Set limit values

      do k=1,mne

         if (dg%dofs(k).gt.1) then

            do ll = 1,dg%dofs(k)

               dg%ZE(ll,k,dg%irk+1) = dg%ZEconst(k,ll)
               dg%QX(ll,k,dg%irk+1) = dg%QXconst(k,ll)
               dg%QY(ll,k,dg%irk+1) = dg%QYconst(k,ll)

#ifdef TRACE
               dg%iota(ll,k,dg%irk+1) = dg%iotaconst(k,ll)
#endif

#ifdef CHEM
               dg%iota(ll,k,dg%irk+1) = dg%iotaconst(k,ll)
               dg%iota2(ll,k,dg%irk+1) = dg%iota2const(k,ll)
#endif

            enddo

         elseif (dg%dofs(k).eq.1) then

            cycle

         endif

      enddo
      
      return
      end subroutine

!******************************************************************************
!     
!     SUBROUTINE SLOPELIMITER9() -- dg%Needs to be updated for dg%ZE,dg%QX,dg%QY and dg%iota2
!     
!     Written 2011
!     
!     This subroutine selects the first adapted Barth--Jespersen limiter over
!     a Taylor Polynomial basis, and is consistent with p-adaptation
!     to arbitrary order p
!     
!     - cem
!     
!******************************************************************************

      SUBROUTINE SLOPELIMITER9()

!.....Use appropriate modules

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

      Allocate ( ZEminel(mne,dg%dofh),ZEmaxel(mne,dg%dofh),QXminel(mne,dg%dofh) )
      Allocate ( QYminel(mne,dg%dofh),QYmaxel(mne,dg%dofh),QXmaxel(mne,dg%dofh) )
      Allocate ( iotaminel(mne,dg%dofh),iotamaxel(mne,dg%dofh) )
      Allocate ( iota2minel(mne,dg%dofh),iota2maxel(mne,dg%dofh) )


!.....We work over the master element
!.....Set initial values

      fd = dg%slope_weight         ! add weight for lower order pieces (fd<1 => stronger limiting)
      

      DO k=1,NE

         if (dg%dofs(k).gt.1) then

            DO ll = 1,dg%dofs(k)

               dg%ZE(ll,k,dg%NRK+2) = dg%ZE(ll,k,dg%IRK+1)
               dg%QX(ll,k,dg%NRK+2) = dg%QX(ll,k,dg%IRK+1)
               dg%QY(ll,k,dg%NRK+2) = dg%QY(ll,k,dg%IRK+1)

#ifdef TRACE
               dg%iota(ll,k,dg%NRK+2) = dg%iota(ll,k,dg%IRK+1)
#endif

#ifdef CHEM
               dg%iota(ll,k,dg%NRK+2) = dg%iota(ll,k,dg%IRK+1)
               dg%iota2(ll,k,dg%NRK+2) = dg%iota2(ll,k,dg%IRK+1)
#endif


            ENDDO

         elseif (dg%dofs(k).eq.1) then

            cycle

         endif

      ENDDO   

      marea = 2.D0              !master elements area

!.....Convert initial values to the Taylor basis (multiply by dg%Nmatrix) on base element


      dg%ZEtaylor = 0.D0 
      dg%QXtaylor = 0.D0
      dg%QYtaylor = 0.D0

#ifdef TRACE
      dg%iotataylor = 0.D0
#endif

#ifdef CHEM
      dg%iotataylor = 0.D0
      dg%iota2taylor = 0.D0
#endif

      do k=1,ne

         if (dg%dofs(k).gt.1) then

            do ll=1,dg%dofs(k)
               
               do ss=1,dg%dofs(k)

                  dg%ZEtaylor(k,ll,1) = dg%ZEtaylor(k,ll,1) + dg%Nmatrix(k,ll,ss,dg%dofs(k)) * dg%ZE(ss,k,dg%nrk+2)
                  dg%QXtaylor(k,ll,1) = dg%QXtaylor(k,ll,1) + dg%Nmatrix(k,ll,ss,dg%dofs(k)) * dg%QX(ss,k,dg%nrk+2)
                  dg%QYtaylor(k,ll,1) = dg%QYtaylor(k,ll,1) + dg%Nmatrix(k,ll,ss,dg%dofs(k)) * dg%QY(ss,k,dg%nrk+2)

#ifdef TRACE
                  dg%iotataylor(k,ll,1) = dg%iotataylor(k,ll,1) + dg%Nmatrix(k,ll,ss,dg%dofs(k))*dg%iota(ss,k,dg%nrk+2)
#endif

#ifdef CHEM
                  dg%iotataylor(k,ll,1) = dg%iotataylor(k,ll,1) + dg%Nmatrix(k,ll,ss,dg%dofs(k))*dg%iota(ss,k,dg%nrk+2)
                  dg%iota2taylor(k,ll,1) = dg%iota2taylor(k,ll,1) + dg%Nmatrix(k,ll,ss,dg%dofs(k))*dg%iota2(ss,k,dg%nrk+2)
#endif
                  
               enddo

            enddo

         elseif (dg%dofs(k).eq.1) then

            cycle

         endif

      enddo

!.....Find values at vertices of base elements and neighbors


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

      CALL UPDATER_ELEM_MOD2(dg%ZEtaylor,dg%QXtaylor,dg%QYtaylor,1,3)

#ifdef TRACE
      CALL UPDATER_ELEM_MOD2(dg%iotataylor,dg%QXtaylor,dg%QYtaylor,1,2 )  
#endif

#ifdef CHEM
      CALL UPDATER_ELEM_MOD2(dg%iotataylor,dg%iota2taylor,dg%QYtaylor,1,2 )  
#endif

#endif

      do k = 1,ne

         do ll = 1,dg%dofs(k) 

            do ell = 1,3        ! Number of edge neighbors for a triangle

!.....Find max and min values over polynomial coefficients
               
               ZEmaxel(k,ll) = max( dg%ZEtaylor(k,ll,1),dg%ZEtaylor(dg%EL_NBORS(ell,k),ll,1), ZEmaxel(k,ll) )
               QXmaxel(k,ll) = max( dg%QXtaylor(k,ll,1),dg%QXtaylor(dg%EL_NBORS(ell,k),ll,1), QXmaxel(k,ll) )
               QYmaxel(k,ll) = max( dg%QYtaylor(k,ll,1),dg%QYtaylor(dg%EL_NBORS(ell,k),ll,1), QYmaxel(k,ll) )

               ZEminel(k,ll) = min( dg%ZEtaylor(k,ll,1),dg%ZEtaylor(dg%EL_NBORS(ell,k),ll,1), ZEminel(k,ll) )
               QXminel(k,ll) = min( dg%QXtaylor(k,ll,1),dg%QXtaylor(dg%EL_NBORS(ell,k),ll,1), QXminel(k,ll) )
               QYminel(k,ll) = min( dg%QYtaylor(k,ll,1),dg%QYtaylor(dg%EL_NBORS(ell,k),ll,1), QYminel(k,ll) )

#ifdef TRACE
               iotamaxel(k,ll) = max( dg%iotataylor(k,ll,1),dg%iotataylor(dg%EL_NBORS(ell,k),ll,1), 
     &              iotamaxel(k,ll) )
               iotaminel(k,ll) = min( dg%iotataylor(k,ll,1),dg%iotataylor(dg%EL_NBORS(ell,k),ll,1), 
     &              iotaminel(k,ll) )
#endif

#ifdef CHEM
               iotamaxel(k,ll) = max( dg%iotataylor(k,ll,1),dg%iotataylor(dg%EL_NBORS(ell,k),ll,1), 
     &              iotamaxel(k,ll) )
               iota2maxel(k,ll) = max( dg%iota2taylor(k,ll,1),dg%iota2taylor(dg%EL_NBORS(ell,k),ll,1),
     &              iota2maxel(k,ll) )
               iotaminel(k,ll) = min( dg%iotataylor(k,ll,1),dg%iotataylor(dg%EL_NBORS(ell,k),ll,1),
     &              iotaminel(k,ll) )
               iota2minel(k,ll) = min( dg%iota2taylor(k,ll,1),dg%iota2taylor(dg%EL_NBORS(ell,k),ll,1),
     &              iota2minel(k,ll) )
#endif
               
            enddo
            
         enddo

      enddo



!.....Must generate linear recostructions at vertices



      dg%ZEtaylorvert = 0.D0
      dg%QXtaylorvert = 0.D0
      dg%Qytaylorvert = 0.D0

#ifdef TRACE
      dg%iotataylorvert = 0.D0
#endif

#ifdef CHEM
      dg%iotataylorvert = 0.D0
      dg%iota2taylorvert = 0.D0
#endif

      do k=1,ne

         if (dg%dofs(k).gt.1) then

            do ll=1,dg%dofs(k)

               do lll=1,3

                  if (ll.eq.1) then

                     dg%ZEtaylorvert(k,ll,lll) = dg%ZEtaylor(k,ll,1) + dg%ZEtaylor(k,ll+1,1)*( dg%xi2vert(k,lll) -dg%xi2BCb(k) )
     &                    + dg%ZEtaylor(k,ll+2,1)*( dg%xi1vert(k,lll) - dg%xi1BCb(k) )
                     dg%QXtaylorvert(k,ll,lll) = dg%QXtaylor(k,ll,1) + dg%QXtaylor(k,ll+1,1)*( dg%xi2vert(k,lll) -dg%xi2BCb(k) )
     &                    + dg%QXtaylor(k,ll+2,1)*( dg%xi1vert(k,lll) - dg%xi1BCb(k) )
                     dg%Qytaylorvert(k,ll,lll) = dg%QYtaylor(k,ll,1) + dg%QYtaylor(k,ll+1,1)*( dg%xi2vert(k,lll) -dg%xi2BCb(k) ) 
     &                    + dg%QYtaylor(k,ll+2,1)*( dg%xi1vert(k,lll) - dg%xi1BCb(k) )

#ifdef TRACE
                     dg%iotataylorvert(k,ll,lll) = dg%iotataylor(k,ll,1) + 
     &                    dg%iotataylor(k,ll+1,1)*( dg%xi2vert(k,lll) - dg%xi2BCb(k) )
     &                    + dg%iotataylor(k,ll+2,1)*( dg%xi1vert(k,lll) - dg%xi1BCb(k) )
#endif

#ifdef CHEM
                     dg%iotataylorvert(k,ll,lll) = dg%iotataylor(k,ll,1) + 
     &                    dg%iotataylor(k,ll+1,1)*( dg%xi2vert(k,lll) - dg%xi2BCb(k) )
     &                    + dg%iotataylor(k,ll+2,1)*( dg%xi1vert(k,lll) - dg%xi1BCb(k) )
                     dg%iota2taylorvert(k,ll,lll) = dg%iota2taylor(k,ll,1) + 
     &                    dg%iota2taylor(k,ll+1,1)*( dg%xi2vert(k,lll) - dg%xi2BCb(k) )
     &                    + dg%iota2taylor(k,ll+2,1)*( dg%xi1vert(k,lll) - dg%xi1BCb(k) )
#endif

                  elseif (ll.gt.1) then

                     dg%ZEtaylorvert(k,ll,lll) = dg%ZEtaylor(k,ll,1) + 
     &                    dg%ZEtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) ),1)*( dg%xi2vert(k,lll) - dg%xi2BCb(k) ) 
     &                    + dg%ZEtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) )+1,1)*( dg%xi1vert(k,lll) - dg%xi1BCb(k) )
                     dg%QXtaylorvert(k,ll,lll) = dg%QXtaylor(k,ll,1) + 
     &                    dg%QXtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) ),1)*( dg%xi2vert(k,lll) - dg%xi2BCb(k) ) 
     &                    + dg%QXtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) )+1,1)*( dg%xi1vert(k,lll) - dg%xi1BCb(k) )
                     dg%Qytaylorvert(k,ll,lll) = dg%QYtaylor(k,ll,1) + 
     &                    dg%QYtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) ),1)*( dg%xi2vert(k,lll) - dg%xi2BCb(k) ) 
     &                    + dg%QYtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) )+1,1)*( dg%xi1vert(k,lll) - dg%xi1BCb(k) )

#ifdef TRACE
                     dg%iotataylorvert(k,ll,lll) = dg%iotataylor(k,ll,1)+
     &                    dg%iotataylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) ),1)*( dg%xi2vert(k,lll) - dg%xi2BCb(k))
     &                    + dg%iotataylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) )+1,1)*( dg%xi1vert(k,lll) - dg%xi1BCb(k) )
#endif

#ifdef CHEM
                     dg%iotataylorvert(k,ll,lll) = dg%iotataylor(k,ll,1)+
     &                    dg%iotataylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) ),1)*( dg%xi2vert(k,lll) - dg%xi2BCb(k))
     &                    + dg%iotataylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) )+1,1)*( dg%xi1vert(k,lll) - dg%xi1BCb(k) )
                     dg%iota2taylorvert(k,ll,lll) = dg%iota2taylor(k,ll,1)+
     &                    dg%iotataylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) ),1)*( dg%xi2vert(k,lll) - dg%xi2BCb(k))
     &                    + dg%iota2taylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) )+1,1)*( dg%xi1vert(k,lll) - dg%xi1BCb(k) )
#endif

                  else

                  endif

               enddo
               
            enddo
            
         elseif (dg%dofs(k).eq.1) then

            cycle

         endif

      enddo

!.....Compute alphas for each variable in each order derivitive


      dg%alphaZE0 = 0.D0
      dg%alphaQX0 = 0.D0
      dg%alphaQY0 = 0.D0

#ifdef TRACE
      dg%alphaiota0 = 0.D0
#endif

#ifdef CHEM
      dg%alphaiota0 = 0.D0
      dg%alphaiota20 = 0.D0
#endif

      do k = 1,ne

         if (dg%dofs(k).gt.1) then

            do ll = 1,dg%dofs(k)

               do lll = 1,3

                  if ( dg%ZEtaylorvert(k,ll,lll).gt.dg%ZEtaylor(k,ll,1).and.
     &                 abs((dg%ZEtaylorvert(k,ll,lll)-dg%ZEtaylor(k,ll,1))).gt.1.0E-15 ) then 

                     dg%alphaZE0(k,ll,lll) = min( 1.D0, ( ZEmaxel(k,ll) - 
     &                    dg%ZEtaylor(k,ll,1) ) / (  dg%ZEtaylorvert(k,ll,lll) - dg%ZEtaylor(k,ll,1) ) )
                     
                  elseif ( dg%ZEtaylorvert(k,ll,lll).eq.dg%ZEtaylor(k,ll,1).or.
     &                    abs((dg%ZEtaylorvert(k,ll,lll)-dg%ZEtaylor(k,ll,1))).le.1.0E-15 ) then

                     dg%alphaZE0(k,ll,lll) = 1.D0

                  elseif (  dg%ZEtaylorvert(k,ll,lll).lt.dg%ZEtaylor(k,ll,1).and.
     &                    abs((dg%ZEtaylorvert(k,ll,lll)-dg%ZEtaylor(k,ll,1))).gt.1.0E-15 ) then

                     dg%alphaZE0(k,ll,lll) = min( 1.D0, ( ZEminel(k,ll)
     &                    - dg%ZEtaylor(k,ll,1) ) / ( dg%ZEtaylorvert(k,ll,lll) - dg%ZEtaylor(k,ll,1) ) )
                     
                  endif

#ifdef TRACE
                  if (iotaminel(k,ll).ne.iotamaxel(k,ll)) then
                     
                     if ( dg%iotataylorvert(k,ll,lll).gt.dg%iotataylor(k,ll,1).and.
     &                    abs((dg%iotataylorvert(k,ll,lll)-dg%iotataylor(k,ll,1))).gt.1.0E-15 ) then 

                        dg%alphaiota0(k,ll,lll) = min(1.D0,( iotamaxel(k,ll)
     &                       -dg%iotataylor(k,ll,1) )/ (dg%iotataylorvert(k,ll,lll) - dg%iotataylor(k,ll,1)))

                                !adapted part

                     elseif ( ( dg%iotataylorvert(k,ll,lll).gt.dg%iotataylor(k,ll,1) ).and.
     &                       ( abs(dg%iotataylorvert(k,ll,lll)-dg%iotataylor(k,ll,1)).gt.1.0E-15 ).and.
     &                       ( iotamaxel(k,ll).eq.dg%iotataylor(k,ll,1) ) ) then 

                        dg%alphaiota0(k,ll,lll) = min(fd, abs( ( iotamaxel(k,ll)
     &                       - iotaminel(k,ll) )/(dg%iotataylorvert(k,ll,lll) - iotamaxel(k,ll)) ) ) 

                                !adapted part

                     elseif ( ( dg%iotataylorvert(k,ll,lll).gt.dg%iotataylor(k,ll,1) ).and.
     &                       ( abs(dg%iotataylorvert(k,ll,lll)-dg%iotataylor(k,ll,1)).gt.1.0E-15 ).and.
     &                       ( iotaminel(k,ll).eq.dg%iotataylor(k,ll,1) ) ) then 

                        dg%alphaiota0(k,ll,lll) = min(fd, abs( ( iotaminel(k,ll)
     &                       - iotamaxel(k,ll) )/(dg%iotataylorvert(k,ll,lll) - iotaminel(k,ll)) ) )
                        
                     elseif (dg%iotataylorvert(k,ll,lll).eq.dg%iotataylor(k,ll,1).or.
     &                       abs((dg%iotataylorvert(k,ll,lll)-dg%iotataylor(k,ll,1))).le.1.0E-15 ) then

                        dg%alphaiota0(k,ll,lll) = 1.D0

                     elseif (dg%iotataylorvert(k,ll,lll).lt.dg%iotataylor(k,ll,1).
     &                       and.abs((dg%iotataylorvert(k,ll,lll)-dg%iotataylor(k,ll,1))).gt.1.0E-15) then

                        dg%alphaiota0(k,ll,lll) = min( 1.D0,( iotaminel(k,ll)
     &                       -dg%iotataylor(k,ll,1) )/( dg%iotataylorvert(k,ll,lll)-dg%iotataylor(k,ll,1)))
                        
                     endif

                  else
                     
                     dg%alphaiota0(k,ll,lll) = 1.D0

                  endif
#endif

#ifdef CHEM                  
                  if ( dg%iotataylorvert(k,ll,lll).gt.dg%iotataylor(k,ll,1).and.
     &                 abs((dg%iotataylorvert(k,ll,lll)-dg%iotataylor(k,ll,1))).gt.1.0E-15 ) then 

                     dg%alphaiota0(k,ll,lll) = min(1.D0,( iotamaxel(k,ll)-
     &                    dg%iotataylor(k,ll,1) )/ (dg%iotataylorvert(k,ll,lll) - dg%iotataylor(k,ll,1)))
                     
                  elseif (dg%iotataylorvert(k,ll,lll).eq.dg%iotataylor(k,ll,1).or.
     &                    abs((dg%iotataylorvert(k,ll,lll)-dg%iotataylor(k,ll,1))).le.1.0E-15 ) then

                     dg%alphaiota0(k,ll,lll) = 1.D0

                  elseif (dg%iotataylorvert(k,ll,lll).lt.dg%iotataylor(k,ll,1).and.
     &                    abs((dg%iotataylorvert(k,ll,lll)-dg%iotataylor(k,ll,1))).gt.1.0E-15) then

                     dg%alphaiota0(k,ll,lll) = min( 1.D0,( iotaminel(k,ll)
     &                    -dg%iotataylor(k,ll,1) )/( dg%iotataylorvert(k,ll,lll)-dg%iotataylor(k,ll,1)))
                     
                  endif

                  if ( dg%iota2taylorvert(k,ll,lll).gt.dg%iota2taylor(k,ll,1).and.
     &                 abs((dg%iota2taylorvert(k,ll,lll)-dg%iota2taylor(k,ll,1))).gt.1.0E-15 ) then 

                     dg%alphaiota20(k,ll,lll) = min(1.D0,( iota2maxel(k,ll)
     &                    -dg%iota2taylor(k,ll,1) )/ (dg%iota2taylorvert(k,ll,lll) - dg%iota2taylor(k,ll,1)))
                     
                  elseif (dg%iota2taylorvert(k,ll,lll).eq.dg%iota2taylor(k,ll,1).or.
     &                    abs((dg%iota2taylorvert(k,ll,lll)-dg%iota2taylor(k,ll,1))).le.1.0E-15 ) then

                     dg%alphaiota20(k,ll,lll) = 1.D0

                  elseif (dg%iota2taylorvert(k,ll,lll).lt.dg%iota2taylor(k,ll,1).and.
     &                    abs((dg%iota2taylorvert(k,ll,lll)-dg%iota2taylor(k,ll,1))).gt.1.0E-15) then

                     dg%alphaiota20(k,ll,lll) = min( 1.D0,( iota2minel(k,ll)
     &                    -dg%iota2taylor(k,ll,1) )/( dg%iota2taylorvert(k,ll,lll)-dg%iota2taylor(k,ll,1)))
                     
                  endif
#endif

                  if ( dg%QXtaylorvert(k,ll,lll).gt.dg%QXtaylor(k,ll,1).and.
     &                 (dg%QXtaylorvert(k,ll,lll)-dg%QXtaylor(k,ll,1)).gt.1.0E-15 ) then !in dg%xi1

                     dg%alphaQX0(k,ll,lll) = min( 1.D0, ( QXmaxel(k,ll) 
     &                    - dg%QXtaylor(k,ll,1) ) / ( dg%QXtaylorvert(k,ll,lll) - dg%QXtaylor(k,ll,1) ) )
                     
                  elseif ( dg%QXtaylorvert(k,ll,lll).eq.dg%QXtaylor(k,ll,1).or.
     &                    (dg%QXtaylorvert(k,ll,lll)-dg%QXtaylor(k,ll,1)).le.1.0E-15  ) then

                     dg%alphaQX0(k,ll,lll) = 1.D0

                  elseif ( dg%QXtaylorvert(k,ll,lll).lt.dg%QXtaylor(k,ll,1).and.
     &                    (dg%QXtaylorvert(k,ll,lll)-dg%QXtaylor(k,ll,1)).gt.1.0E-15 ) then

                     dg%alphaQX0(k,ll,lll) = min( 1.D0, ( QXminel(k,ll)
     &                    - dg%QXtaylor(k,ll,1) ) / ( dg%QXtaylorvert(k,ll,lll) - dg%QXtaylor(k,ll,1) ) )

                  endif

                  if ( dg%QYtaylorvert(k,ll,lll).gt.dg%QYtaylor(k,ll,1).and.
     &                 (dg%QYtaylorvert(k,ll,lll)-dg%QYtaylor(k,ll,1)).gt.1.0E-15 ) then !in dg%xi1

                     dg%alphaQY0(k,ll,lll) = min( 1.D0, ( QYmaxel(k,ll) 
     &                    - dg%QYtaylor(k,ll,1) ) / ( dg%QYtaylorvert(k,ll,lll) - dg%QYtaylor(k,ll,1) ) )
                     
                  elseif ( dg%QYtaylorvert(k,ll,lll).eq.dg%QYtaylor(k,ll,1).or.
     &                    (dg%QYtaylorvert(k,ll,lll)-dg%QYtaylor(k,ll,1)).le.1.0E-15  ) then

                     dg%alphaQY0(k,ll,lll) = 1.D0

                  elseif ( dg%QYtaylorvert(k,ll,lll).lt.dg%QYtaylor(k,ll,1).and.
     &                    (dg%QYtaylorvert(k,ll,lll)-dg%QYtaylor(k,ll,1)).gt.1.0E-15 ) then

                     dg%alphaQY0(k,ll,lll) = min( 1.D0, ( QYminel(k,ll) 
     &                    - dg%QYtaylor(k,ll,1) ) / ( dg%QYtaylorvert(k,ll,lll)  - dg%QYtaylor(k,ll,1) ) )

                  endif

               enddo            !lll

            enddo               !ll

         elseif (dg%dofs(k).eq.1) then

            cycle

         endif

      enddo

!.....Find the prescribed higher limiters by finding smallest local value

      dg%alphaZE = 0.D0
      dg%alphaQX = 0.D0
      dg%alphaQY = 0.D0

#ifdef TRACE
      dg%alphaiota = 0.D0
#endif

#ifdef CHEM
      dg%alphaiota = 0.D0
      dg%alphaiota2 = 0.D0
#endif

      do k=1,ne

         if (dg%dofs(k).gt.1) then

            do ll=1,dg%dofs(k)		

               dg%alphaZE(k,ll) = minval( dg%alphaZE0(k,ll,:) )
               dg%alphaQX(k,ll) = minval( dg%alphaQX0(k,ll,:) )
               dg%alphaQY(k,ll) = minval( dg%alphaQY0(k,ll,:) )

#ifdef TRACE
               dg%alphaiota(k,ll) = minval( dg%alphaiota0(k,ll,:) )
#endif

#ifdef CHEM
               dg%alphaiota(k,ll) = minval( dg%alphaiota0(k,ll,:) )
               dg%alphaiota2(k,ll) = minval( dg%alphaiota20(k,ll,:) )
#endif

            enddo

         elseif (dg%dofs(k).eq.1) then

            cycle

         endif

      enddo

!.... Choose smallest (minimum) alpha for derivative in x or y

      dg%alphaZEm = 0.D0
      dg%alphaQXm = 0.D0
      dg%alphaQYm = 0.D0

#ifdef TRACE
      dg%alphaiotam = 0.D0
#endif

#ifdef CHEM
      dg%alphaiotam = 0.D0
      dg%alphaiota2m = 0.D0
#endif

      do k=1,ne

         if (dg%dofs(k).gt.1) then
            
            do bb = 1,pdg_el(k)

               if( (bb+1)*(bb+2)/2.le.dg%dofs(k) ) then

                  dg%alphaZEm(k,(bb*(bb+1))/2 + 1) = 
     &                 minval( dg%alphaZE(k,(bb*(bb-1))/2 +1:(bb*(bb+1)) / 2  ) )
                  dg%alphaQXm(k,(bb*(bb+1))/2 + 1) = 
     &                 minval( dg%alphaQX(k,(bb*(bb-1))/2 +1:(bb*(bb+1)) / 2  ) )
                  dg%alphaQYm(k,(bb*(bb+1))/2 + 1) = 
     &                 minval( dg%alphaQY(k,(bb*(bb-1))/2 +1:(bb*(bb+1)) / 2  ) )

#ifdef TRACE
                  dg%alphaiotam(k,(bb*(bb+1))/2 + 1) = 
     &                 minval( dg%alphaiota(k,(bb*(bb-1))/2 +1:(bb*(bb+1)) / 2 ) )
#endif

#ifdef CHEM
                  dg%alphaiotam(k,(bb*(bb+1))/2 + 1) = 
     &                 minval( dg%alphaiota(k,(bb*(bb-1))/2 +1:(bb*(bb+1)) / 2  ) )
                  dg%alphaiota2m(k,(bb*(bb+1))/2 + 1) = 
     &                 minval( dg%alphaiota2(k,(bb*(bb-1))/2 +1:(bb*(bb+1)) / 2 ) )
#endif

               endif
               
            enddo

         elseif (dg%dofs(k).eq.1) then

            cycle

         endif

      enddo

!.....Use max higher derivative values for final limiter value

      dg%alphaZE_max = 0.D0
      dg%alphaQX_max = 0.D0
      dg%alphaQY_max = 0.D0

#ifdef TRACE
      dg%alphaiota_max = 0.D0
#endif

#ifdef CHEM
      dg%alphaiota_max = 0.D0
      dg%alphaiota2_max = 0.D0
#endif

      do k=1,ne

         if (dg%dofs(k).gt.1) then

            do bb =1,pdg_el(k)

               if( (bb+1)*(bb+2)/2.le.dg%dofs(k)) then

                  dg%alphaZE_max(k,(bb*(bb+1))/2 + 1) = fd*maxval( 
     &                 dg%alphaZEm(k,(bb*(bb+1))/2 + 1:dg%dofs(k)) )
                  dg%alphaQX_max(k,(bb*(bb+1))/2 + 1) = fd*maxval( 
     &                 dg%alphaQXm(k,(bb*(bb+1))/2 + 1:dg%dofs(k)) )
                  dg%alphaQY_max(k,(bb*(bb+1))/2 + 1) = fd*maxval( 
     &                 dg%alphaQYm(k,(bb*(bb+1))/2 + 1:dg%dofs(k)) )

#ifdef TRACE
                  dg%alphaiota_max(k,(bb*(bb+1))/2 + 1) = fd*maxval( 
     &                 dg%alphaiotam(k,(bb*(bb+1))/2 + 1:dg%dofs(k)) )
#endif

#ifdef CHEM
                  dg%alphaiota_max(k,(bb*(bb+1))/2 + 1) = fd*maxval( 
     &                 dg%alphaiotam(k,(bb*(bb+1))/2 + 1:dg%dofs(k)) )
                  dg%alphaiota2_max(k,(bb*(bb+1))/2 + 1) = fd*maxval( 
     &                 dg%alphaiota2m(k,(bb*(bb+1))/2 + 1:dg%dofs(k)) )
#endif

               endif

            enddo

         elseif (dg%dofs(k).eq.1) then

            cycle

         endif

      enddo

!.....Limit on the Master element in the Taylor basis, via reconstruction 
!.....of unconstrained solutions with alpha constraints


      dg%limitZE = 0.D0
      dg%limitQX = 0.D0
      dg%limitQY = 0.D0

      dg%lim_count_roll = 0

#ifdef TRACE
      dg%limitiota = 0.D0
#endif

#ifdef CHEM
      dg%limitiota = 0.D0
      dg%limitiota2 = 0.D0
#endif

      do k=1,ne

         dg%lim_count = 0

         if (dg%dofs(k).gt.1) then

            do ll=1,dg%dofs(k)

               if ( ll.eq.1 ) then

                  dg%limitZE(k,ll) = dg%ZEtaylor(k,ll,1)
                  dg%limitQX(k,ll) = dg%QXtaylor(k,ll,1)
                  dg%limitQY(k,ll) = dg%QYtaylor(k,ll,1) 

#ifdef TRACE
                  dg%limitiota(k,ll) = dg%iotataylor(k,ll,1)
#endif

#ifdef CHEM
                  dg%limitiota(k,ll) = dg%iotataylor(k,ll,1)
                  dg%limitiota2(k,ll) = dg%iota2taylor(k,ll,1)
#endif

               elseif ( ll.ge.2 ) then
                  
                  do bb=1,pdg_el(k)

                     if ( ll.le.( ( (bb+1)*(bb+2)) / 2.D0 ).and.(ll.gt.
     &                    (bb*(bb+1)/2.D0) ) ) then

                        dg%limitZE(k,ll) = dg%alphaZE_max(k,(bb*(bb+1))/2 + 1) 
     &                       * dg%ZEtaylor(k,ll,1)
                        dg%limitQX(k,ll) = dg%alphaQX_max(k,(bb*(bb+1))/2 + 1) 
     &                       * dg%ZEtaylor(k,ll,1)
                        dg%limitQY(k,ll) = dg%alphaQY_max(k,(bb*(bb+1))/2 + 1) 
     &                       * dg%ZEtaylor(k,ll,1)

#ifdef TRACE
                        dg%limitiota(k,ll) = dg%alphaiota_max(k,(bb*(bb+1))/2 + 1) 
     &                       * dg%iotataylor(k,ll,1)
#endif

#ifdef CHEM
                        dg%limitiota(k,ll) = dg%alphaiota_max(k,(bb*(bb+1))/2 + 1) 
     &                       * dg%iotataylor(k,ll,1)
                        dg%limitiota2(k,ll) = dg%alphaiota2_max(k,(bb*(bb+1))/2 + 1) 
     &                       * dg%iota2taylor(k,ll,1)
#endif


!$$$  ! Make a counter to track limiting
!$$$  
!$$$  if ( ( dg%alphaZE_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
!$$$  &                       dg%alphaQX_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
!$$$  &                       dg%alphaQY_max(k,(bb*(bb+1))/2 + 1).lt.1.D0 ).and.
!$$$  &                       chem_flag.ne.1.and.tracer_flag.ne.1 ) then
!$$$  
!$$$  dg%lim_count = 1  
!$$$  
!$$$  elseif ( (dg%alphaZE_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
!$$$  &                          dg%alphaQX_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
!$$$  &                          dg%alphaQY_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
!$$$  &                          dg%alphaiota_max(k,(bb*(bb+1))/2 + 1).lt.1.D0 ).and.
!$$$  &                          chem_flag.ne.1.and.tracer_flag.eq.1 ) then
!$$$  
!$$$  dg%lim_count = 1   
!$$$  
!$$$  elseif ( (dg%alphaZE_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
!$$$  &                          dg%alphaQX_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
!$$$  &                          dg%alphaQY_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
!$$$  &                          dg%alphaiota_max(k,(bb*(bb+1))/2 + 1).lt.1.D0.or.
!$$$  &                          dg%alphaiota2_max(k,(bb*(bb+1))/2 + 1).lt.1.D0 ).and.
!$$$  &                          chem_flag.eq.1 ) then
!$$$  
!$$$  dg%lim_count = 1 
!$$$  
!$$$  endif
                        
                     endif
                     
                  enddo
                  
               endif

            enddo

         elseif (dg%dofs(k).eq.1) then

            cycle

         endif

                                !dg%lim_count_roll = dg%lim_count_roll + dg%lim_count

      enddo

!.....Transform back to the Dubiner basis (multiply by dg%NmatrixInv),

      dg%ZEconst =  0.D0
      dg%QXconst =  0.D0
      dg%QYconst =  0.D0

#ifdef TRACE
      dg%iotaconst =  0.D0
#endif

#ifdef CHEM
      dg%iotaconst =  0.D0
      dg%iota2const =  0.D0
#endif

      do k=1,ne

         if (dg%dofs(k).gt.1) then
                                !do lll=1,3

            do ll=1,dg%dofs(k)

               do ss=1,dg%dofs(k)

                  dg%ZEconst(k,ll) = dg%ZEconst(k,ll) + dg%NmatrixInv(k,ll,ss,dg%dofs(k)) 
     &                 * dg%limitZE(k,ss)
                  dg%QXconst(k,ll) = dg%QXconst(k,ll) + dg%NmatrixInv(k,ll,ss,dg%dofs(k)) 
     &                 * dg%limitQX(k,ss)
                  dg%QYconst(k,ll) = dg%QYconst(k,ll) + dg%NmatrixInv(k,ll,ss,dg%dofs(k)) 
     &                 * dg%limitQY(k,ss)

#ifdef TRACE
                  dg%iotaconst(k,ll) = dg%iotaconst(k,ll) + 
     &                 dg%NmatrixInv(k,ll,ss,dg%dofs(k)) * dg%limitiota(k,ss)
#endif

#ifdef CHEM
                  dg%iotaconst(k,ll) = dg%iotaconst(k,ll) + 
     &                 dg%NmatrixInv(k,ll,ss,dg%dofs(k)) * dg%limitiota(k,ss)
                  dg%iota2const(k,ll) = dg%iota2const(k,ll) + 
     &                 dg%NmatrixInv(k,ll,ss,dg%dofs(k)) * dg%limitiota2(k,ss)
#endif

               enddo

            enddo

         elseif (dg%dofs(k).eq.1) then

            cycle

         endif

      enddo

!.....Set limit values

      do k=1,ne

         if (dg%dofs(k).gt.1) then

            do ll = 1,dg%dofs(k)

               dg%ZE(ll,k,dg%irk+1) = dg%ZEconst(k,ll)
               dg%QX(ll,k,dg%irk+1) = dg%QXconst(k,ll)
               dg%QY(ll,k,dg%irk+1) = dg%QYconst(k,ll)

#ifdef TRACE
               dg%iota(ll,k,dg%irk+1) = dg%iotaconst(k,ll)
#endif

#ifdef CHEM
               dg%iota(ll,k,dg%irk+1) = dg%iotaconst(k,ll)
               dg%iota2(ll,k,dg%irk+1) = dg%iota2const(k,ll)
#endif

            enddo

         elseif (dg%dofs(k).eq.1) then

            cycle

         endif

      enddo
      
      return
      end subroutine


!***********************************************************************
!     
!     SUBROUTINE SLOPELIMITER7()
!     
!     Written 2011
!     
!     This subroutine selects the hierarchical reconstruction method over
!     a Taylor Polynomial basis, and is consistent with p-adaptation
!     to arbitrary order p
!     
!     --cem
!     
!***********************************************************************

      SUBROUTINE SLOPELIMITER7()

!.....Use appropriate modules

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


      Allocate ( ZE_cand(mne,dg%dofh,1),QX_cand(mne,dg%dofh,1) )
      Allocate ( QY_cand(mne,dg%dofh,1) ) 
      Allocate ( iota_cand(mne,dg%dofh,1),iota2_cand(mne,dg%dofh,1) )
      Allocate ( ZE_minmod(mne,dg%dofh,1),QX_minmod(mne,dg%dofh,1) )
      Allocate ( QY_minmod(mne,dg%dofh,1) ) 
      Allocate ( iota_minmod(mne,dg%dofh,1),iota2_minmod(mne,dg%dofh,1) )
      Allocate ( ZEtaylor0slope(mne,dg%dofh,1),ZEtaylor0(MNE,dg%dofh,1))
      Allocate ( QXtaylor0slope(mne,dg%dofh,1),QXtaylor0(MNE,dg%dofh,1))
      Allocate ( QYtaylor0slope(mne,dg%dofh,1),QYtaylor0(MNE,dg%dofh,1))
      Allocate ( iotataylor0slope(mne,dg%dofh,1),iotataylor0(MNE,dg%dofh,1))
      Allocate ( iota2taylor0slope(mne,dg%dofh,1),iota2taylor0(MNE,dg%dofh,1))





!.....We work over the master element
!.....Set initial values

      MUSCL = 0                 ! add weight for lower order pieces (fd<1 => stronger limiting)
      ENO = 1
      Mixed_MUSCL = 0
      Mixed_ENO = 0
      neigh_flag = 1
      
      DO k=1,NE

         if (dg%dofs(k).gt.1) then

            DO ll = 1,dg%dofs(k)

               dg%ZE(ll,k,dg%NRK+2) = dg%ZE(ll,k,dg%IRK+1)
               dg%QX(ll,k,dg%NRK+2) = dg%QX(ll,k,dg%IRK+1)
               dg%QY(ll,k,dg%NRK+2) = dg%QY(ll,k,dg%IRK+1)

#ifdef TRACE
               dg%iota(ll,k,dg%NRK+2) = dg%iota(ll,k,dg%IRK+1)
#endif

#ifdef CHEM
               dg%iota(ll,k,dg%NRK+2) = dg%iota(ll,k,dg%IRK+1)
               dg%iota2(ll,k,dg%NRK+2) = dg%iota2(ll,k,dg%IRK+1)
#endif

            ENDDO

         elseif (dg%dofs(k).eq.1) then

            cycle

         endif

      ENDDO   

      marea = 2.D0              !master elements area

!.....Convert initial values to the Taylor basis (multiply by dg%Nmatrix) on base element


      dg%ZEtaylor = 0.D0 
      dg%QXtaylor = 0.D0
      dg%QYtaylor = 0.D0

#ifdef TRACE
      dg%iotataylor = 0.D0      
#endif

#ifdef CHEM
      dg%iotataylor = 0.D0
      dg%iota2taylor = 0.D0
#endif

      do k=1,ne

         if (dg%dofs(k).gt.1) then

            do ll=1,dg%dofs(k)
               
               do ss=1,dg%dofs(k)

                  dg%ZEtaylor(k,ll,1) = dg%ZEtaylor(k,ll,1) + dg%Nmatrix(k,ll,ss,dg%dofs(k)) * dg%ZE(ss,k,dg%nrk+2)
                  dg%QXtaylor(k,ll,1) = dg%QXtaylor(k,ll,1) + dg%Nmatrix(k,ll,ss,dg%dofs(k)) * dg%QX(ss,k,dg%nrk+2)
                  dg%QYtaylor(k,ll,1) = dg%QYtaylor(k,ll,1) + dg%Nmatrix(k,ll,ss,dg%dofs(k)) * dg%QY(ss,k,dg%nrk+2)

#ifdef TRACE
                  dg%iotataylor(k,ll,1) = dg%iotataylor(k,ll,1) + dg%Nmatrix(k,ll,ss,dg%dofs(k))*dg%iota(ss,k,dg%nrk+2)
#endif

#ifdef CHEM
                  dg%iotataylor(k,ll,1) = dg%iotataylor(k,ll,1) + dg%Nmatrix(k,ll,ss,dg%dofs(k))*dg%iota(ss,k,dg%nrk+2)
                  dg%iota2taylor(k,ll,1) = dg%iota2taylor(k,ll,1) + dg%Nmatrix(k,ll,ss,dg%dofs(k))*dg%iota2(ss,k,dg%nrk+2)
#endif
                  
               enddo

            enddo

         elseif (dg%dofs(k).eq.1) then

            cycle

         endif

      enddo


!.....Generate linear recostructions at vertices

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

      CALL UPDATER_ELEM_MOD2(dg%ZEtaylor,dg%QXtaylor,dg%QYtaylor,1,3)


#ifdef TRACE
      CALL UPDATER_ELEM_MOD2(dg%iotataylor,dg%QXtaylor,dg%QYtaylor,1,2 )  
#endif

#ifdef CHEM
      CALL UPDATER_ELEM_MOD2(dg%iotataylor,dg%iota2taylor,dg%QYtaylor,1,2 )  
#endif

#endif


      do k=1,ne

         do ll = dg%dofs(k),3,-1

            ZEtaylor0(k,ll,1) = dg%ZEtaylor(k,ll,1)
            ZEtaylor0slope(k,ll,1) = 0.D0
            QXtaylor0(k,ll,1) = dg%QXtaylor(k,ll,1)
            QXtaylor0slope(k,ll,1) = 0.D0
            QYtaylor0(k,ll,1) = dg%QYtaylor(k,ll,1)
            QYtaylor0slope(k,ll,1) = 0.D0
            
#ifdef TRACE
            iotataylor0(k,ll,1) = dg%iotataylor(k,ll,1)
            iotataylor0slope(k,ll,1) = 0.D0
#endif

#ifdef CHEM            
            iotataylor0(k,ll,1) = dg%iotataylor(k,ll,1)
            iotataylor0slope(k,ll,1) = 0.D0
            iota2taylor0(k,ll,1) = dg%iota2taylor(k,ll,1)
            iota2taylor0slope(k,ll,1) = 0.D0
#endif

            if (ll.gt.3) then

                                ! Construct the coeffs in nonlinear reconstruction

               do i = 0,pdg_el(k)

                  do j = 0,pdg_el(k)

                     b_index = (dg%bj(ll)+j) * (dg%bj(ll) + 1 + j) / 2 + ( dg%bj(ll)+j ) * ( dg%bi(ll)+i ) 
     &                    + ( dg%bi(ll) + i )*(dg%bi(ll) + 3 + i)/2 + 1

                     if ( b_index.le.dg%dofs(k).and.(i+j).gt.0) then


                        Call factorial(i,dg%fact(i))
                        Call factorial(j,dg%fact(j))
                        
                        do mm = 1,dg%nagp(pdg_el(k))

                           ZEtaylor0(k,ll,1) =   ZEtaylor0(k,ll,1) + dg%ZEtaylor(k,b_index,1) * 
     &                          (1.D0/dg%fact(i)*dg%fact(j))* ( ( dg%yagp(mm,pdg_el(k)) - dg%xi2BCb(k) )**i
     &                          *( dg%xagp(mm,pdg_el(k)) - dg%xi1BCb(k) )**j  )  * dg%wagp(mm,pdg_el(k))

                           ZEtaylor0slope(k,ll,1) =  ZEtaylor0slope(k,ll,1) + dg%ZEtaylor(k,b_index,1) * 
     &                          (1.D0/dg%fact(i)*dg%fact(j))* ( ( dg%yagp(mm,pdg_el(k)) - dg%xi2BCb(k) )**i
     &                          *( dg%xagp(mm,pdg_el(k)) - dg%xi1BCb(k) )**j  )  * dg%wagp(mm,pdg_el(k))

                           QXtaylor0(k,ll,1) =   QXtaylor0(k,ll,1) + dg%QXtaylor(k,b_index,1) * 
     &                          (1.D0/dg%fact(i)*dg%fact(j))* ( ( dg%yagp(mm,pdg_el(k)) - dg%xi2BCb(k) )**i
     &                          *( dg%xagp(mm,pdg_el(k)) - dg%xi1BCb(k) )**j  )  * dg%wagp(mm,pdg_el(k))

                           QXtaylor0slope(k,ll,1) =  QXtaylor0slope(k,ll,1) + dg%QXtaylor(k,b_index,1) * 
     &                          (1.D0/dg%fact(i)*dg%fact(j))* ( ( dg%yagp(mm,pdg_el(k)) - dg%xi2BCb(k) )**i
     &                          *( dg%xagp(mm,pdg_el(k)) - dg%xi1BCb(k) )**j  )  * dg%wagp(mm,pdg_el(k))

                           QYtaylor0(k,ll,1) =   QYtaylor0(k,ll,1) + dg%QYtaylor(k,b_index,1) * 
     &                          (1.D0/dg%fact(i)*dg%fact(j))* ( ( dg%yagp(mm,pdg_el(k)) - dg%xi2BCb(k) )**i
     &                          *( dg%xagp(mm,pdg_el(k)) - dg%xi1BCb(k) )**j  )  * dg%wagp(mm,pdg_el(k))

                           QYtaylor0slope(k,ll,1) =  QYtaylor0slope(k,ll,1) + dg%QYtaylor(k,b_index,1) * 
     &                          (1.D0/dg%fact(i)*dg%fact(j))* ( ( dg%yagp(mm,pdg_el(k)) - dg%xi2BCb(k) )**i
     &                          *( dg%xagp(mm,pdg_el(k)) - dg%xi1BCb(k) )**j  )  * dg%wagp(mm,pdg_el(k))

#ifdef TRACE
                           iotataylor0(k,ll,1) =   iotataylor0(k,ll,1) + dg%iotataylor(k,b_index,1) * 
     &                          (1.D0/dg%fact(i)*dg%fact(j))* ( ( dg%yagp(mm,pdg_el(k)) - dg%xi2BCb(k) )**i
     &                          *( dg%xagp(mm,pdg_el(k)) - dg%xi1BCb(k) )**j  )  * dg%wagp(mm,pdg_el(k))

                           iotataylor0slope(k,ll,1) =  iotataylor0slope(k,ll,1) + dg%iotataylor(k,b_index,1) * 
     &                          (1.D0/dg%fact(i)*dg%fact(j))* ( ( dg%yagp(mm,pdg_el(k)) - dg%xi2BCb(k) )**i
     &                          *( dg%xagp(mm,pdg_el(k)) - dg%xi1BCb(k) )**j  )  * dg%wagp(mm,pdg_el(k))
#endif

#ifdef CHEM
                           iotataylor0(k,ll,1) =   iotataylor0(k,ll,1) + dg%iotataylor(k,b_index,1) * 
     &                          (1.D0/dg%fact(i)*dg%fact(j))* ( ( dg%yagp(mm,pdg_el(k)) - dg%xi2BCb(k) )**i
     &                          *( dg%xagp(mm,pdg_el(k)) - dg%xi1BCb(k) )**j  )  * dg%wagp(mm,pdg_el(k))

                           iotataylor0slope(k,ll,1) =  iotataylor0slope(k,ll,1) + dg%iotataylor(k,b_index,1) * 
     &                          (1.D0/dg%fact(i)*dg%fact(j))* ( ( dg%yagp(mm,pdg_el(k)) - dg%xi2BCb(k) )**i
     &                          *( dg%xagp(mm,pdg_el(k)) - dg%xi1BCb(k) )**j  )  * dg%wagp(mm,pdg_el(k))

                           iota2taylor0(k,ll,1) =   iota2taylor0(k,ll,1) + dg%iota2taylor(k,b_index,1) * 
     &                          (1.D0/dg%fact(i)*dg%fact(j))* ( ( dg%yagp(mm,pdg_el(k)) - dg%xi2BCb(k) )**i
     &                          *( dg%xagp(mm,pdg_el(k)) - dg%xi1BCb(k) )**j  )  * dg%wagp(mm,pdg_el(k))

                           iota2taylor0slope(k,ll,1) =  iota2taylor0slope(k,ll,1) + dg%iota2taylor(k,b_index,1) * 
     &                          (1.D0/dg%fact(i)*dg%fact(j))* ( ( dg%yagp(mm,pdg_el(k)) - dg%xi2BCb(k) )**i
     &                          *( dg%xagp(mm,pdg_el(k)) - dg%xi1BCb(k) )**j  )  * dg%wagp(mm,pdg_el(k))
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

         if (dg%dofs(k).gt.3) then
            
            do ll=dg%dofs(k),3,-1

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

!$$$  #ifdef CMPI
!$$$  
!$$$  CALL UPDATER_ELEM_MOD2(ZE_minmod,QX_minmod,QY_minmod,1,3)
!$$$  
!$$$  #ifdef TRACE
!$$$  CALL UPDATER_ELEM_MOD2(iota_minmod,iota_minmod,QY_minmod,1,2 )
!$$$  #endif
!$$$  
!$$$  #ifdef CHEM
!$$$  CALL UPDATER_ELEM_MOD2(iota_minmod,iota2_minmod,QY_minmod,1,2 )
!$$$  #endif
!$$$  
!$$$  #endif

                                !Apply the friendly minmod regime
      
      do k=1,mne

         do ll=1,dg%dofs(k)

            ZE_minmod(k,ll,1) = dg%ZEtaylor(k,ll,1)
            QX_minmod(k,ll,1) = dg%QXtaylor(k,ll,1)
            QY_minmod(k,ll,1) = dg%QYtaylor(k,ll,1)

#ifdef TRACE
            iota_minmod(k,ll,1) = dg%iotataylor(k,ll,1)
#endif

#ifdef CHEM
            iota_minmod(k,ll,1) = dg%iotataylor(k,ll,1)
            iota2_minmod(k,ll,1) = dg%iota2taylor(k,ll,1)
#endif
            if (ll.gt.1) then

               if (neigh_flag.eq.1) then !focal neighbors

                  if (MUSCL.eq.1) then !construct MUSCL minmod


                     if (minval(ZE_cand(dg%focal_neigh(k,1:dg%focal_up(k)),ll,1)).gt.0.D0 ) then

                        ZE_minmod(k,ll,1) = minval(ZE_cand(dg%focal_neigh(k,1:dg%focal_up(k)),ll,1))

                     elseif (maxval(ZE_cand(dg%focal_neigh(k,1:dg%focal_up(k)),ll,1)).lt.0.D0 ) then

                        ZE_minmod(k,ll,1) = maxval(ZE_cand(dg%focal_neigh(k,1:dg%focal_up(k)),ll,1))

                     else

                        ZE_minmod(k,ll,1) = 0.D0

                     endif

                     if (minval(QX_cand(dg%focal_neigh(k,1:dg%focal_up(k)),ll,1)).gt.0.D0 ) then

                        QX_minmod(k,ll,1) = minval(QX_cand(dg%focal_neigh(k,1:dg%focal_up(k)),ll,1))

                     elseif (maxval(QX_cand(dg%focal_neigh(k,1:dg%focal_up(k)),ll,1)).lt.0.D0 ) then

                        QX_minmod(k,ll,1) = maxval(QX_cand(dg%focal_neigh(k,1:dg%focal_up(k)),ll,1))

                     else

                        QX_minmod(k,ll,1) = 0.D0

                     endif

                     if (minval(QY_cand(dg%focal_neigh(k,1:dg%focal_up(k)),ll,1)).gt.0.D0 ) then

                        QY_minmod(k,ll,1) = minval(QY_cand(dg%focal_neigh(k,1:dg%focal_up(k)),ll,1))

                     elseif (maxval(QY_cand(dg%focal_neigh(k,1:dg%focal_up(k)),ll,1)).lt.0.D0 ) then

                        QY_minmod(k,ll,1) = maxval(QY_cand(dg%focal_neigh(k,1:dg%focal_up(k)),ll,1))

                     else

                        QY_minmod(k,ll,1) = 0.D0

                     endif
                     
                     
                  elseif (ENO.eq.1) then !construct ENO minmod

                     ZE_minmod(k,ll,1) = minval(abs(ZE_cand(dg%focal_neigh(k,1:dg%focal_up(k)),ll,1)))
                     QX_minmod(k,ll,1) = minval(abs(ZE_cand(dg%focal_neigh(k,1:dg%focal_up(k)),ll,1)))
                     QY_minmod(k,ll,1) = minval(abs(ZE_cand(dg%focal_neigh(k,1:dg%focal_up(k)),ll,1)))

                  elseif (Mixed_MUSCL.eq.1) then !construct Mixed minmod

                  elseif (Mixed_ENO.eq.1) then !construct Mixed minmod

                  endif

#ifdef TRACE
                  if (MUSCL.eq.1) then !construct MUSCL minmod


                     if (minval(iota_cand(dg%focal_neigh(k,1:dg%focal_up(k)),ll,1)).gt.0.D0 ) then

                        iota_minmod(k,ll,1) = minval(iota_cand(dg%focal_neigh(k,1:dg%focal_up(k)),ll,1))

                     elseif (maxval(iota_cand(dg%focal_neigh(k,1:dg%focal_up(k)),ll,1)).lt.0.D0 ) then

                        iota_minmod(k,ll,1) = maxval(iota_cand(dg%focal_neigh(k,1:dg%focal_up(k)),ll,1))

                     else

                        iota_minmod(k,ll,1) = 0.D0

                     endif
                     
                     
                  elseif (ENO.eq.1) then !construct ENO minmod

                     iota_minmod(k,ll,1) = minval(abs(iota_cand(dg%focal_neigh(k,1:dg%focal_up(k)),ll,1)))

                  elseif (Mixed_MUSCL.eq.1) then !construct Mixed minmod

                  elseif (Mixed_ENO.eq.1) then !construct Mixed minmod

                  endif
#endif

#ifdef CHEM
                  if (MUSCL.eq.1) then !construct MUSCL minmod


                     if (minval(iota_cand(dg%focal_neigh(k,1:dg%focal_up(k)),ll,1)).gt.0.D0 ) then

                        iota_minmod(k,ll,1) = minval(iota_cand(dg%focal_neigh(k,1:dg%focal_up(k)),ll,1))

                     elseif (maxval(iota_cand(dg%focal_neigh(k,1:dg%focal_up(k)),ll,1)).lt.0.D0 ) then

                        iota_minmod(k,ll,1) = maxval(iota_cand(dg%focal_neigh(k,1:dg%focal_up(k)),ll,1))

                     else

                        iota_minmod(k,ll,1) = 0.D0

                     endif
                     
                     
                  elseif (ENO.eq.1) then !construct ENO minmod

                     iota_minmod(k,ll,1) = minval(abs(iota_cand(dg%focal_neigh(k,1:dg%focal_up(k)),ll,1)))

                  elseif (Mixed_MUSCL.eq.1) then !construct Mixed minmod

                  elseif (Mixed_ENO.eq.1) then !construct Mixed minmod

                  endif

                  if (MUSCL.eq.1) then !construct MUSCL minmod


                     if (minval(iota2_cand(dg%focal_neigh(k,1:dg%focal_up(k)),ll,1)).gt.0.D0 ) then

                        iota2_minmod(k,ll,1) = minval(iota2_cand(dg%focal_neigh(k,1:dg%focal_up(k)),ll,1))

                     elseif (maxval(iota2_cand(dg%focal_neigh(k,1:dg%focal_up(k)),ll,1)).lt.0.D0 ) then

                        iota2_minmod(k,ll,1) = maxval(iota2_cand(dg%focal_neigh(k,1:dg%focal_up(k)),ll,1))

                     else

                        iota2_minmod(k,ll,1) = 0.D0

                     endif
                     
                     
                  elseif (ENO.eq.1) then !construct ENO minmod

                     iota2_minmod(k,ll,1) = minval(abs(iota2_cand(dg%focal_neigh(k,1:dg%focal_up(k)),ll,1)))

                  elseif (Mixed_MUSCL.eq.1) then !construct Mixed minmod

                  elseif (Mixed_ENO.eq.1) then !construct Mixed minmod

                  endif
#endif

               elseif (neigh_flag.eq.2) then !edge neighbors

                  if (MUSCL.eq.1) then !construct MUSCL minmod

                     if (ZE_cand(k,ll,1).gt.0.D0 ) then

                        ZE_minmod(k,ll,1) = minval(ZE_cand(dg%El_nbors(:,k),ll,1))

                     elseif (ZE_cand(k,ll,1).lt.0.D0 ) then

                        ZE_minmod(k,ll,1) = maxval(ZE_cand(dg%El_nbors(:,k),ll,1))

                     else

                        ZE_minmod(k,ll,1) = 0.D0

                     endif

                     if (QX_cand(k,ll,1).gt.0.D0 ) then

                        QX_minmod(k,ll,1) = minval(QX_cand(dg%El_nbors(:,k),ll,1))

                     elseif (QX_cand(k,ll,1).lt.0.D0 ) then

                        QX_minmod(k,ll,1) = maxval(QX_cand(dg%El_nbors(:,k),ll,1))

                     else

                        QX_minmod(k,ll,1) = 0.D0

                     endif

                     if (QY_cand(k,ll,1).gt.0.D0 ) then

                        QY_minmod(k,ll,1) = minval(QY_cand(dg%El_nbors(:,k),ll,1))

                     elseif (QY_cand(k,ll,1).lt.0.D0 ) then

                        QY_minmod(k,ll,1) = maxval(QY_cand(dg%El_nbors(:,k),ll,1))

                     else

                        QY_minmod(k,ll,1) = 0.D0

                     endif

                  elseif (ENO.eq.1) then !construct ENO minmod

                     ZE_minmod(k,ll,1) = minval(abs(ZE_cand(dg%El_nbors(:,k),ll,1)))
                     QX_minmod(k,ll,1) = minval(abs(QX_cand(dg%El_nbors(:,k),ll,1)))
                     QY_minmod(k,ll,1) = minval(abs(QY_cand(dg%El_nbors(:,k),ll,1)))

                  elseif (Mixed_MUSCL.eq.1) then !construct Mixed minmod

                  elseif (Mixed_ENO.eq.1) then !construct Mixed minmod

                  endif

#ifdef TRACE
                  if (MUSCL.eq.1) then !construct MUSCL minmod

                     if (iota_cand(k,ll,1).gt.0.D0 ) then

                        iota_minmod(k,ll,1) = minval(iota_cand(dg%El_nbors(:,k),ll,1))

                     elseif (iota_cand(k,ll,1).lt.0.D0 ) then

                        iota_minmod(k,ll,1) = maxval(iota_cand(dg%El_nbors(:,k),ll,1))

                     else

                        iota_minmod(k,ll,1) = 0.D0

                     endif

                  elseif (ENO.eq.1) then !construct ENO minmod

                     iota_minmod(k,ll,1) = minval(abs(iota_cand(dg%El_nbors(:,k),ll,1)))

                  elseif (Mixed_MUSCL.eq.1) then !construct Mixed minmod

                  elseif (Mixed_ENO.eq.1) then !construct Mixed minmod

                  endif
#endif

#ifdef CHEM
                  if (MUSCL.eq.1) then !construct MUSCL minmod

                     if (iota_cand(k,ll,1).gt.0.D0 ) then

                        iota_minmod(k,ll,1) = minval(iota_cand(dg%El_nbors(:,k),ll,1))

                     elseif (iota_cand(k,ll,1).lt.0.D0 ) then

                        iota_minmod(k,ll,1) = maxval(iota_cand(dg%El_nbors(:,k),ll,1))

                     else

                        iota_minmod(k,ll,1) = 0.D0

                     endif

                  elseif (ENO.eq.1) then !construct ENO minmod

                     iota_minmod(k,ll,1) = minval(abs(iota_cand(dg%El_nbors(:,k),ll,1)))

                  elseif (Mixed_MUSCL.eq.1) then !construct Mixed minmod

                  elseif (Mixed_ENO.eq.1) then !construct Mixed minmod

                  endif

                  if (MUSCL.eq.1) then !construct MUSCL minmod

                     if (iota2_cand(k,ll,1).gt.0.D0 ) then

                        iota2_minmod(k,ll,1) = minval(iota2_cand(dg%El_nbors(:,k),ll,1))

                     elseif (iota2_cand(k,ll,1).lt.0.D0 ) then

                        iota2_minmod(k,ll,1) = maxval(iota2_cand(dg%El_nbors(:,k),ll,1))

                     else

                        iota2_minmod(k,ll,1) = 0.D0

                     endif

                  elseif (ENO.eq.1) then !construct ENO minmod

                     iota2_minmod(k,ll,1) = minval(abs(iota2_cand(dg%El_nbors(:,k),ll,1)))

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

!.....Transform back to the Dubiner basis (multiply by dg%NmatrixInv),

      dg%ZEconst =  0.D0
      dg%QXconst =  0.D0
      dg%QYconst =  0.D0

#ifdef TRACE
      dg%iotaconst =  0.D0
#endif

#ifdef CHEM
      dg%iotaconst =  0.D0
      dg%iota2const =  0.D0
#endif

      do k=1,ne

         if (dg%dofs(k).gt.1) then
                                !do lll=1,3

            do ll=1,dg%dofs(k)

               do ss=1,dg%dofs(k)

                  dg%ZEconst(k,ll) = dg%ZEconst(k,ll) + dg%NmatrixInv(k,ll,ss,dg%dofs(k)) 
     &                 * ZE_minmod(k,ss,1)
                  dg%QXconst(k,ll) = dg%QXconst(k,ll) + dg%NmatrixInv(k,ll,ss,dg%dofs(k)) 
     &                 * QX_minmod(k,ss,1)
                  dg%QYconst(k,ll) = dg%QYconst(k,ll) + dg%NmatrixInv(k,ll,ss,dg%dofs(k)) 
     &                 * QX_minmod(k,ss,1)

#ifdef TRACE
                  dg%iotaconst(k,ll) = dg%iotaconst(k,ll) + 
     &                 dg%NmatrixInv(k,ll,ss,dg%dofs(k)) * iota_minmod(k,ss,1) 
#endif

#ifdef CHEM
                  dg%iotaconst(k,ll) = dg%iotaconst(k,ll) + 
     &                 dg%NmatrixInv(k,ll,ss,dg%dofs(k)) *  iota_minmod(k,ss,1)
                  dg%iota2const(k,ll) = dg%iota2const(k,ll) + 
     &                 dg%NmatrixInv(k,ll,ss,dg%dofs(k)) * iota2_minmod(k,ss,1)
#endif

               enddo

            enddo

         elseif (dg%dofs(k).eq.1) then

            cycle

         endif

      enddo

!.....Set limit values
      
      CALL SLOPELIMITER5(s)

      do k=1,ne

         if (dg%dofs(k).gt.3) then

            do ll = 3,dg%dofs(k)

               dg%ZE(ll,k,dg%irk+1) = dg%ZEconst(k,ll)
               dg%QX(ll,k,dg%irk+1) = dg%QXconst(k,ll)
               dg%QY(ll,k,dg%irk+1) = dg%QYconst(k,ll)

#ifdef TRACE
               dg%iota(ll,k,dg%irk+1) = dg%iotaconst(k,ll)
#endif

#ifdef CHEM
               dg%iota(ll,k,dg%irk+1) = dg%iotaconst(k,ll)
               dg%iota2(ll,k,dg%irk+1) = dg%iota2const(k,ll)
#endif

            enddo

         elseif (dg%dofs(k).eq.1) then

            cycle

         endif

      enddo
      
      return
      end subroutine


!***********************************************************************
!     
!     SUBROUTINE SLOPELIMITER8()
!     
!     Written 2011
!     
!     The hierarchic recombination method, works for adapting p
!     
!     -cem
!     
!     
!***********************************************************************

      SUBROUTINE SLOPELIMITER8()

!.....Use appropriate modules

      use sizes, only : sz
      use global
      use dg

#ifdef CMPI
      USE MESSENGER_ELEM
      USE MESSENGER
#endif

      IMPLICIT NONE

!.....Declare local variables

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
      Allocate ( iota2_MIN1(NP),iota2_MAX1(NP), temparray1(MNE,dg%dofh,1) )
      Allocate (  temparray2(MNE,dg%dofh,1), temparray3(MNE,dg%dofh,1) )
      Allocate (  temparray4(MNE,dg%dofh,1), temparray5(MNE,dg%dofh,1) )

!     FIND THE MAXIMUM AND MINIMUM OF EACH VARIABLE OVER ALL ELEMENTS 
!     SHARING A NODE

      do k=1,ne

         if (dg%dofs(k).gt.1) then

            do ll = 1,dg%dofs(k)

               dg%ze(ll,k,dg%nrk+2) = dg%ze(ll,k,dg%irk+1)
               dg%qx(ll,k,dg%nrk+2) = dg%qx(ll,k,dg%irk+1)
               dg%qy(ll,k,dg%nrk+2) = dg%qy(ll,k,dg%irk+1)

#ifdef TRACE
               dg%iota(ll,k,dg%nrk+2) = dg%iota(ll,k,dg%irk+1)
#endif

#ifdef CHEM
               dg%iota(ll,k,dg%nrk+2) = dg%iota(ll,k,dg%irk+1)
               dg%iota2(ll,k,dg%nrk+2) = dg%iota2(ll,k,dg%irk+1)
#endif

            enddo

         elseif (dg%dofs(k).eq.1) then

            cycle

         endif

      enddo   

!.....convert initial values to the taylor basis (multiply by dg%nmatrix) on base element


      dg%zetaylor = 0.d0 
      dg%qxtaylor = 0.d0
      dg%qytaylor = 0.d0

#ifdef TRACE
      dg%iotataylor = 0.d0
#endif

#ifdef CHEM
      dg%iotataylor = 0.d0
      dg%iota2taylor = 0.d0
#endif

      do k=1,ne

                                !if (dg%dofs(k).gt.1) then

         do ll=1,dg%dofs(k)
            
            do ss=1,dg%dofs(k)

               dg%zetaylor(k,ll,1) = dg%zetaylor(k,ll,1) + dg%nmatrix(k,ll,ss,dg%dofs(k)) * dg%ze(ss,k,dg%irk+1)
               dg%qxtaylor(k,ll,1) = dg%qxtaylor(k,ll,1) + dg%nmatrix(k,ll,ss,dg%dofs(k)) * dg%qx(ss,k,dg%irk+1)
               dg%qytaylor(k,ll,1) = dg%qytaylor(k,ll,1) + dg%nmatrix(k,ll,ss,dg%dofs(k)) * dg%qy(ss,k,dg%irk+1)

#ifdef TRACE
               dg%iotataylor(k,ll,1) = dg%iotataylor(k,ll,1) + 
     &              dg%nmatrix(k,ll,ss,dg%dofs(k))*dg%iota(ss,k,dg%irk+1)
#endif

#ifdef CHEM
               dg%iotataylor(k,ll,1) = dg%iotataylor(k,ll,1) + dg%nmatrix(k,ll,ss,dg%dofs(k))*dg%iota(ss,k,dg%irk+1)
               dg%iota2taylor(k,ll,1) = dg%iota2taylor(k,ll,1) + dg%nmatrix(k,ll,ss,dg%dofs(k))*dg%iota2(ss,k,dg%irk+1)
#endif
               
            enddo

         enddo

      enddo


#ifdef CMPI

      CALL UPDATER_ELEM_MOD2(dg%ZEtaylor,dg%QXtaylor,dg%QYtaylor,1,3)

#ifdef TRACE
      CALL UPDATER_ELEM_MOD2(dg%iotataylor,dg%iotataylor,dg%QYtaylor,1,2 )
#endif

#ifdef CHEM
      CALL UPDATER_ELEM_MOD2(dg%iotataylor,dg%iota2taylor,dg%QYtaylor,1,2 )
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

      do ll=dg%dofh,1,-1

         if(ll+floor(  0.5D0 + sqrt(real(2*ll)) )+1.le.dg%dofh) then


            Call Slopelim_Aux(dg%ZEtaylor(:,:,1),dg%QXtaylor(:,:,1),dg%QYtaylor(:,:,1),   
     &           dg%iotataylor(:,:,1), dg%iota2taylor(:,:,1),ll,
     &           ll+floor(  0.5D0 + sqrt(real(2*ll)) ),
     &           ll+floor(  0.5D0 + sqrt(real(2*ll)) )+1 )

            if (ll-1+floor(  0.5D0 + sqrt(real(2*(ll-1))) +1).eq.ll+floor(  0.5D0 + sqrt(real(2*(ll))) )) then

               temparray1(:,ll+floor(  0.5D0 + sqrt(real(2*ll))),1) = 
     &              dg%ZEtaylor( :,ll+floor( 0.5D0 + sqrt(real(2*ll)) ),1 )
               temparray2(:,ll+floor(  0.5D0 + sqrt(real(2*ll))),1) = 
     &              dg%QXtaylor( :,ll+floor( 0.5D0 + sqrt(real(2*ll)) ),1 )
               temparray3(:,ll+floor(  0.5D0 + sqrt(real(2*ll))),1) = 
     &              dg%QYtaylor( :,ll+floor( 0.5D0 + sqrt(real(2*ll)) ),1 )

#ifdef TRACE
               temparray4(:,ll+floor(  0.5D0 + sqrt(real(2*ll))),1) = 
     &              dg%iotataylor( :,ll+floor( 0.5D0 + sqrt(real(2*ll)) ),1 )
#endif

#ifdef CHEM
               temparray4(:,ll+floor(  0.5D0 + sqrt(real(2*ll))),1) = 
     &              dg%iotataylor( :,ll+floor( 0.5D0 + sqrt(real(2*ll)) ),1 )
               temparray5(:,ll+floor(  0.5D0 + sqrt(real(2*ll))),1) = 
     &              dg%iota2taylor( :,ll+floor( 0.5D0 + sqrt(real(2*ll)) ),1 )
#endif

            endif

            if (ll+floor(  0.5D0 + sqrt(real(2*(ll))) +1).eq.
     &           ll+1+floor(  0.5D0 + sqrt(real(2*(ll+1))) )) then

               do k = 1,ne

                  if (dg%ZEtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll))+1 ),1 ).gt.0.D0.and.
     &                 temparray1(k,ll+1+floor(  0.5D0 + sqrt(real(2*(ll+1)))),1).gt.0.D0 ) then

                     dg%ZEtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll))+1 ),1 ) = 
     &                    min(dg%ZEtaylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) + 1 ),1),
     &                    temparray1(k,ll+1+floor(  0.5D0 + sqrt(real(2*(ll+1)))),1))

                  elseif (dg%ZEtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) +1 ),1 ).lt.0.D0.and.
     &                    temparray1(k,ll+1+floor(  0.5D0 + sqrt(real(2*(ll+1)))),1).lt.0.D0 ) then
                     
                     dg%ZEtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll))+1 ),1 ) =  
     &                    max(dg%ZEtaylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) + 1),1),
     &                    temparray1(k,ll+1+floor(  0.5D0 + sqrt(real(2*(ll+1)))),1))

                  elseif (temparray1(k,ll+1+floor(  0.5D0 + sqrt(real(2*(ll+1)))),1).ne.0.D0) then

                     dg%ZEtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) +1),1 ) = 0.D0

                  endif

                  if (dg%QXtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll))+1 ),1 ).gt.0.D0.and.
     &                 temparray2(k,ll+1+floor(  0.5D0 + sqrt(real(2*(ll+1)))),1).gt.0.D0 ) then

                     dg%QXtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll))+1 ),1 ) = 
     &                    min(dg%QXtaylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) + 1 ),1),
     &                    temparray2(k,ll+1+floor(  0.5D0 + sqrt(real(2*(ll+1)))),1))

                  elseif (dg%QXtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) +1 ),1 ).lt.0.D0.and.
     &                    temparray2(k,ll+1+floor(  0.5D0 + sqrt(real(2*(ll+1)))),1).lt.0.D0 ) then
                     
                     dg%QXtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll))+1 ),1 ) =  
     &                    max(dg%QXtaylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) + 1),1),
     &                    temparray2(k,ll+1+floor(  0.5D0 + sqrt(real(2*(ll+1)))),1))

                  elseif (temparray2(k,ll+1+floor(  0.5D0 + sqrt(real(2*(ll+1)))),1).ne.0.D0) then

                     dg%QXtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) +1),1 ) = 0.D0

                  endif

                  if (dg%QYtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll))+1 ),1 ).gt.0.D0.and.
     &                 temparray3(k,ll+1+floor(  0.5D0 + sqrt(real(2*(ll+1)))),1).gt.0.D0 ) then

                     dg%QYtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll))+1 ),1 ) = 
     &                    min(dg%QYtaylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) + 1 ),1),
     &                    temparray3(k,ll+1+floor(  0.5D0 + sqrt(real(2*(ll+1)))),1))

                  elseif (dg%QYtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) +1 ),1 ).lt.0.D0.and.
     &                    temparray3(k,ll+1+floor(  0.5D0 + sqrt(real(2*(ll+1)))),1).lt.0.D0 ) then
                     
                     dg%QYtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll))+1 ),1 ) =  
     &                    max(dg%QYtaylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) + 1),1),
     &                    temparray3(k,ll+1+floor(  0.5D0 + sqrt(real(2*(ll+1)))),1))

                  elseif (temparray3(k,ll+1+floor(  0.5D0 + sqrt(real(2*(ll+1)))),1).ne.0.D0) then

                     dg%QYtaylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) +1),1 ) = 0.D0

                  endif

#ifdef TRACE
                  if (dg%iotataylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll))+1 ),1 ).gt.0.D0.and.
     &                 temparray4(k,ll+1+floor(  0.5D0 + sqrt(real(2*(ll+1)))),1).gt.0.D0 ) then

                     dg%iotataylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll))+1 ),1 ) = 
     &                    min(dg%iotataylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) + 1 ),1),
     &                    temparray4(k,ll+1+floor(  0.5D0 + sqrt(real(2*(ll+1)))),1))

                  elseif (dg%iotataylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) +1 ),1 ).lt.0.D0.and.
     &                    temparray4(k,ll+1+floor(  0.5D0 + sqrt(real(2*(ll+1)))),1).lt.0.D0 ) then
                     
                     dg%iotataylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll))+1 ),1 ) =  
     &                    max(dg%iotataylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) + 1),1),
     &                    temparray4(k,ll+1+floor(  0.5D0 + sqrt(real(2*(ll+1)))),1))

                  elseif (temparray4(k,ll+1+floor(  0.5D0 + sqrt(real(2*(ll+1)))),1).ne.0.D0) then

                     dg%iotataylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) +1),1 ) = 0.D0

                  endif
#endif

#ifdef CHEM
                  if (dg%iotataylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll))+1 ),1 ).gt.0.D0.and.
     &                 temparray4(k,ll+1+floor(  0.5D0 + sqrt(real(2*(ll+1)))),1).gt.0.D0 ) then

                     dg%iotataylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll))+1 ),1 ) = 
     &                    min(dg%iotataylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) + 1 ),1),
     &                    temparray4(k,ll+1+floor(  0.5D0 + sqrt(real(2*(ll+1)))),1))

                  elseif (dg%iotataylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) +1 ),1 ).lt.0.D0.and.
     &                    temparray4(k,ll+1+floor(  0.5D0 + sqrt(real(2*(ll+1)))),1).lt.0.D0 ) then
                     
                     dg%iotataylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll))+1 ),1 ) =  
     &                    max(dg%iotataylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) + 1),1),
     &                    temparray4(k,ll+1+floor(  0.5D0 + sqrt(real(2*(ll+1)))),1))

                  elseif (temparray4(k,ll+1+floor(  0.5D0 + sqrt(real(2*(ll+1)))),1).ne.0.D0) then

                     dg%iotataylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) +1),1 ) = 0.D0

                  endif

                  if (dg%iota2taylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll))+1 ),1 ).gt.0.D0.and.
     &                 temparray5(k,ll+1+floor(  0.5D0 + sqrt(real(2*(ll+1)))),1).gt.0.D0 ) then

                     dg%iota2taylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll))+1 ),1 ) = 
     &                    min(dg%iota2taylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) + 1 ),1),
     &                    temparray5(k,ll+1+floor(  0.5D0 + sqrt(real(2*(ll+1)))),1))

                  elseif (dg%iota2taylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) +1 ),1 ).lt.0.D0.and.
     &                    temparray5(k,ll+1+floor(  0.5D0 + sqrt(real(2*(ll+1)))),1).lt.0.D0 ) then
                     
                     dg%iota2taylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll))+1 ),1 ) =  
     &                    max(dg%iota2taylor(k,ll+floor(  0.5D0 + sqrt(real(2*ll)) + 1),1),
     &                    temparray5(k,ll+1+floor(  0.5D0 + sqrt(real(2*(ll+1)))),1))

                  elseif (temparray5(k,ll+1+floor(  0.5D0 + sqrt(real(2*(ll+1)))),1).ne.0.D0) then

                     dg%iota2taylor(k,ll+floor( 0.5D0 + sqrt(real(2*ll)) +1),1 ) = 0.D0

                  endif
#endif

               enddo

            endif

         endif
         
      enddo

#ifdef CMPI

      CALL UPDATER_ELEM_MOD2(dg%ZEtaylor,dg%QXtaylor,dg%QYtaylor,1,3)

#ifdef TRACE
      CALL UPDATER_ELEM_MOD2(dg%iotataylor,dg%iotataylor,dg%QYtaylor,1,2 )
#endif

#ifdef CHEM
      CALL UPDATER_ELEM_MOD2(dg%iotataylor,dg%iota2taylor,dg%QYtaylor,1,2 )
#endif

#endif

!.....Transform back to the Dubiner basis (multiply by dg%NmatrixInv),

      dg%ZEconst =  0.D0
      dg%QXconst =  0.D0
      dg%QYconst =  0.D0

#ifdef TRACE
      dg%iotaconst =  0.D0
#endif

#ifdef CHEM
      dg%iotaconst =  0.D0
      dg%iota2const =  0.D0
#endif

      do k=1,ne

         if (dg%dofs(k).gt.1) then
                                !do lll=1,3

            do ll=1,dg%dofs(k)

               do ss=1,dg%dofs(k)

                  dg%ZEconst(k,ll) = dg%ZEconst(k,ll) + dg%NmatrixInv(k,ll,ss,dg%dofs(k)) 
     &                 * dg%ZEtaylor(k,ss,1)
                  dg%QXconst(k,ll) = dg%QXconst(k,ll) + dg%NmatrixInv(k,ll,ss,dg%dofs(k)) 
     &                 * dg%QXtaylor(k,ss,1)
                  dg%QYconst(k,ll) = dg%QYconst(k,ll) + dg%NmatrixInv(k,ll,ss,dg%dofs(k)) 
     &                 * dg%QXtaylor(k,ss,1)

#ifdef TRACE
                  dg%iotaconst(k,ll) = dg%iotaconst(k,ll) + 
     &                 dg%NmatrixInv(k,ll,ss,dg%dofs(k)) * dg%iotataylor(k,ss,1) 
#endif

#ifdef CHEM
                  dg%iotaconst(k,ll) = dg%iotaconst(k,ll) + 
     &                 dg%NmatrixInv(k,ll,ss,dg%dofs(k)) *  dg%iotataylor(k,ss,1)
                  dg%iota2const(k,ll) = dg%iota2const(k,ll) + 
     &                 dg%NmatrixInv(k,ll,ss,dg%dofs(k)) * dg%iota2taylor(k,ss,1)
#endif

               enddo

            enddo

         elseif (dg%dofs(k).eq.1) then

            cycle

         endif

      enddo

!.....Set limit values

      do k=1,ne

         if (dg%dofs(k).gt.1) then

            do ll = 1,dg%dofs(k)

               dg%ZE(ll,k,dg%irk+1) = dg%ZEconst(k,ll)
               dg%QX(ll,k,dg%irk+1) = dg%QXconst(k,ll)
               dg%QY(ll,k,dg%irk+1) = dg%QYconst(k,ll)

#ifdef TRACE
               dg%iota(ll,k,dg%irk+1) = dg%iotaconst(k,ll)
#endif

#ifdef CHEM
               dg%iota(ll,k,dg%irk+1) = dg%iotaconst(k,ll)
               dg%iota2(ll,k,dg%irk+1) = dg%iota2const(k,ll)
#endif

            enddo

         elseif (dg%dofs(k).eq.1) then

            cycle

         endif

      enddo
      
      return
      end subroutine



!***********************************
!     *
!     SUBROUTINE Slopelim_Aux()    *
!     -cem
!     *
!***********************************

      SUBROUTINE Slopelim_Aux(ZEder,QXder,QYder,
     &     iotader,iota2der,l1,l2,l3)

!.....Use appropriate modules

      USE SIZES, ONLY : SZ
      USE GLOBAL
      USE DG

#ifdef CMPI
      USE MESSENGER
#endif

      IMPLICIT NONE

!.....Declare local variables

      INTEGER L,LL,INC1,INC2,INC3,KDP,NN,IVAR,I,J,kk,k,varnum,ss,lll,l1,l2,l3
      REAL(SZ) ZEC(3),ZEVERTEX(3),DIF(3),SUMLOC,SUMDIF,SIGNDIF,ZEVERTEX2(3),
      REAL(SZ) DIV,REDFAC,REDMAX,zek(MNE,3,1),zecc(MNE,3),zeve(MNE,3)
      Real(SZ), intent(inout) :: ZEder(MNE,dg%dofh,1)
      Real(SZ), intent(inout) :: QXder(MNE,dg%dofh,1)
      Real(SZ), intent(inout) :: QYder(MNE,dg%dofh,1)
      Real(SZ), intent(inout) :: iotader(MNE,dg%dofh,1)
      Real(SZ), intent(inout) :: iota2der(MNE,dg%dofh,1)
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



!     FIND THE MAXIMUM AND MINIMUM OF EACH VARIABLE OVER ALL ELEMENTS 
!     SHARING A NODE


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

!     IF(dg%WDFLG(NBOR_EL).EQ.0) CYCLE ! DON'T COUNT DRY ELEMENTS  sb 02/26/07

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

!     
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
!     
!     LOOP OVER ELEMENTS TO CALCULATE NEW VERTEX VALUES
!     

      zek = 0.D0
      zecc = 0.D0
      zeve = 0.D0

      DO I=1,NE 

         !IF(dg%WDFLG(I).EQ.0) CYCLE ! DON'T COUNT DRY ELEMENTS  sb 02/26/07
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

               Zhi(lll,1) =  dg%var2sigmav(i,lll,1) 
               Zhi(lll,2) =  dg%var2sigmav(i,lll,2) 
               Zhi(lll,3) =  dg%var2sigmav(i,lll,3) 

            enddo

            DO KK=1,3
               ZEVERTEX(1)=ZEVERTEX(1)+ Zhi(1,kk) * ZEC(KK)
               ZEVERTEX(2)=ZEVERTEX(2)+ Zhi(2,kk) * ZEC(KK)
               ZEVERTEX(3)=ZEVERTEX(3)+ Zhi(3,kk) * ZEC(KK)
            ENDDO

            ZEVERTEX2(1)=dg%iota(1,I,dg%IRK+1)
            ZEVERTEX2(2)=dg%iota(1,I,dg%IRK+1)
            ZEVERTEX2(3)=dg%iota(1,I,dg%IRK+1)

            DO KK=2,3
               ZEVERTEX2(1)=ZEVERTEX2(1)+ dg%PHI_CORNER(KK,1,1)*dg%iota(kk,I,dg%IRK+1)
               ZEVERTEX2(2)=ZEVERTEX2(2)+ dg%PHI_CORNER(KK,2,1)*dg%iota(kk,I,dg%IRK+1)
               ZEVERTEX2(3)=ZEVERTEX2(3)+ dg%PHI_CORNER(KK,3,1)*dg%iota(kk,I,dg%IRK+1)
            ENDDO

!     
!     RESET THE VERTEX VALUE TO BE LESS THAN OR EQUAL TO THE MAX AND
!     GREATER THAN OR EQUAL TO THE MIN AT THAT VERTEX
!     
            ZEVERTEX(1)=DMAX1(DMIN1(ZEVERTEX(1),ZEMAX1(1)),ZEMIN1(1))
            ZEVERTEX(2)=DMAX1(DMIN1(ZEVERTEX(2),ZEMAX1(2)),ZEMIN1(2))
            ZEVERTEX(3)=DMAX1(DMIN1(ZEVERTEX(3),ZEMAX1(3)),ZEMIN1(3))

!     LOOP OVER THE VERTICES 3 TIMES
!     IF THE VALUE AT THE VERTEX IS ABOVE (BELOW) THE MAX (MIN) AT THAT
!     VERTEX THEN SUBTRACT OFF THE DIFFERENCE AND ADD IT TO THE OTHER
!     VERTICES
!     
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
!     
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





!.....Subroutine to find the inverse of a square matrix by Guass-Jordan elimination
!.....cem

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
