C***********************************************************************
C     
C     SUBROUTINE NUMERICAL_FLUX
C     
C     This subroutine calls an appropriate numrical flux subroutine
C     based on FLUXTYPE variable.
C     
C     Written by Shintaro Bunya (11-08-2006)
C
C     01-10-2011 - cem - adapted for p_enrichment and multicomponent
C     06-01-2012 - cem - New fluxes, sediment, etc.
C     
C***********************************************************************

      SUBROUTINE NUMERICAL_FLUX(IT,MM)
      
      USE DG,ONLY : FLUXTYPE

      IMPLICIT NONE

      INTEGER IT,MM


      IF(FLUXTYPE.EQ.1) THEN
         CALL ROE_FLUX(IT,MM)
      ELSEIF(FLUXTYPE.EQ.2) THEN
         CALL LLF_FLUX()
      ELSEIF (FLUXTYPE.EQ.3) THEN
         CALL HLLC_FLUX()
      ELSEIF (FLUXTYPE.EQ.4) THEN
         CALL NCP_FLUX()
      ELSE
         STOP 'INVALID FLUXTYPE'
      ENDIF

      RETURN
      END SUBROUTINE
 

***********************************************************************
C     
C     SUBROUTINE ROE_FLUX
C     
C     This subroutine computes the Roe flux for the shallow water
C     equations.
C     
C     Written by Ethan Kubatko (06-11-2004)
C     
C***********************************************************************
      SUBROUTINE ROE_FLUX(IT,MM)
      
C.....Use appropriate modules

      USE GLOBAL
      USE DG
      use sizes, only: layers
      use fparser
      use fparser2

      IMPLICIT NONE
      
C.....Declare local variables
      
      INTEGER II,GED,ll,w,IT,kk,jj,MM,fast
      REAL(SZ) DEN, U_AVG, V_AVG, VEL_NORMAL, q_RoeX, q_RoeY, q_Roe
      Real(SZ) bedX_IN,bedX_EX,dtildeQx(4),dtildeQy(4)
      Real(SZ) bedY_IN,bedY_EX, bed_AVG(layers),lambda4_IN,lambda4_EX
      Real(SZ) rVncU,rVncV,lambda1_IN,lambda1_EX,lambda2_IN,lambda2_EX
      Real(SZ) lambda3_IN,lambda3_EX,Select_IN,Select_EX
      Real T1, T2, T3, T4
      Real(SZ) chi,alpha_0,beta_0,gamma1,gamma2,gamma3,gamma4
      Real(SZ) iot,iot2,iot3,iot4,iot5,J6,J7,J8,bigA,bigB
      Real(SZ) bigD,bigE,bigF,J9,epsilon,test_g,xi
      
C.....Compute the total height of the water column

      HT_IN = ZE_IN*IFNLFA + HB_IN
      HT_EX = ZE_EX*IFNLFA + HB_EX

C.....Compute the velocities in the x and y directions

      U_IN = QX_IN/HT_IN
      U_EX = QX_EX/HT_EX

      V_IN = QY_IN/HT_IN
      V_EX = QY_EX/HT_EX

C.....Compute the Roe averaged variables

      T1    = SQRT(HT_IN)
      T2    = SQRT(HT_EX)
      T3    = SQRT(HT_IN + HT_EX)
      T4    = 1.0D0/(T1 + T2)
      C_ROE = G2ROOT*T3*SQRT(NY**2 + (NX*SFAC_IN)**2) !srb - spherical coordinate correction
      U_ROE = T4*(U_IN*T1 + U_EX*T2)
      V_ROE = T4*(V_IN*T1 + V_EX*T2)

C.....Evaluate the eigenvalues at the Roe averaged variables

      EIGVAL(2) = (U_ROE*NX + V_ROE*NY)*IFNLCT
      EIGVAL(1) = EIGVAL(2) + C_ROE
      EIGVAL(3) = EIGVAL(2) - C_ROE

      !eigmax = max(abs(EIGVAL(1)),abs(EIGVAL(2)),abs(EIGVAL(3)))
      !cfl_max = max(EIGMAX,cfl_max)
      
C.....Evaluate right eigenvectors at Roe averaged variables

      RI(1,1) = 1.D0
      RI(2,1) = U_ROE*IFNLCT + (IFNLCT*C_ROE + (1-IFNLCT)*G/C_ROE)*NX*SFAC_IN !srb - spherical coordinate correction
      RI(3,1) = V_ROE*IFNLCT + (IFNLCT*C_ROE + (1-IFNLCT)*G/C_ROE)*NY

      RI(1,2) = 0.D0
      RI(2,2) = -NY
      RI(3,2) =  NX*SFAC_IN

      RI(1,3) = 1.D0
      RI(2,3) = U_ROE*IFNLCT - (IFNLCT*C_ROE + (1-IFNLCT)*G/C_ROE)*NX*SFAC_IN !srb - spherical coordinate correction
      RI(3,3) = V_ROE*IFNLCT - (IFNLCT*C_ROE + (1-IFNLCT)*G/C_ROE)*NY

C.....Evaluate left eigenvectors at Roe averaged variables

      DEN = 1.0d0/(RI(2,2)*RI(3,3)-RI(2,3)*RI(3,2)+RI(2,1)*RI(3,2)
     &     -RI(3,1)*RI(2,2))

      LE(1,1) =  DEN*( RI(2,2)*RI(3,3)-RI(2,3)*RI(3,2) )
      LE(1,2) =  DEN*RI(3,2)
      LE(1,3) = -DEN*RI(2,2)

      LE(2,1) = -DEN*( RI(2,1)*RI(3,3)-RI(2,3)*RI(3,1) )
      LE(2,2) =  DEN*( RI(3,3)-RI(3,1) )
      LE(2,3) = -DEN*( RI(2,3)-RI(2,1) )

      LE(3,1) =  DEN*( RI(2,1)*RI(3,2) - RI(2,2)*RI(3,1) )
      LE(3,2) = -DEN*RI(3,2)
      LE(3,3) =  DEN*RI(2,2)

C.....Compute the Roe matrix

      DO II = 1,3

         A_ROE(II,1) = RI(II,1)*ABS(EIGVAL(1))*LE(1,1) +
     &                 RI(II,2)*ABS(EIGVAL(2))*LE(2,1) +
     &                 RI(II,3)*ABS(EIGVAL(3))*LE(3,1)

         A_ROE(II,2) = RI(II,1)*ABS(EIGVAL(1))*LE(1,2) +
     &                 RI(II,2)*ABS(EIGVAL(2))*LE(2,2) +
     &                 RI(II,3)*ABS(EIGVAL(3))*LE(3,2)

         A_ROE(II,3) = RI(II,1)*ABS(EIGVAL(1))*LE(1,3) +
     &                 RI(II,2)*ABS(EIGVAL(2))*LE(2,3) +
     &                 RI(II,3)*ABS(EIGVAL(3))*LE(3,3)

      ENDDO

C.....Compute continuity fluxes at exterior state

      FX_EX = (IFNLCT*QX_EX + (1-IFNLCT)*QX_EX*HT_EX)*SFAC_EX
      FY_EX = IFNLCT*QY_EX + (1-IFNLCT)*QY_EX*HT_EX

C.....Compute x momentum fluxes at exterior state

      GX_EX = IFNLCT*(QX_EX*QX_EX/(HB_EX + ZE_EX) +
     &     0.5D0*G*(ZE_EX*ZE_EX + 2.D0*ZE_EX*HB_EX)) + (1-IFNLCT)*G*ZE_EX
      GX_EX = GX_EX*SFAC_EX
      GY_EX = IFNLCT*QX_EX*QY_EX/(HB_EX + ZE_EX)

C.....Compute y momentum fluxes at exterior state

      HX_EX = (IFNLCT*QX_EX*QY_EX/(HB_EX + ZE_EX))*SFAC_EX
      HY_EX = IFNLCT*(QY_EX*QY_EX/(HB_EX + ZE_EX) +
     &     0.5D0*G*(ZE_EX*ZE_EX + 2.D0*ZE_EX*HB_EX)) + (1-IFNLCT)*G*ZE_EX

C.....Compute continuity fluxes at interior state

      FX_IN = (IFNLCT*QX_IN + (1-IFNLCT)*QX_IN*HT_IN)*SFAC_IN
      FY_IN = IFNLCT*QY_IN + (1-IFNLCT)*QY_IN*HT_IN

C.....Compute x momentum fluxes at interior state

      GX_IN = IFNLCT*(QX_IN*QX_IN/(HB_IN + ZE_IN) +
     &     0.5D0*G*(ZE_IN*ZE_IN + 2.D0*ZE_IN*HB_IN)) + (1-IFNLCT)*G*ZE_IN

      GX_IN = GX_IN*SFAC_IN
      GY_IN = IFNLCT*QX_IN*QY_IN/(HB_IN + ZE_IN)

C.....Compute y momentum fluxes at interior state

      HX_IN = (IFNLCT*QX_IN*QY_IN/(HB_IN + ZE_IN))*SFAC_IN
      HY_IN = IFNLCT*(QY_IN*QY_IN/(HB_IN + ZE_IN) +
     &     0.5D0*G*(ZE_IN*ZE_IN + 2.D0*ZE_IN*HB_IN)) + (1-IFNLCT)*G*ZE_IN

C.....Compute the average flux in the normal direction

      F_AVG = 0.5D0*((FX_IN + FX_EX)*NX + (FY_IN + FY_EX)*NY)
      G_AVG = 0.5D0*((GX_IN + GX_EX)*NX + (GY_IN + GY_EX)*NY)
      H_AVG = 0.5D0*((HX_IN + HX_EX)*NX + (HY_IN + HY_EX)*NY)

c$$$      eigmax = max(abs(0.5D0*((FX_IN + FX_EX) + (FY_IN + FY_EX))),
c$$$     &     abs( 0.5D0*((GX_IN + GX_EX) + (GY_IN + GY_EX))),
c$$$     &     abs( 0.5D0*((HX_IN + HX_EX) + (HY_IN + HY_EX))) )
      
c$$$      cfl_max = max(EIGMAX,cfl_max)

C.....Compute the jump in the variables

      JUMP(1) = ZE_EX - ZE_IN
      JUMP(2) = QX_EX - QX_IN
      JUMP(3) = QY_EX - QY_IN

C.....Compute the Roe flux

      F_HAT = F_AVG - 0.5D0*(A_ROE(1,1)*JUMP(1) + A_ROE(1,2)*JUMP(2) +
     &     A_ROE(1,3)*JUMP(3) )
      G_HAT = G_AVG - 0.5D0*(A_ROE(2,1)*JUMP(1) + A_ROE(2,2)*JUMP(2) +
     &     A_ROE(2,3)*JUMP(3) )
      H_HAT = H_AVG - 0.5D0*(A_ROE(3,1)*JUMP(1) + A_ROE(3,2)*JUMP(2) +
     &     A_ROE(3,3)*JUMP(3) )

      !cfl_max = max(abs(F_HAT),abs(G_HAT),abs(H_HAT),cfl_max)

#ifdef SED_LAY

c$$$!.....Not well-balanced, be careful here!
c$$$
c$$$!.....Let us use the fully coupled system, which leads to a different eigenproblem
c$$$
c$$$      GX_IN = 0.D0
c$$$      GY_IN = 0.D0
c$$$      HX_IN = 0.D0
c$$$      HY_IN = 0.D0
c$$$      GX_EX = 0.D0
c$$$      GY_EX = 0.D0
c$$$      HX_EX = 0.D0
c$$$      HY_EX = 0.D0
c$$$      FX_IN = 0.D0
c$$$      FY_IN = 0.D0
c$$$      FX_EX = 0.D0
c$$$      Fy_EX = 0.D0
c$$$
c$$$      f_Hat = 0.D0
c$$$      g_hat = 0.D0
c$$$      h_hat = 0.D0
c$$$      bed_hat = 0.D0
c$$$
c$$$      test_g = 0.D0
c$$$
c$$$C.....Compute continuity fluxes at interior state
c$$$
c$$$      FX_IN = QX_IN*SFAC_IN
c$$$      FY_IN = QY_IN
c$$$
c$$$C.....Compute continuity fluxes at exterior state
c$$$
c$$$      FX_EX = QX_EX*SFAC_EX
c$$$      FY_EX = QY_EX
c$$$
c$$$C.....Compute x momentum fluxes at interior state
c$$$
c$$$      GX_IN = (QX_IN*QX_IN/(HT_IN) 
c$$$     &     + 0.5D0*G*(HT_IN**2))*SFAC_IN 
c$$$      !&     - (HT_IN*bed_IN(1)*G)*SFAC_IN 
c$$$      GY_IN = QX_IN*QY_IN/(HT_IN)
c$$$
c$$$C.....Compute y momentum fluxes at interior state
c$$$
c$$$      HX_IN = (QX_IN*QY_IN/(HT_IN))*SFAC_IN
c$$$      HY_IN = (QY_IN*QY_IN/(HT_IN) +
c$$$     &     0.5D0*G*(HT_IN**2))
c$$$      !&     - (HT_IN*bed_IN(1)*G)
c$$$
c$$$C.....Compute x momentum fluxes at exterior state
c$$$
c$$$      GX_EX = (QX_EX*QX_EX/(HT_EX) 
c$$$     &     + 0.5D0*G*(HT_EX**2))*SFAC_EX
c$$$      !&     - (HT_EX*bed_EX(1)*G)*SFAC_EX
c$$$      GY_EX = QX_EX*QY_EX/(HT_EX)
c$$$
c$$$C.....Compute y momentum fluxes at exterior state
c$$$
c$$$      HX_EX = (QX_EX*QY_EX/(HT_EX))*SFAC_EX
c$$$      HY_EX = (QY_EX*QY_EX/(HT_EX) +
c$$$     &     0.5D0*G*(HT_EX**2))
c$$$      !&     - (HT_EX*bed_EX(1)*G)
c$$$
c$$$
c$$$C.....Compute the average flux in the normal direction
c$$$
c$$$      F_AVG = 0.5D0*((FX_IN + FX_EX)*NX + (FY_IN + FY_EX)*NY)
c$$$      G_AVG = 0.5D0*((GX_IN + GX_EX)*NX + (GY_IN + GY_EX)*NY)
c$$$      H_AVG = 0.5D0*((HX_IN + HX_EX)*NX + (HY_IN + HY_EX)*NY)

      QMag_IN = (QX_IN*QX_IN/(HT_IN**2) + QY_IN*QY_IN/(HT_IN**2) )**(1/2)
      QMag_EX = (QX_EX*QX_EX/(HT_EX**2) + QY_EX*QY_EX/(HT_EX**2) )**(1/2)

      !Note that the choice of linearization can cause this to need updating....
      bedX_IN =  porosity*HT_IN**(-1) * QMag_IN**(2) *QX_IN *SFAC_IN
      bedX_EX =  porosity*HT_EX**(-1) * QMag_EX**(2) *QX_EX *SFAC_EX
      bedY_IN =  porosity*HT_IN**(-1) * QMag_IN**(2) *QY_IN 
      bedY_EX =  porosity*HT_EX**(-1) * QMag_EX**(2) *QY_EX

C.....Compute the remaining Roe averaged variables

      ZE_ROE  = (ZE_IN + ZE_EX) / 2.D0
      QX_ROE  = (QX_IN + QX_EX) / 2.D0
      QY_ROE  = (QY_IN + QY_EX) / 2.D0
      bed_ROE = (bed_IN(1)+bed_EX(1)) / 2.D0


      !if (it.eq.1.and.mm.eq.1) then !do this only once, because it is Slooooow ....
      if (.not.init_parser) then

         open(446, file = "./db_partials_X")
         read(446,'(a)') funcx(1)
         read(446,'(a)') funcx(2)
         read(446,'(a)') funcx(3)
         read(446,'(a)') funcx(4)
         close(446)
         open(447, file = "./db_partials_Y")
         read(447,'(a)') funcy(1)
         read(447,'(a)') funcy(2)
         read(447,'(a)') funcy(3)
         read(447,'(a)') funcy(4)
         close(447)

         !Need to tokenize the sed functions
         CALL initf (4)         ! Initialize function parser for nfunc functions
         CALL initf2 (4)        ! Initialize function parser for nfunc functions
         do w=1,4               ! Loop over functions
            CALL parsef (w, funcx(w), varx) ! Parse and bytecompile w-th function string 
            CALL parsef2 (w, funcy(w), vary) ! Parse and bytecompile w-th function string 
         enddo

         init_parser = .true.

      endif
      
      valx(1) = ZE_ROE
      valx(2) = QX_ROE
      valx(3) = QY_ROE
      valx(4) = bed_ROE

      valy(1) = ZE_ROE
      valy(2) = QX_ROE
      valy(3) = QY_ROE
      valy(4) = bed_ROE

      dtildeQx = 0.0
      dtildeQy = 0.0
      do w=1,4
         dtildeQx(w) = porosity * evalf (w, valx)  *( QMag_IN**(2) + QMag_EX**(2) ) / 2.D0 ! Interprete bytecode representation of w-th function
         dtildeQy(w) = porosity*evalf2 (w, valy) *( QMag_IN**(2) + QMag_EX**(2) ) / 2.D0 ! Interprete bytecode representation of w-th function
c$$$         if (EvalErrType > 0.or.EvalErrType2 > 0) then
c$$$            if(porosity*((3*QX_ROE**2 + QX_ROE*QY_ROE + QY_ROE**2/2)/(ZE_ROE + bed_ROE)**3).ne.dtildeQx(w)) then
c$$$               write(*,*)'That was foolish on timestep', it
c$$$               WRITE(*,*)'*** Error: ',EvalErrMsg ()
c$$$               WRITE(*,*)'*** Error: ',EvalErrMsg2 ()
c$$$               WRITE(*,*) funcx(w),'X=',dtildeQx(w),w
c$$$               WRITE(*,*) funcy(w),'Y=',dtildeQy(w),w
c$$$               WRITE(*,*) 'valx = ', valx(w), w
c$$$               WRITE(*,*) 'valy = ', valy(w), w
c$$$               print*, '((3*QX_ROE**2 + QX_ROE*QY_ROE + QY_ROE**2/2)/(ZE_ROE + bed_ROE)**3)
c$$$     &              actualX=', porosity*((3*QX_ROE**2 + QX_ROE*QY_ROE + QY_ROE**2/2)/(ZE_ROE + bed_ROE)**3)
c$$$               print*, '(QX_ROE*QY_ROE + QY_ROE**2/2)/(ZE_ROE + bed_ROE)**3  
c$$$     &              actualY=', porosity*((QX_ROE*QY_ROE + QY_ROE**2/2)/(ZE_ROE + bed_ROE)**3 ) 
c$$$               write(*,*)'Ya shoulda known bettah!'    
c$$$            endif
c$$$         endif
      enddo


      
      EIGVAL(4) = 0.D0
      do w=1,4
         EIGVAL(4) = EIGVAL(4) + dtildeQx(w)*NX + dtildeQy(w)*NY
      enddo

C.....A bunch of functions

      epsilon = 0.D0 !1.0e-12

      chi = (NY**2) * dtildeQy(2) + dtildeQx(3)*(NY**2) - dtildeQx(3) 
     &     + (dtildeQx(2)-dtildeQy(3))*NY*NX

      xi = C_ROE**2

      alpha_0 = NX*(NY**2)*(dtildeQx(3)*xi + dtildeQy(2)*xi + dtildeQy(2)*(V_ROE**2))
     &     + (NY**3) * (dtildeQy(3)*xi - dtildeQy(2)*U_ROE*V_ROE - dtildeQx(2) * xi )
     &     + NY*(-dtildeQx(2)*(EIGVAL(2)**2)-dtildeQy(3)*(V_ROE**2) + dtildeQx(2)*xi
     &     - dtildeQy(1)*V_ROE + (V_ROE**2) * dtildeQx(2)) + dtildeQx(2)*(-(NY**2) * V_ROE 
     &     + NY*EIGVAL(2))*EIGVAL(1) + (-dtildeQy(2)*V_ROE*EIGVAL(2)-dtildeQy(2)*V_ROE*C_ROE
     &     + dtildeQx(1)*C_ROE)*NY*NX+(-dtildeQx(1)*V_ROE-dtildeQx(3)*(V_ROE**2))*NX
     &     + (dtildeQy(1)*C_ROE+dtildeQy(2)*C_ROE*U_ROE+dtildeQx(2)*V_ROE*EIGVAL(2)
     &     + dtildeQx(2)*V_ROE*C_ROE)*(NY**2) - dtildeQx(2)*V_ROE*EIGVAL(1)

      beta_0 = NX*(NY**2) *(dtildeQx(3)*xi+dtildeQy(2)*xi + dtildeQy(2)*(V_ROE**2))
     &     + (NY**3) * (dtildeQy(3)*xi - dtildeQy(2)*U_ROE*V_ROE - dtildeQx(2) * xi )
     &     + NY*(-dtildeQx(2)*(EIGVAL(2)**2)-dtildeQy(3)*(V_ROE**2) + dtildeQx(2)*xi
     &     - dtildeQy(1)*V_ROE + (V_ROE**2) * dtildeQx(2)) + dtildeQx(2)*(-(NY**2) * V_ROE 
     &     + NY*EIGVAL(2))*EIGVAL(3) + (-dtildeQy(2)*V_ROE*EIGVAL(2)-dtildeQy(2)*V_ROE*C_ROE
     &     + dtildeQx(1)*C_ROE)*NY*NX+(-dtildeQx(1)*V_ROE-dtildeQx(3)*(V_ROE**2))*NX
     &     + (-dtildeQy(1)*C_ROE-dtildeQy(2)*C_ROE*U_ROE+dtildeQx(2)*V_ROE*EIGVAL(2)
     &     - dtildeQx(2)*V_ROE*C_ROE)*(NY**2) + dtildeQx(2)*(-EIGVAL(2)+C_ROE)*V_ROE

      gamma1 = NY * (-(V_ROE**2) * EIGVAL(2) - (V_ROE**2) * C_ROE + (C_ROE**3) - EIGVAL(2)*xi )
     &    + V_ROE * xi + ( -EIGVAL(2) - C_ROE ) * (-NY*(V_ROE**2) + EIGVAL(2)*V_ROE )
     &    + (EIGVAL(2)**2) * V_ROE

      gamma2 = NY * (-(V_ROE**2) * EIGVAL(2) + (V_ROE**2) * C_ROE - (C_ROE**3) - EIGVAL(2)*xi )
     &    + V_ROE * xi + ( -EIGVAL(2) + C_ROE ) * (-NY*(V_ROE**2) + EIGVAL(2)*V_ROE )
     &    + (EIGVAL(2)**2) * V_ROE

      gamma3 = - NX*NY*V_ROE*xi + (NY**2) * U_ROE * xi

      gamma4 = C_ROE*EIGVAL(2)*U_ROE - NX*(C_ROE**3)

      iot = gamma2*EIGVAL(3)*NY+NY*gamma1*EIGVAL(1)+NX*EIGVAL(3)*gamma4+NX*EIGVAL(3)
     &     *gamma3-NX*EIGVAL(1)*gamma4+NX*EIGVAL(1)*gamma3

      iot2 = - (NY**2) * xi + (V_ROE**2)

      iot3 = -chi*(EIGVAL(2)**2) - chi*(EIGVAL(4)**2) + 2.D0*chi*EIGVAL(2)*EIGVAL(4)

      iot4 = beta_0*( (EIGVAL(4)**2) - 2.D0*EIGVAL(4)*EIGVAL(2) - C_ROE*EIGVAL(4)
     &     + (EIGVAL(2)**2) +C_ROE*EIGVAL(2) )

      iot5 = -alpha_0*(EIGVAL(2)**2) + alpha_0*C_ROE*EIGVAL(2) - alpha_0*(EIGVAL(4)**2)
     &     + 2.D0*alpha_0*EIGVAL(2)*EIGVAL(4) - alpha_0*C_ROE*EIGVAL(4)

      J6 = (NY**2) * (xi**2) + iot2*xi

      J7 = iot3*xi + chi * (xi**2)

      J8 = beta_0*(EIGVAL(2)**2) - alpha_0*(EIGVAL(2)**2) - alpha_0*C_ROE*EIGVAL(4)
     &     - 2.D0*beta_0*EIGVAL(2)*EIGVAL(4) + 2.D0*alpha_0*EIGVAL(2)*EIGVAL(4) + alpha_0*
     &     C_ROE*EIGVAL(2) + beta_0*C_ROE*EIGVAL(2) - alpha_0*(EIGVAL(4)**2) + beta_0
     &     *(EIGVAL(4)**2) - beta_0*C_ROE*EIGVAL(4)

      J9 = chi*J6 + iot3*(NY**2)*xi + iot2*iot3

      bigA = EIGVAL(4) + C_ROE - EIGVAL(2)

      bigB = EIGVAL(2) + C_ROE - EIGVAL(4)

      bigD = -gamma2*gamma4*J9+gamma2*gamma3*J9-gamma4*gamma1*J9-gamma3*gamma1*J9
     &     +gamma1*gamma2*NY*J8+(gamma3*gamma1*J7+gamma4*gamma1*J7-gamma2*gamma3*J7
     &     +gamma2*gamma4*J7)*(NY**2)+gamma3*NX*gamma1*iot5+gamma4*NX*gamma1*iot5
     &     +NX*gamma3*gamma2*iot4-NX*gamma4*gamma2*iot4

      bigE = -gamma1*EIGVAL(1)*J9-EIGVAL(3)*gamma2*J9+(-NX*EIGVAL(1)*gamma3*J7-NX*EIGVAL(3)
     &     *gamma3*J7+NX*EIGVAL(1)*gamma4*J7-NX*EIGVAL(3)*gamma4*J7+iot*J7)*NY
     &     +NX*EIGVAL(1)*gamma1*iot5-NX*EIGVAL(3)*gamma2*iot4

      bigF = EIGVAL(1)*gamma3*J9+EIGVAL(3)*gamma3*J9-EIGVAL(1)*gamma4*J9+EIGVAL(3)*gamma4
     &     *J9-gamma2*EIGVAL(3)*NY*J8+(-EIGVAL(1)*gamma3*J7-EIGVAL(3)*gamma3*J7+EIGVAL(1)
     &     *gamma4*J7-EIGVAL(3)*gamma4*J7)*(NY**2)-NX*EIGVAL(1)*gamma3*iot5-NX*EIGVAL(3)
     &     *gamma3*iot5+NX*EIGVAL(1)*gamma4*iot5-NX*EIGVAL(3)*gamma4*iot5+iot*iot5

      RI = 0.D0
      LE = 0.D0
      A_ROE = 0.D0
    
C.....Evaluate right eigenvectors at Roe averaged variables

      if (abs(EIGVAL(1)*alpha_0*gamma1).le.epsilon) then

         RI(1,1) = 0.D0

      else

         RI(1,1) = ( -C_ROE*(EIGVAL(2)**2) + (C_ROE**3) ) * ( EIGVAL(2) + C_ROE - EIGVAL(4) ) 
     &        *((NY**2) * xi - (V_ROE**2) ) / (EIGVAL(1)*alpha_0*gamma1)

      endif

      if (abs(alpha_0*gamma1).le.epsilon) then

         RI(2,1)= 0.D0

      else

         RI(2,1) = ((NY**2) * xi - (V_ROE**2)) * ( EIGVAL(2) + C_ROE - EIGVAL(4) )
     &        *( -gamma4 + gamma3 ) / (alpha_0*gamma1)

      endif

      if (abs(alpha_0).le.epsilon) then

         RI(3,1) = 0.D0
      
      else
         
         RI(3,1) = (EIGVAL(2)+C_ROE-EIGVAL(4))*( (NY*C_ROE)**2 - V_ROE**2 )/ alpha_0

      endif

      RI(4,1) = 1.D0

      RI(1,2) = 0.D0

      if (abs(chi).le.epsilon) then

         RI(2,2) = 0.D0
         RI(3,2) = 0.D0

      else

         RI(2,2) = (EIGVAL(2) - EIGVAL(4)) * NY / chi
         RI(3,2) = (EIGVAL(4) - EIGVAL(2)) * NX / chi

      endif

      RI(4,2) = 1.D0


      if ( abs(EIGVAL(3)*beta_0*gamma2).le.epsilon) then

         RI(1,3) = 0.D0

      else

         RI(1,3) = (C_ROE*EIGVAL(2)**2 - C_ROE**3)*(EIGVAL(2) - C_ROE - EIGVAL(4))
     &        * ((NY**2) * xi - V_ROE**2) / (EIGVAL(3)*beta_0*gamma2)

      endif

      if (abs(beta_0*gamma2).le.epsilon) then

         RI(2,3) = 0.D0

      else

         RI(2,3) = ((NY**2) * xi - V_ROE**2)*(EIGVAL(2)-C_ROE-EIGVAL(4))
     &        *(gamma4 + gamma3 ) / (beta_0*gamma2)

      endif

      if (abs(beta_0).le.epsilon) then
         
         RI(3,3) = 0.D0
      
      else
         
         RI(3,3) = (EIGVAL(2)-C_ROE-EIGVAL(4))*((NY**2) * xi - V_ROE**2) / beta_0
         
      endif

      RI(4,3) = 1.D0

      RI(1,4) = 0.D0
      RI(2,4) = 0.D0
      RI(3,4) = 0.D0
      RI(4,4) = 1.D0

      
C.....Evaluate left eigenvectors (inverse) at Roe averaged variables

      if ( abs(C_ROE*(EIGVAL(2)+C_ROE-EIGVAL(4))*iot*(-(EIGVAL(2)**2)+xi)*iot2)
     &  .le.epsilon) then

         LE(1,1) = 0.D0

      else

         LE(1,1) =  - ( alpha_0*EIGVAL(3)*EIGVAL(1)*gamma1*(NY*gamma2+NX*gamma4+NX*gamma3))
     &        / (C_ROE*bigB*iot*iot2*(-(EIGVAL(2)**2)+xi))


      endif

      if ( abs(C_ROE*(-(EIGVAL(2)**2)+xi)*(EIGVAL(2)-EIGVAL(4))*iot).le.epsilon) then

         LE(2,1) = 0.D0

      else

         LE(2,1) = chi*EIGVAL(1)*EIGVAL(3)*(gamma2*gamma4-gamma2*gamma3+gamma4*gamma1
     &        +gamma3*gamma1) / (C_ROE*(-(EIGVAL(2)**2)+xi)*(EIGVAL(2)-EIGVAL(4))*iot)


      endif

      if ( abs(C_ROE*bigA*iot*(-(EIGVAL(2)**2)+xi)*iot2).le.epsilon) then

         LE(3,1) = 0.D0

      else

         LE(3,1) = -beta_0*EIGVAL(1)*EIGVAL(3)*gamma2*(gamma1*NY-NX*gamma4+NX*gamma3)
     &        / (C_ROE*bigA*iot*(-(EIGVAL(2)**2)+xi)*iot2)

      endif

      if ( abs(C_ROE*iot*bigA*bigB*(EIGVAL(2)-EIGVAL(4))*iot2
     & *(-(EIGVAL(2)**2)+xi)).le.epsilon) then

         LE(4,1) = 0.D0
         
      else
      
         LE(4,1) = (EIGVAL(1)*EIGVAL(3)*bigD) / (C_ROE*iot*iot2*bigA*bigB 
     &        *(EIGVAL(2)-EIGVAL(4))*(-(EIGVAL(2)**2)+xi))

      endif

      if ( abs( bigB*iot*iot2 ).le.epsilon) then

         LE(1,2) = 0.D0

      else

         LE(1,2) = -alpha_0*gamma1*EIGVAL(1)*NX / ( bigB*iot*iot2 )

      endif

      if ( abs((EIGVAL(2)-EIGVAL(4))*iot).le.epsilon) then
      
         LE(2,2) = 0.D0

      else

         LE(2,2) = chi * ( EIGVAL(3)*gamma2 + gamma1*EIGVAL(1) ) / ((EIGVAL(2)-EIGVAL(4))*iot)

      endif

      if ( abs(bigA*iot*iot2).le.epsilon) then

         LE(3,2) = 0.D0
         
      else

         LE(3,2) = beta_0*gamma2*EIGVAL(3)*NX / (bigA*iot*iot2)

      endif

      if ( abs(iot*iot2*bigA*bigB*(EIGVAL(2)-EIGVAL(4))).le.epsilon) then

         LE(4,2) = 0.D0

      else

         LE(4,2) = bigE / (iot*iot2*bigA*bigB*(EIGVAL(2)-EIGVAL(4)))

      endif

      if ( abs(bigB*iot*iot2).le.epsilon) then

         LE(1,3) = 0.D0
         
      else

         LE(1,3) = -alpha_0*(-gamma2*EIGVAL(3)*NY-NX*EIGVAL(3)*gamma4-NX*EIGVAL(3)*gamma3
     &        +NX*EIGVAL(1)*gamma4-NX*EIGVAL(1)*gamma3+iot) / (bigB*iot*iot2)

      endif

      if ( abs((EIGVAL(2)-EIGVAL(4))*iot).le.epsilon) then
         
         LE(2,3) = 0.D0
         
      else

         LE(2,3) = -chi*(EIGVAL(3)*gamma4+EIGVAL(3)*gamma3-EIGVAL(1)*gamma4+EIGVAL(1)*gamma3)
     &        / ((EIGVAL(2)-EIGVAL(4))*iot)


      endif


      if ( abs(bigA*iot*iot2).le.epsilon) then 

         LE(3,3) = 0.D0
         
      else

         LE(3,3) = ( beta_0*gamma2*EIGVAL(3)*NY ) / (bigA*iot*iot2)

      endif

      if ( abs(iot*iot2*bigA*bigB*(EIGVAL(2)-EIGVAL(4))).le.epsilon) then

         LE(4,3) = 0.D0

      else

         LE(4,3) = bigF / (iot*iot2*bigA*bigB*(EIGVAL(2)-EIGVAL(4)))
      
      endif

      LE(1,4) = 0.D0
      LE(2,4) = 0.D0
      LE(3,4) = 0.D0
      LE(4,4) = 1.D0

C.....Compute the Roe matrix


       !Let's get rid of the bad stuff in the limit

      do ii=1,4
         if(abs(EIGVAL(ii)).le.1.0e-14) then

            RI(:,ii) = 0.D0
            LE(ii,:) = 0.D0

         endif
      enddo


      DO II = 1,4

         A_ROE(II,1) = RI(II,1)*ABS(EIGVAL(1))*LE(1,1) +
     &        RI(II,2)*ABS(EIGVAL(2))*LE(2,1) +
     &        RI(II,3)*ABS(EIGVAL(3))*LE(3,1) +
     &        RI(II,4)*ABS(EIGVAL(4))*LE(4,1) 

         A_ROE(II,2) = RI(II,1)*ABS(EIGVAL(1))*LE(1,2) +
     &        RI(II,2)*ABS(EIGVAL(2))*LE(2,2) +
     &        RI(II,3)*ABS(EIGVAL(3))*LE(3,2) +
     &        RI(II,4)*ABS(EIGVAL(4))*LE(4,2)

         A_ROE(II,3) = RI(II,1)*ABS(EIGVAL(1))*LE(1,3) +
     &        RI(II,2)*ABS(EIGVAL(2))*LE(2,3) +
     &        RI(II,3)*ABS(EIGVAL(3))*LE(3,3) +
     &        RI(II,4)*ABS(EIGVAL(4))*LE(4,3)

         A_ROE(II,4) = RI(II,1)*ABS(EIGVAL(1))*LE(1,4) +
     &        RI(II,2)*ABS(EIGVAL(2))*LE(2,4) +
     &        RI(II,3)*ABS(EIGVAL(3))*LE(3,4) + 
     &        RI(II,4)*ABS(EIGVAL(4))*LE(4,4)

      ENDDO

      !Assume the bed load

      bed_AVG(1) = 0.5D0*( (bedX_IN + bedX_EX)*NX + (bedY_IN + bedY_EX)*NY )

      JUMP(4) = bed_EX(1) - bed_IN(1)

      !Let's get rid of the bad stuff in the limit

      do ii=1,4
         if(abs(jump(ii)).le.1.0e-14) then

            A_ROE(:,ii) = 0.D0

         endif
      enddo


! Compute correction (this is not set by default)

c$$$      lambda1_IN = 0.D0
c$$$      lambda1_EX = 0.D0
c$$$      lambda2_IN = 0.D0
c$$$      lambda2_EX = 0.D0
c$$$      lambda3_IN = 0.D0
c$$$      lambda3_EX = 0.D0
c$$$
c$$$      lambda2_IN = (U_IN*NX+V_IN*NY)
c$$$      lambda2_EX = (U_EX*NX+V_EX*NY)
c$$$
c$$$      lambda1_IN = lambda2_IN + sqrt(g*HT_IN)
c$$$      lambda1_EX = lambda2_EX + sqrt(g*HT_EX)
c$$$
c$$$      lambda3_IN = lambda2_IN - sqrt(g*HT_IN)
c$$$      lambda3_EX = lambda2_EX - sqrt(g*HT_EX)
c$$$
c$$$      lambda4_IN = 0.D0
c$$$      lambda4_EX = 0.D0
c$$$
c$$$      lambda4_IN = (-3.D0*QX_IN/(ZE_IN + bed_IN(1))**4 + (ZE_IN + bed_IN(1))**(-3)
c$$$     &     -3.D0*QX_IN/(ZE_IN + bed_IN(1))**4)*QMag_IN*NX 
c$$$     &     + (-3.D0*QY_IN/(ZE_IN + bed_IN(1))**4 + (ZE_IN + bed_IN(1))**(-3)
c$$$     &     -3.D0*QY_IN/(ZE_IN + bed_IN(1))**4)*QMag_IN*NY 
c$$$
c$$$      lambda4_EX = (-3.D0*QX_EX/(ZE_EX + bed_EX(1))**4 + (ZE_EX + bed_EX(1))**(-3)
c$$$     &     -3.D0*QX_EX/(ZE_EX + bed_EX(1))**4)*QMag_EX*NX 
c$$$     &     + (-3.D0*QY_EX/(ZE_EX + bed_EX(1))**4 + (ZE_EX + bed_EX(1))**(-3)
c$$$     &     -3.D0*QY_EX/(ZE_EX + bed_EX(1))**4)*QMag_EX*NY 
c$$$
c$$$
c$$$
c$$$      Select_IN = min(min(lambda1_IN,lambda2_IN,lambda3_IN,lambda4_IN),
c$$$     &                min(lambda1_EX,lambda2_EX,lambda3_EX,lambda4_EX))
c$$$      Select_EX = max(max(lambda1_IN,lambda2_IN,lambda3_IN,lambda4_IN),
c$$$     &                max(lambda1_EX,lambda2_EX,lambda3_EX,lambda4_EX))
c$$$
c$$$      rVncU = 0.D0
c$$$      rVncV = 0.D0
c$$$
c$$$      do ll=1,layers !should be fixed for multiple layers
c$$$
c$$$         rVncU =  rVncU + 0.5D0*G*(HT_IN+HT_EX) 
c$$$     &        *(bed_EX(ll)-bed_IN(ll))*NX
c$$$         rVncV =  rVncV + 0.5D0*G*(HT_IN+HT_EX) 
c$$$     &        *(bed_EX(ll)-bed_IN(ll))*NY
c$$$         
c$$$        ! print*,  rVncU, rVncV
c$$$
c$$$      enddo
c$$$



C.....Compute the Roe flux (only formed for one sediment layer!) !Can add path conservative contribution

      F_HAT = F_AVG  - 0.5D0*(A_ROE(1,1)*JUMP(1) + A_ROE(1,2)*JUMP(2) +
     &     A_ROE(1,3)*JUMP(3) + A_ROE(1,4) * JUMP(4) ) 
      !&        - (Select_IN*Select_EX)*(HT_IN-HT_EX) / (Select_EX-Select_IN)

      !print*, F_Hat

      G_HAT = G_AVG - 0.5D0*(A_ROE(2,1)*JUMP(1) + A_ROE(2,2)*JUMP(2) +
     &     A_ROE(2,3)*JUMP(3) + A_ROE(2,4) * JUMP(4) )
      !&        - 0.5D0*(Select_IN+Select_EX)*rVncU
      !&        / (Select_EX-Select_IN)


      H_HAT = H_AVG - 0.5D0*(A_ROE(3,1)*JUMP(1) + A_ROE(3,2)*JUMP(2) +
     &     A_ROE(3,3)*JUMP(3)+ A_ROE(3,4) * JUMP(4) )
      !&        - 0.5D0*(Select_IN+Select_EX)*rVncV
      !&        / (Select_EX-Select_IN)

      bed_HAT(1) = bed_AVG(1) - 0.5D0*( A_ROE(4,1)*JUMP(1) + A_ROE(4,2)*JUMP(2) +
     &     A_ROE(4,3)*JUMP(3) + A_ROE(4,4) * JUMP(4) )
      !&           - ((Select_IN*Select_EX) * (bed_IN(1) - bed_EX(1)))
      !&           / (Select_EX-Select_IN) 
      
#endif

C.....Compute upwinding scheme for chem and tracer

#ifdef TRACE
      U_AVG = 0.5D0 * ( U_IN + U_EX )
      V_AVG = 0.5D0 * ( V_IN + V_EX )
      VEL_NORMAL = U_AVG*NX + V_AVG*NY
      
      if (VEL_NORMAL.GT.0) THEN
         I_HAT = VEL_NORMAL*iota_IN
      else
         I_HAT = VEL_NORMAL*iota_EX
      endif
#endif

#ifdef CHEM
      U_AVG = 0.5D0 * ( U_IN + U_EX )
      V_AVG = 0.5D0 * ( V_IN + V_EX )
      VEL_NORMAL = U_AVG*NX + V_AVG*NY
      
      if (VEL_NORMAL.GT.0) THEN
         I_HAT = VEL_NORMAL*iota_IN
         J_HAT = VEL_NORMAL*iota2_IN
      else
         I_HAT = VEL_NORMAL*iota_EX
         J_HAT = VEL_NORMAL*iota2_EX
      endif
#endif

#ifdef DYNP
      U_AVG = 0.5D0 * ( U_IN + U_EX )
      V_AVG = 0.5D0 * ( V_IN + V_EX )
      VEL_NORMAL = U_AVG*NX + V_AVG*NY
      
      if (VEL_NORMAL.GT.0) THEN
         K_HAT = VEL_NORMAL*dynP_IN
      else
         K_HAT = VEL_NORMAL*dynP_EX
      endif
#endif

      RETURN
      END SUBROUTINE


C***********************************************************************
C     
C     SUBROUTINE LLF_FLUX
C     
C     This subroutine computes the Local Lax Friedrichs flux for the
C     shallow water equations.
C     
C     Written by Ethan Kubatko (07-15-2005)
C     
C***********************************************************************

      SUBROUTINE LLF_FLUX()

C.....Use appropriate modules

      USE GLOBAL
      USE DG
      use sizes, only: layers
!!!!!!!!!!!!!!!!! ASD !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      USE pressure_module, ONLY: pd_in,pd_ex
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      

      IMPLICIT NONE


C.....Declare local variables

      INTEGER II,l
      REAL(SZ) EIGVALS(6)
      REAL(SZ) EIGMAX, IDEPTH, HUU, HUV, HVV, GH2
      REAL(SZ) C_EX, C_IN, UN_EX, UN_IN, U_AVG, V_AVG, VEL_NORMAL
      REAL(SZ) DEPTH, F1_NL, FU_NL, FV_NL, FG_NL, FH_NL
      REAL(SZ) FX1_IN, FX2_IN, FX3_IN, FY1_IN, FY2_IN, FY3_IN
      REAL(SZ) FX1_EX, FX2_EX, FX3_EX, FY1_EX, FY2_EX, FY3_EX
      REAL(SZ) F1_AVG, F2_AVG, F3_AVG, NX2, NY2, NXY
      REAL(SZ) TOL
      REAL(SZ) DEN, q_RoeX, q_RoeY, q_Roe
      Real(SZ) discharge_modelX_IN,discharge_modelX_EX
      Real(SZ) discharge_modelY_IN,discharge_modelY_EX
      Real(SZ) rVncU,rVncV,lambda1_IN,lambda1_EX,lambda2_IN,lambda2_EX
      Real(SZ) lambda3_IN,lambda3_EX,Select_IN,Select_EX

C.....Compute the jump in the variables

      JUMP(1) = ZE_EX - ZE_IN
      JUMP(2) = QX_EX - QX_IN
      JUMP(3) = QY_EX - QY_IN

C.....Compute the total height of the water column

      HT_IN = ZE_IN*NLEQ + HB_IN

C.....Compute continuity fluxes at interior state

      F1_NL  = NLEQ + LEQ*HB_IN
      FX1_IN = QX_IN*F1_NL*SFAC_IN
      FY1_IN = QY_IN*F1_NL

C.....Compute momentum flux terms at interior state

      FU_NL = NLEQ*QX_IN
      FV_NL = NLEQ*QY_IN
      FG_NL = NLEQG*ZE_IN
      FH_NL = 1.D0/(NLEQ*HT_IN + LEQ)
      U_IN  = QX_IN*FH_NL
      V_IN  = QY_IN*FH_NL

      HUU = FU_NL*U_IN
      HVV = FV_NL*V_IN
      HUV = FU_NL*V_IN
      GH2 = FG_NL*(0.5D0*ZE_IN + HB_IN) + FG_L*ZE_IN

C.....Compute x momentum fluxes at interior state

      FX2_IN = (HUU + GH2)*SFAC_IN
      FY2_IN = HUV

C.....Compute y momentum fluxes at interior state

      FX3_IN = HUV*SFAC_IN
      FY3_IN = HVV + GH2
      
!!!!!!!!!!!!!!!!! ASD !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!          Include pressure terms in momentum flux
      FX2_IN = FX2_IN + pd_in*SFAC_IN
      FY3_IN = FY3_IN + pd_in
      
C.....Compute the total height of the water column

      HT_EX = ZE_EX*NLEQ + HB_EX

C.....Compute continuity fluxes at exterior state

      F1_NL  = NLEQ + LEQ*HB_EX
      FX1_EX = QX_EX*F1_NL*SFAC_EX
      FY1_EX = QY_EX*F1_NL

C.....Compute momentum flux terms at interior state

      FU_NL = NLEQ*QX_EX
      FV_NL = NLEQ*QY_EX
      FG_NL = NLEQG*ZE_EX
      FH_NL = 1.D0/(NLEQ*HT_EX + LEQ)
      U_EX  = QX_EX*FH_NL
      V_EX  = QY_EX*FH_NL

      HUU = FU_NL*U_EX
      HVV = FV_NL*V_EX
      HUV = FU_NL*V_EX
      GH2 = FG_NL*(0.5D0*ZE_EX + HB_EX) + FG_L*ZE_EX

C.....Compute x momentum fluxes at exterior state

      FX2_EX = (HUU + GH2)*SFAC_EX
      FY2_EX = HUV

C.....Compute y momentum fluxes at exterior state

      FX3_EX = HUV*SFAC_EX
      FY3_EX = HVV + GH2
!!!!!!!!!!!!!!!!! ASD !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!          Include pressure terms in momentum flux
      FX2_EX = FX2_EX + pd_ex*SFAC_EX
      FY3_EX = FY3_EX + pd_ex

C.....Compute the average flux function

      F1_AVG = 0.5D0*((FX1_IN + FX1_EX)*NX + (FY1_IN + FY1_EX)*NY)
      F2_AVG = 0.5D0*((FX2_IN + FX2_EX)*NX + (FY2_IN + FY2_EX)*NY)
      F3_AVG = 0.5D0*((FX3_IN + FX3_EX)*NX + (FY3_IN + FY3_EX)*NY)

C.....Evaluate the eigenvalues at the interior and exterior states

      UN_IN = (U_IN*NX + V_IN*NY)*NLEQ
      C_IN = SQRT(G*HT_IN*(NY**2 + (NX*SFAC_IN)**2)) !srb - spherical coordinate correction
      EIGVALS(1) = ABS(UN_IN + C_IN)
      EIGVALS(2) = ABS(UN_IN)
      EIGVALS(3) = ABS(UN_IN - C_IN)

      UN_EX = (U_EX*NX + V_EX*NY)*NLEQ
      C_EX = SQRT(G*HT_EX*(NY**2 + (NX*SFAC_EX)**2)) !srb - spherical coordinate correction
      EIGVALS(4) = ABS(UN_EX + C_EX)
      EIGVALS(5) = ABS(UN_EX)
      EIGVALS(6) = ABS(UN_EX - C_EX)

C.....Find the maximum eigenvalue (in absolute value)

      EIGMAX = MAX( EIGVALS(1), EIGVALS(2), EIGVALS(3),
     &     EIGVALS(4), EIGVALS(5), EIGVALS(6) )

      !cfl_max = max(EIGMAX,cfl_max)
                                !print*,eigmax

C.....Compute the Local Lax Friedrichs Fluxes

      NX2 = NLEQ+LEQ*NX*NX
      NY2 = NLEQ+LEQ*NY*NY
      NXY = NX*NY

      F_HAT = F1_AVG - 0.5D0*EIGMAX*(JUMP(1))
      G_HAT = F2_AVG - 0.5D0*EIGMAX*(JUMP(2)*NX2 + LEQ*JUMP(3)*NXY)
      H_HAT = F3_AVG - 0.5D0*EIGMAX*(JUMP(3)*NY2 + LEQ*JUMP(2)*NXY)

#ifdef SED_LAY !this is not well-tested!

!.....Compute the stabilization term (v_nc in CMM's thesis) 
!.....for the nonconservative product in the momentum eqn

                                !Compute the magnitudes of the velocities for sediment

      QMag_IN = (QX_IN*QX_IN/(HT_IN**2) + QY_IN*QY_IN/(HT_IN**2) )**(1/2)
      QMag_EX = (QX_EX*QX_EX/(HT_EX**2) + QY_EX*QY_EX/(HT_EX**2) )**(1/2)

      rVncU = 0.D0
      rVncV = 0.D0

      do l=1,layers

         rVncU =  rVncU + 0.5D0*G*(HT_IN+HT_EX) 
     &        *(bed_EX(l)-bed_IN(l))*NX
         rVncV =  rVncV + 0.5D0*G*(HT_IN+HT_EX) 
     &        *(bed_EX(l)-bed_IN(l))*NY

      enddo

!.....When not using the NCP formalism
!.....then we enforce the (truncated) HLL wave-speed rendition.
!.....See CMM's thesis, pg 45, substituting the native flux

      lambda1_IN = U_IN*NX + V_IN*NY
      lambda1_EX = U_EX*NX + V_EX*NY

      lambda2_IN = lambda1_IN - (g*HT_IN)**(1/2)
      lambda2_EX = lambda1_EX - (g*HT_EX)**(1/2)

      lambda3_IN = lambda1_IN + (g*HT_IN)**(1/2)
      lambda3_EX = lambda1_EX + (g*HT_EX)**(1/2)

      Select_IN = min(lambda2_IN,lambda2_EX)
      Select_EX = max(lambda3_IN,lambda3_EX)

!.....Computes the sediment discharge rate based on the formulation
!.....developed by Camenen and Larson, so $\tilde{q} = A_g*H^{-3}|q|^{2}q$
!.....which clearly makes sense only for a single layer

      discharge_modelX_IN =  porosity*HT_IN**(-1)*QMag_IN**(2)*QX_IN *SFAC_IN
      discharge_modelX_EX =  porosity*HT_EX**(-1)*QMag_EX**(2)*QX_EX *SFAC_EX
      discharge_modelY_IN =  porosity*HT_IN**(-1)*QMag_IN**(2)*QY_IN 
      discharge_modelY_EX =  porosity*HT_EX**(-1)*QMag_EX**(2)*QY_EX



      q_RoeX = ( QX_IN*(HT_IN)**0.5D0  
     &     + QX_EX*(HT_EX)**0.5D0 ) / ((HT_IN)**0.5D0 + (HT_EX)**0.5D0)
      q_RoeY = ( QY_IN*(HT_IN)**0.5D0  
     &     + QY_EX*(HT_EX)**0.5D0 ) / ((HT_IN)**0.5D0 + (HT_EX)**0.5D0)

      q_Roe = q_RoeX*NX + q_RoeY*NY

      U_AVG = 0.5D0 * ( U_IN + U_EX )
      V_AVG = 0.5D0 * ( V_IN + V_EX )
      VEL_NORMAL = U_AVG*NX + V_AVG*NY

      do l=1,layers

         if (vel_normal.ge.0) THEN   !CAUTION, only makes sense for single layer
            bed_HAT(l) = discharge_modelX_IN*NX + discharge_modelY_IN*NY
         else
            bed_HAT(l) = discharge_modelX_EX*NX + discharge_modelY_EX*NY
         endif
         
      enddo
      
      
#endif

C.....Compute chem and tracer upwinding

#ifdef TRACE
      U_AVG = 0.5D0 * ( U_IN + U_EX )
      V_AVG = 0.5D0 * ( V_IN + V_EX )
      VEL_NORMAL = U_AVG*NX + V_AVG*NY
      
      if (VEL_NORMAL.GT.0) THEN
         I_HAT = VEL_NORMAL*iota_IN
      else
         I_HAT = VEL_NORMAL*iota_EX
      endif
#endif

#ifdef CHEM
      U_AVG = 0.5D0 * ( U_IN + U_EX )
      V_AVG = 0.5D0 * ( V_IN + V_EX )
      VEL_NORMAL = U_AVG*NX + V_AVG*NY
      
      if (VEL_NORMAL.GT.0) THEN
         I_HAT = VEL_NORMAL*iota_IN
         J_HAT = VEL_NORMAL*iota2_IN
      else
         I_HAT = VEL_NORMAL*iota_EX
         J_HAT = VEL_NORMAL*iota2_EX
      endif
#endif

#ifdef DYNP
      U_AVG = 0.5D0 * ( U_IN + U_EX )
      V_AVG = 0.5D0 * ( V_IN + V_EX )
      VEL_NORMAL = U_AVG*NX + V_AVG*NY
      
      if (VEL_NORMAL.GT.0) THEN
         K_HAT = VEL_NORMAL*dynP_IN
      else
         K_HAT = VEL_NORMAL*dynP_EX
      endif
#endif

      RETURN
      END SUBROUTINE
      
C***********************************************************************
C     
C     SUBROUTINE HLLC_FLUX()
C     
C     This subroutine computes the Harten-Lax-van Leer Contact (HLLC)
C     flux for the shallow water equations.
C     
C     Written by Ethan Kubatko
C     
C***********************************************************************

      SUBROUTINE HLLC_FLUX()

C.....Use appropriate modules

      USE GLOBAL
      USE DG
      use sizes, only: layers

      IMPLICIT NONE

C.....Declare local variables

      INTEGER II,l
      REAL(SZ) Q_STAR(3), Q(3)
      REAL(SZ) HUU, HUV, HVV, GH2, U_AVG, V_AVG, VEL_NORMAL
      REAL(SZ) C_EX, C_IN, UN_EX, UN_IN, UT_EX, UT_IN
      REAL(SZ) LS_EX, LS_IN, S_EX, S_IN, S_STAR
      REAL(SZ) H_STAR
      REAL(SZ) DEPTH, F1_NL, FU_NL, FV_NL, FG_NL, FH_NL
      REAL(SZ) FX1_IN, FX2_IN, FX3_IN, FY1_IN, FY2_IN, FY3_IN
      REAL(SZ) FX1_EX, FX2_EX, FX3_EX, FY1_EX, FY2_EX, FY3_EX
      Real(SZ) q_RoeX, q_RoeY, q_Roe
      Real(SZ) discharge_modelX_IN,discharge_modelX_EX
      Real(SZ) discharge_modelY_IN,discharge_modelY_EX
      Real(SZ) rVncU,rVncV,lambda1_IN,lambda1_EX,lambda2_IN,lambda2_EX
      Real(SZ) lambda3_IN,lambda3_EX,Select_IN,Select_EX

C.....Compute the interior and exterior water column heights

      HT_IN = ZE_IN*NLEQ + HB_IN 
      HT_EX = ZE_EX*NLEQ + HB_EX

C.....Compute the interior and exterior normal velocities

      FH_NL = 1.D0/(NLEQ*HT_IN + LEQ)
      U_IN  = QX_IN*FH_NL
      V_IN  = QY_IN*FH_NL
      UN_IN =  U_IN*NX + V_IN*NY
      UT_IN = -U_IN*NY + V_IN*NX

      FH_NL = 1.D0/(NLEQ*HT_EX + LEQ)
      U_EX  = QX_EX*FH_NL
      V_EX  = QY_EX*FH_NL
      UN_EX =  U_EX*NX + V_EX*NY
      UT_EX = -U_EX*NY + V_EX*NX

C.....Compute H*

      C_EX = SQRT(G*HT_EX)
      C_IN = SQRT(G*HT_IN)
      H_STAR = 0.5D0*(C_IN + C_EX) + 0.25D0*(UN_IN - UN_EX)
      H_STAR = 1.D0/G*H_STAR*H_STAR

C.....Compute the wave speeds

      LS_EX = 1.D0
      LS_IN = 1.D0
      IF (H_STAR.GT.HT_EX) LS_EX = SQRT( (H_STAR*H_STAR + H_STAR*HT_EX)
     &     /(2.D0*HT_EX*HT_EX) )
      IF (H_STAR.GT.HT_IN) LS_IN = SQRT( (H_STAR*H_STAR + H_STAR*HT_IN)
     &     /(2.D0*HT_IN*HT_IN) )

      S_IN   = UN_IN - C_IN*LS_IN
      S_EX   = UN_EX + C_EX*LS_EX
      S_STAR = ( S_IN*HT_EX*(UN_EX - S_EX) - S_EX*HT_IN*(UN_IN - S_IN) )
     &     /( HT_EX*(UN_EX - S_EX) - HT_IN*(UN_IN - S_IN) )

C.....Compute the numerical flux based on the wave speeds

      IF (S_IN.GE.0) THEN

C.......Compute continuity fluxes at exterior state

         F1_NL  = NLEQ + LEQ*HB_IN
         FX1_IN = QX_IN*F1_NL
         FY1_IN = QY_IN*F1_NL

C.......Compute momentum flux terms at exterior state

         FU_NL = NLEQ*QX_IN
         FV_NL = NLEQ*QY_IN
         FG_NL = NLEQG*ZE_IN

         HUU = FU_NL*U_IN
         HVV = FV_NL*V_IN
         HUV = FU_NL*V_IN
         GH2 = FG_NL*(0.5D0*ZE_IN + HB_IN) + FG_L*ZE_IN

C.......Compute x momentum fluxes at exterior state

         FX2_IN = HUU + GH2
         FY2_IN = HUV

C.......Compute y momentum fluxes at exterior state

         FX3_IN = HUV
         FY3_IN = HVV + GH2

C.......Set the numerical fluxes equal to the exterior fluxes

         F_HAT = FX1_IN*NX + FY1_IN*NY
         G_HAT = FX2_IN*NX + FY2_IN*NY
         H_HAT = FX3_IN*NX + FY3_IN*NY

      ELSEIF ((S_IN.LE.0).AND.(S_STAR.GE.0)) THEN

C.......Compute continuity fluxes at exterior state

         F1_NL  = NLEQ + LEQ*HB_IN
         FX1_IN = QX_IN*F1_NL
         FY1_IN = QY_IN*F1_NL

C.......Compute momentum flux terms at exterior state

         FU_NL = NLEQ*QX_IN
         FV_NL = NLEQ*QY_IN
         FG_NL = NLEQG*ZE_IN

         HUU = FU_NL*U_IN
         HVV = FV_NL*V_IN
         HUV = FU_NL*V_IN
         GH2 = FG_NL*(0.5D0*ZE_IN + HB_IN) + FG_L*ZE_IN

C.......Compute x momentum fluxes at exterior state

         FX2_IN = HUU + GH2
         FY2_IN = HUV

C.......Compute y momentum fluxes at exterior state

         FX3_IN = HUV
         FY3_IN = HVV + GH2

C.......Compute Q*

         Q_STAR(1) = HT_IN*((S_IN - UN_IN)/(S_IN - S_STAR))
         Q_STAR(2) = Q_STAR(1)*S_STAR
         Q_STAR(3) = Q_STAR(1)*UT_IN
         Q(1) =  HT_IN
         Q(2) =  QX_IN*NX + QY_IN*NY
         Q(3) = -QX_IN*NY + QY_IN*NX

C.......Compute the numerical flux

         F_HAT = FX1_IN*NX + FY1_IN*NY + S_IN*(  Q_STAR(1) - Q(1))
         G_HAT = FX2_IN*NX + FY2_IN*NY + S_IN*( (Q_STAR(2) - Q(2))*NX -
     &        (Q_STAR(3) - Q(3))*NY )
         H_HAT = FX3_IN*NX + FY3_IN*NY + S_IN*( (Q_STAR(2) - Q(2))*NY +
     &        (Q_STAR(3) - Q(3))*NX )

      ELSEIF ((S_STAR.LE.0).AND.(S_EX.GE.0)) THEN

C.......Compute continuity fluxes at exterior state

         F1_NL  = NLEQ + LEQ*HB_EX
         FX1_EX = QX_EX*F1_NL
         FY1_EX = QY_EX*F1_NL

C.......Compute momentum flux terms at exterior state

         FU_NL = NLEQ*QX_EX
         FV_NL = NLEQ*QY_EX
         FG_NL = NLEQG*ZE_EX

         HUU = FU_NL*U_EX
         HVV = FV_NL*V_EX
         HUV = FU_NL*V_EX
         GH2 = FG_NL*(0.5D0*ZE_EX + HB_EX) + FG_L*ZE_EX

C.......Compute x momentum fluxes at exterior state

         FX2_EX = HUU + GH2
         FY2_EX = HUV

C.......Compute y momentum fluxes at exterior state

         FX3_EX = HUV
         FY3_EX = HVV + GH2

C.......Compute Q*

         Q_STAR(1) = HT_EX*((S_EX - UN_EX)/(S_EX - S_STAR))
         Q_STAR(2) = Q_STAR(1)*S_STAR
         Q_STAR(3) = Q_STAR(1)*UT_EX
         Q(1) =  HT_EX
         Q(2) =  QX_EX*NX + QY_EX*NY
         Q(3) = -QX_EX*NY + QY_EX*NX

C.......Compute the numerical flux

         F_HAT = FX1_EX*NX + FY1_EX*NY + S_EX*(  Q_STAR(1) - Q(1))
         G_HAT = FX2_EX*NX + FY2_EX*NY + S_EX*( (Q_STAR(2) - Q(2))*NX -
     &        (Q_STAR(3) - Q(3))*NY )
         H_HAT = FX3_EX*NX + FY3_EX*NY + S_EX*( (Q_STAR(2) - Q(2))*NY +
     &        (Q_STAR(3) - Q(3))*NX )

      ELSEIF (S_EX.LE.0) THEN

C.......Compute continuity fluxes at exterior state

         F1_NL  = NLEQ + LEQ*HB_EX
         FX1_EX = QX_EX*F1_NL
         FY1_EX = QY_EX*F1_NL

C.......Compute momentum flux terms at exterior state

         FU_NL = NLEQ*QX_EX
         FV_NL = NLEQ*QY_EX
         FG_NL = NLEQG*ZE_EX

         HUU = FU_NL*U_EX
         HVV = FV_NL*V_EX
         HUV = FU_NL*V_EX
         GH2 = FG_NL*(0.5D0*ZE_EX + HB_EX) + FG_L*ZE_EX

C.......Compute x momentum fluxes at exterior state

         FX2_EX = HUU + GH2
         FY2_EX = HUV

C.......Compute y momentum fluxes at exterior state

         FX3_EX = HUV
         FY3_EX = HVV + GH2

C.......Set the numerical fluxes equal to the exterior fluxes

         F_HAT = FX1_EX*NX + FY1_EX*NY
         G_HAT = FX2_EX*NX + FY2_EX*NY
         H_HAT = FX3_EX*NX + FY3_EX*NY

      ENDIF

#ifdef SED_LAY !Not thoroughly tested!!

!.....Compute the stabilization term (v_nc in CMM's thesis) 
!.....for the nonconservative product in the momentum eqn

                                !Compute the magnitudes of the velocities for sediment

      QMag_IN = (QX_IN*QX_IN/(HT_IN**2) + QY_IN*QY_IN/(HT_IN**2) )**(1/2)
      QMag_EX = (QX_EX*QX_EX/(HT_EX**2) + QY_EX*QY_EX/(HT_EX**2) )**(1/2)

      rVncU = 0.D0
      rVncV = 0.D0

      do l=1,layers

         rVncU =  rVncU + 0.5D0*G*(HT_IN+HT_EX) 
     &        *(bed_EX(l)-bed_IN(l))*NX
         rVncV =  rVncV + 0.5D0*G*(HT_IN+HT_EX) 
     &        *(bed_EX(l)-bed_IN(l))*NY

      enddo

!.....When not using the NCP formalism
!.....then we enforce the (truncated) HLL wave-speed rendition.
!.....See CMM's thesis, pg 45, substituting for the base flux

      lambda1_IN = U_IN*NX + V_IN*NY
      lambda1_EX = U_EX*NX + V_EX*NY

      lambda2_IN = lambda1_IN - (g*HT_IN)**(1/2)
      lambda2_EX = lambda1_EX - (g*HT_EX)**(1/2)

      lambda3_IN = lambda1_IN + (g*HT_IN)**(1/2)
      lambda3_EX = lambda1_EX + (g*HT_EX)**(1/2)

      Select_IN = min(lambda2_IN,lambda2_EX)
      Select_EX = max(lambda3_IN,lambda3_EX)

      if(Select_IN.ge.0.D0) then 
         
         !F_HAT = FX1_IN*NX + FY1_IN*NY No update
         G_HAT = G_HAT  - 0.5D0*rVncU
         H_HAT = H_HAT  - 0.5D0*rVncV
         
      elseif(Select_IN.lt.0.D0.and.Select_EX.gt.0.D0) then

         !F_HAT = FX1_EX*NX + FY1_EX*NY  No update

         G_HAT = G_HAT - 0.5D0*(Select_IN + Select_EX)/
     &        (Select_EX-Select_IN)*rVncU
         H_HAT = H_HAT - 0.5D0*(Select_IN + Select_EX)/
     &        (Select_EX-Select_IN)*rVncV
         
      else

         !F_HAT = FX1_EX*NX + FY1_EX*NY 
         G_HAT = G_HAT  + 0.5D0*rVncU
         H_HAT = H_HAT  + 0.5D0*rVncV

      endif

!.....Computes the sediment discharge rate based on the formulation
!.....developed by Camenen and Larson, so $\tilde{q} = A_g*H^{-3}|q|^{2}q$
!.....which clearly makes sense only for a single layer

      discharge_modelX_IN =  porosity*HT_IN**(-1)*QMag_IN**(2)*QX_IN*SFAC_IN
      discharge_modelX_EX =  porosity*HT_EX**(-1)*QMag_EX**(2)*QX_EX*SFAC_EX
      discharge_modelY_IN =  porosity*HT_IN**(-1)*QMag_IN**(2)*QY_IN 
      discharge_modelY_EX =  porosity*HT_EX**(-1)*QMag_EX**(2)*QY_EX


      q_RoeX = ( QX_IN*(HT_IN)**0.5D0  
     &     + QX_EX*(HT_EX)**0.5D0 ) / ((HT_IN)**0.5D0 + (HT_EX)**0.5D0)
      q_RoeY = ( QY_IN*(HT_IN)**0.5D0  
     &     + QY_EX*(HT_EX)**0.5D0 ) / ((HT_IN)**0.5D0 + (HT_EX)**0.5D0)

      q_Roe = q_RoeX*NX + q_RoeY*NY

      do l=1,layers

         if (q_Roe.ge.0) THEN   !CAUTION, only makes sense for single layer
            bed_HAT(l) = discharge_modelX_IN*NX + discharge_modelY_IN*NY
         else
            bed_HAT(l) = discharge_modelX_EX*NX + discharge_modelY_EX*NY
         endif
         
      enddo
      
      
#endif

#ifdef TRACE
      U_AVG = 0.5D0 * ( U_IN + U_EX )
      V_AVG = 0.5D0 * ( V_IN + V_EX )
      VEL_NORMAL = U_AVG*NX + V_AVG*NY
      
      if (VEL_NORMAL.GT.0) THEN
         I_HAT = VEL_NORMAL*iota_IN
      else
         I_HAT = VEL_NORMAL*iota_EX
      endif
#endif

#ifdef CHEM
      U_AVG = 0.5D0 * ( U_IN + U_EX )
      V_AVG = 0.5D0 * ( V_IN + V_EX )
      VEL_NORMAL = U_AVG*NX + V_AVG*NY
      
      if (VEL_NORMAL.GT.0) THEN
         I_HAT = VEL_NORMAL*iota_IN
         J_HAT = VEL_NORMAL*iota2_IN
      else
         I_HAT = VEL_NORMAL*iota_EX
         J_HAT = VEL_NORMAL*iota2_EX
      endif
#endif

#ifdef DYNP
      U_AVG = 0.5D0 * ( U_IN + U_EX )
      V_AVG = 0.5D0 * ( V_IN + V_EX )
      VEL_NORMAL = U_AVG*NX + V_AVG*NY
      
      if (VEL_NORMAL.GT.0) THEN
         K_HAT = VEL_NORMAL*dynP_IN
      else
         K_HAT = VEL_NORMAL*dynP_EX
      endif
#endif

      
      Return
      END SUBROUTINE


C***********************************************************************
C     
C     subroutine NCP_flux()
C     
C     This subroutine uses CMM's version of the NCP flux.  
C     See his thesis pages 2.4.3
C     
C     NOTE: THIS FLUX IS ONLY MEANT FOR WHEN SEDIMENT IS ON!!
C     because the formalism is more tightly coupled
C     
C     
C     Written by cem interpreted from CM's code revamp
C     
C***********************************************************************

      SUBROUTINE NCP_FLUX()

C.....Use appropriate modules

      USE GLOBAL
      USE DG
      use sizes, only: layers

      IMPLICIT NONE

C.....Declare local variables

      INTEGER II,ll
      REAL(SZ) Q_STAR(3), Q(3)
      REAL(SZ) HUU, HUV, HVV, GH2, U_AVG, V_AVG, VEL_NORMAL
      REAL(SZ) rB,rC,rD,rQ,rR,rTheta,rAlph,kappa2,GH3
      REAL(SZ) DEPTH, F1_NL, FU_NL, FV_NL, FG_NL, FH_NL_IN,FH_NL_EX
      REAL(SZ) FX1_IN, FX2_IN, FX3_IN, FY1_IN, FY2_IN, FY3_IN
      REAL(SZ) FX1_EX, FX2_EX, FX3_EX, FY1_EX, FY2_EX, FY3_EX
      Real(SZ) q_RoeX, q_RoeY, q_Roe
      Real(SZ) discharge_modelX_IN,discharge_modelX_EX
      Real(SZ) discharge_modelY_IN,discharge_modelY_EX
      Real(SZ) rVncU,rVncV,lambda1_IN,lambda1_EX,lambda2_IN,lambda2_EX
      Real(SZ) lambda3_IN,lambda3_EX,lambda4_IN,lambda4_EX
      Real(SZ) Select_IN,Select_EX

! Beware of well-balancedness issues in this formulation

C.....Compute the interior and exterior water column heights

      HT_IN = ZE_IN*IFNLFA + HB_IN
      HT_EX = ZE_EX*IFNLFA + HB_EX

C.....Compute the interior and exterior normal velocities

      FH_NL_IN = 1.D0/(NLEQ*HT_IN + LEQ)
      U_IN  = QX_IN*FH_NL_IN
      V_IN  = QY_IN*FH_NL_IN
 
      FH_NL_EX = 1.D0/(NLEQ*HT_EX + LEQ)
      U_EX  = QX_EX*FH_NL_EX
      V_EX  = QY_EX*FH_NL_EX
 
C.....Compute the interior and exterior discharge rates

      QMag_IN = (QX_IN*QX_IN/(HT_IN**2) + QY_IN*QY_IN/(HT_IN**2) )**(1/2)
      QMag_EX = (QX_EX*QX_EX/(HT_EX**2) + QY_EX*QY_EX/(HT_EX**2) )**(1/2)

      discharge_modelX_IN =  porosity*HT_IN**(-1)*QMag_IN**(2)*(QX_IN) *SFAC_IN
      discharge_modelX_EX =  porosity*HT_EX**(-1)*QMag_EX**(2)*(QX_EX) *SFAC_EX
      discharge_modelY_IN =  porosity*HT_IN**(-1)*QMag_IN**(2)*(QX_IN) 
      discharge_modelY_EX =  porosity*HT_EX**(-1)*QMag_EX**(2)*(QY_EX)

C.....Compute the eigenvalues and aux variables of interior states

                                !Army of definitions

      lambda1_IN = U_IN*NX + V_IN*NY
 
      rAlph = porosity ! constant
      kappa2 = 1.D0 !kapp2 = 1.D0 - porosity

      rB = -2.D0*lambda1_IN !hardwired for \tilde{q} independent of b!

      rC = (1.D0 - 2.D0*rAlph*G/(kappa2) )* lambda1_IN**2 
     & - G*HT_IN*(HT_IN**(-3)*( QMag_IN * rAlph/kappa2) + 1.D0)


      rD = rAlph*G*lambda1_IN/kappa2  
     & * (3.D0*HT_IN**(-2)*QMag_IN - 2.D0*lambda1_IN**2)


      rQ = (3.D0*rC - rB**2) / 9.D0

      rR = (9.D0*rB*rC - 27.D0*rD - 2.D0*rB**3) / 54.D0

      rTheta = acos(rR / sqrt(-rQ**3))

      lambda2_IN = 2.D0*sqrt(-rQ)*cos(rTheta/3.D0)+ (2.D0/3.D0)*lambda1_IN

      lambda3_IN = 2.D0*sqrt(-rQ)*cos(rTheta/3.D0+ 2.D0*PI/3.D0)+ (2.D0/3.D0)*lambda1_IN

      lambda4_IN = 2.D0*sqrt(-rQ)*cos(rTheta/3.D0+ 4.D0*PI/3.D0)+ (2.D0/3.D0)*lambda1_IN

C.....Compute the eigenvalues and aux vairables of exterior states

                                !Army of definitions

      lambda1_EX = U_EX*NX + V_EX*NY

      rB = -2.D0*lambda1_EX !hardwired for \tilde{q} independent of b!

      rC = (1.D0 - 2.D0*rAlph*G/(kappa2) )* lambda1_EX**2 
     & - G*HT_EX*(rAlph*HT_EX**(-3)*( QMag_EX * rAlph/kappa2) + 1.D0)

      rD = rAlph*G*lambda1_EX/kappa2  
     & * (3.D0*HT_EX**(-2)*QMag_EX - 2.D0*lambda1_EX**2)

      rQ = (3.D0*rC - rB**2) / 9.D0

      rR = (9.D0*rB*rC - 27.D0*rD - 2.D0*rB**3) / 54.D0

      rTheta = acos(rR / sqrt(-rQ**3))

      lambda2_EX = 2.D0*sqrt(-rQ)*cos(rTheta/3.D0)+ (2.D0/3.D0)*lambda1_EX

      lambda3_EX = 2.D0*sqrt(-rQ)*cos(rTheta/3.D0+ 2.D0*PI/3.D0)+ (2.D0/3.D0)*lambda1_EX

      lambda4_EX = 2.D0*sqrt(-rQ)*cos(rTheta/3.D0+ 4.D0*PI/3.D0)+ (2.D0/3.D0)*lambda1_EX

C.....Select the correction displacement stuff for the momentum

      Select_IN = min(min(lambda1_IN,lambda2_IN,lambda3_IN,lambda4_IN),
     &                min(lambda1_EX,lambda2_EX,lambda3_EX,lambda4_EX))
      Select_EX = max(max(lambda1_IN,lambda2_IN,lambda3_IN,lambda4_IN),
     &                max(lambda1_EX,lambda2_EX,lambda3_EX,lambda4_EX))
      rVncU = 0.D0
      rVncV = 0.D0

      do ll=1,layers !should be fixed for multiple layers

         rVncU =  rVncU + 0.5D0*G*(HT_IN+HT_EX) 
     &        *(bed_EX(ll)-bed_IN(ll))*NX
         rVncV =  rVncV + 0.5D0*G*(HT_IN+HT_EX) 
     &        *(bed_EX(ll)-bed_IN(ll))*NY

      enddo

!.....First Case

      if (Select_IN.ge.0) then

C.......Compute continuity fluxes at interior state

         F1_NL  = NLEQ + LEQ*HB_IN
         FX1_IN = QX_IN*F1_NL*SFAC_IN
         FY1_IN = QY_IN*F1_NL

C.......Compute momentum flux terms at interior state

         FU_NL = NLEQ*QX_IN
         FV_NL = NLEQ*QY_IN
         FG_NL = NLEQG*ZE_IN

         HUU = FU_NL*U_IN
         HVV = FV_NL*V_IN
         HUV = FU_NL*V_IN
         GH2 = FG_NL*(0.5D0*ZE_IN + HB_IN) + FG_L*ZE_IN
         GH3 = 0.5D0*G*(HT_IN**2)

C.......Compute x momentum fluxes at interior state

         FX2_IN = (HUU + GH3)*SFAC_IN
         FY2_IN = HUV

C.......Compute y momentum fluxes at interior state

         FX3_IN = HUV*SFAC_IN
         FY3_IN = HVV + GH3

!.......Set the fluxes
         
         F_HAT = FX1_IN*NX + FY1_IN*NY 
         G_HAT = FX2_IN*NX + FY2_IN*NY - 0.5D0*rVncU
         H_HAT = FX3_IN*NX + FY3_IN*NY - 0.5D0*rVncV

         do ll = 1,layers        !CAUTION, only really makes sense for single layer

               bed_HAT(ll) = 0.5D0*(discharge_modelX_IN+discharge_modelX_EX)*NX 
     &           + 0.5D0*(discharge_modelY_IN+discharge_modelY_EX)*NY

         enddo

 
!...Second case
      elseif(Select_IN.lt.0.D0.and.Select_EX.ge.0.D0) then

C.......Compute continuity fluxes

         F1_NL  = NLEQ + LEQ*HB_IN
         FX1_IN = QX_IN*F1_NL*SFAC_IN
         FY1_IN = QY_IN*F1_NL
         F1_NL  = NLEQ + LEQ*HB_EX
         FX1_EX = QX_EX*F1_NL*SFAC_EX
         FY1_EX = QY_EX*F1_NL

C.......Compute momentum fluxes interior

         FU_NL = NLEQ*QX_IN
         FV_NL = NLEQ*QY_IN
         FG_NL = NLEQG*ZE_IN

         HUU = FU_NL*U_IN
         HVV = FV_NL*V_IN
         HUV = FU_NL*V_IN
         GH2 = FG_NL*(0.5D0*ZE_IN + HB_IN) + FG_L*ZE_IN
         GH3 = 0.5D0*G*(HT_IN**2)

C.......Compute x momentum fluxes interior

         FX2_IN = (HUU + GH3)*SFAC_IN
         FY2_IN = HUV

C.......Compute y momentum fluxes interior

         FX3_IN = HUV*SFAC_IN
         FY3_IN = HVV + GH3

C.......Compute momentum fluxes exterior 

         FU_NL = NLEQ*QX_EX
         FV_NL = NLEQ*QY_EX
         FG_NL = NLEQG*ZE_EX

         HUU = FU_NL*U_EX
         HVV = FV_NL*V_EX
         HUV = FU_NL*V_EX
         GH2 = FG_NL*(0.5D0*ZE_EX + HB_EX) + FG_L*ZE_EX
         GH3 = 0.5D0*G*(HT_EX**2)

C.......Compute x momentum fluxes exterior

         FX2_EX = (HUU + GH3)*SFAC_EX
         FY2_EX = HUV

C.......Compute y momentum fluxes at exterior state

         FX3_EX = HUV*SFAC_EX
         FY3_EX = HVV + GH3


!.......Assign NCP/HLL version of fluxes
         
         F_HAT = ((Select_EX*(FX1_IN*NX + FY1_IN*NY) - Select_IN*(FX1_EX*NX + FY1_EX*NY))
     &        - (Select_IN*Select_EX)*(HT_IN-HT_EX) ) / (Select_EX-Select_IN)
     &        
         G_HAT = ((Select_EX*(FX2_IN*NX + FY2_IN*NY) - Select_IN*(FX2_EX*NX + FY2_EX*NY))
     &        - (Select_IN*Select_EX)*(QX_IN-QX_EX) - 0.5D0*(Select_IN+Select_EX)*rVncU)
     &        / (Select_EX-Select_IN)
  
         H_HAT =  ((Select_EX*(FX3_IN*NX + FY3_IN*NY) - Select_IN*(FX3_EX*NX + FY3_EX*NY))
     &        - (Select_IN*Select_EX)*(QY_IN-QY_EX) - 0.5D0*(Select_IN+Select_EX)*rVncV )
     &        / (Select_EX-Select_IN)


         do ll = 1,layers 

            bed_HAT(ll) = ((Select_EX*(discharge_modelX_IN*NX + discharge_modelY_IN*NY)
     &           - Select_IN*(discharge_modelX_EX*NX + discharge_modelY_EX*NY))) 
     &           - ((Select_IN*Select_EX) * (bed_IN(ll) - bed_EX(ll)))
     &           / (Select_EX-Select_IN) 
         enddo


!...Third case 
      else
   

C.......Compute continuity fluxes at exterior state

         F1_NL  = NLEQ + LEQ*HB_EX
         FX1_EX = QX_EX*F1_NL*SFAC_EX
         FY1_EX = QY_EX*F1_NL

C.......Compute momentum flux terms at exterior state

         FU_NL = NLEQ*QX_EX
         FV_NL = NLEQ*QY_EX
         FG_NL = NLEQG*ZE_EX

         HUU = FU_NL*U_EX
         HVV = FV_NL*V_EX
         HUV = FU_NL*V_EX
         GH2 = FG_NL*(0.5D0*ZE_EX + HB_EX) + FG_L*ZE_EX
         GH3 = 0.5D0*G*(HT_EX**2)

C.......Compute x momentum fluxes at exterior state

         FX2_EX = (HUU + GH3)*SFAC_EX
         FY2_EX = HUV

C.......Compute y momentum fluxes at exterior state

         FX3_EX = HUV
         FY3_EX = HVV + GH3

C.......Compute y momentum fluxes at exterior state

         FX3_EX = HUV*SFAC_EX
         FY3_EX = HVV + GH3

!.......Set the fluxes
         
         F_HAT = FX1_EX*NX + FY1_EX*NY 
         G_HAT = FX2_EX*NX + FY2_EX*NY + 0.5D0*rVncU
         H_HAT = FX3_EX*NX + FY3_EX*NY + 0.5D0*rVncV

         do ll = 1,layers        !CAUTION, only really makes sense for single layer

               bed_HAT(ll) = 0.5D0*(discharge_modelX_IN+discharge_modelX_EX)*NX 
     &           + 0.5D0*(discharge_modelY_IN+discharge_modelY_EX)*NY

         enddo

      endif
      
#ifdef TRACE
      U_AVG = 0.5D0 * ( U_IN + U_EX )
      V_AVG = 0.5D0 * ( V_IN + V_EX )
      VEL_NORMAL = U_AVG*NX + V_AVG*NY
      
      if (VEL_NORMAL.GT.0) THEN
         I_HAT = VEL_NORMAL*iota_IN
      else
         I_HAT = VEL_NORMAL*iota_EX
      endif
#endif

#ifdef CHEM
      U_AVG = 0.5D0 * ( U_IN + U_EX )
      V_AVG = 0.5D0 * ( V_IN + V_EX )
      VEL_NORMAL = U_AVG*NX + V_AVG*NY
      
      if (VEL_NORMAL.GT.0) THEN
         I_HAT = VEL_NORMAL*iota_IN
         J_HAT = VEL_NORMAL*iota2_IN
      else
         I_HAT = VEL_NORMAL*iota_EX
         J_HAT = VEL_NORMAL*iota2_EX
      endif
#endif

#ifdef DYNP
      U_AVG = 0.5D0 * ( U_IN + U_EX )
      V_AVG = 0.5D0 * ( V_IN + V_EX )
      VEL_NORMAL = U_AVG*NX + V_AVG*NY
      
      if (VEL_NORMAL.GT.0) THEN
         K_HAT = VEL_NORMAL*dynP_IN
      else
         K_HAT = VEL_NORMAL*dynP_EX
      endif
#endif


      RETURN
      END SUBROUTINE
      
      
