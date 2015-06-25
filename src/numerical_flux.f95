!***********************************************************************
!     
!     SUBROUTINE NUMERICAL_FLUX
!     
!     This subroutine calls an appropriate numrical flux subroutine
!     based on FLUXTYPE variable.
!     
!     Written by Shintaro Bunya (11-08-2006)
!
!     01-10-2011 - cem - adapted for p_enrichment and multicomponent
!     06-01-2012 - cem - New fluxes, sediment, etc.
!     
!***********************************************************************

      SUBROUTINE NUMERICAL_FLUX(s,dg_here,IT,MM)
      
      USE DG
      use sizes

      IMPLICIT NONE

      type (sizes_type) :: s
      type (dg_type) :: dg_here

      INTEGER IT,MM


      IF(dg_here%FLUXTYPE.EQ.1) THEN
         CALL ROE_FLUX(s,dg_here,IT,MM)
      ELSEIF(dg_here%FLUXTYPE.EQ.2) THEN
         CALL LLF_FLUX(s,dg_here)
      ELSEIF (dg_here%FLUXTYPE.EQ.3) THEN
         CALL HLLC_FLUX(s,dg_here)
      ELSEIF (dg_here%FLUXTYPE.EQ.4) THEN
         CALL NCP_FLUX(s,dg_here)
      ELSE
         STOP 'INVALID FLUXTYPE'
      ENDIF

      RETURN
      END SUBROUTINE
 

!***********************************************************************
!     
!     SUBROUTINE ROE_FLUX
!     
!     This subroutine computes the Roe flux for the shallow water
!     equations.
!     
!     Written by Ethan Kubatko (06-11-2004)
!     
!***********************************************************************
      SUBROUTINE ROE_FLUX(s,dg_here,IT,MM)
      
!.....Use appropriate modules

      USE GLOBAL
      USE DG
      use sizes
      use fparser
      use fparser2

      IMPLICIT NONE

      type (sizes_type) :: s
      type (dg_type) :: dg_here
      
!.....Declare local variables
      
      INTEGER II,GED,ll,w,IT,kk,jj,MM,fast
      REAL(SZ) DEN, U_AVG, V_AVG, VEL_NORMAL, q_RoeX, q_RoeY, q_Roe
      Real(SZ) bedX_IN,bedX_EX,dtildeQx(4),dtildeQy(4)
      Real(SZ) bedY_IN,bedY_EX, bed_AVG(s%layers),lambda4_IN,lambda4_EX
      Real(SZ) rVncU,rVncV,lambda1_IN,lambda1_EX,lambda2_IN,lambda2_EX
      Real(SZ) lambda3_IN,lambda3_EX,Select_IN,Select_EX
      Real T1, T2, T3, T4
      Real(SZ) chi,alpha_0,beta_0,gamma1,gamma2,gamma3,gamma4
      Real(SZ) iot,iot2,iot3,iot4,iot5,J6,J7,J8,bigA,bigB
      Real(SZ) bigD,bigE,bigF,J9,epsilon,test_g,xi
      
!.....Compute the total height of the water column

      HT_IN = dg_here%ZE_IN*IFNLFA + dg_here%HB_IN
      HT_EX = dg_here%ZE_EX*IFNLFA + dg_here%HB_EX

!.....Compute the velocities in the x and y directions

      U_IN = dg_here%QX_IN/HT_IN
      U_EX = dg_here%QX_EX/HT_EX

      V_IN = dg_here%QY_IN/HT_IN
      V_EX = dg_here%QY_EX/HT_EX

!.....Compute the Roe averaged variables

      T1    = SQRT(HT_IN)
      T2    = SQRT(HT_EX)
      T3    = SQRT(HT_IN + HT_EX)
      T4    = 1.0D0/(T1 + T2)
      C_ROE = G2ROOT*T3*SQRT(dg_here%NY**2 + (dg_here%NX*dg_here%SFAC_IN)**2) !srb - spherical coordinate correction
      U_ROE = T4*(U_IN*T1 + U_EX*T2)
      V_ROE = T4*(V_IN*T1 + V_EX*T2)

!.....Evaluate the eigenvalues at the Roe averaged variables

      EIGVAL(2) = (U_ROE*dg_here%NX + V_ROE*dg_here%NY)*IFNLCT
      EIGVAL(1) = EIGVAL(2) + C_ROE
      EIGVAL(3) = EIGVAL(2) - C_ROE

      !eigmax = max(abs(EIGVAL(1)),abs(EIGVAL(2)),abs(EIGVAL(3)))
      !cfl_max = max(EIGMAX,cfl_max)
      
!.....Evaluate right eigenvectors at Roe averaged variables

      RI(1,1) = 1.D0
      RI(2,1) = U_ROE*IFNLCT + (IFNLCT*C_ROE + (1-IFNLCT)*G/C_ROE)*dg_here%NX*dg_here%SFAC_IN !srb - spherical coordinate correction
      RI(3,1) = V_ROE*IFNLCT + (IFNLCT*C_ROE + (1-IFNLCT)*G/C_ROE)*dg_here%NY

      RI(1,2) = 0.D0
      RI(2,2) = -dg_here%NY
      RI(3,2) =  dg_here%NX*dg_here%SFAC_IN

      RI(1,3) = 1.D0
      RI(2,3) = U_ROE*IFNLCT - (IFNLCT*C_ROE + (1-IFNLCT)*G/C_ROE)*dg_here%NX*dg_here%SFAC_IN !srb - spherical coordinate correction
      RI(3,3) = V_ROE*IFNLCT - (IFNLCT*C_ROE + (1-IFNLCT)*G/C_ROE)*dg_here%NY

!.....Evaluate left eigenvectors at Roe averaged variables

      DEN = 1.0d0/(RI(2,2)*RI(3,3)-RI(2,3)*RI(3,2)+RI(2,1)*RI(3,2)&
     -RI(3,1)*RI(2,2))

      LE(1,1) =  DEN*( RI(2,2)*RI(3,3)-RI(2,3)*RI(3,2) )
      LE(1,2) =  DEN*RI(3,2)
      LE(1,3) = -DEN*RI(2,2)

      LE(2,1) = -DEN*( RI(2,1)*RI(3,3)-RI(2,3)*RI(3,1) )
      LE(2,2) =  DEN*( RI(3,3)-RI(3,1) )
      LE(2,3) = -DEN*( RI(2,3)-RI(2,1) )

      LE(3,1) =  DEN*( RI(2,1)*RI(3,2) - RI(2,2)*RI(3,1) )
      LE(3,2) = -DEN*RI(3,2)
      LE(3,3) =  DEN*RI(2,2)

!.....Compute the Roe matrix

      DO II = 1,3

         A_ROE(II,1) = RI(II,1)*ABS(EIGVAL(1))*LE(1,1) +&
                 RI(II,2)*ABS(EIGVAL(2))*LE(2,1) +&
                 RI(II,3)*ABS(EIGVAL(3))*LE(3,1)

         A_ROE(II,2) = RI(II,1)*ABS(EIGVAL(1))*LE(1,2) +&
                 RI(II,2)*ABS(EIGVAL(2))*LE(2,2) +&
                 RI(II,3)*ABS(EIGVAL(3))*LE(3,2)

         A_ROE(II,3) = RI(II,1)*ABS(EIGVAL(1))*LE(1,3) +&
                 RI(II,2)*ABS(EIGVAL(2))*LE(2,3) +&
                 RI(II,3)*ABS(EIGVAL(3))*LE(3,3)

      ENDDO

!.....Compute continuity fluxes at exterior state

      FX_EX = (IFNLCT*dg_here%QX_EX + (1-IFNLCT)*dg_here%QX_EX*HT_EX)*dg_here%SFAC_EX
      FY_EX = IFNLCT*dg_here%QY_EX + (1-IFNLCT)*dg_here%QY_EX*HT_EX

!.....Compute x momentum fluxes at exterior state

      GX_EX = IFNLCT*(dg_here%QX_EX*dg_here%QX_EX/(dg_here%HB_EX + dg_here%ZE_EX) +&
     0.5D0*G*(dg_here%ZE_EX*dg_here%ZE_EX + 2.D0*dg_here%ZE_EX*dg_here%HB_EX)) + (1-IFNLCT)*G*dg_here%ZE_EX
      GX_EX = GX_EX*dg_here%SFAC_EX
      GY_EX = IFNLCT*dg_here%QX_EX*dg_here%QY_EX/(dg_here%HB_EX + dg_here%ZE_EX)

!.....Compute y momentum fluxes at exterior state

      HX_EX = (IFNLCT*dg_here%QX_EX*dg_here%QY_EX/(dg_here%HB_EX + dg_here%ZE_EX))*dg_here%SFAC_EX
      HY_EX = IFNLCT*(dg_here%QY_EX*dg_here%QY_EX/(dg_here%HB_EX + dg_here%ZE_EX) +&
     0.5D0*G*(dg_here%ZE_EX*dg_here%ZE_EX + 2.D0*dg_here%ZE_EX*dg_here%HB_EX)) + (1-IFNLCT)*G*dg_here%ZE_EX

!.....Compute continuity fluxes at interior state

      FX_IN = (IFNLCT*dg_here%QX_IN + (1-IFNLCT)*dg_here%QX_IN*HT_IN)*dg_here%SFAC_IN
      FY_IN = IFNLCT*dg_here%QY_IN + (1-IFNLCT)*dg_here%QY_IN*HT_IN

!.....Compute x momentum fluxes at interior state

      GX_IN = IFNLCT*(dg_here%QX_IN*dg_here%QX_IN/(dg_here%HB_IN + dg_here%ZE_IN) +&
     0.5D0*G*(dg_here%ZE_IN*dg_here%ZE_IN + 2.D0*dg_here%ZE_IN*dg_here%HB_IN)) + (1-IFNLCT)*G*dg_here%ZE_IN

      GX_IN = GX_IN*dg_here%SFAC_IN
      GY_IN = IFNLCT*dg_here%QX_IN*dg_here%QY_IN/(dg_here%HB_IN + dg_here%ZE_IN)

!.....Compute y momentum fluxes at interior state

      HX_IN = (IFNLCT*dg_here%QX_IN*dg_here%QY_IN/(dg_here%HB_IN + dg_here%ZE_IN))*dg_here%SFAC_IN
      HY_IN = IFNLCT*(dg_here%QY_IN*dg_here%QY_IN/(dg_here%HB_IN + dg_here%ZE_IN) +&
     0.5D0*G*(dg_here%ZE_IN*dg_here%ZE_IN + 2.D0*dg_here%ZE_IN*dg_here%HB_IN)) + (1-IFNLCT)*G*dg_here%ZE_IN

!.....Compute the average flux in the normal direction

      F_AVG = 0.5D0*((FX_IN + FX_EX)*dg_here%NX + (FY_IN + FY_EX)*dg_here%NY)
      G_AVG = 0.5D0*((GX_IN + GX_EX)*dg_here%NX + (GY_IN + GY_EX)*dg_here%NY)
      H_AVG = 0.5D0*((HX_IN + HX_EX)*dg_here%NX + (HY_IN + HY_EX)*dg_here%NY)

!$$$      eigmax = max(abs(0.5D0*((FX_IN + FX_EX) + (FY_IN + FY_EX))),
!$$$     &     abs( 0.5D0*((GX_IN + GX_EX) + (GY_IN + GY_EX))),
!$$$     &     abs( 0.5D0*((HX_IN + HX_EX) + (HY_IN + HY_EX))) )
      
!$$$      cfl_max = max(EIGMAX,cfl_max)

!.....Compute the jump in the variables

      JUMP(1) = dg_here%ZE_EX - dg_here%ZE_IN
      JUMP(2) = dg_here%QX_EX - dg_here%QX_IN
      JUMP(3) = dg_here%QY_EX - dg_here%QY_IN

!.....Compute the Roe flux

      F_HAT = F_AVG - 0.5D0*(A_ROE(1,1)*JUMP(1) + A_ROE(1,2)*JUMP(2) +&
     A_ROE(1,3)*JUMP(3) )
      G_HAT = G_AVG - 0.5D0*(A_ROE(2,1)*JUMP(1) + A_ROE(2,2)*JUMP(2) +&
     A_ROE(2,3)*JUMP(3) )
      H_HAT = H_AVG - 0.5D0*(A_ROE(3,1)*JUMP(1) + A_ROE(3,2)*JUMP(2) +&
     A_ROE(3,3)*JUMP(3) )

      !cfl_max = max(abs(F_HAT),abs(G_HAT),abs(H_HAT),cfl_max)

#ifdef SED_LAY

!$$$!.....Not well-balanced, be careful here!
!$$$
!$$$!.....Let us use the fully coupled system, which leads to a different eigenproblem
!$$$
!$$$      GX_IN = 0.D0
!$$$      GY_IN = 0.D0
!$$$      HX_IN = 0.D0
!$$$      HY_IN = 0.D0
!$$$      GX_EX = 0.D0
!$$$      GY_EX = 0.D0
!$$$      HX_EX = 0.D0
!$$$      HY_EX = 0.D0
!$$$      FX_IN = 0.D0
!$$$      FY_IN = 0.D0
!$$$      FX_EX = 0.D0
!$$$      Fy_EX = 0.D0
!$$$
!$$$      f_Hat = 0.D0
!$$$      g_hat = 0.D0
!$$$      h_hat = 0.D0
!$$$      dg_here%bed_hat = 0.D0
!$$$
!$$$      test_g = 0.D0
!$$$
!$$$C.....Compute continuity fluxes at interior state
!$$$
!$$$      FX_IN = dg_here%QX_IN*dg_here%SFAC_IN
!$$$      FY_IN = dg_here%QY_IN
!$$$
!$$$C.....Compute continuity fluxes at exterior state
!$$$
!$$$      FX_EX = dg_here%QX_EX*dg_here%SFAC_EX
!$$$      FY_EX = dg_here%QY_EX
!$$$
!$$$C.....Compute x momentum fluxes at interior state
!$$$
!$$$      GX_IN = (dg_here%QX_IN*dg_here%QX_IN/(HT_IN) 
!$$$     &     + 0.5D0*G*(HT_IN**2))*dg_here%SFAC_IN 
!$$$      !&     - (HT_IN*dg_here%bed_IN(1)*G)*dg_here%SFAC_IN 
!$$$      GY_IN = dg_here%QX_IN*dg_here%QY_IN/(HT_IN)
!$$$
!$$$C.....Compute y momentum fluxes at interior state
!$$$
!$$$      HX_IN = (dg_here%QX_IN*dg_here%QY_IN/(HT_IN))*dg_here%SFAC_IN
!$$$      HY_IN = (dg_here%QY_IN*dg_here%QY_IN/(HT_IN) +
!$$$     &     0.5D0*G*(HT_IN**2))
!$$$      !&     - (HT_IN*dg_here%bed_IN(1)*G)
!$$$
!$$$C.....Compute x momentum fluxes at exterior state
!$$$
!$$$      GX_EX = (dg_here%QX_EX*dg_here%QX_EX/(HT_EX) 
!$$$     &     + 0.5D0*G*(HT_EX**2))*dg_here%SFAC_EX
!$$$      !&     - (HT_EX*dg_here%bed_EX(1)*G)*dg_here%SFAC_EX
!$$$      GY_EX = dg_here%QX_EX*dg_here%QY_EX/(HT_EX)
!$$$
!$$$C.....Compute y momentum fluxes at exterior state
!$$$
!$$$      HX_EX = (dg_here%QX_EX*dg_here%QY_EX/(HT_EX))*dg_here%SFAC_EX
!$$$      HY_EX = (dg_here%QY_EX*dg_here%QY_EX/(HT_EX) +
!$$$     &     0.5D0*G*(HT_EX**2))
!$$$      !&     - (HT_EX*dg_here%bed_EX(1)*G)
!$$$
!$$$
!$$$C.....Compute the average flux in the normal direction
!$$$
!$$$      F_AVG = 0.5D0*((FX_IN + FX_EX)*dg_here%NX + (FY_IN + FY_EX)*dg_here%NY)
!$$$      G_AVG = 0.5D0*((GX_IN + GX_EX)*dg_here%NX + (GY_IN + GY_EX)*dg_here%NY)
!$$$      H_AVG = 0.5D0*((HX_IN + HX_EX)*dg_here%NX + (HY_IN + HY_EX)*dg_here%NY)

      dg_here%QMag_IN = (dg_here%QX_IN*dg_here%QX_IN/(HT_IN**2) + dg_here%QY_IN*dg_here%QY_IN/(HT_IN**2) )**(1/2)
      dg_here%QMag_EX = (dg_here%QX_EX*dg_here%QX_EX/(HT_EX**2) + dg_here%QY_EX*dg_here%QY_EX/(HT_EX**2) )**(1/2)

      !Note that the choice of linearization can cause this to need updating....
      bedX_IN =  dg_here%porosity*HT_IN**(-1) * dg_here%QMag_IN**(2) *dg_here%QX_IN *dg_here%SFAC_IN
      bedX_EX =  dg_here%porosity*HT_EX**(-1) * dg_here%QMag_EX**(2) *dg_here%QX_EX *dg_here%SFAC_EX
      bedY_IN =  dg_here%porosity*HT_IN**(-1) * dg_here%QMag_IN**(2) *dg_here%QY_IN 
      bedY_EX =  dg_here%porosity*HT_EX**(-1) * dg_here%QMag_EX**(2) *dg_here%QY_EX

!.....Compute the remaining Roe averaged variables

      ZE_ROE  = (dg_here%ZE_IN + dg_here%ZE_EX) / 2.D0
      QX_ROE  = (dg_here%QX_IN + dg_here%QX_EX) / 2.D0
      QY_ROE  = (dg_here%QY_IN + dg_here%QY_EX) / 2.D0
      bed_ROE = (dg_here%bed_IN(1)+dg_here%bed_EX(1)) / 2.D0


      !if (it.eq.1.and.mm.eq.1) then !do this only once, because it is Slooooow ....
      if (.not.init_parser) then

         open(446, file = "./db_partials_X")
         read(446,'(a)') dg_here%funcx(1)
         read(446,'(a)') dg_here%funcx(2)
         read(446,'(a)') dg_here%funcx(3)
         read(446,'(a)') dg_here%funcx(4)
         close(446)
         open(447, file = "./db_partials_Y")
         read(447,'(a)') dg_here%funcy(1)
         read(447,'(a)') dg_here%funcy(2)
         read(447,'(a)') dg_here%funcy(3)
         read(447,'(a)') dg_here%funcy(4)
         close(447)

         !Need to tokenize the sed functions
         CALL initf (4)         ! Initialize function parser for nfunc functions
         CALL initf2 (4)        ! Initialize function parser for nfunc functions
         do w=1,4               ! Loop over functions
            CALL parsef (w, dg_here%funcx(w), varx) ! Parse and bytecompile w-th function string 
            CALL parsef2 (w, dg_here%funcy(w), vary) ! Parse and bytecompile w-th function string 
         enddo

         init_parser = .true.

      endif
      
      dg_here%valx(1) = ZE_ROE
      dg_here%valx(2) = QX_ROE
      dg_here%valx(3) = QY_ROE
      dg_here%valx(4) = bed_ROE

      dg_here%valy(1) = ZE_ROE
      dg_here%valy(2) = QX_ROE
      dg_here%valy(3) = QY_ROE
      dg_here%valy(4) = bed_ROE

      dtildeQx = 0.0
      dtildeQy = 0.0
      do w=1,4
         dtildeQx(w) = dg_here%porosity * evalf (w, dg_here%valx)  *( dg_here%QMag_IN**(2) + dg_here%QMag_EX**(2) ) / 2.D0 ! Interprete bytecode representation of w-th function
         dtildeQy(w) = dg_here%porosity*evalf2 (w, dg_here%valy) *( dg_here%QMag_IN**(2) + dg_here%QMag_EX**(2) ) / 2.D0 ! Interprete bytecode representation of w-th function
!$$$         if (EvalErrType > 0.or.EvalErrType2 > 0) then
!$$$            if(dg_here%porosity*((3*QX_ROE**2 + QX_ROE*QY_ROE + QY_ROE**2/2)/(ZE_ROE + bed_ROE)**3).ne.dtildeQx(w)) then
!$$$               write(*,*)'That was foolish on timestep', it
!$$$               WRITE(*,*)'*** Error: ',EvalErrMsg ()
!$$$               WRITE(*,*)'*** Error: ',EvalErrMsg2 ()
!$$$               WRITE(*,*) dg_here%funcx(w),'X=',dtildeQx(w),w
!$$$               WRITE(*,*) dg_here%funcy(w),'Y=',dtildeQy(w),w
!$$$               WRITE(*,*) 'dg_here%valx = ', dg_here%valx(w), w
!$$$               WRITE(*,*) 'dg_here%valy = ', dg_here%valy(w), w
!$$$               print*, '((3*QX_ROE**2 + QX_ROE*QY_ROE + QY_ROE**2/2)/(ZE_ROE + bed_ROE)**3)
!$$$     &              actualX=', dg_here%porosity*((3*QX_ROE**2 + QX_ROE*QY_ROE + QY_ROE**2/2)/(ZE_ROE + bed_ROE)**3)
!$$$               print*, '(QX_ROE*QY_ROE + QY_ROE**2/2)/(ZE_ROE + bed_ROE)**3  
!$$$     &              actualY=', dg_here%porosity*((QX_ROE*QY_ROE + QY_ROE**2/2)/(ZE_ROE + bed_ROE)**3 ) 
!$$$               write(*,*)'Ya shoulda known bettah!'    
!$$$            endif
!$$$         endif
      enddo


      
      EIGVAL(4) = 0.D0
      do w=1,4
         EIGVAL(4) = EIGVAL(4) + dtildeQx(w)*dg_here%NX + dtildeQy(w)*dg_here%NY
      enddo

!.....A bunch of functions

      epsilon = 0.D0 !1.0e-12

      chi = (dg_here%NY**2) * dtildeQy(2) + dtildeQx(3)*(dg_here%NY**2) - dtildeQx(3) &
     + (dtildeQx(2)-dtildeQy(3))*dg_here%NY*dg_here%NX

      xi = C_ROE**2

      alpha_0 = dg_here%NX*(dg_here%NY**2)*(dtildeQx(3)*xi + dtildeQy(2)*xi + dtildeQy(2)*(V_ROE**2))&
     + (dg_here%NY**3) * (dtildeQy(3)*xi - dtildeQy(2)*U_ROE*V_ROE - dtildeQx(2) * xi )&
     + dg_here%NY*(-dtildeQx(2)*(EIGVAL(2)**2)-dtildeQy(3)*(V_ROE**2) + dtildeQx(2)*xi&
     - dtildeQy(1)*V_ROE + (V_ROE**2) * dtildeQx(2)) + dtildeQx(2)*(-(dg_here%NY**2) * V_ROE &
     + dg_here%NY*EIGVAL(2))*EIGVAL(1) + (-dtildeQy(2)*V_ROE*EIGVAL(2)-dtildeQy(2)*V_ROE*C_ROE&
     + dtildeQx(1)*C_ROE)*dg_here%NY*dg_here%NX+(-dtildeQx(1)*V_ROE-dtildeQx(3)*(V_ROE**2))*dg_here%NX&
     + (dtildeQy(1)*C_ROE+dtildeQy(2)*C_ROE*U_ROE+dtildeQx(2)*V_ROE*EIGVAL(2)&
     + dtildeQx(2)*V_ROE*C_ROE)*(dg_here%NY**2) - dtildeQx(2)*V_ROE*EIGVAL(1)

      beta_0 = dg_here%NX*(dg_here%NY**2) *(dtildeQx(3)*xi+dtildeQy(2)*xi + dtildeQy(2)*(V_ROE**2))&
     + (dg_here%NY**3) * (dtildeQy(3)*xi - dtildeQy(2)*U_ROE*V_ROE - dtildeQx(2) * xi )&
     + dg_here%NY*(-dtildeQx(2)*(EIGVAL(2)**2)-dtildeQy(3)*(V_ROE**2) + dtildeQx(2)*xi&
     - dtildeQy(1)*V_ROE + (V_ROE**2) * dtildeQx(2)) + dtildeQx(2)*(-(dg_here%NY**2) * V_ROE &
     + dg_here%NY*EIGVAL(2))*EIGVAL(3) + (-dtildeQy(2)*V_ROE*EIGVAL(2)-dtildeQy(2)*V_ROE*C_ROE&
     + dtildeQx(1)*C_ROE)*dg_here%NY*dg_here%NX+(-dtildeQx(1)*V_ROE-dtildeQx(3)*(V_ROE**2))*dg_here%NX&
     + (-dtildeQy(1)*C_ROE-dtildeQy(2)*C_ROE*U_ROE+dtildeQx(2)*V_ROE*EIGVAL(2)&
     - dtildeQx(2)*V_ROE*C_ROE)*(dg_here%NY**2) + dtildeQx(2)*(-EIGVAL(2)+C_ROE)*V_ROE

      gamma1 = dg_here%NY * (-(V_ROE**2) * EIGVAL(2) - (V_ROE**2) * C_ROE + (C_ROE**3) - EIGVAL(2)*xi )&
    + V_ROE * xi + ( -EIGVAL(2) - C_ROE ) * (-dg_here%NY*(V_ROE**2) + EIGVAL(2)*V_ROE )&
    + (EIGVAL(2)**2) * V_ROE

      gamma2 = dg_here%NY * (-(V_ROE**2) * EIGVAL(2) + (V_ROE**2) * C_ROE - (C_ROE**3) - EIGVAL(2)*xi )&
    + V_ROE * xi + ( -EIGVAL(2) + C_ROE ) * (-dg_here%NY*(V_ROE**2) + EIGVAL(2)*V_ROE )&
    + (EIGVAL(2)**2) * V_ROE

      gamma3 = - dg_here%NX*dg_here%NY*V_ROE*xi + (dg_here%NY**2) * U_ROE * xi

      gamma4 = C_ROE*EIGVAL(2)*U_ROE - dg_here%NX*(C_ROE**3)

      iot = gamma2*EIGVAL(3)*dg_here%NY+dg_here%NY*gamma1*EIGVAL(1)+dg_here%NX*EIGVAL(3)*gamma4+dg_here%NX*EIGVAL(3)&
     *gamma3-dg_here%NX*EIGVAL(1)*gamma4+dg_here%NX*EIGVAL(1)*gamma3

      iot2 = - (dg_here%NY**2) * xi + (V_ROE**2)

      iot3 = -chi*(EIGVAL(2)**2) - chi*(EIGVAL(4)**2) + 2.D0*chi*EIGVAL(2)*EIGVAL(4)

      iot4 = beta_0*( (EIGVAL(4)**2) - 2.D0*EIGVAL(4)*EIGVAL(2) - C_ROE*EIGVAL(4)&
     + (EIGVAL(2)**2) +C_ROE*EIGVAL(2) )

      iot5 = -alpha_0*(EIGVAL(2)**2) + alpha_0*C_ROE*EIGVAL(2) - alpha_0*(EIGVAL(4)**2)&
     + 2.D0*alpha_0*EIGVAL(2)*EIGVAL(4) - alpha_0*C_ROE*EIGVAL(4)

      J6 = (dg_here%NY**2) * (xi**2) + iot2*xi

      J7 = iot3*xi + chi * (xi**2)

      J8 = beta_0*(EIGVAL(2)**2) - alpha_0*(EIGVAL(2)**2) - alpha_0*C_ROE*EIGVAL(4)&
     - 2.D0*beta_0*EIGVAL(2)*EIGVAL(4) + 2.D0*alpha_0*EIGVAL(2)*EIGVAL(4) + alpha_0*&
     C_ROE*EIGVAL(2) + beta_0*C_ROE*EIGVAL(2) - alpha_0*(EIGVAL(4)**2) + beta_0&
     *(EIGVAL(4)**2) - beta_0*C_ROE*EIGVAL(4)

      J9 = chi*J6 + iot3*(dg_here%NY**2)*xi + iot2*iot3

      bigA = EIGVAL(4) + C_ROE - EIGVAL(2)

      bigB = EIGVAL(2) + C_ROE - EIGVAL(4)

      bigD = -gamma2*gamma4*J9+gamma2*gamma3*J9-gamma4*gamma1*J9-gamma3*gamma1*J9&
     +gamma1*gamma2*dg_here%NY*J8+(gamma3*gamma1*J7+gamma4*gamma1*J7-gamma2*gamma3*J7&
     +gamma2*gamma4*J7)*(dg_here%NY**2)+gamma3*dg_here%NX*gamma1*iot5+gamma4*dg_here%NX*gamma1*iot5&
     +dg_here%NX*gamma3*gamma2*iot4-dg_here%NX*gamma4*gamma2*iot4

      bigE = -gamma1*EIGVAL(1)*J9-EIGVAL(3)*gamma2*J9+(-dg_here%NX*EIGVAL(1)*gamma3*J7-dg_here%NX*EIGVAL(3)&
     *gamma3*J7+dg_here%NX*EIGVAL(1)*gamma4*J7-dg_here%NX*EIGVAL(3)*gamma4*J7+iot*J7)*dg_here%NY&
     +dg_here%NX*EIGVAL(1)*gamma1*iot5-dg_here%NX*EIGVAL(3)*gamma2*iot4

      bigF = EIGVAL(1)*gamma3*J9+EIGVAL(3)*gamma3*J9-EIGVAL(1)*gamma4*J9+EIGVAL(3)*gamma4&
     *J9-gamma2*EIGVAL(3)*dg_here%NY*J8+(-EIGVAL(1)*gamma3*J7-EIGVAL(3)*gamma3*J7+EIGVAL(1)&
     *gamma4*J7-EIGVAL(3)*gamma4*J7)*(dg_here%NY**2)-dg_here%NX*EIGVAL(1)*gamma3*iot5-dg_here%NX*EIGVAL(3)&
     *gamma3*iot5+dg_here%NX*EIGVAL(1)*gamma4*iot5-dg_here%NX*EIGVAL(3)*gamma4*iot5+iot*iot5

      RI = 0.D0
      LE = 0.D0
      A_ROE = 0.D0
    
!.....Evaluate right eigenvectors at Roe averaged variables

      if (abs(EIGVAL(1)*alpha_0*gamma1).le.epsilon) then

         RI(1,1) = 0.D0

      else

         RI(1,1) = ( -C_ROE*(EIGVAL(2)**2) + (C_ROE**3) ) * ( EIGVAL(2) + C_ROE - EIGVAL(4) ) &
        *((dg_here%NY**2) * xi - (V_ROE**2) ) / (EIGVAL(1)*alpha_0*gamma1)

      endif

      if (abs(alpha_0*gamma1).le.epsilon) then

         RI(2,1)= 0.D0

      else

         RI(2,1) = ((dg_here%NY**2) * xi - (V_ROE**2)) * ( EIGVAL(2) + C_ROE - EIGVAL(4) )&
        *( -gamma4 + gamma3 ) / (alpha_0*gamma1)

      endif

      if (abs(alpha_0).le.epsilon) then

         RI(3,1) = 0.D0
      
      else
         
         RI(3,1) = (EIGVAL(2)+C_ROE-EIGVAL(4))*( (dg_here%NY*C_ROE)**2 - V_ROE**2 )/ alpha_0

      endif

      RI(4,1) = 1.D0

      RI(1,2) = 0.D0

      if (abs(chi).le.epsilon) then

         RI(2,2) = 0.D0
         RI(3,2) = 0.D0

      else

         RI(2,2) = (EIGVAL(2) - EIGVAL(4)) * dg_here%NY / chi
         RI(3,2) = (EIGVAL(4) - EIGVAL(2)) * dg_here%NX / chi

      endif

      RI(4,2) = 1.D0


      if ( abs(EIGVAL(3)*beta_0*gamma2).le.epsilon) then

         RI(1,3) = 0.D0

      else

         RI(1,3) = (C_ROE*EIGVAL(2)**2 - C_ROE**3)*(EIGVAL(2) - C_ROE - EIGVAL(4))&
        * ((dg_here%NY**2) * xi - V_ROE**2) / (EIGVAL(3)*beta_0*gamma2)

      endif

      if (abs(beta_0*gamma2).le.epsilon) then

         RI(2,3) = 0.D0

      else

         RI(2,3) = ((dg_here%NY**2) * xi - V_ROE**2)*(EIGVAL(2)-C_ROE-EIGVAL(4))&
        *(gamma4 + gamma3 ) / (beta_0*gamma2)

      endif

      if (abs(beta_0).le.epsilon) then
         
         RI(3,3) = 0.D0
      
      else
         
         RI(3,3) = (EIGVAL(2)-C_ROE-EIGVAL(4))*((dg_here%NY**2) * xi - V_ROE**2) / beta_0
         
      endif

      RI(4,3) = 1.D0

      RI(1,4) = 0.D0
      RI(2,4) = 0.D0
      RI(3,4) = 0.D0
      RI(4,4) = 1.D0

      
!.....Evaluate left eigenvectors (inverse) at Roe averaged variables

      if ( abs(C_ROE*(EIGVAL(2)+C_ROE-EIGVAL(4))*iot*(-(EIGVAL(2)**2)+xi)*iot2)&
  .le.epsilon) then

         LE(1,1) = 0.D0

      else

         LE(1,1) =  - ( alpha_0*EIGVAL(3)*EIGVAL(1)*gamma1*(dg_here%NY*gamma2+dg_here%NX*gamma4+dg_here%NX*gamma3))&
        / (C_ROE*bigB*iot*iot2*(-(EIGVAL(2)**2)+xi))


      endif

      if ( abs(C_ROE*(-(EIGVAL(2)**2)+xi)*(EIGVAL(2)-EIGVAL(4))*iot).le.epsilon) then

         LE(2,1) = 0.D0

      else

         LE(2,1) = chi*EIGVAL(1)*EIGVAL(3)*(gamma2*gamma4-gamma2*gamma3+gamma4*gamma1&
        +gamma3*gamma1) / (C_ROE*(-(EIGVAL(2)**2)+xi)*(EIGVAL(2)-EIGVAL(4))*iot)


      endif

      if ( abs(C_ROE*bigA*iot*(-(EIGVAL(2)**2)+xi)*iot2).le.epsilon) then

         LE(3,1) = 0.D0

      else

         LE(3,1) = -beta_0*EIGVAL(1)*EIGVAL(3)*gamma2*(gamma1*dg_here%NY-dg_here%NX*gamma4+dg_here%NX*gamma3)&
        / (C_ROE*bigA*iot*(-(EIGVAL(2)**2)+xi)*iot2)

      endif

      if ( abs(C_ROE*iot*bigA*bigB*(EIGVAL(2)-EIGVAL(4))*iot2&
 *(-(EIGVAL(2)**2)+xi)).le.epsilon) then

         LE(4,1) = 0.D0
         
      else
      
         LE(4,1) = (EIGVAL(1)*EIGVAL(3)*bigD) / (C_ROE*iot*iot2*bigA*bigB &
        *(EIGVAL(2)-EIGVAL(4))*(-(EIGVAL(2)**2)+xi))

      endif

      if ( abs( bigB*iot*iot2 ).le.epsilon) then

         LE(1,2) = 0.D0

      else

         LE(1,2) = -alpha_0*gamma1*EIGVAL(1)*dg_here%NX / ( bigB*iot*iot2 )

      endif

      if ( abs((EIGVAL(2)-EIGVAL(4))*iot).le.epsilon) then
      
         LE(2,2) = 0.D0

      else

         LE(2,2) = chi * ( EIGVAL(3)*gamma2 + gamma1*EIGVAL(1) ) / ((EIGVAL(2)-EIGVAL(4))*iot)

      endif

      if ( abs(bigA*iot*iot2).le.epsilon) then

         LE(3,2) = 0.D0
         
      else

         LE(3,2) = beta_0*gamma2*EIGVAL(3)*dg_here%NX / (bigA*iot*iot2)

      endif

      if ( abs(iot*iot2*bigA*bigB*(EIGVAL(2)-EIGVAL(4))).le.epsilon) then

         LE(4,2) = 0.D0

      else

         LE(4,2) = bigE / (iot*iot2*bigA*bigB*(EIGVAL(2)-EIGVAL(4)))

      endif

      if ( abs(bigB*iot*iot2).le.epsilon) then

         LE(1,3) = 0.D0
         
      else

         LE(1,3) = -alpha_0*(-gamma2*EIGVAL(3)*dg_here%NY-dg_here%NX*EIGVAL(3)*gamma4-dg_here%NX*EIGVAL(3)*gamma3&
        +dg_here%NX*EIGVAL(1)*gamma4-dg_here%NX*EIGVAL(1)*gamma3+iot) / (bigB*iot*iot2)

      endif

      if ( abs((EIGVAL(2)-EIGVAL(4))*iot).le.epsilon) then
         
         LE(2,3) = 0.D0
         
      else

         LE(2,3) = -chi*(EIGVAL(3)*gamma4+EIGVAL(3)*gamma3-EIGVAL(1)*gamma4+EIGVAL(1)*gamma3)&
        / ((EIGVAL(2)-EIGVAL(4))*iot)


      endif


      if ( abs(bigA*iot*iot2).le.epsilon) then 

         LE(3,3) = 0.D0
         
      else

         LE(3,3) = ( beta_0*gamma2*EIGVAL(3)*dg_here%NY ) / (bigA*iot*iot2)

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

!.....Compute the Roe matrix


       !Let's get rid of the bad stuff in the limit

      do ii=1,4
         if(abs(EIGVAL(ii)).le.1.0e-14) then

            RI(:,ii) = 0.D0
            LE(ii,:) = 0.D0

         endif
      enddo


      DO II = 1,4

         A_ROE(II,1) = RI(II,1)*ABS(EIGVAL(1))*LE(1,1) +&
        RI(II,2)*ABS(EIGVAL(2))*LE(2,1) +&
        RI(II,3)*ABS(EIGVAL(3))*LE(3,1) +&
        RI(II,4)*ABS(EIGVAL(4))*LE(4,1) 

         A_ROE(II,2) = RI(II,1)*ABS(EIGVAL(1))*LE(1,2) +&
        RI(II,2)*ABS(EIGVAL(2))*LE(2,2) +&
        RI(II,3)*ABS(EIGVAL(3))*LE(3,2) +&
        RI(II,4)*ABS(EIGVAL(4))*LE(4,2)

         A_ROE(II,3) = RI(II,1)*ABS(EIGVAL(1))*LE(1,3) +&
        RI(II,2)*ABS(EIGVAL(2))*LE(2,3) +&
        RI(II,3)*ABS(EIGVAL(3))*LE(3,3) +&
        RI(II,4)*ABS(EIGVAL(4))*LE(4,3)

         A_ROE(II,4) = RI(II,1)*ABS(EIGVAL(1))*LE(1,4) +&
        RI(II,2)*ABS(EIGVAL(2))*LE(2,4) +&
        RI(II,3)*ABS(EIGVAL(3))*LE(3,4) + &
        RI(II,4)*ABS(EIGVAL(4))*LE(4,4)

      ENDDO

      !Assume the dg_here%bed load

      bed_AVG(1) = 0.5D0*( (bedX_IN + bedX_EX)*dg_here%NX + (bedY_IN + bedY_EX)*dg_here%NY )

      JUMP(4) = dg_here%bed_EX(1) - dg_here%bed_IN(1)

      !Let's get rid of the bad stuff in the limit

      do ii=1,4
         if(abs(jump(ii)).le.1.0e-14) then

            A_ROE(:,ii) = 0.D0

         endif
      enddo


! Compute correction (this is not set by default)

!$$$      lambda1_IN = 0.D0
!$$$      lambda1_EX = 0.D0
!$$$      lambda2_IN = 0.D0
!$$$      lambda2_EX = 0.D0
!$$$      lambda3_IN = 0.D0
!$$$      lambda3_EX = 0.D0
!$$$
!$$$      lambda2_IN = (U_IN*dg_here%NX+V_IN*dg_here%NY)
!$$$      lambda2_EX = (U_EX*dg_here%NX+V_EX*dg_here%NY)
!$$$
!$$$      lambda1_IN = lambda2_IN + sqrt(g*HT_IN)
!$$$      lambda1_EX = lambda2_EX + sqrt(g*HT_EX)
!$$$
!$$$      lambda3_IN = lambda2_IN - sqrt(g*HT_IN)
!$$$      lambda3_EX = lambda2_EX - sqrt(g*HT_EX)
!$$$
!$$$      lambda4_IN = 0.D0
!$$$      lambda4_EX = 0.D0
!$$$
!$$$      lambda4_IN = (-3.D0*dg_here%QX_IN/(dg_here%ZE_IN + dg_here%bed_IN(1))**4 + (dg_here%ZE_IN + dg_here%bed_IN(1))**(-3)
!$$$     &     -3.D0*dg_here%QX_IN/(dg_here%ZE_IN + dg_here%bed_IN(1))**4)*dg_here%QMag_IN*dg_here%NX 
!$$$     &     + (-3.D0*dg_here%QY_IN/(dg_here%ZE_IN + dg_here%bed_IN(1))**4 + (dg_here%ZE_IN + dg_here%bed_IN(1))**(-3)
!$$$     &     -3.D0*dg_here%QY_IN/(dg_here%ZE_IN + dg_here%bed_IN(1))**4)*dg_here%QMag_IN*dg_here%NY 
!$$$
!$$$      lambda4_EX = (-3.D0*dg_here%QX_EX/(dg_here%ZE_EX + dg_here%bed_EX(1))**4 + (dg_here%ZE_EX + dg_here%bed_EX(1))**(-3)
!$$$     &     -3.D0*dg_here%QX_EX/(dg_here%ZE_EX + dg_here%bed_EX(1))**4)*dg_here%QMag_EX*dg_here%NX 
!$$$     &     + (-3.D0*dg_here%QY_EX/(dg_here%ZE_EX + dg_here%bed_EX(1))**4 + (dg_here%ZE_EX + dg_here%bed_EX(1))**(-3)
!$$$     &     -3.D0*dg_here%QY_EX/(dg_here%ZE_EX + dg_here%bed_EX(1))**4)*dg_here%QMag_EX*dg_here%NY 
!$$$
!$$$
!$$$
!$$$      Select_IN = min(min(lambda1_IN,lambda2_IN,lambda3_IN,lambda4_IN),
!$$$     &                min(lambda1_EX,lambda2_EX,lambda3_EX,lambda4_EX))
!$$$      Select_EX = max(max(lambda1_IN,lambda2_IN,lambda3_IN,lambda4_IN),
!$$$     &                max(lambda1_EX,lambda2_EX,lambda3_EX,lambda4_EX))
!$$$
!$$$      rVncU = 0.D0
!$$$      rVncV = 0.D0
!$$$
!$$$      do ll=1,layers !should be fixed for multiple layers
!$$$
!$$$         rVncU =  rVncU + 0.5D0*G*(HT_IN+HT_EX) 
!$$$     &        *(dg_here%bed_EX(ll)-dg_here%bed_IN(ll))*dg_here%NX
!$$$         rVncV =  rVncV + 0.5D0*G*(HT_IN+HT_EX) 
!$$$     &        *(dg_here%bed_EX(ll)-dg_here%bed_IN(ll))*dg_here%NY
!$$$         
!$$$        ! print*,  rVncU, rVncV
!$$$
!$$$      enddo
!$$$



!.....Compute the Roe flux (only formed for one sediment layer!) !Can add path conservative contribution

      F_HAT = F_AVG  - 0.5D0*(A_ROE(1,1)*JUMP(1) + A_ROE(1,2)*JUMP(2) +&
     A_ROE(1,3)*JUMP(3) + A_ROE(1,4) * JUMP(4) ) 
      !&        - (Select_IN*Select_EX)*(HT_IN-HT_EX) / (Select_EX-Select_IN)

      !print*, F_Hat

      G_HAT = G_AVG - 0.5D0*(A_ROE(2,1)*JUMP(1) + A_ROE(2,2)*JUMP(2) +&
     A_ROE(2,3)*JUMP(3) + A_ROE(2,4) * JUMP(4) )
      !&        - 0.5D0*(Select_IN+Select_EX)*rVncU
      !&        / (Select_EX-Select_IN)


      H_HAT = H_AVG - 0.5D0*(A_ROE(3,1)*JUMP(1) + A_ROE(3,2)*JUMP(2) +&
     A_ROE(3,3)*JUMP(3)+ A_ROE(3,4) * JUMP(4) )
      !&        - 0.5D0*(Select_IN+Select_EX)*rVncV
      !&        / (Select_EX-Select_IN)

      dg_here%bed_HAT(1) = bed_AVG(1) - 0.5D0*( A_ROE(4,1)*JUMP(1) + A_ROE(4,2)*JUMP(2) +&
     A_ROE(4,3)*JUMP(3) + A_ROE(4,4) * JUMP(4) )
      !&           - ((Select_IN*Select_EX) * (dg_here%bed_IN(1) - dg_here%bed_EX(1)))
      !&           / (Select_EX-Select_IN) 
      
#endif

!.....Compute upwinding scheme for chem and tracer

#ifdef TRACE
      U_AVG = 0.5D0 * ( U_IN + U_EX )
      V_AVG = 0.5D0 * ( V_IN + V_EX )
      VEL_NORMAL = U_AVG*dg_here%NX + V_AVG*dg_here%NY
      
      if (VEL_NORMAL.GT.0) THEN
         I_HAT = VEL_NORMAL*dg_here%iota_IN
      else
         I_HAT = VEL_NORMAL*dg_here%iota_EX
      endif
#endif

#ifdef CHEM
      U_AVG = 0.5D0 * ( U_IN + U_EX )
      V_AVG = 0.5D0 * ( V_IN + V_EX )
      VEL_NORMAL = U_AVG*dg_here%NX + V_AVG*dg_here%NY
      
      if (VEL_NORMAL.GT.0) THEN
         I_HAT = VEL_NORMAL*dg_here%iota_IN
         J_HAT = VEL_NORMAL*dg_here%iota2_IN
      else
         I_HAT = VEL_NORMAL*dg_here%iota_EX
         J_HAT = VEL_NORMAL*dg_here%iota2_EX
      endif
#endif

#ifdef DYNP
      U_AVG = 0.5D0 * ( U_IN + U_EX )
      V_AVG = 0.5D0 * ( V_IN + V_EX )
      VEL_NORMAL = U_AVG*dg_here%NX + V_AVG*dg_here%NY
      
      if (VEL_NORMAL.GT.0) THEN
         K_HAT = VEL_NORMAL*dynP_IN
      else
         K_HAT = VEL_NORMAL*dynP_EX
      endif
#endif

      RETURN
      END SUBROUTINE


!***********************************************************************
!     
!     SUBROUTINE LLF_FLUX
!     
!     This subroutine computes the Local Lax Friedrichs flux for the
!     shallow water equations.
!     
!     Written by Ethan Kubatko (07-15-2005)
!     
!***********************************************************************

      SUBROUTINE LLF_FLUX(s,dg_here)

!.....Use appropriate modules

      USE GLOBAL
      USE DG
      use sizes

      IMPLICIT NONE

      type (sizes_type) :: s
      type (dg_type) :: dg_here

!.....Declare local variables

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

!.....Compute the jump in the variables

      JUMP(1) = dg_here%ZE_EX - dg_here%ZE_IN
      JUMP(2) = dg_here%QX_EX - dg_here%QX_IN
      JUMP(3) = dg_here%QY_EX - dg_here%QY_IN

!.....Compute the total height of the water column

      HT_IN = dg_here%ZE_IN*NLEQ + dg_here%HB_IN

!.....Compute continuity fluxes at interior state

      F1_NL  = NLEQ + LEQ*dg_here%HB_IN
      FX1_IN = dg_here%QX_IN*F1_NL*dg_here%SFAC_IN
      FY1_IN = dg_here%QY_IN*F1_NL

!.....Compute momentum flux terms at interior state

      FU_NL = NLEQ*dg_here%QX_IN
      FV_NL = NLEQ*dg_here%QY_IN
      FG_NL = NLEQG*dg_here%ZE_IN
      FH_NL = 1.D0/(NLEQ*HT_IN + LEQ)
      U_IN  = dg_here%QX_IN*FH_NL
      V_IN  = dg_here%QY_IN*FH_NL

      HUU = FU_NL*U_IN
      HVV = FV_NL*V_IN
      HUV = FU_NL*V_IN
      GH2 = FG_NL*(0.5D0*dg_here%ZE_IN + dg_here%HB_IN) + dg_here%FG_L*dg_here%ZE_IN

!.....Compute x momentum fluxes at interior state

      FX2_IN = (HUU + GH2)*dg_here%SFAC_IN
      FY2_IN = HUV

!.....Compute y momentum fluxes at interior state

      FX3_IN = HUV*dg_here%SFAC_IN
      FY3_IN = HVV + GH2

!.....Compute the total height of the water column

      HT_EX = dg_here%ZE_EX*NLEQ + dg_here%HB_EX

!.....Compute continuity fluxes at exterior state

      F1_NL  = NLEQ + LEQ*dg_here%HB_EX
      FX1_EX = dg_here%QX_EX*F1_NL*dg_here%SFAC_EX
      FY1_EX = dg_here%QY_EX*F1_NL

!.....Compute momentum flux terms at interior state

      FU_NL = NLEQ*dg_here%QX_EX
      FV_NL = NLEQ*dg_here%QY_EX
      FG_NL = NLEQG*dg_here%ZE_EX
      FH_NL = 1.D0/(NLEQ*HT_EX + LEQ)
      U_EX  = dg_here%QX_EX*FH_NL
      V_EX  = dg_here%QY_EX*FH_NL

      HUU = FU_NL*U_EX
      HVV = FV_NL*V_EX
      HUV = FU_NL*V_EX
      GH2 = FG_NL*(0.5D0*dg_here%ZE_EX + dg_here%HB_EX) + dg_here%FG_L*dg_here%ZE_EX

!.....Compute x momentum fluxes at exterior state

      FX2_EX = (HUU + GH2)*dg_here%SFAC_EX
      FY2_EX = HUV

!.....Compute y momentum fluxes at exterior state

      FX3_EX = HUV*dg_here%SFAC_EX
      FY3_EX = HVV + GH2

!.....Compute the average flux function

      F1_AVG = 0.5D0*((FX1_IN + FX1_EX)*dg_here%NX + (FY1_IN + FY1_EX)*dg_here%NY)
      F2_AVG = 0.5D0*((FX2_IN + FX2_EX)*dg_here%NX + (FY2_IN + FY2_EX)*dg_here%NY)
      F3_AVG = 0.5D0*((FX3_IN + FX3_EX)*dg_here%NX + (FY3_IN + FY3_EX)*dg_here%NY)

!.....Evaluate the eigenvalues at the interior and exterior states

      UN_IN = (U_IN*dg_here%NX + V_IN*dg_here%NY)*NLEQ
      C_IN = SQRT(G*HT_IN*(dg_here%NY**2 + (dg_here%NX*dg_here%SFAC_IN)**2)) !srb - spherical coordinate correction
      EIGVALS(1) = ABS(UN_IN + C_IN)
      EIGVALS(2) = ABS(UN_IN)
      EIGVALS(3) = ABS(UN_IN - C_IN)

      UN_EX = (U_EX*dg_here%NX + V_EX*dg_here%NY)*NLEQ
      C_EX = SQRT(G*HT_EX*(dg_here%NY**2 + (dg_here%NX*dg_here%SFAC_EX)**2)) !srb - spherical coordinate correction
      EIGVALS(4) = ABS(UN_EX + C_EX)
      EIGVALS(5) = ABS(UN_EX)
      EIGVALS(6) = ABS(UN_EX - C_EX)

!.....Find the maximum eigenvalue (in absolute value)

      EIGMAX = MAX( EIGVALS(1), EIGVALS(2), EIGVALS(3),&
     EIGVALS(4), EIGVALS(5), EIGVALS(6) )

      !cfl_max = max(EIGMAX,cfl_max)
                                !print*,eigmax

!.....Compute the Local Lax Friedrichs Fluxes

      NX2 = NLEQ+LEQ*dg_here%NX*dg_here%NX
      NY2 = NLEQ+LEQ*dg_here%NY*dg_here%NY
      NXY = dg_here%NX*dg_here%NY

      F_HAT = F1_AVG - 0.5D0*EIGMAX*(JUMP(1))
      G_HAT = F2_AVG - 0.5D0*EIGMAX*(JUMP(2)*NX2 + LEQ*JUMP(3)*NXY)
      H_HAT = F3_AVG - 0.5D0*EIGMAX*(JUMP(3)*NY2 + LEQ*JUMP(2)*NXY)

!this is not well-tested!
#ifdef SED_LAY 

!.....Compute the stabilization term (v_nc in CMM's thesis) 
!.....for the nonconservative product in the momentum eqn

                                !Compute the magnitudes of the velocities for sediment

      dg_here%QMag_IN = (dg_here%QX_IN*dg_here%QX_IN/(HT_IN**2) + dg_here%QY_IN*dg_here%QY_IN/(HT_IN**2) )**(1/2)
      dg_here%QMag_EX = (dg_here%QX_EX*dg_here%QX_EX/(HT_EX**2) + dg_here%QY_EX*dg_here%QY_EX/(HT_EX**2) )**(1/2)

      rVncU = 0.D0
      rVncV = 0.D0

      do l=1,s%layers

         rVncU =  rVncU + 0.5D0*G*(HT_IN+HT_EX) &
        *(dg_here%bed_EX(l)-dg_here%bed_IN(l))*dg_here%NX
         rVncV =  rVncV + 0.5D0*G*(HT_IN+HT_EX) &
        *(dg_here%bed_EX(l)-dg_here%bed_IN(l))*dg_here%NY

      enddo

!.....When not using the NCP formalism
!.....then we enforce the (truncated) HLL wave-speed rendition.
!.....See CMM's thesis, pg 45, substituting the native flux

      lambda1_IN = U_IN*dg_here%NX + V_IN*dg_here%NY
      lambda1_EX = U_EX*dg_here%NX + V_EX*dg_here%NY

      lambda2_IN = lambda1_IN - (g*HT_IN)**(1/2)
      lambda2_EX = lambda1_EX - (g*HT_EX)**(1/2)

      lambda3_IN = lambda1_IN + (g*HT_IN)**(1/2)
      lambda3_EX = lambda1_EX + (g*HT_EX)**(1/2)

      Select_IN = min(lambda2_IN,lambda2_EX)
      Select_EX = max(lambda3_IN,lambda3_EX)

!.....Computes the sediment discharge rate based on the formulation
!.....developed by Camenen and Larson, so $\tilde{q} = A_g*H^{-3}|q|^{2}q$
!.....which clearly makes sense only for a single layer

      discharge_modelX_IN =  dg_here%porosity*HT_IN**(-1)*dg_here%QMag_IN**(2)*dg_here%QX_IN *dg_here%SFAC_IN
      discharge_modelX_EX =  dg_here%porosity*HT_EX**(-1)*dg_here%QMag_EX**(2)*dg_here%QX_EX *dg_here%SFAC_EX
      discharge_modelY_IN =  dg_here%porosity*HT_IN**(-1)*dg_here%QMag_IN**(2)*dg_here%QY_IN 
      discharge_modelY_EX =  dg_here%porosity*HT_EX**(-1)*dg_here%QMag_EX**(2)*dg_here%QY_EX



      q_RoeX = ( dg_here%QX_IN*(HT_IN)**0.5D0  &
     + dg_here%QX_EX*(HT_EX)**0.5D0 ) / ((HT_IN)**0.5D0 + (HT_EX)**0.5D0)
      q_RoeY = ( dg_here%QY_IN*(HT_IN)**0.5D0  &
     + dg_here%QY_EX*(HT_EX)**0.5D0 ) / ((HT_IN)**0.5D0 + (HT_EX)**0.5D0)

      q_Roe = q_RoeX*dg_here%NX + q_RoeY*dg_here%NY

      U_AVG = 0.5D0 * ( U_IN + U_EX )
      V_AVG = 0.5D0 * ( V_IN + V_EX )
      VEL_NORMAL = U_AVG*dg_here%NX + V_AVG*dg_here%NY

      do l=1,s%layers

         if (vel_normal.ge.0) THEN   !CAUTION, only makes sense for single layer
            dg_here%bed_HAT(l) = discharge_modelX_IN*dg_here%NX + discharge_modelY_IN*dg_here%NY
         else
            dg_here%bed_HAT(l) = discharge_modelX_EX*dg_here%NX + discharge_modelY_EX*dg_here%NY
         endif
         
      enddo
      
      
#endif

!.....Compute chem and tracer upwinding

#ifdef TRACE
      U_AVG = 0.5D0 * ( U_IN + U_EX )
      V_AVG = 0.5D0 * ( V_IN + V_EX )
      VEL_NORMAL = U_AVG*dg_here%NX + V_AVG*dg_here%NY
      
      if (VEL_NORMAL.GT.0) THEN
         I_HAT = VEL_NORMAL*dg_here%iota_IN
      else
         I_HAT = VEL_NORMAL*dg_here%iota_EX
      endif
#endif

#ifdef CHEM
      U_AVG = 0.5D0 * ( U_IN + U_EX )
      V_AVG = 0.5D0 * ( V_IN + V_EX )
      VEL_NORMAL = U_AVG*dg_here%NX + V_AVG*dg_here%NY
      
      if (VEL_NORMAL.GT.0) THEN
         I_HAT = VEL_NORMAL*dg_here%iota_IN
         J_HAT = VEL_NORMAL*dg_here%iota2_IN
      else
         I_HAT = VEL_NORMAL*dg_here%iota_EX
         J_HAT = VEL_NORMAL*dg_here%iota2_EX
      endif
#endif

#ifdef DYNP
      U_AVG = 0.5D0 * ( U_IN + U_EX )
      V_AVG = 0.5D0 * ( V_IN + V_EX )
      VEL_NORMAL = U_AVG*dg_here%NX + V_AVG*dg_here%NY
      
      if (VEL_NORMAL.GT.0) THEN
         K_HAT = VEL_NORMAL*dynP_IN
      else
         K_HAT = VEL_NORMAL*dynP_EX
      endif
#endif

      RETURN
      END SUBROUTINE
      
!***********************************************************************
!     
!     SUBROUTINE HLLC_FLUX()
!     
!     This subroutine computes the Harten-Lax-van Leer Contact (HLLC)
!     flux for the shallow water equations.
!     
!     Written by Ethan Kubatko
!     
!***********************************************************************

      SUBROUTINE HLLC_FLUX(s,dg_here)

!.....Use appropriate modules

      USE GLOBAL
      USE DG
      use sizes

      IMPLICIT NONE

      type (sizes_type) :: s
      type (dg_type) :: dg_here

!.....Declare local variables

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

!.....Compute the interior and exterior water column heights

      HT_IN = dg_here%ZE_IN*NLEQ + dg_here%HB_IN 
      HT_EX = dg_here%ZE_EX*NLEQ + dg_here%HB_EX

!.....Compute the interior and exterior normal velocities

      FH_NL = 1.D0/(NLEQ*HT_IN + LEQ)
      U_IN  = dg_here%QX_IN*FH_NL
      V_IN  = dg_here%QY_IN*FH_NL
      UN_IN =  U_IN*dg_here%NX + V_IN*dg_here%NY
      UT_IN = -U_IN*dg_here%NY + V_IN*dg_here%NX

      FH_NL = 1.D0/(NLEQ*HT_EX + LEQ)
      U_EX  = dg_here%QX_EX*FH_NL
      V_EX  = dg_here%QY_EX*FH_NL
      UN_EX =  U_EX*dg_here%NX + V_EX*dg_here%NY
      UT_EX = -U_EX*dg_here%NY + V_EX*dg_here%NX

!.....Compute H*

      C_EX = SQRT(G*HT_EX)
      C_IN = SQRT(G*HT_IN)
      H_STAR = 0.5D0*(C_IN + C_EX) + 0.25D0*(UN_IN - UN_EX)
      H_STAR = 1.D0/G*H_STAR*H_STAR

!.....Compute the wave speeds

      LS_EX = 1.D0
      LS_IN = 1.D0
      IF (H_STAR.GT.HT_EX) LS_EX = SQRT( (H_STAR*H_STAR + H_STAR*HT_EX)&
     /(2.D0*HT_EX*HT_EX) )
      IF (H_STAR.GT.HT_IN) LS_IN = SQRT( (H_STAR*H_STAR + H_STAR*HT_IN)&
     /(2.D0*HT_IN*HT_IN) )

      S_IN   = UN_IN - C_IN*LS_IN
      S_EX   = UN_EX + C_EX*LS_EX
      S_STAR = ( S_IN*HT_EX*(UN_EX - S_EX) - S_EX*HT_IN*(UN_IN - S_IN) )&
     /( HT_EX*(UN_EX - S_EX) - HT_IN*(UN_IN - S_IN) )

!.....Compute the numerical flux based on the wave speeds

      IF (S_IN.GE.0) THEN

!.......Compute continuity fluxes at exterior state

         F1_NL  = NLEQ + LEQ*dg_here%HB_IN
         FX1_IN = dg_here%QX_IN*F1_NL
         FY1_IN = dg_here%QY_IN*F1_NL

!.......Compute momentum flux terms at exterior state

         FU_NL = NLEQ*dg_here%QX_IN
         FV_NL = NLEQ*dg_here%QY_IN
         FG_NL = NLEQG*dg_here%ZE_IN

         HUU = FU_NL*U_IN
         HVV = FV_NL*V_IN
         HUV = FU_NL*V_IN
         GH2 = FG_NL*(0.5D0*dg_here%ZE_IN + dg_here%HB_IN) + dg_here%FG_L*dg_here%ZE_IN

!.......Compute x momentum fluxes at exterior state

         FX2_IN = HUU + GH2
         FY2_IN = HUV

!.......Compute y momentum fluxes at exterior state

         FX3_IN = HUV
         FY3_IN = HVV + GH2

!.......Set the numerical fluxes equal to the exterior fluxes

         F_HAT = FX1_IN*dg_here%NX + FY1_IN*dg_here%NY
         G_HAT = FX2_IN*dg_here%NX + FY2_IN*dg_here%NY
         H_HAT = FX3_IN*dg_here%NX + FY3_IN*dg_here%NY

      ELSEIF ((S_IN.LE.0).AND.(S_STAR.GE.0)) THEN

!.......Compute continuity fluxes at exterior state

         F1_NL  = NLEQ + LEQ*dg_here%HB_IN
         FX1_IN = dg_here%QX_IN*F1_NL
         FY1_IN = dg_here%QY_IN*F1_NL

!.......Compute momentum flux terms at exterior state

         FU_NL = NLEQ*dg_here%QX_IN
         FV_NL = NLEQ*dg_here%QY_IN
         FG_NL = NLEQG*dg_here%ZE_IN

         HUU = FU_NL*U_IN
         HVV = FV_NL*V_IN
         HUV = FU_NL*V_IN
         GH2 = FG_NL*(0.5D0*dg_here%ZE_IN + dg_here%HB_IN) + dg_here%FG_L*dg_here%ZE_IN

!.......Compute x momentum fluxes at exterior state

         FX2_IN = HUU + GH2
         FY2_IN = HUV

!.......Compute y momentum fluxes at exterior state

         FX3_IN = HUV
         FY3_IN = HVV + GH2

!.......Compute Q*

         Q_STAR(1) = HT_IN*((S_IN - UN_IN)/(S_IN - S_STAR))
         Q_STAR(2) = Q_STAR(1)*S_STAR
         Q_STAR(3) = Q_STAR(1)*UT_IN
         Q(1) =  HT_IN
         Q(2) =  dg_here%QX_IN*dg_here%NX + dg_here%QY_IN*dg_here%NY
         Q(3) = -dg_here%QX_IN*dg_here%NY + dg_here%QY_IN*dg_here%NX

!.......Compute the numerical flux

         F_HAT = FX1_IN*dg_here%NX + FY1_IN*dg_here%NY + S_IN*(  Q_STAR(1) - Q(1))
         G_HAT = FX2_IN*dg_here%NX + FY2_IN*dg_here%NY + S_IN*( (Q_STAR(2) - Q(2))*dg_here%NX -&
        (Q_STAR(3) - Q(3))*dg_here%NY )
         H_HAT = FX3_IN*dg_here%NX + FY3_IN*dg_here%NY + S_IN*( (Q_STAR(2) - Q(2))*dg_here%NY +&
        (Q_STAR(3) - Q(3))*dg_here%NX )

      ELSEIF ((S_STAR.LE.0).AND.(S_EX.GE.0)) THEN

!.......Compute continuity fluxes at exterior state

         F1_NL  = NLEQ + LEQ*dg_here%HB_EX
         FX1_EX = dg_here%QX_EX*F1_NL
         FY1_EX = dg_here%QY_EX*F1_NL

!.......Compute momentum flux terms at exterior state

         FU_NL = NLEQ*dg_here%QX_EX
         FV_NL = NLEQ*dg_here%QY_EX
         FG_NL = NLEQG*dg_here%ZE_EX

         HUU = FU_NL*U_EX
         HVV = FV_NL*V_EX
         HUV = FU_NL*V_EX
         GH2 = FG_NL*(0.5D0*dg_here%ZE_EX + dg_here%HB_EX) + dg_here%FG_L*dg_here%ZE_EX

!.......Compute x momentum fluxes at exterior state

         FX2_EX = HUU + GH2
         FY2_EX = HUV

!.......Compute y momentum fluxes at exterior state

         FX3_EX = HUV
         FY3_EX = HVV + GH2

!.......Compute Q*

         Q_STAR(1) = HT_EX*((S_EX - UN_EX)/(S_EX - S_STAR))
         Q_STAR(2) = Q_STAR(1)*S_STAR
         Q_STAR(3) = Q_STAR(1)*UT_EX
         Q(1) =  HT_EX
         Q(2) =  dg_here%QX_EX*dg_here%NX + dg_here%QY_EX*dg_here%NY
         Q(3) = -dg_here%QX_EX*dg_here%NY + dg_here%QY_EX*dg_here%NX

!.......Compute the numerical flux

         F_HAT = FX1_EX*dg_here%NX + FY1_EX*dg_here%NY + S_EX*(  Q_STAR(1) - Q(1))
         G_HAT = FX2_EX*dg_here%NX + FY2_EX*dg_here%NY + S_EX*( (Q_STAR(2) - Q(2))*dg_here%NX -&
        (Q_STAR(3) - Q(3))*dg_here%NY )
         H_HAT = FX3_EX*dg_here%NX + FY3_EX*dg_here%NY + S_EX*( (Q_STAR(2) - Q(2))*dg_here%NY +&
        (Q_STAR(3) - Q(3))*dg_here%NX )

      ELSEIF (S_EX.LE.0) THEN

!.......Compute continuity fluxes at exterior state

         F1_NL  = NLEQ + LEQ*dg_here%HB_EX
         FX1_EX = dg_here%QX_EX*F1_NL
         FY1_EX = dg_here%QY_EX*F1_NL

!.......Compute momentum flux terms at exterior state

         FU_NL = NLEQ*dg_here%QX_EX
         FV_NL = NLEQ*dg_here%QY_EX
         FG_NL = NLEQG*dg_here%ZE_EX

         HUU = FU_NL*U_EX
         HVV = FV_NL*V_EX
         HUV = FU_NL*V_EX
         GH2 = FG_NL*(0.5D0*dg_here%ZE_EX + dg_here%HB_EX) + dg_here%FG_L*dg_here%ZE_EX

!.......Compute x momentum fluxes at exterior state

         FX2_EX = HUU + GH2
         FY2_EX = HUV

!.......Compute y momentum fluxes at exterior state

         FX3_EX = HUV
         FY3_EX = HVV + GH2

!.......Set the numerical fluxes equal to the exterior fluxes

         F_HAT = FX1_EX*dg_here%NX + FY1_EX*dg_here%NY
         G_HAT = FX2_EX*dg_here%NX + FY2_EX*dg_here%NY
         H_HAT = FX3_EX*dg_here%NX + FY3_EX*dg_here%NY

      ENDIF

!Not thoroughly tested!!
#ifdef SED_LAY 

!.....Compute the stabilization term (v_nc in CMM's thesis) 
!.....for the nonconservative product in the momentum eqn

                                !Compute the magnitudes of the velocities for sediment

      dg_here%QMag_IN = (dg_here%QX_IN*dg_here%QX_IN/(HT_IN**2) + dg_here%QY_IN*dg_here%QY_IN/(HT_IN**2) )**(1/2)
      dg_here%QMag_EX = (dg_here%QX_EX*dg_here%QX_EX/(HT_EX**2) + dg_here%QY_EX*dg_here%QY_EX/(HT_EX**2) )**(1/2)

      rVncU = 0.D0
      rVncV = 0.D0

      do l=1,s%layers

         rVncU =  rVncU + 0.5D0*G*(HT_IN+HT_EX) &
        *(dg_here%bed_EX(l)-dg_here%bed_IN(l))*dg_here%NX
         rVncV =  rVncV + 0.5D0*G*(HT_IN+HT_EX) &
        *(dg_here%bed_EX(l)-dg_here%bed_IN(l))*dg_here%NY

      enddo

!.....When not using the NCP formalism
!.....then we enforce the (truncated) HLL wave-speed rendition.
!.....See CMM's thesis, pg 45, substituting for the base flux

      lambda1_IN = U_IN*dg_here%NX + V_IN*dg_here%NY
      lambda1_EX = U_EX*dg_here%NX + V_EX*dg_here%NY

      lambda2_IN = lambda1_IN - (g*HT_IN)**(1/2)
      lambda2_EX = lambda1_EX - (g*HT_EX)**(1/2)

      lambda3_IN = lambda1_IN + (g*HT_IN)**(1/2)
      lambda3_EX = lambda1_EX + (g*HT_EX)**(1/2)

      Select_IN = min(lambda2_IN,lambda2_EX)
      Select_EX = max(lambda3_IN,lambda3_EX)

      if(Select_IN.ge.0.D0) then 
         
         !F_HAT = FX1_IN*dg_here%NX + FY1_IN*dg_here%NY No update
         G_HAT = G_HAT  - 0.5D0*rVncU
         H_HAT = H_HAT  - 0.5D0*rVncV
         
      elseif(Select_IN.lt.0.D0.and.Select_EX.gt.0.D0) then

         !F_HAT = FX1_EX*dg_here%NX + FY1_EX*dg_here%NY  No update

         G_HAT = G_HAT - 0.5D0*(Select_IN + Select_EX)/&
        (Select_EX-Select_IN)*rVncU
         H_HAT = H_HAT - 0.5D0*(Select_IN + Select_EX)/&
        (Select_EX-Select_IN)*rVncV
         
      else

         !F_HAT = FX1_EX*dg_here%NX + FY1_EX*dg_here%NY 
         G_HAT = G_HAT  + 0.5D0*rVncU
         H_HAT = H_HAT  + 0.5D0*rVncV

      endif

!.....Computes the sediment discharge rate based on the formulation
!.....developed by Camenen and Larson, so $\tilde{q} = A_g*H^{-3}|q|^{2}q$
!.....which clearly makes sense only for a single layer

      discharge_modelX_IN =  dg_here%porosity*HT_IN**(-1)*dg_here%QMag_IN**(2)*dg_here%QX_IN*dg_here%SFAC_IN
      discharge_modelX_EX =  dg_here%porosity*HT_EX**(-1)*dg_here%QMag_EX**(2)*dg_here%QX_EX*dg_here%SFAC_EX
      discharge_modelY_IN =  dg_here%porosity*HT_IN**(-1)*dg_here%QMag_IN**(2)*dg_here%QY_IN 
      discharge_modelY_EX =  dg_here%porosity*HT_EX**(-1)*dg_here%QMag_EX**(2)*dg_here%QY_EX


      q_RoeX = ( dg_here%QX_IN*(HT_IN)**0.5D0  &
     + dg_here%QX_EX*(HT_EX)**0.5D0 ) / ((HT_IN)**0.5D0 + (HT_EX)**0.5D0)
      q_RoeY = ( dg_here%QY_IN*(HT_IN)**0.5D0  &
     + dg_here%QY_EX*(HT_EX)**0.5D0 ) / ((HT_IN)**0.5D0 + (HT_EX)**0.5D0)

      q_Roe = q_RoeX*dg_here%NX + q_RoeY*dg_here%NY

      do l=1,s%layers

         if (q_Roe.ge.0) THEN   !CAUTION, only makes sense for single layer
            dg_here%bed_HAT(l) = discharge_modelX_IN*dg_here%NX + discharge_modelY_IN*dg_here%NY
         else
            dg_here%bed_HAT(l) = discharge_modelX_EX*dg_here%NX + discharge_modelY_EX*dg_here%NY
         endif
         
      enddo
      
      
#endif

#ifdef TRACE
      U_AVG = 0.5D0 * ( U_IN + U_EX )
      V_AVG = 0.5D0 * ( V_IN + V_EX )
      VEL_NORMAL = U_AVG*dg_here%NX + V_AVG*dg_here%NY
      
      if (VEL_NORMAL.GT.0) THEN
         I_HAT = VEL_NORMAL*dg_here%iota_IN
      else
         I_HAT = VEL_NORMAL*dg_here%iota_EX
      endif
#endif

#ifdef CHEM
      U_AVG = 0.5D0 * ( U_IN + U_EX )
      V_AVG = 0.5D0 * ( V_IN + V_EX )
      VEL_NORMAL = U_AVG*dg_here%NX + V_AVG*dg_here%NY
      
      if (VEL_NORMAL.GT.0) THEN
         I_HAT = VEL_NORMAL*dg_here%iota_IN
         J_HAT = VEL_NORMAL*dg_here%iota2_IN
      else
         I_HAT = VEL_NORMAL*dg_here%iota_EX
         J_HAT = VEL_NORMAL*dg_here%iota2_EX
      endif
#endif

#ifdef DYNP
      U_AVG = 0.5D0 * ( U_IN + U_EX )
      V_AVG = 0.5D0 * ( V_IN + V_EX )
      VEL_NORMAL = U_AVG*dg_here%NX + V_AVG*dg_here%NY
      
      if (VEL_NORMAL.GT.0) THEN
         K_HAT = VEL_NORMAL*dynP_IN
      else
         K_HAT = VEL_NORMAL*dynP_EX
      endif
#endif

      
      Return
      END SUBROUTINE


!***********************************************************************
!     
!     subroutine NCP_flux()
!     
!     This subroutine uses CMM's version of the NCP flux.  
!     See his thesis pages 2.4.3
!     
!     NOTE: THIS FLUX IS ONLY MEANT FOR WHEN SEDIMENT IS ON!!
!     because the formalism is more tightly coupled
!     
!     
!     Written by cem interpreted from CM's code revamp
!     
!***********************************************************************

      SUBROUTINE NCP_FLUX(s,dg_here)

!.....Use appropriate modules

      USE GLOBAL
      USE DG
      use sizes

      IMPLICIT NONE

      type (sizes_type) :: s
      type (dg_type) :: dg_here

!.....Declare local variables

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

!.....Compute the interior and exterior water column heights

      HT_IN = dg_here%ZE_IN*IFNLFA + dg_here%HB_IN
      HT_EX = dg_here%ZE_EX*IFNLFA + dg_here%HB_EX

!.....Compute the interior and exterior normal velocities

      FH_NL_IN = 1.D0/(NLEQ*HT_IN + LEQ)
      U_IN  = dg_here%QX_IN*FH_NL_IN
      V_IN  = dg_here%QY_IN*FH_NL_IN
 
      FH_NL_EX = 1.D0/(NLEQ*HT_EX + LEQ)
      U_EX  = dg_here%QX_EX*FH_NL_EX
      V_EX  = dg_here%QY_EX*FH_NL_EX
 
!.....Compute the interior and exterior discharge rates

      dg_here%QMag_IN = (dg_here%QX_IN*dg_here%QX_IN/(HT_IN**2) + dg_here%QY_IN*dg_here%QY_IN/(HT_IN**2) )**(1/2)
      dg_here%QMag_EX = (dg_here%QX_EX*dg_here%QX_EX/(HT_EX**2) + dg_here%QY_EX*dg_here%QY_EX/(HT_EX**2) )**(1/2)

      discharge_modelX_IN =  dg_here%porosity*HT_IN**(-1)*dg_here%QMag_IN**(2)*(dg_here%QX_IN) *dg_here%SFAC_IN
      discharge_modelX_EX =  dg_here%porosity*HT_EX**(-1)*dg_here%QMag_EX**(2)*(dg_here%QX_EX) *dg_here%SFAC_EX
      discharge_modelY_IN =  dg_here%porosity*HT_IN**(-1)*dg_here%QMag_IN**(2)*(dg_here%QX_IN) 
      discharge_modelY_EX =  dg_here%porosity*HT_EX**(-1)*dg_here%QMag_EX**(2)*(dg_here%QY_EX)

!.....Compute the eigenvalues and aux variables of interior states

                                !Army of definitions

      lambda1_IN = U_IN*dg_here%NX + V_IN*dg_here%NY
 
      rAlph = dg_here%porosity ! constant
      kappa2 = 1.D0 !kapp2 = 1.D0 - dg_here%porosity

      rB = -2.D0*lambda1_IN !hardwired for \tilde{q} independent of b!

      rC = (1.D0 - 2.D0*rAlph*G/(kappa2) )* lambda1_IN**2 &
 - G*HT_IN*(HT_IN**(-3)*( dg_here%QMag_IN * rAlph/kappa2) + 1.D0)


      rD = rAlph*G*lambda1_IN/kappa2  &
 * (3.D0*HT_IN**(-2)*dg_here%QMag_IN - 2.D0*lambda1_IN**2)


      rQ = (3.D0*rC - rB**2) / 9.D0

      rR = (9.D0*rB*rC - 27.D0*rD - 2.D0*rB**3) / 54.D0

      rTheta = acos(rR / sqrt(-rQ**3))

      lambda2_IN = 2.D0*sqrt(-rQ)*cos(rTheta/3.D0)+ (2.D0/3.D0)*lambda1_IN

      lambda3_IN = 2.D0*sqrt(-rQ)*cos(rTheta/3.D0+ 2.D0*PI/3.D0)+ (2.D0/3.D0)*lambda1_IN

      lambda4_IN = 2.D0*sqrt(-rQ)*cos(rTheta/3.D0+ 4.D0*PI/3.D0)+ (2.D0/3.D0)*lambda1_IN

!.....Compute the eigenvalues and aux vairables of exterior states

                                !Army of definitions

      lambda1_EX = U_EX*dg_here%NX + V_EX*dg_here%NY

      rB = -2.D0*lambda1_EX !hardwired for \tilde{q} independent of b!

      rC = (1.D0 - 2.D0*rAlph*G/(kappa2) )* lambda1_EX**2 &
 - G*HT_EX*(rAlph*HT_EX**(-3)*( dg_here%QMag_EX * rAlph/kappa2) + 1.D0)

      rD = rAlph*G*lambda1_EX/kappa2  &
 * (3.D0*HT_EX**(-2)*dg_here%QMag_EX - 2.D0*lambda1_EX**2)

      rQ = (3.D0*rC - rB**2) / 9.D0

      rR = (9.D0*rB*rC - 27.D0*rD - 2.D0*rB**3) / 54.D0

      rTheta = acos(rR / sqrt(-rQ**3))

      lambda2_EX = 2.D0*sqrt(-rQ)*cos(rTheta/3.D0)+ (2.D0/3.D0)*lambda1_EX

      lambda3_EX = 2.D0*sqrt(-rQ)*cos(rTheta/3.D0+ 2.D0*PI/3.D0)+ (2.D0/3.D0)*lambda1_EX

      lambda4_EX = 2.D0*sqrt(-rQ)*cos(rTheta/3.D0+ 4.D0*PI/3.D0)+ (2.D0/3.D0)*lambda1_EX

!.....Select the correction displacement stuff for the momentum

      Select_IN = min(min(lambda1_IN,lambda2_IN,lambda3_IN,lambda4_IN),&
                min(lambda1_EX,lambda2_EX,lambda3_EX,lambda4_EX))
      Select_EX = max(max(lambda1_IN,lambda2_IN,lambda3_IN,lambda4_IN),&
                max(lambda1_EX,lambda2_EX,lambda3_EX,lambda4_EX))
      rVncU = 0.D0
      rVncV = 0.D0

      do ll=1,s%layers !should be fixed for multiple layers

         rVncU =  rVncU + 0.5D0*G*(HT_IN+HT_EX) &
        *(dg_here%bed_EX(ll)-dg_here%bed_IN(ll))*dg_here%NX
         rVncV =  rVncV + 0.5D0*G*(HT_IN+HT_EX) &
        *(dg_here%bed_EX(ll)-dg_here%bed_IN(ll))*dg_here%NY

      enddo

!.....First Case

      if (Select_IN.ge.0) then

!.......Compute continuity fluxes at interior state

         F1_NL  = NLEQ + LEQ*dg_here%HB_IN
         FX1_IN = dg_here%QX_IN*F1_NL*dg_here%SFAC_IN
         FY1_IN = dg_here%QY_IN*F1_NL

!.......Compute momentum flux terms at interior state

         FU_NL = NLEQ*dg_here%QX_IN
         FV_NL = NLEQ*dg_here%QY_IN
         FG_NL = NLEQG*dg_here%ZE_IN

         HUU = FU_NL*U_IN
         HVV = FV_NL*V_IN
         HUV = FU_NL*V_IN
         GH2 = FG_NL*(0.5D0*dg_here%ZE_IN + dg_here%HB_IN) + dg_here%FG_L*dg_here%ZE_IN
         GH3 = 0.5D0*G*(HT_IN**2)

!.......Compute x momentum fluxes at interior state

         FX2_IN = (HUU + GH3)*dg_here%SFAC_IN
         FY2_IN = HUV

!.......Compute y momentum fluxes at interior state

         FX3_IN = HUV*dg_here%SFAC_IN
         FY3_IN = HVV + GH3

!.......Set the fluxes
         
         F_HAT = FX1_IN*dg_here%NX + FY1_IN*dg_here%NY 
         G_HAT = FX2_IN*dg_here%NX + FY2_IN*dg_here%NY - 0.5D0*rVncU
         H_HAT = FX3_IN*dg_here%NX + FY3_IN*dg_here%NY - 0.5D0*rVncV

         do ll = 1,s%layers        !CAUTION, only really makes sense for single layer

               dg_here%bed_HAT(ll) = 0.5D0*(discharge_modelX_IN+discharge_modelX_EX)*dg_here%NX &
           + 0.5D0*(discharge_modelY_IN+discharge_modelY_EX)*dg_here%NY

         enddo

 
!...Second case
      elseif(Select_IN.lt.0.D0.and.Select_EX.ge.0.D0) then

!.......Compute continuity fluxes

         F1_NL  = NLEQ + LEQ*dg_here%HB_IN
         FX1_IN = dg_here%QX_IN*F1_NL*dg_here%SFAC_IN
         FY1_IN = dg_here%QY_IN*F1_NL
         F1_NL  = NLEQ + LEQ*dg_here%HB_EX
         FX1_EX = dg_here%QX_EX*F1_NL*dg_here%SFAC_EX
         FY1_EX = dg_here%QY_EX*F1_NL

!.......Compute momentum fluxes interior

         FU_NL = NLEQ*dg_here%QX_IN
         FV_NL = NLEQ*dg_here%QY_IN
         FG_NL = NLEQG*dg_here%ZE_IN

         HUU = FU_NL*U_IN
         HVV = FV_NL*V_IN
         HUV = FU_NL*V_IN
         GH2 = FG_NL*(0.5D0*dg_here%ZE_IN + dg_here%HB_IN) + dg_here%FG_L*dg_here%ZE_IN
         GH3 = 0.5D0*G*(HT_IN**2)

!.......Compute x momentum fluxes interior

         FX2_IN = (HUU + GH3)*dg_here%SFAC_IN
         FY2_IN = HUV

!.......Compute y momentum fluxes interior

         FX3_IN = HUV*dg_here%SFAC_IN
         FY3_IN = HVV + GH3

!.......Compute momentum fluxes exterior 

         FU_NL = NLEQ*dg_here%QX_EX
         FV_NL = NLEQ*dg_here%QY_EX
         FG_NL = NLEQG*dg_here%ZE_EX

         HUU = FU_NL*U_EX
         HVV = FV_NL*V_EX
         HUV = FU_NL*V_EX
         GH2 = FG_NL*(0.5D0*dg_here%ZE_EX + dg_here%HB_EX) + dg_here%FG_L*dg_here%ZE_EX
         GH3 = 0.5D0*G*(HT_EX**2)

!.......Compute x momentum fluxes exterior

         FX2_EX = (HUU + GH3)*dg_here%SFAC_EX
         FY2_EX = HUV

!.......Compute y momentum fluxes at exterior state

         FX3_EX = HUV*dg_here%SFAC_EX
         FY3_EX = HVV + GH3


!.......Assign NCP/HLL version of fluxes
         
         F_HAT = ((Select_EX*(FX1_IN*dg_here%NX + FY1_IN*dg_here%NY) - Select_IN*(FX1_EX*dg_here%NX + FY1_EX*dg_here%NY))&
              - (Select_IN*Select_EX)*(HT_IN-HT_EX) ) / (Select_EX-Select_IN)
         
         G_HAT = ((Select_EX*(FX2_IN*dg_here%NX + FY2_IN*dg_here%NY) - Select_IN*(FX2_EX*dg_here%NX + FY2_EX*dg_here%NY))&
              - (Select_IN*Select_EX)*(dg_here%QX_IN-dg_here%QX_EX) - 0.5D0*(Select_IN+Select_EX)*rVncU)&
              / (Select_EX-Select_IN)
         
         H_HAT =  ((Select_EX*(FX3_IN*dg_here%NX + FY3_IN*dg_here%NY) - Select_IN*(FX3_EX*dg_here%NX + FY3_EX*dg_here%NY))&
              - (Select_IN*Select_EX)*(dg_here%QY_IN-dg_here%QY_EX) - 0.5D0*(Select_IN+Select_EX)*rVncV )&
              / (Select_EX-Select_IN)


         do ll = 1,s%layers 

            dg_here%bed_HAT(ll) = ((Select_EX*(discharge_modelX_IN*dg_here%NX + discharge_modelY_IN*dg_here%NY)&
           - Select_IN*(discharge_modelX_EX*dg_here%NX + discharge_modelY_EX*dg_here%NY))) &
           - ((Select_IN*Select_EX) * (dg_here%bed_IN(ll) - dg_here%bed_EX(ll)))&
           / (Select_EX-Select_IN) 
         enddo


!...Third case 
      else
   

!.......Compute continuity fluxes at exterior state

         F1_NL  = NLEQ + LEQ*dg_here%HB_EX
         FX1_EX = dg_here%QX_EX*F1_NL*dg_here%SFAC_EX
         FY1_EX = dg_here%QY_EX*F1_NL

!.......Compute momentum flux terms at exterior state

         FU_NL = NLEQ*dg_here%QX_EX
         FV_NL = NLEQ*dg_here%QY_EX
         FG_NL = NLEQG*dg_here%ZE_EX

         HUU = FU_NL*U_EX
         HVV = FV_NL*V_EX
         HUV = FU_NL*V_EX
         GH2 = FG_NL*(0.5D0*dg_here%ZE_EX + dg_here%HB_EX) + dg_here%FG_L*dg_here%ZE_EX
         GH3 = 0.5D0*G*(HT_EX**2)

!.......Compute x momentum fluxes at exterior state

         FX2_EX = (HUU + GH3)*dg_here%SFAC_EX
         FY2_EX = HUV

!.......Compute y momentum fluxes at exterior state

         FX3_EX = HUV
         FY3_EX = HVV + GH3

!.......Compute y momentum fluxes at exterior state

         FX3_EX = HUV*dg_here%SFAC_EX
         FY3_EX = HVV + GH3

!.......Set the fluxes
         
         F_HAT = FX1_EX*dg_here%NX + FY1_EX*dg_here%NY 
         G_HAT = FX2_EX*dg_here%NX + FY2_EX*dg_here%NY + 0.5D0*rVncU
         H_HAT = FX3_EX*dg_here%NX + FY3_EX*dg_here%NY + 0.5D0*rVncV

         do ll = 1,s%layers        !CAUTION, only really makes sense for single layer

               dg_here%bed_HAT(ll) = 0.5D0*(discharge_modelX_IN+discharge_modelX_EX)*dg_here%NX &
           + 0.5D0*(discharge_modelY_IN+discharge_modelY_EX)*dg_here%NY

         enddo

      endif
      
#ifdef TRACE
      U_AVG = 0.5D0 * ( U_IN + U_EX )
      V_AVG = 0.5D0 * ( V_IN + V_EX )
      VEL_NORMAL = U_AVG*dg_here%NX + V_AVG*dg_here%NY
      
      if (VEL_NORMAL.GT.0) THEN
         I_HAT = VEL_NORMAL*dg_here%iota_IN
      else
         I_HAT = VEL_NORMAL*dg_here%iota_EX
      endif
#endif

#ifdef CHEM
      U_AVG = 0.5D0 * ( U_IN + U_EX )
      V_AVG = 0.5D0 * ( V_IN + V_EX )
      VEL_NORMAL = U_AVG*dg_here%NX + V_AVG*dg_here%NY
      
      if (VEL_NORMAL.GT.0) THEN
         I_HAT = VEL_NORMAL*dg_here%iota_IN
         J_HAT = VEL_NORMAL*dg_here%iota2_IN
      else
         I_HAT = VEL_NORMAL*dg_here%iota_EX
         J_HAT = VEL_NORMAL*dg_here%iota2_EX
      endif
#endif

#ifdef DYNP
      U_AVG = 0.5D0 * ( U_IN + U_EX )
      V_AVG = 0.5D0 * ( V_IN + V_EX )
      VEL_NORMAL = U_AVG*dg_here%NX + V_AVG*dg_here%NY
      
      if (VEL_NORMAL.GT.0) THEN
         K_HAT = VEL_NORMAL*dynP_IN
      else
         K_HAT = VEL_NORMAL*dynP_EX
      endif
#endif


      RETURN
      END SUBROUTINE
      
      
