C***********************************************************************
C
C     SUBROUTINE PRECIPITATION( )
C
C     This subroutine reads the precipitation file and interpolates 
C     over the mesh
C
C     Current version generates only uniform precipitation
C
C     Written by Younghun Kang
C
C***********************************************************************

      MODULE PRECIPITATION
        use global, only : sz, NP, NM
        use DG
        implicit none

        REAL(SZ), ALLOCATABLE :: PREC2(:), PREC3(:) ! Parametric rainfall
        !real(sz) source_r

        contains

      SUBROUTINE alloc_precipitation()
        implicit none
        allocate (prec2(NP), prec3(NP))
        prec3 = 0.0
      END SUBROUTINE

      FUNCTION computeElemTRR(II) RESULT(SOURCE_R1)
        use global, only : sz
        implicit none
        INTEGER, INTENT(IN) :: II ! global index of element
        REAL(SZ) :: SOURCE_R1
        INTEGER :: N1, N2, N3

        N1 = NM(II,1)
        N2 = NM(II,2)
        N3 = NM(II,3)
        SOURCE_R1 =  1.0/3.0*(PREC2(N1) + PREC2(N2) + PREC2(N3))
        if (SOURCE_R1 .LT. 0.0) then
           SOURCE_R1 = 0.0
        ENDIF
      END FUNCTION

      function computeTRR_IPET(dist, LatestRmax, Pn, Pc) result(TRR)
        use global, only : sz
        implicit none
        real(sz) :: TRR
        real(sz), intent(in) :: dist, Pn, Pc
        real(sz), intent(in) :: LatestRmax

        IF (dist <= LatestRmax) THEN
            TRR = ((1.14) + (0.12 * (Pn - Pc))) 
        ELSEIF (dist>LatestRmax) THEN
        IF (dist<500) THEN
            TRR = 1.14 + 0.12 * (Pn - Pc) * (EXP(-0.3 * ((dist - LatestRmax) / LatestRmax))) 
        ELSEIF (dist>500) THEN
            TRR = 0
        END IF
        END IF

        ! Convert mm/hr to m/s
        TRR = TRR * 1e-3 / 3600.0
      END FUNCTION

      function computeTRR_RCLIPER(dist, Vmax) result(TRR)
        use global, only : sz

        implicit none
        real(sz), intent(in) :: dist  ! must be in km
        real(sz), intent(in) :: Vmax
        real(sz) :: TRR

        !...Constants for parametric rainfall model
        REAL(sz) :: a1,a2,a3,a4
        REAL(sz) :: b1,b2,b3,b4
        REAL(sz) :: NMW
        REAL(sz) :: T0
        REAL(sz) :: Tm
        REAL(sz) :: rm
        REAL(sz) :: re

        a1 = -1.10 !inches/day
        a2 = -1.6d0 !inches/day
        a3 =  64.5d0 !kilometers
        a4 =  150.0 !kilometers
        b1 =  3.96d0 !inches/day
        b2 =  4.80 !inches/day
        b3 = -13.0 !kilometers
        b4 = -16.0 !kilometers

        ! Equations 
        NMW=(1.0 + ((Vmax-35.0)/33.0)) !normalized maximum wind
        T0=a1+(b1*NMW) ! rain rate at TC center r=0
        Tm=a2+(b2*NMW) !maximum rain rate
        rm=a3+(b3*NMW) !radius from the center at which maximum rain rate occurs
        re=a4+(b4*NMW) ! curve fit parameter; specifies end behavior of rainfall rate curves
        ! Compute TRR
        IF (dist .GE. 500) THEN
        TRR = 0.0
        ELSEIF (dist .LT. 500) THEN
           IF (dist .LT. rm) THEN
           TRR = T0 + (Tm - T0) * (dist / rm) ! TRR is rain rate in inches/day
           ELSEIF (dist .GE. rm) THEN
           TRR = Tm * EXP(-((dist - rm) / re))
        ENDIF
        ENDIF

        !convert inches/day-> m/sec
        TRR=TRR/3401568.0
      
      end function  


      SUBROUTINE UPDATE_PREC(i, lat, lon, cLon, cLat, Pn, Pc, rmx, Vmax, model_type)
        use global, only : Rearth, deg2rad, DTDP, rad2deg
        use vortex
        use DG, only : IRK
        implicit none
        integer, intent(in) :: i, model_type
        real(sz), intent(in) :: lon, lat, cLon, cLat, Pn, Pc, rmx, Vmax
        real(sz) :: dx, dy, dist, LatestRmax 

        dx = deg2rad * Rearth * (lon - cLon) * COS(deg2rad*cLat)
        dy = deg2rad * Rearth * (lat - cLat)
        dist = SQRT(dx*dx + dy*dy)

        dist = dist / 1000.0  ! convert to km
        
        LatestRmax = rmx*1.852  ! Assign the latest value of rmx to LatestRmax

        
        SELECT CASE (model_type)
          CASE (1)  ! Use RCLIPER
          PREC2(i) = computeTRR_RCLIPER(dist, Vmax)
          CASE (3)  ! Use IPET
          PREC2(i) = computeTRR_IPET(dist,LatestRmax, Pn, Pc)
            END SELECT  
        
        IF (IRK == 1) THEN
            PREC3(i) = PREC3(i) + DTDP * PREC2(i)
        END IF
      END SUBROUTINE

           END MODULE
