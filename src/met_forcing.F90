
!***********************************************************************
!
!     SUBROUTINE MET_FORCING()
!
!     This subroutine handles the meteorological forcing for the DG
!     code
!
!***********************************************************************

      SUBROUTINE MET_FORCING(s,dg_here,global_here,nodalattr_here,IT)
      
      USE GLOBAL
      USE DG
#ifdef HARM
      USE HARM
#endif
      USE SIZES
      USE WIND
#ifdef OWIWIND
      USE OWIWIND,ONLY : NWS12INIT,NWS12GET
#endif
      USE NodalAttributes
#ifdef SWAN
!asey 101118: We need these values from other places.
      USE OWIWIND,     ONLY: WindMultiplier
      USE Couple2Swan, ONLY: COUPWIND,&
                       SWAN_WX2,&
                       SWAN_WY2
#endif
      
      IMPLICIT NONE

      type (sizes_type) :: s
      type (dg_type) :: dg_here
      type (global_type) :: global_here
      type (nodalattr_type) :: nodalattr_here

      REAL(SZ) WindDragLimit
      INTEGER II, IT

!.....Set the wind drag limit

      WindDragLimit = 0.002
      global_here%RampMete = dg_here%rampdg

!asey 130710: Added this section.
      IF(global_here%WTIME1.LT.global_here%ITHS*global_here%DTDP)THEN
         global_here%WTIME1 = global_here%ITHS*global_here%DTDP
         global_here%WTIME2 = global_here%WTIME1 + global_here%WTIMINC
      ENDIF
!-----------------------------------------------------------------------
!
!     global_here%NWS = 1
!
!     Wind stress and atmospheric pressure are read in at all grid nodes
!     at every model time step from the fort.22 file.
!
!-----------------------------------------------------------------------

      IF (global_here%NWS.EQ.1) THEN
         DO II= 1,global_here%NP
         
!..........Read in the data

           READ(s%fort22unit,*) global_here%NHG, global_here%WSX2(II), global_here%WSY2(II), global_here%PR2(II)

!..........Apply the met global_here%ramp function

!           global_here%RampMete = dg_here%RAMPDG
           
           global_here%WSX2(II)    = global_here%RampMete*global_here%WSX2(II)
           global_here%WSY2(II)    = global_here%RampMete*global_here%WSY2(II)
           global_here%PR2(II)     = global_here%RampMete*global_here%PR2(II)
           global_here%WVNXOUT(II) = global_here%WSX2(II)
           global_here%WVNYOUT(II) = global_here%WSY2(II)
         ENDDO
      ENDIF
      
!-----------------------------------------------------------------------
!
!     global_here%NWS = 2
!
!     Wind stress and atmospheric pressure are read in at all grid nodes
!     at a time interval that does not equal the model time step.  In-
!     terpolation in time is used to synchronize the wind and pressure
!     information with the model time step.
!
!-----------------------------------------------------------------------

      IF (ABS(global_here%NWS).EQ.2) THEN
      
!.......Determine if the met file time increment is exceeded
      
        IF (global_here%TIME_A.GT.global_here%WTIME2) THEN
          global_here%WTIME1 = global_here%WTIME2
          global_here%WTIME2 = global_here%WTIME2 + global_here%WTIMINC
          DO II= 1,global_here%NP
          
!...........Shift current data to old
          
            global_here%WVNX1(II) = global_here%WVNX2(II)
            global_here%WVNY1(II) = global_here%WVNY2(II)
            global_here%PRN1(II)  = global_here%PRN2(II)
            
!...........Read in data

            READ(s%fort22unit,*) global_here%NHG, global_here%WVNX2(II), global_here%WVNY2(II), global_here%PRN2(II)
          ENDDO
          PRINT*,'READING IN WIND DATA SET AT TIMESTEP',IT
        ENDIF
        
        global_here%WTRATIO = (global_here%TIME_A - global_here%WTIME1)/global_here%WTIMINC
        DO II= 1,global_here%NP
        
!.........Interpolate in time
        
          global_here%WINDX      = global_here%WVNX1(II) + global_here%WTRATIO*(global_here%WVNX2(II) - global_here%WVNX1(II))
          global_here%WINDY      = global_here%WVNY1(II) + global_here%WTRATIO*(global_here%WVNY2(II) - global_here%WVNY1(II))
          
!.........Apply mete global_here%ramp

!          global_here%RampMete = dg_here%RAMPDG
          
          global_here%WSX2(II)    = global_here%RampMete*global_here%WINDX
          global_here%WSY2(II)    = global_here%RampMete*global_here%WINDY
          global_here%PR2(II)     = global_here%RampMete*(global_here%PRN1(II)+global_here%WTRATIO*(global_here%PRN2(II)-global_here%PRN1(II)))
          global_here%WVNXOUT(II) = global_here%WSX2(II)
          global_here%WVNYOUT(II) = global_here%WSY2(II)
        ENDDO
      ENDIF
      
!-----------------------------------------------------------------------
!
!     global_here%NWS = 3
!
!     Wind velocity in US Navy Fleet Numeric format interpolated in
!     space onto the ADCIRC grid and in time to synchronize the wind and
!     pressure information with the model time step. Garratt's formula
!     is used to compute wind stress from the wind velocity.
!
!-----------------------------------------------------------------------

      IF (global_here%NWS.EQ.3) THEN
      
!.......Determine if the met file time increment is exceeded

        IF (global_here%TIME_A.GT.global_here%WTIME2) THEN
          global_here%WTIME1 = global_here%WTIME2
          global_here%WTIME2 = global_here%WTIME2 + global_here%WTIMINC
          
!.........Shift current data to old
          
          DO II=1,global_here%NP
            global_here%WVNX1(II) = global_here%WVNX2(II)
            global_here%WVNY1(II) = global_here%WVNY2(II)
          ENDDO
          
!.........Obtain the meteorological forcing data
          
          CALL NWS3GET(s, global_here%X, global_here%Y, global_here%SLAM, global_here%SFEA, global_here%WVNX2, global_here%WVNY2, global_here%IWTIME, global_here%IWYR,&
                  global_here%WTIMED, global_here%NP, global_here%NWLON, global_here%NWLAT, global_here%WLATMAX, global_here%WLONMIN,&
                  global_here%WLATINC, global_here%WLONINC, global_here%ICS, global_here%NSCREEN, global_here%ScreenUnit )
         ENDIF

         global_here%WTRATIO = (global_here%TIME_A - global_here%WTIME1)/global_here%WTIMINC
         DO II= 1,global_here%NP
         
!..........Interpolate in time
         
           global_here%WINDX   = global_here%WVNX1(II) + global_here%WTRATIO*(global_here%WVNX2(II) - global_here%WVNX1(II))
           global_here%WINDY   = global_here%WVNY1(II) + global_here%WTRATIO*(global_here%WVNY2(II) - global_here%WVNY1(II))
           
!..........Compute wind drag
           
           global_here%WINDMAG = SQRT( global_here%WINDX*global_here%WINDX + global_here%WINDY*global_here%WINDY )
           global_here%WDRAGCO = WindDrag( global_here%WINDMAG, WindDragLimit, "Garratt   ")
           
!..........Apply directional wind reductions
           
           IF (nodalattr_here%LoadDirEffRLen) THEN
             CALL ApplyDirectionalWindReduction(nodalattr_here, II, global_here%WDRAGCO, global_here%WINDMAG,&
                                          global_here%DP(II), global_here%ETA2(II), global_here%H0, global_here%G,&
                                          global_here%WINDX, global_here%WINDY )
             global_here%WINDMAG = SQRT( global_here%WINDX*global_here%WINDX + global_here%WINDY*global_here%WINDY )
             global_here%WDRAGCO = WindDrag(global_here%WINDMAG, WindDragLimit, "Garratt   " )
           ENDIF
           
!..........Apply met global_here%ramp
!           global_here%RampMete = dg_here%RAMPDG
           global_here%WSX2(II)    = global_here%RampMete*0.001293D0*global_here%WDRAGCO*global_here%WINDX*global_here%WINDMAG
           global_here%WSY2(II)    = global_here%RampMete*0.001293D0*global_here%WDRAGCO*global_here%WINDY*global_here%WINDMAG
           global_here%WVNXOUT(II) = global_here%RampMete*global_here%WINDX
           global_here%WVNYOUT(II) = global_here%RampMete*global_here%WINDY
#ifdef SWAN
!asey 101118: Added these lines for coupling winds to SWAN.
           IF(COUPWIND)THEN
             SWAN_WX2(II,2) = global_here%WINDX
             SWAN_WY2(II,2) = global_here%WINDY
           ENDIF
#endif
        ENDDO
      ENDIF
      
!-----------------------------------------------------------------------
!
!     global_here%NWS = 4
!
!     Wind velocity and atmospheric pressure are read in (PBL/JAG
!     format) at selected ADCIRC grid nodes. Interpolation in time is
!     used to synchronize the wind and pressure information with the
!     model time step. Garratt's formula is used to compute wind stress
!     from wind velocity.
!
!-----------------------------------------------------------------------

      IF (ABS(global_here%NWS).EQ.4) THEN
      
!.......Determine if the met file time increment is exceeded
      
        IF (global_here%TIME_A.GT.global_here%WTIME2) THEN
          global_here%WTIME1 = global_here%WTIME2
          global_here%WTIME2 = global_here%WTIME2 + global_here%WTIMINC

!.........Shift current data to old
          
          DO II = 1,global_here%NP
            global_here%WVNX1(II) = global_here%WVNX2(II)
            global_here%WVNY1(II) = global_here%WVNY2(II)
            global_here%PRN1(II)  = global_here%PRN2(II)
          ENDDO
          
!.........Obtain the meteorological forcing data
          
          CALL NWS4GET( global_here%WVNX2, global_here%WVNY2, global_here%PRN2, global_here%NP, global_here%RHOWAT0, global_here%G )
        ENDIF

        global_here%WTRATIO = (global_here%TIME_A-global_here%WTIME1)/global_here%WTIMINC
        DO II = 1,global_here%NP
         
!.........Interpolate in time
         
          global_here%WINDX = global_here%WVNX1(II) + global_here%WTRATIO*(global_here%WVNX2(II)-global_here%WVNX1(II))
          global_here%WINDY = global_here%WVNY1(II) + global_here%WTRATIO*(global_here%WVNY2(II)-global_here%WVNY1(II))
            
!.........Compute wind drag
            
          global_here%WINDMAG = SQRT( global_here%WINDX*global_here%WINDX + global_here%WINDY*global_here%WINDY )
          global_here%WDRAGCO = WindDrag( global_here%WINDMAG, WindDragLimit, "Garratt   " )
           
!.........Apply directional wind reductions

          IF (nodalattr_here%LoadDirEffRLen) THEN
            CALL ApplyDirectionalWindReduction(nodalattr_here, II, global_here%WDRAGCO, global_here%WINDMAG,&
                                          global_here%DP(II), global_here%ETA2(II), global_here%H0, global_here%G,&
                                          global_here%WINDX, global_here%WINDY )
            global_here%WINDMAG = SQRT(global_here%WINDX*global_here%WINDX+global_here%WINDY*global_here%WINDY)
            global_here%WDRAGCO = WindDrag(global_here%WINDMAG, WindDragLimit, "Garratt   ")
          ENDIF
            
!.........Apply met global_here%ramp
!           global_here%RampMete = dg_here%RAMPDG
          global_here%WSX2(II)    = global_here%RampMete*0.001293d0*global_here%WDRAGCO*global_here%WINDX*global_here%WINDMAG
          global_here%WSY2(II)    = global_here%RampMete*0.001293d0*global_here%WDRAGCO*global_here%WINDY*global_here%WINDMAG
          global_here%PR2(II)     = global_here%RampMete*(global_here%PRN1(II)+global_here%WTRATIO*(global_here%PRN2(II)- global_here%PRN1(II)))
          global_here%WVNXOUT(II) = global_here%RampMete*global_here%WINDX
          global_here%WVNYOUT(II) = global_here%RampMete*global_here%WINDY
#ifdef SWAN
!asey 101118: Added these lines for coupling winds to SWAN.
          IF(COUPWIND)THEN
            SWAN_WX2(II,2) = global_here%WINDX
            SWAN_WY2(II,2) = global_here%WINDY
          ENDIF
#endif
        ENDDO
      ENDIF
      
!-----------------------------------------------------------------------
!
!     global_here%NWS = 5
!
!     Wind velocity and atmospheric pressure are read in at all grid
!     nodes. Interpolation in time is used to synchronize the wind and
!     pressure information with the model time step. Garratt's formula
!     is used to compute wind stress from wind velocity.
!
!-----------------------------------------------------------------------

      IF(ABS(global_here%NWS).EQ.5) THEN
      
!.......Determine if the met file time increment is exceeded
      
        IF(global_here%TIME_A.GT.global_here%WTIME2) THEN
          global_here%WTIME1 = global_here%WTIME2
          global_here%WTIME2 = global_here%WTIME2 + global_here%WTIMINC
          
!.........Shift current data to old
          
          DO II = 1,global_here%NP
            global_here%WVNX1(II) = global_here%WVNX2(II)
            global_here%WVNY1(II) = global_here%WVNY2(II)
            global_here%PRN1(II)  = global_here%PRN2(II)
            
!...........Read in the meteorological forcing data

            READ(s%fort22unit,*) global_here%NHG, global_here%WVNX2(II), global_here%WVNY2(II), global_here%PRN2(II)
          ENDDO
        ENDIF
        
        global_here%WTRATIO = (global_here%TIME_A - global_here%WTIME1)/global_here%WTIMINC
        DO II = 1,global_here%NP
        
!.........Interpolate in time
        
          global_here%WINDX   = global_here%WVNX1(II) + global_here%WTRATIO*(global_here%WVNX2(II) - global_here%WVNX1(II))
          global_here%WINDY   = global_here%WVNY1(II) + global_here%WTRATIO*(global_here%WVNY2(II) - global_here%WVNY1(II))
          
!.........Compute wind drag
          
          global_here%WINDMAG = SQRT( global_here%WINDX*global_here%WINDX + global_here%WINDY*global_here%WINDY )
          global_here%WDRAGCO = WindDrag( global_here%WINDMAG, WindDragLimit, "Garratt   " )
          
!.........Apply directional wind reductions
          
          IF (nodalattr_here%LoadDirEffRLen) THEN
            CALL ApplyDirectionalWindReduction(nodalattr_here, II, global_here%WDRAGCO, global_here%WINDMAG,&
                                          global_here%DP(II), global_here%ETA2(II), global_here%H0, global_here%G,&
                                          global_here%WINDX, global_here%WINDY )
            global_here%WINDMAG = SQRT( global_here%WINDX*global_here%WINDX + global_here%WINDY*global_here%WINDY )
            global_here%WDRAGCO = WindDrag(global_here%WINDMAG, WindDragLimit, "Garratt   ")
          ENDIF
          
!.........Apply met global_here%ramp
          
          global_here%WSX2(II)    = global_here%RampMete*0.001293d0*global_here%WDRAGCO*global_here%WINDX*global_here%WINDMAG
          global_here%WSY2(II)    = global_here%RampMete*0.001293d0*global_here%WDRAGCO*global_here%WINDY*global_here%WINDMAG
          global_here%PR2(II)     = global_here%RampMete*(global_here%PRN1(II)+global_here%WTRATIO*(global_here%PRN2(II)-global_here%PRN1(II)))
          global_here%WVNXOUT(II) = global_here%RampMete*global_here%WINDX
          global_here%WVNYOUT(II) = global_here%RampMete*global_here%WINDY
#ifdef SWAN
!asey 101118: Added these lines for coupling winds to SWAN.
          IF(COUPWIND)THEN
            SWAN_WX2(II,2) = global_here%WINDX
            SWAN_WY2(II,2) = global_here%WINDY
          ENDIF
#endif
        ENDDO
      ENDIF
      
!-----------------------------------------------------------------------
!
!     global_here%NWS = 6
!
!     Wind velocity and atmospheric pressure are read in for a
!     rectangular grid (either in Longitude, Latitude or Cartesian
!     coordinates, consistent with the grid coordinates) and
!     interpolated in space onto the ADCIRC grid and in time to
!     synchronize the wind and pressure information with the model time
!     step. Garratt's formula is used to compute wind stress from the
!     wind velocity.
!
!-----------------------------------------------------------------------

      IF (global_here%NWS.EQ.6) THEN
      
!.......Determine if the met file time increment is exceeded

        IF (global_here%TIME_A.GT.global_here%WTIME2) THEN
          global_here%WTIME1 = global_here%WTIME2
          global_here%WTIME2 = global_here%WTIME2 + global_here%WTIMINC
          
!.........Shift current data to old
          
          DO II= 1,global_here%NP
            global_here%WVNX1(II) = global_here%WVNX2(II)
            global_here%WVNY1(II) = global_here%WVNY2(II)
            global_here%PRN1(II)  = global_here%PRN2(II)
          ENDDO
          global_here%NWSGGWI = global_here%NWSGGWI + 1
          
!.........Obtain meteorological forcing data
          
          CALL NWS6GET( global_here%X, global_here%Y, global_here%SLAM, global_here%SFEA, global_here%WVNX2, global_here%WVNY2, global_here%PRN2, global_here%NP,&
                  global_here%NWLON, global_here%NWLAT, global_here%WLATMAX, global_here%WLONMIN, global_here%WLATINC,&
                  global_here%WLONINC, global_here%ICS, global_here%RHOWAT0, global_here%G )
        ENDIF
        
        global_here%WTRATIO=(global_here%TIME_A-global_here%WTIME1)/global_here%WTIMINC
        DO II= 1,global_here%NP
        
!.........Interpolate in time
        
          global_here%WINDX = global_here%WVNX1(II) + global_here%WTRATIO*(global_here%WVNX2(II)-global_here%WVNX1(II))
          global_here%WINDY = global_here%WVNY1(II) + global_here%WTRATIO*(global_here%WVNY2(II)-global_here%WVNY1(II))
          
!.........Compute wind drag
          
          global_here%WINDMAG = SQRT(global_here%WINDX*global_here%WINDX+global_here%WINDY*global_here%WINDY)
          global_here%WDRAGCO = WindDrag( global_here%WINDMAG, WindDragLimit, "Garratt   " )
          
!.........Apply directional wind reductions

          IF (nodalattr_here%LoadDirEffRLen) THEN
            CALL ApplyDirectionalWindReduction(nodalattr_here, II, global_here%WDRAGCO, global_here%WINDMAG,&
                                          global_here%DP(II), global_here%ETA2(II), global_here%H0, global_here%G,&
                                          global_here%WINDX, global_here%WINDY )
            global_here%WINDMAG = SQRT(global_here%WINDX*global_here%WINDX+global_here%WINDY*global_here%WINDY)
            global_here%WDRAGCO = WindDrag(global_here%WINDMAG, WindDragLimit, "Garratt   ")
          ENDIF
          
!.........Apply met global_here%ramp
          
          global_here%WSX2(II)    = global_here%RampMete*0.001293d0*global_here%WDRAGCO*global_here%WINDX*global_here%WINDMAG
          global_here%WSY2(II)    = global_here%RampMete*0.001293d0*global_here%WDRAGCO*global_here%WINDY*global_here%WINDMAG
          global_here%PR2(II)     = global_here%RampMete*(global_here%PRN1(II)+global_here%WTRATIO*(global_here%PRN2(II)-global_here%PRN1(II)))
          global_here%WVNXOUT(II) = global_here%RampMete*global_here%WINDX
          global_here%WVNYOUT(II) = global_here%RampMete*global_here%WINDY
#ifdef SWAN
!asey 101118: Added these lines for coupling winds to SWAN.
          IF(COUPWIND)THEN
            SWAN_WX2(II,2) = global_here%WINDX
            SWAN_WY2(II,2) = global_here%WINDY
          ENDIF
#endif
        ENDDO
      ENDIF
      
!-----------------------------------------------------------------------
!
!     global_here%NWS = 7
!
!     jgf46.01 New option to read in surface wind stress and atmospheric
!     pressure for a rectangular grid (either in Longitude, Latitude or
!     Cartesian coordinates, consistent with the grid coordinates) and
!     interpolate in space onto the ADCIRC grid. Interpolation in time
!     is used to synchronize the wind and pressure information with the
!     model time step.
!
!-----------------------------------------------------------------------

      IF(ABS(global_here%NWS).EQ.7) THEN
      
!.......Determine if the met file time increment is exceeded
      
        IF (global_here%TIME_A.GT.global_here%WTIME2) THEN
          global_here%WTIME1 = global_here%WTIME2
          global_here%WTIME2 = global_here%WTIME2 + global_here%WTIMINC
          
!.........Shift current data to old
          
          DO II= 1,global_here%NP
            global_here%WVNX1(II) = global_here%WVNX2(II)
            global_here%WVNY1(II) = global_here%WVNY2(II)
            global_here%PRN1(II)  = global_here%PRN2(II)
          ENDDO
          
!.........Obtain the meteorological forcing data
          
          CALL NWS7GET( global_here%X, global_here%Y, global_here%SLAM, global_here%SFEA, global_here%WVNX2, global_here%WVNY2, global_here%PRN2, global_here%NP, global_here%NWLON,&
                  global_here%NWLAT, global_here%WLATMAX, global_here%WLONMIN, global_here%WLATINC, global_here%WLONINC, global_here%ICS,&
                  global_here%RHOWAT0,global_here%G )
        ENDIF

        global_here%WTRATIO=(global_here%TIME_A-global_here%WTIME1)/global_here%WTIMINC
        DO II= 1,global_here%NP
        
!.........Interpolate in time
        
          global_here%WINDX = global_here%WVNX1(II) + global_here%WTRATIO*(global_here%WVNX2(II) - global_here%WVNX1(II))
          global_here%WINDY = global_here%WVNY1(II) + global_here%WTRATIO*(global_here%WVNY2(II) - global_here%WVNY1(II))
          
!.........Apply met global_here%ramp
          
          global_here%WSX2(II)    = global_here%RampMete*global_here%WINDX
          global_here%WSY2(II)    = global_here%RampMete*global_here%WINDY
          global_here%PR2(II)     = global_here%RampMete*(global_here%PRN1(II)+global_here%WTRATIO*(global_here%PRN2(II)-global_here%PRN1(II)))
          global_here%WVNXOUT(II) = global_here%WSX2(II)
          global_here%WVNYOUT(II) = global_here%WSY2(II)
        ENDDO
      ENDIF
      
!-----------------------------------------------------------------------
!
!     global_here%NWS = 8
!
!     jgf46.02 New option to read in hurricane locations and generate
!     hurricane winds from the Holland Wind Model.
!
!-----------------------------------------------------------------------

      IF (ABS(global_here%NWS).EQ.8) THEN
      
!.......Obtain the meteorological forcing data
!         write(*,*) 'calling HollandGet ',global_here%time_a

         CALL HollandGet(s, global_here%X, global_here%Y, global_here%SLAM, global_here%SFEA, global_here%WVNX2, global_here%WVNY2, global_here%PRN2, global_here%NP, global_here%ICS,&
              global_here%RHOWAT0, global_here%G, global_here%TIME_A, global_here%NSCREEN, global_here%ScreenUnit )
         DO II= 1,global_here%NP
            global_here%WINDX = global_here%WVNX2(II)
            global_here%WINDY = global_here%WVNY2(II)
            
            !.........Compute wind drag
            
            global_here%WINDMAG = SQRT(global_here%WINDX*global_here%WINDX+global_here%WINDY*global_here%WINDY)
            global_here%WDRAGCO = WindDrag( global_here%WINDMAG, WindDragLimit, "Garratt   " )
            
            !.........Apply directional wind reductions
            
            IF (nodalattr_here%LoadDirEffRLen) THEN
               CALL ApplyDirectionalWindReduction(nodalattr_here, II, global_here%WDRAGCO, global_here%WINDMAG,&
                    global_here%DP(II), global_here%ETA2(II), global_here%H0, global_here%G,&
                    global_here%WINDX, global_here%WINDY )
               global_here%WINDMAG = SQRT(global_here%WINDX*global_here%WINDX + global_here%WINDY*global_here%WINDY)
               global_here%WDRAGCO = WindDrag( global_here%WINDMAG, WindDragLimit, "Garratt   " )
            ENDIF
            
            !.........Apply met global_here%ramp
            
            global_here%WSX2(II)    = global_here%RampMete*0.001293d0*global_here%WDRAGCO*global_here%WINDX*global_here%WINDMAG
            global_here%WSY2(II)    = global_here%RampMete*0.001293d0*global_here%WDRAGCO*global_here%WINDY*global_here%WINDMAG
            global_here%PR2(II)     = global_here%RampMete*global_here%PRN2(II)
            global_here%WVNXOUT(II) = global_here%RampMete*global_here%WINDX
            global_here%WVNYOUT(II) = global_here%RampMete*global_here%WINDY
         ENDDO
      ENDIF
      
!-----------------------------------------------------------------------
!
!     global_here%NWS = 9
!
!     jgf46.16 Merged:
!     nodalattr_here%cf & cm added global_here%nws = 9: asymmetric hurricane winds
!
!-----------------------------------------------------------------------

!      IF (global_here%NWS.EQ.9) THEN
      
!.......Obtain meteorological forcing data
      
!        CALL NWS9GET(global_here%SLAM,global_here%SFEA,global_here%WVNX2,global_here%WVNY2,global_here%PRN2,global_here%NP,global_here%TIME_A, global_here%ICS)
!        DO II= 1,global_here%NP
!          global_here%WINDX = global_here%WVNX2(II)
!          global_here%WINDY = global_here%WVNY2(II)
          
!.........Compute wind drag
          
!          global_here%WINDMAG = SQRT(global_here%WINDX*global_here%WINDX + global_here%WINDY*global_here%WINDY)
!          global_here%WDRAGCO = WindDrag( global_here%WINDMAG, WindDragLimit, "Garratt   " )
          
!.........Apply directional wind reductions
          
!          IF (nodalattr_here%LoadDirEffRLen) THEN
!            CALL ApplyDirectionalWindReduction( II, global_here%WDRAGCO, global_here%WINDMAG,
!     &                                          global_here%DP(II), global_here%ETA2(II), global_here%H0, global_here%G,
!     &                                          global_here%WINDX, global_here%WINDY )
!            global_here%WINDMAG = SQRT(global_here%WINDX*global_here%WINDX+global_here%WINDY*global_here%WINDY)
!            global_here%WDRAGCO = WindDrag(global_here%WINDMAG, WindDragLimit, "Garratt   ")
!          ENDIF
          
!.........Apply met global_here%ramp

!          global_here%WSX2(II)    = global_here%RampMete*0.001293d0*global_here%WDRAGCO*global_here%WINDX*global_here%WINDMAG
!          global_here%WSY2(II)    = global_here%RampMete*0.001293d0*global_here%WDRAGCO*global_here%WINDY*global_here%WINDMAG
!          global_here%PR2(II)     = global_here%RampMete*global_here%PRN2(II)
!          global_here%WVNXOUT(II) = global_here%RampMete*global_here%WINDX
!          global_here%WVNYOUT(II) = global_here%RampMete*global_here%WINDY
#ifdef SWAN
!asey 101118: Added these lines for coupling winds to SWAN.
!          IF(COUPWIND)THEN
!            SWAN_WX2(II,2) = global_here%WINDX
!            SWAN_WY2(II,2) = global_here%WINDY
!          ENDIF
#endif
!         ENDDO
!      ENDIF
      
!-----------------------------------------------------------------------
!
!     global_here%NWS = 10
!
!     Wind velocity (10 m) and atmospheric pressure are read in from a
!     sequence of National Weather Service (global_here%NWS) Aviation (AVN) model
!     output files. Each AVN file is assumed to contain data on a
!     Gaussian longitude, latitude grid at a single time. Consecutive
!     files in the sequence are separated by N hours in time. Garratt's
!     formula is used to compute wind stress from the wind velocity.
!
!-----------------------------------------------------------------------

      IF (global_here%NWS.EQ.10) THEN
      
!.......Determine if the met file time increment is exceeded
      
        IF (global_here%TIME_A.GT.global_here%WTIME2) THEN
          global_here%WTIME1 = global_here%WTIME2
          global_here%WTIME2 = global_here%WTIME2 + global_here%WTIMINC
          
!.........Shift current data to old
          
          DO II= 1,global_here%NP
            global_here%WVNX1(II) = global_here%WVNX2(II)
            global_here%WVNY1(II) = global_here%WVNY2(II)
            global_here%PRN1(II)  = global_here%PRN2(II)
          ENDDO
          global_here%NWSGGWI = global_here%NWSGGWI + 1
          
!.........Obtain meteorological forcing data
          
          CALL NWS10GET(s,global_here, global_here%NWSGGWI, global_here%SLAM, global_here%SFEA, global_here%WVNX2, global_here%WVNY2, global_here%PRN2, global_here%NP,&
                   global_here%RHOWAT0, global_here%G, global_here%NWLON, global_here%NWLAT, global_here%WTIMINC )
        ENDIF
        
        global_here%WTRATIO = (global_here%TIME_A - global_here%WTIME1)/global_here%WTIMINC
        DO II = 1,global_here%NP
        
!.........Interpolate in time
        
          global_here%WINDX = global_here%WVNX1(II) + global_here%WTRATIO*(global_here%WVNX2(II)-global_here%WVNX1(II))
          global_here%WINDY = global_here%WVNY1(II) + global_here%WTRATIO*(global_here%WVNY2(II)-global_here%WVNY1(II))
           
!.........Compute wind drag
           
          global_here%WINDMAG = SQRT(global_here%WINDX*global_here%WINDX + global_here%WINDY*global_here%WINDY)
          global_here%WDRAGCO = WindDrag( global_here%WINDMAG, WindDragLimit, "Garratt   " )

!.........Apply directional wind reductions

          IF (nodalattr_here%LoadDirEffRLen) THEN
            CALL ApplyDirectionalWindReduction(nodalattr_here, II, global_here%WDRAGCO, global_here%WINDMAG,&
                                          global_here%DP(II), global_here%ETA2(II), global_here%H0, global_here%G,&
                                          global_here%WINDX, global_here%WINDY )
            global_here%WINDMAG = SQRT(global_here%WINDX*global_here%WINDX + global_here%WINDY*global_here%WINDY)
            global_here%WDRAGCO = WindDrag( global_here%WINDMAG, WindDragLimit, "Garratt   " )
          ENDIF
          
!.........Apply met global_here%ramp
          
          global_here%WSX2(II)    = global_here%RampMete*0.001293d0*global_here%WDRAGCO*global_here%WINDX*global_here%WINDMAG
          global_here%WSY2(II)    = global_here%RampMete*0.001293d0*global_here%WDRAGCO*global_here%WINDY*global_here%WINDMAG
          global_here%PR2(II)     = global_here%RampMete*(global_here%PRN1(II)+global_here%WTRATIO*(global_here%PRN2(II)-global_here%PRN1(II)))
          global_here%WVNXOUT(II) = global_here%RampMete*global_here%WINDX
          global_here%WVNYOUT(II) = global_here%RampMete*global_here%WINDY
#ifdef SWAN
!asey 101118: Added these lines for coupling winds to SWAN.
          IF(COUPWIND)THEN
            SWAN_WX2(II,2) = global_here%WINDX
            SWAN_WY2(II,2) = global_here%WINDY
          ENDIF
#endif
        ENDDO
      ENDIF
      
!-----------------------------------------------------------------------
!
!     global_here%NWS = 11
!
!     Wind velocity (10 m) and atmospheric pressure are read in from a
!     sequence of stripped down (?) National Weather Service (global_here%NWS) ETA
!     29km model output files. Each ETA file is assumed to contain data
!     on an E grid for a single global_here%day (8 data sets, one every 3 hours, be-
!     ginning @ 03:00 and continuing through 24:00 of the given global_here%day).
!     The wind data is converted to an east-west, north-south coordinate
!     system inside ADCIRC. Garratt's formula is used to compute wind
!     stress from the wind velocity.
!
!-----------------------------------------------------------------------

      IF(global_here%NWS.EQ.11) THEN

!.......Determine if the met file time increment is exceeded

        IF (global_here%TIME_A.GT.global_here%WTIME2) THEN
          global_here%WTIME1=global_here%WTIME2
          global_here%WTIME2=global_here%WTIME2+global_here%WTIMINC
           
!........Shift current data to old
           
          DO II = 1,global_here%NP
            global_here%WVNX1(II) = global_here%WVNX2(II)
            global_here%WVNY1(II) = global_here%WVNY2(II)
            global_here%PRN1(II)  = global_here%PRN2(II)
          ENDDO
          global_here%IDSETFLG = global_here%IDSETFLG + 1
          IF (global_here%IDSETFLG.GT.8) THEN
            global_here%NWSEGWI = global_here%NWSEGWI + 1
            global_here%IDSETFLG = 1
          ENDIF
            
!.........Obtain meteorological forcing data
            
          CALL NWS11GET(s,global_here, global_here%NWSEGWI, global_here%IDSETFLG, global_here%SLAM, global_here%SFEA, global_here%WVNX2, global_here%WVNY2,&
                   global_here%PRN2, global_here%NP, global_here%RHOWAT0, global_here%G )
        ENDIF

        global_here%WTRATIO=(global_here%TIME_A-global_here%WTIME1)/global_here%WTIMINC
        DO II = 1,global_here%NP
         
!.........Interpolate in time
         
          global_here%WINDX = global_here%WVNX1(II) + global_here%WTRATIO*(global_here%WVNX2(II)-global_here%WVNX1(II))
          global_here%WINDY = global_here%WVNY1(II) + global_here%WTRATIO*(global_here%WVNY2(II)-global_here%WVNY1(II))
            
!.........Compute wind drag
            
          global_here%WINDMAG = SQRT(global_here%WINDX*global_here%WINDX+global_here%WINDY*global_here%WINDY)
          global_here%WDRAGCO = WindDrag(global_here%WINDMAG, WindDragLimit, "Garratt   ")
            
!.........Apply directional wind reductions

          IF (nodalattr_here%LoadDirEffRLen) THEN
            CALL ApplyDirectionalWindReduction(nodalattr_here, II, global_here%WDRAGCO, global_here%WINDMAG,&
                                          global_here%DP(II), global_here%ETA2(II), global_here%H0, global_here%G,&
                                          global_here%WINDX, global_here%WINDY )
            global_here%WINDMAG = SQRT(global_here%WINDX*global_here%WINDX + global_here%WINDY*global_here%WINDY)
            global_here%WDRAGCO = WindDrag( global_here%WINDMAG, WindDragLimit, "Garratt   " )
          ENDIF
            
!.........Apply met global_here%ramp
            
          global_here%WSX2(II)    = global_here%RampMete*0.001293d0*global_here%WDRAGCO*global_here%WINDX*global_here%WINDMAG
          global_here%WSY2(II)    = global_here%RampMete*0.001293d0*global_here%WDRAGCO*global_here%WINDY*global_here%WINDMAG
          global_here%PR2(II)     = global_here%RampMete*(global_here%PRN1(II)+global_here%WTRATIO*(global_here%PRN2(II)-global_here%PRN1(II)))
          global_here%WVNXOUT(II) = global_here%RampMete*global_here%WINDX
          global_here%WVNYOUT(II) = global_here%RampMete*global_here%WINDY
#ifdef SWAN
!asey 101118: Added these lines for coupling winds to SWAN.
          IF(COUPWIND)THEN
            SWAN_WX2(II,2) = global_here%WINDX
            SWAN_WY2(II,2) = global_here%WINDY
          ENDIF
#endif
        ENDDO
      ENDIF
      
!-----------------------------------------------------------------------
!
!     global_here%NWS = 12
!
!     sb46.28sb01 global_here%NWS=12 reads in raw OWI files 09/xx/2006
!
!-----------------------------------------------------------------------
#ifdef OWIWIND
      IF(ABS(global_here%NWS).EQ.12) THEN
      
!.......Determine if the met file time increment is exceeded
      
        IF(global_here%TIME_A.GT.global_here%WTIME2) THEN
          global_here%WTIME1=global_here%WTIME2
          global_here%WTIME2=global_here%WTIME2+global_here%WTIMINC
          
!........Shift current data to old

          DO II =1,global_here%NP
            global_here%WVNX1(II) = global_here%WVNX2(II)
            global_here%WVNY1(II) = global_here%WVNY2(II)
            global_here%PRN1(II)  = global_here%PRN2(II)
          ENDDO
          
!.........Obtain meteorological forcing data
          
          CALL NWS12GET( global_here%WVNX2, global_here%WVNY2, global_here%PRN2, global_here%NP, global_here%RHOWAT0, global_here%G )
        ENDIF

        global_here%WTRATIO=(global_here%TIME_A - global_here%WTIME1)/global_here%WTIMINC
        DO II = 1,global_here%NP
        
!.........Interpolate in time
        
          global_here%WINDX = global_here%WVNX1(II) + global_here%WTRATIO*(global_here%WVNX2(II)-global_here%WVNX1(II))
          global_here%WINDY = global_here%WVNY1(II) + global_here%WTRATIO*(global_here%WVNY2(II)-global_here%WVNY1(II))
          
!.........Compute wind drag
          
          global_here%WINDMAG = SQRT(global_here%WINDX*global_here%WINDX + global_here%WINDY*global_here%WINDY)
          global_here%WDRAGCO = WindDrag( global_here%WINDMAG, WindDragLimit, "Garratt   " )
          
!.........Apply directional wind reductions
          
          IF (nodalattr_here%LoadDirEffRLen) THEN
            CALL ApplyDirectionalWindReduction(nodalattr_here, II, global_here%WDRAGCO, global_here%WINDMAG,&
                                          global_here%DP(II), global_here%ETA2(II), global_here%H0, global_here%G,&
                                          global_here%WINDX, global_here%WINDY )
            global_here%WINDMAG = SQRT(global_here%WINDX*global_here%WINDX + global_here%WINDY*global_here%WINDY)
            global_here%WDRAGCO = WindDrag(global_here%WINDMAG, WindDragLimit, "Garratt   ")
          ENDIF
          
!.........Apply met global_here%ramp
          
          global_here%WSX2(II)    = global_here%RampMete*0.001293d0*global_here%WDRAGCO*global_here%WINDX*global_here%WINDMAG
          global_here%WSY2(II)    = global_here%RampMete*0.001293d0*global_here%WDRAGCO*global_here%WINDY*global_here%WINDMAG
          global_here%PR2(II)     = global_here%RampMete*(global_here%PRN1(II)+global_here%WTRATIO*(global_here%PRN2(II)-global_here%PRN1(II)))
          global_here%WVNXOUT(II) = global_here%RampMete*global_here%WINDX
          global_here%WVNYOUT(II) = global_here%RampMete*global_here%WINDY
#ifdef SWAN
!asey 101118: Added these lines for coupling winds to SWAN.
          IF(COUPWIND)THEN
            SWAN_WX2(II,2) = global_here%WINDX/WindMultiplier
            SWAN_WY2(II,2) = global_here%WINDY/WindMultiplier
          ENDIF
#endif
        ENDDO
      ENDIF
#endif
      
      RETURN
      END SUBROUTINE MET_FORCING
