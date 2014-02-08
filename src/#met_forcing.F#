C***********************************************************************
C
C     SUBROUTINE MET_FORCING()
C
C     This subroutine handles the meteorological forcing for the DG
C     code
C
C***********************************************************************

      SUBROUTINE MET_FORCING(IT)
      
      USE GLOBAL
      USE DG
      USE HARM
      USE SIZES
      USE WIND
      USE OWIWIND,ONLY : NWS12INIT,NWS12GET
      USE NodalAttributes
#ifdef SWAN
Casey 101118: We need these values from other places.
      USE OWIWIND,     ONLY: WindMultiplier
      USE Couple2Swan, ONLY: COUPWIND,
     &                       SWAN_WX2,
     &                       SWAN_WY2
#endif
      
      IMPLICIT NONE
      
      REAL(SZ) WindDragLimit
      INTEGER II, IT

C.....Set the wind drag limit

      WindDragLimit = 0.002
      RampMete = rampdg

Casey 130710: Added this section.
      IF(WTIME1.LT.ITHS*DTDP)THEN
         WTIME1 = ITHS*DTDP
         WTIME2 = WTIME1 + WTIMINC
      ENDIF
C-----------------------------------------------------------------------
C
C     NWS = 1
C
C     Wind stress and atmospheric pressure are read in at all grid nodes
C     at every model time step from the fort.22 file.
C
C-----------------------------------------------------------------------

      IF (NWS.EQ.1) THEN
         DO II= 1,NP
         
C..........Read in the data
         
           READ(22,*) NHG, WSX2(II), WSY2(II), PR2(II)
           
C..........Apply the met ramp function

c           RampMete = RAMPDG
           
           WSX2(II)    = RampMete*WSX2(II)
           WSY2(II)    = RampMete*WSY2(II)
           PR2(II)     = RampMete*PR2(II)
           WVNXOUT(II) = WSX2(II)
           WVNYOUT(II) = WSY2(II)
         ENDDO
      ENDIF
      
C-----------------------------------------------------------------------
C
C     NWS = 2
C
C     Wind stress and atmospheric pressure are read in at all grid nodes
C     at a time interval that does not equal the model time step.  In-
C     terpolation in time is used to synchronize the wind and pressure
C     information with the model time step.
C
C-----------------------------------------------------------------------

      IF (ABS(NWS).EQ.2) THEN
      
C.......Determine if the met file time increment is exceeded
      
        IF (TIME_A.GT.WTIME2) THEN
          WTIME1 = WTIME2
          WTIME2 = WTIME2 + WTIMINC
          DO II= 1,NP
          
C...........Shift current data to old
          
            WVNX1(II) = WVNX2(II)
            WVNY1(II) = WVNY2(II)
            PRN1(II)  = PRN2(II)
            
C...........Read in data
            
            READ(22,*) NHG, WVNX2(II), WVNY2(II), PRN2(II)
          ENDDO
          PRINT*,'READING IN WIND DATA SET AT TIMESTEP',IT
        ENDIF
        
        WTRATIO = (TIME_A - WTIME1)/WTIMINC
        DO II= 1,NP
        
C.........Interpolate in time
        
          WINDX      = WVNX1(II) + WTRATIO*(WVNX2(II) - WVNX1(II))
          WINDY      = WVNY1(II) + WTRATIO*(WVNY2(II) - WVNY1(II))
          
C.........Apply mete ramp

c          RampMete = RAMPDG
          
          WSX2(II)    = RampMete*WINDX
          WSY2(II)    = RampMete*WINDY
          PR2(II)     = RampMete*(PRN1(II)+WTRATIO*(PRN2(II)-PRN1(II)))
          WVNXOUT(II) = WSX2(II)
          WVNYOUT(II) = WSY2(II)
        ENDDO
      ENDIF
      
C-----------------------------------------------------------------------
C
C     NWS = 3
C
C     Wind velocity in US Navy Fleet Numeric format interpolated in
C     space onto the ADCIRC grid and in time to synchronize the wind and
C     pressure information with the model time step. Garratt's formula
C     is used to compute wind stress from the wind velocity.
C
C-----------------------------------------------------------------------

      IF (NWS.EQ.3) THEN
      
C.......Determine if the met file time increment is exceeded

        IF (TIME_A.GT.WTIME2) THEN
          WTIME1 = WTIME2
          WTIME2 = WTIME2 + WTIMINC
          
C.........Shift current data to old
          
          DOII= 1,NP
            WVNX1(II) = WVNX2(II)
            WVNY1(II) = WVNY2(II)
          ENDDO
          
C.........Obtain the meteorological forcing data
          
          CALL NWS3GET( X, Y, SLAM, SFEA, WVNX2, WVNY2, IWTIME, IWYR,
     &                  WTIMED, NP, NWLON, NWLAT, WLATMAX, WLONMIN,
     &                  WLATINC, WLONINC, ICS, NSCREEN, ScreenUnit )
         ENDIF

         WTRATIO = (TIME_A - WTIME1)/WTIMINC
         DO II= 1,NP
         
C..........Interpolate in time
         
           WINDX   = WVNX1(II) + WTRATIO*(WVNX2(II) - WVNX1(II))
           WINDY   = WVNY1(II) + WTRATIO*(WVNY2(II) - WVNY1(II))
           
C..........Compute wind drag
           
           WINDMAG = SQRT( WINDX*WINDX + WINDY*WINDY )
           WDRAGCO = WindDrag( WINDMAG, WindDragLimit, "Garratt   ")
           
C..........Apply directional wind reductions
           
           IF (LoadDirEffRLen) THEN
             CALL ApplyDirectionalWindReduction( II, WDRAGCO, WINDMAG,
     &                                          DP(II), ETA2(II), H0, G,
     &                                          WINDX, WINDY )
             WINDMAG = SQRT( WINDX*WINDX + WINDY*WINDY )
             WDRAGCO = WindDrag(WINDMAG, WindDragLimit, "Garratt   " )
           ENDIF
           
C..........Apply met ramp
c           RampMete = RAMPDG
           WSX2(II)    = RampMete*0.001293D0*WDRAGCO*WINDX*WINDMAG
           WSY2(II)    = RampMete*0.001293D0*WDRAGCO*WINDY*WINDMAG
           WVNXOUT(II) = RampMete*WINDX
           WVNYOUT(II) = RampMete*WINDY
#ifdef SWAN
Casey 101118: Added these lines for coupling winds to SWAN.
           IF(COUPWIND)THEN
             SWAN_WX2(II,2) = WINDX
             SWAN_WY2(II,2) = WINDY
           ENDIF
#endif
        ENDDO
      ENDIF
      
C-----------------------------------------------------------------------
C
C     NWS = 4
C
C     Wind velocity and atmospheric pressure are read in (PBL/JAG
C     format) at selected ADCIRC grid nodes. Interpolation in time is
C     used to synchronize the wind and pressure information with the
C     model time step. Garratt's formula is used to compute wind stress
C     from wind velocity.
C
C-----------------------------------------------------------------------

      IF (ABS(NWS).EQ.4) THEN
      
C.......Determine if the met file time increment is exceeded
      
        IF (TIME_A.GT.WTIME2) THEN
          WTIME1 = WTIME2
          WTIME2 = WTIME2 + WTIMINC

C.........Shift current data to old
          
          DO II = 1,NP
            WVNX1(II) = WVNX2(II)
            WVNY1(II) = WVNY2(II)
            PRN1(II)  = PRN2(II)
          ENDDO
          
C.........Obtain the meteorological forcing data
          
          CALL NWS4GET( WVNX2, WVNY2, PRN2, NP, RHOWAT0, G )
        ENDIF

        WTRATIO = (TIME_A-WTIME1)/WTIMINC
        DO II = 1,NP
         
C.........Interpolate in time
         
          WINDX = WVNX1(II) + WTRATIO*(WVNX2(II)-WVNX1(II))
          WINDY = WVNY1(II) + WTRATIO*(WVNY2(II)-WVNY1(II))
            
C.........Compute wind drag
            
          WINDMAG = SQRT( WINDX*WINDX + WINDY*WINDY )
          WDRAGCO = WindDrag( WINDMAG, WindDragLimit, "Garratt   " )
           
C.........Apply directional wind reductions

          IF (LoadDirEffRLen) THEN
            CALL ApplyDirectionalWindReduction( II, WDRAGCO, WINDMAG,
     &                                          DP(II), ETA2(II), H0, G,
     &                                          WINDX, WINDY )
            WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
            WDRAGCO = WindDrag(WINDMAG, WindDragLimit, "Garratt   ")
          ENDIF
            
C.........Apply met ramp
c           RampMete = RAMPDG
          WSX2(II)    = RampMete*0.001293d0*WDRAGCO*WINDX*WINDMAG
          WSY2(II)    = RampMete*0.001293d0*WDRAGCO*WINDY*WINDMAG
          PR2(II)     = RampMete*(PRN1(II)+WTRATIO*(PRN2(II)- PRN1(II)))
          WVNXOUT(II) = RampMete*WINDX
          WVNYOUT(II) = RampMete*WINDY
#ifdef SWAN
Casey 101118: Added these lines for coupling winds to SWAN.
          IF(COUPWIND)THEN
            SWAN_WX2(II,2) = WINDX
            SWAN_WY2(II,2) = WINDY
          ENDIF
#endif
        ENDDO
      ENDIF
      
C-----------------------------------------------------------------------
C
C     NWS = 5
C
C     Wind velocity and atmospheric pressure are read in at all grid
C     nodes. Interpolation in time is used to synchronize the wind and
C     pressure information with the model time step. Garratt's formula
C     is used to compute wind stress from wind velocity.
C
C-----------------------------------------------------------------------

      IF(ABS(NWS).EQ.5) THEN
      
C.......Determine if the met file time increment is exceeded
      
        IF(TIME_A.GT.WTIME2) THEN
          WTIME1 = WTIME2
          WTIME2 = WTIME2 + WTIMINC
          
C.........Shift current data to old
          
          DO II = 1,NP
            WVNX1(II) = WVNX2(II)
            WVNY1(II) = WVNY2(II)
            PRN1(II)  = PRN2(II)
            
C...........Read in the meteorological forcing data
            
            READ(22,*) NHG, WVNX2(II), WVNY2(II), PRN2(II)
          ENDDO
        ENDIF
        
        WTRATIO = (TIME_A - WTIME1)/WTIMINC
        DO II = 1,NP
        
C.........Interpolate in time
        
          WINDX   = WVNX1(II) + WTRATIO*(WVNX2(II) - WVNX1(II))
          WINDY   = WVNY1(II) + WTRATIO*(WVNY2(II) - WVNY1(II))
          
C.........Compute wind drag
          
          WINDMAG = SQRT( WINDX*WINDX + WINDY*WINDY )
          WDRAGCO = WindDrag( WINDMAG, WindDragLimit, "Garratt   " )
          
C.........Apply directional wind reductions
          
          IF (LoadDirEffRLen) THEN
            CALL ApplyDirectionalWindReduction( II, WDRAGCO, WINDMAG,
     &                                          DP(II), ETA2(II), H0, G,
     &                                          WINDX, WINDY )
            WINDMAG = SQRT( WINDX*WINDX + WINDY*WINDY )
            WDRAGCO = WindDrag(WINDMAG, WindDragLimit, "Garratt   ")
          ENDIF
          
C.........Apply met ramp
          
          WSX2(II)    = RampMete*0.001293d0*WDRAGCO*WINDX*WINDMAG
          WSY2(II)    = RampMete*0.001293d0*WDRAGCO*WINDY*WINDMAG
          PR2(II)     = RampMete*(PRN1(II)+WTRATIO*(PRN2(II)-PRN1(II)))
          WVNXOUT(II) = RampMete*WINDX
          WVNYOUT(II) = RampMete*WINDY
#ifdef SWAN
Casey 101118: Added these lines for coupling winds to SWAN.
          IF(COUPWIND)THEN
            SWAN_WX2(II,2) = WINDX
            SWAN_WY2(II,2) = WINDY
          ENDIF
#endif
        ENDDO
      ENDIF
      
C-----------------------------------------------------------------------
C
C     NWS = 6
C
C     Wind velocity and atmospheric pressure are read in for a
C     rectangular grid (either in Longitude, Latitude or Cartesian
C     coordinates, consistent with the grid coordinates) and
C     interpolated in space onto the ADCIRC grid and in time to
C     synchronize the wind and pressure information with the model time
C     step. Garratt's formula is used to compute wind stress from the
C     wind velocity.
C
C-----------------------------------------------------------------------

      IF (NWS.EQ.6) THEN
      
C.......Determine if the met file time increment is exceeded

        IF (TIME_A.GT.WTIME2) THEN
          WTIME1 = WTIME2
          WTIME2 = WTIME2 + WTIMINC
          
C.........Shift current data to old
          
          DO II= 1,NP
            WVNX1(II) = WVNX2(II)
            WVNY1(II) = WVNY2(II)
            PRN1(II)  = PRN2(II)
          ENDDO
          NWSGGWI = NWSGGWI + 1
          
C.........Obtain meteorological forcing data
          
          CALL NWS6GET( X, Y, SLAM, SFEA, WVNX2, WVNY2, PRN2, NP,
     &                  NWLON, NWLAT, WLATMAX, WLONMIN, WLATINC,
     &                  WLONINC, ICS, RHOWAT0, G )
        ENDIF
        
        WTRATIO=(TIME_A-WTIME1)/WTIMINC
        DO II= 1,NP
        
C.........Interpolate in time
        
          WINDX = WVNX1(II) + WTRATIO*(WVNX2(II)-WVNX1(II))
          WINDY = WVNY1(II) + WTRATIO*(WVNY2(II)-WVNY1(II))
          
C.........Compute wind drag
          
          WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
          WDRAGCO = WindDrag( WINDMAG, WindDragLimit, "Garratt   " )
          
C.........Apply directional wind reductions

          IF (LoadDirEffRLen) THEN
            CALL ApplyDirectionalWindReduction( II, WDRAGCO, WINDMAG,
     &                                          DP(II), ETA2(II), H0, G,
     &                                          WINDX, WINDY )
            WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
            WDRAGCO = WindDrag(WINDMAG, WindDragLimit, "Garratt   ")
          ENDIF
          
C.........Apply met ramp
          
          WSX2(II)    = RampMete*0.001293d0*WDRAGCO*WINDX*WINDMAG
          WSY2(II)    = RampMete*0.001293d0*WDRAGCO*WINDY*WINDMAG
          PR2(II)     = RampMete*(PRN1(II)+WTRATIO*(PRN2(II)-PRN1(II)))
          WVNXOUT(II) = RampMete*WINDX
          WVNYOUT(II) = RampMete*WINDY
#ifdef SWAN
Casey 101118: Added these lines for coupling winds to SWAN.
          IF(COUPWIND)THEN
            SWAN_WX2(II,2) = WINDX
            SWAN_WY2(II,2) = WINDY
          ENDIF
#endif
        ENDDO
      ENDIF
      
C-----------------------------------------------------------------------
C
C     NWS = 7
C
C     jgf46.01 New option to read in surface wind stress and atmospheric
C     pressure for a rectangular grid (either in Longitude, Latitude or
C     Cartesian coordinates, consistent with the grid coordinates) and
C     interpolate in space onto the ADCIRC grid. Interpolation in time
C     is used to synchronize the wind and pressure information with the
C     model time step.
C
C-----------------------------------------------------------------------

      IF(ABS(NWS).EQ.7) THEN
      
C.......Determine if the met file time increment is exceeded
      
        IF (TIME_A.GT.WTIME2) THEN
          WTIME1 = WTIME2
          WTIME2 = WTIME2 + WTIMINC
          
C.........Shift current data to old
          
          DO II= 1,NP
            WVNX1(II) = WVNX2(II)
            WVNY1(II) = WVNY2(II)
            PRN1(II)  = PRN2(II)
          ENDDO
          
C.........Obtain the meteorological forcing data
          
          CALL NWS7GET( X, Y, SLAM, SFEA, WVNX2, WVNY2, PRN2, NP, NWLON,
     &                  NWLAT, WLATMAX, WLONMIN, WLATINC, WLONINC, ICS,
     &                  RHOWAT0,G )
        ENDIF

        WTRATIO=(TIME_A-WTIME1)/WTIMINC
        DO II= 1,NP
        
C.........Interpolate in time
        
          WINDX = WVNX1(II) + WTRATIO*(WVNX2(II) - WVNX1(II))
          WINDY = WVNY1(II) + WTRATIO*(WVNY2(II) - WVNY1(II))
          
C.........Apply met ramp
          
          WSX2(II)    = RampMete*WINDX
          WSY2(II)    = RampMete*WINDY
          PR2(II)     = RampMete*(PRN1(II)+WTRATIO*(PRN2(II)-PRN1(II)))
          WVNXOUT(II) = WSX2(II)
          WVNYOUT(II) = WSY2(II)
        ENDDO
      ENDIF
      
C-----------------------------------------------------------------------
C
C     NWS = 8
C
C     jgf46.02 New option to read in hurricane locations and generate
C     hurricane winds from the Holland Wind Model.
C
C-----------------------------------------------------------------------

      IF (ABS(NWS).EQ.8) THEN
      
C.......Obtain the meteorological forcing data
c         write(*,*) 'calling HollandGet ',time_a

        CALL HollandGet( X, Y, SLAM, SFEA, WVNX2, WVNY2, PRN2, NP, ICS,
     &                   RHOWAT0, G, TIME_A, NSCREEN, ScreenUnit )
        DO II= 1,NP
          WINDX = WVNX2(II)
          WINDY = WVNY2(II)
          
C.........Compute wind drag
          
          WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
          WDRAGCO = WindDrag( WINDMAG, WindDragLimit, "Garratt   " )
          
C.........Apply directional wind reductions
          
          IF (LoadDirEffRLen) THEN
            CALL ApplyDirectionalWindReduction( II, WDRAGCO, WINDMAG,
     &                                          DP(II), ETA2(II), H0, G,
     &                                          WINDX, WINDY )
            WINDMAG = SQRT(WINDX*WINDX + WINDY*WINDY)
            WDRAGCO = WindDrag( WINDMAG, WindDragLimit, "Garratt   " )
          ENDIF
            
C.........Apply met ramp
            
          WSX2(II)    = RampMete*0.001293d0*WDRAGCO*WINDX*WINDMAG
          WSY2(II)    = RampMete*0.001293d0*WDRAGCO*WINDY*WINDMAG
          PR2(II)     = RampMete*PRN2(II)
          WVNXOUT(II) = RampMete*WINDX
          WVNYOUT(II) = RampMete*WINDY
#ifdef SWAN
Casey 101118: Added these lines for coupling winds to SWAN.
          IF(COUPWIND)THEN
            SWAN_WX2(II,2) = WINDX
            SWAN_WY2(II,2) = WINDY
          ENDIF
#endif
        ENDDO
      ENDIF
      
C-----------------------------------------------------------------------
C
C     NWS = 9
C
C     jgf46.16 Merged:
C     cf & cm added nws = 9: asymmetric hurricane winds
C
C-----------------------------------------------------------------------

C      IF (NWS.EQ.9) THEN
      
C.......Obtain meteorological forcing data
      
C        CALL NWS9GET(SLAM,SFEA,WVNX2,WVNY2,PRN2,NP,TIME_A, ICS)
C        DO II= 1,NP
C          WINDX = WVNX2(II)
C          WINDY = WVNY2(II)
          
C.........Compute wind drag
          
C          WINDMAG = SQRT(WINDX*WINDX + WINDY*WINDY)
C          WDRAGCO = WindDrag( WINDMAG, WindDragLimit, "Garratt   " )
          
C.........Apply directional wind reductions
          
C          IF (LoadDirEffRLen) THEN
C            CALL ApplyDirectionalWindReduction( II, WDRAGCO, WINDMAG,
C     &                                          DP(II), ETA2(II), H0, G,
C     &                                          WINDX, WINDY )
C            WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
C            WDRAGCO = WindDrag(WINDMAG, WindDragLimit, "Garratt   ")
C          ENDIF
          
C.........Apply met ramp

C          WSX2(II)    = RampMete*0.001293d0*WDRAGCO*WINDX*WINDMAG
C          WSY2(II)    = RampMete*0.001293d0*WDRAGCO*WINDY*WINDMAG
C          PR2(II)     = RampMete*PRN2(II)
C          WVNXOUT(II) = RampMete*WINDX
C          WVNYOUT(II) = RampMete*WINDY
#ifdef SWAN
Casey 101118: Added these lines for coupling winds to SWAN.
C          IF(COUPWIND)THEN
C            SWAN_WX2(II,2) = WINDX
C            SWAN_WY2(II,2) = WINDY
C          ENDIF
#endif
C         ENDDO
C      ENDIF
      
C-----------------------------------------------------------------------
C
C     NWS = 10
C
C     Wind velocity (10 m) and atmospheric pressure are read in from a
C     sequence of National Weather Service (NWS) Aviation (AVN) model
C     output files. Each AVN file is assumed to contain data on a
C     Gaussian longitude, latitude grid at a single time. Consecutive
C     files in the sequence are separated by N hours in time. Garratt's
C     formula is used to compute wind stress from the wind velocity.
C
C-----------------------------------------------------------------------

      IF (NWS.EQ.10) THEN
      
C.......Determine if the met file time increment is exceeded
      
        IF (TIME_A.GT.WTIME2) THEN
          WTIME1 = WTIME2
          WTIME2 = WTIME2 + WTIMINC
          
C.........Shift current data to old
          
          DO II= 1,NP
            WVNX1(II) = WVNX2(II)
            WVNY1(II) = WVNY2(II)
            PRN1(II)  = PRN2(II)
          ENDDO
          NWSGGWI = NWSGGWI + 1
          
C.........Obtain meteorological forcing data
          
          CALL NWS10GET( NWSGGWI, SLAM, SFEA, WVNX2, WVNY2, PRN2, NP,
     &                   RHOWAT0, G, NWLON, NWLAT, WTIMINC )
        ENDIF
        
        WTRATIO = (TIME_A - WTIME1)/WTIMINC
        DO II = 1,NP
        
C.........Interpolate in time
        
          WINDX = WVNX1(II) + WTRATIO*(WVNX2(II)-WVNX1(II))
          WINDY = WVNY1(II) + WTRATIO*(WVNY2(II)-WVNY1(II))
           
C.........Compute wind drag
           
          WINDMAG = SQRT(WINDX*WINDX + WINDY*WINDY)
          WDRAGCO = WindDrag( WINDMAG, WindDragLimit, "Garratt   " )

C.........Apply directional wind reductions

          IF (LoadDirEffRLen) THEN
            CALL ApplyDirectionalWindReduction( II, WDRAGCO, WINDMAG,
     &                                          DP(II), ETA2(II), H0, G,
     &                                          WINDX, WINDY )
            WINDMAG = SQRT(WINDX*WINDX + WINDY*WINDY)
            WDRAGCO = WindDrag( WINDMAG, WindDragLimit, "Garratt   " )
          ENDIF
          
C.........Apply met ramp
          
          WSX2(II)    = RampMete*0.001293d0*WDRAGCO*WINDX*WINDMAG
          WSY2(II)    = RampMete*0.001293d0*WDRAGCO*WINDY*WINDMAG
          PR2(II)     = RampMete*(PRN1(II)+WTRATIO*(PRN2(II)-PRN1(II)))
          WVNXOUT(II) = RampMete*WINDX
          WVNYOUT(II) = RampMete*WINDY
#ifdef SWAN
Casey 101118: Added these lines for coupling winds to SWAN.
          IF(COUPWIND)THEN
            SWAN_WX2(II,2) = WINDX
            SWAN_WY2(II,2) = WINDY
          ENDIF
#endif
        ENDDO
      ENDIF
      
C-----------------------------------------------------------------------
C
C     NWS = 11
C
C     Wind velocity (10 m) and atmospheric pressure are read in from a
C     sequence of stripped down (?) National Weather Service (NWS) ETA
C     29km model output files. Each ETA file is assumed to contain data
C     on an E grid for a single day (8 data sets, one every 3 hours, be-
C     ginning @ 03:00 and continuing through 24:00 of the given day).
C     The wind data is converted to an east-west, north-south coordinate
C     system inside ADCIRC. Garratt's formula is used to compute wind
C     stress from the wind velocity.
C
C-----------------------------------------------------------------------

      IF(NWS.EQ.11) THEN

C.......Determine if the met file time increment is exceeded

        IF (TIME_A.GT.WTIME2) THEN
          WTIME1=WTIME2
          WTIME2=WTIME2+WTIMINC
           
C........Shift current data to old
           
          DO II = 1,NP
            WVNX1(II) = WVNX2(II)
            WVNY1(II) = WVNY2(II)
            PRN1(II)  = PRN2(II)
          ENDDO
          IDSETFLG = IDSETFLG + 1
          IF (IDSETFLG.GT.8) THEN
            NWSEGWI = NWSEGWI + 1
            IDSETFLG = 1
          ENDIF
            
C.........Obtain meteorological forcing data
            
          CALL NWS11GET( NWSEGWI, IDSETFLG, SLAM, SFEA, WVNX2, WVNY2,
     &                   PRN2, NP, RHOWAT0, G )
        ENDIF

        WTRATIO=(TIME_A-WTIME1)/WTIMINC
        DO II = 1,NP
         
C.........Interpolate in time
         
          WINDX = WVNX1(II) + WTRATIO*(WVNX2(II)-WVNX1(II))
          WINDY = WVNY1(II) + WTRATIO*(WVNY2(II)-WVNY1(II))
            
C.........Compute wind drag
            
          WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
          WDRAGCO = WindDrag(WINDMAG, WindDragLimit, "Garratt   ")
            
C.........Apply directional wind reductions

          IF (LoadDirEffRLen) THEN
            CALL ApplyDirectionalWindReduction( II, WDRAGCO, WINDMAG,
     &                                          DP(II), ETA2(II), H0, G,
     &                                          WINDX, WINDY )
            WINDMAG = SQRT(WINDX*WINDX + WINDY*WINDY)
            WDRAGCO = WindDrag( WINDMAG, WindDragLimit, "Garratt   " )
          ENDIF
            
C.........Apply met ramp
            
          WSX2(II)    = RampMete*0.001293d0*WDRAGCO*WINDX*WINDMAG
          WSY2(II)    = RampMete*0.001293d0*WDRAGCO*WINDY*WINDMAG
          PR2(II)     = RampMete*(PRN1(II)+WTRATIO*(PRN2(II)-PRN1(II)))
          WVNXOUT(II) = RampMete*WINDX
          WVNYOUT(II) = RampMete*WINDY
#ifdef SWAN
Casey 101118: Added these lines for coupling winds to SWAN.
          IF(COUPWIND)THEN
            SWAN_WX2(II,2) = WINDX
            SWAN_WY2(II,2) = WINDY
          ENDIF
#endif
        ENDDO
      ENDIF
      
C-----------------------------------------------------------------------
C
C     NWS = 12
C
C     sb46.28sb01 NWS=12 reads in raw OWI files 09/xx/2006
C
C-----------------------------------------------------------------------

      IF(ABS(NWS).EQ.12) THEN
      
C.......Determine if the met file time increment is exceeded
      
        IF(TIME_A.GT.WTIME2) THEN
          WTIME1=WTIME2
          WTIME2=WTIME2+WTIMINC
          
C........Shift current data to old

          DO II =1,NP
            WVNX1(II) = WVNX2(II)
            WVNY1(II) = WVNY2(II)
            PRN1(II)  = PRN2(II)
          ENDDO
          
C.........Obtain meteorological forcing data
          
          CALL NWS12GET( WVNX2, WVNY2, PRN2, NP, RHOWAT0, G )
        ENDIF

        WTRATIO=(TIME_A - WTIME1)/WTIMINC
        DO II = 1,NP
        
C.........Interpolate in time
        
          WINDX = WVNX1(II) + WTRATIO*(WVNX2(II)-WVNX1(II))
          WINDY = WVNY1(II) + WTRATIO*(WVNY2(II)-WVNY1(II))
          
C.........Compute wind drag
          
          WINDMAG = SQRT(WINDX*WINDX + WINDY*WINDY)
          WDRAGCO = WindDrag( WINDMAG, WindDragLimit, "Garratt   " )
          
C.........Apply directional wind reductions
          
          IF (LoadDirEffRLen) THEN
            CALL ApplyDirectionalWindReduction( II, WDRAGCO, WINDMAG,
     &                                          DP(II), ETA2(II), H0, G,
     &                                          WINDX, WINDY )
            WINDMAG = SQRT(WINDX*WINDX + WINDY*WINDY)
            WDRAGCO = WindDrag(WINDMAG, WindDragLimit, "Garratt   ")
          ENDIF
          
C.........Apply met ramp
          
          WSX2(II)    = RampMete*0.001293d0*WDRAGCO*WINDX*WINDMAG
          WSY2(II)    = RampMete*0.001293d0*WDRAGCO*WINDY*WINDMAG
          PR2(II)     = RampMete*(PRN1(II)+WTRATIO*(PRN2(II)-PRN1(II)))
          WVNXOUT(II) = RampMete*WINDX
          WVNYOUT(II) = RampMete*WINDY
#ifdef SWAN
Casey 101118: Added these lines for coupling winds to SWAN.
          IF(COUPWIND)THEN
            SWAN_WX2(II,2) = WINDX/WindMultiplier
            SWAN_WY2(II,2) = WINDY/WindMultiplier
          ENDIF
#endif
        ENDDO
      ENDIF
      
      RETURN
      END SUBROUTINE MET_FORCING
