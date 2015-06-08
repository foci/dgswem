!***********************************************************************
!
!     SUBROUTINE MET_FORCING()
!
!     This subroutine handles the meteorological forcing for the DG
!     code
!
!***********************************************************************

      SUBROUTINE MET_FORCING(IT)
      
      USE GLOBAL
      USE DG
      USE HARM
      USE SIZES
      USE WIND
      USE OWIWIND,ONLY : NWS12INIT,NWS12GET
      USE NodalAttributes
#ifdef SWAN
!asey 101118: We need these values from other places.
      USE OWIWIND,     ONLY: WindMultiplier
      USE Couple2Swan, ONLY: COUPWIND,&
                       SWAN_WX2,&
                       SWAN_WY2
#endif
      
      IMPLICIT NONE
      
      REAL(SZ) WindDragLimit
      INTEGER II, IT

!.....Set the wind drag limit

      WindDragLimit = 0.002
      RampMete = rampdg

!asey 130710: Added this section.
      IF(WTIME1.LT.ITHS*DTDP)THEN
         WTIME1 = ITHS*DTDP
         WTIME2 = WTIME1 + WTIMINC
      ENDIF
!-----------------------------------------------------------------------
!
!     NWS = 1
!
!     Wind stress and atmospheric pressure are read in at all grid nodes
!     at every model time step from the fort.22 file.
!
!-----------------------------------------------------------------------

      IF (NWS.EQ.1) THEN
         DO II= 1,NP
         
!..........Read in the data
         
           READ(22,*) NHG, WSX2(II), WSY2(II), PR2(II)
           
!..........Apply the met ramp function

!           RampMete = RAMPDG
           
           WSX2(II)    = RampMete*WSX2(II)
           WSY2(II)    = RampMete*WSY2(II)
           PR2(II)     = RampMete*PR2(II)
           WVNXOUT(II) = WSX2(II)
           WVNYOUT(II) = WSY2(II)
         ENDDO
      ENDIF
      
!-----------------------------------------------------------------------
!
!     NWS = 2
!
!     Wind stress and atmospheric pressure are read in at all grid nodes
!     at a time interval that does not equal the model time step.  In-
!     terpolation in time is used to synchronize the wind and pressure
!     information with the model time step.
!
!-----------------------------------------------------------------------

      IF (ABS(NWS).EQ.2) THEN
      
!.......Determine if the met file time increment is exceeded
      
        IF (TIME_A.GT.WTIME2) THEN
          WTIME1 = WTIME2
          WTIME2 = WTIME2 + WTIMINC
          DO II= 1,NP
          
!...........Shift current data to old
          
            WVNX1(II) = WVNX2(II)
            WVNY1(II) = WVNY2(II)
            PRN1(II)  = PRN2(II)
            
!...........Read in data
            
            READ(22,*) NHG, WVNX2(II), WVNY2(II), PRN2(II)
          ENDDO
          PRINT*,'READING IN WIND DATA SET AT TIMESTEP',IT
        ENDIF
        
        WTRATIO = (TIME_A - WTIME1)/WTIMINC
        DO II= 1,NP
        
!.........Interpolate in time
        
          WINDX      = WVNX1(II) + WTRATIO*(WVNX2(II) - WVNX1(II))
          WINDY      = WVNY1(II) + WTRATIO*(WVNY2(II) - WVNY1(II))
          
!.........Apply mete ramp

!          RampMete = RAMPDG
          
          WSX2(II)    = RampMete*WINDX
          WSY2(II)    = RampMete*WINDY
          PR2(II)     = RampMete*(PRN1(II)+WTRATIO*(PRN2(II)-PRN1(II)))
          WVNXOUT(II) = WSX2(II)
          WVNYOUT(II) = WSY2(II)
        ENDDO
      ENDIF
      
!-----------------------------------------------------------------------
!
!     NWS = 3
!
!     Wind velocity in US Navy Fleet Numeric format interpolated in
!     space onto the ADCIRC grid and in time to synchronize the wind and
!     pressure information with the model time step. Garratt's formula
!     is used to compute wind stress from the wind velocity.
!
!-----------------------------------------------------------------------

      IF (NWS.EQ.3) THEN
      
!.......Determine if the met file time increment is exceeded

        IF (TIME_A.GT.WTIME2) THEN
          WTIME1 = WTIME2
          WTIME2 = WTIME2 + WTIMINC
          
!.........Shift current data to old
          
          DO II=1,NP
            WVNX1(II) = WVNX2(II)
            WVNY1(II) = WVNY2(II)
          ENDDO
          
!.........Obtain the meteorological forcing data
          
          CALL NWS3GET( X, Y, SLAM, SFEA, WVNX2, WVNY2, IWTIME, IWYR,&
                  WTIMED, NP, NWLON, NWLAT, WLATMAX, WLONMIN,&
                  WLATINC, WLONINC, ICS, NSCREEN, ScreenUnit )
         ENDIF

         WTRATIO = (TIME_A - WTIME1)/WTIMINC
         DO II= 1,NP
         
!..........Interpolate in time
         
           WINDX   = WVNX1(II) + WTRATIO*(WVNX2(II) - WVNX1(II))
           WINDY   = WVNY1(II) + WTRATIO*(WVNY2(II) - WVNY1(II))
           
!..........Compute wind drag
           
           WINDMAG = SQRT( WINDX*WINDX + WINDY*WINDY )
           WDRAGCO = WindDrag( WINDMAG, WindDragLimit, "Garratt   ")
           
!..........Apply directional wind reductions
           
           IF (LoadDirEffRLen) THEN
             CALL ApplyDirectionalWindReduction( II, WDRAGCO, WINDMAG,&
                                          DP(II), ETA2(II), H0, G,&
                                          WINDX, WINDY )
             WINDMAG = SQRT( WINDX*WINDX + WINDY*WINDY )
             WDRAGCO = WindDrag(WINDMAG, WindDragLimit, "Garratt   " )
           ENDIF
           
!..........Apply met ramp
!           RampMete = RAMPDG
           WSX2(II)    = RampMete*0.001293D0*WDRAGCO*WINDX*WINDMAG
           WSY2(II)    = RampMete*0.001293D0*WDRAGCO*WINDY*WINDMAG
           WVNXOUT(II) = RampMete*WINDX
           WVNYOUT(II) = RampMete*WINDY
#ifdef SWAN
!asey 101118: Added these lines for coupling winds to SWAN.
           IF(COUPWIND)THEN
             SWAN_WX2(II,2) = WINDX
             SWAN_WY2(II,2) = WINDY
           ENDIF
#endif
        ENDDO
      ENDIF
      
!-----------------------------------------------------------------------
!
!     NWS = 4
!
!     Wind velocity and atmospheric pressure are read in (PBL/JAG
!     format) at selected ADCIRC grid nodes. Interpolation in time is
!     used to synchronize the wind and pressure information with the
!     model time step. Garratt's formula is used to compute wind stress
!     from wind velocity.
!
!-----------------------------------------------------------------------

      IF (ABS(NWS).EQ.4) THEN
      
!.......Determine if the met file time increment is exceeded
      
        IF (TIME_A.GT.WTIME2) THEN
          WTIME1 = WTIME2
          WTIME2 = WTIME2 + WTIMINC

!.........Shift current data to old
          
          DO II = 1,NP
            WVNX1(II) = WVNX2(II)
            WVNY1(II) = WVNY2(II)
            PRN1(II)  = PRN2(II)
          ENDDO
          
!.........Obtain the meteorological forcing data
          
          CALL NWS4GET( WVNX2, WVNY2, PRN2, NP, RHOWAT0, G )
        ENDIF

        WTRATIO = (TIME_A-WTIME1)/WTIMINC
        DO II = 1,NP
         
!.........Interpolate in time
         
          WINDX = WVNX1(II) + WTRATIO*(WVNX2(II)-WVNX1(II))
          WINDY = WVNY1(II) + WTRATIO*(WVNY2(II)-WVNY1(II))
            
!.........Compute wind drag
            
          WINDMAG = SQRT( WINDX*WINDX + WINDY*WINDY )
          WDRAGCO = WindDrag( WINDMAG, WindDragLimit, "Garratt   " )
           
!.........Apply directional wind reductions

          IF (LoadDirEffRLen) THEN
            CALL ApplyDirectionalWindReduction( II, WDRAGCO, WINDMAG,&
                                          DP(II), ETA2(II), H0, G,&
                                          WINDX, WINDY )
            WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
            WDRAGCO = WindDrag(WINDMAG, WindDragLimit, "Garratt   ")
          ENDIF
            
!.........Apply met ramp
!           RampMete = RAMPDG
          WSX2(II)    = RampMete*0.001293d0*WDRAGCO*WINDX*WINDMAG
          WSY2(II)    = RampMete*0.001293d0*WDRAGCO*WINDY*WINDMAG
          PR2(II)     = RampMete*(PRN1(II)+WTRATIO*(PRN2(II)- PRN1(II)))
          WVNXOUT(II) = RampMete*WINDX
          WVNYOUT(II) = RampMete*WINDY
#ifdef SWAN
!asey 101118: Added these lines for coupling winds to SWAN.
          IF(COUPWIND)THEN
            SWAN_WX2(II,2) = WINDX
            SWAN_WY2(II,2) = WINDY
          ENDIF
#endif
        ENDDO
      ENDIF
      
!-----------------------------------------------------------------------
!
!     NWS = 5
!
!     Wind velocity and atmospheric pressure are read in at all grid
!     nodes. Interpolation in time is used to synchronize the wind and
!     pressure information with the model time step. Garratt's formula
!     is used to compute wind stress from wind velocity.
!
!-----------------------------------------------------------------------

      IF(ABS(NWS).EQ.5) THEN
      
!.......Determine if the met file time increment is exceeded
      
        IF(TIME_A.GT.WTIME2) THEN
          WTIME1 = WTIME2
          WTIME2 = WTIME2 + WTIMINC
          
!.........Shift current data to old
          
          DO II = 1,NP
            WVNX1(II) = WVNX2(II)
            WVNY1(II) = WVNY2(II)
            PRN1(II)  = PRN2(II)
            
!...........Read in the meteorological forcing data
            
            READ(22,*) NHG, WVNX2(II), WVNY2(II), PRN2(II)
          ENDDO
        ENDIF
        
        WTRATIO = (TIME_A - WTIME1)/WTIMINC
        DO II = 1,NP
        
!.........Interpolate in time
        
          WINDX   = WVNX1(II) + WTRATIO*(WVNX2(II) - WVNX1(II))
          WINDY   = WVNY1(II) + WTRATIO*(WVNY2(II) - WVNY1(II))
          
!.........Compute wind drag
          
          WINDMAG = SQRT( WINDX*WINDX + WINDY*WINDY )
          WDRAGCO = WindDrag( WINDMAG, WindDragLimit, "Garratt   " )
          
!.........Apply directional wind reductions
          
          IF (LoadDirEffRLen) THEN
            CALL ApplyDirectionalWindReduction( II, WDRAGCO, WINDMAG,&
                                          DP(II), ETA2(II), H0, G,&
                                          WINDX, WINDY )
            WINDMAG = SQRT( WINDX*WINDX + WINDY*WINDY )
            WDRAGCO = WindDrag(WINDMAG, WindDragLimit, "Garratt   ")
          ENDIF
          
!.........Apply met ramp
          
          WSX2(II)    = RampMete*0.001293d0*WDRAGCO*WINDX*WINDMAG
          WSY2(II)    = RampMete*0.001293d0*WDRAGCO*WINDY*WINDMAG
          PR2(II)     = RampMete*(PRN1(II)+WTRATIO*(PRN2(II)-PRN1(II)))
          WVNXOUT(II) = RampMete*WINDX
          WVNYOUT(II) = RampMete*WINDY
#ifdef SWAN
!asey 101118: Added these lines for coupling winds to SWAN.
          IF(COUPWIND)THEN
            SWAN_WX2(II,2) = WINDX
            SWAN_WY2(II,2) = WINDY
          ENDIF
#endif
        ENDDO
      ENDIF
      
!-----------------------------------------------------------------------
!
!     NWS = 6
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

      IF (NWS.EQ.6) THEN
      
!.......Determine if the met file time increment is exceeded

        IF (TIME_A.GT.WTIME2) THEN
          WTIME1 = WTIME2
          WTIME2 = WTIME2 + WTIMINC
          
!.........Shift current data to old
          
          DO II= 1,NP
            WVNX1(II) = WVNX2(II)
            WVNY1(II) = WVNY2(II)
            PRN1(II)  = PRN2(II)
          ENDDO
          NWSGGWI = NWSGGWI + 1
          
!.........Obtain meteorological forcing data
          
          CALL NWS6GET( X, Y, SLAM, SFEA, WVNX2, WVNY2, PRN2, NP,&
                  NWLON, NWLAT, WLATMAX, WLONMIN, WLATINC,&
                  WLONINC, ICS, RHOWAT0, G )
        ENDIF
        
        WTRATIO=(TIME_A-WTIME1)/WTIMINC
        DO II= 1,NP
        
!.........Interpolate in time
        
          WINDX = WVNX1(II) + WTRATIO*(WVNX2(II)-WVNX1(II))
          WINDY = WVNY1(II) + WTRATIO*(WVNY2(II)-WVNY1(II))
          
!.........Compute wind drag
          
          WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
          WDRAGCO = WindDrag( WINDMAG, WindDragLimit, "Garratt   " )
          
!.........Apply directional wind reductions

          IF (LoadDirEffRLen) THEN
            CALL ApplyDirectionalWindReduction( II, WDRAGCO, WINDMAG,&
                                          DP(II), ETA2(II), H0, G,&
                                          WINDX, WINDY )
            WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
            WDRAGCO = WindDrag(WINDMAG, WindDragLimit, "Garratt   ")
          ENDIF
          
!.........Apply met ramp
          
          WSX2(II)    = RampMete*0.001293d0*WDRAGCO*WINDX*WINDMAG
          WSY2(II)    = RampMete*0.001293d0*WDRAGCO*WINDY*WINDMAG
          PR2(II)     = RampMete*(PRN1(II)+WTRATIO*(PRN2(II)-PRN1(II)))
          WVNXOUT(II) = RampMete*WINDX
          WVNYOUT(II) = RampMete*WINDY
#ifdef SWAN
!asey 101118: Added these lines for coupling winds to SWAN.
          IF(COUPWIND)THEN
            SWAN_WX2(II,2) = WINDX
            SWAN_WY2(II,2) = WINDY
          ENDIF
#endif
        ENDDO
      ENDIF
      
!-----------------------------------------------------------------------
!
!     NWS = 7
!
!     jgf46.01 New option to read in surface wind stress and atmospheric
!     pressure for a rectangular grid (either in Longitude, Latitude or
!     Cartesian coordinates, consistent with the grid coordinates) and
!     interpolate in space onto the ADCIRC grid. Interpolation in time
!     is used to synchronize the wind and pressure information with the
!     model time step.
!
!-----------------------------------------------------------------------

      IF(ABS(NWS).EQ.7) THEN
      
!.......Determine if the met file time increment is exceeded
      
        IF (TIME_A.GT.WTIME2) THEN
          WTIME1 = WTIME2
          WTIME2 = WTIME2 + WTIMINC
          
!.........Shift current data to old
          
          DO II= 1,NP
            WVNX1(II) = WVNX2(II)
            WVNY1(II) = WVNY2(II)
            PRN1(II)  = PRN2(II)
          ENDDO
          
!.........Obtain the meteorological forcing data
          
          CALL NWS7GET( X, Y, SLAM, SFEA, WVNX2, WVNY2, PRN2, NP, NWLON,&
                  NWLAT, WLATMAX, WLONMIN, WLATINC, WLONINC, ICS,&
                  RHOWAT0,G )
        ENDIF

        WTRATIO=(TIME_A-WTIME1)/WTIMINC
        DO II= 1,NP
        
!.........Interpolate in time
        
          WINDX = WVNX1(II) + WTRATIO*(WVNX2(II) - WVNX1(II))
          WINDY = WVNY1(II) + WTRATIO*(WVNY2(II) - WVNY1(II))
          
!.........Apply met ramp
          
          WSX2(II)    = RampMete*WINDX
          WSY2(II)    = RampMete*WINDY
          PR2(II)     = RampMete*(PRN1(II)+WTRATIO*(PRN2(II)-PRN1(II)))
          WVNXOUT(II) = WSX2(II)
          WVNYOUT(II) = WSY2(II)
        ENDDO
      ENDIF
      
!-----------------------------------------------------------------------
!
!     NWS = 8
!
!     jgf46.02 New option to read in hurricane locations and generate
!     hurricane winds from the Holland Wind Model.
!
!-----------------------------------------------------------------------

      IF (ABS(NWS).EQ.8) THEN
      
!.......Obtain the meteorological forcing data
!         write(*,*) 'calling HollandGet ',time_a

        CALL HollandGet( X, Y, SLAM, SFEA, WVNX2, WVNY2, PRN2, NP, ICS,&
                   RHOWAT0, G, TIME_A, NSCREEN, ScreenUnit )
        DO II= 1,NP
          WINDX = WVNX2(II)
          WINDY = WVNY2(II)
          
!.........Compute wind drag
          
          WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
          WDRAGCO = WindDrag( WINDMAG, WindDragLimit, "Garratt   " )
          
!.........Apply directional wind reductions
          
          IF (LoadDirEffRLen) THEN
            CALL ApplyDirectionalWindReduction( II, WDRAGCO, WINDMAG,&
                                          DP(II), ETA2(II), H0, G,&
                                          WINDX, WINDY )
            WINDMAG = SQRT(WINDX*WINDX + WINDY*WINDY)
            WDRAGCO = WindDrag( WINDMAG, WindDragLimit, "Garratt   " )
          ENDIF
            
!.........Apply met ramp
            
          WSX2(II)    = RampMete*0.001293d0*WDRAGCO*WINDX*WINDMAG
          WSY2(II)    = RampMete*0.001293d0*WDRAGCO*WINDY*WINDMAG
          PR2(II)     = RampMete*PRN2(II)
          WVNXOUT(II) = RampMete*WINDX
          WVNYOUT(II) = RampMete*WINDY
#ifdef SWAN
!asey 101118: Added these lines for coupling winds to SWAN.
          IF(COUPWIND)THEN
            SWAN_WX2(II,2) = WINDX
            SWAN_WY2(II,2) = WINDY
          ENDIF
#endif
        ENDDO
      ENDIF
      
!-----------------------------------------------------------------------
!
!     NWS = 9
!
!     jgf46.16 Merged:
!     cf & cm added nws = 9: asymmetric hurricane winds
!
!-----------------------------------------------------------------------

!      IF (NWS.EQ.9) THEN
      
!.......Obtain meteorological forcing data
      
!        CALL NWS9GET(SLAM,SFEA,WVNX2,WVNY2,PRN2,NP,TIME_A, ICS)
!        DO II= 1,NP
!          WINDX = WVNX2(II)
!          WINDY = WVNY2(II)
          
!.........Compute wind drag
          
!          WINDMAG = SQRT(WINDX*WINDX + WINDY*WINDY)
!          WDRAGCO = WindDrag( WINDMAG, WindDragLimit, "Garratt   " )
          
!.........Apply directional wind reductions
          
!          IF (LoadDirEffRLen) THEN
!            CALL ApplyDirectionalWindReduction( II, WDRAGCO, WINDMAG,
!     &                                          DP(II), ETA2(II), H0, G,
!     &                                          WINDX, WINDY )
!            WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
!            WDRAGCO = WindDrag(WINDMAG, WindDragLimit, "Garratt   ")
!          ENDIF
          
!.........Apply met ramp

!          WSX2(II)    = RampMete*0.001293d0*WDRAGCO*WINDX*WINDMAG
!          WSY2(II)    = RampMete*0.001293d0*WDRAGCO*WINDY*WINDMAG
!          PR2(II)     = RampMete*PRN2(II)
!          WVNXOUT(II) = RampMete*WINDX
!          WVNYOUT(II) = RampMete*WINDY
#ifdef SWAN
!asey 101118: Added these lines for coupling winds to SWAN.
!          IF(COUPWIND)THEN
!            SWAN_WX2(II,2) = WINDX
!            SWAN_WY2(II,2) = WINDY
!          ENDIF
#endif
!         ENDDO
!      ENDIF
      
!-----------------------------------------------------------------------
!
!     NWS = 10
!
!     Wind velocity (10 m) and atmospheric pressure are read in from a
!     sequence of National Weather Service (NWS) Aviation (AVN) model
!     output files. Each AVN file is assumed to contain data on a
!     Gaussian longitude, latitude grid at a single time. Consecutive
!     files in the sequence are separated by N hours in time. Garratt's
!     formula is used to compute wind stress from the wind velocity.
!
!-----------------------------------------------------------------------

      IF (NWS.EQ.10) THEN
      
!.......Determine if the met file time increment is exceeded
      
        IF (TIME_A.GT.WTIME2) THEN
          WTIME1 = WTIME2
          WTIME2 = WTIME2 + WTIMINC
          
!.........Shift current data to old
          
          DO II= 1,NP
            WVNX1(II) = WVNX2(II)
            WVNY1(II) = WVNY2(II)
            PRN1(II)  = PRN2(II)
          ENDDO
          NWSGGWI = NWSGGWI + 1
          
!.........Obtain meteorological forcing data
          
          CALL NWS10GET( NWSGGWI, SLAM, SFEA, WVNX2, WVNY2, PRN2, NP,&
                   RHOWAT0, G, NWLON, NWLAT, WTIMINC )
        ENDIF
        
        WTRATIO = (TIME_A - WTIME1)/WTIMINC
        DO II = 1,NP
        
!.........Interpolate in time
        
          WINDX = WVNX1(II) + WTRATIO*(WVNX2(II)-WVNX1(II))
          WINDY = WVNY1(II) + WTRATIO*(WVNY2(II)-WVNY1(II))
           
!.........Compute wind drag
           
          WINDMAG = SQRT(WINDX*WINDX + WINDY*WINDY)
          WDRAGCO = WindDrag( WINDMAG, WindDragLimit, "Garratt   " )

!.........Apply directional wind reductions

          IF (LoadDirEffRLen) THEN
            CALL ApplyDirectionalWindReduction( II, WDRAGCO, WINDMAG,&
                                          DP(II), ETA2(II), H0, G,&
                                          WINDX, WINDY )
            WINDMAG = SQRT(WINDX*WINDX + WINDY*WINDY)
            WDRAGCO = WindDrag( WINDMAG, WindDragLimit, "Garratt   " )
          ENDIF
          
!.........Apply met ramp
          
          WSX2(II)    = RampMete*0.001293d0*WDRAGCO*WINDX*WINDMAG
          WSY2(II)    = RampMete*0.001293d0*WDRAGCO*WINDY*WINDMAG
          PR2(II)     = RampMete*(PRN1(II)+WTRATIO*(PRN2(II)-PRN1(II)))
          WVNXOUT(II) = RampMete*WINDX
          WVNYOUT(II) = RampMete*WINDY
#ifdef SWAN
!asey 101118: Added these lines for coupling winds to SWAN.
          IF(COUPWIND)THEN
            SWAN_WX2(II,2) = WINDX
            SWAN_WY2(II,2) = WINDY
          ENDIF
#endif
        ENDDO
      ENDIF
      
!-----------------------------------------------------------------------
!
!     NWS = 11
!
!     Wind velocity (10 m) and atmospheric pressure are read in from a
!     sequence of stripped down (?) National Weather Service (NWS) ETA
!     29km model output files. Each ETA file is assumed to contain data
!     on an E grid for a single day (8 data sets, one every 3 hours, be-
!     ginning @ 03:00 and continuing through 24:00 of the given day).
!     The wind data is converted to an east-west, north-south coordinate
!     system inside ADCIRC. Garratt's formula is used to compute wind
!     stress from the wind velocity.
!
!-----------------------------------------------------------------------

      IF(NWS.EQ.11) THEN

!.......Determine if the met file time increment is exceeded

        IF (TIME_A.GT.WTIME2) THEN
          WTIME1=WTIME2
          WTIME2=WTIME2+WTIMINC
           
!........Shift current data to old
           
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
            
!.........Obtain meteorological forcing data
            
          CALL NWS11GET( NWSEGWI, IDSETFLG, SLAM, SFEA, WVNX2, WVNY2,&
                   PRN2, NP, RHOWAT0, G )
        ENDIF

        WTRATIO=(TIME_A-WTIME1)/WTIMINC
        DO II = 1,NP
         
!.........Interpolate in time
         
          WINDX = WVNX1(II) + WTRATIO*(WVNX2(II)-WVNX1(II))
          WINDY = WVNY1(II) + WTRATIO*(WVNY2(II)-WVNY1(II))
            
!.........Compute wind drag
            
          WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
          WDRAGCO = WindDrag(WINDMAG, WindDragLimit, "Garratt   ")
            
!.........Apply directional wind reductions

          IF (LoadDirEffRLen) THEN
            CALL ApplyDirectionalWindReduction( II, WDRAGCO, WINDMAG,&
                                          DP(II), ETA2(II), H0, G,&
                                          WINDX, WINDY )
            WINDMAG = SQRT(WINDX*WINDX + WINDY*WINDY)
            WDRAGCO = WindDrag( WINDMAG, WindDragLimit, "Garratt   " )
          ENDIF
            
!.........Apply met ramp
            
          WSX2(II)    = RampMete*0.001293d0*WDRAGCO*WINDX*WINDMAG
          WSY2(II)    = RampMete*0.001293d0*WDRAGCO*WINDY*WINDMAG
          PR2(II)     = RampMete*(PRN1(II)+WTRATIO*(PRN2(II)-PRN1(II)))
          WVNXOUT(II) = RampMete*WINDX
          WVNYOUT(II) = RampMete*WINDY
#ifdef SWAN
!asey 101118: Added these lines for coupling winds to SWAN.
          IF(COUPWIND)THEN
            SWAN_WX2(II,2) = WINDX
            SWAN_WY2(II,2) = WINDY
          ENDIF
#endif
        ENDDO
      ENDIF
      
!-----------------------------------------------------------------------
!
!     NWS = 12
!
!     sb46.28sb01 NWS=12 reads in raw OWI files 09/xx/2006
!
!-----------------------------------------------------------------------

      IF(ABS(NWS).EQ.12) THEN
      
!.......Determine if the met file time increment is exceeded
      
        IF(TIME_A.GT.WTIME2) THEN
          WTIME1=WTIME2
          WTIME2=WTIME2+WTIMINC
          
!........Shift current data to old

          DO II =1,NP
            WVNX1(II) = WVNX2(II)
            WVNY1(II) = WVNY2(II)
            PRN1(II)  = PRN2(II)
          ENDDO
          
!.........Obtain meteorological forcing data
          
          CALL NWS12GET( WVNX2, WVNY2, PRN2, NP, RHOWAT0, G )
        ENDIF

        WTRATIO=(TIME_A - WTIME1)/WTIMINC
        DO II = 1,NP
        
!.........Interpolate in time
        
          WINDX = WVNX1(II) + WTRATIO*(WVNX2(II)-WVNX1(II))
          WINDY = WVNY1(II) + WTRATIO*(WVNY2(II)-WVNY1(II))
          
!.........Compute wind drag
          
          WINDMAG = SQRT(WINDX*WINDX + WINDY*WINDY)
          WDRAGCO = WindDrag( WINDMAG, WindDragLimit, "Garratt   " )
          
!.........Apply directional wind reductions
          
          IF (LoadDirEffRLen) THEN
            CALL ApplyDirectionalWindReduction( II, WDRAGCO, WINDMAG,&
                                          DP(II), ETA2(II), H0, G,&
                                          WINDX, WINDY )
            WINDMAG = SQRT(WINDX*WINDX + WINDY*WINDY)
            WDRAGCO = WindDrag(WINDMAG, WindDragLimit, "Garratt   ")
          ENDIF
          
!.........Apply met ramp
          
          WSX2(II)    = RampMete*0.001293d0*WDRAGCO*WINDX*WINDMAG
          WSY2(II)    = RampMete*0.001293d0*WDRAGCO*WINDY*WINDMAG
          PR2(II)     = RampMete*(PRN1(II)+WTRATIO*(PRN2(II)-PRN1(II)))
          WVNXOUT(II) = RampMete*WINDX
          WVNYOUT(II) = RampMete*WINDY
#ifdef SWAN
!asey 101118: Added these lines for coupling winds to SWAN.
          IF(COUPWIND)THEN
            SWAN_WX2(II,2) = WINDX/WindMultiplier
            SWAN_WY2(II,2) = WINDY/WindMultiplier
          ENDIF
#endif
        ENDDO
      ENDIF
      
      RETURN
      END SUBROUTINE MET_FORCING
