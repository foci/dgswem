!**************************************************************************
!  mod history
!  v41.06mxxx - date - programmer - describe change 
!                    - mark change in code with  cinitials-mxxx
!
!  v41.11 - 09/14/01 - rl - from 41.09 - modified for global_here%NWS=-2
!  v41.09 - 06/30/01 - jw - from 41.08 - made minor mods as global_here%per vp version 41.05 
!  v41.08 - 06/22/01 - rl - from 41.07 - added 41.05m009 changes to global_here%HABSMIN
!                                        and global_here%ETA2
!  v41.07 - 04/09/01 - rl - from 41.06 - initialized global_here%PRN1(), global_here%PRN2() for global_here%NRS<>0
!**************************************************************************
!
        SUBROUTINE COLDSTART(s)
!
!**************************************************************************
!
!  COLD START PROGRAM SETUP ROUTINE 
!
!**************************************************************************
!
      USE SIZES
      USE GLOBAL
#ifdef HARM
      USE HARM
#endif
      USE WIND
#ifdef OWIWIND
      USE OWIWIND, ONLY : NWS12INIT, NWS12GET
#endif
!      USE NodalAttributes, ONLY : STARTDRY, GeoidOffset, LoadGeoidOffset
      IMPLICIT NONE

      type (sizes_type) :: s

      INTEGER i,j
!
      global_here%ITHS = 0

!...
!....SET AT REST INITIAL CONDITION OVER WHOLE DOMAIN
!....IF BOTTOM IS ABOVE THE GEIOD -> DRY NODE
!....IF BOTTOM IS INITIALLY BELOW THE GEIOD AND STARTDRY=-88888 -> DRY NODE
!...
        global_here%HABSMIN=0.8d0*global_here%H0
        DO I=1,global_here%NP
          global_here%UU1(I) =0.D0
          global_here%VV1(I) =0.D0
          global_here%UU2(I) =0.D0
          global_here%VV2(I) =0.D0
          global_here%ETA2(I)=0.D0
          global_here%NODEREP(I)=MAX0(global_here%NODEWETMIN,global_here%NODEDRYMIN)
          global_here%NNODECODE(I)=1
          IF(global_here%NOLIFA.EQ.2) THEN
            global_here%HTOT=global_here%DP(I)+global_here%ETA2(I)
            IF(global_here%HTOT.LE.global_here%H0) THEN
              global_here%NNODECODE(I)=0
              global_here%ETA2(I)=global_here%H0-global_here%DP(I)
!              ELSE
!              IF(STARTDRY(I).EQ.-88888) THEN
!                global_here%NNODECODE(I)=0
!                global_here%ETA2(I)=global_here%H0-global_here%DP(I)
!                ENDIF
            ENDIF
          ENDIF
          global_here%ETA1(I)=global_here%ETA2(I)
          global_here%ETAS(I)=0.D0
          global_here%CH1(I)=0.d0
          END DO

!...
!....INITIALIZE THE ELEVATION SPECIFIED BOUNDARY CONDITION IF IT REQUIRES THE USE
!....OF THE UNIT 19 FILE.
!...

        IF((global_here%NOPE.GT.0).AND.(global_here%NBFR.EQ.0)) THEN
          IF(global_here%NSCREEN.EQ.1.AND.S%MYPROC.EQ.0) WRITE(6,1112)
          WRITE(16,1112)
          IF(global_here%NSCREEN.EQ.1.AND.S%MYPROC.EQ.0) WRITE(6,1977)
          WRITE(16,1977)
 1977     FORMAT(/,1X,'ELEVATION SPECIFIED INFORMATION READ FROM UNIT ',&
           '19',/)
          OPEN(19,FILE=S%DIRNAME//'/'//'fort.19')
          READ(19,*) global_here%ETIMINC
          DO J=1,global_here%NETA
             READ(19,*) global_here%ESBIN1(J)
          END DO
          DO J=1,global_here%NETA
             READ(19,*) global_here%ESBIN2(J)
          END DO
          global_here%ETIME1 = global_here%STATIM*86400.D0
          global_here%ETIME2 = global_here%ETIME1 + global_here%ETIMINC
        ENDIF

!....INITIALIZE THE NORMAL FLOW BOUNDARY CONDITION

        DO I=1,global_here%NVEL
          global_here%QN2(I)=0.D0
          global_here%QN1(I)=0.D0
          global_here%QN0(I)=0.D0
        END DO

        IF((global_here%NFLUXF.EQ.1).AND.(global_here%NFFR.EQ.0)) THEN
          IF(global_here%NSCREEN.EQ.1.AND.S%MYPROC.EQ.0) WRITE(6,1112)
          WRITE(16,1112)
          IF(global_here%NSCREEN.EQ.1.AND.S%MYPROC.EQ.0) WRITE(6,1979)
          WRITE(16,1979)
 1979     FORMAT(/,1X,'NORMAL FLOW INFORMATION READ FROM UNIT 20',/)
          OPEN(20,FILE=S%DIRNAME//'/'//'fort.20')
          READ(20,*) global_here%FTIMINC
          DO J=1,global_here%NVEL
            global_here%QNIN1(J)=0.D0
            IF((global_here%LBCODEI(J).EQ.2).OR.(global_here%LBCODEI(J).EQ.12)&
                          .OR.(global_here%LBCODEI(J).EQ.22))&
                                          READ(20,*) global_here%QNIN1(J)
            END DO
          DO J=1,global_here%NVEL
            global_here%QNIN2(J)=0.D0
            IF((global_here%LBCODEI(J).EQ.2).OR.(global_here%LBCODEI(J).EQ.12)&
                          .OR.(global_here%LBCODEI(J).EQ.22))&
                                          READ(20,*) global_here%QNIN2(J)
            END DO
          global_here%QTIME1 = global_here%STATIM*86400.D0
          global_here%QTIME2 = global_here%QTIME1 + global_here%FTIMINC
          ENDIF

!...INPUT METEOROLOGICAL INFORMATION FROM UNIT 22 OR UNIT 200 SERIES
!....IF FLEET NUMERIC WIND DATA IS USED, FIND BEGINNING TIME IN FILE,
!....NOTE: CAN'T DEAL WITH WIND THAT STARTS AFTER global_here%WREFTIM!!!!!!!!!!!!
!....READ IN AND INTERPOLATE IN SPACE ONTO THE ADCIRC GRID THE
!....TIME LEVEL 1 AND LEVEL 2 WIND FIELDS

        IF(global_here%NWS.global_here%NE.0) THEN
          DO I=1,global_here%NP
            global_here%WSX1(I)=0.D0
            global_here%WSY1(I)=0.D0
            global_here%PR1(I) =0.D0
            global_here%WSX2(I)=0.D0
            global_here%WSY2(I)=0.D0
            global_here%PR2(I) =0.D0
            ENDDO

          IF(global_here%NSCREEN.EQ.1.AND.S%MYPROC.EQ.0) WRITE(6,1112)
          WRITE(16,1112)
          IF(global_here%NSCREEN.EQ.1.AND.S%MYPROC.EQ.0) WRITE(6,1980)
          WRITE(16,1980)
 1980     FORMAT(/,1X,'WIND (AND PRESSURE) INFORMATION READ.',/)
          ENDIF

        IF(global_here%NWS.EQ.1) THEN
          OPEN(22,FILE=S%DIRNAME//'/'//'fort.22')
          ENDIF

        IF(ABS(global_here%NWS).EQ.2) THEN
          OPEN(22,FILE=S%DIRNAME//'/'//'fort.22')
          READ(22,*) (global_here%NHG,global_here%WVNX1(I),global_here%WVNY1(I),global_here%PRN1(I),I=1,global_here%NP)
          READ(22,*) (global_here%NHG,global_here%WVNX2(I),global_here%WVNY2(I),global_here%PRN2(I),I=1,global_here%NP)
          global_here%WTIME1 = global_here%STATIM*86400.D0
          global_here%WTIME2 = global_here%WTIME1 + global_here%WTIMINC
          ENDIF

        IF(global_here%NWS.EQ.3) THEN
          OPEN(22,FILE=S%DIRNAME//'/'//'fort.22')
 2222     CALL NWS3GET(s, global_here%X, global_here%Y, global_here%SLAM, global_here%SFEA, global_here%WVNX2, global_here%WVNY2, global_here%IWTIME, global_here%IWYR,&
                  global_here%WTIMED, global_here%NP, global_here%NWLON, global_here%NWLAT, global_here%WLATMAX, global_here%WLONMIN,&
                  global_here%WLATINC, global_here%WLONINC, global_here%ICS, global_here%NSCREEN, global_here%ScreenUnit )
          IF(global_here%IWYR.global_here%NE.global_here%IREFYR) THEN
            global_here%IWTIMEP=global_here%IWTIME
            DO I=1,global_here%NP
              global_here%WVNX1(I)=global_here%WVNX2(I)
              global_here%WVNY1(I)=global_here%WVNY2(I)
              END DO
            GOTO 2222
            ENDIF
          IF(global_here%WTIMED.LE.global_here%WREFTIM) THEN
            global_here%IWTIMEP=global_here%IWTIME
            DO I=1,global_here%NP
              global_here%WVNX1(I)=global_here%WVNX2(I)
              global_here%WVNY1(I)=global_here%WVNY2(I)
              END DO
            GOTO 2222
            ENDIF
          IF(global_here%NSCREEN.EQ.1.AND.S%MYPROC.EQ.0) &
         WRITE(6,*)'FOUND WIND DATA AT TIME= ',global_here%IWTIMEP
          WRITE(16,*) 'FOUND WIND DATA AT TIME= ',global_here%IWTIMEP
          IF(global_here%NSCREEN.EQ.1.AND.S%MYPROC.EQ.0) &
         WRITE(6,*)'FOUND WIND DATA AT TIME= ',global_here%IWTIME
          WRITE(16,*) 'FOUND WIND DATA AT TIME= ',global_here%IWTIME
          global_here%WTIME2=global_here%WTIMED-global_here%WREFTIM                  !CAST INTO MODEL TIME REFRENCE
          global_here%WTIME1=global_here%WTIME2-global_here%WTIMINC
          ENDIF

        IF(ABS(global_here%NWS).EQ.4) THEN
          OPEN(22,FILE=S%DIRNAME//'/'//'fort.22')
          global_here%WTIME1 = global_here%STATIM*86400.D0
          global_here%WTIME2=global_here%WTIME1+global_here%WTIMINC
          CALL NWS4GET(global_here%WVNX1,global_here%WVNY1,global_here%PRN1,global_here%NP,global_here%RHOWAT0,global_here%G)
          CALL NWS4GET(global_here%WVNX2,global_here%WVNY2,global_here%PRN2,global_here%NP,global_here%RHOWAT0,global_here%G)
          ENDIF

        IF(ABS(global_here%NWS).EQ.5) THEN
          OPEN(22,FILE=S%DIRNAME//'/'//'fort.22')
          READ(22,*) (global_here%NHG,global_here%WVNX1(I),global_here%WVNY1(I),global_here%PRN1(I),I=1,global_here%NP)
          READ(22,*) (global_here%NHG,global_here%WVNX2(I),global_here%WVNY2(I),global_here%PRN2(I),I=1,global_here%NP)
          global_here%WTIME1 = global_here%STATIM*86400.D0
          global_here%WTIME2 = global_here%WTIME1 + global_here%WTIMINC
          ENDIF

        IF(global_here%NWS.EQ.6) THEN
          OPEN(22,FILE=S%DIRNAME//'/'//'fort.22')
          CALL NWS6GET(global_here%X,global_here%Y,global_here%SLAM,global_here%SFEA,global_here%WVNX1,global_here%WVNY1,global_here%PRN1,global_here%NP,global_here%NWLON,global_here%NWLAT,&
                 global_here%WLATMAX,global_here%WLONMIN,global_here%WLATINC,global_here%WLONINC,global_here%ICS,global_here%RHOWAT0,global_here%G)
          CALL NWS6GET(global_here%X,global_here%Y,global_here%SLAM,global_here%SFEA,global_here%WVNX2,global_here%WVNY2,global_here%PRN2,global_here%NP,global_here%NWLON,global_here%NWLAT,&
                 global_here%WLATMAX,global_here%WLONMIN,global_here%WLATINC,global_here%WLONINC,global_here%ICS,global_here%RHOWAT0,global_here%G)
          global_here%WTIME1 = global_here%STATIM*86400.D0
          global_here%WTIME2 = global_here%WTIME1 + global_here%WTIMINC
          ENDIF

        IF(global_here%NWS.EQ.10) THEN
          global_here%WTIME1=global_here%STATIM*86400.D0
          global_here%WTIME2=global_here%WTIME1+global_here%WTIMINC
          global_here%NWSGGWI=-1
          CALL NWS10GET(s,global_here%NWSGGWI,global_here%SLAM,global_here%SFEA,global_here%WVNX1,global_here%WVNY1,global_here%PRN1,global_here%NP,global_here%RHOWAT0,global_here%G,&
                  global_here%NWLON,global_here%NWLAT,global_here%WTIMINC) !JUST COMPUTE INTERPOLATING FACTORS
          global_here%NWSGGWI=1
          CALL NWS10GET(s,global_here%NWSGGWI,global_here%SLAM,global_here%SFEA,global_here%WVNX2,global_here%WVNY2,global_here%PRN2,global_here%NP,global_here%RHOWAT0,global_here%G,&
                  global_here%NWLON,global_here%NWLAT,global_here%WTIMINC) !NOW INTERPOLATE 1st WIND FIELD
          ENDIF

        IF(global_here%NWS.EQ.11) THEN
          global_here%WTIME1=global_here%STATIM*86400.D0
          global_here%WTIME2=global_here%WTIME1+global_here%WTIMINC
          global_here%NWSEGWI=0
          global_here%IDSETFLG=0
          IF(global_here%NSCREEN.EQ.1.AND.S%MYPROC.EQ.0) WRITE(6,1197)
          WRITE(16,1197)
 1197     FORMAT(/,1X,'THE E29 MET GRID INTERPOLATING FACTORS ARE ',&
                'BEING COMPUTED ')
          CALL NWS11GET(s,global_here%NWSEGWI,global_here%IDSETFLG,global_here%SLAM,global_here%SFEA,global_here%WVNX1,global_here%WVNY1,global_here%PRN1,global_here%NP,&
                  global_here%RHOWAT0,global_here%G)  !JUST COMPUTE INTERPOLATING FACTORS
          IF(global_here%NSCREEN.EQ.1.AND.S%MYPROC.EQ.0) WRITE(6,1198)
          WRITE(16,1198)
 1198     FORMAT(1X,'FINISHED COMPUTING E29 INTERPOLATING FACTORS',/)
          global_here%NWSEGWI=1
          global_here%IDSETFLG=1
          CALL NWS11GET(s,global_here%NWSEGWI,global_here%IDSETFLG,global_here%SLAM,global_here%SFEA,global_here%WVNX2,global_here%WVNY2,global_here%PRN2,global_here%NP,&
                  global_here%RHOWAT0,global_here%G) !NOW INTERPOLATE 1st WIND FIELD
        ENDIF
        
!ek added global_here%nws = 12 for owi winds

#ifdef OWIWIND
        IF(ABS(global_here%NWS).EQ.12) THEN
          CALL NWS12INIT( global_here%WVNX1, global_here%WVNY1, global_here%PRN1, global_here%NP, global_here%RHOWAT0, global_here%G )
          CALL NWS12GET(  global_here%WVNX1, global_here%WVNY1, global_here%PRN1, global_here%NP, global_here%RHOWAT0, global_here%G )
          CALL NWS12GET(  global_here%WVNX2, global_here%WVNY2, global_here%PRN2, global_here%NP, global_here%RHOWAT0, global_here%G )
          global_here%WTIME1 = global_here%STATIM*86400.D0
          global_here%WTIME2 = global_here%WTIME1 + global_here%WTIMINC
        ENDIF
#endif

!...INPUT RADIATION STRESS INFORMATION FROM UNIT 23
!....READ IN THE TIME LEVEL 1 AND LEVEL 2 FIELDS

        IF(global_here%NRS.EQ.1) THEN
          IF(global_here%NWS.EQ.0) THEN
            DO I=1,global_here%NP
              global_here%WSX1(I)=0.D0
              global_here%WSY1(I)=0.D0
              global_here%WSX2(I)=0.D0
              global_here%WSY2(I)=0.D0
              global_here%PRN1(I)=0.D0    !need to be initialized
              global_here%PRN2(I)=0.D0    !even if not used
              ENDDO
            ENDIF
          OPEN(23,FILE=S%DIRNAME//'/'//'fort.23')
          global_here%RSTIME1 = global_here%STATIM*86400.D0
          global_here%RSTIME2 = global_here%RSTIME1+global_here%RSTIMINC
          IF (global_here%FRW.EQ.0) THEN
            CALL RSGET(global_here%RSNX1,global_here%RSNY1,global_here%NP)
            CALL RSGET(global_here%RSNX2,global_here%RSNY2,global_here%NP)
          ENDIF
          IF (global_here%FRW.EQ.1) THEN
!            CALL RSGET_MORE(global_here%RSNX1,global_here%RSNY1,global_here%WAVE_H1,global_here%WAVE_T1,global_here%WAVE_A1,
!     &                                                  global_here%WAVE_D1,global_here%NP)
!            CALL RSGET_MORE(global_here%RSNX2,global_here%RSNY2,global_here%WAVE_H2,global_here%WAVE_T2,global_here%WAVE_A2,
!     &                                                  global_here%WAVE_D2,global_here%NP)
            
          ENDIF
#ifdef SWAN   
!asey 101118: Added this call to initialize the radiation stress gradients.
         IF(global_here%NRS.EQ.3) THEN
           IF(.NOT.ALLOCATED(global_here%RSNX1)) ALLOCATE(global_here%RSNX1(1:global_here%NP))
           IF(.NOT.ALLOCATED(global_here%RSNX2)) ALLOCATE(global_here%RSNX2(1:global_here%NP))
           IF(.NOT.ALLOCATED(global_here%RSNY1)) ALLOCATE(global_here%RSNY1(1:global_here%NP))
           IF(.NOT.ALLOCATED(global_here%RSNY2)) ALLOCATE(global_here%RSNY2(1:global_here%NP))
           DO I=1,global_here%NP
             global_here%RSNX1(I) = 0.D0
             global_here%RSNX2(I) = 0.D0
             global_here%RSNY1(I) = 0.D0
             global_here%RSNY2(I) = 0.D0
           ENDDO
         ENDIF
#endif
          IF(global_here%NSCREEN.EQ.1.AND.S%MYPROC.EQ.0) WRITE(6,1112)
          WRITE(16,1112)
          IF(global_here%NSCREEN.EQ.1.AND.S%MYPROC.EQ.0) WRITE(6,1981)
          WRITE(16,1981)
 1981     FORMAT(/,1X,'RADIATION STRESS INFORMATION READ.',/)
          ENDIF


!...
!...LINES TO USE TIDAL POTENTIAL FORCING
!...
       if (s%CTIP) then
          DO I=1,global_here%NP
             global_here%TIP2(I)=0.0
          END DO
       endif

!WET...
!WET...THE FOLLOWING LINES ARE FOR WETTING AND DRYING
!WET...Dry any landlocked nodes by checking that they are connected to at
!WET...least 1 functioning element.
!WET...
        IF(global_here%NOLIFA.EQ.2) THEN
          DO I=1,global_here%NP
            global_here%MJU(I)=0
            ENDDO
          DO I=1,global_here%NE
            global_here%NM1=global_here%NM(I,1)
            global_here%NM2=global_here%NM(I,2)
            global_here%NM3=global_here%NM(I,3)
            global_here%NC1=global_here%NNODECODE(global_here%NM1)
            global_here%NC2=global_here%NNODECODE(global_here%NM2)
            global_here%NC3=global_here%NNODECODE(global_here%NM3)
            global_here%NCELE=global_here%NC1*global_here%NC2*global_here%NC3
            global_here%MJU(global_here%NM1)=global_here%MJU(global_here%NM1)+global_here%NCELE
            global_here%MJU(global_here%NM2)=global_here%MJU(global_here%NM2)+global_here%NCELE
            global_here%MJU(global_here%NM3)=global_here%MJU(global_here%NM3)+global_here%NCELE
          ENDDO
          DO I=1,global_here%NP
            IF((global_here%NNODECODE(I).EQ.1).AND.(global_here%MJU(I).EQ.0)) THEN
              global_here%NNODECODE(I)=0
              IF(global_here%NSCREEN.EQ.1.AND.S%MYPROC.EQ.0) WRITE(*,9883) I
              WRITE(16,9883) I
              ENDIF
          ENDDO
        ENDIF
!...
!......INITIALIZE 3D SOLUTION
!...

!...LINES TO RUN THE CODE IN 3D VS MODE.

      if (s%C3DVS) then
!       CALL VSSTUP(global_here%DT,global_here%STATIM,NBYTE,global_here%RUNDES,global_here%RUNID,global_here%AGRID,global_here%NT)
      endif

!...LINES TO RUN THE CODE IN 3D DSS MODE

      if (S%C3DDSS) then
!       CALL DSSSTUP(global_here%DT,global_here%STATIM,NBYTE,global_here%RUNDES,global_here%RUNID,global_here%AGRID,global_here%NT)
      endif


!...
!....INITILIZE ELEVATION STATION SPOOL COUNTER
!....OPEN ELEVATION STATION OUTPUT FILE
!....WRITE OUT HEADER INFORMATION INCLUDING global_here%NTRSPE (NO. OF DATA PTS. AT EACH
!....ELEVATION STATION), global_here%NSTAE, global_here%DT*global_here%NSPOOLE, global_here%NSPOOLE, IRTYPE
!...
        global_here%NSCOUE=0
        global_here%IESTP=0

 3220   FORMAT(1X,A32,2X,A24,2X,A24)
 3645   FORMAT(1X,I10,1X,I10,1X,E15.7,1X,I5,1X,I5)

        IF(ABS(global_here%NOUTE).EQ.1) THEN
          OPEN(61,FILE=S%DIRNAME//'/'//'fort.61')
          WRITE(61,3220) global_here%RUNDES,global_here%RUNID,global_here%AGRID
          WRITE(61,3645) global_here%NTRSPE,global_here%NSTAE,global_here%DTDP*global_here%NSPOOLE,global_here%NSPOOLE,1
          global_here%IESTP=2
          ENDIF

          IF(ABS(global_here%NOUTE).EQ.2) THEN
             OPEN(61,FILE=S%DIRNAME//'/'//'fort.61',&
                  ACCESS='DIRECT',RECL=NBYTE)
             IF(NBYTE.EQ.4) THEN
                DO I=1,8
                   WRITE(61,REC=global_here%IESTP+I) global_here%RDES4(I)
                ENDDO
                global_here%IESTP=global_here%IESTP+8
                DO I=1,6
                   WRITE(61,REC=global_here%IESTP+I) global_here%RID4(I)
                ENDDO
                global_here%IESTP=global_here%IESTP+6
                DO I=1,6
                   WRITE(61,REC=global_here%IESTP+I) global_here%AID4(I)
                ENDDO
                global_here%IESTP=global_here%IESTP+6
             ENDIF
             IF(NBYTE.EQ.8) THEN
                DO I=1,4
                   WRITE(61,REC=global_here%IESTP+I) global_here%RDES8(I)
                ENDDO
                global_here%IESTP=global_here%IESTP+4
                DO I=1,3
                   WRITE(61,REC=global_here%IESTP+I) global_here%RID8(I)
                ENDDO
                global_here%IESTP=global_here%IESTP+3
                DO I=1,3
                   WRITE(61,REC=global_here%IESTP+I) global_here%AID8(I)
                ENDDO
                global_here%IESTP=global_here%IESTP+3
             ENDIF
             WRITE(61,REC=global_here%IESTP+1) global_here%NTRSPE
             WRITE(61,REC=global_here%IESTP+2) global_here%NSTAE
             WRITE(61,REC=global_here%IESTP+3) global_here%DT*global_here%NSPOOLE
             WRITE(61,REC=global_here%IESTP+4) global_here%NSPOOLE
             WRITE(61,REC=global_here%IESTP+5) 1
             global_here%IESTP=global_here%IESTP+5
             CLOSE(61)                    ! DO THIS TO FLUSH THE WRITE BUFFER
             OPEN(61,FILE=S%DIRNAME//'/'//'fort.61',&
                  ACCESS='DIRECT',RECL=NBYTE)
          ENDIF

!...
!....INITILIZE VELOCITY STATION SPOOL COUNTER
!....OPEN VELOCITY STATION OUTPUT FILE
!....WRITE OUT HEADER INFORMATION INCLUDING global_here%NTRSPV (NO. OF DATA PTS. AT EACH
!....VELOCITY STATION), global_here%NSTAV, global_here%DT*global_here%NSPOOLV, global_here%NSPOOLV, IRTYPE
!...
        global_here%NSCOUV=0
        global_here%IVSTP=0

        IF(ABS(global_here%NOUTV).EQ.1) THEN
          OPEN(62,FILE=S%DIRNAME//'/'//'fort.62')
          WRITE(62,3220) global_here%RUNDES,global_here%RUNID,global_here%AGRID
          WRITE(62,3645) global_here%NTRSPV,global_here%NSTAV,global_here%DTDP*global_here%NSPOOLV,global_here%NSPOOLV,2
          global_here%IVSTP=2
          ENDIF

        IF(ABS(global_here%NOUTV).EQ.2) THEN
          OPEN(62,FILE=S%DIRNAME//'/'//'fort.62',&
          ACCESS='DIRECT',RECL=NBYTE)
          IF(NBYTE.EQ.4) THEN
            DO I=1,8
              WRITE(62,REC=global_here%IVSTP+I) global_here%RDES4(I)
              ENDDO
            global_here%IVSTP=global_here%IVSTP+8
            DO I=1,6
              WRITE(62,REC=global_here%IVSTP+I) global_here%RID4(I)
              ENDDO
            global_here%IVSTP=global_here%IVSTP+6
            DO I=1,6
              WRITE(62,REC=global_here%IVSTP+I) global_here%AID4(I)
              ENDDO
            global_here%IVSTP=global_here%IVSTP+6
            ENDIF
          IF(NBYTE.EQ.8) THEN
            DO I=1,4
              WRITE(62,REC=global_here%IVSTP+I) global_here%RDES8(I)
              ENDDO
            global_here%IVSTP=global_here%IVSTP+4
            DO I=1,3
              WRITE(62,REC=global_here%IVSTP+I) global_here%RID8(I)
              ENDDO
            global_here%IVSTP=global_here%IVSTP+3
            DO I=1,3
              WRITE(62,REC=global_here%IVSTP+I) global_here%AID8(I)
              ENDDO
            global_here%IVSTP=global_here%IVSTP+3
            ENDIF
          WRITE(62,REC=global_here%IVSTP+1) global_here%NTRSPV
          WRITE(62,REC=global_here%IVSTP+2) global_here%NSTAV
          WRITE(62,REC=global_here%IVSTP+3) global_here%DT*global_here%NSPOOLV
          WRITE(62,REC=global_here%IVSTP+4) global_here%NSPOOLV
          WRITE(62,REC=global_here%IVSTP+5) 2
          global_here%IVSTP=global_here%IVSTP+5
          CLOSE(62)                    ! DO THIS TO FLUSH THE WRITE BUFFER
          OPEN(62,FILE=S%DIRNAME//'/'//'fort.62',&
         ACCESS='DIRECT',RECL=NBYTE)
          ENDIF

!...
!....INITILIZE CONCENTRATION STATION SPOOL COUNTER
!....OPEN ELEVATION STATION OUTPUT FILE
!....WRITE OUT HEADER INFORMATION INCLUDING global_here%NTRSPC (NO. OF DATA PTS. AT EACH
!....CONCENTRATION STATION), global_here%NSTAC, global_here%DT*global_here%NSPOOLC, global_here%NSPOOLC, IRTYPE
!...
        global_here%NSCOUC=0
        global_here%ICSTP=0

        IF(ABS(global_here%NOUTC).EQ.1) THEN
          OPEN(81,FILE=S%DIRNAME//'/'//'fort.81')
          WRITE(81,3220) global_here%RUNDES,global_here%RUNID,global_here%AGRID
          WRITE(81,3645) global_here%NTRSPC,global_here%NSTAC,global_here%DTDP*global_here%NSPOOLC,global_here%NSPOOLC,1
          global_here%ICSTP=2
          ENDIF

        IF(ABS(global_here%NOUTC).EQ.2) THEN
          OPEN(81,FILE=S%DIRNAME//'/'//'fort.81',&
          ACCESS='DIRECT',RECL=NBYTE)
          IF(NBYTE.EQ.4) THEN
            DO I=1,8
              WRITE(81,REC=global_here%ICSTP+I) global_here%RDES4(I)
              ENDDO
            global_here%ICSTP=global_here%ICSTP+8
            DO I=1,6
              WRITE(81,REC=global_here%ICSTP+I) global_here%RID4(I)
              ENDDO
            global_here%ICSTP=global_here%ICSTP+6
            DO I=1,6
              WRITE(81,REC=global_here%ICSTP+I) global_here%AID4(I)
              ENDDO
            global_here%ICSTP=global_here%ICSTP+6
            ENDIF
          IF(NBYTE.EQ.8) THEN
            DO I=1,4
              WRITE(81,REC=global_here%ICSTP+I) global_here%RDES8(I)
              ENDDO
            global_here%ICSTP=global_here%ICSTP+4
            DO I=1,3
              WRITE(81,REC=global_here%ICSTP+I) global_here%RID8(I)
              ENDDO
            global_here%ICSTP=global_here%ICSTP+3
            DO I=1,3
              WRITE(81,REC=global_here%ICSTP+I) global_here%AID8(I)
              ENDDO
            global_here%ICSTP=global_here%ICSTP+3
            ENDIF
          WRITE(81,REC=global_here%ICSTP+1) global_here%NTRSPC
          WRITE(81,REC=global_here%ICSTP+2) global_here%NSTAC
          WRITE(81,REC=global_here%ICSTP+3) global_here%DT*global_here%NSPOOLC
          WRITE(81,REC=global_here%ICSTP+4) global_here%NSPOOLC
          WRITE(81,REC=global_here%ICSTP+5) 1
          global_here%ICSTP=global_here%ICSTP+5
          CLOSE(81)                    ! DO THIS TO FLUSH THE WRITE BUFFER
          OPEN(81,FILE=S%DIRNAME//'/'//'fort.81',&
         ACCESS='DIRECT',RECL=NBYTE)
          ENDIF


!...
!....INITILIZE BATHYMETRY STATION SPOOL COUNTER
!....OPEN ELEVATION STATION OUTPUT FILE
!....WRITE OUT HEADER INFORMATION INCLUDING global_here%NTRSPE (NO. OF DATA PTS. AT EACH
!....ELEVATION STATION), global_here%NSTAE, global_here%DT*global_here%NSPOOLE, global_here%NSPOOLE, IRTYPE
!...
        global_here%NSCOUE=0
        global_here%IESTP=0

        IF(ABS(global_here%NOUTE).EQ.1) THEN
          OPEN(82,FILE=S%DIRNAME//'/'//'fort.82')
          WRITE(82,3220) global_here%RUNDES,global_here%RUNID,global_here%AGRID
          WRITE(82,3645) global_here%NTRSPE,global_here%NSTAE,global_here%DTDP*global_here%NSPOOLE,global_here%NSPOOLE,1
          global_here%IESTP=2
          ENDIF

        IF(ABS(global_here%NOUTE).EQ.2) THEN
          OPEN(82,FILE=S%DIRNAME//'/'//'fort.82',&
          ACCESS='DIRECT',RECL=NBYTE)
          IF(NBYTE.EQ.4) THEN
            DO I=1,8
              WRITE(82,REC=global_here%IESTP+I) global_here%RDES4(I)
              ENDDO
            global_here%IESTP=global_here%IESTP+8
            DO I=1,6
              WRITE(82,REC=global_here%IESTP+I) global_here%RID4(I)
              ENDDO
            global_here%IESTP=global_here%IESTP+6
            DO I=1,6
              WRITE(82,REC=global_here%IESTP+I) global_here%AID4(I)
              ENDDO
            global_here%IESTP=global_here%IESTP+6
            ENDIF
          IF(NBYTE.EQ.8) THEN
            DO I=1,4
              WRITE(82,REC=global_here%IESTP+I) global_here%RDES8(I)
              ENDDO
            global_here%IESTP=global_here%IESTP+4
            DO I=1,3
              WRITE(82,REC=global_here%IESTP+I) global_here%RID8(I)
              ENDDO
            global_here%IESTP=global_here%IESTP+3
            DO I=1,3
              WRITE(82,REC=global_here%IESTP+I) global_here%AID8(I)
              ENDDO
            global_here%IESTP=global_here%IESTP+3
            ENDIF
          WRITE(82,REC=global_here%IESTP+1) global_here%NTRSPE
          WRITE(82,REC=global_here%IESTP+2) global_here%NSTAE
          WRITE(82,REC=global_here%IESTP+3) global_here%DT*global_here%NSPOOLE
          WRITE(82,REC=global_here%IESTP+4) global_here%NSPOOLE
          WRITE(82,REC=global_here%IESTP+5) 1
          global_here%IESTP=global_here%IESTP+5
          CLOSE(82)                    ! DO THIS TO FLUSH THE WRITE BUFFER
          OPEN(82,FILE=S%DIRNAME//'/'//'fort.82',&
         ACCESS='DIRECT',RECL=NBYTE)
          ENDIF


!...
!....INITILIZE METEOROLOGICAL STATION SPOOL COUNTERS
!....OPEN METEOROLOGICAL STATION OUTPUT FILES
!....WRITE OUT HEADER INFORMATION INCLUDING global_here%NTRSPM (NO. OF DATA PTS. AT EACH
!....METEOROLOGICAL STATION), global_here%NSTAM, global_here%DT*global_here%NSPOOLM, global_here%NSPOOLM, IRTYPE
!...
        global_here%NSCOUM=0
        global_here%IPSTP=0
        global_here%IWSTP=0

        IF(ABS(global_here%NOUTM).EQ.1) THEN
          OPEN(71,FILE=S%DIRNAME//'/'//'fort.71')
          WRITE(71,3220) global_here%RUNDES,global_here%RUNID,global_here%AGRID
          WRITE(71,3645) global_here%NTRSPM,global_here%NSTAM,global_here%DTDP*global_here%NSPOOLM,global_here%NSPOOLM,1
          global_here%IPSTP=2
          OPEN(72,FILE=S%DIRNAME//'/'//'fort.72')
          WRITE(72,3220) global_here%RUNDES,global_here%RUNID,global_here%AGRID
          WRITE(72,3645) global_here%NTRSPM,global_here%NSTAM,global_here%DTDP*global_here%NSPOOLM,global_here%NSPOOLM,2
          global_here%IWSTP=2
          ENDIF

        IF(ABS(global_here%NOUTM).EQ.2) THEN
          OPEN(71,FILE=S%DIRNAME//'/'//'fort.71',&
          ACCESS='DIRECT',RECL=NBYTE)
          OPEN(72,FILE=S%DIRNAME//'/'//'fort.72',&
          ACCESS='DIRECT',RECL=NBYTE)
          IF(NBYTE.EQ.4) THEN
            DO I=1,8
              WRITE(71,REC=global_here%IPSTP+I) global_here%RDES4(I)
              WRITE(72,REC=global_here%IWSTP+I) global_here%RDES4(I)
              ENDDO
            global_here%IPSTP=global_here%IPSTP+8
            global_here%IWSTP=global_here%IWSTP+8
            DO I=1,6
              WRITE(71,REC=global_here%IPSTP+I) global_here%RID4(I)
              WRITE(72,REC=global_here%IWSTP+I) global_here%RID4(I)
              ENDDO
            global_here%IPSTP=global_here%IPSTP+6
            global_here%IWSTP=global_here%IWSTP+6
            DO I=1,6
              WRITE(71,REC=global_here%IPSTP+I) global_here%AID4(I)
              WRITE(72,REC=global_here%IWSTP+I) global_here%AID4(I)
              ENDDO
            global_here%IPSTP=global_here%IPSTP+6
            global_here%IWSTP=global_here%IWSTP+6
            ENDIF
          IF(NBYTE.EQ.8) THEN
            DO I=1,4
              WRITE(71,REC=global_here%IPSTP+I) global_here%RDES8(I)
              WRITE(72,REC=global_here%IWSTP+I) global_here%RDES8(I)
              ENDDO
            global_here%IPSTP=global_here%IPSTP+4
            global_here%IWSTP=global_here%IWSTP+4
            DO I=1,3
              WRITE(71,REC=global_here%IPSTP+I) global_here%RID8(I)
              WRITE(72,REC=global_here%IWSTP+I) global_here%RID8(I)
              ENDDO
            global_here%IPSTP=global_here%IPSTP+3
            global_here%IWSTP=global_here%IWSTP+3
            DO I=1,3
              WRITE(71,REC=global_here%IPSTP+I) global_here%AID8(I)
              WRITE(72,REC=global_here%IWSTP+I) global_here%AID8(I)
              ENDDO
            global_here%IPSTP=global_here%IPSTP+3
            global_here%IWSTP=global_here%IWSTP+3
            ENDIF
          WRITE(71,REC=global_here%IPSTP+1) global_here%NTRSPM
          WRITE(71,REC=global_here%IPSTP+2) global_here%NSTAM
          WRITE(71,REC=global_here%IPSTP+3) global_here%DT*global_here%NSPOOLM
          WRITE(71,REC=global_here%IPSTP+4) global_here%NSPOOLM
          WRITE(71,REC=global_here%IPSTP+5) 1
          WRITE(72,REC=global_here%IWSTP+1) global_here%NTRSPM
          WRITE(72,REC=global_here%IWSTP+2) global_here%NSTAM
          WRITE(72,REC=global_here%IWSTP+3) global_here%DT*global_here%NSPOOLM
          WRITE(72,REC=global_here%IWSTP+4) global_here%NSPOOLM
          WRITE(72,REC=global_here%IWSTP+5) 2
          global_here%IPSTP=global_here%IPSTP+5
          global_here%IWSTP=global_here%IWSTP+5
          CLOSE(71)                    ! DO THIS TO FLUSH THE WRITE BUFFER
          CLOSE(72)                    ! DO THIS TO FLUSH THE WRITE BUFFER
          OPEN(71,FILE=S%DIRNAME//'/'//'fort.71',&
         ACCESS='DIRECT',RECL=NBYTE)
          OPEN(72,FILE=S%DIRNAME//'/'//'fort.72',&
         ACCESS='DIRECT',RECL=NBYTE)
          ENDIF

!...
!....INITILIZE GLOBAL ELEVATION SPOOL COUNTER
!....OPEN GLOBAL ELEVATION OUTPUT FILE
!....WRITE OUT HEADER INFORMATION INCLUDING global_here%NDSETSE
!....(NO. OF GLOBAL ELEVATION DATA SETS TO BE SPOOLED),
!....global_here%NP, global_here%DT*global_here%NSPOOLGE, global_here%NSPOOLGE, IRTYPE
!...
        global_here%NSCOUGE=0
        global_here%IGEP=0

        IF(ABS(global_here%NOUTGE).EQ.1) THEN
          OPEN(63,FILE=S%DIRNAME//'/'//'fort.63')
          WRITE(63,3220) global_here%RUNDES,global_here%RUNID,global_here%AGRID
          WRITE(63,3645) global_here%NDSETSE,global_here%NE,global_here%DTDP*global_here%NSPOOLGE,global_here%NSPOOLGE,1
          global_here%IGEP=2
          ENDIF

        IF(ABS(global_here%NOUTGE).EQ.1) THEN
          OPEN(88,FILE=S%DIRNAME//'/'//'fort.88')
          WRITE(88,3220) global_here%RUNDES,global_here%RUNID,global_here%AGRID
          WRITE(88,3645) global_here%NDSETSE,global_here%NE,global_here%DTDP*global_here%NSPOOLGE,global_here%NSPOOLGE,1
          global_here%IGEP=2
          ENDIF

        IF(ABS(global_here%NOUTGE).EQ.1) THEN
          OPEN(89,FILE=S%DIRNAME//'/'//'fort.89')
          WRITE(89,3220) global_here%RUNDES,global_here%RUNID,global_here%AGRID
          WRITE(89,3645) global_here%NDSETSE,global_here%NE,global_here%DTDP*global_here%NSPOOLGE,global_here%NSPOOLGE,1
          global_here%IGEP=2
          ENDIF

        IF(ABS(global_here%NOUTGE).EQ.1) THEN
          OPEN(4441,FILE=S%DIRNAME//'/'//'fort.4l')
          WRITE(4441,3220) global_here%RUNDES,global_here%RUNID,global_here%AGRID
          WRITE(4441,3645) global_here%NDSETSE,global_here%NE,global_here%DTDP*global_here%NSPOOLGE,global_here%NSPOOLGE,1
          global_here%IGEP=2
          ENDIF

        IF(ABS(global_here%NOUTGE).EQ.2) THEN
          OPEN(63,FILE=S%DIRNAME//'/'//'fort.63',&
          ACCESS='DIRECT',RECL=NBYTE)
          IF(NBYTE.EQ.4) THEN
            DO I=1,8
              WRITE(63,REC=global_here%IGEP+I) global_here%RDES4(I)
              ENDDO
            global_here%IGEP=global_here%IGEP+8
            DO I=1,6
              WRITE(63,REC=global_here%IGEP+I) global_here%RID4(I)
              ENDDO
            global_here%IGEP=global_here%IGEP+6
            DO I=1,6
              WRITE(63,REC=global_here%IGEP+I) global_here%AID4(I)
              ENDDO
            global_here%IGEP=global_here%IGEP+6
            ENDIF
          IF(NBYTE.EQ.8) THEN
            DO I=1,4
              WRITE(63,REC=global_here%IGEP+I) global_here%RDES8(I)
              ENDDO
            global_here%IGEP=global_here%IGEP+4
            DO I=1,3
              WRITE(63,REC=global_here%IGEP+I) global_here%RID8(I)
              ENDDO
            global_here%IGEP=global_here%IGEP+3
            DO I=1,3
              WRITE(63,REC=global_here%IGEP+I) global_here%AID8(I)
              ENDDO
            global_here%IGEP=global_here%IGEP+3
            ENDIF
          WRITE(63,REC=global_here%IGEP+1) global_here%NDSETSE
          WRITE(63,REC=global_here%IGEP+2) global_here%NE
          WRITE(63,REC=global_here%IGEP+3) global_here%DT*global_here%NSPOOLGE
          WRITE(63,REC=global_here%IGEP+4) global_here%NSPOOLGE
          WRITE(63,REC=global_here%IGEP+5) 1
          global_here%IGEP=global_here%IGEP+5
          CLOSE(63)                    ! DO THIS TO FLUSH THE WRITE BUFFER
          OPEN(63,FILE=S%DIRNAME//'/'//'fort.63',&
         ACCESS='DIRECT',RECL=NBYTE)
          ENDIF

!...
!....INITILIZE GLOBAL VELOCITY SPOOL COUNTER
!....OPEN GLOBAL VELOCITY OUTPUT FILE
!....WRITE OUT HEADER INFORMATION INCLUDING global_here%NDSETSV
!....(NO. OF GLOBAL VELOCITY DATA SETS TO BE SPOOLED),
!....global_here%NP, global_here%DT*global_here%NSPOOLGV, global_here%NSPOOLGV, IRTYPE
!...
        global_here%NSCOUGV=0
        global_here%IGVP=0

        IF(ABS(global_here%NOUTGV).EQ.1) THEN
          OPEN(64,FILE=S%DIRNAME//'/'//'fort.64')
          WRITE(64,3220) global_here%RUNDES,global_here%RUNID,global_here%AGRID
          WRITE(64,3645) global_here%NDSETSV,global_here%NE,global_here%DTDP*global_here%NSPOOLGV,global_here%NSPOOLGV,2
          global_here%IGVP=2
          ENDIF

        IF(ABS(global_here%NOUTGV).EQ.2) THEN
          OPEN(64,FILE=S%DIRNAME//'/'//'fort.64',&
          ACCESS='DIRECT',RECL=NBYTE)
          IF(NBYTE.EQ.4) THEN
            DO I=1,8
              WRITE(64,REC=global_here%IGVP+I) global_here%RDES4(I)
              ENDDO
            global_here%IGVP=global_here%IGVP+8
            DO I=1,6
              WRITE(64,REC=global_here%IGVP+I) global_here%RID4(I)
              ENDDO
            global_here%IGVP=global_here%IGVP+6
            DO I=1,6
              WRITE(64,REC=global_here%IGVP+I) global_here%AID4(I)
              ENDDO
            global_here%IGVP=global_here%IGVP+6
            ENDIF
          IF(NBYTE.EQ.8) THEN
            DO I=1,4
              WRITE(64,REC=global_here%IGVP+I) global_here%RDES8(I)
              ENDDO
            global_here%IGVP=global_here%IGVP+4
            DO I=1,3
              WRITE(64,REC=global_here%IGVP+I) global_here%RID8(I)
              ENDDO
            global_here%IGVP=global_here%IGVP+3
            DO I=1,3
              WRITE(64,REC=global_here%IGVP+I) global_here%AID8(I)
              ENDDO
            global_here%IGVP=global_here%IGVP+3
            ENDIF
          WRITE(64,REC=global_here%IGVP+1) global_here%NDSETSV
          WRITE(64,REC=global_here%IGVP+2) global_here%NE
          WRITE(64,REC=global_here%IGVP+3) global_here%DT*global_here%NSPOOLGV
          WRITE(64,REC=global_here%IGVP+4) global_here%NSPOOLGV
          WRITE(64,REC=global_here%IGVP+5) 2
          global_here%IGVP=global_here%IGVP+5
          CLOSE(64)                    ! DO THIS TO FLUSH THE WRITE BUFFER
          OPEN(64,FILE=S%DIRNAME//'/'//'fort.64',&
         ACCESS='DIRECT',RECL=NBYTE)
          ENDIF

!...
!....INITILIZE GLOBAL WIND and pressure SPOOL COUNTER
!....OPEN GLOBAL WIND and pressure OUTPUT FILEs
!....WRITE OUT HEADER INFORMATION INCLUDING global_here%NDSETSW
!....(NO. OF GLOBAL WIND DATA SETS TO BE SPOOLED),
!....global_here%NP, global_here%DT*global_here%NSPOOLGW, global_here%NSPOOLGW, IRTYPE
!...
        global_here%NSCOUGW=0
        global_here%IGWP=0
        global_here%igpp=0

        IF(ABS(global_here%NOUTGW).EQ.1) THEN
          open(73,file=s%DIRNAME//'/'//'fort.73')
          write(73,3220) global_here%rundes,global_here%runid,global_here%agrid
          write(73,3645) global_here%ndsetsw,global_here%np,global_here%dtdp*global_here%nspoolgw,global_here%nspoolgw,1
          global_here%igpp=2
          OPEN(74,FILE=S%DIRNAME//'/'//'fort.74')
          WRITE(74,3220) global_here%RUNDES,global_here%RUNID,global_here%AGRID
          WRITE(74,3645) global_here%NDSETSW,global_here%NP,global_here%DTDP*global_here%NSPOOLGW,global_here%NSPOOLGW,2
          global_here%IGWP=2
          ENDIF

        IF(ABS(global_here%NOUTGW).EQ.2) THEN
          open(73,file=s%DIRNAME//'/'//'fort.73',&
          access='direct',recl=nbyte)
          OPEN(74,FILE=S%DIRNAME//'/'//'fort.74',&
          ACCESS='DIRECT',RECL=NBYTE)
          IF(NBYTE.EQ.4) THEN
            DO I=1,8
              write(73,rec=global_here%igpp+i) global_here%rdes4(i)
              WRITE(74,REC=global_here%IGWP+I) global_here%RDES4(I)
              ENDDO
            global_here%igpp=global_here%igpp+8
            global_here%IGWP=global_here%IGWP+8
            DO I=1,6
              write(73,rec=global_here%igpp+i) global_here%rid4(i)
              WRITE(74,REC=global_here%IGWP+I) global_here%RID4(I)
              ENDDO
            global_here%igpp=global_here%igpp+6
            global_here%IGWP=global_here%IGWP+6
            DO I=1,6
              write(73,rec=global_here%igpp+i) global_here%aid4(i)
              WRITE(74,REC=global_here%IGWP+I) global_here%AID4(I)
              ENDDO
            global_here%igpp=global_here%igpp+6
            global_here%IGWP=global_here%IGWP+6
            ENDIF
          IF(NBYTE.EQ.8) THEN
            DO I=1,4
              write(73,rec=global_here%igpp+i) global_here%rdes8(i)
              WRITE(74,REC=global_here%IGWP+I) global_here%RDES8(I)
              ENDDO
            global_here%igpp=global_here%igpp+4
            global_here%IGWP=global_here%IGWP+4
            DO I=1,3
              write(73,rec=global_here%igpp+i) global_here%rid8(i)
              WRITE(74,REC=global_here%IGWP+I) global_here%RID8(I)
              ENDDO
            global_here%igpp=global_here%igpp+3
            global_here%IGWP=global_here%IGWP+3
            DO I=1,3
              write(73,rec=global_here%igpp+i) global_here%aid8(i)
              WRITE(74,REC=global_here%IGWP+I) global_here%AID8(I)
              ENDDO
            global_here%igpp=global_here%igpp+3
            global_here%IGWP=global_here%IGWP+3
            ENDIF
          write(73,rec=global_here%igpp+1) global_here%ndsetsw
          write(73,rec=global_here%igpp+2) global_here%np
          write(73,rec=global_here%igpp+3) global_here%dt*global_here%nspoolgw
          write(73,rec=global_here%igpp+4) global_here%nspoolgw
          write(73,rec=global_here%igpp+5) 2
          global_here%igpp=global_here%igpp+5
          close(73)                    ! DO THIS TO FLUSH THE WRITE BUFFER
          open(73,file=s%DIRNAME//'/'//'fort.73',&
         access='direct',recl=nbyte)
          WRITE(74,REC=global_here%IGWP+1) global_here%NDSETSW
          WRITE(74,REC=global_here%IGWP+2) global_here%NP
          WRITE(74,REC=global_here%IGWP+3) global_here%DT*global_here%NSPOOLGW
          WRITE(74,REC=global_here%IGWP+4) global_here%NSPOOLGW
          WRITE(74,REC=global_here%IGWP+5) 2
          global_here%IGWP=global_here%IGWP+5
          CLOSE(74)                    ! DO THIS TO FLUSH THE WRITE BUFFER
          OPEN(74,FILE=S%DIRNAME//'/'//'fort.74',&
         ACCESS='DIRECT',RECL=NBYTE)
          ENDIF
#ifdef SWAN
!asey 101118: Added the output of radiation stress gradients.
       IGRadS=0
       IF(ABS(global_here%NOUTGW).EQ.1) THEN
          OPEN(164,FILE=S%DIRNAME//'/'//'rads.64')
          WRITE(164,3220) global_here%RUNDES,global_here%RUNID,global_here%AGRID
          WRITE(164,3645) global_here%NDSETSW,global_here%NP,global_here%DTDP*global_here%NSPOOLGW,global_here%NSPOOLGW,2
          IGRadS=2
          ENDIF
       IF(ABS(global_here%NOUTGW).EQ.2) THEN
          OPEN(164,FILE=TRIM(LOCALDIR)//'/'//'rads.64',&
           ACCESS='DIRECT',RECL=NByte)
          IF(NBYTE.EQ.4) THEN
             DO I=1,8
                WRITE(164,REC=IGRadS+I) global_here%RDES4(I)
                ENDDO
             IGRadS=IGRadS+8
             DO I=1,6
                WRITE(164,REC=IGRadS+I) global_here%RID4(I)
                ENDDO
             IGRadS=IGRadS+6
             DO I=1,6
                WRITE(164,REC=IGRadS+I) global_here%AID4(I)
                ENDDO
             IGRadS=IGRadS+6
             ENDIF
          IF(NBYTE.EQ.8) THEN
             DO I=1,4
                WRITE(164,REC=IGRadS+I) global_here%RDES8(I)
                ENDDO
             IGRadS=IGRadS+4
             DO I=1,3
                WRITE(164,REC=IGRadS+I) global_here%RID8(I)
                ENDDO
             IGRadS=IGRadS+3
             DO I=1,3
                WRITE(164,REC=IGRadS+I) global_here%AID8(I)
                ENDDO
             IGRadS=IGRadS+3
             ENDIF
          WRITE(164,REC=IGRadS+1) global_here%NDSETSW
          WRITE(164,REC=IGRadS+2) global_here%NP
          WRITE(164,REC=IGRadS+3) global_here%DT*global_here%NSPOOLGW
          WRITE(164,REC=IGRadS+4) global_here%NSPOOLGW
          WRITE(164,REC=IGRadS+5) 2
          IGRadS=IGRadS+5
          CLOSE(164)
          ENDIF
#endif

!...
!....INITILIZE GLOBAL CONCENTRATION SPOOL COUNTER
!....OPEN GLOBAL CONCENTRATION OUTPUT FILE
!....WRITE OUT HEADER INFORMATION INCLUDING global_here%NDSETSC
!....(NO. OF GLOBAL CONCENTRATION DATA SETS TO BE SPOOLED),
!....global_here%NP, global_here%DT*global_here%NSPOOLGC, global_here%NSPOOLGC, IRTYPE
!...
        global_here%NSCOUGC=0
        global_here%IGCP=0

        IF(ABS(global_here%NOUTGC).EQ.1) THEN
          OPEN(83,FILE=S%DIRNAME//'/'//'fort.83')
          WRITE(83,3220) global_here%RUNDES,global_here%RUNID,global_here%AGRID
          WRITE(83,3645) global_here%NDSETSC,global_here%NP,global_here%DTDP*global_here%NSPOOLGC,global_here%NSPOOLGC,1
          global_here%IGCP=2
          ENDIF

        IF(ABS(global_here%NOUTGC).EQ.2) THEN
          OPEN(83,FILE=S%DIRNAME//'/'//'fort.83',&
          ACCESS='DIRECT',RECL=NBYTE)
          IF(NBYTE.EQ.4) THEN
            DO I=1,8
              WRITE(83,REC=global_here%IGCP+I) global_here%RDES4(I)
              ENDDO
            global_here%IGCP=global_here%IGCP+8
            DO I=1,6
              WRITE(83,REC=global_here%IGCP+I) global_here%RID4(I)
              ENDDO
            global_here%IGCP=global_here%IGCP+6
            DO I=1,6
              WRITE(83,REC=global_here%IGCP+I) global_here%AID4(I)
              ENDDO
            global_here%IGCP=global_here%IGCP+6
            ENDIF
          IF(NBYTE.EQ.8) THEN
            DO I=1,4
              WRITE(83,REC=global_here%IGCP+I) global_here%RDES8(I)
              ENDDO
            global_here%IGCP=global_here%IGCP+4
            DO I=1,3
              WRITE(83,REC=global_here%IGCP+I) global_here%RID8(I)
              ENDDO
            global_here%IGCP=global_here%IGCP+3
            DO I=1,3
              WRITE(83,REC=global_here%IGCP+I) global_here%AID8(I)
              ENDDO
            global_here%IGCP=global_here%IGCP+3
            ENDIF
          WRITE(83,REC=global_here%IGCP+1) global_here%NDSETSC
          WRITE(83,REC=global_here%IGCP+2) global_here%NP
          WRITE(83,REC=global_here%IGCP+3) global_here%DT*global_here%NSPOOLGC
          WRITE(83,REC=global_here%IGCP+4) global_here%NSPOOLGC
          WRITE(83,REC=global_here%IGCP+5) 1
          global_here%IGCP=global_here%IGCP+5
          CLOSE(83)                    ! DO THIS TO FLUSH THE WRITE BUFFER
          OPEN(83,FILE=S%DIRNAME//'/'//'fort.83',&
         ACCESS='DIRECT',RECL=NBYTE)
          ENDIF

#ifdef harm
!...
!....INITIALIZE HARMONIC ANALYSIS MATRICES, MEAN AND SQUARE VECTORS
!...
        IF (global_here%IHARIND.EQ.1) THEN
           global_here%ICHA=0
           CALL HACOLDS(HAFREQ)
           IF(global_here%NHASE.EQ.1) CALL HACOLDSES(global_here%NSTAE)
           IF(global_here%NHASV.EQ.1) CALL HACOLDSVS(global_here%NSTAV)
           IF(global_here%NHAGE.EQ.1) CALL HACOLDSEG(global_here%NP)
           IF(global_here%NHAGV.EQ.1) CALL HACOLDSVG(global_here%NP)
           IF ( CHARMV) THEN
             DO I=1,global_here%NP
                global_here%ELAV(I)=0.D0
                global_here%XVELAV(I)=0.D0
                global_here%YVELAV(I)=0.D0
                global_here%ELVA(I)=0.D0
                global_here%XVELVA(I)=0.D0
                global_here%YVELVA(I)=0.D0
             ENDDO
           ENDIF !  charmv
        ENDIF
#endif
!
 1112 FORMAT(/,1X,79('_'))
 9883 FORMAT(' !!! NODE ',I6,' DRIED (LANDLOCKING)')

!.....Added sediment transport output files (Ethan Kubatko, 8-1-2003)

!.....INITILIZE GLOBAL BATHYMETRY SPOOL COUNTER
!.....OPEN GLOBAL BATHYMETRY OUTPUT FILE
!.....WRITE OUT HEADER INFORMATION INCLUDING global_here%NDSETSE
!.....(NO. OF GLOBAL ELEVATION DATA SETS TO BE SPOOLED),
!.....global_here%NP, global_here%DT*global_here%NSPOOLGE, global_here%NSPOOLGE, IRTYPE

        IF (global_here%SEDFLAG.GE.1) THEN
        
          global_here%NSCOUGE = 0
          global_here%IGEP = 0

          IF (ABS(global_here%NOUTGE).EQ.1) THEN
            OPEN(84,FILE=S%DIRNAME//'/'//'fort.84')
            OPEN(85,FILE=S%DIRNAME//'/'//'fort.85')
            WRITE(84,3220) global_here%RUNDES,global_here%RUNID,global_here%AGRID
            WRITE(84,3645) global_here%NDSETSE,global_here%NP,global_here%DTDP*global_here%NSPOOLGE,global_here%NSPOOLGE,1
            WRITE(85,3220) global_here%RUNDES,global_here%RUNID,global_here%AGRID
            WRITE(85,3645) global_here%NDSETSE,global_here%NP,global_here%DTDP*global_here%NSPOOLGE,global_here%NSPOOLGE,1
            global_here%IGEP=2
          ENDIF

          IF (ABS(global_here%NOUTGE).EQ.2) THEN
            OPEN(84,FILE=S%DIRNAME//'/'//'fort.84',&
          ACCESS='DIRECT',RECL=NBYTE)
            OPEN(85,FILE=S%DIRNAME//'/'//'fort.85',&
          ACCESS='DIRECT',RECL=NBYTE)
            IF(NBYTE.EQ.4) THEN
              DO I=1,8
                WRITE(84,REC=global_here%IGEP+I) global_here%RDES4(I)
                WRITE(85,REC=global_here%IGEP+I) global_here%RDES4(I)
              ENDDO
              global_here%IGEP=global_here%IGEP+8
              DO I=1,6
                WRITE(84,REC=global_here%IGEP+I) global_here%RID4(I)
                WRITE(85,REC=global_here%IGEP+I) global_here%RID4(I)
              ENDDO
              global_here%IGEP=global_here%IGEP+6
              DO I=1,6
                WRITE(84,REC=global_here%IGEP+I) global_here%AID4(I)
                WRITE(85,REC=global_here%IGEP+I) global_here%AID4(I)
              ENDDO
              global_here%IGEP=global_here%IGEP+6
              ENDIF
            IF(NBYTE.EQ.8) THEN
              DO I=1,4
                WRITE(84,REC=global_here%IGEP+I) global_here%RDES8(I)
                WRITE(85,REC=global_here%IGEP+I) global_here%RDES8(I)
              ENDDO
              global_here%IGEP=global_here%IGEP+4
              DO I=1,3
                WRITE(84,REC=global_here%IGEP+I) global_here%RID8(I)
                WRITE(85,REC=global_here%IGEP+I) global_here%RID8(I)
              ENDDO
              global_here%IGEP=global_here%IGEP+3
              DO I=1,3
                WRITE(84,REC=global_here%IGEP+I) global_here%AID8(I)
                WRITE(85,REC=global_here%IGEP+I) global_here%AID8(I)
              ENDDO
              global_here%IGEP=global_here%IGEP+3
            ENDIF
            WRITE(84,REC=global_here%IGEP+1) global_here%NDSETSE
            WRITE(84,REC=global_here%IGEP+2) global_here%NP
            WRITE(84,REC=global_here%IGEP+3) global_here%DT*global_here%NSPOOLGE
            WRITE(84,REC=global_here%IGEP+4) global_here%NSPOOLGE
            WRITE(84,REC=global_here%IGEP+5) 1
            WRITE(85,REC=global_here%IGEP+1) global_here%NDSETSE
            WRITE(85,REC=global_here%IGEP+2) global_here%NP
            WRITE(85,REC=global_here%IGEP+3) global_here%DT*global_here%NSPOOLGE
            WRITE(85,REC=global_here%IGEP+4) global_here%NSPOOLGE
            WRITE(85,REC=global_here%IGEP+5) 1
            global_here%IGEP=global_here%IGEP+5
            CLOSE(84)                    ! DO THIS TO FLUSH THE WRITE BUFFER
            CLOSE(85)                    ! DO THIS TO FLUSH THE WRITE BUFFER
            OPEN(84,FILE=S%DIRNAME//'/'//'fort.84',&
         ACCESS='DIRECT',RECL=NBYTE)
            OPEN(85,FILE=S%DIRNAME//'/'//'fort.85',&
         ACCESS='DIRECT',RECL=NBYTE)
          ENDIF
!...
!.....INITILIZE BATHYMETRY STATION SPOOL COUNTER
!.....OPEN ELEVATION STATION OUTPUT FILE
!.....WRITE OUT HEADER INFORMATION INCLUDING global_here%NTRSPE (NO. OF DATA PTS. AT EACH
!.....ELEVATION STATION), global_here%NSTAE, global_here%DT*global_here%NSPOOLE, global_here%NSPOOLE, IRTYPE
!...
          global_here%NSCOUE=0
          global_here%IESTP=0

          IF (ABS(global_here%NOUTE).EQ.1) THEN
            OPEN(82,FILE=S%DIRNAME//'/'//'fort.82')
            WRITE(82,3220) global_here%RUNDES,global_here%RUNID,global_here%AGRID
            WRITE(82,3645) global_here%NTRSPE,global_here%NSTAE,global_here%DTDP*global_here%NSPOOLE,global_here%NSPOOLE,1
            global_here%IESTP=2
          ENDIF

          IF (ABS(global_here%NOUTE).EQ.2) THEN
            OPEN(82,FILE=S%DIRNAME//'/'//'fort.82',&
          ACCESS='DIRECT',RECL=NBYTE)
            IF(NBYTE.EQ.4) THEN
              DO I=1,8
                WRITE(82,REC=global_here%IESTP+I) global_here%RDES4(I)
              ENDDO
              global_here%IESTP=global_here%IESTP+8
              DO I=1,6
                WRITE(82,REC=global_here%IESTP+I) global_here%RID4(I)
              ENDDO
              global_here%IESTP=global_here%IESTP+6
              DO I=1,6
                WRITE(82,REC=global_here%IESTP+I) global_here%AID4(I)
              ENDDO
              global_here%IESTP=global_here%IESTP+6
            ENDIF
            IF(NBYTE.EQ.8) THEN
              DO I=1,4
                WRITE(82,REC=global_here%IESTP+I) global_here%RDES8(I)
              ENDDO
              global_here%IESTP=global_here%IESTP+4
              DO I=1,3
                WRITE(82,REC=global_here%IESTP+I) global_here%RID8(I)
              ENDDO
              global_here%IESTP=global_here%IESTP+3
              DO I=1,3
                WRITE(82,REC=global_here%IESTP+I) global_here%AID8(I)
              ENDDO
              global_here%IESTP=global_here%IESTP+3
            ENDIF
            WRITE(82,REC=global_here%IESTP+1) global_here%NTRSPE
            WRITE(82,REC=global_here%IESTP+2) global_here%NSTAE
            WRITE(82,REC=global_here%IESTP+3) global_here%DT*global_here%NSPOOLE
            WRITE(82,REC=global_here%IESTP+4) global_here%NSPOOLE
            WRITE(82,REC=global_here%IESTP+5) 1
            global_here%IESTP=global_here%IESTP+5
            CLOSE(82)                    ! DO THIS TO FLUSH THE WRITE BUFFER
            OPEN(82,FILE=S%DIRNAME//'/'//'fort.82',&
         ACCESS='DIRECT',RECL=NBYTE)
          ENDIF

!.....INITILIZE GLOBAL BED LOAD SEDIMENT FLUX SPOOL COUNTER
!.....OPEN GLOBAL SEDIMENT FLUX OUTPUT FILE
!.....WRITE OUT HEADER INFORMATION INCLUDING global_here%NDSETSV
!.....(NO. OF GLOBAL VELOCITY DATA SETS TO BE SPOOLED),
!.....global_here%NP, global_here%DT*global_here%NSPOOLGV, global_here%NSPOOLGV, IRTYPE
!...
          global_here%NSCOUGV=0
          global_here%IGVP=0

          IF(ABS(global_here%NOUTGV).EQ.1) THEN
            OPEN(94,FILE=S%DIRNAME//'/'//'fort.94')
            WRITE(94,3220) global_here%RUNDES,global_here%RUNID,global_here%AGRID
            WRITE(94,3645) global_here%NDSETSV,global_here%NP,global_here%DTDP*global_here%NSPOOLGV,global_here%NSPOOLGV,2
            global_here%IGVP=2
          ENDIF

          IF(ABS(global_here%NOUTGV).EQ.2) THEN
            OPEN(94,FILE=S%DIRNAME//'/'//'fort.94',&
          ACCESS='DIRECT',RECL=NBYTE)
            IF(NBYTE.EQ.4) THEN
              DO I=1,8
                WRITE(94,REC=global_here%IGVP+I) global_here%RDES4(I)
              ENDDO
              global_here%IGVP=global_here%IGVP+8
              DO I=1,6
                WRITE(94,REC=global_here%IGVP+I) global_here%RID4(I)
              ENDDO
              global_here%IGVP=global_here%IGVP+6
              DO I=1,6
                WRITE(94,REC=global_here%IGVP+I) global_here%AID4(I)
              ENDDO
              global_here%IGVP=global_here%IGVP+6
            ENDIF
            IF(NBYTE.EQ.8) THEN
              DO I=1,4
                WRITE(94,REC=global_here%IGVP+I) global_here%RDES8(I)
              ENDDO
              global_here%IGVP=global_here%IGVP+4
              DO I=1,3
                WRITE(94,REC=global_here%IGVP+I) global_here%RID8(I)
              ENDDO
              global_here%IGVP=global_here%IGVP+3
              DO I=1,3
                WRITE(94,REC=global_here%IGVP+I) global_here%AID8(I)
              ENDDO
              global_here%IGVP=global_here%IGVP+3
            ENDIF
            WRITE(94,REC=global_here%IGVP+1) global_here%NDSETSV
            WRITE(94,REC=global_here%IGVP+2) global_here%NP
            WRITE(94,REC=global_here%IGVP+3) global_here%DT*global_here%NSPOOLGV
            WRITE(94,REC=global_here%IGVP+4) global_here%NSPOOLGV
            WRITE(94,REC=global_here%IGVP+5) 2
            global_here%IGVP=global_here%IGVP+5
            CLOSE(94)                    ! DO THIS TO FLUSH THE WRITE BUFFER
            OPEN(94,FILE=S%DIRNAME//'/'//'fort.94',&
         ACCESS='DIRECT',RECL=NBYTE)
          ENDIF
          
!.....INITILIZE GLOBAL SUSPENDED LOAD SEDIMENT FLUX SPOOL COUNTER
!.....OPEN GLOBAL SEDIMENT FLUX OUTPUT FILE
!.....WRITE OUT HEADER INFORMATION INCLUDING global_here%NDSETSV
!.....(NO. OF GLOBAL VELOCITY DATA SETS TO BE SPOOLED),
!.....global_here%NP, global_here%DT*global_here%NSPOOLGV, global_here%NSPOOLGV, IRTYPE
!...
          global_here%NSCOUGV=0
          global_here%IGVP=0

          IF(ABS(global_here%NOUTGV).EQ.1) THEN
            OPEN(96,FILE=S%DIRNAME//'/'//'fort.96')
            WRITE(96,3220) global_here%RUNDES,global_here%RUNID,global_here%AGRID
            WRITE(96,3645) global_here%NDSETSV,global_here%NP,global_here%DTDP*global_here%NSPOOLGV,global_here%NSPOOLGV,2
            global_here%IGVP=2
          ENDIF

          IF(ABS(global_here%NOUTGV).EQ.2) THEN
            OPEN(96,FILE=S%DIRNAME//'/'//'fort.96',&
          ACCESS='DIRECT',RECL=NBYTE)
            IF(NBYTE.EQ.4) THEN
              DO I=1,8
                WRITE(96,REC=global_here%IGVP+I) global_here%RDES4(I)
              ENDDO
              global_here%IGVP=global_here%IGVP+8
              DO I=1,6
                WRITE(96,REC=global_here%IGVP+I) global_here%RID4(I)
              ENDDO
              global_here%IGVP=global_here%IGVP+6
              DO I=1,6
                WRITE(96,REC=global_here%IGVP+I) global_here%AID4(I)
              ENDDO
              global_here%IGVP=global_here%IGVP+6
            ENDIF
            IF(NBYTE.EQ.8) THEN
              DO I=1,4
                WRITE(96,REC=global_here%IGVP+I) global_here%RDES8(I)
              ENDDO
              global_here%IGVP=global_here%IGVP+4
              DO I=1,3
                WRITE(96,REC=global_here%IGVP+I) global_here%RID8(I)
              ENDDO
              global_here%IGVP=global_here%IGVP+3
              DO I=1,3
                WRITE(96,REC=global_here%IGVP+I) global_here%AID8(I)
              ENDDO
              global_here%IGVP=global_here%IGVP+3
            ENDIF
            WRITE(96,REC=global_here%IGVP+1) global_here%NDSETSV
            WRITE(96,REC=global_here%IGVP+2) global_here%NP
            WRITE(96,REC=global_here%IGVP+3) global_here%DT*global_here%NSPOOLGV
            WRITE(96,REC=global_here%IGVP+4) global_here%NSPOOLGV
            WRITE(96,REC=global_here%IGVP+5) 2
            global_here%IGVP=global_here%IGVP+5
            CLOSE(96)                    ! DO THIS TO FLUSH THE WRITE BUFFER
            OPEN(96,FILE=S%DIRNAME//'/'//'fort.96',&
         ACCESS='DIRECT',RECL=NBYTE)
          ENDIF

        ENDIF
          
      RETURN
      END
