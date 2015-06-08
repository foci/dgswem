C**************************************************************************
C  mod history
C  v41.06mxxx - date - programmer - describe change 
C                    - mark change in code with  cinitials-mxxx
C
C  v41.11 - 09/14/01 - rl - from 41.09 - modified for NWS=-2
C  v41.09 - 06/30/01 - jw - from 41.08 - made minor mods as per vp version 41.05 
C  v41.08 - 06/22/01 - rl - from 41.07 - added 41.05m009 changes to HABSMIN
C                                        and ETA2
C  v41.07 - 04/09/01 - rl - from 41.06 - initialized PRN1(), PRN2() for NRS<>0
C**************************************************************************
C
        SUBROUTINE COLDSTART()
C
C**************************************************************************
C
C  COLD START PROGRAM SETUP ROUTINE 
C
C**************************************************************************
C
      USE SIZES
      USE GLOBAL
      USE HARM
      USE WIND
      USE OWIWIND, ONLY : NWS12INIT, NWS12GET
      USE NodalAttributes, ONLY : STARTDRY, GeoidOffset, LoadGeoidOffset
      IMPLICIT NONE

      INTEGER i,j
C
      ITHS = 0

C...
C....SET AT REST INITIAL CONDITION OVER WHOLE DOMAIN
C....IF BOTTOM IS ABOVE THE GEIOD -> DRY NODE
C....IF BOTTOM IS INITIALLY BELOW THE GEIOD AND STARTDRY=-88888 -> DRY NODE
C...
        HABSMIN=0.8d0*H0
        DO I=1,NP
          UU1(I) =0.D0
          VV1(I) =0.D0
          UU2(I) =0.D0
          VV2(I) =0.D0
          ETA2(I)=0.D0
          NODEREP(I)=MAX0(NODEWETMIN,NODEDRYMIN)
          NNODECODE(I)=1
          IF(NOLIFA.EQ.2) THEN
            HTOT=DP(I)+ETA2(I)
            IF(HTOT.LE.H0) THEN
              NNODECODE(I)=0
              ETA2(I)=H0-DP(I)
C              ELSE
C              IF(STARTDRY(I).EQ.-88888) THEN
C                NNODECODE(I)=0
C                ETA2(I)=H0-DP(I)
C                ENDIF
            ENDIF
          ENDIF
          ETA1(I)=ETA2(I)
          ETAS(I)=0.D0
          CH1(I)=0.d0
          END DO

C...
C....INITIALIZE THE ELEVATION SPECIFIED BOUNDARY CONDITION IF IT REQUIRES THE USE
C....OF THE UNIT 19 FILE.
C...

        IF((NOPE.GT.0).AND.(NBFR.EQ.0)) THEN
          IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,1112)
          WRITE(16,1112)
          IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,1977)
          WRITE(16,1977)
 1977     FORMAT(/,1X,'ELEVATION SPECIFIED INFORMATION READ FROM UNIT ',
     &           '19',/)
          OPEN(19,FILE=DIRNAME//'/'//'fort.19')
          READ(19,*) ETIMINC
          DO J=1,NETA
             READ(19,*) ESBIN1(J)
          END DO
          DO J=1,NETA
             READ(19,*) ESBIN2(J)
          END DO
          ETIME1 = STATIM*86400.D0
          ETIME2 = ETIME1 + ETIMINC
        ENDIF

C....INITIALIZE THE NORMAL FLOW BOUNDARY CONDITION

        DO I=1,NVEL
          QN2(I)=0.D0
          QN1(I)=0.D0
          QN0(I)=0.D0
        END DO

        IF((NFLUXF.EQ.1).AND.(NFFR.EQ.0)) THEN
          IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,1112)
          WRITE(16,1112)
          IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,1979)
          WRITE(16,1979)
 1979     FORMAT(/,1X,'NORMAL FLOW INFORMATION READ FROM UNIT 20',/)
          OPEN(20,FILE=DIRNAME//'/'//'fort.20')
          READ(20,*) FTIMINC
          DO J=1,NVEL
            QNIN1(J)=0.D0
            IF((LBCODEI(J).EQ.2).OR.(LBCODEI(J).EQ.12)
     &                          .OR.(LBCODEI(J).EQ.22))
     &                                          READ(20,*) QNIN1(J)
            END DO
          DO J=1,NVEL
            QNIN2(J)=0.D0
            IF((LBCODEI(J).EQ.2).OR.(LBCODEI(J).EQ.12)
     &                          .OR.(LBCODEI(J).EQ.22))
     &                                          READ(20,*) QNIN2(J)
            END DO
          QTIME1 = STATIM*86400.D0
          QTIME2 = QTIME1 + FTIMINC
          ENDIF

C...INPUT METEOROLOGICAL INFORMATION FROM UNIT 22 OR UNIT 200 SERIES
C....IF FLEET NUMERIC WIND DATA IS USED, FIND BEGINNING TIME IN FILE,
C....NOTE: CAN'T DEAL WITH WIND THAT STARTS AFTER WREFTIM!!!!!!!!!!!!
C....READ IN AND INTERPOLATE IN SPACE ONTO THE ADCIRC GRID THE
C....TIME LEVEL 1 AND LEVEL 2 WIND FIELDS

        IF(NWS.NE.0) THEN
          DO I=1,NP
            WSX1(I)=0.D0
            WSY1(I)=0.D0
            PR1(I) =0.D0
            WSX2(I)=0.D0
            WSY2(I)=0.D0
            PR2(I) =0.D0
            ENDDO

          IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,1112)
          WRITE(16,1112)
          IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,1980)
          WRITE(16,1980)
 1980     FORMAT(/,1X,'WIND (AND PRESSURE) INFORMATION READ.',/)
          ENDIF

        IF(NWS.EQ.1) THEN
          OPEN(22,FILE=DIRNAME//'/'//'fort.22')
          ENDIF

        IF(ABS(NWS).EQ.2) THEN
          OPEN(22,FILE=DIRNAME//'/'//'fort.22')
          READ(22,*) (NHG,WVNX1(I),WVNY1(I),PRN1(I),I=1,NP)
          READ(22,*) (NHG,WVNX2(I),WVNY2(I),PRN2(I),I=1,NP)
          WTIME1 = STATIM*86400.D0
          WTIME2 = WTIME1 + WTIMINC
          ENDIF

        IF(NWS.EQ.3) THEN
          OPEN(22,FILE=DIRNAME//'/'//'fort.22')
 2222     CALL NWS3GET( X, Y, SLAM, SFEA, WVNX2, WVNY2, IWTIME, IWYR,
     &                  WTIMED, NP, NWLON, NWLAT, WLATMAX, WLONMIN,
     &                  WLATINC, WLONINC, ICS, NSCREEN, ScreenUnit )
          IF(IWYR.NE.IREFYR) THEN
            IWTIMEP=IWTIME
            DO I=1,NP
              WVNX1(I)=WVNX2(I)
              WVNY1(I)=WVNY2(I)
              END DO
            GOTO 2222
            ENDIF
          IF(WTIMED.LE.WREFTIM) THEN
            IWTIMEP=IWTIME
            DO I=1,NP
              WVNX1(I)=WVNX2(I)
              WVNY1(I)=WVNY2(I)
              END DO
            GOTO 2222
            ENDIF
          IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) 
     &         WRITE(6,*)'FOUND WIND DATA AT TIME= ',IWTIMEP
          WRITE(16,*) 'FOUND WIND DATA AT TIME= ',IWTIMEP
          IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) 
     &         WRITE(6,*)'FOUND WIND DATA AT TIME= ',IWTIME
          WRITE(16,*) 'FOUND WIND DATA AT TIME= ',IWTIME
          WTIME2=WTIMED-WREFTIM                  !CAST INTO MODEL TIME REFRENCE
          WTIME1=WTIME2-WTIMINC
          ENDIF

        IF(ABS(NWS).EQ.4) THEN
          OPEN(22,FILE=DIRNAME//'/'//'fort.22')
          WTIME1 = STATIM*86400.D0
          WTIME2=WTIME1+WTIMINC
          CALL NWS4GET(WVNX1,WVNY1,PRN1,NP,RHOWAT0,G)
          CALL NWS4GET(WVNX2,WVNY2,PRN2,NP,RHOWAT0,G)
          ENDIF

        IF(ABS(NWS).EQ.5) THEN
          OPEN(22,FILE=DIRNAME//'/'//'fort.22')
          READ(22,*) (NHG,WVNX1(I),WVNY1(I),PRN1(I),I=1,NP)
          READ(22,*) (NHG,WVNX2(I),WVNY2(I),PRN2(I),I=1,NP)
          WTIME1 = STATIM*86400.D0
          WTIME2 = WTIME1 + WTIMINC
          ENDIF

        IF(NWS.EQ.6) THEN
          OPEN(22,FILE=DIRNAME//'/'//'fort.22')
          CALL NWS6GET(X,Y,SLAM,SFEA,WVNX1,WVNY1,PRN1,NP,NWLON,NWLAT,
     &                 WLATMAX,WLONMIN,WLATINC,WLONINC,ICS,RHOWAT0,G)
          CALL NWS6GET(X,Y,SLAM,SFEA,WVNX2,WVNY2,PRN2,NP,NWLON,NWLAT,
     &                 WLATMAX,WLONMIN,WLATINC,WLONINC,ICS,RHOWAT0,G)
          WTIME1 = STATIM*86400.D0
          WTIME2 = WTIME1 + WTIMINC
          ENDIF

        IF(NWS.EQ.10) THEN
          WTIME1=STATIM*86400.D0
          WTIME2=WTIME1+WTIMINC
          NWSGGWI=-1
          CALL NWS10GET(NWSGGWI,SLAM,SFEA,WVNX1,WVNY1,PRN1,NP,RHOWAT0,G,
     &                  NWLON,NWLAT,WTIMINC) !JUST COMPUTE INTERPOLATING FACTORS
          NWSGGWI=1
          CALL NWS10GET(NWSGGWI,SLAM,SFEA,WVNX2,WVNY2,PRN2,NP,RHOWAT0,G,
     &                  NWLON,NWLAT,WTIMINC) !NOW INTERPOLATE 1st WIND FIELD
          ENDIF

        IF(NWS.EQ.11) THEN
          WTIME1=STATIM*86400.D0
          WTIME2=WTIME1+WTIMINC
          NWSEGWI=0
          IDSETFLG=0
          IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,1197)
          WRITE(16,1197)
 1197     FORMAT(/,1X,'THE E29 MET GRID INTERPOLATING FACTORS ARE ',
     &                'BEING COMPUTED ')
          CALL NWS11GET(NWSEGWI,IDSETFLG,SLAM,SFEA,WVNX1,WVNY1,PRN1,NP,
     &                  RHOWAT0,G)  !JUST COMPUTE INTERPOLATING FACTORS
          IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,1198)
          WRITE(16,1198)
 1198     FORMAT(1X,'FINISHED COMPUTING E29 INTERPOLATING FACTORS',/)
          NWSEGWI=1
          IDSETFLG=1
          CALL NWS11GET(NWSEGWI,IDSETFLG,SLAM,SFEA,WVNX2,WVNY2,PRN2,NP,
     &                  RHOWAT0,G) !NOW INTERPOLATE 1st WIND FIELD
        ENDIF
        
Cek added nws = 12 for owi winds
        
        IF(ABS(NWS).EQ.12) THEN
          CALL NWS12INIT( WVNX1, WVNY1, PRN1, NP, RHOWAT0, G )
          CALL NWS12GET(  WVNX1, WVNY1, PRN1, NP, RHOWAT0, G )
          CALL NWS12GET(  WVNX2, WVNY2, PRN2, NP, RHOWAT0, G )
          WTIME1 = STATIM*86400.D0
          WTIME2 = WTIME1 + WTIMINC
        ENDIF

C...INPUT RADIATION STRESS INFORMATION FROM UNIT 23
C....READ IN THE TIME LEVEL 1 AND LEVEL 2 FIELDS

        IF(NRS.EQ.1) THEN
          IF(NWS.EQ.0) THEN
            DO I=1,NP
              WSX1(I)=0.D0
              WSY1(I)=0.D0
              WSX2(I)=0.D0
              WSY2(I)=0.D0
              PRN1(I)=0.D0    !need to be initialized
              PRN2(I)=0.D0    !even if not used
              ENDDO
            ENDIF
          OPEN(23,FILE=DIRNAME//'/'//'fort.23')
          RSTIME1 = STATIM*86400.D0
          RSTIME2 = RSTIME1+RSTIMINC
          IF (FRW.EQ.0) THEN
            CALL RSGET(RSNX1,RSNY1,NP)
            CALL RSGET(RSNX2,RSNY2,NP)
          ENDIF
          IF (FRW.EQ.1) THEN
C            CALL RSGET_MORE(RSNX1,RSNY1,WAVE_H1,WAVE_T1,WAVE_A1,
C     &                                                  WAVE_D1,NP)
C            CALL RSGET_MORE(RSNX2,RSNY2,WAVE_H2,WAVE_T2,WAVE_A2,
C     &                                                  WAVE_D2,NP)
            
          ENDIF
#ifdef SWAN   
Casey 101118: Added this call to initialize the radiation stress gradients.
         IF(NRS.EQ.3) THEN
           IF(.NOT.ALLOCATED(RSNX1)) ALLOCATE(RSNX1(1:NP))
           IF(.NOT.ALLOCATED(RSNX2)) ALLOCATE(RSNX2(1:NP))
           IF(.NOT.ALLOCATED(RSNY1)) ALLOCATE(RSNY1(1:NP))
           IF(.NOT.ALLOCATED(RSNY2)) ALLOCATE(RSNY2(1:NP))
           DO I=1,NP
             RSNX1(I) = 0.D0
             RSNX2(I) = 0.D0
             RSNY1(I) = 0.D0
             RSNY2(I) = 0.D0
           ENDDO
         ENDIF
#endif
          IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,1112)
          WRITE(16,1112)
          IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,1981)
          WRITE(16,1981)
 1981     FORMAT(/,1X,'RADIATION STRESS INFORMATION READ.',/)
          ENDIF


C...
C...LINES TO USE TIDAL POTENTIAL FORCING
C...
       if (CTIP) then
          DO I=1,NP
             TIP2(I)=0.0
          END DO
       endif

CWET...
CWET...THE FOLLOWING LINES ARE FOR WETTING AND DRYING
CWET...Dry any landlocked nodes by checking that they are connected to at
CWET...least 1 functioning element.
CWET...
        IF(NOLIFA.EQ.2) THEN
          DO I=1,NP
            MJU(I)=0
            ENDDO
          DO I=1,NE
            NM1=NM(I,1)
            NM2=NM(I,2)
            NM3=NM(I,3)
            NC1=NNODECODE(NM1)
            NC2=NNODECODE(NM2)
            NC3=NNODECODE(NM3)
            NCELE=NC1*NC2*NC3
            MJU(NM1)=MJU(NM1)+NCELE
            MJU(NM2)=MJU(NM2)+NCELE
            MJU(NM3)=MJU(NM3)+NCELE
          ENDDO
          DO I=1,NP
            IF((NNODECODE(I).EQ.1).AND.(MJU(I).EQ.0)) THEN
              NNODECODE(I)=0
              IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(*,9883) I
              WRITE(16,9883) I
              ENDIF
          ENDDO
        ENDIF
C...
C......INITIALIZE 3D SOLUTION
C...

C...LINES TO RUN THE CODE IN 3D VS MODE.

      if (C3DVS) then
c       CALL VSSTUP(DT,STATIM,NBYTE,RUNDES,RUNID,AGRID,NT)
      endif

C...LINES TO RUN THE CODE IN 3D DSS MODE

      if (C3DDSS) then
c       CALL DSSSTUP(DT,STATIM,NBYTE,RUNDES,RUNID,AGRID,NT)
      endif


C...
C....INITILIZE ELEVATION STATION SPOOL COUNTER
C....OPEN ELEVATION STATION OUTPUT FILE
C....WRITE OUT HEADER INFORMATION INCLUDING NTRSPE (NO. OF DATA PTS. AT EACH
C....ELEVATION STATION), NSTAE, DT*NSPOOLE, NSPOOLE, IRTYPE
C...
        NSCOUE=0
        IESTP=0

 3220   FORMAT(1X,A32,2X,A24,2X,A24)
 3645   FORMAT(1X,I10,1X,I10,1X,E15.7,1X,I5,1X,I5)

        IF(ABS(NOUTE).EQ.1) THEN
          OPEN(61,FILE=DIRNAME//'/'//'fort.61')
          WRITE(61,3220) RUNDES,RUNID,AGRID
          WRITE(61,3645) NTRSPE,NSTAE,DTDP*NSPOOLE,NSPOOLE,1
          IESTP=2
          ENDIF

        IF(ABS(NOUTE).EQ.2) THEN
          OPEN(61,FILE=DIRNAME//'/'//'fort.61',
     &          ACCESS='DIRECT',RECL=NBYTE)
          IF(NBYTE.EQ.4) THEN
            DO I=1,8
              WRITE(61,REC=IESTP+I) RDES4(I)
              ENDDO
            IESTP=IESTP+8
            DO I=1,6
              WRITE(61,REC=IESTP+I) RID4(I)
              ENDDO
            IESTP=IESTP+6
            DO I=1,6
              WRITE(61,REC=IESTP+I) AID4(I)
              ENDDO
            IESTP=IESTP+6
            ENDIF
          IF(NBYTE.EQ.8) THEN
            DO I=1,4
              WRITE(61,REC=IESTP+I) RDES8(I)
              ENDDO
            IESTP=IESTP+4
            DO I=1,3
              WRITE(61,REC=IESTP+I) RID8(I)
              ENDDO
            IESTP=IESTP+3
            DO I=1,3
              WRITE(61,REC=IESTP+I) AID8(I)
              ENDDO
            IESTP=IESTP+3
            ENDIF
          WRITE(61,REC=IESTP+1) NTRSPE
          WRITE(61,REC=IESTP+2) NSTAE
          WRITE(61,REC=IESTP+3) DT*NSPOOLE
          WRITE(61,REC=IESTP+4) NSPOOLE
          WRITE(61,REC=IESTP+5) 1
          IESTP=IESTP+5
          CLOSE(61)                    ! DO THIS TO FLUSH THE WRITE BUFFER
          OPEN(61,FILE=DIRNAME//'/'//'fort.61',
     &         ACCESS='DIRECT',RECL=NBYTE)
          ENDIF

C...
C....INITILIZE VELOCITY STATION SPOOL COUNTER
C....OPEN VELOCITY STATION OUTPUT FILE
C....WRITE OUT HEADER INFORMATION INCLUDING NTRSPV (NO. OF DATA PTS. AT EACH
C....VELOCITY STATION), NSTAV, DT*NSPOOLV, NSPOOLV, IRTYPE
C...
        NSCOUV=0
        IVSTP=0

        IF(ABS(NOUTV).EQ.1) THEN
          OPEN(62,FILE=DIRNAME//'/'//'fort.62')
          WRITE(62,3220) RUNDES,RUNID,AGRID
          WRITE(62,3645) NTRSPV,NSTAV,DTDP*NSPOOLV,NSPOOLV,2
          IVSTP=2
          ENDIF

        IF(ABS(NOUTV).EQ.2) THEN
          OPEN(62,FILE=DIRNAME//'/'//'fort.62',
     &          ACCESS='DIRECT',RECL=NBYTE)
          IF(NBYTE.EQ.4) THEN
            DO I=1,8
              WRITE(62,REC=IVSTP+I) RDES4(I)
              ENDDO
            IVSTP=IVSTP+8
            DO I=1,6
              WRITE(62,REC=IVSTP+I) RID4(I)
              ENDDO
            IVSTP=IVSTP+6
            DO I=1,6
              WRITE(62,REC=IVSTP+I) AID4(I)
              ENDDO
            IVSTP=IVSTP+6
            ENDIF
          IF(NBYTE.EQ.8) THEN
            DO I=1,4
              WRITE(62,REC=IVSTP+I) RDES8(I)
              ENDDO
            IVSTP=IVSTP+4
            DO I=1,3
              WRITE(62,REC=IVSTP+I) RID8(I)
              ENDDO
            IVSTP=IVSTP+3
            DO I=1,3
              WRITE(62,REC=IVSTP+I) AID8(I)
              ENDDO
            IVSTP=IVSTP+3
            ENDIF
          WRITE(62,REC=IVSTP+1) NTRSPV
          WRITE(62,REC=IVSTP+2) NSTAV
          WRITE(62,REC=IVSTP+3) DT*NSPOOLV
          WRITE(62,REC=IVSTP+4) NSPOOLV
          WRITE(62,REC=IVSTP+5) 2
          IVSTP=IVSTP+5
          CLOSE(62)                    ! DO THIS TO FLUSH THE WRITE BUFFER
          OPEN(62,FILE=DIRNAME//'/'//'fort.62',
     &         ACCESS='DIRECT',RECL=NBYTE)
          ENDIF

C...
C....INITILIZE CONCENTRATION STATION SPOOL COUNTER
C....OPEN ELEVATION STATION OUTPUT FILE
C....WRITE OUT HEADER INFORMATION INCLUDING NTRSPC (NO. OF DATA PTS. AT EACH
C....CONCENTRATION STATION), NSTAC, DT*NSPOOLC, NSPOOLC, IRTYPE
C...
        NSCOUC=0
        ICSTP=0

        IF(ABS(NOUTC).EQ.1) THEN
          OPEN(81,FILE=DIRNAME//'/'//'fort.81')
          WRITE(81,3220) RUNDES,RUNID,AGRID
          WRITE(81,3645) NTRSPC,NSTAC,DTDP*NSPOOLC,NSPOOLC,1
          ICSTP=2
          ENDIF

        IF(ABS(NOUTC).EQ.2) THEN
          OPEN(81,FILE=DIRNAME//'/'//'fort.81',
     &          ACCESS='DIRECT',RECL=NBYTE)
          IF(NBYTE.EQ.4) THEN
            DO I=1,8
              WRITE(81,REC=ICSTP+I) RDES4(I)
              ENDDO
            ICSTP=ICSTP+8
            DO I=1,6
              WRITE(81,REC=ICSTP+I) RID4(I)
              ENDDO
            ICSTP=ICSTP+6
            DO I=1,6
              WRITE(81,REC=ICSTP+I) AID4(I)
              ENDDO
            ICSTP=ICSTP+6
            ENDIF
          IF(NBYTE.EQ.8) THEN
            DO I=1,4
              WRITE(81,REC=ICSTP+I) RDES8(I)
              ENDDO
            ICSTP=ICSTP+4
            DO I=1,3
              WRITE(81,REC=ICSTP+I) RID8(I)
              ENDDO
            ICSTP=ICSTP+3
            DO I=1,3
              WRITE(81,REC=ICSTP+I) AID8(I)
              ENDDO
            ICSTP=ICSTP+3
            ENDIF
          WRITE(81,REC=ICSTP+1) NTRSPC
          WRITE(81,REC=ICSTP+2) NSTAC
          WRITE(81,REC=ICSTP+3) DT*NSPOOLC
          WRITE(81,REC=ICSTP+4) NSPOOLC
          WRITE(81,REC=ICSTP+5) 1
          ICSTP=ICSTP+5
          CLOSE(81)                    ! DO THIS TO FLUSH THE WRITE BUFFER
          OPEN(81,FILE=DIRNAME//'/'//'fort.81',
     &         ACCESS='DIRECT',RECL=NBYTE)
          ENDIF


C...
C....INITILIZE BATHYMETRY STATION SPOOL COUNTER
C....OPEN ELEVATION STATION OUTPUT FILE
C....WRITE OUT HEADER INFORMATION INCLUDING NTRSPE (NO. OF DATA PTS. AT EACH
C....ELEVATION STATION), NSTAE, DT*NSPOOLE, NSPOOLE, IRTYPE
C...
        NSCOUE=0
        IESTP=0

        IF(ABS(NOUTE).EQ.1) THEN
          OPEN(82,FILE=DIRNAME//'/'//'fort.82')
          WRITE(82,3220) RUNDES,RUNID,AGRID
          WRITE(82,3645) NTRSPE,NSTAE,DTDP*NSPOOLE,NSPOOLE,1
          IESTP=2
          ENDIF

        IF(ABS(NOUTE).EQ.2) THEN
          OPEN(82,FILE=DIRNAME//'/'//'fort.82',
     &          ACCESS='DIRECT',RECL=NBYTE)
          IF(NBYTE.EQ.4) THEN
            DO I=1,8
              WRITE(82,REC=IESTP+I) RDES4(I)
              ENDDO
            IESTP=IESTP+8
            DO I=1,6
              WRITE(82,REC=IESTP+I) RID4(I)
              ENDDO
            IESTP=IESTP+6
            DO I=1,6
              WRITE(82,REC=IESTP+I) AID4(I)
              ENDDO
            IESTP=IESTP+6
            ENDIF
          IF(NBYTE.EQ.8) THEN
            DO I=1,4
              WRITE(82,REC=IESTP+I) RDES8(I)
              ENDDO
            IESTP=IESTP+4
            DO I=1,3
              WRITE(82,REC=IESTP+I) RID8(I)
              ENDDO
            IESTP=IESTP+3
            DO I=1,3
              WRITE(82,REC=IESTP+I) AID8(I)
              ENDDO
            IESTP=IESTP+3
            ENDIF
          WRITE(82,REC=IESTP+1) NTRSPE
          WRITE(82,REC=IESTP+2) NSTAE
          WRITE(82,REC=IESTP+3) DT*NSPOOLE
          WRITE(82,REC=IESTP+4) NSPOOLE
          WRITE(82,REC=IESTP+5) 1
          IESTP=IESTP+5
          CLOSE(82)                    ! DO THIS TO FLUSH THE WRITE BUFFER
          OPEN(82,FILE=DIRNAME//'/'//'fort.82',
     &         ACCESS='DIRECT',RECL=NBYTE)
          ENDIF


C...
C....INITILIZE METEOROLOGICAL STATION SPOOL COUNTERS
C....OPEN METEOROLOGICAL STATION OUTPUT FILES
C....WRITE OUT HEADER INFORMATION INCLUDING NTRSPM (NO. OF DATA PTS. AT EACH
C....METEOROLOGICAL STATION), NSTAM, DT*NSPOOLM, NSPOOLM, IRTYPE
C...
        NSCOUM=0
        IPSTP=0
        IWSTP=0

        IF(ABS(NOUTM).EQ.1) THEN
          OPEN(71,FILE=DIRNAME//'/'//'fort.71')
          WRITE(71,3220) RUNDES,RUNID,AGRID
          WRITE(71,3645) NTRSPM,NSTAM,DTDP*NSPOOLM,NSPOOLM,1
          IPSTP=2
          OPEN(72,FILE=DIRNAME//'/'//'fort.72')
          WRITE(72,3220) RUNDES,RUNID,AGRID
          WRITE(72,3645) NTRSPM,NSTAM,DTDP*NSPOOLM,NSPOOLM,2
          IWSTP=2
          ENDIF

        IF(ABS(NOUTM).EQ.2) THEN
          OPEN(71,FILE=DIRNAME//'/'//'fort.71',
     &          ACCESS='DIRECT',RECL=NBYTE)
          OPEN(72,FILE=DIRNAME//'/'//'fort.72',
     &          ACCESS='DIRECT',RECL=NBYTE)
          IF(NBYTE.EQ.4) THEN
            DO I=1,8
              WRITE(71,REC=IPSTP+I) RDES4(I)
              WRITE(72,REC=IWSTP+I) RDES4(I)
              ENDDO
            IPSTP=IPSTP+8
            IWSTP=IWSTP+8
            DO I=1,6
              WRITE(71,REC=IPSTP+I) RID4(I)
              WRITE(72,REC=IWSTP+I) RID4(I)
              ENDDO
            IPSTP=IPSTP+6
            IWSTP=IWSTP+6
            DO I=1,6
              WRITE(71,REC=IPSTP+I) AID4(I)
              WRITE(72,REC=IWSTP+I) AID4(I)
              ENDDO
            IPSTP=IPSTP+6
            IWSTP=IWSTP+6
            ENDIF
          IF(NBYTE.EQ.8) THEN
            DO I=1,4
              WRITE(71,REC=IPSTP+I) RDES8(I)
              WRITE(72,REC=IWSTP+I) RDES8(I)
              ENDDO
            IPSTP=IPSTP+4
            IWSTP=IWSTP+4
            DO I=1,3
              WRITE(71,REC=IPSTP+I) RID8(I)
              WRITE(72,REC=IWSTP+I) RID8(I)
              ENDDO
            IPSTP=IPSTP+3
            IWSTP=IWSTP+3
            DO I=1,3
              WRITE(71,REC=IPSTP+I) AID8(I)
              WRITE(72,REC=IWSTP+I) AID8(I)
              ENDDO
            IPSTP=IPSTP+3
            IWSTP=IWSTP+3
            ENDIF
          WRITE(71,REC=IPSTP+1) NTRSPM
          WRITE(71,REC=IPSTP+2) NSTAM
          WRITE(71,REC=IPSTP+3) DT*NSPOOLM
          WRITE(71,REC=IPSTP+4) NSPOOLM
          WRITE(71,REC=IPSTP+5) 1
          WRITE(72,REC=IWSTP+1) NTRSPM
          WRITE(72,REC=IWSTP+2) NSTAM
          WRITE(72,REC=IWSTP+3) DT*NSPOOLM
          WRITE(72,REC=IWSTP+4) NSPOOLM
          WRITE(72,REC=IWSTP+5) 2
          IPSTP=IPSTP+5
          IWSTP=IWSTP+5
          CLOSE(71)                    ! DO THIS TO FLUSH THE WRITE BUFFER
          CLOSE(72)                    ! DO THIS TO FLUSH THE WRITE BUFFER
          OPEN(71,FILE=DIRNAME//'/'//'fort.71',
     &         ACCESS='DIRECT',RECL=NBYTE)
          OPEN(72,FILE=DIRNAME//'/'//'fort.72',
     &         ACCESS='DIRECT',RECL=NBYTE)
          ENDIF

C...
C....INITILIZE GLOBAL ELEVATION SPOOL COUNTER
C....OPEN GLOBAL ELEVATION OUTPUT FILE
C....WRITE OUT HEADER INFORMATION INCLUDING NDSETSE
C....(NO. OF GLOBAL ELEVATION DATA SETS TO BE SPOOLED),
C....NP, DT*NSPOOLGE, NSPOOLGE, IRTYPE
C...
        NSCOUGE=0
        IGEP=0

        IF(ABS(NOUTGE).EQ.1) THEN
          OPEN(63,FILE=DIRNAME//'/'//'fort.63')
          WRITE(63,3220) RUNDES,RUNID,AGRID
          WRITE(63,3645) NDSETSE,NE,DTDP*NSPOOLGE,NSPOOLGE,1
          IGEP=2
          ENDIF

        IF(ABS(NOUTGE).EQ.1) THEN
          OPEN(88,FILE=DIRNAME//'/'//'fort.88')
          WRITE(88,3220) RUNDES,RUNID,AGRID
          WRITE(88,3645) NDSETSE,NE,DTDP*NSPOOLGE,NSPOOLGE,1
          IGEP=2
          ENDIF

        IF(ABS(NOUTGE).EQ.1) THEN
          OPEN(89,FILE=DIRNAME//'/'//'fort.89')
          WRITE(89,3220) RUNDES,RUNID,AGRID
          WRITE(89,3645) NDSETSE,NE,DTDP*NSPOOLGE,NSPOOLGE,1
          IGEP=2
          ENDIF

        IF(ABS(NOUTGE).EQ.1) THEN
          OPEN(4441,FILE=DIRNAME//'/'//'fort.4l')
          WRITE(4441,3220) RUNDES,RUNID,AGRID
          WRITE(4441,3645) NDSETSE,NE,DTDP*NSPOOLGE,NSPOOLGE,1
          IGEP=2
          ENDIF

        IF(ABS(NOUTGE).EQ.2) THEN
          OPEN(63,FILE=DIRNAME//'/'//'fort.63',
     &          ACCESS='DIRECT',RECL=NBYTE)
          IF(NBYTE.EQ.4) THEN
            DO I=1,8
              WRITE(63,REC=IGEP+I) RDES4(I)
              ENDDO
            IGEP=IGEP+8
            DO I=1,6
              WRITE(63,REC=IGEP+I) RID4(I)
              ENDDO
            IGEP=IGEP+6
            DO I=1,6
              WRITE(63,REC=IGEP+I) AID4(I)
              ENDDO
            IGEP=IGEP+6
            ENDIF
          IF(NBYTE.EQ.8) THEN
            DO I=1,4
              WRITE(63,REC=IGEP+I) RDES8(I)
              ENDDO
            IGEP=IGEP+4
            DO I=1,3
              WRITE(63,REC=IGEP+I) RID8(I)
              ENDDO
            IGEP=IGEP+3
            DO I=1,3
              WRITE(63,REC=IGEP+I) AID8(I)
              ENDDO
            IGEP=IGEP+3
            ENDIF
          WRITE(63,REC=IGEP+1) NDSETSE
          WRITE(63,REC=IGEP+2) NE
          WRITE(63,REC=IGEP+3) DT*NSPOOLGE
          WRITE(63,REC=IGEP+4) NSPOOLGE
          WRITE(63,REC=IGEP+5) 1
          IGEP=IGEP+5
          CLOSE(63)                    ! DO THIS TO FLUSH THE WRITE BUFFER
          OPEN(63,FILE=DIRNAME//'/'//'fort.63',
     &         ACCESS='DIRECT',RECL=NBYTE)
          ENDIF

C...
C....INITILIZE GLOBAL VELOCITY SPOOL COUNTER
C....OPEN GLOBAL VELOCITY OUTPUT FILE
C....WRITE OUT HEADER INFORMATION INCLUDING NDSETSV
C....(NO. OF GLOBAL VELOCITY DATA SETS TO BE SPOOLED),
C....NP, DT*NSPOOLGV, NSPOOLGV, IRTYPE
C...
        NSCOUGV=0
        IGVP=0

        IF(ABS(NOUTGV).EQ.1) THEN
          OPEN(64,FILE=DIRNAME//'/'//'fort.64')
          WRITE(64,3220) RUNDES,RUNID,AGRID
          WRITE(64,3645) NDSETSV,NE,DTDP*NSPOOLGV,NSPOOLGV,2
          IGVP=2
          ENDIF

        IF(ABS(NOUTGV).EQ.2) THEN
          OPEN(64,FILE=DIRNAME//'/'//'fort.64',
     &          ACCESS='DIRECT',RECL=NBYTE)
          IF(NBYTE.EQ.4) THEN
            DO I=1,8
              WRITE(64,REC=IGVP+I) RDES4(I)
              ENDDO
            IGVP=IGVP+8
            DO I=1,6
              WRITE(64,REC=IGVP+I) RID4(I)
              ENDDO
            IGVP=IGVP+6
            DO I=1,6
              WRITE(64,REC=IGVP+I) AID4(I)
              ENDDO
            IGVP=IGVP+6
            ENDIF
          IF(NBYTE.EQ.8) THEN
            DO I=1,4
              WRITE(64,REC=IGVP+I) RDES8(I)
              ENDDO
            IGVP=IGVP+4
            DO I=1,3
              WRITE(64,REC=IGVP+I) RID8(I)
              ENDDO
            IGVP=IGVP+3
            DO I=1,3
              WRITE(64,REC=IGVP+I) AID8(I)
              ENDDO
            IGVP=IGVP+3
            ENDIF
          WRITE(64,REC=IGVP+1) NDSETSV
          WRITE(64,REC=IGVP+2) NE
          WRITE(64,REC=IGVP+3) DT*NSPOOLGV
          WRITE(64,REC=IGVP+4) NSPOOLGV
          WRITE(64,REC=IGVP+5) 2
          IGVP=IGVP+5
          CLOSE(64)                    ! DO THIS TO FLUSH THE WRITE BUFFER
          OPEN(64,FILE=DIRNAME//'/'//'fort.64',
     &         ACCESS='DIRECT',RECL=NBYTE)
          ENDIF

C...
C....INITILIZE GLOBAL WIND and pressure SPOOL COUNTER
C....OPEN GLOBAL WIND and pressure OUTPUT FILEs
C....WRITE OUT HEADER INFORMATION INCLUDING NDSETSW
C....(NO. OF GLOBAL WIND DATA SETS TO BE SPOOLED),
C....NP, DT*NSPOOLGW, NSPOOLGW, IRTYPE
C...
        NSCOUGW=0
        IGWP=0
        igpp=0

        IF(ABS(NOUTGW).EQ.1) THEN
          open(73,file=dirname//'/'//'fort.73')
          write(73,3220) rundes,runid,agrid
          write(73,3645) ndsetsw,np,dtdp*nspoolgw,nspoolgw,1
          igpp=2
          OPEN(74,FILE=DIRNAME//'/'//'fort.74')
          WRITE(74,3220) RUNDES,RUNID,AGRID
          WRITE(74,3645) NDSETSW,NP,DTDP*NSPOOLGW,NSPOOLGW,2
          IGWP=2
          ENDIF

        IF(ABS(NOUTGW).EQ.2) THEN
          open(73,file=dirname//'/'//'fort.73',
     &          access='direct',recl=nbyte)
          OPEN(74,FILE=DIRNAME//'/'//'fort.74',
     &          ACCESS='DIRECT',RECL=NBYTE)
          IF(NBYTE.EQ.4) THEN
            DO I=1,8
              write(73,rec=igpp+i) rdes4(i)
              WRITE(74,REC=IGWP+I) RDES4(I)
              ENDDO
            igpp=igpp+8
            IGWP=IGWP+8
            DO I=1,6
              write(73,rec=igpp+i) rid4(i)
              WRITE(74,REC=IGWP+I) RID4(I)
              ENDDO
            igpp=igpp+6
            IGWP=IGWP+6
            DO I=1,6
              write(73,rec=igpp+i) aid4(i)
              WRITE(74,REC=IGWP+I) AID4(I)
              ENDDO
            igpp=igpp+6
            IGWP=IGWP+6
            ENDIF
          IF(NBYTE.EQ.8) THEN
            DO I=1,4
              write(73,rec=igpp+i) rdes8(i)
              WRITE(74,REC=IGWP+I) RDES8(I)
              ENDDO
            igpp=igpp+4
            IGWP=IGWP+4
            DO I=1,3
              write(73,rec=igpp+i) rid8(i)
              WRITE(74,REC=IGWP+I) RID8(I)
              ENDDO
            igpp=igpp+3
            IGWP=IGWP+3
            DO I=1,3
              write(73,rec=igpp+i) aid8(i)
              WRITE(74,REC=IGWP+I) AID8(I)
              ENDDO
            igpp=igpp+3
            IGWP=IGWP+3
            ENDIF
          write(73,rec=igpp+1) ndsetsw
          write(73,rec=igpp+2) np
          write(73,rec=igpp+3) dt*nspoolgw
          write(73,rec=igpp+4) nspoolgw
          write(73,rec=igpp+5) 2
          igpp=igpp+5
          close(73)                    ! DO THIS TO FLUSH THE WRITE BUFFER
          open(73,file=dirname//'/'//'fort.73',
     &         access='direct',recl=nbyte)
          WRITE(74,REC=IGWP+1) NDSETSW
          WRITE(74,REC=IGWP+2) NP
          WRITE(74,REC=IGWP+3) DT*NSPOOLGW
          WRITE(74,REC=IGWP+4) NSPOOLGW
          WRITE(74,REC=IGWP+5) 2
          IGWP=IGWP+5
          CLOSE(74)                    ! DO THIS TO FLUSH THE WRITE BUFFER
          OPEN(74,FILE=DIRNAME//'/'//'fort.74',
     &         ACCESS='DIRECT',RECL=NBYTE)
          ENDIF
#ifdef SWAN
Casey 101118: Added the output of radiation stress gradients.
       IGRadS=0
       IF(ABS(NOUTGW).EQ.1) THEN
          OPEN(164,FILE=DIRNAME//'/'//'rads.64')
          WRITE(164,3220) RUNDES,RUNID,AGRID
          WRITE(164,3645) NDSETSW,NP,DTDP*NSPOOLGW,NSPOOLGW,2
          IGRadS=2
          ENDIF
       IF(ABS(NOUTGW).EQ.2) THEN
          OPEN(164,FILE=TRIM(LOCALDIR)//'/'//'rads.64',
     &           ACCESS='DIRECT',RECL=NByte)
          IF(NBYTE.EQ.4) THEN
             DO I=1,8
                WRITE(164,REC=IGRadS+I) RDES4(I)
                ENDDO
             IGRadS=IGRadS+8
             DO I=1,6
                WRITE(164,REC=IGRadS+I) RID4(I)
                ENDDO
             IGRadS=IGRadS+6
             DO I=1,6
                WRITE(164,REC=IGRadS+I) AID4(I)
                ENDDO
             IGRadS=IGRadS+6
             ENDIF
          IF(NBYTE.EQ.8) THEN
             DO I=1,4
                WRITE(164,REC=IGRadS+I) RDES8(I)
                ENDDO
             IGRadS=IGRadS+4
             DO I=1,3
                WRITE(164,REC=IGRadS+I) RID8(I)
                ENDDO
             IGRadS=IGRadS+3
             DO I=1,3
                WRITE(164,REC=IGRadS+I) AID8(I)
                ENDDO
             IGRadS=IGRadS+3
             ENDIF
          WRITE(164,REC=IGRadS+1) NDSETSW
          WRITE(164,REC=IGRadS+2) NP
          WRITE(164,REC=IGRadS+3) DT*NSPOOLGW
          WRITE(164,REC=IGRadS+4) NSPOOLGW
          WRITE(164,REC=IGRadS+5) 2
          IGRadS=IGRadS+5
          CLOSE(164)
          ENDIF
#endif

C...
C....INITILIZE GLOBAL CONCENTRATION SPOOL COUNTER
C....OPEN GLOBAL CONCENTRATION OUTPUT FILE
C....WRITE OUT HEADER INFORMATION INCLUDING NDSETSC
C....(NO. OF GLOBAL CONCENTRATION DATA SETS TO BE SPOOLED),
C....NP, DT*NSPOOLGC, NSPOOLGC, IRTYPE
C...
        NSCOUGC=0
        IGCP=0

        IF(ABS(NOUTGC).EQ.1) THEN
          OPEN(83,FILE=DIRNAME//'/'//'fort.83')
          WRITE(83,3220) RUNDES,RUNID,AGRID
          WRITE(83,3645) NDSETSC,NP,DTDP*NSPOOLGC,NSPOOLGC,1
          IGCP=2
          ENDIF

        IF(ABS(NOUTGC).EQ.2) THEN
          OPEN(83,FILE=DIRNAME//'/'//'fort.83',
     &          ACCESS='DIRECT',RECL=NBYTE)
          IF(NBYTE.EQ.4) THEN
            DO I=1,8
              WRITE(83,REC=IGCP+I) RDES4(I)
              ENDDO
            IGCP=IGCP+8
            DO I=1,6
              WRITE(83,REC=IGCP+I) RID4(I)
              ENDDO
            IGCP=IGCP+6
            DO I=1,6
              WRITE(83,REC=IGCP+I) AID4(I)
              ENDDO
            IGCP=IGCP+6
            ENDIF
          IF(NBYTE.EQ.8) THEN
            DO I=1,4
              WRITE(83,REC=IGCP+I) RDES8(I)
              ENDDO
            IGCP=IGCP+4
            DO I=1,3
              WRITE(83,REC=IGCP+I) RID8(I)
              ENDDO
            IGCP=IGCP+3
            DO I=1,3
              WRITE(83,REC=IGCP+I) AID8(I)
              ENDDO
            IGCP=IGCP+3
            ENDIF
          WRITE(83,REC=IGCP+1) NDSETSC
          WRITE(83,REC=IGCP+2) NP
          WRITE(83,REC=IGCP+3) DT*NSPOOLGC
          WRITE(83,REC=IGCP+4) NSPOOLGC
          WRITE(83,REC=IGCP+5) 1
          IGCP=IGCP+5
          CLOSE(83)                    ! DO THIS TO FLUSH THE WRITE BUFFER
          OPEN(83,FILE=DIRNAME//'/'//'fort.83',
     &         ACCESS='DIRECT',RECL=NBYTE)
          ENDIF

C...
C....INITIALIZE HARMONIC ANALYSIS MATRICES, MEAN AND SQUARE VECTORS
C...
        IF (IHARIND.EQ.1) THEN
           ICHA=0
           CALL HACOLDS(HAFREQ)
           IF(NHASE.EQ.1) CALL HACOLDSES(NSTAE)
           IF(NHASV.EQ.1) CALL HACOLDSVS(NSTAV)
           IF(NHAGE.EQ.1) CALL HACOLDSEG(NP)
           IF(NHAGV.EQ.1) CALL HACOLDSVG(NP)
           IF ( CHARMV) THEN
             DO I=1,NP
                ELAV(I)=0.D0
                XVELAV(I)=0.D0
                YVELAV(I)=0.D0
                ELVA(I)=0.D0
                XVELVA(I)=0.D0
                YVELVA(I)=0.D0
             ENDDO
           ENDIF !  charmv
        ENDIF
C
 1112 FORMAT(/,1X,79('_'))
 9883 FORMAT(' !!! NODE ',I6,' DRIED (LANDLOCKING)')

C.....Added sediment transport output files (Ethan Kubatko, 8-1-2003)

C.....INITILIZE GLOBAL BATHYMETRY SPOOL COUNTER
C.....OPEN GLOBAL BATHYMETRY OUTPUT FILE
C.....WRITE OUT HEADER INFORMATION INCLUDING NDSETSE
C.....(NO. OF GLOBAL ELEVATION DATA SETS TO BE SPOOLED),
C.....NP, DT*NSPOOLGE, NSPOOLGE, IRTYPE

        IF (SEDFLAG.GE.1) THEN
        
          NSCOUGE = 0
          IGEP = 0

          IF (ABS(NOUTGE).EQ.1) THEN
            OPEN(84,FILE=DIRNAME//'/'//'fort.84')
            OPEN(85,FILE=DIRNAME//'/'//'fort.85')
            WRITE(84,3220) RUNDES,RUNID,AGRID
            WRITE(84,3645) NDSETSE,NP,DTDP*NSPOOLGE,NSPOOLGE,1
            WRITE(85,3220) RUNDES,RUNID,AGRID
            WRITE(85,3645) NDSETSE,NP,DTDP*NSPOOLGE,NSPOOLGE,1
            IGEP=2
          ENDIF

          IF (ABS(NOUTGE).EQ.2) THEN
            OPEN(84,FILE=DIRNAME//'/'//'fort.84',
     &          ACCESS='DIRECT',RECL=NBYTE)
            OPEN(85,FILE=DIRNAME//'/'//'fort.85',
     &          ACCESS='DIRECT',RECL=NBYTE)
            IF(NBYTE.EQ.4) THEN
              DO I=1,8
                WRITE(84,REC=IGEP+I) RDES4(I)
                WRITE(85,REC=IGEP+I) RDES4(I)
              ENDDO
              IGEP=IGEP+8
              DO I=1,6
                WRITE(84,REC=IGEP+I) RID4(I)
                WRITE(85,REC=IGEP+I) RID4(I)
              ENDDO
              IGEP=IGEP+6
              DO I=1,6
                WRITE(84,REC=IGEP+I) AID4(I)
                WRITE(85,REC=IGEP+I) AID4(I)
              ENDDO
              IGEP=IGEP+6
              ENDIF
            IF(NBYTE.EQ.8) THEN
              DO I=1,4
                WRITE(84,REC=IGEP+I) RDES8(I)
                WRITE(85,REC=IGEP+I) RDES8(I)
              ENDDO
              IGEP=IGEP+4
              DO I=1,3
                WRITE(84,REC=IGEP+I) RID8(I)
                WRITE(85,REC=IGEP+I) RID8(I)
              ENDDO
              IGEP=IGEP+3
              DO I=1,3
                WRITE(84,REC=IGEP+I) AID8(I)
                WRITE(85,REC=IGEP+I) AID8(I)
              ENDDO
              IGEP=IGEP+3
            ENDIF
            WRITE(84,REC=IGEP+1) NDSETSE
            WRITE(84,REC=IGEP+2) NP
            WRITE(84,REC=IGEP+3) DT*NSPOOLGE
            WRITE(84,REC=IGEP+4) NSPOOLGE
            WRITE(84,REC=IGEP+5) 1
            WRITE(85,REC=IGEP+1) NDSETSE
            WRITE(85,REC=IGEP+2) NP
            WRITE(85,REC=IGEP+3) DT*NSPOOLGE
            WRITE(85,REC=IGEP+4) NSPOOLGE
            WRITE(85,REC=IGEP+5) 1
            IGEP=IGEP+5
            CLOSE(84)                    ! DO THIS TO FLUSH THE WRITE BUFFER
            CLOSE(85)                    ! DO THIS TO FLUSH THE WRITE BUFFER
            OPEN(84,FILE=DIRNAME//'/'//'fort.84',
     &         ACCESS='DIRECT',RECL=NBYTE)
            OPEN(85,FILE=DIRNAME//'/'//'fort.85',
     &         ACCESS='DIRECT',RECL=NBYTE)
          ENDIF
C...
C.....INITILIZE BATHYMETRY STATION SPOOL COUNTER
C.....OPEN ELEVATION STATION OUTPUT FILE
C.....WRITE OUT HEADER INFORMATION INCLUDING NTRSPE (NO. OF DATA PTS. AT EACH
C.....ELEVATION STATION), NSTAE, DT*NSPOOLE, NSPOOLE, IRTYPE
C...
          NSCOUE=0
          IESTP=0

          IF (ABS(NOUTE).EQ.1) THEN
            OPEN(82,FILE=DIRNAME//'/'//'fort.82')
            WRITE(82,3220) RUNDES,RUNID,AGRID
            WRITE(82,3645) NTRSPE,NSTAE,DTDP*NSPOOLE,NSPOOLE,1
            IESTP=2
          ENDIF

          IF (ABS(NOUTE).EQ.2) THEN
            OPEN(82,FILE=DIRNAME//'/'//'fort.82',
     &          ACCESS='DIRECT',RECL=NBYTE)
            IF(NBYTE.EQ.4) THEN
              DO I=1,8
                WRITE(82,REC=IESTP+I) RDES4(I)
              ENDDO
              IESTP=IESTP+8
              DO I=1,6
                WRITE(82,REC=IESTP+I) RID4(I)
              ENDDO
              IESTP=IESTP+6
              DO I=1,6
                WRITE(82,REC=IESTP+I) AID4(I)
              ENDDO
              IESTP=IESTP+6
            ENDIF
            IF(NBYTE.EQ.8) THEN
              DO I=1,4
                WRITE(82,REC=IESTP+I) RDES8(I)
              ENDDO
              IESTP=IESTP+4
              DO I=1,3
                WRITE(82,REC=IESTP+I) RID8(I)
              ENDDO
              IESTP=IESTP+3
              DO I=1,3
                WRITE(82,REC=IESTP+I) AID8(I)
              ENDDO
              IESTP=IESTP+3
            ENDIF
            WRITE(82,REC=IESTP+1) NTRSPE
            WRITE(82,REC=IESTP+2) NSTAE
            WRITE(82,REC=IESTP+3) DT*NSPOOLE
            WRITE(82,REC=IESTP+4) NSPOOLE
            WRITE(82,REC=IESTP+5) 1
            IESTP=IESTP+5
            CLOSE(82)                    ! DO THIS TO FLUSH THE WRITE BUFFER
            OPEN(82,FILE=DIRNAME//'/'//'fort.82',
     &         ACCESS='DIRECT',RECL=NBYTE)
          ENDIF

C.....INITILIZE GLOBAL BED LOAD SEDIMENT FLUX SPOOL COUNTER
C.....OPEN GLOBAL SEDIMENT FLUX OUTPUT FILE
C.....WRITE OUT HEADER INFORMATION INCLUDING NDSETSV
C.....(NO. OF GLOBAL VELOCITY DATA SETS TO BE SPOOLED),
C.....NP, DT*NSPOOLGV, NSPOOLGV, IRTYPE
C...
          NSCOUGV=0
          IGVP=0

          IF(ABS(NOUTGV).EQ.1) THEN
            OPEN(94,FILE=DIRNAME//'/'//'fort.94')
            WRITE(94,3220) RUNDES,RUNID,AGRID
            WRITE(94,3645) NDSETSV,NP,DTDP*NSPOOLGV,NSPOOLGV,2
            IGVP=2
          ENDIF

          IF(ABS(NOUTGV).EQ.2) THEN
            OPEN(94,FILE=DIRNAME//'/'//'fort.94',
     &          ACCESS='DIRECT',RECL=NBYTE)
            IF(NBYTE.EQ.4) THEN
              DO I=1,8
                WRITE(94,REC=IGVP+I) RDES4(I)
              ENDDO
              IGVP=IGVP+8
              DO I=1,6
                WRITE(94,REC=IGVP+I) RID4(I)
              ENDDO
              IGVP=IGVP+6
              DO I=1,6
                WRITE(94,REC=IGVP+I) AID4(I)
              ENDDO
              IGVP=IGVP+6
            ENDIF
            IF(NBYTE.EQ.8) THEN
              DO I=1,4
                WRITE(94,REC=IGVP+I) RDES8(I)
              ENDDO
              IGVP=IGVP+4
              DO I=1,3
                WRITE(94,REC=IGVP+I) RID8(I)
              ENDDO
              IGVP=IGVP+3
              DO I=1,3
                WRITE(94,REC=IGVP+I) AID8(I)
              ENDDO
              IGVP=IGVP+3
            ENDIF
            WRITE(94,REC=IGVP+1) NDSETSV
            WRITE(94,REC=IGVP+2) NP
            WRITE(94,REC=IGVP+3) DT*NSPOOLGV
            WRITE(94,REC=IGVP+4) NSPOOLGV
            WRITE(94,REC=IGVP+5) 2
            IGVP=IGVP+5
            CLOSE(94)                    ! DO THIS TO FLUSH THE WRITE BUFFER
            OPEN(94,FILE=DIRNAME//'/'//'fort.94',
     &         ACCESS='DIRECT',RECL=NBYTE)
          ENDIF
          
C.....INITILIZE GLOBAL SUSPENDED LOAD SEDIMENT FLUX SPOOL COUNTER
C.....OPEN GLOBAL SEDIMENT FLUX OUTPUT FILE
C.....WRITE OUT HEADER INFORMATION INCLUDING NDSETSV
C.....(NO. OF GLOBAL VELOCITY DATA SETS TO BE SPOOLED),
C.....NP, DT*NSPOOLGV, NSPOOLGV, IRTYPE
C...
          NSCOUGV=0
          IGVP=0

          IF(ABS(NOUTGV).EQ.1) THEN
            OPEN(96,FILE=DIRNAME//'/'//'fort.96')
            WRITE(96,3220) RUNDES,RUNID,AGRID
            WRITE(96,3645) NDSETSV,NP,DTDP*NSPOOLGV,NSPOOLGV,2
            IGVP=2
          ENDIF

          IF(ABS(NOUTGV).EQ.2) THEN
            OPEN(96,FILE=DIRNAME//'/'//'fort.96',
     &          ACCESS='DIRECT',RECL=NBYTE)
            IF(NBYTE.EQ.4) THEN
              DO I=1,8
                WRITE(96,REC=IGVP+I) RDES4(I)
              ENDDO
              IGVP=IGVP+8
              DO I=1,6
                WRITE(96,REC=IGVP+I) RID4(I)
              ENDDO
              IGVP=IGVP+6
              DO I=1,6
                WRITE(96,REC=IGVP+I) AID4(I)
              ENDDO
              IGVP=IGVP+6
            ENDIF
            IF(NBYTE.EQ.8) THEN
              DO I=1,4
                WRITE(96,REC=IGVP+I) RDES8(I)
              ENDDO
              IGVP=IGVP+4
              DO I=1,3
                WRITE(96,REC=IGVP+I) RID8(I)
              ENDDO
              IGVP=IGVP+3
              DO I=1,3
                WRITE(96,REC=IGVP+I) AID8(I)
              ENDDO
              IGVP=IGVP+3
            ENDIF
            WRITE(96,REC=IGVP+1) NDSETSV
            WRITE(96,REC=IGVP+2) NP
            WRITE(96,REC=IGVP+3) DT*NSPOOLGV
            WRITE(96,REC=IGVP+4) NSPOOLGV
            WRITE(96,REC=IGVP+5) 2
            IGVP=IGVP+5
            CLOSE(96)                    ! DO THIS TO FLUSH THE WRITE BUFFER
            OPEN(96,FILE=DIRNAME//'/'//'fort.96',
     &         ACCESS='DIRECT',RECL=NBYTE)
          ENDIF

        ENDIF
          
      RETURN
      END
