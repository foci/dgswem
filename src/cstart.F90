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
        SUBROUTINE COLDSTART(s,global_here)
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
      type (global_type) :: global_here

      INTEGER i,j
      
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
          WRITE(s%fort16unit,1112)
          IF(global_here%NSCREEN.EQ.1.AND.S%MYPROC.EQ.0) WRITE(6,1977)
          WRITE(s%fort16unit,1977)
 1977     FORMAT(/,1X,'ELEVATION SPECIFIED INFORMATION READ FROM UNIT ',&
           '19',/)
          OPEN(s%fort19unit,FILE=S%DIRNAME//'/'//'fort.19')
          READ(s%fort19unit,*) global_here%ETIMINC
          DO J=1,global_here%NETA
             READ(s%fort19unit,*) global_here%ESBIN1(J)
          END DO
          DO J=1,global_here%NETA
             READ(s%fort19unit,*) global_here%ESBIN2(J)
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
          WRITE(s%fort16unit,1112)
          IF(global_here%NSCREEN.EQ.1.AND.S%MYPROC.EQ.0) WRITE(6,1979)
          WRITE(s%fort16unit,1979)
 1979     FORMAT(/,1X,'NORMAL FLOW INFORMATION READ FROM UNIT 20',/)
          OPEN(s%fort20unit,FILE=S%DIRNAME//'/'//'fort.20')
          READ(s%fort20unit,*) global_here%FTIMINC
          DO J=1,global_here%NVEL
            global_here%QNIN1(J)=0.D0
            IF((global_here%LBCODEI(J).EQ.2).OR.(global_here%LBCODEI(J).EQ.12)&
                          .OR.(global_here%LBCODEI(J).EQ.22))&
                                          READ(s%fort20unit,*) global_here%QNIN1(J)
            END DO
          DO J=1,global_here%NVEL
            global_here%QNIN2(J)=0.D0
            IF((global_here%LBCODEI(J).EQ.2).OR.(global_here%LBCODEI(J).EQ.12)&
                          .OR.(global_here%LBCODEI(J).EQ.22))&
                                          READ(s%fort20unit,*) global_here%QNIN2(J)
            END DO
          global_here%QTIME1 = global_here%STATIM*86400.D0
          global_here%QTIME2 = global_here%QTIME1 + global_here%FTIMINC
          ENDIF

!...INPUT METEOROLOGICAL INFORMATION FROM UNIT 22 OR UNIT 200 SERIES
!....IF FLEET NUMERIC WIND DATA IS USED, FIND BEGINNING TIME IN FILE,
!....NOTE: CAN'T DEAL WITH WIND THAT STARTS AFTER global_here%WREFTIM!!!!!!!!!!!!
!....READ IN AND INTERPOLATE IN SPACE ONTO THE ADCIRC GRID THE
!....TIME LEVEL 1 AND LEVEL 2 WIND FIELDS

        IF(global_here%NWS.NE.0) THEN
          DO I=1,global_here%NP
            global_here%WSX1(I)=0.D0
            global_here%WSY1(I)=0.D0
            global_here%PR1(I) =0.D0
            global_here%WSX2(I)=0.D0
            global_here%WSY2(I)=0.D0
            global_here%PR2(I) =0.D0
            ENDDO

          IF(global_here%NSCREEN.EQ.1.AND.S%MYPROC.EQ.0) WRITE(6,1112)
          WRITE(s%fort16unit,1112)
          IF(global_here%NSCREEN.EQ.1.AND.S%MYPROC.EQ.0) WRITE(6,1980)
          WRITE(s%fort16unit,1980)
 1980     FORMAT(/,1X,'WIND (AND PRESSURE) INFORMATION READ.',/)
          ENDIF

        IF(global_here%NWS.EQ.1) THEN
          OPEN(s%fort22unit,FILE=S%DIRNAME//'/'//'fort.22')
          ENDIF

        IF(ABS(global_here%NWS).EQ.2) THEN
          OPEN(s%fort22unit,FILE=S%DIRNAME//'/'//'fort.22')
          READ(s%fort22unit,*) (global_here%NHG,global_here%WVNX1(I),global_here%WVNY1(I),global_here%PRN1(I),I=1,global_here%NP)
          READ(s%fort22unit,*) (global_here%NHG,global_here%WVNX2(I),global_here%WVNY2(I),global_here%PRN2(I),I=1,global_here%NP)
          global_here%WTIME1 = global_here%STATIM*86400.D0
          global_here%WTIME2 = global_here%WTIME1 + global_here%WTIMINC
          ENDIF

        IF(global_here%NWS.EQ.3) THEN
          OPEN(s%fort22unit,FILE=S%DIRNAME//'/'//'fort.22')
 2222     CALL NWS3GET(s, global_here%X, global_here%Y, global_here%SLAM, global_here%SFEA, global_here%WVNX2, global_here%WVNY2, global_here%IWTIME, global_here%IWYR,&
                  global_here%WTIMED, global_here%NP, global_here%NWLON, global_here%NWLAT, global_here%WLATMAX, global_here%WLONMIN,&
                  global_here%WLATINC, global_here%WLONINC, global_here%ICS, global_here%NSCREEN, global_here%ScreenUnit )
          IF(global_here%IWYR.NE.global_here%IREFYR) THEN
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
          WRITE(s%fort16unit,*) 'FOUND WIND DATA AT TIME= ',global_here%IWTIMEP
          IF(global_here%NSCREEN.EQ.1.AND.S%MYPROC.EQ.0) &
         WRITE(6,*)'FOUND WIND DATA AT TIME= ',global_here%IWTIME
          WRITE(s%fort16unit,*) 'FOUND WIND DATA AT TIME= ',global_here%IWTIME
          global_here%WTIME2=global_here%WTIMED-global_here%WREFTIM                  !CAST INTO MODEL TIME REFRENCE
          global_here%WTIME1=global_here%WTIME2-global_here%WTIMINC
          ENDIF

        IF(ABS(global_here%NWS).EQ.4) THEN
          OPEN(s%fort22unit,FILE=S%DIRNAME//'/'//'fort.22')
          global_here%WTIME1 = global_here%STATIM*86400.D0
          global_here%WTIME2=global_here%WTIME1+global_here%WTIMINC
          CALL NWS4GET(global_here%WVNX1,global_here%WVNY1,global_here%PRN1,global_here%NP,global_here%RHOWAT0,global_here%G)
          CALL NWS4GET(global_here%WVNX2,global_here%WVNY2,global_here%PRN2,global_here%NP,global_here%RHOWAT0,global_here%G)
          ENDIF

        IF(ABS(global_here%NWS).EQ.5) THEN
          OPEN(s%fort22unit,FILE=S%DIRNAME//'/'//'fort.22')
          READ(s%fort22unit,*) (global_here%NHG,global_here%WVNX1(I),global_here%WVNY1(I),global_here%PRN1(I),I=1,global_here%NP)
          READ(s%fort22unit,*) (global_here%NHG,global_here%WVNX2(I),global_here%WVNY2(I),global_here%PRN2(I),I=1,global_here%NP)
          global_here%WTIME1 = global_here%STATIM*86400.D0
          global_here%WTIME2 = global_here%WTIME1 + global_here%WTIMINC
          ENDIF

        IF(global_here%NWS.EQ.6) THEN
          OPEN(s%fort22unit,FILE=S%DIRNAME//'/'//'fort.22')
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
          CALL NWS10GET(s,global_here,global_here%NWSGGWI,global_here%SLAM,global_here%SFEA,global_here%WVNX1,global_here%WVNY1,global_here%PRN1,global_here%NP,global_here%RHOWAT0,global_here%G,&
                  global_here%NWLON,global_here%NWLAT,global_here%WTIMINC) !JUST COMPUTE INTERPOLATING FACTORS
          global_here%NWSGGWI=1
          CALL NWS10GET(s,global_here,global_here%NWSGGWI,global_here%SLAM,global_here%SFEA,global_here%WVNX2,global_here%WVNY2,global_here%PRN2,global_here%NP,global_here%RHOWAT0,global_here%G,&
                  global_here%NWLON,global_here%NWLAT,global_here%WTIMINC) !NOW INTERPOLATE 1st WIND FIELD
          ENDIF

        IF(global_here%NWS.EQ.11) THEN
          global_here%WTIME1=global_here%STATIM*86400.D0
          global_here%WTIME2=global_here%WTIME1+global_here%WTIMINC
          global_here%NWSEGWI=0
          global_here%IDSETFLG=0
          IF(global_here%NSCREEN.EQ.1.AND.S%MYPROC.EQ.0) WRITE(6,1197)
          WRITE(s%fort16unit,1197)
 1197     FORMAT(/,1X,'THE E29 MET GRID INTERPOLATING FACTORS ARE ',&
                'BEING COMPUTED ')
          CALL NWS11GET(s,global_here,global_here%NWSEGWI,global_here%IDSETFLG,global_here%SLAM,global_here%SFEA,global_here%WVNX1,global_here%WVNY1,global_here%PRN1,global_here%NP,&
                  global_here%RHOWAT0,global_here%G)  !JUST COMPUTE INTERPOLATING FACTORS
          IF(global_here%NSCREEN.EQ.1.AND.S%MYPROC.EQ.0) WRITE(6,1198)
          WRITE(s%fort16unit,1198)
 1198     FORMAT(1X,'FINISHED COMPUTING E29 INTERPOLATING FACTORS',/)
          global_here%NWSEGWI=1
          global_here%IDSETFLG=1
          CALL NWS11GET(s,global_here,global_here%NWSEGWI,global_here%IDSETFLG,global_here%SLAM,global_here%SFEA,global_here%WVNX2,global_here%WVNY2,global_here%PRN2,global_here%NP,&
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
          OPEN(s%fort23unit,FILE=S%DIRNAME//'/'//'fort.23')
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
          WRITE(s%fort16unit,1112)
          IF(global_here%NSCREEN.EQ.1.AND.S%MYPROC.EQ.0) WRITE(6,1981)
          WRITE(s%fort16unit,1981)
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
              WRITE(s%fort16unit,9883) I
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
          OPEN(s%fort61unit,FILE=S%DIRNAME//'/'//'fort.61')
          WRITE(s%fort61unit,3220) global_here%RUNDES,global_here%RUNID,global_here%AGRID
          WRITE(s%fort61unit,3645) global_here%NTRSPE,global_here%NSTAE,global_here%DTDP*global_here%NSPOOLE,global_here%NSPOOLE,1
          global_here%IESTP=2
          ENDIF

          IF(ABS(global_here%NOUTE).EQ.2) THEN
             OPEN(s%fort61unit,FILE=S%DIRNAME//'/'//'fort.61',&
                  ACCESS='DIRECT',RECL=NBYTE)
             IF(NBYTE.EQ.4) THEN
                DO I=1,8
                   WRITE(s%fort61unit,REC=global_here%IESTP+I) global_here%RDES4(I)
                ENDDO
                global_here%IESTP=global_here%IESTP+8
                DO I=1,6
                   WRITE(s%fort61unit,REC=global_here%IESTP+I) global_here%RID4(I)
                ENDDO
                global_here%IESTP=global_here%IESTP+6
                DO I=1,6
                   WRITE(s%fort61unit,REC=global_here%IESTP+I) global_here%AID4(I)
                ENDDO
                global_here%IESTP=global_here%IESTP+6
             ENDIF
             IF(NBYTE.EQ.8) THEN
                DO I=1,4
                   WRITE(s%fort61unit,REC=global_here%IESTP+I) global_here%RDES8(I)
                ENDDO
                global_here%IESTP=global_here%IESTP+4
                DO I=1,3
                   WRITE(s%fort61unit,REC=global_here%IESTP+I) global_here%RID8(I)
                ENDDO
                global_here%IESTP=global_here%IESTP+3
                DO I=1,3
                   WRITE(s%fort61unit,REC=global_here%IESTP+I) global_here%AID8(I)
                ENDDO
                global_here%IESTP=global_here%IESTP+3
             ENDIF
             WRITE(s%fort61unit,REC=global_here%IESTP+1) global_here%NTRSPE
             WRITE(s%fort61unit,REC=global_here%IESTP+2) global_here%NSTAE
             WRITE(s%fort61unit,REC=global_here%IESTP+3) global_here%DT*global_here%NSPOOLE
             WRITE(s%fort61unit,REC=global_here%IESTP+4) global_here%NSPOOLE
             WRITE(s%fort61unit,REC=global_here%IESTP+5) 1
             global_here%IESTP=global_here%IESTP+5
             CLOSE(s%fort61unit)                    ! DO THIS TO FLUSH THE WRITE BUFFER
             OPEN(s%fort61unit,FILE=S%DIRNAME//'/'//'fort.61',&
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
          OPEN(s%fort62unit,FILE=S%DIRNAME//'/'//'fort.62')
          WRITE(s%fort62unit,3220) global_here%RUNDES,global_here%RUNID,global_here%AGRID
          WRITE(s%fort62unit,3645) global_here%NTRSPV,global_here%NSTAV,global_here%DTDP*global_here%NSPOOLV,global_here%NSPOOLV,2
          global_here%IVSTP=2
          ENDIF

        IF(ABS(global_here%NOUTV).EQ.2) THEN
          OPEN(s%fort62unit,FILE=S%DIRNAME//'/'//'fort.62',&
          ACCESS='DIRECT',RECL=NBYTE)
          IF(NBYTE.EQ.4) THEN
            DO I=1,8
              WRITE(s%fort62unit,REC=global_here%IVSTP+I) global_here%RDES4(I)
              ENDDO
            global_here%IVSTP=global_here%IVSTP+8
            DO I=1,6
              WRITE(s%fort62unit,REC=global_here%IVSTP+I) global_here%RID4(I)
              ENDDO
            global_here%IVSTP=global_here%IVSTP+6
            DO I=1,6
              WRITE(s%fort62unit,REC=global_here%IVSTP+I) global_here%AID4(I)
              ENDDO
            global_here%IVSTP=global_here%IVSTP+6
            ENDIF
          IF(NBYTE.EQ.8) THEN
            DO I=1,4
              WRITE(s%fort62unit,REC=global_here%IVSTP+I) global_here%RDES8(I)
              ENDDO
            global_here%IVSTP=global_here%IVSTP+4
            DO I=1,3
              WRITE(s%fort62unit,REC=global_here%IVSTP+I) global_here%RID8(I)
              ENDDO
            global_here%IVSTP=global_here%IVSTP+3
            DO I=1,3
              WRITE(s%fort62unit,REC=global_here%IVSTP+I) global_here%AID8(I)
              ENDDO
            global_here%IVSTP=global_here%IVSTP+3
            ENDIF
          WRITE(s%fort62unit,REC=global_here%IVSTP+1) global_here%NTRSPV
          WRITE(s%fort62unit,REC=global_here%IVSTP+2) global_here%NSTAV
          WRITE(s%fort62unit,REC=global_here%IVSTP+3) global_here%DT*global_here%NSPOOLV
          WRITE(s%fort62unit,REC=global_here%IVSTP+4) global_here%NSPOOLV
          WRITE(s%fort62unit,REC=global_here%IVSTP+5) 2
          global_here%IVSTP=global_here%IVSTP+5
          CLOSE(s%fort62unit)                    ! DO THIS TO FLUSH THE WRITE BUFFER
          OPEN(s%fort62unit,FILE=S%DIRNAME//'/'//'fort.62',&
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
          OPEN(s%fort81unit,FILE=S%DIRNAME//'/'//'fort.81')
          WRITE(s%fort81unit,3220) global_here%RUNDES,global_here%RUNID,global_here%AGRID
          WRITE(s%fort81unit,3645) global_here%NTRSPC,global_here%NSTAC,global_here%DTDP*global_here%NSPOOLC,global_here%NSPOOLC,1
          global_here%ICSTP=2
          ENDIF

        IF(ABS(global_here%NOUTC).EQ.2) THEN
          OPEN(s%fort81unit,FILE=S%DIRNAME//'/'//'fort.81',&
          ACCESS='DIRECT',RECL=NBYTE)
          IF(NBYTE.EQ.4) THEN
            DO I=1,8
              WRITE(s%fort81unit,REC=global_here%ICSTP+I) global_here%RDES4(I)
              ENDDO
            global_here%ICSTP=global_here%ICSTP+8
            DO I=1,6
              WRITE(s%fort81unit,REC=global_here%ICSTP+I) global_here%RID4(I)
              ENDDO
            global_here%ICSTP=global_here%ICSTP+6
            DO I=1,6
              WRITE(s%fort81unit,REC=global_here%ICSTP+I) global_here%AID4(I)
              ENDDO
            global_here%ICSTP=global_here%ICSTP+6
            ENDIF
          IF(NBYTE.EQ.8) THEN
            DO I=1,4
              WRITE(s%fort81unit,REC=global_here%ICSTP+I) global_here%RDES8(I)
              ENDDO
            global_here%ICSTP=global_here%ICSTP+4
            DO I=1,3
              WRITE(s%fort81unit,REC=global_here%ICSTP+I) global_here%RID8(I)
              ENDDO
            global_here%ICSTP=global_here%ICSTP+3
            DO I=1,3
              WRITE(s%fort81unit,REC=global_here%ICSTP+I) global_here%AID8(I)
              ENDDO
            global_here%ICSTP=global_here%ICSTP+3
            ENDIF
          WRITE(s%fort81unit,REC=global_here%ICSTP+1) global_here%NTRSPC
          WRITE(s%fort81unit,REC=global_here%ICSTP+2) global_here%NSTAC
          WRITE(s%fort81unit,REC=global_here%ICSTP+3) global_here%DT*global_here%NSPOOLC
          WRITE(s%fort81unit,REC=global_here%ICSTP+4) global_here%NSPOOLC
          WRITE(s%fort81unit,REC=global_here%ICSTP+5) 1
          global_here%ICSTP=global_here%ICSTP+5
          CLOSE(s%fort81unit)                    ! DO THIS TO FLUSH THE WRITE BUFFER
          OPEN(s%fort81unit,FILE=S%DIRNAME//'/'//'fort.81',&
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
          OPEN(s%fort82unit,FILE=S%DIRNAME//'/'//'fort.82')
          WRITE(s%fort82unit,3220) global_here%RUNDES,global_here%RUNID,global_here%AGRID
          WRITE(s%fort82unit,3645) global_here%NTRSPE,global_here%NSTAE,global_here%DTDP*global_here%NSPOOLE,global_here%NSPOOLE,1
          global_here%IESTP=2
          ENDIF

        IF(ABS(global_here%NOUTE).EQ.2) THEN
          OPEN(s%fort82unit,FILE=S%DIRNAME//'/'//'fort.82',&
          ACCESS='DIRECT',RECL=NBYTE)
          IF(NBYTE.EQ.4) THEN
            DO I=1,8
              WRITE(s%fort82unit,REC=global_here%IESTP+I) global_here%RDES4(I)
              ENDDO
            global_here%IESTP=global_here%IESTP+8
            DO I=1,6
              WRITE(s%fort82unit,REC=global_here%IESTP+I) global_here%RID4(I)
              ENDDO
            global_here%IESTP=global_here%IESTP+6
            DO I=1,6
              WRITE(s%fort82unit,REC=global_here%IESTP+I) global_here%AID4(I)
              ENDDO
            global_here%IESTP=global_here%IESTP+6
            ENDIF
          IF(NBYTE.EQ.8) THEN
            DO I=1,4
              WRITE(s%fort82unit,REC=global_here%IESTP+I) global_here%RDES8(I)
              ENDDO
            global_here%IESTP=global_here%IESTP+4
            DO I=1,3
              WRITE(s%fort82unit,REC=global_here%IESTP+I) global_here%RID8(I)
              ENDDO
            global_here%IESTP=global_here%IESTP+3
            DO I=1,3
              WRITE(s%fort82unit,REC=global_here%IESTP+I) global_here%AID8(I)
              ENDDO
            global_here%IESTP=global_here%IESTP+3
            ENDIF
          WRITE(s%fort82unit,REC=global_here%IESTP+1) global_here%NTRSPE
          WRITE(s%fort82unit,REC=global_here%IESTP+2) global_here%NSTAE
          WRITE(s%fort82unit,REC=global_here%IESTP+3) global_here%DT*global_here%NSPOOLE
          WRITE(s%fort82unit,REC=global_here%IESTP+4) global_here%NSPOOLE
          WRITE(s%fort82unit,REC=global_here%IESTP+5) 1
          global_here%IESTP=global_here%IESTP+5
          CLOSE(s%fort82unit)                    ! DO THIS TO FLUSH THE WRITE BUFFER
          OPEN(s%fort82unit,FILE=S%DIRNAME//'/'//'fort.82',&
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
          OPEN(s%fort17unit,FILE=S%DIRNAME//'/'//'fort.71')
          WRITE(s%fort17unit,3220) global_here%RUNDES,global_here%RUNID,global_here%AGRID
          WRITE(s%fort17unit,3645) global_here%NTRSPM,global_here%NSTAM,global_here%DTDP*global_here%NSPOOLM,global_here%NSPOOLM,1
          global_here%IPSTP=2
          OPEN(s%fort72unit,FILE=S%DIRNAME//'/'//'fort.72')
          WRITE(s%fort72unit,3220) global_here%RUNDES,global_here%RUNID,global_here%AGRID
          WRITE(s%fort72unit,3645) global_here%NTRSPM,global_here%NSTAM,global_here%DTDP*global_here%NSPOOLM,global_here%NSPOOLM,2
          global_here%IWSTP=2
          ENDIF

        IF(ABS(global_here%NOUTM).EQ.2) THEN
          OPEN(s%fort17unit,FILE=S%DIRNAME//'/'//'fort.71',&
          ACCESS='DIRECT',RECL=NBYTE)
          OPEN(s%fort72unit,FILE=S%DIRNAME//'/'//'fort.72',&
          ACCESS='DIRECT',RECL=NBYTE)
          IF(NBYTE.EQ.4) THEN
            DO I=1,8
              WRITE(s%fort17unit,REC=global_here%IPSTP+I) global_here%RDES4(I)
              WRITE(s%fort72unit,REC=global_here%IWSTP+I) global_here%RDES4(I)
              ENDDO
            global_here%IPSTP=global_here%IPSTP+8
            global_here%IWSTP=global_here%IWSTP+8
            DO I=1,6
              WRITE(s%fort17unit,REC=global_here%IPSTP+I) global_here%RID4(I)
              WRITE(s%fort72unit,REC=global_here%IWSTP+I) global_here%RID4(I)
              ENDDO
            global_here%IPSTP=global_here%IPSTP+6
            global_here%IWSTP=global_here%IWSTP+6
            DO I=1,6
              WRITE(s%fort17unit,REC=global_here%IPSTP+I) global_here%AID4(I)
              WRITE(s%fort72unit,REC=global_here%IWSTP+I) global_here%AID4(I)
              ENDDO
            global_here%IPSTP=global_here%IPSTP+6
            global_here%IWSTP=global_here%IWSTP+6
            ENDIF
          IF(NBYTE.EQ.8) THEN
            DO I=1,4
              WRITE(s%fort17unit,REC=global_here%IPSTP+I) global_here%RDES8(I)
              WRITE(s%fort72unit,REC=global_here%IWSTP+I) global_here%RDES8(I)
              ENDDO
            global_here%IPSTP=global_here%IPSTP+4
            global_here%IWSTP=global_here%IWSTP+4
            DO I=1,3
              WRITE(s%fort17unit,REC=global_here%IPSTP+I) global_here%RID8(I)
              WRITE(s%fort72unit,REC=global_here%IWSTP+I) global_here%RID8(I)
              ENDDO
            global_here%IPSTP=global_here%IPSTP+3
            global_here%IWSTP=global_here%IWSTP+3
            DO I=1,3
              WRITE(s%fort17unit,REC=global_here%IPSTP+I) global_here%AID8(I)
              WRITE(s%fort72unit,REC=global_here%IWSTP+I) global_here%AID8(I)
              ENDDO
            global_here%IPSTP=global_here%IPSTP+3
            global_here%IWSTP=global_here%IWSTP+3
            ENDIF
          WRITE(s%fort17unit,REC=global_here%IPSTP+1) global_here%NTRSPM
          WRITE(s%fort17unit,REC=global_here%IPSTP+2) global_here%NSTAM
          WRITE(s%fort17unit,REC=global_here%IPSTP+3) global_here%DT*global_here%NSPOOLM
          WRITE(s%fort17unit,REC=global_here%IPSTP+4) global_here%NSPOOLM
          WRITE(s%fort17unit,REC=global_here%IPSTP+5) 1
          WRITE(s%fort72unit,REC=global_here%IWSTP+1) global_here%NTRSPM
          WRITE(s%fort72unit,REC=global_here%IWSTP+2) global_here%NSTAM
          WRITE(s%fort72unit,REC=global_here%IWSTP+3) global_here%DT*global_here%NSPOOLM
          WRITE(s%fort72unit,REC=global_here%IWSTP+4) global_here%NSPOOLM
          WRITE(s%fort72unit,REC=global_here%IWSTP+5) 2
          global_here%IPSTP=global_here%IPSTP+5
          global_here%IWSTP=global_here%IWSTP+5
          CLOSE(s%fort17unit)                    ! DO THIS TO FLUSH THE WRITE BUFFER
          CLOSE(s%fort72unit)                    ! DO THIS TO FLUSH THE WRITE BUFFER
          OPEN(s%fort17unit,FILE=S%DIRNAME//'/'//'fort.71',&
         ACCESS='DIRECT',RECL=NBYTE)
          OPEN(s%fort72unit,FILE=S%DIRNAME//'/'//'fort.72',&
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
          OPEN(s%fort63unit,FILE=S%DIRNAME//'/'//'fort.63')
          WRITE(s%fort63unit,3220) global_here%RUNDES,global_here%RUNID,global_here%AGRID
          WRITE(s%fort63unit,3645) global_here%NDSETSE,global_here%NE,global_here%DTDP*global_here%NSPOOLGE,global_here%NSPOOLGE,1
          global_here%IGEP=2
          ENDIF

        IF(ABS(global_here%NOUTGE).EQ.1) THEN
          OPEN(s%fort88unit,FILE=S%DIRNAME//'/'//'fort.88')
          WRITE(s%fort88unit,3220) global_here%RUNDES,global_here%RUNID,global_here%AGRID
          WRITE(s%fort88unit,3645) global_here%NDSETSE,global_here%NE,global_here%DTDP*global_here%NSPOOLGE,global_here%NSPOOLGE,1
          global_here%IGEP=2
          ENDIF

        IF(ABS(global_here%NOUTGE).EQ.1) THEN
          OPEN(s%fort89unit,FILE=S%DIRNAME//'/'//'fort.89')
          WRITE(s%fort89unit,3220) global_here%RUNDES,global_here%RUNID,global_here%AGRID
          WRITE(s%fort89unit,3645) global_here%NDSETSE,global_here%NE,global_here%DTDP*global_here%NSPOOLGE,global_here%NSPOOLGE,1
          global_here%IGEP=2
          ENDIF

        IF(ABS(global_here%NOUTGE).EQ.1) THEN
          OPEN(s%fort4lunit,FILE=S%DIRNAME//'/'//'fort.4l')
          WRITE(s%fort4lunit,3220) global_here%RUNDES,global_here%RUNID,global_here%AGRID
          WRITE(s%fort4lunit,3645) global_here%NDSETSE,global_here%NE,global_here%DTDP*global_here%NSPOOLGE,global_here%NSPOOLGE,1
          global_here%IGEP=2
          ENDIF

        IF(ABS(global_here%NOUTGE).EQ.2) THEN
          OPEN(s%fort63unit,FILE=S%DIRNAME//'/'//'fort.63',&
          ACCESS='DIRECT',RECL=NBYTE)
          IF(NBYTE.EQ.4) THEN
            DO I=1,8
              WRITE(s%fort63unit,REC=global_here%IGEP+I) global_here%RDES4(I)
              ENDDO
            global_here%IGEP=global_here%IGEP+8
            DO I=1,6
              WRITE(s%fort63unit,REC=global_here%IGEP+I) global_here%RID4(I)
              ENDDO
            global_here%IGEP=global_here%IGEP+6
            DO I=1,6
              WRITE(s%fort63unit,REC=global_here%IGEP+I) global_here%AID4(I)
              ENDDO
            global_here%IGEP=global_here%IGEP+6
            ENDIF
          IF(NBYTE.EQ.8) THEN
            DO I=1,4
              WRITE(s%fort63unit,REC=global_here%IGEP+I) global_here%RDES8(I)
              ENDDO
            global_here%IGEP=global_here%IGEP+4
            DO I=1,3
              WRITE(s%fort63unit,REC=global_here%IGEP+I) global_here%RID8(I)
              ENDDO
            global_here%IGEP=global_here%IGEP+3
            DO I=1,3
              WRITE(s%fort63unit,REC=global_here%IGEP+I) global_here%AID8(I)
              ENDDO
            global_here%IGEP=global_here%IGEP+3
            ENDIF
          WRITE(s%fort63unit,REC=global_here%IGEP+1) global_here%NDSETSE
          WRITE(s%fort63unit,REC=global_here%IGEP+2) global_here%NE
          WRITE(s%fort63unit,REC=global_here%IGEP+3) global_here%DT*global_here%NSPOOLGE
          WRITE(s%fort63unit,REC=global_here%IGEP+4) global_here%NSPOOLGE
          WRITE(s%fort63unit,REC=global_here%IGEP+5) 1
          global_here%IGEP=global_here%IGEP+5
          CLOSE(s%fort63unit)                    ! DO THIS TO FLUSH THE WRITE BUFFER
          OPEN(s%fort63unit,FILE=S%DIRNAME//'/'//'fort.63',&
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
          OPEN(s%fort64unit,FILE=S%DIRNAME//'/'//'fort.64')
          WRITE(s%fort64unit,3220) global_here%RUNDES,global_here%RUNID,global_here%AGRID
          WRITE(s%fort64unit,3645) global_here%NDSETSV,global_here%NE,global_here%DTDP*global_here%NSPOOLGV,global_here%NSPOOLGV,2
          global_here%IGVP=2
          ENDIF

        IF(ABS(global_here%NOUTGV).EQ.2) THEN
          OPEN(s%fort64unit,FILE=S%DIRNAME//'/'//'fort.64',&
          ACCESS='DIRECT',RECL=NBYTE)
          IF(NBYTE.EQ.4) THEN
            DO I=1,8
              WRITE(s%fort64unit,REC=global_here%IGVP+I) global_here%RDES4(I)
              ENDDO
            global_here%IGVP=global_here%IGVP+8
            DO I=1,6
              WRITE(s%fort64unit,REC=global_here%IGVP+I) global_here%RID4(I)
              ENDDO
            global_here%IGVP=global_here%IGVP+6
            DO I=1,6
              WRITE(s%fort64unit,REC=global_here%IGVP+I) global_here%AID4(I)
              ENDDO
            global_here%IGVP=global_here%IGVP+6
            ENDIF
          IF(NBYTE.EQ.8) THEN
            DO I=1,4
              WRITE(s%fort64unit,REC=global_here%IGVP+I) global_here%RDES8(I)
              ENDDO
            global_here%IGVP=global_here%IGVP+4
            DO I=1,3
              WRITE(s%fort64unit,REC=global_here%IGVP+I) global_here%RID8(I)
              ENDDO
            global_here%IGVP=global_here%IGVP+3
            DO I=1,3
              WRITE(s%fort64unit,REC=global_here%IGVP+I) global_here%AID8(I)
              ENDDO
            global_here%IGVP=global_here%IGVP+3
            ENDIF
          WRITE(s%fort64unit,REC=global_here%IGVP+1) global_here%NDSETSV
          WRITE(s%fort64unit,REC=global_here%IGVP+2) global_here%NE
          WRITE(s%fort64unit,REC=global_here%IGVP+3) global_here%DT*global_here%NSPOOLGV
          WRITE(s%fort64unit,REC=global_here%IGVP+4) global_here%NSPOOLGV
          WRITE(s%fort64unit,REC=global_here%IGVP+5) 2
          global_here%IGVP=global_here%IGVP+5
          CLOSE(s%fort64unit)                    ! DO THIS TO FLUSH THE WRITE BUFFER
          OPEN(s%fort64unit,FILE=S%DIRNAME//'/'//'fort.64',&
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
          open(s%fort73unit,file=s%DIRNAME//'/'//'fort.73')
          write(s%fort73unit,3220) global_here%rundes,global_here%runid,global_here%agrid
          write(s%fort73unit,3645) global_here%ndsetsw,global_here%np,global_here%dtdp*global_here%nspoolgw,global_here%nspoolgw,1
          global_here%igpp=2
          OPEN(s%fort74unit,FILE=S%DIRNAME//'/'//'fort.74')
          WRITE(s%fort74unit,3220) global_here%RUNDES,global_here%RUNID,global_here%AGRID
          WRITE(s%fort74unit,3645) global_here%NDSETSW,global_here%NP,global_here%DTDP*global_here%NSPOOLGW,global_here%NSPOOLGW,2
          global_here%IGWP=2
          ENDIF

        IF(ABS(global_here%NOUTGW).EQ.2) THEN
          open(s%fort73unit,file=s%DIRNAME//'/'//'fort.73',&
          access='direct',recl=nbyte)
          OPEN(s%fort74unit,FILE=S%DIRNAME//'/'//'fort.74',&
          ACCESS='DIRECT',RECL=NBYTE)
          IF(NBYTE.EQ.4) THEN
            DO I=1,8
              write(s%fort73unit,rec=global_here%igpp+i) global_here%rdes4(i)
              WRITE(s%fort74unit,REC=global_here%IGWP+I) global_here%RDES4(I)
              ENDDO
            global_here%igpp=global_here%igpp+8
            global_here%IGWP=global_here%IGWP+8
            DO I=1,6
              write(s%fort73unit,rec=global_here%igpp+i) global_here%rid4(i)
              WRITE(s%fort74unit,REC=global_here%IGWP+I) global_here%RID4(I)
              ENDDO
            global_here%igpp=global_here%igpp+6
            global_here%IGWP=global_here%IGWP+6
            DO I=1,6
              write(s%fort73unit,rec=global_here%igpp+i) global_here%aid4(i)
              WRITE(s%fort74unit,REC=global_here%IGWP+I) global_here%AID4(I)
              ENDDO
            global_here%igpp=global_here%igpp+6
            global_here%IGWP=global_here%IGWP+6
            ENDIF
          IF(NBYTE.EQ.8) THEN
            DO I=1,4
              write(s%fort73unit,rec=global_here%igpp+i) global_here%rdes8(i)
              WRITE(s%fort74unit,REC=global_here%IGWP+I) global_here%RDES8(I)
              ENDDO
            global_here%igpp=global_here%igpp+4
            global_here%IGWP=global_here%IGWP+4
            DO I=1,3
              write(s%fort73unit,rec=global_here%igpp+i) global_here%rid8(i)
              WRITE(s%fort74unit,REC=global_here%IGWP+I) global_here%RID8(I)
              ENDDO
            global_here%igpp=global_here%igpp+3
            global_here%IGWP=global_here%IGWP+3
            DO I=1,3
              write(s%fort73unit,rec=global_here%igpp+i) global_here%aid8(i)
              WRITE(s%fort74unit,REC=global_here%IGWP+I) global_here%AID8(I)
              ENDDO
            global_here%igpp=global_here%igpp+3
            global_here%IGWP=global_here%IGWP+3
            ENDIF
          write(s%fort73unit,rec=global_here%igpp+1) global_here%ndsetsw
          write(s%fort73unit,rec=global_here%igpp+2) global_here%np
          write(s%fort73unit,rec=global_here%igpp+3) global_here%dt*global_here%nspoolgw
          write(s%fort73unit,rec=global_here%igpp+4) global_here%nspoolgw
          write(s%fort73unit,rec=global_here%igpp+5) 2
          global_here%igpp=global_here%igpp+5
          close(s%fort73unit)                    ! DO THIS TO FLUSH THE WRITE BUFFER
          open(s%fort73unit,file=s%DIRNAME//'/'//'fort.73',&
         access='direct',recl=nbyte)
          WRITE(s%fort74unit,REC=global_here%IGWP+1) global_here%NDSETSW
          WRITE(s%fort74unit,REC=global_here%IGWP+2) global_here%NP
          WRITE(s%fort74unit,REC=global_here%IGWP+3) global_here%DT*global_here%NSPOOLGW
          WRITE(s%fort74unit,REC=global_here%IGWP+4) global_here%NSPOOLGW
          WRITE(s%fort74unit,REC=global_here%IGWP+5) 2
          global_here%IGWP=global_here%IGWP+5
          CLOSE(s%fort74unit)                    ! DO THIS TO FLUSH THE WRITE BUFFER
          OPEN(s%fort74unit,FILE=S%DIRNAME//'/'//'fort.74',&
         ACCESS='DIRECT',RECL=NBYTE)
          ENDIF
#ifdef SWAN
!asey 101118: Added the output of radiation stress gradients.
       IGRadS=0
       IF(ABS(global_here%NOUTGW).EQ.1) THEN
          OPEN(rads64unit,FILE=S%DIRNAME//'/'//'rads.64')
          WRITE(rads64unit,3220) global_here%RUNDES,global_here%RUNID,global_here%AGRID
          WRITE(rads64unit,3645) global_here%NDSETSW,global_here%NP,global_here%DTDP*global_here%NSPOOLGW,global_here%NSPOOLGW,2
          IGRadS=2
          ENDIF
       IF(ABS(global_here%NOUTGW).EQ.2) THEN
          OPEN(rads64unit,FILE=TRIM(LOCALDIR)//'/'//'rads.64',&
           ACCESS='DIRECT',RECL=NByte)
          IF(NBYTE.EQ.4) THEN
             DO I=1,8
                WRITE(rads64unit,REC=IGRadS+I) global_here%RDES4(I)
                ENDDO
             IGRadS=IGRadS+8
             DO I=1,6
                WRITE(rads64unit,REC=IGRadS+I) global_here%RID4(I)
                ENDDO
             IGRadS=IGRadS+6
             DO I=1,6
                WRITE(rads64unit,REC=IGRadS+I) global_here%AID4(I)
                ENDDO
             IGRadS=IGRadS+6
             ENDIF
          IF(NBYTE.EQ.8) THEN
             DO I=1,4
                WRITE(rads64unit,REC=IGRadS+I) global_here%RDES8(I)
                ENDDO
             IGRadS=IGRadS+4
             DO I=1,3
                WRITE(rads64unit,REC=IGRadS+I) global_here%RID8(I)
                ENDDO
             IGRadS=IGRadS+3
             DO I=1,3
                WRITE(rads64unit,REC=IGRadS+I) global_here%AID8(I)
                ENDDO
             IGRadS=IGRadS+3
             ENDIF
          WRITE(rads64unit,REC=IGRadS+1) global_here%NDSETSW
          WRITE(rads64unit,REC=IGRadS+2) global_here%NP
          WRITE(rads64unit,REC=IGRadS+3) global_here%DT*global_here%NSPOOLGW
          WRITE(rads64unit,REC=IGRadS+4) global_here%NSPOOLGW
          WRITE(rads64unit,REC=IGRadS+5) 2
          IGRadS=IGRadS+5
          CLOSE(rads64unit)
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
          OPEN(s%fort83unit,FILE=S%DIRNAME//'/'//'fort.83')
          WRITE(s%fort83unit,3220) global_here%RUNDES,global_here%RUNID,global_here%AGRID
          WRITE(s%fort83unit,3645) global_here%NDSETSC,global_here%NP,global_here%DTDP*global_here%NSPOOLGC,global_here%NSPOOLGC,1
          global_here%IGCP=2
          ENDIF

        IF(ABS(global_here%NOUTGC).EQ.2) THEN
          OPEN(s%fort83unit,FILE=S%DIRNAME//'/'//'fort.83',&
          ACCESS='DIRECT',RECL=NBYTE)
          IF(NBYTE.EQ.4) THEN
            DO I=1,8
              WRITE(s%fort83unit,REC=global_here%IGCP+I) global_here%RDES4(I)
              ENDDO
            global_here%IGCP=global_here%IGCP+8
            DO I=1,6
              WRITE(s%fort83unit,REC=global_here%IGCP+I) global_here%RID4(I)
              ENDDO
            global_here%IGCP=global_here%IGCP+6
            DO I=1,6
              WRITE(s%fort83unit,REC=global_here%IGCP+I) global_here%AID4(I)
              ENDDO
            global_here%IGCP=global_here%IGCP+6
            ENDIF
          IF(NBYTE.EQ.8) THEN
            DO I=1,4
              WRITE(s%fort83unit,REC=global_here%IGCP+I) global_here%RDES8(I)
              ENDDO
            global_here%IGCP=global_here%IGCP+4
            DO I=1,3
              WRITE(s%fort83unit,REC=global_here%IGCP+I) global_here%RID8(I)
              ENDDO
            global_here%IGCP=global_here%IGCP+3
            DO I=1,3
              WRITE(s%fort83unit,REC=global_here%IGCP+I) global_here%AID8(I)
              ENDDO
            global_here%IGCP=global_here%IGCP+3
            ENDIF
          WRITE(s%fort83unit,REC=global_here%IGCP+1) global_here%NDSETSC
          WRITE(s%fort83unit,REC=global_here%IGCP+2) global_here%NP
          WRITE(s%fort83unit,REC=global_here%IGCP+3) global_here%DT*global_here%NSPOOLGC
          WRITE(s%fort83unit,REC=global_here%IGCP+4) global_here%NSPOOLGC
          WRITE(s%fort83unit,REC=global_here%IGCP+5) 1
          global_here%IGCP=global_here%IGCP+5
          CLOSE(s%fort83unit)                    ! DO THIS TO FLUSH THE WRITE BUFFER
          OPEN(s%fort83unit,FILE=S%DIRNAME//'/'//'fort.83',&
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
            OPEN(s%fort84unit,FILE=S%DIRNAME//'/'//'fort.84')
            OPEN(s%fort85unit,FILE=S%DIRNAME//'/'//'fort.85')
            WRITE(s%fort84unit,3220) global_here%RUNDES,global_here%RUNID,global_here%AGRID
            WRITE(s%fort84unit,3645) global_here%NDSETSE,global_here%NP,global_here%DTDP*global_here%NSPOOLGE,global_here%NSPOOLGE,1
            WRITE(s%fort85unit,3220) global_here%RUNDES,global_here%RUNID,global_here%AGRID
            WRITE(s%fort85unit,3645) global_here%NDSETSE,global_here%NP,global_here%DTDP*global_here%NSPOOLGE,global_here%NSPOOLGE,1
            global_here%IGEP=2
          ENDIF

          IF (ABS(global_here%NOUTGE).EQ.2) THEN
            OPEN(s%fort84unit,FILE=S%DIRNAME//'/'//'fort.84',&
          ACCESS='DIRECT',RECL=NBYTE)
            OPEN(s%fort85unit,FILE=S%DIRNAME//'/'//'fort.85',&
          ACCESS='DIRECT',RECL=NBYTE)
            IF(NBYTE.EQ.4) THEN
              DO I=1,8
                WRITE(s%fort84unit,REC=global_here%IGEP+I) global_here%RDES4(I)
                WRITE(s%fort85unit,REC=global_here%IGEP+I) global_here%RDES4(I)
              ENDDO
              global_here%IGEP=global_here%IGEP+8
              DO I=1,6
                WRITE(s%fort84unit,REC=global_here%IGEP+I) global_here%RID4(I)
                WRITE(s%fort85unit,REC=global_here%IGEP+I) global_here%RID4(I)
              ENDDO
              global_here%IGEP=global_here%IGEP+6
              DO I=1,6
                WRITE(s%fort84unit,REC=global_here%IGEP+I) global_here%AID4(I)
                WRITE(s%fort85unit,REC=global_here%IGEP+I) global_here%AID4(I)
              ENDDO
              global_here%IGEP=global_here%IGEP+6
              ENDIF
            IF(NBYTE.EQ.8) THEN
              DO I=1,4
                WRITE(s%fort84unit,REC=global_here%IGEP+I) global_here%RDES8(I)
                WRITE(s%fort85unit,REC=global_here%IGEP+I) global_here%RDES8(I)
              ENDDO
              global_here%IGEP=global_here%IGEP+4
              DO I=1,3
                WRITE(s%fort84unit,REC=global_here%IGEP+I) global_here%RID8(I)
                WRITE(s%fort85unit,REC=global_here%IGEP+I) global_here%RID8(I)
              ENDDO
              global_here%IGEP=global_here%IGEP+3
              DO I=1,3
                WRITE(s%fort84unit,REC=global_here%IGEP+I) global_here%AID8(I)
                WRITE(s%fort85unit,REC=global_here%IGEP+I) global_here%AID8(I)
              ENDDO
              global_here%IGEP=global_here%IGEP+3
            ENDIF
            WRITE(s%fort84unit,REC=global_here%IGEP+1) global_here%NDSETSE
            WRITE(s%fort84unit,REC=global_here%IGEP+2) global_here%NP
            WRITE(s%fort84unit,REC=global_here%IGEP+3) global_here%DT*global_here%NSPOOLGE
            WRITE(s%fort84unit,REC=global_here%IGEP+4) global_here%NSPOOLGE
            WRITE(s%fort84unit,REC=global_here%IGEP+5) 1
            WRITE(s%fort85unit,REC=global_here%IGEP+1) global_here%NDSETSE
            WRITE(s%fort85unit,REC=global_here%IGEP+2) global_here%NP
            WRITE(s%fort85unit,REC=global_here%IGEP+3) global_here%DT*global_here%NSPOOLGE
            WRITE(s%fort85unit,REC=global_here%IGEP+4) global_here%NSPOOLGE
            WRITE(s%fort85unit,REC=global_here%IGEP+5) 1
            global_here%IGEP=global_here%IGEP+5
            CLOSE(s%fort84unit)                    ! DO THIS TO FLUSH THE WRITE BUFFER
            CLOSE(s%fort85unit)                    ! DO THIS TO FLUSH THE WRITE BUFFER
            OPEN(s%fort84unit,FILE=S%DIRNAME//'/'//'fort.84',&
         ACCESS='DIRECT',RECL=NBYTE)
            OPEN(s%fort85unit,FILE=S%DIRNAME//'/'//'fort.85',&
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
            OPEN(s%fort82unit,FILE=S%DIRNAME//'/'//'fort.82')
            WRITE(s%fort82unit,3220) global_here%RUNDES,global_here%RUNID,global_here%AGRID
            WRITE(s%fort82unit,3645) global_here%NTRSPE,global_here%NSTAE,global_here%DTDP*global_here%NSPOOLE,global_here%NSPOOLE,1
            global_here%IESTP=2
          ENDIF

          IF (ABS(global_here%NOUTE).EQ.2) THEN
            OPEN(s%fort82unit,FILE=S%DIRNAME//'/'//'fort.82',&
          ACCESS='DIRECT',RECL=NBYTE)
            IF(NBYTE.EQ.4) THEN
              DO I=1,8
                WRITE(s%fort82unit,REC=global_here%IESTP+I) global_here%RDES4(I)
              ENDDO
              global_here%IESTP=global_here%IESTP+8
              DO I=1,6
                WRITE(s%fort82unit,REC=global_here%IESTP+I) global_here%RID4(I)
              ENDDO
              global_here%IESTP=global_here%IESTP+6
              DO I=1,6
                WRITE(s%fort82unit,REC=global_here%IESTP+I) global_here%AID4(I)
              ENDDO
              global_here%IESTP=global_here%IESTP+6
            ENDIF
            IF(NBYTE.EQ.8) THEN
              DO I=1,4
                WRITE(s%fort82unit,REC=global_here%IESTP+I) global_here%RDES8(I)
              ENDDO
              global_here%IESTP=global_here%IESTP+4
              DO I=1,3
                WRITE(s%fort82unit,REC=global_here%IESTP+I) global_here%RID8(I)
              ENDDO
              global_here%IESTP=global_here%IESTP+3
              DO I=1,3
                WRITE(s%fort82unit,REC=global_here%IESTP+I) global_here%AID8(I)
              ENDDO
              global_here%IESTP=global_here%IESTP+3
            ENDIF
            WRITE(s%fort82unit,REC=global_here%IESTP+1) global_here%NTRSPE
            WRITE(s%fort82unit,REC=global_here%IESTP+2) global_here%NSTAE
            WRITE(s%fort82unit,REC=global_here%IESTP+3) global_here%DT*global_here%NSPOOLE
            WRITE(s%fort82unit,REC=global_here%IESTP+4) global_here%NSPOOLE
            WRITE(s%fort82unit,REC=global_here%IESTP+5) 1
            global_here%IESTP=global_here%IESTP+5
            CLOSE(s%fort82unit)                    ! DO THIS TO FLUSH THE WRITE BUFFER
            OPEN(s%fort82unit,FILE=S%DIRNAME//'/'//'fort.82',&
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
            OPEN(s%fort94unit,FILE=S%DIRNAME//'/'//'fort.94')
            WRITE(s%fort94unit,3220) global_here%RUNDES,global_here%RUNID,global_here%AGRID
            WRITE(s%fort94unit,3645) global_here%NDSETSV,global_here%NP,global_here%DTDP*global_here%NSPOOLGV,global_here%NSPOOLGV,2
            global_here%IGVP=2
          ENDIF

          IF(ABS(global_here%NOUTGV).EQ.2) THEN
            OPEN(s%fort94unit,FILE=S%DIRNAME//'/'//'fort.94',&
          ACCESS='DIRECT',RECL=NBYTE)
            IF(NBYTE.EQ.4) THEN
              DO I=1,8
                WRITE(s%fort94unit,REC=global_here%IGVP+I) global_here%RDES4(I)
              ENDDO
              global_here%IGVP=global_here%IGVP+8
              DO I=1,6
                WRITE(s%fort94unit,REC=global_here%IGVP+I) global_here%RID4(I)
              ENDDO
              global_here%IGVP=global_here%IGVP+6
              DO I=1,6
                WRITE(s%fort94unit,REC=global_here%IGVP+I) global_here%AID4(I)
              ENDDO
              global_here%IGVP=global_here%IGVP+6
            ENDIF
            IF(NBYTE.EQ.8) THEN
              DO I=1,4
                WRITE(s%fort94unit,REC=global_here%IGVP+I) global_here%RDES8(I)
              ENDDO
              global_here%IGVP=global_here%IGVP+4
              DO I=1,3
                WRITE(s%fort94unit,REC=global_here%IGVP+I) global_here%RID8(I)
              ENDDO
              global_here%IGVP=global_here%IGVP+3
              DO I=1,3
                WRITE(s%fort94unit,REC=global_here%IGVP+I) global_here%AID8(I)
              ENDDO
              global_here%IGVP=global_here%IGVP+3
            ENDIF
            WRITE(s%fort94unit,REC=global_here%IGVP+1) global_here%NDSETSV
            WRITE(s%fort94unit,REC=global_here%IGVP+2) global_here%NP
            WRITE(s%fort94unit,REC=global_here%IGVP+3) global_here%DT*global_here%NSPOOLGV
            WRITE(s%fort94unit,REC=global_here%IGVP+4) global_here%NSPOOLGV
            WRITE(s%fort94unit,REC=global_here%IGVP+5) 2
            global_here%IGVP=global_here%IGVP+5
            CLOSE(s%fort94unit)                    ! DO THIS TO FLUSH THE WRITE BUFFER
            OPEN(s%fort94unit,FILE=S%DIRNAME//'/'//'fort.94',&
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
            OPEN(s%fort96unit,FILE=S%DIRNAME//'/'//'fort.96')
            WRITE(s%fort96unit,3220) global_here%RUNDES,global_here%RUNID,global_here%AGRID
            WRITE(s%fort96unit,3645) global_here%NDSETSV,global_here%NP,global_here%DTDP*global_here%NSPOOLGV,global_here%NSPOOLGV,2
            global_here%IGVP=2
          ENDIF

          IF(ABS(global_here%NOUTGV).EQ.2) THEN
            OPEN(s%fort96unit,FILE=S%DIRNAME//'/'//'fort.96',&
          ACCESS='DIRECT',RECL=NBYTE)
            IF(NBYTE.EQ.4) THEN
              DO I=1,8
                WRITE(s%fort96unit,REC=global_here%IGVP+I) global_here%RDES4(I)
              ENDDO
              global_here%IGVP=global_here%IGVP+8
              DO I=1,6
                WRITE(s%fort96unit,REC=global_here%IGVP+I) global_here%RID4(I)
              ENDDO
              global_here%IGVP=global_here%IGVP+6
              DO I=1,6
                WRITE(s%fort96unit,REC=global_here%IGVP+I) global_here%AID4(I)
              ENDDO
              global_here%IGVP=global_here%IGVP+6
            ENDIF
            IF(NBYTE.EQ.8) THEN
              DO I=1,4
                WRITE(s%fort96unit,REC=global_here%IGVP+I) global_here%RDES8(I)
              ENDDO
              global_here%IGVP=global_here%IGVP+4
              DO I=1,3
                WRITE(s%fort96unit,REC=global_here%IGVP+I) global_here%RID8(I)
              ENDDO
              global_here%IGVP=global_here%IGVP+3
              DO I=1,3
                WRITE(s%fort96unit,REC=global_here%IGVP+I) global_here%AID8(I)
              ENDDO
              global_here%IGVP=global_here%IGVP+3
            ENDIF
            WRITE(s%fort96unit,REC=global_here%IGVP+1) global_here%NDSETSV
            WRITE(s%fort96unit,REC=global_here%IGVP+2) global_here%NP
            WRITE(s%fort96unit,REC=global_here%IGVP+3) global_here%DT*global_here%NSPOOLGV
            WRITE(s%fort96unit,REC=global_here%IGVP+4) global_here%NSPOOLGV
            WRITE(s%fort96unit,REC=global_here%IGVP+5) 2
            global_here%IGVP=global_here%IGVP+5
            CLOSE(s%fort96unit)                    ! DO THIS TO FLUSH THE WRITE BUFFER
            OPEN(s%fort96unit,FILE=S%DIRNAME//'/'//'fort.96',&
         ACCESS='DIRECT',RECL=NBYTE)
          ENDIF

        ENDIF
          
      RETURN
      END
