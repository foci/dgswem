C**************************************************************************
C PADCIRC RELEASE VERSION 41.11 09/14/2001  
C    last changes in this file VERSION 41.11
C
C  mod history
C  v40.02mxxx - date - programmer - describe change 
C                    - mark change in code with  cinitials-mxxx 
C
C  v41.11 - 09/14/01 - rl - from 41.10 - added NWS = -2 capability
C  v41.09 - 06/30/01 - jw - from 41.08 - made minor mods as per vp version 41.05
C  v40.02m002 - 12/22 - jjw/vjp - Vic suggested this change to avoid compiler conflicts
C  v40.02m001 - 12/21 - jjw - add cross barrier pipes cjjwm001
C************************************************************************** 
C 
        SUBROUTINE HOTSTART()
C
C**************************************************************************
C
C  HOT START PROGRAM SETUP ROUTINE 
C
C**************************************************************************
C
      USE GLOBAL
      USE HARM
      USE WIND
      !USE SED
      USE OWIWIND,ONLY : NWS12INIT,NWS12GET
      USE NodalAttributes, ONLY :
     &     LoadDirEffRLen,
     &     ApplyDirectionalWindReduction
     
      IMPLICIT NONE
      INTEGER i,j, IT
      REAL(SZ) WindDragLimit    !jgf46.01 Cap on max wind drag coefficient
      REAL(SZ) RampExtFlux1   ! Ramp for external flux b.c.s @ITHS-1
      REAL(SZ) RampExtFlux2   ! Ramp for external flux b.c.s @ITHS
      REAL(SZ) RampIntFlux1   ! Ramp for internal flux b.c.s @ITHS-1
      REAL(SZ) RampIntFlux2   ! Ramp for internal flux b.c.s @ITHS
      REAL(SZ) RampTip2    ! Ramp for tidal potential @ITHS
      REAL(SZ) RampMete2   ! Ramp for wind and atmospheric pressure @ITHS
      REAL(SZ) RampWRad2   ! Ramp for wave radiation stress @ITHS
C...
C......READ IN 2DDI HOT START INITIAL CONDITION OVER WHOLE DOMAIN
C......THIS FILE ALWAYS HAS A RECL=8 BECAUSE IT IS ASSUMED THAT THE HARMONIC
C......ANALYSIS IS ALWAYS DONE IN 64 BITS, EVEN ON A WORKSTATION
C...
        IF(IHOT.EQ.67) OPEN(67,FILE=DIRNAME//'/'//'fort.67',
     &        ACCESS='DIRECT',RECL=8)
        IF(IHOT.EQ.68) OPEN(68,FILE=DIRNAME//'/'//'fort.68',
     &       ACCESS='DIRECT',RECL=8)
        IHOTSTP=1
        READ(IHOT,REC=IHOTSTP) IMHS
        IHOTSTP=2
        READ(IHOT,REC=IHOTSTP) TIME_A
        IHOTSTP=3
        READ(IHOT,REC=IHOTSTP) ITHS
        DO I=1,NP
          READ(IHOT,REC=IHOTSTP+1) ETA1(I)
          READ(IHOT,REC=IHOTSTP+2) ETA2(I)
          READ(IHOT,REC=IHOTSTP+3) UU2(I)
          READ(IHOT,REC=IHOTSTP+4) VV2(I)
Casey 130211: Debug.
         !IF (SEDFLAG.GE.1) THEN
         !  READ(IHOT,REC=IHOTSTP+5) DP(I)
         !  READ(IHOT,REC=IHOTSTP+6) DP0(I)
         !  IHOTSTP = IHOTSTP + 6
         !ELSE
            IHOTSTP=IHOTSTP+4
         !ENDIF
          
          IF(IMHS.EQ.10) THEN
            READ(IHOT,REC=IHOTSTP+1) CH1(I)
            IHOTSTP=IHOTSTP+1
            ENDIF
          READ(IHOT,REC=IHOTSTP+1) NNODECODE(I)
          IHOTSTP=IHOTSTP+1
          ETAS(I)=ETA2(I)-ETA1(I)
          NODEREP(I)=MAX0(NODEWETMIN,NODEDRYMIN)
          END DO

        RAMP2=1.0D0
        RAMP1=1.0D0
        IF(NRAMP.EQ.1) THEN
          RAMP1=TANH((2.D0*(ITHS-1)*DTDP/86400.D0)/DRAMP)
          RAMP2=TANH((2.D0*ITHS*DTDP/86400.D0)/DRAMP)
        ENDIF
        
        IF (NRamp.eq.0) THEN
          RampExtFlux1 = 1.0d0
          RampExtFlux2 = 1.0d0
          RampIntFlux1 = 1.0d0
          RampIntFlux2 = 1.0d0
          RampTip2     = 1.0d0
          RampMete2    = 1.0d0
          RampWRad2    = 1.0d0
        ELSE
          RampExtFlux1= TANH((2.D0*(ITHS-1)*DTDP/86400.D0)/DRampExtFlux)
          RampExtFlux2= TANH((2.D0*(ITHS)*DTDP/86400.D0)/DRampExtFlux)
          RampIntFlux1= TANH((2.D0*(ITHS-1)*DTDP/86400.D0)/DRampIntFlux)
          RampIntFlux2= TANH((2.D0*(ITHS)*DTDP/86400.D0)/DRampIntFlux)
          RampTip2    = TANH((2.D0*(ITHS)*DTDP/86400.D0)/DRampTip)
          RampMete2   = TANH((2.D0*(ITHS)*DTDP/86400.D0)/DRampMete)
          RampWRad2   = TANH((2.D0*(ITHS)*DTDP/86400.D0)/DRampWRad)
        ENDIF

C
C....SET POSITIONS IN BOUNDARY CONDITION, WIND AND OUTPUT FILES
C
        WRITE(16,1112)
        WRITE(16,1794)
1794    FORMAT(//,' INFORMATION ABOUT RE-STARTING THE TIME SERIES',
     &            ' OUTPUT FILES (UNITS 61-64,71,74),',
     &    /,' WIND/PRESSURE FILE (UNIT 22) AND FLOW BOUNDARY CONDITION',
     &      ' FILE (UNIT 20)',//)

C......INITIALLY, ZERO OUT THE NORMAL FLOW ON ALL BOUNDARIES

        DO I=1,NVEL
          QN2(I)=0.D0
          QN1(I)=0.D0
          QN0(I)=0.D0
          END DO

C....FIND THE PROPER PLACE IN THE APERIODIC ELEVATION SPECIFIED BOUNDARY CONDITION
C....FILE IF IT IS REQURIED.

        IF((NOPE.GT.0).AND.(NBFR.EQ.0)) THEN
          IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,1112)
          WRITE(16,1112)
          IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,1976)
          WRITE(16,1976)
 1976     FORMAT(/,1X,'LOCATING ELEVATION SPECIFIED INFORMATION IN ',
     &                'UNIT 19',/)
          OPEN(19,FILE=DIRNAME//'/'//'fort.19')
          READ(19,*) ETIMINC
          ETIME1=STATIM*86400.D0
          ETIME2=ETIME1+ETIMINC
          DO J=1,NETA
            READ(19,*) ESBIN1(J)
                     END DO
          DO J=1,NETA
              READ(19,*) ESBIN2(J)
            END DO
          DO IT=1,ITHS-1
            TIMEIT=IT*DTDP + STATIM*86400.D0
            IF(TIMEIT.GT.ETIME2) THEN
              ETIME1=ETIME2
              ETIME2=ETIME1+ETIMINC
              DO J=1,NETA
                ESBIN1(J)=ESBIN2(J)        
                READ(19,*) ESBIN2(J)
                END DO
              ENDIF
            END DO
          IF(TIME_A.GT.ETIME2) THEN
            ETIME1=ETIME2
            ETIME2=ETIME1+ETIMINC
            DO J=1,NETA
              ESBIN1(J)=ESBIN2(J)
              READ(19,*) ESBIN2(J)
              END DO
            ENDIF
          ETRATIO=(TIMEIT-ETIME1)/ETIMINC
          ENDIF

C......FIND PROPER PLACE IN THE APERIODIC NORMAL FLOW BOUNDARY CONDITION FILE IF IT
C......IS REQUIRED

        IF((NFLUXF.EQ.1).AND.(NFFR.EQ.0)) THEN
          IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,1112)
          WRITE(16,1112)
          IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,1978)
          WRITE(16,1978)
 1978     FORMAT(/,1X,'LOCATING NORMAL FLOW INFORMATION IN UNIT 20',/)
          OPEN(20,FILE=DIRNAME//'/'//'fort.20')
            READ(20,*) FTIMINC
          QTIME1=STATIM*86400.D0
          QTIME2=QTIME1+FTIMINC
          DO J=1,NVEL
            QNIN1(J)=0.D0
            IF((LBCODEI(J).EQ.2).OR.(LBCODEI(J).EQ.12)
     &                          .OR.(LBCODEI(J).EQ.22))
     &                                              READ(20,*) QNIN1(J)
            END DO
          DO J=1,NVEL
            QNIN2(J)=0.D0
            IF((LBCODEI(J).EQ.2).OR.(LBCODEI(J).EQ.12)
     &                          .OR.(LBCODEI(J).EQ.22))
     &                                              READ(20,*) QNIN2(J)
            END DO
          DO IT=1,ITHS-1
            TIMEIT=IT*DTDP + STATIM*86400.D0
            IF(TIMEIT.GT.QTIME2) THEN
              QTIME1=QTIME2
              QTIME2=QTIME2+FTIMINC
              DO J=1,NVEL
                IF((LBCODEI(J).EQ.2).OR.(LBCODEI(J).EQ.12)
     &                              .OR.(LBCODEI(J).EQ.22)) THEN
                  QNIN1(J)=QNIN2(J)        
                  READ(20,*) QNIN2(J)
                  ENDIF
                END DO
              ENDIF
            END DO
          QTRATIO=(TIMEIT-QTIME1)/FTIMINC
          DO I=1,NVEL
            QN1(I)=RAMP1*(QNIN1(I)+QTRATIO*(QNIN2(I)-QNIN1(I)))
            END DO
          IF(TIME_A.GT.QTIME2) THEN
            QTIME1=QTIME2
            QTIME2=QTIME1+FTIMINC
            DO J=1,NVEL
              IF((LBCODEI(J).EQ.2).OR.(LBCODEI(J).EQ.12)
     &                            .OR.(LBCODEI(J).EQ.22)) THEN
                QNIN1(J)=QNIN2(J)        
                READ(20,*) QNIN2(J)
                ENDIF
              END DO
            ENDIF
          QTRATIO=(TIME_A-QTIME1)/FTIMINC
          DO I=1,NVEL
            QN2(I)=RAMP2*(QNIN1(I)+QTRATIO*(QNIN2(I)-QNIN1(I)))
            END DO
          ENDIF

C......RESTART THE PERIODIC NORMAL FLOW BOUNDARY CONDITION

        IF((NFLUXF.EQ.1).AND.(NFFR.GT.0)) THEN
          DO J=1,NFFR
            IF(FPER(J).EQ.0.) THEN
              NCYC=0.
              ELSE
              NCYC=INT(TIME_A/FPER(J))
              ENDIF
            ARGJ1=FAMIG(J)*(TIME_A-DTDP-NCYC*FPER(J))+FFACE(J)
            ARGJ2=FAMIG(J)*(TIME_A-NCYC*FPER(J))+FFACE(J)
            RFF1=FFF(J)*RAMP1
            RFF2=FFF(J)*RAMP2
            DO I=1,NVEL
              ARG1=ARGJ1-QNPH(J,I)
              ARG2=ARGJ2-QNPH(J,I)
              QN1(I)=QN1(I)+QNAM(J,I)*RFF1*COS(ARG1)
              QN2(I)=QN2(I)+QNAM(J,I)*RFF2*COS(ARG2)
              END DO
            END DO
            ENDIF

C...
C...RESTART SUPERCRITICAL OUTWARD NORMAL FLOW OVER SPECIFIED
C....EXTERNAL BARRIER BOUNDARY NODES
C...
        IF(NFLUXB.EQ.1) THEN
          DO I=1,NVEL
            IF((LBCODEI(I).EQ.3).OR.(LBCODEI(I).EQ.13)
     &        .OR.(LBCODEI(I).EQ.23)) THEN
              NNBB=NBV(I)
              RBARWL=2.D0*(ETA1(NNBB)-BARLANHT(I))/3.D0
              IF(RBARWL.GT.0.0D0) THEN
                QN1(I)=-RAMP1*BARLANCFSP(I)*RBARWL*(RBARWL*G)**0.5D0
              ENDIF
              RBARWL=2.D0*(ETA2(NNBB)-BARLANHT(I))/3.D0
              IF(RBARWL.GT.0.0D0) THEN
                QN2(I)=-RAMP2*BARLANCFSP(I)*RBARWL*(RBARWL*G)**0.5D0
              ENDIF
            ENDIF
          END DO
        ENDIF

C...
C...RESTART INWARD/OUTWARD NORMAL FLOW OVER SPECIFIED
cjjwm001 - modified/added the following 3 lines
C....INTERNAL BARRIERS AND FOR INTERNAL BARRIER BOUNDARIES
C....WITH CROSS BARRIER PIPES
C....THIS SECTION ONLY RESTARTS THE OVER BARRIER FLOW COMPONENT
C...
cjjwm001 - modified following line
        IF(NFLUXIB.EQ.1) THEN
          DO I=1,NVEL
cjjwm001 - modified following 2 lines
            IF((LBCODEI(I).EQ.4).OR.(LBCODEI(I).EQ.24) 
     &        .OR.(LBCODEI(I).EQ.5).OR.(LBCODEI(I).EQ.25)) THEN
              NNBB1=NBV(I)      ! GLOBAL NODE NUMBER ON THIS SIDE OF BARRIER
              NNBB2=IBCONN(I)   ! GLOBAL NODE NUMBER ON OPPOSITE SIDE OF BARRIER
C.............RESET INFORMATION FOR K-1 TIME LEVEL
              RBARWL1=ETA1(NNBB1)-BARINHT(I)
              RBARWL2=ETA1(NNBB2)-BARINHT(I)
              RBARWL1F=2.0D0*RBARWL1/3.0D0
              RBARWL2F=2.0D0*RBARWL2/3.0D0
              IF((RBARWL1.LT.0.0).AND.(RBARWL2.LT.0.0)) THEN ! WATER LEVEL BELOW BARRIER
                QN1(I)=0.0D0                                   ! NO FLOW
                GOTO 1998
              ENDIF
              IF(RBARWL1.EQ.RBARWL2) THEN ! WATER LEVEL EQUAL ON BOTH SIDES OF BARRIER
                QN1(I)=0.0D0                ! NO FLOW
                GOTO 1998
              ENDIF
              IF(RBARWL1.GT.RBARWL2) THEN ! WATER LEVEL GREATER ON THIS SIDE OF BARRIER
                IF(RBARWL2.GT.RBARWL1F) THEN ! OUTWARD SUBCRITICAL FLOW
                  QN1(I)=-RAMP1*BARINCFSB(I)*RBARWL2*
     &                    (2.d0*G*(RBARWL1-RBARWL2))**0.5D0
                  GOTO 1998
                ELSE                        ! OUTWARD SUPERCRITICAL FLOW
                  QN1(I)=-RAMP1*BARINCFSP(I)*RBARWL1F*
     &                    (RBARWL1F*G)**0.5D0
                  GOTO 1998
                ENDIF
              ENDIF
              IF(RBARWL2.GT.RBARWL1) THEN  ! WATER LEVEL LOWER ON THIS SIDE OF BARRIER
                IF(RBARWL1.GT.RBARWL2F) THEN ! INWARD SUBCRITICAL FLOW
                  QN1(I)=RAMP1*BARINCFSB(I)*RBARWL1*
     &                   (2.d0*G*(RBARWL2-RBARWL1))**0.5D0
                  GOTO 1998
                ELSE                         ! INWARD SUPERCRITICAL FLOW
                  QN1(I)=RAMP1*BARINCFSP(I)*RBARWL2F*(RBARWL2F*G)**0.5D0
                  GOTO 1998
                ENDIF
              ENDIF
1998          CONTINUE
C.............RESET INFORMATION FOR K TIME LEVEL
              RBARWL1=ETA2(NNBB1)-BARINHT(I)
              RBARWL2=ETA2(NNBB2)-BARINHT(I)
              RBARWL1F=2.0D0*RBARWL1/3.0D0
              RBARWL2F=2.0D0*RBARWL2/3.0D0
              IF((RBARWL1.LT.0.0).AND.(RBARWL2.LT.0.0)) THEN ! WATER LEVEL BELOW BARRIER
                QN2(I)=0.0D0                                   ! NO FLOW
                GOTO 1999
              ENDIF
              IF(RBARWL1.EQ.RBARWL2) THEN ! WATER LEVEL EQUAL ON BOTH SIDES OF BARRIER
                QN2(I)=0.0D0                ! NO FLOW
                GOTO 1999
              ENDIF
              IF(RBARWL1.GT.RBARWL2) THEN ! WATER LEVEL GREATER ON THIS SIDE OF BARRIER
                IF(RBARWL2.GT.RBARWL1F) THEN ! OUTWARD SUBCRITICAL FLOW
                  QN2(I)=-RAMP2*BARINCFSB(I)*RBARWL2*
     &                    (2.d0*G*(RBARWL1-RBARWL2))**0.5D0
                  GOTO 1999
                ELSE                         ! OUTWARD SUPERCRITICAL FLOW
                  QN2(I)=-RAMP2*BARINCFSP(I)*RBARWL1F*
     &                    (RBARWL1F*G)**0.5D0
                  GOTO 1999
                ENDIF
              ENDIF
              IF(RBARWL2.GT.RBARWL1) THEN !WATER LEVEL LOWER ON THIS SIDE OF BARRIER
                IF(RBARWL1.GT.RBARWL2F) THEN ! INWARD SUBCRITICAL FLOW
                  QN2(I)=RAMP2*BARINCFSB(I)*RBARWL1*
     &                   (2.d0*G*(RBARWL2-RBARWL1))**0.5D0
                  GOTO 1999
                ELSE                         ! INWARD SUPERCRITICAL FLOW
                  QN2(I)=RAMP2*BARINCFSP(I)*RBARWL2F*(RBARWL2F*G)**0.5D0
                  GOTO 1999
                ENDIF
              ENDIF
1999          CONTINUE
            ENDIF
          END DO
        ENDIF

cjjwm001 - start add
C...
C...RESTART INWARD/OUTWARD NORMAL FLOW OVER SPECIFIED
C....INTERNAL BARRIERS WITH CROSS BARRIER PIPES
C....THIS SECTION RESTARTS THE PIPE FLOW COMPONENT
C....NOTE THAT PIPE FLOW COMPONENT IS ADDED INTO BARRIER FLOW COMPONENT
C....THAT WAS PREVIOUSLY SET
C...
        IF(NFLUXIBP.EQ.1) THEN
          DO I=1,NVEL
            IF((LBCODEI(I).EQ.5).OR.(LBCODEI(I).EQ.25)) THEN
              NNBB1=NBV(I)      ! GLOBAL NODE NUMBER ON THIS SIDE OF BARRIER
              NNBB2=IBCONN(I)   ! GLOBAL NODE NUMBER ON OPPOSITE SIDE OF BARRIER
C.............RESET INFORMATION FOR K-1 TIME LEVEL
              RBARWL1=ETA1(NNBB1)-PIPEHT(I)
              RBARWL2=ETA1(NNBB2)-PIPEHT(I)
              IF((RBARWL1.LT.0.0).AND.(RBARWL2.LT.0.0)) THEN ! WATER LEVEL BELOW PIPE
                QN1(I)=QN1(I)+0.0D0                                 ! NO FLOW
                GOTO 2002
              ENDIF
              IF(RBARWL1.EQ.RBARWL2) THEN ! WATER LEVEL EQUAL ON BOTH SIDES OF PIPE
                QN1(I)=QN1(I)+0.0D0                ! NO FLOW
                GOTO 2002
              ENDIF
              IF(RBARWL1.GT.RBARWL2) THEN ! WATER LEVEL GREATER ON THIS SIDE OF PIPE
                IF(RBARWL2.LE.0) THEN ! OUTWARD FREE DISCHARGE 
                  QN1(I)=QN1(I)-RAMP1*0.25D0*PI*(PIPEDIAM(I))**2
     &                    *(2.D0*G*RBARWL1/(1+PIPECOEF(I)))**0.5D0
                  GOTO 2002
                ELSE                        ! OUTWARD SUBMERGED DISCHARGE
                  QN1(I)=QN1(I)-RAMP1*0.25D0*PI*(PIPEDIAM(I))**2
     &                    *(2.D0*G*(RBARWL1-RBARWL2)/PIPECOEF(I))**0.5D0
                  GOTO 2002
                ENDIF
              ENDIF
              IF(RBARWL2.GT.RBARWL1) THEN  ! WATER LEVEL LOWER ON THIS SIDE OF PIPE
                IF(RBARWL1.LE.0) THEN ! INWARD FREE DISCHARGE 
                  QN1(I)=QN1(I)+RAMP1*0.25D0*PI*(PIPEDIAM(I))**2
     &                    *(2.D0*G*RBARWL2/(1+PIPECOEF(I)))**0.5D0
                  GOTO 2002
                ELSE                         ! INWARD SUBMERGED DISCHARGE
                  QN1(I)=QN1(I)+RAMP1*0.25D0*PI*(PIPEDIAM(I))**2
     &                    *(2.D0*G*(RBARWL2-RBARWL1)/PIPECOEF(I))**0.5D0
                  GOTO 2002
                ENDIF
              ENDIF
2002          CONTINUE
C.............RESET INFORMATION FOR K TIME LEVEL
              RBARWL1=ETA2(NNBB1)-PIPEHT(I)
              RBARWL2=ETA2(NNBB2)-PIPEHT(I)
              IF((RBARWL1.LT.0.0).AND.(RBARWL2.LT.0.0)) THEN ! WATER LEVEL BELOW PIPE
                QN2(I)=QN2(I)+0.0D0                                   ! NO FLOW
                GOTO 2003
              ENDIF
              IF(RBARWL1.EQ.RBARWL2) THEN ! WATER LEVEL EQUAL ON BOTH SIDES OF PIPE
                QN2(I)=QN2(I)+0.0D0                ! NO FLOW
                GOTO 2003
              ENDIF
              IF(RBARWL1.GT.RBARWL2) THEN ! WATER LEVEL GREATER ON THIS SIDE OF PIPE
                IF(RBARWL2.LE.0) THEN ! OUTWARD FREE DISCHARGE
                  QN2(I)=QN2(I)-RAMP2*0.25D0*PI*(PIPEDIAM(I))**2
     &                    *(2.D0*G*RBARWL1/(1+PIPECOEF(I)))**0.5D0
                  GOTO 2003
                ELSE                         ! OUTWARD SUBMERGED DISCHARGE
                  QN2(I)=QN2(I)-RAMP2*0.25D0*PI*(PIPEDIAM(I))**2
     &                    *(2.D0*G*(RBARWL1-RBARWL2)/PIPECOEF(I))**0.5D0
                  GOTO 2003
                ENDIF
              ENDIF
              IF(RBARWL2.GT.RBARWL1) THEN !WATER LEVEL LOWER ON THIS SIDE OF PIPE
                IF(RBARWL1.LE.0) THEN ! INWARD FREE DISCHARGE
                  QN2(I)=QN2(I)+RAMP2*0.25D0*PI*(PIPEDIAM(I))**2
     &                    *(2.D0*G*RBARWL2/(1+PIPECOEF(I)))**0.5D0
                  GOTO 2003
                ELSE                         ! INWARD SUBMERGED DISCHARGE
                  QN2(I)=QN2(I)+RAMP2*0.25D0*PI*(PIPEDIAM(I))**2
     &                    *(2.D0*G*(RBARWL2-RBARWL1)/PIPECOEF(I))**0.5D0
                  GOTO 2003
                ENDIF
              ENDIF
2003          CONTINUE
            ENDIF
          END DO
        ENDIF
cjjwm001 - end add    

C......RESTART WIND AND PRESSURE INFORMATION

        IF(NWS.EQ.1) THEN
          OPEN(22,FILE=DIRNAME//'/'//'fort.22')
          DO J=1,ITHS
            DO I=1,NP
              READ(22,*) NHG,WSX2(I),WSY2(I),PR2(I)
              END DO
            END DO
          DO I=1,NP
            WSX2(I)=RAMP2*WSX2(I)
            WSY2(I)=RAMP2*WSY2(I)
            PR2(I)=RAMP2*PR2(I)
            END DO
          ENDIF

        IF(NWS.EQ.2) THEN
          OPEN(22,FILE=DIRNAME//'/'//'fort.22')
          WTIME1 = STATIM*86400.D0
          WTIME2 = WTIME1 + WTIMINC
          READ(22,*) (NHG,WVNX2(I),WVNY2(I),PRN2(I),I=1,NP)
          DO IT=1,ITHS
            TIMEIT=IT*DTDP + STATIM*86400.D0
            IF(TIMEIT.GT.WTIME2) THEN
              WTIME1=WTIME2
              WTIME2=WTIME2+WTIMINC
              DO I=1,NP
                WVNX1(I)=WVNX2(I)
                WVNY1(I)=WVNY2(I)
                PRN1(I)=PRN2(I)
                READ(22,*) NHG,WVNX2(I),WVNY2(I),PRN2(I)
              END DO
            ENDIF
          END DO
          WTRATIO=(TIME_A-WTIME1)/WTIMINC
          DO I=1,NP
            WINDX = WVNX1(I) + WTRATIO*(WVNX2(I)-WVNX1(I))
            WINDY = WVNY1(I) + WTRATIO*(WVNY2(I)-WVNY1(I))
            WSX2(I) = RAMP2*WINDX
            WSY2(I) = RAMP2*WINDY
            PR2(I)=RAMP2*(PRN1(I)+WTRATIO*(PRN2(I)-PRN1(I)))
            END DO
          ENDIF

        IF(NWS.EQ.-2) THEN
          OPEN(22,FILE=DIRNAME//'/'//'fort.22')
          WTIME1 = TIME_A
          WTIME2 = WTIME1 + WTIMINC
          READ(22,*) (NHG,WVNX1(I),WVNY1(I),PRN1(I),I=1,NP)
          READ(22,*) (NHG,WVNX2(I),WVNY2(I),PRN2(I),I=1,NP)
          DO I=1,NP
            WSX2(I) = RAMP2*WVNX1(I)
            WSY2(I) = RAMP2*WVNY1(I)
            PR2(I)=RAMP2*PRN1(I)
            END DO
          ENDIF

        IF(NWS.EQ.3) THEN
          OPEN(22,FILE=DIRNAME//'/'//'fort.22')
 2223     CALL NWS3GET( X, Y, SLAM, SFEA, WVNX2, WVNY2, IWTIME, IWYR,
     &                  WTIMED, NP, NWLON, NWLAT, WLATMAX, WLONMIN,
     &                  WLATINC, WLONINC, ICS, NSCREEN, ScreenUnit )
          IF(IWYR.NE.IREFYR) THEN
            IWTIMEP=IWTIME
            DO I=1,NP
              WVNX1(I)=WVNX2(I)
              WVNY1(I)=WVNY2(I)
              END DO
            GOTO 2223
            ENDIF
          IF(WTIMED.LE.WREFTIM) THEN
            IWTIMEP=IWTIME
            DO I=1,NP
              WVNX1(I)=WVNX2(I)
              WVNY1(I)=WVNY2(I)
              END DO
            GOTO 2223
            ENDIF
          IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) 
     &         WRITE(6,*)'FOUND WIND DATA AT TIME= ',IWTIMEP
          WRITE(16,*) 'FOUND WIND DATA AT TIME =',IWTIMEP
          IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) 
     &         WRITE(6,*)'FOUND WIND DATA AT TIME= ',IWTIME
          WRITE(16,*) 'FOUND WIND DATA AT TIME =',IWTIME
          WTIME2=WTIMED-WREFTIM                    !CAST INTO MODEL TIME REFERENCE
          WTIME1=WTIME2-WTIMINC
          DO IT=1,ITHS
            TIMEIT=IT*DTDP + STATIM*86400.D0
            IF(TIMEIT.GT.WTIME2) THEN
              WTIME1=WTIME2
              WTIME2=WTIME2+WTIMINC
              DO I=1,NP
                WVNX1(I)=WVNX2(I)
                WVNY1(I)=WVNY2(I)
                END DO
             CALL NWS3GET( X, Y, SLAM, SFEA, WVNX2, WVNY2, IWTIME, IWYR,
     &                  WTIMED, NP, NWLON, NWLAT, WLATMAX, WLONMIN,
     &                  WLATINC, WLONINC, ICS, NSCREEN, ScreenUnit )
              IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) 
     &         WRITE(6,*)'WIND FILE ADVANCED TO TIME',' = ', IWTIME
              WRITE(16,*) 'WIND FILE ADVANCED TO TIME = ',IWTIME
              ENDIF
            END DO
          WTRATIO=(TIME_A-WTIME1)/WTIMINC
          DO I=1,NP                                !INTERPOLATE IN TIME
            WINDX = WVNX1(I) + WTRATIO*(WVNX2(I)-WVNX1(I))
            WINDY = WVNY1(I) + WTRATIO*(WVNY2(I)-WVNY1(I))
            WINDMAG=SQRT(WINDX*WINDX+WINDY*WINDY)
            WDRAGCO = 0.001d0*(0.75d0+0.067d0*WINDMAG)
            IF(WDRAGCO.GT.0.003d0) WDRAGCO=0.003d0
            WSX2(I)=RAMP2*0.001293d0*WDRAGCO*WINDX*WINDMAG
            WSY2(I)=RAMP2*0.001293d0*WDRAGCO*WINDY*WINDMAG
            END DO
          ENDIF

        IF(NWS.EQ.4) THEN
          OPEN(22,FILE=DIRNAME//'/'//'fort.22')
          WTIME1 = STATIM*86400.D0
          WTIME2 = WTIME1 + WTIMINC
          CALL NWS4GET(WVNX2,WVNY2,PRN2,NP,RHOWAT0,G)
          DO IT=1,ITHS
            TIMEIT=IT*DTDP + STATIM*86400.D0
            IF(TIMEIT.GT.WTIME2) THEN
              WTIME1=WTIME2
              WTIME2=WTIME2+WTIMINC
              DO I=1,NP
                WVNX1(I)=WVNX2(I)
                WVNY1(I)=WVNY2(I)
                PRN1(I)=PRN2(I)
                END DO
              CALL NWS4GET(WVNX2,WVNY2,PRN2,NP,RHOWAT0,G)
              ENDIF
            END DO
          WTRATIO=(TIME_A-WTIME1)/WTIMINC
          DO I=1,NP
            WINDX = WVNX1(I) + WTRATIO*(WVNX2(I)-WVNX1(I))
            WINDY = WVNY1(I) + WTRATIO*(WVNY2(I)-WVNY1(I))
            WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
            WDRAGCO = 0.001d0*(0.75d0+0.067d0*WINDMAG)
            IF(WDRAGCO.GT.0.003d0) WDRAGCO=0.003d0
            WSX2(I) = RAMP2*0.001293d0*WDRAGCO*WINDX*WINDMAG
            WSY2(I) = RAMP2*0.001293d0*WDRAGCO*WINDY*WINDMAG
            PR2(I) = RAMP2*(PRN1(I)+WTRATIO*(PRN2(I)-PRN1(I)))
            END DO
          ENDIF

        IF(NWS.EQ.-4) THEN
          OPEN(22,FILE=DIRNAME//'/'//'fort.22')
          WTIME1 = TIME_A
          WTIME2 = WTIME1 + WTIMINC
          CALL NWS4GET(WVNX1,WVNY1,PRN1,NP,RHOWAT0,G)
          CALL NWS4GET(WVNX2,WVNY2,PRN2,NP,RHOWAT0,G)
          DO I=1,NP
            WINDX = WVNX1(I)
            WINDY = WVNY1(I)
            WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
            WDRAGCO = 0.001d0*(0.75d0+0.067d0*WINDMAG)
            IF(WDRAGCO.GT.0.003d0) WDRAGCO=0.003d0
            WSX2(I) = RAMP2*0.001293d0*WDRAGCO*WINDX*WINDMAG
            WSY2(I) = RAMP2*0.001293d0*WDRAGCO*WINDY*WINDMAG
            PR2(I) = RAMP2*PRN1(I)
            END DO
          ENDIF

        IF(NWS.EQ.5) THEN
          OPEN(22,FILE=DIRNAME//'/'//'fort.22')
          WTIME1 = STATIM*86400.D0
          WTIME2 = WTIME1 + WTIMINC
          READ(22,*) (NHG,WVNX2(I),WVNY2(I),PRN2(I),I=1,NP)
          DO IT=1,ITHS
            TIMEIT=IT*DTDP + STATIM*86400.D0
            IF(TIMEIT.GT.WTIME2) THEN
              WTIME1=WTIME2
              WTIME2=WTIME2+WTIMINC
              DO I=1,NP
                WVNX1(I)=WVNX2(I)
                WVNY1(I)=WVNY2(I)
                PRN1(I)=PRN2(I)
                READ(22,*) NHG,WVNX2(I),WVNY2(I),PRN2(I)
              END DO
            ENDIF
          END DO
          WTRATIO=(TIME_A-WTIME1)/WTIMINC
          DO I=1,NP
            WINDX = WVNX1(I) + WTRATIO*(WVNX2(I)-WVNX1(I))
            WINDY = WVNY1(I) + WTRATIO*(WVNY2(I)-WVNY1(I))
            WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
            WDRAGCO = 0.001d0*(0.75d0+0.067d0*WINDMAG)
            IF(WDRAGCO.GT.0.003d0) WDRAGCO=0.003d0
            WSX2(I) = RAMP2*0.001293d0*WDRAGCO*WINDX*WINDMAG
            WSY2(I) = RAMP2*0.001293d0*WDRAGCO*WINDY*WINDMAG
            PR2(I) = RAMP2*(PRN1(I)+WTRATIO*(PRN2(I)-PRN1(I))) 
            END DO
          ENDIF

        IF(NWS.EQ.-5) THEN
          OPEN(22,FILE=DIRNAME//'/'//'fort.22')
          WTIME1 = TIME_A
          WTIME2 = WTIME1 + WTIMINC
          READ(22,*) (NHG,WVNX1(I),WVNY1(I),PRN1(I),I=1,NP)
          READ(22,*) (NHG,WVNX2(I),WVNY2(I),PRN2(I),I=1,NP)
          DO I=1,NP
            WINDX = WVNX1(I)
            WINDY = WVNY1(I)
            WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
            WDRAGCO = 0.001d0*(0.75d0+0.067d0*WINDMAG)
            IF(WDRAGCO.GT.0.003d0) WDRAGCO=0.003d0
            WSX2(I) = RAMP2*0.001293d0*WDRAGCO*WINDX*WINDMAG
            WSY2(I) = RAMP2*0.001293d0*WDRAGCO*WINDY*WINDMAG
            PR2(I) = RAMP2*PRN1(I) 
            END DO
          ENDIF

        IF(NWS.EQ.6) THEN
          OPEN(22,FILE=DIRNAME//'/'//'fort.22')
C   The following 3 lines are a hardwire to allow a non standard met file to be read in at 
C   time zero in a hot start.  They should be eliminated or commented out for normal operation
c         OPEN(199,FILE=DIRNAME//'/'//'fort.199')
c         READ(199,*) (NHG,PRN1(I),WVNX1(I),WVNY1(I),I=1,NP)
c         CLOSE(199)
C   The following CALL statement should be uncommented for normal operation
          CALL NWS6GET(X,Y,SLAM,SFEA,WVNX1,WVNY1,PRN1,NP,NWLON,NWLAT,
     &                 WLATMAX,WLONMIN,WLATINC,WLONINC,ICS,RHOWAT0,G)
          CALL NWS6GET(X,Y,SLAM,SFEA,WVNX2,WVNY2,PRN2,NP,NWLON,NWLAT,
     &                 WLATMAX,WLONMIN,WLATINC,WLONINC,ICS,RHOWAT0,G)
          WTIME1=TIME_A
          WTIME2=WTIME1+WTIMINC
          WTRATIO=(TIME_A-WTIME1)/WTIMINC
          DO I=1,NP
            WINDX = WVNX1(I) + WTRATIO*(WVNX2(I)-WVNX1(I))
            WINDY = WVNY1(I) + WTRATIO*(WVNY2(I)-WVNY1(I))
            WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
            WDRAGCO = 0.001d0*(0.75d0+0.067d0*WINDMAG)
            IF(WDRAGCO.GT.0.003d0) WDRAGCO=0.003d0
            WSX2(I) = RAMP2*0.001293d0*WDRAGCO*WINDX*WINDMAG
            WSY2(I) = RAMP2*0.001293d0*WDRAGCO*WINDY*WINDMAG
            PR2(I) = RAMP2*(PRN1(I)+WTRATIO*(PRN2(I)-PRN1(I)))
            END DO
          ENDIF

        IF(NWS.EQ.10) THEN
          WTIME1=TIME_A
          WTIME2=WTIME1+WTIMINC
          NWSGGWI=-1
          CALL NWS10GET(NWSGGWI,SLAM,SFEA,WVNX1,WVNY1,PRN1,NP,RHOWAT0,G,
     &                  NWLON,NWLAT,WTIMINC)
          NWSGGWI=0
          CALL NWS10GET(NWSGGWI,SLAM,SFEA,WVNX1,WVNY1,PRN1,NP,RHOWAT0,G,
     &                  NWLON,NWLAT,WTIMINC)
          NWSGGWI=1
          CALL NWS10GET(NWSGGWI,SLAM,SFEA,WVNX2,WVNY2,PRN2,NP,RHOWAT0,G,
     &                  NWLON,NWLAT,WTIMINC)
          WTRATIO=(TIME_A-WTIME1)/WTIMINC
          DO I=1,NP
            WINDX = WVNX1(I) + WTRATIO*(WVNX2(I)-WVNX1(I))
            WINDY = WVNY1(I) + WTRATIO*(WVNY2(I)-WVNY1(I))
            WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
            WDRAGCO = 0.001d0*(0.75d0+0.067d0*WINDMAG)
            IF(WDRAGCO.GT.0.003d0) WDRAGCO=0.003d0
            WSX2(I) = RAMP2*0.001293d0*WDRAGCO*WINDX*WINDMAG
            WSY2(I) = RAMP2*0.001293d0*WDRAGCO*WINDY*WINDMAG
            PR2(I) = RAMP2*(PRN1(I)+WTRATIO*(PRN2(I)-PRN1(I)))
            END DO
          ENDIF

        IF(NWS.EQ.11) THEN
          WTIME1=TIME_A
          WTIME2=WTIME1+WTIMINC
          NWSEGWI=0
          IDSETFLG=0
          IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,1197)
          WRITE(16,1197)
          CALL NWS11GET(NWSEGWI,IDSETFLG,SLAM,SFEA,WVNX2,WVNY2,PRN2,NP,
     &                  RHOWAT0,G)  !JUST COMPUTE INTERPOLATING FACTORS
          IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,1198)
          WRITE(16,1198)
          NWSEGWI=0
          IDSETFLG=8
          CALL NWS11GET(NWSEGWI,IDSETFLG,SLAM,SFEA,WVNX1,WVNY1,PRN1,NP,
     &                  RHOWAT0,G)  !NOW COMPUTE HOTSTART WIND FILED
          NWSEGWI=1
          IDSETFLG=1
          CALL NWS11GET(NWSEGWI,IDSETFLG,SLAM,SFEA,WVNX2,WVNY2,PRN2,NP,
     &                  RHOWAT0,G)  !NOW COMPUTE NEXT WIND FIELD
          WTRATIO=(TIME_A-WTIME1)/WTIMINC
          DO I=1,NP
            WINDX = WVNX1(I) + WTRATIO*(WVNX2(I)-WVNX1(I))
            WINDY = WVNY1(I) + WTRATIO*(WVNY2(I)-WVNY1(I))
            WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
            WDRAGCO = 0.001d0*(0.75d0+0.067d0*WINDMAG)
            IF(WDRAGCO.GT.0.003d0) WDRAGCO=0.003d0
            WSX2(I) = RAMP2*0.001293d0*WDRAGCO*WINDX*WINDMAG
            WSY2(I) = RAMP2*0.001293d0*WDRAGCO*WINDY*WINDMAG
            PR2(I) = RAMP2*(PRN1(I)+WTRATIO*(PRN2(I)-PRN1(I)))
            END DO
        ENDIF
          
Cek added NWS = 12 and -12 from version 46
          
        IF (NWS.EQ.12) THEN
          WTIME1 = STATIM*86400.D0
          WTIME2 = WTIME1 + WTIMINC
          CALL NWS12INIT(WVNX1,WVNY1,PRN1,NP,RHOWAT0,G)
          CALL NWS12GET(WVNX2,WVNY2,PRN2,NP,RHOWAT0,G)
          DO IT = 1,ITHS
            TIMEIT = IT*DTDP + STATIM*86400.D0
            IF (TIMEIT.GT.WTIME2) THEN
              WTIME1 = WTIME2
              WTIME2 = WTIME2+WTIMINC
              DO I = 1,NP
                WVNX1(I) = WVNX2(I)
                WVNY1(I) = WVNY2(I)
                PRN1(I)  = PRN2(I)
              ENDDO
              CALL NWS12GET(WVNX2,WVNY2,PRN2,NP,RHOWAT0,G)
            ENDIF
          ENDDO
          WTRATIO = (TIME_A - WTIME1)/WTIMINC
          DO I = 1,NP
            WINDX = WVNX1(I) + WTRATIO*(WVNX2(I)-WVNX1(I))
            WINDY = WVNY1(I) + WTRATIO*(WVNY2(I)-WVNY1(I))
            WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
            WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
            IF (LoadDirEffRLen) THEN
              CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), H0, G, WindX, WindY)
              WindMag = SQRT(WindX*WindX+WindY*WindY)
              WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
            ENDIF
            WSX2(I)    = RampMete2*0.001293d0*WDRAGCO*WINDX*WINDMAG
            WSY2(I)    = RampMete2*0.001293d0*WDRAGCO*WINDY*WINDMAG
            PR2(I)     = RampMete2*PRN1(I)
            WVNXOUT(I) = RampMete2*WINDX
            WVNYOUT(I) = RampMete2*WINDY
          ENDDO
        ENDIF

        IF(NWS.EQ.-12) THEN
          WTIME1 = TIME_A
          WTIME2 = WTIME1 + WTIMINC
          CALL NWS12INIT(WVNX1,WVNY1,PRN1,NP,RHOWAT0,G)
          CALL NWS12GET(WVNX1,WVNY1,PRN1,NP,RHOWAT0,G)
          CALL NWS12GET(WVNX2,WVNY2,PRN2,NP,RHOWAT0,G)
          DO I = 1,NP
            WINDX   = WVNX1(I)
            WINDY   = WVNY1(I)
            WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
            WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
            IF (LoadDirEffRLen) THEN
              CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), H0, G, WindX, WindY)
              WindMag = SQRT(WindX*WindX+WindY*WindY)
              WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
            ENDIF
            WSX2(I)    = RampMete2*0.001293d0*WDRAGCO*WINDX*WINDMAG
            WSY2(I)    = RampMete2*0.001293d0*WDRAGCO*WINDY*WINDMAG
            PR2(I)     = RampMete2*PRN1(I)
            WVNXOUT(I) = RampMete*WINDX
            WVNYOUT(I) = RampMete*WINDY
          ENDDO
       ENDIF

C......RESTART THE WAVE RADIATION STRESS

        IF(NRS.EQ.1) THEN
          OPEN(23,FILE=DIRNAME//'/'//'fort.23')
          RSTIME1 = TIME_A
          RSTIME2 = RSTIME1 + RSTIMINC
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
          DO I=1,NP
            WSX2(I) = WSX2(I)+RAMP2*RSNX1(I)
            WSY2(I) = WSY2(I)+RAMP2*RSNY1(I)
            END DO
          ENDIF

       if (CTIP) then
Cjromo 11-01-00  Initialize TIP2 for HOTSTART
          DO I=1,NP
             TIP2(I)=0.0
          END DO
CTIP  LINES TO USE TIDAL POTENTIAL FORCING
       IF(NTIP.GE.1) THEN
         DO J=1,NTIF
           IF(PERT(J).EQ.0.) THEN
             NCYC=0
             ELSE
             NCYC=INT(TIME_A/PERT(J))
             ENDIF
           ARGT=AMIGT(J)*(TIME_A-NCYC*PERT(J))+FACET(J)
           TPMUL=RAMP2*ETRF(J)*TPK(J)*FFT(J)
           SALTMUL=RAMP2*FFT(J)
           NA=NINT(0.00014/AMIGT(J))
           IF(NA.EQ.1) THEN                        !SEMI-DIURNAL SPECIES
             DO I=1,NP
               ARGTP=ARGT+2.*SLAM(I)
               ARGSALT=ARGT-SALTPHA(J,I)
               CCSFEA=COS(SFEA(I))
               CCSFEA=CCSFEA*CCSFEA
               TIP2(I)=TIP2(I)+TPMUL*CCSFEA*COS(ARGTP)
     &                 +SALTMUL*SALTAMP(J,I)*COS(ARGSALT)
               END DO
             ENDIF
           IF(NA.EQ.2) THEN
             DO I=1,NP
               ARGTP=ARGT+SLAM(I)
               ARGSALT=ARGT-SALTPHA(J,I)
cjjw/vjpm002 - modified/added the following 5 lines
#ifdef REAL8
               S2SFEA=SIN(2.d0*SFEA(I))
#else
               S2SFEA=SIN(2.e0*SFEA(I))
#endif
               TIP2(I)=TIP2(I)+TPMUL*S2SFEA*COS(ARGTP)
     &                +SALTMUL*SALTAMP(J,I)*COS(ARGSALT)
               END DO
             ENDIF
           END DO
         ENDIF
      endif     !   CTIP


C...
C....SET UP TO RESTART TIMESERIES OUTPUT FILES
C....
C...
        IF(NBYTE.EQ.4) ITEMPSTP=20
        IF(NBYTE.EQ.8) ITEMPSTP=10

C...
C....IF RESTARTING THE ELEVATION STATION OUTPUT FILE, GO TO THE PROPER PLACE
C....IN THE FILE.  OTHERWISE ZERO OUT NSCOUE.
C...
        READ(IHOT,REC=IHOTSTP+1) IESTP
        READ(IHOT,REC=IHOTSTP+2) NSCOUE
        IHOTSTP=IHOTSTP+2
        WRITE(16,1040) IESTP,NSCOUE
 1040   FORMAT(//,1X,I6,' LINES OR RECORDS WRITTEN IN ELEVATION ',
     &                  'STATION FILE BY THE TIME OF THE HOT START',
     &             /,8X,'SPOOL COUNTER = ',I6)
        IF(NOUTE.LT.0) THEN
          IESTP=0
          NSCOUE=0        
          IF((NTCYSE.LT.ITHS).AND.(NSPOOLE.GT.0)) THEN
            NTCYSE=NTCYSE+((ITHS-NTCYSE)/NSPOOLE)*NSPOOLE
            IF(NTCYSE.LT.ITHS) NTCYSE=NTCYSE+NSPOOLE
            IF(NSPOOLE.NE.0) NTRSPE=(NTCYFE-NTCYSE)/NSPOOLE
            ENDIF
          WRITE(16,1041)
 1041     FORMAT(//,' A NEW ELEVATION STATION FILE WILL BE STARTED')
          ENDIF

        IF(NOUTE.EQ.-2) THEN
          OPEN(61,FILE=DIRNAME//'/'//'fort.61',
     &          ACCESS='DIRECT',RECL=NBYTE)
          IF (SEDFLAG.GE.1) OPEN(82,FILE=DIRNAME//'/'//'fort.82',
     &          ACCESS='DIRECT',RECL=NBYTE)
          IF(NBYTE.EQ.4) THEN
            DO I=1,8
              WRITE(61,REC=IESTP+I) RDES4(I)
              IF (SEDFLAG.GE.1) WRITE(82,REC=IESTP+I) RDES4(I)
              ENDDO
            IESTP=IESTP+8
            DO I=1,6
              WRITE(61,REC=IESTP+I) RID4(I)
              IF (SEDFLAG.GE.1) WRITE(82,REC=IESTP+I) RID4(I)
              ENDDO
            IESTP=IESTP+6
            DO I=1,6
              WRITE(61,REC=IESTP+I) AID4(I)
              IF (SEDFLAG.GE.1) WRITE(82,REC=IESTP+I) AID4(I)
              ENDDO
            IESTP=IESTP+6
            ENDIF
          IF(NBYTE.EQ.8) THEN
            DO I=1,4
              WRITE(61,REC=IESTP+I) RDES8(I)
              IF (SEDFLAG.GE.1) WRITE(82,REC=IESTP+I) RDES8(I)
              ENDDO
            IESTP=IESTP+4
            DO I=1,3
              WRITE(61,REC=IESTP+I) RID8(I)
              IF (SEDFLAG.GE.1) WRITE(82,REC=IESTP+I) RID8(I)
              ENDDO
            IESTP=IESTP+3
            DO I=1,3
              WRITE(61,REC=IESTP+I) AID8(I)
              IF (SEDFLAG.GE.1) WRITE(82,REC=IESTP+I) AID8(I)
              ENDDO
            IESTP=IESTP+3
            ENDIF
          WRITE(61,REC=IESTP+1) NTRSPE
          IF (SEDFLAG.GE.1) WRITE(82,REC=IESTP+1) NTRSPE
          WRITE(61,REC=IESTP+2) NSTAE
          IF (SEDFLAG.GE.1) WRITE(82,REC=IESTP+2) NSTAE
          WRITE(61,REC=IESTP+3) DT*NSPOOLE
          IF (SEDFLAG.GE.1) WRITE(82,REC=IESTP+3) DT*NSPOOLE
          WRITE(61,REC=IESTP+4) NSPOOLE
          IF (SEDFLAG.GE.1) WRITE(82,REC=IESTP+4) NSPOOLE
          WRITE(61,REC=IESTP+5) 1
          IF (SEDFLAG.GE.1) WRITE(82,REC=IESTP+5) 1
          IESTP=IESTP+5
          CLOSE(61)                    ! DO THIS TO FLUSH THE WRITE BUFFER
          IF (SEDFLAG.GE.1) CLOSE(82)                    ! DO THIS TO FLUSH THE WRITE BUFFER
          OPEN(61,FILE=DIRNAME//'/'//'fort.61',
     &         ACCESS='DIRECT',RECL=NBYTE)
          IF (SEDFLAG.GE.1) OPEN(82,FILE=DIRNAME//'/'//'fort.82',
     &         ACCESS='DIRECT',RECL=NBYTE)
          ENDIF
        IF(NOUTE.EQ.-1) THEN
          OPEN(61,FILE=DIRNAME//'/'//'fort.61')
          IF (SEDFLAG.GE.1) OPEN(82,FILE=DIRNAME//'/'//'fort.82')
          WRITE(61,3220) RUNDES,RUNID,AGRID
          IF (SEDFLAG.GE.1) WRITE(82,3220) RUNDES,RUNID,AGRID
          WRITE(61,3645) NTRSPE,NSTAE,DTDP*NSPOOLE,NSPOOLE,1
          IF (SEDFLAG.GE.1) WRITE(82,3645) NTRSPE,NSTAE,DTDP*NSPOOLE,
     &                                     NSPOOLE,1
          IESTP=2
          ENDIF
        IF(NOUTE.EQ.1) THEN
          OPEN(61,FILE=DIRNAME//'/'//'fort.61')
          IF (SEDFLAG.GE.1) OPEN(82,FILE=DIRNAME//'/'//'fort.82')
          DO I=1,IESTP          !I DON'T KNOW OF A PRACTICAL WAY TO CHANGE NTRSPE
 1050       FORMAT(1X)
            READ(61,1050)
            IF (SEDFLAG.GE.1) READ(82,1050)
            ENDDO
          ENDIF
        IF(NOUTE.EQ.2) THEN
          OPEN(61,FILE=DIRNAME//'/'//'fort.61',
     &          ACCESS='DIRECT',RECL=NBYTE)
          IF (SEDFLAG.GE.1) OPEN(82,FILE=DIRNAME//'/'//'fort.82',
     &          ACCESS='DIRECT',RECL=NBYTE)
          WRITE(61,REC=ITEMPSTP+1) NTRSPE   ! ALLOW ADDITIONAL OUTPUT DATA TO BE WRITTEN
          IF (SEDFLAG.GE.1) WRITE(82,REC=ITEMPSTP+1) NTRSPE   ! ALLOW ADDITIONAL OUTPUT DATA TO BE WRITTEN
          CLOSE(61)                         ! DO THIS TO FLUSH THE WRITE BUFFER
          IF (SEDFLAG.GE.1) CLOSE(82)                         ! DO THIS TO FLUSH THE WRITE BUFFER
          OPEN(61,FILE=DIRNAME//'/'//'fort.61',
     &         ACCESS='DIRECT',RECL=NBYTE)
          IF (SEDFLAG.GE.1) OPEN(82,FILE=DIRNAME//'/'//'fort.82',
     &         ACCESS='DIRECT',RECL=NBYTE)
          ENDIF

C...
C....GO TO THE PROPER PLACE IN THE VELOCITY STATION OUTPUT FILE
C...
        READ(IHOT,REC=IHOTSTP+1) IVSTP
        READ(IHOT,REC=IHOTSTP+2) NSCOUV
        IHOTSTP=IHOTSTP+2
        WRITE(16,1042) IVSTP,NSCOUV
 1042   FORMAT(//,1X,I6,' LINES OR RECORDS WRITTEN IN VELOCITY ',
     &                  'STATION FILE BY THE TIME OF THE HOT START',
     &          /,8X,'SPOOL COUNTER =',I6)
        IF(NOUTV.LT.0) THEN
          IVSTP=0
          NSCOUV=0        
          IF((NTCYSV.LT.ITHS).AND.(NSPOOLV.GT.0)) THEN
            NTCYSV=NTCYSV+((ITHS-NTCYSV)/NSPOOLV)*NSPOOLV
            IF(NTCYSV.LT.ITHS) NTCYSV=NTCYSV+NSPOOLV
            NTRSPV=(NTCYFV-NTCYSV)/NSPOOLV
            ENDIF
          WRITE(16,1043)
 1043     FORMAT(//,' A NEW VELOCITY STATION FILE WILL BE STARTED')
          ENDIF

        IF(NOUTV.EQ.-2) THEN
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
        IF(NOUTV.EQ.-1) THEN
          OPEN(62,FILE=DIRNAME//'/'//'fort.62')
          WRITE(62,3220) RUNDES,RUNID,AGRID
          WRITE(62,3645) NTRSPV,NSTAV,DTDP*NSPOOLV,NSPOOLV,2
          IVSTP=2
          ENDIF
        IF(NOUTV.EQ.1) THEN
          OPEN(62,FILE=DIRNAME//'/'//'fort.62')
          DO I=1,IVSTP          !I DON'T KNOW OF A PRACTICAL WAY TO CHANGE NTRSPV
            READ(62,1050)
            ENDDO
          ENDIF
        IF(NOUTV.EQ.2) THEN
          OPEN(62,FILE=DIRNAME//'/'//'fort.62',
     &          ACCESS='DIRECT',RECL=NBYTE)
          WRITE(62,REC=ITEMPSTP+1) NTRSPV   ! ALLOW ADDITIONAL OUTPUT DATA TO BE WRITTEN
          CLOSE(62)                         ! DO THIS TO FLUSH THE WRITE BUFFER
          OPEN(62,FILE=DIRNAME//'/'//'fort.62',
     &         ACCESS='DIRECT',RECL=NBYTE)
          ENDIF

C...
C....GO TO THE PROPER PLACE IN THE CONCENTRATION STATION OUTPUT FILE
C...
        READ(IHOT,REC=IHOTSTP+1) ICSTP
        READ(IHOT,REC=IHOTSTP+2) NSCOUC
        IHOTSTP=IHOTSTP+2
        WRITE(16,1044) ICSTP,NSCOUC
 1044   FORMAT(//,1X,I6,' LINES OR RECORDS WRITTEN IN CONCENTRATION ',
     &                  'STATION FILE BY THE TIME OF THE HOT START',
     &          /,8X,'SPOOL COUNTER = ',I6)
        IF(NOUTC.LT.0) THEN
          ICSTP=0
          NSCOUC=0        
          IF((NTCYSC.LT.ITHS).AND.(NSPOOLC.GT.0)) THEN
            NTCYSC=NTCYSC+((ITHS-NTCYSC)/NSPOOLC)*NSPOOLC
            IF(NTCYSC.LT.ITHS) NTCYSC=NTCYSC+NSPOOLC
            NTRSPC=(NTCYFC-NTCYSC)/NSPOOLC
            ENDIF
          WRITE(16,1045)
 1045     FORMAT(//,' A NEW CONCENTRATION STATION FILE WILL BE STARTED')
          ENDIF

        IF(NOUTC.EQ.-2) THEN
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
        IF(NOUTC.EQ.-1) THEN
          OPEN(81,FILE='fort.81')
          WRITE(81,3220) RUNDES,RUNID,AGRID
          WRITE(81,3645) NTRSPC,NSTAC,DTDP*NSPOOLC,NSPOOLC,1
          ICSTP=2
          ENDIF
        IF(NOUTC.EQ.1) THEN
          OPEN(81,FILE=DIRNAME//'/'//'fort.81')
          DO I=1,ICSTP          !I DON'T KNOW OF A PRACTICAL WAY TO CHANGE NTRSPC
            READ(81,1050)
            ENDDO
          ENDIF
        IF(NOUTC.EQ.2) THEN
          OPEN(81,FILE=DIRNAME//'/'//'fort.81',
     &          ACCESS='DIRECT',RECL=NBYTE)
crevisit          WRITE(81,REC=ITEMPSTP+1) NTRSPC   ! ALLOW ADDITIONAL OUTPUT DATA TO BE WRITTEN
          CLOSE(81)                         ! DO THIS TO FLUSH THE WRITE BUFFER
          OPEN(81,FILE=DIRNAME//'/'//'fort.81',
     &         ACCESS='DIRECT',RECL=NBYTE)
          ENDIF

C...
C....GO TO THE PROPER PLACE IN THE METEOROLOGICAL STATION OUTPUT FILE
C...
        READ(IHOT,REC=IHOTSTP+1) IPSTP
        READ(IHOT,REC=IHOTSTP+2) IWSTP
        READ(IHOT,REC=IHOTSTP+3) NSCOUM
        IHOTSTP=IHOTSTP+3
        WRITE(16,1038) IWSTP,IPSTP,NSCOUM
 1038   FORMAT(//,1X,I6,' LINES OR RECORDS WRITTEN IN THE WIND STATION',
     &                  ' FILE BY THE TIME OF THE HOT START',
     &          /,1X,I6,' LINES OR RECORDS WRITTIN IN THE PRES STATION',
     &                  ' FILE BY THE TMIE OF THE HOT START',
     &          /,8X,'SPOOL COUNTER = ',I6)
        IF(NOUTM.LT.0) THEN
          IPSTP=0
          IWSTP=0
          NSCOUM=0        
          IF((NTCYSM.LT.ITHS).AND.(NSPOOLM.GT.0)) THEN
            NTCYSM=NTCYSM+((ITHS-NTCYSM)/NSPOOLM)*NSPOOLM
            IF(NTCYSM.LT.ITHS) NTCYSM=NTCYSM+NSPOOLM
            NTRSPM=(NTCYFM-NTCYSM)/NSPOOLM
            ENDIF
          WRITE(16,1039)
 1039    FORMAT(//,' A NEW METEOROLOGICAL STATION FILE WILL BE STARTED')
          ENDIF

        IF(NOUTM.EQ.-2) THEN
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
          CLOSE(72)
          OPEN(71,FILE=DIRNAME//'/'//'fort.71',
     &         ACCESS='DIRECT',RECL=NBYTE)          
          OPEN(72,FILE=DIRNAME//'/'//'fort.72',
     &         ACCESS='DIRECT',RECL=NBYTE)
          ENDIF
        IF(NOUTM.EQ.-1) THEN
          OPEN(71,FILE=DIRNAME//'/'//'fort.71',
     &         ACCESS='DIRECT',RECL=NBYTE)          
          OPEN(72,FILE=DIRNAME//'/'//'fort.72',
     &         ACCESS='DIRECT',RECL=NBYTE)
          WRITE(71,3220) RUNDES,RUNID,AGRID
          WRITE(71,3645) NTRSPM,NSTAM,DTDP*NSPOOLM,NSPOOLM,1
          WRITE(72,3220) RUNDES,RUNID,AGRID
          WRITE(72,3645) NTRSPM,NSTAM,DTDP*NSPOOLM,NSPOOLM,1
          IPSTP=2
          IWSTP=2
          ENDIF
        IF(NOUTM.EQ.1) THEN
          OPEN(71,FILE=DIRNAME//'/'//'fort.71')
          OPEN(72,FILE=DIRNAME//'/'//'fort.72')
          DO I=1,IPSTP          !I DON'T KNOW OF A PRACTICAL WAY TO CHANGE NTRSPM
            READ(71,1050)
            ENDDO
          DO I=1,IWSTP          !I DON'T KNOW OF A PRACTICAL WAY TO CHANGE NTRSPM
            READ(72,1050)
            ENDDO
          ENDIF
        IF(NOUTM.EQ.2) THEN
          OPEN(71,FILE=DIRNAME//'/'//'fort.71',
     &          ACCESS='DIRECT',RECL=NBYTE)
          OPEN(72,FILE=DIRNAME//'/'//'fort.72',
     &          ACCESS='DIRECT',RECL=NBYTE)
crevisit          WRITE(71,REC=ITEMPSTP+1) NTRSPM   ! ALLOW ADDITIONAL OUTPUT DATA TO BE WRITTEN
crevisit          WRITE(72,REC=ITMEPSTP+1) NTRSPM
          CLOSE(71)                         ! DO THIS TO FLUSH THE WRITE BUFFER
          CLOSE(72)
          OPEN(71,FILE=DIRNAME//'/'//'fort.71',
     &         ACCESS='DIRECT',RECL=NBYTE)
          OPEN(72,FILE=DIRNAME//'/'//'fort.72',
     &         ACCESS='DIRECT',RECL=NBYTE)
          ENDIF

C...
C....GO TO THE PROPER PLACE IN THE GLOBAL ELEVATION OUTPUT FILE
C...
        READ(IHOT,REC=IHOTSTP+1) IGEP
        READ(IHOT,REC=IHOTSTP+2) NSCOUGE
        IHOTSTP=IHOTSTP+2
        WRITE(16,1046) IGEP,NSCOUGE
 1046   FORMAT(//,1X,I6,' LINES OR RECORDS WRITTEN IN THE GLOBAL ',
     &                  'ELEVATION FILE BY THE TIME OF THE HOT START',
     &          /,8X,'SPOOL COUNTER =',I6)
        IF(NOUTGE.LT.0) THEN
          IGEP=0
          NSCOUGE=0        
          IF((NTCYSGE.LT.ITHS).AND.(NSPOOLGE.GT.0)) THEN
            NTCYSGE=NTCYSGE+((ITHS-NTCYSGE)/NSPOOLGE)*NSPOOLGE
            IF(NTCYSGE.LT.ITHS) NTCYSGE=NTCYSGE+NSPOOLGE
            NDSETSE=(NTCYFGE-NTCYSGE)/NSPOOLGE
            ENDIF
          WRITE(16,1047)
 1047     FORMAT(//,' A NEW GLOBAL ELEVATION FILE WILL BE STARTED')
          ENDIF

        IF(NOUTGE.EQ.-2) THEN
          OPEN(63,FILE=DIRNAME//'/'//'fort.63',
     &          ACCESS='DIRECT',RECL=NBYTE)
          IF (SEDFLAG.GE.1) OPEN(84,FILE=DIRNAME//'/'//'fort.84',
     &          ACCESS='DIRECT',RECL=NBYTE)
          IF (SEDFLAG.GE.1) OPEN(85,FILE=DIRNAME//'/'//'fort.85',
     &          ACCESS='DIRECT',RECL=NBYTE)
          IF(NBYTE.EQ.4) THEN
            DO I=1,8
              WRITE(63,REC=IGEP+I) RDES4(I)
              IF (SEDFLAG.GE.1) WRITE(84,REC=IGEP+I) RDES4(I)
              IF (SEDFLAG.GE.1) WRITE(85,REC=IGEP+I) RDES4(I)
              ENDDO
            IGEP=IGEP+8
            DO I=1,6
              WRITE(63,REC=IGEP+I) RID4(I)
              IF (SEDFLAG.GE.1) WRITE(84,REC=IGEP+I) RID4(I)
              IF (SEDFLAG.GE.1) WRITE(85,REC=IGEP+I) RID4(I)
              ENDDO
            IGEP=IGEP+6
            DO I=1,6
              WRITE(63,REC=IGEP+I) AID4(I)
              IF (SEDFLAG.GE.1) WRITE(84,REC=IGEP+I) AID4(I)
              IF (SEDFLAG.GE.1) WRITE(85,REC=IGEP+I) AID4(I)
              ENDDO
            IGEP=IGEP+6
            ENDIF
          IF(NBYTE.EQ.8) THEN
            DO I=1,4
              WRITE(63,REC=IGEP+I) RDES8(I)
              IF (SEDFLAG.GE.1) WRITE(84,REC=IGEP+I) RDES8(I)
              IF (SEDFLAG.GE.1) WRITE(85,REC=IGEP+I) RDES8(I)
              ENDDO
            IGEP=IGEP+4
            DO I=1,3
              WRITE(63,REC=IGEP+I) RID8(I)
              IF (SEDFLAG.GE.1) WRITE(84,REC=IGEP+I) RID8(I)
              IF (SEDFLAG.GE.1) WRITE(85,REC=IGEP+I) RID8(I)
              ENDDO
            IGEP=IGEP+3
            DO I=1,3
              WRITE(63,REC=IGEP+I) AID8(I)
              IF (SEDFLAG.GE.1) WRITE(84,REC=IGEP+I) AID8(I)
              IF (SEDFLAG.GE.1) WRITE(85,REC=IGEP+I) AID8(I)
              ENDDO
            IGEP=IGEP+3
            ENDIF
          WRITE(63,REC=IGEP+1) NDSETSE
          IF (SEDFLAG.GE.1) WRITE(84,REC=IGEP+1) NDSETSE
          IF (SEDFLAG.GE.1) WRITE(85,REC=IGEP+1) NDSETSE
          WRITE(63,REC=IGEP+2) NP
          IF (SEDFLAG.GE.1) WRITE(84,REC=IGEP+2) NP
          IF (SEDFLAG.GE.1) WRITE(85,REC=IGEP+2) NP
          WRITE(63,REC=IGEP+3) DT*NSPOOLGE
          IF (SEDFLAG.GE.1) WRITE(84,REC=IGEP+3) DT*NSPOOLGE
          IF (SEDFLAG.GE.1) WRITE(85,REC=IGEP+3) DT*NSPOOLGE
          WRITE(63,REC=IGEP+4) NSPOOLGE
          IF (SEDFLAG.GE.1) WRITE(84,REC=IGEP+4) NSPOOLGE
          IF (SEDFLAG.GE.1) WRITE(85,REC=IGEP+4) NSPOOLGE
          WRITE(63,REC=IGEP+5) 1
          IF (SEDFLAG.GE.1) WRITE(84,REC=IGEP+5) 1
          IF (SEDFLAG.GE.1) WRITE(85,REC=IGEP+5) 1
          IGEP=IGEP+5
          CLOSE(63)                    ! DO THIS TO FLUSH THE WRITE BUFFER
          IF (SEDFLAG.GE.1) CLOSE(84)  ! DO THIS TO FLUSH THE WRITE BUFFER
          IF (SEDFLAG.GE.1) CLOSE(85)  ! DO THIS TO FLUSH THE WRITE BUFFER
          OPEN(63,FILE=DIRNAME//'/'//'fort.63',
     &         ACCESS='DIRECT',RECL=NBYTE)
          IF (SEDFLAG.GE.1) OPEN(84,FILE=DIRNAME//'/'//'fort.84',
     &         ACCESS='DIRECT',RECL=NBYTE)
          IF (SEDFLAG.GE.1) OPEN(85,FILE=DIRNAME//'/'//'fort.85',
     &         ACCESS='DIRECT',RECL=NBYTE)
          ENDIF
        IF(NOUTGE.EQ.-1) THEN
          OPEN(63,FILE=DIRNAME//'/'//'fort.63')
          IF (SEDFLAG.GE.1) OPEN(84,FILE=DIRNAME//'/'//'fort.84')
          IF (SEDFLAG.GE.1) OPEN(85,FILE=DIRNAME//'/'//'fort.85')
          WRITE(63,3220) RUNDES,RUNID,AGRID
          IF (SEDFLAG.GE.1) WRITE(84,3220) RUNDES,RUNID,AGRID
          IF (SEDFLAG.GE.1) WRITE(85,3220) RUNDES,RUNID,AGRID
Casey 130211: Debug.
         !WRITE(63,3645) NDSETSE,NP,DTDP*NSPOOLGE,NSPOOLGE,1
          WRITE(63,3645) NDSETSE,NE,DTDP*NSPOOLGE,NSPOOLGE,1
          IF (SEDFLAG.GE.1) WRITE(84,3645) NDSETSE,NP,DTDP*NSPOOLGE,
     &                                     NSPOOLGE,1
          IF (SEDFLAG.GE.1) WRITE(85,3645) NDSETSE,NP,DTDP*NSPOOLGE,
     &                                     NSPOOLGE,1
          IGEP=2
          ENDIF
        IF(NOUTGE.EQ.1) THEN
          OPEN(63,FILE=DIRNAME//'/'//'fort.63')
          IF (SEDFLAG.GE.1) OPEN(84,FILE=DIRNAME//'/'//'fort.84')
          IF (SEDFLAG.GE.1) OPEN(85,FILE=DIRNAME//'/'//'fort.85')
          DO I=1,IGEP           !I DON'T KNOW OF A PRACTICAL WAY TO CHANGE NDSETSE
            READ(63,1050)
            IF (SEDFLAG.GE.1) READ(84,1050)
            IF (SEDFLAG.GE.1) READ(85,1050)
            ENDDO
          ENDIF
        IF(NOUTGE.EQ.2) THEN
          OPEN(63,FILE=DIRNAME//'/'//'fort.63',
     &          ACCESS='DIRECT',RECL=NBYTE)
          IF (SEDFLAG.GE.1) OPEN(84,FILE=DIRNAME//'/'//'fort.84',
     &          ACCESS='DIRECT',RECL=NBYTE)
          IF (SEDFLAG.GE.1) OPEN(85,FILE=DIRNAME//'/'//'fort.85',
     &          ACCESS='DIRECT',RECL=NBYTE)
          WRITE(63,REC=ITEMPSTP+1) NDSETSE  ! ALLOW ADDITIONAL OUTPUT DATA TO BE WRITTEN
          IF (SEDFLAG.GE.1) WRITE(84,REC=ITEMPSTP+1) NDSETSE  ! ALLOW ADDITIONAL OUTPUT DATA TO BE WRITTEN
          CLOSE(63)                         ! DO THIS TO FLUSH THE WRITE BUFFER
          IF (SEDFLAG.GE.1) CLOSE(84)                         ! DO THIS TO FLUSH THE WRITE BUFFER
          IF (SEDFLAG.GE.1) CLOSE(85)                         ! DO THIS TO FLUSH THE WRITE BUFFER
          OPEN(63,FILE=DIRNAME//'/'//'fort.63',
     &         ACCESS='DIRECT',RECL=NBYTE)
          IF (SEDFLAG.GE.1) OPEN(84,FILE=DIRNAME//'/'//'fort.84',
     &         ACCESS='DIRECT',RECL=NBYTE)
          IF (SEDFLAG.GE.1) OPEN(85,FILE=DIRNAME//'/'//'fort.85',
     &         ACCESS='DIRECT',RECL=NBYTE)
          ENDIF

C...
C....GO TO THE PROPER PLACE IN THE GLOBAL VELOCITY OUTPUT FILE
C...
        READ(IHOT,REC=IHOTSTP+1) IGVP
        READ(IHOT,REC=IHOTSTP+2) NSCOUGV
        IHOTSTP=IHOTSTP+2
        WRITE(16,1048) IGVP,NSCOUGV
 1048   FORMAT(//,1X,I6,' LINES OR RECORDS WRITTEN IN THE GLOBAL ',
     &                  'VELOCITY FILE BY THE TIME OF THE HOT START',
     &          /,8X,'SPOOL COUNTER =',I6)
        IF(NOUTGV.LT.0) THEN
          IGVP=0
          NSCOUGV=0        
          IF((NTCYSGV.LT.ITHS).AND.(NSPOOLGV.GT.0)) THEN
            NTCYSGV=NTCYSGV+((ITHS-NTCYSGV)/NSPOOLGV)*NSPOOLGV
            IF(NTCYSGV.LT.ITHS) NTCYSGV=NTCYSGV+NSPOOLGV
            NDSETSV=(NTCYFGV-NTCYSGV)/NSPOOLGV
            ENDIF
          WRITE(16,1049)
 1049     FORMAT(//,' A NEW GLOBAL VELOCITY FILE WILL BE STARTED')
          ENDIF

        IF(NOUTGV.EQ.-2) THEN
          OPEN(64,FILE=DIRNAME//'/'//'fort.64',
     &          ACCESS='DIRECT',RECL=NBYTE)
          IF (SEDFLAG.GE.1) OPEN(94,FILE=DIRNAME//'/'//'fort.94',
     &          ACCESS='DIRECT',RECL=NBYTE)
          IF (SEDFLAG.GE.1) OPEN(96,FILE=DIRNAME//'/'//'fort.96',
     &          ACCESS='DIRECT',RECL=NBYTE)
          IF(NBYTE.EQ.4) THEN
            DO I=1,8
              WRITE(64,REC=IGVP+I) RDES4(I)
              IF (SEDFLAG.GE.1) WRITE(94,REC=IGVP+I) RDES4(I)
              IF (SEDFLAG.GE.1) WRITE(96,REC=IGVP+I) RDES4(I)
              ENDDO
            IGVP=IGVP+8
            DO I=1,6
              WRITE(64,REC=IGVP+I) RID4(I)
              IF (SEDFLAG.GE.1) WRITE(94,REC=IGVP+I) RID4(I)
              IF (SEDFLAG.GE.1) WRITE(96,REC=IGVP+I) RID4(I)
              ENDDO
            IGVP=IGVP+6
            DO I=1,6
              WRITE(64,REC=IGVP+I) AID4(I)
              IF (SEDFLAG.GE.1) WRITE(94,REC=IGVP+I) AID4(I)
              IF (SEDFLAG.GE.1) WRITE(96,REC=IGVP+I) AID4(I)
              ENDDO
            IGVP=IGVP+6
            ENDIF
          IF(NBYTE.EQ.8) THEN
            DO I=1,4
              WRITE(64,REC=IGVP+I) RDES8(I)
              IF (SEDFLAG.GE.1) WRITE(94,REC=IGVP+I) RDES8(I)
              IF (SEDFLAG.GE.1) WRITE(96,REC=IGVP+I) RDES8(I)
              ENDDO
            IGVP=IGVP+4
            DO I=1,3
              WRITE(64,REC=IGVP+I) RID8(I)
              IF (SEDFLAG.GE.1) WRITE(94,REC=IGVP+I) RID8(I)
              IF (SEDFLAG.GE.1) WRITE(96,REC=IGVP+I) RID8(I)
              ENDDO
            IGVP=IGVP+3
            DO I=1,3
              WRITE(64,REC=IGVP+I) AID8(I)
              IF (SEDFLAG.GE.1) WRITE(94,REC=IGVP+I) AID8(I)
              IF (SEDFLAG.GE.1) WRITE(96,REC=IGVP+I) AID8(I)
              ENDDO
            IGVP=IGVP+3
            ENDIF
          WRITE(64,REC=IGVP+1) NDSETSV
          IF (SEDFLAG.GE.1) WRITE(94,REC=IGVP+1) NDSETSV
          IF (SEDFLAG.GE.1) WRITE(96,REC=IGVP+1) NDSETSV
          WRITE(64,REC=IGVP+2) NP
          IF (SEDFLAG.GE.1) WRITE(94,REC=IGVP+2) NP
          IF (SEDFLAG.GE.1) WRITE(96,REC=IGVP+2) NP
          WRITE(64,REC=IGVP+3) DT*NSPOOLGV
          IF (SEDFLAG.GE.1) WRITE(94,REC=IGVP+3) DT*NSPOOLGV
          IF (SEDFLAG.GE.1) WRITE(96,REC=IGVP+3) DT*NSPOOLGV
          WRITE(64,REC=IGVP+4) NSPOOLGV
          IF (SEDFLAG.GE.1) WRITE(94,REC=IGVP+4) NSPOOLGV
          IF (SEDFLAG.GE.1) WRITE(96,REC=IGVP+4) NSPOOLGV
          WRITE(64,REC=IGVP+5) 2
          IF (SEDFLAG.GE.1) WRITE(94,REC=IGVP+5) 2
          IF (SEDFLAG.GE.1) WRITE(96,REC=IGVP+5) 2
          IGVP=IGVP+5
          CLOSE(64)                    ! DO THIS TO FLUSH THE WRITE BUFFER
          IF (SEDFLAG.GE.1) CLOSE(94)
          IF (SEDFLAG.GE.1) CLOSE(96)
          OPEN(64,FILE=DIRNAME//'/'//'fort.64',
     &         ACCESS='DIRECT',RECL=NBYTE)
          IF (SEDFLAG.GE.1)  OPEN(94,FILE=DIRNAME//'/'//'fort.94',
     &         ACCESS='DIRECT',RECL=NBYTE)
          IF (SEDFLAG.GE.1) OPEN(96,FILE=DIRNAME//'/'//'fort.96',
     &         ACCESS='DIRECT',RECL=NBYTE)
          ENDIF
       IF(NOUTGV.EQ.-1) THEN
          OPEN(64,FILE=DIRNAME//'/'//'fort.64')
          IF (SEDFLAG.GE.1) OPEN(94,FILE=DIRNAME//'/'//'fort.94')
          IF (SEDFLAG.GE.1) OPEN(96,FILE=DIRNAME//'/'//'fort.96')
          WRITE(64,3220) RUNDES,RUNID,AGRID
          IF (SEDFLAG.GE.1) WRITE(94,3220) RUNDES,RUNID,AGRID
          IF (SEDFLAG.GE.1) WRITE(96,3220) RUNDES,RUNID,AGRID
Casey 130211: Debug.
         !WRITE(64,3645) NDSETSV,NP,DTDP*NSPOOLGV,NSPOOLGV,2
          WRITE(64,3645) NDSETSV,NE,DTDP*NSPOOLGV,NSPOOLGV,2
          IF (SEDFLAG.GE.1) WRITE(94,3645) NDSETSV,NP,DTDP*NSPOOLGV,
     &                                     NSPOOLGV,2
          IF (SEDFLAG.GE.1) WRITE(96,3645) NDSETSV,NP,DTDP*NSPOOLGV,
     &                                     NSPOOLGV,2
          IGVP=2
          ENDIF
        IF(NOUTGV.EQ.1) THEN
          OPEN(64,FILE=DIRNAME//'/'//'fort.64')
          IF (SEDFLAG.GE.1) OPEN(94,FILE=DIRNAME//'/'//'fort.94')
          IF (SEDFLAG.GE.1) OPEN(96,FILE=DIRNAME//'/'//'fort.96')
          DO I=1,IGVP           !I DON'T KNOW OF A PRACTICAL WAY TO CHANGE NDSETSV
            READ(64,1050)
            IF (SEDFLAG.GE.1) READ(94,1050)
            IF (SEDFLAG.GE.1) READ(96,1050)
            ENDDO
          ENDIF
        IF(NOUTGV.EQ.2) THEN
          OPEN(64,FILE=DIRNAME//'/'//'fort.64',
     &          ACCESS='DIRECT',RECL=NBYTE)
          IF (SEDFLAG.GE.1) OPEN(94,FILE=DIRNAME//'/'//'fort.94',
     &          ACCESS='DIRECT',RECL=NBYTE)
          IF (SEDFLAG.GE.1) OPEN(96,FILE=DIRNAME//'/'//'fort.96',
     &          ACCESS='DIRECT',RECL=NBYTE)
          WRITE(64,REC=ITEMPSTP+1) NDSETSV  ! ALLOW ADDITIONAL OUTPUT DATA TO BE WRITTEN
          IF (SEDFLAG.GE.1) WRITE(94,REC=ITEMPSTP+1) NDSETSV
          IF (SEDFLAG.GE.1) WRITE(96,REC=ITEMPSTP+1) NDSETSV
          CLOSE(64)
          IF (SEDFLAG.GE.1) CLOSE(94)
          IF (SEDFLAG.GE.1) CLOSE(96)                          ! DO THIS TO FLUSH THE WRITE BUFFER
          OPEN(64,FILE=DIRNAME//'/'//'fort.64',
     &         ACCESS='DIRECT',RECL=NBYTE)
          IF (SEDFLAG.GE.1) OPEN(94,FILE=DIRNAME//'/'//'fort.94',
     &         ACCESS='DIRECT',RECL=NBYTE)
          IF (SEDFLAG.GE.1) OPEN(96,FILE=DIRNAME//'/'//'fort.96',
     &         ACCESS='DIRECT',RECL=NBYTE)
          ENDIF

C...
C....GO TO THE PROPER PLACE IN THE GLOBAL CONCENTRATION OUTPUT FILE
C...
        READ(IHOT,REC=IHOTSTP+1) IGCP
        READ(IHOT,REC=IHOTSTP+2) NSCOUGC
        IHOTSTP=IHOTSTP+2
        WRITE(16,1053) IGCP,NSCOUGC
 1053   FORMAT(//,1X,I6,' LINES OR RECORDS WRITTEN IN THE GLOBAL ',
     &                'CONCENTRATION FILE BY THE TIME OF THE HOT START',
     &          /,8X,'SPOOL COUNTER =',I6)
        IF(NOUTGC.LT.0) THEN
          IGCP=0
          NSCOUGC=0        
          IF((NTCYSGC.LT.ITHS).AND.(NSPOOLGC.GT.0)) THEN
            NTCYSGC=NTCYSGC+((ITHS-NTCYSGC)/NSPOOLGC)*NSPOOLGC
            IF(NTCYSGC.LT.ITHS) NTCYSGC=NTCYSGC+NSPOOLGC
            NDSETSC=(NTCYFGC-NTCYSGC)/NSPOOLGC
            ENDIF
          WRITE(16,1054)
 1054     FORMAT(//,' A NEW GLOBAL CONCENTRATION FILE WILL BE STARTED')
          ENDIF

         IF(NOUTGC.EQ.-2) THEN
          OPEN(83,FILE=DIRNAME//'/'//'fort.83',
     &           ACCESS='DIRECT',RECL=NBYTE)
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
        IF(NOUTGC.EQ.-1) THEN
          OPEN(83,FILE=DIRNAME//'/'//'fort.83')
          WRITE(83,3220) RUNDES,RUNID,AGRID
          WRITE(83,3645) NDSETSC,NP,DTDP*NSPOOLGC,NSPOOLGC,1
          IGCP=2
          ENDIF
        IF(NOUTGC.EQ.1) THEN
          OPEN(83,FILE=DIRNAME//'/'//'fort.83')
          DO I=1,IGCP           !I DON'T KNOW OF A PRACTICAL WAY TO CHANGE NDSETSC
            READ(83,1050)
            ENDDO
          ENDIF
        IF(NOUTGC.EQ.2) THEN
          OPEN(83,FILE=DIRNAME//'/'//'fort.83',
     &          ACCESS='DIRECT',RECL=NBYTE)
          WRITE(83,REC=ITEMPSTP+1) NDSETSC  ! ALLOW ADDITIONAL OUTPUT DATA TO BE WRITTEN
          CLOSE(83)                         ! DO THIS TO FLUSH THE WRITE BUFFER
          OPEN(83,FILE=DIRNAME//'/'//'fort.83',
     &         ACCESS='DIRECT',RECL=NBYTE)
          ENDIF

C...
C....GO TO THE PROPER PLACE IN THE GLOBAL METEOROLOGICAL OUTPUT FILES
C...
        READ(IHOT,REC=IHOTSTP+1) IGPP
        READ(IHOT,REC=IHOTSTP+2) IGWP
        READ(IHOT,REC=IHOTSTP+3) NSCOUGW
        IHOTSTP=IHOTSTP+3
        WRITE(16,1055) IGWP,IGPP,NSCOUGW
 1055   FORMAT(//,1X,I6,' LINES OR RECORDS WRITTEN IN THE GLOBAL ',
     &                  'WIND FILE BY THE TIME OF THE HOT START',
     &          /,1X,I6,'LINES OR RECORDS WRITTEN IN THE GLOBAL ',
     &                  'PRESSURE FILE BY THE TIME OF THE HOT START',   
     &          /,8X,'SPOOL COUNTER =',I6)
        IF(NOUTGW.LT.0) THEN
          igpp=0
          IGWP=0
          NSCOUGW=0        
          IF((NTCYSGW.LT.ITHS).AND.(NSPOOLGW.GT.0)) THEN
            NTCYSGW=NTCYSGW+((ITHS-NTCYSGW)/NSPOOLGW)*NSPOOLGW
            IF(NTCYSGW.LT.ITHS) NTCYSGW=NTCYSGW+NSPOOLGW
            NDSETSW=(NTCYFGW-NTCYSGW)/NSPOOLGW
            ENDIF
          WRITE(16,1056)
 1056     FORMAT(//,' NEW GLOBAL WIND & pressure FILEs WILL BE STARTED')
          ENDIF

        IF(NOUTGW.EQ.-2) THEN
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
        IF(NOUTGW.EQ.-1) THEN
          open(73,file=dirname//'/'//'fort.73')
          write(73,3220) rundes,runid,agird
          write(73,3645) ndsetsw,np,dtdp*nspoolgw,nspoolgw,1
          igpp=2
          OPEN(74,FILE=DIRNAME//'/'//'fort.74')
          WRITE(74,3220) RUNDES,RUNID,AGRID
          WRITE(74,3645) NDSETSW,NP,DTDP*NSPOOLGW,NSPOOLGW,2
          IGWP=2
          ENDIF
        IF(NOUTGW.EQ.1) THEN
          OPEN(73,FILE=DIRNAME//'/'//'fort.73')
          do i=1,igpp           !I DON'T KNOW OF A PRACTICAL WAY TO CHANGE NDSETSW
            read(73,1050)
            enddo          
                OPEN(74,FILE=DIRNAME//'/'//'fort.74')
          DO I=1,IGWP           !I DON'T KNOW OF A PRACTICAL WAY TO CHANGE NDSETSW
            READ(74,1050)
            ENDDO
          ENDIF
        IF(NOUTGW.EQ.2) THEN
          open(73,file=dirname//'/'//'fort.73',
     &          access='direct',recl=nbyte)
          write(73,REC=itempstp+1) ndsetsw  ! ALLOW ADDITIONAL OUTPUT DATA TO BE WRITTEN
          close(73)                         ! DO THIS TO FLUSH THE WRITE BUFFER
          open(73,file=dirname//'/'//'fort.73',
     &         access='direct',recl=nbyte)
          OPEN(74,FILE=DIRNAME//'/'//'fort.74',
     &          ACCESS='DIRECT',RECL=NBYTE)
          WRITE(74,REC=ITEMPSTP+1) NDSETSW  ! ALLOW ADDITIONAL OUTPUT DATA TO BE WRITTEN
          CLOSE(74)                         ! DO THIS TO FLUSH THE WRITE BUFFER
          OPEN(74,FILE=DIRNAME//'/'//'fort.74',
     &         ACCESS='DIRECT',RECL=NBYTE)
          ENDIF
#ifdef SWAN
Casey 101118: Added the output of radiation stress gradients.
        IF(NOUTGW.LT.0) THEN
          IGRadS=0
          ENDIF
        IF(NOUTGW.EQ.-2) THEN
          OPEN(716,FILE=TRIM(LOCALDIR)//'/'//'rads.64',
     &          ACCESS='DIRECT',RECL=NBYTE)
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
          CLOSE(164)                    ! DO THIS TO FLUSH THE WRITE BUFFER
          OPEN(164,FILE=TRIM(LOCALDIR)//'/'//'rads.64',
     &         ACCESS='DIRECT',RECL=NBYTE)
          ENDIF
        IF(NOUTGW.EQ.-1) THEN
          OPEN(164,FILE=DIRNAME//'/'//'rads.64')
          WRITE(164,3220) RUNDES,RUNID,AGRID
          WRITE(164,3645) NDSETSW,NP,DTDP*NSPOOLGW,NSPOOLGW,2
          IGRadS=2
          ENDIF
        IF(NOUTGW.EQ.1) THEN
          OPEN(164,FILE=DIRNAME//'/'//'rads.64')
          DO I=1,IGRadS           !I DON'T KNOW OF A PRACTICAL WAY TO CHANGE NDSETSW
            READ(164,1050)
            ENDDO
          ENDIF
        IF(NOUTGW.EQ.2) THEN
          OPEN(164,FILE=DIRNAME//'/'//'rads.64',
     &          ACCESS='DIRECT',RECL=NBYTE)
          WRITE(164,REC=ITEMPSTP+1) NDSETSW  ! ALLOW ADDITIONAL OUTPUT DATA TO BE WRITTEN
          CLOSE(164)                         ! DO THIS TO FLUSH THE WRITE BUFFER
          OPEN(164,FILE=DIRNAME//'/'//'rads.64',
     &         ACCESS='DIRECT',RECL=NBYTE)
          ENDIF
#endif


C...
C......HOT START INFORMATION FOR HARMONIC ANALYSIS
C...
        IF(IHARIND.EQ.1) THEN
          IHABEG=ITHAS+NHAINC
C...
C........IF HARMONIC ANALYSIS HAS NOT BEGUN, COLD START THE HARMONIC ANALYSIS
C...
          IF(ITHS.LT.IHABEG) THEN
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
              END DO
             ENDIF   !   charmv

           ENDIF

C...
C........IF HARMONIC ANALYSIS HAS ALREADY BEGUN, READ IN HOT START
C........HARMONIC ANALYSIS, MEAN AND SQUARE INFO
C...
          IF(ITHS.GT.ITHAS) THEN
            IHOTSTP=IHOTSTP+1
            READ(IHOT,REC=IHOTSTP) ICHA
            ENDIF
          IF(ITHS.GE.IHABEG) THEN
            CALL HAHOTS(NSTAE,NSTAV,NP,NHASE,NHASV,NHAGE,NHAGV,
     &                  NSCREEN,IHOTSTP,IHOT,MYPROC)
            IF(NHASE.EQ.1) CALL HAHOTSES(NSTAE,IHOTSTP,IHOT)
            IF(NHASV.EQ.1) CALL HAHOTSVS(NSTAV,IHOTSTP,IHOT)
            IF(NHAGE.EQ.1) CALL HAHOTSEG(NP,IHOTSTP,IHOT)
            IF(NHAGV.EQ.1) CALL HAHOTSVG(NP,IHOTSTP,IHOT)
          ENDIF

C..Read in Means and Squares

        if( CHARMV) then
          IF((FMV.NE.0.).AND.(ITHS.GT.ITMV)) THEN
            IHOTSTP=IHOTSTP+1
            READ(IHOT,REC=IHOTSTP) NTSTEPS
            IF(NHAGE.EQ.1) THEN
              DO I=1,NP
                READ(IHOT,REC=IHOTSTP+1) ELAV(I)
                READ(IHOT,REC=IHOTSTP+2) ELVA(I)
                IHOTSTP=IHOTSTP+2
                ENDDO
              ENDIF
            IF(NHAGV.EQ.1) THEN
              DO I=1,NP
                READ(IHOT,REC=IHOTSTP+1) XVELAV(I)
                READ(IHOT,REC=IHOTSTP+2) YVELAV(I)
                READ(IHOT,REC=IHOTSTP+3) XVELVA(I)
                READ(IHOT,REC=IHOTSTP+4) YVELVA(I)
                IHOTSTP=IHOTSTP+4
                ENDDO
              ENDIF
            ENDIF
       endif   !  charmv


          ENDIF


       if(C3DVS) then
C3DVS..UNCOMMENT THE FOLLOWING LINES TO RUN THE CODE IN 3D VS MODE.
C3DVS..COMMENT OUT THE FOLLOWING LINES TO RUN THE CODE IN 2DDI OR 3D DSS MODE.
C3DVS..
c     ! CALL VSSTUP(DT,STATIM,NBYTE,RUNDES,RUNID,AGRID,NT)
C3DVS..
C3DVS..END OF 3D VS STATEMENTS (MORE FOLLOW BELOW)
C3DVS..
       end if



       if(C3DDSS)  then
C3DDSS.UNCOMMENT THE FOLLOWING LINES TO RUN THE CODE IN 3D DSS MODE
C3DDSS.COMMENT OUT THE FOLLOWING LINES TO RUN THE CODE IN 2DDI OR 3D VS MODE.
C3DDSS.
c     ! CALL DSSSTUP(DT,STATIM,NBYTE,RUNDES,RUNID,AGRID,NT)
C3DDSS.
C3DDSS.END OF 3D DSS STATEMENTS (MORE FOLLOW BELOW)
C3DDSS.
       end if

      CLOSE(IHOT)
C
 1112 FORMAT(/,1X,79('_'))
 1197 FORMAT(/,1X,'THE E29 MET GRID INTERPOLATING FACTORS ARE ',
     &                'BEING COMPUTED ')
 1198     FORMAT(1X,'FINISHED COMPUTING E29 INTERPOLATING FACTORS',/)
 3220 FORMAT(1X,A32,2X,A24,2X,A24)
 3645 FORMAT(1X,I10,1X,I10,1X,E15.7,1X,I5,1X,I5)
C

C.....Adjust Morphology timestep counter (defunct)
      

      RETURN
      END
