C***********************************************************************
C
C     SUBROUTINE  CREATE_EDGE_DATA.F
C
C     This program takes the original ADCIRC data sets given in files
C     'fort.14' and 'fort.15' and generates edge based data structures
C
C     Written by: Srinivas Chippada
C
C     Turned into a subroutine by Clint Dawson, May 2002
C
C     Modifications for different boundary types and notational changes
C     made by Ethan Kubatko, March 2005
C
C     Modifications for parallel runs
C     made by Shintaro Bunya, Aug 2005 
C
C     Bug fix
C     made by Shintaro Bunya, Aug 26, 2005 
C
C     Modifications to skip an edge on a barriar overlapping on 
C     the external boundary.  Boundary edge creation in a group 
C     used to be aborted when such an edge is found.
C     made by Shintaro Bunya, Sep  1, 2005 
C
C     Speed up in generating edge pairs
C     made by Shintaro Bunya, Feb 26, 2006
C
C***********************************************************************

      SUBROUTINE CREATE_EDGE_DATA()

C.....Use appropriate modules

      USE GLOBAL
      USE DG     
      use sizes, only:  myproc
#ifdef CMPI
      USE MESSENGER_ELEM
#endif
      IMPLICIT NONE
      
C.....Declare local variables

      INTEGER IED, JED, JJED, IEL, IEL1, IEL2, JJ1, JJ2
Csb--
      INTEGER JEL, JJEL,ED_ID,L,GED,I,J,K
C--
      INTEGER I1, I2, I3, II, LED(2,3)
      INTEGER NED1, NED2, NED3, NED4, NED5, NED6
      INTEGER NERR, NN1, NN2
      INTEGER IC1, IC2, IC3
      INTEGER IERROR
      
      REAL(SZ) mX, mY, diff
      INTEGER FLAG1, FLAG2

C.....Compute maximum number of edges

      MNED = 3*MNE

C.....Allocate the edge data arrays

      CALL ALLOC_EDGES1()

C.....Generate the edge connectivity

      DO J = 1,MNE
        EDFLG(1,J) = 0
        EDFLG(2,J) = 0
        EDFLG(3,J) = 0
        NELED(1,J) = 0
        NELED(2,J) = 0
        NELED(3,J) = 0
      ENDDO


      NEDNO(:,:) = 0
      NEDEL(:,:) = 0

      NEDGES = 0

      IF(MYPROC.EQ.0) THEN
         PRINT *,'CREATING EDGE PAIRS...'
      ENDIF

      DO 30 IEL = 1,NE
        N1 = NM(IEL,1)
        N2 = NM(IEL,2)
        N3 = NM(IEL,3)
        LED(1,1) = N2
        LED(2,1) = N3
        LED(1,2) = N3
        LED(2,2) = N1
        LED(1,3) = N1
        LED(2,3) = N2

        DO 10 IED = 1,3
          IF(EDFLG(IED,IEL).EQ.1) GOTO 10
           
          I1 = LED(1,IED)
          I2 = LED(2,IED)

          NEDGES = NEDGES + 1
          ED_ID = NEDGES
          NELED(IED,IEL) = ED_ID
          NEDNO(1,ED_ID) = I1
          NEDNO(2,ED_ID) = I2
          NEDEL(1,ED_ID) = IEL
          EDFLG(IED,IEL) = 1

          DO 15 JJEL = 1,NNDEL(I1)
            JEL = NDEL(I1,JJEL)
            IF(JEL.EQ.IEL) GOTO 15
            DO JED = 1,3
              J1 = NM(JEL,MOD(JED+0,3)+1)
              J2 = NM(JEL,MOD(JED+1,3)+1)
              IF ( ((J1.EQ.I1).AND.(J2.EQ.I2)).OR.
     &             ((J1.EQ.I2).AND.(J2.EQ.I1)) ) THEN

                IF(EDFLG(JED,JEL).EQ.1) THEN
                  PRINT *,'POSSIBLE DUPLICATE ELEMENT'
                  PRINT *,'MYPROC=',MYPROC
                  PRINT *,'EL=',JEL,', J1=',J1,', J2=',J2
                  PRINT *,'  (CREATE_EDGE_DATA.F)'
                  PRINT *,'EXECUTION WILL BE TERMINATED'
                  PRINT *,'! CHECK THE GRID CAREFULLY !'
                  PRINT*,'--------------------------------------'
                  STOP
                ENDIF
           
                NELED(JED,JEL) = ED_ID
                NEDEL(2,ED_ID) = JEL
                EDFLG(JED,JEL) = 1
                GOTO 10
              ENDIF
            ENDDO
 15       CONTINUE
 10     CONTINUE
 30   CONTINUE

      DEALLOCATE(EDFLG)

      IF(MYPROC.EQ.0) THEN
         PRINT *,'DONE'
         PRINT *,''
      ENDIF

C.....Print out the edge-to-node and edge-to-element connectivity

      WRITE(17,*)NEDGES
      DO IED = 1,NEDGES
        WRITE(17,*)IED,NEDNO(1,IED),NEDNO(2,IED),NEDEL(1,IED),
     &             NEDEL(2,IED)
      ENDDO
        
C.....Print out the element-to-edge connectivity

      DO IEL = 1,NE
        WRITE(17,*)IEL, NELED(1,IEL), NELED(2,IEL), NELED(3,IEL)
      ENDDO

C.....An index to keep track of the edges

      DO I = 1,NEDGES
        NCOUNT(I) = -1
      ENDDO
      
C.....Zero out internal edge counter and array

      NIEDS = 0
      NIEDN = 0
      
C.....Zero out land (no-normal flow) edge counter and array

      NLEDS = 0
      NLEDN = 0
      
C.....Zero out elevation specified open edge counter and array

      NEEDS = 0
      NEEDN = 0
      
C.....Zero out flow specified open edge counter and array

      NFEDS = 0
      NFEDN = 0
      
C.....Zero out the radiation edge counter and array

      NREDS = 0
      NREDN = 0
      
C.....Zero out internal and external barrier counters and arrays

      NIBEDS = 0
      NEBEDS = 0

      NIBSEG = 0
      NEBSEG = 0
      
      NIBEDN  = 0
      NEBEDN  = 0
      NIBSEGN = 0
      NEBSEGN = 0

C.....Younghun: Zero out channel edge counter and array
      NCEDS = 0
      NCEDN = 0
      NCEDNO = 0
      NCEDED = 0
      NEDCEDN = 0

C.....Find node pairs that are not edges of elements (eg. an internal
C.....barrier against a land boundary)

      NOT_AN_EDGE = 0
      WEIR_BUDDY_NODE = 0
      DO 10130 II = 1,NVEL
        IF( (LBCODEI(II).EQ.4).OR.(LBCODEI(II).EQ.24).OR.
     &      (LBCODEI(II).EQ.5).OR.(LBCODEI(II).EQ.25) ) THEN
          J1 = NBV(II)      ! GLOBAL NODE NUMBER ON BACK SIDE OF BARRIER
          J2 = IBCONN(II)   ! GLOBAL NODE NUMBER ON FRONT SIDE OF BARRIER
          DO K = 1,NBOU
            IF ( (SEGTYPE(K).EQ.4 ).OR.(SEGTYPE(K).EQ.24) ) THEN
              IF (WEIR_BUDDY_NODE(J1,1).EQ.0) THEN
                WEIR_BUDDY_NODE(J1,1) = J2
              ELSE
                WEIR_BUDDY_NODE(J1,2) = J2
              ENDIF
            ENDIF
            IF ( (SEGTYPE(K).EQ.0 ).OR.(SEGTYPE(K).EQ.3 ).OR.
     &           (SEGTYPE(K).EQ.2 ).OR.(SEGTYPE(K).EQ.12).OR.
     &           (SEGTYPE(K).EQ.13).OR.(SEGTYPE(K).EQ.20).OR.
     &           (SEGTYPE(K).EQ.22).OR.(SEGTYPE(K).EQ.23) ) THEN
              DO IED = 1,NVELL(K)-1
                N1 = NBVV(K,IED)
                N2 = NBVV(K,IED+1)
                IF ((N1.EQ.J1).OR.(N1.EQ.J2)) THEN
                  IF ((N2.EQ.J1).OR.(N2.EQ.J2)) THEN
                    NOT_AN_EDGE(N1) = 1
                    NOT_AN_EDGE(N2) = 1
                  ENDIF
                ENDIF
              ENDDO
            ELSEIF ( (SEGTYPE(K).EQ.1 ).OR.(SEGTYPE(K).EQ.11).OR.
     &               (SEGTYPE(K).EQ.21) ) THEN
              DO IED = 1,NVELL(K)-1
                N1 = NBVV(K,IED)
                N2 = NBVV(K,IED+1)
                IF ((N1.EQ.J1).OR.(N1.EQ.J2)) THEN
                  IF ((N2.EQ.J1).OR.(N2.EQ.J2)) THEN
                    NOT_AN_EDGE(N1) = 1
                    NOT_AN_EDGE(N2) = 1
                  ENDIF
                ENDIF
              ENDDO
            ENDIF
          ENDDO
        ENDIF
10130 CONTINUE

C.....Find the interior edges

      DO I = 1,NEDGES
        IEL1 = NEDEL(1,I)
        IEL2 = NEDEL(2,I)
        IF((IEL1.NE.0).AND.(IEL2.NE.0))THEN
          NIEDS = NIEDS + 1
          NIEDN(NIEDS) = I
          NCOUNT(I) = 0
        ENDIF
      ENDDO

C.....Find elevation specified boundary edges

      DO 10131 K = 1,NOPE
        DO IED = 1,NVDLL(K)-1
          N1 = NBDV(K,IED)
          N2 = NBDV(K,IED+1)
          IERROR = 0
          DO JED = 1,NEDGES
            J1 = NEDNO(1,JED)
            J2 = NEDNO(2,JED)
            IF ( (J1.EQ.N1).OR.(J1.EQ.N2) ) THEN
              IF ( (J2.EQ.N1).OR.(J2.EQ.N2) ) THEN
                NEEDS = NEEDS + 1
                NEEDN(NEEDS) = JED
                NCOUNT(JED) = 4
                IERROR = 1
              ENDIF
            ENDIF
          ENDDO
Csb-PDG1
#ifdef CMPI 
          IF(IERROR.EQ.0) THEN
            WRITE(16,*)
     &          'ERROR IN PROCESSING OPEN OCEAN BOUNDARY CONDITIONS'
          ENDIF
#else
          IF(IERROR.EQ.0) THEN
            STOP 'ERROR IN PROCESSING OPEN OCEAN BOUNDARY CONDITIONS'
          ENDIF
#endif
C--
        ENDDO
10131 CONTINUE

C.....Find flux specified boundary edges

      JNMM = 0
      ONE_OR_TWO = 1
      DO 10132 K = 1,NBOU
Csb-
        DO 10133 IED = 1,NVELL(K)-1
C--
          N1 = NBVV(K,IED)
          N2 = NBVV(K,IED+1)

          IERROR = 0
          DO JED = 1,NEDGES
            J1 = NEDNO(1,JED)
            J2 = NEDNO(2,JED)
            IF ((NOT_AN_EDGE(N1).EQ.1).AND.(NOT_AN_EDGE(N2).EQ.1)) THEN
C              PRINT*,'NODES ',N1,' AND ',N2,' MAKE UP A BOUNDARY',
C     &                  'SEGMENT THAT IS NOT AN EDGE TO AN ELEMENT; ',
Csb-
C     &                  'THIS EDGE IS SKIPPED'
              GOTO 10133
C--
            ENDIF
            IF ( (N1.EQ.J1).OR.(N1.EQ.J2) ) THEN
              IF ( (N2.EQ.J1).OR.(N2.EQ.J2) ) THEN

                IERROR = 1
                NCOUNT(JED) = 1

C.....Determine the different boundary types

                IF ( (SEGTYPE(K).EQ.0 ).OR.(SEGTYPE(K).EQ.10).OR.
     &               (SEGTYPE(K).EQ.20) ) THEN
                  NLEDS = NLEDS + 1
                  NLEDN(NLEDS) = JED
                ENDIF

                IF ( (SEGTYPE(K).EQ.1 ).OR.(SEGTYPE(K).EQ.11).OR.
     &               (SEGTYPE(K).EQ.21) ) THEN
                  NLEDS = NLEDS + 1
                  NLEDN(NLEDS) = JED
                ENDIF

                IF ( (SEGTYPE(K).EQ.2 ).OR.(SEGTYPE(K).EQ.12).OR.
     &               (SEGTYPE(K).EQ.22) ) THEN
                  NFEDS = NFEDS + 1
                  NFEDN(NFEDS) = JED
                ENDIF
                
                IF ( (SEGTYPE(K).EQ.3 ).OR.(SEGTYPE(K).EQ.13).OR.
     &               (SEGTYPE(K).EQ.23) ) THEN
                  NEBSEG = NEBSEG + 1
                  NEBEDS = NEBEDS + 1
                  NEBEDN(NEBEDS) = JED
                  NEBSEGN(NEBSEG) = JED
                ENDIF
                
                IF ( (SEGTYPE(K).EQ.4 ).OR.(SEGTYPE(K).EQ.14).OR.
     &               (SEGTYPE(K).EQ.24) ) THEN
                  NIBSEG = NIBSEG + 1
                  NIBEDS = NIBEDS + 1
                  NIBEDN(NIBEDS) = JED
                  NIBSEGN(1,NIBSEG) = JED

C.....Find the edge on the opposite side of the barrier

                  NN1 = BACKNODES(1,NIBSEG)
                  NN2 = BACKNODES(2,NIBSEG)
c                  if (myproc.eq.24) then
c                     write(200+myproc,*) jed,nn1,nn2
c                  endif
                  DO JJED = 1,NEDGES
                    JJ1 = NEDNO(1,JJED)
                    JJ2 = NEDNO(2,JJED)
                    IF ( (NN1.EQ.JJ1).OR.(NN1.EQ.JJ2) ) THEN
                      IF ( (NN2.EQ.JJ1).OR.(NN2.EQ.JJ2) ) THEN
                        NIBEDS = NIBEDS + 1
                        NIBEDN(NIBEDS) = JJED
                        NIBSEGN(2,NIBSEG) = JJED
                        NCOUNT(JJED) = 1
                        ONE_OR_TWO(N1) = 2
                        ONE_OR_TWO(N2) = 2
                      ENDIF
                    ENDIF
                  ENDDO
c                  if (nibsegn(2,nibseg).eq.0) then
c                     write(200+myproc,*) 'error in create_edge_data'
c                     write(200+myproc,*) myproc,nibseg,jed,nn1,nn2
c                  endif
                ENDIF
                
                IF ( (SEGTYPE(K).EQ.30) ) THEN
                  NREDS = NREDS + 1
                  NREDN(NREDS) = JED
                ENDIF
                
                ! Younghun added
                IF (SEGTYPE(K).eq.6) THEN
                   NCEDS = NCEDS + 1
                   NCEDN(NCEDS) = JED
                   ! NEDCEDN(JED) = NCEDS
                   
                   NCEDNO(1,NCEDS) = N1
                   NCEDNO(2,NCEDS) = N2
                   
                   ! write(9998,*) K,JED,N1,N2
                   if (J1.eq.N1 .AND. J2.eq.N2) then
                      NCEDNO(3,NCEDS) = 0 ! Flag for nonflip NEDNO direction
                   elseif (J2.eq.N1 .AND. J1.eq.N2) then
                      NCEDNO(3,NCEDS) = 1 ! Flag for flip NEDNO direction
                   endif

                   ! HYDR_INFO(1,N1) = CBWIDTH(IED)
                   ! HYDR_INFO(2,N1) = (CTWIDTH(IED) - CBWIDTH(IED))/(2.D0*CDEPTH(IED))
                   
                   ! HYDR_INFO(1,N2) = CBWIDTH(IED+1)
                   ! HYDR_INFO(2,N2) = (CTWIDTH(IED+1) - CBWIDTH(IED+1))/(2.D0*CDEPTH(IED+1))
                   
                   ! print*, HYDR_INFO(1,N1), HYDR_INFO(2,N1)
                   ! print*, HYDR_INFO(1,N2), HYDR_INFO(2,N2)

                   CDEPTH(1,NCEDS)  = CDEPTHR(K,IED)
                   CDEPTH(2,NCEDS)  = CDEPTHR(K,IED+1)
                   CTWIDTH(1,NCEDS) = CTWIDTHR(K,IED)
                   CTWIDTH(2,NCEDS) = CTWIDTHR(K,IED+1)
                   CBWIDTH(1,NCEDS) = CBWIDTHR(K,IED)
                   CBWIDTH(2,NCEDS) = CBWIDTHR(K,IED+1)
                   ! Compute the side slope
                   do J=1,2
                      if (CTWIDTH(J,NCEDS).eq.CBWIDTH(J,NCEDS)) then
                         CSSLOPE(J,NCEDS) = 1e8
                      else
                         CSSLOPE(J,NCEDS) = 2*CDEPTH(J,NCEDS)/
     &                    (CTWIDTH(J,NCEDS) - CBWIDTH(J,NCEDS))
                      endif
                   enddo
                   CEDSLOPE(NCEDS) = abs( (CDEPTH(1,NCEDS)-DP(N1)) -
     &              (CDEPTH(2,NCEDS)-DP(N2)) )/sqrt((X(N1)-X(N2))**2.D0
     &                                      + (Y(N1)-Y(N2))**2.D0)
                   ! print *, CEDSLOPE(NCEDS)
                ENDIF
                
                
              ENDIF
            ENDIF
Csb-
          ENDDO
C--
Csb-PDG1
#ifdef CMPI 
          IF(IERROR.EQ.0) THEN
            WRITE(16,*) 'NODE PAIR (',N1,',',N2,') IS NOT AN EDGE.'
            WRITE(16,*) 'ERROR IN PROCESSING LAND SEGMENT'
            WRITE(16,*) ''
          ENDIF
#else
           IF(IERROR.EQ.0) then
             WRITE(*,*) 'NODE PAIR (',N1,',',N2,') IS NOT AN EDGE.'
             WRITE(*,*) ''
             WRITE(*,*) 'ERROR IN PROCESSING LAND SEGMENT'
             WRITE(*,*) ''
C             STOP 'ERROR IN PROCESSING LAND SEGMENT'
           endif
#endif
C--
10133   CONTINUE
10132 CONTINUE

C.....2021-05-12 Younghun
#ifdef CMPI
      CALL MSG_TABLE_CHANNEL()
      CALL MESSAGE_START_CHANNEL()
#endif

C.....Younghun added
      ALLOCATE (NCEDIN(MNED,NCEDS+1), NCEDEX(MNED,NCEDS+1))
      NCEDIN = 0
      NCEDEX = 0
      
      DO I = 1,NCEDS
         N1 = NCEDNO(1,I)
         N2 = NCEDNO(2,I)
         
         DO J = 1,NCEDS
            IF (NCEDNO(2,J).eq.N1) THEN
               NCEDED(1,I) = NCEDN(J)
            
               NCEDIN(1,I) = NCEDIN(1,I) + 1
               K = NCEDIN(1,I)
               NCEDIN(K+1,I) = NCEDN(J)
            ENDIF
            IF (NCEDNO(1,J).eq.N2) THEN
               NCEDED(2,I) = NCEDN(J)
               
               NCEDEX(1,I) = NCEDEX(1,I) + 1
               K = NCEDEX(1,I)
               NCEDEX(K+1,I) = NCEDN(J)
               
            ENDIF
         ENDDO
         
         GED = NCEDN(I)
         if (CBWIDTH(1,I).ne.0 .and. CBWIDTH(2,I).ne.0) then
            NEDCEDN(GED) = I
         endif
         ! write(*,*) NCEDED(1,I), NCEDED(2,I)
         ! WRITE(*,*) x(N1), y(N1)
         ! WRITE(*,*) x(N2), y(N2)
         ! print*,NCEDIN(1,I), NCEDEX(1,I)
         
      ENDDO
      
C.....Find the interior edges (exclude channel-edges)
      NIEDS = 0
      NIEDN = 0
      DO I = 1,NEDGES
        IEL1 = NEDEL(1,I)
        IEL2 = NEDEL(2,I)
        DO J = 1, NCEDS
           IF (NCEDN(J).EQ.I) GOTO 40
        ENDDO
        IF((IEL1.NE.0).AND.(IEL2.NE.0))THEN
          NIEDS = NIEDS + 1
          NIEDN(NIEDS) = I
          NCOUNT(I) = 0
        ENDIF
 40   ENDDO

C.....Check the order of the nodes assigned to an edge - important in
C.....the calculation of the unit normal

      DO I = 1,NEDGES
        N1 = NEDNO(1,I)
        N2 = NEDNO(2,I)
        IEL = NEDEL(1,I)
        IF((N1.EQ.NM(IEL,2)).AND.(N2.EQ.NM(IEL,1))) THEN
          WRITE(6,*)'THE ORDER OF NODES ASSIGNED IS WRONG FOR EDGE',I
          WRITE(16,*)'THE ORDER OF NODES ASSIGNED IS WRONG FOR EDGE',I
          STOP
        ENDIF
        IF((N1.EQ.NM(IEL,3)).AND.(N2.EQ.NM(IEL,2)))THEN
          WRITE(6,*)'THE ORDER OF NODES ASSIGNED IS WRONG FOR EDGE',I
          WRITE(16,*)'THE ORDER OF NODES ASSIGNED IS WRONG FOR EDGE',I
          STOP
        ENDIF
        IF((N1.EQ.NM(IEL,1)).AND.(N2.EQ.NM(IEL,3)))THEN
          WRITE(6,*)'THE ORDER OF NODES ASSIGNED IS WRONG FOR EDGE',I
          WRITE(16,*)'THE ORDER OF NODES ASSIGNED IS WRONG FOR EDGE',I
          STOP
        ENDIF
      ENDDO
      
C.....Check for missing edges

Csb-PDG1 modified
      NERR = NEDGES - (NIEDS + NLEDS + NEEDS + NFEDS + NIBEDS + NEBEDS +
     &                 NREDS + NCEDS)

      IF(MYPROC.EQ.0) THEN
         WRITE(6,*) '  '
         WRITE(6,*) 'TOTAL NO. OF EDGES = ', NEDGES
         WRITE(6,*) '  '
         WRITE(6,*) 'NO. OF INTERNAL (NON-BOUNDARY) EDGES = ', NIEDS
         WRITE(6,*) 'NO. OF NO-NORMAL FLOW EDGES = ', NLEDS
         WRITE(6,*) 'NO. OF NON-ZERO NORMAL FLOW EDGES = ', NFEDS
         WRITE(6,*) 'NO. OF ELEVATION SPECIFIED EDGES = ', NEEDS
         WRITE(6,*) 'NO. OF EXTERNAL BARRIER EDGES = ', NEBEDS
         WRITE(6,*) 'NO. OF INTERNAL BARRIER EDGES = ', NIBEDS
         WRITE(6,*) 'NO. OF RADIATION EDGES = ', NREDS
         WRITE(6,*) 'NO. OF OPEN CHANNEL EDGES = ', NCEDS
         WRITE(6,*) 
     &        '-----------------------------------------------------'
         WRITE(6,*) 'NO. OF MISSING EDGES = ',NERR
         WRITE(6,*)  ''
      ENDIF
      WRITE(16,*) '  '
      WRITE(16,*) 'TOTAL NO. OF EDGES = ', NEDGES
      WRITE(16,*) '  '
      WRITE(16,*) 'NO. OF INTERNAL (NON-BOUNDARY) EDGES = ', NIEDS
      WRITE(16,*) 'NO. OF NO-NORMAL FLOW EDGES = ', NLEDS
      WRITE(16,*) 'NO. OF NON-ZERO NORMAL FLOW EDGES = ', NFEDS
      WRITE(16,*) 'NO. OF ELEVATION SPECIFIED EDGES = ', NEEDS
      WRITE(16,*) 'NO. OF EXTERNAL BARRIER EDGES = ', NEBEDS
      WRITE(16,*) 'NO. OF INTERNAL BARRIER EDGES = ', NIBEDS
      WRITE(16,*) 'NO. OF RADIATION EDGES = ', NREDS
      WRITE(16,*) 'NO. OF OPEN CHANNEL EDGES = ', NCEDS
      WRITE(16,*) 
     &        '-----------------------------------------------------'
      WRITE(16,*) 'NO. OF MISSING EDGES = ',NERR
      WRITE(16,*)  ''

      DO I = 1,NEDGES
        IF (NCOUNT(I).LT.0) THEN
          N1 = NEDNO(1,I)
          N2 = NEDNO(2,I)
Csb-PDG1
#ifdef CMPI
C          WRITE(6,*)' '
C          WRITE(6,*)'EDGE ',I,' IS MADE UP OF NODES ',N1,' AND ',N2
C          WRITE(6,*)'EDGE ',I, ' IS NEITHER AN INTERNAL NOR A BOUNDARY',
C     &       'MAKE SURE IF THIS IS DUE TO THE DOMAIN DECOMPOSITION'
#else
          WRITE(16,*)' '
          WRITE(16,*)'EDGE',I,',MADE UP OF NODES',N1,'AND',N2,', IS NOT'
          WRITE(16,*)'AN INTERNAL(NON-BOUNDARY) EDGE OR A BOUNDARY EDGE'
          WRITE(16,*)'ASSUMING EDGE',I,'IS A NO-NORMAL FLOW EDGE !!!'
C          STOP
#endif
C--
          NLEDS = NLEDS + 1
          NLEDN(NLEDS) = I
        ENDIF
      ENDDO
      
C.....Add internal barrier edges to land edge table for wet-dry
C.....post-processing

      DO I = 1,NIBEDS
        NLEDN(NLEDS+I) = NIBEDN(I)
      ENDDO
      
C.....Print out the interior edges

      WRITE(17,*) NIEDS,'       ! NUMBER OF INTERNAL EDGES'
      DO I = 1,NIEDS
        WRITE(17,*) I,NIEDN(I),nedno(1,niedn(i)),nedno(2,niedn(i))
      ENDDO

C.....Prin out land edges info.

      WRITE(17,*) NLEDS,'       ! NUMBER OF NO-NORMAL FLOW EDGES'
      IF (NLEDS.GT.0) THEN
        DO I = 1,NLEDS
          WRITE(17,*) I, NLEDN(I),nedno(1,nledn(i)),
     $          nedno(2,nledn(i))
        ENDDO
      ENDIF
      
C.....Print out elevation specified edge info.

      WRITE(17,*) NEEDS,'       ! NUMBER OF ELEVATION SPECIFIED EDGES'
      IF (NEEDS.GT.0) THEN
        DO I = 1,NEEDS
          WRITE(17,*) I, NEEDN(I)
        ENDDO
      ENDIF

C.....Print out non-zero flow edges info.

      WRITE(17,*) NFEDS,'       ! NUMBER OF FLOW SPECIFIED EDGES'
      IF (NFEDS.GT.0) THEN
        DO I = 1,NFEDS
          WRITE(17,*) I, NFEDN(I)
        ENDDO
      ENDIF
      
C.....Print out external barrier edge info.

      WRITE(17,*) NEBEDS,'       ! NUMBER OF EXTERNAL BARRIER EDGES'
      IF (NEEDS.GT.0) THEN
        DO I = 1,NEBEDS
          WRITE(17,*) I, NEBEDN(I),nedno(1,nebedn(i)),
     $          nedno(2,nebedn(i))
        ENDDO
      ENDIF
      
C.....Print out internal barrier edge info.

      WRITE(17,*) NIBEDS,'       ! NUMBER OF INTERNAL BARRIER EDGES'
      IF (NIBEDS.GT.0) THEN
        DO I = 1,NIBEDS
          WRITE(17,*) I, NIBEDN(I),nedno(1,nibedn(i)),
     $          nedno(2,nibedn(i))
C          PRINT*,'INTERNAL BARRIER EDGE =',I
C          PRINT*,'IS MADE UP OF NODES',NEDNO(1,NIBEDN(I)),'AND',
C     &                                 NEDNO(2,NIBEDN(I))
        ENDDO
      ENDIF

C.....Print out radiation edges info.

      WRITE(17,*) NREDS, '      ! NUMBER OF RADIATION EDGES'
      IF (NREDS.GT.0) THEN
        DO I = 1,NREDS
          WRITE(17,*) I, NREDN(I)
        ENDDO
      ENDIF
      
C.....Construct global edge to local edge (1,2, or 3) table

      NEDSD(:,:) = 0.D0
      DO I = 1,MNED
        DO K = 1,2
          IF (NEDEL(K,I).NE.0) THEN
            N1 = NELED(1,NEDEL(K,I))
            N2 = NELED(2,NEDEL(K,I))
            N3 = NELED(3,NEDEL(K,I))
            IF (N1.EQ.I) NEDSD(K,I) = 1
            IF (N2.EQ.I) NEDSD(K,I) = 2
            IF (N3.EQ.I) NEDSD(K,I) = 3
          ENDIF
        ENDDO
      ENDDO
      
C.....Younghun: Assign Edge type to identify open channel edges
      CALL CALC_TRISLOPE()
      
      DO I = 1,NEDGES
        ET(I) = -1
        flipNEDEL(I) = 0
        
        if (NEDCEDN(I).gt.0) then
           ! ET(I) = -2
        endif
        N1 = NEDNO(1,I)
        N2 = NEDNO(2,I)
        
        mX = (X(N1) + X(N2))/2
        mY = (Y(N1) + Y(N2))/2
        
        IEL1 = NEDEL(1,I)
        IEL2 = NEDEL(2,I)
        
        IF (IEL1 * IEL2.EQ.0) THEN
            CALL pointInTri(mX + S0X(IEL1), mY + S0Y(IEL1),IEL1,FLAG1)
            IF (FLAG1.EQ.-1) THEN
                ET(I) = 1
            ELSE
                ET(I) = 0
            ENDIF 
               
        ELSE
            CALL pointInTri(mX + S0X(IEL1), mY + S0Y(IEL1),IEL1,FLAG1)
            CALL pointInTri(mX + S0X(IEL2), mY + S0Y(IEL2),IEL2,FLAG2)
            
C......Younghun's classification
            
            IF (FLAG1.EQ.-1) THEN
               if (FLAG2.EQ.-1) then
                  ET(I) = -2
               else
                  ET(I) = -1
               endif
               
            ELSEIF (FLAG1.EQ.0) then
               if (FLAG2.EQ.-1) then
                  ET(I) = -1
                  flipNEDEL(I) = 1
               else
                  ET(I) = -1
               endif
               
            ELSEIF (FLAG1.EQ.1) then
               if (FLAG2.EQ.1) then
                  ET(I) = 0
               else
                  ET(I) = -1
                  flipNEDEL(I) = 1
               endif
            ENDIF
                
        ENDIF
        

        ! DO J = 1, NCEDS
           ! IF (NCEDN(J).eq.I .AND. ET(I).ne.-2) THEN
              ! ET(I) = -2
           ! ENDIF
        ! ENDDO 
        
                
        ! DO II = 1, NCEDS
           ! I1 = NCEDN(II)
           ! I2 = NCEDED(1,II)
           ! I3 = NCEDED(2,II)
           
           ! IF (I2*I3.ne.0) THEN
              ! IF (ET(I2).eq.ET(I3)) THEN
                 ! ET(I1) = ET(I2)
              ! ENDIF
           ! ENDIF
        ! ENDDO
        ! WRITE(99999,*) I, mx, my, ET(I)
        
      ENDDO
      
      ! DO I = 1, NCEDS
         ! GED = NCEDN(I)
         ! IF (ET(GED).eq.0 .AND. NEDCEDN(GED).eq.0) THEN
            ! ET(GED) = -1
         ! ENDIF
      ! ENDDO
      
      ! do ii = 1, nceds
         ! i1 = ncedn(ii)
         ! ! write(9999,*) et(i1)
         ! if (et(i1).eq.0) then
            ! et(i1) = -1
         ! endif
      ! enddo
      
      ! DO II = 1, NCEDS
         ! I1 = NCEDN(II)
         ! I2 = NCEDED(1,II)
         ! I3 = NCEDED(2,II)
         
         ! IF (I2*I3.ne.0) THEN
            ! IF (ET(I2).eq.ET(I3) .AND. ET(I2).lt.ET(I1)) THEN
               ! ET(I1) = ET(I2)
            ! ENDIF
         ! ENDIF
      ! ENDDO
      
      ! DO I = 1, NEDGES
         ! READ(164,*) ET(I)
      ! ENDDO
      
      
      ! DO I = 1, NEDGES
         ! N1 = NEDNO(1,I)
         ! N2 = NEDNO(2,I)
         
         ! mX = (X(N1) + X(N2))/2
         ! mY = (Y(N1) + Y(N2))/2
         ! WRITE(9993,*) I, mX, mY, ET(I)
      ! ENDDO
      
      ! DO I = 1, NEDGES
         ! N1 = NEDNO(1,I)
         ! N2 = NEDNO(2,I)
         
         ! mX = (X(N1) + X(N2))/2
         ! mY = (Y(N1) + Y(N2))/2
         ! WRITE(9996,*) I, mX, mY, flipNEDEL(I)
      ! ENDDO
      ! stop
      
C      READ(163,*) I
C      IF (I.ne.NCEDS) THEN
C         print *, '*********** Execution terminated ***********'
C         print *, 'The number of open-channels in fort.163 file does not match with mesh'
C         stop
C      ENDIF
C      
C      ALLOCATE (manningKWE1(NCEDS))
C      
C      DO J = 1,NCEDS
C         READ(163,*) manningKWE1(J)
C      ENDDO
      
      ! manningKWE1 = 0.15
      
      RETURN
      END SUBROUTINE CREATE_EDGE_DATA
