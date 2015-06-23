!***********************************************************************
!
!     SUBROUTINE  CREATE_EDGE_DATA.dg%F
!
!     This program takes the original ADCIRC data sets given in files
!     'fort.14' and 'fort.15' and generates edge based data structures
!
!     Written by: Srinivas Chippada
!
!     Turned into a subroutine by Clint Dawson, May 2002
!
!     Modifications for different boundary types and notational changes
!     made by Ethan Kubatko, March 2005
!
!     Modifications for parallel runs
!     made by Shintaro Bunya, Aug 2005 
!
!     Bug fix
!     made by Shintaro Bunya, Aug 26, 2005 
!
!     Modifications to skip an edge on a barriar overlapping on 
!     the external boundary.  Boundary edge creation in a group 
!     used to be aborted when such an edge is found.
!     made by Shintaro Bunya, Sep  1, 2005 
!
!     Speed up in generating edge pairs
!     made by Shintaro Bunya, Feb 26, 2006
!
!***********************************************************************

      SUBROUTINE CREATE_EDGE_DATA(s)

!.....Use appropriate modules

      USE GLOBAL
      USE DG     
      use sizes
#ifdef CMPI
      USE MESSENGER_ELEM
#endif
      IMPLICIT NONE
      type (sizes_type) :: s
      
!.....Declare local variables

      INTEGER IED, JED, JJED, IEL, IEL1, IEL2, JJ1, JJ2
!sb--
      INTEGER JEL, JJEL,ED_ID,L,GED,I,J,K
!--
      INTEGER I1, I2, I3, II, LED(2,3)
      INTEGER NED1, NED2, NED3, NED4, NED5, NED6
      INTEGER NERR, NN1, NN2
      INTEGER IC1, IC2, IC3
      INTEGER IERROR

!.....Compute maximum number of edges

      dg%MNED = 3*s%MNE

!.....Allocate the edge data arrays

      CALL ALLOC_EDGES1()

!.....Generate the edge connectivity

      DO J = 1,S%MNE
        EDFLG(1,J) = 0
        EDFLG(2,J) = 0
        EDFLG(3,J) = 0
        NELED(1,J) = 0
        NELED(2,J) = 0
        NELED(3,J) = 0
      ENDDO


      dg%NEDNO(:,:) = 0
      dg%NEDEL(:,:) = 0

      dg%NEDGES = 0
#ifdef CMPI
      IF(MYPROC.EQ.0) THEN
         PRINT *,'CREATING EDGE PAIRS...'
      ENDIF
#else
      PRINT *,'CREATING EDGE PAIRS...'
#endif

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

          dg%NEDGES = dg%NEDGES + 1
          ED_ID = dg%NEDGES
          NELED(IED,IEL) = ED_ID
          dg%NEDNO(1,ED_ID) = I1
          dg%NEDNO(2,ED_ID) = I2
          dg%NEDEL(1,ED_ID) = IEL
          EDFLG(IED,IEL) = 1

          DO 15 JJEL = 1,NNDEL(I1)
            JEL = NDEL(I1,JJEL)
            IF(JEL.EQ.IEL) GOTO 15
            DO JED = 1,3
              dg%J1 = NM(JEL,MOD(JED+0,3)+1)
              dg%J2 = NM(JEL,MOD(JED+1,3)+1)
              IF ( ((dg%J1.EQ.I1).AND.(dg%J2.EQ.I2)).OR.&
             ((dg%J1.EQ.I2).AND.(dg%J2.EQ.I1)) ) THEN

                IF(EDFLG(JED,JEL).EQ.1) THEN
                  PRINT *,'POSSIBLE DUPLICATE ELEMENT'
#ifdef CMPI
                  PRINT *,'MYPROC=',MYPROC
#endif
                  PRINT *,'dg%EL=',JEL,', dg%J1=',dg%J1,', dg%J2=',dg%J2
                  PRINT *,'  (CREATE_EDGE_DATA.dg%F)'
                  PRINT *,'EXECUTION WILL BE TERMINATED'
                  PRINT *,'! CHECK THE GRID CAREFULLY !'
                  PRINT*,'--------------------------------------'
                  STOP
                ENDIF
           
                NELED(JED,JEL) = ED_ID
                dg%NEDEL(2,ED_ID) = JEL
                EDFLG(JED,JEL) = 1
                GOTO 10
              ENDIF
            ENDDO
 15       CONTINUE
 10     CONTINUE
 30   CONTINUE

      DEALLOCATE(EDFLG)

#ifdef CMPI
      IF(MYPROC.EQ.0) THEN
         PRINT *,'DONE'
         PRINT *,''
      ENDIF
#else
         PRINT *,'DONE'
         PRINT *,''
#endif

!.....Print out the edge-to-node and edge-to-element connectivity

      WRITE(17,*)dg%NEDGES
      DO IED = 1,dg%NEDGES
        WRITE(17,*)IED,dg%NEDNO(1,IED),dg%NEDNO(2,IED),dg%NEDEL(1,IED),&
             dg%NEDEL(2,IED)
      ENDDO
        
!.....Print out the element-to-edge connectivity

      DO IEL = 1,NE
        WRITE(17,*)IEL, NELED(1,IEL), NELED(2,IEL), NELED(3,IEL)
      ENDDO

!.....An index to keep track of the edges

      DO I = 1,dg%NEDGES
        dg%NCOUNT(I) = -1
      ENDDO
      
!.....Zero out internal edge counter and array

      dg%NIEDS = 0
      dg%NIEDN = 0
      
!.....Zero out land (no-normal flow) edge counter and array

      dg%NLEDS = 0
      dg%NLEDN = 0
      
!.....Zero out elevation specified open edge counter and array

      dg%NEEDS = 0
      dg%NEEDN = 0
      
!.....Zero out flow specified open edge counter and array

      dg%NFEDS = 0
      dg%NFEDN = 0
      
!.....Zero out the radiation edge counter and array

      dg%NREDS = 0
      dg%NREDN = 0
      
!.....Zero out internal and external barrier counters and arrays

      dg%NIBEDS = 0
      dg%NEBEDS = 0

      dg%NIBSEG = 0
      dg%NEBSEG = 0
      
      dg%NIBEDN  = 0
      dg%NEBEDN  = 0
      dg%NIBSEGN = 0
      dg%NEBSEGN = 0

!.....Find node pairs that are not edges of elements (eg. an internal
!.....barrier against a land boundary)

      NOT_AN_EDGE = 0
      WEIR_BUDDY_NODE = 0
      DO 10130 II = 1,NVEL
        IF( (LBCODEI(II).EQ.4).OR.(LBCODEI(II).EQ.24).OR.&
      (LBCODEI(II).EQ.5).OR.(LBCODEI(II).EQ.25) ) THEN
          dg%J1 = NBV(II)      ! GLOBAL NODE NUMBER ON BACK SIDE OF BARRIER
          dg%J2 = IBCONN(II)   ! GLOBAL NODE NUMBER ON FRONT SIDE OF BARRIER
          DO K = 1,NBOU
            IF ( (SEGTYPE(K).EQ.4 ).OR.(SEGTYPE(K).EQ.24) ) THEN
              IF (WEIR_BUDDY_NODE(dg%J1,1).EQ.0) THEN
                WEIR_BUDDY_NODE(dg%J1,1) = dg%J2
              ELSE
                WEIR_BUDDY_NODE(dg%J1,2) = dg%J2
              ENDIF
            ENDIF
            IF ( (SEGTYPE(K).EQ.0 ).OR.(SEGTYPE(K).EQ.3 ).OR.&
           (SEGTYPE(K).EQ.2 ).OR.(SEGTYPE(K).EQ.12).OR.&
           (SEGTYPE(K).EQ.13).OR.(SEGTYPE(K).EQ.20).OR.&
           (SEGTYPE(K).EQ.22).OR.(SEGTYPE(K).EQ.23) ) THEN
              DO IED = 1,NVELL(K)-1
                N1 = NBVV(K,IED)
                N2 = NBVV(K,IED+1)
                IF ((N1.EQ.dg%J1).OR.(N1.EQ.dg%J2)) THEN
                  IF ((N2.EQ.dg%J1).OR.(N2.EQ.dg%J2)) THEN
                    NOT_AN_EDGE(N1) = 1
                    NOT_AN_EDGE(N2) = 1
                  ENDIF
                ENDIF
              ENDDO
            ELSEIF ( (SEGTYPE(K).EQ.1 ).OR.(SEGTYPE(K).EQ.11).OR.&
               (SEGTYPE(K).EQ.21) ) THEN
              DO IED = 1,NVELL(K)-1
                N1 = NBVV(K,IED)
                N2 = NBVV(K,IED+1)
                IF ((N1.EQ.dg%J1).OR.(N1.EQ.dg%J2)) THEN
                  IF ((N2.EQ.dg%J1).OR.(N2.EQ.dg%J2)) THEN
                    NOT_AN_EDGE(N1) = 1
                    NOT_AN_EDGE(N2) = 1
                  ENDIF
                ENDIF
              ENDDO
            ENDIF
          ENDDO
        ENDIF
10130 CONTINUE

!.....Find the interior edges

      DO I = 1,dg%NEDGES
        IEL1 = dg%NEDEL(1,I)
        IEL2 = dg%NEDEL(2,I)
        IF((IEL1.NE.0).AND.(IEL2.NE.0))THEN
          dg%NIEDS = dg%NIEDS + 1
          dg%NIEDN(dg%NIEDS) = I
          dg%NCOUNT(I) = 0
        ENDIF
      ENDDO

!.....Find elevation specified boundary edges

      DO 10131 K = 1,NOPE
        DO IED = 1,NVDLL(K)-1
          N1 = NBDV(K,IED)
          N2 = NBDV(K,IED+1)
          IERROR = 0
          DO JED = 1,dg%NEDGES
            dg%J1 = dg%NEDNO(1,JED)
            dg%J2 = dg%NEDNO(2,JED)
            IF ( (dg%J1.EQ.N1).OR.(dg%J1.EQ.N2) ) THEN
              IF ( (dg%J2.EQ.N1).OR.(dg%J2.EQ.N2) ) THEN
                dg%NEEDS = dg%NEEDS + 1
                dg%NEEDN(dg%NEEDS) = JED
                dg%NCOUNT(JED) = 4
                IERROR = 1
              ENDIF
            ENDIF
          ENDDO
!sb-PDG1
#ifdef CMPI 
          IF(IERROR.EQ.0) THEN
            WRITE(16,*)&
          'ERROR IN PROCESSING OPEN OCEAN BOUNDARY CONDITIONS'
          ENDIF
#else
          IF(IERROR.EQ.0) THEN
            STOP 'ERROR IN PROCESSING OPEN OCEAN BOUNDARY CONDITIONS'
          ENDIF
#endif
!--
        ENDDO
10131 CONTINUE

!.....Find flux specified boundary edges

      JNMM = 0
      ONE_OR_TWO = 1
      DO 10132 K = 1,NBOU
!sb-
        DO 10133 IED = 1,NVELL(K)-1
!--
          N1 = NBVV(K,IED)
          N2 = NBVV(K,IED+1)

          IERROR = 0
          DO JED = 1,dg%NEDGES
            dg%J1 = dg%NEDNO(1,JED)
            dg%J2 = dg%NEDNO(2,JED)
            IF ((NOT_AN_EDGE(N1).EQ.1).AND.(NOT_AN_EDGE(N2).EQ.1)) THEN
!              PRINT*,'NODES ',N1,' AND ',N2,' MAKE UP A BOUNDARY',
!     &                  'SEGMENT THAT IS NOT AN EDGE TO AN ELEMENT; ',
!sb-
!     &                  'THIS EDGE IS SKIPPED'
              GOTO 10133
!--
            ENDIF
            IF ( (N1.EQ.dg%J1).OR.(N1.EQ.dg%J2) ) THEN
              IF ( (N2.EQ.dg%J1).OR.(N2.EQ.dg%J2) ) THEN

                IERROR = 1
                dg%NCOUNT(JED) = 1

!.....Determine the different boundary types

                IF ( (SEGTYPE(K).EQ.0 ).OR.(SEGTYPE(K).EQ.10).OR.&
               (SEGTYPE(K).EQ.20) ) THEN
                  dg%NLEDS = dg%NLEDS + 1
                  dg%NLEDN(dg%NLEDS) = JED
                ENDIF

                IF ( (SEGTYPE(K).EQ.1 ).OR.(SEGTYPE(K).EQ.11).OR.&
               (SEGTYPE(K).EQ.21) ) THEN
                  dg%NLEDS = dg%NLEDS + 1
                  dg%NLEDN(dg%NLEDS) = JED
                ENDIF

                IF ( (SEGTYPE(K).EQ.2 ).OR.(SEGTYPE(K).EQ.12).OR.&
               (SEGTYPE(K).EQ.22) ) THEN
                  dg%NFEDS = dg%NFEDS + 1
                  dg%NFEDN(dg%NFEDS) = JED
                ENDIF
                
                IF ( (SEGTYPE(K).EQ.3 ).OR.(SEGTYPE(K).EQ.13).OR.&
               (SEGTYPE(K).EQ.23) ) THEN
                  dg%NEBSEG = dg%NEBSEG + 1
                  dg%NEBEDS = dg%NEBEDS + 1
                  dg%NEBEDN(dg%NEBEDS) = JED
                  dg%NEBSEGN(dg%NEBSEG) = JED
                ENDIF
                
                IF ( (SEGTYPE(K).EQ.4 ).OR.(SEGTYPE(K).EQ.14).OR.&
               (SEGTYPE(K).EQ.24) ) THEN
                  dg%NIBSEG = dg%NIBSEG + 1
                  dg%NIBEDS = dg%NIBEDS + 1
                  dg%NIBEDN(dg%NIBEDS) = JED
                  dg%NIBSEGN(1,dg%NIBSEG) = JED

!.....Find the edge on the opposite side of the barrier

                  NN1 = dg%BACKNODES(1,dg%NIBSEG)
                  NN2 = dg%BACKNODES(2,dg%NIBSEG)
!                  if (myproc.eq.24) then
!                     write(200+myproc,*) jed,nn1,nn2
!                  endif
                  DO JJED = 1,dg%NEDGES
                    JJ1 = dg%NEDNO(1,JJED)
                    JJ2 = dg%NEDNO(2,JJED)
                    IF ( (NN1.EQ.JJ1).OR.(NN1.EQ.JJ2) ) THEN
                      IF ( (NN2.EQ.JJ1).OR.(NN2.EQ.JJ2) ) THEN
                        dg%NIBEDS = dg%NIBEDS + 1
                        dg%NIBEDN(dg%NIBEDS) = JJED
                        dg%NIBSEGN(2,dg%NIBSEG) = JJED
                        dg%NCOUNT(JJED) = 1
                        ONE_OR_TWO(N1) = 2
                        ONE_OR_TWO(N2) = 2
                      ENDIF
                    ENDIF
                  ENDDO
!                  if (dg%nibsegn(2,dg%nibseg).eq.0) then
!                     write(200+myproc,*) 'error in create_edge_data'
!                     write(200+myproc,*) myproc,dg%nibseg,jed,nn1,nn2
!                  endif
                ENDIF
                
                IF ( (SEGTYPE(K).EQ.30) ) THEN
                  dg%NREDS = dg%NREDS + 1
                  dg%NREDN(dg%NREDS) = JED
                ENDIF
                
              ENDIF
            ENDIF
!sb-
          ENDDO
!--
!sb-PDG1
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
!             STOP 'ERROR IN PROCESSING LAND SEGMENT'
           endif
#endif
!--
10133   CONTINUE
10132 CONTINUE

!.....Check the order of the nodes assigned to an edge - important in
!.....the calculation of the unit normal

      DO I = 1,dg%NEDGES
        N1 = dg%NEDNO(1,I)
        N2 = dg%NEDNO(2,I)
        IEL = dg%NEDEL(1,I)
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
      
!.....Check for missing edges

!sb-PDG1 modified
      NERR = dg%NEDGES - (dg%NIEDS + dg%NLEDS + dg%NEEDS + dg%NFEDS + dg%NIBEDS + dg%NEBEDS +&
                 dg%NREDS)

#ifdef CMPI
      IF(MYPROC.EQ.0) THEN
         WRITE(6,*) '  '
         WRITE(6,*) 'TOTAL NO. OF EDGES = ', dg%NEDGES
         WRITE(6,*) '  '
         WRITE(6,*) 'NO. OF INTERNAL (NON-BOUNDARY) EDGES = ', dg%NIEDS
         WRITE(6,*) 'NO. OF NO-NORMAL FLOW EDGES = ', dg%NLEDS
         WRITE(6,*) 'NO. OF NON-ZERO NORMAL FLOW EDGES = ', dg%NFEDS
         WRITE(6,*) 'NO. OF ELEVATION SPECIFIED EDGES = ', dg%NEEDS
         WRITE(6,*) 'NO. OF EXTERNAL BARRIER EDGES = ', dg%NEBEDS
         WRITE(6,*) 'NO. OF INTERNAL BARRIER EDGES = ', dg%NIBEDS
         WRITE(6,*) 'NO. OF RADIATION EDGES = ', dg%NREDS
         WRITE(6,*) &
        '-----------------------------------------------------'
         WRITE(6,*) 'NO. OF MISSING EDGES = ',NERR
         WRITE(6,*)  ''
      ENDIF
#else
         WRITE(6,*) '  '
         WRITE(6,*) 'TOTAL NO. OF EDGES = ', dg%NEDGES
         WRITE(6,*) '  '
         WRITE(6,*) 'NO. OF INTERNAL (NON-BOUNDARY) EDGES = ', dg%NIEDS
         WRITE(6,*) 'NO. OF NO-NORMAL FLOW EDGES = ', dg%NLEDS
         WRITE(6,*) 'NO. OF NON-ZERO NORMAL FLOW EDGES = ', dg%NFEDS
         WRITE(6,*) 'NO. OF ELEVATION SPECIFIED EDGES = ', dg%NEEDS
         WRITE(6,*) 'NO. OF EXTERNAL BARRIER EDGES = ', dg%NEBEDS
         WRITE(6,*) 'NO. OF INTERNAL BARRIER EDGES = ', dg%NIBEDS
         WRITE(6,*) 'NO. OF RADIATION EDGES = ', dg%NREDS
         WRITE(6,*) &
        '-----------------------------------------------------'
         WRITE(6,*) 'NO. OF MISSING EDGES = ',NERR
         WRITE(6,*)  ''
#endif
      WRITE(16,*) '  '
      WRITE(16,*) 'TOTAL NO. OF EDGES = ', dg%NEDGES
      WRITE(16,*) '  '
      WRITE(16,*) 'NO. OF INTERNAL (NON-BOUNDARY) EDGES = ', dg%NIEDS
      WRITE(16,*) 'NO. OF NO-NORMAL FLOW EDGES = ', dg%NLEDS
      WRITE(16,*) 'NO. OF NON-ZERO NORMAL FLOW EDGES = ', dg%NFEDS
      WRITE(16,*) 'NO. OF ELEVATION SPECIFIED EDGES = ', dg%NEEDS
      WRITE(16,*) 'NO. OF EXTERNAL BARRIER EDGES = ', dg%NEBEDS
      WRITE(16,*) 'NO. OF INTERNAL BARRIER EDGES = ', dg%NIBEDS
      WRITE(16,*) 'NO. OF RADIATION EDGES = ', dg%NREDS
      WRITE(16,*) &
        '-----------------------------------------------------'
      WRITE(16,*) 'NO. OF MISSING EDGES = ',NERR
      WRITE(16,*)  ''

      DO I = 1,dg%NEDGES
        IF (dg%NCOUNT(I).LT.0) THEN
          N1 = dg%NEDNO(1,I)
          N2 = dg%NEDNO(2,I)
!sb-PDG1
#ifdef CMPI
!          WRITE(6,*)' '
!          WRITE(6,*)'EDGE ',I,' IS MADE UP OF NODES ',N1,' AND ',N2
!          WRITE(6,*)'EDGE ',I, ' IS NEITHER AN INTERNAL NOR A BOUNDARY',
!     &       'MAKE SURE IF THIS IS DUE TO THE DOMAIN DECOMPOSITION'
#else
          WRITE(16,*)' '
          WRITE(16,*)'EDGE',I,',MADE UP OF NODES',N1,'AND',N2,', IS NOT'
          WRITE(16,*)'AN INTERNAL(NON-BOUNDARY) EDGE OR A BOUNDARY EDGE'
          WRITE(16,*)'ASSUMING EDGE',I,'IS A NO-NORMAL FLOW EDGE !!!'
!          STOP
#endif
!--
          dg%NLEDS = dg%NLEDS + 1
          dg%NLEDN(dg%NLEDS) = I
        ENDIF
      ENDDO
      
!.....Add internal barrier edges to land edge table for wet-dry
!.....post-processing

      DO I = 1,dg%NIBEDS
        dg%NLEDN(dg%NLEDS+I) = dg%NIBEDN(I)
      ENDDO
      
!.....Print out the interior edges

      WRITE(17,*) dg%NIEDS,'       ! NUMBER OF INTERNAL EDGES'
      DO I = 1,dg%NIEDS
        WRITE(17,*) I,dg%NIEDN(I),dg%nedno(1,dg%niedn(i)),dg%nedno(2,dg%niedn(i))
      ENDDO

!.....Prin out land edges info.

      WRITE(17,*) dg%NLEDS,'       ! NUMBER OF NO-NORMAL FLOW EDGES'
      IF (dg%NLEDS.GT.0) THEN
        DO I = 1,dg%NLEDS
          WRITE(17,*) I, dg%NLEDN(I),dg%nedno(1,dg%nledn(i)),          dg%nedno(2,dg%nledn(i))
        ENDDO
      ENDIF
      
!.....Print out elevation specified edge info.

      WRITE(17,*) dg%NEEDS,'       ! NUMBER OF ELEVATION SPECIFIED EDGES'
      IF (dg%NEEDS.GT.0) THEN
        DO I = 1,dg%NEEDS
          WRITE(17,*) I, dg%NEEDN(I)
        ENDDO
      ENDIF

!.....Print out non-zero flow edges info.

      WRITE(17,*) dg%NFEDS,'       ! NUMBER OF FLOW SPECIFIED EDGES'
      IF (dg%NFEDS.GT.0) THEN
        DO I = 1,dg%NFEDS
          WRITE(17,*) I, dg%NFEDN(I)
        ENDDO
      ENDIF
      
!.....Print out external barrier edge info.

      WRITE(17,*) dg%NEBEDS,'       ! NUMBER OF EXTERNAL BARRIER EDGES'
      IF (dg%NEEDS.GT.0) THEN
        DO I = 1,dg%NEBEDS
          WRITE(17,*) I, dg%NEBEDN(I),dg%nedno(1,dg%nebedn(i)),          dg%nedno(2,dg%nebedn(i))
        ENDDO
      ENDIF
      
!.....Print out internal barrier edge info.

      WRITE(17,*) dg%NIBEDS,'       ! NUMBER OF INTERNAL BARRIER EDGES'
      IF (dg%NIBEDS.GT.0) THEN
        DO I = 1,dg%NIBEDS
          WRITE(17,*) I, dg%NIBEDN(I),dg%nedno(1,dg%nibedn(i)),          dg%nedno(2,dg%nibedn(i))
!          PRINT*,'INTERNAL BARRIER EDGE =',I
!          PRINT*,'IS MADE UP OF NODES',dg%NEDNO(1,dg%NIBEDN(I)),'AND',
!     &                                 dg%NEDNO(2,dg%NIBEDN(I))
        ENDDO
      ENDIF

!.....Print out radiation edges info.

      WRITE(17,*) dg%NREDS, '      ! NUMBER OF RADIATION EDGES'
      IF (dg%NREDS.GT.0) THEN
        DO I = 1,dg%NREDS
          WRITE(17,*) I, dg%NREDN(I)
        ENDDO
      ENDIF
      
!.....Construct global edge to local edge (1,2, or 3) table

      dg%NEDSD(:,:) = 0.D0
      DO I = 1,dg%MNED
        DO K = 1,2
          IF (dg%NEDEL(K,I).NE.0) THEN
            N1 = NELED(1,dg%NEDEL(K,I))
            N2 = NELED(2,dg%NEDEL(K,I))
            N3 = NELED(3,dg%NEDEL(K,I))
            IF (N1.EQ.I) dg%NEDSD(K,I) = 1
            IF (N2.EQ.I) dg%NEDSD(K,I) = 2
            IF (N3.EQ.I) dg%NEDSD(K,I) = 3
          ENDIF
        ENDDO
      ENDDO

      RETURN
      END SUBROUTINE CREATE_EDGE_DATA
