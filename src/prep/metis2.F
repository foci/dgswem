        SUBROUTINE METIS2(MNPROC2,MNP2,MNE2,MNEI2,NWEIR2,
     &                    NNEG2,WEIR2,WEIRD2)
        USE PRE_GLOBAL 
C
C----------------------------------------------------------------------
C  INTERFACE ROUTINE FOR PADCIRC TO USE THE METIS 4.0 LIBRARY 
C  A GRAPH PARTITION TOOL, FOR DOMAIN DECOMPOSITION
C  Version 1.1  vjp  12/15/99
C  vjp 7/28/2001 changed 2*MNEI2 to 3*MNEI2 to prevent out of bounds
C  Redimensioned CO_NODES and IDUAL and added bounds test
C----------------------------------------------------------------------
C
        INTEGER MNED
        INTEGER I, J, IEL, INODE, JNODE, NCOUNT, ITOT
        INTEGER MNPROC2,MNP2,MNE2,MNEI2,NWEIR2
        INTEGER NNEG2(3,MNE2), WEIR2(*), WEIRD2(*)
C
#ifdef CRAY
        INTEGER*8 WEIGHTFLAG, NUMFLAG, NPARTS, OPTIONS(5)
        INTEGER*8 EDGECUT,NNODG8,OPTYPE,NBYTES
#else
        INTEGER  WEIGHTFLAG, NUMFLAG, NPARTS, OPTIONS(5)
        INTEGER  EDGECUT,NNODG8,OPTYPE,NBYTES
#endif
C
        INTEGER,ALLOCATABLE :: IDUAL(:),ITVECT(:)
        INTEGER,ALLOCATABLE :: XADJ(:), ADJNCY(:)
        INTEGER,ALLOCATABLE :: VWGTS(:), EWGTS(:)
        INTEGER,ALLOCATABLE :: CO_NODES(:,:),NEDGES(:)
        INTEGER,ALLOCATABLE :: PROC2(:)
        EXTERNAL metis_partgraphkway,metis_estimatememory
C
        MNED = 3*MNP2
        ALLOCATE ( IDUAL(MNP2),ITVECT(3*MNEI2) )
        ALLOCATE ( XADJ(MNP2+1), ADJNCY(2*MNED) )
        ALLOCATE ( VWGTS(MNP2), EWGTS(2*MNED) )
        ALLOCATE ( CO_NODES(3*MNEI2,MNP2),NEDGES(MNP2) )
        ALLOCATE ( PROC2(MNP2) )
C      
        DO I=1, MNP2
           PROC2(I) = 0
        ENDDO
c       DO I=1, MNE2
c          print *, I, NNEG2(1,I),NNEG2(2,I),NNEG2(3,I)
c       ENDDO

C
C--COMPUTE INDEX OF WEIR DUAL WHICH IS ZERO IF NOT A WEIR NODE
C
        DO INODE=1, MNP2   
           IDUAL(INODE) = 0
        ENDDO
        DO J=1, NWEIR2
           IDUAL(WEIR2(J)) = WEIRD2(J)
           IDUAL(WEIRD2(J)) = WEIR2(J)
        ENDDO
C
C--COMPUTE CO_NODES LISTS AND NUMBER OF EDGES CONTAINING A NODE
C
        DO INODE = 1,MNP2
           NEDGES(INODE) = 0
        ENDDO
C
        DO IEL=1, MNE2
           INODE = NNEG2(1,IEL)
           NCOUNT = NEDGES(INODE)
           CO_NODES(NCOUNT+1,INODE) = NNEG2(2,IEL)
           CO_NODES(NCOUNT+2,INODE) = NNEG2(3,IEL)
           IF (IDUAL(INODE).EQ.0) THEN
             NEDGES(INODE) = NCOUNT+2
           ELSE
             CO_NODES(NCOUNT+3,INODE) = IDUAL(INODE)
             NEDGES(INODE) = NCOUNT+3
           ENDIF
           IF (NEDGES(INODE).GT.3*MNEI2) THEN
             print *, "inode = ", INODE
             print *, "nedges(inode) = ", NEDGES(INODE)
           ENDIF
        ENDDO
C
        DO IEL=1, MNE2
           INODE = NNEG2(2,IEL)
           NCOUNT = NEDGES(INODE)
           CO_NODES(NCOUNT+1,INODE) = NNEG2(1,IEL)
           CO_NODES(NCOUNT+2,INODE) = NNEG2(3,IEL)
           IF (IDUAL(INODE).EQ.0) THEN
             NEDGES(INODE) = NCOUNT+2
           ELSE
             CO_NODES(NCOUNT+3,INODE) = IDUAL(INODE)
             NEDGES(INODE) = NCOUNT+3
           ENDIF
           IF (NEDGES(INODE).GT.3*MNEI2) THEN
             print *, "inode = ", INODE
             print *, "nedges(inode) = ", NEDGES(INODE)
           ENDIF
        ENDDO
C
        DO IEL=1, MNE2
           INODE = NNEG2(3,IEL)
           NCOUNT = NEDGES(INODE)
           CO_NODES(NCOUNT+1,INODE) = NNEG2(1,IEL)
           CO_NODES(NCOUNT+2,INODE) = NNEG2(2,IEL)
           IF (IDUAL(INODE).EQ.0) THEN
             NEDGES(INODE) = NCOUNT+2
           ELSE
             CO_NODES(NCOUNT+3,INODE) = IDUAL(INODE)
             NEDGES(INODE) = NCOUNT+3
           ENDIF
           IF (NEDGES(INODE).GT.3*MNEI2) THEN
             print *, "inode = ", INODE
             print *, "nedges(inode) = ", NEDGES(INODE)
           ENDIF
        ENDDO
C
C  REMOVE REDUNDANCY IN NODE LISTS
C
        DO INODE = 1,MNP2   
           IF (NEDGES(INODE).GT.3*MNEI2) THEN
             print *, "number edges exceeds 2*mnei in metis graph"
             print *, "this occurs at node = ",INODE
             print *, "nedges(inode) = ", NEDGES(INODE)
             stop 'metis interface routine'
           ENDIF
           DO J=1, NEDGES(INODE)
              ITVECT(J) = CO_NODES(J,INODE)
           ENDDO
           IF (NEDGES(INODE).GT.1) THEN
             NCOUNT = NEDGES(INODE)
             CALL SORT(NCOUNT,ITVECT)
             NCOUNT = 1
             JNODE = ITVECT(1)
             CO_NODES(1,INODE) = JNODE
             DO J=2, NEDGES(INODE)
                IF (ITVECT(J).NE.JNODE) THEN
                  NCOUNT = NCOUNT + 1
                  JNODE = ITVECT(J)
                  CO_NODES(NCOUNT,INODE) = JNODE
                ENDIF
             ENDDO
           ENDIF
           NEDGES(INODE) = NCOUNT
        ENDDO
C
C  COMPUTE WEIGHTS OF THE GRAPH'S VERTICES
C
        DO INODE = 1,MNP2   
         VWGTS(INODE) = NEDGES(INODE)
        ENDDO
C
C--COMPUTE ADJACENCY LIST OF GRAPH AND ITS EDGE WEIGHTS
C
        XADJ(1) = 1
        ITOT = 0
        DO INODE = 1,MNP2
        DO J = 1, NEDGES(INODE)
           ITOT = ITOT + 1
           JNODE = CO_NODES(J,INODE)
           ADJNCY(ITOT) = JNODE
           EWGTS(ITOT)  = NEDGES(JNODE)
        ENDDO
        XADJ(INODE+1) = ITOT+1
        ENDDO
C
C--CALL K-WAY METIS FOR PARTITIONING
C
        NUMFLAG  = 1
        NPARTS = MNPROC2
        OPTIONS(1) = 1
        OPTIONS(2) = 3
        OPTIONS(3) = 1
        OPTIONS(4) = 3   !  minimize number of co-domains
        OPTIONS(5) = 0
        WEIGHTFLAG = 3   ! use weights for nodes and edges
C
        NNODG8 = MNP2  
        OPTYPE = 2
        CALL metis_estimatememory(NNODG8,XADJ,ADJNCY,NUMFLAG,
     &        OPTYPE,NBYTES)
C
        print *, ""
        print *, "Grid Partition Data"
        print *, "METIS 4.0 will require approximately ",nbytes," bytes"
C
        CALL metis_partgraphkway( NNODG8,XADJ,ADJNCY,VWGTS,EWGTS,
     &      WEIGHTFLAG,NUMFLAG,NPARTS,OPTIONS,EDGECUT,PROC2)
C
        DO I=1, MNP2
         WRITE(151,*)  PROC2(I)
        ENDDO
        CLOSE(151)
C
        print *, "Total Edges Cut = ",EDGECUT 
C
        RETURN
        END
