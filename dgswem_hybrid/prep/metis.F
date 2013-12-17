        SUBROUTINE METIS()
        USE PRE_GLOBAL 
C
C----------------------------------------------------------------------
C  INTERFACE ROUTINE FOR PADCIRC TO USE THE METIS 4.0 LIBRARY 
C  A GRAPH PARTITION TOOL, FOR DOMAIN DECOMPOSITION
C  Version 1.1  vjp  12/15/99
C  vjp 7/28/2001 changed 2*MNEI to 3*MNEI to prevent out of bounds
C  Redimensioned CO_NODES and IDUAL and added bounds test
C----------------------------------------------------------------------
C
        INTEGER MNED, IDUMP
        INTEGER I, J, IEL, INODE, JNODE, NCOUNT, ITOT, NEDGETOT
C
#ifdef CRAY
Cvjp  redefined as integer*8 to be compabible on Cray machines
        INTEGER*8 WEIGHTFLAG, NUMFLAG, NPARTS, OPTIONS(5)
        INTEGER*8 EDGECUT,NNODG8,OPTYPE,NBYTES
#else
        INTEGER  WEIGHTFLAG, NUMFLAG, NPARTS, OPTIONS(5)
        INTEGER  EDGECUT,NNODG8,OPTYPE,NBYTES
#endif
C
        INTEGER,ALLOCATABLE :: IDUAL(:),ITVECT(:),ITVECT2(:)
        INTEGER,ALLOCATABLE :: XADJ(:), ADJNCY(:)
        INTEGER,ALLOCATABLE :: VWGTS(:), EWGTS(:)
        INTEGER,ALLOCATABLE :: CO_NODES(:,:),NEDGES(:)
        EXTERNAL metis_partgraphkway,metis_estimatememory
C
        MNED = 3*MNP
        ALLOCATE ( ITVECT(3*MNEI),ITVECT2(3*MNEI) )
        ALLOCATE ( IDUAL(MNP) )
        ALLOCATE ( XADJ(MNP+1), ADJNCY(2*MNED) )
        ALLOCATE ( VWGTS(MNP), EWGTS(2*MNED) )
        ALLOCATE ( CO_NODES(3*MNEI,MNP),NEDGES(MNP) )

C
C--COMPUTE INDEX OF WEIR DUAL WHICH IS ZERO IF NOT A WEIR NODE
C
        DO INODE=1, MNP   
           IDUAL(INODE) = 0
        ENDDO
        print *, "This model has ",NWEIR," weir node pairs"
        DO J=1, NWEIR
           IDUAL(WEIR(J)) = WEIRD(J)
           IDUAL(WEIRD(J)) = WEIR(J)
        ENDDO
C
C--COMPUTE CO_NODES LISTS AND NUMBER OF EDGES CONTAINING A NODE
C
        DO INODE = 1,MNP
           NEDGES(INODE) = 0
        ENDDO
C
        DO IEL=1, MNE  
           INODE = NNEG(1,IEL)
           NCOUNT = NEDGES(INODE)
           CO_NODES(NCOUNT+1,INODE) = NNEG(2,IEL)
           CO_NODES(NCOUNT+2,INODE) = NNEG(3,IEL)
           IF (IDUAL(INODE).EQ.0) THEN
             NEDGES(INODE) = NCOUNT+2
           ELSE
             CO_NODES(NCOUNT+3,INODE) = IDUAL(INODE)
             NEDGES(INODE) = NCOUNT+3
           ENDIF
           IF (NEDGES(INODE).GT.3*MNEI) THEN  
             print *, "inode = ", INODE
             print *, "nedges(inode) = ", NEDGES(INODE)
           ENDIF
        ENDDO
C
        DO IEL=1, MNE  
           INODE = NNEG(2,IEL)
           NCOUNT = NEDGES(INODE)
           CO_NODES(NCOUNT+1,INODE) = NNEG(1,IEL)
           CO_NODES(NCOUNT+2,INODE) = NNEG(3,IEL)
           IF (IDUAL(INODE).EQ.0) THEN
             NEDGES(INODE) = NCOUNT+2
           ELSE
             CO_NODES(NCOUNT+3,INODE) = IDUAL(INODE)
             NEDGES(INODE) = NCOUNT+3
           ENDIF
           IF (NEDGES(INODE).GT.3*MNEI) THEN  
             print *, "inode = ", INODE
             print *, "nedges(inode) = ", NEDGES(INODE)
           ENDIF
        ENDDO
C
        DO IEL=1, MNE    
           INODE = NNEG(3,IEL)
           NCOUNT = NEDGES(INODE)
           CO_NODES(NCOUNT+1,INODE) = NNEG(1,IEL)
           CO_NODES(NCOUNT+2,INODE) = NNEG(2,IEL)
           IF (IDUAL(INODE).EQ.0) THEN
             NEDGES(INODE) = NCOUNT+2
           ELSE
             CO_NODES(NCOUNT+3,INODE) = IDUAL(INODE)
             NEDGES(INODE) = NCOUNT+3
           ENDIF
           IF (NEDGES(INODE).GT.3*MNEI) THEN  
             print *, "inode = ", INODE
             print *, "nedges(inode) = ", NEDGES(INODE)
           ENDIF
        ENDDO
C
C  REMOVE REDUNDANCY IN NODE LISTS
C
        NEDGETOT = 0           !  This will be twice number of edges
        DO INODE = 1,MNP   
           IF (NEDGES(INODE).GT.3*MNEI) THEN  
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
             JNODE = ITVECT(1)
             CO_NODES(1,INODE) = JNODE
             NCOUNT = 1
             DO J=2, NEDGES(INODE)
                IF (ITVECT(J).NE.JNODE) THEN
                  NCOUNT = NCOUNT + 1
                  JNODE = ITVECT(J)
                  CO_NODES(NCOUNT,INODE) = JNODE
                ENDIF
             ENDDO
           ELSE
             print *, "node = ",INODE," is isolated"
             stop 'vic'
           ENDIF
           NEDGES(INODE) = NCOUNT
           NEDGETOT = NEDGETOT + NCOUNT
        ENDDO
        NEDGETOT = NEDGETOT/2
        print *, "edge count = ",nedgetot
C
C  COMPUTE WEIGHTS OF THE GRAPH'S VERTICES
C
        DO INODE = 1,MNP   
         VWGTS(INODE) = NEDGES(INODE)
        ENDDO
C
C--COMPUTE ADJACENCY LIST OF GRAPH AND ITS EDGE WEIGHTS
C
        XADJ(1) = 1
        ITOT = 0
        DO INODE = 1,MNP
        DO J = 1, NEDGES(INODE)
           ITOT = ITOT + 1
           JNODE = CO_NODES(J,INODE)
           ADJNCY(ITOT) = JNODE
           EWGTS(ITOT)  = (NEDGES(JNODE)+NEDGES(INODE))
        ENDDO
        XADJ(INODE+1) = ITOT+1
        ENDDO
C
C--CALL K-WAY METIS FOR PARTITIONING
C
        NUMFLAG  = 1
        NPARTS = MNPROC
        OPTIONS(1) = 1
cvjp 3/29/2001 workaround to limit memory usage
c       OPTIONS(2) = 1
        OPTIONS(2) = 3
c
        OPTIONS(3) = 1
c
cvjp 3/29/2001 workaround to limit memory usage
c       OPTIONS(4) = 1   !  
c       OPTIONS(4) = 2   !  
        OPTIONS(4) = 3   !  minimize number of co-domains
        OPTIONS(5) = 0
c
        WEIGHTFLAG = 3   ! use weights for nodes and edges
C
        NNODG8 = MNP  
        OPTYPE = 2
        CALL metis_estimatememory(NNODG8,XADJ,ADJNCY,NUMFLAG,
     &        OPTYPE,NBYTES)
 
        print *, ""
        print *, "Grid Partition Data"
        print *, "METIS 4.0 will require approximately ",nbytes," bytes"
C
        CALL metis_partgraphkway( NNODG8,XADJ,ADJNCY,VWGTS,EWGTS,
     &      WEIGHTFLAG,NUMFLAG,NPARTS,OPTIONS,EDGECUT,PROC)
C
        print *, "Total Edges Cut = ",EDGECUT 
C
C Dump graph to a file for debugging
C
        IDUMP = 1
        IF (IDUMP.EQ.1) THEN
           OPEN(FILE='metis_graph.txt',UNIT=99)
           WRITE(99,100) MNP, NEDGETOT, 11, 1
           DO INODE=1, MNP
              WRITE(99,200) VWGTS(INODE),
     &   (CO_NODES(J,INODE), EWGTS(XADJ(INODE)+J-1),J=1,NEDGES(INODE))
           ENDDO
           CLOSE(99)
c         stop 'vic'
        ENDIF
C
 100    FORMAT(4I10)
 200    FORMAT(100I10)
C
        RETURN
        END
