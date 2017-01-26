      SUBROUTINE metis2(nparts)

      USE globals, ONLY: nn,ne,ned,ged2el,part      
      USE messenger2, ONLY: myrank

      IMPLICIT NONE
      
      INTEGER :: nparts
      
      INTEGER :: ed,nd,j,el
      INTEGER :: el1,el2,nd2
      INTEGER :: tot
      INTEGER :: tmp
      
      INTEGER, ALLOCATABLE, DIMENSION(:) :: nadjels
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: adjels
      
      INTEGER :: numflag,wgtflag
      INTEGER :: edgecut
      INTEGER, ALLOCATABLE, DIMENSION(:) :: xadj
      INTEGER, ALLOCATABLE, DIMENSION(:) :: adjncy
      INTEGER, ALLOCATABLE, DIMENSION(:) :: vwgt
      INTEGER, ALLOCATABLE, DIMENSION(:) :: adjwgt
      INTEGER :: options(5)
      
      edgecut = 0
      
      ! ---------------------------------------
      ! Elements are treated as the graph nodes
      ! ---------------------------------------      
      
      ! Find the nodes adjacent each node
      ALLOCATE(nadjels(ne))
      ALLOCATE(adjels(4,ne))
      
      nadjels = 0
      DO ed = 1,ned
        el1 = ged2el(1,ed) ! find the node numbers on each edge
        el2 = ged2el(2,ed)
        
        IF(el2 /= 0) THEN
          nadjels(el1) = nadjels(el1) + 1 ! count the nodes adjacent to node n1
          nadjels(el2) = nadjels(el2) + 1 ! count the nodes adjacent to node n2
        
          adjels(nadjels(el1),el1) = el2 ! node n2 is adjacent to node n1
          adjels(nadjels(el2),el2) = el1 ! node n1 is adjacent to node n2       
        ENDIF
      ENDDO
      
      ! sort the node numbers adjacent to each node 
      ! (I'm not sure why this matters, but it seems to reduce the edge cut count)
      DO el = 1,ne
        DO ed = 1,nadjels(el)
          DO j = ed,nadjels(el)
            IF(adjels(j,el) < adjels(ed,el)) THEN
              tmp = adjels(j,el)
              adjels(j,el) = adjels(ed,el)
              adjels(ed,el) = tmp
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      
!       DO el = 1,ne
!         PRINT("(I7,5x,I7,5x,10(I5))"), el, nadjels(el), (adjels(ed,el), ed = 1,nadjels(el))
!       ENDDO      
      
      ! Compute the node weights
      ! (ADCPREP uses the number of adjacent nodes)
      ALLOCATE(vwgt(ne))
      
      DO el = 1,ne
        vwgt(el) = nadjels(el)
!         PRINT*, nd,vwgt(nd)
      ENDDO
      
      ! Create the adjacency arrays used by metis
      ! Also compute the edge weight 
      ! (ADCPREP uses the sum of the adjacent nodes for each node on the edge)
      ALLOCATE(xadj(ne+1))
      ALLOCATE(adjncy(2*ned),adjwgt(2*ned))
      
      xadj(1) = 1
      tot = 0
      DO el = 1,ne
        DO ed = 1,nadjels(el)
          tot = tot + 1
          el2 = adjels(ed,el)
          adjncy(tot) = el2 
          adjwgt(tot) = nadjels(el) + nadjels(el2)
!           PRINT*, tot, adjncy(tot), adjwgt(tot)
        ENDDO
        xadj(el+1) = tot+1
      ENDDO      
      
      numflag = 1 ! use fortran numbering
      wgtflag = 3 ! use edge and node weights
      
      options(1) = 1 ! use non-default options
      options(2) = 3 ! use sorted heavy edge matching type
      options(3) = 1 ! use multilevel recursive bisection algorithm during initial partitioning 
      options(4) = 3 ! minimize connectivity among the subdomains 
      options(5) = 0 ! always set to zero

      IF (myrank == 0) THEN
        PRINT "(A)", " "
        PRINT "(A)", "---------------------------------------------"
        PRINT "(A)", "      Metis Partitioning Information         "
        PRINT "(A)", "---------------------------------------------"  
        PRINT "(A)", " "
        PRINT*, "Partitioning into ", nparts, " subdomains"      
      ENDIF
      
      ALLOCATE(part(ne))
      
      IF (nparts > 1) THEN
      
        CALL METIS_PartGraphKway(ne,xadj,adjncy,vwgt,adjwgt,wgtflag, &
                                 numflag,nparts,options,edgecut,part)
                               
      ELSE IF (nparts == 1) THEN
      
        DO el = 1,ne
          part(el) = 1
        ENDDO
      
      ENDIF
      
      IF (myrank == 0) THEN
        PRINT*,"edgecut = ", edgecut
      ENDIF
      
      RETURN 
      END SUBROUTINE metis2