      MODULE read_grid
      
      USE globals, ONLY: solution,coarse,fine,base
      
      CONTAINS
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
      SUBROUTINE read_grids()
      
      USE allocation, ONLY: sizes
      
      IMPLICIT NONE
      
      INTEGER :: i,j
      
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", "           Coarse Grid Information           "
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", " "      
      
      CALL sizes(coarse)
      CALL read_nodes(coarse)      
      
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", "            Fine Grid Information            "
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", " "     
     
      CALL sizes(fine)
      CALL read_nodes(fine)
      
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", "            Base Grid Information            "
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", " "     
      
      CALL sizes(base)
      CALL read_nodes(base)      
      
      RETURN
      END SUBROUTINE read_grids
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE read_nodes(sol)
      
      USE globals, ONLY: pres,solution,exclude_bndel
      USE allocation, ONLY: grid_alloc

      IMPLICIT NONE
      INTEGER :: i,j,k,el,n,nd
      INTEGER :: ne,nn,ctp,nvert,nbseg,btype,nope,nbou
      INTEGER :: mnepn,n1,n2,found            
      INTEGER :: curved_grid
      INTEGER, ALLOCATABLE, DIMENSION(:) :: vflag     
      REAL(pres), ALLOCATABLE, DIMENSION(:,:) :: vxy_temp
      TYPE(solution) :: sol
  
      ! open fort.14 grid file
      OPEN(UNIT = 14, FILE = sol%grid_file)      
      
      PRINT "(A,A)", "Grid file: ", sol%grid_file                          
                       
      ! read in name of grid
      READ(14,"(A)"), sol%grid_name                                         

      PRINT "(A,A)", "Grid name: ", sol%grid_name

      ! read in number of elements and number of nodes
      READ(14,*), sol%ne, sol%nn                                                

      PRINT "(A,I5)", "Number of elements: ", sol%ne
      PRINT "(A,I5)", "Number of nodes: ", sol%nn
      PRINT*, " "
      
      nn = sol%nn
      ne = sol%ne
      ctp = sol%ctp

      CALL grid_alloc(1,sol)

      curved_grid = 0
      
      ! read in node coordinates and depths
      DO i = 1,nn                                                      
        READ(14,*), j, sol%xy(1,j), sol%xy(2,j), sol%depth(j)
      ENDDO
!       PRINT "(A)", "Node coordinates and depth: "
!       DO i = nn-15,nn !nn
!         PRINT "(I5,3(F11.3,3x))", i,sol%xy(1,i), sol%xy(2,i), sol%depth(i)
!       ENDDO
!       PRINT*, " "

      ! read in element connectivity
      sol%ect = 0
      DO i = 1,ne
        READ(14,*) el,sol%nelnds(el),(sol%ect(j,el),j = 1,sol%nelnds(el))
        IF (sol%nelnds(el) == 3) THEN
          sol%el_type(el) = 1
        ELSE IF (sol%nelnds(el) == 4) THEN
          sol%el_type(el) = 2
        ELSE IF (sol%nelnds(el) == (sol%ctp+1)*(sol%ctp+2)/2) THEN
          sol%el_type(el) = 3
          curved_grid = 1
        ELSE IF (sol%nelnds(el) == (sol%ctp+1)*(sol%ctp+1)) THEN
          sol%el_type(el) = 4
          curved_grid = 1
        ELSE
          PRINT*, "Element type not supported or ctp not compatible with grid"
          STOP
        ENDIF 
        
        DO j = 1,sol%nelnds(el)
          sol%elxy(j,el,1) = sol%xy(1,sol%ect(j,el))
          sol%elxy(j,el,2) = sol%xy(2,sol%ect(j,el))
          sol%elhb(j,el)   = sol%depth(sol%ect(j,el))
        ENDDO      
        
      ENDDO
      
      PRINT "(A,I5)", "Curved grid = ",curved_grid
      PRINT*, " "
      
      
      ! find element verticies
      ALLOCATE(vxy_temp(2,sol%nn),vflag(sol%nn))
      vflag = 0
      
      IF(curved_grid == 1) THEN
        DO i = 1,ne
          nvert = sol%nverts(sol%el_type(i))
          DO j = 1,nvert
            sol%vct(j,i) = sol%ect(ctp*(j-1)+1,i)
          ENDDO
        ENDDO        
      ELSE ! in case grid isn't curved but ctp > 1
        DO i = 1,ne
          nvert = sol%nverts(sol%el_type(i))
          DO j = 1,nvert
            sol%vct(j,i) = sol%ect(j,i) 
          ENDDO
        ENDDO           
      ENDIF

      
      n = 0
      DO i = 1,ne
        nvert = sol%nverts(sol%el_type(i))
        DO j = 1,nvert
          nd = sol%vct(j,i)
          IF (vflag(nd) == 0) THEN  
            n = n + 1            
            sol%vxyn(n) = nd
            vxy_temp(1,n) = sol%xy(1,nd)
            vxy_temp(2,n) = sol%xy(2,nd)
            vflag(nd) = 1
          ENDIF
        ENDDO
      ENDDO         

      DO i = 1,n
        sol%vxy(1,i) = vxy_temp(1,i)
        sol%vxy(2,i) = vxy_temp(2,i)
      ENDDO
      
      sol%mnelnds = maxval(sol%nelnds)
      
!       PRINT "(A)", "Element connectivity table: "
!       DO i = 1,15
!         PRINT "(2(I5,3x),8x,9(I5,3x))", i,sol%nelnds(i),(sol%ect(j,i),j=1,sol%nelnds(i))
!       ENDDO
!       PRINT*, " "


      ! Find elements associated with each node            
      ALLOCATE(sol%nepn(sol%nn))
      
      sol%nepn(:) = 0
      DO el = 1,sol%ne
        nvert = sol%nverts(sol%el_type(el))
        DO nd = 1,nvert
          n1 = sol%vct(nd,el)
          sol%nepn(n1) = sol%nepn(n1) + 1
        ENDDO
      ENDDO
      
      mnepn = maxval(sol%nepn)
      
      ALLOCATE(sol%epn(mnepn,sol%nn))
      
      sol%nepn(:) = 0
      DO el = 1,sol%ne
        nvert = sol%nverts(sol%el_type(el))
        DO nd = 1,nvert
          n1 = sol%vct(nd,el)
          sol%nepn(n1) = sol%nepn(n1) + 1
          sol%epn(sol%nepn(n1),n1) = el
        ENDDO
      ENDDO 


      ! Read in open boundaries
      READ(14,*) sol%nope  ! number of open boundaries                                                 
      READ(14,*) sol%neta  ! number of total elevation specified boundary nodes
      
      nope = sol%nope
      
      CALL grid_alloc(2,sol)

      DO i = 1,nope                                                     
        READ(14,*), nbseg  ! read in # of nodes in segment, boundary type
        sol%obseg(i) = nbseg
        DO j = 1,nbseg
          READ(14,*) sol%obnds(j,i) ! read in open boundary node numbers
        ENDDO
      ENDDO
!       PRINT "(A)", "Open boundary segments:"
!       DO i = 1,nope
!         nbseg = obseg(i)
!         PRINT "(A,I5,A,I5,A)", "Open boundary segment ",i," contains ",nbseg," nodes"
!         DO j = 1,nbseg
!           PRINT "(I5)",obnds(j,i)
!         ENDDO
!       ENDDO
!       PRINT*, " "
!       

      ! read in normal flow boundaries
      READ(14,*) sol%nbou  ! number of normal flow boundaries
      READ(14,*) sol%nvel  ! total number of normal flow nodes

      nbou = sol%nbou
       
      CALL grid_alloc(3,sol)

      DO i = 1,nbou
        READ(14,*), nbseg, btype ! read in # of nodes in segment, boundary type
        sol%fbseg(1,i) = nbseg
        sol%fbseg(2,i) = btype
        DO j = 1,nbseg
          READ(14,*), sol%fbnds(j,i)  ! read in normal flow boundary node numbers
        ENDDO
        IF (btype == 1 .OR. btype == 11 .OR. btype == 21) THEN
          IF (sol%fbnds(nbseg,i) /= sol%fbnds(1,i)) THEN
            sol%fbnds(nbseg+1,i) = sol%fbnds(1,i)  ! close island boundaries
            sol%fbseg(1,i) = sol%fbseg(1,i) + 1
          ENDIF
        ENDIF
        
      ENDDO
!       PRINT "(A)", "Normal flow boundary segments: "
!       DO i = 1,nbou
!         nbseg = fbseg(1,i)
!         btype = fbseg(2,i)
!         PRINT "(A,I3,A,I3,A,I5,A)", "Normal flow boundary segment ",i," type ",btype, " contains ",nbseg," nodes"
!         DO j = 1,nbseg
!           PRINT "(I5)", fbnds(j,i)
!         ENDDO
!       ENDDO
!       PRINT*, " "


      CLOSE(14) 
      
      ! Find elements on land boundaries (these are left out of the L2 error calculation in case they're curve)
      ALLOCATE(sol%bndel(sol%ne))
      sol%bndel = 0
      
      IF (exclude_bndel) THEN
      DO i = 1,sol%nbou
        btype = sol%fbseg(2,i)
        nbseg = sol%fbseg(1,i)-1
        IF (btype == 0 .OR. btype == 10 .OR. btype == 20 .OR. &
            btype == 2 .OR. btype == 12 .OR. btype == 22) THEN
          DO j = 1,nbseg
            n1 = sol%fbnds(j,i)
            n2 = sol%fbnds(j+1,i)
      elem: DO el = 1,sol%ne
              found = 0
              DO nd = 1,sol%nelnds(el)
                IF (sol%ect(nd,el) == n1 .OR. sol%ect(nd,el) == n2) THEN
                  found = found + 1
                ENDIF
              ENDDO
              IF (found == 2) THEN
                sol%bndel(el) = 1
                EXIT elem
              ENDIF
            ENDDO elem
          ENDDO
        ENDIF
      ENDDO

      DO i = 1,sol%nope
        nbseg = sol%obseg(i)-1
          DO j = 1,nbseg
            n1 = sol%obnds(j,i)
            n2 = sol%obnds(j+1,i)
     elem2: DO el = 1,sol%ne
              found = 0
              DO nd = 1,sol%nelnds(el)
                IF (sol%ect(nd,el) == n1 .OR. sol%ect(nd,el) == n2) THEN
                  found = found + 1
                ENDIF
              ENDDO
              IF (found == 2) THEN
                sol%bndel(el) = 1
                EXIT elem2
              ENDIF
            ENDDO elem2
          ENDDO
      ENDDO
      ENDIF


      RETURN
      END SUBROUTINE read_nodes
      
      END MODULE read_grid
