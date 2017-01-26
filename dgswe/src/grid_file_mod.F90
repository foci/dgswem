      MODULE grid_file_mod

      USE globals, ONLY: rp,pi,g
      USE quit, ONLY: abort    
   
      
      IMPLICIT NONE
         

      CONTAINS
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      SUBROUTINE read_header(myrank,grid_file,grid_name,ne,nn)

      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: myrank
      CHARACTER(100), INTENT(IN)  :: grid_file
      CHARACTER(100), INTENT(OUT) :: grid_name      
      INTEGER, INTENT(OUT) :: ne
      INTEGER, INTENT(OUT) :: nn   
      
      LOGICAL :: file_exists          
      

      ! open fort.14 grid file
      INQUIRE(FILE=grid_file, EXIST=file_exists)
      IF(file_exists == .FALSE.) THEN
        PRINT*, "grid file does not exist"
        CALL abort()        
      ENDIF
      
      
      IF (myrank == 0 ) PRINT "(A)", "Reading in grid file"       
      
      OPEN(UNIT=14, FILE=grid_file)                 
                       
      ! read in name of grid
      READ(14,"(A)"), grid_name                                         

      ! read in number of elements and number of nodes
      READ(14,*), ne, nn       

      RETURN
      END SUBROUTINE read_header
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  


      SUBROUTINE read_coords(nn,xy,depth,h0)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: nn
      REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: xy
      REAL(rp), DIMENSION(:)  , ALLOCATABLE, INTENT(OUT) :: depth
      REAL(rp), INTENT(IN), OPTIONAL :: h0

      
      INTEGER :: i,j
      INTEGER :: limit_depth
      INTEGER :: alloc_status           
      alloc_status = 0    

      ALLOCATE(xy(2,nn),depth(nn), STAT=alloc_status)
      IF (alloc_status /= 0) THEN
        PRINT*, "Allocation error"
        CALL abort()
      ENDIF
      
      ! read in node coordinates and depths
      DO i = 1,nn                                                      
        READ(14,*), j, xy(1,j), xy(2,j), depth(j)        
      ENDDO
      
      
      limit_depth = 0
      IF (PRESENT(h0)) THEN
        limit_depth = 1
      ENDIF      
      
      IF (limit_depth) THEN
        DO i = 1,nn      
          IF (depth(i) < h0) THEN
            depth(i) = h0
          ENDIF      
        ENDDO
      ENDIF
      
!       PRINT "(A)", "Node coordinates and depth: "
!       DO i = 1,nn
!         PRINT "(I5,3(F11.3,3x))", i,xy(1,i), xy(2,i), depth(i)
!       ENDDO
!       PRINT*, " "

 
      
      RETURN
      END SUBROUTINE read_coords
      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  


      SUBROUTINE read_connectivity(ne,ect,el_type)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: ne
      INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: ect
      INTEGER, DIMENSION(:)  , ALLOCATABLE, INTENT(OUT) :: el_type

      INTEGER :: i,j
      INTEGER :: el
      INTEGER :: nv
      INTEGER :: alloc_status       
      
      
      ALLOCATE(ect(4,ne),el_type(ne), STAT=alloc_status)
      IF (alloc_status /= 0) THEN
        PRINT*, "Allocation error"
        CALL abort()
      ENDIF      
      
      ! read in element connectivity
      DO i = 1,ne
        READ(14,*) el,nv,(ect(j,el),j = 1,nv)
        
        IF (nv == 3) THEN
          el_type(el) = 1
        ELSE IF (nv == 4) THEN
          el_type(el) = 2
        ELSE
          PRINT*, "Element type not supported"
          CALL abort()
        ENDIF 

      ENDDO           
      
!       PRINT "(A)", "Element connectivity table: "
!       DO i = 1,ne
!         PRINT "(2(I5,3x),8x,4(I5,3x))", i,nelnds(i),(ect(j,i),j=1,nelnds(i))
!       ENDDO
!       PRINT*, " "      
      
      RETURN
      END SUBROUTINE read_connectivity 
      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  


      SUBROUTINE read_open_boundaries(nope,neta,obseg,obnds)
      
      IMPLICIT NONE     
      
      INTEGER, INTENT(OUT) :: nope
      INTEGER, INTENT(OUT) :: neta
      INTEGER, DIMENSION(:)  , ALLOCATABLE, INTENT(OUT) :: obseg      
      INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: obnds
      
      INTEGER :: i,j
      INTEGER :: nbseg
      INTEGER :: alloc_status       
      
      
      
      READ(14,*) nope  ! number of open boundaries                                                 
      READ(14,*) neta  ! number of total elevation specified boundary nodes
      
      ALLOCATE(obseg(nope),obnds(neta,nope) ,STAT=alloc_status)
      IF (alloc_status /= 0) THEN
        PRINT*, "Allocation error"
        CALL abort()
      ENDIF          

      DO i = 1,nope                                                     
        READ(14,*), nbseg  ! read in # of nodes in segment, boundary type
        obseg(i) = nbseg
        DO j = 1,nbseg
          READ(14,*) obnds(j,i) ! read in open boundary node numbers
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
      
      RETURN
      END SUBROUTINE read_open_boundaries

      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      SUBROUTINE read_flow_boundaries(nbou,nvel,fbseg,fbnds)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(OUT) :: nbou
      INTEGER, INTENT(OUT) :: nvel
      INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: fbseg
      INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: fbnds
            
      INTEGER :: i,j      
      INTEGER :: nbseg,btype
      INTEGER :: alloc_status       
      
      
      READ(14,*) nbou  ! number of normal flow boundaries
      READ(14,*) nvel  ! total number of normal flow nodes

      ALLOCATE(fbseg(2,nbou),fbnds(nvel,nbou), STAT=alloc_status)
      IF (alloc_status /= 0) THEN
        PRINT*, "Allocation error"
        CALL abort()
      ENDIF        

      DO i = 1,nbou
        READ(14,*), nbseg, btype ! read in # of nodes in segment, boundary type
        fbseg(1,i) = nbseg
        fbseg(2,i) = btype
        DO j = 1,nbseg
          READ(14,*), fbnds(j,i)  ! read in normal flow boundary node numbers
        ENDDO
        IF (btype == 1 .OR. btype == 11 .OR. btype == 21) THEN
          IF (fbnds(nbseg,i) /= fbnds(1,i)) THEN
            fbnds(nbseg+1,i) = fbnds(1,i)  ! close island boundaries
            fbseg(1,i) = fbseg(1,i) + 1
          ENDIF
        ENDIF
      ENDDO
      
      CLOSE(14)        
      
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
      
      RETURN
      END SUBROUTINE read_flow_boundaries
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      

      SUBROUTINE copy_footer(grid_file_from,grid_file_to)
      
      IMPLICIT NONE
      
      CHARACTER(*), INTENT(IN) :: grid_file_from
      CHARACTER(*), INTENT(IN) :: grid_file_to           
      
      CHARACTER(200) :: line
      CHARACTER(10) :: flag
      LOGICAL :: file_exists
      INTEGER :: read_stat
      
      
      INQUIRE(FILE=grid_file_from, EXIST=file_exists)
      IF(file_exists == .FALSE.) THEN
        PRINT*, "grid file does not exist"
        CALL abort()        
      ENDIF   
      
      OPEN(UNIT=141, FILE=grid_file_from)  
      
      DO 
      
        READ(141,*,IOSTAT=read_stat) flag
        IF (flag == "!!!!!!!!!!") THEN
          EXIT
        ELSE IF (read_stat < 0) THEN
          PRINT*, "footer does not exist"
          RETURN
        ENDIF
        
      ENDDO
      
      INQUIRE(FILE=grid_file_to, EXIST=file_exists)
      IF(file_exists == .FALSE.) THEN
        PRINT*, "grid file does not exist"
        CALL abort()        
      ENDIF        
      
      OPEN(UNIT=142, FILE=grid_file_to) 
      
      DO 
      
        READ(142,*,IOSTAT=read_stat) flag
        IF (flag == "!!!!!!!!!!") THEN
          EXIT
        ELSE IF (read_stat < 0) THEN
          PRINT*, "exiting to"
          RETURN
        ENDIF
        
      ENDDO      
           
      DO 
      
        READ(141,"(A)",IOSTAT=read_stat) line
        IF (read_stat < 0) THEN
          EXIT
        ENDIF
      
        WRITE(142,"(A)") line

      ENDDO      
      
      CLOSE(141)
      CLOSE(142)
      
      RETURN
      END SUBROUTINE copy_footer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
      
      SUBROUTINE init_element_coordinates(ne,ctp,el_type,nverts,xy,ect,elxy)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: ne
      INTEGER, INTENT(IN) :: ctp
      INTEGER, DIMENSION(:), INTENT(IN) :: el_type
      INTEGER, DIMENSION(:), INTENT(IN) :: nverts
      REAL(rp), DIMENSION(:,:), INTENT(IN) :: xy
      INTEGER, DIMENSION(:,:), INTENT(IN) :: ect
      REAL(rp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(OUT) :: elxy
            
      INTEGER :: el,nd
      INTEGER :: et,nv
      INTEGER :: mnnds
      INTEGER :: alloc_status
      
      mnnds = (ctp+1)**2
      
      ALLOCATE(elxy(mnnds,ne,2), STAT=alloc_status)
      
      
      DO el = 1,ne        
        et = el_type(el)
        nv = nverts(et)
        DO nd = 1,nv
          elxy(nd,el,1) = xy(1,ect(nd,el))
          elxy(nd,el,2) = xy(2,ect(nd,el))
        ENDDO              
      ENDDO          
      
      RETURN
      END SUBROUTINE init_element_coordinates      

      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE init_element_depth(ne,hbp,el_type,nverts,depth,ect,elhb)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: ne
      INTEGER, INTENT(IN) :: hbp
      INTEGER, DIMENSION(:), INTENT(IN) :: el_type
      INTEGER, DIMENSION(:), INTENT(IN) :: nverts
      REAL(rp), DIMENSION(:), INTENT(IN) :: depth
      INTEGER, DIMENSION(:,:), INTENT(IN) :: ect
      REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: elhb
            
      INTEGER :: el,nd
      INTEGER :: et,nv
      INTEGER :: mnnds
      INTEGER :: alloc_status
      
      mnnds = (hbp+1)**2
      
      ALLOCATE(elhb(mnnds,ne), STAT=alloc_status)              
      
      DO el = 1,ne        
        et = el_type(el)
        nv = nverts(et)
        DO nd = 1,nv
          elhb(nd,el) = depth(ect(nd,el))
        ENDDO              
      ENDDO        
      
      RETURN
      END SUBROUTINE init_element_depth      

      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      SUBROUTINE read_bathy_file(myrank,bathy_file,hbp,ne,el_type,nverts,depth,ect,elhb,hb_file_exists)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: myrank
      CHARACTER(100), INTENT(IN) :: bathy_file
      INTEGER, INTENT(IN) :: hbp
      INTEGER, INTENT(IN) :: ne
      INTEGER, DIMENSION(:), INTENT(IN) :: el_type
      INTEGER, DIMENSION(:), INTENT(IN) :: nverts      
      REAL(rp), DIMENSION(:), INTENT(IN) :: depth
      INTEGER, DIMENSION(:,:), INTENT(IN) :: ect
      REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: elhb      
      LOGICAL, INTENT(OUT), OPTIONAL :: hb_file_exists
      
      INTEGER :: i
      INTEGER :: el,nd
      INTEGER :: et,nv
      INTEGER :: nnds,mnnds
      INTEGER :: ne_check
      INTEGER :: hbp_check
      LOGICAL :: file_exists  
      INTEGER :: alloc_status   
      CHARACTER(200) :: tmp      
      CHARACTER(10) :: flag      
              
      mnnds = (hbp+1)**2
      
      INQUIRE(FILE=bathy_file, EXIST = file_exists)
      IF(file_exists == .FALSE.) THEN
        IF (myrank == 0) PRINT*, "high order bathymetry file does not exist"  
        
        IF (hbp > 1) THEN
          IF(myrank == 0) PRINT*, "high order bathymetry file is required for hbp > 1"
          CALL abort()   
        ELSE  
        
          
          CALL init_element_depth(ne,hbp,el_type,nverts,depth,ect,elhb)          
          
        ENDIF
      ELSE
      
      
        IF (myrank == 0 ) PRINT "(A)", "Reading in high order bathymetry file"  
        
        OPEN(UNIT = 14, FILE = bathy_file)
        
    pre:DO                                !read in the header info
          READ(14,*) flag
          IF (flag == "!!!!!!!!!!") THEN
            EXIT pre
          ENDIF
        ENDDO pre  
        
        READ(14,*) ne_check, hbp_check
        IF (ne_check /= ne .or. hbp_check /= hbp) THEN
          IF (myrank == 0) PRINT*, "incorrect high order bathymetry file"
          CALL abort()
        ENDIF
        
        ALLOCATE(elhb(mnnds,ne), STAT=alloc_status)
        IF (alloc_status /= 0) THEN
          PRINT*, "Allocation error"
          CALL abort()
        ENDIF         
        
        DO i = 1,ne
          READ(14,*) el,nnds,(elhb(nd,el), nd = 1,nnds)
        ENDDO
      
        CLOSE(14)
      ENDIF      
      
      IF (PRESENT(hb_file_exists)) THEN
        hb_file_exists = file_exists
      ENDIF      
      
      
      RETURN 
      END SUBROUTINE read_bathy_file


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      SUBROUTINE read_curve_file(myrank,curve_file,ctp,nbou,xy,bndxy,cb_file_exists)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: myrank
      CHARACTER(100), INTENT(IN) :: curve_file
      INTEGER, INTENT(IN) :: ctp
      INTEGER, INTENT(IN) :: nbou
      REAL(rp), DIMENSION(:,:), INTENT(OUT) :: xy
      REAL(rp), DIMENSION(:,:,:,:), ALLOCATABLE, INTENT(OUT) :: bndxy
      LOGICAL, INTENT(OUT), OPTIONAL :: cb_file_exists 

      INTEGER :: i,j,k
      INTEGER :: nd
      INTEGER :: nbou_check
      INTEGER :: ctp_check
      INTEGER :: nmax
      INTEGER :: nbseg,btype      
      INTEGER :: alloc_status   
      LOGICAL :: file_exists    
      CHARACTER(10) :: flag      
      CHARACTER(200) :: tmp        
      
      IF (ctp == 1) THEN
        RETURN
      ENDIF
      
      
      INQUIRE(FILE=curve_file, EXIST = file_exists)  
      IF (file_exists == .FALSE.) THEN
        IF (myrank == 0) PRINT*, "curved boundary file does not exist"
        
        IF (ctp > 1) THEN
          IF(myrank == 0) PRINT*, "curved boundary file is required for ctp > 1"
          CALL abort()          
        ENDIF        
        
      ELSE
      
      
        IF (myrank == 0) PRINT "(A)", "Reading in curved boundary file"
      
        OPEN(UNIT=14, FILE=curve_file)
        
    pre:DO                                !read in the header info
          READ(14,*) flag
          IF (flag == "!!!!!!!!!!") THEN
            EXIT pre
          ENDIF
        ENDDO pre  
        
        READ(14,*) nbou_check
        READ(14,*) nmax,ctp_check
        
        IF (ctp_check /= ctp) THEN
          PRINT*, "incorrect curved boundary file"
          CALL abort()
        ENDIF
        
!         IF (nbou_check /= nbou) THEN
!           PRINT*, "incorrect curved boundary file"
!           CALL abort()
!         ENDIF        
        
        ALLOCATE(bndxy(2,ctp+1,nmax,nbou_check), STAT = alloc_status)
        IF (alloc_status /= 0) THEN
          PRINT*, "Allocation error"
          CALL abort()
        ENDIF          
        
        DO i = 1,nbou_check
          READ(14,*) nbseg,btype  
          IF(nbseg > 0) THEN
            DO j = 1,nbseg-1
              READ(14,*) nd, xy(1,nd),(bndxy(1,k,j,i), k=1,ctp-1)
              READ(14,*) nd, xy(2,nd),(bndxy(2,k,j,i), k=1,ctp-1)
            ENDDO
            READ(14,*) nd, xy(1,nd)
            READ(14,*) nd, xy(2,nd)
          ENDIF
        ENDDO
        
        CLOSE(14)
      ENDIF      
      
      
      IF (PRESENT(cb_file_exists)) THEN
        cb_file_exists = file_exists
      ENDIF
  
      
      RETURN 
      END SUBROUTINE read_curve_file


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  
      SUBROUTINE print_grid_info(grid_file,grid_name,ne,nn)
      
      IMPLICIT NONE

      CHARACTER(100), INTENT(IN) :: grid_file
      CHARACTER(100), INTENT(IN) :: grid_name
      INTEGER, INTENT(IN) :: ne
      INTEGER, INTENT(IN) :: nn
      
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", "             Grid Information                "
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", " "
      PRINT "(A,A)", "Grid file: ", grid_file                     
      PRINT "(A,A)", "Grid name: ", grid_name      
      PRINT "(A,I15)", "Number of elements: ", ne
      PRINT "(A,I15)", "Number of nodes: ", nn
      PRINT*, " "          
      
      
      RETURN
      END SUBROUTINE print_grid_info


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 


      SUBROUTINE grid_size(ne,el_type,ect,xy,h)
      
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: ne
      INTEGER, DIMENSION(:), INTENT(IN) :: el_type
      INTEGER, DIMENSION(:,:), INTENT(IN) :: ect
      REAL(rp), DIMENSION(:,:), INTENT(IN) :: xy
      REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: h
      
      INTEGER :: el,i,et
      INTEGER :: n1,n2,nv
      REAL(rp) :: x(4),y(4)
      REAL(rp) :: psml,pl
      REAL(rp) :: alpha,gamma
      REAL(rp) :: area,s,r,l(4)

      
      ALLOCATE(h(ne))
                 
      DO el = 1,ne
       
        et = el_type(el) 
        
        IF (mod(et,2) == 1) THEN
          nv = 3
        ELSE IF (mod(et,2) == 0) THEN  
          nv = 4
        ENDIF
        
        ! calculate edge lengths
        DO i = 1,nv
          n1 = mod(i+0,nv) + 1
          n2 = mod(i+1,nv) + 1
            
          x(n1) = xy(1,ect(n1,el))
          x(n2) = xy(1,ect(n2,el))
            
          y(n1) = xy(2,ect(n1,el))
          y(n2) = xy(2,ect(n2,el))
            
          l(i) = sqrt((x(n2)-x(n1))**2 + (y(n2)-y(n1))**2)
        ENDDO
        
        ! calculate semiperimeter
        s = 0
        DO i = 1,nv
          s = s + l(i)
        ENDDO
        s = 0.5*s
        
        
        psml = 1d0  ! product of semiperimeter - edge lengths
        pl = 1d0    ! product of edge lengths
        DO i = 1,nv
         psml = psml*(s-l(i))
         pl = pl*l(i)
        ENDDO
        
        IF (mod(et,2) == 1) THEN
          
          ! calculate radius of smallest inscribed circle for triangles
          r = sqrt(psml/s)       
                               
        ELSE IF (mod(et,2) == 0) THEN  
          
          !calculate radius of circle with equal area for quadrilaterals
          alpha = angle(x(4),y(4),x(1),y(1),x(2),y(2))
          gamma = angle(x(2),y(2),x(3),y(3),x(4),y(4))
          area = sqrt(psml-.5d0*pl*(1d0+cos(alpha+gamma))) ! Bretschneider's formula          
          
          r = sqrt(area/pi)
          
        ENDIF
        
        h(el) = 2d0*r
      ENDDO

      RETURN
      END SUBROUTINE grid_size
      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      REAL(rp) FUNCTION angle(x1,y1,x2,y2,x3,y3)
      
      IMPLICIT NONE  
      
      REAL(rp) :: x1,x2,x3,xy,y1,y2,y3,y4
      REAL(rp) :: xv1,xv2,yv1,yv2
      REAL(rp) :: adotb,la,lb
      
      
      xv1 = x1-x2
      yv1 = y1-y2
      
      xv2 = x3-x2
      yv2 = y3-y2      
      
      adotb = xv1*xv2 + yv1*yv2
      
      la = sqrt(xv1**2 + yv1**2)
      lb = sqrt(xv2**2 + yv2**2)
      
      angle = acos(adotb/(la*lb))*(180d0/pi)
      
      RETURN
      END FUNCTION
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

     
      SUBROUTINE courant(p,ne,u,cfl,el_type,nverts,nnds,elhb,h)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: p
      INTEGER, INTENT(IN) :: ne
      REAL(rp), INTENT(IN) :: u
      REAL(rp), INTENT(IN) :: cfl      
      INTEGER, DIMENSION(:), INTENT(IN) :: el_type
      INTEGER, DIMENSION(:), INTENT(IN) :: nverts
      INTEGER, DIMENSION(:), INTENT(IN) :: nnds
      REAL(rp), DIMENSION(:,:), INTENT(IN) :: elhb
      REAL(rp), DIMENSION(:), INTENT(IN) :: h
      
      INTEGER :: el,nd
      INTEGER :: et,nv,nnd
      REAL(rp) :: hb
      REAL(rp) :: c
      REAL(rp) :: dt,dtmin
      
      dtmin = 9999d0
      
      DO el = 1,ne
        et = el_type(el)
        nv = nverts(et)
        
        IF (mod(et,2) == 1) THEN
          nnd = nnds(5)
        ELSE IF (mod(et,2) == 0) THEN
          nnd = nnds(6)
        ENDIF

        hb = 0d0
        DO nd = 1,nnd
          IF (elhb(nd,el) > hb) THEN
            hb = elhb(nd,el)
          ENDIF
        ENDDO

        
        c = u + sqrt(g*hb)
       
        dt = h(el)/c
        dt = dt/(2d0*real(p,rp)+1d0)
        
        IF (dt < dtmin) THEN
          dtmin = dt
        ENDIF
          
      ENDDO
      
      PRINT ("(A,F6.3,A,F6.3,A,F18.3,A)"), "Rough max timestep based on CFL = ",cfl, " and u = ", u, ": ", dtmin, " sec"
      PRINT*,""
      
      RETURN
      END SUBROUTINE

      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      SUBROUTINE write_header(grid_file,grid_name,ne,nn)

      IMPLICIT NONE
      
      CHARACTER(*), INTENT(IN)  :: grid_file
      CHARACTER(*), INTENT(IN) :: grid_name      
      INTEGER, INTENT(IN) :: ne
      INTEGER, INTENT(IN) :: nn   
      
      ! open fort.14 grid file      
      OPEN(UNIT=14, FILE=TRIM(grid_file))                 
                       
      ! write out name of grid
      WRITE(14,"(A)"), TRIM(grid_name)                                         

      ! read in number of elements and number of nodes
      WRITE(14,"(2(I8,1x))"), ne, nn       

      RETURN
      END SUBROUTINE write_header
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  


      SUBROUTINE write_coords(nn,xy,depth)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: nn
      REAL(rp), DIMENSION(:,:), INTENT(IN) :: xy
      REAL(rp), DIMENSION(:)  , INTENT(IN) :: depth
      
      INTEGER :: i

      ! write out node coordinates and depths
      DO i = 1,nn                                                      
!         WRITE(14,"(I8,1X,3(D24.17,1X))"), i, xy(1,i), xy(2,i), depth(i)
        WRITE(14,"(I8,1X,3(E24.17,1X))"), i, xy(1,i), xy(2,i), depth(i)           
      ENDDO 
      
      RETURN
      END SUBROUTINE write_coords
      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  


      SUBROUTINE write_connectivity(ne,ect,el_type,nverts)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: ne
      INTEGER, DIMENSION(:,:), INTENT(IN) :: ect
      INTEGER, DIMENSION(:)  , INTENT(IN) :: el_type
      INTEGER, DIMENSION(:)  , INTENT(IN) :: nverts

      INTEGER :: i,j
      INTEGER :: et,nv
      
      ! write out element connectivity
      DO i = 1,ne
        et = el_type(i)
        nv = nverts(et)
        WRITE(14,"(6(I8,1x))") i,nv,(ect(j,i),j = 1,nv)        
      ENDDO           
   
      
      RETURN
      END SUBROUTINE write_connectivity 
      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  


      SUBROUTINE write_open_boundaries(nope,neta,obseg,obnds)
      
      IMPLICIT NONE     
      
      INTEGER, INTENT(IN) :: nope
      INTEGER, INTENT(IN) :: neta
      INTEGER, DIMENSION(:)  , INTENT(IN) :: obseg      
      INTEGER, DIMENSION(:,:), INTENT(IN) :: obnds
      
      INTEGER :: i,j
      INTEGER :: nbseg     
                  
      WRITE(14,"(I8,19x,A)") nope, "! number of open boundaries"                                                
      WRITE(14,"(I8,19x,A)") neta, "! number of total elevation specified boundary nodes"

      DO i = 1,nope      
        nbseg = obseg(i)
        WRITE(14,"(I8,10x,A,1x,I8)"), nbseg, "! number of nodes in open boundary", i
        DO j = 1,nbseg
          WRITE(14,"(I8)") obnds(j,i) ! write out open boundary node numbers
        ENDDO
      ENDDO   
      
      RETURN
      END SUBROUTINE write_open_boundaries

      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      SUBROUTINE write_flow_boundaries(nbou,nvel,fbseg,fbnds)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: nbou
      INTEGER, INTENT(IN) :: nvel
      INTEGER, DIMENSION(:,:), INTENT(IN) :: fbseg
      INTEGER, DIMENSION(:,:), INTENT(IN) :: fbnds
            
      INTEGER :: i,j      
      INTEGER :: nbseg,btype      
            
      WRITE(14,"(I8,19x,A)") nbou, "! number of normal flow boundaries"
      WRITE(14,"(I8,19x,A)") nvel, "! total number of normal flow nodes"

      DO i = 1,nbou
        nbseg = fbseg(1,i)
        btype = fbseg(2,i)      
        WRITE(14,"(I8,1x,I8,10x,A,1x,I8)"), nbseg, btype, "! number of nodes in normal flow boundary", i
        DO j = 1,nbseg
          WRITE(14,"(I8)"), fbnds(j,i)  ! read in normal flow boundary node numbers
        ENDDO
      ENDDO

      WRITE(14,"(A)") " "
      WRITE(14,"(A)") "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      
      CLOSE(14)        

      RETURN
      END SUBROUTINE write_flow_boundaries
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      END MODULE grid_file_mod