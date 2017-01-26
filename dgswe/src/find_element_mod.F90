      MODULE find_element

      USE globals, ONLY: rp
      USE kdtree2_module
      USE lapack_interfaces
      USE sort_mod
      
      IMPLICIT NONE
      
      SAVE
      
      INTEGER :: srchdp 
      REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: rsre          
      TYPE(kdtree2), POINTER :: tree_xy      
      TYPE(kdtree2_result), ALLOCATABLE, DIMENSION(:) :: closest   
      INTEGER, DIMENSION(:), ALLOCATABLE :: nnd2el
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: nd2el
      INTEGER, DIMENSION(:), ALLOCATABLE :: nverts
      INTEGER, DIMENSION(:), ALLOCATABLE :: np
      INTEGER, DIMENSION(:), ALLOCATABLE :: nnds      
      REAL(rp), DIMENSION(:), ALLOCATABLE :: elnx,elny
      INTEGER :: mnepn
      INTEGER :: mnnds
      

      CONTAINS
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   

      SUBROUTINE find_element_init(nel_type,nvertex,nporder,nnodes,nn,xy,nepn,epn)
      
      USE transformation, ONLY: init_vandermonde
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: nel_type
      INTEGER, DIMENSION(:), INTENT(IN) :: nvertex
      INTEGER, DIMENSION(:), INTENT(IN) :: nporder
      INTEGER, DIMENSION(:), INTENT(IN) :: nnodes      
      INTEGER, INTENT(IN) :: nn
      REAL(rp), DIMENSION(:,:), INTENT(IN) :: xy      
      INTEGER, DIMENSION(:), INTENT(IN) :: nepn
      INTEGER, DIMENSION(:,:), INTENT(IN) :: epn
            
      INTEGER :: nd,el,et,n,p
      INTEGER :: mnp
      INTEGER :: info
      
      ! initialize module varibles to increase flexibility and decrease number
      ! of subroutine calling arguments
      
      srchdp = 20
      IF (nn < srchdp) THEN
        srchdp = nn
      ENDIF
      
      tree_xy => kdtree2_create(xy(1:2,1:nn), rearrange=.true., sort=.true.)
      ALLOCATE(closest(MAX(nn,srchdp)))      
            
            
      CALL ref_elem_coords(nel_type,rsre)
            
            
      mnepn = MAXVAL(nepn)
      ALLOCATE(nnd2el(nn),nd2el(mnepn,nn))      
      DO nd = 1,nn
        nnd2el(nd) = nepn(nd)
        DO el = 1,nnd2el(nd)
          nd2el(el,nd) = epn(el,nd)
        ENDDO      
      ENDDO
  
  
      ALLOCATE(nverts(nel_type),np(nel_type),nnds(nel_type))      
      DO et = 1,nel_type
        nverts(et) = nvertex(et)
        np(et) = nporder(et)
        nnds(et) = nnodes(et)
      ENDDO
      
      
      mnp = maxval(np)            
      mnnds = (mnp+1)**2
      ALLOCATE(elnx(mnnds),elny(mnnds))
      
      CALL init_vandermonde(nel_type,np)            
          
      
      RETURN
      END SUBROUTINE find_element_init
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      

      SUBROUTINE in_element(xy,el_type,elxy,el_found,rs,closest_eds,closest_vertex,closest_els)

      IMPLICIT NONE
              
      REAL(rp), INTENT(IN) :: xy(2)
      INTEGER, DIMENSION(:), INTENT(IN) :: el_type
      REAL(rp), DIMENSION(:,:,:), INTENT(IN) :: elxy
      INTEGER, INTENT(OUT) :: el_found          
      REAL(rp), INTENT(OUT) :: rs(2)
      INTEGER, INTENT(OUT), OPTIONAL :: closest_eds(4)      
      INTEGER, INTENT(OUT), OPTIONAL :: closest_vertex 
      INTEGER, INTENT(OUT), DIMENSION(:), ALLOCATABLE, OPTIONAL :: closest_els
     
            
      INTEGER :: srch,nd,i,j
      INTEGER :: el,eln,clnd,et,n
      INTEGER :: found,etemp     
      INTEGER :: min_el
      INTEGER :: vert,leds(4)
      INTEGER :: ntested,tested
      INTEGER, ALLOCATABLE, DIMENSION(:) :: el_tested
      REAL(rp) :: diff,min_diff
      REAL(rp) :: tol,dtemp
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: el_diff
            
      ALLOCATE(el_tested(mnepn*srchdp))
      ALLOCATE(el_diff(mnepn*srchdp))
      
      tol = 1d-5 
      CALL kdtree2_n_nearest(tp=tree_xy,qv=xy,nn=srchdp,results=closest) ! find what element xy is in               
        
      ! Test elements to see which element point is located in    
      ntested = 0
      found = 0      
      min_diff = 999d0
search: DO srch = 1,srchdp

          clnd = closest(srch)%idx          
          
!           PRINT("(A,I9)"), "   closest node: ", clnd

    elem: DO el = 1,nnd2el(clnd) 
 
            eln = nd2el(el,clnd)
            
            tested = 0
     check: DO i = 1,ntested
              IF (el_tested(i) == eln) THEN
                tested = 1      
                EXIT check
              ENDIF
            ENDDO check
            
            IF (tested == 1) THEN
              CYCLE elem
            ELSE
              ntested = ntested + 1
              el_tested(ntested) = eln
            ENDIF
            
            
!             PRINT("(A,I9)"), "     testing element: ", eln                               

            ! Compute sum of sub-triangle areas
            CALL sub_element(xy,eln,el_type,elxy,diff,leds,vert,rs)            
            
            IF (diff < min_diff) THEN  ! keep track of element with minimum difference in sum of sub-triangle areas                                  
              min_diff = diff          ! to return if tolerance is not met
              min_el = eln
            ENDIF
            
            el_diff(ntested) = diff
          
            ! The station is in the element if the reference element area and sum of sub triangle are the same
            IF (diff < tol) THEN
!               PRINT("(A,I9)"), "   element found", eln                            
                      
              el_found = eln        
              found = 1                        
              IF (.NOT. PRESENT(closest_els)) THEN ! exit when found unless collecting all closest elements
                EXIT search            
              ENDIF
            ENDIF    
            
!             PRINT*, " "  
          ENDDO elem
        
        ENDDO search    
        

        
        IF (found == 0) THEN
          el_found = min_el        

          PRINT*, "ELEMENT NOT FOUND"  
          PRINT*, "USING ELEMENT ", el_found, "(AREA DIFF = ",min_diff, ")"       

          CALL sub_element(xy,el_found,el_type,elxy,diff,leds,vert,rs)           
          
        ENDIF    
        
        IF (PRESENT(closest_vertex) ) THEN
          closest_vertex = vert
        ENDIF
 
        IF (PRESENT(closest_eds)) THEN
          closest_eds = leds
        ENDIF 
        
        IF (PRESENT(closest_els)) THEN        
          ALLOCATE(closest_els(ntested))                    
          CALL insertion_sort(ntested,el_diff,el_tested)          
          
          DO i = 1,ntested
            closest_els(i) = el_tested(i)
!             PRINT*, el_diff(i),el_tested(i)
          ENDDO                    
        ENDIF

      RETURN
      END SUBROUTINE in_element
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 

      SUBROUTINE sub_element(xy,eln,el_type,elxy,diff,closest_ed,vert,rs)   
      
      USE transformation, ONLY: xy2rs
      
      IMPLICIT NONE
      
      REAL(rp), INTENT(IN) :: xy(2)        
      INTEGER, INTENT(IN) :: eln
      INTEGER, DIMENSION(:), INTENT(IN) :: el_type
      REAL(rp), DIMENSION(:,:,:), INTENT(IN) :: elxy          
      REAL(rp), INTENT(OUT) :: diff      
      INTEGER, INTENT(OUT) :: closest_ed(4)
      INTEGER, INTENT(OUT) :: vert
      REAL(rp), INTENT(OUT) :: rs(2)
      
      INTEGER :: i,j    
      INTEGER :: nv,et,n,npt,nd,p
      INTEGER :: n1,n2
      INTEGER :: etemp
      REAL(rp) :: area,area_sum
      REAL(rp) :: stri_area(4),stri_min
      REAL(rp) :: dist
      REAL(rp) :: x(1),y(1)
      REAL(rp) :: r(1),s(1)
      REAL(rp) :: x1,x2,x3,y1,y2,y3
      REAL(rp) :: atemp
      
      
      et = el_type(eln)
      nv = nverts(et)
      p = np(et)
      n = nnds(et)
      npt = 1

      
      DO nd = 1,n
        elnx(nd) = elxy(nd,eln,1)
        elny(nd) = elxy(nd,eln,2)
      ENDDO
      
      x(1) = xy(1)
      y(1) = xy(2)
                                                        
      ! Compute the local (r,s) coordinates of the (x,y) point
      CALL xy2rs(et,p,elnx,elny,npt,x,y,r,s)
      
      rs(1) = r(1)
      rs(2) = s(1)
          
      ! Find reference element area
      IF (mod(et,2) == 1) THEN
        area = 2d0
      ELSE IF (mod(et,2) == 0) THEN
        area = 4d0
      ENDIF          
          
      ! Compute sum of sub-triangle areas
      area_sum = 0d0
      stri_min = 999d0
      DO i = 1,nv
        n1 = mod(i+0,nv)+1
        n2 = mod(i+1,nv)+1           

        x1 = rsre(1,n1,et)
        y1 = rsre(2,n1,et)
            
        x2 = rsre(1,n2,et)
        y2 = rsre(2,n2,et)
            
        x3 = r(1)
        y3 = s(1)
                        
        stri_area(i) = .5d0*abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1))
        closest_ed(i) = i
                                            
        area_sum = area_sum + stri_area(i)
      ENDDO
          
!     PRINT("(A,F20.15,A,F20.15)"), "   area = ",area, "   area sum = ", area_sum
!     PRINT*, " "
          
          
          
          
          
      diff = abs(area-area_sum) 
      
!       DO i = 1,nv           ! keep track of minimum sub-triangle area to determine
!         DO j = i+1,nv       ! which edge the point lies on, or is closest to
!           IF (stri_area(j) < stri_area(i)) THEN
!             atemp = stri_area(i)
!             stri_area(i) = stri_area(j)
!             stri_area(j) = atemp
!             
!             etemp = closest_ed(i)
!             closest_ed(i) = closest_ed(j)
!             closest_ed(j) = etemp            
!           ENDIF
!         ENDDO
!       ENDDO
      
      CALL insertion_sort(nv,stri_area,closest_ed) ! keep track of minimum sub-triangle area to determine
                                                   ! which edge the point lies on, or is closest to
      
      vert = 0             ! if one sub-triangle area dominates, then the point is close 
      DO i = 1,nv          ! to that  vertex
        IF (abs(stri_area(i)-area_sum) < 1d-8) THEN
          vert = i
        ENDIF
      ENDDO
            
      
      RETURN
      END SUBROUTINE sub_element    
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   


      SUBROUTINE ref_elem_coords(nel_type,rsre)
      
      USE basis, ONLY: element_nodes
      
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nel_type
      REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: rsre
      
      INTEGER :: et,pt,n
      INTEGER :: np,space
      REAL(rp) :: r(4),s(4)
      
      ALLOCATE(rsre(2,4,nel_type))
      
      space = 1
      np = 1
      
      DO et = 1,nel_type      
        
        CALL element_nodes(et,space,np,n,r,s)
        
        DO pt = 1,n
          rsre(1,pt,et) = r(pt)
          rsre(2,pt,et) = s(pt)
        ENDDO
        
      ENDDO
      
      END SUBROUTINE ref_elem_coords
      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   


      SUBROUTINE check_elem(xy,el_found)      
      
      IMPLICIT NONE
      

      REAL(rp), INTENT(IN) :: xy(2)
      INTEGER, INTENT(IN) :: el_found
            
      INTEGER :: el,eln,clnd
      INTEGER :: found      

                  
      CALL kdtree2_n_nearest(tp=tree_xy,qv=xy,nn=1,results=closest)              
        
  
      found = 0      
      clnd = closest(1)%idx          

elem: DO el = 1,nnd2el(clnd) 
 
        eln = nd2el(el,clnd)
            
        IF (eln == el_found) THEN
          found = 1
          EXIT elem
        ENDIF
            
      ENDDO elem
        
        
      IF (found == 0) THEN
        PRINT*, "ERROR FINDING ELEMENT FOR VERTEX NODE"
        STOP
      ENDIF     
 
      
      
      RETURN
      END SUBROUTINE check_elem

      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       

      SUBROUTINE find_element_final()
      
      USE transformation, ONLY: final_vandermonde      
      
      IMPLICIT NONE
      
      NULLIFY(tree_xy)
      DEALLOCATE(closest)
      DEALLOCATE(rsre)
      DEALLOCATE(nnd2el,nd2el)
      DEALLOCATE(nverts,np,nnds)
      DEALLOCATE(elnx,elny)
      
      CALL final_vandermonde
      
      RETURN
      END SUBROUTINE find_element_final

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
      
      END MODULE find_element