      SUBROUTINE find_nesting(coarse,fine,elf2elc)

      USE globals, ONLY:pres,solution
      USE evaluate, ONLY: newton
      USE kdtree2_module

      IMPLICIT NONE
      
      TYPE(solution) :: coarse,fine
      TYPE(kdtree2), POINTER :: tree
      TYPE(kdtree2_result), ALLOCATABLE, DIMENSION(:) :: kdresults
      
      INTEGER :: i,elf,elc,nd,srch
      INTEGER :: etf,etc,nvert,clnd,eln
      INTEGER :: n1,n2
      INTEGER :: el_found,found
      INTEGER :: srchdp
      INTEGER :: elf2elc(fine%ne)
      
      REAL(pres) :: area,sarea,tol
      REAL(pres) :: x(4),y(4)
      REAL(pres) :: xf(1),yf(1),hb(1)     
      REAL(pres) :: r(1),s(1)
      
      srchdp = 4
      tol = 1d-7
      
      

      ! Build kd-tree
      IF (coarse%ne < srchdp) THEN
        srchdp = coarse%ne
      ENDIF      
      ALLOCATE(kdresults(srchdp))            
      tree => kdtree2_create(coarse%vxy, rearrange=.true., sort=.true.)
      
      
      
      
      elf2elc(:) = 0
      el_found = 0 
      
      PRINT("(A)"), TRIM(coarse%sol_name) // " / " // TRIM(fine%sol_name) // " grid"
      
      ! Find coarse element each fine element is located in
elemf:DO elf = 1,fine%ne
      
        IF (fine%bndel(elf) == 1) THEN ! ignore elements on the boundary
          el_found = el_found + 1        
          CYCLE elemf
        ENDIF
      
        etf = fine%el_type(elf)
        
        ! Compute the (x,y) coordinates of the first quadrature point (from fine element)
        xf(1) = 0d0
        yf(1) = 0d0
        DO nd = 1,fine%nnds(etf)
          xf(1) = xf(1) + fine%l(nd,1,etf)*fine%elxy(nd,elf,1)
          yf(1) = yf(1) + fine%l(nd,1,etf)*fine%elxy(nd,elf,2)
        ENDDO
      
        ! Find coarse node closest to quadrature point
        CALL kdtree2_n_nearest(tp=tree,qv=(/xf(1),yf(1)/),nn=srchdp,results=kdresults)
        
search: DO srch = 1,srchdp
          clnd = coarse%vxyn(kdresults(srch)%idx)      
        
!           PRINT("(A,I5,A,I5,A,F11.4)"), "fine element #: ",elf, ",   closest coarse node: ", clnd, ",   distance: ", kdresults(1)%dis      
        
          ! Test elements associated with nearest node to see which coarse element contains the quadrature point
          found = 0  
    elem: DO elc = 1,coarse%nepn(clnd)
            eln = coarse%epn(elc,clnd)
          
            etc = coarse%el_type(eln)
            nvert = coarse%nverts(etc)                
          
            ! Compute the local (r,s) coarse element coordinates of the (x,y) quadrature point
            CALL newton(coarse,xf(1),yf(1),1,eln,r,s,hb)
          
            ! Find reference element area
            IF (mod(etc,2) == 1) THEN
              area = 2d0
            ELSE IF (mod(etc,2) == 0) THEN
              area = 4d0
            ENDIF          
          
            ! Compute sum of sub-triangle areas
            sarea = 0d0
            DO i = 1,nvert
              n1 = mod(i+0,nvert)+1
              n2 = mod(i+1,nvert)+1           

              x(1) = coarse%rsre(1,n1,etc)
              y(1) = coarse%rsre(2,n1,etc)
            
              x(2) = coarse%rsre(1,n2,etc)
              y(2) = coarse% rsre(2,n2,etc)
            
              x(3) = r(1)
              y(3) = s(1)
            
              sarea = sarea + .5d0*abs((x(2)-x(1))*(y(3)-y(1)) - (x(3)-x(1))*(y(2)-y(1)))
            ENDDO
          
!             PRINT("(A,I5,A,F20.15,A,F20.15)"), "   testing: ", eln, "   area = ",area, "   sarea = ", sarea
!             PRINT*, " "
          
            ! The quadrature point is in the element if the reference element area and sum of sub triangle are the same
            IF (abs(area - sarea) < tol) THEN          
!               PRINT("(A,I5)"), achar(27)//'[95m    element found '//achar(27)//'[0m', eln
            
              elf2elc(elf) = eln            
            
              found = 1
              el_found = el_found + 1
                                    
              EXIT elem                        
            ENDIF
          
          ENDDO elem
        
          IF (found == 0) THEN
!             PRINT*, "element not found: going to next node"    ! Try elements around next closest node     
          ELSE   
            EXIT search
          ENDIF
        
!           PRINT*, " " 
        
        ENDDO search
        
        IF (found == 0) THEN
!           PRINT*,  achar(27)//"[31m ERROR: ELEMENT NOT FOUND"//achar(27)//'[0m'
!           STOP
        ENDIF        
      
      ENDDO elemf
      
      PRINT("(A,I5)"), "Missing elements = ", fine%ne-el_found
      PRINT*, " "
      
      IF (fine%ne-el_found /= 0) THEN
        STOP
      ENDIF

      
      
      
      RETURN
      END SUBROUTINE find_nesting