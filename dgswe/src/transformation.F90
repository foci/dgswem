      MODULE transformation
      
      USE globals, ONLY: rp
      USE lapack_interfaces
      
      IMPLICIT NONE
      
      SAVE 
      
      INTEGER :: ldv
      INTEGER :: vandermonde_init = 0
      REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: V
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: ipiv
      
      CONTAINS
      
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   

      SUBROUTINE init_vandermonde(nel_type,np)
      
      USE vandermonde, ONLY: vandermonde_area
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: nel_type
      INTEGER, DIMENSION(:), INTENT(IN) :: np
      
      INTEGER :: et
      INTEGER :: nnd
      INTEGER :: mp,mnnds
      INTEGER :: info
      
      IF (vandermonde_init == 0) THEN
        mp = MAXVAL(np)
        ldv = (mp+1)**2
        ALLOCATE(V(ldv,ldv,nel_type))
        ALLOCATE(ipiv(ldv,nel_type))
      
        DO et = 1,nel_type
          CALL vandermonde_area(et,np(et),nnd,V(:,:,et))
          CALL DGETRF(nnd,nnd,V(1,1,et),ldv,ipiv(1,et),info)           
        ENDDO
      
        vandermonde_init = 1
      ENDIF
      
      RETURN
      END SUBROUTINE      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   

      SUBROUTINE final_vandermonde()
      
      IMPLICIT NONE
      
      IF (vandermonde_init == 1) THEN
      
        DEALLOCATE(V)
        DEALLOCATE(ipiv)     
      
        vandermonde_init = 0
      ENDIF
      
      RETURN
      END SUBROUTINE            
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      


      SUBROUTINE element_transformation(nnd,x,y,l,xpt,ypt,dldr,dlds,drdx,drdy,dsdx,dsdy,detJ)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nnd

      REAL(rp), DIMENSION(:), INTENT(IN) :: x,y
      REAL(rp), DIMENSION(:), INTENT(IN) :: l
      REAL(rp), INTENT(OUT) :: xpt,ypt
      REAL(rp), DIMENSION(:), INTENT(IN), OPTIONAL :: dldr,dlds
      REAL(rp), INTENT(OUT), OPTIONAL :: drdx,drdy,dsdx,dsdy,detJ

      INTEGER :: nd
      REAL(rp) :: dxdr,dxds,dydr,dyds
      INTEGER :: calc_inv


      IF ( PRESENT(dldr) .AND. PRESENT(dlds) .AND. &
           PRESENT(drdx) .AND. PRESENT(drdy) .AND. &
           PRESENT(dsdx) .AND. PRESENT(dsdy) .AND. &
           PRESENT(detJ) ) THEN
        calc_inv = 1
      ELSE
        calc_inv = 0
      ENDIF
     

      xpt = 0d0
      ypt = 0d0
          
      DO nd = 1,nnd        
              
        xpt = xpt + l(nd)*x(nd)
        ypt = ypt + l(nd)*y(nd)
                      
      ENDDO

         
      IF (calc_inv == 1) THEN

        dxdr = 0d0
        dxds = 0d0
        dydr = 0d0
        dyds = 0d0
            
          
        DO nd = 1,nnd
        
          dxdr = dxdr + dldr(nd)*x(nd)
          dxds = dxds + dlds(nd)*x(nd)
          dydr = dydr + dldr(nd)*y(nd)
          dyds = dyds + dlds(nd)*y(nd)              
                      
        ENDDO
            
        detJ = dxdr*dyds - dxds*dydr
            
        drdx =  dyds/detJ
        drdy = -dxds/detJ
        dsdx = -dydr/detJ
        dsdy =  dxdr/detJ

      ENDIF

      END SUBROUTINE element_transformation


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     

      SUBROUTINE xy2rs(et,p,elx,ely,npt,x,y,r,s)

      USE shape_functions_mod, ONLY: shape_functions_area_eval
      USE basis, ONLY: element_basis

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: et
      INTEGER, INTENT(IN) :: p
      REAL(rp), DIMENSION(:), INTENT(IN) :: elx
      REAL(rp), DIMENSION(:), INTENT(IN) :: ely
      INTEGER, INTENT(IN) :: npt
      REAL(rp), DIMENSION(:), INTENT(IN) :: x   
      REAL(rp), DIMENSION(:), INTENT(IN) :: y       
      REAL(rp), DIMENSION(:), INTENT(OUT) :: r   
      REAL(rp), DIMENSION(:), INTENT(OUT) :: s      
      
      INTEGER :: it,n,pt,m
      INTEGER :: maxit
      INTEGER :: mnnds
      INTEGER :: info
      REAL(rp) :: tol    
      REAL(rp) :: f,g
      REAL(rp) :: error(npt),err
      REAL(rp) :: drdf,drdg,dsdf,dsdg,jac
      REAL(rp), DIMENSION(:,:), ALLOCATABLE :: l,dldr,dlds
      REAL(rp), DIMENSION(:,:), ALLOCATABLE :: phi
        
      tol = 1d-10
      maxit = 1000
          
      mnnds = (p+1)**2      
      ALLOCATE(l(mnnds,npt),dldr(mnnds,npt),dlds(mnnds,npt))
      ALLOCATE(phi(mnnds,3*npt))
        
      ! Initial guesses    
      IF (mod(et,2) == 1) THEN
        DO pt = 1,npt
          r(pt) = -1d0/3d0
          s(pt) = -1d0/3d0
        ENDDO
      ELSE IF (mod(et,2) == 0) THEN
        DO pt = 1,npt
          r(pt) = 0d0
          s(pt) = 0d0
        ENDDO
      ENDIF
      

      DO it = 1,maxit     
           
        IF (vandermonde_init) THEN
          ! Evaluate basis functions at r,s coordinates 
          CALL element_basis(et,p,n,npt,r,s,l,dldr,dlds)
          
          DO pt = 1,npt
            DO m = 1,n              
              phi(m,pt) = l(m,pt)
              phi(m,npt+pt) = dldr(m,pt)
              phi(m,2*npt+pt) = dlds(m,pt)          
            ENDDO
          ENDDO                  
        
          ! Solve linear systems to get nodal shape functions/derivatives (l,dldr,dlds) at r,s coordinates
          ! V l(r,s) = phi(r,s), V dldr(r,s) = dpdr(r,s), V dlds(r,s) = dpds(r,s)
          CALL DGETRS("N",n,3*npt,V(1,1,et),ldv,ipiv(1,et),phi,mnnds,info) 
          
          DO pt = 1,npt
            DO m = 1,n              
              l(m,pt) = phi(m,pt)
              dldr(m,pt) = phi(m,npt+pt)
              dlds(m,pt) = phi(m,2*npt+pt)          
            ENDDO
          ENDDO            
          
        ELSE
          CALL shape_functions_area_eval(et,p,n,npt,r,s,l,dldr,dlds) 
        ENDIF

        
        

        
        DO pt = 1,npt
          ! Evaluate transformation function/derivatives at r,s coordinates
          CALL element_transformation(n,elx,ely,l(:,pt),f,g, &
                                    dldr(:,pt),dlds(:,pt),drdf,drdg,dsdf,dsdg,jac)               
        
          ! Newton iteration 
          f = f - x(pt)
          g = g - y(pt)
        
          r(pt) = r(pt) - (drdf*f + drdg*g)
          s(pt) = s(pt) - (dsdf*f + dsdg*g)
          
          error(pt) = MAX(ABS(f),ABS(g))
        ENDDO
              
        err = MAXVAL(error)
        IF ( err < tol) THEN
          EXIT
        ENDIF        
       
        
      ENDDO
      
!       IF (it >= maxit) THEN
!         PRINT("(A,E22.15)"), "   MAX ITERATIONS EXCEEDED, error = ",err
!         PRINT("(2(A,F20.15))"), "   r = ",r(1), "   s = ", s(1)
!       ELSE       
!         PRINT("(A,I7,A,E22.15)"), "   iterations: ",it, "  error = ",err
!         PRINT("(2(A,F20.15))"), "   r = ",r(1), "   s = ", s(1)
!       ENDIF


      RETURN
      END SUBROUTINE xy2rs
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE element_area(ne,nel_type,np,el_type,elxy,area)
      
      USE area_qpts_mod, ONLY: area_qpts
      USE shape_functions_mod, ONLY: shape_functions_area_eval   
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: ne
      INTEGER, INTENT(IN) :: nel_type
      INTEGER, DIMENSION(:), INTENT(IN) :: np
      INTEGER, DIMENSION(:), INTENT(IN) :: el_type
      REAL(rp), DIMENSION(:,:,:), INTENT(IN) :: elxy
      REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: area
      
      INTEGER :: el,pt
      INTEGER :: et,nnd
      INTEGER :: nqpt(nel_type),nnds(nel_type)      
      INTEGER :: p,mnp,mnqpt,mnnds
      REAL(rp) :: xpt,ypt
      REAL(rp) :: drdx,drdy,dsdx,dsdy,detJ
      REAL(rp) :: x1,x2,x3,y1,y2,y3
      REAL(rp), DIMENSION(:,:), ALLOCATABLE :: wpt
      REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: qpt
      REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: psi,dpsidr,dpsids  
      
           
      
      p = 1
      mnp = maxval(np)
      CALL area_qpts(1,p,mnp,nel_type,nqpt,mnqpt,wpt,qpt)    
      
      mnnds = (mnp+1)**2     
      ALLOCATE(psi(mnnds,mnqpt,nel_type),dpsidr(mnnds,mnqpt,nel_type),dpsids(mnnds,mnqpt,nel_type))
      
      DO et = 1,nel_type         
        CALL shape_functions_area_eval(et,np(et),nnds(et),nqpt(et),qpt(:,1,et),qpt(:,2,et), &
                                       psi(:,:,et),dpsidr(:,:,et),dpsids(:,:,et))
      ENDDO
      
      
      ALLOCATE(area(ne))      
      
      DO el = 1,ne
      
        et = el_type(el)
        nnd = nnds(et)
        
        area(el) = 0d0                
        DO pt = 1,nqpt(et)
          CALL element_transformation(nnd,elxy(:,el,1),elxy(:,el,2),psi(:,pt,et),xpt,ypt, &
                                      dpsidr(:,pt,et),dpsids(:,pt,et),drdx,drdy,dsdx,dsdy,detJ)        
        
          area(el) = area(el) + wpt(pt,et)*detJ
        ENDDO
        
!         IF (et == 1) THEN
!           x1 = elxy(1,el,1)
!           x2 = elxy(2,el,1)
!           x3 = elxy(3,el,1)
!           y1 = elxy(1,el,2)
!           y2 = elxy(2,el,2) 
!           y3 = elxy(3,el,2) 
!           
!           PRINT*, area(el), .5d0*((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1))          
!         ENDIF
      ENDDO
      
      RETURN
      END SUBROUTINE element_area
      


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      END MODULE transformation