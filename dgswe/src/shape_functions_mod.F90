      MODULE shape_functions_mod
      
      CONTAINS
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
      
      SUBROUTINE shape_functions_area_eval(nv,p,nnd,npts,r,s,psi,dpdr,dpds)
      
      USE globals, ONLY: rp
      USE basis, ONLY: element_basis
      USE lapack_interfaces
      USE vandermonde, ONLY: vandermonde_area
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: nv,p,npts
      INTEGER, INTENT(OUT) :: nnd
      REAL(rp), DIMENSION(:), INTENT(IN) :: r,s
      REAL(rp), DIMENSION(:,:), INTENT(OUT) :: psi
      REAL(rp), DIMENSION(:,:), OPTIONAL, INTENT(OUT) :: dpdr,dpds      
            
      INTEGER :: info
      INTEGER :: ldv,ldp
      INTEGER :: calc_derv
      INTEGER, DIMENSION(:), ALLOCATABLE :: ipiv
      REAL(rp), DIMENSION(:,:), ALLOCATABLE :: V
      
      IF (PRESENT(dpdr) .AND. PRESENT(dpds)) THEN
        calc_derv = 1
      ELSE
        calc_derv = 0
      ENDIF
            
      ldv = (p+1)**2

      ALLOCATE(V(ldv,ldv))
      ALLOCATE(ipiv(ldv))

      ldp = SIZE(psi,1)

      CALL vandermonde_area(nv,p,nnd,V)
      
      IF (calc_derv == 0) THEN
      
        CALL element_basis(nv,p,nnd,npts,r,s,psi) 
        CALL DGESV(nnd,npts,V,ldv,ipiv,psi,ldp,info)      
        
      ELSE 
      
        CALL element_basis(nv,p,nnd,npts,r,s,psi,dpdr,dpds) 
        CALL DGESV(nnd,npts,V,ldv,ipiv,psi,ldp,info)    
        CALL DGETRS("N",nnd,npts,V,ldv,ipiv,dpdr,ldp,info) 
        CALL DGETRS("N",nnd,npts,V,ldv,ipiv,dpds,ldp,info)         
      
      ENDIF
      
      END SUBROUTINE shape_functions_area_eval
      

      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    
      SUBROUTINE shape_functions_edge_eval(p,n,npts,xi,psi,dpdxi)

      USE globals, ONLY: rp
      USE basis, ONLY: jacobi,djacobi
      USE vandermonde, ONLY: vandermonde_edge
      
      IMPLICIT NONE    
      
      INTEGER, INTENT(IN) :: p
      INTEGER, INTENT(OUT) :: n
      INTEGER, INTENT(IN) :: npts
      REAL(rp), DIMENSION(:), INTENT(IN) :: xi
      REAL(rp), DIMENSION(:,:), INTENT(OUT) :: psi
      REAL(rp), DIMENSION(:,:), OPTIONAL, INTENT(OUT) :: dpdxi
      
      
      
      INTEGER :: pt,i
      INTEGER :: ldp
      INTEGER :: info      
      REAL(rp) :: phi(npts),dphi(npts)  
      INTEGER :: ipiv(p+1)
      REAL(rp) :: V(p+1,p+1)
      INTEGER :: calc_derv
      
      IF (PRESENT(dpdxi)) THEN
        calc_derv = 1
      ELSE 
        calc_derv = 0
      ENDIF
       
      ldp = SIZE(psi,1)
      
      CALL vandermonde_edge(p,n,V)
      
      DO i = 0,p
      
        CALL jacobi(0,0,i,xi,npts,phi)
        
        DO pt = 1,npts          
          psi(i+1,pt) = phi(pt)
        ENDDO
      
      ENDDO

      CALL DGESV(n,npts,V,n,ipiv,psi,ldp,info)  
      
      
      
      IF (calc_derv == 1) THEN
        DO i = 0,p
      
          CALL djacobi(0,0,i,xi,npts,dphi)
        
          DO pt = 1,npts          
            dpdxi(i+1,pt) = dphi(pt)
          ENDDO
      
        ENDDO      
        
        CALL DGETRS("N",n,npts,V,n,ipiv,dpdxi,ldp,info)
      ENDIF 

      RETURN
      END SUBROUTINE shape_functions_edge_eval      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      

      
      END MODULE shape_functions_mod