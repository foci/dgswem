     MODULE vandermonde
     
     CONTAINS
          
     
     SUBROUTINE vandermonde_area(et,p,ndf,V)
     
     USE globals, ONLY: rp
     USE basis, ONLY: element_nodes,element_basis
     
     IMPLICIT NONE
       
     INTEGER, INTENT(IN) :: et
     INTEGER, INTENT(IN) :: p
     INTEGER, INTENT(OUT) :: ndf
     REAL(rp), DIMENSION(:,:), INTENT(OUT) :: V
     
     INTEGER :: npt
     INTEGER :: info
     REAL(rp), DIMENSION((p+1)**2) :: r,s
   
   
     CALL element_nodes(et,1,p,npt,r,s)
     CALL element_basis(et,p,ndf,npt,r,s,V)
     
     
     END SUBROUTINE vandermonde_area       
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!               
      
      
      SUBROUTINE vandermonde_edge(p,n,V)
      
      USE globals, ONLY: rp
      USE basis, ONLY: lglpts,jacobi
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: p
      INTEGER, INTENT(OUT) :: n
      REAL(rp), DIMENSION(:,:), INTENT(OUT) :: V
      
      INTEGER :: pt,i
      REAL(rp) :: xi(p+1)
      REAL(rp) :: phi(p+1)            
      
    
      n = p+1
      
      CALL lglpts(p,xi)       
      
      DO i = 0,p      
      
        CALL jacobi(0,0,i,xi,n,phi)
        
        DO pt = 1,n
          V(i+1,pt) = phi(pt)
        ENDDO       
        
      ENDDO              

                  
      RETURN      
      END SUBROUTINE vandermonde_edge
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
      
      END MODULE vandermonde