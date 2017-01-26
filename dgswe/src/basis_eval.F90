      SUBROUTINE area_basis()

        USE globals, ONLY: rp,nel_type,nqpta,mnqpta,wpta,qpta,ndof,mndof, &
                           phia,dpdr,dpds                           
        USE basis, ONLY: element_basis
        USE allocation, ONLY: alloc_basis_arrays             
        USE messenger2, ONLY: myrank        
        USE read_dginp, ONLY: p
        

        IMPLICIT NONE
        INTEGER :: i,j,pt,et,dof,ndf
        REAL(rp) :: qint
        REAL(rp) :: mm(mndof,mndof)
  
        CALL alloc_basis_arrays()
      
        DO et = 1,nel_type        
          

          CALL element_basis(et,p,ndf,nqpta(et),qpta(:,1,et),qpta(:,2,et),phia(:,:,et),dpdr(:,:,et),dpds(:,:,et))         
          
!           PRINT "(A)", 'Basis functions at quadrature points'
!           DO i = 1,ndf
!             PRINT "(100(F10.3))", (phia(i,j,et),j=1,nqpta(et))
!           ENDDO
!           PRINT "(A)", ' '            
            
          
          ! Compute mass matrix (constant*indentity matrix)
          DO i = 1,ndf
            DO j = 1,ndf
              mm(i,j) = 0d0
              DO pt = 1,nqpta(et)
                mm(i,j) = mm(i,j) + wpta(pt,et)*phia(i,pt,et)*phia(j,pt,et)
              ENDDO
            ENDDO
          ENDDO
        
          
          
        ENDDO                
        
        CALL modal2nodal()
        
        IF (myrank == 0) THEN
          DO et = 1,nel_type
          
            PRINT "(A)", "---------------------------------------------"
            PRINT "(A)", "         Basis Function Information          "
            PRINT "(A)", "---------------------------------------------"
            PRINT "(A)", " "

            PRINT "(A,I5)", "Polynomial order:",p           
          
            PRINT "(A,I5)", "Number of degrees of freedom:",ndf
            PRINT "(A)", " "        

            PRINT "(A)", 'Mass matrix'
            DO i = 1,ndf
              PRINT "(100(F10.3))", (mm(i,j),j=1,ndof(et))
            ENDDO
            PRINT "(A)", ' '

!             PRINT "(A)", 'Basis functions at quadrature points'
!             DO i = 1,ndf
!               PRINT "(100(F10.3))", (phia(i,j,et),j=1,nqpta(et))
!             ENDDO
!             PRINT "(A)", ' '             
            
          ENDDO
        ENDIF


        RETURN
      END SUBROUTINE area_basis
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE edge_basis()

        USE globals, ONLY: rp,ndof,nverts,mndof,nqpte,mnqpte,qpte,wpte,phie,phie_int,nel_type
        USE basis, ONLY: element_basis
        USE read_dginp, ONLY: p

        IMPLICIT NONE

        INTEGER :: dof,pt,et,i
        INTEGER :: tpts,ndf
        REAL(rp) :: phi(mndof*4*mnqpte)       

        DO et = 1,nel_type
          
          tpts = nverts(et)*nqpte(et)
          CALL element_basis(et,p,ndf,tpts,qpte(:,1,et),qpte(:,2,et),phie(:,:,et))
          
          DO pt = 1,tpts
            DO dof = 1,ndf
              phie_int(dof,pt,et) = phie(dof,pt,et)*wpte(pt,et)
            ENDDO
          ENDDO
          
        ENDDO        


        RETURN
      END SUBROUTINE edge_basis      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!