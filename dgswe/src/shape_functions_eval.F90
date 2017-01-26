      SUBROUTINE shape_functions_area_qpts()

      USE globals, ONLY: rp,nel_type,order,nverts, &
                         mnqpta,nqpta,mnqpte,nqpte,np, &
                         qpta,qpte,psia,dpsidr,dpsids
      USE shape_functions_mod, ONLY: shape_functions_area_eval
      
      IMPLICIT NONE      
      
      INTEGER :: pt,i,j,dof
      INTEGER :: typ,et,eo,tpts
      INTEGER :: nv,nnd,nqa,nqe,p
      INTEGER :: info      
      REAL(rp) :: r(mnqpta+4*mnqpte),s(mnqpta+4*mnqpte)

      psia = 0d0
      dpsidr = 0d0
      dpsids = 0d0
      
      
      DO typ = 1,2*nel_type
      
        IF (typ <= 4) THEN
          et = typ
        ELSE
          et = typ - 4
        ENDIF
        
        eo = order(typ)        
        p = np(eo)           
        
        nv = nverts(et)
        nqa = nqpta(et)
        nqe = nv*nqpte(et)   
      
        tpts = nqa+nqe     
      
        DO pt = 1,nqa
          r(pt) = qpta(pt,1,et)
          s(pt) = qpta(pt,2,et)
        ENDDO  
      
        DO i = 1,nqe
          pt = nqa+i
          r(pt) = qpte(i,1,et)
          s(pt) = qpte(i,2,et)
        ENDDO


        CALL shape_functions_area_eval(nv,p,nnd,tpts,r,s,psia(:,:,typ),dpsidr(:,:,typ),dpsids(:,:,typ))                 
!       
!       PRINT*, "psi : "      
!       DO i = 1,nnd
!         PRINT("(300(F27.17))"), (psia(i,j,et), j = 1,tpts)
!       ENDDO         

!       PRINT*, "dpsidr : "      
!       DO i = 1,nnd
!         PRINT("(300(F27.17))"), (dpsidr(i,j,et), j = 1,tpts)
!       ENDDO     
! 
!       PRINT*, "dpsids : "      
!       DO i = 1,nnd
!         PRINT("(300(F27.17))"), (dpsids(i,j,et), j = 1,tpts)
!       ENDDO            

       ENDDO
       
      END SUBROUTINE shape_functions_area_qpts
           
     
     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                 
       
       
      SUBROUTINE shape_functions_edge_qpts()

      USE globals, ONLY: rp,nel_type, &
                         mnqpte,nqpte,np,mnp, &
                         qpte,psie,dpsidxi
      USE shape_functions_mod, ONLY: shape_functions_edge_eval                         
      
      IMPLICIT NONE      
      
      INTEGER :: pt,i
      INTEGER :: et
      INTEGER :: nqe,p,n
      INTEGER :: info      
      REAL(rp) :: xi(mnqpte)
      REAL(rp) :: phi(mnqpte),dphi(mnqpte)
      
      ! Evaulates 1-D edge shape functions at 1-D guass points
      !
      ! Used in computing the Jacobians for the edge integrals.
      ! See edge_transformation.F90
       
      DO et = 1,nel_type
       
        p = np(et)    
        n = p+1
        nqe = nqpte(et)        
       
        DO pt = 1,nqe
          xi(pt) = qpte(pt,2,et) ! These are the 1-D guass points
        ENDDO
      
        CALL shape_functions_edge_eval(p,n,nqe,xi,psie(:,:,et),dpsidxi(:,:,et))
 
      ENDDO

      RETURN
      END SUBROUTINE shape_functions_edge_qpts