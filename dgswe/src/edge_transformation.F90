      SUBROUTINE edge_transformation()
      
      USE globals, ONLY: rp,ned,el_type,nverts,np,nqpte,qpte, &
                         psie,dpsidxi,detJe, &
                         ect,xy,elxy,ged2el,ged2led
      USE basis, ONLY: lglpts,jacobi,djacobi
      
      IMPLICIT NONE
      
      INTEGER :: pt,nd,ed
      INTEGER :: n,nnds,et,nqpts  
      INTEGER :: el_in,el_ex,led_in
      INTEGER :: nv,ind 
      REAL(rp) :: x,y,dxdxi,dydxi
      
      
      
      DO ed = 1,ned
      
        el_in = ged2el(1,ed)        
        led_in = ged2led(1,ed)
        
        et = el_type(el_in)
        nv = nverts(et)
        nqpts = nqpte(et)     
        n = np(et)
          
          DO pt = 1,nqpts
          
            dxdxi = 0d0
            dydxi = 0d0
            
            DO nd = 1,n+1
              IF (led_in == nv-1 .and. nd == n+1) THEN
                ind = 1
              ELSE 
                ind = mod(led_in,nv)*n+nd
              ENDIF
              
              x = elxy(ind,el_in,1)
              y = elxy(ind,el_in,2)
              
              dxdxi = dxdxi + dpsidxi(nd,pt,et)*x
              dydxi = dydxi + dpsidxi(nd,pt,et)*y
              
            ENDDO
            
            detJe(ed,pt) = sqrt(dxdxi*dxdxi + dydxi*dydxi)
          
          ENDDO
           
      ENDDO
          
      
      
      RETURN
      END SUBROUTINE edge_transformation