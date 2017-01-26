      SUBROUTINE area_transformation()

      USE globals, ONLY: rp,ne,el_type,mndof,elxy, &
                         ndof,nnds,nqpta,wpta,dpdr,dpds,phia,mmi_init, &
                         psia,dpsidr,dpsids, &
                         detJa,dpdx_init,dpdy_init,phia,phia_int_init, &
                         r_earth
                         
      USE transformation, ONLY: element_transformation
      USE spherical_mod, ONLY: cpp_factor
      USE read_dginp, ONLY: coord_sys,slam0,sphi0
      USE lapack_interfaces

      IMPLICIT NONE
      
      INTEGER :: i,j,nd,el,pt,m,dof
      INTEGER :: et,nnd,nqa,ndf
      INTEGER :: ipiv(mndof)
      INTEGER ::info
      REAL(rp) :: work(mndof*mndof)      
      REAL(rp) :: mm(mndof,mndof)      
      REAL(rp) :: xpt,ypt
      REAL(rp) :: drdx,drdy,dsdx,dsdy
      REAL(rp) :: Sp

    
          

          
      DO el = 1,ne
        et = el_type(el)
        
        nnd = nnds(et)
        nqa = nqpta(et)
        ndf = ndof(et)
        
        mm = 0d0
        
   pts: DO pt = 1,nqa        

          CALL element_transformation(nnd,elxy(:,el,1),elxy(:,el,2),psia(:,pt,et),xpt,ypt, &
                                      dpsidr(:,pt,et),dpsids(:,pt,et),drdx,drdy,dsdx,dsdy,detJa(el,pt))
                                      
          IF (detJa(el,pt) < 0d0) THEN
            PRINT*, "Negative Jacobian determinant at area quadrature point detected, el = ",el
            PRINT*, detJa(el,pt)
          ENDIF
            
          CALL cpp_factor(coord_sys,r_earth,slam0,sphi0,ypt,Sp)
                     
            
          DO dof = 1,ndf
            i = (dof-1)*nqa + pt        
              
            dpdx_init(el,dof,pt) = wpta(pt,et)*(dpdr(dof,pt,et)*drdx + dpds(dof,pt,et)*dsdx)*detJa(el,pt)*Sp
            dpdy_init(el,dof,pt) = wpta(pt,et)*(dpdr(dof,pt,et)*drdy + dpds(dof,pt,et)*dsdy)*detJa(el,pt)  
              
            phia_int_init(el,dof,pt) = wpta(pt,et)*phia(dof,pt,et)*detJa(el,pt)
           
          ENDDO
            
            
          DO i = 1,ndf
            DO j = 1,ndf
              mm(j,i) = mm(j,i) + wpta(pt,et)*phia(i,pt,et)*phia(j,pt,et)*detJa(el,pt)
            ENDDO
          ENDDO

        ENDDO pts
          
          
        CALL DGETRF(ndf,ndf,mm,mndof,ipiv,info)       
        CALL DGETRI(ndf,mm,mndof,ipiv,work,ndf*ndf,info)

        m = 1
        DO i = 1,ndf
          DO j = 1,ndf
            mmi_init(el,m) = mm(i,j)
            m = m + 1
          ENDDO
        ENDDO          
          
      ENDDO
      
  
      RETURN
      END SUBROUTINE area_transformation
      