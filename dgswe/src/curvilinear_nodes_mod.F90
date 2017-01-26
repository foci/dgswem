      MODULE curvilinear_nodes_mod

      
      USE globals, ONLY: rp,pi
      
      
      IMPLICIT NONE
      
      INTEGER :: blend = 1      

      CONTAINS
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  


      SUBROUTINE shape_functions_linear_at_ctp(nel_type,np,psi)

      USE basis, ONLY: element_nodes
      USE shape_functions_mod, ONLY: shape_functions_area_eval      
      
      IMPLICIT NONE      
      
      INTEGER, INTENT(IN) :: nel_type
      INTEGER, DIMENSION(:), INTENT(IN) :: np
      REAL(rp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(OUT) :: psi
      
      INTEGER :: pt,i,j
      INTEGER :: et,str_typ,crv_typ
      INTEGER :: npts,nnd 
      INTEGER :: mnp,mnnds
      REAL(rp), DIMENSION(:), ALLOCATABLE :: r,s
       

      ! Evaluates linear shape functions at curved element nodal sets
      !
      ! Used to create additional nodes which can then be adjusted to make 
      ! straight elements curved. See edge_coordinates_curved()
      
      mnp = maxval(np)
      mnnds = (mnp+1)**2
      ALLOCATE(r(mnnds),s(mnnds))
      
      ALLOCATE(psi(mnnds,mnnds,2))             
      psi = 0d0 
      
      DO et = 1,2

        IF (mod(et,2) == 1) THEN
          str_typ = 1
          crv_typ = 3
        ELSE IF (mod(et,2) == 0) THEN
          str_typ = 2
          crv_typ = 4
        ENDIF  

        CALL element_nodes(et,1,np(crv_typ),npts,r,s)
        CALL shape_functions_area_eval(et,np(str_typ),nnd,npts,r,s,psi(:,:,et))   
              
!         
!         PRINT*, "psiv"
!         DO i = 1,nnd
!           PRINT "(300(f27.17))", (psi(i,j,et), j = 1,npts)
!         ENDDO
!         PRINT*," "
!         

      ENDDO
      
      END SUBROUTINE shape_functions_linear_at_ctp

      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  


      SUBROUTINE eval_coordinates_curved(ctp,nnds,nverts,el_type,xy,ect,fbseg,fbnds, &
                                         nnfbed,nfbedn,nfbednn,ged2el,ged2led, &
                                         psiv,bndxy,elxy)
     
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: ctp
      INTEGER, DIMENSION(:), INTENT(IN) :: nnds
      INTEGER, DIMENSION(:), INTENT(IN) :: nverts
      INTEGER, DIMENSION(:), INTENT(INOUT) :: el_type
      REAL(rp), DIMENSION(:,:), INTENT(IN) :: xy
      INTEGER, DIMENSION(:,:), INTENT(IN) :: ect
      INTEGER, DIMENSION(:,:), INTENT(IN) :: fbseg
      INTEGER, DIMENSION(:,:), INTENT(IN) :: fbnds
      INTEGER, INTENT(IN) :: nnfbed
      INTEGER, DIMENSION(:), INTENT(IN) :: nfbedn
      INTEGER, DIMENSION(:,:), INTENT(IN) :: nfbednn
      INTEGER, DIMENSION(:,:), INTENT(IN) :: ged2el
      INTEGER, DIMENSION(:,:), INTENT(IN) :: ged2led
      REAL(rp), DIMENSION(:,:,:), INTENT(IN) :: psiv
      REAL(rp), DIMENSION(:,:,:,:), INTENT(IN) :: bndxy
      REAL(rp), DIMENSION(:,:,:), INTENT(INOUT) :: elxy
     
      INTEGER :: i,ed,nd,j
      INTEGER :: ged,seg,n1,el,led
      INTEGER :: n1ind
      REAL(rp), DIMENSION(2,ctp-1) :: segxy     
      
      IF (ctp == 1) THEN
        RETURN
      ENDIF
      
      
      j = 0
      
      DO ed = 1,nnfbed
        ged = nfbedn(ed)
        
        seg = nfbednn(ed,1)
        n1 = nfbednn(ed,2)
        
        el = ged2el(1,ged)
        led = ged2led(1,ged)     
         
        
        n1ind = 0
 search:DO i = 1,fbseg(1,seg)-1
          IF (fbnds(i,seg) == n1) THEN
            n1ind = i
            EXIT search
          ENDIF
        ENDDO search
        
        IF (n1ind == 0) THEN
          PRINT*, "node number not found"
          STOP
        ENDIF
        
        DO nd = 1,ctp-1
          segxy(1,nd) = bndxy(1,nd,n1ind,seg)
          segxy(2,nd) = bndxy(2,nd,n1ind,seg)
        ENDDO
        
        CALL edge_coordinates_curved(el,ctp,led,nnds,nverts,el_type,xy,ect,segxy,psiv,elxy)      
        
      ENDDO     

     
     END SUBROUTINE eval_coordinates_curved


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE edge_coordinates_curved(el,ctp,led,nnds,nverts,el_type,xy,ect,segxy,psiv,elxy)
      
      USE transformation, ONLY: element_transformation            
       
      IMPLICIT NONE
       
      INTEGER, INTENT(IN) :: el
      INTEGER, INTENT(IN) :: ctp
      INTEGER, INTENT(IN) :: led      
      INTEGER, DIMENSION(:), INTENT(IN) :: nnds
      INTEGER, DIMENSION(:), INTENT(IN) :: nverts
      INTEGER, DIMENSION(:), INTENT(INOUT) :: el_type
      REAL(rp), DIMENSION(:,:), INTENT(IN) :: xy
      INTEGER, DIMENSION(:,:), INTENT(IN) :: ect
      REAL(rp), DIMENSION(:,:), INTENT(IN) :: segxy
      REAL(rp), DIMENSION(:,:,:), INTENT(IN) :: psiv
      REAL(rp), DIMENSION(:,:,:), INTENT(INOUT) :: elxy
      
      INTEGER :: pt,nd
      INTEGER :: nnd,nv,et,mnnds
      INTEGER :: pt1,pt2
      INTEGER :: reverse      
      REAL(rp) :: xpt,ypt
      REAL(rp) :: d1,d2
      REAL(rp) :: dx1,dy1,dx2,dy2
      REAL(rp), DIMENSION(:), ALLOCATABLE :: x,y
      REAL(rp), DIMENSION(:,:), ALLOCATABLE :: eddx      
       
      mnnds = MAXVAL(nnds)
      ALLOCATE(x(mnnds),y(mnnds))
      ALLOCATE(eddx(ctp+1,2))
       
       
      et = el_type(el)
      nv = nverts(et)  
        
      IF (et <= 2) THEN            ! only compute extra nodes if element isn't already curved
                                   ! otherwise, just adjust the boundary nodes. 
        IF (mod(et,2) == 1) THEN   ! this is important for elements with two no normal flow boundaries
          el_type(el) = 3 
          nnd = nnds(3)      
        ELSE IF (mod(et,2) == 0) THEN
          el_type(el) = 4
          nnd = nnds(4)                  
        ENDIF         
        
        DO nd = 1,nv
          x(nd) = xy(1,ect(nd,el))
          y(nd) = xy(2,ect(nd,el))
        ENDDO        
        
        DO pt = 1,nnd               

          CALL element_transformation(nv,x,y,psiv(:,pt,et),xpt,ypt)
        
          elxy(pt,el,1) = xpt
          elxy(pt,el,2) = ypt
        ENDDO      
        
      ENDIF        
      
      d1 = 0d0                                   ! Make sure the same boundary orientation is used  
      d2 = 0d0                                   ! between the segxy data and elxy. 
      DO nd = 1,ctp-1                            ! The orientation that gives the minimum total difference 
        pt1 = mod(led,nv)*ctp + 1 + nd           ! should be the correct one
        pt2 = mod(led,nv)*ctp + 1 + ctp - nd

        dx1 = segxy(1,nd)-elxy(pt1,el,1)
        dy1 = segxy(2,nd)-elxy(pt1,el,2)
        
        dx2 = segxy(1,nd)-elxy(pt2,el,1)
        dy2 = segxy(2,nd)-elxy(pt2,el,2)     
        
        d1 = d1 + sqrt(dx1**2+dy1**2)
        d2 = d2 + sqrt(dx2**2+dy2**2)

      ENDDO      
      
      reverse = 0
      IF (d2 < d1) THEN
        reverse = 1
      ENDIF
      
      eddx = 0d0                
      
      DO nd = 1,ctp-1                   ! adjust mid-edge nodes to spline coordinates, 
                                        ! verticies have already been adjusted in read_curve_file
        IF (reverse == 0) THEN                                        
          pt = mod(led,nv)*ctp + 1 + nd                                        
        ELSE IF (reverse == 1) THEN
          pt = mod(led,nv)*ctp + 1 + ctp - nd
        ENDIF
        
!         ytest = elxy(pt,el,2)
!         xpt = elxy(pt,el,1) 
!           
!         IF (ytest < 250d0) THEN
!           ypt = 0d0 + 100d0*(1d0/(COSH(4d0*(xpt-2000d0)/500d0)))
!         ELSE IF (ytest > 250d0) THEN
!           ypt = 500d0 - 100d0*(1d0/(COSH(4d0*(xpt-2000d0)/500d0)))
!         ENDIF
!         
!         eddx(nd+1,2) = ypt - elxy(pt,el,2)        

        eddx(nd+1,1) = segxy(1,nd)-elxy(pt,el,1)
        eddx(nd+1,2) = segxy(2,nd)-elxy(pt,el,2)

      ENDDO           
      
      IF (blend) THEN
        IF (mod(et,2) == 1) THEN
          CALL tri_blend_coordinates(et,el,ctp,led,mnnds,eddx,elxy)
        ELSE IF (mod(et,2) == 0) THEN
          CALL quad_blend_coordinates(et,el,ctp,led,mnnds,eddx,elxy)
!           CALL set_coordinates(nv,el,ctp,led,eddx,elxy)
        ENDIF
      ELSE
        CALL set_coordinates(nv,el,ctp,led,eddx,elxy)
      ENDIF

      RETURN
      END SUBROUTINE edge_coordinates_curved


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE set_coordinates(nv,el,ctp,led,eddx,elxy)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: nv
      INTEGER, INTENT(IN) :: el
      INTEGER, INTENT(IN) :: ctp
      INTEGER, INTENT(IN) :: led
      REAL(rp), DIMENSION(:,:), INTENT(IN) :: eddx
      REAL(rp), DIMENSION(:,:,:), INTENT(INOUT) :: elxy
      
      INTEGER :: nd,pt
      
        DO nd = 1,ctp-1                   
          pt = mod(led,nv)*ctp + 1 + nd                      

          elxy(pt,el,1) = elxy(pt,el,1) + eddx(nd+1,1)
          elxy(pt,el,2) = elxy(pt,el,2) + eddx(nd+1,2) 
        ENDDO        
      
      RETURN
      END SUBROUTINE set_coordinates

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
      
      SUBROUTINE tri_blend_coordinates(et,el,ctp,led,mnnds,eddx,elxy)
      
      USE basis, ONLY: element_nodes,jacobi
      USE transformation, ONLY: xy2rs
      USE lapack_interfaces

      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: et
      INTEGER, INTENT(IN) :: el
      INTEGER, INTENT(IN) :: ctp
      INTEGER, INTENT(IN) :: led
      INTEGER, INTENT(IN) :: mnnds
      REAL(rp), DIMENSION(:,:), INTENT(INOUT) :: eddx
      REAL(rp), DIMENSION(:,:,:), INTENT(INOUT) :: elxy
      
      INTEGER :: i,j,pt
      INTEGER :: nnd,nrhs
      INTEGER :: npt_el,npt_ed
      INTEGER :: ipiv(mnnds),info
      REAL(rp) :: r(mnnds),s(mnnds)
      REAL(rp) :: er(mnnds),es(mnnds)
      REAL(rp) :: elr(mnnds),edr(mnnds)
      REAL(rp) :: Ved(ctp+1,ctp+1),Vel(mnnds,ctp+1)
      REAL(rp) :: phi(mnnds)
      REAL(rp) :: dx(mnnds),dy(mnnds)
      REAL(rp) :: blend,denom
!       REAL(rp) :: elxy_old(mnnds,2)
      
      

        CALL element_nodes(et,1,ctp,npt_el,r,s)
        CALL element_nodes(et,1,ctp,npt_ed,er,es,led)      

        
        SELECT CASE(led)
          CASE (3)
            elr = r  
            edr = er
          CASE (1)
            elr = s
            edr = es
          CASE (2)
            elr = s
            edr = es
        END SELECT    
        


 
        nnd = ctp + 1       
        
        DO i = 1,nnd
          CALL jacobi(0,0,i-1,edr,npt_ed,phi)    
          DO pt = 1,npt_ed              
            Ved(pt,i) = phi(pt)
          ENDDO
        ENDDO
        
        

        CALL DGETRF(nnd,nnd,Ved,nnd,ipiv,info)  

        nrhs = 2
        CALL DGETRS('N',nnd,nrhs,Ved,nnd,ipiv,eddx,nnd,info)
   
        
        DO i = 1,nnd
          CALL jacobi(0,0,i-1,elr,npt_el,phi)    
          DO pt = 1,npt_el              
            Vel(pt,i) = phi(pt)
          ENDDO
        ENDDO
        

        
        dx = 0d0
        dy = 0d0
        DO j = 1,nnd        
          DO i = 1,npt_el
            dx(i) = dx(i) + Vel(i,j)*eddx(j,1)
            dy(i) = dy(i) + Vel(i,j)*eddx(j,2)            
          ENDDO
        ENDDO
     

!         elxy_old = elxy(:,el,:)
 
        DO pt = 1,npt_el
          
          blend = 0d0
          denom = 1d0-elr(pt)
          IF (abs(denom) > 1d-8) THEN
          
            SELECT CASE (led)
              CASE(3)
                blend = -(r(pt)+s(pt))/denom
              CASE(1)
                blend = (r(pt)+1d0)/denom
              CASE(2)
                blend = -(r(pt)+s(pt))/denom
            END SELECT
            
          ENDIF
          
          elxy(pt,el,1) = elxy(pt,el,1) + blend*dx(pt)
          elxy(pt,el,2) = elxy(pt,el,2) + blend*dy(pt)
        ENDDO
        
        
!         CALL xy2rs(et,ctp,elxy_old(:,1),elxy_old(:,2),npt_el,elxy(:,el,1),elxy(:,el,2),r,s)
        
!         DO pt = 1,npt_el
!           PRINT*, r(pt),s(pt)
!         ENDDO
 
      END SUBROUTINE tri_blend_coordinates
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        

      SUBROUTINE quad_blend_coordinates(et,el,ctp,led,mnnds,eddx,elxy)
      
      USE basis, ONLY: element_nodes,jacobi
      USE lapack_interfaces

      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: et
      INTEGER, INTENT(IN) :: el
      INTEGER, INTENT(IN) :: ctp
      INTEGER, INTENT(IN) :: led
      INTEGER, INTENT(IN) :: mnnds
      REAL(rp), DIMENSION(:,:), INTENT(INOUT) :: eddx
      REAL(rp), DIMENSION(:,:,:), INTENT(INOUT) :: elxy
      
      INTEGER :: i,j,pt
      INTEGER :: nnd,nrhs
      INTEGER :: npt_el,npt_ed
      INTEGER :: ipiv(mnnds),info
      REAL(rp) :: r(mnnds),s(mnnds)
      REAL(rp) :: er(mnnds),es(mnnds)
      REAL(rp) :: elr(mnnds),edr(mnnds)
      REAL(rp) :: Ved(ctp+1,ctp+1),Vel(mnnds,ctp+1)
      REAL(rp) :: phi(mnnds)
      REAL(rp) :: dx(mnnds),dy(mnnds)
      REAL(rp) :: blend,denom
      
      

        CALL element_nodes(et,1,ctp,npt_el,r,s)
        CALL element_nodes(et,1,ctp,npt_ed,er,es,led)      

        
        SELECT CASE(led)
          CASE (4)
            elr = r  
            edr = er
          CASE (1)
            elr = s
            edr = es
          CASE (2)
            elr = r
            edr = er
          CASE (3)
            elr = s
            edr = es
        END SELECT    
        


 
        nnd = ctp + 1       
        
        DO i = 1,nnd
          CALL jacobi(0,0,i-1,edr,npt_ed,phi)    
          DO pt = 1,npt_ed              
            Ved(pt,i) = phi(pt)
          ENDDO
        ENDDO
        
        

        CALL DGETRF(nnd,nnd,Ved,nnd,ipiv,info)  

        nrhs = 2
        CALL DGETRS('N',nnd,nrhs,Ved,nnd,ipiv,eddx,nnd,info)
   
        
        DO i = 1,nnd
          CALL jacobi(0,0,i-1,elr,npt_el,phi)    
          DO pt = 1,npt_el              
            Vel(pt,i) = phi(pt)
          ENDDO
        ENDDO
        

        
        dx = 0d0
        dy = 0d0
        DO j = 1,nnd        
          DO i = 1,npt_el
            dx(i) = dx(i) + Vel(i,j)*eddx(j,1)
            dy(i) = dy(i) + Vel(i,j)*eddx(j,2)            
          ENDDO
        ENDDO
     

 
        DO pt = 1,npt_el          
         
            SELECT CASE (led)            
              CASE(4)
                blend = .5d0*(1d0-s(pt))
              CASE(1)
                blend = .5d0*(1d0+r(pt))
              CASE(2)
                blend = .5d0*(1d0+s(pt))
              CASE(3)
                blend = .5d0*(1d0-r(pt))
            END SELECT
          
          elxy(pt,el,1) = elxy(pt,el,1) + blend*dx(pt)
          elxy(pt,el,2) = elxy(pt,el,2) + blend*dy(pt)
        ENDDO
        
        
        
 
 
      END SUBROUTINE quad_blend_coordinates
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
      
      END MODULE curvilinear_nodes_mod