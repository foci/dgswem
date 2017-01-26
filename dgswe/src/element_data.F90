      SUBROUTINE element_data()

      USE globals, ONLY: rp,ect,xy,depth,elxy,ne,nn,ned,ndof,mnnds, &
                         nqpta,mnqpta,qpta,nqpte, &
                         el_type,nel_type,np,nnds,nverts, &
                         area,edlen,edlen_area,normal,ged2nn,ged2el,ged2led, &
                         dhbdx_init,dhbdy_init, &
                         detJa,detJe, &
                         nx_pt,ny_pt, &
                         elhb,hbm
                         
      USE basis, ONLY:       
      USE allocation, ONLY: alloc_trans_arrays
      USE read_dginp, ONLY: p,ctp,hbp
      USE messenger2, ONLY: myrank        
      USE bathymetry_interp_mod, ONLY: bathymetry_nodal2modal

      IMPLICIT NONE
      INTEGER :: el,ed,led,dof,pt,i,nd
      INTEGER :: ind,et
      REAL(rp) :: x1,x2,x3,x4,y1,y2,y3,y4
      REAL(rp) :: dxdr,dxds,dydr,dyds   
      REAL(rp) :: drdx,drdy,dsdx,dsdy      
      REAL(rp) :: hb1,hb2,hb3,hb4,hb,dhbdx1,dhbdy1 
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: r,s
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: jac    
    
      CALL alloc_trans_arrays()                       
      
      IF (ctp > 1) THEN  
        IF (myrank == 0) PRINT "(A)", "Creating curvilinear element information..."
        CALL curvilinear()    
      ENDIF           
      
      IF (myrank == 0) PRINT "(A)", "Calculating area shape functions..."
      CALL shape_functions_area_qpts()                         

      IF (myrank == 0) PRINT "(A)", "Computing area transformations..."
      CALL area_transformation() 
      
      IF (myrank == 0) PRINT "(A)", "Computing element normals..."
      CALL normals()
      
      IF (myrank == 0) PRINT "(A)", "Calculating edge shape functions..."
      CALL shape_functions_edge_qpts()      

      IF (myrank == 0) PRINT "(A)", "Computing edge transformations..."
      CALL edge_transformation()
      
      IF (myrank == 0) PRINT "(A)", "Interpolating bathymety onto area quadrature points..."
      CALL bathymetry_interp_area_qpts()
      
      IF (myrank == 0) PRINT "(A)", "Interpolating bathymety onto edge quadrature points..."      
      CALL bathymetry_interp_edge_qpts()
      
      IF (myrank == 0) PRINT "(A)", "Computing modal bathymetry..."           
      CALL bathymetry_nodal2modal(hbp,mnnds,ne,el_type,elhb,hbm)
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Calculate element area 
      !!!!!!!!!!!!!!!!!!!!!!!!!!
      
      ALLOCATE(r(mnqpta),s(mnqpta),jac(ne,mnqpta))
      DO i = 1,nqpta(2)
        r(i) = qpta(i,1,2)
        s(i) = qpta(i,2,2)
      ENDDO

!       PRINT*, " "
      DO el = 1,ne
        x1 = xy(1,ect(1,el))
        y1 = xy(2,ect(1,el))

        x2 = xy(1,ect(2,el))
        y2 = xy(2,ect(2,el))

        x3 = xy(1,ect(3,el))
        y3 = xy(2,ect(3,el))

        area(el) = .5d0*((x2*y3-x3*y2) + (x3*y1-x1*y3) + (x1*y2-x2*y1))
        IF (el_type(el) == 1) THEN
!           PRINT("(I25,15(e25.16))"), el,area(el),(2d0*detJa(el,pt), pt = 1,nqpta(1))
!           PRINT("(I25,15(e25.16))"), el,(abs(area(el)-2d0*detJa(el,pt)), pt = 1,nqpta(1))
        ELSE IF (el_type(el) == 2) THEN
          x4 = xy(1,ect(4,el))
          y4 = xy(2,ect(4,el))
          
          DO pt = 1,nqpta(2)
          
            dxdr = .25d0*((-1d0+s(pt))*x1 + ( 1d0-s(pt))*x2 + (1d0+s(pt))*x3 + (-1d0-s(pt))*x4)
            dxds = .25d0*((-1d0+r(pt))*x1 + (-1d0-r(pt))*x2 + (1d0+r(pt))*x3 + ( 1d0-r(pt))*x4)
            dydr = .25d0*((-1d0+s(pt))*y1 + ( 1d0-s(pt))*y2 + (1d0+s(pt))*y3 + (-1d0-s(pt))*y4)
            dyds = .25d0*((-1d0+r(pt))*y1 + (-1d0-r(pt))*y2 + (1d0+r(pt))*y3 + ( 1d0-r(pt))*y4)
            
            jac(el,pt) = dxdr*dyds - dxds*dydr
          ENDDO
         
!           PRINT("(I5,4(e23.14),10x,4(e23.14))"), el,(jac(el,pt), pt = 1,nqpta(2)),(detJa(el,pt), pt = 1,nqpta(2))         
!           PRINT("(I5,4(e23.14))"), el,(abs(jac(el,pt)-detJa(el,pt)), pt = 1,nqpta(2))    
!           PRINT("(I5,10x,40(e23.14))"), el,(detJa(el,pt), pt = 1,nqpta(2))     
          
        ELSE IF ( el_type(el) == 4) THEN
        
!           PRINT("(I5,10x,40(e23.14))"), el,(detJa(el,pt), pt = 1,nqpta(2))  

        ENDIF
      ENDDO
     


      !!!!!!!!!!!!!!!!!!!!!!!!!
      ! Calculate edge length 
      !!!!!!!!!!!!!!!!!!!!!!!!!      

      DO ed = 1,ned
        x1 = xy(1,ged2nn(1,ed))
        y1 = xy(2,ged2nn(1,ed))

        x2 = xy(1,ged2nn(2,ed))
        y2 = xy(2,ged2nn(2,ed))

        edlen(ed) = sqrt((x2-x1)**2 + (y2-y1)**2)
        
!         PRINT("(I5,e23.14,10x,4(e23.14))"), ed,0.5d0*edlen(ed),(detJe(ed,pt),pt = 1,nqpte(1))
!         PRINT("(I5,4(e23.14))"), ed,(abs(0.5d0*edlen(ed)-detJe(ed,pt)),pt = 1,nqpte(1))

        edlen_area(1,ed) = edlen(ed)/area(ged2el(1,ed))
        IF (ged2el(2,ed) /= 0) THEN
          edlen_area(2,ed) = edlen(ed)/area(ged2el(2,ed))   
        ENDIF    

! Set up for alternative imlementation where all edlen_area
! multiplications are done right before/during the edge integration instead of after numerical flux calculation
! The idea was to vectorize these multiplications.  Didn't show any improvement.
!         el = ged2el(1,ed)     
!         led = ged2led(1,ed)
!         edlen_area(el,led) = edlen(ed)/area(el)
! 
!         el = ged2el(2,ed)
!         led = ged2led(2,ed)
!         IF (el /= 0) THEN
!           edlen_area(el,led) = edlen(ed)/area(el)   
!         ENDIF 
      ENDDO
        

      !!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Calculate edge normals 
      !!!!!!!!!!!!!!!!!!!!!!!!!!

      DO ed = 1,ned
        x1 = xy(1,ged2nn(1,ed))
        y1 = xy(2,ged2nn(1,ed))

        x2 = xy(1,ged2nn(2,ed))
        y2 = xy(2,ged2nn(2,ed))

        normal(1,ed) =  (y2-y1)/edlen(ed)
        normal(2,ed) = -(x2-x1)/edlen(ed)
        
!         PRINT("(I5,2(e23.14),10x,4(e23.14))"), ed,normal(1,ed),normal(2,ed),(nx_pt(ed,i),i=1,nqpte(1)), (ny_pt(ed,i),i=1,nqpte(1))
!         PRINT("(3(I5),4(e23.14))"), ed,ged2el(1,ed),ged2led(1,ed),(abs(normal(2,ed)-ny_pt(ed,i)),i=1,nqpte(1))
!         PRINT("(3(I5),4(e23.14))"), ed,ged2el(1,ed),ged2led(1,ed),(abs(normal(1,ed)-nx_pt(ed,i)),i=1,nqpte(1))
      ENDDO
            

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Calculate bathymetry derivatives 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!       DO el = 1,nn
!         depth(el) = 10d0
!       ENDDO
      

 

!       DO el = 1,ne
!         x1 = xy(1,ect(1,el))
!         y1 = xy(2,ect(1,el))
! 
!         x2 = xy(1,ect(2,el))
!         y2 = xy(2,ect(2,el))
! 
!         x3 = xy(1,ect(3,el))
!         y3 = xy(2,ect(3,el))
!         
!         hb1 = depth(ect(1,el))
!         hb2 = depth(ect(2,el))
!         hb3 = depth(ect(3,el))
! 
!         dhbdx_init(el) = ( -(.5d0*(y3-y1)+.5d0*(y1-y2))*hb1 + .5d0*(y3-y1)*hb2 + .5d0*(y1-y2)*hb3 )/area(el)
!         dhbdy_init(el) = ( -(.5d0*(x1-x3)+.5d0*(x2-x1))*hb1 + .5d0*(x1-x3)*hb2 + .5d0*(x2-x1)*hb3 )/area(el)
!         
! !         dhbdx(el) = ( -(.5d0*(y3-y1)+.5d0*(y1-y2))*hb1 + .5d0*(y3-y1)*hb2 + .5d0*(y1-y2)*hb3 )/area(el)
! !         dhbdy(el) = ( -(.5d0*(x1-x3)+.5d0*(x2-x1))*hb1 + .5d0*(x1-x3)*hb2 + .5d0*(x2-x1)*hb3 )/area(el)
!  
! !         WRITE(64,*) dhbdx(el),dhbdy(el)
! 
!       ENDDO

      DO el = 1,ne
        et = el_type(el)
        

        
          x1 = xy(1,ect(1,el))
          y1 = xy(2,ect(1,el))

          x2 = xy(1,ect(2,el))
          y2 = xy(2,ect(2,el))

          x3 = xy(1,ect(3,el))
          y3 = xy(2,ect(3,el))
        
          hb1 = depth(ect(1,el))
          hb2 = depth(ect(2,el))
          hb3 = depth(ect(3,el))

        IF (et == 1) THEN          

          dhbdx1 = ( -(.5d0*(y3-y1)+.5d0*(y1-y2))*hb1 + .5d0*(y3-y1)*hb2 + .5d0*(y1-y2)*hb3 )/area(el)
          dhbdy1 = ( -(.5d0*(x1-x3)+.5d0*(x2-x1))*hb1 + .5d0*(x1-x3)*hb2 + .5d0*(x2-x1)*hb3 )/area(el)
          
!           PRINT("(10x,I5,3(e23.14))"), el,(abs(dhbdx1-dhbdx_init(el,i)),i=1,nqpta(1))
!           PRINT("(10x,I5,3(e23.14))"), el,(abs(dhbdy1-dhbdy_init(el,i)),i=1,nqpta(1))          
          
        ELSE IF (et == 2) THEN
        
          x4 = xy(1,ect(4,el))
          y4 = xy(2,ect(4,el))
          
          hb4 = depth(ect(4,el))
          
          DO pt = 1,nqpta(2)
          
            dxdr = .25d0*((-1d0+s(pt))*x1 + ( 1d0-s(pt))*x2 + (1d0+s(pt))*x3 + (-1d0-s(pt))*x4)
            dxds = .25d0*((-1d0+r(pt))*x1 + (-1d0-r(pt))*x2 + (1d0+r(pt))*x3 + ( 1d0-r(pt))*x4)
            dydr = .25d0*((-1d0+s(pt))*y1 + ( 1d0-s(pt))*y2 + (1d0+s(pt))*y3 + (-1d0-s(pt))*y4)
            dyds = .25d0*((-1d0+r(pt))*y1 + (-1d0-r(pt))*y2 + (1d0+r(pt))*y3 + ( 1d0-r(pt))*y4)
            
            jac(el,pt) = dxdr*dyds - dxds*dydr            
            
            drdx =  dyds/jac(el,pt)
            drdy = -dxds/jac(el,pt)
            dsdx = -dydr/jac(el,pt)
            dsdy =  dxdr/jac(el,pt)
            
            dhbdx1 = .25d0*(((-1d0+s(pt))*drdx + (-1d0+r(pt))*dsdx)*hb1 &
                         + (( 1d0-s(pt))*drdx + (-1d0-r(pt))*dsdx)*hb2 &
                         + (( 1d0+s(pt))*drdx + ( 1d0+r(pt))*dsdx)*hb3 &
                         + ((-1d0-s(pt))*drdx + ( 1d0-r(pt))*dsdx)*hb4)
            dhbdy1 = .25d0*(((-1d0+s(pt))*drdy + (-1d0+r(pt))*dsdy)*hb1 &
                         + (( 1d0-s(pt))*drdy + (-1d0-r(pt))*dsdy)*hb2 &
                         + (( 1d0+s(pt))*drdy + ( 1d0+r(pt))*dsdy)*hb3 &
                         + ((-1d0-s(pt))*drdy + ( 1d0-r(pt))*dsdy)*hb4)
                     
!           PRINT("(I5,4(e23.14))"), el,(abs(dhbdx1-dhbdx_init(el,i)),i=1,nqpta(2))
!           PRINT("(I5,4(e23.14))"), el,(abs(dhbdy1-dhbdy_init(el,i)),i=1,nqpta(2))                       
          ENDDO
        
                
        ENDIF
        
!           PRINT("(I5,e23.14,10x,3(e23.14))"), el,dhbdx1,(dhbdx_init(el,i),i=1,nqpta(1))
!           PRINT("(I5,e23.14,10x,3(e23.14))"), el,dhbdy1,(dhbdy_init(el,i),i=1,nqpta(1))               
        
      ENDDO
      

      RETURN
      END SUBROUTINE element_data