      MODULE bathymetry_interp_mod

      USE globals, ONLY: rp,pi
      
      IMPLICIT NONE
      
      CONTAINS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         
      
      
      SUBROUTINE shape_functions_eltypes_at_hbp(space,nel_type,np,psi,dpdr,dpds,ext,nnds)

      USE basis, ONLY: element_nodes
      USE shape_functions_mod, ONLY: shape_functions_area_eval      
      
      IMPLICIT NONE    
      
      INTEGER, INTENT(IN) :: space
      INTEGER, INTENT(IN) :: nel_type
      INTEGER, DIMENSION(:), INTENT(IN) :: np
      REAL(rp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(OUT) :: psi
      REAL(rp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: dpdr
      REAL(rp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: dpds     
      INTEGER, INTENT(IN), OPTIONAL :: ext
      INTEGER, DIMENSION(:), INTENT(INOUT), OPTIONAL :: nnds
      
      
      INTEGER :: pt,i,j
      INTEGER :: et,hbp_type     
      INTEGER :: npts,nnd
      INTEGER :: mnp,mnnds
      INTEGER :: calc_deriv,extract
      REAL(rp), DIMENSION(:), ALLOCATABLE :: r,s
      
      ! Evaluates linear/curvilinear shape functions at high-order batymetry nodal sets
      !
      ! Used to compute function-specified bathymetry at high-order batymetry nodes
      ! for curved elements. See bathy_coordinates()    
       
      mnp = maxval(np)
      mnnds = (mnp+1)**2
      ALLOCATE(r(mnnds),s(mnnds))
      
      ALLOCATE(psi(mnnds,mnnds,nel_type))             
      psi = 0d0 
      
      calc_deriv = 0
      IF (PRESENT(dpdr) .AND. PRESENT(dpds)) THEN
        calc_deriv = 1
        ALLOCATE(dpdr(mnnds,mnnds,nel_type),dpds(mnnds,mnnds,nel_type))
      ENDIF
      
      extract = 0
      IF (PRESENT(ext) .AND. PRESENT(nnds)) THEN
        extract = 1
      ENDIF
      
       
      DO et = 1,nel_type     

        IF (mod(et,2) == 1) THEN
          hbp_type = 5    
        ELSE IF (mod(et,2) == 0) THEN
          hbp_type = 6
        ENDIF   
                 
        IF (extract) THEN
          CALL element_nodes(et,space,np(hbp_type),npts,r,s,ext)
          nnds(hbp_type) = npts
        ELSE
          CALL element_nodes(et,space,np(hbp_type),npts,r,s)                  
        ENDIF
        
        IF (calc_deriv) THEN
          CALL shape_functions_area_eval(et,np(et),nnd,npts,r,s,psi(:,:,et),dpdr(:,:,et),dpds(:,:,et))                    
        ELSE
          CALL shape_functions_area_eval(et,np(et),nnd,npts,r,s,psi(:,:,et))     
        ENDIF
                
!         PRINT*, "psic"        
!         DO i = 1,nnd
!           PRINT "(20(f27.17))", (psi(i,j,et), j = 1,npts)
!         ENDDO
!         PRINT*," "


      ENDDO     
      
      END SUBROUTINE shape_functions_eltypes_at_hbp      
      
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       



      SUBROUTINE bathy_coordinates(el,nnds,nverts,el_type,elxy,psic,xyhb,depth,ect,elhb, &
                                   dpdr,dpds,dhbdx,dhbdy,nhb)
      
      USE transformation, ONLY: element_transformation          
      
      IMPLICIT NONE        
      
!       INTEGER, INTENT(IN) :: ne   
      INTEGER, INTENT(IN) :: el
      INTEGER, DIMENSION(:), INTENT(IN) :: nnds
      INTEGER, DIMENSION(:), INTENT(IN) :: nverts      
      INTEGER, DIMENSION(:), INTENT(IN) :: el_type
      REAL(rp), DIMENSION(:,:,:), INTENT(IN) :: elxy
      REAL(rp), DIMENSION(:,:,:), INTENT(INOUT) :: psic
      REAL(rp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(INOUT) :: xyhb
      REAL(rp), DIMENSION(:), INTENT(IN), OPTIONAL :: depth
      INTEGER, DIMENSION(:,:), INTENT(IN), OPTIONAL :: ect
      REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT), OPTIONAL :: elhb
      REAL(rp), DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: dpdr
      REAL(rp), DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: dpds
      REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT), OPTIONAL :: dhbdx
      REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT), OPTIONAL :: dhbdy
      REAL(rp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(INOUT), OPTIONAL :: nhb
      
      INTEGER :: interp,calc_deriv,calc_norm
      INTEGER :: pt,nd
      INTEGER :: et,npts,mnnds,nv,nnd
      INTEGER :: et_linear
      INTEGER :: nnd_coord,nnd_interp
      REAL(rp) :: xpt,ypt,hb,dhdx,dhdy      
      REAL(rp) :: drdx,drdy,dsdx,dsdy,jac,Sp
      REAL(rp) :: hbvert(4)
      REAL(rp) :: nrm

      
      
!       mnnds = MAXVAL(nnds)                 
      
      interp = 0
      IF ( PRESENT(depth) .AND. PRESENT(ect) .AND. PRESENT(elhb) ) THEN
        interp = 1
!         ALLOCATE(elhb(mnnds,ne))  
      ENDIF
      
      calc_deriv = 0
      IF (PRESENT(dpdr) .AND. PRESENT(dpds) .AND. PRESENT(dhbdx) .AND. PRESENT(dhbdy)) THEN
        calc_deriv = 1
!         ALLOCATE(dhbdx(mnnds,ne),dhbdy(mnnds,ne))
      ENDIF
      
      calc_norm = 0
      IF (PRESENT(nhb)) THEN
        calc_norm = 1
!         ALLOCATE(nhb(mnnds,ne,3))
      ENDIF

      
      
        et = el_type(el)      
      
        nv = nverts(et)
        IF (mod(et,2) == 1) THEN
          npts = nnds(5) 
          et_linear = 1
        ELSE IF (mod(et,2) == 0) THEN
          npts = nnds(6) 
          et_linear = 2
        ENDIF                   
      
        nnd_coord = nnds(et)
        nnd_interp = nnds(et_linear)    
      
        DO pt = 1,npts              

        
        
          IF (calc_deriv) THEN
            CALL element_transformation(nnd_coord,elxy(:,el,1),elxy(:,el,2),psic(:,pt,et),xpt,ypt, &
                                        dpdr(:,pt,et),dpds(:,pt,et),drdx,drdy,dsdx,dsdy,jac)
          ELSE
            CALL element_transformation(nnd_coord,elxy(:,el,1),elxy(:,el,2),psic(:,pt,et),xpt,ypt)
          ENDIF
  
          xyhb(pt,el,1) = xpt
          xyhb(pt,el,2) = ypt 
          
          
        
        
          IF (interp) THEN
            DO nd = 1,nv        
              hbvert(nd) = depth(ect(nd,el))         
            ENDDO        
                         
            Sp = 1d0             
                                              
            IF (calc_deriv) THEN
              CALL bathymetry_interp_eval(nnd_interp,hbvert,psic(:,pt,et_linear),hb, &
                                          dpdr(:,pt,et_linear),dpds(:,pt,et_linear),drdx,drdy,dsdx,dsdy,Sp,dhdx,dhdy)
                                          
              dhbdx(pt,el) = dhdx
              dhbdy(pt,el) = dhdy
            ELSE
              CALL bathymetry_interp_eval(nnd_interp,hbvert,psic(:,pt,et_linear),hb)
            ENDIF 
            
           elhb(pt,el) = hb   
!              elhb(pt,el) = 10d0
!              ypt = 500d0/(f2(xpt)-f1(xpt))*ypt - 500d0/(f2(xpt)-f1(xpt))*f1(xpt)
!              elhb(pt,el) = 10d0 - 5d0*cos(2d0*pi/500d0*ypt) 

          
            IF (calc_norm) THEN
              nrm = sqrt(dhdx**2 + dhdy**2 + 1d0)
              nhb(pt,el,1) = dhdx/nrm
              nhb(pt,el,2) = dhdy/nrm
              nhb(pt,el,3) = -1d0/nrm
            ENDIF
          ENDIF
          
          
        
        ENDDO

      
      END SUBROUTINE bathy_coordinates      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
      SUBROUTINE bathymetry_interp_eval(n,hb,l,h,dldr,dlds,drdx,drdy,dsdx,dsdy,Sp,dhdx,dhdy)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: n
      REAL(rp), DIMENSION(:), INTENT(IN) :: hb
      REAL(rp), DIMENSION(:), INTENT(IN) :: l
      REAL(rp), INTENT(OUT) :: h      
      REAL(rp), DIMENSION(:), INTENT(IN), OPTIONAL :: dldr,dlds
      REAL(rp), INTENT(IN), OPTIONAL :: drdx,drdy,dsdx,dsdy
      REAL(rp), INTENT(IN), OPTIONAL :: Sp
      REAL(rp), INTENT(OUT), OPTIONAL :: dhdx,dhdy
      
      INTEGER :: nd
      INTEGER :: calc_deriv
      
      IF( PRESENT(dldr) .AND. PRESENT(dlds) .AND. &
          PRESENT(drdx) .AND. PRESENT(drdy) .AND. &
          PRESENT(dsdx) .AND. PRESENT(dsdy) .AND. & 
          PRESENT(Sp)   .AND.                     &
          PRESENT(dhdx) .AND. PRESENT(dhdy) ) THEN
          
        calc_deriv = 1
          
      ELSE
      
        calc_deriv = 0
          
      ENDIF

      h = 0d0      
      DO nd = 1,n                
        h =  h + l(nd)*hb(nd)                           
      ENDDO         
      
      IF (calc_deriv == 1) THEN
      
        dhdx = 0d0
        dhdy = 0d0
        DO nd = 1,n  
          dhdx = dhdx + (dldr(nd)*drdx + dlds(nd)*dsdx)*hb(nd)*Sp
          dhdy = dhdy + (dldr(nd)*drdy + dlds(nd)*dsdy)*hb(nd)              
        ENDDO      
      
      ENDIF 
      
      
      RETURN
      END SUBROUTINE bathymetry_interp_eval
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      

      SUBROUTINE bathymetry_nodal2modal(hbp,mnnds,ne,el_type,elhb,hbm)
      
      USE vandermonde, ONLY: vandermonde_area
      USE lapack_interfaces
      
      IMPLICIT NONE
     
      INTEGER, INTENT(IN) :: hbp
      INTEGER, INTENT(IN) :: mnnds
      INTEGER, INTENT(IN) :: ne
      INTEGER, DIMENSION(:), INTENT(IN) :: el_type
      REAL(rp), DIMENSION(:,:), INTENT(IN) :: elhb
      REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: hbm
     
      INTEGER :: et,el,nd,dof
      INTEGER :: n,nnds(2)
      INTEGER :: info
      REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: V
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: ipiv     
     
     ALLOCATE(V(mnnds,mnnds,2))
     ALLOCATE(ipiv(mnnds,2))
     
      DO et = 1,2
        CALL vandermonde_area(et,hbp,n,V(:,:,et))
        CALL DGETRF(n,n,V(1,1,et),mnnds,ipiv(1,et),info)  
        nnds(et) = n
      ENDDO
      
      
      ALLOCATE(hbm(mnnds,ne))
      hbm = 0d0 
      DO el = 1,ne
      
        et = el_type(el)
        IF (mod(et,2) == 1) THEN
          et = 1     
        ELSE IF (mod(et,2) == 0) THEN
          et = 2
        ENDIF
        n = nnds(et)                
                
        DO nd = 1,n
          hbm(nd,el) = elhb(nd,el)
        ENDDO
        
        CALL DGETRS("T",n,1,V(1,1,et),mnnds,ipiv(1,et),hbm(1,el),mnnds,info)                         
        
      ENDDO          
     
     
      RETURN
      END SUBROUTINE bathymetry_nodal2modal      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      SUBROUTINE dgswem_bathymetry_nodal2modal(ne,ect,depth,hbm)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: ne
      INTEGER, DIMENSION(:,:), INTENT(IN) :: ect
      REAL(rp), DIMENSION(:), INTENT(IN) :: depth
      REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: hbm
      
      INTEGER :: el
      INTEGER :: n1,n2,n3      
      
      ALLOCATE(hbm(3,ne))
      
      DO el = 1,ne
         n1 = ect(1,el)
         n2 = ect(2,el)
         n3 = ect(3,el)
         
         hbm(1,el) =  1d0/3d0 * (depth(n1) + depth(n2) + depth(n3))
         hbm(2,el) = -1d0/6d0 * (depth(n1) + depth(n2)) + 1d0/3d0*depth(n3)
         hbm(3,el) = -0.5d0*depth(n1) + 0.5d0*depth(n2)
      ENDDO      
      
      RETURN
      END SUBROUTINE dgswem_bathymetry_nodal2modal      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      FUNCTION f1(x) RESULT(y)

      IMPLICIT NONE
      
      REAL(rp) :: x
      REAL(rp) :: y
      
      y = 0d0 + 100d0*(1d0/(COSH(4d0*(x-2000d0)/500d0)))
      
      END FUNCTION
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      

      FUNCTION f2(x) RESULT(y)

      IMPLICIT NONE
      
      REAL(rp) :: x
      REAL(rp) :: y
      
      y = 500d0 - 100d0*(1d0/(COSH(4d0*(x-2000d0)/500d0)))
      
      END FUNCTION      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      

      END MODULE bathymetry_interp_mod
