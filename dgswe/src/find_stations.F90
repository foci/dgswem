      SUBROUTINE find_stations()

      USE globals, ONLY:rp,ndof,mndof,np,mnnds,nnds,nverts,nel_type,nn,nepn,epn, &                                    
                        xy,el_type,elxy,elhb,elsta,hbsta,nsta,xysta,phi_sta
      USE read_dginp, ONLY:p,hbp                        
      USE basis, ONLY: element_basis      
      USE shape_functions_mod, ONLY: shape_functions_area_eval      
      USE find_element, ONLY: find_element_init,in_element
      USE bathymetry_interp_mod, ONLY: bathymetry_interp_eval

      IMPLICIT NONE

      INTEGER :: sta
      INTEGER :: elin,et,n,npts,nv
      REAL(rp) :: r(1),s(1),hb,rs(2)
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: l,phi           
      

      CALL find_element_init(nel_type,nverts,np,nnds,nn,xy,nepn,epn)         

      ALLOCATE(phi_sta(mndof,nsta))      
      ALLOCATE(elsta(nsta),hbsta(nsta,1))
      ALLOCATE(l(mnnds,1))
      ALLOCATE(phi(mndof,1))         
      
      
      ! Find element each station is located in and determine local coordinates of station
      npts = 1
      DO sta = 1,nsta
      
        CALL in_element(xysta(:,sta),el_type,elxy,elin,rs)       

        et = el_type(elin)                  
        nv = nverts(et)     
        
        r(1) = rs(1)
        s(1) = rs(2)
        
        CALL shape_functions_area_eval(nv,hbp,n,npts,r,s,l)
        CALL bathymetry_interp_eval(n,elhb(:,elin),l(:,1),hb)            
        CALL element_basis(et,p,ndof(et),npts,r,s,phi)           
        
        phi_sta(:,sta) = phi(:,1)           
        elsta(sta) = elin    
        hbsta(sta,1) = hb           
      
      ENDDO

      END SUBROUTINE find_stations