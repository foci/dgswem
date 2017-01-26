      SUBROUTINE read_grid()
      
      USE globals, ONLY: rp,ne,nn,nverts,nnds,ect,xy,depth,elxy,elhb,bndxy,el_size, &
                         nope,neta,obseg,obnds,nvel,nbou,fbseg,fbnds,grid_name, &
                         el_type,mnnds, &
                         r_earth,deg2rad
                             
      USE messenger2, ONLY: myrank
      USE read_dginp, ONLY: grid_file,bathy_file,curve_file, &
                            cb_file_exists,hb_file_exists, &
                            p,ctp,hbp,h0,coord_sys,slam0,sphi0
      USE spherical_mod, ONLY: cpp_transformation                                                        
      USE grid_file_mod                            

      IMPLICIT NONE

      REAL(rp) :: cfl,u      


      CALL read_header(myrank,grid_file,grid_name,ne,nn)  
      
      CALL read_coords(nn,xy,depth,h0)
      
      CALL cpp_transformation(coord_sys,r_earth,slam0,sphi0,nn,xy)

      CALL read_connectivity(ne,ect,el_type)
      
      CALL init_element_coordinates(ne,ctp,el_type,nverts,xy,ect,elxy)            
      
      CALL read_open_boundaries(nope,neta,obseg,obnds)      
      
      CALL read_flow_boundaries(nbou,nvel,fbseg,fbnds)
      
      CALL read_bathy_file(myrank,bathy_file,hbp,ne,el_type,nverts,depth,ect,elhb,hb_file_exists)
                  
      CALL read_curve_file(myrank,curve_file,ctp,nbou,xy,bndxy,cb_file_exists)      
      
      IF (myrank == 0) THEN
        CALL print_grid_info(grid_file,grid_name,ne,nn)
      ENDIF
      
      CALL grid_size(ne,el_type,ect,xy,el_size)
      
      cfl = 1d0
      u = 0d0
      !CALL courant(p,ne,u,cfl,el_type,nverts,nnds,elhb,el_size)
  
      

      RETURN
      END SUBROUTINE read_grid