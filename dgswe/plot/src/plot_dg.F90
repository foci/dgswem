      PROGRAM plot_dg
      
      USE plot_globals
      USE globals, ONLY: ne,nn,nverts,ndof,mndof,np,mnp,nnds,mnnds,nel_type,ect,xy,depth, &
                         nope,neta,obseg,obnds,nvel,nbou,fbseg,fbnds,bndxy,grid_name, &
                         el_type,elxy,elhb, &
                         nepn,epn,mnepn,ned,ged2el,ged2el,ged2led,ged2nn,ed_type,recv_edge, &
                         nied,iedn,nobed,obedn,nfbed,fbedn,nnfbed,nfbedn,nfbednn, &
                         psic,psiv, &
                         deg2rad
      USE grid_file_mod
      USE basis, ONLY: element_nodes,element_basis
      USE read_write_output, ONLY: read_solution_full,read_fort63,read_fort64, &
                                   read_dg63,read_dg64
      USE read_dginp, ONLY: read_input,out_direc,p,ctp,hbp,tf, &
                            grid_file,curve_file,cb_file_exists,bathy_file,hb_file_exists, &
                            sphi0,slam0, &
                            sta_opt,stations_file
      USE plot_mod, ONLY: read_colormap,setup_cbounds,plot_ref_el, &
                          scale_factors,zoom_box,make_plot,make_movie                                                                           
      USE evaluate_mod, ONLY: evaluate_basis
      USE labels_mod, ONLY: latex_axes_labels,run_latex,close_tex, &
                            latex_element_labels,latex_node_labels  
      USE axes_mod, ONLY: write_all_axes                            
                          
      USE triangulation, ONLY: reference_element_delaunay                            
      USE edge_connectivity_mod
      USE curvilinear_nodes_mod
      USE transformation
      USE shape_functions_mod
      USE version
      USE bathymetry_interp_mod, ONLY: bathymetry_nodal2modal,dgswem_bathymetry_nodal2modal
      
      IMPLICIT NONE
      
      INTEGER :: ord
      INTEGER :: start_snap,end_snap
      INTEGER :: ncall
      CHARACTER(3) :: nout
      CHARACTER(25) :: fname
      REAL(rp) :: H
      LOGICAL :: file_exists
      
      space = 0        
      
      CALL version_information(6)
      
      CALL read_plot_input()
      
  
      
      IF (mesh%plot_sol_option == 0 .and. zeta%plot_sol_option == 0 .and. &
          bathy%plot_sol_option == 0 .and. vel%plot_sol_option == 0) THEN
        
        PRINT("(A)"), "No plot options have been specified"
        STOP
      ENDIF    
      
      CALL read_input(0,input_path)
      IF (substitute_path == 1) THEN
        CALL substitute_partial_path(grid_file,replace_path,sub_path)
        CALL substitute_partial_path(curve_file,replace_path,sub_path)    
        CALL substitute_partial_path(stations_file,replace_path,sub_path)       
      ENDIF
      
      slam0 = slam0*deg2rad
      sphi0 = sphi0*deg2rad          
      
      CALL sizes()
      CALL setup_plot_types()

      PRINT("(A)"), grid_file
      CALL read_header(0,grid_file,grid_name,ne,nn)        
      CALL read_coords(nn,xy,depth)
      CALL read_connectivity(ne,ect,el_type) 
      CALL init_element_coordinates(ne,ctp,el_type,nverts,xy,ect,elxy)                  
      CALL read_open_boundaries(nope,neta,obseg,obnds)            
      CALL read_flow_boundaries(nbou,nvel,fbseg,fbnds)                       
      CALL read_curve_file(0,curve_file,ctp,nbou,xy,bndxy,cb_file_exists)      
      CALL print_grid_info(grid_file,grid_name,ne,nn)    
      
      CALL elements_per_node(ne,nn,nverts,el_type,ect,nepn,mnepn,epn)       
      CALL find_edge_pairs(ne,nverts,el_type,ect,nepn,epn,ned,ged2el,ged2nn,ged2led)      
      CALL find_interior_edges(ned,ged2el,nied,iedn,ed_type,recv_edge)      
      CALL find_open_edges(nope,obseg,obnds,ged2nn,nobed,obedn,ed_type,recv_edge)            
      CALL find_flow_edges(nbou,fbseg,fbnds,ged2nn,nnfbed,nfbedn,nfbednn,nfbed,fbedn,recv_edge,ed_type)     
      nred = 0
      CALL print_connect_info(mnepn,ned,nied,nobed,nfbed,nnfbed,nred)
      IF (sta_opt > 0) THEN
        CALL read_stations()
      ENDIF

      PRINT("(A)"), "Calculating curved boundary information..."
      CALL shape_functions_linear_at_ctp(nel_type,np,psiv)                   
      CALL eval_coordinates_curved(ctp,nnds,nverts,el_type,xy,ect,fbseg,fbnds, &
                                   nnfbed,nfbedn,nfbednn,ged2el,ged2led, &
                                   psiv,bndxy,elxy)     
      CALL element_area(ne,nel_type,np,el_type,elxy,el_area)
!       CALL find_element_init(nel_type,nverts,np,nnds,nn,xy,nepn,epn)      
      

      PRINT("(A)"), "Calculating additional ploting point coordinates..."
      nord = (p_high-p_low+1)/p_skip
      
      ALLOCATE(r(mnpp,nel_type*nord),s(mnpp,nel_type*nord))      
      ALLOCATE(psic(mnnds,mnpp,nel_type*nord))
      ALLOCATE(rect(3,3*mnpp,nel_type*nord))      
      ALLOCATE(nptri(nel_type*nord),npplt(nel_type*nord),pplt(nel_type*nord))
      ncall = 0 
      DO et = 1,nel_type 
              
        DO ord = 1,nord
          i = (et-1)*nord+ord
          pplt(i) = (ord-1)*p_skip+p_low
          CALL element_nodes(et,space,pplt(i),npplt(i),r(:,i),s(:,i))                  
          CALL shape_functions_area_eval(et,np(et),nnd,npplt(i),r(:,i),s(:,i),psic(:,:,i))  
          CALL reference_element_delaunay(et,pplt(i),npplt(i),r(:,i),s(:,i),nptri(i),rect(:,:,i))        
          
          PRINT("(4(A,I4))"), "  number of additional nodes/sub-triangles: ", npplt(i),"/",nptri(i) 
          
!           ncall = ncall + 1
!           WRITE(nout,"(I3.3)") ncall 
!           fname = "ref_el_"//nout//".ps"
!           CALL plot_ref_el(fname,figure_width,et,np(et),nptri(i),rect(:,:,i),r(:,i),s(:,i))
          
!           DO el = 1,nptri(i)
!             PRINT "(4(I5))", el,(rect(nd,el,i), nd=1,3)
!           ENDDO
!           PRINT*, "" 
        ENDDO         
        
      ENDDO                                    
           
      ALLOCATE(xyplt(mnpp,ne,2))
      DO el = 1,ne      
        et = el_type(el)                          
        nnd = nnds(et)
        i = (et-1)*nord+nord
        npts = npplt(i)
        DO pt = 1,npts              
          CALL element_transformation(nnd,elxy(:,el,1),elxy(:,el,2),psic(:,pt,i),xpt,ypt)           
          xyplt(pt,el,1) = xpt
          xyplt(pt,el,2) = ypt
        ENDDO
      ENDDO       
             
      
      PRINT("(A)"), "Evaluating reference element coordinate information..."
      CALL evaluate_basis(hbp,nord,mnpp,mndof,nel_type,npplt,r,s,ndof_hb,phi_hb)
      CALL evaluate_basis(p,nord,mnpp,mndof,nel_type,npplt,r,s,ndof_sol,phi_sol)

      
      
      
      CALL read_colormap(cmap_file)
      
      
      
      PRINT("(A)"), "Finding zoom box..."
      CALL zoom_box(ne,nord,npplt,el_type,xyplt,xbox_min,xbox_max,ybox_min,ybox_max, &
                                                     xmin,xmax,ymin,ymax,el_in)

                                                     
      PRINT("(A)"), "Scaling coordinates..."
      CALL scale_factors(figure_width,figure_height,xmin,xmax,ymin,ymax,ax,bx,ay,by)

      
    


      
#ifdef adcirc 

#elif dgswem

#else 
      snap_start = snap_start + 1
      snap_end = snap_end + 1      
#endif      
      
      nsnap_Z = snap_end
      nsnap_Qx = snap_end
      nsnap_Qy = snap_end
      nsnap_hb = 1
      
      t_start = 0d0
      t_end = tf*86400d0 
      t_snap = 0d0
      
      snap_char = "0000"      
      
          
#ifdef adcirc        
      
      IF (zeta%plot_sol_option == 1 .or. vel%plot_sol_option == 1) THEN
        PRINT("(A)"), "Reading zeta solution..."          
        CALL read_fort63(out_direc,t,eta,nsnap_Z) 
        ALLOCATE(Z(3,ne,nsnap_Z))
        DO snap = 1,nsnap_Z
          DO el = 1,ne
            DO nd = 1,3
              Z(nd,el,snap) = eta(ect(nd,el),snap)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
      IF (bathy%plot_sol_option == 1 .or. vel%plot_sol_option == 1) THEN
        PRINT("(A)"), "Reading bathymetry solution..."       
        ALLOCATE(hb(3,ne,1))
        DO el = 1,ne
          DO nd = 1,3
            hb(nd,el,1) = depth(ect(nd,el))
          ENDDO
        ENDDO
      ENDIF        
      IF (vel%plot_sol_option == 1) THEN
        PRINT("(A)"), "Reading u and v solutions..."       
        CALL read_fort64(out_direc,t,uu2,vv2,nsnap_Qx)
        ALLOCATE(Qx(3,ne,nsnap_Qx))
        ALLOCATE(Qy(3,ne,nsnap_Qx))        
        DO snap = 1,nsnap_Qx
          DO el = 1,ne
            DO nd = 1,3
              H = Z(nd,el,snap)+hb(nd,el,1)
              Qx(nd,el,snap) = uu2(ect(nd,el),snap)*H
              Qy(nd,el,snap) = vv2(ect(nd,el),snap)*H
            ENDDO
          ENDDO
        ENDDO                
      ELSE 
        ALLOCATE(Qx(1,1,1),Qy(1,1,1))
      ENDIF             
      
#elif dgswem

      IF (zeta%plot_sol_option == 1 .or. vel%plot_sol_option == 1) THEN
        PRINT("(A)"), "Reading zeta solution..."          
        CALL read_DG63(out_direc,ne,t,Z,nsnap_Z) 
      ENDIF
      IF (vel%plot_sol_option == 1) THEN
        PRINT("(A)"), "Reading Qx and Qy solutions..."       
        CALL read_DG64(out_direc,ne,t,Qx,Qy,nsnap_Qx)        
      ELSE 
        ALLOCATE(Qx(1,1,1),Qy(1,1,1))
      ENDIF  
      IF (bathy%plot_sol_option == 1 .or. vel%plot_sol_option == 1) THEN
        ALLOCATE(hb(3,ne,1))
        CALL dgswem_bathymetry_nodal2modal(ne,ect,depth,hbm)
        hb(:,:,1) = hbm(:,:)
      ENDIF  

#else
 
      IF (zeta%plot_sol_option == 1 .or. vel%plot_sol_option == 1) THEN
        PRINT("(A)"), "Reading zeta solution..."          
        CALL read_solution_full(out_direc,"Z.sol","N",t,Z,nsnap_Z) 
      ENDIF
      IF (vel%plot_sol_option == 1) THEN
        PRINT("(A)"), "Reading Qx and Qy solutions..."       
        CALL read_solution_full(out_direc,"Qx.sol","N",t,Qx,nsnap_Qx)        
        CALL read_solution_full(out_direc,"Qy.sol","N",t,Qy,nsnap_Qy) 
      ELSE 
        ALLOCATE(Qx(1,1,1),Qy(1,1,1))
      ENDIF  
      IF (bathy%plot_sol_option == 1 .or. vel%plot_sol_option == 1) THEN
        INQUIRE(file="hb.sol",exist=file_exists)
        IF (file_exists == .true.) THEN
          PRINT("(A)"), "Reading bathymetry solution..."          
          CALL read_solution_full(out_direc,"hb.sol","N",t,hb,nsnap_hb)  
        ELSE
          CALL read_bathy_file(0,bathy_file,hbp,ne,el_type,nverts,depth,ect,elhb,hb_file_exists)
          ALLOCATE(hb(ndof_hb(4),ne,1))
          CALL bathymetry_nodal2modal(hbp,ndof_hb(4),ne,el_type,elhb,hbm)
          hb(:,:,1) = hbm(:,:)
        ENDIF
      ENDIF   
 
#endif
      

      
      
      CALL setup_cbounds(ne,el_in,el_type,npplt,mesh,1,1)     
      CALL make_plot(1,0d0,mesh)  
      

      IF (zeta%plot_sol_option == 0 .and. vel%plot_sol_option == 0 .and. bathy%plot_sol_option == 0) THEN
        STOP
      ENDIF     
      
      IF (adapt_option == 1) THEN
        OPEN(UNIT = 998, FILE="error.out", STATUS="REPLACE")
        WRITE(998,"(A)") "name     snap     error_total     nptri_total     pplt_max     ne_total"
      ENDIF      
      
      CALL setup_cbounds(ne,el_in,el_type,npplt,bathy,1,1)      
      CALL make_plot(1,t_snap,bathy)                         
      
      
      IF (zeta%plot_sol_option == 0 .and. vel%plot_sol_option == 0) THEN
        STOP
      ENDIF
      
      
             
         
      IF ((snap_start > nsnap_Z) .or. (snap_start > nsnap_Qx) .or. (snap_start > nsnap_Qy)) THEN
        snap_start = MIN(nsnap_Z,nsnap_Qx,nsnap_Qy)
      ENDIF                  
      
      IF ((snap_end > nsnap_Z) .or. (snap_end > nsnap_Qx) .or. (snap_end > nsnap_Qy)) THEN
        snap_end = MIN(nsnap_Z,nsnap_Qx,nsnap_Qy)
      ENDIF
      
      
      


      CALL setup_cbounds(ne,el_in,el_type,npplt,zeta,snap_start,snap_end)
      CALL setup_cbounds(ne,el_in,el_type,npplt,vel,snap_start,snap_end)

      
      PRINT("(A)"), " "
      
      
      DO snap = snap_start,snap_end
      
        
        t_snap = t(snap)      
        PRINT("(A)"), "---------------------------------------------"
        PRINT("(A,I5,A,I5,A,F20.5)"), "Time snap: ", snap-1,"/",snap_end," t = ", t_snap        
     
      
        CALL make_plot(snap,t_snap,zeta)
              
        PRINT("(A)"), " "        

        CALL make_plot(snap,t_snap,vel)
        
      ENDDO
      
      IF (adapt_option == 1) THEN
        CLOSE(998)
      ENDIF
      
      CALL make_movie(zeta,frmt)
      CALL make_movie(vel,frmt)
  

      END PROGRAM plot_dg
      
