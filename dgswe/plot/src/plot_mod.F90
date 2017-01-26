      MODULE plot_mod
      
      USE globals, ONLY: rp,nel_type,np,nnds,psic
      USE plot_globals, ONLY: plot_type,cscale_width,lr_margin,top_margin,dash,fontsize, &
                              axes_width,axes_height, &
                              rmin_page,rmax_page,smin_page,smax_page, &
                              rmin_axes,rmax_axes,smin_axes,smax_axes, &
                              rmin_cbar,rmax_cbar,smin_cbar,smax_cbar, &
                              rmin_tbar,rmax_tbar,smin_tbar,smax_tbar, &
                              xticklabel_pad,yticklabel_pad,cticklabel_pad, &
                              xlabel_pad,ylabel_pad,clabel_pad, &
                              nxtick,nytick,nctick, &
                              nxdec,nydec,ncdec,ntdec, &
                              dr_xlabel,ds_ylabel,ds_clabel,main_font, &
                              ncolors,colors, &
                              mnpp,ps,pc,nord,adapt_option, &
                              phi_hb,phi_sol,ndof_hb,ndof_sol, &
                              ax,ay,bx,by, &
                              xbox_min,xbox_max,ybox_min,ybox_max
      
      IMPLICIT NONE
      
      INTEGER :: unit_count = 99   
       
      
      CONTAINS
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE make_plot(snap,t_snap,fig)
            
      USE globals, ONLY: ne,nn,nverts,el_type,xy,ect,elxy, &
                         nbou,fbseg,fbnds,bndxy, &
                         nnfbed,nfbedn, &
                         nfbed,fbedn, &
                         nobed,obedn, &
                         ged2el,ged2led, &
                         nsta,xysta
      USE read_dginp, ONLY: p,ctp
      USE plot_globals, ONLY: el_in,t_start,t_end,xyplt,pplt,npplt,nptri,rect,r,s, &
                              frmt,density,pc,el_area
      USE labels_mod, ONLY: latex_axes_labels,run_latex,read_latex, & 
                            latex_element_labels,latex_node_labels, &
                            write_latex_ps_body,remove_latex_files, &
                            write_char_array
      USE axes_mod, ONLY: write_all_axes
            
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: snap
      REAL(rp), INTENT(IN) :: t_snap
      TYPE(plot_type), INTENT(INOUT) :: fig      

      INTEGER :: el,et
      INTEGER :: pe
      INTEGER, DIMENSION(:), ALLOCATABLE :: el_in_all
      CHARACTER(:), ALLOCATABLE :: filename
      CHARACTER(4) :: snap_char
      
 
      IF (fig%plot_sol_option < 1) THEN
        RETURN
      ENDIF     
                
      IF (.not. ALLOCATED(fig%el_plt)) THEN
        ALLOCATE(fig%el_plt(ne))
      ENDIF                   
              
      DO el = 1,ne
        et = el_type(el)
        fig%el_plt(el) = (et-1)*nord + nord
      ENDDO              

      fig%sol_min = fig%cscale_vals(snap,2)
      fig%sol_max = fig%cscale_vals(snap,3)                                      
      
                
      PRINT("(2(A,F10.5))"), "  min value = ", fig%sol_min, "  max value = ", fig%sol_max
       
      IF (fig%type_flag == 3 .or. fig%type_flag == 4) THEN 
        WRITE(snap_char,"(I4.4)") snap           
        filename = TRIM(ADJUSTL(fig%name))//"_"//snap_char
        WRITE(fig%cscale_unit,"(3(e24.17,1x))") t_snap,fig%sol_min,fig%sol_max        
      ELSE
        filename = TRIM(ADJUSTL(fig%name))
      ENDIF
      
      
          
      CALL latex_axes_labels(fig,t_snap,t_start,t_end)  
      IF (fig%el_label_option /= "off") THEN
        CALL latex_element_labels(ne,fig%el_label_option,el_type,el_in,nverts,xy,ect,nnfbed,nfbedn,ged2el)
      ENDIF
      IF (fig%nd_label_option /= "off") THEN
        CALL latex_node_labels(nn,fig%nd_label_option,xy,nbou,fbseg,fbnds)
      ENDIF
      CALL run_latex()
      CALL read_latex(fig%latex_header,fig%nline_header,fig%latex_body,fig%nline_body)
      CALL remove_latex_files()      
      
           
      CALL write_psheader(filename//".ps",fig%ps_unit)  
      CALL write_char_array(fig%ps_unit,fig%nline_header,fig%latex_header)  
      CALL plot_background(fig%ps_unit,.75d0,.75d0,.75d0)
      
      
      IF (fig%cbar_flag == 1) THEN
        CALL plot_filled_contours_adapt(fig%ps_unit,ne,el_type,el_in,el_area,elxy,xyplt,pplt,nptri,npplt,rect,r,s,snap,fig)             
      ENDIF
      IF (fig%plot_lines_option == 1) THEN
        CALL plot_line_contours(fig%ps_unit,ne,el_type,el_in,nptri,rect,xyplt,r,s,snap,fig)          
      ENDIF      
      
!         CALL plot_cb_nodes(fig%ps_unit,ctp,nbou,fbseg,fbnds,xy,bndxy)
!         CALL plot_elxy_nodes(fig%ps_unit,ne,el_type,nnds,elxy)

      IF (adapt_option == 1 .and. fig%type_flag > 1) THEN
        CALL SYSTEM("cp "//filename//".ps "//filename//"_pltmesh.ps")
        OPEN(UNIT=999,FILE=filename//"_pltmesh.ps",POSITION="APPEND")
        CALL plot_vis_mesh(999,ne,el_type,el_in,xyplt,nptri,rect,fig)   
        CALL write_all_axes(999,fig%cbar_flag,fig%tbar_flag,t_snap,t_start,t_end)
        CALL write_char_array(999,fig%nline_body,fig%latex_body)        
        CALL convert_ps(filename//"_pltmesh",frmt,density,fig%rm_ps)         
      ENDIF 
      
      

      IF (fig%plot_mesh_option == 1) THEN
        IF (fig%name == "mesh") THEN
          CALL fill_elements(fig%ps_unit,ne,nverts,fig%el_plt,pplt,el_type,el_in,xy,ect,xyplt)   
        ENDIF
        CALL plot_mesh(fig%ps_unit,ne,nverts,fig%el_plt,pplt,el_type,el_in,xy,ect,xyplt)         
        IF (fig%name == "mesh" .and. fig%plot_sta_option == 1) THEN
          CALL plot_stations(fig%ps_unit,fig%sta_start,fig%sta_end,xysta)        
        ENDIF
      ELSE IF (fig%plot_mesh_option == 2) THEN
        pe = (et-1)*nord + nord
        CALL plot_boundaries(fig%ps_unit,nverts,fig%el_plt,pplt,nobed,obedn,ged2el,ged2led,el_type,el_in,ect,xy,xyplt)       
        CALL plot_boundaries(fig%ps_unit,nverts,fig%el_plt,pplt,nnfbed,nfbedn,ged2el,ged2led,el_type,el_in,ect,xy,xyplt)        
        CALL plot_boundaries(fig%ps_unit,nverts,fig%el_plt,pplt,nfbed,fbedn,ged2el,ged2led,el_type,el_in,ect,xy,xyplt)          
      ENDIF 
              
      
      CALL write_all_axes(fig%ps_unit,fig%cbar_flag,fig%tbar_flag,t_snap,t_start,t_end)       
      
!       CALL write_latex_ps_body(fig%ps_unit)
      CALL write_char_array(fig%ps_unit,fig%nline_body,fig%latex_body)  
      
      CALL close_ps(filename,fig%ps_unit)
      CALL convert_ps(filename,frmt,density,fig%rm_ps)      
      
      RETURN
      END SUBROUTINE make_plot


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      SUBROUTINE prep_station_plot(fig,axis_label_option,snap,xmin,xmax,ymin,ymax,filename)
      
      USE plot_globals, ONLY: figure_width,figure_height,frmt,density, &
                              ax,bx,ay,by
      USE labels_mod, ONLY: write_texheader,run_latex,read_latex, & 
                            write_xyaxis_labels,remove_latex_files, &
                            write_char_array
      USE axes_mod, ONLY: write_xyaxis                             
      
      IMPLICIT NONE
            
      TYPE(plot_type), INTENT(INOUT) :: fig   
      CHARACTER(*), INTENT(IN) :: axis_label_option
      INTEGER, INTENT(IN) :: snap
      REAL(rp), INTENT(IN) :: xmax,xmin,ymin,ymax
      
      CHARACTER(:), ALLOCATABLE, INTENT(IN) :: filename     

           
      CALL scale_factors(figure_width,figure_height,xmin,xmax,ymin,ymax,ax,bx,ay,by)      
      
      CALL write_texheader()  
      CALL write_xyaxis_labels(axis_label_option)      
      
      CALL run_latex()
      CALL read_latex(fig%latex_header,fig%nline_header,fig%latex_body,fig%nline_body)
      CALL remove_latex_files()     
                
      CALL write_psheader(filename//".ps",fig%ps_unit)  
      CALL write_char_array(fig%ps_unit,fig%nline_header,fig%latex_header)  
      CALL plot_background(fig%ps_unit,1d0,1d0,1d0)      
               
      
      RETURN
      END SUBROUTINE prep_station_plot
      


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  


      SUBROUTINE finish_station_plot(fig,filename)
      
      USE plot_globals, ONLY: frmt,density
      USE labels_mod, ONLY: write_char_array
      USE axes_mod, ONLY: write_xyaxis           
      
      IMPLICIT NONE
      
      TYPE(plot_type), INTENT(INOUT) :: fig      
      CHARACTER(*), INTENT(IN) :: filename
      
      CALL write_xyaxis(fig%ps_unit)       
      CALL write_char_array(fig%ps_unit,fig%nline_body,fig%latex_body)  
      
      CALL close_ps(filename,fig%ps_unit)
      CALL convert_ps(filename,frmt,density,fig%rm_ps)            
                  
      RETURN
      END SUBROUTINE finish_station_plot
      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      SUBROUTINE zoom_box(ne,nord,npplt,el_type,xyplt,xbox_min,xbox_max,ybox_min,ybox_max,xmin,xmax,ymin,ymax,el_in)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: ne
      INTEGER, INTENT(IN) :: nord
      INTEGER, DIMENSION(:), INTENT(IN) :: npplt      
      INTEGER, DIMENSION(:), INTENT(IN) :: el_type
      REAL(rp), DIMENSION(:,:,:), INTENT(IN) :: xyplt
      REAL(rp), INTENT(IN) :: xbox_min
      REAL(rp), INTENT(IN) :: xbox_max
      REAL(rp), INTENT(IN) :: ybox_min
      REAL(rp), INTENT(IN) :: ybox_max
      REAL(rp), INTENT(OUT) :: xmin
      REAL(rp), INTENT(OUT) :: xmax
      REAL(rp), INTENT(OUT) :: ymin
      REAL(rp), INTENT(OUT) :: ymax
      INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: el_in
      
      INTEGER :: el,pt
      INTEGER :: et,npts
      INTEGER :: mnpts
      INTEGER, DIMENSION(:), ALLOCATABLE :: outside
      INTEGER :: in_flag
      REAL(rp) :: xpt,ypt     
      
      
      xmax = -1d10
      ymax = -1d10
      xmin = 1d10
      ymin = 1d10      
      
      ALLOCATE(outside(maxval(npplt)))
      
      ALLOCATE(el_in(ne))
      el_in = 1
      
      DO el = 1,ne      
        et = el_type(el)                          
        npts = npplt((et-1)*nord+nord)
        outside = 0
        DO pt = 1,npts  
        
          xpt = xyplt(pt,el,1)
          ypt = xyplt(pt,el,2)
        
          IF (xpt > xmax) THEN
            xmax = xpt
          ENDIF
          
          IF (xpt < xmin) THEN
            xmin = xpt
          ENDIF

          IF (ypt > ymax) THEN
            ymax = ypt
          ENDIF
          
          IF (ypt < ymin) THEN
            ymin = ypt
          ENDIF          
          
          IF (xpt < xbox_min .or. xpt > xbox_max) THEN
            outside(pt) = 1
            xmin = xbox_min
            xmax = xbox_max
          ENDIF
          
          IF (ypt < ybox_min .or. ypt > ybox_max) THEN
            outside(pt) = 1
            ymin = ybox_min
            ymax = ybox_max
          ENDIF

        ENDDO
        
        in_flag = 0
        DO pt = 1,npts
          IF (outside(pt) == 0) THEN
            in_flag = 1
          ENDIF
        ENDDO
        
        IF (in_flag == 0) THEN
          el_in(el) = 0
        ENDIF
      ENDDO    
      
!       CALL in_element(xy,el_type,elxy,el_found,rs)
            
      
      RETURN
      END SUBROUTINE zoom_box
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      SUBROUTINE scale_factors(figure_width,figure_height,xmin,xmax,ymin,ymax,ax,bx,ay,by)
      
      IMPLICIT NONE
      
      REAL(rp), INTENT(IN) :: figure_width
      REAL(rp), INTENT(IN) :: figure_height
      REAL(rp), INTENT(IN) :: xmin
      REAL(rp), INTENT(IN) :: xmax
      REAL(rp), INTENT(IN) :: ymin
      REAL(rp), INTENT(IN) :: ymax
      REAL(rp), INTENT(OUT) :: ax
      REAL(rp), INTENT(OUT) :: bx
      REAL(rp), INTENT(OUT) :: ay
      REAL(rp), INTENT(OUT) :: by
!       REAL(rp), DIMENSION(:,:,:), INTENT(INOUT) :: xyplt
!       REAL(rp), DIMENSION(:,:), INTENT(INOUT) :: xy
!       REAL(rp), DIMENSION(:,:,:), INTENT(INOUT) :: elxy      
      
      INTEGER :: el,nd
      INTEGER :: et,npts,nv,nnd
               
      lr_margin = (rmax_page - figure_width)/2d0 
      axes_width = figure_width/1.37d0
      
      cscale_width = axes_width*.05d0
      
      dash = 0.2d0*cscale_width
      xticklabel_pad = 2d0*dash
      yticklabel_pad = 2d0*dash
      xlabel_pad = 10d0*dash
      ylabel_pad = 15d0*dash        
      cticklabel_pad = 2d0*dash
      clabel_pad = 15d0*dash 
      

         
      
      
      rmin_axes = rmin_page + lr_margin + ylabel_pad
      rmax_axes = rmax_page - lr_margin - 2d0*cscale_width - clabel_pad
      smin_axes = smin_page + lr_margin
      smax_axes = smax_page      
      
    
            
      ax = ((rmax_axes-rmin_axes)/(xmax-xmin))
      bx = (rmin_axes*xmax-rmax_axes*xmin)/(xmax-xmin)
      
      IF (figure_height < 0d0) THEN
        ay = ax     ! axis equal            
      ELSE
        axes_height = (rmax_axes-rmin_axes)/1.25d0
        smax_axes = smin_axes+axes_height
        ay = ((smax_axes-smin_axes)/(ymax-ymin))      
      ENDIF
      by = smin_axes-ay*ymin              
      
      smax_axes = ay*ymax + by
      
      IF (smax_axes > smax_page) THEN
        PRINT*, "Error: Equal axis scaling puts y axis off page"
        STOP
      ENDIF
      
      rmin_tbar = rmax_axes + cscale_width
      rmax_tbar = rmin_tbar + clabel_pad
      smin_tbar = smax_axes - 2d0*cscale_width
      smax_tbar = smin_tbar + cscale_width   
      
      RETURN
      END SUBROUTINE scale_factors

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
      
      
      SUBROUTINE write_psheader(file_name,file_unit)
      
      USE labels_mod, ONLY: write_latex_ps_header
      
      IMPLICIT NONE
      
      CHARACTER(*), INTENT(IN) :: file_name
      INTEGER, INTENT(OUT) :: file_unit
      
      CHARACTER(1000) :: line
      INTEGER :: read_stat       

      
      
      unit_count = unit_count + 1
      file_unit = unit_count
      
      OPEN(UNIT=file_unit,FILE=file_name)      
      WRITE(file_unit,"(A)") "%!PS-Adobe-3.0"
!       WRITE(file_unit,"(A)") "%!PS-Adobe-3.0 EPSF-3.0"      
!       WRITE(file_unit,"(A,4(F9.5,1x))") "%%BoundingBox: ",rmin_axes,smin_axes,rmax_axes,smax_axes                
      
      WRITE(file_unit,"(A)") "/draw-line {"
      WRITE(file_unit,"(A)") "newpath"        
      WRITE(file_unit,"(A)") "moveto"
      WRITE(file_unit,"(A)") "lineto"
      WRITE(file_unit,"(A)") ".5 setlinewidth 2 setlinejoin"
      WRITE(file_unit,"(A)") "stroke"      
      WRITE(file_unit,"(A)") "} def"       
      
      WRITE(file_unit,"(A)") "/draw-tri-element {"
      WRITE(file_unit,"(A)") "newpath"        
      WRITE(file_unit,"(A)") "moveto"
      WRITE(file_unit,"(A)") "lineto"
      WRITE(file_unit,"(A)") "lineto"
      WRITE(file_unit,"(A)") "closepath"
      WRITE(file_unit,"(A)") ".5 setlinewidth 2 setlinejoin"
      WRITE(file_unit,"(A)") "stroke"      
      WRITE(file_unit,"(A)") "} def" 
      
      WRITE(file_unit,"(A)") "/draw-quad-element {"
      WRITE(file_unit,"(A)") "newpath"        
      WRITE(file_unit,"(A)") "moveto"
      WRITE(file_unit,"(A)") "lineto"
      WRITE(file_unit,"(A)") "lineto"
      WRITE(file_unit,"(A)") "lineto"      
      WRITE(file_unit,"(A)") "closepath"
      WRITE(file_unit,"(A)") ".5 setlinewidth 2 setlinejoin"
      WRITE(file_unit,"(A)") "stroke"      
      WRITE(file_unit,"(A)") "} def"    
      
      WRITE(file_unit,"(A)") "/fill-tri-element {"
      WRITE(file_unit,"(A)") "newpath"        
      WRITE(file_unit,"(A)") "moveto"
      WRITE(file_unit,"(A)") "lineto"
      WRITE(file_unit,"(A)") "lineto"
      WRITE(file_unit,"(A)") "closepath"
!       WRITE(file_unit,"(A)") "gsave 0 0 1 setrgbcolor fill grestore"     
      WRITE(file_unit,"(A)") "gsave 1 1 1 setrgbcolor fill grestore"         
      WRITE(file_unit,"(A)") "} def" 
      
      WRITE(file_unit,"(A)") "/fill-quad-element {"
      WRITE(file_unit,"(A)") "newpath"        
      WRITE(file_unit,"(A)") "moveto"
      WRITE(file_unit,"(A)") "lineto"
      WRITE(file_unit,"(A)") "lineto"
      WRITE(file_unit,"(A)") "lineto"      
      WRITE(file_unit,"(A)") "closepath"
!       WRITE(file_unit,"(A)") "gsave 0 0 1 setrgbcolor fill grestore"  
      WRITE(file_unit,"(A)") "gsave 1 1 1 setrgbcolor fill grestore"        
      WRITE(file_unit,"(A)") "} def"       
      
      WRITE(file_unit,"(A)") "/draw-box {"
      WRITE(file_unit,"(A)") "newpath"        
      WRITE(file_unit,"(A)") "moveto"
      WRITE(file_unit,"(A)") "lineto"
      WRITE(file_unit,"(A)") "lineto"
      WRITE(file_unit,"(A)") "lineto"      
      WRITE(file_unit,"(A)") "closepath"
      WRITE(file_unit,"(A)") ".5 setlinewidth 2 setlinejoin"
      WRITE(file_unit,"(A)") "stroke"      
      WRITE(file_unit,"(A)") "} def"   
      
      WRITE(file_unit,"(A)") "/draw-line {"
      WRITE(file_unit,"(A)") "newpath"      
      WRITE(file_unit,"(A)") "moveto"
      WRITE(file_unit,"(A)") "lineto"
      WRITE(file_unit,"(A)") ".5 setlinewidth 2 setlinejoin"
      WRITE(file_unit,"(A)") "stroke"      
      WRITE(file_unit,"(A)") "} def"      
      
      WRITE(file_unit,"(A)") "/draw-dot {"
      WRITE(file_unit,"(A)") "newpath"      
      WRITE(file_unit,"(A)") " 2 0 360 arc closepath"
!       WRITE(file_unit,"(A)") "stroke"
!       WRITE(file_unit,"(A)") "gsave 0 0 0 setrgbcolor fill grestore"
      WRITE(file_unit,"(A)") "gsave setrgbcolor fill grestore"
      WRITE(file_unit,"(A)") "} def"
      
      WRITE(file_unit,"(A)") "/trifill {"
      WRITE(file_unit,"(A)") "/A exch def"
      WRITE(file_unit,"(A)") "/B exch def"
      WRITE(file_unit,"(A)") "/C exch def"     
      WRITE(file_unit,"(A)") "/ds ["
      WRITE(file_unit,"(A)") "  0 A aload pop" 
      WRITE(file_unit,"(A)") "  0 B aload pop"
      WRITE(file_unit,"(A)") "  0 C aload pop"
      WRITE(file_unit,"(A)") "] def"
      WRITE(file_unit,"(A)") "newpath"
      WRITE(file_unit,"(A)") "<<"
      WRITE(file_unit,"(A)") "  /ShadingType 4"
      WRITE(file_unit,"(A)") "  /ColorSpace [ /DeviceRGB ]"
      WRITE(file_unit,"(A)") "  /DataSource ds" 
!       WRITE(file_unit,"(A)")   " /AntiAlias true"
      WRITE(file_unit,"(A)") ">>"
      WRITE(file_unit,"(A)") "closepath"      
      WRITE(file_unit,"(A)") "shfill"
      WRITE(file_unit,"(A)") "} def"  
      
      WRITE(file_unit,"(A)") "/recfill {"
      WRITE(file_unit,"(A)") "/A exch def"
      WRITE(file_unit,"(A)") "/B exch def"
      WRITE(file_unit,"(A)") "/C exch def"
      WRITE(file_unit,"(A)") "/D exch def"      
      WRITE(file_unit,"(A)") "<<"
      WRITE(file_unit,"(A)") " /ShadingType 2" 
      WRITE(file_unit,"(A)") " /ColorSpace /DeviceRGB"
      WRITE(file_unit,"(A)") " /BBox A "      
      WRITE(file_unit,"(A)") " /Coords B "
      WRITE(file_unit,"(A)") " /Function"
      WRITE(file_unit,"(A)") " <<"
      WRITE(file_unit,"(A)") "  /FunctionType 2"
      WRITE(file_unit,"(A)") "  /Domain [0 1]"
      WRITE(file_unit,"(A)") "  /C0 C "
      WRITE(file_unit,"(A)") "  /C1 D "
      WRITE(file_unit,"(A)") "  /N 1"
      WRITE(file_unit,"(A)") " >>"
      WRITE(file_unit,"(A)") ">>" 
      WRITE(file_unit,"(A)") "shfill"
      WRITE(file_unit,"(A)") "} def"  
      
      WRITE(file_unit,"(A)") "/scale-box{"    
      WRITE(file_unit,"(A)") "/ymax exch def" 
      WRITE(file_unit,"(A)") "/ymin exch def"   
      WRITE(file_unit,"(A)") "/xmax exch def"    
      WRITE(file_unit,"(A)") "/xmin exch def" 
      WRITE(file_unit,"(A,F20.5,A)") "/ax ",ax, " def"     
      WRITE(file_unit,"(A,F20.5,A)") "/bx ",bx, " def" 
      WRITE(file_unit,"(A,F20.5,A)") "/ay ",ay, " def"  
      WRITE(file_unit,"(A,F20.5,A)") "/by ",by, " def"    
      WRITE(file_unit,"(A)") "xmin ax mul bx add" 
      WRITE(file_unit,"(A)") "ymin ay mul by add" 
      WRITE(file_unit,"(A)") "xmax ax mul bx add"    
      WRITE(file_unit,"(A)") "ymin ay mul by add" 
      WRITE(file_unit,"(A)") "xmax ax mul bx add" 
      WRITE(file_unit,"(A)") "ymax ay mul by add"    
      WRITE(file_unit,"(A)") "xmin ax mul bx add" 
      WRITE(file_unit,"(A)") "ymax ay mul by add"       
      WRITE(file_unit,"(A)") "newpath" 
      WRITE(file_unit,"(A)") "moveto" 
      WRITE(file_unit,"(A)") "lineto" 
      WRITE(file_unit,"(A)") "lineto" 
      WRITE(file_unit,"(A)") "lineto"
      WRITE(file_unit,"(A)") "closepath"   
      WRITE(file_unit,"(A)") "gsave 3 setlinewidth 2 setlinejoin 0 0 0 setrgbcolor stroke grestore"         
      
      WRITE(file_unit,"(A)") "} def"   
        
      
      ! center vertically, left align horizontally
      WRITE(file_unit,"(A)") "/caxis-tick-labels {"
      WRITE(file_unit,"(A)") "moveto gsave dup false charpath pathbbox grestore exch pop sub 2 div 0 exch rmoveto pop show"
      WRITE(file_unit,"(A)") "} def"
      

      ! top align vertically, center horizontally
      WRITE(file_unit,"(A)") "/xaxis-labels {"
      WRITE(file_unit,"(A)") "moveto gsave dup false charpath pathbbox grestore"
      WRITE(file_unit,"(A)") "3 1 roll"      
      WRITE(file_unit,"(A)") "4 1 roll"    
      WRITE(file_unit,"(A)") "sub neg" 
      WRITE(file_unit,"(A)") "3 1 roll" 
      WRITE(file_unit,"(A)") "sub 2 div neg" 
      WRITE(file_unit,"(A)") "exch" 
      WRITE(file_unit,"(A)") "rmoveto" 
      WRITE(file_unit,"(A)") "show" 
      WRITE(file_unit,"(A)") "} def"         
      
      ! center vertically, right align horizontally
      WRITE(file_unit,"(A)") "/yaxis-tick-labels {"
      WRITE(file_unit,"(A)") "moveto gsave dup false charpath pathbbox grestore"
      WRITE(file_unit,"(A)") "3 1 roll"      
      WRITE(file_unit,"(A)") "4 1 roll"    
      WRITE(file_unit,"(A)") "sub 2 div neg" 
      WRITE(file_unit,"(A)") "3 1 roll" 
      WRITE(file_unit,"(A)") "sub neg" 
      WRITE(file_unit,"(A)") "exch" 
      WRITE(file_unit,"(A)") "rmoveto" 
      WRITE(file_unit,"(A)") "show" 
      WRITE(file_unit,"(A)") "} def"  
      
      ! rotate text, center vertically, right align horizontally
      WRITE(file_unit,"(A)") "/yaxis-labels {"   
      WRITE(file_unit,"(A)") "moveto gsave dup false charpath pathbbox grestore"
      WRITE(file_unit,"(A)") "3 1 roll"      
      WRITE(file_unit,"(A)") "4 1 roll"    
      WRITE(file_unit,"(A)") "sub neg" 
      WRITE(file_unit,"(A)") "3 1 roll" 
      WRITE(file_unit,"(A)") "sub 2 div neg" 
      WRITE(file_unit,"(A)") "rmoveto" 
      WRITE(file_unit,"(A)") "90 rotate"        
      WRITE(file_unit,"(A)") "show" 
      WRITE(file_unit,"(A)") "-90 rotate"      
      WRITE(file_unit,"(A)") "} def"  
      
      ! rotate text, center vertically, left align horizontally      
      WRITE(file_unit,"(A)") "/caxis-labels {"   
      WRITE(file_unit,"(A)") "moveto gsave dup false charpath pathbbox grestore"
      WRITE(file_unit,"(A)") "3 1 roll"      
      WRITE(file_unit,"(A)") "4 1 roll"    
      WRITE(file_unit,"(A)") "sub " 
      WRITE(file_unit,"(A)") "3 1 roll" 
      WRITE(file_unit,"(A)") "sub 2 div neg" 
      WRITE(file_unit,"(A)") "rmoveto"  
      WRITE(file_unit,"(A)") "90 rotate"        
      WRITE(file_unit,"(A)") "show" 
      WRITE(file_unit,"(A)") "-90 rotate"      
      WRITE(file_unit,"(A)") "} def"      
      
      WRITE(file_unit,"(A)") "/choosefont {"
      WRITE(file_unit,"(A)") "findfont"      
      WRITE(file_unit,"(I3,A)") fontsize," scalefont"  
      WRITE(file_unit,"(A)") "setfont"
      WRITE(file_unit,"(A)") "} def"   
      
      WRITE(file_unit,"(A)") TRIM(ADJUSTL(main_font))//" choosefont"      
!       WRITE(file_unit,"(A)") "(/usr/share/fonts/type1/gsfonts/n021003l.pfb)     choosefont"    
!       WRITE(file_unit,"(A)") "(/usr/share/fonts/type1/gsfonts/cmr10.pfb)     choosefont"        
      
            
!       CALL write_latex_ps_header(file_unit) 

     
 
  
      RETURN
      END SUBROUTINE write_psheader
      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      

      SUBROUTINE plot_background(file_unit,r,b,g)
     
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: file_unit
      REAL(rp), INTENT(IN) :: r,b,g
     
      WRITE(file_unit,"(A)") "gsave newpath"        
      WRITE(file_unit,"(2(F9.5,1x),A)") rmin_axes,smin_axes,"moveto"
      WRITE(file_unit,"(2(F9.5,1x),A)") rmax_axes,smin_axes,"lineto"
      WRITE(file_unit,"(2(F9.5,1x),A)") rmax_axes,smax_axes,"lineto"      
      WRITE(file_unit,"(2(F9.5,1x),A)") rmin_axes,smax_axes,"lineto"      
      WRITE(file_unit,"(A)") "closepath"     
      WRITE(file_unit,"(3(F9.5,1x),A)") r,g,b," setrgbcolor fill grestore"          
     
      RETURN
      END SUBROUTINE plot_background


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      
      SUBROUTINE close_ps(filename,unit_number)
      
      IMPLICIT NONE
      
      CHARACTER(*), INTENT(IN) :: filename
      INTEGER, INTENT(IN) :: unit_number                    
          
!       WRITE(unit_number,"(A)") "showpage"
      CLOSE(unit_number)            
!       CALL SYSTEM("ps2eps -B -C < "//filename//".ps > "//filename//".eps")
      
      RETURN
      END SUBROUTINE close_ps                           
      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE plot_filled_contours_adapt(file_unit,ne,el_type,el_in,el_area,elxy,xyplt,pplt,nptri,npplt,rect,r,s,snap,fig)
      
      USE transformation, ONLY: init_vandermonde,element_transformation,xy2rs
      USE basis, ONLY: element_basis,linear_basis           
      USE read_dginp, ONLY: p,hbp,ctp   
      USE plot_globals, ONLY: Z,hb,Qx,Qy      
      USE area_qpts_mod, ONLY: tri_cubature,area_qpts
      USE shape_functions_mod, ONLY: shape_functions_area_eval
      USE evaluate_mod, ONLY: evaluate_solution
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: file_unit
      INTEGER, INTENT(IN) :: ne
      INTEGER, DIMENSION(:), INTENT(IN) :: el_type
      INTEGER, DIMENSION(:), INTENT(IN) :: el_in
      REAL(rp), DIMENSION(:), INTENT(IN) :: el_area
      REAL(rp), DIMENSION(:,:,:), INTENT(IN) :: elxy
      REAL(rp), DIMENSION(:,:,:), INTENT(INOUT) :: xyplt 
      INTEGER, DIMENSION(:), INTENT(IN) :: pplt 
      INTEGER, DIMENSION(:), INTENT(IN) :: nptri
      INTEGER, DIMENSION(:), INTENT(IN) :: npplt      
      INTEGER, DIMENSION(:,:,:), INTENT(IN) :: rect 
      REAL(rp), DIMENSION(:,:), INTENT(IN) :: r
      REAL(rp), DIMENSION(:,:), INTENT(IN) :: s
      INTEGER, INTENT(IN) :: snap      
      TYPE(plot_type), INTENT(INOUT) :: fig

      
      INTEGER :: i,j,v
      INTEGER :: el,nd,dof,lev,tri,ord,pt
      INTEGER :: et,nv,pl    
      INTEGER :: nqpt,nnd,mnnds,mnp,ndf,mnqpta,nqpta(nel_type)
      INTEGER :: qpt_order
      INTEGER :: nptri_total
      INTEGER :: ne_total
      INTEGER :: pplt_max
      REAL(rp) :: sol_lev
      REAL(rp) :: dc      
      REAL(rp), DIMENSION(:,:), ALLOCATABLE :: qpt
      REAL(rp), DIMENSION(:), ALLOCATABLE :: wpt
      REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: qpta
      REAL(rp), DIMENSION(:,:), ALLOCATABLE ::  wpta
      REAL(rp), DIMENSION(:,:), ALLOCATABLE :: l
      REAL(rp), DIMENSION(:), ALLOCATABLE :: xpt,ypt
      REAL(rp), DIMENSION(:), ALLOCATABLE :: rpt,spt
      REAL(rp), DIMENSION(:), ALLOCATABLE :: sol_lin,sol_el  
      REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: phia
      REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: psi,dpsidr,dpsids      
      REAL(rp) :: xpta,ypta
      REAL(rp) :: drdx,drdy,dsdx,dsdy      
      REAL(rp) :: qpts(13**2,3)
      REAL(rp) :: xv(3),yv(3)
      REAL(rp) :: detJ
      REAL(rp) :: H,u_vel,v_vel
      REAL(rp) :: sol_avg
      REAL(rp) :: error_total
   
      REAL(rp) :: color_val(3)
      REAL(rp) :: err,max_err
!       REAL(rp) :: rel_tol,abs_tol
      
      REAL(rp), DIMENSION(:), ALLOCATABLE :: hb_val,zeta_val,Qx_val,Qy_val     
      REAL(rp), DIMENSION(:), ALLOCATABLE :: hb_sol,zeta_sol,Qx_sol,Qy_sol    
      
      
      ! quadrature points for elemental solution average
      CALL area_qpts(1,p,ctp,nel_type,nqpta,mnqpta,wpta,qpta)      
      ALLOCATE(phia((p+1)**2,mnqpta,nel_type))

      ! quadrature points for error integral
      qpt_order = 2
      CALL tri_cubature(qpt_order,nqpt,qpts)      
      mnqpta = max(mnqpta,nqpt)
      
      ALLOCATE(hb_val(mnpp),zeta_val(mnpp),Qx_val(mnpp),Qy_val(mnpp))
      ALLOCATE(hb_sol(mnqpta),zeta_sol(mnqpta),Qx_sol(mnqpta),Qy_sol(mnqpta))      
      ALLOCATE(qpt(mnqpta,2),wpt(mnqpta))
      ALLOCATE(l(3,mnqpta))
      ALLOCATE(xpt(mnqpta),ypt(mnqpta))
      ALLOCATE(rpt(mnqpta),spt(mnqpta))
      ALLOCATE(sol_lin(mnqpta),sol_el(mnqpta))
      
      DO pt = 1,nqpt
        qpt(pt,1) = qpts(pt,1)
        qpt(pt,2) = qpts(pt,2)
        wpt(pt) = qpts(pt,3)
      ENDDO
      
      
      ! shape functions for elemental solution average
      mnp = maxval(np)      
      mnnds = (mnp+1)**2     
      ALLOCATE(psi(mnnds,mnqpta,nel_type),dpsidr(mnnds,mnqpta,nel_type),dpsids(mnnds,mnqpta,nel_type))
      
      DO et = 1,nel_type         
        CALL shape_functions_area_eval(et,np(et),nnds(et),nqpta(et),qpta(:,1,et),qpta(:,2,et), &
                                       psi(:,:,et),dpsidr(:,:,et),dpsids(:,:,et))
      ENDDO      
      
      
      
      
      dc = (fig%sol_max-fig%sol_min)/real(ncolors-1,rp)      
      
      CALL init_vandermonde(nel_type,np)      
            
!       rel_tol = 1d-1  
!       abs_tol = 1d-1
      
      

      nptri_total = 0
      error_total = 0d0
      ne_total = 0
      pplt_max = 0
      
 elem:DO el = 1,ne

        IF (el_in(el) == 0) THEN
          CYCLE elem
        ENDIF
        
        et = el_type(el)       
        nnd = nnds(et)        
          
        IF (et > 2) THEN
         pl = pc
        ELSE
          pl = 1
        ENDIF
        
        IF (adapt_option == 1) THEN           
          ! Find element average value for solution
          CALL evaluate_solution(el,et,fig%type_flag,snap,sol_el,nqpta(et),qpta(:,1,et),qpta(:,2,et))
        
          sol_avg = 0d0
          DO pt = 1,nqpta(et)
            CALL element_transformation(nnd,elxy(:,el,1),elxy(:,el,2),psi(:,pt,et),xpta,ypta, &
                                        dpsidr(:,pt,et),dpsids(:,pt,et),drdx,drdy,dsdx,dsdy,detJ)                
            sol_avg = sol_avg + wpta(pt,et)*sol_el(pt)*detJ
          ENDDO
          sol_avg = sol_avg/el_area(el)
             
          PRINT*, ""        
          PRINT "(A,I9,A,I9,A,F15.7)", "ELEMENT: ", el, "    et: ", et, "    sol_avg = ", sol_avg        
        ENDIF
         
 order: DO ord = pl,nord
 
          i = (et-1)*nord+ord
          

          DO pt = 1,npplt(i)              
            CALL element_transformation(nnd,elxy(:,el,1),elxy(:,el,2),psic(:,pt,i),xyplt(pt,el,1),xyplt(pt,el,2))           
          ENDDO                  
     
          ! Evaluate solution at plotting nodes
          CALL evaluate_solution(el,et,fig%type_flag,snap,fig%sol_val(:,el),npplt(i),phi_hb=phi_hb(:,:,i),phi_sol=phi_sol(:,:,i),ndf_hb=ndof_hb(et),ndf_sol=ndof_sol(et))          
!           CALL evaluate_solution(el,et,fig%type_flag,snap,fig%sol_val(:,el),npplt(i),r(:,i),s(:,i))       

          IF (adapt_option == 0) THEN
            EXIT order
          ENDIF
    
          
          CALL linear_basis(ndf,nqpt,qpt(:,1),qpt(:,2),l)          
          
          err = 0d0
          max_err = 0d0
  sub_tri:DO tri = 1,nptri(i)
          
            ! Evaluate linear solution at plotting element quadrature points  
            sol_lin = 0d0
            DO pt = 1,nqpt
              DO v = 1,3              
                 nd = rect(v,tri,i)
                sol_lin(pt) = sol_lin(pt) + l(v,pt)*fig%sol_val(nd,el)
              ENDDO
            ENDDO
            
            DO pt = 1,nqpt
              xpt(pt) = 0d0
              ypt(pt) = 0d0
              DO v = 1,3
                nd = rect(v,tri,i)
                
                xv(v) = xyplt(nd,el,1)
                yv(v) = xyplt(nd,el,2)
                
                xpt(pt) = xpt(pt) + l(v,pt)*xv(v)
                ypt(pt) = ypt(pt) + l(v,pt)*yv(v)
              ENDDO
            ENDDO
            
            detJ = .25d0*((xv(2)-xv(1))*(yv(3)-yv(1))-(xv(3)-xv(1))*(yv(2)-yv(1)))
            
            CALL xy2rs(et,np(et),elxy(:,el,1),elxy(:,el,2),nqpt,xpt,ypt,rpt,spt)   
            
            ! Evaluate DG solution at plotting element quadrature points
            CALL evaluate_solution(el,et,fig%type_flag,snap,sol_el,nqpt,rpt,spt)

            
            
            DO pt = 1,nqpt                        
              err = err + detJ*wpt(pt)*(sol_el(pt)-sol_lin(pt))**2
            ENDDO
 
              max_err = sqrt(err)

!             DO pt = 1,nqpt               
!               err = abs(sol_el(pt)-sol_lin(pt))
!               IF (err > max_err) THEN
!                 max_err = err
!               ENDIF
!             ENDDO
            
            IF (abs(max_err/sol_avg) > fig%rel_tol .and. max_err > fig%abs_tol) THEN
              PRINT "(A,I4,A,I4)", "    Early sub_tri exit: ", tri, "/", nptri(i)              
              EXIT sub_tri
            ENDIF
            
          ENDDO sub_tri
          
          err = sqrt(err)
!           err = max_err
          
          PRINT "(A,I9,A,F15.7,A,F15.7)", "  ORDER: ", pplt(ord), "  Abs. Error = ", err, "  Rel. Error = ", abs(err/sol_avg)         
          
          IF (abs(err/sol_avg) < fig%rel_tol .or. err < fig%abs_tol) THEN
            fig%el_plt(el) = i       
            EXIT order     
          ENDIF          
          
        ENDDO order                 
          
                
        
        ord = fig%el_plt(el)
        
         
        
        DO tri = 1,nptri(ord)

          DO v = 1,3  
   
            nd = rect(v,tri,ord)
          
            sol_lev = fig%sol_min
            lev = ncolors
    levels: DO i = 1,ncolors-1
              IF ((fig%sol_val(nd,el) >= sol_lev) .and. (fig%sol_val(nd,el) < sol_lev+dc)) THEN
                lev = i
                CALL interp_colors(lev,sol_lev,dc,colors,fig%sol_val(nd,el),color_val)
                EXIT levels
              ENDIF
              sol_lev = sol_lev + dc
            ENDDO levels
          
            IF (fig%sol_val(nd,el) <= fig%sol_min) THEN
              lev = 1
              color_val(1) = colors(lev,1)
              color_val(2) = colors(lev,2)
              color_val(3) = colors(lev,3)
            ELSE IF (fig%sol_val(nd,el) > fig%sol_max) THEN
              lev = ncolors
              color_val(1) = colors(lev,1)
              color_val(2) = colors(lev,2)
              color_val(3) = colors(lev,3)           
            ENDIF          

          
            WRITE(file_unit,"(A,5(F9.5,1x),A)") "[",ax*xyplt(nd,el,1)+bx,ay*xyplt(nd,el,2)+by,color_val(1),color_val(2),color_val(3),"]"  
          ENDDO        

          WRITE(file_unit,"(A)") "trifill"        
        
        ENDDO 
        
        error_total = error_total + err**2
        nptri_total = nptri_total + nptri(ord)
        ne_total = ne_total + 1
        
        IF (pplt(ord) > pplt_max) THEN
          pplt_max = pplt(ord)
        ENDIF

      ENDDO elem
      
      IF (adapt_option == 1) THEN
        WRITE(998,"(A25,I25,E24.17,1x,I25,I25,I25)") fig%name, snap, sqrt(error_total), nptri_total, pplt_max, ne_total
      ENDIF
              

      END SUBROUTINE plot_filled_contours_adapt
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   

      SUBROUTINE plot_line_contours(file_unit,ne,el_type,el_in,nptri,rect,xyplt,rre,sre,snap,fig)
      
      USE globals, ONLY: mndof,elxy,xy,ect,np,mnnds
      USE plot_globals, ONLY: Z,hb,Qx,Qy
      USE read_dginp, ONLY: p,hbp
      USE basis, ONLY: element_basis,linear_basis,dgswem_basis
      USE transformation, ONLY: element_transformation
      USE shape_functions_mod, ONLY: shape_functions_area_eval, shape_functions_edge_eval
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: file_unit
      INTEGER, INTENT(IN) :: ne
      INTEGER, DIMENSION(:), INTENT(IN) :: el_type
      INTEGER, DIMENSION(:), INTENT(IN) :: el_in
      INTEGER, DIMENSION(:), INTENT(IN) :: nptri
      INTEGER, DIMENSION(:,:,:), INTENT(IN) :: rect
      REAL(rp), DIMENSION(:,:,:), INTENT(IN) :: xyplt  
      REAL(rp), DIMENSION(:,:), INTENT(IN) :: rre
      REAL(rp), DIMENSION(:,:), INTENT(IN) :: sre
      INTEGER, INTENT(IN) :: snap
      TYPE(plot_type), INTENT(IN) :: fig

      
      INTEGER :: i,j,vrt,it
      INTEGER :: el,nd,dof,lev,tri,ig
      INTEGER :: et,nv,ndf,ntry,nnd,ord    
      REAL(rp) :: sol_lev
      REAL(rp) :: dc
      INTEGER :: above(3)
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: el_below
      INTEGER :: n1,n2
      INTEGER :: nd1,nd2
      REAL(rp) :: x(2),y(2)
      REAL(rp) :: r,s,re(1),se(1)
      REAL(rp) :: r1,r2,s1,s2
      REAL(rp) :: t,f,dfdr,dfds,dfdt,dt,t0
      INTEGER :: max_it
      INTEGER :: fail_flag
      REAL(rp) :: phi(mndof,2),dpdr(mndof,1),dpds(mndof,1)
      REAL(rp) :: l(mnnds,1)
      REAL(rp) :: xv(3),yv(3)
      REAL(rp) :: rv(3),sv(3)
      REAL(rp) :: zeta,bathy,xmom,ymom,H
      REAL(rp) :: dzdr,dzds,dbdr,dbds,dHdr,dHds
      REAL(rp) :: dxmdr,dxmds,dymdr,dymds
      REAL(rp) :: uvel,vvel
      REAL(rp) :: tol
      REAL(rp) :: hbe(2),Ze(2),ve(2),fe(2)
      REAL(rp) :: e,a,b,c,w,u,v
      
      tol = 1d-8
      ntry = 20
      max_it = 1000
      


      ALLOCATE(el_below(maxval(nptri),ne))
      el_below = 0
      
       dc = (fig%sol_max-fig%sol_min)/real(nctick-1,rp)                  
       sol_lev = fig%sol_min      
       
levels:DO lev = 1,nctick
      
    elem:DO el = 1,ne
           et = el_type(el)
           IF (el_in(el) == 0) THEN
             CYCLE elem
           ENDIF
     
           ord = fig%el_plt(el)
           
    subtri:DO tri = 1,nptri(ord)
    
             IF (el_below(tri,el) == 1) THEN
               CYCLE subtri
             ENDIF 

!              PRINT("(3(A,I9))"), "lev = ", lev, " el = ", el, " tri = ",tri
             
             DO vrt = 1,3  
   
               nd = rect(vrt,tri,ord)
          
               IF (fig%sol_val(nd,el) >= sol_lev) THEN
                 above(vrt) = 1
               ELSE
                 above(vrt) = 0
               ENDIF
                
             ENDDO        
 
 
             IF (above(1) == 1 .and. above(2) == 1 .and. above(3) == 1) THEN
               CYCLE subtri
             ENDIF
             
             IF (above(1) == 0 .and. above(2) == 0 .and. above(3) == 0) THEN
               el_below(tri,el) = 1             
               CYCLE subtri
             ENDIF       
             
             
             DO i = 1,3
               nd = rect(i,tri,ord)
               rv(i) = rre(nd,ord)
               sv(i) = sre(nd,ord)
             ENDDO

             
             i = 0             
             
       edge: DO vrt = 1,3
               n1 = mod(vrt+0,3)+1
               n2 = mod(vrt+1,3)+1
                              
               IF ((above(n1) == 1 .and. above(n2) == 1) .or. &
                   (above(n1) == 0 .and. above(n2) == 0)) THEN    
                   CYCLE edge
               ENDIF
                   
               r1 = rv(n1)
               r2 = rv(n2)
               s1 = sv(n1)
               s2 = sv(n2)   
                 
               dt = 1d0/real(ntry,rp)
               t0 = 0d0
               fail_flag = 1
      guesses: DO ig = 1,2*ntry+1
                 
                 t = t0                 
                 it = 0
         newton: DO 
         
                   r = .5d0*((1d0-t)*r1 + (1d0+t)*r2)
                   s = .5d0*((1d0-t)*s1 + (1d0+t)*s2)
                   re(1) = r
                   se(1) = s            
                   
                   IF (fig%type_flag == 2 .or. fig%type_flag == 4) THEN
#ifdef adcirc                   
                     CALL linear_basis(ndf,1,re,se,phi,dpdr,dpds)
#elif dgswem                       
                     CALL dgswem_basis(hbp,ndf,1,re,se,phi,dpdr,dpds)   
#else
                     CALL element_basis(et,hbp,ndf,1,re,se,phi,dpdr,dpds)   
#endif                     
                     bathy = 0d0
                     dbdr = 0d0
                     dbds = 0d0
                     DO dof = 1,ndf
                       bathy = bathy + phi(dof,1)*hb(dof,el,1)
                       dbdr = dbdr + dpdr(dof,1)*hb(dof,el,1)
                       dbds = dbds + dpds(dof,1)*hb(dof,el,1)
                     ENDDO
                     IF (fig%type_flag == 2) THEN
                       f = bathy - sol_lev
                       dfdr = dbdr
                       dfds = dbds
                     ENDIF
                   ENDIF
                   
                   IF (fig%type_flag == 3 .or. fig%type_flag == 4) THEN
#ifdef adcirc
                     CALL linear_basis(ndf,1,re,se,phi,dpdr,dpds)
#elif dgswem
                     CALL dgswem_basis(p,ndf,1,re,se,phi,dpdr,dpds)    
#else                     
                     CALL element_basis(et,p,ndf,1,re,se,phi,dpdr,dpds)    
#endif                     
                     zeta = 0d0
                     dzdr = 0d0
                     dzds = 0d0
                     DO dof = 1,ndf
                       zeta = zeta + phi(dof,1)*Z(dof,el,snap)
                       dzdr = dzdr + dpdr(dof,1)*Z(dof,el,snap)
                       dzds = dzds + dpds(dof,1)*Z(dof,el,snap)
                     ENDDO
                     IF (fig%type_flag == 3) THEN
                       f = zeta - sol_lev
                       dfdr = dzdr
                       dfds = dzds
                       ENDIF
                   ENDIF    
                     
                   IF (fig%type_flag == 4) THEN
                     xmom = 0d0
                     ymom = 0d0
                     dxmdr = 0d0
                     dxmds = 0d0
                     dymdr = 0d0
                     dymds = 0d0
                     DO dof = 1,ndf
                       xmom = xmom + phi(dof,1)*Qx(dof,el,snap)
                       ymom = ymom + phi(dof,1)*Qy(dof,el,snap)
                       dxmdr = dxmdr + dpdr(dof,1)*Qx(dof,el,snap)
                       dxmds = dxmds + dpds(dof,1)*Qx(dof,el,snap)
                       dymdr = dymdr + dpdr(dof,1)*Qy(dof,el,snap)
                       dymds = dymds + dpds(dof,1)*Qy(dof,el,snap)
                     ENDDO
                     H = zeta + bathy          
                     dHdr = dzdr + dbdr
                     dHds = dzds + dbds
                     uvel = xmom/H
                     vvel = ymom/H
                     f = uvel**2 + vvel**2 - sol_lev**2
                     dfdr = 2d0*(uvel*(dxmdr*H-dHdr*xmom)/H**2 + vvel*(dymdr*H-dHdr*ymom)/H**2)
                     dfds = 2d0*(uvel*(dxmds*H-dHds*xmom)/H**2 + vvel*(dymds*H-dHds*ymom)/H**2)
                   ENDIF
                   
                   
                   dfdt = .5d0*(dfdr*(r2-r1) + dfds*(s2-s1))
                 
                   t = t - f/dfdt
                   it = it + 1
                   
                   
                   IF (abs(f) < tol) THEN


                     IF (t <= 1d0+tol .and. t >= -(1d0+tol)) THEN
                       r = .5d0*((1d0-t)*r1 + (1d0+t)*r2)
                       s = .5d0*((1d0-t)*s1 + (1d0+t)*s2)                   
                       fail_flag = 0           
!                        PRINT*, "Iteration successful"
                       EXIT guesses
                     ELSE                 
!                        PRINT("(4(A,I9))"), "  Point outside edge     lev = ", lev, " el = ", el, " tri = ",tri, " vrt = ",vrt
!                        PRINT("(A,F12.6,A,I5,A,E24.17,A,E24.17)"), "      t0 = ", t0, "  it = ", it,"  t = ", t, "  f = ", f
!                        PRINT*, ""     
                       EXIT newton
                     ENDIF
                   ENDIF
                   
                   IF (it > max_it) THEN     
!                    PRINT*, "  Max iterations exceeded"
                     EXIT newton
                   ENDIF
                   
                 ENDDO newton
                 
               t0 = t0*-1d0                 
               IF (mod(ig,2) == 1) THEN
                 t0 = t0 + dt
               ENDIF
                 
               ENDDO guesses
                 
               i = i + 1
               IF (fail_flag == 0) THEN
                 re(1) = r
                 se(1) = s
               ELSE                
                 PRINT*, "  iteration failed ", "r = ",r, " s = ",s 
                 re(1) = .5d0*(r1+r2)
                 se(1) = .5d0*(s1+s2)
               ENDIF
 
               CALL shape_functions_area_eval(et,np(et),nnd,1,re,se,l)
               CALL element_transformation(nnd,elxy(:,el,1),elxy(:,el,2),l(:,1),x(i),y(i))  
                 
 
               !!!! Edge midpoint !!!!
!                  i = i + 1
!                  nd1 = rect(n1,tri,et)
!                  nd2 = rect(n2,tri,et)                 
!                  x(i) = 0.5d0*(xyplt(nd1,el,1) + xyplt(nd2,el,1))
!                  y(i) = 0.5d0*(xyplt(nd1,el,2) + xyplt(nd2,el,2))                   
                 
             ENDDO edge

             IF (i > 2) THEN
               PRINT*, "Error: "
               STOP
             ENDIF 
             
             WRITE(file_unit,"(2(F9.5,1x))") ax*x(1)+bx,ay*y(1)+by
             WRITE(file_unit,"(2(F9.5,1x))") ax*x(2)+bx,ay*y(2)+by     
             WRITE(file_unit,"(A)") "draw-line"  
             
           ENDDO subtri
         ENDDO elem
         sol_lev = sol_lev + dc
      ENDDO levels

      END SUBROUTINE plot_line_contours

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      SUBROUTINE plot_mesh(file_unit,ne,nverts,el_plt,pplt,el_type,el_in,xy,ect,xyplt)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: file_unit
      INTEGER, INTENT(IN) :: ne
      INTEGER, DIMENSION(:), INTENT(IN) :: nverts
      INTEGER, DIMENSION(:), INTENT(IN) :: el_plt
      INTEGER, DIMENSION(:), INTENT(IN) :: pplt
      INTEGER, DIMENSION(:), INTENT(IN) :: el_type
      INTEGER, DIMENSION(:), INTENT(IN) :: el_in      
      REAL(rp), DIMENSION(:,:), INTENT(IN) :: xy
      INTEGER, DIMENSION(:,:), INTENT(IN) :: ect
      REAL(rp), DIMENSION(:,:,:), INTENT(IN) :: xyplt

      INTEGER :: el,nd
      INTEGER :: et

      
 elem:DO el = 1,ne
        
        IF (el_in(el) == 0) THEN
          CYCLE elem
        ENDIF
        
        et = el_type(el)
        
        IF (et == 1) THEN
          DO nd = 1,nverts(et)
            WRITE(file_unit,"(2(F9.5,1x))") ax*xy(1,ect(nd,el))+bx,ay*xy(2,ect(nd,el))+by    
          ENDDO                
          WRITE(file_unit,"(A)") "draw-tri-element"         
        ELSE IF (et == 2) THEN
          DO nd = 1,nverts(et)
            WRITE(file_unit,"(2(F9.5,1x))") ax*xy(1,ect(nd,el))+bx,ay*xy(2,ect(nd,el))+by    
          ENDDO                
          WRITE(file_unit,"(A)") "draw-quad-element"  
        ELSE
          WRITE(file_unit,"(A)") "newpath"
          WRITE(file_unit,"(2(F9.5,1x),A)") ax*xyplt(1,el,1)+bx,ay*xyplt(1,el,2)+by,"moveto" 
          DO nd = 2,nverts(et)*pplt(el_plt(el))
            WRITE(file_unit,"(2(F9.5,1x),A)") ax*xyplt(nd,el,1)+bx,ay*xyplt(nd,el,2)+by, "lineto"
          ENDDO
          WRITE(file_unit,"(A)") "closepath"
          WRITE(file_unit,"(A)") ".5 setlinewidth 2 setlinejoin"          
          WRITE(file_unit,"(A)") "stroke"          
        ENDIF
      ENDDO elem

      END SUBROUTINE plot_mesh
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      SUBROUTINE plot_vis_mesh(file_unit,ne,el_type,el_in,xyplt,nptri,rect,fig)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: file_unit
      INTEGER, INTENT(IN) :: ne
      INTEGER, DIMENSION(:), INTENT(IN) :: el_type
      INTEGER, DIMENSION(:), INTENT(IN) :: el_in      
      REAL(rp), DIMENSION(:,:,:), INTENT(IN) :: xyplt
      INTEGER, DIMENSION(:), INTENT(IN) :: nptri
      INTEGER, DIMENSION(:,:,:), INTENT(IN) :: rect
      TYPE(plot_type), INTENT(IN) :: fig      

      INTEGER :: el,nd,tri,v
      INTEGER :: et,ord      

      
 elem:DO el = 1,ne
        
        IF (el_in(el) == 0) THEN
          CYCLE elem
        ENDIF
        
        et = el_type(el)
        ord = fig%el_plt(el)        

        DO tri = 1,nptri(ord)        
          DO v = 1,3        
            nd = rect(v,tri,ord)
            WRITE(file_unit,"(2(F9.5,1x))") ax*xyplt(nd,el,1)+bx,ay*xyplt(nd,el,2)+by
          ENDDO                
          WRITE(file_unit,"(A)") "draw-tri-element"            
        ENDDO
      ENDDO elem

      END SUBROUTINE plot_vis_mesh      

      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   

     SUBROUTINE plot_boundaries(file_unit,nverts,el_plt,pplt,nbed,bedn,ged2el,ged2led,el_type,el_in,ect,xy,xyplt)
     
     IMPLICIT NONE
     
     INTEGER, INTENT(IN) :: file_unit
     INTEGER, DIMENSION(:), INTENT(IN) :: nverts
     INTEGER, DIMENSION(:), INTENT(IN) :: el_plt
     INTEGER, DIMENSION(:), INTENT(IN) :: pplt
     INTEGER, INTENT(IN) :: nbed
     INTEGER, DIMENSION(:), INTENT(IN) :: bedn
     INTEGER, DIMENSION(:,:), INTENT(IN) :: ged2el
     INTEGER, DIMENSION(:,:), INTENT(IN) :: ged2led
     INTEGER, DIMENSION(:), INTENT(IN) :: el_type
     INTEGER, DIMENSION(:), INTENT(IN) :: el_in
     INTEGER, DIMENSION(:,:), INTENT(IN) :: ect
     REAL(rp), DIMENSION(:,:), INTENT(IN) :: xy
     REAL(rp), DIMENSION(:,:,:), INTENT(IN) :: xyplt
     
     INTEGER :: ed,nd,j
     INTEGER :: el,et,nv
     INTEGER :: ged,led
     INTEGER :: n1,n2
     
edge:DO ed = 1,nbed
       ged = bedn(ed)
       
       el = ged2el(1,ged)
       led = ged2led(1,ged)
       
       et = el_type(el)
       nv = nverts(et)
       
       IF (el_in(el) == 0) THEN
         CYCLE
       ENDIF
       
       IF (et <= 2) THEN
       
          n1 = mod(led+0,nv)+1
          n2 = mod(led+1,nv)+1
          
          WRITE(file_unit,"(2(F9.5,1x))") ax*xy(1,ect(n1,el))+bx,ay*xy(2,ect(n1,el))+by
          WRITE(file_unit,"(2(F9.5,1x))") ax*xy(1,ect(n2,el))+bx,ay*xy(2,ect(n2,el))+by
          WRITE(file_unit,"(A)") "draw-line"          
        ELSE
          
          n1 = mod(led,nv)*pplt(el_plt(el)) + 1
          n2 = n1 + pplt(el_plt(el))
          WRITE(file_unit,"(A)") "newpath"
          WRITE(file_unit,"(2(F9.5,1x),A)") ax*xyplt(n1,el,1)+bx,ay*xyplt(n1,el,2)+by,"moveto" 
          DO j = 1, pplt(el_plt(el))
            nd = n1 + j
            IF (nd == nv*pplt(el_plt(el))+1) THEN
              nd = 1
            ENDIF
            WRITE(file_unit,"(2(F9.5,1x),A)") ax*xyplt(nd,el,1)+bx,ay*xyplt(nd,el,2)+by, "lineto"
          ENDDO
          WRITE(file_unit,"(A)") ".5 setlinewidth 2 setlinejoin"          
          WRITE(file_unit,"(A)") "stroke"              
        ENDIF
       
     ENDDO edge
     
     END SUBROUTINE plot_boundaries
     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      SUBROUTINE plot_cb_nodes(file_unit,ctp,nbou,fbseg,fbnds,xy,bndxy)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: file_unit
      INTEGER, INTENT(IN) :: ctp
      INTEGER, INTENT(IN) :: nbou
      INTEGER, DIMENSION(:,:), INTENT(IN) :: fbseg
      INTEGER, DIMENSION(:,:), INTENT(IN) :: fbnds
      REAL(rp), DIMENSION(:,:), INTENT(IN) :: xy
      REAL(rp), DIMENSION(:,:,:,:), INTENT(IN) :: bndxy
      
      INTEGER :: i,j,k
      INTEGER :: nd
      INTEGER :: nbnds
      INTEGER :: btype
      
      DO i = 1,nbou
        nbnds = fbseg(1,i)
        btype = fbseg(2,i)
        
        IF( btype == 0 .OR. btype == 10 .OR. btype == 20  .OR. &   ! land boundaries
            btype == 1 .OR. btype == 11 .OR. btype == 21 ) THEN    ! island boundaries
            
          DO j = 1,nbnds-1
            
            nd = fbnds(j,i)            
            WRITE(file_unit,"(3(I5,1x))") 0,0,0                    
            WRITE(file_unit,"(2(F9.5,1x))") ax*xy(1,nd)+bx, ay*xy(2,nd)+by
            WRITE(file_unit,"(A)") "draw-dot"
            
            DO k = 1,ctp-1
              WRITE(file_unit,"(3(I5,1x))") 0,0,0             
              WRITE(file_unit,"(2(F9.5,1x))") ax*bndxy(1,k,j,i)+bx, ay*bndxy(2,k,j,i)+by
              WRITE(file_unit,"(A)") "draw-dot"              
            ENDDO
          
          ENDDO
          
          WRITE(file_unit,"(3(I5,1x))") 0,0,0           
          WRITE(file_unit,"(2(F9.5,1x))") ax*xy(1,nbnds)+bx, ay*xy(2,nbnds)+by
          WRITE(file_unit,"(A)") "draw-dot"         
            
        ENDIF    
      ENDDO
      
      RETURN
      END SUBROUTINE plot_cb_nodes
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      SUBROUTINE plot_elxy_nodes(file_unit,ne,el_type,nnds,elxy)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: file_unit    
      INTEGER, INTENT(IN) :: ne
      INTEGER, DIMENSION(:), INTENT(IN) :: el_type
      INTEGER, DIMENSION(:), INTENT(IN) :: nnds
      REAL(rp), DIMENSION(:,:,:), INTENT(IN) :: elxy
      
      INTEGER :: el,nd
      INTEGER :: et
      INTEGER :: nnd

      DO el = 1,ne
        et = el_type(el)
        nnd = nnds(et)
        
        DO nd = 1,nnd
          WRITE(file_unit,"(3(I5,1x))") 0,0,0         
          WRITE(file_unit,"(2(F9.5,1x))") ax*elxy(nd,el,1)+bx, ay*elxy(nd,el,2)+by
          WRITE(file_unit,"(A)") "draw-dot"         
        ENDDO
      ENDDO

      
      RETURN
      END SUBROUTINE plot_elxy_nodes      
      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      SUBROUTINE fill_elements(file_unit,ne,nverts,el_plt,pplt,el_type,el_in,xy,ect,xyplt)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: file_unit
      INTEGER, INTENT(IN) :: ne
      INTEGER, DIMENSION(:), INTENT(IN) :: nverts
      INTEGER, DIMENSION(:), INTENT(IN) :: el_plt
      INTEGER, DIMENSION(:), INTENT(IN) :: pplt
      INTEGER, DIMENSION(:), INTENT(IN) :: el_type
      INTEGER, DIMENSION(:), INTENT(IN) :: el_in      
      REAL(rp), DIMENSION(:,:), INTENT(IN) :: xy
      INTEGER, DIMENSION(:,:), INTENT(IN) :: ect
      REAL(rp), DIMENSION(:,:,:), INTENT(IN) :: xyplt

      INTEGER :: el,nd
      INTEGER :: et
      INTEGER :: i,n
      INTEGER :: fill
      INTEGER, ALLOCATABLE, DIMENSION(:) :: fill_list
      LOGICAL :: file_exists
      
      ALLOCATE(fill_list(ne))
      
!       fill_list = 0
!       INQUIRE(file='element.fill',exist=file_exists)
!       IF (file_exists) THEN
!         OPEN(unit=101,file='element.fill')
!         READ(101,*) n
!         DO i = 1,n
!           READ(101,*) el
!           fill_list(el) = 1
!         ENDDO
!         CLOSE(101)
!       ENDIF

      fill_list = 1
      
 elem:DO el = 1,ne
        
        IF (el_in(el) == 0) THEN
          CYCLE elem
        ENDIF
        
        et = el_type(el)
        fill = fill_list(el)        
        
        IF (et == 1) THEN
          IF (fill == 1) THEN
            DO nd = 1,nverts(et)
              WRITE(file_unit,"(2(F9.5,1x))") ax*xy(1,ect(nd,el))+bx,ay*xy(2,ect(nd,el))+by
            ENDDO                
          WRITE(file_unit,"(A)") "fill-tri-element"            
          ENDIF
        ELSE IF (et == 2) THEN
          IF (fill == 1) THEN
            DO nd = 1,nverts(et)
              WRITE(file_unit,"(2(F9.5,1x))") ax*xy(1,ect(nd,el))+bx,ay*xy(2,ect(nd,el))+by
            ENDDO                
          WRITE(file_unit,"(A)") "fill-quad-element"            
          ENDIF
       ELSE
          IF (fill == 1) THEN       
            WRITE(file_unit,"(A)") "newpath"
            WRITE(file_unit,"(2(F9.5,1x),A)") ax*xyplt(1,el,1)+bx,ay*xyplt(1,el,2)+by,"moveto" 
            DO nd = 2,nverts(et)*pplt(el_plt(el))
              WRITE(file_unit,"(2(F9.5,1x),A)") ax*xyplt(nd,el,1)+bx,ay*xyplt(nd,el,2)+by, "lineto"
            ENDDO
!             WRITE(file_unit,"(A)") "gsave 0 0 1 setrgbcolor fill grestore"
            WRITE(file_unit,"(A)") "gsave 1 1 1 setrgbcolor fill grestore"
          ENDIF        
        ENDIF
      ENDDO elem

      END SUBROUTINE fill_elements     

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   

      SUBROUTINE plot_stations(file_unit,sta_start,sta_end,xysta)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: file_unit
      INTEGER, INTENT(IN) :: sta_start
      INTEGER, INTENT(IN) :: sta_end
      REAL(rp), DIMENSION(:,:), INTENT(IN) :: xysta
      
      INTEGER :: sta
      REAL(rp) :: x,y
      
      DO sta = sta_start,sta_end
            
        x = xysta(1,sta)
        y = xysta(2,sta)
            
        IF (x > xbox_min .AND. x < xbox_max .AND. y > ybox_min .AND. y < ybox_max) THEN  
          WRITE(file_unit,"(3(I5,1x))") 1,0,0        
          WRITE(file_unit,"(2(F9.5,1x))") ax*x+bx,ay*y+by
          WRITE(file_unit,"(A)") "draw-dot"        
        ENDIF
      
      ENDDO
      
      END SUBROUTINE plot_stations

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      SUBROUTINE plot_xy(file_unit,npt,xvec,yvec,rgb,width)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: file_unit
      INTEGER, INTENT(IN) :: npt
      REAL(rp), DIMENSION(:), INTENT(IN) :: xvec
      REAL(rp), DIMENSION(:), INTENT(IN) :: yvec
      REAL(rp), DIMENSION(:), INTENT(IN) :: rgb
      REAL(rp), INTENT(IN) :: width
      
      INTEGER :: pt
      
      WRITE(file_unit,"(3(F9.5,1x),A)") rgb(1),rgb(2),rgb(3)," setrgbcolor"      
      WRITE(file_unit,"(A)") "newpath"        
      WRITE(file_unit,"(2(F9.5,1x),A)") ax*xvec(1)+bx,ay*yvec(1)+by, " moveto"      
      DO pt = 2,npt                
        WRITE(file_unit,"(2(F9.5,1x),A)") ax*xvec(pt)+bx,ay*yvec(pt)+by, " lineto"      
      ENDDO
      WRITE(file_unit,"(F9.5,1x,A)") width," setlinewidth 2 setlinejoin"
      WRITE(file_unit,"(A)") "stroke"            
      WRITE(file_unit,"(A)")  "0 0 0 setrgbcolor"       
      
      RETURN
      END SUBROUTINE plot_xy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      SUBROUTINE interp_colors(lev,sol_lev,dc,colors,sol_val,color_val)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: lev
      REAL(rp), INTENT(IN) :: sol_lev
      REAL(rp), INTENT(IN) :: dc
      REAL(rp), DIMENSION(:,:), INTENT(IN) :: colors
      REAL(rp), INTENT(IN) :: sol_val
      REAL(rp), DIMENSION(:), INTENT(OUT) :: color_val
      
      INTEGER :: i
      REAL(rp) :: l1,l2
      
      l1 = (sol_val-(sol_lev+dc))/(-dc)
      l2 = (sol_val-sol_lev)/dc
      
      DO i = 1,3
        color_val(i) = l1*colors(lev,i) + l2*colors(lev+1,i)
      ENDDO      
      
      
!       DO i = 1,3
!         color_val(i) = colors(lev,i)
!       ENDDO
!       
      RETURN
      END SUBROUTINE interp_colors
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
      
      SUBROUTINE read_colormap(cmap_file)
      
      IMPLICIT NONE
      
      CHARACTER(*), INTENT(IN) :: cmap_file
      
      INTEGER :: j,lev
      REAL(rp) :: cmax      
      LOGICAL :: file_exists
      
      INQUIRE(file=TRIM(ADJUSTL(cmap_file)),exist=file_exists)
      IF (file_exists == .FALSE.) THEN
        PRINT*, "colormap file does not exist"
        STOP
      ENDIF
      
      OPEN(UNIT=101,FILE=TRIM(ADJUSTL(cmap_file)))
      READ(101,*) ncolors
      ALLOCATE(colors(ncolors,3))
      DO lev = 1,ncolors
        READ(101,*) (colors(lev,j), j=1,3)
      ENDDO
      CLOSE(101)
      
      cmax = MAXVAL(colors)
      IF (cmax > 1d0+1d-12) THEN
        colors = colors/255d0
      ENDIF      
      
      RETURN
      END SUBROUTINE read_colormap      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      SUBROUTINE setup_cbounds(ne,el_in,el_type,npplt,fig,snap_start,snap_end)
      
      USE evaluate_mod, ONLY: evaluate_solution
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: ne
      INTEGER, DIMENSION(:), INTENT(IN) :: el_in
      INTEGER, DIMENSION(:), INTENT(IN) :: el_type
      INTEGER, DIMENSION(:), INTENT(iN) :: npplt
      TYPE(plot_type), INTENT(INOUT) :: fig
      INTEGER, INTENT(IN) :: snap_start
      INTEGER, INTENT(IN) :: snap_end         
      
      INTEGER :: i,snap,el,pt
      INTEGER :: et,npts
      INTEGER :: start_snap,end_snap
      CHARACTER(:), ALLOCATABLE :: filename

      
      filename = TRIM(ADJUSTL(fig%name))//".cscale"      
      
      ALLOCATE(fig%cscale_vals(snap_end,3)) 
      ALLOCATE(fig%sol_val(mnpp,ne))      
      
      IF (fig%plot_sol_option == 0 .or. fig%type_flag == 1) THEN
        RETURN
      ENDIF
      
      IF (fig%type_flag > 2) THEN   
        OPEN(unit=fig%cscale_unit,file=filename//".out")  
        WRITE(fig%cscale_unit,"(2I5)") snap_start-1,snap_end-1
      ENDIF      
      
      fig%num_cscale_vals = snap_end - snap_start + 1
      
      IF (fig%cscale_option == "file") THEN
        OPEN(unit=fig%cscale_unit,file=filename)
        READ(fig%cscale_unit,*) start_snap,end_snap
        fig%num_cscale_vals = end_snap - start_snap + 1           
        DO snap = snap_start,snap_end
          READ(fig%cscale_unit,*) fig%cscale_vals(snap,1),fig%cscale_vals(snap,2),fig%cscale_vals(snap,3)
        ENDDO
        CLOSE(fig%cscale_unit)
      ENDIF
      
      IF (fig%cscale_option == "spec") THEN
        DO snap = 1,snap_end
          fig%cscale_vals(snap,2) = fig%cscale_min
          fig%cscale_vals(snap,3) = fig%cscale_max
        ENDDO        
      ENDIF
      
      IF (fig%cscale_option == "auto-snap" .or. fig%cscale_option == "auto-all") THEN
      
        DO snap = snap_start,snap_end        
    
          fig%cscale_min = 1d10
          fig%cscale_max = -1d10      
    
    elem: DO el = 1,ne
    
            IF (el_in(el) == 0) THEN
              CYCLE elem
            ENDIF
            
            et = el_type(el)
            i = (et-1)*nord+nord                 
            npts = npplt(i)       
            
            CALL evaluate_solution(el,et,fig%type_flag,snap,fig%sol_val(:,el),npts,phi_hb=phi_hb(:,:,i),phi_sol=phi_sol(:,:,i),ndf_hb=ndof_hb(et),ndf_sol=ndof_sol(et))  
            
            DO pt = 1,npts
              IF (fig%sol_val(pt,el) > fig%cscale_max) THEN
                fig%cscale_max = fig%sol_val(pt,el)
              ENDIF
              
              IF (fig%sol_val(pt,el) < fig%cscale_min) THEN
                fig%cscale_min = fig%sol_val(pt,el)
              ENDIF
            ENDDO            
          
          ENDDO elem
                
          fig%cscale_vals(snap,2) = fig%cscale_min
          fig%cscale_vals(snap,3) = fig%cscale_max
        ENDDO            
      ENDIF
      
      IF (fig%cscale_option == "auto-all") THEN
        fig%cscale_min = 1d10
        fig%cscale_max = -1d10
        DO snap = snap_start,snap_end
          IF (fig%cscale_vals(snap,3) > fig%cscale_max) THEN
            fig%cscale_max = fig%cscale_vals(snap,3)
          ENDIF
          
          IF (fig%cscale_vals(snap,2) < fig%cscale_min) THEN
            fig%cscale_min = fig%cscale_vals(snap,2)
          ENDIF
        ENDDO
        
        DO snap = 1,snap_end
          fig%cscale_vals(snap,2) = fig%cscale_min
          fig%cscale_vals(snap,3) = fig%cscale_max          
        ENDDO
      ENDIF
      
      RETURN
      END SUBROUTINE setup_cbounds
      
            
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   

      SUBROUTINE convert_ps(filename,frmt,density,rm_ps)
      
      IMPLICIT NONE
      
      CHARACTER(*), INTENT(IN) :: filename
      CHARACTER(*), INTENT(IN) :: frmt
      CHARACTER(*), INTENT(IN) :: density
      INTEGER, INTENT(IN) :: rm_ps
      
      CHARACTER(:), ALLOCATABLE :: command
      
      
      IF (TRIM(ADJUSTL(frmt)) == "png" .or. TRIM(ADJUSTL(frmt)) == "jpg" .or. TRIM(ADJUSTL(frmt)) == "tif") THEN
        command = "convert -trim -density "//density//" "//filename//".ps "//filename//"."//frmt
        PRINT("(A)"), "  Converting to "//frmt//" format..."
!         PRINT*, command
        CALL SYSTEM(command)
        IF (rm_ps == 1) THEN
          PRINT("(A)"), "  Removing PostScript file..."
          CALL SYSTEM("rm "//filename//".ps")
        ENDIF        
      ENDIF
      
      
      RETURN
      END SUBROUTINE convert_ps
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       SUBROUTINE make_movie(fig,frmt)
       
       IMPLICIT NONE
       
       TYPE(plot_type), INTENT(IN) :: fig
       CHARACTER(*), INTENT(IN) :: frmt
       
        IF (fig%plot_sol_option == 1 .and. fig%movie_flag == 1) THEN
          CALL SYSTEM("ffmpeg -i "//TRIM(ADJUSTL(fig%name))//"_%04d."//frmt//" -y "//TRIM(ADJUSTL(fig%name))//".mp4")
        ENDIF          
       
       
       RETURN
       END SUBROUTINE make_movie

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      SUBROUTINE plot_ref_el(filename,figure_width,et,ctp,ntri,ect,r,s)
      
      USE basis, ONLY: element_nodes
      USE axes_mod, ONLY: write_xyaxis    
      USE labels_mod, ONLY: write_xyaxis_labels,write_texheader,run_latex, & 
                            write_latex_ps_body,remove_latex_files      
      
      IMPLICIT NONE
      
      CHARACTER(*), INTENT(IN) :: filename
      REAL(rp), INTENT(IN) :: figure_width
      INTEGER, INTENT(IN) :: et
      INTEGER, INTENT(IN) :: ctp
      INTEGER, INTENT(IN) :: ntri
      INTEGER, DIMENSION(:,:), INTENT(IN) :: ect
      REAL(rp), DIMENSION(:), INTENT(IN) :: r
      REAL(rp), DIMENSION(:), INTENT(IN) :: s
      
      INTEGER :: el,nd
      INTEGER :: file_unit
      INTEGER :: nnd
      REAL(rp) :: r_min,r_max,s_min,s_max
      REAL(rp) :: a_x,b_x,a_y,b_y
      REAL(rp), DIMENSION(:), ALLOCATABLE :: rc,sc
      
      file_unit = 90
      
      r_max = 1d0
      r_min = -1d0
      s_max = 1d0
      s_min = -1d0
      
      lr_margin =(612d0 - figure_width)/2d0 
      axes_width = figure_width/1.37d0      
      cscale_width = axes_width*.05d0      
      dash = 0.2d0*cscale_width
      xticklabel_pad = 2d0*dash
      yticklabel_pad = 2d0*dash
      xlabel_pad = 10d0*dash
      ylabel_pad = 15d0*dash        
      cticklabel_pad = 2d0*dash
      clabel_pad = 15d0*dash 
      
      rmin_axes = rmin_page + lr_margin + ylabel_pad
      rmax_axes = rmax_page - lr_margin - 2d0*cscale_width - clabel_pad
      smin_axes = smin_page + lr_margin
      smax_axes = smax_page     
            
      ax = ((rmax_axes-rmin_axes)/(r_max-r_min))
      bx = (rmin_axes*r_max-rmax_axes*r_min)/(r_max-r_min)
      
      ay = ax    
      by = smin_axes-ax*s_min     
      
      smax_axes = ay*s_max + by      

      CALL write_texheader()              
      CALL write_xyaxis_labels("rs")  
      CALL run_latex()         

      CALL write_psheader(TRIM(filename),file_unit)     
      
      DO el = 1,ntri
        DO nd = 1,3
          WRITE(file_unit,"(2(F9.5,1x))") ax*r(ect(nd,el))+bx,ay*s(ect(nd,el))+by                       
        ENDDO
        WRITE(file_unit,"(A)") "draw-tri-element"           
      ENDDO
      
      ALLOCATE(rc((ctp+1)**2),sc((ctp+1)**2))
      CALL element_nodes(et,1,ctp,nnd,rc,sc)          
    
      CALL write_xyaxis(file_unit)  
      
      DO nd = 1,nnd
        WRITE(file_unit,"(3(I5,1x))") 0,0,0       
        WRITE(file_unit,"(2(F9.5,1x))") ax*rc(nd)+bx,ay*sc(nd)+by
        WRITE(file_unit,"(A)") "draw-dot"            
      ENDDO         
      
      CALL write_latex_ps_body(file_unit)      
    
      WRITE(file_unit,"(A)") "showpage"           
      CLOSE(file_unit)   
      
      CALL remove_latex_files()      

      END SUBROUTINE plot_ref_el
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      END MODULE plot_mod