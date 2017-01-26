      PROGRAM plot_stations


      USE plot_globals
      USE globals, ONLY: nsta
      USE grid_file_mod
      USE read_write_output
      USE plot_mod
      USE read_dginp, ONLY: stations_file     

      IMPLICIT NONE
      
      CHARACTER(100) :: out_direc
      CHARACTER(:), ALLOCATABLE :: filename
      INTEGER :: sta,run
      REAL(rp) :: h,u,v
      REAL(rp) :: x_max,x_min,y_max,y_min
      INTEGER :: npt
      REAL(rp), DIMENSION(:), ALLOCATABLE :: xvec,yvec
      INTEGER :: nrun 
      INTEGER :: nsta_plot
      REAL(rp) :: line_width
      

!       nadc = 2
!       ALLOCATE(adc_sol(nadc))
!       adc_sol(1)%line = "/home/sbrus/data-drive/galveston_spline_flux_fix/galveston_tri/adcirc/ESL0/"
!       adc_sol(2)%line = "/home/sbrus/data-drive/galveston_spline_flux_fix/galveston_tri_x64/adcirc/ESL0/"
!       ndg = 1
!       ALLOCATE(dg_sol(ndg))
!       dg_sol(1)%line = "/home/sbrus/data-drive/galveston_spline_flux_fix/galveston_tri/p3/ctp3/hbp3/"

      nadc = 1
      ALLOCATE(adc_sol(nadc))
      adc_sol(1)%line = "/home/sbrus/data-drive/galveston_spline_flux_fix/galveston_tri_x64/adcirc/ESL0/"
      ndg = 2
      ALLOCATE(dg_sol(ndg))
      dg_sol(1)%line = "/home/sbrus/data-drive/galveston_spline_flux_fix/galveston_tri/p3/ctp3/hbp3/"
      dg_sol(2)%line = "/home/sbrus/data-drive/galveston_spline_flux_fix/galveston_tri/p3/ctp1/hbp1/"      
      
      stations_file = "/home/sbrus/data-drive/galveston_spline_flux_fix/grids/stations.d"
      
      nsta_plot = 480            
      
      nsnap_Z = 576
      nsnap_Qx = 576      
      nsnap_Qy = 576      
      
      out_direc = dg_sol(1)%line    
      CALL read_stations()    
      
      zeta%plot_sol_option = 1
      vel%plot_sol_option = 1

      zeta%name = "zeta"
      vel%name = "vel"            
      
      line_width = 1.25d0
      figure_width = 7d0*72d0
      figure_height = 7d0*72d0
      nxtick = 5
      nytick = 5
      nxdec = 0
      nydec = 2
      fontsize = 12
      font = "sans"
      frmt = "png"
      density = "400"
      zeta%rm_ps = 1
      vel%rm_ps = 1
      
      CALL read_colormap("/home/sbrus/Codes/dgswe/plot/work/lines.cmap")
      
      ALLOCATE(zeta%sta_val(nsta,nsnap_Z+1,ndg+nadc),zeta%t_sta(nsnap_Z+1,ndg+nadc))         
      ALLOCATE(vel%sta_val(nsta,nsnap_Z+1,ndg+nadc),vel%t_sta(nsnap_Z+1,ndg+nadc))           
      
      
      nrun = 0            
      
      DO i = 1,nadc
      
        out_direc = adc_sol(i)%line    
        
        nrun = nrun + 1        
      
        IF (zeta%plot_sol_option == 1 .or. vel%plot_sol_option == 1) THEN
          PRINT("(A)"), "Reading zeta solution..."          
          CALL read_fort6163(out_direc,"61",t,eta,nsnap_Z) 
        ENDIF     
        IF (vel%plot_sol_option == 1) THEN
          PRINT("(A)"), "Reading u and v solutions..."       
          CALL read_fort6264(out_direc,"62",t,uu2,vv2,nsnap_Qx)                   
        ELSE 
          ALLOCATE(Qx(1,1,1),Qy(1,1,1))
        ENDIF   
      
        DO sta = 1,nsta
          zeta%sta_val(sta,1,nrun) = 0d0
        ENDDO            
        zeta%t_sta(1,nrun) = 0d0
        DO snap  = 1,nsnap_Z
          DO sta = 1,nsta
            zeta%sta_val(sta,snap+1,nrun) = eta(sta,snap)
          ENDDO
          zeta%t_sta(snap+1,nrun) = t(snap)/86400d0
        ENDDO
      
      
        DO sta = 1,nsta
          vel%sta_val(sta,1,nrun) = 0d0
        ENDDO      
        vel%t_sta(1,nrun) = 0d0
        DO snap = 1,nsnap_Qx
          DO sta = 1,nsta
            vel%sta_val(sta,snap+1,nrun) = sqrt(uu2(sta,snap)**2 + vv2(sta,snap)**2)
          ENDDO
          vel%t_sta(snap+1,nrun) = t(snap)/86400d0
        ENDDO      
  
      ENDDO 
 
      nsnap_Z = nsnap_Z + 1
      nsnap_Qx = nsnap_Qx + 1
      nsnap_Qy = nsnap_Qy + 1 
      DEALLOCATE(t)
 
 
      DO i = 1,ndg
      
        out_direc = dg_sol(i)%line          
      
        nrun = nrun + 1
        
        IF (zeta%plot_sol_option == 1 .or. vel%plot_sol_option == 1) THEN
          PRINT("(A)"), "Reading zeta solution..."          
          CALL read_solution_full(out_direc,"Z.sta","N",t,Z,nsnap_Z) 
        ENDIF
        IF (vel%plot_sol_option == 1) THEN
          PRINT("(A)"), "Reading Qx and Qy solutions..."       
          CALL read_solution_full(out_direc,"Qx.sta","N",t,Qx,nsnap_Qx)        
          CALL read_solution_full(out_direc,"Qy.sta","N",t,Qy,nsnap_Qy) 
        ELSE 
          ALLOCATE(Qx(1,1,1),Qy(1,1,1))
        ENDIF  
        IF (vel%plot_sol_option == 1) THEN
          PRINT("(A)"), "Reading bathymetry solution..."          
          CALL read_solution_full(out_direc,"hb.sta","N",t,hb,nsnap_hb)  
        ENDIF   
 
   
        DO snap  = 1,nsnap_Z
          DO sta = 1,nsta
            zeta%sta_val(sta,snap,nrun) = Z(sta,1,snap)
          ENDDO
          zeta%t_sta(snap,nrun) = t(snap)/86400d0
        ENDDO
      

        DO snap = 1,nsnap_Qx
          DO sta = 1,nsta
            h = Z(sta,1,snap) + hb(sta,1,1)
            u = Qx(sta,1,snap)/h
            v = Qy(sta,1,snap)/h
            vel%sta_val(sta,snap,nrun) = sqrt(u**2 + v**2)
          ENDDO
          vel%t_sta(snap,nrun) = t(snap)/86400d0
        ENDDO
      
      ENDDO
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      x_min = 1
      x_max = nsta_plot
      
      ALLOCATE(xvec(nsta_plot),yvec(nsta_plot))
      
      DO snap = 2,nsnap_Z
        PRINT*, "zeta snap",snap
      
        y_min = MINVAL(zeta%sta_val(x_min:x_max,snap,1:nrun))  
        y_max = MAXVAL(zeta%sta_val(x_min:x_max,snap,1:nrun)) 
        
        y_min = y_min - abs((y_max-ymin)*.05d0)
        y_max = y_max + abs((y_max-ymin)*.05d0)        
        
        WRITE(snap_char,"(I4.4)") snap           
        filename = TRIM(ADJUSTL(zeta%name))//"_xs_"//snap_char          
        
        CALL prep_station_plot(zeta,"zs",snap,x_min,x_max,y_min,y_max,filename)        
        
        DO run = 1,nrun
          DO pt = 1,nsta_plot
            xvec(pt) = pt
            yvec(pt) = zeta%sta_val(pt,snap,run)
          ENDDO
        
          CALL plot_xy(zeta%ps_unit,nsta_plot,xvec,yvec,colors(run+1,:),line_width)
        ENDDO
        
        CALL finish_station_plot(zeta,filename)
      ENDDO
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
      DO snap = 2,nsnap_Qx
        PRINT*, "vel snap",snap      
      
        y_min = MINVAL(vel%sta_val(x_min:x_max,snap,1:nrun))  
        y_max = MAXVAL(vel%sta_val(x_min:x_max,snap,1:nrun)) 
        
        y_min = 0d0
        y_max = y_max + (y_max-ymin)*.05d0        
        
        WRITE(snap_char,"(I4.4)") snap           
        filename = TRIM(ADJUSTL(vel%name))//"_xs_"//snap_char          
        
        CALL prep_station_plot(vel,"vs",snap,x_min,x_max,y_min,y_max,filename)        
        
        DO run = 1,nrun
          DO pt = 1,nsta_plot
            xvec(pt) = pt
            yvec(pt) = vel%sta_val(pt,snap,run)
          ENDDO
        
          CALL plot_xy(vel%ps_unit,nsta_plot,xvec,yvec,colors(run+1,:),line_width)
        ENDDO
        
        CALL finish_station_plot(vel,filename)
      ENDDO      
      
      DEALLOCATE(xvec,yvec)
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
      
      x_min = 0d0
      x_max = 3d0
      
      nxtick = 4
      nytick = 5
      
      snap_end = INT(((nsnap_Z-1)/4)*3)
      
      ALLOCATE(xvec(snap_end+1),yvec(snap_end+1))
      
      DO sta = 1,nsta_plot
        PRINT*, "zeta station",sta      
      
        y_min = MINVAL(zeta%sta_val(sta,1:snap_end,1:nrun))  
        y_max = MAXVAL(zeta%sta_val(sta,1:snap_end,1:nrun)) 
        
        y_min = y_min - abs((y_max-ymin)*.05d0)
        y_max = y_max + abs((y_max-ymin)*.05d0)
        
        WRITE(snap_char,"(I4.4)") sta           
        filename = TRIM(ADJUSTL(zeta%name))//"_ts_"//snap_char          
        
        CALL prep_station_plot(zeta,"zt",sta,x_min,x_max,y_min,y_max,filename)        
        
        DO run = 1,nrun
          DO pt = 1,snap_end
            xvec(pt) = zeta%t_sta(pt,run)
            yvec(pt) = zeta%sta_val(sta,pt,run)
          ENDDO
        
          CALL plot_xy(zeta%ps_unit,snap_end,xvec,yvec,colors(run+1,:),line_width)
        ENDDO
        
        CALL finish_station_plot(zeta,filename)
      ENDDO
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     

      DO sta = 1,nsta_plot
        PRINT*, "vel station",sta
      
        y_min = MINVAL(vel%sta_val(sta,1:snap_end,1:nrun))  
        y_max = MAXVAL(vel%sta_val(sta,1:snap_end,1:nrun)) 
        
        y_min = 0d0
        y_max = y_max + abs((y_max-ymin)*.05d0)        
        
        WRITE(snap_char,"(I4.4)") sta           
        filename = TRIM(ADJUSTL(vel%name))//"_ts_"//snap_char          
        
        CALL prep_station_plot(vel,"vt",sta,x_min,x_max,y_min,y_max,filename)        
        
        DO run = 1,nrun
          DO pt = 1,snap_end
            xvec(pt) = vel%t_sta(pt,run)
            yvec(pt) = vel%sta_val(sta,pt,run)
          ENDDO
        
          CALL plot_xy(vel%ps_unit,snap_end,xvec,yvec,colors(run+1,:),line_width)
        ENDDO
        
        CALL finish_station_plot(vel,filename)
      ENDDO
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   

      END PROGRAM plot_stations