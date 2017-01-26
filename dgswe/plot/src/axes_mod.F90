      MODULE axes_mod

      USE globals, ONLY: rp
      USE plot_globals, ONLY: cscale_width,dash,fontsize,ax,bx,ay,by, &
                              rmin_page,rmax_page,smin_page,smax_page, &
                              rmin_axes,rmax_axes,smin_axes,smax_axes, &
                              rmin_cbar,rmax_cbar,smin_cbar,smax_cbar, &
                              rmin_tbar,rmax_tbar,smin_tbar,smax_tbar, &
                              nxtick,nytick,nctick, &
                              ncolors,colors

      IMPLICIT NONE

      CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      
      SUBROUTINE write_all_axes(file_unit,color_bar,time_bar,t_snap,t_start,t_end)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: file_unit    
      INTEGER, INTENT(IN) :: color_bar
      INTEGER, INTENT(IN) :: time_bar      
      REAL(rp), INTENT(IN) :: t_snap
      REAL(rp), INTENT(IN) :: t_start
      REAL(rp), INTENT(IN) :: t_end 
      

      CALL write_xyaxis(file_unit)   
      
      IF (color_bar == 1) THEN
        CALL write_colorscale(file_unit)
      ENDIF
      IF (time_bar == 1) THEN
        CALL write_tbar(file_unit,t_snap,t_start,t_end)
      ENDIF            
      
      RETURN
      END SUBROUTINE write_all_axes      


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE write_colorscale(file_unit)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: file_unit

      
      INTEGER :: lev
      INTEGER :: tick         
      REAL(rp) :: r0,r1,s0,s1
      REAL(rp) :: ds
      REAL(rp) :: cval
      REAL(rp) :: dc   
      CHARACTER(20) :: cchar
      

      
      r0 = rmax_axes
      s0 = smin_axes
      r1 = rmax_axes            

      
      ds = (smax_cbar-smin_cbar)/(ncolors-1)
      DO lev = 1,ncolors-1
      
        s1 = s0 + ds
      
        WRITE(file_unit,"(A,3(F9.5,1x),A)") "[",colors(lev+1,1),colors(lev+1,2),colors(lev+1,3),"]"       
        WRITE(file_unit,"(A,3(F9.5,1x),A)") "[",colors(lev,1),colors(lev,2),colors(lev,3),"]"  
        WRITE(file_unit,"(A,4(F9.5,1x),A)") "[",r0,s0,r1,s1,"]" 
        WRITE(file_unit,"(A,4(F9.5,1x),A)") "[",rmin_cbar,s0,rmax_cbar,s1,"]"         
        WRITE(file_unit,"(A)") "recfill" 
        
        s0 = s1
        
      ENDDO
                     
      WRITE(file_unit,"(2(F9.5,1x))") rmin_cbar,s1
      WRITE(file_unit,"(2(F9.5,1x))") rmax_cbar,s1
      WRITE(file_unit,"(2(F9.5,1x))") rmax_cbar,smin_axes
      WRITE(file_unit,"(2(F9.5,1x))") rmin_cbar,smin_axes      
      WRITE(file_unit,"(A)") "draw-box"       
      

      
      s0 = smin_axes 
      ds = (smax_cbar-smin_cbar)/(nctick-1)      
      DO tick = 1,nctick
              
        WRITE(file_unit,"(2(F9.5,1x))") rmin_cbar,s0 
        WRITE(file_unit,"(2(F9.5,1x))") rmin_cbar+dash,s0     
        WRITE(file_unit,"(A)") "draw-line" 
        
        WRITE(file_unit,"(2(F9.5,1x))") rmax_cbar,s0
        WRITE(file_unit,"(2(F9.5,1x))") rmax_cbar-dash,s0     
        WRITE(file_unit,"(A)") "draw-line"  

        s0 = s0 + ds
      ENDDO       
           
      
      RETURN
      END SUBROUTINE write_colorscale
      


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE write_xyaxis(file_unit)
      
      IMPLICIT NONE
            
      INTEGER, INTENT(IN) :: file_unit
      
      INTEGER :: i
      INTEGER :: expnt
      REAL(rp) :: dr,ds
      REAL(rp) :: r0,r1,s0,s1
      REAL(rp) :: xval,yval
    
      CHARACTER(20) :: xchar,ychar
      
      WRITE(file_unit,"(A)")  "newpath"     
      WRITE(file_unit,"(2(F9.5,1x),A)") rmin_page,smin_page," moveto"
      WRITE(file_unit,"(2(F9.5,1x),A)") rmax_page,smin_page," lineto"
      WRITE(file_unit,"(2(F9.5,1x),A)") rmax_page,smax_page," lineto"  
      WRITE(file_unit,"(2(F9.5,1x),A)") rmin_page,smax_page," lineto"    
      WRITE(file_unit,"(2(F9.5,1x),A)") rmin_page,smin_page," lineto"        
      WRITE(file_unit,"(2(F9.5,1x),A)") rmin_axes,smin_axes," moveto"   
      WRITE(file_unit,"(2(F9.5,1x),A)") rmin_axes,smax_axes," lineto"  
      WRITE(file_unit,"(2(F9.5,1x),A)") rmax_axes,smax_axes," lineto"
      WRITE(file_unit,"(2(F9.5,1x),A)") rmax_axes,smin_axes," lineto"
      WRITE(file_unit,"(2(F9.5,1x),A)") rmin_axes,smin_axes," lineto"     
      WRITE(file_unit,"(A)") "gsave 1 1 1 setrgbcolor fill grestore"      
      WRITE(file_unit,"(A)") "clear"      
      
      
      WRITE(file_unit,"(2(F9.5,1x))") rmin_axes,smin_axes
      WRITE(file_unit,"(2(F9.5,1x))") rmax_axes,smin_axes
      WRITE(file_unit,"(2(F9.5,1x))") rmax_axes,smax_axes
      WRITE(file_unit,"(2(F9.5,1x))") rmin_axes,smax_axes      
      WRITE(file_unit,"(A)") "draw-box"  
      
!       ! x-axis line
!       WRITE(file_unit,"(2(F9.5,1x))") rmin_axes,smin_axes
!       WRITE(file_unit,"(2(F9.5,1x))") rmax_axes,smin_axes     
!       WRITE(file_unit,"(A)") "draw-line"  
!       
!       ! y-axis line
!       WRITE(file_unit,"(2(F9.5,1x))") rmin_axes,smin_axes
!       WRITE(file_unit,"(2(F9.5,1x))") rmin_axes,smax_axes     
!       WRITE(file_unit,"(A)") "draw-line"        

      
      dr = (rmax_axes-rmin_axes)/(real(nxtick,rp)-1d0)
      

      r0 = rmin_axes
      
      xval = (r0-bx)/ax
      expnt = INT(LOG10(xval))
      IF (expnt <= 3) THEN
        expnt = 0
      ENDIF
      
      DO i = 1,nxtick        
        WRITE(file_unit,"(2(F9.5,1x))") r0,smin_axes+dash
        WRITE(file_unit,"(2(F9.5,1x))") r0,smin_axes
        WRITE(file_unit,"(A)") "draw-line" 
        
        WRITE(file_unit,"(2(F9.5,1x))") r0,smax_axes
        WRITE(file_unit,"(2(F9.5,1x))") r0,smax_axes-dash
        WRITE(file_unit,"(A)") "draw-line"         

        r0 = r0 + dr
      ENDDO
      
      
      
      
      
      IF (nytick >= 10000) THEN
        ds = dr   
      ELSE
        ds = (smax_axes-smin_axes)/(real(nytick,rp)-1d0)      
      ENDIF      

      s0 = smin_axes      
      
      yval = (s0-by)/ay
      expnt = INT(LOG10(yval))
      IF (expnt <= 3) THEN
        expnt = 0
      ENDIF      
      
      DO i = 1,nytick
      
        IF (s0 > smax_axes) THEN
          EXIT        
        ENDIF
        
        WRITE(file_unit,"(2(F9.5,1x))") rmin_axes,s0
        WRITE(file_unit,"(2(F9.5,1x))") rmin_axes+dash,s0
        WRITE(file_unit,"(A)") "draw-line" 
        
        WRITE(file_unit,"(2(F9.5,1x))") rmax_axes,s0
        WRITE(file_unit,"(2(F9.5,1x))") rmax_axes-dash,s0
        WRITE(file_unit,"(A)") "draw-line"         
 
        s0 = s0 + ds        
      
      ENDDO 
      
      
      RETURN
      END SUBROUTINE write_xyaxis


      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE write_tbar(file_unit,t_snap,t_start,t_end)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: file_unit      
      REAL(rp), INTENT(IN) :: t_snap
      REAL(rp), INTENT(IN) :: t_start
      REAL(rp), INTENT(IN) :: t_end
      
      REAL(rp) ::  rmax_snap
      
      WRITE(file_unit,"(2(F9.5,1x))") rmin_tbar,smin_tbar
      WRITE(file_unit,"(2(F9.5,1x))") rmax_tbar,smin_tbar
      WRITE(file_unit,"(2(F9.5,1x))") rmax_tbar,smax_tbar
      WRITE(file_unit,"(2(F9.5,1x))") rmin_tbar,smax_tbar      
      WRITE(file_unit,"(A)") "draw-box"  
      
      rmax_snap = (t_snap-t_end)/(t_start-t_end)*rmin_tbar + (t_snap-t_start)/(t_end-t_start)*rmax_tbar
      
      WRITE(file_unit,"(A)")  "newpath"     
      WRITE(file_unit,"(2(F9.5,1x),A)") rmin_tbar,smin_tbar," moveto"
      WRITE(file_unit,"(2(F9.5,1x),A)") rmax_snap,smin_tbar," lineto"
      WRITE(file_unit,"(2(F9.5,1x),A)") rmax_snap,smax_tbar," lineto"  
      WRITE(file_unit,"(2(F9.5,1x),A)") rmin_tbar,smax_tbar," lineto"    
      WRITE(file_unit,"(A)") "closepath"            
      WRITE(file_unit,"(A)") "gsave 0 0 0 setrgbcolor fill grestore"        
      
      RETURN
      END SUBROUTINE write_tbar                  
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       


      END MODULE axes_mod
      