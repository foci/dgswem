      MODULE labels_mod

      USE globals, ONLY: rp,r_earth,deg2rad
      USE plot_globals, ONLY: cscale_width,fontsize,font,ax,bx,ay,by, &
                              rmin_page,rmax_page,smin_page,smax_page, &
                              rmin_axes,rmax_axes,smin_axes,smax_axes, &
                              rmin_cbar,rmax_cbar,smin_cbar,smax_cbar, &
                              rmin_tbar,rmax_tbar,smin_tbar,smax_tbar, &
                              xticklabel_pad,yticklabel_pad,cticklabel_pad, &
                              xlabel_pad,ylabel_pad,clabel_pad, &
                              nxtick,nytick,nctick, &
                              nxdec,nydec,ncdec,ntdec, &
                              dr_xlabel,ds_ylabel,ds_clabel, &
                              plot_type,char_array
      USE read_dginp, ONLY: sphi0,slam0                              

      IMPLICIT NONE
      
      INTEGER :: tex_unit = 10   
      INTEGER :: nline_texfile      
      INTEGER :: tex_output_unit = 11            

      CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      SUBROUTINE write_texheader()
      
      IMPLICIT NONE
      

      OPEN(UNIT=tex_unit,FILE="labels.tex")           
      WRITE(tex_unit,"(A,I2,A)") "\documentclass[",fontsize,"pt]{article}"
      WRITE(tex_unit,"(A)") "\pagestyle{empty}"
!       WRITE(tex_unit,"(A)") "\usepackage[absolute,showboxes]{textpos}"
      WRITE(tex_unit,"(A)") "\usepackage[absolute]{textpos}"       
      WRITE(tex_unit,"(A)") "\usepackage{graphicx}"  
      WRITE(tex_unit,"(A)") "\usepackage{xcolor}"            
      WRITE(tex_unit,"(A)") "\setlength{\TPHorizModule}{1pt}"  
      WRITE(tex_unit,"(A)") "\setlength{\TPVertModule}{1pt}"
      WRITE(tex_unit,"(A)") "\textblockorigin{0pt}{0pt}"
      WRITE(tex_unit,"(A)") "\setlength{\parindent}{0pt}"
      IF (TRIM(ADJUSTL(font)) == "times") THEN
        WRITE(tex_unit,"(A)") "\usepackage{mathptmx}"       
      ELSE IF (TRIM(ADJUSTL(font)) == "sans") THEN
        WRITE(tex_unit,"(A)") "\renewcommand{\rmdefault}{\sfdefault}"     
        WRITE(tex_unit,"(A)") "\usepackage{sansmath}" 
        WRITE(tex_unit,"(A)") "\sansmath"        
      ENDIF
      WRITE(tex_unit,"(A)") "\begin{document}"
       
      
      RETURN
      END SUBROUTINE write_texheader
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      SUBROUTINE close_tex()
      
      IMPLICIT NONE 
      
      WRITE(tex_unit,"(A)") "\end{document}"
      CLOSE(tex_unit)
      
      RETURN
      END SUBROUTINE close_tex  
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      SUBROUTINE run_latex()
      
      IMPLICIT NONE
      
      CALL close_tex()
      
!       CALL SYSTEM("latex labels.tex")
!       CALL SYSTEM("dvips -o labels.ps labels.dvi")
      
      CALL SYSTEM("latex labels.tex >/dev/null")
      CALL SYSTEM("dvips -q -o labels.ps labels.dvi")        
      
      RETURN
      END SUBROUTINE run_latex      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      SUBROUTINE latex_axes_labels(fig,t_snap,t_start,t_end)
      
      IMPLICIT NONE
      
      TYPE(plot_type), INTENT(IN) :: fig
      REAL(rp), INTENT(IN) :: t_snap
      REAL(rp), INTENT(IN) :: t_start
      REAL(rp), INTENT(IN) :: t_end 
      
                  
      
      CALL write_texheader()        
      
      IF (abs(sphi0) > 0d0 .and. abs(slam0) > 0d0) THEN
        CALL write_xyaxis_labels("ll")        
      ELSE
        CALL write_xyaxis_labels("xy")  
!         CALL write_xyaxis_labels("rs")          
      ENDIF      
      
      IF (fig%cbar_flag == 1) THEN      
        CALL write_caxis_labels(fig%tbar_flag,fig%sol_min,fig%sol_max,fig%sol_label)
      ENDIF
      IF (fig%tbar_flag == 1) THEN
        CALL write_tbar_labels(t_snap)
      ENDIF                
          
      
      RETURN
      END SUBROUTINE latex_axes_labels        
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE write_caxis_labels(time_bar,sol_min,sol_max,sol_label)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: time_bar
      REAL(rp), INTENT(IN) :: sol_min
      REAL(rp), INTENT(IN) :: sol_max
      CHARACTER(*), INTENT(IN) :: sol_label

      
      INTEGER :: lev
      INTEGER :: tick         
      REAL(rp) :: r0,r1,s0,s1
      REAL(rp) :: ds
      REAL(rp) :: cval
      REAL(rp) :: dc   
      CHARACTER(20) :: cchar
      
      rmin_cbar = rmax_axes + cscale_width
      rmax_cbar = rmin_cbar + cscale_width
      smin_cbar = smin_axes
      IF (time_bar == 1) THEN      
!         smax_cbar = .75d0*smax_axes     
        smax_cbar = smin_tbar - 2d0*cscale_width            
      ELSE 
        smax_cbar = smax_axes 
      ENDIF     
            

      r0 = rmax_axes
      s0 = smin_axes
      r1 = rmax_axes            
      
      cval = sol_min
      dc = (sol_max-sol_min)/(nctick-1)
      ds = (smax_cbar-smin_cbar)/(nctick-1)      
      DO tick = 1,nctick               
        
        CALL format_number(ncdec,cval,0,cchar)         
!         WRITE(tex_unit,"(A,2(1x,F9.5))")  "("//TRIM(ADJUSTL(cchar))//")",rmax_cbar+cticklabel_pad,s0       
!         WRITE(tex_unit,"(A)") "caxis-tick-labels"
        
        WRITE(tex_unit,"(A,F9.5,A,F9.5,A)") "\begin{textblock}{400}[0,0.5](",rmax_cbar+cticklabel_pad,",",smax_page-s0,")"
        WRITE(tex_unit,"(A)") TRIM(ADJUSTL(cchar))        
        WRITE(tex_unit,"(A)") "\end{textblock}"                 
        
        cval = cval + dc
        s0 = s0 + ds
      ENDDO       
      
!         WRITE(tex_unit,"(A,2(1x,F9.5))")  "("//TRIM(ADJUSTL(sol_label))//")",rmax_cbar+clabel_pad,(smin_axes+smax_axes)/2d0       
!         WRITE(tex_unit,"(A)") "caxis-labels"     
        
      WRITE(tex_unit,"(A,F9.5,A,F9.5,A)") "\begin{textblock}{50}[0,0.5](",-50+rmax_cbar+clabel_pad,",",smax_page-(smax_cbar+smin_cbar)/2d0,")"
      WRITE(tex_unit,"(A)") "\hfill \rotatebox[origin=c]{90}{"//TRIM(ADJUSTL(sol_label))//"}\vskip-\TPboxrulesize"       
      WRITE(tex_unit,"(A)") "\end{textblock}"         
      
      RETURN
      END SUBROUTINE write_caxis_labels 
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE write_xyaxis_labels(axis_label_option)
      
      IMPLICIT NONE

      CHARACTER(2), INTENT(IN) :: axis_label_option
      
      INTEGER :: i
      INTEGER :: expnt,expnt_min,expnt_max
      REAL(rp) :: r0,r1,s0,s1
      REAL(rp) :: xval,yval
    
      CHARACTER(20) :: xchar,ychar     
      CHARACTER(25) :: xlabel,ylabel
      
      IF (axis_label_option == "xy") THEN
        xlabel = "$x$"
        ylabel = "$y$"
      ELSE IF (axis_label_option == "rs") THEN
        xlabel = "$r$"
        ylabel = "$s$"
      ELSE IF (axis_label_option == "ll") THEN
        xlabel = "longitude"
        ylabel = "latitude"
      ELSE IF (axis_label_option == "zt") THEN
        xlabel = "time (days)"
        ylabel = "surface elevation (m)"      
      ELSE IF (axis_label_option == "vt") THEN
        xlabel = "time (days)"
        ylabel = "velocity (m/s)"      
      ELSE IF (axis_label_option == "zs") THEN
        xlabel = "station"
        ylabel = "surface elevation (m)"      
      ELSE IF (axis_label_option == "vs") THEN
        xlabel = "station"
        ylabel = "velocity (m/s)"      
      ENDIF
      

      dr_xlabel = (rmax_axes-rmin_axes)/(real(nxtick,rp)-1d0)
      
      
      xval = (rmin_axes-bx)/ax
      expnt_min = INT(LOG10(abs(xval)))
      IF ((expnt_min <= 3 .AND. expnt_min > -2) .OR. xval < 1d-10) THEN
        expnt_min = 0
      ENDIF
      
      xval = (rmax_axes-bx)/ax      
      expnt_max = INT(LOG10(abs(xval)))
      IF ((expnt_max <= 3 .AND. expnt_max > -2) .OR. xval < 1d-10) THEN
        expnt_max = 0
      ENDIF      
      
      IF (expnt_min /= 0) THEN
        expnt = expnt_min
      ELSEIF (expnt_max /= 0) THEN
        expnt = expnt_max
      ENDIF

      
      r0 = rmin_axes
      
      DO i = 1,nxtick             
        
        xval = (r0-bx)/ax
        IF (axis_label_option == "ll") THEN
          xval = xval/(r_earth*cos(sphi0)) + slam0
          xval = xval/deg2rad
          expnt = 0          
        ENDIF
        CALL format_number(nxdec,xval,expnt,xchar)        
      
                  
        WRITE(tex_unit,"(A,F9.5,A,F9.5,A)") "\begin{textblock}{400}[0.5,0](",r0,",",smax_page-smin_axes+xticklabel_pad,")"
        WRITE(tex_unit,"(A)") "\centerline{"//TRIM(ADJUSTL(xchar))//"}"
        WRITE(tex_unit,"(A)") "\end{textblock}"        
        
!         WRITE(tex_unit,"(A,2(1x,F9.5))")  "("//TRIM(ADJUSTL(xchar))//")",r0,smin_axes-xticklabel_pad    
!         WRITE(tex_unit,"(A)") "xaxis-labels"           
        
        r0 = r0 + dr_xlabel
      ENDDO
      
      
      WRITE(tex_unit,"(A,F9.5,A,F9.5,A)") "\begin{textblock}{400}[0.5,0](",(rmax_axes+rmin_axes)/2d0,",",smax_page-smin_axes+xlabel_pad,")"
      WRITE(tex_unit,"(A)") "\centerline{"//TRIM(xlabel)//"}"        
      WRITE(tex_unit,"(A)") "\end{textblock}"  
      
!       WRITE(tex_unit,"(A)") TRIM(ADJUSTL(math_font))//" choosefont"         
!       WRITE(tex_unit,"(A,2(1x,F9.5))")  "(x)",(rmax_axes+rmin_axes)/2d0,smin_axes-xlabel_pad    
!       WRITE(tex_unit,"(A)") "xaxis-labels"      
!       WRITE(tex_unit,"(A)") TRIM(ADJUSTL(main_font))//" choosefont"        
 
      
      IF (expnt /= 0) THEN      
        WRITE(tex_unit,"(A,F9.5,A,F9.5,A)") "\begin{textblock}{400}[0,0](",rmin_axes+(nxtick-1)*dr_xlabel,",",smax_page-smin_axes+xticklabel_pad+xlabel_pad/2d0,")"
        WRITE(tex_unit,"(A,I2,A)") "$\times10^{",expnt,"}$"        
        WRITE(tex_unit,"(A)") "\end{textblock}"         
      ENDIF
      
      
      
      
      
      
      IF (nytick >= 10000) THEN
        ds_ylabel = dr_xlabel      
      ELSE
        ds_ylabel = (smax_axes-smin_axes)/(real(nytick,rp)-1d0)      
      ENDIF
      
      
      
      
      yval = (smin_axes-by)/ay
      expnt_min = INT(LOG10(abs(yval)))
      IF ((expnt_min <= 3 .AND. expnt_min > -2) .OR. yval < 1d-10) THEN
        expnt_min = 0
      ENDIF
      
      yval = (smax_axes-by)/ay      
      expnt_max = INT(LOG10(abs(yval)))
      IF ((expnt_max <= 3 .AND. expnt_max > -2) .OR. yval < 1d-10) THEN
        expnt_max = 0
      ENDIF      
      
      IF (expnt_min /= 0) THEN
        expnt = expnt_min
      ELSEIF (expnt_max /= 0) THEN
        expnt = expnt_max
      ENDIF      
          
      
      s0 = smin_axes         
      
      DO i = 1,nytick
      
        IF (s0 > smax_axes+1d-12) THEN
          EXIT        
        ENDIF
        
        yval = (s0-by)/ay
        IF (axis_label_option == "ll") THEN        
          yval = yval/r_earth
          yval = yval/deg2rad
          expnt = 0
        ENDIF
        
        CALL format_number(nydec,yval,expnt,ychar) 
        
        WRITE(tex_unit,"(A,F9.5,A,F9.5,A)") "\begin{textblock}{400}[1,0.5](",rmin_axes-yticklabel_pad,",",smax_page-s0,")"
        WRITE(tex_unit,"(A)") "\hfill "//TRIM(ADJUSTL(ychar))        
        WRITE(tex_unit,"(A)") "\end{textblock}"  
        
!         WRITE(tex_unit,"(A,2(1x,F9.5))")  "("//TRIM(ADJUSTL(ychar))//")",rmin_axes-yticklabel_pad,s0        
!         WRITE(tex_unit,"(A)") "yaxis-tick-labels"       

        s0 = s0 + ds_ylabel        
      
      ENDDO         
      
      WRITE(tex_unit,"(A,F9.5,A,F9.5,A)") "\begin{textblock}{50}[0,0.5](",rmin_axes-ylabel_pad,",",smax_page-(smax_axes+smin_axes)/2d0,")"
      WRITE(tex_unit,"(A)") "\rotatebox[origin=c]{90}{"//TRIM(ylabel)//"}\vskip-\TPboxrulesize"        
      WRITE(tex_unit,"(A)") "\end{textblock}"  
      
      
!       WRITE(tex_unit,"(A)") TRIM(ADJUSTL(math_font))//" choosefont"        
!       WRITE(tex_unit,"(A,2(1x,F9.5))")  "(y)",rmin_axes-ylabel_pad,(smax_axes+smin_axes)/2d0    
!       WRITE(tex_unit,"(A)") "yaxis-labels"  
!       WRITE(tex_unit,"(A)") TRIM(ADJUSTL(main_font))//" choosefont"        
      
      IF (expnt /= 0) THEN      
        WRITE(tex_unit,"(A,F9.5,A,F9.5,A)") "\begin{textblock}{400}[0,1](",rmin_axes,",",smax_page-(smax_axes + yticklabel_pad),")"
        WRITE(tex_unit,"(A,I2,A)") "$\times10^{",expnt,"}$"        
        WRITE(tex_unit,"(A)") "\end{textblock}"         
      ENDIF
      
      RETURN
      END SUBROUTINE write_xyaxis_labels
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE write_tbar_labels(t_snap)
      
      IMPLICIT NONE

      REAL(rp), INTENT(IN) :: t_snap
      
      REAL(rp) :: tday
      CHARACTER(20) :: tchar
      
      tday = t_snap/86400d0
      
      CALL format_number(ntdec,tday,0,tchar)       
      
      WRITE(tex_unit,"(A,F9.5,A,F9.5,A)") "\begin{textblock}{400}[0.5,0](",(rmin_tbar+rmax_tbar)/2d0,",",smax_page-(smax_axes),")"
      WRITE(tex_unit,"(A,F9.2,A)") "\centerline{$t="//TRIM(ADJUSTL(tchar))//"$ days}"        
      WRITE(tex_unit,"(A)") "\end{textblock}"      
      
      RETURN
      END SUBROUTINE write_tbar_labels    
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE latex_element_labels(ne,label_option,el_type,el_in,nverts,xy,ect,nnfbed,nfbedn,ged2el)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: ne
      CHARACTER(*), INTENT(IN) :: label_option 
      INTEGER, DIMENSION(:), INTENT(IN) :: el_type
      INTEGER, DIMENSION(:), INTENT(IN) :: el_in      
      INTEGER, DIMENSION(:), INTENT(IN) :: nverts
      REAL(rp), DIMENSION(:,:), INTENT(IN) :: xy
      INTEGER, DIMENSION(:,:), INTENT(IN) :: ect
      INTEGER, INTENT(IN) :: nnfbed
      INTEGER, DIMENSION(:), INTENT(IN) :: nfbedn
      INTEGER, DIMENSION(:,:), INTENT(IN) :: ged2el
      
      INTEGER :: i
      INTEGER :: el,nd
      INTEGER :: et,nv 
      INTEGER :: ed,ged
      INTEGER :: label_unit = 20
      INTEGER :: ne_list
      INTEGER, DIMENSION(:), ALLOCATABLE :: el_list
      INTEGER, DIMENSION(:), ALLOCATABLE :: el_flag
      REAL(rp) :: xsum,ysum
      REAL(rp) :: rpos,spos      
      CHARACTER(20) :: elchar
      LOGICAL :: file_exists
      
      IF (TRIM(ADJUSTL(label_option)) == "off") THEN
      
        RETURN
        
      ELSE IF (TRIM(ADJUSTL(label_option)) == "file") THEN
      
        OPEN(unit=label_unit,file="element.label")
        READ(label_unit,*) ne_list 
        ALLOCATE(el_list(ne_list))
        DO i = 1,ne_list
          READ(label_unit,*) el_list(i)
        ENDDO          
        
      ELSE IF (TRIM(ADJUSTL(label_option)) == "all") THEN
      
        ne_list = ne
        ALLOCATE(el_list(ne_list))
        DO i = 1,ne_list
          el_list(i) = i
        ENDDO
        
      ELSE IF (TRIM(ADJUSTL(label_option)) == "bou") THEN
      
        ALLOCATE(el_list(ne))
        ALLOCATE(el_flag(ne))
        el_flag = 0
        ne_list = 0
        DO ed = 1,nnfbed
          ged = nfbedn(ed)
          el = ged2el(1,ged)
          
          IF (el_flag(el) == 0) THEN
            ne_list = ne_list + 1
            el_list(ne_list) = el
            el_flag(el) = 1
          ENDIF
        ENDDO
        
      ELSE
      
        PRINT("(A)"), "Error: Mesh plot element label option not recognized"
        STOP
        
      ENDIF    

      IF (ne_list > 20000) THEN 
        PRINT("(A)"), "Warning: number of element labels is too large" 
        PRINT("(A)"), "         skipping to prevent exceeding latex memory"
        RETURN
      ENDIF

 elem:DO i = 1,ne_list
        el = el_list(i)
        
        IF (el_in(el) == 0) THEN
          CYCLE elem
        ENDIF
        
        et = el_type(el)
        nv = nverts(et)
        
        xsum = 0d0
        ysum = 0d0
        DO nd = 1,nv
          xsum = xsum + ax*xy(1,ect(nd,el)) + bx
          ysum = ysum + ay*xy(2,ect(nd,el)) + by
        ENDDO
        
        rpos = xsum/real(nv,rp)
        spos = ysum/real(nv,rp)
        
        
        WRITE(elchar,"(I20)") el     
        WRITE(tex_unit,"(A,F9.5,A,F9.5,A)") "\begin{textblock}{400}[0,0](",rpos,",",smax_page-spos,")"
        WRITE(tex_unit,"(A)") "{\color{red} \tiny "//TRIM(ADJUSTL(elchar))//"}"
        WRITE(tex_unit,"(A)") "\end{textblock}"        
      ENDDO elem
      
      
      RETURN
      END SUBROUTINE latex_element_labels      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   

      SUBROUTINE latex_node_labels(nn,label_option,xy,nbou,fbseg,fbnds)
      
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nn
      CHARACTER(*), INTENT(IN) :: label_option 
      REAL(rp), DIMENSION(:,:), INTENT(IN) :: xy
      INTEGER, INTENT(IN) :: nbou
      INTEGER, DIMENSION(:,:), INTENT(IN) :: fbseg
      INTEGER, DIMENSION(:,:), INTENT(IN) :: fbnds
      
      INTEGER :: i
      INTEGER :: nd,bou
      INTEGER :: label_unit = 21
      INTEGER :: nn_list
      INTEGER, DIMENSION(:), ALLOCATABLE :: nd_list
      REAL(rp) rpt,spt
      CHARACTER(20) :: ndchar   
      LOGICAL :: file_exists
      
      IF (TRIM(ADJUSTL(label_option)) == "off") THEN
      
        RETURN
        
      ELSE IF (TRIM(ADJUSTL(label_option)) == "file") THEN

        OPEN(unit=label_unit,file="node.label")
        READ(label_unit,*) nn_list 
        ALLOCATE(nd_list(nn_list))
        DO i = 1,nn_list
          READ(label_unit,*) nd_list(i)
        ENDDO          
        
      ELSE IF (TRIM(ADJUSTL(label_option)) == "all") THEN
      
        nn_list = nn
        ALLOCATE(nd_list(nn_list))
        DO i = 1,nn_list
          nd_list(i) = i
        ENDDO
        
      ELSE IF (TRIM(ADJUSTL(label_option)) == "bou") THEN
      
        ALLOCATE(nd_list(nn))
        nn_list = 0
        DO bou = 1,nbou
          DO nd = 1,fbseg(1,bou)          
            nn_list = nn_list + 1
            nd_list(nn_list) = fbnds(nd,bou)          
          ENDDO
        ENDDO
        
      ELSE
      
        PRINT("(A)"), "Error: Mesh plot node label option not recognized"
        STOP
        
      ENDIF
      
      
      IF (nn_list > 20000) THEN 
        PRINT("(A)"), "Warning: number of node labels is too large" 
        PRINT("(A)"), "         skipping to prevent exceeding latex memory"
        RETURN
      ENDIF      
      
      
      DO i = 1,nn_list
        nd = nd_list(i)

        rpt = ax*xy(1,nd)+bx
        spt = ay*xy(2,nd)+by
        
        IF ((rpt > rmin_axes .and. rpt < rmax_axes) .and. &
            (spt > smin_axes .and. spt < smax_axes) ) THEN
            
          WRITE(ndchar,"(I20)") nd     
          WRITE(tex_unit,"(A,F9.5,A,F9.5,A)") "\begin{textblock}{400}[0,0](",rpt,",",smax_page-spt,")"
          WRITE(tex_unit,"(A)") "{\color{blue} \tiny "//TRIM(ADJUSTL(ndchar))//"}"
          WRITE(tex_unit,"(A)") "\end{textblock}"              
            
        ENDIF

      ENDDO 
      
      
      RETURN
      END SUBROUTINE latex_node_labels      
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      SUBROUTINE format_number(ndec,val,expnt,val_char)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: ndec
      REAL(rp), INTENT(IN) :: val
      INTEGER, INTENT(IN) :: expnt
      CHARACTER(*), INTENT(OUT) :: val_char
      
      CHARACTER(:), ALLOCATABLE :: frmt   
      CHARACTER(1) :: ndec_char
      
      IF (ndec > 0) THEN             
        WRITE(ndec_char,"(I1)") ndec      
        frmt = "(F20."//ndec_char//")"      
      ELSE
        frmt = "(I20)"
      ENDIF      

      IF (ndec > 0) THEN
        WRITE(val_char,frmt) val/(10d0**expnt)      
      ELSE
        WRITE(val_char,frmt) NINT(val/(10d0**expnt))
      ENDIF      
      
      RETURN
      END SUBROUTINE format_number      
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE read_latex(header_lines,nhead,body_lines,nbody)
      
      IMPLICIT NONE
            
      TYPE(char_array), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: header_lines
      INTEGER, INTENT(OUT) :: nhead
      TYPE(char_array), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: body_lines      
      INTEGER, INTENT(OUT) :: nbody
      
      INTEGER :: i
      CHARACTER(1000) :: line
      INTEGER :: read_stat       
      
      OPEN(UNIT=tex_output_unit, FILE="labels.ps")
      READ(tex_output_unit,"(A)",IOSTAT=read_stat) line    ! discard first line  
      nhead = 0
head: DO
        READ(tex_output_unit,"(A)",IOSTAT=read_stat) line    
        nhead = nhead + 1         
        IF (read_stat < 0 ) THEN
          PRINT*, "Error in LaTeX file"
          STOP
        ELSEIF (TRIM(ADJUSTL(line)) == "%%EndSetup") THEN
          EXIT head   
        ENDIF
                
      ENDDO head  
      
      nbody = 0
body: DO
        READ(tex_output_unit,"(A)",IOSTAT=read_stat) line       
        IF (read_stat < 0 ) THEN
          EXIT body 
        ELSE 
          nbody = nbody + 1           
        ENDIF                
      ENDDO body              
 
      CLOSE(tex_output_unit)       
      
      
      
      ALLOCATE(header_lines(nhead))
      ALLOCATE(body_lines(nbody))
      OPEN(UNIT=tex_output_unit, FILE="labels.ps",POSITION="REWIND")
      READ(tex_output_unit,"(A)") line                    ! discard first line
      DO i = 1,nhead       
        READ(tex_output_unit,"(A)") header_lines(i)%line      
!         PRINT "(A200)", header_lines(i)%line
      ENDDO
!       PRINT*, "------------------------------------------"
      DO i = 1,nbody
        READ(tex_output_unit,"(A)") body_lines(i)%line   
!         PRINT "(A)", body_lines(i)%line       
      ENDDO
      CLOSE(tex_output_unit)   
      
      RETURN
      END SUBROUTINE
            
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE write_latex_ps_header(file_unit)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: file_unit
      
      CHARACTER(1000) :: line
      INTEGER :: read_stat        
      
      OPEN(UNIT=tex_output_unit, FILE="labels.ps")
      WRITE(file_unit,"(A)") "%!PS-Adobe-3.0"      
!       WRITE(file_unit,"(A)") "%!PS-Adobe-3.0 EPSF-3.0"
!       WRITE(file_unit,"(A,4(F9.5,1x))") "%%BoundingBox: ",rmin_axes,smin_axes,rmax_axes,smax_axes           
      nline_texfile = 0
head: DO
        READ(tex_output_unit,"(A)",IOSTAT=read_stat) line
        nline_texfile = nline_texfile + 1        
        IF (read_stat < 0 ) THEN
          PRINT*, "Error in LaTeX file"
          STOP
        ELSEIF (TRIM(ADJUSTL(line)) == "%%EndSetup") THEN
          WRITE(file_unit,"(A)") line
          EXIT head
        ELSEIF (nline_texfile > 1) THEN
          WRITE(file_unit,"(A)") line                        
        ENDIF
                
      ENDDO head          
      
      
      RETURN
      END SUBROUTINE write_latex_ps_header

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  


      SUBROUTINE write_latex_ps_body(file_unit)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: file_unit
      
      CHARACTER(1000) :: line
      INTEGER :: read_stat   
      
tail: DO
        READ(tex_output_unit,"(A)",IOSTAT=read_stat) line
        nline_texfile = nline_texfile + 1        
        IF (read_stat < 0 ) THEN
          EXIT tail                
        ENDIF
                
        WRITE(file_unit,"(A)") line
                
      ENDDO tail           
   
 
      CLOSE(tex_output_unit)      

      
      RETURN
      END SUBROUTINE write_latex_ps_body           
    
      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE write_char_array(file_unit,nlines,file_lines)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: file_unit
      INTEGER, INTENT(IN) :: nlines      
      TYPE(char_array), DIMENSION(:), INTENT(IN) :: file_lines
      
      INTEGER :: i
      
      DO i = 1,nlines
        WRITE(file_unit,"(A)") file_lines(i)%line      
      ENDDO
      
      RETURN
      END SUBROUTINE write_char_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      SUBROUTINE remove_latex_files()
      
      IMPLICIT NONE
      
      CALL SYSTEM("rm labels.aux labels.dvi labels.log labels.ps labels.tex")
      
      RETURN
      END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      END MODULE labels_mod