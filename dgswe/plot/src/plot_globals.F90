      MODULE plot_globals
      
      USE globals, ONLY: rp
      
      CHARACTER(100) :: input_path
      CHARACTER(100) :: cmap_file
      CHARACTER(100) :: tsnap_spec
      CHARACTER(3) :: frmt
      CHARACTER(3) :: density
      CHARACTER(5) :: font
      
      INTEGER :: substitute_path
      CHARACTER(100) :: replace_path
      CHARACTER(100) :: sub_path
      
!       INTEGER :: coord_sys
!       REAL(rp) :: slam0,sphi0      
      REAL(rp) :: xmin,xmax
      REAL(rp) :: ymin,ymax
      REAL(rp) :: rmin,rmax
      REAL(rp) :: smin,smax
      
      REAL(rp) :: lr_margin 
      REAL(rp) :: top_margin
      REAL(rp) :: cscale_width
      REAL(rp) :: axes_width
      REAL(rp) :: axes_height
      REAL(rp) :: rmin_axes,rmax_axes
      REAL(rp) :: smin_axes,smax_axes    
      REAL(rp) :: rmin_cbar,rmax_cbar
      REAL(rp) :: smin_cbar,smax_cbar   
      REAL(rp) :: rmin_tbar,rmax_tbar
      REAL(rp) :: smin_tbar,smax_tbar        
      REAL(rp) :: ax,bx  
      REAL(rp) :: ay,by      
      
      INTEGER :: ncolors      
      REAL(rp), DIMENSION(:,:), ALLOCATABLE :: colors        
      
      REAL(rp) :: dash      
      REAL(rp) :: xticklabel_pad,yticklabel_pad
      REAL(rp) :: xlabel_pad,ylabel_pad  
      REAL(rp) :: cticklabel_pad
      REAL(rp) :: clabel_pad     
      REAL(rp) :: dr_xlabel,ds_ylabel,ds_clabel
      
      REAL(rp) :: rmin_page = 0d0
      REAL(rp) :: rmax_page = 612d0
      REAL(rp) :: smin_page = 0d0
      REAL(rp) :: smax_page = 792d0      
      
      INTEGER :: nxtick
      INTEGER :: nytick
      INTEGER :: nctick
      
      INTEGER :: nxdec
      INTEGER :: nydec
      INTEGER :: ncdec
      INTEGER :: ntdec
      
      
      CHARACTER(100) :: main_font = "/Times-Roman"
      CHARACTER(100) :: math_font = "/Times-Italic"    
!       CHARACTER(100) :: main_font = "(/usr/share/fonts/type1/gsfonts/cmr10.pfb)"
      INTEGER :: fontsize          
      
      INTEGER :: el,nd,pt
      INTEGER :: et,nv

      
      INTEGER :: adapt_option
      INTEGER :: ps
      INTEGER :: pc
      INTEGER :: p_low
      INTEGER :: p_high
      INTEGER :: p_skip
      INTEGER :: nord
      INTEGER, DIMENSION(:), ALLOCATABLE :: pplt
      INTEGER, DIMENSION(:), ALLOCATABLE :: npplt
      INTEGER, DIMENSION(:), ALLOCATABLE :: nptri
      INTEGER :: mnpp
      INTEGER :: nred 
      
      
      INTEGER, DIMENSION(:), ALLOCATABLE :: ndof_sol      
      REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: phi_sol 
      INTEGER, DIMENSION(:), ALLOCATABLE :: ndof_hb      
      REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: phi_hb

      
      INTEGER :: i,j
      INTEGER :: nsnap
      INTEGER :: snap
      INTEGER :: space
      INTEGER :: n,ndf,npts,nnd,dof
      REAL(rp) :: xpt,ypt
      REAL(rp), TARGET :: t_snap
      REAL(rp) :: t_start,t_end
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: t      
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: r,s  
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: rect        
      REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: xyplt  
      REAL(rp), DIMENSION(:), ALLOCATABLE :: el_area
      INTEGER :: outside
      INTEGER, DIMENSION(:), ALLOCATABLE :: el_in
      REAL(rp) :: xbox_min,xbox_max,ybox_min,ybox_max
      REAL(rp) :: figure_width
      REAL(rp) :: figure_height
      
!       TYPE :: viz
!         INTEGER, DIMENSION(:), ALLOCATABLE :: p
!         INTEGER, DIMENSION(:), ALLOCATABLE :: nnds
!         REAL(rp), DIMENSION(:,:), ALLOCATABLE :: r,s
!         INTEGER, DIMENSION(:), ALLOCATABLE :: ntri        
!         INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: tri
!         REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: xy
!         REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: phi_sol 
!         REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: phi_hb
!       END TYPE
!       
!       TYPE(viz) :: plt
      
      INTEGER :: snap_start
      INTEGER :: snap_end
      
      INTEGER :: nsnap_Z
      INTEGER :: nsnap_Qx
      INTEGER :: nsnap_Qy
      INTEGER :: nsnap_hb
      
      INTEGER :: tex_unit   
           
      CHARACTER(4) :: snap_char     
      CHARACTER(3) :: start_num
      CHARACTER(3) :: nframes
      
      TYPE :: char_array
        CHARACTER(500) :: line
      END TYPE
            
      

      
      REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: Z
      REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: hb
      REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: Qx
      REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: Qy
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: hbm
            
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: eta
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: uu2
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: vv2
      
      TYPE :: plot_type
      
        INTEGER :: plot_sol_option
        INTEGER :: plot_mesh_option
        INTEGER :: plot_lines_option
        INTEGER :: plot_sta_option
        INTEGER :: sta_start
        INTEGER :: sta_end
        
        INTEGER :: cbar_flag
        INTEGER :: tbar_flag
        
        INTEGER :: rm_ps
        INTEGER :: movie_flag        
        
        CHARACTER(100) :: sol_label
        CHARACTER(100) :: name
        
        INTEGER :: ps_unit
        REAL(rp) :: sol_min,sol_max
        REAL(rp) :: snap_min,snap_max
        INTEGER :: nsnap
        INTEGER, DIMENSION(:), ALLOCATABLE :: el_plt
        REAL(rp), DIMENSION(:,:), ALLOCATABLE :: sol_val
        REAL(rp), DIMENSION(:), ALLOCATABLE :: t
        REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: phi 
        INTEGER :: p
        INTEGER :: ndof(4)         
        REAL(rp) :: h0
        
        REAL(rp) :: rel_tol
        REAL(rp) :: abs_tol
      
        INTEGER :: cscale_unit
        CHARACTER(100) :: cscale_option
        INTEGER :: num_cscale_vals
        REAL(rp), DIMENSION(:,:), ALLOCATABLE :: cscale_vals
        REAL(rp) :: cscale_min,cscale_max
        
        CHARACTER(100) :: el_label_option    
        CHARACTER(100) :: nd_label_option 
        
        INTEGER :: type_flag 
        
        REAL(rp), POINTER :: t_snap

        INTEGER :: nline_header
        TYPE(char_array), DIMENSION(:), ALLOCATABLE :: latex_header
        INTEGER :: nline_body
        TYPE(char_array), DIMENSION(:), ALLOCATABLE :: latex_body        
      
      
        REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: sta_val
        REAL(rp), DIMENSION(:,:), ALLOCATABLE :: t_sta        
      END TYPE
      
      TYPE(plot_type) :: zeta
      TYPE(plot_type) :: vel
      TYPE(plot_type) :: bathy
      TYPE(plot_type) :: mesh  
      
!       TYPE(plot_type), DIMENSION(:), ALLOCATABLE :: zeta_sta
!       TYPE(plot_type), DIMENSION(:), ALLOCATABLE :: vel_sta
      
      INTEGER :: nadc
      INTEGER :: ndg
      
      TYPE(char_array), DIMENSION(:), ALLOCATABLE :: adc_sol
      TYPE(char_array), DIMENSION(:), ALLOCATABLE :: dg_sol
      


      
      CONTAINS
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         
      
      SUBROUTINE sizes()
      
      USE globals, ONLY: nverts,ndof,mndof,np,mnp,nnds,mnnds
      USE read_dginp, ONLY: p,ctp,hbp
      
      IMPLICIT NONE
            
      ndof(1) = (p+1)*(p+2)/2
      ndof(2) = (p+1)**2
      ndof(3) = ndof(1)
      ndof(4) = ndof(2)      
      mndof = maxval(ndof)
      
      
      nverts(1) = 3
      nverts(2) = 4
      nverts(3) = 3
      nverts(4) = 4
      
      np(1) = 1
      np(2) = 1
      np(3) = ctp
      np(4) = ctp  
      mnp = maxval(np)+1

      nnds(1) = 3
      nnds(2) = 4
      nnds(3) = (ctp+1)*(ctp+2)/2
      nnds(4) = (ctp+1)*(ctp+1)      
      mnnds = maxval(nnds)   
      
!       pplt(1) = ps
!       pplt(2) = ps
!       pplt(3) = pc
!       pplt(4) = pc

!       npplt(1) = (ps+1)*(ps+2)/2
!       npplt(2) = (ps+1)*(ps+1)
!       npplt(3) = (pc+1)*(pc+2)/2
!       npplt(4) = (pc+1)*(pc+1)
!       mnpp = maxval(npplt)     

      mnpp = (p_high+1)**2
      
      END SUBROUTINE sizes
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      SUBROUTINE setup_plot_types()
      USE read_dginp, ONLY: h0,p,hbp   
      
      IMPLICIT NONE           
      
      mesh%type_flag = 1     
      bathy%type_flag = 2    
      zeta%type_flag = 3     
      vel%type_flag = 4
      
      mesh%cbar_flag = 0
      mesh%tbar_flag = 0
      bathy%cbar_flag = 1
      bathy%tbar_flag = 0
      zeta%cbar_flag = 1
      zeta%tbar_flag = 1
      vel%cbar_flag = 1
      vel%tbar_flag = 1
      
      bathy%el_label_option = "off"
      bathy%nd_label_option = "off"
      zeta%el_label_option = "off"
      zeta%nd_label_option = "off"
      vel%el_label_option = "off"
      vel%nd_label_option = "off"
      
      mesh%sol_label = "mesh"
      mesh%name = "mesh"
      zeta%sol_label = "surface elevation (m)"
      zeta%name = "zeta"      
      vel%sol_label = "velocity (m/s)"
      vel%name = "vel"      
      bathy%sol_label = "bathymetry (m)"
      bathy%name = "bathy"           
      
      mesh%t_snap => t_snap
      bathy%t_snap => t_snap
      zeta%t_snap => t_snap
      vel%t_snap => t_snap            
      
      zeta%cscale_unit = 30
      bathy%cscale_unit = 31
      vel%cscale_unit = 32
      
      bathy%h0 = h0
      
      bathy%p = hbp
      zeta%p = p
      vel%p = p
      
      END SUBROUTINE setup_plot_types

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      END MODULE plot_globals