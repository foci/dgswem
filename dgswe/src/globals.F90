      MODULE globals
      IMPLICIT NONE

      SAVE

      INTEGER, PARAMETER :: rp = kind(1d0) ! double precision 
      INTEGER :: ne ! number of elements
      INTEGER :: nn ! number of nodes
      REAL(rp), PARAMETER :: g = 9.81d0 ! gravitational constant
      REAL(rp), PARAMETER :: pt5g = 0.5d0*g
      REAL(rp), PARAMETER :: pi=3.141592653589793D0
      REAL(rp), PARAMETER :: deg2rad = pi/180d0      
!       REAL(rp), PARAMETER :: r_earth = 6378206.4d0      
      REAL(rp), PARAMETER :: r_earth = 6.340304248283833d6    
      CHARACTER(100) :: grid_name ! name of the grid

      INTEGER, PARAMETER :: nel_type = 4 !(type #s: 1 -> triangles, 2 -> quads, 3 -> curved triangles, 4-> curved quads)
      INTEGER, PARAMETER :: norder = 6 ! # of different orders (straight sided elements for tri/quad = 1, curvilinear tri/quad = ctp, high-order bathymetry = hbp) 
      INTEGER :: nverts(nel_type)
      INTEGER :: np(norder), nnds(norder)
      INTEGER :: order(2*nel_type)
      
      INTEGER, ALLOCATABLE, DIMENSION(:) :: el_type
      
      INTEGER, DIMENSION(nel_type) :: ndof ! number of degrees of freedom
      INTEGER, DIMENSION(nel_type) :: nqpta ! number of area quadrature points
      INTEGER, DIMENSION(nel_type) :: nqpte ! number of edge quadrature points      
      INTEGER :: mnqpta,mnnds,mnqpte,mndof,mnp      
      
      INTEGER :: nblk,nrblk
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: edblk,nfblk,rnfblk
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: iediblk,bediblk
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: elblk
      INTEGER :: mnpartel,mnparted      
      
      INTEGER, ALLOCATABLE, DIMENSION(:) :: part
      INTEGER, ALLOCATABLE, DIMENSION(:) :: npartel
      INTEGER, ALLOCATABLE, DIMENSION(:) :: nparted
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: npartet      
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: lel2gel ! gives the global element number that corresponds to a local partition element number
      INTEGER, ALLOCATABLE, DIMENSION(:) :: gel2lel
      INTEGER, ALLOCATABLE, DIMENSION(:) :: gel2part      
      INTEGER, ALLOCATABLE, DIMENSION(:) :: gel2ael ! gives the aligned element number that corresponds to a global element number
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ael2gel ! gives the global element number that corresponds to an aligned element number      

      REAL(rp) :: t                ! simulation time
      INTEGER :: tskp_sol,tskp_sta ! time steps to skip between output
      INTEGER :: nout_sol,nout_sta ! number of outputs 
      REAL(rp) :: tstage           ! rk stage time
      REAL(rp) :: ramp             ! boundary condition ramp value
      REAL(rp), PARAMETER :: pt3333 = 1d0/3d0 ! 1/3 for 3rd order rk
      
      REAL(rp), DIMENSION(5), PARAMETER :: ark = (/ 0d0, &
                                                   -0.41789047449985196d0, &
                                                   -1.19215169464267693d0, &
                                                   -1.69778469247152784d0, &
                                                   -1.51418344425715578d0 /)
      REAL(rp), DIMENSION(5), PARAMETER :: brk = (/ 0.149659021999229117d0, &
                                                    0.379210312999627281d0, &
                                                    0.822955029386981717d0, &
                                                    0.699450455949122107d0, &
                                                    0.153057247968151993d0 /)
      REAL(rp), DIMENSION(5), PARAMETER :: crk = (/ 0d0, &
                                                    0.149659021999229117d0, &
                                                    0.370400957364204773d0, &
                                                    0.622255763134443168d0, &
                                                    0.958282130674690254d0 /)           

      
      INTEGER :: nsta
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ndsta
      INTEGER, ALLOCATABLE, DIMENSION(:) :: elsta
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: xysta ! x,y coordinates of stations
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: phi_sta      
      
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ect ! element connectivity table
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: xy ! x,y coordinates of nodes
      REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: elxy            
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: depth ! depth at each node
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: elhb
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: hbnodes 
      REAL(rp), ALLOCATABLE, DIMENSION(:,:,:,:) :: bndxy
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: el_size

      INTEGER :: nope ! number of open boundary segents
      INTEGER, ALLOCATABLE, DIMENSION(:) :: obseg ! number of nodes in each open boundary segment
      INTEGER :: neta ! total elevation specified boundary nodes
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: obnds ! open boundary nodes (node,seg)

      INTEGER :: nbou  ! number of normal flow boundary segments
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: fbseg ! number of nodes and type of each normal flow boundary segment (1,seg) = # nodes, (2,seg) = type
      INTEGER :: nvel  ! total number of normal flow boundary nodes
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: fbnds ! normal flow boundary nodes (node,seg)

      INTEGER :: nobfr ! number of periodic forcings on elevation boundaries
      CHARACTER(10), ALLOCATABLE, DIMENSION(:) :: obtag,obtag2 ! constituent name
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: obfreq ! open bounday constituent frequency
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: obper
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: obnfact ! open boundary constituent nodal factor
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: obeq ! open bounday constituent equilibrium argument
      REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: obamp ! open boundary node amplitute
      REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: obamp_qpt !open boundary amplitude interpolated to edge quatrature points 
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: obdepth_qpt
      REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: obph ! open boundary node phase
      REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: obph_qpt ! open boundary phase interpolated to edge quadrature points


      INTEGER :: nfbfr ! number of periodic forcings on flow boundaries
      CHARACTER(10), ALLOCATABLE, DIMENSION(:) :: fbtag,fbtag2 ! constituent name
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: fbfreq ! flow bounday constituent frequency
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: fbper
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: fbnfact ! flow boundary constituent nodal factor
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: fbeq ! flow bounday constituent equilibrium argument
      REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: fbamp ! flow boundary node amplitute  
      REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: fbamp_qpt ! flow boundary amplitude interpolated to edge quadrature points
      REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: fbph ! flow boundary node phase
      REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: fbph_qpt ! flow boundary phase interpolated to edge quadrature points
      
      INTEGER :: nfbsfr ! number of surge forcings on flow boundaries
      CHARACTER(10), ALLOCATABLE, DIMENSION(:) :: fbstag,fbstag2 ! surge name      
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: fbsbgn ! flow boundary surge beginning time
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: fbsend ! flow boundary surge ending time
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: fbssig ! flow boundary surge rise/fall parameter
      REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: fbsamp ! flow boundary nodal surge amplitude 
      REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: fbsamp_qpt ! flow boundary nodal surge amplitude interpolated to edge quadrature points

      REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: Va,Ve
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ipive
      REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: psia,psie,psiv,psic
      REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: dpsidr,dpsids   
      REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: dpsidxi
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: area ! element areas
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: edlen ! element edge lengths
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: edlen_area
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: normal ! element edge normals
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: nx_pt,ny_pt    
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: Spe,cfac
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: dhbdx,dhbdy ! elemental x and y derivatives of linear bathymetry
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: dhbdx_init,dhbdy_init   
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: mmi,mmi_init      
      REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: m2n     
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: ml2,mml

      INTEGER :: ned ! total number of edges
      INTEGER, ALLOCATABLE, DIMENSION(:) :: nepn ! number of elements per node
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: epn ! elements per node 
      INTEGER :: mnepn ! maximum number of elements per node      
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ged2nn ! gives the two node numbers that make up a global edge number
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ged2el ! gives the two element numbers that share a global edge number
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ged2led ! gives the two local edge numbers that share a global edge number
      INTEGER, ALLOCATABLE, DIMENSION(:) :: recv_edge    
      
      INTEGER :: nied ! total number of interior edges
      INTEGER, ALLOCATABLE, DIMENSION(:) :: iedn ! array of interior edge numbers
      INTEGER :: nbed
      INTEGER, ALLOCATABLE, DIMENSION(:) :: bedn
      INTEGER :: nobed ! total number of open boundary edges
      INTEGER, ALLOCATABLE, DIMENSION(:) :: obedn ! array of open boundary edge numbers   
      INTEGER :: nfbed ! total number of flow boundary edges
      INTEGER, ALLOCATABLE, DIMENSION(:) :: fbedn ! array of flow boundary edge numbers
      INTEGER :: nnfbed ! total number of no normal flow boundary edges
      INTEGER, ALLOCATABLE, DIMENSION(:) :: nfbedn ! array of no normal flow boundary edge numbers
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: nfbednn       
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ed_type  ! open/flow boundary edge flags
      
      INTEGER :: check_iedge,check_gedge

      REAL(rp), ALLOCATABLE, TARGET, DIMENSION(:,:) :: H ! degrees of freedom for water column height
      REAL(rp), ALLOCATABLE, TARGET, DIMENSION(:,:) :: Z      
      REAL(rp), ALLOCATABLE, TARGET, DIMENSION(:,:) :: Qx ! degrees of freedom for x momentum
      REAL(rp), ALLOCATABLE, TARGET, DIMENSION(:,:) :: Qy ! degrees of freedom for y momentum
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: hbm                 
      
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: Exx
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: Eyy
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: Exy   
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: Eyx      
      
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: Hold,Hinit ! degrees of freedom for water column height
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: Zold,Zinit 
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: Qxold,Qxinit ! degrees of freedom for x momentum
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: Qyold,Qyinit ! degrees of freedom for y momentum
      
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: hbqpta,hbqpta_init
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: hbqpte,hbqpte_init
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: hbqpted      


      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: rhsZ,rhsH,rhsQx,rhsQy ! right hand side evaluation arrays
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: MirhsZ,MirhsH,MirhsQx,MirhsQy        
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: rhsExx,rhsEyy,rhsExy,rhsEyx

      REAL(rp), ALLOCATABLE, DIMENSION(:) :: Zqpta,Hqpta,Qxqpta,Qyqpta
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: xmoma,ymoma,xymoma
      REAL(rp), ALLOCATABLE, TARGET, DIMENSION(:,:) :: Zqpt,Hqpt,Qxqpt,Qyqpt ! solution quadrature point evaluation arrays
      REAL(rp), ALLOCATABLE, TARGET, DIMENSION(:,:) :: xmom,ymom,xymom ! momentum terms quadrature point evaluation arrays
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: src_x,src_y
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: tau
      
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: Exxqpta,Eyyqpta,Exyqpta,Eyxqpta             
      REAL(rp), ALLOCATABLE, TARGET, DIMENSION(:,:) :: Exxqpt,Eyyqpt,Exyqpt,Eyxqpt         

      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: wpta ! area and edge quadrature weights
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: wpte      
      REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: qpte ! edge quadrature points
      REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: qpta ! area quadrature points

      REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: phia ! basis functions evaluated at area quadrature points
      REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: phia_int,phia_int_init ! basis functions evaluated at area quadrature points multiplied by quadrature weights
      REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: phie ! basis functions evaluated at edge quadrature points
      REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: phie_int 
      REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: phi
      REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: phil ! linear nodal basis functions evaluated at area quadrature points
      REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: dpdr,dpds ! basis function derivatives evaluated at area quadrature points
      REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: dpdx,dpdy            
      REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: dpdx_init,dpdy_init
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: detJa     
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: detJe            

      REAL(rp), ALLOCATABLE, DIMENSION(:) :: recipHa ! temporary variables for momentum term calculation

      REAL(rp), ALLOCATABLE, DIMENSION(:) :: const
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: inx,iny
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: icfac      
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: detJe_in,detJe_ex
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: Zhatv,Hhatv,Qxhatv,Qyhatv
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: rHi,rHe
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: xmomi,xmome,ymomi,ymome,xymomi,xymome     
      REAL(rp), ALLOCATABLE, TARGET, DIMENSION(:,:) :: fbHf,nfbHf,obHf  
      REAL(rp), ALLOCATABLE, TARGET, DIMENSION(:,:) :: fbQxf,nfbQxf,obQxf       
      REAL(rp), ALLOCATABLE, TARGET, DIMENSION(:,:) :: fbQyf,nfbQyf,obQyf  
      REAL(rp), ALLOCATABLE, TARGET, DIMENSION(:,:) :: Hfluxi,Hfluxe     
      REAL(rp), ALLOCATABLE, TARGET, DIMENSION(:,:) :: Zfluxi,Zfluxe      
      REAL(rp), ALLOCATABLE, TARGET, DIMENSION(:,:) :: Qxfluxi,Qxfluxe
      REAL(rp), ALLOCATABLE, TARGET, DIMENSION(:,:) :: Qyfluxi,Qyfluxe      

      TYPE :: edge_ptr_array
        REAL(rp), POINTER :: ptr
      END TYPE edge_ptr_array

      TYPE(edge_ptr_array), ALLOCATABLE, DIMENSION(:,:) :: Hi,He
      TYPE(edge_ptr_array), ALLOCATABLE, DIMENSION(:,:) :: Zi,Ze      
      TYPE(edge_ptr_array), ALLOCATABLE, DIMENSION(:,:) :: Qxi,Qxe
      TYPE(edge_ptr_array), ALLOCATABLE, DIMENSION(:,:) :: Qyi,Qye
      TYPE(edge_ptr_array), ALLOCATABLE, DIMENSION(:,:) :: xmi,xme
      TYPE(edge_ptr_array), ALLOCATABLE, DIMENSION(:,:) :: ymi,yme
      TYPE(edge_ptr_array), ALLOCATABLE, DIMENSION(:,:) :: xymi,xyme   
      
      TYPE(edge_ptr_array), ALLOCATABLE, DIMENSION(:,:) :: Exxi,Exxe
      TYPE(edge_ptr_array), ALLOCATABLE, DIMENSION(:,:) :: Eyyi,Eyye
      TYPE(edge_ptr_array), ALLOCATABLE, DIMENSION(:,:) :: Exyi,Exye            
      TYPE(edge_ptr_array), ALLOCATABLE, DIMENSION(:,:) :: Eyxi,Eyxe
       
      
      TYPE(edge_ptr_array), ALLOCATABLE, DIMENSION(:,:) :: Hwrite
      TYPE(edge_ptr_array), ALLOCATABLE, DIMENSION(:,:) :: Zwrite      
      TYPE(edge_ptr_array), ALLOCATABLE, DIMENSION(:,:) :: Qxwrite
      TYPE(edge_ptr_array), ALLOCATABLE, DIMENSION(:,:) :: Qywrite
      
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: Zout 
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: Qxout 
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: Qyout 
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: Zsta
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: Qxsta
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: Qysta
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: hbsta      
      
      INTEGER :: Zsol_unit
      INTEGER :: Qxsol_unit
      INTEGER :: Qysol_unit
      INTEGER :: Zsta_unit
      INTEGER :: Qxsta_unit
      INTEGER :: Qysta_unit
      
      
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: Qxin,Qyin,Hin,Zin
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: Qxex,Qyex,Hex,Zex
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: xmin,ymin,xymin
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: xmex,ymex,xymex
      
      
      INTEGER, ALLOCATABLE, DIMENSION(:) :: nresel
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: el_g2l,el_l2g

      INTEGER :: nsred
      INTEGER, ALLOCATABLE, DIMENSION(:) :: sredn
      
      INTEGER :: mned_sr
      
      INTEGER, ALLOCATABLE, DIMENSION(:) :: nresnd
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: nd_l2g      
      INTEGER, ALLOCATABLE, DIMENSION(:) :: nd_g2l
       
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: lect
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: lnelnds
      
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: lobseg
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: lobnds
      INTEGER, ALLOCATABLE, DIMENSION(:) :: lneta
      INTEGER, ALLOCATABLE, DIMENSION(:) :: lnope      
      
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: lfbseg
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: lfbnds      
      INTEGER, ALLOCATABLE, DIMENSION(:) :: lnvel
      INTEGER, ALLOCATABLE, DIMENSION(:) :: lnbou     
      REAL(rp), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: lbndxy
      
      REAL(rp), ALLOCATABLE, DIMENSION(:,:,:,:) :: lobamp      
      REAL(rp), ALLOCATABLE, DIMENSION(:,:,:,:) :: lobph

      REAL(rp), ALLOCATABLE, DIMENSION(:,:,:,:) :: lfbamp
      REAL(rp), ALLOCATABLE, DIMENSION(:,:,:,:) :: lfbph 
      REAL(rp), ALLOCATABLE, DIMENSION(:,:,:,:) :: lfbsamp
      INTEGER, ALLOCATABLE, DIMENSION(:) :: lnbouf  
      
      INTEGER, ALLOCATABLE, DIMENSION(:) :: nlsta
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: sta_l2g
      
      
      
      TYPE :: grid
      
        CHARACTER(100) :: grid_name ! name of the grid
        CHARACTER(100) :: grid_file ! name of fort.14 file 
        
        INTEGER :: ne ! number of elements
        INTEGER :: nn ! number of nodes                
      
        INTEGER, ALLOCATABLE, DIMENSION(:) :: el_type      
      
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ect ! element connectivity table
        REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: xy ! x,y coordinates of nodes
        REAL(rp), ALLOCATABLE, DIMENSION(:) :: depth ! depth at each node
 
        INTEGER :: nope ! number of open boundary segents
        INTEGER, ALLOCATABLE, DIMENSION(:) :: obseg ! number of nodes in each open boundary segment
        INTEGER :: neta ! total elevation specified boundary nodes
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: obnds ! open boundary nodes

        INTEGER :: nbou  ! number of normal flow boundary segments
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: fbseg ! number of nodes and type of each normal flow boundary segment
        INTEGER :: nvel  ! total number of normal flow boundary nodes
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: fbnds ! normal flow boundary nodes
      
        INTEGER :: mnepn
        INTEGER, ALLOCATABLE, DIMENSION(:) :: nepn ! number of elements per node
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: epn ! elements per node                      

      END TYPE
         
      
      TYPE(grid) :: base
      TYPE(grid) :: eval
      
#ifdef ALIGN16      
!DIR$ ATTRIBUTES ALIGN:16 :: rhsH,rhsQx,rhsQy
!DIR$ ATTRIBUTES ALIGN:16 :: Hqpt,Qxqpt,Qyqpt
!DIR$ ATTRIBUTES ALIGN:16 :: H,Qx,Qy
!DIR$ ATTRIBUTES ALIGN:16 :: Hold, Qxold,Qyold
!DIR$ ATTRIBUTES ALIGN:16 :: phia,phie
!DIR$ ATTRIBUTES ALIGN:16 :: recipHa
!DIR$ ATTRIBUTES ALIGN:16 :: xmom,ymom,xymom
!DIR$ ATTRIBUTES ALIGN:16 :: tau
!DIR$ ATTRIBUTES ALIGN:16 :: dhbdx,dhbdy
!DIR$ ATTRIBUTES ALIGN:16 :: src_x,src_y
!DIR$ ATTRIBUTES ALIGN:16 :: dpdx,dpdy
!DIR$ ATTRIBUTES ALIGN:16 :: phia_int, phie_int
!DIR$ ATTRIBUTES ALIGN:16 :: Hflux,Qxflux,Qyflux
!DIR$ ATTRIBUTES ALIGN:16 :: const,Hhatv,Qxhatv,Qyhatv
!DIR$ ATTRIBUTES ALIGN:16 :: Hi,He,Qxi,Qxe,Qyi,Qye,xmi,xme,xymi,xyme,ymi,yme
!DIR$ ATTRIBUTES ALIGN:16 :: inx,iny,len_area_ex,len_area_in
!DIR$ ATTRIBUTES ALIGN:16 :: Hfe,Hfi,Qxfi,Qxfe,Qyfi,Qyfe
!DIR$ ATTRIBUTES ALIGN:16 :: Hin,Qxin,Qyin
!DIR$ ATTRIBUTES ALIGN:16 :: Hex,Qxex,Qyex
!DIR$ ATTRIBUTES ALIGN:16 :: xmin,ymin,xymin
!DIR$ ATTRIBUTES ALIGN:16 :: xmex,ymex,xymex      
#endif

#ifdef ALIGN32      
!DIR$ ATTRIBUTES ALIGN:32 :: rhsH,rhsQx,rhsQy
!DIR$ ATTRIBUTES ALIGN:32 :: Hqpt,Qxqpt,Qyqpt
!DIR$ ATTRIBUTES ALIGN:32 :: H,Qx,Qy
!DIR$ ATTRIBUTES ALIGN:32 :: Hold, Qxold,Qyold
!DIR$ ATTRIBUTES ALIGN:32 :: phia,phie
!DIR$ ATTRIBUTES ALIGN:32 :: recipHa
!DIR$ ATTRIBUTES ALIGN:32 :: xmom,ymom,xymom
!DIR$ ATTRIBUTES ALIGN:32 :: tau
!DIR$ ATTRIBUTES ALIGN:32 :: dhbdx,dhbdy
!DIR$ ATTRIBUTES ALIGN:32 :: src_x,src_y
!DIR$ ATTRIBUTES ALIGN:32 :: dpdx,dpdy
!DIR$ ATTRIBUTES ALIGN:32 :: phia_int, phie_int
!DIR$ ATTRIBUTES ALIGN:32 :: Hflux,Qxflux,Qyflux
!DIR$ ATTRIBUTES ALIGN:32 :: const,Hhatv,Qxhatv,Qyhatv
!DIR$ ATTRIBUTES ALIGN:32 :: Hi,He,Qxi,Qxe,Qyi,Qye,xmi,xme,xymi,xyme,ymi,yme
!DIR$ ATTRIBUTES ALIGN:32 :: inx,iny,len_area_ex,len_area_in
!DIR$ ATTRIBUTES ALIGN:32 :: Hfe,Hfi,Qxfi,Qxfe,Qyfi,Qyfe
!DIR$ ATTRIBUTES ALIGN:32 :: Hin,Qxin,Qyin
!DIR$ ATTRIBUTES ALIGN:32 :: Hex,Qxex,Qyex
!DIR$ ATTRIBUTES ALIGN:32 :: xmin,ymin,xymin
!DIR$ ATTRIBUTES ALIGN:32 :: xmex,ymex,xymex
#endif

#ifdef ALIGN64
!DIR$ ATTRIBUTES ALIGN:64 :: rhsH,rhsQx,rhsQy
!DIR$ ATTRIBUTES ALIGN:64 :: Hqpt,Qxqpt,Qyqpt
!DIR$ ATTRIBUTES ALIGN:64 :: H,Qx,Qy
!DIR$ ATTRIBUTES ALIGN:64 :: Hold, Qxold,Qyold
!DIR$ ATTRIBUTES ALIGN:64 :: phia,phie
!DIR$ ATTRIBUTES ALIGN:64 :: recipHa
!DIR$ ATTRIBUTES ALIGN:64 :: xmom,ymom,xymom
!DIR$ ATTRIBUTES ALIGN:64 :: tau
!DIR$ ATTRIBUTES ALIGN:64 :: dhbdx,dhbdy
!DIR$ ATTRIBUTES ALIGN:64 :: src_x,src_y
!DIR$ ATTRIBUTES ALIGN:64 :: dpdx,dpdy
!DIR$ ATTRIBUTES ALIGN:64 :: phia_int, phie_int
!DIR$ ATTRIBUTES ALIGN:64 :: Hflux,Qxflux,Qyflux
!DIR$ ATTRIBUTES ALIGN:64 :: const,Hhatv,Qxhatv,Qyhatv
!DIR$ ATTRIBUTES ALIGN:64 :: Hi,He,Qxi,Qxe,Qyi,Qye,xmi,xme,xymi,xyme,ymi,yme
!DIR$ ATTRIBUTES ALIGN:64 :: inx,iny,len_area_ex,len_area_in
!DIR$ ATTRIBUTES ALIGN:64 :: Hfe,Hfi,Qxfi,Qxfe,Qyfi,Qyfe
!DIR$ ATTRIBUTES ALIGN:64 :: Hin,Qxin,Qyin
!DIR$ ATTRIBUTES ALIGN:64 :: Hex,Qxex,Qyex
!DIR$ ATTRIBUTES ALIGN:64 :: xmin,ymin,xymin
!DIR$ ATTRIBUTES ALIGN:64 :: xmex,ymex,xymex
#endif


!       ark(1) = 0d0
!       ark(2) = - 567301805773d0 / 1357537059087d0
!       ark(3) = -2404267990393d0 / 2016746695238d0
!       ark(4) = -3550918686646d0 / 2091501179385d0
!       ark(5) = -1275806237668d0 / 842570457699d0
      
!       brk(1) = 1432997174477d0 / 9575080441755d0
!       brk(2) = 5161836677717d0 / 13612068292357d0
!       brk(3) = 1720146321549d0 / 2090206949498d0
!       brk(4) = 3134564353537d0 / 4481467310338d0
!       brk(5) = 2277821191437d0 / 14882151754819d0
!       
!       crk(1) = 0d0
!       crk(2) = 1432997174477d0 / 9575080441755d0
!       crk(3) = 2526269341429d0 / 6820363962896d0
!       crk(4) = 2006345519317d0 / 3224310063776d0
!       crk(5) = 2802321613138d0 / 2924317926251d0   

     
      END MODULE globals

