      MODULE globals
      IMPLICIT NONE

      SAVE

      INTEGER, PARAMETER :: pres = kind(1d0) ! double precision 
      INTEGER :: p ! polynomial order
      INTEGER :: ndof ! number of degrees of freedom
      INTEGER :: nqpta ! number of area quadrature points
      INTEGER :: nqpte ! number of edge quadrature points
      INTEGER :: ne ! number of elements
      INTEGER :: nn ! number of nodes
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: split,split2
      INTEGER :: sp,nsp,nsp2

      REAL(pres), PARAMETER :: g = 9.81d0 ! gravitational constant
      REAL(pres), PARAMETER :: pt5g = 0.5d0*g
      REAL(pres), PARAMETER  ::  pi=3.141592653589793D0
      REAL(pres) :: cf ! bottom friction parameter
      CHARACTER(24) :: grid_name ! name of the grid
      CHARACTER(24) :: grid_file ! name of fort.14 file
      CHARACTER(24) :: forcing_file ! name of fort.15 file

      REAL(pres) :: dt ! time step
      REAL(pres) :: t ! simulation time
      REAL(pres) :: tstage ! rk stage time
      REAL(pres) :: ramp ! boundary condition ramp value
      REAL(pres) :: dramp ! numer of ramp days
      REAL(pres), PARAMETER :: pt3333 = 1d0/3d0 ! 1/3 for 3rd order rk
      
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ect ! element connectivity table
      REAL(pres), ALLOCATABLE, DIMENSION(:,:) :: xy ! x,y coordinates of nodes
      REAL(pres), ALLOCATABLE, DIMENSION(:) :: depth ! depth at each node

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
      REAL(pres), ALLOCATABLE, DIMENSION(:) :: obfreq ! open bounday constituent frequency
      REAL(pres), ALLOCATABLE, DIMENSION(:) :: obper
      REAL(pres), ALLOCATABLE, DIMENSION(:) :: obnfact ! open boundary constituent nodal factor
      REAL(pres), ALLOCATABLE, DIMENSION(:) :: obeq ! open bounday constituent equilibrium argument
      REAL(pres), ALLOCATABLE, DIMENSION(:,:) :: obamp ! open boundary node amplitute
      REAL(pres), ALLOCATABLE, DIMENSION(:,:) :: obamp_qpt !open boundary amplitude interpolated to edge quatrature points 
      REAL(pres), ALLOCATABLE, DIMENSION(:,:) :: obdepth_qpt
      REAL(pres), ALLOCATABLE, DIMENSION(:,:) :: obph ! open boundary node phase
      REAL(pres), ALLOCATABLE, DIMENSION(:,:) :: obph_qpt ! open boundary phase interpolated to edge quadrature points


      INTEGER :: nfbfr ! number of periodic forcings on flow boundaries
      CHARACTER(10), ALLOCATABLE, DIMENSION(:) :: fbtag,fbtag2 ! constituent name
      REAL(pres), ALLOCATABLE, DIMENSION(:) :: fbfreq ! flow bounday constituent frequency
      REAL(pres), ALLOCATABLE, DIMENSION(:) :: fbper
      REAL(pres), ALLOCATABLE, DIMENSION(:) :: fbnfact ! flow boundary constituent nodal factor
      REAL(pres), ALLOCATABLE, DIMENSION(:) :: fbeq ! flow bounday constituent equilibrium argument
      REAL(pres), ALLOCATABLE, DIMENSION(:,:) :: fbamp ! flow boundary node amplitute  
      REAL(pres), ALLOCATABLE, DIMENSION(:,:) :: fbamp_qpt ! flow boundary amplitude interpolated to edge quadrature points
      REAL(pres), ALLOCATABLE, DIMENSION(:,:) :: fbph ! flow boundary node phase
      REAL(pres), ALLOCATABLE, DIMENSION(:,:) :: fbph_qpt ! flow boundary phase interpolated to edge quadrature points
      REAL(pres) :: arg
      INTEGER :: bfr

      REAL(pres), ALLOCATABLE, DIMENSION(:) :: area ! element areas
      REAL(pres), ALLOCATABLE, DIMENSION(:) :: edlen ! element edge lengths
      REAL(pres), ALLOCATABLE, DIMENSION(:,:) :: edlen_area
      REAL(pres), ALLOCATABLE, DIMENSION(:,:) :: normal ! element edge normals
      REAL(pres), ALLOCATABLE, DIMENSION(:) :: dhbdx,dhbdy ! elemental x and y derivatives of linear bathymetry

      INTEGER :: ned ! total number of edges
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ged2nn ! gives the two node numbers that make up a global edge number
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ged2el ! gives the two element numbers that share a global edge number
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ged2led ! gives the two local edge numbers that share a global edge number

      INTEGER :: nied ! total number of interior edges
      INTEGER, ALLOCATABLE, DIMENSION(:) :: iedn ! array of interior edge numbers
      INTEGER :: nobed ! total number of open boundary edges
      INTEGER, ALLOCATABLE, DIMENSION(:) :: obedn ! array of open boundary edge numbers   
      INTEGER :: nfbed ! total number of flow boundary edges
      INTEGER, ALLOCATABLE, DIMENSION(:) :: fbedn ! array of flow boundary edge numbers
      INTEGER :: nnfbed ! total number of no normal flow boundary edges
      INTEGER, ALLOCATABLE, DIMENSION(:) :: nfbedn ! array of no normal flow boundary edge numbers

      REAL(pres), ALLOCATABLE, DIMENSION(:,:) :: H,Hold ! degrees of freedom for water column height
      REAL(pres), ALLOCATABLE, DIMENSION(:,:) :: Qx,Qxold ! degrees of freedom for x momentum
      REAL(pres), ALLOCATABLE, DIMENSION(:,:) :: Qy,Qyold ! degrees of freedom for y momentum

      REAL(pres), ALLOCATABLE, DIMENSION(:,:) :: rhsH,rhsQx,rhsQy ! right hand side evaluation arrays

      REAL(pres), ALLOCATABLE, TARGET, DIMENSION(:,:) :: Hqpt,Qxqpt,Qyqpt ! solution quadrature point evaluation arrays
      REAL(pres), ALLOCATABLE, TARGET, DIMENSION(:,:) :: xmom,ymom,xymom ! momentum terms quadrature point evaluation arrays
      REAL(pres), ALLOCATABLE, TARGET, DIMENSION(:,:) :: egnval,Hn,Qxn,Qyn
      REAL(pres), ALLOCATABLE, DIMENSION(:) :: src_x,src_y
      REAL(pres), ALLOCATABLE, DIMENSION(:) :: tau
      REAL(pres) :: u,v

      REAL(pres), ALLOCATABLE, DIMENSION(:) :: wpta,wpte ! area and edge quadrature weights
      REAL(pres), ALLOCATABLE, DIMENSION(:) :: qpte ! edge quadrature points
      REAL(pres), ALLOCATABLE, DIMENSION(:,:) :: qpta ! area quadrature points

      REAL(pres), ALLOCATABLE, DIMENSION(:,:) :: phia ! basis functions evaluated at area quadrature points
      REAL(pres), ALLOCATABLE, DIMENSION(:,:) :: phia_int ! basis functions evaluated at area quadrature points multiplied by quadrature weights
      REAL(pres), ALLOCATABLE, DIMENSION(:,:) :: phie ! basis functions evaluated at edge quadrature points
      REAL(pres), ALLOCATABLE, DIMENSION(:,:) :: phie_int 
      REAL(pres), ALLOCATABLE, DIMENSION(:,:) :: phil ! linear nodal basis functions evaluated at area quadrature points
      REAL(pres), ALLOCATABLE, DIMENSION(:,:) :: dpdx,dpdy ! basis function derivatives evaluated at area quadrature points

      REAL(pres), ALLOCATABLE, DIMENSION(:) :: pressa,recipHa ! temporary variables for momentum term calculation
      REAL(pres) :: press,recipH ! temporary variables for momentum term calculation

      INTEGER :: el,dof,l,pt,ind,ed,led,ged ! do loop variables
      INTEGER :: el_in,el_ex ! interior/exterior element variables
      INTEGER :: led_in,led_ex ! interior/exterior local edge variables
      INTEGER :: gp_in,gp_ex

      REAL(pres) :: alpha,eig_in,eig_ex ! numerical flux constant calculation variables
      REAL(pres), DIMENSION(6) :: eig
      REAL(pres) :: H_in,H_ex,Qx_in,Qx_ex,Qy_in,Qy_ex ! interior/exterior numerical flux solution variables
      REAL(pres) :: xmom_in,xmom_ex,ymom_in,ymom_ex,xymom_in,xymom_ex ! interior/exterior numerical flux momentum variables
      REAL(pres) :: nx,ny,nx2,ny2,nxny,tx,ty ! x and y edge normals
      REAL(pres), ALLOCATABLE, DIMENSION(:,:) :: nxv,nyv
      REAL(pres) :: Hhat,Qxhat,Qyhat ! numerical flux variables
      REAL(pres), ALLOCATABLE, TARGET, DIMENSION(:,:) :: Hflux,Qxflux,Qyflux ! numerical flux arrays
      REAL(pres) :: Qn,Qt

      REAL(pres), ALLOCATABLE, DIMENSION(:) :: const
      REAL(pres), ALLOCATABLE, DIMENSION(:) :: inx,iny,len_area_in,len_area_ex
      REAL(pres), ALLOCATABLE, DIMENSION(:) :: Hhatv,Qxhatv,Qyhatv
      
      REAL(pres), ALLOCATABLE, DIMENSION(:) :: Hin,Hex,Qxin,Qxex,Qyin,Qyex
      REAL(pres), ALLOCATABLE, DIMENSION(:) :: xmin,xmex,ymin,ymex,xymin,xymex   
      
      INTEGER, ALLOCATABLE, DIMENSION(:) :: edpt_in,edpt_ex

      TYPE :: edge_ptr_array
        REAL(pres), POINTER :: ptr
      END TYPE edge_ptr_array

      TYPE(edge_ptr_array), ALLOCATABLE, DIMENSION(:,:) :: Hi,He
      TYPE(edge_ptr_array), ALLOCATABLE, DIMENSION(:,:) :: Qxi,Qxe
      TYPE(edge_ptr_array), ALLOCATABLE, DIMENSION(:,:) :: Qyi,Qye
      TYPE(edge_ptr_array), ALLOCATABLE, DIMENSION(:,:) :: xmi,xme
      TYPE(edge_ptr_array), ALLOCATABLE, DIMENSION(:,:) :: ymi,yme
      TYPE(edge_ptr_array), ALLOCATABLE, DIMENSION(:,:) :: xymi,xyme
      TYPE(edge_ptr_array), ALLOCATABLE, DIMENSION(:,:) :: Hfi,Qxfi,Qyfi
      TYPE(edge_ptr_array), ALLOCATABLE, DIMENSION(:,:) :: Hfe,Qxfe,Qyfe
      TYPE(edge_ptr_array), ALLOCATABLE, DIMENSION(:,:) :: Hni,Qxni,Qyni
      TYPE(edge_ptr_array), ALLOCATABLE, DIMENSION(:,:) :: Hne,Qxne,Qyne      
      TYPE(edge_ptr_array), ALLOCATABLE, DIMENSION(:,:) :: EVi,EVe
!      TYPE(edge_ptr_array), ALLOCATABLE, DIMENSION(:,:) :: rHi,rHe

      TYPE :: edge_ptr_array2
        REAL(pres), POINTER :: in
        REAL(pres), POINTER :: ex
      END TYPE edge_ptr_array2

      TYPE(edge_ptr_array2), ALLOCATABLE, DIMENSION(:,:) :: Hp
      TYPE(edge_ptr_array2), ALLOCATABLE, DIMENSION(:,:) :: Qxp
      TYPE(edge_ptr_array2), ALLOCATABLE, DIMENSION(:,:) :: Qyp
      TYPE(edge_ptr_array2), ALLOCATABLE, DIMENSION(:,:) :: xmp
      TYPE(edge_ptr_array2), ALLOCATABLE, DIMENSION(:,:) :: ymp
      TYPE(edge_ptr_array2), ALLOCATABLE, DIMENSION(:,:) :: xymp
      
      TYPE :: edge_ptr_array3
        REAL(pres), DIMENSION(:),POINTER :: ptr
      END TYPE edge_ptr_array3  
      
      TYPE(edge_ptr_array3), ALLOCATABLE, DIMENSION(:) :: Hai,Hae
      TYPE(edge_ptr_array3), ALLOCATABLE, DIMENSION(:) :: Qxai,Qxae
      TYPE(edge_ptr_array3), ALLOCATABLE, DIMENSION(:) :: Qyai,Qyae
      TYPE(edge_ptr_array3), ALLOCATABLE, DIMENSION(:) :: xmai,xmae
      TYPE(edge_ptr_array3), ALLOCATABLE, DIMENSION(:) :: ymai,ymae
      TYPE(edge_ptr_array3), ALLOCATABLE, DIMENSION(:) :: xymai,xymae
      TYPE(edge_ptr_array3), ALLOCATABLE, DIMENSION(:) :: Hfai,Qxfai,Qyfai
      TYPE(edge_ptr_array3), ALLOCATABLE, DIMENSION(:) :: Hfae,Qxfae,Qyfae      

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
!DIR$ ATTRIBUTES ALIGN:32 :: Hin,Hex,Qxin,Qxex,Qyin,Qyex,xmin,xmex,xymin,xymex,ymin,ymex
!DIR$ ATTRIBUTES ALIGN:32 :: Hfe,Hfi,Qxfi,Qxfe,Qyfi,Qyfe
!DIR$ ATTRIBUTES ALIGN:32 :: Hn,Qxn,Qyn,egnval
!DIR$ ATTRIBUTES ALIGN:32 :: nxv,nyv
!DIR$ ATTRIBUTES ALIGN:32 :: Hni,Hne,Qxni,Qxne,Qyni,Qyne

      END MODULE globals

