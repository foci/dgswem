!******************************************************************************
!     last changes in this file VERSION 46.00
!     
!     Written for ADCIRC v46.00 by Jason G. Fleming.
!******************************************************************************
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     M O D U L E   N O D A L  A T T R I B U T E S
!-----------------------------------------------------------------------
!
!     jgf46.00 This module manages nodal attribute data, including
!     bottom friction, tau0, startdry, directional wind speed reduction,
!     and etc. Will read the Nodal Attributes File (unit 13) and
!     initialize the nodal attribute arrays.
!
!     Handling data by label rather than an integer encoding should
!     result in increased transparency as well as ease the transition to
!     HDF5/NetCDF i/o. The labels were chosen according to the
!     guidelines of the CF Standard. Creating labels according to CF
!     Standard Guidelines should enhance interoperability with other
!     simulation frameworks.
!
!     To use a nodal attribute contained in the fort.13 file, the
!     corresponding attribute name must appear in the fort.15 file.  A
!     list of nodal attributes is read in from the fort.15 file if the
!     fort.15 parameter NWP > 0. This also signals ADCIRC to look for a
!     fort.13 file.
!
!     Summary of the file format for the Nodal Attributes File:
!
!     AGRID                       ! user's comment line - should be a cross
!                                 !    reference to the grid file
!     NumOfNodes                  ! number of nodes, must match NP 
!                                 !    from grid file
!     NAttr                       ! number of attributes contained in this file
!
!     do i=1, NAttr
!        AttrName(i)              ! nodal attribute name (see 
!                                 !    valid names below)
!        Units(i)                 ! physical units (ft, m/s, none)
!        ValuesPerNode(i)         ! number of values at each node for 
!                                 !   a particular attribute
!        DefaultVal(i)            ! default value(s) for the nodal attribute
!     end do
!
!     do i=1, NAttr
!        AttrName(i)              ! label of the attribute, again
!        NumNodesNotDefaultVal(i) ! number of nodes with non-default values
!        do j=1, NumNodesNotDefault(i)
!           n, (Attr(n,k), k=1,ValuesPerNode(i))
!        end do
!     end do
!     
! 
!
!     Valid labels are as follows:
!     
!     ADCIRC Variable:       CF-Style Label:
!      Tau0                  "primitive_weighting_in_continuity_equation"
!      StartDry              "surface_submergence_state"
!      Fric                  "quadratic_friction_coefficient_at_sea_floor"
!      z0Land                "surface_directional_effective_roughness_length"
!                            (z0Land has ValuesPerNode = 12)
!      VCanopy               "surface_canopy_coefficient"
!      BK,BAlpha,BDelX,POAN  "bridge_pilings_friction_parameters"
!                            (bridge_pilings... has ValuesPerNode=4)
!      ManningsN             "mannings_n_at_sea_floor"
!      Chezy                 "chezy_friction_coefficient_at_sea_floor"
!      GeoidOffset           "sea_surface_height_above_geoid"
!      EVM        "average_horizontal_eddy_viscosity_in_sea_water_wrt_depth"
!      EVC        "average_horizontal_eddy_diffusivity_in_sea_water_wrt_depth"
!
!-----------------------------------------------------------------------
      MODULE NodalAttributes
      USE SIZES
      USE GLOBAL
!
!asey 101118: Allow SWAN to handle wave refraction as a nodal attribute.
!             I've placed these changes outside the #ifdef SWAN flags
!             because we want to be able to use the same fort.13 files
!             for both ADCIRC and SWAN+ADCIRC runs.  This way, the new
!             nodal attribute will be processed but only applied when
!             ADCIRC is coupled to SWAN.
      LOGICAL              :: LoadSwanWaveRefrac
      LOGICAL              :: FoundSwanWaveRefrac
      CHARACTER(LEN=80)    :: SwanWaveRefracUnits
      INTEGER              :: SwanWaveRefracNoOfVals
      REAL(SZ)             :: SwanWaveRefracDefVal
      REAL(SZ),ALLOCATABLE :: SwanWaveRefrac(:)
!
!     The following flags are .true. if the corresponding data are
!     required for the run, according to the unit 15 control file
      LOGICAL LoadTau0          
      LOGICAL LoadStartDry      
      LOGICAL LoadDirEffRLen
      LOGICAL LoadCanopyCoef
      LOGICAL LoadQuadraticFric 
      LOGICAL LoadBridgePilings 
      LOGICAL LoadChezy     
      LOGICAL LoadManningsN 
      LOGICAL LoadGeoidOffset
      LOGICAL LoadEVM
      LOGICAL LoadEVC
!
!     The following flags are .true. if there are data with the
!     corresponding label in the unit 13 file.
      LOGICAL FoundTau0          
      LOGICAL FoundStartDry      
      LOGICAL FoundDirEffRLen
      LOGICAL FoundCanopyCoef
      LOGICAL FoundQuadraticFric
      LOGICAL FoundBridgePilings  
      LOGICAL FoundChezy     
      LOGICAL FoundManningsN 
      LOGICAL FoundGeoidOffset
      LOGICAL FoundEVM
      LOGICAL FoundEVC
!
!     These variables hold the strings which describe the attribute's
!     units. These data are loaded from the file, but not used as of
!     v46.00.
      CHARACTER(len=80) Tau0Units   
      CHARACTER(len=80) StartDryUnits      
      CHARACTER(len=80) DirEffRLenUnits  
      CHARACTER(len=80) CanopyCoefUnits
      CHARACTER(len=80) QuadraticFricUnits 
      CHARACTER(len=80) BridgePilingsUnits 
      CHARACTER(len=80) ChezyUnits     
      CHARACTER(len=80) ManningsNUnits 
      CHARACTER(len=80) GeoidOffsetUnits 
      CHARACTER(len=80) EVMUnits 
      CHARACTER(len=80) EVCUnits 
!
!     These variables hold the number of values per node for each
!     attribute.
      INTEGER Tau0NoOfVals
      INTEGER StartDryNoOfVals      
      INTEGER DirEffRLenNoOfVals  
      INTEGER CanopyCoefNoOfVals
      INTEGER QuadraticFricNoOfVals
      INTEGER BridgePilingsNoOfVals
      INTEGER ChezyNoOfVals     
      INTEGER ManningsNNoOfVals 
      INTEGER GeoidOffsetNoOfVals
      INTEGER EVMNoOfVals
      INTEGER EVCNoOfVals
!
!     These variables hold the default values for each attribute.
      REAL(SZ) Tau0DefVal
      REAL(SZ) StartDryDefVal      
      REAL(SZ) DirEffRLenDefVal(12)  
      REAL(SZ) CanopyCoefDefVal
      REAL(SZ) QuadraticFricDefVal
      REAL(SZ) BridgePilingsDefVal(4)
      REAL(SZ) ChezyDefVal     
      REAL(SZ) ManningsNDefVal 
      REAL(SZ) GeoidOffsetDefVal
      REAL(SZ) EVMDefVal
      REAL(SZ) EVCDefVal
!
      INTEGER NumOfNodes    ! number of nodes listed in unit 13 file, cf. NP
      INTEGER NAttr         ! number of nodal attributes in the unit 13 file
!
!     The following variables are inputs from the unit 15 model param. file
      INTEGER NWP     ! number of nodal attributes to read from file
      INTEGER NoLiBF  ! nonlinear bottom friction indicator
      REAL(SZ) Tau0   ! primitive continuity eqn. weight
      REAL(SZ) Tau    ! linear friction coefficient (1/sec)
      REAL(SZ) CF     ! 2DDI bottom fric. coef., effect varies based on NoLiBF
      REAL(SZ) HBreak ! break depth for NOLIBF .eq. 2
      REAL(SZ) FTheta ! dimless param. for NOLIBF .eq. 2
      REAL(SZ) FGamma ! dimless param. for NOLIBF .eq.
      REAL(SZ) ESLM   ! horizontal eddy viscosity (length^2/time)
      REAL(SZ) ESLC   ! horizontal eddy diffusivity (length^2/time)
      INTEGER  IFLINBF! flag to turn on linear bottom friction
      INTEGER  IFNLBF ! flag to turn on nonlinear bottom friction
      INTEGER  IFHYBF ! flag to turn on hybrid bottom friction
!
!     Nodal attributes.
      REAL(SZ), ALLOCATABLE :: STARTDRY(:) ! 1=nodes below geoid initially dry
      REAL(SZ), ALLOCATABLE :: FRIC(:)     ! bottom friction coefficient
      REAL(SZ), ALLOCATABLE :: TAU0VAR(:)  ! primitive equation weighting
      REAL(SZ), ALLOCATABLE :: TAU0BASE(:) ! base (original) primitive equation weighting.
                                           !  Tau0Var may be optimized later based on 
                                           !  Tau0Base and other info if Tau0 value 
                                           !  in Unit 15 is -3. jjw&sb46.39.sb01
      REAL(SZ), ALLOCATABLE :: z0land(:,:) ! directional wind speed red. fac.
      REAL(SZ), ALLOCATABLE :: vcanopy(:)  ! canopy coefficient
!     The following attribute contains BK(I),BALPHA(I),BDELX(I), and POAN(I)
      REAL(SZ), ALLOCATABLE :: BridgePilings(:,:)
      REAL(SZ), ALLOCATABLE :: Chezy(:)
      REAL(SZ), ALLOCATABLE :: ManningsN(:)
      REAL(SZ), ALLOCATABLE :: GeoidOffset(:)
      REAL(SZ), ALLOCATABLE :: EVM(:)
      REAL(SZ), ALLOCATABLE :: EVC(:)
!
      INTEGER i        ! node loop counter
      INTEGER j        ! attribute values loop counter
      INTEGER k        ! attribute loop counter
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      CONTAINS !- - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

!     ----------------------------------------------------------------
!     S U B R O U T I N E     I N I T  N  A  M O D U L E
!     ----------------------------------------------------------------
!
!     jgf46.00 Subroutine to initialize the variables in the nodal
!     attributes module.
!
!     ----------------------------------------------------------------
      SUBROUTINE InitNAModule()
      IMPLICIT NONE
!
!asey 101118: Make changes compact.
      LoadSwanWaveRefrac     = .FALSE.
      FoundSwanWaveRefrac    = .FALSE.
      SwanWaveRefracNoOfVals = 1
      SwanWaveRefracDefVal   = 0.D0
!
      LoadTau0           = .False.
      LoadStartDry       = .False. 
      LoadDirEffRLen     = .False.
      LoadManningsN      = .False. 
      LoadQuadraticFric  = .False.
      LoadChezy          = .False.
      LoadBridgePilings  = .False.
      LoadCanopyCoef     = .False.
      LoadGeoidOffset    = .False.
      LoadEVM            = .False.
      LoadEVC            = .False.
!
      FoundTau0           = .False.
      FoundStartDry       = .False.
      FoundDirEffRLen     = .False.
      FoundManningsN      = .False. 
      FoundQuadraticFric  = .False.
      FoundChezy          = .False.
      FoundBridgePilings  = .False.
      FoundCanopyCoef     = .False.
      FoundGeoidOffset    = .False.
      FoundEVM            = .False.
      FoundEVC            = .False.
!      
      Tau0NoOfVals          = 1
      StartDryNoOfVals      = 1
      DirEffRLenNoOfVals    = 12
      QuadraticFricNoOfVals = 1
      ChezyNoOfVals         = 1
      ManningsNNoOfVals     = 1
      BridgePilingsNoOfVals = 4
      CanopyCoefNoOfVals    = 1
      GeoidOffsetNoOfVals   = 1
      EVMNoOfVals           = 1
      EVCNoOfVals           = 1
!
      Tau0DefVal             = 0.0
      StartDryDefVal         = 0.0
      DO j=1, DirEffRLenNoOfVals
         DirEffRLenDefVal(j) = 0.0
      END DO
      CanopyCoefDefVal       = 0.0
      QuadraticFricDefVal    = 0.0
      DO j=1, BridgePilingsNoOfVals
         BridgePilingsDefVal(j) = 0.0
      END DO
      ChezyDefVal            = 0.0
      ManningsNDefVal        = 0.0
      GeoidOffsetDefVal      = 0.0
      EVMDefVal              = 0.0
      EVCDefVal              = 0.0
!
      HBREAK=1.d0
      FTHETA=1.d0
      FGAMMA=1.d0
!
      ESLM=0.0
      ESLC=0.0
!
      RETURN
!     ----------------------------------------------------------------
      END SUBROUTINE InitNAModule
!     ----------------------------------------------------------------


!     ----------------------------------------------------------------
!     S U B R O U T I N E     R E A D  N O D A L  A T T R
!     ----------------------------------------------------------------
!
!     jgf46.00 Subroutine to read the nodal attributes file (unit 13).
!
!     ----------------------------------------------------------------
      SUBROUTINE ReadNodalAttr(NScreen, ScreenUnit, MyProc, NAbOut)
      IMPLICIT NONE
      INTEGER, intent(in) :: NScreen ! nonzero for debug info to screen
      INTEGER, intent(in) :: ScreenUnit ! i/o for screen 
      INTEGER, intent(in) :: MyProc  ! in parallel, only MyProc=0 i/o to screen
      INTEGER, intent(in) :: NAbOut  ! 1 to abbrev. output to unit 16
!     
      LOGICAL NAFound  ! .true. if Nodal Attributes File (fort.13) exists
      INTEGER ErrorIO  ! zero if file opened successfully
      CHARACTER(len=80) AttrName ! string where the attribute name is stored
      CHARACTER(len=80) header   ! string where alphanumeric file id is stored
      INTEGER NumNodesNotDefault ! number of individual nodes to specify
      LOGICAL SkipDataSet ! .true. if a data set in unit 13 is not needed
      CHARACTER(len=80) Skipped ! data in unit 13 we do not need
      INTEGER L                 ! line counter
!
      NAFound = .False.
      SkipDataSet = .False.
!
!     Check to make sure that NWP is a valid number.
      IF (NWP.LT.0) THEN
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) THEN
            WRITE(ScreenUnit,9972)
            WRITE(ScreenUnit,*) 'NWP =',NWP
            WRITE(ScreenUnit,9728)
            WRITE(ScreenUnit,9973)
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'NWP =',NWP
         WRITE(16,9728)
         WRITE(16,9973)
 9728    FORMAT(/,1X,'Your selection of NWP (a UNIT 15 input ',&
             'parameter) is not an allowable value')
         STOP            ! We're toast.
      ENDIF
!
!     Check to see if there are nodal attributes to be read in. If not,
!     simply return.
      IF (NWP.EQ.0) THEN
         WRITE(16,231) NWP
 231     FORMAT(/,5X,'NWP = ',I2,&
            /,9X,'A Nodal Attributes File (unit 13)',&
             /,9X,'will not be used.')
         RETURN
      ENDIF
!
!     Otherwise, get on with it.
      WRITE(16,232) NWP
 232  FORMAT(/,5X,'NWP = ',I2,&
          /,9X,'Must read Nodal Attributes File (unit 13).')
!
!     Determine if the Nodal Attributes File exists.
      INQUIRE(FILE=TRIM(INPUTDIR)//'/'//'fort.13',EXIST=NAFound)
!
      IF (.not.NAFound) THEN
         WRITE(16,1001)         ! Nodal Attributes file 
         WRITE(16,1011)         ! was not found.
         WRITE(16,9973)         ! execution terminated
         IF (NScreen.ne.0.and.MyProc.eq.0) THEN
            WRITE(ScreenUnit,1001)       
            WRITE(ScreenUnit,1011)       
            WRITE(ScreenUnit,9973)      ! execution terminated
         ENDIF
         STOP
      ENDIF
!
!     Read the unit 15 control file to determine what data must be
!     loaded from nodal attributes file.
      WRITE(16,235) NWP
 235  FORMAT(/,9X,'Need to load ',I2,' nodal attribute(s):')
      DO k=1,NWP
         READ(15,'(A80)') AttrName
         WRITE(16,'(14X,A80)') AttrName
         SELECT CASE (TRIM(ADJUSTL(AttrName)))
         CASE("primitive_weighting_in_continuity_equation")
            LoadTau0 = .True.
         CASE("surface_submergence_state")
            LoadStartDry = .True.
         CASE("quadratic_friction_coefficient_at_sea_floor")
            LoadQuadraticFric = .True.
         CASE("surface_directional_effective_roughness_length") 
            LoadDirEffRLen = .True. 
         CASE("surface_canopy_coefficient") 
            LoadCanopyCoef = .True.
         CASE("bridge_pilings_friction_parameters") 
            LoadBridgePilings = .True.
         CASE("mannings_n_at_sea_floor")
            LoadManningsN = .True. 
         CASE("chezy_friction_coefficient_at_sea_floor") 
            LoadChezy = .True.
         CASE("sea_surface_height_above_geoid") 
            LoadGeoidOffset = .True.
         CASE&
        ("average_horizontal_eddy_viscosity_in_sea_water_wrt_depth")
            LoadEVM = .True.
         CASE&
        ("average_horizontal_eddy_diffusivity_in_sea_water_wrt_depth")
            LoadEVC = .True.
!asey 101118: Allow SWAN to handle wave refraction as a nodal attribute.
         CASE("wave_refraction_in_swan")
            LoadSwanWaveRefrac = .TRUE.
         CASE DEFAULT     
            WRITE(16,1000)          ! unit 15 Model Parameter file
            WRITE(16,1021) AttrName ! contains invalid name
            IF (NScreen.ne.0.and.MyProc.eq.0) THEN
               WRITE(ScreenUnit,1000)       
               WRITE(ScreenUnit,1021) AttrName 
            ENDIF
         END SELECT
      ENDDO
!
!     Now open the nodal attributes (unit 13) file.
      WRITE(16,240) 
 240  FORMAT(/,9X,'Nodal Attributes File (unit 13) was found.',&
          ' Opening file.') 
      OPEN(UNIT=13, FILE=DIRNAME//'/'//'fort.13', &
         IOSTAT=ErrorIO)
      IF ( ErrorIO .GT. 0 ) THEN 
         WRITE(16,1001)         ! Nodal attribute file
         WRITE(16,1005)         ! exists but can't be opened
         WRITE(16,9973)         ! execution terminated
         IF (NScreen.ne.0.and.MyProc.eq.0) THEN
            WRITE(ScreenUnit,1001) 
            WRITE(ScreenUnit,1005)
            WRITE(ScreenUnit,9973) 
         ENDIF
         STOP                   ! We're toast.
      ENDIF
!
!     Read each attribute name, units, number of values, and default value
      READ(13,'(A80)') header
      WRITE(16,250) 
 250  FORMAT(/,9X,'User comment line from unit 13:') 
      WRITE(16,'(14X,A80,/)') header
      READ(13,*) NumOfNodes     ! number of nodes according to unit 13
      READ(13,*) NAttr          ! number of attributes in the unit 13 file
      DO k=1, NAttr
         READ(13,'(A80)') AttrName
         WRITE(16,'(9X,A80)') AttrName
         WRITE(16,260) 
 260     FORMAT(14X,'was found!',/) 
         SELECT CASE (TRIM(ADJUSTL(AttrName)))
         CASE("primitive_weighting_in_continuity_equation")
            FoundTau0 = .True.
            READ(13,'(A80)') Tau0Units   
            READ(13,*) Tau0NoOfVals
            READ(13,*) Tau0DefVal
         CASE("surface_submergence_state")
            FoundStartDry = .True. 
            READ(13,'(A80)') StartDryUnits      
            READ(13,*) StartDryNoOfVals
            READ(13,*) StartDryDefVal            
         CASE("quadratic_friction_coefficient_at_sea_floor")
            FoundQuadraticFric = .True.
            READ(13,'(A80)') QuadraticFricUnits 
            READ(13,*) QuadraticFricNoOfVals
            READ(13,*) QuadraticFricDefVal
         CASE("surface_directional_effective_roughness_length") 
            FoundDirEffRLen = .True. 
            READ(13,'(A80)') DirEffRLenUnits  
            READ(13,*) DirEffRLenNoOfVals
            READ(13,*) (DirEffRLenDefVal(j),j=1,12)    
         CASE("surface_canopy_coefficient") 
            FoundCanopyCoef = .True. 
            READ(13,'(A80)') CanopyCoefUnits  
            READ(13,*) CanopyCoefNoOfVals
            READ(13,*) CanopyCoefDefVal 
         CASE("bridge_pilings_friction_parameters") 
            FoundBridgePilings = .True.
            READ(13,'(A80)') BridgePilingsUnits  
            READ(13,*) BridgePilingsNoOfVals
            READ(13,*) (BridgePilingsDefVal(j),j=1,12)  
         CASE("mannings_n_at_sea_floor")
            FoundManningsN = .True. 
            READ(13,'(A80)') ManningsNUnits 
            READ(13,*) ManningsNNoOfVals
            READ(13,*) ManningsNDefVal 
         CASE("chezy_friction_coefficient_at_sea_floor")
            FoundChezy = .True.
            READ(13,'(A80)') ChezyUnits     
            READ(13,*) ChezyNoOfVals
            READ(13,*) ChezyDefVal          
         CASE("sea_surface_height_above_geoid")
            FoundGeoidOffset = .True.
            READ(13,'(A80)') GeoidOffsetUnits     
            READ(13,*) GeoidOffsetNoOfVals
            READ(13,*) GeoidOffsetDefVal
         CASE&
         ("average_horizontal_eddy_viscosity_in_sea_water_wrt_depth")
            FoundEVM = .True.
            READ(13,'(A80)') EVMUnits     
            READ(13,*) EVMNoOfVals
            READ(13,*) EVMDefVal
         CASE&
        ("average_horizontal_eddy_diffusivity_in_sea_water_wrt_depth")
            READ(13,'(A80)') EVCUnits     
            READ(13,*) EVCNoOfVals
            READ(13,*) EVCDefVal
!asey 101118: Allow SWAN to handle wave refraction as a nodal attribute.
         CASE("wave_refraction_in_swan")
            FoundSwanWaveRefrac = .TRUE.
            READ(13,'(A80)') SwanWaveRefracUnits
            READ(13,*) SwanWaveRefracNoOfVals
            READ(13,*) SwanWaveRefracDefVal
         CASE DEFAULT     
            WRITE(16,1001)          ! Nodal Attributes file
            WRITE(16,1021) AttrName ! contains invalid name
            IF (NScreen.ne.0.and.MyProc.eq.0) THEN
               WRITE(ScreenUnit,1001)       
               WRITE(ScreenUnit,1021) AttrName 
            ENDIF
            READ(13,'(A80)') Skipped  ! skip the Units for the invalid name
            READ(13,'(A80)') Skipped  ! skip the NoOfVals for invalid name
         END SELECT
      END DO
!
!     Determine if there are any attributes required by the fort.15 file
!     that are not in the nodal attributes file.
      IF(((LoadTau0).and.(.not.FoundTau0)).or.&
        ((LoadStartDry).and.(.not.FoundStartDry)).or.&
        ((LoadQuadraticFric).and.&
         (.not.FoundQuadraticFric)).or.&
        ((LoadDirEffRLen).and.&
         (.not.FoundDirEffRLen)).or.&
        ((LoadCanopyCoef).and.&
         (.not.FoundCanopyCoef)).or.&
        ((LoadBridgePilings).and.&
         (.not.FoundBridgePilings)).or.&
        ((LoadManningsN).and.&
         (.not.FoundManningsN)).or.&
        ((LoadGeoidOffset).and.&
         (.not.FoundGeoidOffset)).or.&
        ((LoadChezy).and.(.not.FoundChezy)).or.&
        ((LoadEVM).and.(.not.FoundEVM)).or.&
!asey 101118: Allow SWAN to handle wave refraction as a nodal attribute.&
        ((LoadSwanWaveRefrac).and.(.not.FoundSwanWaveRefrac)).or.&
        ((LoadEVC).and.(.not.FoundEVC)) ) THEN
         WRITE(16,1111)
 1111    FORMAT('ERROR: Nodal Attributes file (unit 13) does '&
             'not contain all the attributes listed in the '&
             /,'model parameter file (unit 15).')
         WRITE(16,9973)         ! execution terminated
         IF (NScreen.ne.0.and.MyProc.eq.0) THEN
            WRITE(ScreenUnit,1111) 
            WRITE(ScreenUnit,9973)       ! execution terminated
         ENDIF
         STOP                   ! We're toast.
      ENDIF

!     Allocate memory to hold our data.
      ALLOCATE(TAU0VAR(NumOfNodes),TAU0BASE(NumOfNodes)) ! jjw&sb46.39sb01
      ALLOCATE(STARTDRY(NumOfNodes)) 
      ALLOCATE(FRIC(NumOfNodes))
      ALLOCATE(z0land(NumOfNodes,DirEffRLenNoOfVals))
      ALLOCATE(vcanopy(NumOfNodes))
      ALLOCATE(BridgePilings(NumOfNodes,BridgePilingsNoOfVals))
      ALLOCATE(GeoidOffset(NumOfNodes))
      ALLOCATE(Chezy(NumOfNodes))
      ALLOCATE(ManningsN(NumOfNodes))
      ALLOCATE(EVM(NumOfNodes))
      ALLOCATE(EVC(NumOfNodes))
!asey 101118: Allow SWAN to handle wave refraction as a nodal attribute.
      ALLOCATE(SwanWaveRefrac(NumOfNodes))
!
!     Now read each of the attributes required by the model parameter
!     (unit 15) file and skip past the others.
      WRITE(16,270) NWP
 270  FORMAT(/,9X,'Now reading ',I2,' nodal attribute(s).')
      DO k=1, NAttr
         WRITE(16,280) k
 280     FORMAT(/,9X,'Attribute ',I2,':')
         READ(13,'(A80)') AttrName
         READ(13,*) NumNodesNotDefault
         WRITE(16,'(14X,A80)') AttrName
         SELECT CASE (TRIM(ADJUSTL(AttrName)))
         CASE("primitive_weighting_in_continuity_equation")
            IF (LoadTau0) THEN
               CALL LoadAttrVec(TAU0VAR, Tau0DefVal,&
                   NumNodesNotDefault, NScreen, MyProc, NAbOut)
            ELSE
               SkipDataSet = .True.
            ENDIF
         CASE("surface_submergence_state")
            IF (LoadStartDry) THEN
               CALL LoadAttrVec(STARTDRY, StartDryDefVal, &
                    NumNodesNotDefault, NScreen, MyProc, NAbOut)
            ELSE
               SkipDataSet = .True.
            ENDIF
         CASE("quadratic_friction_coefficient_at_sea_floor")
            IF (LoadQuadraticFric) THEN
               CALL LoadAttrVec(FRIC, QuadraticFricDefVal, &
                   NumNodesNotDefault, NScreen, MyProc, NAbOut)
            ELSE
               SkipDataSet = .True.
            ENDIF
         CASE("surface_directional_effective_roughness_length") 
            IF (LoadDirEffRLen) THEN
               CALL LoadAttrMat(z0land, DirEffRLenNoOfVals,  &
                   DirEffRLenDefVal, NumNodesNotDefault,  &
                   NScreen, MyProc, NAbOut)
            ELSE
               SkipDataSet = .True.
            ENDIF
         CASE("surface_canopy_coefficient") 
            IF (LoadCanopyCoef) THEN
               CALL LoadAttrVec(vcanopy, CanopyCoefDefVal,&
                    NumNodesNotDefault, NScreen, MyProc, NAbOut)
            ELSE
               SkipDataSet = .True.
            ENDIF
         CASE("bridge_pilings_friction_parameters") 
            IF (LoadBridgePilings) THEN
               CALL LoadAttrMat(BridgePilings, BridgePilingsNoOfVals, &
                   BridgePilingsDefVal,  NumNodesNotDefault, &
                   NScreen, MyProc, NAbOut)
            ELSE
               SkipDataSet = .True.
            ENDIF
         CASE("mannings_n_at_sea_floor")
            IF (LoadManningsN) THEN
               CALL LoadAttrVec(ManningsN, ManningsNDefVal, &
                   NumNodesNotDefault, NScreen, MyProc, NAbOut)
            ELSE
               SkipDataSet = .True.
            ENDIF
         CASE("chezy_friction_coefficient_at_sea_floor")
            IF (LoadChezy) THEN
               CALL LoadAttrVec(Chezy, ChezyDefVal, &
                    NumNodesNotDefault, NScreen, MyProc, NAbOut)
            ELSE
               SkipDataSet = .True.
            ENDIF
         CASE("sea_surface_height_above_geoid")
            IF (LoadGeoidOffset) THEN
               CALL LoadAttrVec(GeoidOffset, GeoidOffsetDefVal, &
                    NumNodesNotDefault, NScreen, MyProc, NAbOut)
            ELSE
               SkipDataSet = .True.
            ENDIF
         CASE&
        ("average_horizontal_eddy_viscosity_in_sea_water_wrt_depth")
            IF (LoadEVM) THEN
               CALL LoadAttrVec(EVM, EVMDefVal, &
                    NumNodesNotDefault, NScreen, MyProc, NAbOut)
            ELSE
               SkipDataSet = .True.
            ENDIF
         CASE&
        ("average_horizontal_eddy_diffusivity_in_sea_water_wrt_depth")
            IF (LoadEVC) THEN
               CALL LoadAttrVec(EVC, EVCDefVal, &
                    NumNodesNotDefault, NScreen, MyProc, NAbOut)
            ELSE
               SkipDataSet = .True.
            ENDIF
!asey 101118: Allow SWAN to handle wave refraction as a nodal attribute.
         CASE("wave_refraction_in_swan")
            IF (LoadSwanWaveRefrac) THEN
               CALL LoadAttrVec(SwanWaveRefrac, SwanWaveRefracDefVal,&
                   NumNodesNotDefault, NScreen, MyProc, NAbOut)
            ELSE
               SkipDataSet = .TRUE.
            ENDIF
         CASE DEFAULT
            SkipDataSet = .True.
            WRITE(16,1001)      ! Nodal Attributes file
            WRITE(16,1021) AttrName ! contains invalid name
            IF (NScreen.ne.0.and.MyProc.eq.0) THEN
               WRITE(ScreenUnit,1001)       
               WRITE(ScreenUnit,1021) AttrName 
            ENDIF
         END SELECT
         IF (SkipDataSet) THEN
            DO L=1, NumNodesNotDefault
               READ(13,*) Skipped
            END DO
            WRITE(16,'(9X,A8)') 'Skipped.'
            SkipDataSet = .False.
         ELSE
            WRITE(16,'(/,9X,A18,A80)') 'Finished loading ', AttrName
         ENDIF
      END DO
!
 1000 FORMAT('ERROR: The Model Parameter File (unit 15)')
 1001 FORMAT('ERROR: The Nodal Attributes File (unit 13)')
 1002 FORMAT('ERROR: The legacy StartDry File (unit 12)')
 1003 FORMAT('ERROR: Spatially Varying Fric. Coeff. File (unit 21)')
!
 1005 FORMAT('exists but cannot be opened.')
 1011 FORMAT('was not found.') 
 1021 FORMAT('contains invalid name: ',A80) 
 9972 FORMAT(////,1X,'!!!!!!!!!! INPUT ERROR !!!!!!!!!',/)
 9973 FORMAT(/,1X,'!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!',//)
!      
      RETURN
!     ----------------------------------------------------------------
      END SUBROUTINE ReadNodalAttr
!     ----------------------------------------------------------------


!     ----------------------------------------------------------------
!         S U B R O U T I N E     L O A D  A T T R  V E C 
!     ----------------------------------------------------------------
!
!     jgf46.00 Subroutine to set a single set of nodal attributes to
!     their user-specified default values, then read the nondefault
!     values from the Nodal Attributes File (unit 13). This subroutine
!     is used for nodal attributes with only one value per node, hence
!     the suffix "vec" in the name.
!
!     ----------------------------------------------------------------
      SUBROUTINE LoadAttrVec(AttributeData, Default, NumNodesNotDef,NScreen, MyProc, NAbOut)
      IMPLICIT NONE
      REAL(SZ), intent(out), dimension(NumOfNodes) :: AttributeData
      REAL(SZ), intent(in):: Default ! default value for all nodes
      INTEGER, intent(in) :: NumNodesNotDef ! number of nodes specified in file
      INTEGER, intent(in) :: NScreen ! 1 for debug info to screen (unit 6)
      INTEGER, intent(in) :: MyProc  ! in parallel, only MyProc=0 i/o to screen
      INTEGER, intent(in) :: NAbOut  ! 1 to abbrev. output to unit 16
!
      INTEGER NodeNum            ! node number listed in the file
!
!     Set all values to user-specified default values.
      IF (NABOUT.EQ.0) WRITE(16,1001) Default 
      DO i=1, NumOfNodes
         AttributeData(i) = Default 
      END DO
!
      IF (NABOUT.EQ.0) WRITE(16,1005) 
      DO i=1, NumNodesNotDef
         READ(13,*) NodeNum, AttributeData(NodeNum)
         IF (NABOUT.EQ.0) WRITE(16,1010) NodeNum, AttributeData(NodeNum)     
      END DO
!
 1001 FORMAT(/,10X,'Set all nodes to the default value of ',E16.8,/)
 1005 FORMAT(/,10X,'Now setting the following nodes to these values:',&
          /,10X,'NODE',5X,'DATA',5X/)
 1010 FORMAT(7X,I6,6X,E16.8)
!
      RETURN
!     ----------------------------------------------------------------
      END SUBROUTINE LoadAttrVec
!     ----------------------------------------------------------------


!     ----------------------------------------------------------------
!         S U B R O U T I N E     L O A D  A T T R  M A T
!     ----------------------------------------------------------------
!
!     jgf46.00 Subroutine to load a single set of nodal attributes from
!     the Nodal Attributes File (unit 13) if there is more than one
!     value per node.
!
!     ----------------------------------------------------------------
      SUBROUTINE LoadAttrMat(AttributeData, NumCol, Default,NumNodesNotDef, NScreen, MyProc, NAbOut)
      IMPLICIT NONE
      INTEGER, intent(in) :: NumCol  ! number of columns in the matrix 
      REAL(SZ), intent(out), &
          dimension(NumOfNodes,NumCol) :: AttributeData
      REAL(SZ), intent(in), dimension(NumCol) :: Default ! default values
      INTEGER, intent(in) :: NumNodesNotDef  ! number of nodes spec. in file
      INTEGER, intent(in) :: NScreen ! 1 for debug info to screen (unit 6)
      INTEGER, intent(in) :: MyProc  ! in parallel, only MyProc=0 i/o to screen
      INTEGER, intent(in) :: NAbOut  ! 1 to abbrev. output to unit 16
!
      INTEGER NodeNum            ! node number listed in the file
!
!     Set all nodes to user-specified default values.
      IF (NABOUT.EQ.0) WRITE(16,1001) 
      DO i=1, NumOfNodes
         DO j=1, NumCol
            AttributeData(i,j)=Default(j)
         END DO
      END DO
!
      IF (NABOUT.EQ.0) WRITE(16,1005) 
      DO i=1, NumNodesNotDef
         READ(13,*) NodeNum, (AttributeData(NodeNum,j),j=1,NumCol)
         IF (NABOUT.EQ.0) WRITE(16,1010) NodeNum, &
             (AttributeData(NodeNum,j),j=1,NumCol)
      END DO
!
 1001 FORMAT(/,10X,'Set all nodes to the default values of ',/,&
          99E16.8,/)
 1005 FORMAT(/,10X,'Now setting the following nodes to these values:',&
          /,10X,'NODE',5X,'DATA',5X/)
 1010 FORMAT(7X,I6,6X,12(1X,E16.8))
!
      RETURN
!     ----------------------------------------------------------------
      END SUBROUTINE LoadAttrMat
!     ----------------------------------------------------------------



!     ----------------------------------------------------------------
!         S U B R O U T I N E     I N I T  N O D A L  A T T R 
!     ----------------------------------------------------------------
!
!     jgf46.00 Subroutine to initialize and error check the nodal 
!     attributes read in from the Nodal Attributes File (unit 13).
!
!     ----------------------------------------------------------------
      SUBROUTINE InitNodalAttr(DP, NP, G, NScreen, ScreenUnit,MyProc, NAbOut)
      IMPLICIT NONE
      INTEGER, intent(in) :: NP ! number of nodes in the grid file 
      REAL(SZ), intent(in), dimension(NP) :: DP ! array of bathymetric depths
      REAL(SZ), intent(in):: G  ! gravitational acceleration
      INTEGER, intent(in) :: NScreen ! nonzero for debug info to screen
      INTEGER, intent(in) :: ScreenUnit ! i/o for debug info to screen
      INTEGER, intent(in) :: MyProc  ! in parallel, only MyProc=0 i/o to screen
      INTEGER, intent(in) :: NAbOut  ! 1 to abbrev. output to unit 16
!
!
!     ERROR CHECK: If a nodal attributes file is being used, check to
!     see that the number of nodes in the nodal attribute file is the
!     same as the number of nodes in the grid file.
      IF (NWP.NE.0.AND.NumOfNodes.NE.NP) THEN
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,9900)
         WRITE(16,9900)
 9900    FORMAT(////,1X,'!!!!!!!!!!  FATAL ERROR  !!!!!!!!!',&
             //,1X,'The number of nodes in the grid file (unit 14) and'&
             /,1X,'the nodal attributes file (unit 13) must match.',&
             //,1X,'!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!',//)
         WRITE(16,*) 'np, numofnodes ',np,NumOfNodes
         STOP                   ! We're toast.
      ENDIF
!
!     ERROR CHECK: If Chezy, Manning's or Quadratic friction was loaded
!     from the nodal attributes file, NOLIBF must be  >= 1. 
      IF ((LoadChezy.or.LoadManningsN.or.LoadQuadraticFric).and.&
          NoLiBF.eq.0) THEN
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,9800)
         WRITE(16,9800) nolibf
 9800    FORMAT(////,1X,'!!!!!!!!!!  FATAL ERROR  !!!!!!!!!',&
             //,1X,'Nonlinear bottom friction coefficients were loaded'&
             /,1X,'from the nodal attributes file (unit 13), so ',&
             /,1X,'NoLiBF must be set to 1. It is set to ',i2,' in',&
             /,1X,'the model parameter (unit 15) file.',&
             //,1X,'!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!',//)
         STOP                   ! We're toast.
      ENDIF
!     
!     I N I T    S T A R T D R Y
      IF (.not.LoadStartDry) THEN
!         IF (NWP.eq.0) THEN
!            ALLOCATE(STARTDRY(NP))
!         ENDIF
         DO I=1, NP
            STARTDRY(I) = 0.0D0
         ENDDO
      ENDIF
!
!     I N I T     T A U 0
      IF (.not.LoadTau0) THEN
!         IF (NWP.eq.0) THEN
!            ALLOCATE(TAU0VAR(NP),TAU0BASE(NP)) ! jjw&sb46.39sb01
!         ENDIF
!
!     jgf46.25 If input tau0 is positive, set all nodes to that value.
         IF (Tau0.ge.0) THEN
            DO I=1,NP
               Tau0Var(I)=Tau0
            END DO
            WRITE(16,7) Tau0
 7          FORMAT(/,5X,&
                'A SPATIALLY CONSTANT WEIGHTING COEFFICIENT (Tau0)'&
                ,/,5X,' WILL BE USED IN THE GENERALIZED WAVE',&
                ' CONTINUITY EQUATION.',&
                /,5X,'Tau0 = ',E15.8,2X,'1/sec',/)
         ELSE
            DO I=1,NP
               Tau0Var(I)=Tau0NodalValue(Tau0,DP(I))
            ENDDO
            
            IF(Tau0.eq.-2) THEN
               WRITE(16,6) ! spatially vary tau0 according to hard coded scheme
               WRITE(16,62) ! description of scheme
 62            FORMAT(/,5X,'IF DEPTH > 200           Tau0 = 0.005',    &
                   /,5X,'IF 200   > DEPTH > 1     Tau0 = 1/DEPTH  ',    &
                   /,5X,'IF 1     > DEPTH         Tau0 = 1.0 ')
            ELSE
               WRITE(16,6) ! spatially vary tau0 according to hard coded scheme
               WRITE(16,61) ! description of scheme
 61            FORMAT(/,5X,' IF DEPTH GE 10           -> TAU0 = 0.005',    &
                   /,5X,' IF DEPTH LT 10           -> TAU0 = 0.020',/)
            ENDIF
         ENDIF
      ENDIF
 6    FORMAT(/,5X,'A SPATIALLY VARIABLE WEIGHTING COEFFICIENT (Tau0)'&
          ,/,5X,' WILL BE USED IN THE GENERALIZED WAVE',&
          ' CONTINUITY EQUATION.',&
          /,5x,'THIS VALUE WILL BE DETERMINED AS FOLLOWS:')

!
!     jgf46.27 If we have already loaded the tau0 values directly from
!     the nodal attributes file, check to see if the default value was
!     negative. If so, this indicates that nodal values of tau0 that
!     were not explicitly set in the nodal attributes file should be set
!     accordiing to one of the hard-coded tau0 schemes.

      IF (LoadTau0.and.Tau0DefVal.lt.0) THEN
         DO I=1,NP
            IF (Tau0Var(I).lt.0) THEN
               Tau0Var(I)=Tau0NodalValue(Tau0DefVal,DP(I))
            ENDIF
         ENDDO
      ENDIF
!
!     jjw&sb46.39.sb01 Save the Tau0Var in Tau0Base for later use
      DO I=1,NP
        Tau0Base(I) = Tau0Var(I)
      ENDDO
!
!     jjw&sb46.38.sb01 If tau0 is loaded from nodal attributes file and 
!     Tau0 is -3, time-varing tau0 optimizer will be applied in timestep.F
      if(LoadTau0.AND.Tau0.eq.-3.d0) then
         WRITE(16,7)
      endif
 8    FORMAT(/,5X,'A SPATIALLY TEMPORALLY VARIABLE OPTIMIZED '&
          ,/,5X,' WEIGHTING COEFFICIENT (Tau0) WILL BE USED '&
          ,/,5X,' IN THE GENERALIZED WAVE CONTINUITY EQUATION.',/)
!
!     I N I T   B O T T O M   F R I C T I O N
      IF(NOLIBF.EQ.0) THEN
         IFNLBF=0
         IFLINBF=1
         IFHYBF=0
      ENDIF
      IF(NOLIBF.EQ.1) THEN
         IFNLBF=1
         IFLINBF=0
         IFHYBF=0
      ENDIF
      IF(NOLIBF.EQ.2) THEN
         IFNLBF=0
         IFLINBF=0
         IFHYBF=1
      ENDIF
!     
!     Initialize bottom friction if it was not loaded from unit 13.
      IF((.not.LoadQuadraticFric).and.(.not.LoadManningsN).and.&
          (.not.LoadChezy)) THEN
         IF (NoLiBF.eq.0) CF=Tau
!     If a nodal attributes file was read, FRIC was allocated there.
!         IF (NWP.eq.0) THEN
!            ALLOCATE(FRIC(NP))
!         ENDIF
         DO I=1,NP
            FRIC(I)=CF
         END DO
      ENDIF

!
!     Initialize bridge pilings.
!
      IF (LoadBridgePilings) THEN
         DO I=1, NP
            IF (BridgePilings(I,1).ne.0) THEN ! only for nodes w/piers
               BridgePilings(I,3) = 4.d0 * &
                   BridgePilings(I,3) / BridgePilings(I,4)
            ENDIF
         END DO
      ENDIF
!
!     I N I T   E D D Y   V I S C O S I T Y  &  D I F F U S I V I T Y
      IF (.not.LoadEVM) THEN
!	 IF (NWP .eq. 0) THEN
!          	ALLOCATE(EVM(NP))
!	 ENDIF
         DO I=1,NP
            EVM(I)=ESLM
         END DO
      ENDIF
      DO I=1,NP
         EVMSUM=EVMSUM+ABS(EVM(I))
      ENDDO
      IF (.not.LoadEVC.and.ESLC.ne.0) THEN
         DO I=1,NP
            EVC(I)=ESLC
         END DO
      ENDIF
!
      RETURN
!     ----------------------------------------------------------------
      END SUBROUTINE InitNodalAttr
!     ----------------------------------------------------------------


!     ----------------------------------------------------------------
!      F U N C T I O N   T A U 0  N O D A L  V A L U E 
!     ----------------------------------------------------------------
!
!     jgf46.27 Function to calculate tau0 based on the scheme selection
!     and the depth. This assumes that Scheme is negative. 
!       
!     ----------------------------------------------------------------
      REAL(SZ) FUNCTION Tau0NodalValue(Scheme, Depth)
      IMPLICIT NONE
      REAL(SZ) Scheme         
      REAL(SZ) Depth
!
      IF (Scheme.eq.-2.d0) THEN 
!     Smoothly varying tau0 with depth.
         IF(Depth.GE.200.) Tau0NodalValue=0.005
         IF((Depth.LT.200.).AND.(Depth.GE.1.)) THEN
            Tau0NodalValue=1./Depth
         ENDIF
         IF(Depth.LT.1.) Tau0NodalValue=1.0
      ELSE
!     Abrupt variation in tau0 with depth.           
         IF(Depth.LE.10.) Tau0NodalValue=0.020d0    
         IF(Depth.GT.10.) Tau0NodalValue=0.005d0
      ENDIF         
!     ----------------------------------------------------------------
      END FUNCTION Tau0NodalValue
!     ----------------------------------------------------------------


!     ----------------------------------------------------------------
!     S U B R O U T I N E     
!     A P P L Y  2 D  B O T T O M  F R I C T I O N
!     ----------------------------------------------------------------
!
!     jgf46.00 Subroutine to apply 2D bottom friction from turbulent
!     viscous effects as well as bridge pilings. This is used in the
!     time stepping loop.
!
!     ----------------------------------------------------------------
!
!     sb46.28sb02 Lower limit of Cd was added as an argument.
! 
!     ----------------------------------------------------------------
      SUBROUTINE Apply2DBottomFriction(UU1, VV1, DP, ETA2, G,&
          IFNLFA, NP, TK, LL)
      USE SIZES
      IMPLICIT NONE
      INTEGER, intent(in) :: NP                   ! number of nodes in grid
      REAL(SZ), intent(in), dimension(NP) :: UU1  ! x-dir velocities
      REAL(SZ), intent(in), dimension(NP) :: VV1  ! y-dir velocities
      REAL(SZ), intent(in), dimension(NP) :: DP   ! bathymetric depths
      REAL(SZ), intent(in), dimension(NP) :: ETA2 ! water surf. elevations
      REAL(SZ), intent(in) :: G                   ! gravitational constant
      INTEGER, intent(in) :: IFNLFA               ! nonlin. finite amp. flag 
      REAL(SZ), intent(inout), dimension(NP) :: TK! depth avg. fric.
      REAL(SZ), intent(in)                :: LL   ! lower limit of Cd sb46.28sb02
!
      REAL(SZ) UV1   ! velocity magnitude (speed)
      REAL(SZ) H1    ! total depth
      REAL(SZ) Fr
      REAL(SZ) FricBP
      REAL(SZ) BK    ! BK(1) is pier shape factor
      REAL(SZ) BALPHA! BALPHA(2) is constriction fraction
      REAL(SZ) BDELX ! BDELX(3) is effective delx
!
!     Step 0. Convert Manning's N to Cd, if necessary.
      IF (LoadManningsN) THEN
         DO I=1, NP
            FRIC(I)=g*ManningsN(I)**2.d0&
                /( ( DP(I)+IFNLFA*ETA2(I) )**(1.d0/3.d0) ) ! sb46.28sb02
            !sb46.28sb02  Lower limit is applied here.
            IF(FRIC(I).LT.LL) THEN
               FRIC(I) = LL
            ENDIF
         ENDDO
      ENDIF
!
!     ... Convert Chezy to Cd, if necessary.
      IF (LoadChezy) THEN
         DO I=1,NP
            FRIC(I)=G*(Chezy(I)**2)
         END DO
      ENDIF
!
!     Step 1. Apply friction arising from turbulent viscous interaction
!     with the sea floor.
      DO I=1, NP
         UV1=SQRT(UU1(I)*UU1(I)+VV1(I)*VV1(I))
         H1=DP(I)+IFNLFA*ETA2(I)
         TK(I)= FRIC(I)*&
             ( IFLINBF +       &! linear&
             (UV1/H1) * (IFNLBF &! nonlinear&
             + IFHYBF*(1+(HBREAK/H1)**FTHETA)**(FGAMMA/FTHETA))) ! hybrid
      END DO
!
!     Step 2. Apply friction arising from flow interaction with bridge
!     pilings, if required.
      IF (LoadBridgePilings) THEN
         DO I=1, NP
            UV1=SQRT(UU1(I)*UU1(I)+VV1(I)*VV1(I))
            H1=DP(I)+IFNLFA*ETA2(I)
            Fr=UV1*UV1/(G*H1)
            BK = BridgePilings(I,1) 
            BALPHA = BridgePilings(I,2) 
            BDELX = BridgePilings(I,3) 
            FricBP=(H1/BDELX)*BK*(BK+5.d0*Fr*Fr-0.6d0)&
                *(BALPHA+15.d0*BALPHA**4)
            TK(I)=TK(I)+FricBP*UV1/H1
         END DO
      ENDIF
!
      RETURN
!     ----------------------------------------------------------------
      END SUBROUTINE Apply2DBottomFriction
!     ----------------------------------------------------------------


!     ----------------------------------------------------------------
!     S U B R O U T I N E     
!     A P P L Y  D I R E C T I O N A L  W I N D  R E D U C T I O N
!     ----------------------------------------------------------------
!
!     jgf46.00 Subroutine to calculate the land wind reduction factor
!     based on a table of directional wind drag values. Originally
!     written into the hstart.F file by jjw in jjw-42.06j. This is used
!     in hstart.F and timestep.F.
!
!     ----------------------------------------------------------------
      SUBROUTINE ApplyDirectionalWindReduction(NodeNumber,WindDragCo,WindMag,BathymetricDepth,Elevation,CutOffDepth,G,WindX,WindY)
      USE SIZES
      IMPLICIT NONE
      INTEGER,  intent(in) :: NodeNumber ! index of node under consideration
      REAL(SZ), intent(in) :: WindDragCo ! wind drag coefficient
      REAL(SZ), intent(in) :: WindMag    ! wind magnitude
      REAL(SZ), intent(in) :: BathymetricDepth ! a.k.a. dp(i),depth below geoid
      REAL(SZ), intent(in) :: Elevation  ! a.k.a. eta2(i)
      REAL(SZ), intent(in) :: CutOffDepth! a.k.a. h0, user-spec. min. depth
      REAL(SZ), intent(in) :: G          ! gravitational constant

      REAL(SZ), intent(inout) :: WindX   ! x-dir component of wind velocity
      REAL(SZ), intent(inout) :: WindY   ! x-dir component of wind velocity

      REAL(SZ) z0m   ! marine roughness coefficient based on Garratt's formula
      REAL(SZ) angle ! direction wind is coming from
      INTEGER idir   ! code for wind direction
      REAL(SZ) z0l   ! drag for a particular node, for particular direction
      REAL(SZ) TotalDepth  ! bathymetric depth + sea surface elevation
      REAL(SZ) fr    ! land wind reduction factor
!
!     compute marine roughness coefficient based on Garratt's formula
      z0m=(0.018d0/G)*WindDragCo*WindMag**2.d0
!
!     compute direction  that the wind is coming from
      if((WindX.eq.0).and.(WindY.eq.0))then
         angle=0.d0
      else
         angle=atan2(WindY,WindX)
      endif
      angle=360.*angle/(2*3.141592654d0)
      idir=0
      if((angle.gt.-15.).and.(angle.le.15))  idir=1
      if((angle.gt.15.).and.(angle.le.45))   idir=2
      if((angle.gt.45.).and.(angle.le.75))   idir=3
      if((angle.gt.75.).and.(angle.le.105))  idir=4
      if((angle.gt.105.).and.(angle.le.135)) idir=5
      if((angle.gt.135.).and.(angle.le.165)) idir=6
      if((angle.gt.165.).and.(angle.le.180)) idir=7
      if((angle.gt.-45.).and.(angle.le.-15)) idir=12
      if((angle.gt.-75.).and.(angle.le.-45)) idir=11
      if((angle.gt.-105.).and.(angle.le.-75)) idir=10
      if((angle.gt.-135.).and.(angle.le.-105)) idir=9
      if((angle.gt.-165.).and.(angle.le.-135)) idir=8
      if((angle.ge.-180.).and.(angle.le.-165)) idir=7
!
!     define land roughness from usace values
      z0l=z0land(NodeNumber,idir)
!
!     reset z0l depending on situation                       
      if(z0l.le.0.006) then    
!     coe set their value to a marine value -> reset to correct marine value
         z0l=z0m    
      else 
!     coe set their value to a land value -> proceed with checking this value
         TotalDepth = BathymetricDepth + Elevation
         if( (TotalDepth.gt.2*CutOffDepth).and.&
             (BathymetricDepth.lt.20)) then
!     compute adjusted z0l to account for overland flooding - do this only 
!     in the case where the water column is greater than twice h0 and 
!     you are not in a river (I assume that rivers are deeper than 20m 
!     and have z0l>0.006) 
            z0l=z0l-TotalDepth/30. ! correction for overland flooding
         endif 
      endif  
!
!     compute land wind reduction factor
      if(z0l.gt.0.0001) then
         fr=(z0m/z0l)**0.0706d0
      else
         fr=1.000d0
      endif
      if(fr.gt.1.0000d0) fr=1.0000d0
!     adjust time interpolated wind field
      WindX = vcanopy(NodeNumber)*fr*WindX 
      WindY = vcanopy(NodeNumber)*fr*WindY 
!
      RETURN
!     ----------------------------------------------------------------
      END SUBROUTINE ApplyDirectionalWindReduction
!     ----------------------------------------------------------------



!     ----------------------------------------------------------------
!                      S U B R O U T I N E
!         R E A D  L E G A C Y  S T A R T  D R Y  F I L E
!     ----------------------------------------------------------------
!
!     jgf46.00 Subroutine to load up the legacy startdry file (unit
!     12). This is just a cut-and-paste from the section of the
!     READ_INPUT subroutine that did the same thing. This subroutine is
!     never called. It is vestigial and listed here purely as reference
!     material.
!
!     ----------------------------------------------------------------
      SUBROUTINE ReadLegacyStartDryFile(NP, NScreen, ScreenUnit, &
          MyProc, NAbOut)
      IMPLICIT NONE
      INTEGER, intent(in) :: NP ! number of nodes in grid file
      INTEGER, intent(in) :: NScreen ! nonzero for debug info to screen
      INTEGER, intent(in) :: ScreenUnit ! i/o for debug info to screen
      INTEGER, intent(in) :: MyProc  ! in parallel, only MyProc=0 i/o to screen
      INTEGER, intent(in) :: NAbOut  ! 1 to abbrev. output to unit 16

      INTEGER JKI          ! node number from file
      INTEGER NE2          ! number of elements, according to fort.12 file
      INTEGER NP2          ! number of nodes, according to fort.12 file

      CHARACTER(len=80) AGRID2 ! users comment/description line
      REAL(SZ) DUM1, DUM2  ! data that we want to skip

      OPEN(12,FILE=TRIM(INPUTDIR)//'/'//'fort.12')
!
!...  READ STARTDRY INFORMATION FROM UNIT 12 
      READ(12,'(A80)') AGRID2
      WRITE(16,2038) AGRID2
2038  FORMAT(5X,'STARTDRY FILE IDENTIFICATION : ',A80,/)
      READ(12,*) NE2,NP2
!
!...  CHECK THAT NE2 AND NP2 MATCH WITH GRID FILE 
!      IF((NE2.NE.NE).OR.(NP2.NE.NP)) THEN
       IF(NP2.NE.NP) THEN
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,9900)
         WRITE(16,9900)
 9900    FORMAT(////,1X,'!!!!!!!!!!  FATAL ERROR  !!!!!!!!!',&
             //,1X,'THE PARAMETER NE2 AND NP2 MUST MATCH NE AND NP ',&
             /,1X,'USER MUST CHECK FORT.12 INPUT FILE ',&
             //,1X,'!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!',//)
         STOP
      ENDIF
!
!...  READ IN STARTDRY CODE VALUES
      DO I=1,NP
         READ(12,*) JKI,DUM1,DUM2,STARTDRY(JKI)
         IF(JKI.NE.I) THEN
            IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,99805)
            WRITE(16,99805)
99805       FORMAT(////,1X,'!!!!!!!!!!  WARNING - NONFATAL ',&
                'INPUT ERROR  !!!!!!!!!',&
                //,1X,'YOUR NODE NUMBERING IS NOT SEQUENTIAL ',&
                'CHECK YOUR UNIT 12 INPUT FILE CAREFULLY',//)
         ENDIF
      END DO
!     
!...  CLOSE UNIT 12 FILE       
      CLOSE(12)     
!
      RETURN
!     ----------------------------------------------------------------
      END SUBROUTINE ReadLegacyStartDryFile
!     ----------------------------------------------------------------


!     ----------------------------------------------------------------
!                       S U B R O U T I N E     
!       R E A D  L E G A C Y  B O T T O M  F R I C T I O N  F I L E
!     ----------------------------------------------------------------
!
!     jgf46.00 Subroutine to load up the legacy Spatially Varying
!     Friction Coefficient File (unit 21). This is just a cut-and-paste
!     from the section of the READ_INPUT subroutine that did the same
!     thing. This subroutine is never called. It is vestigial and listed
!     here purely as reference material.
!
!     ----------------------------------------------------------------
      SUBROUTINE ReadLegacyBottomFrictionFile(NP, NScreen, ScreenUnit,MyProc, NAbOut)
      IMPLICIT NONE
      INTEGER, intent(in) :: NP ! number of nodes in grid file
      INTEGER, intent(in) :: NScreen ! nonzero for debug info to screen
      INTEGER, intent(in) :: ScreenUnit ! i/o for debug info to screen
      INTEGER, intent(in) :: MyProc  ! in parallel, only MyProc=0 i/o to screen
      INTEGER, intent(in) :: NAbOut  ! 1 to abbrev. output to unit 16
      
      CHARACTER(len=80) AFRIC  ! user's comment/description line
      INTEGER NHG    ! node number from file

      OPEN(21,FILE=TRIM(INPUTDIR)//'/'//'fort.21')
      READ(21,'(A80)') AFRIC
      DO I=1,NP
         READ(21,*) NHG,FRIC(NHG)
         IF(NHG.NE.I) THEN
            IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,99803)
            WRITE(16,99803)
99803       FORMAT(////,1X,'!!!!!!!!!!  WARNING - FATAL ',&
                'INPUT ERROR  !!!!!!!!!',//,1X,&
                'YOUR NODAL FRICTION NUMBERING IS NOT SEQUENTIAL ',&
                /,1X,'CHECK YOUR UNIT 21 INPUT FILE CAREFULLY',//,1X,&
                '!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!',//)
            STOP
         ENDIF
      END DO
      WRITE(16,3601) AFRIC
 3601 FORMAT(/,5X,'FRICTION FILE IDENTIFICATN : ',A80,/)
      IF(NABOUT.NE.1) THEN
         WRITE(16,2080)
 2080    FORMAT(/,10X,'NODE',5X,'BOTTOM FRICTION FRIC',5X,/)
         DO I=1,NP
            WRITE(16,2087) I,FRIC(I)
 2087       FORMAT(7X,I6,6X,E17.10)
         END DO
      ELSE
         WRITE(16,3504)
 3504    FORMAT(/,5X,'NODAL BOTTOM FRICTION VALUES ARE AVAILABLE',&
             /,6X,' IN UNIT 21 INPUT FILE')
      ENDIF
!
      RETURN
!     ----------------------------------------------------------------
      END SUBROUTINE ReadLegacyBottomFrictionFile
!     ----------------------------------------------------------------


!-----------------------------------------------------------------------
      END MODULE NodalAttributes
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

