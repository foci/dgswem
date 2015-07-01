!     ----------------------------------------------------------------
!     S U B R O U T I N E     I N I T  N  A  M O D U L E
!     ----------------------------------------------------------------
!
!     jgf46.00 Subroutine to initialize the variables in the nodal
!     attributes module.
!
!     ----------------------------------------------------------------
      SUBROUTINE InitNAModule(nodalattr_here)
        USE NodalAttributes
      IMPLICIT NONE
      type (nodalattr_type) :: nodalattr_here
      integer :: i,j,k
!
!asey 101118: Make changes compact.
      nodalattr_here%LoadSwanWaveRefrac     = .FALSE.
      nodalattr_here%FoundSwanWaveRefrac    = .FALSE.
      nodalattr_here%SwanWaveRefracNoOfVals = 1
      nodalattr_here%SwanWaveRefracDefVal   = 0.D0
!
      nodalattr_here%LoadTau0           = .False.
      nodalattr_here%LoadStartDry       = .False. 
      nodalattr_here%LoadDirEffRLen     = .False.
      nodalattr_here%LoadManningsN      = .False. 
      nodalattr_here%LoadQuadraticFric  = .False.
      nodalattr_here%LoadChezy          = .False.
      nodalattr_here%LoadBridgePilings  = .False.
      nodalattr_here%LoadCanopyCoef     = .False.
      nodalattr_here%LoadGeoidOffset    = .False.
      nodalattr_here%LoadEVM            = .False.
      nodalattr_here%LoadEVC            = .False.
!
      nodalattr_here%FoundTau0           = .False.
      nodalattr_here%FoundStartDry       = .False.
      nodalattr_here%FoundDirEffRLen     = .False.
      nodalattr_here%FoundManningsN      = .False. 
      nodalattr_here%FoundQuadraticFric  = .False.
      nodalattr_here%FoundChezy          = .False.
      nodalattr_here%FoundBridgePilings  = .False.
      nodalattr_here%FoundCanopyCoef     = .False.
      nodalattr_here%FoundGeoidOffset    = .False.
      nodalattr_here%FoundEVM            = .False.
      nodalattr_here%FoundEVC            = .False.
!      
      nodalattr_here%Tau0NoOfVals          = 1
      nodalattr_here%StartDryNoOfVals      = 1
      nodalattr_here%DirEffRLenNoOfVals    = 12
      nodalattr_here%QuadraticFricNoOfVals = 1
      nodalattr_here%ChezyNoOfVals         = 1
      nodalattr_here%ManningsNNoOfVals     = 1
      nodalattr_here%BridgePilingsNoOfVals = 4
      nodalattr_here%CanopyCoefNoOfVals    = 1
      nodalattr_here%GeoidOffsetNoOfVals   = 1
      nodalattr_here%EVMNoOfVals           = 1
      nodalattr_here%EVCNoOfVals           = 1
!
      nodalattr_here%Tau0DefVal             = 0.0
      nodalattr_here%StartDryDefVal         = 0.0
      DO j=1, nodalattr_here%DirEffRLenNoOfVals
         nodalattr_here%DirEffRLenDefVal(j) = 0.0
      END DO
      nodalattr_here%CanopyCoefDefVal       = 0.0
      nodalattr_here%QuadraticFricDefVal    = 0.0
      DO j=1, nodalattr_here%BridgePilingsNoOfVals
         nodalattr_here%BridgePilingsDefVal(j) = 0.0
      END DO
      nodalattr_here%ChezyDefVal            = 0.0
      nodalattr_here%ManningsNDefVal        = 0.0
      nodalattr_here%GeoidOffsetDefVal      = 0.0
      nodalattr_here%EVMDefVal              = 0.0
      nodalattr_here%EVCDefVal              = 0.0
!
      nodalattr_here%HBREAK=1.d0
      nodalattr_here%FTHETA=1.d0
      nodalattr_here%FGAMMA=1.d0
!
      nodalattr_here%ESLM=0.0
      nodalattr_here%ESLC=0.0
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
      SUBROUTINE ReadNodalAttr(s,nodalattr_here, NScreen, ScreenUnit, MyProc, NAbOut)
        use sizes
        use NodalAttributes      
        IMPLICIT NONE
      type (sizes_type) :: s
      type (nodalattr_type) :: nodalattr_here
      integer :: i,j,k

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
!     Check to make sure that nodalattr_here%NWP is a valid number.
      IF (nodalattr_here%NWP.LT.0) THEN
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) THEN
            WRITE(ScreenUnit,9972)
            WRITE(ScreenUnit,*) 'nodalattr_here%NWP =',nodalattr_here%NWP
            WRITE(ScreenUnit,9728)
            WRITE(ScreenUnit,9973)
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'nodalattr_here%NWP =',nodalattr_here%NWP
         WRITE(16,9728)
         WRITE(16,9973)
 9728    FORMAT(/,1X,'Your selection of nodalattr_here%NWP (a UNIT 15 input ',&
             'parameter) is not an allowable value')
         STOP            ! We're toast.
      ENDIF
!
!     Check to see if there are nodal attributes to be read in. If not,
!     simply return.
      IF (nodalattr_here%NWP.EQ.0) THEN
         WRITE(16,231) nodalattr_here%NWP
 231     FORMAT(/,5X,'nodalattr_here%NWP = ',I2,&
            /,9X,'A Nodal Attributes File (unit 13)',&
             /,9X,'will not be used.')
         RETURN
      ENDIF
!
!     Otherwise, get on with it.
      WRITE(16,232) nodalattr_here%NWP
 232  FORMAT(/,5X,'nodalattr_here%NWP = ',I2,&
          /,9X,'Must read Nodal Attributes File (unit 13).')
!
!     Determine if the Nodal Attributes File exists.
      INQUIRE(FILE=TRIM(s%INPUTDIR)//'/'//'fort.13',EXIST=NAFound)
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
      WRITE(16,235) nodalattr_here%NWP
 235  FORMAT(/,9X,'Need to load ',I2,' nodal attribute(s):')
      DO k=1,nodalattr_here%NWP
         READ(15,'(A80)') AttrName
         WRITE(16,'(14X,A80)') AttrName
         SELECT CASE (TRIM(ADJUSTL(AttrName)))
         CASE("primitive_weighting_in_continuity_equation")
            nodalattr_here%LoadTau0 = .True.
         CASE("surface_submergence_state")
            nodalattr_here%LoadStartDry = .True.
         CASE("quadratic_friction_coefficient_at_sea_floor")
            nodalattr_here%LoadQuadraticFric = .True.
         CASE("surface_directional_effective_roughness_length") 
            nodalattr_here%LoadDirEffRLen = .True. 
         CASE("surface_canopy_coefficient") 
            nodalattr_here%LoadCanopyCoef = .True.
         CASE("bridge_pilings_friction_parameters") 
            nodalattr_here%LoadBridgePilings = .True.
         CASE("mannings_n_at_sea_floor")
            nodalattr_here%LoadManningsN = .True. 
         CASE("chezy_friction_coefficient_at_sea_floor") 
            nodalattr_here%LoadChezy = .True.
         CASE("sea_surface_height_above_geoid") 
            nodalattr_here%LoadGeoidOffset = .True.
         CASE&
        ("average_horizontal_eddy_viscosity_in_sea_water_wrt_depth")
            nodalattr_here%LoadEVM = .True.
         CASE&
        ("average_horizontal_eddy_diffusivity_in_sea_water_wrt_depth")
            nodalattr_here%LoadEVC = .True.
!asey 101118: Allow SWAN to handle wave refraction as a nodal attribute.
         CASE("wave_refraction_in_swan")
            nodalattr_here%LoadSwanWaveRefrac = .TRUE.
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
      OPEN(UNIT=13, FILE=s%DIRNAME//'/'//'fort.13', &
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
      READ(13,*) nodalattr_here%NumOfNodes     ! number of nodes according to unit 13
      READ(13,*) nodalattr_here%NAttr          ! number of attributes in the unit 13 file
      DO k=1, nodalattr_here%NAttr
         READ(13,'(A80)') AttrName
         WRITE(16,'(9X,A80)') AttrName
         WRITE(16,260) 
 260     FORMAT(14X,'was found!',/) 
         SELECT CASE (TRIM(ADJUSTL(AttrName)))
         CASE("primitive_weighting_in_continuity_equation")
            nodalattr_here%FoundTau0 = .True.
            READ(13,'(A80)') nodalattr_here%Tau0Units   
            READ(13,*) nodalattr_here%Tau0NoOfVals
            READ(13,*) nodalattr_here%Tau0DefVal
         CASE("surface_submergence_state")
            nodalattr_here%FoundStartDry = .True. 
            READ(13,'(A80)') nodalattr_here%StartDryUnits      
            READ(13,*) nodalattr_here%StartDryNoOfVals
            READ(13,*) nodalattr_here%StartDryDefVal            
         CASE("quadratic_friction_coefficient_at_sea_floor")
            nodalattr_here%FoundQuadraticFric = .True.
            READ(13,'(A80)') nodalattr_here%QuadraticFricUnits 
            READ(13,*) nodalattr_here%QuadraticFricNoOfVals
            READ(13,*) nodalattr_here%QuadraticFricDefVal
         CASE("surface_directional_effective_roughness_length") 
            nodalattr_here%FoundDirEffRLen = .True. 
            READ(13,'(A80)') nodalattr_here%DirEffRLenUnits  
            READ(13,*) nodalattr_here%DirEffRLenNoOfVals
            READ(13,*) (nodalattr_here%DirEffRLenDefVal(j),j=1,12)    
         CASE("surface_canopy_coefficient") 
            nodalattr_here%FoundCanopyCoef = .True. 
            READ(13,'(A80)') nodalattr_here%CanopyCoefUnits  
            READ(13,*) nodalattr_here%CanopyCoefNoOfVals
            READ(13,*) nodalattr_here%CanopyCoefDefVal 
         CASE("bridge_pilings_friction_parameters") 
            nodalattr_here%FoundBridgePilings = .True.
            READ(13,'(A80)') nodalattr_here%BridgePilingsUnits  
            READ(13,*) nodalattr_here%BridgePilingsNoOfVals
            READ(13,*) (nodalattr_here%BridgePilingsDefVal(j),j=1,12)  
         CASE("mannings_n_at_sea_floor")
            nodalattr_here%FoundManningsN = .True. 
            READ(13,'(A80)') nodalattr_here%ManningsNUnits 
            READ(13,*) nodalattr_here%ManningsNNoOfVals
            READ(13,*) nodalattr_here%ManningsNDefVal 
         CASE("chezy_friction_coefficient_at_sea_floor")
            nodalattr_here%FoundChezy = .True.
            READ(13,'(A80)') nodalattr_here%ChezyUnits     
            READ(13,*) nodalattr_here%ChezyNoOfVals
            READ(13,*) nodalattr_here%ChezyDefVal          
         CASE("sea_surface_height_above_geoid")
            nodalattr_here%FoundGeoidOffset = .True.
            READ(13,'(A80)') nodalattr_here%GeoidOffsetUnits     
            READ(13,*) nodalattr_here%GeoidOffsetNoOfVals
            READ(13,*) nodalattr_here%GeoidOffsetDefVal
         CASE&
         ("average_horizontal_eddy_viscosity_in_sea_water_wrt_depth")
            nodalattr_here%FoundEVM = .True.
            READ(13,'(A80)') nodalattr_here%EVMUnits     
            READ(13,*) nodalattr_here%EVMNoOfVals
            READ(13,*) nodalattr_here%EVMDefVal
         CASE&
        ("average_horizontal_eddy_diffusivity_in_sea_water_wrt_depth")
            READ(13,'(A80)') nodalattr_here%EVCUnits     
            READ(13,*) nodalattr_here%EVCNoOfVals
            READ(13,*) nodalattr_here%EVCDefVal
!asey 101118: Allow SWAN to handle wave refraction as a nodal attribute.
         CASE("wave_refraction_in_swan")
            nodalattr_here%FoundSwanWaveRefrac = .TRUE.
            READ(13,'(A80)') nodalattr_here%SwanWaveRefracUnits
            READ(13,*) nodalattr_here%SwanWaveRefracNoOfVals
            READ(13,*) nodalattr_here%SwanWaveRefracDefVal
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
      IF(((nodalattr_here%LoadTau0).and.(.not.nodalattr_here%FoundTau0)).or.&
        ((nodalattr_here%LoadStartDry).and.(.not.nodalattr_here%FoundStartDry)).or.&
        ((nodalattr_here%LoadQuadraticFric).and.&
         (.not.nodalattr_here%FoundQuadraticFric)).or.&
        ((nodalattr_here%LoadDirEffRLen).and.&
         (.not.nodalattr_here%FoundDirEffRLen)).or.&
        ((nodalattr_here%LoadCanopyCoef).and.&
         (.not.nodalattr_here%FoundCanopyCoef)).or.&
        ((nodalattr_here%LoadBridgePilings).and.&
         (.not.nodalattr_here%FoundBridgePilings)).or.&
        ((nodalattr_here%LoadManningsN).and.&
         (.not.nodalattr_here%FoundManningsN)).or.&
        ((nodalattr_here%LoadGeoidOffset).and.&
         (.not.nodalattr_here%FoundGeoidOffset)).or.&
        ((nodalattr_here%LoadChezy).and.(.not.nodalattr_here%FoundChezy)).or.&
        ((nodalattr_here%LoadEVM).and.(.not.nodalattr_here%FoundEVM)).or.&
!asey 101118: Allow SWAN to handle wave refraction as a nodal attribute.&
        ((nodalattr_here%LoadSwanWaveRefrac).and.(.not.nodalattr_here%FoundSwanWaveRefrac)).or.&
        ((nodalattr_here%LoadEVC).and.(.not.nodalattr_here%FoundEVC)) ) THEN
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
      ALLOCATE(nodalattr_here%TAU0VAR(nodalattr_here%NumOfNodes),nodalattr_here%TAU0BASE(nodalattr_here%NumOfNodes)) ! jjw&sb46.39sb01
      ALLOCATE(nodalattr_here%STARTDRY(nodalattr_here%NumOfNodes)) 
      ALLOCATE(nodalattr_here%FRIC(nodalattr_here%NumOfNodes))
      ALLOCATE(nodalattr_here%z0land(nodalattr_here%NumOfNodes,nodalattr_here%DirEffRLenNoOfVals))
      ALLOCATE(nodalattr_here%vcanopy(nodalattr_here%NumOfNodes))
      ALLOCATE(nodalattr_here%BridgePilings(nodalattr_here%NumOfNodes,nodalattr_here%BridgePilingsNoOfVals))
      ALLOCATE(nodalattr_here%GeoidOffset(nodalattr_here%NumOfNodes))
      ALLOCATE(nodalattr_here%Chezy(nodalattr_here%NumOfNodes))
      ALLOCATE(nodalattr_here%ManningsN(nodalattr_here%NumOfNodes))
      ALLOCATE(nodalattr_here%EVM(nodalattr_here%NumOfNodes))
      ALLOCATE(nodalattr_here%EVC(nodalattr_here%NumOfNodes))
!asey 101118: Allow SWAN to handle wave refraction as a nodal attribute.
      ALLOCATE(nodalattr_here%SwanWaveRefrac(nodalattr_here%NumOfNodes))
!
!     Now read each of the attributes required by the model parameter
!     (unit 15) file and skip past the others.
      WRITE(16,270) nodalattr_here%NWP
 270  FORMAT(/,9X,'Now reading ',I2,' nodal attribute(s).')
      DO k=1, nodalattr_here%NAttr
         WRITE(16,280) k
 280     FORMAT(/,9X,'Attribute ',I2,':')
         READ(13,'(A80)') AttrName
         READ(13,*) NumNodesNotDefault
         WRITE(16,'(14X,A80)') AttrName
         SELECT CASE (TRIM(ADJUSTL(AttrName)))
         CASE("primitive_weighting_in_continuity_equation")
            IF (nodalattr_here%LoadTau0) THEN
               CALL LoadAttrVec(nodalattr_here%TAU0VAR, nodalattr_here%Tau0DefVal,&
                   NumNodesNotDefault, NScreen, MyProc, NAbOut,nodalattr_here%NumOfNodes)
            ELSE
               SkipDataSet = .True.
            ENDIF
         CASE("surface_submergence_state")
            IF (nodalattr_here%LoadStartDry) THEN
               CALL LoadAttrVec(nodalattr_here%STARTDRY, nodalattr_here%StartDryDefVal, &
                    NumNodesNotDefault, NScreen, MyProc, NAbOut,nodalattr_here%NumOfNodes)
            ELSE
               SkipDataSet = .True.
            ENDIF
         CASE("quadratic_friction_coefficient_at_sea_floor")
            IF (nodalattr_here%LoadQuadraticFric) THEN
               CALL LoadAttrVec(nodalattr_here%FRIC, nodalattr_here%QuadraticFricDefVal, &
                   NumNodesNotDefault, NScreen, MyProc, NAbOut,nodalattr_here%NumOfNodes)
            ELSE
               SkipDataSet = .True.
            ENDIF
         CASE("surface_directional_effective_roughness_length") 
            IF (nodalattr_here%LoadDirEffRLen) THEN
               CALL LoadAttrMat(nodalattr_here%z0land, nodalattr_here%DirEffRLenNoOfVals,  &
                   nodalattr_here%DirEffRLenDefVal, NumNodesNotDefault,  &
                   NScreen, MyProc, NAbOut,nodalattr_here%NumOfNodes)
            ELSE
               SkipDataSet = .True.
            ENDIF
         CASE("surface_canopy_coefficient") 
            IF (nodalattr_here%LoadCanopyCoef) THEN
               CALL LoadAttrVec(nodalattr_here%vcanopy, nodalattr_here%CanopyCoefDefVal,&
                    NumNodesNotDefault, NScreen, MyProc, NAbOut,nodalattr_here%NumOfNodes)
            ELSE
               SkipDataSet = .True.
            ENDIF
         CASE("bridge_pilings_friction_parameters") 
            IF (nodalattr_here%LoadBridgePilings) THEN
               CALL LoadAttrMat(nodalattr_here%BridgePilings, nodalattr_here%BridgePilingsNoOfVals, &
                   nodalattr_here%BridgePilingsDefVal,  NumNodesNotDefault, &
                   NScreen, MyProc, NAbOut,nodalattr_here%NumOfNodes)
            ELSE
               SkipDataSet = .True.
            ENDIF
         CASE("mannings_n_at_sea_floor")
            IF (nodalattr_here%LoadManningsN) THEN
               CALL LoadAttrVec(nodalattr_here%ManningsN, nodalattr_here%ManningsNDefVal, &
                   NumNodesNotDefault, NScreen, MyProc, NAbOut,nodalattr_here%NumOfNodes)
            ELSE
               SkipDataSet = .True.
            ENDIF
         CASE("chezy_friction_coefficient_at_sea_floor")
            IF (nodalattr_here%LoadChezy) THEN
               CALL LoadAttrVec(nodalattr_here%Chezy, nodalattr_here%ChezyDefVal, &
                    NumNodesNotDefault, NScreen, MyProc, NAbOut,nodalattr_here%NumOfNodes)
            ELSE
               SkipDataSet = .True.
            ENDIF
         CASE("sea_surface_height_above_geoid")
            IF (nodalattr_here%LoadGeoidOffset) THEN
               CALL LoadAttrVec(nodalattr_here%GeoidOffset, nodalattr_here%GeoidOffsetDefVal, &
                    NumNodesNotDefault, NScreen, MyProc, NAbOut,nodalattr_here%NumOfNodes)
            ELSE
               SkipDataSet = .True.
            ENDIF
         CASE&
        ("average_horizontal_eddy_viscosity_in_sea_water_wrt_depth")
            IF (nodalattr_here%LoadEVM) THEN
               CALL LoadAttrVec(nodalattr_here%EVM, nodalattr_here%EVMDefVal, &
                    NumNodesNotDefault, NScreen, MyProc, NAbOut,nodalattr_here%NumOfNodes)
            ELSE
               SkipDataSet = .True.
            ENDIF
         CASE&
        ("average_horizontal_eddy_diffusivity_in_sea_water_wrt_depth")
            IF (nodalattr_here%LoadEVC) THEN
               CALL LoadAttrVec(nodalattr_here%EVC, nodalattr_here%EVCDefVal, &
                    NumNodesNotDefault, NScreen, MyProc, NAbOut,nodalattr_here%NumOfNodes)
            ELSE
               SkipDataSet = .True.
            ENDIF
!asey 101118: Allow SWAN to handle wave refraction as a nodal attribute.
         CASE("wave_refraction_in_swan")
            IF (nodalattr_here%LoadSwanWaveRefrac) THEN
               CALL LoadAttrVec(nodalattr_here%SwanWaveRefrac, nodalattr_here%SwanWaveRefracDefVal,&
                   NumNodesNotDefault, NScreen, MyProc, NAbOut,nodalattr_here%NumOfNodes)
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
 1002 FORMAT('ERROR: The legacy nodalattr_here%StartDry File (unit 12)')
 1003 FORMAT('ERROR: Spatially Varying nodalattr_here%Fric. Coeff. File (unit 21)')
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
      SUBROUTINE ApplyDirectionalWindReduction(nodalattr_here,NodeNumber,WindDragCo,WindMag,BathymetricDepth,Elevation,CutOffDepth,G,WindX,WindY)
      USE SIZES
      USE NodalAttributes
      IMPLICIT NONE
      type (nodalattr_type) :: nodalattr_here
      integer :: i,j,k
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
      z0l=nodalattr_here%z0land(NodeNumber,idir)
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
      WindX = nodalattr_here%vcanopy(NodeNumber)*fr*WindX 
      WindY = nodalattr_here%vcanopy(NodeNumber)*fr*WindY 
!
      RETURN
!     ----------------------------------------------------------------
      END SUBROUTINE ApplyDirectionalWindReduction
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
      SUBROUTINE LoadAttrVec(AttributeData, Default, NumNodesNotDef,NScreen, MyProc, NAbOut,NumOfNodes)
        USE SIZES, ONLY: SZ
      IMPLICIT NONE
      REAL(SZ), intent(out), dimension(NumOfNodes) :: AttributeData
      REAL(SZ), intent(in):: Default ! default value for all nodes
      INTEGER, intent(in) :: NumNodesNotDef ! number of nodes specified in file
      INTEGER, intent(in) :: NScreen ! 1 for debug info to screen (unit 6)
      INTEGER, intent(in) :: MyProc  ! in parallel, only MyProc=0 i/o to screen
      INTEGER, intent(in) :: NAbOut  ! 1 to abbrev. output to unit 16
!
      integer :: i,j,k
      INTEGER NodeNum            ! node number listed in the file
      integer :: NumOfNodes
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
      SUBROUTINE LoadAttrMat(AttributeData, NumCol, Default,NumNodesNotDef, NScreen, MyProc, NAbOut, NumOfNodes)
        USE SIZES, ONLY: SZ
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
      integer :: NumOfNodes

      integer :: i,j,k
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
      SUBROUTINE InitNodalAttr(global_here,nodalattr_here,DP, NP, G, NScreen, ScreenUnit,MyProc, NAbOut)
      USE GLOBAL
      USE SIZES, ONLY: SZ
      USE NodalAttributes
      IMPLICIT NONE

      type (global_type) :: global_here
      type (nodalattr_type) :: nodalattr_here

      INTEGER, intent(in) :: NP ! number of nodes in the grid file 
      REAL(SZ), intent(in), dimension(NP) :: DP ! array of bathymetric depths
      REAL(SZ), intent(in):: G  ! gravitational acceleration
      INTEGER, intent(in) :: NScreen ! nonzero for debug info to screen
      INTEGER, intent(in) :: ScreenUnit ! i/o for debug info to screen
      INTEGER, intent(in) :: MyProc  ! in parallel, only MyProc=0 i/o to screen
      INTEGER, intent(in) :: NAbOut  ! 1 to abbrev. output to unit 16
      integer :: i,j,k
!
!
!     ERROR CHECK: If a nodal attributes file is being used, check to
!     see that the number of nodes in the nodal attribute file is the
!     same as the number of nodes in the grid file.
      IF (nodalattr_here%NWP.NE.0.AND.nodalattr_here%NumOfNodes.NE.NP) THEN
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,9900)
         WRITE(16,9900)
 9900    FORMAT(////,1X,'!!!!!!!!!!  FATAL ERROR  !!!!!!!!!',&
             //,1X,'The number of nodes in the grid file (unit 14) and'&
             /,1X,'the nodal attributes file (unit 13) must match.',&
             //,1X,'!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!',//)
         WRITE(16,*) 'np, nodalattr_here%numofnodes ',np,nodalattr_here%NumOfNodes
         STOP                   ! We're toast.
      ENDIF
!
!     ERROR CHECK: If nodalattr_here%Chezy, Manning's or Quadratic friction was loaded
!     from the nodal attributes file, nodalattr_here%NOLIBF must be  >= 1. 
      IF ((nodalattr_here%LoadChezy.or.nodalattr_here%LoadManningsN.or.nodalattr_here%LoadQuadraticFric).and.&
          nodalattr_here%NoLiBF.eq.0) THEN
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,9800)
         WRITE(16,9800) nodalattr_here%nolibf
 9800    FORMAT(////,1X,'!!!!!!!!!!  FATAL ERROR  !!!!!!!!!',&
             //,1X,'Nonlinear bottom friction coefficients were loaded'&
             /,1X,'from the nodal attributes file (unit 13), so ',&
             /,1X,'nodalattr_here%NoLiBF must be set to 1. It is set to ',i2,' in',&
             /,1X,'the model parameter (unit 15) file.',&
             //,1X,'!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!',//)
         STOP                   ! We're toast.
      ENDIF
!     
!     I N I T    S T A R T D R Y
      IF (.not.nodalattr_here%LoadStartDry) THEN
!         IF (nodalattr_here%NWP.eq.0) THEN
!            ALLOCATE(nodalattr_here%STARTDRY(NP))
!         ENDIF
         DO I=1, NP
            nodalattr_here%STARTDRY(I) = 0.0D0
         ENDDO
      ENDIF
!
!     I N I T     T A U 0
      IF (.not.nodalattr_here%LoadTau0) THEN
!         IF (nodalattr_here%NWP.eq.0) THEN
!            ALLOCATE(nodalattr_here%TAU0VAR(NP),nodalattr_here%TAU0BASE(NP)) ! jjw&sb46.39sb01
!         ENDIF
!
!     jgf46.25 If input nodalattr_here%tau0 is positive, set all nodes to that value.
         IF (nodalattr_here%Tau0.ge.0) THEN
            DO I=1,NP
               nodalattr_here%Tau0Var(I)=nodalattr_here%Tau0
            END DO
            WRITE(16,7) nodalattr_here%Tau0
 7          FORMAT(/,5X,&
                'A SPATIALLY CONSTANT WEIGHTING COEFFICIENT (nodalattr_here%Tau0)'&
                ,/,5X,' WILL BE USED IN THE GENERALIZED WAVE',&
                ' CONTINUITY EQUATION.',&
                /,5X,'nodalattr_here%Tau0 = ',E15.8,2X,'1/sec',/)
         ELSE
            DO I=1,NP
               nodalattr_here%Tau0Var(I)=Tau0NodalValue(nodalattr_here%Tau0,DP(I))
            ENDDO
            
            IF(nodalattr_here%Tau0.eq.-2) THEN
               WRITE(16,6) ! spatially vary nodalattr_here%tau0 according to hard coded scheme
               WRITE(16,62) ! description of scheme
 62            FORMAT(/,5X,'IF DEPTH > 200           nodalattr_here%Tau0 = 0.005',    &
                   /,5X,'IF 200   > DEPTH > 1     nodalattr_here%Tau0 = 1/DEPTH  ',    &
                   /,5X,'IF 1     > DEPTH         nodalattr_here%Tau0 = 1.0 ')
            ELSE
               WRITE(16,6) ! spatially vary nodalattr_here%tau0 according to hard coded scheme
               WRITE(16,61) ! description of scheme
 61            FORMAT(/,5X,' IF DEPTH GE 10           -> nodalattr_here%TAU0 = 0.005',    &
                   /,5X,' IF DEPTH LT 10           -> nodalattr_here%TAU0 = 0.020',/)
            ENDIF
         ENDIF
      ENDIF
 6    FORMAT(/,5X,'A SPATIALLY VARIABLE WEIGHTING COEFFICIENT (nodalattr_here%Tau0)'&
          ,/,5X,' WILL BE USED IN THE GENERALIZED WAVE',&
          ' CONTINUITY EQUATION.',&
          /,5x,'THIS VALUE WILL BE DETERMINED AS FOLLOWS:')

!
!     jgf46.27 If we have already loaded the nodalattr_here%tau0 values directly from
!     the nodal attributes file, check to see if the default value was
!     negative. If so, this indicates that nodal values of nodalattr_here%tau0 that
!     were not explicitly set in the nodal attributes file should be set
!     accordiing to one of the hard-coded nodalattr_here%tau0 schemes.

      IF (nodalattr_here%LoadTau0.and.nodalattr_here%Tau0DefVal.lt.0) THEN
         DO I=1,NP
            IF (nodalattr_here%Tau0Var(I).lt.0) THEN
               nodalattr_here%Tau0Var(I)=Tau0NodalValue(nodalattr_here%Tau0DefVal,DP(I))
            ENDIF
         ENDDO
      ENDIF
!
!     jjw&sb46.39.sb01 Save the nodalattr_here%Tau0Var in nodalattr_here%Tau0Base for later use
      DO I=1,NP
        nodalattr_here%Tau0Base(I) = nodalattr_here%Tau0Var(I)
      ENDDO
!
!     jjw&sb46.38.sb01 If nodalattr_here%tau0 is loaded from nodal attributes file and 
!     nodalattr_here%Tau0 is -3, time-varing nodalattr_here%tau0 optimizer will be applied in timestep.F
      if(nodalattr_here%LoadTau0.AND.nodalattr_here%Tau0.eq.-3.d0) then
         WRITE(16,7)
      endif
 8    FORMAT(/,5X,'A SPATIALLY TEMPORALLY VARIABLE OPTIMIZED '&
          ,/,5X,' WEIGHTING COEFFICIENT (nodalattr_here%Tau0) WILL BE USED '&
          ,/,5X,' IN THE GENERALIZED WAVE CONTINUITY EQUATION.',/)
!
!     I N I T   B O T T O M   F R I C T I O N
      IF(nodalattr_here%NOLIBF.EQ.0) THEN
         nodalattr_here%IFNLBF=0
         nodalattr_here%IFLINBF=1
         nodalattr_here%IFHYBF=0
      ENDIF
      IF(nodalattr_here%NOLIBF.EQ.1) THEN
         nodalattr_here%IFNLBF=1
         nodalattr_here%IFLINBF=0
         nodalattr_here%IFHYBF=0
      ENDIF
      IF(nodalattr_here%NOLIBF.EQ.2) THEN
         nodalattr_here%IFNLBF=0
         nodalattr_here%IFLINBF=0
         nodalattr_here%IFHYBF=1
      ENDIF
!     
!     Initialize bottom friction if it was not loaded from unit 13.
      IF((.not.nodalattr_here%LoadQuadraticFric).and.(.not.nodalattr_here%LoadManningsN).and.&
          (.not.nodalattr_here%LoadChezy)) THEN
         IF (nodalattr_here%NoLiBF.eq.0) nodalattr_here%CF=nodalattr_here%Tau
!     If a nodal attributes file was read, nodalattr_here%FRIC was allocated there.
!         IF (nodalattr_here%NWP.eq.0) THEN
!            ALLOCATE(nodalattr_here%FRIC(NP))
!         ENDIF
         DO I=1,NP
            nodalattr_here%FRIC(I)=nodalattr_here%CF
         END DO
      ENDIF

!
!     Initialize bridge pilings.
!
      IF (nodalattr_here%LoadBridgePilings) THEN
         DO I=1, NP
            IF (nodalattr_here%BridgePilings(I,1).ne.0) THEN ! only for nodes w/piers
               nodalattr_here%BridgePilings(I,3) = 4.d0 * &
                   nodalattr_here%BridgePilings(I,3) / nodalattr_here%BridgePilings(I,4)
            ENDIF
         END DO
      ENDIF
!
!     I N I T   E D D Y   V I S C O S I T Y  &  D I F F U S I V I T Y
      IF (.not.nodalattr_here%LoadEVM) THEN
!	 IF (nodalattr_here%NWP .eq. 0) THEN
!          	ALLOCATE(nodalattr_here%EVM(NP))
!	 ENDIF
         DO I=1,NP
            nodalattr_here%EVM(I)=nodalattr_here%ESLM
         END DO
      ENDIF
      DO I=1,NP
         global_here%EVMSUM=global_here%EVMSUM+ABS(nodalattr_here%EVM(I))
      ENDDO
      IF (.not.nodalattr_here%LoadEVC.and.nodalattr_here%ESLC.ne.0) THEN
         DO I=1,NP
            nodalattr_here%EVC(I)=nodalattr_here%ESLC
         END DO
      ENDIF
!
      RETURN
!     ----------------------------------------------------------------
      END SUBROUTINE InitNodalAttr
!     ----------------------------------------------------------------
