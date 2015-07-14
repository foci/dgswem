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
!
!asey 101118: Allow SWAN to handle wave refraction as a nodal attribute.
!             I've placed these changes outside the #ifdef SWAN flags
!             because we want to be able to use the same fort.13 files
!             for both ADCIRC and SWAN+ADCIRC runs.  This way, the new
!             nodal attribute will be processed but only applied when
!             ADCIRC is coupled to SWAN.
      type nodalattr_type
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
!      INTEGER i        ! node loop counter
!      INTEGER j        ! attribute values loop counter
!      INTEGER k        ! attribute loop counter
!

   end type nodalattr_type
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      CONTAINS !- - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 







!     ----------------------------------------------------------------
!      F U N C T I O N   T A U 0  N O D A L  V A L U E 
!     ----------------------------------------------------------------
!
!     jgf46.27 Function to calculate nodalattr_here%tau0 based on the scheme selection
!     and the depth. This assumes that Scheme is negative. 
!       
!     ----------------------------------------------------------------
      REAL(SZ) FUNCTION Tau0NodalValue(Scheme, Depth)
      IMPLICIT NONE
      REAL(SZ) Scheme         
      REAL(SZ) Depth
!
      IF (Scheme.eq.-2.d0) THEN 
!     Smoothly varying nodalattr_here%tau0 with depth.
         IF(Depth.GE.200.) Tau0NodalValue=0.005
         IF((Depth.LT.200.).AND.(Depth.GE.1.)) THEN
            Tau0NodalValue=1./Depth
         ENDIF
         IF(Depth.LT.1.) Tau0NodalValue=1.0
      ELSE
!     Abrupt variation in nodalattr_here%tau0 with depth.           
         IF(Depth.LE.10.) Tau0NodalValue=0.020d0    
         IF(Depth.GT.10.) Tau0NodalValue=0.005d0
      ENDIF         
!     ----------------------------------------------------------------
      END FUNCTION Tau0NodalValue
!     ----------------------------------------------------------------

#if 0

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
      REAL(SZ), intent(inout), dimension(NP) :: TK! depth avg. nodalattr_here%fric.
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
      IF (nodalattr_here%LoadManningsN) THEN
         DO I=1, NP
            nodalattr_here%FRIC(I)=g*nodalattr_here%ManningsN(I)**2.d0&
                /( ( DP(I)+IFNLFA*ETA2(I) )**(1.d0/3.d0) ) ! sb46.28sb02
            !sb46.28sb02  Lower limit is applied here.
            IF(nodalattr_here%FRIC(I).LT.LL) THEN
               nodalattr_here%FRIC(I) = LL
            ENDIF
         ENDDO
      ENDIF
!
!     ... Convert nodalattr_here%Chezy to Cd, if necessary.
      IF (nodalattr_here%LoadChezy) THEN
         DO I=1,NP
            nodalattr_here%FRIC(I)=G*(nodalattr_here%Chezy(I)**2)
         END DO
      ENDIF
!
!     Step 1. Apply friction arising from turbulent viscous interaction
!     with the sea floor.
      DO I=1, NP
         UV1=SQRT(UU1(I)*UU1(I)+VV1(I)*VV1(I))
         H1=DP(I)+IFNLFA*ETA2(I)
         TK(I)= nodalattr_here%FRIC(I)*&
             ( nodalattr_here%IFLINBF +       &! linear&
             (UV1/H1) * (nodalattr_here%IFNLBF &! nonlinear&
             + nodalattr_here%IFHYBF*(1+(nodalattr_here%HBREAK/H1)**nodalattr_here%FTHETA)**(nodalattr_here%FGAMMA/nodalattr_here%FTHETA))) ! hybrid
      END DO
!
!     Step 2. Apply friction arising from flow interaction with bridge
!     pilings, if required.
      IF (nodalattr_here%LoadBridgePilings) THEN
         DO I=1, NP
            UV1=SQRT(UU1(I)*UU1(I)+VV1(I)*VV1(I))
            H1=DP(I)+IFNLFA*ETA2(I)
            Fr=UV1*UV1/(G*H1)
            BK = nodalattr_here%BridgePilings(I,1) 
            BALPHA = nodalattr_here%BridgePilings(I,2) 
            BDELX = nodalattr_here%BridgePilings(I,3) 
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

#endif


#if 0

!     ----------------------------------------------------------------
!                      S U B R O U T I N E
!         R E A D  L E G A C Y  S T A R T  D R Y  F I L E
!     ----------------------------------------------------------------
!
!     jgf46.00 Subroutine to load up the legacy nodalattr_here%startdry file (unit
!     12). This is just a cut-and-paste from the section of the
!     READ_INPUT subroutine that did the same thing. This subroutine is
!     never called. It is vestigial and listed here purely as reference
!     material.
!
!     ----------------------------------------------------------------
      SUBROUTINE ReadLegacyStartDryFile(s,NP, NScreen, ScreenUnit, &
          MyProc, NAbOut)
      IMPLICIT NONE
      type (sizes_type) :: s
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

      OPEN(s%fort12unit,FILE=TRIM(s%INPUTDIR)//'/'//'fort.12')
!
!...  READ nodalattr_here%STARTDRY INFORMATION FROM UNIT 12 
      READ(s%fort12unit,'(A80)') AGRID2
      WRITE(s%fort16unit,2038) AGRID2
2038  FORMAT(5X,'nodalattr_here%STARTDRY FILE IDENTIFICATION : ',A80,/)
      READ(s%fort12unit,*) NE2,NP2
!
!...  CHECK THAT NE2 AND NP2 MATCH WITH GRID FILE 
!      IF((NE2.NE.NE).OR.(NP2.NE.NP)) THEN
       IF(NP2.NE.NP) THEN
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,9900)
         WRITE(s%fort16unit,9900)
 9900    FORMAT(////,1X,'!!!!!!!!!!  FATAL ERROR  !!!!!!!!!',&
             //,1X,'THE PARAMETER NE2 AND NP2 MUST MATCH NE AND NP ',&
             /,1X,'USER MUST CHECK FORT.12 INPUT FILE ',&
             //,1X,'!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!',//)
         STOP
      ENDIF
!
!...  READ IN nodalattr_here%STARTDRY CODE VALUES
      DO I=1,NP
         READ(s%fort12unit,*) JKI,DUM1,DUM2,nodalattr_here%STARTDRY(JKI)
         IF(JKI.NE.I) THEN
            IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,99805)
            WRITE(s%fort16unit,99805)
99805       FORMAT(////,1X,'!!!!!!!!!!  WARNING - NONFATAL ',&
                'INPUT ERROR  !!!!!!!!!',&
                //,1X,'YOUR NODE NUMBERING IS NOT SEQUENTIAL ',&
                'CHECK YOUR UNIT 12 INPUT FILE CAREFULLY',//)
         ENDIF
      END DO
!     
!...  CLOSE UNIT 12 FILE       
      CLOSE(s%fort12unit)     
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
      SUBROUTINE ReadLegacyBottomFrictionFile(s, NP, NScreen, ScreenUnit,MyProc, NAbOut)
        use sizes
      IMPLICIT NONE
      type (sizes_type) :: s
      INTEGER, intent(in) :: NP ! number of nodes in grid file
      INTEGER, intent(in) :: NScreen ! nonzero for debug info to screen
      INTEGER, intent(in) :: ScreenUnit ! i/o for debug info to screen
      INTEGER, intent(in) :: MyProc  ! in parallel, only MyProc=0 i/o to screen
      INTEGER, intent(in) :: NAbOut  ! 1 to abbrev. output to unit 16
      
      CHARACTER(len=80) AFRIC  ! user's comment/description line
      INTEGER NHG    ! node number from file

      OPEN(s%fort21unit,FILE=TRIM(s%INPUTDIR)//'/'//'fort.21')
      READ(s%fort21unit,'(A80)') AFRIC
      DO I=1,NP
         READ(s%fort21unit,*) NHG,nodalattr_here%FRIC(NHG)
         IF(NHG.NE.I) THEN
            IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,99803)
            WRITE(s%fort16unit,99803)
99803       FORMAT(////,1X,'!!!!!!!!!!  WARNING - FATAL ',&
                'INPUT ERROR  !!!!!!!!!',//,1X,&
                'YOUR NODAL FRICTION NUMBERING IS NOT SEQUENTIAL ',&
                /,1X,'CHECK YOUR UNIT 21 INPUT FILE CAREFULLY',//,1X,&
                '!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!',//)
            STOP
         ENDIF
      END DO
      WRITE(s%fort16unit,3601) AFRIC
 3601 FORMAT(/,5X,'FRICTION FILE IDENTIFICATN : ',A80,/)
      IF(NABOUT.NE.1) THEN
         WRITE(s%fort16unit,2080)
 2080    FORMAT(/,10X,'NODE',5X,'BOTTOM FRICTION nodalattr_here%FRIC',5X,/)
         DO I=1,NP
            WRITE(s%fort16unit,2087) I,nodalattr_here%FRIC(I)
 2087       FORMAT(7X,I6,6X,E17.10)
         END DO
      ELSE
         WRITE(s%fort16unit,3504)
 3504    FORMAT(/,5X,'NODAL BOTTOM FRICTION VALUES ARE AVAILABLE',&
             /,6X,' IN UNIT 21 INPUT FILE')
      ENDIF
!
      RETURN
!     ----------------------------------------------------------------
      END SUBROUTINE ReadLegacyBottomFrictionFile
!     ----------------------------------------------------------------

#endif

!-----------------------------------------------------------------------
      END MODULE NodalAttributes
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------


