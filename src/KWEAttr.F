C***********************************************************************
C     
C     SUBROUTINE ReadKWEAttr()
C     
C     This is a temporaly subroutine reading attributes for the DG rainfall-runoff
C     
C     The KWE attributes include nodal/elemental attributes, 
C     and the format is slightly changed from fort.13
C     
C     Written by Younghun Kang (04-13-2020)
C     
C***********************************************************************

      MODULE KWEAttr
      
      USE SIZES
      USE GLOBAL
      USE DG
      
      CHARACTER(len=80) AttrName ! string where the attribute name is stored
      CHARACTER(len=80) header   ! string where alphanumeric file id is stored
      CHARACTER(len=80) Skipped ! data in unit 13 we do not need
      
      CHARACTER(len=80) HydraulicCondUnits
      CHARACTER(len=80) CapillaryPUnits
      CHARACTER(len=80) MoistDecificUnits
      CHARACTER(len=80) IntDepthUnits
      CHARACTER(len=80) ManningKWE2Units
      CHARACTER(len=80) ManningKWE1Units
      
      INTEGER HydraulicCondNoOfVals
      INTEGER CapillaryPNoOfVals
      INTEGER MoistDecificNoOfVals
      INTEGER IntDepthNoOfVals
      INTEGER ManningKWE2NoOfVals
      INTEGER ManningKWE1NoOfVals
      
      INTEGER HydraulicCondNoOfAttrNodes
      INTEGER CapillaryPNoOfAttrNodes
      INTEGER MoistDecificNoOfAttrNodes
      INTEGER IntDepthNoOfAttrNodes
      INTEGER ManningKWE2NoOfAttrNodes
      INTEGER ManningKWE1NoOfAttrNodes
      
      REAL(SZ) HydraulicCondDefVal
      REAL(SZ) CapillaryPDefVal
      REAL(SZ) MoistDecificDefVal
      REAL(SZ) IntDepthDefVal
      REAL(SZ) ManningKWE2DefVal
      REAL(SZ) ManningKWE1DefVal
      
      ! REAL(SZ), ALLOCATABLE :: KS(:), HF(:), MD(:), IntDp(:)
      ! REAL(SZ), ALLOCATABLE :: manningKWE1(:), manningKWE2(:)
      
      CONTAINS
      
      SUBROUTINE ReadKWEAttr()
      
      READ(131,'(A80)') header
      READ(131,*) NAttr
      Do k=1, NAttr
         READ(131,'(A80)') AttrName
         WRITE(16,'(9X,A80)') AttrName
         WRITE(16,260) 
 260     FORMAT(14X,'was found!',/) 
         SELECT CASE (TRIM(ADJUSTL(AttrName)))
         ! CASE("precipitation_rate")
            ! FoundPrec = .True.
            ! READ(131,'(A80)') PrecUnits
            ! READ(131,*) PrecNoOfVals, PrecNoOfAttrNodes
            ! READ(131,*) PrecDefVal
         CASE("hydraulic_conductivity")
            FoundHydraulicCond = .True.
            READ(131,'(A80)') HydraulicCondUnits
            READ(131,*) HydraulicCondNoOfVals, HydraulicCondNoOfAttrNodes
            ALLOCATE (KS(HydraulicCondNoOfAttrNodes))
            READ(131,*) HydraulicCondDefVal
         CASE("capillary_pressure")
            FoundCapillaryP = .True.
            READ(131,'(A80)') CapillaryPUnits
            READ(131,*) CapillaryPNoOfVals, CapillaryPNoOfAttrNodes
            ALLOCATE (HF(CapillaryPNoOfAttrNodes))
            READ(131,*) CapillaryPDefVal
         CASE("moisture_decific")
            FoundMoistDecific = .True.
            READ(131,'(A80)') MoistDecificUnits
            READ(131,*) MoistDecificNoOfVals, MoistDecificNoOfAttrNodes
            ALLOCATE (MD(MoistDecificNoOfAttrNodes))
            READ(131,*) MoistDecificDefVal
         CASE("intercept_depth")
            FoundIntDepth = .True.
            READ(131,'(A80)') IntDepthUnits
            READ(131,*) IntDepthNoOfVals, IntDepthNoOfAttrNodes
            ALLOCATE (IntDP(IntDepthNoOfAttrNodes))
            READ(131,*) IntDepthDefVal
         CASE("mannings_n_at_element")
            FoundManningKWE2 = .True.
            READ(131,'(A80)') ManningKWE2Units
            READ(131,*) ManningKWE2NoOfVals, ManningKWE2NoOfAttrNodes
            ALLOCATE (ManningKWE2(ManningKWE2NoOfAttrNodes))
            READ(131,*) ManningKWE2DefVal
         CASE("mannings_n_at_open_channel")
            FoundManningKWE1 = .True.
            READ(131,'(A80)') ManningKWE1Units
            READ(131,*) ManningKWE1NoOfVals, ManningKWE1NoOfAttrNodes
            ALLOCATE (ManningKWE1(ManningKWE1NoOfAttrNodes))
            READ(131,*) ManningKWE1DefVal
         CASE DEFAULT     
            ! WRITE(16,1001)          ! Nodal Attributes file
            ! WRITE(16,1021) AttrName ! contains invalid name
            ! IF (NScreen.ne.0.and.MyProc.eq.0) THEN
               ! WRITE(ScreenUnit,1001)       
               ! WRITE(ScreenUnit,1021) AttrName 
            ! ENDIF
            READ(131,'(A80)') Skipped  ! skip the Unitsfor the invalid name
            READ(131,'(A80)') Skipped  ! skip the NoOfVals for invalid name
         END SELECT
      END DO

      DO k=1, NAttr
         WRITE(16,280) k
 280     FORMAT(/,9X,'Attribute ',I2,':')
         READ(131,'(A80)') AttrName
         READ(131,*) NumNodesNotDefault
         WRITE(16,'(14X,A80)') AttrName
         SELECT CASE (TRIM(ADJUSTL(AttrName)))
C         CASE("precipitation_rate")
C            CALL LoadAttrVecKWE(Prec, PrecDefVal, NumNodesNotDefault, 
C     &                   PrecNoOfAttrNodes,NScreen, MyProc, NAbOut)
         CASE("hydraulic_conductivity")
            PRINT *,  HydraulicCondNoOfAttrNodes
            CALL LoadAttrVecKWE(KS, HydraulicCondDefVal, NumNodesNotDefault, 
     &                   HydraulicCondNoOfAttrNodes,NScreen, MyProc, NAbOut)
         CASE("capillary_pressure")
            CALL LoadAttrVecKWE(HF, CapillaryPDefVal, NumNodesNotDefault, 
     &                   CapillaryPNoOfAttrNodes,NScreen, MyProc, NAbOut)
         CASE("moisture_decific")
            CALL LoadAttrVecKWE(MD, MoistDecificDefVal, NumNodesNotDefault, 
     &                   MoistDecificNoOfAttrNodes,NScreen, MyProc, NAbOut)
         CASE("intercept_depth")
            CALL LoadAttrVecKWE(IntDp, IntDepthDefVal, NumNodesNotDefault, 
     &                   IntDepthNoOfAttrNodes,NScreen, MyProc, NAbOut)
         CASE("mannings_n_at_element")
            CALL LoadAttrVecKWE(manningKWE2, ManningKWE2DefVal, NumNodesNotDefault, 
     &                   ManningKWE2NoOfAttrNodes,NScreen, MyProc, NAbOut)
         CASE("mannings_n_at_open_channel")
            CALL LoadAttrVecKWE(manningKWE1, ManningKWE1DefVal, NumNodesNotDefault,
     &                   ManningKWE1NoOfAttrNodes,NScreen, MyProc, NAbOut)
         CASE DEFAULT
            ! SkipDataSet = .True.
            ! WRITE(16,1001)      ! Nodal Attributes file
            ! WRITE(16,1021) AttrName ! contains invalid name
            ! IF (NScreen.ne.0.and.MyProc.eq.0) THEN
               ! WRITE(ScreenUnit,1001)       
               ! WRITE(ScreenUnit,1021) AttrName 
            ! ENDIF
         END SELECT
         ! IF (SkipDataSet) THEN
            ! DO L=1, NumNodesNotDefault
               ! READ(131,*) Skipped
            ! END DO
            ! WRITE(16,'(9X,A8)') 'Skipped.'
            ! SkipDataSet = .False.
         ! ELSE
            ! WRITE(16,'(/,9X,A18,A80)') 'Finished loading ', AttrName
         ! ENDIF
      END DO
      
      END SUBROUTINE ReadKWEAttr
      
C     ----------------------------------------------------------------
C         S U B R O U T I N E     L O A D  A T T R  V E C 
C     ----------------------------------------------------------------
C
C     jgf46.00 Subroutine to set a single set of nodal attributes to
C     their user-specified default values, then read the nondefault
C     values from the Nodal Attributes File (unit 13). This subroutine
C     is used for nodal attributes with only one value per node, hence
C     the suffix "vec" in the name.
C
C     ----------------------------------------------------------------
      SUBROUTINE LoadAttrVecKWE(AttributeData, Default, NumNodesNotDef, 
     &     NumAttrNodes,
     &     NScreen, MyProc, NAbOut)
      IMPLICIT NONE
      
      REAL(SZ), intent(in):: Default ! default value for all nodes
      INTEGER, intent(in) :: NumNodesNotDef ! number of nodes specified in file
      INTEGER, intent(in) :: NScreen ! 1 for debug info to screen (unit 6)
      INTEGER, intent(in) :: MyProc  ! in parallel, only MyProc=0 i/o to screen
      INTEGER, intent(in) :: NAbOut  ! 1 to abbrev. output to unit 16
C
      INTEGER NodeNum            ! node number listed in the file
      INTEGER I
      INTEGER, INTENT(IN) :: NumAttrNodes
      
      REAL(SZ), intent(out), dimension(NumAttrNodes) :: AttributeData
C
C     Set all values to user-specified default values.
      IF (NABOUT.EQ.0) WRITE(16,1001) Default 
      DO i=1, NumAttrNodes
         AttributeData(i) = Default 
      END DO
C
      IF (NABOUT.EQ.0) WRITE(16,1005) 
      DO i=1, NumNodesNotDef
         READ(131,*) NodeNum, AttributeData(NodeNum)
         IF (NABOUT.EQ.0) WRITE(16,1010) NodeNum, AttributeData(NodeNum)     
      END DO
C
 1001 FORMAT(/,10X,'Set all nodes to the default value of ',E16.8,/)
 1005 FORMAT(/,10X,'Now setting the following nodes to these values:',
     &     /,10X,'NODE',5X,'DATA',5X/)
 1010 FORMAT(7X,I6,6X,E16.8)
C
      RETURN
C     ----------------------------------------------------------------
      END SUBROUTINE LoadAttrVecKWE
C     ----------------------------------------------------------------
      
      END MODULE KWEAttr