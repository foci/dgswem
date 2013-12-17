      MODULE DIFF45_41

      USE SIZES, ONLY : SZ

      PUBLIC

C---------------------------------------------------------
C  Imported from sizes.F
C

C     Dimension of vertical FE mesh (To interface 2D & 3D)
      INTEGER :: MNODES

C
C
C---------------------------------------------------------

C---------------------------------------------------------
C  Imported from global.F
C
      INTEGER,ALLOCATABLE ::    NEITAB(:,:),NEITABELE(:,:)
      
C.....for buoyancy forcing
      REAL(SZ),ALLOCATABLE ::   VIDBCPDX1(:), VIDBCPDY1(:), DASigT(:)
      INTEGER,ALLOCATABLE ::    LBArray_Pointer(:)

C.....nominal density of water RHOWAT0
      REAL(SZ), PARAMETER ::  RHOWAT0=1000.D0

C.....Sigma T value of reference density
      REAL(SZ), PARAMETER ::  SigT0=RHOWAT0-1000.D0 

      LOGICAL CBaroclinic

      REAL(SZ),ALLOCATABLE ::   MOM_LV_X(:),MOM_LV_Y(:),GWCE_LV(:)
C
C
C---------------------------------------------------------


C-------------------end of data declarations----------------------------------C



      END MODULE DIFF45_41

