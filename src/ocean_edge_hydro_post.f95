C***********************************************************************
C
C     SUBROUTINE OCEAN_EDGE_HYDRO_POST( )
C
C     It seems elements on open ocean boundaries must be constant 
C     elements to avoid intabilities.  This subroutine simply eliminates
C     high order constituents.
C
C     Written by Shintaro Bunya based on ocean_edge_hydro.F (02-26-2007) 
C
C***********************************************************************

      SUBROUTINE OCEAN_EDGE_HYDRO_POST( )

C.....Use appropriate modules

      USE GLOBAL
      USE DG

      IMPLICIT NONE

C.....Declare local variables

      INTEGER L, LED, GED

      DO 1000 L=1, needs
      
C.....Retrieve the global and local edge number

      GED = NEEDN(L)
      LED = NEDSD(1,GED)

C.....Retrieve the elements which share the edge

      EL_IN = NEDEL(1,GED)
        
      ZE(2:DOF,EL_IN,IRK+1) = 0.D0
      QX(2:DOF,EL_IN,IRK+1) = 0.D0
      QY(2:DOF,EL_IN,IRK+1) = 0.D0

1000  CONTINUE
      RETURN
      END SUBROUTINE
