!***********************************************************************
!
!     SUBROUTINE CALC_NORMAL( )
!
!     This subroutine calculates the unit normal vector for a given edge
!
!     Written by someone at UT
!
!     Comments and slight modifications made by Ethan Kubatko (03-07-05)
!
!***********************************************************************

      SUBROUTINE CALC_NORMAL()

!.....Use appropriate modules

      USE GLOBAL
      USE DG
      
      IMPLICIT NONE
      
!.....Declare local variables

      INTEGER IEL, IED,i

!.....Loop over the edges

      DO I = 1,NEDGES
      
!.....Retrieve the node numbers for the given edge

        N1 = NEDNO(1,I)
        N2 = NEDNO(2,I)

!.....Compute an average SFAC to adjust normal for CPP coordinates

        SAV = (SFAC(N1) + SFAC(N2))/2.0
        
!.....Compute the length of the given egde
        
        XLEN(I) = SQRT((Y(N2) - Y(N1))**2.D0                   + (X(N2) - X(N1))**2.D0)
     
!.....Compute the components of the normal vector

       !COSNX(I) = SAV*(Y(N2) - Y(N1))/XLEN(I)
       COSNX(I) = (Y(N2) - Y(N1))/XLEN(I)
       SINNX(I) = -(X(N2) - X(N1))/XLEN(I)
      ENDDO
        
      RETURN
      END


