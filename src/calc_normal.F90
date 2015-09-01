!***********************************************************************
!
!     SUBROUTINE CALC_NORMAL( )
!
!     This subroutine calculates the unit normal vector for a given global_here%edge
!
!     Written by someone at UT
!
!     Comments and slight modifications made by Ethan Kubatko (03-07-05)
!
!***********************************************************************

      SUBROUTINE CALC_NORMAL(dg_here,global_here)

!.....Use appropriate modules

      USE GLOBAL
      USE DG
      
      IMPLICIT NONE

      type (dg_type) :: dg_here
      type (global_type) :: global_here

!.....Declare local variables

      INTEGER IEL, IED,i

!.....Loop over the edges

      DO I = 1,dg_here%NEDGES
      
!.....Retrieve the node numbers for the given global_here%edge

        global_here%N1 = dg_here%NEDNO(1,I)
        global_here%N2 = dg_here%NEDNO(2,I)

!.....Compute an average global_here%SFAC to adjust normal for CPP coordinates

        dg_here%SAV = (global_here%SFAC(global_here%N1) + global_here%SFAC(global_here%N2))/2.0
        
!.....Compute the length of the given egde
        
        dg_here%XLEN(I) = SQRT((global_here%Y(global_here%N2) - global_here%Y(global_here%N1))**2.D0                   + (global_here%X(global_here%N2) - global_here%X(global_here%N1))**2.D0)
     
!.....Compute the components of the normal vector

       !dg_here%COSNX(I) = dg_here%SAV*(global_here%Y(global_here%N2) - global_here%Y(global_here%N1))/dg_here%XLEN(I)
       dg_here%COSNX(I) = (global_here%Y(global_here%N2) - global_here%Y(global_here%N1))/dg_here%XLEN(I)
       dg_here%SINNX(I) = -(global_here%X(global_here%N2) - global_here%X(global_here%N1))/dg_here%XLEN(I)
      ENDDO
        
      RETURN
      END


