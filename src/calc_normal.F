C***********************************************************************
C
C     SUBROUTINE CALC_NORMAL( )
C
C     This subroutine calculates the unit normal vector for a given edge
C
C     Written by someone at UT
C
C     Comments and slight modifications made by Ethan Kubatko (03-07-05)
C
C***********************************************************************

      SUBROUTINE CALC_NORMAL()

C.....Use appropriate modules

      USE GLOBAL, only : x, y, sfac, MNE, neled
      USE DG, only : xlen, cosnx, sinnx, nedges, nedno,  sav, elem_size
      use sizes, only : sz
      use iso_fortran_env, only : output_unit
      
      IMPLICIT NONE
      
C.....Declare local variables

      INTEGER IEL, IED,i, n1, n2
      real(sz) l1, l2, l3, s

C.....Loop over the edges

      DO I = 1,NEDGES
      
C.....Retrieve the node numbers for the given edge

        N1 = NEDNO(1,I)
        N2 = NEDNO(2,I)

C.....Compute an average SFAC to adjust normal for CPP coordinates

        SAV = (SFAC(N1) + SFAC(N2))/2.0
        
C.....Compute the length of the given egde
        
        XLEN(I) = SQRT((Y(N2) - Y(N1))**2.0
     &                    + (X(N2) - X(N1))**2.0)
     
C.....Compute the components of the normal vector

       !COSNX(I) = SAV*(Y(N2) - Y(N1))/XLEN(I)
       COSNX(I) = (Y(N2) - Y(N1))/XLEN(I)
       SINNX(I) = -(X(N2) - X(N1))/XLEN(I)
      ENDDO

      DO i = 1, MNE
         l1 = XLEN(NELED(1,i))
         l2 = XLEN(NELED(2,i))
         l3 = XLEN(NELED(3,i))
         s = 0.5 * (l1 + l2 + l3)
         elem_size(i) = 2 * sqrt( (s-l1)*(s-l2)*(s-l3) / s ) ! diameter of the incircle
      ENDDO
      write(output_unit, '( "Smallest triangle diameter = ", F9.2, " m")') minval(elem_size)
        
      RETURN
      END


