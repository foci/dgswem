!***********************************************************************
!
!     SUBROUTINE TIDAL_POTENTIAL()
!
!     This subroutine computes the tidal potential forcing terms at each
!     node.  Note that the tidal potential term is defined linearly over
!     an element regardless of the p used in the DG calculation.
!
!     (Taken from timestep.dg%f)
!
!     Comments and general clean-up by Ethan Kubatko (06-02-2005)
!
!***********************************************************************

      SUBROUTINE TIDAL_POTENTIAL(dg_here,global_here)

!.....Use appropriate modules

      USE GLOBAL
      USE DG
      
      IMPLICIT NONE

      type (dg_type) :: dg_here
      type (global_type) :: global_here

      INTEGER j,i
      
      dg_here%TIMEH_DG = global_here%TIMEH - global_here%DTDP + dg_here%DTVD(dg_here%IRK)*global_here%DTDP
      
!.....Loop over the nodes

      DO I = 1,global_here%NP
      
!.....Initialize tidal potential terms

        global_here%TIP2(I) = 0.D0

!.....Loop over the tidal potential constituents

        DO J = 1,global_here%NTIF

          IF (global_here%PERT(J).EQ.0.) THEN
            global_here%NCYC = 0
          ELSE
            global_here%NCYC = INT(dg_here%TIMEH_DG/global_here%PERT(J))
          ENDIF
          global_here%ARGT    = global_here%AMIGT(J)*(dg_here%TIMEH_DG - global_here%NCYC*global_here%PERT(J)) + global_here%FACET(J)
          global_here%TPMUL   = dg_here%RAMPDG*global_here%ETRF(J)*global_here%TPK(J)*global_here%FFT(J)
          global_here%SALTMUL = dg_here%RAMPDG*global_here%FFT(J)
          global_here%NA      = NINT(0.00014/global_here%AMIGT(J))
        
!.....Semi-diurnal species

          IF (global_here%NA.EQ.1) THEN

            global_here%ARGTP    = global_here%ARGT + 2.D0*global_here%SLAM(I)
            global_here%ARGSALT  = global_here%ARGT - global_here%SALTPHA(J,I)
            global_here%CCSFEA   = COS(global_here%SFEA(I))
            global_here%CCSFEA   = global_here%CCSFEA*global_here%CCSFEA
            global_here%TIP2(I) = global_here%TIP2(I) + global_here%TPMUL*global_here%CCSFEA*COS(global_here%ARGTP)&
                             + global_here%SALTMUL*global_here%SALTAMP(J,I)*COS(global_here%ARGSALT)

          ENDIF
        
!.....Diurnal species
        
          IF (global_here%NA.EQ.2) THEN

            global_here%ARGTP    = global_here%ARGT + global_here%SLAM(I)
            global_here%ARGSALT  = global_here%ARGT - global_here%SALTPHA(J,I)
            global_here%S2SFEA   = SIN(2.D0*global_here%SFEA(I))
            global_here%TIP2(I) = global_here%TIP2(I) + global_here%TPMUL*global_here%S2SFEA*COS(global_here%ARGTP)&
                             + global_here%SALTMUL*global_here%SALTAMP(J,I)*COS(global_here%ARGSALT)
     
          ENDIF

        ENDDO

      ENDDO
      
      RETURN
      END SUBROUTINE

