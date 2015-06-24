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

      SUBROUTINE TIDAL_POTENTIAL()

!.....Use appropriate modules

      USE GLOBAL
      USE DG
      
      IMPLICIT NONE

      INTEGER j,i
      
      dg%TIMEH_DG = TIMEH - DTDP + dg%DTVD(dg%IRK)*DTDP
      
!.....Loop over the nodes

      DO I = 1,NP
      
!.....Initialize tidal potential terms

        TIP2(I) = 0.D0

!.....Loop over the tidal potential constituents

        DO J = 1,NTIF

          IF (PERT(J).EQ.0.) THEN
            NCYC = 0
          ELSE
            NCYC = INT(dg%TIMEH_DG/PERT(J))
          ENDIF
          ARGT    = AMIGT(J)*(dg%TIMEH_DG - NCYC*PERT(J)) + FACET(J)
          TPMUL   = dg%RAMPDG*ETRF(J)*TPK(J)*FFT(J)
          SALTMUL = dg%RAMPDG*FFT(J)
          NA      = NINT(0.00014/AMIGT(J))
        
!.....Semi-diurnal species

          IF (NA.EQ.1) THEN

            ARGTP    = ARGT + 2.D0*SLAM(I)
            ARGSALT  = ARGT - SALTPHA(J,I)
            CCSFEA   = COS(SFEA(I))
            CCSFEA   = CCSFEA*CCSFEA
            TIP2(I) = TIP2(I) + TPMUL*CCSFEA*COS(ARGTP)&
                             + SALTMUL*SALTAMP(J,I)*COS(ARGSALT)

          ENDIF
        
!.....Diurnal species
        
          IF (NA.EQ.2) THEN

            ARGTP    = ARGT + SLAM(I)
            ARGSALT  = ARGT - SALTPHA(J,I)
            S2SFEA   = SIN(2.D0*SFEA(I))
            TIP2(I) = TIP2(I) + TPMUL*S2SFEA*COS(ARGTP)&
                             + SALTMUL*SALTAMP(J,I)*COS(ARGSALT)
     
          ENDIF

        ENDDO

      ENDDO
      
      RETURN
      END SUBROUTINE

