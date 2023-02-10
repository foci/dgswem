
C***********************************************************************
C     
C     SUBROUTINE DG_RUNOFF_TIMESTEP(IT)
C     
C     This is the main subroutine for the DG-SAKE.
C     
C     Written by Younghun Kang (01-28-2020)
C     
C***********************************************************************

      SUBROUTINE DG_RUNOFF_TIMESTEP(IT)

C.....Use appropriate modules
      
      USE GLOBAL
      USE DG

#ifdef CMPI
      USE MESSENGER_ELEM 
#endif

#ifdef SWAN
Casey 101118: Add variables for coupling to SWAN.
      USE Couple2Swan, ONLY: ComputeWaveDrivenForces,
     &                       InterpoWeight
#endif

      IMPLICIT NONE
      
C.....Declare local variables


      INTEGER IT,L,GED,NBOREL,NNBORS,NDRYNBORS,NQEDS,Istop,k,j,kk,i,mm
      INTEGER Detected,flagger,store
      REAL(SZ) QBCT1,QBCT2,QBC1,QBC2,ZP(3),QXP(3),QXP_OLD(3),
     $     QYP(3),QYP_OLD(3),ZP_OLD(3),ell_1,ell_2,ell_3
      REAL(SZ) DPAVG, ITDT, ARK, BRK, CRK, DRK, ERK,s,alph,eta,seta
      Real(sz) tempx,tempy,rev,Ox,Oy,time_at,C_0,C_1,sig,sigma

      Real(SZ),allocatable :: XB(:),YB(:),radial(:)

      Allocate ( XB(MNE),YB(MNE),radial(MNE) )

C.....Compute the current time

      ITDT = IT*DTDP
      TIME_A = ITDT + STATIM*86400.D0
      TIMEH  = ITDT + (STATIM - REFTIM)*86400.D0

C.....Update ETA1 (shintaro: do we need this?)

      DO I = 1,MNP
         ETA1(I) = ETA2(I)
      ENDDO

C.....Begin Runge-Kutta time stepper

      DO IRK = 1,NRK
         
C.......Compute the DG time and DG ramp

#ifdef RKSSP
 
         TIMEDG = TIME_A - DTDP + DTVD(IRK)*DTDP
         RAMPDG = 1.D0
         RAMPExtFlux = 1.0d0
         IF (NRAMP.GE.1) THEN
            IF (NRAMP.EQ.1) THEN
               RAMPDG = TANH((2.D0*((IT-1) +
     & DTVD(IRK))*DTDP/86400.D0)/DRAMP)
               RAMPExtFlux = TANH((2.D0*((IT-1) +
     & DTVD(IRK))*DTDP/86400.D0)/DRAMPExtFlux)
            ENDIF
            IF (NRAMP.EQ.2) THEN
               RAMPDG = TANH((2.D0*((IT-1) +
     & DTVD(IRK))*DTDP/86400.D0)/DRAMP)
               RAMPExtFlux = TANH((2.D0*((IT-1) +
     & DTVD(IRK))*DTDP/86400.D0)/DRAMPExtFlux)
            ENDIF
            IF (NRAMP.EQ.3) THEN
               WRITE(*,*) 'NRAMP = 3 not supported '
               STOP
            ENDIF
         ENDIF
         RAMP = RAMPDG
         
C.......Compute no-normal flow edges

         IF (NLEDS.GT.0)  CALL LAND_EDGE_RUNOFF(IT)

C.......Compute radiation edges

         IF (NREDS.GT.0)  CALL RADIATION_EDGE_RUNOFF(IT)
         
C.......Compute internal edges

         CALL INTERNAL_EDGE_RUNOFF(IT)

C.......Compute elements to finish building the rhs
         
         CALL RHS_DG_RUNOFF(IT)

C.......SSP Runge-Kutta Time Scheme

         DO I = 1,IRK
            ARK = ATVD(IRK,I)
            BRK = BTVD(IRK,I)*DT
            DO J = 1,NE
               DO K = 1,DOFS(J)

                  HWD(K,J,IRK+1) = HWD(K,J,irk+1) + ARK*HWD(K,J,I)
     &                 + BRK*RHS_HWD(K,J,I) 
     
               ENDDO
            ENDDO
            
            DO J = 1,NCEDS
               GED = NCEDN(J)
               DO K = 1,ph+1
                  AWD(K,GED,IRK+1) = AWD(K,GED,irk+1) + ARK*AWD(K,GED,I)
     &                 + BRK*RHS_AWD(K,GED,I)
               ENDDO
            ENDDO
         ENDDO
#endif

#ifdef CMPI
         CALL UPDATER_ELEM_MOD(HWD,QX,QY,IRK+1,1)
         CALL UPDATER_CHANNEL(AWD,HWD,IRK+1,1)
#endif
     
#ifdef SLOPEALL
         CALL SLOPELIMITER()
#endif
      ENDDO
      

C.....Update variables for next time step     

      DO J = 1,MNE
         
         DO K = 1,DOFS(J)
            ! ZE(K,J,1) = ZE(K,J,NRK+1)
            ! QX(K,J,1) = QX(K,J,NRK+1)
            ! QY(K,J,1) = QY(K,J,NRK+1)
            
            HWD(K,J,1) = HWD(K,J,NRK+1)

         ENDDO
      ENDDO

      DO KK = 2,NRK+1
         DO J = 1,MNE
            DO K = 1,DOFH
               ! ZE(K,J,KK) = 0.D0
               ! QX(K,J,KK) = 0.D0
               ! QY(K,J,KK) = 0.D0
               ! RHS_ZE(K,J,KK-1) = 0.D0
               ! RHS_QX(K,J,KK-1) = 0.D0
               ! RHS_QY(K,J,KK-1) = 0.D0
               
               HWD(K,J,KK) = 0.D0
               RHS_HWD(K,J,KK-1) = 0.D0

            ENDDO
         ENDDO
      ENDDO
      
      
C.....Update variables for next time step     

      DO J = 1,NCEDS
         GED = NCEDN(J)
         DO K = 1,ph+1
            AWD(K,GED,1) = AWD(K,GED,NRK+1)
         ENDDO
      ENDDO

      DO KK = 2,NRK+1
         DO J = 1,NCEDS
            GED = NCEDN(J)
            DO K = 1,ph+1
               AWD(K,GED,KK) = 0.D0
               RHS_AWD(K,GED,KK-1) = 0.D0
            ENDDO
         ENDDO
      ENDDO

      RETURN
      END SUBROUTINE
