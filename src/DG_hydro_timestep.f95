
!***********************************************************************
!     
!     SUBROUTINE DG_HYDRO_TIMESTEP(IT)
!     
!     This is the main subroutine for the DG hydro
!     
!     Written by Ethan Kubatko (06-11-2004)
!     
!-----------------------------------------------------------------------
!     
!     Modification history on hp_DG_ADCIRC since v9_sb1
!     
!     v??     - 08/05    - sb - parallelized
!     v??     - 08/05    - sb - wetting/drying
!     v10_sb5 - 10/05    - sb - wetting/drying algorithm is improved
!     v10_sb5 - 10/05    - sb - Ethan's slope limiter is consolidated
!     v21     - 11/11    - cem - major update/p_enrich/slope/rkc etc
!     v22     -  6/12    - cem - major update/sediment/diffusion/fluxes
!     Moving to versioning control
!     
!***********************************************************************

      SUBROUTINE DG_HYDRO_TIMESTEP(s,dg_here,IT)

!.....Use appropriate modules
      USE SIZES
      USE GLOBAL
      USE DG

#ifdef CMPI
      USE MESSENGER_ELEM 
#endif

#ifdef SWAN
!asey 101118: Add variables for coupling to SWAN.
      USE Couple2Swan, ONLY: ComputeWaveDrivenForces,&
                            InterpoWeight
#endif

      IMPLICIT NONE
      
!.....Declare local variables

      type (sizes_type) :: s
      type (dg_type) :: dg_here

      integer :: irk ! defining this locally
      INTEGER IT,L,GED,NBOREL,NNBORS,NDRYNBORS,Istop,k,j,kk,i,mm
      INTEGER Detected,flagger,store
      REAL(SZ) QBCT1,QBCT2,QBC1,QBC2,ZP(3),QXP(3),QXP_OLD(3),QYP(3),QYP_OLD(3),ZP_OLD(3),ell_1,ell_2,ell_3
      REAL(SZ) DPAVG, ITDT, ARK, BRK, CRK, DRK, ERK,s_dg,alph,eta,seta
      Real(sz) tempx,tempy,rev,Ox,Oy,time_at,C_0,C_1,sig,sigma

      Real(SZ),allocatable :: XB(:),YB(:),radial(:)

      Allocate ( XB(S%MNE),YB(S%MNE),radial(S%MNE) )

!.....Compute the current time

      ITDT = IT*DTDP
      TIME_A = ITDT + STATIM*86400.D0
      TIMEH  = ITDT + (STATIM - REFTIM)*86400.D0

!.....Update ETA1 (shintaro: do we need this?)

      DO I = 1,S%MNP
         ETA1(I) = ETA2(I)
      ENDDO

!.....Begin Runge-Kutta time stepper

      DO 100 irk = 1,dg_here%NRK
         
!.......Compute the DG time and DG ramp

#ifdef RKSSP

 
         dg_here%TIMEDG = TIME_A - DTDP + dg_here%DTVD(irk)*DTDP
         dg_here%RAMPDG = 1.D0
         RAMPExtFlux = 1.0d0
         IF (NRAMP.GE.1) THEN
            IF (NRAMP.EQ.1) THEN
               dg_here%RAMPDG = TANH((2.D0*((IT-1) + dg_here%DTVD(irk))*DTDP/86400.D0)/DRAMP)
               RAMPExtFlux = &
                   TANH((2.D0*((IT-1) + dg_here%DTVD(irk))*DTDP/86400.D0)/DRAMPExtFlux)
            ENDIF
            IF (NRAMP.EQ.2) THEN
               dg_here%RAMPDG = TANH((2.D0*((IT-1) + dg_here%DTVD(irk))*DTDP/86400.D0)/DRAMP)
               RAMPExtFlux = &
                   TANH((2.D0*((IT-1) + dg_here%DTVD(irk))*DTDP/86400.D0)/DRAMPExtFlux)
            ENDIF
            IF (NRAMP.EQ.3) THEN
               WRITE(*,*) 'NRAMP = 3 not supported '
               STOP
            ENDIF
         ENDIF
         RAMP = dg_here%RAMPDG
         
!.......Obtain the meteorological forcing

         IF (NWS.NE.0) CALL MET_FORCING(s,dg_here,IT)
        
        IF(NRS.GE.1) THEN
          IF(TIME_A.GT.RSTIME2) THEN
            RSTIME1=RSTIME2
            RSTIME2=RSTIME2+RSTIMINC
           !DO I=1,NP
           !  RSNX1(I)=RSNX2(I)
           !  RSNY1(I)=RSNY2(I)
           !END DO
#ifdef SWAN
!asey 090302: Added for coupling to SWAN.
            IF((NRS.EQ.3).AND.(irk.EQ.1)) THEN
              InterpoWeight = 1.0
              CALL ComputeWaveDrivenForces
!asey 090707: We want to extrapolate forward in time.  Load the latest (current) forces
!             into RSNX1/RSNY1, and then load the future forces into RSNX2/RSNY2.
              DO I=1,NP
                RSX = RSNX1(I)
                RSY = RSNY1(I)
                RSNX1(I) = RSNX2(I)
                RSNY1(I) = RSNY2(I)
                RSNX2(I) = RSNX2(I) + (RSNX2(I)-RSX)
                RSNY2(I) = RSNY2(I) + (RSNY2(I)-RSY)
              ENDDO
            ENDIF
#endif
          ENDIF
          RStRatio=(TIME_A-RSTIME1)/RSTIMINC
          DO I=1,NP
            RSX = RampMete*(RSNX1(I) + RStRatio*(RSNX2(I)-RSNX1(I)))
            RSY = RampMete*(RSNY1(I) + RStRatio*(RSNY2(I)-RSNY1(I)))
            WSX2(I) = WSX2(I) + RSX
            WSY2(I) = WSY2(I) + RSY
#ifdef SWAN
!asey 090302: Added these lines for output to the rads.64 file.
            RSNXOUT(I) = RSX
            RSNYOUT(I) = RSY
#endif
         ENDDO
      ENDIF
 
!.......Compute tidal potential terms

         IF (NTIP.NE.0) CALL TIDAL_POTENTIAL(dg_here)

!.......Compute LDG auxiliary equations
         
         IF (EVMSUM.NE.0.D0.or.dg_here%artdif.eq.1) CALL LDG_HYDRO(s,dg_here,IT)
         
!.......Compute elevation specified edges

         IF (dg_here%NEEDS.GT.0)  CALL OCEAN_EDGE_HYDRO(s,dg_here,IT)

!.......Compute no-normal flow edges

         IF (dg_here%NLEDS.GT.0)  CALL LAND_EDGE_HYDRO(s,dg_here,IT)

!.......Compute non-zero flow edges

         IF (dg_here%NFEDS.GT.0)  CALL FLOW_EDGE_HYDRO(s,dg_here,IT)
         
!.......Compute radiation edges

         IF (dg_here%NREDS.GT.0)  CALL RADIATION_EDGE_HYDRO(s,dg_here,IT)

!.......Compute internal barrier edges

         IF (dg_here%NIBEDS.GT.0) CALL IBARRIER_EDGE_HYDRO(s,dg_here,IT)
         
!.......Compute external barrier edges

         IF (dg_here%NIBEDS.GT.0) CALL EBARRIER_EDGE_HYDRO(s,dg_here,IT)
         
!.......Compute internal edges

         CALL INTERNAL_EDGE_HYDRO(s,dg_here,IT)

!.......Compute elements to finish building the rhs
         
         CALL RHS_DG_HYDRO(s,dg_here)

!.......SSP Runge-Kutta Time Scheme

         DO I = 1,irk
            ARK = dg_here%ATVD(irk,I)
            BRK = dg_here%BTVD(irk,I)*DT
            DO J = 1,NE
               DO K = 1,dg_here%DOFS(J)

#ifdef SED_LAY

                  do l = 1,layers
                     
                     dg_here%bed(K,J,irk+1,l) = dg_here%bed(K,J,irk+1,l) + ARK*dg_here%bed(K,J,I,l)&
                         + BRK*dg_here%RHS_bed(K,J,I,l)

                     !adjust for state variable represenation
                     dg_here%ZE(K,J,irk+1) = dg_here%ZE(K,J,irk+1) - BRK*dg_here%RHS_bed(K,J,I,l)


                  enddo
#endif

                  dg_here%ZE(K,J,irk+1) = dg_here%ZE(K,J,irk+1) + ARK*dg_here%ZE(K,J,I)&
                      + BRK*dg_here%RHS_ZE(K,J,I) 
                  dg_here%QX(K,J,irk+1) = dg_here%QX(K,J,irk+1) + ARK*dg_here%QX(K,J,I)&
                      + BRK*dg_here%RHS_QX(K,J,I)
                  dg_here%QY(K,J,irk+1) = dg_here%QY(K,J,irk+1) + dg_here%ATVD(irk,I)*dg_here%QY(K,J,I)&
                      + BRK*dg_here%RHS_QY(K,J,I)


!.......Compute the transported tracer term if flagged

#ifdef TRACE
                  dg_here%iota(K,J,irk+1) = dg_here%iota(K,J,irk+1) + ARK*dg_here%iota(K,J,I)&
                      + BRK*dg_here%RHS_iota(K,J,I)
#endif

!......Compute chemistry transported terms if flagged

#ifdef CHEM
                  dg_here%iota(K,J,irk+1) = dg_here%iota(K,J,irk+1) + ARK*dg_here%iota(K,J,I)&
                      + BRK*dg_here%RHS_iota(K,J,I)
                  dg_here%iota2(K,J,irk+1) = dg_here%iota2(K,J,irk+1) + ARK*dg_here%iota2(K,J,I)&
                      + BRK*dg_here%RHS_iota2(K,J,I)
#endif

!.......Compute the dynamic pressure if flagged

#ifdef DYNP
                  dg_here%dynP(K,J,irk+1) = dg_here%dynP(K,J,irk+1) + ARK*dg_here%dynP(K,J,I)&
                      + BRK*dg_here%RHS_dynP(K,J,I)
#endif


               ENDDO
            ENDDO
         ENDDO

#endif

#ifdef RKC

         print*,"Using RKC"

         dg_here%TIMEDG = TIME_A - DTDP + dg_here%RKC_c(irk)*DTDP
         dg_here%RAMPDG = 1.D0
         RAMPExtFlux = 1.0d0
         IF (NRAMP.GE.1) THEN
            IF (NRAMP.EQ.1) THEN
               dg_here%RAMPDG = TANH((2.D0*(IT + dg_here%RKC_c(irk))*DTDP/86400.D0)/DRAMP)
               RAMPExtFlux = &
                    TANH((2.D0*(IT + dg_here%RKC_c(irk))*DTDP/86400.D0)/DRAMPExtFlux)
            ENDIF
            IF (NRAMP.EQ.2) THEN
               dg_here%RAMPDG = TANH((2.D0*(IT + dg_here%RKC_c(irk))*DTDP/86400.D0)/DRAMP)
               RAMPExtFlux = &
                    TANH((2.D0*(IT + dg_here%RKC_c(irk))*DTDP/86400.D0)/DRAMPExtFlux)
            ENDIF
            IF (NRAMP.EQ.3) THEN
               WRITE(*,*) 'NRAMP = 3 not supported '
               STOP
            ENDIF
         ENDIF
         RAMP = dg_here%RAMPDG
         
!.......Obtain the meteorological forcing

         IF (NWS.NE.0) CALL MET_FORCING(s,dg_here,IT)
         
!.......Compute tidal potential terms

         IF (NTIP.NE.0) CALL TIDAL_POTENTIAL(dg_here)

!.......Compute LDG auxiliary equations
         
         IF (EVMSUM.NE.0.D0.or.dg_here%artdif.eq.1) CALL LDG_HYDRO(s,dg_here,IT)
         
!.......Compute elevation specified edges

         IF (dg_here%NEEDS.GT.0)  CALL OCEAN_EDGE_HYDRO(s,dg_here,IT)

!.......Compute no-normal flow edges

         IF (dg_here%NLEDS.GT.0)  CALL LAND_EDGE_HYDRO(s,dg_here,IT)

!.......Compute non-zero flow edges

         IF (dg_here%NFEDS.GT.0)  CALL FLOW_EDGE_HYDRO(s,dg_here,IT)
         
!.......Compute radiation edges

         IF (dg_here%NREDS.GT.0)  CALL RADIATION_EDGE_HYDRO(s,dg_here,IT)

!.......Compute internal barrier edges

         IF (dg_here%NIBEDS.GT.0) CALL IBARRIER_EDGE_HYDRO(s,dg_here,IT)
         
!.......Compute external barrier edges

         IF (dg_here%NIBEDS.GT.0) CALL EBARRIER_EDGE_HYDRO(s,dg_here,IT)
         
!.......Compute internal edges

         CALL INTERNAL_EDGE_HYDRO(s,dg_here,IT)

!.......Compute elements to finish building the rhs
         
         CALL RHS_DG_HYDRO(s,dg_here)

!.......RKC Time Scheme

         if (irk.eq.1) then

            BRK = dg_here%RKC_tildemu(1)*dt

            DO J = 1,NE
               DO K = 1,dg_here%DOFS(J)

#ifdef SED_LAY
                  do l = 1,layers
                     
                     dg_here%bed(K,J,irk+1,l) = dg_here%bed(K,J,irk+1,l) + dg_here%bed(K,J,1,l) &
                         + BRK*dg_here%RHS_bed(K,J,1,l)

                   !adjust for state variable represenation
                     dg_here%ZE(K,J,irk+1) = dg_here%ZE(K,J,irk+1) - BRK*dg_here%RHS_bed(K,J,1,l)
                     
                  enddo
#endif
                  
                  dg_here%ZE(K,J,irk+1) = dg_here%ZE(K,J,irk+1) + dg_here%ZE(K,J,1)&
                      + BRK*dg_here%RHS_ZE(K,J,1)
                  dg_here%QX(K,J,irk+1) = dg_here%QX(K,J,irk+1) + dg_here%QX(K,J,1)&
                      + BRK*dg_here%RHS_QX(K,J,1)
                  dg_here%QY(K,J,irk+1) = dg_here%QY(K,J,irk+1)  + dg_here%QY(K,J,1)&
                      + BRK*dg_here%RHS_QY(K,J,1)
                  
!.......Compute the transported tracer term if flagged
                  
#ifdef TRACE
                  dg_here%iota(K,J,irk+1) = dg_here%iota(K,J,irk+1) + dg_here%iota(K,J,1)&
                      + BRK*dg_here%RHS_iota(K,J,1)
#endif
                  
!......Compute chemistry transported terms if flagged
                  
#ifdef CHEM
                  dg_here%iota(K,J,irk+1) = dg_here%iota(K,J,irk+1) + dg_here%iota(K,J,1)&
                       + BRK*dg_here%RHS_iota(K,J,1)
                  dg_here%iota2(K,J,irk+1) = dg_here%iota2(K,J,irk+1) + dg_here%iota2(K,J,1)&
                       + BRK*dg_here%RHS_iota2(K,J,1)
#endif

!.......Compute the dynamic pressure if flagged
                  
#ifdef DYNP
                  dg_here%dynP(K,J,irk+1) = dg_here%dynP(K,J,irk+1) + dg_here%dynP(K,J,1)&
                      + BRK*dg_here%RHS_dynP(K,J,1)
#endif
               ENDDO
            ENDDO
            
         else

            ARK = 1.D0 - dg_here%RKC_mu(irk)-dg_here%RKC_nu(irk)
            CRK = dg_here%RKC_mu(irk)
            DRK = dg_here%RKC_nu(irk)

            BRK = dg_here%RKC_tildemu(irk)*dt
            ERK = dg_here%RKC_gamma(irk)*dt

             DO J = 1,NE
               DO K = 1,dg_here%DOFS(J)

#ifdef SED_LAY
                  do l = 1,layers
                     
                     dg_here%bed(K,J,irk+1,l) = dg_here%bed(K,J,irk+1,l) + ARK*dg_here%bed(K,J,1,l) + CRK*dg_here%bed(K,J,irk,l)&
                         + DRK*dg_here%bed(K,J,irk-1,l) + BRK*dg_here%RHS_bed(K,J,irk,l) + ERK*dg_here%RHS_bed(K,J,1,l)


                   !adjust for state variable represenation
                     dg_here%ZE(K,J,irk+1) = dg_here%ZE(K,J,irk+1) - BRK*dg_here%RHS_bed(K,J,irk,l) - ERK*dg_here%RHS_bed(K,J,1,l)
                     
                  enddo
#endif

                  dg_here%ZE(K,J,irk+1) = dg_here%ZE(K,J,irk+1) + ARK*dg_here%ZE(K,J,1) + CRK*dg_here%ZE(K,J,irk)&
                       + DRK*dg_here%ZE(K,J,irk-1) + BRK*dg_here%RHS_ZE(K,J,irk) + ERK*dg_here%RHS_ZE(K,J,1)
                  dg_here%QX(K,J,irk+1) = dg_here%QX(K,J,irk+1) + ARK*dg_here%QX(K,J,1) + CRK*dg_here%QX(K,J,irk)&
                       + DRK*dg_here%QX(K,J,irk-1) + BRK*dg_here%RHS_QX(K,J,irk) + ERK*dg_here%RHS_QX(K,J,1)
                  dg_here%QY(K,J,irk+1) = dg_here%QY(K,J,irk+1) + ARK*dg_here%QY(K,J,1) + CRK*dg_here%QY(K,J,irk)&
                       + DRK*dg_here%QY(K,J,irk-1) + BRK*dg_here%RHS_QY(K,J,irk) + ERK*dg_here%RHS_QY(K,J,1)


!.......Compute the transported tracer term if flagged

#ifdef TRACE
                  dg_here%iota(K,J,irk+1) = dg_here%iota(K,J,irk+1) + ARK*dg_here%iota(K,J,1) + CRK*dg_here%iota(K,J,irk)&
                 + DRK*dg_here%iota(K,J,irk-1) + BRK*dg_here%RHS_iota(K,J,irk) + ERK*dg_here%RHS_iota(K,J,1)

#endif

!......Compute chemistry transported terms if flagged

#ifdef CHEM
                  dg_here%iota(K,J,irk+1) =  dg_here%iota(K,J,irk+1) + ARK*dg_here%iota(K,J,1) + CRK*dg_here%iota(K,J,irk)&
                 + DRK*dg_here%iota(K,J,irk-1) + BRK*dg_here%RHS_iota(K,J,irk) + ERK*dg_here%RHS_iota(K,J,1)
                  dg_here%iota2(K,J,irk+1) = dg_here%iota2(K,J,irk+1) + ARK*dg_here%iota2(K,J,1) + CRK*dg_here%iota2(K,J,irk)&
                 + DRK*dg_here%iota2(K,J,irk-1) + BRK*dg_here%RHS_iota2(K,J,irk) + ERK*dg_here%RHS_iota2(K,J,1)

#endif

!.......Compute the dynamic pressure if flagged

#ifdef DYNP
                  dg_here%dynP(K,J,irk+1) = dg_here%dynP(K,J,irk+1) + ARK*dg_here%dynP(K,J,1) + CRK*dg_here%dynP(K,J,irk)&
                 + DRK*dg_here%dynP(K,J,irk-1) + BRK*dg_here%RHS_dynP(K,J,irk) + ERK*dg_here%RHS_dynP(K,J,1)

#endif

               ENDDO
            ENDDO

         endif

#endif

#ifdef CMPI
         CALL UPDATER_ELEM_MOD(dg_here%ZE,dg_here%QX,dg_here%QY,irk+1,3)
#ifdef TRACE
         CALL UPDATER_ELEM_MOD(dg_here%iota,dg_here%iota2,dg_here%QY,irk+1,2 )
#endif

#ifdef DYNP
         CALL UPDATER_ELEM_MOD(dg_here%dynP,dg_here%iota2,dg_here%QY,irk+1,2 )
#endif

#ifdef CHEM
         CALL UPDATER_ELEM_MOD(dg_here%iota,dg_here%iota2,dg_here%QY,irk+1,2 )
#endif

#ifdef SED_LAY
         do l = 1,layers
            dg_here%arrayfix => dg_here%bed(:,:,:,l) 
            CALL UPDATER_ELEM_MOD(dg_here%arrayfix,dg_here%QX,dg_here%QY,irk+1,1 )
         enddo
#endif

#endif

!.......Apply the slopelimiter if appropriate

#ifdef SLOPEALL
         CALL SLOPELIMITER(s,dg_here)
#endif

#ifdef SLOPE5
         CALL SLOPELIMITER(s,dg_here) 
#endif

#ifdef STBLZR
         if (.not.dg_here%stblzr) then
            CALL SLOPELIMITER(s,dg_here)
            dg_here%stblzr = .true.
         endif
#endif


!.......Apply the wet-dry algorithm if appropriate

      IF (NOLIFA .GE. 2) THEN
         CALL WETDRY(dg_here)
      ENDIF

!$$$C.......Apply the slopelimiter again if auxiliary variable is being used
!$$$
!$$$         if (EVMSUM.NE.0.D0.or.dg_here%artdif.eq.1) then
!$$$#ifdef SLOPEALL
!$$$            CALL SLOPELIMITER()
!$$$#endif
!$$$            
!$$$#ifdef SLOPE5
!$$$            CALL SLOPELIMITER() 
!$$$#endif
!$$$         endif

!.......For parallel run update for dg_here%dofs if slope limiter and/or wetting and
!.......drying is being used

#ifdef CMPI
         if (dg_here%SLOPEFLAG.NE.0.OR.NOLIFA.GE.2) THEN
            CALL UPDATER_ELEM_MOD(dg_here%ZE,dg_here%QX,dg_here%QY,irk+1,3)
            
#ifdef TRACE
            CALL UPDATER_ELEM_MOD(dg_here%iota,dg_here%iota2,dg_here%QY,irk+1,2 )
#endif

#ifdef CHEM
            CALL UPDATER_ELEM_MOD(dg_here%iota,dg_here%iota2,dg_here%QY,irk+1,2 )
#endif

#ifdef DYNP
            CALL UPDATER_ELEM_MOD(dg_here%dynP,dg_here%iota2,dg_here%QY,irk+1,2 )
#endif


#ifdef SED_LAY
            do l = 1,layers
               dg_here%arrayfix => dg_here%bed(:,:,:,l)
               CALL UPDATER_ELEM_MOD(dg_here%arrayfix,dg_here%QX,dg_here%QY,irk+1,1 )
            enddo
#endif
         endif
#endif

!.......p_enrich soln and update

#ifdef P_AD
         if (irk+1.eq.2) then
            CALL p_enrichment(s,dg_here,it,0)
         endif
#ifdef CMPI

         CALL UPDATER_ELEM_MOD(dg_here%ZE,dg_here%QX,dg_here%QY,irk+1,3)

#ifdef TRACE
         CALL UPDATER_ELEM_MOD(dg_here%iota,dg_here%iota2,dg_here%QY,irk+1,2 )
#endif

#ifdef CHEM
         CALL UPDATER_ELEM_MOD(dg_here%iota,dg_here%iota2,dg_here%QY,irk+1,2 )
#endif

#ifdef DYNP
         CALL UPDATER_ELEM_MOD(dg_here%dynP,dg_here%iota2,dg_here%QY,irk+1,2 )
#endif

#ifdef SED_LAY
         do l = 1,layers
            dg_here%arrayfix => dg_here%bed(:,:,:,l)
            CALL UPDATER_ELEM_MOD(dg_here%arrayfix,dg_here%QX,dg_here%QY,irk+1,1 )
         enddo
#endif

#endif 

#endif

 100  CONTINUE

!.....Update variables for next time step

#ifdef FILTER 
      !Stabilization filtering 
      DO J = 1,NE
         DO K = 2,dg_here%DOFS(J)

            alph = 10.D0
            s_dg = dg_here%RK_stage

            if (k.le.3) then
               eta = 1.D0 / ( pdg_el(j)+1.D0 )
            elseif (k.gt.3.and.k.le.6) then
               eta = 2.D0 / ( pdg_el(j)+1.D0 )
            elseif (k.gt.6.and.k.le.10) then
               eta = 3.D0 / ( pdg_el(j)+1.D0 )
            elseif (k.gt.11.and.k.le.15) then
               eta = 4.D0 / ( pdg_el(j)+1.D0 )
            elseif (k.gt.16.and.k.le.21) then
               eta = 5.D0 / ( pdg_el(j)+1.D0 )       
            elseif (k.gt.22.and.k.le.28) then
               eta = 6.D0 / ( pdg_el(j)+1.D0 )   
            elseif (k.gt.29.and.k.le.36) then
               eta = 7.D0 / ( pdg_el(j)+1.D0 )               
            endif

            seta = eta**s_dg
            sigma = exp(-alph*seta)

!.......Filter the layers if flagged

#ifdef SED_LAY
            do l = 1,layers
               
               dg_here%bed(K,J,dg_here%NRK+1,l) = sigma * dg_here%bed(K,J,dg_here%nrk+1,l) 
               
            enddo
#endif

            dg_here%ZE(K,J,dg_here%NRK+1) = sigma*dg_here%ZE(K,J,dg_here%nrk+1) 
            dg_here%QX(K,J,dg_here%NRK+1) = sigma*dg_here%QX(K,J,dg_here%nrk+1)
            dg_here%QY(K,J,dg_here%NRK+1) = sigma*dg_here%QY(K,J,dg_here%nrk+1)


!.......Filter the transported tracer term if flagged

#ifdef TRACE
            dg_here%iota(K,J,dg_here%NRK+1) = sigma*dg_here%iota(K,J,dg_here%nrk+1)

#endif

!......Filter chemistry transported terms if flagged

#ifdef CHEM
            dg_here%iota(K,J,dg_here%NRK+1) =  sigma*dg_here%iota(K,J,dg_here%nrk+1) 
            dg_here%iota2(K,J,dg_here%NRK+1) = sigma*dg_here%iota2(K,J,dg_here%nrk+1) 

#endif

         ENDDO
      ENDDO
#endif
      

      DO J = 1,S%MNE
         
         DO K = 1,dg_here%DOFS(J)
            dg_here%ZE(K,J,1) = dg_here%ZE(K,J,dg_here%NRK+1)
            dg_here%QX(K,J,1) = dg_here%QX(K,J,dg_here%NRK+1)
            dg_here%QY(K,J,1) = dg_here%QY(K,J,dg_here%NRK+1)
#ifdef TRACE
            dg_here%iota(K,J,1) = dg_here%iota(K,J,dg_here%NRK+1)
#endif

#ifdef CHEM
            dg_here%iota(K,J,1) = dg_here%iota(K,J,dg_here%NRK+1)
            dg_here%iota2(K,J,1) = dg_here%iota2(K,J,dg_here%NRK+1)
#endif

#ifdef DYNP
            dg_here%dynP(K,J,1) = dg_here%dynP(K,J,dg_here%NRK+1)
#endif

#ifdef SED_LAY
            dg_here%bed(K,J,1,:) = dg_here%bed(K,J,dg_here%NRK+1,:)
#endif
         ENDDO
      ENDDO

      DO KK = 2,dg_here%NRK+1
         DO J = 1,S%MNE
            DO K = 1,dg_here%DOFH
               dg_here%ZE(K,J,KK) = 0.D0
               dg_here%QX(K,J,KK) = 0.D0
               dg_here%QY(K,J,KK) = 0.D0
               dg_here%RHS_ZE(K,J,KK-1) = 0.D0
               dg_here%RHS_QX(K,J,KK-1) = 0.D0
               dg_here%RHS_QY(K,J,KK-1) = 0.D0
#ifdef TRACE
               dg_here%iota(K,J,KK) = 0.D0
               dg_here%RHS_iota(K,J,KK-1) = 0.D0
#endif

#ifdef CHEM
               dg_here%iota(K,J,KK) = 0.D0
               dg_here%iota2(K,J,KK) = 0.D0
               dg_here%RHS_iota(K,J,KK-1) = 0.D0
               dg_here%RHS_iota2(K,J,KK-1) = 0.D0
#endif

#ifdef DYNP
               dg_here%dynP(K,J,KK) = 0.D0
               dg_here%RHS_dynP(K,J,KK-1) = 0.D0
#endif

#ifdef SED_LAY
               dg_here%bed(K,J,KK,:) = 0.D0
               dg_here%hb(K,J,KK) = 0.D0
               dg_here%RHS_bed(K,J,KK-1,:) = 0.D0
#endif


            ENDDO
         ENDDO
      ENDDO

      RETURN
      END SUBROUTINE
