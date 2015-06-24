
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

      SUBROUTINE DG_HYDRO_TIMESTEP(s,dg,IT)

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
      type (dg_type) :: dg

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

      DO 100 dg%IRK = 1,dg%NRK
         
!.......Compute the DG time and DG ramp

#ifdef RKSSP

 
         dg%TIMEDG = TIME_A - DTDP + dg%DTVD(dg%IRK)*DTDP
         dg%RAMPDG = 1.D0
         RAMPExtFlux = 1.0d0
         IF (NRAMP.GE.1) THEN
            IF (NRAMP.EQ.1) THEN
               dg%RAMPDG = TANH((2.D0*((IT-1) + dg%DTVD(dg%IRK))*DTDP/86400.D0)/DRAMP)
               RAMPExtFlux = &
                   TANH((2.D0*((IT-1) + dg%DTVD(dg%IRK))*DTDP/86400.D0)/DRAMPExtFlux)
            ENDIF
            IF (NRAMP.EQ.2) THEN
               dg%RAMPDG = TANH((2.D0*((IT-1) + dg%DTVD(dg%IRK))*DTDP/86400.D0)/DRAMP)
               RAMPExtFlux = &
                   TANH((2.D0*((IT-1) + dg%DTVD(dg%IRK))*DTDP/86400.D0)/DRAMPExtFlux)
            ENDIF
            IF (NRAMP.EQ.3) THEN
               WRITE(*,*) 'NRAMP = 3 not supported '
               STOP
            ENDIF
         ENDIF
         RAMP = dg%RAMPDG
         
!.......Obtain the meteorological forcing

         IF (NWS.NE.0) CALL MET_FORCING(s,dg,IT)
        
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
            IF((NRS.EQ.3).AND.(dg%IRK.EQ.1)) THEN
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

         IF (NTIP.NE.0) CALL TIDAL_POTENTIAL(dg)

!.......Compute LDG auxiliary equations
         
         IF (EVMSUM.NE.0.D0.or.dg%artdif.eq.1) CALL LDG_HYDRO(s,dg,IT)
         
!.......Compute elevation specified edges

         IF (dg%NEEDS.GT.0)  CALL OCEAN_EDGE_HYDRO(s,dg,IT)

!.......Compute no-normal flow edges

         IF (dg%NLEDS.GT.0)  CALL LAND_EDGE_HYDRO(s,dg,IT)

!.......Compute non-zero flow edges

         IF (dg%NFEDS.GT.0)  CALL FLOW_EDGE_HYDRO(s,dg,IT)
         
!.......Compute radiation edges

         IF (dg%NREDS.GT.0)  CALL RADIATION_EDGE_HYDRO(s,dg,IT)

!.......Compute internal barrier edges

         IF (dg%NIBEDS.GT.0) CALL IBARRIER_EDGE_HYDRO(s,dg,IT)
         
!.......Compute external barrier edges

         IF (dg%NIBEDS.GT.0) CALL EBARRIER_EDGE_HYDRO(s,dg,IT)
         
!.......Compute internal edges

         CALL INTERNAL_EDGE_HYDRO(s,dg,IT)

!.......Compute elements to finish building the rhs
         
         CALL RHS_DG_HYDRO(s,dg)

!.......SSP Runge-Kutta Time Scheme

         DO I = 1,dg%IRK
            ARK = dg%ATVD(dg%IRK,I)
            BRK = dg%BTVD(dg%IRK,I)*DT
            DO J = 1,NE
               DO K = 1,dg%DOFS(J)

#ifdef SED_LAY

                  do l = 1,layers
                     
                     dg%bed(K,J,dg%IRK+1,l) = dg%bed(K,J,dg%irk+1,l) + ARK*dg%bed(K,J,I,l)&
                         + BRK*dg%RHS_bed(K,J,I,l)

                     !adjust for state variable represenation
                     dg%ZE(K,J,dg%IRK+1) = dg%ZE(K,J,dg%IRK+1) - BRK*dg%RHS_bed(K,J,I,l)


                  enddo
#endif

                  dg%ZE(K,J,dg%IRK+1) = dg%ZE(K,J,dg%irk+1) + ARK*dg%ZE(K,J,I)&
                      + BRK*dg%RHS_ZE(K,J,I) 
                  dg%QX(K,J,dg%IRK+1) = dg%QX(K,J,dg%irk+1) + ARK*dg%QX(K,J,I)&
                      + BRK*dg%RHS_QX(K,J,I)
                  dg%QY(K,J,dg%IRK+1) = dg%QY(K,J,dg%irk+1) + dg%ATVD(dg%IRK,I)*dg%QY(K,J,I)&
                      + BRK*dg%RHS_QY(K,J,I)


!.......Compute the transported tracer term if flagged

#ifdef TRACE
                  dg%iota(K,J,dg%IRK+1) = dg%iota(K,J,dg%irk+1) + ARK*dg%iota(K,J,I)&
                      + BRK*dg%RHS_iota(K,J,I)
#endif

!......Compute chemistry transported terms if flagged

#ifdef CHEM
                  dg%iota(K,J,dg%IRK+1) = dg%iota(K,J,dg%irk+1) + ARK*dg%iota(K,J,I)&
                      + BRK*dg%RHS_iota(K,J,I)
                  dg%iota2(K,J,dg%IRK+1) = dg%iota2(K,J,dg%irk+1) + ARK*dg%iota2(K,J,I)&
                      + BRK*dg%RHS_iota2(K,J,I)
#endif

!.......Compute the dynamic pressure if flagged

#ifdef DYNP
                  dg%dynP(K,J,dg%IRK+1) = dg%dynP(K,J,dg%irk+1) + ARK*dg%dynP(K,J,I)&
                      + BRK*dg%RHS_dynP(K,J,I)
#endif


               ENDDO
            ENDDO
         ENDDO

#endif

#ifdef RKC

         print*,"Using RKC"

         dg%TIMEDG = TIME_A - DTDP + dg%RKC_c(dg%IRK)*DTDP
         dg%RAMPDG = 1.D0
         RAMPExtFlux = 1.0d0
         IF (NRAMP.GE.1) THEN
            IF (NRAMP.EQ.1) THEN
               dg%RAMPDG = TANH((2.D0*(IT + dg%RKC_c(dg%IRK))*DTDP/86400.D0)/DRAMP)
               RAMPExtFlux = &
                    TANH((2.D0*(IT + dg%RKC_c(dg%IRK))*DTDP/86400.D0)/DRAMPExtFlux)
            ENDIF
            IF (NRAMP.EQ.2) THEN
               dg%RAMPDG = TANH((2.D0*(IT + dg%RKC_c(dg%IRK))*DTDP/86400.D0)/DRAMP)
               RAMPExtFlux = &
                    TANH((2.D0*(IT + dg%RKC_c(dg%IRK))*DTDP/86400.D0)/DRAMPExtFlux)
            ENDIF
            IF (NRAMP.EQ.3) THEN
               WRITE(*,*) 'NRAMP = 3 not supported '
               STOP
            ENDIF
         ENDIF
         RAMP = dg%RAMPDG
         
!.......Obtain the meteorological forcing

         IF (NWS.NE.0) CALL MET_FORCING(s,dg,IT)
         
!.......Compute tidal potential terms

         IF (NTIP.NE.0) CALL TIDAL_POTENTIAL(dg)

!.......Compute LDG auxiliary equations
         
         IF (EVMSUM.NE.0.D0.or.dg%artdif.eq.1) CALL LDG_HYDRO(s,dg,IT)
         
!.......Compute elevation specified edges

         IF (dg%NEEDS.GT.0)  CALL OCEAN_EDGE_HYDRO(s,dg,IT)

!.......Compute no-normal flow edges

         IF (dg%NLEDS.GT.0)  CALL LAND_EDGE_HYDRO(s,dg,IT)

!.......Compute non-zero flow edges

         IF (dg%NFEDS.GT.0)  CALL FLOW_EDGE_HYDRO(s,dg,IT)
         
!.......Compute radiation edges

         IF (dg%NREDS.GT.0)  CALL RADIATION_EDGE_HYDRO(s,dg,IT)

!.......Compute internal barrier edges

         IF (dg%NIBEDS.GT.0) CALL IBARRIER_EDGE_HYDRO(s,dg,IT)
         
!.......Compute external barrier edges

         IF (dg%NIBEDS.GT.0) CALL EBARRIER_EDGE_HYDRO(s,dg,IT)
         
!.......Compute internal edges

         CALL INTERNAL_EDGE_HYDRO(s,dg,IT)

!.......Compute elements to finish building the rhs
         
         CALL RHS_DG_HYDRO(s,dg)

!.......RKC Time Scheme

         if (dg%irk.eq.1) then

            BRK = dg%RKC_tildemu(1)*dt

            DO J = 1,NE
               DO K = 1,dg%DOFS(J)

#ifdef SED_LAY
                  do l = 1,layers
                     
                     dg%bed(K,J,dg%IRK+1,l) = dg%bed(K,J,dg%irk+1,l) + dg%bed(K,J,1,l) &
                         + BRK*dg%RHS_bed(K,J,1,l)

                   !adjust for state variable represenation
                     dg%ZE(K,J,dg%IRK+1) = dg%ZE(K,J,dg%IRK+1) - BRK*dg%RHS_bed(K,J,1,l)
                     
                  enddo
#endif
                  
                  dg%ZE(K,J,dg%IRK+1) = dg%ZE(K,J,dg%irk+1) + dg%ZE(K,J,1)&
                      + BRK*dg%RHS_ZE(K,J,1)
                  dg%QX(K,J,dg%IRK+1) = dg%QX(K,J,dg%irk+1) + dg%QX(K,J,1)&
                      + BRK*dg%RHS_QX(K,J,1)
                  dg%QY(K,J,dg%IRK+1) = dg%QY(K,J,dg%irk+1)  + dg%QY(K,J,1)&
                      + BRK*dg%RHS_QY(K,J,1)
                  
!.......Compute the transported tracer term if flagged
                  
#ifdef TRACE
                  dg%iota(K,J,dg%IRK+1) = dg%iota(K,J,dg%irk+1) + dg%iota(K,J,1)&
                      + BRK*dg%RHS_iota(K,J,1)
#endif
                  
!......Compute chemistry transported terms if flagged
                  
#ifdef CHEM
                  dg%iota(K,J,dg%IRK+1) = dg%iota(K,J,dg%irk+1) + dg%iota(K,J,1)&
                       + BRK*dg%RHS_iota(K,J,1)
                  dg%iota2(K,J,dg%IRK+1) = dg%iota2(K,J,dg%irk+1) + dg%iota2(K,J,1)&
                       + BRK*dg%RHS_iota2(K,J,1)
#endif

!.......Compute the dynamic pressure if flagged
                  
#ifdef DYNP
                  dg%dynP(K,J,dg%IRK+1) = dg%dynP(K,J,dg%irk+1) + dg%dynP(K,J,1)&
                      + BRK*dg%RHS_dynP(K,J,1)
#endif
               ENDDO
            ENDDO
            
         else

            ARK = 1.D0 - dg%RKC_mu(dg%irk)-dg%RKC_nu(dg%irk)
            CRK = dg%RKC_mu(dg%irk)
            DRK = dg%RKC_nu(dg%irk)

            BRK = dg%RKC_tildemu(dg%irk)*dt
            ERK = dg%RKC_gamma(dg%irk)*dt

             DO J = 1,NE
               DO K = 1,dg%DOFS(J)

#ifdef SED_LAY
                  do l = 1,layers
                     
                     dg%bed(K,J,dg%IRK+1,l) = dg%bed(K,J,dg%irk+1,l) + ARK*dg%bed(K,J,1,l) + CRK*dg%bed(K,J,dg%irk,l)&
                         + DRK*dg%bed(K,J,dg%irk-1,l) + BRK*dg%RHS_bed(K,J,dg%irk,l) + ERK*dg%RHS_bed(K,J,1,l)


                   !adjust for state variable represenation
                     dg%ZE(K,J,dg%IRK+1) = dg%ZE(K,J,dg%IRK+1) - BRK*dg%RHS_bed(K,J,dg%irk,l) - ERK*dg%RHS_bed(K,J,1,l)
                     
                  enddo
#endif

                  dg%ZE(K,J,dg%IRK+1) = dg%ZE(K,J,dg%irk+1) + ARK*dg%ZE(K,J,1) + CRK*dg%ZE(K,J,dg%irk)&
                       + DRK*dg%ZE(K,J,dg%irk-1) + BRK*dg%RHS_ZE(K,J,dg%irk) + ERK*dg%RHS_ZE(K,J,1)
                  dg%QX(K,J,dg%IRK+1) = dg%QX(K,J,dg%irk+1) + ARK*dg%QX(K,J,1) + CRK*dg%QX(K,J,dg%irk)&
                       + DRK*dg%QX(K,J,dg%irk-1) + BRK*dg%RHS_QX(K,J,dg%irk) + ERK*dg%RHS_QX(K,J,1)
                  dg%QY(K,J,dg%IRK+1) = dg%QY(K,J,dg%irk+1) + ARK*dg%QY(K,J,1) + CRK*dg%QY(K,J,dg%irk)&
                       + DRK*dg%QY(K,J,dg%irk-1) + BRK*dg%RHS_QY(K,J,dg%irk) + ERK*dg%RHS_QY(K,J,1)


!.......Compute the transported tracer term if flagged

#ifdef TRACE
                  dg%iota(K,J,dg%IRK+1) = dg%iota(K,J,dg%irk+1) + ARK*dg%iota(K,J,1) + CRK*dg%iota(K,J,dg%irk)&
                 + DRK*dg%iota(K,J,dg%irk-1) + BRK*dg%RHS_iota(K,J,dg%irk) + ERK*dg%RHS_iota(K,J,1)

#endif

!......Compute chemistry transported terms if flagged

#ifdef CHEM
                  dg%iota(K,J,dg%IRK+1) =  dg%iota(K,J,dg%irk+1) + ARK*dg%iota(K,J,1) + CRK*dg%iota(K,J,dg%irk)&
                 + DRK*dg%iota(K,J,dg%irk-1) + BRK*dg%RHS_iota(K,J,dg%irk) + ERK*dg%RHS_iota(K,J,1)
                  dg%iota2(K,J,dg%IRK+1) = dg%iota2(K,J,dg%irk+1) + ARK*dg%iota2(K,J,1) + CRK*dg%iota2(K,J,dg%irk)&
                 + DRK*dg%iota2(K,J,dg%irk-1) + BRK*dg%RHS_iota2(K,J,dg%irk) + ERK*dg%RHS_iota2(K,J,1)

#endif

!.......Compute the dynamic pressure if flagged

#ifdef DYNP
                  dg%dynP(K,J,dg%IRK+1) = dg%dynP(K,J,dg%irk+1) + ARK*dg%dynP(K,J,1) + CRK*dg%dynP(K,J,dg%irk)&
                 + DRK*dg%dynP(K,J,dg%irk-1) + BRK*dg%RHS_dynP(K,J,dg%irk) + ERK*dg%RHS_dynP(K,J,1)

#endif

               ENDDO
            ENDDO

         endif

#endif

#ifdef CMPI
         CALL UPDATER_ELEM_MOD(dg%ZE,dg%QX,dg%QY,dg%IRK+1,3)
#ifdef TRACE
         CALL UPDATER_ELEM_MOD(dg%iota,dg%iota2,dg%QY,dg%IRK+1,2 )
#endif

#ifdef DYNP
         CALL UPDATER_ELEM_MOD(dg%dynP,dg%iota2,dg%QY,dg%IRK+1,2 )
#endif

#ifdef CHEM
         CALL UPDATER_ELEM_MOD(dg%iota,dg%iota2,dg%QY,dg%IRK+1,2 )
#endif

#ifdef SED_LAY
         do l = 1,layers
            dg%arrayfix => dg%bed(:,:,:,l) 
            CALL UPDATER_ELEM_MOD(dg%arrayfix,dg%QX,dg%QY,dg%IRK+1,1 )
         enddo
#endif

#endif

!.......Apply the slopelimiter if appropriate

#ifdef SLOPEALL
         CALL SLOPELIMITER(s,dg)
#endif

#ifdef SLOPE5
         CALL SLOPELIMITER(s,dg) 
#endif

#ifdef STBLZR
         if (.not.dg%stblzr) then
            CALL SLOPELIMITER(s,dg)
            dg%stblzr = .true.
         endif
#endif


!.......Apply the wet-dry algorithm if appropriate

      IF (NOLIFA .GE. 2) THEN
         CALL WETDRY(dg)
      ENDIF

!$$$C.......Apply the slopelimiter again if auxiliary variable is being used
!$$$
!$$$         if (EVMSUM.NE.0.D0.or.dg%artdif.eq.1) then
!$$$#ifdef SLOPEALL
!$$$            CALL SLOPELIMITER()
!$$$#endif
!$$$            
!$$$#ifdef SLOPE5
!$$$            CALL SLOPELIMITER() 
!$$$#endif
!$$$         endif

!.......For parallel run update for dg%dofs if slope limiter and/or wetting and
!.......drying is being used

#ifdef CMPI
         if (dg%SLOPEFLAG.NE.0.OR.NOLIFA.GE.2) THEN
            CALL UPDATER_ELEM_MOD(dg%ZE,dg%QX,dg%QY,dg%IRK+1,3)
            
#ifdef TRACE
            CALL UPDATER_ELEM_MOD(dg%iota,dg%iota2,dg%QY,dg%IRK+1,2 )
#endif

#ifdef CHEM
            CALL UPDATER_ELEM_MOD(dg%iota,dg%iota2,dg%QY,dg%IRK+1,2 )
#endif

#ifdef DYNP
            CALL UPDATER_ELEM_MOD(dg%dynP,dg%iota2,dg%QY,dg%IRK+1,2 )
#endif


#ifdef SED_LAY
            do l = 1,layers
               dg%arrayfix => dg%bed(:,:,:,l)
               CALL UPDATER_ELEM_MOD(dg%arrayfix,dg%QX,dg%QY,dg%IRK+1,1 )
            enddo
#endif
         endif
#endif

!.......p_enrich soln and update

#ifdef P_AD
         if (dg%irk+1.eq.2) then
            CALL p_enrichment(s,dg,it,0)
         endif
#ifdef CMPI

         CALL UPDATER_ELEM_MOD(dg%ZE,dg%QX,dg%QY,dg%IRK+1,3)

#ifdef TRACE
         CALL UPDATER_ELEM_MOD(dg%iota,dg%iota2,dg%QY,dg%IRK+1,2 )
#endif

#ifdef CHEM
         CALL UPDATER_ELEM_MOD(dg%iota,dg%iota2,dg%QY,dg%IRK+1,2 )
#endif

#ifdef DYNP
         CALL UPDATER_ELEM_MOD(dg%dynP,dg%iota2,dg%QY,dg%IRK+1,2 )
#endif

#ifdef SED_LAY
         do l = 1,layers
            dg%arrayfix => dg%bed(:,:,:,l)
            CALL UPDATER_ELEM_MOD(dg%arrayfix,dg%QX,dg%QY,dg%IRK+1,1 )
         enddo
#endif

#endif 

#endif

 100  CONTINUE

!.....Update variables for next time step

#ifdef FILTER 
      !Stabilization filtering 
      DO J = 1,NE
         DO K = 2,dg%DOFS(J)

            alph = 10.D0
            s_dg = dg%RK_stage

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
               
               dg%bed(K,J,dg%NRK+1,l) = sigma * dg%bed(K,J,dg%nrk+1,l) 
               
            enddo
#endif

            dg%ZE(K,J,dg%NRK+1) = sigma*dg%ZE(K,J,dg%nrk+1) 
            dg%QX(K,J,dg%NRK+1) = sigma*dg%QX(K,J,dg%nrk+1)
            dg%QY(K,J,dg%NRK+1) = sigma*dg%QY(K,J,dg%nrk+1)


!.......Filter the transported tracer term if flagged

#ifdef TRACE
            dg%iota(K,J,dg%NRK+1) = sigma*dg%iota(K,J,dg%nrk+1)

#endif

!......Filter chemistry transported terms if flagged

#ifdef CHEM
            dg%iota(K,J,dg%NRK+1) =  sigma*dg%iota(K,J,dg%nrk+1) 
            dg%iota2(K,J,dg%NRK+1) = sigma*dg%iota2(K,J,dg%nrk+1) 

#endif

         ENDDO
      ENDDO
#endif
      

      DO J = 1,S%MNE
         
         DO K = 1,dg%DOFS(J)
            dg%ZE(K,J,1) = dg%ZE(K,J,dg%NRK+1)
            dg%QX(K,J,1) = dg%QX(K,J,dg%NRK+1)
            dg%QY(K,J,1) = dg%QY(K,J,dg%NRK+1)
#ifdef TRACE
            dg%iota(K,J,1) = dg%iota(K,J,dg%NRK+1)
#endif

#ifdef CHEM
            dg%iota(K,J,1) = dg%iota(K,J,dg%NRK+1)
            dg%iota2(K,J,1) = dg%iota2(K,J,dg%NRK+1)
#endif

#ifdef DYNP
            dg%dynP(K,J,1) = dg%dynP(K,J,dg%NRK+1)
#endif

#ifdef SED_LAY
            dg%bed(K,J,1,:) = dg%bed(K,J,dg%NRK+1,:)
#endif
         ENDDO
      ENDDO

      DO KK = 2,dg%NRK+1
         DO J = 1,S%MNE
            DO K = 1,dg%DOFH
               dg%ZE(K,J,KK) = 0.D0
               dg%QX(K,J,KK) = 0.D0
               dg%QY(K,J,KK) = 0.D0
               dg%RHS_ZE(K,J,KK-1) = 0.D0
               dg%RHS_QX(K,J,KK-1) = 0.D0
               dg%RHS_QY(K,J,KK-1) = 0.D0
#ifdef TRACE
               dg%iota(K,J,KK) = 0.D0
               dg%RHS_iota(K,J,KK-1) = 0.D0
#endif

#ifdef CHEM
               dg%iota(K,J,KK) = 0.D0
               dg%iota2(K,J,KK) = 0.D0
               dg%RHS_iota(K,J,KK-1) = 0.D0
               dg%RHS_iota2(K,J,KK-1) = 0.D0
#endif

#ifdef DYNP
               dg%dynP(K,J,KK) = 0.D0
               dg%RHS_dynP(K,J,KK-1) = 0.D0
#endif

#ifdef SED_LAY
               dg%bed(K,J,KK,:) = 0.D0
               dg%hb(K,J,KK) = 0.D0
               dg%RHS_bed(K,J,KK-1,:) = 0.D0
#endif


            ENDDO
         ENDDO
      ENDDO

      RETURN
      END SUBROUTINE
