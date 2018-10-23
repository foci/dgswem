
C***********************************************************************
C     
C     SUBROUTINE DG_HYDRO_TIMESTEP(IT)
C     
C     This is the main subroutine for the DG hydro
C     
C     Written by Ethan Kubatko (06-11-2004)
C     
C-----------------------------------------------------------------------
C     
C     Modification history on hp_DG_ADCIRC since v9_sb1
C     
C     v??     - 08/05    - sb - parallelized
C     v??     - 08/05    - sb - wetting/drying
C     v10_sb5 - 10/05    - sb - wetting/drying algorithm is improved
C     v10_sb5 - 10/05    - sb - Ethan's slope limiter is consolidated
C     v21     - 11/11    - cem - major update/p_enrich/slope/rkc etc
C     v22     -  6/12    - cem - major update/sediment/diffusion/fluxes
C     Moving to versioning control
C     
C***********************************************************************

      SUBROUTINE DG_HYDRO_TIMESTEP(IT)

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

      DO 100 IRK = 1,NRK
         
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
         
C.......Obtain the meteorological forcing

         IF (NWS.NE.0) CALL MET_FORCING(IT)
        
        IF(NRS.GE.1) THEN
          IF(TIME_A.GT.RSTIME2) THEN
            RSTIME1=RSTIME2
            RSTIME2=RSTIME2+RSTIMINC
           !DO I=1,NP
           !  RSNX1(I)=RSNX2(I)
           !  RSNY1(I)=RSNY2(I)
           !END DO
#ifdef SWAN
Casey 090302: Added for coupling to SWAN.
            IF((NRS.EQ.3).AND.(IRK.EQ.1)) THEN
              InterpoWeight = 1.0
              CALL ComputeWaveDrivenForces
Casey 090707: We want to extrapolate forward in time.  Load the latest 
C             (current) forces into RSNX1/RSNY1, and then load the 
C             future forces into RSNX2/RSNY2.
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
Casey 090302: Added these lines for output to the rads.64 file.
            RSNXOUT(I) = RSX
            RSNYOUT(I) = RSY
#endif
         ENDDO
      ENDIF
 
C.......Compute tidal potential terms

         IF (NTIP.NE.0) CALL TIDAL_POTENTIAL()

C.......For non-periodic elevation bcs

         IF ((NBFR.EQ.0).AND.(NOPE.GT.0)) THEN
           IF (TIME_A.GT.ETIME2) THEN
             ETIME1 = ETIME2
             ETIME2 = ETIME1 + ETIMINC
             DO J = 1,NETA
               ESBIN1(J) = ESBIN2(J)
               READ(19,*) ESBIN2(J)
             ENDDO
           ENDIF
           ETRATIO = (TIMEDG - ETIME1)/ETIMINC
           DO I = 1,NETA
             NBDI = NBD(I)
             ETA2(NBDI) = RAMPDG*( ESBIN1(I) + ETRATIO*( ESBIN2(I) -
     &                                                   ESBIN1(I) ) )
           ENDDO
         ENDIF
         
C.......For non-periodic flux bcs

         IF (MNFFR.EQ.1 .AND. NFFR.EQ.0) THEN
           IF (TIME_A.GT.QTIME2) THEN
             QTIME1 = QTIME2
             QTIME2 = QTIME1 + FTIMINC
             DO J=1,NVEL
               QNIN1(J)=QNIN2(J)
               QNIN2(J)=0.D0
               IF((LBCODEI(J).EQ.2).OR.(LBCODEI(J).EQ.12)
     &                             .OR.(LBCODEI(J).EQ.22))    
     &         READ(20,*) QNIN2(J)
             END DO
           ENDIF
           QTRATIO = (TIMEDG - QTIME1)/FTIMINC
           NQEDS = 0
           DO I = 1,NVEL
             IF((LBCODEI(I).EQ.2).OR.(LBCODEI(I).EQ.12)
     &                             .OR.(LBCODEI(I).EQ.22)) THEN 
               NQEDS = NQEDS + 1
               IF (NQEDS .LE. NFEDS) THEN
                 QNAM_DG(1,NQEDS,1) = RAMPDG*( QNIN1(I) + QTRATIO*
     &                                         ( QNIN2(I) - QNIN1(I) ) )
                 QNPH_DG(1,NQEDS,1) = 0.D0
                 QNAM_DG(1,NQEDS,2) = RAMPDG*( QNIN1(I+1) + QTRATIO*
     &                                     ( QNIN2(I+1) - QNIN1(I+1) ) )
                 QNPH_DG(1,NQEDS,2) = 0.D0
               ENDIF
             ENDIF
           ENDDO
         ENDIF       

C.......Compute LDG auxiliary equations
         
         IF (EVMSUM.NE.0.D0.or.artdif.eq.1) CALL LDG_HYDRO(IT)
         
C.......Compute elevation specified edges

         IF (NEEDS.GT.0)  CALL OCEAN_EDGE_HYDRO(IT)

C.......Compute no-normal flow edges

         IF (NLEDS.GT.0)  CALL LAND_EDGE_HYDRO(IT)

C.......Compute non-zero flow edges

         IF (NFEDS.GT.0)  CALL FLOW_EDGE_HYDRO(IT)
         
C.......Compute radiation edges

         IF (NREDS.GT.0)  CALL RADIATION_EDGE_HYDRO(IT)

C.......Compute internal barrier edges

         IF (NIBEDS.GT.0) CALL IBARRIER_EDGE_HYDRO(IT)
         
C.......Compute external barrier edges

         IF (NEBEDS.GT.0) CALL EBARRIER_EDGE_HYDRO(IT)
         
C.......Compute internal edges

         CALL INTERNAL_EDGE_HYDRO(IT)

C.......Compute elements to finish building the rhs
         
         CALL RHS_DG_HYDRO()

C.......SSP Runge-Kutta Time Scheme

         DO I = 1,IRK
            ARK = ATVD(IRK,I)
            BRK = BTVD(IRK,I)*DT
            DO J = 1,NE
               DO K = 1,DOFS(J)

#ifdef SED_LAY

                  do l = 1,layers
                     
                     bed(K,J,IRK+1,l) = bed(K,J,irk+1,l) 
     &                  + ARK*bed(K,J,I,l) + BRK*RHS_bed(K,J,I,l)

                     !adjust for state variable represenation
                     ZE(K,J,IRK+1) = ZE(K,J,IRK+1)-BRK*RHS_bed(K,J,I,l)


                  enddo
#endif

                  ZE(K,J,IRK+1) = ZE(K,J,irk+1) + ARK*ZE(K,J,I)
     &                 + BRK*RHS_ZE(K,J,I) 
                  QX(K,J,IRK+1) = QX(K,J,irk+1) + ARK*QX(K,J,I)
     &                 + BRK*RHS_QX(K,J,I)
                  QY(K,J,IRK+1) = QY(K,J,irk+1) + ATVD(IRK,I)*QY(K,J,I)
     &                 + BRK*RHS_QY(K,J,I)


C.......Compute the transported tracer term if flagged

#ifdef TRACE
                  iota(K,J,IRK+1) = iota(K,J,irk+1) + ARK*iota(K,J,I)
     &                 + BRK*RHS_iota(K,J,I)
#endif

C......Compute chemistry transported terms if flagged

#ifdef CHEM
                  iota(K,J,IRK+1) = iota(K,J,irk+1) + ARK*iota(K,J,I)
     &                 + BRK*RHS_iota(K,J,I)
                  iota2(K,J,IRK+1) = iota2(K,J,irk+1) + ARK*iota2(K,J,I)
     &                 + BRK*RHS_iota2(K,J,I)
#endif

C.......Compute the dynamic pressure if flagged

#ifdef DYNP
                  dynP(K,J,IRK+1) = dynP(K,J,irk+1) + ARK*dynP(K,J,I)
     &                 + BRK*RHS_dynP(K,J,I)
#endif


               ENDDO
            ENDDO
         ENDDO

#endif

#ifdef RKC

         TIMEDG = TIME_A - DTDP + RKC_c(IRK)*DTDP
         RAMPDG = 1.D0
         RAMPExtFlux = 1.0d0
         IF (NRAMP.GE.1) THEN
            IF (NRAMP.EQ.1) THEN
               RAMPDG = TANH((2.D0*(IT + RKC_c(IRK))*DTDP/86400.D0)
     &            /DRAMP)
               RAMPExtFlux = 
     &              TANH((2.D0*(IT + RKC_c(IRK))*DTDP/86400.D0)
     &            /DRAMPExtFlux)
            ENDIF
            IF (NRAMP.EQ.2) THEN
               RAMPDG = TANH((2.D0*(IT + RKC_c(IRK))*DTDP/86400.D0)
     &            /DRAMP)
               RAMPExtFlux = 
     &              TANH((2.D0*(IT + RKC_c(IRK))*DTDP/86400.D0)
     &           /DRAMPExtFlux)
            ENDIF
            IF (NRAMP.EQ.3) THEN
               WRITE(*,*) 'NRAMP = 3 not supported '
               STOP
            ENDIF
         ENDIF
         RAMP = RAMPDG
         
C.......Obtain the meteorological forcing

         IF (NWS.NE.0) CALL MET_FORCING(IT)
         
C.......Compute tidal potential terms

         IF (NTIP.NE.0) CALL TIDAL_POTENTIAL()

C.......Compute LDG auxiliary equations
         
         IF (EVMSUM.NE.0.D0.or.artdif.eq.1) CALL LDG_HYDRO(IT)
         
C.......Compute elevation specified edges

         IF (NEEDS.GT.0)  CALL OCEAN_EDGE_HYDRO(IT)

C.......Compute no-normal flow edges

         IF (NLEDS.GT.0)  CALL LAND_EDGE_HYDRO(IT)

C.......Compute non-zero flow edges

         IF (NFEDS.GT.0)  CALL FLOW_EDGE_HYDRO(IT)
         
C.......Compute radiation edges

         IF (NREDS.GT.0)  CALL RADIATION_EDGE_HYDRO(IT)

C.......Compute internal barrier edges

         IF (NIBEDS.GT.0) CALL IBARRIER_EDGE_HYDRO(IT)
         
C.......Compute external barrier edges

         IF (NEBEDS.GT.0) CALL EBARRIER_EDGE_HYDRO(IT)
         
C.......Compute internal edges

         CALL INTERNAL_EDGE_HYDRO(IT)

C.......Compute elements to finish building the rhs
         
         CALL RHS_DG_HYDRO()

C.......RKC Time Scheme

         if (irk.eq.1) then

            BRK = RKC_tildemu(1)*dt

            DO J = 1,NE
               DO K = 1,DOFS(J)

#ifdef SED_LAY
                  do l = 1,layers
                     
                     bed(K,J,IRK+1,l) = bed(K,J,irk+1,l) + bed(K,J,1,l) 
     &                    + BRK*RHS_bed(K,J,1,l)

                   !adjust for state variable represenation
                     ZE(K,J,IRK+1) = ZE(K,J,IRK+1)-BRK*RHS_bed(K,J,1,l)
                     
                  enddo
#endif
                  
                  ZE(K,J,IRK+1) = ZE(K,J,irk+1) + ZE(K,J,1)
     &                 + BRK*RHS_ZE(K,J,1)
                  QX(K,J,IRK+1) = QX(K,J,irk+1) + QX(K,J,1)
     &                 + BRK*RHS_QX(K,J,1)
                  QY(K,J,IRK+1) = QY(K,J,irk+1)  + QY(K,J,1)
     &                 + BRK*RHS_QY(K,J,1)
                  
C.......Compute the transported tracer term if flagged
                  
#ifdef TRACE
                  iota(K,J,IRK+1) = iota(K,J,irk+1) + iota(K,J,1)
     &                 + BRK*RHS_iota(K,J,1)
#endif
                  
C......Compute chemistry transported terms if flagged
                  
#ifdef CHEM
                  iota(K,J,IRK+1) = iota(K,J,irk+1) + iota(K,J,1)
     &                 + BRK*RHS_iota(K,J,1)
                  iota2(K,J,IRK+1) = iota2(K,J,irk+1) + iota2(K,J,1)
     &                 + BRK*RHS_iota2(K,J,1)
#endif

C.......Compute the dynamic pressure if flagged
                  
#ifdef DYNP
                  dynP(K,J,IRK+1) = dynP(K,J,irk+1) + dynP(K,J,1)
     &                 + BRK*RHS_dynP(K,J,1)
#endif
               ENDDO
            ENDDO
            
         else

            ARK = 1.D0 - RKC_mu(irk)-RKC_nu(irk)
            CRK = RKC_mu(irk)
            DRK = RKC_nu(irk)

            BRK = RKC_tildemu(irk)*dt
            ERK = RKC_gamma(irk)*dt

             DO J = 1,NE
               DO K = 1,DOFS(J)

#ifdef SED_LAY
                  do l = 1,layers
                     
                     bed(K,J,IRK+1,l) = bed(K,J,irk+1,l) 
     &                   + ARK*bed(K,J,1,l) + CRK*bed(K,J,irk,l)
     &                   + DRK*bed(K,J,irk-1,l) + BRK*RHS_bed(K,J,irk,l)
     &                   + ERK*RHS_bed(K,J,1,l)


                   !adjust for state variable represenation
                     ZE(K,J,IRK+1) = ZE(K,J,IRK+1) 
     &                   - BRK*RHS_bed(K,J,irk,l) - ERK*RHS_bed(K,J,1,l)
                     
                  enddo
#endif

                  ZE(K,J,IRK+1) = ZE(K,J,irk+1) + ARK*ZE(K,J,1) 
     &                + CRK*ZE(K,J,irk) + DRK*ZE(K,J,irk-1) 
     &                + BRK*RHS_ZE(K,J,irk) + ERK*RHS_ZE(K,J,1)
                  QX(K,J,IRK+1) = QX(K,J,irk+1) + ARK*QX(K,J,1)
     &                + CRK*QX(K,J,irk) + DRK*QX(K,J,irk-1) 
     &                + BRK*RHS_QX(K,J,irk) + ERK*RHS_QX(K,J,1)
                  QY(K,J,IRK+1) = QY(K,J,irk+1) + ARK*QY(K,J,1) 
     &                + CRK*QY(K,J,irk) + DRK*QY(K,J,irk-1)
     &                + BRK*RHS_QY(K,J,irk) + ERK*RHS_QY(K,J,1)


C.......Compute the transported tracer term if flagged

#ifdef TRACE
                  iota(K,J,IRK+1) = iota(K,J,irk+1) + ARK*iota(K,J,1) 
     &                + CRK*iota(K,J,irk) + DRK*iota(K,J,irk-1) 
     &                + BRK*RHS_iota(K,J,irk) + ERK*RHS_iota(K,J,1)

#endif

C......Compute chemistry transported terms if flagged

#ifdef CHEM
                  iota(K,J,IRK+1) =  iota(K,J,irk+1) + ARK*iota(K,J,1) 
     &                 + CRK*iota(K,J,irk) + DRK*iota(K,J,irk-1) 
     &                 + BRK*RHS_iota(K,J,irk) + ERK*RHS_iota(K,J,1)
                  iota2(K,J,IRK+1) = iota2(K,J,irk+1) + ARK*iota2(K,J,1)
     &                 + CRK*iota2(K,J,irk) + DRK*iota2(K,J,irk-1) 
     &                 + BRK*RHS_iota2(K,J,irk) + ERK*RHS_iota2(K,J,1)

#endif

C.......Compute the dynamic pressure if flagged

#ifdef DYNP
                  dynP(K,J,IRK+1) = dynP(K,J,irk+1) + ARK*dynP(K,J,1) 
     &                 + CRK*dynP(K,J,irk) + DRK*dynP(K,J,irk-1) 
     &                 + BRK*RHS_dynP(K,J,irk) + ERK*RHS_dynP(K,J,1)

#endif

               ENDDO
            ENDDO

         endif

#endif

#ifdef CMPI
         CALL UPDATER_ELEM_MOD(ZE,QX,QY,IRK+1,3)
#ifdef TRACE
         CALL UPDATER_ELEM_MOD(iota,iota2,QY,IRK+1,2 )
#endif

#ifdef DYNP
         CALL UPDATER_ELEM_MOD(dynP,iota2,QY,IRK+1,2 )
#endif

#ifdef CHEM
         CALL UPDATER_ELEM_MOD(iota,iota2,QY,IRK+1,2 )
#endif

#ifdef SED_LAY
         do l = 1,layers
            arrayfix => bed(:,:,:,l) 
            CALL UPDATER_ELEM_MOD(arrayfix,QX,QY,IRK+1,1 )
         enddo
#endif

#endif

C.......Apply the slopelimiter if appropriate

#ifdef SLOPEALL
         CALL SLOPELIMITER()
#endif

#ifdef SLOPE5
         CALL SLOPELIMITER() 
#endif

#ifdef STBLZR
         if (.not.stblzr) then
            CALL SLOPELIMITER()
            stblzr = .true.
         endif
#endif


C.......Apply the wet-dry algorithm if appropriate

      IF (NOLIFA .GE. 2) THEN
         CALL WETDRY()
      ENDIF

c$$$C...Apply the slopelimiter again if auxiliary variable is being used
c$$$
c$$$         if (EVMSUM.NE.0.D0.or.artdif.eq.1) then
c$$$#ifdef SLOPEALL
c$$$            CALL SLOPELIMITER()
c$$$#endif
c$$$            
c$$$#ifdef SLOPE5
c$$$            CALL SLOPELIMITER() 
c$$$#endif
c$$$         endif

C...For parallel run update for dofs if slope limiter and/or wetting and
C...drying is being used

#ifdef CMPI
         if (SLOPEFLAG.NE.0.OR.NOLIFA.GE.2) THEN
            CALL UPDATER_ELEM_MOD(ZE,QX,QY,IRK+1,3)
            
#ifdef TRACE
            CALL UPDATER_ELEM_MOD(iota,iota2,QY,IRK+1,2 )
#endif

#ifdef CHEM
            CALL UPDATER_ELEM_MOD(iota,iota2,QY,IRK+1,2 )
#endif

#ifdef DYNP
            CALL UPDATER_ELEM_MOD(dynP,iota2,QY,IRK+1,2 )
#endif


#ifdef SED_LAY
            do l = 1,layers
               arrayfix => bed(:,:,:,l)
               CALL UPDATER_ELEM_MOD(arrayfix,QX,QY,IRK+1,1 )
            enddo
#endif
         endif
#endif

C.......p_enrich soln and update

#ifdef P_AD
         if (irk+1.eq.2) then
            CALL p_enrichment(it,0)
         endif
#ifdef CMPI

         CALL UPDATER_ELEM_MOD(ZE,QX,QY,IRK+1,3)

#ifdef TRACE
         CALL UPDATER_ELEM_MOD(iota,iota2,QY,IRK+1,2 )
#endif

#ifdef CHEM
         CALL UPDATER_ELEM_MOD(iota,iota2,QY,IRK+1,2 )
#endif

#ifdef DYNP
         CALL UPDATER_ELEM_MOD(dynP,iota2,QY,IRK+1,2 )
#endif

#ifdef SED_LAY
         do l = 1,layers
            arrayfix => bed(:,:,:,l)
            CALL UPDATER_ELEM_MOD(arrayfix,QX,QY,IRK+1,1 )
         enddo
#endif

#endif 

#endif

 100  CONTINUE

C.....Update variables for next time step

#ifdef FILTER 
      !Stabilization filtering 
      DO J = 1,NE
         DO K = 2,DOFS(J)

            alph = 10.D0
            s = RK_stage

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

            seta = eta**s
            sigma = exp(-alph*seta)

C.......Filter the layers if flagged

#ifdef SED_LAY
            do l = 1,layers
               
               bed(K,J,NRK+1,l) = sigma * bed(K,J,nrk+1,l) 
               
            enddo
#endif

            ZE(K,J,NRK+1) = sigma*ZE(K,J,nrk+1) 
            QX(K,J,NRK+1) = sigma*QX(K,J,nrk+1)
            QY(K,J,NRK+1) = sigma*QY(K,J,nrk+1)


C.......Filter the transported tracer term if flagged

#ifdef TRACE
            iota(K,J,NRK+1) = sigma*iota(K,J,nrk+1)

#endif

C......Filter chemistry transported terms if flagged

#ifdef CHEM
            iota(K,J,NRK+1) =  sigma*iota(K,J,nrk+1) 
            iota2(K,J,NRK+1) = sigma*iota2(K,J,nrk+1) 

#endif

         ENDDO
      ENDDO
#endif
      

      DO J = 1,MNE
         
         DO K = 1,DOFS(J)
            ZE(K,J,1) = ZE(K,J,NRK+1)
            QX(K,J,1) = QX(K,J,NRK+1)
            QY(K,J,1) = QY(K,J,NRK+1)
#ifdef TRACE
            iota(K,J,1) = iota(K,J,NRK+1)
#endif

#ifdef CHEM
            iota(K,J,1) = iota(K,J,NRK+1)
            iota2(K,J,1) = iota2(K,J,NRK+1)
#endif

#ifdef DYNP
            dynP(K,J,1) = dynP(K,J,NRK+1)
#endif

#ifdef SED_LAY
            bed(K,J,1,:) = bed(K,J,NRK+1,:)
#endif
         ENDDO
      ENDDO

      DO KK = 2,NRK+1
         DO J = 1,MNE
            DO K = 1,DOFH
               ZE(K,J,KK) = 0.D0
               QX(K,J,KK) = 0.D0
               QY(K,J,KK) = 0.D0
               RHS_ZE(K,J,KK-1) = 0.D0
               RHS_QX(K,J,KK-1) = 0.D0
               RHS_QY(K,J,KK-1) = 0.D0
#ifdef TRACE
               iota(K,J,KK) = 0.D0
               RHS_iota(K,J,KK-1) = 0.D0
#endif

#ifdef CHEM
               iota(K,J,KK) = 0.D0
               iota2(K,J,KK) = 0.D0
               RHS_iota(K,J,KK-1) = 0.D0
               RHS_iota2(K,J,KK-1) = 0.D0
#endif

#ifdef DYNP
               dynP(K,J,KK) = 0.D0
               RHS_dynP(K,J,KK-1) = 0.D0
#endif

#ifdef SED_LAY
               bed(K,J,KK,:) = 0.D0
               hb(K,J,KK) = 0.D0
               RHS_bed(K,J,KK-1,:) = 0.D0
#endif


            ENDDO
         ENDDO
      ENDDO

      RETURN
      END SUBROUTINE
