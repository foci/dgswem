C***********************************************************************
C     
C     SUBROUTINE LDG_HYDRO(IT)
C     
C     Compute variable Z to be used in the LDG terms.
C     
C     Written by Shintaro Bunya (01-01-2007)
C
C     01-10-2011 - cem - adapted for p_enrichment and multicomponent
C     06-01-2012 - cem - sediment diffusion added
C     
C***********************************************************************

      SUBROUTINE LDG_HYDRO(IT)

C.....Use appropriate modules
      
      USE SIZES,ONLY : SZ, myproc
      USE GLOBAL
      USE DG
      USE NodalAttributes, ONLY : EVM

#ifdef CMPI
      USE MESSENGER_ELEM,ONLY : UPDATELZ_ELEM,UPDATEMZ_ELEM
#endif

      IMPLICIT NONE
      
C.....Declare local variables

      INTEGER IT,L,GED,NBOREL,NNBORS,NDRYNBORS,k,i,ll
      INTEGER Detected,ELT
      REAL(SZ) QBCT1,QBCT2,QBC1,QBC2,ZP(3)
      REAL(SZ) DPAVG,lim_by
      Real(SZ),Allocatable :: tmp_mz(:,:,:,:)

      Allocate ( tmp_mz(dofh,2,1,MNE) )

C.....Initialize for viscosity/diffusion
#ifdef WAVE_DIF
      HZ = 0.D0
#endif 
      LZ = 0.D0
#ifdef TRACE
      TZ = 0.d0
#endif
#ifdef SED_LAY
      MZ = 0.d0
#endif
     
C.....Compute elevation specified edges

      IF (NEEDS.GT.0)  CALL OCEAN_EDGE_LDG_HYDRO()

C.....Compute no-normal flow edges

      IF (NLEDS.GT.0)  CALL LAND_EDGE_LDG_HYDRO()

C.....Compute non-zero flow edges

      IF (NFEDS.GT.0)  CALL FLOW_EDGE_LDG_HYDRO()
      
C.....Compute radiation edges

      IF (NREDS.GT.0)  CALL RADIATION_EDGE_LDG_HYDRO()

C.....Compute internal edges

      CALL INTERNAL_EDGE_LDG_HYDRO()

C.....Loop over interior elements

      CALL RHS_LDG_HYDRO()

      do L=1,NE

         N1 = NM(L,1)
         N2 = NM(L,2)
         N3 = NM(L,3)
         if (EDN(L) .eq. 4) then
            N4 = NM(L,4)
         endif 
         if (EVMSUM.NE.0.D0) then
            if (EDN(L) .eq. 3) then
               EVMAvg = (EVM(N1)+EVM(N2)+EVM(N3))/3.d0
            elseif (EDN(L) .eq. 4) then
               EVMAvg = (EVM(N1)+EVM(N2)+EVM(N3)+EVM(N4))/4.d0
            endif
         else
            EVMAvg = 0.D0
         endif

         bg_dif = 0.D0 
         trc_dif = 0.D0 
          
#ifdef SED_LAY
         SEVDMAvg = SEVDM 
#endif
 
#ifdef ARTDIF

         lim_by = 1.D0
#endif
         
#ifdef WAVE_DIF
         if (entrop(1,L).gt.s0+kappa) then
            
            HZ(:,:,:,L) = (bg_dif+e1(1)) *
     &           (abs(lim_by+balance(1)))**pa * HZ(:,:,:,L)

         else

            HZ(:,:,:,L) = bg_dif*HZ(:,:,:,L)

         endif
#endif

         if (entrop(2,L).gt.s0+kappa) then
            
            LZ(:,1,1,L) = (EVMAvg+e1(2) * 
     &           (abs(lim_by+balance(2)))**pa)  * LZ(:,1,1,L)
            LZ(:,1,2,L) = (EVMAvg+e1(2) *
     &           (abs(lim_by+balance(2)))**pa ) * LZ(:,1,2,L)
            LZ(:,2,1,L) = (EVMAvg+e1(2) *
     &           (abs(lim_by+balance(2)))**pa ) * LZ(:,2,1,L)
     
         else

            LZ(:,1,:,L) = EVMAvg* LZ(:,1,:,L)
            LZ(:,2,1,L) = EVMAvg* LZ(:,2,1,L)
            
         endif

         if (entrop(3,L).gt.s0+kappa) then
            
            LZ(:,2,2,L) = (EVMAvg+e1(3) * 
     &           (abs(lim_by+balance(3)))**pa) * LZ(:,2,2,L)

            if (entrop(2,L).le.s0+kappa) then
               
               LZ(:,2,1,L) = (EVMAvg+e1(3) *
     &              (abs(lim_by+balance(3)))**pa ) * LZ(:,2,1,L)
               LZ(:,1,2,L) = (EVMAvg+e1(3) *
     &              (abs(lim_by+balance(3)))**pa ) * LZ(:,1,2,L)
               
            endif
            
         else
            
            LZ(:,2,2,L) = EVMAvg* LZ(:,2,2,L)
            
            if (entrop(2,L).le.s0+kappa) then
               
               LZ(:,1,2,L) = EVMAvg* LZ(:,1,2,L)
               LZ(:,2,1,L) = EVMAvg* LZ(:,2,1,L)
               
            endif
            
         endif

#ifdef TRACE
         if (entrop(4,L).gt.s0+kappa) then
            
            TZ(:,:,:,L) = (trc_dif+e1(4)) *
     &           (abs(lim_by+balance(4)))**pa * TZ(:,:,:,L)

         else

            TZ(:,:,:,L) = trc_dif*TZ(:,:,:,L)

         endif
#endif

#ifdef SED_LAY

         if (entrop(5,L).gt.s0+kappa) then
            
            MZ(:,:,:,L) = (SEVDMAvg+e1(5)) *
     &           (abs(lim_by+balance(5)))**pa * MZ(:,:,:,L)

         else

            MZ(:,:,:,L) = MZ(:,:,:,L)*SEVDMAvg
            
         endif
#endif
         
#ifdef WAVE_DIF
         HZ(:,:,:,L) = bg_dif*HZ(:,:,:,L)
#endif
         LZ(:,:,:,L) = EVMAvg*LZ(:,:,:,L)
               
#ifdef TRACE
         TZ(:,:,:,L) = trc_dif*TZ(:,:,:,L)
#endif

#ifdef SED_LAY
         MZ(:,:,:,L) = MZ(:,:,:,L)*SEVDMAvg
#endif

      enddo
C********************************************

#ifdef CMPI
#ifdef WAVE_DIF
      CALL UPDATELZ_ELEM(HZ)
#endif
      CALL UPDATELZ_ELEM(LZ)
#ifdef TRACE
      CALL UPDATELZ_ELEM(TZ)
#endif
#ifdef SED_LAY
      do ll=1,layers
         tmp_mz(:,:,1,:) = MZ(:,:,ll,:) 
         CALL UPDATEMZ_ELEM(tmp_mz)
         MZ(:,:,ll,:) = tmp_mz(:,:,1,:)
      enddo
#endif
#endif

      RETURN
      END SUBROUTINE

C***********************************************************************
C     
C     SUBROUTINE OCEAN_EDGE_LDG_HYDRO( )
C     
C     This subroutine does the following:
C     
C     1.  Calculates the values of the necessary variables at the edge
C     gauss points for NO-NORMAL FLOW edges
C     2.  Compute the LDG flux at these points (sediment diffusion?).
C     3.  Compute the boundary integrals.
C     
C     Written by Shintaro Bunya (01-04-2007)
C     01-10-2011 - cem - adapted for p_enrichment and multicomponent
C     
C***********************************************************************
      SUBROUTINE OCEAN_EDGE_LDG_HYDRO()

C.....Use appropriate modules

      USE SIZES,ONLY : layers
      USE GLOBAL
      USE DG

      IMPLICIT NONE

C.....Declare local variables

      INTEGER L,LED,GED,k,I,ll,j
      REAL(SZ) QX_AVG, QY_AVG,ZE_AVG,bed_AVG(layers)
      REAL(SZ) iota_AVG

      DO 1000 L = 1, needs
         
C.....Retrieve the global and local edge number

         GED = NEEDN(L)
         LED = NEDSD(1,GED)

C.....Retrieve the elements which share the edge

         EL_IN = NEDEL(1,GED)

         J = EDN(EL_IN) - 2
         pa = PDG_EL(EL_IN)

#ifdef P0
         if (pa.eq.0) then
            pa = 1
         endif
#endif
         
C.....Retrieve the components of the normal vector to the edge
         
         NX = COSNX(GED)
         NY = SINNX(GED)
C.....Retrieve the nodes of the edge
         
         N1 = NEDNO(1,GED)
         N2 = NEDNO(2,GED)
         
C.....Compute ZE, QX, QY, and HB at each edge Gauss quadrature point

         DO I = 1,NEGP(pa)

#ifdef WAVE_DIF
            ZE_IN = 0.D0
            ZE_EX = 0.D0
#endif
            QX_IN = 0.D0
            QY_IN = 0.D0

            QX_EX = 0.D0
            QY_EX = 0.D0

#ifdef TRACE
            iota_IN = 0.D0
            iota_EX = 0.D0
#endif

            !deal with sediment
#ifdef SED_LAY
             do ll=1,layers
               bed_IN(ll) = 0.D0
               bed_EX(ll) = 0.D0
            enddo
#endif
           
C.....Compute the solution at the interior state

            DO K = 1,DOFS(EL_IN)
#ifdef WAVE_DIF
               ZE_IN = ZE_IN + ZE(K,EL_IN,IRK)*PHI_EDGE(K,I,LED,pa,J)
#endif
               QX_IN = QX_IN + QX(K,EL_IN,IRK)*PHI_EDGE(K,I,LED,pa,J)
               QY_IN = QY_IN + QY(K,EL_IN,IRK)*PHI_EDGE(K,I,LED,pa,J)
#ifdef TRACE
            iota_IN = iota_IN + iota(K,EL_IN,IRK)*PHI_EDGE(K,I,LED,pa,J)
#endif
#ifdef SED_LAY
               do ll=1,layers
                  bed_IN(ll) = bed_IN(ll) 
     &                 + bed(K,EL_IN,IRK,ll)*PHI_EDGE(K,I,LED,pa,J)
               enddo
#endif
            ENDDO
            
C.....Set the exterior state flows equal to the interior state flows

            QX_EX = QX_IN
            QY_EX = QY_IN

#ifdef SED_LAY
            do ll=1,layers
               bed_EX(ll) = bed_IN(ll)
            enddo
#endif

C.....Take the average
#ifdef WAVE_DIF
            ZE_AVG = ZE_IN
#endif
            QX_AVG = 0.5D0*(QX_IN + QX_EX)*SFACED(I,LED,EL_IN,pa)
            QY_AVG = 0.5D0*(QY_IN + QY_EX)
#ifdef TRACE
            iota_AVG = iota_IN
#endif
#ifdef SED_LAY
            do ll=1,layers
               bed_AVG(ll) = 0.5D0*( bed_IN(ll) + bed_EX(ll) )
            enddo
#endif

C.....Compute the edge integral
            
            DO K = 1,DOFS(EL_IN)
               CALL EDGE_INT_LDG_HYDRO
     &        (K,EL_IN,LED,GED,I,iota_AVG,ZE_AVG,QX_AVG,QY_AVG,NX,NY,pa)
#ifdef SED_LAY
               do ll=1,layers
                  CALL EDGE_INT_LDG_sediment
     &                       (K,EL_IN,LED,GED,I,bed_AVG(ll),NX,NY,pa,ll)
               enddo
#endif
            ENDDO
            
         ENDDO

 1000 CONTINUE
      
      RETURN
      END SUBROUTINE


C***********************************************************************
C     
C     SUBROUTINE LAND_EDGE_HYDRO_LDG( )
C     
C     This subroutine does the following:
C     
C     1.  Calculates the values of the necessary variables at the edge
C     gauss points for NO-NORMAL FLOW edges
C     2.  Compute the LDG flux at these points.
C     3.  Compute the boundary integrals.
C     
C     Written by Shintaro Bunya (01-02-2007)
C     01-10-2011 - cem - adapted for p_enrichment and multicomponent
C     
C***********************************************************************

      SUBROUTINE LAND_EDGE_LDG_HYDRO()

C.....Use appropriate modules

      USE SIZES,ONLY : SZ,layers
      USE GLOBAL
      USE DG
      IMPLICIT NONE

C.....Declare local variables


      INTEGER K,L,LED,GED,GP,i,kk,ll,J
      REAL(SZ) AREA, IMASS
      REAL(SZ) TX, TY, QX_AVG, QY_AVG, bed_AVG(layers)
      REAL(SZ) ZE_AVG,iota_AVG
      
      DO 1000 L = 1,NLEDS
         
C.....Retrieve the global and local edge number
         
         GED = NLEDN(L)
         LED = NEDSD(1,GED)
         
C.....Retrieve the elements which share the edge
         
         EL_IN = NEDEL(1,GED)
         
         pa = PDG_EL(EL_IN)
         
#ifdef P0
         if (pa.eq.0) then
            pa = 1
         endif
#endif

C.....Retrieve the components of the normal vector to the edge
         
         NX = COSNX(GED)
         NY = SINNX(GED)
         
C.....Set the components for the tangential vector to the edge
         
         TX = -NY
         TY =  NX

C.....Compute ZE, QX, QY, and HB at each Gauss point
         J = EDN(EL_IN) - 2
         DO I = 1,NEGP(pa)
#ifdef WAVE_DIF
            ZE_IN = 0.D0
#endif
            QX_IN = 0.D0
            QY_IN = 0.D0
#ifdef TRACE
            iota_IN = 0.d0
#endif
#ifdef SED_LAY
            do ll=1,layers
               bed_IN(ll) = 0.D0
            enddo
#endif

C.....Compute the solution at the interior state

            DO K = 1,DOFS(EL_IN)
#ifdef WAVE_DIF
               ZE_IN = ZE_IN + ZE(K,EL_IN,IRK)*PHI_EDGE(K,I,LED,pa,J)
#endif
               QX_IN = QX_IN + QX(K,EL_IN,IRK)*PHI_EDGE(K,I,LED,pa,J)
               QY_IN = QY_IN + QY(K,EL_IN,IRK)*PHI_EDGE(K,I,LED,pa,J)
#ifdef TRACE
            iota_IN = iota_IN + iota(K,EL_IN,IRK)*PHI_EDGE(K,I,LED,pa,J)
#endif
#ifdef SED_LAY
               do ll=1,layers
                  bed_IN(ll) = bed_IN(ll) + bed(K,EL_IN,IRK,ll)*
     &                                            PHI_EDGE(K,I,LED,pa,J)
               enddo
#endif

            ENDDO
            
C.....Compute the velocity in the normal and tangental direction
            
            Q_N_INT = QX_IN*NX + QY_IN*NY
            Q_T_INT = QX_IN*TX + QY_IN*TY

C.....Reflect the velocity in the normal direction

            Q_N_EXT = -Q_N_INT
            Q_T_EXT =  Q_T_INT
            
C.....Compute the x and y components of the external state flow

            QX_EX = ( TY*Q_N_EXT - NY*Q_T_EXT)/(NX*TY - NY*TX)
            QY_EX = (-TX*Q_N_EXT + NX*Q_T_EXT)/(NX*TY - NY*TX)

#ifdef SED_LAY
            do ll=1,layers !maybe this should be changed?
               bed_EX(ll) =  bed_IN(ll)
            enddo
#endif            


C.....Take the average
#ifdef WAVE_DIF
            ZE_AVG = ZE_IN
#endif
            QX_AVG = 0.5D0*(QX_IN + QX_EX)*SFACED(I,LED,EL_IN,pa)
            QY_AVG = 0.5D0*(QY_IN + QY_EX)
#ifdef TRACE
            iota_AVG = iota_IN
#endif
#ifdef SED_LAY
            do ll=1,layers !take the Fronenius norm
               bed_AVG(ll) = 0.5D0*( bed_IN(ll) + bed_EX(ll) )
            enddo
#endif

C.....Compute the edge integral

            DO K = 1,DOFS(EL_IN)
               CALL EDGE_INT_LDG_HYDRO
     &        (K,EL_IN,LED,GED,I,iota_AVG,ZE_AVG,QX_AVG,QY_AVG,NX,NY,pa)
#ifdef SED_LAY
               do ll=1,layers
                  CALL EDGE_INT_LDG_sediment
     &                 (K,EL_IN,LED,GED,I,bed_AVG(ll),NX,NY,pa,ll)
               enddo
#endif
            ENDDO

         ENDDO
            
 1000 CONTINUE
      
      RETURN
      END SUBROUTINE


C***********************************************************************
C     
C     SUBROUTINE FLOW_EDGE_LDG_HYDRO( )
C     
C     This subroutine does the following:
C     
C     1.  Calculates the values of the necessary variables at the edge
C     gauss points for INTERNAL edges
C     2.  Take the average of these values (sediment diffusion)
C     3.  Perform boundary integration
C     
C     Written by Shintaro Bunya (01-05-2007)
C     01-10-2011 - cem - adapted for p_enrichment and multicomponent
C     
C***********************************************************************

      SUBROUTINE FLOW_EDGE_LDG_HYDRO()

C.....Use appropriate modules

      USE SIZES,ONLY : layers
      USE GLOBAL
      USE DG

      IMPLICIT NONE

C.....Declare local variables

      INTEGER L, LED, GED,k,i,jj,ll,J
      REAL(SZ) TX, TY, QX_AVG, QY_AVG, bed_AVG(layers)
      REAL(SZ) iota_AVG,ZE_AVG

      DO 1000 L = 1,NFEDS
         
C.....Retrieve the global and local edge number

         GED = NFEDN(L)
         LED = NEDSD(1,GED)

C.....Retrieve the element to which the edge belongs

         EL_IN = NEDEL(1,GED)

         PA = PDG_EL(EL_IN)
         J = EDN(EL_IN) - 2
#ifdef P0
         if (pa.eq.0) then
            pa = 1
         endif
#endif
         
C.....Retrieve the components of the normal vector to the edge
         
         NX = COSNX(GED)
         NY = SINNX(GED)
         
C.....Set the components for the tangential vector to the edge

         TX = -NY
         TY =  NX
         
C.....Compute ZE, QX, QY, and HB at each edge Gauss quadrature point

         DO I = 1,NEGP(pa)

#ifdef WAVE_DIF
            ZE_IN = 0.D0
            ZE_EX = 0.D0
#endif

            QX_IN = 0.D0
            QY_IN = 0.D0

            QX_EX = 0.D0
            QY_EX = 0.D0
#ifdef TRACE
            iota_IN = 0.D0
            iota_EX = 0.D0
#endif
#ifdef SED_LAY
             do ll=1,layers
               bed_IN(ll) = 0.D0
               bed_EX(ll) = 0.D0
            enddo
#endif
            
C.....Compute the specified flow boundaries for the exterior state

            Q_N_EXT = 0.D0
            DO JJ = 1,MNFFR
               IF (NFFR.EQ.0) THEN
                  NCYC = 0.D0
                  ARGJ = 0.D0
                  RFF = RAMPExtFlux
               ELSEIF (FPER(JJ).EQ.0.D0) THEN
                  NCYC = 0.D0
                  ARGJ = FAMIG(JJ)*(TIMEDG - NCYC*FPER(JJ)) + FFACE(JJ)
                  RFF  = FFF(JJ)*RAMPExtFlux
               ELSE
                  NCYC = INT(TIMEDG/FPER(JJ))
                  ARGJ = FAMIG(JJ)*(TIMEDG - NCYC*FPER(JJ)) + FFACE(JJ)
                  RFF  = FFF(JJ)*RAMPExtFlux
               ENDIF
               
               QNAM_GP = 0.5D0*(QNAM_DG(JJ,L,1) + QNAM_DG(JJ,L,2))
     &            + 0.5D0*(QNAM_DG(JJ,L,2) - QNAM_DG(JJ,L,1))*XEGP(I,pa)
               QNPH_GP = 0.5D0*(QNPH_DG(JJ,L,1) + QNPH_DG(JJ,L,2))
     &            + 0.5D0*(QNPH_DG(JJ,L,2) - QNPH_DG(JJ,L,1))*XEGP(I,pa)
               
               ARG = ARGJ - QNPH_GP
               
               Q_N_EXT = Q_N_EXT + QNAM_GP*RFF*COS(ARG)
               Q_T_EXT =  0.D0


            ENDDO
            
C.....Compute the solution at the interior state

            DO K = 1,DOFS(EL_IN)

#ifdef WAVE_DIF
               ZE_IN = ZE_IN + ZE(K,EL_IN,IRK)*PHI_EDGE(K,I,LED,pa,J)
#endif

               QX_IN = QX_IN + QX(K,EL_IN,IRK)*PHI_EDGE(K,I,LED,pa,J)
               QY_IN = QY_IN + QY(K,EL_IN,IRK)*PHI_EDGE(K,I,LED,pa,J)

#ifdef SED_LAY
               do ll=1,layers
                  bed_in(ll) = bed_in(ll) +
     &                        bed(K,EL_IN,IRK,ll)*PHI_EDGE(K,I,LED,pa,J)
                  bed_ex(ll) = bed_ex(ll) + 
     &                        bed(K,EL_IN,IRK,ll)*PHI_EDGE(K,I,LED,pa,J)
               enddo
#endif
            ENDDO


               QX_EX = -( TY*Q_N_EXT - NY*Q_T_EXT)/(NX*TY - NY*TX)
               QY_EX = -(-TX*Q_N_EXT + NX*Q_T_EXT)/(NX*TY - NY*TX)
            
C.....Take the average

#ifdef WAVE_DIF
               ZE_AVG = ZE_IN
#endif
               QX_AVG = 0.5*(QX_IN + QX_EX)*SFACED(I,LED,EL_IN,pa)
               QY_AVG = 0.5*(QY_IN + QY_EX)
#ifdef TRACE
               iota_AVG = iota_IN
#endif

#ifdef SED_LAY
               do ll=1,layers
                  bed_AVG(ll) = 0.5D0*( bed_IN(ll) + bed_EX(ll) )
               enddo
#endif

C.....Compute the edge integral

            DO K = 1,DOFS(EL_IN)
               CALL EDGE_INT_LDG_HYDRO
     &        (K,EL_IN,LED,GED,I,iota_AVG,ZE_AVG,QX_AVG,QY_AVG,NX,NY,pa)
#ifdef SED_LAY
               do ll=1,layers
                  CALL EDGE_INT_LDG_sediment
     &                 (K,EL_IN,LED,GED,I,bed_AVG(ll),NX,NY,pa,ll)
               enddo
#endif
            ENDDO
         ENDDO

 1000 CONTINUE
      
      RETURN
      END SUBROUTINE


C***********************************************************************
C     
C     SUBROUTINE RADIATION_EDGE_LDG_HYDRO( )
C     
C     This subroutine does the following:
C     
C     1.  Calculates the values of the necessary variables at the edge
C     gauss points for INTERNAL edges (sediment diffusion?)
C     2.  Take the average of these values
C     3.  Perform boundary integration
C     
C     Written by Shintaro Bunya (01-10-2007)
C     01-10-2011 - cem - adapted for p_enrichment and multicomponent
C     
C***********************************************************************

      SUBROUTINE RADIATION_EDGE_LDG_HYDRO()

C.....Use appropriate modules

      USE SIZES,ONLY : layers
      USE GLOBAL
      USE DG

      IMPLICIT NONE

C.....Declare local variables

      INTEGER L, LED, GED,k,i,ll,J
      REAL(SZ) TX, TY

      DO 1000 L = 1,NREDS
         
C.....Retrieve the global and local edge number

         GED = NREDN(L)
         LED = NEDSD(1,GED)

C.....Retrieve the elements which share the edge

         EL_IN = NEDEL(1,GED)
         J = EDN(EL_IN) - 2
         PA = PDG_EL(EL_IN)

#ifdef P0
         if (pa.eq.0) then
            pa = 1
         endif
#endif

C.....Retrieve the components of the normal vector to the edge
         
         NX = COSNX(GED)
         NY = SINNX(GED)

C.....Compute ZE, QX, QY, and HB at each Gauss point

         DO I = 1,NEGP(pa)
            
#ifdef WAVE_DIF
            ZE_IN = 0.D0
#endif
            QX_IN = 0.D0
            QY_IN = 0.D0
#ifdef TRACE
            iota_IN = 0.d0
#endif
#ifdef SED_LAY
            do ll=1,layers
               bed_IN(ll) = 0.D0
            enddo
#endif

C.....Compute the solution at the interior state

            DO K = 1,DOFS(EL_IN)
               
#ifdef WAVE_DIF
               ZE_IN = ZE_IN + ZE(K,EL_IN,IRK)*PHI_EDGE(K,I,LED,pa,J)
#endif
               QX_IN = QX_IN + QX(K,EL_IN,IRK)*PHI_EDGE(K,I,LED,pa,J)
               QY_IN = QY_IN + QY(K,EL_IN,IRK)*PHI_EDGE(K,I,LED,pa,J)
#ifdef TRACE
            iota_IN = iota_IN + iota(K,EL_IN,IRK)*PHI_EDGE(K,I,LED,pa,J)
#endif

#ifdef SED_LAY
               do ll=1,layers
                  bed_in(ll) = bed_in(ll) + 
     &                        bed(K,EL_IN,IRK,ll)*PHI_EDGE(K,I,LED,pa,J)
               enddo
#endif
        

            ENDDO

            QX_IN = QX_IN*SFACED(I,LED,EL_IN,pa)

C.....Compute the edge integral
           
            DO K = 1,DOFS(EL_IN)
               CALL EDGE_INT_LDG_HYDRO
     &            (K,EL_IN,LED,GED,I,iota_IN,ZE_IN,QX_IN,QY_IN,NX,NY,pa)
#ifdef SED_LAY
               do ll=1,layers
                  CALL EDGE_INT_LDG_sediment
     &                 (K,EL_IN,LED,GED,I,bed_IN(ll),NX,NY,pa,ll)
               enddo
#endif
            ENDDO

         ENDDO

 1000 CONTINUE      
      RETURN
      END SUBROUTINE

C***********************************************************************
C     
C     SUBROUTINE INTERNAL_EDGE_LDG_HYDRO( )
C     
C     This subroutine does the following:
C     
C     1.  Calculates the values of the necessary variables at the edge
C     gauss points for INTERNAL edges (sediment diffusion?)
C     2.  Take the average of these values
C     3.  Perform boundary integration
C     
C     Written by Shintaro Bunya (01-02-2007)
C     01-10-2011 - cem - adapted for p_enrichment and multicomponent
C     
C***********************************************************************

      SUBROUTINE INTERNAL_EDGE_LDG_HYDRO()

C.....Use appropriate modules
      
      USE SIZES,ONLY : SZ,layers
      USE GLOBAL
      USE DG

      IMPLICIT NONE

C.....Declare local variables

      INTEGER L, LED_IN, LED_EX, GED, GP_IN, GP_EX,k,i,ll,JI,JX,kk
      REAL(SZ) ZE_AVG,QX_AVG,QY_AVG,bed_AVG(layers),W_IN,W_EX
      real(sz) iota_AVG
      real(sz), allocatable :: vtemp1(:),vtemp2(:)


      DO 1000 L = 1,NIEDS

C.....Retrieve the global and local edge number

         GED = NIEDN(L)
         LED_IN = NEDSD(1,GED)
         LED_EX = NEDSD(2,GED)

C.....Retrieve the elements which share the edge

         EL_IN = NEDEL(1,GED)
         EL_EX = NEDEL(2,GED)

         EL = EL_IN
         IF (DOFS(EL_EX).LT.DOFS(EL_IN)) then
            EL = EL_EX
         endif
C.....Determine whether each element is a triangle or quad; value is 1 for tri, 2 for quad
         JI = EDN(EL_IN) - 2
         JX = EDN(EL_EX) - 2
         
         pa = PDG_EL(EL)
#ifdef P0
         if (pa.eq.0) then
            pa = 1
         endif
#endif   
         if (JI .eq. 2 .or. JX .eq. 2) then
            allocate( vtemp1(DOFS(EL)),vtemp2(DOFS(EL)) )
         endif
C.....Retrieve the components of the normal vector to the edge

         NX = COSNX(GED)
         NY = SINNX(GED)

C.....Compute ZE, QX, QY, and HB at each edge Gauss quadrature point

         DO I = 1,NEGP(pa)

            GP_IN = I
            GP_EX = NEGP(pa) - I + 1
            
#ifdef WAVE_DIF
            ZE_IN = 0.D0
#endif
            QX_IN = 0.D0
            QY_IN = 0.D0
#ifdef WAVE_DIF
            ZE_EX = 0.D0
#endif
            QX_EX = 0.D0
            QY_EX = 0.D0
#ifdef TRACE
            iota_IN = 0.D0
            iota_EX = 0.D0
#endif

!deal with sediment
#ifdef SED_LAY
            do ll=1,layers
               bed_IN(ll) = 0.D0
               bed_EX(ll) = 0.D0
            enddo
#endif

            DO K = 1,DOFS(EL)

#ifdef WAVE_DIF
                  ZE_IN = ZE_IN + 
     &                    ZE(K,EL_IN,IRK)*PHI_EDGE(K,GP_IN,LED_IN,pa,JI)
#endif
                  QX_IN = QX_IN + 
     &                    QX(K,EL_IN,IRK)*PHI_EDGE(K,GP_IN,LED_IN,pa,JI)
                  QY_IN = QY_IN + 
     &                    QY(K,EL_IN,IRK)*PHI_EDGE(K,GP_IN,LED_IN,pa,JI)
#ifdef WAVE_DIF
                  ZE_EX = ZE_EX + 
     &                    ZE(K,EL_EX,IRK)*PHI_EDGE(K,GP_EX,LED_EX,pa,JX)
#endif
                  QX_EX = QX_EX + 
     &                    QX(K,EL_EX,IRK)*PHI_EDGE(K,GP_EX,LED_EX,pa,JX)
                  QY_EX = QY_EX + 
     &                    QY(K,EL_EX,IRK)*PHI_EDGE(K,GP_EX,LED_EX,pa,JX)
#ifdef TRACE
      iota_IN = iota_IN+iota(K,EL_IN,IRK)*PHI_EDGE(K,GP_IN,LED_IN,pa,JI)
      iota_EX = iota_EX+iota(K,EL_EX,IRK)*PHI_EDGE(K,GP_EX,LED_EX,pa,JX)
#endif

#ifdef SED_LAY
               do ll=1,layers
                  bed_IN(ll) = bed_IN(ll) +
     &                bed(K,EL_IN,IRK,ll)*PHI_EDGE(K,GP_IN,LED_IN,pa,JI)
                  bed_EX(ll) = bed_EX(ll) +
     &                bed(K,EL_EX,IRK,ll)*PHI_EDGE(K,GP_EX,LED_EX,pa,JX)
               enddo
#endif

            ENDDO

C.....Take the average

#ifdef WAVE_DIF
            ZE_AVG = 0.5D0*(ZE_IN + ZE_EX)
#endif

            QX_AVG = 0.5D0*(QX_IN*SFACED(GP_IN,LED_IN,EL_IN,pa) 
     &                            + QX_EX*SFACED(GP_EX,LED_EX,EL_EX,pa))
            QY_AVG = 0.5D0*(QY_IN + QY_EX)

#ifdef TRACE
            iota_AVG = 0.5D0*(iota_IN + iota_EX)
#endif

#ifdef SED_LAY
            do ll=1,layers
               bed_AVG(ll) = 0.5*( bed_IN(ll) + bed_EX(ll) )
            enddo
#endif

C.....Compute the edge integral

            DO K = 1,DOFS(EL)

               if (EDN(EL_IN) .eq. 3) then
                  W_IN = 2.0d0*M_INV(K,pa)/AREAS(EL_IN)*XLEN(GED)*
     &                     PHI_EDGE(K,GP_IN,LED_IN,pa,JI)*WEGP(GP_IN,pa)
               elseif (EDN(EL_IN) .eq. 4) then
                  do kk = 1,DOFS(EL)
                     vtemp1(kk) = M_INVM(EL_IN,K,kk,pa)
                     vtemp2(kk) = PHI_EDGE(kk,GP_IN,LED_IN,pa,JI)
                  enddo
                  W_IN = dot_product(vtemp1,vtemp2)
                  W_IN = W_IN*0.5d0*XLEN(GED)*WEGP(GP_IN,pa)
               endif
               
#ifdef WAVE_DIF
               HZ(K,1,1,EL_IN) = HZ(K,1,1,EL_IN) - ZE_AVG*NX*W_IN
               HZ(K,2,2,EL_IN) = HZ(K,2,2,EL_IN) - ZE_AVG*NY*W_IN
#endif
               LZ(K,1,1,EL_IN) = LZ(K,1,1,EL_IN) - QX_AVG*NX*W_IN
               LZ(K,1,2,EL_IN) = LZ(K,1,2,EL_IN) - QX_AVG*NY*W_IN
               LZ(K,2,1,EL_IN) = LZ(K,2,1,EL_IN) - QY_AVG*NX*W_IN
               LZ(K,2,2,EL_IN) = LZ(K,2,2,EL_IN) - QY_AVG*NY*W_IN
            
               if (EDN(EL_EX) .eq. 3) then
                  W_EX = 2.0d0*M_INV(K,pa)/AREAS(EL_EX)*XLEN(GED)*
     &                     PHI_EDGE(K,GP_EX,LED_EX,pa,JX)*WEGP(GP_EX,pa)
               elseif (EDN(EL_EX) .eq. 4) then
                  do kk = 1,DOFS(EL)
                     vtemp1(kk) = M_INVM(EL_EX,K,kk,pa)
                     vtemp2(kk) = PHI_EDGE(kk,GP_EX,LED_EX,pa,JX)
                  enddo
                  W_EX = dot_product(vtemp1,vtemp2)
                  W_EX = W_EX*0.5d0*XLEN(GED)*WEGP(GP_EX,pa)
               endif

#ifdef WAVE_DIF
               HZ(K,1,1,EL_EX) = HZ(K,1,1,EL_EX) + ZE_AVG*NX*W_EX
               HZ(K,2,2,EL_EX) = HZ(K,2,2,EL_EX) + ZE_AVG*NY*W_EX
#endif
               LZ(K,1,1,EL_EX) = LZ(K,1,1,EL_EX) + QX_AVG*NX*W_EX
               LZ(K,1,2,EL_EX) = LZ(K,1,2,EL_EX) + QX_AVG*NY*W_EX
               LZ(K,2,1,EL_EX) = LZ(K,2,1,EL_EX) + QY_AVG*NX*W_EX
               LZ(K,2,2,EL_EX) = LZ(K,2,2,EL_EX) + QY_AVG*NY*W_EX

#ifdef TRACE
               TZ(K,1,1,EL_IN) = TZ(K,1,1,EL_IN) - iota_AVG*NX*W_IN
               TZ(K,2,2,EL_IN) = TZ(K,2,2,EL_IN) - iota_AVG*NY*W_IN
               TZ(K,1,1,EL_EX) = TZ(K,1,1,EL_EX) + iota_AVG*NX*W_EX
               TZ(K,2,2,EL_EX) = TZ(K,2,2,EL_EX) + iota_AVG*NY*W_EX
#endif

#ifdef SED_LAY
               do ll=1,layers
                  MZ(K,1,ll,EL_IN) = MZ(K,1,ll,EL_IN) 
     &                 - bed_AVG(ll)*NX*W_IN
                  MZ(K,2,ll,EL_IN) = MZ(K,2,ll,EL_IN) 
     &                 - bed_AVG(ll)*NY*W_IN
                  MZ(K,1,ll,EL_EX) = MZ(K,1,ll,EL_EX) 
     &                 + bed_AVG(ll)*NX*W_EX
                  MZ(K,2,ll,EL_EX) = MZ(K,2,ll,EL_EX) 
     &                 + bed_AVG(ll)*NY*W_EX
               enddo
#endif

            enddo
         enddo
         if (JI .eq. 2 .or. JX .eq. 2) deallocate( vtemp1,vtemp2 )
 1000 CONTINUE
      RETURN
      END SUBROUTINE


C***********************************************************************
C     
C     SUBROUTINE EDGE_INT_LDG_HYDRO()
C     
C     This subroutine computes the edge integrals for the LDG boundary
C     terms using Gauss quadrature and adds them to the LZ. (sediment diffusion?)
C     
C     Written by Shintaro Bunya (01-02-2007)
C     01-10-2011 - cem - adapted for p_enrichment and multicomponent
C     
C***********************************************************************

      SUBROUTINE EDGE_INT_LDG_HYDRO(K,EL,LED,GED,GP,iota_Avg,ZE_Avg,
     &     QX_Avg, QY_Avg,NX,NY,pa)

C.....Use appropriate modules

      USE SIZES,ONLY : SZ,layers
      USE GLOBAL,ONLY : AREAS,nm,EDN
      USE DG,ONLY : M_INV,XLEN,PHI_EDGE,WEGP,LZ,HZ,TZ,M_INVM
 
      IMPLICIT NONE
      
C.....Declare local variables

      INTEGER K,EL,LED,GED,GP,i,pa,pp,J
      REAL(SZ) AREA,IMASS,n1,n2,n3,QX_Avg,QY_Avg
      Real(SZ) ZE_Avg,iota_Avg
      real(sz) NX,NY,vtemp1((pa+1)*(pa+2)/2),vtemp2((pa+1)*(pa+2)/2)
      
C.....Retrieve the element area
      
      AREA = 0.5D0*AREAS(EL)

C.....Comput the edge integral
      J = EDN(EL) - 2
      if (EDN(EL) .eq. 3) then
         IMASS = M_INV(K,pa)/(0.5D0*AREA)
#ifdef WAVE_DIF
         HZ(K,1,1,EL) = HZ(K,1,1,EL) - IMASS*XLEN(GED)*
     &               0.5D0*ZE_Avg*PHI_EDGE(K,GP,LED,pa,J)*NX*WEGP(GP,pa)
         HZ(K,2,2,EL) = HZ(K,2,2,EL) - IMASS*XLEN(GED)*
     &               0.5D0*ZE_Avg*PHI_EDGE(K,GP,LED,pa,J)*NY*WEGP(GP,pa)
#endif
         LZ(K,1,1,EL) = LZ(K,1,1,EL) - IMASS*XLEN(GED)*
     &               0.5D0*QX_Avg*PHI_EDGE(K,GP,LED,pa,J)*NX*WEGP(GP,pa)
         LZ(K,1,2,EL) = LZ(K,1,2,EL) - IMASS*XLEN(GED)*
     &               0.5D0*QX_Avg*PHI_EDGE(K,GP,LED,pa,J)*NY*WEGP(GP,pa)
         LZ(K,2,1,EL) = LZ(K,2,1,EL) - IMASS*XLEN(GED)*
     &               0.5D0*QY_Avg*PHI_EDGE(K,GP,LED,pa,J)*NX*WEGP(GP,pa)
         LZ(K,2,2,EL) = LZ(K,2,2,EL) - IMASS*XLEN(GED)*
     &               0.5D0*QY_Avg*PHI_EDGE(K,GP,LED,pa,J)*NY*WEGP(GP,pa)

#ifdef TRACE
         TZ(K,1,1,EL) = TZ(K,1,1,EL)  - IMASS*XLEN(GED)*
     &             0.5D0*iota_Avg*PHI_EDGE(K,GP,LED,pa,J)*NX*WEGP(GP,pa)
         TZ(K,2,2,EL) = TZ(K,2,2,EL) - IMASS*XLEN(GED)*
     &             0.5D0*iota_Avg*PHI_EDGE(K,GP,LED,pa,J)*NY*WEGP(GP,pa)
#endif
      elseif (EDN(EL) .eq. 4) then
         pp = (pa+1)*(pa+2)/2
         do i = 1,pp
            vtemp1(i) = M_INVM(EL,K,i,pa)
            vtemp2(i) = PHI_EDGE(i,GP,LED,pa,J)
         enddo
         IMASS = dot_product(vtemp1,vtemp2)
#ifdef WAVE_DIF
         HZ(K,1,1,EL) = HZ(K,1,1,EL) - IMASS*XLEN(GED)*
     &                                       0.5D0*ZE_Avg*NX*WEGP(GP,pa)
         HZ(K,2,2,EL) = HZ(K,2,2,EL) - IMASS*XLEN(GED)*
     &                                       0.5D0*ZE_Avg*NY*WEGP(GP,pa)
#endif
         LZ(K,1,1,EL) = LZ(K,1,1,EL) - IMASS*XLEN(GED)*
     &                                       0.5D0*QX_Avg*NX*WEGP(GP,pa)
         LZ(K,1,2,EL) = LZ(K,1,2,EL) - IMASS*XLEN(GED)*
     &                                       0.5D0*QX_Avg*NY*WEGP(GP,pa)
         LZ(K,2,1,EL) = LZ(K,2,1,EL) - IMASS*XLEN(GED)*
     &                                       0.5D0*QY_Avg*NX*WEGP(GP,pa)
         LZ(K,2,2,EL) = LZ(K,2,2,EL) - IMASS*XLEN(GED)*
     &                                       0.5D0*QY_Avg*NY*WEGP(GP,pa)

#ifdef TRACE
         TZ(K,1,1,EL) = TZ(K,1,1,EL)  - IMASS*XLEN(GED)*
     &                                     0.5D0*iota_Avg*NX*WEGP(GP,pa)
         TZ(K,2,2,EL) = TZ(K,2,2,EL) - IMASS*XLEN(GED)*
     &                                     0.5D0*iota_Avg*NY*WEGP(GP,pa)
#endif
      endif
         
      RETURN
      END SUBROUTINE

C***********************************************************************
C     
C     SUBROUTINE EDGE_INT_LDG_sediment()
C     
C     This subroutine computes the edge integrals for the LDG boundary
C     terms for sediment 
C     
C     2012 - cem
C     
C***********************************************************************

      SUBROUTINE EDGE_INT_LDG_sediment
     &     (K,EL,LED,GED,GP,bed_avg,NX,NY,pa,ll)
                                ! <ezpp-noinst>
      
C.....Use appropriate modules

      USE SIZES,ONLY : SZ,layers
      USE GLOBAL,ONLY : AREAS, EDN
      USE DG,ONLY : M_INV,XLEN,PHI_EDGE,WEGP,MZ,M_INVM

      IMPLICIT NONE
      
C.....Declare local variables

      INTEGER K, EL, LED, GED, GP,i,pa,ll,pp,J
      REAL(SZ) AREA, IMASS, bed_avg, NX, NY
      REAL(SZ) vtemp1((pa+1)*(pa+2)/2),vtemp2((pa+1)*(pa+2)/2)
      
C.....Retrieve the element area
      
      AREA = 0.5D0*AREAS(EL)
      J = EDN(EL) - 2
C.....Compute the edge integral
      if (EDN(EL) .eq. 3) then
         IMASS = M_INV(K,pa)/(0.5D0*AREA)
         MZ(K,1,ll,EL) = MZ(K,1,ll,EL) - IMASS*
     &    XLEN(GED)*0.5d0*bed_avg*PHI_EDGE(K,GP,LED,pa,J)*NX*WEGP(GP,pa)
         MZ(K,2,ll,EL) = MZ(K,2,ll,EL) - IMASS*
     &    XLEN(GED)*0.5d0*bed_avg*PHI_EDGE(K,GP,LED,pa,J)*NY*WEGP(GP,pa)
      
      elseif (EDN(EL) .eq. 4) then
         pp = (pa+1)*(pa+2)/2
         do i = 1,pp
            vtemp1(i) = M_INVM(EL,K,i,pa)
            vtemp2(i) = PHI_EDGE(i,GP,LED,pa,J)
         enddo
         IMASS = dot_product(vtemp1,vtemp2)
         MZ(K,1,ll,EL) = MZ(K,1,ll,EL) - IMASS*
     &                            XLEN(GED)*0.5d0*bed_avg*NX*WEGP(GP,pa)
         MZ(K,2,ll,EL) = MZ(K,2,ll,EL) - IMASS*
     &                            XLEN(GED)*0.5d0*bed_avg*NY*WEGP(GP,pa)
      endif
      
      RETURN
      END SUBROUTINE

C***********************************************************************
C     
C     SUBROUTINE RHS_LDG_HYDRO()
C     
C     This subroutine computes the area integrals for the LDG hydro and
C     adds them into the LZ (+sediment diffusion)
C     
C     Written by Shintaro Bunya (01-02-2007)
C     01-10-2011 - cem - adapted for p_enrichment and multicomponent
C     2012 - cem - added sediment layers
C***********************************************************************

      SUBROUTINE RHS_LDG_HYDRO()
      
C.....Use appropriate modules

      USE SIZES,ONLY : layers,SZ
      USE GLOBAL,ONLY : NE,NM,N1,N2,N3,N4,pdg_el,entrop,EDN
      USE DG,ONLY : QX_IN,QY_IN,LZ,M_INV,SFAC_ELEM,ZE_IN,
     &     DRPHI,DSPHI,DRDX,DRDY,DSDX,DSDY,mz,hz,M_INVM,
     &     WAGP,NAGP,DOFS,QX,QY,PHI_AREA,IRK,bed_IN,iota_IN,
     &     slimit1,slimit2,slimit5,balance,DRDXQ,DRDYQ,DSDXQ,DSDYQ,
     &     ze,iota,iota2,slimit3,slimit4,TZ,XAGP,YAGP,Q1,Q2,Q3

      IMPLICIT NONE
      
C.....Declare local variables
      INTEGER L,K,I,pa,ll,kk,jj
      real(sz) ze_sensor1, qx_sensor1, qy_sensor1
      real(sz) ze_sensor2, qx_sensor2, qy_sensor2
      real(sz) iota_sensor1,iota_sensor2,Jac
      logical fl1,fl2,fl3,fl4,fl5
      real(sz), allocatable :: vtemp1(:), vtemp2(:), MTEMPI(:,:)
#ifdef SED_LAY
      real(sz) tbed_sensor1,tbed_sensor2
#endif

      DO L=1,NE

         pa = PDG_EL(L)

#ifdef P0
         if (pa.eq.0) then
            pa = 1
         endif
#endif
         
C.....Retrieve the global node numbers for the element
         
         N1 = NM(L,1)
         N2 = NM(L,2)
         N3 = NM(L,3)
         if (EDN(L) .eq. 4) then
            allocate( MTEMPI(DOFS(L),DOFS(L)), vtemp1(DOFS(L)), 
     &                                                 vtemp2(DOFS(L)) )
            N4 = NM(L,4)
         endif
         
         jj = EDN(L) - 2
         
C.....Compute ZE, QX, QY, and HB at each area Gauss quadrature point

         DO I = 1,NAGP(pa,jj)
            
#ifdef WAVE_DIF
            ZE_IN = 0.D0
#endif
            QX_IN = 0.D0
            QY_IN = 0.D0
#ifdef TRACE
            iota_IN = 0.D0
#endif

#ifdef SED_LAY
            do ll=1,layers
               bed_IN(ll) = 0.D0
            enddo
#endif
            
            DO K = 1,DOFS(L)
#ifdef WAVE_DIF
               ZE_IN = ZE_IN + ZE(K,L,IRK)*PHI_AREA(K,I,pa,jj)
#endif
               QX_IN = QX_IN + QX(K,L,IRK)*PHI_AREA(K,I,pa,jj)
               QY_IN = QY_IN + QY(K,L,IRK)*PHI_AREA(K,I,pa,jj)
#ifdef TRACE
               iota_IN = iota_IN + iota(K,L,IRK)*PHI_AREA(K,I,pa,jj)
#endif

#ifdef SED_LAY
               do ll=1,layers
               bed_IN(ll)=bed_IN(ll)+bed(K,L,IRK,ll)*PHI_AREA(K,I,pa,jj)
               enddo
#endif

            ENDDO

C.....Build the rhs
C........Perform matrix vector multiplication for quads
            if (EDN(L) .eq. 4) then
               do k = 1,dofs(L)
                  do kk = 1,dofs(L)

                     Jac = Q1(L) + Q2(L)*XAGP(I,pa,2) +
     &                                              Q3(L)*YAGP(I,pa,2)
                     MTEMPI(kk,k) = Jac*M_INVM(L,kk,k,pa)
                  enddo
                  vtemp1(k) = DRPHI(k,I,pa,2)
                  vtemp2(k) = DSPHI(k,I,pa,2)
               enddo
               vtemp1 = matmul(MTEMPI,vtemp1)
               vtemp2 = matmul(MTEMPI,vtemp2)
            endif
            
            DO K = 1,DOFS(L)
               if (EDN(L) .eq. 3) then
#ifdef WAVE_DIF
                  HZ(K,1,1,L) = HZ(K,1,1,L)
     &            + M_INV(K,pa)*ZE_IN*SFAC_ELEM(I,L,pa)* ! <--- ZE/Mk
     &            ( DRPHI(K,I,pa,1)*DRDX(L) + DSPHI(K,I,pa,1)*DSDX(L) )* ! <--- dphi/dx
     &            WAGP(I,pa,1)  ! <--- weight
                  HZ(K,2,2,L) = HZ(K,2,2,L)
     &            + M_INV(K,pa)*ZE_IN* ! <--- ZE/Mk
     &            ( DRPHI(K,I,pa,1)*DRDY(L) + DSPHI(K,I,pa,1)*DSDY(L) )* ! <--- dphi/dy
     &            WAGP(I,pa,1)  ! <--- weight
#endif

                  LZ(K,1,1,L) = LZ(K,1,1,L)
     &            + M_INV(K,pa)*QX_IN*SFAC_ELEM(I,L,pa)* ! <--- QX/Mk
     &            ( DRPHI(K,I,pa,1)*DRDX(L) + DSPHI(K,I,pa,1)*DSDX(L) )* ! <--- dphi/dx
     &            WAGP(I,pa,1)  ! <--- weight

                  LZ(K,1,2,L) = LZ(K,1,2,L)
     &            + M_INV(K,pa)*QX_IN* ! <--- QX/Mk
     &            ( DRPHI(K,I,pa,1)*DRDY(L) + DSPHI(K,I,pa,1)*DSDY(L) )* ! <--- dphi/dy
     &            WAGP(I,pa,1)  ! <--- weight

                  LZ(K,2,1,L) = LZ(K,2,1,L)
     &            + M_INV(K,pa)*QY_IN*SFAC_ELEM(I,L,pa)* ! <--- QY/Mk
     &            ( DRPHI(K,I,pa,1)*DRDX(L) + DSPHI(K,I,pa,1)*DSDX(L) )* ! <--- dphi/dx
     &            WAGP(I,pa,1)  ! <--- weight

                  LZ(K,2,2,L) = LZ(K,2,2,L)
     &            + M_INV(K,pa)*QY_IN* ! <--- QY/Mk
     &            ( DRPHI(K,I,pa,1)*DRDY(L) + DSPHI(K,I,pa,1)*DSDY(L) )* ! <--- dphi/dy
     &            WAGP(I,pa,1)  ! <--- weight
#ifdef TRACE
                  TZ(K,1,1,L) = TZ(K,1,1,L)
     &            + M_INV(K,pa)*iota_IN*SFAC_ELEM(I,L,pa)* ! <--- iota*H/Mk
     &            ( DRPHI(K,I,pa,1)*DRDX(L) + DSPHI(K,I,pa,1)*DSDX(L) )* ! <--- dphi/dx
     &            WAGP(I,pa,1)  ! <--- weight
                  TZ(K,2,2,L) = TZ(K,2,2,L)
     &            + M_INV(K,pa)*iota_IN* ! <--- iota*H/Mk
     &            ( DRPHI(K,I,pa,1)*DRDY(L) + DSPHI(K,I,pa,1)*DSDY(L) )* ! <--- dphi/dy
     &            WAGP(I,pa,1)  ! <--- weight
#endif

#ifdef SED_LAY
                  do ll=1,layers
                     MZ(K,1,ll,L) =  MZ(K,1,ll,L)
     &            + M_INV(K,pa)*bed_IN(ll)*SFAC_ELEM(I,L,pa)* 
     &            ( DRPHI(K,I,pa,1)*DRDX(L) + DSPHI(K,I,pa,1)*DSDX(L) )*
     &            WAGP(I,pa,1) 
                     MZ(K,2,ll,L) =  MZ(K,2,ll,L)
     &            + M_INV(K,pa)*bed_IN(ll)* 
     &            ( DRPHI(K,I,pa,1)*DRDY(L) + DSPHI(K,I,pa,1)*DSDY(L) )*
     &            WAGP(I,pa,1) 
                  enddo
#endif
               elseif (EDN(L) .eq. 4) then
#ifdef WAVE_DIF
                  HZ(K,1,1,L) = HZ(K,1,1,L)
     &               + ZE_IN*SFAC_ELEM(I,L,pa)*WAGP(I,pa,2)*
     &               (vtemp1(k)*DRDXQ(L,I,pa) + vtemp2(k)*DSDXQ(L,I,pa))
     
                  HZ(K,2,2,L) = HZ(K,2,2,L)
     &               + ZE_IN*WAGP(I,pa,2)*
     &               (vtemp1(k)*DRDYQ(L,I,pa) + vtemp2(k)*DSDYQ(L,I,pa))
#endif

                  LZ(K,1,1,L) = LZ(K,1,1,L)
     &               + QX_IN*SFAC_ELEM(I,L,pa)*WAGP(I,pa,2)*
     &               (vtemp1(k)*DRDXQ(L,I,pa) + vtemp2(k)*DSDXQ(L,I,pa))

                  LZ(K,1,2,L) = LZ(K,1,2,L)
     &               + QX_IN*WAGP(I,pa,2)* 
     &               (vtemp1(k)*DRDYQ(L,I,pa) + vtemp2(k)*DSDYQ(L,I,pa))

                  LZ(K,2,1,L) = LZ(K,2,1,L)
     &               + QY_IN*SFAC_ELEM(I,L,pa)*WAGP(I,pa,2)* 
     &               (vtemp1(k)*DRDXQ(L,I,pa) + vtemp2(k)*DSDXQ(L,I,pa))

                  LZ(K,2,2,L) = LZ(K,2,2,L)
     &               + QY_IN*WAGP(I,pa,2)* 
     &               (vtemp1(k)*DRDYQ(L,I,pa) + vtemp2(k)*DSDYQ(L,I,pa))
#ifdef TRACE
                  TZ(K,1,1,L) = TZ(K,1,1,L)
     &               + iota_IN*SFAC_ELEM(I,L,pa)*WAGP(I,pa,2)*
     &               (vtemp1(k)*DRDXQ(L,I,pa) + vtemp2(k)*DSDXQ(L,I,pa))
     
                  TZ(K,2,2,L) = TZ(K,2,2,L)
     &               + iota_IN*WAGP(I,pa,2)*
     &               (vtemp1(k)*DRDYQ(L,I,pa) + vtemp2(k)*DSDYQ(L,I,pa))
#endif

#ifdef SED_LAY
                  do ll=1,layers
                     MZ(K,1,ll,L) =  MZ(K,1,ll,L)
     &               + bed_IN(ll)*SFAC_ELEM(I,L,pa)*WAGP(I,pa,2)*
     &               (vtemp1(k)*DRDXQ(L,I,pa) + vtemp2(k)*DSDXQ(L,I,pa))
                      
                     MZ(K,2,ll,L) =  MZ(K,2,ll,L)
     &               + bed_IN(ll)*WAGP(I,pa,2)*
     &               (vtemp1(k)*DRDYQ(L,I,pa) + vtemp2(k)*DSDYQ(L,I,pa))
                  enddo
#endif
               endif
            enddo
         enddo

#ifdef ARTDIF

         slimit1 = -100.D0
         slimit2 = -100.D0
         slimit3 = -100.D0

         entrop(:,L) = -100.D0
#endif
           
#ifdef WAVE_DIF 
         ze_sensor1 = 0.d0
         ze_sensor2 = 0.d0  
#endif             
         qx_sensor1 = 0.d0
         qy_sensor1 = 0.d0
           
         qx_sensor2 = 0.d0
         qy_sensor2 = 0.d0

#ifdef TRACE
         slimit4 = -100.D0 
         iota_sensor1 = 0.d0
         iota_sensor2 = 0.d0
#endif
         
#ifdef SED_LAY
         slimit5 = -100.D0
         tbed_sensor1 = 0.d0
         tbed_sensor2 = 0.d0
#endif
         
         do I = 1,NAGP(pa,jj)
            
            do kk = 2, dofs(L)  !Compute the first sensor
#ifdef WAVE_DIF 
               ze_sensor1 = ze_sensor1 + (ze(kk,L,irk)* 
     &              phi_area(kk,I,pa,jj))**2.D0 * wagp(I,pa,jj)
#endif
               qx_sensor1 = qx_sensor1 + (qx(kk,L,irk)* 
     &              phi_area(kk,I,pa,jj))**2.D0 * wagp(I,pa,jj)
               qy_sensor1 = qy_sensor1 + (qy(kk,L,irk)* 
     &              phi_area(kk,I,pa,jj))**2.D0 * wagp(I,pa,jj)
#ifdef TRACE
               iota_sensor1 = iota_sensor1 + (iota(kk,L,irk)
     &              * phi_area(kk,I,pa,jj))**2.D0 * wagp(I,pa,jj)
#endif
               
#ifdef SED_LAY       
               do ll=1,layers
                  tbed_sensor1 = tbed_sensor1 + (bed(kk,L,irk,ll)* 
     &                 phi_area(kk,I,pa,jj))**2.D0 * wagp(I,pa,jj)
               enddo
#endif

            enddo
            
            do kk = 1,dofs(L) !Compute the second sensor
#ifdef WAVE_DIF 
               ze_sensor2 = ze_sensor2 + (ze(kk,L,irk)* 
     &              phi_area(kk,I,pa,jj))**2.D0 * wagp(I,pa,jj)
#endif
               qx_sensor2 = qx_sensor2 + (qx(kk,L,irk)* 
     &              phi_area(kk,I,pa,jj))**2.D0 * wagp(I,pa,jj)
               qy_sensor2 = qy_sensor2 + (qy(kk,L,irk)* 
     &              phi_area(kk,I,pa,jj))**2.D0 * wagp(I,pa,jj)
#ifdef TRACE
               iota_sensor2 = iota_sensor2 + (iota(kk,L,irk)
     &              * phi_area(kk,I,pa,jj))**2.D0 * wagp(I,pa,jj)
#endif
               
#ifdef SED_LAY       
               do ll=1,layers
                  tbed_sensor2 = tbed_sensor2 + (bed(kk,L,irk,ll)* !Adjust sensor for multilayers!
     &                 phi_area(kk,I,pa,jj))**2.D0 * wagp(I,pa,jj)
               enddo
#endif

            enddo
            
         enddo
#ifdef WAVE_DIF 
         if  (ze_sensor2.gt.1.0e-4.and.ze_sensor1.gt.1.0e-4 ) then
            slimit1 = log10( ze_sensor1/ze_sensor2 ) + balance(1)
            entrop(1,L) = slimit1
         endif
#endif
         if (qx_sensor2.gt.1.0e-4.and.qx_sensor1.gt.1.0e-4 ) then
            slimit2 = log10( qx_sensor1/qx_sensor2 ) + balance(2)
            entrop(2,L) = slimit2
         endif
         if ( qy_sensor2.gt.1.0e-4.and.qy_sensor1.gt.1.0e-4 ) then
            slimit3 = log10( qy_sensor1/qy_sensor2 ) + balance(3)
            entrop(3,L) = slimit3
         endif
#ifdef TRACE
          if ( iota_sensor2.gt.1.0e-4.and.iota_sensor1.gt.1.0e-4 ) then
            slimit4 = log10( iota_sensor1/iota_sensor2 ) + balance(4)
            entrop(4,L) = slimit4
          endif
#endif
         
#ifdef SED_LAY
         if ( tbed_sensor2.gt.1.0e-4.and.tbed_sensor1.gt.1.0e-4 ) then
            slimit5 = log10( tbed_sensor1/tbed_sensor2 ) + balance(5)
            entrop(5,L) = slimit5 
        endif
#endif                
 
         if (EDN(L) .eq. 4) deallocate( MTEMPI, vtemp1, vtemp2 )

      ENDDO

      RETURN
      END SUBROUTINE

