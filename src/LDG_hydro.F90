!***********************************************************************
!     
!     SUBROUTINE LDG_HYDRO(IT)
!     
!     Compute variable Z to be used in the LDG terms.
!     
!     Written by Shintaro Bunya (01-01-2007)
!
!     01-10-2011 - cem - adapted for p_enrichment and multicomponent
!     06-01-2012 - cem - sediment diffusion added
!     
!***********************************************************************

      SUBROUTINE LDG_HYDRO(s,dg_here,global_here,nodalattr_here,IT)

!.....Use appropriate modules
      
      USE SIZES
      USE GLOBAL
      USE DG
      USE NodalAttributes

#ifdef CMPI
      USE MESSENGER_ELEM,ONLY : UPDATELZ_ELEM,UPDATEMZ_ELEM
#endif

      IMPLICIT NONE
      
      type (sizes_type) :: s
      type (dg_type) :: dg_here
      type (global_type) :: global_here
      type (nodalattr_type) :: nodalattr_here
      
!.....Declare local variables

      INTEGER IT,L,GED,NBOREL,NNBORS,NDRYNBORS,k,i,ll
      INTEGER Detected
      REAL(SZ) QBCT1,QBCT2,QBC1,QBC2,ZP(3)
      REAL(SZ) DPAVG,lim_by
      Real(SZ),Allocatable :: tmp_mz(:,:,:,:)

      Allocate ( tmp_mz(dg_here%dofh,2,1,s%MNE) )

!.....Initialize for viscosity/diffusion
#ifdef WAVE_DIF
      dg_here%HZ = 0.D0
#endif 
      dg_here%LZ = 0.D0
#ifdef TRACE
      dg_here%TZ = 0.d0
#endif
#ifdef SED_LAY
      dg_here%MZ = 0.d0
#endif
     
!.....Compute elevation specified edges

      IF (dg_here%NEEDS.GT.0)  CALL OCEAN_EDGE_LDG_HYDRO(s,dg_here,global_here)

!.....Compute no-normal flow edges

      IF (dg_here%NLEDS.GT.0)  CALL LAND_EDGE_LDG_HYDRO(s,dg_here,global_here)

!.....Compute non-zero flow edges

      IF (dg_here%NFEDS.GT.0)  CALL FLOW_EDGE_LDG_HYDRO(s,dg_here,global_here)
      
!.....Compute radiation edges

      IF (dg_here%NREDS.GT.0)  CALL RADIATION_EDGE_LDG_HYDRO(s,dg_here,global_here)

!.....Compute internal edges

      CALL INTERNAL_EDGE_LDG_HYDRO(s,dg_here,global_here)

!.....Loop over interior elements

      CALL RHS_LDG_HYDRO(s,dg_here,global_here)

      do L=1,global_here%NE

         global_here%N1 = global_here%NM(L,1)
         global_here%N2 = global_here%NM(L,2)
         global_here%N3 = global_here%NM(L,3)
         if (global_here%EVMSUM.NE.0.D0) then
            dg_here%EVMAvg = (NODALATTR_HERE%EVM(global_here%N1)+NODALATTR_HERE%EVM(global_here%N2)+NODALATTR_HERE%EVM(global_here%N3))/3.d0
         else
            dg_here%EVMAvg = 0.D0
         endif

         dg_here%bg_dif = 0.D0 
         dg_here%trc_dif = 0.D0 
          
#ifdef SED_LAY
         dg_here%SEVDMAvg = dg_here%SEVDM 
#endif
 
#ifdef ARTDIF

         lim_by = 1.D0
         
#ifdef WAVE_DIF
         if (global_here%entrop(1,L).gt.dg_here%s0+dg_here%kappa) then
            
            dg_here%HZ(:,:,:,L) = (dg_here%bg_dif+dg_here%e1(1)) *&
           (abs(lim_by+dg_here%balance(1)))**dg_here%pa * dg_here%HZ(:,:,:,L)

         else

            dg_here%HZ(:,:,:,L) = dg_here%bg_dif*dg_here%HZ(:,:,:,L)

         endif
#endif

         if (global_here%entrop(2,L).gt.dg_here%s0+dg_here%kappa) then
            
            dg_here%LZ(:,1,1,L) = (dg_here%EVMAvg+dg_here%e1(2) * &
           (abs(lim_by+dg_here%balance(2)))**dg_here%pa)  * dg_here%LZ(:,1,1,L)
            dg_here%LZ(:,1,2,L) = (dg_here%EVMAvg+dg_here%e1(2) *&
           (abs(lim_by+dg_here%balance(2)))**dg_here%pa ) * dg_here%LZ(:,1,2,L)
            dg_here%LZ(:,2,1,L) = (dg_here%EVMAvg+dg_here%e1(2) *&
           (abs(lim_by+dg_here%balance(2)))**dg_here%pa ) * dg_here%LZ(:,2,1,L)
     
         else

            dg_here%LZ(:,1,:,L) = dg_here%EVMAvg* dg_here%LZ(:,1,:,L)
            dg_here%LZ(:,2,1,L) = dg_here%EVMAvg* dg_here%LZ(:,2,1,L)
            
         endif

         if (global_here%entrop(3,L).gt.dg_here%s0+dg_here%kappa) then
            
            dg_here%LZ(:,2,2,L) = (dg_here%EVMAvg+dg_here%e1(3) * &
           (abs(lim_by+dg_here%balance(3)))**dg_here%pa) * dg_here%LZ(:,2,2,L)

            if (global_here%entrop(2,L).global_here%le.dg_here%s0+dg_here%kappa) then
               
               dg_here%LZ(:,2,1,L) = (dg_here%EVMAvg+dg_here%e1(3) *&
              (abs(lim_by+dg_here%balance(3)))**dg_here%pa ) * dg_here%LZ(:,2,1,L)
               dg_here%LZ(:,1,2,L) = (dg_here%EVMAvg+dg_here%e1(3) *&
              (abs(lim_by+dg_here%balance(3)))**dg_here%pa ) * dg_here%LZ(:,1,2,L)
               
            endif
            
         else
            
            dg_here%LZ(:,2,2,L) = dg_here%EVMAvg* dg_here%LZ(:,2,2,L)
            
            if (global_here%entrop(2,L).global_here%le.dg_here%s0+dg_here%kappa) then
               
               dg_here%LZ(:,1,2,L) = dg_here%EVMAvg* dg_here%LZ(:,1,2,L)
               dg_here%LZ(:,2,1,L) = dg_here%EVMAvg* dg_here%LZ(:,2,1,L)
               
            endif
            
         endif

#ifdef TRACE
         if (global_here%entrop(4,L).gt.dg_here%s0+dg_here%kappa) then
            
            dg_here%TZ(:,:,:,L) = (dg_here%trc_dif+dg_here%e1(4)) *&
           (abs(lim_by+dg_here%balance(4)))**dg_here%pa * dg_here%TZ(:,:,:,L)

         else

            dg_here%TZ(:,:,:,L) = dg_here%trc_dif*dg_here%TZ(:,:,:,L)

         endif
#endif

#ifdef SED_LAY

         if (global_here%entrop(5,L).gt.dg_here%s0+dg_here%kappa) then
            
            dg_here%MZ(:,:,:,L) = (dg_here%SEVDMAvg+dg_here%e1(5)) *&
           (abs(lim_by+dg_here%balance(5)))**dg_here%pa * dg_here%MZ(:,:,:,L)

         else

            dg_here%MZ(:,:,:,L) = dg_here%MZ(:,:,:,L)*dg_here%SEVDMAvg
            
         endif
#endif
         
#else
#ifdef WAVE_DIF
         dg_here%HZ(:,:,:,L) = dg_here%bg_dif*dg_here%HZ(:,:,:,L)
#endif
         dg_here%LZ(:,:,:,L) = dg_here%EVMAvg*dg_here%LZ(:,:,:,L)
               
#ifdef TRACE
         dg_here%TZ(:,:,:,L) = dg_here%trc_dif*dg_here%TZ(:,:,:,L)
#endif

#ifdef SED_LAY
         dg_here%MZ(:,:,:,L) = dg_here%MZ(:,:,:,L)*dg_here%SEVDMAvg
#endif


#endif

      enddo
!********************************************

#ifdef CMPI
#ifdef WAVE_DIF
      CALL UPDATELZ_ELEM(DG_HERE,dg_here%HZ)
#endif
      CALL UPDATELZ_ELEM(DG_HERE,dg_here%LZ)
#ifdef TRACE
      CALL UPDATELZ_ELEM(DG_HERE,dg_here%TZ)
#endif
#ifdef SED_LAY
      do ll=1,layers
         tmp_mz(:,:,1,:) = dg_here%MZ(:,:,ll,:) 
         CALL UPDATEMZ_ELEM(tmp_mz)
         dg_here%MZ(:,:,ll,:) = tmp_mz(:,:,1,:)
      enddo
#endif
#endif

      RETURN
    END SUBROUTINE LDG_HYDRO

!***********************************************************************
!     
!     SUBROUTINE OCEAN_EDGE_LDG_HYDRO( )
!     
!     This subroutine does the following:
!     
!     1.  Calculates the values of the necessary variables at the global_here%edge
!     gauss points for NO-NORMAL FLOW edges
!     2.  Compute the LDG flux at these points (sediment diffusion?).
!     3.  Compute the boundary integrals.
!     
!     Written by Shintaro Bunya (01-04-2007)
!     01-10-2011 - cem - adapted for p_enrichment and multicomponent
!     
!***********************************************************************
      SUBROUTINE OCEAN_EDGE_LDG_HYDRO(s,dg_here,global_here)

!.....Use appropriate modules

      USE SIZES
      USE GLOBAL
      USE DG

      IMPLICIT NONE

      type (sizes_type) :: s
      type (dg_type) :: dg_here
      type (global_type) :: global_here

!.....Declare local variables

      INTEGER L,LED,GED,k,I,ll
      REAL(SZ) QX_AVG, QY_AVG,ZE_AVG,bed_AVG(s%layers)
      REAL(SZ) iota_AVG

      DO 1000 L = 1, dg_here%needs
         
!.....Retrieve the global and local global_here%edge number

         GED = dg_here%NEEDN(L)
         LED = dg_here%NEDSD(1,GED)

!.....Retrieve the elements which share the global_here%edge

         global_here%EL_IN = dg_here%NEDEL(1,GED)

         dg_here%pa = global_here%PDG_EL(global_here%EL_IN)

#ifdef P0
         if (dg_here%pa.eq.0) then
            dg_here%pa = 1
         endif
#endif
         
!.....Retrieve the components of the normal vector to the global_here%edge
         
         dg_here%NX = dg_here%COSNX(GED)
         dg_here%NY = dg_here%SINNX(GED)
         
!.....Retrieve the nodes of the global_here%edge
         
         global_here%N1 = dg_here%NEDNO(1,GED)
         global_here%N2 = dg_here%NEDNO(2,GED)
         
!.....Compute dg_here%ZE, dg_here%QX, dg_here%QY, and dg_here%HB at each global_here%edge Gauss quadrature point

         DO I = 1,dg_here%NEGP(dg_here%pa)

#ifdef WAVE_DIF
            dg_here%ZE_IN = 0.D0
            dg_here%ZE_EX = 0.D0
#endif
            dg_here%QX_IN = 0.D0
            dg_here%QY_IN = 0.D0

            dg_here%QX_EX = 0.D0
            dg_here%QY_EX = 0.D0

#ifdef TRACE
            dg_here%iota_IN = 0.D0
            dg_here%iota_EX = 0.D0
#endif

            !deal with sediment
#ifdef SED_LAY
             do ll=1,s%layers
               dg_here%bed_IN(ll) = 0.D0
               dg_here%bed_EX(ll) = 0.D0
            enddo
#endif
           
!.....Compute the solution at the interior state

            DO K = 1,dg_here%DOFS(global_here%EL_IN)
#ifdef WAVE_DIF
               dg_here%ZE_IN = dg_here%ZE_IN + dg_here%ZE(K,global_here%EL_IN,dg_here%IRK)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
#endif
               dg_here%QX_IN = dg_here%QX_IN + dg_here%QX(K,global_here%EL_IN,dg_here%IRK)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
               dg_here%QY_IN = dg_here%QY_IN + dg_here%QY(K,global_here%EL_IN,dg_here%IRK)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
#ifdef TRACE
               dg_here%iota_IN = dg_here%iota_IN + dg_here%iota(K,global_here%EL_IN,dg_here%IRK)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
#endif
#ifdef SED_LAY
               do ll=1,s%layers
                  dg_here%bed_IN(ll) = dg_here%bed_IN(ll) &
                 + dg_here%bed(K,global_here%EL_IN,dg_here%IRK,ll)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
               enddo
#endif
            ENDDO
            
!.....Set the exterior state flows equal to the interior state flows

            dg_here%QX_EX = dg_here%QX_IN
            dg_here%QY_EX = dg_here%QY_IN

#ifdef SED_LAY
            do ll=1,s%layers
               dg_here%bed_EX(ll) = dg_here%bed_IN(ll)
            enddo
#endif

!.....Take the average
#ifdef WAVE_DIF
            ZE_AVG = dg_here%ZE_IN
#endif
            QX_AVG = 0.5D0*(dg_here%QX_IN + dg_here%QX_EX)*dg_here%SFACED(I,LED,global_here%EL_IN,dg_here%pa)
            QY_AVG = 0.5D0*(dg_here%QY_IN + dg_here%QY_EX)
#ifdef TRACE
            iota_AVG = dg_here%iota_IN
#endif
#ifdef SED_LAY
            do ll=1,s%layers
               bed_AVG(ll) = 0.5D0*( dg_here%bed_IN(ll) + dg_here%bed_EX(ll) )
            enddo
#endif

!.....Compute the global_here%edge integral
            
            DO K = 1,dg_here%DOFS(global_here%EL_IN)
               CALL EDGE_INT_LDG_HYDRO&
              (s,dg_here,global_here,K,global_here%EL_IN,LED,GED,I,iota_AVG,ZE_AVG,QX_AVG,QY_AVG,dg_here%NX,dg_here%NY,dg_here%pa)
#ifdef SED_LAY
               do ll=1,s%layers
                  CALL EDGE_INT_LDG_sediment&
                (K,global_here%EL_IN,LED,GED,I,bed_AVG(ll),dg_here%NX,dg_here%NY,dg_here%pa,ll)
               enddo
#endif
            ENDDO
            
         ENDDO

 1000 CONTINUE
      
      RETURN
      END SUBROUTINE


!***********************************************************************
!     
!     SUBROUTINE LAND_EDGE_HYDRO_LDG( )
!     
!     This subroutine does the following:
!     
!     1.  Calculates the values of the necessary variables at the global_here%edge
!     gauss points for NO-NORMAL FLOW edges
!     2.  Compute the LDG flux at these points.
!     3.  Compute the boundary integrals.
!     
!     Written by Shintaro Bunya (01-02-2007)
!     01-10-2011 - cem - adapted for p_enrichment and multicomponent
!     
!***********************************************************************

      SUBROUTINE LAND_EDGE_LDG_HYDRO(s,dg_here,global_here)

!.....Use appropriate modules

      USE SIZES
      USE GLOBAL
      USE DG
      IMPLICIT NONE

      type (sizes_type) :: s
      type (dg_type) :: dg_here
      type (global_type) :: global_here

!.....Declare local variables

      INTEGER K,L,LED,GED,GP,i,kk,ll
      REAL(SZ) AREA, IMASS
      REAL(SZ) TX, TY, QX_AVG, QY_AVG, bed_AVG(s%layers)
      REAL(SZ) ZE_AVG,iota_AVG
      
      DO 1000 L = 1,dg_here%NLEDS
         
!.....Retrieve the global and local global_here%edge number
         
         GED = dg_here%NLEDN(L)
         LED = dg_here%NEDSD(1,GED)
         
!.....Retrieve the elements which share the global_here%edge
         
         global_here%EL_IN = dg_here%NEDEL(1,GED)
         
         dg_here%pa = global_here%PDG_EL(global_here%EL_IN)
         
#ifdef P0
         if (dg_here%pa.eq.0) then
            dg_here%pa = 1
         endif
#endif

!.....Retrieve the components of the normal vector to the global_here%edge
         
         dg_here%NX = dg_here%COSNX(GED)
         dg_here%NY = dg_here%SINNX(GED)
         
!.....Set the components for the tangential vector to the global_here%edge
         
         TX = -dg_here%NY
         TY =  dg_here%NX

!.....Compute dg_here%ZE, dg_here%QX, dg_here%QY, and dg_here%HB at each Gauss point

         DO I = 1,dg_here%NEGP(dg_here%pa)
#ifdef WAVE_DIF
            dg_here%ZE_IN = 0.D0
#endif
            dg_here%QX_IN = 0.D0
            dg_here%QY_IN = 0.D0
#ifdef TRACE
            dg_here%iota_IN = 0.d0
#endif
#ifdef SED_LAY
            do ll=1,s%layers
               dg_here%bed_IN(ll) = 0.D0
            enddo
#endif

!.....Compute the solution at the interior state

            DO K = 1,dg_here%DOFS(global_here%EL_IN)
#ifdef WAVE_DIF
               dg_here%ZE_IN = dg_here%ZE_IN + dg_here%ZE(K,global_here%EL_IN,dg_here%IRK)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
#endif
               dg_here%QX_IN = dg_here%QX_IN + dg_here%QX(K,global_here%EL_IN,dg_here%IRK)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
               dg_here%QY_IN = dg_here%QY_IN + dg_here%QY(K,global_here%EL_IN,dg_here%IRK)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
#ifdef TRACE
               dg_here%iota_IN = dg_here%iota_IN + dg_here%iota(K,global_here%EL_IN,dg_here%IRK)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
#endif
#ifdef SED_LAY
               do ll=1,s%layers
                  dg_here%bed_IN(ll) = dg_here%bed_IN(ll) + dg_here%bed(K,global_here%EL_IN,dg_here%IRK,ll)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
               enddo
#endif

            ENDDO
            
!.....Compute the velocity in the normal and tangental direction
            
            global_here%Q_N_INT = dg_here%QX_IN*dg_here%NX + dg_here%QY_IN*dg_here%NY
            global_here%Q_T_INT = dg_here%QX_IN*TX + dg_here%QY_IN*TY

!.....Reflect the velocity in the normal direction

            global_here%Q_N_EXT = -global_here%Q_N_INT
            global_here%Q_T_EXT =  global_here%Q_T_INT
            
!.....Compute the global_here%x and global_here%y components of the external state flow

            dg_here%QX_EX = ( TY*global_here%Q_N_EXT - dg_here%NY*global_here%Q_T_EXT)/(dg_here%NX*TY - dg_here%NY*TX)
            dg_here%QY_EX = (-TX*global_here%Q_N_EXT + dg_here%NX*global_here%Q_T_EXT)/(dg_here%NX*TY - dg_here%NY*TX)

#ifdef SED_LAY
            do ll=1,s%layers !maybe this should be changed?
               dg_here%bed_EX(ll) =  dg_here%bed_IN(ll)
            enddo
#endif            


!.....Take the average
#ifdef WAVE_DIF
            ZE_AVG = dg_here%ZE_IN
#endif
            QX_AVG = 0.5D0*(dg_here%QX_IN + dg_here%QX_EX)*dg_here%SFACED(I,LED,global_here%EL_IN,dg_here%pa)
            QY_AVG = 0.5D0*(dg_here%QY_IN + dg_here%QY_EX)
#ifdef TRACE
            iota_AVG = dg_here%iota_IN
#endif
#ifdef SED_LAY
            do ll=1,s%layers !take the Fronenius norm
               bed_AVG(ll) = 0.5D0*( dg_here%bed_IN(ll) + dg_here%bed_EX(ll) )
            enddo
#endif

!.....Compute the global_here%edge integral

            DO K = 1,dg_here%DOFS(global_here%EL_IN)
               CALL EDGE_INT_LDG_HYDRO&
              (s,dg_here,global_here,K,global_here%EL_IN,LED,GED,I,iota_AVG,ZE_AVG,QX_AVG,QY_AVG,dg_here%NX,dg_here%NY,dg_here%pa)
#ifdef SED_LAY
               do ll=1,s%layers
                  CALL EDGE_INT_LDG_sediment&
                 (K,global_here%EL_IN,LED,GED,I,bed_AVG(ll),dg_here%NX,dg_here%NY,dg_here%pa,ll)
               enddo
#endif
            ENDDO

         ENDDO
            
 1000 CONTINUE
      
      RETURN
      END SUBROUTINE


!***********************************************************************
!     
!     SUBROUTINE FLOW_EDGE_LDG_HYDRO( )
!     
!     This subroutine does the following:
!     
!     1.  Calculates the values of the necessary variables at the global_here%edge
!     gauss points for INTERNAL edges
!     2.  Take the average of these values (sediment diffusion)
!     3.  Perform boundary integration
!     
!     Written by Shintaro Bunya (01-05-2007)
!     01-10-2011 - cem - adapted for p_enrichment and multicomponent
!     
!***********************************************************************

      SUBROUTINE FLOW_EDGE_LDG_HYDRO(s,dg_here,global_here)

!.....Use appropriate modules

      USE SIZES
      USE GLOBAL
      USE DG

      IMPLICIT NONE

      type (sizes_type) :: s
      type (dg_type) :: dg_here
      type (global_type) :: global_here

!.....Declare local variables

      INTEGER L, LED, GED,k,i,jj,ll
      REAL(SZ) TX, TY, QX_AVG, QY_AVG, bed_AVG(s%layers)
      REAL(SZ) iota_AVG,ZE_AVG

      DO 1000 L = 1,dg_here%NFEDS
         
!.....Retrieve the global and local global_here%edge number

         GED = dg_here%NFEDN(L)
         LED = dg_here%NEDSD(1,GED)

!.....Retrieve the element to which the global_here%edge belongs

         global_here%EL_IN = dg_here%NEDEL(1,GED)

         dg_here%PA = global_here%PDG_EL(global_here%EL_IN)

#ifdef P0
         if (dg_here%pa.eq.0) then
            dg_here%pa = 1
         endif
#endif
         
!.....Retrieve the components of the normal vector to the global_here%edge
         
         dg_here%NX = dg_here%COSNX(GED)
         dg_here%NY = dg_here%SINNX(GED)
         
!.....Set the components for the tangential vector to the global_here%edge

         TX = -dg_here%NY
         TY =  dg_here%NX
         
!.....Compute dg_here%ZE, dg_here%QX, dg_here%QY, and dg_here%HB at each global_here%edge Gauss quadrature point

         DO I = 1,dg_here%NEGP(dg_here%pa)

#ifdef WAVE_DIF
            dg_here%ZE_IN = 0.D0
            dg_here%ZE_EX = 0.D0
#endif

            dg_here%QX_IN = 0.D0
            dg_here%QY_IN = 0.D0

            dg_here%QX_EX = 0.D0
            dg_here%QY_EX = 0.D0
#ifdef TRACE
            dg_here%iota_IN = 0.D0
            dg_here%iota_EX = 0.D0
#endif
#ifdef SED_LAY
             do ll=1,s%layers
               dg_here%bed_IN(ll) = 0.D0
               dg_here%bed_EX(ll) = 0.D0
            enddo
#endif
            
!.....Compute the specified flow boundaries for the exterior state

            global_here%Q_N_EXT = 0.D0
            DO JJ = 1,global_here%NFFR
               
               IF(global_here%FPER(JJ).EQ.0.D0) THEN
                  global_here%NCYC = 0.D0
               ELSE
                  global_here%NCYC = INT(dg_here%TIMEDG/global_here%FPER(JJ))
               ENDIF

               global_here%ARGJ = global_here%FAMIG(JJ)*(dg_here%TIMEDG - global_here%NCYC*global_here%FPER(JJ)) + global_here%FFACE(JJ)
               global_here%RFF  = global_here%FFF(JJ)*global_here%RAMPExtFlux
               
               dg_here%QNAM_GP = 0.5D0*(dg_here%QNAM_DG(JJ,L,1) + dg_here%QNAM_DG(JJ,L,2))&
              + 0.5D0*(dg_here%QNAM_DG(JJ,L,2) - dg_here%QNAM_DG(JJ,L,1))*dg_here%XEGP(I,dg_here%pa)
               dg_here%QNPH_GP = 0.5D0*(dg_here%QNPH_DG(JJ,L,1) + dg_here%QNPH_DG(JJ,L,2))&
              + 0.5D0*(dg_here%QNPH_DG(JJ,L,2) - dg_here%QNPH_DG(JJ,L,1))*dg_here%XEGP(I,dg_here%pa)
               
               global_here%ARG = global_here%ARGJ - dg_here%QNPH_GP
               
               global_here%Q_N_EXT = global_here%Q_N_EXT + dg_here%QNAM_GP*global_here%RFF*COS(global_here%ARG)
               global_here%Q_T_EXT =  0.D0


            ENDDO
            
!.....Compute the solution at the interior state

            DO K = 1,dg_here%DOFS(global_here%EL_IN)

#ifdef WAVE_DIF
               dg_here%ZE_IN = dg_here%ZE_IN + dg_here%ZE(K,global_here%EL_IN,dg_here%IRK)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
#endif

               dg_here%QX_IN = dg_here%QX_IN + dg_here%QX(K,global_here%EL_IN,dg_here%IRK)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
               dg_here%QY_IN = dg_here%QY_IN + dg_here%QY(K,global_here%EL_IN,dg_here%IRK)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)

#ifdef SED_LAY
               do ll=1,s%layers
                  dg_here%bed_in(ll) = dg_here%bed_in(ll) + dg_here%bed(K,global_here%EL_IN,dg_here%IRK,ll)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
                  dg_here%bed_ex(ll) = dg_here%bed_ex(ll) + dg_here%bed(K,global_here%EL_IN,dg_here%IRK,ll)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
               enddo
#endif
            ENDDO


               dg_here%QX_EX = -( TY*global_here%Q_N_EXT - dg_here%NY*global_here%Q_T_EXT)/(dg_here%NX*TY - dg_here%NY*TX)
               dg_here%QY_EX = -(-TX*global_here%Q_N_EXT + dg_here%NX*global_here%Q_T_EXT)/(dg_here%NX*TY - dg_here%NY*TX)
            
!.....Take the average

#ifdef WAVE_DIF
               ZE_AVG = dg_here%ZE_IN
#endif
               QX_AVG = 0.5*(dg_here%QX_IN + dg_here%QX_EX)*dg_here%SFACED(I,LED,global_here%EL_IN,dg_here%pa)
               QY_AVG = 0.5*(dg_here%QY_IN + dg_here%QY_EX)
#ifdef TRACE
               iota_AVG = dg_here%iota_IN
#endif

#ifdef SED_LAY
               do ll=1,s%layers
                  bed_AVG(ll) = 0.5D0*( dg_here%bed_IN(ll) + dg_here%bed_EX(ll) )
               enddo
#endif

!.....Compute the global_here%edge integral

            DO K = 1,dg_here%DOFS(global_here%EL_IN)
               CALL EDGE_INT_LDG_HYDRO&
              (s,dg_here,global_here,K,global_here%EL_IN,LED,GED,I,iota_AVG,ZE_AVG,QX_AVG,QY_AVG,dg_here%NX,dg_here%NY,dg_here%pa)
#ifdef SED_LAY
               do ll=1,s%layers
                  CALL EDGE_INT_LDG_sediment&
                 (K,global_here%EL_IN,LED,GED,I,bed_AVG(ll),dg_here%NX,dg_here%NY,dg_here%pa,ll)
               enddo
#endif
            ENDDO
         ENDDO

 1000 CONTINUE
      
      RETURN
      END SUBROUTINE


!***********************************************************************
!     
!     SUBROUTINE RADIATION_EDGE_LDG_HYDRO( )
!     
!     This subroutine does the following:
!     
!     1.  Calculates the values of the necessary variables at the global_here%edge
!     gauss points for INTERNAL edges (sediment diffusion?)
!     2.  Take the average of these values
!     3.  Perform boundary integration
!     
!     Written by Shintaro Bunya (01-10-2007)
!     01-10-2011 - cem - adapted for p_enrichment and multicomponent
!     
!***********************************************************************

      SUBROUTINE RADIATION_EDGE_LDG_HYDRO(s,dg_here,global_here)

!.....Use appropriate modules

      USE SIZES
      USE GLOBAL
      USE DG

      IMPLICIT NONE

      type (sizes_type) :: s
      type (dg_type) :: dg_here
      type (global_type) :: global_here

!.....Declare local variables

      INTEGER L, LED, GED,k,i,ll
      REAL(SZ) TX, TY

      DO 1000 L = 1,dg_here%NREDS
         
!.....Retrieve the global and local global_here%edge number

         GED = dg_here%NREDN(L)
         LED = dg_here%NEDSD(1,GED)

!.....Retrieve the elements which share the global_here%edge

         global_here%EL_IN = dg_here%NEDEL(1,GED)

         dg_here%PA = global_here%PDG_EL(global_here%EL_IN)

#ifdef P0
         if (dg_here%pa.eq.0) then
            dg_here%pa = 1
         endif
#endif

!.....Retrieve the components of the normal vector to the global_here%edge
         
         dg_here%NX = dg_here%COSNX(GED)
         dg_here%NY = dg_here%SINNX(GED)

!.....Compute dg_here%ZE, dg_here%QX, dg_here%QY, and dg_here%HB at each Gauss point

         DO I = 1,dg_here%NEGP(dg_here%pa)
            
#ifdef WAVE_DIF
            dg_here%ZE_IN = 0.D0
#endif
            dg_here%QX_IN = 0.D0
            dg_here%QY_IN = 0.D0
#ifdef TRACE
            dg_here%iota_IN = 0.d0
#endif
#ifdef SED_LAY
            do ll=1,s%layers
               dg_here%bed_IN(ll) = 0.D0
            enddo
#endif

!.....Compute the solution at the interior state

            DO K = 1,dg_here%DOFS(global_here%EL_IN)
               
#ifdef WAVE_DIF
               dg_here%ZE_IN = dg_here%ZE_IN + dg_here%ZE(K,global_here%EL_IN,dg_here%IRK)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
#endif
               dg_here%QX_IN = dg_here%QX_IN + dg_here%QX(K,global_here%EL_IN,dg_here%IRK)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
               dg_here%QY_IN = dg_here%QY_IN + dg_here%QY(K,global_here%EL_IN,dg_here%IRK)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
#ifdef TRACE
               dg_here%iota_IN = dg_here%iota_IN + dg_here%iota(K,global_here%EL_IN,dg_here%IRK)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
#endif

#ifdef SED_LAY
               do ll=1,s%layers
                  dg_here%bed_in(ll) = dg_here%bed_in(ll) + dg_here%bed(K,global_here%EL_IN,dg_here%IRK,ll)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
               enddo
#endif
        

            ENDDO

            dg_here%QX_IN = dg_here%QX_IN*dg_here%SFACED(I,LED,global_here%EL_IN,dg_here%pa)

!.....Compute the global_here%edge integral
           
            DO K = 1,dg_here%DOFS(global_here%EL_IN)
               CALL EDGE_INT_LDG_HYDRO&
              (s,dg_here,global_here,K,global_here%EL_IN,LED,GED,I,dg_here%iota_IN,dg_here%ZE_IN,dg_here%QX_IN,dg_here%QY_IN,dg_here%NX,dg_here%NY,dg_here%pa)
#ifdef SED_LAY
               do ll=1,s%layers
                  CALL EDGE_INT_LDG_sediment&
                 (K,global_here%EL_IN,LED,GED,I,dg_here%bed_IN(ll),dg_here%NX,dg_here%NY,dg_here%pa,ll)
               enddo
#endif
            ENDDO

         ENDDO

 1000 CONTINUE      
      RETURN
      END SUBROUTINE

!***********************************************************************
!     
!     SUBROUTINE INTERNAL_EDGE_LDG_HYDRO( )
!     
!     This subroutine does the following:
!     
!     1.  Calculates the values of the necessary variables at the global_here%edge
!     gauss points for INTERNAL edges (sediment diffusion?)
!     2.  Take the average of these values
!     3.  Perform boundary integration
!     
!     Written by Shintaro Bunya (01-02-2007)
!     01-10-2011 - cem - adapted for p_enrichment and multicomponent
!     
!***********************************************************************

      SUBROUTINE INTERNAL_EDGE_LDG_HYDRO(s,dg_here,global_here)

!.....Use appropriate modules
      
      USE SIZES
      USE GLOBAL
      USE DG

      IMPLICIT NONE

      type (sizes_type) :: s
      type (dg_type) :: dg_here
      type (global_type) :: global_here

!.....Declare local variables

      INTEGER L, LED_IN, LED_EX, GED, GP_IN, GP_EX,k,i,ll
      REAL(SZ) ZE_AVG,QX_AVG,QY_AVG,bed_AVG(s%layers),W_IN,W_EX
      real(sz) iota_AVG


      DO 1000 L = 1,dg_here%NIEDS

!.....Retrieve the global and local global_here%edge number

         GED = dg_here%NIEDN(L)
         LED_IN = dg_here%NEDSD(1,GED)
         LED_EX = dg_here%NEDSD(2,GED)

!.....Retrieve the elements which share the global_here%edge

         global_here%EL_IN = dg_here%NEDEL(1,GED)
         global_here%EL_EX = dg_here%NEDEL(2,GED)

         dg_here%EL = global_here%EL_IN

         IF (dg_here%DOFS(global_here%EL_EX).LT.dg_here%DOFS(global_here%EL_IN)) then
            dg_here%EL = global_here%EL_EX
         endif

         dg_here%pa = global_here%PDG_EL(dg_here%EL)

#ifdef P0
         if (dg_here%pa.eq.0) then
            dg_here%pa = 1
         endif
#endif

!.....Retrieve the components of the normal vector to the global_here%edge

         dg_here%NX = dg_here%COSNX(GED)
         dg_here%NY = dg_here%SINNX(GED)

!.....Compute dg_here%ZE, dg_here%QX, dg_here%QY, and dg_here%HB at each global_here%edge Gauss quadrature point

         DO I = 1,dg_here%NEGP(dg_here%pa)

            GP_IN = I
            GP_EX = dg_here%NEGP(dg_here%pa) - I + 1
            
#ifdef WAVE_DIF
            dg_here%ZE_IN = 0.D0
#endif
            dg_here%QX_IN = 0.D0
            dg_here%QY_IN = 0.D0
#ifdef WAVE_DIF
            dg_here%ZE_EX = 0.D0
#endif
            dg_here%QX_EX = 0.D0
            dg_here%QY_EX = 0.D0
#ifdef TRACE
            dg_here%iota_IN = 0.D0
            dg_here%iota_EX = 0.D0
#endif

!deal with sediment
#ifdef SED_LAY
            do ll=1,s%layers
               dg_here%bed_IN(ll) = 0.D0
               dg_here%bed_EX(ll) = 0.D0
            enddo
#endif

            DO K = 1,dg_here%DOFS(dg_here%EL)

#ifdef WAVE_DIF
               dg_here%ZE_IN = dg_here%ZE_IN + dg_here%ZE(K,global_here%EL_IN,dg_here%IRK)*dg_here%PHI_EDGE(K,GP_IN,LED_IN,dg_here%pa)
#endif
               dg_here%QX_IN = dg_here%QX_IN + dg_here%QX(K,global_here%EL_IN,dg_here%IRK)*dg_here%PHI_EDGE(K,GP_IN,LED_IN,dg_here%pa)
               dg_here%QY_IN = dg_here%QY_IN + dg_here%QY(K,global_here%EL_IN,dg_here%IRK)*dg_here%PHI_EDGE(K,GP_IN,LED_IN,dg_here%pa)
#ifdef WAVE_DIF
               dg_here%ZE_EX = dg_here%ZE_EX + dg_here%ZE(K,global_here%EL_EX,dg_here%IRK)*dg_here%PHI_EDGE(K,GP_EX,LED_EX,dg_here%pa)
#endif
               dg_here%QX_EX = dg_here%QX_EX + dg_here%QX(K,global_here%EL_EX,dg_here%IRK)*dg_here%PHI_EDGE(K,GP_EX,LED_EX,dg_here%pa)
               dg_here%QY_EX = dg_here%QY_EX + dg_here%QY(K,global_here%EL_EX,dg_here%IRK)*dg_here%PHI_EDGE(K,GP_EX,LED_EX,dg_here%pa)
#ifdef TRACE
               dg_here%iota_IN = dg_here%iota_IN + dg_here%iota(K,global_here%EL_IN,dg_here%IRK)*dg_here%PHI_EDGE(K,GP_IN,LED_IN,dg_here%pa)
               dg_here%iota_EX = dg_here%iota_EX + dg_here%iota(K,global_here%EL_EX,dg_here%IRK)*dg_here%PHI_EDGE(K,GP_EX,LED_EX,dg_here%pa)
#endif

#ifdef SED_LAY
               do ll=1,s%layers
                  dg_here%bed_IN(ll) = dg_here%bed_IN(ll) &
                 + dg_here%bed(K,global_here%EL_IN,dg_here%IRK,ll)*dg_here%PHI_EDGE(K,GP_IN,LED_IN,dg_here%pa)
                  dg_here%bed_EX(ll) = dg_here%bed_EX(ll) &
                 + dg_here%bed(K,global_here%EL_EX,dg_here%IRK,ll)*dg_here%PHI_EDGE(K,GP_EX,LED_EX,dg_here%pa)
               enddo
#endif

            ENDDO

!.....Take the average

#ifdef WAVE_DIF
            ZE_AVG = 0.5D0*(dg_here%ZE_IN + dg_here%ZE_EX)
#endif

            QX_AVG = 0.5D0*(dg_here%QX_IN*dg_here%SFACED(GP_IN,LED_IN,global_here%EL_IN,dg_here%pa) &
           + dg_here%QX_EX*dg_here%SFACED(GP_EX,LED_EX,global_here%EL_EX,dg_here%pa))
            QY_AVG = 0.5D0*(dg_here%QY_IN + dg_here%QY_EX)

#ifdef TRACE
            iota_AVG = 0.5D0*(dg_here%iota_IN + dg_here%iota_EX)
#endif

#ifdef SED_LAY
            do ll=1,s%layers
               bed_AVG(ll) = 0.5*( dg_here%bed_IN(ll) + dg_here%bed_EX(ll) )
            enddo
#endif

!.....Compute the global_here%edge integral

            DO K = 1,dg_here%DOFS(dg_here%EL)

               W_IN = 2.0*dg_here%M_INV(K,dg_here%pa)/global_here%AREAS(global_here%EL_IN)*dg_here%XLEN(GED)*&
              dg_here%PHI_EDGE(K,GP_IN,LED_IN,dg_here%pa)*dg_here%WEGP(GP_IN,dg_here%pa)
               W_EX = 2.0*dg_here%M_INV(K,dg_here%pa)/global_here%AREAS(global_here%EL_EX)*dg_here%XLEN(GED)*&
              dg_here%PHI_EDGE(K,GP_EX,LED_EX,dg_here%pa)*dg_here%WEGP(GP_EX,dg_here%pa)
               
#ifdef WAVE_DIF
               dg_here%HZ(K,1,1,global_here%EL_IN) = dg_here%HZ(K,1,1,global_here%EL_IN) - ZE_AVG*dg_here%NX*W_IN
               dg_here%HZ(K,2,2,global_here%EL_IN) = dg_here%HZ(K,2,2,global_here%EL_IN) - ZE_AVG*dg_here%NY*W_IN
#endif
               dg_here%LZ(K,1,1,global_here%EL_IN) = dg_here%LZ(K,1,1,global_here%EL_IN) - QX_AVG*dg_here%NX*W_IN
               dg_here%LZ(K,1,2,global_here%EL_IN) = dg_here%LZ(K,1,2,global_here%EL_IN) - QX_AVG*dg_here%NY*W_IN
               dg_here%LZ(K,2,1,global_here%EL_IN) = dg_here%LZ(K,2,1,global_here%EL_IN) - QY_AVG*dg_here%NX*W_IN
               dg_here%LZ(K,2,2,global_here%EL_IN) = dg_here%LZ(K,2,2,global_here%EL_IN) - QY_AVG*dg_here%NY*W_IN
               
#ifdef WAVE_DIF
               dg_here%HZ(K,1,1,global_here%EL_EX) = dg_here%HZ(K,1,1,global_here%EL_EX) + ZE_AVG*dg_here%NX*W_EX
               dg_here%HZ(K,2,2,global_here%EL_EX) = dg_here%HZ(K,2,2,global_here%EL_EX) + ZE_AVG*dg_here%NY*W_EX
#endif
               dg_here%LZ(K,1,1,global_here%EL_EX) = dg_here%LZ(K,1,1,global_here%EL_EX) + QX_AVG*dg_here%NX*W_EX
               dg_here%LZ(K,1,2,global_here%EL_EX) = dg_here%LZ(K,1,2,global_here%EL_EX) + QX_AVG*dg_here%NY*W_EX
               dg_here%LZ(K,2,1,global_here%EL_EX) = dg_here%LZ(K,2,1,global_here%EL_EX) + QY_AVG*dg_here%NX*W_EX
               dg_here%LZ(K,2,2,global_here%EL_EX) = dg_here%LZ(K,2,2,global_here%EL_EX) + QY_AVG*dg_here%NY*W_EX
#ifdef TRACE
               dg_here%TZ(K,1,1,global_here%EL_IN) = dg_here%TZ(K,1,1,global_here%EL_IN) - iota_AVG*dg_here%NX*W_IN
               dg_here%TZ(K,2,2,global_here%EL_IN) = dg_here%TZ(K,2,2,global_here%EL_IN) - iota_AVG*dg_here%NY*W_IN

               dg_here%TZ(K,1,1,global_here%EL_EX) = dg_here%TZ(K,1,1,global_here%EL_EX) + iota_AVG*dg_here%NX*W_EX
               dg_here%TZ(K,2,2,global_here%EL_EX) = dg_here%TZ(K,2,2,global_here%EL_EX) + iota_AVG*dg_here%NY*W_EX
#endif

#ifdef SED_LAY
               do ll=1,s%layers
                  dg_here%MZ(K,1,ll,global_here%EL_IN) = dg_here%MZ(K,1,ll,global_here%EL_IN) &
                 - bed_AVG(ll)*dg_here%NX*W_IN
                  dg_here%MZ(K,2,ll,global_here%EL_IN) = dg_here%MZ(K,2,ll,global_here%EL_IN) &
                 - bed_AVG(ll)*dg_here%NY*W_IN
                  dg_here%MZ(K,1,ll,global_here%EL_EX) = dg_here%MZ(K,1,ll,global_here%EL_EX) &
                 + bed_AVG(ll)*dg_here%NX*W_EX
                  dg_here%MZ(K,2,ll,global_here%EL_EX) = dg_here%MZ(K,2,ll,global_here%EL_EX) &
                 + bed_AVG(ll)*dg_here%NY*W_EX
               enddo
#endif

            enddo
         enddo
         
 1000 CONTINUE
      RETURN
      END SUBROUTINE


!***********************************************************************
!     
!     SUBROUTINE EDGE_INT_LDG_HYDRO()
!     
!     This subroutine computes the edge integrals for the LDG boundary
!     terms using Gauss quadrature and adds them to the LZ. (sediment diffusion?)
!     
!     Written by Shintaro Bunya (01-02-2007)
!     01-10-2011 - cem - adapted for p_enrichment and multicomponent
!     
!***********************************************************************

      SUBROUTINE EDGE_INT_LDG_HYDRO(s,dg_here,global_here,K,EL,LED,GED,GP,iota_Avg,ZE_Avg,QX_Avg,&
     QY_Avg,NX,NY,pa)

!.....Use appropriate modules

      USE SIZES
      USE GLOBAL
      USE DG
 
      IMPLICIT NONE
      
      type (sizes_type) :: s
      type (dg_type) :: dg_here
      type (global_type) :: global_here
      
!.....Declare local variables

      INTEGER K,EL,LED,GED,GP,i,pa
      REAL(SZ) AREA, IMASS,n1,n2,n3,QX_Avg,QY_Avg
      Real(SZ) ZE_Avg,iota_Avg
      real(sz) NX,NY
      
!.....Retrieve the element area
      
      AREA = 0.5D0*global_here%AREAS(EL)

!.....Comput the edge integral

      IMASS = dg_here%M_INV(K,pa)/(0.5D0*AREA)

#ifdef WAVE_DIF
      dg_here%HZ(K,1,1,EL) = dg_here%HZ(K,1,1,EL)&
     - IMASS*dg_here%XLEN(GED)*0.5*ZE_Avg*dg_here%PHI_EDGE(K,GP,LED,pa)*NX*dg_here%WEGP(GP,pa)
      dg_here%HZ(K,2,2,EL) = dg_here%HZ(K,2,2,EL)&
     - IMASS*dg_here%XLEN(GED)*0.5*ZE_Avg*dg_here%PHI_EDGE(K,GP,LED,pa)*NY*dg_here%WEGP(GP,pa)
#endif
      dg_here%LZ(K,1,1,EL) = dg_here%LZ(K,1,1,EL)&
     - IMASS*dg_here%XLEN(GED)*0.5*QX_Avg*dg_here%PHI_EDGE(K,GP,LED,pa)*NX*dg_here%WEGP(GP,pa)
      dg_here%LZ(K,1,2,EL) = dg_here%LZ(K,1,2,EL)&
     - IMASS*dg_here%XLEN(GED)*0.5*QX_Avg*dg_here%PHI_EDGE(K,GP,LED,pa)*NY*dg_here%WEGP(GP,pa)
      dg_here%LZ(K,2,1,EL) = dg_here%LZ(K,2,1,EL)&
     - IMASS*dg_here%XLEN(GED)*0.5*QY_Avg*dg_here%PHI_EDGE(K,GP,LED,pa)*NX*dg_here%WEGP(GP,pa)
      dg_here%LZ(K,2,2,EL) = dg_here%LZ(K,2,2,EL)&
     - IMASS*dg_here%XLEN(GED)*0.5*QY_Avg*dg_here%PHI_EDGE(K,GP,LED,pa)*NY*dg_here%WEGP(GP,pa)

#ifdef TRACE
      dg_here%TZ(K,1,1,EL) = dg_here%TZ(K,1,1,EL)&
     - IMASS*dg_here%XLEN(GED)*0.5D0*iota_Avg*dg_here%PHI_EDGE(K,GP,LED,pa)*NX*dg_here%WEGP(GP,pa)
      TZ(K,2,2,EL) = TZ(K,2,2,EL)&
     - IMASS*XLEN(GED)*0.5D0*iota_Avg*PHI_EDGE(K,GP,LED,pa)*NY*WEGP(GP,pa)

#endif
      
      RETURN
      END SUBROUTINE

!***********************************************************************
!     
!     SUBROUTINE EDGE_INT_LDG_sediment()
!     
!     This subroutine computes the edge integrals for the LDG boundary
!     terms for sediment 
!     
!     2012 - cem
!     
!***********************************************************************

      SUBROUTINE EDGE_INT_LDG_sediment&
     (s,dg_here,global_here,K,EL,LED,GED,GP,bed_avg,NX,NY,pa,ll)
                                ! <ezpp-noinst>
      
!.....Use appropriate modules

      USE SIZES
      USE GLOBAL
      USE DG

      IMPLICIT NONE
      
      type (sizes_type) :: s
      type (dg_type) :: dg_here
      type (global_type) :: global_here

!.....Declare local variables

      INTEGER K, EL, LED, GED, GP,i,pa,ll
      REAL(SZ) AREA, IMASS, bed_avg, NX, NY
      
!.....Retrieve the element area
      
      AREA = 0.5D0*global_here%AREAS(EL)

!.....Comput the edge integral

      IMASS = dg_here%M_INV(K,pa)/(0.5D0*AREA)

      DG_HERE%MZ(K,1,ll,EL) = DG_HERE%MZ(K,1,ll,EL)&
     - IMASS*dg_here%XLEN(GED)*0.5*bed_avg*dg_here%PHI_EDGE(K,GP,LED,pa)*NX*dg_here%WEGP(GP,pa)
      DG_HERE%MZ(K,2,ll,EL) = DG_HERE%MZ(K,2,ll,EL)&
     - IMASS*dg_here%XLEN(GED)*0.5*bed_avg*dg_here%PHI_EDGE(K,GP,LED,pa)*NY*dg_here%WEGP(GP,pa)
      
      RETURN
      END SUBROUTINE

!***********************************************************************
!     
!     SUBROUTINE RHS_LDG_HYDRO()
!     
!     This subroutine computes the area integrals for the LDG hydro and
!     adds them into the dg_here%LZ (+sediment diffusion)
!     
!     Written by Shintaro Bunya (01-02-2007)
!     01-10-2011 - cem - adapted for p_enrichment and multicomponent
!     2012 - cem - added sediment layers
!***********************************************************************

      SUBROUTINE RHS_LDG_HYDRO(s,dg_here,global_here)
      
!.....Use appropriate modules

      USE SIZES
      USE GLOBAL
      USE DG

      IMPLICIT NONE

      type (sizes_type) :: s
      type (dg_type) :: dg_here
      type (global_type) :: global_here
      
!.....Declare local variables
      INTEGER L,K,I,pa,ll,kk
      real(sz) ze_sensor1, qx_sensor1, qy_sensor1
      real(sz) ze_sensor2, qx_sensor2, qy_sensor2
      real(sz) iota_sensor1,iota_sensor2
      logical fl1,fl2,fl3,fl4,fl5
#ifdef SED_LAY
      real(sz) tbed_sensor1,tbed_sensor2
#endif

      DO L=1,global_here%NE

         pa = global_here%PDG_EL(L)

#ifdef P0
         if (pa.eq.0) then
            pa = 1
         endif
#endif
         
!.....Retrieve the global node numbers for the element
         
         global_here%N1 = global_here%NM(L,1)
         global_here%N2 = global_here%NM(L,2)
         global_here%N3 = global_here%NM(L,3)
         
!.....Compute dg_here%ZE, dg_here%QX, dg_here%QY, and HB at each area Gauss quadrature point

         DO I = 1,dg_here%NAGP(pa)
            
#ifdef WAVE_DIF
            dg_here%ZE_IN = 0.D0
#endif
            dg_here%QX_IN = 0.D0
            dg_here%QY_IN = 0.D0
#ifdef TRACE
            dg_here%iota_IN = 0.D0
#endif

#ifdef SED_LAY
            do ll=1,s%layers
               dg_here%bed_IN(ll) = 0.D0
            enddo
#endif
            
            DO K = 1,dg_here%DOFS(L)
#ifdef WAVE_DIF
               dg_here%ZE_IN = dg_here%ZE_IN + dg_here%ZE(K,L,dg_here%IRK)*dg_here%PHI_AREA(K,I,pa)
#endif
               dg_here%QX_IN = dg_here%QX_IN + dg_here%QX(K,L,dg_here%IRK)*dg_here%PHI_AREA(K,I,pa)
               dg_here%QY_IN = dg_here%QY_IN + dg_here%QY(K,L,dg_here%IRK)*dg_here%PHI_AREA(K,I,pa)
#ifdef TRACE
               dg_here%iota_IN = dg_here%iota_IN + dg_here%iota(K,L,dg_here%IRK)*dg_here%PHI_AREA(K,I,pa)
#endif

#ifdef SED_LAY
               do ll=1,s%layers
                  dg_here%bed_IN(ll) = dg_here%bed_IN(ll) + dg_here%bed(K,L,dg_here%IRK,ll)*dg_here%PHI_AREA(K,I,pa)
               enddo
#endif

            ENDDO

!.....Build the rhs
                
            DO K = 1,dg_here%DOFS(L)
#ifdef WAVE_DIF
               dg_here%HZ(K,1,1,L) = dg_here%HZ(K,1,1,L)&
                    + dg_here%M_INV(K,pa)*dg_here%ZE_IN*dg_here%SFAC_ELEM(I,L,pa)*& ! <--- dg_here%ZE/Mk
                    ( dg_here%DRPHI(K,I,pa)*dg_here%DRDX(L) + dg_here%DSPHI(K,I,pa)*dg_here%DSDX(L) )*& ! <--- dphi/dx
                    dg_here%WAGP(I,pa)  ! <--- weight
               dg_here%HZ(K,2,2,L) = dg_here%HZ(K,2,2,L)&
                    + dg_here%M_INV(K,pa)*dg_here%ZE_IN*& ! <--- dg_here%ZE/Mk
                    ( dg_here%DRPHI(K,I,pa)*dg_here%DRDY(L) + dg_here%DSPHI(K,I,pa)*dg_here%DSDY(L) )*& ! <--- dphi/dy
                    dg_here%WAGP(I,pa)  ! <--- weight
#endif

               dg_here%LZ(K,1,1,L) = dg_here%LZ(K,1,1,L)&
                   + dg_here%M_INV(K,pa)*dg_here%QX_IN*dg_here%SFAC_ELEM(I,L,pa)*& ! <--- dg_here%QX/Mk
                   ( dg_here%DRPHI(K,I,pa)*dg_here%DRDX(L) + dg_here%DSPHI(K,I,pa)*dg_here%DSDX(L) )*& ! <--- dphi/dx
                   dg_here%WAGP(I,pa)  ! <--- weight
               
               dg_here%LZ(K,1,2,L) = dg_here%LZ(K,1,2,L)&
                    + dg_here%M_INV(K,pa)*dg_here%QX_IN*& ! <--- dg_here%QX/Mk
                    ( dg_here%DRPHI(K,I,pa)*dg_here%DRDY(L) + dg_here%DSPHI(K,I,pa)*dg_here%DSDY(L) )*& ! <--- dphi/dy
                    dg_here%WAGP(I,pa)  ! <--- weight

               dg_here%LZ(K,2,1,L) = dg_here%LZ(K,2,1,L)&
                    + dg_here%M_INV(K,pa)*dg_here%QY_IN*dg_here%SFAC_ELEM(I,L,pa)*& ! <--- dg_here%QY/Mk
                    ( dg_here%DRPHI(K,I,pa)*dg_here%DRDX(L) + dg_here%DSPHI(K,I,pa)*dg_here%DSDX(L) )*& ! <--- dphi/dx
                    dg_here%WAGP(I,pa)  ! <--- weight

               dg_here%LZ(K,2,2,L) = dg_here%LZ(K,2,2,L)&
                    + dg_here%M_INV(K,pa)*dg_here%QY_IN*& ! <--- dg_here%QY/Mk
                    ( dg_here%DRPHI(K,I,pa)*dg_here%DRDY(L) + dg_here%DSPHI(K,I,pa)*dg_here%DSDY(L) )*& ! <--- dphi/dy
                    dg_here%WAGP(I,pa)  ! <--- weight
#ifdef TRACE
               dg_here%TZ(K,1,1,L) = dg_here%TZ(K,1,1,L)&
                    + dg_here%M_INV(K,pa)*dg_here%iota_IN*dg_here%SFAC_ELEM(I,L,pa)*& ! <--- dg_here%iota*H/Mk
                    ( dg_here%DRPHI(K,I,pa)*dg_here%DRDX(L) + dg_here%DSPHI(K,I,pa)*dg_here%DSDX(L) )*& ! <--- dphi/dx
                    dg_here%WAGP(I,pa)  ! <--- weight
               dg_here%TZ(K,2,2,L) = dg_here%TZ(K,2,2,L)&
                    + dg_here%M_INV(K,pa)*dg_here%iota_IN*& ! <--- dg_here%iota*H/Mk
                    ( dg_here%DRPHI(K,I,pa)*dg_here%DRDY(L) + dg_here%DSPHI(K,I,pa)*dg_here%DSDY(L) )*& ! <--- dphi/dy
                    dg_here%WAGP(I,pa)  ! <--- weight
#endif

#ifdef SED_LAY
               do ll=1,s%layers
                  dg_here%MZ(K,1,ll,L) =  dg_here%MZ(K,1,ll,L)&
                 + dg_here%M_INV(K,pa)*dg_here%bed_IN(ll)*dg_here%SFAC_ELEM(I,L,pa)* &
                 ( dg_here%DRPHI(K,I,pa)*dg_here%DRDX(L) + dg_here%DSPHI(K,I,pa)*dg_here%DSDX(L) )* &
                 dg_here%WAGP(I,pa) 
                  dg_here%MZ(K,2,ll,L) =  dg_here%MZ(K,2,ll,L)&
                 + dg_here%M_INV(K,pa)*dg_here%bed_IN(ll)* &
                 ( dg_here%DRPHI(K,I,pa)*dg_here%DRDY(L) + dg_here%DSPHI(K,I,pa)*dg_here%DSDY(L) )* &
                 dg_here%WAGP(I,pa) 
               enddo
#endif

            enddo
         enddo

#ifdef ARTDIF

         dg_here%slimit1 = -100.D0
         dg_here%slimit2 = -100.D0
         dg_here%slimit3 = -100.D0

         global_here%entrop(:,L) = -100.D0
           
#ifdef WAVE_DIF 
         ze_sensor1 = 0.d0
         ze_sensor2 = 0.d0  
#endif             
         qx_sensor1 = 0.d0
         qy_sensor1 = 0.d0
           
         qx_sensor2 = 0.d0
         qy_sensor2 = 0.d0

#ifdef TRACE
         dg_here%slimit4 = -100.D0 
         iota_sensor1 = 0.d0
         iota_sensor2 = 0.d0
#endif
         
#ifdef SED_LAY
         dg_here%slimit5 = -100.D0
         tbed_sensor1 = 0.d0
         tbed_sensor2 = 0.d0
#endif
         
         do I = 1,dg_here%NAGP(pa)
            
            do kk = 2, dg_here%dofs(L)  !Compute the first sensor
#ifdef WAVE_DIF 
               ze_sensor1 = ze_sensor1 + (dg_here%ze(kk,L,dg_here%irk)* &
              dg_here%phi_area(kk,I,pa))**2.D0 * dg_here%wagp(I,pa)
#endif
               qx_sensor1 = qx_sensor1 + (dg_here%qx(kk,L,dg_here%irk)* &
              dg_here%phi_area(kk,I,pa))**2.D0 * dg_here%wagp(I,pa)
               qy_sensor1 = qy_sensor1 + (dg_here%qy(kk,L,dg_here%irk)* &
              dg_here%phi_area(kk,I,pa))**2.D0 * dg_here%wagp(I,pa)
#ifdef TRACE
               iota_sensor1 = iota_sensor1 + (dg_here%iota(kk,L,dg_here%irk)&
              * dg_here%phi_area(kk,I,pa))**2.D0 * dg_here%wagp(I,pa)
#endif
               
#ifdef SED_LAY       
               do ll=1,s%layers
                  tbed_sensor1 = tbed_sensor1 + (dg_here%bed(kk,L,dg_here%irk,ll)* &
                 dg_here%phi_area(kk,I,pa))**2.D0 * dg_here%wagp(I,pa)
               enddo
#endif

            enddo
            
            do kk = 1,dg_here%dofs(L) !Compute the second sensor
#ifdef WAVE_DIF 
               ze_sensor2 = ze_sensor2 + (dg_here%ze(kk,L,dg_here%irk)* &
              dg_here%phi_area(kk,I,pa))**2.D0 * dg_here%wagp(I,pa)
#endif
               qx_sensor2 = qx_sensor2 + (dg_here%qx(kk,L,dg_here%irk)* &
              dg_here%phi_area(kk,I,pa))**2.D0 * dg_here%wagp(I,pa)
               qy_sensor2 = qy_sensor2 + (dg_here%qy(kk,L,dg_here%irk)* &
              dg_here%phi_area(kk,I,pa))**2.D0 * dg_here%wagp(I,pa)
#ifdef TRACE
               iota_sensor2 = iota_sensor2 + (dg_here%iota(kk,L,dg_here%irk)&
              * dg_here%phi_area(kk,I,pa))**2.D0 * dg_here%wagp(I,pa)
#endif
               
#ifdef SED_LAY       
               do ll=1,s%layers
                  tbed_sensor2 = tbed_sensor2 + (dg_here%bed(kk,L,dg_here%irk,ll)*& !Adjust sensor for multilayers!&
                 dg_here%phi_area(kk,I,pa))**2.D0 * dg_here%wagp(I,pa)
               enddo
#endif

            enddo
            
         enddo
#ifdef WAVE_DIF 
         if  (ze_sensor2.gt.1.0e-4.and.ze_sensor1.gt.1.0e-4 ) then
            dg_here%slimit1 = log10( ze_sensor1/ze_sensor2 ) + dg_here%balance(1)
            global_here%entrop(1,L) = dg_here%slimit1
         endif
#endif
         if (qx_sensor2.gt.1.0e-4.and.qx_sensor1.gt.1.0e-4 ) then
            dg_here%slimit2 = log10( qx_sensor1/qx_sensor2 ) + dg_here%balance(2)
            global_here%entrop(2,L) = dg_here%slimit2
         endif
         if ( qy_sensor2.gt.1.0e-4.and.qy_sensor1.gt.1.0e-4 ) then
            dg_here%slimit3 = log10( qy_sensor1/qy_sensor2 ) + dg_here%balance(3)
            global_here%entrop(3,L) = dg_here%slimit3
         endif
#ifdef TRACE
          if ( iota_sensor2.gt.1.0e-4.and.iota_sensor1.gt.1.0e-4 ) then
            dg_here%slimit4 = log10( iota_sensor1/iota_sensor2 ) + dg_here%balance(4)
            global_here%entrop(4,L) = dg_here%slimit4
          endif
#endif
         
#ifdef SED_LAY
         if ( tbed_sensor2.gt.1.0e-4.and.tbed_sensor1.gt.1.0e-4 ) then
            dg_here%slimit5 = log10( tbed_sensor1/tbed_sensor2 ) + dg_here%balance(5)
            global_here%entrop(5,L) = dg_here%slimit5 
        endif
#endif
                  
#endif                  
 
      ENDDO

      RETURN
      END SUBROUTINE

