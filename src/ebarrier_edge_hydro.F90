      SUBROUTINE EBARRIER_EDGE_HYDRO(s,dg_here,global_here,IT)
      
!.....Use appropriate modules

      use sizes
      USE GLOBAL
      USE DG

      IMPLICIT NONE

      type (sizes_type) :: s
      type (dg_type) :: dg_here
      type (global_type) :: global_here

!.....Declare local variables

      INTEGER GED, LED, L,i,k,ll,IT
      REAL(SZ) ABOVE, SUBSUP, W_IN, TX, TY

!.....Loop over the external barrier segments

      DO 1000 L = 1,dg_here%NEBSEG
         
!.......Obtain the global and local global_here%edge numbers

         GED = dg_here%NEBSEGN(L)
         LED = dg_here%NEDSD(1,GED)

!.......Obtain the element

         dg_here%EL = dg_here%NEDEL(1,GED)
         dg_here%PA = global_here%PDG_EL(global_here%EL_IN)

#ifdef P0
         if (dg_here%pa.eq.0) then
            dg_here%pa = 1
         endif
#endif
         
         IF (dg_here%WDFLG(dg_here%EL).EQ.0) GOTO 1000

!.....Retrieve the components of the normal vector to the global_here%edge

         dg_here%NX = dg_here%COSNX(GED)
         dg_here%NY = dg_here%SINNX(GED)
         
!.....Set the components for the tangential vector to the global_here%edge

         TX = -dg_here%NY
         TY =  dg_here%NX
         
!.......Compute the variables at the quadrature points

         DO I = 1,dg_here%NEGP(dg_here%pa)
            
!.........Obtain the height of the barrier at the quadrature point

            dg_here%ZE_IN = dg_here%ZE(1,dg_here%EL,dg_here%IRK)
            dg_here%QX_IN = dg_here%QX(1,dg_here%EL,dg_here%IRK)
            dg_here%QY_IN = dg_here%QY(1,dg_here%EL,dg_here%IRK)
            dg_here%HB_IN = dg_here%BATHED(I,LED,global_here%EL_IN,dg_here%pa)

#ifdef TRACE
            dg_here%iota_IN = dg_here%iota(1,dg_here%EL,dg_here%IRK)
#endif

#ifdef CHEM
            dg_here%iota_IN = dg_here%iota(1,dg_here%EL,dg_here%IRK)
            dg_here%iota2_IN = dg_here%iota2(1,dg_here%EL,dg_here%IRK)
#endif

#ifdef DYNP
            dynP_IN = dg_here%dynP(1,dg_here%EL,dg_here%IRK)
#endif

#ifdef SED_LAY
            dg_here%HB(:,dg_here%EL,dg_here%irk) = 0.D0
            do ll=1,layers
               dg_here%HB(1,dg_here%EL,dg_here%irk) = dg_here%HB(1,dg_here%EL,dg_here%irk) + dg_here%bed(1,dg_here%EL,dg_here%irk,ll)
            enddo
            dg_here%bed_IN(:) = dg_here%bed(1,dg_here%EL,dg_here%irk,:)
            dg_here%HB_IN = dg_here%HB(1,global_here%EL_IN,dg_here%irk)
#endif


            DO K = 2,dg_here%DOFS(dg_here%EL)
               dg_here%ZE_IN = dg_here%ZE_IN + dg_here%ZE(K,dg_here%EL,dg_here%IRK)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
               dg_here%QX_IN = dg_here%QX_IN + dg_here%QX(K,dg_here%EL,dg_here%IRK)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
               dg_here%QY_IN = dg_here%QY_IN + dg_here%QY(K,dg_here%EL,dg_here%IRK)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
               !dg_here%HB_IN = dg_here%HB_IN + dg_here%HB(K,dg_here%EL,dg_here%IRK)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)

#ifdef TRACE
               dg_here%iota_IN = dg_here%iota_IN + dg_here%iota(K,dg_here%EL,dg_here%IRK)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
#endif

#ifdef CHEM
               dg_here%iota_IN = dg_here%iota_IN + dg_here%iota(K,dg_here%EL,dg_here%IRK)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
               dg_here%iota2_IN = dg_here%iota2_IN + dg_here%iota2(K,dg_here%EL,dg_here%IRK)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
#endif

#ifdef DYNP
               dynP_IN = dynP_IN + dg_here%dynP(K,dg_here%EL,dg_here%IRK)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
#endif

#ifdef SED_LAY

            do ll = 1,layers
               
               dg_here%bed_IN(ll) = dg_here%bed_IN(ll) + dg_here%bed(K,dg_here%EL,dg_here%IRK,ll)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
               dg_here%HB_IN = dg_here%HB_IN + dg_here%bed(k,dg_here%EL,dg_here%irk,ll)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
               
            enddo
#endif

            ENDDO

            dg_here%SFAC_IN = dg_here%SFACED(I,LED,dg_here%EL,dg_here%pa)

            ABOVE  = dg_here%ZE_IN - dg_here%EBHT(L)
            SUBSUP = 2.D0*ABOVE/3.D0
            
!.........Case 1:  Water is below barrier
!     ---------------------------------------------
            
            IF (ABOVE.LT.BARMIN) THEN
               global_here%Q_N_EXT = -(dg_here%QX_IN*dg_here%NX + dg_here%QY_IN*dg_here%NY)
               global_here%Q_T_EXT =   dg_here%QX_IN*TX + dg_here%QY_IN*TY
               
!.........Case 2:  Water is above barrier
!     ------------------------------------------------------------
               
            ELSE
               global_here%Q_N_EXT = dg_here%RAMPDG*dg_here%EBCFSP(L)*SUBSUP*SQRT(SUBSUP*global_here%G)
               global_here%Q_T_EXT = 0.D0
            ENDIF
            
!.....Set the exterior dg_here%bed and surface elevation equal to the interior

            dg_here%ZE_EX = dg_here%ZE_IN
            dg_here%HB_EX = dg_here%HB_IN
            dg_here%SFAC_EX = dg_here%SFAC_IN

#ifdef TRACE
            dg_here%iota_EX = dg_here%iota_IN
#endif

#ifdef CHEM
            dg_here%iota_EX = dg_here%iota_IN
            dg_here%iota2_EX = dg_here%iota2_IN
#endif

#ifdef DYNP
            dynP_EX = dynP_IN
#endif

#ifdef SED_LAY

            dg_here%bed_EX(:) = dg_here%bed_IN(:)

#endif

!.........Compute the global_here%x and global_here%y components of the external state flow

            dg_here%QX_EX = ( TY*global_here%Q_N_EXT - dg_here%NY*global_here%Q_T_EXT)/(dg_here%NX*TY - dg_here%NY*TX)
            dg_here%QY_EX = (-TX*global_here%Q_N_EXT + dg_here%NX*global_here%Q_T_EXT)/(dg_here%NX*TY - dg_here%NY*TX)
            
!.........Compute numerical flux
            
            CALL NUMERICAL_FLUX(s,dg_here,global_here,IT)

            DO K = 1,dg_here%DOFS(dg_here%EL)

               W_IN = 2.0*dg_here%M_INV(K,dg_here%pa)/global_here%AREAS(dg_here%EL)*dg_here%XLEN(GED)&
              *dg_here%PHI_EDGE(K,I,LED,dg_here%pa)*dg_here%WEGP(I,dg_here%pa)

               dg_here%RHS_ZE(K,dg_here%EL,dg_here%IRK) = dg_here%RHS_ZE(K,dg_here%EL,dg_here%IRK) - W_IN*global_here%F_HAT
               dg_here%RHS_QX(K,dg_here%EL,dg_here%IRK) = dg_here%RHS_QX(K,dg_here%EL,dg_here%IRK) - W_IN*global_here%G_HAT
               dg_here%RHS_QY(K,dg_here%EL,dg_here%IRK) = dg_here%RHS_QY(K,dg_here%EL,dg_here%IRK) - W_IN*global_here%H_HAT

#ifdef TRACE
               dg_here%RHS_iota(K,dg_here%EL,dg_here%IRK) = dg_here%RHS_iota(K,dg_here%EL,dg_here%IRK) - W_IN*global_here%I_HAT
#endif

#ifdef CHEM
               dg_here%RHS_iota(K,dg_here%EL,dg_here%IRK) = dg_here%RHS_iota(K,dg_here%EL,dg_here%IRK) - W_IN*global_here%I_HAT
               dg_here%RHS_iota2(K,dg_here%EL,dg_here%IRK) = dg_here%RHS_iota2(K,dg_here%EL,dg_here%IRK) - W_IN*global_here%J_HAT
#endif

#ifdef DYNP
               dg_here%RHS_dynP(K,dg_here%EL,dg_here%IRK) = dg_here%RHS_dynP(K,dg_here%EL,dg_here%IRK) - W_IN*global_here%I_HAT
#endif

#ifdef SED_LAY

               do ll=1,layers

                  dg_here%RHS_bed(K,dg_here%EL,dg_here%IRK,ll) = dg_here%RHS_bed(K,dg_here%EL,dg_here%IRK,ll) - W_IN*dg_here%Bed_HAT(ll)

               enddo

#endif

            ENDDO
         ENDDO
 1000 CONTINUE

      RETURN
      END SUBROUTINE EBARRIER_EDGE_HYDRO
