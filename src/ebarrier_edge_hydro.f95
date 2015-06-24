      SUBROUTINE EBARRIER_EDGE_HYDRO(s,dg,IT)
      
!.....Use appropriate modules

      use sizes
      USE GLOBAL
      USE DG

      IMPLICIT NONE

      type (sizes_type) :: s
      type (dg_type) :: dg

!.....Declare local variables

      INTEGER GED, LED, L,i,k,ll,IT
      REAL(SZ) ABOVE, SUBSUP, W_IN, TX, TY

!.....Loop over the external barrier segments

      DO 1000 L = 1,dg%NEBSEG
         
!.......Obtain the global and local edge numbers

         GED = dg%NEBSEGN(L)
         LED = dg%NEDSD(1,GED)

!.......Obtain the element

         dg%EL = dg%NEDEL(1,GED)
         dg%PA = PDG_EL(EL_IN)

#ifdef P0
         if (dg%pa.eq.0) then
            dg%pa = 1
         endif
#endif
         
         IF (dg%WDFLG(dg%EL).EQ.0) GOTO 1000

!.....Retrieve the components of the normal vector to the edge

         dg%NX = dg%COSNX(GED)
         dg%NY = dg%SINNX(GED)
         
!.....Set the components for the tangential vector to the edge

         TX = -dg%NY
         TY =  dg%NX
         
!.......Compute the variables at the quadrature points

         DO I = 1,dg%NEGP(dg%pa)
            
!.........Obtain the height of the barrier at the quadrature point

            dg%ZE_IN = dg%ZE(1,dg%EL,dg%IRK)
            dg%QX_IN = dg%QX(1,dg%EL,dg%IRK)
            dg%QY_IN = dg%QY(1,dg%EL,dg%IRK)
            dg%HB_IN = dg%BATHED(I,LED,EL_IN,dg%pa)

#ifdef TRACE
            dg%iota_IN = dg%iota(1,dg%EL,dg%IRK)
#endif

#ifdef CHEM
            dg%iota_IN = dg%iota(1,dg%EL,dg%IRK)
            dg%iota2_IN = dg%iota2(1,dg%EL,dg%IRK)
#endif

#ifdef DYNP
            dynP_IN = dg%dynP(1,dg%EL,dg%IRK)
#endif

#ifdef SED_LAY
            dg%HB(:,dg%EL,dg%irk) = 0.D0
            do ll=1,layers
               dg%HB(1,dg%EL,dg%irk) = dg%HB(1,dg%EL,dg%irk) + dg%bed(1,dg%EL,dg%irk,ll)
            enddo
            dg%bed_IN(:) = dg%bed(1,dg%EL,dg%irk,:)
            dg%HB_IN = dg%HB(1,EL_IN,dg%irk)
#endif


            DO K = 2,dg%DOFS(dg%EL)
               dg%ZE_IN = dg%ZE_IN + dg%ZE(K,dg%EL,dg%IRK)*dg%PHI_EDGE(K,I,LED,dg%pa)
               dg%QX_IN = dg%QX_IN + dg%QX(K,dg%EL,dg%IRK)*dg%PHI_EDGE(K,I,LED,dg%pa)
               dg%QY_IN = dg%QY_IN + dg%QY(K,dg%EL,dg%IRK)*dg%PHI_EDGE(K,I,LED,dg%pa)
               !dg%HB_IN = dg%HB_IN + dg%HB(K,dg%EL,dg%IRK)*dg%PHI_EDGE(K,I,LED,dg%pa)

#ifdef TRACE
               dg%iota_IN = dg%iota_IN + dg%iota(K,dg%EL,dg%IRK)*dg%PHI_EDGE(K,I,LED,dg%pa)
#endif

#ifdef CHEM
               dg%iota_IN = dg%iota_IN + dg%iota(K,dg%EL,dg%IRK)*dg%PHI_EDGE(K,I,LED,dg%pa)
               dg%iota2_IN = dg%iota2_IN + dg%iota2(K,dg%EL,dg%IRK)*dg%PHI_EDGE(K,I,LED,dg%pa)
#endif

#ifdef DYNP
               dynP_IN = dynP_IN + dg%dynP(K,dg%EL,dg%IRK)*dg%PHI_EDGE(K,I,LED,dg%pa)
#endif

#ifdef SED_LAY

            do ll = 1,layers
               
               dg%bed_IN(ll) = dg%bed_IN(ll) + dg%bed(K,dg%EL,dg%IRK,ll)*dg%PHI_EDGE(K,I,LED,dg%pa)
               dg%HB_IN = dg%HB_IN + dg%bed(k,dg%EL,dg%irk,ll)*dg%PHI_EDGE(K,I,LED,dg%pa)
               
            enddo
#endif

            ENDDO

            dg%SFAC_IN = dg%SFACED(I,LED,dg%EL,dg%pa)

            ABOVE  = dg%ZE_IN - dg%EBHT(L)
            SUBSUP = 2.D0*ABOVE/3.D0
            
!.........Case 1:  Water is below barrier
!     ---------------------------------------------
            
            IF (ABOVE.LT.BARMIN) THEN
               Q_N_EXT = -(dg%QX_IN*dg%NX + dg%QY_IN*dg%NY)
               Q_T_EXT =   dg%QX_IN*TX + dg%QY_IN*TY
               
!.........Case 2:  Water is above barrier
!     ------------------------------------------------------------
               
            ELSE
               Q_N_EXT = dg%RAMPDG*dg%EBCFSP(L)*SUBSUP*SQRT(SUBSUP*G)
               Q_T_EXT = 0.D0
            ENDIF
            
!.....Set the exterior dg%bed and surface elevation equal to the interior

            dg%ZE_EX = dg%ZE_IN
            dg%HB_EX = dg%HB_IN
            dg%SFAC_EX = dg%SFAC_IN

#ifdef TRACE
            dg%iota_EX = dg%iota_IN
#endif

#ifdef CHEM
            dg%iota_EX = dg%iota_IN
            dg%iota2_EX = dg%iota2_IN
#endif

#ifdef DYNP
            dynP_EX = dynP_IN
#endif

#ifdef SED_LAY

            dg%bed_EX(:) = dg%bed_IN(:)

#endif

!.........Compute the x and y components of the external state flow

            dg%QX_EX = ( TY*Q_N_EXT - dg%NY*Q_T_EXT)/(dg%NX*TY - dg%NY*TX)
            dg%QY_EX = (-TX*Q_N_EXT + dg%NX*Q_T_EXT)/(dg%NX*TY - dg%NY*TX)
            
!.........Compute numerical flux
            
            CALL NUMERICAL_FLUX(s,dg,IT)

            DO K = 1,dg%DOFS(dg%EL)

               W_IN = 2.0*dg%M_INV(K,dg%pa)/AREAS(dg%EL)*dg%XLEN(GED)&
              *dg%PHI_EDGE(K,I,LED,dg%pa)*dg%WEGP(I,dg%pa)

               dg%RHS_ZE(K,dg%EL,dg%IRK) = dg%RHS_ZE(K,dg%EL,dg%IRK) - W_IN*F_HAT
               dg%RHS_QX(K,dg%EL,dg%IRK) = dg%RHS_QX(K,dg%EL,dg%IRK) - W_IN*G_HAT
               dg%RHS_QY(K,dg%EL,dg%IRK) = dg%RHS_QY(K,dg%EL,dg%IRK) - W_IN*H_HAT

#ifdef TRACE
               dg%RHS_iota(K,dg%EL,dg%IRK) = dg%RHS_iota(K,dg%EL,dg%IRK) - W_IN*I_HAT
#endif

#ifdef CHEM
               dg%RHS_iota(K,dg%EL,dg%IRK) = dg%RHS_iota(K,dg%EL,dg%IRK) - W_IN*I_HAT
               dg%RHS_iota2(K,dg%EL,dg%IRK) = dg%RHS_iota2(K,dg%EL,dg%IRK) - W_IN*J_HAT
#endif

#ifdef DYNP
               dg%RHS_dynP(K,dg%EL,dg%IRK) = dg%RHS_dynP(K,dg%EL,dg%IRK) - W_IN*I_HAT
#endif

#ifdef SED_LAY

               do ll=1,layers

                  dg%RHS_bed(K,dg%EL,dg%IRK,ll) = dg%RHS_bed(K,dg%EL,dg%IRK,ll) - W_IN*dg%Bed_HAT(ll)

               enddo

#endif

            ENDDO
         ENDDO
 1000 CONTINUE

      RETURN
      END SUBROUTINE EBARRIER_EDGE_HYDRO
