!***********************************************************************
!     
!     SUBROUTINE WETDRY()
!     
!     This subroutine performs wetting and drying -- for details see:
!     
!     "A Wetting and Drying Treatment for the Runge-Kutta Discontinuous
!     Galerkin Solution to the Shallow Water Equations", Shintaro Bunya
!     et al, Computer Methods in Applied Mechanics and Eng., in review.
!     
!     Written by Shintaro Bunya and Ethan Kubatko
!     01-10-2011 - cem - adapted for p_enrichment and multicomponent
!     
!     
!***********************************************************************

      SUBROUTINE WETDRY()
      
!.....Use appropriate modules
      
      USE SIZES
      USE GLOBAL
      USE DG
      
      IMPLICIT NONE

!.....Declare local variables

      REAL(SZ), PARAMETER :: ZERO = 1.D-6
      REAL(SZ), PARAMETER :: HUGE = 1.D+6

      INTEGER I, J, K1, K2, K, No1, No2, No3, ISTOP,DOFD,l
      INTEGER EL1, EL2, EL3, DRYALL, PWD, DOFWD
      INTEGER ELEMENT_CHECK, HT_MAX_I, HT_MIN_I, NDRYNODE, TMP
      INTEGER L2H(3)
      
      REAL(SZ) AREA
      REAL(SZ) DEN, DP_AVG, DP_MIN
      REAL(SZ) HT_AVG, HT_MAX, HT_MIN, HT_VOL, H0C
      REAL(SZ) NUM1, NUM2
      REAL(SZ) QX_SUMs, QY_SUMs
      REAL(SZ) ZE_VOL

      REAL(SZ) DH_HAT(3)
      REAL(SZ) HT_HAT(dg%NCHECK(dg%ph)), HT_MODE(42)
      REAL(SZ) HT_NODE(dg%NCHECK(dg%ph))
      REAL(SZ) QX_HAT(3), QX_NODE(3), QY_HAT(3), QY_NODE(3)
      REAL(SZ) ZE_HAT(dg%NCHECK(dg%ph),dg%ph), ZE_NODE(dg%NCHECK(dg%ph),dg%ph)

      REAL(SZ), PARAMETER:: HUGEVEL = 1.D+5  
      REAL(SZ):: HS0, U_HAT(3), V_HAT(3)
!.....Initialize node codes to dry

      DO I = 1,NP
         NODECODE(I) = 0
      ENDDO
      
!.....Loop over the elements

      DO 100 J = 1,NE

         dg%pa = pdg_el(j)
         
#ifdef P0
         if (dg%pa.eq.0) then
            dg%pa = 1
         endif
#endif
         
!.......Retrieve vertex  numbers for element 
!.......(the vertex numbers are the node #'s if p=1)

         No1 = NM(J,1)
         No2 = NM(J,2)
         No3 = NM(J,3)
         
!.......Compute value of dg%ZE and HT at each node, find HT max and min,
!.......and count number of dry nodes

         NDRYNODE = 0
         HT_MIN = HUGE
         HT_MAX = ZERO
         DO I = 1,dg%NCHECK(dg%pa)
            !
            ZE_NODE(I,dg%pa) = dg%ZE(1,J,dg%IRK+1)
            DO K = 2,dg%DOFS(J)
               ZE_NODE(I,dg%pa) = ZE_NODE(I,dg%pa) + dg%PHI_CHECK(K,I,dg%pa)*dg%ZE(K,J,dg%IRK+1)
                                ! print*,dg%phi_check(k,i,dg%pa),k,i,dg%pa
            ENDDO
            HT_NODE(I) = ZE_NODE(I,dg%pa) + dg%DP_NODE(I,J,dg%pa)
                                ! print*,dg%DP_node(i,j,dg%pa), i,j
#ifdef SED_LAY   
            HT_NODE(I) = 0.D0
            do l=1,layers !notice we are summing over layers at nodes
               dg%DP_NODE(I,J,dg%pa) = dg%bed(1,J,dg%IRK+1,l)
               DO K = 2,dg%DOFS(J)
                  dg%DP_NODE(I,J,dg%pa) = dg%DP_NODE(I,J,dg%pa) + dg%PHI_CHECK(K,I,dg%pa)*dg%bed(K,J,dg%IRK+1,l)
               ENDDO
            enddo
#endif
            HT_NODE(I) = dg%DP_NODE(I,J,dg%pa) + ZE_NODE(I,dg%pa)

            IF (HT_NODE(I).LT.HT_MIN) THEN
               HT_MIN_I = I
               HT_MIN = HT_NODE(I)
            ENDIF
            IF (HT_NODE(I).GT.HT_MAX) THEN
               HT_MAX_I = I
               HT_MAX = HT_NODE(I)
            ENDIF
            !Note, this can be greater than 3 for element J if p>1
            IF (HT_NODE(I).LE.(H0+ZERO)) NDRYNODE = NDRYNODE + 1
         ENDDO
         
!.......Set "average" bathymetry based on p

         
         !Notice that the average is the average over the sum of the layers
         !DP_AVG = 1.D0/dg%ncheck(dg%pa) * (sum(dg%DP_NODE(:,J,dg%pa)))        
         DP_AVG = dg%C13*(dg%DP_NODE(1,J,dg%pa) + dg%DP_NODE(2,J,dg%pa) + dg%DP_NODE(3,J,dg%pa))

#ifdef P0
         IF (pdg_el(J).EQ.0) THEN
            DP_AVG = 0.D0
            DP_AVG = dg%DPE_MIN(J)
         ENDIF
#endif

!-----------------------------------------------------------------------
!.......Case 1: All nodes are wet
!-----------------------------------------------------------------------

         IF (NDRYNODE.EQ.0) THEN

!.........Case 1A:  Element was previously wet
!----------------------------------------------

            IF (dg%WDFLG(J).EQ.1) THEN
               
!...........Make no adjustments and set element vertex codes to wet

               ELEMENT_CHECK = 0
               NODECODE(No1) = 1
               NODECODE(No2) = 1
               NODECODE(No3) = 1
               
!.........Case 1B:  Element was previously dry
!----------------------------------------------
               
            ELSE
               
!...........Make no adjustments now

               ZE_HAT(:,dg%pa) = ZE_NODE(:,dg%pa)

!...........But set flag to do element check

               ELEMENT_CHECK = 1

            ENDIF
            
!-----------------------------------------------------------------------
!.......Case 2: All nodes are dry
!-----------------------------------------------------------------------

         ELSEIF (NDRYNODE.EQ.dg%NCHECK(dg%pa)) THEN
            
  
!.........Set surface elevation parallel to bathymetry such that total
!.........water depth at all points = average water depth
            IF (dg%pa.EQ.0) THEN
               HT_AVG = H0
            ELSEIF (dg%pa.EQ.1) THEN
               HT_AVG = dg%ZE(1,J,dg%IRK+1) + DP_AVG
            ELSE
               ZE_VOL = 0.D0
               DO K = 1,dg%DOFS(J)
                  ZE_VOL = ZE_VOL + dg%ZE(K,J,dg%IRK+1)*dg%PHI_INTEGRATED(K,dg%pa)
               ENDDO
               HT_VOL = dg%DP_VOL(J,dg%pa) + 0.25D0*AREAS(J)*ZE_VOL
               HT_AVG = HT_VOL/(0.5D0*AREAS(J))
            ENDIF
            
            dg%ZE(:,J,dg%IRK+1) = 0.D0
            dg%ZE(1,J,dg%IRK+1) = HT_AVG - DP_AVG ! for dg%pa = 0, i.e FVM
            IF ( dg%pa > 0 ) THEN 
               ! Prevent water from becoming dangeously low by            !
               ! bumping the water depth back to ZERO, i.e. a small value !
               IF ( HT_AVG < ZERO ) dg%ZE(1,J,dg%IRK+1) = ZERO - DP_AVG 
            END IF
            dg%ZE(2,J,dg%IRK+1) = -dg%HB(2,J,1)
            dg%ZE(3,J,dg%IRK+1) = -dg%HB(3,J,1)

!.........Zero out the fluxes
            dg%QX(:,J,dg%IRK+1) = 0.D0
            dg%QY(:,J,dg%IRK+1) = 0.D0
            
!.........Set flags to dry and skip element check
            dg%WDFLG(J) = 0
            ELEMENT_CHECK = 0

!          
!            !zero out higher order stuff
!            do i = 4,dg%dofh
!               
!               dg%ZE(i,J,dg%IRK+1) = 0.D0
!               dg%QX(i,J,dg%IRK+1) = 0.D0
!               dg%QY(i,J,dg%IRK+1) = 0.D0
!
!$$$#ifdef SED_LAY !is this more stable?
!$$$               do l=1,layers 
!$$$
!$$$                  dg%bed(i,J,dg%IRK+1,l) = 0.D0
!$$$
!$$$               enddo
!$$$#endif
!               
!            enddo
!

            !transported quantities should still be stable
!$$$            if (tracer_flag.eq.1) then
!$$$               
!$$$               dg%iota(:,J,dg%irk+1) = 0.D0
!$$$               
!$$$            endif
!$$$            
!$$$            if (chem_flag.eq.1) then
!$$$               
!$$$               dg%iota(:,J,dg%irk+1) = 0.D0
!$$$               dg%iota2(:,j,dg%irk+1) = 0.D0
!$$$               
!$$$            endif            
            
!----------------------------------------------------------------------------
!.......Case 4: At least one node (at GP) is dry and at least one node is wet
!----------------------------------------------------------------------------

         ELSE

!.........Set averages based on p
            IF (dg%pa.EQ.0) THEN
               HT_AVG = H0
            ELSEIF (dg%pa.EQ.1) THEN
               HT_AVG = dg%ZE(1,J,dg%IRK+1) + DP_AVG
            ELSE
               ZE_VOL = 0.D0
               DO K = 1,dg%DOFS(J)
                  ZE_VOL = ZE_VOL + dg%ZE(K,J,dg%IRK+1)*dg%PHI_INTEGRATED(K,dg%pa)
               ENDDO
               HT_VOL = dg%DP_VOL(J,dg%pa) + 0.25D0*AREAS(J)*ZE_VOL
               HT_AVG = HT_VOL/(0.5D0*AREAS(J))
            ENDIF
            
!     --------------------------------------------------------------
!.........Case 3A: Average water depth < H0
!     --------------------------------------------------------------

            IF (HT_AVG.LE.(H0+ZERO)) THEN

!...........Set surface elevation parallel to bathymetry such that total
!...........water depth at all points = average water depth

               dg%ZE(:,J,dg%IRK+1) = 0.D0
               dg%ZE(1,J,dg%IRK+1) = HT_AVG - DP_AVG 
               IF ( dg%pa > 0 ) THEN
                  ! dw: just in case, it is not likely to happen 
                  IF ( HT_AVG < ZERO ) dg%ZE(1,J,dg%IRK+1) = ZERO - DP_AVG

               END IF
               dg%ZE(2,J,dg%IRK+1) = -dg%HB(2,J,1)
               dg%ZE(3,J,dg%IRK+1) = -dg%HB(3,J,1)
               ZE_HAT(:,dg%pa) = dg%ZE(1,J,dg%IRK+1)

!...........Zero out the fluxes
               dg%QX(:,J,dg%IRK+1) = 0.D0
               dg%QY(:,J,dg%IRK+1) = 0.D0

!     --------------------------------------------------------------
!.........Case 3B: Aaverage water depth > H0
!     --------------------------------------------------------------
               
            ELSE

!...........If p>1 then find a linear function that has the same average
!...........water depth and the same linear components

               IF (dg%pa.GT.1) THEN
                                !Notice the trick here.  Nodes = vertices
                                !since we NOW restrict to a linear regime
                  !c--- dw: should it be (?)
                  ! HT_MODE(2) = -(dg%ZE(2,J,1) + dg%HB(2,J,1))
                  ! HT_MODE(3) = -(dg%ZE(3,J,1) + dg%HB(3,J,1))
                  HT_MODE(1) = HT_AVG
                  HT_MODE(2) = -(dg%ZE(2,J,1) + dg%HB(2,J,1))
                  HT_MODE(3) = -(dg%ZE(3,J,1) + dg%HB(3,J,1))
                  HT_NODE = 0.D0
                  DO K = 1,3
                     HT_NODE(1) = HT_NODE(1) + dg%PHI_CORNER(K,1,dg%pa)*HT_MODE(K)
                     HT_NODE(2) = HT_NODE(2) + dg%PHI_CORNER(K,2,dg%pa)*HT_MODE(K)
                     HT_NODE(3) = HT_NODE(3) + dg%PHI_CORNER(K,3,dg%pa)*HT_MODE(K)
                  ENDDO

               ENDIF

!...........Reorder the nodes=vertices from shallowest to deepest

               L2H(1) = 1
               L2H(2) = 2
               L2H(3) = 3
               DO K1 = 1,3
                  DO K2 = K1,3
                     IF ( HT_NODE(L2H(K1)).GT.HT_NODE(L2H(K2)) ) THEN
                        TMP = L2H(K1)
                        L2H(K1) = L2H(K2)
                        L2H(K2) = TMP
                     ENDIF
                  ENDDO
               ENDDO

!...........Compute the factors needed for nodal/vertex ajustment

               DH_HAT(L2H(1)) = MAX(H0, HT_NODE(L2H(1))) - HT_NODE(L2H(1))
               DH_HAT(L2H(2)) = MAX(H0, HT_NODE(L2H(2))-0.5*DH_HAT(L2H(1)))&
                    - HT_NODE(L2H(2))
               DH_HAT(L2H(3)) = - DH_HAT(L2H(1)) - DH_HAT(L2H(2))

!...........Adjust nodal/vertex heights, find location of maximum, sum nodal/vertex
!...........fluxes on dry nodes, and count new # of dry nodes/vertices

               NDRYNODE = 0
               QX_SUMs = 0
               QY_SUMs = 0
               HT_MAX = 0
               DO I = 1,3
                  HT_HAT(I)  = HT_NODE(I) + DH_HAT(I)
                  QX_NODE(I) = dg%QX(1,J,dg%IRK+1)
                  QY_NODE(I) = dg%QY(1,J,dg%IRK+1)
                  IF (HT_HAT(I).GT.HT_MAX) THEN
                     HT_MAX_I = I
                     HT_MAX = HT_HAT(I)
                  ENDIF
                  DO K = 2,dg%DOFS(J) !NOTE: We look at integrated value AT the vertex!
                     QX_NODE(I) = QX_NODE(I) + dg%PHI_CORNER(K,I,dg%pa)*dg%QX(K,J,dg%IRK+1)
                     QY_NODE(I) = QY_NODE(I) + dg%PHI_CORNER(K,I,dg%pa)*dg%QY(K,J,dg%IRK+1)
                  ENDDO
                  IF (HT_HAT(I).LE.(H0+1.D-10)) THEN
                     NDRYNODE = NDRYNODE + 1
                     QX_SUMs = QX_SUMs + QX_NODE(I)
                     QY_SUMs = QY_SUMs + QY_NODE(I)
                  ENDIF
               ENDDO

!...........If all adjusted nodes/vertices are dry then zero out fluxes

               IF (NDRYNODE.EQ.3) THEN
                  
                  dg%QX(:,J,dg%IRK+1) = 0.D0
                  dg%QY(:,J,dg%IRK+1) = 0.D0

!$$$  if (tracer_flag.eq.1) then
!$$$  
!$$$  dg%iota(:,J,dg%irk+1) = 0.D0
!$$$  
!$$$  endif
!$$$
!$$$  if (tracer_flag.eq.1) then
!$$$  
!$$$  dg%iota(:,J,dg%irk+1) = 0.D0
!$$$  
!$$$  endif
!$$$  
!$$$  if (chem_flag.eq.1) then
!$$$  
!$$$  dg%iota(:,J,dg%irk+1) = 0.D0
!$$$  dg%iota2(:,j,dg%irk+1) = 0.D0
!$$$  
!$$$  endif 

!...........Else re-distribute nodal/vertex fluxes

               ELSE
                  !
                  QX_SUMs = QX_SUMs/(3.D0-REAL(NDRYNODE))
                  QY_SUMs = QY_SUMs/(3.D0-REAL(NDRYNODE))
                  DO I = 1,3
                     IF(HT_HAT(I).LE.(H0+1.D-10)) THEN
                        QX_HAT(I) = 0.D0
                        QY_HAT(I) = 0.D0
                     ELSE
                        !c
                        !c... dw: seem to not agree with sb's paper (??)
                        !c       (momemtum could be increased)
                        !IF(QX_SUMs*QX_NODE(I).LT.0.D0) THEN
                        !   QX_HAT(I) = QX_NODE(I) + QX_SUMs
                        !ELSE
                        !   QX_HAT(I) = QX_NODE(I)
                        !ENDIF
                        !IF(QY_SUMs*QY_NODE(I).LT.0.D0) THEN
                        !   QY_HAT(I) = QY_NODE(I) + QY_SUMs
                        !ELSE
                        !   QY_HAT(I) = QY_NODE(I)
                        !ENDIF

                        !c
                        !c--- dw: make it consistant with sb's paper
                        !c
                        QX_HAT(I) = QX_NODE(I) + QX_SUMs 
                        QY_HAT(I) = QY_NODE(I) + QY_SUMs  
                     ENDIF
                  ENDDO
                  !
                  
!c                  !c............. dw, consider cliping velocity greater than HUGEVEL 
!c                  !c               (experimental)
!c                  U_HAT(1:3) = QX_HAT(1:3)/HT_HAT(1:3) ;
!c                  V_HAT(1:3) = QY_HAT(1:3)/HT_HAT(1:3) ;
!c                  DO I = 1, 3
!c                     IF ( abs(U_HAT(I)) > HUGEVEL ) U_HAT(I) = HUGEVEL ;
!c                     IF ( abs(V_HAT(I)) > HUGEVEL ) V_HAT(I) = HUGEVEL ;
!c                  END 
!c                  QX_HAT(1:3) = U_HAT(1:3)*HT_HAT(1:3) ;
!c                  QY_HAT(1:3) = U_HAT(1:3)*HT_HAT(1:3) ;
!c                  !c............



!.............Compute new linear modal dg%dofs for x-direction flux

                  dg%QX(1,J,dg%IRK+1) =  dg%C13*(QX_HAT(1) + QX_HAT(2) + QX_HAT(3))
                  dg%QX(2,J,dg%IRK+1) = -dg%C16*(QX_HAT(1) + QX_HAT(2))+dg%C13*QX_HAT(3)
                  dg%QX(3,J,dg%IRK+1) = -0.5D0*QX_HAT(1) + 0.5D0*QX_HAT(2)

!.............Compute new linear modal dg%dofs for y-direction flux

                  dg%QY(1,J,dg%IRK+1) =  dg%C13*(QY_HAT(1) + QY_HAT(2) + QY_HAT(3))
                  dg%QY(2,J,dg%IRK+1) = -dg%C16*(QY_HAT(1) + QY_HAT(2))+dg%C13*QY_HAT(3)
                  dg%QY(3,J,dg%IRK+1) = -0.5D0*QY_HAT(1) + 0.5D0*QY_HAT(2)
               ENDIF

!...........Compute new nodal surface elevations

               ZE_HAT(1:3,dg%pa) = HT_HAT(1:3) - dg%DP_NODE(1:3,J,dg%pa)

!...........Compute the new linear modal dg%dofs for the surface elevation

               dg%ZE(2,J,dg%IRK+1) = -dg%C16*(ZE_HAT(1,dg%pa) + ZE_HAT(2,dg%pa)) + dg%C13*ZE_HAT(3,dg%pa)
               dg%ZE(3,J,dg%IRK+1) = -0.5D0*ZE_HAT(1,dg%pa) + 0.5D0*ZE_HAT(2,dg%pa)

            ENDIF
            

!.........If applicable zero out higher dg%dofs

            do i = 4,dg%dofh
               
               dg%ZE(i,J,dg%IRK+1) = 0.D0
               dg%QX(i,J,dg%IRK+1) = 0.D0
               dg%QY(i,J,dg%IRK+1) = 0.D0

!Do we want to zero this out?
#ifdef SED_LAY 
               do l=1,layers
                  dg%bed(i,J,dg%irk+1,l) = 0.D0
               enddo
#endif
               
#ifdef TRACE
               dg%iota(i,J,dg%irk+1) = 0.D0
#endif
               
#ifdef CHEM         
               dg%iota(i,J,dg%irk+1) = 0.D0
               dg%iota2(i,j,dg%irk+1) = 0.D0
#endif

#ifdef DYNP
               dg%dynP(i,J,dg%irk+1) = 0.D0
#endif

            enddo


!.........Set flag to do element check

            ELEMENT_CHECK = 1
            
!!!   If element is dry, or partially dry, p is forced down to linears !!!
!!!   NOTE: This effects the global order of the solution if wetdry is on !!!

                                !if (dg%PADAPT.EQ.1) THEN
#ifdef P_AD
            if (pdg_el(j).gt.1) then

               pdg_el(j) = 1
               dg%dofs(J) = 3

            endif
#endif


            
            
         ENDIF

!-----------------------------------------------------------------------
!     End of cases
!-----------------------------------------------------------------------

!.......Element check:
!.......Not clear to me what Shintaro is checking here (EJK)
!.......Setting the lower bound
         
         IF (ELEMENT_CHECK.EQ.1) THEN

!Now test to see if the water column at the min vertex is wet
            IF (MINVAL(ZE_HAT(:,dg%pa)+dg%DP_NODE(:,J,dg%pa)).GT.H0 ) THEN  ! cem fix
!            IF (ZE_HAT(HT_MAX_I,dg%pa).GT.(-dg%DPE_MIN(J)+H0)) THEN     ! Shintaro's criteria 

!.........Make no adjustments and set element and node flags to wet

               dg%WDFLG(J) = 1
               NODECODE(No1) = 1
               NODECODE(No2) = 1
               NODECODE(No3) = 1

            ELSE

!...........Make adjustments and set element node codes to dry

               dg%QX(:,J,dg%IRK+1) = 0.D0
               dg%QY(:,J,dg%IRK+1) = 0.D0

               dg%WDFLG(J) = 0

            ENDIF
         ENDIF

 100  CONTINUE

      RETURN
      END SUBROUTINE WETDRY
