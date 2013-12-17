C***********************************************************************
C     
C     SUBROUTINE WETDRY()
C     
C     This subroutine performs wetting and drying -- for details see:
C     
C     "A Wetting and Drying Treatment for the Runge-Kutta Discontinuous
C     Galerkin Solution to the Shallow Water Equations", Shintaro Bunya
C     et al, Computer Methods in Applied Mechanics and Eng., in review.
C     
C     Written by Shintaro Bunya and Ethan Kubatko
C     01-10-2011 - cem - adapted for p_enrichment and multicomponent
C     
C     
C***********************************************************************

      SUBROUTINE WETDRY()
      
C.....Use appropriate modules
      
      USE SIZES
      USE GLOBAL
      USE DG
      
      IMPLICIT NONE

C.....Declare local variables

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
      REAL(SZ) HT_HAT(NCHECK(ph)), HT_MODE(42)
      REAL(SZ) HT_NODE(NCHECK(ph))
      REAL(SZ) QX_HAT(3), QX_NODE(3), QY_HAT(3), QY_NODE(3)
      REAL(SZ) ZE_HAT(NCHECK(ph),ph), ZE_NODE(NCHECK(ph),ph)
      
C.....Initialize node codes to dry

      DO I = 1,NP
         NODECODE(I) = 0
      ENDDO
      
C.....Loop over the elements

      DO 100 J = 1,NE

         pa = pdg_el(j)
         
#ifdef P0
         if (pa.eq.0) then
            pa = 1
         endif
#endif
         
C.......Retrieve vertex  numbers for element 
C.......(the vertex numbers are the node #'s if p=1)

         No1 = NM(J,1)
         No2 = NM(J,2)
         No3 = NM(J,3)
         
C.......Compute value of ZE and HT at each node, find HT max and min,
C.......and count number of dry nodes

         NDRYNODE = 0
         HT_MIN = HUGE
         HT_MAX = ZERO
         DO I = 1,NCHECK(pa)
            ZE_NODE(I,pa) = ZE(1,J,IRK+1)
            DO K = 2,DOFS(J)
               ZE_NODE(I,pa) = ZE_NODE(I,pa) + PHI_CHECK(K,I,pa)*ZE(K,J,IRK+1)
                                !print*,phi_check(k,i,pa),k,i,pa
            ENDDO
            HT_NODE(I) = ZE_NODE(I,pa) + DP_NODE(I,J,pa)
                                !print*,DP_node(i,j,pa), i,j
#ifdef SED_LAY   
            HT_NODE(I) = 0.D0
            do l=1,layers !notice we are summing over layers at nodes
               DP_NODE(I,J,pa) = bed(1,J,IRK+1,l)
               DO K = 2,DOFS(J)
                  DP_NODE(I,J,pa) = DP_NODE(I,J,pa) + PHI_CHECK(K,I,pa)*bed(K,J,IRK+1,l)
               ENDDO
            enddo
#endif
            HT_NODE(I) = DP_NODE(I,J,pa) + ZE_NODE(I,pa)

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
         
C.......Set "average" bathymetry based on p

         

         !Notice that the average is the average over the sum of the layers
         !DP_AVG = 1.D0/ncheck(pa) * (sum(DP_NODE(:,J,pa)))        
         DP_AVG = C13*(DP_NODE(1,J,pa) + DP_NODE(2,J,pa) + DP_NODE(3,J,pa))

#ifdef P0
         IF (pdg_el(J).EQ.0) THEN
            DP_AVG = 0.D0
            DP_AVG = DPE_MIN(J)
         ENDIF
#endif

C-----------------------------------------------------------------------
C.......Case 1: All nodes are wet
C-----------------------------------------------------------------------

                                !print*,ndrynode

         IF (NDRYNODE.EQ.0) THEN

C.........Case 1A:  Element was previously wet
C----------------------------------------------

            IF (WDFLG(J).EQ.1) THEN
               
C...........Make no adjustments and set element vertex codes to wet

               ELEMENT_CHECK = 0
               NODECODE(No1) = 1
               NODECODE(No2) = 1
               NODECODE(No3) = 1
               
C.........Case 1B:  Element was previously dry
C----------------------------------------------
               
            ELSE
               
C...........Make no adjustments now

               ZE_HAT(:,pa) = ZE_NODE(:,pa)

C...........But set flag to do element check

               ELEMENT_CHECK = 1

            ENDIF
            
C-----------------------------------------------------------------------
C.......Case 2: All nodes are dry
C-----------------------------------------------------------------------

         ELSEIF (NDRYNODE.EQ.NCHECK(pa)) THEN
            
C.........Set surface elevation parallel to bathymetry such that total
C..........water depth at all points = average water depth

            ZE(:,J,IRK+1) = 0.D0
            ZE(1,J,IRK+1) = H0 - DP_AVG
            ZE(2,J,IRK+1) = -HB(2,J,1)
            ZE(3,J,IRK+1) = -HB(3,J,1)

C.........Zero out the fluxes

            QX(:,J,IRK+1) = 0.D0
            QY(:,J,IRK+1) = 0.D0
            
C.........Set flags to dry and skip element check

            WDFLG(J) = 0
            ELEMENT_CHECK = 0

            !zero out higher order stuff
            do i = 4,dofh
               
               ZE(i,J,IRK+1) = 0.D0
               QX(i,J,IRK+1) = 0.D0
               QY(i,J,IRK+1) = 0.D0

c$$$#ifdef SED_LAY !is this more stable?
c$$$               do l=1,layers 
c$$$
c$$$                  bed(i,J,IRK+1,l) = 0.D0
c$$$
c$$$               enddo
c$$$#endif
               
            enddo

            !transported quantities should still be stable
c$$$            if (tracer_flag.eq.1) then
c$$$               
c$$$               iota(:,J,irk+1) = 0.D0
c$$$               
c$$$            endif
c$$$            
c$$$            if (chem_flag.eq.1) then
c$$$               
c$$$               iota(:,J,irk+1) = 0.D0
c$$$               iota2(:,j,irk+1) = 0.D0
c$$$               
c$$$            endif            
            
C----------------------------------------------------------------------------
C.......Case 4: At least one node (at GP) is dry and at least one node is wet
C----------------------------------------------------------------------------

         ELSE

C.........Set averages based on p

            IF (pa.EQ.0) THEN
               HT_AVG = H0
            ELSEIF (pa.EQ.1) THEN
               HT_AVG = ZE(1,J,IRK+1) + DP_AVG
            ELSE
               ZE_VOL = 0.D0
               DO K = 1,DOFS(J)
                  ZE_VOL = ZE_VOL + ZE(K,J,IRK+1)*PHI_INTEGRATED(K,pa)
               ENDDO
               HT_VOL = DP_VOL(J,pa) + 0.25D0*AREAS(J)*ZE_VOL
               HT_AVG = HT_VOL/(0.5D0*AREAS(J))
            ENDIF
            
C     --------------------------------------------------------------
C.........Case 3A: Average water depth < H0
C     --------------------------------------------------------------

            IF (HT_AVG.LE.(H0+ZERO)) THEN

C...........Set surface elevation parallel to bathymetry such that total
C...........water depth at all points = average water depth

               ZE(:,J,IRK+1) = 0.D0
               ZE(1,J,IRK+1) = HT_AVG - DP_AVG
               ZE(2,J,IRK+1) = -HB(2,J,1)
               ZE(3,J,IRK+1) = -HB(3,J,1)
               ZE_HAT(:,pa) = ZE(1,J,IRK+1)

C...........Zero out the fluxes

               QX(:,J,IRK+1) = 0.D0
               QY(:,J,IRK+1) = 0.D0

C     --------------------------------------------------------------
C.........Case 3B: Aaverage water depth > H0
C     --------------------------------------------------------------
               
            ELSE

C...........If p>1 then find a linear function that has the same average
C...........water depth and the same linear components

               IF (pa.GT.1) THEN
                                !Notice the trick here.  Nodes = vertices
                                !since we NOW restrict to a linear regime
                  HT_MODE(1) = HT_AVG
                  HT_MODE(2) = -(ZE(2,J,1) + HB(2,J,1))
                  HT_MODE(3) = -(ZE(3,J,1) + HB(3,J,1))
                  HT_NODE = 0.D0
                  DO K = 1,3
                     HT_NODE(1) = HT_NODE(1) + PHI_CORNER(K,1,pa)*HT_MODE(K)
                     HT_NODE(2) = HT_NODE(2) + PHI_CORNER(K,2,pa)*HT_MODE(K)
                     HT_NODE(3) = HT_NODE(3) + PHI_CORNER(K,3,pa)*HT_MODE(K)
                  ENDDO

               ENDIF

C...........Reorder the nodes=vertices from shallowest to deepest

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

C...........Compute the factors needed for nodal/vertex ajustment

               DH_HAT(L2H(1)) = MAX(H0, HT_NODE(L2H(1))) - HT_NODE(L2H(1))
               DH_HAT(L2H(2)) = MAX(H0, HT_NODE(L2H(2))-0.5*DH_HAT(L2H(1)))
     &              - HT_NODE(L2H(2))
               DH_HAT(L2H(3)) = - DH_HAT(L2H(1)) - DH_HAT(L2H(2))

C...........Adjust nodal/vertex heights, find location of maximum, sum nodal/vertex
C...........fluxes on dry nodes, and count new # of dry nodes/vertices

               NDRYNODE = 0
               QX_SUMs = 0
               QY_SUMs = 0
               HT_MAX = 0
               DO I = 1,3
                  HT_HAT(I)  = HT_NODE(I) + DH_HAT(I)
                  QX_NODE(I) = QX(1,J,IRK+1)
                  QY_NODE(I) = QY(1,J,IRK+1)
                  IF (HT_HAT(I).GT.HT_MAX) THEN
                     HT_MAX_I = I
                     HT_MAX = HT_HAT(I)
                  ENDIF
                  DO K = 2,DOFS(J) !NOTE: We look at integrated value AT the vertex!
                     QX_NODE(I) = QX_NODE(I) + PHI_CORNER(K,I,pa)*QX(K,J,IRK+1)
                     QY_NODE(I) = QY_NODE(I) + PHI_CORNER(K,I,pa)*QY(K,J,IRK+1)
                  ENDDO
                  IF (HT_HAT(I).LE.(H0+1.D-10)) THEN
                     NDRYNODE = NDRYNODE + 1
                     QX_SUMs = QX_SUMs + QX_NODE(I)
                     QY_SUMs = QY_SUMs + QY_NODE(I)
                  ENDIF
               ENDDO

C...........If all adjusted nodes/vertices are dry then zero out fluxes

               IF (NDRYNODE.EQ.3) THEN
                  QX(:,J,IRK+1) = 0.D0
                  QY(:,J,IRK+1) = 0.D0

c$$$  if (tracer_flag.eq.1) then
c$$$  
c$$$  iota(:,J,irk+1) = 0.D0
c$$$  
c$$$  endif
c$$$
c$$$  if (tracer_flag.eq.1) then
c$$$  
c$$$  iota(:,J,irk+1) = 0.D0
c$$$  
c$$$  endif
c$$$  
c$$$  if (chem_flag.eq.1) then
c$$$  
c$$$  iota(:,J,irk+1) = 0.D0
c$$$  iota2(:,j,irk+1) = 0.D0
c$$$  
c$$$  endif 

C...........Else re-distribute nodal/vertex fluxes

               ELSE
                  QX_SUMs = QX_SUMs/(3.D0-REAL(NDRYNODE))
                  QY_SUMs = QY_SUMs/(3.D0-REAL(NDRYNODE))
                  DO I = 1,3
                     IF(HT_HAT(I).LE.(H0+1.D-10)) THEN
                        QX_HAT(I) = 0.D0
                        QY_HAT(I) = 0.D0
                     ELSE
                        IF(QX_SUMs*QX_NODE(I).LT.0.D0) THEN
                           QX_HAT(I) = QX_NODE(I) + QX_SUMs
                        ELSE
                           QX_HAT(I) = QX_NODE(I)
                        ENDIF
                        IF(QY_SUMs*QY_NODE(I).LT.0.D0) THEN
                           QY_HAT(I) = QY_NODE(I) + QY_SUMs
                        ELSE
                           QY_HAT(I) = QY_NODE(I)
                        ENDIF
                     ENDIF
                  ENDDO

C.............Compute new linear modal dofs for x-direction flux

                  QX(1,J,IRK+1) =  C13*(QX_HAT(1) + QX_HAT(2) + QX_HAT(3))
                  QX(2,J,IRK+1) = -C16*(QX_HAT(1) + QX_HAT(2))+C13*QX_HAT(3)
                  QX(3,J,IRK+1) = -0.5D0*QX_HAT(1) + 0.5D0*QX_HAT(2)

C.............Compute new linear modal dofs for y-direction flux

                  QY(1,J,IRK+1) =  C13*(QY_HAT(1) + QY_HAT(2) + QY_HAT(3))
                  QY(2,J,IRK+1) = -C16*(QY_HAT(1) + QY_HAT(2))+C13*QY_HAT(3)
                  QY(3,J,IRK+1) = -0.5D0*QY_HAT(1) + 0.5D0*QY_HAT(2)

               ENDIF

C...........Compute new nodal surface elevations

               ZE_HAT(1:3,pa) = HT_HAT(1:3) - DP_NODE(1:3,J,pa)

C...........Compute the new linear modal dofs for the surface elevation

               ZE(2,J,IRK+1) = -C16*(ZE_HAT(1,pa) + ZE_HAT(2,pa)) + C13*ZE_HAT(3,pa)
               ZE(3,J,IRK+1) = -0.5D0*ZE_HAT(1,pa) + 0.5D0*ZE_HAT(2,pa)

            ENDIF
            

C.........If applicable zero out higher dofs

            do i = 4,dofh
               
               ZE(i,J,IRK+1) = 0.D0
               QX(i,J,IRK+1) = 0.D0
               QY(i,J,IRK+1) = 0.D0

#ifdef SED_LAY !Do we want to zero this out?
               do l=1,layers
                  bed(i,J,irk+1,l) = 0.D0
               enddo
#endif
               
#ifdef TRACE
               iota(i,J,irk+1) = 0.D0
#endif
               
#ifdef CHEM         
               iota(i,J,irk+1) = 0.D0
               iota2(i,j,irk+1) = 0.D0
#endif

#ifdef DYNP
               dynP(i,J,irk+1) = 0.D0
#endif

            enddo


C.........Set flag to do element check

            ELEMENT_CHECK = 1
            
!!!   If element is dry, or partially dry, p is forced down to linears !!!
!!!   NOTE: This effects the global order of the solution if wetdry is on !!!

                                !if (PADAPT.EQ.1) THEN
#ifdef P_AD
            if (pdg_el(j).gt.1) then

               pdg_el(j) = 1
               dofs(J) = 3

            endif
#endif


            
            
         ENDIF

C-----------------------------------------------------------------------
C     End of cases
C-----------------------------------------------------------------------

C.......Element check:
C.......Not clear to me what Shintaro is checking here (EJK)
C.......Setting the lower bound
         
         IF (ELEMENT_CHECK.EQ.1) THEN
            IF (ZE_HAT(HT_MAX_I,pa).GT.(-DPE_MIN(J)+H0)) THEN

C.........Make no adjustments and set element and node flags to wet

               WDFLG(J) = 1
               NODECODE(No1) = 1
               NODECODE(No2) = 1
               NODECODE(No3) = 1

            ELSE

C...........Make adjustments and set element node codes to dry

               QX(:,J,IRK+1) = 0.D0
               QY(:,J,IRK+1) = 0.D0

               WDFLG(J) = 0

            ENDIF
         ENDIF

 100  CONTINUE

      RETURN
      END SUBROUTINE WETDRY
