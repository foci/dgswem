!***********************************************************************
!     
!     SUBROUTINE DG_PREP()
!     
!     This subroutine does preparatory stuff for DG
!     
!     Written by Ethan Kubatko (03-07-2005)
!     
!-----------------------------------------------------------------------
!     
!     mod history of hp_ADCIRC since v9
!     
!     v9_sb1       - Aug 05 - sb - parallelized
!     v9_sb2       - Aug 05 - sb - wetting/drying
!     v10_sb1      - Aug 05 - sb - reflect fort.12as an initial surface elevation
!     v10_sb5      - Oct 03 - sb - consolidate Ethan's slope limiter
!     v10_sb5      - Oct 24 - sb - DG.65 output is added
!     - Aug 29, 2007 - Fixed bug for startdry = 1
!     01-10-2011 - cem - adapted for p_enrichment and multicomponent
!     06-01-2012 - cem -sediment added
!     
!***********************************************************************

      SUBROUTINE PREP_DG(s)

!.....Use appropriate modules
      USE SIZES
      USE GLOBAL
      USE DG
      USE NodalAttributes, ONLY : STARTDRY, FRIC, GeoidOffset,LoadGeoidOffset,LoadManningsN,ManningsN
#ifdef CMPI
      USE MESSENGER_ELEM
      USE MESSENGER
#endif

      IMPLICIT NONE
      
      type (sizes_type) :: s
      

!.....Declare local variables

      INTEGER II, l, P_0, DOF_0,j,k,kk,jj,i,chi,ll,mm
      REAL(SZ) AREA, ANGLE_SUM, HBB(3), CASUM, DP_MIN,temp_lay
      REAL(SZ) XI, YI, ZE1, ZE2, ZE3, l2er,l2erh2,xcen,ycen,epsl,pi_n
      REAL(SZ) ZP(3), DHBX, ell_1,ell_2,ell_3,int_hb,int_ze,int_yd

      real(sz) checkarea,arint(2,2),rhsint(2),edgeint,dpsdx,psimid
      real(sz) determ,sfacdub2max,sfacdub3max,R
      integer i1,i2,sfac_flag,led
      integer ifac2max,ifac3max

      real(sz) xmid,ymid,Ox,Oy,rev,C_0,sig,C_1

      Real(SZ),allocatable :: XBCbt(:),YBCbt(:),radial(:),XB(:),YB(:),l2e(:)
      Real(SZ),allocatable :: iota_check(:),iota_check2(:),hbo(:,:,:), ydubo(:,:)
      Real(SZ),allocatable :: YELEM(:),YED(:),HB1(:,:,:,:), zeo(:,:,:)

      integer :: myproc_here

      Allocate ( XBCbt(S%MNE),YBCbt(S%MNE),radial(S%MNE),XB(S%MNE),YB(S%MNE),l2e(S%MNE) )
      Allocate ( iota_check(S%MNE),iota_check2(S%MNE),hbo(36,S%MNE,1),ydubo(36,s%mne) ) 
      Allocate ( YELEM(dg%ph),YED(dg%ph),hb1(36,s%mne,1,dg%ph), zeo(36,s%mne,1) )


      C13 = 1.D0/3.D0
      C16 = 1.D0/6.D0
      R = 6378206.4d0

#ifdef CMPI
      myproc_here = myproc
#else
      myproc_here = 0
#endif


!     sb-PDG1 moved from other places

!.....Obtain RK time scheme parameters

      CALL RK_TIME()
      
!.....Compute the degrees of freedom per element

      DG%DOF = (dg%pl+1)*(dg%pl+2)/2
      dg%dofx = (dg%px+1)*(dg%px+2)/2    ! dg%dofx for variable functions f=f(x) 
      P_0 = dg%pl
      DOF_0 = (dg%pl+1)*(dg%pl+2)/2   ! dof at lowest order when p!=0
      dg%dofh = (dg%ph + 1)*(dg%ph + 2)/2


!.....Allocate some DG stuff

      IF (DG%PADAPT.EQ.1) THEN

         dg%dofh = (dg%ph + 1)*(dg%ph + 2)/2
         dg%dofl = (dg%pl + 1)*(dg%pl + 2)/2
         pa = dg%pl
           
      elseif (dg%padapt.eq.0) then 

         dg%dofh = dg%dofh
         dg%dofl = DOF_0
         pa = dg%pl

      endif

#ifdef SED_LAY
!.....Initialize funtion parser for sediment types
      init_parser = .false.
#endif

#ifdef SED_LAY
!.....Initialize stabilizer sweep for A.D.
      stblzr = .false.
#endif


!.....Compute the number of gauss points needed for the edge integrals

      CALL ALLOC_DG4(s)          !moved here 6.28.10, for p_adapt because of messenger_elem      

         dofs(:) = dg%dofl
         PDG_EL(:) = dg%pl
         PDG(:) = dg%pl  
         PCOUNT(:) = 0
         pa = dg%pl

      do chi=dg%pl,dg%ph
         NEGP(chi) = chi + 1
      enddo

      IF (dg%pl.eq.0) THEN

         PDG_EL(:) = 1
         PDG(:) = 1     
         DG%DOF    = 3
         dg%pl     = 1
         dg%dofl   = 3
         P_0    = 0
         DOF_0  = 1 
         NEGP(dg%pl) = 2

      ENDIF

!     cnd
!     iwrite=0
      
!.....Initilization for parallel DG run

#ifdef CMPI

      CALL MSG_TYPES_ELEM()     ! Determine Word Sizes for Message-Passing
      CALL MSG_TABLE_ELEM()     ! Read Message-Passing Tables

      IF (DG%SLOPEFLAG.ge.4) THEN
         CALL MSG_TYPES()
         CALL MSG_TABLE()
      ENDIF

#endif

!.....Create the edge based data

      IF(MYPROC_HERE.EQ.0) THEN
         PRINT*, 'CREATING EDGE DATA...'
         PRINT*, ''
      ENDIF
      CALL CREATE_EDGE_DATA(s)
      IF(MYPROC_HERE.EQ.0) THEN
         print *, 'CREATING EDGE DATA DONE'
         print *, ''
      ENDIF

#ifdef CMPI
      CALL MESSAGE_START_ELEM() ! Startup persistent message passing
      IF (DG%SLOPEFLAG.ge.4) CALL MESSAGE_START()
#endif

!.....Re-arrange elevation specified boundary segment data for DG

      IF (NEEDS.GT.0) THEN
         CALL ALLOC_DG1(s%MNBFR)
         II = 1
         JJ = 1
         DO I = 1,NBFR
            DO J = 1,NOPE
               DO K = 1,NVDLL(J)-1
                  EMO_DG(I,II,1) = EMO(I,JJ)
                  EMO_DG(I,II,2) = EMO(I,JJ+1)
                  EFA_DG(I,II,1) = EFA(I,JJ)
                  EFA_DG(I,II,2) = EFA(I,JJ+1)
                  UMO_DG(I,II,1) = UMO(I,JJ)
                  UMO_DG(I,II,2) = UMO(I,JJ+1)
                  UFA_DG(I,II,1) = UFA(I,JJ)
                  UFA_DG(I,II,2) = UFA(I,JJ+1)
                  VMO_DG(I,II,1) = VMO(I,JJ)
                  VMO_DG(I,II,2) = VMO(I,JJ+1)
                  VFA_DG(I,II,1) = VFA(I,JJ)
                  VFA_DG(I,II,2) = VFA(I,JJ+1)
                  II = II + 1
                  JJ = JJ + 1
               ENDDO
               JJ = JJ + 1
            ENDDO
            II = 1
            JJ = 1
         ENDDO
      ENDIF

!.....Re-arrange non-zero flow specified boundary segment data for DG

      IF (NFEDS.GT.0) THEN
         CALL ALLOC_DG2(s%MNFFR)
         II = 1
         JJ = 1
         DO I = 1,NFFR
            DO J = 1,NBOU
               IF ( (SEGTYPE(J).EQ.2 ).OR.(SEGTYPE(J).EQ.12)&
              .OR.(SEGTYPE(J).EQ.22) ) THEN
                  DO K = 1,NVELL(J)-1
                     QNAM_DG(I,II,1) = QNAM(I,JJ)
                     QNAM_DG(I,II,2) = QNAM(I,JJ+1)
                     QNPH_DG(I,II,1) = QNPH(I,JJ)
                     QNPH_DG(I,II,2) = QNPH(I,JJ+1)
                     II = II + 1
                     JJ = JJ + 1
                  ENDDO
                  JJ = JJ + 1
               ENDIF
            ENDDO
            II = 1
            JJ = 1
         ENDDO
      ENDIF
      
!.....If there are internal barriers allocate some stuff

      IF (NIBEDS.NE.0) CALL ALLOC_DG3(S%MNP)

!.....Allocate the array for node to element table

      CALL ALLOC_NNOEL1(S)

!.....Determine the number of elements connected at each node

      EL_COUNT = 0
      MAXEL = 1
      DO K = 1,3
         DO J = 1,S%MNE
            N1 = NM(J,K)
            EL_COUNT(N1) = EL_COUNT(N1) + 1
         ENDDO
      ENDDO
      MAXEL = MAXVAL(EL_COUNT)

!.....Allocate the array for the node to element table

      CALL ALLOC_NNOEL2(S,MAXEL)
      
!.....Construct node to element table

      EL_COUNT = 0
      DO K = 1,3
         DO J = 1,S%MNE
            N1 = NM(J,K)
            NNOEL(N1,1+EL_COUNT(N1)) = J
            EL_COUNT(N1) = EL_COUNT(N1) + 1
         ENDDO
      ENDDO

!.....Construct node to element angle table

      DO I = 1,S%MNP
         ETAMAX(I) = -99999
         KK = 1
         ELETAB(I,1) = I
         S1 = SFAC(I)
         J1 = NEITAB(I,1)
         DO 111 K = 1,NNEIGH(I)-1
            J2 = NEITAB(I,1+K)
            IF (K.LT.(NNEIGH(I)-1)) THEN
               J3 = NEITAB(I,2+K)
            ELSE
               J3 = NEITAB(I,2)
            ENDIF
            DO J = 1,EL_COUNT(I)
               DG%EL = NNOEL(I,J)
               N1 = NM(DG%EL,1)
               N2 = NM(DG%EL,2)
               N3 = NM(DG%EL,3)
               IF ((J1.EQ.N1).OR.(J1.EQ.N2).OR.(J1.EQ.N3)) THEN
                  IF ((J2.EQ.N1).OR.(J2.EQ.N2).OR.(J2.EQ.N3)) THEN
                     IF ((J3.EQ.N1).OR.(J3.EQ.N2).OR.(J3.EQ.N3)) THEN
                        ELETAB(I,1+KK) = DG%EL
                        S2  = SFAC(J2)
                        SAV = (S1 + S2)/2.D0
                        VEC1(1) =      X(J1) - X(J2)
                        VEC1(2) = SAV*(Y(J1) - Y(J2))
                        VEC2(1) =      X(J1) - X(J3)
                        VEC2(2) = SAV*(Y(J1) - Y(J3))
                        MAG1 = SQRT(VEC1(1)**2 + VEC1(2)**2)
                        MAG2 = SQRT(VEC2(1)**2 + VEC2(2)**2)
                        DOT = DOT_PRODUCT(VEC1,VEC2)
                        EL_ANG  = ACOS(DOT/(MAG1*MAG2))
                        ANGTAB(I,KK+1) = RAD2DEG*EL_ANG
                        KK = KK + 1
                        GOTO 111
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO
 111     CONTINUE
      ENDDO


!     C.....Allocate some DG stuff

!     CALL ALLOC_DG4()
      
!.....Initialize the DG arrays

      ZE = 0.D0
      zeo = 0.D0
      QX = 0.D0
      QY = 0.D0
      hbo = 0.D0
      hb = 0.D0
      LZ = 0.d0
      HZ = 0.D0
      TZ = 0.D0
      MZ = 0.D0
      iota = 0.D0
      iota2 = 0.D0
      iotaa = 0.D0
      iotaa2 = 0.D0 
      iotaa3 = 0.D0       
      MassMax = 0.D0
      MARK = 0

      bed = 0.D0

      RHS_ZE = 0.D0
      RHS_QX = 0.D0
      RHS_QY = 0.D0
      RHS_iota = 0.D0
      RHS_iota2 = 0.D0
      
      WSX2(:) = 0
      WSY2(:) = 0

!.....If using modal initial conditions transform the bathymetry from
!.....nodal coordinates to modal dof

      DO J = 1,S%MNE
         N1 = NM(J,1)
         N2 = NM(J,2)
         N3 = NM(J,3)
         hbo(1,J,1) =  1.D0/3.D0 * (DP(N1) + DP(N2) + DP(N3))
         hbo(2,J,1) = -1.D0/6.D0 * (DP(N1) + DP(N2)) + 1.D0/3.D0*DP(N3)
         hbo(3,J,1) = -0.5D0*DP(N1) + 0.5D0*DP(N2)

         ydubo(1,J)= 1.D0/3.D0*(Y(N1) + Y(N2) + Y(N3))
         ydubo(2,J) = -1.D0/6.D0*(Y(N1) + Y(N2))&
        + 1.D0/3.D0*Y(N3)
         ydubo(3,J) = -0.5D0*Y(N1) + 0.5D0*Y(N2)

      ENDDO


      if (LoadManningsN) then
         DO J = 1,NE
            N1 = NM(J,1)
            N2 = NM(J,2)
            N3 = NM(J,3)
            MANN(1,J) =  1.D0/3.D0*(ManningsN(N1)&
           + ManningsN(N2) + ManningsN(N3))
            MANN(2,J) = -1.D0/6.D0*(ManningsN(N1) &
           + ManningsN(N2)) + 1.D0/3.D0*ManningsN(N3)
            MANN(3,J) = -0.5D0*ManningsN(N1) + 0.5D0*ManningsN(N2)
         ENDDO
      endif

      IF (DG%MODAL_IC.EQ.0) THEN
!     this assumes a cold start
         if (LoadGeoidOffset) then
            DO J = 1,NE
               N1 = NM(J,1)
               N2 = NM(J,2)
               N3 = NM(J,3)
               zeo(1,J,1)=1.d0/3.d0*(GeoidOffset(N1)+GeoidOffset(N2)+&
              GeoidOffset(N3))
               IF (dof_0.NE.1) THEN
                  zeo(2,J,1)=-1.d0/6.d0*(GeoidOffset(N1)+GeoidOffset(N2))&
                 +1.d0/3.d0*GeoidOffset(N3)
                  zeo(3,J,1)=-.5d0*GeoidOffset(N1)+.5d0*GeoidOffset(N2)
               ENDIF
            ENDDO
         endif
      ENDIF

!--   

!As part of initializing the system, let us determine the partials
!of the sediment discharge equation, fed in by fort.dg

#ifdef SED_LAY
      IF(MYPROC_HERE.EQ.0)THEN
         print*, 'Parsing the following sediment discharge equations:'
         print *, ''
         print*, 'In X we have: ', sed_equationX
         print*, 'In Y we have: ', sed_equationY
         open(444, file = "./sedlaw.X")
         write(444,'(a)') sed_equationX
         close(444)
         open(445, file = "./sedlaw.Y")
         write(445,*) sed_equationY
         close(445)
         CALL SYSTEM('python py_scriptX') !this writes db_partials_X file
         CALL SYSTEM('python py_scriptY') !this writes db_partials_Y file
         print *, ''
      ENDIF
#endif

      IF(MYPROC_HERE.EQ.0) THEN
         print *, 'PREP FOR WET/DRY BEGINS...'
      ENDIF

!.....1. Set initial surface elevation above the bed elevation
!.....if wetting-and-drying is enabled and the initial water depth is
!.....not specified by fort.12
!.....2. Set wet-and-dry elemental flags
!.....3. Set the DOF at dry elements = 1

      NCHECK(1) = 3
      if (dg%ph.gt.1) then
         do chi = 2,dg%ph
            NCHECK(chi) = NCHECK(1) + 3*negp(chi)
         enddo
      endif
      
      CALL ALLOC_DG_WETDRY(s)
      PHI_CHECK = 0.D0
      PSI_CHECK = 0.D0
      H0L = H0
      H0H = H0 * 1.0
      HABSMIN = H0 * 1.0
            
!.....Retrieve the normals to the edges

      CALL CALC_NORMAL()

!.....Retrieve the area integral gauss quadrature points
      
      do j=1,dg%ph
         
         CALL QUAD_PTS_AREA(s,2*j,j)

      enddo

!.....Retrieve the edge integral gauss quadrature points
      
      do j=1,dg%ph

         CALL QUAD_PTS_EDGE(s,j,j)

      enddo

!.....Evaluate the orthogonal basis and its derivatives at the area
!.....gauss quadrature points

      do j=1,dg%ph

         CALL ORTHOBASIS_AREA(j,j)

      enddo

!.....Evaluate the orthogonal basis at the edge gauss quadrature points

      do j=1,dg%ph

         CALL ORTHOBASIS_EDGE(j,j)

      enddo

!.....Do the L2-projection of the initial conditions

      hb = 0.D0
      qx = 0.D0
      qy = 0.D0
      ze = 0.D0
      ydub = 0.d0
      hb1 = 0.D0

 
      hb(1:dg%dofh,:,1) = hbo(1:dg%dofh,:,1)
      ze(1:dg%dofh,:,1) = zeo(1:dg%dofh,:,1)
      ydub(1:dg%dofh,:,1) = ydubo(1:dg%dofh,:)

      do chi=dg%pl,dg%ph
         hb1(1:dg%dofh,:,1,chi) = hbo(1:dg%dofh,:,1)
      enddo

      !if layers are on, distribute them evenly across the total bed load
      !other approaches are clearly available, this is a simple choice
      !adapt for higher order initial data

#ifdef SED_LAY

      do ll=1,layers
         
         bed(:,:,1,ll) = hb(:,:,1) / layers
         
      enddo

#endif      

#ifdef ARTDIF
      !Set up the artififical diffusion stuff
      e1(:) = 0.D0
      balance(:) = 0.D0
      entrop(:,:) = -100.D0

      if (dg%tune_by_hand.eq.1) then

         balance(4) = 0.D0     
         
         s0    =  0.0D0       
         kappa =   -1.D0      
         
         e1(1) = 0.D0
         e1(2) = 0.D0
         e1(3) = 0.D0
         e1(4) = 2.5e-6       
         e1(5) = 0.D0
         
      else
         
         e1 = uniform_dif

      endif

#endif

!Update DPE_MIN

      DO J = 1,NE
         DPE_MIN(J) = MIN(DP(NM(J,1)),DP(NM(J,2)),DP(NM(J,3)))
      ENDDO

!.....Compute the values of the nodal basis functions at the
!.....area gauss quadrature points, at every p level chi

      do chi = 1,dg%ph
         do I = 1,NAGP(chi)
            PSI1(I,chi) = -1.D0/2.D0*(XAGP(I,chi) + YAGP(I,chi))
            PSI2(I,chi) =  1.D0/2.D0*(XAGP(I,chi) + 1.D0)
            PSI3(I,chi) =  1.D0/2.D0*(YAGP(I,chi) + 1.D0)
         enddo
      enddo
      
!.....Store the derivatives of the (linear) nodal basis functions

      DRPSI(1) = -1.D0/2.D0
      DRPSI(2) =  1.D0/2.D0
      DRPSI(3) =  0.D0
      DSPSI(1) = -1.D0/2.D0
      DSPSI(2) =  0.D0
      DSPSI(3) =  1.D0/2.D0
      

!.....Pre-compute the derivatives of the coordinate transformation for
!.....each element

      DO J = 1,S%MNE

!.....Retrieve the global node numbers for the element

         N1 = NM(J,1)
         N2 = NM(J,2)
         N3 = NM(J,3)
         x1=x(n1)
         y1=y(n1)
         x2=x(n2)
         y2=y(n2)
         x3=x(n3)
         y3=y(n3)
         AREA = (X1 - X3)*(Y2 - Y3) + (X3 - X2)*(Y1 - Y3)
         area=area*.5d0

!.....Compute the derivatives of the coordinate transformation

         DRDX(J) = 1.D0/AREA*(Y(N3) - Y(N1))
         DSDX(J) = 1.D0/AREA*(Y(N1) - Y(N2))

         DRDY(J) = 1.D0/AREA*(X(N1) - X(N3))
         DSDY(J) = 1.D0/AREA*(X(N2) - X(N1))
         
!.......Compute elemental Coriolis and friction terms

         CORI_EL(J) = (CORIF(N1) + CORIF(N2) + CORIF(N3))/3.D0
         FRIC_EL(J) = (FRIC(N1) + FRIC(N2) + FRIC(N3))/3.D0

!.......Pre-compute the bathymetry and the gradient of the bathymetry at
!.......the quadrature points and compute volume of water

         DP_VOL(J,:) = 0.D0
         SFAC_ELEM(:,J,:) = 0.D0

         do chi = 1,dg%ph

            if (chi.gt.1) then

               do I = 1,NAGP(chi) ! Forced for wetting and drying's sake

                  BATH(I,J,chi) = 0.D0
                  DBATHDX(I,J,chi) = 0.D0
                  DBATHDY(I,J,chi) = 0.D0
                  YELEM(chi) = 0.d0

                  DO K = 1,(chi+1)*(chi+2)/2 

                     DPHIDX = DRPHI(K,I,chi)*DRDX(J) + DSPHI(K,I,chi)*DSDX(J)
                     DPHIDY = DRPHI(K,I,chi)*DRDY(J) + DSPHI(K,I,chi)*DSDY(J)
                     XFAC(K,I,J,chi)  = M_INV(K,chi)*WAGP(I,chi)*DPHIDX
                     YFAC(K,I,J,chi)  = M_INV(K,chi)*WAGP(I,chi)*DPHIDY
                     SRFAC(K,I,J,chi) = M_INV(K,chi)*WAGP(I,chi)*PHI_AREA(K,I,chi)
                     BATH(I,J,chi) = BATH(I,J,chi) + HB(K,J,1)*PHI_AREA(K,I,chi)
                     YELEM(chi) = YELEM(chi) + YDUB(K,J,chi)*PHI_AREA(K,I,chi)
                     
                     IF (ICS.EQ.1) THEN
                        SFAC_ELEM(I,J,chi)=1.0D0
                     ELSE
                        SFAC_ELEM(I,J,chi)=COS(SFEA0)/COS(YELEM(chi)/R)
                     ENDIF

                     DBATHDX(I,J,chi) = DBATHDX(I,J,chi) + HB(K,J,1)*DPHIDX
                     DBATHDY(I,J,chi) = DBATHDY(I,J,chi) + HB(K,J,1)*DPHIDY
                     DP_VOL(J,chi) = DP_VOL(J,chi) + WAGP(I,chi)*HB(K,J,1)*PHI_AREA(K,I,chi)

                  ENDDO

               enddo

            else

               do I = 1,NAGP(chi) ! Area quadrature points

                  BATH(I,J,chi) = 0.D0
                  DBATHDX(I,J,chi) = 0.D0
                  DBATHDY(I,J,chi) = 0.D0
                  YELEM(chi) = 0.d0

                  DO K = 1,(chi+1)*(chi+2)/2 

                     DPHIDX = DRPHI(K,I,chi)*DRDX(J) + DSPHI(K,I,chi)*DSDX(J)
                     DPHIDY = DRPHI(K,I,chi)*DRDY(J) + DSPHI(K,I,chi)*DSDY(J)
                     XFAC(K,I,J,chi)  = M_INV(K,chi)*WAGP(I,chi)*DPHIDX
                     YFAC(K,I,J,chi)  = M_INV(K,chi)*WAGP(I,chi)*DPHIDY
                     SRFAC(K,I,J,chi) = M_INV(K,chi)*WAGP(I,chi)*PHI_AREA(K,I,chi)
                     BATH(I,J,chi) = BATH(I,J,chi) + HB(K,J,1)*PHI_AREA(K,I,chi)
                     YELEM(chi) = YELEM(chi) + YDUB(K,J,chi)*PHI_AREA(K,I,chi)
                     
                     IF (ICS.EQ.1) THEN
                        SFAC_ELEM(I,J,chi)=1.0D0
                     ELSE
                        SFAC_ELEM(I,J,chi)=COS(SFEA0)/COS(YELEM(chi)/R)
                     ENDIF

                     DBATHDX(I,J,chi) = DBATHDX(I,J,chi) + HB(K,J,1)*DPHIDX
                     DBATHDY(I,J,chi) = DBATHDY(I,J,chi) + HB(K,J,1)*DPHIDY
                     DP_VOL(J,chi) = DP_VOL(J,chi) + WAGP(I,chi)*HB(K,J,1)*PHI_AREA(K,I,chi)

                  ENDDO

               enddo     

            endif

         enddo 

         do chi = 1,dg%ph
            DP_VOL(J,chi) = 0.25D0*AREAS(J)*DP_VOL(J,chi)
         enddo

         do chi = 1,dg%ph

            if (chi.ge.1) then

               DO L = 1,3
                  do I = 1,NEGP(chi) ! Edge quadrature points

                     BATHED(I,L,J,chi) = 0.D0
                     yed(chi) = 0.d0

                     DO K = 1,(chi+1)*(chi+2)/2
                        
                        BATHED(I,L,J,chi) = BATHED(I,L,J,chi)+HB(K,J,1)*PHI_EDGE(K,I,L,chi)

                        YED(chi) = YED(chi) + YDUB(K,J,chi)*PHI_EDGE(K,I,L,chi)

                        IF (ICS.EQ.1) THEN
                           SFACED(I,L,J,chi)=1.0d0
                        ELSE
                           SFACED(I,L,J,chi)=COS(SFEA0)/COS(YED(chi)/R)
                        ENDIF

                        EDGEQ(K,I,L,chi) = 2.0*M_INV(K,chi)*PHI_EDGE(K,I,L,chi)*WEGP(I,chi)

                     ENDDO
                  enddo
               ENDDO

            else

               DO L = 1,3
                  do I = 1,NEGP(chi) ! Edge quadrature points

                     BATHED(I,L,J,chi) = 0.D0
                     yed(chi) = 0.d0

                     DO K = 1,(chi+1)*(chi+2)/2
                        
                        BATHED(I,L,J,chi) = BATHED(I,L,J,chi)+HB(K,J,1)*PHI_EDGE(K,I,L,chi)

                        YED(chi) = YED(chi) + YDUB(K,J,chi)*PHI_EDGE(K,I,L,chi)

                        IF (ICS.EQ.1) THEN
                           SFACED(I,L,J,chi)=1.0d0
                        ELSE
                           SFACED(I,L,J,chi)=COS(SFEA0)/COS(YED(chi)/R)
                        ENDIF

                        EDGEQ(K,I,L,chi) = 2.0*M_INV(K,chi)*PHI_EDGE(K,I,L,chi)*WEGP(I,chi)

                     ENDDO
                  enddo
               ENDDO

            endif

         enddo

!........Store bathymetry at triangular vertices and edge gauss points for wet-dry
!........


         do chi = 1,dg%ph

            DO I = 1,3
               DP_NODE(I,J,chi) = DP(NM(J,I))
            ENDDO


            IF (NCHECK(chi).GT.3) THEN
               II = 4
               DO L = 1,3
                  DO I = 1,NEGP(chi)
                     DP_NODE(II,J,chi) = BATHED(I,L,J,chi)
                     II = II + 1
                  ENDDO
               ENDDO
            ENDIF
         enddo
      ENDDO
      
      DO I = 1,3
         IF (I.EQ.1) THEN
            XI = -1.D0
            YI = -1.D0
         ELSEIF (I.EQ.2) THEN
            XI =  1.D0 
            YI = -1.D0
         ELSE
            XI = -1.D0
            YI =  1.D0
         ENDIF
         PSI_CHECK(1,I) = -1.D0/2.D0*(XI + YI)
         PSI_CHECK(2,I) =  1.D0/2.D0*(XI + 1.D0)
         PSI_CHECK(3,I) =  1.D0/2.D0*(YI + 1.D0)
         do chi = 1,dg%ph
            DO K = 1,(chi+1)*(chi+2)/2 
               PHI_CHECK(K,I,chi) = PHI_CORNER(K,I,chi)
            ENDDO
         enddo
      ENDDO

      do chi =1,dg%ph
         IF (NCHECK(chi).GT.3) THEN
            II = 4
            DO L = 1,3
               DO I = 1,NEGP(chi)
                  DO K = 1,(chi+1)*(chi+2)/2
                     
                     PHI_CHECK(K,II,chi) = PHI_EDGE(K,I,L,chi)
                     
                  ENDDO
                  II = II + 1
               ENDDO
            ENDDO
         ENDIF
      enddo

!.....Integrate the basis functions

      PHI_INTEGRATED = 0.D0
      do chi = 1,dg%ph
         DO I = 1,NAGP(chi)
            DO K = 1,(chi+1)*(chi+2)/2
               PHI_INTEGRATED(K,chi) = PHI_INTEGRATED(K,chi) + WAGP(I,chi)*PHI_AREA(K,I,chi)
            enddo
         ENDDO
      ENDDO

!.....Wetting and drying is not turned on

      IF(NOLIFA.EQ.0.OR.NOLIFA.EQ.1) THEN
         DO J = 1,S%MNE
            WDFLG(J) = 1
                                !DOFS(J) = 3
         ENDDO
         
!.....Wetting and drying is turned on but there are no dry nodes below
!.....geoid

      ELSEIF (NOLIFA.EQ.2.AND.NSTARTDRY.EQ.0) THEN

         DO J = 1,S%MNE
            
!.........Check to see if there are initially any dry nodes

                                !ZE1 = 0.D0
                                !ZE2 = 0.D0
                                !ZE3 = 0.D0 

            ZE1 = ze(1,J,1)
            ZE2 = ze(1,J,1)
            ZE3 = ze(1,J,1)
            IF (DP(NM(J,1)).LT.H0) ZE1 = max(ze1,H0 - DP(NM(J,1)))
            IF (DP(NM(J,2)).LT.H0) ZE2 = max(ze2,H0 - DP(NM(J,2)))
            IF (DP(NM(J,3)).LT.H0) ZE3 = max(ze3,H0 - DP(NM(J,3)))
            
!.........If so set initial surface elevation values
            
            IF ((ZE1 + ZE2 + ZE3)/3.D0.NE.ze(1,J,1)) THEN
               IF (p_0.EQ.0) THEN
                  DP_MIN = MIN(DP(NM(J,1)),DP(NM(J,2)),DP(NM(J,3)))
                  ze(1,J,1) = max(ze(1,j,1),H0 - DP_MIN)
               ELSE
                  IF (ze(1,J,1).GT.(ZE1+ZE2+ZE3)/3.d0) THEN
                     IF (dof_0.NE.1) THEN
                        ze(2,J,1)=0.d0
                        ze(3,J,1)=0.d0
                        ze(4:dg%dofh,J,1) = 0.D0 ! forced again for transparency
                     ENDIF
                  ELSE
                     ze(1,J,1)=(ZE1+ZE2+ZE3)/3.D0
                     IF (DOF_0.NE.1) THEN
                        ze(2,J,1) = -1.D0/6.D0*(ZE1 + ZE2) + 1.D0/3.D0*ZE3
                        ze(3,J,1) = -0.5D0*ZE1 + 0.5D0*ZE2
                        ze(4:dg%dofh,J,1) = 0.D0 ! forced again for transparency
                     ENDIF
                  ENDIF
               ENDIF

!............and set wet/dry flag (0 = dry, 1 = wet)

               WDFLG(J) = 0
               DOFS(J) = dof_0
            ELSE
               WDFLG(J) = 1
               DOFS(J) = dofs(J)
            ENDIF

         ENDDO


!.....If there are dry nodes below geoid

      ELSEIF (NOLIFA.EQ.2.AND.NSTARTDRY.EQ.1) THEN
         
!.......Loop over elements
         
         DO J = 1,NE

!........Retrieve global node numbers for element

            N1 = NM(J,1)
            N2 = NM(J,2)
            N3 = NM(J,3)

!........Check to see if nodes are initially dry

            ZE1 = 0
            ZE2 = 0
            ZE3 = 0
            IF (STARTDRY(N1).EQ.1) ZE1 = H0 - DP(N1)
            IF (DP(N1).LT.H0) ZE1 = H0 - DP(N1)
            IF (STARTDRY(N2).EQ.1) ZE2 = H0 - DP(N2)
            IF (DP(N2).LT.H0) ZE2 = H0 - DP(N2)
            IF (STARTDRY(N3).EQ.1) ZE3 = H0 - DP(N3)
            IF (DP(N3).LT.H0) ZE3 = H0 - DP(N3)

            IF (DG%MODAL_IC.EQ.3) THEN
               IF (STARTDRY(N1).EQ.-88888) then
                  ZE1 = H0 - DP(N1)
               else
                  ZE1 = STARTDRY(N1)
               endif
               IF (STARTDRY(N2).EQ.-88888) then
                  ZE2 = H0 - DP(N2)
               else
                  ZE2 = STARTDRY(N2)
               endif
               IF (STARTDRY(N3).EQ.-88888) then
                  ZE3 = H0 - DP(N3)
               else
                  ZE3 = STARTDRY(N3)
               endif
            ENDIF

!.........If so set initial surface elevation values

            IF ((ZE1 + ZE2 + ZE3).NE.0) THEN
               IF (P_0.EQ.0) THEN
                  DP_MIN = MIN(DP(NM(J,1)),DP(NM(J,2)),DP(NM(J,3)))
                  ze(1,J,1) = max(ze(1,j,1),H0 - DP_MIN)
               ELSE
                  IF (ze(1,J,1).GT.(ZE1+ZE2+ZE3)/3.d0) THEN
                     IF (DOF_0.NE.1) THEN
                        ze(2,J,1)=0.d0
                        ze(3,J,1)=0.d0
                        ze(4:dg%dofh,J,1) = 0.D0 ! forced again for transparency
                     ENDIF
                  ELSE
                     ze(1,J,1)=(ZE1+ZE2+ZE3)/3.D0
                     IF (DOF_0.NE.1) THEN
                        ze(2,J,1) = -1.D0/6.D0*(ZE1 + ZE2) + 1.D0/3.D0*ZE3
                        ze(3,J,1) = -0.5D0*ZE1 + 0.5D0*ZE2
                        ze(4:dg%dofh,J,1) = 0.D0 ! forced again for transparency
                     ENDIF
                  ENDIF
               ENDIF

!...........Set wet/dry flag (0 = dry, 1 = wet)

               WDFLG(J) = 0
               DOFS(J) = DOF_0
            ELSE
               WDFLG(J) = 1
               DOFS(J) = DOFS(J)
            ENDIF
         ENDDO
      ENDIF

      IF(MYPROC_HERE.EQ.0) THEN
         print *, 'DONE'
         print *, ''
      ENDIF

!.....Read in modal dof for initial conditions
!Asserts error if you project onto lower order basis
!That is, do not expect convergence

      IF (DG%MODAL_IC.EQ.1) THEN
         OPEN(163,FILE=S%DIRNAME//'/'//'Initial_Conditions.163')
         OPEN(164,FILE=S%DIRNAME//'/'//'Initial_Conditions.164')
         OPEN(114,FILE=S%DIRNAME//'/'//'Initial_Bathymetry.114')
         READ(163,*) P_READ
         READ(164,*) P_READ,P_READ
         READ(114,*) P_READ
         IF (P_READ.NE.dg%ph) THEN
            PRINT*,'INCONSISTENCY IN P -- CHECK INPUT FILES'
            STOP
         ENDIF
         DO J = 1,S%MNE
            DO K = 1,DG%DOFH
               READ(163,*) ze(K,J,1)
               READ(164,*) QX(K,J,1), QY(K,J,1)
               READ(114,*) HB(K,J,1)
            ENDDO
         ENDDO
         H_TRI = SQRT((X(1)-X(2))**2 + (Y(1)-Y(2))**2)
         CLOSE(163)
         CLOSE(164)
         CLOSE(114)
      ENDIF

!.....Read in modal dof for hot start conditions

      IF (DG%MODAL_IC.EQ.2) THEN
         OPEN(263,FILE=S%DIRNAME//'/'//'Hot_start.263')
         OPEN(264,FILE=S%DIRNAME//'/'//'Hot_start.264')
         OPEN(214,FILE=S%DIRNAME//'/'//'Hot_start.214')
#ifdef TRACE
         OPEN(288,FILE=S%DIRNAME//'/'//'Hot_start.288')
#endif
#ifdef CHEM
         OPEN(289,FILE=S%DIRNAME//'/'//'Hot_start.289')
#endif
#ifdef DYNP
         OPEN(291,FILE=S%DIRNAME//'/'//'Hot_start.291')
#endif
!$$$#ifdef SED_LAY
!$$$         OPEN(290,FILE=S%DIRNAME//'/'//'Hot_start.290')
!$$$#endif
         READ(263,*) P_READ
         READ(264,*) P_READ,P_READ
         READ(214,*) ITHS
         IF (P_READ.NE.DG%PH) THEN
            PRINT*,'INCONSISTENCY IN P -- CHECK INPUT FILES'
            STOP
         ENDIF
         DO J = 1,S%MNE
            DO K = 1,DG%DOFH
               READ(263,*) ze(K,J,1)
               READ(264,*) QX(K,J,1), QY(K,J,1)
               READ(214,*) HB(K,J,1), WDFLG(J)
#ifdef TRACE
               READ(288,*) iota(K,J,1)
#endif
#ifdef CHEM
               READ(289,*) iota(K,J,1),iota2(K,J,1)
#endif
#ifdef DYNP
               READ(291,*) dynP(K,J,1)
#endif
!$$$#ifdef SED_LAY
!$$$               do ll=1,layers
!$$$                  READ(290,*) bed(K,J,1,ll)
!$$$               enddo
!$$$#endif
            ENDDO
         ENDDO
         CLOSE(263)
         CLOSE(264)
         CLOSE(214)
#ifdef TRACE
         CLOSE(288)
#endif
#ifdef CHEM
         CLOSE(289)
#endif
!$$$#ifdef SED_LAY
!$$$         CLOSE(290)
!$$$#endif
#ifdef DYNP
         CLOSE(291)
#endif
      ENDIF
      
!.....Initialize the DG.63 output file

      IF (ABS(NOUTGE).EQ.1) THEN
         OPEN(631,FILE=S%DIRNAME//'/'//'DG.63')
         WRITE(631,3220) RUNDES, RUNID, AGRID
         WRITE(631,3645) NDSETSE, dg%dofh, DTDP*NSPOOLGE, NSPOOLGE, 1
      ENDIF

!.....Initialize the DG.64 output file

      IF (ABS(NOUTGV).EQ.1) THEN
         OPEN(641,FILE=S%DIRNAME//'/'//'DG.64')
         WRITE(641,3220) RUNDES, RUNID, AGRID
         WRITE(641,3645) NDSETSV, dg%dofh, DTDP*NSPOOLGV, NSPOOLGV, 2
      ENDIF

!.....Initialize the DG.65 output file (contains elemental statuses such
!.....as the wet/dry status.

      IF ((ABS(NOUTGE).EQ.1).AND.(NOLIFA.GE.2)) THEN
         OPEN(651,FILE=S%DIRNAME//'/'//'DG.65')
         WRITE(651,3220) RUNDES, RUNID, AGRID
         WRITE(651,3645) NDSETSE, dg%dofh, DTDP*NSPOOLGE, NSPOOLGE, 1
      ENDIF
 3220 FORMAT(1X,A32,2X,A24,2X,A24)
 3645 FORMAT(1X,I10,1X,I10,1X,E15.7,1X,I5,1X,I5)

!.....Set p back to original value if p = 0

      IF (P_0.NE.dg%pl) THEN
         PDG_EL(:) = 0
         PDG(:) = 0
         DG%DOF = 1
         DG%DOFL = 1
         DOFS(:) = 1
         DOF_0 = 1
         dg%pl = 0
                                !p = 0
         pa = 0
      ENDIF
      
!.....Compute basis functions at stations

      IF (NSTAE.GT.0) THEN      ! Elevation stations
         CALL ALLOC_STAE( NSTAE )
         DO I = 1,NSTAE
            CALL STA_BASIS( XEL(I), YEL(I),  NNE(I), PHI_STAE(:,I) )
         ENDDO
      ENDIF
      
      IF (NSTAV.GT.0) THEN      ! Velocity Stations
         CALL ALLOC_STAV( NSTAV )
         DO I = 1,NSTAV
            CALL STA_BASIS( XEV(I), YEV(I),  NNV(I), PHI_STAV(:,I) )
         ENDDO
      ENDIF

!.....Prep the slopelimiter

      IF (DG%SLOPEFLAG.NE.0) THEN
         IF(MYPROC_HERE.EQ.0)THEN
            print *, 'Slope limiting prep begins, "kshanti"'
         ENDIF
         CALL ALLOC_SLOPELIM(s)
         CALL PREP_SLOPELIM(s)
         IF(MYPROC_HERE.EQ.0)THEN
            print *, 'Finished'
         ENDIF
      ENDIF
!--   

!.....Close files

      CLOSE(115)
      CLOSE(25)

      RETURN
      END SUBROUTINE PREP_DG

!***********************************************************************
!     
!     SUBROUTINE RK_TIME()
!     
!     This subroutine does preparatory stuff for DG
!     
!     Written by Ethan Kubatko (03-07-2005)
!     
!***********************************************************************

      SUBROUTINE RK_TIME()

      USE GLOBAL
      USE DG

      IMPLICIT NONE

      INTEGER L,i,j,k
      REAL(SZ) ARK, BRK, CASUM, MAX_BOA
      Real(SZ) eps_const,RKC_omega0,RKC_omega1

!.....Allocate the time stepping arrays

      NRK = DG%RK_STAGE
      CALL ALLOC_RK()

#ifdef RKSSP
!.....The forward Euler method

      IF ((DG%RK_STAGE.EQ.1).AND.(DG%RK_ORDER.EQ.1)) THEN

         ATVD(:,:) = 0.D0
         BTVD(:,:) = 0.D0
         CTVD(:,:) = 0.D0
         DTVD(:)   = 0.D0

         ATVD(1,1) = 1.D0
         BTVD(1,1) = 1.D0

!.....SSP(s,2) schemes

      ELSEIF (DG%RK_ORDER.EQ.2) THEN

         ATVD(:,:) = 0.D0
         BTVD(:,:) = 0.D0
         CTVD(:,:) = 0.D0
         DTVD(:)   = 0.D0

         DO I = 1,NRK
            DO J = 0,NRK-1

               IF ((J.EQ.(I-1)).AND.(I.LT.NRK)) THEN
                  ATVD(I,J+1) = 1.D0
                  BTVD(I,J+1) = 1.D0/(NRK-1)
               ELSEIF ((J.EQ.0).AND.(I.EQ.NRK)) THEN
                  ATVD(I,J+1) = 1.D0/NRK
               ELSEIF ((J.EQ.(NRK-1)).AND.(I.EQ.NRK)) THEN
                  ATVD(I,J+1) = (NRK-1.D0)/NRK
                  BTVD(I,J+1) = 1.D0/NRK
               ENDIF

            ENDDO
         ENDDO

!.....SSP(3,3) scheme

      ELSEIF ((DG%RK_STAGE.EQ.3).AND.(DG%RK_ORDER.EQ.3)) THEN

         ATVD(:,:) = 0.D0
         BTVD(:,:) = 0.D0
         CTVD(:,:) = 0.D0
         DTVD(:)   = 0.D0

         ATVD(1,1) = 1.D0
         ATVD(2,1) = 3.D0/4.D0
         ATVD(2,2) = 1.D0/4.D0
         ATVD(3,1) = 1.D0/3.D0
         ATVD(3,3) = 2.D0/3.D0

         BTVD(1,1) = 1.D0
         BTVD(2,2) = 1.D0/4.D0
         BTVD(3,3) = 2.D0/3.D0

!.....SSP(4,3) scheme

      ELSEIF ((DG%RK_STAGE.EQ.4).AND.(DG%RK_ORDER.EQ.3)) THEN

         ATVD(:,:) = 0.D0
         BTVD(:,:) = 0.D0
         CTVD(:,:) = 0.D0
         DTVD(:)   = 0.D0

         ATVD(1,1) = 1.D0
         ATVD(2,2) = 1.D0
         ATVD(3,1) = 2.D0/3.D0
         ATVD(3,3) = 1.D0/3.D0
         ATVD(4,4) = 1.D0

         BTVD(1,1) = 1.D0/2.D0
         BTVD(2,2) = 1.D0/2.D0
         BTVD(3,3) = 1.D0/6.D0
         BTVD(4,4) = 1.D0/2.D0

!.....SSP(5,3) scheme

      ELSEIF ((DG%RK_STAGE.EQ.5).AND.(DG%RK_ORDER.EQ.3)) THEN

         ATVD(:,:) = 0.D0
         BTVD(:,:) = 0.D0
         CTVD(:,:) = 0.D0
         DTVD(:)   = 0.D0

         ATVD(1,1) = 1.D0
         ATVD(2,2) = 1.D0
         ATVD(3,1) = 0.355909775063327D0
         ATVD(3,3) = 0.644090224936674D0
         ATVD(4,1) = 0.367933791638137D0
         ATVD(4,4) = 0.632066208361863D0
         ATVD(5,3) = 0.237593836598569D0
         ATVD(5,5) = 0.762406163401431D0

         BTVD(1,1) = 0.377268915331368D0
         BTVD(2,2) = 0.377268915331368D0
         BTVD(3,3) = 0.242995220537396D0
         BTVD(4,4) = 0.238458932846290D0
         BTVD(5,5) = 0.287632146308408D0

!.....SSP(6,3) scheme

      ELSEIF ((DG%RK_STAGE.EQ.6).AND.(DG%RK_ORDER.EQ.3)) THEN

         ATVD(:,:) = 0.D0
         BTVD(:,:) = 0.D0
         CTVD(:,:) = 0.D0
         DTVD(:)   = 0.D0

         ATVD(1,1) = 1.D0
         ATVD(2,2) = 1.D0
         ATVD(3,3) = 1.D0
         ATVD(4,1) = 0.476769811285196D0
         ATVD(4,2) = 0.098511733286064D0
         ATVD(4,4) = 0.424718455428740D0
         ATVD(5,5) = 1.D0
         ATVD(6,3) = 0.155221702560091D0
         ATVD(6,6) = 0.844778297439909D0

         BTVD(1,1) = 0.284220721334261D0
         BTVD(2,2) = 0.284220721334261D0
         BTVD(3,3) = 0.284220721334261D0
         BTVD(4,4) = 0.120713785765930D0
         BTVD(5,5) = 0.284220721334261D0
         BTVD(6,6) = 0.240103497065900D0

!.....SSP(7,3) scheme

      ELSEIF ((DG%RK_STAGE.EQ.7).AND.(DG%RK_ORDER.EQ.3)) THEN

         ATVD(:,:) = 0.D0
         BTVD(:,:) = 0.D0
         CTVD(:,:) = 0.D0
         DTVD(:)   = 0.D0

         ATVD(1,1) = 1.D0
         ATVD(2,2) = 1.D0
         ATVD(3,3) = 1.D0
         ATVD(4,1) = 0.184962588071072D0
         ATVD(4,4) = 0.815037411928928D0
         ATVD(5,1) = 0.180718656570380D0
         ATVD(5,2) = 0.314831034403793D0
         ATVD(5,5) = 0.504450309025826D0
         ATVD(6,6) = 1.D0
         ATVD(7,4) = 0.120199000000000D0
         ATVD(7,7) = 0.879801000000000D0

         BTVD(1,1) = 0.233213863663009D0
         BTVD(2,2) = 0.233213863663009D0
         BTVD(3,3) = 0.233213863663009D0
         BTVD(4,4) = 0.190078023865845D0
         BTVD(5,5) = 0.117644805593912D0
         BTVD(6,6) = 0.233213863663009D0
         BTVD(7,7) = 0.205181790464579D0

!.....SSP(8,3) scheme

      ELSEIF ((DG%RK_STAGE.EQ.8).AND.(DG%RK_ORDER.EQ.3)) THEN

         ATVD(:,:) = 0.D0
         BTVD(:,:) = 0.D0
         CTVD(:,:) = 0.D0
         DTVD(:)   = 0.D0

         ATVD(1,1) = 1.D0
         ATVD(2,2) = 1.D0
         ATVD(3,3) = 1.D0
         ATVD(4,4) = 1.D0
         ATVD(5,1) = 0.421366967085359D0
         ATVD(5,2) = 0.005949401107575D0
         ATVD(5,5) = 0.572683631807067D0
         ATVD(6,2) = 0.004254010666365D0
         ATVD(6,6) = 0.995745989333635D0
         ATVD(7,3) = 0.104380143093325D0
         ATVD(7,4) = 0.243265240906726D0
         ATVD(7,7) = 0.652354615999950D0
         ATVD(8,8) = 1.D0

         BTVD(1,1) = 0.195804015330143D0
         BTVD(2,2) = 0.195804015330143D0
         BTVD(3,3) = 0.195804015330143D0
         BTVD(4,4) = 0.195804015330143D0
         BTVD(5,5) = 0.112133754621673D0
         BTVD(6,6) = 0.194971062960412D0
         BTVD(7,7) = 0.127733653231944D0
         BTVD(8,8) = 0.195804015330143D0

!.....SSP(5,4) scheme

      ELSEIF ((DG%RK_STAGE.EQ.5).AND.(DG%RK_ORDER.EQ.4)) THEN

         ATVD(:,:) = 0.D0
         BTVD(:,:) = 0.D0
         CTVD(:,:) = 0.D0
         DTVD(:)   = 0.D0

         ATVD(1,1) = 1.D0
         ATVD(2,1) = 0.44437049406734D0
         ATVD(2,2) = 0.55562950593266D0
         ATVD(3,1) = 0.62010185138540D0
         ATVD(3,3) = 0.37989814861460D0
         ATVD(4,1) = 0.17807995410773D0
         ATVD(4,4) = 0.82192004589227D0
         ATVD(5,1) = 0.00683325884039D0
         ATVD(5,3) = 0.51723167208978D0
         ATVD(5,4) = 0.12759831133288D0
         ATVD(5,5) = 0.34833675773694D0

         BTVD(1,1) = 0.39175222700392D0
         BTVD(2,2) = 0.36841059262959D0
         BTVD(3,3) = 0.25189177424738D0
         BTVD(4,4) = 0.54497475021237D0
         BTVD(5,4) = 0.08460416338212D0
         BTVD(5,5) = 0.22600748319395D0

!.....SSP(6,4) scheme

      ELSEIF ((DG%RK_STAGE.EQ.6).AND.(DG%RK_ORDER.EQ.4)) THEN

         ATVD(:,:) = 0.D0
         BTVD(:,:) = 0.D0
         CTVD(:,:) = 0.D0
         DTVD(:)   = 0.D0

         ATVD(1,1) = 1.00000000000000D0
         ATVD(2,1) = 0.30948026455053D0
         ATVD(2,2) = 0.69051973544947D0
         ATVD(3,1) = 0.54205244285557D0
         ATVD(3,3) = 0.45794755714443D0
         ATVD(4,1) = 0.35984960863377D0
         ATVD(4,4) = 0.64015039136623D0
         ATVD(5,5) = 1.00000000000000D0
         ATVD(6,1) = 0.05776282890116D0
         ATVD(6,3) = 0.44216432622405D0
         ATVD(6,5) = 0.10115567086469D0
         ATVD(6,6) = 0.39891717401009D0

         BTVD(1,1) = 0.39270746575722D0
         BTVD(2,2) = 0.30154043149172D0
         BTVD(3,3) = 0.19997937335132D0
         BTVD(4,4) = 0.27954483459696D0
         BTVD(5,5) = 0.43668618869443D0
         BTVD(6,3) = 0.09150931531680D0
         BTVD(6,5) = 0.04417328437472D0
         BTVD(6,6) = 0.14911300530736D0

!.....SSP(7,4) scheme

      ELSEIF ((DG%RK_STAGE.EQ.7).AND.(DG%RK_ORDER.EQ.4)) THEN

         ATVD(:,:) = 0.D0
         BTVD(:,:) = 0.D0
         CTVD(:,:) = 0.D0
         DTVD(:)   = 0.D0

         ATVD(1,1) = 1.D0
         ATVD(2,1) = 0.20161507213829D0
         ATVD(2,2) = 0.79838492786171D0
         ATVD(3,1) = 0.19469598207921D0
         ATVD(3,3) = 0.80530401792079D0
         ATVD(4,1) = 0.58143386885601D0
         ATVD(4,4) = 0.41856613114399D0
         ATVD(5,1) = 0.01934367892154D0
         ATVD(5,5) = 0.98065632107846D0
         ATVD(6,6) = 1.D0
         ATVD(7,1) = 0.06006304558847D0
         ATVD(7,3) = 0.30152730794242D0
         ATVD(7,4) = 0.10518998496676D0
         ATVD(7,5) = 0.01483791154585D0
         ATVD(7,7) = 0.51838174995650D0

         BTVD(1,1) = 0.30111872706068D0
         BTVD(2,2) = 0.24040865318216D0
         BTVD(3,3) = 0.24249212077315D0
         BTVD(4,4) = 0.12603810060080D0
         BTVD(5,5) = 0.29529398308716D0
         BTVD(6,6) = 0.30111872706068D0
         BTVD(7,3) = 0.09079551914158D0
         BTVD(7,4) = 0.02888359354880D0
         BTVD(7,7) = 0.15609445267839D0

!.....SSP(8,4) scheme

      ELSEIF ((DG%RK_STAGE.EQ.8).AND.(DG%RK_ORDER.EQ.4)) THEN

         ATVD(:,:) = 0.D0
         BTVD(:,:) = 0.D0
         CTVD(:,:) = 0.D0
         DTVD(:)   = 0.D0

         ATVD(1,1) = 1.D0
         ATVD(2,1) = 0.10645325745007D0
         ATVD(2,2) = 0.89354674254993D0
         ATVD(3,3) = 1.D0
         ATVD(4,1) = 0.57175518477257D0
         ATVD(4,4) = 0.42824481522743D0
         ATVD(5,1) = 0.19161667219044D0
         ATVD(5,5) = 0.80838332780956D0
         ATVD(6,6) = 1.D0
         ATVD(7,7) = 1.D0
         ATVD(8,1) = 0.02580435327923D0
         ATVD(8,3) = 0.03629901341774D0
         ATVD(8,4) = 0.31859181340256D0
         ATVD(8,5) = 0.05186768980103D0
         ATVD(8,6) = 0.03944076217320D0
         ATVD(8,7) = 0.00511633747411D0
         ATVD(8,8) = 0.52288003045213D0

         BTVD(1,1) = 0.24120020561311D0
         BTVD(2,2) = 0.21552365802797D0
         BTVD(3,3) = 0.24120020561311D0
         BTVD(4,4) = 0.10329273748560D0
         BTVD(5,5) = 0.19498222488188D0
         BTVD(6,6) = 0.24120020561311D0
         BTVD(7,7) = 0.24120020561311D0
         BTVD(8,3) = 0.00875532949991D0
         BTVD(8,4) = 0.06195575835101D0
         BTVD(8,6) = 0.00951311994571D0
         BTVD(8,8) = 0.12611877085604D0

      ENDIF

!.....Compute the time dependent parameters

      DO K = 0,DG%RK_STAGE-1
         DO I = 1,DG%RK_STAGE
            CASUM = 0.D0
            DO L = K+1,I-1
               CASUM = CASUM + CTVD(L,K+1)*ATVD(I,L+1)
            ENDDO
            CTVD(I,K+1) = BTVD(I,K+1) + CASUM
         ENDDO
      ENDDO

      DO K = 1,DG%RK_STAGE-1
         DTVD(K+1) = 0.D0
         DO L = 0,K-1
            DTVD(K+1) = DTVD(K+1) + CTVD(K,L+1)
         ENDDO
      ENDDO

!.....Compute the maximum beta over alpha ratio at each stage

      DO IRK = 1,DG%RK_STAGE
         MAX_BOA = 0.D0
         DO I = 1,IRK
            ARK = ATVD(IRK,I)
            BRK = BTVD(IRK,I)
            IF (ARK.NE.0.D0) THEN
               IF (MAX_BOA.LT.BRK/ARK) MAX_BOA = BRK/ARK
            ENDIF
         ENDDO
         MAX_BOA_DT(IRK) = MAX_BOA*DT
      ENDDO

#endif

!-----------------------------------------------------------
!.... Compute the Runge-Kutta Chebyshev (RKC) version

#ifdef RKC
         
         RKC_omega0 = 0.D0
         RKC_omega1 = 0.D0

         eps_const = 0.153846154D0 

                                !Compute the RKC version

         RKC_T(0) = 1.D0
         RKC_U(0) = 1.D0
         
         RKC_omega0 = 1.D0 + eps_const*nrk**(-2.D0)

         RKC_c(0) = 0.D0
         RKC_b(0) = 0.D0
         RKC_b(1) = 1.D0/RKC_omega0

         RKC_T(1) = RKC_omega0
         RKC_U(1) = 2.D0*RKC_omega0

         RKC_Tprime(1) = nrk**(2.D0)
         RKC_Tdprime(1) = 1.D0/3.D0 * (nrk**(2.D0)*(nrk**(2.D0) - 1.D0))


         do j = 0,nrk

                                !Chebyshev polynomial stuff

            if(j.gt.1) then

               RKC_T(j) = 2.0 * RKC_omega0 * RKC_T(j-1) - RKC_T(j-2)
               RKC_U(j) = 2.0 * RKC_omega0 * RKC_U(j-1) - RKC_U(j-2)
               
            endif

            if(j.gt.1) then

               RKC_Tprime(j) = j*RKC_U(j-1)
               RKC_Tdprime(j) = (j*(j*RKC_T(j) - RKC_omega0*RKC_U(j-1))/( RKC_omega0**(2.D0) - 1.D0 ) )
               RKC_b(j) = RKC_Tdprime(j)*RKC_Tprime(j)**(-2.D0)

            endif
            
            RKC_a(j)  = 1.D0 - RKC_b(j)*RKC_T(j)
            
         enddo

         if(Dg%rk_stage.gt.1) then

            RKC_b(0) = RKC_b(2)

         endif



         RKC_omega1 =  RKC_Tprime(nrk)/RKC_Tdprime(nrk)

         RKC_tildemu(1)= nrk**(-1.D0) * RKC_omega1
         
         do j=2,nrk

            RKC_mu(j) = 2.D0*RKC_b(j)*RKC_omega0/RKC_b(j-1)
            
            RKC_tildemu(j) = 2.D0*RKC_b(j)*RKC_omega1/RKC_b(j-1)
            
            RKC_nu(j) = - RKC_b(j)/RKC_b(j-2)
            
            RKC_gamma(j) = -RKC_a(j-1)*RKC_tildemu(j)
            
            RKC_c(j) = (j**(2.D0)-1.0) / ( nrk**(2.D0) - 1.D0 )

         enddo

         RKC_c(1) = (RKC_c(2)) / (4.D0*RKC_omega0)

         RKC_c(Dg%rk_stage) = 0.D0
         RKC_c(Dg%rk_stage) = 1.D0

#endif


      DEALLOCATE(CTVD)
      RETURN
      END SUBROUTINE RK_TIME
