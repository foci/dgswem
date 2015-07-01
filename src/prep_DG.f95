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

      SUBROUTINE PREP_DG(s,dg_here,global_here)

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
      type (dg_type) :: dg_here
      type (global_type) :: global_here

!.....Declare local variables

      INTEGER II, l, P_0, DOF_0,j,k,kk,jj,i,chi,ll,mm,irk
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
      Allocate ( YELEM(dg_here%ph),YED(dg_here%ph),hb1(36,s%mne,1,dg_here%ph), zeo(36,s%mne,1) )


      dg_here%C13 = 1.D0/3.D0
      dg_here%C16 = 1.D0/6.D0
      R = 6378206.4d0

#ifdef CMPI
      myproc_here = myproc
#else
      myproc_here = 0
#endif


!     sb-PDG1 moved from other places

!.....Obtain RK time scheme parameters

      CALL RK_TIME(dg_here,global_here)
      
!.....Compute the degrees of freedom global_here%per element

      dg_here%DOF = (dg_here%pl+1)*(dg_here%pl+2)/2
      dg_here%dofx = (dg_here%px+1)*(dg_here%px+2)/2    ! dg_here%dofx for variable functions dg_here%f=dg_here%f(global_here%x) 
      P_0 = dg_here%pl
      DOF_0 = (dg_here%pl+1)*(dg_here%pl+2)/2   ! dg_here%dof at lowest order when p!=0
      dg_here%dofh = (dg_here%ph + 1)*(dg_here%ph + 2)/2


!.....Allocate some DG stuff

      IF (dg_here%PADAPT.EQ.1) THEN

         dg_here%dofh = (dg_here%ph + 1)*(dg_here%ph + 2)/2
         dg_here%dofl = (dg_here%pl + 1)*(dg_here%pl + 2)/2
         dg_here%pa = dg_here%pl
           
      elseif (dg_here%padapt.eq.0) then 

         dg_here%dofh = dg_here%dofh
         dg_here%dofl = DOF_0
         dg_here%pa = dg_here%pl

      endif

#ifdef SED_LAY
!.....Initialize funtion parser for sediment types
      init_parser = .false.
#endif

#ifdef SED_LAY
!.....Initialize stabilizer sweep for A.D.
      dg_here%stblzr = .false.
#endif


!.....Compute the number of gauss points needed for the global_here%edge integrals

      CALL ALLOC_DG4(s,dg_here)          !moved here 6.28.10, for p_adapt because of messenger_elem      

         dg_here%dofs(:) = dg_here%dofl
         global_here%PDG_EL(:) = dg_here%pl
         dg_here%PDG(:) = dg_here%pl  
         dg_here%PCOUNT(:) = 0
         dg_here%pa = dg_here%pl

      do chi=dg_here%pl,dg_here%ph
         dg_here%NEGP(chi) = chi + 1
      enddo

      IF (dg_here%pl.eq.0) THEN

         global_here%PDG_EL(:) = 1
         dg_here%PDG(:) = 1     
         dg_here%DOF    = 3
         dg_here%pl     = 1
         dg_here%dofl   = 3
         P_0    = 0
         DOF_0  = 1 
         dg_here%NEGP(dg_here%pl) = 2

      ENDIF

!     cnd
!     dg_here%iwrite=0
      
!.....Initilization for parallel DG run

#ifdef CMPI

      CALL MSG_TYPES_ELEM()     ! Determine Word Sizes for Message-Passing
      CALL MSG_TABLE_ELEM()     ! Read Message-Passing Tables

      IF (dg_here%SLOPEFLAG.ge.4) THEN
         CALL MSG_TYPES()
         CALL MSG_TABLE()
      ENDIF

#endif

!.....Create the global_here%edge based data

      IF(MYPROC_HERE.EQ.0) THEN
         PRINT*, 'CREATING global_here%EDGE DATA...'
         PRINT*, ''
      ENDIF
      CALL CREATE_EDGE_DATA(s,dg_here,global_here)
      IF(MYPROC_HERE.EQ.0) THEN
         print *, 'CREATING global_here%EDGE DATA DONE'
         print *, ''
      ENDIF

#ifdef CMPI
      CALL MESSAGE_START_ELEM() ! Startup persistent message passing
      IF (dg_here%SLOPEFLAG.ge.4) CALL MESSAGE_START()
#endif

!.....Re-arrange elevation specified boundary segment data for DG

      IF (dg_here%NEEDS.GT.0) THEN
         CALL ALLOC_DG1(dg_here,s%MNBFR)
         II = 1
         JJ = 1
         DO I = 1,global_here%NBFR
            DO J = 1,global_here%NOPE
               DO K = 1,global_here%NVDLL(J)-1
                  dg_here%EMO_DG(I,II,1) = global_here%EMO(I,JJ)
                  dg_here%EMO_DG(I,II,2) = global_here%EMO(I,JJ+1)
                  dg_here%EFA_DG(I,II,1) = global_here%EFA(I,JJ)
                  dg_here%EFA_DG(I,II,2) = global_here%EFA(I,JJ+1)
                  dg_here%UMO_DG(I,II,1) = global_here%UMO(I,JJ)
                  dg_here%UMO_DG(I,II,2) = global_here%UMO(I,JJ+1)
                  dg_here%UFA_DG(I,II,1) = global_here%UFA(I,JJ)
                  dg_here%UFA_DG(I,II,2) = global_here%UFA(I,JJ+1)
                  dg_here%VMO_DG(I,II,1) = global_here%VMO(I,JJ)
                  dg_here%VMO_DG(I,II,2) = global_here%VMO(I,JJ+1)
                  dg_here%VFA_DG(I,II,1) = global_here%VFA(I,JJ)
                  dg_here%VFA_DG(I,II,2) = global_here%VFA(I,JJ+1)
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

      IF (dg_here%NFEDS.GT.0) THEN
         CALL ALLOC_DG2(dg_here,s%MNFFR)
         II = 1
         JJ = 1
         DO I = 1,global_here%NFFR
            DO J = 1,global_here%NBOU
               IF ( (global_here%SEGTYPE(J).EQ.2 ).OR.(global_here%SEGTYPE(J).EQ.12)&
              .OR.(global_here%SEGTYPE(J).EQ.22) ) THEN
                  DO K = 1,global_here%NVELL(J)-1
                     dg_here%QNAM_DG(I,II,1) = global_here%QNAM(I,JJ)
                     dg_here%QNAM_DG(I,II,2) = global_here%QNAM(I,JJ+1)
                     dg_here%QNPH_DG(I,II,1) = global_here%QNPH(I,JJ)
                     dg_here%QNPH_DG(I,II,2) = global_here%QNPH(I,JJ+1)
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

      IF (dg_here%NIBEDS.ne.0) CALL ALLOC_DG3(dg_here,S%MNP)

!.....Allocate the array for node to element table

      CALL ALLOC_NNOEL1(S,global_here)

!.....Determine the number of elements connected at each node

      global_here%EL_COUNT = 0
      global_here%MAXEL = 1
      DO K = 1,3
         DO J = 1,S%MNE
            global_here%N1 = global_here%NM(J,K)
            global_here%EL_COUNT(global_here%N1) = global_here%EL_COUNT(global_here%N1) + 1
         ENDDO
      ENDDO
      global_here%MAXEL = MAXVAL(global_here%EL_COUNT)

!.....Allocate the array for the node to element table

      CALL ALLOC_NNOEL2(S,global_here,global_here%MAXEL)
      
!.....Construct node to element table

      global_here%EL_COUNT = 0
      DO K = 1,3
         DO J = 1,S%MNE
            global_here%N1 = global_here%NM(J,K)
            global_here%NNOEL(global_here%N1,1+global_here%EL_COUNT(global_here%N1)) = J
            global_here%EL_COUNT(global_here%N1) = global_here%EL_COUNT(global_here%N1) + 1
         ENDDO
      ENDDO

!.....Construct node to element angle table

      DO I = 1,S%MNP
         global_here%ETAMAX(I) = -99999
         KK = 1
         global_here%ELETAB(I,1) = I
         dg_here%S1 = global_here%SFAC(I)
         dg_here%J1 = global_here%NEITAB(I,1)
         DO 111 K = 1,global_here%NNEIGH(I)-1
            dg_here%J2 = global_here%NEITAB(I,1+K)
            IF (K.LT.(global_here%NNEIGH(I)-1)) THEN
               dg_here%J3 = global_here%NEITAB(I,2+K)
            ELSE
               dg_here%J3 = global_here%NEITAB(I,2)
            ENDIF
            DO J = 1,global_here%EL_COUNT(I)
               dg_here%EL = global_here%NNOEL(I,J)
               global_here%N1 = global_here%NM(dg_here%EL,1)
               global_here%N2 = global_here%NM(dg_here%EL,2)
               global_here%N3 = global_here%NM(dg_here%EL,3)
               IF ((dg_here%J1.EQ.global_here%N1).OR.(dg_here%J1.EQ.global_here%N2).OR.(dg_here%J1.EQ.global_here%N3)) THEN
                  IF ((dg_here%J2.EQ.global_here%N1).OR.(dg_here%J2.EQ.global_here%N2).OR.(dg_here%J2.EQ.global_here%N3)) THEN
                     IF ((dg_here%J3.EQ.global_here%N1).OR.(dg_here%J3.EQ.global_here%N2).OR.(dg_here%J3.EQ.global_here%N3)) THEN
                        global_here%ELETAB(I,1+KK) = dg_here%EL
                        dg_here%S2  = global_here%SFAC(dg_here%J2)
                        dg_here%SAV = (dg_here%S1 + dg_here%S2)/2.D0
                        dg_here%VEC1(1) =      global_here%X(dg_here%J1) - global_here%X(dg_here%J2)
                        dg_here%VEC1(2) = dg_here%SAV*(global_here%Y(dg_here%J1) - global_here%Y(dg_here%J2))
                        dg_here%VEC2(1) =      global_here%X(dg_here%J1) - global_here%X(dg_here%J3)
                        dg_here%VEC2(2) = dg_here%SAV*(global_here%Y(dg_here%J1) - global_here%Y(dg_here%J3))
                        dg_here%MAG1 = SQRT(dg_here%VEC1(1)**2 + dg_here%VEC1(2)**2)
                        dg_here%MAG2 = SQRT(dg_here%VEC2(1)**2 + dg_here%VEC2(2)**2)
                        dg_here%DOT = DOT_PRODUCT(dg_here%VEC1,dg_here%VEC2)
                        dg_here%EL_ANG  = ACOS(dg_here%DOT/(dg_here%MAG1*dg_here%MAG2))
                        global_here%ANGTAB(I,KK+1) = RAD2DEG*dg_here%EL_ANG
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

      dg_here%ZE = 0.D0
      zeo = 0.D0
      dg_here%QX = 0.D0
      dg_here%QY = 0.D0
      hbo = 0.D0
      dg_here%hb = 0.D0
      dg_here%LZ = 0.d0
      dg_here%HZ = 0.D0
      dg_here%TZ = 0.D0
      dg_here%MZ = 0.D0
      dg_here%iota = 0.D0
      dg_here%iota2 = 0.D0
      dg_here%iotaa = 0.D0
      dg_here%iotaa2 = 0.D0 
      dg_here%iotaa3 = 0.D0       
      global_here%MassMax = 0.D0
      dg_here%MARK = 0
#ifdef SEDLAY
      dg_here%bed = 0.D0
#endif
      dg_here%RHS_ZE = 0.D0
      dg_here%RHS_QX = 0.D0
      dg_here%RHS_QY = 0.D0
      dg_here%RHS_iota = 0.D0
      dg_here%RHS_iota2 = 0.D0
      
      global_here%WSX2(:) = 0
      global_here%WSY2(:) = 0

!.....If using modal initial conditions transform the bathymetry from
!.....nodal coordinates to modal dg_here%dof

      DO J = 1,S%MNE
         global_here%N1 = global_here%NM(J,1)
         global_here%N2 = global_here%NM(J,2)
         global_here%N3 = global_here%NM(J,3)
         hbo(1,J,1) =  1.D0/3.D0 * (global_here%DP(global_here%N1) + global_here%DP(global_here%N2) + global_here%DP(global_here%N3))
         hbo(2,J,1) = -1.D0/6.D0 * (global_here%DP(global_here%N1) + global_here%DP(global_here%N2)) + 1.D0/3.D0*global_here%DP(global_here%N3)
         hbo(3,J,1) = -0.5D0*global_here%DP(global_here%N1) + 0.5D0*global_here%DP(global_here%N2)

         ydubo(1,J)= 1.D0/3.D0*(global_here%Y(global_here%N1) + global_here%Y(global_here%N2) + global_here%Y(global_here%N3))
         ydubo(2,J) = -1.D0/6.D0*(global_here%Y(global_here%N1) + global_here%Y(global_here%N2))&
        + 1.D0/3.D0*global_here%Y(global_here%N3)
         ydubo(3,J) = -0.5D0*global_here%Y(global_here%N1) + 0.5D0*global_here%Y(global_here%N2)

      ENDDO


      if (LoadManningsN) then
         DO J = 1,global_here%NE
            global_here%N1 = global_here%NM(J,1)
            global_here%N2 = global_here%NM(J,2)
            global_here%N3 = global_here%NM(J,3)
            dg_here%MANN(1,J) =  1.D0/3.D0*(ManningsN(global_here%N1)&
           + ManningsN(global_here%N2) + ManningsN(global_here%N3))
            dg_here%MANN(2,J) = -1.D0/6.D0*(ManningsN(global_here%N1) &
           + ManningsN(global_here%N2)) + 1.D0/3.D0*ManningsN(global_here%N3)
            dg_here%MANN(3,J) = -0.5D0*ManningsN(global_here%N1) + 0.5D0*ManningsN(global_here%N2)
         ENDDO
      endif

      IF (dg_here%MODAL_IC.EQ.0) THEN
!     this assumes a cold start
         if (LoadGeoidOffset) then
            DO J = 1,global_here%NE
               global_here%N1 = global_here%NM(J,1)
               global_here%N2 = global_here%NM(J,2)
               global_here%N3 = global_here%NM(J,3)
               zeo(1,J,1)=1.d0/3.d0*(GeoidOffset(global_here%N1)+GeoidOffset(global_here%N2)+&
              GeoidOffset(global_here%N3))
               IF (dof_0.ne.1) THEN
                  zeo(2,J,1)=-1.d0/6.d0*(GeoidOffset(global_here%N1)+GeoidOffset(global_here%N2))&
                 +1.d0/3.d0*GeoidOffset(global_here%N3)
                  zeo(3,J,1)=-.5d0*GeoidOffset(global_here%N1)+.5d0*GeoidOffset(global_here%N2)
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
         print*, 'In global_here%X we have: ', global_here%sed_equationX
         print*, 'In global_here%Y we have: ', global_here%sed_equationY
         open(444, file = "./sedlaw.global_here%X")
         write(444,'(a)') global_here%sed_equationX
         close(444)
         open(445, file = "./sedlaw.global_here%Y")
         write(445,*) global_here%sed_equationY
         close(445)
         CALL SYSTEM('python py_scriptX') !this writes db_partials_X file
         CALL SYSTEM('python py_scriptY') !this writes db_partials_Y file
         print *, ''
      ENDIF
#endif

      IF(MYPROC_HERE.EQ.0) THEN
         print *, 'PREP FOR WET/DRY BEGINS...'
      ENDIF

!.....1. Set initial surface elevation above the dg_here%bed elevation
!.....if wetting-and-drying is enabled and the initial water depth is
!.....not specified by fort.12
!.....2. Set wet-and-dry elemental flags
!.....3. Set the dg_here%DOF at dry elements = 1

      dg_here%NCHECK(1) = 3
      if (dg_here%ph.gt.1) then
         do chi = 2,dg_here%ph
            dg_here%NCHECK(chi) = dg_here%NCHECK(1) + 3*dg_here%negp(chi)
         enddo
      endif
      
      CALL ALLOC_DG_WETDRY(s,dg_here)
      dg_here%PHI_CHECK = 0.D0
      dg_here%PSI_CHECK = 0.D0
      global_here%H0L = global_here%H0
      global_here%H0H = global_here%H0 * 1.0
      global_here%HABSMIN = global_here%H0 * 1.0
            
!.....Retrieve the normals to the edges

      CALL CALC_NORMAL(dg_here,global_here)

!.....Retrieve the area integral gauss quadrature points
      
      do j=1,dg_here%ph
         
         CALL QUAD_PTS_AREA(s,dg_here,global_here,2*j,j)

      enddo

!.....Retrieve the global_here%edge integral gauss quadrature points
      
      do j=1,dg_here%ph

         CALL QUAD_PTS_EDGE(s,dg_here,global_here,j,j)

      enddo

!.....Evaluate the orthogonal basis and its derivatives at the area
!.....gauss quadrature points

      do j=1,dg_here%ph

         CALL ORTHOBASIS_AREA(dg_here,global_here,j,j)

      enddo

!.....Evaluate the orthogonal basis at the global_here%edge gauss quadrature points

      do j=1,dg_here%ph

         CALL ORTHOBASIS_EDGE(dg_here,global_here,j,j)

      enddo

!.....Do the L2-projection of the initial conditions

      dg_here%hb = 0.D0
      dg_here%qx = 0.D0
      dg_here%qy = 0.D0
      dg_here%ze = 0.D0
      global_here%ydub = 0.d0
      hb1 = 0.D0

!$$$      do k = 1,S%MNE
!$$$
!$$$         
!$$$         global_here%n1 = global_here%NM(k,1)
!$$$         global_here%n2 = global_here%NM(k,2)
!$$$         global_here%n3 = global_here%NM(k,3)
!$$$                                !Define lagrange transform
!$$$
!$$$         do mm = 1,dg_here%nagp(dg_here%ph)     !global_here%ICs should not have higher order than dg_here%ph
!$$$
!$$$            ell_1 = -0.5D0 * ( dg_here%xagp(mm,dg_here%ph) + dg_here%yagp(mm,dg_here%ph) )
!$$$            ell_2 =  0.5D0 * ( dg_here%xagp(mm,dg_here%ph) + 1.D0 )
!$$$            ell_3 =  0.5D0 * ( dg_here%yagp(mm,dg_here%ph) + 1.D0 )
!$$$
!$$$            XBCbt(k) = global_here%x(global_here%n1)*ell_1 + global_here%x(global_here%n2)*ell_2 + global_here%x(global_here%n3)*ell_3
!$$$            YBCbt(k) = global_here%y(global_here%n1)*ell_1 + global_here%y(global_here%n2)*ell_2 + global_here%y(global_here%n3)*ell_3
!$$$
!$$$
!$$$            rev =  3.141592653589793D0 / 4.D0
!$$$
!$$$            Ox = XBCbt(k)*cos(rev) - YBCbt(k)*sin(rev)
!$$$            Oy = YBCbt(k)*cos(rev) + XBCbt(k)*sin(rev)
!$$$
!$$$            radial(k) = min(sqrt( Ox**2 + (Oy + 0.25D0)**2 ), 0.18D0)/0.18D0
!$$$
!$$$
!$$$            do j = 1,dg_here%dofh
!$$$
!$$$               dg_here%QX(j,k,1) = dg_here%QX(j,k,1) + YBCbt(k) * dg_here%wagp(mm,dg_here%ph) * dg_here%phi_area(j,mm,dg_here%ph)
!$$$
!$$$               dg_here%QY(j,k,1) = dg_here%QY(j,k,1) - XBCbt(k) * dg_here%wagp(mm,dg_here%ph) * dg_here%phi_area(j,mm,dg_here%ph)        
!$$$
!$$$               dg_here%iota(j,k,1) = dg_here%iota(j,k,1) + 0.25D0 *( 1.0D0 + cos(3.141592653589793D0*radial(k)) ) 
!$$$     &              * dg_here%wagp(mm,dg_here%ph) * dg_here%phi_area(j,mm,dg_here%ph)  
!$$$               
!$$$c$$$               dg_here%iota(j,k,1) = dg_here%iota(j,k,1) + 0.5D0 *( exp ( - ( (XBCbt(k)+.05D0)**2 + (YBCbt(k)+.05D0)**2  ) /0.001D0 ) ) 
!$$$c$$$     &              * dg_here%wagp(mm,dg_here%ph) * dg_here%phi_area(j,mm,dg_here%ph) 
!$$$
!$$$               if( ( sqrt((Ox + 0.25D0 )**2 + Oy**2)).le.0.18D0.and.(sqrt((Ox + 0.25D0 )**2 + 
!$$$     &              Oy**2)).ge.0.025D0.and.(Ox.le.(-0.23D0) )) then
!$$$
!$$$                  dg_here%iota(j,k,1) = dg_here%iota(j,k,1) + 1.0D0 * dg_here%wagp(mm,dg_here%ph) * dg_here%phi_area(j,mm,dg_here%ph) 
!$$$
!$$$               elseif(sqrt( (Ox -.25D0)**2 + Oy**2 ).le.0.18D0 ) then
!$$$
!$$$                  dg_here%iota(j,k,1) =  dg_here%iota(j,k,1) + (1.D0 - ( 1.D0 / 0.18D0 ) * sqrt((Ox -0.25D0)**2 + Oy**2 ) ) 
!$$$     &                 * dg_here%wagp(mm,dg_here%ph) * dg_here%phi_area(j,mm,dg_here%ph) 
!$$$
!$$$               endif
!$$$
!$$$               dg_here%hb(j,k,1) = dg_here%hb(j,k,1) + 1.D0*dg_here%wagp(mm,dg_here%ph) * dg_here%phi_area(j,mm,dg_here%ph)
!$$$
!$$$            enddo
!$$$
!$$$
!$$$         enddo
!$$$
!$$$                                ! get back the coeffs in each component
!$$$
!$$$         do j= 1,dg_here%dofh
!$$$
!$$$            dg_here%QX(j,k,1) =  dg_here%QX(j,k,1) * dg_here%M_inv(j,dg_here%ph)
!$$$            dg_here%QY(j,k,1) =  dg_here%QY(j,k,1) * dg_here%M_inv(j,dg_here%ph)
!$$$            dg_here%iota(j,k,1) =  dg_here%iota(j,k,1) * dg_here%M_inv(j,dg_here%ph)
!$$$            dg_here%hb(j,k,1) =  dg_here%hb(j,k,1) * dg_here%M_inv(j,dg_here%ph)
!$$$
!$$$         enddo
!$$$
!$$$      enddo

      !dg_here%iotaa = dg_here%iota
 
      dg_here%hb(1:dg_here%dofh,:,1) = hbo(1:dg_here%dofh,:,1)
      dg_here%ze(1:dg_here%dofh,:,1) = zeo(1:dg_here%dofh,:,1)
      global_here%ydub(1:dg_here%dofh,:,1) = ydubo(1:dg_here%dofh,:)

      do chi=dg_here%pl,dg_here%ph
         hb1(1:dg_here%dofh,:,1,chi) = hbo(1:dg_here%dofh,:,1)
      enddo

      !if layers are on, distribute them evenly across the total dg_here%bed load
      !other approaches are clearly available, this is a simple choice
      !adapt for higher order initial data

#ifdef SED_LAY

      do ll=1,layers
         
         dg_here%bed(:,:,1,ll) = dg_here%hb(:,:,1) / layers
         
      enddo

#endif      

#ifdef ARTDIF
      !Set up the artififical diffusion stuff
      dg_here%e1(:) = 0.D0
      dg_here%balance(:) = 0.D0
      global_here%entrop(:,:) = -100.D0

      if (dg_here%tune_by_hand.eq.1) then

         dg_here%balance(4) = 0.D0     
         
         dg_here%s0    =  0.0D0       
         dg_here%kappa =   -1.D0      
         
         dg_here%e1(1) = 0.D0
         dg_here%e1(2) = 0.D0
         dg_here%e1(3) = 0.D0
         dg_here%e1(4) = 2.5e-6       
         dg_here%e1(5) = 0.D0
         
      else
         
         dg_here%e1 = dg_here%uniform_dif

      endif

#endif

!Update dg_here%DPE_MIN

      DO J = 1,global_here%NE
         dg_here%DPE_MIN(J) = MIN(global_here%DP(global_here%NM(J,1)),global_here%DP(global_here%NM(J,2)),global_here%DP(global_here%NM(J,3)))
      ENDDO

!.....Compute the values of the nodal basis functions at the
!.....area gauss quadrature points, at every p level chi

      do chi = 1,dg_here%ph
         do I = 1,dg_here%NAGP(chi)
            dg_here%PSI1(I,chi) = -1.D0/2.D0*(dg_here%XAGP(I,chi) + dg_here%YAGP(I,chi))
            dg_here%PSI2(I,chi) =  1.D0/2.D0*(dg_here%XAGP(I,chi) + 1.D0)
            dg_here%PSI3(I,chi) =  1.D0/2.D0*(dg_here%YAGP(I,chi) + 1.D0)
         enddo
      enddo
      
!.....Store the derivatives of the (linear) nodal basis functions

      dg_here%DRPSI(1) = -1.D0/2.D0
      dg_here%DRPSI(2) =  1.D0/2.D0
      dg_here%DRPSI(3) =  0.D0
      dg_here%DSPSI(1) = -1.D0/2.D0
      dg_here%DSPSI(2) =  0.D0
      dg_here%DSPSI(3) =  1.D0/2.D0
      

!.....Pre-compute the derivatives of the coordinate transformation for
!.....each element

      DO J = 1,S%MNE

!.....Retrieve the global node numbers for the element

         global_here%N1 = global_here%NM(J,1)
         global_here%N2 = global_here%NM(J,2)
         global_here%N3 = global_here%NM(J,3)
         global_here%x1=global_here%x(global_here%n1)
         global_here%y1=global_here%y(global_here%n1)
         global_here%x2=global_here%x(global_here%n2)
         global_here%y2=global_here%y(global_here%n2)
         global_here%x3=global_here%x(global_here%n3)
         global_here%y3=global_here%y(global_here%n3)
         AREA = (global_here%X1 - global_here%X3)*(global_here%Y2 - global_here%Y3) + (global_here%X3 - global_here%X2)*(global_here%Y1 - global_here%Y3)
         area=area*.5d0

!.....Compute the derivatives of the coordinate transformation

         dg_here%DRDX(J) = 1.D0/AREA*(global_here%Y(global_here%N3) - global_here%Y(global_here%N1))
         dg_here%DSDX(J) = 1.D0/AREA*(global_here%Y(global_here%N1) - global_here%Y(global_here%N2))

         dg_here%DRDY(J) = 1.D0/AREA*(global_here%X(global_here%N1) - global_here%X(global_here%N3))
         dg_here%DSDY(J) = 1.D0/AREA*(global_here%X(global_here%N2) - global_here%X(global_here%N1))
         
!.......Compute elemental Coriolis and friction terms

         dg_here%CORI_EL(J) = (global_here%CORIF(global_here%N1) + global_here%CORIF(global_here%N2) + global_here%CORIF(global_here%N3))/3.D0
         dg_here%FRIC_EL(J) = (FRIC(global_here%N1) + FRIC(global_here%N2) + FRIC(global_here%N3))/3.D0

!.......Pre-compute the bathymetry and the gradient of the bathymetry at
!.......the quadrature points and compute volume of water

         dg_here%DP_VOL(J,:) = 0.D0
         dg_here%SFAC_ELEM(:,J,:) = 0.D0

         do chi = 1,dg_here%ph

            if (chi.gt.1) then

               do I = 1,dg_here%NAGP(chi) ! Forced for wetting and drying's sake

                  dg_here%BATH(I,J,chi) = 0.D0
                  dg_here%DBATHDX(I,J,chi) = 0.D0
                  dg_here%DBATHDY(I,J,chi) = 0.D0
                  YELEM(chi) = 0.d0

                  DO K = 1,(chi+1)*(chi+2)/2 

                     dg_here%DPHIDX = dg_here%DRPHI(K,I,chi)*dg_here%DRDX(J) + dg_here%DSPHI(K,I,chi)*dg_here%DSDX(J)
                     dg_here%DPHIDY = dg_here%DRPHI(K,I,chi)*dg_here%DRDY(J) + dg_here%DSPHI(K,I,chi)*dg_here%DSDY(J)
                     dg_here%XFAC(K,I,J,chi)  = dg_here%M_INV(K,chi)*dg_here%WAGP(I,chi)*dg_here%DPHIDX
                     dg_here%YFAC(K,I,J,chi)  = dg_here%M_INV(K,chi)*dg_here%WAGP(I,chi)*dg_here%DPHIDY
                     dg_here%SRFAC(K,I,J,chi) = dg_here%M_INV(K,chi)*dg_here%WAGP(I,chi)*dg_here%PHI_AREA(K,I,chi)
                     dg_here%BATH(I,J,chi) = dg_here%BATH(I,J,chi) + dg_here%HB(K,J,1)*dg_here%PHI_AREA(K,I,chi)
                     YELEM(chi) = YELEM(chi) + global_here%YDUB(K,J,chi)*dg_here%PHI_AREA(K,I,chi)
                     
                     IF (global_here%ICS.EQ.1) THEN
                        dg_here%SFAC_ELEM(I,J,chi)=1.0D0
                     ELSE
                        dg_here%SFAC_ELEM(I,J,chi)=COS(global_here%SFEA0)/COS(YELEM(chi)/R)
                     ENDIF

                     dg_here%DBATHDX(I,J,chi) = dg_here%DBATHDX(I,J,chi) + dg_here%HB(K,J,1)*dg_here%DPHIDX
                     dg_here%DBATHDY(I,J,chi) = dg_here%DBATHDY(I,J,chi) + dg_here%HB(K,J,1)*dg_here%DPHIDY
                     dg_here%DP_VOL(J,chi) = dg_here%DP_VOL(J,chi) + dg_here%WAGP(I,chi)*dg_here%HB(K,J,1)*dg_here%PHI_AREA(K,I,chi)

                  ENDDO

               enddo

            else

               do I = 1,dg_here%NAGP(chi) ! Area quadrature points

                  dg_here%BATH(I,J,chi) = 0.D0
                  dg_here%DBATHDX(I,J,chi) = 0.D0
                  dg_here%DBATHDY(I,J,chi) = 0.D0
                  YELEM(chi) = 0.d0

                  DO K = 1,(chi+1)*(chi+2)/2 

                     dg_here%DPHIDX = dg_here%DRPHI(K,I,chi)*dg_here%DRDX(J) + dg_here%DSPHI(K,I,chi)*dg_here%DSDX(J)
                     dg_here%DPHIDY = dg_here%DRPHI(K,I,chi)*dg_here%DRDY(J) + dg_here%DSPHI(K,I,chi)*dg_here%DSDY(J)
                     dg_here%XFAC(K,I,J,chi)  = dg_here%M_INV(K,chi)*dg_here%WAGP(I,chi)*dg_here%DPHIDX
                     dg_here%YFAC(K,I,J,chi)  = dg_here%M_INV(K,chi)*dg_here%WAGP(I,chi)*dg_here%DPHIDY
                     dg_here%SRFAC(K,I,J,chi) = dg_here%M_INV(K,chi)*dg_here%WAGP(I,chi)*dg_here%PHI_AREA(K,I,chi)
                     dg_here%BATH(I,J,chi) = dg_here%BATH(I,J,chi) + dg_here%HB(K,J,1)*dg_here%PHI_AREA(K,I,chi)
                     YELEM(chi) = YELEM(chi) + global_here%YDUB(K,J,chi)*dg_here%PHI_AREA(K,I,chi)
                     
                     IF (global_here%ICS.EQ.1) THEN
                        dg_here%SFAC_ELEM(I,J,chi)=1.0D0
                     ELSE
                        dg_here%SFAC_ELEM(I,J,chi)=COS(global_here%SFEA0)/COS(YELEM(chi)/R)
                     ENDIF

                     dg_here%DBATHDX(I,J,chi) = dg_here%DBATHDX(I,J,chi) + dg_here%HB(K,J,1)*dg_here%DPHIDX
                     dg_here%DBATHDY(I,J,chi) = dg_here%DBATHDY(I,J,chi) + dg_here%HB(K,J,1)*dg_here%DPHIDY
                     dg_here%DP_VOL(J,chi) = dg_here%DP_VOL(J,chi) + dg_here%WAGP(I,chi)*dg_here%HB(K,J,1)*dg_here%PHI_AREA(K,I,chi)

                  ENDDO

               enddo     

            endif

         enddo 

         do chi = 1,dg_here%ph
            dg_here%DP_VOL(J,chi) = 0.25D0*global_here%AREAS(J)*dg_here%DP_VOL(J,chi)
         enddo

         do chi = 1,dg_here%ph

            if (chi.ge.1) then

               DO L = 1,3
                  do I = 1,dg_here%NEGP(chi) ! global_here%Edge quadrature points

                     dg_here%BATHED(I,L,J,chi) = 0.D0
                     yed(chi) = 0.d0

                     DO K = 1,(chi+1)*(chi+2)/2
                        
                        dg_here%BATHED(I,L,J,chi) = dg_here%BATHED(I,L,J,chi)+dg_here%HB(K,J,1)*dg_here%PHI_EDGE(K,I,L,chi)

                        YED(chi) = YED(chi) + global_here%YDUB(K,J,chi)*dg_here%PHI_EDGE(K,I,L,chi)

                        IF (global_here%ICS.EQ.1) THEN
                           dg_here%SFACED(I,L,J,chi)=1.0d0
                        ELSE
                           dg_here%SFACED(I,L,J,chi)=COS(global_here%SFEA0)/COS(YED(chi)/R)
                        ENDIF

                        dg_here%EDGEQ(K,I,L,chi) = 2.0*dg_here%M_INV(K,chi)*dg_here%PHI_EDGE(K,I,L,chi)*dg_here%WEGP(I,chi)

                     ENDDO
                  enddo
               ENDDO

            else

               DO L = 1,3
                  do I = 1,dg_here%NEGP(chi) ! global_here%Edge quadrature points

                     dg_here%BATHED(I,L,J,chi) = 0.D0
                     yed(chi) = 0.d0

                     DO K = 1,(chi+1)*(chi+2)/2
                        
                        dg_here%BATHED(I,L,J,chi) = dg_here%BATHED(I,L,J,chi)+dg_here%HB(K,J,1)*dg_here%PHI_EDGE(K,I,L,chi)

                        YED(chi) = YED(chi) + global_here%YDUB(K,J,chi)*dg_here%PHI_EDGE(K,I,L,chi)

                        IF (global_here%ICS.EQ.1) THEN
                           dg_here%SFACED(I,L,J,chi)=1.0d0
                        ELSE
                           dg_here%SFACED(I,L,J,chi)=COS(global_here%SFEA0)/COS(YED(chi)/R)
                        ENDIF

                        dg_here%EDGEQ(K,I,L,chi) = 2.0*dg_here%M_INV(K,chi)*dg_here%PHI_EDGE(K,I,L,chi)*dg_here%WEGP(I,chi)

                     ENDDO
                  enddo
               ENDDO

            endif

         enddo

!........Store bathymetry at triangular vertices and global_here%edge gauss points for wet-dry
!........


         do chi = 1,dg_here%ph

            DO I = 1,3
               dg_here%DP_NODE(I,J,chi) = global_here%DP(global_here%NM(J,I))
            ENDDO


            IF (dg_here%NCHECK(chi).GT.3) THEN
               II = 4
               DO L = 1,3
                  DO I = 1,dg_here%NEGP(chi)
                     dg_here%DP_NODE(II,J,chi) = dg_here%BATHED(I,L,J,chi)
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
         dg_here%PSI_CHECK(1,I) = -1.D0/2.D0*(XI + YI)
         dg_here%PSI_CHECK(2,I) =  1.D0/2.D0*(XI + 1.D0)
         dg_here%PSI_CHECK(3,I) =  1.D0/2.D0*(YI + 1.D0)
         do chi = 1,dg_here%ph
            DO K = 1,(chi+1)*(chi+2)/2 
               dg_here%PHI_CHECK(K,I,chi) = dg_here%PHI_CORNER(K,I,chi)
            ENDDO
         enddo
      ENDDO

      do chi =1,dg_here%ph
         IF (dg_here%NCHECK(chi).GT.3) THEN
            II = 4
            DO L = 1,3
               DO I = 1,dg_here%NEGP(chi)
                  DO K = 1,(chi+1)*(chi+2)/2
                     
                     dg_here%PHI_CHECK(K,II,chi) = dg_here%PHI_EDGE(K,I,L,chi)
                     
                  ENDDO
                  II = II + 1
               ENDDO
            ENDDO
         ENDIF
      enddo

!.....Integrate the basis functions

      dg_here%PHI_INTEGRATED = 0.D0
      do chi = 1,dg_here%ph
         DO I = 1,dg_here%NAGP(chi)
            DO K = 1,(chi+1)*(chi+2)/2
               dg_here%PHI_INTEGRATED(K,chi) = dg_here%PHI_INTEGRATED(K,chi) + dg_here%WAGP(I,chi)*dg_here%PHI_AREA(K,I,chi)
            enddo
         ENDDO
      ENDDO

!.....Wetting and drying is not turned on

      IF(global_here%NOLIFA.EQ.0.OR.global_here%NOLIFA.EQ.1) THEN
         DO J = 1,S%MNE
            dg_here%WDFLG(J) = 1
                                !dg_here%DOFS(J) = 3
         ENDDO
         
!.....Wetting and drying is turned on but there are no dry nodes below
!.....geoid

      ELSEIF (global_here%NOLIFA.EQ.2.AND.global_here%NSTARTDRY.EQ.0) THEN

         DO J = 1,S%MNE
            
!.........Check to see if there are initially any dry nodes

                                !ZE1 = 0.D0
                                !ZE2 = 0.D0
                                !ZE3 = 0.D0 

            ZE1 = dg_here%ze(1,J,1)
            ZE2 = dg_here%ze(1,J,1)
            ZE3 = dg_here%ze(1,J,1)
            IF (global_here%DP(global_here%NM(J,1)).LT.global_here%H0) ZE1 = max(ze1,global_here%H0 - global_here%DP(global_here%NM(J,1)))
            IF (global_here%DP(global_here%NM(J,2)).LT.global_here%H0) ZE2 = max(ze2,global_here%H0 - global_here%DP(global_here%NM(J,2)))
            IF (global_here%DP(global_here%NM(J,3)).LT.global_here%H0) ZE3 = max(ze3,global_here%H0 - global_here%DP(global_here%NM(J,3)))
            
!.........If so set initial surface elevation values
            
            IF ((ZE1 + ZE2 + ZE3)/3.D0.ne.dg_here%ze(1,J,1)) THEN
               IF (p_0.EQ.0) THEN
                  DP_MIN = MIN(global_here%DP(global_here%NM(J,1)),global_here%DP(global_here%NM(J,2)),global_here%DP(global_here%NM(J,3)))
                  dg_here%ze(1,J,1) = max(dg_here%ze(1,j,1),global_here%H0 - DP_MIN)
               ELSE
                  IF (dg_here%ze(1,J,1).GT.(ZE1+ZE2+ZE3)/3.d0) THEN
                     IF (dof_0.ne.1) THEN
                        dg_here%ze(2,J,1)=0.d0
                        dg_here%ze(3,J,1)=0.d0
                        dg_here%ze(4:dg_here%dofh,J,1) = 0.D0 ! forced again for transparency
                     ENDIF
                  ELSE
                     dg_here%ze(1,J,1)=(ZE1+ZE2+ZE3)/3.D0
                     IF (DOF_0.ne.1) THEN
                        dg_here%ze(2,J,1) = -1.D0/6.D0*(ZE1 + ZE2) + 1.D0/3.D0*ZE3
                        dg_here%ze(3,J,1) = -0.5D0*ZE1 + 0.5D0*ZE2
                        dg_here%ze(4:dg_here%dofh,J,1) = 0.D0 ! forced again for transparency
                     ENDIF
                  ENDIF
               ENDIF

!............and set wet/dry flag (0 = dry, 1 = wet)

               dg_here%WDFLG(J) = 0
               dg_here%DOFS(J) = dof_0
            ELSE
               dg_here%WDFLG(J) = 1
               dg_here%DOFS(J) = dg_here%dofs(J)
            ENDIF

         ENDDO


!.....If there are dry nodes below geoid

      ELSEIF (global_here%NOLIFA.EQ.2.AND.global_here%NSTARTDRY.EQ.1) THEN
         
!.......Loop over elements
         
         DO J = 1,global_here%NE

!........Retrieve global node numbers for element

            global_here%N1 = global_here%NM(J,1)
            global_here%N2 = global_here%NM(J,2)
            global_here%N3 = global_here%NM(J,3)

!........Check to see if nodes are initially dry

            ZE1 = 0
            ZE2 = 0
            ZE3 = 0
            IF (STARTDRY(global_here%N1).EQ.1) ZE1 = global_here%H0 - global_here%DP(global_here%N1)
            IF (global_here%DP(global_here%N1).LT.global_here%H0) ZE1 = global_here%H0 - global_here%DP(global_here%N1)
            IF (STARTDRY(global_here%N2).EQ.1) ZE2 = global_here%H0 - global_here%DP(global_here%N2)
            IF (global_here%DP(global_here%N2).LT.global_here%H0) ZE2 = global_here%H0 - global_here%DP(global_here%N2)
            IF (STARTDRY(global_here%N3).EQ.1) ZE3 = global_here%H0 - global_here%DP(global_here%N3)
            IF (global_here%DP(global_here%N3).LT.global_here%H0) ZE3 = global_here%H0 - global_here%DP(global_here%N3)

            IF (dg_here%MODAL_IC.EQ.3) THEN
               IF (STARTDRY(global_here%N1).EQ.-88888) then
                  ZE1 = global_here%H0 - global_here%DP(global_here%N1)
               else
                  ZE1 = STARTDRY(global_here%N1)
               endif
               IF (STARTDRY(global_here%N2).EQ.-88888) then
                  ZE2 = global_here%H0 - global_here%DP(global_here%N2)
               else
                  ZE2 = STARTDRY(global_here%N2)
               endif
               IF (STARTDRY(global_here%N3).EQ.-88888) then
                  ZE3 = global_here%H0 - global_here%DP(global_here%N3)
               else
                  ZE3 = STARTDRY(global_here%N3)
               endif
            ENDIF

!.........If so set initial surface elevation values

            IF ((ZE1 + ZE2 + ZE3).ne.0) THEN
               IF (P_0.EQ.0) THEN
                  DP_MIN = MIN(global_here%DP(global_here%NM(J,1)),global_here%DP(global_here%NM(J,2)),global_here%DP(global_here%NM(J,3)))
                  dg_here%ze(1,J,1) = max(dg_here%ze(1,j,1),global_here%H0 - DP_MIN)
               ELSE
                  IF (dg_here%ze(1,J,1).GT.(ZE1+ZE2+ZE3)/3.d0) THEN
                     IF (DOF_0.ne.1) THEN
                        dg_here%ze(2,J,1)=0.d0
                        dg_here%ze(3,J,1)=0.d0
                        dg_here%ze(4:dg_here%dofh,J,1) = 0.D0 ! forced again for transparency
                     ENDIF
                  ELSE
                     dg_here%ze(1,J,1)=(ZE1+ZE2+ZE3)/3.D0
                     IF (DOF_0.ne.1) THEN
                        dg_here%ze(2,J,1) = -1.D0/6.D0*(ZE1 + ZE2) + 1.D0/3.D0*ZE3
                        dg_here%ze(3,J,1) = -0.5D0*ZE1 + 0.5D0*ZE2
                        dg_here%ze(4:dg_here%dofh,J,1) = 0.D0 ! forced again for transparency
                     ENDIF
                  ENDIF
               ENDIF

!...........Set wet/dry flag (0 = dry, 1 = wet)

               dg_here%WDFLG(J) = 0
               dg_here%DOFS(J) = DOF_0
            ELSE
               dg_here%WDFLG(J) = 1
               dg_here%DOFS(J) = dg_here%DOFS(J)
            ENDIF
         ENDDO
      ENDIF

      IF(MYPROC_HERE.EQ.0) THEN
         print *, 'DONE'
         print *, ''
      ENDIF

!.....Read in modal dg_here%dof for initial conditions
!Asserts error if you project onto lower order basis
!That is, do not expect convergence

      IF (dg_here%MODAL_IC.EQ.1) THEN
         OPEN(163,FILE=S%DIRNAME//'/'//'Initial_Conditions.163')
         OPEN(164,FILE=S%DIRNAME//'/'//'Initial_Conditions.164')
         OPEN(114,FILE=S%DIRNAME//'/'//'Initial_Bathymetry.114')
         READ(163,*) dg_here%P_READ
         READ(164,*) dg_here%P_READ,dg_here%P_READ
         READ(114,*) dg_here%P_READ
         IF (dg_here%P_READ.ne.dg_here%ph) THEN
            PRINT*,'INCONSISTENCY IN P -- CHECK INPUT FILES'
            STOP
         ENDIF
         DO J = 1,S%MNE
            DO K = 1,dg_here%DOFH
               READ(163,*) dg_here%ze(K,J,1)
               READ(164,*) dg_here%QX(K,J,1), dg_here%QY(K,J,1)
               READ(114,*) dg_here%HB(K,J,1)
            ENDDO
         ENDDO
         dg_here%H_TRI = SQRT((global_here%X(1)-global_here%X(2))**2 + (global_here%Y(1)-global_here%Y(2))**2)
         CLOSE(163)
         CLOSE(164)
         CLOSE(114)
      ENDIF

!.....Read in modal dg_here%dof for hot start conditions

      IF (dg_here%MODAL_IC.EQ.2) THEN
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
         READ(263,*) dg_here%P_READ
         READ(264,*) dg_here%P_READ,dg_here%P_READ
         READ(214,*) global_here%ITHS
         IF (dg_here%P_READ.ne.dg_here%PH) THEN
            PRINT*,'INCONSISTENCY IN P -- CHECK INPUT FILES'
            STOP
         ENDIF
         DO J = 1,S%MNE
            DO K = 1,dg_here%DOFH
               READ(263,*) dg_here%ze(K,J,1)
               READ(264,*) dg_here%QX(K,J,1), dg_here%QY(K,J,1)
               READ(214,*) dg_here%HB(K,J,1), dg_here%WDFLG(J)
#ifdef TRACE
               READ(288,*) dg_here%iota(K,J,1)
#endif
#ifdef CHEM
               READ(289,*) dg_here%iota(K,J,1),dg_here%iota2(K,J,1)
#endif
#ifdef DYNP
               READ(291,*) dg_here%dynP(K,J,1)
#endif
!$$$#ifdef SED_LAY
!$$$               do ll=1,layers
!$$$                  READ(290,*) dg_here%bed(K,J,1,ll)
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

      IF (ABS(global_here%NOUTGE).EQ.1) THEN
         OPEN(631,FILE=S%DIRNAME//'/'//'DG.63')
         WRITE(631,3220) global_here%RUNDES, global_here%RUNID, global_here%AGRID
         WRITE(631,3645) global_here%NDSETSE, dg_here%dofh, global_here%DTDP*global_here%NSPOOLGE, global_here%NSPOOLGE, 1
      ENDIF

!.....Initialize the DG.64 output file

      IF (ABS(global_here%NOUTGV).EQ.1) THEN
         OPEN(641,FILE=S%DIRNAME//'/'//'DG.64')
         WRITE(641,3220) global_here%RUNDES, global_here%RUNID, global_here%AGRID
         WRITE(641,3645) global_here%NDSETSV, dg_here%dofh, global_here%DTDP*global_here%NSPOOLGV, global_here%NSPOOLGV, 2
      ENDIF

!.....Initialize the DG.65 output file (contains elemental statuses such
!.....as the wet/dry status.

      IF ((ABS(global_here%NOUTGE).EQ.1).AND.(global_here%NOLIFA.GE.2)) THEN
         OPEN(651,FILE=S%DIRNAME//'/'//'DG.65')
         WRITE(651,3220) global_here%RUNDES, global_here%RUNID, global_here%AGRID
         WRITE(651,3645) global_here%NDSETSE, dg_here%dofh, global_here%DTDP*global_here%NSPOOLGE, global_here%NSPOOLGE, 1
      ENDIF
 3220 FORMAT(1X,A32,2X,A24,2X,A24)
 3645 FORMAT(1X,I10,1X,I10,1X,E15.7,1X,I5,1X,I5)

!.....Set p back to original value if p = 0

      IF (P_0.ne.dg_here%pl) THEN
         global_here%PDG_EL(:) = 0
         dg_here%PDG(:) = 0
         dg_here%DOF = 1
         dg_here%DOFL = 1
         dg_here%DOFS(:) = 1
         DOF_0 = 1
         dg_here%pl = 0
                                !p = 0
         dg_here%pa = 0
      ENDIF
      
!.....Compute basis functions at stations

      IF (global_here%NSTAE.GT.0) THEN      ! Elevation stations
         CALL ALLOC_STAE(dg_here, global_here%NSTAE )
         DO I = 1,global_here%NSTAE
            CALL STA_BASIS(dg_here, global_here, global_here%XEL(I), global_here%YEL(I),  global_here%NNE(I), dg_here%PHI_STAE(:,I) )
         ENDDO
      ENDIF
      
      IF (global_here%NSTAV.GT.0) THEN      ! Velocity Stations
         CALL ALLOC_STAV(dg_here, global_here%NSTAV )
         DO I = 1,global_here%NSTAV
            CALL STA_BASIS(dg_here, global_here, global_here%XEV(I), global_here%YEV(I),  global_here%NNV(I), dg_here%PHI_STAV(:,I) )
         ENDDO
      ENDIF

!.....Prep the slopelimiter

      IF (dg_here%SLOPEFLAG.ne.0) THEN
         IF(MYPROC_HERE.EQ.0)THEN
            print *, 'Slope limiting prep begins, "kshanti"'
         ENDIF
         CALL ALLOC_SLOPELIM(s,dg_here)
         CALL PREP_SLOPELIM(s,dg_here,global_here)
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

      SUBROUTINE RK_TIME(dg_here,global_here)

      USE GLOBAL
      USE DG

      IMPLICIT NONE

      type (dg_type) :: dg_here
      type (global_type) :: global_here
      
      INTEGER L,i,j,k,irk
      REAL(SZ) ARK, BRK, CASUM, MAX_BOA
      Real(SZ) eps_const,RKC_omega0,RKC_omega1

!.....Allocate the time stepping arrays

      dg_here%NRK = dg_here%RK_STAGE
      CALL ALLOC_RK(dg_here)

#ifdef RKSSP
!.....The forward Euler method

      IF ((dg_here%RK_STAGE.EQ.1).AND.(dg_here%RK_ORDER.EQ.1)) THEN

         dg_here%ATVD(:,:) = 0.D0
         dg_here%BTVD(:,:) = 0.D0
         dg_here%CTVD(:,:) = 0.D0
         dg_here%DTVD(:)   = 0.D0

         dg_here%ATVD(1,1) = 1.D0
         dg_here%BTVD(1,1) = 1.D0

!.....SSP(s,2) schemes

      ELSEIF (dg_here%RK_ORDER.EQ.2) THEN

         dg_here%ATVD(:,:) = 0.D0
         dg_here%BTVD(:,:) = 0.D0
         dg_here%CTVD(:,:) = 0.D0
         dg_here%DTVD(:)   = 0.D0

         DO I = 1,dg_here%NRK
            DO J = 0,dg_here%NRK-1

               IF ((J.EQ.(I-1)).AND.(I.LT.dg_here%NRK)) THEN
                  dg_here%ATVD(I,J+1) = 1.D0
                  dg_here%BTVD(I,J+1) = 1.D0/(dg_here%NRK-1)
               ELSEIF ((J.EQ.0).AND.(I.EQ.dg_here%NRK)) THEN
                  dg_here%ATVD(I,J+1) = 1.D0/dg_here%NRK
               ELSEIF ((J.EQ.(dg_here%NRK-1)).AND.(I.EQ.dg_here%NRK)) THEN
                  dg_here%ATVD(I,J+1) = (dg_here%NRK-1.D0)/dg_here%NRK
                  dg_here%BTVD(I,J+1) = 1.D0/dg_here%NRK
               ENDIF

            ENDDO
         ENDDO

!.....SSP(3,3) scheme

      ELSEIF ((dg_here%RK_STAGE.EQ.3).AND.(dg_here%RK_ORDER.EQ.3)) THEN

         dg_here%ATVD(:,:) = 0.D0
         dg_here%BTVD(:,:) = 0.D0
         dg_here%CTVD(:,:) = 0.D0
         dg_here%DTVD(:)   = 0.D0

         dg_here%ATVD(1,1) = 1.D0
         dg_here%ATVD(2,1) = 3.D0/4.D0
         dg_here%ATVD(2,2) = 1.D0/4.D0
         dg_here%ATVD(3,1) = 1.D0/3.D0
         dg_here%ATVD(3,3) = 2.D0/3.D0

         dg_here%BTVD(1,1) = 1.D0
         dg_here%BTVD(2,2) = 1.D0/4.D0
         dg_here%BTVD(3,3) = 2.D0/3.D0

!.....SSP(4,3) scheme

      ELSEIF ((dg_here%RK_STAGE.EQ.4).AND.(dg_here%RK_ORDER.EQ.3)) THEN

         dg_here%ATVD(:,:) = 0.D0
         dg_here%BTVD(:,:) = 0.D0
         dg_here%CTVD(:,:) = 0.D0
         dg_here%DTVD(:)   = 0.D0

         dg_here%ATVD(1,1) = 1.D0
         dg_here%ATVD(2,2) = 1.D0
         dg_here%ATVD(3,1) = 2.D0/3.D0
         dg_here%ATVD(3,3) = 1.D0/3.D0
         dg_here%ATVD(4,4) = 1.D0

         dg_here%BTVD(1,1) = 1.D0/2.D0
         dg_here%BTVD(2,2) = 1.D0/2.D0
         dg_here%BTVD(3,3) = 1.D0/6.D0
         dg_here%BTVD(4,4) = 1.D0/2.D0

!.....SSP(5,3) scheme

      ELSEIF ((dg_here%RK_STAGE.EQ.5).AND.(dg_here%RK_ORDER.EQ.3)) THEN

         dg_here%ATVD(:,:) = 0.D0
         dg_here%BTVD(:,:) = 0.D0
         dg_here%CTVD(:,:) = 0.D0
         dg_here%DTVD(:)   = 0.D0

         dg_here%ATVD(1,1) = 1.D0
         dg_here%ATVD(2,2) = 1.D0
         dg_here%ATVD(3,1) = 0.355909775063327D0
         dg_here%ATVD(3,3) = 0.644090224936674D0
         dg_here%ATVD(4,1) = 0.367933791638137D0
         dg_here%ATVD(4,4) = 0.632066208361863D0
         dg_here%ATVD(5,3) = 0.237593836598569D0
         dg_here%ATVD(5,5) = 0.762406163401431D0

         dg_here%BTVD(1,1) = 0.377268915331368D0
         dg_here%BTVD(2,2) = 0.377268915331368D0
         dg_here%BTVD(3,3) = 0.242995220537396D0
         dg_here%BTVD(4,4) = 0.238458932846290D0
         dg_here%BTVD(5,5) = 0.287632146308408D0

!.....SSP(6,3) scheme

      ELSEIF ((dg_here%RK_STAGE.EQ.6).AND.(dg_here%RK_ORDER.EQ.3)) THEN

         dg_here%ATVD(:,:) = 0.D0
         dg_here%BTVD(:,:) = 0.D0
         dg_here%CTVD(:,:) = 0.D0
         dg_here%DTVD(:)   = 0.D0

         dg_here%ATVD(1,1) = 1.D0
         dg_here%ATVD(2,2) = 1.D0
         dg_here%ATVD(3,3) = 1.D0
         dg_here%ATVD(4,1) = 0.476769811285196D0
         dg_here%ATVD(4,2) = 0.098511733286064D0
         dg_here%ATVD(4,4) = 0.424718455428740D0
         dg_here%ATVD(5,5) = 1.D0
         dg_here%ATVD(6,3) = 0.155221702560091D0
         dg_here%ATVD(6,6) = 0.844778297439909D0

         dg_here%BTVD(1,1) = 0.284220721334261D0
         dg_here%BTVD(2,2) = 0.284220721334261D0
         dg_here%BTVD(3,3) = 0.284220721334261D0
         dg_here%BTVD(4,4) = 0.120713785765930D0
         dg_here%BTVD(5,5) = 0.284220721334261D0
         dg_here%BTVD(6,6) = 0.240103497065900D0

!.....SSP(7,3) scheme

      ELSEIF ((dg_here%RK_STAGE.EQ.7).AND.(dg_here%RK_ORDER.EQ.3)) THEN

         dg_here%ATVD(:,:) = 0.D0
         dg_here%BTVD(:,:) = 0.D0
         dg_here%CTVD(:,:) = 0.D0
         dg_here%DTVD(:)   = 0.D0

         dg_here%ATVD(1,1) = 1.D0
         dg_here%ATVD(2,2) = 1.D0
         dg_here%ATVD(3,3) = 1.D0
         dg_here%ATVD(4,1) = 0.184962588071072D0
         dg_here%ATVD(4,4) = 0.815037411928928D0
         dg_here%ATVD(5,1) = 0.180718656570380D0
         dg_here%ATVD(5,2) = 0.314831034403793D0
         dg_here%ATVD(5,5) = 0.504450309025826D0
         dg_here%ATVD(6,6) = 1.D0
         dg_here%ATVD(7,4) = 0.120199000000000D0
         dg_here%ATVD(7,7) = 0.879801000000000D0

         dg_here%BTVD(1,1) = 0.233213863663009D0
         dg_here%BTVD(2,2) = 0.233213863663009D0
         dg_here%BTVD(3,3) = 0.233213863663009D0
         dg_here%BTVD(4,4) = 0.190078023865845D0
         dg_here%BTVD(5,5) = 0.117644805593912D0
         dg_here%BTVD(6,6) = 0.233213863663009D0
         dg_here%BTVD(7,7) = 0.205181790464579D0

!.....SSP(8,3) scheme

      ELSEIF ((dg_here%RK_STAGE.EQ.8).AND.(dg_here%RK_ORDER.EQ.3)) THEN

         dg_here%ATVD(:,:) = 0.D0
         dg_here%BTVD(:,:) = 0.D0
         dg_here%CTVD(:,:) = 0.D0
         dg_here%DTVD(:)   = 0.D0

         dg_here%ATVD(1,1) = 1.D0
         dg_here%ATVD(2,2) = 1.D0
         dg_here%ATVD(3,3) = 1.D0
         dg_here%ATVD(4,4) = 1.D0
         dg_here%ATVD(5,1) = 0.421366967085359D0
         dg_here%ATVD(5,2) = 0.005949401107575D0
         dg_here%ATVD(5,5) = 0.572683631807067D0
         dg_here%ATVD(6,2) = 0.004254010666365D0
         dg_here%ATVD(6,6) = 0.995745989333635D0
         dg_here%ATVD(7,3) = 0.104380143093325D0
         dg_here%ATVD(7,4) = 0.243265240906726D0
         dg_here%ATVD(7,7) = 0.652354615999950D0
         dg_here%ATVD(8,8) = 1.D0

         dg_here%BTVD(1,1) = 0.195804015330143D0
         dg_here%BTVD(2,2) = 0.195804015330143D0
         dg_here%BTVD(3,3) = 0.195804015330143D0
         dg_here%BTVD(4,4) = 0.195804015330143D0
         dg_here%BTVD(5,5) = 0.112133754621673D0
         dg_here%BTVD(6,6) = 0.194971062960412D0
         dg_here%BTVD(7,7) = 0.127733653231944D0
         dg_here%BTVD(8,8) = 0.195804015330143D0

!.....SSP(5,4) scheme

      ELSEIF ((dg_here%RK_STAGE.EQ.5).AND.(dg_here%RK_ORDER.EQ.4)) THEN

         dg_here%ATVD(:,:) = 0.D0
         dg_here%BTVD(:,:) = 0.D0
         dg_here%CTVD(:,:) = 0.D0
         dg_here%DTVD(:)   = 0.D0

         dg_here%ATVD(1,1) = 1.D0
         dg_here%ATVD(2,1) = 0.44437049406734D0
         dg_here%ATVD(2,2) = 0.55562950593266D0
         dg_here%ATVD(3,1) = 0.62010185138540D0
         dg_here%ATVD(3,3) = 0.37989814861460D0
         dg_here%ATVD(4,1) = 0.17807995410773D0
         dg_here%ATVD(4,4) = 0.82192004589227D0
         dg_here%ATVD(5,1) = 0.00683325884039D0
         dg_here%ATVD(5,3) = 0.51723167208978D0
         dg_here%ATVD(5,4) = 0.12759831133288D0
         dg_here%ATVD(5,5) = 0.34833675773694D0

         dg_here%BTVD(1,1) = 0.39175222700392D0
         dg_here%BTVD(2,2) = 0.36841059262959D0
         dg_here%BTVD(3,3) = 0.25189177424738D0
         dg_here%BTVD(4,4) = 0.54497475021237D0
         dg_here%BTVD(5,4) = 0.08460416338212D0
         dg_here%BTVD(5,5) = 0.22600748319395D0

!.....SSP(6,4) scheme

      ELSEIF ((dg_here%RK_STAGE.EQ.6).AND.(dg_here%RK_ORDER.EQ.4)) THEN

         dg_here%ATVD(:,:) = 0.D0
         dg_here%BTVD(:,:) = 0.D0
         dg_here%CTVD(:,:) = 0.D0
         dg_here%DTVD(:)   = 0.D0

         dg_here%ATVD(1,1) = 1.00000000000000D0
         dg_here%ATVD(2,1) = 0.30948026455053D0
         dg_here%ATVD(2,2) = 0.69051973544947D0
         dg_here%ATVD(3,1) = 0.54205244285557D0
         dg_here%ATVD(3,3) = 0.45794755714443D0
         dg_here%ATVD(4,1) = 0.35984960863377D0
         dg_here%ATVD(4,4) = 0.64015039136623D0
         dg_here%ATVD(5,5) = 1.00000000000000D0
         dg_here%ATVD(6,1) = 0.05776282890116D0
         dg_here%ATVD(6,3) = 0.44216432622405D0
         dg_here%ATVD(6,5) = 0.10115567086469D0
         dg_here%ATVD(6,6) = 0.39891717401009D0

         dg_here%BTVD(1,1) = 0.39270746575722D0
         dg_here%BTVD(2,2) = 0.30154043149172D0
         dg_here%BTVD(3,3) = 0.19997937335132D0
         dg_here%BTVD(4,4) = 0.27954483459696D0
         dg_here%BTVD(5,5) = 0.43668618869443D0
         dg_here%BTVD(6,3) = 0.09150931531680D0
         dg_here%BTVD(6,5) = 0.04417328437472D0
         dg_here%BTVD(6,6) = 0.14911300530736D0

!.....SSP(7,4) scheme

      ELSEIF ((dg_here%RK_STAGE.EQ.7).AND.(dg_here%RK_ORDER.EQ.4)) THEN

         dg_here%ATVD(:,:) = 0.D0
         dg_here%BTVD(:,:) = 0.D0
         dg_here%CTVD(:,:) = 0.D0
         dg_here%DTVD(:)   = 0.D0

         dg_here%ATVD(1,1) = 1.D0
         dg_here%ATVD(2,1) = 0.20161507213829D0
         dg_here%ATVD(2,2) = 0.79838492786171D0
         dg_here%ATVD(3,1) = 0.19469598207921D0
         dg_here%ATVD(3,3) = 0.80530401792079D0
         dg_here%ATVD(4,1) = 0.58143386885601D0
         dg_here%ATVD(4,4) = 0.41856613114399D0
         dg_here%ATVD(5,1) = 0.01934367892154D0
         dg_here%ATVD(5,5) = 0.98065632107846D0
         dg_here%ATVD(6,6) = 1.D0
         dg_here%ATVD(7,1) = 0.06006304558847D0
         dg_here%ATVD(7,3) = 0.30152730794242D0
         dg_here%ATVD(7,4) = 0.10518998496676D0
         dg_here%ATVD(7,5) = 0.01483791154585D0
         dg_here%ATVD(7,7) = 0.51838174995650D0

         dg_here%BTVD(1,1) = 0.30111872706068D0
         dg_here%BTVD(2,2) = 0.24040865318216D0
         dg_here%BTVD(3,3) = 0.24249212077315D0
         dg_here%BTVD(4,4) = 0.12603810060080D0
         dg_here%BTVD(5,5) = 0.29529398308716D0
         dg_here%BTVD(6,6) = 0.30111872706068D0
         dg_here%BTVD(7,3) = 0.09079551914158D0
         dg_here%BTVD(7,4) = 0.02888359354880D0
         dg_here%BTVD(7,7) = 0.15609445267839D0

!.....SSP(8,4) scheme

      ELSEIF ((dg_here%RK_STAGE.EQ.8).AND.(dg_here%RK_ORDER.EQ.4)) THEN

         dg_here%ATVD(:,:) = 0.D0
         dg_here%BTVD(:,:) = 0.D0
         dg_here%CTVD(:,:) = 0.D0
         dg_here%DTVD(:)   = 0.D0

         dg_here%ATVD(1,1) = 1.D0
         dg_here%ATVD(2,1) = 0.10645325745007D0
         dg_here%ATVD(2,2) = 0.89354674254993D0
         dg_here%ATVD(3,3) = 1.D0
         dg_here%ATVD(4,1) = 0.57175518477257D0
         dg_here%ATVD(4,4) = 0.42824481522743D0
         dg_here%ATVD(5,1) = 0.19161667219044D0
         dg_here%ATVD(5,5) = 0.80838332780956D0
         dg_here%ATVD(6,6) = 1.D0
         dg_here%ATVD(7,7) = 1.D0
         dg_here%ATVD(8,1) = 0.02580435327923D0
         dg_here%ATVD(8,3) = 0.03629901341774D0
         dg_here%ATVD(8,4) = 0.31859181340256D0
         dg_here%ATVD(8,5) = 0.05186768980103D0
         dg_here%ATVD(8,6) = 0.03944076217320D0
         dg_here%ATVD(8,7) = 0.00511633747411D0
         dg_here%ATVD(8,8) = 0.52288003045213D0

         dg_here%BTVD(1,1) = 0.24120020561311D0
         dg_here%BTVD(2,2) = 0.21552365802797D0
         dg_here%BTVD(3,3) = 0.24120020561311D0
         dg_here%BTVD(4,4) = 0.10329273748560D0
         dg_here%BTVD(5,5) = 0.19498222488188D0
         dg_here%BTVD(6,6) = 0.24120020561311D0
         dg_here%BTVD(7,7) = 0.24120020561311D0
         dg_here%BTVD(8,3) = 0.00875532949991D0
         dg_here%BTVD(8,4) = 0.06195575835101D0
         dg_here%BTVD(8,6) = 0.00951311994571D0
         dg_here%BTVD(8,8) = 0.12611877085604D0

      ENDIF

!.....Compute the time dependent parameters

      DO K = 0,dg_here%RK_STAGE-1
         DO I = 1,dg_here%RK_STAGE
            CASUM = 0.D0
            DO L = K+1,I-1
               CASUM = CASUM + dg_here%CTVD(L,K+1)*dg_here%ATVD(I,L+1)
            ENDDO
            dg_here%CTVD(I,K+1) = dg_here%BTVD(I,K+1) + CASUM
         ENDDO
      ENDDO

      DO K = 1,dg_here%RK_STAGE-1
         dg_here%DTVD(K+1) = 0.D0
         DO L = 0,K-1
            dg_here%DTVD(K+1) = dg_here%DTVD(K+1) + dg_here%CTVD(K,L+1)
         ENDDO
      ENDDO

!.....Compute the maximum beta over global_here%alpha ratio at each stage

      DO irk = 1,dg_here%RK_STAGE
         MAX_BOA = 0.D0
         DO I = 1,irk
            ARK = dg_here%ATVD(irk,I)
            BRK = dg_here%BTVD(irk,I)
            IF (ARK.ne.0.D0) THEN
               IF (MAX_BOA.LT.BRK/ARK) MAX_BOA = BRK/ARK
            ENDIF
         ENDDO
         dg_here%MAX_BOA_DT(irk) = MAX_BOA*global_here%DT
      ENDDO

#endif

!-----------------------------------------------------------
!.... Compute the Runge-Kutta Chebyshev (RKC) version

#ifdef RKC
         
         RKC_omega0 = 0.D0
         RKC_omega1 = 0.D0

         eps_const = 0.153846154D0 

                                !Compute the RKC version

         dg_here%RKC_T(0) = 1.D0
         dg_here%RKC_U(0) = 1.D0
         
         RKC_omega0 = 1.D0 + eps_const*dg_here%nrk**(-2.D0)

         dg_here%RKC_c(0) = 0.D0
         dg_here%RKC_b(0) = 0.D0
         dg_here%RKC_b(1) = 1.D0/RKC_omega0

         dg_here%RKC_T(1) = RKC_omega0
         dg_here%RKC_U(1) = 2.D0*RKC_omega0

         dg_here%RKC_Tprime(1) = dg_here%nrk**(2.D0)
         dg_here%RKC_Tdprime(1) = 1.D0/3.D0 * (dg_here%nrk**(2.D0)*(dg_here%nrk**(2.D0) - 1.D0))


         do j = 0,dg_here%nrk

                                !Chebyshev polynomial stuff

            if(j.gt.1) then

               dg_here%RKC_T(j) = 2.0 * RKC_omega0 * dg_here%RKC_T(j-1) - dg_here%RKC_T(j-2)
               dg_here%RKC_U(j) = 2.0 * RKC_omega0 * dg_here%RKC_U(j-1) - dg_here%RKC_U(j-2)
               
            endif

            if(j.gt.1) then

               dg_here%RKC_Tprime(j) = j*dg_here%RKC_U(j-1)
               dg_here%RKC_Tdprime(j) = (j*(j*dg_here%RKC_T(j) - RKC_omega0*dg_here%RKC_U(j-1))/( RKC_omega0**(2.D0) - 1.D0 ) )
               dg_here%RKC_b(j) = dg_here%RKC_Tdprime(j)*dg_here%RKC_Tprime(j)**(-2.D0)

            endif
            
            dg_here%RKC_a(j)  = 1.D0 - dg_here%RKC_b(j)*dg_here%RKC_T(j)
            
         enddo

         if(dg_here%RK_stage.gt.1) then

            dg_here%RKC_b(0) = dg_here%RKC_b(2)

         endif



         RKC_omega1 =  dg_here%RKC_Tprime(dg_here%nrk)/dg_here%RKC_Tdprime(dg_here%nrk)

         dg_here%RKC_tildemu(1)= dg_here%nrk**(-1.D0) * RKC_omega1
         
         do j=2,dg_here%nrk

            dg_here%RKC_mu(j) = 2.D0*dg_here%RKC_b(j)*RKC_omega0/dg_here%RKC_b(j-1)
            
            dg_here%RKC_tildemu(j) = 2.D0*dg_here%RKC_b(j)*RKC_omega1/dg_here%RKC_b(j-1)
            
            dg_here%RKC_nu(j) = - dg_here%RKC_b(j)/dg_here%RKC_b(j-2)
            
            dg_here%RKC_gamma(j) = -dg_here%RKC_a(j-1)*dg_here%RKC_tildemu(j)
            
            dg_here%RKC_c(j) = (j**(2.D0)-1.0) / ( dg_here%nrk**(2.D0) - 1.D0 )

         enddo

         dg_here%RKC_c(1) = (dg_here%RKC_c(2)) / (4.D0*RKC_omega0)

         dg_here%RKC_c(dg_here%RK_stage) = 0.D0
         dg_here%RKC_c(dg_here%RK_stage) = 1.D0

#endif


      DEALLOCATE(dg_here%CTVD)
      RETURN
      END SUBROUTINE RK_TIME
