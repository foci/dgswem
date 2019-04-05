C***********************************************************************
C     
C     SUBROUTINE DG_PREP()
C     
C     This subroutine does preparatory stuff for DG
C     
C     Written by Ethan Kubatko (03-07-2005)
C     
C-----------------------------------------------------------------------
C     
C     mod history of hp_ADCIRC since v9
C     
C     v9_sb1       - Aug 05 - sb - parallelized
C     v9_sb2       - Aug 05 - sb - wetting/drying
C     v10_sb1      - Aug 05 - sb - reflect fort.12as an initial surface 
C                                  elevation
C     v10_sb5      - Oct 03 - sb - consolidate Ethan's slope limiter
C     v10_sb5      - Oct 24 - sb - DG.65 output is added
C     - Aug 29, 2007 - Fixed bug for startdry = 1
C     01-10-2011 - cem - adapted for p_enrichment and multicomponent
C     06-01-2012 - cem -sediment added
C     
C***********************************************************************

      SUBROUTINE PREP_DG()

C.....Use appropriate modules

      USE GLOBAL
      USE DG
      USE NodalAttributes, ONLY : STARTDRY, FRIC, GeoidOffset,
     $     LoadGeoidOffset,LoadManningsN,ManningsN
#ifdef CMPI
      USE MESSENGER_ELEM
      USE MESSENGER
#endif

      IMPLICIT NONE
      
C.....Declare local variables

      INTEGER II,l,P_0,DOF_0,j,k,kk,jj,jjj,i,chi,ll,mm,Q,M,P,SZ2,w,III
      CHARACTER(LEN=8) REGION
      REAL(SZ) AREA, ANGLE_SUM, HBB(MXEDGS),CASUM, DP_MIN,temp_lay,XP,YP
      REAL(SZ) XI,YI,ZE1,ZE2,ZE3,ZE4,l2er,l2erh2,xcen,ycen,epsl,pi_n
      REAL(SZ) ZP(MXEDGS), DHBX, ell_1,ell_2,ell_3,int_hb,int_ze,int_yd
      REAL(SZ), Allocatable :: PTS(:,:), VERT(:,:,:), DBASIS(:,:)
      REAL(SZ), Allocatable :: BARY(:,:), WTS(:), BASIS(:), PT(:)
      REAL(SZ), Allocatable :: PHI_AREAM(:,:), PHI_AREAMT(:,:), WM(:,:)
      REAL(SZ), Allocatable :: MTEMP1(:,:),MTEMP2(:,:),MTEMPI(:,:)
      REAL(SZ), Allocatable :: vtemp1(:),vtemp2(:),PHI_AREAV(:)

      real(sz) checkarea,arint(2,2),rhsint(2),edgeint,dpsdx,psimid,
     $     determ,sfacdub2max,sfacdub3max,R,Jac,XX,YY
      integer i1,i2,sfac_flag,led,ADDGP(2),NQEDS
      integer ifac2max,ifac3max,phh,ELT,pp,NN,KI(MNP)
      logical Logi

      real(sz) xmid,ymid,Ox,Oy,rev,C_0,sig,C_1,Aq,Atr1,Atr2

      Real(SZ),allocatable :: XBCbt(:),YBCbt(:),radial(:),XB(:),YB(:),
     & l2e(:)
      Real(SZ),allocatable :: iota_check(:),iota_check2(:),hbo(:,:,:), 
     & ydubo(:,:)
      Real(SZ),allocatable :: YELEM(:),YED(:),HB1(:,:,:,:), zeo(:,:,:)

      Allocate ( XBCbt(MNE),YBCbt(MNE),radial(MNE),XB(MNE),YB(MNE),
     & l2e(MNE) )
      Allocate ( iota_check(MNE),iota_check2(MNE),hbo(36,MNE,1),
     & ydubo(36,mne) ) 
      Allocate ( YELEM(ph),YED(ph),hb1(36,mne,1,ph), zeo(36,mne,1), 
     & BARY(DIM,2), VERT(4,DIM,2) )
      Allocate ( PT(DIM) )
      
      C13 = 1.D0/3.D0
      C16 = 1.D0/6.D0
      R = 6378206.4d0

C     sb-PDG1 moved from other places

C.....Obtain RK time scheme parameters

      CALL RK_TIME()
      
C.....Compute the degrees of freedom per element

      DOF = (pl+1)*(pl+2)/2
      dofx = (px+1)*(px+2)/2    ! dofx for variable functions f=f(x) 
      P_0 = pl
      DOF_0 = (pl+1)*(pl+2)/2   ! dof at lowest order when p!=0
      dofh = (ph + 1)*(ph + 2)/2


C.....Allocate some DG stuff

      IF (PADAPT.EQ.1) THEN

         dofh = (ph + 1)*(ph + 2)/2
         dofl = (pl + 1)*(pl + 2)/2
         pa = pl
           
      elseif (padapt.eq.0) then 

         dofh = dofh
         dofl = DOF_0
         pa = pl

      endif

#ifdef SED_LAY
C.....Initialize funtion parser for sediment types
      init_parser = .false.
#endif

#ifdef SED_LAY
C.....Initialize stabilizer sweep for A.D.
      stblzr = .false.
#endif


C.....Compute the number of gauss points needed for the edge integrals

      if (MXEDGS .eq. 4) THEN
         Logi = .true.
      elseif (MXEDGS .eq. 3) THEN
         Logi = .false.
      endif
      CALL ALLOC_DG4(Logi)          !moved here 6.28.10, for p_adapt because of messenger_elem      

         dofs(:) = dofl
         PDG_EL(:) = pl
         PDG(:) = pl  
         PCOUNT(:) = 0
         pa = pl

      do chi=pl,ph
         NEGP(chi) = chi + 1
      enddo

      IF (pl.eq.0) THEN

         PDG_EL(:) = 1
         PDG(:) = 1     
         DOF    = 3
         pl     = 1
         dofl   = 3
         P_0    = 0
         DOF_0  = 1 
         NEGP(pl) = 2

      ENDIF

!     cnd
!     iwrite=0
      
C.....Initilization for parallel DG run

#ifdef CMPI

      CALL MSG_TYPES_ELEM()     ! Determine Word Sizes for Message-Passing
      CALL MSG_TABLE_ELEM()     ! Read Message-Passing Tables

      IF (SLOPEFLAG.ge.4) THEN
         CALL MSG_TYPES()
         CALL MSG_TABLE()
      ENDIF

#endif

C.....Create the edge based data

      IF(MYPROC.EQ.0) THEN
         PRINT*, 'CREATING EDGE DATA...'
         PRINT*, ''
      ENDIF
      CALL CREATE_EDGE_DATA()
      IF(MYPROC.EQ.0) THEN
         print *, 'CREATING EDGE DATA DONE'
         print *, ''
      ENDIF

#ifdef CMPI
      CALL MESSAGE_START_ELEM() ! Startup persistent message passing
      IF (SLOPEFLAG.ge.4) CALL MESSAGE_START()
#endif

C.....Re-arrange elevation specified boundary segment data for DG

      IF (NEEDS.GT.0) THEN
         CALL ALLOC_DG1(MNBFR)
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

C.....Re-arrange non-zero flow specified boundary segment data for DG

      IF (NFEDS.GT.0) THEN
         CALL ALLOC_DG2(MNFFR)
         II = 1
         JJ = 1
         DO I = 1,MNFFR
            IF (NFFR .EQ. 0) THEN
               QTRATIO = (TIMEDG - QTIME1)/FTIMINC
               NQEDS = 0
               DO J = 1,NVEL
                  IF((LBCODEI(J).EQ.2).OR.(LBCODEI(J).EQ.12)
     &                                .OR.(LBCODEI(J).EQ.22)) THEN 
                     NQEDS = NQEDS + 1
                     IF (NQEDS .LE. NFEDS) THEN
                        QNAM_DG(1,NQEDS,1) = RAMPDG*( QNIN1(J) + 
     &                                 QTRATIO*( QNIN2(J) - QNIN1(J) ) )
                        QNPH_DG(1,NQEDS,1) = 0.D0
                        QNAM_DG(1,NQEDS,2) = RAMPDG*( QNIN1(J+1) + 
     &                             QTRATIO*( QNIN2(J+1) - QNIN1(J+1) ) )
                        QNPH_DG(1,NQEDS,2) = 0.D0
                     ENDIF
                  ENDIF
               ENDDO
            ELSE
               DO J = 1,NBOU
                  IF ( (SEGTYPE(J).EQ.2 ).OR.(SEGTYPE(J).EQ.12)
     &                                  .OR.(SEGTYPE(J).EQ.22) ) THEN
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
            ENDIF
         ENDDO
      ENDIF
      
C.....If there are internal barriers allocate some stuff

      IF (NIBEDS.NE.0) CALL ALLOC_DG3(MNP)

C.....Allocate the array for node to element table

      CALL ALLOC_NNOEL1(MNP)

C.....Determine the number of elements connected at each node

      EL_COUNT = 0
      MAXEL = 1

      DO J = 1,MNE
         DO K = 1,EDN(J)
            N1 = NM(J,K)
            EL_COUNT(N1) = EL_COUNT(N1) + 1
         ENDDO
      ENDDO
      MAXEL = MAXVAL(EL_COUNT)

C.....Allocate the array for the node to element table

      CALL ALLOC_NNOEL2(MNP,MAXEL)
      
C.....Construct node to element table

      EL_COUNT = 0
      DO J = 1,MNE
         DO K = 1,EDN(J)
            N1 = NM(J,K)
            NNOEL(N1,1+EL_COUNT(N1)) = J
            EL_COUNT(N1) = EL_COUNT(N1) + 1
         ENDDO
      ENDDO

C.....Construct node to element angle table
      DO I = 1,MNP
         ETAMAX(I) = -99999
         KK = 1
         ELETAB(I,1) = I
         S1 = SFAC(I)
         J1 = NEITAB(I,1)

         DO 111 K = 1,NNEIGH(I)-1
            J2 = NEITAB(I,1+K)
            DO J = 1,EL_COUNT(I)
               EL = NNOEL(I,J)
               JJ = EDN(EL) - 2
               J2 = NEITAB(I,1+K)
               IF (K.LT.(NNEIGH(I)-JJ)) THEN
                  J3 = NEITAB(I,JJ+1+K)
               ELSE
                  J3 = NEITAB(I,2)
               ENDIF
               IF (EDN(EL) .EQ. 3) THEN
                  N1 = NM(EL,1)
                  N2 = NM(EL,2)
                  N3 = NM(EL,3)
                  IF ((J1.EQ.N1).OR.(J1.EQ.N2).OR.(J1.EQ.N3)) THEN
                     IF ((J2.EQ.N1).OR.(J2.EQ.N2).OR.(J2.EQ.N3)) THEN
                        IF ((J3.EQ.N1).OR.(J3.EQ.N2).OR.(J3.EQ.N3)) THEN
                           ELETAB(I,1+KK) = EL
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
               ELSEIF (EDN(EL) .EQ. 4) THEN
                  N1 = NM(EL,1)
                  N2 = NM(EL,2)
                  N3 = NM(EL,3)
                  N4 = NM(EL,4)
                  IF (J1.EQ.N1) THEN
                     IF ((J2.EQ.N1).OR.(J2.EQ.N2).OR.(J2.EQ.N4)) THEN
                        IF((J3.EQ.N1).OR.(J3.EQ.N2).OR.(J3.EQ.N4)) THEN
                           ELETAB(I,1+KK) = EL
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
                  ELSEIF (J1.EQ.N2) THEN
                     IF((J2.EQ.N1).OR.(J2.EQ.N2).OR.(J2.EQ.N3)) THEN
                        IF((J3.EQ.N1).OR.(J3.EQ.N2).OR.(J3.EQ.N3)) THEN
                           ELETAB(I,1+KK) = EL
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
                  ELSEIF (J1.EQ.N3) THEN
                     IF((J2.EQ.N2).OR.(J2.EQ.N3).OR.(J2.EQ.N4)) THEN
                        IF((J3.EQ.N2).OR.(J3.EQ.N3).OR.(J3.EQ.N4)) THEN
                           ELETAB(I,1+KK) = EL
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
                  ELSEIF (J1.EQ.N4) THEN
                     IF((J2.EQ.N1).OR.(J2.EQ.N3).OR.(J2.EQ.N4)) THEN
                        IF((J3.EQ.N1).OR.(J3.EQ.N3).OR.(J3.EQ.N4)) THEN
                           ELETAB(I,1+KK) = EL
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
                  !IF ((J1.EQ.N1).OR.(J1.EQ.N2).OR.(J1.EQ.N3).OR.(J1.EQ.N4)) THEN
                  !   IF((J2.EQ.N1).OR.(J2.EQ.N2).OR.(J2.EQ.N3).OR.(J2.EQ.N4)) THEN
                  !      IF((J3.EQ.N1).OR.(J3.EQ.N2).OR.(J3.EQ.N3).OR.(J3.EQ.N4)) THEN
                  !         ELETAB(I,1+KK) = EL
                  !         S2  = SFAC(J2)
                  !         SAV = (S1 + S2)/2.D0
                  !         VEC1(1) =      X(J1) - X(J2)
                  !         VEC1(2) = SAV*(Y(J1) - Y(J2))
                  !         VEC2(1) =      X(J1) - X(J3)
                  !         VEC2(2) = SAV*(Y(J1) - Y(J3))
                  !         MAG1 = SQRT(VEC1(1)**2 + VEC1(2)**2)
                  !         MAG2 = SQRT(VEC2(1)**2 + VEC2(2)**2)
                  !         DOT = DOT_PRODUCT(VEC1,VEC2)
                  !         EL_ANG  = ACOS(DOT/(MAG1*MAG2))
                  !         ANGTAB(I,KK+1) = RAD2DEG*EL_ANG
                  !         KK = KK + 1
                  !         GOTO 111
                  !      ENDIF
                  !   ENDIF
                  !ENDIF
               ENDIF
            ENDDO
 111     CONTINUE
         KI(I) = KK
      ENDDO

      
C.....Delete any duplicate entries in ELETAB and ANGTAB
C.....Typically only occurs for nodes with just 1 element connection.
      DO I = 1,MNP
         DO j = 2,KI(I)-1
            do while (ELETAB(i,j) .eq. ELETAB(i,j+1) .and. ELETAB(i,j) .ne. 0)
               do jj = j+1,KI(I)
                  if (jj .ne. KI(I)) then
                      ELETAB(i,jj) = ELETAB(i,jj+1)
                      ANGTAB(i,jj) = ANGTAB(i,jj+1)
                  else
                      ELETAB(i,jj) = 0
                      ANGTAB(i,jj) = 0
                  endif
               enddo
            enddo
         ENDDO
      ENDDO
 
!     C.....Allocate some DG stuff

!     CALL ALLOC_DG4()
      
C.....Initialize the DG arrays

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

C.....If using modal initial conditions transform the bathymetry from
C.....nodal coordinates to modal dof
      ! COME BACK TO THIS!
      DO J = 1,MNE
         N1 = NM(J,1)
         N2 = NM(J,2)
         N3 = NM(J,3)
         hbo(1,J,1) =  1.D0/3.D0 * (DP(N1) + DP(N2) + DP(N3))
         hbo(2,J,1) = -1.D0/6.D0 * (DP(N1) + DP(N2)) + 1.D0/3.D0*DP(N3)
         hbo(3,J,1) = -0.5D0*DP(N1) + 0.5D0*DP(N2)

         ydubo(1,J)= 1.D0/3.D0*(Y(N1) + Y(N2) + Y(N3))
         ydubo(2,J) = -1.D0/6.D0*(Y(N1) + Y(N2))
     $        + 1.D0/3.D0*Y(N3)
         ydubo(3,J) = -0.5D0*Y(N1) + 0.5D0*Y(N2)

      ENDDO

      ! COME BACK TO THIS!
      if (LoadManningsN) then
         DO J = 1,NE
            N1 = NM(J,1)
            N2 = NM(J,2)
            N3 = NM(J,3)
            MANN(1,J) =  1.D0/3.D0*(ManningsN(N1)
     $           + ManningsN(N2) + ManningsN(N3))
            MANN(2,J) = -1.D0/6.D0*(ManningsN(N1) 
     $           + ManningsN(N2)) + 1.D0/3.D0*ManningsN(N3)
            MANN(3,J) = -0.5D0*ManningsN(N1) + 0.5D0*ManningsN(N2)
         ENDDO
      endif
      ! COME BACK TO THIS!
      IF (MODAL_IC.EQ.0) THEN
C     this assumes a cold start
         if (LoadGeoidOffset) then
            DO J = 1,NE
               N1 = NM(J,1)
               N2 = NM(J,2)
               N3 = NM(J,3)
               zeo(1,J,1)=1.d0/3.d0*(GeoidOffset(N1)+GeoidOffset(N2)+
     $              GeoidOffset(N3))
               IF (dof_0.NE.1) THEN
                 zeo(2,J,1)=-1.d0/6.d0*(GeoidOffset(N1)+GeoidOffset(N2))
     $                 +1.d0/3.d0*GeoidOffset(N3)
                  zeo(3,J,1)=-.5d0*GeoidOffset(N1)+.5d0*GeoidOffset(N2)
               ENDIF
            ENDDO
         endif
      ENDIF

C--   

!As part of initializing the system, let us determine the partials
!of the sediment discharge equation, fed in by fort.dg

#ifdef SED_LAY
      IF(MYPROC.EQ.0)THEN
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

      IF(MYPROC.EQ.0) THEN
         print *, 'PREP FOR WET/DRY BEGINS...'
      ENDIF

C.....1. Set initial surface elevation above the bed elevation
C.....if wetting-and-drying is enabled and the initial water depth is
C.....not specified by fort.12
C.....2. Set wet-and-dry elemental flags
C.....3. Set the DOF at dry elements = 1
      do j=1,ph
          jj = 2*j
          NEGP(j)  = CEILING((jj+3)/2.0d0)
      enddo

      NCHECK(1,1) = MNEDGS
      NCHECK(1,2) = MXEDGS
      if (ph.gt.1) then
         do chi = 2,ph
            NCHECK(chi,1) = NCHECK(1,1) + MNEDGS*negp(chi)
            NCHECK(chi,2) = NCHECK(1,2) + MXEDGS*negp(chi)
         enddo
      endif
      
      CALL ALLOC_DG_WETDRY()
      PHI_CHECK = 0.D0
      PSI_CHECK = 0.D0
      H0L = H0
      H0H = H0 * 1.0
      HABSMIN = H0 * 1.0
            
C.....Retrieve the normals to the edges

      CALL CALC_NORMAL()

C.....Retrieve the area integral gauss quadrature points
      
      do j=1,ph
         if (j.eq.1) then
            phh = 2*ph
            do ELT = MNEDGS,MXEDGS
               if (ELT.eq.3 .and. DIM.eq.2) then
                  if (phh.le.0) then
                     NAGP(ph,1) = 1
                  elseif (phh.le.2) then
                     NAGP(ph,1) = 3
                  elseif (phh.le.4) then
                     NAGP(ph,1) = 6
                  elseif (phh.le.6) then
                     NAGP(ph,1) = 12
                  elseif (phh.le.8) then
                     NAGP(ph,1) = 16
                  elseif (phh.le.10) then
                     NAGP(ph,1) = 25
                  elseif (phh.le.12) then
                     NAGP(ph,1) = 33
                  elseif (phh.le.14) then
                     NAGP(ph,1) = 42
                  elseif (phh.le.16) then
                     NAGP(ph,1) = 55
                  elseif (phh.le.18) then
                     NAGP(ph,1) = 72
                  elseif (phh.le.20) then
                     NAGP(ph,1) = 88
                  endif
               elseif (ELT.eq.4 .and. DIM.eq.2) then
                  if (phh.le.0) then
                     NAGP(ph,2) = 1
                  elseif (phh.le.2) then
                     NAGP(ph,2) = 3
                  elseif (phh.le.4) then
                     NAGP(ph,2) = 6
                  elseif (phh.le.6) then
                     NAGP(ph,2) = 10
                  elseif (phh.le.8) then
                     NAGP(ph,2) = 16
                  elseif (phh.le.10) then
                     NAGP(ph,2) = 22
                  elseif (phh.le.12) then
                     NAGP(ph,2) = 31
                  elseif (phh.le.14) then
                     NAGP(ph,2) = 44
                  elseif (phh.le.16) then
                     NAGP(ph,2) = 56
                  elseif (phh.le.18) then
                     NAGP(ph,2) = 68
                  elseif (phh.le.20) then
                     NAGP(ph,2) = 81
                  elseif (phh.le.22) then
                     NAGP(ph,2) = 100
                  endif
               endif
            enddo
C...........Allocate XaGP, YaGP, and WaGP
            MXNAGP = max(NAGP(ph,1),NAGP(ph,2))
            call ALLOC_AREA_GAUSS(Logi)
         endif
         do ELT = MNEDGS,MXEDGS
			   jj = 2*j
            if (ELT.eq.3 .and. DIM.eq.2) then
               REGION = 'TRIANGLE'
               if (jj.le.0) then
                  SZ2 = 1
               elseif (jj.le.2) then
                  SZ2 = 3
               elseif (jj.le.4) then
                  SZ2 = 6
               elseif (jj.le.6) then
                  SZ2 = 12
               elseif (jj.le.8) then
                  SZ2 = 16
               elseif (jj.le.10) then
                  SZ2 = 25
               elseif (jj.le.12) then
                  SZ2 = 33
               elseif (jj.le.14) then
                  SZ2 = 42
               elseif (jj.le.16) then
                  SZ2 = 55
               elseif (jj.le.18) then
                  SZ2 = 72
               elseif (jj.le.20) then
                  SZ2 = 88
               endif 
            elseif (ELT.eq.4 .and. DIM.eq.2) then
               REGION = 'SQUARE  '
               if (jj.le.0) then
                  SZ2 = 1
               elseif (jj.le.2) then
                  SZ2 = 3
               elseif (jj.le.4) then
                  SZ2 = 6
               elseif (jj.le.6) then
                  SZ2 = 10
               elseif (jj.le.8) then
                  SZ2 = 16
               elseif (jj.le.10) then
                  SZ2 = 22
               elseif (jj.le.12) then
                  SZ2 = 31
               elseif (jj.le.14) then
                  SZ2 = 44
               elseif (jj.le.16) then
                  SZ2 = 56
               elseif (jj.le.18) then
                  SZ2 = 68
               elseif (jj.le.20) then
                  SZ2 = 81
               elseif (jj.le.22) then
                  SZ2 = 100
               endif
            endif

            Allocate( PTS(SZ2,DIM), WTS(SZ2) )
            NN = 0

            call QUADRATURE( jj,REGION,PTS,WTS,NN,SZ2,DIM )
            
            NAGP(j,ELT-2) = NN
            do jj = 1,SZ2 
               XAGP(jj,j,ELT-2) = PTS(jj,1)
               YAGP(jj,j,ELT-2) = PTS(jj,2)
               WAGP(jj,j,ELT-2) = WTS(jj)
            enddo
            Deallocate( PTS, WTS )
         enddo
      enddo

C.....Retrieve the edge integral gauss quadrature points
      
      do j=1,ph

          if (j.eq.1) then

             NEGP(ph) = CEILING((phh+3)/2.0d0)
             CALL ALLOC_EDGE_GAUSS(Logi)

          endif
          jj = 2*j
          SZ2 = CEILING((jj+3)/2.0d0)
          NEGP(j) = SZ2

          Allocate( PTS(SZ2,1), WTS(SZ2) )
          CALL quad_rules_general(NEGP(j),.true.,0.0d0,0.0d0,PTS,WTS)
     
          XEGP(1:NEGP(j),j) = PTS(1:NEGP(j),1)
          WEGP(1:NEGP(j),j) = WTS(1:NEGP(j))
          Deallocate( PTS, WTS )
      enddo

C.....Evaluate the orthogonal basis and its derivatives at the area
C.....gauss quadrature points

C.....Determine vertices and element barycenter by elements type

      BARY(:,1) = (/ -1.D0/3.D0, -1.D0/3.D0/)
      VERT(1,:,1) = (/ -1.D0, -1.D0 /)
      VERT(2,:,1) = (/ 1.D0, -1.D0 /)
      VERT(3,:,1) = (/ -1.D0, 1.D0 /)
      ADDGP(1) = 3
      if (MXEDGS .eq. 4) then
         BARY(:,2) = (/ 0.D0, 0.D0 /)
         VERT(1,:,2) = (/ -1.D0, -1.D0 /)
         VERT(2,:,2) = (/ 1.D0, -1.D0 /)
         VERT(3,:,2) = (/ 1.D0, 1.D0 /)
         VERT(4,:,2) = (/ -1.D0, 1.D0 /)
         ADDGP(2) = 4
      endif

C.....Loop over orders of polynomials required
      do w=1,ph
         P = w
         L = w
         SZ2 = (P+2)*(P+1)**(DIM-1)/2
         Allocate ( BASIS(SZ2),DBASIS(SZ2,DIM) )

C........Allocate JACOBI (orthobasis routine is no longer used so
C........ALLOC_JACOBI is no longer needed)

C$$$        if (L.eq.1 .and. p.eq.1) then
C$$$ 
C$$$           call ALLOC_JACOBI()
C$$$
C$$$        endif

C.......Calculate orthogonal basis functions and derivatives at
C.......gauss area integral points, element barycenter and vertices

C........Loop over area gauss points, then element barycenter, then vertices
         do ELT = MNEDGS,MXEDGS
            J = ELT - 2
            do Q=1,NAGP(L,J)+1+ADDGP(J)
               if (Q.eq.NAGP(L,J)+1) then
                  PT(1) = BARY(1,J)
                  PT(2) = BARY(2,J)
               elseif (Q.gt.NAGP(L,J)+1) then
                  PT(1) = VERT(Q - NAGP(L,J)-1,1,J)
                  PT(2) = VERT(Q - NAGP(L,J)-1,2,J)
               else
                  PT(1) = XAGP(Q,L,J)
                  PT(2) = YAGP(Q,L,J)
               endif
   
c...........Evaluate basis functions for specified element at current point, order basis functions hierarchically         
               call ORTHOGONAL_BASIS(J,PT,P,DIM,BASIS,DBASIS)
               if (Q.le.NAGP(L,J)) then
                  PHI_AREA(:,Q,P,J) = BASIS
                  DRPHI(:,Q,P,J) = DBASIS(:,1)
                  DSPHI(:,Q,P,J) = DBASIS(:,2)
               elseif (Q.eq.NAGP(L,J)+1) then
                  PHI_CENTER(:,P,J) = BASIS
               else
                  PHI_CORNER(:,Q - NAGP(L,J)-1,P,J) = BASIS
               endif
            enddo
         enddo
			Deallocate ( BASIS, DBASIS )
      enddo
      Deallocate( BARY, VERT )


C.....Evaluate the orthogonal basis at the edge gauss quadrature points

      do w=1,ph
         P = w
         L = w
         SZ2 = (P+2)*(P+1)**(DIM-1)/2
         Allocate ( BASIS(SZ2),DBASIS(SZ2,DIM) )

c........Calculate entries for inverse mass stiffness matrix
          
         do ELT = MNEDGS,MXEDGS
            JJJ = ELT - 2
            if (JJJ .eq. 1 .and. DIM .eq. 2) then
c........For triangular elements
               M = 1
               do J = 0,P
                  do I = 0,J
                     JJ = J - I
                     M_INV(M,P) =((2.D0*I+1.D0)*(2.D0*JJ+2.D0*I+2.D0)
     &                                                            /4.D0)
                     M = M+1
                  enddo
               enddo
            elseif (JJJ .eq. 2 .and. DIM .eq. 2) then
c........Below is for parallelogram elements; if we can determine  
c........whether an element is parallelogram, then we may use this  
c........approximation for M_INV; until we implement that capability,
C........I am leaving this commented out, -DW. 
C........l(This should not be hard to implement though, so coming soon.) 
c               M = 1
c               do J = 0,P
c                  do I = 0,J
c                     JJ = J - I
c                     M_INV(M,P) =((2.D0*I+1.D0)*(2.D0*JJ+1.D0)/4.D0)
c                     M = M+1
c                  enddo
c               enddo
C........For the general quadrilateral element, orthogonality is lost,
C........and so is our 'matrix free' implementation, so we must compute
C........the entire mass-stiffness matrix and its inverse, for each element.
               Allocate( PHI_AREAM(NAGP(P,2),(P+1)*(P+2)/2), 
     &                             WM(NAGP(P,2),NAGP(P,2)), 
     &                             PHI_AREAMT((P+1)*(P+2)/2,NAGP(P,2)),
     &                             MTEMP1((P+1)*(P+2)/2,NAGP(P,2)),
     &                             MTEMP2((P+1)*(P+2)/2,(P+1)*(P+2)/2),
     &                             MTEMPI((P+1)*(P+2)/2,(P+1)*(P+2)/2) )
               do jj = 1,NE
                  if (EDN(jj) .eq. 4) then
C.................Compute element Jacobian for each area gauss point
                     X1 = X(NM(jj,1))
                     X2 = X(NM(jj,2))
                     X3 = X(NM(jj,3))
                     X4 = X(NM(jj,4)) 
                     Y1 = Y(NM(jj,1))
                     Y2 = Y(NM(jj,2))
                     Y3 = Y(NM(jj,3))
                     Y4 = Y(NM(jj,4))
                     Aq = 0.5d0*(X1*Y2 + X2*Y3 + X3*Y4 + X4*Y1  
     &                                  - X2*Y1 - X3*Y2 - X4*Y3 - X1*Y4)
                     Atr1 = 0.5d0*((X1 - X3)*(Y2 - Y3) 
     &                                            + (X3 - X2)*(Y1 - Y3))
                     Atr2 = 0.5d0*((X1 - X4)*(Y2 - Y4) 
     &                                            + (X4 - X2)*(Y1 - Y4))
C.................Instead of storing the Jacobian, we store these,
C.................so that we may later easily compute the Jacobian for each integral.
                     Q1(jj) = 0.25d0*Aq
                     Q2(jj) = 0.25d0*(Atr2 - Atr1)
                     Q3(jj) = 0.25d0*(Aq - Atr1 - Atr2)
                     WM = 0
                     do j = 1,NAGP(P,2)
                        Jac = Q1(jj) + Q2(jj)*XAGP(j,P,2) +
     &                                                Q3(jj)*YAGP(j,P,2)
                        do k = 1,(P+1)*(P+2)/2
                           PHI_AREAM(j,k) = PHI_AREA(k,j,P,2)
                           PHI_AREAMT(k,j) = PHI_AREA(k,j,P,2)
                        enddo
                        WM(j,j) = Jac*WAGP(j,P,2)
                     enddo
                     MTEMP1 = matmul(PHI_AREAMT,WM)
                     MTEMP2 = matmul(MTEMP1,PHI_AREAM)
                     call inv(MTEMP2,MTEMPI,(P+1)*(P+2)/2)
                     do k = 1,(P+1)*(P+2)/2
                        do j = 1,(P+1)*(P+2)/2
                           M_INVM(jj,j,k,P) = MTEMPI(j,k)
                        enddo
                     enddo
                  endif
               enddo
               Deallocate(PHI_AREAM,WM,PHI_AREAMT,MTEMP1,MTEMP2,MTEMPI )
            elseif (JJJ .eq. 1 .and. DIM .eq. 3) then
c........For triangular prism elements
               M = 1 
               do k = 0,P
                  do J = 0,P
                     do I = 0,J
                        JJ = J - I
                        M_INV(M,P) =((2.D0*I+1.D0)*(2.D0*JJ+2.D0
     &                                      *I+2.D0)*(2.D0*k+1.D0)/8.D0)
                        M = M+1
                     enddo
                  enddo
               enddo
            elseif (JJJ .eq. 2 .and. DIM .eq. 3) then
C........This is only for parallelepiped elements should change this
C........later to reflect the above approach for 2D quads . 
C               M = 1 
C               do k = 0,P
C                  do J = 0,P
C                     do I = 0,J
C                        JJ = J - I                   
C                        M_INV(M,P) = ((2.D0*I+1.D0)*(2.D0*JJ+1.D0)*
C     &                                               (2.D0*k+1.D0)/8.D0)
C                        M = M+1
C                     enddo
C                  enddo
C               enddo
            endif
            
c........Loop over edges
            do 100 II = 1,ELT
c..........Loop over edge gauss points
               do Q=1,NEGP(L)+1
c..............Get edge gauss points (and midpoint) for triangle edges
                  if (JJJ .eq. 1 .and. DIM .eq. 2) then
                     if (II.eq.1) then
                        if (Q.lt.(NEGP(L)+1)) then
                           PT(1) = -XEGP(Q,L)
                           PT(2) =  XEGP(Q,L)
                        else
                           PT(1) = 0.D0
                           PT(2) = 0.D0
                        endif
                     elseif (II.eq.2) then
                        if (Q.lt.(NEGP(L)+1)) then
                           PT(1) = -1.D0
                           PT(2) = -XEGP(Q,L)
                        else
                           PT(1) = -1.D0
                           PT(2) =  0.D0
                        endif
                     elseif (II.eq.3) then
                        if (Q.lt.(NEGP(L)+1)) then
                           PT(1) = XEGP(Q,L)
                           PT(2) = -1.D0
                        else
                           PT(1) =  0.D0
                           PT(2) = -1.D0
                        endif
                     endif
c..............Get edge gauss points (and midpoint) for rectangle edges
                  elseif (JJJ .eq. 2 .and. DIM .eq. 2) then
                     if (II.eq.1) then
                        if (Q.lt.(NEGP(L)+1)) then
                           PT(1) = 1.D0
                           PT(2) = XEGP(Q,L)
                        else
                           PT(1) = 1.D0
                           PT(2) = 0.D0
                        endif
                     elseif (II.eq.2) then
                        if (Q.lt.(NEGP(L)+1)) then
                           PT(1) = -XEGP(Q,L)
                           PT(2) = 1.D0
                        else
                           PT(1) = 0.D0
                           PT(2) = 1.D0
                        endif
                     elseif (II.eq.3) then
                        if (Q.lt.(NEGP(L)+1)) then
                           PT(1) = -1.D0
                           PT(2) = -XEGP(Q,L)
                        else
                           PT(1) = -1.D0
                           PT(2) =  0.D0
                        endif
                     elseif (II.eq.4) then
                        if (Q.lt.(NEGP(L)+1)) then
                           PT(1) = XEGP(Q,L)
                           PT(2) = -1.D0
                        else
                           PT(1) =  0.D0
                           PT(2) = -1.D0
                        endif
                     endif
                  endif

c..............Evaluate basis functions for specified element at current point, order basis functions hierarchically   
                  call ORTHOGONAL_BASIS(JJJ,PT,P,DIM,BASIS,DBASIS)
                  if (Q.le.NEGP(L)+1) then
                     PHI_EDGE(:,Q,II,P,JJJ) = BASIS
                  else
                     PHI_MID(:,II,P,JJJ) = BASIS
                  endif
               enddo

               if (JJJ .eq. 1) then !.and. DIM .eq. 2) then
                  XP = PT(1)
                  YP = PT(2)
                  III = 3
                  do Q = 1,NEGP(L)+1
                     if (P.gt.1) then
                        if (Q.lt.(NEGP(L)+1)) then
                           PSI_CHECK(1,III,JJJ) = -1.D0/2.D0*(XP + YP )
                           PSI_CHECK(2,III,JJJ) =  1.D0/2.D0*(XP + 1.D0)
                           PSI_CHECK(3,III,JJJ) =  1.D0/2.D0*(YP + 1.D0)
                           III = III + 1
                        endif
                     endif
                  enddo
               elseif (JJJ .eq. 2) then ! .and. DIM .eq. 2) then
                  XP = PT(1)
                  YP = PT(2)
                  III = 4
                  do Q = 1,NEGP(L)+1
                     if (P.gt.1) then
                        if (Q.lt.(NEGP(L)+1)) then
                           PSI_CHECK(1,III,JJJ) =  1.D0/4.D0*
     &                                         (1.0D0 - XP)*(1.0D0 - YP)
                           PSI_CHECK(2,III,JJJ) =  1.D0/4.D0*
     &                                         (1.0D0 + XP)*(1.0D0 - YP)
                           PSI_CHECK(3,III,JJJ) =  1.D0/4.D0*
     &                                         (1.0D0 + XP)*(1.0D0 + YP)
                           PSI_CHECK(4,III,JJJ) =  1.D0/4.D0*
     &                                         (1.0D0 - XP)*(1.0D0 + YP)
                           III = III + 1
                        endif
                     endif
                  enddo
               endif
 100        continue
         enddo
         Deallocate ( BASIS, DBASIS )
      enddo
      
C.....Do the L2-projection of the initial conditions

      hb = 0.D0
      qx = 0.D0
      qy = 0.D0
      ze = 0.D0
      ydub = 0.d0
      hb1 = 0.D0

c$$$      do k = 1,MNE
c$$$
c$$$         
c$$$         n1 = NM(k,1)
c$$$         n2 = NM(k,2)
c$$$         n3 = NM(k,3)
c$$$                                !Define lagrange transform
c$$$
c$$$         do mm = 1,nagp(ph)     !ICs should not have higher order than ph
c$$$
c$$$            ell_1 = -0.5D0 * ( xagp(mm,ph) + yagp(mm,ph) )
c$$$            ell_2 =  0.5D0 * ( xagp(mm,ph) + 1.D0 )
c$$$            ell_3 =  0.5D0 * ( yagp(mm,ph) + 1.D0 )
c$$$
c$$$            XBCbt(k) = x(n1)*ell_1 + x(n2)*ell_2 + x(n3)*ell_3
c$$$            YBCbt(k) = y(n1)*ell_1 + y(n2)*ell_2 + y(n3)*ell_3
c$$$
c$$$
c$$$            rev =  3.141592653589793D0 / 4.D0
c$$$
c$$$            Ox = XBCbt(k)*cos(rev) - YBCbt(k)*sin(rev)
c$$$            Oy = YBCbt(k)*cos(rev) + XBCbt(k)*sin(rev)
c$$$
c$$$            radial(k) = min(sqrt( Ox**2 + (Oy + 0.25D0)**2 ), 0.18D0)/0.18D0
c$$$
c$$$
c$$$            do j = 1,dofh
c$$$
c$$$               QX(j,k,1) = QX(j,k,1) + YBCbt(k) * wagp(mm,ph) * phi_area(j,mm,ph)
c$$$
c$$$               QY(j,k,1) = QY(j,k,1) - XBCbt(k) * wagp(mm,ph) * phi_area(j,mm,ph)        
c$$$
c$$$               iota(j,k,1) = iota(j,k,1) + 0.25D0 *( 1.0D0 + cos(3.141592653589793D0*radial(k)) ) 
c$$$     &              * wagp(mm,ph) * phi_area(j,mm,ph)  
c$$$               
c$$$c$$$               iota(j,k,1) = iota(j,k,1) + 0.5D0 *( exp ( - ( (XBCbt(k)+.05D0)**2 + (YBCbt(k)+.05D0)**2  ) /0.001D0 ) ) 
c$$$c$$$     &              * wagp(mm,ph) * phi_area(j,mm,ph) 
c$$$
c$$$               if( ( sqrt((Ox + 0.25D0 )**2 + Oy**2)).le.0.18D0.and.(sqrt((Ox + 0.25D0 )**2 + 
c$$$     &              Oy**2)).ge.0.025D0.and.(Ox.le.(-0.23D0) )) then
c$$$
c$$$                  iota(j,k,1) = iota(j,k,1) + 1.0D0 * wagp(mm,ph) * phi_area(j,mm,ph) 
c$$$
c$$$               elseif(sqrt( (Ox -.25D0)**2 + Oy**2 ).le.0.18D0 ) then
c$$$
c$$$                  iota(j,k,1) =  iota(j,k,1) + (1.D0 - ( 1.D0 / 0.18D0 ) * sqrt((Ox -0.25D0)**2 + Oy**2 ) ) 
c$$$     &                 * wagp(mm,ph) * phi_area(j,mm,ph) 
c$$$
c$$$               endif
c$$$
c$$$               hb(j,k,1) = hb(j,k,1) + 1.D0*wagp(mm,ph) * phi_area(j,mm,ph)
c$$$
c$$$            enddo
c$$$
c$$$
c$$$         enddo
c$$$
c$$$                                ! get back the coeffs in each component
c$$$
c$$$         do j= 1,dofh
c$$$
c$$$            QX(j,k,1) =  QX(j,k,1) * M_inv(j,ph)
c$$$            QY(j,k,1) =  QY(j,k,1) * M_inv(j,ph)
c$$$            iota(j,k,1) =  iota(j,k,1) * M_inv(j,ph)
c$$$            hb(j,k,1) =  hb(j,k,1) * M_inv(j,ph)
c$$$
c$$$         enddo
c$$$
c$$$      enddo

      !iotaa = iota
 
      hb(1:dofh,:,1) = hbo(1:dofh,:,1)
      ze(1:dofh,:,1) = zeo(1:dofh,:,1)
      ydub(1:dofh,:,1) = ydubo(1:dofh,:)

      do chi=pl,ph
         hb1(1:dofh,:,1,chi) = hbo(1:dofh,:,1)
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
      !Set up the artificial diffusion stuff
      e1(:) = 0.D0
      balance(:) = 0.D0
      entrop(:,:) = -100.D0

      if (tune_by_hand.eq.1) then

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
         IF (EDN(J).EQ.3) THEN
            DPE_MIN(J) = MIN(DP(NM(J,1)),DP(NM(J,2)),DP(NM(J,3)))
         ELSEIF (EDN(J).EQ.4) THEN
            DPE_MIN(J) = MIN(DP(NM(J,1)),DP(NM(J,2)),DP(NM(J,3)),
     &                                                      DP(NM(J,4)))
         ENDIF
      ENDDO

C.....Compute the values of the nodal basis functions at the
C.....area gauss quadrature points, at every p level chi

      do ELT = MNEDGS,MXEDGS
         J = ELT - 2
         if (J .eq. 1) then
            do chi = 1,ph
               do I = 1,NAGP(chi,J)
                  PSI1(I,chi,J) = -1.D0/2.D0*(XAGP(I,chi,J) 
     &                                                  + YAGP(I,chi,J))
                  PSI2(I,chi,J) =  1.D0/2.D0*(XAGP(I,chi,J) + 1.D0)
                  PSI3(I,chi,J) =  1.D0/2.D0*(YAGP(I,chi,J) + 1.D0)
               enddo
            enddo
         elseif (J .eq. 2) THEN
            do chi = 1,ph
               do I = 1,NAGP(chi,J)
                  PSI1(I,chi,J) = 0.25D0*(1.D0 - XAGP(I,chi,J))*
     &                                            (1.D0 - YAGP(I,chi,J))
                  PSI2(I,chi,J) = 0.25D0*(1.D0 + XAGP(I,chi,J))*
     &                                            (1.D0 - YAGP(I,chi,J))
                  PSI3(I,chi,J) = 0.25D0*(1.D0 + XAGP(I,chi,J))*
     &                                            (1.D0 + YAGP(I,chi,J))
                  PSI4(I,chi) = 0.25D0*(1.D0 - XAGP(I,chi,J))*
     &                                            (1.D0 + YAGP(I,chi,J))
               enddo
            enddo
         endif
      enddo
      
      
C.....Store the derivatives of the (linear) nodal basis functions
      if (MXEDGS .eq. 3) then
         call ALLOC_DPSI_TRI()
      elseif (MXEDGS .eq. 4) then
         call ALLOC_DPSI_QUAD()
      endif
      do ELT = MNEDGS,MXEDGS
         if (ELT .eq. 3) then
            DRPSI(1,1,1,1) = -1.D0/2.D0
            DRPSI(2,1,1,1) =  1.D0/2.D0
            DRPSI(3,1,1,1) =  0.D0
            DSPSI(1,1,1,1) = -1.D0/2.D0
            DSPSI(2,1,1,1) =  0.D0
            DSPSI(3,1,1,1) =  1.D0/2.D0
         elseif (ELT .eq. 4) then
            do chi = 1,ph
               do I = 1,NAGP(chi,2)
                  DRPSI(1,I,chi,2) = -0.25D0*(1.D0 - YAGP(I,chi,2))
                  DRPSI(2,I,chi,2) =  0.25D0*(1.D0 - YAGP(I,chi,2))
                  DRPSI(3,I,chi,2) =  0.25D0*(1.D0 + YAGP(I,chi,2))
                  DRPSI(4,I,chi,2) = -0.25D0*(1.D0 + YAGP(I,chi,2))
                  DSPSI(1,I,chi,2) = -0.25D0*(1.D0 - XAGP(I,chi,2))
                  DSPSI(2,I,chi,2) = -0.25D0*(1.D0 + XAGP(I,chi,2))
                  DSPSI(3,I,chi,2) =  0.25D0*(1.D0 + XAGP(I,chi,2))
                  DSPSI(4,I,chi,2) =  0.25D0*(1.D0 - XAGP(I,chi,2))
               enddo
            enddo 
         endif
      enddo
      

C.....Pre-compute the derivatives of the coordinate transformation for
C.....each element

      DO J = 1,MNE

C.....Retrieve the global node numbers for the element
         IF (EDN(J).EQ.3) THEN
            N1 = NM(J,1)
            N2 = NM(J,2)
            N3 = NM(J,3)
            X1=X(N1)
            Y1=Y(N1)
            X2=X(N2)
            Y2=Y(N2)
            X3=X(N3)
            Y3=Y(N3)
            AREA = (X1 - X3)*(Y2 - Y3) + (X3 - X2)*(Y1 - Y3)
            AREA=AREA*0.5D0
C..........Compute the derivatives of the coordinate transformation

            DRDX(J) = 1.D0/AREA*(Y3 - Y1)
            DSDX(J) = 1.D0/AREA*(Y1 - Y2)

            DRDY(J) = 1.D0/AREA*(X1 - X3)
            DSDY(J) = 1.D0/AREA*(X2 - X1)
C..........Compute elemental Coriolis and friction terms

            CORI_EL(J) = (CORIF(N1) + CORIF(N2) + CORIF(N3))/3.D0
            FRIC_EL(J) = (FRIC(N1) + FRIC(N2) + FRIC(N3))/3.D0
         ELSEIF (EDN(J).EQ.4) THEN
            N1 = NM(J,1)
            N2 = NM(J,2)
            N3 = NM(J,3)
            N4 = NM(J,4)
            X1=X(N1)
            Y1=Y(N1)
            X2=X(N2)
            Y2=Y(N2)
            X3=X(N3)
            Y3=Y(N3)
            X4=X(N4)
            Y4=Y(N4)
C........Compute the derivatives of the coordinate transformation.
C........The commented expressions are for PARALLELOGRAMS ONLY; they are
C........commented out until we implement a method to identify parallelogram elements
C            DRDX(J) = 2.D0/(X2 - X1)
C            DSDX(J) = 0.D0       
C            DRDY(J) = 0.D0
C            DSDY(J) = 2.D0/(Y4 - Y1)
C........Compute derivatives for each order of rule and for each area gauss point.
            do chi = 1,ph
               do I = 1,NAGP(chi,2)
                  XX = XAGP(I,chi,2)
                  YY = YAGP(I,chi,2)
                  Jac = Q1(J) + Q2(J)*XX + Q3(J)*YY
                  DRDXQ(J,I,chi) = -0.25d0*1.0d0/Jac*((1-XX)*Y1 + 
     &                                (1+XX)*Y2 - (1+XX)*Y3 - (1-XX)*Y4)
                  DSDXQ(J,I,chi) = 0.25d0*1.0d0/Jac*((1-YY)*Y1 - 
     &                                (1-YY)*Y2 - (1+YY)*Y3 + (1+YY)*Y4)
                  DRDYQ(J,I,chi) = 0.25d0*1.0d0/Jac*((1-XX)*X1 + 
     &                                (1+XX)*X2 - (1+XX)*X3 - (1-XX)*X4)
                  DSDYQ(J,I,chi) = -0.25d0*1.0d0/Jac*((1-YY)*X1 - 
     &                                (1-YY)*X2 - (1+YY)*X3 + (1+YY)*X4)
               enddo
            end do

C..........Compute elemental Coriolis and friction terms

            CORI_EL(J)=(CORIF(N1)+CORIF(N2)+CORIF(N3)+CORIF(N4))/4.D0
            FRIC_EL(J)=(FRIC(N1) + FRIC(N2) + FRIC(N3) + FRIC(N4))/4.D0
         ENDIF

C.......Pre-compute the bathymetry and the gradient of the bathymetry at
C.......the quadrature points and compute volume of water

         DP_VOL(J,:) = 0.D0
         SFAC_ELEM(:,J,:) = 0.D0
         JJ = EDN(J) - 2
         do chi = 1,ph
            pp = (chi+1)*(chi+2)/2
            if (EDN(J) .eq. 4) then
               allocate(vtemp1(pp),vtemp2(pp), MTEMPI(pp,pp),
     &                                                    PHI_AREAV(pp))
            endif

            if (chi.gt.1) then

               do I = 1,NAGP(chi,JJ) ! Forced for wetting and drying's sake

                  BATH(I,J,chi) = 0.D0
                  DBATHDX(I,J,chi) = 0.D0
                  DBATHDY(I,J,chi) = 0.D0
                  YELEM(chi) = 0.d0

                  DO K = 1,pp 

                     if (EDN(J) .eq. 3) then
                        DPHIDX=DRPHI(K,I,chi,JJ)*DRDX(J)+
     &                                         DSPHI(K,I,chi,JJ)*DSDX(J)
                        DPHIDY=DRPHI(K,I,chi,JJ)*DRDY(J)+
     &                                         DSPHI(K,I,chi,JJ)*DSDY(J)
                        XFAC(K,I,J,chi)=M_INV(K,chi)*WAGP(I,chi,JJ)
     &                                                           *DPHIDX
                        YFAC(K,I,J,chi)=M_INV(K,chi)*WAGP(I,chi,JJ)
     &                                                           *DPHIDY
                        SRFAC(K,I,J,chi)=M_INV(K,chi)*WAGP(I,chi,JJ)*
     &                                              PHI_AREA(K,I,chi,JJ)
                     elseif (EDN(J) .eq. 4) then
                        DPHIDX=DRPHI(K,I,chi,JJ)*DRDXQ(J,I,chi)+
     &                                  DSPHI(K,I,chi,JJ)*DSDXQ(J,I,chi)
                        DPHIDY=DRPHI(K,I,chi,JJ)*DRDYQ(J,I,chi)+
     &                                  DSPHI(K,I,chi,JJ)*DSDYQ(J,I,chi)
                        vtemp1(K) = DPHIDX
                        vtemp2(K) = DPHIDY
                        PHI_AREAV(K) = PHI_AREA(K,I,chi,JJ)
                     endif
                     
                     BATH(I,J,chi) = BATH(I,J,chi) + HB(K,J,1)*
     &                                              PHI_AREA(K,I,chi,JJ)
                     YELEM(chi) = YELEM(chi) + YDUB(K,J,chi)*
     &                                              PHI_AREA(K,I,chi,JJ)
                     IF (ICS.EQ.1) THEN
                        SFAC_ELEM(I,J,chi)=1.0D0
                     ELSE
                        SFAC_ELEM(I,J,chi)=COS(SFEA0)/COS(YELEM(chi)/R)
                     ENDIF

                     DBATHDX(I,J,chi) =DBATHDX(I,J,chi)+HB(K,J,1)*DPHIDX
                     DBATHDY(I,J,chi) =DBATHDY(I,J,chi)+HB(K,J,1)*DPHIDY
                     if (EDN(J) .eq. 3) then 
                        DP_VOL(J,chi)=DP_VOL(J,chi)+WAGP(I,chi,JJ)*
     &                                    HB(K,J,1)*PHI_AREA(K,I,chi,JJ)
                     elseif (EDN(J) .eq. 4) then
                        Jac = Q1(J) + Q2(J)*XAGP(I,chi,JJ) 
     &                                            + Q3(J)*YAGP(I,chi,JJ)
                        DP_VOL(J,chi)=DP_VOL(J,chi)+Jac*WAGP(I,chi,JJ)*
     &                                    HB(K,J,1)*PHI_AREA(K,I,chi,JJ)
                     endif


                  ENDDO
                  if (EDN(J) .eq. 4) then
                     do k = 1,pp
                        do ii = 1,pp
                           Jac = Q1(J) + Q2(J)*XAGP(I,chi,JJ) +
     &                                              Q3(J)*YAGP(I,chi,JJ)
                           MTEMPI(ii,k) = Jac*M_INVM(J,ii,k,chi)
                        enddo
                     enddo
                     vtemp1 = matmul(MTEMPI,vtemp1)
                     vtemp1 = vtemp1*WAGP(I,chi,JJ)
                     vtemp2 = matmul(MTEMPI,vtemp2)
                     vtemp2 = vtemp2*WAGP(I,chi,JJ)
                     PHI_AREAV = matmul(MTEMPI,PHI_AREAV)
                     PHI_AREAV = PHI_AREAV*WAGP(I,chi,JJ)
                     do k = 1,pp
                        XFAC(k,I,J,chi) = vtemp1(k)
                        YFAC(k,I,J,chi) = vtemp2(k)
                        SRFAC(k,I,J,chi) = PHI_AREAV(k)
                     enddo
                  endif

               enddo
               if (EDN(J) .eq. 4) then
                  deallocate(vtemp1,vtemp2,MTEMPI,PHI_AREAV)
               endif
            else

               do I = 1,NAGP(chi,JJ) ! Area quadrature points

                  BATH(I,J,chi) = 0.D0
                  DBATHDX(I,J,chi) = 0.D0
                  DBATHDY(I,J,chi) = 0.D0
                  YELEM(chi) = 0.d0

                  DO K = 1,pp 

                     if (EDN(J) .eq. 3) then
                        DPHIDX=DRPHI(K,I,chi,JJ)*DRDX(J)+
     &                                         DSPHI(K,I,chi,JJ)*DSDX(J)
                        DPHIDY=DRPHI(K,I,chi,JJ)*DRDY(J)+
     &                                         DSPHI(K,I,chi,JJ)*DSDY(J)
                        XFAC(K,I,J,chi)=M_INV(K,chi)*WAGP(I,chi,JJ)
     &                                                           *DPHIDX
                        YFAC(K,I,J,chi)=M_INV(K,chi)*WAGP(I,chi,JJ)
     &                                                           *DPHIDY
                        SRFAC(K,I,J,chi)=M_INV(K,chi)*WAGP(I,chi,JJ)*
     &                                              PHI_AREA(K,I,chi,JJ)
                     elseif (EDN(J) .eq. 4) then
                        DPHIDX=DRPHI(K,I,chi,JJ)*DRDXQ(J,I,chi)+
     &                                  DSPHI(K,I,chi,JJ)*DSDXQ(J,I,chi)
                        DPHIDY=DRPHI(K,I,chi,JJ)*DRDYQ(J,I,chi)+
     &                                  DSPHI(K,I,chi,JJ)*DSDYQ(J,I,chi)
                        vtemp1(K) = DPHIDX
                        vtemp2(K) = DPHIDY
                        PHI_AREAV(K) = PHI_AREA(K,I,chi,JJ)
                     endif
                     
                     BATH(I,J,chi) = BATH(I,J,chi) + HB(K,J,1)*
     &                                              PHI_AREA(K,I,chi,JJ)
                     YELEM(chi) = YELEM(chi) + YDUB(K,J,chi)*
     &                                              PHI_AREA(K,I,chi,JJ)
                     IF (ICS.EQ.1) THEN
                        SFAC_ELEM(I,J,chi)=1.0D0
                     ELSE
                        SFAC_ELEM(I,J,chi)=COS(SFEA0)/COS(YELEM(chi)/R)
                     ENDIF

                     DBATHDX(I,J,chi) =DBATHDX(I,J,chi)+HB(K,J,1)*DPHIDX
                     DBATHDY(I,J,chi) =DBATHDY(I,J,chi)+HB(K,J,1)*DPHIDY
                     if (EDN(J) .eq. 3) then 
                        DP_VOL(J,chi)=DP_VOL(J,chi)+WAGP(I,chi,JJ)*
     &                                    HB(K,J,1)*PHI_AREA(K,I,chi,JJ)
                     elseif (EDN(J) .eq. 4) then
                        Jac = Q1(J) + Q2(J)*XAGP(I,chi,JJ) 
     &                                            + Q3(J)*YAGP(I,chi,JJ)
                        DP_VOL(J,chi)=DP_VOL(J,chi)+Jac*WAGP(I,chi,JJ)*
     &                                    HB(K,J,1)*PHI_AREA(K,I,chi,JJ)
                     endif

                  ENDDO
                  if (EDN(J) .eq. 4) then
                     do k = 1,pp
                        do ii = 1,pp
                           Jac = Q1(J) + Q2(J)*XAGP(I,chi,JJ) +
     &                                              Q3(J)*YAGP(I,chi,JJ)
                           MTEMPI(ii,k) = Jac*M_INVM(J,ii,k,chi)
                        enddo
                     enddo
                     vtemp1 = matmul(MTEMPI,vtemp1)
                     vtemp1 = vtemp1*WAGP(I,chi,JJ)
                     vtemp2 = matmul(MTEMPI,vtemp2)
                     vtemp2 = vtemp2*WAGP(I,chi,JJ)
                     PHI_AREAV = matmul(MTEMPI,PHI_AREAV)
                     PHI_AREAV = PHI_AREAV*WAGP(I,chi,JJ)
                     do k = 1,pp
                        XFAC(k,I,J,chi) = vtemp1(k)
                        YFAC(k,I,J,chi) = vtemp2(k)
                        SRFAC(k,I,J,chi) = PHI_AREAV(k)
                     enddo
                  endif

               enddo
               if (EDN(J) .eq. 4) then
                  deallocate(vtemp1,vtemp2,MTEMPI,PHI_AREAV)
               endif
            endif

         enddo 
         if (EDN(J).EQ.3) then
            do chi = 1,ph
               DP_VOL(J,chi) = 0.25D0*AREAS(J)*DP_VOL(J,chi)
            enddo
         endif

         do chi = 1,ph
            pp = (chi+1)*(chi+2)/2
            if (chi.ge.1) then

               DO L = 1,EDN(J)
                  if (EDN(J) .EQ. 4) THEN
                     allocate(MTEMP1(pp,NEGP(chi)),MTEMPI(pp,pp),
     &                                         WM(NEGP(chi),NEGP(chi)) )
                  endif
                  do I = 1,NEGP(chi) ! Forced for wetting and drying's sake

                     BATHED(I,L,J,chi) = 0.D0
                     yed(chi) = 0.d0
                     
                     DO K = 1,pp
                        
                        BATHED(I,L,J,chi) = BATHED(I,L,J,chi)+HB(K,J,1)*
     &                                            PHI_EDGE(K,I,L,chi,JJ)

                        YED(chi) = YED(chi) + YDUB(K,J,chi)*
     &                                            PHI_EDGE(K,I,L,chi,JJ)

                        IF (ICS.EQ.1) THEN
                           SFACED(I,L,J,chi)=1.0d0
                        ELSE
                           SFACED(I,L,J,chi)=COS(SFEA0)/COS(YED(chi)/R)
                        ENDIF
                        IF (EDN(J) .EQ. 3) THEN

                           EDGEQ(K,I,L,chi,JJ) = 2.0d0*M_INV(K,chi)*
     &                                PHI_EDGE(K,I,L,chi,JJ)*WEGP(I,chi)
                        ENDIF
                     ENDDO
                  enddo
                  
                  if (EDN(J) .eq. 4) then
                     WM = 0
                     do I = 1,NEGP(chi)
                        WM(I,I) = WEGP(I,chi)
                        do k = 1,pp
                           MTEMP1(k,I) = PHI_EDGE(k,I,L,chi,JJ)
                        enddo
                     enddo
                     do k = 1,pp
                        do I = 1,pp
                           MTEMPI(I,k) = M_INVM(J,I,k,chi)
                        enddo
                     enddo
                     MTEMP1 = matmul(MTEMP1,WM)
                     MTEMP1 = matmul(MTEMPI,MTEMP1)
                     do I = 1,NEGP(chi)
                        do k = 1,pp
                           EDGEQ(k,I,L,chi,JJ) = MTEMP1(k,I)
                        enddo
                     enddo
                     deallocate(MTEMP1,MTEMPI,WM)
                  endif
                  
               ENDDO

            else

               DO L = 1,EDN(J)
                  if (EDN(J) .EQ. 4) THEN
                     allocate(MTEMP1(pp,NEGP(chi)),MTEMPI(pp,pp),
     &                                         WM(NEGP(chi),NEGP(chi)) )
                  endif
                  do I = 1,NEGP(chi) ! Edge quadrature points

                     BATHED(I,L,J,chi) = 0.D0
                     yed(chi) = 0.d0
                     
                     DO K = 1,pp
                        
                        BATHED(I,L,J,chi) = BATHED(I,L,J,chi)+HB(K,J,1)*
     &                                            PHI_EDGE(K,I,L,chi,JJ)

                        YED(chi) = YED(chi) + YDUB(K,J,chi)*
     &                                            PHI_EDGE(K,I,L,chi,JJ)

                        IF (ICS.EQ.1) THEN
                           SFACED(I,L,J,chi)=1.0d0
                        ELSE
                           SFACED(I,L,J,chi)=COS(SFEA0)/COS(YED(chi)/R)
                        ENDIF
                        IF (EDN(J) .EQ. 3) THEN

                           EDGEQ(K,I,L,chi,JJ) = 2.0d0*M_INV(K,chi)*
     &                                PHI_EDGE(K,I,L,chi,JJ)*WEGP(I,chi)
                        ENDIF
                     ENDDO
                  enddo
                  
                  if (EDN(J) .eq. 4) then
                     WM = 0
                     do I = 1,NEGP(chi)
                        WM(I,I) = WEGP(I,chi)
                        do k = 1,pp
                           MTEMP1(k,I) = PHI_EDGE(k,I,L,chi,JJ)
                        enddo
                     enddo
                     do k = 1,pp
                        do I = 1,pp
                           MTEMPI(I,k) = M_INVM(J,I,k,chi)
                        enddo
                     enddo
                     MTEMP1 = matmul(MTEMPI,MTEMP1)
                     MTEMP1 = matmul(MTEMP1,WM)
                     do I = 1,NEGP(chi)
                        do k = 1,pp
                           EDGEQ(k,I,L,chi,JJ) = MTEMP1(k,I)
                        enddo
                     enddo
                     deallocate(MTEMP1,MTEMPI,WM)
                  endif
                  
               ENDDO

            endif

         enddo

C........Store bathymetry at triangular vertices and edge gauss points for wet-dry
C........


         do chi = 1,ph

            DO I = 1,EDN(J)
               DP_NODE(I,J,chi) = DP(NM(J,I))
            ENDDO


            IF (NCHECK(chi,EDN(J)-2).GT.EDN(J)) THEN
               II = EDN(J) + 1
               DO L = 1,EDN(J)
                  DO I = 1,NEGP(chi)
                     DP_NODE(II,J,chi) = BATHED(I,L,J,chi)
                     II = II + 1
                  ENDDO
               ENDDO
            ENDIF
         enddo
      ENDDO
      
      DO ELT = MNEDGS,MXEDGS
        J = ELT - 2
        IF (J .eq. 1) THEN
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
            PSI_CHECK(1,I,1) = -1.D0/2.D0*(XI + YI)
            PSI_CHECK(2,I,1) =  1.D0/2.D0*(XI + 1.D0)
            PSI_CHECK(3,I,1) =  1.D0/2.D0*(YI + 1.D0)
            do chi = 1,ph
              DO K = 1,(chi+1)*(chi+2)/2 
                PHI_CHECK(K,I,chi,1) = PHI_CORNER(K,I,chi,1)
              ENDDO
            enddo
          ENDDO
        ELSEIF (J .eq. 2) THEN
          DO I = 1,4
            IF (I.EQ.1) THEN
              XI = -1.D0
              YI = -1.D0
            ELSEIF (I.EQ.2) THEN
              XI =  1.D0 
              YI = -1.D0
            ELSEIF (I.EQ.3) THEN
              XI =  1.D0
              YI =  1.D0
            ELSE
              XI = -1.0D0
              YI =  1.0D0
            ENDIF
            PSI_CHECK(1,I,2) =  0.25d0*(1.0D0 - XP)*(1.0D0 - YP)                        
            PSI_CHECK(2,I,2) =  0.25d0*(1.0D0 + XP)*(1.0D0 - YP)                        
            PSI_CHECK(3,I,2) =  0.25d0*(1.0D0 + XP)*(1.0D0 + YP)                       
            PSI_CHECK(4,I,2) =  0.25d0*(1.0D0 - XP)*(1.0D0 + YP)                   
            do chi = 1,ph
              DO K = 1,(chi+1)*(chi+2)/2 
                PHI_CHECK(K,I,chi,2) = PHI_CORNER(K,I,chi,2)
              ENDDO
            enddo
          ENDDO
        ENDIF

        do chi =1,ph
           IF (NCHECK(chi,J).GT.ELT) THEN
              II = ELT + 1
              DO L = 1,ELT
                 DO I = 1,NEGP(chi)
                    DO K = 1,(chi+1)*(chi+2)/2
                       
                       PHI_CHECK(K,II,chi,J) = PHI_EDGE(K,I,L,chi,J)
                       
                    ENDDO
                    II = II + 1
                 ENDDO
              ENDDO
           ENDIF
        enddo
      
C.....Integrate the basis functions

        PHI_INTEGRATED = 0.D0
        do chi = 1,ph
           DO I = 1,NAGP(chi,J)
              DO K = 1,(chi+1)*(chi+2)/2
                 PHI_INTEGRATED(K,chi,J) =PHI_INTEGRATED(K,chi,J)+
     &                                 WAGP(I,chi,J)*PHI_AREA(K,I,chi,J)
              enddo
           ENDDO
         ENDDO
      ENDDO

C.....Wetting and drying is not turned on

      IF (NOLIFA.EQ.0 .OR. NOLIFA.EQ.1) THEN
         DO J = 1,MNE
            WDFLG(J) = 1
                                !DOFS(J) = 3
         ENDDO
         
C.....Wetting and drying is turned on but there are no dry nodes below
C.....geoid

      ELSEIF (NOLIFA.EQ.2 .AND. NSTARTDRY.EQ.0) THEN

         DO J = 1,MNE
            
C.........Check to see if there are initially any dry nodes

                                !ZE1 = 0.D0
                                !ZE2 = 0.D0
                                !ZE3 = 0.D0 

            ZE1 = ze(1,J,1)
            ZE2 = ze(1,J,1)
            ZE3 = ze(1,J,1)
            
            IF (DP(NM(J,1)).LT.H0) ZE1 = max(ze1,H0 - DP(NM(J,1)))
            IF (DP(NM(J,2)).LT.H0) ZE2 = max(ze2,H0 - DP(NM(J,2)))
            IF (DP(NM(J,3)).LT.H0) ZE3 = max(ze3,H0 - DP(NM(J,3)))
            
            IF (EDN(J) .EQ. 4) THEN
               ZE4 = ze(1,J,1)
               IF (DP(NM(J,4)).LT.H0) ZE4 = max(ze4,H0 - DP(NM(J,4)))
            ENDIF
C.........If so set initial surface elevation values
            
            IF (EDN(J) .EQ. 3) THEN
               IF ((ZE1 + ZE2 + ZE3)/3.D0 .NE. ze(1,J,1)) THEN
                  IF (p_0.EQ.0) THEN
                     DP_MIN = MIN(DP(NM(J,1)),DP(NM(J,2)),DP(NM(J,3)))
                     ze(1,J,1) = max(ze(1,j,1),H0 - DP_MIN)
                  ELSE
                     IF (ze(1,J,1).GT.(ZE1+ZE2+ZE3)/3.d0) THEN
                        IF (dof_0.NE.1) THEN
                           ze(2,J,1)=0.d0
                           ze(3,J,1)=0.d0
                           ze(4:dofh,J,1) = 0.D0 ! forced again for transparency
                        ENDIF
                     ELSE
                        ze(1,J,1)=(ZE1+ZE2+ZE3)/3.D0
                        IF (DOF_0.NE.1) THEN
                           ze(2,J,1)=-1.D0/6.D0*(ZE1+ZE2)+1.D0/3.D0*ZE3
                           ze(3,J,1) = -0.5D0*ZE1 + 0.5D0*ZE2
                           ze(4:dofh,J,1) = 0.D0 ! forced again for transparency
                        ENDIF
                     ENDIF
                  ENDIF
               
C...........   .and set wet/dry flag (0 = dry, 1 = wet)
               
                  WDFLG(J) = 0
                  DOFS(J) = dof_0
               ELSE
                  WDFLG(J) = 1
                  DOFS(J) = dofs(J)
               ENDIF
            ELSEIF (EDN(J) .EQ. 4) THEN
               IF ((ZE1 + ZE2 + ZE3 + ZE4)/4.D0 .NE. ze(1,J,1)) THEN
                  IF (p_0.EQ.0) THEN
                     DP_MIN = MIN(DP(NM(J,1)),DP(NM(J,2)),DP(NM(J,3)),
     &                                                      DP(NM(J,4)))
                     ze(1,J,1) = max(ze(1,j,1),H0 - DP_MIN)
                  ELSE
                     IF (ze(1,J,1).GT.(ZE1+ZE2+ZE3+ZE4)/4.d0) THEN
                        IF (dof_0.NE.1) THEN
                           ze(2,J,1)=0.d0 
                           ze(3,J,1)=0.d0
                           ze(4:dofh,J,1) = 0.D0 ! forced again for transparency
                        ENDIF
                     ELSE
                        ze(1,J,1)=(ZE1+ZE2+ZE3+ZE4)/4.D0
                        IF (DOF_0.NE.1) THEN
                           ! COME BACK TO THIS!
                           ze(2,J,1) = -1.D0/6.D0*(ZE1 + ZE2)+1.D0/3.D0*ZE3
                           ze(3,J,1) = -0.5D0*ZE1 + 0.5D0*ZE2
                           ze(4:dofh,J,1) = 0.D0 ! forced again for transparency
                        ENDIF
                     ENDIF
                  ENDIF
               
C...........   .and set wet/dry flag (0 = dry, 1 = wet)
               
                  WDFLG(J) = 0
                  DOFS(J) = dof_0
               ELSE
                  WDFLG(J) = 1
                  DOFS(J) = dofs(J)
               ENDIF
            ENDIF
         ENDDO


C.....If there are dry nodes below geoid

      ELSEIF (NOLIFA.EQ.2 .AND. NSTARTDRY.EQ.1) THEN
         
C.......Loop over elements
         
         DO J = 1,NE

C........Retrieve global node numbers for element

            N1 = NM(J,1)
            N2 = NM(J,2)
            N3 = NM(J,3)

C........Check to see if nodes are initially dry

            ZE1 = 0.0d0
            ZE2 = 0.0d0
            ZE3 = 0.0d0
            IF (STARTDRY(N1).EQ.1) ZE1 = H0 - DP(N1)
            IF (DP(N1).LT.H0) ZE1 = H0 - DP(N1)
            IF (STARTDRY(N2).EQ.1) ZE2 = H0 - DP(N2)
            IF (DP(N2).LT.H0) ZE2 = H0 - DP(N2)
            IF (STARTDRY(N3).EQ.1) ZE3 = H0 - DP(N3)
            IF (DP(N3).LT.H0) ZE3 = H0 - DP(N3)

            IF (MODAL_IC.EQ.3) THEN
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
            IF (EDN(J) .EQ. 4) THEN
               N4 = NM(J,4)
               ZE4 = 0.0d0
               IF (STARTDRY(N4).EQ.1) ZE4 = H0 - DP(N4)
               IF (DP(N4).LT.H0) ZE4 = H0 - DP(N4)
               IF (MODAL_IC.EQ.3) THEN
                  IF (STARTDRY(N4).EQ.-88888) THEN
                     ZE4 = H0 - DP(N4)
                  ELSE
                     ZE4 = STARTDRY(N4)
                  ENDIF
               ENDIF
            ENDIF

C.........If so set initial surface elevation values
            IF (EDN(J) .EQ. 3) THEN
               IF ((ZE1 + ZE2 + ZE3).NE.0.0d0) THEN
                  IF (P_0.EQ.0) THEN
                     DP_MIN = MIN(DP(NM(J,1)),DP(NM(J,2)),DP(NM(J,3)))
                     ze(1,J,1) = max(ze(1,j,1),H0 - DP_MIN)
                  ELSE
                     IF (ze(1,J,1).GT.(ZE1+ZE2+ZE3)/3.d0) THEN
                        IF (DOF_0.NE.1) THEN
                           ze(2,J,1)=0.d0
                           ze(3,J,1)=0.d0
                           ze(4:dofh,J,1) = 0.D0 ! forced again for transparency
                        ENDIF
                     ELSE
                        ze(1,J,1)=(ZE1+ZE2+ZE3)/3.D0
                        IF (DOF_0.NE.1) THEN
                           ze(2,J,1)=-1.D0/6.D0*(ZE1+ZE2)+1.D0/3.D0*ZE3
                           ze(3,J,1) = -0.5D0*ZE1 + 0.5D0*ZE2
                           ze(4:dofh,J,1) = 0.D0 ! forced again for transparency
                        ENDIF
                     ENDIF
                  ENDIF
               
C...........   Set wet/dry flag (0 = dry, 1 = wet)
               
                  WDFLG(J) = 0
                  DOFS(J) = DOF_0
               ELSE
                  WDFLG(J) = 1
                  DOFS(J) = DOFS(J)
               ENDIF
            ELSEIF (EDN(J) .EQ. 4) THEN
               IF ((ZE1 + ZE2 + ZE3 + ZE4).NE.0.0d0) THEN
                  IF (P_0.EQ.0) THEN
                     DP_MIN = MIN(DP(NM(J,1)),DP(NM(J,2)),DP(NM(J,3)),
     &                                                      DP(NM(J,4)))
                     ze(1,J,1) = max(ze(1,j,1),H0 - DP_MIN)
                  ELSE
                     IF (ze(1,J,1).GT.(ZE1+ZE2+ZE3+ZE4)/4.d0) THEN
                        IF (DOF_0.NE.1) THEN
                           ze(2,J,1)=0.d0
                           ze(3,J,1)=0.d0
                           ze(4:dofh,J,1) = 0.D0 ! forced again for transparency
                        ENDIF
                     ELSE
                        ze(1,J,1)=(ZE1+ZE2+ZE3+ZE4)/4.D0
                        IF (DOF_0.NE.1) THEN
                           ! COME BACK TO THIS!
                           ze(2,J,1) = -1.D0/6.D0*(ZE1 + ZE2)+1.D0/3.D0*ZE3
                           ze(3,J,1) = -0.5D0*ZE1 + 0.5D0*ZE2
                           ze(4:dofh,J,1) = 0.D0 ! forced again for transparency
                        ENDIF
                     ENDIF
                  ENDIF
               
C...........   Set wet/dry flag (0 = dry, 1 = wet)
               
                  WDFLG(J) = 0
                  DOFS(J) = DOF_0
               ELSE
                  WDFLG(J) = 1
                  DOFS(J) = DOFS(J)
               ENDIF
            ENDIF
         ENDDO
      ENDIF

      IF(MYPROC.EQ.0) THEN
         print *, 'DONE'
         print *, ''
      ENDIF

C.....Read in modal dof for initial conditions
!Asserts error if you project onto lower order basis
!That is, do not expect convergence

      IF (MODAL_IC.EQ.1) THEN
         OPEN(163,FILE=DIRNAME//'/'//'Initial_Conditions.163')
         OPEN(164,FILE=DIRNAME//'/'//'Initial_Conditions.164')
         OPEN(114,FILE=DIRNAME//'/'//'Initial_Bathymetry.114')
         READ(163,*) P_READ
         IF (P_READ.NE.ph) THEN
            PRINT*,'INCONSISTENCY IN P -- CHECK INPUT FILES'
            STOP
         ENDIF
         READ(164,*) P_READ,P_READ2
         IF (P_READ.NE.ph) THEN
            PRINT*,'INCONSISTENCY IN P -- CHECK INPUT FILES'
            STOP
         ENDIF
         IF (P_READ2.NE.ph) THEN
            PRINT*,'INCONSISTENCY IN P -- CHECK INPUT FILES'
            STOP
         ENDIF
         READ(114,*) P_READ
         IF (P_READ.NE.ph) THEN
            PRINT*,'INCONSISTENCY IN P -- CHECK INPUT FILES'
            STOP
         ENDIF
         DO J = 1,MNE
            DO K = 1,DOFH
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

C.....Read in modal dof for hot start conditions

      IF (MODAL_IC.EQ.2) THEN
         OPEN(263,FILE=DIRNAME//'/'//'Hot_start.263')
         OPEN(264,FILE=DIRNAME//'/'//'Hot_start.264')
         OPEN(214,FILE=DIRNAME//'/'//'Hot_start.214')
#ifdef TRACE
         OPEN(288,FILE=DIRNAME//'/'//'Hot_start.288')
#endif
#ifdef CHEM
         OPEN(289,FILE=DIRNAME//'/'//'Hot_start.289')
#endif
#ifdef DYNP
         OPEN(291,FILE=DIRNAME//'/'//'Hot_start.291')
#endif
c$$$#ifdef SED_LAY
c$$$         OPEN(290,FILE=DIRNAME//'/'//'Hot_start.290')
c$$$#endif
         READ(263,*) P_READ
         IF (P_READ.NE.PH) THEN
            PRINT*,'INCONSISTENCY IN P -- CHECK INPUT FILES'
            STOP
         ENDIF
         READ(264,*) P_READ,P_READ2
         IF (P_READ.NE.PH) THEN
            PRINT*,'INCONSISTENCY IN P -- CHECK INPUT FILES'
            STOP
         ENDIF
         IF (P_READ2.NE.PH) THEN
            PRINT*,'INCONSISTENCY IN P -- CHECK INPUT FILES'
            STOP
         ENDIF
         READ(214,*) P_READ
         IF (P_READ.NE.PH) THEN
            PRINT*,'INCONSISTENCY IN P -- CHECK INPUT FILES'
            STOP
         ENDIF
         DO J = 1,MNE
            DO K = 1,DOFH
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
c$$$#ifdef SED_LAY
c$$$               do ll=1,layers
c$$$                  READ(290,*) bed(K,J,1,ll)
c$$$               enddo
c$$$#endif
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
c$$$#ifdef SED_LAY
c$$$         CLOSE(290)
c$$$#endif
#ifdef DYNP
         CLOSE(291)
#endif
      ENDIF
      
C.....Initialize the DG.63 output file

      IF (ABS(NOUTGE).EQ.1) THEN
         OPEN(631,FILE=DIRNAME//'/'//'DG.63')
         WRITE(631,3220) RUNDES, RUNID, AGRID
         WRITE(631,3645) NDSETSE, dofh, DTDP*NSPOOLGE, NSPOOLGE, 1
      ENDIF

C.....Initialize the DG.64 output file

      IF (ABS(NOUTGV).EQ.1) THEN
         OPEN(641,FILE=DIRNAME//'/'//'DG.64')
         WRITE(641,3220) RUNDES, RUNID, AGRID
         WRITE(641,3645) NDSETSV, dofh, DTDP*NSPOOLGV, NSPOOLGV, 2
      ENDIF

C.....Initialize the DG.65 output file (contains elemental statuses such
C.....as the wet/dry status.

      IF ((ABS(NOUTGE).EQ.1).AND.(NOLIFA.GE.2)) THEN
         OPEN(651,FILE=DIRNAME//'/'//'DG.65')
         WRITE(651,3220) RUNDES, RUNID, AGRID
         WRITE(651,3645) NDSETSE, dofh, DTDP*NSPOOLGE, NSPOOLGE, 1
      ENDIF
 3220 FORMAT(1X,A32,2X,A24,2X,A24)
 3645 FORMAT(1X,I10,1X,I10,1X,E15.7,1X,I5,1X,I5)

C.....Set p back to original value if p = 0

      IF (P_0.NE.pl) THEN
         PDG_EL(:) = 0
         PDG(:) = 0
         DOF = 1
         DOFL = 1
         DOFS(:) = 1
         DOF_0 = 1
         pl = 0
                                !p = 0
         pa = 0
      ENDIF
      
C.....Compute basis functions at stations

      IF (NSTAE.GT.0) THEN      ! Elevation stations
         CALL ALLOC_STAE( NSTAE )
         DO I = 1,NSTAE
           CALL STA_BASIS( DIM, XEL(I), YEL(I),  NNE(I), PHI_STAE(:,I) )
         ENDDO
      ENDIF
      
      IF (NSTAV.GT.0) THEN      ! Velocity Stations
         CALL ALLOC_STAV( NSTAV )
         DO I = 1,NSTAV
           CALL STA_BASIS( DIM, XEV(I), YEV(I),  NNV(I), PHI_STAV(:,I) )
         ENDDO
      ENDIF

C.....Prep the slopelimiter

      IF (SLOPEFLAG.NE.0) THEN
         IF(MYPROC.EQ.0)THEN
            print *, 'Slope limiting prep begins, "kshanti"'
         ENDIF
         CALL ALLOC_SLOPELIM()
         CALL PREP_SLOPELIM()
         IF(MYPROC.EQ.0)THEN
            print *, 'Finished'
         ENDIF
      ENDIF
C--   

C.....Close files

!      CLOSE(115)
!      CLOSE(25)
      
      RETURN
      END SUBROUTINE PREP_DG

C***********************************************************************
C     
C     SUBROUTINE RK_TIME()
C     
C     This subroutine does preparatory stuff for DG
C     
C     Written by Ethan Kubatko (03-07-2005)
C     
C***********************************************************************

      SUBROUTINE RK_TIME()

      USE GLOBAL
      USE DG

      IMPLICIT NONE

      INTEGER L,i,j,k
      REAL(SZ) ARK, BRK, CASUM, MAX_BOA
      Real(SZ) eps_const,RKC_omega0,RKC_omega1

C.....Allocate the time stepping arrays

      NRK = RK_STAGE
      CALL ALLOC_RK()

#ifdef RKSSP
C.....The forward Euler method

      IF ((RK_STAGE.EQ.1).AND.(RK_ORDER.EQ.1)) THEN

         ATVD(:,:) = 0.D0
         BTVD(:,:) = 0.D0
         CTVD(:,:) = 0.D0
         DTVD(:)   = 0.D0

         ATVD(1,1) = 1.D0
         BTVD(1,1) = 1.D0

C.....SSP(s,2) schemes
C.....THESE ARE DG-OPTIMIZED METHODS!
C..... Kubatko, Ethan J., Benjamin A. Yeager, and David I. Ketcheson. 
C....."Optimal strong-stability-preserving Runge–Kutta time 
C..... discretizations for discontinuous Galerkin methods." 
C..... Journal of Scientific Computing 60.2 (2014): 313-344.
      ELSEIF (RK_ORDER.EQ.2) THEN
      
        ATVD(:,:) = 0.D0
        BTVD(:,:) = 0.D0
        CTVD(:,:) = 0.D0
        DTVD(:)   = 0.D0
        
        IF (RK_STAGE.EQ.3) THEN

          ATVD(1,1) =  1.000000000000000D0
          ATVD(2,1) =  0.087353119859156D0
          ATVD(3,1) =  0.344956917166841D0
          ATVD(2,2) =  0.912646880140844D0
          ATVD(3,3) =  0.655043082833159D0

          BTVD(1,1) =  0.528005024856522D0
          BTVD(3,1) =  0.022826837460491D0
          BTVD(2,2) =  0.481882138633993D0
          BTVD(3,3) =  0.345866039233415D0

        ELSEIF (RK_STAGE.EQ.4) THEN

          ATVD(1,1) =  1.000000000000000D0
          ATVD(2,1) =  0.394806441339829D0
          ATVD(3,1) =  0.002797307087390D0
          ATVD(4,1) =  0.252860909354373D0
          ATVD(2,2) =  0.605193558660171D0
          ATVD(3,3) =  0.997202692912610D0
          ATVD(4,4) =  0.747139090645627D0

          BTVD(1,1) =  0.406584463657504D0
          BTVD(3,1) =  0.013637216641451D0
          BTVD(4,1) =  0.016453567333598D0
          BTVD(2,2) =  0.246062298456822D0
          BTVD(3,3) =  0.405447122055692D0
          BTVD(4,4) =  0.303775146447707D0

        ELSEIF (RK_STAGE.EQ.5) THEN

          ATVD(1,1) = 1.000000000000000D0
          ATVD(2,1) = 0.235593265061659D0
          ATVD(3,1) = 0.174017972351526D0
          ATVD(4,1) = 0.235264368870758D0
          ATVD(5,1) = 0.141720372339803D0
          ATVD(2,2) = 0.764406734938341D0
          ATVD(4,2) = 0.000058643383967D0
          ATVD(5,2) = 0.095374613155521D0
          ATVD(3,3) = 0.825982027648475D0
          ATVD(5,3) = 0.000311763705780D0
          ATVD(4,4) = 0.764676987745275D0
          ATVD(5,5) = 0.762593250798895D0

          BTVD(1,1) = 0.324840618151514D0
          BTVD(3,1) = 0.108822380501601D0
          BTVD(4,1) = 0.054392262422093D0
          BTVD(5,1) = 0.000000180291569D0
          BTVD(2,2) = 0.248310356296551D0
          BTVD(4,2) = 0.000019049753098D0
          BTVD(5,2) = 0.030981548293401D0
          BTVD(3,3) = 0.268312512443371D0
          BTVD(5,3) = 0.000101273514903D0
          BTVD(4,4) = 0.248398145385413D0
          BTVD(5,5) = 0.247721262987686D0
        ELSE


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
        ENDIF
C.....SSP(3,3) scheme

      ELSEIF ((RK_STAGE.EQ.3).AND.(RK_ORDER.EQ.3)) THEN

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

C.....SSP(4,3) scheme

      ELSEIF ((RK_STAGE.EQ.4).AND.(RK_ORDER.EQ.3)) THEN

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

C.....SSP(5,3) scheme

      ELSEIF ((RK_STAGE.EQ.5).AND.(RK_ORDER.EQ.3)) THEN

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

C.....SSP(6,3) scheme

      ELSEIF ((RK_STAGE.EQ.6).AND.(RK_ORDER.EQ.3)) THEN

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

C.....SSP(7,3) scheme

      ELSEIF ((RK_STAGE.EQ.7).AND.(RK_ORDER.EQ.3)) THEN

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

C.....SSP(8,3) scheme

      ELSEIF ((RK_STAGE.EQ.8).AND.(RK_ORDER.EQ.3)) THEN

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

C.....SSP(5,4) scheme

      ELSEIF ((RK_STAGE.EQ.5).AND.(RK_ORDER.EQ.4)) THEN

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

C.....SSP(6,4) scheme

      ELSEIF ((RK_STAGE.EQ.6).AND.(RK_ORDER.EQ.4)) THEN

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

C.....SSP(7,4) scheme

      ELSEIF ((RK_STAGE.EQ.7).AND.(RK_ORDER.EQ.4)) THEN

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

C.....SSP(8,4) scheme

      ELSEIF ((RK_STAGE.EQ.8).AND.(RK_ORDER.EQ.4)) THEN

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

C.....Compute the time dependent parameters

      DO K = 0,RK_STAGE-1
         DO I = 1,RK_STAGE
            CASUM = 0.D0
            DO L = K+1,I-1
               CASUM = CASUM + CTVD(L,K+1)*ATVD(I,L+1)
            ENDDO
            CTVD(I,K+1) = BTVD(I,K+1) + CASUM
         ENDDO
      ENDDO

      DO K = 1,RK_STAGE-1
         DTVD(K+1) = 0.D0
         DO L = 0,K-1
            DTVD(K+1) = DTVD(K+1) + CTVD(K,L+1)
         ENDDO
      ENDDO

C.....Compute the maximum beta over alpha ratio at each stage

      DO IRK = 1,RK_STAGE
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

C-----------------------------------------------------------
C.... Compute the Runge-Kutta Chebyshev (RKC) version

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
               RKC_Tdprime(j) = (j*(j*RKC_T(j) - RKC_omega0*RKC_U(j-1))/
     &            ( RKC_omega0**(2.D0) - 1.D0 ) )
               RKC_b(j) = RKC_Tdprime(j)*RKC_Tprime(j)**(-2.D0)

            endif
            
            RKC_a(j)  = 1.D0 - RKC_b(j)*RKC_T(j)
            
         enddo

         if(RK_stage.gt.1) then

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

         RKC_c(RK_stage) = 0.D0
         RKC_c(RK_stage) = 1.D0

#endif


      DEALLOCATE(CTVD)
      RETURN
      END SUBROUTINE RK_TIME

C.....Subroutine to find the inverse of a square matrix by Guass-Jordan elimination

      subroutine Inv(matrix, inverse, n)
      Use sizes, only : sz

      implicit none
      integer n
      real(sz), dimension(n,n) :: matrix
      real(sz), dimension(n,n) :: inverse
      
      integer :: i, j, k, l
      real(sz) :: m
      real(sz), dimension(n,2*n) :: augmatrix !augmented matrix
      
                                !Augment input matrix with an identity matrix

      do i = 1, n

         do j = 1, 2*n

            if ( j.le.n ) then

               augmatrix(i,j) = matrix(i,j)

            else if ((i+n) == j) then

               augmatrix(i,j) = 1

            else

               augmatrix(i,j) = 0

            endif

         enddo

      enddo
      
                                !Reduce augmented matrix to upper traingular form

      do k =1, n-1

         if (augmatrix(k,k) == 0) then


            do i = k+1, n

               if (augmatrix(i,k) /= 0) then

                  do j = 1,2*n

                     augmatrix(k,j) = augmatrix(k,j)+augmatrix(i,j)

                  enddo

               endif


            enddo

         endif

         do j = k+1, n

            m = augmatrix(j,k)/augmatrix(k,k)

            do i = k, 2*n

               augmatrix(j,i) = augmatrix(j,i) - m*augmatrix(k,i)

            enddo

         enddo

      enddo 
      
                                !Test for invertibility

      do i = 1, n

         if (augmatrix(i,i) == 0) then

            inverse = 0

            return

         endif

      enddo
      
                                !Make diagonal elements as 1

      do i = 1 , n

         m = augmatrix(i,i)

         do j = i, (2 * n)

            augmatrix(i,j) = (augmatrix(i,j) / m)

         enddo

      enddo
      
                                !Reduced right side half of augmented matrix to identity matrix

      do k = n-1, 1, -1

         do i =1, k

            m = augmatrix(i,k+1)

            do j = k, (2*n)

               augmatrix(i,j) = augmatrix(i,j) -augmatrix(k+1,j) * m

            enddo

         enddo

      enddo
      
                                !Compute answer

      do i =1, n

         do j = 1, n

            inverse(i,j) = augmatrix(i,j+n)

         enddo

      enddo

      end subroutine Inv
