C******************************************************************************
C PADCIRC VERSION 45.01 09/30/2004                                            *
C  last changes in this file VERSION 44.15                                    *
C                                                                             *
C-----------------------------------------------------------------------------*   
C                                                                             *
C  Modifications for parallel runs, Shintaro Bunya, Aug 2005                  *
C                                                                             *
C******************************************************************************   
C                                                                             *
C  GLOBAL DATA FOR 3D VS ROUTINES                                             *
C                                                                             *
C****************************************************************************** 
C 
      MODULE GLOBAL_3DVS

C...
C...BRING IN NECESSARY EXTERNAL MODE (2DDI) DATA
C...
Csb-PDG1 added some of variables in padcirc45
      USE GLOBAL,
     &ONLY : SZ, NBYTE, MNP, MNE, MNEI, MNVEL,
     &       MYPROC,LNAME,DIRNAME,
     &       X, Y, DP, NNEIGH, NEITAB, AREAS,
     &       WSX1, WSY1, WSX2, WSY2,
     &       BSX => BSX1, BSY => BSY1,
     &       UU => UU2, VV => VV2,
     &       DUU => DUU1, DUV => DUV1, DVV => DVV1,
     &       ETA1, ETA2, CORIF,
     &       LBCODEI, CSII, SIII, QNORMSP1 => QN2,
     &       NP, NOLICA, NOLIFA, NSCREEN, IHOT, IHSFIL, IHOTSTP, NWS,
     &       RHOWAT0, G,
     &       RUNDES, RDES4, RDES8, RUNID, RID4, RID8, AGRID, AID4, AID8

      USE DIFF45_41,
     &ONLY : MNODES, NEITABELE, DASigT, LBARRAY_POINTER, SIGT0,
     &       CBaroclinic,
     &       VIDBCPDX => VIDBCPDX1, VIDBCPDY => VIDBCPDY1,
     &       BTP => MOM_LV_X

      PUBLIC
C--

C...
C...DESCRIPTION OF EXTERNAL MODE DATA
C...
C   INTEGER  NSCREEN              : flag to suppress or allow screen output
C   INTEGER  NEITAB(MNP,MNEI)     : table of neighbor nodes for each node
C   INTEGER  NNEIGH(MNP)          : number of neighboring nodes for each node
C   INTEGER  NEITABELE(MNP,MNei)  : table of neighboring elements for each node 
C   INTEGER  NOLIFA               : nonlinear finite amplitude flag (1=yes,0=no)
C   INTEGER  NOLICA               : nonlinear advection flag (1=yes, 0=no)
C   INTEGER  NP                   : number of horizontal nodes
C   INTEGER  LBARRAY_POINTER(MNP) : pointer into array of land or flux boundary
C   INTEGER  LBCODEI(MNVel)       : array of land or flux boundary codes
C   REAL TAUSX1(MNP),TAUSY1(MNP)     : Wind stress components at time level s
C   REAL TAUSX2(MNP),TAUSY2(MNP)     : Wind stress components at time level s+1
C   REAL CORIF(MNP)                  : nodal values of Coriolis parameter
C   REAL EVM(MNP)                    : lateral eddy viscosity for momentum
C   REAL VIDBCPDX(MNP),VIDBCPDY(MNP) : (X,Y) derivatives of vertically integrated pressure grad
C   REAL UBAR(MNP),VBAR(MNP)         : vertically averaged velocity components 
C   REAL DUU(MNP),DVV(MNP),DUV(MNP)  : velocity dispersion terms
C   REAL TAUBX(MNP),TAUBY(MNP)       : bottom stress computed after velocity solution
C   REAL QNORMSP1(MNVel)             : specified normal flux boundary condition at time level s+1
C   REAL CSII(MNVel),SIII(MNVel)     : cosine and sine of normal flux boundary node
C   REAL ETA1(MNP),ETA2(MNP)         : water surf elev at time levels s, s+1
C   REAL DP(MNP)                     : still water depth
C   REAL BTP(MNP)                    : barotropic pressure (incl TP & wl) at time levels s+1/2
C   REAL*8  X(MNP),Y(MNP)            : nodal coordinates
C   REAL*8  AREAS(MNE)               : 2*Element Area

C...
C...DECLARE INTERNAL MODE GLOBAL ARRAYS
C...
      COMPLEX,ALLOCATABLE :: GAMMA(:), Q(:,:)

      REAL(SZ),ALLOCATABLE :: SIGMA(:), EVTOT(:)
      REAL(SZ),ALLOCATABLE :: INM(:,:), LVN(:)
      REAL(SZ),ALLOCATABLE :: WZ(:,:)
      REAL(SZ),ALLOCATABLE :: SIGT(:,:), TEMP(:,:), SAL(:,:), BCP(:,:)
      REAL(SZ),ALLOCATABLE :: Q20(:,:), L(:,:)

      INTEGER,ALLOCATABLE :: ISDHOUT(:), ISVHOUT(:), ISTHOUT(:)

C...
C...  DECLARE INTERNAL MODE GLOBAL SCALARS
C...
      COMPLEX :: I, IDTALP1, IDT1MALP1

      REAL(SZ) :: A, B, AMB, GORHO, GORHOOAMB
      REAL(SZ) :: KP, EVMIN, EVCON, Z0S
      REAL(SZ) :: Z0B, DTALP3, DT1MALP3, DTALP2, DT1MALP2
      REAL(SZ) :: THETA1, THETA2

      INTEGER :: NFEN, IDIAG, ISLIP, IEVC
      INTEGER :: ISTART, IDEN
      INTEGER :: I3DSD, NTSSSD, NTSFSD, NSPO3DSD, NSSD, NHN3DSD, ISDREC
      INTEGER :: I3DSV, NTSSSV, NTSFSV, NSPO3DSV, NSSV, NHN3DSV, ISVREC
      INTEGER :: I3DST, NTSSST, NTSFST, NSPO3DST, NSST, NHN3DST, ISTREC
      INTEGER :: I3DGD, NTSSGD, NTSFGD, NSPO3DGD, NSGD, IGDREC
      INTEGER :: I3DGV, NTSSGV, NTSFGV, NSPO3DGV, NSGV, IGVREC
      INTEGER :: I3DGT, NTSSGT, NTSFGT, NSPO3DGT, NSGT, IGTREC

      logical :: turb_allocated = .false.

C...
C...DESCRIPTION OF INTERNAL MODE GLOBAL DATA
C...
C   ISTART                            : counter of # of timesteps VSSOL has been called
C   IDIAG                             : flag to specify amount of diagnostic output
C   NFEN                              : number of vertical nodes
C   ISLIP                             : slip coefficient flag (0=no slip, 1=linear, 2=quadratic)
C   IDEN                              : flag, IDEN=0: barotropic run, IDEN=1: baroclinic run
C   NTSSSD,NTSFSD,NSPO3DSD,NSSD       : parameters controlling station density output
C   I3DSD,NHN3DSD,ISDHOUT(MNP),ISDREC : parameters controlling station density output
C   NTSSSV,NTSFSV,NSPO3DSV,NSSV       : parameters controlling station velocity output
C   I3DSV,NHN3DSV,ISVHOUT(MNP),ISVREC : parameters controlling station velocity output
C   NTSSST,NTSFST,NSPO3DST,NSST       : parameters controlling station velocity output
C   I3DST,NHN3DST,ISTHOUT(MNP),ISTREC : parameters controlling station turbulence output
C   NTSSGD,NTSFGD,NSPO3DGD,NSGD       : parameters controlling global density output
C   I3DGD,IGDREC                      : parameters controlling global density output
C   NTSSGV,NTSFGV,NSPO3DGV,NSGV       : parameters controlling global velocity output     
C   I3DGV,IGVREC                      : parameters controlling global velocity output
C   NTSSGT,NTSFGT,NSPO3DGT,NSGT       : parameters controlling global turbulence output   
C   I3DGT,IGTREC                      : parameters controlling global turbulence output
C   KP                                : input bottom friction coefficient
C   DTALP2,DT1MALP2                   : DelT*alpha2, DelT*(1-alpha2)
C   DTALP3,DT1MALP3                   : DelT*alpha3, DelT*(1-alpha3)
C   WZ(MNP,NFEN)                    : "z" vertical velocity
C   A,B                             : top and bottom sigma values
C   AMB                             : (a-b) difference between top and bottom sigma values
C   SIGMA(NFEN)                     : Sigma levels of vertical nodes
C   Q20(MNP,NFEN)                   : turbulent kinetic energy computed by MY closure
C   L(MNP,NFEN)                     : turbulent length scale computed by MY closure
C   EVTOT(NFEN)                     : vertical eddy viscosity     
C   BCP(MNP,NFEN)                   : baroclinic pressure (integrated down from surface)
C   GORhoOAMB                       :       gravity/(reference density)/(a-b)
C   SIGT(MNP,NFEN)                  : sigma T values
C   TEMP(MNP,NFEN)                  : temperature values
C   SAL(MNP,NFEN)                   : salinity values
C   Inm(NFEN,3)                     : Integral used in vertical assembly
C   LVn(NFEN)                       : Integral used in vertical assembly
C   Q(MNP,NFEN)                     : horizontal velocity in the complex form u + iv
C   GAMMA(NFEN)                     : horizontal velocity soln in the complex form u + iv
C   IDTALP1,IDT1MALP1                 : i*DelT*alpha1, i*DelT*(1-alpha1)
C   I                                 : square root of (-1)


C-------------------end of data declarations----------------------------------C


      CONTAINS

      SUBROUTINE ALLOC_3DVS()
C
C     Allocate space for arrays used in 3D VS routines
C
      ALLOCATE( SIGMA(NFEN), EVTOT(NFEN) )
      ALLOCATE( GAMMA(NFEN), INM(NFEN,3), LVN(NFEN) )
      ALLOCATE( Q(MNP,NFEN), WZ(MNP,NFEN) )
      ALLOCATE( SIGT(MNP,NFEN), TEMP(MNP,NFEN) )
      ALLOCATE( SAL(MNP,NFEN), BCP(MNP,NFEN) )
      ALLOCATE( Q20(MNP,NFEN), L(MNP,NFEN) )
      ALLOCATE( ISDHOUT(MNP), ISVHOUT(MNP), ISTHOUT(MNP) )

      RETURN
      END SUBROUTINE ALLOC_3DVS

      END MODULE GLOBAL_3DVS
