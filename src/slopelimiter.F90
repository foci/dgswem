#ifdef SLOPE5

!***********************************************************************
!     
!     SUBROUTINE SLOPELIMITER5()
!     
!     Written by Clint Dawson - 30 June 2010
!     01-10-2011 - cem - adapted for p_enrichment and multicomponent
!     
!     
!***********************************************************************

      SUBROUTINE SLOPELIMITER5_PARTA(s,dg_here,global_here)

!.....Use appropriate modules

      USE SIZES
      USE GLOBAL
      USE DG

#ifdef CMPI
      USE MESSENGER
#endif

      IMPLICIT NONE

      type (sizes_type) :: s
      type (dg_type) :: dg_here
      type (global_type) :: global_here

!.....Declare local variables

      INTEGER L, LL, INC1,INC2,INC3,KDP,NN,IVAR,I,J,kk,k,varnum,bb,varnum_prev
      REAL(SZ) ZEC(3),ZEVERTEX(3),DIF(3),SUMLOC,SUMDIF,SIGNDIF
      REAL(SZ)    DIV,REDFAC,REDMAX,tmp1,tmp2,tmp3,bound
      Real(SZ) ZEMIN1(3),ZEMAX1(3),QXMIN1(3),QXMAX1(3)
      Real(SZ) QYMIN1(3),QYMAX1(3)
      Real(SZ) iotaMIN1(3),iotaMAX1(3)
      Real(SZ) iota2MIN1(3),iota2MAX1(3)
      Real(SZ), pointer:: arraymin(:),arraymax(:)

!     FIND THE MAXIMUM AND MINIMUM OF EACH VARIABLE OVER ALL ELEMENTS 
!     SHARING A NODE

      bound = 0.0D0

      DO I = 1,global_here%NP
         ZE_MIN1(I)=99999.
         ZE_MAX1(I)=-99999.
         QX_MIN1(I)=99999.
         QX_MAX1(I)=-99999.
         QY_MIN1(I)=99999.
         QY_MAX1(I)=-99999.


         global_here%NO_NBORS = global_here%EL_COUNT(I)

         DO J = 1,global_here%NO_NBORS
            global_here%NBOR_EL = global_here%ELETAB(I,1+J)

!     IF(dg_here%WDFLG(global_here%NBOR_EL).EQ.0) CYCLE ! DON'T COUNT DRY ELEMENTS  sb 02/26/07

            global_here%ZE_DG(J) = dg_here%ZE(1,global_here%NBOR_EL,dg_here%IRK+1)
            global_here%QX_DG(J) = dg_here%QX(1,global_here%NBOR_EL,dg_here%IRK+1)
            global_here%QY_DG(J) = dg_here%QY(1,global_here%NBOR_EL,dg_here%IRK+1)



!     
            IF (global_here%ZE_DG(J).LT.ZE_MIN1(I))THEN
               ZE_MIN1(I)=global_here%ZE_DG(J)
            ENDIF
            IF (global_here%ZE_DG(J).GT.ZE_MAX1(I)) THEN
               ZE_MAX1(I)=global_here%ZE_DG(J)
            ENDIF
            IF (global_here%QX_DG(J).LT.QX_MIN1(I))THEN
               QX_MIN1(I)=global_here%QX_DG(J)
            ENDIF
            IF (global_here%QX_DG(J).GT.QX_MAX1(I)) THEN
               QX_MAX1(I)=global_here%QX_DG(J)
            ENDIF
            IF (global_here%QY_DG(J).LT.QY_MIN1(I))THEN
               QY_MIN1(I)=global_here%QY_DG(J)
            ENDIF
            IF (global_here%QY_DG(J).GT.QY_MAX1(I)) THEN
               QY_MAX1(I)=global_here%QY_DG(J)
            ENDIF


         ENDDO
      ENDDO
! #ifdef CMPI
! 
!       CALL UPDATER(ZE_MIN1,ZE_MAX1,QX_MIN1,3)
!       CALL UPDATER(QX_MAX1,QY_MIN1,QY_MAX1,3)
! 
! #endif

      RETURN
      END SUBROUTINE 

      SUBROUTINE SLOPELIMITER5_PARTB(s,dg_here,global_here)

!.....Use appropriate modules

      USE SIZES
      USE GLOBAL
      USE DG

#ifdef CMPI
      USE MESSENGER
#endif

      IMPLICIT NONE

      type (sizes_type) :: s
      type (dg_type) :: dg_here
      type (global_type) :: global_here

!.....Declare local variables

      INTEGER L, LL, INC1,INC2,INC3,KDP,NN,IVAR,I,J,kk,k,varnum,bb,varnum_prev
      REAL(SZ) ZEC(3),ZEVERTEX(3),DIF(3),SUMLOC,SUMDIF,SIGNDIF
      REAL(SZ)    DIV,REDFAC,REDMAX,tmp1,tmp2,tmp3,bound
      Real(SZ) ZEMIN1(3),ZEMAX1(3),QXMIN1(3),QXMAX1(3)
      Real(SZ) QYMIN1(3),QYMAX1(3)
      Real(SZ) iotaMIN1(3),iotaMAX1(3)
      Real(SZ) iota2MIN1(3),iota2MAX1(3)
      Real(SZ), pointer:: arraymin(:),arraymax(:)


!     
!     LOOP OVER ELEMENTS TO CALCULATE NEW VERTEX VALUES
!     

      bb = 1

      DO I=1,global_here%NE 
         !IF(dg_here%WDFLG(I).EQ.0) CYCLE ! DON'T COUNT DRY ELEMENTS  sb 02/26/07
         global_here%N1=global_here%NM(I,1)
         global_here%N2=global_here%NM(I,2)
         global_here%N3=global_here%NM(I,3) 
         
         varnum = 3



         DO IVAR=1,varnum

            IF (IVAR.EQ.1) THEN
               ZEC(1)=dg_here%ZE(1,I,dg_here%IRK+1)
               ZEC(2)=dg_here%ZE(2,I,dg_here%IRK+1)
               ZEC(3)=dg_here%ZE(3,I,dg_here%IRK+1)
               ZEMAX1(1)=ZE_MAX1(global_here%N1)
               ZEMIN1(1)=ZE_MIN1(global_here%N1)
               ZEMAX1(2)=ZE_MAX1(global_here%N2)
               ZEMIN1(2)=ZE_MIN1(global_here%N2)
               ZEMAX1(3)=ZE_MAX1(global_here%N3)
               ZEMIN1(3)=ZE_MIN1(global_here%N3)
            ENDIF

            IF (IVAR.EQ.2) THEN
               ZEC(1)=dg_here%QX(1,I,dg_here%IRK+1)
               ZEC(2)=dg_here%QX(2,I,dg_here%IRK+1)
               ZEC(3)=dg_here%QX(3,I,dg_here%IRK+1)
               ZEMAX1(1)=QX_MAX1(global_here%N1)
               ZEMIN1(1)=QX_MIN1(global_here%N1)
               ZEMAX1(2)=QX_MAX1(global_here%N2)
               ZEMIN1(2)=QX_MIN1(global_here%N2)
               ZEMAX1(3)=QX_MAX1(global_here%N3)
               ZEMIN1(3)=QX_MIN1(global_here%N3)
            ENDIF

            IF (IVAR.EQ.3) THEN
               ZEC(1)=dg_here%QY(1,I,dg_here%IRK+1)
               ZEC(2)=dg_here%QY(2,I,dg_here%IRK+1)
               ZEC(3)=dg_here%QY(3,I,dg_here%IRK+1)
               ZEMAX1(1)=QY_MAX1(global_here%N1)
               ZEMIN1(1)=QY_MIN1(global_here%N1)
               ZEMAX1(2)=QY_MAX1(global_here%N2)
               ZEMIN1(2)=QY_MIN1(global_here%N2)
               ZEMAX1(3)=QY_MAX1(global_here%N3)
               ZEMIN1(3)=QY_MIN1(global_here%N3)
            ENDIF

            

!     COMPUTE THE VERTEX VALUES

            ZEVERTEX(1)=ZEC(1)
            ZEVERTEX(2)=ZEC(1)
            ZEVERTEX(3)=ZEC(1)
            DO KK=2,3
               ZEVERTEX(1)=ZEVERTEX(1)+ dg_here%PHI_CORNER(KK,1,1)*ZEC(KK)
               ZEVERTEX(2)=ZEVERTEX(2)+ dg_here%PHI_CORNER(KK,2,1)*ZEC(KK)
               ZEVERTEX(3)=ZEVERTEX(3)+ dg_here%PHI_CORNER(KK,3,1)*ZEC(KK)
            ENDDO

            
!     RESET THE VERTEX VALUE TO BE LESS THAN OR EQUAL TO THE MAX AND
!     GREATER THAN OR EQUAL TO THE MIN AT THAT VERTEX
!     
            ZEVERTEX(1)=DMAX1(DMIN1(ZEVERTEX(1),ZEMAX1(1)),ZEMIN1(1))
            ZEVERTEX(2)=DMAX1(DMIN1(ZEVERTEX(2),ZEMAX1(2)),ZEMIN1(2))
            ZEVERTEX(3)=DMAX1(DMIN1(ZEVERTEX(3),ZEMAX1(3)),ZEMIN1(3))

            tmp1 = ZEVERTEX(1)
            tmp2 = ZEVERTEX(2)
            tmp3 = ZEVERTEX(3)


!     LOOP OVER THE VERTICES 3 TIMES
!     IF THE VALUE AT THE VERTEX IS ABOVE (BELOW) THE MAX (MIN) AT THAT
!     VERTEX THEN SUBTRACT OFF THE DIFFERENCE AND ADD IT TO THE OTHER
!     VERTICES
!     
            DO LL=1,3
               SUMLOC=(ZEVERTEX(1)+ZEVERTEX(2)+ZEVERTEX(3))/3.0D0
               SUMDIF=(SUMLOC-ZEC(1))*3.0D0
               SIGNDIF=DSIGN(1.D0,SUMDIF)
               DIF(1)=(ZEVERTEX(1)-ZEC(1))*SIGNDIF
               DIF(2)=(ZEVERTEX(2)-ZEC(1))*SIGNDIF
               DIF(3)=(ZEVERTEX(3)-ZEC(1))*SIGNDIF
               INC1=0
               IF (DIF(1).GT.0) INC1=1
               INC2=0
               IF (DIF(2).GT.0) INC2=1
               INC3=0
               IF (DIF(3).GT.0) INC3=1
               KDP=INC1+INC2+INC3
!     
               DO K=1,3
                  DIV=DMAX1(1.D0,DFLOAT(KDP))
                  IF (DIF(K).GT.0) THEN
                     REDFAC=SUMDIF*SIGNDIF/DIV
                     KDP=KDP-1
                  ELSE
                     REDFAC=0
                  ENDIF
                  IF (SIGNDIF.GT.0) THEN
                     REDMAX=ZEVERTEX(K)-ZEMIN1(K)
                  ELSE
                     REDMAX=ZEMAX1(K)-ZEVERTEX(K)
                  ENDIF
                  REDFAC=DMIN1(REDFAC,REDMAX)
                  SUMDIF=SUMDIF-REDFAC*SIGNDIF
                  ZEVERTEX(K)=ZEVERTEX(K)-REDFAC*SIGNDIF
               ENDDO
            ENDDO
            IF (IVAR.EQ.1) THEN


               dg_here%ZE(2,I,dg_here%IRK+1)=-1.d0/6.d0*(ZEVERTEX(1)+ZEVERTEX(2))&
                   +1.d0/3.d0*ZEVERTEX(3)
               dg_here%ZE(3,I,dg_here%IRK+1)=-.5d0*ZEVERTEX(1)+.5d0*ZEVERTEX(2)
            ENDIF
            IF (IVAR.EQ.2) THEN


               dg_here%QX(2,I,dg_here%IRK+1)=-1.d0/6.d0*(ZEVERTEX(1)+ZEVERTEX(2))&
              +1.d0/3.d0*ZEVERTEX(3)
               dg_here%QX(3,I,dg_here%IRK+1)=-.5d0*ZEVERTEX(1)+.5d0*ZEVERTEX(2)
            ENDIF
            IF (IVAR.EQ.3) THEN


               dg_here%QY(2,I,dg_here%IRK+1)=-1.d0/6.d0*(ZEVERTEX(1)+ZEVERTEX(2))&
                   +1.d0/3.d0*ZEVERTEX(3)
               dg_here%QY(3,I,dg_here%IRK+1)=-.5d0*ZEVERTEX(1)+.5d0*ZEVERTEX(2)
            ENDIF



      ENDDO
      RETURN
      END SUBROUTINE 


#endif
