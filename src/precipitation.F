C***********************************************************************
C
C     SUBROUTINE PRECIPITATION( )
C
C     This subroutine reads the precipitation file and interpolate 
C     over the mesh
C
C     Current version generate only uniform precipitation
C
C     Written by Younghun Kang
C
C***********************************************************************

      SUBROUTINE PRECIPITATION(ITDT,II,SOURCE_R1)

C.....Use appropriate modules

      USE GLOBAL
      USE DG
      
      IMPLICIT NONE
      
C.....Declare local variables
      INTEGER I, PrecFlag
      INTEGER, INTENT(IN) :: II ! global index of element
      REAL(SZ), INTENT(IN) :: ITDT
      REAL(SZ), INTENT(OUT) :: SOURCE_R1
      REAL(SZ) :: DT1, I1, I2, infilt
      REAL(SZ) :: SOURCE_R11,SOURCE_R22,infilt1,infilt2,TotInf2
      
      ! Linearly interpolate fort.21 file
      IF (RAINFALL.EQ.1 .OR. RAINFALL.EQ.2) THEN
         DO I = 1, NTPREC-1
            IF (ITDT.ge.ITPREC(I) .AND. ITDT.le.ITPREC(I+1)) THEN
               SOURCE_R1 = (PREC(II,I+1) - PREC(II,I))/
     &          (ITPREC(I+1) - ITPREC(I))
     &                        *(ITDT - ITPREC(I)) + PREC(II,I)
               EXIT
            ENDIF
         ENDDO
      ELSE ! Compute element-averaged rainfall for parametric models
         N1 = NM(II,1)
         N2 = NM(II,2)
         N3 = NM(II,3)
         SOURCE_R1 =  1.D0/3.D0*(PREC2(N1)
     &                       + PREC2(N2) + PREC2(N3))
      ENDIF
      
      ! Linearly interpolate fort.prec file
      IF (RAINFALL.EQ.1 .OR. RAINFALL.EQ.2) THEN
      SOURCE_R1 = 0.d0
      DO I = 1, NTPREC-1
         IF (TIMEDG.ge.ITPREC(I) .AND. TIMEDG.le.ITPREC(I+1)) THEN
            SOURCE_R1 = (PREC(II,I+1) - PREC(II,I))/
     &       (ITPREC(I+1) - ITPREC(I))
     &                     *(TIMEDG - ITPREC(I)) + PREC(II,I)
            EXIT
         ENDIF
      ENDDO
      
      SOURCE_R11 = 0.d0
      DO I = 1, NTPREC-1
         IF (ITDT.ge.ITPREC(I) .AND. ITDT.le.ITPREC(I+1)) THEN
            SOURCE_R11 = (PREC(II,I+1) - PREC(II,I))/
     &       (ITPREC(I+1) - ITPREC(I))
     &                     *(ITDT - ITPREC(I)) + PREC(II,I)
            EXIT
         ENDIF
      ENDDO
      
      SOURCE_R22 = 0.d0
      DO I = 1, NTPREC-1
         IF (ITDT+DTDP.ge.ITPREC(I) .AND. ITDT+DTDP.le.ITPREC(I+1)) THEN
            SOURCE_R22 = (PREC(II,I+1) - PREC(II,I))/
     &       (ITPREC(I+1) - ITPREC(I))
     &                     *(ITDT+DTDP - ITPREC(I)) + PREC(II,I)
            EXIT
         ENDIF
      ENDDO
      ENDIF
      
      if (RAINFALL .EQ. 1) then ! no infiltration
         SOURCE_R1 = SOURCE_R11 + (SOURCE_R22 - SOURCE_R11)*
     &    (TIMEDG - ITDT)/DTDP
         
      elseif (RAINFALL .EQ. 2) then ! apply Green-Ampt model
         SOURCE_R1 = SOURCE_R11 + (SOURCE_R22 - SOURCE_R11)*
     &    (TIMEDG - ITDT)/DTDP
         ! I = I + 1
         DT1 = DTVD(IRK)*DTDP
         ! IF (IRK.EQ.1) THEN
            ! DT1 = DTVD(IRK)*DTDP
         ! ELSE
            ! DT1 = (DTVD(IRK) - DTVD(IRK-1))*DTDP
         ! ENDIF
         ! DT1 = DTDP/NRK
         
         I1 = (Ks(II)*DTDP - 2.D0*TotInf(II));
         I2 = (Ks(II)*TotInf(II) + Ks(II)*Hf(II)*Md(II));
         infilt1 = (1.D0/(2.D0*DTDP))*
     &    (I1 + sqrt(I1**2 + 8.D0*I2*DTDP));
         infilt1 = Ks(II)*(1+Hf(II)*Md(II)/TotInf(II))
         if (SOURCE_R11 .le. infilt1) then
             infilt1 = SOURCE_R11;
         endif
         
         TotInf2 = TotInf(II)  + infilt1*DTDP;
         I1 = (Ks(II)*DTDP - 2.D0*TotInf2);
         I2 = (Ks(II)*TotInf2 + Ks(II)*Hf(II)*Md(II));
         infilt2 = (1.D0/(2.D0*DTDP))*
     &    (I1 + sqrt(I1**2 + 8.D0*I2*DTDP));
         infilt2 = Ks(II)*(1+Hf(II)*Md(II)/TotInf2)
         if (SOURCE_R22 .le. infilt2) then
             infilt2 = SOURCE_R22;
         endif
         
         infilt = infilt1 + (infilt2 - infilt1)*(TIMEDG - ITDT)/DTDP
         if (SOURCE_R1 .le. infilt) then
             infilt = SOURCE_R1;
             SOURCE_R1 = 0.D0;
         else
             SOURCE_R1 = SOURCE_R1 - infilt;
         endif
         
         if (irk.eq.nrk) then
         TotInf(II) = TotInf(II)  + infilt1*DTDP;
         endif
         
         !! SOURCE_R1 = PREC(II,I)
         !if (SOURCE_R1*DT1.gt.IntDp(ii)) then
         !   SOURCE_R1 = SOURCE_R1 - IntDp(II)/DT1;
         !   IntDp(II) = 0.D0;
         !   I1 = (Ks(II)*DT1 - 2.D0*TotInf(II));
         !   I2 = (Ks(II)*TotInf(II) + Ks(II)*Hf(II)*Md(II));
         !   infilt = (1.D0/(2.D0*DT1))*
     &   !    (I1 + sqrt(I1**2 + 8.D0*I2*DT1));
         !   infilt = Ks(II)*(1+Hf(II)*Md(II)/TotInf(II))
         !   if (SOURCE_R1 .le. infilt) then
         !       infilt = SOURCE_R1;
         !       SOURCE_R1 = 0.D0;
         !   else
         !       SOURCE_R1 = SOURCE_R1 - infilt;
         !   endif
         !   TotInf(II) = TotInf(II)  + infilt*DT1;
         !else
         !   IntDp(II) = IntDp(II) - SOURCE_R1*DT1;
         !   SOURCE_R1 = 0.D0;
         !endif
      endif
      
      ! stop
      ! DO I = 1,NP
        ! PREC(I) = rainfall_rate
      ! ENDDO
      
      
      
        
      RETURN
      END


