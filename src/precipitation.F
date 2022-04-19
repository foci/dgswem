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
      
      IF (RAINFALL.EQ.3) THEN
         N1 = NM(II,1)
         N2 = NM(II,2)
         N3 = NM(II,3)
         SOURCE_R1 =  1.D0/3.D0*(PREC2(N1)
     &                       + PREC2(N2) + PREC2(N3))
      ENDIF
      
      ! Linearly interpolate fort.prec file
      IF (RAINFALL.EQ.1 .OR. RAINFALL.EQ.2) THEN
      DO I = 1, NTPREC-1
         IF (ITDT.ge.ITPREC(I) .AND. ITDT.le.ITPREC(I+1)) THEN
            SOURCE_R1 = (PREC(II,I+1) - PREC(II,I))/
     &       (ITPREC(I+1) - ITPREC(I))
     &                     *(ITDT - ITPREC(I)) + PREC(II,I)
            EXIT
         ENDIF
      ENDDO
      ENDIF
      
      if (RAINFALL .EQ. 1) then ! no infiltration
         
      elseif (RAINFALL .EQ. 2) then ! apply Green-Ampt model
         if (ITDT.gt.ITPREC(NTPREC)) then
            SOURCE_R1 = 0.D0
         else
            I = I + 1
            DT1 = DTVD(IRK)*DTDP
          
            SOURCE_R1 = PREC(II,I)
            if (SOURCE_R1*DT1.gt.IntDp(ii)) then
               SOURCE_R1 = SOURCE_R1 - IntDp(II)/DT1;
               IntDp(II) = 0.D0;
               I1 = (Ks(II)*DT1 - 2.D0*TotInf(II));
               I2 = (Ks(II)*TotInf(II) + Ks(II)*Hf(II)*Md(II));
               infilt = (1.D0/(2.D0*DT1))*
     &          (I1 + sqrt(I1**2 + 8.D0*I2*DT1));
             
               if (SOURCE_R1 .le. infilt) then
                   infilt = SOURCE_R1;
                   SOURCE_R1 = 0.D0;
               else
                   SOURCE_R1 = SOURCE_R1 - infilt;
               endif
               TotInf(II) = TotInf(II)  + infilt*DT1;
            else
               IntDp(II) = IntDp(II) - SOURCE_R1*DT1;
               SOURCE_R1 = 0.D0;
            endif
         endif
      endif
      
      ! stop
      ! DO I = 1,NP
        ! PREC(I) = rainfall_rate
      ! ENDDO
      
      
      
        
      RETURN
      END


