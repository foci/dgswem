!***********************************************************************
!     
!     SUBROUTINE:  WRITE_RESULTS
!     
!     Taken from the timestep subroutine.  Modified to include
!     interpolation of DG modal degrees of freedom to the nodes.  These
!     multiple nodal values are then averaged to a single value.
!     
!     Aug ??, 2005, sb, Modifications for parallel runs
!     Jan 01, 2007, sb, Files are forced to be written if FORCE_WRITE = .TRUE.
!     
!***********************************************************************

      SUBROUTINE WRITE_RESULTS(s,dg,IT,FORCE_WRITE)

!.....Use appropriate modules

      USE GLOBAL
      USE DG
#ifdef HARM
      USE HARM
#endif
      use sizes

#ifdef CMPI
      USE MESSENGER_ELEM 
#endif
      
      type (sizes_type) :: s
      type (dg_type) :: dg

!.....Declare local variables

      INTEGER IT,I,K,Maxp(0:8),Minp(0:8),J,KK
      REAL(SZ) AREA, DEPTH, ANGLE_SUM, FH_NL, DP1, DP2, DP3, DP00, ZE00
      LOGICAL FORCE_WRITE
      real(sz) qmaxe,elmaxe,tempx,tempy,rev,Ox,Oy,iota_error,iota_int
      real(sz) iota2_int,iota3_int,iota2_error,iota3_error
      integer imaxze,imaxq, ModetoNode

      Real(SZ),allocatable :: XBCbt(:),YBCbt(:),radial(:),temper(:)

      Allocate ( XBCbt(s%MNE),YBCbt(s%MNE),radial(s%MNE) )
      Allocate ( temper(s%mne) )
      
!.....Transform from modal coordinates to nodal coordinates and average
!.....to single nodal values 

      ModetoNode = 0

      if (ModetoNode.eq.1) then

         DO I = 1,S%MNP

            NO_NBORS = EL_COUNT(I)

            AREA_SUM = 0
            ANGLE_SUM = 0
            CEN_SUM = 0
            DO 333 J = 1,NO_NBORS
               NBOR_EL = ELETAB(I,1+J)

               IF(dg%WDFLG(NBOR_EL).EQ.0) CYCLE ! DON'T COUNT DRY ELEMENTS  sb 02/26/07

               DO K = 1,3
                  N1 = NM(NBOR_EL,K)
                  IF (N1.EQ.I) THEN
                     ZE_DG(J) = dg%ZE(1,NBOR_EL,1)
                     QX_DG(J) = dg%QX(1,NBOR_EL,1)
                     QY_DG(J) = dg%QY(1,NBOR_EL,1)
                     HB_DG(J) = dg%HB(1,NBOR_EL,1)

#ifdef TRACE
                     iota_DG(J) = dg%iota(1,NBOR_EL,1)
                     iotaa_DG(J) = dg%iotaa(1,NBOR_EL,1)
#endif

#ifdef CHEM
                     iota_DG(J) = dg%iota(1,NBOR_EL,1)
                     iota2_DG(J) = dg%iota2(1,NBOR_EL,1)
#endif

#ifdef DYNP
                     dynP_DG(J) = dg%dynP(1,NBOR_EL,1)
#endif

#ifdef SED_LAY
                     HB_DG(J) = 0.d0
                     bed_DG(J,:) = dg%bed(1,NBOR_EL,1,:)
                     HB_DG(J) = sum(bed_DG(J,:))
#endif


                     DO KK = 2,dg%DOFH
                        ZE_DG(J) = ZE_DG(J) + dg%PHI_CORNER(KK,K,dg%ph)*dg%ZE(KK,NBOR_EL,1)
                        QX_DG(J) = QX_DG(J) + dg%PHI_CORNER(KK,K,dg%ph)*dg%QX(KK,NBOR_EL,1)
                        QY_DG(J) = QY_DG(J) + dg%PHI_CORNER(KK,K,dg%ph)*dg%QY(KK,NBOR_EL,1)

#ifdef TRACE
                        iota_DG(J) = iota_DG(J) + dg%PHI_CORNER(KK,K,dg%ph)*dg%iota(KK,NBOR_EL,1)
                        iotaa_DG(J) = iotaa_DG(J) + dg%PHI_CORNER(KK,K,dg%ph)*dg%iotaa(KK,NBOR_EL,1)
#endif

#ifdef CHEM
                        iota_DG(J) = iota_DG(J) + dg%PHI_CORNER(KK,K,dg%ph)*dg%iota(KK,NBOR_EL,1)
                        iota2_DG(J) = iota2_DG(J) + dg%PHI_CORNER(KK,K,dg%ph)*dg%iota2(KK,NBOR_EL,1)
#endif

#ifdef DYNP
                        dynP_DG(J) = dynP_DG(J) + dg%PHI_CORNER(KK,K,dg%ph)*dg%dynP(KK,NBOR_EL,1)
#endif

#ifdef SED_LAY
                        do l=1,s%layers
                           bed_DG(J,l) = bed_DG(J,l) + dg%PHI_CORNER(KK,K,dg%ph)*dg%bed(KK,NBOR_EL,1,l)
                           HB_DG(J) = HB_DG(J) + dg%PHI_CORNER(KK,K,dg%ph)*dg%bed(KK,NBOR_EL,1,l)
                           iotaa_DG(J) = dg%iotaa(1,NBOR_EL,1)
                        enddo
#else

                        HB_DG(J) = HB_DG(J) + dg%PHI_CORNER(KK,K,dg%ph)*dg%HB(KK,NBOR_EL,1)
#endif
                        
                     ENDDO
                     AREA = 0.5*AREAS(NBOR_EL)
                     AREA_SUM = AREA_SUM + AREA
                     GOTO 333
                  ENDIF
               ENDDO
 333        CONTINUE

            ETA2(I) = 0.D0
            tracer(I) = 0.D0
            tracer2(I) = 0.D0
            UU2(I)  = 0.D0
            VV2(I)  = 0.D0

            IF(SEDFLAG.GE.2) DP(I) = 0.D0

            DO J = 1,NO_NBORS
               NBOR_EL = ELETAB(I,1+J)

               IF(dg%WDFLG(NBOR_EL).EQ.0) CYCLE ! DON'T COUNT DRY ELEMENTS  sb 02/26/07

               AREA = 0.5*AREAS(NBOR_EL)/AREA_SUM
               DEPTH = ZE_DG(J) + HB_DG(J)
               FH_NL = 1.D0/(NLEQ*DEPTH + LEQ)
               ETA2(I) = ETA2(I) + AREA*ZE_DG(J)

#ifdef TRACE
               tracer(I)  = tracer(I)  + AREA*iota_DG(J)*FH_NL
               tracer2(I) = tracer2(I) + AREA*iotaa_DG(J)*FH_NL
#endif

#ifdef CHEM
               tracer(I) = tracer(I) + AREA*iota_DG(J)*FH_NL
               tracer2(I) = tracer2(I) + AREA*iota2_DG(J)*FH_NL
#endif 

#ifdef DYNP
               dyn_P(I) = dyn_P(I) + AREA*dynP_DG(J)*FH_NL
#endif    
               
               if (etamax(i).lt.eta2(i)) etamax(i)=eta2(i)
               UU2(I)  = UU2(I)  + AREA*QX_DG(J)*FH_NL
               VV2(I)  = VV2(I)  + AREA*QY_DG(J)*FH_NL
               IF(SEDFLAG.GE.2) DP(I) = DP(I) + (AREA/AREA_SUM)*HB_DG(J)
            ENDDO
         ENDDO

      else

         ETA2 = 0.D0
         tracer = 0.D0
         tracer2 = 0.D0
         UU2  = 0.D0
         VV2  = 0.D0
         DPe = 0.D0

         do j=1,ne

            dg%pa = pdg_el(j)
            
            if (dg%pa.eq.0) then

               dg%pa = 1

            endif

            do I = 1,dg%NAGP(dg%pa)

               do k = 1,dg%dofs(j)

                  eta2(j) = eta2(j)+ dg%WAGP(I,dg%pa) * dg%PHI_AREA(K,I,dg%pa) * dg%ZE(K,j,1) * 0.5D0 *AREAS(j)
                  UU2(j)  = UU2(j) + dg%WAGP(I,dg%pa) * dg%PHI_AREA(K,I,dg%pa) * dg%QX(K,j,1) * 0.5D0 *AREAS(j)
                  VV2(j)  = VV2(j) + dg%WAGP(I,dg%pa) * dg%PHI_AREA(K,I,dg%pa) * dg%QY(K,j,1) * 0.5D0 *AREAS(j)


#ifdef TRACE
                  tracer(J)  = tracer(J)  +  dg%WAGP(I,dg%pa) * dg%PHI_AREA(K,I,dg%pa) * dg%iota(K,j,1)  * 0.5D0*AREAS(j)
                  tracer2(J) = tracer2(J) +  dg%WAGP(I,dg%pa) * dg%PHI_AREA(K,I,dg%pa) * dg%iotaa(K,j,1) * 0.5D0*AREAS(j) 
#endif

#ifdef CHEM
                  tracer(J)  = tracer(J)  +  dg%WAGP(I,dg%pa) * dg%PHI_AREA(K,I,dg%pa) * dg%iota(K,j,1) * 0.5D0*AREAS(j)
                  tracer2(J) = tracer2(J) +  dg%WAGP(I,dg%pa) * dg%PHI_AREA(K,I,dg%pa) * dg%iota2(K,j,1) * 0.5D0*AREAS(j)
#endif

#ifdef DYNP
                  dyn_P(J)  = dyn_P(J)  +  dg%WAGP(I,dg%pa) * dg%PHI_AREA(K,I,dg%pa) * dg%dynP(K,j,1) * 0.5D0*AREAS(j)
#endif

#ifdef SED_LAY
                  do l=1,s%layers

                     bed_int(J,l) = bed_int(J,l) + dg%WAGP(I,dg%pa) 
&                    *dg%PHI_AREA(K,I,dg%pa)*dg%bed(K,J,1,l)*0.5D0*AREAS(j)

                     DPe(j)  = DPe(j) + dg%WAGP(I,dg%pa)*dg%PHI_AREA(K,I,dg%pa)*dg%bed(K,J,1,l) 
&                    * 0.5D0 *AREAS(j)
                     
                  enddo

                  tracer2(J) = tracer2(J) +  dg%WAGP(I,dg%pa) * dg%PHI_AREA(K,I,dg%pa) * dg%iotaa2(K,j,1) * 0.5D0*AREAS(j)
#else

                  DPe(j)  = DPe(j) + dg%WAGP(I,dg%pa) * dg%PHI_AREA(K,I,dg%pa) * dg%HB(K,j,1) * 0.5D0 *AREAS(j)
                  tracer2(J) = tracer2(J) +  dg%WAGP(I,dg%pa) * dg%PHI_AREA(K,I,dg%pa) * dg%iotaa2(K,j,1) * 0.5D0*AREAS(j)

#endif
                  
               enddo

            enddo

            DEPTH =  (eta2(j) + dpe(j)) 
            FH_NL = 1.D0/(NLEQ*DEPTH + LEQ)
            eta2(j) = eta2(j)/(AREAS(J))
            UU2(j)  = UU2(j) *FH_NL
            VV2(j)  = VV2(j) *FH_NL
            dpe(j)  = dpe(j)/(AREAS(J))

#ifdef SED_LAY
            do l=1,s%layers
               bed_int(j,l) = bed_int(j,l)/(AREAS(J))
            enddo
#endif

#ifdef TRACE            
            tracer(J)  = tracer(J)  * FH_NL
            tracer2(J) = tracer2(J) * FH_NL  
#endif

#ifdef CHEM
            tracer(J)  = tracer(J)  * FH_NL
            tracer2(J) = tracer2(J) * FH_NL
#endif

#ifdef DYNP            
            dyn_P(J)  = dyn_P(J)  * FH_NL
#endif
            
         enddo

      endif

      iota_int = 0.D0
      
      do j = 1,ne
         
         DO I = 1,dg%NAGP(dg%pa)
            
            iota_error = 0.D0
            
            do k = 1,dg%dofs(j)
               
               iota_error = iota_error + (dg%iota(k,j,1)- dg%iotaa(k,j,1)) * dg%PHI_AREA(k,I,dg%pa) 
               
            ENDDO
            
            iota_int = iota_int + iota_error**2.D0 * dg%WAGP(I,dg%pa)*AREAS(J)
            
         ENDDO
         
      enddo

!...  
!...  OUTPUT ELEVATION RECORDING STATION INFORMATION IF NOUTE<>0 AND THE
!.... TIME STEP FALLS WITHIN THE SPECIFIED WINDOW
!...  CALCULATE ELEVATION SOLUTIONS AT STATIONS USING INTERPOLATION
!...  
      
      IF(NOUTE.NE.0) THEN
         IF((IT.GT.NTCYSE).AND.(IT.LE.NTCYFE).OR.FORCE_WRITE) THEN
            NSCOUE=NSCOUE+1
            IF(NSCOUE.EQ.NSPOOLE.OR.FORCE_WRITE) THEN
               DO I=1,NSTAE           
                                !EE1=ETA2(NM(NNE(I),1))
                                !EE2=ETA2(NM(NNE(I),2))
                                !EE3=ETA2(NM(NNE(I),3))

                  ee1=eta2(nne(i))

                  NC1=NODECODE(NM(NNE(I),1))
                  NC2=NODECODE(NM(NNE(I),2))
                  NC3=NODECODE(NM(NNE(I),3))
                  NCELE=NC1*NC2*NC3
!     IF(NCELE.EQ.1) ET00(I)=EE1*STAIE1(I)+EE2*STAIE2(I)
!     &                                              +EE3*STAIE3(I)
                  if (ncele.eq.1) then
                     ET00(I)=ee1
                  endif
                  IF(NCELE.EQ.0) ET00(I)=-99999.
                  
                  ee1 = dpe(nne(i))

                                !EE1=DP(NM(NNE(I),1))
                                !EE2=DP(NM(NNE(I),2))
                                !EE3=DP(NM(NNE(I),3))
                  IF(NCELE.EQ.1) BT00(I)= ee1

                                !EE1*STAIE1(I)+EE2*STAIE2(I)
                                !&                 +EE3*STAIE3(I)
                  IF(NCELE.EQ.0) BT00(I)=-99999.
               END DO
               IF(ABS(NOUTE).EQ.1) THEN
                  WRITE(61,2120) time_A,IT
                  WRITE(82,2120) time_A,IT
                  DO I=1,NSTAE
                     WRITE(61,2453) I,ET00(I)
                     WRITE(82,2453) I,BT00(I)
                  END DO
                  IESTP = IESTP+1+NSTAE
               ENDIF
               IF(ABS(NOUTE).EQ.2) THEN
                  WRITE(61,REC=IESTP+1) time_A
                  WRITE(61,REC=IESTP+2) IT
                  WRITE(82,REC=IESTP+1) time_A
                  WRITE(82,REC=IESTP+2) IT
                  IESTP = IESTP + 2
                  DO I=1,NSTAE
                     WRITE(61,REC=IESTP+I) ET00(I)
                     WRITE(82,REC=IESTP+I) BT00(I)
                  END DO
                  IESTP = IESTP + NSTAE
               ENDIF
               NSCOUE=0
            ENDIF
         ENDIF
         IF(IT.EQ.NTCYFE) CLOSE(61)
         IF(IT.EQ.NTCYFE) CLOSE(82)
      ENDIF

!...  
!...  OUTPUT VELOCITY RECORDING STATION TIME SERIES INFORMATION IF NOUTV<>0
!.... AND THE TIME STEP FALLS WITHIN THE SPECIFIED WINDOW
!...  CALCULATE VELOCITY SOLUTIONS AT STATIONS USING INTERPOLATION
!...  
      IF(NOUTV.NE.0) THEN
         IF((IT.GT.NTCYSV).AND.(IT.LE.NTCYFV).OR.FORCE_WRITE) THEN
            NSCOUV=NSCOUV+1
            IF(NSCOUV.EQ.NSPOOLV.OR.FORCE_WRITE) THEN
               DO I=1,NSTAV
                  U11=UU2(nnv(I))
                                !U22=UU2(NM(NNV(I),2))
                                !U33=UU2(NM(NNV(I),3))
                  V11=VV2(NNV(I))
                                !V22=VV2(NM(NNV(I),2))
                                !V33=VV2(NM(NNV(I),3))
                  UU00(I)= u11  !U11*STAIV1(I)+U22*STAIV2(I)+U33*STAIV3(I)
                  VV00(I)= v11  !V11*STAIV1(I)+V22*STAIV2(I)+V33*STAIV3(I)
               END DO
               IF(ABS(NOUTV).EQ.1) THEN
                  WRITE(62,2120) time_A,IT
                  DO I=1,NSTAV
                     WRITE(62,2454) I,UU00(I),VV00(I)
                  END DO
                  IVSTP = IVSTP+1+NSTAV
               ENDIF
               IF(ABS(NOUTV).EQ.2) THEN
                  WRITE(62,REC=IVSTP+1) time_A
                  WRITE(62,REC=IVSTP+2) IT
                  IVSTP = IVSTP + 2
                  DO I=1,NSTAV
                     WRITE(62,REC=IVSTP+2*I-1) UU00(I)
                     WRITE(62,REC=IVSTP+2*I) VV00(I)
                  END DO
                  IVSTP = IVSTP + 2*NSTAV
               ENDIF
               NSCOUV=0
            ENDIF
         ENDIF
         IF(IT.EQ.NTCYFV) CLOSE(62)
      ENDIF

!...  
!...  OUTPUT CONCENTRATION RECORDING STATION INFORMATION IF NOUTC<>0 AND THE
!.... TIME STEP FALLS WITHIN THE SPECIFIED WINDOW
!...  CALCULATE CONCENTRATION SOLUTIONS AT STATIONS USING INTERPOLATION
!...  

      IF(NOUTC.NE.0) THEN
         IF((IT.GT.NTCYSC).AND.(IT.LE.NTCYFC).OR.FORCE_WRITE) THEN
            NSCOUC=NSCOUC+1
            IF(NSCOUC.EQ.NSPOOLC.OR.FORCE_WRITE) THEN
               DO I=1,NSTAC
                  NM1=NM(NNC(I),1)
                  NM2=NM(NNC(I),2)
                  NM3=NM(NNC(I),3)
                  HH2N1=dpe(nnc(i))+IFNLFA*ETA2(nnv(i))  !DP(NM1)+IFNLFA*ETA2(NM1)
                  !HH2N2=  DP(NM2)+IFNLFA*ETA2(NM2)
                  !HH2N3=  DP(NM3)+IFNLFA*ETA2(NM3)
                  C1=CH1(NM1)/HH2N1
                  C2=0.D0
                  C3=0.D0
                  !C2=CH1(NM2)/HH2N2
                  !C3=CH1(NM3)/HH2N3
                  NC1=NODECODE(NM1)
                  NC2=NODECODE(NM2)
                  NC3=NODECODE(NM3)
                  NCELE=NC1*NC2*NC3
                  IF(NCELE.EQ.1) CC00(I)=C1*STAIC1(I)+C2*STAIC2(I)&
                       +C3*STAIC3(I)
                  IF(NCELE.EQ.0) CC00(I)=-99999.
               END DO
               IF(ABS(NOUTC).EQ.1) THEN
                  WRITE(81,2120) time_A,IT
                  DO I=1,NSTAC
                     WRITE(81,2453) I,CC00(I)
                  END DO
                  ICSTP = ICSTP+1+NSTAC
               ENDIF
               IF(ABS(NOUTC).EQ.2) THEN
                  WRITE(81,REC=ICSTP+1) time_A
                  WRITE(81,REC=ICSTP+2) IT
                  ICSTP = ICSTP + 2
                  DO I=1,NSTAC
                     WRITE(81,REC=ICSTP+I) CC00(I)
                  END DO
                  ICSTP = ICSTP + NSTAC
               ENDIF
               NSCOUC=0
            ENDIF
         ENDIF
         IF(IT.EQ.NTCYFC) CLOSE(81)
      ENDIF

!...  
!...  OUTPUT METEOROLOGICAL RECORDING STATION INFORMATION IF NWS>0 AND THE
!.... TIME STEP FALLS WITHIN THE SPECIFIED WINDOW
!...  CALCULATE METEOROLOGICAL SOLUTIONS AT STATIONS USING INTERPOLATION
!...  

      IF((NWS.NE.0).AND.(NOUTM.NE.0)) THEN
         IF((IT.GT.NTCYSM).AND.(IT.LE.NTCYFM).OR.FORCE_WRITE) THEN
            NSCOUM=NSCOUM+1
            IF(NSCOUM.EQ.NSPOOLM.OR.FORCE_WRITE) THEN
               DO I=1,NSTAM
                  NM1=NM(NNM(I),1)
                  NM2=NM(NNM(I),2)
                  NM3=NM(NNM(I),3)
                  U11=wvnxout(NM1)
                  U22=wvnxout(NM2)
                  U33=wvnxout(NM3)
                  V11=wvnyout(NM1)
                  V22=wvnyout(NM2)
                  V33=wvnyout(NM3)
                  P11=PR2(NM1)
                  P22=PR2(NM2)
                  P33=PR2(NM3)
                  RMU00(I)=U11*STAIM1(I)+U22*STAIM2(I)+U33*STAIM3(I)
                  RMV00(I)=V11*STAIM1(I)+V22*STAIM2(I)+V33*STAIM3(I)
                  RMP00(I)=P11*STAIM1(I)+P22*STAIM2(I)+P33*STAIM3(I)
               END DO
               IF(ABS(NOUTM).EQ.1) THEN
                  WRITE(71,2120) time_A,IT
                  WRITE(72,2120) time_A,IT
                  DO I=1,NSTAM
                     WRITE(71,2453) I,RMP00(I)
                     WRITE(72,2454) I,RMU00(I),RMV00(I)
                  END DO
                  IPSTP=IPSTP+1+NSTAM
                  IWSTP=IWSTP+1+NSTAM
               ENDIF
               IF(ABS(NOUTM).EQ.2) THEN
                  WRITE(71,REC=IPSTP+1) time_A
                  WRITE(71,REC=IPSTP+2) IT
                  WRITE(72,REC=IWSTP+1) time_A
                  WRITE(72,REC=IWSTP+2) IT
                  IPSTP=IPSTP+2
                  IWSTP=IWSTP+2
                  DO I=1,NSTAM
                     WRITE(71,REC=IPSTP+I) RMP00(I)
                     WRITE(72,REC=IWSTP+2*I-1) RMU00(I)
                     WRITE(72,REC=IWSTP+2*I) RMV00(I)
                  END DO
                  IPSTP=IPSTP+NSTAM
                  IWSTP=IWSTP+2*NSTAM
               ENDIF
               NSCOUM=0
            ENDIF
         ENDIF
         IF(IT.EQ.NTCYFM) THEN
            CLOSE(71)
            CLOSE(72)
         ENDIF
      ENDIF

!.....Output the gloabl elevation data if NOUTGE ~= 0 and the
!.....time step falls within the specified window

      IF (NOUTGE.NE.0) THEN
         IF ((IT.GT.NTCYSGE).AND.(IT.LE.NTCYFGE).OR.FORCE_WRITE) THEN
            NSCOUGE = NSCOUGE + 1
            IF (NSCOUGE.EQ.NSPOOLGE.OR.FORCE_WRITE) THEN
               IF (ABS(NOUTGE).EQ.1) THEN
                  WRITE(63,2120) TIME_A, IT
                  WRITE(88,2120) TIME_A, IT
                  WRITE(89,2120) TIME_A, IT
                  WRITE(4441,2120) TIME_A, IT
                                !WRITE(895,21 20) TIME_A, IT
                  WRITE(631,2120) TIME_A, IT
 2120             FORMAT(2X,E20.10,5X,I10)
                  DO I=1,NE
                     IF (ABS(ETA2(I)).LE.(10.0**(-30))) ETA2(I) = 0.D0
                     WRITE(63,2453) I,ETA2(I)
                                !IF (NODECODE(I).EQ.0) WRITE(63,2453) I,-99999.
                     IF (ABS(tracer(I)).LE.(10.0**(-30))) tracer(I) = 0.D0
                     WRITE(88,2453) I,tracer(I)
                                !IF (NODECODE(I).EQ.0) WRITE(88,2453) I,.
                     IF (ABS(tracer2(I)).LE.(10.0**(-30))) tracer2(I) = 0.D0
                     WRITE(89,2453) I,tracer2(I)
                                !IF (NODECODE(I).EQ.0) WRITE(89,2453) I,0.\
                     IF (ABS(dpe(I)).LE.(10.0**(-30))) dpe(I) = 0.D0
                     write(4441,2453) I,dpe(I)
 2453                FORMAT(2X,I8,2X,E16.8E3)
                  ENDDO
 
!     cnd...for tecplot results
!     cmm modified 4/6/09
!     cem...modified for cell-centered data, nodes give incorrect results for p>1
                  
                  if (mod(it,NSPOOLGE).eq.0) then
!Casey 120813: Begin the OUT_TEC conditional.

#ifdef OUT_TEC
                     
                     if (ModetoNode.eq.0) then

           write(777,*) 'ZONE ZONETYPE=FETRIANGLE ',
     $     'NODES=', np, 
     $     ' ELEMENTS=', ne, 
     $     ' DATAPACKING=BLOCK ',
     $     'SOLUTIONTIME=',time_a,
     $     'VARLOCATION=([3,4,5,6,7,8,10,11,12,13,14]=CELLCENTERED)'
                     do i=1,np
                        write(777,7777)  x(i)
                     enddo
                     do i=1,np
                        write(777,7777)  y(i)
                     enddo
                     do i=1,ne 
                        write(777,7777)  dpe(i)
                     enddo
                     do i=1,ne 
                        write(777,7777)  eta2(i)
                     enddo
                     do i=1,ne 
                        write(777,7777)  eta2(i)+dpe(i)
                     enddo
                     do i=1,ne 
                        write(777,7777)  uu2(i)
                     enddo
                     do i=1,ne 
                        write(777,7777)  vv2(i)
                     enddo
                     do i=1,ne 
                        write(777,7777)  sqrt(uu2(i)**2+vv2(i)**2)
                     enddo
                     do i=1,np 
                        write(777,7777)  sqrt(wsx2(i)**2+wsy2(i)**2)
                     enddo
                     do i=1,ne 
                        write(777,7777)  entrop(4,i) !DBLE(pdg_el(i))
                     enddo
                     do i=1,ne 
                        write(777,7777)  tracer(i)
                     enddo
                     do i=1,ne 
                        write(777,7777)  tracer2(i)
                     enddo
                     do i=1,ne 
                        write(777,7777)  abs(tracer(i)+tracer2(i))
                     enddo
                     do i=1,ne 
                        write(777,7777)  abs(tracer(i)-tracer2(i))    
                     enddo

                  else

                 write(777,*) 'ZONE ZONETYPE=FETRIANGLE ',
     $                'NODES=', np, 
     $                ' ELEMENTS=', ne, 
     $                ' DATAPACKING=POINT ','SOLUTIONTIME=',time_a
                 do i=1,np
                    if (ics.eq.2) then
                       write(777,7777) slam(i)/deg2rad, sfea(i)/deg2rad, 
!                 write(777,7777) x(i), y(i), 
     $                   dp(i), eta2(i), eta2(i)+dp(i),uu2(i),vv2(i),
     $           sqrt(uu2(i)**2+vv2(i)**2),sqrt(wsx2(i)**2+wsy2(i)**2),
     $                   myproc
                    else
                       write(777,7777) x(i), y(i), 
     $                   dp(i), eta2(i), eta2(i)+dp(i),uu2(i),vv2(i),
     $           sqrt(uu2(i)**2+vv2(i)**2),sqrt(wsx2(i)**2+wsy2(i)**2),
     $                   myproc
                    endif

                 enddo
              endif


 7777         format(9f20.8,i10)
              do i=1,ne
                 write(777,"(3i12)") nm(i,1), nm(i,2), nm(i,3)
              enddo


#endif

!.....Write DG.63 results
                     
                     DO J = 1,NE
                        DO K = 1,dg%DOFS(J)
                           WRITE(631,*) J, dg%ZE(K,J,1)
                        ENDDO
                     ENDDO
                     
!.....Write DG.65 results, elemental statuses (wet/dry)
                     
                     IF (NOLIFA.GE.2) THEN
                        WRITE(651,2120) TIME_A, IT
                        DO J = 1,NE
                           WRITE(651,2455) J,dg%WDFLG(J) 
 2455                      FORMAT(2X,I8,2X,I2)
                        ENDDO
                     ENDIF
                     
                     IGEP = IGEP + 1 + NP
                  ENDIF
                  IF (ABS(NOUTGE).EQ.2) THEN
                     WRITE(63,REC=IGEP+1) TIME_A
                     WRITE(63,REC=IGEP+2) IT
                     IGEP = IGEP + 2
                     DO I = 1,NP
                        WRITE(63,REC=IGEP+I)ETA2(I)
                                !IF (NODECODE(I).EQ.0) WRITE(63,REC=IGEP+I) -99999.
                     ENDDO
                     IF(SEDFLAG.GE.2) THEN
                        WRITE(84,REC=IGEP+1) TIME_A
                        WRITE(84,REC=IGEP+2) IT
                        DO I = 1,NP
                           WRITE(84,REC=IGEP+I) DP(I)
                        ENDDO
                     ENDIF
                     IGEP = IGEP + NP
                  ENDIF
                  NSCOUGE = 0
               ENDIF
            ENDIF
            IF (IT.EQ.NTCYFGE) CLOSE(63)
            IF (IT.EQ.NTCYFGE) CLOSE(631)
            IF (IT.EQ.NTCYFGE.AND.NOLIFA.GE.2) CLOSE(651)
            IF (IT.EQ.NTCYFGE.AND.SEDFLAG.GE.1) CLOSE(84)

         endif
      ENDIF

!...  
!...  OUTPUT GLOBAL VELOCITY DATA IF NOUTGV<>0 AND THE
!.... TIME STEP FALLS WITHIN THE SPECIFIED WINDOW
!...  
      IF (NOUTGV.NE.0) THEN
         IF ((IT.GT.NTCYSGV).AND.(IT.LE.NTCYFGV).OR.FORCE_WRITE) THEN
            NSCOUGV=NSCOUGV+1
            IF (NSCOUGV.EQ.NSPOOLGV.OR.FORCE_WRITE) THEN
               IF(ABS(NOUTGV).EQ.1) THEN
                  WRITE(64,2120) TIME_A,IT
                  WRITE(641,2120) TIME_A,IT
                  DO I = 1,NE
!.....Had trouble writing out numbers with less than E-99 (EJK)
                     IF (ABS(UU2(I)).LE.(10.0**(-30))) UU2(I) = 0.D0
                     IF (ABS(VV2(I)).LE.(10.0**(-30))) VV2(I) = 0.D0
                     WRITE(64,2454) I,UU2(I),VV2(I)
 2454                FORMAT(2X,I8,2(2X,E16.8E3))
                  ENDDO
                  
!.....Write DG.64 results

                  DO J = 1,NE
                     DO K = 1,dg%DOFS(J)
                        WRITE(641,*) J, dg%QX(K,J,1), dg%QY(K,J,1)
                     ENDDO
                  ENDDO
                  
                  IGVP = IGVP + 1 + NP
               ENDIF
               IF (ABS(NOUTGV).EQ.2) THEN
                  WRITE(64,REC=IGVP+1) time_A
                  WRITE(64,REC=IGVP+2) IT
                  IGVP = IGVP + 2
                  DO I=1,NP
                     WRITE(64,REC=IGVP+2*I-1) UU2(I)
                     WRITE(64,REC=IGVP+2*I) VV2(I)
                  ENDDO
                  IGVP = IGVP + 2*NP
               ENDIF
               NSCOUGV=0
            ENDIF
         ENDIF
         IF(IT.EQ.NTCYFGV) CLOSE(64)
         IF(IT.EQ.NTCYFGV) CLOSE(641)
      ENDIF

!...  
!...  OUTPUT GLOBAL WIND STRESS and atmospheric pressure data IF NOUTGW<>0 AND THE
!.... TIME STEP FALLS WITHIN THE SPECIFIED WINDOW
!...  
      IF((NWS.NE.0).AND.(NOUTGW.NE.0)) THEN
         IF((IT.GT.NTCYSGW).AND.(IT.LE.NTCYFGW).OR.FORCE_WRITE) THEN
            NSCOUGW=NSCOUGW+1
            IF(NSCOUGW.EQ.NSPOOLGW.OR.FORCE_WRITE) THEN
               IF(ABS(NOUTGW).EQ.1) THEN
                  write(73,2120) time_A,it
                  WRITE(74,2120) time_A,IT
                  DO I=1,NP
                     write(73,2453) i,pr2(i)
                     WRITE(74,2454) i,wvnxout(i),wvnyout(i)
                  ENDDO
                  igpp = igpp+1+np
                  IGWP = IGWP+1+NP
               ENDIF
               IF(ABS(NOUTGW).EQ.2) THEN
                  WRITE(73,REC=igpp+1) time_A
                  WRITE(73,REC=igpp+2) IT
                  WRITE(74,REC=IGWP+1) time_A
                  WRITE(74,REC=IGWP+2) IT
                  igpp = igpp + 2
                  IGWP = IGWP + 2
                  DO I=1,NP
                     write(73,rec=igpp+i) pr2(i)
                     WRITE(74,REC=IGWP+2*I-1) wvnxout(i)
                     WRITE(74,REC=IGWP+2*I) wvnyout(i)
                  END DO
                  igpp = igpp + np
                  IGWP = IGWP + 2*NP
               ENDIF
               NSCOUGW=0
            ENDIF
         ENDIF
         IF(IT.EQ.NTCYFGW) then
            close(73)
            close(74)
         ENDIF
      endif

      if (it.eq.nt) then
         open(963,FILE=S%DIRNAME//'/'//'maxele.63')
         write(963,*) np
         do i=1,np
            write(963,9633) etamax(i)
         enddo
 9633    format(f20.8)
         close(963)
#ifdef OUT_TEC
         if (ModetoNode.eq.1) then
            write(778,*) 'TITLE = "dgswem output"'
            write(778,*) 
     $           'VARIABLES = "x", "y", "maxeta"'
            write(778,*) 'ZONE ZONETYPE=FETRIANGLE ',
     $           'NODES=', np, 
     $           ' ELEMENTS=', ne, 
     $           ' DATAPACKING=POINT ','SOLUTIONTIME=',time_a
            do i=1,np
               write(778,7778) slam(i)/deg2rad, sfea(i)/deg2rad, 
     $              etamax(i)
            enddo
 7778       format(3f20.8)
            do i=1,ne
               write(778,"(3i12)") nm(i,1), nm(i,2), nm(i,3)
            enddo
            close(778)
         endif
#endif
      endif

!...  
!...  OUTPUT GLOBAL CONCENTRATION DATA IF NOUTGC<>0 AND THE
!.... TIME STEP FALLS WITHIN THE SPECIFIED WINDOW
!...  
      IF(NOUTGC.NE.0) THEN
         IF((IT.GT.NTCYSGC).AND.(IT.LE.NTCYFGC).OR.FORCE_WRITE) THEN
            NSCOUGC=NSCOUGC+1
            IF(NSCOUGC.EQ.NSPOOLGC.OR.FORCE_WRITE) THEN
               IF(ABS(NOUTGC).EQ.1) THEN
                  WRITE(83,2120) time_A,IT
                  DO I=1,NP
                     HH2=DP(I)+IFNLFA*ETA2(I)
                     C1=CH1(I)/HH2
                     IF(NODECODE(I).EQ.1) WRITE(83,2453) I,C1
                     IF(NODECODE(I).EQ.0) WRITE(83,2453) I,-99999.
                  ENDDO
                  IGCP=IGCP+1+NP
               ENDIF
               IF(ABS(NOUTGC).EQ.2) THEN
                  WRITE(83,REC=IGEP+1) time_A
                  WRITE(83,REC=IGEP+2) IT
                  IGCP = IGCP + 2
                  DO I=1,NP
                     HH2=DP(I)+IFNLFA*ETA2(I)
                     C1=CH1(I)/HH2
                     IF(NODECODE(I).EQ.1) WRITE(83,REC=IGCP+I) C1
                     IF(NODECODE(I).EQ.0) WRITE(83,REC=IGCP+I) -99999.
                  ENDDO
                  IGCP=IGCP+NP
               ENDIF
               NSCOUGC=0
            ENDIF
         ENDIF
         IF(IT.EQ.NTCYFGC) CLOSE(83)
      ENDIF

!...  
!...  IF IHARIND=1 AND THE TIME STEP FALLS WITHIN THE SPECIFIED WINDOW AND
!...  ON THE SPECIFIED INCREMENT, USE MODEL RESULTS TO UPDATE HARMONIC
!...  ANALYSIS MATRIX AND LOAD VECTORS.  NOTE: AN 8 BYTE RECORD SHOULD BE
!...  USED THROUGHOUT THE HARMONIC ANALYSIS SUBROUTINES, EVEN ON 32 BIT
!...  WORKSTATIONS, SINCE IN THAT CASE THE HARMONIC ANALYSIS IS DONE IN
!...  DOUBLE PRECISION.
!...  
#ifdef HARM
      IF(IHARIND.EQ.1) THEN
         IF((IT.GT.ITHAS).AND.(IT.LE.ITHAF)) THEN
            IF(ICHA.EQ.NHAINC) ICHA=0
            ICHA=ICHA+1
            IF(ICHA.EQ.NHAINC) THEN
!...  
!.....UPDATE THE LHS MATRIX
!...  
               CALL LSQUPDLHS(timeh,IT)
!...  
!.....IF DESIRED COMPUTE ELEVATION STATION INFORMATION AND UPDATE LOAD VECTOR
!...  
               IF(NHASE.EQ.1) THEN
                  DO I=1,NSTAE
!     EE1=ETA2(NM(NNE(I),1))
!     EE2=ETA2(NM(NNE(I),2))
!     EE3=ETA2(NM(NNE(I),3))
!     ET00(I)=EE1*STAIE1(I)+EE2*STAIE2(I)+EE3*STAIE3(I)
                     ET00(I) = dg%ZE(1,NNE(I),1)
                     DO K = 2,dg%DOFH
                        ET00(I) = ET00(I) + dg%ZE(K,NNE(I),1)*dg%PHI_STAE(K,I)
                     ENDDO
                  ENDDO
                  CALL LSQUPDES(ET00,NSTAE)
               ENDIF
!...  
!.....IF DESIRED COMPUTE VELOCITY STATION INFORMATION AND UPDATE LOAD VECTOR
!...  
               IF(NHASV.EQ.1) THEN
                  DO I=1,NSTAV
!     U11=UU2(NM(NNV(I),1))
!     U22=UU2(NM(NNV(I),2))
!     U33=UU2(NM(NNV(I),3))
!     V11=VV2(NM(NNV(I),1))
!     V22=VV2(NM(NNV(I),2))
!     V33=VV2(NM(NNV(I),3))
!     UU00(I)=U11*STAIV1(I)+U22*STAIV2(I)+U33*STAIV3(I)
!     VV00(I)=V11*STAIV1(I)+V22*STAIV2(I)+V33*STAIV3(I)
                     DP1     = DPE(nnv(i)) !DP(NM(NNV(I),1))
                     !DP2     = DPE(nnv(i)) !DP(NM(NNV(I),2))
                     !DP3     = DPE(nnv(i)) !DP(NM(NNV(I),3))
                     !DP00    = DP1*STAIV1(I) + DP2*STAIV2(I) +DP3*STAIV3(I)
                     ZE00    = dg%ZE(1,NNV(I),1)
                     UU00(I) = dg%QX(1,NNV(I),1)
                     VV00(I) = dg%QY(1,NNV(I),1)
                     DO K = 2,dg%DOFH
                        ZE00    = ZE00    + dg%ZE(K,NNV(I),1)*dg%PHI_STAV(K,I)
                        UU00(I) = UU00(I) + dg%QX(K,NNV(I),1)*dg%PHI_STAV(K,I)
                        VV00(I) = VV00(I) + dg%QY(K,NNV(I),1)*dg%PHI_STAV(K,I)
                     ENDDO
                     DEPTH   = ZE00 + DP1 !DP00
                     FH_NL   = 1.D0/(NLEQ*DEPTH + LEQ)
                     UU00(I) = UU00(I)*FH_NL
                     VV00(I) = VV00(I)*FH_NL
                  ENDDO
                  CALL LSQUPDVS(UU00,VV00,NSTAV)
               ENDIF
!...  
!.....IF DESIRED UPDATE GLOBAL ELEVATION LOAD VECTOR
!...  
               IF(NHAGE.EQ.1) CALL LSQUPDEG(ETA2,NP)
!...  
!.....IF DESIRED UPDATE GLOBAL VELOCITY LOAD VECTOR
!...  
               IF(NHAGV.EQ.1) CALL LSQUPDVG(UU2,VV2,NP)

            ENDIF
         ENDIF

!...  LINES TO COMPUTE MEANS AND VARIANCES

         if (CHARMV) then
            IF(IT.GT.ITMV) THEN
               NTSTEPS=NTSTEPS+1
               DO I=1,NP
                  ELAV(I)=ELAV(I)+ETA2(I)
                  XVELAV(I)=XVELAV(I)+UU2(I)
                  YVELAV(I)=YVELAV(I)+VV2(I)
                  ELVA(I)=ELVA(I)+ETA2(I)*ETA2(I)
                  XVELVA(I)=XVELVA(I)+UU2(I)*UU2(I)
                  YVELVA(I)=YVELVA(I)+VV2(I)*VV2(I)
               END DO
            ENDIF
         endif                  !   charmv


      ENDIF
#endif

!...  
!...  WRITE OUT HOT START INFORMATION IF NHSTAR=1 AND AT CORRECT TIME STEP
!...  NOTE: THE HOT START FILES USE A RECORD LENGTH OF 8 ON BOTH 32 BIT
!.... WORKSTATIONS AND THE 64 BIT CRAY.  THIS IS BECAUSE THE HARMONIC
!.... ANALYSIS IS DONE IN DOUBLE PRECISION (64 BITS) ON WORKSTATIONS.
!...  
#ifdef HOTSTART
      IF(NHSTAR.EQ.1) THEN
         ITEST=(IT/NHSINC)*NHSINC  
         IF(ITEST.EQ.IT) THEN
            IF(IHSFIL.EQ.67) OPEN(67,FILE=S%DIRNAME//'/'//'fort.67',&
            ACCESS='DIRECT',RECL=8)
            IF(IHSFIL.EQ.68) OPEN(68,FILE=S%DIRNAME//'/'//'fort.68',&
           ACCESS='DIRECT',RECL=8)
            IHOTSTP=1
            WRITE(IHSFIL,REC=IHOTSTP) IM
            IHOTSTP=2
            WRITE(IHSFIL,REC=IHOTSTP) TIME_A
            IHOTSTP=3
            WRITE(IHSFIL,REC=IHOTSTP) IT
            DO I=1,NP
               WRITE(IHSFIL,REC=IHOTSTP+1) ETA1(I)
               WRITE(IHSFIL,REC=IHOTSTP+2) ETA2(I)
               WRITE(IHSFIL,REC=IHOTSTP+3) UU2(I)
               WRITE(IHSFIL,REC=IHOTSTP+4) VV2(I)
               IHOTSTP = IHOTSTP + 4
               IF(IM.EQ.10) THEN
                  WRITE(IHSFIL,REC=IHOTSTP+1) CH1(I)
                  IHOTSTP=IHOTSTP+1
               ENDIF
               WRITE(IHSFIL,REC=IHOTSTP+1) NODECODE(I)
               IHOTSTP=IHOTSTP+1
            END DO
            WRITE(IHSFIL,REC=IHOTSTP+1) IESTP
            WRITE(IHSFIL,REC=IHOTSTP+2) NSCOUE
            IHOTSTP=IHOTSTP+2
            WRITE(IHSFIL,REC=IHOTSTP+1) IVSTP
            WRITE(IHSFIL,REC=IHOTSTP+2) NSCOUV
            IHOTSTP=IHOTSTP+2
            WRITE(IHSFIL,REC=IHOTSTP+1) ICSTP
            WRITE(IHSFIL,REC=IHOTSTP+2) NSCOUC
            IHOTSTP=IHOTSTP+2
            WRITE(IHSFIL,REC=IHOTSTP+1) IPSTP
            WRITE(IHSFIL,REC=IHOTSTP+2) IWSTP
            WRITE(IHSFIL,REC=IHOTSTP+2) NSCOUM
            IHOTSTP=IHOTSTP+3
            WRITE(IHSFIL,REC=IHOTSTP+1) IGEP
            WRITE(IHSFIL,REC=IHOTSTP+2) NSCOUGE
            IHOTSTP=IHOTSTP+2
            WRITE(IHSFIL,REC=IHOTSTP+1) IGVP
            WRITE(IHSFIL,REC=IHOTSTP+2) NSCOUGV
            IHOTSTP=IHOTSTP+2
            WRITE(IHSFIL,REC=IHOTSTP+1) IGCP
            WRITE(IHSFIL,REC=IHOTSTP+2) NSCOUGC
            IHOTSTP=IHOTSTP+2
            WRITE(IHSFIL,REC=IHOTSTP+1) IGPP
            WRITE(IHSFIL,REC=IHOTSTP+2) IGWP
            WRITE(IHSFIL,REC=IHOTSTP+3) NSCOUGW
            IHOTSTP=IHOTSTP+3
!...  
!...  IF APPROPRIATE ADD HARMONIC ANALYSIS INFORMATION TO HOT START FILE
!...  
            IF((IHARIND.EQ.1).AND.(IT.GT.ITHAS)) THEN
               WRITE(IHSFIL,REC=IHOTSTP+1) ICHA
               IHOTSTP = IHOTSTP + 1
               CALL HAHOUT(NP,NSTAE,NSTAV,NHASE,NHASV,NHAGE,NHAGV,&
              IHSFIL,IHOTSTP)
!     
               IF(NHASE.EQ.1) CALL HAHOUTES(NSTAE,IHSFIL,IHOTSTP)
               IF(NHASV.EQ.1) CALL HAHOUTVS(NSTAV,IHSFIL,IHOTSTP)
               IF(NHAGE.EQ.1) CALL HAHOUTEG(NP,IHSFIL,IHOTSTP)
               IF(NHAGV.EQ.1) CALL HAHOUTVG(NP,IHSFIL,IHOTSTP)
            ENDIF

            if( CHARMV) then
               IF((IHARIND.EQ.1).AND.(IT.GT.ITMV)) THEN
                  IHOTSTP=IHOTSTP+1
                  WRITE(IHSFIL,REC=IHOTSTP) NTSTEPS
                  IF(NHAGE.EQ.1) THEN
                     DO I=1,NP
                        WRITE(IHSFIL,REC=IHOTSTP+1) ELAV(I)
                        WRITE(IHSFIL,REC=IHOTSTP+2) ELVA(I)
                        IHOTSTP=IHOTSTP+2
                     END DO
                  ENDIF
                  IF(NHAGV.EQ.1) THEN
                     DO I=1,NP
                        WRITE(IHSFIL,REC=IHOTSTP+1) XVELAV(I)
                        WRITE(IHSFIL,REC=IHOTSTP+2) YVELAV(I)
                        WRITE(IHSFIL,REC=IHOTSTP+3) XVELVA(I)
                        WRITE(IHSFIL,REC=IHOTSTP+4) YVELVA(I)
                        IHOTSTP=IHOTSTP+4
                     END DO
                  ENDIF
               ENDIF
            endif               !  charmv
!     
            
            IF (C3D) THEN
!     CALL HSTART3D_OUT()
            ENDIF

!...  
!...  CLOSE THE HOT START OUTPUT FILE
!...  
            CLOSE(IHSFIL)
            IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) THEN
               WRITE(6,24541) IHSFIL,IT,TIME_A
            ENDIF
            WRITE(16,24541) IHSFIL,IT,TIME_A
24541       FORMAT(1X,'HOT START OUTPUT WRITTEN TO UNIT ',I2,&
           ' AT TIME STEP = ',I9,' TIME = ',E15.8)
            IF(IHSFIL.EQ.67) THEN
               IHSFIL=68
            ELSE
               IHSFIL=67
            ENDIF
         ENDIF
      ENDIF
#endif      


!...  FIND AND PRINT TO UNIT 6, THE MAXIMUM ELEVATION, THE MAXIMUM VELOCITY
!.... AND THE NODE NUMBERS AT WHICH THEY OCCUR ON MYPROC=0
!.... IF ELMAX EXCEEDS THRESHOLD, PRINT INFORMATION ON ALL PROCESSORS WHERE
!.... THIS OCCURS 
!...  
      IF(NSCREEN.EQ.1) THEN
         ELMAX=0.0
         VELMAX=0.0
         KEMAX = 0
         KVMAX = 0
!     
         imaxze=0.0
         elmaxe=0.0
         qmaxe=0.0
         imaxq = 0.0
         DO I=1,NP
            IF((NODECODE(I).EQ.1).AND.(ABS(ETA2(I)).GT.ELMAX))THEN
               ELMAX=ABS(ETA2(I))
               KEMAX=I
            ENDIF
            VELABS=UU2(I)*UU2(I)+VV2(I)*VV2(I)
            IF (VELABS.GT.VELMAX) THEN
               VELMAX=VELABS
               KVMAX=I
            ENDIF
         END DO
!     nd
         do i=1,ne
            if (dg%wdflg(i).ne.0) then
               if (abs(dg%ze(1,i,1)).gt.elmaxe) then
                  imaxze=i
                  elmaxe=abs(dg%ze(1,i,1))
               endif
               if (sqrt(dg%qx(1,i,1)**2+dg%qy(1,i,1)**2).gt.qmaxe) then
                  qmaxe=sqrt(dg%qx(1,i,1)**2+dg%qy(1,i,1)**2)
                  imaxq=i
               endif
               
            endif
         enddo
         
         VELMAX=VELMAX**0.5d0
         
#ifdef CMPI
!     sb
!     IF(MYPROC.EQ.0.AND.ELMAX.LT.200.0.AND.KEMAX.GT.0) THEN
         IF( (MYPROC.EQ.0).AND.(ELMAX.LT.200.0).AND.&
        (MOD(IT,NSCREEN_INC).EQ.0) ) THEN
            if (kemax.eq.0) then 
               kemax = 1 
            endif
       print*,'________________________________',&
       '__________________________________________'
            print*,''

         WRITE(6,1991) IT,NUMITR,ETA2(KEMAX),KEMAX,VELMAX,KVMAX,&
           MYPROC
 1991   FORMAT(1X,'TIME STEP =',I8,5X,'ITERATIONS =',I5,&
           /,2X,'ELMAX = ', E11.4,' AT NODE',I7,&
           2X,'SPEEDMAX = ',E11.4,' AT NODE',I7,&
           2X,'ON MYPRnOC = ',I4)
   !write(6,*) 'elmaxe, imaxze, qmaxe, imaxq ',elmaxe,imaxze,
                                !$           qmaxe,imaxq


!$$$  if (chem_flag.eq.1) then
!$$$  print*,'__________________________________________________________________________'
!$$$  print*,'|                       |                                                |'
!$$$  print*,'| Chemistry turned "ON" |',' Maximum mass action =',-minval(MassMax(:)),' |'
!$$$  print*,'|_______________________|________________________________________________|'
!$$$  endif


            if (dg%padapt.eq.1) then

               do k=dg%pl,dg%ph
                  Minp(k) = 0
                  do i = 1,ne
                     if (pdg_el(i).eq.k) then
                        Minp(k) = Minp(k) + 1
                     endif
                  enddo
               enddo

             print*,'_________________________________________________'
             print*,'|                      |                        |'
             print*,'|  Polynomial order    |     # of Elements      |'
             print*,'|______________________|________________________|'
              do k=dg%pl,dg%ph
             print*,'|',k,'         |',Minp(k),'           |'
             print*,'|______________________|________________________|'
              enddo
               print*,''
               print*,'With',ne,'total elements'               
            endif

!$$$  if (dg%slopeflag.eq.4.or.dg%slopeflag.ge.6) then
!$$$  
!$$$  print*,''
!$$$  print*,'p-adaptive slope limiting "ACTIVE."'
!$$$  print*,''
!$$$  print*,dg%lim_count_roll,'elements of a total of'
!$$$  print*,ne,'elements were VERTEX LIMITED'
!$$$  print*,''
!$$$  
!$$$  endif

            print*,'___________________________________',&
          '_______________________________________'

         ENDIF
!     sb
!     IF(ELMAX.GT.200.0.AND.KEMAX.GT.0) THEN
         IF(ELMAX.GT.200.0) THEN
            if (kemax.eq.0) then 
               kemax = 1 
            endif
            WRITE(6,1993) IT,NUMITR,ETA2(KEMAX),KEMAX,VELMAX,KVMAX,&
           MYPROC
            WRITE(16,1993) IT,NUMITR,ETA2(KEMAX),KEMAX,VELMAX,KVMAX,&
           MYPROC
 1993       FORMAT(1X,'TIME STEP =',I8,6X,'ITERATIONS =',I5,&
      /,2X,'ELMAX = ', E11.4,' AT NODE',I7,&
      2X,'SPEEDMAX = ',E11.4,' AT NODE',I7,&
      2X,'ON MYPROC = ',I4,' !!! WARNING - HIGH ELEVATION !!!')
            STOP
         ENDIF


#else
!     sb
!     IF(ELMAX.LT.200.0.AND.KEMAX.GT.0) THEN

         IF(ELMAX.LT.200.0.and.(MOD(IT,NSCREEN_INC).EQ.0) ) THEN
            if (kemax.eq.0) then 
               kemax = 1 
            endif
            print*,'______________________________________',&
            '____________________________________'
            print*,''
            WRITE(6,1992) IT,NUMITR,ETA2(KEMAX),KEMAX,VELMAX,KVMAX
 1992       FORMAT(1X,'TIME STEP =',I8,5X,'ITERATIONS =',I5,&
           /,2X,'ELMAX = ', E11.4,' AT NODE',I7,&
           2X,'SPEEDMAX = ',E11.4,' AT NODE',I7)
                                !write(6,*) 'elmaxe, imaxze, qmaxe, imaxq ',elmaxe,imaxze,
                                !$           qmaxe,imaxq

!$$$  if (chem_flag.eq.1) then
!$$$  print*,'__________________________________________________________________________'
!$$$  print*,'|                       |                                                |'
!$$$  print*,'| Chemistry turned "ON" |',' Maximum mass action =',-minval(MassMax(:)),' |'
!$$$  print*,'|_______________________|________________________________________________|'
!$$$  endif


            if (dg%padapt.eq.1) then
               
               do k=dg%pl,dg%ph
                  Minp(k) = 0
                  do i = 1,ne
                     if (pdg_el(i).eq.k) then
                        Minp(k) = Minp(k) + 1
                     endif
                  enddo
               enddo

             print*,'_________________________________________________'
             print*,'|                      |                        |'
             print*,'|  Polynomial order    |     # of Elements      |'
             print*,'|______________________|________________________|'
               do k=dg%pl,dg%ph
             print*,'|',k,'         |',Minp(k),'           |'
             print*,'|______________________|________________________|'
               enddo
               print*,''
               print*,'With',ne,'total elements'               
            endif

!$$$  if (dg%slopeflag.eq.4.or.dg%slopeflag.ge.6) then
!$$$  
!$$$  print*,''
!$$$  print*,'p-adaptive slope limiting "ACTIVE."'
!$$$  print*,''
!$$$  print*,dg%lim_count_roll,'elements of a total of'
!$$$  print*,ne,'elements were VERTEX LIMITED'
!$$$  print*,''
!$$$  
!$$$  endif

            print*,'___________________________________',&
          '_______________________________________'

         ENDIF
!     sb
!     IF(ELMAX.GT.200.0.AND.KEMAX.GT.0) THEN
         IF(ELMAX.GT.200.0) THEN
          WRITE(6,1994) IT,NUMITR,ETA2(KEMAX),KEMAX,VELMAX,KVMAX
          WRITE(16,1994) IT,NUMITR,ETA2(KEMAX),KEMAX,VELMAX,KVMAX
 1994       FORMAT(1X,'TIME STEP =',I8,6X,'ITERATIONS =',I5,&
           /,2X,'ELMAX = ', E11.4,' AT NODE',I7,&
           2X,'SPEEDMAX = ',E11.4,' AT NODE',I7,&
           2X,' !!! WARNING - HIGH ELEVATION !!!')
            STOP
         ENDIF
#endif
      ENDIF
      
!.....If applicable write out a DG hot start

      IF ((dg%DGHOT.EQ.1).AND.(MOD(IT,dg%DGHOTSPOOL).EQ.0)) THEN
         
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
#ifdef SED_LAY
         OPEN(290,FILE=S%DIRNAME//'/'//'Hot_start.290')
#endif

         WRITE(263,*) dg%PH
         WRITE(264,*) dg%PH, dg%PH
         WRITE(214,*) IT

         DO J = 1,s%MNE
            DO K = 1,dg%DOFH
               WRITE(263,*) dg%ZE(K,J,1)
               WRITE(264,*) dg%QX(K,J,1), dg%QY(K,J,1)
               WRITE(214,*) dg%HB(K,J,1), dg%WDFLG(J)
#ifdef TRACE
               WRITE(288,*) dg%iota(K,J,1)
#endif
#ifdef CHEM
               WRITE(289,*) dg%iota(K,J,1),dg%iota2(K,J,1)
#endif
#ifdef DYNP
               WRITE(291,*) dg%dynP(K,J,1)
#endif
#ifdef SED_LAY
               do ll=1,s%layers
                  WRITE(290,*) dg%bed(K,J,1,ll)
               enddo
#endif
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
#ifdef DYNP
         CLOSE(291)
#endif
#ifdef SED_LAY
         CLOSE(290)
#endif
      ENDIF

!.....If end of run write out DG data

!     IF (IT.EQ.NT) THEN
!     dg%H_TRI = SQRT((X(1)-X(2))**2 + (Y(1)-Y(2))**2)
!     WRITE(88,*) dg%H_TRI,P,0
!     WRITE(89,*) dg%H_TRI,P,0
!     IF (MYPROC == 0) THEN
!     PRINT*,'T FINAL =',NT*DTDP
!     PRINT*,'DT =',DTDP
!     PRINT*,'NT =',NT
!     ENDIF

!     DO J=1,s%MNE
!     ZE_DG(1) = 0.D0
!     QX_DG(1) = 0.D0
!     QY_DG(1) = 0.D0
!     
!     DO K=1,dg%DOF
!     ZE_DG(1) = ZE_DG(1) + dg%PHI_CENTER(K)*dg%ZE(K,J,1)
!     QX_DG(1) = QX_DG(1) + dg%PHI_CENTER(K)*dg%QX(K,J,1)
!     QY_DG(1) = QY_DG(1) + dg%PHI_CENTER(K)*dg%QY(K,J,1)
!     
!     WRITE(88,*) dg%ZE(K,J,1),dg%QX(K,J,1),dg%QY(K,J,1)
!     
!     ENDDO
!     
!     WRITE(89,*) ZE_DG(1),QX_DG(1),QY_DG(1)
!     
!     ENDDO
!     ENDIF
!...  
!...  ****************** TIME STEPPING LOOP ENDS HERE ********************
!...  
      RETURN
      END

!***********************************************************************
!
!     SUBROUTINE:  WRITE_DG_IC
!
!     Outputs DG initial conditions
!
!     Written by Shintaro Bunya, Nov. 2005 
!
!***********************************************************************

      SUBROUTINE WRITE_DG_IC(dg)

      USE GLOBAL
      USE DG

      type (dg_type) :: dg

!.....DG.63.IC = Record of the initial surface elevation
!      OPEN(632,FILE=S%DIRNAME//'/'//'DG.63.IC')
!      WRITE(632,3220) RUNDES, RUNID, AGRID
!      WRITE(632,3645) 1, dg%DOF, DTDP*NSPOOLGE, NSPOOLGE, 1

!      WRITE(632,2120) 0.0, 1
!      DO J = 1,NE
!        DO K = 1,dg%DOF
!          WRITE(632,2453) J, dg%ZE(K,J,1)
!        ENDDO
!      ENDDO
      
!      CLOSE(632)

!.....DG.64.IC = Record of the initial discharge
!      OPEN(642,FILE=S%DIRNAME//'/'//'DG.64.IC')
!      WRITE(642,3220) RUNDES, RUNID, AGRID
!      WRITE(642,3645) 1, dg%DOF, DTDP*NSPOOLGV, NSPOOLGV, 2
!
!      WRITE(642,2120) 0.0, 1
!      DO J = 1,NE
!        DO K = 1,dg%DOF
!          WRITE(642,2454) J, dg%QX(K,J,1), dg%QY(K,J,1)
!        ENDDO
!      ENDDO
!
!      CLOSE(642)

!.....DG.65.IC = Record of the initial wet/dry flags
!      IF (NOLIFA.GE.2) THEN
!        OPEN(652,FILE=S%DIRNAME//'/'//'DG.65.IC')
!        WRITE(652,3220) RUNDES, RUNID, AGRID
!        WRITE(652,3645) 1, dg%DOF, DTDP*NSPOOLGE, NSPOOLGE, 1
!
!        WRITE(652,2120) 0.0, 1
!       DO J = 1,NE
!          WRITE(652,2455) J,dg%WDFLG(J)
!       ENDDO
!
!        CLOSE(652)
!      ENDIF

2120  FORMAT(2X,E20.10,5X,I10)
2453  FORMAT(2X,I8,2X,E16.8E3)
2454  FORMAT(2X,I8,2(2X,E16.8E3))
2455  FORMAT(2X,I8,2X,I2)
3220  FORMAT(1X,A32,2X,A24,2X,A24)
3645  FORMAT(1X,I10,1X,I10,1X,E15.7,1X,I5,1X,I5)
        
      RETURN
      END
