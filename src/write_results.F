C***********************************************************************
C     
C     SUBROUTINE:  WRITE_RESULTS
C     
C     Taken from the timestep subroutine.  Modified to include
C     interpolation of DG modal degrees of freedom to the nodes.  These
C     multiple nodal values are then averaged to a single value.
C     
C     Aug ??, 2005, sb, Modifications for parallel runs
C     Jan 01, 2007, sb, Files are forced to be written if FORCE_WRITE = .TRUE.
C     
C***********************************************************************

      SUBROUTINE WRITE_RESULTS(IT,FORCE_WRITE)

C.....Use appropriate modules

      USE GLOBAL
      USE DG
      USE HARM
      use sizes, only:  myproc,layers

#ifdef CMPI
      USE MESSENGER_ELEM 
#endif
      
C.....Declare local variables

      INTEGER IT,I,K,Maxp(0:8),Minp(0:8)
      REAL(SZ) AREA, DEPTH, ANGLE_SUM, FH_NL, DP1, DP2, DP3, DP00, ZE00
      LOGICAL FORCE_WRITE
      real(sz) qmaxe,elmaxe,tempx,tempy,rev,Ox,Oy,iota_error,iota_int
      real(sz) iota2_int,iota3_int,iota2_error,iota3_error
      integer imaxze,imaxq, ModetoNode

      Real(SZ),allocatable :: XBCbt(:),YBCbt(:),radial(:),temper(:)

      Allocate ( XBCbt(MNE),YBCbt(MNE),radial(MNE) )
      Allocate ( temper(mne) )
      
C.....Transform from modal coordinates to nodal coordinates and average
C.....to single nodal values 

      ModetoNode = 0

      if (ModetoNode.eq.1) then

         DO I = 1,MNP

            NO_NBORS = EL_COUNT(I)

            AREA_SUM = 0
            ANGLE_SUM = 0
            CEN_SUM = 0
            DO 333 J = 1,NO_NBORS
               NBOR_EL = ELETAB(I,1+J)

               IF(WDFLG(NBOR_EL).EQ.0) CYCLE ! DON'T COUNT DRY ELEMENTS  sb 02/26/07

               DO K = 1,3
                  N1 = NM(NBOR_EL,K)
                  IF (N1.EQ.I) THEN
                     ZE_DG(J) = ZE(1,NBOR_EL,1)
                     QX_DG(J) = QX(1,NBOR_EL,1)
                     QY_DG(J) = QY(1,NBOR_EL,1)
                     HB_DG(J) = HB(1,NBOR_EL,1)

#ifdef TRACE
                     iota_DG(J) = iota(1,NBOR_EL,1)
                     iotaa_DG(J) = iotaa(1,NBOR_EL,1)
#endif

#ifdef CHEM
                     iota_DG(J) = iota(1,NBOR_EL,1)
                     iota2_DG(J) = iota2(1,NBOR_EL,1)
#endif

#ifdef DYNP
                     dynP_DG(J) = dynP(1,NBOR_EL,1)
#endif

#ifdef SED_LAY
                     HB_DG(J) = 0.d0
                     bed_DG(J,:) = bed(1,NBOR_EL,1,:)
                     HB_DG(J) = sum(bed_DG(J,:))
#endif


                     DO KK = 2,DOFH
                        ZE_DG(J) = ZE_DG(J) + PHI_CORNER(KK,K,ph)*ZE(KK,NBOR_EL,1)
                        QX_DG(J) = QX_DG(J) + PHI_CORNER(KK,K,ph)*QX(KK,NBOR_EL,1)
                        QY_DG(J) = QY_DG(J) + PHI_CORNER(KK,K,ph)*QY(KK,NBOR_EL,1)

#ifdef TRACE
                        iota_DG(J) = iota_DG(J) + PHI_CORNER(KK,K,ph)*iota(KK,NBOR_EL,1)
                        iotaa_DG(J) = iotaa_DG(J) + PHI_CORNER(KK,K,ph)*iotaa(KK,NBOR_EL,1)
#endif

#ifdef CHEM
                        iota_DG(J) = iota_DG(J) + PHI_CORNER(KK,K,ph)*iota(KK,NBOR_EL,1)
                        iota2_DG(J) = iota2_DG(J) + PHI_CORNER(KK,K,ph)*iota2(KK,NBOR_EL,1)
#endif

#ifdef DYNP
                        dynP_DG(J) = dynP_DG(J) + PHI_CORNER(KK,K,ph)*dynP(KK,NBOR_EL,1)
#endif

#ifdef SED_LAY
                        do l=1,layers
                           bed_DG(J,l) = bed_DG(J,l) + PHI_CORNER(KK,K,ph)*bed(KK,NBOR_EL,1,l)
                           HB_DG(J) = HB_DG(J) + PHI_CORNER(KK,K,ph)*bed(KK,NBOR_EL,1,l)
                           iotaa_DG(J) = iotaa(1,NBOR_EL,1)
                        enddo
#else

                        HB_DG(J) = HB_DG(J) + PHI_CORNER(KK,K,ph)*HB(KK,NBOR_EL,1)
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

            IF(SEDFLAG.GE.1) DP(I) = 0.D0

            DO J = 1,NO_NBORS
               NBOR_EL = ELETAB(I,1+J)

               IF(WDFLG(NBOR_EL).EQ.0) CYCLE ! DON'T COUNT DRY ELEMENTS  sb 02/26/07

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
               IF(SEDFLAG.GE.1) DP(I) = DP(I) + (AREA/AREA_SUM)*HB_DG(J)
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

            pa = pdg_el(j)
            
            if (pa.eq.0) then

               pa = 1

            endif

            do I = 1,NAGP(pa)

               do k = 1,dofs(j)

                  eta2(j) = eta2(j)+ WAGP(I,pa) * PHI_AREA(K,I,pa) * ZE(K,j,1) * 0.5D0 *AREAS(j)
                  UU2(j)  = UU2(j) + WAGP(I,pa) * PHI_AREA(K,I,pa) * QX(K,j,1) * 0.5D0 *AREAS(j)
                  VV2(j)  = VV2(j) + WAGP(I,pa) * PHI_AREA(K,I,pa) * QY(K,j,1) * 0.5D0 *AREAS(j)


#ifdef TRACE
                  tracer(J)  = tracer(J)  +  WAGP(I,pa) * PHI_AREA(K,I,pa) * iota(K,j,1)  * 0.5D0*AREAS(j)
                  tracer2(J) = tracer2(J) +  WAGP(I,pa) * PHI_AREA(K,I,pa) * iotaa(K,j,1) * 0.5D0*AREAS(j) 
#endif

#ifdef CHEM
                  tracer(J)  = tracer(J)  +  WAGP(I,pa) * PHI_AREA(K,I,pa) * iota(K,j,1) * 0.5D0*AREAS(j)
                  tracer2(J) = tracer2(J) +  WAGP(I,pa) * PHI_AREA(K,I,pa) * iota2(K,j,1) * 0.5D0*AREAS(j)
#endif

#ifdef DYNP
                  dyn_P(J)  = dyn_P(J)  +  WAGP(I,pa) * PHI_AREA(K,I,pa) * dynP(K,j,1) * 0.5D0*AREAS(j)
#endif

#ifdef SED_LAY
                  do l=1,layers

                     bed_int(J,l) = bed_int(J,l) + WAGP(I,pa) 
     &                    *PHI_AREA(K,I,pa)*bed(K,J,1,l)*0.5D0*AREAS(j)

                     DPe(j)  = DPe(j) + WAGP(I,pa)*PHI_AREA(K,I,pa)*bed(K,J,1,l) 
     &                    * 0.5D0 *AREAS(j)
                     
                  enddo

                  tracer2(J) = tracer2(J) +  WAGP(I,pa) * PHI_AREA(K,I,pa) * iotaa2(K,j,1) * 0.5D0*AREAS(j)
#else

                  DPe(j)  = DPe(j) + WAGP(I,pa) * PHI_AREA(K,I,pa) * HB(K,j,1) * 0.5D0 *AREAS(j)
                  tracer2(J) = tracer2(J) +  WAGP(I,pa) * PHI_AREA(K,I,pa) * iotaa2(K,j,1) * 0.5D0*AREAS(j)

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
            do l=1,layers
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
         
         DO I = 1,NAGP(pa)
            
            iota_error = 0.D0
            
            do k = 1,dofs(j)
               
               iota_error = iota_error + (iota(k,j,1)- iotaa(k,j,1)) * PHI_AREA(k,I,pa) 
               
            ENDDO
            
            iota_int = iota_int + iota_error**2.D0 * WAGP(I,pa)*AREAS(J)
            
         ENDDO
         
      enddo
      
      !temperg = 0.D0
      !l2e = 0.D0

      temper =  tracer - tracer2

      temperg = maxval(temper(:))
      l2er_global = iota_int**(1.D0/2.D0)

#ifdef CMPI
 
      Call para_max(temperg)
      Call para_sum(l2er_global)
      
      !call MPI_Allreduce(maxval(temper(:)), temperg, 1, MPI_DOUBLE_PRECISION, MPI_MAX, COMM, ierr);
      !call MPI_Allreduce( iota_int**(1.D0/2.D0), l2er_global, 1, MPI_DOUBLE_PRECISION, MPI_SUM, COMM, ierr);

#endif
      IF  (MYPROC.EQ.0) then
         write(*,*) 'L2nerror('') = ', l2er_global,it,myproc 
         write(*,*) 'linferror('') = ', temperg,it,myproc
      endif

C...  
C...  OUTPUT ELEVATION RECORDING STATION INFORMATION IF NOUTE<>0 AND THE
C.... TIME STEP FALLS WITHIN THE SPECIFIED WINDOW
C...  CALCULATE ELEVATION SOLUTIONS AT STATIONS USING INTERPOLATION
C...  
      
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
c     IF(NCELE.EQ.1) ET00(I)=EE1*STAIE1(I)+EE2*STAIE2(I)
c     &                                              +EE3*STAIE3(I)
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

C...  
C...  OUTPUT VELOCITY RECORDING STATION TIME SERIES INFORMATION IF NOUTV<>0
C.... AND THE TIME STEP FALLS WITHIN THE SPECIFIED WINDOW
C...  CALCULATE VELOCITY SOLUTIONS AT STATIONS USING INTERPOLATION
C...  
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

C...  
C...  OUTPUT CONCENTRATION RECORDING STATION INFORMATION IF NOUTC<>0 AND THE
C.... TIME STEP FALLS WITHIN THE SPECIFIED WINDOW
C...  CALCULATE CONCENTRATION SOLUTIONS AT STATIONS USING INTERPOLATION
C...  

      IF(NOUTC.NE.0) THEN
         IF((IT.GT.NTCYSC).AND.(IT.LE.NTCYFC).OR.FORCE_WRITE) THEN
            NSCOUC=NSCOUC+1
            IF(NSCOUC.EQ.NSPOOLC.OR.FORCE_WRITE) THEN
               DO I=1,NSTAC
                  NM1=NM(NNC(I),1)
                  NM2=NM(NNC(I),2)
                  NM3=NM(NNC(I),3)
                  HH2N1=DP(NM1)+IFNLFA*ETA2(NM1)
                  HH2N2=DP(NM2)+IFNLFA*ETA2(NM2)
                  HH2N3=DP(NM3)+IFNLFA*ETA2(NM3)
                  C1=CH1(NM1)/HH2N1
                  C2=CH1(NM2)/HH2N2
                  C3=CH1(NM3)/HH2N3
                  NC1=NODECODE(NM1)
                  NC2=NODECODE(NM2)
                  NC3=NODECODE(NM3)
                  NCELE=NC1*NC2*NC3
                  IF(NCELE.EQ.1) CC00(I)=C1*STAIC1(I)+C2*STAIC2(I)
     &                 +C3*STAIC3(I)
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

C...  
C...  OUTPUT METEOROLOGICAL RECORDING STATION INFORMATION IF NWS>0 AND THE
C.... TIME STEP FALLS WITHIN THE SPECIFIED WINDOW
C...  CALCULATE METEOROLOGICAL SOLUTIONS AT STATIONS USING INTERPOLATION
C...  

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

C.....Output the gloabl elevation data if NOUTGE ~= 0 and the
C.....time step falls within the specified window

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
 2453                FORMAT(2X,I8,2X,E15.8E3)
                  ENDDO
 
!     cnd...for tecplot results
!     cmm modified 4/6/09
!     cem...modified for cell-centered data, nodes give incorrect results for p>1
                  
                  if (mod(it,NSPOOLGE).eq.0) then
!Casey 120813: Begin the OUT_TEC conditional.
#ifdef OUT_TEC
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
                        write(777,7777)  DBLE(pdg_el(i))
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

 7777                format(9f20.8,i10)
                     do i=1,ne
                        write(777,"(3i12)") nm(i,1), nm(i,2), nm(i,3)
                     enddo
#endif

C.....Write DG.63 results
                     
                     DO J = 1,NE
                        DO K = 1,DOFS(J)
                           WRITE(631,*) J, ZE(K,J,1)
                        ENDDO
                     ENDDO
                     
C.....Write DG.65 results, elemental statuses (wet/dry)
                     
                     IF (NOLIFA.GE.2) THEN
                        WRITE(651,2120) TIME_A, IT
                        DO J = 1,NE
                           WRITE(651,2455) J,WDFLG(J) 
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
                     IF(SEDFLAG.GE.1) THEN
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

C...  
C...  OUTPUT GLOBAL VELOCITY DATA IF NOUTGV<>0 AND THE
C.... TIME STEP FALLS WITHIN THE SPECIFIED WINDOW
C...  
      IF (NOUTGV.NE.0) THEN
         IF ((IT.GT.NTCYSGV).AND.(IT.LE.NTCYFGV).OR.FORCE_WRITE) THEN
            NSCOUGV=NSCOUGV+1
            IF (NSCOUGV.EQ.NSPOOLGV.OR.FORCE_WRITE) THEN
               IF(ABS(NOUTGV).EQ.1) THEN
                  WRITE(64,2120) TIME_A,IT
                  WRITE(641,2120) TIME_A,IT
                  DO I = 1,NE
C.....Had trouble writing out numbers with less than E-99 (EJK)
                     IF (ABS(UU2(I)).LE.(10.0**(-30))) UU2(I) = 0.D0
                     IF (ABS(VV2(I)).LE.(10.0**(-30))) VV2(I) = 0.D0
                     WRITE(64,2454) I,UU2(I),VV2(I)
 2454                FORMAT(2X,I8,2(2X,E15.8E3))
                  ENDDO
                  
C.....Write DG.64 results

                  DO J = 1,NE
                     DO K = 1,DOFS(J)
                        WRITE(641,*) J, QX(K,J,1), QY(K,J,1)
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

C...  
C...  OUTPUT GLOBAL WIND STRESS and atmospheric pressure data IF NOUTGW<>0 AND THE
C.... TIME STEP FALLS WITHIN THE SPECIFIED WINDOW
C...  
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
         open(963,FILE=DIRNAME//'/'//'maxele.63')
         write(963,*) np
         do i=1,np
            write(963,9633) etamax(i)
         enddo
 9633    format(f20.8)
         close(963)
#ifdef OUT_TEC
         write(778,*) 'TITLE = "PADCIRC-DG output"'
         write(778,*) 
     $        'VARIABLES = "x", "y", "maxeta"'
         write(778,*) 'ZONE ZONETYPE=FETRIANGLE ',
     $        'NODES=', np, 
     $        ' ELEMENTS=', ne, 
     $        ' DATAPACKING=POINT ','SOLUTIONTIME=',time_a
         do i=1,np
            write(778,7778) slam(i)/deg2rad, sfea(i)/deg2rad, 
     $           etamax(i)
         enddo
 7778    format(3f20.8)
         do i=1,ne
            write(778,"(3i12)") nm(i,1), nm(i,2), nm(i,3)
         enddo
         close(778)
#endif
      endif

C...  
C...  OUTPUT GLOBAL CONCENTRATION DATA IF NOUTGC<>0 AND THE
C.... TIME STEP FALLS WITHIN THE SPECIFIED WINDOW
C...  
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

C...  
C...  IF IHARIND=1 AND THE TIME STEP FALLS WITHIN THE SPECIFIED WINDOW AND
C...  ON THE SPECIFIED INCREMENT, USE MODEL RESULTS TO UPDATE HARMONIC
C...  ANALYSIS MATRIX AND LOAD VECTORS.  NOTE: AN 8 BYTE RECORD SHOULD BE
C...  USED THROUGHOUT THE HARMONIC ANALYSIS SUBROUTINES, EVEN ON 32 BIT
C...  WORKSTATIONS, SINCE IN THAT CASE THE HARMONIC ANALYSIS IS DONE IN
C...  DOUBLE PRECISION.
C...  
      IF(IHARIND.EQ.1) THEN
         IF((IT.GT.ITHAS).AND.(IT.LE.ITHAF)) THEN
            IF(ICHA.EQ.NHAINC) ICHA=0
            ICHA=ICHA+1
            IF(ICHA.EQ.NHAINC) THEN
C...  
C.....UPDATE THE LHS MATRIX
C...  
               CALL LSQUPDLHS(timeh,IT)
C...  
C.....IF DESIRED COMPUTE ELEVATION STATION INFORMATION AND UPDATE LOAD VECTOR
C...  
               IF(NHASE.EQ.1) THEN
                  DO I=1,NSTAE
C     EE1=ETA2(NM(NNE(I),1))
C     EE2=ETA2(NM(NNE(I),2))
C     EE3=ETA2(NM(NNE(I),3))
C     ET00(I)=EE1*STAIE1(I)+EE2*STAIE2(I)+EE3*STAIE3(I)
                     ET00(I) = ZE(1,NNE(I),1)
                     DO K = 2,DOFH
                        ET00(I) = ET00(I) + ZE(K,NNE(I),1)*PHI_STAE(K,I)
                     ENDDO
                  ENDDO
                  CALL LSQUPDES(ET00,NSTAE)
               ENDIF
C...  
C.....IF DESIRED COMPUTE VELOCITY STATION INFORMATION AND UPDATE LOAD VECTOR
C...  
               IF(NHASV.EQ.1) THEN
                  DO I=1,NSTAV
C     U11=UU2(NM(NNV(I),1))
C     U22=UU2(NM(NNV(I),2))
C     U33=UU2(NM(NNV(I),3))
C     V11=VV2(NM(NNV(I),1))
C     V22=VV2(NM(NNV(I),2))
C     V33=VV2(NM(NNV(I),3))
C     UU00(I)=U11*STAIV1(I)+U22*STAIV2(I)+U33*STAIV3(I)
C     VV00(I)=V11*STAIV1(I)+V22*STAIV2(I)+V33*STAIV3(I)
                     DP1     = DP(NM(NNV(I),1))
                     DP2     = DP(NM(NNV(I),2))
                     DP3     = DP(NM(NNV(I),3))
                     DP00    = DP1*STAIV1(I) + DP2*STAIV2(I) +DP3*STAIV3(I)
                     ZE00    = ZE(1,NNV(I),1)
                     UU00(I) = QX(1,NNV(I),1)
                     VV00(I) = QY(1,NNV(I),1)
                     DO K = 2,DOFH
                        ZE00    = ZE00    + ZE(K,NNV(I),1)*PHI_STAV(K,I)
                        UU00(I) = UU00(I) + QX(K,NNV(I),1)*PHI_STAV(K,I)
                        VV00(I) = VV00(I) + QY(K,NNV(I),1)*PHI_STAV(K,I)
                     ENDDO
                     DEPTH   = ZE00 + DP00
                     FH_NL   = 1.D0/(NLEQ*DEPTH + LEQ)
                     UU00(I) = UU00(I)*FH_NL
                     VV00(I) = VV00(I)*FH_NL
                  ENDDO
                  CALL LSQUPDVS(UU00,VV00,NSTAV)
               ENDIF
C...  
C.....IF DESIRED UPDATE GLOBAL ELEVATION LOAD VECTOR
C...  
               IF(NHAGE.EQ.1) CALL LSQUPDEG(ETA2,NP)
C...  
C.....IF DESIRED UPDATE GLOBAL VELOCITY LOAD VECTOR
C...  
               IF(NHAGV.EQ.1) CALL LSQUPDVG(UU2,VV2,NP)

            ENDIF
         ENDIF

C...  LINES TO COMPUTE MEANS AND VARIANCES

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

C...  
C...  WRITE OUT HOT START INFORMATION IF NHSTAR=1 AND AT CORRECT TIME STEP
C...  NOTE: THE HOT START FILES USE A RECORD LENGTH OF 8 ON BOTH 32 BIT
C.... WORKSTATIONS AND THE 64 BIT CRAY.  THIS IS BECAUSE THE HARMONIC
C.... ANALYSIS IS DONE IN DOUBLE PRECISION (64 BITS) ON WORKSTATIONS.
C...  
      IF(NHSTAR.EQ.1) THEN
         ITEST=(IT/NHSINC)*NHSINC  
         IF(ITEST.EQ.IT) THEN
            IF(IHSFIL.EQ.67) OPEN(67,FILE=DIRNAME//'/'//'fort.67',
     &           ACCESS='DIRECT',RECL=8)
            IF(IHSFIL.EQ.68) OPEN(68,FILE=DIRNAME//'/'//'fort.68',
     &           ACCESS='DIRECT',RECL=8)
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
C...  
C...  IF APPROPRIATE ADD HARMONIC ANALYSIS INFORMATION TO HOT START FILE
C...  
            IF((IHARIND.EQ.1).AND.(IT.GT.ITHAS)) THEN
               WRITE(IHSFIL,REC=IHOTSTP+1) ICHA
               IHOTSTP = IHOTSTP + 1
               CALL HAHOUT(NP,NSTAE,NSTAV,NHASE,NHASV,NHAGE,NHAGV,
     &              IHSFIL,IHOTSTP)
C     
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
c     
            
            IF (C3D) THEN
c     CALL HSTART3D_OUT()
            ENDIF

C...  
C...  CLOSE THE HOT START OUTPUT FILE
C...  
            CLOSE(IHSFIL)
            IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) THEN
               WRITE(6,24541) IHSFIL,IT,TIME_A
            ENDIF
            WRITE(16,24541) IHSFIL,IT,TIME_A
24541       FORMAT(1X,'HOT START OUTPUT WRITTEN TO UNIT ',I2,
     &           ' AT TIME STEP = ',I9,' TIME = ',E15.8)
            IF(IHSFIL.EQ.67) THEN
               IHSFIL=68
            ELSE
               IHSFIL=67
            ENDIF
         ENDIF
      ENDIF
      

C...  FIND AND PRINT TO UNIT 6, THE MAXIMUM ELEVATION, THE MAXIMUM VELOCITY
C.... AND THE NODE NUMBERS AT WHICH THEY OCCUR ON MYPROC=0
C.... IF ELMAX EXCEEDS THRESHOLD, PRINT INFORMATION ON ALL PROCESSORS WHERE
C.... THIS OCCURS 
C...  
      IF(NSCREEN.EQ.1) THEN
         ELMAX=0.0
         VELMAX=0.0
         KEMAX = 0
         KVMAX = 0
c     
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
c     nd
         do i=1,ne
            if (wdflg(i).ne.0) then
               if (abs(ze(1,i,1)).gt.elmaxe) then
                  imaxze=i
                  elmaxe=abs(ze(1,i,1))
               endif
               if (sqrt(qx(1,i,1)**2+qy(1,i,1)**2).gt.qmaxe) then
                  qmaxe=sqrt(qx(1,i,1)**2+qy(1,i,1)**2)
                  imaxq=i
               endif
               
            endif
         enddo
         
         VELMAX=VELMAX**0.5d0
         
#ifdef CMPI
C     sb
C     IF(MYPROC.EQ.0.AND.ELMAX.LT.200.0.AND.KEMAX.GT.0) THEN
         IF( (MYPROC.EQ.0).AND.(ELMAX.LT.200.0).AND.
     &        (MOD(IT,NSCREEN_INC).EQ.0) ) THEN
            if (kemax.eq.0) then 
               kemax = 1 
            endif
       print*,'________________________________',
     &       '__________________________________________'
            print*,''

         WRITE(6,1991) IT,NUMITR,ETA2(KEMAX),KEMAX,VELMAX,KVMAX,
     &           MYPROC
 1991   FORMAT(1X,'TIME STEP =',I8,5X,'ITERATIONS =',I5,
     &           /,2X,'ELMAX = ', E10.4,' AT NODE',I7,
     &           2X,'SPEEDMAX = ',E10.4,' AT NODE',I7,
     &           2X,'ON MYPRnOC = ',I4)
   !write(6,*) 'elmaxe, imaxze, qmaxe, imaxq ',elmaxe,imaxze,
                                !$           qmaxe,imaxq


c$$$  if (chem_flag.eq.1) then
c$$$  print*,'__________________________________________________________________________'
c$$$  print*,'|                       |                                                |'
c$$$  print*,'| Chemistry turned "ON" |',' Maximum mass action =',-minval(MassMax(:)),' |'
c$$$  print*,'|_______________________|________________________________________________|'
c$$$  endif


            if (padapt.eq.1) then

               do k=pl,ph
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
              do k=pl,ph
             print*,'|',k,'         |',Minp(k),'           |'
             print*,'|______________________|________________________|'
              enddo
               print*,''
               print*,'With',ne,'total elements'               
            endif

c$$$  if (slopeflag.eq.4.or.slopeflag.ge.6) then
c$$$  
c$$$  print*,''
c$$$  print*,'p-adaptive slope limiting "ACTIVE."'
c$$$  print*,''
c$$$  print*,lim_count_roll,'elements of a total of'
c$$$  print*,ne,'elements were VERTEX LIMITED'
c$$$  print*,''
c$$$  
c$$$  endif

            print*,'___________________________________',
     &          '_______________________________________'

         ENDIF
C     sb
C     IF(ELMAX.GT.200.0.AND.KEMAX.GT.0) THEN
         IF(ELMAX.GT.200.0) THEN
            if (kemax.eq.0) then 
               kemax = 1 
            endif
            WRITE(6,1993) IT,NUMITR,ETA2(KEMAX),KEMAX,VELMAX,KVMAX,
     &           MYPROC
            WRITE(16,1993) IT,NUMITR,ETA2(KEMAX),KEMAX,VELMAX,KVMAX,
     &           MYPROC
 1993       FORMAT(1X,'TIME STEP =',I8,6X,'ITERATIONS =',I5,
     &      /,2X,'ELMAX = ', E10.4,' AT NODE',I7,
     &      2X,'SPEEDMAX = ',E10.4,' AT NODE',I7,
     &      2X,'ON MYPROC = ',I4,' !!! WARNING - HIGH ELEVATION !!!')
            STOP
         ENDIF


#else
C     sb
C     IF(ELMAX.LT.200.0.AND.KEMAX.GT.0) THEN

         IF(ELMAX.LT.200.0.and.(MOD(IT,NSCREEN_INC).EQ.0) ) THEN
            if (kemax.eq.0) then 
               kemax = 1 
            endif
            print*,'______________________________________',
     &            '____________________________________'
            print*,''
            WRITE(6,1992) IT,NUMITR,ETA2(KEMAX),KEMAX,VELMAX,KVMAX
 1992       FORMAT(1X,'TIME STEP =',I8,5X,'ITERATIONS =',I5,
     &           /,2X,'ELMAX = ', E10.4,' AT NODE',I7,
     &           2X,'SPEEDMAX = ',E10.4,' AT NODE',I7)
                                !write(6,*) 'elmaxe, imaxze, qmaxe, imaxq ',elmaxe,imaxze,
                                !$           qmaxe,imaxq

c$$$  if (chem_flag.eq.1) then
c$$$  print*,'__________________________________________________________________________'
c$$$  print*,'|                       |                                                |'
c$$$  print*,'| Chemistry turned "ON" |',' Maximum mass action =',-minval(MassMax(:)),' |'
c$$$  print*,'|_______________________|________________________________________________|'
c$$$  endif


            if (padapt.eq.1) then
               
               do k=pl,ph
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
               do k=pl,ph
             print*,'|',k,'         |',Minp(k),'           |'
             print*,'|______________________|________________________|'
               enddo
               print*,''
               print*,'With',ne,'total elements'               
            endif

c$$$  if (slopeflag.eq.4.or.slopeflag.ge.6) then
c$$$  
c$$$  print*,''
c$$$  print*,'p-adaptive slope limiting "ACTIVE."'
c$$$  print*,''
c$$$  print*,lim_count_roll,'elements of a total of'
c$$$  print*,ne,'elements were VERTEX LIMITED'
c$$$  print*,''
c$$$  
c$$$  endif

            print*,'___________________________________',
     &          '_______________________________________'

         ENDIF
C     sb
C     IF(ELMAX.GT.200.0.AND.KEMAX.GT.0) THEN
         IF(ELMAX.GT.200.0) THEN
          WRITE(6,1994) IT,NUMITR,ETA2(KEMAX),KEMAX,VELMAX,KVMAX
          WRITE(16,1994) IT,NUMITR,ETA2(KEMAX),KEMAX,VELMAX,KVMAX
 1994       FORMAT(1X,'TIME STEP =',I8,6X,'ITERATIONS =',I5,
     &           /,2X,'ELMAX = ', E10.4,' AT NODE',I7,
     &           2X,'SPEEDMAX = ',E10.4,' AT NODE',I7,
     &           2X,' !!! WARNING - HIGH ELEVATION !!!')
            STOP
         ENDIF
#endif
      ENDIF
      
C.....If applicable write out a DG hot start

      IF ((DGHOT.EQ.1).AND.(MOD(IT,DGHOTSPOOL).EQ.0)) THEN
         
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
#ifdef SED_LAY
         OPEN(290,FILE=DIRNAME//'/'//'Hot_start.290')
#endif

         WRITE(263,*) PH
         WRITE(264,*) PH, PH
         WRITE(214,*) IT

         DO J = 1,MNE
            DO K = 1,DOFH
               WRITE(263,*) ZE(K,J,1)
               WRITE(264,*) QX(K,J,1), QY(K,J,1)
               WRITE(214,*) HB(K,J,1), WDFLG(J)
#ifdef TRACE
               WRITE(288,*) iota(K,J,1)
#endif
#ifdef CHEM
               WRITE(289,*) iota(K,J,1),iota2(K,J,1)
#endif
#ifdef DYNP
               WRITE(291,*) dynP(K,J,1)
#endif
#ifdef SED_LAY
               do ll=1,layers
                  WRITE(290,*) bed(K,J,1,ll)
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

C.....If end of run write out DG data

C     IF (IT.EQ.NT) THEN
C     H_TRI = SQRT((X(1)-X(2))**2 + (Y(1)-Y(2))**2)
C     WRITE(88,*) H_TRI,P,0
C     WRITE(89,*) H_TRI,P,0
C     IF (MYPROC == 0) THEN
C     PRINT*,'T FINAL =',NT*DTDP
C     PRINT*,'DT =',DTDP
C     PRINT*,'NT =',NT
C     ENDIF

C     DO J=1,MNE
C     ZE_DG(1) = 0.D0
C     QX_DG(1) = 0.D0
C     QY_DG(1) = 0.D0
C     
C     DO K=1,DOF
C     ZE_DG(1) = ZE_DG(1) + PHI_CENTER(K)*ZE(K,J,1)
C     QX_DG(1) = QX_DG(1) + PHI_CENTER(K)*QX(K,J,1)
C     QY_DG(1) = QY_DG(1) + PHI_CENTER(K)*QY(K,J,1)
C     
C     WRITE(88,*) ZE(K,J,1),QX(K,J,1),QY(K,J,1)
C     
C     ENDDO
C     
C     WRITE(89,*) ZE_DG(1),QX_DG(1),QY_DG(1)
C     
C     ENDDO
C     ENDIF
C...  
C...  ****************** TIME STEPPING LOOP ENDS HERE ********************
C...  
      RETURN
      END

C***********************************************************************
C
C     SUBROUTINE:  WRITE_DG_IC
C
C     Outputs DG initial conditions
C
C     Written by Shintaro Bunya, Nov. 2005 
C
C***********************************************************************

      SUBROUTINE WRITE_DG_IC()

      USE GLOBAL
      USE DG

C.....DG.63.IC = Record of the initial surface elevation
C      OPEN(632,FILE=DIRNAME//'/'//'DG.63.IC')
C      WRITE(632,3220) RUNDES, RUNID, AGRID
C      WRITE(632,3645) 1, DOF, DTDP*NSPOOLGE, NSPOOLGE, 1

C      WRITE(632,2120) 0.0, 1
C      DO J = 1,NE
C        DO K = 1,DOF
C          WRITE(632,2453) J, ZE(K,J,1)
C        ENDDO
C      ENDDO
      
C      CLOSE(632)

C.....DG.64.IC = Record of the initial discharge
C      OPEN(642,FILE=DIRNAME//'/'//'DG.64.IC')
C      WRITE(642,3220) RUNDES, RUNID, AGRID
C      WRITE(642,3645) 1, DOF, DTDP*NSPOOLGV, NSPOOLGV, 2
C
C      WRITE(642,2120) 0.0, 1
C      DO J = 1,NE
C        DO K = 1,DOF
C          WRITE(642,2454) J, QX(K,J,1), QY(K,J,1)
C        ENDDO
C      ENDDO
C
C      CLOSE(642)

C.....DG.65.IC = Record of the initial wet/dry flags
C      IF (NOLIFA.GE.2) THEN
C        OPEN(652,FILE=DIRNAME//'/'//'DG.65.IC')
C        WRITE(652,3220) RUNDES, RUNID, AGRID
C        WRITE(652,3645) 1, DOF, DTDP*NSPOOLGE, NSPOOLGE, 1
C
C        WRITE(652,2120) 0.0, 1
C       DO J = 1,NE
C          WRITE(652,2455) J,WDFLG(J)
C       ENDDO
C
C        CLOSE(652)
C      ENDIF

2120  FORMAT(2X,E20.10,5X,I10)
2453  FORMAT(2X,I8,2X,E15.8E3)
2454  FORMAT(2X,I8,2(2X,E15.8E3))
2455  FORMAT(2X,I8,2X,I2)
3220  FORMAT(1X,A32,2X,A24,2X,A24)
3645  FORMAT(1X,I10,1X,I10,1X,E15.7,1X,I5,1X,I5)
        
      RETURN
      END
