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

      SUBROUTINE WRITE_RESULTS(s,dg_here,global_here,IT,FORCE_WRITE)

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
      type (dg_type) :: dg_here
      type (global_type) :: global_here


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



      ModetoNode = 1

      if (ModetoNode.eq.1) then
         
         !print*, 'modetonode is 1, entering those loops'

         DO I = 1,S%MNP

            global_here%NO_NBORS = global_here%EL_COUNT(I)

            global_here%AREA_SUM = 0
            ANGLE_SUM = 0
            global_here%CEN_SUM = 0
            DO 333 J = 1,global_here%NO_NBORS
               global_here%NBOR_EL = global_here%ELETAB(I,1+J)

               IF(dg_here%WDFLG(global_here%NBOR_EL).EQ.0) CYCLE ! DON'T COUNT DRY ELEMENTS  sb 02/26/07

               DO K = 1,3
                  global_here%N1 = global_here%NM(global_here%NBOR_EL,K)
                  IF (global_here%N1.EQ.I) THEN
                     global_here%ZE_DG(J) = dg_here%ZE(1,global_here%NBOR_EL,1)
                     global_here%QX_DG(J) = dg_here%QX(1,global_here%NBOR_EL,1)
                     global_here%QY_DG(J) = dg_here%QY(1,global_here%NBOR_EL,1)
                     global_here%HB_DG(J) = dg_here%HB(1,global_here%NBOR_EL,1)

#ifdef TRACE
                     global_here%iota_DG(J) = dg_here%iota(1,global_here%NBOR_EL,1)
                     global_here%iotaa_DG(J) = dg_here%iotaa(1,global_here%NBOR_EL,1)
#endif

#ifdef CHEM
                     global_here%iota_DG(J) = dg_here%iota(1,global_here%NBOR_EL,1)
                     global_here%iota2_DG(J) = dg_here%iota2(1,global_here%NBOR_EL,1)
#endif

#ifdef DYNP
                     global_here%dynP_DG(J) = dg_here%dynP(1,global_here%NBOR_EL,1)
#endif

#ifdef SED_LAY
                     global_here%HB_DG(J) = 0.d0
                     global_here%bed_DG(J,:) = dg_here%bed(1,global_here%NBOR_EL,1,:)
                     global_here%HB_DG(J) = sum(global_here%bed_DG(J,:))
#endif


                     DO KK = 2,dg_here%DOFH
                        global_here%ZE_DG(J) = global_here%ZE_DG(J) + dg_here%PHI_CORNER(KK,K,dg_here%ph)*dg_here%ZE(KK,global_here%NBOR_EL,1)
                        global_here%QX_DG(J) = global_here%QX_DG(J) + dg_here%PHI_CORNER(KK,K,dg_here%ph)*dg_here%QX(KK,global_here%NBOR_EL,1)
                        global_here%QY_DG(J) = global_here%QY_DG(J) + dg_here%PHI_CORNER(KK,K,dg_here%ph)*dg_here%QY(KK,global_here%NBOR_EL,1)

#ifdef TRACE
                        global_here%iota_DG(J) = global_here%iota_DG(J) + dg_here%PHI_CORNER(KK,K,dg_here%ph)*dg_here%iota(KK,global_here%NBOR_EL,1)
                        global_here%iotaa_DG(J) = global_here%iotaa_DG(J) + dg_here%PHI_CORNER(KK,K,dg_here%ph)*dg_here%iotaa(KK,global_here%NBOR_EL,1)
#endif

#ifdef CHEM
                        global_here%iota_DG(J) = global_here%iota_DG(J) + dg_here%PHI_CORNER(KK,K,dg_here%ph)*dg_here%iota(KK,global_here%NBOR_EL,1)
                        global_here%iota2_DG(J) = global_here%iota2_DG(J) + dg_here%PHI_CORNER(KK,K,dg_here%ph)*dg_here%iota2(KK,global_here%NBOR_EL,1)
#endif

#ifdef DYNP
                        global_here%dynP_DG(J) = global_here%dynP_DG(J) + dg_here%PHI_CORNER(KK,K,dg_here%ph)*dg_here%dynP(KK,global_here%NBOR_EL,1)
#endif

#ifdef SED_LAY
                        do l=1,s%layers
                           global_here%bed_DG(J,l) = global_here%bed_DG(J,l) + dg_here%PHI_CORNER(KK,K,dg_here%ph)*dg_here%bed(KK,global_here%NBOR_EL,1,l)
                           global_here%HB_DG(J) = global_here%HB_DG(J) + dg_here%PHI_CORNER(KK,K,dg_here%ph)*dg_here%bed(KK,global_here%NBOR_EL,1,l)
                           global_here%iotaa_DG(J) = dg_here%iotaa(1,global_here%NBOR_EL,1)
                        enddo
#else

                        global_here%HB_DG(J) = global_here%HB_DG(J) + dg_here%PHI_CORNER(KK,K,dg_here%ph)*dg_here%HB(KK,global_here%NBOR_EL,1)
#endif
                        
                     ENDDO
                     AREA = 0.5*global_here%AREAS(global_here%NBOR_EL)
                     global_here%AREA_SUM = global_here%AREA_SUM + AREA
                     GOTO 333
                  ENDIF
               ENDDO
 333        CONTINUE

            global_here%ETA2(I) = 0.D0
            global_here%tracer(I) = 0.D0
            global_here%tracer2(I) = 0.D0
            global_here%UU2(I)  = 0.D0
            global_here%VV2(I)  = 0.D0

            IF(global_here%SEDFLAG.GE.2) global_here%DP(I) = 0.D0

            DO J = 1,global_here%NO_NBORS
               global_here%NBOR_EL = global_here%ELETAB(I,1+J)

               IF(dg_here%WDFLG(global_here%NBOR_EL).EQ.0) CYCLE ! DON'T COUNT DRY ELEMENTS  sb 02/26/07

               AREA = 0.5*global_here%AREAS(global_here%NBOR_EL)/global_here%AREA_SUM
               DEPTH = global_here%ZE_DG(J) + global_here%HB_DG(J)
               FH_NL = 1.D0/(global_here%NLEQ*DEPTH + global_here%LEQ)
               global_here%ETA2(I) = global_here%ETA2(I) + AREA*global_here%ZE_DG(J)

#ifdef TRACE
               global_here%tracer(I)  = global_here%tracer(I)  + AREA*global_here%iota_DG(J)*FH_NL
               global_here%tracer2(I) = global_here%tracer2(I) + AREA*global_here%iotaa_DG(J)*FH_NL
#endif

#ifdef CHEM
               global_here%tracer(I) = global_here%tracer(I) + AREA*global_here%iota_DG(J)*FH_NL
               global_here%tracer2(I) = global_here%tracer2(I) + AREA*global_here%iota2_DG(J)*FH_NL
#endif 

#ifdef DYNP
               global_here%dyn_P(I) = global_here%dyn_P(I) + AREA*global_here%dynP_DG(J)*FH_NL
#endif    
               
               if (global_here%etamax(i).lt.global_here%eta2(i)) global_here%etamax(i)=global_here%eta2(i)
               global_here%UU2(I)  = global_here%UU2(I)  + AREA*global_here%QX_DG(J)*FH_NL
               global_here%VV2(I)  = global_here%VV2(I)  + AREA*global_here%QY_DG(J)*FH_NL
               IF(global_here%SEDFLAG.GE.2) global_here%DP(I) = global_here%DP(I) + (AREA/global_here%AREA_SUM)*global_here%HB_DG(J)
            ENDDO
         ENDDO

      else

         global_here%ETA2 = 0.D0
         global_here%tracer = 0.D0
         global_here%tracer2 = 0.D0
         global_here%UU2  = 0.D0
         global_here%VV2  = 0.D0
         global_here%DPe = 0.D0

         do j=1,global_here%ne

            dg_here%pa = global_here%pdg_el(j)
            
            if (dg_here%pa.eq.0) then

               dg_here%pa = 1

            endif

            do I = 1,dg_here%NAGP(dg_here%pa)

               do k = 1,dg_here%dofs(j)

                  global_here%eta2(j) = global_here%eta2(j)+ dg_here%WAGP(I,dg_here%pa) * dg_here%PHI_AREA(K,I,dg_here%pa) * dg_here%ZE(K,j,1) * 0.5D0 *global_here%AREAS(j)
                  global_here%UU2(j)  = global_here%UU2(j) + dg_here%WAGP(I,dg_here%pa) * dg_here%PHI_AREA(K,I,dg_here%pa) * dg_here%QX(K,j,1) * 0.5D0 *global_here%AREAS(j)
                  global_here%VV2(j)  = global_here%VV2(j) + dg_here%WAGP(I,dg_here%pa) * dg_here%PHI_AREA(K,I,dg_here%pa) * dg_here%QY(K,j,1) * 0.5D0 *global_here%AREAS(j)


#ifdef TRACE
                  global_here%tracer(J)  = global_here%tracer(J)  +  dg_here%WAGP(I,dg_here%pa) * dg_here%PHI_AREA(K,I,dg_here%pa) * dg_here%iota(K,j,1)  * 0.5D0*global_here%AREAS(j)
                  global_here%tracer2(J) = global_here%tracer2(J) +  dg_here%WAGP(I,dg_here%pa) * dg_here%PHI_AREA(K,I,dg_here%pa) * dg_here%iotaa(K,j,1) * 0.5D0*global_here%AREAS(j) 
#endif

#ifdef CHEM
                  global_here%tracer(J)  = global_here%tracer(J)  +  dg_here%WAGP(I,dg_here%pa) * dg_here%PHI_AREA(K,I,dg_here%pa) * dg_here%iota(K,j,1) * 0.5D0*global_here%AREAS(j)
                  global_here%tracer2(J) = global_here%tracer2(J) +  dg_here%WAGP(I,dg_here%pa) * dg_here%PHI_AREA(K,I,dg_here%pa) * dg_here%iota2(K,j,1) * 0.5D0*global_here%AREAS(j)
#endif

#ifdef DYNP
                  global_here%dyn_P(J)  = global_here%dyn_P(J)  +  dg_here%WAGP(I,dg_here%pa) * dg_here%PHI_AREA(K,I,dg_here%pa) * dg_here%dynP(K,j,1) * 0.5D0*global_here%AREAS(j)
#endif

#ifdef SED_LAY
                  do l=1,s%layers

                     global_here%bed_int(J,l) = global_here%bed_int(J,l) + dg_here%WAGP(I,dg_here%pa) 
&                    *dg_here%PHI_AREA(K,I,dg_here%pa)*dg_here%bed(K,J,1,l)*0.5D0*global_here%AREAS(j)

                     global_here%DPe(j)  = global_here%DPe(j) + dg_here%WAGP(I,dg_here%pa)*dg_here%PHI_AREA(K,I,dg_here%pa)*dg_here%bed(K,J,1,l) 
&                    * 0.5D0 *global_here%AREAS(j)
                     
                  enddo

                  global_here%tracer2(J) = global_here%tracer2(J) +  dg_here%WAGP(I,dg_here%pa) * dg_here%PHI_AREA(K,I,dg_here%pa) * dg_here%iotaa2(K,j,1) * 0.5D0*global_here%AREAS(j)
#else

                  global_here%DPe(j)  = global_here%DPe(j) + dg_here%WAGP(I,dg_here%pa) * dg_here%PHI_AREA(K,I,dg_here%pa) * dg_here%HB(K,j,1) * 0.5D0 *global_here%AREAS(j)
                  global_here%tracer2(J) = global_here%tracer2(J) +  dg_here%WAGP(I,dg_here%pa) * dg_here%PHI_AREA(K,I,dg_here%pa) * dg_here%iotaa2(K,j,1) * 0.5D0*global_here%AREAS(j)

#endif
                  
               enddo

            enddo

            DEPTH =  (global_here%eta2(j) + global_here%dpe(j)) 
            FH_NL = 1.D0/(global_here%NLEQ*DEPTH + global_here%LEQ)
            global_here%eta2(j) = global_here%eta2(j)/(global_here%AREAS(J))
            global_here%UU2(j)  = global_here%UU2(j) *FH_NL
            global_here%VV2(j)  = global_here%VV2(j) *FH_NL
            global_here%dpe(j)  = global_here%dpe(j)/(global_here%AREAS(J))

#ifdef SED_LAY
            do l=1,s%layers
               global_here%bed_int(j,l) = global_here%bed_int(j,l)/(global_here%AREAS(J))
            enddo
#endif

#ifdef TRACE            
            global_here%tracer(J)  = global_here%tracer(J)  * FH_NL
            global_here%tracer2(J) = global_here%tracer2(J) * FH_NL  
#endif

#ifdef CHEM
            global_here%tracer(J)  = global_here%tracer(J)  * FH_NL
            global_here%tracer2(J) = global_here%tracer2(J) * FH_NL
#endif

#ifdef DYNP            
            global_here%dyn_P(J)  = global_here%dyn_P(J)  * FH_NL
#endif
            
         enddo

      endif

      iota_int = 0.D0
      
      do j = 1,global_here%ne
         
         DO I = 1,dg_here%NAGP(dg_here%pa)
            
            iota_error = 0.D0
            
            do k = 1,dg_here%dofs(j)
               
               iota_error = iota_error + (dg_here%iota(k,j,1)- dg_here%iotaa(k,j,1)) * dg_here%PHI_AREA(k,I,dg_here%pa) 
               
            ENDDO
            
            iota_int = iota_int + iota_error**2.D0 * dg_here%WAGP(I,dg_here%pa)*global_here%AREAS(J)
            
         ENDDO
         
      enddo

!...  
!...  OUTPUT ELEVATION RECORDING STATION INFORMATION IF global_here%NOUTE<>0 AND THE
!.... TIME STEP FALLS WITHIN THE SPECIFIED WINDOW
!...  CALCULATE ELEVATION SOLUTIONS AT STATIONS USING INTERPOLATION
!...  
      
      IF(global_here%NOUTE.ne.0) THEN
         IF((IT.GT.global_here%NTCYSE).AND.(IT.LE.global_here%NTCYFE).OR.FORCE_WRITE) THEN
            global_here%NSCOUE=global_here%NSCOUE+1
            IF(global_here%NSCOUE.EQ.global_here%NSPOOLE.OR.FORCE_WRITE) THEN
               DO I=1,global_here%NSTAE           
                                !global_here%EE1=global_here%ETA2(global_here%NM(global_here%NNE(I),1))
                                !global_here%EE2=global_here%ETA2(global_here%NM(global_here%NNE(I),2))
                                !global_here%EE3=global_here%ETA2(global_here%NM(global_here%NNE(I),3))

                  global_here%ee1=global_here%eta2(global_here%nne(i))

                  global_here%NC1=global_here%NODECODE(global_here%NM(global_here%NNE(I),1))
                  global_here%NC2=global_here%NODECODE(global_here%NM(global_here%NNE(I),2))
                  global_here%NC3=global_here%NODECODE(global_here%NM(global_here%NNE(I),3))
                  global_here%NCELE=global_here%NC1*global_here%NC2*global_here%NC3
!     IF(global_here%NCELE.EQ.1) global_here%ET00(I)=global_here%EE1*global_here%STAIE1(I)+global_here%EE2*global_here%STAIE2(I)
!     &                                              +global_here%EE3*global_here%STAIE3(I)
                  if (global_here%ncele.eq.1) then
                     global_here%ET00(I)=global_here%ee1
                  endif
                  IF(global_here%NCELE.EQ.0) global_here%ET00(I)=-99999.
                  
                  global_here%ee1 = global_here%dpe(global_here%nne(i))

                                !global_here%EE1=global_here%DP(global_here%NM(global_here%NNE(I),1))
                                !global_here%EE2=global_here%DP(global_here%NM(global_here%NNE(I),2))
                                !global_here%EE3=global_here%DP(global_here%NM(global_here%NNE(I),3))
                  IF(global_here%NCELE.EQ.1) global_here%BT00(I)= global_here%ee1

                                !global_here%EE1*global_here%STAIE1(I)+global_here%EE2*global_here%STAIE2(I)
                                !&                 +global_here%EE3*global_here%STAIE3(I)
                  IF(global_here%NCELE.EQ.0) global_here%BT00(I)=-99999.
               END DO
               IF(ABS(global_here%NOUTE).EQ.1) THEN
                  WRITE(s%fort61unit,2120) global_here%time_A,IT
                  WRITE(S%FORT82UNIT,2120) global_here%time_A,IT
                  DO I=1,global_here%NSTAE
                     WRITE(s%fort61unit,2453) I,global_here%ET00(I)
                     WRITE(S%FORT82UNIT,2453) I,global_here%BT00(I)
                  END DO
                  global_here%IESTP = global_here%IESTP+1+global_here%NSTAE
               ENDIF
               IF(ABS(global_here%NOUTE).EQ.2) THEN
                  WRITE(s%fort61unit,REC=global_here%IESTP+1) global_here%time_A
                  WRITE(s%fort61unit,REC=global_here%IESTP+2) IT
                  WRITE(S%FORT82UNIT,REC=global_here%IESTP+1) global_here%time_A
                  WRITE(S%FORT82UNIT,REC=global_here%IESTP+2) IT
                  global_here%IESTP = global_here%IESTP + 2
                  DO I=1,global_here%NSTAE
                     WRITE(s%fort61unit,REC=global_here%IESTP+I) global_here%ET00(I)
                     WRITE(S%FORT82UNIT,REC=global_here%IESTP+I) global_here%BT00(I)
                  END DO
                  global_here%IESTP = global_here%IESTP + global_here%NSTAE
               ENDIF
               global_here%NSCOUE=0
            ENDIF
         ENDIF
         IF(IT.EQ.global_here%NTCYFE) CLOSE(s%fort61unit)
         IF(IT.EQ.global_here%NTCYFE) CLOSE(s%fort82unit)
      ENDIF

!...  
!...  OUTPUT VELOCITY RECORDING STATION TIME SERIES INFORMATION IF global_here%NOUTV<>0
!.... AND THE TIME STEP FALLS WITHIN THE SPECIFIED WINDOW
!...  CALCULATE VELOCITY SOLUTIONS AT STATIONS USING INTERPOLATION
!...  
      IF(global_here%NOUTV.ne.0) THEN
         IF((IT.GT.global_here%NTCYSV).AND.(IT.LE.global_here%NTCYFV).OR.FORCE_WRITE) THEN
            global_here%NSCOUV=global_here%NSCOUV+1
            IF(global_here%NSCOUV.EQ.global_here%NSPOOLV.OR.FORCE_WRITE) THEN
               DO I=1,global_here%NSTAV
                  global_here%U11=global_here%UU2(global_here%nnv(I))
                                !global_here%U22=global_here%UU2(global_here%NM(global_here%NNV(I),2))
                                !global_here%U33=global_here%UU2(global_here%NM(global_here%NNV(I),3))
                  global_here%V11=global_here%VV2(global_here%NNV(I))
                                !global_here%V22=global_here%VV2(global_here%NM(global_here%NNV(I),2))
                                !global_here%V33=global_here%VV2(global_here%NM(global_here%NNV(I),3))
                  global_here%UU00(I)= global_here%u11  !global_here%U11*global_here%STAIV1(I)+global_here%U22*global_here%STAIV2(I)+global_here%U33*global_here%STAIV3(I)
                  global_here%VV00(I)= global_here%v11  !global_here%V11*global_here%STAIV1(I)+global_here%V22*global_here%STAIV2(I)+global_here%V33*global_here%STAIV3(I)
               END DO
               IF(ABS(global_here%NOUTV).EQ.1) THEN
                  WRITE(s%fort62unit,2120) global_here%time_A,IT
                  DO I=1,global_here%NSTAV
                     WRITE(s%fort62unit,2454) I,global_here%UU00(I),global_here%VV00(I)
                  END DO
                  global_here%IVSTP = global_here%IVSTP+1+global_here%NSTAV
               ENDIF
               IF(ABS(global_here%NOUTV).EQ.2) THEN
                  WRITE(s%fort62unit,REC=global_here%IVSTP+1) global_here%time_A
                  WRITE(s%fort62unit,REC=global_here%IVSTP+2) IT
                  global_here%IVSTP = global_here%IVSTP + 2
                  DO I=1,global_here%NSTAV
                     WRITE(s%fort62unit,REC=global_here%IVSTP+2*I-1) global_here%UU00(I)
                     WRITE(s%fort62unit,REC=global_here%IVSTP+2*I) global_here%VV00(I)
                  END DO
                  global_here%IVSTP = global_here%IVSTP + 2*global_here%NSTAV
               ENDIF
               global_here%NSCOUV=0
            ENDIF
         ENDIF
         IF(IT.EQ.global_here%NTCYFV) CLOSE(s%fort62unit)
      ENDIF

!...  
!...  OUTPUT CONCENTRATION RECORDING STATION INFORMATION IF global_here%NOUTC<>0 AND THE
!.... TIME STEP FALLS WITHIN THE SPECIFIED WINDOW
!...  CALCULATE CONCENTRATION SOLUTIONS AT STATIONS USING INTERPOLATION
!...  

      IF(global_here%NOUTC.ne.0) THEN
         IF((IT.GT.global_here%NTCYSC).AND.(IT.LE.global_here%NTCYFC).OR.FORCE_WRITE) THEN
            global_here%NSCOUC=global_here%NSCOUC+1
            IF(global_here%NSCOUC.EQ.global_here%NSPOOLC.OR.FORCE_WRITE) THEN
               DO I=1,global_here%NSTAC
                  global_here%NM1=global_here%NM(global_here%NNC(I),1)
                  global_here%NM2=global_here%NM(global_here%NNC(I),2)
                  global_here%NM3=global_here%NM(global_here%NNC(I),3)
                  global_here%HH2N1=global_here%dpe(global_here%nnc(i))+global_here%IFNLFA*global_here%ETA2(global_here%nnv(i))  !global_here%DP(global_here%NM1)+global_here%IFNLFA*global_here%ETA2(global_here%NM1)
                  !global_here%HH2N2=  global_here%DP(global_here%NM2)+global_here%IFNLFA*global_here%ETA2(global_here%NM2)
                  !global_here%HH2N3=  global_here%DP(global_here%NM3)+global_here%IFNLFA*global_here%ETA2(global_here%NM3)
                  global_here%C1=global_here%CH1(global_here%NM1)/global_here%HH2N1
                  global_here%C2=0.D0
                  global_here%C3=0.D0
                  !global_here%C2=global_here%CH1(global_here%NM2)/global_here%HH2N2
                  !global_here%C3=global_here%CH1(global_here%NM3)/global_here%HH2N3
                  global_here%NC1=global_here%NODECODE(global_here%NM1)
                  global_here%NC2=global_here%NODECODE(global_here%NM2)
                  global_here%NC3=global_here%NODECODE(global_here%NM3)
                  global_here%NCELE=global_here%NC1*global_here%NC2*global_here%NC3
                  IF(global_here%NCELE.EQ.1) global_here%CC00(I)=global_here%C1*global_here%STAIC1(I)+global_here%C2*global_here%STAIC2(I)&
                       +global_here%C3*global_here%STAIC3(I)
                  IF(global_here%NCELE.EQ.0) global_here%CC00(I)=-99999.
               END DO
               IF(ABS(global_here%NOUTC).EQ.1) THEN
                  WRITE(S%FORT81UNIT,2120) global_here%time_A,IT
                  DO I=1,global_here%NSTAC
                     WRITE(S%FORT81UNIT,2453) I,global_here%CC00(I)
                  END DO
                  global_here%ICSTP = global_here%ICSTP+1+global_here%NSTAC
               ENDIF
               IF(ABS(global_here%NOUTC).EQ.2) THEN
                  WRITE(S%FORT81UNIT,REC=global_here%ICSTP+1) global_here%time_A
                  WRITE(S%FORT81UNIT,REC=global_here%ICSTP+2) IT
                  global_here%ICSTP = global_here%ICSTP + 2
                  DO I=1,global_here%NSTAC
                     WRITE(S%FORT81UNIT,REC=global_here%ICSTP+I) global_here%CC00(I)
                  END DO
                  global_here%ICSTP = global_here%ICSTP + global_here%NSTAC
               ENDIF
               global_here%NSCOUC=0
            ENDIF
         ENDIF
         IF(IT.EQ.global_here%NTCYFC) CLOSE(s%fort81unit)
      ENDIF

!...  
!...  OUTPUT METEOROLOGICAL RECORDING STATION INFORMATION IF global_here%NWS>0 AND THE
!.... TIME STEP FALLS WITHIN THE SPECIFIED WINDOW
!...  CALCULATE METEOROLOGICAL SOLUTIONS AT STATIONS USING INTERPOLATION
!...  

      IF((global_here%NWS.ne.0).AND.(global_here%NOUTM.ne.0)) THEN
         IF((IT.GT.global_here%NTCYSM).AND.(IT.LE.global_here%NTCYFM).OR.FORCE_WRITE) THEN
            global_here%NSCOUM=global_here%NSCOUM+1
            IF(global_here%NSCOUM.EQ.global_here%NSPOOLM.OR.FORCE_WRITE) THEN
               DO I=1,global_here%NSTAM
                  global_here%NM1=global_here%NM(global_here%NNM(I),1)
                  global_here%NM2=global_here%NM(global_here%NNM(I),2)
                  global_here%NM3=global_here%NM(global_here%NNM(I),3)
                  global_here%U11=global_here%wvnxout(global_here%NM1)
                  global_here%U22=global_here%wvnxout(global_here%NM2)
                  global_here%U33=global_here%wvnxout(global_here%NM3)
                  global_here%V11=global_here%wvnyout(global_here%NM1)
                  global_here%V22=global_here%wvnyout(global_here%NM2)
                  global_here%V33=global_here%wvnyout(global_here%NM3)
                  global_here%P11=global_here%PR2(global_here%NM1)
                  global_here%P22=global_here%PR2(global_here%NM2)
                  global_here%P33=global_here%PR2(global_here%NM3)
                  global_here%RMU00(I)=global_here%U11*global_here%STAIM1(I)+global_here%U22*global_here%STAIM2(I)+global_here%U33*global_here%STAIM3(I)
                  global_here%RMV00(I)=global_here%V11*global_here%STAIM1(I)+global_here%V22*global_here%STAIM2(I)+global_here%V33*global_here%STAIM3(I)
                  global_here%RMP00(I)=global_here%P11*global_here%STAIM1(I)+global_here%P22*global_here%STAIM2(I)+global_here%P33*global_here%STAIM3(I)
               END DO
               IF(ABS(global_here%NOUTM).EQ.1) THEN
                  WRITE(S%FORT71UNIT,2120) global_here%time_A,IT
                  WRITE(S%FORT72UNIT,2120) global_here%time_A,IT
                  DO I=1,global_here%NSTAM
                     WRITE(S%FORT71UNIT,2453) I,global_here%RMP00(I)
                     WRITE(S%FORT72UNIT,2454) I,global_here%RMU00(I),global_here%RMV00(I)
                  END DO
                  global_here%IPSTP=global_here%IPSTP+1+global_here%NSTAM
                  global_here%IWSTP=global_here%IWSTP+1+global_here%NSTAM
               ENDIF
               IF(ABS(global_here%NOUTM).EQ.2) THEN
                  WRITE(S%FORT71UNIT,REC=global_here%IPSTP+1) global_here%time_A
                  WRITE(S%FORT71UNIT,REC=global_here%IPSTP+2) IT
                  WRITE(S%FORT72UNIT,REC=global_here%IWSTP+1) global_here%time_A
                  WRITE(S%FORT72UNIT,REC=global_here%IWSTP+2) IT
                  global_here%IPSTP=global_here%IPSTP+2
                  global_here%IWSTP=global_here%IWSTP+2
                  DO I=1,global_here%NSTAM
                     WRITE(S%FORT71UNIT,REC=global_here%IPSTP+I) global_here%RMP00(I)
                     WRITE(S%FORT72UNIT,REC=global_here%IWSTP+2*I-1) global_here%RMU00(I)
                     WRITE(S%FORT72UNIT,REC=global_here%IWSTP+2*I) global_here%RMV00(I)
                  END DO
                  global_here%IPSTP=global_here%IPSTP+global_here%NSTAM
                  global_here%IWSTP=global_here%IWSTP+2*global_here%NSTAM
               ENDIF
               global_here%NSCOUM=0
            ENDIF
         ENDIF
         IF(IT.EQ.global_here%NTCYFM) THEN
            CLOSE(s%fort71unit)
            CLOSE(s%fort72unit)
         ENDIF
      ENDIF

!.....Output the gloabl elevation data if global_here%NOUTGE ~= 0 and the
!.....time step falls within the specified window

      IF (global_here%NOUTGE.ne.0) THEN
         IF ((IT.GT.global_here%NTCYSGE).AND.(IT.LE.global_here%NTCYFGE).OR.FORCE_WRITE) THEN
            global_here%NSCOUGE = global_here%NSCOUGE + 1
            IF (global_here%NSCOUGE.EQ.global_here%NSPOOLGE.OR.FORCE_WRITE) THEN
               IF (ABS(global_here%NOUTGE).EQ.1) THEN
                  WRITE(s%fort63unit,2120) global_here%TIME_A, IT
                  WRITE(S%FORT88UNIT,2120) global_here%TIME_A, IT
                  WRITE(S%FORT89UNIT,2120) global_here%TIME_A, IT
!                  WRITE(4441,2120) global_here%TIME_A, IT
                                !WRITE(895,21 20) global_here%TIME_A, IT
                  WRITE(s%dg63unit,2120) global_here%TIME_A, IT
 2120             FORMAT(2X,E20.10,5X,I10)
                  DO I=1,global_here%NE
                     IF (ABS(global_here%ETA2(I)).LE.(10.0**(-30))) global_here%ETA2(I) = 0.D0
                     WRITE(s%fort63unit,2453) I,global_here%ETA2(I)
                                !IF (global_here%NODECODE(I).EQ.0) WRITE(63,2453) I,-99999.
                     IF (ABS(global_here%tracer(I)).LE.(10.0**(-30))) global_here%tracer(I) = 0.D0
                     WRITE(S%FORT88UNIT,2453) I,global_here%tracer(I)
                                !IF (global_here%NODECODE(I).EQ.0) WRITE(S%FORT88UNIT,2453) I,.
                     IF (ABS(global_here%tracer2(I)).LE.(10.0**(-30))) global_here%tracer2(I) = 0.D0
                     WRITE(S%FORT89UNIT,2453) I,global_here%tracer2(I)
                                !IF (global_here%NODECODE(I).EQ.0) WRITE(S%FORT89UNIT,2453) I,0.\
                     IF (ABS(global_here%dpe(I)).LE.(10.0**(-30))) global_here%dpe(I) = 0.D0
 !                    write(4441,2453) I,global_here%dpe(I)
 2453                FORMAT(2X,I8,2X,E16.8E3)
                  ENDDO
 
!     cnd...for tecplot results
!     cmm modified 4/6/09
!     cem...modified for cell-centered data, nodes give incorrect results for p>1
                  
                  if (mod(it,global_here%NSPOOLGE).eq.0) then
!Casey 120813: Begin the OUT_TEC conditional.

#ifdef OUT_TEC
                     
                     if (ModetoNode.eq.0) then

           write(777,*) 'ZONE ZONETYPE=FETRIANGLE ',
     $     'NODES=', global_here%np, 
     $     ' ELEMENTS=', global_here%ne, 
     $     ' DATAPACKING=BLOCK ',
     $     'SOLUTIONTIME=',global_here%time_a,
     $     'VARLOCATION=([3,4,5,6,7,8,10,11,12,13,14]=CELLCENTERED)'
                     do i=1,global_here%np
                        write(777,7777)  global_here%x(i)
                     enddo
                     do i=1,global_here%np
                        write(777,7777)  global_here%y(i)
                     enddo
                     do i=1,global_here%ne 
                        write(777,7777)  global_here%dpe(i)
                     enddo
                     do i=1,global_here%ne 
                        write(777,7777)  global_here%eta2(i)
                     enddo
                     do i=1,global_here%ne 
                        write(777,7777)  global_here%eta2(i)+global_here%dpe(i)
                     enddo
                     do i=1,global_here%ne 
                        write(777,7777)  global_here%uu2(i)
                     enddo
                     do i=1,global_here%ne 
                        write(777,7777)  global_here%vv2(i)
                     enddo
                     do i=1,global_here%ne 
                        write(777,7777)  sqrt(global_here%uu2(i)**2+global_here%vv2(i)**2)
                     enddo
                     do i=1,global_here%np 
                        write(777,7777)  sqrt(global_here%wsx2(i)**2+global_here%wsy2(i)**2)
                     enddo
                     do i=1,global_here%ne 
                        write(777,7777)  global_here%entrop(4,i) !DBLE(global_here%pdg_el(i))
                     enddo
                     do i=1,global_here%ne 
                        write(777,7777)  global_here%tracer(i)
                     enddo
                     do i=1,global_here%ne 
                        write(777,7777)  global_here%tracer2(i)
                     enddo
                     do i=1,global_here%ne 
                        write(777,7777)  abs(global_here%tracer(i)+global_here%tracer2(i))
                     enddo
                     do i=1,global_here%ne 
                        write(777,7777)  abs(global_here%tracer(i)-global_here%tracer2(i))    
                     enddo

                  else

                 write(777,*) 'ZONE ZONETYPE=FETRIANGLE ',
     $                'NODES=', global_here%np, 
     $                ' ELEMENTS=', global_here%ne, 
     $                ' DATAPACKING=POINT ','SOLUTIONTIME=',global_here%time_a
                 do i=1,global_here%np
                    if (global_here%ics.eq.2) then
                       write(777,7777) global_here%slam(i)/deg2rad, global_here%sfea(i)/deg2rad, 
!                 write(777,7777) global_here%x(i), global_here%y(i), 
     $                   global_here%dp(i), global_here%eta2(i), global_here%eta2(i)+global_here%dp(i),global_here%uu2(i),global_here%vv2(i),
     $           sqrt(global_here%uu2(i)**2+global_here%vv2(i)**2),sqrt(global_here%wsx2(i)**2+global_here%wsy2(i)**2),
     $                   s%myproc
                    else
                       write(777,7777) global_here%x(i), global_here%y(i), 
     $                   global_here%dp(i), global_here%eta2(i), global_here%eta2(i)+global_here%dp(i),global_here%uu2(i),global_here%vv2(i),
     $           sqrt(global_here%uu2(i)**2+global_here%vv2(i)**2),sqrt(global_here%wsx2(i)**2+global_here%wsy2(i)**2),
     $                   s%myproc
                    endif

                 enddo
              endif


 7777         format(9f20.8,i10)
              do i=1,global_here%ne
                 write(777,"(3i12)") global_here%nm(i,1), global_here%nm(i,2), global_here%nm(i,3)
              enddo


#endif

!.....Write DG.63 results
                     
                     DO J = 1,global_here%NE
                        DO K = 1,dg_here%DOFS(J)
                           !srb: add bathymetry output to work with Dam's plotting scripts
                           WRITE(s%dg63unit,*) J, dg_here%ZE(K,J,1), dg_here%HB(K,J,1)  
                        ENDDO
                     ENDDO
                     
!.....Write DG.65 results, elemental statuses (wet/dry)
                     
                     IF (global_here%NOLIFA.GE.2) THEN
                        WRITE(s%dg65unit,2120) global_here%TIME_A, IT
                        DO J = 1,global_here%NE
                           WRITE(s%dg65unit,2455) J,dg_here%WDFLG(J) 
 2455                      FORMAT(2X,I8,2X,I2)
                        ENDDO
                     ENDIF
                     
                     global_here%IGEP = global_here%IGEP + 1 + global_here%NP
                  ENDIF
                  IF (ABS(global_here%NOUTGE).EQ.2) THEN
                     WRITE(s%fort63unit,REC=global_here%IGEP+1) global_here%TIME_A
                     WRITE(s%fort63unit,REC=global_here%IGEP+2) IT
                     global_here%IGEP = global_here%IGEP + 2
                     DO I = 1,global_here%NP
                        WRITE(s%fort63unit,REC=global_here%IGEP+I)global_here%ETA2(I)
                                !IF (global_here%NODECODE(I).EQ.0) WRITE(63,REC=global_here%IGEP+I) -99999.
                     ENDDO
                     IF(global_here%SEDFLAG.GE.2) THEN
                        WRITE(S%FORT84UNIT,REC=global_here%IGEP+1) global_here%TIME_A
                        WRITE(S%FORT84UNIT,REC=global_here%IGEP+2) IT
                        DO I = 1,global_here%NP
                           WRITE(S%FORT84UNIT,REC=global_here%IGEP+I) global_here%DP(I)
                        ENDDO
                     ENDIF
                     global_here%IGEP = global_here%IGEP + global_here%NP
                  ENDIF
                  global_here%NSCOUGE = 0
               ENDIF
            ENDIF
            IF (IT.EQ.global_here%NTCYFGE) CLOSE(s%fort63unit)
            IF (IT.EQ.global_here%NTCYFGE) CLOSE(s%dg63unit)
            IF (IT.EQ.global_here%NTCYFGE.AND.global_here%NOLIFA.GE.2) CLOSE(s%dg65unit)
            IF (IT.EQ.global_here%NTCYFGE.AND.global_here%SEDFLAG.GE.1) CLOSE(s%fort84unit)

         endif
      ENDIF

!...  
!...  OUTPUT GLOBAL VELOCITY DATA IF global_here%NOUTGV<>0 AND THE
!.... TIME STEP FALLS WITHIN THE SPECIFIED WINDOW
!...  
      IF (global_here%NOUTGV.ne.0) THEN
         IF ((IT.GT.global_here%NTCYSGV).AND.(IT.LE.global_here%NTCYFGV).OR.FORCE_WRITE) THEN
            global_here%NSCOUGV=global_here%NSCOUGV+1
            IF (global_here%NSCOUGV.EQ.global_here%NSPOOLGV.OR.FORCE_WRITE) THEN
               IF(ABS(global_here%NOUTGV).EQ.1) THEN
                  WRITE(s%fort64unit,2120) global_here%TIME_A,IT
                  WRITE(s%dg64unit,2120) global_here%TIME_A,IT
                  DO I = 1,global_here%NE
!.....Had trouble writing out numbers with less than E-99 (EJK)
                     IF (ABS(global_here%UU2(I)).LE.(10.0**(-30))) global_here%UU2(I) = 0.D0
                     IF (ABS(global_here%VV2(I)).LE.(10.0**(-30))) global_here%VV2(I) = 0.D0
                     WRITE(s%fort64unit,2454) I,global_here%UU2(I),global_here%VV2(I)
 2454                FORMAT(2X,I8,2(2X,E16.8E3))
                  ENDDO
                  
!.....Write DG.64 results

                  DO J = 1,global_here%NE
                     DO K = 1,dg_here%DOFS(J)
                        WRITE(s%dg64unit,*) J, dg_here%QX(K,J,1), dg_here%QY(K,J,1)
                     ENDDO
                  ENDDO
                  
                  global_here%IGVP = global_here%IGVP + 1 + global_here%NP
               ENDIF
               IF (ABS(global_here%NOUTGV).EQ.2) THEN
                  WRITE(s%fort64unit,REC=global_here%IGVP+1) global_here%time_A
                  WRITE(s%fort64unit,REC=global_here%IGVP+2) IT
                  global_here%IGVP = global_here%IGVP + 2
                  DO I=1,global_here%NP
                     WRITE(s%fort64unit,REC=global_here%IGVP+2*I-1) global_here%UU2(I)
                     WRITE(s%fort64unit,REC=global_here%IGVP+2*I) global_here%VV2(I)
                  ENDDO
                  global_here%IGVP = global_here%IGVP + 2*global_here%NP
               ENDIF
               global_here%NSCOUGV=0
            ENDIF
         ENDIF
         IF(IT.EQ.global_here%NTCYFGV) CLOSE(s%fort64unit)
         IF(IT.EQ.global_here%NTCYFGV) CLOSE(s%dg64unit)
      ENDIF

!...  
!...  OUTPUT GLOBAL WIND STRESS and atmospheric pressure data IF global_here%NOUTGW<>0 AND THE
!.... TIME STEP FALLS WITHIN THE SPECIFIED WINDOW
!...  
      IF((global_here%NWS.ne.0).AND.(global_here%NOUTGW.ne.0)) THEN
         IF((IT.GT.global_here%NTCYSGW).AND.(IT.LE.global_here%NTCYFGW).OR.FORCE_WRITE) THEN
            global_here%NSCOUGW=global_here%NSCOUGW+1
            IF(global_here%NSCOUGW.EQ.global_here%NSPOOLGW.OR.FORCE_WRITE) THEN
               IF(ABS(global_here%NOUTGW).EQ.1) THEN
                  write(s%fort73unit,2120) global_here%time_A,it
                  WRITE(S%FORT74UNIT,2120) global_here%time_A,IT
                  DO I=1,global_here%NP
                     write(s%fort73unit,2453) i,global_here%pr2(i)
                     WRITE(S%FORT74UNIT,2454) i,global_here%wvnxout(i),global_here%wvnyout(i)
                  ENDDO
                  global_here%igpp = global_here%igpp+1+global_here%np
                  global_here%IGWP = global_here%IGWP+1+global_here%NP
               ENDIF
               IF(ABS(global_here%NOUTGW).EQ.2) THEN
                  WRITE(S%FORT73UNIT,REC=global_here%igpp+1) global_here%time_A
                  WRITE(S%FORT73UNIT,REC=global_here%igpp+2) IT
                  WRITE(S%FORT74UNIT,REC=global_here%IGWP+1) global_here%time_A
                  WRITE(S%FORT74UNIT,REC=global_here%IGWP+2) IT
                  global_here%igpp = global_here%igpp + 2
                  global_here%IGWP = global_here%IGWP + 2
                  DO I=1,global_here%NP
                     write(s%fort73unit,rec=global_here%igpp+i) global_here%pr2(i)
                     WRITE(S%FORT74UNIT,REC=global_here%IGWP+2*I-1) global_here%wvnxout(i)
                     WRITE(S%FORT74UNIT,REC=global_here%IGWP+2*I) global_here%wvnyout(i)
                  END DO
                  global_here%igpp = global_here%igpp + global_here%np
                  global_here%IGWP = global_here%IGWP + 2*global_here%NP
               ENDIF
               global_here%NSCOUGW=0
            ENDIF
         ENDIF
         IF(IT.EQ.global_here%NTCYFGW) then
            close(s%fort73unit)
            close(s%fort74unit)
         ENDIF
      endif

      if (it.eq.global_here%nt) then
         open(s%maxele63unit,FILE=S%DIRNAME//'/'//'maxele.63')
         write(s%maxele63unit,*) global_here%np
         do i=1,global_here%np
            write(s%maxele63unit,9633) global_here%etamax(i)
         enddo
 9633    format(f20.8)
         close(s%maxele63unit)
#ifdef OUT_TEC
         if (ModetoNode.eq.1) then
            write(778,*) 'TITLE = "dgswem output"'
            write(778,*) 
     $           'VARIABLES = "global_here%x", "global_here%y", "maxeta"'
            write(778,*) 'ZONE ZONETYPE=FETRIANGLE ',
     $           'NODES=', global_here%np, 
     $           ' ELEMENTS=', global_here%ne, 
     $           ' DATAPACKING=POINT ','SOLUTIONTIME=',global_here%time_a
            do i=1,global_here%np
               write(778,7778) global_here%slam(i)/deg2rad, global_here%sfea(i)/deg2rad, 
     $              global_here%etamax(i)
            enddo
 7778       format(3f20.8)
            do i=1,global_here%ne
               write(778,"(3i12)") global_here%nm(i,1), global_here%nm(i,2), global_here%nm(i,3)
            enddo
            close(778)
         endif
#endif
      endif

!...  
!...  OUTPUT GLOBAL CONCENTRATION DATA IF global_here%NOUTGC<>0 AND THE
!.... TIME STEP FALLS WITHIN THE SPECIFIED WINDOW
!...  
      IF(global_here%NOUTGC.ne.0) THEN
         IF((IT.GT.global_here%NTCYSGC).AND.(IT.LE.global_here%NTCYFGC).OR.FORCE_WRITE) THEN
            global_here%NSCOUGC=global_here%NSCOUGC+1
            IF(global_here%NSCOUGC.EQ.global_here%NSPOOLGC.OR.FORCE_WRITE) THEN
               IF(ABS(global_here%NOUTGC).EQ.1) THEN
                  WRITE(s%fort83unit,2120) global_here%time_A,IT
                  DO I=1,global_here%NP
                     global_here%HH2=global_here%DP(I)+global_here%IFNLFA*global_here%ETA2(I)
                     global_here%C1=global_here%CH1(I)/global_here%HH2
                     IF(global_here%NODECODE(I).EQ.1) WRITE(s%fort83unit,2453) I,global_here%C1
                     IF(global_here%NODECODE(I).EQ.0) WRITE(s%fort83unit,2453) I,-99999.
                  ENDDO
                  global_here%IGCP=global_here%IGCP+1+global_here%NP
               ENDIF
               IF(ABS(global_here%NOUTGC).EQ.2) THEN
                  WRITE(s%fort83unit,REC=global_here%IGEP+1) global_here%time_A
                  WRITE(s%fort83unit,REC=global_here%IGEP+2) IT
                  global_here%IGCP = global_here%IGCP + 2
                  DO I=1,global_here%NP
                     global_here%HH2=global_here%DP(I)+global_here%IFNLFA*global_here%ETA2(I)
                     global_here%C1=global_here%CH1(I)/global_here%HH2
                     IF(global_here%NODECODE(I).EQ.1) WRITE(s%fort83unit,REC=global_here%IGCP+I) global_here%C1
                     IF(global_here%NODECODE(I).EQ.0) WRITE(s%fort83unit,REC=global_here%IGCP+I) -99999.
                  ENDDO
                  global_here%IGCP=global_here%IGCP+global_here%NP
               ENDIF
               global_here%NSCOUGC=0
            ENDIF
         ENDIF
         IF(IT.EQ.global_here%NTCYFGC) CLOSE(s%fort83unit)
      ENDIF

!...  
!...  IF global_here%IHARIND=1 AND THE TIME STEP FALLS WITHIN THE SPECIFIED WINDOW AND
!...  ON THE SPECIFIED INCREMENT, USE MODEL RESULTS TO UPDATE HARMONIC
!...  ANALYSIS MATRIX AND LOAD VECTORS.  NOTE: AN 8 BYTE RECORD SHOULD BE
!...  USED THROUGHOUT THE HARMONIC ANALYSIS SUBROUTINES, EVEN ON 32 BIT
!...  WORKSTATIONS, SINCE IN THAT CASE THE HARMONIC ANALYSIS IS DONE IN
!...  DOUBLE PRECISION.
!...  
#ifdef HARM
      IF(global_here%IHARIND.EQ.1) THEN
         IF((IT.GT.global_here%ITHAS).AND.(IT.LE.global_here%ITHAF)) THEN
            IF(global_here%ICHA.EQ.global_here%NHAINC) global_here%ICHA=0
            global_here%ICHA=global_here%ICHA+1
            IF(global_here%ICHA.EQ.global_here%NHAINC) THEN
!...  
!.....UPDATE THE LHS MATRIX
!...  
               CALL LSQUPDLHS(global_here%timeh,IT)
!...  
!.....IF DESIRED COMPUTE ELEVATION STATION INFORMATION AND UPDATE LOAD VECTOR
!...  
               IF(global_here%NHASE.EQ.1) THEN
                  DO I=1,global_here%NSTAE
!     global_here%EE1=global_here%ETA2(global_here%NM(global_here%NNE(I),1))
!     global_here%EE2=global_here%ETA2(global_here%NM(global_here%NNE(I),2))
!     global_here%EE3=global_here%ETA2(global_here%NM(global_here%NNE(I),3))
!     global_here%ET00(I)=global_here%EE1*global_here%STAIE1(I)+global_here%EE2*global_here%STAIE2(I)+global_here%EE3*global_here%STAIE3(I)
                     global_here%ET00(I) = dg_here%ZE(1,global_here%NNE(I),1)
                     DO K = 2,dg_here%DOFH
                        global_here%ET00(I) = global_here%ET00(I) + dg_here%ZE(K,global_here%NNE(I),1)*dg_here%PHI_STAE(K,I)
                     ENDDO
                  ENDDO
                  CALL LSQUPDES(global_here%ET00,global_here%NSTAE)
               ENDIF
!...  
!.....IF DESIRED COMPUTE VELOCITY STATION INFORMATION AND UPDATE LOAD VECTOR
!...  
               IF(global_here%NHASV.EQ.1) THEN
                  DO I=1,global_here%NSTAV
!     global_here%U11=global_here%UU2(global_here%NM(global_here%NNV(I),1))
!     global_here%U22=global_here%UU2(global_here%NM(global_here%NNV(I),2))
!     global_here%U33=global_here%UU2(global_here%NM(global_here%NNV(I),3))
!     global_here%V11=global_here%VV2(global_here%NM(global_here%NNV(I),1))
!     global_here%V22=global_here%VV2(global_here%NM(global_here%NNV(I),2))
!     global_here%V33=global_here%VV2(global_here%NM(global_here%NNV(I),3))
!     global_here%UU00(I)=global_here%U11*global_here%STAIV1(I)+global_here%U22*global_here%STAIV2(I)+global_here%U33*global_here%STAIV3(I)
!     global_here%VV00(I)=global_here%V11*global_here%STAIV1(I)+global_here%V22*global_here%STAIV2(I)+global_here%V33*global_here%STAIV3(I)
                     DP1     = global_here%DPE(global_here%nnv(i)) !global_here%DP(global_here%NM(global_here%NNV(I),1))
                     !DP2     = global_here%DPE(global_here%nnv(i)) !global_here%DP(global_here%NM(global_here%NNV(I),2))
                     !DP3     = global_here%DPE(global_here%nnv(i)) !global_here%DP(global_here%NM(global_here%NNV(I),3))
                     !DP00    = DP1*global_here%STAIV1(I) + DP2*global_here%STAIV2(I) +DP3*global_here%STAIV3(I)
                     ZE00    = dg_here%ZE(1,global_here%NNV(I),1)
                     global_here%UU00(I) = dg_here%QX(1,global_here%NNV(I),1)
                     global_here%VV00(I) = dg_here%QY(1,global_here%NNV(I),1)
                     DO K = 2,dg_here%DOFH
                        ZE00    = ZE00    + dg_here%ZE(K,global_here%NNV(I),1)*dg_here%PHI_STAV(K,I)
                        global_here%UU00(I) = global_here%UU00(I) + dg_here%QX(K,global_here%NNV(I),1)*dg_here%PHI_STAV(K,I)
                        global_here%VV00(I) = global_here%VV00(I) + dg_here%QY(K,global_here%NNV(I),1)*dg_here%PHI_STAV(K,I)
                     ENDDO
                     DEPTH   = ZE00 + DP1 !DP00
                     FH_NL   = 1.D0/(global_here%NLEQ*DEPTH + global_here%LEQ)
                     global_here%UU00(I) = global_here%UU00(I)*FH_NL
                     global_here%VV00(I) = global_here%VV00(I)*FH_NL
                  ENDDO
                  CALL LSQUPDVS(global_here%UU00,global_here%VV00,global_here%NSTAV)
               ENDIF
!...  
!.....IF DESIRED UPDATE GLOBAL ELEVATION LOAD VECTOR
!...  
               IF(global_here%NHAGE.EQ.1) CALL LSQUPDEG(global_here%ETA2,global_here%NP)
!...  
!.....IF DESIRED UPDATE GLOBAL VELOCITY LOAD VECTOR
!...  
               IF(global_here%NHAGV.EQ.1) CALL LSQUPDVG(global_here%UU2,global_here%VV2,global_here%NP)

            ENDIF
         ENDIF

!...  LINES TO COMPUTE MEANS AND VARIANCES

         if (CHARMV) then
            IF(IT.GT.global_here%ITMV) THEN
               global_here%NTSTEPS=global_here%NTSTEPS+1
               DO I=1,global_here%NP
                  global_here%ELAV(I)=global_here%ELAV(I)+global_here%ETA2(I)
                  global_here%XVELAV(I)=global_here%XVELAV(I)+global_here%UU2(I)
                  global_here%YVELAV(I)=global_here%YVELAV(I)+global_here%VV2(I)
                  global_here%ELVA(I)=global_here%ELVA(I)+global_here%ETA2(I)*global_here%ETA2(I)
                  global_here%XVELVA(I)=global_here%XVELVA(I)+global_here%UU2(I)*global_here%UU2(I)
                  global_here%YVELVA(I)=global_here%YVELVA(I)+global_here%VV2(I)*global_here%VV2(I)
               END DO
            ENDIF
         endif                  !   charmv


      ENDIF
#endif

!...  
!...  WRITE OUT HOT START INFORMATION IF global_here%NHSTAR=1 AND AT CORRECT TIME STEP
!...  NOTE: THE HOT START FILES USE A RECORD LENGTH OF 8 ON BOTH 32 BIT
!.... WORKSTATIONS AND THE 64 BIT CRAY.  THIS IS BECAUSE THE HARMONIC
!.... ANALYSIS IS DONE IN DOUBLE PRECISION (64 BITS) ON WORKSTATIONS.
!...  
#ifdef HOTSTART
      IF(global_here%NHSTAR.EQ.1) THEN
         global_here%ITEST=(IT/global_here%NHSINC)*global_here%NHSINC  
         IF(global_here%ITEST.EQ.IT) THEN
            IF(global_here%IHSFIL.EQ.67) OPEN(67,FILE=S%DIRNAME//'/'//'fort.67',&
            ACCESS='DIRECT',RECL=8)
            IF(global_here%IHSFIL.EQ.68) OPEN(68,FILE=S%DIRNAME//'/'//'fort.68',&
           ACCESS='DIRECT',RECL=8)
            global_here%IHOTSTP=1
            WRITE(global_here%IHSFIL,REC=global_here%IHOTSTP) global_here%IM
            global_here%IHOTSTP=2
            WRITE(global_here%IHSFIL,REC=global_here%IHOTSTP) global_here%TIME_A
            global_here%IHOTSTP=3
            WRITE(global_here%IHSFIL,REC=global_here%IHOTSTP) IT
            DO I=1,global_here%NP
               WRITE(global_here%IHSFIL,REC=global_here%IHOTSTP+1) global_here%ETA1(I)
               WRITE(global_here%IHSFIL,REC=global_here%IHOTSTP+2) global_here%ETA2(I)
               WRITE(global_here%IHSFIL,REC=global_here%IHOTSTP+3) global_here%UU2(I)
               WRITE(global_here%IHSFIL,REC=global_here%IHOTSTP+4) global_here%VV2(I)
               global_here%IHOTSTP = global_here%IHOTSTP + 4
               IF(global_here%IM.EQ.10) THEN
                  WRITE(global_here%IHSFIL,REC=global_here%IHOTSTP+1) global_here%CH1(I)
                  global_here%IHOTSTP=global_here%IHOTSTP+1
               ENDIF
               WRITE(global_here%IHSFIL,REC=global_here%IHOTSTP+1) global_here%NODECODE(I)
               global_here%IHOTSTP=global_here%IHOTSTP+1
            END DO
            WRITE(global_here%IHSFIL,REC=global_here%IHOTSTP+1) global_here%IESTP
            WRITE(global_here%IHSFIL,REC=global_here%IHOTSTP+2) global_here%NSCOUE
            global_here%IHOTSTP=global_here%IHOTSTP+2
            WRITE(global_here%IHSFIL,REC=global_here%IHOTSTP+1) global_here%IVSTP
            WRITE(global_here%IHSFIL,REC=global_here%IHOTSTP+2) global_here%NSCOUV
            global_here%IHOTSTP=global_here%IHOTSTP+2
            WRITE(global_here%IHSFIL,REC=global_here%IHOTSTP+1) global_here%ICSTP
            WRITE(global_here%IHSFIL,REC=global_here%IHOTSTP+2) global_here%NSCOUC
            global_here%IHOTSTP=global_here%IHOTSTP+2
            WRITE(global_here%IHSFIL,REC=global_here%IHOTSTP+1) global_here%IPSTP
            WRITE(global_here%IHSFIL,REC=global_here%IHOTSTP+2) global_here%IWSTP
            WRITE(global_here%IHSFIL,REC=global_here%IHOTSTP+2) global_here%NSCOUM
            global_here%IHOTSTP=global_here%IHOTSTP+3
            WRITE(global_here%IHSFIL,REC=global_here%IHOTSTP+1) global_here%IGEP
            WRITE(global_here%IHSFIL,REC=global_here%IHOTSTP+2) global_here%NSCOUGE
            global_here%IHOTSTP=global_here%IHOTSTP+2
            WRITE(global_here%IHSFIL,REC=global_here%IHOTSTP+1) global_here%IGVP
            WRITE(global_here%IHSFIL,REC=global_here%IHOTSTP+2) global_here%NSCOUGV
            global_here%IHOTSTP=global_here%IHOTSTP+2
            WRITE(global_here%IHSFIL,REC=global_here%IHOTSTP+1) global_here%IGCP
            WRITE(global_here%IHSFIL,REC=global_here%IHOTSTP+2) global_here%NSCOUGC
            global_here%IHOTSTP=global_here%IHOTSTP+2
            WRITE(global_here%IHSFIL,REC=global_here%IHOTSTP+1) global_here%IGPP
            WRITE(global_here%IHSFIL,REC=global_here%IHOTSTP+2) global_here%IGWP
            WRITE(global_here%IHSFIL,REC=global_here%IHOTSTP+3) global_here%NSCOUGW
            global_here%IHOTSTP=global_here%IHOTSTP+3
!...  
!...  IF APPROPRIATE ADD HARMONIC ANALYSIS INFORMATION TO HOT START FILE
!...  
            IF((global_here%IHARIND.EQ.1).AND.(IT.GT.global_here%ITHAS)) THEN
               WRITE(global_here%IHSFIL,REC=global_here%IHOTSTP+1) global_here%ICHA
               global_here%IHOTSTP = global_here%IHOTSTP + 1
               CALL HAHOUT(global_here%NP,global_here%NSTAE,global_here%NSTAV,global_here%NHASE,global_here%NHASV,global_here%NHAGE,global_here%NHAGV,&
              global_here%IHSFIL,global_here%IHOTSTP)
!     
               IF(global_here%NHASE.EQ.1) CALL HAHOUTES(global_here%NSTAE,global_here%IHSFIL,global_here%IHOTSTP)
               IF(global_here%NHASV.EQ.1) CALL HAHOUTVS(global_here%NSTAV,global_here%IHSFIL,global_here%IHOTSTP)
               IF(global_here%NHAGE.EQ.1) CALL HAHOUTEG(global_here%NP,global_here%IHSFIL,global_here%IHOTSTP)
               IF(global_here%NHAGV.EQ.1) CALL HAHOUTVG(global_here%NP,global_here%IHSFIL,global_here%IHOTSTP)
            ENDIF

            if( CHARMV) then
               IF((global_here%IHARIND.EQ.1).AND.(IT.GT.global_here%ITMV)) THEN
                  global_here%IHOTSTP=global_here%IHOTSTP+1
                  WRITE(global_here%IHSFIL,REC=global_here%IHOTSTP) global_here%NTSTEPS
                  IF(global_here%NHAGE.EQ.1) THEN
                     DO I=1,global_here%NP
                        WRITE(global_here%IHSFIL,REC=global_here%IHOTSTP+1) global_here%ELAV(I)
                        WRITE(global_here%IHSFIL,REC=global_here%IHOTSTP+2) global_here%ELVA(I)
                        global_here%IHOTSTP=global_here%IHOTSTP+2
                     END DO
                  ENDIF
                  IF(global_here%NHAGV.EQ.1) THEN
                     DO I=1,global_here%NP
                        WRITE(global_here%IHSFIL,REC=global_here%IHOTSTP+1) global_here%XVELAV(I)
                        WRITE(global_here%IHSFIL,REC=global_here%IHOTSTP+2) global_here%YVELAV(I)
                        WRITE(global_here%IHSFIL,REC=global_here%IHOTSTP+3) global_here%XVELVA(I)
                        WRITE(global_here%IHSFIL,REC=global_here%IHOTSTP+4) global_here%YVELVA(I)
                        global_here%IHOTSTP=global_here%IHOTSTP+4
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
            CLOSE(global_here%IHSFIL)
            IF(global_here%NSCREEN.EQ.1.AND.s%MYPROC.EQ.0) THEN
               WRITE(6,24541) global_here%IHSFIL,IT,global_here%TIME_A
            ENDIF
            WRITE(s%fort16unit,24541) global_here%IHSFIL,IT,global_here%TIME_A
24541       FORMAT(1X,'HOT START OUTPUT WRITTEN TO UNIT ',I2,&
           ' AT TIME STEP = ',I9,' TIME = ',E15.8)
            IF(global_here%IHSFIL.EQ.67) THEN
               global_here%IHSFIL=68
            ELSE
               global_here%IHSFIL=67
            ENDIF
         ENDIF
      ENDIF
#endif      


!...  FIND AND PRINT TO UNIT 6, THE MAXIMUM ELEVATION, THE MAXIMUM VELOCITY
!.... AND THE NODE NUMBERS AT WHICH THEY OCCUR ON MYPROC=0
!.... IF global_here%ELMAX EXCEEDS THRESHOLD, PRINT INFORMATION ON ALL PROCESSORS WHERE
!.... THIS OCCURS 
!...  
      IF(global_here%NSCREEN.EQ.1) THEN
         global_here%ELMAX=0.0
         global_here%VELMAX=0.0
         global_here%KEMAX = 0
         global_here%KVMAX = 0
!     
         imaxze=0.0
         elmaxe=0.0
         qmaxe=0.0
         imaxq = 0.0
         DO I=1,global_here%NP
            IF((global_here%NODECODE(I).EQ.1).AND.(ABS(global_here%ETA2(I)).GT.global_here%ELMAX))THEN
               global_here%ELMAX=ABS(global_here%ETA2(I))
               global_here%KEMAX=I
            ENDIF
            global_here%VELABS=global_here%UU2(I)*global_here%UU2(I)+global_here%VV2(I)*global_here%VV2(I)
            IF (global_here%VELABS.GT.global_here%VELMAX) THEN
               global_here%VELMAX=global_here%VELABS
               global_here%KVMAX=I
            ENDIF
         END DO
!     nd
         do i=1,global_here%ne
            if (dg_here%wdflg(i).ne.0) then
               if (abs(dg_here%ze(1,i,1)).gt.elmaxe) then
                  imaxze=i
                  elmaxe=abs(dg_here%ze(1,i,1))
               endif
               if (sqrt(dg_here%qx(1,i,1)**2+dg_here%qy(1,i,1)**2).gt.qmaxe) then
                  qmaxe=sqrt(dg_here%qx(1,i,1)**2+dg_here%qy(1,i,1)**2)
                  imaxq=i
               endif
               
            endif
         enddo
         
!         global_here%VELMAX=global_here%VELMAX**0.5d0
         
#ifdef CMPI
!     sb
!     IF(MYPROC.EQ.0.AND.global_here%ELMAX.LT.200.0.AND.global_here%KEMAX.GT.0) THEN
         IF( (s%MYPROC.EQ.0).AND.(global_here%ELMAX.LT.200.0).AND.&
        (MOD(IT,global_here%NSCREEN_INC).EQ.0) ) THEN
            if (global_here%kemax.eq.0) then 
               global_here%kemax = 1 
            endif
       print*,'________________________________',&
       '__________________________________________'
            print*,''

         WRITE(6,1991) IT,global_here%NUMITR,global_here%ETA2(global_here%KEMAX),global_here%KEMAX,global_here%VELMAX,global_here%KVMAX,&
           s%MYPROC
 1991   FORMAT(1X,'TIME STEP =',I8,5X,'ITERATIONS =',I5,&
           /,2X,'global_here%ELMAX = ', E11.4,' AT NODE',I7,&
           2X,'SPEEDMAX = ',E11.4,' AT NODE',I7,&
           2X,'ON MYPRnOC = ',I4)
   !write(6,*) 'elmaxe, imaxze, qmaxe, imaxq ',elmaxe,imaxze,
                                !$           qmaxe,imaxq


!$$$  if (global_here%chem_flag.eq.1) then
!$$$  print*,'__________________________________________________________________________'
!$$$  print*,'|                       |                                                |'
!$$$  print*,'| Chemistry turned "ON" |',' Maximum mass action =',-minval(global_here%MassMax(:)),' |'
!$$$  print*,'|_______________________|________________________________________________|'
!$$$  endif


            if (dg_here%padapt.eq.1) then

               do k=dg_here%pl,dg_here%ph
                  Minp(k) = 0
                  do i = 1,global_here%ne
                     if (global_here%pdg_el(i).eq.k) then
                        Minp(k) = Minp(k) + 1
                     endif
                  enddo
               enddo

             print*,'_________________________________________________'
             print*,'|                      |                        |'
             print*,'|  Polynomial order    |     # of Elements      |'
             print*,'|______________________|________________________|'
              do k=dg_here%pl,dg_here%ph
             print*,'|',k,'         |',Minp(k),'           |'
             print*,'|______________________|________________________|'
              enddo
               print*,''
               print*,'With',global_here%ne,'total elements'               
            endif

!$$$  if (dg_here%slopeflag.eq.4.or.dg_here%slopeflag.ge.6) then
!$$$  
!$$$  print*,''
!$$$  print*,'p-adaptive slope limiting "ACTIVE."'
!$$$  print*,''
!$$$  print*,dg_here%lim_count_roll,'elements of a total of'
!$$$  print*,global_here%ne,'elements were VERTEX LIMITED'
!$$$  print*,''
!$$$  
!$$$  endif

            print*,'___________________________________',&
          '_______________________________________'

         ENDIF
!     sb
!     IF(global_here%ELMAX.GT.200.0.AND.global_here%KEMAX.GT.0) THEN
         IF(global_here%ELMAX.GT.200.0) THEN
            if (global_here%kemax.eq.0) then 
               global_here%kemax = 1 
            endif
            WRITE(6,1993) IT,global_here%NUMITR,global_here%ETA2(global_here%KEMAX),global_here%KEMAX,global_here%VELMAX,global_here%KVMAX,&
           s%MYPROC
            WRITE(s%fort16unit,1993) IT,global_here%NUMITR,global_here%ETA2(global_here%KEMAX),global_here%KEMAX,global_here%VELMAX,global_here%KVMAX,&
           s%MYPROC
 1993       FORMAT(1X,'TIME STEP =',I8,6X,'ITERATIONS =',I5,&
      /,2X,'global_here%ELMAX = ', E11.4,' AT NODE',I7,&
      2X,'SPEEDMAX = ',E11.4,' AT NODE',I7,&
      2X,'ON MYPROC = ',I4,' !!! WARNING - HIGH ELEVATION !!!')
            STOP
         ENDIF


#else
!     sb
!     IF(global_here%ELMAX.LT.200.0.AND.global_here%KEMAX.GT.0) THEN

         IF(global_here%ELMAX.LT.200.0.and.(MOD(IT,global_here%NSCREEN_INC).EQ.0).AND.&
            (s%MYPROC.EQ.0)) THEN
            if (global_here%kemax.eq.0) then 
               global_here%kemax = 1 
            endif
            print*,'______________________________________',&
            '____________________________________'
            print*,''
            WRITE(6,1992) IT,global_here%NUMITR,global_here%ETA2(global_here%KEMAX),global_here%KEMAX,global_here%VELMAX,global_here%KVMAX
 1992       FORMAT(1X,'TIME STEP =',I8,5X,'ITERATIONS =',I5,&
           /,2X,'global_here%ELMAX = ', E11.4,' AT NODE',I7,&
           2X,'SPEEDMAX = ',E11.4,' AT NODE',I7)
                                !write(6,*) 'elmaxe, imaxze, qmaxe, imaxq ',elmaxe,imaxze,
                                !$           qmaxe,imaxq

!$$$  if (global_here%chem_flag.eq.1) then
!$$$  print*,'__________________________________________________________________________'
!$$$  print*,'|                       |                                                |'
!$$$  print*,'| Chemistry turned "ON" |',' Maximum mass action =',-minval(global_here%MassMax(:)),' |'
!$$$  print*,'|_______________________|________________________________________________|'
!$$$  endif


            if (dg_here%padapt.eq.1) then
               
               do k=dg_here%pl,dg_here%ph
                  Minp(k) = 0
                  do i = 1,global_here%ne
                     if (global_here%pdg_el(i).eq.k) then
                        Minp(k) = Minp(k) + 1
                     endif
                  enddo
               enddo

             print*,'_________________________________________________'
             print*,'|                      |                        |'
             print*,'|  Polynomial order    |     # of Elements      |'
             print*,'|______________________|________________________|'
               do k=dg_here%pl,dg_here%ph
             print*,'|',k,'         |',Minp(k),'           |'
             print*,'|______________________|________________________|'
               enddo
               print*,''
               print*,'With',global_here%ne,'total elements'               
            endif

!$$$  if (dg_here%slopeflag.eq.4.or.dg_here%slopeflag.ge.6) then
!$$$  
!$$$  print*,''
!$$$  print*,'p-adaptive slope limiting "ACTIVE."'
!$$$  print*,''
!$$$  print*,dg_here%lim_count_roll,'elements of a total of'
!$$$  print*,global_here%ne,'elements were VERTEX LIMITED'
!$$$  print*,''
!$$$  
!$$$  endif

            print*,'___________________________________',&
          '_______________________________________'

         ENDIF
!     sb
!     IF(global_here%ELMAX.GT.200.0.AND.global_here%KEMAX.GT.0) THEN
         IF(global_here%ELMAX.GT.200.0) THEN
          WRITE(6,1994) IT,global_here%NUMITR,global_here%ETA2(global_here%KEMAX),global_here%KEMAX,global_here%VELMAX,global_here%KVMAX
          WRITE(s%fort16unit,1994) IT,global_here%NUMITR,global_here%ETA2(global_here%KEMAX),global_here%KEMAX,global_here%VELMAX,global_here%KVMAX
 1994       FORMAT(1X,'TIME STEP =',I8,6X,'ITERATIONS =',I5,&
           /,2X,'global_here%ELMAX = ', E11.4,' AT NODE',I7,&
           2X,'SPEEDMAX = ',E11.4,' AT NODE',I7,&
           2X,' !!! WARNING - HIGH ELEVATION !!!')
            STOP
         ENDIF
#endif
      ENDIF
      
!.....If applicable write out a DG hot start

      IF ((dg_here%DGHOT.EQ.1).AND.(MOD(IT,dg_here%DGHOTSPOOL).EQ.0)) THEN
         
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

         WRITE(263,*) dg_here%PH
         WRITE(264,*) dg_here%PH, dg_here%PH
         WRITE(214,*) IT

         DO J = 1,s%MNE
            DO K = 1,dg_here%DOFH
               WRITE(263,*) dg_here%ZE(K,J,1)
               WRITE(264,*) dg_here%QX(K,J,1), dg_here%QY(K,J,1)
               WRITE(214,*) dg_here%HB(K,J,1), dg_here%WDFLG(J)
#ifdef TRACE
               WRITE(288,*) dg_here%iota(K,J,1)
#endif
#ifdef CHEM
               WRITE(289,*) dg_here%iota(K,J,1),dg_here%iota2(K,J,1)
#endif
#ifdef DYNP
               WRITE(291,*) dg_here%dynP(K,J,1)
#endif
#ifdef SED_LAY
               do ll=1,s%layers
                  WRITE(290,*) dg_here%bed(K,J,1,ll)
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

!     IF (IT.EQ.global_here%NT) THEN
!     dg_here%H_TRI = SQRT((global_here%X(1)-global_here%X(2))**2 + (global_here%Y(1)-global_here%Y(2))**2)
!     WRITE(S%FORT88UNIT,*) dg_here%H_TRI,P,0
!     WRITE(S%FORT89UNIT,*) dg_here%H_TRI,P,0
!     IF (MYPROC == 0) THEN
!     PRINT*,'T FINAL =',global_here%NT*global_here%DTDP
!     PRINT*,'global_here%DT =',global_here%DTDP
!     PRINT*,'global_here%NT =',global_here%NT
!     ENDIF

!     DO J=1,s%MNE
!     global_here%ZE_DG(1) = 0.D0
!     global_here%QX_DG(1) = 0.D0
!     global_here%QY_DG(1) = 0.D0
!     
!     DO K=1,dg_here%DOF
!     global_here%ZE_DG(1) = global_here%ZE_DG(1) + dg_here%PHI_CENTER(K)*dg_here%ZE(K,J,1)
!     global_here%QX_DG(1) = global_here%QX_DG(1) + dg_here%PHI_CENTER(K)*dg_here%QX(K,J,1)
!     global_here%QY_DG(1) = global_here%QY_DG(1) + dg_here%PHI_CENTER(K)*dg_here%QY(K,J,1)
!     
!     WRITE(S%FORT88UNIT,*) dg_here%ZE(K,J,1),dg_here%QX(K,J,1),dg_here%QY(K,J,1)
!     
!     ENDDO
!     
!     WRITE(S%FORT89UNIT,*) global_here%ZE_DG(1),global_here%QX_DG(1),global_here%QY_DG(1)
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

      SUBROUTINE WRITE_DG_IC(dg_here)

      USE GLOBAL
      USE DG

      type (dg_type) :: dg_here

!.....DG.63.IC = Record of the initial surface elevation
!      OPEN(632,FILE=S%DIRNAME//'/'//'DG.63.IC')
!      WRITE(632,3220) global_here%RUNDES, global_here%RUNID, global_here%AGRID
!      WRITE(632,3645) 1, dg_here%DOF, global_here%DTDP*global_here%NSPOOLGE, global_here%NSPOOLGE, 1

!      WRITE(632,2120) 0.0, 1
!      DO J = 1,global_here%NE
!        DO K = 1,dg_here%DOF
!          WRITE(632,2453) J, dg_here%ZE(K,J,1)
!        ENDDO
!      ENDDO
      
!      CLOSE(632)

!.....DG.64.IC = Record of the initial discharge
!      OPEN(642,FILE=S%DIRNAME//'/'//'DG.64.IC')
!      WRITE(642,3220) global_here%RUNDES, global_here%RUNID, global_here%AGRID
!      WRITE(642,3645) 1, dg_here%DOF, global_here%DTDP*global_here%NSPOOLGV, global_here%NSPOOLGV, 2
!
!      WRITE(642,2120) 0.0, 1
!      DO J = 1,global_here%NE
!        DO K = 1,dg_here%DOF
!          WRITE(642,2454) J, dg_here%QX(K,J,1), dg_here%QY(K,J,1)
!        ENDDO
!      ENDDO
!
!      CLOSE(642)

!.....DG.65.IC = Record of the initial wet/dry flags
!      IF (global_here%NOLIFA.GE.2) THEN
!        OPEN(652,FILE=S%DIRNAME//'/'//'DG.65.IC')
!        WRITE(652,3220) global_here%RUNDES, global_here%RUNID, global_here%AGRID
!        WRITE(652,3645) 1, dg_here%DOF, global_here%DTDP*global_here%NSPOOLGE, global_here%NSPOOLGE, 1
!
!        WRITE(652,2120) 0.0, 1
!       DO J = 1,global_here%NE
!          WRITE(652,2455) J,dg_here%WDFLG(J)
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
