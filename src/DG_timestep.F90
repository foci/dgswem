!***********************************************************************
!     
!     SUBROUTINE DG_TIMESTEP(IT)
!     
!     This is the main subroutine for the DG hydro
!     
!     Written by Shintaro Bunya (01-01-2007)
!     
!***********************************************************************

      SUBROUTINE DG_TIMESTEP(s,dg_here,global_here,nodalattr_here,IT)

!.....Use appropriate modules
      
      USE SIZES
      USE GLOBAL
      USE NodalAttributes

#ifdef CMPI
      USE MESSENGER_ELEM
#endif

      USE DG

      IMPLICIT NONE

      type (sizes_type) :: s
      type (dg_type) :: dg_here
      type (global_type) :: global_here
      type (nodalattr_type) :: nodalattr_here

!.....Declare local variables

      INTEGER IT, IRK
      
!.....Hydrodynamics
#ifdef HPX
      CALL DG_HYDRO_TIMESTEP(s,dg_here,global_here,nodalattr_here,IT,IRK)
#else


      DO IRK = 1,dg_here%NRK
        CALL DG_HYDRO_TIMESTEP(s,dg_here,global_here,nodalattr_here,IT,IRK)
      ENDDO
      
      CALL DG_TIMESTEP_ADVANCE(s,dg_here,global_here,nodalattr_here,IT)
      
#endif

!.....Write out results

! #ifdef CMPI
! #ifdef BLKOUT
!       CALL MSG_BLOCKSYNC_START()
! #endif
! #endif
!       CALL WRITE_RESULTS(s,dg_here,global_here,IT,.FALSE.)
! #ifdef CMPI
! #ifdef BLKOUT
!       CALL MSG_BLOCKSYNC_FINISH()
! #endif
! #endif
! 
!       CALL SCRUTINIZE_SOLUTION(s,dg_here,global_here,IT)

      RETURN
      END SUBROUTINE


!***********************************************************************
!     
!     SUBROUTINE SCRUTINIZE_SOLUTION(IT)
!     
!     Detect a negative water depth and a NaN in solution arrays and 
!     terminate if either of them is found.
!     
!     Written by Shintaro Bunya (01-01-2007)
!     
!***********************************************************************

      SUBROUTINE SCRUTINIZE_SOLUTION(s,dg_here,global_here,IT)

!.....Use appropriate modules
      
#ifdef CMPI
      USE MESSENGER_ELEM, ONLY: MESSAGE_FINI, ErrorElevSum
#else
      USE SIZES
#endif
      USE GLOBAL
      USE DG

      IMPLICIT NONE
      
      type (sizes_type) :: s
      type (dg_type) :: dg_here
      type (global_type) :: global_here

!.....Declare local variables

      INTEGER IT,J,K,l
      LOGICAL Detected
      REAL(SZ) DPAVG

      INTEGER ErrorElevExceeded

!.....Detect negative mass
      ErrorElevExceeded=0
      DO J = 1,global_here%NE
         DPAVG =  1.D0/dg_here%ncheck(global_here%pdg_el(j)) * (sum(dg_here%DP_NODE(:,J,global_here%pdg_el(j)))) !(global_here%DP(global_here%NM(J,1))+global_here%DP(global_here%NM(J,2))+global_here%DP(global_here%NM(J,3)))/3.D0
         Detected = .FALSE.
         IF((dg_here%ZE(1,J,1)+DPAVG).LE.0.D0) THEN
            PRINT *, ''
#ifdef CMPI
            PRINT *, 'IN SUBDOMAIN ',s%MYPROC
#endif
            PRINT *, '  dg_here%ZE(',1,',',J,',1) = ',dg_here%ZE(1,J,1)
            PRINT *, '  HE(',1,',',J,',1) = ',dg_here%ZE(1,J,1)+DPAVG
            PRINT *, '  global_here%DP = ',global_here%DP(global_here%NM(J,1)),global_here%DP(global_here%NM(J,2)),global_here%DP(global_here%NM(J,3))
            PRINT *, '  DPAVG = ',DPAVG
            PRINT *, ' global_here%x,  global_here%y ',global_here%SLAM(global_here%NM(J,1))/deg2rad,global_here%SFEA(global_here%NM(J,1))/DEG2RAD
            Detected = .TRUE.
            ErrorElevExceeded=1    
         ENDIF
      ENDDO
#ifdef CMPI
      CALL ErrorElevSum(ErrorElevExceeded)
      IF(ErrorElevExceeded.NE.0) THEN
         if (s%myproc.eq.0) then
            PRINT *, ''
            PRINT *,'  PROGRAM WILL BE TERMINATED in DG_timestep'
            PRINT *, ''
         endif

         call message_fini (s)
         STOP
      endif
#else
!$$$      IF(ErrorElevExceeded.global_here%NE.0) THEN
!$$$         PRINT *, ''
!$$$         PRINT *,'  PROGRAM WILL BE TERMINATED in DG_timestep'
!$$$         PRINT *, ''
!$$$         stop
!$$$      ENDIF
#endif


!$$$
!$$$         IF(Detected) THEN
!$$$            PRINT *, ''
!$$$            PRINT *,'  WRITING OUT THE LAST-MOMENT RESULTS...'
!$$$            PRINT *, ''
!$$$            CALL WRITE_RESULTS(IT,.TRUE.)
!$$$            PRINT *, ''
!$$$            PRINT *,'  PROGRAM WILL BE TERMINATED.'
!$$$            PRINT *, ''
!$$$#ifdef CMPI
!$$$            call message_fini ()
!$$$            STOP
!$$$#else
!$$$            stop
!$$$#endif
!$$$         ENDIF
!$$$      ENDDO

!.....Detect NaN
      
      DO J = 1,global_here%NE
         DO K = 1,dg_here%DOF

            Detected = .FALSE.
            IF(dg_here%ZE(K,J,1).NE.dg_here%ZE(K,J,1)) THEN
               PRINT *, ''
#ifdef CMPI
               PRINT *, 'IN SUBDOMAIN ',s%MYPROC
#endif
               PRINT *, '  dg_here%ZE(',K,',',J,',1) IS ',dg_here%ZE(K,J,1)
               PRINT *, ' global_here%x,  global_here%y ',global_here%SLAM(global_here%NM(J,1))/DEG2RAD,global_here%SFEA(global_here%NM(J,1))/DEG2RAD
               Detected = .TRUE.
               ErrorElevExceeded = 1
            ENDIF
            IF(dg_here%QX(K,J,1).NE.dg_here%QX(K,J,1)) THEN
#ifdef CMPI
               PRINT *, 'IN SUBDOMAIN ',s%MYPROC
#endif
               PRINT *, ''
               PRINT *, '  dg_here%QX(',K,',',J,',1) IS ',dg_here%QX(K,J,1)
               PRINT *, ' global_here%x,  global_here%y ',global_here%SLAM(global_here%NM(J,1))/DEG2RAD,global_here%SFEA(global_here%NM(J,1))/DEG2RAD
               Detected = .TRUE.
               ErrorElevExceeded = 1
            ENDIF
            IF(dg_here%QY(K,J,1).NE.dg_here%QY(K,J,1)) THEN
#ifdef CMPI
               PRINT *, 'IN SUBDOMAIN ',s%MYPROC
#endif
               PRINT *, ''
               PRINT *, '  dg_here%QY(',K,',',J,',1) IS ',dg_here%QY(K,J,1)
               PRINT *, ' global_here%x,  global_here%y ',global_here%SLAM(global_here%NM(J,1))/deg2rad,global_here%SFEA(global_here%NM(J,1))/deg2rad
               Detected = .TRUE.
               ErrorElevExceeded = 1
            ENDIF

#ifdef TRACE
            IF(dg_here%iota(K,J,1).NE.dg_here%iota(K,J,1)) THEN
               PRINT *, ''

#ifdef CMPI
               PRINT *, 'IN SUBDOMAIN ',s%MYPROC
#endif
               PRINT *, '  dg_here%iota(',K,',',J,',1) IS ',dg_here%iota(K,J,1)
               PRINT *, ' global_here%x,  global_here%y ',global_here%X(global_here%NM(J,1)),global_here%Y(global_here%NM(J,1))
               Detected = .TRUE.
               ErrorElevExceeded = 1
            ENDIF
#endif
            
#ifdef CHEM
            IF(dg_here%iota(K,J,1).NE.dg_here%iota(K,J,1)) THEN
               PRINT *, ''
#ifdef CMPI
               PRINT *, 'IN SUBDOMAIN ',s%MYPROC
#endif
               PRINT *, '  dg_here%iota(',K,',',J,',1) IS ',dg_here%iota(K,J,1)
               PRINT *, ' global_here%x,  global_here%y ',global_here%X(global_here%NM(J,1)),global_here%Y(global_here%NM(J,1))
               Detected = .TRUE.
            ENDIF
            IF(dg_here%iota2(K,J,1).NE.dg_here%iota2(K,J,1)) THEN
               PRINT *, ''
#ifdef CMPI
               PRINT *, 'IN SUBDOMAIN ',s%MYPROC
#endif
               PRINT *, '  dg_here%iota2(',K,',',J,',1) IS ',dg_here%iota2(K,J,1)
               PRINT *, ' global_here%x,  global_here%y ',global_here%X(global_here%NM(J,1)),global_here%Y(global_here%NM(J,1))
               Detected = .TRUE.
               ErrorElevExceeded = 1
            ENDIF
#endif

#ifdef DYNP
            IF(dg_here%dynP(K,J,1).NE.dg_here%dynP(K,J,1)) THEN
               PRINT *, ''

#ifdef CMPI
               PRINT *, 'IN SUBDOMAIN ',MYPROC
#endif
               PRINT *, '  dg_here%dynP(',K,',',J,',1) IS ',dg_here%dynP(K,J,1)
               PRINT *, ' global_here%x,  global_here%y ',global_here%X(global_here%NM(J,1)),global_here%Y(global_here%NM(J,1))
               Detected = .TRUE.
               ErrorElevExceeded = 1
            ENDIF
#endif
#ifdef SED_LAY
            
            do l = 1,s%layers
               IF(dg_here%bed(K,J,1,l).NE.dg_here%bed(K,J,1,l)) THEN
                  PRINT *, ''
                  
#ifdef CMPI
                  PRINT *, 'IN SUBDOMAIN ',MYPROC
#endif
                  PRINT *, '  dg_here%bed(',K,',',J,',1,',l,') IS ', dg_here%bed(K,J,1,l)
                  PRINT *, ' global_here%x,  global_here%y ',global_here%X(global_here%NM(J,1)),global_here%Y(global_here%NM(J,1))
                  Detected = .TRUE.
                  ErrorElevExceeded = 1
               ENDIF
            enddo
#endif

!$$$            IF(Detected) THEN
!$$$               PRINT *, ''
!$$$               PRINT *,'  WRITING OUT THE LAST-MOMENT RESULTS...'
!$$$               PRINT *, ''
!$$$               CALL WRITE_RESULTS(IT,.TRUE.)
!$$$               PRINT *, ''
!$$$               PRINT *,'  PROGRAM WILL BE TERMINATED.'
!$$$               PRINT *, ''
!$$$#ifdef CMPI
!$$$               call message_fini()
!$$$               STOP
!$$$#else
!$$$               stop
!$$$#endif
!$$$            ENDIF
         ENDDO
      ENDDO
#ifdef CMPI
      CALL ErrorElevSum(ErrorElevExceeded)
!      write(*,*) myproc,ErrorElevExceeded
      IF (ErrorElevExceeded.NE.0) THEN
         CALL MESSAGE_FINI(s)
         stop
      ENDIF
#else
!      write(*,*)ErrorElevExceeded
!$$$      IF (ErrorElevExceeded.global_here%NE.0) THEN
!$$$         STOP
!$$$      ENDIF
#endif

      RETURN
      END SUBROUTINE
