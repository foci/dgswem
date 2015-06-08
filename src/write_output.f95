C-----------------------------------------------------------------------
C     S U B R O U T I N E   W R I T E  O U T  A R R A Y  
C-----------------------------------------------------------------------
C     jgf48.03 This subroutine was created to write out a column 
C     vector (i.e., nodal data such as water surface elevation or 
C     pressure) to a file.
C-----------------------------------------------------------------------
Casey 090302: Added the filename as the last argument.
      SUBROUTINE writeOutArray(lun, time, it, descript, pack_cmd,
     &                         unpack_cmd, filepos, fn)

#ifdef SWAN

      USE SIZES
      USE GLOBAL
      USE GLOBAL_IO, ONLY : collectFullDomainArray
      IMPLICIT NONE
C     args
      INTEGER, intent(in) :: lun ! logical unit number of file to write to
      REAL(8), intent(in) :: time ! seconds since cold start
      INTEGER, intent(in) :: it   ! number of time steps since cold start
      type(OutputDataDescript_t), intent(in) :: descript !describes output data
      EXTERNAL :: pack_cmd   ! subroutine used to pack data on subdomain
      EXTERNAL :: unpack_cmd ! subroutine used to unpack data on proc 0
      INTEGER, intent(inout) :: filepos  ! current position in the output file
C     local vars
Casey 090302: Increase the length of the filename from 7 to 15.
      CHARACTER(15) :: fn     ! outfile name (valid for lun between 10 and 99)
      INTEGER :: I           ! loop counter

C     initialize output file name
Casey 090302: The filename is now passed as an input argument.
C     fn(1:5) = 'fort.'
C     WRITE(fn(6:7),2) lun

C     collect up the data from subdomains if running in parallel
      IF ((MNPROC.gt.1).and.(.not.WRITE_LOCAL_FILES)) THEN
C         write(16,*) 'About to collectFullDomainArray' !jgfdebug48.03
         CALL collectFullDomainArray(descript, pack_cmd, unpack_cmd)
      ENDIF

C     write data according to format specifier from fort.15 (e.g., NOUTE)
      SELECT CASE (ABS(descript % specifier))

      CASE(1) !ascii text

C         write(16,*) 'About to open globalio text file.' !jgfdebug48.03
         IF ( (MNPROC.gt.1).and.(MyProc.eq.0)
     &        .and.(.not.WRITE_LOCAL_FILES)) THEN
            OPEN(lun,FILE=TRIM(GLOBALDIR)//'/'//fn,
     &           ACCESS='SEQUENTIAL',POSITION='APPEND')
            WRITE(lun,2120) time,IT
            IF (descript % num_items_per_record .eq. 1) THEN
               DO I=1, descript % num_fd_records
C                  WRITE(16,*) I ! jgfdebug48.03
                  WRITE(lun,2453) I, descript % array_g(I)
               ENDDO
            ENDIF
            IF (descript % num_items_per_record .eq. 2) THEN
               DO I=1, descript % num_fd_records
C                  WRITE(16,*) I ! jgfdebug48.03
                  WRITE(lun,2454) I, descript % array_g(I),
     &                               descript % array2_g(I)
               ENDDO
            ENDIF
C            WRITE(16,*) 'About to close globalio text file.'!jgfdebug48.03
            CLOSE(lun)
         ENDIF

         IF ((MNPROC.eq.1).or.(WRITE_LOCAL_FILES)) THEN
            OPEN(lun,FILE=TRIM(LOCALDIR)//'/'//fn,
     &           ACCESS='SEQUENTIAL',POSITION='APPEND')
            WRITE(lun,2120) time,IT
            IF (descript % num_items_per_record .eq. 1) THEN
!......TCM - v48.4618 -- Fixed dry node output for serial run case
               IF ((trim(descript % field_name) .eq. 'Elev').and.
     &             (descript % ConsiderWetDry .EQV. .TRUE.)) THEN
                  DO I=1, descript % num_records_this
                     if(NODECODE(I).EQ.1) THEN
                        WRITE(lun,2453) I, descript % array(I)
                     ELSE
                        WRITE(lun,2453) I, descript % alternate_value  !-99999.0 for dry nodes
                     ENDIF
                  END DO
               ELSE
                  DO I=1, descript % num_records_this
                     WRITE(lun,2453) I, descript % array(I)
                  END DO
               ENDIF
            ENDIF
            IF (descript % num_items_per_record .eq. 2) THEN
               DO I=1, descript % num_records_this
                  WRITE(lun,2454) I, descript % array(I),
     &                               descript % array2(I)
               END DO
            ENDIF

         ENDIF
         filepos = filepos+1+descript % num_records_this

      CASE(2) !binary (nonportable)

         IF ( (MNPROC.gt.1).and.(MyProc.eq.0)
     &         .and.(.not.WRITE_LOCAL_FILES)) THEN
            OPEN(lun,FILE=TRIM(GLOBALDIR)//'/'//fn,
     &           ACCESS='DIRECT',RECL=NBYTE)
            WRITE(lun,REC=filepos+1) time
            WRITE(lun,REC=filepos+2) IT
            filepos = filepos + 2
            IF ( descript % num_items_per_record .eq. 1 ) THEN
               DO I=1, descript % num_fd_records
                  WRITE(lun,REC=filepos+I) descript % array_g(I)
               END DO
            ENDIF
            IF ( descript % num_items_per_record .eq. 2 ) THEN
               DO I=1, descript % num_fd_records
                  WRITE(lun,REC=filepos+2*I-1) descript % array_g(I)
                  WRITE(lun,REC=filepos+2*I)   descript % array2_g(I)
               END DO
            ENDIF
            CLOSE(lun)
         ENDIF

         IF ((MNPROC.eq.1).or.(WRITE_LOCAL_FILES)) THEN
            OPEN(lun,FILE=TRIM(LOCALDIR)//'/'//fn,
     &           ACCESS='DIRECT',RECL=NBYTE)
            WRITE(lun,REC=filepos+1) time
            WRITE(lun,REC=filepos+2) IT
            filepos = filepos + 2
            IF ( descript % num_items_per_record .eq. 1 ) THEN
               IF ((trim(descript % field_name) .eq. 'Elev').and.
     &             (descript % ConsiderWetDry .EQV. .TRUE.)) THEN
                  DO I=1, descript % num_records_this
                     if(NODECODE(I).EQ.1) THEN
                        WRITE(lun,REC=filepos+I) descript % array(I)
                     ELSE
                        WRITE(lun,REC=filepos+I)
     &                    descript % alternate_value !-99999.0 for dry nodes
                     ENDIF
                  END DO
               ELSE
                  DO I=1, descript % num_records_this
                     WRITE(lun,REC=filepos+I) descript % array(I)
                  END DO
               ENDIF
            ENDIF
            IF ( descript % num_items_per_record .eq. 2 ) THEN
               DO I=1, descript % num_records_this
                  !tcmv48.4618 -- changed from array_g to array
                  WRITE(lun,REC=filepos+2*I-1) descript % array(I)
                  !tcmv48.4618 -- changed from array2_g to array2
                  WRITE(lun,REC=filepos+2*I)   descript % array2(I)
               END DO
            ENDIF
            CLOSE(lun)
         ENDIF
         filepos = filepos + descript % num_records_this

      CASE(3) !netcdf (portable)
#ifdef NETCDF
!       IF(MYPROC.EQ.0) PRINT *,"W: BEFORE WRITING OUT ARRAY "
         CALL writeOutArrayNetCDF(lun, time, it, descript)
!       IF(MYPROC.EQ.0) PRINT *,"W: AFTER WRITING OUT ARRAY "
#else
         WRITE(ScreenUnit,*) 'ERROR: NetCDF is not available.'
         WRITE(16,*) 'ERROR: NetCDF is not available.'
#endif

      CASE DEFAULT
         WRITE(ScreenUnit,*) 'ERROR: Invalid output specifier.'
      END SELECT

 2    FORMAT(I2)
 2120 FORMAT(2X,E20.10,5X,I10)
 2453 FORMAT(2x, i8, 2x, E20.10, E20.10, E20.10, E20.10)
 2454 FORMAT(2X,I8,2(2X,E15.8))

#endif

C-----------------------------------------------------------------------
       END SUBROUTINE writeOutArray
C-----------------------------------------------------------------------


