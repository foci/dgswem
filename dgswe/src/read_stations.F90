      SUBROUTINE read_stations()

      USE globals, ONLY: nsta,xysta
      USE read_dginp, ONLY: sta_opt,stations_file
      USE messenger2, ONLY: myrank
      USE quit, ONLY: abort

      IMPLICIT NONE
      
      INTEGER :: sta
      LOGICAL :: file_exists
      

      
      file_exists = .FALSE.
      INQUIRE(FILE=trim(stations_file),EXIST=file_exists)
      
      IF(file_exists) THEN
        
        IF (myrank == 0 ) PRINT "(A)", "Reading in stations file"
      
        OPEN(UNIT=15,FILE=trim(stations_file))
        READ(15,*) nsta        
        ALLOCATE(xysta(2,nsta))
        DO sta = 1,nsta
          READ(15,*) xysta(1,sta), xysta(2,sta)
        ENDDO
        CLOSE(15)      

      ELSE
     
        IF (sta_opt > 0) THEN
          IF (myrank == 0) PRINT*, "stations file required for sta_opt > 0"
          CALL abort()
        ENDIF

      ENDIF
        
        
        
      RETURN
      END SUBROUTINE read_stations