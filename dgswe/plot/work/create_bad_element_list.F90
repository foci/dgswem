      PROGRAM create_bad_elem_list

      IMPLICIT NONE

      
      INTEGER :: read_stat
      INTEGER :: el
      INTEGER :: i
      INTEGER :: bad_list(10000)
      INTEGER :: nbad
      INTEGER :: flag
      CHARACTER(100) :: line
      CHARACTER(7) :: el_char

      
      OPEN(UNIT=12,FILE='bad_elems.txt')
      
      nbad = 0
      
      DO 
      
        READ(12,"(A)",IOSTAT=read_stat) line   
        PRINT "(A)", line
        
        IF (read_stat < 0) THEN
          EXIT
        ENDIF
        
        IF (INDEX(line,":",back=.true.) > 35) THEN
!           PRINT*, line
          el_char = line(39:45)
          READ(el_char,*) el
          
          flag = 0
          DO i = 1,nbad
            IF (bad_list(i) == el) THEN
              flag = 1
              EXIT
            ENDIF
          ENDDO
          
          IF (flag == 0) THEN
            nbad = nbad + 1
            bad_list(nbad) = el
          ENDIF
          
        ENDIF        
      
      ENDDO
      
      OPEN(UNIT=13,FILE='element.fill')
      
      WRITE(13,*) nbad
      DO i = 1,nbad
        WRITE(13,*) bad_list(i)
        PRINT*, bad_list(i)
      ENDDO
      CLOSE(13)
      
      

      END PROGRAM create_bad_elem_list