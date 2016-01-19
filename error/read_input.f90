      SUBROUTINE read_input()

      USE globals, ONLY: coarse,fine,base,lines,tf,exclude_bndel

      IMPLICIT NONE
      
      INTEGER, PARAMETER :: ninp = 14
      INTEGER :: inp_read,skipped
      CHARACTER(100) :: temp
      
      OPEN(unit=15,file='error.inp')
      
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", "             Input Information               "
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", " "
      
      inp_read = 0
      skipped = 0
      DO WHILE (inp_read < ninp)
      
        READ(15,"(A100)") temp
                    
        IF ( INDEX(temp,"!") == 1 .or. INDEX(temp,"          ") == 1) THEN
!           PRINT*, "input ignored"
            skipped = skipped + 1
        ELSE
          inp_read = inp_read + 1
          SELECT CASE (inp_read)
            CASE (1)
              coarse%grid_file = TRIM(temp)
              PRINT("(A,A)"), "coarse grid_file = ", coarse%grid_file
            CASE (2)
              coarse%out_direc = TRIM(temp)
              PRINT("(A,A)"), "coarse output directory = ", coarse%out_direc
            CASE (3)
              READ(temp,*) coarse%p
              PRINT("(A,I3)"), "coarse p = ", coarse%p
            CASE (4)
              READ(temp,*) coarse%ctp
              PRINT("(A,I3)"), "coarse ctp = ", coarse%ctp
            CASE (5)
              READ(temp,*) coarse%dt
              PRINT("(A,F15.7)"), "coarse dt = ", coarse%dt
              PRINT*, " " 
            CASE (6)
              fine%grid_file = TRIM(temp)
              PRINT("(A,A)"), "fine grid_file = ", fine%grid_file
            CASE (7)
              fine%out_direc = TRIM(temp)
              PRINT("(A,A)"), "fine output directory = ", fine%out_direc
            CASE (8)
              READ(temp,*) fine%p
              PRINT("(A,I3)"), "fine p = ", fine%p           
            CASE (9)
              READ(temp,*) fine%ctp
              PRINT("(A,I3)"), "fine ctp = ", fine%ctp
            CASE (10)
              READ(temp,*) fine%dt
              PRINT("(A,F15.7)"), "fine dt = ", fine%dt 
              PRINT*, " " 
            CASE (11)
              READ(temp,*) tf
              tf = tf*86400d0
              PRINT("(A,F15.7)"), "tf = ", tf     
            CASE (12)
              READ(temp,*) lines
              PRINT("(A,I5)"), "lines = ", lines  
              PRINT*, " " 
            CASE (13)
              base%grid_file = TRIM(temp)
              PRINT("(A,A)"), "base grid directory = ", base%grid_file    
            CASE (14)
              READ(temp,*) base%ctp
              PRINT("(A,I3)"), "base ctp = ", base%ctp                        
          END SELECT
            
        ENDIF
      
      
      ENDDO
      
      PRINT*, " "
      PRINT("(A,I5)"), "Lines skipped: ", skipped
      PRINT*, " "  
      
      IF(coarse%p /= fine%p) THEN
        PRINT("(A)"), "Warning: inconsistent p"
      ENDIF
      
      CLOSE(15)
      
      fine%sol_name = "fine"
      coarse%sol_name = "coarse"
      base%sol_name = "base"
      
      exclude_bndel = .false.
      IF (base%ctp > 1 .or. coarse%ctp > 1 .or. fine%ctp > 1) THEN
        exclude_bndel = .true.
      ENDIF

      END SUBROUTINE  read_input