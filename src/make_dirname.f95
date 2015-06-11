SUBROUTINE MAKE_DIRNAME(g)
  implicit none
  type (global_data_type) :: g
  INTEGER :: LNAME, IARGC, ARGCOUNT, I, iprefix
  CHARACTER(2048) :: CMDLINEARG
  CHARACTER(8)    :: PREFIX(2) = (/ '/PE0000 ', '/DOM0000' /)
  logical         :: fileFound
  
  g%INPUTDIR  = ""
  g%GLOBALDIR = ""
  g%LOCALDIR  = ""
  g%ARGCOUNT  = g%IARGC()
      
      !      WRITE_LOCAL_FILES = MNPROC == 1
      !asey 121128: The DG output will be written to local files anyway.
      !             This setting only affects the SWAN files.
      !     WRITE_LOCAL_FILES = 1
      g%WRITE_LOCAL_FILES = .FALSE.
      g%WRITE_LOCAL_HOT_START_FILES = .TRUE.
      
#ifdef CMPI
      WRITE(g%DIRNAME(3:6),'(I4.4)') g%MYPROC
#else
      g%MYPROC=0
#endif
      
      IF (g%ARGCOUNT > 0) THEN
         I = 0
         DO WHILE (I < g%ARGCOUNT)
            I = I + 1
            CALL GETARG( I, g%CMDLINEARG )
            IF (CMDLINEARG(1:2) == "-I") THEN
               I = I + 1
               CALL GETARG( I, INPUTDIR )
            ELSEIF (CMDLINEARG(1:2) == "-O") THEN
               I = I + 1
               CALL GETARG(I,GLOBALDIR)
          ELSEIF (CMDLINEARG(1:2) == "-L") THEN
             g%WRITE_LOCAL_FILES = .TRUE.
          ENDIF
       ENDDO
    ENDIF
    
    !.....Default root working directory
    
    IF (g%INPUTDIR == "") THEN
       g%ROOTDIR = '.'
       g%INPUTDIR = '.'
    ELSE
       g%ROOTDIR = g%INPUTDIR
    ENDIF
    
    !.....Set the global input directory
    
    g%GBLINPUTDIR = g%ROOTDIR
    
    !asey 121128: Uncommented these lines.
    if (g%GLOBALDIR == "") then
       g%ROOTDIR = '.'
    else
       g%ROOTDIR = g%GLOBALDIR
    endif
    
    WRITE(g%GLOBALDIR,'(A)') TRIM(g%ROOTDIR)
    
    RETURN
  END SUBROUTINE MAKE_DIRNAME
  
