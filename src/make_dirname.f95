SUBROUTINE MAKE_DIRNAME()
  INTEGER :: LNAME, IARGC, ARGCOUNT, I, iprefix
  CHARACTER(2048) :: CMDLINEARG
  CHARACTER(8)    :: PREFIX(2) = (/ '/PE0000 ', '/DOM0000' /)
  logical         :: fileFound
  
  INPUTDIR  = ""
  GLOBALDIR = ""
  LOCALDIR  = ""
      ARGCOUNT  = IARGC()
      
      !      WRITE_LOCAL_FILES = MNPROC == 1
      !asey 121128: The DG output will be written to local files anyway.
      !             This setting only affects the SWAN files.
      !     WRITE_LOCAL_FILES = 1
      WRITE_LOCAL_FILES = .FALSE.
      WRITE_LOCAL_HOT_START_FILES = .TRUE.
      
#ifdef CMPI
      WRITE(DIRNAME(3:6),'(I4.4)') MYPROC
#else
      MYPROC=0
#endif
      
      IF (ARGCOUNT > 0) THEN
         I = 0
         DO WHILE (I < ARGCOUNT)
            I = I + 1
            CALL GETARG( I, CMDLINEARG )
            IF (CMDLINEARG(1:2) == "-I") THEN
               I = I + 1
               CALL GETARG( I, INPUTDIR )
            ELSEIF (CMDLINEARG(1:2) == "-O") THEN
               I = I + 1
               CALL GETARG(I,GLOBALDIR)
          ELSEIF (CMDLINEARG(1:2) == "-L") THEN
             WRITE_LOCAL_FILES = .TRUE.
          ENDIF
       ENDDO
    ENDIF
    
    !.....Default root working directory
    
    IF (INPUTDIR == "") THEN
       ROOTDIR = '.'
       INPUTDIR = '.'
    ELSE
       ROOTDIR = INPUTDIR
    ENDIF
    
    !.....Set the global input directory
    
    GBLINPUTDIR = ROOTDIR
    
    !asey 121128: Uncommented these lines.
    if (GLOBALDIR == "") then
       ROOTDIR = '.'
    else
       ROOTDIR = GLOBALDIR
    endif
    
    WRITE(GLOBALDIR,'(A)') TRIM(ROOTDIR)
    
    RETURN
  END SUBROUTINE MAKE_DIRNAME
  
