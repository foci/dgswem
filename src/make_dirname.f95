SUBROUTINE MAKE_DIRNAME(s)
  use sizes
  implicit none
  type (sizes_type) :: s
  INTEGER :: IARGC, ARGCOUNT, I, iprefix
  CHARACTER(2048) :: CMDLINEARG
  CHARACTER(8)    :: PREFIX(2) = (/ '/PE0000 ', '/DOM0000' /)
  logical         :: fileFound
  
  s%INPUTDIR  = ""
  s%GLOBALDIR = ""
  s%LOCALDIR  = ""
  ARGCOUNT  = IARGC()
  
  !      s%WRITE_LOCAL_FILES = s%MNPROC == 1
  !asey 121128: The DG output will be written to local files anyway.
  !             This setting only affects the SWAN files.
  !     s%WRITE_LOCAL_FILES = 1
  s%WRITE_LOCAL_FILES = .FALSE.
  s%WRITE_LOCAL_HOT_START_FILES = .TRUE.
  
#ifdef CMPI
  WRITE(s%DIRNAME(3:6),'(I4.4)') s%MYPROC
#else
  s%MYPROC=0
#endif
  
  IF (ARGCOUNT > 0) THEN
     I = 0
     DO WHILE (I < ARGCOUNT)
        I = I + 1
        CALL GETARG( I, CMDLINEARG )
        IF (CMDLINEARG(1:2) == "-I") THEN
           I = I + 1
           CALL GETARG( I, s%INPUTDIR )
        ELSEIF (CMDLINEARG(1:2) == "-O") THEN
           I = I + 1
           CALL GETARG(I,s%GLOBALDIR)
        ELSEIF (CMDLINEARG(1:2) == "-L") THEN
           s%WRITE_LOCAL_FILES = .TRUE.
        ENDIF
     ENDDO
  ENDIF
  
  !.....Default root working directory
  
  IF (s%INPUTDIR == "") THEN
     s%ROOTDIR = '.'
     s%INPUTDIR = '.'
  ELSE
     s%ROOTDIR = s%INPUTDIR
  ENDIF
  
  !.....Set the global input directory
  
  s%GBLINPUTDIR = s%ROOTDIR
  
#ifdef CMPI
  !      iprefix = 0
  !      do i = 1, 2
  !        INQUIRE(file=TRIM(s%ROOTDIR)//TRIM(PREFIX(i))//'/'//'fort.14',
  !     &    exist=fileFound)
  !       if (fileFound) then
  !          iprefix = i
  !          exit
  !        end if
  !      end do
  !      if (.not. fileFound) then
  !        print *, "Failed to find prefix directory"
  !        call msg_abort()
  !      end if
  
  !      WRITE(s%INPUTDIR,'(2A)') TRIM(s%ROOTDIR),PREFIX(iprefix)
  !      s%LNAME = LEN_TRIM(s%INPUTDIR)
  !      WRITE(s%INPUTDIR(s%LNAME-3:s%LNAME),'(I4.4)') s%MYPROC
#else
  !      WRITE(s%INPUTDIR,'(A)') TRIM(s%ROOTDIR)
#endif
  
  !asey 121128: Uncommented these lines.
  if (s%GLOBALDIR == "") then
     s%ROOTDIR = '.'
  else
     s%ROOTDIR = s%GLOBALDIR
  endif
  
  WRITE(s%GLOBALDIR,'(A)') TRIM(s%ROOTDIR)
  

#ifdef CMPI
  !      WRITE(s%LOCALDIR,'(2A)') TRIM(s%ROOTDIR),TRIM(PREFIX(iprefix))
  !      s%LNAME = LEN_TRIM(s%LOCALDIR)
  !      WRITE(s%LOCALDIR(s%LNAME-3:s%LNAME),'(I4.4)') s%MYPROC
  !      call MAKEDIR(trim(s%LOCALDIR))
#else
  !      WRITE(s%LOCALDIR,'(A)') TRIM(s%ROOTDIR)
#endif
  !      if (s%WRITE_LOCAL_FILES) s%GLOBALDIR = s%LOCALDIR
  !      s%HOTSTARTDIR = s%LOCALDIR
  
  RETURN
END SUBROUTINE MAKE_DIRNAME
