      MODULE read_dginp
      
      ! Module containing subroutines to read the fort.dg file
      !
      !  - Supports both fixed and keyword formats
      !
      !  - Keyword format is intended to increase flexibility in adding/depreciating features
      !    while maintaining forward compatibility and some degree of backward compatibility 
      !    for the fort.dg file.  
      !
      !      * The keywords are configured in the FORT_DG_SETUP subroutine.
      !        (See subrotine header for details.)
      !
      !      * keyword fort.dg file format rules are:
      !
      !          1) options are assigned in keyword = value format (e.g. fluxtype = 1)
      !          2) one option per line
      !          3) options can be specified in any order
      !          4) line beginning with ! indicates a comment
      !          5) blank lines are skipped
      !          6) comments may follow an assignment (e.g. fluxtype = 1 !comment)
      !          7) unrecognized keyword assignments are skipped      
      !          8) unassigned options that are not required will use 
      !             default values specified in FORT_DG_SETUP
      !          
      !  - Subroutines contained are:
      !
      !      1) READ_FIXED_FORT_DG
      !         
      !         * reads old fixed format fort.dg used in commit 6ada0d60 and eariler
      !
      !      2) READ_KEYWORD_FORT_DG
      !
      !         * reads keyword format fort.dg described above
      !
      !      3) CHECK_ERRORS
      !
      !         * Handles missing options 
      !         * Terminates if required options are missing
      !         * Warns that default values are used for missing optional options and continues      
      !
      !      4) FORT_DG_SETUP
      !
      !         * Responsible for configuring fort.dg options
      !         * MODIFICATIONS FOR ADDITION/REMOVAL OF FORT.DG OPTIONS SHOULD BE DONE HERE
      
      USE globals, ONLY: rp
      USE quit, ONLY: abort
      
      TYPE :: key_val
        CHARACTER(15) :: key            ! keyword
        REAL(rp), POINTER :: rptr       ! pointer to real target
        INTEGER, POINTER :: iptr        ! pointer to integer target
        CHARACTER(100), POINTER :: cptr ! pointer to character target      
        
        INTEGER :: vartype              ! target type indicator: 1=integer, 2=real, 3=character
        
        LOGICAL :: required             ! required/optional flag
        
        INTEGER :: flag                 ! successful read flag
      END TYPE key_val
      
      INTEGER, PARAMETER :: maxopt = 100          ! maximum allowable fort.dg options
      TYPE(key_val), DIMENSION(maxopt) :: dginp      
      
      INTEGER :: nopt                             ! number of valid options in dginp structure
      INTEGER, DIMENSION(maxopt) :: dginp_ind    ! indicies of valid options in dginp structure
      
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Input file variables
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CHARACTER(100), TARGET :: grid_file
      CHARACTER(100), TARGET :: forcing_file
      CHARACTER(100), TARGET :: out_direc 
      CHARACTER(100), TARGET :: bathy_file
      CHARACTER(100), TARGET :: curve_file
      CHARACTER(100), TARGET :: stations_file
      INTEGER, TARGET :: p
      INTEGER, TARGET :: ctp
      INTEGER, TARGET :: hbp
      INTEGER, TARGET :: rk_type 
      INTEGER, TARGET  :: npart 
      INTEGER, TARGET :: coord_sys      
      INTEGER, TARGET :: sol_opt
      INTEGER, TARGET :: sta_opt
      REAL(rp), TARGET :: sol_snap
      REAL(rp), TARGET :: sta_snap
      REAL(rp), TARGET  :: dt
      REAL(rp), TARGET  :: tf
      REAL(rp), TARGET  :: dramp
      REAL(rp), TARGET  :: cf
      REAL(rp), TARGET  :: lines
      REAL(rp), TARGET :: slam0,sphi0
      REAL(rp), TARGET :: h0   
      REAL(rp), TARGET :: esl
      LOGICAL :: hb_file_exists
      LOGICAL :: cb_file_exists
 
          
 
      CONTAINS 
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
      SUBROUTINE read_input(myrank,dirname)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: myrank
      CHARACTER(*), INTENT(IN) :: dirname
      
      INTEGER :: opt_read
      INTEGER :: read_stat
      INTEGER :: skipped
      INTEGER :: read_file
      LOGICAL :: file_exists
      CHARACTER(100) :: temp
      
      ! initialize the dginp option structure
      CALL dginp_setup(myrank)      

      INQUIRE(FILE=TRIM(ADJUSTL(dirname))//'/'//'dgswe.inp', EXIST = file_exists)
      IF(file_exists == .FALSE.) THEN
        PRINT*, "dgswe.inp file does not exist"
        CALL abort()
      ENDIF      
      
      OPEN(UNIT=15,FILE=TRIM(ADJUSTL(dirname))//'/'//'dgswe.inp') 
      
      opt_read = 0
      skipped = 0 
      read_file = 0

      DO WHILE (opt_read < nopt)
      
        READ(15,"(A100)",IOSTAT=read_stat) temp
        IF(read_stat /= 0) THEN                    ! check for end-of-file
          EXIT
        ENDIF
        
        IF ( INDEX(temp,"!") == 1 .or. LEN(TRIM(temp)) == 0) THEN
            skipped = skipped + 1
        ELSE  
          read_file = 1
          EXIT
        ENDIF
        
      ENDDO
      
      CLOSE(15)
      

      IF (myrank == 0) THEN
        PRINT "(A)", "---------------------------------------------"
        PRINT "(A)", "             Input Information               "
        PRINT "(A)", "---------------------------------------------"
        PRINT "(A)", " "            
      ENDIF        
      

      IF (read_file == 1) THEN
        IF (INDEX(temp,"=") > 0) THEN
          IF (myrank == 0) THEN
            PRINT*, "Reading keyword format input file"
          ENDIF
          CALL read_keyword_dginp(myrank,dirname)
        ELSE
          IF (myrank == 0) THEN
            PRINT*, "Reading fixed format input file"
          ENDIF
          CALL read_fixed_dginp(myrank,dirname)
        ENDIF
      ELSE 
        PRINT*, "ERROR: dgswe.inp does not contain any information"
        CALL abort()
      ENDIF
          
           
      
      
      END SUBROUTINE read_input
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      

      SUBROUTINE read_fixed_dginp(myrank,dirname)

      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: myrank
      CHARACTER(*), INTENT(IN) :: dirname      
      
      INTEGER, PARAMETER :: ninp = 12
      INTEGER :: i
      INTEGER :: exind
      INTEGER :: inp_read,skipped
      CHARACTER(100) :: temp
      CHARACTER(100) :: line
      LOGICAL :: file_exists
      
!       OPEN(UNIT=16, FILE=dirname//'/'//'fort.16')      
      
      OPEN(UNIT=15,FILE=TRIM(ADJUSTL(dirname))//'/'//'dgswe.inp',POSITION="rewind")      
      
      inp_read = 0
      skipped = 0
      DO WHILE (inp_read < ninp)
      
        READ(15,"(A100)") temp
                    
        IF ( INDEX(temp,"!") == 1 .or. LEN(TRIM(temp)) == 0) THEN
            skipped = skipped + 1
        ELSE
        
          exind = INDEX(temp,"!")
          IF (exind > 0) THEN
            line = ADJUSTL(temp(1:exind-1)) 
          ELSE
            line = ADJUSTL(temp)          
          ENDIF
          
          inp_read = inp_read + 1
          SELECT CASE (inp_read)
            CASE (1)
              grid_file = TRIM(line)
            CASE (2)
              forcing_file = TRIM(line)
            CASE (3)
              READ(line,*) p
            CASE (4)
              READ(line,*) ctp
            CASE (5)
              READ(line,*) hbp            
            CASE (6)
              READ(line,*) dt
            CASE (7)
              READ(line,*) tf
!               tf = tf*86400d0
            CASE (8)
              READ(line,*) dramp
            CASE (9)
              READ(line,*) cf
            CASE (10)
              READ(line,*) lines
            CASE (11)
              out_direc = TRIM(line)
            CASE (12)
              READ(line,*) npart
          END SELECT
            
        ENDIF
      
      
      ENDDO
           
      CLOSE(15)
      
 
      
      IF (myrank == 0) THEN     
      
      ! print inputs
      DO i = 1,maxopt        
        IF (ASSOCIATED(dginp(i)%iptr)) THEN  
          PRINT("(A,A,I8)"), dginp(i)%key," = ",dginp(i)%iptr
        ENDIF
        
        IF (ASSOCIATED(dginp(i)%rptr)) THEN 
          PRINT("(A,A,E21.8)"), dginp(i)%key," = ",dginp(i)%rptr
        ENDIF
        
        IF (ASSOCIATED(dginp(i)%cptr)) THEN 
          PRINT("(A,A,A)"), dginp(i)%key," = ",dginp(i)%cptr 
        ENDIF
      ENDDO      
      
      PRINT*, " "
      
      ENDIF

      RETURN
      END SUBROUTINE read_fixed_dginp
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 
      SUBROUTINE read_keyword_dginp(myrank,dirname)    

      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: myrank
      CHARACTER(*), INTENT(IN) :: dirname      
      
      INTEGER :: i,j,opt
      INTEGER :: read_stat
      INTEGER :: opt_read
      INTEGER :: comment,blank
      INTEGER :: eqind,exind
      LOGICAL :: found      
      CHARACTER(100) :: temp,line,empty
      CHARACTER(15) :: test_opt
      CHARACTER(100) :: test_val     
      
      opt_read = 0
      comment = 0 
      blank = 0
      
      
      OPEN(UNIT=25,FILE=TRIM(ADJUSTL(dirname))//'/'//'dgswe.inp',POSITION="rewind")        
     
      
      DO WHILE (opt_read < nopt)
      
        READ(25,"(A100)",IOSTAT=read_stat) temp
        IF(read_stat /= 0) THEN                    ! check for end-of-file
          EXIT
        ENDIF
        
        line = ADJUSTL(temp)
        
        IF(INDEX(line,"!") == 1) THEN              ! lines beginning with ! are skipped
        
          comment = comment + 1
          
        ELSE IF (LEN(TRIM(line)) == 0) THEN        ! blank lines are skipped
        
          blank = blank + 1
          
        ELSE
  
          ! determine keyword and assignment value
          eqind = INDEX(line,"=")
          exind = INDEX(line,"!")
          test_opt = line(1:eqind-1)
          IF (exind > 0) THEN                          ! handle trailing comment 
            test_val = ADJUSTL(line(eqind+1:exind-1))  ! (only necessary if there is no space between value and the !)
          ELSE
            test_val = ADJUSTL(line(eqind+1:))
          ENDIF         
          
          ! Look for a match for the keyword
          found = .false.
    test: DO opt = 1,nopt
    
            i = dginp_ind(opt)    
    
            IF (test_opt == dginp(i)%key) THEN
              
              ! Set variables equal to value from fort.dg through pointer using an internal read
              SELECT CASE (dginp(i)%vartype) 
                CASE(1)
                  READ(test_val,*) dginp(i)%iptr
                  IF (myrank == 0) PRINT("(A,A,I8)"), test_opt," = ",dginp(i)%iptr
                CASE(2)
                  READ(test_val,*) dginp(i)%rptr
                  IF (myrank == 0) PRINT("(A,A,E21.8)"), test_opt," = ",dginp(i)%rptr                  
                CASE(3)
                  dginp(i)%cptr = TRIM(test_val)
                  IF (myrank == 0) PRINT("(A,A,A)"), test_opt," = ",dginp(i)%cptr                  
              END SELECT

              found = .true.          ! flag match
              opt_read = opt_read + 1
              dginp(i)%flag = 1      ! flag option as found
              
              EXIT test
              
            ENDIF
          ENDDO test
                    
          IF (found == .false. .and. eqind > 0) THEN
            ! unmatched lines with an equal sign are either incorrect or no longer supported
            IF (myrank == 0) PRINT("(3A)"),"*** WARNING: ",test_opt, " is an incorrect or depreciated value ***"            
          ELSE IF (found == .false.) THEN
            ! unmatched lines without an equal sign are ignored
            IF (myrank == 0) PRINT("(A)"), "*** WARNING: non-comment line does not contain a keyword assignment***"           
          ENDIF
          
        ENDIF
      ENDDO 
      
      IF (myrank == 0) PRINT*, ""
     
      CALL check_errors(myrank,opt_read)
      
      IF (myrank == 0) PRINT*, ""
      CLOSE(25)
            
      END SUBROUTINE read_keyword_dginp
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     SUBROUTINE check_errors(myrank,opt_read)

     IMPLICIT NONE
     
     INTEGER, INTENT(IN) :: myrank
     INTEGER, INTENT(IN) :: opt_read     
     
     INTEGER :: i,j,opt
     INTEGER :: quit
     INTEGER :: l
     
     IF(opt_read /= nopt) THEN

       ! check for required options that are unspecifed 
       quit = 0
       DO opt = 1,nopt
         i = dginp_ind(opt)
         IF (dginp(i)%flag == 0 .and. dginp(i)%required == .true.) THEN
           quit = 1   ! flag fatal error
         ENDIF
       ENDDO
        
       IF (quit == 1) THEN
        
          IF (myrank == 0) PRINT("(A)"), "*** ERROR: There are missing required options in the fort.dg file ***"  
          IF (myrank == 0) PRINT("(A)"), "           The following options must be specified: "      
          j = 0        
          DO opt = 1,nopt
            i = dginp_ind(opt)
            IF (dginp(i)%flag == 0 .and. dginp(i)%required == .true.) THEN
              j = j+1
              IF (myrank == 0) PRINT "(A,I3,2A)", "              ",j,") ",dginp(i)%key
            ENDIF
          ENDDO          
          
          IF (myrank == 0) PRINT("(A)"), "!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!"
          STOP
          
       ELSE
        
          IF (myrank == 0) PRINT("(A)"), "*** WARNING: There are missing optional options in the fort.dg file ***"
          IF (myrank == 0) PRINT("(A)"), "             The following default values will be used: "    
          j = 0        
          DO opt = 1,nopt
            i = dginp_ind(opt)
            IF (dginp(i)%flag == 0 .and. dginp(i)%required == .false.) THEN
              
              j = j+1
              SELECT CASE (dginp(i)%vartype) 
                CASE(1)
                  IF (myrank == 0) PRINT("(A,I3,A,A,A,I8)"),     "              ",j,") ",dginp(i)%key," = ",dginp(i)%iptr
                CASE(2)
                  IF (myrank == 0) PRINT("(A,I3,A,A,A,E21.8)"),  "              ",j,") ",dginp(i)%key," = ",dginp(i)%rptr                  
                CASE(3)
                  IF (myrank == 0) PRINT("(A,I3,A,A,A,A)"),      "              ",j,") ",dginp(i)%key," = ",dginp(i)%cptr                  
              END SELECT
              
            ENDIF
          ENDDO 
          
          IF (myrank == 0) PRINT("(A)"), '!!!!!! EXECUTION WILL CONTINUE !!!!!!!!'
          
       ENDIF       
                  
     ENDIF        
     
     ! make sure out_direc ends with /
     l = len(trim(out_direc))
     IF (out_direc(l:l) /= "/") THEN
       out_direc = out_direc(1:l) // "/"
       PRINT*, "Added / to out_direc: ", out_direc
       
     ENDIF
     
     
     RETURN
     END SUBROUTINE check_errors
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!           
      
      SUBROUTINE dginp_setup(myrank)
      
      ! Subroutine that configures the fort.dg options
      !
      !   This subroutine is meant to add flexibility in adding/depreciating  
      !   features while maintaining forward (and some degree of backward) compatibility 
      !
      !   - Options can be added to the fort.dg file by:
      !     
      !     At the top of this (read_dginp) module:
      !       1) Declaring the variable with the TARGET attribute 
      !       
      !     In the dginp structure (found in this subroutine):
      !       2) Specifying a keyword in a unused index (<= maxopt) 
      !       3) Associating the appropriate pointer with the corresponding variable
      !          Note: pointer must agree with the associated variable type 
      !                (iptr=integer, rptr=real, cptr=character)      
      !       4) Specifying whether the variable is required 
      !       5) Providing a default value
      ! 
      !   - Options can be removed from the fort.dg file by:
      !       1) Commenting out or deleting an existing entry in the dginp structure
      !          Note: re-indexing subsequent entries is not necessary 
      !       
      !       OR
      !
      !       2) Setting the dginp(i)%required variable to .false.
      ! 
      !   - New features should be added as dginp(i)%required = .false. as much as possible 
      !     to maintain backward compatibility, older fort.dg files not containing these 
      !     options will cause provided default values to be used (these should be set so 
      !     the feature is turned off)
      !
      !   - fort.dg files containing new feature options can still be used for previous  
      !     versions of the code because the new options will be ignored
          
      
      IMPLICIT NONE        
      
      INTEGER, INTENT(IN) :: myrank
      
      INTEGER :: i
      INTEGER :: ncheck
      CHARACTER(15) :: empty
      
      ! initialize dginp structure
      DO i = 1,maxopt
        NULLIFY(dginp(i)%iptr)
        NULLIFY(dginp(i)%rptr)
        NULLIFY(dginp(i)%cptr)
        
        dginp(i)%key = empty
        dginp(i)%flag = 0        
      ENDDO
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Configure fort.dg options here:
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      !    keywords                         target variables                      requirement                 default values
      dginp(1)%key = "grid_file";       dginp(1)%cptr => grid_file ;      dginp(1)%required = .true.;      dginp(1)%cptr = ""
      dginp(2)%key = "forcing_file";    dginp(2)%cptr => forcing_file;    dginp(2)%required = .true.;      dginp(2)%cptr = ""
      dginp(3)%key = "p";               dginp(3)%iptr => p;               dginp(3)%required = .true.;      dginp(3)%iptr = 1
      dginp(4)%key = "ctp";             dginp(4)%iptr => ctp;             dginp(4)%required = .false.;     dginp(4)%iptr = 1
      dginp(5)%key = "hbp";             dginp(5)%iptr => hbp;             dginp(5)%required = .false.;     dginp(5)%iptr = 1
      dginp(6)%key = "rk";              dginp(6)%iptr => rk_type;         dginp(6)%required = .false.;     dginp(6)%iptr = 22      
      dginp(7)%key = "dt";              dginp(7)%rptr => dt;              dginp(7)%required = .true.;      dginp(7)%rptr = 0.5d0
      dginp(8)%key = "tf";              dginp(8)%rptr => tf;              dginp(8)%required = .true.;      dginp(8)%rptr = 1
      dginp(9)%key = "dramp";           dginp(9)%rptr => dramp;           dginp(9)%required = .true.;      dginp(9)%rptr = 0.08d0
      dginp(10)%key = "cf";             dginp(10)%rptr => cf;             dginp(10)%required = .true.;     dginp(10)%rptr = 0.0025d0
!       dginp(11)%key = "lines";          dginp(11)%rptr => lines;          dginp(11)%required = .true.;     dginp(11)%rptr = 10d0
      dginp(12)%key = "out_direc";      dginp(12)%cptr => out_direc;      dginp(12)%required = .false.;    dginp(12)%cptr = "./"
      dginp(13)%key = "npart";          dginp(13)%iptr => npart;          dginp(13)%required = .false.;     dginp(13)%iptr = 1
      dginp(14)%key = "bathy_file";     dginp(14)%cptr => bathy_file;     dginp(14)%required = .false.;    dginp(14)%cptr = "./elem_nodes.d"
      dginp(15)%key = "coord_sys";      dginp(15)%iptr => coord_sys;      dginp(15)%required = .false.;    dginp(15)%iptr = 1
      dginp(16)%key = "slam0";          dginp(16)%rptr => slam0;          dginp(16)%required = .false.;    dginp(16)%rptr = 0d0 
      dginp(17)%key = "sphi0";          dginp(17)%rptr => sphi0;          dginp(17)%required = .false.;    dginp(17)%rptr = 0d0 
      dginp(18)%key = "h0";             dginp(18)%rptr => h0;             dginp(18)%required = .false.;    dginp(18)%rptr = 1d0   
      dginp(19)%key = "curve_file";     dginp(19)%cptr => curve_file;     dginp(19)%required = .false.;    dginp(19)%cptr = "./cl_nodes.cb"
      dginp(20)%key = "sol_opt";        dginp(20)%iptr => sol_opt;        dginp(20)%required = .false.;    dginp(20)%iptr = 0  
      dginp(21)%key = "sta_opt";        dginp(21)%iptr => sta_opt;        dginp(21)%required = .false.;    dginp(21)%iptr = 0  
      dginp(22)%key = "sol_snap";       dginp(22)%rptr => sol_snap;       dginp(22)%required = .false.;    dginp(22)%rptr = 0  
      dginp(23)%key = "sta_snap";       dginp(23)%rptr => sta_snap;       dginp(23)%required = .false.;    dginp(23)%rptr = 0  
      dginp(24)%key = "sta_file";       dginp(24)%cptr => stations_file;  dginp(24)%required = .false.;    dginp(24)%cptr = "./stations.d"  
      dginp(25)%key = "esl";            dginp(25)%rptr => esl;            dginp(25)%required = .false.;    dginp(25)%rptr = 0d0       
!       dginp(25)%key = "??";             dginp(25)%?ptr => ??;             dginp(25)%required = .false.;    dginp(25)%?ptr = ??  
      

      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! End configuration
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      nopt = 0
      ncheck = 0
      DO i = 1,maxopt
      
        ! find and keep track of populated indicies
        IF (dginp(i)%key .ne. empty) THEN      
          nopt = nopt + 1      
          dginp_ind(nopt) = i
        ENDIF
        
        ! determine target variable type by checking association status
        dginp(i)%vartype = 0    
        
        IF (ASSOCIATED(dginp(i)%iptr)) THEN  ! integer
          ncheck = ncheck + 1   
          dginp(i)%vartype = 1
        ENDIF
        
        IF (ASSOCIATED(dginp(i)%rptr)) THEN ! real
          ncheck = ncheck + 1
          dginp(i)%vartype = 2
        ENDIF
        
        IF (ASSOCIATED(dginp(i)%cptr)) THEN ! character
          ncheck = ncheck + 1        
          dginp(i)%vartype = 3        
        ENDIF
      ENDDO
      
!       PRINT*, "Number of options = ", nopt
!       PRINT*, "Number of pointer associations = ", ncheck
      
      ! ensure user has associated each keyword pointer
      IF (nopt /= ncheck) THEN
        IF (myrank == 0) PRINT("(A)"), "*** ERROR: fort.dg option pointer association error ***"
        IF (myrank == 0) PRINT("(A)"), "           check keyword configuration in dginp_setup subroutine"
        STOP
      ENDIF
          
      
      RETURN
      END SUBROUTINE dginp_setup
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     

      SUBROUTINE write_local(pe)
      
      IMPLICIT NONE
            
      INTEGER, INTENT(IN) :: pe
      CHARACTER(6) :: dirname
      INTEGER,PARAMETER :: lname = 6
      
      dirname = "PE0000"      
      WRITE(dirname(3:lname),"(I4.4)") pe-1          
      
      grid_file = dirname(1:lname)//'/'//"fort.14"
      forcing_file = dirname(1:lname)//'/'//"fort.15"      
      out_direc = './' // dirname(1:lname) //'/'
      bathy_file = './' // dirname(1:lname) //'/'//"fort.hb"
      curve_file = './' // dirname(1:lname) //'/'//"fort.cb"      
      stations_file = './' // dirname(1:lname) //'/'//"fort.sta"         
      
            
      OPEN(UNIT=10,FILE=TRIM(ADJUSTL(dirname(1:lname)))//'/'//'dgswe.inp')
      CALL write_input(file_unit=10)
      CLOSE(10)
      
      END SUBROUTINE write_local

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        

      SUBROUTINE write_input(file_unit)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: file_unit
      INTEGER :: i
      

      DO i = 1,maxopt        
        IF (ASSOCIATED(dginp(i)%iptr)) THEN  
          WRITE(file_unit,"(A,A,I8)") dginp(i)%key," = ",dginp(i)%iptr
        ENDIF
        
        IF (ASSOCIATED(dginp(i)%rptr)) THEN 
          WRITE(file_unit,"(A,A,E21.8)") dginp(i)%key," = ",dginp(i)%rptr
        ENDIF
        
        IF (ASSOCIATED(dginp(i)%cptr)) THEN 
          WRITE(file_unit,"(A,A,A)") dginp(i)%key," = ",dginp(i)%cptr 
        ENDIF
      ENDDO      
                
      
      RETURN
      END SUBROUTINE write_input

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      SUBROUTINE write_file_SHAs(file_unit,dirname)
           
      IMPLICIT NONE
           
      INTEGER, INTENT(IN) :: file_unit
      CHARACTER(*) :: dirname
      
      CHARACTER(40) :: sha1
      
      CALL directed_output(" ",file_unit)          
      CALL directed_output("grid file SHA: "//sha1(grid_file,dirname)         ,file_unit)    
      CALL directed_output("bathy file SHA: "//sha1(bathy_file,dirname)       ,file_unit)  
      CALL directed_output("curve file SHA: "//sha1(curve_file,dirname)       ,file_unit)   
      CALL directed_output("forcing file SHA: "//sha1(forcing_file,dirname)   ,file_unit)  
      CALL directed_output("stations file SHA: "//sha1(stations_file,dirname) ,file_unit)        
      
      RETURN
      END SUBROUTINE write_file_SHAs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      END MODULE read_dginp