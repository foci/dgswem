      MODULE FORT_DG
      
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
      !         * reads old fixed format fort.dg used in dgswem v11.13/dg-adcirc v22
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
      
      USE sizes, ONLY: sz
      
      TYPE :: key_val
        CHARACTER(15) :: key            ! keyword
        REAL(SZ), POINTER :: rptr       ! pointer to real target
        INTEGER, POINTER :: iptr        ! pointer to integer target
        CHARACTER(100), POINTER :: cptr ! pointer to character target      
        
        INTEGER :: vartype              ! target type indicator: 1=integer, 2=real, 3=character
        
        INTEGER :: required             ! required/optional flag
        
        INTEGER :: flag                 ! successful read flag
      END TYPE key_val
      
      INTEGER, PARAMETER :: maxopt = 100          ! maximum allowable fort.dg options
      TYPE(key_val), DIMENSION(maxopt) :: fortdg      
      
      INTEGER :: nopt                             ! number of valid options in fortdg structure
      INTEGER, DIMENSION(maxopt) :: fortdg_ind    ! indicies of valid options in fortdg structure
 
      CONTAINS            
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      

      SUBROUTINE READ_FIXED_FORT_DG()
      
      USE global, ONLY: dgswe,dg_to_cg,sedflag,reaction_rate,
     & sed_equationX,sed_equationY, rhowat0,vertexslope
      USE sizes, ONLY: myproc,layers,dirname
      USE dg, ONLY: padapt,pflag,gflag,diorism,pl,ph,px,slimit,plimit,
     &      pflag2con1,pflag2con2,lebesgueP,fluxtype,rk_stage,rk_order,
     &      modal_ic,dghot,dghotspool,slopeflag,slope_weight,porosity,
     &      sevdm,mnes,artdif,kappa,s0,uniform_dif,tune_by_hand,
     &      sl2_m,sl2_nyu,sl3_md
            
      
      IMPLICIT NONE
      
      INTEGER :: i
      CHARACTER(256) :: LINE
      
      CALL FORT_DG_SETUP()      

      OPEN(25,FILE=DIRNAME//'/'//'fort.dg',POSITION="rewind")  
      
      IF (myproc == 0) THEN
        PRINT*, ""
        PRINT("(A)"), "READING FIXED FORMAT FORT.DG..."
        PRINT*, ""      
      ENDIF
      
      READ(25,*) DGSWE
      READ(25,*) padapt,pflag
      READ(25,*) gflag,diorism
      READ(25,*) pl,ph,px
      READ(25,*) slimit
      READ(25,*) plimit
      READ(25,*) pflag2con1,pflag2con2,lebesgueP 
      READ(25,*) FLUXTYPE
      READ(25,*) RK_STAGE, RK_ORDER
      READ(25,*) DG_TO_CG
      READ(25,*) MODAL_IC
      READ(25,*) DGHOT, DGHOTSPOOL
      READ(25,"(A256)") LINE
      READ(LINE,*) SLOPEFLAG
      IF(SLOPEFLAG.EQ.2) THEN
         READ(LINE,*) SLOPEFLAG, SL2_M, SL2_NYU
      ENDIF
      IF(SLOPEFLAG.EQ.3) THEN
         READ(LINE,*) SLOPEFLAG, SL2_M, SL2_NYU, SL3_MD
      ENDIF
      IF(SLOPEFLAG.EQ.4) THEN
         READ(LINE,*) SLOPEFLAG,slope_weight
         vertexslope = .True.
      ENDIF
      IF(SLOPEFLAG.EQ.5) THEN
         READ(LINE,*) SLOPEFLAG
         vertexslope = .True.
      ENDIF
      IF(SLOPEFLAG.EQ.6) THEN
         READ(LINE,*) SLOPEFLAG,slope_weight
         vertexslope = .True.
      ENDIF
      IF(SLOPEFLAG.EQ.7) THEN
         READ(LINE,*) SLOPEFLAG,slope_weight
         vertexslope = .True.
      ENDIF
      IF(SLOPEFLAG.EQ.8) THEN
         READ(LINE,*) SLOPEFLAG,slope_weight
         vertexslope = .True.
      ENDIF
      IF(SLOPEFLAG.EQ.9) THEN
         READ(LINE,*) SLOPEFLAG,slope_weight
         vertexslope = .True.
      ENDIF
      IF(SLOPEFLAG.EQ.10) THEN
         READ(LINE,*) SLOPEFLAG,slope_weight
         vertexslope = .True.
      ENDIF
      READ(25,*) SEDFLAG,porosity,SEVDM,layers
      READ(25,*) reaction_rate
      READ(25,*) MNES
      READ(25,*) artdif,kappa,s0,uniform_dif,tune_by_hand
      READ(25,'(a)') sed_equationX
      READ(25,'(a)') sed_equationY
      
      IF(FLUXTYPE.NE.1.AND.FLUXTYPE.NE.2.AND.FLUXTYPE.NE.3.AND.
     & FLUXTYPE.NE.4) THEN
         IF (myproc == 0) THEN 
           PRINT *,'SPECIFIED FLUXTYPE (=', FLUXTYPE,') IS NOT ALLOWED.'
           PRINT *,'EXECUTION WILL BE TERMINATED.'
         ENDIF 
         
         STOP 
      ENDIF     
      
      ! print inputs
      IF (myproc == 0) THEN
        DO i = 1,maxopt        
          IF (ASSOCIATED(fortdg(i)%iptr)) THEN  
            PRINT("(A,A,I8)"), fortdg(i)%key," = ",fortdg(i)%iptr
          ENDIF
        
          IF (ASSOCIATED(fortdg(i)%rptr)) THEN 
            PRINT("(A,A,E21.8)"), fortdg(i)%key," = ",fortdg(i)%rptr
          ENDIF
        
          IF (ASSOCIATED(fortdg(i)%cptr)) THEN 
            PRINT("(A,A,A)"), fortdg(i)%key," = ",fortdg(i)%cptr 
          ENDIF
        ENDDO      
        PRINT*, " "
      ENDIF 
      
      CLOSE(25)
      
      RETURN
      END SUBROUTINE READ_FIXED_FORT_DG
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 
      SUBROUTINE READ_KEYWORD_FORT_DG()
      
      USE sizes, ONLY: myproc,dirname
      USE global, ONLY: nfover

      IMPLICIT NONE
      
      INTEGER :: i,j,opt
      INTEGER :: read_stat
      INTEGER :: opt_read
      INTEGER :: comment,blank
      INTEGER :: eqind,exind
      INTEGER :: found      
      CHARACTER(100) :: temp,line,empty
      CHARACTER(15) :: test_opt
      CHARACTER(100) :: test_val
      
      ! initialize the fortdg option structure
      CALL FORT_DG_SETUP()
      
      opt_read = 0
      comment = 0 
      blank = 0
      
      
      OPEN(25,FILE=DIRNAME//'/'//'fort.dg',POSITION="rewind")   
      IF (myproc == 0) THEN
        PRINT*, ""
        PRINT("(A)"), "READING KEYWORD FORMAT FORT.DG..."
        PRINT*, ""
      ENDIF
     
      
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
          found = 0
       test: DO opt = 1,nopt
    
            i = fortdg_ind(opt)    
    
            IF (test_opt == fortdg(i)%key) THEN
              
              ! Set variables equal to value from fort.dg through pointer using an internal read
              SELECT CASE (fortdg(i)%vartype) 
                CASE(1)
                  READ(test_val,*) fortdg(i)%iptr
                  IF (myproc == 0) PRINT("(A,A,I8)"),
     &              test_opt," = ",fortdg(i)%iptr
                CASE(2)
                  READ(test_val,*) fortdg(i)%rptr
                  IF (myproc == 0) PRINT("(A,A,E21.8)"),
     &             test_opt," = ",fortdg(i)%rptr
                CASE(3)
                  fortdg(i)%cptr = TRIM(test_val)
                  IF (myproc == 0) PRINT("(A,A,A)"),
     &             test_opt," = ",fortdg(i)%cptr
              END SELECT

              found = 1          ! flag match
              opt_read = opt_read + 1
              fortdg(i)%flag = 1      ! flag option as found
              
              EXIT test
              
            ENDIF
          ENDDO test
               
          IF (myproc == 0 ) THEN     
            IF (found == 0 .and. eqind > 0) THEN
              ! unmatched lines with an equal sign are either incorrect or no longer supported
              PRINT("(3A)"),"*** WARNING: ",test_opt,
     &         " is an incorrect or depreciated value ***"
            ELSE IF (found == 0) THEN
              ! unmatched lines without an equal sign are ignored
              PRINT("(A)"),"*** WARNING ***"
            ENDIF
          ENDIF 
           
        ENDIF
      ENDDO 
      
      IF (myproc == 0) PRINT*, ""
     
      CALL CHECK_ERRORS(opt_read)
      
      IF (myproc == 0) PRINT*, ""
      CLOSE(25)
            
      END SUBROUTINE READ_KEYWORD_FORT_DG
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE CHECK_ERRORS(opt_read)
     
      USE sizes, ONLY: myproc     
     
      IMPLICIT NONE
     
      INTEGER :: i,j,opt
      INTEGER :: opt_read
      INTEGER :: quit
     
      IF(opt_read /= nopt) THEN

       ! check for required options that are unspecifed 
       quit = 0
       DO opt = 1,nopt
         i = fortdg_ind(opt)
         IF ((fortdg(i)%flag == 0) .and. (fortdg(i)%required == 1)) THEN
           quit = 1   ! flag fatal error
         ENDIF
       ENDDO
        
       IF (quit == 1) THEN
        
          IF (myproc == 0) THEN
            PRINT("(A)"), "*** ERROR ***"
            PRINT("(A)"), "the following options must be specified: "
            j = 0        
            DO opt = 1,nopt
              i = fortdg_ind(opt)
              IF ((fortdg(i)%flag == 0) .and. (fortdg(i)%required == 1))
     &         THEN
                j = j+1
                PRINT "(A,I3,2A)", "              ",j,") ",fortdg(i)%key
              ENDIF
            ENDDO          
          
            PRINT("(A)"),
     &            "!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!"
          ENDIF
          
          STOP
          
       ELSE
        
          IF (myproc == 0) THEN
            PRINT("(A)"), "*** WARNING: There are missing optional options in the fort.dg file ***"
            PRINT("(A)"), "             The following default values will be used: "    
            j = 0        
            DO opt = 1,nopt
              i = fortdg_ind(opt)
              IF ((fortdg(i)%flag == 0) .and. (fortdg(i)%required == 0)) THEN
              
                j = j+1
                SELECT CASE (fortdg(i)%vartype) 
                  CASE(1)
                    PRINT("(A,I3,A,A,A,I8)"),     "              ",j,") ",fortdg(i)%key," = ",fortdg(i)%iptr
                  CASE(2)
                    PRINT("(A,I3,A,A,A,E21.8)"),  "              ",j,") ",fortdg(i)%key," = ",fortdg(i)%rptr                  
                  CASE(3)
                    PRINT("(A,I3,A,A,A,A)"),      "              ",j,") ",fortdg(i)%key," = ",fortdg(i)%cptr                  
                END SELECT
              
              ENDIF
            ENDDO 
          
            PRINT("(A)"), '!!!!!! EXECUTION WILL CONTINUE !!!!!!!!'
          ENDIF
          
       ENDIF       
                  
      ENDIF
     
     
      RETURN
      END SUBROUTINE CHECK_ERRORS
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!           
      
      SUBROUTINE FORT_DG_SETUP()
      
      ! Subroutine that configures the fort.dg options
      !
      !   This subroutine is meant to add flexibility in adding/depreciating  
      !   features while maintaining forward (and some degree of backward) compatibility 
      !
      !   - Options can be added to the fort.dg file by:
      !       1) Specifying a keyword in a unused index (<= maxopt) of the fortdg structure
      !       2) Associating the appropriate pointer with the corresponding variable
      !          Note: pointer must agree with the associated variable type 
      !                (iptr=integer, rptr=real, cptr=character)      
      !          Note: the associated variable must be declared using the TARGET attribute
      !       3) Specifying whether the variable is required (1 = yes, 0 = no)
      !       4) Providing a default value
      ! 
      !   - Options can be removed from the fort.dg file by:
      !       1) Commenting out or deleting an existing entry in the fortdg structure
      !          Note: re-indexing subsequent entries is not necessary (see fortdg(17) below)
      !       
      !       OR
      !
      !       2) Setting the fortdg(i)%required variable to 0
      ! 
      !   - New features should be added as fortdg(i)%required = 0 as much as possible 
      !     to maintain backward compatibility, older fort.dg files not containing these 
      !     options will cause provided default values to be used (these should be set so 
      !     the feature is turned off)
      !
      !   - fort.dg files containing new feature options can still be used for previous  
      !     versions of the code because the new options will be ignored
      
      
      USE global, ONLY: dgswe,dg_to_cg,sedflag,reaction_rate,
     & sed_equationX,sed_equationY,dgsake,rainfall
      USE sizes, ONLY: myproc,layers
      USE dg, ONLY: padapt,pflag,gflag,diorism,pl,ph,px,slimit,plimit,
     &     pflag2con1,pflag2con2,lebesgueP,fluxtype,rk_stage,rk_order,
     &     modal_ic,dghot,dghotspool,slopeflag,slope_weight,porosity,
     &     sevdm,mnes,artdif,kappa,s0,uniform_dif,tune_by_hand
      
      IMPLICIT NONE        
      
      INTEGER :: i
      INTEGER :: ncheck
      CHARACTER(15) :: empty
      CHARACTER(28) :: sedXdef,sedYdef
      
      ! initialize fortdg structure
      DO i = 1,maxopt
        NULLIFY(fortdg(i)%iptr)
        NULLIFY(fortdg(i)%rptr)
        NULLIFY(fortdg(i)%cptr)
        
        fortdg(i)%key = empty
        fortdg(i)%flag = 0        
      ENDDO
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Configure fort.dg options here:
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      sedXdef = "(ZE_ROE+bed_ROE)**-1 *QX_ROE"
      sedYdef = "(ZE_ROE+bed_ROE)**-1 *QY_ROE"
      
      !    keywords                         target variables                      requirement                 default values
      fortdg(1)%key = "dgswe";          fortdg(1)%iptr => dgswe ;          fortdg(1)%required = 1;   fortdg(1)%iptr = 1
      fortdg(2)%key = "padapt";         fortdg(2)%iptr => padapt;          fortdg(2)%required = 1;   fortdg(2)%iptr = 0
      fortdg(3)%key = "pflag";          fortdg(3)%iptr => pflag;           fortdg(3)%required = 1;   fortdg(3)%iptr = 2
      fortdg(4)%key = "gflag";          fortdg(4)%iptr => gflag;           fortdg(4)%required = 1;   fortdg(4)%iptr = 1
      fortdg(5)%key = "dis_tol";        fortdg(5)%rptr => diorism;         fortdg(5)%required = 1;   fortdg(5)%rptr = 8
      fortdg(6)%key = "pl";             fortdg(6)%iptr => pl;              fortdg(6)%required = 1;   fortdg(6)%iptr = 1
      fortdg(7)%key = "ph";             fortdg(7)%iptr => ph;              fortdg(7)%required = 1;   fortdg(7)%iptr = 1
      fortdg(8)%key = "px";             fortdg(8)%iptr => px;              fortdg(8)%required = 1;   fortdg(8)%iptr = 1
      fortdg(9)%key = "slimit";         fortdg(9)%rptr => slimit;          fortdg(9)%required = 1;   fortdg(9)%rptr = 0.00005
      fortdg(10)%key = "plimit";        fortdg(10)%rptr => plimit;         fortdg(10)%required = 1;  fortdg(10)%rptr = 10
      fortdg(11)%key = "k";             fortdg(11)%rptr => pflag2con1;     fortdg(11)%required = 1;  fortdg(11)%rptr = 1
      fortdg(12)%key = "ks";            fortdg(12)%rptr => pflag2con2;     fortdg(12)%required = 1;  fortdg(12)%rptr = 0.5
      fortdg(13)%key = "L";             fortdg(13)%iptr => lebesgueP;      fortdg(13)%required = 1;  fortdg(13)%iptr = 2
      fortdg(14)%key = "fluxtype";      fortdg(14)%iptr => fluxtype;       fortdg(14)%required = 1;  fortdg(14)%iptr = 1
      fortdg(15)%key = "rk_stage";      fortdg(15)%iptr => rk_stage;       fortdg(15)%required = 1;  fortdg(15)%iptr = 2
      fortdg(16)%key = "rk_order";      fortdg(16)%iptr => rk_order;       fortdg(16)%required = 1;  fortdg(16)%iptr = 2
!       fortdg(17)%key = "dg_to_cg";      fortdg(17)%iptr => dg_to_cg;     fortdg(17)%required = 1;  fortdg(17)%iptr = 1
      fortdg(18)%key = "modal_ic";      fortdg(18)%iptr => modal_ic;       fortdg(18)%required = 1;  fortdg(18)%iptr = 0
      fortdg(19)%key = "dghot";         fortdg(19)%iptr => dghot;          fortdg(19)%required = 1;  fortdg(19)%iptr = 0
      fortdg(20)%key = "dghotspool";    fortdg(20)%iptr => dghotspool;     fortdg(20)%required = 1;  fortdg(20)%iptr = 86400
      fortdg(21)%key = "slopeflag";     fortdg(21)%iptr => slopeflag;      fortdg(21)%required = 1;  fortdg(21)%iptr = 5
      fortdg(22)%key = "weight";        fortdg(22)%rptr => slope_weight;   fortdg(22)%required = 1;  fortdg(22)%rptr = 1
      fortdg(23)%key = "sedflag";       fortdg(23)%iptr => sedflag;        fortdg(23)%required = 1;  fortdg(23)%iptr = 0
      fortdg(24)%key = "porosity";      fortdg(24)%rptr => porosity;       fortdg(24)%required = 1;  fortdg(24)%rptr = 0.0001
      fortdg(25)%key = "sevdm";         fortdg(25)%rptr => sevdm;          fortdg(25)%required = 1;  fortdg(25)%rptr = 0.00001
      fortdg(26)%key = "layers";        fortdg(26)%iptr => layers;         fortdg(26)%required = 0;  fortdg(26)%iptr = 1
      fortdg(27)%key = "rxn_rate";      fortdg(27)%rptr => reaction_rate;  fortdg(27)%required = 1;  fortdg(27)%rptr = 1.0
      fortdg(28)%key = "nelem";         fortdg(28)%iptr => mnes;           fortdg(28)%required = 1;  fortdg(28)%iptr = 23556
      fortdg(29)%key = "artdif";        fortdg(29)%iptr => artdif;         fortdg(29)%required = 1;  fortdg(29)%iptr = 0
      fortdg(30)%key = "kappa";         fortdg(30)%rptr => kappa;          fortdg(30)%required = 1;  fortdg(30)%rptr = -1.0
      fortdg(31)%key = "s0";            fortdg(31)%rptr => s0;             fortdg(31)%required = 1;  fortdg(31)%rptr = 0.0
      fortdg(32)%key = "uniform_dif";   fortdg(32)%rptr => uniform_dif;    fortdg(32)%required = 1;  fortdg(32)%rptr = 2.5e-6
      fortdg(33)%key = "tune_by_hand";  fortdg(33)%iptr => tune_by_hand;   fortdg(33)%required = 1;  fortdg(33)%iptr = 0
      fortdg(34)%key = "sed_equationX"; fortdg(34)%cptr => sed_equationX;  fortdg(34)%required = 0;  fortdg(34)%cptr = sedXdef
      fortdg(35)%key = "sed_equationY"; fortdg(35)%cptr => sed_equationY;  fortdg(35)%required = 0;  fortdg(35)%cptr = sedYdef
      ! 2021-04-13 Younghun: Added for dgsake
      fortdg(36)%key = "dgsake";        fortdg(36)%iptr => dgsake;         fortdg(36)%required = 0;  fortdg(36)%iptr = 0
      fortdg(37)%key = "rainfall";    fortdg(37)%iptr => rainfall;     fortdg(37)%required = 0;  fortdg(37)%iptr = 1
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! End configuration
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      nopt = 0
      ncheck = 0
      DO i = 1,maxopt
      
        ! find and keep track of populated indicies
        IF (fortdg(i)%key .ne. empty) THEN      
          nopt = nopt + 1      
          fortdg_ind(nopt) = i
        ENDIF
        
        ! determine target variable type by checking association status
        fortdg(i)%vartype = 0    
        
        IF (ASSOCIATED(fortdg(i)%iptr)) THEN  ! integer
          ncheck = ncheck + 1   
          fortdg(i)%vartype = 1
        ENDIF
        
        IF (ASSOCIATED(fortdg(i)%rptr)) THEN ! real
          ncheck = ncheck + 1
          fortdg(i)%vartype = 2
        ENDIF
        
        IF (ASSOCIATED(fortdg(i)%cptr)) THEN ! character
          ncheck = ncheck + 1        
          fortdg(i)%vartype = 3        
        ENDIF
      ENDDO
      
!       PRINT*, "Number of options = ", nopt
!       PRINT*, "Number of pointer associations = ", ncheck
      
      ! ensure user has associated each keyword pointer
      IF (nopt /= ncheck) THEN
        IF (myproc == 0) THEN 
          PRINT("(A)"), "*** ERROR: fort.dg option pointer association error ***"
          PRINT("(A)"), "           check keyword configuration in fort_dg_setup subroutine"
        ENDIF 
        
        STOP
      ENDIF
          
      
      RETURN
      END SUBROUTINE FORT_DG_SETUP
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!            
      
      END MODULE FORT_DG
