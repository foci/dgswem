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
        
        LOGICAL :: required             ! required/optional flag
        
        INTEGER :: flag                 ! successful read flag
      END TYPE key_val
      
      INTEGER, PARAMETER :: maxopt = 100          ! maximum allowable fort.dg options
      TYPE(key_val), DIMENSION(maxopt) :: fortdg      
      
      INTEGER :: nopt                             ! number of valid options in fortdg structure
      INTEGER, DIMENSION(maxopt) :: fortdg_ind    ! indicies of valid options in fortdg structure
 
      CONTAINS            
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      

      SUBROUTINE READ_FIXED_FORT_DG(s)
      
      USE global, ONLY: dgswe,dg_to_cg,sedflag,reaction_rate,sed_equationX,sed_equationY, &
                        rhowat0,vertexslope
      USE sizes
      USE dg, ONLY: padapt,pflag,gflag,diorism,pl,ph,px,slimit,plimit, &
                    pflag2con1,pflag2con2,lebesgueP,fluxtype,dg%rk_stage,rk_order, &
                    modal_ic,dghot,dghotspool,slopeflag,slope_weight,porosity, &
                    sevdm,mnes,artdif,kappa,s0,uniform_dif,tune_by_hand, &
                    sl2_m,sl2_nyu,sl3_md
            
      
      IMPLICIT NONE
      type (sizes_type) :: s

      INTEGER :: i
      CHARACTER(256) :: LINE
      
      CALL FORT_DG_SETUP()      

      OPEN(25,FILE=s%DIRNAME//'/'//'fort.dg',POSITION="rewind")  
      
      PRINT*, ""
      PRINT("(A)"), "READING FIXED FORMAT FORT.DG..."
      PRINT*, ""      
      
      READ(25,*) DGSWE
      READ(25,*) dg%padapt,dg%pflag
      READ(25,*) dg%gflag,diorism
      READ(25,*) dg%pl,dg%ph,dg%px
      READ(25,*) slimit
      READ(25,*) plimit
      READ(25,*) pflag2con1,pflag2con2,dg%lebesgueP 
      READ(25,*) DG%FLUXTYPE
      READ(25,*) DG%RK_STAGE, DG%RK_ORDER
      READ(25,*) DG_TO_CG
      READ(25,*) DG%MODAL_IC
      READ(25,*) DG%DGHOT, DG%DGHOTSPOOL
      READ(25,"(A256)") LINE
      READ(LINE,*) DG%SLOPEFLAG
      IF(DG%SLOPEFLAG.EQ.2) THEN
         READ(LINE,*) DG%SLOPEFLAG, SL2_M, SL2_NYU
      ENDIF
      IF(DG%SLOPEFLAG.EQ.3) THEN
         READ(LINE,*) DG%SLOPEFLAG, SL2_M, SL2_NYU, SL3_MD
      ENDIF
      IF(DG%SLOPEFLAG.EQ.4) THEN
         READ(LINE,*) DG%SLOPEFLAG,slope_weight
         vertexslope = .True.
      ENDIF
      IF(DG%SLOPEFLAG.EQ.5) THEN
         READ(LINE,*) DG%SLOPEFLAG
         vertexslope = .True.
      ENDIF
      IF(DG%SLOPEFLAG.EQ.6) THEN
         READ(LINE,*) DG%SLOPEFLAG,slope_weight
         vertexslope = .True.
      ENDIF
      IF(DG%SLOPEFLAG.EQ.7) THEN
         READ(LINE,*) DG%SLOPEFLAG,slope_weight
         vertexslope = .True.
      ENDIF
      IF(DG%SLOPEFLAG.EQ.8) THEN
         READ(LINE,*) DG%SLOPEFLAG,slope_weight
         vertexslope = .True.
      ENDIF
      IF(DG%SLOPEFLAG.EQ.9) THEN
         READ(LINE,*) DG%SLOPEFLAG,slope_weight
         vertexslope = .True.
      ENDIF
      IF(DG%SLOPEFLAG.EQ.10) THEN
         READ(LINE,*) DG%SLOPEFLAG,slope_weight
         vertexslope = .True.
      ENDIF
      READ(25,*) SEDFLAG,porosity,SEVDM,s%layers
      READ(25,*) reaction_rate
      READ(25,*) DG%MNES
      READ(25,*) dg%artdif,kappa,s0,uniform_dif,dg%tune_by_hand
      READ(25,'(a)') sed_equationX
      READ(25,'(a)') sed_equationY
      
      IF(DG%FLUXTYPE.NE.1.AND.DG%FLUXTYPE.NE.2.AND.DG%FLUXTYPE.NE.3.AND.DG%FLUXTYPE.NE.4) THEN
         PRINT *, 'SPECIFIED FLUXTYPE (=', DG%FLUXTYPE,') IS NOT ALLOWED.'
         PRINT *, 'EXECUTION WILL BE TERMINATED.'
         STOP 'SPECIFIED FLUXTYPE IS NOT ALLOWED.'
      ENDIF     
      
      ! print inputs
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
      CLOSE(25)
      
      RETURN
      END SUBROUTINE READ_FIXED_FORT_DG
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 
      SUBROUTINE READ_KEYWORD_FORT_DG(s)
      
      USE sizes
      USE global, ONLY: nfover

      IMPLICIT NONE
      
      type (sizes_type) :: s
      
      INTEGER :: i,j,opt
      INTEGER :: read_stat
      INTEGER :: opt_read
      INTEGER :: comment,blank
      INTEGER :: eqind,exind
      LOGICAL :: found      
      CHARACTER(100) :: temp,line
      CHARACTER(15) :: test_opt
      CHARACTER(100) :: test_val
      
      ! initialize the fortdg option structure
      CALL FORT_DG_SETUP()
      
      opt_read = 0
      comment = 0 
      blank = 0
      
      
      OPEN(25,FILE=s%DIRNAME//'/'//'fort.dg',POSITION="rewind")   
      
      PRINT*, ""
      PRINT("(A)"), "READING KEYWORD FORMAT FORT.DG..."
      PRINT*, ""
     
      
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
    
            i = fortdg_ind(opt)    
    
            IF (test_opt == fortdg(i)%key) THEN
              
              ! Set variables equal to value from fort.dg through pointer using an internal read
              SELECT CASE (fortdg(i)%vartype) 
                CASE(1)
                  READ(test_val,*) fortdg(i)%iptr
                  PRINT("(A,A,I8)"), test_opt," = ",fortdg(i)%iptr
                CASE(2)
                  READ(test_val,*) fortdg(i)%rptr
                  PRINT("(A,A,E21.8)"), test_opt," = ",fortdg(i)%rptr                  
                CASE(3)
                  fortdg(i)%cptr = TRIM(test_val)
                  PRINT("(A,A,A)"), test_opt," = ",fortdg(i)%cptr                  
              END SELECT

              found = .true.          ! flag match
              opt_read = opt_read + 1
              fortdg(i)%flag = 1      ! flag option as found
              
              EXIT test
              
            ENDIF
          ENDDO test
                    
          IF (found .eqv. .false. .and. eqind > 0) THEN
            ! unmatched lines with an equal sign are either incorrect or no longer supported
            PRINT("(3A)"),"*** WARNING: ",test_opt, " is an incorrect or depreciated value ***"            
          ELSE IF (found .eqv. .false.) THEN
            ! unmatched lines without an equal sign are ignored
            PRINT("(A)"), "*** WARNING: non-comment line does not contain a keyword assignment***"           
          ENDIF
          
        ENDIF
      ENDDO 
      
      PRINT*, ""
     
      CALL CHECK_ERRORS(opt_read)
      
      PRINT*, ""
      CLOSE(25)
            
      END SUBROUTINE READ_KEYWORD_FORT_DG
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     SUBROUTINE CHECK_ERRORS(opt_read)
     
     IMPLICIT NONE
     
     INTEGER :: i,j,opt
     INTEGER :: opt_read
     INTEGER :: quit
     
     IF(opt_read /= nopt) THEN

       ! check for required options that are unspecifed 
       quit = 0
       DO opt = 1,nopt
         i = fortdg_ind(opt)
         IF (fortdg(i)%flag == 0 .and. fortdg(i)%required) THEN
           quit = 1   ! flag fatal error
         ENDIF
       ENDDO
        
       IF (quit == 1) THEN
        
          PRINT("(A)"), "*** ERROR: There are missing required options in the fort.dg file ***"  
          PRINT("(A)"), "           The following options must be specified: "      
          j = 0        
          DO opt = 1,nopt
            i = fortdg_ind(opt)
            IF (fortdg(i)%flag == 0 .and. fortdg(i)%required) THEN
              j = j+1
              PRINT "(A,I3,2A)", "              ",j,") ",fortdg(i)%key
            ENDIF
          ENDDO          
          
          PRINT("(A)"), "!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!"
          STOP
          
       ELSE
        
          PRINT("(A)"), "*** WARNING: There are missing optional options in the fort.dg file ***"
          PRINT("(A)"), "             The following default values will be used: "    
          j = 0        
          DO opt = 1,nopt
            i = fortdg_ind(opt)
            IF (fortdg(i)%flag == 0 .and. fortdg(i)%required .eqv. .false.) THEN
              
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
      !       3) Specifying whether the variable is required 
      !       4) Providing a default value
      ! 
      !   - Options can be removed from the fort.dg file by:
      !       1) Commenting out or deleting an existing entry in the fortdg structure
      !          Note: re-indexing subsequent entries is not necessary (see fortdg(17) below)
      !       
      !       OR
      !
      !       2) Setting the fortdg(i)%required variable to .false.
      ! 
      !   - New features should be added as fortdg(i)%required = .false. as much as possible 
      !     to maintain backward compatibility, older fort.dg files not containing these 
      !     options will cause provided default values to be used (these should be set so 
      !     the feature is turned off)
      !
      !   - fort.dg files containing new feature options can still be used for previous  
      !     versions of the code because the new options will be ignored
      
      
      USE global, ONLY: dgswe,dg_to_cg,sedflag,reaction_rate,sed_equationX,sed_equationY
      USE sizes
      USE dg, ONLY: padapt,pflag,gflag,diorism,pl,ph,px,slimit,plimit, &
                    pflag2con1,pflag2con2,lebesgueP,fluxtype,rk_stage,rk_order, &
                    modal_ic,dghot,dghotspool,slopeflag,slope_weight,porosity, &
                    sevdm,mnes,artdif,kappa,s0,uniform_dif,tune_by_hand
      
      IMPLICIT NONE        

      type (sizes_type) :: s
      
      INTEGER :: i
      INTEGER :: ncheck

      !declare targets because we can't have targets in the data types
      integer, target :: layers_target
!      integer, target :: dgflag_target
      integer, target :: dghot_target
      integer, target :: dghotspool_target
      integer, target :: mnes_target
      integer, target :: artdif_target
      integer, target :: tune_by_hand_target
      integer, target :: modal_ic_target
      integer, target :: slopeflag_target
      integer, target :: fluxtype_target
      integer, target :: rk_stage_target
      integer, target :: rk_order_target
      integer, target :: padapt_target
      integer, target :: pflag_target
      integer, target :: pl_target
      integer, target :: ph_target
      integer, target :: px_target
      integer, target :: lebesgueP_target
      integer, target :: gflag_target
      
      real(sz), target :: diorism_target
      real(sz), target :: porosity_target
      real(sz), target :: sevdm_target
      real(sz), target :: slimit_target
      real(sz), target :: plimit_target
      real(sz), target :: pflag2con1_target
      real(sz), target :: pflag2con2_target
      real(sz), target :: slope_weight_target
      real(sz), target :: kappa_target
      real(sz), target :: s0_target
      real(sz), target :: uniform_dif_target
            


      ! initialize fortdg structure
      DO i = 1,maxopt
        NULLIFY(fortdg(i)%iptr)
        NULLIFY(fortdg(i)%rptr)
        NULLIFY(fortdg(i)%cptr)
        
        fortdg(i)%key = "empty"
        fortdg(i)%flag = 0        
      ENDDO
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Configure fort.dg options here:
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      !    keywords                         target variables                      requirement                 default values
      fortdg(1)%key = "dgswe";          fortdg(1)%iptr => dgswe;          fortdg(1)%required = .true.;      fortdg(1)%iptr = 1
      fortdg(2)%key = "padapt";         fortdg(2)%iptr => padapt_target;          fortdg(2)%required = .true.;      fortdg(2)%iptr = 0
      fortdg(3)%key = "pflag";          fortdg(3)%iptr => pflag_target;           fortdg(3)%required = .true.;      fortdg(3)%iptr = 2
      fortdg(4)%key = "gflag";          fortdg(4)%iptr => gflag_target;           fortdg(4)%required = .true.;      fortdg(4)%iptr = 1
      fortdg(5)%key = "dis_tol";        fortdg(5)%rptr => diorism_target;         fortdg(5)%required = .true.;      fortdg(5)%rptr = 8
      fortdg(6)%key = "pl";             fortdg(6)%iptr => pl_target;              fortdg(6)%required = .true.;      fortdg(6)%iptr = 1
      fortdg(7)%key = "ph";             fortdg(7)%iptr => ph_target;              fortdg(7)%required = .true.;      fortdg(7)%iptr = 1
      fortdg(8)%key = "px";             fortdg(8)%iptr => px_target;              fortdg(8)%required = .true.;      fortdg(8)%iptr = 1
      fortdg(9)%key = "slimit";         fortdg(9)%rptr => slimit_target;          fortdg(9)%required = .true.;      fortdg(9)%rptr = 0.00005
      fortdg(10)%key = "plimit";        fortdg(10)%rptr => plimit_target;         fortdg(10)%required = .true.;     fortdg(10)%rptr = 10
      fortdg(11)%key = "k";             fortdg(11)%rptr => pflag2con1_target;     fortdg(11)%required = .true.;     fortdg(11)%rptr = 1
      fortdg(12)%key = "ks";            fortdg(12)%rptr => pflag2con2_target;     fortdg(12)%required = .true.;     fortdg(12)%rptr = 0.5
      fortdg(13)%key = "L";             fortdg(13)%iptr => lebesgueP_target;      fortdg(13)%required = .true.;     fortdg(13)%iptr = 2
      fortdg(14)%key = "fluxtype";      fortdg(14)%iptr => fluxtype_target;       fortdg(14)%required = .true.;     fortdg(14)%iptr = 1
      fortdg(15)%key = "rk_stage";      fortdg(15)%iptr => rk_stage_target;       fortdg(15)%required = .true.;     fortdg(15)%iptr = 2
      fortdg(16)%key = "rk_order";      fortdg(16)%iptr => rk_order_target;       fortdg(16)%required = .true.;     fortdg(16)%iptr = 2
!       fortdg(17)%key = "dg_to_cg";      fortdg(17)%iptr => dg_to_cg;     fortdg(17)%required = .true.;     fortdg(17)%iptr = 1
      fortdg(18)%key = "modal_ic";      fortdg(18)%iptr => modal_ic_target;       fortdg(18)%required = .true.;     fortdg(18)%iptr = 0
      fortdg(19)%key = "dghot";         fortdg(19)%iptr => dghot_target;          fortdg(19)%required = .true.;     fortdg(19)%iptr = 0
      fortdg(20)%key = "dghotspool";    fortdg(20)%iptr => dghotspool_target;     fortdg(20)%required = .true.;     fortdg(20)%iptr = 86400
      fortdg(21)%key = "slopeflag";     fortdg(21)%iptr => slopeflag_target;      fortdg(21)%required = .true.;     fortdg(21)%iptr = 5
      fortdg(22)%key = "weight";        fortdg(22)%rptr => slope_weight_target;   fortdg(22)%required = .true.;     fortdg(22)%rptr = 1
      fortdg(23)%key = "sedflag";       fortdg(23)%iptr => sedflag;        fortdg(23)%required = .true.;     fortdg(23)%iptr = 0
      fortdg(24)%key = "porosity";      fortdg(24)%rptr => porosity_target;       fortdg(24)%required = .true.;     fortdg(24)%rptr = 0.0001
      fortdg(25)%key = "sevdm";         fortdg(25)%rptr => sevdm_target;          fortdg(25)%required = .true.;     fortdg(25)%rptr = 0.00001
      fortdg(26)%key = "layers";        fortdg(26)%iptr => layers_target;         fortdg(26)%required = .false.;    fortdg(26)%iptr = 1
      fortdg(27)%key = "rxn_rate";      fortdg(27)%rptr => reaction_rate;  fortdg(27)%required = .true.;     fortdg(27)%rptr = 1.0
      fortdg(28)%key = "nelem";         fortdg(28)%iptr => mnes_target;           fortdg(28)%required = .true.;     fortdg(28)%iptr = 23556
      fortdg(29)%key = "artdif";        fortdg(29)%iptr => artdif_target;         fortdg(29)%required = .true.;     fortdg(29)%iptr = 0
      fortdg(30)%key = "kappa";         fortdg(30)%rptr => kappa_target;          fortdg(30)%required = .true.;     fortdg(30)%rptr = -1.0
      fortdg(31)%key = "s0";            fortdg(31)%rptr => s0_target;             fortdg(31)%required = .true.;     fortdg(31)%rptr = 0.0;
      fortdg(32)%key = "uniform_dif";   fortdg(32)%rptr => uniform_dif_target;    fortdg(32)%required = .true.;     fortdg(32)%rptr = 2.5e-6;
      fortdg(33)%key = "tune_by_hand";  fortdg(33)%iptr => tune_by_hand_target;   fortdg(33)%required = .true.;     fortdg(33)%iptr = 0;
      fortdg(34)%key = "sed_equationX"; fortdg(34)%cptr => sed_equationX;  fortdg(34)%required = .false.;    fortdg(34)%cptr = "(ZE_ROE+bed_ROE)**-1 *QX_ROE";
      fortdg(35)%key = "sed_equationY"; fortdg(35)%cptr => sed_equationY;  fortdg(35)%required = .false.;    fortdg(35)%cptr = "(ZE_ROE+bed_ROE)**-1 *QY_ROE";
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! End configuration
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! assign targets to the real variables
      s%layers = layers_target
!      dgflag_target
      dg%dghot = dghot_target
      dg%dghotspool = dghotspool_target
      dg%mnes = mnes_target
      dg%artdif = artdif_target
      dg%tune_by_hand = tune_by_hand_target
      dg%modal_ic = modal_ic_target
      dg%slopeflag = slopeflag_target
      dg%fluxtype = fluxtype_target
      dg%rk_stage = rk_stage_target
      dg%rk_order = rk_order_target
      dg%padapt = padapt_target
      dg%pflag = pflag_target
      dg%pl = pl_target
      dg%ph = ph_target
      dg%px = px_target
      dg%lebesegueP = lebesgueP_target
      dg%gflag = gflag_target
      
      dg%diorism = diorism_target
      dg%porosity = porosity_target
      dg%sevdm = sevdm_target
      dg%slimit = slimit_target
      dg%plimit = plimit_target
      dg%pflag2con1 = pflag2con1_target
      dg%pflag2con2 = pflag2con2_target
      dg%slope_weight = slope_weight_target
      dg%kappa = kappa_target
      dg%s0 = s0_target
      dg%uniform_dif = uniform_dif_target

      
      nopt = 0
      ncheck = 0
      DO i = 1,maxopt
      
        ! find and keep track of populated indicies         
        ! zach temp debugging
        IF (fortdg(i)%key .ne. "empty") THEN      
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
      
       PRINT*, "Number of options = ", nopt
       PRINT*, "Number of pointer associations = ", ncheck
      
      ! ensure user has associated each keyword pointer
      IF (nopt /= ncheck) THEN
        PRINT("(A)"), "*** ERROR: fort.dg option pointer association error ***"
        PRINT("(A)"), "           check keyword configuration in fort_dg_setup subroutine"
        STOP
      ENDIF
          
      
      RETURN
      END SUBROUTINE FORT_DG_SETUP
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!            
      
      END MODULE FORT_DG
