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
  
END MODULE FORT_DG
