PROGRAM DAGSWEM
    USE SIZES
    USE GLOBAL
    USE DG
    USE DAGSWEM_STATE
    USE DAGSWEM_CL
    USE FSTARPU_MOD
    USE ISO_C_BINDING, ONLY: C_LOC, C_SIZEOF, C_PTR, C_NULL_PTR

    IMPLICIT NONE
    
    ! Orchestrator variables
    INTEGER :: I, ret
    TYPE SubdomainData
        TYPE(C_PTR) :: handles(NUM_STATE_HANDLES)
    END TYPE
    TYPE(SubdomainData), ALLOCATABLE, TARGET :: subdomains(:)
    LOGICAL :: fileFound
    
    ! Find total number of subdomains (MNPROC)
    MNPROC = 0
    DO
        WRITE(DIRNAME(3:6),'(I4.4)') MNPROC
        INQUIRE(FILE=TRIM(DIRNAME)//'/fort.14', EXIST=fileFound)
        IF (.NOT. fileFound) EXIT
        MNPROC = MNPROC + 1
    END DO
    
    IF (MNPROC == 0) THEN
        PRINT *, "No PE**** directories found."
        STOP
    END IF
    
    PRINT *, "Found ", MNPROC, " subdomains."
    
    ALLOCATE(subdomains(0:MNPROC-1))
    
    ! Initialize each subdomain
    DO I = 0, MNPROC - 1
        MYPROC = I
        CALL MAKE_DIRNAME()
        
        ! Call DG-SWEM initialization routines for this subdomain
        ! This populates pointers in GLOBAL and DG
        PRINT *, "Initializing Subdomain ", I
        CALL COLDSTART()
        
        ! Save context (pointers are saved, and nullified)
        CALL DGSWEM_STATE_REGISTER(subdomains(I)%handles)
    END DO
    
    ! StarPU initialization
    ret = fstarpu_init(C_NULL_PTR)
    
    ! TODO: register handles, run task graph, wait, shutdown
    
    CALL fstarpu_shutdown()
    
END PROGRAM DAGSWEM
