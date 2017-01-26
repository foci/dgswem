       PROGRAM nc_test

         USE netcdf
          USE globals

          IMPLICIT NONE
          INTEGER :: nc_id
          INTEGER :: NC_DimID_time
          INTEGER :: NC_DimID_ndof
          INTEGER :: NC_DimID_ne

          CALL check(NF90_CREATE('solution_H.nc',NF90_CLOBBER,nc_id))

          CALL check(NF90_DEF_DIM(nc_id,'time',NF90_UNLIMITED,NC_DimID_time))
          CALL check(NF90_DEF_DIM(nc_id,'dof',ndof,NC_DimID_ndof))
          CALL check(NF90_DEF_DIM(nc_id,'ne',ne,NC_DimID_ne))


          
        END PROGRAM nc_test

        SUBROUTINE check(status)
          USE netcdf

          IMPLICIT NONE
          INTEGER, INTENT(IN) :: status

          IF(status /= NF90_NOERR) THEN
            PRINT("(A,A)"), "fatal error from ", TRIM(NF90_STRERROR(status))  
          ENDIF
               
        END SUBROUTINE check
