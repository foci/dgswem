      SUBROUTINE read_solution(sol,t)

      USE globals, ONLY:pres,solution,lines

      IMPLICIT NONE
      
      TYPE(solution) :: sol
      
      INTEGER :: i,dof,el,line,et,n
      REAL(pres) :: t,hb
      
      PRINT("(A)"), TRIM(sol%sol_name) // " solution"
      
#ifndef adcirc      
      
      OPEN(UNIT=11, FILE=trim(sol%out_direc) //"solution_H.d")
      OPEN(UNIT=12, FILE=trim(sol%out_direc) //"solution_Qx.d")     
      OPEN(UNIT=13, FILE=trim(sol%out_direc) //"solution_Qy.d")      
      
      READ(11,*) 
      READ(12,*)
      READ(13,*)
      
      DO line = 1,lines+1
        READ(11,*) t
        READ(12,*) t
        READ(13,*) t 
        
        DO dof = 1,sol%mndof
          READ(11,*) (sol%H(el,dof), el = 1,sol%ne)
          READ(12,*) (sol%Qx(el,dof), el = 1,sol%ne)
          READ(13,*) (sol%Qy(el,dof), el = 1,sol%ne)
        ENDDO    
        
      ENDDO
      
      CLOSE(11)
      CLOSE(12)
      CLOSE(13)     
      
#else      
      
      OPEN(UNIT=11, FILE=trim(sol%out_direc) //"DG.63")
      OPEN(UNIT=12, FILE=trim(sol%out_direc) //"DG.64")     
      
      READ(11,*) 
      READ(11,*)
      READ(12,*)
      READ(12,*)
      
      DO line = 1,lines
        READ(11,*) t,n
        READ(12,*) t,n 
        
        DO el = 1,sol%ne
          et = sol%el_type(el)
          DO dof = 1,sol%ndof(et)
            READ(11,*) i,sol%H(el,dof) !, hb
            READ(12,*) i,sol%Qx(el,dof), sol%Qy(el,dof)
          ENDDO    
        ENDDO        
        
      ENDDO
      
      CLOSE(11)
      CLOSE(12)

#endif      

      PRINT("(A)"), "Done"
      PRINT*, " "

      RETURN 
      END SUBROUTINE read_solution
