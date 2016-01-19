      MODULE allocation

      USE globals, ONLY: solution
      
      CONTAINS
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
      SUBROUTINE sizes(sol)
      
      IMPLICIT NONE   
      
      INTEGER :: p,ctp
      TYPE(solution) :: sol
      
      p = sol%p
      ctp = sol%ctp

      sol%nverts(1) = 3
      sol%nverts(2) = 4
      sol%nverts(3) = 3
      sol%nverts(4) = 4      
      
      sol%ndof(1) = (p+1)*(p+2)/2
      sol%ndof(2) = (p+1)*(p+1)
      sol%ndof(3) = sol%ndof(1)
      sol%ndof(4) = sol%ndof(2)
      sol%mndof = maxval(sol%ndof)
      
      sol%np(1) = 1
      sol%np(2) = 1
      sol%np(3) = ctp
      sol%np(4) = ctp        
      
      sol%nnds(1) = 3
      sol%nnds(2) = 4
      sol%nnds(3) = (ctp+1)*(ctp+2)/2
      sol%nnds(4) = (ctp+1)*(ctp+1) 
      sol%mnnds = maxval(sol%nnds)         

      RETURN
      END SUBROUTINE sizes      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      

      SUBROUTINE grid_alloc(stage,sol)

      IMPLICIT NONE
      
      INTEGER :: stage
      TYPE(solution) :: sol
      INTEGER :: i      
      INTEGER :: n 
      INTEGER :: alloc_status(10)

      alloc_status = 0
      
      IF (stage == 1) THEN
      
        n = 10
        
        ALLOCATE(sol%ect(sol%mnnds,sol%ne),STAT = alloc_status(1))
        ALLOCATE(sol%vct(4,sol%ne),STAT = alloc_status(2))
        ALLOCATE(sol%xy(2,sol%nn),STAT = alloc_status(3))
        ALLOCATE(sol%depth(sol%nn),STAT = alloc_status(4))
        ALLOCATE(sol%nelnds(sol%ne),STAT = alloc_status(5))
        ALLOCATE(sol%el_type(sol%ne),STAT = alloc_status(6))
        ALLOCATE(sol%elxy(sol%mnnds,sol%ne,2),STAT = alloc_status(7))
        ALLOCATE(sol%elhb(sol%mnnds,sol%ne),STAT = alloc_status(8))
        ALLOCATE(sol%vxyn(sol%nn),STAT = alloc_status(9))
        ALLOCATE(sol%vxy(2,sol%nn),STAT = alloc_status(10))
      
      ELSE IF (stage == 2) THEN
      
        n = 2     
        
        ALLOCATE(sol%obseg(sol%nope),STAT = alloc_status(1))
        ALLOCATE(sol%obnds(sol%neta,sol%nope),STAT = alloc_status(2))
      
      ELSE IF (stage == 3) THEN
      
        n = 2
        
        ALLOCATE(sol%fbseg(2,sol%nbou),STAT = alloc_status(1))
        ALLOCATE(sol%fbnds(sol%nvel,sol%nbou),STAT = alloc_status(2))        
      
      ENDIF
      
      DO i = 1,n
        IF (alloc_status(i) /= 0) THEN
          PRINT*, "Allocation error: grid_alloc"
          PRINT*, "Stage = ", stage
        ENDIF
      ENDDO

      RETURN
      END SUBROUTINE grid_alloc
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      SUBROUTINE nest_alloc()
      
      USE globals, ONLY: fine,elf2elc,elf2elb 
      
      IMPLICIT NONE
      
      ALLOCATE(elf2elc(fine%ne),elf2elb(fine%ne))    
      
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", "          Finding Element Nesting            "
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", " "        
      
      RETURN
      END SUBROUTINE nest_alloc

      SUBROUTINE error_alloc()
      
      USE globals, ONLY: mnqpta,coarse,fine, &
                         Hc,Qxc,Qyc,Hf,Qxf,Qyf, &
                         xf,yf,r,s,hb,phi
      
      IMPLICIT NONE
      
      INTEGER :: ne,mndof
      
      ALLOCATE(xf(mnqpta),yf(mnqpta))
      ALLOCATE(r(mnqpta),s(mnqpta),hb(mnqpta))
      
      mndof = coarse%mndof
      ne = coarse%ne
      
      ALLOCATE(coarse%H(ne,mndof),coarse%Qx(ne,mndof),coarse%Qy(ne,mndof))
      ALLOCATE(coarse%phi(mndof,mnqpta,1))
      ALLOCATE(phi(mndof*mnqpta))
 
      mndof = fine%mndof
      ne = fine%ne
      
      ALLOCATE(fine%H(ne,mndof),fine%Qx(ne,mndof),fine%Qy(ne,mndof)) 
      
      ALLOCATE(Hc(mnqpta),Qxc(mnqpta),Qyc(mnqpta))
      ALLOCATE(Hf(mnqpta),Qxf(mnqpta),Qyf(mnqpta))    
      
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", "            Reading Solutions                "
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", " "      
      
      RETURN
      END SUBROUTINE  error_alloc

      END MODULE allocation

