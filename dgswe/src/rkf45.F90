      SUBROUTINE rkf45()

      USE globals, ONLY: rp,dof,ndof,el,ne,nel_type,blk,npart,npartet,elblk, &
                         t,tstage,dt,pt3333,ramp,dramp, &
                         Hold,H,MirhsH,Qxold,Qx,MirhsQx,Qyold,Qy,MirhsQy

      IMPLICIT NONE
      
      INTEGER :: et
      INTEGER :: flag
      INTEGER :: maxit
      INTEGER :: k
      
      REAL(rp) :: maxdt
      REAL(rp) :: a(6,6),ta(6),b(6),c(5)
      
      
      maxit = 100
      maxdt = 10d0
      
      a = 0d0
      a(1,1) = 1d0/4d0
      
      a(2,1) = 3d0/32d0
      a(2,2) = 9d0/32d0
      
      a(3,1) = 1923d0/2197d0
      a(3,2) = 7200d0/2197d0
      a(3,3) = 7296d0/2198d0
      
      a(4,1) = 439d0/216d0
      a(4,2) = -8d0
      a(4,3) = 3680d0/513d0
      a(4,4) = 845d0/4104d0
      
      a(5,1) = -8d0/27d0     
      a(5,2) = 2d0
      a(5,3) = -3544d0/2565d0
      a(5,4) = 1859d0/4104d0
      a(5,5) = -11d0/40d0
      
      ta(1) = 1d0
      ta(2) = 1d0/4d0
      ta(3) = 3d0/8d0
      ta(4) = 12d0/13d0
      ta(5) = 1d0
      ta(6) = 1d0/2d0
      
      b(1) = 16d0/135d0
      b(2) = 0d0
      b(3) = 6656d0/12825d0 
      b(4) = 28561d0/56430d0
      b(5) = 9d0/5d0
      b(6) = 2d0/55d0
      
      c(1) = 25d0/216d0
      c(2) = 0d0
      c(3) = 1408d0/2565d0
      c(4) = 2197d0/4104d0
      c(5) = 1d0/5d0
      
      


      ! Save previous solution
      DO blk = 1,npart
        DO et = 1,nel_type
          IF (npartet(et,blk) > 0) THEN
          
            DO dof = 1,ndof(et)
!!DIR$ VECTOR ALIGNED
              DO el = elblk(1,blk,et),elblk(2,blk,et)
                Hold(el,dof) = H(el,dof)
              ENDDO
!!DIR$ VECTOR ALIGNED
              DO el = elblk(1,blk,et),elblk(2,blk,et)
                Qxold(el,dof) = Qx(el,dof)
              ENDDO
!!DIR$ VECTOR ALIGNED
              DO el = elblk(1,blk,et),elblk(2,blk,et)
                Qyold(el,dof) = Qy(el,dof)
              ENDDO
            ENDDO
            
          ENDIF
        ENDDO
      ENDDO   
      

      tstage = t
      ramp = TANH((2d0*tstage)/(86400d0*dramp))
      
      CALL rhs2()

      DO blk = 1,npart
        DO et  = 1,nel_type
          IF (npartet(et,blk) > 0) THEN
          
            DO dof = 1,ndof(et)
!!DIR$ VECTOR ALIGNED
              DO el = elblk(1,blk,et),elblk(2,blk,et)
                Hst(el,dof,1) = dt*MirhsH(el,dof)
              ENDDO
!!DIR$ VECTOR ALIGNED
              DO el = elblk(1,blk,et),elblk(2,blk,et)
                Qxst(el,dof,1) = dt*MirhsQx(el,dof)
              ENDDO
!!DIR$ VECTOR ALIGNED
              DO el = elblk(1,blk,et),elblk(2,blk,et)
                Qyst(el,dof,1) = dt*MirhsQy(el,dof)
              ENDDO
            ENDDO
            
          ENDIF
        ENDDO
      ENDDO
            
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      

      DO st = 2,6
      
      tstage = t + ta(st)*dt
      ramp = TANH((2d0*tstage)/(86400d0*dramp))
      
      
      DO blk = 1,npart
        DO et  = 1,nel_type
          IF (npartet(et,blk) > 0) THEN
          
            DO dof = 1,ndof(et)
              DO el = elblk(1,blk,et),elblk(2,blk,et)
                H(el,dof) = Hold(el,dof) 
              ENDDO

              DO el = elblk(1,blk,et),elblk(2,blk,et)
                Qx(el,dof) = Qxold(el,dof) 
              ENDDO

              DO el = elblk(1,blk,et),elblk(2,blk,et)
                Qy(el,dof) = Qyold(el,dof)
              ENDDO
            ENDDO
             
      
          DO k = 1,st-1
          
            DO dof = 1,ndof(et)
              DO el = elblk(1,blk,et),elblk(2,blk,et)
                H(el,dof) = H(el,dof) + a(k)*Hst(el,dof,k)
              ENDDO

              DO el = elblk(1,blk,et),elblk(2,blk,et)
                Qx(el,dof) = Qx(el,dof) + a(k)*Qxst(el,dof,k)
              ENDDO

              DO el = elblk(1,blk,et),elblk(2,blk,et)
                Qy(el,dof) = Qy(el,dof) + a(k)*Qyst(el,dof,k)
              ENDDO
            ENDDO
            
          ENDDO    
          
        ENDIF

      ENDDO         
      ENDDO
      
      
      CALL rhs2()

      DO blk = 1,npart
        DO et  = 1,nel_type
          IF (npartet(et,blk) > 0) THEN
          
            DO dof = 1,ndof(et)
              DO el = elblk(1,blk,et),elblk(2,blk,et)
                Hstage(el,dof,st) = dt*MirhsH(el,dof)
              ENDDO

              DO el = elblk(1,blk,et),elblk(2,blk,et)
                Qxstage(el,dof,st) = dt*MirhsQx(el,dof)
              ENDDO

              DO el = elblk(1,blk,et),elblk(2,blk,et)
                Qystage(el,dof,st) = dt*MirhsQy(el,dof)
              ENDDO
            ENDDO
            
          ENDIF
          
        ENDDO
      ENDDO      
      
      ENDDO
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
       DO blk = 1,npart
        DO et  = 1,nel_type
          IF (npartet(et,blk) > 0) THEN
          
            DO dof = 1,ndof(et)
              DO el = elblk(1,blk,et),elblk(2,blk,et)
                H5(el,dof) = Hold(el,dof) + b(1)*Hst(el,dof,1)
                H4(el,dof) = Hold(el,dof) + c(1)*Hst(el,dof,1)
              ENDDO

              DO el = elblk(1,blk,et),elblk(2,blk,et)
                Qx5(el,dof) = Qxold(el,dof) + b(1)*Qxst(el,dof,1)
                Qx4(el,dof) = Qxold(el,dof) + c(1)*Qxst(el,dof,1)
              ENDDO

              DO el = elblk(1,blk,et),elblk(2,blk,et)
                Qy5(el,dof) = Qyold(el,dof) + b(1)*Qyst(el,dof,1)
                Qy4(el,dof) = Qyold(el,dof) + c(1)*Qyst(el,dof,1)
              ENDDO
            ENDDO
            
            DO k = 3,5
            DO dof = 1,ndof(et)
              DO el = elblk(1,blk,et),elblk(2,blk,et)
                H5(el,dof) = H5(el,dof) + b(k)*Hst(el,dof,k)
                H4(el,dof) = H4(el,dof) + c(k)*Hst(el,dof,k)
              ENDDO

              DO el = elblk(1,blk,et),elblk(2,blk,et)
                Qx5(el,dof) = Qx5(el,dof) + b(k)*Qxst(el,dof,k)
                Qx4(el,dof) = Qx4(el,dof) + c(k)*Qxst(el,dof,k)
              ENDDO

              DO el = elblk(1,blk,et),elblk(2,blk,et)
                Qy5(el,dof) = Qy5(el,dof) + b(k)*Qyst(el,dof,k)
                Qy4(el,dof) = Qy4(el,dof) + c(k)*Qyst(el,dof,k)
              ENDDO
            ENDDO
            ENDDO
            
            DO dof = 1,ndof(et)
              DO el = elblk(1,blk,et),elblk(2,blk,et)
                H5(el,dof) = H5(el,dof) + b(6)*Hst(el,dof,6)
              ENDDO

              DO el = elblk(1,blk,et),elblk(2,blk,et)
                Qx5(el,dof) = Qx5(el,dof) + b(6)*Qxst(el,dof,6)
              ENDDO

              DO el = elblk(1,blk,et),elblk(2,blk,et)
                Qy5(el,dof) = Qy5(el,dof) + b(6)*Qyst(el,dof,6)
              ENDDO
            ENDDO            
            
            
            
          ENDIF
        ENDDO
      ENDDO          




     CALL nan_check()

      RETURN
      END SUBROUTINE RK
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
      
      SUBROUTINE nan_check()
      
      USE globals, ONLY: dof,ndof,el,ne, &
                         H,Qx,Qy, &
                         blk,npart,nel_type,elblk,npartet
      
      IMPLICIT NONE
      
      INTEGER :: et
      
      DO blk = 1,npart
        DO et = 1,nel_type
          IF (npartet(et,blk) > 0) THEN
      
            DO dof = 1,ndof(et)
      
              DO el = elblk(1,blk,et),elblk(2,blk,et)
                IF (H(el,dof) /= H(el,dof)) THEN
                  PRINT*, "NaN detected in H solution"
                  STOP
                ENDIF
              ENDDO
        
              DO el = elblk(1,blk,et),elblk(2,blk,et)
                IF (Qx(el,dof) /= Qx(el,dof)) THEN
                  PRINT*, "NaN detected in Qx solution"
                  STOP
                ENDIF
              ENDDO
        
              DO el = elblk(1,blk,et),elblk(2,blk,et)
                IF (Qy(el,dof) /= Qy(el,dof)) THEN
                  PRINT*, "NaN detected in Qy solution"
                  STOP
                ENDIF
              ENDDO
        
            ENDDO
      
          ENDIF
        ENDDO
      ENDDO
      
      RETURN
      END SUBROUTINE nan_check
