      SUBROUTINE rk()

      USE globals, ONLY: ndof,ne,nel_type,npartet,elblk, &
                         t,tstage,pt3333,ramp,ark,brk,crk, &
                         Hold,H,MirhsH,Zold,Z,MirhsZ,Qxold,Qx,MirhsQx,Qyold,Qy,MirhsQy  
                         
      USE read_dginp, ONLY: npart,dt,dramp,rk_type
      USE quit, ONLY: abort

      IMPLICIT NONE
      
      INTEGER :: et,stg,dof,el,blk
         
      
   
!       CALL nan_check()
!       PRINT*, "NaN checked initial condition"      

      SELECT CASE (rk_type)
      
        CASE(11)
        
          CALL swap()
          CALL forward_euler()
          
        CASE(22)
        
          CALL swap()
          CALL forward_euler()

          ! Evaluate RHS
          tstage = t + dt
          ramp = TANH((2d0*tstage)/(86400d0*dramp))
          CALL rhs2()
!           CALL rhs3()
!           CALL rhs4()
      
          ! Second RK stage
          DO blk = 1,npart+1
            DO et  = 1,nel_type
              IF (npartet(et,blk) > 0) THEN   
          
                DO dof = 1,ndof(et)
!DIR$ IVDEP                
!!DIR$ VECTOR ALIGNED
                  DO el = elblk(1,blk,et),elblk(2,blk,et)
!                     H(el,dof) = .5d0*(Hold(el,dof) + H(el,dof) + dt*MirhsH(el,dof))
                    Z(el,dof) = .5d0*(Zold(el,dof) + Z(el,dof) + dt*MirhsZ(el,dof))
                    MirhsZ(el,dof) = 0d0
                  ENDDO
!DIR$ IVDEP                  
!!DIR$ VECTOR ALIGNED
                  DO el = elblk(1,blk,et),elblk(2,blk,et)
                    Qx(el,dof) = .5d0*(Qxold(el,dof) + Qx(el,dof) + dt*MirhsQx(el,dof))
                    MirhsQx(el,dof) = 0d0
                  ENDDO
!DIR$ IVDEP                  
!!DIR$ VECTOR ALIGNED
                  DO el = elblk(1,blk,et),elblk(2,blk,et)
                    Qy(el,dof) = .5d0*(Qyold(el,dof) + Qy(el,dof) + dt*MirhsQy(el,dof))
                    MirhsQy(el,dof) = 0d0
                  ENDDO
                ENDDO
      
              ENDIF
            ENDDO
          ENDDO      

        CASE(33)
        
          CALL swap()
          CALL forward_euler()

          ! Evaluate RHS
          tstage = t + dt
          ramp = TANH((2d0*tstage)/(86400d0*dramp))
          CALL rhs2()
!          CALL rhs3()
!          CALL rhs4()
      
          ! Second RK stage
          DO blk = 1,npart+1
            DO et  = 1,nel_type
              IF (npartet(et,blk) > 0) THEN    
          
                DO dof = 1,ndof(et)
!DIR$ IVDEP
!!DIR$ VECTOR ALIGNED
                  DO el = elblk(1,blk,et),elblk(2,blk,et)
!                     H(el,dof) = .25d0*(3d0*Hold(el,dof) + H(el,dof) + dt*MirhsH(el,dof))
                    Z(el,dof) = .25d0*(3d0*Zold(el,dof) + Z(el,dof) + dt*MirhsZ(el,dof))
                    MirhsZ(el,dof) = 0d0
                  ENDDO
!DIR$ IVDEP
!!DIR$ VECTOR ALIGNED
                  DO el = elblk(1,blk,et),elblk(2,blk,et)
                    Qx(el,dof) = .25d0*(3d0*Qxold(el,dof) + Qx(el,dof) + dt*MirhsQx(el,dof))
                    MirhsQx(el,dof) = 0d0
                  ENDDO
!DIR$ IVDEP        
!!DIR$ VECTOR ALIGNED
                  DO el = elblk(1,blk,et),elblk(2,blk,et)
                    Qy(el,dof) = .25d0*(3d0*Qyold(el,dof) + Qy(el,dof) + dt*MirhsQy(el,dof))
                    MirhsQy(el,dof) = 0d0
                  ENDDO
                ENDDO
      
              ENDIF
            ENDDO
          ENDDO   
      
!          CALL nan_check()      

          ! Evaluate RHS
          tstage = t + .5d0*dt
          ramp = TANH((2d0*tstage)/(86400d0*dramp))
          CALL rhs2()
!           CALL rhs3()
!           CALL rhs4()
      
          ! Third RK stage
          DO blk = 1,npart+1
            DO et  = 1,nel_type
              IF (npartet(et,blk) > 0) THEN      
          
                DO dof = 1,ndof(et)
!DIR$ IVDEP      
!!DIR$ VECTOR ALIGNED
                  DO el = elblk(1,blk,et),elblk(2,blk,et)
!                     H(el,dof) = pt3333*(Hold(el,dof) + 2d0*(H(el,dof) + dt*MirhsH(el,dof)))
                    Z(el,dof) = pt3333*(Zold(el,dof) + 2d0*(Z(el,dof) + dt*MirhsZ(el,dof)))
                    MirhsZ(el,dof) = 0d0
                  ENDDO
!DIR$ IVDEP        
!!DIR$ VECTOR ALIGNED
                  DO el = elblk(1,blk,et),elblk(2,blk,et)
                    Qx(el,dof) = pt3333*(Qxold(el,dof) + 2d0*(Qx(el,dof) + dt*MirhsQx(el,dof)))
                    MirhsQx(el,dof) = 0d0
                  ENDDO
!DIR$ IVDEP        
!!DIR$ VECTOR ALIGNED
                  DO el = elblk(1,blk,et),elblk(2,blk,et)
                    Qy(el,dof) = pt3333*(Qyold(el,dof) + 2d0*(Qy(el,dof) + dt*MirhsQy(el,dof)))
                    MirhsQy(el,dof) = 0d0
                  ENDDO
                ENDDO
      
              ENDIF
            ENDDO
          ENDDO      

        CASE(45)
                    
          DO stg = 1,5

          tstage = t + crk(stg)*dt 
          ramp = TANH((2d0*tstage)/(86400d0*dramp))
          CALL rhs2()
      
          DO blk = 1,npart+1
            DO et  = 1,nel_type
              IF (npartet(et,blk) > 0) THEN    
          
                DO dof = 1,ndof(et)
!DIR$ IVDEP
                  DO el = elblk(1,blk,et),elblk(2,blk,et)
                    Zold(el,dof) = ark(stg)*Zold(el,dof) + dt*MirhsZ(el,dof)
                    Z(el,dof) = Z(el,dof) + brk(stg)*Zold(el,dof)
                    MirhsZ(el,dof) = 0d0
                  ENDDO
!DIR$ IVDEP
                  DO el = elblk(1,blk,et),elblk(2,blk,et)
                    Qxold(el,dof) = ark(stg)*Qxold(el,dof) + dt*MirhsQx(el,dof)
                    Qx(el,dof) = Qx(el,dof) + brk(stg)*Qxold(el,dof)
                    MirhsQx(el,dof) = 0d0
                  ENDDO
!DIR$ IVDEP        
                  DO el = elblk(1,blk,et),elblk(2,blk,et)
                    Qyold(el,dof) = ark(stg)*Qyold(el,dof) + dt*MirhsQy(el,dof)
                    Qy(el,dof) = Qy(el,dof) + brk(stg)*Qyold(el,dof)
                    MirhsQy(el,dof) = 0d0
                  ENDDO
                ENDDO
      
              ENDIF
            ENDDO
          ENDDO   
      
          ENDDO
          
        CASE DEFAULT
        
          PRINT*, "Time-stepping option not availiable: ", rk_type
          CALL abort()
          
      END SELECT          
      
      

      CALL nan_check()
         

      RETURN
      END SUBROUTINE RK
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      SUBROUTINE swap()
      
      USE globals, ONLY: ndof,nel_type,npartet,elblk, &
                         Hold,H,Zold,Z,Qxold,Qx,Qyold,Qy   
      USE read_dginp, ONLY: npart                         
                       
      IMPLICIT NONE                       
                       
      INTEGER :: et,dof,el,blk                         
                         
!       CALL nan_check()
!       PRINT*, "NaN checked initial condition"      

      ! Save previous solution
      DO blk = 1,npart+1
        DO et = 1,nel_type
          IF (npartet(et,blk) > 0) THEN
          
            DO dof = 1,ndof(et)
!DIR$ IVDEP            
!!DIR$ VECTOR ALIGNED
              DO el = elblk(1,blk,et),elblk(2,blk,et)
!                 Hold(el,dof) = H(el,dof)
                Zold(el,dof) = Z(el,dof)
              ENDDO
!DIR$ IVDEP              
!!DIR$ VECTOR ALIGNED
              DO el = elblk(1,blk,et),elblk(2,blk,et)
                Qxold(el,dof) = Qx(el,dof)
              ENDDO
!DIR$ IVDEP              
!!DIR$ VECTOR ALIGNED
              DO el = elblk(1,blk,et),elblk(2,blk,et)
                Qyold(el,dof) = Qy(el,dof)
              ENDDO
            ENDDO
            
          ENDIF
        ENDDO
      ENDDO                         
      
      END SUBROUTINE swap
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      SUBROUTINE forward_euler()
      
      USE globals, ONLY: ndof,ne,nel_type,npartet,elblk, &
                         t,tstage,ramp, &
                         Hold,H,MirhsH,Zold,Z,MirhsZ,Qxold,Qx,MirhsQx,Qyold,Qy,MirhsQy    
      USE read_dginp, ONLY: dt,npart,dramp
      
      IMPLICIT NONE
      
      INTEGER :: et,dof,el,blk
      
      ! Evalutate right hand side  
      tstage = t
      ramp = TANH((2d0*tstage)/(86400d0*dramp))
      CALL rhs2()
!       CALL rhs3()
!        CALL rhs4()
      
      ! First RK stage
      DO blk = 1,npart+1
        DO et  = 1,nel_type
          IF (npartet(et,blk) > 0) THEN
          
            DO dof = 1,ndof(et)
!DIR$ IVDEP            
!!DIR$ VECTOR ALIGNED
              DO el = elblk(1,blk,et),elblk(2,blk,et)
!                 H(el,dof) = Hold(el,dof) + dt*MirhsH(el,dof)
                Z(el,dof) = Zold(el,dof) + dt*MirhsZ(el,dof)
                MirhsZ(el,dof) = 0d0
              ENDDO
!DIR$ IVDEP
!!DIR$ VECTOR ALIGNED
              DO el = elblk(1,blk,et),elblk(2,blk,et)
                Qx(el,dof) = Qxold(el,dof) + dt*MirhsQx(el,dof)
                MirhsQx(el,dof) = 0d0
              ENDDO
!DIR$ IVDEP              
!!DIR$ VECTOR ALIGNED
              DO el = elblk(1,blk,et),elblk(2,blk,et)
                Qy(el,dof) = Qyold(el,dof) + dt*MirhsQy(el,dof)
                MirhsQy(el,dof) = 0d0
              ENDDO
            ENDDO
            
          ENDIF
        ENDDO
      ENDDO
           
      
      END SUBROUTINE forward_euler

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
      
      
      SUBROUTINE nan_check()
      
      USE globals, ONLY: ndof,ne, &
                         H,Z,Qx,Qy, &
                         nel_type,elblk,npartet, &
                         t
                         
      USE quit, ONLY: finish 
      USE output, ONLY: output_solution,close_output
      USE read_dginp, ONLY: npart     
      USE messenger2, ONLY: myrank,lel2gel
      
#ifdef CMPI       
      USE mpi                            
#endif        
      
      IMPLICIT NONE
      
      INTEGER :: et,dof,el,blk
      INTEGER :: nan_error,nan_error_sum,cnt,ierr
      
      nan_error = 0
      nan_error_sum = 0
      
      DO blk = 1,npart+1
        DO et = 1,nel_type
          IF (npartet(et,blk) > 0) THEN
      
            DO dof = 1,ndof(et)
      
              DO el = elblk(1,blk,et),elblk(2,blk,et)
                IF (Z(el,dof) /= Z(el,dof)) THEN
                  nan_error = 1
                  PRINT "(2(A,I9))", "NaN detected in Z solution for (global) element : ",lel2gel(el), " on PE: ", myrank
                  PRINT("(A,e15.8)"), 't = ', t                      
                ENDIF
              ENDDO
        
              DO el = elblk(1,blk,et),elblk(2,blk,et)
                IF (Qx(el,dof) /= Qx(el,dof)) THEN
                  nan_error = 1
                  PRINT "(2(A,I9))", "NaN detected in Qx solution for (global) element : ",lel2gel(el), " on PE: ", myrank
                  PRINT("(A,e15.8)"), 't = ', t
                ENDIF
              ENDDO
        
              DO el = elblk(1,blk,et),elblk(2,blk,et)
                IF (Qy(el,dof) /= Qy(el,dof)) THEN
                  nan_error = 1
                  PRINT "(2(A,I9))", "NaN detected in Qy solution for (global) element : ",lel2gel(el), " on PE: ", myrank
                  PRINT("(A,e15.8)"), 't = ', t
                ENDIF
              ENDDO
        
            ENDDO
      
          ENDIF
        ENDDO
      ENDDO
      
#ifdef CMPI
      cnt = 1
      CALL MPI_ALLREDUCE(nan_error,nan_error_sum,cnt,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)      
#else          
      nan_error_sum = nan_error    
#endif
      
      IF (nan_error_sum > 0) THEN
        CALL output_solution(.false.)
        CALL close_output()                  
        CALL finish(myrank)      
      ENDIF
      
      RETURN
      END SUBROUTINE nan_check
