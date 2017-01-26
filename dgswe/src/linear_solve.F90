      SUBROUTINE linear_solve(et,sel,eel,ndof)

      USE globals, ONLY: rhsZ,rhsQx,rhsQy, &
                         MirhsZ,MirhsQx,MirhsQy, &
                         mmi

      IMPLICIT NONE

      INTEGER :: et,sel,eel,ndof
      INTEGER :: i,j,l,el,m

            SELECT CASE(et)
               
              CASE(1)
                DO l = 1,ndof
                  DO el = sel,eel
                    MirhsZ(el,l)  = MirhsZ(el,l)  + mmi(el,1)*rhsZ(el,l) 
                    MirhsQx(el,l) = MirhsQx(el,l) + mmi(el,1)*rhsQx(el,l) 
                    MirhsQy(el,l) = MirhsQy(el,l) + mmi(el,1)*rhsQy(el,l) 
                    
                    rhsZ(el,l) = 0d0
                    rhsQx(el,l) = 0d0
                    rhsQy(el,l) = 0d0
                  ENDDO
                ENDDO
 
              CASE DEFAULT
              
                m = 1
                DO i = 1,ndof
                  DO j = 1,ndof
                    DO el = sel,eel
                      MirhsZ(el,i)  = MirhsZ(el,i)  + mmi(el,m)*rhsZ(el,j) 
                      MirhsQx(el,i) = MirhsQx(el,i) + mmi(el,m)*rhsQx(el,j) 
                      MirhsQy(el,i) = MirhsQy(el,i) + mmi(el,m)*rhsQy(el,j)                                       
                    ENDDO
                    m = m + 1
                  ENDDO
                ENDDO 
                
                DO i = 1,ndof
                  DO el = sel,eel
                    rhsZ(el,i) = 0d0
                    rhsQx(el,i) = 0d0
                    rhsQy(el,i) = 0d0     
                  ENDDO
                ENDDO
                
            END SELECT      

      END SUBROUTINE linear_solve
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      


      SUBROUTINE linear_solve_ldg(et,sel,eel,ndof)
      
      USE globals, ONLY: Exx,Eyy,Exy,Eyx, &
                         rhsExx,rhsEyy,rhsExy,rhsEyx, &
                         mmi
      
      IMPLICIT NONE
      
      INTEGER :: et
      INTEGER :: sel,eel
      INTEGER :: ndof
      
      INTEGER :: l,el,i,j
      INTEGER :: m
      
      DO l = 1,ndof
        DO el = sel,eel
          Exx(el,l) = 0d0
          Eyy(el,l) = 0d0
          Exy(el,l) = 0d0
!           Eyx(el,l) = 0d0
        ENDDO
      ENDDO
            
      
      SELECT CASE(et)
      
        CASE(1)
        
          DO l = 1,ndof
            DO el = sel,eel
              Exx(el,l) = Exx(el,l) + mmi(el,1)*rhsExx(el,l)
              Eyy(el,l) = Eyy(el,l) + mmi(el,1)*rhsEyy(el,l)
              Exy(el,l) = Exy(el,l) + mmi(el,1)*rhsExy(el,l)
!               Eyx(el,l) = Eyx(el,l) + mmi(el,1)*rhsEyx(el,l)              
              
              rhsExx(el,l) = 0d0
              rhsEyy(el,l) = 0d0
              rhsExy(el,l) = 0d0
!               rhsEyx(el,l) = 0d0              
            ENDDO
          ENDDO
        
        CASE DEFAULT
        
          m = 1
          DO i = 1,ndof
            DO j = 1,ndof
              DO el = sel,eel
                Exx(el,i) = Exx(el,i) + mmi(el,m)*rhsExx(el,j)
                Eyy(el,i) = Eyy(el,i) + mmi(el,m)*rhsEyy(el,j)
                Exy(el,i) = Exy(el,i) + mmi(el,m)*rhsExy(el,j)
!                 Eyx(el,i) = Eyx(el,i) + mmi(el,m)*rhsEyx(el,j)                
              ENDDO
              m = m + 1
            ENDDO
          ENDDO
          
          DO i = 1,ndof
            DO el = sel,eel                
              rhsExx(el,i) = 0d0
              rhsEyy(el,i) = 0d0
              rhsExy(el,i) = 0d0
!               rhsEyx(el,i) = 0d0                
            ENDDO
          ENDDO 
          
      END SELECT
      
      RETURN
      END SUBROUTINE linear_solve_ldg