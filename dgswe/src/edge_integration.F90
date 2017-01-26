       SUBROUTINE interior_edge_eval(et,sel,eel,ndof,tnqpte)

       USE globals, ONLY: Zqpt,Hqpt,Qxqpt,Qyqpt, &
                          Z,Qx,Qy, &
                          phie,hbqpte, &
                          g,pt5g,recipHa,xmom,ymom,xymom
                                               

       IMPLICIT NONE
       
       INTEGER :: pt,el,dof,et
       INTEGER :: tnqpte,ndof
       INTEGER :: sel,eel
       
ed_points: DO pt = 1,tnqpte

!!DIR$ VECTOR ALIGNED     
!DIR$ SIMD
              DO el = sel,eel
                Zqpt(el,pt)  = Z(el,1)
                Qxqpt(el,pt) = Qx(el,1)
                Qyqpt(el,pt) = Qy(el,1)
              ENDDO

    ed_basis: DO dof = 2,ndof     
!!DIR$ VECTOR ALIGNED
!DIR$ SIMD
                DO el = sel,eel   ! Compute solutions at edge quadrature points                
                  Zqpt(el,pt)  = Zqpt(el,pt)  + Z(el,dof)*phie(dof,pt,et)            
                  Qxqpt(el,pt) = Qxqpt(el,pt) + Qx(el,dof)*phie(dof,pt,et)            
                  Qyqpt(el,pt) = Qyqpt(el,pt) + Qy(el,dof)*phie(dof,pt,et)                              
                ENDDO

              ENDDO ed_basis

!!DIR$ VECTOR ALIGNED
!DIR$ SIMD
              DO el = sel,eel  ! Compute momentum terms   
                Hqpt(el,pt) = Zqpt(el,pt) + hbqpte(el,pt)    
                recipHa(el) = 1d0/Hqpt(el,pt)
                
                xmom(el,pt) = pt5g*(Hqpt(el,pt)*Hqpt(el,pt) - hbqpte(el,pt)*hbqpte(el,pt))+ Qxqpt(el,pt)*Qxqpt(el,pt)*recipHa(el) 
                ymom(el,pt) = pt5g*(Hqpt(el,pt)*Hqpt(el,pt) - hbqpte(el,pt)*hbqpte(el,pt))+ Qyqpt(el,pt)*Qyqpt(el,pt)*recipHa(el) 
                xymom(el,pt) = Qxqpt(el,pt)*Qyqpt(el,pt)*recipHa(el)
              ENDDO

          ENDDO ed_points

      END SUBROUTINE interior_edge_eval
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      SUBROUTINE interior_edge_eval_ldg_Q(et,sel,eel,ndof,tnqpte)
      
      USE globals, ONLY: Qxqpt,Qyqpt, &
                         Qx,Qy, &
                         phie
      
      IMPLICIT NONE
      
      INTEGER :: et
      INTEGER :: tnqpte,ndof
      INTEGER :: sel,eel
       
      INTEGER :: pt,el,dof
      
      DO pt = 1,tnqpte
        DO el = sel,eel
          Qxqpt(el,pt) = Qx(el,1)
          Qyqpt(el,pt) = Qy(el,1)          
        ENDDO
        
        DO dof = 2,ndof
          DO el = sel,eel
            Qxqpt(el,pt) = Qxqpt(el,pt) + Qx(el,dof)*phie(dof,pt,et)            
            Qyqpt(el,pt) = Qyqpt(el,pt) + Qy(el,dof)*phie(dof,pt,et)            
          ENDDO
        ENDDO
      ENDDO
            
      RETURN
      END SUBROUTINE interior_edge_eval_ldg_Q


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     


      SUBROUTINE interior_edge_eval_ldg_E(et,sel,eel,ndof,tnqpte)
      
      USE globals, ONLY: Exxqpt,Eyyqpt,Exyqpt,Eyxqpt, &
                         Exx,Eyy,Exy,Eyx, &
                         phie
      
      IMPLICIT NONE
      
      INTEGER :: et
      INTEGER :: tnqpte,ndof
      INTEGER :: sel,eel
       
      INTEGER :: pt,el,dof
      
      DO pt = 1,tnqpte
        DO el = sel,eel
          Exxqpt(el,pt) = Exx(el,1)
          Eyyqpt(el,pt) = Eyy(el,1)          
          Exyqpt(el,pt) = Exy(el,1)           
!           Eyxqpt(el,pt) = Eyx(el,1)          
        ENDDO
        
        DO dof = 2,ndof
          DO el = sel,eel
            Exxqpt(el,pt) = Exxqpt(el,pt) + Exx(el,dof)*phie(dof,pt,et)            
            Eyyqpt(el,pt) = Eyyqpt(el,pt) + Eyy(el,dof)*phie(dof,pt,et)  
            Exyqpt(el,pt) = Exyqpt(el,pt) + Exy(el,dof)*phie(dof,pt,et)
!             Eyxqpt(el,pt) = Eyxqpt(el,pt) + Eyx(el,dof)*phie(dof,pt,et)            
          ENDDO
        ENDDO
      ENDDO
            
      RETURN
      END SUBROUTINE interior_edge_eval_ldg_E
      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
       
       SUBROUTINE edge_integration(et,sel,eel,ndof,tnqpte)
       
       USE globals, ONLY: rhsZ,rhsQx,rhsQy, &
                          Zqpt,Qxqpt,Qyqpt, &
                          phie_int
       
       IMPLICIT NONE
       
       INTEGER :: pt,l,el
       INTEGER :: et
       INTEGER :: sel,eel
       INTEGER :: ndof,tnqpte
       
       DO pt = 1,tnqpte
         DO l = 1,ndof
!!DIR$ VECTOR ALIGNED
           DO el = sel,eel
             rhsZ(el,l) = rhsZ(el,l) - Zqpt(el,pt)*phie_int(l,pt,et)
             rhsQx(el,l) = rhsQx(el,l) - Qxqpt(el,pt)*phie_int(l,pt,et)
             rhsQy(el,l) = rhsQy(el,l) - Qyqpt(el,pt)*phie_int(l,pt,et)                   
           ENDDO
         ENDDO                                    
       ENDDO       
       
       RETURN
       END SUBROUTINE edge_integration            

       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      SUBROUTINE edge_integration_ldg(et,sel,eel,ndof,tnqpte)
      
      USE globals, ONLY: rhsExx,rhsEyy,rhsExy,rhsEyx, &
                         Exxqpt,Eyyqpt,Exyqpt,Eyxqpt, &
                         phie_int
                         
      USE read_dginp, ONLY: esl                          
      
      IMPLICIT NONE
      
      INTEGER :: et
      INTEGER :: sel,eel
      INTEGER :: ndof,tnqpte      
      
      INTEGER :: pt,l,el
      
      DO pt = 1,tnqpte
        DO l = 1,ndof
          DO el = sel,eel
            rhsExx(el,l) = rhsExx(el,l) + Exxqpt(el,pt)*phie_int(l,pt,et)
            rhsEyy(el,l) = rhsEyy(el,l) + Eyyqpt(el,pt)*phie_int(l,pt,et)
            rhsExy(el,l) = rhsExy(el,l) + Exyqpt(el,pt)*phie_int(l,pt,et)
!             rhsEyx(el,l) = rhsEyx(el,l) + Eyxqpt(el,pt)*phie_int(l,pt,et)            
          ENDDO
        ENDDO
      ENDDO
      
      DO l = 1,ndof
        DO el = sel,eel
          rhsExx(el,l) = 2d0*esl*rhsExx(el,l)
          rhsEyy(el,l) = 2d0*esl*rhsEyy(el,l)
          rhsExy(el,l) = esl*rhsExy(el,l)

!           rhsExx(el,l) = esl*rhsExx(el,l)
!           rhsEyy(el,l) = esl*rhsEyy(el,l)
!           rhsExy(el,l) = esl*rhsExy(el,l)
!           rhsEyx(el,l) = esl*rhsEyx(el,l)
        ENDDO
      ENDDO
      
      RETURN
      END SUBROUTINE edge_integration_ldg
      
      


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       