      SUBROUTINE area_integration(et,sel,eel,ndof,nqpta)
      
      USE globals, ONLY: Zqpta,Hqpta,Qxqpta,Qyqpta, &
                         Z,Qx,Qy, &
                         phia,phia_int,dpdx,dpdy, &
                         recipHa,xmoma,ymoma,xymoma, &
                         tau,src_x,src_y, &
                         dhbdx,dhbdy, &
                         rhsZ,rhsQx,rhsQy,&
                         g,pt5g,hbqpta, &
                         Exx,Eyy,Exy,Eyx, &
                         Exxqpta,Eyyqpta,Exyqpta,Eyxqpta
                         
      USE read_dginp, ONLY: cf                         

      IMPLICIT NONE

      INTEGER :: nqpta,ndof,et
      INTEGER :: sel,eel

      INTEGER :: pt,dof,el,l


a_points: DO pt = 1,nqpta

!!DIR$ VECTOR ALIGNED
!DIR$ SIMD
            DO el = sel,eel  ! First basis function is 1
              Zqpta(el)  = Z(el,1)
              Qxqpta(el) = Qx(el,1)
              Qyqpta(el) = Qy(el,1)
              
              Exxqpta(el) = Exx(el,1)
              Eyyqpta(el) = Eyy(el,1) 
              Exyqpta(el) = Exy(el,1)
!               Eyxqpta(el) = Eyx(el,1)              
            ENDDO

   a_basis: DO dof = 2,ndof
!!DIR$ VECTOR ALIGNED
!DIR$ SIMD
              DO el = sel,eel   ! Evaluate solution at area quadrature point
                Zqpta(el)  = Zqpta(el)  + Z(el,dof)*phia(dof,pt,et)
                Qxqpta(el) = Qxqpta(el) + Qx(el,dof)*phia(dof,pt,et) 
                Qyqpta(el) = Qyqpta(el) + Qy(el,dof)*phia(dof,pt,et)
                
                Exxqpta(el) = Exxqpta(el) + Exx(el,dof)*phia(dof,pt,et)
                Eyyqpta(el) = Eyyqpta(el) + Eyy(el,dof)*phia(dof,pt,et)
                Exyqpta(el) = Exyqpta(el) + Exy(el,dof)*phia(dof,pt,et)
!                 Eyxqpta(el) = Eyxqpta(el) + Eyx(el,dof)*phia(dof,pt,et)                
              ENDDO

            ENDDO a_basis

!!DIR$ VECTOR ALIGNED
!DIR$ SIMD
            DO el = sel,eel   ! Compute momentum terms
              Hqpta(el) = Zqpta(el) + hbqpta(el,pt)
              recipHa(el) = 1d0/Hqpta(el)

              xmoma(el) = pt5g*(Hqpta(el)*Hqpta(el)-hbqpta(el,pt)*hbqpta(el,pt)) + Qxqpta(el)*Qxqpta(el)*recipHa(el) - Exxqpta(el)
              ymoma(el) = pt5g*(Hqpta(el)*Hqpta(el)-hbqpta(el,pt)*hbqpta(el,pt)) + Qyqpta(el)*Qyqpta(el)*recipHa(el) - Eyyqpta(el)
              xymoma(el) = Qxqpta(el)*Qyqpta(el)*recipHa(el) - Exyqpta(el)             
!               xymoma(el) = Qxqpta(el)*Qyqpta(el)*recipHa(el) 
            ENDDO 

!!DIR$ VECTOR ALIGNED
!DIR$ SIMD
            DO el =  sel,eel   ! Compute source terms
              tau(el) = cf*sqrt((Qxqpta(el)*recipHa(el))**2 + (Qyqpta(el)*recipHa(el))**2)*recipHa(el)
              src_x(el) = g*Zqpta(el)*dhbdx(el,pt) - tau(el)*Qxqpta(el) 
              src_y(el) = g*Zqpta(el)*dhbdy(el,pt) - tau(el)*Qyqpta(el)
            ENDDO

!!DIR$ VECTOR ALIGNED
!DIR$ SIMD
            DO el = sel,eel   ! Derivatives are 0 for first dof
              rhsQx(el,1) = rhsQx(el,1) + src_x(el)*phia_int(el,1,pt)
              rhsQy(el,1) = rhsQy(el,1) + src_y(el)*phia_int(el,1,pt)
            ENDDO

      test: DO l = 2,ndof 
!!DIR$ VECTOR ALIGNED          
!DIR$ SIMD
              DO el = sel,eel
                rhsZ(el,l)  = rhsZ(el,l)  + Qxqpta(el)*dpdx(el,l,pt) + Qyqpta(el)*dpdy(el,l,pt)

                rhsQx(el,l) = rhsQx(el,l) + xmoma(el)*dpdx(el,l,pt)  + xymoma(el)*dpdy(el,l,pt) + src_x(el)*phia_int(el,l,pt)           
                rhsQy(el,l) = rhsQy(el,l) + xymoma(el)*dpdx(el,l,pt) + ymoma(el)*dpdy(el,l,pt)  + src_y(el)*phia_int(el,l,pt)                
                
!                 rhsQx(el,l) = rhsQx(el,l) + xmoma(el)*dpdx(el,l,pt)  + (xymoma(el) - Exyqpta(el))*dpdy(el,l,pt) + src_x(el)*phia_int(el,l,pt)           
!                 rhsQy(el,l) = rhsQy(el,l) + (xymoma(el) - Eyxqpta(el))*dpdx(el,l,pt) + ymoma(el)*dpdy(el,l,pt)  + src_y(el)*phia_int(el,l,pt)
              ENDDO

            ENDDO test 

          ENDDO a_points
          
      RETURN 
      END SUBROUTINE area_integration
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      SUBROUTINE area_integration_ldg(et,sel,eel,ndof,nqpta)
      
      USE globals, ONLY: Qxqpta,Qyqpta, &
                         Qx,Qy, &
                         rhsExx,rhsEyy,rhsExy,rhsEyx, &
                         phia,dpdx,dpdy
      
      IMPLICIT NONE
      
      INTEGER :: nqpta,ndof,et
      INTEGER :: sel,eel

      INTEGER :: pt,dof,el,l
      
      
      DO pt = 1,nqpta
        DO el = sel,eel
          Qxqpta(el) = Qx(el,1)
          Qyqpta(el) = Qy(el,1)        
        ENDDO
        
        DO dof = 2,ndof
          DO el = sel,eel
            Qxqpta(el) = Qxqpta(el) + Qx(el,dof)*phia(dof,pt,et) 
            Qyqpta(el) = Qyqpta(el) + Qy(el,dof)*phia(dof,pt,et)          
          ENDDO
        ENDDO 
        
        DO l = 2,ndof
          DO el = sel,eel
            rhsExx(el,l) = rhsExx(el,l) - Qxqpta(el)*dpdx(el,l,pt)
            rhsEyy(el,l) = rhsEyy(el,l) - Qyqpta(el)*dpdy(el,l,pt)
            rhsExy(el,l) = rhsExy(el,l) - (Qxqpta(el)*dpdy(el,l,pt) + Qyqpta(el)*dpdx(el,l,pt))            
!             rhsExy(el,l) = rhsExy(el,l) - Qxqpta(el)*dpdy(el,l,pt)
!             rhsEyx(el,l) = rhsEyx(el,l) - Qyqpta(el)*dpdx(el,l,pt)
          ENDDO
        ENDDO
      ENDDO
     
      
      RETURN
      END SUBROUTINE area_integration_ldg
