      SUBROUTINE rhs()

      USE globals,  ONLY: l,ind,el,ne,pt,led,ied,nqpta,nqpte,dof,ndof,g, & 
                          H,Qx,Qy, &
                          rhsH,rhsQx,rhsQy, &
                          Hqpt,Qxqpt,Qyqpt, &
                          wpta,phia,dpdx,dpdy,wpte,phie, &
                          press,recipH,xmom,ymom,xymom,tau, &
                          area,edlen,normal, &
                          el_in,el_ex,led_in,led_ex, &
                          eig_in,eig_ex,alpha,nx,ny,nx2,ny2,nxny, &
                          H_in,H_ex,Qx_in,Qx_ex,Qy_in,Qy_ex, &
                          xmom_in,xmom_ex,ymom_in,ymom_ex,xymom_in,xymom_ex, &
                          Hhat,Qxhat,Qyhat, &
                          Hflux,Qxflux,Qyflux, &
                          nied,iedn,nbed,bedn, &
                          ged2led,ged2el,ged
                   
      IMPLICIT NONE

!     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     c Area Integrals
!     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  

      rhsH(:,:) = 0d0
      rhsQx(:,:) = 0d0
      rhsQy(:,:) = 0d0    

a_points: DO pt = 1,nqpta

            Hqpt(:,1) = 0d0
            Qxqpt(:,1) = 0d0
            Qyqpt(:,1) = 0d0

   a_basis: DO dof = 1,ndof

              DO el = 1,ne
                Hqpt(el,1) = Hqpt(el,1) + H(el,dof)*phia(dof,pt)
              ENDDO
              DO el = 1,ne
                Qxqpt(el,1) = Qxqpt(el,1) + Qx(el,dof)*phia(dof,pt) 
              ENDDO
              DO el = 1,ne
                Qyqpt(el,1) = Qyqpt(el,1) + Qy(el,dof)*phia(dof,pt)
              ENDDO

            ENDDO a_basis

            DO el = 1,ne
              press = .5d0*g*Hqpt(el,1)*Hqpt(el,1)
              recipH = 1d0/Hqpt(el,1)

              xmom(el,1) = Qxqpt(el,1)*Qxqpt(el,1)*recipH + press
              ymom(el,1) = Qyqpt(el,1)*Qyqpt(el,1)*recipH + press
              xymom(el,1) = Qxqpt(el,1)*Qyqpt(el,1)*recipH
            ENDDO

      test: DO l = 1,ndof

              ind = (l-1)*nqpta+pt
              DO el = 1,ne
                rhsH(el,l) = rhsH(el,l) + wpta(pt)*(Qxqpt(el,1)*dpdx(el,ind) + Qyqpt(el,1)*dpdy(el,ind))
              ENDDO
              DO el = 1,ne
                rhsQx(el,l) = rhsQx(el,l) + wpta(pt)*(xmom(el,1)*dpdx(el,ind) + xymom(el,1)*dpdy(el,ind) - area(el)*tau*Qxqpt(el,1)*phia(l,pt))
                rhsQy(el,l) = rhsQy(el,l) + wpta(pt)*(xymom(el,1)*dpdx(el,ind) + ymom(el,1)*dpdy(el,ind) - area(el)*tau*Qyqpt(el,1)*phia(l,pt))
              ENDDO

            ENDDO test 

          ENDDO a_points

      
!     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     c Edge Integrals
!     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
     

ed_points: DO pt = 1,nqpte

            Hqpt(:,:) = 0d0
            Qxqpt(:,:) = 0d0
            Qyqpt(:,:) = 0d0

    l_edge: DO led = 1,3
              ind = (led-1)*nqpte + pt
    ed_basis: DO dof = 1,ndof
                 
                DO el = 1,ne
                  Hqpt(el,led) = Hqpt(el,led) + H(el,dof)*phie(dof,ind)            
                ENDDO
                DO el = 1,ne
                  Qxqpt(el,led) = Qxqpt(el,led) + Qx(el,dof)*phie(dof,ind)            
                ENDDO
                DO el = 1,ne
                  Qyqpt(el,led) = Qyqpt(el,led) + Qy(el,dof)*phie(dof,ind)            
                ENDDO

              ENDDO ed_basis

              DO el = 1,ne
                press = .5d0*g*Hqpt(el,led)*Hqpt(el,led)
                recipH = 1d0/Hqpt(el,led)
 
                xmom(el,led) = Qxqpt(el,led)*Qxqpt(el,led)*recipH + press
                ymom(el,led) = Qyqpt(el,led)*Qyqpt(el,led)*recipH + press
                xymom(el,led) = Qxqpt(el,led)*Qyqpt(el,led)*recipH
              ENDDO

            ENDDO l_edge

            ! Interior edges
           DO ied = 1,nied

              ged = iedn(ied)

              led_in = ged2led(1,ged)
              led_ex = ged2led(2,ged)

              el_in = ged2el(1,ged)
              el_ex = ged2el(2,ged)

              H_in = Hqpt(el_in,led_in)
              H_ex = Hqpt(el_ex,led_ex)

              Qx_in = Qxqpt(el_in,led_in)
              Qx_ex = Qxqpt(el_ex,led_ex)

              Qy_in = Qyqpt(el_in,led_in)
              Qy_ex = Qyqpt(el_ex,led_ex)

              nx = normal(1,ged)
              ny = normal(2,ged)

              eig_in = abs(Qx_in*nx + Qy_in*ny)/H_in+sqrt(g*H_in) 
              eig_ex = abs(Qx_ex*nx + Qy_ex*ny)/H_ex+sqrt(g*H_ex)

              alpha = max(eig_in,eig_ex)

              xmom_in = xmom(el_in,led_in)
              xmom_ex = xmom(el_ex,led_ex)

              ymom_in = ymom(el_in,led_in)
              ymom_ex = ymom(el_ex,led_ex)

              xymom_in = xymom(el_in,led_in)
              xymom_ex = xymom(el_ex,led_ex)

              CALL numerical_flux()

              Hflux(el_in,led_in) = edlen(ged)*Hhat
              Hflux(el_ex,led_ex) = -Hflux(el_in,led_in)

              Qxflux(el_in,led_in) = edlen(ged)*Qxhat
              Qxflux(el_ex,led_ex) = -Qxflux(el_in,led_in)

              Qyflux(el_in,led_in) = edlen(ged)*Qyhat
              Qyflux(el_ex,led_ex) = -Qyflux(el_in,led_in)

            ENDDO
 
            ! No normal flow boundary condition 
            DO ied = 1,nbed

              ged = bedn(ied)

              led_in = ged2led(1,ged)

              el_in = ged2el(1,ged)

              nx = normal(1,ged)
              ny = normal(2,ged)

              nx2 = nx*nx
              ny2 = ny*ny
              nxny = nx*ny

              H_in = Hqpt(el_in,led_in)
              H_ex = H_in

              Qx_in = Qxqpt(el_in,led_in)
              Qy_in = Qyqpt(el_in,led_in)

              Qx_ex = Qx_in*(ny2-nx2) - 2d0*nxny*Qy_in
              Qy_ex = Qy_in*(nx2-ny2) - 2d0*nxny*Qx_in

              alpha = abs(Qx_in*nx + Qy_in*ny)/H_in+sqrt(g*H_in) 
 
              press = .5d0*g*H_ex*H_ex
              recipH = 1d0/H_ex

              xmom_in = xmom(el_in,led_in)
              xmom_ex = (Qx_ex*Qx_ex)*recipH + press

              ymom_in = ymom(el_in,led_in)
              ymom_ex = (Qy_ex*Qy_ex)*recipH + press

              xymom_in = xymom(el_in,led_in)
              xymom_ex = (Qx_ex*Qy_ex)*recipH

              CALL numerical_flux()

              Hflux(el_in,led_in) = edlen(ged)*Hhat

              Qxflux(el_in,led_in) = edlen(ged)*Qxhat

              Qyflux(el_in,led_in) = edlen(ged)*Qyhat
  
             ENDDO

             
             DO led = 1,3
               ind = (led-1)*nqpte+pt
               DO l = 1,ndof
                 DO el = 1,ne
                   rhsH(el,l) = rhsH(el,l) - wpte(pt)*Hflux(el,led)*phie(l,ind)
                   rhsQx(el,l) = rhsQx(el,l) - wpte(pt)*Qxflux(el,led)*phie(l,ind)
                   rhsQy(el,l) = rhsQy(el,l) - wpte(pt)*Qyflux(el,led)*phie(l,ind)
                 ENDDO
               ENDDO                                    
             ENDDO


          ENDDO ed_points

          ! Multiply by inverse of mass matrix
          DO dof = 1,ndof
            DO el = 1,ne
              rhsH(el,dof) = rhsH(el,dof)/area(el)
            ENDDO

            DO el = 1,ne
              rhsQx(el,dof) = rhsQx(el,dof)/area(el)
            ENDDO

            DO el = 1,ne
              rhsQy(el,dof) = rhsQy(el,dof)/area(el)
            ENDDO
          ENDDO


      RETURN
      END SUBROUTINE rhs