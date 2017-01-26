      SUBROUTINE rhs2()

      USE globals,  ONLY: l,ind,el,ne,pt,led,ed,dof,ndof,g,tstage,ramp,arg,bfr, & 
                          H,Qx,Qy, &
                          rhsH,rhsQx,rhsQy, &
                          Hqpt,Qxqpt,Qyqpt, &
                          nqpta,wpta,phia,dpdx,dpdy,nqpte,wpte,phie, &
                          press,recipH,xmom,ymom,xymom,cf,tau,u,v, &
                          area,edlen,normal,dhbdx,dhbdy, &
                          el_in,el_ex,led_in,led_ex,gp_in,gp_ex, &
                          nx,ny,nx2,ny2,nxny,tx,ty, &
                          H_in,H_ex,Qx_in,Qx_ex,Qy_in,Qy_ex,Qn,Qt, &
                          xmom_in,xmom_ex,ymom_in,ymom_ex,xymom_in,xymom_ex, &
                          Hhat,Qxhat,Qyhat, &
                          Hflux,Qxflux,Qyflux, &
                          nied,iedn,nobed,obedn,nnfbed,nfbedn,nfbed,fbedn, &
                          nobfr,obfreq,obper,obnfact,obamp_qpt,obph_qpt,obeq,obdepth_qpt, &
                          nfbfr,fbfreq,fbper,fbnfact,fbamp_qpt,fbph_qpt,fbeq, &
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

            DO el = 1,ne
              u = Qxqpt(el,1)/Hqpt(el,1)
              v = Qyqpt(el,1)/Hqpt(el,1)
              tau(el) = cf*sqrt(u*u + v*v)/Hqpt(el,1)
            ENDDO

      test: DO l = 1,ndof

              ind = (l-1)*nqpta+pt
              DO el = 1,ne
                rhsH(el,l) = rhsH(el,l) + wpta(pt)*(Qxqpt(el,1)*dpdx(el,ind) + Qyqpt(el,1)*dpdy(el,ind))
              ENDDO
              DO el = 1,ne
                rhsQx(el,l) = rhsQx(el,l) + wpta(pt)*(xmom(el,1)*dpdx(el,ind) + xymom(el,1)*dpdy(el,ind) + area(el)*( g*Hqpt(el,1)*dhbdx(el) - tau(el)*Qxqpt(el,1) )*phia(l,pt))

                rhsQy(el,l) = rhsQy(el,l) + wpta(pt)*(xymom(el,1)*dpdx(el,ind) + ymom(el,1)*dpdy(el,ind) + area(el)*( g*Hqpt(el,1)*dhbdy(el) - tau(el)*Qyqpt(el,1) )*phia(l,pt))
              ENDDO

            ENDDO test 

          ENDDO a_points

      
!     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     c Edge Integrals
!     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
     

          Hqpt(:,:) = 0d0
          Qxqpt(:,:) = 0d0
          Qyqpt(:,:) = 0d0

ed_points: DO pt = 1,nqpte

    l_edge: DO led = 1,3
              ind = (led-1)*nqpte + pt
    ed_basis: DO dof = 1,ndof
                 
                DO el = 1,ne
                  Hqpt(el,ind) = Hqpt(el,ind) + H(el,dof)*phie(dof,ind)            
                ENDDO
                DO el = 1,ne
                  Qxqpt(el,ind) = Qxqpt(el,ind) + Qx(el,dof)*phie(dof,ind)            
                ENDDO
                DO el = 1,ne
                  Qyqpt(el,ind) = Qyqpt(el,ind) + Qy(el,dof)*phie(dof,ind)            
                ENDDO

              ENDDO ed_basis

              DO el = 1,ne
                press = .5d0*g*Hqpt(el,ind)*Hqpt(el,ind)
                recipH = 1d0/Hqpt(el,ind)
 
                xmom(el,ind) = Qxqpt(el,ind)*Qxqpt(el,ind)*recipH + press
                ymom(el,ind) = Qyqpt(el,ind)*Qyqpt(el,ind)*recipH + press
                xymom(el,ind) = Qxqpt(el,ind)*Qyqpt(el,ind)*recipH
              ENDDO

            ENDDO l_edge

          ENDDO ed_points

ed_points2: DO pt = 1,nqpte

            ! Interior edges
            DO ed = 1,nied

              ged = iedn(ed)

              led_in = ged2led(1,ged)
              led_ex = ged2led(2,ged)

              gp_in = (led_in-1)*nqpte + pt
              gp_ex = (led_ex-1)*nqpte + nqpte - pt + 1

              el_in = ged2el(1,ged)
              el_ex = ged2el(2,ged)

              H_in = Hqpt(el_in,gp_in)
              H_ex = Hqpt(el_ex,gp_ex)

              Qx_in = Qxqpt(el_in,gp_in)
              Qx_ex = Qxqpt(el_ex,gp_ex)

              Qy_in = Qyqpt(el_in,gp_in)
              Qy_ex = Qyqpt(el_ex,gp_ex)

              nx = normal(1,ged)
              ny = normal(2,ged)

              xmom_in = xmom(el_in,gp_in)
              xmom_ex = xmom(el_ex,gp_ex)

              ymom_in = ymom(el_in,gp_in)
              ymom_ex = ymom(el_ex,gp_ex)

              xymom_in = xymom(el_in,gp_in)
              xymom_ex = xymom(el_ex,gp_ex)

              CALL numerical_flux()

              Hflux(el_in,gp_in) = edlen(ged)*Hhat
              Hflux(el_ex,gp_ex) = -Hflux(el_in,gp_in)

              Qxflux(el_in,gp_in) = edlen(ged)*Qxhat
              Qxflux(el_ex,gp_ex) = -Qxflux(el_in,gp_in)

              Qyflux(el_in,gp_in) = edlen(ged)*Qyhat
              Qyflux(el_ex,gp_ex) = -Qyflux(el_in,gp_in) 

            ENDDO
 
            ! No normal flow boundary condition 
            DO ed = 1,nnfbed

              ged = nfbedn(ed)

              led_in = ged2led(1,ged)

              gp_in = (led_in-1)*nqpte + pt

              el_in = ged2el(1,ged)

              nx = normal(1,ged)
              ny = normal(2,ged)

              nx2 = nx*nx
              ny2 = ny*ny
              nxny = nx*ny

              H_in = Hqpt(el_in,gp_in)
              H_ex = H_in

              Qx_in = Qxqpt(el_in,gp_in)
              Qy_in = Qyqpt(el_in,gp_in)

              Qx_ex = Qx_in*(ny2-nx2) - 2d0*nxny*Qy_in
              Qy_ex = Qy_in*(nx2-ny2) - 2d0*nxny*Qx_in

              press = .5d0*g*H_ex*H_ex
              recipH = 1d0/H_ex

              xmom_in = xmom(el_in,gp_in)
              xmom_ex = (Qx_ex*Qx_ex)*recipH + press

              ymom_in = ymom(el_in,gp_in)
              ymom_ex = (Qy_ex*Qy_ex)*recipH + press

              xymom_in = xymom(el_in,gp_in)
              xymom_ex = (Qx_ex*Qy_ex)*recipH

              CALL numerical_flux()

              Hflux(el_in,gp_in) = edlen(ged)*Hhat

              Qxflux(el_in,gp_in) = edlen(ged)*Qxhat

              Qyflux(el_in,gp_in) = edlen(ged)*Qyhat
   
            ENDDO

            ! Flow specified boundary edges
            DO ed = 1,nfbed

              ged = fbedn(ed)

              led_in = ged2led(1,ged)

              gp_in = (led_in-1)*nqpte + pt

              el_in = ged2el(1,ged)

              nx = normal(1,ged)
              ny = normal(2,ged)

              tx = -ny
              ty = nx

              H_in = Hqpt(el_in,gp_in)
              H_ex = H_in

              Qx_in = Qxqpt(el_in,gp_in)
              Qy_in = Qyqpt(el_in,gp_in)

              Qn = 0d0
              Qt= 0d0
              DO bfr = 1,nfbfr
                ind = (pt-1)*nfbfr + bfr
                arg = fbfreq(bfr)*(tstage - INT(tstage/fbper(bfr))*fbper(bfr)) + fbeq(bfr)
                Qn = Qn + fbamp_qpt(ind,ed)*fbnfact(bfr)*ramp*COS(arg-fbph_qpt(ind,ed))
              ENDDO

!               Qn = 5d0

              Qn = -Qn

              Qx_ex = ( ty*Qn - ny*Qt)/(nx*ty-ny*tx)
              Qy_ex = (-tx*Qn + nx*Qt)/(nx*ty-ny*tx)
 
              press = .5d0*g*H_ex*H_ex
              recipH = 1d0/H_ex

              xmom_in = xmom(el_in,gp_in)
              xmom_ex = (Qx_ex*Qx_ex)*recipH + press

              ymom_in = ymom(el_in,gp_in)
              ymom_ex = (Qy_ex*Qy_ex)*recipH + press

              xymom_in = xymom(el_in,gp_in)
              xymom_ex = (Qx_ex*Qy_ex)*recipH

              CALL numerical_flux()

              Hflux(el_in,gp_in) = edlen(ged)*Hhat

              Qxflux(el_in,gp_in) = edlen(ged)*Qxhat

              Qyflux(el_in,gp_in) = edlen(ged)*Qyhat

            ENDDO

             ! Open boundary edges (elevation specified)
            DO ed = 1,nobed

              ged = obedn(ed)

              led_in = ged2led(1,ged)

              gp_in = (led_in-1)*nqpte + pt

              el_in = ged2el(1,ged)

              nx = normal(1,ged)
              ny = normal(2,ged)

              H_in = Hqpt(el_in,gp_in)

              H_ex = 0d0
              DO bfr = 1,nobfr
                ind = (pt-1)*nobfr + bfr
                arg = obfreq(bfr)*(tstage-INT(tstage/obper(bfr))*obper(bfr)) + obeq(bfr)
                H_ex = H_ex + obamp_qpt(ind,ed)*obnfact(bfr)*ramp*COS(arg-obph_qpt(ind,ed))
              ENDDO
              H_ex = H_ex + obdepth_qpt(ed,pt)

!               H_ex = 10d0

              Qx_in = Qxqpt(el_in,gp_in)
              Qy_in = Qyqpt(el_in,gp_in)

              Qx_ex = Qx_in
              Qy_ex = Qy_in
 
              press = .5d0*g*H_ex*H_ex
              recipH = 1d0/H_ex

              xmom_in = xmom(el_in,gp_in)
              xmom_ex = (Qx_ex*Qx_ex)*recipH + press

              ymom_in = ymom(el_in,gp_in)
              ymom_ex = (Qy_ex*Qy_ex)*recipH + press

              xymom_in = xymom(el_in,gp_in)
              xymom_ex = (Qx_ex*Qy_ex)*recipH

              CALL numerical_flux()

              Hflux(el_in,gp_in) = edlen(ged)*Hhat

              Qxflux(el_in,gp_in) = edlen(ged)*Qxhat

              Qyflux(el_in,gp_in) = edlen(ged)*Qyhat                                 

             ENDDO

          ENDDO ed_points2

             
           DO pt = 1,nqpte
             DO led = 1,3
               ind = (led-1)*nqpte+pt
               DO l = 1,ndof
                 DO el = 1,ne
                   rhsH(el,l) = rhsH(el,l) - wpte(pt)*Hflux(el,ind)*phie(l,ind)
                   rhsQx(el,l) = rhsQx(el,l) - wpte(pt)*Qxflux(el,ind)*phie(l,ind)
                   rhsQy(el,l) = rhsQy(el,l) - wpte(pt)*Qyflux(el,ind)*phie(l,ind)
                 ENDDO
               ENDDO                                    
             ENDDO
           ENDDO

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
      END SUBROUTINE rhs2