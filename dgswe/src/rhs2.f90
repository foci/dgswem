      SUBROUTINE rhs2()

      USE globals,  ONLY: l,ind,el,ne,pt,led,ed,dof,ndof,split,sp,nsp,split2,nsp2, &
                          g,pt5g,tstage,ramp,arg,bfr, & 
                          H,Qx,Qy, &
                          rhsH,rhsQx,rhsQy, &
                          Hqpt,Qxqpt,Qyqpt, &
                          nqpta,wpta,phia,phia_int,dpdx,dpdy,nqpte,wpte,phie,phie_int, &
                          press,recipH,xmom,ymom,xymom,cf,tau,u,v,src_x,src_y, &
                          area,edlen_area,normal,dhbdx,dhbdy, &
                          el_in,el_ex,led_in,led_ex,gp_in,gp_ex, &
                          nx,ny,nx2,ny2,nxny,tx,ty, &
                          H_in,H_ex,Qx_in,Qx_ex,Qy_in,Qy_ex,Qn,Qt, &
                          xmom_in,xmom_ex,ymom_in,ymom_ex,xymom_in,xymom_ex, &
                          Hhat,Qxhat,Qyhat, &
                          Hflux,Qxflux,Qyflux, &
                          nied,iedn,nobed,obedn,nnfbed,nfbedn,nfbed,fbedn, &
                          nobfr,obfreq,obper,obnfact,obamp_qpt,obph_qpt,obeq,obdepth_qpt, &
                          nfbfr,fbfreq,fbper,fbnfact,fbamp_qpt,fbph_qpt,fbeq, &
                          ged2led,ged2el,ged, &
                          pressa,recipHa, & 
                          Hi,He,Qxi,Qxe,Qyi,Qye, &
                          xmi,xme,ymi,yme,xymi,xyme, &
                          Hfi,Hfe,Qxfi,Qxfe,Qyfi,Qyfe, &
                          const,inx,iny,len_area_in,len_area_ex, &
                          Hhatv,Qxhatv,Qyhatv, &
                          nxv,nyv,Hn,Qxn,Qyn,egnval, &
                          Hni,Hne,Qxni,Qxne,Qyni,Qyne,EVi,EVe, &
                          Hp,Qxp,Qyp,xmp,ymp,xymp, &
                          Hin,Hex,Qxin,Qxex,Qyin,Qyex, &
                          xmin,xmex,ymin,ymex,xymin,xymex, &
                          Hai,Hae,Qxai,Qxae,Qyai,Qyae, &
                          xmai,xmae,ymai,ymae,xymai,xymae, &
                          Hfai,Hfae,Qxfai,Qxfae,Qyfai,Qyfae                          
                          
                   
      IMPLICIT NONE

!     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     c Area Integrals
!     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  

!DIR$ VECTOR ALIGNED
      rhsH(:,:) = 0d0
!DIR$ VECTOR ALIGNED
      rhsQx(:,:) = 0d0
!DIR$ VECTOR ALIGNED
      rhsQy(:,:) = 0d0   

      DO sp = 1,nsp 

a_points: DO pt = 1,nqpta

!DIR$ VECTOR ALIGNED
            DO el = split(1,sp),split(2,sp)
              Hqpt(el,1) = H(el,1)
              Qxqpt(el,1) = Qx(el,1)
              Qyqpt(el,1) = Qy(el,1)
            ENDDO

   a_basis: DO dof = 2,ndof
!DIR$ VECTOR ALIGNED
              DO el = split(1,sp),split(2,sp) ! Evaluate solution at area quadrature point
                Hqpt(el,1) = Hqpt(el,1) + H(el,dof)*phia(dof,pt)
                Qxqpt(el,1) = Qxqpt(el,1) + Qx(el,dof)*phia(dof,pt) 
                Qyqpt(el,1) = Qyqpt(el,1) + Qy(el,dof)*phia(dof,pt)
              ENDDO

            ENDDO a_basis


!DIR$ VECTOR ALIGNED
            DO el = split(1,sp),split(2,sp) ! Compute momentum terms
              recipHa(el) = 1d0/Hqpt(el,1)

              xmom(el,1) = pt5g*Hqpt(el,1)*Hqpt(el,1) + Qxqpt(el,1)*Qxqpt(el,1)*recipHa(el)
              ymom(el,1) = pt5g*Hqpt(el,1)*Hqpt(el,1) + Qyqpt(el,1)*Qyqpt(el,1)*recipHa(el) 
              xymom(el,1) = Qxqpt(el,1)*Qyqpt(el,1)*recipHa(el)

            ENDDO 

!DIR$ VECTOR ALIGNED
            DO el = split(1,sp),split(2,sp) ! Compute source terms
              tau(el) = cf*sqrt((Qxqpt(el,1)*recipHa(el))**2 + (Qyqpt(el,1)*recipHa(el))**2)*recipHa(el)
              src_x(el) = g*Hqpt(el,1)*dhbdx(el) - tau(el)*Qxqpt(el,1) 
              src_y(el) = g*Hqpt(el,1)*dhbdy(el) - tau(el)*Qyqpt(el,1)
            ENDDO

!DIR$ VECTOR ALIGNED
            DO el = split(1,sp),split(2,sp) ! Derivatives are 0 for first dof
              rhsQx(el,1) = rhsQx(el,1) + src_x(el)*phia_int(1,pt)
              rhsQy(el,1) = rhsQy(el,1) + src_y(el)*phia_int(1,pt)
            ENDDO

      test: DO l = 2,ndof 
              ind = (l-1)*nqpta+pt
!DIR$ VECTOR ALIGNED          
              DO el = split(1,sp),split(2,sp)
                rhsH(el,l) = rhsH(el,l) + Qxqpt(el,1)*dpdx(el,ind) + Qyqpt(el,1)*dpdy(el,ind)

                rhsQx(el,l) = rhsQx(el,l) + xmom(el,1)*dpdx(el,ind) + xymom(el,1)*dpdy(el,ind) + src_x(el)*phia_int(l,pt)
                rhsQy(el,l) = rhsQy(el,l) + xymom(el,1)*dpdx(el,ind) + ymom(el,1)*dpdy(el,ind) + src_y(el)*phia_int(l,pt)
              ENDDO

            ENDDO test 

          ENDDO a_points

       ENDDO

      
!     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     c Edge Integrals
!     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
     

!           Hqpt(:,:) = 0d0
!           Qxqpt(:,:) = 0d0
!           Qyqpt(:,:) = 0d0

       DO sp = 1,nsp

ed_points: DO pt = 1,3*nqpte

!     l_edge: DO led = 1,3
!               ind = (led-1)*nqpte + pt
!DIR$ VECTOR ALIGNED               
              DO el = split(1,sp),split(2,sp)
                Hqpt(el,pt) = H(el,1)
                Qxqpt(el,pt) = Qx(el,1)
                Qyqpt(el,pt) = Qy(el,1)
              ENDDO

    ed_basis: DO dof = 2,ndof     
!DIR$ VECTOR ALIGNED
                DO el = split(1,sp),split(2,sp) ! Compute solutions at edge quadrature points
                  Hqpt(el,pt) = Hqpt(el,pt) + H(el,dof)*phie(dof,pt)            
                  Qxqpt(el,pt) = Qxqpt(el,pt) + Qx(el,dof)*phie(dof,pt)            
                  Qyqpt(el,pt) = Qyqpt(el,pt) + Qy(el,dof)*phie(dof,pt)            
                ENDDO

              ENDDO ed_basis

!DIR$ VECTOR ALIGNED
              DO el = split(1,sp),split(2,sp) ! Compute momentum terms
                recipHa(el) = 1d0/Hqpt(el,pt)
                
                xmom(el,pt) = pt5g*Hqpt(el,pt)*Hqpt(el,pt) + Qxqpt(el,pt)*Qxqpt(el,pt)*recipHa(el) 
                ymom(el,pt) = pt5g*Hqpt(el,pt)*Hqpt(el,pt) + Qyqpt(el,pt)*Qyqpt(el,pt)*recipHa(el) 
                xymom(el,pt) = Qxqpt(el,pt)*Qyqpt(el,pt)*recipHa(el)
              ENDDO
              
!               led = (pt-1)/nqpte + 1
! !DIR$ IVDEP
!               DO el = split(1,sp),split(2,sp)
!                 Hn(el,pt) = nxv(el,led)*Qxqpt(el,pt) + nyv(el,led)*Qyqpt(el,pt)
!                 Qxn(el,pt) = nxv(el,led)*xmom(el,pt) + nyv(el,led)*xymom(el,pt)
!                 Qyn(el,pt) = nxv(el,led)*xymom(el,pt) + nyv(el,led)*ymom(el,pt)
!               ENDDO
! 
!               DO el = split(1,sp),split(2,sp)
!                 egnval(el,pt) = ABS(Hn(el,pt))*recipHa(el) + SQRT(g*Hqpt(el,pt))
!               ENDDO

!             ENDDO l_edge

          ENDDO ed_points

       ENDDO
       
!        DO ed = 1,nied
!          DO pt = 1,nqpte
!            const(pt) = max(abs(Qxai(ed)%ptr(pt)*inx(ed) + Qyai(ed)%ptr(pt)*iny(ed))/Hai(ed)%ptr(pt) + sqrt(g*Hai(ed)%ptr(pt)), &
!                            abs(Qxae(ed)%ptr(pt)*inx(ed) + Qyae(ed)%ptr(pt)*iny(ed))/Hae(ed)%ptr(pt) + sqrt(g*Hae(ed)%ptr(pt)))
!          ENDDO
!          
!          DO pt = 1,nqpte
!            Hhatv(pt) = .5d0*(inx(ed)*(Qxai(ed)%ptr(pt) + Qxae(ed)%ptr(pt)) + iny(ed)*(Qyai(ed)%ptr(pt) + Qyae(ed)%ptr(pt)) &
!                                      - const(pt)*(Hae(ed)%ptr(pt) - Hai(ed)%ptr(pt)))         
!          ENDDO
!          
!          DO pt = 1,nqpte
!            Qxhatv(pt) = .5d0*(inx(ed)*(xmai(ed)%ptr(pt) + xmae(ed)%ptr(pt)) + iny(ed)*(xymai(ed)%ptr(pt) + xymae(ed)%ptr(pt))  &
!                                      - const(pt)*(Qxae(ed)%ptr(pt) - Qxai(ed)%ptr(pt)))         
!          ENDDO
!          
!          DO pt = 1,nqpte
!            Qyhatv(pt) = .5d0*(inx(ed)*(xymai(ed)%ptr(pt) + xymae(ed)%ptr(pt)) + iny(ed)*(ymai(ed)%ptr(pt) + ymae(ed)%ptr(pt))  &
!                                      - const(pt)*(Qyae(ed)%ptr(pt) - Qyai(ed)%ptr(pt)))         
!          ENDDO
!          
!          DO pt = 1,nqpte
!            Hfae(ed)%ptr(pt) = -len_area_ex(ed)*Hhatv(pt)
!            Hfai(ed)%ptr(pt) =  len_area_in(ed)*Hhatv(pt)         
!          ENDDO
!          
!          DO pt = 1,nqpte
!            Qxfae(ed)%ptr(pt) = -len_area_ex(ed)*Qxhatv(pt)
!            Qxfai(ed)%ptr(pt) =  len_area_in(ed)*Qxhatv(pt)         
!          ENDDO
!          
!          DO pt = 1,nqpte
!            Qyfae(ed)%ptr(pt) = -len_area_ex(ed)*Qyhatv(pt)
!            Qyfai(ed)%ptr(pt) =  len_area_in(ed)*Qyhatv(pt)         
!          ENDDO         
!        ENDDO
     
     
          DO sp = 1,nsp2
ed_points2: DO pt = 1,nqpte ! Compute numerical fluxes for all edges

!               ! Interior edges
! !DIR$ VECTOR ALIGNED
!               DO ed = split2(1,sp),split2(2,sp)
!                 const(ed) = max(EVi(ed,pt)%ptr, EVe(ed,pt)%ptr)
!               ENDDO
!               
! !DIR$ VECTOR ALIGNED
! !DIR$ IVDEP
!               DO ed = split2(1,sp),split2(2,sp)
!                 Hhatv(ed) = .5d0*(Hni(ed,pt)%ptr + Hne(ed,pt)%ptr &
!                                         - const(ed)*(He(ed,pt)%ptr - Hi(ed,pt)%ptr))
!                 Hfe(ed,pt)%ptr = -len_area_ex(ed)*Hhatv(ed)
!                 Hfi(ed,pt)%ptr =  len_area_in(ed)*Hhatv(ed)
!               ENDDO
! 
! !DIR$ VECTOR ALIGNED
! !DIR$ IVDEP
!               DO ed = split2(1,sp),split2(2,sp)
!                 Qxhatv(ed) = .5d0*(Qxni(ed,pt)%ptr + Qxne(ed,pt)%ptr  &
!                                         - const(ed)*(Qxe(ed,pt)%ptr - Qxi(ed,pt)%ptr))
!                 Qxfe(ed,pt)%ptr = -len_area_ex(ed)*Qxhatv(ed)
!                 Qxfi(ed,pt)%ptr =  len_area_in(ed)*Qxhatv(ed)
!               ENDDO
!               
! !DIR$ VECTOR ALIGNED
! !DIR$ IVDEP
!               DO ed = split2(1,sp),split2(2,sp)
!                 Qyhatv(ed) = .5d0*(Qyni(ed,pt)%ptr + Qyne(ed,pt)%ptr  &
!                                         - const(ed)*(Qye(ed,pt)%ptr - Qyi(ed,pt)%ptr))
!                 Qyfe(ed,pt)%ptr = -len_area_ex(ed)*Qyhatv(ed)
!                 Qyfi(ed,pt)%ptr =  len_area_in(ed)*Qyhatv(ed)
!               ENDDO
!               
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!               ! Interior edges
! !DIR$ VECTOR ALIGNED
!               DO ed = split2(1,sp),split2(2,sp)
!                 const(ed) = max(EVi(ed,pt)%ptr, EVe(ed,pt)%ptr)
! 
!                 Hhatv(ed) = .5d0*(Hni(ed,pt)%ptr + Hne(ed,pt)%ptr &
!                                         - const(ed)*(He(ed,pt)%ptr - Hi(ed,pt)%ptr))
!                 Qxhatv(ed) = .5d0*(Qxni(ed,pt)%ptr + Qxne(ed,pt)%ptr  &
!                                         - const(ed)*(Qxe(ed,pt)%ptr - Qxi(ed,pt)%ptr))
!                 Qyhatv(ed) = .5d0*(Qyni(ed,pt)%ptr + Qyne(ed,pt)%ptr  &
!                                         - const(ed)*(Qye(ed,pt)%ptr - Qyi(ed,pt)%ptr))
!               ENDDO
!               
! !DIR$ VECTOR ALIGNED
! !DIR$ IVDEP
!               DO ed = split2(1,sp),split2(2,sp)
!                 Hfe(ed,pt)%ptr = -len_area_ex(ed)*Hhatv(ed)
!                 Hfi(ed,pt)%ptr =  len_area_in(ed)*Hhatv(ed)
!               ENDDO
! 
! !DIR$ VECTOR ALIGNED
! !DIR$ IVDEP
!               DO ed = split2(1,sp),split2(2,sp)
!                 Qxfe(ed,pt)%ptr = -len_area_ex(ed)*Qxhatv(ed)
!                 Qxfi(ed,pt)%ptr =  len_area_in(ed)*Qxhatv(ed)
!               ENDDO
!               
! !DIR$ VECTOR ALIGNED
! !DIR$ IVDEP
!               DO ed = split2(1,sp),split2(2,sp)
!                 Qyfe(ed,pt)%ptr = -len_area_ex(ed)*Qyhatv(ed)
!                 Qyfi(ed,pt)%ptr =  len_area_in(ed)*Qyhatv(ed)
!               ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

              
!DIR$ VECTOR ALIGNED
              DO ed = split2(1,sp),split2(2,sp)
                const(ed) = max(abs(Qxi(ed,pt)%ptr*inx(ed) + Qyi(ed,pt)%ptr*iny(ed))/Hi(ed,pt)%ptr + sqrt(g*Hi(ed,pt)%ptr), &
                                abs(Qxe(ed,pt)%ptr*inx(ed) + Qye(ed,pt)%ptr*iny(ed))/He(ed,pt)%ptr + sqrt(g*He(ed,pt)%ptr))
              ENDDO
!DIR$ IVDEP
!DIR$ VECTOR ALIGNED
              DO ed = split2(1,sp),split2(2,sp)
                Hhatv(ed) = .5d0*(inx(ed)*(Qxi(ed,pt)%ptr + Qxe(ed,pt)%ptr) + iny(ed)*(Qyi(ed,pt)%ptr + Qye(ed,pt)%ptr) &
                                        - const(ed)*(He(ed,pt)%ptr - Hi(ed,pt)%ptr))
                                        
                Hfe(ed,pt)%ptr = -len_area_ex(ed)*Hhatv(ed)
                Hfi(ed,pt)%ptr =  len_area_in(ed)*Hhatv(ed)
              ENDDO      
!DIR$ IVDEP
!DIR$ VECTOR ALIGNED
              DO ed = split2(1,sp),split2(2,sp)
                Qxhatv(ed) = .5d0*(inx(ed)*(xmi(ed,pt)%ptr + xme(ed,pt)%ptr) + iny(ed)*(xymi(ed,pt)%ptr + xyme(ed,pt)%ptr)  &
                                        - const(ed)*(Qxe(ed,pt)%ptr - Qxi(ed,pt)%ptr))
                                        
                Qxfe(ed,pt)%ptr = -len_area_ex(ed)*Qxhatv(ed)
                Qxfi(ed,pt)%ptr =  len_area_in(ed)*Qxhatv(ed)
              ENDDO   
!DIR$ IVDEP
!DIR$ VECTOR ALIGNED
              DO ed = split2(1,sp),split2(2,sp)
                Qyhatv(ed) = .5d0*(inx(ed)*(xymi(ed,pt)%ptr + xyme(ed,pt)%ptr) + iny(ed)*(ymi(ed,pt)%ptr + yme(ed,pt)%ptr)  &
                                        - const(ed)*(Qye(ed,pt)%ptr - Qyi(ed,pt)%ptr))
                Qyfe(ed,pt)%ptr = -len_area_ex(ed)*Qyhatv(ed)
                Qyfi(ed,pt)%ptr =  len_area_in(ed)*Qyhatv(ed)
              ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              
!              DO ed = split2(1,sp),split2(2,sp)
!               Hin(ed) = Hi(ed,pt)%ptr
!               Qxin(ed) = Qxi(ed,pt)%ptr
!               Qyin(ed) = Qyi(ed,pt)%ptr
!               xmin(ed) = xmi(ed,pt)%ptr
!               ymin(ed) = ymi(ed,pt)%ptr
!               xymin(ed) = xymi(ed,pt)%ptr
! 
!               Hex(ed) = He(ed,pt)%ptr
!               Qxex(ed) = Qxe(ed,pt)%ptr
!               Qyex(ed) = Qye(ed,pt)%ptr
!               xmex(ed) = xme(ed,pt)%ptr
!               ymex(ed) = yme(ed,pt)%ptr
!               xymex(ed) = xyme(ed,pt)%ptr
!             ENDDO
! 
!               ! Interior edges
! !DIR$ VECTOR ALIGNED              
!               DO ed = split2(1,sp),split2(2,sp)
!                 const(ed) = max(abs(Qxin(ed)*inx(ed) + Qyin(ed)*iny(ed))/Hin(ed) + sqrt(g*Hin(ed)), &
!                                 abs(Qxex(ed)*inx(ed) + Qyex(ed)*iny(ed))/Hex(ed) + sqrt(g*Hex(ed)))
!               ENDDO
! !DIR$ VECTOR ALIGNED              
!               DO ed = split2(1,sp),split2(2,sp)     
!                 Hhatv(ed) = .5d0*(inx(ed)*(Qxin(ed) + Qxex(ed)) + iny(ed)*(Qyin(ed) + Qyex(ed)) &
!                                         - const(ed)*(Hex(ed) - Hin(ed)))
!               ENDDO 
! !DIR$ VECTOR ALIGNED              
!               DO ed = split2(1,sp),split2(2,sp)              
!                 Qxhatv(ed) = .5d0*(inx(ed)*(xmin(ed) + xmex(ed)) + iny(ed)*(xymin(ed) + xymex(ed))  &
!                                         - const(ed)*(Qxex(ed) - Qxin(ed)))
!               ENDDO 
! !DIR$ VECTOR ALIGNED                                         
!               DO ed = split2(1,sp),split2(2,sp)              
!                 Qyhatv(ed) = .5d0*(inx(ed)*(xymin(ed) + xymex(ed)) + iny(ed)*(ymin(ed) + ymex(ed))  &
!                                         - const(ed)*(Qyex(ed) - Qyin(ed)))
!               ENDDO
! 
! !DIR$ IVDEP
! !DIR$ VECTOR ALIGNED
!               DO ed = split2(1,sp),split2(2,sp)
!                 Hfe(ed,pt)%ptr = -len_area_ex(ed)*Hhatv(ed)
!                 Hfi(ed,pt)%ptr =  len_area_in(ed)*Hhatv(ed)
!               ENDDO
! !DIR$ IVDEP
! !DIR$ VECTOR ALIGNED
!               DO ed = split2(1,sp),split2(2,sp)
!                 Qxfe(ed,pt)%ptr = -len_area_ex(ed)*Qxhatv(ed)
!                 Qxfi(ed,pt)%ptr =  len_area_in(ed)*Qxhatv(ed)
!               ENDDO
! !DIR$ IVDEP
! !DIR$ VECTOR ALIGNED
!               DO ed = split2(1,sp),split2(2,sp)
!                 Qyfe(ed,pt)%ptr = -len_area_ex(ed)*Qyhatv(ed)
!                 Qyfi(ed,pt)%ptr =  len_area_in(ed)*Qyhatv(ed)
!               ENDDO

     ENDDO ed_points2
   ENDDO
          
          DO pt = 1,nqpte
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

              Hflux(el_in,gp_in) = edlen_area(1,ged)*Hhat

              Qxflux(el_in,gp_in) = edlen_area(1,ged)*Qxhat

              Qyflux(el_in,gp_in) = edlen_area(1,ged)*Qyhat
 
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

              Hflux(el_in,gp_in) = edlen_area(1,ged)*Hhat

              Qxflux(el_in,gp_in) = edlen_area(1,ged)*Qxhat

              Qyflux(el_in,gp_in) = edlen_area(1,ged)*Qyhat

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

              Hflux(el_in,gp_in) = edlen_area(1,ged)*Hhat

              Qxflux(el_in,gp_in) = edlen_area(1,ged)*Qxhat

              Qyflux(el_in,gp_in) = edlen_area(1,ged)*Qyhat                            

             ENDDO

          ENDDO 


          DO sp = 1,nsp             
            DO pt = 1,3*nqpte
               DO l = 1,ndof
!DIR$ VECTOR ALIGNED
                 DO el = split(1,sp),split(2,sp)
                   rhsH(el,l) = rhsH(el,l) - Hflux(el,pt)*phie_int(l,pt)
                   rhsQx(el,l) = rhsQx(el,l) - Qxflux(el,pt)*phie_int(l,pt)
                   rhsQy(el,l) = rhsQy(el,l) - Qyflux(el,pt)*phie_int(l,pt)
                 ENDDO
               ENDDO                                    
           ENDDO
         ENDDO


      RETURN
      END SUBROUTINE rhs2
