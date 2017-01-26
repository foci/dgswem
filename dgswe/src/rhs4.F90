      SUBROUTINE rhs4()

      USE globals,  ONLY: rp,l,ind,el,ne,pt,led,ed,dof,ndof, &
                          blk,nblk,elblk,edblk,nfblk,rnfblk,nrblk,npart, &
                          g,pt5g,tstage,ramp,arg,bfr, & 
                          H,Qx,Qy, &
                          rhsH,rhsQx,rhsQy, &
                          Hqpt,Qxqpt,Qyqpt, &
                          nqpta,wpta,phi,phia,phia_int,dpdx,dpdy,nqpte,wpte,phie,phie_int, &
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
                          ged2led,ged2el,gel2ael,ged, &
                          pressa,recipHa, & 
                          Hi,He,Qxi,Qxe,Qyi,Qye, &
                          xmi,xme,ymi,yme,xymi,xyme, &
                          Hfi,Hfe,Qxfi,Qxfe,Qyfi,Qyfe, &
                          const,inx,iny,len_area_in,len_area_ex, &
                          Hhatv,Qxhatv,Qyhatv, &
                          rHi,rHe,xmomi,xmome,ymomi,ymome,xymomi,xymome, &
                          Hin,Qxin,Qyin, &
                          Hex,Qxex,Qyex, &
                          xmin,ymin,xymin, &
                          xmex,ymex,xymex
!$    USE omp_lib                          
                          
                   
      IMPLICIT NONE
      INTEGER :: i,j

!     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     c Area Integrals
!     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  

!!DIR$ VECTOR ALIGNED
      rhsH(:,:) = 0d0
!!DIR$ VECTOR ALIGNED
      rhsQx(:,:) = 0d0
!!DIR$ VECTOR ALIGNED
      rhsQy(:,:) = 0d0   

      
!$OMP parallel default(none)  &
!$OMP             private(i,ind,blk,pt,dof,el,l,ed,ged,led_in,gp_in,el_in, &
!$OMP                     nx,ny,tx,ty,nx2,ny2,nxny,H_in,H_ex,Qx_in,Qx_ex,Qy_in,Qy_ex, &
!$OMP                     press,recipH,xmom_in,xmom_ex,ymom_in,ymom_ex,xymom_in,xymom_ex, &
!$OMP                     Hhat,Qxhat,Qyhat, &
!$OMP                     Qn,Qt,arg,bfr) &
!$OMP             shared(nqpta,nqpte,ndof,H,Qx,Qy,Hqpt,Qxqpt,Qyqpt, & 
!$OMP                    phi,phia,phie,recipHa,xmom,ymom,xymom,tau,src_x,src_y, &
!$OMP                    dhbdx,dhbdy,cf,phia_int,phie_int,dpdx,dpdy,rhsH,rhsQx,rhsQy, &
!$OMP                    const,Qxi,Qxe,Qyi,Qye,Hi,He,xmi,xme,ymi,yme,xymi,xyme, &
!$OMP                    Hhatv,Qxhatv,Qyhatv,Hfi,Hfe,Qxfi,Qxfe,Qyfi,Qyfe, &
!$OMP                    Hflux,Qxflux,Qyflux,inx,iny,len_area_in,len_area_ex, &
!$OMP                    nnfbed,nfbedn,nfbed,fbedn,nobed,obedn, &
!$OMP                    nfbfr,fbfreq,fbper,fbeq,fbamp_qpt,fbnfact,fbph_qpt, &
!$OMP                    nobfr,obfreq,obper,obeq,obamp_qpt,obnfact,obph_qpt,obdepth_qpt, &
!$OMP                    ramp,tstage, &
!$OMP                    edlen_area,normal, &
!$OMP                    ged2led,gel2ael,ged2el, &
!$OMP                    nblk,npart,nrblk,elblk,edblk,nfblk,rnfblk)

!$OMP do
   

      DO blk = 1,npart
          i = 1
a_points: DO pt = 1,nqpta+3*nqpte

!!DIR$ VECTOR ALIGNED
            DO el = edblk(1,blk),edblk(2,blk) ! First basis function is 1
              Hqpt(el,i) = H(el,1)
              Qxqpt(el,i) = Qx(el,1)
              Qyqpt(el,i) = Qy(el,1)
            ENDDO

   a_basis: DO dof = 2,ndof
!DIR$ IVDEP   
!!DIR$ VECTOR ALIGNED
              DO el = edblk(1,blk),edblk(2,blk) ! Evaluate solution at area quadrature point
                Hqpt(el,i) = Hqpt(el,i) + H(el,dof)*phi(dof,pt)
                Qxqpt(el,i) = Qxqpt(el,i) + Qx(el,dof)*phi(dof,pt) 
                Qyqpt(el,i) = Qyqpt(el,i) + Qy(el,dof)*phi(dof,pt)
              ENDDO

            ENDDO a_basis

!!DIR$ VECTOR ALIGNED
            DO el = edblk(1,blk),edblk(2,blk) ! Compute momentum terms
              recipHa(el) = 1d0/Hqpt(el,i)

              xmom(el,i) = pt5g*Hqpt(el,i)*Hqpt(el,i) + Qxqpt(el,i)*Qxqpt(el,i)*recipHa(el)
              ymom(el,i) = pt5g*Hqpt(el,i)*Hqpt(el,i) + Qyqpt(el,i)*Qyqpt(el,i)*recipHa(el) 
              xymom(el,i) = Qxqpt(el,i)*Qyqpt(el,i)*recipHa(el)
            ENDDO 
            
          IF (pt <= nqpta) THEN

!!DIR$ VECTOR ALIGNED
            DO el = edblk(1,blk),edblk(2,blk) ! Compute source terms
              tau(el) = cf*sqrt((Qxqpt(el,1)*recipHa(el))**2 + (Qyqpt(el,1)*recipHa(el))**2)*recipHa(el)
              src_x(el) = g*Hqpt(el,1)*dhbdx(el) - tau(el)*Qxqpt(el,1) 
              src_y(el) = g*Hqpt(el,1)*dhbdy(el) - tau(el)*Qyqpt(el,1)
            ENDDO

!!DIR$ VECTOR ALIGNED
            DO el = edblk(1,blk),edblk(2,blk) ! Derivatives are 0 for first dof
              rhsQx(el,1) = rhsQx(el,1) + src_x(el)*phia_int(1,pt)
              rhsQy(el,1) = rhsQy(el,1) + src_y(el)*phia_int(1,pt)
            ENDDO

      test: DO l = 2,ndof 
              ind = (l-1)*nqpta+pt
!!DIR$ VECTOR ALIGNED          
              DO el = edblk(1,blk),edblk(2,blk)
                rhsH(el,l) = rhsH(el,l) + Qxqpt(el,1)*dpdx(el,ind) + Qyqpt(el,1)*dpdy(el,ind)

                rhsQx(el,l) = rhsQx(el,l) + xmom(el,1)*dpdx(el,ind) + xymom(el,1)*dpdy(el,ind) + src_x(el)*phia_int(l,pt)           
                rhsQy(el,l) = rhsQy(el,l) + xymom(el,1)*dpdx(el,ind) + ymom(el,1)*dpdy(el,ind) + src_y(el)*phia_int(l,pt)
              ENDDO

            ENDDO test 
            
          ELSEIF (pt > nqpta .and. pt < nqpta+3*nqpte) THEN
          
            i = i+1         
          
          ELSE 

ed_points2: DO j = 1,nqpte ! Compute numerical fluxes for all edges
              
!!DIR$ VECTOR ALIGNED
              DO ed = nfblk(1,blk),nfblk(2,blk) 
                
                const(ed) = max(abs(Qxi(ed,j)%ptr*inx(ed) + Qyi(ed,j)%ptr*iny(ed))/Hi(ed,j)%ptr + sqrt(g*Hi(ed,j)%ptr), &
                                abs(Qxe(ed,j)%ptr*inx(ed) + Qye(ed,j)%ptr*iny(ed))/He(ed,j)%ptr + sqrt(g*He(ed,j)%ptr))
              ENDDO
                                         
!DIR$ IVDEP
!!DIR$ VECTOR ALIGNED
              DO ed = nfblk(1,blk),nfblk(2,blk)             
                Hhatv(ed) = .5d0*(inx(ed)*(Qxi(ed,j)%ptr + Qxe(ed,j)%ptr) + iny(ed)*(Qyi(ed,j)%ptr + Qye(ed,j)%ptr) &
                                        - const(ed)*(He(ed,j)%ptr - Hi(ed,j)%ptr))

                Hfe(ed,j)%ptr = -len_area_ex(ed)*Hhatv(ed)
                Hfi(ed,j)%ptr =  len_area_in(ed)*Hhatv(ed)
              ENDDO    
                        
!DIR$ IVDEP
!!DIR$ VECTOR ALIGNED
              DO ed = nfblk(1,blk),nfblk(2,blk)          
                Qxhatv(ed) = .5d0*(inx(ed)*(xmi(ed,j)%ptr + xme(ed,j)%ptr) + iny(ed)*(xymi(ed,j)%ptr + xyme(ed,j)%ptr)  &
                                        - const(ed)*(Qxe(ed,j)%ptr - Qxi(ed,j)%ptr))
                                 
                Qxfe(ed,j)%ptr = -len_area_ex(ed)*Qxhatv(ed)
                Qxfi(ed,j)%ptr =  len_area_in(ed)*Qxhatv(ed)
              ENDDO   
                           
!DIR$ IVDEP
!!DIR$ VECTOR ALIGNED
              DO ed = nfblk(1,blk),nfblk(2,blk)
                Qyhatv(ed) = .5d0*(inx(ed)*(xymi(ed,j)%ptr + xyme(ed,j)%ptr) + iny(ed)*(ymi(ed,j)%ptr + yme(ed,j)%ptr)  &
                                        - const(ed)*(Qye(ed,j)%ptr - Qyi(ed,j)%ptr))
                                     
                Qyfe(ed,j)%ptr = -len_area_ex(ed)*Qyhatv(ed)
                Qyfi(ed,j)%ptr =  len_area_in(ed)*Qyhatv(ed)
              ENDDO                           

            ENDDO ed_points2            
          
          ENDIF

        ENDDO a_points

       ENDDO
      
            
!$OMP do         
      DO blk = 1,nrblk
       DO pt = 1,nqpte
!!DIR$ VECTOR ALIGNED
              DO ed = rnfblk(1,blk),rnfblk(2,blk)
!               DO ed = nfblk(1,npart+1),nfblk(2,npart+1)
                const(ed) = max(abs(Qxi(ed,pt)%ptr*inx(ed) + Qyi(ed,pt)%ptr*iny(ed))/Hi(ed,pt)%ptr + sqrt(g*Hi(ed,pt)%ptr), &
                                abs(Qxe(ed,pt)%ptr*inx(ed) + Qye(ed,pt)%ptr*iny(ed))/He(ed,pt)%ptr + sqrt(g*He(ed,pt)%ptr))
              ENDDO
              
                            
!DIR$ IVDEP
!!DIR$ VECTOR ALIGNED
              DO ed = rnfblk(1,blk),rnfblk(2,blk)
!               DO ed = nfblk(1,npart+1),nfblk(2,npart+1)
                Hhatv(ed) = .5d0*(inx(ed)*(Qxi(ed,pt)%ptr + Qxe(ed,pt)%ptr) + iny(ed)*(Qyi(ed,pt)%ptr + Qye(ed,pt)%ptr) &
                                        - const(ed)*(He(ed,pt)%ptr - Hi(ed,pt)%ptr))
                                     
                Hfe(ed,pt)%ptr = -len_area_ex(ed)*Hhatv(ed)
                Hfi(ed,pt)%ptr =  len_area_in(ed)*Hhatv(ed)               
              ENDDO     
              
              
              
!DIR$ IVDEP
!!DIR$ VECTOR ALIGNED
              DO ed = rnfblk(1,blk),rnfblk(2,blk)
!               DO ed = nfblk(1,npart+1),nfblk(2,npart+1)
                Qxhatv(ed) = .5d0*(inx(ed)*(xmi(ed,pt)%ptr + xme(ed,pt)%ptr) + iny(ed)*(xymi(ed,pt)%ptr + xyme(ed,pt)%ptr)  &
                                        - const(ed)*(Qxe(ed,pt)%ptr - Qxi(ed,pt)%ptr))
                                   
                Qxfe(ed,pt)%ptr = -len_area_ex(ed)*Qxhatv(ed)
                Qxfi(ed,pt)%ptr =  len_area_in(ed)*Qxhatv(ed)
              ENDDO   
              
              
              
!DIR$ IVDEP
!!DIR$ VECTOR ALIGNED
              DO ed = rnfblk(1,blk),rnfblk(2,blk)
!               DO ed = nfblk(1,npart+1),nfblk(2,npart+1)
                Qyhatv(ed) = .5d0*(inx(ed)*(xymi(ed,pt)%ptr + xyme(ed,pt)%ptr) + iny(ed)*(ymi(ed,pt)%ptr + yme(ed,pt)%ptr)  &
                                        - const(ed)*(Qye(ed,pt)%ptr - Qyi(ed,pt)%ptr))
                                      
                Qyfe(ed,pt)%ptr = -len_area_ex(ed)*Qyhatv(ed)
                Qyfi(ed,pt)%ptr =  len_area_in(ed)*Qyhatv(ed)
              ENDDO

        ENDDO       
        
      ENDDO
!$OMP end do


!$OMP do     
!           DO pt = 1,nqpte
            ! No normal flow boundary condition 
            DO ed = 1,nnfbed
              DO pt = 1,nqpte

              ged = nfbedn(ed)

              led_in = ged2led(1,ged)

              gp_in = (led_in-1)*nqpte + pt

              el_in = gel2ael(ged2el(1,ged))

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

              CALL numerical_flux(Qx_in,Qy_in,H_in,Qx_ex,Qy_ex,H_ex,nx,ny,Qxhat,Qyhat,Hhat)

              Hflux(el_in,gp_in) = edlen_area(1,ged)*Hhat

              Qxflux(el_in,gp_in) = edlen_area(1,ged)*Qxhat

              Qyflux(el_in,gp_in) = edlen_area(1,ged)*Qyhat
              
              ENDDO
            ENDDO
!$OMP end do

!$OMP do

            ! Flow specified boundary edges
            DO ed = 1,nfbed
              DO pt = 1,nqpte

              ged = fbedn(ed)

              led_in = ged2led(1,ged)

              gp_in = (led_in-1)*nqpte + pt

              el_in = gel2ael(ged2el(1,ged))
              
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
 
              CALL numerical_flux(Qx_in,Qy_in,H_in,Qx_ex,Qy_ex,H_ex,nx,ny,Qxhat,Qyhat,Hhat)

              Hflux(el_in,gp_in) = edlen_area(1,ged)*Hhat
              
              Qxflux(el_in,gp_in) = edlen_area(1,ged)*Qxhat

              Qyflux(el_in,gp_in) = edlen_area(1,ged)*Qyhat
              
              ENDDO
            ENDDO
            
!$OMP end do

!$OMP do

             ! Open boundary edges (elevation specified)
            DO ed = 1,nobed
              DO pt = 1,nqpte

              ged = obedn(ed)

              led_in = ged2led(1,ged)

              gp_in = (led_in-1)*nqpte + pt

              el_in = gel2ael(ged2el(1,ged))

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
 
              CALL numerical_flux(Qx_in,Qy_in,H_in,Qx_ex,Qy_ex,H_ex,nx,ny,Qxhat,Qyhat,Hhat)

              Hflux(el_in,gp_in) = edlen_area(1,ged)*Hhat

              Qxflux(el_in,gp_in) = edlen_area(1,ged)*Qxhat

              Qyflux(el_in,gp_in) = edlen_area(1,ged)*Qyhat                            
              ENDDO
             ENDDO

!           ENDDO 
!$OMP end do        


!$OMP do          
          DO blk = 1,nblk             
            DO pt = 1,3*nqpte
               DO l = 1,ndof
!!DIR$ VECTOR ALIGNED
                 DO el = elblk(1,blk),elblk(2,blk)
                   rhsH(el,l) = rhsH(el,l) - Hflux(el,pt)*phie_int(l,pt)
                   rhsQx(el,l) = rhsQx(el,l) - Qxflux(el,pt)*phie_int(l,pt)
                   rhsQy(el,l) = rhsQy(el,l) - Qyflux(el,pt)*phie_int(l,pt)                   
                 ENDDO
               ENDDO                                    
           ENDDO
         ENDDO
!$OMP end do  



!$OMP end parallel         

      RETURN
      END SUBROUTINE rhs4
