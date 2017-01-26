      SUBROUTINE boundary_edge_flow()

      USE globals, ONLY: rp,el_type,ndof,nqpte,pi, &
                         nfbed,fbedn,nfbfr, &
                         tstage,ramp, &
                         fbfreq,fbper,fbeq, &
                         fbnfact,fbamp_qpt,fbph_qpt, &
                         ged2led,ged2el,gel2ael, &
                         nx_pt,ny_pt,Spe,detJe,hbqpted, &
                         Zqpt,Qxqpt,Qyqpt, &
                         nfbsfr,fbsamp_qpt,fbsbgn,fbsend,fbssig, &
                         Exxqpt,Eyyqpt,Exyqpt,Eyxqpt

      IMPLICIT NONE

      INTEGER :: ed,pt,bfr
      INTEGER :: ged
      INTEGER :: led_in,el_in,gp_in
      REAL(rp) :: Z_in,Z_ex,Qx_in,Qx_ex,Qy_in,Qy_ex
      REAL(rp) :: Zhat,Qxhat,Qyhat
      REAL(rp) :: tx,ty,nx,ny,nx2,ny2,nxny,sp,hb
      REAL(rp) :: arg,Qn,Qt
      REAL(rp) :: Exx_in,Eyy_in,Exy_in,Eyx_in
      

!$OMP do

      ! Flow specified boundary edges
      DO ed = 1,nfbed
            
        ged = fbedn(ed)
        led_in = ged2led(1,ged)  
        el_in = gel2ael(ged2el(1,ged))      
              
        DO pt = 1,nqpte(1)

          gp_in = (led_in-1)*nqpte(1) + pt
       
          nx = nx_pt(ged,pt)
          ny = ny_pt(ged,pt)
          sp = Spe(ged,pt)      
              
          hb = hbqpted(ged,pt)              

          tx = -ny
          ty = nx

          Z_in = Zqpt(el_in,gp_in)
          Z_ex = Z_in

          Qx_in = Qxqpt(el_in,gp_in)
          Qy_in = Qyqpt(el_in,gp_in)

          Qn = 0d0
          Qt = 0d0
          DO bfr = 1,nfbfr
            arg = fbfreq(bfr)*(tstage - real(INT(tstage/fbper(bfr)),rp)*fbper(bfr)) + fbeq(bfr)
            Qn = Qn + fbamp_qpt(bfr,pt,ed)*fbnfact(bfr)*ramp*COS(arg-fbph_qpt(bfr,pt,ed))
          ENDDO
          
          DO bfr = 1,nfbsfr
              Qn = Qn + .5d0*(fbsamp_qpt(bfr,pt,ed)*exp(-((tstage-fbsbgn(bfr))/fbssig(bfr))**2) &
                             -fbsamp_qpt(bfr,pt,ed)*exp(-((tstage-fbsend(bfr))/fbssig(bfr))**2))
          ENDDO          

          Qn = -Qn
          
          Qx_ex = ( ty*Qn - ny*Qt)/(nx*ty-ny*tx)
          Qy_ex = (-tx*Qn + nx*Qt)/(nx*ty-ny*tx)
 
          CALL numerical_flux(Qx_in,Qy_in,Z_in,Qx_ex,Qy_ex,Z_ex,hb,nx,ny,sp,Qxhat,Qyhat,Zhat)
          
          Exx_in = Exxqpt(el_in,gp_in)
          Eyy_in = Eyyqpt(el_in,gp_in)
          Exy_in = Exyqpt(el_in,gp_in)
!           Eyx_in = Eyxqpt(el_in,gp_in)          

          Zqpt(el_in,gp_in) = detJe(ged,pt)*Zhat              
          Qxqpt(el_in,gp_in) = detJe(ged,pt)*(Qxhat - nx*Exx_in - ny*Exy_in)
          Qyqpt(el_in,gp_in) = detJe(ged,pt)*(Qyhat - nx*Exy_in - ny*Eyy_in)
          
!           Qxqpt(el_in,gp_in) = detJe(ged,pt)*(Qxhat - nx*Exx_in - ny*Exy_in)
!           Qyqpt(el_in,gp_in) = detJe(ged,pt)*(Qyhat - nx*Eyx_in - ny*Eyy_in)

!           Qxqpt(el_in,gp_in) = detJe(ged,pt)*Qxhat
!           Qyqpt(el_in,gp_in) = detJe(ged,pt)*Qyhat
              
        ENDDO
      ENDDO
            
!$OMP end do

      END SUBROUTINE boundary_edge_flow

      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      SUBROUTINE boundary_edge_flow_ldg()
      
      USE globals, ONLY: rp,nfbed,fbedn, &
                         nqpte, &
                         ged2led,gel2ael,ged2el, &
                         detJe,nx_pt,ny_pt,Spe, &
                         tstage,ramp, &
                         nfbfr,fbfreq,fbper,fbeq,fbnfact, &
                         fbamp_qpt,fbph_qpt, &
                         nfbsfr,fbsamp_qpt,fbsbgn,fbsend,fbssig, &
                         Qxqpt,Qyqpt, &
                         Exxqpt,Eyyqpt,Exyqpt,Eyxqpt

      IMPLICIT NONE
      
      INTEGER :: ed,pt,bfr
      INTEGER :: ged
      INTEGER :: led_in,el_in,gp_in
      REAL(rp) :: nx,ny,sp
      REAL(rp) :: tx,ty
      REAL(rp) :: Qn,Qt
      REAL(rp) :: arg
      REAL(rp) :: Qx_in,Qy_in,Qx_ex,Qy_ex
      
      DO ed = 1,nfbed
        
        ged = fbedn(ed)
        led_in = ged2led(1,ged)
        el_in = gel2ael(ged2el(1,ged))
      
        DO pt = 1,nqpte(1)
        
          gp_in = (led_in-1)*nqpte(1) + pt
          
          nx = nx_pt(ged,pt)
          ny = ny_pt(ged,pt)
          sp = Spe(ged,pt)
          
          tx = -ny
          ty = nx          
          
          Qx_in = Qxqpt(el_in,gp_in)
          Qy_in = Qyqpt(el_in,gp_in)          
          
          Qn = 0d0
          Qt = 0d0
          DO bfr = 1,nfbfr
            arg = fbfreq(bfr)*(tstage - real(INT(tstage/fbper(bfr)),rp)*fbper(bfr)) + fbeq(bfr)
            Qn = Qn + fbamp_qpt(bfr,pt,ed)*fbnfact(bfr)*ramp*COS(arg-fbph_qpt(bfr,pt,ed))
          ENDDO
          
          DO bfr = 1,nfbsfr
              Qn = Qn + .5d0*(fbsamp_qpt(bfr,pt,ed)*exp(-((tstage-fbsbgn(bfr))/fbssig(bfr))**2) &
                             -fbsamp_qpt(bfr,pt,ed)*exp(-((tstage-fbsend(bfr))/fbssig(bfr))**2))
          ENDDO          

          Qn = -Qn
                    
          Qx_ex = ( ty*Qn - ny*Qt)/(nx*ty-ny*tx)
          Qy_ex = (-tx*Qn + nx*Qt)/(nx*ty-ny*tx)
          
!           Exxqpt(el_in,gp_in) = detJe(ged,pt)*nx*Qx_ex
!           Eyyqpt(el_in,gp_in) = detJe(ged,pt)*ny*Qy_ex
! !           Exyqpt(el_in,gp_in) = detJe(ged,pt)*(ny*Qx_ex + nx*Qy_ex)
  

          Exxqpt(el_in,gp_in) = detJe(ged,pt)*nx*.5d0*(Qx_ex+Qx_in)
          Eyyqpt(el_in,gp_in) = detJe(ged,pt)*ny*.5d0*(Qy_ex+Qy_in)
          Exyqpt(el_in,gp_in) = detJe(ged,pt)*.5d0*(ny*(Qx_ex+Qx_in)+nx*(Qy_ex+Qy_in))          
!           Exyqpt(el_in,gp_in) = detJe(ged,pt)*ny*.5d0*(Qx_ex+Qx_in)
!           Eyxqpt(el_in,gp_in) = detJe(ged,pt)*nx*.5d0*(Qy_ex+Qy_in)
        
        ENDDO
      
      ENDDO
      
      
      RETURN
      END SUBROUTINE boundary_edge_flow_ldg