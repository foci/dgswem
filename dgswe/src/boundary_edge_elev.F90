      SUBROUTINE boundary_edge_elev()

      USE globals, ONLY: rp,el_type,ndof,nqpte, &
                         nobed,obedn,nobfr, &
                         tstage,ramp, &
                         obfreq,obper,obeq, &
                         obnfact,obamp_qpt,obph_qpt, &
                         ged2led,ged2el,gel2ael, &
                         nx_pt,ny_pt,Spe,detJe,hbqpted, &
                         Zqpt,Qxqpt,Qyqpt, &
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

      ! Open boundary edges (elevation specified)
      DO ed = 1,nobed
            
        ged = obedn(ed)
        led_in = ged2led(1,ged)
        el_in = gel2ael(ged2el(1,ged))
              
        DO pt = 1,nqpte(1)

          gp_in = (led_in-1)*nqpte(1) + pt

          nx = nx_pt(ged,pt)
          ny = ny_pt(ged,pt)
          sp = Spe(ged,pt)
              
          hb = hbqpted(ged,pt)                 

          Z_in = Zqpt(el_in,gp_in)

          Z_ex = 0d0
          DO bfr = 1,nobfr
            arg = obfreq(bfr)*(tstage-real(INT(tstage/obper(bfr)),rp)*obper(bfr)) + obeq(bfr)
            Z_ex = Z_ex + obamp_qpt(bfr,pt,ed)*obnfact(bfr)*ramp*COS(arg-obph_qpt(bfr,pt,ed))
          ENDDO

          Qx_in = Qxqpt(el_in,gp_in)
          Qy_in = Qyqpt(el_in,gp_in)

          Qx_ex = Qx_in
          Qy_ex = Qy_in
 
          CALL numerical_flux(Qx_in,Qy_in,Z_in,Qx_ex,Qy_ex,Z_ex,hb,nx,ny,sp,Qxhat,Qyhat,Zhat)
          
          Exx_in = Exxqpt(el_in,gp_in)
          Eyy_in = Eyyqpt(el_in,gp_in)
          Exy_in = Exyqpt(el_in,gp_in)

          Zqpt(el_in,gp_in) = detJe(ged,pt)*Zhat
          Qxqpt(el_in,gp_in) = detJe(ged,pt)*(Qxhat - nx*Exx_in - ny*Exy_in)
          Qyqpt(el_in,gp_in) = detJe(ged,pt)*(Qyhat - nx*Exy_in - ny*Eyy_in)

!           Qxqpt(el_in,gp_in) = detJe(ged,pt)*(Qxhat - nx*Exx_in - ny*Exy_in)
!           Qyqpt(el_in,gp_in) = detJe(ged,pt)*(Qyhat - nx*Eyx_in - ny*Eyy_in)
                          
          ENDDO
        ENDDO
!$OMP end do  

      END SUBROUTINE boundary_edge_elev
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      SUBROUTINE boundary_edge_elev_ldg()
      
      USE globals, ONLY: rp,nobed,obedn, &
                         ged2led,ged2el,gel2ael, &
                         nqpte, &
                         nx_pt,ny_pt,Spe, &
                         Qxqpt,Qyqpt, &
                         Exxqpt,Eyyqpt,Exyqpt,Eyxqpt, &
                         detJe
      
      IMPLICIT NONE
      
      INTEGER :: ed,pt
      INTEGER :: ged,led_in,gp_in,el_in
      REAL(rp) :: nx,ny,sp
      REAL(rp) :: Qx_ex,Qy_ex
      
      DO ed = 1,nobed
      
         ged = obedn(ed)    
         led_in = ged2led(1,ged)
         el_in =  gel2ael(ged2el(1,ged))
      
         DO pt = 1,nqpte(1)
           
           gp_in = (led_in-1)*nqpte(1) + pt           
           
           nx = nx_pt(ged,pt)
           ny = ny_pt(ged,pt)
           sp = Spe(ged,pt)
           
           Qx_ex = Qxqpt(el_in,gp_in)
           Qy_ex = Qyqpt(el_in,gp_in)  
           
           Exxqpt(el_in,gp_in) = detJe(ged,pt)*nx*Qx_ex
           Eyyqpt(el_in,gp_in) = detJe(ged,pt)*ny*Qy_ex
           Exyqpt(el_in,gp_in) = detJe(ged,pt)*(ny*Qx_ex+nx*Qy_ex)          
!            Exyqpt(el_in,gp_in) = detJe(ged,pt)*ny*Qx_ex
!            Eyxqpt(el_in,gp_in) = detJe(ged,pt)*nx*Qy_ex
           
         ENDDO
      ENDDO
      
      RETURN
      END SUBROUTINE boundary_edge_elev_ldg
      