      SUBROUTINE boundary_edge_land()

      USE globals, ONLY: rp,el_type,ndof,nqpte,nnfbed,nfbedn, &
                         ged2led,ged2el,gel2ael, &
                         nx_pt,ny_pt,Spe,detJe,hbqpted, &
                         Z,Qx,Qy, &
                         Zqpt,Qxqpt,Qyqpt, &
                         rhsZ,rhsQx,rhsQy, &
                         phie,phie_int, &
                         Exxqpt,Eyyqpt,Exyqpt,Eyxqpt
                         
      USE read_dginp, ONLY: ctp                         
                         

      IMPLICIT NONE

      INTEGER :: ed,pt,l,dof
      INTEGER :: ged,et
      INTEGER :: led_in,el_in,gp_in
      REAL(rp) :: Z_in,Z_ex,Qx_in,Qx_ex,Qy_in,Qy_ex
      REAL(rp) :: Zhat,Qxhat,Qyhat
      REAL(rp) :: nx,ny,nx2,ny2,nxny,sp,hb
      REAL(rp) :: Exx_in,Eyy_in,Exy_in,Eyx_in

      IF (ctp == 1) THEN
!$OMP do     
        ! No normal flow boundary condition 
        DO ed = 1,nnfbed
            
          ged = nfbedn(ed)
          led_in = ged2led(1,ged)        
          el_in = gel2ael(ged2el(1,ged))  
              
          DO pt = 1,nqpte(1)

            gp_in = (led_in-1)*nqpte(1) + pt

            nx = nx_pt(ged,pt)
            ny = ny_pt(ged,pt)
            sp = Spe(ged,pt)
              
            hb = hbqpted(ged,pt)         

            nx2 = nx*nx
            ny2 = ny*ny
            nxny = nx*ny

            Z_in = Zqpt(el_in,gp_in)
            Z_ex = Z_in              

            Qx_in = Qxqpt(el_in,gp_in)
            Qy_in = Qyqpt(el_in,gp_in)

            Qx_ex = Qx_in*(ny2-nx2) - 2d0*nxny*Qy_in
            Qy_ex = Qy_in*(nx2-ny2) - 2d0*nxny*Qx_in

            CALL numerical_flux(Qx_in,Qy_in,Z_in,Qx_ex,Qy_ex,Z_ex,hb,nx,ny,sp,Qxhat,Qyhat,Zhat)
            
            Exx_in = Exxqpt(el_in,gp_in)
            Eyy_in = Eyyqpt(el_in,gp_in)
            Exy_in = Exyqpt(el_in,gp_in)
!             Eyx_in = Eyxqpt(el_in,gp_in)            

            Zqpt(el_in,gp_in) = detJe(ged,pt)*Zhat
            Qxqpt(el_in,gp_in) = detJe(ged,pt)*(Qxhat - nx*Exx_in - ny*Exy_in)
            Qyqpt(el_in,gp_in) = detJe(ged,pt)*(Qyhat - nx*Exy_in - ny*Eyy_in)
!             Qxqpt(el_in,gp_in) = detJe(ged,pt)*(Qxhat - nx*Exx_in - ny*Exy_in)
!             Qyqpt(el_in,gp_in) = detJe(ged,pt)*(Qyhat - nx*Eyx_in - ny*Eyy_in)

!             Qxqpt(el_in,gp_in) = detJe(ged,pt)*Qxhat
!             Qyqpt(el_in,gp_in) = detJe(ged,pt)*Qyhat
              
          ENDDO
        ENDDO
!$OMP end do
      ELSE
!$OMP do     
        ! No normal flow boundary condition 
        DO ed = 1,nnfbed
            
          ged = nfbedn(ed)
          led_in = ged2led(1,ged)            
          el_in = gel2ael(ged2el(1,ged))
          et = el_type(ged2el(1,ged)) 
              
          DO pt = 1,nqpte(et)

            gp_in = (led_in-1)*nqpte(et) + pt

            nx = nx_pt(ged,pt)
            ny = ny_pt(ged,pt)
            sp = Spe(ged,pt) 
                
            hb = hbqpted(ged,pt)                

            nx2 = nx*nx
            ny2 = ny*ny
            nxny = nx*ny

            Z_in = Z(el_in,1)
            Qx_in = Qx(el_in,1)
            Qy_in = Qy(el_in,1)
            DO dof = 2,ndof(et)
              Z_in  = Z_in  + Z(el_in,dof)*phie(dof,gp_in,et)
              Qx_in = Qx_in + Qx(el_in,dof)*phie(dof,gp_in,et)
              Qy_in = Qy_in + Qy(el_in,dof)*phie(dof,gp_in,et)
            ENDDO            
              
            Z_ex = Z_in
            Qx_ex = Qx_in*(ny2-nx2) - 2d0*nxny*Qy_in
            Qy_ex = Qy_in*(nx2-ny2) - 2d0*nxny*Qx_in

            CALL numerical_flux(Qx_in,Qy_in,Z_in,Qx_ex,Qy_ex,Z_ex,hb,nx,ny,sp,Qxhat,Qyhat,Zhat)
            
            Exx_in = Exxqpt(el_in,gp_in)
            Eyy_in = Eyyqpt(el_in,gp_in)
            Exy_in = Exyqpt(el_in,gp_in)
!             Eyx_in = Eyxqpt(el_in,gp_in)            
           
            Zhat  = detJe(ged,pt)*Zhat              
            Qxhat = detJe(ged,pt)*(Qxhat - nx*Exx_in - ny*Exy_in)
            Qyhat = detJe(ged,pt)*(Qyhat - nx*Exy_in - ny*Eyy_in)

!             Qxhat = detJe(ged,pt)*(Qxhat - nx*Exx_in - ny*Exy_in)
!             Qyhat = detJe(ged,pt)*(Qyhat - nx*Eyx_in - ny*Eyy_in)
              
            DO l = 1,ndof(et)
              rhsZ(el_in,l)  = rhsZ(el_in,l)  - Zhat*phie_int(l,gp_in,et)
              rhsQx(el_in,l) = rhsQx(el_in,l) - Qxhat*phie_int(l,gp_in,et)
              rhsQy(el_in,l) = rhsQy(el_in,l) - Qyhat*phie_int(l,gp_in,et)                   
            ENDDO  
            
            Zqpt(el_in,gp_in) = 0d0
            Qxqpt(el_in,gp_in) = 0d0
            Qyqpt(el_in,gp_in) = 0d0          
              
          ENDDO
        ENDDO
!$OMP end do
      ENDIF     


      END SUBROUTINE boundary_edge_land
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      SUBROUTINE boundary_edge_land_ldg()
      
      USE globals, ONLY: rp,nnfbed,nfbedn, &
                         ged2led,gel2ael,ged2el, &
                         el_type,nqpte,ndof, &
                         nx_pt,ny_pt,Spe,detJe, &
                         Qx,Qy,phie,phie_int, &                         
                         Qxqpt,Qyqpt, &
                         Exxqpt,Eyyqpt,Exyqpt,Eyxqpt, &
                         rhsExx,rhsEyy,rhsExy,rhsEyx
                         
      USE read_dginp, ONLY: ctp                          
      
      IMPLICIT NONE
      
      INTEGER :: ed,pt,l,dof
      INTEGER :: ged,et
      INTEGER :: led_in,el_in,gp_in
      REAL(rp) :: nx,ny,sp
      REAL(rp) :: nx2,ny2,nxny
      REAL(rp) :: Qx_in,Qy_in
      REAL(rp) :: Qx_ex,Qy_ex
      REAL(rp) :: Exx_hat,Eyy_hat,Exy_hat,Eyx_hat
      
      
      IF (ctp == 1) THEN
      
        DO ed = 1,nnfbed
        
          ged = nfbedn(ed) 
          led_in = ged2led(1,ged)
          el_in = gel2ael(ged2el(1,ged))
          
          DO pt = 1,nqpte(1)
          
            gp_in = (led_in-1)*nqpte(1) + pt          
            
            nx = nx_pt(ged,pt)
            ny = ny_pt(ged,pt)
            sp = Spe(ged,pt)          
            
            nx2 = nx*nx
            ny2 = ny*ny
            nxny = nx*ny     
            
            Qx_in = Qxqpt(el_in,gp_in)
            Qy_in = Qyqpt(el_in,gp_in)

            Qx_ex = Qx_in*(ny2-nx2) - 2d0*nxny*Qy_in
            Qy_ex = Qy_in*(nx2-ny2) - 2d0*nxny*Qx_in
            
!             Exxqpt(el_in,gp_in) = detJe(ged,pt)*nx*Qx_ex
!             Eyyqpt(el_in,gp_in) = detJe(ged,pt)*ny*Qy_ex
! !             Exyqpt(el_in,gp_in) = detJe(ged,pt)*(ny*Qx_ex + nx*Qy_ex)  
            
            Exxqpt(el_in,gp_in) = detJe(ged,pt)*nx*.5d0*(Qx_ex+Qx_in)
            Eyyqpt(el_in,gp_in) = detJe(ged,pt)*ny*.5d0*(Qy_ex+Qy_in)
            Exyqpt(el_in,gp_in) = detJe(ged,pt)*.5d0*(ny*(Qx_ex+Qx_in)+nx*(Qy_ex+Qy_in))            
!             Exyqpt(el_in,gp_in) = detJe(ged,pt)*ny*.5d0*(Qx_ex+Qx_in)
!             Eyxqpt(el_in,gp_in) = detJe(ged,pt)*nx*.5d0*(Qy_ex+Qy_in)             
          
          ENDDO
          
        ENDDO
      
      ELSE
       
        DO ed = 1,nnfbed
        
          ged = nfbedn(ed) 
          led_in = ged2led(1,ged)
          el_in = gel2ael(ged2el(1,ged))  
          et = el_type(ged2el(1,ged))
        
          DO pt = 1,nqpte(et)
          
            gp_in = (led_in-1)*nqpte(et) + pt

            nx = nx_pt(ged,pt)
            ny = ny_pt(ged,pt)
            sp = Spe(ged,pt)           
            
            nx2 = nx*nx
            ny2 = ny*ny
            nxny = nx*ny
            
            Qx_in = Qx(el_in,1)
            Qy_in = Qy(el_in,1)
            DO dof = 2,ndof(et)
              Qx_in = Qx_in + Qx(el_in,dof)*phie(dof,gp_in,et)
              Qy_in = Qy_in + Qy(el_in,dof)*phie(dof,gp_in,et)
            ENDDO               
            
            Qx_ex = Qx_in*(ny2-nx2) - 2d0*nxny*Qy_in
            Qy_ex = Qy_in*(nx2-ny2) - 2d0*nxny*Qx_in    
            
            Exx_hat = detJe(ged,pt)*nx*.5d0*(Qx_ex+Qx_in)
            Eyy_hat = detJe(ged,pt)*ny*.5d0*(Qy_ex+Qy_in)
            Exy_hat = detJe(ged,pt)*.5d0*(ny*(Qx_ex+Qx_in)+nx*(Qy_ex+Qy_in))            
!             Exy_hat = detJe(ged,pt)*ny*.5d0*(Qx_ex+Qx_in)
!             Eyx_hat = detJe(ged,pt)*nx*.5d0*(Qy_ex+Qy_in)
            
            DO l = 1,ndof(et)
              rhsExx(el_in,l) = rhsExx(el_in,l) + Exx_hat*phie_int(l,gp_in,et)
              rhsEyy(el_in,l) = rhsEyy(el_in,l) + Eyy_hat*phie_int(l,gp_in,et)
              rhsExy(el_in,l) = rhsExy(el_in,l) + Exy_hat*phie_int(l,gp_in,et)   
              rhsEyx(el_in,l) = rhsEyx(el_in,l) + Eyx_hat*phie_int(l,gp_in,et)                
            ENDDO  
            
            Exxqpt(el_in,gp_in) = 0d0
            Eyyqpt(el_in,gp_in) = 0d0         
            Exyqpt(el_in,gp_in) = 0d0  
!             Eyxqpt(el_in,gp_in) = 0d0              
          
          ENDDO
        ENDDO
      
      ENDIF
      
      RETURN
      END SUBROUTINE boundary_edge_land_ldg