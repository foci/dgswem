
      SUBROUTINE interior_edge_nflux(sed,eed,nqpte)

      USE globals, ONLY:  const,inx,iny,icfac,g, &
                          Zi,Ze,Hi,He,Qxi,Qxe,Qyi,Qye, &
                          xmi,xme,ymi,yme,xymi,xyme, &
                          Zhatv,Qxhatv,Qyhatv, &
                          detJe_in,detJe_ex, &
                          Exxi,Exxe,Eyyi,Eyye,Exyi,Exye,Eyxi,Eyxe
   
      IMPLICIT NONE

      INTEGER :: sed,eed,nqpte
      INTEGER :: pt,ed
                  
        
ed_points: DO pt = 1,nqpte ! Compute numerical fluxes for all edges
              
!DIR$ SIMD              
!!DIR$ VECTOR ALIGNED
              DO ed = sed,eed
                
                const(ed) = max(abs(Qxi(ed,pt)%ptr*inx(ed,pt) + Qyi(ed,pt)%ptr*iny(ed,pt))/Hi(ed,pt)%ptr + sqrt(g*Hi(ed,pt)%ptr*icfac(ed,pt)), &
                                abs(Qxe(ed,pt)%ptr*inx(ed,pt) + Qye(ed,pt)%ptr*iny(ed,pt))/He(ed,pt)%ptr + sqrt(g*He(ed,pt)%ptr*icfac(ed,pt)))
              ENDDO
                                        
!DIR$ SIMD              
!!DIR$ VECTOR ALIGNED
              DO ed = sed,eed                  
                Zhatv(ed) = .5d0*(inx(ed,pt)*(Qxi(ed,pt)%ptr + Qxe(ed,pt)%ptr) + iny(ed,pt)*(Qyi(ed,pt)%ptr + Qye(ed,pt)%ptr) &
                                        - const(ed)*(Ze(ed,pt)%ptr - Zi(ed,pt)%ptr))           
              ENDDO
              
!DIR$ SIMD
!!DIR$ VECTOR ALIGNED
              DO ed = sed,eed     
                Qxhatv(ed) = .5d0*(inx(ed,pt)*(xmi(ed,pt)%ptr + xme(ed,pt)%ptr) + iny(ed,pt)*(xymi(ed,pt)%ptr + xyme(ed,pt)%ptr)  &
                                        - const(ed)*(Qxe(ed,pt)%ptr - Qxi(ed,pt)%ptr))                                     
              ENDDO
              
!DIR$ SIMD
!!DIR$ VECTOR ALIGNED
              DO ed = sed,eed
                Qyhatv(ed) = .5d0*(inx(ed,pt)*(xymi(ed,pt)%ptr + xyme(ed,pt)%ptr) + iny(ed,pt)*(ymi(ed,pt)%ptr + yme(ed,pt)%ptr)  &
                                        - const(ed)*(Qye(ed,pt)%ptr - Qyi(ed,pt)%ptr))                                     
              ENDDO
              
              DO ed = sed,eed
                Qxhatv(ed) = Qxhatv(ed) - inx(ed,pt)*.5d0*(Exxe(ed,pt)%ptr+Exxi(ed,pt)%ptr) - iny(ed,pt)*.5d0*(Exye(ed,pt)%ptr+Exyi(ed,pt)%ptr)               
              ENDDO
              
              DO ed = sed,eed
!                 Qyhatv(ed) = Qyhatv(ed) - inx(ed,pt)*.5d0*(Eyxe(ed,pt)%ptr+Eyxi(ed,pt)%ptr) - iny(ed,pt)*.5d0*(Eyye(ed,pt)%ptr+Eyyi(ed,pt)%ptr)
                Qyhatv(ed) = Qyhatv(ed) - inx(ed,pt)*.5d0*(Exye(ed,pt)%ptr+Exyi(ed,pt)%ptr) - iny(ed,pt)*.5d0*(Eyye(ed,pt)%ptr+Eyyi(ed,pt)%ptr)                
              ENDDO
              
              
!DIR$ IVDEP              
              DO ed = sed,eed
                Ze(ed,pt)%ptr = -detJe_ex(ed,pt)*Zhatv(ed)
                Zi(ed,pt)%ptr =  detJe_in(ed,pt)*Zhatv(ed)                
              ENDDO          
!DIR$ IVDEP              
              DO ed = sed,eed             
!                 Qxe(ed,pt)%ptr = -detJe_ex(ed,pt)*(Qxhatv(ed) - inx(ed,pt)*Exxe(ed,pt)%ptr - iny(ed,pt)*Exye(ed,pt)%ptr) 
!                 Qxi(ed,pt)%ptr =  detJe_in(ed,pt)*(Qxhatv(ed) - inx(ed,pt)*Exxi(ed,pt)%ptr - iny(ed,pt)*Exyi(ed,pt)%ptr)               

                Qxe(ed,pt)%ptr = -detJe_ex(ed,pt)*Qxhatv(ed)
                Qxi(ed,pt)%ptr =  detJe_in(ed,pt)*Qxhatv(ed)
              ENDDO   
!DIR$ IVDEP              
              DO ed = sed,eed                      
!                 Qye(ed,pt)%ptr = -detJe_ex(ed,pt)*(Qyhatv(ed) - inx(ed,pt)*Exye(ed,pt)%ptr - iny(ed,pt)*Eyye(ed,pt)%ptr) 
!                 Qyi(ed,pt)%ptr =  detJe_in(ed,pt)*(Qyhatv(ed) - inx(ed,pt)*Exyi(ed,pt)%ptr - iny(ed,pt)*Eyyi(ed,pt)%ptr)               

                Qye(ed,pt)%ptr = -detJe_ex(ed,pt)*Qyhatv(ed)
                Qyi(ed,pt)%ptr =  detJe_in(ed,pt)*Qyhatv(ed)
              ENDDO               

        ENDDO ed_points  


       RETURN 
       END SUBROUTINE interior_edge_nflux
       
       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  


      SUBROUTINE interior_edge_nflux_ldg(sed,eed,nqpte)  
      
      USE globals, ONLY: Exxi,Exxe,Eyyi,Eyye,Exyi,Exye,Eyxi,Eyxe, &
                         Qxi,Qxe,Qyi,Qye, &
                         inx,iny, &
                         detJe_in,detJe_ex
      
      IMPLICIT NONE
      
      INTEGER :: sed,eed
      INTEGER :: nqpte
      
      INTEGER :: pt,ed
                  
      DO pt = 1,nqpte
        DO ed = sed,eed

!           Exxe(ed,pt)%ptr = -detJe_ex(ed,pt)*inx(ed,pt)*Qxi(ed,pt)%ptr        
!           Exxi(ed,pt)%ptr =  detJe_in(ed,pt)*inx(ed,pt)*Qxe(ed,pt)%ptr
! 
!           Eyye(ed,pt)%ptr = -detJe_ex(ed,pt)*iny(ed,pt)*Qyi(ed,pt)%ptr          
!           Eyyi(ed,pt)%ptr =  detJe_in(ed,pt)*iny(ed,pt)*Qye(ed,pt)%ptr
! 
!           Exye(ed,pt)%ptr = -detJe_ex(ed,pt)*(iny(ed,pt)*Qxi(ed,pt)%ptr + inx(ed,pt)*Qyi(ed,pt)%ptr)          
!           Exyi(ed,pt)%ptr =  detJe_in(ed,pt)*(iny(ed,pt)*Qxe(ed,pt)%ptr + inx(ed,pt)*Qye(ed,pt)%ptr)
! 

          Exxe(ed,pt)%ptr = -detJe_ex(ed,pt)*inx(ed,pt)*.5d0*(Qxi(ed,pt)%ptr+Qxe(ed,pt)%ptr)
          Exxi(ed,pt)%ptr =  detJe_in(ed,pt)*inx(ed,pt)*.5d0*(Qxe(ed,pt)%ptr+Qxi(ed,pt)%ptr)

          Eyye(ed,pt)%ptr = -detJe_ex(ed,pt)*iny(ed,pt)*.5d0*(Qyi(ed,pt)%ptr+Qye(ed,pt)%ptr)           
          Eyyi(ed,pt)%ptr =  detJe_in(ed,pt)*iny(ed,pt)*.5d0*(Qye(ed,pt)%ptr+Qyi(ed,pt)%ptr)
          
          Exye(ed,pt)%ptr = -detJe_ex(ed,pt)*.5d0*(iny(ed,pt)*(Qxi(ed,pt)%ptr+Qxe(ed,pt)%ptr)+inx(ed,pt)*(Qyi(ed,pt)%ptr+Qye(ed,pt)%ptr))
          Exyi(ed,pt)%ptr =  detJe_in(ed,pt)*.5d0*(iny(ed,pt)*(Qxe(ed,pt)%ptr+Qxi(ed,pt)%ptr)+inx(ed,pt)*(Qye(ed,pt)%ptr+Qyi(ed,pt)%ptr))          
          
!           Exye(ed,pt)%ptr = -detJe_ex(ed,pt)*iny(ed,pt)*.5d0*(Qxi(ed,pt)%ptr+Qxe(ed,pt)%ptr)
!           Exyi(ed,pt)%ptr =  detJe_in(ed,pt)*iny(ed,pt)*.5d0*(Qxe(ed,pt)%ptr+Qxi(ed,pt)%ptr)
! 
!           Eyxe(ed,pt)%ptr = -detJe_ex(ed,pt)*inx(ed,pt)*.5d0*(Qyi(ed,pt)%ptr+Qye(ed,pt)%ptr)            
!           Eyxi(ed,pt)%ptr =  detJe_in(ed,pt)*inx(ed,pt)*.5d0*(Qye(ed,pt)%ptr+Qyi(ed,pt)%ptr)
 
          
        ENDDO
      ENDDO
      
      RETURN
      END SUBROUTINE interior_edge_nflux_ldg
                            
                            
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
       
       
      SUBROUTINE recieve_edge_nflux(nred,nqpte)

      USE globals, ONLY: const,g,pt5g,recipHa, &
                         Zhatv,Qxhatv,Qyhatv
         
      USE messenger2, ONLY: Zri,Zre,Hri,Hre,Qxri,Qxre,Qyri,Qyre, &
                            xmri,xmre,ymri,ymre,xymri,xymre, &
                            detJe_recv,rcfac,hbr,rnx,rny, &
                            Exxri,Exxre,Eyyri,Eyyre,Exyri,Exyre,Eyxri,Eyxre
   
      IMPLICIT NONE

      INTEGER :: nred,nqpte
      INTEGER :: pt,ed
                  
        
      DO pt = 1,nqpte
      
        DO ed = 1,nred
          Hre(ed) = Zre(ed,pt)%ptr + hbr(ed,pt)
        ENDDO
      
!!DIR$ VECTOR ALIGNED      
        DO ed = 1,nred
          const(ed) = max(abs(Qxri(ed,pt)%ptr*rnx(ed,pt) + Qyri(ed,pt)%ptr*rny(ed,pt))/Hri(ed,pt)%ptr + sqrt(g*Hri(ed,pt)%ptr*rcfac(ed,pt)), &
                          abs(Qxre(ed,pt)%ptr*rnx(ed,pt) + Qyre(ed,pt)%ptr*rny(ed,pt))/Hre(ed)        + sqrt(g*Hre(ed)*rcfac(ed,pt)))          
        ENDDO
        

!DIR$ IVDEP
!!DIR$ VECTOR ALIGNED
        DO ed = 1,nred      
          Zhatv(ed) = .5d0*(rnx(ed,pt)*(Qxri(ed,pt)%ptr + Qxre(ed,pt)%ptr) + rny(ed,pt)*(Qyri(ed,pt)%ptr + Qyre(ed,pt)%ptr) &
                                        - const(ed)*(Zre(ed,pt)%ptr - Zri(ed,pt)%ptr))
        ENDDO

!DIR$ IVDEP
!!DIR$ VECTOR ALIGNED        
        DO ed = 1,nred
          recipHa(ed) = 1d0/Hre(ed)
          
          xmre(ed) = pt5g*(Hre(ed)*Hre(ed) - hbr(ed,pt)*hbr(ed,pt)) + Qxre(ed,pt)%ptr*Qxre(ed,pt)%ptr*recipHa(ed)
          ymre(ed) = pt5g*(Hre(ed)*Hre(ed) - hbr(ed,pt)*hbr(ed,pt)) + Qyre(ed,pt)%ptr*Qyre(ed,pt)%ptr*recipHa(ed)
          xymre(ed) = Qxre(ed,pt)%ptr*Qyre(ed,pt)%ptr*recipHa(ed)
        ENDDO
 
!DIR$ IVDEP
!!DIR$ VECTOR ALIGNED 
        DO ed = 1,nred
          Qxhatv(ed) = .5d0*(rnx(ed,pt)*(xmri(ed,pt)%ptr + xmre(ed)) + rny(ed,pt)*(xymri(ed,pt)%ptr + xymre(ed))  &
                                        - const(ed)*(Qxre(ed,pt)%ptr - Qxri(ed,pt)%ptr))
        ENDDO
  
!DIR$ IVDEP
!!DIR$ VECTOR ALIGNED  
        DO ed = 1,nred
          Qyhatv(ed) = .5d0*(rnx(ed,pt)*(xymri(ed,pt)%ptr + xymre(ed)) + rny(ed,pt)*(ymri(ed,pt)%ptr + ymre(ed))  &
                                        - const(ed)*(Qyre(ed,pt)%ptr - Qyri(ed,pt)%ptr))
        ENDDO
        
        DO ed = 1,nred
          Qxhatv(ed) = Qxhatv(ed) - rnx(ed,pt)*.5d0*(Exxri(ed,pt)%ptr+Exxre(ed,pt)%ptr) - rny(ed,pt)*.5d0*(Exyri(ed,pt)%ptr+Exyre(ed,pt)%ptr)
        ENDDO
        
        DO ed = 1,nred
!           Qyhatv(ed) = Qyhatv(ed) - rnx(ed,pt)*.5d0*(Eyxri(ed,pt)%ptr+Eyxre(ed,pt)%ptr) - rny(ed,pt)*.5d0*(Eyyri(ed,pt)%ptr+Eyyre(ed,pt)%ptr)
          Qyhatv(ed) = Qyhatv(ed) - rnx(ed,pt)*.5d0*(Exyri(ed,pt)%ptr+Exyre(ed,pt)%ptr) - rny(ed,pt)*.5d0*(Eyyri(ed,pt)%ptr+Eyyre(ed,pt)%ptr)          
        ENDDO

        
        
        DO ed = 1,nred                                       
          Zri(ed,pt)%ptr =  detJe_recv(ed,pt)*Zhatv(ed)          
        ENDDO   
        DO ed = 1,nred                                         
!           Qxri(ed,pt)%ptr =  detJe_recv(ed,pt)*(Qxhatv(ed) - rnx(ed,pt)*Exxri(ed,pt)%ptr - rny(ed,pt)*Exyri(ed,pt)%ptr)      
          Qxri(ed,pt)%ptr =  detJe_recv(ed,pt)*Qxhatv(ed)
        ENDDO        
        DO ed = 1,nred 
!           Qyri(ed,pt)%ptr =  detJe_recv(ed,pt)*(Qyhatv(ed) - rnx(ed,pt)*Exyri(ed,pt)%ptr - rny(ed,pt)*Eyyri(ed,pt)%ptr)
          Qyri(ed,pt)%ptr =  detJe_recv(ed,pt)*Qyhatv(ed)
        ENDDO
      ENDDO

      RETURN 
      END SUBROUTINE recieve_edge_nflux
       
       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      SUBROUTINE recieve_edge_nflux_ldg(nred,nqpte)

      USE messenger2, ONLY: Qxre,Qyre,Qxri,Qyri, &
                            Exxri,Eyyri,Exyri,Eyxri, &
                            detJe_recv,rnx,rny
      
      IMPLICIT NONE
      
      INTEGER :: nred,nqpte
      
      INTEGER :: pt,ed
      
      DO pt = 1,nqpte
        DO ed = 1,nred
        
!           Exxri(ed,pt)%ptr = detJe_recv(ed,pt)*rnx(ed,pt)*Qxre(ed,pt)%ptr
!           Eyyri(ed,pt)%ptr = detJe_recv(ed,pt)*rny(ed,pt)*Qyre(ed,pt)%ptr
!           Exyri(ed,pt)%ptr = detJe_recv(ed,pt)*(rny(ed,pt)*Qxre(ed,pt)%ptr + rnx(ed,pt)*Qyre(ed,pt)%ptr)

          Exxri(ed,pt)%ptr = detJe_recv(ed,pt)*rnx(ed,pt)*.5d0*(Qxre(ed,pt)%ptr+Qxri(ed,pt)%ptr)
          Eyyri(ed,pt)%ptr = detJe_recv(ed,pt)*rny(ed,pt)*.5d0*(Qyre(ed,pt)%ptr+Qyri(ed,pt)%ptr)
          Exyri(ed,pt)%ptr = detJe_recv(ed,pt)*.5d0*(rny(ed,pt)*(Qxre(ed,pt)%ptr+Qxri(ed,pt)%ptr)+rnx(ed,pt)*(Qyre(ed,pt)%ptr+Qyri(ed,pt)%ptr))
!           Exyri(ed,pt)%ptr = detJe_recv(ed,pt)*rny(ed,pt)*.5d0*(Qxre(ed,pt)%ptr+Qxri(ed,pt)%ptr)
!           Eyxri(ed,pt)%ptr = detJe_recv(ed,pt)*rnx(ed,pt)*.5d0*(Qyre(ed,pt)%ptr+Qyri(ed,pt)%ptr)
          
        ENDDO
      ENDDO
                  
      RETURN
      END SUBROUTINE recieve_edge_nflux_ldg

      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      SUBROUTINE numerical_flux(Qx_in,Qy_in,Z_in,Qx_ex,Qy_ex,Z_ex,hb,nx,ny,sp,Qxhat,Qyhat,Zhat)
      
!       USE globals, ONLY:rp,nx,ny,g,eig_in,eig_ex,alpha, &
!                         H_in,H_ex,Qx_in,Qx_ex,Qy_in,Qy_ex, &
!                         xmom_in,xmom_ex,ymom_in,ymom_ex,xymom_in,xymom_ex, &
!                         Hhat,Qxhat,Qyhat

      USE globals, ONLY:rp,g


      IMPLICIT NONE
      
      REAL(rp) :: eig_in,eig_ex,alpha
      REAL(rp) :: Qx_in,Qy_in,Z_in
      REAL(rp) :: Qx_ex,Qy_ex,Z_ex
      REAL(rp) :: nx,ny,sp,hb
      REAL(rp) :: Qxhat,Qyhat,Zhat
      REAL(rp) :: H_in,H_ex,hb_in,hb_ex
      REAL(rp) :: xmom_in,ymom_in,xymom_in
      REAL(rp) :: xmom_ex,ymom_ex,xymom_ex
      REAL(rp) :: rHin,rHex
      REAL(rp) :: nxsp,cfac
      
      hb_in = hb
      hb_ex = hb
      
      H_in = Z_in + hb_in
      H_ex = Z_ex + hb_ex
            
      rHin = 1d0/H_in
      rHex = 1d0/H_ex
      
      nxsp = nx*sp
      cfac = ny**2 + (nxsp)**2
      
      eig_in = abs(Qx_in*nxsp + Qy_in*ny)*rHin + sqrt(g*H_in*cfac) 
      eig_ex = abs(Qx_ex*nxsp + Qy_ex*ny)*rHex + sqrt(g*H_ex*cfac)

      alpha = max(eig_in,eig_ex)

!       eig(1) = abs((Qx_in*nxsp + Qy_in*ny)/H_in + sqrt(g*H_in*cfac))
!       eig(2) = abs((Qx_in*nxsp + Qy_in*ny)/H_in)
!       eig(3) = abs((Qx_in*nxsp + Qy_in*ny)/H_in - sqrt(g*H_in*cfac))
! 
!       eig(4) = abs((Qx_ex*nxsp + Qy_ex*ny)/H_ex + sqrt(g*H_ex*cfac))
!       eig(5) = abs((Qx_ex*nxsp + Qy_ex*ny)/H_ex)
!       eig(6) = abs((Qx_ex*nxsp + Qy_ex*ny)/H_ex - sqrt(g*H_ex*cfac))
! 
!       alpha = maxval(eig,1)

      xmom_in = Qx_in*Qx_in*rHin + .5d0*g*(H_in*H_in - hb_in*hb_in)
      xmom_ex = Qx_ex*Qx_ex*rHex + .5d0*g*(H_ex*H_ex - hb_ex*hb_ex)
      
      ymom_in = Qy_in*Qy_in*rHin + .5d0*g*(H_in*H_in - hb_in*hb_in)
      ymom_ex = Qy_ex*Qy_ex*rHex + .5d0*g*(H_ex*H_ex - hb_ex*hb_ex)
       
      xymom_in = Qx_in*Qy_in*rHin
      xymom_ex = Qx_ex*Qy_ex*rHex

      Zhat  = .5d0*( nxsp*(Qx_in+Qx_ex) + ny*(Qy_in+Qy_ex) - alpha*(Z_ex-Z_in) )
      Qxhat = .5d0*( nxsp*(xmom_in+xmom_ex) + ny*(xymom_in+xymom_ex) - alpha*(Qx_ex-Qx_in) )
      Qyhat = .5d0*( nxsp*(xymom_in+xymom_ex) + ny*(ymom_in+ymom_ex) - alpha*(Qy_ex-Qy_in) )

      RETURN
      END SUBROUTINE numerical_flux