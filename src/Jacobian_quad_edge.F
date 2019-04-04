C..... Turns out the jacobian for edge integrals is still length/2 for quads so this code isn't necessary/done.
      subroutine Jacobian_quad_edge(LED,EL,NX,NY,XEGP,JM)

      use global,only : X1,X2,X3,X4,Y1,Y2,Y3,Y4,nm
      
      implicit none
      
      integer, intent(in) LED, N1, N2, EL
      real(8), intent(in) XEGP
      real(8), intent(out) JM
      real(8) J1, J2, J3, J4, Xi, Yi, ksi, eta
      
      N1 = NM(EL,1)
      N2 = NM(EL,2)
      N3 = NM(EL,3)
      N4 = NM(EL,4)
      
      X1 = X(N1)
      X2 = X(N2)
      X3 = X(N3)
      X4 = X(N4)
      Y1 = Y(N1)
      Y2 = Y(N2)
      Y3 = Y(N3)
      Y4 = Y(N4)
      
      Xi = X1 + (X2 - X1)*(XEGP + 1.d0)*abs(NX)*0.5d0
      Yi = Y1 + (Y2 - Y1)*(XEGP + 1.d0)*abs(NX)*0.5d0
      if (LED .eq. 1) then
         eta = -1.d0
         J1 = -.25d0*(1.d0-eta)*X1 + .25d0*(1.d0-eta)*X2 +
     &                         .25d0*(1.d0+eta)*X3 - .25d0*(1.d0+eta)*X4
         J3 = -.25d0*(1.d0-eta)*Y1 + .25d0*(1.d0-eta)*Y2 +
     &                         .25d0*(1.d0+eta)*Y3 - .25d0*(1.d0+eta)*Y4
         JM = sqrt(J1**2.d0 + J3**2.d0)
      elseif (LED .eq. 2) then
         ksi = 1.d0
         eta = XEGP
         
      endif 
      
      end subroutine