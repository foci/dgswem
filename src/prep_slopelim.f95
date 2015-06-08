C***********************************************************************
C     
C     SUBROUTINE PREP_SLOPELIM()
C     
C     This subroutine does preparatory stuff for the slope limiter
C     
C     Written by Ethan Kubatko (07-13-2005)
C     01-10-2011 - cem - adapted for slew of new limiters
C     NOTE: This is a very expensive step and can be saved in a file 
C     for multiple runs; needs to be optimized still
C     
C***********************************************************************

      SUBROUTINE PREP_SLOPELIM()
      
C.....Use appropriate modules

      USE GLOBAL
      USE DG

      IMPLICIT NONE
      
C.....Declare local variables

      INTEGER GED,i,j,ll,k,lll,mm,ell,nin,bbb,bmm
      Real(SZ) areau,xmax,xmin,ymax,ymin
      Real(SZ) dxdxi1,dydxi1,dxdxi2,dydxi2,dxi1dx,dxi2dx
      Real(SZ) dxi1dy,dxi2dy,ell_1,ell_2,ell_3
      Real(SZ) ZEVERTEX2(3),ZEVERTEX(3)

      Real(SZ),Allocatable :: tempmat(:,:),tempInv(:,:),tempag(:,:)
      Real(SZ),Allocatable :: AreaV_integral(:,:,:,:), A(:,:)
      Real(SZ),Allocatable :: Full_M_inv(:,:), temp_p(:,:),temp_t(:,:)
      Real(SZ),Allocatable :: Taylor_mass(:,:,:)

      Allocate ( tempmat(dofh,dofh),tempInv(dofh,dofh), tempag(dofh,dofh) )
      Allocate ( AreaV_integral(MNE,0:ph,0:ph,3), A(dofh,dofh) )
      Allocate ( Full_M_inv(dofh,dofh), temp_p(dofh,dofh) )
      Allocate ( Taylor_mass(MNE,dofh,dofh),temp_t(dofh,dofh))
      
C.....Initialize SL3
      
      SL3 = 0

C.....Loop over the elements

      DO J = 1,MNE
         
C.....Retrieve the nodal coordinates of the given element

         N1 = NM(J,1)
         N2 = NM(J,2)
         N3 = NM(J,3)

C.....Compute the barycenter coordinates of the element and store

         XBC(J) = 1.D0/3.D0*(X(N1) + X(N2) + X(N3))
         YBC(J) = 1.D0/3.D0*(Y(N1) + Y(N2) + Y(N3))

         X1 = XBC(J)
         Y1 = YBC(J)
         
C.....Loop over the edges to fnd the neighboring elements

         DO I = 1,3
            
C.....Retrieve the global edge number of the element for the given edge

            GED = NELED(I,J)

C.....Retrieve the neighboring element number and store into an array

            EL_NBORS(I,J) = NEDEL(1,GED)

            IF (EL_NBORS(I,J).EQ.J) EL_NBORS(I,J) = NEDEL(2,GED)

C.....If the element has an edge that is on boundary go to next element
C     sb-2007/07/27 commented out

C     IF (EL_NBORS(I,J).EQ.0) GOTO 111
            
         ENDDO
         
C.....Set the 4th neighbroing element equal to the first

         EL_NBORS(4,J) = EL_NBORS(1,J)

C.....Now loop over the three edges of the element again

         DO I = 1,3

            IF(EL_NBORS(I,J).EQ.0.OR.EL_NBORS(I+1,J).EQ.0) GOTO 110
            
C.....Compute the barycenter coordinates of two neighboring elements

            N1 = NM(EL_NBORS(I,J),1)
            N2 = NM(EL_NBORS(I,J),2)
            N3 = NM(EL_NBORS(I,J),3)
            
            X2 = 1.D0/3.D0*(X(N1) + X(N2) + X(N3))
            Y2 = 1.D0/3.D0*(Y(N1) + Y(N2) + Y(N3))
            
            N1 = NM(EL_NBORS(I+1,J),1)
            N2 = NM(EL_NBORS(I+1,J),2)
            N3 = NM(EL_NBORS(I+1,J),3)

            X3 = 1.D0/3.D0*(X(N1) + X(N2) + X(N3))
            Y3 = 1.D0/3.D0*(Y(N1) + Y(N2) + Y(N3))
            
C.....Compute the time independent planar constant

            SL3(I,J) = X1*(Y2 - Y3) + X2*(Y3 - Y1) + X3*(Y1 - Y2)

            IF (SL3(I,J).LE.0.AND.SLOPEFLAG.NE.0) then
           WRITE(16,*) 'WARNING. SL3(',I,',',J,') =',SL3(I,J),' <= 0.',
     &              '    ELEMENT ',J,
     &              ' WILL NOT BE CONSIDERED IN SLOPE LIMITING.'
               WRITE(16,*)
            WRITE(*,*) 'WARNING. SL3(',I,',',J,') =',SL3(I,J),' <= 0.',
     &              '    ELEMENT ',J,
     &              ' WILL NOT BE CONSIDERED IN SLOPE LIMITING.'
               WRITE(*,*)
               SL3(I,J) = 0.D0
            ENDIF
            
 110        CONTINUE

         ENDDO
         
      ENDDO


C******************************************************************************
C.....Vertex-based slope limiter (need the following stuff for integration)
C.....and must fill array for all possible p (ie. dofl:dofh)

      
#ifdef SLOPEALL

      XBCb = 0.D0
      YBCb = 0.D0
      Deltx = 0.D0
      Delty = 0.D0
      xtransform = 0.D0
      ytransform = 0.D0
      xi1BCb = 0.D0
      xi2BCb = 0.D0
      xi1vert = 0.D0
      xi2vert = 0.D0
      xtransformv = 0.D0
      ytransformv = 0.D0
      YBCv = 0.D0
      XBCv = 0.D0
      xi1BCv = 0.D0
      xi2BCv = 0.D0
      f = 0.D0
      g0 = 0.D0 
      varsigma0 = 0.D0
      fv = 0.D0
      g0v = 0.D0 
      varsigma0v = 0.D0
      var2sigmag = 0.D0
      NmatrixInv = 0.D0
      Nmatrix = 0.D0


C.....Loop over p, and allocate for each order

      do ll = 1,ph

         if (ll.gt.0) then
            
C.....Loop over the elements

            do k = 1,MNE

C.....Set areas over the physical elements
               
               areau = 0.5D0*AREAS(k)

C.....Retrieve the nodal coordinates (above) of the given element

               N1 = NM(k,1)
               N2 = NM(k,2)
               N3 = NM(k,3)

                                !areat = 0.5*(x(n2)*y(n3)-x(n3)*y(n2)+x(n3)*y(n1)-x(n1)*y(n3) + x(n1)*y(n2)-x(n2)*y(n1))
                                !areat = areau

C.....Find cell conditioners

               xmax = max( x(n1),x(n2),x(n3) )
               xmin = min( x(n1),x(n2),x(n3) )
               ymax = max( y(n1),y(n2),y(n3) )
               ymin = min( y(n1),y(n2),y(n3) )

               if (ll.le.2) then

                  Deltx(k) = ( xmax - xmin ) / 2.D0 
                  Delty(k) = ( ymax - ymin ) / 2.D0 

               else 

                  Deltx(k) = ( xmax - xmin ) / ll 
                  Delty(k) = ( ymax - ymin ) / ll

               endif

               

C.....Compute the centroid coordinates of the base element in physical space

               XBCb(k) = ( x(n1) + x(n2) + x(n3) ) / 3.D0
               YBCb(k) = ( y(n1) + y(n2) + y(n3) ) / 3.D0

C.....Transform quad points to physical space (for integration) xi --> x

               do mm=1,nagp(ll)

                  ell_1 = -0.5D0 * ( xagp(mm,ph) + yagp(mm,ph) )
                  ell_2 =  0.5D0 * ( xagp(mm,ph) + 1.D0 )
                  ell_3 =  0.5D0 * ( yagp(mm,ph) + 1.D0 )

                  xtransform(k,mm) = x(n1)*ell_1 + x(n2)*ell_2 + x(n3)*ell_3
                  ytransform(k,mm) = y(n1)*ell_1 + y(n2)*ell_2 + y(n3)*ell_3

               enddo

C.....Find centroid coordinates in the master element frame

               xi1BCb(k) =  ( (y(N3)-y(N1))*( XBCb(k) -0.5D0 * 
     &              (x(N2) + x(N3))) + (x(N1) - x(N3))*(YBCb(k)
     &              - 0.5D0*(y(N2) + y(N3)) ) ) / areau
               xi2BCb(k) =  ( (y(N1)-y(N2))*( XBCb(k) -0.5D0 * (x(N2) 
     &              + x(N3))) + (x(N2) - x(N1))*(YBCb(k)
     &              - 0.5D0*(y(N2) + y(N3)) ) ) / areau

C.....Find vertices in the master element frame

               do lll=1,3

                  xi1vert(k,lll) =  ( (y(N3)-y(N1))*( x(nm(k,lll)) -0.5D0 * 
     &                 (x(N2) + x(N3))) + (x(N1) - x(N3))* (y(nm(k,lll))
     &                 - 0.5D0*(y(N2) + y(N3)) ) ) / areau
                  xi2vert(k,lll) =  ( (y(N1)-y(N2))*( x(nm(k,lll)) -0.5D0 * 
     &                 (x(N2) + x(N3))) + (x(N2) - x(N1))* (y(nm(k,lll))
     &                 - 0.5D0*(y(N2) + y(N3)) ) ) / areau

               enddo

C.....Find all neighbors of shared vertex

               do mm =1,MNE     !number of elements

                  do lll =1,3   !number of vertices

                     do nin =1,3 !number of vertices

                        if( NM(k,lll).eq.NM(mm,nin).and.k.ne.mm ) then !find common vertices of "nearby" elements
                           
C.....Compute the centroids of all conterminous (of codimension 2) elements (by vertex) of base element k in physical space

                           XBCv(k,mm) = 1.D0/3.D0*( X(NM(mm,1)) + 
     &                          X(NM(mm,2)) + X(NM(mm,3)) )
                           YBCv(k,mm) = 1.D0/3.D0*( Y(NM(mm,1)) + 
     &                          Y(NM(mm,2)) + Y(NM(mm,3)) )

                                !Stored_neighbors(k,mm,lll) = k*mm*lll !store neighboring elements

C.....Convert the centroid coordinates of all conterminous (of codimension 2) elements (by vertex) 
C.....of base element k to the master element space

                           xi1BCv(k,mm) = ( (y(NM(mm,3))-y(NM(mm,1)))*
     &                          ( XBCv(k,mm) 
     &                          -0.5D0 * (x(Nm(mm,2)) + x(NM(mm,3)))) + 
     &                          (x(NM(mm,1)) - x(NM(mm,3)))*(YBCv(k,mm) 
     &                          - 0.5D0*(y(NM(mm,2)) + y(NM(mm,3))) ) ) / areau
                           xi2BCv(k,mm) = ( (y(NM(mm,1))-y(NM(mm,2)))*
     &                          ( XBCv(k,mm) 
     &                          -0.5D0 * (x(Nm(mm,2)) + x(NM(mm,3)))) + 
     &                          (x(NM(mm,2)) - x(NM(mm,1)))*(YBCv(k,mm)
     &                          - 0.5D0*(y(NM(mm,2)) + y(NM(mm,3))) ) ) / areau
                           
                        endif
                     enddo

                  enddo

               enddo


C.....Now compute the derivatives of the Taylor basis with respect to the physical 
C.....basis using the transformation rules from the paper (e.g. Leibniz and Faa' di Bruno formulas)
               
               dxdxi1 = 0.5D0 * ( x(N2) - x(N1) )
               dydxi1 = 0.5D0 * ( y(N2) - y(N1) )
               dxdxi2 = 0.5D0 * ( x(N3) - x(N1) )
               dydxi2 = 0.5D0 * ( y(N3) - y(N1) )

               dxi1dx = ( y(N3) - y(N1) ) / areau
               dxi2dx = ( y(N1) - y(N2) ) / areau
               dxi1dy = ( x(N1) - x(N3) ) / areau 
               dxi2dy = ( x(N2) - x(N1) ) / areau

C.....Write the generalized Taylor basis of order p in physical 
C.....coordinates (x(xi1,xi2), y(xi1,xi2)) and integrate over elements
C.....using the physical to master transformation, e.g. T^-1:x-->xi

               do i = 0,ll      !max polynomial degree in x

                  do j = 0,ll   !max polynomial degree in y

                     Call factorial(i,fact(i))
                     Call factorial(j,fact(j))

                     Area_integral(k,i,j) = 0.D0


                     do mm = 1,NAGP(ll) !number of quad points

                        Area_integral(k,i,j) = Area_integral(k,i,j) + 
     &                       ( (  xtransform(k,mm) - XBCb(k) )**i 
     &                       * (  ytransform(k,mm)- YBCb(k))**j * ( wagp(mm,ll) ) )
     &                       * abs( dxdxi1*dydxi2 - dxdxi2*dydxi1 )
     &                       / ( fact(i)*fact(j)*Deltx(k)**i * Delty(k)**j )		


                     enddo
                     
                     do mm =1,nagp(ll) !at quad points

                        f(k,mm,i,j) = (  xtransform(k,mm) - XBCb(k)  )**i 
     &                       / ( fact(i) * Deltx(k)**i )
                        g0(k,mm,i,j) = (  ytransform(k,mm) - YBCb(k)  )**j 
     &                       / ( fact(j) * Delty(k)**j )		

                        if (i.eq.0.and.j.eq.0) then

                           varsigma0(k,mm,i,j) = 1

                        else
                           
                           varsigma0(k,mm,i,j) = ( f(k,mm,i,j) * g0(k,mm,i,j) ) 
     &                          - Area_integral(k,i,j)/areau

                        endif

                     enddo


                     do lll = 1,3 !at vertices

                        AreaV_integral(k,i,j,lll) = 0.D0

                        do mm = 1,nagp(ll) !number of quad points

                           AreaV_integral(k,i,j,lll) = AreaV_integral(k,i,j,lll) + 
     &                          ( (  x(nm(k,lll)) - XBCb(k) )**i 
     &                          * (  y(nm(k,lll)) - YBCb(k) )**j * ( wagp(mm,ll) ) )
     &                          * abs( dxdxi1*dydxi2 - dxdxi2*dydxi1 )
     &                          / ( fact(i)*fact(j)*Deltx(k)**i * Delty(k)**j )

                        enddo

                        fv(k,lll,i,j) = (  x(nm(k,lll)) - XBCb(k) )**i / 
     &                       ( fact(i) * Deltx(k)**i )
                        
                        g0v(k,lll,i,j) =  ( y(nm(k,lll)) - YBCb(k) )**j / 
     &                       ( fact(j) * Delty(k)**j )
                        
                        if (i.eq.0.and.j.eq.0) then

                           varsigma0v(k,lll,i,j) = 1

                        else
                           
                           varsigma0v(k,lll,i,j) = ( fv(k,lll,i,j) * 
     &                          g0v(k,lll,i,j) ) 
     &                          - Area_integral(k,i,j)/areau

                        endif

                     enddo

                  enddo

               enddo



C.....Re-order the Taylor basis functions and componentwise derivatives into hierarchical order

               bbb = 1 
               do j = 0,ll

                  do i = 0,j

                     do mm = 1,NAGP(ll)


                        var2sigmag(k,mm,bbb) = varsigma0(k,mm,i,j-i)

                        
                        if ( abs(var2sigmag(k,mm,bbb)).lt.1.0E-15 ) then 

                           var2sigmag(k,mm,bbb) = 0.D0 
                           
                        endif


                     enddo      !mm

                     do lll = 1,3

                        var2sigmav(k,lll,bbb) = varsigma0v(k,lll,i,j-i)

                        if ( abs(var2sigmav(k,lll,bbb)).lt.1.0E-15 ) then 

                                !var2sigmav(k,lll,bbb) = 0.D0 
                           
                        endif	

                     enddo      !lll

                     
                     bi(bbb) = i
                     bj(bbb) = j

                     bbb = bbb + 1

                  enddo         !i

               enddo            !j
               


C.....Compute the inner product matrix Pmatrix, of the Taylor 
C.....basis with the Dubiner basis, and compute the transformation 
C.....matrix Nmatrix=Pmatrix*M(-1), using the mass matrix inverse M_inv,

               ell = (ll+1)*(ll+2)/2

               A = 0.D0
               Full_M_inv = 0.D0
               temp_p = 0.D0



               do i = 1,ell

                  do j = 1,ell
                     
                     pmatrix(k,i,j)=0.D0

                     do mm=1,NAGP(ll)
                        
                        pmatrix(k,i,j) = pmatrix(k,i,j)  + wagp(mm,ll) * 
     &                       var2sigmag(k,mm,i) * phi_area(j,mm,ll) 
                        

                        temp_p(i,j) = pmatrix(k,i,j)

                        Taylor_mass(k,i,j) = Taylor_mass(k,i,j) +  wagp(mm,ll) 
     &                       *  var2sigmag(k,mm,i) *  var2sigmag(k,mm,j)

                        temp_t(i,j) = Taylor_mass(k,i,j)

                     enddo
                     
                  enddo
                  
                  Full_M_inv(i,i) =   M_INV(i,ll)
                  
               enddo

               Call Inv(temp_t(1:ell,1:ell), tempInv(1:ell,1:ell), ell)
               
               Nmatrix(k,1:ell,1:ell,ell) = matmul(tempInv(1:ell,1:ell), temp_p(1:ell,1:ell))

               tempmat(1:ell,1:ell) = Nmatrix(k,1:ell,1:ell,ell)

C.....Invert the transformation matrix Nmatrix^(-1), 

               Call Inv(tempmat(1:ell,1:ell), tempInv(1:ell,1:ell), ell) 
               
               NmatrixInv(k,1:ell,1:ell,ell) = tempInv(1:ell,1:ell)

            enddo               !k-elements
            
         endif                  !ll loop

      enddo                     !p_adapt

C.....Construct focal neighbors for non-vertex based limiters
      
      focal_neigh = 0
      focal_up = 0

      do j = 1,mne

         bmm = 1

         do lll = 1,3

            do ell = 1,nneigh_elem(nm(j,lll))

               focal_neigh(j,bmm) = neigh_elem(nm(j,lll), ell)  

               focal_up(j) = bmm

               bmm = bmm + 1

            enddo

         enddo

      enddo

#endif

      RETURN
      END SUBROUTINE PREP_SLOPELIM

C.....Need factorial function for generalization

#ifdef SLOPEALL

      subroutine factorial(n,p)
      
      implicit none
      integer n,p,i

      p = 1

      do i = 1, n

         p = p * i

      enddo
      

      end subroutine factorial

#endif

C.....Subroutine to find the inverse of a square matrix by Guass-Jordan elimination

#ifdef SLOPEALL

      subroutine Inv(matrix, inverse, n)
      Use sizes, only : sz

      implicit none
      integer n
      real(sz), dimension(n,n) :: matrix
      real(sz), dimension(n,n) :: inverse
      
      integer :: i, j, k, l
      real(sz) :: m
      real(sz), dimension(n,2*n) :: augmatrix !augmented matrix
      
                                !Augment input matrix with an identity matrix

      do i = 1, n

         do j = 1, 2*n

            if ( j.le.n ) then

               augmatrix(i,j) = matrix(i,j)

            else if ((i+n) == j) then

               augmatrix(i,j) = 1

            else

               augmatrix(i,j) = 0

            endif

         enddo

      enddo
      
                                !Reduce augmented matrix to upper traingular form

      do k =1, n-1

         if (augmatrix(k,k) == 0) then


            do i = k+1, n

               if (augmatrix(i,k) /= 0) then

                  do j = 1,2*n

                     augmatrix(k,j) = augmatrix(k,j)+augmatrix(i,j)

                  enddo

               endif


            enddo

         endif

         do j = k+1, n			

            m = augmatrix(j,k)/augmatrix(k,k)

            do i = k, 2*n

               augmatrix(j,i) = augmatrix(j,i) - m*augmatrix(k,i)

            enddo

         enddo

      enddo 
      
                                !Test for invertibility

      do i = 1, n

         if (augmatrix(i,i) == 0) then

            inverse = 0

            return

         endif

      enddo
      
                                !Make diagonal elements as 1

      do i = 1 , n

         m = augmatrix(i,i)

         do j = i, (2 * n)				

            augmatrix(i,j) = (augmatrix(i,j) / m)

         enddo

      enddo
      
                                !Reduced right side half of augmented matrix to identity matrix

      do k = n-1, 1, -1

         do i =1, k

            m = augmatrix(i,k+1)

            do j = k, (2*n)

               augmatrix(i,j) = augmatrix(i,j) -augmatrix(k+1,j) * m

            enddo

         enddo

      enddo				
      
                                !Compute answer

      do i =1, n

         do j = 1, n

            inverse(i,j) = augmatrix(i,j+n)

         enddo

      enddo

      end subroutine Inv

#endif

