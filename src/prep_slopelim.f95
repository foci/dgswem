!***********************************************************************
!     
!     SUBROUTINE PREP_SLOPELIM()
!     
!     This subroutine does preparatory stuff for the slope limiter
!     
!     Written by Ethan Kubatko (07-13-2005)
!     01-10-2011 - cem - adapted for slew of new limiters
!     NOTE: This is a very expensive step and can be saved in a file 
!     for multiple runs; dg%needs to be optimized still
!     
!***********************************************************************

      SUBROUTINE PREP_SLOPELIM(s)
      
!.....Use appropriate modules
      use sizes
      USE GLOBAL
      USE DG

      IMPLICIT NONE
      
      type (sizes_type) :: s

!.....Declare local variables

      INTEGER GED,i,j,ll,k,lll,mm,ell,nin,bbb,bmm
      Real(SZ) areau,xmax,xmin,ymax,ymin
      Real(SZ) dxdxi1,dydxi1,dxdxi2,dydxi2,dxi1dx,dxi2dx
      Real(SZ) dxi1dy,dxi2dy,ell_1,ell_2,ell_3
      Real(SZ) ZEVERTEX2(3),ZEVERTEX(3)

      Real(SZ),Allocatable :: tempmat(:,:),tempInv(:,:),tempag(:,:)
      Real(SZ),Allocatable :: AreaV_integral(:,:,:,:), A(:,:)
      Real(SZ),Allocatable :: Full_M_inv(:,:), temp_p(:,:),temp_t(:,:)
      Real(SZ),Allocatable :: Taylor_mass(:,:,:)

      Allocate ( tempmat(dg%dofh,dg%dofh),tempInv(dg%dofh,dg%dofh), tempag(dg%dofh,dg%dofh) )
      Allocate ( AreaV_integral(s%MNE,0:dg%ph,0:dg%ph,3), A(dg%dofh,dg%dofh) )
      Allocate ( Full_M_inv(dg%dofh,dg%dofh), temp_p(dg%dofh,dg%dofh) )
      Allocate ( Taylor_mass(s%MNE,dg%dofh,dg%dofh),temp_t(dg%dofh,dg%dofh))
      
!.....Initialize dg%SL3
      
      dg%SL3 = 0

!.....Loop over the elements

      DO J = 1,s%MNE
         
!.....Retrieve the nodal coordinates of the given element

         N1 = NM(J,1)
         N2 = NM(J,2)
         N3 = NM(J,3)

!.....Compute the barycenter coordinates of the element and store

         dg%XBC(J) = 1.D0/3.D0*(X(N1) + X(N2) + X(N3))
         dg%YBC(J) = 1.D0/3.D0*(Y(N1) + Y(N2) + Y(N3))

         X1 = dg%XBC(J)
         Y1 = dg%YBC(J)
         
!.....Loop over the edges to fnd the neighboring elements

         DO I = 1,3
            
!.....Retrieve the global edge number of the element for the given edge

            GED = NELED(I,J)

!.....Retrieve the neighboring element number and store into an array

            dg%EL_NBORS(I,J) = dg%NEDEL(1,GED)

            IF (dg%EL_NBORS(I,J).EQ.J) dg%EL_NBORS(I,J) = dg%NEDEL(2,GED)

!.....If the element has an edge that is on boundary go to next element
!     sb-2007/07/27 commented out

!     IF (dg%EL_NBORS(I,J).EQ.0) GOTO 111
            
         ENDDO
         
!.....Set the 4th neighbroing element equal to the first

         dg%EL_NBORS(4,J) = dg%EL_NBORS(1,J)

!.....Now loop over the three edges of the element again

         DO I = 1,3

            IF(dg%EL_NBORS(I,J).EQ.0.OR.dg%EL_NBORS(I+1,J).EQ.0) GOTO 110
            
!.....Compute the barycenter coordinates of two neighboring elements

            N1 = NM(dg%EL_NBORS(I,J),1)
            N2 = NM(dg%EL_NBORS(I,J),2)
            N3 = NM(dg%EL_NBORS(I,J),3)
            
            X2 = 1.D0/3.D0*(X(N1) + X(N2) + X(N3))
            Y2 = 1.D0/3.D0*(Y(N1) + Y(N2) + Y(N3))
            
            N1 = NM(dg%EL_NBORS(I+1,J),1)
            N2 = NM(dg%EL_NBORS(I+1,J),2)
            N3 = NM(dg%EL_NBORS(I+1,J),3)

            X3 = 1.D0/3.D0*(X(N1) + X(N2) + X(N3))
            Y3 = 1.D0/3.D0*(Y(N1) + Y(N2) + Y(N3))
            
!.....Compute the time independent planar constant

            dg%SL3(I,J) = X1*(Y2 - Y3) + X2*(Y3 - Y1) + X3*(Y1 - Y2)

            IF (dg%SL3(I,J).LE.0.AND.dg%SLOPEFLAG.NE.0) then
           WRITE(16,*) 'WARNING. dg%SL3(',I,',',J,') =',dg%SL3(I,J),' <= 0.',&
                   '    ELEMENT ',J,&
                   ' WILL NOT BE CONSIDERED IN SLOPE LIMITING.'
               WRITE(16,*)
            WRITE(*,*) 'WARNING. dg%SL3(',I,',',J,') =',dg%SL3(I,J),' <= 0.',&
                   '    ELEMENT ',J,&
                   ' WILL NOT BE CONSIDERED IN SLOPE LIMITING.'
               WRITE(*,*)
               dg%SL3(I,J) = 0.D0
            ENDIF
            
 110        CONTINUE

         ENDDO
         
      ENDDO


!******************************************************************************
!.....Vertex-based slope limiter (need the following stuff for integration)
!.....and must fill array for all possible p (ie. dg%dofl:dg%dofh)

      
#ifdef SLOPEALL

      dg%XBCb = 0.D0
      dg%YBCb = 0.D0
      dg%Deltx = 0.D0
      dg%Delty = 0.D0
      dg%xtransform = 0.D0
      dg%ytransform = 0.D0
      dg%xi1BCb = 0.D0
      dg%xi2BCb = 0.D0
      dg%xi1vert = 0.D0
      dg%xi2vert = 0.D0
      dg%xtransformv = 0.D0
      dg%ytransformv = 0.D0
      dg%YBCv = 0.D0
      dg%XBCv = 0.D0
      dg%xi1BCv = 0.D0
      dg%xi2BCv = 0.D0
      dg%f = 0.D0
      dg%g0 = 0.D0 
      dg%varsigma0 = 0.D0
      dg%fv = 0.D0
      dg%g0v = 0.D0 
      dg%varsigma0v = 0.D0
      dg%var2sigmag = 0.D0
      dg%NmatrixInv = 0.D0
      dg%Nmatrix = 0.D0


!.....Loop over p, and allocate for each order

      do ll = 1,dg%ph

         if (ll.gt.0) then
            
!.....Loop over the elements

            do k = 1,s%MNE

!.....Set areas over the physical elements
               
               areau = 0.5D0*AREAS(k)

!.....Retrieve the nodal coordinates (above) of the given element

               N1 = NM(k,1)
               N2 = NM(k,2)
               N3 = NM(k,3)

                                !areat = 0.5*(x(n2)*y(n3)-x(n3)*y(n2)+x(n3)*y(n1)-x(n1)*y(n3) + x(n1)*y(n2)-x(n2)*y(n1))
                                !areat = areau

!.....Find cell conditioners

               xmax = max( x(n1),x(n2),x(n3) )
               xmin = min( x(n1),x(n2),x(n3) )
               ymax = max( y(n1),y(n2),y(n3) )
               ymin = min( y(n1),y(n2),y(n3) )

               if (ll.le.2) then

                  dg%Deltx(k) = ( xmax - xmin ) / 2.D0 
                  dg%Delty(k) = ( ymax - ymin ) / 2.D0 

               else 

                  dg%Deltx(k) = ( xmax - xmin ) / ll 
                  dg%Delty(k) = ( ymax - ymin ) / ll

               endif

               

!.....Compute the centroid coordinates of the base element in physical space

               dg%XBCb(k) = ( x(n1) + x(n2) + x(n3) ) / 3.D0
               dg%YBCb(k) = ( y(n1) + y(n2) + y(n3) ) / 3.D0

!.....Transform quad points to physical space (for integration) xi --> x

               do mm=1,dg%nagp(ll)

                  ell_1 = -0.5D0 * ( dg%xagp(mm,dg%ph) + dg%yagp(mm,dg%ph) )
                  ell_2 =  0.5D0 * ( dg%xagp(mm,dg%ph) + 1.D0 )
                  ell_3 =  0.5D0 * ( dg%yagp(mm,dg%ph) + 1.D0 )

                  dg%xtransform(k,mm) = x(n1)*ell_1 + x(n2)*ell_2 + x(n3)*ell_3
                  dg%ytransform(k,mm) = y(n1)*ell_1 + y(n2)*ell_2 + y(n3)*ell_3

               enddo

!.....Find centroid coordinates in the master element frame

               dg%xi1BCb(k) =  ( (y(N3)-y(N1))*( dg%XBCb(k) -0.5D0 * &
              (x(N2) + x(N3))) + (x(N1) - x(N3))*(dg%YBCb(k)&
              - 0.5D0*(y(N2) + y(N3)) ) ) / areau
               dg%xi2BCb(k) =  ( (y(N1)-y(N2))*( dg%XBCb(k) -0.5D0 * (x(N2) &
              + x(N3))) + (x(N2) - x(N1))*(dg%YBCb(k)&
              - 0.5D0*(y(N2) + y(N3)) ) ) / areau

!.....Find vertices in the master element frame

               do lll=1,3

                  dg%xi1vert(k,lll) =  ( (y(N3)-y(N1))*( x(nm(k,lll)) -0.5D0 * &
                 (x(N2) + x(N3))) + (x(N1) - x(N3))* (y(nm(k,lll))&
                 - 0.5D0*(y(N2) + y(N3)) ) ) / areau
                  dg%xi2vert(k,lll) =  ( (y(N1)-y(N2))*( x(nm(k,lll)) -0.5D0 * &
                 (x(N2) + x(N3))) + (x(N2) - x(N1))* (y(nm(k,lll))&
                 - 0.5D0*(y(N2) + y(N3)) ) ) / areau

               enddo

!.....Find all neighbors of shared vertex

               do mm =1,S%MNE     !number of elements

                  do lll =1,3   !number of vertices

                     do nin =1,3 !number of vertices

                        if( NM(k,lll).eq.NM(mm,nin).and.k.ne.mm ) then !find common vertices of "nearby" elements
                           
!.....Compute the centroids of all conterminous (of codimension 2) elements (by vertex) of base element k in physical space

                           dg%XBCv(k,mm) = 1.D0/3.D0*( X(NM(mm,1)) + &
                          X(NM(mm,2)) + X(NM(mm,3)) )
                           dg%YBCv(k,mm) = 1.D0/3.D0*( Y(NM(mm,1)) + &
                          Y(NM(mm,2)) + Y(NM(mm,3)) )

                                !Stored_neighbors(k,mm,lll) = k*mm*lll !store neighboring elements

!.....Convert the centroid coordinates of all conterminous (of codimension 2) elements (by vertex) 
!.....of base element k to the master element space

                           dg%xi1BCv(k,mm) = ( (y(NM(mm,3))-y(NM(mm,1)))*&
                          ( dg%XBCv(k,mm) &
                          -0.5D0 * (x(Nm(mm,2)) + x(NM(mm,3)))) + &
                          (x(NM(mm,1)) - x(NM(mm,3)))*(dg%YBCv(k,mm) &
                          - 0.5D0*(y(NM(mm,2)) + y(NM(mm,3))) ) ) / areau
                           dg%xi2BCv(k,mm) = ( (y(NM(mm,1))-y(NM(mm,2)))*&
                          ( dg%XBCv(k,mm) &
                          -0.5D0 * (x(Nm(mm,2)) + x(NM(mm,3)))) + &
                          (x(NM(mm,2)) - x(NM(mm,1)))*(dg%YBCv(k,mm)&
                          - 0.5D0*(y(NM(mm,2)) + y(NM(mm,3))) ) ) / areau
                           
                        endif
                     enddo

                  enddo

               enddo


!.....Now compute the derivatives of the Taylor basis with respect to the physical 
!.....basis using the transformation rules from the paper (e.g. Leibniz and Faa' di Bruno formulas)
               
               dxdxi1 = 0.5D0 * ( x(N2) - x(N1) )
               dydxi1 = 0.5D0 * ( y(N2) - y(N1) )
               dxdxi2 = 0.5D0 * ( x(N3) - x(N1) )
               dydxi2 = 0.5D0 * ( y(N3) - y(N1) )

               dxi1dx = ( y(N3) - y(N1) ) / areau
               dxi2dx = ( y(N1) - y(N2) ) / areau
               dxi1dy = ( x(N1) - x(N3) ) / areau 
               dxi2dy = ( x(N2) - x(N1) ) / areau

!.....Write the generalized Taylor basis of order p in physical 
!.....coordinates (x(dg%xi1,dg%xi2), y(dg%xi1,dg%xi2)) and integrate over elements
!.....using the physical to master transformation, e.g. T^-1:x-->xi

               do i = 0,ll      !max polynomial degree in x

                  do j = 0,ll   !max polynomial degree in y

                     Call factorial(i,dg%fact(i))
                     Call factorial(j,dg%fact(j))

                     dg%Area_integral(k,i,j) = 0.D0


                     do mm = 1,dg%NAGP(ll) !number of quad points

                        dg%Area_integral(k,i,j) = dg%Area_integral(k,i,j) + &
                       ( (  dg%xtransform(k,mm) - dg%XBCb(k) )**i &
                       * (  dg%ytransform(k,mm)- dg%YBCb(k))**j * ( dg%wagp(mm,ll) ) )&
                       * abs( dxdxi1*dydxi2 - dxdxi2*dydxi1 )&
                       / ( dg%fact(i)*dg%fact(j)*dg%Deltx(k)**i * dg%Delty(k)**j )		


                     enddo
                     
                     do mm =1,dg%nagp(ll) !at quad points

                        dg%f(k,mm,i,j) = (  dg%xtransform(k,mm) - dg%XBCb(k)  )**i &
                       / ( dg%fact(i) * dg%Deltx(k)**i )
                        dg%g0(k,mm,i,j) = (  dg%ytransform(k,mm) - dg%YBCb(k)  )**j &
                       / ( dg%fact(j) * dg%Delty(k)**j )		

                        if (i.eq.0.and.j.eq.0) then

                           dg%varsigma0(k,mm,i,j) = 1

                        else
                           
                           dg%varsigma0(k,mm,i,j) = ( dg%f(k,mm,i,j) * dg%g0(k,mm,i,j) ) &
                          - dg%Area_integral(k,i,j)/areau

                        endif

                     enddo


                     do lll = 1,3 !at vertices

                        AreaV_integral(k,i,j,lll) = 0.D0

                        do mm = 1,dg%nagp(ll) !number of quad points

                           AreaV_integral(k,i,j,lll) = AreaV_integral(k,i,j,lll) + &
                          ( (  x(nm(k,lll)) - dg%XBCb(k) )**i &
                          * (  y(nm(k,lll)) - dg%YBCb(k) )**j * ( dg%wagp(mm,ll) ) )&
                          * abs( dxdxi1*dydxi2 - dxdxi2*dydxi1 )&
                          / ( dg%fact(i)*dg%fact(j)*dg%Deltx(k)**i * dg%Delty(k)**j )

                        enddo

                        dg%fv(k,lll,i,j) = (  x(nm(k,lll)) - dg%XBCb(k) )**i / &
                       ( dg%fact(i) * dg%Deltx(k)**i )
                        
                        dg%g0v(k,lll,i,j) =  ( y(nm(k,lll)) - dg%YBCb(k) )**j / &
                       ( dg%fact(j) * dg%Delty(k)**j )
                        
                        if (i.eq.0.and.j.eq.0) then

                           dg%varsigma0v(k,lll,i,j) = 1

                        else
                           
                           dg%varsigma0v(k,lll,i,j) = ( dg%fv(k,lll,i,j) * &
                          dg%g0v(k,lll,i,j) ) &
                          - dg%Area_integral(k,i,j)/areau

                        endif

                     enddo

                  enddo

               enddo



!.....Re-order the Taylor basis functions and componentwise derivatives into hierarchical order

               bbb = 1 
               do j = 0,ll

                  do i = 0,j

                     do mm = 1,dg%NAGP(ll)


                        dg%var2sigmag(k,mm,bbb) = dg%varsigma0(k,mm,i,j-i)

                        
                        if ( abs(dg%var2sigmag(k,mm,bbb)).lt.1.0E-15 ) then 

                           dg%var2sigmag(k,mm,bbb) = 0.D0 
                           
                        endif


                     enddo      !mm

                     do lll = 1,3

                        dg%var2sigmav(k,lll,bbb) = dg%varsigma0v(k,lll,i,j-i)

                        if ( abs(dg%var2sigmav(k,lll,bbb)).lt.1.0E-15 ) then 

                                !dg%var2sigmav(k,lll,bbb) = 0.D0 
                           
                        endif	

                     enddo      !lll

                     
                     dg%bi(bbb) = i
                     dg%bj(bbb) = j

                     bbb = bbb + 1

                  enddo         !i

               enddo            !j
               


!.....Compute the inner product matrix dg%Pmatrix, of the Taylor 
!.....basis with the Dubiner basis, and compute the transformation 
!.....matrix dg%Nmatrix=dg%Pmatrix*M(-1), using the mass matrix inverse dg%M_inv,

               ell = (ll+1)*(ll+2)/2

               A = 0.D0
               Full_M_inv = 0.D0
               temp_p = 0.D0



               do i = 1,ell

                  do j = 1,ell
                     
                     dg%pmatrix(k,i,j)=0.D0

                     do mm=1,dg%NAGP(ll)
                        
                        dg%pmatrix(k,i,j) = dg%pmatrix(k,i,j)  + dg%wagp(mm,ll) * &
                       dg%var2sigmag(k,mm,i) * dg%phi_area(j,mm,ll) 
                        

                        temp_p(i,j) = dg%pmatrix(k,i,j)

                        Taylor_mass(k,i,j) = Taylor_mass(k,i,j) +  dg%wagp(mm,ll) &
                       *  dg%var2sigmag(k,mm,i) *  dg%var2sigmag(k,mm,j)

                        temp_t(i,j) = Taylor_mass(k,i,j)

                     enddo
                     
                  enddo
                  
                  Full_M_inv(i,i) =   dg%M_INV(i,ll)
                  
               enddo

               Call Inv(temp_t(1:ell,1:ell), tempInv(1:ell,1:ell), ell)
               
               dg%Nmatrix(k,1:ell,1:ell,ell) = matmul(tempInv(1:ell,1:ell), temp_p(1:ell,1:ell))

               tempmat(1:ell,1:ell) = dg%Nmatrix(k,1:ell,1:ell,ell)

!.....Invert the transformation matrix dg%Nmatrix^(-1), 

               Call Inv(tempmat(1:ell,1:ell), tempInv(1:ell,1:ell), ell) 
               
               dg%NmatrixInv(k,1:ell,1:ell,ell) = tempInv(1:ell,1:ell)

            enddo               !k-elements
            
         endif                  !ll loop

      enddo                     !p_adapt

!.....Construct focal neighbors for non-vertex based limiters
      
      dg%focal_neigh = 0
      dg%focal_up = 0

      do j = 1,s%MNE

         bmm = 1

         do lll = 1,3

            do ell = 1,nneigh_elem(nm(j,lll))

               dg%focal_neigh(j,bmm) = neigh_elem(nm(j,lll), ell)  

               dg%focal_up(j) = bmm

               bmm = bmm + 1

            enddo

         enddo

      enddo

#endif

      RETURN
      END SUBROUTINE PREP_SLOPELIM

!.....Need factorial function for generalization

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

!.....Subroutine to find the inverse of a square matrix by Guass-Jordan elimination

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

