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

      SUBROUTINE PREP_SLOPELIM(s,dg_here)
      
!.....Use appropriate modules
      use sizes
      USE GLOBAL
      USE DG

      IMPLICIT NONE
      
      type (sizes_type) :: s
      type (dg_type) :: dg_here

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

      Allocate ( tempmat(dg_here%dofh,dg_here%dofh),tempInv(dg_here%dofh,dg_here%dofh), tempag(dg_here%dofh,dg_here%dofh) )
      Allocate ( AreaV_integral(s%MNE,0:dg_here%ph,0:dg_here%ph,3), A(dg_here%dofh,dg_here%dofh) )
      Allocate ( Full_M_inv(dg_here%dofh,dg_here%dofh), temp_p(dg_here%dofh,dg_here%dofh) )
      Allocate ( Taylor_mass(s%MNE,dg_here%dofh,dg_here%dofh),temp_t(dg_here%dofh,dg_here%dofh))
      
!.....Initialize dg_here%SL3
      
      dg_here%SL3 = 0

!.....Loop over the elements

      DO J = 1,s%MNE
         
!.....Retrieve the nodal coordinates of the given element

         N1 = NM(J,1)
         N2 = NM(J,2)
         N3 = NM(J,3)

!.....Compute the barycenter coordinates of the element and store

         dg_here%XBC(J) = 1.D0/3.D0*(X(N1) + X(N2) + X(N3))
         dg_here%YBC(J) = 1.D0/3.D0*(Y(N1) + Y(N2) + Y(N3))

         X1 = dg_here%XBC(J)
         Y1 = dg_here%YBC(J)
         
!.....Loop over the edges to fnd the neighboring elements

         DO I = 1,3
            
!.....Retrieve the global edge number of the element for the given edge

            GED = NELED(I,J)

!.....Retrieve the neighboring element number and store into an array

            dg_here%EL_NBORS(I,J) = dg_here%NEDEL(1,GED)

            IF (dg_here%EL_NBORS(I,J).EQ.J) dg_here%EL_NBORS(I,J) = dg_here%NEDEL(2,GED)

!.....If the element has an edge that is on boundary go to next element
!     sb-2007/07/27 commented out

!     IF (dg_here%EL_NBORS(I,J).EQ.0) GOTO 111
            
         ENDDO
         
!.....Set the 4th neighbroing element equal to the first

         dg_here%EL_NBORS(4,J) = dg_here%EL_NBORS(1,J)

!.....Now loop over the three edges of the element again

         DO I = 1,3

            IF(dg_here%EL_NBORS(I,J).EQ.0.OR.dg_here%EL_NBORS(I+1,J).EQ.0) GOTO 110
            
!.....Compute the barycenter coordinates of two neighboring elements

            N1 = NM(dg_here%EL_NBORS(I,J),1)
            N2 = NM(dg_here%EL_NBORS(I,J),2)
            N3 = NM(dg_here%EL_NBORS(I,J),3)
            
            X2 = 1.D0/3.D0*(X(N1) + X(N2) + X(N3))
            Y2 = 1.D0/3.D0*(Y(N1) + Y(N2) + Y(N3))
            
            N1 = NM(dg_here%EL_NBORS(I+1,J),1)
            N2 = NM(dg_here%EL_NBORS(I+1,J),2)
            N3 = NM(dg_here%EL_NBORS(I+1,J),3)

            X3 = 1.D0/3.D0*(X(N1) + X(N2) + X(N3))
            Y3 = 1.D0/3.D0*(Y(N1) + Y(N2) + Y(N3))
            
!.....Compute the time independent planar constant

            dg_here%SL3(I,J) = X1*(Y2 - Y3) + X2*(Y3 - Y1) + X3*(Y1 - Y2)

            IF (dg_here%SL3(I,J).LE.0.AND.dg_here%SLOPEFLAG.NE.0) then
           WRITE(16,*) 'WARNING. dg_here%SL3(',I,',',J,') =',dg_here%SL3(I,J),' <= 0.',&
                   '    ELEMENT ',J,&
                   ' WILL NOT BE CONSIDERED IN SLOPE LIMITING.'
               WRITE(16,*)
            WRITE(*,*) 'WARNING. dg_here%SL3(',I,',',J,') =',dg_here%SL3(I,J),' <= 0.',&
                   '    ELEMENT ',J,&
                   ' WILL NOT BE CONSIDERED IN SLOPE LIMITING.'
               WRITE(*,*)
               dg_here%SL3(I,J) = 0.D0
            ENDIF
            
 110        CONTINUE

         ENDDO
         
      ENDDO


!******************************************************************************
!.....Vertex-based slope limiter (need the following stuff for integration)
!.....and must fill array for all possible p (ie. dg_here%dofl:dg_here%dofh)

      
#ifdef SLOPEALL

      dg_here%XBCb = 0.D0
      dg_here%YBCb = 0.D0
      dg_here%Deltx = 0.D0
      dg_here%Delty = 0.D0
      dg_here%xtransform = 0.D0
      dg_here%ytransform = 0.D0
      dg_here%xi1BCb = 0.D0
      dg_here%xi2BCb = 0.D0
      dg_here%xi1vert = 0.D0
      dg_here%xi2vert = 0.D0
      dg_here%xtransformv = 0.D0
      dg_here%ytransformv = 0.D0
      dg_here%YBCv = 0.D0
      dg_here%XBCv = 0.D0
      dg_here%xi1BCv = 0.D0
      dg_here%xi2BCv = 0.D0
      dg_here%f = 0.D0
      dg_here%g0 = 0.D0 
      dg_here%varsigma0 = 0.D0
      dg_here%fv = 0.D0
      dg_here%g0v = 0.D0 
      dg_here%varsigma0v = 0.D0
      dg_here%var2sigmag = 0.D0
      dg_here%NmatrixInv = 0.D0
      dg_here%Nmatrix = 0.D0


!.....Loop over p, and allocate for each order

      do ll = 1,dg_here%ph

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

                  dg_here%Deltx(k) = ( xmax - xmin ) / 2.D0 
                  dg_here%Delty(k) = ( ymax - ymin ) / 2.D0 

               else 

                  dg_here%Deltx(k) = ( xmax - xmin ) / ll 
                  dg_here%Delty(k) = ( ymax - ymin ) / ll

               endif

               

!.....Compute the centroid coordinates of the base element in physical space

               dg_here%XBCb(k) = ( x(n1) + x(n2) + x(n3) ) / 3.D0
               dg_here%YBCb(k) = ( y(n1) + y(n2) + y(n3) ) / 3.D0

!.....Transform quad points to physical space (for integration) xi --> x

               do mm=1,dg_here%nagp(ll)

                  ell_1 = -0.5D0 * ( dg_here%xagp(mm,dg_here%ph) + dg_here%yagp(mm,dg_here%ph) )
                  ell_2 =  0.5D0 * ( dg_here%xagp(mm,dg_here%ph) + 1.D0 )
                  ell_3 =  0.5D0 * ( dg_here%yagp(mm,dg_here%ph) + 1.D0 )

                  dg_here%xtransform(k,mm) = x(n1)*ell_1 + x(n2)*ell_2 + x(n3)*ell_3
                  dg_here%ytransform(k,mm) = y(n1)*ell_1 + y(n2)*ell_2 + y(n3)*ell_3

               enddo

!.....Find centroid coordinates in the master element frame

               dg_here%xi1BCb(k) =  ( (y(N3)-y(N1))*( dg_here%XBCb(k) -0.5D0 * &
              (x(N2) + x(N3))) + (x(N1) - x(N3))*(dg_here%YBCb(k)&
              - 0.5D0*(y(N2) + y(N3)) ) ) / areau
               dg_here%xi2BCb(k) =  ( (y(N1)-y(N2))*( dg_here%XBCb(k) -0.5D0 * (x(N2) &
              + x(N3))) + (x(N2) - x(N1))*(dg_here%YBCb(k)&
              - 0.5D0*(y(N2) + y(N3)) ) ) / areau

!.....Find vertices in the master element frame

               do lll=1,3

                  dg_here%xi1vert(k,lll) =  ( (y(N3)-y(N1))*( x(nm(k,lll)) -0.5D0 * &
                 (x(N2) + x(N3))) + (x(N1) - x(N3))* (y(nm(k,lll))&
                 - 0.5D0*(y(N2) + y(N3)) ) ) / areau
                  dg_here%xi2vert(k,lll) =  ( (y(N1)-y(N2))*( x(nm(k,lll)) -0.5D0 * &
                 (x(N2) + x(N3))) + (x(N2) - x(N1))* (y(nm(k,lll))&
                 - 0.5D0*(y(N2) + y(N3)) ) ) / areau

               enddo

!.....Find all neighbors of shared vertex

               do mm =1,S%MNE     !number of elements

                  do lll =1,3   !number of vertices

                     do nin =1,3 !number of vertices

                        if( NM(k,lll).eq.NM(mm,nin).and.k.ne.mm ) then !find common vertices of "nearby" elements
                           
!.....Compute the centroids of all conterminous (of codimension 2) elements (by vertex) of base element k in physical space

                           dg_here%XBCv(k,mm) = 1.D0/3.D0*( X(NM(mm,1)) + &
                          X(NM(mm,2)) + X(NM(mm,3)) )
                           dg_here%YBCv(k,mm) = 1.D0/3.D0*( Y(NM(mm,1)) + &
                          Y(NM(mm,2)) + Y(NM(mm,3)) )

                                !Stored_neighbors(k,mm,lll) = k*mm*lll !store neighboring elements

!.....Convert the centroid coordinates of all conterminous (of codimension 2) elements (by vertex) 
!.....of base element k to the master element space

                           dg_here%xi1BCv(k,mm) = ( (y(NM(mm,3))-y(NM(mm,1)))*&
                          ( dg_here%XBCv(k,mm) &
                          -0.5D0 * (x(Nm(mm,2)) + x(NM(mm,3)))) + &
                          (x(NM(mm,1)) - x(NM(mm,3)))*(dg_here%YBCv(k,mm) &
                          - 0.5D0*(y(NM(mm,2)) + y(NM(mm,3))) ) ) / areau
                           dg_here%xi2BCv(k,mm) = ( (y(NM(mm,1))-y(NM(mm,2)))*&
                          ( dg_here%XBCv(k,mm) &
                          -0.5D0 * (x(Nm(mm,2)) + x(NM(mm,3)))) + &
                          (x(NM(mm,2)) - x(NM(mm,1)))*(dg_here%YBCv(k,mm)&
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
!.....coordinates (x(dg_here%xi1,dg_here%xi2), y(dg_here%xi1,dg_here%xi2)) and integrate over elements
!.....using the physical to master transformation, e.g. T^-1:x-->xi

               do i = 0,ll      !max polynomial degree in x

                  do j = 0,ll   !max polynomial degree in y

                     Call factorial(i,dg_here%fact(i))
                     Call factorial(j,dg_here%fact(j))

                     dg_here%Area_integral(k,i,j) = 0.D0


                     do mm = 1,dg_here%NAGP(ll) !number of quad points

                        dg_here%Area_integral(k,i,j) = dg_here%Area_integral(k,i,j) + &
                       ( (  dg_here%xtransform(k,mm) - dg_here%XBCb(k) )**i &
                       * (  dg_here%ytransform(k,mm)- dg_here%YBCb(k))**j * ( dg_here%wagp(mm,ll) ) )&
                       * abs( dxdxi1*dydxi2 - dxdxi2*dydxi1 )&
                       / ( dg_here%fact(i)*dg_here%fact(j)*dg_here%Deltx(k)**i * dg_here%Delty(k)**j )		


                     enddo
                     
                     do mm =1,dg_here%nagp(ll) !at quad points

                        dg_here%f(k,mm,i,j) = (  dg_here%xtransform(k,mm) - dg_here%XBCb(k)  )**i &
                       / ( dg_here%fact(i) * dg_here%Deltx(k)**i )
                        dg_here%g0(k,mm,i,j) = (  dg_here%ytransform(k,mm) - dg_here%YBCb(k)  )**j &
                       / ( dg_here%fact(j) * dg_here%Delty(k)**j )		

                        if (i.eq.0.and.j.eq.0) then

                           dg_here%varsigma0(k,mm,i,j) = 1

                        else
                           
                           dg_here%varsigma0(k,mm,i,j) = ( dg_here%f(k,mm,i,j) * dg_here%g0(k,mm,i,j) ) &
                          - dg_here%Area_integral(k,i,j)/areau

                        endif

                     enddo


                     do lll = 1,3 !at vertices

                        AreaV_integral(k,i,j,lll) = 0.D0

                        do mm = 1,dg_here%nagp(ll) !number of quad points

                           AreaV_integral(k,i,j,lll) = AreaV_integral(k,i,j,lll) + &
                          ( (  x(nm(k,lll)) - dg_here%XBCb(k) )**i &
                          * (  y(nm(k,lll)) - dg_here%YBCb(k) )**j * ( dg_here%wagp(mm,ll) ) )&
                          * abs( dxdxi1*dydxi2 - dxdxi2*dydxi1 )&
                          / ( dg_here%fact(i)*dg_here%fact(j)*dg_here%Deltx(k)**i * dg_here%Delty(k)**j )

                        enddo

                        dg_here%fv(k,lll,i,j) = (  x(nm(k,lll)) - dg_here%XBCb(k) )**i / &
                       ( dg_here%fact(i) * dg_here%Deltx(k)**i )
                        
                        dg_here%g0v(k,lll,i,j) =  ( y(nm(k,lll)) - dg_here%YBCb(k) )**j / &
                       ( dg_here%fact(j) * dg_here%Delty(k)**j )
                        
                        if (i.eq.0.and.j.eq.0) then

                           dg_here%varsigma0v(k,lll,i,j) = 1

                        else
                           
                           dg_here%varsigma0v(k,lll,i,j) = ( dg_here%fv(k,lll,i,j) * &
                          dg_here%g0v(k,lll,i,j) ) &
                          - dg_here%Area_integral(k,i,j)/areau

                        endif

                     enddo

                  enddo

               enddo



!.....Re-order the Taylor basis functions and componentwise derivatives into hierarchical order

               bbb = 1 
               do j = 0,ll

                  do i = 0,j

                     do mm = 1,dg_here%NAGP(ll)


                        dg_here%var2sigmag(k,mm,bbb) = dg_here%varsigma0(k,mm,i,j-i)

                        
                        if ( abs(dg_here%var2sigmag(k,mm,bbb)).lt.1.0E-15 ) then 

                           dg_here%var2sigmag(k,mm,bbb) = 0.D0 
                           
                        endif


                     enddo      !mm

                     do lll = 1,3

                        dg_here%var2sigmav(k,lll,bbb) = dg_here%varsigma0v(k,lll,i,j-i)

                        if ( abs(dg_here%var2sigmav(k,lll,bbb)).lt.1.0E-15 ) then 

                                !dg_here%var2sigmav(k,lll,bbb) = 0.D0 
                           
                        endif	

                     enddo      !lll

                     
                     dg_here%bi(bbb) = i
                     dg_here%bj(bbb) = j

                     bbb = bbb + 1

                  enddo         !i

               enddo            !j
               


!.....Compute the inner product matrix dg_here%Pmatrix, of the Taylor 
!.....basis with the Dubiner basis, and compute the transformation 
!.....matrix dg_here%Nmatrix=dg_here%Pmatrix*M(-1), using the mass matrix inverse dg_here%M_inv,

               ell = (ll+1)*(ll+2)/2

               A = 0.D0
               Full_M_inv = 0.D0
               temp_p = 0.D0



               do i = 1,ell

                  do j = 1,ell
                     
                     dg_here%pmatrix(k,i,j)=0.D0

                     do mm=1,dg_here%NAGP(ll)
                        
                        dg_here%pmatrix(k,i,j) = dg_here%pmatrix(k,i,j)  + dg_here%wagp(mm,ll) * &
                       dg_here%var2sigmag(k,mm,i) * dg_here%phi_area(j,mm,ll) 
                        

                        temp_p(i,j) = dg_here%pmatrix(k,i,j)

                        Taylor_mass(k,i,j) = Taylor_mass(k,i,j) +  dg_here%wagp(mm,ll) &
                       *  dg_here%var2sigmag(k,mm,i) *  dg_here%var2sigmag(k,mm,j)

                        temp_t(i,j) = Taylor_mass(k,i,j)

                     enddo
                     
                  enddo
                  
                  Full_M_inv(i,i) =   dg_here%M_INV(i,ll)
                  
               enddo

               Call Inv(temp_t(1:ell,1:ell), tempInv(1:ell,1:ell), ell)
               
               dg_here%Nmatrix(k,1:ell,1:ell,ell) = matmul(tempInv(1:ell,1:ell), temp_p(1:ell,1:ell))

               tempmat(1:ell,1:ell) = dg_here%Nmatrix(k,1:ell,1:ell,ell)

!.....Invert the transformation matrix dg_here%Nmatrix^(-1), 

               Call Inv(tempmat(1:ell,1:ell), tempInv(1:ell,1:ell), ell) 
               
               dg_here%NmatrixInv(k,1:ell,1:ell,ell) = tempInv(1:ell,1:ell)

            enddo               !k-elements
            
         endif                  !ll loop

      enddo                     !p_adapt

!.....Construct focal neighbors for non-vertex based limiters
      
      dg_here%focal_neigh = 0
      dg_here%focal_up = 0

      do j = 1,s%MNE

         bmm = 1

         do lll = 1,3

            do ell = 1,nneigh_elem(nm(j,lll))

               dg_here%focal_neigh(j,bmm) = neigh_elem(nm(j,lll), ell)  

               dg_here%focal_up(j) = bmm

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

