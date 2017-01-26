      SUBROUTINE initial()

      USE globals, ONLY: rp,ne,nn,mndof,ndof,nel_type,mnnds,nnds,np,nverts, &
                         ect,xy,depth,elxy,elhb,hbnodes, &
                         nqpta,qpta,wpta,nqpte,qpte, &
                         H,Qx,Qy, &
                         Zinit,Hinit,Qxinit,Qyinit, &
                         phia,phil, &
                         el_type,psia, &
                         detJa,mmi_init
                         
      USE allocation, ONLY: alloc_sol_arrays   
      USE basis, ONLY: tri_nodes,tri_basis,quad_nodes,quad_basis
      USE read_dginp, ONLY: out_direc,grid_file
      USE messenger2, ONLY: myrank      

      IMPLICIT NONE
      INTEGER :: i,j,el,l,pt,dof,ed,nd,k,m,n,p
      INTEGER :: et,nv
      INTEGER :: mnqpte
      INTEGER :: ipiv(mnnds,nel_type),info
      REAL(rp) :: x1,x2,x3,y1,y2,y3,x,y
      REAL(rp) :: sigma,xc,yc,h0,h00
      REAL(rp) :: qint,f,qint2
      REAL(rp) :: xn,yn
      REAL(rp),ALLOCATABLE :: rhsH(:,:),rhsH2(:,:)      
      REAL(rp) :: r(mnnds),s(mnnds),phi(mnnds*mnnds),hb(mnnds)
      REAL(rp) :: V(mnnds,mnnds,nel_type)

      sigma = 10d0
      xc = 0d0 
      yc = 0d0
      
      IF (myrank == 0) PRINT "(A)", "Computing initial condition..."
      CALL alloc_sol_arrays()

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Compute initial condition
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ALLOCATE(rhsH(ne,mndof))
      ALLOCATE(rhsH2(ne,mndof))
      
      Hinit = 0d0
      
!       DO el = 1,ne
!         et = el_type(el)
!         nv = nverts(et)
!         DO l = 1,ndof(et)
!           qint = 0d0
!           qint2 = 0d0
!           DO pt = 1,nqpta(et)
!           
!             r = qpta(pt,1,et)
!             s = qpta(pt,2,et)
!             
!             x = 0d0
!             y = 0d0
!             DO nd = 1,nv
! 
!               xn = elxy(nd,el,1)
!               yn = elxy(nd,el,2)            
!               
!               x = x + psia(nd,pt,et)*xn
!               y = y + psia(nd,pt,et)*yn
!             ENDDO
!             
! !             x1 = xy(1,ect(1,el))
! !             x2 = xy(1,ect(2,el))
! !             x3 = xy(1,ect(3,el))
! ! 
! !             y1 = xy(2,ect(1,el))
! !             y2 = xy(2,ect(2,el))
! !             y3 = xy(2,ect(3,el))
! 
! !             x = -.5d0*((r+s)*x1 - (1d0+r)*x2 - (1d0+s)*x3)
! !             y = -.5d0*((r+s)*y1 - (1d0+r)*y2 - (1d0+s)*y3)
! 
!             h0 = 0d0
!             DO nd = 1,nv   ! This assumes there is an equal order representation between the bathymetry and the coordinate transformation
!               h0 = h0 + psia(nd,pt,et)*elhb(nd,el)
!             ENDDO
! 
! !             h00 = depth(ect(1,el))*phil(1,pt,1) + depth(ect(2,el))*phil(2,pt,1) + depth(ect(3,el))*phil(3,pt,1)
! 
! 
! !             f = exp(-sigma*((x-xc)**2d0+(y-yc)**2d0))+h0 
! !               f = .0002*x + h0
!             f = h0
!             
!             qint = qint + wpta(pt,et)*f*phia(l,pt,et)*detJa(el,pt)
!             qint2 = qint2 + 0.5d0*wpta(pt,et)*f*phia(l,pt,et)
!           ENDDO
!           rhsH(el,l) = qint
!           rhsH2(el,l) = qint2
! !           Hinit(el,1) = el ! for debugging
!         ENDDO
!         
!         m = 1
!         DO i = 1,ndof(et)        
!           DO j = 1,ndof(et)
!             Hinit(el,i) = Hinit(el,i) + mmi_init(el,m)*rhsH(el,j)
!             m = m + 1
!           ENDDO
!         ENDDO
!       ENDDO
      




      Qxinit(:,:) = 0d0
      Qyinit(:,:) = 0d0
      Zinit(:,:) = 0d0      



      RETURN
      END SUBROUTINE
