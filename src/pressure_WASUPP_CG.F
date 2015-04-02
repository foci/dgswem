!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    SUBROUTINE pressure_WASUPP_CG 
!    Aaron Donahue and Jessica Meixner 
!    March 2015
!    Calculates dynamic pressure 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      
      
      SUBROUTINE pressure_WASUPP_CG
      USE AddSPARSKIT
      USE sizes,    ONLY: SZ 
      USE global,   ONLY: G,DT, ne, areas, mnp, x, y, nm, eta2, uu2, vv2, dp
      USE pressure_module, ONLY: NumCols, NumRows, 
     &                            Eta_pp,U_pp,V_pp,Depth_pp, 
     &                            Etax_pp,Etay_pp,Ux_pp,Uy_pp,Vx_pp,Vy_pp,
     &                            Etaxx_pp,Etayy_pp,Depx_pp,Depy_pp,Depxx_pp,Depyy_pp,Depxy_pp,
     &                            p1_pp,p1x_pp,p1y_pp,p1old,p1oldest, 
     &                            totaldepth_pp, wetdry_pp,
     &                            DX,DY,pd_pp,pb_pp, !!! ASD Need to add these
     &                            totalnz, pp_dg_cg_mat, pp_dg_cg_ilumat,
     &                            phinodal,csrloc,
     &                            pp_p1,pp_p2,pp_pd,pp_pb,
     &                            IWASUPP
      USE DG,      ONLY: wagp,nagp
      
      IMPLICIT NONE
      
      INTEGER                       :: VTYPE,NX,NY,NT,Nghost 
      LOGICAL                       :: PERIODIC=.FALSE. 
      REAL(SZ),ALLOCATABLE           :: rhs(:),p1sol(:)
      REAL(SZ)                       :: DropTol,RESID,infor(4)
      INTEGER                       :: ifil,initmem,itmax,infoi(4)
      !c FD matrix
      TYPE (spmat):: fdmat
      TYPE (spmat), pointer:: ilumat
      
      INTEGER                 :: l,pad,i,j,pt
      REAL(SZ),DIMENSION(3,3) :: p1xcf,p1ycf,p1cf
      REAL(SZ),DIMENSION(3)   :: prhs,pdrhs,pbrhs,phix,phiy,xloc,yloc
      REAL(SZ)                :: cgp1(mnp),cgpd(mnp),cgpb(mnp)!,pp_p2(mnp),pp_p1(mnp)
      REAL(SZ)	              :: zi,zx,zy
      REAL(SZ)	              :: hi,hx,hy,hxy
      REAL(SZ)	              :: ui,ux,uy
      REAL(SZ)	              :: vi,vx,vy
      REAL(SZ)	              :: pi,px,py
      REAL(SZ)	              :: arealoc
      TYPE (spmat):: pp_p1_cg_mat
      TYPE (spmat), pointer:: pp_p1_cg_ilumat
      
      
      allocate(pp_p1_cg_mat.idv1(totalnz+2),pp_p1_cg_mat.idv2(mnp+1),pp_p1_cg_mat.val(totalnz+2))
      pp_p1_cg_mat.str = 'csr'
      pp_p1_cg_mat.ndim = mnp
      pp_p1_cg_mat.val(:) = 0.d0
      cgp1(:) = 0.d0
      
      pp_p1(:) = 0.d0
      pp_p2(:) = 0.d0
      pp_pd(:) = 0.d0
      pp_pb(:) = 0.d0
      
      if (IWASUPP.eq.0) then
        goto 100
        
      else
      
      pad = 1
!       open( unit = 60, file = 'p1coeff.crs' ) ; 
!       open( unit = 70, file = 'phixphiy.crs' ) ;
      do l = 1,ne
        arealoc = areas(l)/2.d0
        do i = 1,3
          xloc(i) = x(nm(l,i))
          yloc(i) = y(nm(l,i))
        end do
        do i = 1,3
          phix(i) = 1.d0/2.d0/arealoc*(yloc(mod(i,3)+1)-yloc(mod(i+1,3)+1))
          phiy(i) = 1.d0/2.d0/arealoc*(xloc(mod(i+1,3)+1)-xloc(mod(i,3)+1))
!           write(70,'(2i,2f32.16)') l, i, phix(i), phiy(i)
        end do
        
        ! integrate using guass quadrature
        p1xcf(:,:) = 0.d0
        p1ycf(:,:) = 0.d0
        p1cf(:,:)  = 0.d0
        prhs(:)    = 0.d0
        do pt = 1,nagp(pad)
          zi = 0.d0
          zx = 0.d0
          zy = 0.d0
          
          hi = 0.d0
          hx = 0.d0
          hy = 0.d0
          
          ui = 0.d0
          ux = 0.d0
          uy = 0.d0
          
          vi = 0.d0
          vx = 0.d0
          vy = 0.d0
          
          hxy = 0.d0
          
          do j = 1,3
            zi = zi + eta2(nm(l,j))*phinodal(j,pt)
            zx = zx + eta2(nm(l,j))*phix(j)
            zy = zy + eta2(nm(l,j))*phiy(j)
            
            hi = hi + dp(nm(l,j))*phinodal(j,pt)
            hx = hx + dp(nm(l,j))*phix(j)
            hy = hy + dp(nm(l,j))*phiy(j)
            
            ui = ui + uu2(nm(l,j))*phinodal(j,pt)
            ux = ux + uu2(nm(l,j))*phix(j)
            uy = uy + uu2(nm(l,j))*phiy(j)
            
            vi = vi + vv2(nm(l,j))*phinodal(j,pt)
            vx = vx + vv2(nm(l,j))*phix(j)
            vy = vy + vv2(nm(l,j))*phiy(j)
            
!             print*, l,pt,j,ui
!             if (isnan(ui)) then
!               print*, phinodal(j,pt), uu2(nm(l,j))
!               call exit
!             end if
            
          end do          
          
          do i = 1,3
            do j = 1,3
              p1xcf(i,j) = p1xcf(i,j) + arealoc/2.d0*wagp(pt,pad)*(
     &                        (1.d0/2.d0*hi**2*phix(i))                       ! Linear Flat
     &                       +(1.d0/4.d0*hi*hx*phinodal(i,pt))                ! Linear Vary
     &                       +(1.d0/6.d0*(3*zi**2*phix(i)+4*hi*zx*phinodal(i,pt) ! Nonlin Flat
     &                          +zi*(6*hi*phix(i)+4*zx*phinodal(i,pt))))
     &                       +(1.d0/4.d0*zi*hx*phinodal(i,pt))                ! Nonlin Vary
     &                                                           )*phix(j)
              p1ycf(i,j) = p1ycf(i,j) + arealoc/2.d0*wagp(pt,pad)*(
     &                        (1.d0/2.d0*hi**2*phiy(i))                       ! Linear Flat
     &                       +(1.d0/4.d0*hi*hy*phinodal(i,pt))                ! Linear Vary
     &                       +(1.d0/6.d0*(3*zi**2*phiy(i)+4*hi*zy*phinodal(i,pt) ! Nonlin Flat
     &                          +zi*(6*hi*phiy(i)+4*zy*phinodal(i,pt))))      
     &                       +(1.d0/4.d0*zi*hy*phinodal(i,pt))                ! Nonlin Vary
     &                                                           )*phiy(j)
              p1cf(i,j)  = p1cf(i,j)  + arealoc/2.d0*wagp(pt,pad)*(
     &                        (5.d0/4.d0*phinodal(i,pt))                      ! Linear Flat
     &                       +(-1.d0/4.d0*((hx**2+hy**2)*phinodal(i,pt)+2*hi  ! Linear Vary
     &                          *(hy*phiy(i)+hx*phix(i))))
     &                       +(1.d0/3.d0*((hi+zi)*(zx*phix(i)+zy*phiy(i))     ! Nonlin Flat
     &                          +3.d0*phinodal(i,pt)*(zx**2+zy**2)))
     &                       +(-1.d0/2.d0*(zi*(hx*phix(i)+hy*phiy(i))         ! Nonlin Vary
     &                          +phinodal(i,pt)*(hy*zy+hx*zx)))
     &                                                           )*phinodal(j,pt)     
            end do
            
            
            prhs(i) = prhs(i) - arealoc/2.d0*wagp(pt,pad)*(
     &                   (1.d0/12.d0*g*(-15.d0*zi*phinodal(i,pt)+4.d0*hi**2   ! Linear Flat
     &                      *(zx*phix(i)+zy*phiy(i))))
     &                  +(1.d0/12.d0*g*(3.d0*zi*(phinodal(i,pt)*(hy**2+hx**2) ! Linear Vary
     &                      +2.d0*hi*(hy*phiy(i)+hx*phix(i)))
     &                      +2.d0*hi*phinodal(i,pt)*(hy*zy+hx*zx)))
     &                  +(1.d0/6.d0*(10.d0*zi**2*phinodal(i,pt)               ! Nonlin Flat
     &                      *(vy**2+ux*vy*ux**2+uy*vx)+hi*phinodal(i,pt)
     &                      *(15.d0*ui*vi*hxy+2.d0*(3.d0*g*(zy**2+zx**2)
     &                      +5.d0*hi*(vy**2+vy*ux+ux**2+uy*vx)))
     &                      +zi*(15.d0*ui*vi*hxy*phinodal(i,pt)
     &                      +2.d0*hi*g*(zx*phix(i)+zy*phiy(i))
     &                      +10.d0*phinodal(i,pt)*(vy**2+vy*ux+ux**2+uy*vx)))) 
     &                  +(-1.d0/12.d0*(6.d0*g*zi**2*(hy*phiy(i)+hx*phix(i))  ! Nonlin Vary
     &                      +15.d0*vi**2*hy*((zi+hi)*phiy(i)+(zy+hy)*phinodal(i,pt))
     &                      +30.d0*(zi+hi)*vi*vy*hy+15*ui*hx*(
     &                       ui*(hi*phix(i)+phinodal(i,pt)*(hx+zx))+2.d0*hi*ux*phinodal(i,pt))
     &                      +zi*(15.d0*ui**2*hx*phix(i)-8.d0*g*phinodal(i,pt)*(hy*zy+hx*zx)
     &                      +30.d0*ui*phinodal(i,pt)*hx*ux)))
     &                                                     )
     
!           write(*,'(7f16.8)') prhs(i), ui, vi, ux ,vx, uy, vy
          end do          
          
          
          
          ! Now construct LHS and RHS CG matrices through Gauss integration
          do j = 1,3
            cgp1(nm(l,j)) = cgp1(nm(l,j)) + prhs(j)
            do i = 1,3
!               write(60,'(3i,4f32.16)') l, i, j, p1xcf(i,j), p1ycf(i,j), p1cf(i,j), prhs(j)
              pp_p1_cg_mat.val(csrloc(l,j,i)) = pp_p1_cg_mat.val(csrloc(l,j,i))+
     &                                               (p1xcf(j,i)+p1ycf(j,i)+p1cf(j,i))
            end do
          end do                   
        end do        
        
      end do
!       close(60)
!       close(70)
      pp_p1_cg_mat.idv1(:) = pp_dg_cg_mat.idv1(:)
      pp_p1_cg_mat.idv2(:) = pp_dg_cg_mat.idv2(:)
      
!       open( unit = 10, file = 'valp1.crs' ) ; 
!       open( unit = 20, file = 'colp1.crs' ) ; 
!       open( unit = 30, file = 'ptrp1.crs' ) ; 
!       open( unit = 40, file = 'rhsp1.crs' ) ; 
!       open( unit = 50, file = 'p1areas.crs' ) ; 
!       DO I = 1, ubound(pp_p1_cg_mat.val,1) 
! 	    WRITE(10,'(2f32.16)') pp_p1_cg_mat.val(i), pp_dg_cg_mat.val(i) ;
! 	    WRITE(20,*) pp_p1_cg_mat.idv1(i), pp_dg_cg_mat.idv1(i) ;
! 	    WRITE(10,*) pp_p1_cg_mat.val(i) ;
! 	    WRITE(20,*) pp_p1_cg_mat.idv1(i) ;
!       END DO
!       close(10) ;
!       close(20) ;
      
!       DO I = 1, mnp + 1
! 	    WRITE(30,*) pp_p1_cg_mat.idv2(i), pp_dg_cg_mat.idv2(i) ; 
! 	    if (i.le.mnp) then
! 	      WRITE(40,'(f32.16)') cgp1(i)
!               print*, i, cgp1(i)
! 	    end if
!       END DO
!       do l = 1,ne
!         write(50,'(4F32.16)') areas(l)/2.d0, dp(nm(l,1)), dp(nm(l,2)), dp(nm(l,3))
!       end do
!       close(30)
!       close(40)
!       close(50)

      
      
      ! Now find p1 solution
      DropTol = 0.00001D0 ; 
      ifil = 21 ; !c ifil -- 3*bandwidth
      initmem = 20 ; 
      pp_p1_cg_ilumat => WrpILUT( pp_p1_cg_mat, DROPTOL, IFILOP = IFIL, INITMEMOP = initmem ) ;
      
      pp_p1 = WrpSPBCCGSTAB0( cgp1, pp_p1, pp_p1_cg_mat, pp_p1_cg_ilumat, 1.d-8, mnp/4, INFOR, INFOI ) ;
      
      ! Find P2, pd, and pb for coupling with SWE
      cgp1(:) = 0.d0
      cgpd(:) = 0.d0
      cgpb(:) = 0.d0
      
      do l = 1,ne
        arealoc = areas(l)/2.d0
        do i = 1,3
          xloc(i) = x(nm(l,i))
          yloc(i) = y(nm(l,i))
        end do
        do i = 1,3
          phix(i) = 1.d0/2.d0/arealoc*(yloc(mod(i,3)+1)-yloc(mod(i+1,3)+1))
          phiy(i) = 1.d0/2.d0/arealoc*(xloc(mod(i+1,3)+1)-xloc(mod(i,3)+1))
        end do
        
        prhs(:)    = 0.d0
        pdrhs(:)   = 0.d0
        pbrhs(:)   = 0.d0
        
        do pt = 1,nagp(pad)
          pi = 0.d0
          px = 0.d0
          py = 0.d0
          
          zi = 0.d0
          
          hi = 0.d0
          hx = 0.d0
          hy = 0.d0
          
          ui = 0.d0
          
          vi = 0.d0

          
          do j = 1,3
            zi = zi + eta2(nm(l,j))*phinodal(j,pt)
            
            hi = hi + dp(nm(l,j))*phinodal(j,pt)
            hx = hx + dp(nm(l,j))*phix(j)
            hy = hy + dp(nm(l,j))*phiy(j)
            
            pi = pi + pp_p1(nm(l,j))*phinodal(j,pt)
            px = px + pp_p1(nm(l,j))*phix(j)
            py = py + pp_p1(nm(l,j))*phiy(j)
            
            ui = ui + uu2(nm(l,j))*phinodal(j,pt)
            
            vi = vi + vv2(nm(l,j))*phinodal(j,pt)
          end do
          
          do j = 1,3
            prhs(j) = prhs(j) - arealoc/2.d0*wagp(pt,pad)*3.d0/4.d0*(
     &                  (1.d0+hx**2+hy**2)*(g*zi-pi)
     &                 +(hi+zi)*(hx*px+hy*py)    ! ASD
     &                                                              )*phinodal(j,pt)/4.d0  ! ASD I have no idea why /4!!!!!
            pdrhs(j) = pdrhs(j) + arealoc/2.d0*wagp(pt,pad)/2.d0*(zi+hi)*(
     &                               pi-g*zi )*phinodal(j,pt)/4.d0
            pbrhs(j) = pbrhs(j) - arealoc/2.d0*wagp(pt,pad)/3.d0*(
     &                               3.d0*g*zi-3.d0*pi)*phinodal(j,pt)/4.d0

          end do

          do j = 1,3
            cgp1(nm(l,j)) = cgp1(nm(l,j)) + prhs(j)
            cgpd(nm(l,j)) = cgpd(nm(l,j)) + pdrhs(j)
            cgpb(nm(l,j)) = cgpb(nm(l,j)) + pbrhs(j)
          end do
        
        end do        
      end do
      
      pp_p2 = WrpSPBCCGSTAB0( cgp1, pp_p2, pp_dg_cg_mat, pp_dg_cg_ilumat, 1.d-8, mnp/4, INFOR, INFOI )
      
      ! Just add in the contribution from p2
      do l = 1,ne
        arealoc = areas(l)/2.d0
        do i = 1,3
          xloc(i) = x(nm(l,i))
          yloc(i) = y(nm(l,i))
        end do
        do i = 1,3
          phix(i) = 1.d0/2.d0/arealoc*(yloc(mod(i,3)+1)-yloc(mod(i+1,3)+1))
          phiy(i) = 1.d0/2.d0/arealoc*(xloc(mod(i+1,3)+1)-xloc(mod(i,3)+1))
        end do
        
        pbrhs(:)   = 0.d0
        
        do pt = 1,nagp(pad)
          pi = 0.d0
          do j = 1,3
            pi = pi + pp_p2(nm(l,j))*phinodal(j,pt)
          end do
          do j = 1,3
            pbrhs(j) = pbrhs(j) - arealoc/2.d0*wagp(pt,pad)/3.d0*(pi)*phinodal(j,pt)/4.d0 !!! HERE TOO
          end do
        
          do j = 1,3
            cgpb(nm(l,j)) = cgpb(nm(l,j)) + pbrhs(j)
          end do
        
        end do        
      end do
      
      pp_pd = WrpSPBCCGSTAB0( cgpd, pp_pd, pp_dg_cg_mat, pp_dg_cg_ilumat, 1.d-8, mnp/4, INFOR, INFOI )
      pp_pb = WrpSPBCCGSTAB0( cgpb, pp_pb, pp_dg_cg_mat, pp_dg_cg_ilumat, 1.d-8, mnp/4, INFOR, INFOI )

      end if
      
100   continue
      
      CALL EmptySpmat(pp_p1_cg_mat)
      CALL EmptySpmat(pp_p1_cg_ilumat)
      deallocate(pp_p1_cg_ilumat)
      
      
      
      
!       call exit
! 
!       !.Allocate local variables
!             
!        NX = NumRows 
!        NY = NumCols 
!        Nghost = 0
!        NT = NX*NY
!        ALLOCATE(rhs(NT),p1sol(NT))
!       
!         fdmat.str = 'csr'  ;
!         fdmat.ndim = NT ; !c number of rows 
!       
!         allocate(fdmat.val(5*NT)) ;  !c estimate nnz for second order matrix
!         allocate(fdmat.idv1(5*NT)) ; !c column index (col_ind)
!         allocate(fdmat.idv2(NT+1)) ; !c row index (row_prt)
!       !.First Derivatives  
!         CALL DERIVATIVE_X(NX,NY,1,NX,1,NY,wetdry_pp,DX,U_pp,Ux_pp)
!         CALL DERIVATIVE_X(NX,NY,1,NX,1,NY,wetdry_pp,DX,V_pp,Vx_pp)
!         CALL DERIVATIVE_X(NX,NY,1,NX,1,NY,wetdry_pp,DX,Eta_pp,Etax_pp)
!         CALL DERIVATIVE_X(NX,NY,1,NX,1,NY,wetdry_pp,DX,Depth_pp,Depx_pp)
!       
!         
!         CALL DERIVATIVE_Y(NX,NY,1,NX,1,NY,wetdry_pp,DY,U_pp,Uy_pp)
!         CALL DERIVATIVE_Y(NX,NY,1,NX,1,NY,wetdry_pp,DY,V_pp,Vy_pp)
!         CALL DERIVATIVE_Y(NX,NY,1,NX,1,NY,wetdry_pp,DY,Eta_pp,Etay_pp)
!         CALL DERIVATIVE_Y(NX,NY,1,NX,1,NY,wetdry_pp,DY,Depth_pp,Depy_pp)
!       
!       !.Second Derivatives
!         CALL DERIVATIVE_XX(NX,NY,1,NX,1,NY,wetdry_pp,DX,Eta_pp,Etaxx_pp)
!         CALL DERIVATIVE_YY(NX,NY,1,NX,1,NY,wetdry_pp,DY,Eta_pp,Etayy_pp)
!         CALL DERIVATIVE_XX(NX,NY,1,NX,1,NY,wetdry_pp,DX,Depth_pp,Depxx_pp)
!         CALL DERIVATIVE_YY(NX,NY,1,NX,1,NY,wetdry_pp,DY,Depth_pp,Depyy_pp)
!         CALL DERIVATIVE_XY(NX,NY,1,NX,1,NY,wetdry_pp,DX,DY,Depth_pp,Depxy_pp)
!       
!         
!       !.Determine RHS vector of Pressure Poisson problem for P1
!         CALL GETRHS()
!       
!       !.Construct the left hand side Matrix for solving P1
!         CALL FDDISCRETE()   
!         
! !         CAll testilut(totalnz,fdmat)
! !         call exit
!       
!       !.Precondition Matrix  
!         DropTol = 0.0001D0 ; 
!         ifil = 20 ; !c ifil -- 3*bandwidth
!         initmem = 20 ; 
!         ilumat => WrpILUT( fdmat, DROPTOL, IFILOP = IFIL, INITMEMOP = initmem ) ;
!       
!       !.Solve Matrix using BiCGSTAB 
!       !   p1sol = 0.0D0
!         RESID = 1.d-4 ; 
!         ITMAX = NT/4 ;
!         p1sol = WrpSPBCCGSTAB0( rhs, p1sol, fdmat, ilumat, RESID, ITMAX, INFOR, INFOI ) ;
!         p1oldest = p1old
!         p1old = p1sol
!         
!       !   print*, 'Iter = ',infoi(2)
!       
!       !.Find Pd and Pb from p1
!         p1_pp = 0.0D0
!         p1_pp(1:NX,1:NY) = reshape(p1sol,(/NX,NY/))
!         VTYPE = 1
!       !   CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,p1,VTYPE,PERIODIC)
!         CALL DERIVATIVE_X(NX,NY,1,NX,1,NY,wetdry_pp,DX,p1_pp,p1x_pp)
!         CALL DERIVATIVE_Y(NX,NY,1,NX,1,NY,wetdry_pp,DY,p1_pp,p1y_pp)
!         CALL GETPDPB()
!         
!         ! Assign Ghost Node Values
!         VTYPE = 1
!       !   CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Pd,VTYPE,PERIODIC)
!       !   CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Pb,VTYPE,PERIODIC)
!       
!         ! Deallocate arrays
!         CALL EmptySpmat( fdmat );
!         CALL EmptySpmat( ilumat ) ;
!         deallocate( ilumat ) ;
!         deallocate(rhs,p1sol) ;
!       CONTAINS
!       !==================================================================================!
!         SUBROUTINE GETRHS( )
!           IMPLICIT NONE
!           INTEGER  :: I,J,id
!           
!           REAL(SZ) :: rhsx, rhsy
!       
!       !     open(unit=99,file='p1sol.txt',status='replace')
!           DO J = 1,NY
!             DO I = 1,NX
!               
!               id = (J-Nghost-1)*(NX)+(I-Nghost)
!       
!               IF (wetdry_pp(I,J).EQ.0) THEN
!                 rhs(id) = G*Eta_pp(I,J)
!               ELSE
!       !..................Linear Flatbed
!                rhs(id) = 1.0D0/12.0D0*G*( 15.0D0*Eta_pp(I,J)+4.0D0*Depth_pp(I,J)**2* 
!      &                    (Etaxx_pp(I,J)+Etayy_pp(I,J)) )
!                rhsx = 1.0D0/12.0D0*G*( 15.0D0*Eta_pp(I,J)/2.d0+4.0D0*Depth_pp(I,J)**2* 
!      &                    (Etaxx_pp(I,J)) )
!                rhsy = 1.0D0/12.0D0*G*( 15.0D0*Eta_pp(I,J)/2.d0+4.0D0*Depth_pp(I,J)**2* 
!      &                    (Etayy_pp(I,J)) )
!       !..................NonLinear Flatbed
!                rhs(id) = rhs(id) + 1.0D0/3.0D0*G*Depth_pp(I,J)*(-2.0D0*(Etax_pp(I,J)**2+Etay_pp(I,J)**2) 
!      &                   +Eta_pp(I,J)*(Etaxx_pp(I,J)+Etayy_pp(I,J)))
!                rhsx    = rhsx + 1.0D0/3.0D0*G*Depth_pp(I,J)*(-2.0D0*(Etax_pp(I,J)**2) 
!      &                   +Eta_pp(I,J)*(Etaxx_pp(I,J)))
!                rhsy    = rhsy + 1.0D0/3.0D0*G*Depth_pp(I,J)*(-2.0D0*(Etay_pp(I,J)**2) 
!      &                   +Eta_pp(I,J)*(Etayy_pp(I,J)))
!       !..................Linear Varying Bathy
!                rhs(id) = rhs(id) + 1.0D0/4.0D0*G*( 4.0D0*Depth_pp(I,J)*(Depy_pp(I,J)*Etay_pp(I,J) 
!      &                   +Depx_pp(I,J)*Etax_pp(I,J))+Eta_pp(I,J)*(Depy_pp(I,J)**2+Depx_pp(I,J)**2 
!      &                   +2.0D0*Depth_pp(I,J)*(Depxx_pp(I,J)+Depyy_pp(I,J))) )
!                rhsx   = rhsx + 1.0D0/4.0D0*G*( 4.0D0*Depth_pp(I,J)*( 
!      &                   +Depx_pp(I,J)*Etax_pp(I,J))+Eta_pp(I,J)*(Depx_pp(I,J)**2 
!      &                   +2.0D0*Depth_pp(I,J)*(Depxx_pp(I,J))) )
!                rhsy   = rhsy + 1.0D0/4.0D0*G*( 4.0D0*Depth_pp(I,J)*(Depy_pp(I,J)*Etay_pp(I,J) 
!      &                   )+Eta_pp(I,J)*(Depy_pp(I,J)**2 
!      &                   +2.0D0*Depth_pp(I,J)*(Depyy_pp(I,J))) )
!       !..................Nonlinear Varying Bathy
!                rhs(id) = rhs(id) + 1.0D0/6.0D0*G*Eta_pp(I,J)*( 4.0D0*Depy_pp(I,J)*Etay_pp(I,J) 
!      &                   +4.0D0*Depx_pp(I,J)*Etax_pp(I,J)+3.0D0*Eta_pp(I,J)*(Depxx_pp(I,J)+Depyy_pp(I,J)) )
!                rhsx    = rhsx + 1.0D0/6.0D0*G*Eta_pp(I,J)*(  
!      &                   +4.0D0*Depx_pp(I,J)*Etax_pp(I,J)+3.0D0*Eta_pp(I,J)*(Depxx_pp(I,J)) )
!                rhsy    = rhsy + 1.0D0/6.0D0*G*Eta_pp(I,J)*( 4.0D0*Depy_pp(I,J)*Etay_pp(I,J) 
!      &                  +3.0D0*Eta_pp(I,J)*(Depyy_pp(I,J)) )                  
!       !..................Velocity Flatbed
!                rhs(id) = rhs(id) + 5.0D0/3.0D0*(Depth_pp(I,J)+Eta_pp(I,J))**2*(Ux_pp(I,J)**2+Vy_pp(I,J)**2) 
!      &                   +5.0D0/3.0D0*(Depth_pp(I,J)+Eta_pp(I,J))**2*(Ux_pp(I,J)*Vy_pp(I,J)+Uy_pp(I,J)*Vx_pp(I,J))
!                rhsx    = rhsx + 5.0D0/3.0D0*(Depth_pp(I,J)+Eta_pp(I,J))**2*(Ux_pp(I,J)**2) 
!      &                   +5.0D0/3.0D0*(Depth_pp(I,J)+Eta_pp(I,J))**2*(Ux_pp(I,J)*Vy_pp(I,J)+Uy_pp(I,J)*Vx_pp(I,J))/2.d0
!                rhsy    = rhsy + 5.0D0/3.0D0*(Depth_pp(I,J)+Eta_pp(I,J))**2*(Vy_pp(I,J)**2) 
!      &                   +5.0D0/3.0D0*(Depth_pp(I,J)+Eta_pp(I,J))**2*(Ux_pp(I,J)*Vy_pp(I,J)+Uy_pp(I,J)*Vx_pp(I,J))/2.d0                  
!       !..................Velocity Bathy
!                rhs(id) = rhs(id) - 5.0D0/4.0D0*(Depth_pp(I,J)+Eta_pp(I,J))*( Depxx_pp(I,J)*U_pp(I,J)**2
!      &                       + Depyy_pp(I,J)*V_pp(I,J)**2+2.0D0*Depxy_pp(I,J)*V_pp(I,J)*U_pp(I,J) )
!                rhsx    = rhsx - 5.0D0/4.0D0*(Depth_pp(I,J)+Eta_pp(I,J))*( Depxx_pp(I,J)*U_pp(I,J)**2 
!      &                       + 2.0D0*Depxy_pp(I,J)*V_pp(I,J)*U_pp(I,J)/2.d0 )
!                rhsy    = rhsy - 5.0D0/4.0D0*(Depth_pp(I,J)+Eta_pp(I,J))*(  
!      &                       + Depyy_pp(I,J)*V_pp(I,J)**2+2.0D0*Depxy_pp(I,J)*V_pp(I,J)*U_pp(I,J)/2.d0 )                      
!               ENDIF
!             
!             p1sol(id) = p1old(id) + (p1old(id)-p1oldest(id))*DT
!       !       write(99,'(i,4f32.16)') id, p1sol(id), rhs(id), rhsx, rhsy
!             ENDDO
!           ENDDO
!       !     close(99)
!           
!       !     p1sol = p1old
!           
!           RETURN ;
!         END SUBROUTINE GETRHS
!       !==================================================================================!
!         !c Second order FD
!         !c   Generate CSR ---> fdmat
!         SUBROUTINE FDDISCRETE(  )
!           IMPLICIT NONE
!           
!           INTEGER  :: I,J
!           REAL(SZ) :: val(5)
!           INTEGER  :: col_ind(5),iarr(5)
!           INTEGER  :: nz
!           REAL(SZ) :: mydx2,mydy2,mydx1,mydy1
!           INTEGER  :: ind_beg,ind_end,id
!           REAL(SZ) :: mP1yy,mP1xx,mP1y,mP1x,mP1
!           
!           mydx1 = 1.0D0/DX
!           mydy1 = 1.0D0/DY
!           mydx2 = 1.0D0/DX/DX
!           mydy2 = 1.0D0/DY/DY
!           
!           
!       !     open(unit=98,file='pp_mat_coeff.txt',status='replace')
!           ind_beg = 1
!           DO J = 1,NY
!             DO I = 1,NX
!               
!               id = (J-Nghost-1)*NX+(I-Nghost)
!       !.......Column indices
!               col_ind(1) = ((J-1)-Nghost-1)*NX + (I-Nghost)
!               col_ind(2) = (J-Nghost-1)*NX + ((I-1)-Nghost)
!               col_ind(3) = (J-Nghost-1)*NX + (I-Nghost)
!               col_ind(4) = (J-Nghost-1)*NX + ((I+1)-Nghost)
!               col_ind(5) = ((J+1)-Nghost-1)*NX   + (I-Nghost)
!       !.......Index values (based on Poisson equation)
!       !.......Linear Flatbed
!               mP1xx = -1.0D0/2.0D0*Depth_pp(I,J)**2
!               mP1yy = -1.0D0/2.0D0*Depth_pp(I,J)**2
!               mP1x  = 0.0D0
!               mP1y  = 0.0D0
!               mP1   = 5.0D0/4.0D0
!       !.......NonLinear Flatbed
!               mP1xx = mP1xx - 1.0D0/2.0D0*Eta_pp(I,J)*(2.0D0*Depth_pp(I,J)+Eta_pp(I,J))
!               mP1yy = mP1yy - 1.0D0/2.0D0*Eta_pp(I,J)*(2.0D0*Depth_pp(I,J)+Eta_pp(I,J))
!               mP1x  = mP1x  - 2.0D0/3.0D0*Etax_pp(I,J)*(Depth_pp(I,J)+Eta_pp(I,J))
!               mP1y  = mP1y  - 2.0D0/3.0D0*Etay_pp(I,J)*(Depth_pp(I,J)+Eta_pp(I,J))
!               mP1   = mP1   + 1.0D0/3.0D0*( 
!      &                        2.0D0*Etay_pp(I,J)**2-Depth_pp(I,J)*(Etaxx_pp(I,J)+Etayy_pp(I,J)) 
!      &                       +2.0D0*Etax_pp(I,J)**2-Eta_pp(I,J) * (Etaxx_pp(I,J)+Etayy_pp(I,J)) 
!      &                                       )
!       !.......Linear Bathy
!       !         mP1xx = mP1xx
!       !         mP1yy = mP1yy
!               mP1x  = mP1x  - 1.0D0/4.0D0*Depx_pp(I,J)*Depth_pp(I,J)
!               mP1y  = mP1y  - 1.0D0/4.0D0*Depy_pp(I,J)*Depth_pp(I,J)
!               mP1   = mP1   + 1.0D0/4.0D0*(Depx_pp(I,J)**2+Depy_pp(I,J)**2 
!      &                          +2.0D0*Depth_pp(I,J)*(Depxx_pp(I,J)+Depyy_pp(I,J)))
!       !.......NonLinear Bathy
!       !         mP1xx = mP1xx
!       !         mP1yy = mP1yy
!               mP1x  = mP1x  - 1.0D0/4.0D0*Depx_pp(I,J)*Eta_pp(I,J)
!               mP1y  = mP1y  - 1.0D0/4.0D0*Depy_pp(I,J)*Eta_pp(I,J)
!               mP1   = mP1   + 1.0D0/6.0D0*( -2.0D0*(Depx_pp(I,J)*Etax_pp(I,J)+Depy_pp(I,J)*Etay_pp(I,J)) 
!      &                           +3.0D0*Eta_pp(I,J)*(Depxx_pp(I,J)+Depyy_pp(I,J)) )
!               
!       !.......Finalize        
!               val(1) = mP1yy*mydy2 - mP1y*mydy1/2.0D0
!               
!               val(2) = mP1xx*mydx2 - mP1x*mydx1/2.0D0
!               
!               val(3) = -2.0D0*mP1xx*mydx2 
!      &                  -2.0D0*mP1yy*mydy2 
!      &                  + mP1
!              
!               val(4) = mP1xx*mydx2 + mP1x*mydx1/2.0D0
!               
!               val(5) = mP1yy*mydy2 + mP1y*mydy1/2.0D0
!               
!       !         write(98,'(5f32.16)') mP1xx, mP1yy, mP1x, mP1y, mP1
!       !.......Populate the index values
!               IF (PERIODIC) THEN
!                 IF (wetdry_pp(I,J).EQ.0) THEN
!                    nz = 1
!                    val(3) = 1.0D0
!                    iarr(1:nz) = (/ 3 /) ;
!                 ELSEIF ( I > 1 .AND. I < NX .AND. J > 1 .AND. J < NY ) THEN
!                    iarr(1:5) = (/ 1, 2, 3, 4, 5 /) ; 
!                    nz = 5 ; 
!                 ELSEIF ( I > 1 .AND. I < NX ) THEN
!                    nz = 5 ; 
!                    IF ( J == 1 ) THEN
!       !                 val(5) = val(5)+val(1) ;
!                       col_ind(1) = ((NY)-Nghost-1)*NX + (I-Nghost)
!                       iarr(1:nz) = (/ 1, 2, 3, 4, 5 /) ;
!                    ELSE
!       !                 val(1) = val(5)+val(1) ;
!                       col_ind(5) = ((1)-Nghost-1)*NX   + (I-Nghost)
!                       iarr(1:nz) = (/ 1, 2, 3, 4, 5 /) ; 
!                    END IF
!                 ELSEIF (  J > 1 .AND. J < NY ) THEN
!                    nz = 4 ; 
!                    IF ( I == 1 ) THEN
!                       val(4) = val(4)+val(2) ;
!                       iarr(1:nz) = (/ 1, 3, 4, 5 /) ;
!                    ELSE
!                       val(2) = val(4)+val(2) ;
!                       iarr(1:nz) = (/ 1, 2, 3, 5 /) ; 
!                    END IF
!                 ELSE
!                    nz = 4 ;
!                    IF ( I == 1  ) THEN
!                       IF ( J == 1 ) THEN
!       !                    val(5) = val(5)+val(1) ;
!                          val(4) = val(4)+val(2) ;
!                          col_ind(1) = ((NY)-Nghost-1)*NX + (I-Nghost)
!                          iarr(1:nz) = (/ 1, 3, 4, 5 /) ;
!                       ELSE
!       !                    val(1) = val(5)+val(1) ;
!                          val(4) = val(4)+val(2) ;
!                          col_ind(5) = ((1)-Nghost-1)*NX   + (I-Nghost)
!                          iarr(1:nz) = (/ 1, 3, 4, 5 /) ; 
!                       END IF
!                    ELSEIF ( I == NX ) THEN
!                       IF ( J == 1 ) THEN
!       !                    val(5) = val(5)+val(1) ;
!                          val(2) = val(4)+val(2) ;
!                          col_ind(1) = ((NY)-Nghost-1)*NX + (I-Nghost)
!                          iarr(1:nz) = (/ 1, 2, 3, 5 /) ;
!                       ELSE
!       !                    val(1) = val(5)+val(1) ;
!                          val(2) = val(4)+val(2) ;
!                          col_ind(5) = ((1)-Nghost-1)*NX   + (I-Nghost)
!                          iarr(1:nz) = (/ 1, 2, 3, 5 /) ; 
!                       END IF
!                    END IF 
!                 END IF
!               ELSE
!                 IF (wetdry_pp(I,J).EQ.0) THEN
!                    nz = 1
!                    val(3) = 1.0D0
!                    iarr(1:nz) = (/ 3 /) ;
!                 ELSEIF ( I > 1 .AND. I < NX .AND. J > 1 .AND. J < NY ) THEN
!                    iarr(1:5) = (/ 1, 2, 3, 4, 5 /) ; 
!                    nz = 5 ; 
!                 ELSEIF ( I > 1 .AND. I < NX ) THEN
!                    nz = 4 ; 
!                    IF ( J == 1 ) THEN
!                       val(5) = val(5)+val(1) ;
!                       iarr(1:nz) = (/ 2, 3, 4, 5 /) ;
!                    ELSE
!                       val(1) = val(5)+val(1) ;
!                       iarr(1:nz) = (/ 1, 2, 3, 4 /) ; 
!                    END IF
!                 ELSEIF (  J > 1 .AND. J < NY ) THEN
!                    nz = 4 ; 
!                    IF ( I == 1 ) THEN
!                       val(4) = val(4)+val(2) ;
!                       iarr(1:nz) = (/ 1, 3, 4, 5 /) ;
!                    ELSE
!                       val(2) = val(4)+val(2) ;
!                       iarr(1:nz) = (/ 1, 2, 3, 5 /) ; 
!                    END IF
!                 ELSE
!                    nz = 3 ;
!                    IF ( I == 1  ) THEN
!                       IF ( J == 1 ) THEN
!                          val(5) = val(5)+val(1) ;
!                          val(4) = val(4)+val(2) ;
!                          iarr(1:nz) = (/ 3, 4, 5 /) ;
!                       ELSE
!                          val(1) = val(5)+val(1) ;
!                          val(4) = val(4)+val(2) ;
!                          iarr(1:nz) = (/ 1, 3, 4 /) ; 
!                       END IF
!                    ELSEIF ( I == NX ) THEN
!                       IF ( J == 1 ) THEN
!                          val(5) = val(5)+val(1) ;
!                          val(2) = val(4)+val(2) ;
!                          iarr(1:nz) = (/ 2, 3, 5 /) ;
!                       ELSE
!                          val(1) = val(5)+val(1) ;
!                          val(2) = val(4)+val(2) ;
!                          iarr(1:nz) = (/ 1, 2, 3 /) ; 
!                       END IF
!                    END IF 
!                 END IF
!               END IF    
!               ind_end = ind_beg + nz
!                 
!               fdmat.val(ind_beg:(ind_end-1)) = val(iarr(1:nz)) ;
!               fdmat.idv1(ind_beg:(ind_end-1)) = col_ind(iarr(1:nz)) ;
!       
!               fdmat.idv2(id) = ind_beg ;
!                 
!               ind_beg = ind_end ;
!               
!             ENDDO
!           ENDDO
!       !     close(98)
!           
!           fdmat.idv2(id+1) = ind_beg ;
!           
!           open( unit = 10, file = 'valpp.crs' ) ; 
!           open( unit = 20, file = 'colpp.crs' ) ; 
!           open( unit = 30, file = 'ptrpp.crs' ) ; 
!           DO I = 1, ubound(fdmat.val,1) 
!              WRITE(10,*) fdmat.val(i) ;
!              WRITE(20,*) fdmat.idv1(i) ; 
!           END DO
!           close(10) ;
!           close(20) ;
!         
!           DO I = 1, id + 1
!              WRITE(30,*) fdmat.idv2(i) ; 
!           END DO
!           close(30)
!           
!           RETURN ;
!         END SUBROUTINE FDDISCRETE
!       !==================================================================================!
!         SUBROUTINE GETPDPB( )
!           IMPLICIT NONE
!           INTEGER :: I,J
!           REAL(SZ) :: P2,Hloc
!           
!           DO J = 1,NY
!             DO I = 1,NX
!               IF (wetdry_pp(I,J).EQ.0) THEN
!                 pd_pp(I,J) = 0.0D0
!                 pb_pp(I,J) = 0.0D0
!               ELSE
!       !............Linear Flatbed
!               P2 = 3.0D0/4.0D0*(p1_pp(I,J)-G*Eta_pp(I,J))
!       !............Linear Bathy
!               P2 = P2 +3.0D0/4.0D0*(Depx_pp(I,J)**2+Depy_pp(I,J)**2)*(p1_pp(I,J)-G*Eta_pp(I,J)) 
!      &             -3.0D0/4.0D0*Depth_pp(I,J)*(Depx_pp(I,J)*p1x_pp(I,J)+Depy_pp(I,J)*p1y_pp(I,J))
!       !............NonLinear Bathy
!               P2 = P2 - 3.0D0/4.0D0*Eta_pp(I,J)*(Depx_pp(I,J)*p1x_pp(I,J)+Depy_pp(I,J)*p1y_pp(I,J))
!       ! !............Velocity Bathy
!               P2 = P2 + 3.0D0/4.0D0*(Depth_pp(I,J)+Eta_pp(I,J))*( Depxx_pp(I,J)*U_pp(I,J)**2 
!      &                +Depyy_pp(I,J)*V_pp(I,J)**2+2.0D0*U_pp(I,J)*V_pp(I,J)*Depxy_pp(I,J) )
!               
!               
!       !.......Pd
!       !................Linear
!               pd_pp(I,J) = 1.0D0/2.0D0*Depth_pp(I,J)*(p1_pp(I,J)-G*Eta_pp(I,J))
!       !................Nonlinear
!               pd_pp(I,J) = pd_pp(I,J) + 1.0D0/2.0D0*Eta_pp(I,J)*(p1_pp(I,J)-G*Eta_pp(I,J))
!       !.......and Pb 
!       !................Linear
!               pb_pp(I,J) = p1_pp(I,J)-1.0D0/3.0D0*P2-G*Eta_pp(I,J)
!               ENDIF
!             END DO
!           END DO    
!           
!           RETURN ;
!         END SUBROUTINE GETPDPB
!       !==================================================================================!
!         SUBROUTINE DERIVATIVE_X(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,WETDRY,DX,F,Fx)
!           IMPLICIT NONE
!           INTEGER :: I,J
!           INTEGER,INTENT(IN) :: Mloc,Nloc,Ibeg,Iend,Jbeg,Jend
!           INTEGER,INTENT(IN) :: WETDRY(Mloc,Nloc)
!           REAL(SZ),INTENT(IN) :: F(Mloc,Nloc),DX
!           REAL(SZ),INTENT(OUT) :: Fx(Mloc,Nloc)
!           
!           DO J = Jbeg,Jend
!             DO I = Ibeg,Iend      
!               IF (WETDRY(I,J).EQ.0) THEN
!                 FX(I,J) = 0.D0
!               ELSEIF ( I > Ibeg .AND. I < Iend .AND. J > Jbeg .AND. J < Jend ) THEN  ! INTERIOR
!                 FX(I,J) = (F(I+1,J) - F(I-1,J) )/2.D0/DX
!               ELSEIF ( I > Ibeg .AND. I < Iend ) THEN 
!                 IF ( J == Jbeg ) THEN ! SOUTH
!                   FX(I,J) = (F(I+1,J) - F(I-1,J) )/2.D0/DX
!                 ELSE ! NORTH
!                   FX(I,J) = (F(I+1,J) - F(I-1,J) )/2.D0/DX
!                 END IF
!               ELSEIF (  J > Jbeg .AND. J < Jend ) THEN
!                 IF ( I == Ibeg ) THEN ! EAST
!                   FX(I,J) = 0.D0!(F(I+1,J) - F(I+1,J) )/2.D0/DX
!                 ELSE ! WEST
!                   FX(I,J) = 0.D0!(F(I-1,J) - F(I-1,J) )/2.D0/DX
!                 END IF
!               ELSE ! CORNER
!                 IF ( I == Ibeg  ) THEN
!                   IF ( J == Jbeg ) THEN
!                     FX(I,J) = 0.D0!(F(I+1,J) - F(I+1,J) )/2.D0/DX
!                   ELSE
!                     FX(I,J) = 0.D0!(F(I+1,J) - F(I+1,J) )/2.D0/DX
!                   END IF
!                 ELSEIF ( I == Iend ) THEN
!                   IF ( J == Jbeg ) THEN
!                     FX(I,J) = 0.D0!(F(I-1,J) - F(I-1,J) )/2.D0/DX
!                   ELSE
!                     FX(I,J) = 0.D0!(F(I-1,J) - F(I-1,J) )/2.D0/DX 
!                   END IF
!                 END IF 
!               END IF
!             END DO
!           END DO
!           
!           
!           RETURN 
!         END SUBROUTINE DERIVATIVE_X
!       !==================================================================================!
!         SUBROUTINE DERIVATIVE_Y(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,WETDRY,DY,F,FY)
!           IMPLICIT NONE
!           INTEGER :: I,J
!           INTEGER,INTENT(IN) :: Mloc,Nloc,Ibeg,Iend,Jbeg,Jend
!           INTEGER,INTENT(IN) :: WETDRY(Mloc,Nloc)
!           REAL(SZ),INTENT(IN) :: F(Mloc,Nloc),DY
!           REAL(SZ),INTENT(OUT) :: FY(Mloc,Nloc)
!           
!           DO J = Jbeg,Jend
!             DO I = Ibeg,Iend      
!               IF (WETDRY(I,J).EQ.0) THEN
!                 FY(I,J) = 0.D0
!               ELSEIF ( I > Ibeg .AND. I < Iend .AND. J > Jbeg .AND. J < Jend ) THEN  ! INTERIOR
!                 FY(I,J) = (F(I,J+1) - F(I,J-1) )/2.D0/DY
!               ELSEIF ( I > Ibeg .AND. I < Iend ) THEN 
!                 IF ( J == Jbeg ) THEN ! SOUTH
!                   FY(I,J) = 0.D0!(F(I,J+1) - F(I,J+1) )/2.D0/DX
!                 ELSE ! NORTH
!                   FY(I,J) = 0.D0!(F(I,J-1) - F(I,J-1) )/2.D0/DX
!                 END IF
!               ELSEIF (  J > Jbeg .AND. J < Jend ) THEN
!                 IF ( I == Ibeg ) THEN ! EAST
!                   FY(I,J) = (F(I,J+1) - F(I,J-1) )/2.D0/DY
!                 ELSE ! WEST
!                   FY(I,J) = (F(I,J+1) - F(I,J-1) )/2.D0/DY
!                 END IF
!               ELSE ! CORNER
!                 IF ( I == Ibeg  ) THEN
!                   IF ( J == Jbeg ) THEN
!                     FY(I,J) = 0.D0!(F(I,J+1) - F(I,J+1) )/2.D0/DX
!                   ELSE
!                     FY(I,J) = 0.D0!(F(I,J-1) - F(I,J-1) )/2.D0/DX
!                   END IF
!                 ELSEIF ( I == Iend ) THEN
!                   IF ( J == Jbeg ) THEN
!                     FY(I,J) = 0.D0!(F(I,J+1) - F(I,J+1) )/2.D0/DX
!                   ELSE
!                     FY(I,J) = 0.D0!(F(I,J-1) - F(I,J-1) )/2.D0/DX 
!                   END IF
!                 END IF 
!               END IF
!             END DO
!           END DO
!           
!           
!           RETURN 
!         END SUBROUTINE DERIVATIVE_Y
!       !==================================================================================!
!         SUBROUTINE DERIVATIVE_XX(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,WETDRY,DX,F,Fxx)
!           IMPLICIT NONE
!           INTEGER :: I,J
!           INTEGER,INTENT(IN) :: Mloc,Nloc,Ibeg,Iend,Jbeg,Jend
!           INTEGER,INTENT(IN) :: WETDRY(Mloc,Nloc)
!           REAL(SZ),INTENT(IN) :: F(Mloc,Nloc),DX
!           REAL(SZ),INTENT(OUT) :: Fxx(Mloc,Nloc)
!           
!           DO J = Jbeg,Jend
!             DO I = Ibeg,Iend      
!               IF (WETDRY(I,J).EQ.0) THEN
!                 FXX(I,J) = 0.D0
!               ELSEIF ( I > Ibeg .AND. I < Iend .AND. J > Jbeg .AND. J < Jend ) THEN  ! INTERIOR
!                 FXX(I,J) = (F(I+1,J) + F(I-1,J) - 2.D0*F(I,J) )/DX/DX
!               ELSEIF ( I > Ibeg .AND. I < Iend ) THEN 
!                 IF ( J == Jbeg ) THEN ! SOUTH
!                   FXX(I,J) = (F(I+1,J) + F(I-1,J) - 2.D0*F(I,J) )/DX/DX
!                 ELSE ! NORTH
!                   FXX(I,J) = (F(I+1,J) + F(I-1,J) - 2.D0*F(I,J) )/DX/DX
!                 END IF
!               ELSEIF (  J > Jbeg .AND. J < Jend ) THEN
!                 IF ( I == Ibeg ) THEN ! EAST
!                   FXX(I,J) = (F(I+1,J) + F(I+1,J) - 2.D0*F(I,J) )/DX/DX
!                 ELSE ! WEST
!                   FXX(I,J) = (F(I-1,J) + F(I-1,J) - 2.D0*F(I,J) )/DX/DX
!                 END IF
!               ELSE ! CORNER
!                 IF ( I == Ibeg  ) THEN
!                   IF ( J == Jbeg ) THEN
!                     FXX(I,J) = (F(I+1,J) + F(I+1,J) - 2.D0*F(I,J) )/DX/DX
!                   ELSE
!                     FXX(I,J) = (F(I+1,J) + F(I+1,J) - 2.D0*F(I,J) )/DX/DX
!                   END IF
!                 ELSEIF ( I == Iend ) THEN
!                   IF ( J == Jbeg ) THEN
!                     FXX(I,J) = (F(I-1,J) + F(I-1,J) - 2.D0*F(I,J) )/DX/DX
!                   ELSE
!                     FXX(I,J) = (F(I-1,J) + F(I-1,J) - 2.D0*F(I,J) )/DX/DX 
!                   END IF
!                 END IF 
!               END IF
!             END DO
!           END DO
!           
!           
!           RETURN 
!         END SUBROUTINE DERIVATIVE_XX
!       !==================================================================================!
!         SUBROUTINE DERIVATIVE_YY(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,WETDRY,DY,F,Fyy)
!           IMPLICIT NONE
!           INTEGER :: I,J
!           INTEGER,INTENT(IN) :: Mloc,Nloc,Ibeg,Iend,Jbeg,Jend
!           INTEGER,INTENT(IN) :: WETDRY(Mloc,Nloc)
!           REAL(SZ),INTENT(IN) :: F(Mloc,Nloc),DY
!           REAL(SZ),INTENT(OUT) :: Fyy(Mloc,Nloc)
!           
!           DO J = Jbeg,Jend
!             DO I = Ibeg,Iend      
!               IF (WETDRY(I,J).EQ.0) THEN
!                 FYY(I,J) = 0.D0
!               ELSEIF ( I > Ibeg .AND. I < Iend .AND. J > Jbeg .AND. J < Jend ) THEN  ! INTERIOR
!                 FYY(I,J) = ( F(I,J+1) + F(I,J-1) - 2.D0*F(I,J) )/DY/DY
!               ELSEIF ( I > Ibeg .AND. I < Iend ) THEN 
!                 IF ( J == Jbeg ) THEN ! SOUTH
!                   FYY(I,J) = ( F(I,J+1) + F(I,J+1) - 2.D0*F(I,J) )/DY/DY
!                 ELSE ! NORTH
!                   FYY(I,J) = ( F(I,J-1) + F(I,J-1) - 2.D0*F(I,J) )/DY/DY
!                 END IF
!               ELSEIF (  J > Jbeg .AND. J < Jend ) THEN
!                 IF ( I == Ibeg ) THEN ! EAST
!                   FYY(I,J) = ( F(I,J+1) + F(I,J-1) - 2.D0*F(I,J) )/DY/DY
!                 ELSE ! WEST
!                   FYY(I,J) = ( F(I,J+1) + F(I,J-1) - 2.D0*F(I,J) )/DY/DY
!                 END IF
!               ELSE ! CORNER
!                 IF ( I == Ibeg  ) THEN
!                   IF ( J == Jbeg ) THEN
!                     FYY(I,J) = ( F(I,J+1) + F(I,J+1) - 2.D0*F(I,J) )/DY/DY
!                   ELSE
!                     FYY(I,J) = ( F(I,J-1) + F(I,J-1) - 2.D0*F(I,J) )/DY/DY
!                   END IF
!                 ELSEIF ( I == Iend ) THEN
!                   IF ( J == Jbeg ) THEN
!                     FYY(I,J) = ( F(I,J+1) + F(I,J+1) - 2.D0*F(I,J) )/DY/DY
!                   ELSE
!                     FYY(I,J) = ( F(I,J-1) + F(I,J-1) - 2.D0*F(I,J) )/DY/DY 
!                   END IF
!                 END IF 
!               END IF
!             END DO
!           END DO
!           
!           
!           RETURN 
!         END SUBROUTINE DERIVATIVE_YY
!       !==================================================================================!
!         SUBROUTINE DERIVATIVE_XY(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,WETDRY,DX,DY,F,FXY)
!           IMPLICIT NONE
!           INTEGER :: I,J
!           INTEGER,INTENT(IN) :: Mloc,Nloc,Ibeg,Iend,Jbeg,Jend
!           INTEGER,INTENT(IN) :: WETDRY(Mloc,Nloc)
!           REAL(SZ),INTENT(IN) :: F(Mloc,Nloc),DX,DY
!           REAL(SZ),INTENT(OUT) :: FXY(Mloc,Nloc)
!           
!           DO J = Jbeg,Jend
!             DO I = Ibeg,Iend      
!               IF (WETDRY(I,J).EQ.0) THEN
!                 FXY(I,J) = 0.D0
!               ELSEIF ( I > Ibeg .AND. I < Iend .AND. J > Jbeg .AND. J < Jend ) THEN  ! INTERIOR
!                 FXY(I,J) = (F(I+1,J) + F(I-1,J) - 2.D0*F(I,J) )/DX/DX
!               ELSE ! BOUNDARY
!                 FXY(I,J) = 0.D0 
!               END IF
!             END DO
!           END DO
!           
!           
!           RETURN 
!         END SUBROUTINE DERIVATIVE_XY
! !       !==================================================================================!
! !       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !       SUBROUTINE TESTILUT(nnz,fdmat)
! !       
! !       USE AddSparskit
! !       
! !       USE sizes, ONLY: SZ
! !       USE global, ONLY: mnp
! !       
! !       IMPLICIT NONE
! !       TYPE (spmat),intent(in) :: fdmat
! !       TYPE (spmat):: fdmat1
! !       TYPE (spmat), pointer:: ilumat1
! !       
! ! !       REAL(SZ) :: DropTol
! ! !       INTEGER :: ifil,initmem
! !       INTEGER :: i
! !       INTEGER,INTENT(in) :: nnz
! !       
! !       
! !       fdmat1.str = 'csr'  ;
! !       fdmat1.ndim = mnp ; !c number of rows
! !         
! !       allocate(fdmat1.val(7*mnp)) ;
! !       allocate(fdmat1.idv1(7*mnp)) ;
! !       allocate(fdmat1.idv2(mnp+1)) ;
! !       
! !       open( unit = 10, file = 'valnew.crs' , action='READ') ; 
! !       open( unit = 20, file = 'colnew.crs' , action='READ') ; 
! !       open( unit = 30, file = 'ptr.crs' , action='READ') ;
! !       do i = 1,nnz
! !         read(10,*) fdmat1.val(i)
! !         read(20,'(i)') fdmat1.idv1(i)
! !       end do
! ! 
! !       do i = 1,mnp+1
! !         read(30,'(i)') fdmat1.idv2(i)
! !       end do
! !       close(10)
! !       close(20)
! !       close(30)
! !       
! !       do i = 1,10
! ! !         write(*,*) fdmat.val(i), fdmat1.val(i)
! ! !         write(*,*) fdmat.idv1(i), fdmat1.idv1(i)
! ! !         write(*,*) fdmat.idv2(i), fdmat1.idv2(i)
! !           write(*,*) fdmat1.idv2(mnp+2-i)
! !       end do
! !       
! !       !.Precondition Matrix  
! !       DropTol = 0.00001D0 ; 
! !       ifil = 20 ; !c ifil -- 3*bandwidth
! !       initmem = 20 ; 
! !       ilumat1 => WrpILUT( fdmat1, DROPTOL, IFILOP = IFIL, INITMEMOP = initmem ) ;
! !       print*, 'goodbye'
! !       
! !       END SUBROUTINE TESTILUT
      END SUBROUTINE pressure_WASUPP_CG
