!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    SUBROUTINE pressure_WASSUP 
!    Aaron Donahue and Jessica Meixner 
!    March 2015
!    Calculates dynamic pressure 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      
      
      SUBROUTINE pressure_WASUPP
      USE AddSPARSKIT
      USE sizes,    ONLY: SZ 
      USE global,   ONLY: G,DT 
      USE pressure_module, ONLY: NumCols, NumRows, 
     &                            Eta_pp,U_pp,V_pp,Depth_pp, 
     &                            Etax_pp,Etay_pp,Ux_pp,Uy_pp,Vx_pp,Vy_pp,
     &                            Etaxx_pp,Etayy_pp,Depx_pp,Depy_pp,Depxx_pp,Depyy_pp,Depxy_pp,
     &                            p1_pp,p1x_pp,p1y_pp,p1old,p1oldest, 
     &                            totaldepth_pp, wetdry_pp,
     &                            DX,DY,pd_pp,pb_pp !!! ASD Need to add these
      
      IMPLICIT NONE
      
      INTEGER                       :: VTYPE,NX,NY,NT,Nghost 
      LOGICAL                       :: PERIODIC=.FALSE. 
      REAL(SZ),ALLOCATABLE           :: rhs(:),p1sol(:)
      REAL(SZ)                       :: DropTol,RESID,infor(4)
      INTEGER                       :: ifil,initmem,itmax,infoi(4)
      !c FD matrix
      TYPE (spmat):: fdmat
      TYPE (spmat), pointer:: ilumat
       
      
        
      !.Allocate local variables
      
       NX = NumRows 
       NY = NumCols 
       Nghost = 0
       NT = NX*NY
       ALLOCATE(rhs(NT),p1sol(NT))
      
        fdmat.str = 'csr'  ;
        fdmat.ndim = NT ; !c number of rows 
      
        allocate(fdmat.val(5*NT)) ;  !c estimate nnz for second order matrix
        allocate(fdmat.idv1(5*NT)) ; !c column index (col_ind)
        allocate(fdmat.idv2(NT+1)) ; !c row index (row_prt)
      !.First Derivatives  
        CALL DERIVATIVE_X(NX,NY,1,NX,1,NY,wetdry_pp,DX,U_pp,Ux_pp)
        CALL DERIVATIVE_X(NX,NY,1,NX,1,NY,wetdry_pp,DX,V_pp,Vx_pp)
        CALL DERIVATIVE_X(NX,NY,1,NX,1,NY,wetdry_pp,DX,Eta_pp,Etax_pp)
        CALL DERIVATIVE_X(NX,NY,1,NX,1,NY,wetdry_pp,DX,Depth_pp,Depx_pp)
      
        
        CALL DERIVATIVE_Y(NX,NY,1,NX,1,NY,wetdry_pp,DY,U_pp,Uy_pp)
        CALL DERIVATIVE_Y(NX,NY,1,NX,1,NY,wetdry_pp,DY,V_pp,Vy_pp)
        CALL DERIVATIVE_Y(NX,NY,1,NX,1,NY,wetdry_pp,DY,Eta_pp,Etay_pp)
        CALL DERIVATIVE_Y(NX,NY,1,NX,1,NY,wetdry_pp,DY,Depth_pp,Depy_pp)
      
      !.Second Derivatives
        CALL DERIVATIVE_XX(NX,NY,1,NX,1,NY,wetdry_pp,DX,Eta_pp,Etaxx_pp)
        CALL DERIVATIVE_YY(NX,NY,1,NX,1,NY,wetdry_pp,DY,Eta_pp,Etayy_pp)
        CALL DERIVATIVE_XX(NX,NY,1,NX,1,NY,wetdry_pp,DX,Depth_pp,Depxx_pp)
        CALL DERIVATIVE_YY(NX,NY,1,NX,1,NY,wetdry_pp,DY,Depth_pp,Depyy_pp)
        CALL DERIVATIVE_XY(NX,NY,1,NX,1,NY,wetdry_pp,DX,DY,Depth_pp,Depxy_pp)
      
        
      !.Determine RHS vector of Pressure Poisson problem for P1
        CALL GETRHS()
      
      !.Construct the left hand side Matrix for solving P1
        CALL FDDISCRETE()
      
      !.Precondition Matrix  
        DropTol = 0.0001D0 ; 
        ifil = 20 ; !c ifil -- 3*bandwidth
        initmem = 20 ; 
        ilumat => WrpILUT( fdmat, DROPTOL, IFILOP = IFIL, INITMEMOP = initmem ) ;
      
      !.Solve Matrix using BiCGSTAB 
      !   p1sol = 0.0D0
        RESID = 1.d-4 ; 
        ITMAX = NT/4 ;
        p1sol = WrpSPBCCGSTAB0( rhs, p1sol, fdmat, ilumat,RESID, ITMAX, INFOR, INFOI ) ;
        p1oldest = p1old
        p1old = p1sol
        
      !   print*, 'Iter = ',infoi(2)
      
      !.Find Pd and Pb from p1
        p1_pp = 0.0D0
        p1_pp(1:NX,1:NY) = reshape(p1sol,(/NX,NY/))
        VTYPE = 1
      !   CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,p1,VTYPE,PERIODIC)
        CALL DERIVATIVE_X(NX,NY,1,NX,1,NY,wetdry_pp,DX,p1_pp,p1x_pp)
        CALL DERIVATIVE_Y(NX,NY,1,NX,1,NY,wetdry_pp,DY,p1_pp,p1y_pp)
        CALL GETPDPB()
        
        ! Assign Ghost Node Values
        VTYPE = 1
      !   CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Pd,VTYPE,PERIODIC)
      !   CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Pb,VTYPE,PERIODIC)
      
        ! Deallocate arrays
        CALL EmptySpmat( fdmat );
        CALL EmptySpmat( ilumat ) ;
        deallocate( ilumat ) ;
        deallocate(rhs,p1sol) ;
      CONTAINS
      !==================================================================================!
        SUBROUTINE GETRHS( )
          IMPLICIT NONE
          INTEGER  :: I,J,id
          
          REAL(SZ) :: rhsx, rhsy
      
      !     open(unit=99,file='p1sol.txt',status='replace')
          DO J = 1,NY
            DO I = 1,NX
              
              id = (J-Nghost-1)*(NX)+(I-Nghost)
      
              IF (wetdry_pp(I,J).EQ.0) THEN
                rhs(id) = G*Eta_pp(I,J)
              ELSE
      !..................Linear Flatbed
               rhs(id) = 1.0D0/12.0D0*G*( 15.0D0*Eta_pp(I,J)+4.0D0*Depth_pp(I,J)**2* 
     &                    (Etaxx_pp(I,J)+Etayy_pp(I,J)) )
               rhsx = 1.0D0/12.0D0*G*( 15.0D0*Eta_pp(I,J)/2.d0+4.0D0*Depth_pp(I,J)**2* 
     &                    (Etaxx_pp(I,J)) )
               rhsy = 1.0D0/12.0D0*G*( 15.0D0*Eta_pp(I,J)/2.d0+4.0D0*Depth_pp(I,J)**2* 
     &                    (Etayy_pp(I,J)) )
      !..................NonLinear Flatbed
               rhs(id) = rhs(id) + 1.0D0/3.0D0*G*Depth_pp(I,J)*(-2.0D0*(Etax_pp(I,J)**2+Etay_pp(I,J)**2) 
     &                   +Eta_pp(I,J)*(Etaxx_pp(I,J)+Etayy_pp(I,J)))
               rhsx    = rhsx + 1.0D0/3.0D0*G*Depth_pp(I,J)*(-2.0D0*(Etax_pp(I,J)**2) 
     &                   +Eta_pp(I,J)*(Etaxx_pp(I,J)))
               rhsy    = rhsy + 1.0D0/3.0D0*G*Depth_pp(I,J)*(-2.0D0*(Etay_pp(I,J)**2) 
     &                   +Eta_pp(I,J)*(Etayy_pp(I,J)))
      !..................Linear Varying Bathy
               rhs(id) = rhs(id) + 1.0D0/4.0D0*G*( 4.0D0*Depth_pp(I,J)*(Depy_pp(I,J)*Etay_pp(I,J) 
     &                   +Depx_pp(I,J)*Etax_pp(I,J))+Eta_pp(I,J)*(Depy_pp(I,J)**2+Depx_pp(I,J)**2 
     &                   +2.0D0*Depth_pp(I,J)*(Depxx_pp(I,J)+Depyy_pp(I,J))) )
               rhsx   = rhsx + 1.0D0/4.0D0*G*( 4.0D0*Depth_pp(I,J)*( 
     &                   +Depx_pp(I,J)*Etax_pp(I,J))+Eta_pp(I,J)*(Depx_pp(I,J)**2 
     &                   +2.0D0*Depth_pp(I,J)*(Depxx_pp(I,J))) )
               rhsy   = rhsy + 1.0D0/4.0D0*G*( 4.0D0*Depth_pp(I,J)*(Depy_pp(I,J)*Etay_pp(I,J) 
     &                   )+Eta_pp(I,J)*(Depy_pp(I,J)**2 
     &                   +2.0D0*Depth_pp(I,J)*(Depyy_pp(I,J))) )
      !..................Nonlinear Varying Bathy
               rhs(id) = rhs(id) + 1.0D0/6.0D0*G*Eta_pp(I,J)*( 4.0D0*Depy_pp(I,J)*Etay_pp(I,J) 
     &                   +4.0D0*Depx_pp(I,J)*Etax_pp(I,J)+3.0D0*Eta_pp(I,J)*(Depxx_pp(I,J)+Depyy_pp(I,J)) )
               rhsx    = rhsx + 1.0D0/6.0D0*G*Eta_pp(I,J)*(  
     &                   +4.0D0*Depx_pp(I,J)*Etax_pp(I,J)+3.0D0*Eta_pp(I,J)*(Depxx_pp(I,J)) )
               rhsy    = rhsy + 1.0D0/6.0D0*G*Eta_pp(I,J)*( 4.0D0*Depy_pp(I,J)*Etay_pp(I,J) 
     &                  +3.0D0*Eta_pp(I,J)*(Depyy_pp(I,J)) )                  
      !..................Velocity Flatbed
               rhs(id) = rhs(id) + 5.0D0/3.0D0*(Depth_pp(I,J)+Eta_pp(I,J))**2*(Ux_pp(I,J)**2+Vy_pp(I,J)**2) 
     &                   +5.0D0/3.0D0*(Depth_pp(I,J)+Eta_pp(I,J))**2*(Ux_pp(I,J)*Vy_pp(I,J)+Uy_pp(I,J)*Vx_pp(I,J))
               rhsx    = rhsx + 5.0D0/3.0D0*(Depth_pp(I,J)+Eta_pp(I,J))**2*(Ux_pp(I,J)**2) 
     &                   +5.0D0/3.0D0*(Depth_pp(I,J)+Eta_pp(I,J))**2*(Ux_pp(I,J)*Vy_pp(I,J)+Uy_pp(I,J)*Vx_pp(I,J))/2.d0
               rhsy    = rhsy + 5.0D0/3.0D0*(Depth_pp(I,J)+Eta_pp(I,J))**2*(Vy_pp(I,J)**2) 
     &                   +5.0D0/3.0D0*(Depth_pp(I,J)+Eta_pp(I,J))**2*(Ux_pp(I,J)*Vy_pp(I,J)+Uy_pp(I,J)*Vx_pp(I,J))/2.d0                  
      !..................Velocity Bathy
               rhs(id) = rhs(id) - 5.0D0/4.0D0*(Depth_pp(I,J)+Eta_pp(I,J))*( Depxx_pp(I,J)*U_pp(I,J)**2
     &                       + Depyy_pp(I,J)*V_pp(I,J)**2+2.0D0*Depxy_pp(I,J)*V_pp(I,J)*U_pp(I,J) )
               rhsx    = rhsx - 5.0D0/4.0D0*(Depth_pp(I,J)+Eta_pp(I,J))*( Depxx_pp(I,J)*U_pp(I,J)**2 
     &                       + 2.0D0*Depxy_pp(I,J)*V_pp(I,J)*U_pp(I,J)/2.d0 )
               rhsy    = rhsy - 5.0D0/4.0D0*(Depth_pp(I,J)+Eta_pp(I,J))*(  
     &                       + Depyy_pp(I,J)*V_pp(I,J)**2+2.0D0*Depxy_pp(I,J)*V_pp(I,J)*U_pp(I,J)/2.d0 )                      
              ENDIF
            
            p1sol(id) = p1old(id) + (p1old(id)-p1oldest(id))*DT
      !       write(99,'(i,4f32.16)') id, p1sol(id), rhs(id), rhsx, rhsy
            ENDDO
          ENDDO
      !     close(99)
          
      !     p1sol = p1old
          
          RETURN ;
        END SUBROUTINE GETRHS
      !==================================================================================!
        !c Second order FD
        !c   Generate CSR ---> fdmat
        SUBROUTINE FDDISCRETE(  )
          IMPLICIT NONE
          
          INTEGER  :: I,J
          REAL(SZ) :: val(5)
          INTEGER  :: col_ind(5),iarr(5)
          INTEGER  :: nz
          REAL(SZ) :: mydx2,mydy2,mydx1,mydy1
          INTEGER  :: ind_beg,ind_end,id
          REAL(SZ) :: mP1yy,mP1xx,mP1y,mP1x,mP1
          
          mydx1 = 1.0D0/DX
          mydy1 = 1.0D0/DY
          mydx2 = 1.0D0/DX/DX
          mydy2 = 1.0D0/DY/DY
          
          
      !     open(unit=98,file='pp_mat_coeff.txt',status='replace')
          ind_beg = 1
          DO J = 1,NY
            DO I = 1,NX
              
              id = (J-Nghost-1)*NX+(I-Nghost)
      !.......Column indices
              col_ind(1) = ((J-1)-Nghost-1)*NX + (I-Nghost)
              col_ind(2) = (J-Nghost-1)*NX + ((I-1)-Nghost)
              col_ind(3) = (J-Nghost-1)*NX + (I-Nghost)
              col_ind(4) = (J-Nghost-1)*NX + ((I+1)-Nghost)
              col_ind(5) = ((J+1)-Nghost-1)*NX   + (I-Nghost)
      !.......Index values (based on Poisson equation)
      !.......Linear Flatbed
              mP1xx = -1.0D0/2.0D0*Depth_pp(I,J)**2
              mP1yy = -1.0D0/2.0D0*Depth_pp(I,J)**2
              mP1x  = 0.0D0
              mP1y  = 0.0D0
              mP1   = 5.0D0/4.0D0
      !.......NonLinear Flatbed
              mP1xx = mP1xx - 1.0D0/2.0D0*Eta_pp(I,J)*(2.0D0*Depth_pp(I,J)+Eta_pp(I,J))
              mP1yy = mP1yy - 1.0D0/2.0D0*Eta_pp(I,J)*(2.0D0*Depth_pp(I,J)+Eta_pp(I,J))
              mP1x  = mP1x  - 2.0D0/3.0D0*Etax_pp(I,J)*(Depth_pp(I,J)+Eta_pp(I,J))
              mP1y  = mP1y  - 2.0D0/3.0D0*Etay_pp(I,J)*(Depth_pp(I,J)+Eta_pp(I,J))
              mP1   = mP1   + 1.0D0/3.0D0*( 
     &                        2.0D0*Etay_pp(I,J)**2-Depth_pp(I,J)*(Etaxx_pp(I,J)+Etayy_pp(I,J)) 
     &                       +2.0D0*Etax_pp(I,J)**2-Eta_pp(I,J) * (Etaxx_pp(I,J)+Etayy_pp(I,J)) 
     &                                       )
      !.......Linear Bathy
      !         mP1xx = mP1xx
      !         mP1yy = mP1yy
              mP1x  = mP1x  - 1.0D0/4.0D0*Depx_pp(I,J)*Depth_pp(I,J)
              mP1y  = mP1y  - 1.0D0/4.0D0*Depy_pp(I,J)*Depth_pp(I,J)
              mP1   = mP1   + 1.0D0/4.0D0*(Depx_pp(I,J)**2+Depy_pp(I,J)**2 
     &                          +2.0D0*Depth_pp(I,J)*(Depxx_pp(I,J)+Depyy_pp(I,J)))
      !.......NonLinear Bathy
      !         mP1xx = mP1xx
      !         mP1yy = mP1yy
              mP1x  = mP1x  - 1.0D0/4.0D0*Depx_pp(I,J)*Eta_pp(I,J)
              mP1y  = mP1y  - 1.0D0/4.0D0*Depy_pp(I,J)*Eta_pp(I,J)
              mP1   = mP1   + 1.0D0/6.0D0*( -2.0D0*(Depx_pp(I,J)*Etax_pp(I,J)+Depy_pp(I,J)*Etay_pp(I,J)) 
     &                           +3.0D0*Eta_pp(I,J)*(Depxx_pp(I,J)+Depyy_pp(I,J)) )
              
      !.......Finalize        
              val(1) = mP1yy*mydy2 - mP1y*mydy1/2.0D0
              
              val(2) = mP1xx*mydx2 - mP1x*mydx1/2.0D0
              
              val(3) = -2.0D0*mP1xx*mydx2 
     &                  -2.0D0*mP1yy*mydy2 
     &                  + mP1
             
              val(4) = mP1xx*mydx2 + mP1x*mydx1/2.0D0
              
              val(5) = mP1yy*mydy2 + mP1y*mydy1/2.0D0
              
      !         write(98,'(5f32.16)') mP1xx, mP1yy, mP1x, mP1y, mP1
      !.......Populate the index values
              IF (PERIODIC) THEN
                IF (wetdry_pp(I,J).EQ.0) THEN
                   nz = 1
                   val(3) = 1.0D0
                   iarr(1:nz) = (/ 3 /) ;
                ELSEIF ( I > 1 .AND. I < NX .AND. J > 1 .AND. J < NY ) THEN
                   iarr(1:5) = (/ 1, 2, 3, 4, 5 /) ; 
                   nz = 5 ; 
                ELSEIF ( I > 1 .AND. I < NX ) THEN
                   nz = 5 ; 
                   IF ( J == 1 ) THEN
      !                 val(5) = val(5)+val(1) ;
                      col_ind(1) = ((NY)-Nghost-1)*NX + (I-Nghost)
                      iarr(1:nz) = (/ 1, 2, 3, 4, 5 /) ;
                   ELSE
      !                 val(1) = val(5)+val(1) ;
                      col_ind(5) = ((1)-Nghost-1)*NX   + (I-Nghost)
                      iarr(1:nz) = (/ 1, 2, 3, 4, 5 /) ; 
                   END IF
                ELSEIF (  J > 1 .AND. J < NY ) THEN
                   nz = 4 ; 
                   IF ( I == 1 ) THEN
                      val(4) = val(4)+val(2) ;
                      iarr(1:nz) = (/ 1, 3, 4, 5 /) ;
                   ELSE
                      val(2) = val(4)+val(2) ;
                      iarr(1:nz) = (/ 1, 2, 3, 5 /) ; 
                   END IF
                ELSE
                   nz = 4 ;
                   IF ( I == 1  ) THEN
                      IF ( J == 1 ) THEN
      !                    val(5) = val(5)+val(1) ;
                         val(4) = val(4)+val(2) ;
                         col_ind(1) = ((NY)-Nghost-1)*NX + (I-Nghost)
                         iarr(1:nz) = (/ 1, 3, 4, 5 /) ;
                      ELSE
      !                    val(1) = val(5)+val(1) ;
                         val(4) = val(4)+val(2) ;
                         col_ind(5) = ((1)-Nghost-1)*NX   + (I-Nghost)
                         iarr(1:nz) = (/ 1, 3, 4, 5 /) ; 
                      END IF
                   ELSEIF ( I == NX ) THEN
                      IF ( J == 1 ) THEN
      !                    val(5) = val(5)+val(1) ;
                         val(2) = val(4)+val(2) ;
                         col_ind(1) = ((NY)-Nghost-1)*NX + (I-Nghost)
                         iarr(1:nz) = (/ 1, 2, 3, 5 /) ;
                      ELSE
      !                    val(1) = val(5)+val(1) ;
                         val(2) = val(4)+val(2) ;
                         col_ind(5) = ((1)-Nghost-1)*NX   + (I-Nghost)
                         iarr(1:nz) = (/ 1, 2, 3, 5 /) ; 
                      END IF
                   END IF 
                END IF
              ELSE
                IF (wetdry_pp(I,J).EQ.0) THEN
                   nz = 1
                   val(3) = 1.0D0
                   iarr(1:nz) = (/ 3 /) ;
                ELSEIF ( I > 1 .AND. I < NX .AND. J > 1 .AND. J < NY ) THEN
                   iarr(1:5) = (/ 1, 2, 3, 4, 5 /) ; 
                   nz = 5 ; 
                ELSEIF ( I > 1 .AND. I < NX ) THEN
                   nz = 4 ; 
                   IF ( J == 1 ) THEN
                      val(5) = val(5)+val(1) ;
                      iarr(1:nz) = (/ 2, 3, 4, 5 /) ;
                   ELSE
                      val(1) = val(5)+val(1) ;
                      iarr(1:nz) = (/ 1, 2, 3, 4 /) ; 
                   END IF
                ELSEIF (  J > 1 .AND. J < NY ) THEN
                   nz = 4 ; 
                   IF ( I == 1 ) THEN
                      val(4) = val(4)+val(2) ;
                      iarr(1:nz) = (/ 1, 3, 4, 5 /) ;
                   ELSE
                      val(2) = val(4)+val(2) ;
                      iarr(1:nz) = (/ 1, 2, 3, 5 /) ; 
                   END IF
                ELSE
                   nz = 3 ;
                   IF ( I == 1  ) THEN
                      IF ( J == 1 ) THEN
                         val(5) = val(5)+val(1) ;
                         val(4) = val(4)+val(2) ;
                         iarr(1:nz) = (/ 3, 4, 5 /) ;
                      ELSE
                         val(1) = val(5)+val(1) ;
                         val(4) = val(4)+val(2) ;
                         iarr(1:nz) = (/ 1, 3, 4 /) ; 
                      END IF
                   ELSEIF ( I == NX ) THEN
                      IF ( J == 1 ) THEN
                         val(5) = val(5)+val(1) ;
                         val(2) = val(4)+val(2) ;
                         iarr(1:nz) = (/ 2, 3, 5 /) ;
                      ELSE
                         val(1) = val(5)+val(1) ;
                         val(2) = val(4)+val(2) ;
                         iarr(1:nz) = (/ 1, 2, 3 /) ; 
                      END IF
                   END IF 
                END IF
              END IF    
              ind_end = ind_beg + nz
                
              fdmat.val(ind_beg:(ind_end-1)) = val(iarr(1:nz)) ;
              fdmat.idv1(ind_beg:(ind_end-1)) = col_ind(iarr(1:nz)) ;
      
              fdmat.idv2(id) = ind_beg ;
                
              ind_beg = ind_end ;
              
            ENDDO
          ENDDO
      !     close(98)
          
          fdmat.idv2(id+1) = ind_beg ;
          
      !     open( unit = 10, file = 'val.crs' ) ; 
      !     open( unit = 20, file = 'col.crs' ) ; 
      !     open( unit = 30, file = 'ptr.crs' ) ; 
      !     DO I = 1, ubound(fdmat.val,1) 
      !        WRITE(10,*) fdmat.val(i) ;
      !        WRITE(20,*) fdmat.idv1(i) ; 
      !     END DO
      !     close(10) ;
      !     close(20) ;
      !   
      !     DO I = 1, id + 1
      !        WRITE(30,*) fdmat.idv2(i) ; 
      !     END DO
      !     close(30)
          
          RETURN ;
        END SUBROUTINE FDDISCRETE
      !==================================================================================!
        SUBROUTINE GETPDPB( )
          IMPLICIT NONE
          INTEGER :: I,J
          REAL(SZ) :: P2,Hloc
          
          DO J = 1,NY
            DO I = 1,NX
              IF (wetdry_pp(I,J).EQ.0) THEN
                pd_pp(I,J) = 0.0D0
                pb_pp(I,J) = 0.0D0
              ELSE
      !............Linear Flatbed
              P2 = 3.0D0/4.0D0*(p1_pp(I,J)-G*Eta_pp(I,J))
      !............Linear Bathy
              P2 = P2 +3.0D0/4.0D0*(Depx_pp(I,J)**2+Depy_pp(I,J)**2)*(p1_pp(I,J)-G*Eta_pp(I,J)) 
     &             -3.0D0/4.0D0*Depth_pp(I,J)*(Depx_pp(I,J)*p1x_pp(I,J)+Depy_pp(I,J)*p1y_pp(I,J))
      !............NonLinear Bathy
              P2 = P2 - 3.0D0/4.0D0*Eta_pp(I,J)*(Depx_pp(I,J)*p1x_pp(I,J)+Depy_pp(I,J)*p1y_pp(I,J))
      ! !............Velocity Bathy
              P2 = P2 + 3.0D0/4.0D0*(Depth_pp(I,J)+Eta_pp(I,J))*( Depxx_pp(I,J)*U_pp(I,J)**2 
     &                +Depyy_pp(I,J)*V_pp(I,J)**2+2.0D0*U_pp(I,J)*V_pp(I,J)*Depxy_pp(I,J) )
              
              
      !.......Pd
      !................Linear
              pd_pp(I,J) = 1.0D0/2.0D0*Depth_pp(I,J)*(p1_pp(I,J)-G*Eta_pp(I,J))
      !................Nonlinear
              pd_pp(I,J) = pd_pp(I,J) + 1.0D0/2.0D0*Eta_pp(I,J)*(p1_pp(I,J)-G*Eta_pp(I,J))
      !.......and Pb 
      !................Linear
              pb_pp(I,J) = p1_pp(I,J)-1.0D0/3.0D0*P2-G*Eta_pp(I,J)
              ENDIF
            END DO
          END DO    
          
          RETURN ;
        END SUBROUTINE GETPDPB
      !==================================================================================!
        SUBROUTINE DERIVATIVE_X(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,WETDRY,DX,F,Fx)
          IMPLICIT NONE
          INTEGER :: I,J
          INTEGER,INTENT(IN) :: Mloc,Nloc,Ibeg,Iend,Jbeg,Jend
          INTEGER,INTENT(IN) :: WETDRY(Mloc,Nloc)
          REAL(SZ),INTENT(IN) :: F(Mloc,Nloc),DX
          REAL(SZ),INTENT(OUT) :: Fx(Mloc,Nloc)
          
          DO J = Jbeg,Jend
            DO I = Ibeg,Iend      
              IF (WETDRY(I,J).EQ.0) THEN
                FX(I,J) = 0.D0
              ELSEIF ( I > Ibeg .AND. I < Iend .AND. J > Jbeg .AND. J < Jend ) THEN  ! INTERIOR
                FX(I,J) = (F(I+1,J) - F(I-1,J) )/2.D0/DX
              ELSEIF ( I > Ibeg .AND. I < Iend ) THEN 
                IF ( J == Jbeg ) THEN ! SOUTH
                  FX(I,J) = (F(I+1,J) - F(I-1,J) )/2.D0/DX
                ELSE ! NORTH
                  FX(I,J) = (F(I+1,J) - F(I-1,J) )/2.D0/DX
                END IF
              ELSEIF (  J > Jbeg .AND. J < Jend ) THEN
                IF ( I == Ibeg ) THEN ! EAST
                  FX(I,J) = 0.D0!(F(I+1,J) - F(I+1,J) )/2.D0/DX
                ELSE ! WEST
                  FX(I,J) = 0.D0!(F(I-1,J) - F(I-1,J) )/2.D0/DX
                END IF
              ELSE ! CORNER
                IF ( I == Ibeg  ) THEN
                  IF ( J == Jbeg ) THEN
                    FX(I,J) = 0.D0!(F(I+1,J) - F(I+1,J) )/2.D0/DX
                  ELSE
                    FX(I,J) = 0.D0!(F(I+1,J) - F(I+1,J) )/2.D0/DX
                  END IF
                ELSEIF ( I == Iend ) THEN
                  IF ( J == Jbeg ) THEN
                    FX(I,J) = 0.D0!(F(I-1,J) - F(I-1,J) )/2.D0/DX
                  ELSE
                    FX(I,J) = 0.D0!(F(I-1,J) - F(I-1,J) )/2.D0/DX 
                  END IF
                END IF 
              END IF
            END DO
          END DO
          
          
          RETURN 
        END SUBROUTINE DERIVATIVE_X
      !==================================================================================!
        SUBROUTINE DERIVATIVE_Y(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,WETDRY,DY,F,FY)
          IMPLICIT NONE
          INTEGER :: I,J
          INTEGER,INTENT(IN) :: Mloc,Nloc,Ibeg,Iend,Jbeg,Jend
          INTEGER,INTENT(IN) :: WETDRY(Mloc,Nloc)
          REAL(SZ),INTENT(IN) :: F(Mloc,Nloc),DY
          REAL(SZ),INTENT(OUT) :: FY(Mloc,Nloc)
          
          DO J = Jbeg,Jend
            DO I = Ibeg,Iend      
              IF (WETDRY(I,J).EQ.0) THEN
                FY(I,J) = 0.D0
              ELSEIF ( I > Ibeg .AND. I < Iend .AND. J > Jbeg .AND. J < Jend ) THEN  ! INTERIOR
                FY(I,J) = (F(I,J+1) - F(I,J-1) )/2.D0/DY
              ELSEIF ( I > Ibeg .AND. I < Iend ) THEN 
                IF ( J == Jbeg ) THEN ! SOUTH
                  FY(I,J) = 0.D0!(F(I,J+1) - F(I,J+1) )/2.D0/DX
                ELSE ! NORTH
                  FY(I,J) = 0.D0!(F(I,J-1) - F(I,J-1) )/2.D0/DX
                END IF
              ELSEIF (  J > Jbeg .AND. J < Jend ) THEN
                IF ( I == Ibeg ) THEN ! EAST
                  FY(I,J) = (F(I,J+1) - F(I,J-1) )/2.D0/DY
                ELSE ! WEST
                  FY(I,J) = (F(I,J+1) - F(I,J-1) )/2.D0/DY
                END IF
              ELSE ! CORNER
                IF ( I == Ibeg  ) THEN
                  IF ( J == Jbeg ) THEN
                    FY(I,J) = 0.D0!(F(I,J+1) - F(I,J+1) )/2.D0/DX
                  ELSE
                    FY(I,J) = 0.D0!(F(I,J-1) - F(I,J-1) )/2.D0/DX
                  END IF
                ELSEIF ( I == Iend ) THEN
                  IF ( J == Jbeg ) THEN
                    FY(I,J) = 0.D0!(F(I,J+1) - F(I,J+1) )/2.D0/DX
                  ELSE
                    FY(I,J) = 0.D0!(F(I,J-1) - F(I,J-1) )/2.D0/DX 
                  END IF
                END IF 
              END IF
            END DO
          END DO
          
          
          RETURN 
        END SUBROUTINE DERIVATIVE_Y
      !==================================================================================!
        SUBROUTINE DERIVATIVE_XX(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,WETDRY,DX,F,Fxx)
          IMPLICIT NONE
          INTEGER :: I,J
          INTEGER,INTENT(IN) :: Mloc,Nloc,Ibeg,Iend,Jbeg,Jend
          INTEGER,INTENT(IN) :: WETDRY(Mloc,Nloc)
          REAL(SZ),INTENT(IN) :: F(Mloc,Nloc),DX
          REAL(SZ),INTENT(OUT) :: Fxx(Mloc,Nloc)
          
          DO J = Jbeg,Jend
            DO I = Ibeg,Iend      
              IF (WETDRY(I,J).EQ.0) THEN
                FXX(I,J) = 0.D0
              ELSEIF ( I > Ibeg .AND. I < Iend .AND. J > Jbeg .AND. J < Jend ) THEN  ! INTERIOR
                FXX(I,J) = (F(I+1,J) + F(I-1,J) - 2.D0*F(I,J) )/DX/DX
              ELSEIF ( I > Ibeg .AND. I < Iend ) THEN 
                IF ( J == Jbeg ) THEN ! SOUTH
                  FXX(I,J) = (F(I+1,J) + F(I-1,J) - 2.D0*F(I,J) )/DX/DX
                ELSE ! NORTH
                  FXX(I,J) = (F(I+1,J) + F(I-1,J) - 2.D0*F(I,J) )/DX/DX
                END IF
              ELSEIF (  J > Jbeg .AND. J < Jend ) THEN
                IF ( I == Ibeg ) THEN ! EAST
                  FXX(I,J) = (F(I+1,J) + F(I+1,J) - 2.D0*F(I,J) )/DX/DX
                ELSE ! WEST
                  FXX(I,J) = (F(I-1,J) + F(I-1,J) - 2.D0*F(I,J) )/DX/DX
                END IF
              ELSE ! CORNER
                IF ( I == Ibeg  ) THEN
                  IF ( J == Jbeg ) THEN
                    FXX(I,J) = (F(I+1,J) + F(I+1,J) - 2.D0*F(I,J) )/DX/DX
                  ELSE
                    FXX(I,J) = (F(I+1,J) + F(I+1,J) - 2.D0*F(I,J) )/DX/DX
                  END IF
                ELSEIF ( I == Iend ) THEN
                  IF ( J == Jbeg ) THEN
                    FXX(I,J) = (F(I-1,J) + F(I-1,J) - 2.D0*F(I,J) )/DX/DX
                  ELSE
                    FXX(I,J) = (F(I-1,J) + F(I-1,J) - 2.D0*F(I,J) )/DX/DX 
                  END IF
                END IF 
              END IF
            END DO
          END DO
          
          
          RETURN 
        END SUBROUTINE DERIVATIVE_XX
      !==================================================================================!
        SUBROUTINE DERIVATIVE_YY(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,WETDRY,DY,F,Fyy)
          IMPLICIT NONE
          INTEGER :: I,J
          INTEGER,INTENT(IN) :: Mloc,Nloc,Ibeg,Iend,Jbeg,Jend
          INTEGER,INTENT(IN) :: WETDRY(Mloc,Nloc)
          REAL(SZ),INTENT(IN) :: F(Mloc,Nloc),DY
          REAL(SZ),INTENT(OUT) :: Fyy(Mloc,Nloc)
          
          DO J = Jbeg,Jend
            DO I = Ibeg,Iend      
              IF (WETDRY(I,J).EQ.0) THEN
                FYY(I,J) = 0.D0
              ELSEIF ( I > Ibeg .AND. I < Iend .AND. J > Jbeg .AND. J < Jend ) THEN  ! INTERIOR
                FYY(I,J) = ( F(I,J+1) + F(I,J-1) - 2.D0*F(I,J) )/DY/DY
              ELSEIF ( I > Ibeg .AND. I < Iend ) THEN 
                IF ( J == Jbeg ) THEN ! SOUTH
                  FYY(I,J) = ( F(I,J+1) + F(I,J+1) - 2.D0*F(I,J) )/DY/DY
                ELSE ! NORTH
                  FYY(I,J) = ( F(I,J-1) + F(I,J-1) - 2.D0*F(I,J) )/DY/DY
                END IF
              ELSEIF (  J > Jbeg .AND. J < Jend ) THEN
                IF ( I == Ibeg ) THEN ! EAST
                  FYY(I,J) = ( F(I,J+1) + F(I,J-1) - 2.D0*F(I,J) )/DY/DY
                ELSE ! WEST
                  FYY(I,J) = ( F(I,J+1) + F(I,J-1) - 2.D0*F(I,J) )/DY/DY
                END IF
              ELSE ! CORNER
                IF ( I == Ibeg  ) THEN
                  IF ( J == Jbeg ) THEN
                    FYY(I,J) = ( F(I,J+1) + F(I,J+1) - 2.D0*F(I,J) )/DY/DY
                  ELSE
                    FYY(I,J) = ( F(I,J-1) + F(I,J-1) - 2.D0*F(I,J) )/DY/DY
                  END IF
                ELSEIF ( I == Iend ) THEN
                  IF ( J == Jbeg ) THEN
                    FYY(I,J) = ( F(I,J+1) + F(I,J+1) - 2.D0*F(I,J) )/DY/DY
                  ELSE
                    FYY(I,J) = ( F(I,J-1) + F(I,J-1) - 2.D0*F(I,J) )/DY/DY 
                  END IF
                END IF 
              END IF
            END DO
          END DO
          
          
          RETURN 
        END SUBROUTINE DERIVATIVE_YY
      !==================================================================================!
        SUBROUTINE DERIVATIVE_XY(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,WETDRY,DX,DY,F,FXY)
          IMPLICIT NONE
          INTEGER :: I,J
          INTEGER,INTENT(IN) :: Mloc,Nloc,Ibeg,Iend,Jbeg,Jend
          INTEGER,INTENT(IN) :: WETDRY(Mloc,Nloc)
          REAL(SZ),INTENT(IN) :: F(Mloc,Nloc),DX,DY
          REAL(SZ),INTENT(OUT) :: FXY(Mloc,Nloc)
          
          DO J = Jbeg,Jend
            DO I = Ibeg,Iend      
              IF (WETDRY(I,J).EQ.0) THEN
                FXY(I,J) = 0.D0
              ELSEIF ( I > Ibeg .AND. I < Iend .AND. J > Jbeg .AND. J < Jend ) THEN  ! INTERIOR
                FXY(I,J) = (F(I+1,J) + F(I-1,J) - 2.D0*F(I,J) )/DX/DX
              ELSE ! BOUNDARY
                FXY(I,J) = 0.D0 
              END IF
            END DO
          END DO
          
          
          RETURN 
        END SUBROUTINE DERIVATIVE_XY
      !==================================================================================!
      END SUBROUTINE pressure_WASUPP
