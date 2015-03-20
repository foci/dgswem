!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!      This module is for adding non-hydrostatic pressure to the calculation using 
!        Green-Nhagdhi Pressure Poisson Finite Difference formulations
!
!       
!      pressure_module.F
!
!      Jessica Meixner and Aaron Donahue  Feb 2015 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! ASD - Added projection of DG solution to CG solution
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      MODULE pressure_module

      USE SIZES,    ONLY: SZ
      USE global,   ONLY: NP,mnp
      use dg,       ONLY: nagp
      USE AddSPARSKIT

!      INTEGER LengthX, LengthY 
      INTEGER NumCols, NumRows 


      !Variables that will be passed from DG-SWEM to Pressure Poisson 
      REAL(SZ), ALLOCATABLE :: totaldepth_pp(:,:) !total depth	
      REAL(SZ), ALLOCATABLE :: eta_pp(:,:)        !water surface elevation
      REAL(SZ), ALLOCATABLE :: u_pp(:,:)          !u- current in x direction
      REAL(SZ), ALLOCATABLE :: v_pp(:,:)          !v- current in y direction
      REAL(SZ), ALLOCATABLE :: depth_pp(:,:)       !bathymetry 
      INTEGER,  ALLOCATABLE :: wetdry_pp(:,:)     !wet (1) dry (0) flag 

      REAL(SZ), ALLOCATABLE :: Etax_pp(:,:),Etay_pp(:,:),Etaxx_pp(:,:),Etayy_pp(:,:) !derivatives of surface elevation 
      REAL(SZ), ALLOCATABLE :: Ux_pp(:,:),Uy_pp(:,:),Vx_pp(:,:),Vy_pp(:,:)        !derivatives of current 
      REAL(SZ), ALLOCATABLE :: Depx_pp(:,:),Depy_pp(:,:),Depxx_pp(:,:),Depyy_pp(:,:),Depxy_pp(:,:)   !derivatevs of bathymetry 
      REAL(SZ), ALLOCATABLE :: p1_pp(:,:),p1x_pp(:,:),p1y_pp(:,:)

      REAL(SZ)              :: DX,DY
      REAL(SZ), ALLOCATABLE :: pd_pp(:,:),pb_pp(:,:)
      REAL(SZ), ALLOCATABLE :: p1old(:),p1oldest(:) 
      
      INTEGER,ALLOCATABLE   :: csrloc(:,:,:)
      TYPE (spmat)          :: pp_dg_cg_mat
      TYPE (spmat), pointer :: pp_dg_cg_ilumat
      
      REAL(SZ),allocatable  :: phinodal(:,:)
      
      INTEGER               :: totalnz
      
      REAL(SZ),allocatable  :: pp_p1(:),pp_p2(:),pp_pd(:),pp_pb(:)
      
            
     
      CONTAINS 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE PP_DG_TO_CG
      ! This subroutine calculate the CG Global Mass matrix to be used in the
      ! projection of any DG solution onto a CG framework.
      USE global, ONLY: MNP,NM,NE,neimin,neimax,neitab,nneigh,areas
      USE sizes,  ONLY: SZ
      USE dg,     ONLY: XAGP,YAGP,WAGP,NAGP
      USE AddSPARSKIT
      
      IMPLICIT NONE            
      INTEGER               :: i,j,l,k,pt,pad,a(3)
      INTEGER	            :: ifil,initmem
      REAL(SZ)              :: DropTol
      REAL(SZ),ALLOCATABLE  :: rhs(:),MCG(:,:)
      
      REAL(SZ),ALLOCATABLE  :: matv(:,:)
      
      
      
      ! Determine the value of the linear nodal basis functions at the quadrature points.
      pad = 1      
      allocate(rhs(np),MCG(3,3))
      allocate(phinodal(3,NAGP(pad))) 
      
      allocate(pp_p1(mnp),pp_p2(mnp))
      allocate(pp_pd(mnp),pp_pb(mnp))
            
      do pt = 1,NAGP(pad)
        phinodal(1,pt) = -1.d0/2.d0*(xagp(pt,pad)+yagp(pt,pad))
        phinodal(2,pt) =  1.d0/2.d0*(1.d0+xagp(pt,pad))
        phinodal(3,pt) =  1.d0/2.d0*(1.d0+yagp(pt,pad))
      end do
      
      ! Determine the elemental Mass Matrix
      MCG(:,:) = 0.d0
      do j = 1,3
        do i = 1,3
          do pt = 1,NAGP(pad)
            MCG(i,j) = MCG(i,j) + phinodal(i,pt)*phinodal(j,pt)*wagp(pt,pad)
          end do
        end do  
      end do  
      
      print*, 'MCG Matrix:'
      open(unit=11,file='MCGlocal.txt',status='replace')
      do i = 1,3
        write(*,'(3F32.16)') mcg(i,:)
        write(11,'(3F32.16)') mcg(i,:)
      end do
      
      ! Construct the global CG Mass matrix
      allocate(matv(mnp,neimax))
      matv(:,:) = 0.d0 ! start all values at 0
      
      ! Determine the csr starting point for each node "row" in domain
      pp_dg_cg_mat.str = 'csr'  ;
      pp_dg_cg_mat.ndim = mnp ; !c number of rows
      allocate(pp_dg_cg_mat.idv2(mnp+1)) ; !c row index (row_prt)
      pp_dg_cg_mat.idv2(:) = 0
      pp_dg_cg_mat.idv2(1) = 1
      totalnz = nneigh(1);
      do i = 2,mnp
        pp_dg_cg_mat.idv2(i) =  pp_dg_cg_mat.idv2(i-1)+nneigh(i-1)
        totalnz = totalnz + nneigh(i-1)
      end do
      pp_dg_cg_mat.idv2(mnp+1) = totalnz+1
      
!       allocate(pp_dg_cg_mat.idv1(totalnz)) ; !c column index (col_ind)      
      allocate(pp_dg_cg_mat.idv1(totalnz+2)) ; !c column index (col_ind)
      pp_dg_cg_mat.idv1(:) = 0
      allocate(csrloc(ne,3,3)) ! Determine the connectivity of each node with respects to CSR format
      do l = 1,ne
        ! Node numbers associated eith element
        a(1) = NM(l,1)
        a(2) = NM(l,2)
        a(3) = NM(l,3)
        ! Connectivity associated with node a(j)
        do j = 1,3
        do i = 1,3
          do k = 1,neimax ! cycle through all connected nodes
            if (neitab(a(j),k).eq.a(i)) then
              csrloc(l,j,i) = k+pp_dg_cg_mat.idv2(a(j))-1 ! designate the csr column format location for this connection
              pp_dg_cg_mat.idv1(csrloc(l,j,i)) = a(i)
              ! Note that when building the CSR matrix
              !  l: Designates the current element
              !  j: Designates the row
              !  i: Designates the csr row
              goto 900
            end if           
          end do
          print*, 'Error establishing connectivity table for CSR format'
          print*, 'Error encountered for element', l
          call exit          
900       continue
        end do
        end do
      end do
      
      open(unit=10,file='csrloc_test.txt',status='replace')
      open(unit=20,file='nneigh_elem.txt',status='replace')
      open(unit=30,file='elem_table.txt',status='replace')
      do l = 1,ne
        write(30,'(3i)') NM(l,1:3)
        do i = 1,3
          do j = 1,3
          write(10,'(4i)') l, i, j, csrloc(l,j,i)
          end do
        end do        
      end do
      do i = 1,mnp
        write(20,'(3i)') i, nneigh(i), pp_dg_cg_mat.idv2(i)
      end do
      close(10)
      close(20)
      close(30)
      close(11)
      
      
      ! Build global mass matrix in CSR format using csrloc and csrnode
      
      ! Allocate memory for the global CG mass matrix to be used by Sparsekit      
!       allocate(pp_dg_cg_mat.val(totalnz)) ;  !c estimate nnz for second order matrix
      allocate(pp_dg_cg_mat.val(totalnz+2)) ;  !c estimate nnz for second order matrix
      pp_dg_cg_mat.val(:) = 0.d0
      do l = 1,ne
!         a(1) = NM(l,1)
!         a(2) = NM(l,2)
!         a(3) = NM(l,3)
        do j = 1,3
          do i = 1,3
            pp_dg_cg_mat.val(csrloc(l,j,i)) = pp_dg_cg_mat.val(csrloc(l,j,i)) + AREAS(l)/4.d0*MCG(j,i)            
          end do
        end do
      end do
      
      open( unit = 10, file = 'val.crs' ) ; 
      open( unit = 20, file = 'col.crs' ) ; 
      open( unit = 30, file = 'ptr.crs' ) ; 
      DO I = 1, ubound(pp_dg_cg_mat.val,1) 
	    WRITE(10,'(f32.16)') pp_dg_cg_mat.val(i) ;
	    WRITE(20,*) pp_dg_cg_mat.idv1(i) ; 
      END DO
      close(10) ;
      close(20) ;
      
      DO I = 1, mnp + 1
	    WRITE(30,*) pp_dg_cg_mat.idv2(i) ; 
      END DO
      close(30)

      !! .Precondition Matrix  
      DropTol = 0.00001D0 ; 
      ifil = 21 ; !c ifil -- 3*bandwidth
      initmem = 20 ; 
      pp_dg_cg_ilumat => WrpILUT( pp_dg_cg_mat, DROPTOL, IFILOP = IFIL, INITMEMOP = initmem ) ;
      
      deallocate(rhs,MCG)
      deallocate(matv)
      
      
      END SUBROUTINE PP_DG_TO_CG

            
      !==================================================================================!      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE  read_and_initialize_PPvalues 

      USE global,   ONLY: DP 

      IMPLICIT NONE 

      INTEGER i, j, inode


      !read in what is LengthInX, LengthInY 
      !read in what is NumCols, NumRows
      NumCols = 101 
      NumRows = 6
      DX = .1d0 
      DY = 1d0
  
      !read in what are the different boundary types: 
         !only 1 (reflective) boundary type for now. 

      ALLOCATE(totaldepth_pp(NumRows,NumCols)) 
      ALLOCATE(eta_pp(NumRows,NumCols))
      ALLOCATE(u_pp(NumRows,NumCols))
      ALLOCATE(v_pp(NumRows,NumCols))
      ALLOCATE(depth_pp(NumRows,NumCols))
      ALLOCATE(wetdry_pp(NumRows,NumCols))
      ALLOCATE(Etax_pp(NumRows,NumCols),Etay_pp(NumRows,NumCols))
      ALLOCATE(Etaxx_pp(NumRows,NumCols),Etayy_pp(NumRows,NumCols))
      ALLOCATE(Ux_pp(NumRows,NumCols),Uy_pp(NumRows,NumCols))
      ALLOCATE(Vx_pp(NumRows,NumCols),Vy_pp(NumRows,NumCols))      
      ALLOCATE(Depx_pp(NumRows,NumCols),Depy_pp(NumRows,NumCols))
      ALLOCATE(Depxx_pp(NumRows,NumCols),Depyy_pp(NumRows,NumCols))
      ALLOCATE(Depxy_pp(NumRows,NumCols))   
      ALLOCATE(p1_pp(NumRows,NumCols),p1x_pp(NumRows,NumCols),p1y_pp(NumRows,NumCols))
      ALLOCATE(pd_pp(NumRows,NumCols),pb_pp(NumRows,NumCols))
      ALLOCATE(p1old(NumRows*NumCols),p1oldest(NumRows*NumCols))  

 
      !Inititalize values:     
      totaldepth_pp = 0.d0 
      eta_pp = 0.d0 
      u_pp = 0.d0 
      v_pp = 0.d0 
      wetdry_pp = 1
      Etax_pp=0.d0; Etay_pp=0.d0 
      Etaxx_pp=0.d0; Etayy_pp=0.d0
      Ux_pp=0.d0;Uy_pp=0.d0
      Vx_pp=0.d0;Vy_pp=0.d0 
      Depx_pp=0.d0;Depy_pp=0.d0
      Depxx_pp=0.d0; Depyy_pp=0.d0
      Depxy_pp=0.d0
      p1_pp=0.d0;p1x_pp=0.d0;p1y_pp=0.d0
      pd_pp=0.d0;pb_pp=0.d0
      p1old=0.d0;p1oldest=0.d0 
 
 
      inode = 1
      DO i = 1, NumRows
        DO j = 1, NumCols
           depth_pp(i,j) = DP(inode)
           inode = inode + 1
        END DO !i 
      END DO !j 

      END SUBROUTINE  read_and_initialize_PPvalues 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE  from_DG_to_Pressure(jjRK) 

      !Need to obtain u,v,eta,depth and wet/dry at each node.  then need to pass to (i,j) type format 

      USE GLOBAL, ONLY: MNP, EL_COUNT, ELETAB, NM,
     &                  ZE_DG, QX_DG, QY_DG, HB_DG, 
     &                  AREAS, ETA2, UU2, VV2, 
     &                  NLEQ, LEQ, ne, x, y, dp
      USE DG,     ONLY: WDFLG, ZE, QX, QY, HB, DOFH, 
     &                  PHI_CORNER, ph, wagp, nagp, phi_area    

      IMPLICIT NONE 

C.....Declare local variables

      INTEGER, INTENT(In) :: jjRK
      INTEGER I, NO_NBORS, J,K,N1,KK, NBOR_EL, inode 
      REAL(SZ) AREA_SUM, AREA, DEPTH, FH_NL
      REAL(SZ), ALLOCATABLE :: TOTALDEPTH(:)
      INTEGER, ALLOCATABLE :: WDFLG_VERT(:) 
      
      REAL(SZ),ALLOCATABLE :: dgcgeta(:),dgcgu(:),dgcgv(:),dgcgloc(:)
      REAL(SZ)	           :: ze_in,ue_in,ve_in,INFOR(4)
      INTEGER	           :: INFOI(4)
      INTEGER              :: pad, l
      

      ALLOCATE(TOTALDEPTH(MNP))
      ALLOCATE(WDFLG_VERT(MNP))
      WDFLG_VERT = 1 !Start as all being wet, will set to dry if any element connected to vertex is dry 

      !Obtain u,v,eta, depth and wet/dry flag at each node 

      DO I = 1,MNP

        NO_NBORS = EL_COUNT(I)

        AREA_SUM = 0
        DO 333 J = 1,NO_NBORS
          NBOR_EL = ELETAB(I,1+J)

          IF(WDFLG(NBOR_EL).EQ.0) THEN 
            WDFLG_VERT(I) = 0 !if any neighboring elem is dry, consider vertex dry 
            CYCLE ! DON'T COUNT DRY ELEMENTS  sb 02/26/07
          END IF 

          DO K = 1,3
            N1 = NM(NBOR_EL,K)
            IF (N1.EQ.I) THEN
              ZE_DG(J) = ZE(1,NBOR_EL,jjRK)
              QX_DG(J) = QX(1,NBOR_EL,jjRK)
              QY_DG(J) = QY(1,NBOR_EL,jjRK)
              HB_DG(J) = HB(1,NBOR_EL,jjRK)
              DO KK = 2,DOFH
                ZE_DG(J) = ZE_DG(J) + PHI_CORNER(KK,K,ph)*ZE(KK,NBOR_EL,jjRK)
                QX_DG(J) = QX_DG(J) + PHI_CORNER(KK,K,ph)*QX(KK,NBOR_EL,jjRK)
                QY_DG(J) = QY_DG(J) + PHI_CORNER(KK,K,ph)*QY(KK,NBOR_EL,jjRK)
                HB_DG(J) = HB_DG(J) + PHI_CORNER(KK,K,ph)*HB(KK,NBOR_EL,jjRK)
              ENDDO
              AREA = 0.5*AREAS(NBOR_EL)
              AREA_SUM = AREA_SUM + AREA
              GOTO 333
            ENDIF
          ENDDO
 333    CONTINUE

        ETA2(I) = 0.D0
        UU2(I)  = 0.D0
        VV2(I)  = 0.D0
        TOTALDEPTH(I) = 0.d0 


        DO J = 1,NO_NBORS
          NBOR_EL = ELETAB(I,1+J)

          IF(WDFLG(NBOR_EL).EQ.0) CYCLE ! DON'T COUNT DRY ELEMENTS  sb 02/26/07

          AREA = 0.5*AREAS(NBOR_EL)/AREA_SUM
          DEPTH = ZE_DG(J) + HB_DG(J)
          FH_NL = 1.D0/(NLEQ*DEPTH + LEQ)
          ETA2(I) = ETA2(I) + AREA*ZE_DG(J)
          UU2(I)  = UU2(I)  + AREA*QX_DG(J)*FH_NL
          VV2(I)  = VV2(I)  + AREA*QY_DG(J)*FH_NL
          TOTALDEPTH(I) = TOTALDEPTH(I) + AREA*DEPTH  
        ENDDO
      ENDDO



      !now that we have ETA2,UU2,VV2,TOTALDEPTH, and WDFLG_VERT 
      !need to go from a node numbering to a (i,j) vertex numbering PP FD can use.  
      open(unit=10,file='dgfd_ver1.txt',status='replace',action='write')
      inode = 1 
      DO i = 1, NumRows 
        DO j = 1, NumCols 
           totaldepth_pp(i,j) = TOTALDEPTH(inode) 
           eta_pp(i,j) = ETA2(inode) 
           u_pp(i,j) = UU2(inode) 
           v_pp(i,j) = VV2(inode) 
           wetdry_pp(i,j) = WDFLG_VERT(inode)
           inode = inode + 1
           write(10,'(2i,f)') i, j, eta_pp(i,j)
        END DO !i 
      END DO !j                   
      close(10)

            
!       CALL PP_DG_TO_CG
      
      !!!!! ASD ASD ASD !!!!
      pad = 1
      ! Build DG right hand side
      allocate(dgcgeta(mnp),dgcgu(mnp),dgcgv(mnp),dgcgloc(3))
      dgcgeta(:) = 0.d0
      dgcgu(:) = 0.d0
      dgcgv(:) = 0.d0
      do l = 1,ne
        dgcgloc(:) = 0.d0
        ! determine ze at each quad point        
        do i = 1,nagp(pad)
          ze_in = ze(1,l,jjrk)
          ue_in = qx(1,l,jjrk)
          ve_in = qy(1,l,jjrk)
          do kk = 2,dofh
            ze_in = ze_in + ze(kk,l,jjrk)*phi_area(kk,i,pad)
            ue_in = ue_in + qx(kk,l,jjrk)*phi_area(kk,i,pad)
            ve_in = ve_in + qy(kk,l,jjrk)*phi_area(kk,i,pad)
          end do
          
          do kk = 1,3
            dgcgeta(NM(l,kk)) = dgcgeta(NM(l,kk)) + wagp(i,pad)*phinodal(kk,i)*ze_in*areas(l)/4.d0
            dgcgu(NM(l,kk))   = dgcgu(NM(l,kk))   + wagp(i,pad)*phinodal(kk,i)*ue_in*areas(l)/4.d0
            dgcgv(NM(l,kk))   = dgcgv(NM(l,kk))   + wagp(i,pad)*phinodal(kk,i)*ve_in*areas(l)/4.d0
          end do
        end do
      end do
      
      ETA2 = WrpSPBCCGSTAB0( dgcgeta, ETA2, pp_dg_cg_mat, pp_dg_cg_ilumat, 1.d-12, mnp/4, INFOR, INFOI ) ;
      UU2 = WrpSPBCCGSTAB0( dgcgu, UU2, pp_dg_cg_mat, pp_dg_cg_ilumat, 1.d-12, mnp/4, INFOR, INFOI ) ;
      VV2 = WrpSPBCCGSTAB0( dgcgv, VV2, pp_dg_cg_mat, pp_dg_cg_ilumat, 1.d-12, mnp/4, INFOR, INFOI ) ;
!       print*, 'infoi'
!       do i = 1,4
!         print*, INFOI(i)
!       end do
!       print*, 'Iter = ',infoi(2)

      CALL pressure_WASUPP_CG
            
      open(unit=20,file='dgfd_ver2.txt',status='replace',action='write')
      open(unit=30,file='dgcgsoln.txt',status='replace',action='write')
      do l = 1,ne
        write(30,'(4i,30f32.16)') l, nm(l,1), nm(l,2), nm(l,3), 
     &                 x(nm(l,1)), x(nm(l,2)), x(nm(l,3)),
     &                 y(nm(l,1)), y(nm(l,2)), y(nm(l,3)),
     &                 eta2(nm(l,1)), eta2(nm(l,2)), eta2(nm(l,3)),
     &                 uu2(nm(l,1)), uu2(nm(l,2)), uu2(nm(l,3)),
     &                 vv2(nm(l,1)), vv2(nm(l,2)), vv2(nm(l,3)),
     &                 dp(nm(l,1)), dp(nm(l,2)), dp(nm(l,3)),
     &                 pp_p1(nm(l,1)), pp_p1(nm(l,2)), pp_p1(nm(l,3)),
     &                 pp_p2(nm(l,1)), pp_p2(nm(l,2)), pp_p2(nm(l,3)),
     &                 pp_pd(nm(l,1)), pp_pd(nm(l,2)), pp_pd(nm(l,3)),
     &                 pp_pb(nm(l,1)), pp_pb(nm(l,2)), pp_pb(nm(l,3))
      end do
      inode = 1 
      DO i = 1, NumRows 
        DO j = 1, NumCols 
           eta_pp(i,j) = ETA2(inode) 
           inode = inode + 1
           write(20,'(2i,f)') i, j, eta_pp(i,j)
        END DO !i 
      END DO !j                   
      close(20)
      close(30)
      
!       CALL pressure_WASUPP
      
      deallocate(dgcgeta,dgcgloc)
!       call exit
      
      END SUBROUTINE from_DG_to_Pressure      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      END MODULE

