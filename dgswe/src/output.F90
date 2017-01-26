      MODULE output
      
#ifdef NETCDF
      USE netcdf       
#endif      

      USE globals, ONLY: rp
      USE messenger2, ONLY: myrank,dirname,lname
      USE read_dginp, ONLY: out_direc,grid_file
      USE read_write_output, ONLY: file_init,write_solution_snap,time_snaps
       
      IMPLICIT NONE
      INTEGER :: ncid
      INTEGER :: dimID_time,dimID_mndof,dimID_el
      INTEGER :: dimID_ldof,dimID_tdof,dimID_mnnds,dimID_neltype
      INTEGER :: zid,qxid,qyid,tid,hbid   
      INTEGER :: ml2id,mmlid,m2nid
      INTEGER, DIMENSION(3) :: sol_dim
      INTEGER, DIMENSION(2) :: hb_dim        
      INTEGER, DIMENSION(2) :: ml2_dim          
      INTEGER, DIMENSION(2) :: mml_dim
      INTEGER, DIMENSION(3) :: m2n_dim
!       INTEGER, DIMENSION(0) :: scal_dim

      INTEGER :: nsnap
       
      CONTAINS
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      

      SUBROUTINE nc_setup()

      USE globals, ONLY: ne,ndof,mndof,mnnds,ml2,mml,nel_type

      IMPLICIT NONE

#ifdef NETCDF
      CALL check(NF90_CREATE(dirname(1:lname)//'/'//'solution.nc',NF90_CLOBBER,ncid))

      CALL check(NF90_DEF_DIM(ncid,'time',NF90_UNLIMITED,dimID_time))
      CALL check(NF90_DEF_DIM(ncid,'tri_dof',ndof(1),dimID_tdof))
      CALL check(NF90_DEF_DIM(ncid,'mndof',mndof,dimID_mndof))   
      CALL check(NF90_DEF_DIM(ncid,'mnnds',mnnds,dimID_mnnds))            
      CALL check(NF90_DEF_DIM(ncid,'el',ne,dimID_el))
      CALL check(NF90_DEF_DIM(ncid,'linear_dof',3,dimID_ldof))
      CALL check(NF90_DEF_DIM(ncid,'nel_type',nel_type,dimID_neltype))          
          
      sol_dim = (/ dimID_el, dimID_mndof, dimID_time /)
        
      CALL check(NF90_DEF_VAR(ncid,'Z',NF90_DOUBLE,sol_dim,zid))
      CALL check(NF90_DEF_VAR(ncid,'Qx',NF90_DOUBLE,sol_dim,qxid))
      CALL check(NF90_DEF_VAR(ncid,'Qy',NF90_DOUBLE,sol_dim,qyid))
      CALL check(NF90_DEF_VAR(ncid,'t',NF90_DOUBLE,dimID_time,tid))
          
      hb_dim = (/ dimID_el, dimID_mnnds /)
      CALL check(NF90_DEF_VAR(ncid,'hb',NF90_DOUBLE,hb_dim,hbid))          
          
      m2n_dim = (/ dimID_mnnds, dimID_mndof, dimID_neltype/)
      CALL check(NF90_DEF_VAR(ncid,'m2n',NF90_DOUBLE,m2n_dim,m2nid))
          
      ml2_dim = (/ dimID_ldof, dimID_tdof /)        
      CALL check(NF90_DEF_VAR(ncid,'ml2',NF90_DOUBLE,ml2_dim,ml2id))
          
      mml_dim = (/ dimID_ldof, dimID_ldof /)          
      CALL check(NF90_DEF_VAR(ncid,'mml',NF90_DOUBLE,mml_dim,mmlid))
          
      CALL check(NF90_ENDDEF(ncid))
          
          
      CALL check(NF90_PUT_VAR(ncid,ml2id,ml2))
      CALL check(NF90_PUT_VAR(ncid,mmlid,mml))       
          
      nsnap = 1
#endif      
          
          
      END SUBROUTINE nc_setup
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
        
      SUBROUTINE output_solution(init)
      
      USE globals, ONLY: t,mndof,ne,mnnds,tskp_sol,nout_sol, &
                         Hwrite,Zwrite,Qxwrite,Qywrite,hbm, &
                         Zout,Qxout,Qyout, &
                         Zsol_unit,Qxsol_unit,Qysol_unit
      USE read_dginp, ONLY: tf,dt,sol_opt,sol_snap,hbp                         

      IMPLICIT NONE

      INTEGER :: dof,el
      LOGICAL :: init
      INTEGER :: nout_hb
      INTEGER :: hb_unit
      INTEGER :: ndof_hb
      INTEGER, DIMENSION(3) :: var_start,var_end
      REAL(rp), DIMENSION(1) :: t_tmp
     

      IF (init) THEN
      
        IF (myrank == 0) THEN 
          PRINT "(A)", "Initializing solution output files..."         
        ENDIF
        
        ALLOCATE(Zout(ne,mndof),Qxout(ne,mndof),Qyout(ne,mndof))          
        
        ! Set number of timesteps between output
        CALL time_snaps(sol_opt,sol_snap,tf,dt,tskp_sol,nout_sol) 
        nout_hb = 1
      
        ! Initialize files and write initial condition                    
        CALL file_init(out_direc,"Z.sol",mndof,ne,nout_sol+1,Zsol_unit)
        CALL file_init(out_direc,"Qx.sol",mndof,ne,nout_sol+1,Qxsol_unit)
        CALL file_init(out_direc,"Qy.sol",mndof,ne,nout_sol+1,Qysol_unit) 
        
        ndof_hb = (hbp+1)**2
        CALL file_init(out_direc,"hb.sol",ndof_hb,ne,nout_hb,hb_unit)
        CALL write_solution_snap(hb_unit,ndof_hb,ne,"N",t,hbm)         
        CLOSE(hb_unit)
        
        ! Set up netcdf output files
        CALL nc_setup()       

      ELSE

        IF(myrank == 0) THEN
          PRINT("(A,e15.8)"), 't = ', t
        ENDIF
        
                  
      ENDIF

      
      DO dof = 1,mndof
        DO el = 1,ne
          Zout(el,dof)  = Zwrite(el,dof)%ptr
          Qxout(el,dof) = Qxwrite(el,dof)%ptr
          Qyout(el,dof) = Qywrite(el,dof)%ptr
        ENDDO
      ENDDO
      
      CALL write_solution_snap(Zsol_unit,mndof,ne,"T",t,Zout) 
      CALL write_solution_snap(Qxsol_unit,mndof,ne,"T",t,Qxout)      
      CALL write_solution_snap(Qysol_unit,mndof,ne,"T",t,Qyout)            
      
      
#ifdef NETCDF      

      
      
      
      var_start = (/ 1, 1, nsnap /)
      var_end = (/ ne, mndof, 1 /)
      
      t_tmp(1) = t
      CALL check(NF90_PUT_VAR(ncid,tid,t_tmp,(/nsnap/),(/1/)))
      CALL check(NF90_PUT_VAR(ncid,zid,Znc,var_start,var_end ))
      CALL check(NF90_PUT_VAR(ncid,qxid,Qxnc,var_start,var_end ))
      CALL check(NF90_PUT_VAR(ncid,qyid,Qync,var_start,var_end ))

      nsnap = nsnap + 1
#endif      
      

      RETURN
      END SUBROUTINE output_solution
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        

      SUBROUTINE output_stations(init)
      
      USE globals, ONLY: t,ndof,mndof,nsta,xysta,elsta,el_type,hbsta,phi_sta,tskp_sta,nout_sta, &
                         Hwrite,Zwrite,Qxwrite,Qywrite, &
                         Zsta,Qxsta,Qysta, &
                         Zsta_unit,Qxsta_unit,Qysta_unit
      USE read_dginp, ONLY: tf,dt,sta_opt,sta_snap                         

      IMPLICIT NONE

      INTEGER :: dof,sta
      INTEGER :: elin,et
      INTEGER :: ncol
      INTEGER :: nout_hb
      INTEGER :: hb_unit      
      LOGICAL :: init

      INTEGER, DIMENSION(3) :: var_start,var_end
      REAL(rp), DIMENSION(1) :: t_tmp
     

      IF (init) THEN      
      
        IF (sta_opt > 0) THEN
          CALL read_stations()         
          CALL find_stations()
        ELSE 
          nsta = 0
          RETURN
        ENDIF      
      
        IF (myrank == 0) THEN 
          PRINT "(A)", "Initializing station output files..."      
        ENDIF      
        
        ALLOCATE(Zsta(nsta,1),Qxsta(nsta,1),Qysta(nsta,1))
      
        ! Set number of timesteps between output
        CALL time_snaps(sta_opt,sta_snap,tf,dt,tskp_sta,nout_sta)    
        nout_hb = 1
        ncol = 1
        
        ! Initialize files and write initial condition                     
        CALL file_init(out_direc,"Z.sta",nsta,ncol,nout_sta+1,Zsta_unit)   
        CALL file_init(out_direc,"Qx.sta",nsta,ncol,nout_sta+1,Qxsta_unit)    
        CALL file_init(out_direc,"Qy.sta",nsta,ncol,nout_sta+1,Qysta_unit)  
        
        CALL file_init(out_direc,"hb.sta",nsta,ncol,nout_hb,hb_unit)         
        CALL write_solution_snap(hb_unit,nsta,ncol,"N",t,hbsta)         
        CLOSE(hb_unit)                
        
        ! Set up netcdf output files
!         CALL nc_setup()    

      ENDIF


                    
      Zsta = 0d0
      Qxsta = 0d0
      Qysta = 0d0   
                
      DO sta = 1,nsta
               
        elin = elsta(sta)                 
        et = el_type(elin)                 

        DO dof = 1,ndof(et)
          Zsta(sta,1)  = Zsta(sta,1)  + Zwrite(elin,dof)%ptr*phi_sta(dof,sta)
          Qxsta(sta,1) = Qxsta(sta,1) + Qxwrite(elin,dof)%ptr*phi_sta(dof,sta)            
          Qysta(sta,1) = Qysta(sta,1) + Qywrite(elin,dof)%ptr*phi_sta(dof,sta)
        ENDDO                
      ENDDO  
      
      
      
      CALL write_solution_snap(Zsta_unit,nsta,1,"N",t,Zsta) 
      CALL write_solution_snap(Qxsta_unit,nsta,1,"N",t,Qxsta)    
      CALL write_solution_snap(Qysta_unit,nsta,1,"N",t,Qysta)          
      
#ifdef NETCDF      
!       
!       var_start = (/ 1, 1, nsnap /)
!       var_end = (/ ne, mndof, 1 /)
!       
!       t_tmp(1) = t
!       CALL check(NF90_PUT_VAR(ncid,tid,t_tmp,(/nsnap/),(/1/)))
!       CALL check(NF90_PUT_VAR(ncid,zid,Znc,var_start,var_end ))
!       CALL check(NF90_PUT_VAR(ncid,qxid,Qxnc,var_start,var_end ))
!       CALL check(NF90_PUT_VAR(ncid,qyid,Qync,var_start,var_end ))
! 
!       nsnap = nsnap + 1
#endif      
      

      RETURN
      END SUBROUTINE output_stations
      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    

      SUBROUTINE close_output()     
      
      USE globals, ONLY: Zsol_unit,Qxsol_unit,Qysol_unit, &
                         Zsta_unit,Qxsta_unit,Qysta_unit
      
      IMPLICIT NONE      
      
      INTEGER :: i
      

      CLOSE(Zsol_unit)
      CLOSE(Qxsol_unit)
      CLOSE(Qysol_unit)
      
      CLOSE(Zsta_unit)
      CLOSE(Qxsta_unit)
      CLOSE(Qysta_unit)

#ifdef NETCDF      
      CALL check(NF90_CLOSE(ncid))
#endif      
      
      RETURN
      END SUBROUTINE close_output       
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
                
      SUBROUTINE check(status)

      IMPLICIT NONE
      INTEGER :: status

#ifdef NETCDF      
      IF(status /= NF90_NOERR) THEN
        PRINT("(A,A)"), "fatal error from ", TRIM(NF90_STRERROR(status))  
      ENDIF
#endif      
               
      END SUBROUTINE check

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      
      END MODULE output

      
      
      
      

 