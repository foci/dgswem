      MODULE messenger2

#ifdef CMPI      
      USE mpi
#endif      

!$    USE omp_lib 

      USE globals, ONLY: rp,nrblk
      USE quit, ONLY: abort

      IMPLICIT NONE
      
      INTEGER :: ierr
      INTEGER :: nproc
      INTEGER :: nthreads
      INTEGER :: myrank
      INTEGER :: myid
      INTEGER :: world_group
      INTEGER :: comp_group
      INTEGER :: comp_comm
      INTEGER :: comm_dist_graph
      
      INTEGER :: lname
      CHARACTER(6) :: dirname
       
      INTEGER nproc_sr
      INTEGER, ALLOCATABLE, DIMENSION(:) :: proc_sr
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ned_sr
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: pe_sr      
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: el_sr
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: led_sr
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: nqpte_sr
      REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: nx_sr
      REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: ny_sr
      REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: hb_sr
      REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: detJe_sr
      
      INTEGER, ALLOCATABLE, DIMENSION(:) :: proc_group
      
      INTEGER, ALLOCATABLE, DIMENSION(:) :: lel2gel
      
      INTEGER :: nred
      INTEGER, ALLOCATABLE, DIMENSION(:) :: redn
      
      INTEGER :: mnired
      INTEGER :: mnelred
      
      TYPE :: edge_ptr_array
        REAL(rp), POINTER :: ptr
      END TYPE edge_ptr_array
      
      TYPE(edge_ptr_array), ALLOCATABLE, DIMENSION(:,:) :: Hri !,Hre
      TYPE(edge_ptr_array), ALLOCATABLE, DIMENSION(:,:) :: Zri,Zre      
      TYPE(edge_ptr_array), ALLOCATABLE, DIMENSION(:,:) :: Qxri,Qxre
      TYPE(edge_ptr_array), ALLOCATABLE, DIMENSION(:,:) :: Qyri,Qyre   
      TYPE(edge_ptr_array), ALLOCATABLE, DIMENSION(:,:) :: xmri,ymri,xymri    
      TYPE(edge_ptr_array), ALLOCATABLE, DIMENSION(:,:) :: Exxri,Eyyri,Exyri,Eyxri
      TYPE(edge_ptr_array), ALLOCATABLE, DIMENSION(:,:) :: Exxre,Eyyre,Exyre,Eyxre      

      REAL(rp), ALLOCATABLE, DIMENSION(:) :: xmre,ymre,xymre
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: Hre
      
      REAL(rp), ALLOCATABLE, TARGET, DIMENSION(:,:) :: sol_recv
      REAL(rp), ALLOCATABLE, TARGET, DIMENSION(:,:) :: sol_send      
      
      REAL(rp), ALLOCATABLE, TARGET, DIMENSION(:,:) :: sol_recv_ldg
      REAL(rp), ALLOCATABLE, TARGET, DIMENSION(:,:) :: sol_send_ldg        

      TYPE(edge_ptr_array), ALLOCATABLE, DIMENSION(:,:) :: send_ptr    
      TYPE(edge_ptr_array), ALLOCATABLE, DIMENSION(:,:) :: send_ptr_ldg         
      
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: hbr      
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: rnx,rny
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: rcfac      
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: detJe_recv      
      
      INTEGER, ALLOCATABLE, DIMENSION(:) :: solreq_recv
      INTEGER, ALLOCATABLE, DIMENSION(:) :: solreq_send
      INTEGER, ALLOCATABLE, DIMENSION(:) :: solreq  
      INTEGER, ALLOCATABLE, DIMENSION(:) :: solreq_recv_ldg
      INTEGER, ALLOCATABLE, DIMENSION(:) :: solreq_send_ldg      
      INTEGER, ALLOCATABLE, DIMENSION(:) :: solreq_ldg        
      
      INTEGER, ALLOCATABLE, DIMENSION(:) :: win    
      
      INTEGER :: match_edge
      
#ifdef ALIGN64
!DIR$ ATTRIBUTES ALIGN:64 :: Zri,Zre,Hri,Hre,Qxri,Qxre,Qyri,Qyre
!DIR$ ATTRIBUTES ALIGN:64 :: xmri,xmre,ymri,ymre,xymri,xymre
!DIR$ ATTRIBUTES ALIGN:64 :: Zfri,Hfri,Qxfri,Qyfri
!DIR$ ATTRIBUTES ALIGN:64 :: rnx,rny,len_area_recv
!DIR$ ATTRIBUTES ALIGN:64 :: sol_send,sol_recv
!DIR$ ATTRIBUTES ALIGN:64 :: send_ptr
#endif

#ifdef ALIGN32
!DIR$ ATTRIBUTES ALIGN:32 :: Zri,Zre,Hri,Hre,Qxri,Qxre,Qyri,Qyre
!DIR$ ATTRIBUTES ALIGN:32 :: xmri,xmre,ymri,ymre,xymri,xymre
!DIR$ ATTRIBUTES ALIGN:32 :: Zfri,Hfri,Qxfri,Qyfri
!DIR$ ATTRIBUTES ALIGN:32 :: rnx,rny,len_area_recv
!DIR$ ATTRIBUTES ALIGN:32 :: sol_send,sol_recv
!DIR$ ATTRIBUTES ALIGN:32 :: send_ptr
#endif 

#ifdef ALIGN16
!DIR$ ATTRIBUTES ALIGN:16 :: Zri,Zre,Hri,Hre,Qxri,Qxre,Qyri,Qyre
!DIR$ ATTRIBUTES ALIGN:16 :: xmri,xmre,ymri,ymre,xymri,xymre
!DIR$ ATTRIBUTES ALIGN:16 :: Zfri,Hfri,Qxfri,Qyfri
!DIR$ ATTRIBUTES ALIGN:16 :: rnx,rny,len_area_recv
!DIR$ ATTRIBUTES ALIGN:16 :: sol_send,sol_recv
!DIR$ ATTRIBUTES ALIGN:16 :: send_ptr
#endif

      CONTAINS
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE directory_name()
       
      IMPLICIT NONE

#ifdef CMPI       
      dirname = "PE0000"
      lname = 6
       
      WRITE(dirname(3:lname),"(I4.4)") myrank
       
!      PRINT*, "dirname = ", dirname
#else
      dirname = "."
      lname = 1
#endif       
       
      RETURN
      END SUBROUTINE directory_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
      SUBROUTINE message_init()
      
      USE version, ONLY: version_information
         
      IMPLICIT NONE
      
      CHARACTER(10) :: date,time
      
          
       
#ifdef CMPI       
      CALL MPI_INIT(ierr)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
!      CALL MPI_COMM_GROUP(MPI_COMM_WORLD,world_group,ierr)


      IF(myrank == 0) THEN       
        PRINT*, "Processor rank is ", myrank, "/", nproc
      ENDIF
       
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)       
      
      nrblk = 1
#else
      myrank = 0
      nrblk = 1      
      nproc = 1
#endif

#ifdef openmp      
      nthreads = omp_get_max_threads()        
      
      PRINT*, "Thread numbers :"
!$OMP  parallel private(myid)
      myid = omp_get_thread_num()
      PRINT*, myid        
!$OMP end parallel

      nrblk = nthreads  
#endif

      CALL directory_name()
      
      
      CALL DATE_AND_TIME(date,time)     
      IF (myrank == 0) THEN
        PRINT*, ' '
        PRINT*, 'Started: ',date,time
        PRINT*, ' '        
        
        CALL version_information(unit=6)
      ENDIF           
       
      RETURN           
      END SUBROUTINE message_init
       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      


      SUBROUTINE read_message_files()
      
      USE globals, ONLY: rp,ne,mndof,nout_sol,mnqpte
      
      IMPLICIT NONE
      
   
      INTEGER :: pe,ed,el,j,pt
      LOGICAL :: file_exists
      INTEGER :: npe,tne,mne,ne2,ndof2
      REAL(rp) :: nout2
      INTEGER :: reorder,info

#ifdef CMPI         
      IF (myrank == 0) THEN 
        PRINT "(A)", "Reading message passing files..."
      ENDIF
      
      INQUIRE(FILE=dirname(1:lname)//'/'//'fort.18', EXIST = file_exists)
      IF(file_exists == .FALSE.) THEN
        PRINT*, "fort.18 file does not exist for process", myrank
        CALL abort()
      ENDIF     
      
      OPEN(UNIT=18,FILE=dirname(1:lname)//'/'//'fort.18')
      
      READ(18,182) nproc_sr
      
      ALLOCATE(proc_sr(nproc_sr),ned_sr(nproc_sr))
      ALLOCATE(el_sr(ne,nproc_sr),led_sr(ne,nproc_sr))
      ALLOCATE(nx_sr(mnqpte,nred,nproc_sr),ny_sr(mnqpte,nred,nproc_sr))
      ALLOCATE(hb_sr(mnqpte,nred,nproc_sr))
      ALLOCATE(detJe_sr(mnqpte,nred,nproc_sr))
      ALLOCATE(nqpte_sr(nred,nproc_sr))
      
      DO pe = 1,nproc_sr
!        j = mod(pe-1+myrank,nproc_sr)+1 ! DG ADCIRC method
!        READ(18,181) proc_sr(j),ned_sr(j) 
!        READ(18,180) (el_sr(ed,j), ed = 1,ned_sr(j))
!        READ(18,180) (led_sr(ed,j), ed = 1,ned_sr(j))
        

        READ(18,181) proc_sr(pe),ned_sr(pe) !DG ADCIRC may try to do some straggering here, not really sure about what they are doing 
        READ(18,180) (el_sr(ed,pe), ed = 1,ned_sr(pe))
        READ(18,180) (led_sr(ed,pe), ed = 1,ned_sr(pe))
        DO ed = 1,ned_sr(pe)
          READ(18,182) nqpte_sr(ed,pe)
          DO pt = 1,nqpte_sr(ed,pe)
            READ(18,183) nx_sr(pt,ed,pe),ny_sr(pt,ed,pe),hb_sr(pt,ed,pe),detJe_sr(pt,ed,pe)
          ENDDO
        ENDDO
      ENDDO
 
      CLOSE(18)

!      ALLOCATE(proc_group(nproc_sr+1))
!      DO pe = 1,nproc_sr
!        proc_group(pe) = proc_sr(pe)
!      ENDDO      
!      proc_group(nproc_sr+1) = myrank

!      CALL MPI_GROUP_INCL(world_group,nproc_sr+1,proc_group,comp_group,ierr)
!      CALL MPI_COMM_CREATE(MPI_COMM_WORLD,comp_group,comp_comm,ierr)

      reorder = 1
      info = MPI_INFO_NULL

      CALL MPI_Dist_graph_create_adjacent(MPI_COMM_WORLD,nproc_sr,proc_sr,ned_sr, &
                                                         nproc_sr,proc_sr,ned_sr, &
                                                         info,reorder,comm_dist_graph,ierr)      
      
      INQUIRE(FILE=dirname(1:lname)//'/'//'fort.80', EXIST = file_exists)
      IF(file_exists == .FALSE.) THEN
        PRINT*, "fort.80 file does not exist for process", myrank
        CALL abort()
      ENDIF   
      
      OPEN(UNIT=80,FILE=dirname(1:lname)//'/'//'fort.80')  
      
      READ(80,*) npe
      READ(80,*) tne
      READ(80,*) mne
      READ(80,*) ne2
      READ(80,*) ndof2
      READ(80,*) nout2
      
      ! Check the post processing information against run input 
      IF(npe /= nproc) THEN
        PRINT*, "number of processors does not agree with files"
        CALL abort()
      ENDIF
      
      IF(ne /= ne2) THEN
        PRINT*, "number of elements does not agree with files"
        CALL abort()
      ENDIF
      
      IF(ndof2 /= mndof) THEN
        PRINT*, "number of degrees of freedom does not agree with files"
        CALL abort()
      ENDIF
      
!       IF(nout2 /= nout_sol) THEN
!         PRINT*, "number of output lines does not agree with files"
!         CALL abort()
!       ENDIF
      
      ALLOCATE(lel2gel(ne))
      
      DO el = 1,ne
        READ(80,*) j, lel2gel(el)
      ENDDO
      
      CLOSE(80)
      
      
 180  FORMAT(8X,9I8)
 181  FORMAT(8X,2I8) 
 182  FORMAT(8X,I8) 
 183  FORMAT(8X,4(D24.17,1X))
#else

      ALLOCATE(lel2gel(ne))
      
      DO el = 1,ne
        lel2gel(el) = el
      ENDDO
 
#endif   
      RETURN      
      END SUBROUTINE read_message_files
    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      SUBROUTINE message_setup()
      
      USE globals, ONLY: ned,mnqpte,nqpte,detJe,nx_pt,ny_pt,hbqpte, &
                         Zqpt,Hqpt,Qxqpt,Qyqpt, &
                         xmom,ymom,xymom, &
                         Exxqpt,Eyyqpt,Exyqpt,Eyxqpt, &                         
                         ged2el,ged2led,gel2ael, &
                         Spe,cfac,xy,ect
      
      IMPLICIT NONE
      
#ifdef CMPI
      INTEGER :: pe,ed,pt
      INTEGER :: el,led,gpt,ged,el_in
      INTEGER :: i,j,edcnt
      INTEGER :: gp_in,gp_ex
      INTEGER :: mnedsr,nqpt
      INTEGER :: match
      
      IF (myrank == 0) THEN 
        PRINT "(A)", "Seting up message passing arrays..."      
      ENDIF
        
      mnedsr = MAXVAL(ned_sr)
      
      ALLOCATE(send_ptr(3*mnedsr*mnqpte,nproc_sr))
      ALLOCATE(sol_send(3*mnedsr*mnqpte,nproc_sr))
      
      ALLOCATE(send_ptr_ldg(4*mnedsr*mnqpte,nproc_sr))
      ALLOCATE(sol_send_ldg(4*mnedsr*mnqpte,nproc_sr))      
      
      DO pe = 1,nproc_sr
        i = 0
        DO ed = 1,ned_sr(pe)
          el = el_sr(ed,pe)
          el_in = gel2ael(el)
          led = led_sr(ed,pe)          
          DO pt = 1,nqpte(1)
            i = i+1
            
!             gpt = (led-1)*nqpte(1) + pt
            gpt = (led-1)*nqpte(1) + nqpte(1) - pt + 1            
            
!             send_ptr(i,pe)%ptr => Hqpt(el_in,gpt)
            send_ptr(i,pe)%ptr => Zqpt(el_in,gpt)
            send_ptr(ned_sr(pe)*nqpte(1) + i,pe)%ptr => Qxqpt(el_in,gpt)
            send_ptr(2*ned_sr(pe)*nqpte(1) + i,pe)%ptr => Qyqpt(el_in,gpt)
            
            send_ptr_ldg(i,pe)%ptr => Exxqpt(el_in,gpt)
            send_ptr_ldg(ned_sr(pe)*nqpte(1) + i,pe)%ptr => Eyyqpt(el_in,gpt)
            send_ptr_ldg(2*ned_sr(pe)*nqpte(1) + i,pe)%ptr => Exyqpt(el_in,gpt)
            send_ptr_ldg(3*ned_sr(pe)*nqpte(1) + i,pe)%ptr => Eyxqpt(el_in,gpt)
          ENDDO
        ENDDO
      ENDDO
      
      
      ALLOCATE(sol_recv(3*mnedsr*mnqpte,nproc_sr))
      ALLOCATE(sol_recv_ldg(4*mnedsr*mnqpte,nproc_sr))      
!       ALLOCATE(Hre(nred,mnqpte))
      ALLOCATE(Zre(nred,mnqpte),Qxre(nred,mnqpte),Qyre(nred,mnqpte))
      ALLOCATE(Zri(nred,mnqpte),Hri(nred,mnqpte),Qxri(nred,mnqpte),Qyri(nred,mnqpte))
      ALLOCATE(xmri(nred,mnqpte),ymri(nred,mnqpte),xymri(nred,mnqpte))
      ALLOCATE(Exxri(nred,mnqpte),Eyyri(nred,mnqpte),Exyri(nred,mnqpte),Eyxri(nred,mnqpte))
      ALLOCATE(Exxre(nred,mnqpte),Eyyre(nred,mnqpte),Exyre(nred,mnqpte),Eyxre(nred,mnqpte))      
      ALLOCATE(xmre(nred),ymre(nred),xymre(nred))
      ALLOCATE(Hre(nred))      
      ALLOCATE(hbr(nred,mnqpte))
      
      match_edge = 1
      edcnt = 0     
      DO pe = 1,nproc_sr
        i = 0
        DO ed = 1,ned_sr(pe)
          edcnt = edcnt + 1
          el = el_sr(ed,pe)
          led = led_sr(ed,pe)                            
          
          DO pt = 1,nqpte(1)
            i = i+1
            
!             gp_ex = nqpte(1) - pt + 1
            gp_ex = pt
            
!             Hre(edcnt,gp_ex)%ptr => sol_recv(i,pe)
            Zre(edcnt,gp_ex)%ptr => sol_recv(i,pe)            
            Qxre(edcnt,gp_ex)%ptr => sol_recv(ned_sr(pe)*nqpte(1) + i,pe)
            Qyre(edcnt,gp_ex)%ptr => sol_recv(2*ned_sr(pe)*nqpte(1) + i,pe)
            
            Exxre(edcnt,gp_ex)%ptr => sol_recv_ldg(i,pe)
            Eyyre(edcnt,gp_ex)%ptr => sol_recv_ldg(ned_sr(pe)*nqpte(1) + i,pe) 
            Exyre(edcnt,gp_ex)%ptr => sol_recv_ldg(2*ned_sr(pe)*nqpte(1) + i,pe)
            Eyxre(edcnt,gp_ex)%ptr => sol_recv_ldg(3*ned_sr(pe)*nqpte(1) + i,pe)
            
            gp_in = (led-1)*nqpte(1) + pt
            el_in = gel2ael(el)            
            
            Hri(edcnt,pt)%ptr  => Hqpt(el_in,gp_in)
            Zri(edcnt,pt)%ptr  => Zqpt(el_in,gp_in)
            Qxri(edcnt,pt)%ptr => Qxqpt(el_in,gp_in)
            Qyri(edcnt,pt)%ptr => Qyqpt(el_in,gp_in)
            
            xmri(edcnt,pt)%ptr  => xmom(el_in,gp_in)
            ymri(edcnt,pt)%ptr  => ymom(el_in,gp_in)
            xymri(edcnt,pt)%ptr => xymom(el_in,gp_in)

            Exxri(edcnt,pt)%ptr => Exxqpt(el_in,gp_in)
            Eyyri(edcnt,pt)%ptr => Eyyqpt(el_in,gp_in)
            Exyri(edcnt,pt)%ptr => Exyqpt(el_in,gp_in)
            Eyxri(edcnt,pt)%ptr => Eyxqpt(el_in,gp_in)            
            
          ENDDO
        ENDDO
      ENDDO
      
      ALLOCATE(rnx(nred,mnqpte),rny(nred,mnqpte),rcfac(nred,mnqpte),detJe_recv(nred,mnqpte))
      
      edcnt = 0
      DO pe = 1,nproc_sr
        DO ed = 1,ned_sr(pe)
          edcnt = edcnt+1
          
          el = el_sr(ed,pe)
          led = led_sr(ed,pe)
          nqpt = nqpte_sr(ed,pe)
          el_in = gel2ael(el)
          
          match = 0
          DO i = 1,nred
            ged = redn(i)
            IF(el == ged2el(1,ged) .and. led == ged2led(1,ged)) THEN
            
              DO pt = 1,nqpte(1)
              
                gp_in = (led-1)*nqpte(1) + pt
                
!                 rnx(edcnt,pt) = nx_pt(ged,pt)*Spe(ged,pt)
!                 rny(edcnt,pt) = ny_pt(ged,pt)
                rnx(edcnt,pt) = nx_sr(pt,ed,pe)*Spe(ged,pt)
                rny(edcnt,pt) = ny_sr(pt,ed,pe)
                                       
                hbr(edcnt,pt) = hb_sr(pt,ed,pe) ! assumes continuous bathymetry     
                hbqpte(el_in,gp_in) = hb_sr(pt,ed,pe)

!                 rcfac(edcnt,pt) = cfac(ed,pt)
                rcfac(edcnt,pt) = ny_sr(pt,ed,pe)**2+(nx_sr(pt,ed,pe)*Spe(ged,pt))**2
              
                detJe_recv(edcnt,pt) = detJe_sr(pt,ed,pe)
              ENDDO
              
              match = 1
              
            ENDIF
          ENDDO
          
          IF(match == 0) THEN
            PRINT*, "global edge not found for recieve edge"
            CALL abort()
          ENDIF
          
          
        ENDDO
      ENDDO
#endif       
      
      RETURN
      END SUBROUTINE message_setup
     

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      SUBROUTINE communication_setup()
      
      USE globals, ONLY: nqpte,rp
      
      IMPLICIT NONE
#ifdef CMPI      
      
      INTEGER :: pe,tag
      
      IF (myrank == 0) THEN 
        PRINT "(A)", "Initializing message passing..."       
      ENDIF
      
      ALLOCATE(solreq_recv(nproc_sr))
      ALLOCATE(solreq_send(nproc_sr))
      ALLOCATE(solreq(2*nproc_sr))
      ALLOCATE(solreq_recv_ldg(nproc_sr))
      ALLOCATE(solreq_send_ldg(nproc_sr))
      ALLOCATE(solreq_ldg(2*nproc_sr))
      
      tag = 1
      
      DO pe = 1,nproc_sr
        CALL MPI_RECV_INIT(sol_recv(1,pe),3*ned_sr(pe)*nqpte(1),MPI_DOUBLE_PRECISION,proc_sr(pe), &
                           tag,comm_dist_graph,solreq_recv(pe),ierr)

!        CALL MPI_RECV_INIT(sol_recv(1,pe),3*ned_sr(pe)*nqpte(1),MPI_DOUBLE_PRECISION,proc_sr(pe), &
!                           tag,MPI_COMM_WORLD,solreq_recv(pe),ierr)

!        CALL MPI_RECV_INIT(sol_recv(1,pe),3*ned_sr(pe)*nqpte(1),MPI_DOUBLE_PRECISION,pe, &
!                           tag,comp_comm,solreq_recv(pe),ierr)
      ENDDO
      
      DO pe = 1,nproc_sr
        CALL MPI_SEND_INIT(sol_send(1,pe),3*ned_sr(pe)*nqpte(1),MPI_DOUBLE_PRECISION,proc_sr(pe), &
                           tag,comm_dist_graph,solreq_send(pe),ierr)

!        CALL MPI_SEND_INIT(sol_send(1,pe),3*ned_sr(pe)*nqpte(1),MPI_DOUBLE_PRECISION,proc_sr(pe), &
!                           tag,MPI_COMM_WORLD,solreq_send(pe),ierr)                           

!        CALL MPI_SEND_INIT(sol_send(1,pe),3*ned_sr(pe)*nqpte(1),MPI_DOUBLE_PRECISION,pe, &
!                           tag,comp_comm,solreq_send(pe),ierr)                                 
      ENDDO

      DO pe = 1,nproc_sr
        solreq(pe) = solreq_recv(pe)
        solreq(pe+nproc_sr) = solreq_send(pe)
      ENDDO
      
!       ALLOCATE(win(nproc_sr))
!       
!       DO pe = 1,nproc_sr
!         CALL MPI_Win_create(sol_send(1,pe),3*ned_sr(pe)*nqpte(1)*rp,rp, &
!                             MPI_INFO_NULL,MPI_COMM_WORLD,win(pe))
!       ENDDO
      
      
      DO pe = 1,nproc_sr
        CALL MPI_RECV_INIT(sol_recv_ldg(1,pe),3*ned_sr(pe)*nqpte(1),MPI_DOUBLE_PRECISION,proc_sr(pe), &
                           tag,comm_dist_graph,solreq_recv_ldg(pe),ierr)
      ENDDO
      
      DO pe = 1,nproc_sr
        CALL MPI_SEND_INIT(sol_send_ldg(1,pe),3*ned_sr(pe)*nqpte(1),MPI_DOUBLE_PRECISION,proc_sr(pe), &
                           tag,comm_dist_graph,solreq_send_ldg(pe),ierr)                          
      ENDDO
      
      DO pe = 1,nproc_sr
        solreq_ldg(pe) = solreq_recv_ldg(pe)
        solreq_ldg(pe+nproc_sr) = solreq_send_ldg(pe)
      ENDDO
      
#endif        
      RETURN
      END SUBROUTINE communication_setup
    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef CMPI
      SUBROUTINE message_recieve()
      
      IMPLICIT NONE

      CALL MPI_STARTALL(nproc_sr,solreq_recv,ierr)
      
      RETURN
      END SUBROUTINE message_recieve
#endif      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef CMPI
      SUBROUTINE message_send()
      
      USE globals, ONLY: nqpte
      
      IMPLICIT NONE
      
      INTEGER :: pe,i
      
      DO pe = 1,nproc_sr
        DO i = 1,3*ned_sr(pe)*nqpte(1)
        
          sol_send(i,pe) = send_ptr(i,pe)%ptr
          
        ENDDO
      ENDDO
      
!      CALL MPI_STARTALL(nproc_sr,solreq_send,ierr)
      CALL MPI_STARTALL(2*nproc_sr,solreq,ierr)
      
      RETURN
      END SUBROUTINE message_send
#endif      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef CMPI
      SUBROUTINE message_send_ldg()
      
      USE globals, ONLY: nqpte
      
      IMPLICIT NONE
      
      INTEGER :: pe,i
      
      DO pe = 1,nproc_sr
        DO i = 1,4*ned_sr(pe)*nqpte(1)
        
          sol_send_ldg(i,pe) = send_ptr_ldg(i,pe)%ptr
          
        ENDDO
      ENDDO
      
      CALL MPI_STARTALL(2*nproc_sr,solreq_ldg,ierr)
      
      RETURN
      END SUBROUTINE message_send_ldg
#endif     

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       SUBROUTINE end_time(t_start,nproc)

       IMPLICIT NONE
       
       INTEGER :: i
       INTEGER :: ierr
       INTEGER :: nproc
       REAL(rp) :: t_start,t_end       
       REAL(rp) :: t_max,t_min,t_avg
       REAL(rp) :: cpu_times(nproc)
       

#ifdef openmp      
      t_end = omp_get_wtime()
#else
      CALL CPU_TIME(t_end)     
#endif      
       
#ifdef CMPI       
      CALL MPI_GATHER(t_end-t_start,1,MPI_DOUBLE_PRECISION,cpu_times,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      
      IF (myrank == 0) THEN
        t_avg = 0d0
        DO i = 1,nproc
          t_avg = t_avg + cpu_times(i)
        ENDDO
        t_avg = t_avg/nproc
      
        t_max = MAXVAL(cpu_times)
        t_min = MINVAL(cpu_times)
      
        PRINT*, ' '      
        PRINT("(A,F25.5,A)"), "Average CPU time = ",t_avg," seconds"
        PRINT("(A,F25.5,A)"), "Minimum CPU time = ",t_min," seconds"
        PRINT("(A,F25.5,A)"), "Maximum CPU time = ",t_max," seconds" 
      ENDIF   
      
#else      
      PRINT*, ' '      
      PRINT("(A,F25.5,A)"), "CPU time = ",t_end-t_start," seconds"      
#endif




       RETURN
       END SUBROUTINE end_time

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       

      END MODULE messenger2
