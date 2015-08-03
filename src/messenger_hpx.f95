#ifdef HPX
       subroutine HPX_GET_ELEMS(dg_here,neighbor,volume,sendbuf)     
       
       use dg
       use sizes
       
       implicit none       

       type (dg_type) :: dg_here
       
       logical :: neighbor_found
       integer :: i,index
       integer :: el,dof
       integer :: ncount
       integer :: neighbor
       integer :: volume
       
!       real(sz) :: sendbuf(volume)
       real(sz) :: sendbuf(MAX_BUFFER_SIZE)   
       
       
       
       neighbor_found = .false.
       DO i=1,dg_here%NEIGHPROC_S

         IF (dg_here%IPROC_S(i) == neighbor) THEN
           index = i
           neighbor_found = .true.
           EXIT
         ENDIF
      
       ENDDO     
       
       IF (neighbor_found) THEN
       
         PRINT*, " PROC: ", dg_here%IDPROC, " SENDING ", dg_here%NELEMSEND(index) ," ELEMENTS TO: ", dg_here%IPROC_S(index)       
         PRINT 180, (dg_here%ISENDLOC(el,index), el = 1,dg_here%NELEMSEND(index))
         PRINT*, "" 
      
         ncount = 0
         DO el=1,dg_here%NELEMSEND(index)
            DO dof=1,dg_here%DOFH
              ncount = ncount+1
              sendbuf(ncount)=dg_here%ZE(dof,dg_here%ISENDLOC(el,index),dg_here%IRK+1)
!              IF (abs(sendbuf(ncount) ) > 1d-14) THEN              
!               PRINT 181, dg_here%ISENDLOC(el,index),sendbuf(ncount)
!              ENDIF
            ENDDO
         ENDDO
       
         DO el=1,dg_here%NELEMSEND(index)
           DO dof=1,dg_here%DOFH
             ncount = ncount+1
             sendbuf(ncount)=dg_here%QX(dof,dg_here%ISENDLOC(el,index),dg_here%IRK+1)
             IF (abs(sendbuf(ncount) ) > 1d-14) THEN
               PRINT 181, dg_here%ISENDLOC(el,index),sendbuf(ncount)
             ENDIF
           ENDDO
         ENDDO
          
         DO el=1,dg_here%NELEMSEND(index)
           DO dof=1,dg_here%DOFH
             ncount = ncount+1
             sendbuf(ncount)=dg_here%QY(dof,dg_here%ISENDLOC(el,index),dg_here%IRK+1)
!              IF (abs(sendbuf(ncount) ) > 1d-14) THEN                      
!                PRINT 181, dg_here%ISENDLOC(el,index),sendbuf(ncount)
!              ENDIF
           ENDDO
         ENDDO
         
         PRINT*, "----------------------------------------------------------------------------------"                  
       
       ELSE 
       
         PRINT*, "FORTRAN ERROR: neighbor not found"
         STOP             
    
       ENDIF
       
 180  FORMAT(8X,9I8)    
 181  FORMAT(I5,2X,ES24.17)
       
       end subroutine HPX_GET_ELEMS
       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!           
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       

       subroutine HPX_PUT_ELEMS(dg_here,neighbor,volume,recvbuf)
       
       use dg
       use sizes
       
       implicit none
       
       type (dg_type) :: dg_here
       
       logical :: neighbor_found
       integer :: i,index       
       integer :: el,dof
       integer :: ncount
       integer :: neighbor
       integer :: volume
       
!       real(sz) :: recvbuf(volume)
       real(sz) :: recvbuf(MAX_BUFFER_SIZE)
        
       
       neighbor_found = .false.
       DO i=1,dg_here%NEIGHPROC_R

         IF (dg_here%IPROC_R(i) == neighbor) THEN
           index = i
           neighbor_found = .true.
           EXIT
         ENDIF
      
       ENDDO           
             
       IF (neighbor_found) THEN  
       
         PRINT*, " PROC: ", dg_here%IDPROC, " RECEIVING ",dg_here%NELEMRECV(index), " ELEMENTS FROM: ", dg_here%IPROC_R(index)
         PRINT 180, (dg_here%IRECVLOC(el,index), el = 1,dg_here%NELEMRECV(index))
         PRINT*, ""
         
         ncount = 0
         DO el=1,dg_here%NELEMRECV(index)
           DO dof=1,dg_here%DOFH
             ncount = ncount+1
             dg_here%ZE(dof,dg_here%IRECVLOC(el,index),dg_here%IRK+1) = recvbuf(ncount)
!              IF (abs(recvbuf(ncount) ) > 1d-14) THEN             
!                PRINT 181, dg_here%IRECVLOC(el,index),recvbuf(ncount) 
!              ENDIF
           ENDDO
         ENDDO
       
         DO el=1,dg_here%NELEMRECV(index)
           DO dof=1,dg_here%DOFH
             ncount = ncount+1
             dg_here%QX(dof,dg_here%IRECVLOC(el,index),dg_here%IRK+1) = recvbuf(ncount)
             IF (abs(recvbuf(ncount) ) > 1d-14) THEN
               PRINT 181, dg_here%IRECVLOC(el,index),recvbuf(ncount)
             ENDIF
           ENDDO
         ENDDO
       
         DO el=1,dg_here%NELEMRECV(index)
           DO dof=1,dg_here%DOFH
             ncount = ncount+1
             dg_here%QY(dof,dg_here%IRECVLOC(el,index),dg_here%IRK+1) = recvbuf(ncount)
!              IF (abs(recvbuf(ncount) ) > 1d-14) THEN                
!                PRINT 181, dg_here%IRECVLOC(el,index),recvbuf(ncount)
!              ENDIF
           ENDDO
         ENDDO    
         
         PRINT*, "----------------------------------------------------------------------------------"         
       
       ELSE
       
         PRINT*, "FORTRAN ERROR: neighbor not found"
         STOP       
         
       ENDIF
         

 180  FORMAT(8X,9I8)     
 181  FORMAT(I5,2X,ES24.17) 
       
       end subroutine HPX_PUT_ELEMS
       
       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!           
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       

       subroutine HPX_SWAP_ELEMS(domain,neighbor)
       
       use dg
       use sizes
       
       implicit none
       
       type (dg_type) :: domain
       type (dg_type) :: neighbor
       
       logical :: neighbor_found
       logical :: domain_found
       integer :: i
       integer :: dindex,nindex       
       integer :: el,dof
       integer :: ncount

       
        
       
       neighbor_found = .false.
       DO i=1,domain%NEIGHPROC_R

         IF (domain%IPROC_R(i) == neighbor%IDPROC) THEN
           nindex = i
           neighbor_found = .true.
           EXIT
         ENDIF
      
       ENDDO   
       
       domain_found = .false.
       DO i = 1,neighbor%NEIGHPROC_S
         
         IF (neighbor%IPROC_S(i) == domain%IDPROC) THEN
           dindex = i
           domain_found = .true.
           EXIT
         ENDIF
       
       ENDDO
             
       IF (neighbor_found .and. domain_found) THEN  
       
!          PRINT*, " PROC: ", domain%IDPROC, " RECEIVING ",domain%NELEMRECV(nindex), " ELEMENTS FROM: ", neighbor%IDPROC
!          PRINT*, (domain%IRECVLOC(el,nindex), el = 1,domain%NELEMRECV(nindex))
!          PRINT*, ""
!          PRINT*, " PROC: ", neighbor%IDPROC, " SENDING ",neighbor%NELEMSEND(dindex), " ELEMENTS TO: ", domain%IDPROC
!          PRINT*, (neighbor%ISENDLOC(el,dindex), el = 1,neighbor%NELEMSEND(dindex))
!          PRINT*, ""         
         
         DO el=1,domain%NELEMRECV(nindex)
           DO dof=1,domain%DOFH
             domain%ZE  (dof,domain%IRECVLOC(el,nindex)  ,domain%IRK+1) =  &
             neighbor%ZE(dof,neighbor%ISENDLOC(el,dindex),domain%IRK+1)
           ENDDO
         ENDDO
       
         DO el=1,domain%NELEMRECV(nindex)
           DO dof=1,domain%DOFH
             domain%QX  (dof,domain%IRECVLOC(el,nindex)  ,domain%IRK+1) = &
             neighbor%QX(dof,neighbor%ISENDLOC(el,dindex),domain%IRK+1)             
           ENDDO
         ENDDO
       
         DO el=1,domain%NELEMRECV(nindex)
           DO dof=1,domain%DOFH
             domain%QY  (dof,domain%IRECVLOC(el,nindex)  ,domain%IRK+1) = &
             neighbor%QY(dof,neighbor%ISENDLOC(el,dindex),domain%IRK+1)             
           ENDDO
         ENDDO    
       
       ELSE
       
         PRINT*, "FORTRAN ERROR: neighbor not found"
         STOP       
         
       ENDIF
         

     
       
       end subroutine HPX_SWAP_ELEMS
              
       
       
#endif
