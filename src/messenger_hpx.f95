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
      
         ncount = 0
         DO el=1,dg_here%NELEMSEND(index)
            DO dof=1,dg_here%DOFH
              ncount = ncount+1
              sendbuf(ncount)=dg_here%ZE(dof,dg_here%ISENDLOC(el,index),dg_here%IRK)
            ENDDO
         ENDDO
       
         DO el=1,dg_here%NELEMSEND(index)
           DO dof=1,dg_here%DOFH
             ncount = ncount+1
             sendbuf(ncount)=dg_here%QX(dof,dg_here%ISENDLOC(el,index),dg_here%IRK)
           ENDDO
         ENDDO
          
         DO el=1,dg_here%NELEMSEND(index)
           DO dof=1,dg_here%DOFH
             ncount = ncount+1
             sendbuf(ncount)=dg_here%QY(dof,dg_here%ISENDLOC(el,index),dg_here%IRK)
           ENDDO
         ENDDO
       
       ELSE 
       
         PRINT*, "FORTRAN ERROR: neighbor not found"
         STOP             
    
       ENDIF
       
       
       
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
       
         ncount = 0
         DO el=1,dg_here%NELEMRECV(index)
           DO dof=1,dg_here%DOFH
             ncount = ncount+1
             dg_here%ZE(dof,dg_here%IRECVLOC(el,index),dg_here%IRK) = recvbuf(ncount)
           ENDDO
         ENDDO
       
         DO el=1,dg_here%NELEMRECV(index)
           DO dof=1,dg_here%DOFH
             ncount = ncount+1
             dg_here%QX(dof,dg_here%IRECVLOC(el,index),dg_here%IRK) = recvbuf(ncount)
           ENDDO
         ENDDO
       
         DO el=1,dg_here%NELEMRECV(index)
           DO dof=1,dg_here%DOFH
             ncount = ncount+1
             dg_here%QY(dof,dg_here%IRECVLOC(el,index),dg_here%IRK) = recvbuf(ncount)
           ENDDO
         ENDDO    
       
       ELSE
       
         PRINT*, "FORTRAN ERROR: neighbor not found"
         STOP       
         
       ENDIF
         

     
       
       end subroutine HPX_PUT_ELEMS
#endif
