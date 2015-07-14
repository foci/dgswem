#ifdef HPX
       subroutine HPX_GET_ELEMS(dg_here,neighbor,volume,sendbuf)     
       
       use global
       use dg
       use sizes
       
       implicit none       

       type (sizes_type) :: s
       type (dg_type) :: dg_here
       type (global_type) :: global_here
       
       integer :: el,k
       integer :: ncount
       
       real(sz) :: sendbuf(volume)
       
      
       ncount = 0
       DO el=1,dg_here%NELEMSEND(neighbor)
          DO dof=1,dg_here%DOFH
            ncount = ncount+1
            sendbuf(ncount)=dg_here%ZE(dof,dg_here%ISENDLOC(el,neighbor),dg_here%IRK)
          ENDDO
       ENDDO
       
       DO el=1,dg_here%NELEMSEND(neighbor)
         DO dof=1,dg_here%DOFH
           ncount = ncount+1
           sendbuf(ncount)=dg_here%QX(dof,dg_here%ISENDLOC(el,neighbor),dg_here%IRK)
         ENDDO
       ENDDO
          
       DO el=1,dg_here%NELEMSEND(neighbor)
         DO dof=1,dg_here%DOFH
           ncount = ncount+1
           sendbuf(ncount)=dg_here%QY(dof,dg_here%ISENDLOC(el,neighbor),dg_here%IRK)
         ENDDO
       ENDDO
    
       
       
       end subroutine HPX_GET_ELEMS
       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!           
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       

       subroutine HPX_PUT_ELEMS(dg_here,neighbor,volume,recvbuf)
       
       use global
       use dg
       use sizes
       
       implicit none
       
       type (dg_type) :: dg_here
       
       integer :: i,k
       integer :: ncount
       integer :: neighbor
       integer :: volume
       
       real(sz) :: recvbuf(volume)      
       
             
       
       ncount = 0
       DO el=1,dg_here%NELEMRECV(neighbor)
         DO dof=1,dg_here%DOFH
           ncount = ncount+1
           dg_here%ZE(dof,dg_here%IRECVLOC(el,neighbor),dg_here%IRK) = recvbuf(ncount)
         ENDDO
       ENDDO
       
       DO el=1,dg_here%NELEMRECV(neighbor)
         DO dof=1,dg_here%DOFH
           ncount = ncount+1
           dg_here%QX(dof,dg_here%IRECVLOC(el,neighbor),dg_here%IRK) = recvbuf(ncount)
         ENDDO
       ENDDO
       
       DO el=1,dg_here%NELEMRECV(neighbor)
         DO dof=1,dg_here%DOFH
           ncount = ncount+1
           dg_here%QY(dof,dg_here%IRECVLOC(el,neighbor),dg_here%IRK) = recvbuf(ncount)
         ENDDO
       ENDDO         
         

     
       
       end subroutine HPX_PUT_ELEMS
#endif
