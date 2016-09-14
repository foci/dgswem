#ifdef HPX
       subroutine HPX_GET_ELEMS(dg_here,domain,volume,sendbuf,rkindex)     
       
       use dg
       use sizes
       
       implicit none       

       type (dg_type) :: dg_here
       
       logical :: neighbor_found
       integer :: i,index
       integer :: el,dof
       integer :: ncount
       integer :: domain
       integer :: volume
       integer :: rkindex
       
!       real(sz) :: sendbuf(volume)
       real(sz) :: sendbuf(MAX_BUFFER_SIZE)   

       if (rkindex.eq.0) then
          rkindex = dg_here%irk+1
       endif
       
!       PRINT*, "HPX_GET_ELEMS, rkindex = ", rkindex
       
       neighbor_found = .false.
       DO i=1,dg_here%NEIGHPROC_S

         IF (dg_here%IPROC_S(i) == domain) THEN
           index = i
           neighbor_found = .true.
           EXIT
         ENDIF
      
       ENDDO     
       
       IF (neighbor_found) THEN
       
!          PRINT*, " PROC: ", dg_here%IDPROC, " SENDING ", dg_here%NELEMSEND(index) ," ELEMENTS TO: ", dg_here%IPROC_S(index)       
!          PRINT 180, (dg_here%ISENDLOC_ELEM(el,index), el = 1,dg_here%NELEMSEND(index))
!          PRINT*, "" 
      
         ncount = 0
         DO el=1,dg_here%NELEMSEND(index)
            DO dof=1,dg_here%DOFH
              ncount = ncount+1
              sendbuf(ncount)=dg_here%ZE(dof,dg_here%ISENDLOC_ELEM(el,index),rkindex)
!              IF (abs(sendbuf(ncount) ) > 1d-14) THEN              
!               PRINT 181, dg_here%ISENDLOC_ELEM(el,index),sendbuf(ncount)
!              ENDIF
            ENDDO
         ENDDO
       
         DO el=1,dg_here%NELEMSEND(index)
           DO dof=1,dg_here%DOFH
             ncount = ncount+1
             sendbuf(ncount)=dg_here%QX(dof,dg_here%ISENDLOC_ELEM(el,index),rkindex)
!              IF (abs(sendbuf(ncount) ) > 1d-14) THEN
!                PRINT 181, dg_here%ISENDLOC_ELEM(el,index),sendbuf(ncount)
!              ENDIF
           ENDDO
         ENDDO
          
         DO el=1,dg_here%NELEMSEND(index)
           DO dof=1,dg_here%DOFH
             ncount = ncount+1
             sendbuf(ncount)=dg_here%QY(dof,dg_here%ISENDLOC_ELEM(el,index),rkindex)
!              IF (abs(sendbuf(ncount) ) > 1d-14) THEN                      
!                PRINT 181, dg_here%ISENDLOC_ELEM(el,index),sendbuf(ncount)
!              ENDIF
           ENDDO
         ENDDO
         
         volume = ncount

!          PRINT*, "----------------------------------------------------------------------------------"                  
       
       ELSE 
       
!          PRINT*, "FORTRAN ERROR: send neighbor not found"   
!          PRINT*, "SEND PROC: ", dg_here%IDPROC, " DOMAIN: ", domain                 
!          STOP             
    
       ENDIF
       
 180  FORMAT(8X,9I8)    
 181  FORMAT(I5,2X,ES24.17)
       
       end subroutine HPX_GET_ELEMS
       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!           
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       

       subroutine HPX_PUT_ELEMS(dg_here,neighbor,volume,recvbuf,rkindex)
       
       use dg
       use sizes
       
       implicit none
       
       type (dg_type) :: dg_here
       
       logical :: neighbor_found
       integer :: i,l,index       
       integer :: el,dof
       integer :: ncount
       integer :: neighbor
       integer :: volume
       integer :: rkindex
       
!       real(sz) :: recvbuf(volume)
       real(sz) :: recvbuf(MAX_BUFFER_SIZE)
        
       if (rkindex.eq.0) then
          rkindex = dg_here%irk+1
       endif

!       PRINT*, "HPX_PUT_ELEMS, rkindex = ", rkindex
       
       neighbor_found = .false.
       DO i=1,dg_here%NEIGHPROC_R

         IF (dg_here%IPROC_R(i) == neighbor) THEN
           index = i
           neighbor_found = .true.
           EXIT
         ENDIF
      
       ENDDO           
             
       IF (neighbor_found) THEN  
       
!          PRINT*, " PROC: ", dg_here%IDPROC, " RECEIVING ",dg_here%NELEMRECV(index), " ELEMENTS FROM: ", dg_here%IPROC_R(index)
!          PRINT 180, (dg_here%IRECVLOC_ELEM(el,index), el = 1,dg_here%NELEMRECV(index))
!          PRINT*, ""
         
         ncount = 0
         DO el=1,dg_here%NELEMRECV(index)
           DO dof=1,dg_here%DOFH
             ncount = ncount+1
             dg_here%ZE(dof,dg_here%IRECVLOC_ELEM(el,index),rkindex) = recvbuf(ncount)
!              IF (abs(recvbuf(ncount) ) > 1d-14) THEN             
!                PRINT 181, dg_here%IRECVLOC_ELEM(el,index),recvbuf(ncount) 
!              ENDIF
           ENDDO
         ENDDO
       
         DO el=1,dg_here%NELEMRECV(index)
           DO dof=1,dg_here%DOFH
             ncount = ncount+1
             dg_here%QX(dof,dg_here%IRECVLOC_ELEM(el,index),rkindex) = recvbuf(ncount)
!              IF (abs(recvbuf(ncount) ) > 1d-14) THEN
!                PRINT 181, dg_here%IRECVLOC_ELEM(el,index),recvbuf(ncount)
!              ENDIF
           ENDDO
         ENDDO
       
         DO el=1,dg_here%NELEMRECV(index)
           DO dof=1,dg_here%DOFH
             ncount = ncount+1
             dg_here%QY(dof,dg_here%IRECVLOC_ELEM(el,index),rkindex) = recvbuf(ncount)
!              IF (abs(recvbuf(ncount) ) > 1d-14) THEN                
!                PRINT 181, dg_here%IRECVLOC_ELEM(el,index),recvbuf(ncount)
!              ENDIF
           ENDDO
         ENDDO    
         
!          PRINT*, "----------------------------------------------------------------------------------"         
       
       ELSE
       
!          PRINT*, "FORTRAN ERROR: recieve neighbor not found"
!          PRINT*, "RECV PROC: ", dg_here%IDPROC, " NEIGHBOR: ", neighbor 
!          STOP       
         
       ENDIF
         

 180  FORMAT(8X,9I8)     
 181  FORMAT(I5,2X,ES24.17) 
       
       end subroutine HPX_PUT_ELEMS
       
       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!           
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      subroutine HPX_GET_NODES(dg_here,domain,volume,sendbuf)
       
       use dg
       use sizes
       
       implicit none       

       type (dg_type) :: dg_here
       
       logical :: neighbor_found
       integer :: i,index
       integer :: nd
       integer :: ncount
       integer :: domain
       integer :: volume
       
!       real(sz) :: sendbuf(volume)
       real(sz) :: sendbuf(MAX_BUFFER_SIZE)   

       
       neighbor_found = .false.
       DO i=1,dg_here%NEIGHPROC

         IF (dg_here%IPROC(i) == domain) THEN
           index = i
           neighbor_found = .true.
           EXIT
         ENDIF
      
       ENDDO     
       
       IF (neighbor_found) THEN
       
!          PRINT*, " PROC: ", dg_here%IDPROC, " SENDING ", dg_here%NNODESEND(index) ," NODES TO: ", dg_here%IPROC(index)       
!          PRINT 180, (dg_here%ISENDLOC_NODE(nd,index), nd = 1,dg_here%NNODESEND(index))
!          PRINT*, "" 
      
         ncount = 0
         DO nd=1,dg_here%NNODSEND(index)
            ncount = ncount+1
            sendbuf(ncount)=dg_here%ZE_MIN1(dg_here%ISENDLOC_NODE(nd,index))
         ENDDO

         DO nd=1,dg_here%NNODSEND(index)
            ncount = ncount+1
            sendbuf(ncount)=dg_here%ZE_MAX1(dg_here%ISENDLOC_NODE(nd,index))
         ENDDO

         DO nd=1,dg_here%NNODSEND(index)
            ncount = ncount+1
            sendbuf(ncount)=dg_here%QX_MIN1(dg_here%ISENDLOC_NODE(nd,index))
         ENDDO

         DO nd=1,dg_here%NNODSEND(index)
            ncount = ncount+1
            sendbuf(ncount)=dg_here%QX_MAX1(dg_here%ISENDLOC_NODE(nd,index))
         ENDDO

         DO nd=1,dg_here%NNODSEND(index)
            ncount = ncount+1
            sendbuf(ncount)=dg_here%QY_MIN1(dg_here%ISENDLOC_NODE(nd,index))
         ENDDO

         DO nd=1,dg_here%NNODSEND(index)
            ncount = ncount+1
            sendbuf(ncount)=dg_here%QY_MAX1(dg_here%ISENDLOC_NODE(nd,index))
         ENDDO
       
         
         volume = ncount

!          PRINT*, "----------------------------------------------------------------------------------"                  
       
       ELSE 
       
!          PRINT*, "FORTRAN ERROR: send neighbor not found"   
!          PRINT*, "SEND PROC: ", dg_here%IDPROC, " DOMAIN: ", domain                 
!          STOP             
    
       ENDIF
       
 180  FORMAT(8X,9I8)    
 181  FORMAT(I5,2X,ES24.17)

      return
      end subroutine HPX_GET_NODES

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!           
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

       subroutine HPX_PUT_NODES(dg_here,neighbor,volume,recvbuf)

       use dg
       use sizes
       
       implicit none
       
       type (dg_type) :: dg_here
       
       logical :: neighbor_found
       integer :: i,index       
       integer :: nd
       integer :: ncount
       integer :: neighbor
       integer :: volume
       
!       real(sz) :: recvbuf(volume)
       real(sz) :: recvbuf(MAX_BUFFER_SIZE)
        

       neighbor_found = .false.
       DO i=1,dg_here%NEIGHPROC

         IF (dg_here%IPROC(i) == neighbor) THEN
           index = i
           neighbor_found = .true.
           EXIT
         ENDIF
      
       ENDDO           
             
       IF (neighbor_found) THEN  
       
!          PRINT*, " PROC: ", dg_here%IDPROC, " RECEIVING ",dg_here%NNODERECV(index), " NODES FROM: ", dg_here%IPROC(index)
!          PRINT 180, (dg_here%IRECVLOC_NODE(el,index), el = 1,dg_here%NNODERECV(index))
!          PRINT*, ""
         
         ncount = 0
         DO nd=1,dg_here%NNODRECV(index)
           ncount = ncount+1
           dg_here%ZE_MIN1(dg_here%IRECVLOC_NODE(nd,index)) = recvbuf(ncount)
         ENDDO

         DO nd=1,dg_here%NNODRECV(index)
           ncount = ncount+1
           dg_here%ZE_MAX1(dg_here%IRECVLOC_NODE(nd,index)) = recvbuf(ncount)
         ENDDO

         DO nd=1,dg_here%NNODRECV(index)
           ncount = ncount+1
           dg_here%QX_MIN1(dg_here%IRECVLOC_NODE(nd,index)) = recvbuf(ncount)
         ENDDO

         DO nd=1,dg_here%NNODRECV(index)
           ncount = ncount+1
           dg_here%QX_MAX1(dg_here%IRECVLOC_NODE(nd,index)) = recvbuf(ncount)
         ENDDO

         DO nd=1,dg_here%NNODRECV(index)
           ncount = ncount+1
           dg_here%QY_MIN1(dg_here%IRECVLOC_NODE(nd,index)) = recvbuf(ncount)
         ENDDO

         DO nd=1,dg_here%NNODRECV(index)
           ncount = ncount+1
           dg_here%QY_MAX1(dg_here%IRECVLOC_NODE(nd,index)) = recvbuf(ncount)
         ENDDO

       
        
!          PRINT*, "----------------------------------------------------------------------------------"         
       
       ELSE
       
!          PRINT*, "FORTRAN ERROR: recieve neighbor not found"
!          PRINT*, "RECV PROC: ", dg_here%IDPROC, " NEIGHBOR: ", neighbor 
!          STOP       
         
       ENDIF
         

 180  FORMAT(8X,9I8)     
 181  FORMAT(I5,2X,ES24.17) 

       return
       end subroutine HPX_PUT_NODES


#endif
