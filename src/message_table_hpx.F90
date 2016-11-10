#ifdef HPX
      SUBROUTINE ELEM_MSG_TABLE_HPX (s, dg_here) 
!
!--------------------------------------------------------------------------
!  Routine preforms following steps:
!
!   (1) Read Message-Passing Information from file "DG.18"
!   (2) Determine resident nodes: RESNODE(I) is true  if I is resident node
!   (3) Determine ghost nodes:    RESNODE(I) is false if I is ghost node    
!   (4) Determine number of neighbor subdomains
!   (5) MPI rank of each neighbor and number of ghosts nodes to receive
!   (6) Read Message-Passing Receive List
!   (7) MPI rank of each neighbor and number of ghosts nodes to send
!   (8) Read Message-Passing Send List
!
!  vjp  8/06/1999
!--------------------------------------------------------------------------
!
      use sizes
      use dg
      IMPLICIT NONE
      type (sizes_type) :: s
      type (dg_type) :: dg_here
  
      INTEGER :: IDPROC,NLOCAL,I,J,JJ,NEIGH
      INTEGER :: NCOMMELEM_R,NCOMMELEM_S
    !  INTEGER :: RDIM ! can this be local?
  
      OPEN(s%dg18unit,FILE=s%DIRNAME(1:s%LNAME)//'/'//'DG.18')
      READ(s%dg18unit,3010) IDPROC,NLOCAL
  
      dg_here%IDPROC = IDPROC
  
      ALLOCATE ( DG_HERE%NELEMLOC(NLOCAL) )
      READ(s%dg18unit,1130) (DG_HERE%NELEMLOC(I), I=1,NLOCAL)
      ALLOCATE ( DG_HERE%IBELONGTO_ELEM(S%MNE),DG_HERE%RESELEM(S%MNE) )
      DO I=1,S%MNE
         DG_HERE%IBELONGTO_ELEM(I) = 0
      ENDDO
      DO I=1,NLOCAL
         DG_HERE%IBELONGTO_ELEM(DG_HERE%NELEMLOC(I)) = IDPROC + 1
      ENDDO
      DO I=1, S%MNE
         IF (DG_HERE%IBELONGTO_ELEM(I)-1.EQ.s%MYPROC) THEN
            DG_HERE%RESELEM(I) = .TRUE.
         ELSE 
            DG_HERE%RESELEM(I) = .FALSE.
         ENDIF
      ENDDO
  
      READ(s%dg18unit,3010) DG_HERE%NEIGHPROC_R,DG_HERE%NEIGHPROC_S
      ! check if neighproc_r or neighproc_s are larger than MAX_DOMAIN_NEIGHBOR
      IF(DG_HERE%NEIGHPROC_R.GT.MAX_DOMAIN_NEIGHBORS) THEN
         WRITE(*,*) "NEIGHPROC_R GREATER THAN MAX_DOMAIN_NEIGHBOR"
      END IF

      IF(DG_HERE%NEIGHPROC_S.GT.MAX_DOMAIN_NEIGHBORS) THEN
         WRITE(*,*) "NEIGHPROC_S GREATER THAN MAX_DOMAIN_NEIGHBOR"
      END IF
  
      dg_here%RDIM_ELEM = DG_HERE%NEIGHPROC_R + DG_HERE%NEIGHPROC_S
      ALLOCATE( dg_here%INDEX_ELEM(dg_here%RDIM_ELEM) )
  
      ALLOCATE( DG_HERE%IPROC_R(DG_HERE%NEIGHPROC_R),DG_HERE%NELEMRECV(DG_HERE%NEIGHPROC_R) )
      ALLOCATE( DG_HERE%IRECVLOC_ELEM(S%MNE,DG_HERE%NEIGHPROC_R) )
  
      DO JJ=1,DG_HERE%NEIGHPROC_R
         J = MOD(JJ-1+s%MYPROC,DG_HERE%NEIGHPROC_R)+1
         READ(s%dg18unit,3010) DG_HERE%IPROC_R(J),DG_HERE%NELEMRECV(J)
         READ(s%dg18unit,1130) (DG_HERE%IRECVLOC_ELEM(I,J), I=1,DG_HERE%NELEMRECV(J))
      ENDDO
  
      ALLOCATE( DG_HERE%IPROC_S(DG_HERE%NEIGHPROC_S),DG_HERE%NELEMSEND(DG_HERE%NEIGHPROC_S) )
      ALLOCATE( DG_HERE%ISENDLOC_ELEM(S%MNE,DG_HERE%NEIGHPROC_S) )
  
      DO JJ=1,DG_HERE%NEIGHPROC_S
         J = MOD(JJ-1+s%MYPROC,DG_HERE%NEIGHPROC_S)+1
         READ(s%dg18unit,3010) DG_HERE%IPROC_S(J),DG_HERE%NELEMSEND(J)
         READ(s%dg18unit,1130) (DG_HERE%ISENDLOC_ELEM(I,J), I=1,DG_HERE%NELEMSEND(J))
      ENDDO
  

      CLOSE(s%dg18unit)

  
  
      NCOMMELEM_R = 0
      NCOMMELEM_S = 0
      DO J=1,DG_HERE%NEIGHPROC_R
        NCOMMELEM_R = NCOMMELEM_R + DG_HERE%NELEMRECV(J)
      ENDDO
      DO J=1,DG_HERE%NEIGHPROC_S
        NCOMMELEM_S = NCOMMELEM_S + DG_HERE%NELEMSEND(J)
      ENDDO  
  
      DG_HERE%SEND_VOL = NCOMMELEM_S*DG_HERE%DOFH*3  
      DG_HERE%RECV_VOL = NCOMMELEM_R*DG_HERE%DOFH*3    
  
  
     if (dg_here%send_vol .gt. MAX_BUFFER_SIZE) then
         print*, "FORTRAN ERROR: send_vol greater than MAX_BUFFER_SIZE"
         stop
      endif
      if (dg_here%recv_vol .gt. MAX_BUFFER_SIZE) then
         print*, "FORTRAN ERROR: recv_vol greater than MAX_BUFFER_SIZE"
         stop
      endif


      RETURN
  
1130  FORMAT(8X,9I8)
3010  FORMAT(8X,2I8)

      END SUBROUTINE ELEM_MSG_TABLE_HPX


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!           
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       


      SUBROUTINE NODE_MSG_TABLE_HPX (s, dg_here) 
!
!--------------------------------------------------------------------------
!  Routine preforms following steps:
!
!   (1) Read Message-Passing Information from file "fort.18"
!   (2) Determine resident nodes: RESNODE(I) is true  if I is resident node
!   (3) Determine ghost nodes:    RESNODE(I) is false if I is ghost node    
!   (4) Determine number of neighbor subdomains
!   (5) MPI rank of each neighbor and number of ghosts nodes to receive
!   (6) Read Message-Passing Receive List
!   (7) MPI rank of each neighbor and number of ghosts nodes to send
!   (8) Read Message-Passing Send List
!
!  vjp  8/06/1999
!--------------------------------------------------------------------------
!

      USE sizes
      USE dg

      IMPLICIT NONE

      type (sizes_type) :: s
      type (dg_type) :: dg_here
      INTEGER :: IDPROC,NLOCAL,I,J

      OPEN(s%dg18unit,FILE=s%DIRNAME(1:s%LNAME)//'/'//'fort.18')



      READ(s%dg18unit,3010) IDPROC,NLOCAL    

      ALLOCATE ( dg_here%NNODELOC(NLOCAL) )

      READ(s%dg18unit,1130) (dg_here%NNODELOC(I), I=1,NLOCAL)

      ALLOCATE ( dg_here%IBELONGTO_NODE(s%MNP),dg_here%RESNODE(s%MNP) )

      DO I=1,s%MNP
         dg_here%IBELONGTO_NODE(I) = 0
      ENDDO
      DO I=1,NLOCAL
         dg_here%IBELONGTO_NODE(dg_here%NNODELOC(I)) = IDPROC + 1
      ENDDO
      DO I=1, s%MNP
         IF (dg_here%IBELONGTO_NODE(I)-1.EQ.s%MYPROC) THEN
           dg_here%RESNODE(I) = .TRUE.
         ELSE 
           dg_here%RESNODE(I) = .FALSE.
         ENDIF
      ENDDO

      READ(s%dg18unit,3015) dg_here%NEIGHPROC

      IF(DG_HERE%NEIGHPROC.GT.MAX_DOMAIN_NEIGHBORS) THEN
         WRITE(*,*) "NEIGHPROC GREATER THAN MAX_DOMAIN_NEIGHBOR"
      END IF      

      dg_here%RDIM_NODE = 2*dg_here%NEIGHPROC
      ALLOCATE( dg_here%INDEX_NODE(dg_here%RDIM_NODE) )

      ALLOCATE( dg_here%IPROC(dg_here%NEIGHPROC),dg_here%NNODRECV(dg_here%NEIGHPROC) )
      ALLOCATE( dg_here%IRECVLOC_NODE(s%MNP,dg_here%NEIGHPROC) )

      DO J=1,dg_here%NEIGHPROC
         READ(s%dg18unit,3010) dg_here%IPROC(J),dg_here%NNODRECV(J)
         READ(s%dg18unit,1130) (dg_here%IRECVLOC_NODE(I,J), I=1,dg_here%NNODRECV(J))
      ENDDO

      ALLOCATE( dg_here%NNODSEND(dg_here%NEIGHPROC) )
      ALLOCATE( dg_here%ISENDLOC_NODE(s%MNP,dg_here%NEIGHPROC) )

      DO J=1,dg_here%NEIGHPROC
         READ(s%dg18unit,3010) dg_here%IPROC(J),dg_here%NNODSEND(J)
         READ(s%dg18unit,1130) (dg_here%ISENDLOC_NODE(I,J), I=1,dg_here%NNODSEND(J))
      ENDDO

      CLOSE(s%dg18unit)
      RETURN

1130  FORMAT(8X,9I8)
3010  FORMAT(8X,2I8)
3015  FORMAT(8X,I8)
      END SUBROUTINE
#endif
