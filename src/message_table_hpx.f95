#ifdef HPX
SUBROUTINE MSG_TABLE_HPX (s, dg_here) 
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
  INTEGER :: RDIM ! can this be local?
  
  OPEN(dg18unit,FILE=s%DIRNAME(1:s%LNAME)//'/'//'DG.18')
  READ(dg18unit,3010) IDPROC,NLOCAL
  ALLOCATE ( DG_HERE%NELEMLOC(NLOCAL) )
  READ(dg18unit,1130) (DG_HERE%NELEMLOC(I), I=1,NLOCAL)
  ALLOCATE ( DG_HERE%IBELONGTO(S%MNE),DG_HERE%RESELEM(S%MNE) )
  DO I=1,S%MNE
     DG_HERE%IBELONGTO(I) = 0
  ENDDO
  DO I=1,NLOCAL
     DG_HERE%IBELONGTO(DG_HERE%NELEMLOC(I)) = IDPROC + 1
  ENDDO
  DO I=1, S%MNE
     IF (DG_HERE%IBELONGTO(I)-1.EQ.s%MYPROC) THEN
        DG_HERE%RESELEM(I) = .TRUE.
     ELSE 
        DG_HERE%RESELEM(I) = .FALSE.
     ENDIF
  ENDDO
  
  READ(dg18unit,3010) DG_HERE%NEIGHPROC_R,DG_HERE%NEIGHPROC_S
  
  RDIM = DG_HERE%NEIGHPROC_R + DG_HERE%NEIGHPROC_S
  ALLOCATE( INDEX(RDIM) )
  
  ALLOCATE( DG_HERE%IPROC_R(DG_HERE%NEIGHPROC_R),DG_HERE%NELEMRECV(DG_HERE%NEIGHPROC_R) )
  ALLOCATE( DG_HERE%IRECVLOC(S%MNE,DG_HERE%NEIGHPROC_R) )
  
  DO JJ=1,DG_HERE%NEIGHPROC_R
     J = MOD(JJ-1+s%MYPROC,DG_HERE%NEIGHPROC_R)+1
     READ(dg18unit,3010) DG_HERE%IPROC_R(J),DG_HERE%NELEMRECV(J)
     READ(dg18unit,1130) (DG_HERE%IRECVLOC(I,J), I=1,DG_HERE%NELEMRECV(J))
  ENDDO
  
  ALLOCATE( DG_HERE%IPROC_S(DG_HERE%NEIGHPROC_S),DG_HERE%NELEMSEND(DG_HERE%NEIGHPROC_S) )
  ALLOCATE( DG_HERE%ISENDLOC(S%MNE,DG_HERE%NEIGHPROC_S) )
  
  DO JJ=1,DG_HERE%NEIGHPROC_S
     J = MOD(JJ-1+s%MYPROC,DG_HERE%NEIGHPROC_S)+1
     READ(dg18unit,3010) DG_HERE%IPROC_S(J),DG_HERE%NELEMSEND(J)
     READ(dg18unit,1130) (DG_HERE%ISENDLOC(I,J), I=1,DG_HERE%NELEMSEND(J))
  ENDDO
  
  CLOSE(dg18unit)
  RETURN
  
1130 FORMAT(8X,9I8)
3010 FORMAT(8X,2I8)
END SUBROUTINE MSG_TABLE_HPX
#endif
