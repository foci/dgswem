subroutine hpx_init(g)
  use global_data
  implicit none

  type (global_data_type) :: g
  
  integer :: i,ie

  CALL MAKE_DIRNAME(g)       ! Establish Working Directory Name
  CALL READ_INPUT(g)         ! Establish sizes by reading fort.14 and fort.15

  CALL COLDSTART(g)
  
  !...  
  !...  DETERMINE THE NUMBER OF ACTIVE ELEMENTS (MJU) and total number of 
  !...  elements (NODELE) ATTACHED TO EACH NODE
  !...  

  DO I=1,g%NP
     g%MJU(I)=0
     g%NODELE(I)=0
     g%NODECODE(I)=g%NNODECODE(I)
  END DO
  
  DO IE=1,g%NE
     g%NM1=g%NM(IE,1)
     g%NM2=g%NM(IE,2)
     g%NM3=g%NM(IE,3)
     g%NCELE=g%NODECODE(g%NM1)*g%NODECODE(g%NM2)*g%NODECODE(g%NM3)
     g%MJU(g%NM1)=g%MJU(g%NM1)+g%NCELE
     g%MJU(g%NM2)=g%MJU(g%NM2)+g%NCELE
     g%MJU(g%NM3)=g%MJU(g%NM3)+g%NCELE
     g%NODELE(g%NM1)=g%NODELE(g%NM1)+1
     g%NODELE(g%NM2)=g%NODELE(g%NM2)+1
     g%NODELE(g%NM3)=g%NODELE(g%NM3)+1
  END DO
  
  DO I=1,g%NP
     IF(g%MJU(I).EQ.0) g%MJU(I)=1
  END DO
  
  
  !...  
  !...  ************* SET FLAGS AND COEFFICIENTS USED IN TIME STEPPING *******
  !...  
  
  !...  NONLINEAR FLAGS
  
  IF(g%NOLIBF.EQ.0) THEN
     g%IFNLBF=0
     g%IFLINBF=1
     g%IFHYBF=0
  ENDIF
  IF(g%NOLIBF.EQ.1) THEN
     g%IFNLBF=1
     g%IFLINBF=0
     g%IFHYBF=0
  ENDIF
  IF(g%NOLIBF.EQ.2) THEN
     g%IFNLBF=0
     g%IFLINBF=0
     g%IFHYBF=1
  ENDIF
  IF(g%NOLIFA.EQ.0) THEN
     g%IFNLFA=0
  ELSE
     g%IFNLFA=1
  ENDIF
  IF(g%NOLICA.EQ.0) THEN
     g%IFNLCT=0
     g%NLEQ = 0.D0
     g%LEQ = 1.D0
  ELSE
     g%IFNLCT=1
     g%NLEQ = 1.D0
     g%LEQ = 0.D0
  ENDIF
  IF(g%NOLICAT.EQ.0) THEN
     g%IFNLCAT=0
     g%NLEQ = 0.D0
     g%LEQ = 1.D0
  ELSE
     g%IFNLCAT=1
     g%NLEQ = 1.D0
     g%LEQ = 0.D0
  ENDIF
  g%NLEQG = g%NLEQ*g%G
  g%FG_L = g%LEQ*g%G
  
  g%IFWIND=1
  IF(g%IM.EQ.1) g%IFWIND=0
  
  !...  CONSTANT COEFFICIENTS
  g%GA00=g%G*g%A00
  g%GC00=g%G*g%C00
  g%TADVODT=g%IFNLCAT/g%DT
  g%GB00A00=g%G*(g%B00+g%A00)
  g%GFAO2=g%G*g%IFNLFA/2.D0
  g%GO3=g%G/3.D0
  g%DTO2=g%DT/2.D0
  g%DT2=g%DT*2.D0
  g%GDTO2=g%G*g%DT/2.D0
  g%SADVDTO3=g%IFNLCT*g%DT/3.D0  
  
  CALL PREP_DG(g)

  CALL WRITE_RESULTS(g,0,.FALSE.)

!  CALL WRITE_DG_IC()

end subroutine hpx_init
