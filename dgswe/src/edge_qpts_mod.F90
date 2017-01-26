      MODULE edge_qpts_mod
      
      USE globals, ONLY: rp
      
      IMPLICIT NONE
      
      CONTAINS
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
      
      SUBROUTINE edge_qpts(myrank,p,ctp,nel_type,nqpte,mnqpte,wpte,qpte)

      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: myrank
      INTEGER, INTENT(IN) :: p
      INTEGER, INTENT(IN) :: ctp
      INTEGER, INTENT(IN) :: nel_type
      INTEGER, DIMENSION(:), INTENT(OUT) :: nqpte
      INTEGER, INTENT(OUT) :: mnqpte
      REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: wpte
      REAL(rp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(OUT) :: qpte      
      
      INTEGER :: pt,et,led,i
      INTEGER :: order(nel_type),npt(nel_type)
      REAL(rp) :: w(13,nel_type),r(13,nel_type)
      
      order(1) = p+1
      order(2) = order(1)
      order(3) = p+ctp+1
      order(4) = order(3)

!       order(1) = p+ctp+1
!       order(2) = order(1)
!       order(3) = order(1)
!       order(4) = order(1)
      
      DO et = 1,nel_type
        CALL gauss_qpts(order(et),npt(et),w(:,et),r(:,et))
      ENDDO
      
      
      mnqpte = maxval(npt)
      
      
      ALLOCATE(qpte(4*mnqpte,2,nel_type))        
      ALLOCATE(wpte(4*mnqpte,nel_type))
   
   
      DO et = 1,nel_type 
      
        nqpte(et) = npt(et)              
        
        IF (mod(et,2) == 1) THEN
          
          DO led = 1,3
            DO i = 1,nqpte(et)
              pt = (led-1)*nqpte(et) + i
              SELECT CASE(led)
                CASE(1)
                  wpte(pt,et)   = w(i,et)
                  qpte(pt,1,et) = -r(i,et)
                  qpte(pt,2,et) =  r(i,et)
                CASE(2)
                  wpte(pt,et)   = w(i,et)
                  qpte(pt,1,et) = -1d0
                  qpte(pt,2,et) = -r(i,et)
                CASE(3)
                  wpte(pt,et)   = w(i,et)
                  qpte(pt,1,et) = r(i,et)
                  qpte(pt,2,et) = -1d0
              END SELECT
            ENDDO
          ENDDO          
        
        ELSE IF (mod(et,2) == 0) THEN

          DO led = 1,4
            DO i = 1,nqpte(et)
              pt = (led-1)*nqpte(et) + i
              SELECT CASE(led)
                CASE(1)
                  wpte(pt,et)   = w(i,et)
                  qpte(pt,1,et) = 1d0
                  qpte(pt,2,et) = r(i,et)
                CASE(2)
                  wpte(pt,et)   = w(i,et)
                  qpte(pt,1,et) = -r(i,et)
                  qpte(pt,2,et) = 1d0
                CASE(3)
                  wpte(pt,et)   = w(i,et)
                  qpte(pt,1,et) = -1d0
                  qpte(pt,2,et) = -r(i,et)
                CASE(4)
                  wpte(pt,et)   = w(i,et)
                  qpte(pt,1,et) = r(i,et)
                  qpte(pt,2,et) = -1d0 
              END SELECT
            ENDDO
          ENDDO          
        
        ENDIF
        
      ENDDO


      IF (myrank == 0) THEN
        PRINT "(A)", "---------------------------------------------"
        PRINT "(A)", "        Edge Integration Information         "
        PRINT "(A)", "---------------------------------------------"
        PRINT "(A)", " "

!         PRINT "(A)", 'Edge quadrature weights and points'
      
        DO i = 1,nel_type
          PRINT "(A,I3)", "Number of edge quadrature points",nqpte(i)      
!           DO pt = 1,nqpte(i)
!             PRINT "(2(F10.3))", wpte(pt,i),qpte(pt,2,i)
!           ENDDO
!           PRINT*, ' '
        ENDDO
        PRINT*, " " 
      ENDIF

      RETURN
      END SUBROUTINE edge_qpts
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      SUBROUTINE gauss_qpts(p,nqpte,wpte,qpte)

      USE globals, ONLY: rp
      
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: p
      INTEGER, INTENT(OUT) :: nqpte
      REAL(rp), DIMENSION(:), INTENT(OUT) :: wpte
      REAL(rp), DIMENSION(:), INTENT(OUT) :: qpte

      wpte = 0d0
      qpte = 0d0
      
      SELECT CASE(p)

      CASE(2)
         nqpte = 2

!          qpte(1) = -0.57735026918963D0
!          qpte(2) =  0.57735026918963D0
         qpte(1) = -0.577350269189625765D0
         qpte(2) =  0.577350269189625765D0
         
         wpte(1) = 1.D0
         wpte(2) = 1.D0
         
      CASE(3)
         nqpte = 3
         
!          qpte(1) = -0.77459666924148D0
!          qpte(2) =  0.D0
!          qpte(3) =  0.77459666924148D0
         qpte(1) = -0.774596669241483377D0
         qpte(2) =  0.D0
         qpte(3) =  0.774596669241483377D0
         
!          wpte(1) =  0.55555555555556D0
!          wpte(2) =  0.88888888888888D0
!          wpte(3) =  0.55555555555556D0
         wpte(1) =  0.555555555555555556D0
         wpte(2) =  0.888888888888888889D0
         wpte(3) =  0.555555555555555556D0         
         
      CASE(4)
         nqpte = 4
         
!          qpte(1) = -0.86113631159405D0
!          qpte(2) = -0.33998104358486D0
!          qpte(3) =  0.33998104358486D0
!          qpte(4) =  0.86113631159405D0
         qpte(1) = -0.861136311594052575D0
         qpte(2) = -0.339981043584856265D0
         qpte(3) =  0.339981043584856265D0
         qpte(4) =  0.861136311594052575D0
         
!          wpte(1) =  0.34785484513745D0
!          wpte(2) =  0.65214515486255D0
!          wpte(3) =  0.65214515486255D0
!          wpte(4) =  0.34785484513745D0
         wpte(1) =  0.347854845137453857D0
         wpte(2) =  0.652145154862546143D0
         wpte(3) =  0.652145154862546143D0
         wpte(4) =  0.347854845137453857D0         

      CASE(5)
         nqpte = 5

         qpte(1) = -0.90617984593866D0
         qpte(2) = -0.53846931010568D0
         qpte(3) =  0.D0
         qpte(4) =  0.53846931010568D0
         qpte(5) =  0.90617984593866D0
         
         wpte(1) =  0.23692688505619D0
         wpte(2) =  0.47862867049937D0
         wpte(3) =  0.56888888888889D0
         wpte(4) =  0.47862867049937D0
         wpte(5) =  0.23692688505619D0

      CASE(6)
         nqpte = 6

         qpte(1) = -0.93246951420315D0
         qpte(2) = -0.66120938646626D0
         qpte(3) = -0.23861918608320D0
         qpte(4) =  0.23861918608320D0
         qpte(5) =  0.66120938646626D0
         qpte(6) =  0.93246951420315D0
         
         wpte(1) =  0.17132449237917D0
         wpte(2) =  0.36076157304814D0
         wpte(3) =  0.46791393457269D0
         wpte(4) =  0.46791393457269D0
         wpte(5) =  0.36076157304814D0
         wpte(6) =  0.17132449237917D0

      CASE(7)
         nqpte = 7
         
         qpte(1) = -0.94910791234276D0
         qpte(2) = -0.74153118559939D0
         qpte(3) = -0.40584515137740D0
         qpte(4) =  0.D0
         qpte(5) =  0.40584515137740D0
         qpte(6) =  0.74153118559939D0
         qpte(7) =  0.94910791234276D0
         
         wpte(1) =  0.12948496616887D0
         wpte(2) =  0.27970539148928D0
         wpte(3) =  0.38183005050512D0
         wpte(4) =  0.41795918367347D0
         wpte(5) =  0.38183005050512D0
         wpte(6) =  0.27970539148928D0
         wpte(7) =  0.12948496616887D0

      CASE(8)
         nqpte = 8
         
         qpte(1) = -0.96028985649754D0
         qpte(2) = -0.79666647741363D0
         qpte(3) = -0.52553240991633D0
         qpte(4) = -0.18343464249565D0
         qpte(5) =  0.18343464249565D0
         qpte(6) =  0.52553240991633D0
         qpte(7) =  0.79666647741363D0
         qpte(8) =  0.96028985649754D0
         
         wpte(1) =  0.10122853629038D0
         wpte(2) =  0.22238103445337D0
         wpte(3) =  0.31370664587789D0
         wpte(4) =  0.36268378337836D0
         wpte(5) =  0.36268378337836D0
         wpte(6) =  0.31370664587789D0
         wpte(7) =  0.22238103445337D0
         wpte(8) =  0.10122853629038D0

      CASE(9)
         nqpte = 9
         
         qpte(1) = -0.96816023950763D0
         qpte(2) = -0.83603110732664D0
         qpte(3) = -0.61337143270059D0
         qpte(4) = -0.32425342340381D0
         qpte(5) =  0.D0
         qpte(6) =  0.32425342340381D0
         qpte(7) =  0.61337143270059D0
         qpte(8) =  0.83603110732664D0
         qpte(9) =  0.96816023950763D0
         
         wpte(1) = 0.08127438836163D0
         wpte(2) = 0.18064816069483D0
         wpte(3) = 0.26061069640294D0
         wpte(4) = 0.31234707704000D0
         wpte(5) = 0.33023935500126D0
         wpte(6) = 0.31234707704000D0
         wpte(7) = 0.26061069640294D0
         wpte(8) = 0.18064816069483D0
         wpte(9) = 0.08127438836163D0
         
      CASE(10)
         nqpte = 10
         
         qpte(1)  = -0.97390652851717D0
         qpte(2)  = -0.86506336668898D0
         qpte(3)  = -0.67940956829902D0
         qpte(4)  = -0.43339539412925D0
         qpte(5)  = -0.14887433898163D0
         qpte(6)  =  0.14887433898163D0
         qpte(7)  =  0.43339539412925D0
         qpte(8)  =  0.67940956829902D0
         qpte(9)  =  0.86506336668898D0
         qpte(10) =  0.97390652851717D0

         wpte(1)  =  0.06667134430869D0
         wpte(2)  =  0.14945134915058D0
         wpte(3)  =  0.21908636251598D0
         wpte(4)  =  0.26926671931000D0
         wpte(5)  =  0.29552422471475D0
         wpte(6)  =  0.29552422471475D0
         wpte(7)  =  0.26926671931000D0
         wpte(8)  =  0.21908636251598D0
         wpte(9)  =  0.14945134915058D0
         wpte(10) =  0.06667134430869D0

      CASE(11)
         nqpte = 11
         
         qpte(1)  = -0.97822865814606D0
         qpte(2)  = -0.88706259976810D0
         qpte(3)  = -0.73015200557405D0
         qpte(4)  = -0.51909612920681D0
         qpte(5)  = -0.26954315595234D0
         qpte(6)  =  0.D0
         qpte(7)  =  0.26954315595234D0
         qpte(8)  =  0.51909612920681D0
         qpte(9)  =  0.73015200557405D0
         qpte(10) =  0.88706259976810D0
         qpte(11) =  0.97822865814606D0
         
         wpte(1)  =  0.05566856711627D0
         wpte(2)  =  0.12558036946485D0
         wpte(3)  =  0.18629021092774D0
         wpte(4)  =  0.23319376459199D0
         wpte(5)  =  0.26280454451025D0
         wpte(6)  =  0.27292508677790D0
         wpte(7)  =  0.26280454451025D0
         wpte(8)  =  0.23319376459199D0
         wpte(9)  =  0.18629021092774D0
         wpte(10) =  0.12558036946485D0
         wpte(11) =  0.05566856711627D0
         
      CASE(12)
         nqpte = 12
         
         qpte(1)  = -0.98156063424672D0
         qpte(2)  = -0.90411725637047D0
         qpte(3)  = -0.76990267419430D0
         qpte(4)  = -0.58731795428662D0
         qpte(5)  = -0.36783149899818D0
         qpte(6)  = -0.12523340851147D0
         qpte(7)  =  0.12523340851147D0
         qpte(8)  =  0.36783149899818D0
         qpte(9)  =  0.58731795428662D0
         qpte(10) =  0.76990267419430D0
         qpte(11) =  0.90411725637047D0
         qpte(12) =  0.98156063424672D0
         
         wpte(1)  =  0.04717533638677D0
         wpte(2)  =  0.10693932599520D0
         wpte(3)  =  0.16007832854334D0
         wpte(4)  =  0.20316742672308D0
         wpte(5)  =  0.23349253653835D0
         wpte(6)  =  0.24914704581340D0
         wpte(7)  =  0.24914704581340D0
         wpte(8)  =  0.23349253653835D0
         wpte(9)  =  0.20316742672308D0
         wpte(10) =  0.16007832854334D0
         wpte(11) =  0.10693932599520D0
         wpte(12) =  0.04717533638677D0

      CASE(13)
         nqpte = 13
         
         qpte(1)  = -0.98418305471859D0
         qpte(2)  = -0.91759839922298D0
         qpte(3)  = -0.80157809073331D0
         qpte(4)  = -0.64234933944034D0
         qpte(5)  = -0.44849275103645D0
         qpte(6)  = -0.23045831595513D0
         qpte(7)  =  0.D0
         qpte(8)  =  0.23045831595513D0
         qpte(9)  =  0.44849275103645D0
         qpte(10) =  0.64234933944034D0
         qpte(11) =  0.80157809073331D0
         qpte(12) =  0.91759839922298D0
         qpte(13) =  0.98418305471859D0

         wpte(1)  =  0.04048400476532D0
         wpte(2)  =  0.09212149983773D0
         wpte(3)  =  0.13887351021979D0
         wpte(4)  =  0.17814598076195D0
         wpte(5)  =  0.20781604753689D0
         wpte(6)  =  0.22628318026290D0
         wpte(7)  =  0.23255155323087D0
         wpte(8)  =  0.22628318026290D0
         wpte(9)  =  0.20781604753689D0
         wpte(10) =  0.17814598076195D0
         wpte(11) =  0.13887351021979D0
         wpte(12) =  0.09212149983773D0
         wpte(13) =  0.04048400476532D0

       CASE DEFAULT

         PRINT*, 'Error: invalid p'
         PRINT*, 'Edge quadrature points not available'
         STOP

      END SELECT


      RETURN
      END SUBROUTINE gauss_qpts      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
      END MODULE edge_qpts_mod