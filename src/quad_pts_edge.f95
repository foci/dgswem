C***********************************************************************
C     
C     SUBROUTINE QUAD_PTS_EDGE()
C     
C     Stores the Gauss quadrature points and weights for up to 13 point
C     quadrature rules on the interval (-1,1) (an n point Legendre-Gauss
C     quadrature integrates a 2*n-1 degree polynomial exactly)
C     
C     Written by Ethan Kubatko (06-08-2004)
C     
C     Modifications for parallel runs, Shintaro Bunya, Aug 2005 
C     
C***********************************************************************

      SUBROUTINE QUAD_PTS_EDGE(P,pad)

C.....Use appropriate modules

      USE GLOBAL
      USE DG

      IMPLICIT NONE

C.....Declare local variables

      INTEGER P,pad

C.....Allocate the gauss points and weights arrays for the edge integral
      
      if (pad.eq.1.and.p.eq.1) then

         NEGP(ph) = ph + 1

         CALL ALLOC_EDGE_GAUSS()

      endif

C     sb moved to prep_DG.F
C     C.....Compute the number of gauss points needed for the edge integrals
C     
      NEGP(pad) = P + 1
      
C.....Retrieve the correct gauss quadrature points for the edge integral
      
      IF (NEGP(pad).EQ.2) THEN

         XEGP(1,pad) = -0.57735026918963D0
         XEGP(2,pad) =  0.57735026918963D0
         
         WEGP(1,pad) = 1.D0
         WEGP(2,pad) = 1.D0
         
      ELSEIF (NEGP(pad).EQ.3) THEN
         
         XEGP(1,pad) = -0.77459666924148D0
         XEGP(2,pad) =  0.D0
         XEGP(3,pad) =  0.77459666924148D0
         
         WEGP(1,pad) =  0.55555555555556D0
         WEGP(2,pad) =  0.88888888888888D0
         WEGP(3,pad) =  0.55555555555556D0
         
      ELSEIF (NEGP(pad).EQ.4) THEN

         XEGP(1,pad) = -0.86113631159405D0
         XEGP(2,pad) = -0.33998104358486D0
         XEGP(3,pad) =  0.33998104358486D0
         XEGP(4,pad) =  0.86113631159405D0
         
         WEGP(1,pad) =  0.34785484513745D0
         WEGP(2,pad) =  0.65214515486255D0
         WEGP(3,pad) =  0.65214515486255D0
         WEGP(4,pad) =  0.34785484513745D0

      ELSEIF (NEGP(pad).EQ.5) THEN

         XEGP(1,pad) = -0.90617984593866D0
         XEGP(2,pad) = -0.53846931010568D0
         XEGP(3,pad) =  0.D0
         XEGP(4,pad) =  0.53846931010568D0
         XEGP(5,pad) =  0.90617984593866D0
         
         WEGP(1,pad) =  0.23692688505619D0
         WEGP(2,pad) =  0.47862867049937D0
         WEGP(3,pad) =  0.56888888888889D0
         WEGP(4,pad) =  0.47862867049937D0
         WEGP(5,pad) =  0.23692688505619D0

      ELSEIF (NEGP(pad).EQ.6) THEN

         XEGP(1,pad) = -0.93246951420315D0
         XEGP(2,pad) = -0.66120938646626D0
         XEGP(3,pad) = -0.23861918608320D0
         XEGP(4,pad) =  0.23861918608320D0
         XEGP(5,pad) =  0.66120938646626D0
         XEGP(6,pad) =  0.93246951420315D0
         
         WEGP(1,pad) =  0.17132449237917D0
         WEGP(2,pad) =  0.36076157304814D0
         WEGP(3,pad) =  0.46791393457269D0
         WEGP(4,pad) =  0.46791393457269D0
         WEGP(5,pad) =  0.36076157304814D0
         WEGP(6,pad) =  0.17132449237917D0

      ELSEIF (NEGP(pad).EQ.7) THEN
         
         XEGP(1,pad) = -0.94910791234276D0
         XEGP(2,pad) = -0.74153118559939D0
         XEGP(3,pad) = -0.40584515137740D0
         XEGP(4,pad) =  0.D0
         XEGP(5,pad) =  0.40584515137740D0
         XEGP(6,pad) =  0.74153118559939D0
         XEGP(7,pad) =  0.94910791234276D0
         
         WEGP(1,pad) =  0.12948496616887D0
         WEGP(2,pad) =  0.27970539148928D0
         WEGP(3,pad) =  0.38183005050512D0
         WEGP(4,pad) =  0.41795918367347D0
         WEGP(5,pad) =  0.38183005050512D0
         WEGP(6,pad) =  0.27970539148928D0
         WEGP(7,pad) =  0.12948496616887D0

      ELSEIF (NEGP(pad).EQ.8) THEN
         
         XEGP(1,pad) = -0.96028985649754D0
         XEGP(2,pad) = -0.79666647741363D0
         XEGP(3,pad) = -0.52553240991633D0
         XEGP(4,pad) = -0.18343464249565D0
         XEGP(5,pad) =  0.18343464249565D0
         XEGP(6,pad) =  0.52553240991633D0
         XEGP(7,pad) =  0.79666647741363D0
         XEGP(8,pad) =  0.96028985649754D0
         
         WEGP(1,pad) =  0.10122853629038D0
         WEGP(2,pad) =  0.22238103445337D0
         WEGP(3,pad) =  0.31370664587789D0
         WEGP(4,pad) =  0.36268378337836D0
         WEGP(5,pad) =  0.36268378337836D0
         WEGP(6,pad) =  0.31370664587789D0
         WEGP(7,pad) =  0.22238103445337D0
         WEGP(8,pad) =  0.10122853629038D0

      ELSEIF (NEGP(pad).EQ.9) THEN
         
         XEGP(1,pad) = -0.96816023950763D0
         XEGP(2,pad) = -0.83603110732664D0
         XEGP(3,pad) = -0.61337143270059D0
         XEGP(4,pad) = -0.32425342340381D0
         XEGP(5,pad) =  0.D0
         XEGP(6,pad) =  0.32425342340381D0
         XEGP(7,pad) =  0.61337143270059D0
         XEGP(8,pad) =  0.83603110732664D0
         XEGP(9,pad) =  0.96816023950763D0
         
         WEGP(1,pad) = 0.08127438836163D0
         WEGP(2,pad) = 0.18064816069483D0
         WEGP(3,pad) = 0.26061069640294D0
         WEGP(4,pad) = 0.31234707704000D0
         WEGP(5,pad) = 0.33023935500126D0
         WEGP(6,pad) = 0.31234707704000D0
         WEGP(7,pad) = 0.26061069640294D0
         WEGP(8,pad) = 0.18064816069483D0
         WEGP(9,pad) = 0.08127438836163D0
         
      ELSEIF (NEGP(pad).EQ.10) THEN
         
         XEGP(1,pad)  = -0.97390652851717D0
         XEGP(2,pad)  = -0.86506336668898D0
         XEGP(3,pad)  = -0.67940956829902D0
         XEGP(4,pad)  = -0.43339539412925D0
         XEGP(5,pad)  = -0.14887433898163D0
         XEGP(6,pad)  =  0.14887433898163D0
         XEGP(7,pad)  =  0.43339539412925D0
         XEGP(8,pad)  =  0.67940956829902D0
         XEGP(9,pad)  =  0.86506336668898D0
         XEGP(10,pad) =  0.97390652851717D0

         WEGP(1,pad)  =  0.06667134430869D0
         WEGP(2,pad)  =  0.14945134915058D0
         WEGP(3,pad)  =  0.21908636251598D0
         WEGP(4,pad)  =  0.26926671931000D0
         WEGP(5,pad)  =  0.29552422471475D0
         WEGP(6,pad)  =  0.29552422471475D0
         WEGP(7,pad)  =  0.26926671931000D0
         WEGP(8,pad)  =  0.21908636251598D0
         WEGP(9,pad)  =  0.14945134915058D0
         WEGP(10,pad) =  0.06667134430869D0

      ELSEIF (NEGP(pad).EQ.11) THEN
         
         XEGP(1,pad)  = -0.97822865814606D0
         XEGP(2,pad)  = -0.88706259976810D0
         XEGP(3,pad)  = -0.73015200557405D0
         XEGP(4,pad)  = -0.51909612920681D0
         XEGP(5,pad)  = -0.26954315595234D0
         XEGP(6,pad)  =  0.D0
         XEGP(7,pad)  =  0.26954315595234D0
         XEGP(8,pad)  =  0.51909612920681D0
         XEGP(9,pad)  =  0.73015200557405D0
         XEGP(10,pad) =  0.88706259976810D0
         XEGP(11,pad) =  0.97822865814606D0
         
         WEGP(1,pad)  =  0.05566856711627D0
         WEGP(2,pad)  =  0.12558036946485D0
         WEGP(3,pad)  =  0.18629021092774D0
         WEGP(4,pad)  =  0.23319376459199D0
         WEGP(5,pad)  =  0.26280454451025D0
         WEGP(6,pad)  =  0.27292508677790D0
         WEGP(7,pad)  =  0.26280454451025D0
         WEGP(8,pad)  =  0.23319376459199D0
         WEGP(9,pad)  =  0.18629021092774D0
         WEGP(10,pad) =  0.12558036946485D0
         WEGP(11,pad) =  0.05566856711627D0
         
      ELSEIF (NEGP(pad).EQ.12) THEN
         
         XEGP(1,pad)  = -0.98156063424672D0
         XEGP(2,pad)  = -0.90411725637047D0
         XEGP(3,pad)  = -0.76990267419430D0
         XEGP(4,pad)  = -0.58731795428662D0
         XEGP(5,pad)  = -0.36783149899818D0
         XEGP(6,pad)  = -0.12523340851147D0
         XEGP(7,pad)  =  0.12523340851147D0
         XEGP(8,pad)  =  0.36783149899818D0
         XEGP(9,pad)  =  0.58731795428662D0
         XEGP(10,pad) =  0.76990267419430D0
         XEGP(11,pad) =  0.90411725637047D0
         XEGP(12,pad) =  0.98156063424672D0
         
         WEGP(1,pad)  =  0.04717533638677D0
         WEGP(2,pad)  =  0.10693932599520D0
         WEGP(3,pad)  =  0.16007832854334D0
         WEGP(4,pad)  =  0.20316742672308D0
         WEGP(5,pad)  =  0.23349253653835D0
         WEGP(6,pad)  =  0.24914704581340D0
         WEGP(7,pad)  =  0.24914704581340D0
         WEGP(8,pad)  =  0.23349253653835D0
         WEGP(9,pad)  =  0.20316742672308D0
         WEGP(10,pad) =  0.16007832854334D0
         WEGP(11,pad) =  0.10693932599520D0
         WEGP(12,pad) =  0.04717533638677D0

      ELSEIF (NEGP(pad).EQ.13) THEN
         
         XEGP(1,pad)  = -0.98418305471859D0
         XEGP(2,pad)  = -0.91759839922298D0
         XEGP(3,pad)  = -0.80157809073331D0
         XEGP(4,pad)  = -0.64234933944034D0
         XEGP(5,pad)  = -0.44849275103645D0
         XEGP(6,pad)  = -0.23045831595513D0
         XEGP(7,pad)  =  0.D0
         XEGP(8,pad)  =  0.23045831595513D0
         XEGP(9,pad)  =  0.44849275103645D0
         XEGP(10,pad) =  0.64234933944034D0
         XEGP(11,pad) =  0.80157809073331D0
         XEGP(12,pad) =  0.91759839922298D0
         XEGP(13,pad) =  0.98418305471859D0

         WEGP(1,pad)  =  0.04048400476532D0
         WEGP(2,pad)  =  0.09212149983773D0
         WEGP(3,pad)  =  0.13887351021979D0
         WEGP(4,pad)  =  0.17814598076195D0
         WEGP(5,pad)  =  0.20781604753689D0
         WEGP(6,pad)  =  0.22628318026290D0
         WEGP(7,pad)  =  0.23255155323087D0
         WEGP(8,pad)  =  0.22628318026290D0
         WEGP(9,pad)  =  0.20781604753689D0
         WEGP(10,pad) =  0.17814598076195D0
         WEGP(11,pad) =  0.13887351021979D0
         WEGP(12,pad) =  0.09212149983773D0
         WEGP(13,pad) =  0.04048400476532D0

      ENDIF

      RETURN
      END SUBROUTINE
