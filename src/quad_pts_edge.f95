!***********************************************************************
!     
!     SUBROUTINE QUAD_PTS_EDGE()
!     
!     Stores the Gauss quadrature points and weights for up to 13 point
!     quadrature rules on the interval (-1,1) (an n point Legendre-Gauss
!     quadrature integrates a 2*n-1 degree polynomial exactly)
!     
!     Written by Ethan Kubatko (06-08-2004)
!     
!     Modifications for parallel runs, Shintaro Bunya, Aug 2005 
!     
!***********************************************************************

      SUBROUTINE QUAD_PTS_EDGE(s,dg,P,pad)

!.....Use appropriate modules

      USE GLOBAL
      USE DG
      use sizes

      IMPLICIT NONE
      
      type (sizes_type) :: s
      type (dg_type) :: dg

!.....Declare local variables

      INTEGER P,pad

!.....Allocate the gauss points and weights arrays for the edge integral
      
      if (pad.eq.1.and.p.eq.1) then

         dg%NEGP(dg%ph) = dg%ph + 1

         CALL ALLOC_EDGE_GAUSS(s)

      endif

!     sb moved to prep_DG.dg%F
!     C.....Compute the number of gauss points needed for the edge integrals
!     
      dg%NEGP(pad) = P + 1
      
!.....Retrieve the correct gauss quadrature points for the edge integral
      
      IF (dg%NEGP(pad).EQ.2) THEN

         dg%XEGP(1,pad) = -0.57735026918963D0
         dg%XEGP(2,pad) =  0.57735026918963D0
         
         dg%WEGP(1,pad) = 1.D0
         dg%WEGP(2,pad) = 1.D0
         
      ELSEIF (dg%NEGP(pad).EQ.3) THEN
         
         dg%XEGP(1,pad) = -0.77459666924148D0
         dg%XEGP(2,pad) =  0.D0
         dg%XEGP(3,pad) =  0.77459666924148D0
         
         dg%WEGP(1,pad) =  0.55555555555556D0
         dg%WEGP(2,pad) =  0.88888888888888D0
         dg%WEGP(3,pad) =  0.55555555555556D0
         
      ELSEIF (dg%NEGP(pad).EQ.4) THEN

         dg%XEGP(1,pad) = -0.86113631159405D0
         dg%XEGP(2,pad) = -0.33998104358486D0
         dg%XEGP(3,pad) =  0.33998104358486D0
         dg%XEGP(4,pad) =  0.86113631159405D0
         
         dg%WEGP(1,pad) =  0.34785484513745D0
         dg%WEGP(2,pad) =  0.65214515486255D0
         dg%WEGP(3,pad) =  0.65214515486255D0
         dg%WEGP(4,pad) =  0.34785484513745D0

      ELSEIF (dg%NEGP(pad).EQ.5) THEN

         dg%XEGP(1,pad) = -0.90617984593866D0
         dg%XEGP(2,pad) = -0.53846931010568D0
         dg%XEGP(3,pad) =  0.D0
         dg%XEGP(4,pad) =  0.53846931010568D0
         dg%XEGP(5,pad) =  0.90617984593866D0
         
         dg%WEGP(1,pad) =  0.23692688505619D0
         dg%WEGP(2,pad) =  0.47862867049937D0
         dg%WEGP(3,pad) =  0.56888888888889D0
         dg%WEGP(4,pad) =  0.47862867049937D0
         dg%WEGP(5,pad) =  0.23692688505619D0

      ELSEIF (dg%NEGP(pad).EQ.6) THEN

         dg%XEGP(1,pad) = -0.93246951420315D0
         dg%XEGP(2,pad) = -0.66120938646626D0
         dg%XEGP(3,pad) = -0.23861918608320D0
         dg%XEGP(4,pad) =  0.23861918608320D0
         dg%XEGP(5,pad) =  0.66120938646626D0
         dg%XEGP(6,pad) =  0.93246951420315D0
         
         dg%WEGP(1,pad) =  0.17132449237917D0
         dg%WEGP(2,pad) =  0.36076157304814D0
         dg%WEGP(3,pad) =  0.46791393457269D0
         dg%WEGP(4,pad) =  0.46791393457269D0
         dg%WEGP(5,pad) =  0.36076157304814D0
         dg%WEGP(6,pad) =  0.17132449237917D0

      ELSEIF (dg%NEGP(pad).EQ.7) THEN
         
         dg%XEGP(1,pad) = -0.94910791234276D0
         dg%XEGP(2,pad) = -0.74153118559939D0
         dg%XEGP(3,pad) = -0.40584515137740D0
         dg%XEGP(4,pad) =  0.D0
         dg%XEGP(5,pad) =  0.40584515137740D0
         dg%XEGP(6,pad) =  0.74153118559939D0
         dg%XEGP(7,pad) =  0.94910791234276D0
         
         dg%WEGP(1,pad) =  0.12948496616887D0
         dg%WEGP(2,pad) =  0.27970539148928D0
         dg%WEGP(3,pad) =  0.38183005050512D0
         dg%WEGP(4,pad) =  0.41795918367347D0
         dg%WEGP(5,pad) =  0.38183005050512D0
         dg%WEGP(6,pad) =  0.27970539148928D0
         dg%WEGP(7,pad) =  0.12948496616887D0

      ELSEIF (dg%NEGP(pad).EQ.8) THEN
         
         dg%XEGP(1,pad) = -0.96028985649754D0
         dg%XEGP(2,pad) = -0.79666647741363D0
         dg%XEGP(3,pad) = -0.52553240991633D0
         dg%XEGP(4,pad) = -0.18343464249565D0
         dg%XEGP(5,pad) =  0.18343464249565D0
         dg%XEGP(6,pad) =  0.52553240991633D0
         dg%XEGP(7,pad) =  0.79666647741363D0
         dg%XEGP(8,pad) =  0.96028985649754D0
         
         dg%WEGP(1,pad) =  0.10122853629038D0
         dg%WEGP(2,pad) =  0.22238103445337D0
         dg%WEGP(3,pad) =  0.31370664587789D0
         dg%WEGP(4,pad) =  0.36268378337836D0
         dg%WEGP(5,pad) =  0.36268378337836D0
         dg%WEGP(6,pad) =  0.31370664587789D0
         dg%WEGP(7,pad) =  0.22238103445337D0
         dg%WEGP(8,pad) =  0.10122853629038D0

      ELSEIF (dg%NEGP(pad).EQ.9) THEN
         
         dg%XEGP(1,pad) = -0.96816023950763D0
         dg%XEGP(2,pad) = -0.83603110732664D0
         dg%XEGP(3,pad) = -0.61337143270059D0
         dg%XEGP(4,pad) = -0.32425342340381D0
         dg%XEGP(5,pad) =  0.D0
         dg%XEGP(6,pad) =  0.32425342340381D0
         dg%XEGP(7,pad) =  0.61337143270059D0
         dg%XEGP(8,pad) =  0.83603110732664D0
         dg%XEGP(9,pad) =  0.96816023950763D0
         
         dg%WEGP(1,pad) = 0.08127438836163D0
         dg%WEGP(2,pad) = 0.18064816069483D0
         dg%WEGP(3,pad) = 0.26061069640294D0
         dg%WEGP(4,pad) = 0.31234707704000D0
         dg%WEGP(5,pad) = 0.33023935500126D0
         dg%WEGP(6,pad) = 0.31234707704000D0
         dg%WEGP(7,pad) = 0.26061069640294D0
         dg%WEGP(8,pad) = 0.18064816069483D0
         dg%WEGP(9,pad) = 0.08127438836163D0
         
      ELSEIF (dg%NEGP(pad).EQ.10) THEN
         
         dg%XEGP(1,pad)  = -0.97390652851717D0
         dg%XEGP(2,pad)  = -0.86506336668898D0
         dg%XEGP(3,pad)  = -0.67940956829902D0
         dg%XEGP(4,pad)  = -0.43339539412925D0
         dg%XEGP(5,pad)  = -0.14887433898163D0
         dg%XEGP(6,pad)  =  0.14887433898163D0
         dg%XEGP(7,pad)  =  0.43339539412925D0
         dg%XEGP(8,pad)  =  0.67940956829902D0
         dg%XEGP(9,pad)  =  0.86506336668898D0
         dg%XEGP(10,pad) =  0.97390652851717D0

         dg%WEGP(1,pad)  =  0.06667134430869D0
         dg%WEGP(2,pad)  =  0.14945134915058D0
         dg%WEGP(3,pad)  =  0.21908636251598D0
         dg%WEGP(4,pad)  =  0.26926671931000D0
         dg%WEGP(5,pad)  =  0.29552422471475D0
         dg%WEGP(6,pad)  =  0.29552422471475D0
         dg%WEGP(7,pad)  =  0.26926671931000D0
         dg%WEGP(8,pad)  =  0.21908636251598D0
         dg%WEGP(9,pad)  =  0.14945134915058D0
         dg%WEGP(10,pad) =  0.06667134430869D0

      ELSEIF (dg%NEGP(pad).EQ.11) THEN
         
         dg%XEGP(1,pad)  = -0.97822865814606D0
         dg%XEGP(2,pad)  = -0.88706259976810D0
         dg%XEGP(3,pad)  = -0.73015200557405D0
         dg%XEGP(4,pad)  = -0.51909612920681D0
         dg%XEGP(5,pad)  = -0.26954315595234D0
         dg%XEGP(6,pad)  =  0.D0
         dg%XEGP(7,pad)  =  0.26954315595234D0
         dg%XEGP(8,pad)  =  0.51909612920681D0
         dg%XEGP(9,pad)  =  0.73015200557405D0
         dg%XEGP(10,pad) =  0.88706259976810D0
         dg%XEGP(11,pad) =  0.97822865814606D0
         
         dg%WEGP(1,pad)  =  0.05566856711627D0
         dg%WEGP(2,pad)  =  0.12558036946485D0
         dg%WEGP(3,pad)  =  0.18629021092774D0
         dg%WEGP(4,pad)  =  0.23319376459199D0
         dg%WEGP(5,pad)  =  0.26280454451025D0
         dg%WEGP(6,pad)  =  0.27292508677790D0
         dg%WEGP(7,pad)  =  0.26280454451025D0
         dg%WEGP(8,pad)  =  0.23319376459199D0
         dg%WEGP(9,pad)  =  0.18629021092774D0
         dg%WEGP(10,pad) =  0.12558036946485D0
         dg%WEGP(11,pad) =  0.05566856711627D0
         
      ELSEIF (dg%NEGP(pad).EQ.12) THEN
         
         dg%XEGP(1,pad)  = -0.98156063424672D0
         dg%XEGP(2,pad)  = -0.90411725637047D0
         dg%XEGP(3,pad)  = -0.76990267419430D0
         dg%XEGP(4,pad)  = -0.58731795428662D0
         dg%XEGP(5,pad)  = -0.36783149899818D0
         dg%XEGP(6,pad)  = -0.12523340851147D0
         dg%XEGP(7,pad)  =  0.12523340851147D0
         dg%XEGP(8,pad)  =  0.36783149899818D0
         dg%XEGP(9,pad)  =  0.58731795428662D0
         dg%XEGP(10,pad) =  0.76990267419430D0
         dg%XEGP(11,pad) =  0.90411725637047D0
         dg%XEGP(12,pad) =  0.98156063424672D0
         
         dg%WEGP(1,pad)  =  0.04717533638677D0
         dg%WEGP(2,pad)  =  0.10693932599520D0
         dg%WEGP(3,pad)  =  0.16007832854334D0
         dg%WEGP(4,pad)  =  0.20316742672308D0
         dg%WEGP(5,pad)  =  0.23349253653835D0
         dg%WEGP(6,pad)  =  0.24914704581340D0
         dg%WEGP(7,pad)  =  0.24914704581340D0
         dg%WEGP(8,pad)  =  0.23349253653835D0
         dg%WEGP(9,pad)  =  0.20316742672308D0
         dg%WEGP(10,pad) =  0.16007832854334D0
         dg%WEGP(11,pad) =  0.10693932599520D0
         dg%WEGP(12,pad) =  0.04717533638677D0

      ELSEIF (dg%NEGP(pad).EQ.13) THEN
         
         dg%XEGP(1,pad)  = -0.98418305471859D0
         dg%XEGP(2,pad)  = -0.91759839922298D0
         dg%XEGP(3,pad)  = -0.80157809073331D0
         dg%XEGP(4,pad)  = -0.64234933944034D0
         dg%XEGP(5,pad)  = -0.44849275103645D0
         dg%XEGP(6,pad)  = -0.23045831595513D0
         dg%XEGP(7,pad)  =  0.D0
         dg%XEGP(8,pad)  =  0.23045831595513D0
         dg%XEGP(9,pad)  =  0.44849275103645D0
         dg%XEGP(10,pad) =  0.64234933944034D0
         dg%XEGP(11,pad) =  0.80157809073331D0
         dg%XEGP(12,pad) =  0.91759839922298D0
         dg%XEGP(13,pad) =  0.98418305471859D0

         dg%WEGP(1,pad)  =  0.04048400476532D0
         dg%WEGP(2,pad)  =  0.09212149983773D0
         dg%WEGP(3,pad)  =  0.13887351021979D0
         dg%WEGP(4,pad)  =  0.17814598076195D0
         dg%WEGP(5,pad)  =  0.20781604753689D0
         dg%WEGP(6,pad)  =  0.22628318026290D0
         dg%WEGP(7,pad)  =  0.23255155323087D0
         dg%WEGP(8,pad)  =  0.22628318026290D0
         dg%WEGP(9,pad)  =  0.20781604753689D0
         dg%WEGP(10,pad) =  0.17814598076195D0
         dg%WEGP(11,pad) =  0.13887351021979D0
         dg%WEGP(12,pad) =  0.09212149983773D0
         dg%WEGP(13,pad) =  0.04048400476532D0

      ENDIF

      RETURN
      END SUBROUTINE
