C***********************************************************************
C     
C     MODULE CONT_INT_TEST
C     
C     Module contains the functions required to run test_suite
C     
C     Written by Max Bremer
C     
C***********************************************************************
C
      MODULE CONT_INT_TEST
C
      USE SIZES, ONLY : SZ

      PUBLIC
      REAL(SZ),allocatable :: XBCbt(:),YBCbt(:),radial(:)
      REAL(SZ) :: s_sed
#ifdef WETDR 
C.....Set wetting and drying parameters for parabolic bowl
      REAL(SZ),parameter :: Xpb = 1.d0 !X
      REAL(SZ),parameter :: Ypb = -0.41884d0 !Y
      REAL(SZ),parameter :: alphapb = 1.6D-7
      REAL(SZ) :: omegapb = sqrt( 80.d0 * alphapb ) !g=10 is hardcoded here
#endif
      CONTAINS

      SUBROUTINE INTIALIZE_MAN_SOLN(hbo,hb1,ydubo)
C--------------------------------------------------------------------------
C  Routine sets up intial conditions for manufactured solutions
C--------------------------------------------------------------------------
C.....Use appropriate modules

      USE GLOBAL
      USE DG
      USE SIZES

#ifdef CMPI
      USE MESSENGER_ELEM 
#endif       

      IMPLICIT NONE
      REAL(SZ), allocatable, INTENT(INOUT) :: HB1(:,:,:,:)
      REAL(SZ), allocatable, INTENT(IN) :: ydubo(:,:),hbo(:,:,:)
      REAL(SZ), allocatable :: local_inst(:)

C.....Declare local variables
      REAL(SZ) pi_n,rev,Ox,Oy,epsl
      REAL(SZ) int_hb,int_ze,int_yd
      REAL(SZ) C_0,C_1,ell_1,ell_2,ell_3

      INTEGER k,mm,j,chi


      Allocate ( XBCbt(MNE),YBCbt(MNE),radial(MNE) )

C.....Compute initial conditions
      pi_n = 3.1415926535897932384626433832795
      rev = pi_n/4.d0
      
#ifdef SED_LAY
      bed = 0.D0
      epsl = 0.0000001D0
#endif
      
      do k = 1,MNE

         n1 = NM(k,1)
         n2 = NM(k,2)
         n3 = NM(k,3)
                                !Define lagrange transform
#ifdef SED_LAY
         int_hb = 0.D0
         int_ze = 0.D0
         int_yd = 0.D0

         C_0 = 0.000001D0
#endif

         do mm = 1,nagp(pdg_el(k)) !ICs should not have higher order than ph
            

            ell_1 = -0.5D0 * ( xagp(mm,pdg_el(k)) + yagp(mm,pdg_el(k)) )
            ell_2 =  0.5D0 * ( xagp(mm,pdg_el(k)) + 1.D0 )
            ell_3 =  0.5D0 * ( yagp(mm,pdg_el(k)) + 1.D0 )

            XBCbt(k) = x(n1)*ell_1 + x(n2)*ell_2 + x(n3)*ell_3
            YBCbt(k) = y(n1)*ell_1 + y(n2)*ell_2 + y(n3)*ell_3

 
#ifdef SED_LAY
            do j = 1,dofh
C.....Initialize initial conditions for sediment test case
               QX(j,k,1) = QX(j,k,1) +   0.001D0*YBCbt(k) * wagp(mm,ph) * phi_area(j,mm,ph)
     &              *0.5D0 *AREAS(K)
               QY(j,k,1) = QY(j,k,1) -   0.001D0*XBCbt(k) * wagp(mm,ph) * phi_area(j,mm,ph)
     &              *0.5D0 *AREAS(K)

c$$$               hb(j,k,1) = hb(j,k,1) + XBCbt(k) * wagp(mm,ph) * phi_area(j,mm,ph)
c$$$     &              *0.5D0 *AREAS(K)

               hb(j,k,1) = hb(j,k,1) + (1.0D0 + 0.2D0
     &              *exp(-( ( XBCbt(k))**2 / 0.01D0 + ( YBCbt(k))**2 / 0.01D0 ) )) * wagp(mm,ph) * phi_area(j,mm,ph)
     &              *0.5D0 *AREAS(K)

               ze(j,k,1) = ze(j,k,1) + (2.D0 - (1.0D0 + 0.2D0
     &              *exp(-( ( XBCbt(k))**2 / 0.01D0 + ( YBCbt(k))**2 / 0.01D0 ) )) ) * wagp(mm,ph) * phi_area(j,mm,ph)
     &              *0.5D0 *AREAS(K)
            enddo
#elif TRACE
            Ox = XBCbt(k)*cos(rev) - YBCbt(k)*sin(rev)
            Oy = YBCbt(k)*cos(rev) + XBCbt(k)*sin(rev)

            radial(k) = min(sqrt( Ox**2 + (Oy + 0.25D0)**2 ), 0.18D0)/0.18D0
            do j = 1,dofh
C.....Initialize test conditions for rotating blob test case
               QX(j,k,1) = QX(j,k,1) + YBCbt(k) * wagp(mm,ph) * phi_area(j,mm,ph)

               QY(j,k,1) = QY(j,k,1) - XBCbt(k) * wagp(mm,ph) * phi_area(j,mm,ph)   
               iota(j,k,1) = iota(j,k,1) + 0.5D0 *( exp ( - ( (XBCbt(k)+.05D0)**2 + (YBCbt(k)+.05D0)**2  ) /0.001D0 ) ) 
     &              * wagp(mm,ph) * phi_area(j,mm,ph) 

               hb(j,k,1) = hb(j,k,1) + 1.D0*wagp(mm,ph) * phi_area(j,mm,ph)
            enddo

#endif

         enddo

#ifdef SED_LAY
         do chi=1,ph            ! Need every level data for functions independent of time

            do mm = 1,nagp(chi)            

               do j=1,(chi+1)*(chi+2)/2

                  ydub(j,k,chi) = ydub(j,k,chi) + ydubo(j,k) * wagp(mm,chi) * phi_area(j,mm,chi) * 0.5D0 * AREAS(K)
                  hb1(j,k,1,chi) = hb1(j,k,1,chi) + hbo(j,k,1) * wagp(mm,chi) * phi_area(j,mm,chi)* 0.5D0 * AREAS(K)

               enddo

            enddo
            
         enddo
                                ! get back the coeffs in each component

         do j= 1,dofs(k)

            qx(j,k,1) =  qx(j,k,1) * M_inv(j,pdg_el(k))/ (0.5D0 *AREAS(K))
            qy(j,k,1) =  qy(j,k,1) * M_inv(j,pdg_el(k))/ (0.5D0 *AREAS(K))
            hb(j,k,1) =  hb(j,k,1) * M_inv(j,pdg_el(k))/ (0.5D0 *AREAS(K))
            ze(j,k,1) =  ze(j,k,1) * M_inv(j,pdg_el(k))/ (0.5D0 *AREAS(K))

         enddo

         do chi = 1,ph
            do j=1,(chi+1)*(chi+2)/2              
  !These must agree with the active degree of the computation in p-enrichment
               ydub(j,k,chi)  =  ydub(j,k,chi) * M_inv(j,chi)/ (0.5D0 *AREAS(K))
               hb1(j,k,1,chi) =  hb1(j,k,1,chi) * M_inv(j,chi)/ (0.5D0 *AREAS(K))

            enddo
         enddo
#elif TRACE

      do j= 1,dofh

            QX(j,k,1) =  QX(j,k,1) * M_inv(j,pdg_el(k))
            QY(j,k,1) =  QY(j,k,1) * M_inv(j,pdg_el(k))
            iota(j,k,1) =  iota(j,k,1) * M_inv(j,pdg_el(k))
            hb(j,k,1) =  hb(j,k,1) * M_inv(j,pdg_el(k))

      enddo
#endif


      enddo !Close k

 
#ifdef SED_LAY
      iotaa2 = ze
      iotaa3 = qx
#elif TRACE
      iotaa = iota
#endif

      END SUBROUTINE

      SUBROUTINE RHS_MAN_SED_LAY(i,l)
C--------------------------------------------------------------------------
C  Routine constructs the right hand side for the manufactured solution
C  for the sediment test case
C--------------------------------------------------------------------------
C.....Use appropriate modules
      USE GLOBAL
      USE DG     

      IMPLICIT NONE

C.....Declare local variables
      REAL(SZ) s_mass,fx,fy
      REAL(SZ) fgauss,sig,b_0,H_0
      REAL(SZ) C_0,C_1,ell_1,ell_2,ell_3

      INTEGER i,l

      s_sed = 0.D0
      s_mass = 0.D0

      XBCbt(l) =0.D0
      YBCbt(l)  = 0.D0

      C_0 = 0.1D0
      

      ell_1 = -0.5D0 * ( xagp(i,pdg_el(l)) + yagp(i,pdg_el(l)) )
      ell_2 =  0.5D0 * ( xagp(i,pdg_el(l)) + 1.D0 )
      ell_3 =  0.5D0 * ( yagp(i,pdg_el(l)) + 1.D0 )

      XBCbt(l) = x(n1)*ell_1 + x(n2)*ell_2 + x(n3)*ell_3
      YBCbt(l) = y(n1)*ell_1 + y(n2)*ell_2 + y(n3)*ell_3 


      sig = 0.01D0
      fx = exp(-( ( XBCbt(l))**2 / sig ) )
      fy = exp(-( ( YBCbt(l))**2 / sig ) )
      fgauss = exp(-( ( XBCbt(l))**2 / sig + ( YBCbt(l))**2 / sig ) )

      H_0 = 2.0D0
      C_1 = 0.2D0
      b_0 = 1.0D0

      s_sed = -sin(timedg)*C_1*fgauss

      RHS_ZE(1,L,IRK)    = RHS_ZE(1,L,IRK)    + SRFAC(1,I,L,pa)*s_mass
      RHS_bed(1,L,IRK,1) = RHS_bed(1,L,IRK,1) + SRFAC(1,I,L,pa)*s_sed

      END SUBROUTINE RHS_MAN_SED_LAY

      SUBROUTINE TRUE_MAN_SOLN()
C--------------------------------------------------------------------------
C  Routine computes true state for manufactured solution at given point
C  in time
C--------------------------------------------------------------------------
C.....Use appropriate modules
      USE GLOBAL
      USE DG     

      IMPLICIT NONE

C.....Declare local variables
      REAL(SZ) pi_n,rev,Ox,Oy
      REAL(SZ) C_0,C_1,ell_1,ell_2,ell_3
      REAL(SZ) tempx,tempy,rex2

      INTEGER k,mm,j

#ifdef WETDR 
      rex2 = (Xpb + Ypb * cos(omegapb * time_a) ) / (alphapb * (Xpb**2 - Ypb**2))
#endif

      do k = 1,MNE

         iotaa(:,k,1) = 0.D0
#if defined (SED_LAY) || defined (WETDR)
         iotaa2(:,k,1) = 0.D0
         iotaa3(:,k,1) = 0.D0

         C_0 = 0.000001D0
#endif

         pi_n = 3.1415926535897932384626433832795
         
         n1 = NM(k,1)
         n2 = NM(k,2)
         n3 = NM(k,3)
         

                                !Define lagrange transform
         
         do mm = 1,nagp(ph)     !ICs should not have higher order than ph

            ell_1 = -0.5D0 * ( xagp(mm,ph) + yagp(mm,ph) )
            ell_2 =  0.5D0 * ( xagp(mm,ph) + 1.D0 )
            ell_3 =  0.5D0 * ( yagp(mm,ph) + 1.D0 )

            XBCbt(k) = x(n1)*ell_1 + x(n2)*ell_2 + x(n3)*ell_3
            YBCbt(k) = y(n1)*ell_1 + y(n2)*ell_2 + y(n3)*ell_3
  
#ifdef SED_LAY
            do j = 1,dofh
               

               iotaa(j,k,1) = iotaa(j,k,1) +  (1.0D0 + ( 0.2D0
     &              *exp(-( ( XBCbt(k))**2 / 0.01D0 + ( YBCbt(k))**2 / 0.01D0 ) ) * cos(time_a) ) )
     &              * wagp(mm,ph) * phi_area(j,mm,ph)
     &              *0.5D0 *AREAS(K)

               iotaa2(j,k,1) = iotaa2(j,k,1) +  (2.D0 - (1.0D0 + ( 0.2D0
     &              *exp(-( ( XBCbt(k))**2 / 0.01D0 + ( YBCbt(k))**2 / 0.01D0 ) ) * cos(time_a) ) ))
     &              * wagp(mm,ph) * phi_area(j,mm,ph)
     &              *0.5D0 *AREAS(K)


               iotaa3(j,k,1) = iotaa3(j,k,1) +  0.001D0* YBCbt(k)
     &              * wagp(mm,ph) * phi_area(j,mm,ph)
     &              *0.5D0 *AREAS(K)
               
            enddo
#endif

#ifdef TRACE
            rev =  3.141592653589793D0 / 4.D0

            Ox = XBCbt(k)*cos(rev) - YBCbt(k)*sin(rev)
            Oy = YBCbt(k)*cos(rev) + XBCbt(k)*sin(rev)

            tempx = XBCbt(k)*cos(time_a) - YBCbt(k)*sin(time_a)
            tempy = YBCbt(k)*cos(time_a) + XBCbt(k)*sin(time_a)

            radial(k) = min(sqrt( tempx**2 + (tempy + 0.25D0)**2 ), 0.18D0)/0.18D0

          do j = 1,dofh
             iotaa(j,k,1) = iotaa(j,k,1) + 0.5D0 *( exp ( - ( (tempx+.05D0)**2 + (tempy+.05D0)**2  )
     &            /0.001D0 ) ) * wagp(mm,ph) * phi_area(j,mm,ph) *  (0.5D0 *AREAS(K))
            enddo    
#endif

#ifdef WETDR
            if (XBCbt(k)**2 + YBCbt(k)**2 .LE. rex2) then
            do j = 1,dofh
C.....Compute ze
               iotaa(j,k,1) = iotaa(j,k,1) +  (1.d0/(Xpb + Ypb * cos(omegapb* time_a))
     &              + alphapb *(Ypb**2 - Xpb**2)* (XBCbt(k)**2 + YBCbt(k)**2)
     &              /(Xpb + Ypb* cos(omegapb*time_a)) - alphapb * (XBCbt(k)**2 + YBCbt(k)**2))  
     &              * wagp(mm,ph) * phi_area(j,mm,ph)
     &              *0.5D0 *AREAS(K)
C.....Compute uH
               iotaa2(j,k,1) = iotaa2(j,k,1) +  (-Ypb * omegapb * sin(omegapb *time_a)
     &              /(Xpb + Ypb * cos( omegapb * time_a)) * XBCbt(k) * 0.5d0 )
     &              * (1.d0/(Xpb + Ypb * cos(omegapb* time_a))                                                                                                                
     &              + alphapb *(Ypb**2 - Xpb**2)* (XBCbt(k)**2 + YBCbt(k)**2)                
     &              /(Xpb + Ypb* cos(omegapb*time_a)))
     &              * wagp(mm,ph) * phi_area(j,mm,ph)
     &              *0.5D0 *AREAS(K)
C.....Compute vH
               iotaa3(j,k,1) = iotaa3(j,k,1) +  (-Ypb * omegapb * sin(omegapb *time_a)
     &              /(Xpb + Ypb * cos( omegapb * time_a)) * YBCbt(k) * 0.5d0 )
     &              * (1.d0/(Xpb + Ypb * cos(omegapb* time_a))
     &              + alphapb *(Ypb**2 - Xpb**2)* (XBCbt(k)**2 + YBCbt(k)**2)
     &              /(Xpb + Ypb* cos(omegapb*time_a)))
     &              * wagp(mm,ph) * phi_area(j,mm,ph)
     &              *0.5D0 *AREAS(K)
            enddo
            else
               do j = 1,dofh
C.....Compute ze                                                                                                                                                                                     
               iotaa(j,k,1) = iotaa(j,k,1) +  H0                                                                                                       
     &              * wagp(mm,ph) * phi_area(j,mm,ph)
     &              *0.5D0 *AREAS(K)
C.....Compute uH                                                                                                                                                  
               iotaa2(j,k,1) = 0.d0
C.....Compute vH                                                                                                                                                                                    
               iotaa3(j,k,1) = 0.d0
               enddo
            endif

#endif

#ifdef SED_LAY
            do j = 1,dofh


               iotaa(j,k,1) = iotaa(j,k,1) +  (1.0D0 + ( 0.2D0
     &              *exp(-( ( XBCbt(k))**2 / 0.01D0 + ( YBCbt(k))**2 / 0.01D0 ) ) * cos(time_a) \
) )                                                                                              
     &              * wagp(mm,ph) * phi_area(j,mm,ph)
     &              *0.5D0 *AREAS(K)

               iotaa2(j,k,1) = iotaa2(j,k,1) +  (2.D0 - (1.0D0 + ( 0.2D0
     &              *exp(-( ( XBCbt(k))**2 / 0.01D0 + ( YBCbt(k))**2 / 0.01D0 ) ) * cos(time_a) \
) ))                                                                                             
     &              * wagp(mm,ph) * phi_area(j,mm,ph)
     &              *0.5D0 *AREAS(K)


               iotaa3(j,k,1) = iotaa3(j,k,1) +  0.001D0* YBCbt(k)
     &              * wagp(mm,ph) * phi_area(j,mm,ph)
     &              *0.5D0 *AREAS(K)

            enddo
#endif

         enddo

                                ! get back the coeffs in each component

         do j= 1,dofh

            iotaa(j,k,1)  =  iotaa(j,k,1) * M_inv(j,ph) / (0.5D0 *AREAS(K))
#if defined (SED_LAY) || defined (WETDR)
            iotaa2(j,k,1) =  iotaa2(j,k,1) * M_inv(j,ph) / (0.5D0 *AREAS(K))
            iotaa3(j,k,1) =  iotaa3(j,k,1) * M_inv(j,ph) / (0.5D0 *AREAS(K))
#endif

         enddo

      enddo
      END SUBROUTINE TRUE_MAN_SOLN

      SUBROUTINE COMPUTE_ERROR(IT)
C--------------------------------------------------------------------------
C  Compute Error of numerical solution for manufactured solution
C--------------------------------------------------------------------------
C.....Use appropriate modules
      USE GLOBAL
      USE DG

#ifdef CMPI
      USE MESSENGER_ELEM 
#endif     

      IMPLICIT NONE

C.....Declare local variables
      Real(sz) iota_error,iota_int,iota2_int,iota3_int,iota2_error,iota3_error,testsum
      Real(SZ), allocatable :: temper(:)

      Integer i,j,k,it

      Allocate ( temper(mne) )

      iota_int = 0.D0

      do j = 1,ne         
         DO I = 1,NAGP(pa)

            iota_error = 0.D0           
            do k = 1,dofs(j)
#ifdef SED_LAY
               iota_error = iota_error + (bed(k,j,1,1)- iotaa(k,j,1)) * PHI_AREA(k,I,pa)*AREAS(J)*0.5D0
#elif TRACE
               iota_error = iota_error + (iota(k,j,1)- iotaa(k,j,1)) * PHI_AREA(k,I,pa)
#elif WETDR
               iota_error = iota_error + (ze(k,j,1) - iotaa(k,j,1)) * PHI_AREA(k,I,pa)
#endif
            enddo

            iota_int = iota_int + iota_error**2 * WAGP(I,pa)*AREAS(J)

         ENDDO
      enddo

#ifdef WETDR
      iota2_int = 0.D0
      iota3_int = 0.D0

      do j = 1,ne
         DO I = 1,NAGP(pa)
            iota2_error = 0.D0
            iota3_error = 0.D0
            do k = 1,dofs(j)
               iota2_error = iota2_error + (qx(k,J,1)- iotaa2(k,j,1)) * PHI_AREA(k,I,pa)*AREAS(J)*0.5D0                                                                                               
               iota3_error = iota3_error + (qy(k,J,1)- iotaa3(k,j,1)) * PHI_AREA(k,I,pa)*AREAS(J)*0.5D0                                                                                               
            enddo

            iota2_int = iota2_int + iota2_error**2 * WAGP(I,pa)*AREAS(J)
            iota3_int = iota3_int + iota3_error**2 * WAGP(I,pa)*AREAS(J)

         ENDDO

      enddo

      temper =  0.d0 !abs(dpe - tracer2)
#endif

#ifdef SED_LAY
      iota2_int = 0.D0
      iota3_int = 0.D0

      do j = 1,ne
         DO I = 1,NAGP(pa)
            iota2_error = 0.D0
            iota3_error = 0.D0
            do k = 1,dofs(j)
               iota2_error = iota2_error + (ze(k,J,1)- iotaa2(k,j,1)) * PHI_AREA(k,I,pa)*AREAS(J)*0.5D0
               iota3_error = iota3_error + (qx(K,J,1)- iotaa3(k,j,1)) * PHI_AREA(k,I,pa)*AREAS(J)*0.5D0
            enddo
 
            iota2_int = iota2_int + iota2_error**2 * WAGP(I,pa)*AREAS(J)
            iota3_int = iota3_int + iota3_error**2 * WAGP(I,pa)*AREAS(J)

         ENDDO
         
      enddo

      temper =  abs(dpe - tracer2)
#endif
   

#ifdef TRACE
      temper = abs(tracer - tracer2)
#endif

      temperg = maxval(temper)
         
C.....Write errors into console
#ifdef CMPI
      Call para_max(temperg)
      Call para_sum(iota_int)
#endif
      IF (MYPROC .EQ. 0) THEN
         write(*,*) 'L2nerror('') = ',  sqrt( iota_int ) , it 
         write(*,*) 'linferror('') = ',  temperg , it
      ENDIF
#if defined (SED_LAY) || defined (WETDR)
#ifdef CMPI
      Call para_sum(iota2_int)
      Call para_sum(iota3_int)
#endif
      IF (MYPROC .EQ. 0) THEN
         write(*,*) 'L2nerror2('') = ',  sqrt( iota2_int ) , it 
         write(*,*) 'L2nerror3('') = ',  sqrt( iota3_int ) , it 
      ENDIF
#endif

      END SUBROUTINE COMPUTE_ERROR
      END MODULE CONT_INT_TEST
