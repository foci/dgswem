      SUBROUTINE cfl_prep(ne,no,ics,X,Y,Hb,EC,SLAM0,SFEA0,dt_local,H0,slopeflag)
      USE PRE_GLOBAL, ONLY: CENTROID_X, CENTROID_Y, cflfac, Z0
      IMPLICIT NONE

      REAL(8), PARAMETER :: G=9.80665D0,PI=3.141592653589793D0,
     $     DEG2RAD=PI/180.D0, R=6378206.4d0

      INTEGER, INTENT(IN) :: ne,no,ics,EC(3,ne),slopeflag
      REAL(8) X(no),Y(no),Hb(no),SLAM0,SFEA0,H0
      REAL(8), INTENT(OUT) :: dt_local(ne)

      INTEGER(4)   :: I,J,K
      REAL(8)      :: DH_MIN,Hb_Average,C,Stab_MIN
      REAL(8)      :: CEDGE_X(3),CEDGE_Y(3),DX(3),DY(3),DH(3),Stab(3)
      REAL(8)      :: HBEL(3),EDGELEN(3)


      dt_local(:)=1d6; 
      if (ics.eq.2) then
         SLAM0=SLAM0*DEG2RAD
         SFEA0=SFEA0*DEG2RAD
         do i=1,no
            X(i)=R*(X(i)*DEG2RAD-SLAM0)*COS(SFEA0)
            Y(i)=Y(i)*R*DEG2RAD
         enddo
      endif

      DO I=1,ne

      	! Average Bathymetry
        do j=1,3
c           HBEL(j)=max(abs(HB(EC(j,i))),H0)
           HBEL(j)=max(abs(HB(EC(j,i))),2.)
        enddo
      	Hb_Average = (1.D0/3.D0)*(HBEL(1)+HBEL(2)+HBEL(3))
c  account for initial water height
     $       +Z0
	!Hb_Average = 1.d0/G;
      	C = Sqrt(G*Hb_Average)

      	! EDGE 1 - Node 1 <--> Node 2
      	IF(X(EC(1,I))<X(EC(2,I))) THEN
      		CEDGE_X(1) = (X(EC(2,I))-X(EC(1,I)))/2.D0 + X(EC(1,I))
      	ELSE
      		CEDGE_X(1) = (X(EC(1,I))-X(EC(2,I)))/2.D0 + X(EC(2,I))
      	ENDIF
      	IF(Y(EC(1,I))<Y(EC(2,I))) THEN
      		CEDGE_Y(1) = (Y(EC(2,I))-Y(EC(1,I)))/2.D0 + Y(EC(1,I))
      	ELSE
      		CEDGE_Y(1) = (Y(EC(1,I))-Y(EC(2,I)))/2.D0 + Y(EC(2,I))
      	ENDIF
        EDGELEN(1)=SQRT((X(EC(2,I))-X(EC(1,I)))**2 
     $       + (Y(EC(1,I))-Y(EC(2,I)))**2)
      	! EDGE 2 - Node 2 <--> Node 3
      	IF(X(EC(2,I))<X(EC(3,I))) THEN
      		CEDGE_X(2) = (X(EC(3,I))-X(EC(2,I)))/2.D0 + X(EC(2,I))
      	ELSE
      		CEDGE_X(2) = (X(EC(2,I))-X(EC(3,I)))/2.D0 + X(EC(3,I))
      	ENDIF
      	IF(Y(EC(2,I))<Y(EC(3,I))) THEN
      		CEDGE_Y(2) = (Y(EC(3,I))-Y(EC(2,I)))/2.D0 + Y(EC(2,I))
      	ELSE
      		CEDGE_Y(2) = (Y(EC(2,I))-Y(EC(3,I)))/2.D0 + Y(EC(3,I))
      	ENDIF
        EDGELEN(2)=SQRT((X(EC(2,I))-X(EC(3,I)))**2 
     $       + (Y(EC(3,I))-Y(EC(2,I)))**2)
      	! EDGE 3 - Node 3 <--> Node 1
      	IF(X(EC(1,I))<X(EC(3,I))) THEN
      		CEDGE_X(3) = (X(EC(3,I))-X(EC(1,I)))/2.D0 + X(EC(1,I))
      	ELSE
      		CEDGE_X(3) = (X(EC(1,I))-X(EC(3,I)))/2.D0 + X(EC(3,I))
      	ENDIF
      	IF(Y(EC(1,I))<Y(EC(3,I))) THEN
      		CEDGE_Y(3) = (Y(EC(3,I))-Y(EC(1,I)))/2.D0 + Y(EC(1,I))
      	ELSE
      		CEDGE_Y(3) = (Y(EC(1,I))-Y(EC(3,I)))/2.D0 + Y(EC(3,I))
      	ENDIF
        EDGELEN(3)=SQRT((X(EC(3,I))-X(EC(1,I)))**2 
     $       + (Y(EC(1,I))-Y(EC(3,I)))**2)

      	CENTROID_X(i) =(1.D0/3.D0)*(X(EC(1,I)) + X(EC(2,I)) + X(EC(3,I)))
      	CENTROID_Y(i) =(1.D0/3.D0)*(Y(EC(1,I)) + Y(EC(2,I)) + Y(EC(3,I)))

      	DO J=1,3
      		DX(J)=ABS(CENTROID_X(i) - CEDGE_X(J))
      		DY(J)=ABS(CENTROID_Y(i) - CEDGE_Y(J))
      		DH(J)=SQRT(DX(J)**2 + DY(J)**2)
                DH(J)=MIN(DH(J),EDGELEN(J))	
      
      		! For RK2/1
                if (slopeflag.eq.0) then
                   Stab(J)=cflfac*(1.D0/3.D0)*(1.D0/Sqrt(2.D0))*DH(J)/C
                endif
      		! For RK3/2
                if (slopeflag.ne.0) then
                   Stab(J)=cflfac*(1.D0/Sqrt(2.D0))*DH(J)/C
                endif

      		IF (Stab(J)<dt_local(I)) dt_local(i)=Stab(J)

      	ENDDO
c        write(700,*) i,dt_local(i)
      ENDDO
c
      if (ics.eq.2) then
         do i=1,no
            X(i)=(X(i)/(R*COS(SFEA0))+SLAM0)/DEG2RAD
            Y(i)=(Y(i)/R)/DEG2RAD
         enddo
         SLAM0=SLAM0/DEG2RAD
         SFEA0=SFEA0/DEG2RAD
      endif

      RETURN
      END SUBROUTINE

