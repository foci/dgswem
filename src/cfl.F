      SUBROUTINE cfl(ne,no,ics,X,Y,Hb,EC,SLAM0,SFEA0,dt_local,slopeflag
     $     ,myproc)
      USE GLOBAL, ONLY:  H0
      USE DG, ONLY: CENTROID_X, CENTROID_Y, UFLUX, ZE,cflfac, Z0
      IMPLICIT NONE

      REAL(8), PARAMETER :: G=9.80665D0,PI=3.141592653589793D0,DEG2RAD=PI/180.D0

      INTEGER(4), INTENT(IN) :: ne,no,ics,EC(ne,3),slopeflag,myproc
      REAL(8), INTENT(IN) :: X(no),Y(no),Hb(no),SLAM0,SFEA0
      REAL(8), INTENT(OUT) :: dt_local(ne)

      INTEGER(4)         :: I,J,K
      REAL(8)            :: DH_MIN,Hb_Average,C,Stab_MIN
      REAL(8)            :: CEDGE_X(3),CEDGE_Y(3),DX(3),DY(3),DH(3),Stab(3)
      REAL(8)            :: HBEL(3)


      dt_local(:)=1d6; 

      DO I=1,ne

      	! Average Bathymetry
        do j=1,3
c           HBEL(j)=max(abs(HB(EC(i,j))),H0)
           HBEL(j)=HB(EC(i,j))
c        write(*,*) i,j,hbel(j)
        enddo
      	Hb_Average = (1.D0/3.D0)*(HBEL(1)+HBEL(2)+HBEL(3))
c      	IF (Hb_Average<1D-6) Hb_Average=1D-6
	!Hb_Average = 1.d0/G;
      	C = abs(uflux(i))+Sqrt(G*Hb_Average+Z0)
c      	C = Sqrt(G*Hb_Average)
		!WRITE(99,*) I,C
c        if (myproc.eq.141.and.i.eq.14176) write(987+myproc,*) C
      	! EDGE 1 - Node 1 <--> Node 2
      	IF(X(EC(I,1))<X(EC(I,2))) THEN
      		CEDGE_X(1) = (X(EC(I,2))-X(EC(I,1)))/2.D0 + X(EC(I,1))
      	ELSE
      		CEDGE_X(1) = (X(EC(I,1))-X(EC(I,2)))/2.D0 + X(EC(I,2))
      	ENDIF
      	IF(Y(EC(I,1))<Y(EC(I,2))) THEN
      		CEDGE_Y(1) = (Y(EC(I,2))-Y(EC(I,1)))/2.D0 + Y(EC(I,1))
      	ELSE
      		CEDGE_Y(1) = (Y(EC(I,1))-Y(EC(I,2)))/2.D0 + Y(EC(I,2))
      	ENDIF

      	! EDGE 2 - Node 2 <--> Node 3
      	IF(X(EC(I,2))<X(EC(I,3))) THEN
      		CEDGE_X(2) = (X(EC(I,3))-X(EC(I,2)))/2.D0 + X(EC(I,2))
      	ELSE
      		CEDGE_X(2) = (X(EC(I,2))-X(EC(I,3)))/2.D0 + X(EC(I,3))
      	ENDIF
      	IF(Y(EC(I,2))<Y(EC(I,3))) THEN
      		CEDGE_Y(2) = (Y(EC(I,3))-Y(EC(I,2)))/2.D0 + Y(EC(I,2))
      	ELSE
      		CEDGE_Y(2) = (Y(EC(I,2))-Y(EC(I,3)))/2.D0 + Y(EC(I,3))
      	ENDIF

      	! EDGE 3 - Node 3 <--> Node 1
      	IF(X(EC(I,1))<X(EC(I,3))) THEN
      		CEDGE_X(3) = (X(EC(I,3))-X(EC(I,1)))/2.D0 + X(EC(I,1))
      	ELSE
      		CEDGE_X(3) = (X(EC(I,1))-X(EC(I,3)))/2.D0 + X(EC(I,3))
      	ENDIF
      	IF(Y(EC(I,1))<Y(EC(I,3))) THEN
      		CEDGE_Y(3) = (Y(EC(I,3))-Y(EC(I,1)))/2.D0 + Y(EC(I,1))
      	ELSE
      		CEDGE_Y(3) = (Y(EC(I,1))-Y(EC(I,3)))/2.D0 + Y(EC(I,3))
      	ENDIF

      	CENTROID_X(i) = (1.D0/3.D0)*(X(EC(I,1)) + X(EC(I,2)) + X(EC(I,3)))
      	CENTROID_Y(i) = (1.D0/3.D0)*(Y(EC(I,1)) + Y(EC(I,2)) + Y(EC(I,3)))

      	DO J=1,3
      		DX(J)=ABS(CENTROID_X(i) - CEDGE_X(J))
      		DY(J)=ABS(CENTROID_Y(i) - CEDGE_Y(J))
      		DH(J)=SQRT(DX(J)**2 + DY(J)**2)
	
		!WRITE(*,*) i,j,DH(J)
      
                if (slopeflag.eq.0) then
                   Stab(J)=cflfac*(1.D0/3.D0)*(1.D0/Sqrt(2.D0))*DH(J)/C
                endif
                if (slopeflag.ne.0) then
                   Stab(J)=cflfac*(1.D0/Sqrt(2.D0))*DH(J)/C
                endif

 

      		!Stab(J)=(0.5882D0)*(1.D0/Sqrt(2.D0))*DH(J)/C

      		IF (Stab(J)<dt_local(I)) dt_local(i)=Stab(J)

      	ENDDO
               IF (HB_Average <0.D0) THEN
                 dt_local(i) = 0.0000001D0
              ENDIF
      ENDDO


      RETURN
      END SUBROUTINE

