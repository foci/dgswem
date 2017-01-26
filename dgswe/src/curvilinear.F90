      SUBROUTINE curvilinear()

      USE globals, ONLY: rp,fbseg,fbnds,nnfbed,nfbedn,nfbednn, &
                         ged2el,ged2led,ged2nn, &
                         nverts,el_type,elxy,xy,bndxy,ect,elhb,ne,mnnds,nnds, &
                         nel_type,np,nverts,psiv,psic
      USE read_dginp, ONLY: ctp
      USE curvilinear_nodes_mod
      USE bathymetry_interp_mod

      IMPLICIT NONE
      
      INTEGER :: ed,led,nd,pt,i
      INTEGER :: ged,el,ind,seg
      INTEGER :: et,typ,eo
      INTEGER :: nvert,nnd,npts,nt,nq
      INTEGER :: n1,n1ind
      INTEGER :: space
      REAL(rp) :: xpt,ypt,ytest,hb
      REAL(rp) :: x(mnnds),y(mnnds)      
      REAL(rp) :: rq(mnnds),sq(mnnds)
      REAL(rp) :: rt(mnnds),st(mnnds)  
      REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: xyhb
      
      
      
      CALL shape_functions_linear_at_ctp(nel_type,np,psiv)             
      
      CALL eval_coordinates_curved(ctp,nnds,nverts,el_type,xy,ect,fbseg,fbnds, &
                                   nnfbed,nfbedn,nfbednn,ged2el,ged2led, &
                                   psiv,bndxy,elxy)     
                                   
      space = 1                                   
      CALL shape_functions_eltypes_at_hbp(space,nel_type,np,psic)                                      
      
      ALLOCATE(xyhb(mnnds,ne,2))
      DO el = 1,ne
        CALL bathy_coordinates(el,nnds,nverts,el_type,elxy,psic,xyhb)      
      ENDDO 
      

      
!       OPEN(unit=242,file='bathy.d')
!       DO el = 1,ne
!       
!         et = el_type(el)     
!         IF (mod(et,2) == 1) THEN
!           npts = nnds(5)
!         ELSE IF (mod(et,2) == 0) THEN
!           npts = nnds(6)
!         ENDIF
!         
!         DO pt = 1,npts
!           WRITE(242,"(3(e24.17,1x))") xyhb(pt,el,1),xyhb(pt,el,2),elhb(pt,el)
!         ENDDO
!       ENDDO
!       CLOSE(242)
! 
!       OPEN(unit=242,file='elxy.d')     
!       WRITE(242,*) ne,mnnds
!       DO el = 1,ne     
!         WRITE(242,"(10(e24.17,1x))") (elxy(pt,el,1), pt = 1,mnnds)
!       ENDDO
!       DO el = 1,ne     
!         WRITE(242,"(10(e24.17,1x))") (elxy(pt,el,2), pt = 1,mnnds)
!       ENDDO      
!       CLOSE(242)      
      
      RETURN
      END SUBROUTINE curvilinear
