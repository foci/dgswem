      SUBROUTINE interp_forcing()

      USE globals, ONLY: rp,nqpte,qpte, &
                         nope,obseg,nbou,fbseg, &
                         nobed,nfbed,nobfr,nfbfr, &
                         obamp,obph,obeq,fbamp,fbph,fbeq, &
                         obfreq,obper,fbfreq,fbper, &
                         obamp_qpt,obph_qpt,obdepth_qpt,fbamp_qpt,fbph_qpt, &
                         depth,obnds, &
                         deg2rad,pi, &
                         nfbsfr,fbsamp,fbsamp_qpt,fbsbgn,fbsend,fbssig
                         
      USE allocation, ONLY: alloc_forcing_arrays 
      USE messenger2, ONLY: myrank        

      IMPLICIT NONE

      INTEGER :: pt,bfr,seg,k
      INTEGER :: nd,ed
      INTEGER :: segtype
      REAL(rp) :: xi,l1,l2
      LOGICAL :: any_nfb
      
      IF (myrank == 0) PRINT "(A)", "Interoplating boundary forcing..."
      
      ! Convert degrees to radians and frequencies to periods
      
      DO bfr = 1,nobfr
        obeq(bfr) = obeq(bfr)*deg2rad
        IF(obfreq(bfr) == 0d0) THEN
          obper(bfr) = 1d0
        ELSE
          obper(bfr) = 2d0*pi/obfreq(bfr)
        ENDIF
      ENDDO

      DO bfr = 1,nobfr
        DO seg = 1,nope
          DO nd = 1,obseg(seg)
            obph(nd,seg,bfr) = obph(nd,seg,bfr)*deg2rad
          ENDDO
        ENDDO
      ENDDO


      
      any_nfb = .false. ! determine if there are normal flow boundaries
      DO seg = 1,nbou
        segtype = fbseg(2,seg)
        IF(segtype == 2 .OR. segtype == 12 .OR. segtype == 22)THEN
          any_nfb = .true.
        ENDIF        
      ENDDO      

      IF (any_nfb) THEN

        DO bfr = 1,nfbfr
          fbeq(bfr) = fbeq(bfr)*deg2rad
          IF(fbfreq(bfr) == 0d0) THEN
            fbper(bfr) = 1d0
          ELSE
            fbper(bfr) = 2d0*pi/fbfreq(bfr)        
          ENDIF
        ENDDO

        DO bfr = 1,nfbfr
          DO seg = 1,nbou
            segtype = fbseg(2,seg)
            IF(segtype == 2 .OR. segtype == 12 .OR. segtype == 22)THEN
              DO nd = 1,fbseg(1,seg)
                fbph(nd,seg,bfr) = fbph(nd,seg,bfr)*deg2rad
              ENDDO
            ENDIF
          ENDDO
        ENDDO
        
        DO bfr = 1,nfbsfr
          fbssig(bfr) = fbssig(bfr)*86400d0
          fbsbgn(bfr) = fbsbgn(bfr)*86400d0 + 3d0*fbssig(bfr)
          fbsend(bfr) = fbsend(bfr)*86400d0 - 3d0*fbssig(bfr)
        ENDDO        
        
      ENDIF
             

      ! Interpolate boundary forcing data 
      ! qpte(1:nqpte,2,1) are the 1-D edge quadrature points 
      
      CALL alloc_forcing_arrays(4)

      DO pt = 1,nqpte(1)
        xi = qpte(pt,2,1)
        l1 = .5d0*(1d0-xi)
        l2 = .5d0*(1d0+xi)
        DO bfr = 1,nobfr
          ed = 1
          DO seg = 1,nope
            DO nd = 1,obseg(seg)-1
              obamp_qpt(bfr,pt,ed) = obamp(nd,seg,bfr)*l1 + obamp(nd+1,seg,bfr)*l2
              obph_qpt(bfr,pt,ed)  =  obph(nd,seg,bfr)*l1 +  obph(nd+1,seg,bfr)*l2
!               obdepth_qpt(ed,pt) = depth(obnds(k,seg))*l1 + depth(obnds(k+1,seg))*l2
              ed = ed + 1
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      DO pt = 1,nqpte(1)
        xi = qpte(pt,2,1)
        l1 = .5d0*(1d0-xi)
        l2 = .5d0*(1d0+xi)
        DO bfr = 1,nfbfr
          ed = 1
          DO seg = 1,nbou
            segtype = fbseg(2,seg)
            IF(segtype == 2 .OR. segtype == 12 .OR. segtype == 22)THEN
              DO nd = 1,fbseg(1,seg)-1
                fbamp_qpt(bfr,pt,ed) = fbamp(nd,seg,bfr)*l1 + fbamp(nd+1,seg,bfr)*l2
                fbph_qpt(bfr,pt,ed)  =  fbph(nd,seg,bfr)*l1 +  fbph(nd+1,seg,bfr)*l2
                ed = ed + 1
              ENDDO
            ENDIF
          ENDDO
        ENDDO
        
        DO bfr = 1,nfbsfr
          ed = 1
          DO seg = 1,nbou
            segtype = fbseg(2,seg)
            IF(segtype == 2 .OR. segtype == 12 .OR. segtype == 22)THEN
              DO nd = 1,fbseg(1,seg)-1
                fbsamp_qpt(bfr,pt,ed) = fbsamp(nd,seg,bfr)*l1 + fbsamp(nd+1,seg,bfr)*l2
                ed = ed + 1
              ENDDO
            ENDIF
          ENDDO
        ENDDO        
      ENDDO
      
    

      RETURN
      END SUBROUTINE