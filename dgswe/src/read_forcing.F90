      SUBROUTINE read_forcing()

      USE globals, ONLY: rp,neta,nvel,nope,nbou,obseg,fbseg, &
                         nobfr,obtag,obtag2,obfreq,obnfact,obeq,obamp,obph, &
                         nfbfr,fbtag,fbtag2,fbfreq,fbnfact,fbeq,fbamp,fbph, &
                         fbper,obper,pi,deg2rad, &
                         nfbsfr,fbsbgn,fbsend,fbssig,fbsamp,fbstag,fbstag2
                         
      USE allocation, ONLY: alloc_forcing_arrays   
      USE messenger2, ONLY: myrank
      USE quit, ONLY: abort
      USE read_dginp, ONLY: forcing_file

      IMPLICIT NONE
      INTEGER :: bfr,nd,seg,segtype
      LOGICAL :: file_exists,any_nfb         
      
      
      INQUIRE(FILE=forcing_file, EXIST = file_exists)
      IF(file_exists == .FALSE.) THEN
        PRINT*, "fort.15 file does not exist"
        CALL abort()
      ENDIF
      
         
      
      OPEN(UNIT=15, FILE=forcing_file)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Read in open boundary periodic forcing data
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      READ(15,*) nobfr

      CALL alloc_forcing_arrays(1)

      DO bfr = 1,nobfr
        READ(15,*) obtag(bfr)
        READ(15,*) obfreq(bfr),obnfact(bfr),obeq(bfr)
      ENDDO

      DO bfr = 1,nobfr
        READ(15,*) obtag2(bfr) 
        DO seg = 1,nope
          DO nd = 1,obseg(seg)
            READ(15,*) obamp(nd,seg,bfr),obph(nd,seg,bfr)
          ENDDO
        ENDDO
      ENDDO

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Read in flow boundary periodic forcing data
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      any_nfb = .false. ! determine if there are normal flow boundaries
      DO seg = 1,nbou
        segtype = fbseg(2,seg)
        IF(segtype == 2 .OR. segtype == 12 .OR. segtype == 22)THEN
          any_nfb = .true.
        ENDIF        
      ENDDO      

      nfbfr = 0
      IF (any_nfb) THEN
        READ(15,*) nfbfr

        CALL alloc_forcing_arrays(2)

        DO bfr = 1,nfbfr
          READ(15,*) fbtag(bfr)     
          READ(15,*) fbfreq(bfr),fbnfact(bfr),fbeq(bfr)
        ENDDO


        DO bfr = 1,nfbfr
          READ(15,*) fbtag2(bfr) 
          DO seg = 1,nbou
            segtype = fbseg(2,seg)
            IF(segtype == 2 .OR. segtype == 12 .OR. segtype == 22)THEN
              DO nd = 1,fbseg(1,seg)
                READ(15,*) fbamp(nd,seg,bfr),fbph(nd,seg,bfr)
              ENDDO
            ENDIF
          ENDDO
        ENDDO
        
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Read in flow boundary surge forcing data
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
        
        READ(15,*) nfbsfr
        
        CALL alloc_forcing_arrays(3)
        
        DO bfr = 1,nfbsfr
          READ(15,*) fbstag(bfr)
          READ(15,*) fbsbgn(bfr),fbsend(bfr),fbssig(bfr)
        ENDDO
        
        DO bfr = 1,nfbsfr
          READ(15,*) fbstag2(bfr)
          DO seg = 1,nbou
            segtype = fbseg(2,seg)
            IF(segtype == 2 .OR. segtype == 12 .OR. segtype == 22)THEN
              DO nd = 1,fbseg(1,seg)
                READ(15,*) fbsamp(nd,seg,bfr)
              ENDDO
            ENDIF          
          ENDDO
        ENDDO
        
        
      ENDIF
      
      CLOSE(15)
      
      
      
      
      
      
      
      IF (myrank == 0) THEN
        PRINT "(A,I5)", 'Number of open boundary forcings',nobfr  
        IF (nobfr > 0) THEN
          DO bfr = 1,nobfr
            PRINT "(A,A)", "  ",obtag(bfr)
          ENDDO      
        ENDIF
        PRINT "(A,I5)", 'Number of flow boundary forcings',nfbfr      
        IF (nfbfr > 0) THEN
          DO bfr = 1,nfbfr
            PRINT "(A,A)", "  ",fbtag(bfr)        
          ENDDO      
        ENDIF
        PRINT "(A,I5)", 'Number of flow boundary surge forcings',nfbsfr      
        IF (nfbfr > 0) THEN
          DO bfr = 1,nfbsfr
            PRINT "(A,A)", "  ",fbstag(bfr)        
          ENDDO      
        ENDIF        
      ENDIF


      RETURN
      END SUBROUTINE read_forcing