      SUBROUTINE make_partitions()
      
      USE globals, ONLY:npart,grid_file
      
      IMPLICIT NONE
      
      OPEN(UNIT=10,FILE='partition.d')
      
      WRITE(10,"(I4)") npart
      WRITE(10,*) grid_file
      CLOSE(10)
      
      
      PRINT*, "Sarting ADCPREP"
      CALL SYSTEM("./adcprep < partition.d")
      
      
      END SUBROUTINE make_partitions
      
            
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
            
      
      SUBROUTINE read_partitions()
      
      USE globals, ONLY: ne,npart,tnpel,mnelpp,peln,nprel,preln
      
      IMPLICIT NONE
      
      INTEGER :: part,nd,el,i
      INTEGER :: nelg,nnodg
      INTEGER :: mnnpp
      INTEGER :: ignore1,ignore2
      INTEGER :: pn,nnp,nres
      INTEGER :: neighp_r,neighp_s
      INTEGER :: nrel,nsel
      INTEGER :: el_cnt
      INTEGER :: alloc_status
      
      CHARACTER(80) :: ignore
      
      INTEGER, ALLOCATABLE, DIMENSION(:) :: npresel
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: preseln
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: pnodes
      INTEGER, ALLOCATABLE, DIMENSION(:) :: npsel
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: pseln
      
      OPEN(UNIT=80,FILE="fort.80")
      
      READ(80,*) ignore
      READ(80,*) ignore
      READ(80,"(A80)") ignore
      
      READ(80,*) nelg,nnodg ! total number of elements and nodes
      READ(80,*) npart ! number of element partitions
      READ(80,*) mnnpp ! max. number of nodes per partition
      READ(80,*) mnelpp ! max. number of elements per partition
      
      READ(80,*) ignore1
      READ(80,*) ignore1
      READ(80,*) ignore1
      READ(80,*) ignore1,ignore2
      
      IF(nelg /= ne) THEN
        PRINT*, "Number of elements from parition file /= number of elements from grid file"
        PRINT*, "   Number of elements from partition file: ", nelg
        PRINT*, "   Number of elements from grid file: ", ne
        STOP
      ENDIF
      
      
      ALLOCATE(pnodes(mnnpp,npart),STAT = alloc_status) 
      IF(alloc_status /= 0) THEN
        PRINT*, "Allocation error: pnodes"
      ENDIF
      
      DO part = 1,npart ! read in partion nodes (these are not needed)
        READ(80,*) pn,nnp,nres
        READ(80,1130) (pnodes(nd,part), nd = 1,nnp)
      ENDDO
      
      ALLOCATE(tnpel(npart),peln(mnelpp,npart),STAT = alloc_status) 
      IF(alloc_status /= 0) THEN
        PRINT*, "Allocation error: tnpel,peln"
      ENDIF
      
      DO part = 1,npart ! read in partition elements
        READ(80,*) pn,tnpel(part)
        IF (pn /= part-1) THEN
          PRINT*, "Error reading fort.80"
          STOP
        ENDIF
        READ(80,1130) (peln(el,part), el=1,tnpel(part))
      ENDDO
      
      CLOSE(80)
      
!       PRINT*,"Partition elements"
!       DO part = 1,npart
!         PRINT*, part,tnpel(part)
!         PRINT 1130, (peln(el,part), el=1,tnpel(part))
!       ENDDO      
!       PRINT*, " "
      
      OPEN(UNIT=18,FILE="DG.18")
      
      ALLOCATE(npresel(npart),preseln(mnelpp,npart),STAT = alloc_status)
      IF(alloc_status /= 0) THEN
        PRINT*, "Allocation error: npresel,preseln"
      ENDIF      
      
      ALLOCATE(nprel(npart),preln(mnelpp,npart),STAT = alloc_status)
      IF(alloc_status /= 0) THEN
        PRINT*, "Allocation error: nprel,preln"
      ENDIF  
      nprel = 0
      
      ALLOCATE(npsel(npart),pseln(mnelpp,npart),STAT = alloc_status)
      IF(alloc_status /= 0) THEN
        PRINT*, "Allocation error: npsel,pseln"
      ENDIF  
      npsel = 0      
      
      DO part = 1,npart
      
        READ(18,3010) pn,npresel(part)
        IF (pn /= part-1) THEN
          PRINT*, "Error reading DG.18"
          STOP
        ENDIF
        
        READ(18,1130) (preseln(el,part), el=1,npresel(part))
        READ(18,3010) neighp_r,neighp_s
        
        DO i = 1,neighp_r
          READ(18,3010) pn,nrel
          READ(18,1130) (preln(el+nprel(part),part), el = 1,nrel)
          nprel(part) = nprel(part) + nrel
        ENDDO
        
        DO i = 1,neighp_s
          READ(18,3010)pn,nsel
          READ(18,1130) (pseln(el+npsel(part),part), el = 1,nsel)
          npsel(part) = npsel(part) + nsel
        ENDDO
        
      ENDDO
      
!       PRINT*, "Partition recieve elements"
!       DO part = 1,npart
!         PRINT*, part,nprel(part)
!         PRINT 1130, (preln(el,part), el = 1,nprel(part))
!       ENDDO
!       PRINT*, " "
      
      CLOSE(18)
      
 1130 FORMAT(8X,9I8)  
 3010 FORMAT(8X,2I8)
      
      RETURN
      END SUBROUTINE read_partitions
      
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
      
      
      SUBROUTINE align_partitions()

      USE globals, ONLY: ne,ndof,nqpta,npart,tnpel,npel,peln,mnelpp,preln,psplit, &
                         lel2gel,gel2ael,ael2gel,gel2part, &
                         H,Hinit,Qx,Qxinit,Qy,Qyinit, &
                         dpdx,dpdx_init,dpdy,dpdy_init, &
                         dhbdx,dhbdx_init,dhbdy,dhbdy_init, &
                         Hwrite,Qxwrite,Qywrite, &
                         sp,nsp,split,mnpel
      
      IMPLICIT NONE
      
      INTEGER :: part,el,dof
      INTEGER :: el_cnt,tel
      INTEGER :: yes_rel,j
      INTEGER :: alloc_status
      INTEGER, ALLOCATABLE, DIMENSION(:) :: elflag
      
      ALLOCATE(elflag(ne))
      elflag = 0
      
      ALLOCATE(lel2gel(mnelpp,npart),npel(npart),gel2ael(ne),ael2gel(ne),gel2part(ne),STAT = alloc_status)
      IF(alloc_status /= 0) THEN
        PRINT*, "Allocation error: npel,lel2gel,gel2ael,ael2gel"
      ENDIF       
      
      ! throw out recv (ghost or halo) elements and generate local element to global element table
      DO part = 1,npart 
      
        el_cnt = 0
        DO el = 1,tnpel(part)
        
          yes_rel = 0
          DO j = 1,tnpel(part)
            IF (el == preln(j,part)) THEN
              yes_rel = 1
            ENDIF
          ENDDO
        
!           IF( ANY(preln(:,part).eq.el) ) THEN
          IF (yes_rel == 1) THEN
            ! ignore if element is recv
          ELSE
            IF (elflag(peln(el,part)) < 1) THEN
              el_cnt = el_cnt + 1
              lel2gel(el_cnt,part) = peln(el,part)
              gel2part(lel2gel(el_cnt,part)) = part
              elflag(lel2gel(el_cnt,part)) = elflag(lel2gel(el_cnt,part)) + 1
            ENDIF
          ENDIF
          
        ENDDO
        npel(part) = el_cnt
        
      ENDDO
      
!       PRINT*, "Total partition elements, res+send only"
!       DO part = 1,npart
!         PRINT*, tnpel(part),npel(part)
!         DO el = 1,npel(part)
!           PRINT*, lel2gel(el,part)
!         ENDDO
!       ENDDO

!       DO el = 1,ne
!         PRINT*, el,elflag(el),gel2part(el)
!       ENDDO


      
      ! Check if sum of partition elements = number of elements in mesh
      tel = 0
      DO part = 1,npart
        tel = tel + npel(part)
      ENDDO
      
      IF (tel /= ne) THEN
        PRINT*, "Error: Sum of partition elements /= total elements"
        PRINT*, "Total partition element count ", tel
        PRINT*, "Total number of elements", ne
        STOP
      ENDIF
      
      DO el = 1,ne
        IF (elflag(el) /= 1) THEN
          PRINT*, "ELEMENT ERROR"
          STOP
        ENDIF
      ENDDO
      
      ! rearrange dof arrays and genrate aligned element to global element table and global to aligned
      el_cnt = 1
      DO part = 1,npart
        DO el = 1,npel(part)
          H(el_cnt,1:ndof) = Hinit(lel2gel(el,part),1:ndof)
          Qx(el_cnt,1:ndof) = Qxinit(lel2gel(el,part),1:ndof)
          Qy(el_cnt,1:ndof) = Qyinit(lel2gel(el,part),1:ndof)
          
          dpdx(el_cnt,1:ndof*nqpta) = dpdx_init(lel2gel(el,part),1:ndof*nqpta)
          dpdy(el_cnt,1:ndof*nqpta) = dpdy_init(lel2gel(el,part),1:ndof*nqpta)      
          
          dhbdx(el_cnt) = dhbdx_init(lel2gel(el,part))
          dhbdy(el_cnt) = dhbdy_init(lel2gel(el,part))
          
          ael2gel(el_cnt) = lel2gel(el,part)
          gel2ael(lel2gel(el,part)) = el_cnt
          
          el_cnt = el_cnt + 1
        ENDDO
      ENDDO
      
      DEALLOCATE(Hinit,Qxinit,Qyinit)
      DEALLOCATE(dpdx_init,dpdy_init)
      DEALLOCATE(dhbdx_init,dhbdy_init)
      
!       PRINT*, " "
!       PRINT*, "Aligned element, global element"
!       DO el = 1,ne
!         PRINT*, el, ael2gel(el) !, gel2ael(el)
!       ENDDO
!       
!       PRINT*, " "
!       PRINT*, "Element, H"
!       DO el = 1,ne
!         PRINT*, el, H(el,1) !, gel2ael(el)
!       ENDDO
      
!       PRINT*, " "
!       PRINT*, "global element, aligned element"
!       DO el = 1,ne
!         PRINT*, el, gel2ael(el)
!       ENDDO
      
      ALLOCATE(Hwrite(ne,ndof),Qxwrite(ne,ndof),Qywrite(ne,ndof),STAT = alloc_status)
      IF(alloc_status /= 0) THEN
        PRINT*, "Allocation error: Hwrite,Qxwrite,Qywrite"
      ENDIF      
      
      ! set up write pointer arrays
      DO el = 1,ne
        DO dof = 1,ndof
          Hwrite(el,dof)%ptr => H(gel2ael(el),dof)
          Qxwrite(el,dof)%ptr => Qx(gel2ael(el),dof)
          Qywrite(el,dof)%ptr => Qy(gel2ael(el),dof)          
        ENDDO
      ENDDO
      
!       PRINT*, " "
!       PRINT*, "Element, Hwrite"
!       DO el = 1,ne
!         PRINT*, el, Hwrite(el,1)%ptr !, gel2ael(el)
!       ENDDO

      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", "           Loop Blocking Information         "
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", " "

      ! determine integration element loop blocking info
      PRINT*,"split: "
      DO sp = 1,nsp
        split(1,sp) = (sp-1)*(ne/nsp) + 1
        split(2,sp) = sp*(ne/nsp)
      ENDDO
      split(2,nsp) = ne

      DO sp = 1,nsp
        PRINT*, split(1,sp),split(2,sp) 
      ENDDO
      PRINT*, " "
      
      PRINT*, "Number of partitions: ", npart
      
      ! determine edge integration element loop blocking bounds
      ALLOCATE(psplit(2,npart),STAT = alloc_status)
      IF(alloc_status /= 0) THEN
        PRINT*, "Allocation error: psplit"
      ENDIF         
      
      tel = 0
      psplit(1,1) = 1
      DO part = 1,npart-1      
        psplit(2,part) = tel + npel(part)
        psplit(1,part+1) = tel + npel(part) + 1 
        tel = tel + npel(part)
      ENDDO
      psplit(2,part) = tel + npel(npart)
      
      PRINT*, " "
      PRINT*, "psplit: "
      DO part = 1,npart
        PRINT*, psplit(1,part), psplit(2,part), npel(part)
      ENDDO
      mnpel = MAXVAL(npel(:))
      PRINT*, "Max elements per partition: ", mnpel
      PRINT*, " "
      
      RETURN
      END SUBROUTINE align_partitions
      
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
            
      
      SUBROUTINE edge_partition()
      
      USE globals, ONLY: npart,nied,esplit,npel, &
                         iedn,ged2el,lel2gel,piedn,mnpied
                         
      
      IMPLICIT NONE
      
      INTEGER :: part,ed,el
      INTEGER :: ged,el_in,el_ex
      INTEGER :: yes_in,yes_ex
      INTEGER :: edcnt,ted
      INTEGER :: alloc_status
      
      INTEGER, ALLOCATABLE, DIMENSION(:) :: edflag, npied     
      
      ALLOCATE(edflag(nied),npied(npart+1),piedn(nied))
      
      edflag = 0
      npied = 0
      
      edcnt = 0
      
      ! loop through each partition and find interior edges associated with those elements
      DO part = 1,npart 

        DO ed = 1,nied

          ged = iedn(ed)
          el_in = ged2el(1,ged)
          el_ex = ged2el(2,ged)
          
          yes_in = 0
          yes_ex = 0
          
          DO el = 1,npel(part)
            IF (lel2gel(el,part) == el_in) THEN
              yes_in = 1
            ENDIF
            IF (lel2gel(el,part) == el_ex) THEN
              yes_ex = 1
            ENDIF
          ENDDO
            
          ! check if both elements are in the partition
!           IF(ANY(lel2gel(:,part).eq.el_in) .and. ANY(lel2gel(:,part).eq.el_ex)) THEN
          IF (yes_in == 1 .and. yes_ex == 1) THEN
!             IF (edflag(ed) < 1) THEN
            
              edflag(ed) = edflag(ed) + 1
              edcnt = edcnt + 1
              
              npied(part) = npied(part) + 1
              piedn(edcnt) = ged
                
              CALL point_to_el(edcnt,ged,el_in,el_ex)
!             ENDIF    
          ELSE
            ! ignore edges that contain elements from two different partitions (for now)
          ENDIF

        ENDDO
          
      ENDDO
      
!       PRINT*, "edge count after partitioning = ", edcnt
      
!       PRINT*, " "
!       DO ed = 1,nied
!         PRINT*, ed, edflag(ed)
!       ENDDO

      DO ed = 1,nied
        IF(edflag(ed) > 1) THEN
          PRINT*, "Edge partition error: edflag > 1"
          STOP
        ENDIF
      ENDDO
      
      DO ed = 1,nied
        IF(edflag(ed) == 0) THEN
        
          edflag(ed) = edflag(ed) + 1
          edcnt = edcnt + 1
          
          npied(npart+1) = npied(npart+1) + 1
          
          ged = iedn(ed)
          el_in = ged2el(1,ged)
          el_ex = ged2el(2,ged)
          
          piedn(edcnt) = ged
          
          CALL point_to_el(edcnt,ged,el_in,el_ex)
          
        ENDIF
      ENDDO
      
!       PRINT*, " "
!       DO ed = 1,nied
!         PRINT*, ed, edflag(ed)
!       ENDDO

      ! check for errors
      DO ed = 1,nied
        IF(edflag(ed) /= 1) THEN
          PRINT*, "Edge partition error: edflag /= 1"
          STOP
        ENDIF
      ENDDO
      
      ted = 0
      DO part = 1,npart+1
        ted = ted + npied(part)
      ENDDO
      
      IF(ted /= nied .or. edcnt /= nied) THEN
        PRINT*, "Error: sum of partition interior edges /= total interior edges"
        PRINT*, "ted = ", ted
        PRINT*, "nied = ", nied
        PRINT*, "edcnt = ", edcnt
        STOP
      ENDIF
      
      ! determine edge integration edge loop blocking bounds      
      ALLOCATE(esplit(2,npart+1),STAT = alloc_status)
      IF(alloc_status /= 0) THEN
        PRINT*, "Allocation error: esplit"
      ENDIF     
      
      ted = 0
      esplit(1,1) = 1
      DO part = 1,npart      
        esplit(2,part) = ted + npied(part)
        esplit(1,part+1) = ted + npied(part) + 1 
        ted = ted + npied(part)
      ENDDO
      esplit(2,npart+1) = ted + npied(npart+1)
      
      PRINT*, " "
      PRINT*, "esplit: "
      DO part = 1,npart+1
        PRINT*, esplit(1,part), esplit(2,part), npied(part)
      ENDDO  
      mnpied = MAXVAL(npied(1:npart)) 
      PRINT*, "Max interior edges per partition: ", mnpied
      
      
      RETURN
      END SUBROUTINE edge_partition
      
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!            
      
      
      
      SUBROUTINE point_to_el(ed,ged,el_in,el_ex)
      
      USE globals, ONLY: nqpte,ged2led,gel2ael,gel2part, &
                         Hqpt,Qxqpt,Qyqpt, &
                         xmom,ymom,xymom, &
                         Hflux,Qxflux,Qyflux, &
                         Hi,He,Qxi,Qxe,Qyi,Qye, &
                         xmi,xme,ymi,yme,xymi,xyme, &
                         Hfi,Hfe,Qxfi,Qxfe,Qyfi,Qyfe, &
                         normal,edlen_area, &
                         inx,iny,len_area_in,len_area_ex
                      
      
      IMPLICIT NONE
      
      INTEGER :: pt
      INTEGER :: ged
      INTEGER :: ed
      INTEGER :: el_in,el_ex
      INTEGER :: led_in,led_ex
      INTEGER :: gp_in,gp_ex
      INTEGER :: ael_in,ael_ex
      
      ael_in = gel2ael(el_in)
      ael_ex = gel2ael(el_ex)
      
!       PRINT "(I8,8x,3(I8),8X,3(I8))", ed,el_in,ael_in,gel2part(el_in), el_ex,ael_ex,gel2part(el_ex)
      
      DO pt = 1,nqpte
      
        led_in = ged2led(1,ged)
        led_ex = ged2led(2,ged)

        gp_in = (led_in-1)*nqpte + pt
        gp_ex = (led_ex-1)*nqpte + nqpte - pt + 1

        Hi(ed,pt)%ptr => Hqpt(ael_in,gp_in)
        He(ed,pt)%ptr => Hqpt(ael_ex,gp_ex)

        Qxi(ed,pt)%ptr => Qxqpt(ael_in,gp_in)
        Qxe(ed,pt)%ptr => Qxqpt(ael_ex,gp_ex)

        Qyi(ed,pt)%ptr => Qyqpt(ael_in,gp_in)
        Qye(ed,pt)%ptr => Qyqpt(ael_ex,gp_ex)

        xmi(ed,pt)%ptr => xmom(ael_in,gp_in)
        xme(ed,pt)%ptr => xmom(ael_ex,gp_ex)

        ymi(ed,pt)%ptr => ymom(ael_in,gp_in)
        yme(ed,pt)%ptr => ymom(ael_ex,gp_ex)

        xymi(ed,pt)%ptr => xymom(ael_in,gp_in)
        xyme(ed,pt)%ptr => xymom(ael_ex,gp_ex)

        Hfi(ed,pt)%ptr => Hflux(ael_in,gp_in)
        Hfe(ed,pt)%ptr => Hflux(ael_ex,gp_ex)

        Qxfi(ed,pt)%ptr => Qxflux(ael_in,gp_in)
        Qxfe(ed,pt)%ptr => Qxflux(ael_ex,gp_ex)

        Qyfi(ed,pt)%ptr => Qyflux(ael_in,gp_in)
        Qyfe(ed,pt)%ptr => Qyflux(ael_ex,gp_ex)
     
      
      ENDDO
      
      inx(ed) = normal(1,ged)
      iny(ed) = normal(2,ged)
          
      len_area_in(ed) = edlen_area(1,ged)
      len_area_ex(ed) = edlen_area(2,ged)
      
      RETURN
      END SUBROUTINE point_to_el
