      SUBROUTINE edge_partition2()
      
      USE globals, ONLY: nn,ne,ndof,mndof,nqpta,mnqpta,mnqpte,part,nel_type,el_type, &
                         ged2el,&
                         npartel,nparted,npartet, &
                         gel2part,gel2lel,lel2gel, &
                         ael2gel,gel2ael, &
                         nied,iedn, &
                         H,Z,Qx,Qy,Hinit,Zinit,Qxinit,Qyinit, &
                         dpdx,dpdy,dpdx_init,dpdy_init, &
                         dhbdx,dhbdy,dhbdx_init,dhbdy_init, &
                         hbqpta,hbqpta_init,hbqpte,hbqpte_init, &
                         phia_int,phia_int_init, &
                         mmi,mmi_init, &
                         Hwrite,Zwrite,Qxwrite,Qywrite, &
                         nblk,elblk,edblk,nfblk,nrblk,rnfblk, &
                         mnpartel,mnparted   
                         
      USE allocation, ONLY: alloc_ptr_arrays,alloc_blk_arrays,dealloc_init_arrays                
      USE messenger2, ONLY: myrank,nred,redn
      USE read_dginp, ONLY: npart
      


      IMPLICIT NONE

      INTEGER :: el,pe,ed,dof,blk,et
      INTEGER :: ged
      INTEGER :: el_in,el_ex
      INTEGER :: pe_in,pe_ex
      INTEGER :: elcnt,edcnt
      INTEGER :: ted,tel
      
      INTEGER, ALLOCATABLE, DIMENSION(:) :: edflag      
      INTEGER, ALLOCATABLE, DIMENSION(:) :: elflag
      
      
      ! Find elements that have send/recieve edges
      ALLOCATE(elflag(ne))
      elflag = 0
      
      DO ed = 1,nred
        ged = redn(ed)
        el = ged2el(1,ged)
        
        IF(elflag(el) == 0) THEN
          elflag(el) = 1
        ENDIF
      ENDDO
      
      
      CALL alloc_ptr_arrays()      
      CALL alloc_blk_arrays()
      
      
      ! Create aligned element sections
      npartel = 0
      npartet = 0
      DO et = 1,nel_type
        DO el = 1,ne
          IF (el_type(el) == et) THEN
            pe = part(el) 
            
            IF (elflag(el) == 1) THEN                 ! First aligned section is for elements with send/recieve edges
              npartel(1) = npartel(1) + 1 
              npartet(et,1) = npartet(et,1) + 1  
            
              gel2lel(el) = npartel(1) 
              lel2gel(npartel(1),1) = el 
            ELSE
              npartel(pe+1) = npartel(pe+1) + 1       ! count elements in partition pe
              npartet(et,pe+1) = npartet(et,pe+1) + 1  
            
              gel2lel(el) = npartel(pe+1)             ! global element el has local element number npartel(pe) in partition pe
              lel2gel(npartel(pe+1),pe+1) = el        ! local element npartel(pe) on partition pe has global element number el            
            ENDIF
            
            gel2part(el) = pe                         ! global element el is in parition pe            
          ENDIF
        ENDDO
      ENDDO                    
      

      ! Map initial data to aligned sections
      elcnt = 0 
      DO pe = 1,npart+1
!       PRINT*, " "
        DO el = 1,npartel(pe)
          elcnt = elcnt + 1
        
          H(elcnt,1:mndof) = Hinit(lel2gel(el,pe),1:mndof)
          Z(elcnt,1:mndof) = Zinit(lel2gel(el,pe),1:mndof)
          Qx(elcnt,1:mndof) = Qxinit(lel2gel(el,pe),1:mndof)
          Qy(elcnt,1:mndof) = Qyinit(lel2gel(el,pe),1:mndof)  
          
          hbqpta(elcnt,1:mnqpta) = hbqpta_init(lel2gel(el,pe),1:mnqpta)
          hbqpte(elcnt,1:4*mnqpte) = hbqpte_init(lel2gel(el,pe),1:4*mnqpte)
          
          dpdx(elcnt,1:mndof,1:mnqpta) = dpdx_init(lel2gel(el,pe),1:mndof,1:mnqpta)
          dpdy(elcnt,1:mndof,1:mnqpta) = dpdy_init(lel2gel(el,pe),1:mndof,1:mnqpta)
          
          phia_int(elcnt,1:mndof,1:mnqpta) = phia_int_init(lel2gel(el,pe),1:mndof,1:mnqpta)
          
          dhbdx(elcnt,1:mnqpta) = dhbdx_init(lel2gel(el,pe),1:mnqpta)
          dhbdy(elcnt,1:mnqpta) = dhbdy_init(lel2gel(el,pe),1:mnqpta)
          
          mmi(elcnt,1:mndof*mndof) = mmi_init(lel2gel(el,pe),1:mndof*mndof)
          
          ael2gel(elcnt) = lel2gel(el,pe)
          gel2ael(lel2gel(el,pe)) = elcnt
          
!           PRINT*, elcnt,ael2gel(elcnt))
        ENDDO
      ENDDO
      
      CALL dealloc_init_arrays()
      
      
      ! Setup pointer arrays for output
      DO el = 1,ne
        DO dof = 1,mndof
          Hwrite(el,dof)%ptr => H(gel2ael(el),dof)
          Zwrite(el,dof)%ptr => Z(gel2ael(el),dof)
          Qxwrite(el,dof)%ptr => Qx(gel2ael(el),dof)
          Qywrite(el,dof)%ptr => Qy(gel2ael(el),dof)
        ENDDO
      ENDDO
      
      
      ! Identify all partition edges and setup pointer arrays
      ALLOCATE(edflag(nied))
      
      edflag = 0
      edcnt = 0
      nparted = 0
      DO pe = 1,npart
        DO ed = 1,nied
        
          ged = iedn(ed)
          el_in = ged2el(1,ged)
          el_ex = ged2el(2,ged)
          
          pe_in = gel2part(el_in) 
          pe_ex = gel2part(el_ex)
          
          IF (pe_in == pe .and. pe_ex == pe) THEN
            
            edflag(ed) = edflag(ed) + 1
            edcnt = edcnt + 1
            
            nparted(pe) = nparted(pe) + 1
            
            CALL point_to_el(edcnt,ged,el_in,el_ex)
            
          ENDIF
        
        ENDDO
      ENDDO
      
      ! Identify remainder edges and setup pointer arrays
      DO ed = 1,nied
        IF (edflag(ed) == 0) THEN
          edflag(ed) = edflag(ed) + 1
          edcnt = edcnt + 1
          
          ged = iedn(ed)
          el_in = ged2el(1,ged)
          el_ex = ged2el(2,ged)
          
          nparted(npart+1) = nparted(npart+1) + 1
          
          CALL point_to_el(edcnt,ged,el_in,el_ex)          
          
        ENDIF
      ENDDO
      
      
      
      
      ! check for errors
      DO ed = 1,nied
        IF(edflag(ed) /= 1) THEN
          PRINT*, "Edge partition error: edflag /= 1"
          STOP
        ENDIF
      ENDDO    
      
      ted = 0
      DO pe = 1,npart+1
        ted = ted + nparted(pe)
      ENDDO
      
      IF(ted /= nied .or. edcnt /= nied) THEN
        PRINT*, "Error: sum of partition interior edges /= total interior edges"
        PRINT*, "ted = ", ted
        PRINT*, "nied = ", nied
        PRINT*, "edcnt = ", edcnt
        STOP
      ENDIF
      
      

      ! Calculate blocking arrays
      
      tel = 0                               ! element blocking arrays
      DO blk = 1,npart+1
        elblk(1,blk,1) = tel + 1
        DO et = 1,nel_type-1
          tel = tel + npartet(et,blk)
          elblk(2,blk,et) = tel
          elblk(1,blk,et+1) = tel+1
        ENDDO
        elblk(2,blk,nel_type) = tel + npartet(nel_type,blk)        
      ENDDO
      
      ted = 0                               ! partition edge blocking arrays              
      nfblk(1,1) = 1
      DO blk = 1,npart
        ted = ted + nparted(blk)
        nfblk(2,blk) = ted
        nfblk(1,blk+1) = ted + 1
      ENDDO
      nfblk(2,npart+1) = ted + nparted(npart+1)
      
      rnfblk(1,1) = ted+1                   ! remainder edge blocking arrays
      DO blk = 1,nrblk-1
        ted = ted + (nparted(npart+1)/nrblk)
        rnfblk(2,blk) = ted  
        rnfblk(1,blk+1) = ted + 1
      ENDDO
      rnfblk(2,nrblk) = nied
      
      
      IF (myrank == 0) THEN
        PRINT "(A)", ""
        PRINT "(A)", "---------------------------------------------"
        PRINT "(A)", "           Loop Blocking Information         "
        PRINT "(A)", "---------------------------------------------"
        PRINT "(A)", " "            
          
      
        PRINT*, "Number of partitions: ", npart      
        PRINT*, " "
        PRINT*, "elblk: "
        DO blk = 1,npart+1
          PRINT*, elblk(1,blk,1), elblk(2,blk,1), npartet(1,blk)
          PRINT*, elblk(1,blk,2), elblk(2,blk,2), npartet(2,blk)
          PRINT*, elblk(1,blk,3), elblk(2,blk,3), npartet(3,blk)
          PRINT*, elblk(1,blk,4), elblk(2,blk,4), npartet(4,blk)    
          PRINT*, " "        
        ENDDO
        mnpartel = MAXVAL(npartel)
        PRINT*, "Max elements per partition: ", mnpartel      
      
        PRINT*, " "
        PRINT*, "nfblk: "
        DO blk = 1,npart+1
          PRINT*, nfblk(1,blk), nfblk(2,blk), nparted(blk)
        ENDDO
        mnparted = MAXVAL(nparted(1:npart))
        PRINT*, "Max edges per partition: ", mnparted
        PRINT*, " " 
      
        PRINT*, " "
        PRINT*, "rnfblk: "
        DO blk = 1,nrblk
          PRINT*, rnfblk(1,blk), rnfblk(2,blk)
        ENDDO      
        PRINT*, " "
      ENDIF
      
      END SUBROUTINE edge_partition2
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
      SUBROUTINE point_to_el(ed,ged,el_in,el_ex)
      
      USE globals, ONLY: nqpte,ged2led,gel2ael,gel2part, &
                         Hqpt,Zqpt,Qxqpt,Qyqpt, &
                         xmom,ymom,xymom, &
                         Hi,He,Zi,Ze,Qxi,Qxe,Qyi,Qye, &
                         xmi,xme,ymi,yme,xymi,xyme, &
                         nx_pt,ny_pt,detJe,Spe,cfac, &
                         inx,iny,detJe_in,detJe_ex,icfac, &
                         Exxi,Exxe,Eyyi,Eyye,Exyi,Exye,Eyxi,Eyxe, &
                         Exxqpt,Eyyqpt,Exyqpt,Eyxqpt
                      
      
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
      
      DO pt = 1,nqpte(1)
      
        led_in = ged2led(1,ged)
        led_ex = ged2led(2,ged)

        gp_in = (led_in-1)*nqpte(1) + pt
        gp_ex = (led_ex-1)*nqpte(1) + nqpte(1) - pt + 1

        Hi(ed,pt)%ptr => Hqpt(ael_in,gp_in)
        He(ed,pt)%ptr => Hqpt(ael_ex,gp_ex)
        
        Zi(ed,pt)%ptr => Zqpt(ael_in,gp_in)
        Ze(ed,pt)%ptr => Zqpt(ael_ex,gp_ex)        

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
        
        Exxi(ed,pt)%ptr => Exxqpt(ael_in,gp_in)
        Exxe(ed,pt)%ptr => Exxqpt(ael_ex,gp_ex)
        
        Eyyi(ed,pt)%ptr => Eyyqpt(ael_in,gp_in)
        Eyye(ed,pt)%ptr => Eyyqpt(ael_ex,gp_ex)
        
        Exyi(ed,pt)%ptr => Exyqpt(ael_in,gp_in)
        Exye(ed,pt)%ptr => Exyqpt(ael_ex,gp_ex)        
        
        Eyxi(ed,pt)%ptr => Eyxqpt(ael_in,gp_in)
        Eyxe(ed,pt)%ptr => Eyxqpt(ael_ex,gp_ex)  
        
     
        
        inx(ed,pt) = nx_pt(ged,pt)*Spe(ged,pt)
        iny(ed,pt) = ny_pt(ged,pt)
        
        icfac(ed,pt) = cfac(ged,pt)
          
        detJe_in(ed,pt) = detJe(ged,pt)
        detJe_ex(ed,pt) = detJe(ged,pt)        
     
      
      ENDDO
      
!       ! straight sided edges => constant normals
!       inx(ed) = nx_pt(ged,1)
!       iny(ed) = ny_pt(ged,1)
!           
!       ! straight sided edges => constant Jacobian
!       detJe_in(ed) = detJe(ged,1)
!       detJe_ex(ed) = detJe(ged,1)
      
      RETURN
      END SUBROUTINE point_to_el
      
