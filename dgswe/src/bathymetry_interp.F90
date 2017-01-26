
      SUBROUTINE bathymetry_interp_area_qpts()


      USE globals, ONLY: rp,el_type,nnds,nqpta,order, &
                         ne,elxy,elhb, &
                         hbqpta_init,dhbdx_init,dhbdy_init, &
                         psia,dpsidr,dpsids, &
                         r_earth
      USE transformation, ONLY: element_transformation
      USE spherical_mod, ONLY: cpp_factor
      USE read_dginp, ONLY: coord_sys,slam0,sphi0
      USE bathymetry_interp_mod, ONLY: bathymetry_interp_eval
                         

      IMPLICIT NONE
      
      INTEGER :: el,pt,nd
      INTEGER :: typ,et,eo
      REAL(rp) :: hb
      REAL(rp) :: xpt,ypt,detJ,Sp
      REAL(rp) :: drdx,drdy,dsdx,dsdy



!       OPEN(unit=45, file='dhb.d')
!       OPEN(unit=46, file='hb.d')
!       WRITE(45,*) ne,16
!       WRITE(46,*) ne, 16
    
      
      
      DO el = 1,ne
      
        et = el_type(el)
        typ = et + 4
        eo = order(typ)        
        
!         PRINT "(I5,200(f10.4))", et,(elhb(nd,el),   nd = 1,nnds(eo))          
        
     pts: DO pt = 1,nqpta(et)        
    

            CALL element_transformation(nnds(et),elxy(:,el,1),elxy(:,el,2),psia(:,pt,et),xpt,ypt, &
                                        dpsidr(:,pt,et),dpsids(:,pt,et),drdx,drdy,dsdx,dsdy,detJ)
            
            CALL cpp_factor(coord_sys,r_earth,slam0,sphi0,ypt,Sp)                                                                                                                                              
            
            CALL bathymetry_interp_eval(nnds(eo),elhb(:,el),psia(:,pt,typ),hbqpta_init(el,pt), &
                                        dpsidr(:,pt,typ),dpsids(:,pt,typ),drdx,drdy,dsdx,dsdy,Sp, &
                                        dhbdx_init(el,pt),dhbdy_init(el,pt))
                     
!              IF (et == 1) THEN        
!                WRITE(45,"(I7,4(e24.17,1x))") el,xpt,ypt,dhbdx_init(el,pt),dhbdy_init(el,pt)      
!                WRITE(46,"(I7,4(e24.17,1x))") el,xpt,ypt,hbqpta_init(el,pt)
!              ENDIF
             
          ENDDO pts          
          
                 
                  
      ENDDO
      
!       CLOSE(45)
!       CLOSE(46)      
      
      RETURN

      END SUBROUTINE bathymetry_interp_area_qpts

      
           
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!              
      
      SUBROUTINE bathymetry_interp_edge_qpts()
      
      USE globals, ONLY: rp,el_type,ed_type,nnds,nqpta,nqpte,order, &
                         ned,elhb,ged2el,ged2led, &
                         hbqpte_init,hbqpted, &
                         psia, &
                         recv_edge  
      USE bathymetry_interp_mod, ONLY: bathymetry_interp_eval                         


      IMPLICIT NONE
      
      INTEGER :: i,j,pt,nd,ed
      INTEGER :: typ,et,eo,edpt
      INTEGER :: el_in,el_ex,led_in,led_ex,gp_in,gp_ex
      INTEGER :: nqa,nqe
      REAL(rp) :: hb
      
!       OPEN(unit=47, file='hbe.d')     
      
      DO ed = 1,ned
      
        el_in = ged2el(1,ed)
        led_in = ged2led(1,ed)          
        el_ex = ged2el(2,ed)
        led_ex = ged2led(2,ed)           
          
        IF (ed_type(ed) == 0 .OR. ed_type(ed) == -1) THEN    

          et = el_type(el_in)
          
          ! Interior edge or recieve edge => assumed to be straight
          IF (mod(et,2) == 1) THEN
            et = 1
          ELSE IF (mod(et,2) == 0) THEN
            et = 2
          ENDIF
          
          nqa = nqpta(et)
          nqe = nqpte(et)                    
        
          typ = et + 4
          eo = order(typ)
           
          DO j = 1,nqe
            pt = nqa + (led_in-1)*nqe + j
            gp_in = (led_in-1)*nqe + j  
            gp_ex = (led_ex-1)*nqe + nqe - j + 1
            
            CALL bathymetry_interp_eval(nnds(eo),elhb(:,el_in),psia(:,pt,typ),hb)
              
            hbqpte_init(el_in,gp_in) = hb
            
            IF (ed_type(ed) /= -1) THEN
              hbqpte_init(el_ex,gp_ex) = hb
            ENDIF

            hbqpted(ed,j) = hb
            
          ENDDO
            
!           PRINT("(6(E24.17))"), (hbqpte_init(el_in,gp_in), gp_in = (led_in-1)*nqe+1,(led_in-1)*nqe+nqe)   
!           PRINT("(6(E24.17))"), (hbqpte_init(el_ex,gp_ex), gp_ex = (led_ex-1)*nqe+1,(led_ex-1)*nqe+nqe)  
!           PRINT*, "  "

        ELSE IF (ed_type(ed) == 10) THEN
        
          ! No normal flow boundary treated as curved
        
          et = el_type(el_in) 
          nqa = nqpta(et)
          nqe = nqpte(et)   

          typ = et + 4          
          eo = order(typ)      
        
          DO i = 1,nqe        
            pt = nqa + (led_in-1)*nqe+i                  
            
            CALL bathymetry_interp_eval(nnds(eo),elhb(:,el_in),psia(:,pt,typ),hbqpted(ed,i))
          ENDDO                            

!             PRINT("(4(F24.17))"), (hbqpted(ed,i), i = 1,nqe)   
            
        ELSE
        
          et = el_type(el_in) 
          
          ! Specified flow/open edge assumed to be straight
          IF (mod(et,2) == 1) THEN
            et = 1
          ELSE IF (mod(et,2) == 0) THEN
            et = 2
          ENDIF      
          
          nqa = nqpta(et)
          nqe = nqpte(et)           

          typ = et + 4          
          eo = order(typ)      
        
          DO i = 1,nqe        
            pt = nqa + (led_in-1)*nqe+i            
          
            CALL bathymetry_interp_eval(nnds(eo),elhb(:,el_in),psia(:,pt,typ),hbqpted(ed,i))
          ENDDO                            
        
          
        ENDIF
        
      ENDDO         
      
!       edpt = 0
!       DO ed = 1,ned
!       
!         el_in = ged2el(1,ed)
!         led_in = ged2led(1,ed)          
!         el_ex = ged2el(2,ed)
!         led_ex = ged2led(2,ed)  
!         
!         DO pt = 1,nqpte(1)
!         
!           gp_in = (led_in-1)*nqpte(1) + pt
!           gp_ex = (led_ex-1)*nqpte(1) + nqpte(1) - pt + 1
!           
!           IF (abs(hbqpte_init(el_in,gp_in)- hbqpte_init(el_ex,gp_ex)) > 1d-13 .and. recv_edge(ed) == 0 ) THEN
!             PRINT*, "edge batymetry values don't match"
!             PRINT*, "recv_edge =", recv_edge(ed)
!             PRINT*, hbqpte_init(el_in,gp_in),hbqpte_init(el_ex,gp_ex)   
!             edpt = edpt + 1
!           ENDIF
!           
!         
!         ENDDO
! 
!       ENDDO      
!       PRINT*, "# of unmatched points", edpt
      
!       CLOSE(47)


      RETURN
      END SUBROUTINE bathymetry_interp_edge_qpts

