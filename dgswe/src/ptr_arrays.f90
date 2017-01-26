     SUBROUTINE ptr_arrays()
     
     USE globals, ONLY: ne,ned,nied,nqpte,iedn,ged2led,ged2el,edlen_area,normal, &
                       Hi,He,Qxi,Qxe,Qyi,Qye, &
                       xmi,xme,ymi,yme,xymi,xyme, &
                       Hfi,Hfe,Qxfi,Qxfe,Qyfi,Qyfe, &
                       Hqpt,Qxqpt,Qyqpt, &
                       xmom,ymom,xymom, &
                       Hflux,Qxflux,Qyflux, &
                       inx,iny,len_area_in,len_area_ex,const, &
                       Hhatv,Qxhatv,Qyhatv, &
                       nxv,nyv, & 
                       Hn,Qxn,Qyn,egnval, &
                       Hni,Qxni,Qyni,Hne,Qxne,Qyne,EVi,EVe, &
                       Hp,Qxp,Qyp,xmp,ymp,xymp, &
                       Hai,Hae,Qxai,Qxae,Qyai,Qyae, &
                       xmai,xmae,ymai,ymae,xymai,xymae, &
                       Hfai,Hfae,Qxfai,Qxfae,Qyfai,Qyfae

     IMPLICIT NONE
     INTEGER :: alloc_status
     INTEGER :: ed,ged,pt
     INTEGER :: gp_in,gp_ex,led_in,led_ex,el_in,el_ex

     ALLOCATE(Hi(nied,nqpte),He(nied,nqpte),Qxi(nied,nqpte),Qxe(nied,nqpte),Qyi(nied,nqpte),Qye(nied,nqpte),STAT=alloc_status)
     IF(alloc_status /= 0) THEN
       PRINT*, "Allocation error: Hi,He,Qxi,Qxe,Qyi,Qye"
     ENDIF

     ALLOCATE(xmi(nied,nqpte),xme(nied,nqpte),ymi(nied,nqpte),yme(nied,nqpte),xymi(nied,nqpte),xyme(nied,nqpte),STAT=alloc_status)
     IF(alloc_status /= 0) THEN
       PRINT*, "Allocation error: xmi,xme,ymi,yme,xymi,xyme"
     ENDIF

!      ALLOCATE(Hin(nied,nqpte),Hex(nied,nqpte),Qxin(nied,nqpte),Qxex(nied,nqpte),Qyin(nied,nqpte),Qyex(nied,nqpte),STAT=alloc_status)
!      IF(alloc_status /= 0) THEN
!        PRINT*, "Allocation error: Hin,Hex,Qxin,Qxex,Qyin,Qyex"
!      ENDIF
! 
!      ALLOCATE(xmin(nied,nqpte),xmex(nied,nqpte),ymin(nied,nqpte),ymex(nied,nqpte),xymin(nied,nqpte),xymex(nied,nqpte),STAT=alloc_status)
!      IF(alloc_status /= 0) THEN
!        PRINT*, "Allocation error: xmin,xmex,ymin,ymex,xymin,xymex"
!      ENDIF

     ALLOCATE(Hfi(nied,nqpte),Hfe(nied,nqpte),Qxfi(nied,nqpte),Qxfe(nied,nqpte),Qyfi(nied,nqpte),Qyfe(nied,nqpte),STAT=alloc_status)
     IF(alloc_status /= 0) THEN
       PRINT*, "Allocation error: Hfi,Hfe,Qxfi,Qxfe,Qyfi,Qyfe"
     ENDIF

     ALLOCATE(const(nied),inx(nied),iny(nied),len_area_in(nied),len_area_ex(nied),STAT=alloc_status)
     IF(alloc_status /= 0) THEN
       PRINT*, "Allocation error: const,inx,iny,len_area_in,len_area_ex"
     ENDIF

     ALLOCATE(Hhatv(nied),Qxhatv(nied),Qyhatv(nied),STAT=alloc_status)
     IF(alloc_status /= 0) THEN
       PRINT*, "Allocation error: Hhatv,Qxhatv,Qyhatv"
     ENDIF

     ALLOCATE(nxv(ne,3),nyv(ne,3),STAT=alloc_status)
     IF(alloc_status /= 0) THEN
       PRINT*, "Allocation error: nxv,nyv"
     ENDIF

     ALLOCATE(Hni(nied,nqpte),Hne(nied,nqpte),Qxni(nied,nqpte),Qxne(nied,nqpte),Qyni(nied,nqpte),Qyne(nied,nqpte),STAT=alloc_status)
     IF(alloc_status /= 0) THEN
       PRINT*, "Allocation error: Hni,Hne,Qxni,Qxne,Qyni,Qyne"
     ENDIF

     ALLOCATE(EVi(nied,nqpte),EVe(nied,nqpte),STAT=alloc_status)
     IF(alloc_status /= 0) THEN
       PRINT*, "Allocation error: EVi,EVe"
     ENDIF

!     ALLOCATE(rHi(nied,nqpte),rHe(nied,nqpte),STAT=alloc_status)
!     IF(alloc_status /= 0) THEN
!       PRINT*, "Allocation error: rHi,rHe"
!     ENDIF

     ALLOCATE(Hp(nied,nqpte),Qxp(nied,nqpte),Qyp(nied,nqpte),STAT=alloc_status)
     IF(alloc_status /= 0) THEN
       PRINT*, "Allocation error: Hp,Qxp,Qyp"
     ENDIF

     ALLOCATE(xmp(nied,nqpte),ymp(nied,nqpte),xymp(nied,nqpte),STAT=alloc_status)
     IF(alloc_status /= 0) THEN
       PRINT*, "Allocation error: xmp,ymp,xymp"
     ENDIF
     
     ALLOCATE(Hai(nied),Hae(nied),Qxai(nied),Qxae(nied),Qyai(nied),Qyae(nied))
     ALLOCATE(xmai(nied),xmae(nied),ymai(nied),ymae(nied),xymai(nied),xymae(nied))
     ALLOCATE(Hfai(nied),Hfae(nied),Qxfai(nied),Qxfae(nied),Qyfai(nied),Qyfae(nied))

     DO pt = 1,nqpte
       DO ed = 1,nied
       
         ged = iedn(ed)

         led_in = ged2led(1,ged)
         led_ex = ged2led(2,ged)

         gp_in = (led_in-1)*nqpte + pt
         gp_ex = (led_ex-1)*nqpte + nqpte - pt + 1
  
         el_in = ged2el(1,ged)
         el_ex = ged2el(2,ged)

         Hi(ed,pt)%ptr => Hqpt(el_in,gp_in)
         He(ed,pt)%ptr => Hqpt(el_ex,gp_ex)

         Qxi(ed,pt)%ptr => Qxqpt(el_in,gp_in)
         Qxe(ed,pt)%ptr => Qxqpt(el_ex,gp_ex)

         Qyi(ed,pt)%ptr => Qyqpt(el_in,gp_in)
         Qye(ed,pt)%ptr => Qyqpt(el_ex,gp_ex)

         xmi(ed,pt)%ptr => xmom(el_in,gp_in)
         xme(ed,pt)%ptr => xmom(el_ex,gp_ex)

         ymi(ed,pt)%ptr => ymom(el_in,gp_in)
         yme(ed,pt)%ptr => ymom(el_ex,gp_ex)

         xymi(ed,pt)%ptr => xymom(el_in,gp_in)
         xyme(ed,pt)%ptr => xymom(el_ex,gp_ex)

         Hfi(ed,pt)%ptr => Hflux(el_in,gp_in)
         Hfe(ed,pt)%ptr => Hflux(el_ex,gp_ex)

         Qxfi(ed,pt)%ptr => Qxflux(el_in,gp_in)
         Qxfe(ed,pt)%ptr => Qxflux(el_ex,gp_ex)

         Qyfi(ed,pt)%ptr => Qyflux(el_in,gp_in)
         Qyfe(ed,pt)%ptr => Qyflux(el_ex,gp_ex)

         Hni(ed,pt)%ptr => Hn(el_in,gp_in)
         Hne(ed,pt)%ptr => Hn(el_ex,gp_ex)

         Qxni(ed,pt)%ptr => Qxn(el_in,gp_in)
         Qxne(ed,pt)%ptr => Qxn(el_ex,gp_ex)

         Qyni(ed,pt)%ptr => Qyn(el_in,gp_in)
         Qyne(ed,pt)%ptr => Qyn(el_ex,gp_ex)

         EVi(ed,pt)%ptr => egnval(el_in,gp_in)
         EVe(ed,pt)%ptr => egnval(el_ex,gp_ex)

!        rHi(ed,pt)%ptr => recipHa(el_in,gp_in)
!        rHe(ed,pt)%ptr => recipHa(el_ex,gp_ex)



         Hp(ed,pt)%in => Hqpt(el_in,gp_in)
         Hp(ed,pt)%ex => Hqpt(el_ex,gp_ex)

         Qxp(ed,pt)%in => Qxqpt(el_in,gp_in)
         Qxp(ed,pt)%ex => Qxqpt(el_ex,gp_ex)

         Qyp(ed,pt)%in => Qyqpt(el_in,gp_in)
         Qyp(ed,pt)%ex => Qyqpt(el_ex,gp_ex)

         xmp(ed,pt)%in => xmom(el_in,gp_in)
         xmp(ed,pt)%ex => xmom(el_ex,gp_ex)

         ymp(ed,pt)%in => ymom(el_in,gp_in)
         ymp(ed,pt)%ex => ymom(el_ex,gp_ex)

         xymp(ed,pt)%in => xymom(el_in,gp_in)
         xymp(ed,pt)%ex => xymom(el_ex,gp_ex)

       ENDDO
     ENDDO
     
     DO ed = 1,nied
       ged = iedn(ed)

       led_in = ged2led(1,ged)
       led_ex = ged2led(2,ged)    
       
       gp_in = (led_in-1)*nqpte 
       gp_ex = (led_ex-1)*nqpte       
       
       el_in = ged2el(1,ged)
       el_ex = ged2el(2,ged)         
     
       Hai(ed)%ptr => Hqpt(el_in,gp_in+1:gp_in+nqpte)
       Hae(ed)%ptr => Hqpt(el_ex,gp_ex+nqpte:gp_in+1:-1)
       
       Qxai(ed)%ptr => Qxqpt(el_in,gp_in+1:gp_in+nqpte)
       Qxae(ed)%ptr => Qxqpt(el_ex,gp_ex+nqpte:gp_in+1:-1)   
       
       Qyai(ed)%ptr => Qyqpt(el_in,gp_in+1:gp_in+nqpte)
       Qyae(ed)%ptr => Qyqpt(el_ex,gp_ex+nqpte:gp_in+1:-1)    
       
       xmai(ed)%ptr => xmom(el_in,gp_in+1:gp_in+nqpte)
       xmae(ed)%ptr => xmom(el_ex,gp_ex+nqpte:gp_in+1:-1)    
       
       ymai(ed)%ptr => ymom(el_in,gp_in+1:gp_in+nqpte)
       ymae(ed)%ptr => ymom(el_ex,gp_ex+nqpte:gp_in+1:-1)        
       
       xymai(ed)%ptr => xymom(el_in,gp_in+1:gp_in+nqpte)
       xymae(ed)%ptr => xymom(el_ex,gp_ex+nqpte:gp_in+1:-1) 
       
       Hfai(ed)%ptr => Hflux(el_in,gp_in+1:gp_in+nqpte)
       Hfae(ed)%ptr => Hflux(el_ex,gp_ex+nqpte:gp_in+1:-1)  
       
       Qxfai(ed)%ptr => Qxflux(el_in,gp_in+1:gp_in+nqpte)
       Qxfae(ed)%ptr => Qxflux(el_ex,gp_ex+nqpte:gp_in+1:-1) 
       
       Qyfai(ed)%ptr => Qyflux(el_in,gp_in+1:gp_in+nqpte)
       Qyfae(ed)%ptr => Qyflux(el_ex,gp_ex+nqpte:gp_in+1:-1)         
     ENDDO

     DO ed = 1,nied

      ged = iedn(ed)

      inx(ed) = normal(1,ged)
      iny(ed) = normal(2,ged)

      len_area_in(ed) = edlen_area(1,ged)
      len_area_ex(ed) = edlen_area(2,ged)   

     ENDDO

     DO ged = 1,ned

       el_in = ged2el(1,ged)
       el_ex = ged2el(2,ged)

       led_in = ged2led(1,ged)
       led_ex = ged2led(2,ged)

       nxv(el_in,led_in) = normal(1,ged)
       nyv(el_in,led_in) = normal(2,ged)

       IF(el_ex.gt.0 .and. led_ex.gt.0) THEN
         nxv(el_ex,led_ex) = normal(1,ged)
         nyv(el_ex,led_ex) = normal(2,ged)
       ENDIF
       
     ENDDO
          

     RETURN
     END SUBROUTINE ptr_arrays
