       DO sp = 1,nsp

ed_points: DO pt = 1,nqpte

    l_edge: DO led = 1,3
              ind = (led-1)*nqpte + pt
!DIR$ VECTOR ALIGNED               
!               DO el = 1,ne
              DO el = split(1,sp),split(2,sp)
                Hqpt(el,ind) = H(el,1)
                Qxqpt(el,ind) = Qx(el,1)
                Qyqpt(el,ind) = Qy(el,1)
              ENDDO

    ed_basis: DO dof = 2,ndof     
!DIR$ VECTOR ALIGNED
!                 DO el = 1,ne ! Compute solutions at edge quadrature points
                DO el = split(1,sp),split(2,sp)
                  Hqpt(el,ind) = Hqpt(el,ind) + H(el,dof)*phie(dof,ind)            
                  Qxqpt(el,ind) = Qxqpt(el,ind) + Qx(el,dof)*phie(dof,ind)            
                  Qyqpt(el,ind) = Qyqpt(el,ind) + Qy(el,dof)*phie(dof,ind)            
                ENDDO

              ENDDO ed_basis

!               DO el = 1,ne ! Compute momentum terms
!DIR$ VECTOR ALIGNED
!!DIR$ SIMD
              DO el = split(1,sp),split(2,sp)
!               pressa(el) = .5d0*g*Hqpt(el,ind)*Hqpt(el,ind)
                recipHa(el) = 1d0/Hqpt(el,ind)
                

                xmom(el,ind) = pt5g*Hqpt(el,ind)*Hqpt(el,ind) + Qxqpt(el,ind)*Qxqpt(el,ind)*recipHa(el) 
                ymom(el,ind) = pt5g*Hqpt(el,ind)*Hqpt(el,ind) + Qyqpt(el,ind)*Qyqpt(el,ind)*recipHa(el) 
                xymom(el,ind) = Qxqpt(el,ind)*Qyqpt(el,ind)*recipHa(el)

              ENDDO
! !DIR$ IVDEP
!               DO el = split(1,sp),split(2,sp)
!                 Hn(el,ind) = nxv(el,led)*Qxqpt(el,ind) + nyv(el,led)*Qyqpt(el,ind)
!                 Qxn(el,ind) = nxv(el,led)*xmom(el,ind) + nyv(el,led)*xymom(el,ind)
!                 Qyn(el,ind) = nxv(el,led)*xymom(el,ind) + nyv(el,led)*ymom(el,ind)
!               ENDDO
! 
!               DO el = split(1,sp),split(2,sp)
!                 egnval(el,ind) = ABS(Hn(el,ind))*recipHa(el) + SQRT(g*Hqpt(el,ind))
!               ENDDO

            ENDDO l_edge

          ENDDO ed_points

       ENDDO
     
     
          DO sp = 1,nsp2
ed_points2: DO pt = 1,nqpte ! Compute numerical fluxes for all edges

!               ! Interior edges
!               DO ed = split2(1,sp),split2(2,sp)
!                 const(ed) = max(EVi(ed,pt)%ptr, EVe(ed,pt)%ptr)
! 
!                 Hhatv(ed) = .5d0*(Hni(ed,pt)%ptr + Hne(ed,pt)%ptr &
!                                         - const(ed)*(He(ed,pt)%ptr - Hi(ed,pt)%ptr))
!                 Qxhatv(ed) = .5d0*(Qxni(ed,pt)%ptr + Qxne(ed,pt)%ptr  &
!                                         - const(ed)*(Qxe(ed,pt)%ptr - Qxi(ed,pt)%ptr))
!                 Qyhatv(ed) = .5d0*(Qyni(ed,pt)%ptr + Qyne(ed,pt)%ptr  &
!                                         - const(ed)*(Qye(ed,pt)%ptr - Qyi(ed,pt)%ptr))
!               ENDDO

!DIR$ VECTOR ALIGNED
              DO ed = split2(1,sp),split2(2,sp)
                const(ed) = max(abs(Qxi(ed,pt)%ptr*inx(ed) + Qyi(ed,pt)%ptr*iny(ed))/Hi(ed,pt)%ptr + sqrt(g*Hi(ed,pt)%ptr), &
                                abs(Qxe(ed,pt)%ptr*inx(ed) + Qye(ed,pt)%ptr*iny(ed))/He(ed,pt)%ptr + sqrt(g*He(ed,pt)%ptr))
              ENDDO
!DIR$ IVDEP
              DO ed = split2(1,sp),split2(2,sp)
                Hhatv(ed) = .5d0*(inx(ed)*(Qxi(ed,pt)%ptr + Qxe(ed,pt)%ptr) + iny(ed)*(Qyi(ed,pt)%ptr + Qye(ed,pt)%ptr) &
                                        - const(ed)*(He(ed,pt)%ptr - Hi(ed,pt)%ptr))
                                        
                Hfe(ed,pt)%ptr = -len_area_ex(ed)*Hhatv(ed)
                Hfi(ed,pt)%ptr =  len_area_in(ed)*Hhatv(ed)
              ENDDO      
!DIR$ IVDEP
              DO ed = split2(1,sp),split2(2,sp)
                Qxhatv(ed) = .5d0*(inx(ed)*(xmi(ed,pt)%ptr + xme(ed,pt)%ptr) + iny(ed)*(xymi(ed,pt)%ptr + xyme(ed,pt)%ptr)  &
                                        - const(ed)*(Qxe(ed,pt)%ptr - Qxi(ed,pt)%ptr))
                                        
                Qxfe(ed,pt)%ptr = -len_area_ex(ed)*Qxhatv(ed)
                Qxfi(ed,pt)%ptr =  len_area_in(ed)*Qxhatv(ed)
              ENDDO   
!DIR$ IVDEP
              DO ed = split2(1,sp),split2(2,sp)
                Qyhatv(ed) = .5d0*(inx(ed)*(xymi(ed,pt)%ptr + xyme(ed,pt)%ptr) + iny(ed)*(ymi(ed,pt)%ptr + yme(ed,pt)%ptr)  &
                                        - const(ed)*(Qye(ed,pt)%ptr - Qyi(ed,pt)%ptr))
                Qyfe(ed,pt)%ptr = -len_area_ex(ed)*Qyhatv(ed)
                Qyfi(ed,pt)%ptr =  len_area_in(ed)*Qyhatv(ed)
              ENDDO

!               DO ed = 1,nied
!                 const(ed) = max(abs(Qxp(ed,pt)%in*inx(ed) + Qyp(ed,pt)%in*iny(ed))/Hp(ed,pt)%in + sqrt(g*Hp(ed,pt)%in), &
!                                 abs(Qxp(ed,pt)%ex*inx(ed) + Qyp(ed,pt)%ex*iny(ed))/Hp(ed,pt)%ex + sqrt(g*Hp(ed,pt)%ex))
! 
! 
! 
!                 Hhatv(ed) = .5d0*(inx(ed)*(Qxp(ed,pt)%in + Qxp(ed,pt)%ex) + iny(ed)*(Qyp(ed,pt)%in + Qyp(ed,pt)%ex) &
!                                         - const(ed)*(Hp(ed,pt)%ex - Hp(ed,pt)%in))
!                 Qxhatv(ed) = .5d0*(inx(ed)*(xmp(ed,pt)%in + xmp(ed,pt)%ex) + iny(ed)*(xymp(ed,pt)%in + xymp(ed,pt)%ex)  &
!                                         - const(ed)*(Qxp(ed,pt)%ex - Qxp(ed,pt)%in))
!                 Qyhatv(ed) = .5d0*(inx(ed)*(xymp(ed,pt)%in + xymp(ed,pt)%ex) + iny(ed)*(ymp(ed,pt)%in + ymp(ed,pt)%ex)  &
!                                         - const(ed)*(Qyp(ed,pt)%ex - Qyp(ed,pt)%in))
!               ENDDO


! ! !DIR$ SIMD
! !DIR$ IVDEP
!               DO ed = split2(1,sp),split2(2,sp)
!                 Hfe(ed,pt)%ptr = -len_area_ex(ed)*Hhatv(ed)
!                 Hfi(ed,pt)%ptr =  len_area_in(ed)*Hhatv(ed)
!               ENDDO
! ! !DIR$ SIMD
! !DIR$ IVDEP
!               DO ed = split2(1,sp),split2(2,sp)
!                 Qxfe(ed,pt)%ptr = -len_area_ex(ed)*Qxhatv(ed)
!                 Qxfi(ed,pt)%ptr =  len_area_in(ed)*Qxhatv(ed)
!               ENDDO
! ! !DIR$ SIMD
! !DIR$ IVDEP
!               DO ed = split2(1,sp),split2(2,sp)
!                 Qyfe(ed,pt)%ptr = -len_area_ex(ed)*Qyhatv(ed)
!                 Qyfi(ed,pt)%ptr =  len_area_in(ed)*Qyhatv(ed)
!               ENDDO

!DIR$ VECTOR ALIGNED
            DO ed = split2(1,sp),split2(2,sp)
              Hin(ed) = Hi(ed,pt)%ptr
              Qxin(ed) = Qxi(ed,pt)%ptr
              Qyin(ed) = Qyi(ed,pt)%ptr
              xmin(ed) = xmi(ed,pt)%ptr
              ymin(ed) = ymi(ed,pt)%ptr
              xymin(ed) = xymi(ed,pt)%ptr

              Hex(ed) = He(ed,pt)%ptr
              Qxex(ed) = Qxe(ed,pt)%ptr
              Qyex(ed) = Qye(ed,pt)%ptr
              xmex(ed) = xme(ed,pt)%ptr
              ymex(ed) = yme(ed,pt)%ptr
              xymex(ed) = xyme(ed,pt)%ptr
            ENDDO

              ! Interior edges
!DIR$ VECTOR ALIGNED              
              DO ed = split2(1,sp),split2(2,sp)
                const(ed) = max(abs(Qxin(ed)*inx(ed) + Qyin(ed)*iny(ed))/Hin(ed) + sqrt(g*Hin(ed)), &
                                abs(Qxex(ed)*inx(ed) + Qyex(ed)*iny(ed))/Hex(ed) + sqrt(g*Hex(ed)))
              ENDDO
!DIR$ VECTOR ALIGNED              
              DO ed = split2(1,sp),split2(2,sp)     
                Hhatv(ed) = .5d0*(inx(ed)*(Qxin(ed) + Qxex(ed)) + iny(ed)*(Qyin(ed) + Qyex(ed)) &
                                        - const(ed)*(Hex(ed) - Hin(ed)))
              ENDDO 
!DIR$ VECTOR ALIGNED              
              DO ed = split2(1,sp),split2(2,sp)              
                Qxhatv(ed) = .5d0*(inx(ed)*(xmin(ed) + xmex(ed)) + iny(ed)*(xymin(ed) + xymex(ed))  &
                                        - const(ed)*(Qxex(ed) - Qxin(ed)))
              ENDDO 
!DIR$ VECTOR ALIGNED                                         
              DO ed = split2(1,sp),split2(2,sp)              
                Qyhatv(ed) = .5d0*(inx(ed)*(xymin(ed) + xymex(ed)) + iny(ed)*(ymin(ed) + ymex(ed))  &
                                        - const(ed)*(Qyex(ed) - Qyin(ed)))
              ENDDO

!DIR$ IVDEP
!DIR$ VECTOR ALIGNED
              DO ed = split2(1,sp),split2(2,sp)
                Hfe(ed,pt)%ptr = -len_area_ex(ed)*Hhatv(ed)
                Hfi(ed,pt)%ptr =  len_area_in(ed)*Hhatv(ed)
              ENDDO
!DIR$ IVDEP
!DIR$ VECTOR ALIGNED
              DO ed = split2(1,sp),split2(2,sp)
                Qxfe(ed,pt)%ptr = -len_area_ex(ed)*Qxhatv(ed)
                Qxfi(ed,pt)%ptr =  len_area_in(ed)*Qxhatv(ed)
              ENDDO
!DIR$ IVDEP
!DIR$ VECTOR ALIGNED
              DO ed = split2(1,sp),split2(2,sp)
                Qyfe(ed,pt)%ptr = -len_area_ex(ed)*Qyhatv(ed)
                Qyfi(ed,pt)%ptr =  len_area_in(ed)*Qyhatv(ed)
              ENDDO
              
! !DIR$ IVDEP
! !DIR$ VECTOR ALIGNED
!               DO ed = nfblk(1,blk),nfblk(2,blk)             
!                 Hhatv(ed) = .5d0*(inx(ed)*(Qxi(ed,pt)%ptr + Qxe(ed,pt)%ptr) + iny(ed)*(Qyi(ed,pt)%ptr + Qye(ed,pt)%ptr) &
!                                         - const(ed)*(He(ed,pt)%ptr - Hi(ed,pt)%ptr))
! 
!                 Hfe(ed,pt)%ptr = -len_area_ex(ed)*Hhatv(ed)
!                 Hfi(ed,pt)%ptr =  len_area_in(ed)*Hhatv(ed)
!               ENDDO    
!               
!               DO ed = nfblk(1,blk),nfblk(2,blk)
!                 rHi(ed) = 1d0/Hi(ed,pt)%ptr
!                 rHe(ed) = 1d0/He(ed,pt)%ptr
!               
!                 xmomi(ed) = pt5g*Hi(ed,pt)%ptr*Hi(ed,pt)%ptr + Qxi(ed,pt)%ptr*Qxi(ed,pt)%ptr*rHi(ed)
!                 xmome(ed) = pt5g*He(ed,pt)%ptr*He(ed,pt)%ptr + Qxe(ed,pt)%ptr*Qxe(ed,pt)%ptr*rHe(ed)
!                 
!                 ymomi(ed) = pt5g*Hi(ed,pt)%ptr*Hi(ed,pt)%ptr + Qyi(ed,pt)%ptr*Qyi(ed,pt)%ptr*rHi(ed)
!                 ymome(ed) = pt5g*He(ed,pt)%ptr*He(ed,pt)%ptr + Qye(ed,pt)%ptr*Qye(ed,pt)%ptr*rHe(ed)
!                 
!                 xymomi(ed) = Qxi(ed,pt)%ptr*Qyi(ed,pt)%ptr*rHi(ed)
!                 xymome(ed) = Qxe(ed,pt)%ptr*Qye(ed,pt)%ptr*rHe(ed)                
!               ENDDO
!               
!               
! !DIR$ IVDEP
! !DIR$ VECTOR ALIGNED
!               DO ed = nfblk(1,blk),nfblk(2,blk)          
!                 Qxhatv(ed) = .5d0*(inx(ed)*(xmomi(ed) + xmome(ed)) + iny(ed)*(xymomi(ed) + xymome(ed))  &
!                                         - const(ed)*(Qxe(ed,pt)%ptr - Qxi(ed,pt)%ptr))
!                                  
!                 Qxfe(ed,pt)%ptr = -len_area_ex(ed)*Qxhatv(ed)
!                 Qxfi(ed,pt)%ptr =  len_area_in(ed)*Qxhatv(ed)
!               ENDDO   
!               
!               
!               
! !DIR$ IVDEP
! !DIR$ VECTOR ALIGNED
!               DO ed = nfblk(1,blk),nfblk(2,blk)
!                 Qyhatv(ed) = .5d0*(inx(ed)*(xymomi(ed) + xymome(ed)) + iny(ed)*(ymomi(ed) + ymome(ed))  &
!                                         - const(ed)*(Qye(ed,pt)%ptr - Qyi(ed,pt)%ptr))
!                                      
!                 Qyfe(ed,pt)%ptr = -len_area_ex(ed)*Qyhatv(ed)
!                 Qyfi(ed,pt)%ptr =  len_area_in(ed)*Qyhatv(ed)
!               ENDDO                     

     ENDDO ed_points2
   ENDDO