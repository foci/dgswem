!***********************************************************************
!     
!     subroutine p_enrichment()
!     
!     In concert, was born and evolved -- cem
!     
!***********************************************************************
      subroutine p_enrichment(s,it,irkp)
      
!.....use appropriate modules

      use sizes
      use global
      use dg

#ifdef CMPI
      USE MESSENGER_ELEM
#endif

      implicit none

      type (sizes_type) :: s

!.....declare local variables

      integer l, it,k,i,mm,j,irkp,m,ll
      real(sz) x_center, y_center, x_mid(3), y_mid(3)
      real(sz) ze_center, qx_center, qy_center
      real(sz) ze_sensor1, qx_sensor1, qy_sensor1
      real(sz) ze_sensor2, qx_sensor2, qy_sensor2
      real(sz) ze_mid, qx_mid, qy_mid,lebesgue,roo,kons
      real(sz) kon
      real(sz) maxze_sensor,maxqx_sensor,maxqy_sensor
      real(sz) deltai,dio
      real(sz) minze_sensor,minqx_sensor,minqy_sensor
      real(sz) ze_delta,qx_delta,qy_delta
      real(sz) avg_zesen,avg_qxsen,avg_qysen
      real(sz) temp_ze,temp_qx,temp_qy
      real(sz) ze_sensor(s%mne) 
      real(sz) qx_sensor(s%mne)
      real(sz) qy_sensor(s%mne)
#ifdef TRACE
      real(sz) iota_sensor1, iota_sensor2, iota_center, iota_mid
      real(sz) maxiota_sensor, miniota_sensor, iota_delta
      real(sz) avg_iotasen, temp_iota, iota_sensor(s%mne)
#endif
#ifdef CHEM
      real(sz) iota_sensor1, iota_sensor2, iota_center, iota_mid
      real(sz) maxiota_sensor, miniota_sensor, iota_delta
      real(sz) avg_iotasen, temp_iota, iota_sensor(s%mne)
      real(sz) iota2_sensor1, iota2_sensor2, iota2_center 
      real(sz) iota2_mid, slimit5, maxiota2_sensor, miniota2_sensor
      ! slimit5 is part of module DG but also being declared as a local variable
      real(sz) iota2_delta, avg_iota2sen, temp_iota2, iota2_sensor(s%mne)
      
#endif
#ifdef SED_LAY
      real(sz) tbed_sensor1, tbed_sensor2
      real(sz) tbed_center, tbed_mid
      real(sz) mintbed_sensor,maxtbed_sensor
      real(sz) avg_tbedsen, temp_tbed
      real(sz) tbed_sensor(s%mne)
      real(sz) slimit6, tbed_delta
#endif


!......Set the initial arrays, note we work over the total dg%bed load
!......not the individual layers

      if (dg%pflag.eq.1) then
         
         maxze_sensor = -100.d0 
         maxqx_sensor = -100.d0
         maxqy_sensor = -100.d0
         
         minze_sensor = 100.d0 
         minqx_sensor = 100.d0
         minqy_sensor = 100.d0
         
#ifdef TRACE
         maxiota_sensor = -100.d0
         miniota_sensor = 100.d0
#endif
         
#ifdef CHEM   
         maxiota_sensor = -100.d0
         maxiota2_sensor = -100.d0
         
         miniota_sensor = 100.d0
         miniota2_sensor = 100.d0
#endif

#ifdef SED_LAY
         maxtbed_sensor = -100.d0
         mintbed_sensor = 100.d0
#endif
         
         do l = 1,ne
            
            dg%pa = pdg_el(l)
            
!..........................
!.... type 1 p-enrichment |
!..........................
            
            ze_sensor(l) = -100.d0 
            qx_sensor(l) = -100.d0
            qy_sensor(l) = -100.d0
            
            ze_center = 0.d0
            qx_center = 0.d0
            qy_center = 0.d0
            
#ifdef TRACE
            iota_sensor(l) = -100.d0
            iota_center = 0.d0
#endif

#ifdef CHEM            
            iota_sensor(l) = -100.d0
            iota2_sensor(l) = -100.d0
            iota_center = 0.d0
            iota2_center = 0.d0
#endif
      
#ifdef SED_LAY
            tbed_sensor(l) = -100.d0
            tbed_center = 0.d0
#endif
            
            n1 = nm(l,1)
            n2 = nm(l,2)
            n3 = nm(l,3)
            
            x_center = 1.d0/3.d0*(x(n1) + x(n2) + x(n3))
            y_center = 1.d0/3.d0*(y(n1) + y(n2) + y(n3))
            
            x_mid(1) = 1.d0/2.d0*(x(n2) + x(n3))
            x_mid(2) = 1.d0/2.d0*(x(n3) + x(n1))
            x_mid(3) = 1.d0/2.d0*(x(n1) + x(n2))
            
            y_mid(1) = 1.d0/2.d0*(y(n2) + y(n3))
            y_mid(2) = 1.d0/2.d0*(y(n3) + y(n1))
            y_mid(3) = 1.d0/2.d0*(y(n1) + y(n2))
            
            if (dg%pa.ge.1) then   ! since the dg%pa=0 case is vacuous
               
               do k = 1,dg%dofs(l)
                  
                  ze_center = ze_center + dg%ze(k,l,dg%irk+1)*dg%phi_center(k,dg%pa)
                  qx_center = qx_center + dg%qx(k,l,dg%irk+1)*dg%phi_center(k,dg%pa)
                  qy_center = qy_center + dg%qy(k,l,dg%irk+1)*dg%phi_center(k,dg%pa)
                  
#ifdef TRACE
                  iota_center = iota_center + dg%iota(k,l,dg%irk+1)*dg%phi_center(k,dg%pa)
#endif
                  
#ifdef CHEM            
                  iota_center = iota_center + dg%iota(k,l,dg%irk+1)*dg%phi_center(k,dg%pa)
                  iota2_center = iota2_center + dg%iota2(k,l,dg%irk+1)*dg%phi_center(k,dg%pa)
#endif

#ifdef SED_LAY
                  do ll=1,layers
                     tbed_center = tbed_center + dg%bed(k,l,dg%irk+1,ll)*dg%phi_center(k,dg%pa)
                  enddo
#endif
                  
               enddo
               
            endif
            
!......Compute the sensor using edges in each variable
            
            do i = 1,3
               
               ze_mid = 0.d0
               qx_mid = 0.d0
               qy_mid = 0.d0
               
#ifdef TRACE
               iota_mid = 0.d0
#endif

#ifdef CHEM               
               iota_mid = 0.d0
               iota2_mid = 0.d0
#endif

#ifdef SED_LAY
               tbed_mid = 0.d0
#endif
               
               dist = sqrt( (x_mid(i) - x_center)**2.d0 +&
              (y_mid(i) - y_center)**2.d0 )
               
               
               if (dg%pa.ge.1) then
                  
                  do k = 1,dg%dofs(l)
                     
                     ze_mid = ze_mid + dg%ze(k,l,dg%irk+1)*dg%phi_mid(k,i,dg%pa)
                     qx_mid = qx_mid + dg%qx(k,l,dg%irk+1)*dg%phi_mid(k,i,dg%pa)
                     qy_mid = qy_mid + dg%qy(k,l,dg%irk+1)*dg%phi_mid(k,i,dg%pa)
                     
#ifdef TRACE
                     iota_mid = iota_mid + dg%iota(k,l,dg%irk+1)*dg%phi_mid(k,i,dg%pa)
#endif
                     
#ifdef CHEM               
                     iota_mid = iota_mid + dg%iota(k,l,dg%irk+1)*dg%phi_mid(k,i,dg%pa)
                     iota2_mid = iota2_mid + dg%iota2(k,l,dg%irk+1)*dg%phi_mid(k,i,dg%pa)
#endif

#ifdef SED_LAY
                     do ll=1,layers
                        tbed_mid = tbed_mid + dg%bed(k,l,dg%irk+1,ll)*dg%phi_mid(k,i,dg%pa)
                     enddo
#endif
                     
                  enddo
                  
                  ze_sensor(l) = max(abs( (ze_mid - ze_center)/dist ),ze_sensor(l))
                  qx_sensor(l) = max(abs( (qx_mid - qx_center)/dist ),qx_sensor(l))
                  qy_sensor(l) = max(abs( (qy_mid - qy_center)/dist ),qy_sensor(l))
                  
#ifdef TRACE
                  iota_sensor(l) = max(abs( (iota_mid - iota_center)/dist ),iota_sensor(l))
#endif

#ifdef CHEM                  
                  iota_sensor(l) = max(abs( (iota_mid - iota_center)/dist ),iota_sensor(l))
                  iota2_sensor(l) = max(abs( (iota2_mid - iota2_center)/dist ),iota_sensor(l))
#endif

#ifdef SED_LAY
                  tbed_sensor(l) = max(abs( (tbed_mid - tbed_center)/dist ),tbed_sensor(l))
#endif
                  
               endif
               
            enddo
            
            
            if (dg%gflag.eq.0) then ! dioristic algorithm is OFF
               
#ifdef TRACE

#elif CHEM

#elif SED_LAY

#else
               
!.....if the sensor is greater than the limit and the p is low increase
!.....the order of p

               if ( ( (qy_sensor(l).ge.dg%slimit).or.&
              (qx_sensor(l).ge.dg%slimit).or.&
              (ze_sensor(l).ge.dg%slimit) ).and.&
              (pdg_el(l).lt.dg%ph) ) then
                  
                  pdg_el(l) = pdg_el(l) + 1 
                  dg%dofs(l) = (pdg_el(l) +1 )*(pdg_el(l) +2 ) / 2
                  
                  dg%pcount(l) = it
                  
!.....if the sensor is less than the limit and the p is high decrease
!.....the order of p
                  
               elseif ( ((ze_sensor(l).lt.dg%slimit).and.&
                 (qx_sensor(l).lt.dg%slimit).and.&
                 (qy_sensor(l).lt.dg%slimit) ).and.&
                 (pdg_el(l).gt.dg%pl).and.&
                 ( (it-dg%pcount(l)).gt.dg%plimit ) ) then
                  
                  pdg_el(l) = pdg_el(l) - 1
                  dg%dofs(l) = (pdg_el(l) + 1 )*(pdg_el(l) + 2 ) / 2
                  
                  dg%ze(dg%dofs(l)+1:dg%dofh,l,dg%irk+1) = 0.d0
                  dg%qx(dg%dofs(l)+1:dg%dofh,l,dg%irk+1) = 0.d0
                  dg%qy(dg%dofs(l)+1:dg%dofh,l,dg%irk+1) = 0.d0
                  
                  dg%pcount(l) = it
                  
                  
               elseif (pdg_el(l).eq.0.and.( (it-dg%pcount(l)).gt.dg%plimit ) ) then
                  
                  pdg_el(l) = pdg_el(l) + 1 
                  dg%dofs(l) = (pdg_el(l) + 1 )*(pdg_el(l) + 2 ) / 2
                  
                  dg%pcount(l) = it
                  
               endif

#endif

#ifdef TRACE       
               
!.....if the sensor is greater than the limit and the p is low increase
!.....the order of p

               if ( ( (qy_sensor(l).ge.dg%slimit).or.&
              (qx_sensor(l).ge.dg%slimit).or.&
              (ze_sensor(l).ge.dg%slimit).or.&
              ( iota_sensor(l).ge.dg%slimit) ).and.&
              (pdg_el(l).lt.dg%ph).and.&
              ( (it-dg%pcount(l)).gt.dg%plimit ) ) then
                  
                  pdg_el(l) = pdg_el(l) + 1 
                  dg%dofs(l) = (pdg_el(l) +1 )*(pdg_el(l) +2 ) / 2
                  
                  dg%pcount(l) = it
                  
                  
!.....if the sensor is less than the limit and the p is high decrease
!.....the order of p
                  
                  
               elseif ( ((ze_sensor(l).lt.dg%slimit).and.&
                 (qx_sensor(l).lt.dg%slimit).and.&
                 (qy_sensor(l).lt.dg%slimit).and.&
                 (iota_sensor(l).lt.dg%slimit) ).and.&
                 (pdg_el(l).gt.dg%pl).and.&
                 ( (it-dg%pcount(l)).gt.dg%plimit ) ) then
                  
                  pdg_el(l) = pdg_el(l) - 1
                  dg%dofs(l) = (pdg_el(l) + 1 )*(pdg_el(l) + 2 ) / 2
                  
                  dg%ze(dg%dofs(l)+1:dg%dofh,l,dg%irk+1) = 0.d0
                  dg%iota(dg%dofs(l)+1:dg%dofh,l,dg%irk+1) = 0.d0
                  dg%qx(dg%dofs(l)+1:dg%dofh,l,dg%irk+1) = 0.d0
                  dg%qy(dg%dofs(l)+1:dg%dofh,l,dg%irk+1) = 0.d0
                  
                  dg%pcount(l) = it
                  
                  
               elseif (pdg_el(l).eq.0.and.( (it-dg%pcount(l)).gt.dg%plimit ) ) then
                  
                  pdg_el(l) = pdg_el(l) + 1 
                  dg%dofs(l) = (pdg_el(l) + 1 )*(pdg_el(l) + 2 ) / 2
                  
                  dg%pcount(l) = it
                  
               endif

#endif
               
#ifdef CHEM       
               
!.....if the sensor is greater than the limit and the p is low increase
!.....the order of p
               
               if ( ( (qy_sensor(l).ge.dg%slimit).or.&
              (qx_sensor(l).ge.dg%slimit).or.&
              (ze_sensor(l).ge.dg%slimit).or.&
              (iota_sensor(l).ge.dg%slimit).or.&
              (iota2_sensor(l).ge.dg%slimit) ).and.&
              (pdg_el(l).lt.dg%ph) ) then
                  
                  pdg_el(l) = pdg_el(l) + 1 
                  dg%dofs(l) = (pdg_el(l) +1 )*(pdg_el(l) +2 ) / 2
                  
                  dg%pcount(l) = it
                  
                  
!.....if the sensor is less than the limit and the p is high decrease
!.....the order of p
                  
               elseif ( ((ze_sensor(l).lt.dg%slimit).and.&
                 (qx_sensor(l).lt.dg%slimit).and.&
                 (qy_sensor(l).lt.dg%slimit).and.&
                 (iota_sensor(l).lt.dg%slimit).and.&
                 (iota2_sensor(l).lt.dg%slimit)).and.&
                 (pdg_el(l).gt.dg%pl).and.&
                 ( (it-dg%pcount(l)).gt.dg%plimit ) ) then
                  
                  pdg_el(l) = pdg_el(l) - 1
                  dg%dofs(l) = (pdg_el(l) + 1 )*(pdg_el(l) + 2 ) / 2
                  
                  dg%ze(dg%dofs(l)+1:dg%dofh,l,dg%irk+1) = 0.d0
                  dg%iota(dg%dofs(l)+1:dg%dofh,l,dg%irk+1) = 0.d0
                  dg%iota2(dg%dofs(l)+1:dg%dofh,l,dg%irk+1) = 0.d0
                  dg%qx(dg%dofs(l)+1:dg%dofh,l,dg%irk+1) = 0.d0
                  dg%qy(dg%dofs(l)+1:dg%dofh,l,dg%irk+1) = 0.d0
                  
                  dg%pcount(l) = it
                  
                  
               elseif (pdg_el(l).eq.0.and.( (it-dg%pcount(l)).gt.dg%plimit ) ) then
                  
                  pdg_el(l) = pdg_el(l) + 1 
                  dg%dofs(l) = (pdg_el(l) + 1 )*(pdg_el(l) + 2 ) / 2
                  
                  dg%pcount(l) = it
                  
               endif
#endif

!when chem or tracers are on, this coupling could be made tighter      
#ifdef SED_LAY 
               
!.....if the sensor is greater than the limit and the p is low increase
!.....the order of p

               if ( ( (qy_sensor(l).ge.dg%slimit).or.&
              (qx_sensor(l).ge.dg%slimit).or.&
              (ze_sensor(l).ge.dg%slimit).or.&
              (tbed_sensor(l).ge.dg%slimit) ).and.&
              (pdg_el(l).lt.dg%ph).and.&
              ( (it-dg%pcount(l)).gt.dg%plimit ) ) then
                  
                  pdg_el(l) = pdg_el(l) + 1 
                  dg%dofs(l) = (pdg_el(l) +1 )*(pdg_el(l) +2 ) / 2
                  
                  dg%pcount(l) = it
                  
                  
!.....if the sensor is less than the limit and the p is high decrease
!.....the order of p
                  
                  
               elseif ( ((ze_sensor(l).lt.dg%slimit).and.&
                 (qx_sensor(l).lt.dg%slimit).and.&
                 (qy_sensor(l).lt.dg%slimit).and.&
                 (tbed_sensor(l).lt.dg%slimit) ).and.&
                 (pdg_el(l).gt.dg%pl).and.&
                 ( (it-dg%pcount(l)).gt.dg%plimit ) ) then
                  
                  pdg_el(l) = pdg_el(l) - 1
                  dg%dofs(l) = (pdg_el(l) + 1 )*(pdg_el(l) + 2 ) / 2
                  
                  dg%ze(dg%dofs(l)+1:dg%dofh,l,dg%irk+1) = 0.d0
                  dg%iota(dg%dofs(l)+1:dg%dofh,l,dg%irk+1) = 0.d0
                  dg%qx(dg%dofs(l)+1:dg%dofh,l,dg%irk+1) = 0.d0
                  dg%qy(dg%dofs(l)+1:dg%dofh,l,dg%irk+1) = 0.d0
                  dg%bed(dg%dofs(l)+1:dg%dofh,l,dg%irk+1,:)=0.d0
                  
                  dg%pcount(l) = it
                  
                  
               elseif (pdg_el(l).eq.0.and.( (it-dg%pcount(l)).gt.dg%plimit ) ) then
                  
                  pdg_el(l) = pdg_el(l) + 1 
                  dg%dofs(l) = (pdg_el(l) + 1 )*(pdg_el(l) + 2 ) / 2
                  
                  dg%pcount(l) = it
                  
               endif

#endif
               
            endif

         enddo

         if (dg%gflag.eq.1) then   ! dioristic algorithm is ON

            avg_zesen = 0.d0
            avg_qxsen = 0.d0
            avg_qysen = 0.d0

#ifdef TRACE
            avg_iotasen = 0.d0
#endif

#ifdef CHEM
            avg_iotasen = 0.d0
            avg_iota2sen = 0.d0
#endif

#ifdef SED_LAY
            avg_tbedsen = 0.d0
#endif
            
            do l=1,ne
               
#ifdef CMPI

               if (reselem(l)) then
                  
#endif

                  if (ze_sensor(l).ne.0.d0) then

                     minze_sensor = min(ze_sensor(l),minze_sensor)
                     maxze_sensor = max(ze_sensor(l),maxze_sensor)
                     
                  endif

                  if (qx_sensor(l).ne.0.d0) then

                     minqx_sensor = min(qx_sensor(l),minqx_sensor)
                     maxqx_sensor = max(qx_sensor(l),maxqx_sensor)

                  endif

                  if (qy_sensor(l).ne.0.d0) then

                     minqy_sensor = min(qy_sensor(l),minqy_sensor)
                     maxqy_sensor = max(qy_sensor(l),maxqy_sensor)
                     
                  endif

#ifdef TRACE
                  if (iota_sensor(l).ne.0.d0) then

                     miniota_sensor = min(iota_sensor(l),miniota_sensor)
                     maxiota_sensor = max(iota_sensor(l),maxiota_sensor)

                  endif
#endif

#ifdef CHEM
                  if (iota_sensor(l).ne.0.d0) then

                     miniota_sensor = min(iota_sensor(l),miniota_sensor)
                     maxiota_sensor = max(iota_sensor(l),maxiota_sensor)

                  endif

                  if (iota2_sensor(l).ne.0.d0) then
                     
                     miniota2_sensor = min(iota2_sensor(l),miniota2_sensor)
                     maxiota2_sensor = max(iota2_sensor(l),maxiota2_sensor)
                     
                  endif
#endif

#ifdef SED_LAY
                  if (tbed_sensor(l).ne.0.d0) then

                     mintbed_sensor = min(tbed_sensor(l),mintbed_sensor)
                     maxtbed_sensor = max(tbed_sensor(l),maxtbed_sensor)

                  endif
#endif
                  
                  avg_zesen = sum(ze_sensor)/dg%MNES
                  avg_qxsen = sum(qx_sensor)/dg%MNES
                  avg_qysen = sum(qy_sensor)/dg%MNES
                  
#ifdef TRACE
                  avg_iotasen = sum(iota_sensor)/dg%MNES
#endif

#ifdef CHEM                  
                  avg_iotasen = sum(iota_sensor)/dg%MNES
                  avg_iota2sen = sum(iota2_sensor)/dg%MNES
#endif

#ifdef SED_LAY
                  avg_tbedsen = sum(tbed_sensor)/dg%MNES
#endif

#ifdef CMPI

               endif

#endif
               
            enddo

#ifdef CMPI
            CALL PARA_MAX(MAXZE_SENSOR)
            CALL PARA_MIN(MINZE_SENSOR)
            CALL PARA_SUM(AVG_ZESEN)
            CALL PARA_SUM(AVG_QXSEN)
            CALL PARA_MAX(MAXQX_SENSOR)
            CALL PARA_MIN(MINQX_SENSOR)
            CALL PARA_SUM(AVG_QYSEN)
            CALL PARA_MAX(MAXQY_SENSOR)
            CALL PARA_MIN(MINQY_SENSOR)
#endif
#ifdef CMPI_TRACE
            CALL PARA_MAX(MAXiota_SENSOR)
            CALL PARA_MIN(MINiota_SENSOR)
            CALL PARA_SUM(AVG_iotaSEN)
#endif

#ifdef CMPI_CHEM
            CALL PARA_MAX(MAXiota_SENSOR)
            CALL PARA_MIN(MINiota_SENSOR)
            CALL PARA_SUM(AVG_iotaSEN)
            CALL PARA_MAX(MAXiota2_SENSOR)
            CALL PARA_MIN(MINiota2_SENSOR)
            CALL PARA_SUM(AVG_iota2SEN)
#endif


#ifdef CMPI_SED_LAY
            CALL PARA_MAX(MAXtbed_SENSOR)
            CALL PARA_MIN(MINtbed_SENSOR)
            CALL PARA_SUM(AVG_tbedSEN)
#endif

            ze_delta = 0.d0
            qx_delta = 0.d0
            qy_delta = 0.d0

#ifdef TRACE
            iota_delta = 0.d0
#endif

#ifdef CHEM
            iota_delta = 0.d0
            iota2_delta = 0.d0
#endif

#ifdef SED_LAY
            tbed_delta = 0.d0
#endif
            

            ze_delta = ( dg%diorism / 100.d0 )* (abs(maxze_sensor - minze_sensor) )
            qx_delta = ( dg%diorism / 100.d0 )* (abs(maxqx_sensor - minqx_sensor) )
            qy_delta = ( dg%diorism / 100.d0 )* (abs(maxqy_sensor - minqy_sensor) )   
            
#ifdef TRACE
            iota_delta = ( dg%diorism / 100.d0 )* (abs(maxiota_sensor - miniota_sensor) )
#endif
            
#ifdef CHEM
            iota_delta = ( dg%diorism / 100.d0 )* (abs(maxiota_sensor - miniota_sensor) )
            iota2_delta = ( dg%diorism / 100.d0 )* (abs(maxiota2_sensor - miniota2_sensor) )
#endif

#ifdef SED_LAY    
            tbed_delta = ( dg%diorism / 100.d0 )* (abs(maxtbed_sensor - mintbed_sensor) )
#endif

            
            do l = 1,ne
               
               temp_ze = 0.d0  
               temp_qx = 0.d0
               temp_qy = 0.d0
               
#ifdef TRACE
               temp_iota = 0.d0
#endif

#ifdef CHEM
               temp_iota = 0.d0
               temp_iota2 = 0.d0
#endif

#ifdef SED_LAY
               temp_tbed = 0.d0
#endif

               temp_ze = abs(ze_sensor(l) - avg_zesen)
               temp_qx = abs(qx_sensor(l) - avg_qxsen)
               temp_qy = abs(qy_sensor(l) - avg_qysen)

#ifdef TRACE
               temp_iota = abs(iota_sensor(l) - avg_iotasen)
#endif

#ifdef CHEM
               temp_iota = abs(iota_sensor(l) - avg_iotasen)
               temp_iota2 = abs(iota2_sensor(l) - avg_iota2sen)
#endif

#ifdef SED_LAY
               temp_tbed = abs(tbed_sensor(l) - avg_tbedsen)
#endif
               
#ifdef TRACE

#elif CHEM

#elif SED_LAY

#else 
               
!.....if the sensor is more than the limit and the p is low increase
!.....the order of p
               
               if ( ( (qy_sensor(l).ge.qy_delta).or.&
              (qx_sensor(l).ge.qx_delta).or.&
              (ze_sensor(l).ge.ze_delta) ).and.&
              (pdg_el(l).lt.dg%ph) ) then
                  
                  pdg_el(l) = pdg_el(l) + 1 
                  dg%dofs(l) = (pdg_el(l) +1 )*(pdg_el(l) +2 ) / 2
                  
                  dg%pcount(l) = it
                  
!.....if the sensor is less than the limit and the p is high decrease
!.....the order of p
                  
               elseif ( ((ze_sensor(l).lt.ze_delta).and.&
                 (qx_sensor(l).lt.qx_delta).and.&
                 (qy_sensor(l).lt.qy_delta) ).and.&
                 (pdg_el(l).gt.dg%pl).and.&
                 ( (it-dg%pcount(l)).gt.dg%plimit ) ) then
                  
                  pdg_el(l) = pdg_el(l) - 1
                  dg%dofs(l) = (pdg_el(l) + 1 )*(pdg_el(l) + 2 ) / 2
                  
                  dg%ze(dg%dofs(l)+1:dg%dofh,l,dg%irk+1) = 0.d0
                  dg%qx(dg%dofs(l)+1:dg%dofh,l,dg%irk+1) = 0.d0
                  dg%qy(dg%dofs(l)+1:dg%dofh,l,dg%irk+1) = 0.d0
                  
                  dg%pcount(l) = it
                  
               elseif (pdg_el(l).eq.0.and.&
                 ( (it-dg%pcount(l)).gt.dg%plimit ) ) then
                  
                  pdg_el(l) = pdg_el(l) + 1 
                  dg%dofs(l) = (pdg_el(l) + 1 )*(pdg_el(l) + 2 ) / 2
                  
                  dg%pcount(l) = it
                  
               endif
#endif

#ifdef TRACE

!.....if the sensor is more than the limit and the p is low increase
!.....the order of p
               
               
               if ( ( (temp_qy.ge.qy_delta).or.&
              (temp_qx.ge.qx_delta).or.&
              (temp_ze.ge.ze_delta).or.&
              (temp_iota.ge.iota_delta)).and.&
              (pdg_el(l).lt.dg%ph).and.&
              ( (it-dg%pcount(l)).gt.dg%plimit ) ) then
                  
                  pdg_el(l) = pdg_el(l) + 1 
                  dg%dofs(l) = (pdg_el(l) +1 )*(pdg_el(l) +2 ) / 2
                  
                  dg%pcount(l) = it
                  
!.....if the sensor is less than the limit and the p is high decrease
!.....the order of p
                  
               elseif ( ((temp_ze.lt.ze_delta).and.&
                 (temp_qx.lt.qx_delta).and.&
                 (temp_qy.lt.qy_delta).and.&
                 (temp_iota.lt.iota_delta)).and.&
                 (pdg_el(l).gt.dg%pl).and.&
                 ( (it-dg%pcount(l)).gt.dg%plimit ) ) then
                  
                  pdg_el(l) = pdg_el(l) - 1
                  dg%dofs(l) = (pdg_el(l) + 1 )*(pdg_el(l) + 2 ) / 2
                  
                  dg%ze(dg%dofs(l)+1:dg%dofh,l,dg%irk+1) = 0.d0
                  dg%qx(dg%dofs(l)+1:dg%dofh,l,dg%irk+1) = 0.d0
                  dg%qy(dg%dofs(l)+1:dg%dofh,l,dg%irk+1) = 0.d0
                  dg%iota(dg%dofs(l)+1:dg%dofh,l,dg%irk+1) = 0.d0
                  
                  dg%pcount(l) = it
                  
               elseif (pdg_el(l).eq.0.and.&
                 ( (it-dg%pcount(l)).gt.dg%plimit ) ) then
                  
                  pdg_el(l) = pdg_el(l) + 1 
                  dg%dofs(l) = (pdg_el(l) + 1 )*(pdg_el(l) + 2 ) / 2
                  
                  dg%pcount(l) = it
                  
               endif
#endif
               
#ifdef CHEM

!.....if the sensor is more than the limit and the p is low increase
!.....the order of p

               if ( ( (qy_sensor(l).ge.qy_delta).or.&
              (qx_sensor(l).ge.qx_delta).or.&
              (ze_sensor(l).ge.ze_delta).or.&
              iota_sensor(l).ge.iota_delta.or.&
              iota2_sensor(l).ge.iota2_delta).and.&
              (pdg_el(l).lt.dg%ph) ) then
                  
                  pdg_el(l) = pdg_el(l) + 1 
                  dg%dofs(l) = (pdg_el(l) +1 )*(pdg_el(l) +2 ) / 2
                  
                  dg%pcount(l) = it
                  
!.....if the sensor is less than the limit and the p is high decrease
!.....the order of p
                  
               elseif ( ((ze_sensor(l).lt.ze_delta).and.&
                 (qx_sensor(l).lt.qx_delta).and.&
                 (qy_sensor(l).lt.qy_delta) ).and.&
                 (iota_sensor(l).lt.iota_delta).and.&
                 (iota2_sensor(l).lt.iota2_delta).and.&
                 (pdg_el(l).gt.dg%pl).and.&
                 ( (it-dg%pcount(l)).gt.dg%plimit ) ) then
                  
                  pdg_el(l) = pdg_el(l) - 1
                  dg%dofs(l) = (pdg_el(l) + 1 )*(pdg_el(l) + 2 ) / 2
                  
                  dg%ze(dg%dofs(l)+1:dg%dofh,l,dg%irk+1) = 0.d0
                  dg%qx(dg%dofs(l)+1:dg%dofh,l,dg%irk+1) = 0.d0
                  dg%qy(dg%dofs(l)+1:dg%dofh,l,dg%irk+1) = 0.d0
                  dg%iota(dg%dofs(l)+1:dg%dofh,l,dg%irk+1) = 0.d0
                  dg%iota2(dg%dofs(l)+1:dg%dofh,l,dg%irk+1) = 0.d0
                  
                  dg%pcount(l) = it
                  
               elseif (pdg_el(l).eq.0.and.&
                 ( (it-dg%pcount(l)).gt.dg%plimit ) ) then
                  
                  pdg_el(l) = pdg_el(l) + 1 
                  dg%dofs(l) = (pdg_el(l) + 1 )*(pdg_el(l) + 2 ) / 2
                  
                  dg%pcount(l) = it
                  
               endif
#endif


#ifdef SED_LAY

!.....if the sensor is more than the limit and the p is low increase
!.....the order of p
               
               
               if ( ( (temp_qy.ge.qy_delta).or.&
              (temp_qx.ge.qx_delta).or.&
              (temp_ze.ge.ze_delta).or.&
              (temp_tbed.ge.tbed_delta)).and.&
              (pdg_el(l).lt.dg%ph).and.&
              ( (it-dg%pcount(l)).gt.dg%plimit ) ) then
                  
                  pdg_el(l) = pdg_el(l) + 1 
                  dg%dofs(l) = (pdg_el(l) +1 )*(pdg_el(l) +2 ) / 2
                  
                  dg%pcount(l) = it
                  
!.....if the sensor is less than the limit and the p is high decrease
!.....the order of p
                  
               elseif ( ((temp_ze.lt.ze_delta).and.&
                 (temp_qx.lt.qx_delta).and.&
                 (temp_qy.lt.qy_delta).and.&
                 (temp_tbed.lt.tbed_delta)).and.&
                 (pdg_el(l).gt.dg%pl).and.&
                 ( (it-dg%pcount(l)).gt.dg%plimit ) ) then
                  
                  pdg_el(l) = pdg_el(l) - 1
                  dg%dofs(l) = (pdg_el(l) + 1 )*(pdg_el(l) + 2 ) / 2
                  
                  dg%ze(dg%dofs(l)+1:dg%dofh,l,dg%irk+1) = 0.d0
                  dg%qx(dg%dofs(l)+1:dg%dofh,l,dg%irk+1) = 0.d0
                  dg%qy(dg%dofs(l)+1:dg%dofh,l,dg%irk+1) = 0.d0
                  dg%iota(dg%dofs(l)+1:dg%dofh,l,dg%irk+1) = 0.d0
                  dg%bed(dg%dofs(l)+1:dg%dofh,l,dg%irk+1,:)=0.d0
                  
                  dg%pcount(l) = it
                  
               elseif (pdg_el(l).eq.0.and.&
                 ( (it-dg%pcount(l)).gt.dg%plimit ) ) then
                  
                  pdg_el(l) = pdg_el(l) + 1 
                  dg%dofs(l) = (pdg_el(l) + 1 )*(pdg_el(l) + 2 ) / 2
                  
                  dg%pcount(l) = it
                  
               endif
#endif
               
            enddo
            
         endif
         
      endif

!..........................
!.... type 2 p-enrichment |
!..........................
      
      if (dg%pflag.eq.2) then

         maxze_sensor = -100.d0 
         maxqx_sensor = -100.d0
         maxqy_sensor = -100.d0

         minze_sensor = 100.d0 
         minqx_sensor = 100.d0
         minqy_sensor = 100.d0

#ifdef TRACE
         maxiota_sensor = -100.d0
         miniota_sensor = 100.d0               
#endif

#ifdef CHEM
         maxiota_sensor = -100.d0
         maxiota2_sensor = -100.d0

         miniota_sensor = 100.d0
         miniota2_sensor = 100.d0
#endif

#ifdef SED_LAY
         maxtbed_sensor = -100.d0
         mintbed_sensor =  100.d0 
#endif

         
         kon = dg%pflag2con1       ! usually in [0,1]
         kons = dg%pflag2con2      ! usually 6 works
         lebesgue = dg%lebesguep
         roo = 1.d0/lebesgue      

         do l=1,ne

            ze_sensor(l) = 0.d0             
            qx_sensor(l) = 0.d0
            qy_sensor(l) = 0.d0
            
#ifdef TRACE
            iota_sensor(l) = 0.d0            
#endif
            
#ifdef CHEM      
            iota_sensor(l) = 0.d0
            iota2_sensor(l) = 0.d0
#endif

#ifdef SED_LA
!$$$            do ll=1,layers
!$$$               tbed_center = tbed_center + dg%bed(k,l,dg%irk+1,ll)*dg%phi_center(k,dg%pa)
!$$$            enddo
            tbed_sensor(l) = 0.d0            
#endif
            
            ze_sensor1 = 0.d0             
            qx_sensor1 = 0.d0
            qy_sensor1 = 0.d0
            
            ze_sensor2 = 0.d0            
            qx_sensor2 = 0.d0
            qy_sensor2 = 0.d0
            
            dg%slimit1 = 0.d0 
            dg%slimit2 = 0.d0
            dg%slimit3 = 0.d0

#ifdef TRACE
            iota_sensor1 = 0.d0

            iota_sensor2 = 0.d0
            
            dg%slimit4 = 0.d0
#endif
            
#ifdef CHEM      
            iota_sensor1 = 0.d0
            iota2_sensor1 = 0.d0
            
            iota_sensor2 = 0.d0
            iota2_sensor2 = 0.d0
            
            dg%slimit4 = 0.d0
            dg%slimit5 = 0.d0
#endif

#ifdef SED_LAY
            tbed_sensor1 = 0.d0

            tbed_sensor2 = 0.d0
            
            slimit6 = 0.d0
#endif

            dg%pa = pdg_el(l)

                                !compute the first sensor type          
            if (dg%pa.eq.0) then
               
               cycle
               
            else

               do k = (dg%pa*(dg%pa+1))/2 + 1 ,dg%dofs(l)
                  do mm = 1,dg%nagp(dg%pa)
                     
                     ze_sensor1 = ze_sensor1 + abs(dg%ze(k,l,irkp+1)* &
                    dg%phi_area(k,mm,dg%pa))**lebesgue * dg%wagp(mm,dg%pa)
                     qx_sensor1 = qx_sensor1 + abs(dg%qx(k,l,irkp+1)* &
                    dg%phi_area(k,mm,dg%pa))**lebesgue * dg%wagp(mm,dg%pa)
                     qy_sensor1 = qy_sensor1 + abs(dg%qy(k,l,irkp+1)* &
                    dg%phi_area(k,mm,dg%pa))**lebesgue * dg%wagp(mm,dg%pa)
                     
#ifdef TRACE       
                     iota_sensor1 = iota_sensor1 + abs(dg%iota(k,l,irkp+1)* &
                    dg%phi_area(k,mm,dg%pa))**lebesgue * dg%wagp(mm,dg%pa)
#endif
                     
#ifdef CHEM                
                     iota_sensor1  = iota_sensor1 + abs(dg%iota(k,l,irkp+1)* &
                    dg%phi_area(k,mm,dg%pa))**lebesgue * dg%wagp(mm,dg%pa)
                     iota2_sensor1 = iota2_sensor1 + abs(dg%iota2(k,l,irkp+1)* &
                    dg%phi_area(k,mm,dg%pa))**lebesgue * dg%wagp(mm,dg%pa)
#endif

#ifdef SED_LAY       

                     do ll=1,layers
                        tbed_sensor1 = tbed_sensor1 + abs(dg%bed(k,l,irkp+1,ll)* &
                       dg%phi_area(k,mm,dg%pa))**lebesgue * dg%wagp(mm,dg%pa)
                     enddo
#endif
                     
                  enddo
               enddo
               
            endif
            
                                !compute the second sensor type

            if (dg%pa.eq.0) then
               
               cycle
               
            else
               
               do k = 1, dg%dofs(l)

                  do mm=1,dg%nagp(dg%pa)
                     
                     ze_sensor2 = ze_sensor2 + abs(dg%ze(k,l,irkp+1)* &
                    dg%phi_area(k,mm,dg%pa))**lebesgue * dg%wagp(mm,dg%pa)
                     qx_sensor2 = qx_sensor2 + abs(dg%qx(k,l,irkp+1)* &
                    dg%phi_area(k,mm,dg%pa))**lebesgue * dg%wagp(mm,dg%pa)
                     qy_sensor2 = qy_sensor2 + abs(dg%qy(k,l,irkp+1)* &
                    dg%phi_area(k,mm,dg%pa))**lebesgue * dg%wagp(mm,dg%pa)
                     
#ifdef TRACE
                     iota_sensor2 = iota_sensor2 + abs(dg%iota(k,l,irkp+1)* &
                    dg%phi_area(k,mm,dg%pa))**lebesgue * dg%wagp(mm,dg%pa)
#endif
                     
#ifdef CHEM
                     iota_sensor2 = iota_sensor2 + abs(dg%iota(k,l,irkp+1)* &
                    dg%phi_area(k,mm,dg%pa))**lebesgue * dg%wagp(mm,dg%pa)
                     iota2_sensor2 = iota2_sensor2 + abs(dg%iota2(k,l,irkp+1)* &
                    dg%phi_area(k,mm,dg%pa))**lebesgue * dg%wagp(mm,dg%pa)
#endif

                     
#ifdef SED_LAY
                     do ll=1,layers
                        tbed_sensor2 = tbed_sensor2 + abs(dg%bed(k,l,irkp+1,ll)* &
                       dg%phi_area(k,mm,dg%pa))**lebesgue * dg%wagp(mm,dg%pa)
                     enddo
#endif
                     
                  enddo

               enddo
               
            endif

!................................................................

            if (dg%gflag.eq.0) then !dioristic algorithm OFF
               
               if (dg%pa.gt.0) then
                  
                  if  (ze_sensor2.gt.1.0e-12.and.ze_sensor1.gt.1.0e-12 ) then
                     
                     ze_sensor(l) = log10( ( (ze_sensor1 )**roo) / (ze_sensor2**roo) )
                     dg%slimit1 = log10(kon*real(pdg_el(l))**(- lebesgue**2)) - kons
                     
                  else
                     
                     ze_sensor(l) = 1.d0
                     dg%slimit1 = 2.d0
                     
                  endif

                  if (qx_sensor2.gt.1.0e-12.and.qx_sensor1.gt.1.0e-12 ) then
                     
                     qx_sensor(l) = log10( ( (qx_sensor1)**roo) / (qx_sensor2**roo) )
                     dg%slimit2 = log10(kon*real(pdg_el(l))**(- lebesgue**2)) - kons
                     
                  else
                     
                     qx_sensor(l) = 1.d0
                     dg%slimit2 = 2.d0
                     
                  endif
                  
                  if ( qy_sensor2.gt.1.0e-12.and.qy_sensor1.gt.1.0e-12 ) then
                     
                     qy_sensor(l) = log10( ((qy_sensor1)**roo) / (qy_sensor2**roo) )
                     dg%slimit3 = log10(kon*real(pdg_el(l))**(- lebesgue**2)) - kons
                     
                  else
                     
                     qy_sensor(l) = 1.d0
                     dg%slimit3 = 2.d0
                     
                  endif
                  
#ifdef TRACE
                  if ( iota_sensor2.gt.1.0e-12.and.iota_sensor1.gt.1.0e-12 ) then
                     
                     iota_sensor(l) = log10( ( (iota_sensor1)**roo) / (iota_sensor2**roo) )
                     dg%slimit4 = log10(kon*real(pdg_el(l))**(- lebesgue**2)) - kons
                     
                  else
                     
                     iota_sensor(l) = 1.d0
                     dg%slimit4 = 2.d0
                     
                     
                  endif
#endif
                  
#ifdef CHEM            
                  if ( iota_sensor2.gt.1.0e-12.and.iota_sensor1.gt.1.0e-12) then
                     
                     iota_sensor(l) = log10( ( (iota_sensor1 )**roo) / (iota_sensor2**roo) )
                     dg%slimit4 = log10(kon*real(pdg_el(l))**(- lebesgue**2)) - kons
                     
                  else
                     
                     iota_sensor(l) = 1.d0
                     dg%slimit4 = 2.d0
                     
                  endif
                  
                  if ( iota2_sensor2.gt.1.0e-12.and.iota2_sensor1.gt.1.0e-12 ) then
                     
                     iota2_sensor(l) = log10( ( (iota2_sensor1)**roo) / (iota2_sensor2**roo) )
                     dg%slimit5 = log10(kon*real(pdg_el(l))**(- lebesgue**2)) - kons
                     
                  else
                     
                     iota2_sensor(l) = 1.d0
                     dg%slimit5 = 2.d0
                     
                  endif
#endif

#ifdef SED_LAY
                  if ( tbed_sensor2.gt.1.0e-12.and.tbed_sensor1.gt.1.0e-12 ) then
                     
                     tbed_sensor(l) = log10( ( (tbed_sensor1)**roo) / (tbed_sensor2**roo) )
                     slimit6 = log10(kon*real(pdg_el(l))**(- lebesgue**2)) - kons
                     
                  else
                     
                     tbed_sensor(l) = 1.d0
                     slimit6 = 2.d0
                     
                     
                  endif
#endif
                  
               else 
                  
                  ze_sensor(l) = dg%slimit1
                  qx_sensor(l) = dg%slimit2
                  qy_sensor(l) = dg%slimit3
                  
#ifdef TRACE
                  iota_sensor(l) = dg%slimit4
#endif
                  
#ifdef CHEM            
                  iota_sensor(l) = dg%slimit4
                  iota2_sensor(l) = dg%slimit5
#endif

#ifdef SED_LAY
                  tbed_sensor(l) = slimit6
#endif
                  
               endif

#ifdef TRACE

#elif CHEM

#elif SED_LAY

#else  

!.....if the sensor is less than the limit and the p is low increase
!.....the order of p

               
               if ( ((ze_sensor(l).lt.dg%slimit1).and.&
              (qx_sensor(l).lt.dg%slimit2).and.&
              (qy_sensor(l).lt.dg%slimit3)).and.&
              (pdg_el(l).lt.dg%ph).and.&
              (it-dg%pcount(l)).ge.dg%plimit )  then
                  
                  pdg_el(l) = pdg_el(l) + 1
                  dg%dofs(l) = (pdg_el(l) + 1 )*(pdg_el(l) +2 ) / 2
                  
                  dg%pcount(l) = it
                  
!.....if the sensor is more than the limit and the p is high decrease
!.....the order of p
                  
               elseif ( ((ze_sensor(l).ge.dg%slimit1).or.&
                 (qx_sensor(l).ge.dg%slimit2).or.&
                 (qy_sensor(l).ge.dg%slimit3)).and.&
                 (pdg_el(l).gt.dg%pl) ) then !.and.
                                !&                    (it-dg%pcount(l)).ge.dg%plimit ) then
                  
                  pdg_el(l) = pdg_el(l) - 1
                  dg%dofs(l) = (pdg_el(l) + 1 )*(pdg_el(l) + 2 ) / 2
                  
                  dg%ze(dg%dofs(l)+1:dg%dofh,l,irkp+1) = 0.d0
                  dg%qx(dg%dofs(l)+1:dg%dofh,l,irkp+1) = 0.d0
                  dg%qy(dg%dofs(l)+1:dg%dofh,l,irkp+1) = 0.d0  

                                !dg%pcount(l) = it 

               endif
#endif
               
#ifdef TRACE  

!.....if the sensor is less than the limit and the p is low increase
!.....the order of p
               
               if ( ( (ze_sensor(l).lt.dg%slimit1).and.&
              (qx_sensor(l).lt.dg%slimit2).and.&
              (qy_sensor(l).lt.dg%slimit3).and.&
              (iota_sensor(l).lt.dg%slimit4)).and.&
              (pdg_el(l).lt.dg%ph).and.&
              (it-dg%pcount(l)).ge.dg%plimit ) then
                  
                  pdg_el(l) = pdg_el(l) + 1
                  dg%dofs(l) = (pdg_el(l) + 1 )*(pdg_el(l) +2 ) / 2
                  
                  dg%pcount(l) = it
                  
!.....if the sensor is more than the limit and the p is high decrease
!.....the order of p
                  
               elseif ( ((ze_sensor(l).ge.dg%slimit1).or.&
                 (qx_sensor(l).ge.dg%slimit2).or.&
                 (qy_sensor(l).ge.dg%slimit3).or.&
                 (iota_sensor(l).ge.dg%slimit4)).and.&
                 (pdg_el(l).gt.dg%pl) ) then !.and.
                                !&                    (it-dg%pcount(l)).ge.dg%plimit ) then
                  
                  pdg_el(l) = pdg_el(l) - 1
                  dg%dofs(l) = (pdg_el(l) + 1 )*(pdg_el(l) + 2 ) / 2
                  
                  dg%ze(dg%dofs(l)+1:dg%dofh,l,irkp+1) = 0.d0
                  dg%iota(dg%dofs(l)+1:dg%dofh,l,irkp+1) = 0.d0
                  dg%qx(dg%dofs(l)+1:dg%dofh,l,irkp+1) = 0.d0
                  dg%qy(dg%dofs(l)+1:dg%dofh,l,irkp+1) = 0.d0
                  
                                !dg%pcount(l) = it
                  
               endif
#endif

#ifdef CHEM 
               
!.....if the sensor is less than the limit and the p is low increase
!.....the order of p
               
               
               if ( ( (ze_sensor(l).lt.dg%slimit1).and.&
              (qx_sensor(l).lt.dg%slimit2).and.&
              (qy_sensor(l).lt.dg%slimit3).and.&
              (iota_sensor(l).lt.dg%slimit4).and.&
              (iota2_sensor(l).lt.dg%slimit5)).and. &
              (pdg_el(l).lt.dg%ph).and. &
              (it-dg%pcount(l)).ge.dg%plimit) then
                  
                  pdg_el(l) = pdg_el(l) + 1
                  dg%dofs(l) = (pdg_el(l) + 1 )*(pdg_el(l) +2 ) / 2
                  
                  dg%pcount(l) = it
                  
!.....if the sensor is more than the limit and the p is high decrease
!.....the order of p
                  
               elseif ( ((ze_sensor(l).ge.dg%slimit1).or.&
                 (qx_sensor(l).ge.dg%slimit2).or.&
                 (qy_sensor(l).ge.dg%slimit3).or.&
                 (iota_sensor(l).ge.dg%slimit4).or.&
                 (iota2_sensor(l).ge.dg%slimit5)).and.&
                 (pdg_el(l).gt.dg%pl) ) then !.and.
                                !&                    (it-dg%pcount(l)).ge.dg%plimit ) then
                  
                  pdg_el(l) = pdg_el(l) - 1
                  dg%dofs(l) = (pdg_el(l) + 1 )*(pdg_el(l) + 2 ) / 2
                  
                  dg%ze(dg%dofs(l)+1:dg%dofh,l,irkp+1) = 0.d0
                  dg%iota(dg%dofs(l)+1:dg%dofh,l,irkp+1) = 0.d0
                  dg%iota2(dg%dofs(l)+1:dg%dofh,l,irkp+1) = 0.d0
                  dg%qx(dg%dofs(l)+1:dg%dofh,l,irkp+1) = 0.d0
                  dg%qy(dg%dofs(l)+1:dg%dofh,l,irkp+1) = 0.d0

                                !dg%pcount(l) = it
                  
               endif
#endif

#ifdef SED_LAY 

!.....if the sensor is less than the limit and the p is low increase
!.....the order of p
               
               if ( ( (ze_sensor(l).lt.dg%slimit1).and.&
              (qx_sensor(l).lt.dg%slimit2).and.&
              (qy_sensor(l).lt.dg%slimit3).and.&
              (tbed_sensor(l).lt.slimit6)).and.&
              (pdg_el(l).lt.dg%ph).and.&
              (it-dg%pcount(l)).ge.dg%plimit ) then
                  
                  pdg_el(l) = pdg_el(l) + 1
                  dg%dofs(l) = (pdg_el(l) + 1 )*(pdg_el(l) +2 ) / 2
                  
                  dg%pcount(l) = it
                  
!.....if the sensor is more than the limit and the p is high decrease
!.....the order of p
                  
               elseif ( ((ze_sensor(l).ge.dg%slimit1).or.&
                 (qx_sensor(l).ge.dg%slimit2).or.&
                 (qy_sensor(l).ge.dg%slimit3).or.&
                 (tbed_sensor(l).ge.slimit6)).and.&
                 (pdg_el(l).gt.dg%pl) ) then !.and.
                                !&                    (it-dg%pcount(l)).ge.dg%plimit ) then
                  
                  pdg_el(l) = pdg_el(l) - 1
                  dg%dofs(l) = (pdg_el(l) + 1 )*(pdg_el(l) + 2 ) / 2
                  
                  dg%ze(dg%dofs(l)+1:dg%dofh,l,irkp+1) = 0.d0
                  dg%bed(dg%dofs(l)+1:dg%dofh,l,irkp+1,:) = 0.d0
                  dg%qx(dg%dofs(l)+1:dg%dofh,l,irkp+1) = 0.d0
                  dg%qy(dg%dofs(l)+1:dg%dofh,l,irkp+1) = 0.d0
                  
                                !dg%pcount(l) = it
                  
               endif
#endif
               
            endif
!................................................................

            if (dg%gflag.eq.1) then ! dioristic algorithm ON

               if (ze_sensor2.gt.1.0e-12.and.ze_sensor1.gt.1.0e-12 ) then
                  
                  ze_sensor(l) = log10(( ze_sensor1**roo) / (ze_sensor2**roo))
                  
               else

                  ze_sensor(l) = 0.d0

               endif 
               
               if  (qx_sensor2.gt.1.0e-12.and.qx_sensor1.gt.1.0e-12 ) then
                  
                  qx_sensor(l) = log10(( qx_sensor1**roo) / (qx_sensor2**roo))
                  
               else
                  
                  qx_sensor(l) = 0.d0

               endif

               if  (qy_sensor2.gt.1.0e-12.and.qy_sensor1.gt.1.0e-12 ) then
                  
                  qy_sensor(l) = log10(( qy_sensor1**roo) / (qy_sensor2**roo))

               else

                  qy_sensor(l) = 0.d0

               endif

               
#ifdef TRACE
               if  (iota_sensor2.gt.1.0e-12.and.iota_sensor1.gt.1.0e-12 ) then
                  
                  iota_sensor(l) = log10(( iota_sensor1**roo) / (iota_sensor2**roo))
                  
               else
                  
                  iota_sensor(l) = 0.d0
                  
               endif
#endif
               
#ifdef CHEM         
               if  (iota_sensor2.gt.1.0e-12.and.iota_sensor1.gt.1.0e-12 ) then
                  
                  iota_sensor(l) = log10(( iota_sensor1**roo) / (iota_sensor2**roo))
                  
               else

                  iota_sensor(l) = 0.d0

               endif

               if  (iota2_sensor2.gt.1.0e-12.and.iota2_sensor1.gt.1.0e-12 ) then
                  
                  iota2_sensor(l) = log10(( iota2_sensor1**roo) / (iota2_sensor2**roo)) 
                  
               else

                  iota_sensor(l) = 0.d0

               endif
#endif

             
#ifdef SED_LAY
               if  (tbed_sensor2.gt.1.0e-12.and.tbed_sensor1.gt.1.0e-12 ) then
                  
                  tbed_sensor(l) = log10(( tbed_sensor1**roo) / (tbed_sensor2**roo))
                  
               else
                  
                  tbed_sensor(l) = 0.d0
                  
               endif
#endif
               
            endif

         enddo


         if (dg%gflag.eq.1) then   !if dioristic algorithm ON

            avg_zesen = 0.d0
            avg_qxsen = 0.d0
            avg_qysen = 0.d0

#ifdef TRACE
            avg_iotasen = 0.d0
#endif

#ifdef CHEM
            avg_iotasen = 0.d0
            avg_iota2sen = 0.d0
#endif

#ifdef SED_LAY
            avg_tbedsen = 0.d0
#endif
            
            do l=1,ne
               
#ifdef CMPI

               if (reselem(l)) then
                  
#endif

                  if (ze_sensor(l).ne.0.d0) then

                     minze_sensor = min(ze_sensor(l),minze_sensor)
                     maxze_sensor = max(ze_sensor(l),maxze_sensor)
                     
                  endif

                  if (qx_sensor(l).ne.0.d0) then

                     minqx_sensor = min(qx_sensor(l),minqx_sensor)
                     maxqx_sensor = max(qx_sensor(l),maxqx_sensor)

                  endif

                  if (qy_sensor(l).ne.0.d0) then

                     minqy_sensor = min(qy_sensor(l),minqy_sensor)
                     maxqy_sensor = max(qy_sensor(l),maxqy_sensor)
                     
                  endif
#ifdef TRACE
                  if (iota_sensor(l).ne.0.d0) then

                     miniota_sensor = min(iota_sensor(l),miniota_sensor)
                     maxiota_sensor = max(iota_sensor(l),maxiota_sensor)

                  endif
#endif

#ifdef CHEM
                  if (iota_sensor(l).ne.0.d0) then

                     miniota_sensor = min(iota_sensor(l),miniota_sensor)
                     maxiota_sensor = max(iota_sensor(l),maxiota_sensor)

                  endif

                  if (iota2_sensor(l).ne.0.d0) then
                     
                     miniota2_sensor = min(iota2_sensor(l),miniota2_sensor)
                     maxiota2_sensor = max(iota2_sensor(l),maxiota2_sensor)
                     
                  endif
#endif

#ifdef SED_LAY
                  if (tbed_sensor(l).ne.0.d0) then

                     mintbed_sensor = min(tbed_sensor(l),mintbed_sensor)
                     maxtbed_sensor = max(tbed_sensor(l),maxtbed_sensor)

                  endif
#endif
                  
                  avg_zesen = sum(ze_sensor)/dg%MNES
                  avg_qxsen = sum(qx_sensor)/dg%MNES
                  avg_qysen = sum(qy_sensor)/dg%MNES
                  
#ifdef TRACE
                  avg_iotasen = sum(iota_sensor)/dg%MNES
#endif
                  
#ifdef CHEM             
                  avg_iotasen = sum(iota_sensor)/dg%MNES
                  avg_iota2sen = sum(iota2_sensor)/dg%MNES
#endif

#ifdef SED_LAY
                  avg_tbedsen = sum(tbed_sensor)/dg%MNES
#endif

#ifdef CMPI

               endif

#endif
               
            enddo

#ifdef CMPI
            CALL PARA_MAX(MAXZE_SENSOR)
            CALL PARA_MIN(MINZE_SENSOR)
            CALL PARA_SUM(AVG_ZESEN)
            CALL PARA_SUM(AVG_QXSEN)
            CALL PARA_MAX(MAXQX_SENSOR)
            CALL PARA_MIN(MINQX_SENSOR)
            CALL PARA_SUM(AVG_QYSEN)
            CALL PARA_MAX(MAXQY_SENSOR)
            CALL PARA_MIN(MINQY_SENSOR)
#endif
#ifdef CMPI_TRACE
            CALL PARA_MAX(MAXiota_SENSOR)
            CALL PARA_MIN(MINiota_SENSOR)
            CALL PARA_SUM(AVG_iotaSEN)
#endif

#ifdef CMPI_CHEM
            CALL PARA_MAX(MAXiota_SENSOR)
            CALL PARA_MIN(MINiota_SENSOR)
            CALL PARA_SUM(AVG_iotaSEN)
            CALL PARA_MAX(MAXiota2_SENSOR)
            CALL PARA_MIN(MINiota2_SENSOR)
            CALL PARA_SUM(AVG_iota2SEN)
#endif


#ifdef CMPI_SED_LAY
            CALL PARA_MAX(MAXtbed_SENSOR)
            CALL PARA_MIN(MINtbed_SENSOR)
            CALL PARA_SUM(AVG_tbedSEN)
#endif

            ze_delta = 0.d0
            qx_delta = 0.d0
            qy_delta = 0.d0

#ifdef TRACE
            iota_delta = 0.d0
#endif

#ifdef CHEM
            iota_delta = 0.d0
            iota2_delta = 0.d0
#endif

#ifdef SED_LAY
            tbed_delta = 0.d0
#endif
            
            ze_delta = ( dg%diorism / 100.d0 )* (abs(maxze_sensor - minze_sensor) )
            qx_delta = ( dg%diorism / 100.d0 )* (abs(maxqx_sensor - minqx_sensor) )
            qy_delta = ( dg%diorism / 100.d0 )* (abs(maxqy_sensor - minqy_sensor) )   
            
#ifdef TRACE
            iota_delta = ( dg%diorism / 100.d0 )* (abs(maxiota_sensor - miniota_sensor) )
#endif
            
#ifdef CHEM
            iota_delta = ( dg%diorism / 100.d0 )* (abs(maxiota_sensor - miniota_sensor) )
            iota2_delta = ( dg%diorism / 100.d0 )* (abs(maxiota2_sensor - miniota2_sensor) )
#endif

#ifdef SED_LAY
            tbed_delta = ( dg%diorism / 100.d0 )* (abs(maxtbed_sensor - mintbed_sensor) )
#endif
            
            do l = 1,ne
               
               temp_ze = 0.d0  
               temp_qx = 0.d0
               temp_qy = 0.d0
               
#ifdef TRACE
               temp_iota = 0.d0
#endif

#ifdef CHEM
               temp_iota = 0.d0
               temp_iota2 = 0.d0
#endif

#ifdef SED_LAY
               temp_tbed = 0.d0
#endif

               temp_ze = abs(ze_sensor(l) - avg_zesen)
               temp_qx = abs(qx_sensor(l) - avg_qxsen)
               temp_qy = abs(qy_sensor(l) - avg_qysen)

#ifdef TRACE
               temp_iota = abs(iota_sensor(l) - avg_iotasen)
#endif

#ifdef CHEM
               temp_iota = abs(iota_sensor(l) - avg_iotasen)
               temp_iota2 = abs(iota2_sensor(l) - avg_iota2sen)
#endif

#ifdef SED_LAY
               temp_tbed = abs(tbed_sensor(l) - avg_tbedsen)
#endif

               
#ifdef TRACE

#elif CHEM

#elif SED_LAY

#else


               !print*,'temp=', temp_ze,temp_qx,temp_qy
               !print*,'delta=', ze_delta,qx_delta,qy_delta
               
!.....if the sensor is less than delta and the p is low increase
!.....the order of p
               
               if ( ( (temp_qy.lt.qy_delta).and.&
              (temp_qx.lt.qx_delta).and.&
              (temp_ze.lt.ze_delta) ).and.&
              (pdg_el(l).lt.dg%ph).and.&
              (it-dg%pcount(l)).ge.dg%plimit ) then
                  
                  pdg_el(l) = pdg_el(l) + 1 
                  dg%dofs(l) = (pdg_el(l) +1 )*(pdg_el(l) +2 ) / 2
                  
                  dg%pcount(l) = it
                  
!.....if the sensor is more than delta and the p is high decrease
!.....the order of p
                  
               elseif ( ((temp_ze.ge.ze_delta).or.&
                 (temp_qx.ge.qx_delta).or.&
                 (temp_qy.ge.qy_delta) ).and.&
                 (pdg_el(l).gt.dg%pl) ) then !.and.
                                !&                    (it-dg%pcount(l)).ge.dg%plimit ) then
                  
                  pdg_el(l) = pdg_el(l) - 1
                  dg%dofs(l) = (pdg_el(l) + 1 )*(pdg_el(l) + 2 ) / 2
                  
                  dg%ze(dg%dofs(l)+1:dg%dofh,l,irkp+1) = 0.d0
                  dg%qx(dg%dofs(l)+1:dg%dofh,l,irkp+1) = 0.d0
                  dg%qy(dg%dofs(l)+1:dg%dofh,l,irkp+1) = 0.d0

                                !dg%pcount(l) = it
                  
               endif
#endif

#ifdef TRACE

!.....if the sensor is less than delta and the p is low increase
!.....the order of p

               if ( ((temp_qy.lt.qy_delta).and.&
              (temp_qx.lt.qx_delta).and.&
              (temp_ze.lt.ze_delta).and.&
              (temp_iota.lt.iota_delta)).and.&
              (pdg_el(l).lt.dg%ph).and.&
              (it-dg%pcount(l)).ge.dg%plimit ) then

                  pdg_el(l) = pdg_el(l) + 1 
                  dg%dofs(l) = (pdg_el(l) +1 )*(pdg_el(l) +2 ) / 2

                  dg%pcount(l) = it
                  
!.....if the sensor is less than delta and the p is high decrease
!.....the order of p

               elseif (((temp_ze.ge.ze_delta).or.&
                 (temp_qx.ge.qx_delta).or.&
                 (temp_qy.ge.qy_delta).or.&
                 (temp_iota.ge.iota_delta)).and.&
                 (pdg_el(l).gt.dg%pl) ) then !.and.
                                !&                    (it-dg%pcount(l)).ge.dg%plimit) then

                  pdg_el(l) = pdg_el(l) - 1
                  dg%dofs(l) = (pdg_el(l) + 1 )*(pdg_el(l) + 2 ) / 2
                  
                  dg%ze(dg%dofs(l)+1:dg%dofh,l,irkp+1) = 0.d0
                  dg%qx(dg%dofs(l)+1:dg%dofh,l,irkp+1) = 0.d0
                  dg%qy(dg%dofs(l)+1:dg%dofh,l,irkp+1) = 0.d0
                  dg%iota(dg%dofs(l)+1:dg%dofh,l,irkp+1) = 0.d0

                                !dg%pcount(l) = it
                  
               endif
#endif

#ifdef CHEM           

!.....if the sensor is less than delta and the p is low increase
!.....the order of p

               if ( ((temp_qy.lt.qy_delta).and.&
              (temp_qx.lt.qx_delta).and.&
              (temp_ze.lt.ze_delta).and.&
              (temp_iota.lt.iota_delta).and.&
              (temp_iota2.lt.iota2_delta)).and.&
              (pdg_el(l).lt.dg%ph).and.&
              (it-dg%pcount(l)).ge.dg%plimit ) then

                  pdg_el(l) = pdg_el(l) + 1 
                  dg%dofs(l) = (pdg_el(l) +1 )*(pdg_el(l) +2 ) / 2

                  dg%pcount(l) = it
                  
!.....if the sensor is less than delta and the p is high decrease
!.....the order of p

               elseif (((temp_ze.ge.ze_delta).or.&
                 (temp_qx.ge.qx_delta).or.&
                 (temp_qy.ge.qy_delta).or.&
                 (temp_iota.ge.iota_delta).or.&
                 (temp_iota2.ge.iota2_delta)).and.&
                 (pdg_el(l).gt.dg%pl) ) then !.and.
                                !&                    (it-dg%pcount(l)).ge.dg%plimit) then

                  pdg_el(l) = pdg_el(l) - 1
                  dg%dofs(l) = (pdg_el(l) + 1 )*(pdg_el(l) + 2 ) / 2
                  
                  dg%ze(dg%dofs(l)+1:dg%dofh,l,irkp+1) = 0.d0
                  dg%qx(dg%dofs(l)+1:dg%dofh,l,irkp+1) = 0.d0
                  dg%qy(dg%dofs(l)+1:dg%dofh,l,irkp+1) = 0.d0
                  dg%iota(dg%dofs(l)+1:dg%dofh,l,irkp+1) = 0.d0
                  dg%iota2(dg%dofs(l)+1:dg%dofh,l,irkp+1) = 0.d0

                                !dg%pcount(l) = it
                  
               endif         
#endif

#ifdef SED_LAY

!.....if the sensor is less than delta and the p is low increase
!.....the order of p

               if ( ((temp_qy.lt.qy_delta).and.&
              (temp_qx.lt.qx_delta).and.&
              (temp_ze.lt.ze_delta).and.&
              (temp_tbed.lt.tbed_delta)).and.&
              (pdg_el(l).lt.dg%ph).and.&
              (it-dg%pcount(l)).ge.dg%plimit ) then

                  pdg_el(l) = pdg_el(l) + 1 
                  dg%dofs(l) = (pdg_el(l) +1 )*(pdg_el(l) +2 ) / 2

                  dg%pcount(l) = it
                  
!.....if the sensor is less than delta and the p is high decrease
!.....the order of p

               elseif (((temp_ze.ge.ze_delta).or.&
                 (temp_qx.ge.qx_delta).or.&
                 (temp_qy.ge.qy_delta).or.&
                 (temp_tbed.ge.tbed_delta)).and.&
                 (pdg_el(l).gt.dg%pl) ) then !.and.
                                !&                    (it-dg%pcount(l)).ge.dg%plimit) then

                  pdg_el(l) = pdg_el(l) - 1
                  dg%dofs(l) = (pdg_el(l) + 1 )*(pdg_el(l) + 2 ) / 2
                  
                  dg%ze(dg%dofs(l)+1:dg%dofh,l,irkp+1) = 0.d0
                  dg%qx(dg%dofs(l)+1:dg%dofh,l,irkp+1) = 0.d0
                  dg%qy(dg%dofs(l)+1:dg%dofh,l,irkp+1) = 0.d0
                  dg%bed(dg%dofs(l)+1:dg%dofh,l,irkp+1,:) = 0.d0

                                !dg%pcount(l) = it
                  
               endif
#endif


            enddo
            
         endif

      endif

!.....deal with the interior barrier problem, by making the front and back elements agree
!.....at the min of the union for stability

      if (dg%nibseg.gt.0) then

         do j = 1,dg%nibseg

            if ( pdg_el(dg%nedel(1,dg%nibsegn(1,j))).ne.pdg_el(dg%nedel(1,dg%nibsegn(2,j))) ) then

               pdg_el(dg%nedel(1,dg%nibsegn(1,j))) = max(pdg_el(dg%nedel(1,dg%nibsegn(1,j))), pdg_el(dg%nedel(1,dg%nibsegn(2,j))) )
               pdg_el(dg%nedel(2,dg%nibsegn(1,j))) = max(pdg_el(dg%nedel(1,dg%nibsegn(1,j))), pdg_el(dg%nedel(1,dg%nibsegn(2,j))) )

               dg%dofs(dg%nedel(1,dg%nibsegn(1,j))) = max(dg%dofs(dg%nedel(1,dg%nibsegn(1,j))), dg%dofs(dg%nedel(1,dg%nibsegn(2,j))) )
               dg%dofs(dg%nedel(2,dg%nibsegn(1,j))) = max(dg%dofs(dg%nedel(1,dg%nibsegn(1,j))), dg%dofs(dg%nedel(1,dg%nibsegn(2,j))) )

            endif

         enddo

      endif
      

      end subroutine p_enrichment


  

