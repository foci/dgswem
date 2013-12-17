c***********************************************************************
c     
c     subroutine p_enrichment()
c     
c     In concert, was born and evolved -- cem
c     
c***********************************************************************
      subroutine p_enrichment(it,irkp)
      
c.....use appropriate modules

      use global
      use dg

#ifdef CMPI
      USE MESSENGER_ELEM
#endif

      implicit none

c.....declare local variables

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
      real(sz) ze_sensor(mne) 
      real(sz) qx_sensor(mne)
      real(sz) qy_sensor(mne)
#ifdef TRACE
      real(sz) iota_sensor1, iota_sensor2, iota_center, iota_mid
      real(sz) maxiota_sensor, miniota_sensor, iota_delta
      real(sz) avg_iotasen, temp_iota, iota_sensor(mne)
#endif
#ifdef CHEM
      real(sz) iota_sensor1, iota_sensor2, iota_center, iota_mid
      real(sz) maxiota_sensor, miniota_sensor, iota_delta
      real(sz) avg_iotasen, temp_iota, iota_sensor(mne)
      real(sz) iota2_sensor1, iota2_sensor2, iota2_center 
      real(sz) iota2_mid, slimit5, maxiota2_sensor, miniota2_sensor
      real(sz) iota2_delta, avg_iota2sen, temp_iota2, iota2_sensor(mne)
      
#endif
#ifdef SED_LAY
      real(sz) tbed_sensor1, tbed_sensor2
      real(sz) tbed_center, tbed_mid
      real(sz) mintbed_sensor,maxtbed_sensor
      real(sz) avg_tbedsen, temp_tbed
      real(sz) tbed_sensor(mne)
      real(sz) slimit6, tbed_delta
#endif


c......Set the initial arrays, note we work over the total bed load
c......not the individual layers	

      if (pflag.eq.1) then
         
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
            
            pa = pdg_el(l)
            
c..........................
c.... type 1 p-enrichment |
c..........................
            
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
            
            if (pa.ge.1) then   ! since the pa=0 case is vacuous
               
               do k = 1,dofs(l)
                  
                  ze_center = ze_center + ze(k,l,irk+1)*phi_center(k,pa)
                  qx_center = qx_center + qx(k,l,irk+1)*phi_center(k,pa)
                  qy_center = qy_center + qy(k,l,irk+1)*phi_center(k,pa)
                  
#ifdef TRACE
                  iota_center = iota_center + iota(k,l,irk+1)*phi_center(k,pa)
#endif
                  
#ifdef CHEM            
                  iota_center = iota_center + iota(k,l,irk+1)*phi_center(k,pa)
                  iota2_center = iota2_center + iota2(k,l,irk+1)*phi_center(k,pa)
#endif

#ifdef SED_LAY
                  do ll=1,layers
                     tbed_center = tbed_center + bed(k,l,irk+1,ll)*phi_center(k,pa)
                  enddo
#endif
                  
               enddo
               
            endif
            
c......Compute the sensor using edges in each variable
            
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
               
               dist = sqrt( (x_mid(i) - x_center)**2.d0 +
     &              (y_mid(i) - y_center)**2.d0 )
               
               
               if (pa.ge.1) then
                  
                  do k = 1,dofs(l)
                     
                     ze_mid = ze_mid + ze(k,l,irk+1)*phi_mid(k,i,pa)
                     qx_mid = qx_mid + qx(k,l,irk+1)*phi_mid(k,i,pa)
                     qy_mid = qy_mid + qy(k,l,irk+1)*phi_mid(k,i,pa)
                     
#ifdef TRACE
                     iota_mid = iota_mid + iota(k,l,irk+1)*phi_mid(k,i,pa)
#endif
                     
#ifdef CHEM               
                     iota_mid = iota_mid + iota(k,l,irk+1)*phi_mid(k,i,pa)
                     iota2_mid = iota2_mid + iota2(k,l,irk+1)*phi_mid(k,i,pa)
#endif

#ifdef SED_LAY
                     do ll=1,layers
                        tbed_mid = tbed_mid + bed(k,l,irk+1,ll)*phi_mid(k,i,pa)
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
            
            
            if (gflag.eq.0) then ! dioristic algorithm is OFF
               
#ifdef TRACE

#elif CHEM

#elif SED_LAY

#else
               
c.....if the sensor is greater than the limit and the p is low increase
c.....the order of p

               if ( ( (qy_sensor(l).ge.slimit).or.
     &              (qx_sensor(l).ge.slimit).or.
     &              (ze_sensor(l).ge.slimit) ).and.
     &              (pdg_el(l).lt.ph) ) then
                  
                  pdg_el(l) = pdg_el(l) + 1 
                  dofs(l) = (pdg_el(l) +1 )*(pdg_el(l) +2 ) / 2
                  
                  pcount(l) = it
                  
c.....if the sensor is less than the limit and the p is high decrease
c.....the order of p
                  
               elseif ( ((ze_sensor(l).lt.slimit).and.
     &                 (qx_sensor(l).lt.slimit).and.
     &                 (qy_sensor(l).lt.slimit) ).and.
     &                 (pdg_el(l).gt.pl).and.
     &                 ( (it-pcount(l)).gt.plimit ) ) then
                  
                  pdg_el(l) = pdg_el(l) - 1
                  dofs(l) = (pdg_el(l) + 1 )*(pdg_el(l) + 2 ) / 2
                  
                  ze(dofs(l)+1:dofh,l,irk+1) = 0.d0
                  qx(dofs(l)+1:dofh,l,irk+1) = 0.d0
                  qy(dofs(l)+1:dofh,l,irk+1) = 0.d0
                  
                  pcount(l) = it
                  
                  
               elseif (pdg_el(l).eq.0.and.( (it-pcount(l)).gt.plimit ) ) then
                  
                  pdg_el(l) = pdg_el(l) + 1 
                  dofs(l) = (pdg_el(l) + 1 )*(pdg_el(l) + 2 ) / 2
                  
                  pcount(l) = it
                  
               endif

#endif

#ifdef TRACE       
               
c.....if the sensor is greater than the limit and the p is low increase
c.....the order of p

               if ( ( (qy_sensor(l).ge.slimit).or.
     &              (qx_sensor(l).ge.slimit).or.
     &              (ze_sensor(l).ge.slimit).or.
     &              ( iota_sensor(l).ge.slimit) ).and.
     &              (pdg_el(l).lt.ph).and.
     &              ( (it-pcount(l)).gt.plimit ) ) then
                  
                  pdg_el(l) = pdg_el(l) + 1 
                  dofs(l) = (pdg_el(l) +1 )*(pdg_el(l) +2 ) / 2
                  
                  pcount(l) = it
                  
                  
c.....if the sensor is less than the limit and the p is high decrease
c.....the order of p
                  
                  
               elseif ( ((ze_sensor(l).lt.slimit).and.
     &                 (qx_sensor(l).lt.slimit).and.
     &                 (qy_sensor(l).lt.slimit).and.
     &                 (iota_sensor(l).lt.slimit) ).and.
     &                 (pdg_el(l).gt.pl).and.
     &                 ( (it-pcount(l)).gt.plimit ) ) then
                  
                  pdg_el(l) = pdg_el(l) - 1
                  dofs(l) = (pdg_el(l) + 1 )*(pdg_el(l) + 2 ) / 2
                  
                  ze(dofs(l)+1:dofh,l,irk+1) = 0.d0
                  iota(dofs(l)+1:dofh,l,irk+1) = 0.d0
                  qx(dofs(l)+1:dofh,l,irk+1) = 0.d0
                  qy(dofs(l)+1:dofh,l,irk+1) = 0.d0
                  
                  pcount(l) = it
                  
                  
               elseif (pdg_el(l).eq.0.and.( (it-pcount(l)).gt.plimit ) ) then
                  
                  pdg_el(l) = pdg_el(l) + 1 
                  dofs(l) = (pdg_el(l) + 1 )*(pdg_el(l) + 2 ) / 2
                  
                  pcount(l) = it
                  
               endif

#endif
               
#ifdef CHEM       
               
c.....if the sensor is greater than the limit and the p is low increase
c.....the order of p
               
               if ( ( (qy_sensor(l).ge.slimit).or.
     &              (qx_sensor(l).ge.slimit).or.
     &              (ze_sensor(l).ge.slimit).or.
     &              (iota_sensor(l).ge.slimit).or.
     &              (iota2_sensor(l).ge.slimit) ).and.
     &              (pdg_el(l).lt.ph) ) then
                  
                  pdg_el(l) = pdg_el(l) + 1 
                  dofs(l) = (pdg_el(l) +1 )*(pdg_el(l) +2 ) / 2
                  
                  pcount(l) = it
                  
                  
c.....if the sensor is less than the limit and the p is high decrease
c.....the order of p
                  
               elseif ( ((ze_sensor(l).lt.slimit).and.
     &                 (qx_sensor(l).lt.slimit).and.
     &                 (qy_sensor(l).lt.slimit).and.
     &                 (iota_sensor(l).lt.slimit).and.
     &                 (iota2_sensor(l).lt.slimit)).and.
     &                 (pdg_el(l).gt.pl).and.
     &                 ( (it-pcount(l)).gt.plimit ) ) then
                  
                  pdg_el(l) = pdg_el(l) - 1
                  dofs(l) = (pdg_el(l) + 1 )*(pdg_el(l) + 2 ) / 2
                  
                  ze(dofs(l)+1:dofh,l,irk+1) = 0.d0
                  iota(dofs(l)+1:dofh,l,irk+1) = 0.d0
                  iota2(dofs(l)+1:dofh,l,irk+1) = 0.d0
                  qx(dofs(l)+1:dofh,l,irk+1) = 0.d0
                  qy(dofs(l)+1:dofh,l,irk+1) = 0.d0
                  
                  pcount(l) = it
                  
                  
               elseif (pdg_el(l).eq.0.and.( (it-pcount(l)).gt.plimit ) ) then
                  
                  pdg_el(l) = pdg_el(l) + 1 
                  dofs(l) = (pdg_el(l) + 1 )*(pdg_el(l) + 2 ) / 2
                  
                  pcount(l) = it
                  
               endif
#endif

#ifdef SED_LAY !when chem or tracers are on, this coupling could be made tighter      
               
c.....if the sensor is greater than the limit and the p is low increase
c.....the order of p

               if ( ( (qy_sensor(l).ge.slimit).or.
     &              (qx_sensor(l).ge.slimit).or.
     &              (ze_sensor(l).ge.slimit).or.
     &              (tbed_sensor(l).ge.slimit) ).and.
     &              (pdg_el(l).lt.ph).and.
     &              ( (it-pcount(l)).gt.plimit ) ) then
                  
                  pdg_el(l) = pdg_el(l) + 1 
                  dofs(l) = (pdg_el(l) +1 )*(pdg_el(l) +2 ) / 2
                  
                  pcount(l) = it
                  
                  
c.....if the sensor is less than the limit and the p is high decrease
c.....the order of p
                  
                  
               elseif ( ((ze_sensor(l).lt.slimit).and.
     &                 (qx_sensor(l).lt.slimit).and.
     &                 (qy_sensor(l).lt.slimit).and.
     &                 (tbed_sensor(l).lt.slimit) ).and.
     &                 (pdg_el(l).gt.pl).and.
     &                 ( (it-pcount(l)).gt.plimit ) ) then
                  
                  pdg_el(l) = pdg_el(l) - 1
                  dofs(l) = (pdg_el(l) + 1 )*(pdg_el(l) + 2 ) / 2
                  
                  ze(dofs(l)+1:dofh,l,irk+1) = 0.d0
                  iota(dofs(l)+1:dofh,l,irk+1) = 0.d0
                  qx(dofs(l)+1:dofh,l,irk+1) = 0.d0
                  qy(dofs(l)+1:dofh,l,irk+1) = 0.d0
                  bed(dofs(l)+1:dofh,l,irk+1,:)=0.d0
                  
                  pcount(l) = it
                  
                  
               elseif (pdg_el(l).eq.0.and.( (it-pcount(l)).gt.plimit ) ) then
                  
                  pdg_el(l) = pdg_el(l) + 1 
                  dofs(l) = (pdg_el(l) + 1 )*(pdg_el(l) + 2 ) / 2
                  
                  pcount(l) = it
                  
               endif

#endif
               
            endif

         enddo

         if (gflag.eq.1) then   ! dioristic algorithm is ON

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
                  
                  avg_zesen = sum(ze_sensor)/MNES
                  avg_qxsen = sum(qx_sensor)/MNES
                  avg_qysen = sum(qy_sensor)/MNES
                  
#ifdef TRACE
                  avg_iotasen = sum(iota_sensor)/MNES
#endif

#ifdef CHEM                  
                  avg_iotasen = sum(iota_sensor)/MNES
                  avg_iota2sen = sum(iota2_sensor)/MNES
#endif

#ifdef SED_LAY
                  avg_tbedsen = sum(tbed_sensor)/MNES
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
            

            ze_delta = ( diorism / 100.d0 )* (abs(maxze_sensor - minze_sensor) )
            qx_delta = ( diorism / 100.d0 )* (abs(maxqx_sensor - minqx_sensor) )
            qy_delta = ( diorism / 100.d0 )* (abs(maxqy_sensor - minqy_sensor) )   
            
#ifdef TRACE
            iota_delta = ( diorism / 100.d0 )* (abs(maxiota_sensor - miniota_sensor) )
#endif
            
#ifdef CHEM
            iota_delta = ( diorism / 100.d0 )* (abs(maxiota_sensor - miniota_sensor) )
            iota2_delta = ( diorism / 100.d0 )* (abs(maxiota2_sensor - miniota2_sensor) )
#endif

#ifdef SED_LAY    
            tbed_delta = ( diorism / 100.d0 )* (abs(maxtbed_sensor - mintbed_sensor) )
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
               
c.....if the sensor is more than the limit and the p is low increase
c.....the order of p
               
               if ( ( (qy_sensor(l).ge.qy_delta).or.
     &              (qx_sensor(l).ge.qx_delta).or.
     &              (ze_sensor(l).ge.ze_delta) ).and.
     &              (pdg_el(l).lt.ph) ) then
                  
                  pdg_el(l) = pdg_el(l) + 1 
                  dofs(l) = (pdg_el(l) +1 )*(pdg_el(l) +2 ) / 2
                  
                  pcount(l) = it
                  
c.....if the sensor is less than the limit and the p is high decrease
c.....the order of p
                  
               elseif ( ((ze_sensor(l).lt.ze_delta).and.
     &                 (qx_sensor(l).lt.qx_delta).and.
     &                 (qy_sensor(l).lt.qy_delta) ).and.
     &                 (pdg_el(l).gt.pl).and.
     &                 ( (it-pcount(l)).gt.plimit ) ) then
                  
                  pdg_el(l) = pdg_el(l) - 1
                  dofs(l) = (pdg_el(l) + 1 )*(pdg_el(l) + 2 ) / 2
                  
                  ze(dofs(l)+1:dofh,l,irk+1) = 0.d0
                  qx(dofs(l)+1:dofh,l,irk+1) = 0.d0
                  qy(dofs(l)+1:dofh,l,irk+1) = 0.d0
                  
                  pcount(l) = it
                  
               elseif (pdg_el(l).eq.0.and.
     &                 ( (it-pcount(l)).gt.plimit ) ) then
                  
                  pdg_el(l) = pdg_el(l) + 1 
                  dofs(l) = (pdg_el(l) + 1 )*(pdg_el(l) + 2 ) / 2
                  
                  pcount(l) = it
                  
               endif
#endif

#ifdef TRACE

c.....if the sensor is more than the limit and the p is low increase
c.....the order of p
               
               
               if ( ( (temp_qy.ge.qy_delta).or.
     &              (temp_qx.ge.qx_delta).or.
     &              (temp_ze.ge.ze_delta).or.
     &              (temp_iota.ge.iota_delta)).and.
     &              (pdg_el(l).lt.ph).and.
     &              ( (it-pcount(l)).gt.plimit ) ) then
                  
                  pdg_el(l) = pdg_el(l) + 1 
                  dofs(l) = (pdg_el(l) +1 )*(pdg_el(l) +2 ) / 2
                  
                  pcount(l) = it
                  
c.....if the sensor is less than the limit and the p is high decrease
c.....the order of p
                  
               elseif ( ((temp_ze.lt.ze_delta).and.
     &                 (temp_qx.lt.qx_delta).and.
     &                 (temp_qy.lt.qy_delta).and.
     &                 (temp_iota.lt.iota_delta)).and.
     &                 (pdg_el(l).gt.pl).and.
     &                 ( (it-pcount(l)).gt.plimit ) ) then
                  
                  pdg_el(l) = pdg_el(l) - 1
                  dofs(l) = (pdg_el(l) + 1 )*(pdg_el(l) + 2 ) / 2
                  
                  ze(dofs(l)+1:dofh,l,irk+1) = 0.d0
                  qx(dofs(l)+1:dofh,l,irk+1) = 0.d0
                  qy(dofs(l)+1:dofh,l,irk+1) = 0.d0
                  iota(dofs(l)+1:dofh,l,irk+1) = 0.d0
                  
                  pcount(l) = it
                  
               elseif (pdg_el(l).eq.0.and.
     &                 ( (it-pcount(l)).gt.plimit ) ) then
                  
                  pdg_el(l) = pdg_el(l) + 1 
                  dofs(l) = (pdg_el(l) + 1 )*(pdg_el(l) + 2 ) / 2
                  
                  pcount(l) = it
                  
               endif
#endif
               
#ifdef CHEM

c.....if the sensor is more than the limit and the p is low increase
c.....the order of p

               if ( ( (qy_sensor(l).ge.qy_delta).or.
     &              (qx_sensor(l).ge.qx_delta).or.
     &              (ze_sensor(l).ge.ze_delta).or.
     &              iota_sensor(l).ge.iota_delta.or.
     &              iota2_sensor(l).ge.iota2_delta).and.
     &              (pdg_el(l).lt.ph) ) then
                  
                  pdg_el(l) = pdg_el(l) + 1 
                  dofs(l) = (pdg_el(l) +1 )*(pdg_el(l) +2 ) / 2
                  
                  pcount(l) = it
                  
c.....if the sensor is less than the limit and the p is high decrease
c.....the order of p
                  
               elseif ( ((ze_sensor(l).lt.ze_delta).and.
     &                 (qx_sensor(l).lt.qx_delta).and.
     &                 (qy_sensor(l).lt.qy_delta) ).and.
     &                 (iota_sensor(l).lt.iota_delta).and.
     &                 (iota2_sensor(l).lt.iota2_delta).and.
     &                 (pdg_el(l).gt.pl).and.
     &                 ( (it-pcount(l)).gt.plimit ) ) then
                  
                  pdg_el(l) = pdg_el(l) - 1
                  dofs(l) = (pdg_el(l) + 1 )*(pdg_el(l) + 2 ) / 2
                  
                  ze(dofs(l)+1:dofh,l,irk+1) = 0.d0
                  qx(dofs(l)+1:dofh,l,irk+1) = 0.d0
                  qy(dofs(l)+1:dofh,l,irk+1) = 0.d0
                  iota(dofs(l)+1:dofh,l,irk+1) = 0.d0
                  iota2(dofs(l)+1:dofh,l,irk+1) = 0.d0
                  
                  pcount(l) = it
                  
               elseif (pdg_el(l).eq.0.and.
     &                 ( (it-pcount(l)).gt.plimit ) ) then
                  
                  pdg_el(l) = pdg_el(l) + 1 
                  dofs(l) = (pdg_el(l) + 1 )*(pdg_el(l) + 2 ) / 2
                  
                  pcount(l) = it
                  
               endif
#endif


#ifdef SED_LAY

c.....if the sensor is more than the limit and the p is low increase
c.....the order of p
               
               
               if ( ( (temp_qy.ge.qy_delta).or.
     &              (temp_qx.ge.qx_delta).or.
     &              (temp_ze.ge.ze_delta).or.
     &              (temp_tbed.ge.tbed_delta)).and.
     &              (pdg_el(l).lt.ph).and.
     &              ( (it-pcount(l)).gt.plimit ) ) then
                  
                  pdg_el(l) = pdg_el(l) + 1 
                  dofs(l) = (pdg_el(l) +1 )*(pdg_el(l) +2 ) / 2
                  
                  pcount(l) = it
                  
c.....if the sensor is less than the limit and the p is high decrease
c.....the order of p
                  
               elseif ( ((temp_ze.lt.ze_delta).and.
     &                 (temp_qx.lt.qx_delta).and.
     &                 (temp_qy.lt.qy_delta).and.
     &                 (temp_tbed.lt.tbed_delta)).and.
     &                 (pdg_el(l).gt.pl).and.
     &                 ( (it-pcount(l)).gt.plimit ) ) then
                  
                  pdg_el(l) = pdg_el(l) - 1
                  dofs(l) = (pdg_el(l) + 1 )*(pdg_el(l) + 2 ) / 2
                  
                  ze(dofs(l)+1:dofh,l,irk+1) = 0.d0
                  qx(dofs(l)+1:dofh,l,irk+1) = 0.d0
                  qy(dofs(l)+1:dofh,l,irk+1) = 0.d0
                  iota(dofs(l)+1:dofh,l,irk+1) = 0.d0
                  bed(dofs(l)+1:dofh,l,irk+1,:)=0.d0
                  
                  pcount(l) = it
                  
               elseif (pdg_el(l).eq.0.and.
     &                 ( (it-pcount(l)).gt.plimit ) ) then
                  
                  pdg_el(l) = pdg_el(l) + 1 
                  dofs(l) = (pdg_el(l) + 1 )*(pdg_el(l) + 2 ) / 2
                  
                  pcount(l) = it
                  
               endif
#endif
               
            enddo
            
         endif
         
      endif

c..........................
c.... type 2 p-enrichment |
c..........................
      
      if (pflag.eq.2) then

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

         
         kon = pflag2con1       ! usually in [0,1]
         kons = pflag2con2      ! usually 6 works
         lebesgue = lebesguep
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
c$$$            do ll=1,layers
c$$$               tbed_center = tbed_center + bed(k,l,irk+1,ll)*phi_center(k,pa)
c$$$            enddo
            tbed_sensor(l) = 0.d0            
#endif
            
            ze_sensor1 = 0.d0             
            qx_sensor1 = 0.d0
            qy_sensor1 = 0.d0
            
            ze_sensor2 = 0.d0            
            qx_sensor2 = 0.d0
            qy_sensor2 = 0.d0
            
            slimit1 = 0.d0 
            slimit2 = 0.d0
            slimit3 = 0.d0

#ifdef TRACE
            iota_sensor1 = 0.d0

            iota_sensor2 = 0.d0
            
            slimit4 = 0.d0
#endif
            
#ifdef CHEM      
            iota_sensor1 = 0.d0
            iota2_sensor1 = 0.d0
            
            iota_sensor2 = 0.d0
            iota2_sensor2 = 0.d0
            
            slimit4 = 0.d0
            slimit5 = 0.d0
#endif

#ifdef SED_LAY
            tbed_sensor1 = 0.d0

            tbed_sensor2 = 0.d0
            
            slimit6 = 0.d0
#endif

            pa = pdg_el(l)

                                !compute the first sensor type          
            if (pa.eq.0) then
               
               cycle
               
            else

               do k = (pa*(pa+1))/2 + 1 ,dofs(l)
                  do mm = 1,nagp(pa)
                     
                     ze_sensor1 = ze_sensor1 + abs(ze(k,l,irkp+1)* 
     &                    phi_area(k,mm,pa))**lebesgue * wagp(mm,pa)
                     qx_sensor1 = qx_sensor1 + abs(qx(k,l,irkp+1)* 
     &                    phi_area(k,mm,pa))**lebesgue * wagp(mm,pa)
                     qy_sensor1 = qy_sensor1 + abs(qy(k,l,irkp+1)* 
     &                    phi_area(k,mm,pa))**lebesgue * wagp(mm,pa)
                     
#ifdef TRACE       
                     iota_sensor1 = iota_sensor1 + abs(iota(k,l,irkp+1)* 
     &                    phi_area(k,mm,pa))**lebesgue * wagp(mm,pa)
#endif
                     
#ifdef CHEM                
                     iota_sensor1  = iota_sensor1 + abs(iota(k,l,irkp+1)* 
     &                    phi_area(k,mm,pa))**lebesgue * wagp(mm,pa)
                     iota2_sensor1 = iota2_sensor1 + abs(iota2(k,l,irkp+1)* 
     &                    phi_area(k,mm,pa))**lebesgue * wagp(mm,pa)
#endif

#ifdef SED_LAY       

                     do ll=1,layers
                        tbed_sensor1 = tbed_sensor1 + abs(bed(k,l,irkp+1,ll)* 
     &                       phi_area(k,mm,pa))**lebesgue * wagp(mm,pa)
                     enddo
#endif
                     
                  enddo
               enddo
               
            endif
            
                                !compute the second sensor type

            if (pa.eq.0) then
               
               cycle
               
            else
               
               do k = 1, dofs(l)

                  do mm=1,nagp(pa)
                     
                     ze_sensor2 = ze_sensor2 + abs(ze(k,l,irkp+1)* 
     &                    phi_area(k,mm,pa))**lebesgue * wagp(mm,pa)
                     qx_sensor2 = qx_sensor2 + abs(qx(k,l,irkp+1)* 
     &                    phi_area(k,mm,pa))**lebesgue * wagp(mm,pa)
                     qy_sensor2 = qy_sensor2 + abs(qy(k,l,irkp+1)* 
     &                    phi_area(k,mm,pa))**lebesgue * wagp(mm,pa)
                     
#ifdef TRACE
                     iota_sensor2 = iota_sensor2 + abs(iota(k,l,irkp+1)* 
     &                    phi_area(k,mm,pa))**lebesgue * wagp(mm,pa)
#endif
                     
#ifdef CHEM
                     iota_sensor2 = iota_sensor2 + abs(iota(k,l,irkp+1)* 
     &                    phi_area(k,mm,pa))**lebesgue * wagp(mm,pa)
                     iota2_sensor2 = iota2_sensor2 + abs(iota2(k,l,irkp+1)* 
     &                    phi_area(k,mm,pa))**lebesgue * wagp(mm,pa)
#endif

                     
#ifdef SED_LAY
                     do ll=1,layers
                        tbed_sensor2 = tbed_sensor2 + abs(bed(k,l,irkp+1,ll)* 
     &                       phi_area(k,mm,pa))**lebesgue * wagp(mm,pa)
                     enddo
#endif
                     
                  enddo

               enddo
               
            endif

c................................................................

            if (gflag.eq.0) then !dioristic algorithm OFF
               
               if (pa.gt.0) then
                  
                  if  (ze_sensor2.gt.1.0e-12.and.ze_sensor1.gt.1.0e-12 ) then
                     
                     ze_sensor(l) = log10( ( (ze_sensor1 )**roo) / (ze_sensor2**roo) )
                     slimit1 = log10(kon*real(pdg_el(l))**(- lebesgue**2)) - kons
                     
                  else
                     
                     ze_sensor(l) = 1.d0
                     slimit1 = 2.d0	
                     
                  endif

                  if (qx_sensor2.gt.1.0e-12.and.qx_sensor1.gt.1.0e-12 ) then
                     
                     qx_sensor(l) = log10( ( (qx_sensor1)**roo) / (qx_sensor2**roo) )
                     slimit2 = log10(kon*real(pdg_el(l))**(- lebesgue**2)) - kons
                     
                  else
                     
                     qx_sensor(l) = 1.d0
                     slimit2 = 2.d0	
                     
                  endif
                  
                  if ( qy_sensor2.gt.1.0e-12.and.qy_sensor1.gt.1.0e-12 ) then
                     
                     qy_sensor(l) = log10( ((qy_sensor1)**roo) / (qy_sensor2**roo) )
                     slimit3 = log10(kon*real(pdg_el(l))**(- lebesgue**2)) - kons
                     
                  else
                     
                     qy_sensor(l) = 1.d0
                     slimit3 = 2.d0	
                     
                  endif
                  
#ifdef TRACE
                  if ( iota_sensor2.gt.1.0e-12.and.iota_sensor1.gt.1.0e-12 ) then
                     
                     iota_sensor(l) = log10( ( (iota_sensor1)**roo) / (iota_sensor2**roo) )
                     slimit4 = log10(kon*real(pdg_el(l))**(- lebesgue**2)) - kons
                     
                  else
                     
                     iota_sensor(l) = 1.d0
                     slimit4 = 2.d0
                     
                     
                  endif
#endif
                  
#ifdef CHEM            
                  if ( iota_sensor2.gt.1.0e-12.and.iota_sensor1.gt.1.0e-12) then
                     
                     iota_sensor(l) = log10( ( (iota_sensor1 )**roo) / (iota_sensor2**roo) )
                     slimit4 = log10(kon*real(pdg_el(l))**(- lebesgue**2)) - kons
                     
                  else
                     
                     iota_sensor(l) = 1.d0
                     slimit4 = 2.d0
                     
                  endif
                  
                  if ( iota2_sensor2.gt.1.0e-12.and.iota2_sensor1.gt.1.0e-12 ) then
                     
                     iota2_sensor(l) = log10( ( (iota2_sensor1)**roo) / (iota2_sensor2**roo) )
                     slimit5 = log10(kon*real(pdg_el(l))**(- lebesgue**2)) - kons
                     
                  else
                     
                     iota2_sensor(l) = 1.d0
                     slimit5 = 2.d0
                     
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
                  
                  ze_sensor(l) = slimit1
                  qx_sensor(l) = slimit2
                  qy_sensor(l) = slimit3
                  
#ifdef TRACE
                  iota_sensor(l) = slimit4
#endif
                  
#ifdef CHEM            
                  iota_sensor(l) = slimit4
                  iota2_sensor(l) = slimit5
#endif

#ifdef SED_LAY
                  tbed_sensor(l) = slimit6
#endif
                  
               endif

#ifdef TRACE

#elif CHEM

#elif SED_LAY

#else  

c.....if the sensor is less than the limit and the p is low increase
c.....the order of p

               
               if ( ((ze_sensor(l).lt.slimit1).and.
     &              (qx_sensor(l).lt.slimit2).and.
     &              (qy_sensor(l).lt.slimit3)).and.
     &              (pdg_el(l).lt.ph).and.
     &              (it-pcount(l)).ge.plimit )  then
                  
                  pdg_el(l) = pdg_el(l) + 1
                  dofs(l) = (pdg_el(l) + 1 )*(pdg_el(l) +2 ) / 2
                  
                  pcount(l) = it
                  
c.....if the sensor is more than the limit and the p is high decrease
c.....the order of p
                  
               elseif ( ((ze_sensor(l).ge.slimit1).or.
     &                 (qx_sensor(l).ge.slimit2).or.
     &                 (qy_sensor(l).ge.slimit3)).and.
     &                 (pdg_el(l).gt.pl) ) then !.and.
                                !&                    (it-pcount(l)).ge.plimit ) then
                  
                  pdg_el(l) = pdg_el(l) - 1
                  dofs(l) = (pdg_el(l) + 1 )*(pdg_el(l) + 2 ) / 2
                  
                  ze(dofs(l)+1:dofh,l,irkp+1) = 0.d0
                  qx(dofs(l)+1:dofh,l,irkp+1) = 0.d0
                  qy(dofs(l)+1:dofh,l,irkp+1) = 0.d0  

                                !pcount(l) = it 

               endif
#endif
               
#ifdef TRACE  

c.....if the sensor is less than the limit and the p is low increase
c.....the order of p
               
               if ( ( (ze_sensor(l).lt.slimit1).and.
     &              (qx_sensor(l).lt.slimit2).and.
     &              (qy_sensor(l).lt.slimit3).and.
     &              (iota_sensor(l).lt.slimit4)).and.
     &              (pdg_el(l).lt.ph).and.
     &              (it-pcount(l)).ge.plimit ) then
                  
                  pdg_el(l) = pdg_el(l) + 1
                  dofs(l) = (pdg_el(l) + 1 )*(pdg_el(l) +2 ) / 2
                  
                  pcount(l) = it
                  
c.....if the sensor is more than the limit and the p is high decrease
c.....the order of p
                  
               elseif ( ((ze_sensor(l).ge.slimit1).or.
     &                 (qx_sensor(l).ge.slimit2).or.
     &                 (qy_sensor(l).ge.slimit3).or.
     &                 (iota_sensor(l).ge.slimit4)).and.
     &                 (pdg_el(l).gt.pl) ) then !.and.
                                !&                    (it-pcount(l)).ge.plimit ) then
                  
                  pdg_el(l) = pdg_el(l) - 1
                  dofs(l) = (pdg_el(l) + 1 )*(pdg_el(l) + 2 ) / 2
                  
                  ze(dofs(l)+1:dofh,l,irkp+1) = 0.d0
                  iota(dofs(l)+1:dofh,l,irkp+1) = 0.d0
                  qx(dofs(l)+1:dofh,l,irkp+1) = 0.d0
                  qy(dofs(l)+1:dofh,l,irkp+1) = 0.d0
                  
                                !pcount(l) = it
                  
               endif
#endif

#ifdef CHEM 
               
c.....if the sensor is less than the limit and the p is low increase
c.....the order of p
               
               
               if ( ( (ze_sensor(l).lt.slimit1).and.
     &              (qx_sensor(l).lt.slimit2).and.
     &              (qy_sensor(l).lt.slimit3).and.
     &              (iota_sensor(l).lt.slimit4).and.
     &              (iota2_sensor(l).lt.slimit5)).and. 
     &              (pdg_el(l).lt.ph).and. 
     &              (it-pcount(l)).ge.plimit) then
                  
                  pdg_el(l) = pdg_el(l) + 1
                  dofs(l) = (pdg_el(l) + 1 )*(pdg_el(l) +2 ) / 2
                  
                  pcount(l) = it
                  
c.....if the sensor is more than the limit and the p is high decrease
c.....the order of p
                  
               elseif ( ((ze_sensor(l).ge.slimit1).or.
     &                 (qx_sensor(l).ge.slimit2).or.
     &                 (qy_sensor(l).ge.slimit3).or.
     &                 (iota_sensor(l).ge.slimit4).or.
     &                 (iota2_sensor(l).ge.slimit5)).and.
     &                 (pdg_el(l).gt.pl) ) then !.and.
                                !&                    (it-pcount(l)).ge.plimit ) then
                  
                  pdg_el(l) = pdg_el(l) - 1
                  dofs(l) = (pdg_el(l) + 1 )*(pdg_el(l) + 2 ) / 2
                  
                  ze(dofs(l)+1:dofh,l,irkp+1) = 0.d0
                  iota(dofs(l)+1:dofh,l,irkp+1) = 0.d0
                  iota2(dofs(l)+1:dofh,l,irkp+1) = 0.d0
                  qx(dofs(l)+1:dofh,l,irkp+1) = 0.d0
                  qy(dofs(l)+1:dofh,l,irkp+1) = 0.d0

                                !pcount(l) = it
                  
               endif
#endif

#ifdef SED_LAY 

c.....if the sensor is less than the limit and the p is low increase
c.....the order of p
               
               if ( ( (ze_sensor(l).lt.slimit1).and.
     &              (qx_sensor(l).lt.slimit2).and.
     &              (qy_sensor(l).lt.slimit3).and.
     &              (tbed_sensor(l).lt.slimit6)).and.
     &              (pdg_el(l).lt.ph).and.
     &              (it-pcount(l)).ge.plimit ) then
                  
                  pdg_el(l) = pdg_el(l) + 1
                  dofs(l) = (pdg_el(l) + 1 )*(pdg_el(l) +2 ) / 2
                  
                  pcount(l) = it
                  
c.....if the sensor is more than the limit and the p is high decrease
c.....the order of p
                  
               elseif ( ((ze_sensor(l).ge.slimit1).or.
     &                 (qx_sensor(l).ge.slimit2).or.
     &                 (qy_sensor(l).ge.slimit3).or.
     &                 (tbed_sensor(l).ge.slimit6)).and.
     &                 (pdg_el(l).gt.pl) ) then !.and.
                                !&                    (it-pcount(l)).ge.plimit ) then
                  
                  pdg_el(l) = pdg_el(l) - 1
                  dofs(l) = (pdg_el(l) + 1 )*(pdg_el(l) + 2 ) / 2
                  
                  ze(dofs(l)+1:dofh,l,irkp+1) = 0.d0
                  bed(dofs(l)+1:dofh,l,irkp+1,:) = 0.d0
                  qx(dofs(l)+1:dofh,l,irkp+1) = 0.d0
                  qy(dofs(l)+1:dofh,l,irkp+1) = 0.d0
                  
                                !pcount(l) = it
                  
               endif
#endif
               
            endif
c................................................................

            if (gflag.eq.1) then ! dioristic algorithm ON

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


         if (gflag.eq.1) then   !if dioristic algorithm ON

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
                  
                  avg_zesen = sum(ze_sensor)/MNES
                  avg_qxsen = sum(qx_sensor)/MNES
                  avg_qysen = sum(qy_sensor)/MNES
                  
#ifdef TRACE
                  avg_iotasen = sum(iota_sensor)/MNES
#endif
                  
#ifdef CHEM             
                  avg_iotasen = sum(iota_sensor)/MNES
                  avg_iota2sen = sum(iota2_sensor)/MNES
#endif

#ifdef SED_LAY
                  avg_tbedsen = sum(tbed_sensor)/MNES
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
            
            ze_delta = ( diorism / 100.d0 )* (abs(maxze_sensor - minze_sensor) )
            qx_delta = ( diorism / 100.d0 )* (abs(maxqx_sensor - minqx_sensor) )
            qy_delta = ( diorism / 100.d0 )* (abs(maxqy_sensor - minqy_sensor) )   
            
#ifdef TRACE
            iota_delta = ( diorism / 100.d0 )* (abs(maxiota_sensor - miniota_sensor) )
#endif
            
#ifdef CHEM
            iota_delta = ( diorism / 100.d0 )* (abs(maxiota_sensor - miniota_sensor) )
            iota2_delta = ( diorism / 100.d0 )* (abs(maxiota2_sensor - miniota2_sensor) )
#endif

#ifdef SED_LAY
            tbed_delta = ( diorism / 100.d0 )* (abs(maxtbed_sensor - mintbed_sensor) )
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
               
c.....if the sensor is less than delta and the p is low increase
c.....the order of p
               
               if ( ( (temp_qy.lt.qy_delta).and.
     &              (temp_qx.lt.qx_delta).and.
     &              (temp_ze.lt.ze_delta) ).and.
     &              (pdg_el(l).lt.ph).and.
     &              (it-pcount(l)).ge.plimit ) then
                  
                  pdg_el(l) = pdg_el(l) + 1 
                  dofs(l) = (pdg_el(l) +1 )*(pdg_el(l) +2 ) / 2
                  
                  pcount(l) = it
                  
c.....if the sensor is more than delta and the p is high decrease
c.....the order of p
                  
               elseif ( ((temp_ze.ge.ze_delta).or.
     &                 (temp_qx.ge.qx_delta).or.
     &                 (temp_qy.ge.qy_delta) ).and.
     &                 (pdg_el(l).gt.pl) ) then !.and.
                                !&                    (it-pcount(l)).ge.plimit ) then
                  
                  pdg_el(l) = pdg_el(l) - 1
                  dofs(l) = (pdg_el(l) + 1 )*(pdg_el(l) + 2 ) / 2
                  
                  ze(dofs(l)+1:dofh,l,irkp+1) = 0.d0
                  qx(dofs(l)+1:dofh,l,irkp+1) = 0.d0
                  qy(dofs(l)+1:dofh,l,irkp+1) = 0.d0

                                !pcount(l) = it
                  
               endif
#endif

#ifdef TRACE

c.....if the sensor is less than delta and the p is low increase
c.....the order of p

               if ( ((temp_qy.lt.qy_delta).and.
     &              (temp_qx.lt.qx_delta).and.
     &              (temp_ze.lt.ze_delta).and.
     &              (temp_iota.lt.iota_delta)).and.
     &              (pdg_el(l).lt.ph).and.
     &              (it-pcount(l)).ge.plimit ) then

                  pdg_el(l) = pdg_el(l) + 1 
                  dofs(l) = (pdg_el(l) +1 )*(pdg_el(l) +2 ) / 2

                  pcount(l) = it
                  
c.....if the sensor is less than delta and the p is high decrease
c.....the order of p

               elseif (((temp_ze.ge.ze_delta).or.
     &                 (temp_qx.ge.qx_delta).or.
     &                 (temp_qy.ge.qy_delta).or.
     &                 (temp_iota.ge.iota_delta)).and.
     &                 (pdg_el(l).gt.pl) ) then !.and.
                                !&                    (it-pcount(l)).ge.plimit) then

                  pdg_el(l) = pdg_el(l) - 1
                  dofs(l) = (pdg_el(l) + 1 )*(pdg_el(l) + 2 ) / 2
                  
                  ze(dofs(l)+1:dofh,l,irkp+1) = 0.d0
                  qx(dofs(l)+1:dofh,l,irkp+1) = 0.d0
                  qy(dofs(l)+1:dofh,l,irkp+1) = 0.d0
                  iota(dofs(l)+1:dofh,l,irkp+1) = 0.d0

                                !pcount(l) = it
                  
               endif
#endif

#ifdef CHEM           

c.....if the sensor is less than delta and the p is low increase
c.....the order of p

               if ( ((temp_qy.lt.qy_delta).and.
     &              (temp_qx.lt.qx_delta).and.
     &              (temp_ze.lt.ze_delta).and.
     &              (temp_iota.lt.iota_delta).and.
     &              (temp_iota2.lt.iota2_delta)).and.
     &              (pdg_el(l).lt.ph).and.
     &              (it-pcount(l)).ge.plimit ) then

                  pdg_el(l) = pdg_el(l) + 1 
                  dofs(l) = (pdg_el(l) +1 )*(pdg_el(l) +2 ) / 2

                  pcount(l) = it
                  
c.....if the sensor is less than delta and the p is high decrease
c.....the order of p

               elseif (((temp_ze.ge.ze_delta).or.
     &                 (temp_qx.ge.qx_delta).or.
     &                 (temp_qy.ge.qy_delta).or.
     &                 (temp_iota.ge.iota_delta).or.
     &                 (temp_iota2.ge.iota2_delta)).and.
     &                 (pdg_el(l).gt.pl) ) then !.and.
                                !&                    (it-pcount(l)).ge.plimit) then

                  pdg_el(l) = pdg_el(l) - 1
                  dofs(l) = (pdg_el(l) + 1 )*(pdg_el(l) + 2 ) / 2
                  
                  ze(dofs(l)+1:dofh,l,irkp+1) = 0.d0
                  qx(dofs(l)+1:dofh,l,irkp+1) = 0.d0
                  qy(dofs(l)+1:dofh,l,irkp+1) = 0.d0
                  iota(dofs(l)+1:dofh,l,irkp+1) = 0.d0
                  iota2(dofs(l)+1:dofh,l,irkp+1) = 0.d0

                                !pcount(l) = it
                  
               endif         
#endif

#ifdef SED_LAY

c.....if the sensor is less than delta and the p is low increase
c.....the order of p

               if ( ((temp_qy.lt.qy_delta).and.
     &              (temp_qx.lt.qx_delta).and.
     &              (temp_ze.lt.ze_delta).and.
     &              (temp_tbed.lt.tbed_delta)).and.
     &              (pdg_el(l).lt.ph).and.
     &              (it-pcount(l)).ge.plimit ) then

                  pdg_el(l) = pdg_el(l) + 1 
                  dofs(l) = (pdg_el(l) +1 )*(pdg_el(l) +2 ) / 2

                  pcount(l) = it
                  
c.....if the sensor is less than delta and the p is high decrease
c.....the order of p

               elseif (((temp_ze.ge.ze_delta).or.
     &                 (temp_qx.ge.qx_delta).or.
     &                 (temp_qy.ge.qy_delta).or.
     &                 (temp_tbed.ge.tbed_delta)).and.
     &                 (pdg_el(l).gt.pl) ) then !.and.
                                !&                    (it-pcount(l)).ge.plimit) then

                  pdg_el(l) = pdg_el(l) - 1
                  dofs(l) = (pdg_el(l) + 1 )*(pdg_el(l) + 2 ) / 2
                  
                  ze(dofs(l)+1:dofh,l,irkp+1) = 0.d0
                  qx(dofs(l)+1:dofh,l,irkp+1) = 0.d0
                  qy(dofs(l)+1:dofh,l,irkp+1) = 0.d0
                  bed(dofs(l)+1:dofh,l,irkp+1,:) = 0.d0

                                !pcount(l) = it
                  
               endif
#endif


            enddo
            
         endif

      endif

c.....deal with the interior barrier problem, by making the front and back elements agree
c.....at the min of the union for stability

      if (nibseg.gt.0) then

         do j = 1,nibseg

            if ( pdg_el(nedel(1,nibsegn(1,j))).ne.pdg_el(nedel(1,nibsegn(2,j))) ) then

               pdg_el(nedel(1,nibsegn(1,j))) = max(pdg_el(nedel(1,nibsegn(1,j))), pdg_el(nedel(1,nibsegn(2,j))) )
               pdg_el(nedel(2,nibsegn(1,j))) = max(pdg_el(nedel(1,nibsegn(1,j))), pdg_el(nedel(1,nibsegn(2,j))) )

               dofs(nedel(1,nibsegn(1,j))) = max(dofs(nedel(1,nibsegn(1,j))), dofs(nedel(1,nibsegn(2,j))) )
               dofs(nedel(2,nibsegn(1,j))) = max(dofs(nedel(1,nibsegn(1,j))), dofs(nedel(1,nibsegn(2,j))) )

            endif

         enddo

      endif
      

      end subroutine p_enrichment


  

