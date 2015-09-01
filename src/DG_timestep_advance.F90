      SUBROUTINE DG_TIMESTEP_ADVANCE(s,dg_here,global_here,nodalattr_here,IT)
      
!.....Use appropriate modules
      USE SIZES
      USE GLOBAL
      USE DG
      USE NodalAttributes      
      
      type (sizes_type) :: s
      type (dg_type) :: dg_here
      type (global_type) :: global_here
      type (nodalattr_type) :: nodalattr_here

      INTEGER IT,K,J,KK,I,MM 
      REAL(SZ) s_dg,alph,eta,seta
      Real(sz) sigma      

!.......Update variables for next time step

#ifdef FILTER 
        !Stabilization filtering 
        DO J = 1,global_here%NE
           DO K = 2,dg_here%DOFS(J)

              alph = 10.D0
              s_dg = dg_here%RK_stage

              if (k.le.3) then
                 eta = 1.D0 / ( global_here%pdg_el(j)+1.D0 )
              elseif (k.gt.3.and.k.le.6) then
                 eta = 2.D0 / ( global_here%pdg_el(j)+1.D0 )
              elseif (k.gt.6.and.k.le.10) then
                 eta = 3.D0 / ( global_here%pdg_el(j)+1.D0 )
              elseif (k.gt.11.and.k.le.15) then
                 eta = 4.D0 / ( global_here%pdg_el(j)+1.D0 )
              elseif (k.gt.16.and.k.le.21) then
                 eta = 5.D0 / ( global_here%pdg_el(j)+1.D0 )       
              elseif (k.gt.22.and.k.le.28) then
                 eta = 6.D0 / ( global_here%pdg_el(j)+1.D0 )   
              elseif (k.gt.29.and.k.le.36) then
                 eta = 7.D0 / ( global_here%pdg_el(j)+1.D0 )               
              endif

              seta = eta**s_dg
              sigma = exp(-alph*seta)

!.........Filter the layers if flagged

#ifdef SED_LAY
              do l = 1,layers
                
                 dg_here%bed(K,J,dg_here%NRK+1,l) = sigma * dg_here%bed(K,J,dg_here%nrk+1,l) 
               
              enddo
#endif

              dg_here%ZE(K,J,dg_here%NRK+1) = sigma*dg_here%ZE(K,J,dg_here%nrk+1) 
              dg_here%QX(K,J,dg_here%NRK+1) = sigma*dg_here%QX(K,J,dg_here%nrk+1)
              dg_here%QY(K,J,dg_here%NRK+1) = sigma*dg_here%QY(K,J,dg_here%nrk+1)


!.........Filter the transported global_here%tracer term if flagged

#ifdef TRACE
              dg_here%iota(K,J,dg_here%NRK+1) = sigma*dg_here%iota(K,J,dg_here%nrk+1)

#endif

!........Filter chemistry transported terms if flagged

#ifdef CHEM
              dg_here%iota(K,J,dg_here%NRK+1) =  sigma*dg_here%iota(K,J,dg_here%nrk+1) 
              dg_here%iota2(K,J,dg_here%NRK+1) = sigma*dg_here%iota2(K,J,dg_here%nrk+1) 

#endif

           ENDDO
        ENDDO
#endif
      

        DO J = 1,S%MNE
           
           DO K = 1,dg_here%DOFS(J)
              dg_here%ZE(K,J,1) = dg_here%ZE(K,J,dg_here%NRK+1)
              dg_here%QX(K,J,1) = dg_here%QX(K,J,dg_here%NRK+1)
              dg_here%QY(K,J,1) = dg_here%QY(K,J,dg_here%NRK+1)
#ifdef TRACE
              dg_here%iota(K,J,1) = dg_here%iota(K,J,dg_here%NRK+1)
#endif

#ifdef CHEM
              dg_here%iota(K,J,1) = dg_here%iota(K,J,dg_here%NRK+1)
              dg_here%iota2(K,J,1) = dg_here%iota2(K,J,dg_here%NRK+1)
#endif

#ifdef DYNP
              dg_here%dynP(K,J,1) = dg_here%dynP(K,J,dg_here%NRK+1)
#endif

#ifdef SED_LAY
              dg_here%bed(K,J,1,:) = dg_here%bed(K,J,dg_here%NRK+1,:)
#endif
           ENDDO
        ENDDO

        DO KK = 2,dg_here%NRK+1
           DO J = 1,S%MNE
              DO K = 1,dg_here%DOFH
                 dg_here%ZE(K,J,KK) = 0.D0
                 dg_here%QX(K,J,KK) = 0.D0
                 dg_here%QY(K,J,KK) = 0.D0
                 dg_here%RHS_ZE(K,J,KK-1) = 0.D0
                 dg_here%RHS_QX(K,J,KK-1) = 0.D0
                 dg_here%RHS_QY(K,J,KK-1) = 0.D0
#ifdef TRACE
                 dg_here%iota(K,J,KK) = 0.D0
                 dg_here%RHS_iota(K,J,KK-1) = 0.D0
#endif

#ifdef CHEM
                 dg_here%iota(K,J,KK) = 0.D0
                 dg_here%iota2(K,J,KK) = 0.D0
                 dg_here%RHS_iota(K,J,KK-1) = 0.D0
                 dg_here%RHS_iota2(K,J,KK-1) = 0.D0
#endif

#ifdef DYNP
                 dg_here%dynP(K,J,KK) = 0.D0
                 dg_here%RHS_dynP(K,J,KK-1) = 0.D0
#endif

#ifdef SED_LAY
                 dg_here%bed(K,J,KK,:) = 0.D0
                 dg_here%hb(K,J,KK) = 0.D0
                 dg_here%RHS_bed(K,J,KK-1,:) = 0.D0
#endif


              ENDDO
           ENDDO
        ENDDO
        
#ifdef CMPI
#ifdef BLKOUT
        CALL MSG_BLOCKSYNC_START()
#endif
#endif

        CALL WRITE_RESULTS(s,dg_here,global_here,IT,.FALSE.)
        
#ifdef CMPI
#ifdef BLKOUT
        CALL MSG_BLOCKSYNC_FINISH()
#endif
#endif

        CALL SCRUTINIZE_SOLUTION(s,dg_here,global_here,IT)        
      
       END SUBROUTINE
