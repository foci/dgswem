      SUBROUTINE connect()

      USE globals, ONLY: rp,nn,ne,ned,ect,el_type,nverts, &
                         mnepn,epn,nepn, &
                         ged2nn,ged2el,ged2led, &
                         nbed,bedn,nied,iedn, &
                         nobed,obedn,nfbed,fbedn,nnfbed,nfbedn,nfbednn, &
                         nope,obseg,obnds,nbou,fbseg,fbnds, &
                         ed_type,recv_edge                      
                         
      USE messenger2, ONLY: dirname,lname,myrank,nred,redn      
      USE read_dginp, ONLY: out_direc
      USE edge_connectivity_mod, ONLY: elements_per_node,find_edge_pairs,find_interior_edges, &
                                       find_open_edges,find_flow_edges,find_recieve_edges, &
                                       print_connect_info
                         

      IMPLICIT NONE
      INTEGER :: i
      
                
      CALL elements_per_node(ne,nn,nverts,el_type,ect,nepn,mnepn,epn) 
      
      CALL find_edge_pairs(ne,nverts,el_type,ect,nepn,epn,ned,ged2el,ged2nn,ged2led)
      
      CALL find_interior_edges(ned,ged2el,nied,iedn,ed_type,recv_edge,nbed,bedn)
      
      CALL find_open_edges(nope,obseg,obnds,ged2nn,nobed,obedn,ed_type,recv_edge)      
      
      CALL find_flow_edges(nbou,fbseg,fbnds,ged2nn,nnfbed,nfbedn,nfbednn,nfbed,fbedn,recv_edge,ed_type)
      
      CALL find_recieve_edges(ned,recv_edge,nred,redn,ed_type)
      
      IF (myrank == 0) THEN      
        CALL print_connect_info(mnepn,ned,nied,nobed,nfbed,nnfbed,nred)
      ENDIF

      
!       OPEN(unit=17,file=trim(out_direc) // 'connect.d')      
! 
!       
!       WRITE(17,"(A)") "---------------------------------------------"
!       WRITE(17,"(A)") "       Edge Connectivity Information         "
!       WRITE(17,"(A)") "---------------------------------------------"
!       WRITE(17,"(A)") " "        
!       WRITE(17,"(A,I7)") '   maximum elements per node:', mnepn
!       WRITE(17,"(A)") ' '            
!       WRITE(17,"(A,I7)") '   number of total edges:', ned
!       WRITE(17,"(A)") ' '  
!       WRITE(17,"(A,I7)") '   number of interior edges:', nied
!       WRITE(17,"(A)") ' '        
!       WRITE(17,"(A,I7)") '   number of open boundary edges:', nobed
!       WRITE(17,"(A)") ' '                
!       WRITE(17,"(A,I7)") '   number of specified normal boundary edges:', nfbed
!       WRITE(17,"(A,I7)") '   number of no normal flow boundary edges:', nnfbed
!       WRITE(17,"(A)") ' '        
!       WRITE(17,"(A,I7)") '   number of recieve edges:', nred
!       WRITE(17,"(A)") ' ' 
!       WRITE(17,"(A,I7)") '   number of missing edges:',ned-(nied+nobed+nfbed+nnfbed+nred)
!       WRITE(17,"(A)") ' '      
!       
!       CLOSE(17)

! !       write edge connectivity information in similar format to fort.17
!       OPEN(UNIT=17,FILE=dirname(1:lname)//'/'//'fort.17')
!       
!       WRITE(17,*) ned
!       DO i = 1,ned
!         WRITE(17,*) i,ged2nn(1,i),ged2nn(2,i),ged2el(1,i),ged2el(2,i)
!       ENDDO
! 
!       WRITE(17,*) 'number of interior edges:', nied
!       DO i = 1,nied
!         WRITE(17,*) i,iedn(i),ged2nn(1,iedn(i)),ged2nn(2,iedn(i))
!       ENDDO
! 
!       WRITE(17,*) 'number of no normal flow boundary edges:', nnfbed
!       DO i = 1,nnfbed
!         WRITE(17,*) i,nfbedn(i),ged2nn(1,nfbedn(i)),ged2nn(2,nfbedn(i))
!       ENDDO
! 
!       WRITE(17,*) 'number of open boundary edges:', nobed
!       DO i = 1,nobed
!         WRITE(17,*) i,obedn(i),ged2nn(1,obedn(i)),ged2nn(2,obedn(i))
!       ENDDO
! 
!       WRITE(17,*) 'number of flow specified boundary edges:', nfbed
!       DO i = 1,nfbed
!         WRITE(17,*) i,fbedn(i),ged2nn(1,fbedn(i)),ged2nn(2,fbedn(i))
!       ENDDO
! 
!       WRITE(17,*) "global to local edge table"
!       DO i = 1,ned
!         WRITE(17,*) i,ged2led(1,i),ged2led(2,i)
!       ENDDO
!       
!       CLOSE(17)

      RETURN
      END SUBROUTINE connect
