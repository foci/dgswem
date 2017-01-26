     SUBROUTINE ptr_arrays()
     
     USE globals, ONLY: nied,mnqpte, &
                       Hi,He,Qxi,Qxe,Qyi,Qye, &
                       xmi,xme,ymi,yme,xymi,xyme, &
                       Hfi,Hfe,Qxfi,Qxfe,Qyfi,Qyfe, &
                       inx,iny,detJe_in,detJe_ex,const, &
                       Hhatv,Qxhatv,Qyhatv, &

     IMPLICIT NONE
     INTEGER :: alloc_status

     ALLOCATE(Hi(nied,mnqpte),He(nied,mnqpte),Qxi(nied,mnqpte),Qxe(nied,mnqpte),Qyi(nied,mnqpte),Qye(nied,mnqpte),STAT=alloc_status)
     IF(alloc_status /= 0) PRINT*, "Allocation error: Hi,He,Qxi,Qxe,Qyi,Qye"

     ALLOCATE(xmi(nied,mnqpte),xme(nied,mnqpte),ymi(nied,mnqpte),yme(nied,mnqpte),xymi(nied,mnqpte),xyme(nied,mnqpte),STAT=alloc_status)
     IF(alloc_status /= 0) PRINT*, "Allocation error: xmi,xme,ymi,yme,xymi,xyme"       

     ALLOCATE(Hfi(nied,mnqpte),Hfe(nied,mnqpte),Qxfi(nied,mnqpte),Qxfe(nied,mnqpte),Qyfi(nied,mnqpte),Qyfe(nied,mnqpte),STAT=alloc_status)
     IF(alloc_status /= 0) PRINT*, "Allocation error: Hfi,Hfe,Qxfi,Qxfe,Qyfi,Qyfe"
       
     ALLOCATE(const(nied),inx(nied,mnqpte),iny(nied,mnqpte),detJe_in(nied,mnqpte),detJe_ex(nied,mnqpte),STAT=alloc_status)
     IF(alloc_status /= 0) PRINT*, "Allocation error: const,inx,iny,detJe_in,detJe_ex"       

     ALLOCATE(Hhatv(nied),Qxhatv(nied),Qyhatv(nied),STAT=alloc_status)
     IF(alloc_status /= 0) PRINT*, "Allocation error: Hhatv,Qxhatv,Qyhatv"


     RETURN
     END SUBROUTINE ptr_arrays
