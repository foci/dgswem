      subroutine reset_lts()

      USE GLOBAL
      USE DG
#ifdef CMPI
      USE MESSENGER_NEW
      USE MESSENGER_NODE_NEW
#endif
c
      implicit none
      integer ged,led,gedb,gedf,elb,elf,i,j,k,kk,l,idtemp
      real*8 qxedge,qnormal,qyedge,zedge,hedge,dtemp(2),ddummy(2),
     $     ufluxmax,iprocmax,slamdeg,sfeadeg,windmax
      logical dtmod

c      slammin=-99.
c      slammax=-85.
c      sfeamin=17.
c      sfeamax=30.
      dtmod=.false.
      Nelem(:) = 0
      if (callresetflag.eq.1) go to 9000
      if (dtcheck.eq.-1) go to 7777
      ufluxmax=0.

      DO j=1,NE

C.....Obtain the global and local edge numbers

         pa = pdg_el(j)
         

#ifdef P0
         if (pa.eq.0) then
            pa = 1
         endif
#endif

c         ufluxold(j)=uflux(j)
         uflux(j)=0
         do kk=1,3
            ged=neled(kk,j)
            led=nedsd(1,ged)
            qnormal=0
            zedge=0
            do i=1,negp(pa)
               qxedge=0
               qyedge=0
               do k=1,dofs(j)
                  qxedge=qxedge+qx(k,j,1)*phi_edge(k,i,led,pa)
                  qyedge=qyedge+qy(k,j,1)*phi_edge(k,i,led,pa)
                  zedge=zedge+ze(k,j,1)*phi_edge(k,i,led,pa)
               enddo
               qnormal=qnormal+qxedge*cosnx(ged)+qyedge*sinnx(ged)
            enddo
            n1=nedno(1,ged)
            n2=nedno(2,ged)
            hedge=.5*(dp(n1)+dp(n2))+zedge/negp(pa)
            hedge=max(hedge,h0)
            uflux(j)=uflux(j)+(qnormal/negp(pa))/hedge
         enddo
         uflux(j)=uflux(j)/3
         if (abs(uflux(j)).gt.ufluxmax) then
            ufluxmax=abs(uflux(j))
         endif
c         if (myproc.eq.22.and.j.eq.8801)  then
c            write(5555,*) 'uflux ',j,uflux(j) 
c         endif
      enddo
c
      if (LTSmodflag.eq.1) then
         windmax=0.0
         do i=1,np
            WINDX = WSX2(I)
            WINDY = WSY2(I)
            WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
            windmax=max(windmax,windmag)
         enddo
#ifdef CMPI      
c      dtemp(1)=ufluxmax
c      call FindMax(dtemp(1),ddummy(1))
c      ufluxmax=ddummy(1)
         dtemp(1)=windmax
         call FindMax(dtemp(1),ddummy(1))
         windmax=ddummy(1)
#endif
         if (myproc.eq.0) 
     $        write(*,*) 'maximums',windmax
      endif
      CALL cfl(NE,NP,ics,X,Y,dp,NM,SLAM0,SFEA0,DtLoc,slopeflag,myproc)
c      if (myproc.eq.22) write(5555,*) 'dtloc',dtloc(8801)
c      dt_max = MAXVAL(Dtloc)
c      dt_min = MINVAL(Dtloc)
c      dt_min = max(dt_min,dt_min0)
c      Mint = dt_max/dt_min;
c      WRITE(*,*) "Before Scale ---------------------------------------"
c      WRITE(*,*) "Largest stable time-step over grid from CFL: ",dt_max
c      WRITE(*,*) "Smallest stable time-step over grid from CFL: ",dt_min
c      WRITE(*,*) "Largest to smallest ratio: ",Mint


      ! Scale Local time-steps to user input
      !IF (dtmax > DTDP) THEN
c      write(925+myproc,*) 'updating time step'
      do i=1,NE
         TsGroupSave(i)=TsGroup(i)
      enddo
c
      dt_min0=max(dt_min0,DtDom(1))
      DO i=1,NE
!        DtLoc(i) = min((DtLoc(i)/dt_max)*DT,dtloc(i))
c        DtLoc(i) = min(DTDP,dtloc(i))
        DtLoc(i)=max(dt_min0,dtloc(i))
      ENDDO
c      dt_max = MAXVAL(Dtloc);
c      dt_min = MINVAL(Dtloc);
c      dt_min = max(dt_min,dt_min0)
c      Mint = dt_max/dt_min;
c      WRITE(*,*) "After Scale ---------------------------------------"
c      WRITE(*,*) "Largest stable time-step over grid from CFL: ",dt_max
c      WRITE(*,*) "Smallest stable time-step over grid from CFL: ",dt_min
c      WRITE(*,*) "Largest to smallest ratio: ",Mint,Mbar

!  CND
c      if (lts_flag /= 1) then
c         Mint=1
c      endif
!  
c      if (Mint==1) then
c         Ndom=1
c         DtDom(1)=dt_min
c         Msteps(1)=1
c      else
!  hardwired for at most 2 domains
c         Ndom=MaxNdom
!         Msteps(1) = dt_max/dt_min;
!         Msteps(2) = 1;
!         DtDom(1) = dt_max/Msteps(1);
!         DtDom(2) = dt_max;
c         do i=1,Ndom-1
c            Msteps(i)=Mint**(1./(Ndom-1))
c            Msteps(i)=min(Msteps(i),Mbar)
c         enddo
c         Msteps(Ndom)=1
c         write(*,*) (Msteps(i),i=1,Ndom)
c         if (Msteps(1)<2) then
c            Ndom=1
c            Dtdom(1)=dt_min
c         else
c            DtDom(1) = dt_min
c            do i=2,Ndom+1
c               DtDom(i)=Msteps(i-1)*dtdom(i-1)
c            enddo
c         endif
c      endif
!!!!
c      Nelem(:) = 0
      TsGroup(:)=1
      if (nDom > 1) then
         DO i=1,NE
            DO j=2,nDom
               IF ( (DtLoc(i)>=DtDom(j-1)) .and. (DtLoc(i)<DtDom(j)) ) THEN
                  DtLoc(i) = DtDom(j-1)
                  TsGroup(i) = j-1
c                  Nelem(j-1) = Nelem(j-1) + 1
c                  ElemNum(j-1,Nelem(j-1)) = i
               Endif
            enddo
            if (Dtloc(i) >= DtDom(Ndom)) then
               dtloc(i)=dtdom(Ndom)
               TsGroup(i)=Ndom
c               Nelem(Ndom) = Nelem(Ndom) + 1
c               ElemNum(Ndom,Nelem(Ndom)) = i
            endif
         ENDDO
      else
         DO i=1,NE
            DtLoc(i) = DtDom(1)
            TsGroup(i) = 1
          enddo
      endif
c
c      do kk=1,3
c         do l=1,nibseg
c            GEDB = NIBSEGN(1,L)
c            GEDF = NIBSEGN(2,L)
c            ELB = NEDEL(1,GEDB)
c            ELF = NEDEL(1,GEDF)
c            if (tsgroup(elb).ne.tsgroup(elf) ) then
c               tsgroup(elb)=min(tsgroup(elb),tsgroup(elf))
c               tsgroup(elf)=min(tsgroup(elb),tsgroup(elf))
c         tsgroup(elb)=1
c         tsgroup(elf)=1
c            endif
c         enddo
c      enddo
c  
c  try some new stuff
c      do i=1,ne
c         if (TsGroup(i).gt.TsGroupSave(i)+1) then
c            TsGroup(i)=TsGroupSave(i)+1
c         endif
c         if (abs(uflux(i)-ufluxold(i)).gt..1*abs(ufluxold(i))) then
c            TsGroup(i)=1
c         endif
c      enddo
c      do i=1,ne
c         do j=1,3
c            n1=nm(i,j)
c            slamdeg=slam(n1)/deg2rad
c            sfeadeg=sfea(n1)/deg2rad
c            if (slamdeg.ge.slammin.and.slamdeg.lt.slammax) then
c               if (sfeadeg.ge.sfeamin.and.sfeadeg.lt.sfeamax) then
c                  tsgroup(i)=1
c               endif
c            endif
c         enddo
c      enddo
      if (LTSmodflag.eq.1) then
         do k=1,nmodtimes
            if (myproc.eq.0) write(*,*) ' in reset ',time_a,k,tmod1(k),tmod2(k)
            if (time_a .ge. tmod1(k) .and. time_a .lt. tmod2(k)) go to 4802
         enddo
      endif
 4802 continue

      do i=1,ne
         if (LTSmodflag.eq.1) then
            windmag=0.0
            do j=1,3
               n1=nm(i,j)
               WINDX = WSX2(n1)
               WINDY = WSY2(n1)
               WINDMAG = windmag+SQRT(WINDX*WINDX+WINDY*WINDY)
            enddo
            windmag=windmag/3
            if (windmag.gt..2*windmax) then
               if (TsGroup(i).gt.TsGroupmod(k)) TsGroup(i)=TsGroupmod(k)
            endif
         endif
         if (TsGroup(i).gt.TsGroupSave(i)+1) then
            TsGroup(i)=TsGroupSave(i)+1
         endif
      enddo
c  check internal barriers
      do l=1,nibseg
         GEDB = NIBSEGN(1,L)
         GEDF = NIBSEGN(2,L)
c
         if (gedf.eq.0) go to 8899
         ELB = NEDEL(1,GEDB)
         ELF = NEDEL(1,GEDF)
c         
c         tsgroup(elb)=min(tsgroup(elb),tsgroupsave(elb))
c         tsgroup(elf)=min(tsgroup(elf),tsgroupsave(elf))
         tsgroup(elb)=1
         tsgroup(elf)=1
 8899    continue
      enddo
 9000 continue
      do i=1,ne
cnd try something
         if (TSGROUPFLAG(i)) then
            TsGroup(i)=1
         endif
      enddo
cnd
         

#ifdef CMPI
      call updatei_elem(tsgroup,tsgroup,1)
#endif

c
c      do i=1,ne
c         if (TsGroupSave(i).ne.TsGroup(i)) then
c            write(987,*) myproc,i,
c     $           TsGroup(i),TsGroupSave(i)
c            dtmod=.true.
c         endif
c      enddo
c      if (dtmod) then
c         write(987,*) time_a/86400.
c         write(987,*)
c      endif

c
      go to 7788
 7777 do i=1,ne
         TsGroup(i)=1
      enddo
      Nelem(:) = 0
      if (myproc.eq.0) write(*,*) 'reset Tsgroup to 1'
 7788 continue
      if (slopeflag.ne.0) then
         TsGroupNode(:,:)=0
         do i=1,ne
            do j=1,3
               TsGroupNode(nm(i,j),TsGroup(i))=TsGroup(i)
c               if (myproc.eq.1.and.nm(i,j).eq.193) then
c                  write(911,*) i,Tsgroup(i)
c               endif
            enddo
         enddo
      endif
c
c
#ifdef CMPI
c
      call update_msgpass_elem
c
      if (slopeflag.ne.0) then
         do kk=1,Ndom
            do i=1,np
               TsGroupNodetemp(i)=TsGroupNode(i,kk)
            enddo
            call UPDATEI(TsGroupNodetemp,TsGroupNodetemp,1)
            do i=1,np
               TsGroupNode(i,kk)=TsGroupNodetemp(i)
            enddo
         enddo
c
         call update_msgpass_node
      endif

#endif               
c
      do i=1,ne
         do j=1,Ndom
            if (tsgroup(i).eq.j) then
               nelem(j)=nelem(j)+1
               elemnum(j,nelem(j))=i
            endif
         enddo
      enddo
c      if (myproc.eq.0) write(801,*) time_a/86400.
c      do j=1,Ndom
c         idtemp=nelem(j)
c         call ErrorElevSum(idtemp)
c         if (myproc.eq.0) then
c            write(801,*) j,idtemp
c         endif
c      enddo
c      OPEN(69,FILE=DIRNAME//'/'//'Hot_start.299')
c      write(69,*) Ndom,(Dtdom(i),i=1,Ndom)
c      do i=1,ne
c         write(69,*) TsGroup(i)
c      enddo
c      close(69)

c
      return
      end subroutine
