########################################################################
# Intel-Linux computers using Portland group compilers

########################################################################
# Linux computers using Intel or PGI compiler
# These flags were received from Tim Campbell
#
# ***NOTE*** User must select between various Linux setups
#            by commenting/uncommenting these sections
#
ifeq ($(MACHINE)-$(OS),i686-linux-gnu)
compiler=intel
#compiler=pgi
ifeq ($(compiler),intel)   # Intel
  PPFC	        :=  ifort
  FC	        :=  ifort
  PFC	        :=  mpif77
#  FFLAGS1	:=  -w -O3 -axP -FI -132
#  FFLAGS2	:=  -w -O3 -axP -FI -132
#  FFLAGS3	:=  -w -O3 -axP -FI -132
#  FFLAGS4	:=  -w -O3 -axP -FI -132
  FFLAGS1	:=  -w -O2 -FI -Vaxlib
  FFLAGS2	:=  -w -O2 -FI -Vaxlib
  FFLAGS3	:=  -w -O2 -FI -Vaxlib
  FFLAGS4	:=  -w -O2 -FI -Vaxlib
  DA  	        :=  -DREAL8 -DLINUX -DCSCA
  DP  	        :=  -DREAL8 -DLINUX -DCSCA -DCMPI
  DPRE	        :=  -DREAL8 -DLINUX
  DPRE2	        :=  -DREAL8 -DLINUX -DCMPI
  IMODS		:=  -I
  CC            :=  icc
  CFLAGS        :=  -I. -O2 -DLINUX
  CLIBS         := 
  LIBS  	:=  -L ../metis  -lmetis
  MSGLIBS	:=  
endif
ifeq ($(compiler),pgi)   # PGI
  PPFC	        :=  pgf90
  FC	        :=  pgf90
  FC	        :=  ifort
  PFC	        :=  mpif90
  FFLAGS1	:=  -fastsse
  FFLAGS2	:=  -fastsse
  FFLAGS3	:=  -fastsse
  FFLAGS4	:=  -fastsse
  DA  	        :=  -DREAL8 -DLINUX -DCSCA
  DP  	        :=  -DREAL8 -DLINUX -DCSCA -DCMPI
  DPRE	        :=  -DREAL8 -DLINUX
  DPRE2         :=  -DREAL8 -DLINUX -DCMPI
  IMODS 	:=  -I
  CC            :=  pgcc
  CFLAGS        :=  -I. -fastsse -DLINUX
  CLIBS         := 
  LIBS  	:=  -L ../metis  -lmetis
  MSGLIBS	:=  
endif
endif

########################################################################
# Apple/IBM 
# xlf/xlc compilers or absoft compilers
# powerpc-apple-darwin8.*.0

ifneq (,$(findstring powerpc-darwin,$(MACHINE)-$(OS)))
#compiler=ibmxl
compiler=absoft
#compiler=gfortran
ifeq ($(compiler),ibmxl)   # IBM XL compiler
  PPFC	        := mpif77
  FC	        := mpif77
  PFC	        := mpif77
  FFLAGS1	:=  -w -O2 -w -qstrict -qfixed -qarch=auto -qcache=auto -I .
  FFLAGS2	:=  -w -O2 -w -qstrict -qfixed -qarch=auto -qcache=auto -I .
  FFLAGS3	:=  -w -O2 -w -qstrict -qfixed -qarch=auto -qcache=auto -I .
  FFLAGS4	:=  -w -O2 -w -qstrict -qfixed -qarch=auto -qcache=auto -I .
  DA  	   	:=  -WF,"-DREAL8 -DCSCA -DIBM -DIBMXL"
  DP  	   	:=  -WF,"-DREAL8 -DCSCA -DCMPI -DIBM -DIBMXL"
  DPRE	   	:=  -WF,"-DREAL8 -DIBM -DIBMXL"
  DPRE2	   	:=  -WF,"-DREAL8 -DCMPI -DLINUX -DIBM -DIBMXL"
  IMODS  	:=  -I 
  PMODS  	:=  -p 
  CC            :=  mpicc
# CFLAGS	:=  -I. -O2 -bmaxdata:0x30000000 -DIBMXL
# LDFLAGS	:=  -bmaxdata:0x30000000
  CFLAGS	:=  -I. -O2 -w -qstrict -qfixed -qarch=auto -qcache=auto
  LDFLAGS	:=  
# LIBS	        :=  -L ../metis  -lmetis -L /Applications/Absoft/lib -lU77
  LIBS	        :=  -L ../metis  -lmetis 
  MSGLIBS	:=  -lm
endif

ifeq ($(compiler),absoft)   # absoft
  PPFC          := f90
  FC            := f90
  PFC           := mpif77
  FFLAGS1       :=  -w -O3 -cpu:g5 -f fixed -DIBM -I .
  FFLAGS2       :=  -w -O3 -cpu:g5 -f fixed -I . -DIBM
  FFLAGS3       :=  -w -O3 -cpu:g5 -f fixed -I . -N11
  FFLAGS4       :=  -w -O3 -cpu:g5 -f fixed -I . -N11
  DA            :=  -DREAL8 -DCSCA -DIBM -W132
  DP            :=  -DREAL8 -DCSCA -DIBM -W132 -DCMPI -DBLKOUT
#  DP            :=  -DREAL8 -DCSCA -DIBM -W132 -DCMPI -DBLKOUT
  DPRE          :=  -DREAL4 -DIBM
  DPRE2         :=  -DREAL4 -DIBM
  IMODS         :=  -p
  CC            :=  gcc
  CFLAGS        :=  -I. -O2 -DIBM
  LDFLAGS       :=
  LIBS          :=  -L ../metis  -lmetis -lU77
  PERFLIBS      :=
# PERFLIBS      :=  -L$(HOME)/bin -lparaperf  -lezstub 
  MSGLIBS       :=  -lm
endif

ifeq ($(compiler),gfortran)   # gfortran
  PPFC          := gfortran
  FC            := gfortran
  PFC           := not-specified
  FFLAGS1       :=  -O2 -DLINUX
  FFLAGS2       :=  -O2 -DLINUX
  FFLAGS3       :=  -O2 -DLINUX
  FFLAGS4       :=  -O2 -DLINUX
  DA            :=  -DREAL8 -DCSCA -DLINUX
  DP            :=  -DREAL8 -DCSCA -DLINUX -DCMPI -DBLKOUT
#  DP            :=  -DREAL8 -DCSCA -DIBM -W132 -DCMPI -DBLKOUT
  DPRE          :=  -DREAL4 -DLINUX
  DPRE2         :=  -DREAL4 -DLINUX
  IMODS         :=  -I 
  CC            :=  mpicc
  CFLAGS        :=  -I. -O2 -DIBM
  LDFLAGS       :=
  LIBS          :=  -L ../metis  -lmetis -lU77
  PERFLIBS      :=
# PERFLIBS      :=  -L$(HOME)/bin -lparaperf  -lezstub 
  MSGLIBS       :=  -lm
endif

endif

########################################################################
# Windows computers using gfortran or PGI compiler on Cygwin
#
# ***NOTE*** User must select between the provided compiler options
#
ifeq ($(MACHINE)-$(OS),i686-pc-cygwin)
#compiler=gfortran
compiler=pgi
ifeq ($(compiler),gfortran)   # gfortran
  PPFC	        :=  gfc
  FC	        :=  gfc
  PFC	        :=  mpif77
  FFLAGS1	:=  -w -O3 -mtune=opteron
  FFLAGS2	:=  -w -O3 -mtune=opteron
  FFLAGS3	:=  -w -O3 -mtune=opteron
  FFLAGS4	:=  -w -O3 -mtune=opteron
  DA  	        :=  -DREAL8 -DLINUX -DGFORTRAN -DCSCA
  DP  	        :=  -DREAL8 -DLINUX -DGFORTRAN -DCSCA -DCMPI
  DPRE	        :=  -DREAL8 -DLINUX -DGFORTRAN
  DPRE2	        :=  -DREAL8 -DLINUX -DGFORTRAN -DCMPI
  IMODS		:=  -I
  CC            :=  gcc
  CFLAGS        :=  -I. -O2 -DLINUX
  CLIBS         := 
  LIBS  	:=  -L ../metis  -lmetis
  MSGLIBS	:=  
endif
ifeq ($(compiler),pgi)   # pgi
  PPFC	        :=  pgf90
  FC	        :=  pgf90
  PFC	        :=  mpif77
  FFLAGS1	:=  -O2 -tp athlonxp
  FFLAGS2	:=  -O2 -tp athlonxp
  FFLAGS3	:=  -O2 -tp athlonxp
  FFLAGS4	:=  -O2 -tp athlonxp
#  FFLAGS1	:=  
# FFLAGS2	:= 
# FFLAGS3	:=
# FFLAGS4	:=
  DA  	        :=  -DREAL8 -DLINUX -DCSCA
  DP  	        :=  -DREAL8 -DLINUX -DCSCA -DCMPI
  DPRE	        :=  -DREAL8 -DLINUX
  DPRE2	        :=  -DREAL8 -DLINUX -DCMPI
  IMODS		:=  -I
  CC            :=  pgcc
  CFLAGS        :=  -I. -O2 -DLINUX
  CLIBS         := 
  LIBS  	:=  -L ../metis  -lmetis
  MSGLIBS	:=  
endif
endif

########################################################################
# SGI Origin - SGI using standard compilers

ifneq (,$(findstring mips-irix,$(MACHINE)-$(OS)))
  PPFC              := f90 
  FC              := f90 
  PFC             := f90 
  FFLAGS1	  := -O2 -OPT:Olimit=4257 -OPT:reorg_common=ON
  FFLAGS2	  := -O2 -OPT:Olimit=4257 -OPT:reorg_common=ON
  FFLAGS3	  := -O2 -OPT:Olimit=4257 -OPT:reorg_common=ON
  FFLAGS4	  := -O2 -OPT:Olimit=4257 -OPT:reorg_common=ON
  DA	          :=  -DREAL8 -DCSCA 
  DP	          :=  -DREAL8 -DCMPI  -DSGI -DCSCA 
  DPRE	          :=  -DREAL4 -DSGI
  DPRE2	          :=  -DREAL4 -DCMPI -DSGI
  IMODS   	  :=  -I
  CC              :=  cc   
  CFLAGS          :=  -I. -O2
  LIBS		  :=  -L ../metis  -lmetis
  MSGLIBS	  := -lmpi
endif

########################################################################
# IBM SP - IBM using standard compilers

ifneq (,$(findstring powerpc-aix,$(MACHINE)-$(OS)))
  PPFC          := xlf90_r -q64
  FC            := xlf90_r -q64
  PFC           := mpxlf90_r -q64
  CC            := mpcc_r  -q64
  FFLAGS1       := -w -O2 -qfixed -qarch=auto -qcache=auto
  FFLAGS2       := -w -O2 -qfixed -qarch=auto -qcache=auto
  FFLAGS3       := -w -O2 -qhot -qstrict -qfixed=132 -qarch=auto -qcache=auto
  FFLAGS4       := -w -O3 -qstrict -qfixed=132 -qarch=auto -qcache=auto -qinitauto -qipa=inline=roe_flux:edge_int_hydro
  DA            := -WF,"-DREAL8,-DIBM,-DCSCA"
  DP            := -tF -WF,"-DREAL8,-DIBM,-DCSCA,-DCMPI"
# DP            := -tF -WF,"-DREAL8,-DIBM,-DCSCA,-DCMPI,-DUSE_PARAPERF"
  DPRE          := -tF -WF,"-DREAL4,-DIBM"
  DPRE2         := -tF -WF,"-DREAL4,-DIBM,-DCMPI"
  IMODS         := -I
  CFLAGS        := -I. -O2 -DIBM
  ARFLAGS       := -X64
  LIBS          := -L ../metis  -lmetis
  PERFLIBS      :=
# PERFLIBS      :=  -L ../metis -lmetis -L$(HOME)/bin -lparaperf -L$(LOCAL_LIBRARY) -lpapi64 -lpmapi
  MSGLIBS       := -lm
endif

########################################################################
# Sun Solaris Using standard compilers

ifneq (,$(findstring sparc-solaris,$(MACHINE)-$(OS)))
  PPFC	        := tmf90
  FC	        := tmf90
  PFC	        := tmf90
  ARCH	        := v8plusa
  LD	        := $(FC)
  FFLAGS1	:= -fast -xO4 -xdepend -fsimple=1 -f -dalign -xtarget=ultra -xarch=$(ARCH)
  FFLAGS2	:= -fast -xO4 -xdepend -fsimple=1 -f -dalign -xtarget=ultra -xarch=$(ARCH)
  FFLAGS3	:= -fast -xO4 -xdepend -fsimple=1 -f -dalign -xtarget=ultra -xarch=$(ARCH)
  FFLAGS4	:= -fast -xO4 -xdepend -fsimple=1 -f -dalign -xtarget=ultra -xarch=$(ARCH)
  DA	        := -DREAL8 -DMACHSUN -DCSCA 
  DP	        := -DREAL8 -DCMPI -DMACHSUN -DCSCA 
  DPRE          := -DREAL4 -DMACHSUN
  DPRE2         := -DREAL4 -DCMPI -DMACHSUN
  IMODS		:= -M
  CC       	:= tmcc 
  CFLAGS   	:= -I. 
  LIBS     	:= -L../metis  -lmetis
  MSGLIBS  	:= -lmpi
endif
########################################################################
# Intel-Linux computers using Portland group compilers

ifneq (,$(findstring alphaev6-linux,$(MACHINE)-$(OS)))
  PPFC	        :=  fort
  FC	        :=  fort
  PFC	        :=  mpif90
  FFLAGS1	:=  -fast -fixed 
  FFLAGS2	:=  -fast -fixed
  FFLAGS3	:=  -fast -fixed 
  FFLAGS4	:=  -fast -fixed 
  DA  	        :=  -DREAL8 -DCSCA   
  DP  	        :=  -DREAL8 -DCMPI  -DCSCA 
  DPRE	        :=  -DREAL4 -DLINUX 
  DPRE2	        :=  -DREAL4 -DCMPI -DLINUX
  IMODS 	:=  -I
  CC            := ccc
  CFLAGS        := -I. -DLINUX -O2
  CLIBS         := 
  LIBS  	:=  -L ../metis  -lmetis
  MSGLIBS	:=  
endif

########################################################################
# Compiler flags for Linux operating system on 64bit x86 CPU
#
ifneq (,$(findstring x86_64-linux-gnu,$(MACHINE)-$(OS)))
#
# ***NOTE*** User must select between various Linux setups
#            by commenting/uncommenting the appropriate compiler
#
compiler=intel-lonestar
#compiler=cray_xt3
#compiler=kraken
#
#
ifeq ($(compiler),cray_xt3)
  PPFC          :=  pgf90
  FC            :=  ftn
  PFC           :=  ftn
  CC            :=  pgcc
  FFLAGS1       :=  -O2
  FFLAGS2       :=  -O2
  FFLAGS3       :=  -g -O2 -Mextend
  FFLAGS4       :=  -tp k8-64 -fastsse -Mextend -Minline=name:roe_flux,name:edge_int_hydro
  DA            :=  -DREAL8 -DLINUX -DCSCA
  DP            :=  -DREAL8 -DCMPI  -DLINUX -DCSCA
# DP            :=  -DREAL8 -DCMPI  -DLINUX -DCSCA -DUSE_PARAPERF
  DPRE          :=  -DREAL4 -DLINUX
  DPRE2         :=  -DREAL4 -DLINUX -DCMPI
  CFLAGS        :=  -I. -DLINUX
  IMODS         :=  -module 
  LIBS          :=  -L ../metis -lmetis
  PERFLIBS      :=
# PERFLIBS      :=  -L$(HOME)/bin -lparaperf -L/opt/xt-tools/papi/3.2.1/lib/cnos64 -lpapi -lperfctr
  MSGLIBS       :=
endif
#
# sb46.50.02 These flags work on the UT Austin Lonstar cluster.
ifeq ($(compiler),intel-lonestar)
  PPFC          :=  ifort	
  FC            :=  ifort
  PFC           :=  mpif90
  CC            :=  icc
  FFLAGS1       :=  $(INCDIRS) -O0 -msse3 -132  #-debug -traceback -g -check all #-prof-gen -prof-dir /bevo2/michoski/v21/work -pg -prof-use
  FFLAGS2       :=  $(FFLAGS1)
  FFLAGS3       :=  $(FFLAGS1)
  FFLAGS4       :=  $(FFLAGS1)
  DA            :=  -DREAL8 -DLINUX -DCSCA -DRKSSP -DWETDR -DOUT_TEC -DSLOPE5 #-DARTDIF  #  -DWETDR -DSLOPE5 # -DSED_LAY -DOUT_TEC #-DRKC -DTRACE -DSED_LAY -DCHEM -DP0 -DP_AD -DSLOPEALL -DFILTER
  DP            :=  -DREAL8 -DLINUX -DCSCA -DCMPI -DRKSSP -DWETDR -DOUT_TEC -DSLOPE5 #-DARTDIF #-DWETDR -DSLOPE5 # -DSED_LAY -DRKC -DTRACE -DSED_LAY -DCHEM -DP0 -DP_AD -DSLOPEALL
  DPRE          :=  -DREAL8 -DLINUX -DRKSSP -DOUT_TEC -DWETDR -DSLOPE5 #-DSWAN #-DARTDIF  # -DWETDR -DSLOPE5 #-DSED_LAY -DSWAN #-DOUT_TEC #-DSWAN #-DSLOPE5 -DRKC -DTRACE -DSED_LAY -DCHEM -DP0 -DP_AD -DSLOPEALL
  DPRE2         :=  -DREAL8 -DLINUX -DCMPI 
  CFLAGS        :=  -O3 -xSSSE3 -I.
  IMODS         :=  -I
  LIBS          :=  -L ../metis -lmetis
  MSGLIBS       :=
endif

endif
