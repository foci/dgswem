compiler=intel-lonestar
#compiler=nvhpc
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
  CC            :=  icx
  FFLAGS1       :=  -r8 $(INCDIRS) -O2 -msse3 -132  #-traceback -check all #-prof-gen -prof-dir /bevo2/michoski/v21/work -pg -prof-use
  FFLAGS2       :=  $(FFLAGS1)
  FFLAGS3       :=  $(FFLAGS1)
  FFLAGS4       :=  $(FFLAGS1)
  DA            :=  -DREAL8 -DLINUX -DCSCA -DRKSSP -DSLOPE5 #-DDGOUT #-DOUT_TEC #-DARTDIF  #  -DWETDR # -DSED_LAY -DOUT_TEC #-DRKC -DTRACE -DSED_LAY -DCHEM -DP0 -DP_AD -DSLOPEALL -DFILTER
  DP            :=  -DREAL8 -DLINUX -DCSCA -DCMPI -DRKSSP -DSLOPE5 #-DDGOUT #-DOUT_TEC #-DARTDIF #-DWETDR  # -DSED_LAY -DRKC -DTRACE -DSED_LAY -DCHEM -DP0 -DP_AD -DSLOPEALL
  DPRE          :=  -DREAL8 -DLINUX -DRKSSP -DSLOPE5 #-DOUT_TEC #-DSWAN #-DARTDIF  # -DWETDR #-DSED_LAY -DSWAN #-DOUT_TEC #-DSWAN -DRKC -DTRACE -DSED_LAY -DCHEM -DP0 -DP_AD -DSLOPEALL
  DPRE2         :=  -DREAL8 -DLINUX -DCMPI 
  CFLAGS        :=  -O3 -xSSSE3 -I. -Wno-implicit-function-declaration
  IMODS         :=  -I
  LIBS          :=  -L ../metis -lmetis
  MSGLIBS       :=
endif

ifeq ($(compiler),nvhpc)   # NVIDIA
  sz            := 8
  PPFC	        :=  nvfortran
  FC	        :=  nvfortran
  PFC	        :=  mpif90
  FFLAGS1	:=  -r$(sz) -Mextend -traceback -g -O3
  FFLAGS2	:=  $(FFLAGS1)
  FFLAGS3	:=  $(FFLAGS1)
  FFLAGS4	:=  $(FFLAGS1)
  DA  	        :=  -DREAL$(sz) -DLINUX -DCSCA -DRKSSP -DSLOPE5
  DP  	        :=  -DREAL$(sz) -DLINUX -DCSCA -DCMPI -DRKSSP -DSLOPE5
  DPRE	        :=  -DREAL$(sz) -DLINUX -DRKSSP -DSLOPE5
  DPRE2         :=  -DREAL$(sz) -DLINUX -DCMPI
  IMODS 	:=  -I
  CC            :=  nvc
  CXX    := nvc++
  CXXFLAGS := -O3 -g
  CFLAGS        :=  -O3 -I. -fastsse -DLINUX
  LIBS  	:=  -L ../metis  -lmetis
  MSGLIBS	:=
endif
