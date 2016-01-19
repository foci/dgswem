SHELL = /bin/sh

VPATH := ../src/
ODIR = odir/


#FC = ifort
FC = gfortran
#FFLAGS = -O3 -132 -I$(ODIR) -Dadcirc
FFLAGS = -extend-source -I$(ODIR) -Dadcirc
# FFLAGS = -132 -C -g -traceback -I$(ODIR) #-Dadcirc

FULLNAME := $(shell hostname)
NAME := $(strip $(subst ., , $(FULLNAME)))
#FC_VERSION := $(shell $(FC) -v 2>&1)
FC_VERSION := gfortran4.9.3
DATE := $(shell date)
GIT_SHA := $(shell git rev-parse HEAD)
GIT_BRANCH := $(shell git rev-parse --abbrev-ref HEAD )
GIT_MOD := $(shell git diff-index --quiet HEAD; echo $$?)
GIT_MOD := $(strip $(subst 1,+,$(GIT_MOD)))
GIT_MOD := $(strip $(subst 0,,$(GIT_MOD)))
GIT_MFILES := $(strip $(shell git ls-files ../* -m))

########################################################################
#                  Library Links (Platform Specific)                   #
########################################################################

ifeq ($(NAME),chl-tilos)
  LAPACK_LIB = -llapack
endif

ifeq ($(NAME),sbrus-laptop)
  LAPACK_LIB = -llapack
endif

ifeq ($(NAME),zachtop)
  LAPACK_LIB = -llapack
endif

ifeq ($(filter crc nd edu,$(NAME)),crc nd edu)
 LAPACK_LIB = -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a  $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm
endif

ifeq ($(filter tacc utexas edu,$(NAME)),tacc utexas edu)
  LAPACK_LIB =  -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm
endif

########################################################################
#                               Objects                                #
########################################################################

objects = globals.o version.o allocation.o basis.o evaluate.o read_input.o read_grid.o read_solution.o kdtree2.o find_nesting.o area_qpts.o error.o

obj := $(patsubst %.o, $(ODIR)%.o,$(objects))

########################################################################
#                          Executable Target                           #
########################################################################
.PHONY :  version $(ODIR) error 

error:  version $(ODIR) $(obj)
	@echo "\n\n"
	@echo "Git SHA: $(GIT_SHA)"
	@echo "Modified: $(GIT_MOD) \n"
	$(FC) $(FFLAGS) -o error $(obj) $(LAPACK_LIB)

########################################################################
#                        Premilinary Targets                           #
########################################################################	

$(ODIR) :
	mkdir -p $@

#version:
#	@sed -i '10 c \      gitBranch = "$(GIT_BRANCH)" ' $(VPATH)version.F90
#	@sed -i '11 c \      gitSHA = "$(GIT_SHA) $(GIT_MOD)" ' $(VPATH)version.F90
	@sed -i '12 c \      compiler_version = "$(FC_VERSION)" ' $(VPATH)version.F90
#	@sed -i '13 c \      compiler_flags = "$(FFLAGS)" ' $(VPATH)version.F90
#	@sed -i '14 c \      modified_files = "$(GIT_MFILES)" ' $(VPATH)version.F90
#	@sed -i '15 c \      compile_date = "$(DATE)" ' $(VPATH)version.F90
#	@sed -i '16 c \      host = "$(NAME)" ' $(VPATH)version.F90	

########################################################################
#                            File Targets                              #
########################################################################	

$(ODIR)globals.o: globals.f90
	$(FC) $(FFLAGS) -c $< -o $@
	mv globals.mod $(ODIR)

$(ODIR)allocation.o: allocation.f90                             $(ODIR)globals.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv allocation.mod $(ODIR)

$(ODIR)basis.o: basis.f90                                       $(ODIR)globals.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv basis.mod $(ODIR)

$(ODIR)evaluate.o: evaluate.F90                                 $(ODIR)globals.o  $(ODIR)basis.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv evaluate.mod $(ODIR)

$(ODIR)read_input.o: read_input.f90                             $(ODIR)globals.o
	$(FC) $(FFLAGS) -c $< -o $@

$(ODIR)read_grid.o: read_grid.f90                               $(ODIR)globals.o  $(ODIR)allocation.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv read_grid.mod $(ODIR)

$(ODIR)read_solution.o: read_solution.F90                       $(ODIR)globals.o
	$(FC) $(FFLAGS) -c $< -o $@

$(ODIR)kdtree2.o: kdtree2.F
	$(FC) $(FFLAGS) -c $< -o $@
	mv kdtree2_*.mod  $(ODIR)

$(ODIR)find_nesting.o: find_nesting.f90                         $(ODIR)globals.o  $(ODIR)kdtree2.o  $(ODIR)evaluate.o  
	$(FC) $(FFLAGS) -c $< -o $@

$(ODIR)area_qpts.o: area_qpts.f90
	$(FC) $(FFLAGS) -c $< -o $@

$(ODIR)version.o: version.F90
	$(FC) $(FFLAGS) -c $< -o $@	

$(ODIR)error.o: error.F90                                       $(ODIR)globals.o  $(ODIR)allocation.o  $(ODIR)basis.o  $(ODIR)read_grid.o  $(ODIR)evaluate.o  
	$(FC) $(FFLAGS) -c $< -o $@

########################################################################
#                            Clean Targets                             #
########################################################################

.PHONY : clean
clean : 
	rm -r error $(ODIR)
