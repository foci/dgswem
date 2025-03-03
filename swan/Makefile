# ------------------------------------------------------------------------------
#                      Makefile for building SWAN program and documentation
# ------------------------------------------------------------------------------
#
# Before compilation, type "make config" first!
#
# To compile the serial executable type "make ser"
# To compile the OpenMP executable type "make omp"
# To compile the MPI executable type "make mpi"
# To compile the PunSWAN executable "make punswan"
#
# To remove compiled objects and modules: type "make clean"
#
# To remove compiled objects, modules and executable: type "make clobber"
#
# To compile the SWAN documentation type "make doc"
#
# To remove the PDF and HTML documents type "make cleandoc"
#
# Please do not change anything below, unless you're very sure what you're doing
# ------------------------------------------------------------------------------

include macros.inc

SWAN_EXE = swan.exe

SWAN_OBJS = \
swmod1.$(EXTO) \
swmod2.$(EXTO) \
m_constants.$(EXTO) \
m_fileio.$(EXTO) \
serv_xnl4v5.$(EXTO) \
mod_xnl4v5.$(EXTO) \
SwanGriddata.$(EXTO) \
SwanGridobjects.$(EXTO) \
SwanCompdata.$(EXTO) \
swanmain.$(EXTO) \
swanpre1.$(EXTO) \
swanpre2.$(EXTO) \
swancom1.$(EXTO) \
swancom2.$(EXTO) \
swancom3.$(EXTO) \
swancom4.$(EXTO) \
swancom5.$(EXTO) \
swanout1.$(EXTO) \
swanout2.$(EXTO) \
swanser.$(EXTO) \
swanparll.$(EXTO) \
SwanReadGrid.$(EXTO) \
SwanReadADCGrid.$(EXTO) \
SwanReadTriangleGrid.$(EXTO) \
SwanReadEasymeshGrid.$(EXTO) \
SwanInitCompGrid.$(EXTO) \
SwanCheckGrid.$(EXTO) \
SwanCreateEdges.$(EXTO) \
SwanGridTopology.$(EXTO) \
SwanGridVert.$(EXTO) \
SwanGridCell.$(EXTO) \
SwanGridFace.$(EXTO) \
SwanPrintGridInfo.$(EXTO) \
SwanFindPoint.$(EXTO) \
SwanPointinMesh.$(EXTO) \
SwanBpntlist.$(EXTO) \
SwanPrepComp.$(EXTO) \
SwanVertlist.$(EXTO) \
SwanCompUnstruc.$(EXTO) \
SwanDispParm.$(EXTO) \
SwanPropvelX.$(EXTO) \
SwanSweepSel.$(EXTO) \
SwanPropvelS.$(EXTO) \
SwanTranspAc.$(EXTO) \
SwanTranspX.$(EXTO) \
SwanDiffPar.$(EXTO) \
SwanGSECorr.$(EXTO) \
SwanInterpolatePoint.$(EXTO) \
SwanInterpolateAc.$(EXTO) \
SwanInterpolateOutput.$(EXTO) \
SwanConvAccur.$(EXTO) \
SwanConvStopc.$(EXTO) \
SwanFindObstacles.$(EXTO) \
SwanCrossObstacle.$(EXTO) \
SwanComputeForce.$(EXTO) \
SwanIntgratSpc.$(EXTO) \
SwanReadfort18.$(EXTO) \
SwanSumOverNodes.$(EXTO) \
SwanMinOverNodes.$(EXTO) \
SwanMaxOverNodes.$(EXTO) \
ocpids.$(EXTO) \
ocpcre.$(EXTO) \
ocpmix.$(EXTO)

HCAT_EXE = hcat.exe
HCAT_OBJS = swanhcat.$(EXTO)

MSG_OBJS = \
$(O_DIR)sizes.$(EXTO) \
$(O_DIR)global.$(EXTO) \
$(O_DIR)messenger.$(EXTO) \
$(O_DIR)dg.$(EXTO) \
$(O_DIR)diff45_41.$(EXTO)

.SUFFIXES: .f .F .for .f90 .F90

.PHONEY: help

help:
	@echo "This Makefile supports the following:"
	@echo "make config    -- makes machine-dependent macros include file"
	@echo "make ser       -- makes the serial $(SWAN_EXE) executable"
	@echo "make omp       -- makes the OpenMP $(SWAN_EXE) executable"
	@echo "make mpi       -- makes the    MPI $(SWAN_EXE) executable"
	@echo "make punswan   -- makes the parallel un$(SWAN_EXE) executable"
	@echo "make doc       -- makes the SWAN documentation (PDF)"
	@echo "make clean     -- removes compiled objects and modules"
	@echo "make clobber   -- removes compiled objects, modules and $(SWAN_EXE)"
	@echo "make cleandoc  -- removes all SWAN documents"

config:
	@perl platform.pl

install:
	@perl platform.pl

ser:
#Casey 080108: Added the -swan flag to the following line.
	@perl switch.pl $(swch) -swan *.ftn *.ftn90
	$(MAKE) FOR=$(F90_SER) FFLAGS="$(FLAGS_OPT) $(FLAGS_MSC) $(FLAGS_SER)" \
                FFLAGS90="$(FLAGS_OPT) $(FLAGS_MSC) $(FLAGS_SER) $(FLAGS_F90)" \
                INCS="$(INCS_SER)" LIBS="$(LIBS_SER)" OBJS="$(SWAN_OBJS)" $(SWAN_EXE)
	$(MAKE) FOR="$(F90_SER)" FFLAGS="$(FLAGS_OPT) $(FLAGS_MSC) $(FLAGS_SER)" \
                FFLAGS90="$(FLAGS_OPT) $(FLAGS_MSC) $(FLAGS_SER) $(FLAGS_F90)" $(HCAT_EXE)

omp:
#Casey 090826: Added the -swan flag to the following line.
	@perl switch.pl $(swch) -omp -swan *.ftn *.ftn90
	$(MAKE) FOR=$(F90_OMP) FFLAGS="$(FLAGS_OPT) $(FLAGS_MSC) $(FLAGS_OMP)" \
                FFLAGS90="$(FLAGS_OPT) $(FLAGS_MSC) $(FLAGS_OMP) $(FLAGS_F90)" \
                INCS="$(INCS_OMP)" LIBS="$(LIBS_OMP)" OBJS="$(SWAN_OBJS)" $(SWAN_EXE)
	$(MAKE) FOR="$(F90_SER)" FFLAGS="$(FLAGS_OPT) $(FLAGS_MSC) $(FLAGS_SER)" \
                FFLAGS90="$(FLAGS_OPT) $(FLAGS_MSC) $(FLAGS_SER) $(FLAGS_F90)" $(HCAT_EXE)

mpi:
#Casey 090826: Added the -swan flag to the following line.
	@perl switch.pl $(swch) -mpi -swan *.ftn *.ftn90
	$(MAKE) FOR=$(F90_MPI) FFLAGS="$(FLAGS_OPT) $(FLAGS_MSC) $(FLAGS_MPI)" \
                FFLAGS90="$(FLAGS_OPT) $(FLAGS_MSC) $(FLAGS_MPI) $(FLAGS_F90)" \
                INCS="$(INCS_MPI)" LIBS="$(LIBS_MPI)" OBJS="$(SWAN_OBJS)" $(SWAN_EXE)
	$(MAKE) FOR="$(F90_SER)" FFLAGS="$(FLAGS_OPT) $(FLAGS_MSC) $(FLAGS_SER)" \
                FFLAGS90="$(FLAGS_OPT) $(FLAGS_MSC) $(FLAGS_SER) $(FLAGS_F90)" $(HCAT_EXE)

punswan:
#Casey 080118: Added the -swan flag to the following line.
	@perl switch.pl $(swch) -pun -swan *.ftn *.ftn90
	$(MAKE) FOR=$(F90_MPI) FFLAGS="$(FLAGS_OPT) $(FLAGS_MSC) $(FLAGS_MPI)" \
                FFLAGS90="$(FLAGS_OPT) $(FLAGS_MSC) $(FLAGS_MPI) $(FLAGS_F90)" \
                INCS="$(INCS_MPI) -I$(O_DIR)" LIBS="$(LIBS_MPI)" \
                OBJS="$(MSG_OBJS) $(SWAN_OBJS)" $(SWAN_EXE)
	$(MAKE) FOR="$(F90_SER)" FFLAGS="$(FLAGS_OPT) $(FLAGS_MSC) $(FLAGS_SER)" \
                FFLAGS90="$(FLAGS_OPT) $(FLAGS_MSC) $(FLAGS_SER) $(FLAGS_F90)" $(HCAT_EXE)

doc:
	$(MAKE) -f Makefile.latex TARGET=swanuse doc
	$(MAKE) -f Makefile.latex TARGET=swantech doc
	$(MAKE) -f Makefile.latex TARGET=swanimp doc
	$(MAKE) -f Makefile.latex TARGET=swanpgr doc
	$(MAKE) -f Makefile.latex TARGET=latexfordummies doc

$(HCAT_EXE): $(HCAT_OBJS)
	$(FOR) $(HCAT_OBJS) $(OUT)$(HCAT_EXE)

$(SWAN_EXE): $(SWAN_OBJS)
	$(FOR) $(OBJS) $(OUT)$(SWAN_EXE) $(INCS) $(LIBS)

.f.o:
	$(FOR) $< -c $(FFLAGS) $(INCS)

.f90.o:
	$(FOR) $< -c $(FFLAGS90) $(INCS)

.F.o:
	$(FOR) $< -c $(FFLAGS) $(INCS)

.F90.o:
	$(FOR) $< -c $(FFLAGS90) $(INCS)

.for.o:
	$(FOR) $< -c $(FFLAGS) $(INCS)

.for.obj:
	$(FOR) $< -c $(FFLAGS) $(INCS)

.f90.obj:
	$(FOR) $< -c $(FFLAGS90) $(INCS)

clean:
	$(RM) *.$(EXTO) *.mod

clobber:
	$(RM) *.$(EXTO) *.mod *.f *.F *.F90 *.for *.f90 $(SWAN_EXE) $(HCAT_EXE)

allclean:
	$(RM) *.$(EXTO) *.mod *.f *.F *.F90 *.for *.f90 $(SWAN_EXE) $(HCAT_EXE)

cleandoc:
	$(MAKE) -f Makefile.latex TARGET=swanuse cleandoc
	$(MAKE) -f Makefile.latex TARGET=swantech cleandoc
	$(MAKE) -f Makefile.latex TARGET=swanimp cleandoc
	$(MAKE) -f Makefile.latex TARGET=swanpgr cleandoc
	$(MAKE) -f Makefile.latex TARGET=latexfordummies cleandoc
