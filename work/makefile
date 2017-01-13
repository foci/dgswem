SHELL:=/bin/sh
#
#  Makefile to Build DGSWEM and its pre/post-processor 

########################################################################
#  Get Canonical Machine NAME from config.guess
#1
#NAME     := x86_64-unknown-linux-gnu #$(shell ./config.guess)
#LIST     := $(subst -, ,$(NAME))
#MACHINE  := $(word 1, $(LIST))
#VENDOR   := $(word 2, $(LIST))
#OS       := $(subst  $(MACHINE)-$(VENDOR)-,,$(strip $(NAME)))

#include cmplrflags.mk
#Casey 121126: Include variables/flags from SWAN.
#include ../swan/macros.inc


  PPFC          :=  ifort
  FC            :=  ifort
  PFC           :=  mpif90
  CC            :=  icc
  FFLAGS1       :=  $(INCDIRS) -O2 -132  #-traceback -g -check all #-prof-gen -prof-dir /bevo2/michoski/v21/work -pg -prof-use                                                                     
  FFLAGS2       :=  $(FFLAGS1)
  FFLAGS3       :=  $(FFLAGS1)
  FFLAGS4       :=  $(FFLAGS1)
  DA            :=  -DREAL8 -DLINUX -DCSCA -DRKSSP -DSLOPE5 #-DARTDIF  #  -DWETDR -DSLOPE5 # -DSED_LAY -DOUT_TEC #-DRKC -DTRACE -DSED_LAY -DCHEM -DP0 -DP_AD -DSLOPEALL -DFILTER                            
  DP            :=  -DREAL8 -DLINUX -DCSCA -DCMPI -DRKSSP -DSLOPE5 #-DARTDIF #-DWETDR -DSLOPE5 # -DSED_LAY -DRKC -DTRACE -DSED_LAY -DCHEM -DP0 -DP_AD -DSLOPEALL                                            
  DPRE          :=  -DREAL8 -DLINUX -DRKSSP -DSLOPE5 #-DSWAN #-DARTDIF  # -DWETDR -DSLOPE5 #-DSED_LAY -DSWAN #-DOUT_TEC #-DSWAN #-DSLOPE5 -DRKC -DTRACE -DSED_LAY -DCHEM -DP0 -DP_AD -DSLOPEALL   
  DPRE2         :=  -DREAL8 -DLINUX -DCMPI
  CFLAGS        :=  -O3 -xSSSE3 -I.
  IMODS         :=  -I
  LIBS          :=  -L ../metis -lmetis
  MSGLIBS       :=


######################## Target Specific Rules ###################################

#                                      adcprep   
ifeq ($(BUILDTYPE),adcprep)
  CF:= $(PPFC)
  O_DIR:=odir1/
  FFLAGS:= $(FFLAGS1) $(DPRE) $(IMODS)$(O_DIR) 
  VPATH :=  ../prep 
  PMSG_OBJ:=
endif
#                                      adcpost  
ifeq ($(BUILDTYPE),adcpost)
  CF:= $(PPFC)
  O_DIR:=odir3/
  FFLAGS:= $(FFLAGS1) $(DPRE) $(IMODS)$(O_DIR) 
  VPATH :=  ../prep
endif
#                                      adcirc   
ifeq ($(BUILDTYPE),dgswem_serial)
  CF:= $(PPFC)
  O_DIR:=odir4/
  FFLAGS:= $(FFLAGS3) $(DA) $(IMODS)$(O_DIR)
  VPATH:=  ../src 
  MSG_MOBJ:=
endif
#Casey 121126: Add rules for the serial SWAN+ADCIRC.
#                                      adcswan
ifeq ($(BUILDTYPE),swemswan)
  CF:= $(PPFC)
  O_DIR:=odir44/
  FFLAGS:= $(FFLAGS3) -DSWAN $(DA) $(IMODS) $(O_DIR)
  VPATH:=  ../src:../swan
  MSG_MOBJ:= $(O_DIR)globalio.o
endif
#                                      padcirc   
ifeq ($(BUILDTYPE),dgswem)
  CF:= $(PFC)
  O_DIR:=odir5/
  FFLAGS:= $(FFLAGS4) $(DP) $(IMODS)$(O_DIR)
  VPATH :=  ../src 
  MSG_MOBJ:= $(O_DIR)messenger_elem.o $(O_DIR)messenger.o
  MSG_OBJ:= 
endif
#Casey 121126: Add rules for the parallel SWAN+ADCIRC.
#                                      padcswan
ifeq ($(BUILDTYPE),dgswemswan)
  CF:= $(PFC)
  O_DIR:=odir55/
  FFLAGS:= $(FFLAGS4) -DSWAN $(DP) $(IMODS)$(O_DIR)
  VPATH :=  ../src:../swan
  MSG_MOBJ:= $(O_DIR)messenger_elem.o $(O_DIR)messenger.o $(O_DIR)globalio.o
  MSG_OBJ:=
endif

$(O_DIR):
	mkdir -p $@

######################### Module Source, Object, & Mod Files ######################

PREP_MSRC  =  presizes.F 
PGLO_MSRC  =  pre_global.F  
POST_MSRC  =  post_global.F
ADC_MSRC   =  sizes.F diff45_41.F global.F dg.F global_3dvs.F nodalattr.F fparser.F90 fparser2.F90 \
              read_fort_dg.F90 timer.F90              
HARM_MSRC  =  harm.F  
WIND_MSRC  =  wind.F  
#Casey 121126: Add file for coupling to unstructured SWAN.
COUP_MSRC  =  couple2swan.F

PREP_MOBJ:= $(patsubst %.F, $(O_DIR)%.o, $(PREP_MSRC) )
PGLO_MOBJ:= $(patsubst %.F, $(O_DIR)%.o, $(PGLO_MSRC) )
POST_MOBJ:= $(patsubst %.F, $(O_DIR)%.o, $(POST_MSRC) )
ADC_MOBJ := $(patsubst %.F, $(O_DIR)%.o, $(ADC_MSRC)  )
ADC_MOBJ := $(patsubst %.F90, $(O_DIR)%.o, $(ADC_MOBJ)  )
HARM_MOBJ := $(patsubst %.F, $(O_DIR)%.o, $(HARM_MSRC)  )
WIND_MOBJ := $(patsubst %.F, $(O_DIR)%.o, $(WIND_MSRC)  )
#Casey 121126: Add rules for coupling to unstructured SWAN.
COUP_MOBJ := $(patsubst %.F, $(O_DIR)%.o, $(COUP_MSRC)  )


############################# Source & Object Files ##############################

METIS1_SRC =  metis.F
METIS2_SRC =  metis2.F
PREP_SRC   =  adcprep.F  decomp.F   read_global.F   prep.F   interp.F  machdep.F
POST_SRC   =  adcpost.F    post.F  compare.F  diffmerge.F
ADC_SRC    =  dgswem.F owiwind.F read_input.F cstart.F hstart.F radiation_edge_hydro.F \
	ibarrier_edge_hydro.F ebarrier_edge_hydro.F  \
	calc_normal.F create_edge_data.F detector.F \
	DG_timestep.F DG_hydro_timestep.F \
	edge_int_hydro.F flow_edge_hydro.F \
	internal_edge_hydro.F land_edge_hydro.F \
	ocean_edge_hydro.F  \
	orthobasis_area.F orthobasis_edge.F p_enrichment.F prep_DG.F quad_pts_area.F \
	quad_pts_edge.F radiation_edge_hydro.F read_input.F rhs_dg_hydro.F \
	numerical_flux.F tidal_potential.F write_results.F \
	slopelimiter.F prep_slopelim.F ocean_edge_hydro_post.F ocean_flow_edge_hydro.F \
	wetdry.F LDG_hydro.F met_forcing.F sta_basis.F modal2nodal.F write_output.F
#Casey 121126: Added modal2nodal.F and write_output.F to the previous line.

METIS1_OBJ:= $(patsubst %.F, $(O_DIR)%.o, $(METIS1_SRC) )
METIS2_OBJ:= $(patsubst %.F, $(O_DIR)%.o, $(METIS2_SRC) )
PREP_OBJ:= $(patsubst %.F, $(O_DIR)%.o, $(PREP_SRC) )
POST_OBJ:= $(patsubst %.F, $(O_DIR)%.o, $(POST_SRC) )
ADC_OBJ:= $(patsubst %.F, $(O_DIR)%.o, $(ADC_SRC) )

#################### SWAN Source & Object Files #########################
#Casey 121126: Added this section for coupling to unstructured SWAN.

SWAN_MSRC  = swmod1.f swmod2.f m_constants.f90 m_fileio.f90 \
serv_xnl4v5.f90 mod_xnl4v5.f90 \
SwanGriddata.f90 SwanGridobjects.f90 SwanCompdata.f90 \
couple2adcirc.f90

SWAN_MOBJ := $(patsubst %.f,   $(O_DIR)%.o, $(SWAN_MSRC)  )
SWAN_MOBJ := $(patsubst %.f90, $(O_DIR)%.o, $(SWAN_MOBJ)  )

SWAN_SRC   = swanmain.f swanpre1.f swanpre2.f swancom1.f swancom2.f swancom3.f swancom4.f \
swancom5.f swanout1.f swanout2.f swanser.f swanparll.f SwanReadGrid.f90 \
SwanReadADCGrid.f90 SwanReadTriangleGrid.f90 SwanReadEasymeshGrid.f90 \
SwanInitCompGrid.f90 SwanCheckGrid.f90 SwanCreateEdges.f90 SwanGridTopology.f90 SwanGridVert.f90 \
SwanGridCell.f90 SwanGridFace.f90 SwanPrintGridInfo.f90 SwanFindPoint.f90 \
SwanPointinMesh.f90 SwanBpntlist.f90 SwanPrepComp.f90 SwanVertlist.f90 SwanCompUnstruc.f90 \
SwanDispParm.f90 SwanPropvelX.f90 SwanSweepSel.f90 SwanPropvelS.f90 \
SwanTranspAc.f90 SwanTranspX.f90 SwanDiffPar.f90 SwanGSECorr.f90 SwanInterpolatePoint.f90 \
SwanInterpolateAc.f90 SwanInterpolateOutput.f90 SwanConvAccur.f90 SwanConvStopc.f90 \
SwanFindObstacles.f90 SwanCrossObstacle.f90 \
SwanComputeForce.f90 SwanIntgratSpc.f90 SwanReadfort18.f90 \
SwanSumOverNodes.f90 SwanMinOverNodes.f90 SwanMaxOverNodes.f90 \
ocpids.f ocpcre.f ocpmix.f

SWAN_OBJ:= $(patsubst %.f,   $(O_DIR)%.o, $(SWAN_SRC) )
SWAN_OBJ:= $(patsubst %.f90, $(O_DIR)%.o, $(SWAN_OBJ) )

######################## compilation rules #####################################

$(O_DIR)%.o  : %.F
#  @echo depend $<
#  @echo target $@
	$(CF) -c $(FFLAGS) -o $@  $<
	if [ "`echo *.mod`" != '*.mod' ]; then mv *.mod $(O_DIR); fi

#Casey 121126: Changed the f90 files to F90 extension to avoid conflict with SWAN.
$(O_DIR)%.o  : %.F90
	$(CF) -c $(FFLAGS) -o $@  $<
	if [ "`echo *.mod`" != '*.mod' ]; then mv *.mod $(O_DIR); fi

#Casey 121126: Added the following rules for the unstructured SWAN.
$(O_DIR)%.o  : %.f
	$(CF) -c $(FLAGS_OPT) $(FLAGS_MSC) $(FLAGS_SER) -module $(O_DIR) -o $@  $<
	if [ "`echo *.mod`" != '*.mod' ]; then mv *.mod $(O_DIR); fi
$(O_DIR)%.o  : %.f90
	$(CF) -c $(FLAGS_OPT) $(FLAGS_MSC) $(FLAGS_SER) -module $(O_DIR) -o $@  $<
	if [ "`echo *.mod`" != '*.mod' ]; then mv *.mod $(O_DIR); fi

########################## Executable Targets ##################################

.PHONY: all metis adcprep adcpost dgswem_serial dgswem winds

all :  metis adcprep adcpost dgswem_serial dgswem winds
 
winds : 
	python ./get_winds.py 

ifeq ($(MAKELEVEL),0)
   metis:
	$(MAKE) -C ../metis/Lib/ CC="$(CC)"  CFLAGS="$(CFLAGS)"
   graphchk:
	$(MAKE) -C ../metis/Programs/ CC="$(CC)"  CFLAGS="$(CFLAGS)" 
	mv ../metis/graphchk ../work
   adcprep:
	$(MAKE) BUILDTYPE=adcprep  $@            
   adcprep2:
	$(MAKE) BUILDTYPE=adcprep2 $@ 
   adcpost:
	$(MAKE) BUILDTYPE=adcpost $@            
   dgswem_serial:
	$(MAKE) BUILDTYPE=dgswem_serial  $@
   dgswem:
	$(MAKE) BUILDTYPE=dgswem $@
#Casey 121126: Added the following lines for coupling to unstructured SWAN.
   swemswan:
	@perl ../swan/switch.pl -unix -f95 -adcirc ../swan/*.ftn ../swan/*.ftn90
	$(MAKE) BUILDTYPE=swemswan $@
   dgswemswan:
	@perl ../swan/switch.pl -unix -f95 -pun -adcirc ../swan/*.ftn ../swan/*.ftn90
	$(MAKE) BUILDTYPE=dgswemswan $@
else
   adcprep::  $(O_DIR)
   adcprep2:: $(O_DIR)
   adcpost::  $(O_DIR)
   dgswem_serial::   $(O_DIR)
   dgswem::  $(O_DIR)

   adcprep ::  $(METIS1_OBJ) $(PREP_OBJ) 
	$(CF) $(FFLAGS) -o $@  $(O_DIR)*.o  $(LIBS) $(LDFLAGS)

   adcprep2 :: $(PMSG_OBJ) $(METIS2_OBJ) $(PREP_OBJ) 
	$(CF) $(FFLAGS) -o $@  $(O_DIR)*.o  $(LIBS) $(MSGLIBS) $(LDFLAGS)

   adcpost ::  $(POST_OBJ) 
	$(CF) $(FFLAGS) -o $@  $(O_DIR)*.o

   dgswem_serial ::  $(ADC_OBJ) 
	$(CF) $(FFLAGS) -o $@  $(O_DIR)*.o

   dgswem ::  $(MSG_MOBJ) $(MSG_OBJ) $(ADC_OBJ) 
	$(CF) $(FFLAGS) -o $@  $(O_DIR)*.o  $(MSGLIBS)

#Casey 121126: Added the following lines for coupling to unstructured SWAN.
   swemswan :: $(O_DIR)
   swemswan :: $(MSG_MOBJ) $(SWAN_MOBJ) $(COUP_MOBJ) $(ADC_OBJ) $(SWAN_OBJ)
	$(CF) $(FFLAGS) -o $@ $(O_DIR)*.o

   dgswemswan :: $(O_DIR)
   dgswemswan :: $(MSG_MOBJ) $(SWAN_MOBJ) $(COUP_MOBJ) $(MSG_OBJ) $(ADC_OBJ) $(SWAN_OBJ)
	$(CF) $(FFLAGS) -o $@ $(O_DIR)*.o $(MSGLIBS)
endif

########################## Misc Commands ####################################

clean:
	 rm -f odir*/*.o  odir*/*.mod 
clobber:
	rm -r -f odir*
	rm -f  graphchk adcprep adcprep2 adcpost dgswem_serial dgswem swemswan dgswemswan \
                ../metis/Lib/*.o  ../metis/libmetis.a ../metis/Programs/*.o
help:
	@echo "This makefile supports the following:"
	@echo "make all      - makes all six targets"
	@echo "make adcprep  - makes the adcprep  executable"
	@echo "make adcprep2 - makes the adcprep2  executable"
	@echo "make adcpost  - makes the adcpost  executable"
	@echo "make dgswem_serial   - makes the serial dgswem executable"
	@echo "make dgswem  - makes the dgswem executable"
#Casey 121126: Add text for the SWAN+ADCIRC targets.
	@echo "make swemswan  - makes the serial SWAN+ADCIRC executable"
	@echo "make dgswemswan - makes the parallel SWAN+ADCIRC executable"


########################## Defining the DAG  ####################################

#  adcprep & adcprep2 modules

$(O_DIR)presize.o     :  presize.F
$(O_DIR)pre_global.o  :  pre_global.F  $(PREP_MOBJ)

#  adcprep & adcprep2 

$(O_DIR)adcprep.o     :  adcprep.F  $(PGLO_MOBJ) $(PMSG_OBJ)
$(O_DIR)decomp.o      :  decomp.F   $(PGLO_MOBJ)
$(O_DIR)read_global.o :  read_global.F  $(PGLO_MOBJ)
$(O_DIR)prep.o        :  prep.F   $(PGLO_MOBJ)
$(O_DIR)interp.o      :  interp.F  
$(O_DIR)machdep.o     :  machdep.F
$(O_DIR)metis.o       :  metis.F $(PGLO_MOBJ)
$(O_DIR)metis2.o      :  metis2.F $(PGLO_MOBJ)
$(O_DIR)picomsg.o     :  picomsg.F

#  adpost modules

$(O_DIR)post_global.o :  post_global.F 

#  adpost                

$(O_DIR)adcpost.o     :  adcpost.F $(POST_MOBJ)
$(O_DIR)post.o        :  post.F  $(POST_MOBJ)
$(O_DIR)compare.o     :  compare.F  
$(O_DIR)diffmerge.o   :  diffmerge.F

#  adcirc & padcirc modules

$(O_DIR)sizes.o       			:  sizes.F
$(O_DIR)fparser.o	                :  fparser.F90 $(O_DIR)sizes.o
$(O_DIR)fparser2.o	                :  fparser2.F90 $(O_DIR)sizes.o
$(O_DIR)global.o      			:  global.F  $(O_DIR)sizes.o
$(O_DIR)messenger.o   			:  messenger.F  $(ADC_MOBJ)
$(O_DIR)messenger_elem.o		:  messenger_elem.F  $(ADC_MOBJ)
$(O_DIR)harm.o        			:  harm.F
$(O_DIR)wind.o        			:  wind.F
$(O_DIR)itpackv.o     			:  itpackv.F   $(ADC_MOBJ)
$(O_DIR)dg.o          			:  dg.F sizes.F global.F
$(O_DIR)diff45_41.o   			:  diff45_41.F
$(O_DIR)nodalattr.o   			:  nodalattr.F
#Casey 121126: Added the next line.
$(O_DIR)globalio.o                      :  globalio.F $(ADC_MOBJ)
$(O_DIR)couple2swan.o                   :  couple2swan.F $(ADC_MOBJ)
$(O_DIR)timer.o	                        :  timer.F90 $(ADC_MOBJ)
$(ODIR)read_fort_dg.o                   :  read_fort_dg.F90 $(O_DIR)sizes.o $(O_DIR)dg.o $(O_DIR)global.o

#  adcirc & padcirc

$(O_DIR)owiwind.o   			:  owiwind.F $(ADC_MOBJ)
$(O_DIR)read_input.o  			:  read_input.F $(ADC_MOBJ) $(HARM_MOBJ) $(WIND_MOBJ)
$(O_DIR)cstart.o      			:  cstart.F $(ADC_MOBJ)  $(HARM_MOBJ)  $(WIND_MOBJ)
$(O_DIR)hstart.o      			:  hstart.F $(ADC_MOBJ)  $(HARM_MOBJ)  $(WIND_MOBJ)
$(O_DIR)dgswem.o      			:  dgswem.F $(ADC_MOBJ)  $(HARM_MOBJ) $(WIND_MOBJ) 
$(O_DIR)calc_normal.o         	:  calc_normal.F $(ADC_MOBJ)
$(O_DIR)create_edge_data.o    	:  create_edge_data.F $(ADC_MOBJ)
$(O_DIR)detector.o	      	:  detector.F $(ADC_MOBJ)
$(O_DIR)DG_timestep.o         	:  DG_timestep.F $(ADC_MOBJ)
$(O_DIR)DG_hydro_timestep.o   	:  DG_hydro_timestep.F $(ADC_MOBJ)
$(O_DIR)edge_int_hydro.o      	:  edge_int_hydro.F $(ADC_MOBJ)
$(O_DIR)flow_edge_hydro.o     	:  flow_edge_hydro.F $(ADC_MOBJ)
$(O_DIR)internal_edge_hydro.o 	:  internal_edge_hydro.F $(ADC_MOBJ)
$(O_DIR)land_edge_hydro.o     	:  land_edge_hydro.F $(ADC_MOBJ)
$(O_DIR)ocean_edge_hydro.o    	:  ocean_edge_hydro.F $(ADC_MOBJ)
$(O_DIR)orthobasis_area.o     	:  orthobasis_area.F $(ADC_MOBJ)
$(O_DIR)orthobasis_edge.o     	:  orthobasis_edge.F $(ADC_MOBJ)
$(O_DIR)p_enrichment.o     	:  p_enrichment.F $(ADC_MOBJ)
$(O_DIR)prep_DG.o             	:  prep_DG.F $(ADC_MOBJ)
$(O_DIR)quad_pts_area.o       	:  quad_pts_area.F $(ADC_MOBJ)
$(O_DIR)quad_pts_edge.o       	:  quad_pts_edge.F $(ADC_MOBJ)
$(O_DIR)read_input.o          	:  read_input.F $(ADC_MOBJ)
$(O_DIR)rhs_dg_hydro.o        	:  rhs_dg_hydro.F $(ADC_MOBJ)
$(O_DIR)numerical_flux.o      	:  numerical_flux.F $(ADC_MOBJ)
$(O_DIR)radiation_edge_hydro.o	:  radiation_edge_hydro.F $(ADC_MOBJ)
$(O_DIR)tidal_potential.o     	:  tidal_potential.F $(ADC_MOBJ)
$(O_DIR)ibarrier_edge_hydro.o 	:  ibarrier_edge_hydro.F $(ADC_MOBJ)
$(O_DIR)ebarrier_edge_hydro.o       :  ebarrier_edge_hydro.F $(ADC_MOBJ)
$(O_DIR)write_results.o             :  write_results.F $(ADC_MOBJ)
$(O_DIR)pdg_debug.o                 :  pdg_debug.F  $(ADC_MOBJ) $(MSG_MOBJ)
$(O_DIR)ghost_internal_edge_hydro.o :  ghost_internal_edge_hydro.F  $(ADC_MOBJ) $(MSG_MOBJ)
$(O_DIR)slopelimiter.o              :  slopelimiter.F  $(ADC_MOBJ)
$(O_DIR)prep_slopelim.o             :  prep_slopelim.F  $(ADC_MOBJ)
$(O_DIR)ocean_edge_hydro_post.o     :  ocean_edge_hydro_post.F $(ADC_MOBJ)
$(O_DIR)ocean_flow_edge_hydro.o     :  ocean_flow_edge_hydro.F $(ADC_MOBJ)
$(O_DIR)wetdry.o                    :  wetdry.F  $(ADC_MOBJ)
$(O_DIR)LDG_hydro.o                 :  LDG_hydro.F  $(ADC_MOBJ)
$(O_DIR)met_forcing.o               :  met_forcing.F  $(ADC_MOBJ)
$(O_DIR)sta_basis.o                 :  sta_basis.F  $(ADC_MOBJ)
$(O_DIR)modal2nodal.o               :  modal2nodal.F  $(ADC_MOBJ)
$(O_DIR)write_output.o              :  write_output.F
#Casey 121126: Added the previous lines for modal2nodal.F and write_output.F.

debug:
	@echo $(ADC_MSRC) 
	@echo $(ADC_MOBJ) 

# graphchk

$(O_DIR)io.o		:  ../metis/Lib

#Casey 121126: Added the following lines for SWAN.
$(O_DIR)swmod1.o                : swmod1.f
$(O_DIR)swmod2.o                : swmod2.f
$(O_DIR)m_constants.o           : m_constants.f90
$(O_DIR)m_fileio.o              : m_fileio.f90
$(O_DIR)serv_xnl4v5.o           : serv_xnl4v5.f90
$(O_DIR)mod_xnl4v5.0            : mod_xnl4v5.f90
$(O_DIR)couple2adcirc.o         : couple2adcirc.f90
$(O_DIR)SwanGriddata.o          : SwanGriddata.f90
$(O_DIR)SwanGridobjects.o       : SwanGridobjects.f90
$(O_DIR)SwanCompdata.o          : SwanCompdata.f90
$(O_DIR)swanmain.o              : swanmain.f
$(O_DIR)swanpre1.o              : swanpre1.f
$(O_DIR)swanpre2.o              : swanpre2.f
$(O_DIR)swancom1.o              : swancom1.f
$(O_DIR)swancom2.o              : swancom2.f
$(O_DIR)swancom3.o              : swancom3.f
$(O_DIR)swancom4.o              : swancom4.f
$(O_DIR)swancom5.o              : swancom5.f
$(O_DIR)swanout1.o              : swanout1.f
$(O_DIR)swanout2.o              : swanout2.f
$(O_DIR)swanser.o               : swanser.f
$(O_DIR)swanparll.o             : swanparll.f
$(O_DIR)SwanReadGrid.o          : SwanReadGrid.f90
$(O_DIR)SwanReadADCGrid.o       : SwanReadADCGrid.f90
$(O_DIR)SwanReadTriangleGrid.o  : SwanReadTriangleGrid.f90
$(O_DIR)SwanReadEasymeshGrid.o  : SwanReadEasymeshGrid.f90
$(O_DIR)SwanInitCompGrid.o      : SwanInitCompGrid.f90
$(O_DIR)SwanCheckGrid.o         : SwanCheckGrid.f90
$(O_DIR)SwanCreateEdges.o       : SwanCreateEdges.f90
$(O_DIR)SwanGridTopology.o      : SwanGridTopology.f90
$(O_DIR)SwanGridVert.o          : SwanGridVert.f90
$(O_DIR)SwanGridCell.o          : SwanGridCell.f90
$(O_DIR)SwanGridFace.o          : SwanGridFace.f90
$(O_DIR)SwanPrintGridInfo.o     : SwanPrintGridInfo.f90
$(O_DIR)SwanFindPoint.o         : SwanFindPoint.f90
$(O_DIR)SwanPointinMesh.o       : SwanPointinMesh.f90
$(O_DIR)SwanBpntlist.o          : SwanBpntlist.f90
$(O_DIR)SwanPrepComp.o          : SwanPrepComp.f90
$(O_DIR)SwanVertlist.o          : SwanVertlist.f90
$(O_DIR)SwanCompUnstruc.o       : SwanCompUnstruc.f90
$(O_DIR)SwanDispParm.o          : SwanDispParm.f90
$(O_DIR)SwanPropvelX.o          : SwanPropvelX.f90
$(O_DIR)SwanSweepSel.o          : SwanSweepSel.f90
$(O_DIR)SwanPropvelS.o          : SwanPropvelS.f90
$(O_DIR)SwanTranspAc.o          : SwanTranspAc.f90
$(O_DIR)SwanTranspX.o           : SwanTranspX.f90
$(O_DIR)SwanDiffPar.o           : SwanDiffPar.f90
$(O_DIR)SwanGSECorr.o           : SwanGSECorr.f90
$(O_DIR)SwanInterpolatePoint.o  : SwanInterpolatePoint.f90
$(O_DIR)SwanInterpolateAc.o     : SwanInterpolateAc.f90
$(O_DIR)SwanInterpolateOutput.o : SwanInterpolateOutput.f90
$(O_DIR)SwanConvAccur.o         : SwanConvAccur.f90
$(O_DIR)SwanConvStopc.o         : SwanConvStopc.f90
$(O_DIR)SwanFindObstacles.o     : SwanFindObstacles.f90
$(O_DIR)SwanCrossObstacle.o     : SwanCrossObstacle.f90
$(O_DIR)SwanComputeForce.o      : SwanComputeForce.f90
$(O_DIR)SwanIntgratSpc.o        : SwanIntgratSpc.f90
$(O_DIR)SwanReadfort18.o        : SwanReadfort18.f90
$(O_DIR)SwanSumOverNodes.o      : SwanSumOverNodes.f90
$(O_DIR)SwanMinOverNodes.o      : SwanMinOverNodes.f90
$(O_DIR)SwanMaxOverNodes.o      : SwanMaxOverNodes.f90
$(O_DIR)ocpids.o                : ocpids.f
$(O_DIR)ocpcre.o                : ocpcre.f
$(O_DIR)ocpmix.o                : ocpmix.f

