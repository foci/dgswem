SHELL:=/bin/sh
#
#  Makefile to Build DGSWEM and its pre/post-processor 

########################################################################
#  Get Canonical Machine NAME from config.guess
#1
NAME     := $(shell ./config.guess)
LIST     := $(subst -, ,$(NAME))
MACHINE  := $(word 1, $(LIST))
VENDOR   := $(word 2, $(LIST))
OS       := $(subst  $(MACHINE)-$(VENDOR)-,,$(strip $(NAME)))

#include cmplrflags.mk
include cmplrflags.zachtop.mk

# Debugging 
$(info $$NAME is [${NAME}])
$(info $$FC is [${FC}])
#
######################## Target Specific Rules ###################################

ifeq ($(BUILDTYPE),dgswem_serial)
  CF:= $(PPFC)
  O_DIR:=odir4/
  FFLAGS:= $(FFLAGS3) $(DA) $(IMODS)$(O_DIR)
  VPATH:=  ../src 
  MSG_MOBJ:=
endif

$(O_DIR):
	mkdir -p $@

######################### Module Source, Object, & Mod Files ######################

PREP_MSRC  =  presizes.F 
PGLO_MSRC  =  pre_global.F  
POST_MSRC  =  post_global.F
ADC_MSRC   =  sizes.F diff45_41.F global.f95 dg.f95 global_3dvs.F nodalattr.F fparser.F90 fparser2.F90
HARM_MSRC  =  harm.F  
WIND_MSRC  =  wind.F  
#Casey 121126: Add file for coupling to unstructured SWAN.
COUP_MSRC  =  couple2swan.F

PREP_MOBJ:= $(patsubst %.F, $(O_DIR)%.o, $(PREP_MSRC) )
PGLO_MOBJ:= $(patsubst %.F, $(O_DIR)%.o, $(PGLO_MSRC) )
POST_MOBJ:= $(patsubst %.F, $(O_DIR)%.o, $(POST_MSRC) )
ADC_MOBJ := $(patsubst %.F, $(O_DIR)%.o, $(ADC_MSRC)  )
ADC_MOBJ := $(patsubst %.F90, $(O_DIR)%.o, $(ADC_MOBJ)  )
ADC_MOBJ := $(patsubst %.f95, $(O_DIR)%.o, $(ADC_MOBJ)  )
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

######################## compilation rules #####################################

$(O_DIR)%.o  : %.f95
#  @echo depend $<
#  @echo target $@
	$(CF) -c $(FFLAGS) -o $@  $<
	if [ "`echo *.mod`" != '*.mod' ]; then mv *.mod $(O_DIR); fi

$(O_DIR)%.o  : %.F
#  @echo depend $<
#  @echo target $@
	$(CF) -c $(FFLAGS) -o $@  $<
	if [ "`echo *.mod`" != '*.mod' ]; then mv *.mod $(O_DIR); fi

#Casey 121126: Changed the f90 files to F90 extension to avoid conflict with SWAN.
$(O_DIR)%.o  : %.F90
	$(CF) -c $(FFLAGS) -o $@  $<
	if [ "`echo *.mod`" != '*.mod' ]; then mv *.mod $(O_DIR); fi

########################## Executable Targets ##################################

.PHONY: all metis adcprep adcpost dgswem_serial dgswem winds

all :  adcprep adcpost dgswem_serial dgswem winds
 
winds : 
	python ./get_winds.py 

ifeq ($(MAKELEVEL),0)
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
	@echo "make dgswem_serial   - makes the serial dgswem executable"

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
$(O_DIR)global.o      			:  global.f95  $(O_DIR)sizes.o
$(O_DIR)messenger.o   			:  messenger.F  $(ADC_MOBJ)
$(O_DIR)messenger_elem.o		:  messenger_elem.F  $(ADC_MOBJ)
$(O_DIR)harm.o        			:  harm.F
$(O_DIR)wind.o        			:  wind.F
$(O_DIR)itpackv.o     			:  itpackv.F   $(ADC_MOBJ)
$(O_DIR)dg.o          			:  dg.f95 sizes.F global.f95
$(O_DIR)diff45_41.o   			:  diff45_41.F
$(O_DIR)nodalattr.o   			:  nodalattr.F
#Casey 121126: Added the next line.
$(O_DIR)globalio.o                      :  globalio.F $(ADC_MOBJ)
$(O_DIR)couple2swan.o                   :  couple2swan.F $(ADC_MOBJ)

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

