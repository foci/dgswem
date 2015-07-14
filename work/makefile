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
$(info $$PFC is [${PFC}])
$(info MAKELEVEL is [${MAKELEVEL}])
#
######################## Target Specific Rules ###################################

ifeq ($(BUILDTYPE),dgswem_serial)
  CF:= $(PPFC)
  O_DIR:=odir4/
  FFLAGS:= $(FFLAGS3) $(DA) $(IMODS)$(O_DIR)
  VPATH:=  ../src 
  MSG_MOBJ:=
endif

ifeq ($(BUILDTYPE),dgswem)
  CF:= $(PFC)
  O_DIR:=odir5/
  FFLAGS:= $(FFLAGS4) $(DP) $(IMODS)$(O_DIR)
  VPATH :=  ../src 
  MSG_MOBJ:= $(O_DIR)messenger_elem.o $(O_DIR)messenger.o
  MSG_OBJ:= 
endif

ifeq ($(BUILDTYPE),dgswem_hpx)
  CF:= $(PPFC)
  O_DIR:=odir4/
  FFLAGS:= -DHPX $(FFLAGS3) $(DA) $(IMODS)$(O_DIR)
  VPATH:=  ../src 
  MSG_MOBJ:=
endif


$(O_DIR):
	mkdir -p $@

######################### Module Source, Object, & Mod Files ######################

PREP_MSRC  =  presizes.f95 
PGLO_MSRC  =  pre_global.f95  
POST_MSRC  =  post_global.f95
ADC_MSRC   =  sizes.f95 diff45_41.f95 global.f95 dg.f95 nodalattr.f95 fort_dg_module.f95 
HARM_MSRC  =  harm.f95  
WIND_MSRC  =  wind.f95  
#Casey 121126: Add file for coupling to unstructured SWAN.
COUP_MSRC  =  couple2swan.f95

PREP_MOBJ:= $(patsubst %.f95, $(O_DIR)%.o, $(PREP_MSRC) )
PGLO_MOBJ:= $(patsubst %.f95, $(O_DIR)%.o, $(PGLO_MSRC) )
POST_MOBJ:= $(patsubst %.f95, $(O_DIR)%.o, $(POST_MSRC) )
ADC_MOBJ := $(patsubst %.f95, $(O_DIR)%.o, $(ADC_MSRC)  )
ADC_MOBJ := $(patsubst %.f9590, $(O_DIR)%.o, $(ADC_MOBJ)  )
ADC_MOBJ := $(patsubst %.f95, $(O_DIR)%.o, $(ADC_MOBJ)  )
HARM_MOBJ := $(patsubst %.f95, $(O_DIR)%.o, $(HARM_MSRC)  )
WIND_MOBJ := $(patsubst %.f95, $(O_DIR)%.o, $(WIND_MSRC)  )
#Casey 121126: Add rules for coupling to unstructured SWAN.
COUP_MOBJ := $(patsubst %.f95, $(O_DIR)%.o, $(COUP_MSRC)  )


############################# Source & Object Files ##############################

METIS1_SRC =  metis.f95
METIS2_SRC =  metis2.f95
PREP_SRC   =  adcprep.f95  decomp.f95   read_global.f95   prep.f95   interp.f95  machdep.f95
POST_SRC   =  adcpost.f95    post.f95  compare.f95  diffmerge.f95
ADC_SRC    =  dgswem.f95 read_input.f95 cstart.f95 hstart.f95 radiation_edge_hydro.f95 \
	ibarrier_edge_hydro.f95 ebarrier_edge_hydro.f95  \
	calc_normal.f95 create_edge_data.f95 \
	DG_timestep.f95 DG_hydro_timestep.f95 \
	flow_edge_hydro.f95 \
	internal_edge_hydro.f95 land_edge_hydro.f95 \
	ocean_edge_hydro.f95  \
	orthobasis_area.f95 orthobasis_edge.f95 p_enrichment.f95 prep_DG.f95 quad_pts_area.f95 \
	quad_pts_edge.f95 radiation_edge_hydro.f95 read_input.f95 rhs_dg_hydro.f95 \
	numerical_flux.f95 tidal_potential.f95 write_results.f95 \
	slopelimiter.f95 prep_slopelim.f95 \
	wetdry.f95 LDG_hydro.f95 met_forcing.f95 sta_basis.f95 make_dirname.f95 nodalattr_subs.f95 dgswem_init.f95 cwrappers.f95 make_dirname.f95 read_fort_dg.f95 message_table_hpx.f95 init_fileunits.f95 get_neighbors.f95
#Casey 121126: Added modal2nodal.f95 and write_output.f95 to the previous line.

METIS1_OBJ:= $(patsubst %.f95, $(O_DIR)%.o, $(METIS1_SRC) )
METIS2_OBJ:= $(patsubst %.f95, $(O_DIR)%.o, $(METIS2_SRC) )
PREP_OBJ:= $(patsubst %.f95, $(O_DIR)%.o, $(PREP_SRC) )
POST_OBJ:= $(patsubst %.f95, $(O_DIR)%.o, $(POST_SRC) )
ADC_OBJ:= $(patsubst %.f95, $(O_DIR)%.o, $(ADC_SRC) )

######################## compilation rules #####################################

$(O_DIR)%.o  : %.f95
#  @echo depend $<
#  @echo target $@
	$(CF) -c -cpp $(FFLAGS) -o $@  $<
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

.PHONY: all metis adcprep adcpost dgswem_serial dgswem_hpx dgswem winds

all :  adcprep adcpost dgswem_serial dgswem_hpx dgswem winds
 
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
   dgswem_hpx:
	$(MAKE) BUILDTYPE=dgswem_hpx  $@
   dgswem:
	$(MAKE) BUILDTYPE=dgswem $@

else
   adcprep::  $(O_DIR)
   adcprep2:: $(O_DIR)
   adcpost::  $(O_DIR)
   dgswem_serial::   $(O_DIR)
   dgswem_hpx::   $(O_DIR)
   dgswem::  $(O_DIR)

   adcprep ::  $(METIS1_OBJ) $(PREP_OBJ) 
	$(CF) $(FFLAGS) -o $@  $(O_DIR)*.o  $(LIBS) $(LDFLAGS)

   adcprep2 :: $(PMSG_OBJ) $(METIS2_OBJ) $(PREP_OBJ) 
	$(CF) $(FFLAGS) -o $@  $(O_DIR)*.o  $(LIBS) $(MSGLIBS) $(LDFLAGS)

   adcpost ::  $(POST_OBJ) 
	$(CF) $(FFLAGS) -o $@  $(O_DIR)*.o

   dgswem_serial ::  $(ADC_OBJ) 
	$(CF) $(FFLAGS) -o $@  $(O_DIR)*.o

   dgswem_hpx ::  $(ADC_OBJ) 
	$(CF) $(FFLAGS) -o $@  $(O_DIR)*.o

   dgswem ::  $(MSG_MOBJ) $(MSG_OBJ) $(ADC_OBJ) 
	$(CF) $(FFLAGS) -o $@  $(O_DIR)*.o  $(MSGLIBS)
endif


########################## Misc Commands ####################################

clean:
	 rm -f odir*/*.o  odir*/*.mod 
clobber:
	rm -r -f odir*
	rm -f  graphchk adcprep adcprep2 adcpost dgswem_serial dgswem swemswan dgswemswan dgswem_hpx \
                ../metis/Lib/*.o  ../metis/libmetis.a ../metis/Programs/*.o
help:
	@echo "This makefile supports the following:"
	@echo "make dgswem_serial   - makes the serial dgswem executable"

########################## Defining the DAG  ####################################

#  adcprep & adcprep2 modules

$(O_DIR)presize.o     :  presize.f95
$(O_DIR)pre_global.o  :  pre_global.f95  $(PREP_MOBJ)

#  adcprep & adcprep2 

$(O_DIR)adcprep.o     :  adcprep.f95  $(PGLO_MOBJ) $(PMSG_OBJ)
$(O_DIR)decomp.o      :  decomp.f95   $(PGLO_MOBJ)
$(O_DIR)read_global.o :  read_global.f95  $(PGLO_MOBJ)
$(O_DIR)prep.o        :  prep.f95   $(PGLO_MOBJ)
$(O_DIR)interp.o      :  interp.f95  
$(O_DIR)machdep.o     :  machdep.f95
$(O_DIR)metis.o       :  metis.f95 $(PGLO_MOBJ)
$(O_DIR)metis2.o      :  metis2.f95 $(PGLO_MOBJ)
$(O_DIR)picomsg.o     :  picomsg.f95

#  adpost modules

$(O_DIR)post_global.o :  post_global.f95 

#  adpost                

$(O_DIR)adcpost.o     :  adcpost.f95 $(POST_MOBJ)
$(O_DIR)post.o        :  post.f95  $(POST_MOBJ)
$(O_DIR)compare.o     :  compare.f95  
$(O_DIR)diffmerge.o   :  diffmerge.f95

#  adcirc & padcirc modules

$(O_DIR)sizes.o       			:  sizes.f95
$(O_DIR)global.o      			:  global.f95  $(O_DIR)sizes.o
$(O_DIR)messenger.o   			:  messenger.f95  $(ADC_MOBJ)
$(O_DIR)messenger_elem.o		:  messenger_elem.f95  $(ADC_MOBJ)
#$(O_DIR)harm.o        			:  harm.f95
$(O_DIR)wind.o        			:  wind.f95
$(O_DIR)itpackv.o     			:  itpackv.f95   $(ADC_MOBJ)
$(O_DIR)dg.o          			:  dg.f95 sizes.f95 global.f95
$(O_DIR)diff45_41.o   			:  diff45_41.f95
$(O_DIR)nodalattr.o   			:  nodalattr.f95
$(O_DIR)fort_dg_module.o                :  fort_dg_module.f95
#Casey 121126: Added the next line.
$(O_DIR)globalio.o                      :  globalio.f95 $(ADC_MOBJ)
$(O_DIR)couple2swan.o                   :  couple2swan.f95 $(ADC_MOBJ)

$(ODIR)read_fort_dg.o                   :  read_fort_dg.f95 $(O_DIR)sizes.o $(O_DIR)dg.o $(O_DIR)global.o $(O_DIR) fort_dg_module.o

#  adcirc & padcirc

$(O_DIR)make_dirname.o              :  make_dirname.f95 $(ADC_MOBJ)
#$(O_DIR)owiwind.o   			:  owiwind.f95 $(ADC_MOBJ)
$(O_DIR)read_input.o  			:  read_input.f95 $(ADC_MOBJ) $(WIND_MOBJ)
$(O_DIR)cstart.o      			:  cstart.f95 $(ADC_MOBJ)  $(WIND_MOBJ)
$(O_DIR)hstart.o      			:  hstart.f95 $(ADC_MOBJ)  $(WIND_MOBJ)
$(O_DIR)dgswem.o      			:  dgswem.f95 $(ADC_MOBJ) $(WIND_MOBJ) 
$(O_DIR)calc_normal.o         	:  calc_normal.f95 $(ADC_MOBJ)
$(O_DIR)create_edge_data.o    	:  create_edge_data.f95 $(ADC_MOBJ)
$(O_DIR)DG_timestep.o         	:  DG_timestep.f95 $(ADC_MOBJ)
$(O_DIR)DG_hydro_timestep.o   	:  DG_hydro_timestep.f95 $(ADC_MOBJ)
$(O_DIR)flow_edge_hydro.o     	:  flow_edge_hydro.f95 $(ADC_MOBJ)
$(O_DIR)internal_edge_hydro.o 	:  internal_edge_hydro.f95 $(ADC_MOBJ)
$(O_DIR)land_edge_hydro.o     	:  land_edge_hydro.f95 $(ADC_MOBJ)
$(O_DIR)ocean_edge_hydro.o    	:  ocean_edge_hydro.f95 $(ADC_MOBJ)
$(O_DIR)orthobasis_area.o     	:  orthobasis_area.f95 $(ADC_MOBJ)
$(O_DIR)orthobasis_edge.o     	:  orthobasis_edge.f95 $(ADC_MOBJ)
$(O_DIR)p_enrichment.o     	:  p_enrichment.f95 $(ADC_MOBJ)
$(O_DIR)prep_DG.o             	:  prep_DG.f95 $(ADC_MOBJ)
$(O_DIR)quad_pts_area.o       	:  quad_pts_area.f95 $(ADC_MOBJ)
$(O_DIR)quad_pts_edge.o       	:  quad_pts_edge.f95 $(ADC_MOBJ)
$(O_DIR)read_input.o          	:  read_input.f95 $(ADC_MOBJ)
$(O_DIR)rhs_dg_hydro.o        	:  rhs_dg_hydro.f95 $(ADC_MOBJ)
$(O_DIR)numerical_flux.o      	:  numerical_flux.f95 $(ADC_MOBJ)
$(O_DIR)radiation_edge_hydro.o	:  radiation_edge_hydro.f95 $(ADC_MOBJ)
$(O_DIR)tidal_potential.o     	:  tidal_potential.f95 $(ADC_MOBJ)
$(O_DIR)ibarrier_edge_hydro.o 	:  ibarrier_edge_hydro.f95 $(ADC_MOBJ)
$(O_DIR)ebarrier_edge_hydro.o       :  ebarrier_edge_hydro.f95 $(ADC_MOBJ)
$(O_DIR)write_results.o             :  write_results.f95 $(ADC_MOBJ)
$(O_DIR)pdg_debug.o                 :  pdg_debug.f95  $(ADC_MOBJ) $(MSG_MOBJ)
$(O_DIR)ghost_internal_edge_hydro.o :  ghost_internal_edge_hydro.f95  $(ADC_MOBJ) $(MSG_MOBJ)
$(O_DIR)slopelimiter.o              :  slopelimiter.f95  $(ADC_MOBJ)
$(O_DIR)prep_slopelim.o             :  prep_slopelim.f95  $(ADC_MOBJ)
$(O_DIR)wetdry.o                    :  wetdry.f95  $(ADC_MOBJ)
$(O_DIR)LDG_hydro.o                 :  LDG_hydro.f95  $(ADC_MOBJ)
$(O_DIR)met_forcing.o               :  met_forcing.f95  $(ADC_MOBJ)
$(O_DIR)sta_basis.o                 :  sta_basis.f95  $(ADC_MOBJ)
$(O_DIR)nodalattr_subs.o            : nodalattr_subs.f95 $(ADC_MOBJ)
$(0_DIR)dgswem_init.o               : dgswem_init.f95 $(ADC_MOBJ)
$(0_DIR)cwrappers.o               : cwrappers.f95 $(ADC_MOBJ)


#Casey 121126: Added the previous lines for modal2nodal.F and write_output.F.

debug:
	@echo $(ADC_MSRC) 
	@echo $(ADC_MOBJ) 

# graphchk

$(O_DIR)io.o		:  ../metis/Lib

