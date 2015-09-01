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

include cmplrflags.mk
#include cmplrflags.zachtop.mk

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

ifeq ($(BUILDTYPE),adcpost)
  CF:= $(PPFC)
  O_DIR:=odir3/
  FFLAGS:= $(FFLAGS1) $(DPRE) $(IMODS)$(O_DIR) 
  VPATH :=  ../prep
endif


$(O_DIR):
	mkdir -p $@

######################### Module Source, Object, & Mod Files ######################

PREP_MSRC  =  presizes.F90 
PGLO_MSRC  =  pre_global.F90  
POST_MSRC  =  post_global.F
ADC_MSRC   =  sizes.F90 diff45_41.F90 global.F90 dg.F90 nodalattr.F90 fort_dg_module.F90 
HARM_MSRC  =  harm.F90  
WIND_MSRC  =  wind.F90  
#Casey 121126: Add file for coupling to unstructured SWAN.
COUP_MSRC  =  couple2swan.F90

PREP_MOBJ:= $(patsubst %.F90, $(O_DIR)%.o, $(PREP_MSRC) )
PGLO_MOBJ:= $(patsubst %.F90, $(O_DIR)%.o, $(PGLO_MSRC) )
POST_MOBJ:= $(patsubst %.F, $(O_DIR)%.o, $(POST_MSRC) )
ADC_MOBJ := $(patsubst %.F90, $(O_DIR)%.o, $(ADC_MSRC)  )
ADC_MOBJ := $(patsubst %.F9090, $(O_DIR)%.o, $(ADC_MOBJ)  )
ADC_MOBJ := $(patsubst %.F90, $(O_DIR)%.o, $(ADC_MOBJ)  )
HARM_MOBJ := $(patsubst %.F90, $(O_DIR)%.o, $(HARM_MSRC)  )
WIND_MOBJ := $(patsubst %.F90, $(O_DIR)%.o, $(WIND_MSRC)  )
#Casey 121126: Add rules for coupling to unstructured SWAN.
COUP_MOBJ := $(patsubst %.F90, $(O_DIR)%.o, $(COUP_MSRC)  )


############################# Source & Object Files ##############################

METIS1_SRC =  metis.F90
METIS2_SRC =  metis2.F90
PREP_SRC   =  adcprep.F90  decomp.F90   read_global.F90   prep.F90   interp.F90  machdep.F90
POST_SRC   =  adcpost.F    post.F  compare.F  diffmerge.F
ADC_SRC    =  dgswem.F90 read_input.F90 cstart.F90 hstart.F90 radiation_edge_hydro.F90 \
	ibarrier_edge_hydro.F90 ebarrier_edge_hydro.F90  \
	calc_normal.F90 create_edge_data.F90 \
	DG_timestep.F90 DG_hydro_timestep.F90 DG_timestep_advance.F90\
	flow_edge_hydro.F90 \
	internal_edge_hydro.F90 land_edge_hydro.F90 \
	ocean_edge_hydro.F90  \
	orthobasis_area.F90 orthobasis_edge.F90 p_enrichment.F90 prep_DG.F90 quad_pts_area.F90 \
	quad_pts_edge.F90 radiation_edge_hydro.F90 read_input.F90 rhs_dg_hydro.F90 \
	numerical_flux.F90 tidal_potential.F90 write_results.F90 \
	slopelimiter.F90 prep_slopelim.F90 \
	wetdry.F90 LDG_hydro.F90 met_forcing.F90 sta_basis.F90 make_dirname.F90 nodalattr_subs.F90 dgswem_init.F90 cwrappers.F90 make_dirname.F90 read_fort_dg.F90 message_table_hpx.F90 init_fileunits.F90 get_neighbors.F90 messenger_hpx.F90 hpx_read_n_domains.F90 lgd_yield_subdomain_coord.F90
#Casey 121126: Added modal2nodal.F90 and write_output.F90 to the previous line.

METIS1_OBJ:= $(patsubst %.F90, $(O_DIR)%.o, $(METIS1_SRC) )
METIS2_OBJ:= $(patsubst %.F90, $(O_DIR)%.o, $(METIS2_SRC) )
PREP_OBJ:= $(patsubst %.F90, $(O_DIR)%.o, $(PREP_SRC) )
POST_OBJ:= $(patsubst %.F, $(O_DIR)%.o, $(POST_SRC) )
ADC_OBJ:= $(patsubst %.F90, $(O_DIR)%.o, $(ADC_SRC) )

######################## compilation rules #####################################

$(O_DIR)%.o  : %.F90
#  @echo depend $<
#  @echo target $@
	$(CF) -c  $(FFLAGS) -o $@  $<
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

$(O_DIR)presize.o     :  presize.F90
$(O_DIR)pre_global.o  :  pre_global.F90  $(PREP_MOBJ)

#  adcprep & adcprep2 

$(O_DIR)adcprep.o     :  adcprep.F90  $(PGLO_MOBJ) $(PMSG_OBJ)
$(O_DIR)decomp.o      :  decomp.F90   $(PGLO_MOBJ)
$(O_DIR)read_global.o :  read_global.F90  $(PGLO_MOBJ)
$(O_DIR)prep.o        :  prep.F90   $(PGLO_MOBJ)
$(O_DIR)interp.o      :  interp.F90  
$(O_DIR)machdep.o     :  machdep.F90
$(O_DIR)metis.o       :  metis.F90 $(PGLO_MOBJ)
$(O_DIR)metis2.o      :  metis2.F90 $(PGLO_MOBJ)
$(O_DIR)picomsg.o     :  picomsg.F90

#  adpost modules

$(O_DIR)post_global.o :  post_global.F 

#  adpost                

$(O_DIR)adcpost.o     :  adcpost.F $(POST_MOBJ)
$(O_DIR)post.o        :  post.F  $(POST_MOBJ)
$(O_DIR)compare.o     :  compare.F  
$(O_DIR)diffmerge.o   :  diffmerge.F

#  adcirc & padcirc modules

$(O_DIR)sizes.o       			:  sizes.F90
$(O_DIR)global.o      			:  global.F90  $(O_DIR)sizes.o
$(O_DIR)messenger.o   			:  messenger.F90  $(ADC_MOBJ)
$(O_DIR)messenger_elem.o		:  messenger_elem.F90  $(ADC_MOBJ)
#$(O_DIR)harm.o        			:  harm.F90
$(O_DIR)wind.o        			:  wind.F90
$(O_DIR)itpackv.o     			:  itpackv.F90   $(ADC_MOBJ)
$(O_DIR)dg.o          			:  dg.F90 sizes.F90 global.F90
$(O_DIR)diff45_41.o   			:  diff45_41.F90
$(O_DIR)nodalattr.o   			:  nodalattr.F90
$(O_DIR)fort_dg_module.o                :  fort_dg_module.F90
#Casey 121126: Added the next line.
$(O_DIR)globalio.o                      :  globalio.F90 $(ADC_MOBJ)
$(O_DIR)couple2swan.o                   :  couple2swan.F90 $(ADC_MOBJ)

$(ODIR)read_fort_dg.o                   :  read_fort_dg.F90 $(O_DIR)sizes.o $(O_DIR)dg.o $(O_DIR)global.o $(O_DIR) fort_dg_module.o

#  adcirc & padcirc

$(O_DIR)make_dirname.o              :  make_dirname.F90 $(ADC_MOBJ)
#$(O_DIR)owiwind.o   			:  owiwind.F90 $(ADC_MOBJ)
$(O_DIR)read_input.o  			:  read_input.F90 $(ADC_MOBJ) $(WIND_MOBJ)
$(O_DIR)cstart.o      			:  cstart.F90 $(ADC_MOBJ)  $(WIND_MOBJ)
$(O_DIR)hstart.o      			:  hstart.F90 $(ADC_MOBJ)  $(WIND_MOBJ)
$(O_DIR)dgswem.o      			:  dgswem.F90 $(ADC_MOBJ) $(WIND_MOBJ) 
$(O_DIR)calc_normal.o         	:  calc_normal.F90 $(ADC_MOBJ)
$(O_DIR)create_edge_data.o    	:  create_edge_data.F90 $(ADC_MOBJ)
$(O_DIR)DG_timestep.o         	:  DG_timestep.F90 $(ADC_MOBJ)
$(O_DIR)DG_hydro_timestep.o   	:  DG_hydro_timestep.F90 $(ADC_MOBJ)
$(O_DIR)DG_timestep_advance.o  :  DG_timestep_advance.F90 $(ADC_MOBJ)
$(O_DIR)flow_edge_hydro.o     	:  flow_edge_hydro.F90 $(ADC_MOBJ)
$(O_DIR)internal_edge_hydro.o 	:  internal_edge_hydro.F90 $(ADC_MOBJ)
$(O_DIR)land_edge_hydro.o     	:  land_edge_hydro.F90 $(ADC_MOBJ)
$(O_DIR)ocean_edge_hydro.o    	:  ocean_edge_hydro.F90 $(ADC_MOBJ)
$(O_DIR)orthobasis_area.o     	:  orthobasis_area.F90 $(ADC_MOBJ)
$(O_DIR)orthobasis_edge.o     	:  orthobasis_edge.F90 $(ADC_MOBJ)
$(O_DIR)p_enrichment.o     	:  p_enrichment.F90 $(ADC_MOBJ)
$(O_DIR)prep_DG.o             	:  prep_DG.F90 $(ADC_MOBJ)
$(O_DIR)quad_pts_area.o       	:  quad_pts_area.F90 $(ADC_MOBJ)
$(O_DIR)quad_pts_edge.o       	:  quad_pts_edge.F90 $(ADC_MOBJ)
$(O_DIR)read_input.o          	:  read_input.F90 $(ADC_MOBJ)
$(O_DIR)rhs_dg_hydro.o        	:  rhs_dg_hydro.F90 $(ADC_MOBJ)
$(O_DIR)numerical_flux.o      	:  numerical_flux.F90 $(ADC_MOBJ)
$(O_DIR)radiation_edge_hydro.o	:  radiation_edge_hydro.F90 $(ADC_MOBJ)
$(O_DIR)tidal_potential.o     	:  tidal_potential.F90 $(ADC_MOBJ)
$(O_DIR)ibarrier_edge_hydro.o 	:  ibarrier_edge_hydro.F90 $(ADC_MOBJ)
$(O_DIR)ebarrier_edge_hydro.o       :  ebarrier_edge_hydro.F90 $(ADC_MOBJ)
$(O_DIR)write_results.o             :  write_results.F90 $(ADC_MOBJ)
$(O_DIR)pdg_debug.o                 :  pdg_debug.F90  $(ADC_MOBJ) $(MSG_MOBJ)
$(O_DIR)ghost_internal_edge_hydro.o :  ghost_internal_edge_hydro.F90  $(ADC_MOBJ) $(MSG_MOBJ)
$(O_DIR)slopelimiter.o              :  slopelimiter.F90  $(ADC_MOBJ)
$(O_DIR)prep_slopelim.o             :  prep_slopelim.F90  $(ADC_MOBJ)
$(O_DIR)wetdry.o                    :  wetdry.F90  $(ADC_MOBJ)
$(O_DIR)LDG_hydro.o                 :  LDG_hydro.F90  $(ADC_MOBJ)
$(O_DIR)met_forcing.o               :  met_forcing.F90  $(ADC_MOBJ)
$(O_DIR)sta_basis.o                 :  sta_basis.F90  $(ADC_MOBJ)
$(O_DIR)nodalattr_subs.o            : nodalattr_subs.F90 $(ADC_MOBJ)
$(0_DIR)dgswem_init.o               : dgswem_init.F90 $(ADC_MOBJ)
$(0_DIR)cwrappers.o               : cwrappers.F90 $(ADC_MOBJ)


#Casey 121126: Added the previous lines for modal2nodal.F and write_output.F.

debug:
	@echo $(ADC_MSRC) 
	@echo $(ADC_MOBJ) 

# graphchk

$(O_DIR)io.o		:  ../metis/Lib

