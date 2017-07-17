SHELL=/bin/sh
#
#  PC Makefile to build a single processor version of ADCIRC
#  Known to work with Visual Fortran 6.6
#
#  To use this file, rename it to: makefile
#
#  Options:
#
#  To create a dgswem executable, type:        nmake dgswem
#  To create preprocessed module files, type:  nmake preprocessed_mods
#  To create preprocessed source files, type:  nmake preprocessed_srcs
#  To eliminate all partial files, type:       nmake clean
#


##################### Preprocessor and Compiler Flags #####################
 
#DA   =  -P  -DREAL4  -DCSCA         #Preprocessor flags for a single precision, scalar code
DA   =  -P  -DREAL8  -DCSCA  -DRKSSP         #Preprocessor flags for a double presision, scalar code
#DA   =  -P  -DREAL4  -DCVEC         #Preprocessor flags for a single presision, vector code
#DA   =  -P  -DREAL8  -DCVEC         #Preprocessor flags for a double presision, vector code

FC         = f90

FFLAGS     = /inline:all /fast /cxml #Compiler flags  -DVF:Visual Fortran
#FFLAGS    = -g -C -f -dalign        #These are left over Sun compiler flags

####################### Source File Dependencies ############################
#
# this assigns a preprocessed .f file in the current (work) directory with the
# corresponding target name
#

dgswem.f: ../src/dgswem.F
	fpp $(DA) ../src/dgswem.F  ../work/dgswem.f

calc_normal.f: ../src/calc_normal.F
	fpp $(DA)  ../src/calc_normal.F   ../work/calc_normal.f

create_edge_data.f: ../src/create_edge_data.F
	fpp $(DA) ../src/create_edge_data.F  ../work/create_edge_data.f

cstart.f: ../src/cstart.F
	fpp $(DA)  ../src/cstart.F   ../work/cstart.f

detector.f: ../src/detector.F
	fpp $(DA)  ../src/detector.F   ../work/detector.f

dg.f: ../src/dg.F
	fpp $(DA)  ../src/dg.F   ../work/dg.f

DG_hydro_timestep.f: ../src/DG_hydro_timestep.F
	fpp $(DA)  ../src/DG_hydro_timestep.F   ../work/DG_hydro_timestep.f

edge_int_hydro.f: ../src/edge_int_hydro.F
	fpp $(DA)  ../src/edge_int_hydro.F   ../work/edge_int_hydro.f

flow_edge_hydro.f: ../src/flow_edge_hydro.F
	fpp $(DA)  ../src/flow_edge_hydro.F   ../work/flow_edge_hydro.f

global.f: ../src/global.F    
	fpp $(DA)  ../src/global.F ../work/global.f

harm.f: ../src/harm.F
	fpp $(DA)  ../src/harm.F  ../work/harm.f

hstart.f: ../src/hstart.F
	fpp $(DA)  ../src/hstart.F   ../work/hstart.f

ibarrier_edge_hydro.f: ../src/ibarrier_edge_hydro.F
	fpp $(DA) ../src/ibarrier_edge_hydro.F   ../work/ibarrier_edge_hydro.f

ebarrier_edge_hydro.f: ../src/ebarrier_edge_hydro.F
	fpp $(DA) ../src/ebarrier_edge_hydro.F   ../work/ebarrier_edge_hydro.f

internal_edge_hydro.f: ../src/internal_edge_hydro.F
	fpp $(DA) ../src/internal_edge_hydro.F   ../work/internal_edge_hydro.f

land_edge_hydro.f: ../src/land_edge_hydro.F
	fpp $(DA) ../src/land_edge_hydro.F   ../work/land_edge_hydro.f

LDG_hydro.f: ../src/LDG_hydro.F
	fpp $(DA) ../src/LDG_hydro.F   ../work/LDG_hydro.f

modal2nodal.f: ../src/modal2nodal.F
	fpp $(DA) ../src/modal2nodal.F   ../work/modal2nodal.f

ocean_edge_hydro.f: ../src/ocean_edge_hydro.F
	fpp $(DA)  ../src/ocean_edge_hydro.F   ../work/ocean_edge_hydro.f

ocean_edge_hydro_post.f: ../src/ocean_edge_hydro_post.F
	fpp $(DA)  ../src/ocean_edge_hydro_post.F   ../work/ocean_edge_hydro_post.f

orthogonal_basis_v1.f: ../src/orthogonal_basis_v1.F
	fpp $(DA)  ../src/orthogonal_basis_v1.F   ../work/orthogonal_basis_v1.f

p_enrichment.f: ../src/p_enrichment.F    
	fpp $(DA)  ../src/p_enrichment.F ../work/p_enrichment.f

prep_dg.f: ../src/prep_dg.F    
	fpp $(DA)  ../src/prep_dg.F ../work/prep_dg.f

quadrature_v1.f: ../src/quadrature_v1.F
	fpp $(DA)  ../src/quadrature_v1.F   ../work/quadrature_v1.f

quad_rules_general.f: ../src/quad_rules_general.F
	fpp $(DA)  ../src/quad_rules_general.F   ../work/quad_rules_general.f

radiation_edge_hydro.f: ../src/radiation_edge_hydro.F
	fpp $(DA) ../src/radiation_edge_hydro.F   ../work/radiation_edge_hydro.f

read_fort_dg.f: ../src/read_fort_dg.F
	fpp $(DA)  ../src/read_fort_dg.F   ../work/read_fort_dg.f

read_input.f: ../src/read_input.F
	fpp $(DA) ../src/read_input.F   ../work/read_input.f

rhs_dg_hydro.f: ../src/rhs_dg_hydro.F
	fpp $(DA)  ../src/rhs_dg_hydro.F   ../work/rhs_dg_hydro.f

sta_basis.f: ../src/sta_basis.F
	fpp $(DA)  ../src/sta_basis.F   ../work/sta_basis.f

global_3dvs.f: ../src/global_3dvs.F
	fpp $(DA)  ../src/global_3dvs.F   ../work/global_3dvs.f

diff45_41.f: ../src/diff45_41.F
	fpp $(DA)  ../src/diff45_41.F   ../work/diff45_41.f

sizes.f: ../src/sizes.F   
	fpp $(DA)  ../src/sizes.F  ../work/sizes.f

tidal_potential.f: ../src/tidal_potential.F
      fpp $(DA) ../src/tidal_potential.F   ../work/tidal_potential.f

wind.f: ../src/wind.F   
	fpp $(DA)  ../src/wind.F ../work/wind.f

write_results.f: ../src/write_results.F
	fpp $(DA)  ../src/write_results.F   ../work/write_results.f

slopelimiter.f: ../src/slopelimiter.F
	fpp $(DA)  ../src/slopelimiter.F   ../work/slopelimiter.f

prep_slopelim.f: ../src/prep_slopelim.F
	fpp $(DA)  ../src/prep_slopelim.F   ../work/prep_slopelim.f

numerical_flux.f: ../src/numerical_flux.F
	fpp $(DA)  ../src/numerical_flux.F   ../work/numerical_flux.f

DG_timestep.f: ../src/DG_timestep.F
	fpp $(DA)  ../src/DG_timestep.F   ../work/DG_timestep.f

wetdry.f: ../src/wetdry.F
	fpp $(DA)  ../src/wetdry.F   ../work/wetdry.f

write_output.f: ../src/write_output.F
	fpp $(DA)  ../src/write_output.F   ../work/write_output.f

met_forcing.f: ../src/met_forcing.F
	fpp $(DA)  ../src/met_forcing.F   ../work/met_forcing.f
owiwind.f: ../src/owiwind.F
	fpp $(DA)  ../src/owiwind.F   ../work/owiwind.f

nodalattr.f: ../src/nodalattr.F
	fpp $(DA)  ../src/nodalattr.F   ../work/nodalattr.f

####################### Processed File Dependencies ##########################
#
#  this splits the proprocessed code into "modules" and "mains" and assigns
#  a general variable name to both
#

DGSWEM_MODS   =  sizes.f global.f harm.f wind.f dg.f diff45_41.f global_3dvs.f owiwind.f nodalattr.f read_fort_dg.f
DGSWEM_MANS   =  dgswem.f calc_normal.f create_edge_data.f cstart.f detector.f DG_hydro_timestep.f edge_int_hydro.f flow_edge_hydro.f hstart.f ibarrier_edge_hydro.f internal_edge_hydro.f land_edge_hydro.f ocean_edge_hydro.f orthogonal_basis_v1.f p_enrichment.f prep_DG.f quadrature_v1.f quad_rules_general.f radiation_edge_hydro.f read_input.f rhs_dg_hydro.f tidal_potential.f wind.f write_results.f slopelimiter.f prep_slopelim.f sta_basis.f LDG_hydro.f wetdry.f DG_timestep.f numerical_flux.f ocean_edge_hydro_post.f ebarrier_edge_hydro.f met_forcing.f write_output.f



########################## Executable Targets ###############################
#
#  this associates a target with each variable name
#

preprocessed_mods:  $(DGSWEM_MODS)
preprocessed_mans:  $(DGSWEM_MANS)


####################### Object File Dependencies ############################
#
#  this assigns a general variable name to the modules and mains object files
#

DGSWEM_MODS_OBJS   =  sizes.obj global.obj harm.obj wind.obj dg.obj diff45_41.obj global_3dvs.obj owiwind.obj nodalattr.obj read_fort_dg.obj
DGSWEM_MANS_OBJS   =  dgswem.obj calc_normal.obj create_edge_data.obj cstart.obj detector.obj DG_hydro_timestep.obj edge_int_hydro.obj flow_edge_hydro.obj hstart.obj ibarrier_edge_hydro.obj internal_edge_hydro.obj land_edge_hydro.obj ocean_edge_hydro.obj orthogonal_basis_v1.obj p_enrichment.obj prep_DG.obj quadrature_v1.obj quad_rules_general.obj radiation_edge_hydro.obj read_input.obj rhs_dg_hydro.obj tidal_potential.obj write_results.obj slopelimiter.obj prep_slopelim.obj sta_basis.obj LDG_hydro.obj wetdry.obj DG_timestep.obj numerical_flux.obj ocean_edge_hydro_post.obj ebarrier_edge_hydro.obj met_forcing.obj write_output.obj



######################## compilation rules #################################

$(DGSWEM_MODS_OBJS):  preprocessed_mods
	$(FC) -c $(FFLAGS)  $*.f

$(DGSWEM_MANS_OBJS):  preprocessed_mans
	$(FC) -c  $(FFLAGS) $*.f

########################## Executable Targets ###############################

dgswem: $(DGSWEM_MODS_OBJS) $(DGSWEM_MANS_OBJS)
	$(FC) /exe:dgswem  $(DGSWEM_MODS_OBJS) $(DGSWEM_MANS_OBJS)

########################## Misc Commands ####################################

clean:
	del *.f  *.mod  *.obj