# Makefile created by mkmf $Id: mkmf,v 18.0 2010/03/02 23:26:08 fms Exp $ 

SRCROOT = /Users/eweng/lm3ppa/src/



include /Users/eweng/lm3ppa/bin/mkmf.template.Mac


.DEFAULT:
	-echo $@ does not exist.
all: fms.x
astronomy.o: $(SRCROOT)shared/astronomy/astronomy.F90 utilities.o time_manager.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/astronomy/astronomy.F90
atmos_model.o: $(SRCROOT)atmos_solo_land/atmos_model.F90 utilities.o time_manager.o get_atmos_data.o astronomy.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_solo_land/atmos_model.F90
cana_tile.o: $(SRCROOT)land_lad2/canopy_air/cana_tile.F90 land_tile_diag_sel.o constants.o land_constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)land_lad2/canopy_air/cana_tile.F90
canopy_air.o: $(SRCROOT)land_lad2/canopy_air/canopy_air.F90 /Users/eweng/lm3ppa/src/land_lad2/canopy_air/../shared/debug.inc utilities.o time_manager.o constants.o sphum.o nf_utils.o land_constants.o cana_tile.o land_tile.o land_data.o land_tile_io.o land_debug.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/Users/eweng/lm3ppa/src/land_lad2/canopy_air/../shared	$(SRCROOT)land_lad2/canopy_air/canopy_air.F90
constants.o: $(SRCROOT)shared/constants/constants.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/constants/constants.F90
coupler_main.o: $(SRCROOT)coupler/coupler_main.F90 utilities.o constants.o time_manager.o atmos_model.o land_model.o flux_exchange.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)coupler/coupler_main.F90
flux_exchange.o: $(SRCROOT)coupler/flux_exchange.F90 utilities.o atmos_model.o land_model.o surface_flux.o monin_obukhov.o time_manager.o sat_vapor_pres.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)coupler/flux_exchange.F90
get_atmos_data.o: $(SRCROOT)atmos_solo_land/get_atmos_data.F90 utilities.o time_manager.o get_cal_time.o time_interp.o astronomy.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_solo_land/get_atmos_data.F90
get_cal_time.o: $(SRCROOT)shared/time_manager/get_cal_time.F90 utilities.o time_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/time_manager/get_cal_time.F90
glac_tile.o: $(SRCROOT)land_lad2/glacier/glac_tile.F90 /Users/eweng/lm3ppa/src/shared/include/fms_platform.h utilities.o constants.o land_constants.o land_tile_diag_sel.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/Users/eweng/lm3ppa/src/shared/include	$(SRCROOT)land_lad2/glacier/glac_tile.F90
glacier.o: $(SRCROOT)land_lad2/glacier/glacier.F90 utilities.o time_manager.o constants.o glac_tile.o land_constants.o land_tile.o land_tile_diag.o land_data.o land_io.o land_tile_io.o nf_utils.o land_debug.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)land_lad2/glacier/glacier.F90
lake.o: $(SRCROOT)land_lad2/lake/lake.F90 utilities.o time_manager.o constants.o land_constants.o lake_tile.o land_tile.o land_tile_diag.o land_data.o land_tile_io.o nf_utils.o land_debug.o land_utils.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)land_lad2/lake/lake.F90
lake_tile.o: $(SRCROOT)land_lad2/lake/lake_tile.F90 /Users/eweng/lm3ppa/src/shared/include/fms_platform.h utilities.o constants.o land_constants.o land_tile_diag_sel.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/Users/eweng/lm3ppa/src/shared/include	$(SRCROOT)land_lad2/lake/lake_tile.F90
land_constants.o: $(SRCROOT)land_lad2/land_constants.F90 constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)land_lad2/land_constants.F90
land_data.o: $(SRCROOT)land_lad2/land_data.F90 utilities.o constants.o time_manager.o land_tile.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)land_lad2/land_data.F90
land_debug.o: $(SRCROOT)land_lad2/shared/land_debug.F90 utilities.o time_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)land_lad2/shared/land_debug.F90
land_io.o: $(SRCROOT)land_lad2/shared/land_io.F90 utilities.o constants.o nf_utils.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)land_lad2/shared/land_io.F90
land_model.o: $(SRCROOT)land_lad2/land_model.F90 /Users/eweng/lm3ppa/src/land_lad2/shared/debug.inc utilities.o time_manager.o constants.o astronomy.o sphum.o land_constants.o glacier.o lake.o soil.o snow.o vegetation.o cana_tile.o canopy_air.o river.o topo_rough.o soil_tile.o vegn_tile.o lake_tile.o glac_tile.o snow_tile.o land_numerics.o land_tile.o land_data.o nf_utils.o land_utils.o land_tile_io.o land_tile_diag.o land_debug.o vegn_static_override.o transitions.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/Users/eweng/lm3ppa/src/land_lad2/shared	$(SRCROOT)land_lad2/land_model.F90
land_numerics.o: $(SRCROOT)land_lad2/shared/land_numerics.F90 utilities.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)land_lad2/shared/land_numerics.F90
land_tile.o: $(SRCROOT)land_lad2/land_tile.F90 utilities.o land_constants.o glac_tile.o lake_tile.o soil_tile.o cana_tile.o vegn_tile.o snow_tile.o land_tile_diag_sel.o land_tile_diag_buff.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)land_lad2/land_tile.F90
land_tile_diag.o: $(SRCROOT)land_lad2/shared/land_tile_diag.F90 utilities.o time_manager.o land_tile_diag_sel.o land_tile.o land_data.o land_tile_diag_buff.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)land_lad2/shared/land_tile_diag.F90
land_tile_diag_buff.o: $(SRCROOT)land_lad2/shared/land_tile_diag_buff.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)land_lad2/shared/land_tile_diag_buff.F90
land_tile_diag_sel.o: $(SRCROOT)land_lad2/shared/land_tile_diag_sel.F90 utilities.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)land_lad2/shared/land_tile_diag_sel.F90
land_tile_io.o: $(SRCROOT)land_lad2/shared/land_tile_io.F90 utilities.o time_manager.o nf_utils.o land_io.o land_tile.o land_data.o land_utils.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)land_lad2/shared/land_tile_io.F90
land_utils.o: $(SRCROOT)land_lad2/shared/land_utils.F90 land_tile.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)land_lad2/shared/land_utils.F90
monin_obukhov.o: $(SRCROOT)atmos_param/monin_obukhov/monin_obukhov.F90 utilities.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_param/monin_obukhov/monin_obukhov.F90
monin_obukhov_kernel.o: $(SRCROOT)atmos_param/monin_obukhov/monin_obukhov_kernel.F90 /Users/eweng/lm3ppa/src/shared/include/fms_platform.h /Users/eweng/lm3ppa/src/atmos_param/monin_obukhov/monin_obukhov_interfaces.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/Users/eweng/lm3ppa/src/shared/include -I/Users/eweng/lm3ppa/src/atmos_param/monin_obukhov	$(SRCROOT)atmos_param/monin_obukhov/monin_obukhov_kernel.F90
nf_utils.o: $(SRCROOT)land_lad2/shared/nf_utils/nf_utils.F90 nfu.o nfc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)land_lad2/shared/nf_utils/nf_utils.F90
nfc.o: $(SRCROOT)land_lad2/shared/nf_utils/nfc.F90 /Users/eweng/lm3ppa/src/land_lad2/shared/nf_utils/getput_compressed.inc nfu.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/Users/eweng/lm3ppa/src/land_lad2/shared/nf_utils	$(SRCROOT)land_lad2/shared/nf_utils/nfc.F90
nfu.o: $(SRCROOT)land_lad2/shared/nf_utils/nfu.F90 /Users/eweng/lm3ppa/src/land_lad2/shared/nf_utils/getput.inc
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/Users/eweng/lm3ppa/src/land_lad2/shared/nf_utils	$(SRCROOT)land_lad2/shared/nf_utils/nfu.F90
platform.o: $(SRCROOT)shared/platform/platform.F90 /Users/eweng/lm3ppa/src/shared/include/fms_platform.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/Users/eweng/lm3ppa/src/shared/include	$(SRCROOT)shared/platform/platform.F90
river.o: $(SRCROOT)land_lad2/river/river.F90 utilities.o time_manager.o river_type.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)land_lad2/river/river.F90
river_type.o: $(SRCROOT)land_lad2/river/river_type.F90 time_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)land_lad2/river/river_type.F90
sat_vapor_pres.o: $(SRCROOT)shared/sat_vapor_pres/sat_vapor_pres.F90 utilities.o constants.o sat_vapor_pres_k.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/sat_vapor_pres/sat_vapor_pres.F90
sat_vapor_pres_k.o: $(SRCROOT)shared/sat_vapor_pres/sat_vapor_pres_k.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/sat_vapor_pres/sat_vapor_pres_k.F90
snow.o: $(SRCROOT)land_lad2/snow/snow.F90 /Users/eweng/lm3ppa/src/land_lad2/canopy_air/../shared/debug.inc utilities.o time_manager.o constants.o land_constants.o snow_tile.o land_tile.o land_data.o land_tile_io.o nf_utils.o land_debug.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/Users/eweng/lm3ppa/src/land_lad2/canopy_air/../shared	$(SRCROOT)land_lad2/snow/snow.F90
snow_tile.o: $(SRCROOT)land_lad2/snow/snow_tile.F90 /Users/eweng/lm3ppa/src/shared/include/fms_platform.h utilities.o constants.o land_constants.o land_tile_diag_sel.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/Users/eweng/lm3ppa/src/shared/include	$(SRCROOT)land_lad2/snow/snow_tile.F90
soil.o: $(SRCROOT)land_lad2/soil/soil.F90 /Users/eweng/lm3ppa/src/land_lad2/canopy_air/../shared/debug.inc utilities.o time_manager.o constants.o land_constants.o soil_tile.o land_tile.o land_utils.o land_tile_diag.o land_data.o land_tile_io.o nf_utils.o vegn_tile.o vegn_cohort.o land_debug.o uptake.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/Users/eweng/lm3ppa/src/land_lad2/canopy_air/../shared	$(SRCROOT)land_lad2/soil/soil.F90
soil_tile.o: $(SRCROOT)land_lad2/soil/soil_tile.F90 /Users/eweng/lm3ppa/src/shared/include/fms_platform.h utilities.o constants.o land_constants.o land_io.o land_tile_diag_sel.o nf_utils.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/Users/eweng/lm3ppa/src/shared/include	$(SRCROOT)land_lad2/soil/soil_tile.F90
sphum.o: $(SRCROOT)land_lad2/shared/sphum.F90 constants.o sat_vapor_pres.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)land_lad2/shared/sphum.F90
surface_flux.o: $(SRCROOT)coupler/surface_flux.F90 utilities.o monin_obukhov.o sat_vapor_pres.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)coupler/surface_flux.F90
table_printer.o: $(SRCROOT)land_lad2/shared/table_printer.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)land_lad2/shared/table_printer.F90
time_interp.o: $(SRCROOT)shared/time_interp/time_interp.F90 utilities.o time_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/time_interp/time_interp.F90
time_manager.o: $(SRCROOT)shared/time_manager/time_manager.F90 constants.o utilities.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/time_manager/time_manager.F90
topo_rough.o: $(SRCROOT)land_lad2/topo_rough/topo_rough.F90 utilities.o time_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)land_lad2/topo_rough/topo_rough.F90
transitions.o: $(SRCROOT)land_lad2/transitions/transitions.F90 /Users/eweng/lm3ppa/src/shared/include/fms_platform.h /Users/eweng/lm3ppa/src/land_lad2/canopy_air/../shared/debug.inc utilities.o constants.o time_manager.o get_cal_time.o time_interp.o nfu.o vegn_data.o cana_tile.o snow_tile.o vegn_tile.o soil_tile.o land_tile.o land_tile_io.o land_data.o vegn_harvesting.o land_debug.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/Users/eweng/lm3ppa/src/shared/include -I/Users/eweng/lm3ppa/src/land_lad2/canopy_air/../shared	$(SRCROOT)land_lad2/transitions/transitions.F90
uptake.o: $(SRCROOT)land_lad2/soil/uptake.F90 /Users/eweng/lm3ppa/src/land_lad2/canopy_air/../shared/debug.inc utilities.o constants.o soil_tile.o land_debug.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/Users/eweng/lm3ppa/src/land_lad2/canopy_air/../shared	$(SRCROOT)land_lad2/soil/uptake.F90
utilities.o: $(SRCROOT)shared/utilities/utilities.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/utilities/utilities.F90
vegetation.o: $(SRCROOT)land_lad2/vegetation/vegetation.F90 /Users/eweng/lm3ppa/src/land_lad2/canopy_air/../shared/debug.inc utilities.o time_manager.o constants.o sphum.o nf_utils.o vegn_tile.o soil_tile.o land_constants.o land_tile.o land_tile_diag.o land_data.o land_tile_io.o vegn_data.o vegn_cohort.o canopy_air.o soil.o vegn_cohort_io.o land_debug.o vegn_radiation.o vegn_photosynthesis.o vegn_static_override.o vegn_dynamics.o vegn_disturbance.o vegn_harvesting.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/Users/eweng/lm3ppa/src/land_lad2/canopy_air/../shared	$(SRCROOT)land_lad2/vegetation/vegetation.F90
vegn_cohort.o: $(SRCROOT)land_lad2/vegetation/vegn_cohort.F90 constants.o land_constants.o vegn_data.o soil_tile.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)land_lad2/vegetation/vegn_cohort.F90
vegn_cohort_io.o: $(SRCROOT)land_lad2/vegetation/vegn_cohort_io.F90 /Users/eweng/lm3ppa/src/land_lad2/vegetation/vegn_cohort_io.inc utilities.o nf_utils.o land_io.o land_tile.o land_tile_io.o vegn_cohort.o land_data.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/Users/eweng/lm3ppa/src/land_lad2/vegetation	$(SRCROOT)land_lad2/vegetation/vegn_cohort_io.F90
vegn_data.o: $(SRCROOT)land_lad2/vegetation/vegn_data.F90 utilities.o constants.o table_printer.o land_constants.o land_tile_diag_sel.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)land_lad2/vegetation/vegn_data.F90
vegn_disturbance.o: $(SRCROOT)land_lad2/vegetation/vegn_disturbance.F90 /Users/eweng/lm3ppa/src/land_lad2/canopy_air/../shared/debug.inc utilities.o land_constants.o land_debug.o vegn_data.o vegn_tile.o vegn_cohort.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/Users/eweng/lm3ppa/src/land_lad2/canopy_air/../shared	$(SRCROOT)land_lad2/vegetation/vegn_disturbance.F90
vegn_dynamics.o: $(SRCROOT)land_lad2/vegetation/vegn_dynamics.F90 /Users/eweng/lm3ppa/src/land_lad2/canopy_air/../shared/debug.inc utilities.o time_manager.o constants.o land_constants.o land_tile_diag.o land_data.o vegn_data.o vegn_tile.o soil_tile.o vegn_cohort.o land_debug.o land_numerics.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/Users/eweng/lm3ppa/src/land_lad2/canopy_air/../shared	$(SRCROOT)land_lad2/vegetation/vegn_dynamics.F90
vegn_harvesting.o: $(SRCROOT)land_lad2/vegetation/vegn_harvesting.F90 utilities.o vegn_data.o vegn_tile.o vegn_cohort.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)land_lad2/vegetation/vegn_harvesting.F90
vegn_photosynthesis.o: $(SRCROOT)land_lad2/vegetation/vegn_photosynthesis.F90 /Users/eweng/lm3ppa/src/land_lad2/canopy_air/../shared/debug.inc utilities.o constants.o sphum.o land_constants.o land_debug.o vegn_data.o vegn_cohort.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/Users/eweng/lm3ppa/src/land_lad2/canopy_air/../shared	$(SRCROOT)land_lad2/vegetation/vegn_photosynthesis.F90
vegn_radiation.o: $(SRCROOT)land_lad2/vegetation/vegn_radiation.F90 utilities.o constants.o land_constants.o vegn_data.o vegn_tile.o vegn_cohort.o snow.o land_debug.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)land_lad2/vegetation/vegn_radiation.F90
vegn_static_override.o: $(SRCROOT)land_lad2/vegetation/vegn_static_override.F90 /Users/eweng/lm3ppa/src/land_lad2/vegetation/read_remap_cohort_data.inc utilities.o constants.o time_manager.o get_cal_time.o time_interp.o nf_utils.o land_data.o land_io.o land_numerics.o land_tile_io.o land_tile.o vegn_cohort.o vegn_cohort_io.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/Users/eweng/lm3ppa/src/land_lad2/vegetation	$(SRCROOT)land_lad2/vegetation/vegn_static_override.F90
vegn_tile.o: $(SRCROOT)land_lad2/vegetation/vegn_tile.F90 utilities.o constants.o land_constants.o land_tile_diag_sel.o vegn_data.o vegn_cohort.o soil_tile.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)land_lad2/vegetation/vegn_tile.F90
./glacier.F90: $(SRCROOT)land_lad2/glacier/glacier.F90
	cp $(SRCROOT)land_lad2/glacier/glacier.F90 .
./soil.F90: $(SRCROOT)land_lad2/soil/soil.F90
	cp $(SRCROOT)land_lad2/soil/soil.F90 .
./time_interp.F90: $(SRCROOT)shared/time_interp/time_interp.F90
	cp $(SRCROOT)shared/time_interp/time_interp.F90 .
./vegn_data.F90: $(SRCROOT)land_lad2/vegetation/vegn_data.F90
	cp $(SRCROOT)land_lad2/vegetation/vegn_data.F90 .
./vegetation.F90: $(SRCROOT)land_lad2/vegetation/vegetation.F90
	cp $(SRCROOT)land_lad2/vegetation/vegetation.F90 .
./uptake.F90: $(SRCROOT)land_lad2/soil/uptake.F90
	cp $(SRCROOT)land_lad2/soil/uptake.F90 .
./river.F90: $(SRCROOT)land_lad2/river/river.F90
	cp $(SRCROOT)land_lad2/river/river.F90 .
./land_tile_io.F90: $(SRCROOT)land_lad2/shared/land_tile_io.F90
	cp $(SRCROOT)land_lad2/shared/land_tile_io.F90 .
./vegn_photosynthesis.F90: $(SRCROOT)land_lad2/vegetation/vegn_photosynthesis.F90
	cp $(SRCROOT)land_lad2/vegetation/vegn_photosynthesis.F90 .
./land_constants.F90: $(SRCROOT)land_lad2/land_constants.F90
	cp $(SRCROOT)land_lad2/land_constants.F90 .
./land_io.F90: $(SRCROOT)land_lad2/shared/land_io.F90
	cp $(SRCROOT)land_lad2/shared/land_io.F90 .
./river_type.F90: $(SRCROOT)land_lad2/river/river_type.F90
	cp $(SRCROOT)land_lad2/river/river_type.F90 .
./vegn_cohort_io.inc: /Users/eweng/lm3ppa/src/land_lad2/vegetation/vegn_cohort_io.inc
	cp /Users/eweng/lm3ppa/src/land_lad2/vegetation/vegn_cohort_io.inc .
./vegn_cohort.F90: $(SRCROOT)land_lad2/vegetation/vegn_cohort.F90
	cp $(SRCROOT)land_lad2/vegetation/vegn_cohort.F90 .
./vegn_harvesting.F90: $(SRCROOT)land_lad2/vegetation/vegn_harvesting.F90
	cp $(SRCROOT)land_lad2/vegetation/vegn_harvesting.F90 .
./vegn_disturbance.F90: $(SRCROOT)land_lad2/vegetation/vegn_disturbance.F90
	cp $(SRCROOT)land_lad2/vegetation/vegn_disturbance.F90 .
./getput.inc: /Users/eweng/lm3ppa/src/land_lad2/shared/nf_utils/getput.inc
	cp /Users/eweng/lm3ppa/src/land_lad2/shared/nf_utils/getput.inc .
./table_printer.F90: $(SRCROOT)land_lad2/shared/table_printer.F90
	cp $(SRCROOT)land_lad2/shared/table_printer.F90 .
./cana_tile.F90: $(SRCROOT)land_lad2/canopy_air/cana_tile.F90
	cp $(SRCROOT)land_lad2/canopy_air/cana_tile.F90 .
./monin_obukhov_interfaces.h: /Users/eweng/lm3ppa/src/atmos_param/monin_obukhov/monin_obukhov_interfaces.h
	cp /Users/eweng/lm3ppa/src/atmos_param/monin_obukhov/monin_obukhov_interfaces.h .
./constants.F90: $(SRCROOT)shared/constants/constants.F90
	cp $(SRCROOT)shared/constants/constants.F90 .
./get_atmos_data.F90: $(SRCROOT)atmos_solo_land/get_atmos_data.F90
	cp $(SRCROOT)atmos_solo_land/get_atmos_data.F90 .
./utilities.F90: $(SRCROOT)shared/utilities/utilities.F90
	cp $(SRCROOT)shared/utilities/utilities.F90 .
./land_debug.F90: $(SRCROOT)land_lad2/shared/land_debug.F90
	cp $(SRCROOT)land_lad2/shared/land_debug.F90 .
./snow_tile.F90: $(SRCROOT)land_lad2/snow/snow_tile.F90
	cp $(SRCROOT)land_lad2/snow/snow_tile.F90 .
./nf_utils.F90: $(SRCROOT)land_lad2/shared/nf_utils/nf_utils.F90
	cp $(SRCROOT)land_lad2/shared/nf_utils/nf_utils.F90 .
./lake_tile.F90: $(SRCROOT)land_lad2/lake/lake_tile.F90
	cp $(SRCROOT)land_lad2/lake/lake_tile.F90 .
./platform.F90: $(SRCROOT)shared/platform/platform.F90
	cp $(SRCROOT)shared/platform/platform.F90 .
./land_tile_diag.F90: $(SRCROOT)land_lad2/shared/land_tile_diag.F90
	cp $(SRCROOT)land_lad2/shared/land_tile_diag.F90 .
./snow.F90: $(SRCROOT)land_lad2/snow/snow.F90
	cp $(SRCROOT)land_lad2/snow/snow.F90 .
./land_utils.F90: $(SRCROOT)land_lad2/shared/land_utils.F90
	cp $(SRCROOT)land_lad2/shared/land_utils.F90 .
./flux_exchange.F90: $(SRCROOT)coupler/flux_exchange.F90
	cp $(SRCROOT)coupler/flux_exchange.F90 .
./land_tile_diag_sel.F90: $(SRCROOT)land_lad2/shared/land_tile_diag_sel.F90
	cp $(SRCROOT)land_lad2/shared/land_tile_diag_sel.F90 .
./land_tile.F90: $(SRCROOT)land_lad2/land_tile.F90
	cp $(SRCROOT)land_lad2/land_tile.F90 .
./sphum.F90: $(SRCROOT)land_lad2/shared/sphum.F90
	cp $(SRCROOT)land_lad2/shared/sphum.F90 .
./debug.inc: /Users/eweng/lm3ppa/src/land_lad2/shared/debug.inc
	cp /Users/eweng/lm3ppa/src/land_lad2/shared/debug.inc .
./vegn_dynamics.F90: $(SRCROOT)land_lad2/vegetation/vegn_dynamics.F90
	cp $(SRCROOT)land_lad2/vegetation/vegn_dynamics.F90 .
./land_model.F90: $(SRCROOT)land_lad2/land_model.F90
	cp $(SRCROOT)land_lad2/land_model.F90 .
./lake.F90: $(SRCROOT)land_lad2/lake/lake.F90
	cp $(SRCROOT)land_lad2/lake/lake.F90 .
./nfu.F90: $(SRCROOT)land_lad2/shared/nf_utils/nfu.F90
	cp $(SRCROOT)land_lad2/shared/nf_utils/nfu.F90 .
./sat_vapor_pres.F90: $(SRCROOT)shared/sat_vapor_pres/sat_vapor_pres.F90
	cp $(SRCROOT)shared/sat_vapor_pres/sat_vapor_pres.F90 .
./sat_vapor_pres_k.F90: $(SRCROOT)shared/sat_vapor_pres/sat_vapor_pres_k.F90
	cp $(SRCROOT)shared/sat_vapor_pres/sat_vapor_pres_k.F90 .
./fms_platform.h: /Users/eweng/lm3ppa/src/shared/include/fms_platform.h
	cp /Users/eweng/lm3ppa/src/shared/include/fms_platform.h .
./nfc.F90: $(SRCROOT)land_lad2/shared/nf_utils/nfc.F90
	cp $(SRCROOT)land_lad2/shared/nf_utils/nfc.F90 .
./canopy_air.F90: $(SRCROOT)land_lad2/canopy_air/canopy_air.F90
	cp $(SRCROOT)land_lad2/canopy_air/canopy_air.F90 .
./glac_tile.F90: $(SRCROOT)land_lad2/glacier/glac_tile.F90
	cp $(SRCROOT)land_lad2/glacier/glac_tile.F90 .
./vegn_static_override.F90: $(SRCROOT)land_lad2/vegetation/vegn_static_override.F90
	cp $(SRCROOT)land_lad2/vegetation/vegn_static_override.F90 .
./vegn_tile.F90: $(SRCROOT)land_lad2/vegetation/vegn_tile.F90
	cp $(SRCROOT)land_lad2/vegetation/vegn_tile.F90 .
./get_cal_time.F90: $(SRCROOT)shared/time_manager/get_cal_time.F90
	cp $(SRCROOT)shared/time_manager/get_cal_time.F90 .
./land_data.F90: $(SRCROOT)land_lad2/land_data.F90
	cp $(SRCROOT)land_lad2/land_data.F90 .
./read_remap_cohort_data.inc: /Users/eweng/lm3ppa/src/land_lad2/vegetation/read_remap_cohort_data.inc
	cp /Users/eweng/lm3ppa/src/land_lad2/vegetation/read_remap_cohort_data.inc .
./vegn_radiation.F90: $(SRCROOT)land_lad2/vegetation/vegn_radiation.F90
	cp $(SRCROOT)land_lad2/vegetation/vegn_radiation.F90 .
./debug.inc: /Users/eweng/lm3ppa/src/land_lad2/canopy_air/../shared/debug.inc
	cp /Users/eweng/lm3ppa/src/land_lad2/canopy_air/../shared/debug.inc .
./soil_tile.F90: $(SRCROOT)land_lad2/soil/soil_tile.F90
	cp $(SRCROOT)land_lad2/soil/soil_tile.F90 .
./astronomy.F90: $(SRCROOT)shared/astronomy/astronomy.F90
	cp $(SRCROOT)shared/astronomy/astronomy.F90 .
./monin_obukhov.F90: $(SRCROOT)atmos_param/monin_obukhov/monin_obukhov.F90
	cp $(SRCROOT)atmos_param/monin_obukhov/monin_obukhov.F90 .
./land_numerics.F90: $(SRCROOT)land_lad2/shared/land_numerics.F90
	cp $(SRCROOT)land_lad2/shared/land_numerics.F90 .
./vegn_cohort_io.F90: $(SRCROOT)land_lad2/vegetation/vegn_cohort_io.F90
	cp $(SRCROOT)land_lad2/vegetation/vegn_cohort_io.F90 .
./transitions.F90: $(SRCROOT)land_lad2/transitions/transitions.F90
	cp $(SRCROOT)land_lad2/transitions/transitions.F90 .
./monin_obukhov_kernel.F90: $(SRCROOT)atmos_param/monin_obukhov/monin_obukhov_kernel.F90
	cp $(SRCROOT)atmos_param/monin_obukhov/monin_obukhov_kernel.F90 .
./time_manager.F90: $(SRCROOT)shared/time_manager/time_manager.F90
	cp $(SRCROOT)shared/time_manager/time_manager.F90 .
./coupler_main.F90: $(SRCROOT)coupler/coupler_main.F90
	cp $(SRCROOT)coupler/coupler_main.F90 .
./topo_rough.F90: $(SRCROOT)land_lad2/topo_rough/topo_rough.F90
	cp $(SRCROOT)land_lad2/topo_rough/topo_rough.F90 .
./surface_flux.F90: $(SRCROOT)coupler/surface_flux.F90
	cp $(SRCROOT)coupler/surface_flux.F90 .
./atmos_model.F90: $(SRCROOT)atmos_solo_land/atmos_model.F90
	cp $(SRCROOT)atmos_solo_land/atmos_model.F90 .
./getput_compressed.inc: /Users/eweng/lm3ppa/src/land_lad2/shared/nf_utils/getput_compressed.inc
	cp /Users/eweng/lm3ppa/src/land_lad2/shared/nf_utils/getput_compressed.inc .
./land_tile_diag_buff.F90: $(SRCROOT)land_lad2/shared/land_tile_diag_buff.F90
	cp $(SRCROOT)land_lad2/shared/land_tile_diag_buff.F90 .
SRC = $(SRCROOT)coupler/surface_flux.F90 $(SRCROOT)land_lad2/soil/soil.F90 $(SRCROOT)land_lad2/soil/uptake.F90 $(SRCROOT)land_lad2/vegetation/vegetation.F90 $(SRCROOT)atmos_param/monin_obukhov/monin_obukhov_kernel.F90 $(SRCROOT)land_lad2/vegetation/vegn_dynamics.F90 $(SRCROOT)land_lad2/shared/land_debug.F90 $(SRCROOT)land_lad2/soil/soil_tile.F90 $(SRCROOT)land_lad2/land_model.F90 $(SRCROOT)coupler/flux_exchange.F90 $(SRCROOT)land_lad2/shared/sphum.F90 $(SRCROOT)land_lad2/river/river.F90 $(SRCROOT)land_lad2/vegetation/vegn_cohort.F90 $(SRCROOT)shared/sat_vapor_pres/sat_vapor_pres_k.F90 $(SRCROOT)land_lad2/shared/land_tile_diag_buff.F90 $(SRCROOT)land_lad2/vegetation/vegn_radiation.F90 $(SRCROOT)atmos_solo_land/get_atmos_data.F90 $(SRCROOT)land_lad2/vegetation/vegn_harvesting.F90 $(SRCROOT)land_lad2/land_constants.F90 $(SRCROOT)shared/constants/constants.F90 $(SRCROOT)land_lad2/shared/nf_utils/nfu.F90 $(SRCROOT)land_lad2/glacier/glacier.F90 $(SRCROOT)shared/astronomy/astronomy.F90 $(SRCROOT)atmos_solo_land/atmos_model.F90 $(SRCROOT)shared/time_manager/time_manager.F90 $(SRCROOT)land_lad2/shared/table_printer.F90 $(SRCROOT)land_lad2/shared/land_tile_diag.F90 $(SRCROOT)land_lad2/shared/land_tile_io.F90 $(SRCROOT)land_lad2/shared/land_utils.F90 $(SRCROOT)land_lad2/topo_rough/topo_rough.F90 $(SRCROOT)land_lad2/shared/land_io.F90 $(SRCROOT)coupler/coupler_main.F90 $(SRCROOT)land_lad2/canopy_air/cana_tile.F90 $(SRCROOT)land_lad2/vegetation/vegn_static_override.F90 $(SRCROOT)land_lad2/lake/lake.F90 $(SRCROOT)land_lad2/transitions/transitions.F90 $(SRCROOT)shared/time_manager/get_cal_time.F90 $(SRCROOT)land_lad2/lake/lake_tile.F90 $(SRCROOT)land_lad2/vegetation/vegn_cohort_io.F90 $(SRCROOT)shared/platform/platform.F90 $(SRCROOT)land_lad2/shared/land_numerics.F90 $(SRCROOT)land_lad2/vegetation/vegn_data.F90 $(SRCROOT)land_lad2/land_tile.F90 $(SRCROOT)land_lad2/vegetation/vegn_photosynthesis.F90 $(SRCROOT)land_lad2/glacier/glac_tile.F90 $(SRCROOT)land_lad2/land_data.F90 $(SRCROOT)land_lad2/vegetation/vegn_tile.F90 $(SRCROOT)land_lad2/river/river_type.F90 $(SRCROOT)land_lad2/shared/nf_utils/nf_utils.F90 $(SRCROOT)land_lad2/snow/snow.F90 $(SRCROOT)shared/time_interp/time_interp.F90 $(SRCROOT)land_lad2/shared/land_tile_diag_sel.F90 $(SRCROOT)land_lad2/canopy_air/canopy_air.F90 $(SRCROOT)land_lad2/snow/snow_tile.F90 $(SRCROOT)shared/utilities/utilities.F90 $(SRCROOT)land_lad2/shared/nf_utils/nfc.F90 $(SRCROOT)shared/sat_vapor_pres/sat_vapor_pres.F90 $(SRCROOT)land_lad2/vegetation/vegn_disturbance.F90 $(SRCROOT)atmos_param/monin_obukhov/monin_obukhov.F90 /Users/eweng/lm3ppa/src/land_lad2/vegetation/vegn_cohort_io.inc /Users/eweng/lm3ppa/src/atmos_param/monin_obukhov/monin_obukhov_interfaces.h /Users/eweng/lm3ppa/src/land_lad2/canopy_air/../shared/debug.inc /Users/eweng/lm3ppa/src/land_lad2/shared/nf_utils/getput_compressed.inc /Users/eweng/lm3ppa/src/land_lad2/shared/debug.inc /Users/eweng/lm3ppa/src/land_lad2/shared/nf_utils/getput.inc /Users/eweng/lm3ppa/src/shared/include/fms_platform.h /Users/eweng/lm3ppa/src/land_lad2/vegetation/read_remap_cohort_data.inc
OBJ = surface_flux.o soil.o uptake.o vegetation.o monin_obukhov_kernel.o vegn_dynamics.o land_debug.o soil_tile.o land_model.o flux_exchange.o sphum.o river.o vegn_cohort.o sat_vapor_pres_k.o land_tile_diag_buff.o vegn_radiation.o get_atmos_data.o vegn_harvesting.o land_constants.o constants.o nfu.o glacier.o astronomy.o atmos_model.o time_manager.o table_printer.o land_tile_diag.o land_tile_io.o land_utils.o topo_rough.o land_io.o coupler_main.o cana_tile.o vegn_static_override.o lake.o transitions.o get_cal_time.o lake_tile.o vegn_cohort_io.o platform.o land_numerics.o vegn_data.o land_tile.o vegn_photosynthesis.o glac_tile.o land_data.o vegn_tile.o river_type.o nf_utils.o snow.o time_interp.o land_tile_diag_sel.o canopy_air.o snow_tile.o utilities.o nfc.o sat_vapor_pres.o vegn_disturbance.o monin_obukhov.o
OFF = $(SRCROOT)land_lad2/glacier/glacier.F90 $(SRCROOT)land_lad2/soil/soil.F90 $(SRCROOT)shared/time_interp/time_interp.F90 $(SRCROOT)land_lad2/vegetation/vegn_data.F90 $(SRCROOT)land_lad2/vegetation/vegetation.F90 $(SRCROOT)land_lad2/soil/uptake.F90 $(SRCROOT)land_lad2/river/river.F90 $(SRCROOT)land_lad2/shared/land_tile_io.F90 $(SRCROOT)land_lad2/vegetation/vegn_photosynthesis.F90 $(SRCROOT)land_lad2/land_constants.F90 $(SRCROOT)land_lad2/shared/land_io.F90 $(SRCROOT)land_lad2/river/river_type.F90 /Users/eweng/lm3ppa/src/land_lad2/vegetation/vegn_cohort_io.inc $(SRCROOT)land_lad2/vegetation/vegn_cohort.F90 $(SRCROOT)land_lad2/vegetation/vegn_harvesting.F90 $(SRCROOT)land_lad2/vegetation/vegn_disturbance.F90 /Users/eweng/lm3ppa/src/land_lad2/shared/nf_utils/getput.inc $(SRCROOT)land_lad2/shared/table_printer.F90 $(SRCROOT)land_lad2/canopy_air/cana_tile.F90 /Users/eweng/lm3ppa/src/atmos_param/monin_obukhov/monin_obukhov_interfaces.h $(SRCROOT)shared/constants/constants.F90 $(SRCROOT)atmos_solo_land/get_atmos_data.F90 $(SRCROOT)shared/utilities/utilities.F90 $(SRCROOT)land_lad2/shared/land_debug.F90 $(SRCROOT)land_lad2/snow/snow_tile.F90 $(SRCROOT)land_lad2/shared/nf_utils/nf_utils.F90 $(SRCROOT)land_lad2/lake/lake_tile.F90 $(SRCROOT)shared/platform/platform.F90 $(SRCROOT)land_lad2/shared/land_tile_diag.F90 $(SRCROOT)land_lad2/snow/snow.F90 $(SRCROOT)land_lad2/shared/land_utils.F90 $(SRCROOT)coupler/flux_exchange.F90 $(SRCROOT)land_lad2/shared/land_tile_diag_sel.F90 $(SRCROOT)land_lad2/land_tile.F90 $(SRCROOT)land_lad2/shared/sphum.F90 /Users/eweng/lm3ppa/src/land_lad2/shared/debug.inc $(SRCROOT)land_lad2/vegetation/vegn_dynamics.F90 $(SRCROOT)land_lad2/land_model.F90 $(SRCROOT)land_lad2/lake/lake.F90 $(SRCROOT)land_lad2/shared/nf_utils/nfu.F90 $(SRCROOT)shared/sat_vapor_pres/sat_vapor_pres.F90 $(SRCROOT)shared/sat_vapor_pres/sat_vapor_pres_k.F90 /Users/eweng/lm3ppa/src/shared/include/fms_platform.h $(SRCROOT)land_lad2/shared/nf_utils/nfc.F90 $(SRCROOT)land_lad2/canopy_air/canopy_air.F90 $(SRCROOT)land_lad2/glacier/glac_tile.F90 $(SRCROOT)land_lad2/vegetation/vegn_static_override.F90 $(SRCROOT)land_lad2/vegetation/vegn_tile.F90 $(SRCROOT)shared/time_manager/get_cal_time.F90 $(SRCROOT)land_lad2/land_data.F90 /Users/eweng/lm3ppa/src/land_lad2/vegetation/read_remap_cohort_data.inc $(SRCROOT)land_lad2/vegetation/vegn_radiation.F90 /Users/eweng/lm3ppa/src/land_lad2/canopy_air/../shared/debug.inc $(SRCROOT)land_lad2/soil/soil_tile.F90 $(SRCROOT)shared/astronomy/astronomy.F90 $(SRCROOT)atmos_param/monin_obukhov/monin_obukhov.F90 $(SRCROOT)land_lad2/shared/land_numerics.F90 $(SRCROOT)land_lad2/vegetation/vegn_cohort_io.F90 $(SRCROOT)land_lad2/transitions/transitions.F90 $(SRCROOT)atmos_param/monin_obukhov/monin_obukhov_kernel.F90 $(SRCROOT)shared/time_manager/time_manager.F90 $(SRCROOT)coupler/coupler_main.F90 $(SRCROOT)land_lad2/topo_rough/topo_rough.F90 $(SRCROOT)coupler/surface_flux.F90 $(SRCROOT)atmos_solo_land/atmos_model.F90 /Users/eweng/lm3ppa/src/land_lad2/shared/nf_utils/getput_compressed.inc $(SRCROOT)land_lad2/shared/land_tile_diag_buff.F90
clean: neat
	-rm -f .fms.x.cppdefs $(OBJ) fms.x
neat:
	-rm -f $(TMPFILES)
localize: $(OFF)
	cp $(OFF) .
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
fms.x: $(OBJ) 
	$(LD) $(OBJ) -o fms.x  $(LDFLAGS)
