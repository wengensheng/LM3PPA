module flux_exchange_mod
!-----------------------------------------------------------------------
!                   GNU General Public License                        !                                                                      
! This program is free software; you can redistribute it and/or modify it and  
! are expected to follow the terms of the GNU General Public License  
! as published by the Free Software Foundation; either version 2 of   
! the License, or (at your option) any later version.                 
!                                                                      
! MOM is distributed in the hope that it will be useful, but WITHOUT    
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  
! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public    
! License for more details.                                           
!                                                                      
! For the full text of the GNU General Public License,                
! write to: Free Software Foundation, Inc.,                           
!           675 Mass Ave, Cambridge, MA 02139, USA.                   
! or see:   http://www.gnu.org/licenses/gpl.html                      
!-----------------------------------------------------------------------
! <CONTACT EMAIL="Bruce.Wyman@noaa.gov"> Bruce Wyman </CONTACT>
! <CONTACT EMAIL="V.Balaji@noaa.gov"> V. Balaji </CONTACT>
! <CONTACT EMAIL="Sergey.Malyshev@noaa.gov"> Sergey Malyshev </CONTACT>

! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>

! <OVERVIEW>
!   The flux_exchange module provides interfaces to couple
!   the following component models: atmosphere and land.
! </OVERVIEW>

! <DESCRIPTION>
!  <PRE>
!  1.This version of flux_exchange_mod allows the definition of physically independent
!    grids for atmosphere and land.
!
!  2.Each component model must have a public defined data type containing specific 
!    boundary fields. A list of these quantities is located in the NOTES of this document. 
!
!  3.The surface flux of sensible heat and surface evaporation can be implicit functions
!    of surface temperature. As a consequence, the parts of the land model
!    that updates the surface temperature must be called on the atmospheric time step 
!
!  4.The surface fluxes of all other tracers and of momentum are assumed to be explicit
!    functions of all surface parameters 
!
!  5.While no explicit reference is made within this module to the implicit treatment 
!    of vertical diffusion in the atmosphere and land models, the 
!    module is designed to allow for simultaneous implicit time integration on both 
!    sides of the surface interface. 
!
!  6.Due to #5, the diffusion part of the land model must be called on the 
!    atmospheric time step.
  
!  </PRE>
! </DESCRIPTION>

  use utilities_mod,  only: logunit, get_unit, atmos_tracer_longnames, atmos_tracer_units, &
                            atmos_tracer_names, num_atmos_tracers, land_tracer_names, num_land_tracers, NO_TRACER, &
                            error_mesg, mpp_error, FATAL, NOTE, check_nml_error

!model_boundary_data_type contains all model fields at the boundary.
!model1_model2_boundary_type contains fields that model2 gets
!from model1, may also include fluxes. These are declared by
!flux_exchange_mod and have private components. All model fields in
!model_boundary_data_type may not be exchanged.
  use atmos_model_mod, only: atmos_data_type, land_ice_atmos_boundary_type
  use  land_model_mod, only: land_data_type, atmos_land_boundary_type

  use  surface_flux_mod, only: surface_flux
  use monin_obukhov_mod, only: mo_profile     

  use  time_manager_mod, only: time_type, get_date

  use sat_vapor_pres_mod, only: compute_qs

  use      constants_mod, only: rdgas, rvgas, cp_air, stefan, WTMAIR, HLV, HLF, Radius, PI, CP_OCEAN, &
                                WTMCO2, WTMC

  implicit none
  include 'netcdf.inc'
private

  character(len=48), parameter :: module_name = 'flux_exchange_mod'
  integer, parameter :: ind_u10=2, ind_psurf=3

  public :: flux_exchange_init,   &
     sfc_boundary_layer,   &
     flux_down_from_atmos, &
     flux_up_to_atmos

!-----------------------------------------------------------------------
  character(len=128) :: version = '$Id: flux_exchange.F90,v 1.1.2.1.2.19 2012/06/19 18:34:53 pjp Exp $'
  character(len=128) :: tagname = '$Name: no_fms_b_pjp $'
!-----------------------------------------------------------------------

integer :: n_xgrid_sfc
real, parameter :: area=2.42256098020317e-05

!-----------------------------------------------------------------------
!-------- namelist (for diagnostics) ------

character(len=4), parameter :: mod_name = 'flux'

  integer :: id_drag_moist,  id_drag_heat,  id_drag_mom,     &
     id_rough_moist, id_rough_heat, id_rough_mom, id_land_mask, &
     id_u_star, id_b_star, id_q_star, id_u_flux, id_v_flux,   &
     id_t_surf, id_t_flux, id_r_flux, id_q_flux, id_slp,      &
     id_t_atm,  id_u_atm,  id_v_atm,  id_wind,                &
     id_t_ref,  id_rh_ref, id_u_ref,  id_v_ref, id_wind_ref,  &
     id_del_h,  id_del_m,  id_del_q,  id_rough_scale,         &
     id_t_ca,   id_q_surf, id_q_atm, id_z_atm, id_p_atm, id_gust, &
     id_t_ref_land, id_rh_ref_land, id_u_ref_land, id_v_ref_land, &
     id_q_ref,  id_q_ref_land, id_q_flux_land, id_rh_ref_cmip

integer :: id_co2_atm_dvmr, id_co2_surf_dvmr

integer, allocatable :: id_tr_atm(:), id_tr_surf(:), id_tr_flux(:), id_tr_mol_flux(:)

logical :: first_static = .true.
logical :: do_init = .true.
integer :: remap_method = 1

real, parameter :: bound_tol = 1e-7

real, parameter :: d622 = rdgas/rvgas
real, parameter :: d378 = 1.0-d622
character(*), parameter ::      diag_format='(i4,"-",i2.2,"-",i2.2,1x,i2.2,":",i2.2,":",i2.2,", ",1pe24.16)'
character(*), parameter :: tile_diag_format='(i4,"-",i2.2,"-",i2.2,1x,i2.2,":",i2.2,":",i2.2," tile=",i2.2,", ",1pe24.16)'
integer, parameter :: fm_field_name_len=48, fm_string_len=128

type, public :: coupler_1d_values_type
 character(len=fm_field_name_len) :: name = ' '
 real, pointer, dimension(:) :: values => NULL()
 logical :: mean = .true.
 logical :: override = .false.
 integer :: id_diag = 0
 character(len=fm_string_len) :: long_name = ' '
 character(len=fm_string_len) :: units = ' '
end type coupler_1d_values_type

type, public :: coupler_1d_field_type
 character(len=fm_field_name_len) :: name = ' '
 integer :: num_fields = 0
 type(coupler_1d_values_type), pointer, dimension(:) :: field => NULL()
 character(len=fm_string_len) :: flux_type = ' '
 character(len=fm_string_len) :: implementation = ' '
 real, pointer, dimension(:) :: param => NULL()
 logical, pointer, dimension(:) :: flag => NULL()
 integer :: atm_tr_index = 0
 character(len=fm_string_len) :: ice_restart_file = ' '
 character(len=fm_string_len) :: ocean_restart_file = ' '
 logical :: use_atm_pressure
 logical :: use_10m_wind_speed
 logical :: pass_through_ice
 real :: mol_wt = 0.0
end type coupler_1d_field_type

type, public :: coupler_1d_bc_type
 integer :: num_bcs = 0
 type(coupler_1d_field_type), pointer, dimension(:) :: bc => NULL()
end type coupler_1d_bc_type

!--- namelist interface ------------------------------------------------------
! <NAMELIST NAME="flux_exchange_nml">
!   <DATA NAME="z_ref_heat"  TYPE="real"  DEFAULT="2.0">
!    eference height (meters) for temperature and relative humidity 
!    diagnostics (t_ref,rh_ref,del_h,del_q)
!   </DATA>
!   <DATA NAME="z_ref_mom"  TYPE="real"  DEFAULT="10.0">
!    reference height (meters) for momentum diagnostics (u_ref,v_ref,del_m)
!   </DATA>
!   <DATA NAME="ex_u_star_smooth_bug"  TYPE="logical"  DEFAULT="false">
!    By default, the global exchange grid u_star will not be interpolated from 
!    atmospheric grid, this is different from Jakarta behavior and will
!    change answers. So to perserve Jakarta behavior and reproduce answers
!    explicitly set this namelist variable to .true. in input.nml.
!    Talk to mw, ens for details.
!   </DATA>
!   <DATA NAME="do_runoff"  TYPE="logical"  DEFAULT=".TRUE.">
!    Turns on/off the land runoff interpolation to the ocean.
!   </DATA>


  real ::  z_ref_heat =  2.,  &
           z_ref_mom  = 10.
  logical :: ex_u_star_smooth_bug = .false.
  logical :: do_area_weighted_flux = .FALSE.
  logical :: do_runoff = .TRUE.
  logical :: do_forecast = .false.

namelist /flux_exchange_nml/ z_ref_heat, z_ref_mom, ex_u_star_smooth_bug, &
         do_area_weighted_flux, do_runoff, do_forecast
! </NAMELIST>

! ---- allocatable module storage --------------------------------------------
real, allocatable, dimension(:) :: &
     ! NOTE: T canopy is only differet from t_surf over vegetated land
     ex_t_surf,    &   ! surface temperature for radiation calc, degK
     ex_t_surf_miz,&   ! miz
     ex_t_ca,      &   ! near-surface (canopy) air temperature, degK
     ex_p_surf,    &   ! surface pressure
     ex_slp,       &   ! surface pressure

     ex_flux_t,    &   ! sens heat flux
     ex_flux_lw,   &   ! longwave radiation flux

     ex_dhdt_surf, &   ! d(sens.heat.flux)/d(T canopy)
     ex_dedt_surf, &   ! d(water.vap.flux)/d(T canopy)
     ex_dqsatdt_surf, &   ! d(water.vap.flux)/d(q canopy)
     ex_e_q_n,     &
     ex_drdt_surf, &   ! d(LW flux)/d(T surf)
     ex_dhdt_atm,  &   ! d(sens.heat.flux)/d(T atm)
     ex_flux_u,    &   ! u stress on atmosphere
     ex_flux_v,    &   ! v stress on atmosphere
     ex_dtaudu_atm,&   ! d(stress)/d(u)
     ex_dtaudv_atm,&   ! d(stress)/d(v)
     ex_albedo_fix,&
     ex_albedo_vis_dir_fix,&
     ex_albedo_nir_dir_fix,&
     ex_albedo_vis_dif_fix,&
     ex_albedo_nir_dif_fix,&
     ex_old_albedo,&   ! old value of albedo for downward flux calculations
     ex_drag_q,    &   ! q drag.coeff.
     ex_cd_t,      &
     ex_cd_m,      &
     ex_b_star,    &
     ex_u_star,    &
     ex_wind,      &
     ex_z_atm

real, allocatable, dimension(:,:) :: &
     ex_tr_surf,    & ! near-surface tracer fields
     ex_flux_tr,    & ! tracer fluxes
     ex_dfdtr_surf, & ! d(tracer flux)/d(surf tracer)
     ex_dfdtr_atm,  & ! d(tracer flux)/d(atm tracer)
     ex_e_tr_n,     & ! coefficient in implicit scheme 
     ex_f_tr_delt_n   ! coefficient in implicit scheme

logical, allocatable, dimension(:) :: &
     ex_avail       ! true where data on exchange grid are available
real, allocatable, dimension(:) :: ex_e_t_n, ex_f_t_delt_n

integer :: n_atm_tr      ! number of prognostic tracers in the atmos model
integer :: n_atm_tr_tot  ! number of prognostic tracers in the atmos model
integer :: n_lnd_tr      ! number of prognostic tracers in the land model 
integer :: n_lnd_tr_tot  ! number of prognostic tracers in the land model
integer :: n_exch_tr     ! number of tracers exchanged between models

type :: tracer_ind_type
   integer :: atm, lnd ! indices of the tracer in the respective models
end type 
type(tracer_ind_type), allocatable :: tr_table(:) ! table of tracer indices
type :: tracer_exch_ind_type
   integer :: exch = 0  ! exchange grid index
   integer :: lnd = 0   ! land model index
end type tracer_exch_ind_type
type(tracer_exch_ind_type), allocatable :: tr_table_map(:)
integer :: isphum = NO_TRACER       ! index of specific humidity tracer in tracer table
integer :: ico2   = NO_TRACER       ! index of co2 tracer in tracer table

type(coupler_1d_bc_type), save        :: ex_gas_fields_atm  ! gas fields in atm
                     ! Place holder for various atmospheric fields.
type(coupler_1d_bc_type), save        :: ex_gas_fluxes      ! gas flux
                     ! Place holder of intermediate calculations, such as
                     ! piston velocities etc.

real, dimension(3) :: ccc ! for conservation checks

integer ::  runoff_id_diag =-1 
integer, parameter ::  outunit = 6 
contains

!#######################################################################
! <SUBROUTINE NAME="flux_exchange_init">
!  <OVERVIEW>
!   Initialization routine.
!  </OVERVIEW>
!  <DESCRIPTION>
!   Initializes the interpolation routines,diagnostics and boundary data
!  </DESCRIPTION>
!  <TEMPLATE>
!   call flux_exchange_init ( Time, Atm, Land, land_ice_atmos_boundary, dt_atmos, dt_cpld )
!
!  </TEMPLATE>
!  <IN NAME=" Time" TYPE="time_type">
!   current time
!  </IN>
!  <IN NAME="Atm" TYPE="atmos_data_type">
!   A derived data type to specify atmosphere boundary data.
!  </IN>
!  <IN NAME="Land" TYPE="land_data_type">
!   A derived data type to specify land boundary data.
!  </IN>
!  <INOUT NAME="land_ice_atmos_boundary" TYPE="land_ice_atmos_boundary_type">
!   A derived data type to specify properties and fluxes passed from exchange grid to
!   the atmosphere and land.
!  </INOUT>
!  <IN NAME="dt_atmos" TYPE="integer">
!  Atmos time step in secs.
!  </IN>
!  <IN NAME="dt_cpld" TYPE="integer">
!  Coupled time step in secs.
!  </IN>

! Here are some notes on the stuff in Land, which is a variable of type land_data_type
! Dimensioning:
! real, dimension(1,1,1) :: tile_size, t_surf, t_ca, albedo, albedo_vis_dir, 
! real, dimension(1,1,1) :: albedo_vis_dir, albedo_nir_dir, albedo_vis_dif, albedo_nir_dif
! real, dimension(1,1,1) :: rough_mom, rough_heat, rough_scale
! real, dimension(1,1,1,2) :: tr
! real, dimension(1,1) :: discharge, discharge_heat, discharge_snow, discharge_snow_heat
! logical, dimension(1,1,1) :: mask
! integer, dimension(2) :: axes

! Here are some notes on the stuff in Atm
! real, dimension(2,2) :: glon_bnd, glat_bnd, lon_bnd, lat_bnd
! real, dimension(1,1) :: t_bot, z_bot, p_bot, u_bot, v_bot, p_surf, slp, gust, coszen
! real, dimension(1,1) :: flux_sw_dir, flux_sw_dif, flux_sw_down_vis_dir, flux_sw_down_vis_dif
! real, dimension(1,1) :: flux_sw_down_total_dir, flux_sw_down_total_dif
! real, dimension(1,1) :: flux_sw_vis, flux_sw_vis_dir, flux_sw_vis_dif, lprec, fprec
! real, dimension(1,1,5) :: tr_bot

! Here are some notes on the stuff in Land_Ice_Atmos_Boundary
! real, dimension(1,1) :: albedo, albedo_vis_dir, albedo_nir_dir, albedo_vis_dif, albedo_nir_dif
! real, dimension(1,1) :: dt_t, u_flux, v_flux, dtaudu, dtaudv, u_star, b_star, q_star, rough_mom
! real, dimension(1,1,5) :: dt_tr


subroutine flux_exchange_init ( Time, Atm, Land, land_ice_atmos_boundary, dt_atmos, dt_cpld )

  type(time_type),                   intent(in)  :: Time
  type(atmos_data_type),             intent(inout)  :: Atm
  type(land_data_type),              intent(in)  :: Land
! All intent(OUT) derived types with pointer components must be 
! COMPLETELY allocated here and in subroutines called from here;
! NO pointer components should have been allocated before entry if the
! derived type has intent(OUT) otherwise they may be lost.
  type(land_ice_atmos_boundary_type),intent(inout) :: land_ice_atmos_boundary
  integer, optional,                 intent(in)    :: dt_atmos, dt_cpld

  character(len=64), parameter    :: sub_name = 'flux_exchange_init'
  character(len=256), parameter   :: error_header = '==>Error from ' // trim(module_name) //   &
                                                    '(' // trim(sub_name) // '):'
  character(len=256), parameter   :: warn_header = '==>Warning from ' // trim(module_name) //  &
                                                   '(' // trim(sub_name) // '):'
  character(len=256), parameter   :: note_header = '==>Note from ' // trim(module_name) //     &
                                                   '(' // trim(sub_name) // '):'

  integer        :: isg, ieg, jsg, jeg
  integer        :: isc, iec, jsc, jec
  integer        :: isd, ied, jsd, jed
  integer        :: isc2, iec2, jsc2, jec2
  integer        :: nml_unit, io,  i, j
  integer        :: nlon, nlat
  real, dimension(:,:), allocatable :: tmpx(:,:), tmpy(:,:)
  integer :: is, ie, js, je, kd
  character(32) :: tr_name
  logical       :: found

  integer :: n, npes_atm, npes_ocn, npes_all

!-----------------------------------------------------------------------

!
!       initialize atmos_ocean_fluxes
! Setting up flux types, allocates the arrays.
!

!
!       ocean_tracer_flux_init is called first since it has the meaningful value to set
!       for the input/output file names for the tracer flux values used in restarts. These
!       values could be set in the field table, and this ordering allows this.
!       atmos_tracer_flux_init is called last since it will use the values set in 
!       ocean_tracer_flux_init with the exception of atm_tr_index, which can only
!       be meaningfully set from the atmospheric model (not from the field table)
!

!-----------------------------------------------------------------------
!----- read namelist -------
    nml_unit = get_unit()
    open(nml_unit, file='input.nml', form='formatted', action='read', status='old')
    do 
       read (nml_unit, nml=flux_exchange_nml, iostat=io, end=10)
       if (check_nml_error (io, 'flux_exchange_nml')==0) exit ! from loop
    enddo
10  close (nml_unit)

!----- write namelist to logfile -----
    write (logunit,'(/,80("="),/(a))') trim(version), trim(tagname)
    write( logunit, nml=flux_exchange_nml )

!----- find out number of atmospheric prognostic tracers and index of specific 
!      humidity in the tracer table
  n_atm_tr_tot = num_atmos_tracers
  n_atm_tr     = num_atmos_tracers
  n_lnd_tr_tot = num_land_tracers
  n_lnd_tr     = num_land_tracers

  ! assemble the table of tracer number translation by matching names of
  ! prognostic tracers in the atmosphere and surface models; skip all atmospheric
  ! tracers that have no corresponding surface tracers.
  allocate(tr_table(n_atm_tr))
  allocate(tr_table_map(n_atm_tr))
  n = 1
  do i = 1,n_atm_tr
     tr_table(n)%atm = i
     tr_table(n)%lnd = NO_TRACER
     do j = 1,num_land_tracers
       if(land_tracer_names(j) == atmos_tracer_names(i)) then
         tr_table(n)%lnd = j
         exit
       endif
     enddo
     tr_table_map(i)%lnd = tr_table(n)%lnd
     if(tr_table(n)%lnd /= NO_TRACER) then
       tr_table_map(i)%exch = n
       n = n + 1
     endif
  enddo
  n_exch_tr = n - 1
  !
  !     Set up tracer table entries for ocean-atm gas fluxes where the names of tracers in the
  !     atmosphere and ocean may not be equal
  do n = 1, ex_gas_fluxes%num_bcs  !{
    if (ex_gas_fluxes%bc(n)%atm_tr_index .gt. 0) then  !{
      found = .false.
      do i = 1, n_exch_tr  !{
        if (ex_gas_fluxes%bc(n)%atm_tr_index .eq. tr_table(i)%atm) then
          found = .true.
          exit
        endif
      enddo  !} i
      if (.not. found) then
        n_exch_tr = n_exch_tr + 1
        tr_table(n_exch_tr)%atm = ex_gas_fluxes%bc(n)%atm_tr_index
        tr_table(n_exch_tr)%lnd = NO_TRACER ! because this would have been found above
        tr_table_map(n_exch_tr)%exch = n_exch_tr
        tr_table_map(n_exch_tr)%lnd = tr_table(n_exch_tr)%lnd
      endif
    endif  !}
  enddo  !} n
  write(outunit,*) trim(note_header), ' Number of exchanged tracers = ', n_exch_tr
  write(logunit,*) trim(note_header), ' Number of exchanged tracers = ', n_exch_tr
  do i = 1,n_exch_tr
     tr_name = atmos_tracer_names(tr_table(i)%atm)
     write(outunit,*)'Tracer field name :'//trim(tr_name)
     write(logunit,*)'Tracer field name :'//trim(tr_name)
  enddo

  ! find out which tracer is specific humidity

  ! +fix-me-slm+ specific humidity may not be present if we are running with
  ! dry atmosphere. Besides, model may use mixing ratio ('mix_rat') (?). However,
  ! some atmos code also assumes 'sphum' is present, so for now the following
  ! code may be good enough.

  do i = 1,n_exch_tr
     tr_name = atmos_tracer_names(tr_table(i)%atm)
     if(trim(tr_name)=='sphum') then
        isphum = i
     endif
  ! jgj: find out which exchange tracer is co2
     if(trim(tr_name)=='co2') then
        ico2 = i
        write(outunit,*)'Exchange tracer index for '//trim(tr_name),' : ',ico2
     endif
  enddo

  if (isphum==NO_TRACER) then
     call error_mesg('flux_exchange_mod',&
          'tracer "sphum" must be present in the atmosphere', FATAL )
  endif

  if (ico2==NO_TRACER) then
     call error_mesg('flux_exchange_mod',&
          'tracer "co2" not present in the atmosphere', NOTE )
  endif

!-----------------------------------------------------------------------

    if( Atm%pe )then
        nlon = 1
        nlat = 1          
        n_xgrid_sfc = 1
        Land%tile_size(1,1,1) = 1.0

!-----------------------------------------------------------------------
!----- initialize diagnostic fields -----
!----- all fields will be output on the atmospheric grid -----

        call diag_field_init ( Time, Atm%axes(1:2), Land%axes )

!allocate land_ice_atmos_boundary
        is=1; ie=1; js=1; je=1
        allocate( land_ice_atmos_boundary%t(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%albedo(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%albedo_vis_dir(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%albedo_nir_dir(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%albedo_vis_dif(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%albedo_nir_dif(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%land_frac(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%dt_t(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%dt_tr(is:ie,js:je,n_atm_tr) )
        allocate( land_ice_atmos_boundary%u_flux(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%v_flux(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%dtaudu(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%dtaudv(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%u_star(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%b_star(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%q_star(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%rough_mom(is:ie,js:je) )
! initialize boundary values
        land_ice_atmos_boundary%t=273.0
        land_ice_atmos_boundary%albedo=0.0
        land_ice_atmos_boundary%albedo_vis_dir=0.0
        land_ice_atmos_boundary%albedo_nir_dir=0.0
        land_ice_atmos_boundary%albedo_vis_dif=0.0
        land_ice_atmos_boundary%albedo_nir_dif=0.0
        land_ice_atmos_boundary%land_frac=0.0
        land_ice_atmos_boundary%dt_t=0.0
        land_ice_atmos_boundary%dt_tr=0.0
        land_ice_atmos_boundary%u_flux=0.0
        land_ice_atmos_boundary%v_flux=0.0
        land_ice_atmos_boundary%dtaudu=0.0
        land_ice_atmos_boundary%dtaudv=0.0
        land_ice_atmos_boundary%u_star=0.0
        land_ice_atmos_boundary%b_star=0.0
        land_ice_atmos_boundary%q_star=0.0
        land_ice_atmos_boundary%rough_mom=0.01
    end if

    do_init = .false.

  end subroutine flux_exchange_init
! </SUBROUTINE>

!#######################################################################
! <SUBROUTINE NAME="sfc_boundary_layer">
!  <OVERVIEW>
!   Computes explicit fluxes as well as derivatives that will be used to compute an implicit flux correction. 
!  </OVERVIEW>
!  <DESCRIPTION>
!  <PRE>
!  The following quantities in the land_ice_atmos_boundary_type are computed:
!
!     
!         t_surf_atm = surface temperature (used for radiation)    (K)
!         albedo_atm = surface albedo      (used for radiation)    (nondimensional)
!      rough_mom_atm = surface roughness for momentum (m)
!      land_frac_atm = fractional area of land beneath an atmospheric
!                      grid box 
!         dtaudu_atm, dtaudv_atm = derivatives of wind stress w.r.t. the
!                      lowest level wind speed  (Pa/(m/s))
!         flux_u_atm = zonal wind stress  (Pa)
!         flux_v_atm = meridional wind stress (Pa)
!         u_star_atm = friction velocity (m/s)
!         b_star_atm = buoyancy scale    (m2/s)
!
!         (u_star and b_star are defined so that u_star**2 = magnitude
!           of surface stress divided by density of air at the surface, 
!           and u_star*b_star = buoyancy flux at the surface)
!
!   </PRE>
!  </DESCRIPTION>

!  <TEMPLATE>
!   call sfc_boundary_layer ( dt, Time, Atm, Land, Boundary )
!
!  </TEMPLATE>
!  <IN NAME=" dt" TYPE="real">
!   time step. 
!  </IN>
!  <IN NAME="Time" TYPE="time_type">
!   current time
!  </IN>
!  <INOUT NAME="Atm" TYPE="atmos_data_type">
!   A derived data type to specify atmosphere boundary data.
!  </INOUT>
!  <INOUT NAME="Land" TYPE="land_data_type">
!   A derived data type to specify land boundary data.
!  </INOUT>
!  <INOUT NAME="Boundary" TYPE="land_ice_atmos_boundary_type">
!   A derived data type to specify properties and fluxes passed from exchange grid to
!   the atmosphere and land.
!  </INOUT>
!
subroutine sfc_boundary_layer ( dt, Time, Atm, Land, Land_Ice_Atmos_Boundary )

  real,                  intent(in)  :: dt
  type(time_type),       intent(in)  :: Time
  type(atmos_data_type), intent(inout)  :: Atm
  type(land_data_type),  intent(inout)  :: Land
  type(land_ice_atmos_boundary_type), intent(inout) :: Land_Ice_Atmos_Boundary

  ! ---- local vars ----------------------------------------------------------
  real, dimension(n_xgrid_sfc) :: &
       ex_albedo,     &
       ex_albedo_vis_dir,     &
       ex_albedo_nir_dir,     &
       ex_albedo_vis_dif,     &
       ex_albedo_nir_dif,     &
       ex_land_frac,  &
       ex_t_atm,      & 
       ex_p_atm,      &
       ex_u_atm, ex_v_atm,    &
       ex_gust,       &
       ex_t_surf4,    &
       ex_u_surf, ex_v_surf,  &
       ex_rough_mom, ex_rough_heat, ex_rough_moist, &
       ex_rough_scale,&
       ex_q_star,     &
       ex_cd_q,       &
       ex_ref, ex_ref_u, ex_ref_v, ex_u10, &
       ex_ref2,       &
       ex_t_ref,      &
       ex_qs_ref,     &
       ex_qs_ref_cmip,     &
       ex_del_m,      &
       ex_del_h,      &
       ex_del_q

  real, dimension(n_xgrid_sfc,n_exch_tr) :: ex_tr_atm
! jgj: added for co2_atm diagnostic
  real, dimension(n_xgrid_sfc)           :: ex_co2_atm_dvmr
  real, dimension(size(Land_Ice_Atmos_Boundary%t,1),size(Land_Ice_Atmos_Boundary%t,2)) :: diag_atm
  real, dimension(size(Land%t_ca, 1),size(Land%t_ca,2), size(Land%t_ca,3)) :: diag_land
  real    :: zrefm, zrefh
  logical :: used
  character(32) :: tr_name ! tracer name
  integer :: tr, n, m ! tracer indices
  integer :: i
  integer :: year, month, day, hour, minute, second

  ! [1] check that the module was initialized
!   <ERROR MSG="must call flux_exchange_init first " STATUS="FATAL">
!      flux_exchange_init has not been called before calling sfc_boundary_layer.
!   </ERROR>
  if (do_init) call error_mesg ('flux_exchange_mod',  &
       'must call flux_exchange_init first', FATAL)
  ! [2] allocate storage for variables that are also used in flux_up_to_atmos
  allocate ( &
       ex_t_surf   (n_xgrid_sfc),  &
       ex_t_surf_miz(n_xgrid_sfc), &
       ex_p_surf   (n_xgrid_sfc),  &
       ex_slp      (n_xgrid_sfc),  &
       ex_t_ca     (n_xgrid_sfc),  &
       ex_dhdt_surf(n_xgrid_sfc),  &
       ex_dedt_surf(n_xgrid_sfc),  &
       ex_dqsatdt_surf(n_xgrid_sfc),  &
       ex_drdt_surf(n_xgrid_sfc),  &
       ex_dhdt_atm (n_xgrid_sfc),  &
       ex_flux_t   (n_xgrid_sfc),  &
       ex_flux_lw  (n_xgrid_sfc),  &
       ex_drag_q   (n_xgrid_sfc),  &
       ex_avail    (n_xgrid_sfc),  &
       ex_f_t_delt_n(n_xgrid_sfc), &

       ex_tr_surf     (n_xgrid_sfc, n_exch_tr), &
       ex_dfdtr_surf  (n_xgrid_sfc, n_exch_tr), &
       ex_dfdtr_atm   (n_xgrid_sfc, n_exch_tr), &
       ex_flux_tr     (n_xgrid_sfc, n_exch_tr), &
       ex_f_tr_delt_n (n_xgrid_sfc, n_exch_tr), &
       ex_e_tr_n      (n_xgrid_sfc, n_exch_tr), &

! MOD these were moved from local ! so they can be passed to flux down
       ex_flux_u(n_xgrid_sfc),    &
       ex_flux_v(n_xgrid_sfc),    &
       ex_dtaudu_atm(n_xgrid_sfc),&
       ex_dtaudv_atm(n_xgrid_sfc),&

! values added for LM3
       ex_cd_t     (n_xgrid_sfc),  &
       ex_cd_m     (n_xgrid_sfc),  &
       ex_b_star   (n_xgrid_sfc),  &
       ex_u_star   (n_xgrid_sfc),  &
       ex_wind     (n_xgrid_sfc),  &
       ex_z_atm    (n_xgrid_sfc),  &

       ex_e_t_n    (n_xgrid_sfc),  &
       ex_e_q_n    (n_xgrid_sfc)   )
  do n = 1, ex_gas_fields_atm%num_bcs  !{
    do m = 1, ex_gas_fields_atm%bc(n)%num_fields  !{
      if (associated(ex_gas_fields_atm%bc(n)%field(m)%values)) then  !{
        call mpp_error( FATAL, 'sfc_boundary_layer: ex_gas_fields_atm already allocated.' )
      endif  !}
      allocate ( ex_gas_fields_atm%bc(n)%field(m)%values(n_xgrid_sfc) )
      ex_gas_fields_atm%bc(n)%field(m)%values = 0.0
    enddo  !} m
  enddo  !} n

  do n = 1, ex_gas_fluxes%num_bcs  !{
    do m = 1, ex_gas_fluxes%bc(n)%num_fields  !{
      if (associated(ex_gas_fluxes%bc(n)%field(m)%values)) then  !{
        call mpp_error( FATAL, 'sfc_boundary_layer: ex_gas_fluxes already allocated.' )
      endif  !}
      allocate ( ex_gas_fluxes%bc(n)%field(m)%values(n_xgrid_sfc) )
      ex_gas_fluxes%bc(n)%field(m)%values = 0.0
    enddo  !} m
  enddo  !} n

!
!       Call the atmosphere tracer driver to gather the data needed for extra gas tracers
! For ocean only model

!  call atmos_get_fields_for_flux(Atm)

  ! [3] initialize some values on exchange grid: this is actually a safeguard
  ! against using undefined values
  ex_t_surf   = 200.
  ex_u_surf   =   0.
  ex_v_surf   =   0.
  ex_albedo = 0. ! bw 
  ex_albedo_vis_dir = 0.
  ex_albedo_nir_dir = 0.
  ex_albedo_vis_dif = 0.
  ex_albedo_nir_dif = 0.

  !---- do not use if relax time /= 0 ----
  ex_cd_t = 0.0
  ex_cd_m = 0.0
  ex_cd_q = 0.0

!---- put atmosphere quantities onto exchange grid ----

  ! [4] put all the qantities we need onto exchange grid
  ! [4.1] put atmosphere quantities onto exchange grid
  if (do_forecast) then
    ex_t_surf_miz(1) = Atm%Surf_diff%sst_miz(1,1)
  endif

  ex_gust(1)   = Atm%gust(1,1)
  ex_t_atm(1)  = Atm%t_bot(1,1)
  ex_z_atm(1)  = Atm%z_bot(1,1)
  ex_p_atm(1)  = Atm%p_bot(1,1)
  ex_u_atm(1)  = Atm%u_bot(1,1)
  ex_v_atm(1)  = Atm%v_bot(1,1)
  ex_p_surf(1) = Atm%p_surf(1,1)
  ex_slp(1)    = Atm%slp(1,1)
  ex_gust(1)   = Atm%gust(1,1)

! put atmosphere bottom layer tracer data onto exchange grid
  do tr = 1,n_exch_tr
     ex_tr_atm(1,tr) = Atm%tr_bot(1,1,tr_table(tr)%atm)
  enddo

  ex_tr_surf = ex_tr_atm
  ex_t_ca = ex_t_surf

  if (do_forecast) then
    ex_t_surf_miz(1) = Land%t_surf(1,1,1)
    ex_t_ca(:) = ex_t_surf_miz(:)
  end if

  ex_t_surf(1) = Land%t_surf(1,1,1)
  ex_t_ca(1) = Land%t_ca(1,1,1)
  ex_rough_mom(1) = Land%rough_mom(1,1,1)
  ex_rough_heat(1) = Land%rough_heat(1,1,1)
  ex_rough_moist(1) = Land%rough_heat(1,1,1)
  ex_albedo(1) = Land%albedo(1,1,1)
  ex_albedo_vis_dir(1) = Land%albedo_vis_dir(1,1,1)
  ex_albedo_nir_dir(1) = Land%albedo_nir_dir(1,1,1)
  ex_albedo_vis_dif(1) = Land%albedo_vis_dif(1,1,1)
  ex_albedo_nir_dif(1) = Land%albedo_nir_dif(1,1,1)
  ex_rough_scale = ex_rough_mom
  ex_rough_scale(1) = Land%rough_scale(1,1,1)
 
  do tr = 1,n_exch_tr
     n = tr_table(tr)%lnd
     if(n /= NO_TRACER ) then
        ex_tr_surf(1,tr) = Land%tr(1,1,1,n)
     else
        ! do nothing, since ex_tr_surf is prefilled with ex_tr_atm, and therefore
        ! fluxes will be 0
     endif
  enddo

  ex_land_frac = 0.0
  call put_logical_to_real (Land%mask, ex_land_frac)

  if (do_forecast) then
     ex_t_surf = ex_t_surf_miz
  end if

  ! [5] compute explicit fluxes and tendencies at all available points ---
  call surface_flux (&
       ex_t_atm, ex_tr_atm(:,isphum),  ex_u_atm, ex_v_atm,  ex_p_atm,  ex_z_atm,  &
       ex_p_surf,ex_t_surf, ex_t_ca,  ex_tr_surf(:,isphum),                       &
       ex_u_surf, ex_v_surf,                                                      &
       ex_rough_mom, ex_rough_heat, ex_rough_moist, ex_rough_scale, ex_gust,      &
       ex_flux_t, ex_flux_tr(:,isphum), ex_flux_lw, ex_flux_u, ex_flux_v,         &
       ex_cd_m,   ex_cd_t, ex_cd_q,                                               &
       ex_wind,   ex_u_star, ex_b_star, ex_q_star,                                &
       ex_dhdt_surf, ex_dedt_surf, ex_dfdtr_surf(:,isphum),  ex_drdt_surf,        &
       ex_dhdt_atm,  ex_dfdtr_atm(:,isphum),  ex_dtaudu_atm, ex_dtaudv_atm, dt)

  zrefm = 10.0
  zrefh = z_ref_heat
  ex_avail = .true.
  !      ---- optimize calculation ----
  call mo_profile ( zrefm, zrefh, ex_z_atm,   ex_rough_mom, &
       ex_rough_heat, ex_rough_moist,          &
       ex_u_star, ex_b_star, ex_q_star,        &
       ex_del_m, ex_del_h, ex_del_q, ex_avail  )
  ex_u10 = 0.
  ex_ref_u = ex_u_surf + (ex_u_atm-ex_u_surf) * ex_del_m 
  ex_ref_v = ex_v_surf + (ex_v_atm-ex_v_surf) * ex_del_m
  ex_u10 = sqrt(ex_ref_u**2 + ex_ref_v**2)

  ! fill derivatives for all tracers
  ! F = C0*u*rho*delta_q, C0*u*rho is the same for all tracers, copy from sphum
  do tr = 1,n_exch_tr
     if (tr==isphum) cycle
     ex_dfdtr_atm  (:,tr) = ex_dfdtr_atm  (:,isphum)
     ex_dfdtr_surf (:,tr) = ex_dfdtr_surf (:,isphum)
     ex_flux_tr    (:,tr) = ex_dfdtr_surf(:,tr)*(ex_tr_surf(:,tr)-ex_tr_atm(:,tr))
  enddo

! Combine explicit ocean flux and implicit land flux of extra flux fields.

  ! NB: names of the fields are constructed using tracer name and certain 
  ! prefixes / suffixes. For example, for the tracer named "sphum" (specific humidity) they will be:
  ! "ex_flux_sphum", "ex_dfdsphum_surf", and "ex_dfdsphum_atm".
  ! 
  ! For sensible heat flux names are "ex_flux_t", "ex_dhdt_surf", and "ex_dhdt_atm"; 
  ! despite the name those are actually in energy units, W/m2, W/(m2 degK), and
  ! W/(m2 degK) respectively

  ex_drag_q = ex_wind*ex_cd_q
  ! [6] get mean quantities on atmosphere grid
  ! [6.1] compute t surf for radiation
  ex_t_surf4 = ex_t_surf ** 4

  ! [6.2] put relevant quantities onto atmospheric boundary
  Land_Ice_Atmos_Boundary%t(1,1) = ex_t_surf4(1)
  Land_Ice_Atmos_Boundary%albedo(1,1) = ex_albedo(1)
  Land_Ice_Atmos_Boundary%albedo_vis_dir(1,1) = ex_albedo_vis_dir(1)
  Land_Ice_Atmos_Boundary%albedo_nir_dir(1,1) = ex_albedo_nir_dir(1)
  Land_Ice_Atmos_Boundary%albedo_vis_dif(1,1) = ex_albedo_vis_dif(1)
  Land_Ice_Atmos_Boundary%albedo_nir_dif(1,1) = ex_albedo_nir_dif(1)
  Land_Ice_Atmos_Boundary%rough_mom(1,1) = ex_rough_mom(1)
  Land_Ice_Atmos_Boundary%land_frac(1,1) = ex_land_frac(1)

  Land_Ice_Atmos_Boundary%u_flux(1,1) = ex_flux_u(1)
  Land_Ice_Atmos_Boundary%v_flux(1,1) = ex_flux_v(1)
  Land_Ice_Atmos_Boundary%dtaudu(1,1) = ex_dtaudu_atm(1)
  Land_Ice_Atmos_Boundary%dtaudv(1,1) = ex_dtaudv_atm(1)
  Land_Ice_Atmos_Boundary%u_star(1,1) = ex_u_star(1)
  Land_Ice_Atmos_Boundary%b_star(1,1) = ex_b_star(1)
  Land_Ice_Atmos_Boundary%q_star(1,1) = ex_q_star(1)

  Land_Ice_Atmos_Boundary%t = Land_Ice_Atmos_Boundary%t ** 0.25

  ! [6.3] save atmos albedo fix and old albedo (for downward SW flux calculations)
  ! on exchange grid
  ! allocate ( ex_old_albedo(n_xgrid_sfc)  )
  ! ex_old_albedo = ex_albedo
  
!!  STILL NEEDED   ????
!! IS THIS CORRECT ??
  allocate ( ex_albedo_fix(n_xgrid_sfc) )
  ex_albedo_fix = 0.
  ex_albedo_fix(1) = Land_Ice_Atmos_Boundary%albedo(1,1)
  ex_albedo_fix = (1.0-ex_albedo) / (1.0-ex_albedo_fix)

  allocate ( ex_albedo_vis_dir_fix(n_xgrid_sfc) )
  ex_albedo_vis_dir_fix = 0.
  ex_albedo_vis_dir_fix(1) = Land_Ice_Atmos_Boundary%albedo_vis_dir(1,1)
  ex_albedo_vis_dir_fix = (1.0-ex_albedo_vis_dir) /  (1.0-ex_albedo_vis_dir_fix)
  allocate ( ex_albedo_nir_dir_fix(n_xgrid_sfc) )
  ex_albedo_nir_dir_fix = 0.
  ex_albedo_nir_dir_fix(1) = Land_Ice_Atmos_Boundary%albedo_nir_dir(1,1)
  ex_albedo_nir_dir_fix = (1.0-ex_albedo_nir_dir) /  (1.0-ex_albedo_nir_dir_fix)
  allocate ( ex_albedo_vis_dif_fix(n_xgrid_sfc) )
  ex_albedo_vis_dif_fix = 0.
  ex_albedo_vis_dif_fix(1) = Land_Ice_Atmos_Boundary%albedo_vis_dif(1,1)
  ex_albedo_vis_dif_fix = (1.0-ex_albedo_vis_dif) /   (1.0-ex_albedo_vis_dif_fix)
  allocate ( ex_albedo_nir_dif_fix(n_xgrid_sfc) )
  ex_albedo_nir_dif_fix = 0.
  ex_albedo_nir_dif_fix(1) = Land_Ice_Atmos_Boundary%albedo_nir_dif(1,1)
  ex_albedo_nir_dif_fix = (1.0-ex_albedo_nir_dif) /   (1.0-ex_albedo_nir_dif_fix)

  !=======================================================================
  ! [7] diagnostics section

  call get_date(Time, year, month, day, hour, minute, second)
  !------- save static fields first time only ------
  if (first_static) then

     !------- land fraction ------
     write(id_land_mask, diag_format) year, month, day, hour, minute, second, Land_Ice_Atmos_Boundary%land_frac

     first_static = .false.
  endif

  !------- drag coeff moisture -----------
  if ( id_wind > 0 ) then
     diag_atm(1,1) = ex_wind(1)
     write(id_wind, diag_format) year, month, day, hour, minute, second, diag_atm
  endif
  !------- drag coeff moisture -----------
  if ( id_drag_moist > 0 ) then
     diag_atm(1,1) = ex_cd_q(1)
     write(id_drag_moist, diag_format) year, month, day, hour, minute, second, diag_atm
  endif

  !------- drag coeff heat -----------
  if ( id_drag_heat > 0 ) then
     diag_atm(1,1) = ex_cd_t(1)
     write(id_drag_heat, diag_format) year, month, day, hour, minute, second, diag_atm
  endif
  
  !------- drag coeff momemtum -----------
  if ( id_drag_mom > 0 ) then
     diag_atm(1,1) = ex_cd_m(1)
     write(id_drag_mom, diag_format) year, month, day, hour, minute, second, diag_atm
  endif
  
  !------- roughness moisture -----------
  if ( id_rough_moist > 0 ) then
     diag_atm(1,1) = ex_rough_moist(1)
     write(id_rough_moist, diag_format) year, month, day, hour, minute, second, diag_atm
  endif
  
  !------- roughness heat -----------
  if ( id_rough_heat > 0 ) then
     diag_atm(1,1) = ex_rough_heat(1)
     write(id_rough_heat, diag_format) year, month, day, hour, minute, second, diag_atm
  endif
  
  !------- roughness momemtum -----------
!  write(id_rough_mom, diag_format) year, month, day, hour, minute, second, Land_Ice_Atmos_Boundary%rough_mom
  
  !------- friction velocity -----------
!  write(id_u_star, diag_format) year, month, day, hour, minute, second, Land_Ice_Atmos_Boundary%u_star
  
  !------- bouyancy -----------
!  write(id_b_star, diag_format) year, month, day, hour, minute, second, Land_Ice_Atmos_Boundary%b_star

  !------- moisture scale -----------
!  write(id_q_star, diag_format) year, month, day, hour, minute, second, Land_Ice_Atmos_Boundary%q_star

  !-----------------------------------------------------------------------
  !------ diagnostics for fields at bottom atmospheric level ------
  
  if ( id_t_atm > 0 ) then
     diag_atm(1,1) = ex_t_atm(1)
!     write(id_t_atm, diag_format) year, month, day, hour, minute, second, diag_atm
  endif
  
  if ( id_u_atm > 0 ) then
     diag_atm(1,1) = ex_u_atm(1)
!     write(id_u_atm, diag_format) year, month, day, hour, minute, second, diag_atm
  endif
  
  if ( id_v_atm > 0 ) then
     diag_atm(1,1) = ex_v_atm(1)
!     write(id_v_atm, diag_format) year, month, day, hour, minute, second, diag_atm
  endif

  do tr = 1,n_exch_tr
     tr_name = atmos_tracer_names(tr_table(tr)%atm)
     if ( id_tr_atm(tr) > 0 ) then
        diag_atm(1,1) = ex_tr_atm(1,tr)
!        write(id_tr_atm(tr), diag_format) year, month, day, hour, minute, second, diag_atm
     endif
!!jgj: add dryvmr co2_atm
! - slm Mar 25 2010: moved to resolve interdependence of diagnostic fields
     if ( id_co2_atm_dvmr > 0 .and. trim(tr_name)=='co2') then
        ex_co2_atm_dvmr = (ex_tr_atm(:,tr) / (1.0 - ex_tr_atm(:,isphum))) * WTMAIR/WTMCO2
        diag_atm(1,1) = ex_co2_atm_dvmr(1)
!        write(id_co2_atm_dvmr, diag_format) year, month, day, hour, minute, second, diag_atm
     endif
  enddo

  ! - slm, Mar 25, 2002
  if ( id_p_atm > 0 ) then
     diag_atm(1,1) = ex_p_atm(1)
!     write(id_p_atm, diag_format) year, month, day, hour, minute, second, diag_atm
  endif
  if ( id_z_atm > 0 ) then
     diag_atm(1,1) = ex_z_atm(1)
!     write(id_z_atm, diag_format) year, month, day, hour, minute, second, diag_atm
  endif
  if ( id_gust > 0 ) then
     diag_atm(1,1) = ex_gust(1)
!     write(id_gust, diag_format) year, month, day, hour, minute, second, diag_atm
  endif

  ! - bw, Sep 17, 2007
  if ( id_slp > 0 ) then
     diag_atm(1,1) = ex_slp(1)
!     write(id_slp, diag_format) year, month, day, hour, minute, second, diag_atm
  endif

  !-----------------------------------------------------------------------
  !--------- diagnostics for fields at reference level ---------
  
  if ( id_t_ref > 0 .or. id_rh_ref > 0 .or. &
       id_u_ref > 0 .or. id_v_ref  > 0 .or. id_wind_ref > 0 .or. &
       id_q_ref > 0 .or. id_q_ref_land > 0 .or. &
       id_t_ref_land > 0 .or. id_rh_ref_land > 0 .or. &
       id_rh_ref_cmip >0 .or. &
       id_u_ref_land > 0 .or. id_v_ref_land  > 0 ) then
     
     zrefm = z_ref_mom
     zrefh = z_ref_heat
     !      ---- optimize calculation ----
     if ( id_t_ref <= 0 ) zrefh = zrefm
     
     call mo_profile ( zrefm, zrefh, ex_z_atm,   ex_rough_mom, &
          ex_rough_heat, ex_rough_moist,          &
          ex_u_star, ex_b_star, ex_q_star,        &
          ex_del_m, ex_del_h, ex_del_q, ex_avail  )

     !    ------- reference relative humidity -----------
     if ( id_rh_ref > 0 .or. id_rh_ref_land > 0 .or. &
          id_rh_ref_cmip > 0 .or. &
          id_q_ref > 0 .or. id_q_ref_land >0 ) then
          ex_ref = 1.0e-06
          ex_ref   = ex_tr_surf(:,isphum) + (ex_tr_atm(:,isphum)-ex_tr_surf(:,isphum)) * ex_del_q
        if(id_q_ref > 0) then
           diag_atm(1,1) = ex_ref(1)
           write(id_q_ref, diag_format) year, month, day, hour, minute, second, diag_atm
        endif
        if(id_q_ref_land > 0) then
           diag_land(1,1,1) = ex_ref(1)
           do n=1,size(diag_land,3)
             write(id_q_ref_land, tile_diag_format) year, month, day, hour, minute, second, n, diag_land(:,:,n)
           enddo
        endif
        ex_t_ref = 200.
        ex_t_ref = ex_t_ca + (ex_t_atm-ex_t_ca) * ex_del_h
        call compute_qs (ex_t_ref, ex_p_surf, ex_qs_ref, q = ex_ref)
        call compute_qs (ex_t_ref, ex_p_surf, ex_qs_ref_cmip,  &
                         q = ex_ref, es_over_liq_and_ice = .true.)
! remove cap on relative humidity -- this mod requested by cjg, ljd
!RSH    ex_ref    = MIN(100.,100.*ex_ref/ex_qs_ref)
        ex_ref2   = 100.*ex_ref/ex_qs_ref_cmip
        ex_ref    = 100.*ex_ref/ex_qs_ref

        if ( id_rh_ref_land > 0 ) then
           diag_land(1,1,1) = ex_ref(1)
           do n=1,size(diag_land,3)
             write(id_rh_ref_land, tile_diag_format) year, month, day, hour, minute, second, n, diag_land(:,:,n)
           enddo
        endif
        if(id_rh_ref > 0) then
           diag_atm(1,1) = ex_ref(1)
           write(id_rh_ref, diag_format) year, month, day, hour, minute, second, diag_atm
        endif
        if(id_rh_ref_cmip > 0) then
           diag_atm(1,1) = ex_ref2(1)
           write(id_rh_ref_cmip, diag_format) year, month, day, hour, minute, second, diag_atm
        endif
     endif

     !    ------- reference temp -----------
     if ( id_t_ref > 0 .or. id_t_ref_land > 0 ) then
        ex_ref = ex_t_ca + (ex_t_atm-ex_t_ca) * ex_del_h
        if (id_t_ref_land > 0) then
           diag_land(1,1,1) = ex_ref(1)
           do n=1,size(diag_land,3)
!             write(id_t_ref_land, tile_diag_format) year, month, day, hour, minute, second, n, diag_land(:,:,n)
           enddo
        endif
        if ( id_t_ref > 0 ) then
           diag_atm(1,1) = ex_ref(1)
!           write(id_t_ref, diag_format) year, month, day, hour, minute, second, diag_atm
        endif
     endif

     !    ------- reference u comp -----------
     if ( id_u_ref > 0 .or. id_u_ref_land > 0) then
        ex_ref = ex_u_surf + (ex_u_atm-ex_u_surf) * ex_del_m
        if ( id_u_ref_land > 0 ) then
           diag_land(1,1,1) = ex_ref(1)
           do n=1,size(diag_land,3)
!             write(id_u_ref_land, tile_diag_format) year, month, day, hour, minute, second, n, diag_land(:,:,n)
           enddo
        endif
        if ( id_u_ref > 0 ) then
           diag_atm(1,1) = ex_ref(1)
!           write(id_u_ref, diag_format) year, month, day, hour, minute, second, diag_atm
        endif
     endif

     !    ------- reference v comp -----------
     if ( id_v_ref > 0 .or. id_v_ref_land > 0 ) then
        ex_ref = ex_v_surf + (ex_v_atm-ex_v_surf) * ex_del_m
        if ( id_v_ref_land > 0 ) then
           diag_land(1,1,1) = ex_ref(1)
           do n=1,size(diag_land,3)
!             write(id_v_ref_land, tile_diag_format) year, month, day, hour, minute, second, n, diag_land(:,:,n)
           enddo
        endif
        if ( id_v_ref > 0 ) then
           diag_atm(1,1) = ex_ref(1)
!           write(id_v_ref, diag_format) year, month, day, hour, minute, second, diag_atm
        endif
     endif

     !    ------- reference-level absolute wind -----------
     if ( id_wind_ref > 0 ) then
        ex_ref = sqrt((ex_u_surf + (ex_u_atm-ex_u_surf) * ex_del_m)**2 &
                        +(ex_v_surf + (ex_v_atm-ex_v_surf) * ex_del_m)**2)
        diag_atm(1,1) = ex_ref(1)
!        write(id_wind_ref, diag_format) year, month, day, hour, minute, second, diag_atm
     endif

     !    ------- interp factor for heat ------
     if ( id_del_h > 0 ) then
        diag_atm(1,1) = ex_del_h(1)
!        write(id_del_h, diag_format) year, month, day, hour, minute, second, diag_atm
     endif

     !    ------- interp factor for momentum ------
     if ( id_del_m > 0 ) then
        diag_atm(1,1) = ex_del_m(1)
!        write(id_del_m, diag_format) year, month, day, hour, minute, second, diag_atm
     endif

     !    ------- interp factor for moisture ------
     if ( id_del_q > 0 ) then
        diag_atm(1,1) = ex_del_q(1)
!        write(id_del_q, diag_format) year, month, day, hour, minute, second, diag_atm
     endif

  endif
  ! topographic roughness scale
  if(id_rough_scale>0) then
     diag_atm(1,1) = (log(ex_z_atm(1)/ex_rough_mom(1)+1)/log(ex_z_atm(1)/ex_rough_scale(1)+1))**2
!     write(id_rough_scale, diag_format) year, month, day, hour, minute, second, diag_atm
  endif

!=======================================================================

end subroutine sfc_boundary_layer
! </SUBROUTINE>

!#######################################################################

! <SUBROUTINE NAME="flux_down_from_atmos">
!  <OVERVIEW>
!   Returns fluxes and derivatives corrected for the implicit treatment of atmospheric 
!   diffusive fluxes, as well as the increments in the temperature and specific humidity 
!   of the lowest atmospheric layer due to all explicit processes as well as the diffusive 
!   fluxes through the top of this layer. 
!  </OVERVIEW>
!  <DESCRIPTION>
!  <PRE>
!    The following elements from Atmos_boundary are used as input: 
!
!        flux_u_atm = zonal wind stress (Pa)  
!        flux_v_atm = meridional wind stress (Pa)
!
!
!    The following elements of Land_boundary are output: 
!
!       flux_t_land = sensible heat flux (W/m2)
!       flux_q_land = specific humidity flux (Kg/(m2 s)
!      flux_lw_land = net longwave flux (W/m2), uncorrected for
!                     changes in surface temperature
!      flux_sw_land = net shortwave flux (W/m2)
!         dhdt_land = derivative of sensible heat flux w.r.t.
!                     surface temperature (on land model grid)  (W/(m2 K)
!         dedt_land = derivative of specific humidity flux w.r.t.
!                     surface temperature (on land model grid)  (Kg/(m2 s K)
!         drdt_land = derivative of upward longwave flux w.r.t.
!                     surface temperature (on land model grid) (W/(m2 K)
!        lprec_land = liquid precipitation, mass for one time step
!                      (Kg/m2)
!        fprec_land = frozen precipitation, mass for one time step
!                      (Kg/m2)
!
!   </PRE>
!  </DESCRIPTION>
!  <TEMPLATE>
!   call flux_down_from_atmos (Time, Atm, Land, Atmos_boundary, Land_boundary )
!
!  </TEMPLATE>
!  <IN NAME="Time" TYPE="time_type">
!   current time
!  </IN>
!  <INOUT NAME="Atm" TYPE="atmos_data_type">
!   A derived data type to specify atmosphere boundary data.
!  </INOUT>
!  <IN NAME="Land" TYPE="land_data_type">
!   A derived data type to specify land boundary data.
!  </IN>
!  <IN NAME="Atmos_boundary" TYPE="land_ice_atmos_boundary_type">
!   A derived data type to specify properties and fluxes passed from exchange grid to
!   the atmosphere and land.
!  </IN>
!  <INOUT NAME="Land_boundary" TYPE="atmos_land_boundary_type">
!   A derived data type to specify properties and fluxes passed from atmosphere to land.
!  </INOUT>
!
subroutine flux_down_from_atmos (Time, Atm, Land, Atmos_boundary, Land_boundary )

  type(time_type),       intent(in) :: Time
  type(atmos_data_type), intent(inout) :: Atm
  type(land_data_type),  intent(in) :: Land
  type(land_ice_atmos_boundary_type),intent(in) :: Atmos_boundary
  type(atmos_land_boundary_type),    intent(inout):: Land_boundary

  real, dimension(n_xgrid_sfc) :: ex_flux_sw, ex_flux_lwd, &
       ex_flux_sw_dir, ex_flux_sw_dif,  &
       ex_flux_sw_down_vis_dir, ex_flux_sw_down_total_dir,  &
       ex_flux_sw_down_vis_dif, ex_flux_sw_down_total_dif,  &
       ex_flux_sw_vis, &
       ex_flux_sw_vis_dir, &
       ex_flux_sw_vis_dif, &
       ex_lprec, ex_fprec,      &
       ex_tprec, & ! temperature of precipitation, currently equal to atm T
       ex_u_star_smooth,        &
       ex_coszen

  real, dimension(n_xgrid_sfc) :: ex_gamma  , ex_dtmass,  &
       ex_delta_t, ex_delta_u, ex_delta_v, ex_dflux_t

  real, dimension(n_xgrid_sfc,n_exch_tr) :: &
       ex_delta_tr, & ! tracer tendencies
       ex_dflux_tr    ! fracer flux change

  real    :: cp_inv
  logical :: used
  logical :: ov
  integer :: ier, year, month, day, hour, minute, second

  character(32) :: tr_name ! name of the tracer
  integer :: tr, n, m ! tracer indices

  ov = .FALSE.

!---- put atmosphere quantities onto exchange grid ----

  ex_flux_sw_dir     = 0.0
  ex_flux_sw_vis_dir = 0.0
  ex_flux_sw_dif     = 0.0
  ex_flux_sw_vis_dif = 0.0
  ex_flux_lwd        = 0.0                           
  ex_flux_sw_dir(1) = Atm%flux_sw_dir(1,1)
  ex_flux_sw_vis_dir(1) = Atm%flux_sw_vis_dir(1,1)
  ex_flux_sw_dif(1) = Atm%flux_sw_dif(1,1)
  ex_flux_sw_vis_dif(1) = Atm%flux_sw_vis_dif(1,1)
  ex_flux_sw_down_vis_dir(1) = Atm%flux_sw_down_vis_dir(1,1)
  ex_flux_sw_down_total_dir(1) = Atm%flux_sw_down_total_dir(1,1)
  ex_flux_sw_down_vis_dif(1) = Atm%flux_sw_down_vis_dif(1,1)
  ex_flux_sw_down_total_dif(1) = Atm%flux_sw_down_total_dif(1,1)
  ex_flux_lwd(1) = Atm%flux_lw(1,1)
  ex_lprec(1) = Atm%lprec(1,1)
  ex_fprec(1) = Atm%fprec(1,1)
  ex_tprec(1) = Atm%t_bot(1,1)

  ex_coszen(1) = Atm%coszen(1,1)

! MOD changed the following two lines to put Atmos%surf_diff%delta_u and v
! on exchange grid instead of the stresses themselves so that only the 
! implicit corrections are filtered through the atmospheric grid not the
! stresses themselves
  ex_delta_u = 0.0; ex_delta_v = 0.0
  ex_delta_u(1) = Atm%Surf_Diff%delta_u(1,1)
  ex_delta_v(1) = Atm%Surf_Diff%delta_v(1,1)

  ! MOD update stresses using atmos delta's but derivatives on exchange grid
  ex_flux_u = ex_flux_u + ex_delta_u*ex_dtaudu_atm
  ex_flux_v = ex_flux_v + ex_delta_v*ex_dtaudv_atm

!-----------------------------------------------------------------------
!---- adjust sw flux for albedo variations on exch grid ----
!---- adjust 4 categories (vis/nir dir/dif) separately  ----
!-----------------------------------------------------------------------

  ex_flux_sw_dir = ex_flux_sw_dir - ex_flux_sw_vis_dir     ! temporarily nir/dir
  ex_flux_sw_dir = ex_flux_sw_dir * ex_albedo_nir_dir_fix  ! fix nir/dir
  ex_flux_sw_vis_dir = ex_flux_sw_vis_dir * ex_albedo_vis_dir_fix ! fix vis/dir
  ex_flux_sw_dir = ex_flux_sw_dir + ex_flux_sw_vis_dir     ! back to total dir

  ex_flux_sw_dif = ex_flux_sw_dif - ex_flux_sw_vis_dif     ! temporarily nir/dif
  ex_flux_sw_dif = ex_flux_sw_dif * ex_albedo_nir_dif_fix  ! fix nir/dif
  ex_flux_sw_vis_dif = ex_flux_sw_vis_dif * ex_albedo_vis_dif_fix ! fix vis/dif
  ex_flux_sw_dif = ex_flux_sw_dif + ex_flux_sw_vis_dif     ! back to total dif

  ex_flux_sw_vis = ex_flux_sw_vis_dir + ex_flux_sw_vis_dif ! legacy, remove later
  ex_flux_sw     = ex_flux_sw_dir     + ex_flux_sw_dif     ! legacy, remove later

  deallocate ( ex_albedo_fix )
  deallocate ( ex_albedo_vis_dir_fix )
  deallocate ( ex_albedo_nir_dir_fix )
  deallocate ( ex_albedo_vis_dif_fix )
  deallocate ( ex_albedo_nir_dif_fix )
!----- compute net longwave flux (down-up) -----
  ! (note: lw up already in ex_flux_lw)

  ex_flux_lw = ex_flux_lwd - ex_flux_lw

!-----------------------------------------------------------------------
!----- adjust fluxes for implicit dependence on atmosphere ----


  ex_dtmass(1) = Atm%Surf_Diff%dtmass(1,1)
  ex_delta_t(1) = Atm%Surf_Diff%delta_t(1,1)
  ex_dflux_t(1) = Atm%Surf_Diff%dflux_t(1,1)
  do tr = 1,n_exch_tr
     n = tr_table(tr)%atm
     ex_delta_tr(1,tr) = Atm%Surf_Diff%delta_tr(1,1,n)
     ex_dflux_tr(1,tr) = Atm%Surf_Diff%dflux_tr(1,1,n)
  enddo

  cp_inv = 1.0/cp_air


  ! temperature

  ex_gamma      =  1./ (1.0 - ex_dtmass*(ex_dflux_t + ex_dhdt_atm*cp_inv))
  ex_e_t_n      =  ex_dtmass*ex_dhdt_surf*cp_inv*ex_gamma
  ex_f_t_delt_n = (ex_delta_t + ex_dtmass * ex_flux_t*cp_inv) * ex_gamma    
     
  ex_flux_t     =  ex_flux_t        + ex_dhdt_atm * ex_f_t_delt_n 
  ex_dhdt_surf  =  ex_dhdt_surf     + ex_dhdt_atm * ex_e_t_n   

  ! moisture vs. surface temperture, assuming saturation
  ex_gamma   =  1.0 / (1.0 - ex_dtmass*(ex_dflux_tr(:,isphum) + ex_dfdtr_atm(:,isphum)))
  ex_e_q_n      =  ex_dtmass * ex_dedt_surf * ex_gamma
  ex_dedt_surf  =  ex_dedt_surf + ex_dfdtr_atm(:,isphum) * ex_e_q_n
  do tr = 1,n_exch_tr
     ex_gamma   =  1.0 / (1.0 - ex_dtmass*(ex_dflux_tr(:,tr) + ex_dfdtr_atm(:,tr)))

     ex_e_tr_n(:,tr)      =  ex_dtmass*ex_dfdtr_surf(:,tr)*ex_gamma
     ex_f_tr_delt_n(:,tr) = (ex_delta_tr(:,tr)+ex_dtmass*ex_flux_tr(:,tr))*ex_gamma    

     ex_flux_tr(:,tr)     =  ex_flux_tr(:,tr) + ex_dfdtr_atm(:,tr)*ex_f_tr_delt_n(:,tr) 
     ex_dfdtr_surf(:,tr)  =  ex_dfdtr_surf(:,tr) + ex_dfdtr_atm(:,tr)*ex_e_tr_n(:,tr)
  enddo
!-----------------------------------------------------------------------
!---- output fields on the land grid -------

  Land_boundary%t_flux(1,1,1) = ex_flux_t(1)
  Land_boundary%sw_flux(1,1,1) = ex_flux_sw(1)
  Land_boundary%sw_flux_down_vis_dir(1,1,1) = ex_flux_sw_down_vis_dir(1)
  Land_boundary%sw_flux_down_total_dir(1,1,1) = ex_flux_sw_down_total_dir(1)
  Land_boundary%sw_flux_down_vis_dif(1,1,1) = ex_flux_sw_down_vis_dif(1)
  Land_boundary%sw_flux_down_total_dif(1,1,1) = ex_flux_sw_down_total_dif(1)
  Land_boundary%lw_flux(1,1,1) = ex_flux_lw(1)
  Land_boundary%dhdt(1,1,1) = ex_dhdt_surf(1)
  Land_boundary%drdt(1,1,1) = ex_drdt_surf(1)
  Land_boundary%p_surf(1,1,1) = ex_p_surf(1)

  Land_boundary%lprec(1,1,1) = ex_lprec(1)
  Land_boundary%fprec(1,1,1) = ex_fprec(1)
  Land_boundary%tprec(1,1,1) = ex_tprec(1)

  if(associated(Land_boundary%drag_q)) then
     Land_boundary%drag_q(1,1,1) = ex_drag_q(1)
  endif
  if(associated(Land_boundary%lwdn_flux)) then
     Land_boundary%lwdn_flux(1,1,1) = ex_flux_lwd(1)
  endif
  if(associated(Land_boundary%cd_m)) then
     Land_boundary%cd_m(1,1,1) = ex_cd_m(1)
  endif
  if(associated(Land_boundary%cd_t)) then
     Land_boundary%cd_t(1,1,1) = ex_cd_t(1)
  endif
  if(associated(Land_boundary%bstar)) then
     Land_boundary%bstar(1,1,1) = ex_b_star(1)
  endif
  if(associated(Land_boundary%ustar)) then
     Land_boundary%ustar(1,1,1) = ex_u_star(1)
  endif
  if(associated(Land_boundary%wind)) then
     Land_boundary%wind(1,1,1) = ex_wind(1)
  endif
  if(associated(Land_boundary%z_bot)) then
     Land_boundary%z_bot(1,1,1) = ex_z_atm(1)
  endif

  Land_boundary%tr_flux(:,:,:,:) = 0.0
  Land_boundary%dfdtr(:,:,:,:) = 0.0
  do tr = 1,n_exch_tr
     n = tr_table(tr)%lnd
     if(n /= NO_TRACER ) then
        Land_boundary%tr_flux(1,1,1,n) = ex_flux_tr(1,tr)
        Land_boundary%dfdtr(1,1,1,n) = ex_dfdtr_surf(1,tr)
     endif
  enddo

  deallocate ( ex_flux_u, ex_flux_v, ex_dtaudu_atm, ex_dtaudv_atm)

  !=======================================================================
  !-------------------- diagnostics section ------------------------------
  call get_date(Time, year, month, day, hour, minute, second)

  !------- zonal wind stress -----------
!  write(id_u_flux, diag_format) year, month, day, hour, minute, second, Atmos_boundary%u_flux

  !------- meridional wind stress -----------
!  write(id_v_flux, diag_format) year, month, day, hour, minute, second, Atmos_boundary%v_flux

!=======================================================================

  end subroutine flux_down_from_atmos
! </SUBROUTINE>

!#######################################################################
! <SUBROUTINE NAME="flux_up_to_atmos">
!  <OVERVIEW>
!   Corrects the fluxes for consistency with the new surface temperatures in land model.
!  </OVERVIEW>
!  <DESCRIPTION>
!   Corrects the fluxes for consistency with the new surface temperatures in land model.
!   Final increments for temperature and specific humidity in the 
!   lowest atmospheric layer are computed and returned to the atmospheric model
!   so that it can finalize the increments in the rest of the atmosphere. 
!  <PRE>
!
!   The following elements of the land_ice_atmos_boundary_type are computed:
!        dt_t  = temperature change at the lowest
!                 atmospheric level (deg k)
!        dt_q  = specific humidity change at the lowest
!                 atmospheric level (kg/kg)
!  </PRE>
!  </DESCRIPTION>
!  <TEMPLATE>
!   call flux_up_to_atmos ( Time, Land, Land_Ice_Atmos_Boundary, Land_boundary )
!
!  </TEMPLATE>
!  <IN NAME=" Time" TYPE="time_type">
!   Current time.
!  </IN>
!  <INOUT NAME="Land" TYPE="land_data_type">
!   A derived data type to specify land boundary data.
!  </INOUT>
!  <INOUT NAME="Land_Ice_Atmos_Boundary" TYPE="land_ice_atmos_boundary_type">
!   A derived data type to specify properties and fluxes passed from exchange grid to
!   the atmosphere and land.
!  </INOUT>
!
subroutine flux_up_to_atmos ( Time, Land, Land_Ice_Atmos_Boundary, Land_boundary )

  type(time_type),      intent(in)  :: Time
  type(land_data_type), intent(inout)  :: Land
  type(land_ice_atmos_boundary_type), intent(inout) :: Land_Ice_Atmos_Boundary
  type(atmos_land_boundary_type) :: Land_boundary

  real, dimension(n_xgrid_sfc) ::  &
       ex_t_surf_new, &
       ex_dt_t_surf,  &
       ex_delta_t_n,  &
       ex_t_ca_new,   &
       ex_dt_t_ca
  real, dimension(n_xgrid_sfc,n_exch_tr) :: &
       ex_tr_surf_new,    & ! updated tracer values at the surface
       ex_dt_tr_surf,     & ! tendency of tracers at the surface
       ex_delta_tr_n
! jgj: added for co2_surf diagnostic 
  real, dimension(n_xgrid_sfc) :: &
       ex_co2_surf_dvmr   ! updated CO2 tracer values at the surface (dry vmr)

  real, dimension(size(Land_Ice_Atmos_Boundary%dt_t,1),size(Land_Ice_Atmos_Boundary%dt_t,2)) :: diag_atm, &
       evap_atm
  real, dimension(size(Land_boundary%lprec,1), size(Land_boundary%lprec,2), size(Land_boundary%lprec,3)) :: data_lnd, diag_land
  logical :: used

  integer :: tr       ! tracer index
  character(32) :: tr_name ! tracer name
  integer :: n, i, m, ier, year, month, day, hour, minute, second

  !-----------------------------------------------------------------------

  !----- compute surface temperature change -----

  ex_t_surf_new = 200.0

  ex_t_ca_new(1) = Land%t_ca(1,1,1)
  ex_t_surf_new(1) = Land%t_surf(1,1,1)

  ex_dt_t_ca   = ex_t_ca_new   - ex_t_ca   ! changes in near-surface T
  ex_dt_t_surf = ex_t_surf_new - ex_t_surf ! changes in radiative T

  !-----------------------------------------------------------------------
  !-----  adjust fluxes and atmospheric increments for 
  !-----  implicit dependence on surface temperature -----
  do tr = 1,n_exch_tr
     ! set up updated surface tracer field so that flux to atmos for absent
     ! tracers is zero
     if (ex_dfdtr_surf(1,tr)/=0) then
        ex_dt_tr_surf(1,tr) = -ex_flux_tr(1,tr)/ex_dfdtr_surf(1,tr)
     else
        ex_dt_tr_surf(1,tr) = 0
     endif
     ex_tr_surf_new(1,tr) = ex_tr_surf(1,tr)+ex_dt_tr_surf(1,tr)
     ! get all tracers available from land, and calculate changes in near-tracer field
     n = tr_table(tr)%lnd
     if(n /= NO_TRACER ) then
        ex_tr_surf_new(:,tr) = Land%tr(1,1,1,n)
     endif

     ! get all tracers available from ocean here 

     ! update tracer tendencies in the atmosphere
     ex_dt_tr_surf(:,tr) = ex_tr_surf_new(:,tr) - ex_tr_surf(:,tr)
     ex_delta_tr_n(:,tr) = ex_f_tr_delt_n(:,tr) + ex_dt_tr_surf(:,tr) * ex_e_tr_n(:,tr)
     ex_flux_tr(:,tr)    = ex_flux_tr(:,tr)     + ex_dt_tr_surf(:,tr) * ex_dfdtr_surf(:,tr)
  enddo

  do tr=1,n_exch_tr
     ! get updated tracer tendency on the atmospheic grid
     n=tr_table(tr)%atm
     Land_Ice_Atmos_Boundary%dt_tr(1,1,n) = ex_delta_tr_n(1,tr)
  enddo

  ex_delta_t_n = 0.0

  ex_flux_t     = ex_flux_t  + ex_dt_t_ca   * ex_dhdt_surf
  ex_flux_lw    = ex_flux_lw - ex_dt_t_surf * ex_drdt_surf
  ex_delta_t_n  = ex_f_t_delt_n  + ex_dt_t_ca*ex_e_t_n

  !-----------------------------------------------------------------------
  !---- get mean quantites on atmospheric grid ----

  Land_Ice_Atmos_Boundary%dt_t(1,1) = ex_delta_t_n(1)

  !=======================================================================
  !-------------------- diagnostics section ------------------------------
  call get_date(Time, year, month, day, hour, minute, second)

  !------- new surface temperature -----------
  if ( id_t_surf > 0 ) then
     diag_atm(1,1) = ex_t_surf_new(1)
!     write(id_t_surf, diag_format) year, month, day, hour, minute, second, diag_atm
  endif


  ! + slm, Mar 27 2002
  ! ------ new canopy temperature --------
  !   NOTE, that in the particular case of LM2 t_ca is identical to t_surf,
  !   but this will be changed in future version of the land madel
  if ( id_t_ca > 0 ) then
     diag_atm(1,1) = ex_t_ca_new(1)
  endif

  !------- updated surface tracer fields ------
  do tr=1,n_exch_tr
     tr_name = atmos_tracer_names(tr_table(tr)%atm)
     if ( id_tr_surf(tr) > 0 ) then
        diag_atm(1,1) = ex_tr_surf_new(1,tr)
!        write(id_tr_surf(tr), diag_format) year, month, day, hour, minute, second, diag_atm
     endif
!!jgj:  add dryvmr co2_surf
! - slm Mar 25, 2010: moved to resolve interdependence of diagnostic fields
     if ( id_co2_surf_dvmr > 0 .and. trim(tr_name)=='co2') then
       ex_co2_surf_dvmr = (ex_tr_surf_new(:,tr) / (1.0 - ex_tr_surf_new(:,isphum))) * WTMAIR/WTMCO2
       diag_atm(1,1) = ex_co2_surf_dvmr(1)
!       write(id_co2_surf_dvmr, diag_format) year, month, day, hour, minute, second, diag_atm
     endif
  enddo

  !------- sensible heat flux -----------
  if ( id_t_flux > 0 ) then
     diag_atm(1,1) = ex_flux_t(1)
!     write(id_t_flux, diag_format) year, month, day, hour, minute, second, diag_atm
  endif

  !------- net longwave flux -----------
  if ( id_r_flux > 0 ) then
     diag_atm(1,1) = ex_flux_lw(1)
 !    write(id_r_flux, diag_format) year, month, day, hour, minute, second, diag_atm
  endif

  !------- tracer fluxes ------------
  ! tr_mol_flux diagnostic will be correct for co2 tracer only. 
  ! will need update code to use correct molar mass for tracers other than co2
  do tr=1,n_exch_tr
     if ( id_tr_flux(tr) > 0 .or. id_tr_mol_flux(tr) > 0 ) then
        diag_atm(1,1) = ex_flux_tr(1,tr)
        if(id_tr_flux(tr) > 0) then
!          write(id_tr_flux(tr), diag_format) year, month, day, hour, minute, second, diag_atm
        endif
        if(id_tr_mol_flux(tr) > 0) then
!          write(id_tr_mol_flux(tr), diag_format) year, month, day, hour, minute, second, diag_atm*1000./WTMCO2
        endif
     endif
  enddo

  !-----------------------------------------------------------------------
  !---- accumulate global integral of evaporation (mm/day) -----
  evap_atm(1,1) = ex_flux_tr(1,isphum)
!  write(id_q_flux, diag_format) year, month, day, hour, minute, second, evap_atm
  if( id_q_flux_land > 0 ) then
     diag_land(1,1,1) = ex_flux_tr(1,isphum)
     do n=1,size(diag_land,3)
!       write(id_q_flux_land, tile_diag_format) year, month, day, hour, minute, second, n, diag_land(:,:,n)
     enddo
  endif

  data_lnd(1,1,1) = ex_flux_tr(1,isphum)

  !=======================================================================
  !---- deallocate module storage ----
  deallocate ( &
       ex_t_surf   ,  &
       ex_t_surf_miz, &
       ex_p_surf   ,  &
       ex_slp      ,  &
       ex_t_ca     ,  &
       ex_dhdt_surf,  &
       ex_dedt_surf,  &
       ex_dqsatdt_surf,  &
       ex_drdt_surf,  &
       ex_dhdt_atm ,  &
       ex_flux_t   ,  &
       ex_flux_lw  ,  &
       ex_drag_q   ,  &
       ex_avail    ,  &
       ex_f_t_delt_n, &
       ex_tr_surf  ,  &
       
  ex_dfdtr_surf  , &
       ex_dfdtr_atm   , &
       ex_flux_tr     , &
       ex_f_tr_delt_n , &
       ex_e_tr_n      , &
       
  ex_e_t_n    ,  &
       ex_e_q_n    ,  &
       ! values added for LM3
       ex_cd_t     ,  &
       ex_cd_m     ,  &
       ex_b_star   ,  &
       ex_u_star   ,  &
       ex_wind     ,  &
       ex_z_atm )

! Extra fluxes
  do n = 1, ex_gas_fields_atm%num_bcs  !{
     do m = 1, ex_gas_fields_atm%bc(n)%num_fields  !{
        deallocate ( ex_gas_fields_atm%bc(n)%field(m)%values )
        nullify ( ex_gas_fields_atm%bc(n)%field(m)%values )
     enddo  !} m
  enddo  !} n
  do n = 1, ex_gas_fluxes%num_bcs  !{
     do m = 1, ex_gas_fluxes%bc(n)%num_fields  !{
        deallocate ( ex_gas_fluxes%bc(n)%field(m)%values )
        nullify ( ex_gas_fluxes%bc(n)%field(m)%values )
     enddo  !} m
  enddo  !} n

!-----------------------------------------------------------------------

end subroutine flux_up_to_atmos
! </SUBROUTINE>

!#######################################################################

subroutine put_logical_to_real (mask, ex_mask)

  logical, intent(in)    :: mask(:,:,:)
  real   , intent(inout) :: ex_mask(:)

  !-----------------------------------------------------------------------
  !    puts land model masks (with partitions) onto the
  !    exchange grid as a real array (1.=true, 0.=false)
  !-----------------------------------------------------------------------

  real, dimension(size(mask,1),size(mask,2),size(mask,3)) :: rmask
  
  where (mask)
     rmask = 1.0
  elsewhere
     rmask = 0.0
  endwhere

  ex_mask(1) = rmask(1,1,1)

end subroutine put_logical_to_real

!#######################################################################

subroutine diag_field_init ( Time, atmos_axes, land_axes )

  type(time_type), intent(in) :: Time
  integer,         intent(in) :: atmos_axes(2)
  integer,         intent(in) :: land_axes(2)

  integer :: iref
  character(len=6) :: label_zm, label_zh
  real, dimension(2) :: trange = (/  100., 400. /), &
       vrange = (/ -400., 400. /), &
       frange = (/ -0.01, 1.01 /)
  character(len=32)  :: name     ! name of the tracer
  character(len=128) :: longname ! long name of the tracer
  integer            :: tr       ! tracer index
!-----------------------------------------------------------------------
!  initializes diagnostic fields that may be output from this module
!  (the id numbers may be referenced anywhere in this module)
!-----------------------------------------------------------------------

  !------ labels for diagnostics -------
  !  (z_ref_mom, z_ref_heat are namelist variables)

  iref = int(z_ref_mom+0.5)
  if ( real(iref) == z_ref_mom ) then
     write (label_zm,105) iref
     if (iref < 10) write (label_zm,100) iref
  else
     write (label_zm,110) z_ref_mom
  endif

  iref = int(z_ref_heat+0.5)
  if ( real(iref) == z_ref_heat ) then
     write (label_zh,105) iref
     if (iref < 10) write (label_zh,100) iref
  else
     write (label_zh,110) z_ref_heat
  endif

100 format (i1,' m',3x)
105 format (i2,' m',2x)
110 format (f4.1,' m')

  !--------- initialize static diagnostic fields --------------------

  id_land_mask = get_unit()
  open(unit=id_land_mask, file='DIAGNOSTICS/'//mod_name//'_land_mask', form='formatted', action='write', position='rewind')
  write(id_land_mask,'(a)') 'fractional amount of land  unit=none'
  
  !--------- initialize diagnostic fields --------------------

  id_wind = get_unit()
  open(unit=id_wind, file='DIAGNOSTICS/'//mod_name//'_wind', form='formatted', action='write', position='rewind')
  write(id_wind,'(a)') 'wind speed for flux calculations   units=m/s'
  
  id_drag_moist = get_unit()
  open(unit=id_drag_moist, file='DIAGNOSTICS/'//mod_name//'_drag_moist', form='formatted', action='write', position='rewind')
  write(id_drag_moist,'(a)') 'drag coeff for moisture   units=none'
  
  id_drag_heat = get_unit()
  open(unit=id_drag_heat, file='DIAGNOSTICS/'//mod_name//'_drag_heat', form='formatted', action='write', position='rewind')
  write(id_drag_heat,'(a)') 'drag coeff for heat   units=none'
  
  id_drag_mom = get_unit()
  open(unit=id_drag_mom, file='DIAGNOSTICS/'//mod_name//'_drag_mom', form='formatted', action='write', position='rewind')
  write(id_drag_mom,'(a)') 'drag coeff for momentum   units=none'
  
  id_rough_moist = get_unit()
  open(unit=id_rough_moist, file='DIAGNOSTICS/'//mod_name//'_rough_moist', form='formatted', action='write', position='rewind')
  write(id_rough_moist,'(a)') 'surface roughness for moisture   units=m'

  id_rough_heat = get_unit()
  open(unit=id_rough_heat, file='DIAGNOSTICS/'//mod_name//'_rough_heat', form='formatted', action='write', position='rewind')
  write(id_rough_heat,'(a)') 'surface roughness for heat   units=m'

  id_rough_mom = get_unit()
  open(unit=id_rough_mom, file='DIAGNOSTICS/'//mod_name//'_rough_mom', form='formatted', action='write', position='rewind')
  write(id_rough_mom,'(a)') 'surface roughness for momentum   units=m'

  id_u_star = get_unit()
  open(unit=id_u_star, file='DIAGNOSTICS/'//mod_name//'_u_star', form='formatted', action='write', position='rewind')
  write(id_u_star,'(a)') 'friction velocity   units=m/s'

  id_b_star = get_unit()
  open(unit=id_b_star, file='DIAGNOSTICS/'//mod_name//'_b_star', form='formatted', action='write', position='rewind')
  write(id_b_star,'(a)') 'buoyancy scale   units=m/s2'

  id_q_star = get_unit()
  open(unit=id_q_star, file='DIAGNOSTICS/'//mod_name//'_q_star', form='formatted', action='write', position='rewind')
  write(id_q_star,'(a)') 'moisture scale   units=kg water/kg air'

  id_u_flux = get_unit()
  open(unit=id_u_flux, file='DIAGNOSTICS/'//mod_name//'_tau_x', form='formatted', action='write', position='rewind')
  write(id_u_flux,'(a)') 'zonal wind stress   units=pa'

  id_v_flux = get_unit()
  open(unit=id_v_flux, file='DIAGNOSTICS/'//mod_name//'_tau_y', form='formatted', action='write', position='rewind')
  write(id_v_flux,'(a)') 'meridional wind stress   units=pa'

  id_t_surf = get_unit()
  open(unit=id_t_surf, file='DIAGNOSTICS/'//mod_name//'_t_surf', form='formatted', action='write', position='rewind')
  write(id_t_surf,'(a)') 'surface temperature   units=deg_k'

  ! + slm, Mar 25, 2002 -- add diagnositcs for t_ca, q_ca, and q_atm
  id_t_ca = get_unit()
  open(unit=id_t_ca, file='DIAGNOSTICS/'//mod_name//'_t_ca', form='formatted', action='write', position='rewind')
  write(id_t_ca,'(a)') 'canopy air temperature   units=deg_k'

  ! - slm, Mar 25, 2002
  id_z_atm = get_unit()
  open(unit=id_z_atm, file='DIAGNOSTICS/'//mod_name//'_z_atm', form='formatted', action='write', position='rewind')
  write(id_z_atm,'(a)') 'height of btm level   units=m'

  id_p_atm = get_unit()
  open(unit=id_p_atm, file='DIAGNOSTICS/'//mod_name//'_p_atm', form='formatted', action='write', position='rewind')
  write(id_p_atm,'(a)') 'pressure at btm level   units=pa'

  ! - bw, Mar 25, 2002 -- added diagnostic slp
  id_slp = get_unit()
  open(unit=id_slp, file='DIAGNOSTICS/'//mod_name//'_slp', form='formatted', action='write', position='rewind')
  write(id_slp,'(a)') 'sea level pressure   units=pa'

  id_gust = get_unit()
  open(unit=id_gust, file='DIAGNOSTICS/'//mod_name//'_gust', form='formatted', action='write', position='rewind')
  write(id_gust,'(a)') 'gust scale   units=m/s'

  id_t_flux = get_unit()
  open(unit=id_t_flux, file='DIAGNOSTICS/'//mod_name//'_shflx', form='formatted', action='write', position='rewind')
  write(id_t_flux,'(a)') 'sensible heat flux   units=w/m2'

  id_r_flux = get_unit()
  open(unit=id_r_flux, file='DIAGNOSTICS/'//mod_name//'_lwflx', form='formatted', action='write', position='rewind')
  write(id_r_flux,'(a)') 'net (down-up) longwave flux   units=w/m2'

  id_t_atm = get_unit()
  open(unit=id_t_atm, file='DIAGNOSTICS/'//mod_name//'_t_atm', form='formatted', action='write', position='rewind')
  write(id_t_atm,'(a)') 'temperature at btm level   units=deg_k'

  id_u_atm = get_unit()
  open(unit=id_u_atm, file='DIAGNOSTICS/'//mod_name//'_u_atm', form='formatted', action='write', position='rewind')
  write(id_u_atm,'(a)') 'u wind component at btm level   units=m/s'

  id_v_atm = get_unit()
  open(unit=id_v_atm, file='DIAGNOSTICS/'//mod_name//'_v_atm', form='formatted', action='write', position='rewind')
  write(id_v_atm,'(a)') 'v wind component at btm level   units=m/s'

  id_t_ref = get_unit()
  open(unit=id_t_ref, file='DIAGNOSTICS/'//mod_name//'_t_ref', form='formatted', action='write', position='rewind')
  write(id_t_ref,'(a)') 'temperature at '//label_zh//'   units=deg_k'

  id_rh_ref = get_unit()
  open(unit=id_rh_ref, file='DIAGNOSTICS/'//mod_name//'_rh_ref', form='formatted', action='write', position='rewind')
  write(id_rh_ref,'(a)') 'relative humidity at '//label_zh//'   units=percent'

  id_rh_ref_cmip = get_unit()
  open(unit=id_rh_ref_cmip, file='DIAGNOSTICS/'//mod_name//'_rh_ref_cmip', form='formatted', action='write', position='rewind')
  write(id_rh_ref_cmip,'(a)') 'relative humidity at '//label_zh//'   units=percent'

  id_u_ref = get_unit()
  open(unit=id_u_ref, file='DIAGNOSTICS/'//mod_name//'_u_ref', form='formatted', action='write', position='rewind')
  write(id_u_ref,'(a)') 'zonal wind component at '//label_zm//'   units=m/s'

  id_v_ref = get_unit()
  open(unit=id_v_ref, file='DIAGNOSTICS/'//mod_name//'_v_ref', form='formatted', action='write', position='rewind')
  write(id_v_ref,'(a)') 'meridional wind component at '//label_zm//'   units=m/s'

  id_wind_ref = get_unit()
  open(unit=id_wind_ref, file='DIAGNOSTICS/'//mod_name//'_wind_ref', form='formatted', action='write', position='rewind')
  write(id_wind_ref,'(a)') 'absolute value of wind at '//label_zm//'   units=m/s'

  id_del_h = get_unit()
  open(unit=id_del_h, file='DIAGNOSTICS/'//mod_name//'_del_h', form='formatted', action='write', position='rewind')
  write(id_del_h,'(a)') 'ref height interp factor for heat   units=none'

  id_del_m = get_unit()
  open(unit=id_del_m, file='DIAGNOSTICS/'//mod_name//'_del_m', form='formatted', action='write', position='rewind')
  write(id_del_m,'(a)') 'ref height interp factor for momentum   units=none'

  id_del_q = get_unit()
  open(unit=id_del_q, file='DIAGNOSTICS/'//mod_name//'_del_q', form='formatted', action='write', position='rewind')
  write(id_del_q,'(a)') 'ref height interp factor for moisture   units=none'

  ! + slm Jun 02, 2002 -- diagnostics of reference values over the land
  id_t_ref_land = get_unit()
  open(unit=id_t_ref_land, file='DIAGNOSTICS/'//mod_name//'_t_ref_land', form='formatted', action='write', position='rewind')
  write(id_t_ref_land,'(a)') 'temperature at '//label_zh//'   units=deg_k'

  id_rh_ref_land = get_unit()
  open(unit=id_rh_ref_land, file='DIAGNOSTICS/'//mod_name//'_rh_ref_land', form='formatted', action='write', position='rewind')
  write(id_rh_ref_land,'(a)') 'relative humidity at '//label_zh//'   units=percent'

  id_u_ref_land = get_unit()
  open(unit=id_u_ref_land, file='DIAGNOSTICS/'//mod_name//'_u_ref_land', form='formatted', action='write', position='rewind')
  write(id_u_ref_land,'(a)') 'zonal wind component at '//label_zh//'   units=m/s'

  id_v_ref_land = get_unit()
  open(unit=id_v_ref_land, file='DIAGNOSTICS/'//mod_name//'_v_ref_land', form='formatted', action='write', position='rewind')
  write(id_v_ref_land,'(a)') 'meridional wind component at '//label_zh//'   units=m/s'

  ! - slm Jun 02, 2002
  id_q_ref = get_unit()
  open(unit=id_q_ref, file='DIAGNOSTICS/'//mod_name//'_q_ref', form='formatted', action='write', position='rewind')
  write(id_q_ref,'(a)') 'specific humidity at '//label_zh//'   units=kg/kg'

  id_q_ref_land = get_unit()
  open(unit=id_q_ref_land, file='DIAGNOSTICS/'//mod_name//'_q_ref_land', form='formatted', action='write', position='rewind')
  write(id_q_ref_land,'(a)') 'specific humidity at '//label_zh//'   units=kg/kg'

  id_rough_scale = get_unit()
  open(unit=id_rough_scale, file='DIAGNOSTICS/'//mod_name//'_rough_scale', form='formatted', action='write', position='rewind')
  write(id_rough_scale,'(a)') 'topographic scaling factor for momentum drag   units=1'
!-----------------------------------------------------------------------

  allocate(id_tr_atm(n_exch_tr))
  allocate(id_tr_surf(n_exch_tr))
  allocate(id_tr_flux(n_exch_tr))
  allocate(id_tr_mol_flux(n_exch_tr))

  do tr = 1, n_exch_tr
     name = atmos_tracer_names(tr_table(tr)%atm)
     longname = atmos_tracer_longnames(tr_table(tr)%atm)
     id_tr_atm(tr) = get_unit()
     open(unit=id_tr_atm(tr), file='DIAGNOSTICS/'//mod_name//'_'//trim(name), form='formatted', action='write', position='rewind')
     write(id_tr_atm(tr),'(a)') trim(longname)//'  units='//trim(atmos_tracer_units(tr))

     id_tr_surf(tr) = get_unit()
     open(unit=id_tr_surf(tr), file='DIAGNOSTICS/'//mod_name//'_'//trim(name)//'_surf', form='formatted', action='write', position='rewind')
     write(id_tr_surf(tr),'(a)') trim(longname)//' at the surface   units='//trim(atmos_tracer_units(tr))

     id_tr_flux(tr) = get_unit()
     open(unit=id_tr_flux(tr), file='DIAGNOSTICS/'//mod_name//'_'//trim(name)//'_flux', form='formatted', action='write', position='rewind')
     write(id_tr_flux(tr),'(a)') 'flux of '//trim(longname)//'   units='//trim(atmos_tracer_units(tr))//' kg air/(m2 s)'
!! add dryvmr co2_surf and co2_atm
     if ( trim(name)=='co2') then
! - slm Mar 25, 2010: moved registration of mol_flux inside 'if' to disable 
! saving incorrect results (mol fluxes for other tracers computed with CO2 molar 
! mass)
       id_tr_mol_flux(tr) = get_unit()
       open(unit=id_tr_mol_flux(tr), file='DIAGNOSTICS/'//mod_name//'_'//trim(name)//'_mol_flux', form='formatted', action='write', position='rewind')
       write(id_tr_mol_flux(tr),'(a)') 'flux of '//trim(longname)//'   units=mol CO2/(m2 s)'

       id_co2_atm_dvmr = get_unit()
       open(unit=id_co2_atm_dvmr, file='DIAGNOSTICS/'//mod_name//'_'//trim(name)//'_atm_dvmr', form='formatted', action='write', position='rewind')
       write(id_co2_atm_dvmr,'(a)') trim(longname)//' at btm level   units=mol CO2 /mol air'

       id_co2_surf_dvmr = get_unit()
       open(unit=id_co2_surf_dvmr, file='DIAGNOSTICS/'//mod_name//'_'//trim(name)//'_surf_dvmr', form='formatted', action='write', position='rewind')
       write(id_co2_surf_dvmr,'(a)') trim(longname)//' at the surface   units=mol CO2 /mol air'
     else
       id_tr_mol_flux(tr) = -1
     endif
  enddo

  id_q_flux = get_unit()
  open(unit=id_q_flux, file='DIAGNOSTICS/'//mod_name//'_evap', form='formatted', action='write', position='rewind')
  write(id_q_flux,'(a)') 'evaporation rate   units=kg/m2/s'

  id_q_flux_land = get_unit()
  open(unit=id_q_flux_land, file='DIAGNOSTICS/'//mod_name//'_evap_land', form='formatted', action='write', position='rewind')
  write(id_q_flux_land,'(a)') 'evaporation rate over land   units=kg/m2/s'

  end subroutine diag_field_init

!######################################################################################
! Divide data by area while avoiding zero area elements
  subroutine divide_by_area(data, area)
    real, intent(inout) :: data(:,:)
    real, intent(in)    :: area(:,:)

    if(size(data, dim=1) /= size(area, dim=1) .or. size(data, dim=2) /= size(area, dim=2)) then
       ! no op
       return
    endif

    where(area /= 0) 
       data = data / area
    end where

  end subroutine divide_by_area
!#######################################################################

! <DIAGFIELDS>
!   <NETCDF NAME="land_mask" UNITS="none">
!     fractional amount of land
!   </NETCDF>
!   <NETCDF NAME="wind" UNITS="m/s">
!     wind speed for flux calculations
!   </NETCDF>
!   <NETCDF NAME="drag_moist" UNITS="none">
!     drag coeff for moisture
!   </NETCDF>
!   <NETCDF NAME="drag_heat" UNITS="none">
!     drag coeff for heat
!   </NETCDF>
!   <NETCDF NAME="drag_mom" UNITS="none">
!     drag coeff for momentum
!   </NETCDF>
!   <NETCDF NAME="rough_moist" UNITS="m">
!     surface roughness for moisture
!   </NETCDF>
!   <NETCDF NAME="rough_heat" UNITS="m">
!     surface roughness for heat
!   </NETCDF>
!   <NETCDF NAME="rough_mom" UNITS="m">
!     surface roughness for momentum
!   </NETCDF>
!   <NETCDF NAME="u_star" UNITS="m/s">
!     friction velocity
!   </NETCDF>
!   <NETCDF NAME="b_star" UNITS="m/s">
!     buoyancy scale
!   </NETCDF>
!   <NETCDF NAME="q_star" UNITS="kg water/kg air">
!     moisture scale
!   </NETCDF>
!   <NETCDF NAME="t_atm" UNITS="deg_k">
!     temperature at btm level
!   </NETCDF>
!   <NETCDF NAME="u_atm" UNITS="m/s">
!     u wind component at btm level
!   </NETCDF>
!   <NETCDF NAME="v_atm" UNITS="m/s">
!     v wind component at btm level
!   </NETCDF>
!   <NETCDF NAME="q_atm" UNITS="kg/kg">
!     specific humidity at btm level
!   </NETCDF>
!   <NETCDF NAME="p_atm" UNITS="pa">
!     pressure at btm level
!   </NETCDF>
!   <NETCDF NAME="z_atm" UNITS="m">
!     height of btm level
!   </NETCDF>
!   <NETCDF NAME="gust" UNITS="m/s">
!     gust scale 
!   </NETCDF>
!   <NETCDF NAME="rh_ref" UNITS="percent">
!     relative humidity at ref height
!   </NETCDF>
!   <NETCDF NAME="t_ref" UNITS="deg_k">
!    temperature at ref height
!   </NETCDF>
!   <NETCDF NAME="u_ref" UNITS="m/s">
!    zonal wind component at ref height
!   </NETCDF>
!   <NETCDF NAME="v_ref" UNITS="m/s">
!    meridional wind component at ref height 
!   </NETCDF>
!   <NETCDF NAME="del_h" UNITS="none">
!    ref height interp factor for heat 
!   </NETCDF>
!   <NETCDF NAME="del_m" UNITS="none">
!    ref height interp factor for momentum 
!   </NETCDF>
!   <NETCDF NAME="del_q" UNITS="none">
!    ref height interp factor for moisture
!   </NETCDF>
!   <NETCDF NAME="tau_x" UNITS="pa">
!    zonal wind stress
!   </NETCDF>
!   <NETCDF NAME="tau_y" UNITS="pa">
!    meridional wind stress
!   </NETCDF>
!   <NETCDF NAME="t_surf" UNITS="deg_k">
!     surface temperature
!   </NETCDF>
!   <NETCDF NAME="t_ca" UNITS="deg_k">
!     canopy air temperature
!   </NETCDF>
!   <NETCDF NAME="q_surf" UNITS="kg/kg">
!     surface specific humidity 
!   </NETCDF>
!   <NETCDF NAME="shflx" UNITS="w/m2">
!     sensible heat flux
!   </NETCDF>
!   <NETCDF NAME="evap" UNITS="kg/m2/s">
!     evaporation rate 
!   </NETCDF>
!   <NETCDF NAME="lwflx" UNITS="w/m2">
!    net (down-up) longwave flux 
!   </NETCDF>

! </DIAGFIELDS>

! <INFO>


!   <NOTE>
!   <PRE>
!
!  MAIN PROGRAM EXAMPLE
!  --------------------
!
!       DO slow time steps (ocean)
!
!           call ICE_SLOW_UP
!
!           DO fast time steps (atmos)
!
!                call sfc_boundary_layer
!
!                call ATMOS_DOWN
!
!                call flux_down_from_atmos
!
!                call LAND_FAST
!
!                call ICE_FAST
!
!                call flux_up_to_atmos
!
!                call ATMOS_UP
!
!           END DO
!
!           call ICE_SLOW_DN
!
!           call OCEAN
!
!      END DO
!
!   LAND_FAST and ICE_FAST must update the surface temperature
!
! =======================================================================
!
! REQUIRED VARIABLES IN DEFINED DATA TYPES FOR COMPONENT MODELS
! --------------------------------------------------------------
!
! type (atmos_boundary_data_type) :: Atm
! type (surf_diff_type) :: Atm%Surf_Diff
!
! real, dimension(:)
!
!    Atm%lon_bnd   longitude axis grid box boundaries in radians
!                  must be monotonic
!    Atm%lat_bnd   latitude axis grid box boundaries in radians
!                  must be monotonic
!
! real, dimension(:,:)
!
!    Atm%t_bot     temperature at lowest model level
!    Atm%q_bot     specific humidity at lowest model level
!    Atm%z_bot     height above the surface for the lowest model level (m)
!    Atm%p_bot     pressure at lowest model level (pa)
!    Atm%u_bot     zonal wind component at lowest model level (m/s)
!    Atm%v_bot     meridional wind component at lowest model level (m/s)
!    Atm%p_surf    surface pressure (pa)
!    Atm%slp       sea level pressure (pa)
!    Atm%gust      gustiness factor (m/s)
!    Atm%flux_sw   net shortwave flux at the surface
!    Atm%flux_lw   downward longwave flux at the surface
!    Atm%lprec     liquid precipitation (kg/m2)
!    Atm%fprec     water equivalent frozen precipitation (kg/m2)
!    Atm%coszen    cosine of the zenith angle
!
!   (the following five fields are gathered into a data type for convenience in passing
!   this information through the different levels of the atmospheric model --
!   these fields are rlated to the simultaneous implicit time steps in the
!   atmosphere and surface models -- they are described more fully in
!   flux_exchange.tech.ps and
!   in the documntation for vert_diff_mod
!
!
!    Atm%Surf_Diff%dtmass   = dt/mass where dt = atmospheric time step ((i+1) = (i-1) for leapfrog) (s)
!                           mass = mass per unit area of lowest atmosphehic layer  (Kg/m2))
!    Atm%Surf_Diff%delta_t  increment ((i+1) = (i-1) for leapfrog) in temperature of
!                           lowest atmospheric layer  (K)
!    Atm%Surf_Diff%delta_q  increment ((i+1) = (i-1) for leapfrog) in specific humidity of
!                           lowest atmospheric layer (nondimensional -- Kg/Kg)
!    Atm%Surf_Diff%dflux_t  derivative of implicit part of downward temperature flux at top of lowest
!                           atmospheric layer with respect to temperature
!                           of lowest atmospheric layer (Kg/(m2 s))
!    Atm%Surf_Diff%dflux_q  derivative of implicit part of downward moisture flux at top of lowest
!                           atmospheric layer with respect to specific humidity of
!                           of lowest atmospheric layer (Kg/(m2 s))
!
!
! integer, dimension(4)
!
!    Atm%axes      Axis identifiers returned by diag_axis_init for the
!                  atmospheric model axes: X, Y, Z_full, Z_half.
!
! -----------------------------------------------
!
! type (land_boundary_data_type) :: Land
!
! real, dimension(:)
!
!    Land%lon_bnd     longitude axis grid box boundaries in radians
!                     must be monotonic
!    Land%lat_bnd     latitude axis grid box boundaries in radians
!                     must be monotonic
!
! logical, dimension(:,:,:)
!
!    Land%mask        land/sea mask (true for land)
!    Land%glacier     glacier mask  (true for glacier)
!
! real, dimension(:,:,:)
!
!    Land%tile_size   fractional area of each tile (partition)
!
!    Land%t_surf      surface temperature (deg k)
!    Land%albedo      surface albedo (fraction)
!    Land%rough_mom   surface roughness for momentum (m)
!    Land%rough_heat  surface roughness for heat/moisture (m)
!    Land%stomatal    stomatal resistance
!    Land%snow        snow depth (water equivalent) (kg/m2)
!    Land%water       water depth of the uppermost bucket (kg/m2)
!    Land%max_water   maximum water depth allowed in the uppermost bucket (kg/m2)
!
! -----------------------------------------------
!
!   </PRE>
!   </NOTE>
! </INFO>

end module flux_exchange_mod
