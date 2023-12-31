! ============================================================================
! canopy air
! ============================================================================
module canopy_air_mod

#include "../shared/debug.inc"


use utilities_mod, only : logunit, get_unit, error_mesg, FATAL, NOTE, check_nml_error
use time_manager_mod, only : time_type, time_type_to_real
use constants_mod, only : rdgas, rvgas, cp_air, PI, VONKARM
use sphum_mod, only : qscomp

use nf_utils_mod, only : nfu_inq_var
use land_constants_mod, only : NBANDS,d608,mol_CO2,mol_air
use cana_tile_mod, only : cana_tile_type, cana_prog_type, &
     canopy_air_mass, canopy_air_mass_for_tracers, cpw
use land_tile_mod, only : land_tile_type, land_tile_enum_type, &
     first_elmt, tail_elmt, next_elmt, current_tile, operator(/=)
use land_data_mod,      only : land_state_type, lnd
use land_tile_io_mod, only : create_tile_out_file, read_tile_data_r0d_fptr, write_tile_data_r0d_fptr, &
     get_input_restart_name, print_netcdf_error
use land_debug_mod, only : is_watch_point, check_temp_range

implicit none
private

! ==== public interfaces =====================================================
public :: read_cana_namelist
public :: cana_init
public :: cana_end
public :: save_cana_restart
public :: cana_turbulence
public :: cana_roughness
public :: cana_state
public :: cana_step_1
public :: cana_step_2
! ==== end of public interfaces ==============================================

! ==== module constants ======================================================
character(len=*), private, parameter :: &
  version = '$Id: canopy_air.F90,v 1.1.2.7 2012/06/19 18:34:53 pjp Exp $', &
  tagname = '$Name: no_fms_b_pjp $', &
  module_name = 'canopy_air_mod'

! options for turbulence parameter calculations
integer, parameter :: TURB_LM3W = 1, TURB_LM3V = 2
! ==== module variables ======================================================

!---- namelist ---------------------------------------------------------------
real :: init_T           = 288.
real :: init_T_cold      = 260.
real :: init_q           = 0.
real :: init_co2         = 350.0e-6 ! ppmv = mol co2/mol of dry air
real :: rav_lit_vi       = 0.       ! litter resistance to vapor per v_idx
character(len=32) :: turbulence_to_use = 'lm3w' ! or lm3v
logical :: save_qco2     = .TRUE.
logical :: sfc_dir_albedo_bug = .FALSE. ! if true, reverts to buggy behavior
! where direct albedo was mistakenly used for part of sub-canopy diffuse light
logical :: bare_ground_z0_bug = .FALSE. ! bug switch for bare z0m,d: if TRUE,
! z0m can't be less then "height", and that has lower limit of 0.1 m
namelist /cana_nml/ &
  init_T, init_T_cold, init_q, init_co2, turbulence_to_use, &
  canopy_air_mass, canopy_air_mass_for_tracers, cpw, rav_lit_vi, save_qco2, &
  sfc_dir_albedo_bug, bare_ground_z0_bug
!---- end of namelist --------------------------------------------------------

logical            :: module_is_initialized =.FALSE.
type(time_type)    :: time ! *** NOT YET USED
real               :: delta_time      ! fast time step
integer :: turbulence_option ! selected option of turbulence parameters 
     ! calculations

! ==== NetCDF declarations ===================================================
include 'netcdf.inc'
#define __NF_ASRT__(x) call print_netcdf_error((x),module_name,__LINE__)


contains


! ============================================================================
subroutine read_cana_namelist()
  ! ---- local vars
  integer :: io           ! i/o status for the namelist
  integer :: nml_unit

  write (logunit,'(/,80("="),/(a))') trim(version), trim(tagname)
  nml_unit = get_unit()
  open(nml_unit, file='input.nml', form='formatted', action='read', status='old')
  do 
     read (nml_unit, nml=cana_nml, iostat=io, end=10)
     if (check_nml_error (io, 'cana_nml')==0) exit ! from loop
  enddo
10  close (nml_unit)
  write (logunit, nml=cana_nml)
end subroutine read_cana_namelist

! ============================================================================
! initialize canopy air
subroutine cana_init ( id_lon, id_lat )
  integer, intent(in)          :: id_lon  ! ID of land longitude (X) axis  
  integer, intent(in)          :: id_lat  ! ID of land latitude (Y) axis

  ! ---- local vars ----------------------------------------------------------
  integer :: unit         ! unit for various i/o
  type(land_tile_enum_type)     :: te,ce ! last and current tile
  type(land_tile_type), pointer :: tile   ! pointer to current tile
  character(len=256) :: restart_file_name
  logical :: restart_exists

  module_is_initialized = .TRUE.

  ! ---- make module copy of time --------------------------------------------
  time       = lnd%time
  delta_time = time_type_to_real(lnd%dt_fast)

  ! ---- initialize cana state -----------------------------------------------
  ! first, set the initial values
  te = tail_elmt (lnd%tile_map)
  ce = first_elmt(lnd%tile_map)
  do while(ce /= te)
     tile=>current_tile(ce)  ! get pointer to current tile
     ce=next_elmt(ce)       ! advance position to the next tile
     
     if (.not.associated(tile%cana)) cycle
     
     if (associated(tile%glac)) then
        tile%cana%prog%T = init_T_cold
     else
        tile%cana%prog%T = init_T
     endif
     tile%cana%prog%q = init_q
     ! convert to kg CO2/kg wet air
     tile%cana%prog%co2 = init_co2*mol_CO2/mol_air*(1-tile%cana%prog%q) 
  enddo

  ! then read the restart if it exists
  call get_input_restart_name('INPUT/cana.res.nc',restart_exists,restart_file_name)
  if (restart_exists) then
     call error_mesg('cana_init',&
          'reading NetCDF restart "'//trim(restart_file_name)//'"',&
          NOTE)
     __NF_ASRT__(nf_open(restart_file_name,NF_NOWRITE,unit))
     call read_tile_data_r0d_fptr(unit, 'temp'  , cana_T_ptr  )
     call read_tile_data_r0d_fptr(unit, 'sphum' , cana_q_ptr )
     if(nfu_inq_var(unit,'co2')==NF_NOERR) then
        call read_tile_data_r0d_fptr(unit, 'co2', cana_co2_ptr )
     endif
     __NF_ASRT__(nf_close(unit))     
  else
     call error_mesg('cana_init',&
          'cold-starting canopy air',&
          NOTE)
  endif

  ! initialize options, to avoid expensive character comparisons during 
  ! run-time
  if (trim(turbulence_to_use)=='lm3v') then
     turbulence_option = TURB_LM3V
  else if (trim(turbulence_to_use)=='lm3w') then
     turbulence_option = TURB_LM3W
  else
     call error_mesg('cana_init', 'canopy air turbulence option turbulence_to_use="'// &
          trim(turbulence_to_use)//'" is invalid, use "lm3w" or "lm3v"', FATAL)
  endif
  
end subroutine cana_init


! ============================================================================
! release memory
subroutine cana_end ()

  module_is_initialized =.FALSE.

end subroutine cana_end


! ============================================================================
subroutine save_cana_restart (tile_dim_length, timestamp)
  integer, intent(in) :: tile_dim_length ! length of tile dim. in the output file
  character(*), intent(in) :: timestamp ! timestamp to add to the file name

  ! ---- local vars ----------------------------------------------------------
  integer :: unit            ! restart file i/o unit

!  call error_mesg('cana_end','writing NetCDF restart',NOTE)
  call create_tile_out_file(unit,'RESTART/'//trim(timestamp)//'cana.res.nc',&
          lnd%coord_glon, lnd%coord_glat, cana_tile_exists, tile_dim_length)

     ! write fields
     call write_tile_data_r0d_fptr(unit,'temp' ,cana_T_ptr,'canopy air temperature','degrees_K')
     call write_tile_data_r0d_fptr(unit,'sphum',cana_q_ptr,'canopy air specific humidity','kg/kg')
     if (save_qco2) then
        call write_tile_data_r0d_fptr(unit,'co2'  ,cana_co2_ptr,'canopy air co2 concentration','(kg CO2)/(kg wet air)')
     endif
     ! close output file
     __NF_ASRT__(nf_close(unit))
end subroutine save_cana_restart


! ============================================================================
subroutine cana_turbulence (u_star,&
     vegn_cover, vegn_layerfrac, vegn_height, vegn_bottom, vegn_lai, vegn_sai, vegn_d_leaf, &
     land_d, land_z0m, land_z0s, grnd_z0s, &
     con_v_h, con_v_v, con_g_h, con_g_v )
  real, intent(in) :: &
       u_star, & ! friction velocity, m/s
       land_d, land_z0m, land_z0s, grnd_z0s, & 
       vegn_cover, vegn_height(:), vegn_layerfrac(:), &
       vegn_bottom(:), & ! height of the bottom of the canopy, m
       vegn_lai(:), vegn_sai(:), vegn_d_leaf(:)
  real, intent(out) :: &
       con_v_h(:), con_v_v(:), & ! one-sided foliage-CAS conductance per unit ground area
       con_g_h   , con_g_v       ! ground-CAS conductance per unit ground area

  !---- local constants
  real, parameter :: a_max = 3
  real, parameter :: leaf_co = 0.01 ! meters per second^(1/2)
                                    ! leaf_co = g_b(z)/sqrt(wind(z)/d_leaf)
  real, parameter :: min_thickness = 0.01 ! thickness for switching to thin-canopy approximation, m
  real, parameter :: min_height = 0.1 ! min height of the canopy in TURB_LM3V case, m
  ! ---- local vars 
  real :: a        ! parameter of exponential wind profile within canopy:
                   ! u = u(ztop)*exp(-a*(1-z/ztop))
  real :: height   ! height of the current vegetation cohort, m
  real :: ztop     ! height of the tallest vegetation, m
  real :: wind     ! normalized wind on top of canopy, m/s
  real :: Kh_top   ! turbulent exchange coefficient on top of the canopy
  real :: vegn_idx ! total vegetation index = LAI+SAI, sum over cohorts
  real :: rah_sca  ! ground-SCA resistance
  real :: rav_lit  ! additional resistance of litter to vapor transport
  real :: h0       ! height of the canopy bottom, m

  integer :: i

  ! TODO: check array sizes

  vegn_idx = sum((vegn_lai+vegn_sai)*vegn_layerfrac)  ! total vegetation index
  
  select case(turbulence_option)
  case(TURB_LM3W)
     if(vegn_cover > 0) then
        ztop   = vegn_height(1)
        
        wind  = u_star/VONKARM*log((ztop-land_d)/land_z0m) ! normalized wind on top of the canopy
        a     = vegn_cover*a_max
        do i = 1,size(vegn_lai)
           height = vegn_height(i) ! effective height of the vegetation
           h0     = vegn_bottom(i) ! height of the bottom of the canopy
           if(height-h0>min_thickness) then
              con_v_h(i) = 2*vegn_lai(i)*leaf_co*sqrt(wind/vegn_d_leaf(i))*ztop/(height-h0)&
                 *(exp(-a/2*(ztop-height)/ztop)-exp(-a/2*(ztop-h0)/ztop))/a
           else
              ! thin cohort canopy limit
              con_v_h(i) = vegn_lai(i)*leaf_co*sqrt(wind/vegn_d_leaf(i))&
                 *exp(-a/2*(ztop-height)/ztop)
           endif
        enddo
        con_g_h = u_star*a*VONKARM*(1-land_d/ztop) &
             / (exp(a*(1-grnd_z0s/ztop)) - exp(a*(1-(land_z0s+land_d)/ztop)))
     else
        con_v_h = 0
        con_g_h = 0
     endif
  case(TURB_LM3V)
     ztop = max(vegn_height(1),min_height)
     
     a = a_max
     wind=u_star/VONKARM*log((ztop-land_d)/land_z0m) ! normalized wind on top of the canopy
  
     do i = 1,size(vegn_lai)
        height = max(vegn_height(i),min_height) ! effective height of the vegetation
        h0     = vegn_bottom(i) ! height of the canopy bottom above ground
        if(height-h0>min_thickness) then
           con_v_h(i) = 2*vegn_lai(i)*leaf_co*sqrt(wind/vegn_d_leaf(i))*ztop/(height-h0)&
              *(exp(-a/2*(ztop-height)/ztop)-exp(-a/2*(ztop-h0)/ztop))/a
        else
           ! thin cohort canopy limit
           con_v_h(i) = vegn_lai(i)*leaf_co*sqrt(wind/vegn_d_leaf(i))&
              *exp(-a/2*(ztop-height)/ztop)
        endif
     enddo

     if (land_d > 0.06 .and. vegn_idx > 0.25) then
        Kh_top = VONKARM*u_star*(ztop-land_d)
        rah_sca = ztop/a/Kh_top * &
             (exp(a*(1-grnd_z0s/ztop)) - exp(a*(1-(land_z0m+land_d)/ztop)))
        rah_sca = min(rah_sca,1250.0)
     else
        rah_sca=0.01
     endif
     con_g_h = 1.0/rah_sca
  end select
! not a good parameterization, but just using for sensitivity analyses now.
! ignores differing biomass and litter turnover rates.
  rav_lit = rav_lit_vi * vegn_idx
  con_g_v = con_g_h/(1.+rav_lit*con_g_h)
  con_v_v = con_v_h  
end subroutine cana_turbulence

! ============================================================================
! update effective surface roughness lengths for CAS-to-atmosphere fluxes
! and conductances for canopy-to-CAS and ground-to-CAS fluxes
!
! Strategy: Always define a canopy present. Non-vegetated situation is simply
! a limit as vegetation density approaches (but isn't allowed to reach) zero.
! Create expressions for the outputs that reduce to the special
! cases of full canopy cover and no canopy. Full canopy solution is that
! from Bonan (NCAR/TN-417+STR, 1996, p. 63). Thus, setting cover=1 in
! recovers Bonan. Letting cover approach 0 makes con_v_coef go to zero,
! preventing exchange with canopy, and makes con_g_coef go infinite,
! removing sub-canopy resistance and putting all resistance above the
! canopy, where it can be affected by stability adjustments.
!
! ** However, there is still a problem with this formulation when a
! canopy is present, because surface flux (I think) is not told to
! subtract out the resistances associated with con_v_coef and con_g_coef,
! which thus seem to be double-counted. For testing LM2, we should set them
! to zero anyway.
subroutine cana_roughness(lm2, &
     subs_z0m, subs_z0s, &
     snow_z0m, snow_z0s, snow_area, &
     vegn_cover, vegn_height, vegn_lai, vegn_sai, &
     land_d, land_z0m, land_z0s )
  logical, intent(in) :: lm2
  real, intent(in) :: &
       subs_z0m, subs_z0s, snow_z0m, snow_z0s, snow_area, vegn_cover, vegn_height, &
       vegn_lai, vegn_sai
  real, intent(out) :: &
       land_d    ,&
       land_z0m  ,&
       land_z0s

  !---- local constants
  real, parameter :: d_h_max = 2./3.
  real, parameter :: z0m_h_max = 1/7.35

  ! ---- local vars 
  real :: d_h      ! ratio of displacement height to vegetation height
  real :: z0m_h    ! ratio of roughness length to vegetation height
  real :: grnd_z0m, grnd_z0s
  real :: z0s_h, z0s_h_max
  real :: vegn_idx ! total vegetation index = LAI+SAI
  real :: height   ! effective vegetation height

  grnd_z0m = exp( (1-snow_area)*log(subs_z0m) + snow_area*log(snow_z0m))
  grnd_z0s = exp( (1-snow_area)*log(subs_z0s) + snow_area*log(snow_z0s))

  select case(turbulence_option)
  case(TURB_LM3W)
     if(vegn_cover > 0) then
        z0s_h_max = z0m_h_max*grnd_z0s/grnd_z0m ! to ensure cover->0 limit works
        d_h = vegn_cover*d_h_max
        if (lm2) then
           if (vegn_lai.gt.1) then  ! TEMP ***
              z0m_h = z0m_h_max
              z0s_h = z0s_h_max
           else
              z0m_h = grnd_z0m/vegn_height
              z0s_h = grnd_z0s/vegn_height
           endif
        else
           z0m_h = exp( vegn_cover*log(z0m_h_max) + (1-vegn_cover)*log(grnd_z0m/vegn_height))
           z0s_h = exp( vegn_cover*log(z0s_h_max) + (1-vegn_cover)*log(grnd_z0s/vegn_height))
        endif
        land_d   = d_h*vegn_height
        land_z0m = z0m_h*vegn_height
        land_z0s = z0s_h*vegn_height
     else
        land_d   = 0
        land_z0m = grnd_z0m
        land_z0s = grnd_z0s
     endif
     
  case(TURB_LM3V)
     height = max(vegn_height,0.1) ! effective height of the vegetation
     vegn_idx = vegn_lai+vegn_sai  ! total vegetation index
     if(vegn_idx>1e-4) then
        land_d = 1.1*height*log(1+(0.07*vegn_idx)**0.25)
        if(vegn_idx>2.85) then
           land_z0m = 0.3*(height-land_d)
        else
           land_z0m = grnd_z0m + 0.3*height*sqrt(0.07*vegn_idx)
        endif
     else 
        if (bare_ground_z0_bug) then
           ! bare soil or leaf off
           land_z0m = 0.1 *height
           land_d   = 0.66*height
        else
           land_d   = 0
           land_z0m = grnd_z0m
        endif
     endif
     land_z0s = land_z0m*exp(-2.0) 

  end select

end subroutine cana_roughness

! ============================================================================
subroutine cana_state ( cana, cana_T, cana_q, cana_co2 )
  type(cana_tile_type), intent(in)  :: cana
  real, optional      , intent(out) :: cana_T, cana_q, cana_co2

  if (present(cana_T))   cana_T   = cana%prog%T
  if (present(cana_q))   cana_q   = cana%prog%q
  if (present(cana_co2)) cana_co2 = cana%prog%co2
  
end subroutine

! ============================================================================
subroutine cana_step_1 ( cana,&
     p_surf, con_g_h, con_g_v, grnd_T, grnd_rh, grnd_rh_psi, &
     Hge,  DHgDTg, DHgDTc,    &
     Ege,  DEgDTg, DEgDqc, DEgDpsig     )
  type(cana_tile_type), intent(in) :: cana
  real, intent(in) :: &
     p_surf,  & ! surface pressure, Pa
     con_g_h, & ! conductivity between ground and CAS for heat
     con_g_v, & ! conductivity between ground and CAS for vapor
     grnd_T,  & ! ground temperature, degK
     grnd_rh, & ! ground relative humidity
     grnd_rh_psi ! psi derivative of ground relative humidity
  real, intent(out) ::   &
     Hge,  DHgDTg, DHgDTc, & ! linearization of the sensible heat flux from ground
     Ege,  DEgDTg, DEgDqc, DEgDpsig    ! linearization of evaporation from ground

  ! ---- local vars
  real :: rho, grnd_q, qsat, DqsatDTg

  call check_temp_range(grnd_T,'cana_step_1','grnd_T',lnd%time)

  call qscomp(grnd_T,p_surf,qsat,DqsatDTg)
  grnd_q = grnd_rh * qsat

  rho      =  p_surf/(rdgas*cana%prog%T*(1+d608*cana%prog%q))
  Hge      =  rho*cp_air*con_g_h*(grnd_T - cana%prog%T)
  DHgDTg   =  rho*cp_air*con_g_h
  DHgDTc   = -rho*cp_air*con_g_h
  Ege      =  rho*con_g_v*(grnd_q  - cana%prog%q)
  DEgDTg   =  rho*con_g_v*DqsatDTg*grnd_rh
  DEgDqc   = -rho*con_g_v
  DEgDpsig =  rho*con_g_v*qsat*grnd_rh_psi
  if(is_watch_point())then
     write(*,*)'#### cana_step_1 input ####'
     __DEBUG1__(p_surf)
     __DEBUG2__(con_g_h,con_g_v)
     __DEBUG2__(grnd_T,grnd_rh)
     write(*,*)'#### cana_step_1 internals ####'
     __DEBUG4__(rho, grnd_q, qsat, DqsatDTg)
     __DEBUG2__(cana%prog%T,cana%prog%q)
     write(*,*)'#### cana_step_1 output ####'
     __DEBUG3__(Hge,  DHgDTg, DHgDTc)
     __DEBUG4__(Ege,  DEgDTg, DEgDqc, DEgDpsig)
  endif
end subroutine 


! ============================================================================
subroutine cana_step_2 ( cana, delta_Tc, delta_qc )
  type(cana_tile_type), intent(inout) :: cana
  real, intent(in) ::  &
     delta_Tc, & ! change in canopy air temperature
     delta_qc    ! change in canopy air humidity

  cana%prog%T = cana%prog%T + delta_Tc
  cana%prog%q = cana%prog%q + delta_qc
end subroutine cana_step_2

! ============================================================================
! tile existence detector: returns a logical value indicating wether component
! model tile exists or not
logical function cana_tile_exists(tile)
   type(land_tile_type), pointer :: tile
   cana_tile_exists = associated(tile%cana)
end function cana_tile_exists

! ============================================================================
! accessor functions: given a pointer to a land tile, they return pointer
! to the desired member of the land tile, of NULL if this member does not
! exist.
#define DEFINE_CANA_ACCESSOR_0D(xtype,x) subroutine cana_ ## x ## _ptr(t,p);\
type(land_tile_type),pointer::t;xtype,pointer::p;p=>NULL();if(associated(t))then;if(associated(t%cana))p=>t%cana%prog%x;endif;end subroutine

DEFINE_CANA_ACCESSOR_0D(real,T)
DEFINE_CANA_ACCESSOR_0D(real,q)
DEFINE_CANA_ACCESSOR_0D(real,co2)

end module canopy_air_mod
