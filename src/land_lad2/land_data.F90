module land_data_mod

use utilities_mod     , only : num_land_tracers, land_tracer_names, logunit
use constants_mod     , only : PI
use time_manager_mod  , only : time_type
use land_tile_mod     , only : land_tile_type, land_tile_list_type, &
     land_tile_list_init, land_tile_list_end, nitems

implicit none
private

! ==== public interfaces =====================================================
public :: land_data_init
public :: land_data_end
public :: lnd            ! global data 

public :: atmos_land_boundary_type ! container for information passed from the 
                         ! atmosphere to land
public :: land_data_type ! container for information passed from land to 
                         ! the atmosphere
! both hold information on land grid (that is, after flux exchange translated 
! it from the atmosphere)

public :: dealloc_land2cplr ! deallocates a land_data_type structure
public :: realloc_land2cplr ! allocates a land_data_type members for current 
                            ! number of tiles
public :: dealloc_cplr2land ! deallocates an atmos_land_boundary_type structure
public :: realloc_cplr2land ! allocates an atmos_land_boundary_type members 
                            ! for current number of tiles
! NOTE: realloc_* procedures can be called regardless of the current state
! of the argument data structures, since they deallocate data first.

public :: land_state_type
! ==== end of public interfaces ==============================================

! ---- module constants ------------------------------------------------------
character(len=*), parameter :: &
     module_name = 'land_data_mod', &
     version     = '$Id: land_data.F90,v 1.1.2.9 2012/06/19 18:34:54 pjp Exp $', &
     tagname     = '$Name: no_fms_b_pjp $'

! init_value is used to fill most of the allocated boundary condition arrays.
! It is supposed to be double-precision signaling NaN, to triger a trap when
! the program is compiled with trapping uninitialized values.  
! See http://ftp.uniovi.es/~antonio/uned/ieee754/IEEE-754references.html
real, parameter :: init_value = Z'FFF0000000000001'

! ---- types -----------------------------------------------------------------
type :: atmos_land_boundary_type
   ! data passed from the coupler to the surface
   real, dimension(:,:,:), pointer :: & ! (lon, lat, tile)
        t_flux    => NULL(), &   ! sensible heat flux, W/m2
        lw_flux   => NULL(), &   ! net longwave radiation flux, W/m2
        lwdn_flux => NULL(), &   ! downward longwave radiation flux, W/m2
        sw_flux   => NULL(), &   ! net shortwave radiation flux, W/m2
        swdn_flux => NULL(), &   ! downward shortwave radiation flux, W/m2
        lprec     => NULL(), &   ! liquid precipitation rate, kg/(m2 s)
        fprec     => NULL(), &   ! frozen precipitation rate, kg/(m2 s)
        tprec     => NULL(), &   ! temperature of precipitation, degK
   ! components of downward shortwave flux, W/m2  
        sw_flux_down_vis_dir   => NULL(), & ! visible direct 
        sw_flux_down_total_dir => NULL(), & ! total direct
        sw_flux_down_vis_dif   => NULL(), & ! visible diffuse
        sw_flux_down_total_dif => NULL(), & ! total diffuse
   ! derivatives of the fluxes
        dhdt      => NULL(), &   ! sensible w.r.t. surface temperature
        dhdq      => NULL(), &   ! sensible w.r.t. surface humidity
        drdt      => NULL(), &   ! longwave w.r.t. surface radiative temperature 
   !
        cd_m      => NULL(), &   ! drag coefficient for momentum, dimensionless
        cd_t      => NULL(), &   ! drag coefficient for tracers, dimensionless
        ustar     => NULL(), &   ! turbulent wind scale, m/s
        bstar     => NULL(), &   ! turbulent buoyancy scale, m/s
        wind      => NULL(), &   ! abs wind speed at the bottom of the atmos, m/s
        z_bot     => NULL(), &   ! height of the bottom atmospheric layer above the surface, m
        drag_q    => NULL(), &   ! product of cd_q by wind
        p_surf    => NULL()      ! surface pressure, Pa

   real, dimension(:,:,:,:), pointer :: & ! (lon, lat, tile, tracer)
        tr_flux => NULL(),   &   ! tracer flux, including water vapor flux
        dfdtr   => NULL()        ! derivative of the flux w.r.t. tracer surface value, 
                                 ! including evap over surface specific humidity

   integer :: xtype             !REGRID, REDIST or DIRECT
end type atmos_land_boundary_type


type :: land_data_type
   ! data passed from the surface to the coupler
   logical :: pe ! data presence indicator for stock calculations
   real, pointer, dimension(:,:,:)   :: &  ! (lon, lat, tile)
        tile_size      => NULL(),  & ! fractional coverage of cell by tile, dimensionless
        t_surf         => NULL(),  & ! ground surface temperature, degK
        t_ca           => NULL(),  & ! canopy air temperature, degK
        albedo         => NULL(),  & ! broadband land albedo [unused?]
        albedo_vis_dir => NULL(),  & ! albedo for direct visible radiation
        albedo_nir_dir => NULL(),  & ! albedo for direct NIR radiation 
        albedo_vis_dif => NULL(),  & ! albedo for diffuse visible radiation 
        albedo_nir_dif => NULL(),  & ! albedo for diffuse NIR radiation
        rough_mom      => NULL(),  & ! surface roughness length for momentum, m
        rough_heat     => NULL(),  & ! roughness length for tracers and heat, m
        rough_scale    => NULL()     ! topographic scaler for momentum drag, m

   real, pointer, dimension(:,:,:,:)   :: &  ! (lon, lat, tile, tracer)
        tr    => NULL()              ! tracers, including canopy air specific humidity

   ! NOTE that in contrast to most of the other fields in this structure, the discharges
   ! hold data per-gridcell, rather than per-tile basis. This, and the order of updates,
   ! have implications for the data reallocation procedure.
   real, pointer, dimension(:,:) :: &  ! (lon, lat)
     discharge           => NULL(),  & ! liquid water flux from land to ocean
     discharge_heat      => NULL(),  & ! sensible heat of discharge (0 C datum)
     discharge_snow      => NULL(),  & ! solid water flux from land to ocean
     discharge_snow_heat => NULL()     ! sensible heat of discharge_snow (0 C datum)

   logical, pointer, dimension(:,:,:):: &
        mask => NULL()               ! true if land

   integer :: axes(2)        ! IDs of diagnostic axes
   logical, pointer :: maskmap(:,:) 
end type land_data_type


! land_state_type combines the general information about state of the land model:
! coordinates, time steps, etc. There is only one variable of this type,
! and it is public in this module.
type :: land_state_type
   integer        :: is,ie,js,je
   integer        :: nlon,nlat
   integer        :: ntprog      ! number of prognostic tracers
   integer        :: isphum=0    ! index of specific humidity in tracer array
   integer        :: ico2=0      ! index of carbon dioxide in tracer array
   type(time_type):: dt_fast     ! fast (physical) time step
   type(time_type):: dt_slow     ! slow time step
   type(time_type):: time        ! current time

   real, pointer  :: lon (:,:), lat (:,:) ! grid center coordinates, radian
   real, pointer  :: lonb(:,:), latb(:,:) ! grid vertices, radian
   real, pointer  :: area(:,:)      ! land area per grid cell, m2
   real, pointer  :: cellarea(:,:)  ! grid cell area, m2
   real, pointer  :: coord_glon(:), coord_glonb(:) ! longitudes, degrees East
   real, pointer  :: coord_glat(:), coord_glatb(:) ! latitudes, degrees North

   ! map of tiles
   type(land_tile_list_type), pointer :: tile_map(:,:)

   integer :: nfaces ! number of mosaic faces
   integer :: face  ! the current mosaic face
   integer, allocatable :: pelist(:) ! list of processors that run land model
   integer, allocatable :: io_pelist(:)
   integer :: io_id     ! suffix in the distributed files.
end type land_state_type

! ---- public module variables -----------------------------------------------
type(land_state_type),save :: lnd


! ---- private module variables ----------------------------------------------
logical :: module_is_initialized =.FALSE.


#define __DEALLOC__(x) if (associated(x)) deallocate(x)

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-



! ============================================================================
subroutine land_data_init(layout, io_layout, time, dt_fast, dt_slow)
  integer, intent(inout) :: layout(2)
  integer, intent(inout) :: io_layout(2) ! layout for land model io
  type(time_type), intent(in) :: &
       time,    & ! current model time
       dt_fast, & ! fast (physical) time step
       dt_slow    ! slow time step

  ! ---- local vars
  integer :: nlon, nlat ! size of global grid in lon and lat directions
  integer :: ntiles     ! number of tiles in the mosaic grid 
  integer :: i,j,tr
  integer :: n_io_pes(2) ! number of PEs along x and y
  integer :: io_id(1)

  ! write the version and tag name to the logfile
  write (logunit,'(/,80("="),/(a))') trim(version), trim(tagname)

  ! define the processor layout information according to the global grid size 
  ntiles = 1
  lnd%nfaces = ntiles
  nlon = 1
  nlat = 1
  lnd%nlon = nlon
  lnd%nlat = nlat
  layout(1) = 1
  layout(2) = 1
  io_layout(1) = 1
  io_layout(2) = 1

  allocate(lnd%io_pelist(1))
  lnd%io_pelist(1) = 0
  lnd%io_id        = 0
  lnd%is=1; lnd%ie=1; lnd%js=1; lnd%je=1
  lnd%face = 1

  allocate(lnd%tile_map(lnd%is:lnd%ie,   lnd%js:lnd%je))
  allocate(lnd%lonb    (lnd%is:lnd%ie+1, lnd%js:lnd%je+1))
  allocate(lnd%latb    (lnd%is:lnd%ie+1, lnd%js:lnd%je+1))
  allocate(lnd%lon     (lnd%is:lnd%ie,   lnd%js:lnd%je))
  allocate(lnd%lat     (lnd%is:lnd%ie,   lnd%js:lnd%je))
  allocate(lnd%area    (lnd%is:lnd%ie,   lnd%js:lnd%je))
  allocate(lnd%cellarea(lnd%is:lnd%ie,   lnd%js:lnd%je))
  allocate(lnd%coord_glon(nlon), lnd%coord_glonb(nlon+1))
  allocate(lnd%coord_glat(nlat), lnd%coord_glatb(nlat+1))

  ! initialize coordinates
  lnd%latb(1,1) =  35.40
  lnd%latb(2,1) =  35.40
  lnd%latb(1,2) =  36.40
  lnd%latb(2,2) =  36.40
  lnd%lonb(1,1) = 275.17
  lnd%lonb(1,2) = 275.17
  lnd%lonb(2,1) = 276.17
  lnd%lonb(2,2) = 276.17
  lnd%lon(1,1)  = 275.67
  lnd%lat(1,1)  =  35.90
  lnd%cellarea  = 1.0015480220222198E+10
  lnd%area      = 1.0015480220222198E+10
  lnd%coord_glon(1) = 275.67
  lnd%coord_glat(1) = 35.9
  lnd%coord_glonb(:) = (/275.17,276.17/)
  lnd%coord_glatb(:) = (/ 35.40, 36.40/)

  ! convert coordinates to radian
  lnd%lonb = lnd%lonb*pi/180.0 ; lnd%lon = lnd%lon*pi/180.0
  lnd%latb = lnd%latb*pi/180.0 ; lnd%lat = lnd%lat*pi/180.0

  ! initialize land tile map
  do j = lnd%js,lnd%je
  do i = lnd%is,lnd%ie
     call land_tile_list_init(lnd%tile_map(i,j))
  enddo
  enddo

  ! initialize land model tracers, if necessary
  do tr=1,num_land_tracers
    if(trim(land_tracer_names(tr)) == 'sphum') exit
  enddo
  lnd%isphum = tr
  lnd%ntprog = num_land_tracers

  do tr=1,num_land_tracers
    if(trim(land_tracer_names(tr)) == 'co2') exit
  enddo
  lnd%ico2 = tr

  ! initialize model's time-related parameters
  lnd%time    = time
  lnd%dt_fast = dt_fast
  lnd%dt_slow = dt_slow

  ! initialize the land model processor list
  allocate(lnd%pelist(0:0))
  lnd%pelist = 0

end subroutine land_data_init

! ============================================================================
subroutine land_data_end()

  integer :: i,j
  
  module_is_initialized = .FALSE.

  ! deallocate land tile map here. 
  do j = lnd%js,lnd%je
  do i = lnd%is,lnd%ie
     call land_tile_list_end(lnd%tile_map(i,j))
  enddo
  enddo

  ! deallocate grid data
  deallocate(lnd%lonb, lnd%latb, lnd%lon, lnd%lat, lnd%area, lnd%cellarea, &
       lnd%coord_glonb, lnd%coord_glon, &
       lnd%coord_glatb, lnd%coord_glat, &
       lnd%tile_map, lnd%pelist, lnd%io_pelist)

end subroutine land_data_end


! ============================================================================
subroutine realloc_land2cplr ( bnd )
  type(land_data_type), intent(inout) :: bnd     ! data to allocate

  ! ---- local vars
  integer :: n_tiles

  call dealloc_land2cplr(bnd, dealloc_discharges=.FALSE.)

  n_tiles = max_n_tiles()

  allocate( bnd%mask(lnd%is:lnd%ie,lnd%js:lnd%je,n_tiles) )

  allocate( bnd%tile_size(lnd%is:lnd%ie,lnd%js:lnd%je,n_tiles) )
  allocate( bnd%t_surf(lnd%is:lnd%ie,lnd%js:lnd%je,n_tiles) )
  allocate( bnd%t_ca(lnd%is:lnd%ie,lnd%js:lnd%je,n_tiles) )
  allocate( bnd%tr(lnd%is:lnd%ie,lnd%js:lnd%je,n_tiles,lnd%ntprog) )
  allocate( bnd%albedo(lnd%is:lnd%ie,lnd%js:lnd%je,n_tiles) )
  allocate( bnd%albedo_vis_dir(lnd%is:lnd%ie,lnd%js:lnd%je,n_tiles) )
  allocate( bnd%albedo_nir_dir(lnd%is:lnd%ie,lnd%js:lnd%je,n_tiles) )
  allocate( bnd%albedo_vis_dif(lnd%is:lnd%ie,lnd%js:lnd%je,n_tiles) )
  allocate( bnd%albedo_nir_dif(lnd%is:lnd%ie,lnd%js:lnd%je,n_tiles) )
  allocate( bnd%rough_mom(lnd%is:lnd%ie,lnd%js:lnd%je,n_tiles) )
  allocate( bnd%rough_heat(lnd%is:lnd%ie,lnd%js:lnd%je,n_tiles) )
  allocate( bnd%rough_scale(lnd%is:lnd%ie,lnd%js:lnd%je,n_tiles) )

  bnd%mask              = .FALSE.
  bnd%tile_size         = init_value
  bnd%t_surf            = init_value
  bnd%t_ca              = init_value
  bnd%tr                = init_value
  bnd%albedo            = init_value
  bnd%albedo_vis_dir    = init_value
  bnd%albedo_nir_dir    = init_value
  bnd%albedo_vis_dif    = init_value
  bnd%albedo_nir_dif    = init_value
  bnd%rough_mom         = init_value
  bnd%rough_heat        = init_value
  bnd%rough_scale       = init_value

  ! in contrast to the rest of the land boundary condition fields, discharges 
  ! are specified per grid cell, not per tile; therefore they should not be 
  ! re-allocated when the number of tiles changes. In fact, they must not be
  ! changed at all here because their values are assigned in update_land_model_fast,
  ! not in update_land_bc_*, and therefore would be lost if re-allocated.
  if (.not.associated(bnd%discharge)) then
     allocate( bnd%discharge          (lnd%is:lnd%ie,lnd%js:lnd%je) )
     allocate( bnd%discharge_heat     (lnd%is:lnd%ie,lnd%js:lnd%je) )
     allocate( bnd%discharge_snow     (lnd%is:lnd%ie,lnd%js:lnd%je) )
     allocate( bnd%discharge_snow_heat(lnd%is:lnd%ie,lnd%js:lnd%je) )

     ! discharge and discaharge_snow must be, in contrast to the rest of the boundary
     ! values, filled with zeroes. The reason is because not all of the usable elements
     ! are updated by the land model (only coastal points are).
     bnd%discharge           = 0.0
     bnd%discharge_heat      = 0.0
     bnd%discharge_snow      = 0.0
     bnd%discharge_snow_heat = 0.0
  endif
end subroutine realloc_land2cplr


! ============================================================================
! deallocates boundary data memory
! NOTE that the discharges should be deallocated only at the final clean-up
! stage; during the model run they should be preserved unchanged even when
! other fields are reallocated.
subroutine dealloc_land2cplr ( bnd, dealloc_discharges )
  type(land_data_type), intent(inout) :: bnd  ! data to de-allocate
  logical, intent(in) :: dealloc_discharges

  __DEALLOC__( bnd%tile_size )
  __DEALLOC__( bnd%tile_size )
  __DEALLOC__( bnd%t_surf )
  __DEALLOC__( bnd%t_ca )
  __DEALLOC__( bnd%tr )
  __DEALLOC__( bnd%albedo )
  __DEALLOC__( bnd%albedo_vis_dir )
  __DEALLOC__( bnd%albedo_nir_dir )
  __DEALLOC__( bnd%albedo_vis_dif )
  __DEALLOC__( bnd%albedo_nir_dif )
  __DEALLOC__( bnd%rough_mom )
  __DEALLOC__( bnd%rough_heat )
  __DEALLOC__( bnd%rough_scale )
  __DEALLOC__( bnd%mask )

  if (dealloc_discharges) then
     __DEALLOC__( bnd%discharge           )
     __DEALLOC__( bnd%discharge_heat      )
     __DEALLOC__( bnd%discharge_snow      )
     __DEALLOC__( bnd%discharge_snow_heat )
  end if

end subroutine dealloc_land2cplr


! ============================================================================
! initializes data for data override.
! NOTE: previously the body of the procedure was in the flux_exchange_init,
! currently it is called from land_model_init
subroutine realloc_cplr2land( bnd )
  type(atmos_land_boundary_type), intent(inout) :: bnd

  ! ---- local vars
  integer :: kd

  call dealloc_cplr2land(bnd)

  kd = max_n_tiles()

  allocate( bnd%t_flux(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )
  allocate( bnd%lw_flux(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )
  allocate( bnd%sw_flux(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )
  allocate( bnd%lprec(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )
  allocate( bnd%fprec(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )
  allocate( bnd%tprec(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )
  allocate( bnd%dhdt(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )
  allocate( bnd%dhdq(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )
  allocate( bnd%drdt(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )
  allocate( bnd%p_surf(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )
  allocate( bnd%tr_flux(lnd%is:lnd%ie,lnd%js:lnd%je,kd,lnd%ntprog) )
  allocate( bnd%dfdtr(lnd%is:lnd%ie,lnd%js:lnd%je,kd,lnd%ntprog) )

  allocate( bnd%lwdn_flux(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )
  allocate( bnd%swdn_flux(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )
  allocate( bnd%sw_flux_down_vis_dir(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )
  allocate( bnd%sw_flux_down_total_dir(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )
  allocate( bnd%sw_flux_down_vis_dif(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )
  allocate( bnd%sw_flux_down_total_dif(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )
  allocate( bnd%cd_t(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )
  allocate( bnd%cd_m(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )
  allocate( bnd%bstar(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )
  allocate( bnd%ustar(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )
  allocate( bnd%wind(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )
  allocate( bnd%z_bot(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )

  allocate( bnd%drag_q(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )

  bnd%t_flux                 = init_value
  bnd%lw_flux                = init_value
  bnd%sw_flux                = init_value
  bnd%lprec                  = init_value
  bnd%fprec                  = init_value
  bnd%tprec                  = init_value
  bnd%dhdt                   = init_value
  bnd%dhdq                   = init_value
  bnd%drdt                   = init_value
  bnd%p_surf                 = init_value
  bnd%tr_flux                = init_value
  bnd%dfdtr                  = init_value

  bnd%lwdn_flux              = init_value
  bnd%swdn_flux              = init_value
  bnd%sw_flux_down_vis_dir   = init_value
  bnd%sw_flux_down_total_dir = init_value
  bnd%sw_flux_down_vis_dif   = init_value
  bnd%sw_flux_down_total_dif = init_value
  bnd%cd_t                   = init_value
  bnd%cd_m                   = init_value
  bnd%bstar                  = init_value
  bnd%ustar                  = init_value
  bnd%wind                   = init_value
  bnd%z_bot                  = init_value

  bnd%drag_q                 = init_value

end subroutine realloc_cplr2land


! ============================================================================
subroutine dealloc_cplr2land( bnd )
  type(atmos_land_boundary_type), intent(inout) :: bnd

  __DEALLOC__( bnd%t_flux )
  __DEALLOC__( bnd%lw_flux )
  __DEALLOC__( bnd%sw_flux )
  __DEALLOC__( bnd%lprec )
  __DEALLOC__( bnd%fprec )
  __DEALLOC__( bnd%dhdt )
  __DEALLOC__( bnd%dhdq )
  __DEALLOC__( bnd%drdt )
  __DEALLOC__( bnd%p_surf )
  __DEALLOC__( bnd%lwdn_flux )
  __DEALLOC__( bnd%swdn_flux )
  __DEALLOC__( bnd%sw_flux_down_vis_dir )
  __DEALLOC__( bnd%sw_flux_down_total_dir )
  __DEALLOC__( bnd%sw_flux_down_vis_dif )
  __DEALLOC__( bnd%sw_flux_down_total_dif )
  __DEALLOC__( bnd%cd_t )
  __DEALLOC__( bnd%cd_m )
  __DEALLOC__( bnd%bstar )
  __DEALLOC__( bnd%ustar )
  __DEALLOC__( bnd%wind )
  __DEALLOC__( bnd%z_bot )
  __DEALLOC__( bnd%tr_flux )
  __DEALLOC__( bnd%dfdtr )

end subroutine dealloc_cplr2land

! ============================================================================
function max_n_tiles() result(n)
  integer :: n
  integer :: i,j

  n=1
  do j=lnd%js,lnd%je
  do i=lnd%is,lnd%ie
     n=max(n, nitems(lnd%tile_map(i,j)))
  enddo
  enddo

end function 

end module land_data_mod
