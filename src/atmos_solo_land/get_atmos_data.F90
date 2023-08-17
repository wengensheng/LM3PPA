! ============================================================================
! this module provides prescribed atmospheric forcing
! ============================================================================

#define __NF_ASRT__(err_code) call cdf_err(err_code, __FILE__, __LINE__, FATAL)

module prescr_forcing_mod

use utilities_mod,      only : logunit, get_unit, error_mesg, FATAL, NOTE, check_nml_error
use time_manager_mod,   only : time_type, get_date, get_time, set_time, &
     operator(-), operator(+), operator(/), print_date, set_date
use get_cal_time_mod,   only : get_cal_time
use time_interp_mod ,   only : time_interp, YEAR, NONE
use astronomy_mod,      only : diurnal_solar
use constants_mod,      only : PI

implicit none
private

#include "netcdf.inc"

! ==== public interfaces =====================================================
public:: prescr_forcing_type ! forcing information data type

public:: prescr_forcing_init ! initializes the forcing calculations
public:: prescr_forcing_end  ! finishes using the forcing

public:: get_prescr_forcing  ! calculates the forcing at a given time
! ==== end of public interfaces ==============================================


type :: prescr_forcing_type
   private
   type(time_type)       :: dt      ! time step of atmospheric model
end type prescr_forcing_type


! ---- module constants ------------------------------------------------------
character(len=*),   parameter :: mod_name = 'prescr_forcing'
character(len=128), parameter :: version  = '$Id: get_atmos_data.F90,v 1.1.2.1.2.1.2.5 2011/03/28 16:24:14 pjp Exp $'
character(len=128), parameter :: tagname  = '$Name: no_fms_b_pjp $'
integer, parameter :: &
     iprec  = 1,      &
     itemp  = 2,      &
     iq     = 3,      &
     iwind  = 4,      &
     ips    = 5,      &
     ilwdn  = 6,      &
     iswdn  = 7,      &
     itime  = 8

! ---- module types ----------------------------------------------------------
type :: field_type
   character(len=32) :: name    = ''
   real, pointer     :: data(:) => NULL() 
end type field_type

! ---- module variables ------------------------------------------------------
type(field_type) :: field(8)
data field(:)%name / &
     'prcp',       'tas',        'shum',        'wind',         'pres', &
     'dlwrf',      'dswrf',      'time'/

type(time_type), allocatable :: time_line(:)
type(time_type) :: ts, te
real :: lat, lon ! coordinates of our point

logical :: module_is_initialized = .false.

! ---- namelist --------------------------------------------------------------
character(256) :: filelist = ''            ! name of the file containing list of the data file names
character(10)  :: timeline  = 'normal'     ! type of the timeline ('normal' or 'loop')
integer, dimension(6)  :: &
     start_loop=(/1982,1,1,0,0,0/),  &     ! beginning of the time loop
     end_loop  =(/1999,1,1,0,0,0/)         ! end of the time loop
real    :: t_fprec = 273.16                ! precip is solid if temperature below this
logical :: use_diurnal_solar = .false.     ! if true, calculated diurnal cycle is 
   ! used for SW instead of linear interpolation
real    :: swdn_dir_part = 0.5 ! fraction of direct (beam) light in the sw downward flux
real    :: swdn_vis_part = 0.5 ! fraction of the visible (PAR) light in the sw downward flux
logical :: limit_swdn_by_TOA_flux = .false.! if true, scaled swdn is limited byt the
   ! calculated TOA flux
real    :: solar = 1360.0                  ! solar constant for calculations of TOA swdn, W/m2
real    :: min_lwdn = 0.0    ! lower bound for downward lw flux, W/m2

namelist /prescr_forcing_nml/ &
     use_diurnal_solar, & ! if true, calculated diurnal cycle is used for SW 
     t_fprec,   &  ! precip is solid if temperature below this threshold
     swdn_dir_part, & ! fraction of direct (beam) light in the sw downward flux
     swdn_vis_part, & ! fraction of the visible (PAR) light in the sw downward flux
     limit_swdn_by_TOA_flux, & ! if true, scaled swdn is limited by the calculated TOA flux
     solar, &      ! solar constant for calculations of TOA swdn, W/m2
     min_lwdn, & ! lower bound for the longwave downward radiation
     filelist,  &  ! name of the file containing the list of input data files
     timeline,  &  ! type of the timeline ('normal' or 'loop')
     start_loop, end_loop  ! beginning and end of the time loop

contains

! ============================================================================
! Assumptions:
! (1) all variables are on the same time line
! (2) all variables are in a single netcdf file
 subroutine prescr_forcing_init ( Forcing, glon, glat, Time, dt)
! initializes the forcing calculations
  type(prescr_forcing_type), intent(out) :: &
       Forcing                        ! data to initialize
  real, intent(in) :: glon(:,:)       ! longitudes of global model grid
  real, intent(in) :: glat(:,:)       ! latitudes of global model grid
  type(time_type), intent(in) :: Time ! current time
  type(time_type), intent(in) :: dt   ! fast time step of atmos model 
  
  ! ---- local vars ----------------------------------------------------------
  integer :: unit      ! unit number for i/o
  integer :: io, ierr  ! error codes for various purposes
  logical :: fexist
  integer :: i, j
  character(len=1024) :: filename ! name of the input files
  integer :: ncid, varid, timeid, unlimdim
  integer :: len, time_len, time_idx
  character(len=256) :: units    ! units ot time in the file
  character(len=256) :: calendar ! calendar of the data
  
  ! read namelist
  unit = get_unit()
  open(unit, file='input.nml', form='formatted', action='read', status='old')
  do 
     read (unit, nml=prescr_forcing_nml, iostat=io, end=10)
     if (check_nml_error (io, 'prescr_forcing_nml')==0) exit ! from loop
  enddo
10   close (unit)
  
  ! write the namelist to the log file
  write (logunit,'(/,80("="),/(a))') trim(version), trim(tagname)
  write (logunit, nml=prescr_forcing_nml)
  
  Forcing%dt     = dt

  ! initialize input data set
  ts = set_date(start_loop(1),start_loop(2),start_loop(3),start_loop(4),start_loop(5),start_loop(6))
  te = set_date(end_loop(1),end_loop(2),end_loop(3), end_loop(4),end_loop(5),end_loop(6))

  ! calculate the total size of the time line
  time_len = 0
  unit = get_unit()
  open(unit, file=trim(filelist), form='formatted', action='read')
  do while (.true.)
     filename=''
     read(unit,'(a)',end=100) filename
     inquire (file=trim(filename), exist=fexist)
     if (.not.fexist) then
       call error_mesg ( mod_name, 'Requested data file "'//trim(filename)// '" does not exist', FATAL)
     endif
     __NF_ASRT__(nf_open(trim(filename),NF_NOWRITE,ncid))
     __NF_ASRT__(nf_inq_unlimdim(ncid,unlimdim))
     __NF_ASRT__(nf_inq_dimlen(ncid,unlimdim,len))
     time_len = time_len + len
     __NF_ASRT__(nf_close(ncid))
  enddo
100 continue
  close(unit)

  ! allocate storage for the input data
  do i = 1, size(field)
     allocate(field(i)%data(time_len))
  enddo
  allocate(time_line(time_len))

  ! read the input data
  unit = get_unit()
  open(unit, file=trim(filelist), form='formatted', action='read')
  time_idx = 1
  do while (.true.)
     filename=''
     read(unit,'(a)',end=101) filename
     __NF_ASRT__(nf_open(trim(filename),NF_NOWRITE,ncid))
     
     __NF_ASRT__(nf_inq_varid(ncid,'longitude',varid))
     __NF_ASRT__(nf_get_var_double(ncid,varid,lon))
     __NF_ASRT__(nf_inq_varid(ncid,'latitude',varid))
     __NF_ASRT__(nf_get_var_double(ncid,varid,lat))

     do i = 1, size(field)
        __NF_ASRT__(nf_inq_varid(ncid,field(i)%name,varid))
        __NF_ASRT__(nf_get_var_double(ncid,varid,field(i)%data(time_idx)))
     enddo

     ! ADVANCE THE CURRENT INDEX IN THE TIMELINE
     __NF_ASRT__(nf_inq_unlimdim(ncid,unlimdim))
     __NF_ASRT__(nf_inq_dimlen(ncid,unlimdim,len))
     __NF_ASRT__(nf_inq_varid(ncid,field(itime)%name,timeid))

     ! GET UNITS OF THE TIME
     units = ' '
     __NF_ASRT__(nf_get_att_text(ncid,timeid,'units',units))

     ! GET CALENDAR OF THE DATA
     calendar = ' '
     ierr = nf_get_att_text(ncid,timeid,'calendar',calendar)
     if(ierr/=NF_NOERR) &
          ierr = nf_get_att_text(ncid,timeid,'calendar_type',calendar)
     if(ierr/=NF_NOERR) &
          calendar='JULIAN' ! use model calendar? how to get the name of the model calendar?

     ! CONVERT TIME
     do j = time_idx, time_idx+len-1
        time_line(j) = get_cal_time(field(itime)%data(j),units,calendar)
     enddo
     time_idx = time_idx+len

     __NF_ASRT__(nf_close(ncid))
  enddo
101 continue
  close(unit) 

  ! CONVERT COORDINATES TO RADIAN
  lon = lon*PI/180
  lat = lat*PI/180

  module_is_initialized = .true.

end subroutine prescr_forcing_init


! ============================================================================
subroutine prescr_forcing_end ( Forcing )
! ============================================================================
! terminates forcing data
  ! ---- arguments -----------------------------------------------------------
  type(prescr_forcing_type), intent(inout) :: Forcing

  module_is_initialized = .false.

end subroutine prescr_forcing_end

! ============================================================================
subroutine get_prescr_forcing ( Forcing, Time, &
     swdn_vis_dir, swdn_vis_dif, swdn_tot_dir, swdn_tot_dif, lwdn, &
     Ps, temp, sphum, U, V, lprec, fprec )
! ============================================================================
! calculates prescribed forcing for a given time
  ! ---- arguments -----------------------------------------------------------
  type(prescr_forcing_type), intent(inout) :: &
       Forcing                          ! data to use and possibly update
  type(time_type), intent(in) :: &
       Time                             ! current time
  ! the following arg are inout because they only updated where the data exist
  real, intent(inout) :: &
       swdn_vis_dir(:,:), &
       swdn_vis_dif(:,:), &
       swdn_tot_dir(:,:), &
       swdn_tot_dif(:,:), &
       lwdn(:,:),  & ! downward LW radiation @ sfc, W/m2
       Ps(:,:),    & ! surface pressure, N/m2
       temp(:,:),  & ! temperature, degK
       sphum(:,:), & ! specific humidity, kg/kg
       U (:,:),    & ! zonal wind, m/s
       V (:,:),    & ! meridional wind, m/s
       lprec(:,:), & ! liquid precip, kg/m2/s
       fprec(:,:)    ! frozen precip, kg/m2/s

  ! ---- local vars ---------------------------------------------------------
  real :: &
       sw_ave,   &  ! incident SW at TOA, averaged over input data time interval
       sw_inst,  &  ! incident SW at TOA, averaged over model time step
       swdn,     &  ! total shortwave downward radiation
       cosz,     &  ! cosine of zenith angle
       fracday,  &  ! fraction of day
       rrsun,    &  ! insolation at the Earth orbit
       sw_factor    ! scaling factor for solar radiation
  integer :: t0,t1
  real :: weight

  ! ---- preprocessor definition - shorthand for time interpolation 
#define __INTERP__(id) field(id)%data(t0)*(1-weight)+field(id)%data(t1)*weight 

  ! interpolate in time : t0 and t1 are the indices of the lower and upper boundaries 
  ! of the time axis interval that contains "time"
  if (timeline == 'loop') then
     call time_interp(time, ts, te, time_line, weight, t0, t1,&
          correct_leap_year_inconsistency=.TRUE.)
  else
     call time_interp(time, time_line, weight, t0, t1, YEAR)
  endif

  lprec = __INTERP__(iprec)
  temp  = __INTERP__(itemp)
  sphum = __INTERP__(iq)
  U     = __INTERP__(iwind)
  V     =   0 ! set meridional wind to zero
  ps    = __INTERP__(ips)
  lwdn  = __INTERP__(ilwdn)

  ! distribute precipitation between liquid and solid phase
  fprec = 0
  where (temp < t_fprec)
     fprec = lprec
     lprec = 0
  endwhere

  sw_ave = 0 ; sw_inst = 0 ; cosz = 0 ; fracday = 0

  if(use_diurnal_solar) then
     ! calculate the proper time for the sw fluxes: for diurnal solar 
     ! recalculations we assume that the data are averages over the entire
     ! time interval, so first we calculate the 
     if (timeline == 'loop') then
        call time_interp(time+forcing%dt/2, ts, te, time_line, weight, t0, t1,&
             correct_leap_year_inconsistency=.TRUE.)
     else
        call time_interp(time+forcing%dt/2, time_line, weight, t0, t1, YEAR)
     endif

     swdn = max(0.0,field(iswdn)%data(t0))
     swdn_tot_dir = swdn*swdn_dir_part
     swdn_tot_dif = swdn-swdn_tot_dir
     swdn_vis_dir = swdn_tot_dir*swdn_vis_part
     swdn_vis_dif = swdn_tot_dif*swdn_vis_part

     ! calculate average TOA downward SW for data time interval
     call diurnal_solar( lat, lon, time_line(t0), cosz, fracday, rrsun, time_line(2)-time_line(1))
     sw_ave = cosz*fracday*rrsun*solar
     ! calculate instant TOA downward SW for model time step
     call diurnal_solar( lat, lon, time, cosz, fracday, rrsun, Forcing%dt )
     sw_inst = cosz*fracday*rrsun*solar
     if(sw_ave>0) then
        sw_factor = sw_inst/sw_ave
     else
        sw_factor = 0.0
     endif
     if(limit_swdn_by_toa_flux.and.swdn*sw_factor>sw_inst) &
          sw_factor = sw_inst/swdn
     swdn_vis_dir = swdn_vis_dir*sw_factor
     swdn_vis_dif = swdn_vis_dif*sw_factor
     swdn_tot_dir = swdn_tot_dir*sw_factor
     swdn_tot_dif = swdn_tot_dif*sw_factor
  else
     swdn = max(0.0,__INTERP__(iswdn))
     swdn_tot_dir = swdn*swdn_dir_part
     swdn_tot_dif = swdn-swdn_tot_dir
     swdn_vis_dir = swdn_tot_dir*swdn_vis_part
     swdn_vis_dif = swdn_tot_dif*swdn_vis_part
  endif

end subroutine get_prescr_forcing

! ============================================================================
subroutine cdf_err(err, file, line, level)
! ============================================================================
! prints netcdf error message
  ! ---- arguments -----------------------------------------------------------
  integer, intent(in)           :: err   ! error code
  character (len=*), intent(in) :: file  ! file name
  integer, intent(in)           :: line  ! line number
  integer, intent(in)           :: level ! error level

  ! ---- local vars ----------------------------------------------------------
  character(len=2048) :: info

  if (err /= NF_NOERR) then
     write(info,'("File ::",a," Line ::",i4.4)')trim(file), line
     call error_mesg(info, nf_strerror(err), level)
  endif
end subroutine cdf_err

end module prescr_forcing_mod
