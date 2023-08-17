!
!  coupler_main couples component models and controls the time integration
!
program coupler_main
!-----------------------------------------------------------------------
!                   GNU General Public License                        
!                                                                      
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

! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>

! <OVERVIEW>
!  A main program that couples component models for atmosphere and land.
! </OVERVIEW>

! <DESCRIPTION>
!  This version couples model components representing atmosphere and land. 
!  Each model component is represented by a 
!  data type giving the instantaneous model state.
!
!  The component models are coupled to allow implicit vertical diffusion of 
!  heat and moisture at the interfaces of the atmosphere and land models.
!  As a result, the atmosphere and land models use the same time step. 
!  The atmospheric model has been separated into down and up calls that 
!  correspond to the down and up sweeps of the standard tridiagonal elimination.
!
!  This program contains the model's main time loop. Each iteration of the 
!  main time loop is one coupled (slow) time step. Within this slow time step 
!  loop is a fast time step loop, using the atmospheric time step, where the
!  tridiagonal vertical diffusion equations are solved.
!
! <PRE>
!      MAIN PROGRAM EXAMPLE
!      --------------------
!
!         DO slow time steps
!
!              DO fast time steps (atmos)
!
!                   call flux_calculation
!
!                   call ATMOS_DOWN
!
!                   call flux_down_from_atmos
!
!                   call LAND_FAST
!
!                   call flux_up_to_atmos
!
!                   call ATMOS_UP
!
!              END DO
!
!         END DO

!  </PRE>

! </DESCRIPTION>
! <INFO>
!   <NOTE>
!     <PRE>
!   1.If no value is set for current_date, start_date, or calendar (or default value 
!     specified) then the value from restart file "INPUT/coupler.res" will be used. 
!     If neither a namelist value or restart file value exist the program will fail. 
!   2.The actual run length will be the sum of months, days, hours, minutes, and 
!     seconds. A run length of zero is not a valid option. 
!   3.The run length must be an intergal multiple of the coupling timestep dt_cpld. 
!     </PRE>
!   </NOTE>

!   <ERROR MSG="no namelist value for current_date " STATUS="FATAL">
!     A namelist value for current_date must be given if no restart file for
!     coupler_main (INPUT/coupler.res) is found. 
!   </ERROR>
!   <ERROR MSG="invalid namelist value for calendar" STATUS="FATAL">
!     The value of calendar must be 'julian', 'noleap', or 'thirty_day'. 
!     See the namelist documentation. 
!   </ERROR>
!   <ERROR MSG="no namelist value for calendar" STATUS="FATAL">
!     If no restart file is present, then a namelist value for calendar 
!     must be specified. 
!   </ERROR>
!   <ERROR MSG="initial time is greater than current time" STATUS="FATAL">
!     If a restart file is present, then the namelist value for either 
!     current_date or start_date was incorrectly set. 
!   </ERROR>
!   <ERROR MSG="run length must be multiple of slow time step " STATUS="FATAL">
!     There must be an even number of slow time steps for the requested run length. 
!   </ERROR>
!   <ERROR MSG="final time does not match expected ending time " STATUS="WARNING">
!     This error should probably not occur because of checks done at initialization time. 
!   </ERROR>

! </INFO>

  use utilities_mod,           only: logunit, get_unit, error_mesg, NOTE, FATAL, WARNING, &
                                     check_nml_error
  use constants_mod,           only: constants_init

  use time_manager_mod,        only: time_type, set_calendar_type, set_time
  use time_manager_mod,        only: set_date, get_date, days_in_month, month_name
  use time_manager_mod,        only: operator(+), operator(-), operator (<)
  use time_manager_mod,        only: operator (>), operator ( /= ), operator ( / )
  use time_manager_mod,        only: operator (*), THIRTY_DAY_MONTHS, JULIAN
  use time_manager_mod,        only: NOLEAP, NO_CALENDAR, INVALID_CALENDAR
  use time_manager_mod,        only: date_to_string, increment_date
  use time_manager_mod,        only: operator(>=), operator(<=), operator(==)

 !use fms_io_mod,              only: fms_io_exit

!
! model interfaces used to couple the component models: atmosphere and land
!

  use atmos_model_mod,         only: atmos_model_init, atmos_model_end
  use atmos_model_mod,         only: update_atmos_model_down
  use atmos_model_mod,         only: update_atmos_model_up
  use atmos_model_mod,         only: atmos_data_type
  use atmos_model_mod,         only: land_ice_atmos_boundary_type
  use atmos_model_mod,         only: atmos_model_restart

  use land_model_mod,          only: land_model_init, land_model_end
  use land_model_mod,          only: land_data_type, atmos_land_boundary_type
  use land_model_mod,          only: update_land_model_fast, update_land_model_slow
  use land_model_mod,          only: land_model_restart

!
! flux_ calls translate information between model grids - see flux_exchange.f90
!

  use flux_exchange_mod,       only: flux_exchange_init
  use flux_exchange_mod,       only: sfc_boundary_layer
  use flux_exchange_mod,       only: flux_down_from_atmos
  use flux_exchange_mod,       only: flux_up_to_atmos

  implicit none

!-----------------------------------------------------------------------

  character(len=128) :: version = '$Id: coupler_main.F90,v 1.1.2.1.2.16 2011/03/31 13:46:18 Peter.Phillipps Exp $'
  character(len=128) :: tagname = '$Name: no_fms_b_pjp $'

!-----------------------------------------------------------------------
!---- model defined-types ----

  type (atmos_data_type) :: Atm
  type  (land_data_type) :: Land

  type(atmos_land_boundary_type)     :: Atmos_land_boundary
  type(land_ice_atmos_boundary_type) :: Land_ice_atmos_boundary

!-----------------------------------------------------------------------
! ----- coupled model time -----

  type (time_type) :: Time, Time_init, Time_end, &
                      Time_step_atmos, Time_step_cpld
  type(time_type) :: Time_atmos
  integer :: num_atmos_calls, na
  integer :: num_cpld_calls, nc

!------ for intermediate restart
  type(time_type)                      :: Time_restart, Time_restart_current, Time_start
  character(len=32)                    :: timestamp

! ----- coupled model initial date -----

  integer :: date_init(6) = (/ 0, 0, 0, 0, 0, 0 /)
  integer :: calendar_type = INVALID_CALENDAR

!-----------------------------------------------------------------------
!------ namelist interface -------

! <NAMELIST NAME="coupler_nml">
!   <DATA NAME="current_date"  TYPE="integer, dimension(6)"  DEFAULT="0">
!     The date that the current integration starts with. 
!   </DATA>
!   <DATA NAME="force_date_from_namelist"  TYPE="logical"  DEFAULT=".false.">
!     Flag that determines whether the namelist variable current_date should 
!     override the date in the restart file INPUT/coupler.res. If the restart 
!     file does not exist then force_date_from_namelist has not effect, the value of current_date 
!     will be used.
!   </DATA>
!   <DATA NAME="calendar"  TYPE="character(maxlen=17)"  DEFAULT="''">
!     The calendar type used by the current integration. Valid values are consistent 
!     with the time_manager module: 'julian', 'noleap', or 'thirty_day'. The value 
!     'no_calendar' can not be used because the time_manager's date  function are used. 
!   </DATA>
!   <DATA NAME="months "  TYPE="integer"  DEFAULT="0">
!     The number of months that the current integration will be run for. 
!   </DATA>
!   <DATA NAME="days "  TYPE="integer"  DEFAULT="0">
!     The number of days that the current integration will be run for. 
!   </DATA>
!   <DATA NAME="hours"  TYPE="integer"  DEFAULT="0">
!     The number of hours that the current integration will be run for. 
!   </DATA>
!   <DATA NAME="minutes "  TYPE="integer"  DEFAULT="0">
!     The number of minutes that the current integration will be run for. 
!   </DATA>
!   <DATA NAME="seconds"  TYPE="integer"  DEFAULT="0">
!     The number of seconds that the current integration will be run for. 
!   </DATA>
!   <DATA NAME="dt_atmos"  TYPE="integer"  DEFAULT="0">
!     Atmospheric model time step in seconds, including the fast coupling with land.
!   </DATA>
!   <DATA NAME="dt_cpld"  TYPE="integer"  DEFAULT="0">
!     Time step in seconds.
!     Must be an integral multiple of dt_atmos. This is the "slow" timestep.
!   </DATA>
!  <DATA NAME="do_atmos, do_land, do_flux" TYPE="logical">
!  If true (default), that particular model component (atmos, etc.) is run.
!  If false, the execution of that component is skipped. This is used when
!  ALL the output fields sent by that component to the coupler have been
!  overridden using the data_override feature. For advanced users only:
!  if you're not sure, you should leave these values at TRUE.
!  </DATA> 
!  <DATA NAME="atmos_npes, land_npes" TYPE="integer">
!  Both must be 1
!  </DATA> 
!  </DATA>
!  <DATA NAME="n_mask" TYPE="integer">
!    number of region to be masked out. Its value should be less than MAX_PES.
!  </DATA>
!  <DATA NAME="mask_list(2,MAXPES)" TYPE="integer, dimension(2,MAX_MASK_REGION)">
!    The position of the region to be masked out. mask_list(1,:) is the x-layout position
!    and mask_list(2,:) is y-layout position.  
!  </DATA>
!  <DATA NAME="layout_mask" TYPE="integer, dimension(2)">
!   Processor domain layout for all the component model. layout_mask need to be set when and only 
!   when n_mask is greater than 0 ( some domain region is masked out ). When this namelist is set,
!   it will overload the layout in each component model. The default value is (0,0).
!   Currently we require all the component model has the same layout and same grid size.
!  </DATA>
!  <DATA NAME="restart_interval" TYPE="integer, dimension(6)"  DEFAULT="0">
!     The time interval that write out intermediate restart file. The format is (yr,mo,day,hr,min,sec).
!     When restart_interval is all zero, no intermediate restart file will be written out.
!   </DATA>
!   <NOTE>
!     <PRE>
!     1.If no value is set for current_date, start_date, or calendar (or default value specified) then the value from restart
!       file "INPUT/coupler.res" will be used. If neither a namelist value or restart file value exist the program will fail. 
!     2.The actual run length will be the sum of months, days, hours, minutes, and seconds. A run length of zero is not a
!       valid option. 
!     3.The run length must be an intergal multiple of the coupling timestep dt_cpld. 
!     </PRE>
!   </NOTE>
! </NAMELIST>

  integer, parameter :: MAXPES=1
  integer, dimension(6) :: restart_interval = (/ 0, 0, 0, 0, 0, 0/)
  integer, dimension(6) :: current_date     = (/ 0, 0, 0, 0, 0, 0 /)
  character(len=17) :: calendar = '                 '
  logical :: force_date_from_namelist = .false.  ! override restart values for date
  integer :: months=0, days=0, hours=0, minutes=0, seconds=0
  integer :: dt_atmos = 0  ! fluxes passed between atmosphere & land
  integer :: dt_cpld  = 0


  integer :: atmos_npes=1, land_npes=1
  logical :: do_atmos =.true., do_land =.true.
  logical :: do_flux =.true.
  integer :: layout_mask(2) = (/0 , 0/)
  integer :: n_mask = 0
  integer :: mask_list(2, MAXPES), n, m 
  integer, parameter :: mp = 2*MAXPES
  data ((mask_list(n,m),n=1, 2),m=1,MAXPES) /mp*0/

  namelist /coupler_nml/ current_date, calendar, force_date_from_namelist, months, days, hours, &
                         minutes, seconds, dt_cpld, dt_atmos, do_atmos, do_land, do_flux, &
                         atmos_npes, land_npes, n_mask, layout_mask, mask_list, restart_interval

  character(len=80) :: text
  character(len=48), parameter                    :: mod_name = 'coupler_main_mod'
 
  integer, allocatable :: ensemble_pelist(:, :) 
  integer, parameter :: outunit = 6
!
!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------
!

character(len=64), parameter    :: sub_name = 'coupler_main'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: warn_header =                                &
     '==>Warning from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

!#######################################################################

  call constants_init

  call coupler_init

!-----------------------------------------------------------------------
!------ slow integration loop ------

  do nc = 1, num_cpld_calls

     ! To print the value of frazil heat flux at the right time the following block
     ! needs to sit here rather than at the end of the coupler loop.

     if( Atm%pe )then

        !-----------------------------------------------------------------------
        !   ------ atmos/fast-land/fast integration loop -------

        do na = 1, num_atmos_calls

           Time_atmos = Time_atmos + Time_step_atmos

           if (do_flux) then
              call sfc_boundary_layer( REAL(dt_atmos), Time_atmos, Atm, Land, Land_ice_atmos_boundary )
           end if

           !      ---- atmosphere down ----

           if (do_atmos) call update_atmos_model_down( Land_ice_atmos_boundary, Atm )

           call flux_down_from_atmos(Time_atmos, Atm, Land, Land_ice_atmos_boundary, Atmos_land_boundary )

           !      ---- land model ----

           if (do_land)  call update_land_model_fast( Atmos_land_boundary, Land )

           !      ---- atmosphere up ----

           call flux_up_to_atmos( Time_atmos, Land, Land_ice_atmos_boundary, Atmos_land_boundary )

           if (do_atmos) call update_atmos_model_up( Land_ice_atmos_boundary, Atm )

           !--------------

        enddo

        !   ------ end of atmospheric time step loop -----
        if (do_land) call update_land_model_slow(Atmos_land_boundary,Land)
        !-----------------------------------------------------------------------

        Time = Time_atmos
     end if                     !Atm%pe block

     !--- write out intermediate restart file when needed.
     if( Time >= Time_restart ) then
        Time_restart_current = Time
        Time_restart = increment_date(Time, restart_interval(1), restart_interval(2), &
             restart_interval(3), restart_interval(4), restart_interval(5), restart_interval(6) )
        timestamp = date_to_string(time_restart_current)

        write(outunit,*) '=> NOTE from program coupler: intermediate restart file is written and ', &
             trim(timestamp),' is appended as prefix to each restart file name'
        if( Atm%pe )then        
           call atmos_model_restart(Atm, timestamp)
           call land_model_restart(timestamp)
        endif
        call coupler_restart(Time, Time_restart_current, timestamp)
     end if

     !--------------
     write( text,'(a,i4)' )'Main loop at coupling timestep=', nc


  enddo

!-----------------------------------------------------------------------

  call coupler_end

!-----------------------------------------------------------------------

contains

!#######################################################################

  subroutine coupler_init

!-----------------------------------------------------------------------
!   initialize all defined exchange grids and all boundary maps
!-----------------------------------------------------------------------
 
!
!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------
!

    character(len=64), parameter    :: sub_name = 'coupler_init'
    character(len=256), parameter   :: error_header =                               &
         '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
    character(len=256), parameter   :: warn_header =                                &
         '==>Warning from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
    character(len=256), parameter   :: note_header =                                &
         '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

    integer :: unit, io, m, i
    integer :: date(6)
    type (time_type) :: Run_length
    character(len=9) :: month
    integer :: pe, npes

    integer :: ens_siz(4), ensemble_size

    integer :: atmos_pe_start=0, atmos_pe_end=0
    integer :: n
    logical :: other_fields_exist, fexist
    logical, allocatable :: maskmap(:,:)
    character(len=256) :: err_msg
    integer :: date_restart(6)
    character(len=64)  :: filename, fieldname
    integer :: id_restart, l
!-----------------------------------------------------------------------

!----- write version to logfile -------
    open(logunit, file='logfile.out', form='formatted', action='write')
    write (logunit,'(/,80("="),/(a))') trim(version), trim(tagname)

!----- read namelist -------

     unit = get_unit()
     open(unit, file='input.nml', form='formatted', action='read', status='old')
     do 
        read (unit, nml=coupler_nml, iostat=io, end=10)
        if (check_nml_error (io, 'coupler_nml')==0) exit ! from loop
     enddo
10   close (unit)
    
!----- read date and calendar type from restart file -----

    inquire (file='INPUT/coupler.res', exist=fexist)
    if( fexist )then
        unit = get_unit()
        open( unit, file='INPUT/coupler.res', form='formatted', action='read')
        read( unit,* ) calendar_type
        read( unit,* ) date_init
        read( unit,* ) date
        close(unit)
    else
        force_date_from_namelist = .true.
    endif

!----- use namelist value (either no restart or override flag on) ---

    if ( force_date_from_namelist ) then

        if ( sum(current_date) <= 0 ) then
            call error_mesg ('program coupler',  &
                 'no namelist value for base_date or current_date', FATAL)
        else
            date      = current_date
        endif

!----- override calendar type with namelist value -----

        select case( trim(calendar) )
        case( 'JULIAN' )
            calendar_type = JULIAN
        case( 'NOLEAP' )
            calendar_type = NOLEAP
        case( 'THIRTY_DAY' )
            calendar_type = THIRTY_DAY_MONTHS
        case( 'NO_CALENDAR' )
            calendar_type = NO_CALENDAR
        case( 'julian' )
            calendar_type = JULIAN
        case( 'noleap' )
            calendar_type = NOLEAP
        case( 'thirty_day' )
            calendar_type = THIRTY_DAY_MONTHS
        case( 'no_calendar' )
            calendar_type = NO_CALENDAR
        end select

    endif

    call set_calendar_type (calendar_type, err_msg)
    if(err_msg /= '') then
      call error_mesg('coupler_init',trim(err_msg), FATAL)
    endif

    npes = 1

    allocate( Atm%pelist(atmos_npes) )
    Atm%pelist(1) = 0
    Atm%pe  = .true.
    Land%pe = Atm%pe
 
    !Write out messages on root PEs
    write( text,'(a,2i6)' )'Atmos PE range: ', Atm%pelist(1), Atm%pelist(atmos_npes)
    call error_mesg('coupler_init', trim(text), NOTE )

!----- write namelist to logfile -----
    write( logunit, nml=coupler_nml )

!----- write current/initial date actually used to logfile file -----

    write( logunit, 16 )date(1),trim(month_name(date(2))),date(3:6)
16  format ('  current date used = ',i4,1x,a,2i3,2(':',i2.2),' gmt') 

!----- check the value of layout and setup the maskmap for domain layout.
    if( n_mask > 0 ) then
       if(do_atmos .OR. do_land) call error_mesg( &
            'program coupler', 'do_atmos and do_land should be false when n_mask > 0', FATAL)

       if( layout_mask(1)*layout_mask(2) - n_mask .NE. npes ) call error_mesg( &
            'program coupler', 'layout(1)*layout(2) - n_mask should equal to npes when n_mask>0', FATAL)
       call error_mesg('program coupler', 'layout_mask and mask_list is set in coupler_nml, ' // &
                      'the value of layout_mask will override the layout specified in each component model', NOTE)

       allocate(maskmap(layout_mask(1), layout_mask(2)) )
       maskmap = .TRUE.
       do n=1, n_mask
          if (mask_list(1,n) .gt. layout_mask(1) ) &
             call error_mesg( 'program coupler', 'mask_list elements outside layout defines.', FATAL )
          if (mask_list(2,n) .gt. layout_mask(2) ) &
             call error_mesg( 'program coupler', 'mask_list elements outside layout defines.', FATAL )
          maskmap(mask_list(1,n),mask_list(2,n)) = .false.
       enddo
       !--- copy maskmap value to each model data type
       allocate(Atm%maskmap(layout_mask(1), layout_mask(2)), Land%maskmap(layout_mask(1), layout_mask(2)) )
       Atm%maskmap = maskmap;  Land%maskmap = maskmap
       deallocate(maskmap)
    else
       if( layout_mask(1)*layout_mask(2) .NE. 0 ) call error_mesg( &
            'program coupler', 'when no region is masked out, layout_mask need not be set', NOTE )
    end if

!-----------------------------------------------------------------------

!----- use current date if no base date ------

    if ( date_init(1) == 0 ) date_init = date

!----- set initial and current time types ------

    Time_init = set_date (date_init(1), date_init(2), date_init(3), &
         date_init(4), date_init(5), date_init(6))

    Time      = set_date (date(1), date(2), date(3),  &
         date(4), date(5), date(6))

    Time_start = Time

!----- compute the ending time -----

    Time_end = Time
    do m=1,months
       Time_end = Time_end + set_time(0,days_in_month(Time_end))
    end do
    Time_end   = Time_end + set_time(hours*3600+minutes*60+seconds, days)
    Run_length = Time_end - Time

!--- get the time that last intermediate restart file was written out.
    inquire (file='INPUT/coupler.intermediate.res', exist=fexist)
    if( fexist )then
       unit = get_unit()
       open( unit, file='INPUT/coupler.intermediate.res', form='formatted', action='read')
       read(unit,*) date_restart
       close(unit)
    else
       date_restart = date
    endif

    Time_restart_current = Time
    if(ALL(restart_interval ==0)) then
       Time_restart = increment_date(Time_end, 1, 0, 0, 0, 0, 0)   ! no intermediate restart
    else
       Time_restart = set_date(date_restart(1), date_restart(2), date_restart(3),  &
                               date_restart(4), date_restart(5), date_restart(6) )
       Time_restart = increment_date(Time_restart, restart_interval(1), restart_interval(2), &
            restart_interval(3), restart_interval(4), restart_interval(5), restart_interval(6) )
       if(Time_restart <= Time) call error_mesg( &
           'program coupler', 'The first intermediate restart time is no larger than the start time', FATAL)
    end if

!-----------------------------------------------------------------------
!----- write time stamps (for start time and end time) ------

    unit = get_unit()
    open( unit, file='time_stamp.out', form='formatted', action='write')

    month = month_name(date(2))
    write (unit,20) date, month(1:3)

    call get_date (Time_end, date(1), date(2), date(3),  &
         date(4), date(5), date(6))
    month = month_name(date(2))
    write (unit,20) date, month(1:3)

    close(unit)

20  format (6i4,2x,a3)

!-----------------------------------------------------------------------
!----- compute the time steps ------

    Time_step_cpld  = set_time (dt_cpld ,0)
    Time_step_atmos = set_time (dt_atmos,0)

!----- determine maximum number of iterations per loop ------

    num_cpld_calls  = Run_length      / Time_step_cpld
    num_atmos_calls = Time_step_cpld  / Time_step_atmos

!-----------------------------------------------------------------------
!------------------- some error checks ---------------------------------

!----- initial time cannot be greater than current time -------

    if ( Time_init > Time ) call error_mesg ('program coupler',  &
         'initial time is greater than current time', FATAL)

!----- make sure run length is a multiple of slow time step ------

    if ( num_cpld_calls * Time_step_cpld  /= Run_length )  &
         call error_mesg ('program coupler',  &
         'run length must be multiple of coupled time step', FATAL)

! ---- make sure cpld time step is a multiple of atmos time step ----

    if ( num_atmos_calls * Time_step_atmos /= Time_step_cpld )  &
         call error_mesg ('program coupler',   &
         'cpld time step is not a multiple of the atmos time step', FATAL)

!
!       Initialize the tracer manager. This needs to be done on all PEs,
!       before the individual models are initialized.
!

!-----------------------------------------------------------------------
!------ initialize component models ------

    if( Atm%pe )then
!---- atmosphere ----
        call atmos_model_init( Atm, Time_init, Time, Time_step_atmos )

!---- land ----------
        call land_model_init( Atmos_land_boundary, Land, Time_init, Time, &
             Time_step_atmos, Time_step_cpld )
    end if

!-----------------------------------------------------------------------
!---- initialize flux exchange module ----
    call flux_exchange_init ( Time, Atm, Land, land_ice_atmos_boundary, dt_atmos=dt_atmos, dt_cpld=dt_cpld)

    Time_atmos = Time

!
!       read in extra fields for the air-sea gas fluxes
!

!-----------------------------------------------------------------------
!---- open and close dummy file in restart dir to check if dir exists --

    unit = get_unit()
    open( unit, file='RESTART/file', form='formatted', action='write')
    close(unit)

!-----------------------------------------------------------------------
  end subroutine coupler_init

!#######################################################################

  subroutine coupler_end

!-----------------------------------------------------------------------

!----- check time versus expected ending time ----

    if (Time /= Time_end) call error_mesg ('program coupler',  &
         'final time does not match expected ending time', WARNING)

!-----------------------------------------------------------------------
!the call to fms_io_exit has been moved here
    if( Atm%pe )then
        call atmos_model_end (Atm)
        call  land_model_end (Atmos_land_boundary, Land)
    end if

    !----- write restart file ------
    call coupler_restart(Time, Time_restart_current)

   !call fms_io_exit

!-----------------------------------------------------------------------

  end subroutine coupler_end

  !--- writing restart file that contains running time and restart file writing time.
  subroutine coupler_restart(Time_run, Time_res, time_stamp)
    type(time_type),   intent(in)           :: Time_run, Time_res
    character(len=*), intent(in),  optional :: time_stamp
    character(len=128)                      :: file_run, file_res
    integer :: yr, mon, day, hr, min, sec, date(6), unit

    ! write restart file
    if(present(time_stamp)) then
       file_run = 'RESTART/'//trim(time_stamp)//'.coupler.res'
       file_res = 'RESTART/'//trim(time_stamp)//'.coupler.intermediate.res'
    else
       file_run = 'RESTART/coupler.res'
       file_res = 'RESTART/coupler.intermediate.res'
    endif

    !----- compute current date ------
    call get_date (Time_run, date(1), date(2), date(3),  &
                   date(4), date(5), date(6))
    unit = get_unit()
    open( unit, file=trim(file_run), form='formatted', action='write')
    write( unit, '(i6,8x,a)' )calendar_type, &
         '(Calendar: no_calendar=0, thirty_day_months=1, julian=2, gregorian=3, noleap=4)'

    write( unit, '(6i6,8x,a)' )date_init, &
         'Model start time:   year, month, day, hour, minute, second'
    write( unit, '(6i6,8x,a)' )date, &
         'Current model time: year, month, day, hour, minute, second'

    if(Time_res > Time_start) then
       call get_date(Time_res ,yr,mon,day,hr,min,sec)
       write( unit, '(6i6,8x,a)' )yr,mon,day,hr,min,sec, &
            'Current intermediate restart time: year, month, day, hour, minute, second'
    end if
    close(unit)

  end subroutine coupler_restart

!--------------------------------------------------------------------------

  !#######################################################################

  end program coupler_main

