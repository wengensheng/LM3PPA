module utilities_mod

  implicit none
  private

  character(len=128) :: version='$Id  $'
  character(len=128) :: tagname='$Name: no_fms_b_pjp $'

  integer, parameter :: logunit=10
  integer, parameter :: num_atmos_tracers=5, num_land_tracers=2, NO_TRACER=-1
  integer, parameter :: NOTE=0, WARNING=1, FATAL=2

  character(len=8), parameter, dimension(num_atmos_tracers) :: &
    atmos_tracer_names=(/'sphum   ','co2     ','liq_wat ','ice_wat ','cld_amt '/)

  character(len=40), parameter, dimension(num_atmos_tracers) ::          &
    atmos_tracer_longnames=(/'specific humidity                       ', &
                             'carbon dioxide                          ', &
                             'cloud liquid specific humidity          ', &
                             'cloud ice water specific humidity       ', &
                             'cloud fraction                          '/)

  character(len=8), parameter, dimension(num_atmos_tracers) :: &
    atmos_tracer_units=(/'kg/kg   ','kg/kg   ','kg/kg   ','kg/kg   ','none    '/)

  character(len=8), parameter, dimension(num_land_tracers) :: land_tracer_names=(/'sphum   ','co2     '/)

  public :: version, tagname, logunit, get_unit, NO_TRACER, error_mesg, mpp_error, NOTE, WARNING, FATAL
  public :: num_atmos_tracers, atmos_tracer_names, num_land_tracers, land_tracer_names, atmos_tracer_longnames, atmos_tracer_units
  public :: string
  public :: check_nml_error

  !--- public interface ---
  interface string
     module procedure string_from_integer
     module procedure string_from_real
  end interface

  integer :: num_nml_error_codes, nml_error_codes(20)
  logical :: do_nml_error_init = .TRUE.

contains


!======================================================================
function get_unit()
  integer :: get_unit
  logical :: is_open

  do get_unit=103,1024
    inquire(get_unit, OPENED=is_open)
    if(.not.is_open) exit
  enddo
  if(is_open) then
    call error_mesg('get_unit','Too many units open',FATAL) ! this is temporary
  endif
  return
end function get_unit


!======================================================================
subroutine error_mesg (routine, message, level)
  character(len=*), intent(in) :: routine, message
  integer, intent(in)          :: level
  character(len=512)           :: text

  call  mpp_error(level,trim(routine)//': '//trim(message))
end subroutine error_mesg


!======================================================================
subroutine mpp_error(level, message)
  character(len=*), intent(in) :: message
  integer, intent(in)          :: level
  character(len=512)           :: text
  integer                      :: error

  select case( level )
  case(NOTE)
     text = 'NOTE'     !just FYI
  case(WARNING)
     text = 'WARNING'  !probable error
  case(FATAL)
     text = 'FATAL'    !fatal error
  case default
     text = 'WARNING: non-existent errortype (must be NOTE|WARNING|FATAL)'
  end select

  text = trim(text)//': '//trim(message)
  write(6,'(a)') trim(text)

  select case( level )
  case(NOTE)
    ! carry on
  case(WARNING)
    ! carry on
  case(FATAL)
    call FLUSH(6)
    stop
  case default
    call FLUSH(6)
    stop
  end select
end subroutine mpp_error


!======================================================================
function string_from_integer(n)
  integer, intent(in) :: n
  character(len=16) :: string_from_integer

  if(n<0) then
     call mpp_error(FATAL, 'fms_io_mod: n should be non-negative integer, contact developer')
  else if( n<10 ) then
     write(string_from_integer,'(i1)') n
  else if( n<100 ) then
     write(string_from_integer,'(i2)') n
  else if( n<1000 ) then
     write(string_from_integer,'(i3)') n
  else if( n<10000 ) then
     write(string_from_integer,'(i4)') n
  else if( n<100000 ) then
     write(string_from_integer,'(i5)') n
  else if( n<1000000 ) then
     write(string_from_integer,'(i6)') n
  else if( n<10000000 ) then
     write(string_from_integer,'(i7)') n
  else if( n<100000000 ) then
     write(string_from_integer,'(i8)') n
  else
     call mpp_error(FATAL, 'fms_io_mod: n is too big, contact developer')
  end if
end function string_from_integer

!======================================================================
function string_from_real(a)
  real, intent(in) :: a
  character(len=32) :: string_from_real

  write(string_from_real,*) a
end function string_from_real


! ==============================================================================
!   private routine for initializing allowable error codes

subroutine nml_error_init
! some compilers return non-zero iostat values while
! reading through files with multiple namelist records
! this routines "attempts" to identify the iostat values associated
! with records not belonging to the requested namelist

  integer  unit, io, ir
  real    ::  a=1.
  integer ::  b=1
  logical ::  c=.true.
  character(len=8) ::  d='testing'
  namelist /b_nml/  a,b,c,d
  
  nml_error_codes(1) = 0

! create dummy namelist file that resembles actual
  unit = get_unit()
  open (unit, file='_read_error.nml', form='FORMATTED', status='REPLACE')
! due to namelist bug this will not always work
  write (unit, 1)
1 format ('    ', &
    /' &a_nml  a=1.  /',    &
    /'#------------------', &
    /' &b_nml  a=5., b=0, c=.false., d=''test'',  &end')
  close (unit)

! read namelist files and save error codes
  open (unit, file='_read_error.nml', form='FORMATTED', status='OLD')
  ir=1; io=1; do
     read  (unit, nml=b_nml, iostat=io, end=2)
     if (io == 0) exit
     ir=ir+1; nml_error_codes(ir)=io
  enddo
2 close (unit, STATUS='DELETE')

  num_nml_error_codes = ir
  write(*,*) 'nml_error_codes=',nml_error_codes(1:ir)
  do_nml_error_init = .false.
end subroutine nml_error_init


! ==============================================================================
function check_nml_error (iostat, nml_name) result (error_code)
  integer,          intent(in) :: iostat
  character(len=*), intent(in) :: nml_name

  integer   error_code, i
  character(len=128) :: err_str

  if (do_nml_error_init) call nml_error_init()

  error_code = iostat

  do i = 1, num_nml_error_codes
    if (error_code == nml_error_codes(i)) return
  enddo

  call mpp_error ( FATAL, &
    'error reading namelist "'//trim(nml_name)//'", iostat='//string(error_code))
end function check_nml_error

end module utilities_mod
