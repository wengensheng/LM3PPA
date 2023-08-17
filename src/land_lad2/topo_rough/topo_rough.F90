module topo_rough_mod
! <CONTACT EMAIL="slm@gfdl.noaa.gov">
!   Sergey Malyshev
! </CONTACT>

  use utilities_mod,      only : logunit, get_unit, check_nml_error
  use time_manager_mod,   only : time_type, get_date

implicit none
private
! ==== public interface ======================================================
public :: topo_rough_init
public :: topo_rough_end
public :: update_topo_rough
! ==== end of public interface ===============================================


! <NAMELIST NAME = "topo_rough_nml">
!   <DATA NAME="use_topo_rough" TYPE="logical" DEFAULT="false">
!     If true, the topographic momentum drag scaling scheme is used
!   </DATA>
!   <DATA NAME="max_topo_rough" TYPE="real" DEFAULT="100" UNITS="m">
!     Maximum of topographic "roughness length" used for momentum drag scaling
!   </DATA>
!   <DATA NAME="topo_rough_factor" TYPE="real" DEFAULT="1.0">
!     Scaling factor to convert topography variance to topographic 
!     "roughness length"
!   </DATA>
!   <DATA NAME="topo_rough_file" TYPE="character(len=256)" DEFAULT="INPUT/mg_drag.data.nc">
!     Name of the file to be used as an input for sub-grid topography variance data. 
!     The file can be either NetCDF (in this case variable name can also be specified), or
!     IEEE.
!   </DATA>
!   <DATA NAME="topo_rough_var" TYPE="character(len=128)" DEFAULT="ghprime">
!     Name of the NetCDF variable to be used as a topography variance field. Ignored if
!     the file specified in topo_rough_file is not NetCDF file.
!   </DATA>
! </NAMELIST>

logical     :: use_topo_rough    = .false.
real        :: max_topo_rough    = 100 ! m
real        :: topo_rough_factor = 1.0
real        :: topo_stdev = 96.591287022675942
character(len=8), parameter :: topo_rough_source = 'computed'
character(len=256):: topo_rough_file   = 'INPUT/mg_drag.data.nc'
character(len=128):: topo_rough_var    = 'ghprime'

namelist/topo_rough_nml/ use_topo_rough, topo_rough_factor, max_topo_rough, &
         topo_rough_file, topo_rough_var, topo_stdev

! ==== module constants ======================================================
character(len=*), parameter :: &
     module_name   = 'she_topo_rough', &
     diag_mod_name = 'topo_rough', &
     version       = '$Id: topo_rough.F90,v 1.1.2.8 2012/06/19 18:34:54 pjp Exp $', &
     tagname       = '$Name: no_fms_b_pjp $'

! ==== module private data ===================================================
!real, allocatable, save ::topo_stdev(:,:)
logical :: module_is_initialized = .FALSE.

! ==== NetCDF declarations ===================================================
include 'netcdf.inc'

contains ! ###################################################################

subroutine topo_rough_init(time, lonb, latb, id_lon,id_lat)
  type(time_type), intent(in) :: time            ! current time
  real           , intent(in) :: latb(:,:),lonb(:,:) ! boundaries of the grid cells
  integer        , intent(in) :: id_lon,id_lat   ! IDs of diagnostic axes
!   <ERROR MSG="... is not a valid value for topo_rough_source" STATUS="FATAL">
!     specified value of namelist parameter topo_rough_source is invalid; 
!     valid values are 'computed' or 'input'.
!   </ERROR>
  ! --- local vars
  integer :: io
  integer :: id
  logical :: used, got_stdev
  integer :: nml_unit, year, month, day, hour, minute, second

  ! write the version and tagname to the logfile
  write (logunit,'(/,80("="),/(a))') trim(version), trim(tagname)

  ! read and write (to logfile) namelist variables
  nml_unit = get_unit()
  open(nml_unit, file='input.nml', form='formatted', action='read', status='old')
  do 
     read (nml_unit, nml=topo_rough_nml, iostat=io, end=10)
     if (check_nml_error (io, 'topo_rough_nml')==0) exit ! from loop
  enddo
10  close (nml_unit)

  write (logunit, nml=topo_rough_nml)

  ! allocate topo_stdev
 !allocate(topo_stdev(size(lonb,1)-1, size(lonb,2)-1))

  if (use_topo_rough) then
     topo_stdev = min(topo_stdev*topo_rough_factor,max_topo_rough)
  else
     topo_stdev = 0.0
  endif

  ! diag output : send topo_stdev to diagnostics
  id = get_unit()
  open(unit=id, file='DIAGNOSTICS/'//trim(diag_mod_name)//'_topo_rough', form='formatted', action='write', position='rewind')
  write(id,'(a)') 'units = m'
  call get_date(time, year, month, day, hour, minute, second)
  write(id,100) year, month, day, hour, minute, second, topo_stdev
100 format(i4,'-',i2.2,'-',i2.2,1x,i2.2,':',i2.2,':',i2.2,2x,1pe24.16)

  module_is_initialized = .TRUE.
end subroutine topo_rough_init

! ============================================================================
subroutine topo_rough_end()
 !deallocate(topo_stdev)
  module_is_initialized = .FALSE.
end subroutine

! ============================================================================
subroutine update_topo_rough(topo_rough)
  real, intent(out) :: topo_rough(:,:,:)

  ! ---- local vars
  integer :: k

  ! just assign standard deviation (scaled and trimmed according to namelist 
  ! parameters) to the output field 
  do k = 1, size(topo_rough,3)
     topo_rough(:,:,k) = topo_stdev
  enddo
end subroutine

end module topo_rough_mod
