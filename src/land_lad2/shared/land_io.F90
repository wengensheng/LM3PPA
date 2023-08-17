module land_io_mod

use utilities_mod,     only : logunit, error_mesg, FATAL
use constants_mod,     only : PI
!use horiz_interp_mod,  only : horiz_interp_type, &
!     horiz_interp_new, horiz_interp_del, horiz_interp

use nf_utils_mod,      only : nfu_validtype, nfu_get_dim, nfu_get_dim_bounds, &
     nfu_get_valid_range, nfu_is_valid, nfu_inq_var, nfu_get_var

implicit none
private

! ==== public interface ======================================================
public :: print_netcdf_error
! ==== end of public interface ===============================================

!interface read_field
!   module procedure read_field_N_2D, read_field_N_3D
!   module procedure read_field_I_2D, read_field_I_3D
!end interface

! ==== NetCDF declarations ===================================================
include 'netcdf.inc'
#define __NF_ASRT__(x) call print_netcdf_error((x),__FILE__,__LINE__)

! ==== module constants ======================================================
character(len=*), parameter :: &
     module_name = 'land_io_mod', &
     version     = '$Id: land_io.F90,v 1.1.2.5 2012/06/19 18:34:54 pjp Exp $', &
     tagname     = '$Name: no_fms_b_pjp $'


contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

! ============================================================================
!subroutine read_field_N_2D(filename, varname, lon, lat, data, interp)
!  character(len=*), intent(in) :: filename
!  character(len=*), intent(in) :: varname
!  real, intent(in)  :: lon(:,:),lat(:,:)
!  real, intent(out) :: data(:,:)
!  character(len=*), intent(in), optional :: interp

   ! ---- local vars ----------------------------------------------------------
!  real    :: data3(size(data,1),size(data,2),1)

!  call read_field_N_3D(filename, varname, lon, lat, data3, interp)
!  data = data3(:,:,1)

!end subroutine read_field_N_2D

! ============================================================================
!subroutine read_field_N_3D(filename, varname, lon, lat, data, interp)
!  character(len=*), intent(in) :: filename
!  character(len=*), intent(in) :: varname
!  real, intent(in)  :: lon(:,:),lat(:,:)
!  real, intent(out) :: data(:,:,:)
!  character(len=*), intent(in), optional :: interp

   ! ---- local vars ----------------------------------------------------------
!  integer :: ncid
!  integer :: iret

!  iret = nf_open(filename,NF_NOWRITE,ncid)
!  if(iret/=NF_NOERR) then
!     call error_mesg('read_field','Can''t open netcdf file "'//trim(filename)//'"',FATAL)
!  endif
!  call read_field_I_3D(ncid, varname, lon, lat, data, interp)
!  __NF_ASRT__( nf_close(ncid) )

!end subroutine read_field_N_3D

! ============================================================================
!subroutine read_field_I_2D(ncid, varname, lon, lat, data, interp)
!  integer, intent(in) :: ncid
!  character(len=*), intent(in) :: varname
!  real, intent(in) :: lon(:,:),lat(:,:)
!  real, intent(out) :: data(:,:)
!  character(len=*), intent(in), optional  :: interp
  ! ---- local vars
!  real    :: data3(size(data,1),size(data,2),1)

!  call read_field_I_3D(ncid, varname, lon, lat, data3, interp)
!  data = data3(:,:,1)

!end subroutine read_field_I_2D

! ============================================================================
!subroutine read_field_I_3D(ncid, varname, lon, lat, data, interp)
!  integer, intent(in) :: ncid
!  character(len=*), intent(in) :: varname
!  real, intent(in) :: lon(:,:),lat(:,:)
!  real, intent(out) :: data(:,:,:)
!  character(len=*), intent(in), optional  :: interp

  ! ---- local vars ----------------------------------------------------------
!  integer :: nlon, nlat, nlev ! size of input grid
!  integer :: varndims ! number of variable dimension
!  integer :: vardims(NF_MAX_VAR_DIMS) ! IDs of variable dimension
!  integer :: dimlens(NF_MAX_VAR_DIMS) ! sizes of respective dimensions
!  real,    allocatable :: in_lonb(:), in_latb(:), in_lon(:), in_lat(:)
!  real,    allocatable :: x(:,:,:) ! input buffer
!  logical, allocatable :: mask(:,:,:) ! mask of valid values
!  real,    allocatable :: rmask(:,:,:) ! real mask for interpolator
!  character(len=20) :: interpolation 
!  integer :: i,j,k,imap,jmap !
!  type(nfu_validtype) :: v
!  character(len=4) :: string1, string2

! interpolation = "bilinear"
! if(present(interp)) interpolation = interp
  
  ! get the dimensions of our variable
! __NF_ASRT__( nfu_inq_var(ncid,varname,ndims=varndims,dimids=vardims,dimlens=dimlens) )
! if(varndims<2.or.varndims>3) then
!    write(string1,'(i4)') varndims
!    call error_mesg('read_field','variable "'//trim(varname)//'" is'//string1//&
!         'D, but only reading 2D or 3D variables is supported', FATAL)
! endif
! nlon = dimlens(1) ; nlat = dimlens(2)
! nlev = 1; 
! if (varndims==3) nlev=dimlens(3)
! if(nlev/=size(data,3)) then
!    write(string1,'(i4)') nlev
!    write(string2,'(i4)') size(data,3)
!    call error_mesg('read_field','3rd dimension length of the variable "'&
!         //trim(varname)//'" ('//trim(string1)//') is different from the expected size of data ('// &
!         trim(string2)//')', FATAL)
! endif

! allocate (                 &
!      in_lon  (nlon),   in_lat  (nlat),   &
!      in_lonb (nlon+1), in_latb (nlat+1), &
!      x       (nlon, nlat, nlev) ,&
!      mask    (nlon, nlat, nlev) , rmask(nlon, nlat, nlev) )

  ! read boundaries of the grid cells in longitudinal direction
! __NF_ASRT__(nfu_get_dim(ncid, vardims(1), in_lon))
! __NF_ASRT__(nfu_get_dim(ncid, vardims(2), in_lat))
! in_lon = in_lon*PI/180.0; in_lat = in_lat*PI/180.0
! __NF_ASRT__(nfu_get_dim_bounds(ncid, vardims(1), in_lonb))
! __NF_ASRT__(nfu_get_dim_bounds(ncid, vardims(2), in_latb))
! in_lonb = in_lonb*PI/180.0; in_latb = in_latb*PI/180.0
! __NF_ASRT__(nfu_get_valid_range(ncid,varname,v))
  ! read input data
! __NF_ASRT__( nfu_get_var(ncid,varname,x) ) ! assuming real is real*8
! mask = nfu_is_valid(x,v)
! rmask = 1.0
! where(.not.mask) rmask = 0.0

! select case(trim(interpolation))
! case ("bilinear")
!    do k = 1,size(data,3)
!       call horiz_interp(x(:,:,k), in_lonb, in_latb, lon,lat, data(:,:,k), mask_in=rmask(:,:,k), &
!            interp_method='bilinear')
!    enddo
! case ("nearest")
!    do k = 1,size(data,3)
!    do j = 1,size(data,2)
!    do i = 1,size(data,1)
!       call nearest (mask(:,:,k), in_lon, in_lat, lon(i,j), lat(i,j), imap, jmap)
!       data(i,j,k) = x(imap,jmap,k)
!    enddo
!    enddo
!    enddo
! case default
!    call error_mesg(module_name, interpolation//" is not a valid interpolation method",FATAL)
! end select

! deallocate(in_lonb, in_latb, in_lon, in_lat, x, mask, rmask)

!end subroutine read_field_I_3D

! ============================================================================
subroutine print_netcdf_error(ierr, file, line)
  ! prints out NetCDF library error message, including file name and line number
  integer,          intent(in) :: ierr ! error code
  character(len=*), intent(in) :: file ! name of the file
  integer,          intent(in) :: line ! number of line in the file

  ! ---- local vars
  character(len=1024) :: mesg

  if (ierr.ne.NF_NOERR) then
     write(mesg, "('File ',a,' Line ',i4.4,' :: ',a)") &
          trim(file),line,trim(NF_STRERROR(ierr))
     call error_mesg('NetCDF', mesg, FATAL)
  endif
end subroutine print_netcdf_error

end module
