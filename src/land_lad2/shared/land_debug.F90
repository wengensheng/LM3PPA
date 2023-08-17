module land_debug_mod

use    utilities_mod, only: logunit, get_unit, error_mesg, FATAL, check_nml_error
use time_manager_mod, only: time_type, get_date

implicit none
private

! ==== public interfaces =====================================================
public :: land_debug_init
public :: land_debug_end

public :: set_current_point
public :: get_current_point
public :: current_i, current_j, current_k, current_face
public :: is_watch_point
public :: get_watch_point

public :: check_temp_range
public :: check_conservation

public :: dpri

interface dpri
   module procedure debug_printout_r0d
   module procedure debug_printout_i0d
   module procedure debug_printout_l0d
   module procedure debug_printout_r1d
   module procedure debug_printout_i1d
   module procedure debug_printout_r2d
end interface dpri

! conservation tolerances for use across the code. This module doesn't use
! them, just serves as a convenient place to share them across all land code
public :: water_cons_tol
public :: carbon_cons_tol

! ==== module constants ======================================================
character(len=*), parameter, private   :: &
    module_name = 'land_debug',&
    version     = '$Id: land_debug.F90,v 1.1.2.8 2012/06/19 18:34:54 pjp Exp $',&
    tagname     = '$Name: no_fms_b_pjp $'

! ==== module variables ======================================================
integer :: current_debug_level = 0
integer :: mosaic_tile = 0
integer :: curr_i, curr_j, curr_k
character(9) :: fixed_format='(a12,99g)'

!---- namelist ---------------------------------------------------------------
integer :: watch_point(4)=(/0,0,0,1/) ! coordinates of the point of interest, 
           ! i,j,tile,mosaic_tile
real    :: temp_lo = 120.0 ! lower limit of "reasonable" temperature range, deg K
real    :: temp_hi = 373.0 ! upper limit of "reasonable" temperature range, deg K
logical :: print_hex_debug = .FALSE. ! if TRUE, hex representation of debug 
           ! values is also printed
integer :: label_len = 12  ! minimum length of text labels for debug output
logical :: trim_labels = .FALSE. ! if TRUE, the length of text labels in debug 
           ! printout is never allowed to exceed label_len, resulting in 
           ! trimming of the labels. Set it to TRUE to match earlier debug 
           ! printout
real    :: water_cons_tol  = 1e-11 ! tolerance of water conservation checks 
real    :: carbon_cons_tol = 1e-13 ! tolerance of carbon conservation checks  
namelist/land_debug_nml/ watch_point, temp_lo, temp_hi, &
   print_hex_debug, label_len, trim_labels, &
   water_cons_tol, carbon_cons_tol


contains

! ============================================================================
subroutine land_debug_init()
  ! ---- local vars
  integer :: io, nml_unit

  write (logunit,'(/,80("="),/(a))') trim(version), trim(tagname)
  nml_unit = get_unit()
  open(nml_unit, file='input.nml', form='formatted', action='read', status='old')
  do 
     read (nml_unit, nml=land_debug_nml, iostat=io, end=10)
     if (check_nml_error (io, 'land_debug_nml')==0) exit ! from loop
  enddo
10  close (nml_unit)
  write(logunit, nml=land_debug_nml)
  ! set number of our mosaic tile 
  mosaic_tile = 1

end subroutine land_debug_init

! ============================================================================
subroutine land_debug_end()
end subroutine

! ============================================================================
subroutine set_current_point(i,j,k)
  integer, intent(in) :: i,j,k

  curr_i = i ; curr_j = j ; curr_k = k

  current_debug_level = 0
  if ( watch_point(1)==i.and. &
       watch_point(2)==j.and. &
       watch_point(3)==k.and. &
       watch_point(4)==mosaic_tile) then
     current_debug_level = 1
  endif
end subroutine set_current_point

! ============================================================================
subroutine get_current_point(i,j,k,face)
  integer, intent(out), optional :: i,j,k,face
  if (present(i)) i = curr_i
  if (present(j)) j = curr_j
  if (present(k)) k = curr_k
  if (present(face)) face = mosaic_tile
end subroutine get_current_point

! ============================================================================
integer function current_i() ; current_i = curr_i ; end function
integer function current_j() ; current_j = curr_j ; end function
integer function current_k() ; current_k = curr_k ; end function
integer function current_face() ; current_face = mosaic_tile ; end function

! ============================================================================
function is_watch_point()
  logical :: is_watch_point
  is_watch_point = (current_debug_level > 0)
end function is_watch_point

! ============================================================================
subroutine get_watch_point(i,j,k,face)
  integer, intent(out), optional :: i,j,k,face
  if (present(i)) i = watch_point(1)
  if (present(j)) j = watch_point(2)
  if (present(k)) k = watch_point(3)
  if (present(face)) face = watch_point(4)
end subroutine get_watch_point

! ============================================================================
! checks if the temperature within reasonable range, and prints a message
! if it isn't
subroutine check_temp_range(temp, tag, varname, time)
  real, intent(in) :: temp ! temperature to check
  character(*), intent(in) :: tag ! tag to print
  character(*), intent(in) :: varname ! name of the variable for printout
  type(time_type), intent(in) :: time ! current time

  ! ---- local vars
  integer :: y,mo,d,h,m,s ! components of date

  if(temp_lo<temp.and.temp<temp_hi) then
     return
  else
     call get_date(time,y,mo,d,h,m,s)
     write(*,'(a," : ",a,g,4(x,a,i4),x,a,i4.4,2("-",i2.2),x,i2.2,2(":",i2.2))')&
          trim(tag), trim(varname)//' out of range: value=', &
         temp,'at i=',curr_i,'j=',curr_j,'tile=',curr_k,'face=',mosaic_tile, &
         'time=',y,mo,d,h,m,s
  endif
end subroutine 

! ============================================================================
! debug printout procedures
subroutine debug_printout_r0d(description,value)
  character(*), intent(in) :: description
  real        , intent(in) :: value
  
  if (trim_labels.or.len_trim(description)<label_len) then
     write(*,fixed_format,advance='NO')trim(description)
  else
     write(*,'(x,a)',advance='NO')trim(description)
  endif
  write(*,'(g)',advance='NO')value
  if(print_hex_debug) write(*,'(z17)',advance='NO')value
end subroutine


subroutine debug_printout_i0d(description,value)
  character(*), intent(in) :: description
  integer     , intent(in) :: value
  
  if (trim_labels.or.len_trim(description)<label_len) then
     write(*,fixed_format,advance='NO')trim(description),value
  else
     write(*,'(x,a,g)',advance='NO')trim(description),value
  endif
end subroutine


subroutine debug_printout_l0d(description,value)
  character(*), intent(in) :: description
  logical     , intent(in) :: value
  
  if (trim_labels.or.len_trim(description)<label_len) then
     write(*,fixed_format,advance='NO')trim(description),value
  else
     write(*,'(x,a,g)',advance='NO')trim(description),value
  endif
end subroutine


subroutine debug_printout_r1d(description,values)
  character(*), intent(in) :: description
  real        , intent(in) :: values(:)
  
  integer :: i

  if (trim_labels.or.len_trim(description)<label_len) then
     write(*,fixed_format,advance='NO')trim(description)
  else
     write(*,'(x,a,99g)',advance='NO')trim(description)
  endif
  do i = 1,size(values)
     write(*,'(g)',advance='NO')values(i)
     if(print_hex_debug) write(*,'(z17)',advance='NO')values(i)
  enddo
end subroutine

subroutine debug_printout_i1d(description,values)
  character(*), intent(in) :: description
  integer     , intent(in) :: values(:)
  
  integer :: i

  if (trim_labels.or.len_trim(description)<label_len) then
     write(*,fixed_format,advance='NO')trim(description),values
  else
     write(*,'(x,a,99g)',advance='NO')trim(description),values
  endif
end subroutine 

subroutine debug_printout_r2d(description,values)
  character(*), intent(in) :: description
  real        , intent(in) :: values(:,:)
  
  if (trim_labels.or.len_trim(description)<label_len) then
     write(*,fixed_format,advance='NO')trim(description),values
  else
     write(*,'(x,a,99g)',advance='NO')trim(description),values
  endif
  ! TODO: print values as a matrix
end subroutine

! ============================================================================
! checks the conservation of a substance and issues a message with specified
! severity if the difference is not within tolerance.
subroutine check_conservation(tag, substance, d1, d2, tolerance, time, severity)
  character(*), intent(in) :: tag ! message tag (subroutine name or some such)
  character(*), intent(in) :: substance ! name of the substance for printout
  real, intent(in) :: d1,d2 ! values to check
  real, intent(in) :: tolerance ! tolerance of the test
  type(time_type), intent(in) :: time ! current time
  integer, intent(in), optional :: severity ! severity of the non-conservation error:
         ! Can be WARNING, FATAL, or negative. Negative means check is not done.

  ! ---- local vars
  integer :: y,mo,d,h,m,s ! components of date
  character(512) :: message
  integer :: severity_
  
  severity_=FATAL
  if (present(severity))severity_=severity

  if (severity_<0)return

  if (abs(d2-d1)<tolerance) then
     if (is_watch_point()) then
     write(*,'(3(x,a,g))')&
          trim(tag)//': conservation of '//trim(substance)//'; before=', d1, 'after=', d2, 'diff=',d2-d1
     endif
  else
     call get_date(time,y,mo,d,h,m,s)
     write(message,'(3(x,a,g),4(x,a,i4),x,a,i4.4,2("-",i2.2),x,i2.2,2(":",i2.2))')&
          'conservation of '//trim(substance)//' is violated; before=', d1, 'after=', d2, 'diff=',d2-d1,&
          'at i=',curr_i,'j=',curr_j,'tile=',curr_k,'face=',mosaic_tile, &
          'time=',y,mo,d,h,m,s
     call error_mesg(tag,message,severity_)
  endif
end subroutine 

end module land_debug_mod
