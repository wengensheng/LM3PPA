module land_tile_diag_mod

use utilities_mod,      only : logunit, get_unit
use time_manager_mod,   only : time_type, get_date

use land_tile_selectors_mod, only : tile_selectors_init, tile_selectors_end, &
     tile_selector_type, register_tile_selector, selector_suffix, &
     get_n_selectors, get_selector
use land_tile_mod,      only : land_tile_type, diag_buff_type, &
     land_tile_list_type, first_elmt, tail_elmt, next_elmt, get_elmt_indices, &
     land_tile_enum_type, operator(/=), current_tile, &
     tile_is_selected
use land_data_mod,      only : lnd
use tile_diag_buff_mod, only : diag_buff_type, realloc_diag_buff

implicit none
private


! ==== public interface ======================================================
public :: tile_diag_init
public :: tile_diag_end

public :: diag_buff_type

public :: register_tiled_diag_field
public :: register_tiled_static_field
public :: send_tile_data
public :: send_tile_data_r0d_fptr, send_tile_data_r1d_fptr
public :: send_tile_data_i0d_fptr

public :: dump_tile_diag_fields

! codes of tile aggregaton operations
public :: OP_AVERAGE, OP_SUM

interface send_tile_data
   module procedure send_tile_data_0d
   module procedure send_tile_data_1d
end interface
! ==== end of public interface ===============================================


! ==== module constants ======================================================
character(len=*), parameter :: &
     module_name = 'lan_tile_diag_mod', &
     version     = '$Id: land_tile_diag.F90,v 1.1.2.9 2011/03/28 18:36:39 pjp Exp $', &
     tagname     = '$Name: no_fms_b_pjp $'

! operations used for tile data aggregation
integer, parameter :: OP_AVERAGE = 0 ! weighted average of tile values
integer, parameter :: OP_SUM     = 1 ! sum of all tile values


! ==== module data ===========================================================
logical :: module_is_initialized = .false.



contains



! ============================================================================
subroutine tile_diag_init()

  if (module_is_initialized) return

  module_is_initialized = .true.
  write (logunit,'(/,80("="),/(a))') trim(version), trim(tagname)

  ! initialize diag selectors
  call tile_selectors_init()
  call register_tile_selector('')
end subroutine tile_diag_init



! ============================================================================
subroutine tile_diag_end()

  ! destroy selectors
  call tile_selectors_end()

  module_is_initialized = .false.
end subroutine tile_diag_end


! ============================================================================
function register_tiled_diag_field(mod_name, field_name, long_name, units, op) result (id)

  integer :: id

  character(len=*), intent(in) :: mod_name, field_name
  character(len=*), intent(in), optional :: long_name
  character(len=*), intent(in), optional :: units
  integer,          intent(in), optional :: op ! aggregation operation code
  
  id = reg_field(.false., mod_name, field_name, long_name, units, op=op)

end function register_tiled_diag_field

! ============================================================================
function register_tiled_static_field(mod_name, field_name, long_name, units, op) result (id)

  integer :: id

  character(len=*), intent(in) :: mod_name, field_name
  character(len=*), intent(in), optional :: long_name
  character(len=*), intent(in), optional :: units
  integer,          intent(in), optional :: op ! aggregation operation code
  
  ! --- local vars

  id = reg_field(.true., mod_name, field_name, long_name, units, op)

end function register_tiled_static_field

! ============================================================================
! provides unified interface for registering a diagnostic field with full set
! of selectors
function reg_field(static, mod_name, field_name, long_name, units, op) result(id)
  integer :: id

  logical,          intent(in) :: static
  character(len=*), intent(in) :: mod_name, field_name
  character(len=*), intent(in), optional :: long_name
  character(len=*), intent(in), optional :: units
  integer,          intent(in), optional :: op

  id = get_unit()
  open(unit=id, file='DIAGNOSTICS/'//trim(mod_name)//'_'//trim(field_name), form='formatted', action='write', position='rewind')
  if(present(long_name)) then
     write(id,'(a20,a)') 'long_name = ',trim(long_name)
  else
     write(id,'(a20,a)') 'long_name = ',''
  endif
  if(present(units)) then
     write(id,'(a20,a)') 'units = ',trim(units)
  else
     write(id,'(a20,a)') 'units = ',''
  endif
end function reg_field

! ============================================================================
subroutine send_tile_data_0d(id, x, buffer)
  integer, intent(in) :: id
  real   , intent(in) :: x
  type(diag_buff_type), intent(inout) :: buffer ! buffer is unused in this version

  integer :: year, month, day, hour, minute, second
  character(len=*), parameter :: &
       diag_format='(i4,"-",i2.2,"-",i2.2,1x,i2.2,":",i2.2,":",i2.2,2x,", ",99e16.8)'

  call get_date(lnd%time, year, month, day, hour, minute, second)
  write(id,diag_format) year,month,day,hour,minute,second,x
  
end subroutine send_tile_data_0d

! ============================================================================
subroutine send_tile_data_1d(id, x, buffer)
  integer, intent(in) :: id
  real   , intent(in) :: x(:)
  type(diag_buff_type), intent(inout) :: buffer ! buffer is unused in this version

   integer :: year, month, day, hour, minute, second
  character(len=*), parameter :: &
       diag_format='(i4,"-",i2.2,"-",i2.2,1x,i2.2,":",i2.2,":",i2.2,2x,", ",99e16.8)'

  integer :: idx, i

  call get_date(lnd%time, year, month, day, hour, minute, second)
  write(id,diag_format) year,month,day,hour,minute,second,x
end subroutine send_tile_data_1d

! NOTE: 2-d fields can be handled similarly to 1-d with reshape

! ============================================================================
subroutine send_tile_data_r0d_fptr(id, tile_map, fptr)
  integer, intent(in) :: id
  type(land_tile_list_type), intent(inout) :: tile_map(:,:)
  ! subroutine returning the pointer to the tile data
  interface
     subroutine fptr(tile, ptr)
       use land_tile_mod, only : land_tile_type
       type(land_tile_type), pointer :: tile ! input
       real                , pointer :: ptr  ! returned pointer to the data
     end subroutine fptr 
  end interface

  type(land_tile_enum_type)     :: te,ce   ! tail and current tile list elements
  type(land_tile_type), pointer :: tileptr ! pointer to tile   
  real                , pointer :: ptr     ! pointer to the data element within a tile

  if(id <= 0) return
  ce = first_elmt( tile_map )
  te = tail_elmt ( tile_map )
  do while(ce /= te)
     tileptr => current_tile(ce)
     call fptr(tileptr,ptr)
     if(associated(ptr)) call send_tile_data(id,ptr,tileptr%diag)
     ce=next_elmt(ce)
  enddo
end subroutine send_tile_data_r0d_fptr


! ============================================================================
subroutine send_tile_data_r1d_fptr(id, tile_map, fptr)
  integer, intent(in) :: id
  type(land_tile_list_type), intent(inout) :: tile_map(:,:)
  ! subroutine returning the pointer to the tile data
  interface
     subroutine fptr(tile, ptr)
       use land_tile_mod, only : land_tile_type
       type(land_tile_type), pointer :: tile ! input
       real                , pointer :: ptr(:)  ! returned pointer to the data
     end subroutine fptr 
  end interface

  type(land_tile_enum_type)     :: te,ce   ! tail and current tile list elements
  type(land_tile_type), pointer :: tileptr ! pointer to tile   
  real                , pointer :: ptr(:)     ! pointer to the data element within a tile

  if(id <= 0) return
  ce = first_elmt( tile_map )
  te = tail_elmt ( tile_map )
  do while(ce /= te)
     tileptr => current_tile(ce)
     call fptr(tileptr,ptr)
     if(associated(ptr)) call send_tile_data(id,ptr,tileptr%diag)
     ce=next_elmt(ce)
  enddo
end subroutine send_tile_data_r1d_fptr


! ============================================================================
subroutine send_tile_data_i0d_fptr(id, tile_map, fptr)
  integer, intent(in) :: id
  type(land_tile_list_type), intent(inout) :: tile_map(:,:)
  ! subroutine returning the pointer to the tile data
  interface
     subroutine fptr(tile, ptr)
       use land_tile_mod, only : land_tile_type
       type(land_tile_type), pointer :: tile ! input
       integer             , pointer :: ptr  ! returned pointer to the data
     end subroutine fptr 
  end interface

  type(land_tile_enum_type)     :: te,ce   ! tail and current tile list elements
  type(land_tile_type), pointer :: tileptr ! pointer to tile   
  integer             , pointer :: ptr     ! pointer to the data element within a tile

  if(id <= 0) return
  ce = first_elmt( tile_map )
  te = tail_elmt ( tile_map )
  do while(ce /= te)
     tileptr => current_tile(ce)
     call fptr(tileptr,ptr)
     if(associated(ptr)) call send_tile_data(id,real(ptr),tileptr%diag)
     ce=next_elmt(ce)
  enddo
end subroutine send_tile_data_i0d_fptr


! ============================================================================
subroutine dump_tile_diag_fields(tiles, time) ! This is a dummy routine
  type(land_tile_list_type), intent(in) :: tiles(:,:) ! 
  type(time_type)          , intent(in) :: time       ! current time
end subroutine dump_tile_diag_fields

end module land_tile_diag_mod
