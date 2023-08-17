module river_mod
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
! <CONTACT EMAIL="Kirsten.Findell@@noaa.gov"> Kirsten Findell </CONTACT> 
! <CONTACT EMAIL="Zhi.Liang@@noaa.gov"> Zhi Liang </CONTACT> 

  use utilities_mod,       only : logunit, get_unit
  use time_manager_mod,    only : time_type, increment_time, get_time
  use river_type_mod,      only : river_type
  use constants_mod,       only : PI, RADIAN, tfreeze, DENS_H2O, hlf

  implicit none
  private

!--- version information ---------------------------------------------
  character(len=128) :: version = 'Dummy version of river.F90'
  character(len=128) :: tagname = '  '

!--- public interface ------------------------------------------------
  public :: river_init, river_end, river_type, update_river
  public :: save_river_restart

!---------------------------------------------------------------------

  integer, parameter :: outunit = 6
  real,    parameter :: CONST_OMEAN = 80000

  real,  allocatable, dimension(:,:)   :: discharge2ocean_next   ! store discharge value
  real,  allocatable, dimension(:,:,:) :: discharge2ocean_next_c ! store discharge value
  integer                :: unit
  type(river_type), save :: River

contains

!#####################################################################
  subroutine river_init( land_lon, land_lat, time, dt_fast, &
                         land_frac, id_lon, id_lat, river_land_mask )
    real,            intent(in) :: land_lon(:,:)     ! geographical lontitude of cell center
    real,            intent(in) :: land_lat(:,:)     ! geographical lattitude of cell center
    type(time_type), intent(in) :: time              ! current time
    type(time_type), intent(in) :: dt_fast           ! fast time step
    real,            intent(in) :: land_frac(:,:)    ! land area fraction from land model
    integer,         intent(in) :: id_lon, id_lat    ! IDs of diagnostic axes
    logical,         intent(out):: river_land_mask(:,:) ! land mask seen by rivers

    character(len=128)   :: filename='INPUT/river.res'
    logical              :: fexist

!--- write version and namelist info to logfile --------------------
    write (logunit,'(/,80("="),/(a))') trim(version), trim(tagname)

    return ! Dummy routine

!--- read restart file ---

    inquire (file=trim(filename), exist=fexist)
    if(fexist) then
        unit = get_unit()
        open(unit, file=trim(filename), form='unformatted', action='read')
        read(unit) River%storage
        read(unit) River%storage_c
        read(unit) discharge2ocean_next
        read(unit) discharge2ocean_next_c
        read(unit) River%outflowmean
        close(unit)
        write(outunit,*) 'Read restart files INPUT/river.res'
    else
        River%storage    = 0.0
        River%storage_c  = 0.0
        discharge2ocean_next   = 0.0
        discharge2ocean_next_c = 0.0
        River%outflowmean = CONST_OMEAN
        write(outunit,*) 'cold restart, set data to 0 '
    endif

  end subroutine river_init

!#####################################################################
  subroutine update_river ( runoff, runoff_c, discharge2ocean, discharge2ocean_c )
    real, dimension(:,:),   intent(in)  :: runoff
    real, dimension(:,:,:), intent(in)  :: runoff_c
    real, dimension(:,:),   intent(out) :: discharge2ocean
    real, dimension(:,:,:), intent(out) :: discharge2ocean_c

    integer, save :: n = 0  ! fast time step with each slow time step

    discharge2ocean = 0
    discharge2ocean_c = 0
    return ! Dummy routine

  end subroutine update_river

!#####################################################################

  subroutine river_end

    return ! Dummy routine

  end subroutine river_end

!#####################################################################
  !--- write to restart file
  subroutine save_river_restart(timestamp)
    character(*), intent(in) :: timestamp

    return ! Dummy routine

  end subroutine save_river_restart

!#####################################################################

end module river_mod
