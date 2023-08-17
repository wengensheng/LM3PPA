! ============================================================================
! updates carbon pools and rates on the fast time scale
! ============================================================================
module vegn_dynamics_mod

#include "../shared/debug.inc"

use utilities_mod,    only: logunit,error_mesg, WARNING

use time_manager_mod, only: time_type, get_date
use constants_mod, only : PI,tfreeze
use land_constants_mod, only : seconds_per_year, mol_C
use land_tile_diag_mod, only : &
     register_tiled_diag_field, send_tile_data, diag_buff_type
use land_data_mod,      only : land_state_type, lnd
use vegn_data_mod, only : spdata, NSPECIES,                                      &
     CMPT_NSC, CMPT_SAPWOOD, CMPT_LEAF, CMPT_ROOT, CMPT_VLEAF, CMPT_WOOD,        &
     LEAF_ON, LEAF_OFF, &
     fsc_liv, fsc_wood, K1, K2, K_nitrogen, MLmixRatio,   &
     soil_carbon_depth_scale, C2B, agf_bs, &
     l_fract, mcv_min, mcv_lai, do_ppa, b0_growth, tau_seed, &
     understory_lai_factor, bseed_distr,wood_fract_min,      &
     CN0fastSOM,CN0slowSOM ! Nitrogen

use vegn_tile_mod, only: vegn_tile_type, vegn_tile_carbon
use soil_tile_mod, only: soil_tile_type, soil_ave_temp, soil_ave_theta, clw, csw
use vegn_cohort_mod, only : vegn_cohort_type, &
     update_biomass_pools, update_bio_living_fraction, update_species, &
     leaf_area_from_biomass, init_cohort_allometry_ppa

use land_debug_mod, only : is_watch_point
use land_numerics_mod, only : rank_descending

implicit none
private

! ==== public interfaces =====================================================
public :: vegn_dynamics_init

public :: vegn_carbon_int_lm3   ! fast time-scale integrator of carbon balance
public :: vegn_carbon_int_ppa   ! fast time-scale integrator of carbon balance
public :: vegn_growth           ! slow time-scale redistributor of accumulated carbon
public :: vegn_phenology_lm3 
public :: vegn_phenology_ppa
public :: vegn_biogeography

public :: vegn_starvation_ppa   !
public :: vegn_reproduction_ppa ! reproduction for PPA case
public :: relayer_cohorts
public :: vegn_mergecohorts_ppa ! merge cohorts
public :: merge_noncanopy_cohorts
public :: kill_lowdensity_cohorts
public :: vegn_annualLAImax
! ==== end of public interfaces ==============================================

! ==== module constants ======================================================
character(len=*), private, parameter :: &
   version = '$Id: vegn_dynamics.F90,v 1.2.2.6 2012/06/19 18:34:55 pjp Exp $', &
   tagname = '$Name: no_fms_b_pjp $' ,&
   module_name = 'vegn'
real, parameter :: GROWTH_RESP=0.333  ! fraction of npp lost as growth respiration


! ==== module data ===========================================================
real    :: dt_fast_yr ! fast (physical) time step, yr (year is defined as 365 days)
integer :: i
! diagnostic field IDs
integer :: id_npp, id_nep, id_gpp, id_rsoil, id_rsoil_fast
integer :: id_N_loss,id_QNmineral
! integer :: id_fsc, id_ssc
integer :: id_resp, id_resl, id_resr, id_resg, id_asoil
integer :: id_soilt, id_theta, id_litter
integer :: id_gpp_c1,id_npp_c1,id_resp_c1,id_nsc_c1 ! Weng 2013-01-21
integer :: id_gpp_under,id_npp_under, id_resp_under
integer :: id_cc_An_op,id_cc_An_cl,id_cc_gpp,id_cc_npp,id_cc_resp,id_cc_ca,id_cc_la
integer :: id_cc_PFT, id_CC_Nup

contains

! ============================================================================
subroutine vegn_dynamics_init(id_lon, id_lat, time, delta_time)
  integer        , intent(in) :: id_lon ! ID of land longitude (X) axis 
  integer        , intent(in) :: id_lat ! ID of land latitude (Y) axis
  type(time_type), intent(in) :: time       ! initial time for diagnostic fields
  real           , intent(in) :: delta_time ! fast time step, s

  write (logunit,'(/,80("="),/(a))') trim(version), trim(tagname)

  ! set up global variables
  dt_fast_yr = delta_time/seconds_per_year

  ! register diagnostic fields
  id_N_loss    = register_tiled_diag_field( trim(module_name), 'N_loss',  'Nitrogen loss rate', 'kg N/(m2 year)')
  id_QNmineral = register_tiled_diag_field ( trim(module_name),'QNmineral','Soil mineral nitrogen', 'kg N/m2' )
  id_gpp = register_tiled_diag_field ( trim(module_name), 'gpp',  'gross primary productivity', 'kg C/(m2 year)')
  id_npp = register_tiled_diag_field ( trim(module_name), 'npp',  'net primary productivity', 'kg C/(m2 year)')
  id_nep = register_tiled_diag_field ( trim(module_name), 'nep',  'net ecosystem productivity', 'kg C/(m2 year)')
  id_litter = register_tiled_diag_field (trim(module_name),  'litter',   'litter productivity', 'kg C/(m2 year)')
 ! id_fsc = register_tiled_diag_field ( trim(module_name), 'fsc',   'fast soil carbon', 'kg C/m2')
 ! id_ssc = register_tiled_diag_field ( trim(module_name), 'ssc',   'slow soil carbon', 'kg C/m2')
  id_resp = register_tiled_diag_field ( trim(module_name), 'resp', 'respiration', 'kg C/(m2 year)')
  id_resl = register_tiled_diag_field ( trim(module_name), 'resl', 'leaf respiration', 'kg C/(m2 year)')
  id_resr = register_tiled_diag_field ( trim(module_name), 'resr', 'root respiration', 'kg C/(m2 year)')
  id_resg = register_tiled_diag_field ( trim(module_name), 'resg', 'growth respiration', 'kg C/(m2 year)')
  id_rsoil = register_tiled_diag_field ( trim(module_name), 'rsoil',  'soil respiration', 'kg C/(m2 year)')
  id_rsoil_fast = register_tiled_diag_field( trim(module_name), 'rsoil_fast',   'fast soil carbon respiration', 'kg C/(m2 year)')
  id_asoil = register_tiled_diag_field( trim(module_name), 'asoil',   'aerobic activity modifier')
  id_soilt = register_tiled_diag_field( trim(module_name), 'tsoil_av',   'average soil temperature for carbon decomposition', 'degK')
  id_theta = register_tiled_diag_field( trim(module_name), 'theta',   'average soil wetness for carbon decomposition', 'm3/m3')

! Weng, for output the first cohort's variables
  id_gpp_c1  = register_tiled_diag_field ( trim(module_name), 'gpp_c1',  'gross primary productivity of cohort1', 'kg C/(tree year)')
  id_npp_c1  = register_tiled_diag_field ( trim(module_name), 'npp_c1',  'net primary productivity of cohort1',   'kg C/(tree year)')
  id_resp_c1 = register_tiled_diag_field ( trim(module_name), 'resp_c1', 'respiration of cohort1', 'kg C/(tree year)')
  id_nsc_c1  = register_tiled_diag_field ( trim(module_name), 'nsc_c1',  'non-structural carbon of cohort1', 'kg C/tree')

  id_gpp_under  = register_tiled_diag_field ( trim(module_name), 'gpp_under',  'understory gross primary productivity', 'kg C/(m2 year)')
  id_npp_under  = register_tiled_diag_field ( trim(module_name), 'npp_under',  'understory net primary productivity',   'kg C/(m2 year)')
  id_resp_under = register_tiled_diag_field ( trim(module_name), 'resp_under', 'understory respiration', 'kg C/(m2 year)')
! Weng, for all cohorts
  id_cc_PFT    = register_tiled_diag_field ( trim(module_name), 'cc_PFT','species')
  id_cc_An_op  = register_tiled_diag_field ( trim(module_name), 'cc_An_op','An_op,leaf photosynthesis', 'molC/(m2L yr)')
  id_cc_An_cl  = register_tiled_diag_field ( trim(module_name), 'cc_An_cl','An_cl,leaf respiration', 'molC/(m2L yr)')
  id_cc_gpp  = register_tiled_diag_field ( trim(module_name),   'cc_gpp',  'gross primary productivity', 'kg C/(tree year)')
  id_cc_npp = register_tiled_diag_field ( trim(module_name),    'cc_npp',  'net primary productivity', 'kg C/(tree year)')
  id_cc_resp = register_tiled_diag_field ( trim(module_name),   'cc_resp','respiration', 'kg C/(tree year)')
  id_cc_ca = register_tiled_diag_field ( trim(module_name),   'cc_crownarea','crown area', 'm2 tree-1')
  id_cc_la = register_tiled_diag_field ( trim(module_name),   'cc_leafarea','leaf area', 'm2 tree-1')
  id_cc_Nup = register_tiled_diag_field ( trim(module_name),   'cc_N_uptake','N_uptake', 'kg N tree-1')
end subroutine vegn_dynamics_init

! =============================================================================
! Added by Weng 2012-02-29
subroutine vegn_annualLAImax(vegn)
  type(vegn_tile_type), intent(inout) :: vegn

  ! ---- local vars
  real   :: currLAImax,nextLAImax
  real   :: LAI_Nitrogen,ccLAImax
  logical:: mineral_N_based

  mineral_N_based = .TRUE.

 if(mineral_N_based)then
    do i=0,NSPECIES
        associate (sp => spdata(i) )
!       LAI_nitrogen = 0.5*vegn%annualN  * sp%CNleaf0 * sp%leafLS / sp%LMA
        LAI_nitrogen = 0.5*vegn%previousN * sp%CNleaf0 * sp%leafLS / sp%LMA
        sp%LAImax = MAX(0.05, MIN(LAI_nitrogen,sp%LAI_light))
        sp%underLAImax = MIN(sp%LAImax,1.2)

        end associate
    enddo
 else ! cc%NSN based
!   Reset sp%LAImax
    do i = 1,vegn%n_cohorts
       associate ( cc => vegn%cohorts(i), &
                   sp => spdata(cc%species) )
       if(cc%layer == 1 .and. cc%crownarea>3. .and. cc%NSN>0.)then
           sp%LAImax    = 0.0  ! max(sp%LAImax,ccLAImax)
           sp%layerfrac = 0.0
           sp%n_cc      = 0
       endif
       end associate
    enddo
!   Sum ccLAImax in the first layer
    do i = 1,vegn%n_cohorts
       associate ( cc => vegn%cohorts(i), &
                   sp => spdata(cc%species) )
       if(cc%layer == 1 .and. cc%crownarea>3. .and. cc%NSN>0.)then
           LAI_nitrogen = 0.15*cc%NSN* sp%CNleaf0 * sp%leafLS / sp%LMA / cc%crownarea
           ccLAImax     = MIN(LAI_nitrogen,sp%LAI_light)
           sp%LAImax    = sp%LAImax + ccLAImax*cc%layerfrac  ! max(sp%LAImax,ccLAImax)
           sp%n_cc      = sp%n_cc + 1
           sp%layerfrac = sp%layerfrac + cc%layerfrac
       endif
       end associate
    enddo
!   Mean LAImax
    do i=0,NSPECIES
       if(spdata(i)%n_cc>=1) &
!          spdata(i)%LAImax = MAX(1.0, spdata(i)%LAImax/spdata(i)%n_cc)
           spdata(i)%LAImax = MAX(1.0, spdata(i)%LAImax/spdata(i)%layerfrac)
           spdata(i)%underLAImax = MIN(spdata(i)%LAImax,1.2)
    enddo
 endif   
end subroutine vegn_annualLAImax
! =============================================================================
! Added by Weng 2012-02-29
subroutine vegn_phenology_ppa(vegn)
  type(vegn_tile_type), intent(inout) :: vegn

  ! ---- local vars
  integer :: i
  real    :: leaf_litter,litterN
  real    :: stem_fall, stem_litter, grassdensity   ! for grasses only
  real    :: leaf_add, root_add       ! per day
  real    :: leaf_fall, leaf_fall_rate ! per day
  real    :: root_mortality, root_mort_rate
  real    :: BL_u,BL_c
  real    :: retransN  ! retranslocation coefficient of Nitrogen
  logical :: cc_firstday = .false.
  logical :: growingseason

  retransN = 0.5
  leaf_fall_rate = 0.075
  root_mort_rate = 0.0
  vegn%litter = 0   ! daily litter
! ON and OFF of phenology: change the indicator of growing season for deciduous
  do i = 1,vegn%n_cohorts   
     associate ( cc => vegn%cohorts(i), &
                 sp => spdata(cc%species) )

!    EVERGREEN: status = LEAF_ON
     if(sp%phenotype==1)then
         if(cc%status==LEAF_OFF)cc%status=LEAF_ON
         BL_u = sp%LMA*cc%crownarea*(1.0-sp%internal_gap_frac)* &
                sp%underLAImax
!               MAX(understory_lai_factor*sp%LAImax,1.0)

         BL_c = sp%LMA * sp%LAImax * cc%crownarea * &
                (1.0-sp%internal_gap_frac)
         if (cc%layer > 1 .and. cc%firstlayer == 0) then ! changed back, Weng 2014-01-23
            cc%topyear = 0.0
            cc%bl_max = BL_u
            cc%br_max = 0.8*cc%bl_max/(sp%LMA*sp%SRA)
         else
!           update the years this cohort is in the canopy layer, Weng 2014-01-06
            if(cc%layer == 1)cc%topyear = cc%topyear + 1./365.  ! daily step
            if(cc%layer>1)cc%firstlayer = 0 ! Just for the first year, those who were
                                            ! pushed to understory have the
                                            ! characteristics of canopy trees
            cc%bl_max = BL_u + min(cc%topyear/5.0,1.0)*(BL_c - BL_u)
            cc%br_max = sp%phiRL*cc%bl_max/(sp%LMA*sp%SRA)

         endif
         ! update NSNmax
!         cc%NSNmax = (cc%bl_max +  cc%br_max)/8.0
         cc%NSNmax = 5.0*(cc%bl_max/(sp%CNleaf0*sp%leafLS)+cc%br_max/sp%CNroot0)
     endif ! sp%phenotype==1
!     write(*,'(5(f9.4))')cc%nindivs, cc%nsc, cc%bsw, cc%bwood, cc%bl
     ! if drought-deciduous or cold-deciduous species
     ! temp=10 degrees C gives good growing season pattern        
     ! spp=0 is c3 grass,1 c3 grass,2 deciduous, 3 evergreen
     ! assumption is that no difference between drought and cold deciduous
!     cc%gdd = vegn%gdd  ! turned off by Weng 2014-07-15
     cc_firstday = .false.
     if((sp%lifeform /=0 .OR. (sp%lifeform ==0 .and. cc%layer==1)).and. & 
         sp%phenotype==0 .and. cc%status==LEAF_OFF .and.                &
         cc%gdd>sp%gdd_crit .and. vegn%tc_pheno>sp%tc_crit_on)then

         cc%status = LEAF_ON
         cc_firstday = .True.

     elseif(sp%phenotype==0.and.cc%status==LEAF_ON.and.vegn%tc_pheno<sp%tc_crit)then
        cc%status = LEAF_OFF
        cc%gdd   = 0.0
        vegn%gdd = 0.0
     endif

!    for grasses' first day of a growing season
     if(cc_firstday .and. sp%lifeform ==0 .and. &
        cc%age>3 .and.(cc%nsc>2*sp%seedlingsize .or. cc%nsc<1*sp%seedlingsize) ) then
!        reset grass density and size for perenials
         grassdensity = cc%nindivs
         cc%nindivs = grassdensity*cc%nsc /(2*sp%seedlingsize)
         cc%br  = grassdensity/cc%nindivs * cc%br
         cc%nsc = grassdensity/cc%nindivs * cc%nsc
         cc%bsw = bseed_distr(CMPT_SAPWOOD)*sp%seedlingsize  ! for setting up a initial size
         cc%nsc = cc%nsc - cc%bsw

!!        Nitrogen pools
         cc%QNroot = grassdensity/cc%nindivs * cc%QNroot
         cc%NSN    = grassdensity/cc%nindivs * cc%NSN
         cc%QNsw   = cc%bsw/sp%CNsw0
         cc%NSN    = cc%NSN - cc%QNsw

         call init_cohort_allometry_ppa(cc);
         cc%bliving     = cc%br + cc%bl + cc%bsw + cc%blv
         cc%DBH_ys      = cc%DBH
         cc%BM_ys       = cc%bsw + cc%bwood
         cc%leaf_age = 0.0

         if (cc%nindivs<=0.0) call error_mesg('vegn_phenology_ppa','initial density too low',WARNING)
     endif

!    calculate targe ammounts of leaves and fine roots of this growing season
!    for deciduous species
     if(cc_firstday) then
         BL_u = sp%LMA*cc%crownarea*(1.0-sp%internal_gap_frac)* &
                sp%underLAImax
            !   MAX(understory_lai_factor*sp%LAImax,1.0)

         BL_c = sp%LMA * sp%LAImax * cc%crownarea * &
              (1.0-sp%internal_gap_frac)
         if (cc%layer > 1 .and. cc%firstlayer == 0) then ! changed back, Weng 2014-01-23
!           update the years this cohort is in the canopy layer, Weng 2014-01-06
            cc%topyear = 0.0 
!           The factor (1-spdata(sp)%internal_gap_frac) is added by Weng 2013-12-19
!           to counteract LAI calculation
!           So, the design of internal_gap actually reduce effective crown area
            cc%bl_max = BL_u
!            cc%br_max = sp%phiRL * sp%LAImax/sp%SRA * cc%crownarea * understory_lai_factor ! Weng 2013-09-04
!           Keep understory tree's root low and constant
            cc%br_max = 0.8*cc%bl_max/(sp%LMA*sp%SRA)
         else
!           update the years this cohort is in the canopy layer, Weng 2014-01-06
            if (cc%layer == 1) cc%topyear = cc%topyear + 1.0 
            cc%bl_max = BL_u + min(cc%topyear/5.0,1.0)*(BL_c - BL_u)
            cc%br_max = sp%phiRL*cc%bl_max/(sp%LMA*sp%SRA)
            if(cc%layer>1)cc%firstlayer = 0 ! Just for the first year, those who were
                                            ! pushed to understory have the characteristics of canopy trees
         endif
         ! update NSNmax
!         cc%NSNmax = (cc%bl_max +  cc%br_max)/8.0
         cc%NSNmax = 5.0*(cc%bl_max/(sp%CNleaf0*sp%leafLS)+cc%br_max/sp%CNroot0)
     endif  ! first day


!    Ending a growing season: leaves fall for deciduous
     if(cc%status == LEAF_OFF .AND. cc%bl > 0.)then
        leaf_fall = min(leaf_fall_rate * cc%bl_max, cc%bl)
        if(sp%lifeform==0)then  ! grasses
            stem_fall = MIN(1.0,leaf_fall/cc%bl) * cc%bsw 
        else                    ! trees
            stem_fall = 0.0
        endif

        root_mortality = min( root_mort_rate * cc%br_max, cc%br)  ! Just for test: keeps roots
!        root_mortality = min( leaf_fall_rate * cc%br_max, cc%br)  ! Just for test

        cc%nsc     = cc%nsc + l_fract * (leaf_fall+ root_mortality+stem_fall)
        cc%bl      = cc%bl  - leaf_fall
        cc%br      = cc%br  - root_mortality
        cc%bsw     = cc%bsw - stem_fall
        ! cc%bwood   = cc%bwood*(1.0-min(1.0,stem_fall/(cc%bsw+cc%bwood)))

!       update NPP for leaves, fine roots, and wood
        cc%NPPleaf = cc%NPPleaf - l_fract * leaf_fall
        cc%NPProot = cc%NPProot - l_fract * root_mortality
        cc%NPPwood = cc%NPPwood - l_fract * stem_fall
        cc%lai     = leaf_area_from_biomass(cc%bl,cc%species,cc%layer,cc%firstlayer)/ &
                                    (cc%crownarea *(1.0-sp%internal_gap_frac))
        if(cc%bl == 0.)cc%leaf_age = 0.0
        cc%bliving = cc%blv + cc%br + cc%bl + cc%bsw

        leaf_litter = (1.-l_fract) * (leaf_fall+root_mortality+stem_fall) * cc%nindivs
        vegn%litter = vegn%litter + leaf_litter
        vegn%fast_soil_C = vegn%fast_soil_C +        fsc_liv *leaf_litter
        vegn%slow_soil_C = vegn%slow_soil_C + (1.0 - fsc_liv)*leaf_litter
        vegn%fsc_in  = vegn%fsc_in  + leaf_litter
        vegn%veg_out = vegn%veg_out + leaf_litter

!!       Nitrogen retransloaction/resorption
        if(cc%QNleaf+cc%QNroot>0.0)cc%NSN = cc%NSN +  &
                  retransN * (cc%QNleaf+cc%QNroot)
!                 retransN * (leaf_fall    / cc%bl * cc%QNleaf &
!                          + root_mortality/ cc%br * cc%QNroot  &
!                          + stem_fall     / cc%bsw* cc%QNsw)
!        cc%QNleaf = 0.0 ! cc%QNleaf - leaf_fall    / cc%bl * cc%QNleaf
!        cc%QNroot = 0.0 ! cc%QNroot - root_mortality/ cc%br * cc%QNroot
!        cc%QNsw   = cc%QNsw   - stem_fall     / cc%bsw* cc%QNsw

        litterN = (1.-retransN) * cc%nindivs * (cc%QNleaf+cc%QNroot)
!                         (  leaf_fall     / cc%bl * cc%QNleaf &
!                          + root_mortality/ cc%br * cc%QNroot  &
!                          + stem_fall     / cc%bsw* cc%QNsw) 
!        vegn%litterN   = vegn%litterN + litterN
        vegn%QNfastSOM = vegn%QNfastSOM +        fsc_liv  * litterN
        vegn%QNslowSOM = vegn%QNslowSOM + (1.0 - fsc_liv) * litterN

        cc%QNleaf = 0.0 ! cc%QNleaf - leaf_fall    / cc%bl * cc%QNleaf
        cc%QNroot = 0.0 ! cc%QNroot - root_mortality/ cc%br * cc%QNroot

     endif


!    Grass leaf and stem fall when a cohort is not in the canopy layer
     if(sp%lifeform ==0 .AND. cc%layer >1 .AND. cc%status == LEAF_ON .AND. cc%bl > 0.)then
        leaf_fall = min(leaf_fall_rate * cc%bl_max, cc%bl)
        stem_fall = MIN(1.0,leaf_fall/cc%bl) * cc%bsw 

!       root_mortality = min( root_mort_rate * cc%br_max, cc%br)  ! Just for test
        root_mortality = min( leaf_fall_rate * cc%br_max, cc%br)  ! Just for test

        cc%nsc     = cc%nsc + l_fract * (leaf_fall+ root_mortality+stem_fall)
        cc%bl      = cc%bl  - leaf_fall
        cc%br      = cc%br  - root_mortality
        cc%bsw     = cc%bsw - stem_fall
        ! cc%bwood   = cc%bwood*(1.0-min(1.0,stem_fall/(cc%bsw+cc%bwood)))
        cc%NPPleaf = cc%NPPleaf - l_fract * leaf_fall
        cc%NPProot = cc%NPProot - l_fract * root_mortality
        cc%NPPwood = cc%NPPwood - l_fract * stem_fall
        cc%lai     = leaf_area_from_biomass(cc%bl,cc%species,cc%layer,cc%firstlayer)/ &
                                    (cc%crownarea *(1.0-sp%internal_gap_frac))
        cc%bliving = cc%blv + cc%br + cc%bl + cc%bsw

        leaf_litter = (1.-l_fract) * (leaf_fall+root_mortality+stem_fall) * cc%nindivs
        vegn%litter = vegn%litter + leaf_litter
        vegn%fast_soil_C = vegn%fast_soil_C +        fsc_liv *leaf_litter
        vegn%slow_soil_C = vegn%slow_soil_C + (1.0 - fsc_liv)*leaf_litter
        vegn%fsc_in  = vegn%fsc_in  + leaf_litter
        vegn%veg_out = vegn%veg_out + leaf_litter

        ! reset cc%status and DBH, CA, Hight.
        if(cc%bsw < 0.00000001) then 
            cc%status = LEAF_OFF
            cc%leaf_age = 0.0
            cc%bsw      = bseed_distr(CMPT_SAPWOOD)* cc%nsc
            call init_cohort_allometry_ppa(cc)
        endif

!!       Nitrogen cycle
        cc%NSN = cc%NSN + retransN * (cc%QNleaf+cc%QNroot+cc%QNsw)
!                 retransN * (leaf_fall    / cc%bl * cc%QNleaf &
!                          + root_mortality/ cc%br * cc%QNroot  &
 !                         + stem_fall     / cc%bsw* cc%QNsw)

        cc%QNleaf = 0.0 !cc%QNleaf - leaf_fall    / cc%bl * cc%QNleaf
        cc%QNroot = 0.0 !cc%QNroot - root_mortality/ cc%br * cc%QNroot
        cc%QNsw   = 0.0 !cc%QNsw   - stem_fall     / cc%bsw* cc%QNsw

        litterN = (1.-retransN) * cc%nindivs *  (cc%QNleaf+cc%QNroot+cc%QNsw)
!                         (  leaf_fall     / cc%bl * cc%QNleaf &
!                          + root_mortality/ cc%br * cc%QNroot  &
!                          + stem_fall     / cc%bsw* cc%QNsw) 
!        vegn%litterN   = vegn%litterN + litterN
        vegn%QNfastSOM = vegn%QNfastSOM +        fsc_liv  * litterN
        vegn%QNslowSOM = vegn%QNslowSOM + (1.0 - fsc_liv) * litterN

     endif

     end associate
  enddo

! TODO: fix the mistake below
!  if(vegn%cohorts(1)%gdd == 0.0)vegn%gdd = 0.0 ! incorrect, must not be inside cohort loop
!  vegn%gdd = vegn%cohorts(1)%gdd  ! turned off by Weng 2014-07-15
end subroutine vegn_phenology_ppa

! ============================================================================
! updates cohort biomass pools, LAI, SAI, and height using accumulated 
! carbon_gain and bwood_gain
subroutine vegn_growth (vegn)
  type(vegn_tile_type), intent(inout) :: vegn

  ! ---- local vars
  type(vegn_cohort_type), pointer :: cc    ! current cohort
  real :: cmass0,cmass1 ! for debug only
  integer :: i

  if (is_watch_point()) then
     cmass0 = vegn_tile_carbon(vegn)
  endif

  do i = 1, vegn%n_cohorts   
     cc => vegn%cohorts(i)

     if (do_ppa) then
        call biomass_allocation_ppa(cc)
     else
        cc%bwood   = cc%bwood   + cc%bwood_gain
        cc%bliving = cc%bliving + cc%carbon_gain
        
        if(cc%bliving < 0) then
           cc%bwood    = cc%bwood+cc%bliving
           cc%bliving  = 0
           if (cc%bwood < 0) &
                cc%bwood = 0 ! in principle, that's not conserving carbon
        endif
        
        call update_biomass_pools(cc)
        
        ! reset carbon acculmulation terms
        cc%carbon_gain = 0
        cc%carbon_loss = 0
        cc%bwood_gain  = 0
     endif

     if (cc%status == LEAF_ON) then
        cc%leaf_age = cc%leaf_age + 1.0
        ! limit the maximum leaf age by the leaf time span (reciprocal of leaf 
        ! turnover rate alpha) for given species. alpha is in 1/year, factor of
        ! 365 converts the result to days.
        if (spdata(cc%species)%alpha(CMPT_LEAF) > 0) &
             cc%leaf_age = min(cc%leaf_age,365.0/spdata(cc%species)%alpha(CMPT_LEAF))
     endif
  end do

  if (is_watch_point()) then
     cmass1 = vegn_tile_carbon(vegn)
     write(*,*)"############## vegn_growth #################"
     __DEBUG3__(cmass1-cmass0, cmass1, cmass0)
  endif
end subroutine vegn_growth



! ============================================================================
subroutine plant_respiration(cc, tsoil)
  type(vegn_cohort_type), intent(inout) :: cc
  real, intent(in) :: tsoil
  
  real :: tf,tfs ! thermal inhibition factors for above- and below-ground biomass
  real :: r_leaf, r_vleaf, r_stem, r_root
  real :: Acambium  ! cambium area, m2
  
  integer :: sp ! shorthand for cohort species
  sp = cc%species

  tf  = thermal_inhibition(cc%prog%Tv)
!  tfs = thermal_inhibition(tsoil)  ! original
  tfs = Rm_T_response_function(tsoil) ! Weng 2014-01-14
! With nitrogen model, leaf respiration should be a function of leaf nitrogen
  r_leaf = -mol_C*cc%An_cl*cc%leafarea;
  r_vleaf = spdata(sp)%beta(CMPT_VLEAF) * cc%blv*tf;
  if (do_ppa) then
     ! Stem auto-respiration is proportional to cambium area, not sapwood biomass
     Acambium = PI * cc%DBH * cc%height * 1.2
     r_stem   = spdata(sp)%beta(CMPT_SAPWOOD) * Acambium * tf
  else
     r_stem   = spdata(sp)%beta(CMPT_SAPWOOD) * cc%bsw * tf
  endif
  r_root  = spdata(sp)%beta(CMPT_ROOT) * cc%br*tfs;

  cc%resp = r_leaf + r_vleaf + r_stem + r_root;
  cc%resl = r_leaf;
  cc%resr = r_root;
end subroutine plant_respiration

! ============================================================================
! hourly
! include Nitrogen uptake and partitioning
subroutine vegn_carbon_int_ppa (vegn, tsoil, theta, diag)
  type(vegn_tile_type), intent(inout) :: vegn
  real, intent(in) :: tsoil ! average temperature of soil for soil carbon decomposition, deg K
  real, intent(in) :: theta ! average soil wetness, unitless
  type(diag_buff_type), intent(inout) :: diag

  real :: resp, resl, resr, resg ! respiration terms accumulated for each tree 
  real :: md_alive, md_wood;
  real :: gpp ! gross primary productivity per tile
  real :: gpp_under, npp_under, resp_under
  real :: cgrowth ! carbon spent during this time step for growth, kg C/individual
  real :: deltaBL, deltaBR, deltaBSW ! leaf and fine root carbon tendencies
  real :: b  ! temporary var for several calculations
  integer :: i
  real :: cmass0, cmass1 ! for debug only
  real :: leaf_add, root_add, sapwood_add
  real :: NSC_supply,LR_demand,LR_deficit
  real :: LeafGrowthMin, RootGrowthMin,NSCtarget,v
  real :: LR_growth,WS_growth
  real :: R_days,fNSC,fLFR,fStem,fSeed
! Nitrogen
  logical:: NSN_not_full
  real :: avgNup,turnoverN,pho_N_up
  real :: N_roots0,N_roots
! for cohort output
  integer :: year, month, day, hour, minute, second, maxoutcc
  character(len=100) :: diag_format, diag_format2

  if(is_watch_point()) then
     write(*,*)'#### vegn_carbon_int_ppa input ####'
     __DEBUG2__(tsoil,theta)
     associate(c=>vegn%cohorts(1:vegn%n_cohorts))
     __DEBUG1__(c%bl)
     __DEBUG1__(c%br)
     __DEBUG1__(c%bsw)
     __DEBUG1__(c%bwood)
     __DEBUG1__(c%nsc)
     __DEBUG1__(c%carbon_gain)
     __DEBUG1__(c%carbon_loss)
     __DEBUG1__(c%bwood_gain)
     end associate
     cmass0 = vegn_tile_carbon(vegn) 
  endif

!! Nitrogen uptake
  pho_N_up = 0.5 ! * 8760 * dt_fast_yr  ! need to find out how roots uptake nitrogen
  N_roots0 = 0.3 ! kg C m-2 , roots
  !N_roots0 = 0.2 ! m2/m2, for crown area weighted N uptake 
  N_Roots  = 0.0
  if(vegn%QNmineral > 0.0)then
     do i = 1, vegn%n_cohorts
        associate ( cc => vegn%cohorts(i), &
                    sp => spdata(cc%species) )
!!        There is a problem for the commented scheme: 
!!        the deciduous absorbs too little N in
!!        competing with the evergreen.
!        NSN_not_full = (cc%NSN < 0.08 * cc%nsc).and.(cc%status == LEAF_ON) ! 3*(cc%bl_max/sp%CNleaf0 + cc%br_max/sp%CNroot0)
!        if(NSN_not_full)N_Roots = N_Roots + cc%br * cc%nindivs

!!       A scheme for deciduous to get enough N:
        cc%NSNmax = 0.2 * cc%crownarea  ! just for test
        NSN_not_full = (cc%NSN < cc%NSNmax) ! 5*(cc%bl_max/sp%CNleaf0 + cc%br_max/sp%CNroot0)) !
        if(NSN_not_full)N_Roots = N_Roots + cc%br_max * cc%nindivs
!        if(NSN_not_full)N_Roots = N_Roots + sp%phiRL*sp%LAImax*cc%crownarea/sp%SRA* cc%nindivs   ! just for test
        end associate
     enddo
     ! M-M equation for Nitrogen absoption, check literature
     ! McMurtrie et al. 2012, Ecology & Evolution
     avgNup = pho_N_up /(N_roots0 + N_roots) * vegn%QNmineral ! kg N /kg root C
     vegn%QNmineral = vegn%QNmineral - avgNup * N_roots

     ! Nitrogen uptaken by each cohort, N_uptake
     do i = 1, vegn%n_cohorts
        associate ( cc => vegn%cohorts(i), &
                    sp => spdata(cc%species) )
!        NSN_not_full =(cc%NSN < 0.08 * cc%nsc).and.(cc%status == LEAF_ON) ! 3*(cc%bl_max/sp%CNleaf0 + cc%br_max/sp%CNroot0)
!        if(NSN_not_full) cc%nsn = cc%nsn + cc%br*avgNup

        NSN_not_full = (cc%NSN < cc%NSNmax) ! 5*(cc%bl_max/sp%CNleaf0 + cc%br_max/sp%CNroot0)) !
        if(NSN_not_full)then
            cc%N_uptake=cc%br_max*avgNup
            cc%nsn = cc%nsn + cc%N_uptake ! cc%br_max*avgNup
            cc%N_up_yr=cc%N_up_yr +  cc%N_uptake
        endif

        end associate
     enddo
  endif

  ! update plant carbon for all cohorts
  vegn%npp = 0
  resp = 0 ; resl = 0 ; resr = 0 ; resg = 0 ; gpp = 0
  ! understory GPP, NPP, RESP
  gpp_under = 0; npp_under = 0; resp_under = 0
  do i = 1, vegn%n_cohorts
     associate ( cc => vegn%cohorts(i), &
                 sp => spdata(cc%species) )

     ! that was in eddy_npp_PPA
     call plant_respiration(cc,tsoil)
     cc%gpp  = (cc%An_op - cc%An_cl)*mol_C*cc%leafarea
     cc%nsc  = cc%nsc +  (cc%gpp-cc%resp)*dt_fast_yr 
     cgrowth = 0.0

!    A new scheme for plant growth, modified 9/3/2013 based on Steve's suggestions
     R_days = 5.0
     fNSC = 0.2 * 365. * dt_fast_yr ! 0.2(daily) -->0.2/24 (hourly) 2014-10-22
     fLFR = 0.2 * 365. * dt_fast_yr 
     fStem= dt_fast_yr/sp%tauNSC  ! 0.05
     fSeed= dt_fast_yr/tau_seed   ! 0.05
     NSCtarget = 3.0*(cc%bl_max + cc%br_max)    ! 1.5*(cc%bl_max + cc%br_max)
     
     if(cc%nsc > NSCtarget)then
         v=1.
     else
         v=0.
     endif
!    '365.*dt_fast_yr/R_days' is '1/hours', eg. 1/120
!    the fraction filled per hour
     leaf_add = 0.0; root_add = 0.0; sapwood_add = 0.0
     LR_demand = 0.0; NSC_supply= 0.0
     IF(cc%nsc > 0. .AND. cc%status == LEAF_ON)then
         LR_deficit=max(cc%bl_max+cc%br_max-cc%bl-cc%br,0.0)
         LR_demand =min(fLFR*LR_deficit, fNSC*cc%nsc)
         NSC_supply=max((cc%nsc-0.5*NSCtarget)*fStem,0.0) ! Weng 2014-01-23 for
!                                                           smoothing deltaDBH
     ENDIF
     cc%nsc  = cc%nsc - (LR_demand + NSC_supply)

     cc%resg    = GROWTH_RESP/(1.+GROWTH_RESP) * (LR_demand + NSC_supply)/dt_fast_yr   !     + &
     LR_growth  = LR_demand  /(1.+GROWTH_RESP) 
     WS_growth  = NSC_supply /(1.+GROWTH_RESP)
     cc%carbon_gain = cc%carbon_gain + (LR_growth + WS_growth) !

     cc%resp = cc%resp + cc%resg
     cc%npp  = cc%gpp - cc%resp

!    Turnover of leaves and roots regardless of STATUS according to leaf
!    longevity. Deciduous: 0; Evergreen 0.035/LMa
     deltaBL = cc%bl * sp%alpha(CMPT_LEAF) * dt_fast_yr
     deltaBR = cc%br * sp%alpha(CMPT_ROOT) * dt_fast_yr
     md_alive = deltaBL + deltaBR
     cc%carbon_loss = cc%carbon_loss + md_alive  ! used in diagnostics only
!    update QNleaf and QNroot
     turnoverN = 0.0
     if(cc%bl>0.0)then
         cc%QNleaf = cc%QNleaf - deltaBL/cc%bl*cc%QNleaf
         turnoverN = deltaBL/cc%bl*cc%QNleaf
     endif
     if(cc%br>0.0)then
         cc%QNroot = cc%QNroot - deltaBR/cc%br*cc%QNroot
         turnoverN = turnoverN + deltaBR/cc%br*cc%QNroot
     endif
     cc%bl = cc%bl - deltaBL
     cc%br = cc%br - deltaBR

     ! compute branch and coarse wood losses for tree types
     md_wood = 0.0
     if (cc%species > 1) then ! TODO: add input data selector for woody/grass species
        ! md_wood = 0.6 *cc%bwood * sp%alpha(CMPT_WOOD)*dt_fast_yr
        ! set turnoverable wood biomass as a linear funcion of bl_max (max.foliage
        ! biomass) (Wang, Chuankuan 2006)
        md_wood = MIN(cc%bwood,cc%bl_max) * sp%alpha(CMPT_LEAF)*dt_fast_yr
     endif
     ! Why md_wood is set to 0?
     md_wood = 0.0     
!     cc%bwood = cc%bwood - md_wood
     cc%bwood_gain  = 0.0 ! cc%bwood_gain - md_wood

!    It's not necessary to calculate wood C gain here with the new allocation scheme
!     cc%bwood_gain = cc%bwood_gain + cc%Psw_alphasw * cc%bliving * dt_fast_yr ! - md_wood  !
!    In the new scheme, bwood_gain is driven by the turnover of sapwood, which is driven by NPP allocated to it in turn
!     if (cc%bwood_gain < 0.0) cc%bwood_gain=0.0; ! potential non-conservation ?

     ! add maintenance demand from leaf and root pools to fast soil carbon
     vegn%fast_soil_C = vegn%fast_soil_C + (   fsc_liv *md_alive +    fsc_wood *md_wood)*cc%nindivs
     vegn%slow_soil_C = vegn%slow_soil_C + ((1-fsc_liv)*md_alive + (1-fsc_wood)*md_wood)*cc%nindivs

!!    Nitrogen cycle
     vegn%QNfastSOM = vegn%QNfastSOM +        fsc_liv  * turnoverN * cc%nindivs
     vegn%QNslowSOM = vegn%QNslowSOM + (1.0 - fsc_liv) * turnoverN * cc%nindivs

     ! accumulate tile-level NPP and GPP
     vegn%npp = vegn%npp + cc%npp * cc%nindivs
     gpp      = gpp      + cc%gpp * cc%nindivs

!    Weng 2015-09-18
     cc%annualGPP = cc%annualGPP + cc%gpp *dt_fast_yr
     cc%annualNPP = cc%annualNPP + cc%npp *dt_fast_yr

     ! accumulate respiration terms for tile-level reporting
     resp = resp + cc%resp*cc%nindivs ; resl = resl + cc%resl*cc%nindivs
     resr = resr + cc%resr*cc%nindivs ; resg = resg + cc%resg*cc%nindivs

     ! accumulate tile-level understory NPP, GPP, and RESP
     if(cc%layer > 1) then
         npp_under  = npp_under  + cc%npp * cc%nindivs
         gpp_under  = gpp_under  + cc%gpp * cc%nindivs
         resp_under = resp_under + cc%resp* cc%nindivs
     endif
     ! increment tha cohort age
     cc%age = cc%age + dt_fast_yr

     cc%resg = 0.0
     end associate
  enddo 

  if(is_watch_point()) then
     associate(c=>vegn%cohorts(1:vegn%n_cohorts))
     write(*,*)'#### vegn_carbon_int_ppa output ####'
     __DEBUG1__(c%species)
     __DEBUG5__(c%bl, c%br, c%bsw, c%bwood, c%nsc)
     __DEBUG2__(c%An_op, c%An_cl)
     __DEBUG2__(c%npp, c%gpp)
     __DEBUG2__(c%resp, c%resg)
     __DEBUG2__(c%carbon_gain, c%bwood_gain)
     end associate
     cmass1 = vegn_tile_carbon(vegn)
     __DEBUG3__(cmass1-cmass0-(gpp-resp)*dt_fast_yr,cmass1,cmass0)
     write(*,*)'#### end of vegn_carbon_int_ppa output ####'
  endif

  ! update soil carbon
!  call Dsdt(vegn, diag, tsoil, theta)
  call Dsdtmicrobes(vegn, diag, tsoil, theta)
  ! NEP is equal to NNP minus soil respiration
  vegn%nep = vegn%npp - vegn%rh

  call update_soil_pools(vegn)
  vegn%age = vegn%age + dt_fast_yr;

! ------ diagnostic section
  call send_tile_data(id_gpp,gpp,diag)
  call send_tile_data(id_npp,vegn%npp,diag)
  call send_tile_data(id_nep,vegn%nep,diag)
  call send_tile_data(id_litter,vegn%litter,diag)
  call send_tile_data(id_resp, resp, diag)
  call send_tile_data(id_resl, resl, diag)
  call send_tile_data(id_resr, resr, diag)
  call send_tile_data(id_resg, resg, diag)
  call send_tile_data(id_soilt,tsoil,diag)
  call send_tile_data(id_theta,theta,diag)
!  cc => vegn%cohorts(1)
!  call send_tile_data(id_nsc_c1,cc%nsc,diag)
  associate (cc => vegn%cohorts(1))
  call send_tile_data(id_gpp_c1, cc%gpp,diag)
  call send_tile_data(id_npp_c1, cc%npp,diag)
  call send_tile_data(id_resp_c1,cc%resp, diag)
  call send_tile_data(id_nsc_c1, cc%nsc, diag)
  end associate
! understory gpp, npp, resp
  call send_tile_data(id_gpp_under, gpp_under,  diag)
  call send_tile_data(id_npp_under, npp_under,  diag)
  call send_tile_data(id_resp_under,resp_under, diag)
!! cohorts GPP,An_op,An_cl, NPP, Resp
  diag_format  ='(i4,"-",i2.2,"-",i2.2,1x,i2.2,":",i2.2,":",i2.2,2x,", ", i4,",", 60(E10.4,","))'
  diag_format2 ='(i4,"-",i2.2,"-",i2.2,1x,i2.2,":",i2.2,":",i2.2,2x,", ", i4,",", 60(i4,","))'
  call get_date(lnd%time, year, month, day, hour, minute, second)
  maxoutcc = min(60,vegn%n_cohorts)
  associate (cc => vegn%cohorts(1:vegn%n_cohorts))
  if(year >1940 .AND. year<2251)then
!     id_cc_An_op,id_cc_An_cl,id_cc_gpp,id_cc_npp,id_cc_resp,id_cc_ca,id_cc_la
!     cc%gpp  = (cc%An_op - cc%An_cl)*mol_C*cc%leafarea
     write(id_cc_PFT,diag_format2) year,month,day,hour,minute,second, &
                                 maxoutcc, (cc(i)%species, i=1, maxoutcc)

!     write(id_cc_An_op,diag_format) year,month,day,hour,minute,second, &
!                                 maxoutcc, (cc(i)%An_op, i=1, maxoutcc)

!     write(id_cc_An_cl,diag_format) year,month,day,hour,minute,second, &
!                                 maxoutcc, (cc(i)%An_cl, i=1, maxoutcc)

     write(id_cc_gpp,diag_format)   year,month,day,hour,minute,second, &
                                 maxoutcc, (cc(i)%gpp, i=1, maxoutcc)

     write(id_cc_npp,diag_format)   year,month,day,hour,minute,second, &
                                 maxoutcc, (cc(i)%npp, i=1, maxoutcc)

!     write(id_cc_resp,diag_format)  year,month,day,hour,minute,second, &
!                                maxoutcc, (cc(i)%resp, i=1, maxoutcc)

     write(id_cc_ca,diag_format)    year,month,day,hour,minute,second, &
                                maxoutcc, (cc(i)%crownarea, i=1, maxoutcc)

     write(id_cc_la,diag_format)    year,month,day,hour,minute,second, &
                                maxoutcc, (cc(i)%leafarea, i=1,maxoutcc)

     write(id_cc_Nup,diag_format)    year,month,day,hour,minute,second, &
                                maxoutcc, (cc(i)%N_uptake, i=1,maxoutcc)
  endif
  end associate
  
end subroutine vegn_carbon_int_ppa

!========================================================================
! Starvation due to low NSC
subroutine vegn_starvation_ppa (vegn)
  type(vegn_tile_type), intent(inout) :: vegn

  ! ---- local vars
  real :: loss_alive,loss_wood
  real :: deathrate ! mortality rate, 1/year
  real :: deadtrees ! number of trees that died over the time step
  real :: lossN_wood, lossN_alive
  integer :: i, k
  type(vegn_cohort_type), pointer :: cc(:) ! array to hold new cohorts

  do i = 1, vegn%n_cohorts
     associate ( cc => vegn%cohorts(i)   , &
                 sp => spdata(cc%species)  )

     ! Mortality due to starvation
    deathrate = 0.0
!    if (cc%bsw<0 .or. cc%nsc < 0.00001*cc%bl_max .OR.(cc%layer >1 .and. sp%lifeform ==0)) then
    if (cc%bsw<0 .or. cc%nsc < 0.00001*cc%bl_max) then
         deathrate = 1.0
         deadtrees = cc%nindivs * deathrate !individuals / m2
!        add dead C from leaf and root pools to fast soil carbon
         loss_wood  = deadtrees * cc%bwood
         loss_alive = deadtrees * (cc%bl + cc%br + cc%bsw + cc%bseed + cc%nsc)

         vegn%fast_soil_C = vegn%fast_soil_C +    fsc_liv *loss_alive +   &
                            fsc_wood*loss_wood
         vegn%slow_soil_C = vegn%slow_soil_C + (1-fsc_liv)*loss_alive +   &
                            (1-fsc_wood)*loss_wood

!!        Nitrogen to soil pools
         lossN_wood = deadtrees * cc%QNwood
         lossN_alive= deadtrees * (cc%QNleaf+cc%QNroot+cc%QNsw+cc%QNseed+cc%nsn)
         vegn%QNfastSOM = vegn%QNfastSOM +    fsc_liv  *lossN_alive +   &
                                              fsc_wood *lossN_wood
         vegn%QNslowSOM = vegn%QNslowSOM +(1.-fsc_liv) *lossN_alive +   &
                                          (1.-fsc_wood)*lossN_wood

!        update cohort individuals
         cc%nindivs = cc%nindivs*(1.0-deathrate)

!        water goes to other cohorts
!         if(i>1)then
!            vegn%cohorts(i-1)%prog%Wl=vegn%cohorts(i-1)%prog%Wl + &
!                                  cc%prog%Wl*deadtrees/vegn%cohorts(i-1)%nindivs
!            vegn%cohorts(i-1)%prog%Ws=vegn%cohorts(i-1)%prog%Ws + &
!                                  cc%prog%Ws*deadtrees/vegn%cohorts(i-1)%nindivs
!         elseif(i==1 .and. vegn%n_cohorts>1)then
!            vegn%cohorts(i+1)%prog%Wl=vegn%cohorts(i+1)%prog%Wl + &
!                                  cc%prog%Wl*deadtrees/vegn%cohorts(i+1)%nindivs
!            vegn%cohorts(i+1)%prog%Ws=vegn%cohorts(i+1)%prog%Ws + &
!                                  cc%prog%Ws*deadtrees/vegn%cohorts(i+1)%nindivs
!         endif


     ! for budget tracking - temporary
     vegn%ssc_in = vegn%ssc_in + fsc_liv*loss_alive + fsc_wood *loss_wood
     vegn%veg_out = vegn%veg_out + loss_alive + loss_wood

     else
         deathrate = 0.0
     endif

     end associate
  enddo

  ! calculate the number of remaining cohorts
  k = 0
  do i = 1, vegn%n_cohorts
     if (vegn%cohorts(i)%nindivs > 0.0) k=k+1
  enddo

  if (k==0) call error_mesg('vegn_starvation_ppa','All cohorts died',WARNING)
  if (is_watch_point()) then
     write(*,*)'###### vegn_nat_mortality_ppa #######'
     __DEBUG1__(vegn%cohorts%nindivs)
     __DEBUG1__(k)
  endif

  ! exclude cohorts that have zero individuals
  if (k < vegn%n_cohorts) then
     allocate(cc(k))
     k=0
     do i = 1,vegn%n_cohorts
        if (vegn%cohorts(i)%nindivs > 0.0) then
           k=k+1
           cc(k) = vegn%cohorts(i)
        endif
     enddo

     vegn%n_cohorts = k
     deallocate (vegn%cohorts)
     vegn%cohorts=>cc

!! count canopy layers
!     vegn%n_canopycc = 0
!     do i = 1,vegn%n_cohorts
!        if (vegn%cohorts(i)%layer ==1) then
!           vegn%n_canopycc=vegn%n_canopycc+1
!        else
!           exit
!        endif
!     enddo

  endif

end subroutine vegn_starvation_ppa

! ==============================================================================
! updates cohort vegetation structure, biomass pools, LAI, SAI, and height  
! using accumulated carbon_gain
subroutine biomass_allocation_ppa(cc)
  type(vegn_cohort_type), intent(inout) :: cc
  
  real :: CSAtot ! total cross section area, m2
  real :: CSAsw  ! Sapwood cross sectional area, m2
  real :: CSAwd  ! Heartwood cross sectional area, m2
  real :: DBHwd  ! diameter of heartwood at breast height, m
  real :: BSWmax ! max sapwood biomass, kg C/individual
  real :: G_LFR  ! amount of carbon spent on leaf and root growth
  real :: deltaSeed
  real :: deltaBL, deltaBR ! tendencies of leaf and root biomass, kgC/individual
  real :: deltaBSW ! tendency of sapwood biomass, kgC/individual
  real :: deltaBwood ! tendency of wood biomass, kgC/individual
  real :: deltaDBH ! tendency of breast height diameter, m
  real :: deltaCA ! tendency of crown area, m2/individual
  real :: deltaHeight ! tendency of vegetation height
  real :: deltaNsw    ! Nitrogen from SW to HW
  real :: sw2nsc = 0.0 ! conversion of sapwood to non-structural carbon
  real :: b,BL_u,BL_c
  real :: alphaBL, alphaBR
  real :: DBHtp
! for nitrogen
  real :: N_supply, N_demand,fNr,Nsupplyratio,extraQNsw

  logical :: do_differ_allometry = .False.

  DBHtp = 0.8
  fNr   = 0.25
!! for testing H2 in LM3-PPA paper
!  do_differ_allometry = .TRUE.

  associate (sp => spdata(cc%species)) ! F2003

  ! TODO: what if carbon_gain is not 0, but leaves are OFF (marginal case? or
  ! typical in lm3?)
  if (cc%status == LEAF_ON) then
     ! calculate the carbon spent on growth of leaves and roots
     G_LFR    = max(0.0, min(cc%bl_max+cc%br_max-cc%bl-cc%br,  &
                            (1.-Wood_fract_min)*cc%carbon_gain))
     ! and distribute it between roots and leaves
     deltaBL  = min(G_LFR, max(0.0, &
          (G_LFR*cc%bl_max + cc%bl_max*cc%br - cc%br_max*cc%bl)/(cc%bl_max + cc%br_max) &
          ))
     deltaBR  = G_LFR - deltaBL
     ! calculate carbon spent on growth of sapwood growth
     if(cc%layer == 1.AND. cc%age > sp%maturalage)then
         deltaSeed=      sp%v_seed * (cc%carbon_gain - G_LFR)
         deltaBSW = (1.0-sp%v_seed)* (cc%carbon_gain - G_LFR)
     else
         deltaSeed= 0.0 
         deltaBSW = cc%carbon_gain - G_LFR
     endif
!    Specially for grasses, temporary
     if(sp%lifeform ==0 )then
         deltaSeed = deltaSeed + 0.15*G_LFR
         G_LFR     = 0.85 * G_LFR
         deltaBR  =  0.85 * deltaBR  
         deltaBL  =  0.85 * deltaBL
     endif
!!    Nitrogen effects on allocations between wood and leaves+roots
!
!!    Nitrogen demand by leaves, roots, and seeds (Their C/N ratios are fixed.)
     N_demand = deltaBL/sp%CNleaf0 + deltaBR/sp%CNroot0 + deltaSeed/sp%CNseed0 ! + deltaBSW/sp%CNsw0
!!    Nitrogen available for all tisues, including wood
     N_supply= fNr*cc%NSN     
!!    same ratio reduction for leaf, root, and seed if(N_supply < N_demand)
     IF(N_demand > 0.0)then
         Nsupplyratio = MIN(1.0, N_supply / N_demand)
     ELSE
         Nsupplyratio = 1.0
     ENDIF
     deltaBSW = deltaBSW + (1.0 - Nsupplyratio)*(deltaBL+deltaBR+deltaSeed)
     deltaBR  =  Nsupplyratio * deltaBR  
     deltaBL  =  Nsupplyratio * deltaBL
     deltaSeed=  Nsupplyratio * deltaSeed

!    update biomass pools
     cc%bl     = cc%bl    + deltaBL  ! updated in vegn_int_ppa
     cc%br     = cc%br    + deltaBR
     cc%bsw    = cc%bsw   + deltaBSW
     cc%bseed  = cc%bseed + deltaSeed
!!    update nitrogen pools, Nitrogen allocation
     cc%NSN    = cc%NSN - N_supply
     cc%QNleaf = cc%QNleaf + deltaBL   /sp%CNleaf0
     cc%QNroot = cc%QNroot + deltaBR   /sp%CNroot0
     cc%QNseed = cc%QNseed + deltaSeed/sp%CNseed0
     cc%QNsw   = cc%QNsw   + MAX((N_supply - N_demand),0.0)
!    Return excessiive Nitrogen in SW back to NSN
     extraQNsw = MAX(0.0, cc%QNsw - cc%bsw/sp%CNsw0)
     cc%NSN    = cc%NSN  + extraQNsw ! MAX(0.0, cc%QNsw - cc%bsw/sp%CNsw0)
     cc%QNsw   = cc%QNsw - extraQNsw ! MAX(0.0, cc%QNsw - cc%bsw/sp%CNsw0)

!    accumulated C allocated to leaf, root, and wood
     cc%NPPleaf = cc%NPPleaf + deltaBL
     cc%NPProot = cc%NPProot + deltaBR
     cc%NPPwood = cc%NPPwood + deltaBSW

!    update breast height diameter given increase of bsw
     deltaDBH     = deltaBSW / (sp%thetaBM * sp%alphaBM * cc%DBH**(sp%thetaBM-1))
     deltaHeight  = sp%thetaHT * sp%alphaHT * cc%DBH**(sp%thetaHT-1) * deltaDBH
     if(do_differ_allometry) then
         deltaCA  = sp%thetaCA * sp%alphaCA *     &
               (1./(exp(15.*(cc%DBH-DBHtp))+1.))  *     &     
               cc%DBH**(sp%thetaCA-1) * deltaDBH
     else
         deltaCA  = sp%thetaCA * sp%alphaCA * cc%DBH**(sp%thetaCA-1) * deltaDBH
     endif

     cc%DBH       = cc%DBH       + deltaDBH
     cc%height    = cc%height    + deltaHeight
     cc%crownarea = cc%crownarea + deltaCA

     ! calculate DBH, BLmax, BRmax, BSWmax, bwood_gain using allometric relationships
     ! Weng 2012-01-31 update_bio_living_fraction
!    conversion of sapwood to heartwood
     if(sp%lifeform>0)then
        CSAsw    = sp%alphaCSASW * cc%DBH**sp%thetaCSASW
        CSAtot   = PI * (cc%DBH/2.0)**2
        CSAwd    = max(0.0, CSAtot - CSAsw)
        DBHwd    = 2*sqrt(CSAwd/PI)
        BSWmax   = sp%alphaBM * (cc%DBH**sp%thetaBM - DBHwd**sp%thetaBM)
        deltaBwood = max(cc%bsw - BSWmax, 0.0)
        deltaNsw   = deltaBwood/cc%bsw *cc%QNsw
        cc%bwood   = cc%bwood + deltaBwood
        cc%bsw     = cc%bsw   - deltaBwood
!        write(*,*)'Conversion sw->wood',cc%bsw, deltaBwood
!!       Nitrogen from sapwood to heart wood
        cc%QNsw   = cc%QNsw - deltaNsw
        cc%QNwood = cc%QNwood + deltaNsw
     endif

!    update bl_max and br_max daily
     BL_u = sp%LMA*cc%crownarea*(1.0-sp%internal_gap_frac)* &
                 sp%underLAImax
!                MAX(understory_lai_factor*sp%LAImax,1.0)

     BL_c = sp%LMA * sp%LAImax * cc%crownarea * &
            (1.0-sp%internal_gap_frac)
     if (cc%layer > 1 .and. cc%firstlayer == 0) then ! changed back, Weng 2014-01-23
        cc%bl_max = BL_u
        cc%br_max = 0.8*cc%bl_max/(sp%LMA*sp%SRA)
     else
        if(sp%lifeform == 0) then
           cc%bl_max = BL_c
        else
           cc%bl_max = BL_u + min(cc%topyear/5.0,1.0)*(BL_c - BL_u)
        endif
        cc%br_max = sp%phiRL*cc%bl_max/(sp%LMA*sp%SRA)
     endif

     if(is_watch_point()) then
        __DEBUG4__(cc%bl, cc%br, cc%bl_max, cc%br_max)
        __DEBUG5__(cc%carbon_gain, G_LFR, deltaBL, deltaBR, deltaBSW)
     endif

  else  ! if(cc%status == LEAF_OFf)then
     cc%nsc = cc%nsc + cc%carbon_gain

  endif ! "cc%status == LEAF_ON"

  ! reset carbon acculmulation terms
  cc%carbon_gain = 0
  cc%carbon_loss = 0
  cc%bwood_gain  = 0

  end associate ! F2003
end subroutine biomass_allocation_ppa

! ============================================================================
! the reproduction of each canopy cohort, yearly time step
! calculate the new cohorts added in this step and states:
! tree density, DBH, woddy and living biomass
subroutine vegn_reproduction_ppa (vegn)
  type(vegn_tile_type), intent(inout) :: vegn

! ---- local vars
  type(vegn_cohort_type), pointer :: parent, cc ! parent and child cohort pointers
  type(vegn_cohort_type), pointer :: ccold(:)   ! pointer to old cohort array
  real :: failed_seeds,grassdensity, N_failedseed !, prob_g, prob_e
  logical :: invasion = .FALSE.
  integer :: newcohorts,mutantcc,mutantyr ! number of new cohorts to be created
  integer :: i, k ! cohort indices

!! parameter values

!  invasion = .TRUE.
! Check if reproduction happens
  newcohorts = 0
  do i = 1, vegn%n_cohorts
     if (cohort_can_reproduce(vegn%cohorts(i))) newcohorts=newcohorts + 1
  enddo
  if (newcohorts == 0) return ! do nothing if no cohorts are ready for reproduction

  ccold => vegn%cohorts ! keep old cohort information
  allocate(vegn%cohorts(vegn%n_cohorts+newcohorts))
  vegn%cohorts(1:vegn%n_cohorts) = ccold(1:vegn%n_cohorts) ! copy old cohort information
  deallocate (ccold)

  ! set up new cohorts
  k = vegn%n_cohorts
  do i = 1,vegn%n_cohorts
    if (.not.cohort_can_reproduce(vegn%cohorts(i))) cycle ! nothing to do

    k = k+1 ! increment new cohort index
    
    ! update child cohort parameters
    associate (cc     => vegn%cohorts(k), &  ! F2003
               parent => vegn%cohorts(i), &
               sp     => spdata(parent%species))
    ! copy all parent values into the new cohort then change some of them below
    cc         = parent

    cc%status  = LEAF_OFF
    cc%firstlayer = 0
    cc%topyear = 0.0
    cc%age     = 0.0
    cc%bl      = sp%seedlingsize * bseed_distr(CMPT_LEAF)
    cc%br      = sp%seedlingsize * bseed_distr(CMPT_ROOT)
    cc%bsw     = sp%seedlingsize * bseed_distr(CMPT_SAPWOOD) ! sp%seedlingsize*0.1
    cc%bwood   = sp%seedlingsize * bseed_distr(CMPT_WOOD)    ! sp%seedlingsize*0.05
    cc%nsc     = sp%seedlingsize - cc%bsw ! sp%seedlingsize * bseed_distr(CMPT_NSC)     ! sp%seedlingsize*0.85
    cc%blv     = sp%seedlingsize * bseed_distr(CMPT_VLEAF)
    cc%bseed   = 0.0
!!  for Evergreen, fixing the error of tempertature out of range
!    if(sp%phenotype==1)then
!        cc%status  = LEAF_ON
!        cc%bl      = 0.05*cc%bsw
!        cc%br      = 0.05*cc%bsw
!        cc%bsw     = 0.9 *cc%bsw
!    endif

!!   Nitrogen pools
    cc%QNleaf  = cc%bl/sp%CNleaf0
    cc%QNroot  = cc%br/sp%CNroot0
    cc%QNsw    = cc%bsw/sp%CNsw0
    cc%QNwood  = cc%bwood/sp%CNwood0
    cc%NSN     = MAX(0.001,cc%QNseed*(sp%seedlingsize/parent%bseed) - &
                 (cc%QNleaf+cc%QNroot+cc%QNsw+cc%QNwood))

    cc%QNseed  = 0.0

    cc%nindivs = parent%bseed*parent%nindivs * sp%prob_g * sp%prob_e   &
                 /(sp%seedlingsize*sum(bseed_distr(:)))
!    write(*,*)cc%nindivs
    failed_seeds = (1.-sp%prob_g*sp%prob_e) * parent%bseed * parent%nindivs

    vegn%litter = vegn%litter + failed_seeds
    vegn%fast_soil_C = vegn%fast_soil_C +        fsc_liv *failed_seeds
    vegn%slow_soil_C = vegn%slow_soil_C + (1.0 - fsc_liv)*failed_seeds
    vegn%fsc_in  = vegn%fsc_in  + failed_seeds
    vegn%veg_out = vegn%veg_out + failed_seeds

!!   Nitrogen of seeds to soil SOMs
    N_failedseed= (1.-sp%prob_g*sp%prob_e) * parent%QNseed * parent%nindivs
    vegn%QNfastSOM   = vegn%QNfastSOM   +        fsc_liv * N_failedseed
    vegn%QNslowSOM   = vegn%QNslowSOM   + (1.0 - fsc_liv)* N_failedseed


!   reset bseed
    parent%bseed = 0.0
    parent%QNseed= 0.0

    call init_cohort_allometry_ppa(cc)
    cc%carbon_gain = 0.0
    cc%carbon_loss = 0.0
    cc%bwood_gain  = 0.0
    cc%bliving     = cc%br + cc%bl + cc%bsw + cc%blv
    cc%DBH_ys      = cc%DBH
    cc%BM_ys       = cc%bsw + cc%bwood
    cc%gpp          = 0.0
    cc%npp          = 0.0  
    cc%resp         = 0.0
    cc%resl         = 0.0
    cc%resr         = 0.0
    cc%resg         = 0.0
    cc%npp_previous_day     = 0.0
    cc%npp_previous_day_tmp = 0.0
    cc%leaf_age     = 0.0
   
    cc%annualGPP = 0.0
    cc%annualNPP = 0.0 
    cc%NPPleaf = 0.0
    cc%NPProot = 0.0
    cc%NPPwood = 0.0

    cc%N_up_yr = 0.0 ! annual cohort N uptake

    ! we assume that the newborn cohort is dry; since nindivs of the parent
    ! doesn't change we don't need to do anything with its Wl and Ws to 
    ! conserve water (since Wl and Ws are per individual)
    cc%prog%Wl = 0 ; cc%prog%Ws = 0
    ! TODO: make sure that energy is conserved in reproduction
    cc%prog%Tv = 278.0 ! parent%prog%Tv
    
    end associate   ! F2003
  enddo
  
  vegn%n_cohorts = k

!! reset grass initial size and density
!  do i = 1,vegn%n_cohorts
!     associate (cc     => vegn%cohorts(i), &  ! F2003
!                sp     => spdata(cc%species))
!!    set stem diameter and crown area for grasses at the begining of the growing season
!     if(sp%lifeform == 0 .and. cc%age>0.0)then
!         grassdensity = cc%nindivs
!         cc%nindivs = grassdensity*cc%nsc /sp%seedlingsize
!         cc%br  = grassdensity/cc%nindivs * cc%br
!         cc%nsc = grassdensity*cc%nsc/cc%nindivs
!         cc%bsw = bseed_distr(CMPT_SAPWOOD)*cc%nsc  ! for setting up a initial size
!         cc%nsc = cc%nsc - cc%bsw
!         call init_cohort_allometry_ppa(cc);
!         cc%bliving     = cc%br + cc%bl + cc%bsw + cc%blv
!         cc%DBH_ys      = cc%DBH
!         cc%BM_ys       = cc%bsw + cc%bwood
!         cc%leaf_age = 0.0
!     endif
!     end associate
!  enddo

! switch to other species, Weng 2013-11-01
  if(invasion)then
     mutantyr = 400
     mutantcc = 3
     if(vegn%n_years == mutantyr)then
         newcohorts = 1
         ccold => vegn%cohorts ! keep old cohort information
         allocate(vegn%cohorts(vegn%n_cohorts+newcohorts))
         ccold(mutantcc)%nindivs = ccold(mutantcc)%nindivs/(1+newcohorts)
         vegn%cohorts(1:vegn%n_cohorts) = ccold(1:vegn%n_cohorts) ! copy old cohort information
         do i=1,newcohorts
             vegn%cohorts(vegn%n_cohorts+i)         = ccold(mutantcc)
             vegn%cohorts(vegn%n_cohorts+i)%species = ccold(mutantcc)%species - i
         enddo
         deallocate (ccold)
         vegn%n_cohorts = vegn%n_cohorts + newcohorts
     endif
  endif  

end subroutine vegn_reproduction_ppa

! ============================================================================
function cohort_can_reproduce(cc); logical cohort_can_reproduce
  type(vegn_cohort_type), intent(in) :: cc
  
  cohort_can_reproduce = (cc%layer == 1 .and. &
                          cc%age   > spdata(cc%species)%maturalage .and. &
                          cc%bseed > 0.0)
end function 


! ============================================================================
! Merge similar cohorts in a tile
subroutine vegn_mergecohorts_ppa(vegn)
  type(vegn_tile_type), intent(inout) :: vegn

! ---- local vars
  type(vegn_cohort_type), pointer :: cc(:) ! array to hold new cohorts
  logical :: merged(vegn%n_cohorts)        ! mask to skip cohorts that were already merged
  real, parameter :: mindensity = 1.0E-6
  real :: loss_alive,loss_wood
!  real            :: loss_wood, loss_alive
  integer :: i,j,k

  allocate(cc(vegn%n_cohorts))

  merged(:)=.FALSE. ; k = 0
  do i = 1, vegn%n_cohorts 
     if(merged(i)) cycle ! skip cohorts that were already merged
     k = k+1
     cc(k) = vegn%cohorts(i)
     ! try merging the rest of the cohorts into current one
     do j = i+1, vegn%n_cohorts
        if (merged(j)) cycle ! skip cohorts that are already merged
        if (cohorts_can_be_merged(vegn%cohorts(j),cc(k))) then
           call merge_cohorts(vegn%cohorts(j),cc(k))
           merged(j) = .TRUE.
        endif
     enddo
  enddo
  ! at this point, k is the number of new cohorts
  vegn%n_cohorts = k
  deallocate(vegn%cohorts)
  vegn%cohorts=>cc
!  write(*,*),1,cc(1)%bsw,cc(1)%nsc,cc(1)%bl_max,cc(1)%layer
  ! note that the size of the vegn%cohorts may be larger than vegn%n_cohorts

!! count canopy layers
!     vegn%n_canopycc = 0
!     do i = 1,vegn%n_cohorts
!        if (vegn%cohorts(i)%layer ==1) then
!           vegn%n_canopycc=vegn%n_canopycc+1
!        else
!           exit
!        endif
!     enddo

end subroutine vegn_mergecohorts_ppa

! ============================================================================
! kill low density cohorts, a new function seperated from vegn_mergecohorts_ppa
! Weng, 2014-07-22
subroutine kill_lowdensity_cohorts(vegn)
  type(vegn_tile_type), intent(inout) :: vegn

! ---- local vars
  type(vegn_cohort_type), pointer :: cc(:) ! array to hold new cohorts
  logical :: merged(vegn%n_cohorts)        ! mask to skip cohorts that were already merged
  real, parameter :: mindensity = 5.0E-5
  real :: loss_alive,loss_wood
  real :: lossN_alive,lossN_wood
!  real            :: loss_wood, loss_alive
  integer :: i,j,k


 ! calculate the number of cohorts with indivs>mindensity
  k = 0
  do i = 1, vegn%n_cohorts
     if (vegn%cohorts(i)%nindivs >  mindensity) k=k+1
  enddo
  if (k==0) call error_mesg('kill_lowdensity_cohorts','All cohorts died',WARNING)

  ! exclude cohorts that have low individuals
  if (k < vegn%n_cohorts) then
     allocate(cc(k))
     k=0
     do i = 1,vegn%n_cohorts
        if (vegn%cohorts(i)%nindivs > mindensity) then
           k=k+1
           cc(k) = vegn%cohorts(i)
        else
!          add dead C from leaf and root pools to fast soil carbon
           loss_wood  = vegn%cohorts(i)%nindivs * vegn%cohorts(i)%bwood
           loss_alive = vegn%cohorts(i)%nindivs * &
                     (vegn%cohorts(i)%bl   + &
                      vegn%cohorts(i)%br   + &
                      vegn%cohorts(i)%bsw  + &
                      vegn%cohorts(i)%bseed+ &
                      vegn%cohorts(i)%nsc)
           vegn%fast_soil_C = vegn%fast_soil_C +    fsc_liv *loss_alive +   &
                            fsc_wood*loss_wood
           vegn%slow_soil_C = vegn%slow_soil_C + (1-fsc_liv)*loss_alive +   &
                            (1-fsc_wood)*loss_wood

!!          Nitrogen to soil SOMs
           lossN_wood  = vegn%cohorts(i)%nindivs * vegn%cohorts(i)%QNwood
           lossN_alive = vegn%cohorts(i)%nindivs * &
                     (vegn%cohorts(i)%QNleaf + &
                      vegn%cohorts(i)%QNroot + &
                      vegn%cohorts(i)%QNsw   + &
                      vegn%cohorts(i)%QNseed + &
                      vegn%cohorts(i)%NSN)
           vegn%QNfastSOM = vegn%QNfastSOM +    fsc_liv  * lossN_alive +   &
                                                fsc_wood * lossN_wood
           vegn%QNslowSOM = vegn%QNslowSOM + (1.-fsc_liv)* lossN_alive +   &
                                          (1. - fsc_wood)* lossN_wood

        endif
     enddo
     vegn%n_cohorts = k
     deallocate (vegn%cohorts)
     vegn%cohorts=>cc
  endif

end subroutine kill_lowdensity_cohorts

! ============================================================================
! Merge non-canopy cohorts in a tile
subroutine merge_noncanopy_cohorts(vegn)
  type(vegn_tile_type), intent(inout) :: vegn

! ---- local vars
  type(vegn_cohort_type), pointer :: cc(:) ! array to hold new cohorts
  logical :: merged(vegn%n_cohorts)        ! mask to skip cohorts that were already merged
  real, parameter :: mindensity = 1.0E-10
  real :: loss_alive,loss_wood
!  real            :: loss_wood, loss_alive
  integer :: i,j,k

  allocate(cc(vegn%n_cohorts))

  merged(:)=.FALSE. ; k = 0
  do i = 1, vegn%n_cohorts 
     if(merged(i)) cycle ! skip cohorts that were already merged
     k = k+1
     cc(k) = vegn%cohorts(i)
     ! try merging the rest of the cohorts into current one
     if(cc(k)%layer==1)cycle
     do j = i+1, vegn%n_cohorts
        if (merged(j)) cycle ! skip cohorts that are already merged
        if (cohorts_can_be_merged(vegn%cohorts(j),cc(k))) then
           call merge_cohorts(vegn%cohorts(j),cc(k))
           merged(j) = .TRUE.
        endif
     enddo
  enddo
  ! at this point, k is the number of new cohorts
  vegn%n_cohorts = k
  deallocate(vegn%cohorts)
  vegn%cohorts=>cc

end subroutine merge_noncanopy_cohorts

! ============================================================================
subroutine merge_cohorts(c1,c2)
  type(vegn_cohort_type), intent(in) :: c1
  type(vegn_cohort_type), intent(inout) :: c2
  
  real :: x1, x2 ! normalized relative weights
  real :: HEAT1, HEAT2 ! heat stored in respective canopies

  x1 = c1%nindivs/(c1%nindivs+c2%nindivs)
  x2 = 1.0-x1
  ! update number of individuals in merged cohort
  c2%nindivs = c1%nindivs+c2%nindivs
#define __MERGE__(field) c2%field = x1*c1%field + x2*c2%field
  HEAT1 = (clw*c1%prog%Wl + csw*c1%prog%Ws + c1%mcv_dry)*(c1%prog%Tv-tfreeze)
  HEAT2 = (clw*c2%prog%Wl + csw*c2%prog%Ws + c2%mcv_dry)*(c2%prog%Tv-tfreeze)

  __MERGE__(prog%Wl)
  __MERGE__(prog%Ws)
 
  __MERGE__(bl)      ! biomass of leaves, kg C/indiv
  __MERGE__(blv)     ! biomass of virtual leaves (labile store), kg C/indiv
  __MERGE__(br)      ! biomass of fine roots, kg C/indiv
  __MERGE__(bsw)     ! biomass of sapwood, kg C/indiv
  __MERGE__(bwood)   ! biomass of heartwood, kg C/indiv
  __MERGE__(bseed)   ! future progeny, kgC/indiv
  __MERGE__(nsc)     ! non-structural carbon, kgC/indiv
  __MERGE__(bliving) ! leaves, fine roots, and sapwood biomass, kgC/indiv
  __MERGE__(dbh)     ! diameter at breast height
  __MERGE__(height)  ! cohort height
  __MERGE__(crownarea)   ! area of cohort crown

  __MERGE__(age)     ! age of individual
  __MERGE__(carbon_gain) ! carbon gain during a day, kg C/indiv
  __MERGE__(carbon_loss) ! carbon loss during a day, kg C/indiv [diag only]
  __MERGE__(bwood_gain)  ! heartwood gain during a day, kg C/indiv
  __MERGE__(topyear)     ! the years that a cohort in canopy layer

!! Merge Nitrogen pools
  __MERGE__(QNleaf)
  __MERGE__(QNroot)
  __MERGE__(QNsw)
  __MERGE__(QNwood)
  __MERGE__(QNseed)
  __MERGE__(NSN)
!  write(*,*)c2%NSN
  ! calculate the resulting dry heat capacity
  c2%leafarea = leaf_area_from_biomass(c2%bl, c2%species, c2%layer, c2%firstlayer)
  c2%mcv_dry = max(mcv_min,mcv_lai*c2%leafarea)
  ! update canopy temperature -- just merge it based on area weights if the heat 
  ! capacities are zero, or merge it based on the heat content if the heat contents
  ! are non-zero
  if(HEAT1==0.and.HEAT2==0) then
     __MERGE__(prog%Tv)
  else
     c2%prog%Tv = (HEAT1*x1+HEAT2*x2) / &
          (clw*c2%prog%Wl + csw*c2%prog%Ws + c2%mcv_dry) + tfreeze
  endif
#undef  __MERGE__
  
end subroutine merge_cohorts

! ============================================================================
function cohorts_can_be_merged(c1,c2); logical cohorts_can_be_merged
   type(vegn_cohort_type), intent(in) :: c1,c2

   real, parameter :: mindensity = 1.0E-4
   logical :: sameSpecies, sameLayer, sameSize, sameSizeTree, sameSizeGrass, lowDensity

   sameSpecies  = c1%species == c2%species
   sameLayer    = (c1%layer == c2%layer) ! .and. (c1%firstlayer == c2%firstlayer)
   sameSizeTree = (spdata(c1%species)%lifeform > 0).and.  &
                  (spdata(c2%species)%lifeform > 0).and.  &
                 ((abs(c1%DBH - c2%DBH)/c2%DBH < 0.2 ) .or.  &
                  (abs(c1%DBH - c2%DBH)        < 0.001))  ! it'll be always true for grasses
   sameSizeGrass= (spdata(c1%species)%lifeform ==0) .and. &
                  (spdata(c2%species)%lifeform ==0) .and. &
                 (((c1%DBH == c2%DBH) .And.(c1%nsc == c2%nsc)) .OR. &
                    c1%topyear==c2%topyear)  ! it'll be always true for grasses
   sameSize = sameSizeTree .OR. sameSizeGrass
   lowDensity  = .FALSE. ! c1%nindivs < mindensity 
                         ! Weng, 2014-01-27, turned off
   cohorts_can_be_merged = sameSpecies .and. sameLayer .and. sameSize
end function 

!=================================================================================
! given an array of cohorts, create a new array with old cohorts re-arranged
! in layers according to their height and crown areas.
subroutine relayer_cohorts (vegn)
  type(vegn_tile_type), intent(inout) :: vegn ! input cohorts

  ! ---- local constants
  real, parameter :: tolerance = 1e-6 
  real, parameter :: layer_vegn_cover = 1.0   
  ! ---- local vars
  integer :: idx(vegn%n_cohorts) ! indices of cohorts in decreasing height order
  integer :: i ! new cohort index
  integer :: k ! old cohort index
  integer :: L ! layer index (top-down)
  integer :: N0,N1 ! initial and final number of cohorts 
  real    :: frac ! fraction of the layer covered so far by the canopies
  type(vegn_cohort_type), pointer :: cc(:),new(:)
  real    :: nindivs
  logical :: rand_sorting = .False.

!  rand_sorting = .TRUE. ! .False.
  
  ! rank cohorts in descending order by height. For now, assume that they are 
  ! in order
  if(rand_sorting)then
      call random_sorting(vegn)
      do k=1,vegn%n_cohorts
          idx(k)=k
      enddo
      N0 = vegn%n_cohorts; cc=>vegn%cohorts
  else
      N0 = vegn%n_cohorts; cc=>vegn%cohorts
      call rank_descending(cc(1:N0)%height,idx)
  endif
  
  ! calculate max possible number of new cohorts : it is equal to the number of
  ! old cohorts, plus the number of layers -- since the number of full layers is 
  ! equal to the maximum number of times an input cohort can be split by a layer 
  ! boundary.
  N1 = vegn%n_cohorts + int(sum(cc(1:N0)%nindivs*cc(1:N0)%crownarea))
  allocate(new(N1))

  ! copy cohort information to the new cohorts, splitting the old cohorts that 
  ! stride the layer boundaries
  i = 1 ; k = 1 ; L = 1 ; frac = 0.0 ; nindivs = cc(idx(k))%nindivs
  do 
     new(i)         = cc(idx(k))
     new(i)%nindivs = min(nindivs,(layer_vegn_cover-frac)/cc(idx(k))%crownarea)
     new(i)%layer   = L
     if (L==1) new(i)%firstlayer = 1
!    if (L>1)  new(i)%firstlayer = 0  ! switch off "push-down effects"
     frac = frac+new(i)%nindivs*new(i)%crownarea
     nindivs = nindivs - new(i)%nindivs
     
     if (abs(nindivs*cc(idx(k))%crownarea)<tolerance) then
       new(i)%nindivs = new(i)%nindivs + nindivs ! allocate the remainder of individuals to the last cohort
       if (k==N0) exit ! end of loop
       k = k+1 ; nindivs = cc(idx(k))%nindivs  ! go to the next input cohort
     endif
     
     if (abs(layer_vegn_cover-frac)<tolerance) then
       L = L+1 ; frac = 0.0              ! start new layer
     endif

     i = i+1
  enddo
  
  ! replace the array of cohorts
  deallocate(vegn%cohorts)
  vegn%cohorts => new ; vegn%n_cohorts = i

! Merge those pushed into the second layer
!  call merge_noncanopy_cohorts(vegn)

end subroutine relayer_cohorts


! =============================================================================
subroutine random_sorting (vegn)
  type(vegn_tile_type), intent(inout) :: vegn

  ! ---- local vars
  type(vegn_cohort_type), pointer :: cc(:)   ! current cohort
  type(vegn_cohort_type)         :: ckey
  integer :: i,j
  real,dimension(1):: harvest 

  call random_seed

!  write(*,*)"TEST, Layering:"
! sort cohorts
  associate(cc=>vegn%cohorts)

  do i = 1, vegn%n_cohorts -1
     do j=i+1, vegn%n_cohorts
!       sort first layer cohorts accordint to their heights
        if(cc(j)%height > cc(i)%height .AND. cc(j)%layer == 1)then
           ckey  = cc(i)
           cc(i) = cc(j)
           cc(j) = ckey
        endif
!       sort the 2nd layer cohorts randomly
        call random_number(harvest)
        if(cc(j)%layer /= 1 .AND.cc(i)%layer /= 1 .AND. harvest(1) >0.5)then
           ckey  = cc(i)
           cc(i) = cc(j)
           cc(j) = ckey
        endif
     enddo    
  enddo

  end associate

!! count canopy layers
!     vegn%n_canopycc = 0
!     do i = 1,vegn%n_cohorts
!        if (vegn%cohorts(i)%layer ==1) then
!           vegn%n_canopycc=vegn%n_canopycc+1
!        else
!           exit
!        endif
!     enddo

end subroutine random_sorting


! ============================================================================
! ============================================================================
!!!! LM3 subrountines

! =============================================================================
subroutine vegn_phenology_lm3(vegn, wilt)
  type(vegn_tile_type), intent(inout) :: vegn
  real, intent(in) :: wilt ! ratio of wilting to saturated water content

  ! ---- local vars
  type(vegn_cohort_type), pointer :: cc
  real :: leaf_litter,root_litter;    
  real :: theta_crit; ! critical ratio of average soil water to sat. water
  integer :: i
  
  vegn%litter = 0

  do i = 1,vegn%n_cohorts   
     cc => vegn%cohorts(i)
     
     if(is_watch_point())then
        write(*,*)'####### vegn_phenology_lm3 #######'
        __DEBUG4__(vegn%theta_av, wilt, spdata(cc%species)%cnst_crit_phen, spdata(cc%species)%fact_crit_phen)
        __DEBUG1__(cc%species)
        __DEBUG2__(vegn%tc_av,spdata(cc%species)%tc_crit)
     endif
     ! if drought-deciduous or cold-deciduous species
     ! temp=10 degrees C gives good growing season pattern        
     ! spp=0 is c3 grass,1 c3 grass,2 deciduous, 3 evergreen
     ! assumption is that no difference between drought and cold deciduous
     cc%status = LEAF_ON; ! set status to indicate no leaf drop
      
     if(cc%species < 4 )then! deciduous species
        ! actually either fact_crit_phen or cnst_crit_phen is zero, enforced
        ! by logic in the vegn_data.F90
        theta_crit = spdata(cc%species)%cnst_crit_phen &
              + wilt*spdata(cc%species)%fact_crit_phen
        theta_crit = max(0.0,min(1.0, theta_crit))
        if ( (vegn%theta_av < theta_crit) &
             .or.(vegn%tc_av < spdata(cc%species)%tc_crit) ) then
           cc%status = LEAF_OFF; ! set status to indicate leaf drop 
           cc%leaf_age = 0;
           
           leaf_litter = (1.0-l_fract)*cc%bl*cc%nindivs;
           root_litter = (1.0-l_fract)*cc%br*cc%nindivs;
           
           ! add to patch litter flux terms
           vegn%litter = vegn%litter + leaf_litter + root_litter;
           
           vegn%fast_soil_C = vegn%fast_soil_C +    fsc_liv *(leaf_litter+root_litter);
           vegn%slow_soil_C = vegn%slow_soil_C + (1-fsc_liv)*(leaf_litter+root_litter);
     
           ! vegn%fsc_in+=data->fsc_liv*(leaf_litter+root_litter);
           ! vegn%ssc_in+=(1.0-data->fsc_liv)*(leaf_litter+root_litter);
           vegn%fsc_in  = vegn%fsc_in  + leaf_litter+root_litter;
           vegn%veg_out = vegn%veg_out + leaf_litter+root_litter;
           
           cc%blv = cc%blv + l_fract*(cc%bl+cc%br);
           cc%bl  = 0.0;
           cc%br  = 0.0;
           cc%lai = 0.0;
           
           ! update state
           cc%bliving = cc%blv + cc%br + cc%bl + cc%bsw;
           call update_bio_living_fraction(cc);   
        endif
     endif
  enddo
end subroutine vegn_phenology_lm3


! ============================================================================
subroutine vegn_carbon_int_lm3(vegn, soilt, theta, diag)
  type(vegn_tile_type), intent(inout) :: vegn
  real, intent(in) :: soilt ! average temperature of soil for soil carbon decomposition, deg K
  real, intent(in) :: theta ! average soil wetness, unitless
  type(diag_buff_type), intent(inout) :: diag

  type(vegn_cohort_type), pointer :: cc
  type(vegn_cohort_type), pointer :: c(:) ! for debug only
  real :: resp, resl, resr, resg ! respiration terms accumulated for all cohorts 
  real :: md_alive, md_wood;
  real :: md ! plant tissue maintenance, kg C/timestep
  real :: gpp ! gross primary productivity per tile
  integer :: sp ! shorthand for current cohort specie
  integer :: i

  c=>vegn%cohorts(1:vegn%n_cohorts)
  if(is_watch_point()) then
     write(*,*)'#### vegn_carbon_int_lm3 input ####'
     __DEBUG2__(soilt,theta)
     __DEBUG1__(c%npp_previous_day)
     __DEBUG1__(c%carbon_gain)
     __DEBUG1__(c%carbon_loss)
     __DEBUG1__(c%bwood_gain)
  endif

  !  update plant carbon
  vegn%npp = 0
  resp = 0 ; resl = 0 ; resr = 0 ; resg = 0 ; gpp = 0
  do i = 1, vegn%n_cohorts   
     cc => vegn%cohorts(i)
     sp = cc%species

     call plant_respiration(cc,soilt);
   
     cc%gpp = (cc%An_op - cc%An_cl)*mol_C*cc%leafarea;
     cc%npp = cc%gpp - cc%resp;
   
     if(cc%npp_previous_day > 0) then
        cc%resg = GROWTH_RESP*cc%npp_previous_day;
        cc%npp  = cc%npp  - GROWTH_RESP*cc%npp_previous_day;
        cc%resp = cc%resp + GROWTH_RESP*cc%npp_previous_day;
     else
        cc%resg = 0;
     endif
     ! accumulate npp for the current day
     cc%npp_previous_day_tmp = cc%npp_previous_day_tmp + cc%npp; 
     
     cc%carbon_gain = cc%carbon_gain + cc%npp*dt_fast_yr;
     
     ! check if leaves/roots are present and need to be accounted in maintenance
     if(cc%status == LEAF_ON) then
        md_alive = (cc%Pl * spdata(sp)%alpha(CMPT_LEAF) + &
                    cc%Pr * spdata(sp)%alpha(CMPT_ROOT))* &
              cc%bliving*dt_fast_yr;    
     else
        md_alive = 0
     endif
     
     ! compute branch and coarse wood losses for tree types
     md_wood =0;
     if (sp > 1) then
        md_wood = 0.6 *cc%bwood * spdata(sp)%alpha(CMPT_WOOD)*dt_fast_yr;
     endif
        
     md = md_alive + cc%Psw_alphasw * cc%bliving * dt_fast_yr;
     cc%bwood_gain = cc%bwood_gain + cc%Psw_alphasw * cc%bliving * dt_fast_yr;
     cc%bwood_gain = cc%bwood_gain - md_wood;
!     if (cc%bwood_gain < 0.0) cc%bwood_gain=0.0; ! potential non-conservation ?
     cc%carbon_gain = cc%carbon_gain - md;
     cc%carbon_loss = cc%carbon_loss + md; ! used in diagnostics only

     ! add maintenance demand from leaf and root pools to fast soil carbon
     vegn%fast_soil_C = vegn%fast_soil_C + (   fsc_liv *md_alive +    fsc_wood *md_wood)*cc%nindivs;
     vegn%slow_soil_C = vegn%slow_soil_C + ((1-fsc_liv)*md_alive + (1-fsc_wood)*md_wood)*cc%nindivs;

     ! for budget tracking
!/*     cp->fsc_in+= data->fsc_liv*md_alive+data->fsc_wood*md_wood; */
!/*     cp->ssc_in+= (1.- data->fsc_liv)*md_alive+(1-data->fsc_wood)*md_wood; */
     vegn%fsc_in  = vegn%fsc_in + (     1 *md_alive+   0 *md_wood)*cc%nindivs;
     vegn%ssc_in  = vegn%ssc_in + ((1.- 1)*md_alive+(1-0)*md_wood)*cc%nindivs;

     vegn%veg_in  = vegn%veg_in  + cc%npp*cc%nindivs*dt_fast_yr;
     vegn%veg_out = vegn%veg_out + (md_alive+md_wood)*cc%nindivs;

     ! accumulate tile-level NPP and GPP
     vegn%npp = vegn%npp + cc%npp*cc%nindivs
     gpp = gpp + cc%gpp*cc%nindivs
     ! accumulate respiration terms for tile-level reporting
     resp = resp + cc%resp*cc%nindivs ; resl = resl + cc%resl*cc%nindivs
     resr = resr + cc%resr*cc%nindivs ; resg = resg + cc%resg*cc%nindivs
     ! update cohort age
     cc%age = cc%age + dt_fast_yr
  enddo
  if(is_watch_point()) then
     write(*,*)'#### vegn_carbon_int_lm3 output ####'
     __DEBUG1__(c%species)
     __DEBUG1__(c%bl)
     __DEBUG1__(c%br)
     __DEBUG1__(c%bsw)
     __DEBUG1__(c%bwood)
     __DEBUG1__(c%An_op)
     __DEBUG1__(c%An_cl)
     __DEBUG1__(c%lai)
     __DEBUG1__(c%npp)
     __DEBUG1__(c%gpp)
     __DEBUG1__(c%resp)
     __DEBUG1__(c%resl)
     __DEBUG1__(c%resr)
     __DEBUG1__(c%resg)
     __DEBUG1__(c%carbon_gain)
     __DEBUG1__(c%carbon_loss)
     __DEBUG1__(c%bwood_gain)
  endif

  ! update soil carbon
  call Dsdt(vegn, diag, soilt, theta)

  ! NEP is equal to NNP minus soil respiration
  vegn%nep = vegn%npp - vegn%rh

  call update_soil_pools(vegn)
  vegn%age = vegn%age + dt_fast_yr;


  ! ---- diagnostic section
  call send_tile_data(id_gpp,gpp,diag)
  call send_tile_data(id_npp,vegn%npp,diag)
  call send_tile_data(id_nep,vegn%nep,diag)
  call send_tile_data(id_litter,vegn%litter,diag)
  call send_tile_data(id_resp, resp, diag)
  call send_tile_data(id_resl, resl, diag)
  call send_tile_data(id_resr, resr, diag)
  call send_tile_data(id_resg, resg, diag)
  call send_tile_data(id_soilt,soilt,diag)
  call send_tile_data(id_theta,theta,diag)
  
end subroutine vegn_carbon_int_lm3

! ============================================================================
! Nitrogen mineralization and immoblization
subroutine Dsdt(vegn, diag, soilt, theta)
  type(vegn_tile_type), intent(inout) :: vegn
  type(diag_buff_type), intent(inout) :: diag
  real                , intent(in)    :: soilt ! soil temperature, deg K 
  real                , intent(in)    :: theta

  real :: fast_C_loss,slow_C_loss
  real :: runoff ! Runoff in one time step (kg water)
  real :: N_loss
  real :: fast_N_free = 0.0
  real :: slow_N_free = 0.0
  real :: fCNfast,fCNslow, CNfast, CNslow
  real :: A  ! decomp rate reduction due to moisture and temperature
  
  A=A_function(soilt,theta)
  runoff = vegn%Wrunoff * 365*24*3600 *dt_fast_yr ! kgH2O m-2 s-1 -->kgH2O m-2/time step
  fast_C_loss = vegn%fast_soil_C*A*K1*dt_fast_yr
  slow_C_loss = vegn%slow_soil_C*A*K2*dt_fast_yr
! N_loss      = MAX(0.,vegn%QNmineral-0.0006) * A * K_nitrogen * dt_fast_yr
! N_loss      = MAX(0.,vegn%QNmineral)        * A * K_nitrogen * dt_fast_yr
  N_loss      = MAX(0.,vegn%QNmineral)*(1.0-exp(-0.05*runoff)+ A * K_nitrogen) * dt_fast_yr

! Assume half Nitrogen goes out to mineral N pool
  CNfast = vegn%fast_soil_C / vegn%QNfastSOM
  CNslow = vegn%slow_soil_C / vegn%QNslowSOM
  fCNfast = MAX(1.0,(CN0fastSOM/CNfast)**2)
  fCNslow = MAX(1.0,(CN0slowSOM/CNslow)**2)
  fast_N_free = fCNfast * fast_C_loss/CNfast
  slow_N_free = fCNslow * slow_C_loss/CNslow
  vegn%QNmineral = vegn%QNmineral + fast_N_free + slow_N_free 
  vegn%annualN   = vegn%annualN + fast_N_free + slow_N_free
! update soil N pool
  vegn%QNfastSOM = vegn%QNfastSOM - fast_N_free
  vegn%QNslowSOM = vegn%QNslowSOM - slow_N_free

! update soil carbon pool  
  vegn%fast_soil_C = vegn%fast_soil_C - fast_C_loss
  vegn%slow_soil_C = vegn%slow_soil_C - slow_C_loss

  ! for budget check
  vegn%fsc_out = vegn%fsc_out + fast_C_loss
  vegn%ssc_out = vegn%ssc_out + slow_C_loss

  ! loss of C to atmosphere and leaching
  vegn%rh =   (fast_C_loss+slow_C_loss)/dt_fast_yr

  ! accumulate decomposition rate reduction for the soil carbon restart output
  vegn%asoil_in = vegn%asoil_in + A

!! Nitrogen mineralization, immoblization, and deposition
! Check if soil C/N is above CN0
  fast_N_free = MAX(0.0,vegn%QNfastSOM - vegn%fast_soil_C/CN0fastSOM)
  slow_N_free = MAX(0.0,vegn%QNslowSOM - vegn%slow_soil_C/CN0slowSOM)
  vegn%QNfastSOM = vegn%QNfastSOM - fast_N_free
  vegn%QNslowSOM = vegn%QNslowSOM - slow_N_free ! + vegn%N_input * dt_fast_yr !
!  Weng 2015-09-17: put N_input to QNslowSOM, not the mineral N pool
  vegn%QNmineral = vegn%QNmineral + fast_N_free + slow_N_free &
                 + vegn%N_input * dt_fast_yr                  &
                 - N_loss
  vegn%annualN   = vegn%annualN + fast_N_free + slow_N_free   &
                 + vegn%N_input * dt_fast_yr                  &
                 - N_loss
  ! ---- diagnostic section
!  call send_tile_data(id_fsc, vegn%fast_soil_C, diag) ! in vegetation.F90
!  call send_tile_data(id_ssc, vegn%slow_soil_C, diag) ! in vegetation.F90
  call send_tile_data(id_N_loss,   N_loss/dt_fast_yr, diag)
  call send_tile_data(id_QNmineral,vegn%QNmineral,    diag)
  call send_tile_data(id_rsoil_fast, fast_C_loss/dt_fast_yr, diag)
  call send_tile_data(id_rsoil, vegn%rh, diag)
  call send_tile_data(id_asoil, A, diag)

end subroutine Dsdt

! ============================================================================
! Nitrogen mineralization and immoblization with microbial C & N pools
! it's a new decomposition model with coupled C & N pools and variable 
! carbon use efficiency
subroutine Dsdtmicrobes(vegn, diag, soilt, theta)
  type(vegn_tile_type), intent(inout) :: vegn
  type(diag_buff_type), intent(inout) :: diag
  real                , intent(in)    :: soilt ! soil temperature, deg K 
  real                , intent(in)    :: theta

  real :: CUE0=0.4  ! default microbial CUE
  real :: phoMicrobial = 2.5 ! turnover rate of microbes (yr-1)
  real :: CUEfast,CUEslow
  real :: CNm = 10.0  ! Microbial C/N ratio
  real :: NforM, fNM=0.0  ! mineral N available for microbes
  real :: micr_C_loss, fast_C_loss, slow_C_loss
  real :: runoff ! kg m-2 /step
  real :: etaN = 0.05! loss rate of Nmineral with runoff
  real :: N_loss
  real :: DON_fast,DON_slow,DON_loss ! Dissolved organic N loss, kg N m-2 step-1
  real :: fDON=0.0   ! 0.02     ! fractio of DON production in decomposition
  real :: fast_N_free 
  real :: slow_N_free 
  real :: CNfast, CNslow
  real :: A  ! decomp rate reduction due to moisture and temperature
  
  runoff = vegn%Wrunoff * 365*24*3600 *dt_fast_yr ! kgH2O m-2 s-1 -->kgH2O m-2/time step
! CN ratios of soil C pools
  CNfast = vegn%fast_soil_C/vegn%QNfastSOM
  CNslow = vegn%slow_soil_C/vegn%QNslowSOM

! C decomposition
  A=A_function(soilt,theta)
  micr_C_loss = vegn%microbialC *A*phoMicrobial* dt_fast_yr
  fast_C_loss = vegn%fast_soil_C*A*K1          * dt_fast_yr
  slow_C_loss = vegn%slow_soil_C*A*K2          * dt_fast_yr

! Carbon use efficiencies of microbes
  NforM = fNM * vegn%QNmineral
  CUEfast = MIN(CUE0,CNm*(fast_C_loss/CNfast + NforM)/fast_C_loss)
  CUEslow = MIN(CUE0,CNm*(slow_C_loss/CNslow + NforM)/slow_C_loss)

! update C and N pools
! Carbon pools
  vegn%microbialC  = vegn%microbialC - micr_C_loss &
                    + fast_C_loss * CUEfast &
                    + slow_C_loss * CUEslow
  vegn%fast_soil_C = vegn%fast_soil_C - fast_C_loss
  vegn%slow_soil_C = vegn%slow_soil_C - slow_C_loss

! DON loss, revised 2016-03-03
  fDON=0.0
  DON_fast    = fDON*fast_C_loss/CNfast*(1.0-exp(-etaN*runoff))
  DON_slow    = fDON*slow_C_loss/CNslow*(1.0-exp(-etaN*runoff))
  DON_loss    = DON_fast + DON_slow

! Nitrogen pools
  vegn%microbialN= vegn%microbialC/CNm
  vegn%QNfastSOM = vegn%QNfastSOM - fast_C_loss/CNfast - DON_fast
  vegn%QNslowSOM = vegn%QNslowSOM - slow_C_loss/CNslow - DON_slow

! Mixing of microbes to litters
  vegn%fast_soil_C = vegn%fast_soil_C + MLmixRatio*fast_C_loss * CUEfast
  vegn%QNfastSOM   = vegn%QNfastSOM   + MLmixRatio*fast_C_loss * CUEfast/CNm

  vegn%slow_soil_C = vegn%slow_soil_C + MLmixRatio*slow_C_loss * CUEslow
  vegn%QNslowSOM   = vegn%QNslowSOM   + MLmixRatio*slow_C_loss * CUEslow/CNm

  vegn%microbialC  = vegn%microbialC  - MLmixRatio*(fast_C_loss*CUEfast+slow_C_loss*CUEslow)
  vegn%microbialN  = vegn%microbialC/CNm
    
! update mineral N pool (QNmineral)
  fast_N_free = MAX(0.0, fast_C_loss*(1./CNfast - CUEfast/CNm))
  slow_N_free = MAX(0.0, slow_C_loss*(1./CNslow - CUEslow/CNm))

! N_loss      = MAX(0.,vegn%QNmineral)        * A * K_nitrogen * dt_fast_yr
  N_loss      = MAX(0.,vegn%QNmineral)* &
                (1.0-exp(-etaN*runoff)+ A * K_nitrogen) *dt_fast_yr
 ! put N_input to QNslowSOM, Weng 2015-09-17
  vegn%QNmineral = vegn%QNmineral - N_loss     &
                  + vegn%N_input * dt_fast_yr  &
                  + fast_N_free + slow_N_free  &
                  + micr_C_loss/CNm
  vegn%annualN   = vegn%annualN - N_loss     &
                  + vegn%N_input * dt_fast_yr  &
                  + fast_N_free + slow_N_free  &
                  + micr_C_loss/CNm
 ! Weng 2015-09-17
!  vegn%QNslowSOM = vegn%QNslowSOM + vegn%N_input * dt_fast_yr 

! Check if soil C/N is above CN0
  fast_N_free = MAX(0. ,vegn%QNfastSOM - vegn%fast_soil_C/CN0fastSOM)
  slow_N_free = MAX(0. ,vegn%QNslowSOM - vegn%slow_soil_C/CN0slowSOM)
  vegn%QNfastSOM = vegn%QNfastSOM - fast_N_free
  vegn%QNslowSOM = vegn%QNslowSOM - slow_N_free
  vegn%QNmineral = vegn%QNmineral + fast_N_free + slow_N_free
  vegn%annualN   = vegn%annualN + fast_N_free + slow_N_free

! for budget check
  vegn%fsc_out = vegn%fsc_out + fast_C_loss
  vegn%ssc_out = vegn%ssc_out + slow_C_loss

! loss of C to atmosphere and leaching
  vegn%rh =   (micr_C_loss + &
               fast_C_loss*(1.-CUEfast)+ &
               slow_C_loss*(1.-CUEslow)) &
              /dt_fast_yr

! accumulate decomposition rate reduction for the soil carbon restart output
  vegn%asoil_in = vegn%asoil_in + A

  ! ---- diagnostic section
!  call send_tile_data(id_fsc, vegn%fast_soil_C, diag) ! in vegetation.F90
!  call send_tile_data(id_ssc, vegn%slow_soil_C, diag) ! in vegetation.F90
  call send_tile_data(id_N_loss,   N_loss/dt_fast_yr, diag)
  call send_tile_data(id_QNmineral,vegn%QNmineral,    diag)
  call send_tile_data(id_rsoil_fast, fast_C_loss/dt_fast_yr, diag)
  call send_tile_data(id_rsoil, vegn%rh, diag)
  call send_tile_data(id_asoil, A, diag)

end subroutine Dsdtmicrobes

! ============================================================================
! The combined reduction in decomposition rate as a funciton of TEMP and MOIST
! Based on CENTURY Parton et al 1993 GBC 7(4):785-809 and Bolker's copy of
! CENTURY code
function A_function(soilt, theta) result(A)
  real :: A                 ! return value, resulting reduction in decomposition rate
  real, intent(in) :: soilt ! effective temperature for soil carbon decomposition
  real, intent(in) :: theta 

  real :: soil_temp; ! temperature of the soil, deg C
  real :: Td; ! rate multiplier due to temp
  real :: Wd; ! rate reduction due to mositure

  ! coefficeints and terms used in temperaturex term
  real :: Topt,Tmax,t1,t2,tshl,tshr;

  soil_temp = soilt-273.16;

  ! EFFECT OF TEMPERATURE
  ! from Bolker's century code
  Tmax=45.0;
  if (soil_temp > Tmax) soil_temp = Tmax;
  Topt=35.0;
  tshr=0.2; tshl=2.63;
  t1=(Tmax-soil_temp)/(Tmax-Topt);
  t2=exp((tshr/tshl)*(1.-t1**tshl));
  Td=t1**tshr*t2;

  if (soil_temp > -10) Td=Td+0.05;
  if (Td > 1.) Td=1.;

  ! EFFECT OF MOISTURE
  ! Linn and Doran, 1984, Soil Sci. Amer. J. 48:1267-1272
  ! This differs from the Century Wd
  ! was modified by slm/ens based on the figures from the above paper 
  !     (not the reported function)

  if(theta <= 0.3) then
     Wd = 0.2;
  else if(theta <= 0.6) then
     Wd = 0.2+0.8*(theta-0.3)/0.3;
  else 
     Wd = exp(2.3*(0.6-theta));
  endif

  A = (Td*Wd); ! the combined (multiplicative) effect of temp and water
               ! on decomposition rates
end function A_function


! ============================================================================
! calculated thermal inhibition factor depending on temperature
function thermal_inhibition(T) result(tfs); real tfs
  real, intent(in) :: T ! demperature, degK
  
  tfs = exp(3000.0*(1.0/288.16-1.0/T));
  tfs = tfs / ( &
              (1.0+exp(0.4*(5.0-T+273.16)))* &
              (1.0+exp(0.4*(T - 273.16-45.0)))&
              )
end function

! need a new equation for respiration, Weng 2014-01-13
function Rm_T_response_function(T) result(tfs); real tfs
  real, intent(in) :: T ! demperature, degK
!  According to Jarvi & Burton 2013, Tree Physiology, sugar maple root resp
!  R0=2.87 (kg C of resp/ Kg C of root biomass)
!  E0= 9000
!  R(T)=R0*exp[E0*(1/288.16-1/T)

  tfs = exp(9000.0*(1.0/288.16-1.0/T));

end function


! =============================================================================
subroutine vegn_biogeography(vegn)
  type(vegn_tile_type), intent(inout) :: vegn

  ! ---- local vars
  integer :: i

  do i = 1, vegn%n_cohorts   
     call update_species(vegn%cohorts(i), vegn%t_ann, vegn%t_cold, &
          vegn%p_ann*seconds_per_year, vegn%ncm, vegn%landuse)
  enddo
end subroutine

! =============================================================================
! The stuff below comes from she_update.c -- it looks like it belongs here, 
! since it is essentially a part of the carbon integration (update_patch_fast
! is only called immediately after carbon_int in lm3v)
! =============================================================================


! =============================================================================
subroutine update_soil_pools(vegn)
  type(vegn_tile_type), intent(inout) :: vegn
  
  ! ---- local vars
  real :: delta;

  ! update fsc input rate so that intermediate fsc pool is never
  ! depleted below zero; on the other hand the pool can be only 
  ! depleted, never increased
  vegn%fsc_rate = MAX( 0.0, MIN(vegn%fsc_rate, vegn%fsc_pool/dt_fast_yr));
  delta = vegn%fsc_rate*dt_fast_yr;
  vegn%fast_soil_C = vegn%fast_soil_C + delta;
  vegn%fsc_pool    = vegn%fsc_pool    - delta;

  ! update ssc input rate so that intermediate ssc pool is never
  ! depleted below zero; on the other hand the pool can be only 
  ! depleted, never increased
  vegn%ssc_rate = MAX(0.0, MIN(vegn%ssc_rate, vegn%ssc_pool/dt_fast_yr));
  delta = vegn%ssc_rate*dt_fast_yr;
  vegn%slow_soil_C = vegn%slow_soil_C + delta;
  vegn%ssc_pool    = vegn%ssc_pool    - delta;
end subroutine update_soil_pools


end module vegn_dynamics_mod
