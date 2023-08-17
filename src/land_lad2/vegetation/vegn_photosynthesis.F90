module vegn_photosynthesis_mod

#include "../shared/debug.inc"

use utilities_mod,      only : logunit, error_mesg, FATAL
use constants_mod,      only : TFREEZE 
use sphum_mod,          only : qscomp

use land_constants_mod, only : Rugas, seconds_per_year, mol_h2o, mol_air
use land_debug_mod,     only : is_watch_point
use vegn_data_mod,      only : PT_C4, spdata
use vegn_cohort_mod,    only : vegn_cohort_type, get_vegn_wet_frac

implicit none
private

! ==== public interfaces =====================================================
public :: vegn_photosynthesis_init
public :: vegn_photosynthesis
! ==== end of public interfaces ==============================================

! ==== module constants ======================================================
character(len=*), private, parameter :: &
   version = '$Id: vegn_photosynthesis.F90,v 1.1.2.3 2012/06/19 18:34:55 pjp Exp $', &
   tagname = '$Name: no_fms_b_pjp $', &
   module_name = 'vegn_photosynthesis'
! values for internal vegetation photosynthesis option selector
integer, parameter :: VEGN_PHOT_SIMPLE  = 1 ! zero photosynthesis
integer, parameter :: VEGN_PHOT_LEUNING = 2 ! photosynthesis according to simplified Leuning model

! ==== module variables ======================================================
integer :: vegn_phot_option = -1 ! selector of the photosynthesis option


contains


! ============================================================================
subroutine vegn_photosynthesis_init(photosynthesis_to_use)
  character(*), intent(in) :: photosynthesis_to_use

  write (logunit,'(/,80("="),/(a))') trim(version), trim(tagname)

  ! convert symbolic names of photosynthesis options into numeric IDs to
  ! speed up selection during run-time
  if (trim(photosynthesis_to_use)=='simple') then
     vegn_phot_option = VEGN_PHOT_SIMPLE
  else if (trim(photosynthesis_to_use)=='leuning') then
     vegn_phot_option = VEGN_PHOT_LEUNING
  else
     call error_mesg('vegn_photosynthesis_init',&
          'vegetation photosynthesis option photosynthesis_to_use="'//&
          trim(photosynthesis_to_use)//'" is invalid, use "simple" or "leuning"',&
          FATAL)
  endif

end subroutine vegn_photosynthesis_init


! ============================================================================
! compute stomatal conductance, photosynthesis and respiration
! updates cohort%An_op and cohort%An_cl
subroutine vegn_photosynthesis ( cohort, &
     PAR_dn, PAR_net, cana_q, cana_co2, p_surf, drag_q, &
     soil_beta, soil_water_supply,  stomatal_cond )
  type(vegn_cohort_type), intent(inout) :: cohort
  real, intent(in)  :: PAR_dn   ! downward PAR at the top of the canopy, W/m2 
  real, intent(in)  :: PAR_net  ! net PAR absorbed by the canopy, W/m2
  real, intent(in)  :: cana_q   ! specific humidity in canopy air space, kg/kg
  real, intent(in)  :: cana_co2 ! co2 concentration in canopy air space, mol CO2/mol dry air
  real, intent(in)  :: p_surf   ! surface pressure
  real, intent(in)  :: drag_q   ! drag coefficient for specific humidity
  real, intent(in)  :: soil_beta
  real, intent(in)  :: soil_water_supply ! max supply of water to roots per unit
                                ! active root biomass per second, kg/(indiv s)
  real, intent(out) :: stomatal_cond ! stomatal conductance, m/s(?)

  ! ---- local vars
  real    :: water_supply ! water supply per m2 of leaves
  real    :: fw, fs ! wet and snow-covered fraction of leaves
  real    :: psyn   ! net photosynthesis, mol C/(m2 of leaves s)
  real    :: resp   ! leaf respiration, mol C/(m2 of leaves s)
  real    :: w_scale2

  select case (vegn_phot_option)

  case(VEGN_PHOT_SIMPLE)
     ! beta non-unity only for "beta" models
     stomatal_cond = soil_beta / (cohort%rs_min  + (1-soil_beta)/drag_q)
     cohort%An_op  = 0
     cohort%An_cl  = 0
     psyn = 0
     resp = 0

  case(VEGN_PHOT_LEUNING)
     if(cohort%lai > 0) then
        ! recalculate the water supply to mol H20 per m2 of leaf per second
        water_supply = soil_water_supply/(mol_h2o*cohort%leafarea)
      
        call get_vegn_wet_frac (cohort, fw=fw, fs=fs)
        call gs_Leuning(PAR_dn, PAR_net, cohort%prog%Tv, cana_q, cohort%lai, &
             cohort%leaf_age, p_surf, water_supply, cohort%species, cohort%pt, &
             cana_co2, cohort%extinct, fs+fw, cohort%layer, &
             ! output:
             stomatal_cond, psyn, resp,w_scale2 )
        ! store the calculated photosynthesis and photorespiration for future use
        ! in carbon_int
        cohort%An_op  = psyn * seconds_per_year
        cohort%An_cl  = resp * seconds_per_year
        cohort%w_scale  = w_scale2
        ! convert stomatal conductance from units per unit area of leaf to the 
        ! units per unit area of cohort
        stomatal_cond = stomatal_cond*cohort%lai
     else
        ! no leaves means no photosynthesis and no stomatal conductance either
        cohort%An_op  = 0
        cohort%An_cl  = 0
        cohort%w_scale  = -9999
        stomatal_cond = 0
     endif

  case default
     call error_mesg('vegn_stomatal_cond', &
          'invalid vegetation photosynthesis option', FATAL)
  end select

end subroutine vegn_photosynthesis


! ============================================================================
subroutine gs_Leuning(rad_top, rad_net, tl, ea, lai, leaf_age, &
                   p_surf, ws, pft, pt, ca, kappa, leaf_wet, layer, &
                   gs, apot, acl,w_scale2)
  real,    intent(in)    :: rad_top ! PAR dn on top of the canopy, w/m2
  real,    intent(in)    :: rad_net ! PAR net on top of the canopy, w/m2
  real,    intent(in)    :: tl   ! leaf temperature, degK
  real,    intent(in)    :: ea   ! specific humidity in the canopy air (?), kg/kg
  real,    intent(in)    :: lai  ! leaf area index
  real,    intent(in)    :: leaf_age ! age of leaf since budburst (deciduos), days
  real,    intent(in)    :: p_surf ! surface pressure, Pa
  real,    intent(in)    :: ws   ! water supply, mol H20/(m2 of leaf s)
  integer, intent(in)    :: pft  ! species
  integer, intent(in)    :: pt   ! physiology type (C3 or C4)
  real,    intent(in)    :: ca   ! concentartion of CO2 in the canopy air space, mol CO2/mol dry air
  real,    intent(in)    :: kappa! canopy extinction coefficient (move inside f(pft))
  real,    intent(in)    :: leaf_wet ! fraction of leaf that's wet or snow-covered
  integer, intent(in)    :: layer  ! the layer of this canopy
  ! note that the output is per area of leaf; to get the quantities per area of
  ! land, multiply them by LAI
  real,    intent(out)   :: gs   ! stomatal conductance, m/s
  real,    intent(out)   :: apot ! net photosynthesis, mol C/(m2 s)
  real,    intent(out)   :: acl  ! leaf respiration, mol C/(m2 s)
  real,    intent(out)   :: w_scale2  ! water limitation 

  ! ---- local vars     
  ! photosynthesis
  real :: vm;
  real :: kc,ko; ! Michaelis-Menten constants for CO2 and O2, respectively
  real :: ci;
  real :: capgam; ! CO2 compensation point
  real :: f2,f3;
  real :: coef0,coef1;

  real :: Resp;

  ! conductance related
  real :: b;
  real :: ds;  ! humidity deficit, kg/kg
  real :: hl;  ! saturated specific humidity at the leaf temperature, kg/kg
  real :: do1;
  
  ! misceleneous
  real :: dum2;
  real, parameter :: light_crit = 0;
  real, parameter :: gs_lim = 0.25;
  real, parameter :: Rgas = 8.314 ! J mol-1 K-1, universal gas constant
  ! new average computations
  real :: lai_eq;
  real, parameter :: rad_phot = 0.0000046 ! PAR conversion factor of J -> mol of quanta 
  real :: light_top;
  real :: par_net;
  real :: Ag;
  real :: An;
  real :: Ag_l;
  real :: Ag_rb;
  real :: anbar;
  real :: gsbar;
  real :: w_scale;
  real, parameter :: p_sea = 1.0e5 ! sea level pressure, Pa
  ! soil water stress
  real :: Ed,an_w,gs_w;

  if (is_watch_point()) then
     write(*,*) '####### gs_leuning input #######'
     __DEBUG2__(rad_top, rad_net)
     __DEBUG1__(tl)
     __DEBUG1__(ea)
     __DEBUG1__(lai)
     __DEBUG1__(leaf_age)
     __DEBUG1__(p_surf)
     __DEBUG1__(ws)
     __DEBUG1__(pft)
     __DEBUG1__(ca)
     __DEBUG1__(kappa)
     __DEBUG1__(leaf_wet)
     __DEBUG1__(pt)
     write(*,*) '####### end of ### gs_leuning input #######'
  endif

  b=0.01;
  do1=0.09 ; ! kg/kg
  if (pft < 2) do1=0.15;


  ! Convert Solar influx from W/(m^2s) to mol_of_quanta/(m^2s) PAR,
  ! empirical relationship from McCree is light=rn*0.0000046
  light_top = rad_top*rad_phot;
  par_net   = rad_net*rad_phot;
  
  ! calculate humidity deficit, kg/kg
  call qscomp(tl, p_surf, hl)
  ds = max(hl-ea,0.0)

  ! capgam=0.209/(9000.0*exp(-5000.0*(1.0/288.2-1.0/tl))); - Foley formulation, 1986

!  ko=0.25   *exp(1400.0*(1.0/288.2-1.0/tl))*p_sea/p_surf;
!  kc=0.00015*exp(6000.0*(1.0/288.2-1.0/tl))*p_sea/p_surf;
!  vm=spdata(pft)%Vmax*exp(3000.0*(1.0/288.2-1.0/tl));

! corrected by Weng, 2013-01-17
! Weng, 2013-01-10
  ko=0.248    * exp(35948/Rgas*(1.0/298.2-1.0/tl))*p_sea/p_surf ! Weng, 2013-01-10
  kc=0.000404 * exp(59356/Rgas*(1.0/298.2-1.0/tl))*p_sea/p_surf ! Weng, 2013-01-10
  vm=spdata(pft)%Vmax*exp(24920/Rgas*(1.0/298.2-1.0/tl)) ! / ((layer-1)*1.0+1.0) ! Ea = 33920

  !decrease Vmax due to aging of temperate deciduous leaves 
  !(based on Wilson, Baldocchi and Hanson (2001)."Plant,Cell, and Environment", vol 24, 571-583)
!! Turned off by Weng, 2013-02-01, since we can't trace new leaves
!  if (spdata(pft)%leaf_age_tau>0 .and. leaf_age>spdata(pft)%leaf_age_onset) then
!     vm=vm*exp(-(leaf_age-spdata(pft)%leaf_age_onset)/spdata(pft)%leaf_age_tau)
!  endif

  capgam=0.5*kc/ko*0.21*0.209; ! Farquhar & Caemmerer 1982

  ! Find respiration for the whole canopy layer
  
!  Resp=spdata(pft)%gamma_resp*vm*lai /((layer-1)*1.0+1.0)  ! Weng, 2013-01-17 add '/ ((layer-1)*1.0+1.0)'

! 2014-09-03, for Nitrogen model: resp = D*(A + B*LMA)
! (A+B*LMA) = LNA, D=Vmax/LNA = 25E-6/0.0012 = 0.02 for a standard deciduous species
!! Leaf resp as a function of nitrogen
!  Resp=spdata(pft)%gamma_resp*0.04*spdata(pft)%LNA  & ! basal rate, mol m-2 s-1
!       * exp(24920/Rgas*(1.0/298.2-1.0/tl))         & ! temperature scaled
!       * lai                                        & ! whole canopy
!       /((layer-1)*1.0+1.0)                         !
!! as a function of LMA
!  Resp=(spdata(pft)%gamma_LNbase*spdata(pft)%LNbase+spdata(pft)%gamma_LMA*spdata(pft)%LMA)  & ! basal rate, mol m-2 s-1
!  Resp=spdata(pft)%gamma_LNbase*(2.5*spdata(pft)%LNA-1.5*spdata(pft)%LNbase)     & ! basal rate, mol m-2 s-1
  Resp=spdata(pft)%gamma_LNbase * spdata(pft)%LNA     &
       * exp(24920/Rgas*(1.0/298.2-1.0/tl))           & ! temperature scaled
       * lai        !                                  & ! whole canopy
!       /((layer-1)*1.0+1.0)
! We don't need this equation:
   Resp=Resp/((1.0+exp(0.4*(5.0-tl+TFREEZE)))*(1.0+exp(0.4*(tl-45.0-TFREEZE))));
  
  
  ! ignore the difference in concentrations of CO2 near
  !  the leaf and in the canopy air, rb=0.
 
  Ag_l=0.;
  Ag_rb=0.;
  Ag=0.;
  anbar=-Resp/lai;
  gsbar=b;

 
  ! find the LAI level at which gross photosynthesis rates are equal
  ! only if PAR is positive
  if ( light_top > light_crit ) then
     if (pt==PT_C4) then ! C4 species
        coef0=(1+ds/do1)/spdata(pft)%m_cond;
        ci=(ca+1.6*coef0*capgam)/(1+1.6*coef0); 
!        ci=(ca+1.6*coef0*capgam)/(1+1.6*coef0) * 0.5  ! Weng, to reduce ci
        if (ci>capgam) then
           f2=vm;
           f3=18000.0*vm*ci; ! 18000 or 1800?
!           f3=1800.0 *vm*ci 
       
           dum2=min(f2,f3)
           
           ! find LAI level at which rubisco limited rate is equal to light limited rate
           lai_eq = -log(dum2/(kappa*spdata(pft)%alpha_phot*light_top))/kappa;
           lai_eq = min(max(0.0,lai_eq),lai) ! limit lai_eq to physically possible range

           ! gross photosynthesis for light-limited part of the canopy
           Ag_l   = spdata(pft)%alpha_phot * par_net &
                * (exp(-lai_eq*kappa)-exp(-lai*kappa))/(1-exp(-lai*kappa))
           ! gross photosynthesis for rubisco-limited part of the canopy
           Ag_rb  = dum2*lai_eq

           Ag=(Ag_l+Ag_rb)/ &
             ((1.0+exp(0.4*(5.0-tl+TFREEZE)))*(1.0+exp(0.4*(tl-45.0-TFREEZE))));
           An=Ag-Resp;
           anbar=An/lai;
     
           if(anbar>0.0) then
               gsbar=anbar/(ci-capgam)/coef0;
           endif
        endif ! ci>capgam
     else ! C3 species
        coef0=(1+ds/do1)/spdata(pft)%m_cond;
        coef1=kc*(1.0+0.209/ko);
        ci=(ca+1.6*coef0*capgam)/(1+1.6*coef0);
        f2=vm*(ci-capgam)/(ci+coef1);
        f3=vm/2.;
        dum2=min(f2,f3);
        if (ci>capgam) then
           ! find LAI level at which rubisco limited rate is equal to light limited rate
           lai_eq=-log(dum2*(ci+2.*capgam)/(ci-capgam)/ &
                       (spdata(pft)%alpha_phot*light_top*kappa))/kappa;
           lai_eq = min(max(0.0,lai_eq),lai) ! limit lai_eq to physically possible range

           ! gross photosynthesis for light-limited part of the canopy
           Ag_l   = spdata(pft)%alpha_phot * (ci-capgam)/(ci+2.*capgam) * par_net &
                * (exp(-lai_eq*kappa)-exp(-lai*kappa))/(1-exp(-lai*kappa))
           ! gross photosynthesis for rubisco-limited part of the canopy
           Ag_rb  = dum2*lai_eq

           Ag=(Ag_l+Ag_rb) /((1.0+exp(0.4*(5.0-tl+TFREEZE)))*(1.0+exp(0.4*(tl-45.0-TFREEZE))));
           An=Ag-Resp;
           anbar=An/lai;

           if(anbar>0.0) then
               gsbar=anbar/(ci-capgam)/coef0;
           endif
        endif ! ci>capgam
     endif
  endif ! light is available for photosynthesis
!  write(*,'(1(I4,","),10(E10.4,","))') &
!       layer, light_top, par_net, kappa, lai, lai_eq, ci, capgam, Ag_l, Ag_rb, Ag
  
  an_w=anbar;
  if (an_w > 0.) then
     an_w=an_w*(1-spdata(pft)%wet_leaf_dreg*leaf_wet);
  endif
  
  gs_w=gsbar*(1-spdata(pft)%wet_leaf_dreg*leaf_wet);

  if (gs_w > gs_lim) then
      if(an_w > 0.) an_w = an_w*gs_lim/gs_w;
      gs_w = gs_lim;
  endif

#if 1
  ! find water availability
  ! diagnostic demand

  Ed=gs_w*ds*mol_air/mol_h2o;
  ! the factor mol_air/mol_h2o makes units of gs_w and humidity deficit ds compatible:
  ! ds*mol_air/mol_h2o is the humidity deficit in [mol_h2o/mol_air]

  if (Ed>ws) then
     w_scale=ws/Ed;
     gs_w=w_scale*gs_w;
     if(an_w > 0.0) an_w = an_w*w_scale;
     if(an_w < 0.0.and.gs_w >b) gs_w=b;
     if (is_watch_point()) then
        write(*,*)'#### gs is water-limited'
        __DEBUG1__(w_scale)
        __DEBUG3__(gs_w, an_w, b)
     endif
  endif
  gs=gs_w;
  apot=an_w;
  acl=-Resp/lai;
 
! just for reporting
  w_scale2=min(1.0,ws/Ed) 

#else
! no water limitation on stomata
   gs=gsbar;  
   apot=anbar; 
   acl=-Resp/lai; 
#endif 

   ! finally, convert units of stomatal conductance to m/s from mol/(m2 s) by
   ! multiplying it by a volume of a mole of gas
   gs = gs * Rugas * Tl / p_surf

   if (is_watch_point()) then
      __DEBUG3__(gs, apot, acl)
   endif
end subroutine gs_Leuning

! Add another photosynthesis model to replace 'gs_Leuning'?
! Just for testing responses of gs to elevated CO2

end module vegn_photosynthesis_mod
