module CanopyFluxesMod

#include "shr_assert.h"

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Performs calculation of leaf temperature and surface fluxes.
  ! SoilFluxes then determines soil/snow and ground temperatures and updates the surface
  ! fluxes for the new ground temperature.
  !
  ! !USES:
  use shr_sys_mod           , only : shr_sys_flush
  use shr_kind_mod          , only : r8 => shr_kind_r8
  use shr_log_mod           , only : errMsg => shr_log_errMsg
  use abortutils            , only : endrun
  use elm_varctl            , only : iulog, use_cn, use_lch4, use_c13, use_c14, use_fates
  use elm_varctl            , only : use_hydrstress
  use elm_varpar            , only : nlevgrnd, nlevsno
  use elm_varcon            , only : namep
  use elm_varcon            , only : mm_epsilon
  use elm_varcon            , only : pa_to_kpa
  use pftvarcon             , only : crop, nfixer
  use decompMod             , only : bounds_type
  use PhotosynthesisMod     , only : Photosynthesis, PhotosynthesisTotal, Fractionation, PhotoSynthesisHydraulicStress
  use SoilMoistStressMod    , only : calc_effective_soilporosity, calc_volumetric_h2oliq
  use SoilMoistStressMod    , only : calc_root_moist_stress, set_perchroot_opt
  use SurfaceResistanceMod  , only : do_soilevap_beta
  use VegetationPropertiesType        , only : veg_vp
  use atm2lndType           , only : atm2lnd_type
  use CanopyStateType       , only : canopystate_type, perchroot, perchroot_alt
  use CNStateType           , only : cnstate_type
  use EnergyFluxType        , only : energyflux_type
  use FrictionvelocityType  , only : frictionvel_type
  use SoilStateType         , only : soilstate_type
  use SolarAbsorbedType     , only : solarabs_type
  use SurfaceAlbedoType     , only : surfalb_type
  use CH4Mod                , only : ch4_type
  use PhotosynthesisType    , only : photosyns_type
  use GridcellType          , only : grc_pp
  use TopounitDataType      , only : top_as, top_af
  use ColumnType            , only : col_pp
  use ColumnDataType        , only : col_es, col_ef, col_ws
  use VegetationType        , only : veg_pp
  use VegetationDataType    , only : veg_es, veg_ef, veg_ws, veg_wf

!!! using elm_instMod messes with the compilation order
  use FrictionVelocityMod, only : FrictionVelocity
  use elm_instMod           , only : alm_fates, soil_water_retention_curve
  use TemperatureType , only : temperature_vars
  use perfMod_GPU
  use timeinfoMod
  use spmdmod          , only: masterproc
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: CanopyFluxes

  !------------------------------------------------------------------------------
  ! Derived type for patch-level workspace variables in CanopyFluxes
  ! This groups all per-patch temporary variables to improve cache locality
  ! and reduce memory footprint in the parallel loop
  !------------------------------------------------------------------------------
  type :: worker_type
     ! Atmospheric and surface properties
     real(r8) :: zldis              ! reference height minus zero displacement height [m]
     real(r8) :: ugust_total        ! gustiness including convective velocity [m/s]
     real(r8) :: dth                ! diff of virtual temp. between ref. height and surface
     real(r8) :: dthv               ! diff of vir. poten. temp. between ref. height and surface
     real(r8) :: dqh                ! diff of humidity between ref. height and surface
     real(r8) :: ur                 ! wind speed at reference height [m/s]
     real(r8) :: wc                 ! convective velocity [m/s]
     real(r8) :: dayl_factor        ! scalar (0-1) for daylength effect on Vcmax

     ! Temperature and humidity profiles
     real(r8) :: temp1              ! relation for potential temperature profile
     real(r8) :: temp12m            ! relation for potential temperature profile applied at 2-m
     real(r8) :: temp2              ! relation for specific humidity profile
     real(r8) :: temp22m            ! relation for specific humidity profile applied at 2-m
     real(r8) :: tstar              ! temperature scaling parameter
     real(r8) :: qstar              ! moisture scaling parameter
     real(r8) :: thvstar            ! virtual potential temperature scaling parameter
     
     ! Canopy air properties
     real(r8) :: taf                ! air temperature within canopy space [K]
     real(r8) :: qaf                ! humidity of canopy air [kg/kg]
     
     ! Resistances
     real(r8) :: rb                 ! leaf boundary layer resistance [s/m]
     real(r8) :: rah(2)             ! thermal resistance [s/m] (above/below canopy)
     real(r8) :: raw(2)             ! moisture resistance [s/m] (above/below canopy)
     
     ! Heat conductances
     real(r8) :: wta                ! heat conductance for air [m/s]
     real(r8) :: wtl                ! heat conductance for leaf [m/s]
     real(r8) :: wtg                ! heat conductance for ground [m/s]
     real(r8) :: wta0               ! normalized heat conductance for air [-]
     real(r8) :: wtl0               ! normalized heat conductance for leaf [-]
     real(r8) :: wtg0               ! normalized heat conductance for ground [-]
     real(r8) :: wtal               ! normalized heat conductance for air and leaf [-]
     real(r8) :: wtga               ! normalized heat cond. for air and ground [-]
     real(r8) :: wtshi              ! sensible heat resistance for air, grnd and leaf [-]
     
     ! Latent heat conductances
     real(r8) :: wtaq               ! latent heat conductance for air [m/s]
     real(r8) :: wtlq               ! latent heat conductance for leaf [m/s]
     real(r8) :: wtgq               ! latent heat conductance for ground [m/s]
     real(r8) :: wtaq0              ! normalized latent heat conductance for air [-]
     real(r8) :: wtlq0              ! normalized latent heat conductance for leaf [-]
     real(r8) :: wtgq0              ! normalized heat conductance for ground [-]
     real(r8) :: wtalq              ! normalized latent heat cond. for air and leaf [-]
     real(r8) :: wtgaq              ! normalized latent heat cond. for air and ground [-]
     real(r8) :: wtsqi              ! latent heat resistance for air, grnd and leaf [-]
     
     ! Saturation and vapor pressure
     real(r8) :: el                 ! vapor pressure on leaf surface [pa]
     real(r8) :: deldT              ! derivative of "el" on "t_veg" [pa/K]
     real(r8) :: qsatl              ! leaf specific humidity [kg/kg]
     real(r8) :: qsatldT            ! derivative of "qsatl" on "t_veg"
     real(r8) :: e_ref2m            ! 2 m height surface saturated vapor pressure [Pa]
     real(r8) :: de2mdT             ! derivative of 2 m height surface saturated vapor pressure on t_ref2m
     real(r8) :: qsat_ref2m         ! 2 m height surface saturated specific humidity [kg/kg]
     real(r8) :: dqsat2mdT          ! derivative of 2 m height surface saturated specific humidity on t_ref2m
     real(r8) :: svpts              ! saturation vapor pressure at t_veg (pa)
     real(r8) :: eah                ! canopy air vapor pressure (pa)

     ! Radiation terms
     real(r8) :: air                ! atmos. radiation temporary set
     real(r8) :: bir                ! atmos. radiation temporary set
     real(r8) :: cir                ! atmos. radiation temporary set
     real(r8) :: lw_grnd            ! ground longwave radiation
     
     ! Energy flux derivatives and temporaries
     real(r8) :: dc1                ! derivative of energy flux [W/m2/K]
     real(r8) :: dc2                ! derivative of energy flux [W/m2/K]
     real(r8) :: delt               ! temporary
     real(r8) :: delq               ! temporary
     real(r8) :: delt_snow          ! temperature difference for snow
     real(r8) :: delt_soil          ! temperature difference for soil
     real(r8) :: delt_h2osfc        ! temperature difference for surface water
     real(r8) :: delq_snow          ! humidity difference for snow
     real(r8) :: delq_soil          ! humidity difference for soil
     real(r8) :: delq_h2osfc        ! humidity difference for surface water
     
     ! Iteration variables
     real(r8) :: del                ! absolute change in leaf temp in current iteration [K]
     real(r8) :: del2               ! change in leaf temperature in previous iteration [K]
     real(r8) :: dele               ! change in latent heat flux from leaf [K]
     real(r8) :: dels               ! change in leaf temperature in current iteration [K]
     real(r8) :: det                ! maximum leaf temp. change in two consecutive iter [K]
     real(r8) :: tlbef              ! leaf temperature from previous iteration [K]
     
     ! Evaporation and transpiration
     real(r8) :: rpp                ! fraction of potential evaporation from leaf [-]
     real(r8) :: rppdry             ! fraction of potential evaporation through transp [-]
     real(r8) :: efeb               ! latent heat flux from leaf (previous iter) [mm/s]
     real(r8) :: efeold             ! latent heat flux from leaf (previous iter) [mm/s]
     real(r8) :: efpot              ! potential latent energy flux [kg/m2/s]
     real(r8) :: efe                ! water flux from leaf [mm/s]
     real(r8) :: efsh               ! sensible heat from leaf [mm/s]
     real(r8) :: ecidif             ! excess energies [W/m2]
     real(r8) :: erre               ! balance error
     
     ! Stomatal resistance (for convergence checking)
     real(r8) :: rssun_old          ! previous sunlit stomatal resistance (s/m)
     real(r8) :: rssha_old          ! previous shaded stomatal resistance (s/m)
     
     ! Stability and turbulence
     real(r8) :: obuold             ! Obukhov length scale from previous iteration
     integer  :: nmozsgn            ! number of times stability changes sign
     real(r8) :: fm                 ! needed for BGC only to diagnose 10m wind speed
     
     ! Canopy properties
     real(r8) :: cf                 ! heat transfer coefficient from leaves [-]
     real(r8) :: exp_ai             ! exp(-LSAI)
     real(r8) :: egvf               ! effective green vegetation fraction
     real(r8) :: lt                 ! elai+esai
     
     ! Soil and litter properties
     real(r8) :: csoilb             ! turbulent transfer coefficient over bare soil (unitless)
     real(r8) :: csoilcn            ! interpolated csoilc for less than dense canopies
     real(r8) :: ri                 ! stability parameter for under canopy air (unitless)
     real(r8) :: ricsoilc           ! modified transfer coefficient under dense canopy (unitless)
     real(r8) :: snow_depth_c       ! critical snow depth to cover plant litter (m)
     real(r8) :: rdl                ! dry litter layer resistance for water vapor (s/m)
     real(r8) :: elai_dl            ! exposed (dry) plant litter area index
     real(r8) :: fsno_dl            ! effective snow cover over plant litter
     
     ! Wind stress (implicit coupling)
     real(r8) :: wind_speed0        ! wind speed from atmosphere at start of iteration
     real(r8) :: wind_speed_adj     ! adjusted wind speed for iteration
     real(r8) :: tau                ! stress used in iteration
     real(r8) :: tau_diff           ! difference from previous iteration tau
     real(r8) :: prev_tau           ! previous iteration tau
     real(r8) :: prev_tau_diff      ! previous difference in iteration tau
     
     ! Iteration counters and convergence flags
     integer  :: itlef              ! counter for leaf temperature iteration [-]
     integer  :: itstoma            ! counter for stoma iteration [-]
     integer  :: iter_final         ! number of iterations used
     logical  :: converge_stoma     ! stomatal convergence flag
     logical  :: converge_tveg      ! temperature convergence flag
     real(r8) :: del_gs             ! max difference in stomatal conductance [m/s]

     real(r8) :: err                ! balance error
  end type worker_type

contains

  !------------------------------------------------------------------------------
  subroutine CanopyFluxes(bounds,  num_nolakeurbanp, filter_nolakeurbanp, &
       atm2lnd_vars, canopystate_vars, cnstate_vars, energyflux_vars, &
       frictionvel_vars, soilstate_vars, solarabs_vars, surfalb_vars, &
       ch4_vars, photosyns_vars)
    ! !DESCRIPTION:
    ! 1. Calculates the leaf temperature:
    ! 2. Calculates the leaf fluxes, transpiration, photosynthesis and
    !    updates the dew accumulation due to evaporation.
    !
    ! Method:
    ! Use the Newton-Raphson iteration to solve for the foliage
    ! temperature that balances the surface energy budget:
    !
    ! f(t_veg) = Net radiation - Sensible - Latent = 0
    ! f(t_veg) + d(f)/d(t_veg) * dt_veg = 0     (*)
    !
    ! Note:
    ! (1) In solving for t_veg, t_grnd is given from the previous timestep.
    ! (2) The partial derivatives of aerodynamical resistances, which cannot
    !     be determined analytically, are ignored for d(H)/dT and d(LE)/dT
    ! (3) The weighted stomatal resistance of sunlit and shaded foliage is used
    ! (4) Canopy air temperature and humidity are derived from => Hc + Hg = Ha
    !                                                          => Ec + Eg = Ea
    ! (5) Energy loss is due to: numerical truncation of energy budget equation
    !     (*); and "ecidif" (see the code) which is dropped into the sensible
    !     heat
    ! (6) The convergence criteria: the difference, del = t_veg(n+1)-t_veg(n)
    !     and del2 = t_veg(n)-t_veg(n-1) less than 0.01 K, and the difference
    !     of water flux from the leaf between the iteration step (n+1) and (n)
    !     less than 0.1 W/m2; or the iterative steps over 40.
    !
    ! !USES:
      !$acc routine seq
    use shr_const_mod      , only : SHR_CONST_TKFRZ, SHR_CONST_RGAS
    use shr_flux_mod       , only : shr_flux_update_stress
    use elm_varcon         , only : sb, cpair, hvap, vkc, grav, denice
    use elm_varcon         , only : denh2o, tfrz, csoilc, tlsai_crit, alpha_aero
    use elm_varcon         , only : isecspday, degpsec
    use pftvarcon          , only : irrigated
    use elm_varcon         , only : c14ratio

    !NEW
    use elm_varsur         , only : firrig
    use TopounitType       , only : top_pp
    use QSatMod            , only : QSat
    use FrictionVelocityMod, only : MoninObukIni, &
         implicit_stress, atm_gustiness, force_land_gustiness
    use SoilWaterRetentionCurveMod, only : soil_water_retention_curve_type
    use SurfaceResistanceMod, only : getlblcef
    use PhotosynthesisType, only : photosyns_vars_TimeStepInit
    !
    ! !ARGUMENTS:
    type(bounds_type)         , intent(in)    :: bounds
    integer                   , intent(in)    :: num_nolakeurbanp       ! number of column non-lake, non-urban points in pft filter
    integer                   , intent(in)    :: filter_nolakeurbanp(:) ! patch filter for non-lake, non-urban points
    type(atm2lnd_type)        , intent(inout) :: atm2lnd_vars
    type(canopystate_type)    , intent(inout) :: canopystate_vars
    type(cnstate_type)        , intent(inout) :: cnstate_vars
    type(energyflux_type)     , intent(inout) :: energyflux_vars
    type(frictionvel_type)    , intent(inout) :: frictionvel_vars
    type(solarabs_type)       , intent(inout) :: solarabs_vars
    type(surfalb_type)        , intent(inout) :: surfalb_vars
    type(soilstate_type)      , intent(inout) :: soilstate_vars
    type(ch4_type)            , intent(inout) :: ch4_vars
    type(photosyns_type)      , intent(inout) :: photosyns_vars
    real(r8) :: dtime
    integer  ::  time
    !
    ! !LOCAL VARIABLES:
    real(r8), parameter :: btran0 = 0.0_r8  ! initial value
    real(r8), parameter :: zii = 1000.0_r8  ! convective boundary layer height [m]
    real(r8), parameter :: beta = 1.0_r8    ! coefficient of convective velocity [-]
    real(r8), parameter :: delmax = 1.0_r8  ! maxchange in  leaf temperature [K]
    real(r8), parameter :: dlemin = 0.1_r8  ! max limit for energy flux convergence [w/m2]
    real(r8), parameter :: dtmin = 0.01_r8  ! max limit for temperature convergence [K]
    real(r8), parameter :: dtaumin = 0.01_r8! max limit for stress convergence [Pa]
    integer , parameter :: itmax = 41       ! maximum number of iteration [-]
    integer , parameter :: itmin = 3        ! minimum number of iteration [-]
    real(r8), parameter :: irrig_min_lai = 0.0_r8           ! Minimum LAI for irrigation
    real(r8), parameter :: irrig_btran_thresh = 0.999999_r8 ! Irrigate when btran falls below 0.999999 rather than 1 to allow for round-off error
    integer , parameter :: irrig_start_time = isecspday/4   ! (6AM) Time of day to check whether we need irrigation, seconds (0 = midnight).

    ! We start applying the irrigation in the time step FOLLOWING this time,
    ! since we won't begin irrigating until the next call to CanopyHydrology
    ! Desired amount of time to irrigate per day (sec). Actual time may
    ! differ if this is not a multiple of dtime. Irrigation won't work properly
    ! if dtime > secsperday

    integer , parameter :: irrig_length = isecspday/6     ! 4 hours irrigation

    ! Determines target soil moisture level for irrigation. If h2osoi_liq_so
    ! is the soil moisture level at which stomata are fully open and
    ! h2osoi_liq_sat is the soil moisture level at saturation (eff_porosity),
    ! then the target soil moisture level is
    !     (h2osoi_liq_so + irrig_factor*(h2osoi_liq_sat - h2osoi_liq_so)).
    ! A value of 0 means that the target soil moisture level is h2osoi_liq_so.

    ! A value of 1 means that the target soil moisture level is h2osoi_liq_sat
    real(r8), parameter :: irrig_factor = 0.7_r8

    !added by K.Sakaguchi for litter resistance
    real(r8), parameter :: lai_dl = 0.5_r8           ! placeholder for (dry) plant litter area index (m2/m2)
    real(r8), parameter :: z_dl = 0.05_r8            ! placeholder for (dry) litter layer thickness (m)

    !added by K.Sakaguchi for stability formulation
    real(r8), parameter :: ria  = 0.5_r8             ! free parameter for stable formulation (currently = 0.5, "gamma" in Sakaguchi&Zeng,2008)

    ! Arrays that remain (used before/after loop or need to persist)
    real(r8) :: zldis(bounds%begp:bounds%endp)       ! reference height "minus" zero displacement height [m]
    real(r8) :: err(bounds%begp:bounds%endp)         ! balance error
    
    real(r8) :: dt_veg(bounds%begp:bounds%endp)      ! change in t_veg, last iteration (Kelvin)
    logical  :: check_for_irrig(bounds%begp:bounds%endp) ! where do we need to check soil moisture to see if we need to irrigate?
    logical  :: frozen_soil(bounds%begp:bounds%endp)     ! set to true if we have encountered a frozen soil layer

    ! Scalars and indices
    real(r8) :: cf_bare                              ! heat transfer coefficient from bare ground [-]
    integer  :: j                                    ! soil/snow level index
    integer  :: p                                    ! patch index
    integer  :: c                                    ! column index
    integer  :: l                                    ! landunit index
    integer  :: t                                    ! topounit index
    integer  :: tpu_ind                              ! index of topounit to grid
    integer  :: g                                    ! gridcell index
    integer  :: fp                                   ! lake filter pft index
    integer  :: fn_noveg                             ! number of values in bare ground pft filter
    integer  :: filterp_noveg(bounds%endp-bounds%begp+1) ! bare ground pft filter
    integer  :: fn                                   ! number of values in vegetated pft filter
    integer  :: filterp(bounds%endp-bounds%begp+1)   ! vegetated pft filter
    integer  :: f                                    ! filter index
    logical  :: found                                ! error flag for canopy above forcing hgt
    integer  :: index                                ! patch index for error
    integer  :: local_time                           ! local time at start of time step (seconds after solar midnight)
    integer  :: seconds_since_irrig_start_time
    integer  :: irrig_nsteps_per_day                 ! number of time steps per day in which we irrigate
    integer  :: jtop(bounds%begc:bounds%endc)        ! lbning
    integer  :: filterc_tmp(bounds%endp-bounds%begp+1)   ! temporary variable
    integer  :: ft                                   ! plant functional type index
    real(r8) :: temprootr
    integer  :: iv

    real(r8) :: deficit            ! soil moisture deficit [kg/m2]
    real(r8) :: vol_liq_so         ! partial volume of liquid water for smp_node = smpso
    real(r8) :: h2osoi_liq_so      ! liquid water corresponding to vol_liq_so [kg/m2]
    real(r8) :: h2osoi_liq_sat     ! liquid water at eff_porosity [kg/m2]

    
    character(len=64) :: event !! timing event
    
    ! Worker type for patch-level temporary variables
    type(worker_type) :: w                           ! Patch workspace variables
    
    ! Indices for raw and rah
    integer, parameter :: above_canopy = 1         ! Above canopy
    integer, parameter :: below_canopy = 2         ! Below canopy

    ! Lower bound for VPD (based on CLM)
    real(r8), parameter :: vpd_min = 50._r8

    ! We set the minum allowable difference in the conductance iteration
    ! to be equal to the maximum allowable stomatal resistance (this number is from fates)
    real(r8),parameter :: max_del_gs =  1._r8/2.e8_r8   ! [m/s]

    integer, parameter  :: itmax_stomata = 3

    logical, parameter :: do_b4b = .false.  ! Set this true to reproduce results before
                                           ! refactoring the patch-loops
    
    !------------------------------------------------------------------------------

    associate(                                                               &
         snl                  => col_pp%snl                                   , & ! Input:  [integer  (:)   ]  number of snow layers
         dayl                 => grc_pp%dayl                                  , & ! Input:  [real(r8) (:)   ]  daylength (s)
         max_dayl             => grc_pp%max_dayl                              , & ! Input:  [real(r8) (:)   ]  maximum daylength for this grid cell (s)

         forc_lwrad           => top_af%lwrad                              , & ! Input:  [real(r8) (:)   ]  downward infrared (longwave) radiation (W/m**2)                       
         forc_q               => top_as%qbot                               , & ! Input:  [real(r8) (:)   ]  atmospheric specific humidity (kg/kg)                                 
         forc_pbot            => top_as%pbot                               , & ! Input:  [real(r8) (:)   ]  atmospheric pressure (Pa)                                             
         forc_th              => top_as%thbot                              , & ! Input:  [real(r8) (:)   ]  atmospheric potential temperature (Kelvin)                            
         forc_rho             => top_as%rhobot                             , & ! Input:  [real(r8) (:)   ]  air density (kg/m**3)                                                     
         forc_t               => top_as%tbot                               , & ! Input:  [real(r8) (:)   ]  atmospheric temperature (Kelvin)                                      
         forc_u               => top_as%ubot                               , & ! Input:  [real(r8) (:)   ]  atmospheric wind speed in east direction (m/s)                        
         forc_v               => top_as%vbot                               , & ! Input:  [real(r8) (:)   ]  atmospheric wind speed in north direction (m/s)                       
         wsresp               => top_as%wsresp                             , & ! Input:  [real(r8) (:)   ]  response of wind to surface stress (m/s/Pa)
         tau_est              => top_as%tau_est                            , & ! Input:  [real(r8) (:)   ]  approximate atmosphere change to zonal wind (m/s)
         ugust                => top_as%ugust                              , & ! Input:  [real(r8) (:)   ]  gustiness from atmosphere (m/s)
         forc_pco2            => top_as%pco2bot                            , & ! Input:  [real(r8) (:)   ]  partial pressure co2 (Pa)                                             
         forc_po2             => top_as%po2bot                             , & ! Input:  [real(r8) (:)   ]  partial pressure o2 (Pa)                                              

         dleaf                => veg_vp%dleaf                          , & ! Input:  [real(r8) (:)   ]  characteristic leaf dimension (m)
         smpso                => veg_vp%smpso                          , & ! Input:  [real(r8) (:)   ]  soil water potential at full stomatal opening (mm)
         smpsc                => veg_vp%smpsc                          , & ! Input:  [real(r8) (:)   ]  soil water potential at full stomatal closure (mm)

         htvp                 => col_ef%htvp                  , & ! Input:  [real(r8) (:)   ]  latent heat of evaporation (/sublimation) [J/kg] (constant)

         sabv                 => solarabs_vars%sabv_patch                  , & ! Input:  [real(r8) (:)   ]  solar radiation absorbed by vegetation (W/m**2)

         lbl_rsc_h2o          => canopystate_vars%lbl_rsc_h2o_patch        , & ! Output: [real(r8) (:)   ] laminar boundary layer resistance for h2o
         frac_veg_nosno       => canopystate_vars%frac_veg_nosno_patch     , & ! Input:  [integer  (:)   ]  fraction of vegetation not covered by snow (0 OR 1) [-]
         elai                 => canopystate_vars%elai_patch               , & ! Input:  [real(r8) (:)   ]  one-sided leaf area index with burying by snow
         esai                 => canopystate_vars%esai_patch               , & ! Input:  [real(r8) (:)   ]  one-sided stem area index with burying by snow
         laisun               => canopystate_vars%laisun_patch             , & ! Input:  [real(r8) (:)   ]  sunlit leaf area
         laisha               => canopystate_vars%laisha_patch             , & ! Input:  [real(r8) (:)   ]  shaded leaf area
         displa               => canopystate_vars%displa_patch             , & ! Input:  [real(r8) (:)   ]  displacement height (m)
         htop                 => canopystate_vars%htop_patch               , & ! Input:  [real(r8) (:)   ]  canopy top(m)
         altmax_lastyear_indx => canopystate_vars%altmax_lastyear_indx_col , & ! Input:  [integer  (:)   ]  prior year maximum annual depth of thaw
         altmax_indx          => canopystate_vars%altmax_indx_col          , & ! Input:  [integer  (:)   ]  maximum annual depth of thaw

         dleaf_patch          => canopystate_vars%dleaf_patch                 , & ! Output: [real(r8) (:)   ]  mean leaf diameter for this patch/pft
         watsat               => soilstate_vars%watsat_col                 , & ! Input:  [real(r8) (:,:) ]  volumetric soil water at saturation (porosity)   (constant)
         watdry               => soilstate_vars%watdry_col                 , & ! Input:  [real(r8) (:,:) ]  btran parameter for btran=0                      (constant)
         watopt               => soilstate_vars%watopt_col                 , & ! Input:  [real(r8) (:,:) ]  btran parameter for btran=1                      (constant)
         eff_porosity         => soilstate_vars%eff_porosity_col           , & ! Output: [real(r8) (:,:) ]  effective soil porosity

         sucsat               => soilstate_vars%sucsat_col                 , & ! Input:  [real(r8) (:,:) ]  minimum soil suction (mm)                        (constant)
         bsw                  => soilstate_vars%bsw_col                    , & ! Input:  [real(r8) (:,:) ]  Clapp and Hornberger "b"                         (constant)
         rootfr               => soilstate_vars%rootfr_patch               , & ! Input:  [real(r8) (:,:) ]  fraction of roots in each soil layer
         soilbeta             => soilstate_vars%soilbeta_col               , & ! Input:  [real(r8) (:)   ]  soil wetness relative to field capacity
         rootr                => soilstate_vars%rootr_patch                , & ! Output: [real(r8) (:,:) ]  effective fraction of roots in each soil layer

         forc_hgt_u_patch     => frictionvel_vars%forc_hgt_u_patch         , & ! Input:  [real(r8) (:)   ]  observational height of wind at pft level [m]
         z0mg                 => frictionvel_vars%z0mg_col                 , & ! Input:  [real(r8) (:)   ]  roughness length of ground, momentum [m]
         ram1                 => frictionvel_vars%ram1_patch               , & ! Output: [real(r8) (:)   ]  aerodynamical resistance (s/m)
         z0mv                 => frictionvel_vars%z0mv_patch               , & ! Output: [real(r8) (:)   ]  roughness length over vegetation, momentum [m]
         z0hv                 => frictionvel_vars%z0hv_patch               , & ! Output: [real(r8) (:)   ]  roughness length over vegetation, sensible heat [m]
         z0qv                 => frictionvel_vars%z0qv_patch               , & ! Output: [real(r8) (:)   ]  roughness length over vegetation, latent heat [m]
         rb1                  => frictionvel_vars%rb1_patch                , & ! Output: [real(r8) (:)   ]  boundary layer resistance (s/m)
         num_iter             => frictionvel_vars%num_iter_patch           , & ! Output: number of iterations required

         t_h2osfc             => col_es%t_h2osfc             , & ! Input:  [real(r8) (:)   ]  surface water temperature
         t_soisno             => col_es%t_soisno             , & ! Input:  [real(r8) (:,:) ]  soil temperature (Kelvin)
         t_grnd               => col_es%t_grnd               , & ! Input:  [real(r8) (:)   ]  ground surface temperature [K]
         thv                  => col_es%thv                  , & ! Input:  [real(r8) (:)   ]  virtual potential temperature (kelvin)
         thm                  => veg_es%thm                  , & ! Input:  [real(r8) (:)   ]  intermediate variable (forc_t+0.0098*forc_hgt_t_patch)
         emv                  => veg_es%emv                  , & ! Input:  [real(r8) (:)   ]  vegetation emissivity
         emg                  => col_es%emg                  , & ! Input:  [real(r8) (:)   ]  vegetation emissivity
         t_veg                => veg_es%t_veg                , & ! Output: [real(r8) (:)   ]  vegetation temperature (Kelvin)
         t_ref2m              => veg_es%t_ref2m              , & ! Output: [real(r8) (:)   ]  2 m height surface air temperature (Kelvin)
         t_ref2m_r            => veg_es%t_ref2m_r            , & ! Output: [real(r8) (:)   ]  Rural 2 m height surface air temperature (Kelvin)

         frac_h2osfc          => col_ws%frac_h2osfc           , & ! Input:  [real(r8) (:)   ]  fraction of surface water
         fwet                 => veg_ws%fwet                , & ! Input:  [real(r8) (:)   ]  fraction of canopy that is wet (0 to 1)
         fdry                 => veg_ws%fdry                , & ! Input:  [real(r8) (:)   ]  fraction of foliage that is green and dry [-]
         frac_sno             => col_ws%frac_sno_eff          , & ! Input:  [real(r8) (:)   ]  fraction of ground covered by snow (0 to 1)
         snow_depth           => col_ws%snow_depth            , & ! Input:  [real(r8) (:)   ]  snow height (m)
         qg_snow              => col_ws%qg_snow               , & ! Input:  [real(r8) (:)   ]  specific humidity at snow surface [kg/kg]
         qg_soil              => col_ws%qg_soil               , & ! Input:  [real(r8) (:)   ]  specific humidity at soil surface [kg/kg]
         qg_h2osfc            => col_ws%qg_h2osfc             , & ! Input:  [real(r8) (:)   ]  specific humidity at h2osfc surface [kg/kg]
         qg                   => col_ws%qg                    , & ! Input:  [real(r8) (:)   ]  specific humidity at ground surface [kg/kg]
         dqgdT                => col_ws%dqgdT                 , & ! Input:  [real(r8) (:)   ]  temperature derivative of "qg"
         h2osoi_ice           => col_ws%h2osoi_ice            , & ! Input:  [real(r8) (:,:) ]  ice lens (kg/m2)
         h2osoi_vol           => col_ws%h2osoi_vol            , & ! Input:  [real(r8) (:,:) ]  volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3] by F. Li and S. Levis
         h2osoi_liq           => col_ws%h2osoi_liq            , & ! Input:  [real(r8) (:,:) ]  liquid water (kg/m2)
         h2osoi_liqvol        => col_ws%h2osoi_liqvol         , & ! Output: [real(r8) (:,:) ]  volumetric liquid water (v/v)

         h2ocan               => veg_ws%h2ocan              , & ! Output: [real(r8) (:)   ]  canopy water (mm H2O)
         q_ref2m              => veg_ws%q_ref2m             , & ! Output: [real(r8) (:)   ]  2 m height surface specific humidity (kg/kg)
         rh_ref2m_r           => veg_ws%rh_ref2m_r          , & ! Output: [real(r8) (:)   ]  Rural 2 m height surface relative humidity (%)
         rh_ref2m             => veg_ws%rh_ref2m            , & ! Output: [real(r8) (:)   ]  2 m height surface relative humidity (%)
         rhaf                 => veg_ws%rh_af               , & ! Output: [real(r8) (:)   ]  fractional humidity of canopy air [dimensionless]

         !pgwgt                => veg_pp%wtgcell              , & ! Input:  [integer  (:)   ]  pft's weight in gridcell
         n_irrig_steps_left   => veg_wf%n_irrig_steps_left   , & ! Output: [integer  (:)   ]  number of time steps for which we still need to irrigate today
         irrig_rate           => veg_wf%irrig_rate           , & ! Output: [real(r8) (:)   ]  current irrigation rate [mm/s]
         qflx_tran_veg        => veg_wf%qflx_tran_veg        , & ! Output: [real(r8) (:)   ]  vegetation transpiration (mm H2O/s) (+ = to atm)
         qflx_evap_veg        => veg_wf%qflx_evap_veg        , & ! Output: [real(r8) (:)   ]  vegetation evaporation (mm H2O/s) (+ = to atm)
         qflx_evap_soi        => veg_wf%qflx_evap_soi        , & ! Output: [real(r8) (:)   ]  soil evaporation (mm H2O/s) (+ = to atm)
         qflx_ev_snow         => veg_wf%qflx_ev_snow         , & ! Output: [real(r8) (:)   ]  evaporation flux from snow (W/m**2) [+ to atm]
         qflx_ev_soil         => veg_wf%qflx_ev_soil         , & ! Output: [real(r8) (:)   ]  evaporation flux from soil (W/m**2) [+ to atm]
         qflx_ev_h2osfc       => veg_wf%qflx_ev_h2osfc       , & ! Output: [real(r8) (:)   ]  evaporation flux from h2osfc (W/m**2) [+ to atm]

         rssun                => photosyns_vars%rssun_patch                , & ! Output: [real(r8) (:)   ]  leaf sunlit stomatal resistance (s/m) (output from Photosynthesis)
         rssha                => photosyns_vars%rssha_patch                , & ! Output: [real(r8) (:)   ]  leaf shaded stomatal resistance (s/m) (output from Photosynthesis)
         !rssun_old            => photosyns_vars%rssun_old_patch            , & ! Output: [real(r8) (:)   ]  previous leaf sunlit stomatal resistance (s/m) (output from Photosynthesis)
         !rssha_old            => photosyns_vars%rssha_old_patch            , & ! Output: [real(r8) (:)   ]  previous leaf shaded stomatal resistance (s/m) (output from Photosynthesis)

         grnd_ch4_cond        => ch4_vars%grnd_ch4_cond_patch              , & ! Output: [real(r8) (:)   ]  tracer conductance for boundary layer [m/s]

         btran2               => energyflux_vars%btran2_patch              , & ! Output: [real(r8) (:)   ]  F. Li and S. Levis
         btran                => energyflux_vars%btran_patch               , & ! Output: [real(r8) (:)   ]  transpiration wetness factor (0 to 1)
         rresis               => energyflux_vars%rresis_patch              , & ! Output: [real(r8) (:,:) ]  root resistance by layer (0-1)  (nlevgrnd)
         taux                 => veg_ef%taux                , & ! Output: [real(r8) (:)   ]  wind (shear) stress: e-w (kg/m/s**2)
         tauy                 => veg_ef%tauy                , & ! Output: [real(r8) (:)   ]  wind (shear) stress: n-s (kg/m/s**2)
         canopy_cond          => energyflux_vars%canopy_cond_patch         , & ! Output: [real(r8) (:)   ]  tracer conductance for canopy [m/s]
         cgrnds               => veg_ef%cgrnds              , & ! Output: [real(r8) (:)   ]  deriv. of soil sensible heat flux wrt soil temp [w/m2/k]
         cgrndl               => veg_ef%cgrndl              , & ! Output: [real(r8) (:)   ]  deriv. of soil latent heat flux wrt soil temp [w/m**2/k]
         dlrad                => veg_ef%dlrad               , & ! Output: [real(r8) (:)   ]  downward longwave radiation below the canopy [W/m2]
         ulrad                => veg_ef%ulrad               , & ! Output: [real(r8) (:)   ]  upward longwave radiation above the canopy [W/m2]
         cgrnd                => veg_ef%cgrnd               , & ! Output: [real(r8) (:)   ]  deriv. of soil energy flux wrt to soil temp [w/m2/k]
         eflx_sh_snow         => veg_ef%eflx_sh_snow        , & ! Output: [real(r8) (:)   ]  sensible heat flux from snow (W/m**2) [+ to atm]
         eflx_sh_h2osfc       => veg_ef%eflx_sh_h2osfc      , & ! Output: [real(r8) (:)   ]  sensible heat flux from soil (W/m**2) [+ to atm]
         eflx_sh_soil         => veg_ef%eflx_sh_soil        , & ! Output: [real(r8) (:)   ]  sensible heat flux from soil (W/m**2) [+ to atm]
         eflx_sh_veg          => veg_ef%eflx_sh_veg         , & ! Output: [real(r8) (:)   ]  sensible heat flux from leaves (W/m**2) [+ to atm]
         eflx_sh_grnd         => veg_ef%eflx_sh_grnd        , & ! Output: [real(r8) (:)   ]  sensible heat flux from ground (W/m**2) [+ to atm]
         rah_above            => frictionvel_vars%rah_above_patch , & ! Output: [real(r8) (:)   ]  above-canopy sensible heat flux resistance [s/m]
         rah_below            => frictionvel_vars%rah_above_patch , & ! Output: [real(r8) (:)   ]  below-canopy sensible heat flux resistance [s/m]
         raw_above            => frictionvel_vars%raw_below_patch , & ! Output: [real(r8) (:)   ]  above-canopy water vapour flux resistance [s/m]
         raw_below            => frictionvel_vars%raw_below_patch , & ! Output: [real(r8) (:)   ]  below-canopy water vapour flux resistance [s/m]
         ustar                => frictionvel_vars%ustar_patch     , & ! Output: [real(r8) (:)   ]  friction velocity [m/s]
         um                   => frictionvel_vars%um_patch        , & ! Output: [real(r8) (:)   ]  wind speed including the stablity effect [m/s]
         uaf                  => frictionvel_vars%uaf_patch       , & ! Output: [real(r8) (:)   ]  canopy air wind speed [m/s]
         taf                  => frictionvel_vars%taf_patch       , & ! Output: [real(r8) (:)   ]  canopy air temperature [K]
         qaf                  => frictionvel_vars%qaf_patch       , & ! Output: [real(r8) (:)   ]  canopy air specific humidity [kg/kg]
         obu                  => frictionvel_vars%obu_patch       , & ! Output: [real(r8) (:)   ]  Obukhov length scale [m]
         zeta                 => frictionvel_vars%zeta_patch      , & ! Output: [real(r8) (:)   ]  dimensionless stability parameter 
         vpd                  => frictionvel_vars%vpd_patch       , & ! Output: [real(r8) (:)   ]  vapour pressure deficit [kPa]
         begp                 => bounds%begp                               , &
         endp                 => bounds%endp                                 &
         )

      ! Determine step size
      dtime = dtime_mod
      !yr = year_curr; mon = mon_curr; day = day_curr;
      time = secs_curr;

      irrig_nsteps_per_day = ((irrig_length + (dtime - 1))/dtime)  ! round up
      ! First - set the following values over points where frac vegetation NOT covered by snow is zero
      ! (e.g. btran, t_veg, rootr, rresis)
      do fp = 1,num_nolakeurbanp
         p = filter_nolakeurbanp(fp)
         c = veg_pp%column(p)
         t = veg_pp%topounit(p)
         if (frac_veg_nosno(p) == 0) then
            btran(p) = 0._r8
            t_veg(p) = forc_t(t)
            cf_bare  = forc_pbot(t)/(SHR_CONST_RGAS*0.001_r8*thm(p))*1.e06_r8
            rssun(p) = 1._r8/1.e15_r8 * cf_bare
            rssha(p) = 1._r8/1.e15_r8 * cf_bare
            lbl_rsc_h2o(p)=0._r8
            do j = 1, nlevgrnd
               rootr(p,j)  = 0._r8
               rresis(p,j) = 0._r8
            end do
         end if
      end do
      ! -----------------------------------------------------------------
      ! Time step initialization of photosynthesis variables
      ! -----------------------------------------------------------------

      call photosyns_vars_TimeStepInit(photosyns_vars,bounds)

      ! -----------------------------------------------------------------
      ! Filter patches where frac_veg_nosno IS NON-ZERO
      ! -----------------------------------------------------------------
      fn = 0
      do fp = 1,num_nolakeurbanp
         p = filter_nolakeurbanp(fp)
         if (frac_veg_nosno(p) /= 0) then
            fn = fn + 1
            filterp(fn) = p
         end if
      end do

#ifndef _OPENACC
      if (use_fates) then
         call alm_fates%prep_canopyfluxes( bounds )
      end if
#endif

      rb1(begp:endp) = 0._r8
      !assign the temporary filter
      do f = 1, fn
         p = filterp(f)
         btran(p)   = btran0
         btran2(p)  = btran0
         filterc_tmp(f)=veg_pp%column(p)
      enddo

      !compute effective soil porosity
      call calc_effective_soilporosity(bounds,                          &
           ubj = nlevgrnd,                                              &
           numf = fn,                                                   &
           filter = filterc_tmp(1:fn),                                  &
           watsat = watsat(bounds%begc:bounds%endc, 1:nlevgrnd),        &
           h2osoi_ice = h2osoi_ice(bounds%begc:bounds%endc,1:nlevgrnd), &
           denice = denice,                                             &
           eff_por=eff_porosity(bounds%begc:bounds%endc, 1:nlevgrnd) )

      !compute volumetric liquid water content
      jtop(bounds%begc:bounds%endc) = 1

      call calc_volumetric_h2oliq(bounds,                                    &
           jtop = jtop(bounds%begc:bounds%endc),                             &
           lbj = 1,                                                          &
           ubj = nlevgrnd,                                                   &
           numf = fn,                                                        &
           filter = filterc_tmp(1:fn),                                       &
           eff_porosity = eff_porosity(bounds%begc:bounds%endc, 1:nlevgrnd), &
           h2osoi_liq = h2osoi_liq(bounds%begc:bounds%endc, 1:nlevgrnd),     &
           denh2o = denh2o,                                                  &
           vol_liq = h2osoi_liqvol(bounds%begc:bounds%endc, 1:nlevgrnd) )

      !set up perchroot options
      call set_perchroot_opt(perchroot, perchroot_alt)
      ! --------------------------------------------------------------------------
      ! if this is a FATES simulation
      ! ask fates to calculate btran functions and distribution of uptake
      ! this will require boundary conditions from CLM, boundary conditions which
      ! may only be available from a smaller subset of patches that meet the
      ! exposed veg.
      ! calc_root_moist_stress already calculated root soil water stress 'rresis'
      ! this is the input boundary condition to calculate the transpiration
      ! wetness factor btran and the root weighting factors for FATES.  These
      ! values require knowledge of the belowground root structure.
      ! --------------------------------------------------------------------------

      if(use_fates)then
#ifndef _OPENACC
         call alm_fates%wrap_btran(bounds, fn, filterc_tmp(1:fn), soilstate_vars, &
               energyflux_vars, soil_water_retention_curve)
#endif
      else
         !calculate root moisture stress
         call calc_root_moist_stress(bounds,     &
              nlevgrnd = nlevgrnd,               &
              fn = fn,                           &
              filterp = filterp,                 &
              canopystate_vars=canopystate_vars, &
              energyflux_vars=energyflux_vars,   &
              soilstate_vars=soilstate_vars      &
              )

      end if !use_fates
      ! Determine if irrigation is needed (over irrigated soil columns)

      ! First, determine in what grid cells we need to bother 'measuring' soil water, to see if we need irrigation
      ! Also set n_irrig_steps_left for these grid cells
      ! n_irrig_steps_left(p) > 0 is ok even if irrig_rate(p) ends up = 0
      ! in this case, we'll irrigate by 0 for the given number of time steps

      do f = 1, fn
         p = filterp(f)
         c = veg_pp%column(p)
         g = veg_pp%gridcell(p)
         if ( .not.veg_pp%is_fates(p)             .and. &
              irrigated(veg_pp%itype(p)) == 1._r8 .and. &
              elai(p) > irrig_min_lai          .and. &
              btran(p) < irrig_btran_thresh ) then

            ! see if it's the right time of day to start irrigating:
            local_time = modulo(time + nint(grc_pp%londeg(g)/degpsec), isecspday)
            seconds_since_irrig_start_time = modulo(local_time - irrig_start_time, isecspday)
            if (seconds_since_irrig_start_time < dtime) then
               ! it's time to start irrigating
               check_for_irrig(p)    = .true.
               n_irrig_steps_left(p) = irrig_nsteps_per_day
               irrig_rate(p)         = 0._r8  ! reset; we'll add to this later
            else
               check_for_irrig(p)    = .false.
            end if
         else  ! non-irrig pft or elai<=irrig_min_lai or btran>irrig_btran_thresh
            check_for_irrig(p)       = .false.
         end if

      end do


      ! Now 'measure' soil water for the grid cells identified above and see if the
      ! soil is dry enough to warrant irrigation
      ! (Note: frozen_soil could probably be a column-level variable, but that would be
      ! slightly less robust to potential future modifications)
      ! This should not be operating on FATES patches (see is_fates filter above, pushes
      ! check_for_irrig = false
      frozen_soil(bounds%begp : bounds%endp) = .false.
      do j = 1,nlevgrnd
         do f = 1, fn
            p = filterp(f)
            c = veg_pp%column(p)
            t = veg_pp%topounit(p)
            tpu_ind = top_pp%topo_grc_ind(t)  !Get topounit index on the grid
            g = veg_pp%gridcell(p)
            if (check_for_irrig(p) .and. .not. frozen_soil(p)) then
               ! if level L was frozen, then we don't look at any levels below L
               if (t_soisno(c,j) <= SHR_CONST_TKFRZ) then
                  frozen_soil(p) = .true.
               else if (rootfr(p,j) > 0._r8) then
                  ! determine soil water deficit in this layer:

                  ! Calculate vol_liq_so - i.e., vol_liq at which smp_node = smpso - by inverting the above equations
                  ! for the root resistance factors
                  vol_liq_so   = eff_porosity(c,j) * (-smpso(veg_pp%itype(p))/sucsat(c,j))**(-1/bsw(c,j))

                  ! Translate vol_liq_so and eff_porosity into h2osoi_liq_so and h2osoi_liq_sat and calculate deficit
                  h2osoi_liq_so  = vol_liq_so * denh2o * col_pp%dz(c,j)
                  h2osoi_liq_sat = eff_porosity(c,j) * denh2o * col_pp%dz(c,j)
                  deficit        = max((h2osoi_liq_so + firrig(g,tpu_ind)*(h2osoi_liq_sat - h2osoi_liq_so)) - h2osoi_liq(c,j), 0._r8)

                  ! Add deficit to irrig_rate, converting units from mm to mm/sec
                  irrig_rate(p)  = irrig_rate(p) + deficit/(dtime*irrig_nsteps_per_day)

               end if  ! else if (rootfr(p,j) > 0)
            end if     ! if (check_for_irrig(p) .and. .not. frozen_soil(p))
         end do        ! do f
      end do           ! do j

      
      event = 'can_iter'
      call t_start_lnd(event)
      
      ! Use adaptive scheduling based on workload variability:
      ! - FATES: dynamic(1) for highly variable photosynthesis costs
      ! - Non-FATES: guided for lower overhead with more uniform work
      !if (use_fates) then
      !$OMP PARALLEL DO PRIVATE (f,p,c,t,g,w) SCHEDULE(DYNAMIC, 1) if (use_fates)
      !else
      !  OMP PARALLEL DO PRIVATE (f,p,c,t,g,w) SCHEDULE(GUIDED)
      !end if

      do_patch: do f = 1, fn
         p = filterp(f)
         c = veg_pp%column(p)
         t = veg_pp%topounit(p)
         g = veg_pp%gridcell(p)

         ! Initialize worker variables
         w%del = 0._r8
         w%efeb = 0._r8
         w%wtlq0 = 0._r8
         w%wtalq = 0._r8
         w%wtgq = 0._r8
         w%wtaq0 = 0._r8
         w%obuold = 0._r8
         w%nmozsgn = 0

         ! Modify aerodynamic parameters for sparse/dense canopy (X. Zeng)
         w%lt = min(elai(p)+esai(p), tlsai_crit)
         w%egvf =(1._r8 - alpha_aero * exp(-w%lt)) / (1._r8 - alpha_aero * exp(-tlsai_crit))
         displa(p) = w%egvf * displa(p)
         z0mv(p)   = exp(w%egvf * log(z0mv(p)) + (1._r8 - w%egvf) * log(z0mg(c)))
         z0hv(p)   = z0mv(p)
         z0qv(p)   = z0mv(p)

         ! Net absorbed longwave radiation by canopy and ground
         ! =air+bir*t_veg**4+cir*t_grnd(c)**4

         w%air =   emv(p) * (1._r8+(1._r8-emv(p))*(1._r8-emg(c))) * forc_lwrad(t)
         w%bir = - (2._r8-emv(p)*(1._r8-emg(c))) * emv(p) * sb
         w%cir =   emv(p)*emg(c)*sb

         ! Saturated vapor pressure, specific humidity, and their derivatives
         ! at the leaf surface

         call QSat (t_veg(p), forc_pbot(t), w%el, w%deldT, w%qsatl, w%qsatldT)

         ! Initialize flux profile

         w%taf = (t_grnd(c) + thm(p))/2._r8
         w%qaf = (forc_q(t)+qg(c))/2._r8

         ! Initialize winds for iteration.
         if (implicit_stress) then
            w%wind_speed0 = max(0.01_r8, hypot(forc_u(t), forc_v(t)))
            w%wind_speed_adj = w%wind_speed0
            w%ur = max(1.0_r8, sqrt(w%wind_speed_adj**2 + ugust(t)**2))

            w%prev_tau = tau_est(t)
         else
            w%ur = max(1.0_r8,sqrt(forc_u(t)*forc_u(t)+forc_v(t)*forc_v(t)+ugust(t)*ugust(t)))
         end if
         w%tau_diff = 1.e100_r8
         w%ugust_total = ugust(t)

         w%dth = thm(p)-w%taf
         w%dqh = forc_q(t)-w%qaf
         w%delq = qg(c) - w%qaf
         w%dthv = w%dth*(1._r8+0.61_r8*forc_q(t))+0.61_r8*forc_th(t)*w%dqh
         zldis(p) = forc_hgt_u_patch(p) - displa(p)

         ! Initialize Obukhov length scale and wind speed
         call MoninObukIni(w%ur, thv(c), w%dthv, zldis(p), z0mv(p), um(p), obu(p))
         
         num_iter(p) = 0
         !rssun(p) is carried over from previous time-step
         !rssha(p) is carried over from previous time-step
         w%rssun_old = -100._r8
         w%rssha_old = -100._r8
        

         ! Begin stability iterations
         ! We evaluate temperature convergence inside
         ! of stomatal conductance convergence. We nest these
         ! loops to minimize stomatal conductance calculations

         w%itstoma = 0
         w%converge_stoma = .false.
         iterate_stoma: do while(.not.w%converge_stoma) 

            ! Set counter for leaf temperature iteration (itlef)
            w%itlef = 1
            w%converge_tveg = .false.
            iterate_tveg: do while(.not.w%converge_tveg)
            
               ! Determine friction velocity, and potential temperature and humidity
               ! profiles of the surface boundary layer
               call FrictionVelocityPatch (begp, endp, p, &
                    displa(p), z0mv(p), z0hv(p), z0qv(p), &
                    obu(p), w%itlef, w%ur, um(p), w%ugust_total, ustar(p), &
                    w%temp1, w%temp2, w%temp12m, w%temp22m, w%fm, &
                    frictionvel_vars)
               
               w%tlbef = t_veg(p)
               w%del2 = w%del

               ! Determine aerodynamic resistances
               ram1(p)  = 1._r8/(ustar(p)*ustar(p)/um(p))
               w%rah(above_canopy) = 1._r8/(w%temp1*ustar(p))
               w%raw(above_canopy) = 1._r8/(w%temp2*ustar(p))

               ! Forbid removing more than 99% of wind speed in a time step.
               ! This is mainly to avoid convergence issues since this is such a
               ! basic form of iteration in this loop...
               if (implicit_stress) then
                  w%tau = forc_rho(t)*w%wind_speed_adj/ram1(p)
                  call shr_flux_update_stress(w%wind_speed0, wsresp(t), tau_est(t), &
                       w%tau, w%prev_tau, w%tau_diff, w%prev_tau_diff, &
                       w%wind_speed_adj)
                  w%ur = max(1.0_r8, sqrt(w%wind_speed_adj**2 + ugust(t)**2))
               end if

               ! Bulk boundary layer resistance of leaves

               uaf(p) = um(p)*sqrt( 1._r8/(ram1(p)*um(p)) )

               ! Use pft parameter for leaf characteristic width
               ! dleaf_patch if this is not an ed patch.
               ! Otherwise, the value has already been loaded
               ! during the FATES dynamics and/or initialization call
               if(.not.veg_pp%is_fates(p)) then
                  dleaf_patch(p) = dleaf(veg_pp%itype(p))
               end if
               
               w%cf  = 0.01_r8/(sqrt(uaf(p))*sqrt( dleaf_patch(p) ))
               
               w%rb  = 1._r8/(w%cf*uaf(p))
               rb1(p) = w%rb

               ! Parameterization for variation of csoilc with canopy density from
               ! X. Zeng, University of Arizona
               
               w%exp_ai = exp(-(elai(p)+esai(p)))
               
               ! changed by K.Sakaguchi from here
               ! transfer coefficient over bare soil is changed to a local variable
               ! just for readability of the code (from line 680)
               w%csoilb = (vkc/(0.13_r8*(z0mg(c)*uaf(p)/1.5e-5_r8)**0.45_r8))
               
               !compute the stability parameter for ricsoilc  ("S" in Sakaguchi&Zeng,2008)

               w%ri = ( grav*htop(p) * (w%taf - t_grnd(c)) ) / (w%taf * uaf(p) **2.00_r8)
               
               !! modify csoilc value (0.004) if the under-canopy is in stable condition
               
               if ( (w%taf - t_grnd(c) ) > 0._r8) then
                  ! decrease the value of csoilc by dividing it with (1+gamma*min(S, 10.0))
                  ! ria ("gmanna" in Sakaguchi&Zeng, 2008) is a constant (=0.5)
                  w%ricsoilc = csoilc / (1.00_r8 + ria*min( w%ri, 10.0_r8) )
                  w%csoilcn = w%csoilb*w%exp_ai + w%ricsoilc*(1._r8-w%exp_ai)
               else
                  w%csoilcn = w%csoilb*w%exp_ai + csoilc*(1._r8-w%exp_ai)
               end if

               !! Sakaguchi changes for stability formulation ends here

               w%rah(below_canopy) = 1._r8/(w%csoilcn*uaf(p))
               w%raw(below_canopy) = w%rah(below_canopy)
               if (use_lch4) then
                  grnd_ch4_cond(p) = 1._r8/(w%raw(above_canopy)+w%raw(below_canopy))
               end if

               ! Stomatal resistances for sunlit and shaded fractions of canopy.
               ! Done each iteration to account for differences in eah, tv.
               
               w%svpts = w%el                         ! Pa
               w%eah = forc_pbot(t) * w%qaf / mm_epsilon   ! Pa
               rhaf(p) = w%eah/w%svpts
               
               ! variables for history fields
               rah_above(p)  = w%rah(above_canopy)
               raw_above(p)  = w%raw(above_canopy)
               rah_below(p)  = w%rah(below_canopy)
               raw_below(p)  = w%raw(below_canopy)
               vpd(p)        = max((w%svpts - w%eah), vpd_min) * pa_to_kpa ! kPa
               
               ! Modification for shrubs proposed by X.D.Z
               ! Equivalent modification for soy following AgroIBIS
               ! NOTE: the following block of code was moved out of Photosynthesis subroutine and
               ! into here by M. Vertenstein on 4/6/2014 as part of making the photosynthesis
               ! routine a separate module. This move was also suggested by S. Levis in the previous
               ! version of the code.
               ! BUG MV 4/7/2014 - is this the correct place to have it in the iteration?
               ! THIS SHOULD BE MOVED OUT OF THE ITERATION but will change answers -

               if(.not.veg_pp%is_fates(p)) then
                  ! soybean (crop with N fixation)
                  if (crop(veg_pp%itype(p)) >= 1 .and. nfixer(veg_pp%itype(p)) == 1) then
                     btran(p) = min(1._r8, btran(p) * 1.25_r8)
                  end if
               end if

               w%dayl_factor =min(1._r8,max(0.01_r8,(dayl(g)*dayl(g))/(max_dayl(g)*max_dayl(g))))
               
               if(do_b4b)then
                  call WrapPhotosynthesis(bounds,p,w%svpts,w%eah,forc_po2(t),forc_pco2(t),w%rb,w%dayl_factor, &
                       btran(p),w%qsatl,w%qaf,atm2lnd_vars,canopystate_vars,photosyns_vars, &
                       soilstate_vars, surfalb_vars,solarabs_vars,cnstate_vars,energyflux_vars)
               end if
                  
               ! Sensible heat conductance for air, leaf and ground
               ! Moved the original subroutine in-line...

               w%wta    = 1._r8/w%rah(above_canopy)  ! air
               w%wtl    = (elai(p)+esai(p))/w%rb    ! leaf
               w%wtg = 1._r8/w%rah(below_canopy)  ! ground
               w%wtshi  = 1._r8/(w%wta+w%wtl+w%wtg)
               w%wtl0 = w%wtl*w%wtshi         ! leaf
               w%wtg0    = w%wtg*w%wtshi      ! ground
               w%wta0 = w%wta*w%wtshi         ! air
               
               w%wtga    = w%wta0+w%wtg0      ! ground + air
               w%wtal = w%wta0+w%wtl0   ! air + leaf
               
               ! Fraction of potential evaporation from leaf

               if (fdry(p) > 0._r8) then
                  w%rppdry  = fdry(p)*w%rb*(laisun(p)/(w%rb+rssun(p)) + &
                       laisha(p)/(w%rb+rssha(p)))/elai(p)
               else
                  w%rppdry = 0._r8
               end if

               ! Calculate canopy conductance for methane / oxygen
               ! (e.g. stomatal conductance & leaf bdy cond)
               if (use_lch4) then
                  canopy_cond(p) = (laisun(p)/(w%rb+rssun(p)) + laisha(p) / &
                       (w%rb+rssha(p)))/max(elai(p), 0.01_r8)
               end if

               w%efpot = forc_rho(t)*w%wtl*(w%qsatl-w%qaf)
               ! When the hydraulic stress parameterization is active calculate rpp
               ! but not transpiration
               if ( use_hydrstress ) then
                  if (w%efpot > 0._r8) then
                     if (btran(p) > btran0) then
                        w%rpp = w%rppdry + fwet(p)
                     else
                        w%rpp = fwet(p)
                     end if
                     !Check total evapotranspiration from leaves
                     w%rpp = min(w%rpp, (qflx_tran_veg(p)+h2ocan(p)/dtime)/w%efpot)
                  else
                     w%rpp = 1._r8
                  end if
               else
                  if (w%efpot > 0._r8) then
                     if (btran(p) > btran0) then
                        qflx_tran_veg(p) = w%efpot*w%rppdry
                        w%rpp = w%rppdry + fwet(p)
                     else
                        !No transpiration if btran below 1.e-10
                        w%rpp = fwet(p)
                        qflx_tran_veg(p) = 0._r8
                     end if
                     !Check total evapotranspiration from leaves
                     w%rpp = min(w%rpp, (qflx_tran_veg(p)+h2ocan(p)/dtime)/w%efpot)
                  else
                     !No transpiration if potential evaporation less than zero
                     w%rpp = 1._r8
                     qflx_tran_veg(p) = 0._r8
                  end if
               end if
               
               ! Update conductances for changes in rpp
               ! Latent heat conductances for ground and leaf.
               ! Air has same conductance for both sensible and latent heat.
               ! Moved the original subroutine in-line...

               w%wtaq    = frac_veg_nosno(p)/w%raw(above_canopy)             ! air
               w%wtlq    = frac_veg_nosno(p)*(elai(p)+esai(p))/w%rb * w%rpp   ! leaf

               !Litter layer resistance. Added by K.Sakaguchi
               w%snow_depth_c = z_dl ! critical depth for 100% litter burial by snow (=litter thickness)
               w%fsno_dl = snow_depth(c)/w%snow_depth_c    ! effective snow cover for (dry)plant litter
               w%elai_dl = lai_dl*(1._r8 - min(w%fsno_dl,1._r8)) ! exposed (dry)litter area index
               w%rdl = ( 1._r8 - exp(-w%elai_dl) ) / ( 0.004_r8*uaf(p)) ! dry litter layer resistance
               
               ! add litter resistance and Lee and Pielke 1992 beta
               if (w%delq < 0._r8) then  !dew. Do not apply beta for negative flux (follow old rsoil)
                  w%wtgq = frac_veg_nosno(p)/(w%raw(below_canopy)+w%rdl)
               else
                  if (do_soilevap_beta()) then
                     w%wtgq = soilbeta(c)*frac_veg_nosno(p)/(w%raw(below_canopy)+w%rdl)
                  endif
               end if
               
               w%wtsqi   = 1._r8/(w%wtaq+w%wtlq+w%wtgq)
               
               w%wtgq0    = w%wtgq*w%wtsqi      ! ground
               w%wtlq0 = w%wtlq*w%wtsqi         ! leaf
               w%wtaq0 = w%wtaq*w%wtsqi         ! air
               
               w%wtgaq    = w%wtaq0+w%wtgq0     ! air + ground
               w%wtalq = w%wtaq0+w%wtlq0  ! air + leaf
               
               w%dc1 = forc_rho(t)*cpair*w%wtl
               w%dc2 = hvap*forc_rho(t)*w%wtlq
               
               w%efsh   = w%dc1*(w%wtga*t_veg(p)-w%wtg0*t_grnd(c)-w%wta0*thm(p))
               w%efe = w%dc2*(w%wtgaq*w%qsatl-w%wtgq0*qg(c)-w%wtaq0*forc_q(t))

               ! Evaporation flux from foliage
               w%erre = 0._r8
               if (w%efe*w%efeb < 0._r8) then
                  w%efeold = w%efe
                  w%efe  = 0.1_r8*w%efeold
                  w%erre = w%efe - w%efeold
               end if

               ! fractionate ground emitted longwave
               w%lw_grnd=(frac_sno(c)*t_soisno(c,snl(c)+1)**4 &
                    +(1._r8-frac_sno(c)-frac_h2osfc(c))*t_soisno(c,1)**4 &
                    +frac_h2osfc(c)*t_h2osfc(c)**4)
               
               dt_veg(p) = (sabv(p) + w%air + w%bir*t_veg(p)**4 + &
                    w%cir*w%lw_grnd - w%efsh - w%efe) / &
                    (- 4._r8*w%bir*t_veg(p)**3 +w%dc1*w%wtga +w%dc2*w%wtgaq*w%qsatldT)
               t_veg(p) = w%tlbef + dt_veg(p)
               w%dels = dt_veg(p)
               w%del  = abs(w%dels)
               w%err = 0._r8
               if (w%del > delmax) then
                  dt_veg(p) = delmax*w%dels/w%del
                  t_veg(p) = w%tlbef + dt_veg(p)
                  w%err = sabv(p) + w%air + w%bir*w%tlbef**3*(w%tlbef + &
                       4._r8*dt_veg(p)) + w%cir*w%lw_grnd - &
                       (w%efsh + w%dc1*w%wtga*dt_veg(p)) - (w%efe + &
                       w%dc2*w%wtgaq*w%qsatldT*dt_veg(p))
               end if
               
               ! Fluxes from leaves to canopy space
               ! "efe" was limited as its sign changes frequently.  This limit may
               ! result in an imbalance in "hvap*qflx_evap_veg" and
               ! "efe + dc2*wtgaq*qsatdt_veg"
               
               w%efpot = forc_rho(t)*w%wtl*(w%wtgaq*(w%qsatl+w%qsatldT*dt_veg(p)) &
                    -w%wtgq0*qg(c)-w%wtaq0*forc_q(t))
               qflx_evap_veg(p) = w%rpp*w%efpot
               
               ! Calculation of evaporative potentials (efpot) and
               ! interception losses; flux in kg m**-2 s-1.  ecidif
               ! holds the excess energy if all intercepted water is evaporated
               ! during the timestep.  This energy is later added to the
               ! sensible heat flux.
               if ( use_hydrstress ) then
                  w%ecidif = max(0._r8,qflx_evap_veg(p)-qflx_tran_veg(p)-h2ocan(p)/dtime)
                  qflx_evap_veg(p) = min(qflx_evap_veg(p),qflx_tran_veg(p)+h2ocan(p)/dtime)
               else
                  
                  w%ecidif = 0._r8
                  if (w%efpot > 0._r8 .and. btran(p) > btran0) then
                     qflx_tran_veg(p) = w%efpot*w%rppdry
                  else
                     qflx_tran_veg(p) = 0._r8
                  end if
                  w%ecidif = max(0._r8, qflx_evap_veg(p)-qflx_tran_veg(p)-h2ocan(p)/dtime)
                  qflx_evap_veg(p) = min(qflx_evap_veg(p),qflx_tran_veg(p)+h2ocan(p)/dtime)
               end if

               ! The energy loss due to above two limits is added to
               ! the sensible heat flux.
               eflx_sh_veg(p) = w%efsh + w%dc1*w%wtga*dt_veg(p) + w%err + w%erre + hvap*w%ecidif
               
               ! Re-calculate saturated vapor pressure, specific humidity, and their
               ! derivatives at the leaf surface
               
               call QSat(t_veg(p), forc_pbot(t), w%el, w%deldT, w%qsatl, w%qsatldT)
               
               ! Update vegetation/ground surface temperature, canopy air
               ! temperature, canopy vapor pressure, aerodynamic temperature, and
               ! Monin-Obukhov stability parameter for next iteration.
               
               w%taf = w%wtg0*t_grnd(c) + w%wta0*thm(p) + w%wtl0*t_veg(p)
               w%qaf = w%wtlq0*w%qsatl + w%wtgq0*qg(c) + forc_q(t)*w%wtaq0
               
               ! Update Obukhov length scale and wind speed including the
               ! stability effect
               
               w%dth = thm(p)-w%taf
               w%dqh = forc_q(t)-w%qaf
               w%delq = w%wtalq*qg(c)-w%wtlq0*w%qsatl-w%wtaq0*forc_q(t)
               
               w%tstar = w%temp1*w%dth
               w%qstar = w%temp2*w%dqh
               
               w%thvstar = w%tstar*(1._r8+0.61_r8*forc_q(t)) + 0.61_r8*forc_th(t)*w%qstar
               zeta(p) = zldis(p)*vkc*grav*w%thvstar/(ustar(p)**2*thv(c))
               
               if (zeta(p) >= 0._r8) then     !stable
                  zeta(p) = min(2._r8,max(zeta(p),0.01_r8))
                  um(p) = max(w%ur,0.1_r8)
               else                     !unstable
                  zeta(p) = max(-100._r8,min(zeta(p),-0.01_r8))
                  if ((.not. atm_gustiness) .or. force_land_gustiness) then
                     w%wc = beta*(-grav*ustar(p)*w%thvstar*zii/thv(c))**0.333_r8
                     w%ugust_total = sqrt(ugust(t)**2 + w%wc**2)
                     um(p) = sqrt(w%ur*w%ur+w%wc*w%wc)
                  else
                     um(p) = max(w%ur,0.1_r8)
                  end if
               end if
               obu(p) = zldis(p)/zeta(p)

               if (w%obuold*obu(p) < 0._r8) w%nmozsgn = w%nmozsgn+1
               if (w%nmozsgn >= 4) obu(p) = zldis(p)/(-0.01_r8)
               w%obuold = obu(p)

               ! laminar boundary resistance for h2o over leaf,
               ! should I make this consistent for latent heat calculation?
               lbl_rsc_h2o(p) = getlblcef(forc_rho(t),t_veg(p))*uaf(p)/(uaf(p)**2._r8+1.e-10_r8)   

               ! Test for convergence
               w%iter_final = w%itlef
               w%itlef = w%itlef+1
               
               w%dele = abs(w%efe-w%efeb)
               w%efeb = w%efe
               w%det  = max(w%del,w%del2)
               num_iter(p) = real(w%itlef,r8)
               
               if ( (.not. (w%det < dtmin .and. w%dele < dlemin) .or. &
                    (implicit_stress .and. abs(w%tau_diff) >= dtaumin)) .and. &
                    (w%itlef < itmax)) then
                  w%converge_tveg = .false.
               else
                  w%converge_tveg = .true.
               end if
            end do iterate_tveg

            ! Evaluate quality of conductance solution
            !
            ! Criteria for finding a solution to the outer loop
            !
            ! 1) Always make sure that at least one photosynthesis call
            !    is made. (ie itstoma>0)
            ! 2) Calculate the change in resistance that was made on the
            !    last solution. If the difference is negligable, and
            !    condition 1 is satisfied, then you have a solution
            ! 3) Exit if too many attempts and accept what you have
            !    (ie. itstoma>itmax_stomata

            !reldel_rs = 2._r8*max( abs(rssun(p)-rssun_old(p))/(rssun(p)+rssun_old(p)), &
            !     abs(rssha(p)-rssha_old(p))/(rssha(p)+rssha_old(p)) )
            
            w%del_gs = max( abs(1._r8/rssun(p)-1._r8/w%rssun_old), &
                          abs(1._r8/rssha(p)-1._r8/w%rssha_old) )

            istoma_converge_if: if( do_b4b .or. &
                                    (w%del_gs < max_del_gs ) .or.  &
                                    (w%itstoma>=itmax_stomata) ) then
               w%converge_stoma = .true.
               
            else
               
               ! Update the outer (stomata c) counter
               w%itstoma = w%itstoma + 1
               num_iter(p) = num_iter(p) + 1
               
               ! Update the previous resistances
               w%rssun_old = rssun(p)
               w%rssha_old = rssha(p)

               ! Call photosynthesis and retrieve
               ! updated stomatal conductances
               
               ! Instead of updating stomatal conductances
               ! we hold the value calculated in the outer loop
               ! as constant during this inner loop (tveg) iteration

               call WrapPhotosynthesis(bounds,p,w%svpts,w%eah,forc_po2(t),forc_pco2(t),w%rb,w%dayl_factor, &
                    btran(p),w%qsatl,w%qaf,atm2lnd_vars,canopystate_vars,photosyns_vars, &
                    soilstate_vars, surfalb_vars,solarabs_vars,cnstate_vars,energyflux_vars)
               
            end if istoma_converge_if

         end do iterate_stoma
         

         ! Energy balance check in canopy

         w%lw_grnd=(frac_sno(c)*t_soisno(c,snl(c)+1)**4 &
              +(1._r8-frac_sno(c)-frac_h2osfc(c))*t_soisno(c,1)**4 &
              +frac_h2osfc(c)*t_h2osfc(c)**4)

         w%err = sabv(p) + w%air + w%bir*w%tlbef**3*(w%tlbef + 4._r8*dt_veg(p)) &
                                !+ cir(p)*t_grnd(c)**4 - eflx_sh_veg(p) - hvap*qflx_evap_veg(p)
              + w%cir*w%lw_grnd - eflx_sh_veg(p) - hvap*qflx_evap_veg(p)

         ! Fluxes from ground to canopy space

         w%delt    = w%wtal*t_grnd(c)-w%wtl0*t_veg(p)-w%wta0*thm(p)
         taux(p) = -forc_rho(t)*forc_u(t)/ram1(p)
         tauy(p) = -forc_rho(t)*forc_v(t)/ram1(p)
         if (implicit_stress) then
            taux(p) = taux(p) * (w%wind_speed_adj / w%wind_speed0)
            tauy(p) = tauy(p) * (w%wind_speed_adj / w%wind_speed0)
         end if
         eflx_sh_grnd(p) = cpair*forc_rho(t)*w%wtg*w%delt

         ! compute individual sensible heat fluxes
         w%delt_snow = w%wtal*t_soisno(c,snl(c)+1)-w%wtl0*t_veg(p)-w%wta0*thm(p)
         eflx_sh_snow(p) = cpair*forc_rho(t)*w%wtg*w%delt_snow

         w%delt_soil  = w%wtal*t_soisno(c,1)-w%wtl0*t_veg(p)-w%wta0*thm(p)
         eflx_sh_soil(p) = cpair*forc_rho(t)*w%wtg*w%delt_soil

         w%delt_h2osfc  = w%wtal*t_h2osfc(c)-w%wtl0*t_veg(p)-w%wta0*thm(p)
         eflx_sh_h2osfc(p) = cpair*forc_rho(t)*w%wtg*w%delt_h2osfc
         qflx_evap_soi(p) = forc_rho(t)*w%wtgq*w%delq

         ! compute individual latent heat fluxes
         w%delq_snow = w%wtalq*qg_snow(c)-w%wtlq0*w%qsatl-w%wtaq0*forc_q(t)
         qflx_ev_snow(p) = forc_rho(t)*w%wtgq*w%delq_snow

         w%delq_soil = w%wtalq*qg_soil(c)-w%wtlq0*w%qsatl-w%wtaq0*forc_q(t)
         qflx_ev_soil(p) = forc_rho(t)*w%wtgq*w%delq_soil

         w%delq_h2osfc = w%wtalq*qg_h2osfc(c)-w%wtlq0*w%qsatl-w%wtaq0*forc_q(t)
         qflx_ev_h2osfc(p) = forc_rho(t)*w%wtgq*w%delq_h2osfc

         ! 2 m height air temperature

         t_ref2m(p) = thm(p) + w%temp1*w%dth*(1._r8/w%temp12m - 1._r8/w%temp1)
         t_ref2m_r(p) = t_ref2m(p)

         ! 2 m height specific humidity

         q_ref2m(p) = forc_q(t) + w%temp2*w%dqh*(1._r8/w%temp22m - 1._r8/w%temp2)

         ! 2 m height relative humidity

         call QSat(t_ref2m(p), forc_pbot(t), w%e_ref2m, w%de2mdT, w%qsat_ref2m, w%dqsat2mdT)
         rh_ref2m(p) = min(100._r8, q_ref2m(p) / w%qsat_ref2m * 100._r8)
         rh_ref2m_r(p) = rh_ref2m(p)

         ! Downward longwave radiation below the canopy

         dlrad(p) = (1._r8-emv(p))*emg(c)*forc_lwrad(t) + &
              emv(p)*emg(c)*sb*w%tlbef**3*(w%tlbef + 4._r8*dt_veg(p))

         ! Upward longwave radiation above the canopy

         ulrad(p) = ((1._r8-emg(c))*(1._r8-emv(p))*(1._r8-emv(p))*forc_lwrad(t) &
              + emv(p)*(1._r8+(1._r8-emg(c))*(1._r8-emv(p)))*sb*w%tlbef**3*(w%tlbef + &
              4._r8*dt_veg(p)) + emg(c)*(1._r8-emv(p))*sb*w%lw_grnd)

         ! Derivative of soil energy flux with respect to soil temperature

         cgrnds(p) = cgrnds(p) + cpair*forc_rho(t)*w%wtg*w%wtal
         cgrndl(p) = cgrndl(p) + forc_rho(t)*w%wtgq*w%wtalq*dqgdT(c)
         cgrnd(p)  = cgrnds(p) + cgrndl(p)*htvp(c)

         ! Update dew accumulation (kg/m2)

         h2ocan(p) = max(0._r8,h2ocan(p)+(qflx_tran_veg(p)-qflx_evap_veg(p))*dtime)

         ! Check for convergence of stress.
         if (implicit_stress .and. abs(w%tau_diff) > dtaumin) then
            if (nstep_mod > 0) then ! Suppress common warnings on the first time step.
               write(iulog,*)'WARNING: Stress did not converge for canopy ',&
                    ' nstep = ',nstep_mod,' p= ',p,' prev_tau_diff= ',w%prev_tau_diff,&
                    ' tau_diff= ',w%tau_diff,' tau= ',w%tau,&
                    ' wind_speed_adj= ',w%wind_speed_adj,' iter_final= ',w%iter_final
            end if
         end if

         ! Copy back values needed after loop
         err(p) = w%err
         zldis(p) = w%zldis

      end do do_patch
      !$OMP END PARALLEL DO

      ! Check if forcing height is below canopy height for any patch
      if ( .not. use_fates ) then
         do f = 1, fn
            p = filterp(f)
            if (zldis(p) < 0._r8) then
#ifndef _OPENACC
               write(iulog,*)'Error: Forcing height is below canopy height for pft index ',p
               call endrun(decomp_index=index, elmlevel=namep, msg=errmsg(__FILE__, __LINE__))
#endif
            end if
         end do
      end if
      
      call t_stop_lnd(event)

      if ( use_fates ) then

#ifndef _OPENACC
        call alm_fates%WrapAccumulateFluxes(bounds,fn,filterp(1:fn))
        call alm_fates%wrap_hydraulics_drive(bounds,fn,filterp(1:fn),soilstate_vars, &
             solarabs_vars,energyflux_vars)
#endif
      else

         ! Determine total photosynthesis
         call PhotosynthesisTotal(fn, filterp, &
              atm2lnd_vars, cnstate_vars, canopystate_vars, photosyns_vars)
      end if
      
      ! Report high energy balance errors
      do f = 1, fn
         p = filterp(f)
         if (abs(err(p)) > 0.1_r8) then
            write(iulog,*) 'energy balance in canopy ',p,', err=',err(p)
         end if
      end do

    end associate
  end subroutine CanopyFluxes

  ! =========================================================================
  
  subroutine WrapPhotosynthesis(bounds,p,svpts,eah,o2,co2,rb,dayl_factor, &
       btran,qsatl,qaf,atm2lnd_vars,canopystate_vars,photosyns_vars, &
       soilstate_vars, surfalb_vars,solarabs_vars,cnstate_vars,energyflux_vars)

    
    type(bounds_type)         , intent(in)    :: bounds
    integer  :: p                                        ! patch index from begp:endp
    real(r8) :: svpts
    real(r8) :: eah
    real(r8) :: o2
    real(r8) :: co2
    real(r8) :: rb
    real(r8) :: dayl_factor
    real(r8) :: btran
    real(r8) :: qsatl
    real(r8) :: qaf
    
    type(atm2lnd_type)        , intent(inout) :: atm2lnd_vars
    type(canopystate_type)    , intent(inout) :: canopystate_vars
    type(photosyns_type)      , intent(inout) :: photosyns_vars
    type(soilstate_type)      , intent(inout) :: soilstate_vars
    type(surfalb_type)        , intent(inout) :: surfalb_vars
    type(solarabs_type)       , intent(inout) :: solarabs_vars
    type(cnstate_type)        , intent(inout) :: cnstate_vars
    type(energyflux_type)     , intent(inout) :: energyflux_vars

    real(r8) :: svpts_a(bounds%begp:bounds%endp)
    real(r8) :: eah_a(bounds%begp:bounds%endp)
    real(r8) :: o2_a(bounds%begp:bounds%endp)
    real(r8) :: co2_a(bounds%begp:bounds%endp)
    real(r8) :: rb_a(bounds%begp:bounds%endp)
    real(r8) :: dayl_factor_a(bounds%begp:bounds%endp)
    real(r8) :: btran_a(bounds%begp:bounds%endp)
    real(r8) :: qsatl_a(bounds%begp:bounds%endp)
    real(r8) :: qaf_a(bounds%begp:bounds%endp)
    
    integer :: begp,endp

    begp = bounds%begp
    endp = bounds%endp
    
    ! Instead of updating stomatal conductances
    ! we hold the value calculated in the outer loop
    ! as constant during this inner loop (tveg) iteration
    if_fates: if ( use_fates ) then
#ifndef _OPENACC

       call alm_fates%WrapPatchPhotosynthesis(bounds, p, &
            svpts, eah, o2, &
            co2, rb, dayl_factor, &
            atm2lnd_vars, canopystate_vars, photosyns_vars)
#endif
    else ! not use_fates

       svpts_a(p) = svpts
       eah_a(p) = eah
       o2_a(p) = o2
       co2_a(p) = co2
       rb_a(p) = rb
       dayl_factor_a(p) = dayl_factor
       btran_a(p) = btran
       qsatl_a(p) = qsatl
       qaf_a(p) = qaf
       
       if ( use_hydrstress ) then
          call PhotosynthesisHydraulicStress (bounds, 1, [p], &
               svpts_a(begp:endp), eah_a, o2_a, co2_a, rb_a, &
               energyflux_vars%bsun_patch(begp:endp), energyflux_vars%bsha_patch(begp:endp), &
               btran_a, dayl_factor_a, &
               qsatl_a, qaf_a,     &
               atm2lnd_vars, soilstate_vars, surfalb_vars, solarabs_vars,    &
               canopystate_vars, photosyns_vars)
       else
          call Photosynthesis (bounds, 1, [p], &
               svpts_a, eah_a, o2_a, co2_a, rb_a, btran_a, &
               dayl_factor_a, atm2lnd_vars,  surfalb_vars, solarabs_vars, &
               canopystate_vars, photosyns_vars, 'sun')
       end if
       
       if ( use_c13 ) then
          call Fractionation (bounds, 1, [p], &
               cnstate_vars, solarabs_vars, surfalb_vars, photosyns_vars, 1)
       endif
       
       ! soybean (crop with N fixation)
       if (crop(veg_pp%itype(p)) >= 1 .and. nfixer(veg_pp%itype(p)) == 1) then
          btran_a(p) = min(1._r8, btran_a(p) * 1.25_r8)
       end if
       
       if ( .not. use_hydrstress ) then
          call Photosynthesis (bounds, 1, [p], &
               svpts_a, eah_a, o2_a, co2_a, rb_a, btran_a, &
               dayl_factor_a, atm2lnd_vars,surfalb_vars, solarabs_vars, &
               canopystate_vars, photosyns_vars, 'sha')
      end if
      
       if ( use_c13 ) then
          call Fractionation (bounds, 1, [p],  &
               cnstate_vars, solarabs_vars, surfalb_vars, photosyns_vars, 0)
       end if
    end if if_fates ! end of if use_fates

  end subroutine WrapPhotosynthesis

  ! ===========================================================================
  
  subroutine FrictionVelocityPatch(lbn,ubn, p, displa, z0m, z0h, z0q, &
       obu, iter, ur, um, ugust, ustar, &
       temp1, temp2, temp12m, temp22m, fm, &
       frictionvel_vars)

    implicit none
    !
    ! !ARGUMENTS:
    integer  , intent(in)  :: lbn
    integer  , intent(in)  :: ubn
    integer  , intent(in)  :: p
    real(r8) , intent(in)  :: displa
    real(r8) , intent(in)  :: z0m
    real(r8) , intent(in)  :: z0h
    real(r8) , intent(in)  :: z0q
    real(r8) , intent(in)  :: obu
    integer  , intent(in)  :: iter
    real(r8) , intent(in)  :: ur
    real(r8) , intent(in)  :: um
    real(r8) , intent(in)  :: ugust
    real(r8) , intent(out) :: ustar
    real(r8) , intent(out) :: temp1
    real(r8) , intent(out) :: temp12m
    real(r8) , intent(out) :: temp2
    real(r8) , intent(out) :: temp22m
    real(r8) , intent(out) :: fm
    type(frictionvel_type) , intent(inout) :: frictionvel_vars


    real(r8) :: displa_a( lbn:ubn ) 
    real(r8) :: z0m_a( lbn:ubn ) 
    real(r8) :: z0h_a( lbn:ubn ) 
    real(r8) :: z0q_a( lbn:ubn ) 
    real(r8) :: obu_a( lbn:ubn ) 
    real(r8) :: ur_a( lbn:ubn ) 
    real(r8) :: um_a( lbn:ubn ) 
    real(r8) :: ugust_a( lbn:ubn ) 
    real(r8) :: ustar_a   ( lbn:ubn )         ! friction velocity [m/s] [lbn:ubn]
    real(r8) :: temp1_a   ( lbn:ubn )         ! relation for potential temperature profile [lbn:ubn]
    real(r8) :: temp12m_a ( lbn:ubn )         ! relation for potential temperature profile applied at 2-m [lbn:ubn]
    real(r8) :: temp2_a   ( lbn:ubn )         ! relation for specific humidity profile [lbn:ubn]
    real(r8) :: temp22m_a ( lbn:ubn )         ! relation for specific humidity profile applied at 2-m [lbn:ubn]
    real(r8) :: fm_a      ( lbn:ubn )         ! diagnose 10m wind (DUST only) [lbn:ubn]


    displa_a(p) = displa
    z0m_a(p)  = z0m
    z0h_a(p)  = z0h
    z0q_a(p)  = z0q
    obu_a(p)  = obu
    ur_a(p)   = ur
    um_a(p)   = um
    ugust_a(p) = ugust
    
    call FrictionVelocity (lbn, ubn, 1, [p], &
         displa_a, z0m_a, z0h_a, z0q_a, &
         obu_a, iter, ur_a, um_a, ugust_a, ustar_a, &
         temp1_a, temp2_a, temp12m_a, temp22m_a, fm_a, &
         frictionvel_vars)

    ustar   = ustar_a(p)
    temp1   = temp1_a(p)
    temp2   = temp2_a(p)
    temp12m = temp12m_a(p)
    temp22m = temp22m_a(p)
    fm      = fm_a(p)

    return
  end subroutine FrictionVelocityPatch

end module CanopyFluxesMod
