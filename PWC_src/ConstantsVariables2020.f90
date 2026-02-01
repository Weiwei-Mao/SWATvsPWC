MODULE constants_and_variables

    IMPLICIT NONE
    SAVE

    ! cummulatives for research:
    REAL :: hold_precip
    REAL :: hold_irrig
    REAL :: hold_runoff

    INTEGER :: xx= 89  ! temporary dummy get rid of this

    REAL, PARAMETER :: Version_Number = 1.003 !version 1

    LOGICAL :: use_bidiagonal
    ! tracks the appication widow for placement into a storage array for median calculation of concentrations
    INTEGER :: app_window_counter

    INTEGER, PARAMETER :: number_medians = 11

    REAL :: hold_for_medians_wb(366, number_medians)   ! holds app window values for median determination
    REAL :: hold_for_medians_WPEZ(366,number_medians)  ! holds app window values for median determination
    REAL :: hold_for_medians_TPEZ(366,1)               ! Holds TPEZ acute value, make as a 2d array for the smedian subroutine call

    REAL :: hold_for_medians_daughter(366,number_medians)      ! holds app window values for median determination
    REAL :: hold_for_medians_WPEZ_daughter(366,number_medians) ! holds app window values for median determination
    REAL :: hold_for_medians_TPEZ_daughter(366,1)              ! Holds TPEZ acute value

    REAL :: hold_for_medians_grandaughter(366,number_medians)       ! holds app window values for median determination
    REAL :: hold_for_medians_WPEZ_grandaughter(366,number_medians)  ! holds app window values for median determination
    REAL :: hold_for_medians_TPEZ_grandaughter(366,1)               ! Holds TPEZ acute value

    ! REAL :: medians_conc(number_medians)  !keeps final medians for main (current) waterbody simulations

    ! ********************************************************************************************************************
    ! File Names
    CHARACTER (LEN=512) PRZMVVWMinputfile
    CHARACTER (LEN=512) inputfile ! new input file to replace przmvvwminputfile. this file is the same as the pwc input file

    ! ********************************************************************************************************************
    ! File Unit Numbers
    INTEGER, PARAMETER ::  waterbody_timeseries_unit  = 12
    INTEGER, PARAMETER ::  wpez_timeseries_unit       = 13

    INTEGER, PARAMETER ::  summary_output_unit      = 44
    INTEGER, PARAMETER ::  summary_output_unit_deg1 = 45
    INTEGER, PARAMETER ::  summary_output_unit_deg2 = 46

    INTEGER, PARAMETER ::  summary_output_unit_tpez      = 47
    INTEGER, PARAMETER ::  summary_output_unit_tpez_deg1 = 48
    INTEGER, PARAMETER ::  summary_output_unit_tpez_deg2 = 49


    INTEGER, PARAMETER ::  summary_wpez_unit      = 50
    INTEGER, PARAMETER ::  summary_wpez_unit_deg1 = 51
    INTEGER, PARAMETER ::  summary_wpez_unit_deg2 = 52

    INTEGER, PARAMETER ::  inputfile_unit_number  = 53
    INTEGER, PARAMETER ::  PRZMinputUnit          = 54
    INTEGER, PARAMETER ::  MetFileUnit            = 55
    INTEGER, PARAMETER ::  ScenarioFileUnit       = 56
    INTEGER, PARAMETER ::  BatchFileUnit          = 57

    INTEGER, PARAMETER ::  TimeSeriesUnit2         = 58
    INTEGER, PARAMETER ::  waterbody_file_unit     = 59

    INTEGER, PARAMETER ::  median_unit                   = 60
    INTEGER, PARAMETER ::  median_daughter_unit          = 61
    INTEGER, PARAMETER ::  median_grandaughter_unit      = 62

    INTEGER, PARAMETER ::  median_unit_wpez              = 63
    INTEGER, PARAMETER ::  median_daughter_unit_wpez     = 64
    INTEGER, PARAMETER ::  median_grandaughter_unit_wpez = 65

    INTEGER, PARAMETER ::  median_unit_tpez              = 66
    INTEGER, PARAMETER ::  median_daughter_unit_tpez     = 67
    INTEGER, PARAMETER ::  median_grandaughter_unit_tpez = 68


    CHARACTER(LEN = 500), PARAMETER :: summary_outputfile      =  "summary.txt"
    CHARACTER(LEN = 500), PARAMETER :: summary_outputfile_deg1 =  "summary_Deg1.txt"
    CHARACTER(LEN = 500), PARAMETER :: summary_outputfile_deg2 =  "summary_Deg2.txt"

    CHARACTER(LEN = 500), PARAMETER :: summary_outputfile_tpez      =  "summary_tpez.txt"
    CHARACTER(LEN = 500), PARAMETER :: summary_outputfile_tpez_deg1 =  "summary_tpez_deg1.txt"
    CHARACTER(LEN = 500), PARAMETER :: summary_outputfile_tpez_deg2 =  "summary_tpez_deg2.txt"

    CHARACTER(LEN = 500), PARAMETER :: summary_WPEZoutputfile      =  "summary_WPEZ.txt"
    CHARACTER(LEN = 500), PARAMETER :: summary_WPEZoutputfile_deg1 =  "summary_WPEZ_deg1.txt"
    CHARACTER(LEN = 500), PARAMETER :: summary_WPEZoutputfile_deg2 =  "summary_WPEZ_deg2.txt"


    ! scenario_unit_number :: 80 local number in READ scenario routine

    ! ********************************************************************************************************************
    INTEGER :: number_of_schemes
    ! CHARACTER (LEN=512) scheme_name
    ! INTEGER :: scheme_number

    ! ! Maximum nuber of crop periods that the scenario file holds, Currently the PWC 2020 it is 7 periods oer year
    INTEGER, PARAMETER :: max_number_crop_periods  = 7

    CHARACTER(LEN=10)  :: waterbodytext  ! text for water body type (not an input)

    CHARACTER(LEN= 256) :: run_id
    CHARACTER(LEN= 512) :: full_run_identification  ! run id plus all path information

    ! GW --------------------------------------

    ! output
    REAL :: effective_washout, effective_watercol_metab, effective_hydrolysis, effective_photolysis, &
        effective_volatization, effective_total_deg1, effective_burial, effective_benthic_metab, &
        effective_benthic_hydrolysis, effective_total_deg2


    ! Scheme inputs
    INTEGER, ALLOCATABLE, DIMENSION(:) :: number_of_scenarios
    INTEGER, ALLOCATABLE, DIMENSION(:) :: app_reference_point_schemes  ! 0=absolute date, 1=emergence, 2=maturity, 3=removal
    INTEGER, ALLOCATABLE, DIMENSION(:) :: num_apps_in_schemes

    LOGICAL, ALLOCATABLE, DIMENSION(:) :: is_adjust_for_rain_schemes
    REAL,    ALLOCATABLE, DIMENSION(:) :: rain_limit_schemes
    INTEGER, ALLOCATABLE, DIMENSION(:) :: optimum_application_window_schemes
    INTEGER, ALLOCATABLE, DIMENSION(:) :: intolerable_rain_window_schemes
    INTEGER, ALLOCATABLE, DIMENSION(:) ::	min_days_between_apps_schemes

    REAL,    ALLOCATABLE, DIMENSION(:) :: runoff_mitigation_schemes  ! scheme array of mass reducers for runoff
    REAL,    ALLOCATABLE, DIMENSION(:) :: erosion_mitigation_schemes
    REAL,    ALLOCATABLE, DIMENSION(:) :: drift_mitigation_schemes
    REAL                               :: runoff_mitigation          ! these are scheme-level proportional mass reductions
    REAL                               :: erosion_mitigation
    REAL                               :: drift_mitigation


    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: method_schemes
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: days_until_applied_schemes   ! (scheme #, application # max of 366)
    LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: is_absolute_year_schemes     ! are applications put on during specific years?
    LOGICAL, ALLOCATABLE, DIMENSION(:)   :: is_absolute_year

    ! CHARACTER(LEN=10), ALLOCATABLE,DIMENSION(:,:) absolute_date_option
    ! this is an INTEGER corresponding to the drift in the waterbody file, inexed: scheme number, app number
    INTEGER,ALLOCATABLE,DIMENSION(:,:)      :: drift_schemes
    !distance (ft) of spray buffers or whatever distance unit the spray table is in
    REAL,ALLOCATABLE,DIMENSION(:,:)         :: driftfactor_schemes

    REAL,ALLOCATABLE,DIMENSION(:,:)  :: application_rate_schemes
    REAL,ALLOCATABLE,DIMENSION(:,:)  :: depth_schemes
    REAL,ALLOCATABLE,DIMENSION(:,:)  :: split_schemes

    INTEGER,ALLOCATABLE,DIMENSION(:,:) :: lag_schemes
    INTEGER,ALLOCATABLE,DIMENSION(:,:) :: periodicity_schemes

    INTEGER :: app_reference_point   ! 0=absolute date, 1=emergence, 2=maturity, 3=removal;   scalar sent to przm-vvwm

    CHARACTER (LEN=512), ALLOCATABLE, DIMENSION(:,:) :: scenario_names      ! (i,j) i is scheme #, j is scenario #

    CHARACTER (LEN=512), ALLOCATABLE, DIMENSION(:)   :: scenario_batchfile  ! name of batch scenario file if used
    LOGICAL, ALLOCATABLE, DIMENSION(:)               :: is_batch_scenario


    ! ********************************************************************************************************************
    LOGICAL :: is_koc
    LOGICAL :: is_needs_poundkg_conversion

    ! The following are new additionas and need to be populated and converted as appropriate

    ! molar conversions daughter granddaughter
    REAL :: xAerobic(2)
    REAL :: xBenthic(2)
    REAL :: xPhoto(2)
    REAL ::  xHydro(2)

    CHARACTER (LEN=256) :: outputfiLEName = "TestThis"      !output file name

    !REAL:: D_over_dx
    !REAL:: benthic_depth
    !REAL:: porosity
    !REAL:: bulk_density  ! dry mass/total volume  g/ml
    !REAL:: FROC2         ! fraction oc in benthic sed
    !REAL:: DOC2
    !REAL:: BNMAS
    !REAL:: DFAC
    !REAL:: SUSED   !mg/L
    !REAL:: CHL

    !REAL:: FROC1
    !REAL:: DOC1
    !REAL:: PLMAS
    !REAL :: area_waterbody ! water body area (m2)
    !REAL :: depth_0        ! initial water body depth (m)
    !REAL :: depth_max      ! maximum water body depth (m)
    !Ldays used to average flow for special case reservoir, 0 indicates complete sim average
    !INTEGER :: flow_averaging
    !REAL :: baseflow

    REAL :: cropped_fraction ! fractional area used for crop: only used here to display in output
    REAL, PARAMETER :: CLOUD = 0.

    REAL, PARAMETER :: wind_height= 6.     !height at which wind measurements are reported (m)
    REAL, PARAMETER :: DELT_vvwm = 86400.  !seconds, simulation TIME INTERVAL IS ONE DAY. used for vvwm only
    REAL, PARAMETER :: minimum_depth = 0.00001     !minimum water body depth for vvwm to operate stability wise

    LOGICAL :: is_hed_files_made
    LOGICAL :: is_add_return_frequency
    REAL    :: additional_return_frequency
    LOGICAL :: ConstVolnoflow_type, ConstVolflow_type, pond_type, reservoir_type
    ! LOGICAL ::vvwm_type,


    CHARACTER (LEN = 10) waterbody_id

    ! ********************************************************************************************************************
    LOGICAL :: is_irrigation
    LOGICAL :: is_allyear_irrigation
    LOGICAL :: is_above_crop_irrigation

    ! Time series trasfer from PRZM to VVWM
    ! REAL, ALLOCATABLE, DIMENSION(:) :: runoff_series   !runoff in cm
    REAL, ALLOCATABLE, DIMENSION(:) :: erosion_save  !eroded solids in (metric tons)

    REAL, ALLOCATABLE, DIMENSION(:) :: runoff_save   !runoff in cm
    !REAL, ALLOCATABLE, DIMENSION(:) :: erosion_save

    REAL, ALLOCATABLE, DIMENSION(:) :: et_save

    REAL, ALLOCATABLE, DIMENSION(:,:) :: conc_last_horizon_save


    ! ********************************************************************************************************************
    ! PARAMETERs used for calibration Routines
    INTEGER :: data_number = 0
    INTEGER :: data_date(100) = 0
    REAL    :: data_erosion(100) = 0.0
    REAL    :: data_runoff(100) = 0.0
    REAL    :: data_mgsed(100) = 0.0
    REAL    :: data_mgh2o(100) = 0.0

    REAL    :: model_erosion(100) = 0.0
    REAL    :: model_runoff(100) = -0.0
    REAL    :: model_mg_ro(100) = 0.0
    REAL    :: model_mg_er(100) = 0.0

    REAL    :: cn_moisture(100) = 0.0

    ! ********************************************************************************************************************
    ! LOGICAL :: is_TC_lag_method

    INTEGER :: day_number_chemtrans  !track days for przm chem transport loop

    INTEGER, PARAMETER :: maxFileLENgth = 300  ! Maximum path and File name LENgth  probably should limit to 256 due to windows limit
    INTEGER, PARAMETER :: max_horizons = 25    ! maximum number of soil horizons
    INTEGER, PARAMETER :: max_number_plots = 100

    CHARACTER (LEN = maxFileLENgth) :: working_directory
    CHARACTER (LEN = maxFileLENgth) :: family_name
    CHARACTER (LEN = maxFileLENgth) :: weatherFile_directory


    !Scenario PARAMETERs
    CHARACTER (LEN = 100)             :: scenario_id  ! used in batch scenario file
    CHARACTER (LEN = maxFileLENgth)   :: weatherfiLEName
    CHARACTER (LEN = maxFileLENgth)   :: weatherfiledirectory
    REAL                              :: latitude


    INTEGER, PARAMETER :: Num_hydro_factors = 100
    ! 1 day time step, this value is not adjustable unless considerations are made for runoff which is daily
    REAL, PARAMETER    :: DELT = 1.0
    ! number of sub time steps in tridigonal calculation (used to investigate effect on nonlinear isotherms)
    INTEGER     :: number_subdelt
    REAL        :: subdelt         !delt/number_subdelt

    INTEGER :: NCOM2               ! actual number of total compartments
    ! a user specified node that is used for specific output of flux at a depth (node), not hooked up yet
    INTEGER :: user_node =1


    ! ********************************************************************************************************************
    ! Soil Profile Dependent Input PARAMETERs
    INTEGER :: NHORIZ                           !input value, number of data horizons
    INTEGER :: Num_delx(max_horizons)           !Input Number of discretizations in horizon
    REAL    :: bd_input(max_horizons)           !bulk density of each horizon input
    REAL    :: oc_input(max_horizons)           !organic carbon input Percent %

    REAL    :: fc_input(max_horizons)           !field capacity input
    REAL    :: wp_input(max_horizons)           !wilting point input
    REAL    :: soil_temp_input(max_horizons)    !input of initial soil temperatures
    REAL    :: thickness(max_horizons)          !thickness of each horizon, data input not reformulated by auto
    REAL    :: dispersion_input(max_horizons)   !dispersion


    !Optional PARAMETERs for discretization
    INTEGER, PARAMETER :: max_discretized_layers = 100 !new user defined  discretized layering

    ! use the specified discetization and let program calcuate properties, ignore delta in soil profile
    LOGICAL :: is_auto_profile
    INTEGER :: number_of_discrete_layers  ! number of discrete layers with different discretizations
    REAL    :: profile_thick(max_discretized_layers)              ! thickness of uniquely discretized layer
    INTEGER :: profile_number_increments(max_discretized_layers)  ! number of increments in a layer

    REAL ::  hydrolysis_halflife_input(3)       ! halflife in per days
    REAL ::  hydrolysis_rate(3)                 ! rate in per sec

    REAL :: water_column_halflife_input(3)      ! Degradation halflife in per day
    REAL :: water_column_rate(3)                ! converted to per sec
    REAL :: water_column_ref_temp(3)

    REAL :: benthic_halflife_input(3)     ! Degradation halflife in per day
    REAL :: benthic_rate(3)               ! converted to per sec
    REAL :: benthic_ref_temp(3)

    REAL :: photo_halflife_input(3)    ! Degradation halflife in days
    REAL :: photo_rate(3)              ! Rate in per sec
    REAL :: rflat(3)                   ! input latitude for photolysis study

    LOGICAL :: is_hydrolysis_override  ! if true THEN hydrolysis will overide the aquatic phase portion of soil meatbolism if
    ! hydrolyis is greater than soil metabolism, solid phase metabolism remains the same

    LOGICAL :: is_total_degradation    ! True idf total soil degradation, False if aqueous phase only degradatiuon
    REAL    :: soil_degradation_halflife_input(3), soil_ref_temp(3), xsoil(3)

    REAL    :: Aq_rate_input(3)    ! not REALly an input, soil aqueous degradation rate parent, daughter, grand daughter, per day
    REAL    :: Sorb_rate_input(3)  ! not REALly input soil sorbed degradation rate parent, daughter, grand daughter
    REAL    :: Gas_rate_input(3)   ! not REALly input soil as degradation rate parent, daughter, grand daughter


    REAL    :: k_f_input(3)        ! input Freundlich Kf for each horizon
    REAL    :: N_f_input(3)        ! input Freundlich N for each horizon

    !nonequilibrium PARAMETERs
    REAL :: k_f_2_input(3)  ! input Freundlich Kf for each horizon Region 2
    REAL :: N_f_2_input(3)  ! input Freundlich N for each horizon Region 2


    !Soil Degradation Conversion Factors---Changed to Mass Conversion Factors on 1/13/2024 by multiplying MWT ratios
    REAL :: MolarConvert_aq12_input(max_horizons),MolarConvert_aq13_input(max_horizons),MolarConvert_aq23_input(max_horizons)
    REAL :: MolarConvert_s12_input(max_horizons), MolarConvert_s13_input(max_horizons), MolarConvert_s23_input(max_horizons)
    ! Below this concentration (mg/L), Freundlich isotherm is approximated to be linear to prevent numerical issues
    REAL :: lowest_conc

    !Compartment-Specific PARAMETERs

    REAL,ALLOCATABLE,DIMENSION(:,:) :: k_freundlich   ! Freundich Coefficient in equilibrium Region
    REAL,ALLOCATABLE,DIMENSION(:,:) :: N_freundlich   ! Freundich Exponent in equilibrium Region
    REAL,ALLOCATABLE,DIMENSION(:,:) :: k_freundlich_2 ! Freundich Coefficient in Nonequilibrium Region 2
    REAL,ALLOCATABLE,DIMENSION(:,:) :: N_freundlich_2 ! Freundich Exponent in Nonequilibrium Region 2

    REAL,ALLOCATABLE,DIMENSION(:) :: soil_depth       ! vector of soil depths
    REAL,ALLOCATABLE,DIMENSION(:) :: wiltpoint_water  ! wilting point content weighted by size of compartment (absolute water content)
    REAL,ALLOCATABLE,DIMENSION(:) :: fieldcap_water   ! field capacity content weighted by size of compartment (absolute water content)
    REAL,ALLOCATABLE,DIMENSION(:) :: delx
    REAL,ALLOCATABLE,DIMENSION(:) :: bulkdensity

    REAL,ALLOCATABLE,DIMENSION(:) :: orgcarb
    REAL,ALLOCATABLE,DIMENSION(:) :: theta_zero  ! beginning day water content
    REAL,ALLOCATABLE,DIMENSION(:) :: theta_end   ! end day water content
    REAL,ALLOCATABLE,DIMENSION(:) :: theta_sat   ! stauration water content (fractional)
    REAL,ALLOCATABLE,DIMENSION(:) :: theta_fc    ! field capacity water content (fractional)
    REAL,ALLOCATABLE,DIMENSION(:) :: theta_wp    ! wilting point water content (fractional)

    REAL,ALLOCATABLE,DIMENSION(:) :: soil_temp          ! soil temperature local day only
    REAL,ALLOCATABLE,DIMENSION(:,:) :: soil_temp_save   ! soil temperature all days all compartments
    REAL,ALLOCATABLE,DIMENSION(:) :: soilwater


    ! Note: the soil degradation rates below have a correction for PRZMs implicit routine,
    ! Thus they are not applicable for other routines without implict calcs
    ! this is problematic for TPEZ so TPEZ needs separate soil degradation calcs
    ! current soil degradation rates for water solid gas, CORRECTED FOR IMPLICIT ROUTINE, per day
    REAL,ALLOCATABLE,DIMENSION(:,:) :: dwrate, dsrate, dgrate
    REAL,ALLOCATABLE,DIMENSION(:,:) :: dwrate_atRefTemp, dsrate_atRefTemp, dgrate_atRefTemp ! input degradation rates at ref temp

    ! For TPEZ calculation of soil degrdataion, these will be used to "uncorrect" the implicit correction
    ! So they are now in this global module rather than local
    REAL :: aq_rate_corrected(3)  ! degradation inputs corrected for implicit routine
    REAL :: sorb_rate_corrected(3)
    REAL :: gas_rate_corrected(3)


    REAL,ALLOCATABLE,DIMENSION(:) :: MolarConvert_aq12,MolarConvert_aq13,MolarConvert_aq23
    REAL,ALLOCATABLE,DIMENSION(:) :: MolarConvert_s12, MolarConvert_s13, MolarConvert_s23

    REAL,ALLOCATABLE,DIMENSION(:) :: dispersion  ! formerly disp
    REAL,ALLOCATABLE,DIMENSION(:) :: thair_old
    REAL,ALLOCATABLE,DIMENSION(:) :: thair_new

    REAL,ALLOCATABLE,DIMENSION(:,:) :: Kd_new
    REAL,ALLOCATABLE,DIMENSION(:,:) :: Kd_old    ! current time step Kd (applies to nnlinear isothrems)
    REAL,ALLOCATABLE,DIMENSION(:,:) :: Kd_2

    REAL,ALLOCATABLE,DIMENSION(:,:) :: new_henry  ! henry constant for current time step (end of current step)
    REAL,ALLOCATABLE,DIMENSION(:,:) :: old_Henry  ! Henry constant previous time step (start of current step)

    REAL,ALLOCATABLE,DIMENSION(:,:) :: Sorbed2  ! nonequilib phase conc (internal units of g/g to correspond with mg/L aqueous)

    REAL,ALLOCATABLE,DIMENSION(:,:) :: conc_porewater
    REAL,ALLOCATABLE,DIMENSION(:,:) :: conc_total_per_water     ! pest conc in Total Mass in all phases divided by water content)
    REAL,ALLOCATABLE,DIMENSION(:,:) :: mass_in_compartment      ! New experimetal variable to replace "conc_total_per_water"
    REAL,ALLOCATABLE,DIMENSION(:,:) :: mass_in_compartment2     ! mass in nonequilibrium compartment
    ! REAL,ALLOCATABLE,DIMENSION(:,:) :: SOILAP  !g/cm2  Mass of applied pesticide in soil

    REAL,ALLOCATABLE,DIMENSION(:) :: DGAIR !effective air diffusion coefficient through profile for chemical in the loop


    ! ********************************************************************************************************************
    !****Pesticide Flux**************

    REAL,ALLOCATABLE,DIMENSION(:,:) :: SRCFLX  !Production of chemical. i.e., from degrdation
    !*******STORED OUTPUT VARIABLES**************************************

    REAL :: ERFLUX(3), WOFLUX(3), ROFLUX(3)
    REAL,ALLOCATABLE,DIMENSION(:,:) :: DKFLUX
    REAL,ALLOCATABLE,DIMENSION(:,:) :: PVFLUX
    REAL,ALLOCATABLE,DIMENSION(:,:) :: UPFLUX
    REAL :: SDKFLX(3)     = 0.0
    REAL :: DCOFLX(3)     = 0.0   !Pesticide outflow below soil core
    REAL,ALLOCATABLE,DIMENSION(:,:) :: soil_applied_washoff


    REAL, PARAMETER :: washoff_incorp_depth  = 2.0 !depth at which foliar washoff pesticide is incorporated into soil
    !previously this value was fixed to the runoff extraction depth
    INTEGER :: number_washoff_nodes  !number of nodes coreresponding to washoff_incorp_depth

    ! INTEGER :: NCOM1  !holds the lowest relavant evaopration node for each day
    ! INTEGER :: NCOM0  !holds the node of corresponding to anetd

    INTEGER :: min_evap_node  !holds the node of corresponding to anetd (min_evap_depth)

    REAL, PARAMETER :: cam123_soil_depth = 4.0   ! cm  Default pesticide incorporation depth
    REAL, PARAMETER :: CN_moisture_depth = 10.0  ! depth of soil moisture for CN calculations

    REAL :: r1
    REAL :: R0MIN  = Tiny(r1)

    LOGICAL :: REALly_not_thrufl  !flag to prevent washoff of undercanopy irrigation in PRZM3 (not used in PRZM5)
    REAL    :: canopy_flow        !flow through the canopy
    REAL    :: effective_rain     !total water, without consideration for canopy holdup
    REAL    :: THRUFL             !total water minus canopy holdup

    !***Simulation Timne **************************************************
    INTEGER  :: startday      !julian day start date referenced to 1900 (note that we use two digit date, makes no difference)
    INTEGER  :: julday1900    !the current day referenced to Jan 1, 1900
    INTEGER  :: first_year
    INTEGER  :: last_year
    INTEGER  :: num_years    !number of caLENdar years in simulation (number of different year numbers, can be incomplete years)
    INTEGER  :: first_mon  !first month of simulation
    INTEGER  :: first_day    !first day of simulation


    !***PRZM5 Advanced Analysis OPtipon
    LOGICAL :: is_Freundlich !True if Freundlich
    LOGICAL :: is_nonequilibrium
    LOGICAL :: calibrationFlag =.FALSE.

    LOGICAL :: ADJUST_CN

    !*****Pesticide Applications ******************

    ! These PARAMETERs are arrays that have each application for the entire simulation, day by day
    ! juilan date (1900 references) array of application dates
    INTEGER,ALLOCATABLE,DIMENSION(:) :: application_date 
    !holds the original application date so that "application_date" can be modified during runs
    INTEGER,ALLOCATABLE,DIMENSION(:) :: application_date_original

    ! 1 is default 4cm; 2 default above 4cm; 3 uniform; 4 Tband 5 custm 6 custom
    INTEGER,ALLOCATABLE,DIMENSION(:)  :: pest_app_method
    REAL,ALLOCATABLE,DIMENSION(:)     :: DEPI
    ! READ in application_rate_in as KG/HA but converted to G/CM**2 in initialization
    REAL,ALLOCATABLE,DIMENSION(:)     :: TAPP
    REAL,ALLOCATABLE,DIMENSION(:)     :: APPEFF
    REAL,ALLOCATABLE,DIMENSION(:)     :: Tband_top
    REAL,ALLOCATABLE,DIMENSION(:)     :: drift_kg_per_m2          ! the drift application rate to pond
    ! REAL,ALLOCATABLE,DIMENSION(:)   :: tpez_drift_kg_per_m2     ! the drift application rate to tpez

    ! array of indices to order the applications chronoLOGICALly, needed only if the rain fast option is to be used.
    INTEGER,ALLOCATABLE,DIMENSION(:) :: application_order 

    INTEGER :: total_applications      ! Total applications for simulations


    ! These PARAMETERs are the application PARAMETERs that are READ in from the input file.
    ! They represent a short-hand of the actual applications above

    ! number of apps reported in the input file (shorthand, does not include all individual apps ocurring every year)
    INTEGER :: num_applications_input

    INTEGER,ALLOCATABLE,DIMENSION(:) :: pest_app_method_in
    REAL,ALLOCATABLE,DIMENSION(:)    :: DEPI_in
    REAL,ALLOCATABLE,DIMENSION(:)    :: application_rate_in

    REAL,ALLOCATABLE,DIMENSION(:)    :: APPEFF_in  ! no longer an input, calculated based on spray losses
    REAL,ALLOCATABLE,DIMENSION(:)    :: Tband_top_in

    ! ********************************************************************************************************************
    REAL,ALLOCATABLE,DIMENSION(:)    :: drift_value  !no longer an input----to be modified once spray table completed
    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    !spray table coordinates:
    INTEGER,ALLOCATABLE,DIMENSION(:) :: drift_row
    REAL,ALLOCATABLE,DIMENSION(:)    :: drift_buffer !column

    INTEGER,ALLOCATABLE,DIMENSION(:) :: lag_app_in
    INTEGER,ALLOCATABLE,DIMENSION(:) :: repeat_app_in

    INTEGER,ALLOCATABLE,DIMENSION(:) :: days_until_applied          !for relative dates


    ! g/cm2, keeps track of the cumulative mass applied, used in output CHARACTERizations, zero it for every run
    REAL :: applied_mass_sum_gram_per_cm2
    REAL :: fraction_off_field             !fraction of applied mass that leaves field by runoff, erosion, and drift

    !*** Application Optimizing for Rainfall *********
    LOGICAL :: is_adjust_for_rain       !True means to adjust the applications dates
    REAL    :: rain_limit               !cm
    INTEGER :: optimum_application_window
    INTEGER :: intolerable_rain_window
    INTEGER :: min_days_between_apps



    REAL    :: plant_app     !store the plant applied pesticide
    LOGICAL :: some_applications_were_foliar

    !*****Pesticide Properties ******************
    REAL    :: foliar_halflife_input(3)           !plant decay halflife input

    REAL    :: plant_pesticide_degrade_rate(3)    !plant decay rate  --convert to this in initialization
    REAL    :: plant_washoff_coeff(3)             !plant washoff coefficient (per cm rainfall)




    REAL :: mwt(3)
    REAL :: solubilty(3)
    REAL :: vapor_press(3)
    REAL :: Henry_unitless(3)     !Henry's Constant ConcAir/ConcWater volumetric   now called HenryK same as PRZM
    REAL :: Heat_of_Henry(3)      !Enthalpy of air water exchange (J/mol)
    REAL :: ENPY(3)             ! convert enthalpy from Joules to Kcal for PRZM  4184 J/kcal



    INTEGER :: PCMC                               !signal to tell if its Koc to use
    INTEGER :: KDFLAG                             !Select Kd calculation
    INTEGER :: THFLAG                             !computw Wp and FC based on texture
    REAL    :: k2(3)                              !noneq mass transfer coeff


    REAL :: foliar_formation_ratio_12, foliar_formation_ratio_23 !foliar transformation rates
    REAL :: plant_volatilization_rate(3)   !plant pesticide volatilization rate (does not produce degradate)


    REAL :: UPTKF(3)         !PLANT UPTAKE FACTOR, CONVERTED TO GAMMA1 IN PLANTGROW

    !********Volatilization
    REAL           :: Height_stagnant_air_layer_cm = 0.5  !(can be changed by input file)
    REAL           :: uWind_Reference_Height = 10.0       !(can be changed by input file)
    REAL,PARAMETER :: vonKarman = 0.4

    REAL :: CONDUC(3)        !Some canopy volatilization PARAMETER
    REAL :: DAIR(3)
    !    REAL :: HENRYK(3)

    REAL,PARAMETER :: henry_ref_temp  = 25.0

    REAL    :: Q_10   !temperature basis for formerly QFAC
    !    REAL    :: TBASE  !temperature reference for Q_10 degradation

    INTEGER :: NCHEM


    !*******Heat transfer variables*********************
    REAL,PARAMETER    :: emmiss =0.97
    REAL    :: UBT                            !used in volatilization and heat routines
    REAL    :: Albedo(13), bottom_bc(13) ! albedo and bottom boundary condition constabt temperature
    LOGICAL :: is_temperature_simulated
    REAL    :: enriched_eroded_solids  !mass of eroded solids bumped up be an enrichment factor,  g/cm2, (analogous to runoff)
    REAL    :: TCNC(3)     !used for output only

    !********* Groundwater Related Variables
    REAL :: retardation_factor(3)   !retardation factor of entire soil profile (parent, daughter, grand); exclude last 2 nodes
    REAL :: maxcap_volume           !summation of max capacities, this is the effective aqueous transport pore volume
    INTEGER :: top_node_last_horizon, bottom_node_last_horizon
    REAL :: gw_peak(3)
    REAL :: post_bt_avg(3)
    REAL :: throughputs(3)
    REAL :: simulation_avg(3)

    !******Weather Related Variables******************
    REAL :: precip_rain     !rain-only component, minus snow
    REAL,PARAMETER :: SFAC = 0.274  !USDA value
    REAL :: SNOWFL
    REAL,PARAMETER :: PFAC = 1.0   !originally used to adjust pan evap to field evaporation
    REAL :: open_water_adj         !allows tthe evaporation in waterbody to be adjusted (like pfac for field)
    INTEGER :: num_records         !number of records in weather file

    REAL, ALLOCATABLE, DIMENSION(:) :: precip, pet_evap,air_temperature, wind_speed, solar_radiation  !from weather file
    REAL, ALLOCATABLE, DIMENSION(:) :: precip_m,evap_m , temp_avg, wind_m  ! with VVWM units, all meters, sec, temp is 30-day avg



    !Weather Scalers

    REAL :: wind
    REAL :: precipitation,air_TEMP,PEVP,SOLRAD !old ones

    !*****Irrigation PARAMETERs
    INTEGER :: IRFLAG
    INTEGER :: IRTYPE

    REAL    :: PCDEPL,FLEACH      !inputs
    REAL    :: max_irrig  !maximum  irrigation water over 24 hour period
    REAL    :: under_canopy_irrig, over_canopy_irrig

    REAL    :: IRRR                      !daily irrigation water, used for time series reporting only


    REAL, ALLOCATABLE, DIMENSION(:)   :: enriched_erosion_save

    REAL, ALLOCATABLE, DIMENSION(:)   :: irrigation_save
    REAL, ALLOCATABLE, DIMENSION(:)   :: canopy_flow_save
    REAL, ALLOCATABLE, DIMENSION(:,:) :: thair_save   !(ncom2, num_records)
    REAL, ALLOCATABLE, DIMENSION(:,:) :: theta_end_save
    REAL, ALLOCATABLE, DIMENSION(:,:) :: soilwater_save
    REAL, ALLOCATABLE, DIMENSION(:,:) :: velocity_save
    REAL, ALLOCATABLE, DIMENSION(:,:) :: theta_zero_save

    REAL, ALLOCATABLE, DIMENSION(:,:) :: thair_old_save

    REAL, ALLOCATABLE,DIMENSION(:,:)   ::  infiltration_save


    LOGICAL :: UserSpecifiesDepth        !True if custom depth, False if auto with root depth
    REAL    :: user_irrig_depth          !user-specified depth of irrigation, input value
    INTEGER :: user_irrig_depth_node     !node for user-specified depth of irrigation, calc in INITL subroutine

    !*****Crop PARAMETERs *********

    REAL    :: min_evap_depth  ! the minimum depth that is used for PET satisfaction
    !VECTORS ODF SIZE NUM_RECORDS
    REAL,   ALLOCATABLE, DIMENSION(:) :: crop_fraction_of_max  !vector of daily fration of crop growth
    REAL,   ALLOCATABLE, DIMENSION(:) :: canopy_cover          !Fractional Coverage (but in PWC interface it is %)
    REAL,   ALLOCATABLE, DIMENSION(:) :: canopy_height
    REAL,   ALLOCATABLE, DIMENSION(:) :: canopy_holdup
    REAL,   ALLOCATABLE, DIMENSION(:) :: root_depth
    INTEGER,ALLOCATABLE, DIMENSION(:) :: evapo_root_node  !Node that represents the root depth or evap zone (whichever greatest) Ncom1
    INTEGER,ALLOCATABLE, DIMENSION(:) :: root_node        !Node for the root (formerly nrzomp)
    INTEGER,ALLOCATABLE, DIMENSION(:) :: atharvest_pest_app

    LOGICAL,ALLOCATABLE, DIMENSION(:) :: is_harvest_day !True if the current day is a harvest day

    !VECTORS OF SIZE OF INPUT READ

    INTEGER,DIMENSION(max_number_crop_periods) :: emd   !these are the input values for day month
    INTEGER,DIMENSION(max_number_crop_periods) :: emm
    INTEGER,DIMENSION(max_number_crop_periods) :: mad
    INTEGER,DIMENSION(max_number_crop_periods) :: mam
    INTEGER,DIMENSION(max_number_crop_periods) :: had
    INTEGER,DIMENSION(max_number_crop_periods) :: ham
    INTEGER,DIMENSION(max_number_crop_periods) :: foliar_disposition       !1 = pesticide is surface applied at harvest, 2 = complete removal, 3 = left alone, 4 = bypass, 0 = ????
    REAL   ,DIMENSION(max_number_crop_periods) :: max_root_depth
    REAL   ,DIMENSION(max_number_crop_periods) :: max_canopy_cover
    REAL   ,DIMENSION(max_number_crop_periods) :: max_canopy_holdup
    REAL   ,DIMENSION(max_number_crop_periods) :: max_canopy_height
    INTEGER,DIMENSION(max_number_crop_periods) :: crop_lag
    INTEGER,DIMENSION(max_number_crop_periods) :: crop_periodicity

    INTEGER,ALLOCATABLE,DIMENSION(:,:) :: emergence_date  !(Crop#), these are the processed julian1900 dates
    INTEGER,ALLOCATABLE,DIMENSION(:,:) :: maturity_date
    INTEGER,ALLOCATABLE,DIMENSION(:,:) :: harvest_date
    LOGICAL :: evergreen


    !Scalers for core przm run
    INTEGER :: root_node_daily
    INTEGER :: evapo_root_node_daily       !the nodes to the bottom of evap/root zone
    REAL    :: cover                       !Fractional Coverage (but in PWC interface it is %)
    REAL    :: height
    INTEGER :: harvest_placement
    LOGICAL :: harvest_day

    INTEGER :: num_crop_periods_input !number of cropping periods READ from the input file,not the actual number in the simulation
    !   INTEGER :: total_crop_periods !total number of crop periods in the simulation, considering lag and periodicity and years

    !*****Erosion PARAMETERs *************************
    !   REAL :: AFIELD           !INPUT, field area input, m2---NOW IN WATERBODY PARAMETERs


    REAL :: AFIELD_ha        !   Area in ha, used in erosion routine

    REAL :: USLEK            !INPUT,
    REAL :: USLEP            !INPUT,
    REAL :: USLELS           !INPUT,
    !  REAL :: hydro_LENgth     !INPUT,
    REAL :: SLP

    INTEGER :: erflag         !erosion flag

    INTEGER :: CN_index                    !Current index for ereosion and runoff PARAMETERs, current one being used
    !  REAL    :: N1                          !The current mannings n value
    REAL    :: CFAC                        !The current c factor
    REAL    :: USLEC(Num_hydro_factors)    !inpuit values for c facator
    REAL    :: MNGN(Num_hydro_factors)

    LOGICAL :: use_usleyears               !True means usle values are year specific and do not repeat
    INTEGER :: NUSLEC                      !number of c factors
    INTEGER :: GDUSLEC(Num_hydro_factors)  !day
    INTEGER :: GMUSLEC(Num_hydro_factors)  !month
    INTEGER :: GYUSLEC(Num_hydro_factors)  !year when this option is used



    INTEGER :: JUSLEC(Num_hydro_factors)   !julian day for each of the erosion days converted from gduslec and gmuslec


    REAL :: CN_2(Num_hydro_factors) !Curve Number 2 READ from inputs, Now a REAL number ( was an INTEGER because its used as index in lookup table )



    INTEGER :: IREG        !erosion/rain intensity map region


    !******Hydrology*********
    REAL :: STTDET              !evaporation amount in the top 5 cm, used for temperature calculations
    REAL :: CEVAP               !HOLDS ACTUAL PET
    REAL :: SMELT               !snow melt
    REAL :: CN_moisture_ref     !water content halfway between wp and fc averaged over cn depth
    REAL :: INABS
    REAL :: Infiltration_Out   !Infiltration out of last compartment, this info was previously not captured by PRZM
    INTEGER :: cn_moist_node   !number of compartments that make up runff moisture depth


    !REAL :: THETH    dfy renamed CN_moisture_ref  8/24/17          !water content halfway between wp and fc averaged over cn depth (not sure of purpose)

    !***Accumulators***************

    REAL, ALLOCATABLE,DIMENSION(:)   :: ainf        !infiltration INTO compartment(needs to be DIMENSIONed as ncom2+1 to capture outflow)
    REAL, ALLOCATABLE,DIMENSION(:)   :: vel         !velocity out of compartment = ainf(i+1)/theta
    REAL, ALLOCATABLE,DIMENSION(:)   :: EvapoTran

    REAL, ALLOCATABLE,DIMENSION(:,:) :: GAMMA1
    REAL :: snow              !accumulated snow for each day, calls itself do not initialize here
    REAL :: CINT              !the actual water held on the canopy during the current day, calls itself do not initialize here

    REAL :: runoff_on_day
    REAL :: foliar_degrade_loss(3)                  !foliar loss by degradation
    REAL :: FOLPST(3)                               !foliar storage
    REAL :: SUPFLX(3)                               !plant uptake

    REAL :: Foliar_volatile_loss(3)                 !for output only volatilization, ONLY defined for parent,i never did calcs for degradates

    REAL :: potential_canopy_holdup                 !daily maximum water that the canopy can hold

    REAL :: TDET                                    !evapotranspiration in the top ncom1 compartments

    REAL :: SEDL




    !      INTEGER :: MTR1  !some kind of comaprtment count = ncom2


    !****Miscelaneous******


    !CHARACTER(LEN=78) ::  TITLE     !INPUTused in EXAMS output files and time series
    CHARACTER(LEN=20) ::  PSTNAM(3) !INPUT used in EXAMS

    !INTEGER :: current_day        !day of month
    !INTEGER :: current_month
    !INTEGER :: current_year       !taken from metfile

    INTEGER  :: CFLAG                 !INPUT initial pesticide concentration conversion flag

    CHARACTER(LEN=4)  :: MODE(max_number_plots)      !INPUT time series plot modes
    CHARACTER(LEN=4)  :: PLNAME(max_number_plots)    !INPUT Plot IDs for time series
    INTEGER           :: ARG(max_number_plots)      !Formerly IARG, INPUT arguments for time series plots
    INTEGER           :: ARG2(max_number_plots)      !INPUT arguments for time series plots

    INTEGER           :: NPLOTS                       !INPUT number of time series plots
    REAL              :: CONST(max_number_plots)             !INPUT, the user defined adjustment to the time series output
    REAL              :: OUTPUJ(max_number_plots)
    REAL              :: OUTPJJ(max_number_plots)     !output accumulator
    INTEGER           :: chem_id(max_number_plots)    !plot CHARACTERizer chem 1 2 3


    CHARACTER(LEN=4)  :: temp_MODE(max_number_plots)      !INPUT time series plot modes
    CHARACTER(LEN=4)  :: temp_PLNAME(max_number_plots)    !INPUT Plot IDs for time series
    INTEGER           :: temp_chem_id(max_number_plots)   !plot CHARACTERizer chem 1 2 3
    INTEGER           :: temp_ARG(max_number_plots)       !Formerly IARG, INPUT arguments for time series plots
    INTEGER           :: temp_ARG2(max_number_plots)      !INPUT arguments for time series plots
    REAL              :: temp_CONST(max_number_plots)     !INPUT, the user defined adjustment to the time series output
    INTEGER           :: extra_plots                      !number of user-specified extra plots
    LOGICAL           :: is_timeseriesfile                !create a time series file




    REAL   ::  curve_number_daily   !holds curve number for plotting

    !***** PRZM5 runoff and erosion extraction variables ****************
    REAL :: runoff_effic              !amount of runoff bypassing surface soil
    REAL :: runoff_decline            !exponential factor for interaction decline with depth
    REAL :: runoff_extr_depth         !depth of runoff interaction
    REAL, ALLOCATABLE,DIMENSION(:) :: runoff_intensity  !PRZM5 uses this instead of DRI, fraction of runoff per cm depth

    REAL :: erosion_effic             !propotional reduction in erosion intensity
    REAL :: erosion_decline           !exponential factor for interaction decline with depth
    REAL :: erosion_depth             !erosion depth
    REAL, ALLOCATABLE,DIMENSION(:) :: erosion_intensity

    INTEGER :: RNCMPT             !the number of compartments corresponding to the runoff extraction depth
    INTEGER :: erosion_compt      !number of compartments for rerosion extraction


    CHARACTER(LEN= 256) :: outputfile_parent_daily
    CHARACTER(LEN= 256) :: outputfile_deg1_daily
    CHARACTER(LEN= 256) :: outputfile_deg2_daily

    !  CHARACTER(LEN= 256) :: outputfile_parent_analysis
    !  CHARACTER(LEN= 256) :: outputfile_deg1_analysis
    CHARACTER(LEN= 256) :: outputfile_deg2_analysis

    CHARACTER(LEN= 256) :: outputfile_parent_deem
    CHARACTER(LEN= 256) :: outputfile_deg1_deem
    CHARACTER(LEN= 256) :: outputfile_deg2_deem

    CHARACTER(LEN= 256) :: outputfile_parent_caLENdex
    CHARACTER(LEN= 256) :: outputfile_deg1_caLENdex
    CHARACTER(LEN= 256) :: outputfile_deg2_caLENdex

    CHARACTER(LEN= 256) :: outputfile_parent_esa
    CHARACTER(LEN= 256) :: outputfile_deg1_esa
    CHARACTER(LEN= 256) :: outputfile_deg2_esa


    !VVWM PARAMETERs
    REAL,ALLOCATABLE,DIMENSION(:) :: fraction_to_benthic
    REAL,ALLOCATABLE,DIMENSION(:) :: eroded_solids_mass   !READ in from PRZMV output as tonnes/day but immediately converted to kg/day
    REAL,ALLOCATABLE,DIMENSION(:) :: burial               !(kg/sec) array of daily burial rate
    REAL,ALLOCATABLE,DIMENSION(:) :: flowthru_the_body    !(m3/sec) array of daily flow
    REAL,ALLOCATABLE,DIMENSION(:) :: edge_of_field        !(kg/m3) array of daily runoff concentrations (zero if no runoff) includes base flow dilution if any



    REAL,ALLOCATABLE,DIMENSION(:) :: k_aer_aq     !aqueous-phase aerobic rate (per sec)
    REAL,ALLOCATABLE,DIMENSION(:) :: k_anaer_aq   !aqueous-phase anaerobic rate (per sec)
    REAL,ALLOCATABLE,DIMENSION(:) :: k_aer_s      !sorbed-phase aerobic rate (per sec)
    REAL,ALLOCATABLE,DIMENSION(:) :: k_anaer_s    !sorbed-phase anaerobic rate (per sec)
    REAL,ALLOCATABLE,DIMENSION(:) :: k_volatile   !first order volatilization rate (per sec)
    REAL,ALLOCATABLE,DIMENSION(:) :: k_photo      !photolysis rate (1/sec)
    REAL,ALLOCATABLE,DIMENSION(:) :: daily_depth  !daily water body depths
    REAL,ALLOCATABLE,DIMENSION(:) :: k_hydro      !hydrolysis rate (per sec)
    REAL,ALLOCATABLE,DIMENSION(:) :: k_flow       !array of daily wash out rates (per second)
    REAL,ALLOCATABLE,DIMENSION(:) :: k_burial
    REAL,ALLOCATABLE,DIMENSION(:) :: theta        !solute holding capacity ratio [--]
    REAL,ALLOCATABLE,DIMENSION(:) :: capacity_1   !solute holding capacity of region 1 [m3]
    REAL,ALLOCATABLE,DIMENSION(:) :: degradateProduced1,degradateProduced2
    REAL,ALLOCATABLE,DIMENSION(:) :: fw1          !fraction of region 1 solute in aqueous phase
    REAL,ALLOCATABLE,DIMENSION(:) :: gamma_1              !effective littoral degradation
    REAL,ALLOCATABLE,DIMENSION(:) :: gamma_2              !effective benthic degradation

    REAL                          :: capacity_2   !solute holding capacity of region 2 [m3]
    REAL                          :: fw2          !fraction of region 2 solute in aqueous phase
    REAL                          :: kd_sed_1     !local Kd of ss (m3/kg)
    REAL                          :: omega        !mass transfer coefficient

    REAL                              :: m_sed_1                           !base amount of suspended solids mass (kg)
    REAL,ALLOCATABLE,DIMENSION(:)     :: volume1                           !array of daily water column volumes  [m3]
    REAL:: v2                                                              !aqueous volume of pore water in sediment
    REAL, ALLOCATABLE, DIMENSION(:)    :: A,B,E,F                          !final coefficients for the 2 simultaneous equations
    REAL, ALLOCATABLE, DIMENSION(:)    :: m1_input, m2_input               !at start of time step: the mass input in litt and benthic
    REAL, ALLOCATABLE, DIMENSION(:)    :: m1_store,m2_store,mavg1_store    !array of daily peak/avg
    REAL, ALLOCATABLE, DIMENSION(:)    :: m_total                          !total average daily mass in system: mavg1 + mavg2
    REAL, ALLOCATABLE, DIMENSION(:)    :: aq1_store, aq2_store             !begining day concentrations, after mass inputs

    REAL,ALLOCATABLE,DIMENSION(:,:,:) :: mass_off_field            !daily mass loading from runoff (column 1) & erosion (column 2) (kg)
    REAL,ALLOCATABLE,DIMENSION(:)     :: spray_additions           !daily mass of spray to be added to water column (kg)
    REAL,ALLOCATABLE,DIMENSION(:)     :: aqconc_avg1 ,aqconc_avg2  !daily average aqueous concentrations


    REAL:: spray_total(3)
    REAL:: runoff_total(3)
    REAL:: erosion_total(3)

    REAL :: Sediment_conversion_factor(3) !Converts pore water to concentration (total mass)/(dry solid mass)
    REAL :: Daily_Avg_Runoff
    REAL :: Daily_avg_flow_out


    REAL :: runoff_fraction  !fraction of chemical transport due to runoff
    REAL :: erosion_fraction
    REAL :: drift_fraction
    LOGICAL :: First_time_through_wb(3)        !used for batch READer
    LOGICAL :: First_time_through_wpez         !used for batch READer
    LOGICAL :: First_time_through_tpez         !used for batch READer
    LOGICAL :: First_time_through_PRZM         !used to WRITE output przm time series file headers, so can keep all output WRITEs in one place
    LOGICAL :: First_time_through_medians      ! median output file
    LOGICAL :: First_time_through_medians_wpez
    LOGICAL :: First_time_through_medians_tpez



    !OUTPUT FLAGS
    LOGICAL :: is_runoff_output
    LOGICAL :: is_erosion_output
    LOGICAL :: is_runoff_chem_output
    LOGICAL :: is_erosion_chem_output

    LOGICAL :: is_conc_bottom_output
    LOGICAL :: is_daily_volatilized_output
    LOGICAL :: is_cumulative_volatilized_output

    LOGICAL :: is_daily_chem_leached_output
    REAL    :: leachdepth

    LOGICAL :: is_chem_decayed_part_of_soil_output
    REAL :: decay_start, decay_end

    LOGICAL :: is_chem_in_all_soil_output
    LOGICAL :: is_chem_in_part_soil_output
    REAL :: fieldmass_start, fieldmass_end

    LOGICAL :: is_chem_on_foliage_output

    LOGICAL :: is_precipitation_output
    LOGICAL :: is_evapotranspiration_output
    LOGICAL :: is_soil_water_output
    LOGICAL :: is_irrigation_output

    LOGICAL :: is_infiltration_at_depth_output
    REAL    :: infiltration_point
    LOGICAL :: is_infiltrated_bottom_output

    LOGICAL :: is_waterbody_info_output
    !LOGICAL :: is_waterbody_depth_output
    !LOGICAL :: is_waterbody_concen_output
    !LOGICAL :: is_waterbody_porewater_output


    LOGICAL :: is_output_spraydrift
    LOGICAL :: is_gw_btc



    LOGICAL::  is_constant_profile, is_ramp_profile, is_exp_profile
    REAL :: ramp1, ramp2, ramp3, exp_profile1, exp_profile2

    LOGICAL, ALLOCATABLE, DIMENSION(:) :: is_app_window
    INTEGER, ALLOCATABLE, DIMENSION(:) :: app_window_span, app_window_step


END MODULE constants_and_variables