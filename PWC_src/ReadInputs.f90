MODULE READinputs
    IMPLICIT NONE

CONTAINS

    ! ******************************************************************
    SUBROUTINE READ_inputfile
    ! NEW READ for new input file, READing same file as the PWC Input file ******
    USE constants_and_variables, ONLY: inputfile, inputfile_unit_number, is_koc, is_freundlich, is_nonequilibrium,             &
        k_f_input, N_f_input, k_f_2_input, N_f_2_input, lowest_conc, number_subdelt, k2, water_column_ref_temp, xAerobic,      &
        benthic_halflife_input, benthic_ref_temp, xBenthic, photo_halflife_input, xPhoto, hydrolysis_halflife_input, rflat,    &
        xhydro, soil_degradation_halflife_input, soil_ref_temp, xsoil, foliar_halflife_input, foliar_formation_ratio_12,       &
        foliar_formation_ratio_23, plant_washoff_coeff, mwt, vapor_press, solubilty, dair,Henry_unitless, Heat_of_Henry, q_10, &
        number_of_schemes, num_apps_in_schemes, number_of_scenarios, app_reference_point_schemes, application_rate_schemes,    &
        depth_schemes, split_schemes, drift_schemes, lag_schemes, periodicity_schemes, driftfactor_schemes, method_schemes,    &
        days_until_applied_schemes, scenario_names, working_directory, family_name, weatherfiledirectory, erflag, nchem,       &
        ADJUST_CN, water_column_halflife_input, is_runoff_output, is_erosion_output, is_runoff_chem_output,                    &
        is_erosion_chem_output, is_conc_bottom_output, is_daily_volatilized_output, is_daily_chem_leached_output, leachdepth,  &
        is_chem_decayed_part_of_soil_output, decay_start, decay_end, is_chem_in_part_soil_output, fieldmass_start,             &
        fieldmass_end, is_chem_on_foliage_output, is_precipitation_output, is_evapotranspiration_output, is_soil_water_output, &
        is_irrigation_output, is_chem_in_all_soil_output, is_infiltration_at_depth_output, infiltration_point,                 &
        is_infiltrated_bottom_output, infiltration_point, is_infiltrated_bottom_output, extra_plots, temp_PLNAME, temp_chem_id,&
        temp_MODE, temp_ARG, temp_ARG2, temp_CONST, is_constant_profile, is_ramp_profile, ramp1, ramp2, ramp3, is_exp_profile, &
        exp_profile1, exp_profile2, is_total_degradation, is_app_window, app_window_span, app_window_step, is_timeseriesfile,  &
        is_waterbody_info_output, is_adjust_for_rain_schemes, rain_limit_schemes, optimum_application_window_schemes,          &
        intolerable_rain_window_schemes, min_days_between_apps_schemes, is_batch_scenario, scenario_batchfile,                 &
        is_needs_poundkg_conversion, open_water_adj, is_hydrolysis_override, is_output_spraydrift, is_absolute_year_schemes,   &
        is_gw_btc, runoff_mitigation_schemes, erosion_mitigation_schemes, drift_mitigation_schemes


    USE waterbody_PARAMETERs, ONLY: itsapond, itsareservoir, itsother, itstpezwpez, waterbody_names, USEPA_reservoir,          &
        USEPA_pond, use_tpezbuffer
    USE utilities

    INTEGER :: num_waterbodies, num_special_waterbodies

    INTEGER :: status, i, j
    INTEGER, PARAMETER :: Max_READ_line = 200
    CHARACTER (LEN=Max_READ_line) ::  wholeline
    INTEGER :: absolute_app_month
    INTEGER :: absolute_app_day
    INTEGER :: absolute_app_year
    INTEGER :: comma_1, comma_2, comma_3, comma_4, comma_5, comma_6, comma_7, comma_8
    INTEGER :: first_slash, second_slash
    INTEGER :: start_wb

    CHARACTER (LEN=256) :: scheme_name
    INTEGER             :: scheme_number_READin
    INTEGER             :: comma_place

    CHARACTER (LEN=10) :: dummy

    OPEN(Unit=inputfile_unit_number, FILE=(inputfile), STATUS='OLD', IOSTAT=status)
    IF (status .NE. 0) THEN
        WRITE(*,*)'Problem with PRZMVVWM input file: ', inputfile
        STOP
    ENDIF

    READ(inputfile_unit_number,*) !Version info                                                                     ! Line 1

    READ(inputfile_unit_number,'(A)')  working_directory                                                            ! Line 2
    WRITE(*,'(A23, A300)') ' working directory is: ', adjustl(working_directory)
    READ(inputfile_unit_number,'(A)') family_name                                                                   ! Line 3
    READ(inputfile_unit_number,'(A)') weatherfiledirectory                                                          ! Line 4
    WRITE(*,'(A23, A300)') ' weather directory is: ', adjustl(weatherfiledirectory)
    READ(inputfile_unit_number,*) open_water_adj                                                                    ! Line 5


    READ(inputfile_unit_number,*) is_koc, is_freundlich, is_nonequilibrium, is_needs_poundkg_conversion, is_hydrolysis_override ! Line 6
    READ(inputfile_unit_number,*) nchem                                                                             ! Line 7
    READ(inputfile_unit_number,*) k_f_input(1), k_f_input(2), k_f_input(3)                                          ! Line 8
    READ(inputfile_unit_number,*) N_f_input(1), N_f_input(2), N_f_input(3)                                          ! Line 9
    READ(inputfile_unit_number,*) k_f_2_input(1), k_f_2_input(2), k_f_2_input(3)                                    ! Line 10
    READ(inputfile_unit_number,*) N_f_2_input(1), N_f_2_input(2), N_f_2_input(3)                                    ! Line 11
    READ(inputfile_unit_number,*) K2(1), K2(2), K2(3)                                                               ! Line 12
    READ(inputfile_unit_number,*) lowest_conc, number_subdelt                                                       ! Line 13
    READ(inputfile_unit_number,*) water_column_halflife_input(1), &                                                 ! Line 14
        water_column_halflife_input(2), water_column_halflife_input(3), xAerobic(1), xAerobic(2)
    READ(inputfile_unit_number,*) water_column_ref_temp(1),  water_column_ref_temp(2),  water_column_ref_temp(3)    ! Line 15
    READ(inputfile_unit_number,*) benthic_halflife_input(1), benthic_halflife_input(2), benthic_halflife_input(3),& ! Line 16
        xBenthic(1), xbenthic(2)
    READ(inputfile_unit_number,*) benthic_ref_temp(1),     benthic_ref_temp(2),     benthic_ref_temp(3)             ! Line 17
    READ(inputfile_unit_number,*) photo_halflife_input(1), photo_halflife_input(2), photo_halflife_input(3),&       ! Line 18
        xPhoto(1), xphoto(2)
    READ(inputfile_unit_number,*) rflat(1), rflat(2), rflat(3)                                                      ! Line 19
    READ(inputfile_unit_number,*) hydrolysis_halflife_input(1), hydrolysis_halflife_input(2), &                     ! Line 20
        hydrolysis_halflife_input(3), xhydro(1), xhydro(2)


    READ(inputfile_unit_number,*) soil_degradation_halflife_input(1), soil_degradation_halflife_input(2),&          ! Line 21
        soil_degradation_halflife_input(3), xsoil(1), xsoil(2), is_total_degradation
    READ(inputfile_unit_number,*) soil_ref_temp(1), soil_ref_temp(2), soil_ref_temp(3)                              ! Line 22

    READ(inputfile_unit_number,*) foliar_halflife_input(1),  foliar_halflife_input(2),&                             ! Line 23
        foliar_halflife_input(3), foliar_formation_ratio_12, foliar_formation_ratio_23

    READ(inputfile_unit_number,*) plant_washoff_coeff(1), plant_washoff_coeff(2), plant_washoff_coeff(3)            ! Line 24
    READ(inputfile_unit_number,*) mwt(1), mwt(2), mwt(3)                                                            ! Line 25
    READ(inputfile_unit_number,*) vapor_press(1), vapor_press(2), vapor_press(3)                                    ! Line 26
    READ(inputfile_unit_number,*) solubilty(1), solubilty(2), solubilty(3)                                          ! Line 27
    READ(inputfile_unit_number,*) Henry_unitless(1), Henry_unitless(2), Henry_unitless(3)                           ! Line 28
    READ(inputfile_unit_number,*) DAIR(1), DAIR(2), DAIR(3)                                                         ! Line 29
    READ(inputfile_unit_number,*) Heat_of_Henry(1), Heat_of_Henry(2), Heat_of_Henry(3)                              ! Line 30
    READ(inputfile_unit_number,*) Q_10                                                                              ! Line 31
    READ(inputfile_unit_number,*) is_constant_profile                                                               ! Line 32
    READ(inputfile_unit_number,*) is_ramp_profile, ramp1, ramp2, ramp3                                              ! Line 33
    READ(inputfile_unit_number,*) is_exp_profile , exp_profile1, exp_profile2                                       ! Line 34
    READ(inputfile_unit_number,*) number_of_schemes                                                                 ! Line 35

    WRITE(*,*) "Number of schemes = ", number_of_schemes

    ! Lines become variable from here
    ALLOCATE (num_apps_in_schemes(number_of_schemes))
    ALLOCATE (app_reference_point_schemes(number_of_schemes))
    ALLOCATE (days_until_applied_schemes(number_of_schemes,366))
    ALLOCATE (is_absolute_year_schemes(number_of_schemes,366))

    ALLOCATE (application_rate_schemes(number_of_schemes,366))
    ALLOCATE (method_schemes(number_of_schemes,366))
    ALLOCATE (depth_schemes(number_of_schemes,366))
    ALLOCATE (split_schemes(number_of_schemes,366))
    ALLOCATE (drift_schemes(number_of_schemes,366))
    ALLOCATE (driftfactor_schemes(number_of_schemes,366))

    ALLOCATE (lag_schemes(number_of_schemes,366))
    ALLOCATE (periodicity_schemes(number_of_schemes,366))
    ALLOCATE (scenario_names(number_of_schemes,1000))
    ALLOCATE (number_of_scenarios(number_of_schemes))

    ALLOCATE (is_app_window(number_of_schemes))
    ALLOCATE (app_window_span(number_of_schemes))
    ALLOCATE (app_window_step(number_of_schemes))

    ALLOCATE (is_adjust_for_rain_schemes(number_of_schemes))
    ALLOCATE (rain_limit_schemes(number_of_schemes))
    ALLOCATE (optimum_application_window_schemes(number_of_schemes))
    ALLOCATE (intolerable_rain_window_schemes(number_of_schemes))
    ALLOCATE (min_days_between_apps_schemes(number_of_schemes))

    ALLOCATE (runoff_mitigation_schemes(number_of_schemes))
    ALLOCATE (erosion_mitigation_schemes(number_of_schemes))
    ALLOCATE (drift_mitigation_schemes(number_of_schemes))

    ALLOCATE (is_batch_scenario(number_of_schemes))
    ALLOCATE (scenario_batchfile(number_of_schemes))

    is_absolute_year_schemes = .FALSE.
    DO i=1, number_of_schemes
        READ(inputfile_unit_number,*) scheme_number_READin, scheme_name                                  ! scheme line 1  36
        !      WRITE(*,'(A22,I5,1X,A256)') " Scheme Number & Name ",scheme_number_READin, adjustl(scheme_name)

        READ(inputfile_unit_number,*) app_reference_point_schemes(i)                                     ! scheme line 2  37
        READ(inputfile_unit_number,*) num_apps_in_schemes(i)                                             ! scheme line 3  38

        DO j=1, num_apps_in_schemes(i)

            SELECT CASE (app_reference_point_schemes(i))
                CASE (0)
                    READ(inputfile_unit_number,'(A)')wholeline

                    first_slash  = index(wholeline, '/')
                    second_slash = index(wholeline, '/', .TRUE.)
                    comma_1      = index(wholeline, ',')

                    ! Scheme Line Group 4
                    IF (first_slash == second_slash) THEN  !this checks to see if only one "/" exists
                        ! only a day and month are given, so this is repeating every year

                        ! convert string dates into INTEGERs
                        READ(wholeline(1:first_slash-1),*)  absolute_app_month

                        READ(wholeline((first_slash+1): (comma_1 -1)),*)  absolute_app_day

                        ! Convert to julian day. Treat these similar to REALtive application, using Jan 1 as reference date

                        days_until_applied_schemes(i,j) = jd (1900, absolute_app_month,absolute_app_day)
                    ELSE
                        is_absolute_year_schemes(i,j) = .TRUE.

                        READ(wholeline(1:(first_slash-1)),*)  absolute_app_month
                        READ(wholeline((first_slash+1):(second_slash-1)),*)  absolute_app_day
                        READ(wholeline((second_slash+1): (comma_1 -1)),*)  absolute_app_year

                        !   WRITE(*,*) absolute_app_month, absolute_app_day, absolute_app_year

                        ! For now, treating just like absolute month /day.
                        ! THe absolute year can be used as a flag that is different than 1900 to indiucate absolute year

                        days_until_applied_schemes(i,j) = jd (absolute_app_year, absolute_app_month,absolute_app_day)


                        !the year is include, so this is year specific
                    END IF

                    comma_2 = comma_1 + index(wholeline((comma_1+1):LEN(wholeline)), ',')
                    comma_3 = comma_2 + index(wholeline((comma_2+1):LEN(wholeline)), ',')
                    comma_4 = comma_3 + index(wholeline((comma_3+1):LEN(wholeline)), ',')
                    comma_5 = comma_4 + index(wholeline((comma_4+1):LEN(wholeline)), ',')
                    comma_6 = comma_5 + index(wholeline((comma_5+1):LEN(wholeline)), ',')
                    comma_7 = comma_6 + index(wholeline((comma_6+1):LEN(wholeline)), ',')
                    comma_8 = comma_7 + index(wholeline((comma_7+1):LEN(wholeline)), ',')


                    READ(wholeline((comma_1+1):(comma_2-1)),*)         application_rate_schemes(i,j)
                    READ(wholeline((comma_2+1):(comma_3-1)),*)         method_schemes(i,j)
                    READ(wholeline((comma_3+1):(comma_4-1)),*)         depth_schemes(i,j)
                    READ(wholeline((comma_4+1):(comma_5-1)),*)         split_schemes(i,j)

                    READ(wholeline((comma_5+1):(comma_6-1)),*)         drift_schemes(i,j)

                    READ(wholeline((comma_6+1):(comma_7-1)),*)		   driftfactor_schemes(i,j)

                    READ(wholeline((comma_7+1):(comma_8-1)),*)         periodicity_schemes(i,j)

                    READ(wholeline((comma_8+1):LEN(wholeline)),*)      lag_schemes(i,j)
                    !READ(wholeline((comma_8+1):LEN(wholeline)),*)

                CASE DEFAULT
                    READ(inputfile_unit_number,*) days_until_applied_schemes(i,j), application_rate_schemes(i,j), &
                        method_schemes(i,j), depth_schemes(i,j), split_schemes(i,j),   &
                        drift_schemes(i,j), driftfactor_schemes(i,j) ,periodicity_schemes(i,j) , lag_schemes(i,j)    ! 39
            END SELECT
        END DO

        READ(inputfile_unit_number,*) is_app_window(i), app_window_span(i), app_window_step(i)                       ! Scheme Line 5  40
        IF (not(is_app_window(i) )) THEN
            app_window_span(i) =0  !the stepping starts with zero, so the span iz zero for only one iteration
            app_window_step(i)= 1
        END IF

        READ(inputfile_unit_number,*) is_adjust_for_rain_schemes(i),rain_limit_schemes(i), &                         ! Scheme Line 6  41
            optimum_application_window_schemes(i),intolerable_rain_window_schemes(i),min_days_between_apps_schemes(i)

        READ(inputfile_unit_number,*) number_of_scenarios(i)                                                         ! Scheme Line 7  42

        DO j=1, number_of_scenarios(i)
            READ(inputfile_unit_number,'(A512)')  scenario_names(i,j) !(i,j) i is scheme #, j is scenario #          ! Scheme Group of Lines 8  43
        END DO

        READ(inputfile_unit_number,*) is_batch_scenario(i)                                                           ! Scheme Line 9  44
        READ(inputfile_unit_number,'(A)') scenario_batchfile(i)                                                      ! Scheme Line 10  45

        !remove comma from end of file name, it will be last CHARACTER
        comma_place = index(scenario_batchfile(i),",",.TRUE.)
        IF (comma_place >0) THEN
            scenario_batchfile(i) = scenario_batchfile(i)(1: comma_place -1)
        END IF

        READ(inputfile_unit_number,*) ! "Mitigations"                                                                              Line 46
        READ(inputfile_unit_number,*)  runoff_mitigation_schemes(i) , erosion_mitigation_schemes(i),drift_mitigation_schemes(i)  ! Line 47

    ENDDO


    !Convert application rates to kg/ha if they are in lb/acre, 0.892179 kg/ha per lb/acre
    IF (is_needs_poundkg_conversion) THEN
        !   application_rate_schemes = 0.892179*application_rate_schemes          oops!
        application_rate_schemes = 1.12085*application_rate_schemes
    END IF


    READ(inputfile_unit_number,*) erflag                                                              ! WS Line 1  48
    READ(inputfile_unit_number,*)                                                                     ! WS Line 2  49
    READ(inputfile_unit_number,*)                                                                     ! WS Line 3  50
    READ(inputfile_unit_number,*)                                                                     ! WS Line 4  51
    READ(inputfile_unit_number,*)                                                                     ! WS Line 5  52
    READ(inputfile_unit_number,*)                                                                     ! WS Line 6  53
    READ(inputfile_unit_number,*) adjust_cn                                                           ! WS Line 7  54
    READ(inputfile_unit_number,*) itsapond,itsareservoir,itsother, itstpezwpez, use_tpezbuffer        ! WS Line 8  55
    READ(inputfile_unit_number,*) num_special_waterbodies                                             ! WS Line 9  56

    IF (.NOT. (itsapond .AND. itsareservoir)) start_wb = 0
    IF (itsapond .NEQV. itsareservoir) start_wb = 1
    IF (itsapond .AND. itsareservoir) start_wb = 2

    !The size of the array is the exact number of waterbodies, this SIZE is used later for setting  the run loop
    IF (itsother) THEN
        num_waterbodies = num_special_waterbodies + start_wb
        ALLOCATE (waterbody_names( num_waterbodies )) !add possibility of USEPA reservoir and pond.
    ELSE
        num_waterbodies = start_wb
        ALLOCATE (waterbody_names( num_waterbodies )) ! no special water bodies
    END IF


    IF (itsapond .AND. itsareservoir) THEN
        waterbody_names(1) = USEPA_pond
        waterbody_names(2) = USEPA_reservoir
    ELSE IF (itsapond .AND. .NOT. itsareservoir) THEN
        waterbody_names(1) = USEPA_pond
    ELSE IF (.NOT. itsapond .AND. itsareservoir) THEN
        waterbody_names(1) = USEPA_reservoir
    END IF


    DO i = 1 + start_wB,  num_special_waterbodies + start_wB
        IF (itsother) THEN
            READ(inputfile_unit_number,'(A)') waterbody_names(i)                                              ! WS Line 10  57
        ELSE !skip over the special water bodies that nmight be listed, and its other not checked
            READ(inputfile_unit_number,*)                                                                     ! WS Line 10 (alt)
        END IF
    END DO

    !READ in optional output to zts file
1   is_timeseriesfile = .FALSE.

    READ(inputfile_unit_number,*) is_runoff_output                                                ! OUTPUT Line 1  58
    IF (is_runoff_output ) is_timeseriesfile = .TRUE.
    READ(inputfile_unit_number,*) is_erosion_output                                               ! OUTPUT Line 2  59
    IF (is_erosion_output ) is_timeseriesfile = .TRUE.
    READ(inputfile_unit_number,*) is_runoff_chem_output                                           ! OUTPUT Line 3  60
    IF (is_runoff_chem_output) is_timeseriesfile = .TRUE.
    READ(inputfile_unit_number,*) is_erosion_chem_output                                          ! OUTPUT Line 4  61
    IF (is_erosion_chem_output) is_timeseriesfile = .TRUE.
    READ(inputfile_unit_number,*) is_conc_bottom_output                                           ! OUTPUT Line 5  62
    IF (is_conc_bottom_output) is_timeseriesfile = .TRUE.
    READ(inputfile_unit_number,*) is_daily_volatilized_output                                     ! OUTPUT Line 6  63
    IF (is_daily_volatilized_output) is_timeseriesfile = .TRUE.
    READ(inputfile_unit_number,*) is_daily_chem_leached_output, leachdepth                        ! OUTPUT Line 7  64
    IF (is_daily_chem_leached_output) is_timeseriesfile = .TRUE.
    READ(inputfile_unit_number,*) is_chem_decayed_part_of_soil_output, decay_start, decay_end     ! OUTPUT Line 8  65
    IF (is_chem_decayed_part_of_soil_output) is_timeseriesfile = .TRUE.
    READ(inputfile_unit_number,*) is_chem_in_all_soil_output                                      ! OUTPUT Line 9  66
    IF (is_chem_in_all_soil_output) is_timeseriesfile = .TRUE.
    READ(inputfile_unit_number,*) is_chem_in_part_soil_output, fieldmass_start, fieldmass_end     ! OUTPUT Line 10  67
    IF (is_chem_in_part_soil_output) is_timeseriesfile = .TRUE.
    READ(inputfile_unit_number,*)  is_chem_on_foliage_output                                      ! OUTPUT Line 11  68
    IF (is_chem_on_foliage_output) is_timeseriesfile = .TRUE.
    READ(inputfile_unit_number,*)  is_precipitation_output                                        ! OUTPUT Line 12  69
    IF (is_precipitation_output) is_timeseriesfile = .TRUE.
    READ(inputfile_unit_number,*)  is_evapotranspiration_output                                   ! OUTPUT Line 13  70
    IF (is_evapotranspiration_output) is_timeseriesfile = .TRUE.
    READ(inputfile_unit_number,*)  is_soil_water_output                                           ! OUTPUT Line 14  71
    IF (is_soil_water_output) is_timeseriesfile = .TRUE.
    READ(inputfile_unit_number,*)  is_irrigation_output                                           ! OUTPUT Line 15  72
    IF (is_irrigation_output) is_timeseriesfile = .TRUE.
    READ(inputfile_unit_number,*)  is_infiltration_at_depth_output,infiltration_point             ! OUTPUT Line 71  73
    IF (is_infiltration_at_depth_output) is_timeseriesfile = .TRUE.
    READ(inputfile_unit_number,*) is_infiltrated_bottom_output                                    ! OUTPUT Line 17  74
    IF (is_infiltrated_bottom_output) is_timeseriesfile = .TRUE.

    !these next 3 are for output from waterbody, this goes into a separate output file
    READ(inputfile_unit_number,*) is_waterbody_info_output   !waterbody depth, conc and benthic   ! OUTPUT Line 18  75

    READ(inputfile_unit_number,*)   is_output_spraydrift                                          ! OUTPUT Line 19  76
    READ(inputfile_unit_number,*)	is_gw_btc                                                     !             20  77
    IF (is_gw_btc) is_timeseriesfile = .TRUE.
    READ(inputfile_unit_number,*)	                                                              ! OUTPUT Line 21  78
    READ(inputfile_unit_number,*)                                                                 ! OUTPUT Line 22  79
    READ(inputfile_unit_number,*)	                                                              ! OUTPUT Line 23  80
    READ(inputfile_unit_number,*)	                                                              ! OUTPUT Line 24  81
    READ(inputfile_unit_number,*)	                                                              ! OUTPUT Line 25  82


    READ(inputfile_unit_number,*) extra_plots                                                     ! OUTPUT Line 26  83
    IF (extra_plots > 0) is_timeseriesfile = .TRUE.

    DO i = 1, extra_plots
        READ(inputfile_unit_number,*)  temp_PLNAME(i),  temp_chem_id(i), temp_MODE(i),temp_ARG(i),temp_ARG2(i),temp_CONST(i)   !OUTPUT Line 27
    END DO

    END SUBROUTINE READ_inputfile


    subroutine READ_scenario_file(schemenumber,scenarionumber, error)
    !READs scn2x files
    use constants_and_variables, ONLY:ScenarioFileUnit, scenario_names, &
        weatherfiLEName, latitude, num_crop_periods_input, &
        emd,emm,mad,mam,had,ham,max_root_depth,max_canopy_cover, max_number_crop_periods, &
        max_canopy_height, max_canopy_holdup,foliar_disposition, crop_lag, crop_periodicity  , PFAC,SFAC, min_evap_depth,&
        FLEACH,PCDEPL,max_irrig ,UserSpecifiesDepth, user_irrig_depth,irtype, USLEK,USLELS,USLEP, IREG,SLP, &
        nhoriz,thickness,bd_input,fc_input, wp_input, oc_input, bd_input, Num_delx,dispersion_input, &
        is_temperature_simulated , albedo, NUSLEC,GDUSLEC,GMUSLEC,cn_2, uslec, &
        runoff_extr_depth,runoff_decline,runoff_effic,erosion_depth, erosion_decline, erosion_effic,use_usleyears,Height_stagnant_air_layer_cm, &
        is_auto_profile,number_of_discrete_layers,  profile_thick, profile_number_increments,  evergreen,soil_temp_input, scenario_id, bottom_bc

    INTEGER :: eof
    LOGICAL, intent(out) :: error
    INTEGER, intent(in) :: schemenumber,scenarionumber
    CHARACTER (LEN=512) fiLEName
    INTEGER :: i,status


    !local that will likely need to go to module
    !LOGICAL :: evergreen

    CHARACTER(LEN= 50) ::dummy, dummy2
    REAL :: scalar_albedo, scaler_soil_temp  !these values are arrays in the program, but they are initialesd with constants

    !LOGICAL :: checkopen


    thickness  = 0.0
    bd_input   = 0.0
    fc_input   = 0.0
    wp_input   = 0.0
    oc_input   = 0.0

    error = .FALSE.
    fiLEName = trim(scenario_names(schemenumber,scenarionumber))

    !WRITE(*,*) "Operating on scenario: ", trim(fiLEName)

    !inquire(53, OPENED=checkopen)

    OPEN(Unit = ScenarioFileUnit, FILE=(fiLEName),STATUS='OLD', IOSTAT=status  )
    IF (status .NE. 0) THEN
        WRITE(*,*)'Problem opening scenario file: ', status, trim(fiLEName)
        Error = .TRUE.
        return
    ENDIF

    READ(ScenarioFileUnit,'(A)') scenario_id                     ! Line 1 Scenario !D
    READ(ScenarioFileUnit,'(A)') weatherfiLEName                 ! Line 2 Weather file name
    READ(ScenarioFileUnit,*,  IOSTAT=status ) latitude           ! Line 3 Latitude
    IF (status .NE. 0) THEN
        call scenario_error(error)  !set the error to true
        return
    END IF

    !the following comes from elsewhere:
    READ(ScenarioFileUnit,*)  ! Line 4 ignore get water type somewhere else  !A1,A2,A3,A4,A5, A6,A7
    READ(ScenarioFileUnit,*)  ! Line 5 UserSpecifiedFlowAvg.Checked, ReservoirFlowAvgDays.Text, CustomFlowAvgDays.Text)
    READ(ScenarioFileUnit,*)  ! Line 6 BurialButton.Checked default is always on in pwc 2020
    READ(ScenarioFileUnit,*)  ! Line 7 AFIELD
    READ(ScenarioFileUnit,*)  ! Line 8 waterAreaBox.Text)
    READ(ScenarioFileUnit,*)  ! Line 9 initialDepthBox.Text)
    READ(ScenarioFileUnit,*)  ! Line 10 maxDepthBox.Text)
    READ(ScenarioFileUnit,*)  ! Line 11 massXferBox.Text)
    READ(ScenarioFileUnit,*)  ! Line 12 calculate_prben.Checked, prbenBox.Text)
    READ(ScenarioFileUnit,*)  ! Line 13 benthicdepthBox.Text)
    READ(ScenarioFileUnit,*)  ! Line 14 porosityBox.Text)
    READ(ScenarioFileUnit,*)  ! Line 15 bdBox.Text)
    READ(ScenarioFileUnit,*)  ! Line 16 foc2Box.Text)
    READ(ScenarioFileUnit,*)  ! Line 17 DOC2Box.Text)
    READ(ScenarioFileUnit,*)  ! Line 18 biomass2Box.Text)
    READ(ScenarioFileUnit,*)  ! Line 19 dfacBox.Text)
    READ(ScenarioFileUnit,*)  ! Line 20 ssBox.Text)
    READ(ScenarioFileUnit,*)  ! Line 21 ChlorophyllBox.Text)
    READ(ScenarioFileUnit,*)  ! Line 22 foc1Box.Text)
    READ(ScenarioFileUnit,*)  ! Line 23 DOC1Box.Text)
    READ(ScenarioFileUnit,*)  ! Line 24 Biomass1Box.Text)
    READ(ScenarioFileUnit,*)  ! Line 25 vbNewLine & EpaDefaultsCheck.Checked                                                            'line 53
    READ(ScenarioFileUnit,*)  ! Line 26 String.Format("{0}{1},{2}", vbNewLine, ReservoirCroppedAreaBox.Text, CustomCroppedAreaBox.Text) 'Line 54
    READ(ScenarioFileUnit,*)  ! Line 27 vbNewLine

    READ(ScenarioFileUnit,'(A)') dummy !   msg = "******** start of PRZM information ******************" & vbNewLine

    READ(ScenarioFileUnit,*, IOSTAT=status) dummy, evergreen         ! Line 29, evergreen

    IF (status .NE. 0) THEN
        call scenario_error(error)
        return
    END IF
    READ(ScenarioFileUnit,*,  IOSTAT=status ) num_crop_periods_input ! Line 30, number of crop periods input
    IF (status .NE. 0) THEN
        call scenario_error(error)
        return
    END IF

    READ(ScenarioFileUnit,*) ! msg & String.Format("{0}{1},", vbNewLine, simpleRB.Checked)  ! Line 31 unused

    ! LInes 32 to 38 --------------------------------
    do  i=1,num_crop_periods_input
        READ(ScenarioFileUnit,*, IOSTAT = status) emd(i),emm(i) ,mad(i),mam(i),had(i),ham(i),max_root_depth(i),max_canopy_cover(i),  &
            max_canopy_height(i), max_canopy_holdup(i),foliar_disposition(i), crop_periodicity(i), crop_lag(i)   ! Line 32
        IF (status .NE. 0) THEN
            call scenario_error(error)
            return
        END IF

        !PWC saves canopy cover it as a percent, but przm needs fraction
        max_canopy_cover(i)=max_canopy_cover(i)/100.
    END DO

    ! the maximum rows are 7. If it is less than 7, use this do iteration to ignore these lines
    do i=1, max_number_crop_periods - num_crop_periods_input
        READ(ScenarioFileUnit,'(A)') dummy
    END DO
    !----------------------------------------------
    READ(ScenarioFileUnit,*)    ! Line 39 msg = msg & String.Format("{0}{1},{2},{3},{4}", vbNewLine, altRootDepth.Text, altCanopyCover.Text, altCanopyHeight.Text, altCanopyHoldup.Text)
    READ(ScenarioFileUnit,*)    ! Line 40 msg = msg & String.Format("{0}{1},{2},{3},{4}", vbNewLine, altEmergeDay.Text, altEmergeMon.Text, altDaysToMaturity.Text, altDaysToHarvest.Text)

    READ(ScenarioFileUnit,*, IOSTAT= status)  dummy,dummy, min_evap_depth ! Line 41, depth of soil evaporation, soil profile, first PARAMETER
    IF (status .NE. 0) THEN
        call scenario_error(error)
        return
    END IF

    READ(ScenarioFileUnit,*)                          ! Line 42  "*** irrigation information start ***"
    READ(ScenarioFileUnit,*) irtype                   ! Line 43 0 = none, 1 = overcanopy 2 = under canopy, Irrigation
    READ(ScenarioFileUnit,*) FLEACH,PCDEPL,max_irrig  ! Line 44, extra water fraction, allowed depletion, 0-0.6; max rate, cm

    READ(ScenarioFileUnit,*) UserSpecifiesDepth !, user_irrig_depth  ! UserSpecifiesIrrigDepth.Checked, IrrigationDepthUserSpec.Text)

    if (UserSpecifiesDepth) THEN
        backspace(ScenarioFileUnit)  !if true the userirrig depth will be blank sometinmes
        READ(ScenarioFileUnit,*) UserSpecifiesDepth, user_irrig_depth        ! Line 45 soil irrigation depth, if user specified or use root zone
    endif

    READ(ScenarioFileUnit,'(A)') !dummy "*** spare line for expansion"  ! Line 46

    READ(ScenarioFileUnit,'(A)') !dummy "*** spare line for expansion"  ! Line 47

    READ(ScenarioFileUnit,'(A)') !dummy "*** spare line for expansion"  ! Line 48

    READ(ScenarioFileUnit,*) USLEK,USLELS,USLEP            ! Line 49, runoff and erosion, this line and next line, the whole right 5 PARAMETERs
    READ(ScenarioFileUnit,*) IREG,SLP                      ! Line 50, runoff and erosion,
    READ(ScenarioFileUnit,*)                               ! Line 51   *** Horizon Info *******
    READ(ScenarioFileUnit,*) NHORIZ                        ! Line 52, number of soil layers
    READ(ScenarioFileUnit,*) (thickness(i), i=1, nhoriz)   ! Line 53, thickness
    READ(ScenarioFileUnit,*) (bd_input(i), i=1, nhoriz)    ! Line 54, bulk density
    READ(ScenarioFileUnit,*) (fc_input(i), i=1, nhoriz)    ! Line 55, max water fraction
    READ(ScenarioFileUnit,*) (wp_input(i), i=1, nhoriz)    ! Line 56, minimum water fraction
    READ(ScenarioFileUnit,*) (oc_input(i), i=1, nhoriz)    ! Line 57, organic carbon %
    READ(ScenarioFileUnit,*) (Num_delx(i), i=1, nhoriz)    ! Line 58, number of increment
    READ(ScenarioFileUnit,*) !(sand_input(i), i=1, nhoriz) ! Line 59
    READ(ScenarioFileUnit,*) !(clay_input(i), i=1, nhoriz) ! Line 60

    dispersion_input = 0.0


    READ(ScenarioFileUnit,*)   !*** Horizon End, Temperature Start ********  Line 61

    READ(ScenarioFileUnit,"(A50)", IOSTAT= status) dummy

    !msg = msg & String.Format("{0}{1},{2}", vbNewLine, albedoBox.Text, bcTemp.Text)
    READ(dummy,*, IOSTAT= status) scalar_albedo,scaler_soil_temp                     ! This sentence is not actually enforced
    IF (status .NE. 0) THEN  !this isn't populated for standard scnarios in 2023
        scalar_albedo = 0.2
        scaler_soil_temp = 15.
        !           WRITE(*,*) "Using Default albedo and gw temperature"
    END IF

    ALBEDO = scalar_albedo                 !albedo is monthly in przm
    soil_temp_input = scaler_soil_temp

    !need to set bottom boundary condition
    bottom_bc= scaler_soil_temp

    READ(ScenarioFileUnit,*) is_temperature_simulated !msg = msg & vbNewLine & simTemperature.Checked    ! Line 63

    READ(ScenarioFileUnit,*) !msg = msg & vbNewLine & "***spare line for expansion"                        Line 64
    READ(ScenarioFileUnit,*) !msg = msg & vbNewLine & "***spare line for expansion"                        Line 65
    READ(ScenarioFileUnit,*) !msg = msg & vbNewLine & "*** Erosion & Curve Number Info **********"         Line 66

    READ(ScenarioFileUnit,*) NUSLEC  !msg = msg & vbNewLine & NumberOfFactors.Text                         Line 67 runoff and erosion, row numbers
    READ(ScenarioFileUnit,*) (GDUSLEC(i), i=1,nuslec)                               ! Line 68 day
    READ(ScenarioFileUnit,*) (GMUSLEC(i), i=1,nuslec)                               ! Line 69 month
    READ(ScenarioFileUnit,*) (CN_2(i), i=1,nuslec)                                  ! Line 70 cn
    READ(ScenarioFileUnit,*) (USLEC(i), i=1,nuslec)                                 ! Line 71 usle-c
    READ(ScenarioFileUnit,*) !no longer use manning n mngn                          ! Line 72

    READ(ScenarioFileUnit,*) runoff_extr_depth,runoff_decline,runoff_effic            ! Line 73 runoff and erosion, lastly, left
    READ(ScenarioFileUnit,*)  erosion_depth, erosion_decline, erosion_effic           ! Line 74 runoff and erosion, lastly, right

    READ(ScenarioFileUnit,*)  use_usleyears                                           ! Line 75 if erosion is not changing every year

    READ(ScenarioFileUnit,*) !years of uslec need to add later if we want years       ! Line 76
    READ(ScenarioFileUnit,*) Height_stagnant_air_layer_cm                             ! Line 77 volatilization boundary layer (cm)
    !msg = msg & vbNewLine & volatilizationBoundaryBox.Text

    profile_thick= 0.0
    profile_number_increments=0
    READ(ScenarioFileUnit,*, IOSTAT=eof)  is_auto_profile                             ! Line 78 is autoprofile using PARAMETERs below and bypass value


    if (eof >= 0) THEN             !provides for the possibilty of older scenarios without this feature
        if (is_auto_profile) THEN
            READ(ScenarioFileUnit,*) number_of_discrete_layers
            do i = 1,  number_of_discrete_layers
                READ(ScenarioFileUnit,*) profile_thick(i), profile_number_increments(i)

            END DO
        END IF
    else
        is_auto_profile= .FALSE.  !for older files set default to No aurto profile
    END IF


    close (ScenarioFileUnit)
    end subroutine READ_scenario_file


    subroutine  scenario_error(error)
    implicit none
    LOGICAL, intent(out) :: error
    WRITE(*,*) 'Problem READing scenario file, skipping scenario ?????????????????????????????'
    error = .TRUE.
    end subroutine  scenario_error


    subroutine READ_batch_scenarios(batchfileunit, end_of_file, error_on_READ)
    use utilities_1
    use constants_and_variables, ONLY: scenario_id, weatherfiLEName,latitude, min_evap_depth, IREG, irtype,max_irrig, PCDEPL, fleach, USLEP, USLEK,USLELS,SLP, &
        num_crop_periods_input,NHORIZ, thickness,bd_input,fc_input,wp_input,oc_input, evergreen,emm, emd, mad, mam, had, ham, &
        PFAC,SFAC, max_canopy_cover, max_canopy_holdup, max_root_depth, crop_periodicity, crop_lag,UserSpecifiesDepth, is_temperature_simulated, &
        ALBEDO, soil_temp_input, dispersion_input, nuslec,cn_2, USLEC, gmuslec, gduslec, use_usleyears, &
        runoff_extr_depth,runoff_decline,runoff_effic,erosion_depth,erosion_decline,erosion_effic, Height_stagnant_air_layer_cm, &
        profile_thick, profile_number_increments, is_auto_profile,  number_of_discrete_layers, foliar_disposition, UserSpecifiesDepth,  &
        user_irrig_depth , bottom_bc

    LOGICAL, intent(out) :: end_of_file, error_on_READ
    INTEGER, intent(in)  :: batchfileunit

    CHARACTER(LEN=5) :: dummy
    INTEGER :: julian_emerg, julian_matur, julian_harv  !need conversion to days & months and put into emd(i),emm(i) ,mad(i),mam(i),had(i),ham(i
    REAL :: canopy_holdup      !needs to be put into array max_canopy_holdup(i)
    REAL :: canopy_coverage    !needs to be put into array max_canopy_cover(i)
    REAL :: root_depth         !needs to be put into array max_root_depth(i)
    REAL :: cn_cov, cn_fal, usle_c_cov, usle_c_fal
    INTEGER year !dummy not used but for sub op
    CHARACTER(LEN=50) :: weather_grid
    INTEGER :: iostatus
    INTEGER :: i
    REAL ::  gw_depth, gw_temp  !meters and C
    CHARACTER(LEN=1000) :: input_string
    CHARACTER(LEN=3) ::cropgroup,region

    thickness  = 0.0
    bd_input   = 0.0
    fc_input   = 0.0
    wp_input   = 0.0
    oc_input   = 0.0

    !Presets
    foliar_disposition = 1
    UserSpecifiesDepth = .false.
    user_irrig_depth = 0.0
    crop_periodicity(1) =1
    crop_lag(1) = 0
    !PFAC = 1.0                       !always using PET files now, No need to adjust
    !SFAC = 0.274                     !USDA value  now a PARAMETER
    ALBEDO = 0.2                      !albedo is monthly in przm
    UserSpecifiesDepth = .FALSE.      !only root irrigation
    is_temperature_simulated = .TRUE.
    dispersion_input = 0.0
    nuslec = 2
    use_usleyears = .FALSE.
    runoff_extr_depth = 8.0
    runoff_decline    = 1.4
    runoff_effic      = 0.19
    erosion_depth     = 0.1
    erosion_decline   = 0.0
    erosion_effic     = 1.0
    Height_stagnant_air_layer_cm = 5.0

    is_auto_profile = .TRUE.  !we will use discretizations specified independent of horizon info
    number_of_discrete_layers	=   6
    profile_thick(1) =  3.0
    profile_thick(2) =  7.0
    profile_thick(3) = 10.0
    profile_thick(4) = 80.0
    ! profile_thick(5) = gw_depth - 100.  ! gw_depth is depth to aquifer surface, subtract the 1 meter from above
    profile_thick(6) = 100.

    profile_number_increments(1) = 30
    profile_number_increments(2) =  7
    profile_number_increments(3) =  2
    profile_number_increments(4) =  4
    !  profile_number_increments(5) = int(gw_depth)/50 -2
    profile_number_increments(6) = 2

    end_of_file    = .FALSE.
    error_on_READ  = .FALSE.


    !READ in as string for better list directed error checking. Prevents READs to subsequent lines in the event of missing data
    READ(BatchFileUnit, '(A)',IOSTAT=iostatus ) input_string
    !Check for end of file

    if(IS_IOSTAT_END(iostatus)) THEN
        end_of_file = .TRUE.
        return
    else
        !   WRITE(*,'(A)') input_string
    END IF


    READ(input_string, *, IOSTAT=iostatus)   scenario_id, dummy, weather_grid, dummy ,dummy, dummy, dummy, dummy, dummy, &  !enter the long string of values
        latitude, dummy, min_evap_depth , IREG, irtype, max_irrig,  PCDEPL, FLEACH,dummy, julian_emerg, julian_matur, julian_harv, &
        dummy, dummy, dummy, dummy, canopy_holdup, canopy_coverage, root_depth, cn_cov, cn_fal, usle_c_cov, usle_c_fal, USLEP, USLEK,USLELS,SLP, &
        NHORIZ, thickness(1), thickness(2), thickness(3), thickness(4), thickness(5), thickness(6), thickness(7), thickness(8), &
        bd_input(1),bd_input(2),bd_input(3),bd_input(4),bd_input(5),bd_input(6),bd_input(7),bd_input(8), &
        fc_input(1),fc_input(2),fc_input(3),fc_input(4),fc_input(5),fc_input(6),fc_input(7),fc_input(8), &
        wp_input(1),wp_input(2),wp_input(3),wp_input(4),wp_input(5),wp_input(6),wp_input(7),wp_input(8), &
        oc_input(1),oc_input(2),oc_input(3),oc_input(4),oc_input(5),oc_input(6),oc_input(7),oc_input(8), &
        dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy, &
        dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,cropgroup,region, gw_depth, gw_temp


    if (iostatus /= 0 ) THEN  !there is a problem
        error_on_READ = .TRUE.
        WRITE (*,*) 'ERROR IN THIS SCENARIO (ignore if this is just the header)'
        return
    END IF


    scenario_id =  trim(scenario_id) //"_C" // trim(cropgroup) //"_R" //trim(region)

    !**************************************************
    !Foliar Disposition seems to be undefined in the batch file
    !canopy height is undefined, set to 1 meter by default put into: max_canopy_height(i)


    !added the following to enable complete fiLENames in weather_grid field in csv file, per ian
    if( INDEX(weather_grid, ".wea") > 0) THEN ! assume any file with this string is a complete weather file name
        weatherfiLEName = trim(adjustl(weather_grid))
    else
        weatherfiLEName = trim(adjustl(weather_grid)) // '_grid.wea'
    END IF


    num_crop_periods_input = 1
    evergreen = .FALSE.
    if (julian_emerg==0 .AND.  julian_matur==1 .AND. julian_harv==364) THEN
        evergreen = .TRUE.
    END IF

    max_canopy_cover(1) = canopy_coverage/100.  ! READ in as percent
    max_canopy_holdup(1)   = canopy_holdup
    max_root_depth(1)      = root_depth


    call get_date(julian_emerg, year, emm(1), emd(1) )

    call get_date(julian_matur, year, mam(1), mad(1) )

    call get_date(julian_harv , year,  ham(1), had(1) )


    !nuscle = 2 so only 2 curve numbers
    CN_2(1)  = cn_cov
    CN_2(2)  = cn_fal
    USLEC(1) = usle_c_cov
    USLEC(2) = usle_c_fal

    !only 1 crop
    GDUSLEC(1) = emd(1)  !emergence
    GMUSLEC(1) = emm(1)

    GDUSLEC(2) = had(1)  !harvest
    GMUSLEC(2) = ham(1)

    soil_temp_input  = gw_temp         !array for horizons, set as constant for all horizons as initial condition.
    bottom_bc = gw_temp                !bottom boundary is constant temperature



    profile_thick(5) = gw_depth - 100.  ! gw_depth is depth to aquifer surface, subtract the 1 meter (thickness of the horizons 1 to 4)
    profile_number_increments(5) = nint (profile_thick(5)/50.)

    WRITE(*,*) gw_depth, profile_thick(5), profile_number_increments(5)

    end subroutine READ_batch_scenarios


    subroutine READ_Weatherfile
    !open weather file, count file, allocate the weather variables, READ in dat and place in vectors to store in constant/variable mod.
    use utilities_1
    use constants_and_Variables, ONLY:precip, pet_evap,air_temperature, wind_speed, solar_radiation, num_records, &
        metfileunit,weatherfiLEName, weatherfiledirectory, &
        startday, first_year, last_year, num_years,first_mon, first_day





    INTEGER :: dummy, status,i, year
    CHARACTER(LEN=1024) :: weatherfile_pathandname


    weatherfile_pathandname = trim(adjustl(weatherfiledirectory)) // trim(adjustl(weatherfiLEName))

    OPEN(Unit = metfileunit,FILE=trim(adjustl(weatherfile_pathandname)),STATUS='OLD', IOSTAT=status)
    if (status .NE. 0) THEN
        WRITE(*,*)'Problem opening weather file. Name is ', trim(adjustl(weatherfiLEName))
        stop
    endif

    num_records = 0

    !READ first Date of Simulation
    READ (metfileunit,*, IOSTAT=status) first_mon, first_day, first_year

    if (status /= 0) THEN
        WRITE(*,*) "No data or other problem in Weather File"
        stop
    END IF

    startday = jd(first_year, first_mon, first_day)

    !Count the records in the weather file
    num_records=1
    do
        READ (metfileunit,*, IOSTAT=status) dummy
        if (status /= 0)exit
        num_records=num_records+1
    END DO

    !Allocate the weather PARAMETERs
    allocate (precip(num_records))
    allocate (pet_evap(num_records))
    allocate (air_temperature(num_records))
    allocate (wind_speed(num_records))
    allocate (solar_radiation(num_records))

    !ONLY READS WEA FORMAT
    !rewind and READ in weather data
    rewind(metfileunit)
    do i = 1, num_records
        READ(MetFileUnit,*, IOSTAT=status) dummy,dummy,year,precip(i),pet_evap(i),air_temperature(i),wind_speed(i),solar_radiation(i)
        if (status /= 0)THEN
            WRITE (*,*) "weather file problem on line ", i, status
            exit
        END IF
    END DO

    last_year = year

    num_years = last_year-first_year+1

    close(metfileunit)
    end subroutine READ_Weatherfile

    end module READinputs