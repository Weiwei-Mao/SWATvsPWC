PROGRAM PRZMVVWM

    USE allocations
    USE READinputs
    USE constants_and_variables, ONLY: maxFileLENgth, inputfile,number_of_schemes,  &
        number_of_scenarios,  First_time_through_wb, First_time_through_wpez,       &
        First_time_through_tpez,First_time_through_medians, app_window_span,        &
        app_window_step, application_date, application_date_original,               &
        is_adjust_for_rain, is_batch_scenario, scenario_batchfile ,                 &
        BatchFileUnit, run_id, app_window_counter, First_time_through_medians_wpez, &
        First_time_through_medians_tpez

    USE waterbody_PARAMETERs, ONLY: READ_waterbodyfile, get_pond_PARAMETERs,        &
        get_reservoir_PARAMETERs,waterbody_names,USEPA_reservoir,USEPA_pond,        &
        spraytable,itstpezwpez
    USE clock_variables

    USE PRZM_VERSION
    ! use PRZM_part
    USE initialization
    USE VVWM_solution_setup
    USE schemeload
    USE utilities_1
    USE plantgrowth
    USE field_hydrology
    USE chemical_transport
    USE Output_From_Field
    USE Pesticide_Applications
    !   use READbatchscenario
    USE TPEZ_WPEZ
    USE process_medians, ONLY: calculate_medians
    USE TPEZ_spray_initialization, ONLY:tpez_drift_finalize, set_tpez_spray
    IMPLICIT NONE

    INTEGER :: LENgth               ! LENgth of input file CHARACTERs
    INTEGER :: hh, i ,jj, kk,  iostatus
    LOGICAL :: error
    LOGICAL :: end_of_file, error_on_READ
    ! TRUE if USEPA_pond has just been simulated AND TPEZ/WPEZ run has been requested (i.e., istpezwpez==treue)
    LOGICAL :: run_tpez_wpez
    CHARACTER :: dummy
    CHARACTER(LEN=20) ::message

    !################################################
    CALL CPU_TIME (cputime_begin)

    CALL SYSTEM_CLOCK(c_count, c_rate, c_max)
    clock_time_0 = REAL(c_count)/REAL(c_rate)
    Cumulative_cpu_3 = 0.0

    !################################################

    First_time_through_wb   = .TRUE.  !ARRAY of 3 parent, deg1, deg2, used to WRITE output file headers
    First_time_through_tpez = .TRUE.
    First_time_through_wpez = .TRUE.

    First_time_through_medians_wpez = .TRUE.
    First_time_through_medians_tpez = .TRUE.

    CALL get_command_argument(1,inputfile,LENgth)
    inputfile = 'PRZMVVWM.txt'                       ! Local used for testing, modified 2/24/2025
    CALL przm_id                                     ! Stamp the runstatus file

    CALL READ_inputfile                              ! READ input file, PRZMVVWM.txt
    CALL chemical_manipulations                      ! daughter and granddaughter mole transfer

    ! LOOP Description
    ! Outermost (hh) is for info describing watershed and waterbody, At present, there is a dependency
    ! on the watershed size and the PER/AREA output from the field. This means that a field output
    ! cannot be used for multiple waterbodies.
    ! A field must be run for each waterbody. Future work should look at removing this
    ! over-PARAMETERization of the erosion routine,
    ! and that would allow a single PRZM run for all waterbodies
    ! The dependency occurs in the hydraulic LENgth which is area dependent AND in the MUSS, MUSL equations.
    ! Middle Loop (i) is for Application Schemes
    ! Inner Loop (j) is for scenarios

    DO hh = 1, size(waterbody_names)             ! loop based on waterbody numbers

        run_tpez_wpez = .FALSE.
        ! we want a separate file for each median (move this uop if all in one file, but file has waterbody name
        First_time_through_medians = .TRUE.

        ! note: waterbody spray drift table is loaded here (tpez will need different table)
        SELECT CASE   (waterbody_names(hh))
        CASE (USEPA_reservoir)
            CALL get_reservoir_PARAMETERs
        CASE (USEPA_pond)
            CALL get_pond_PARAMETERs
            IF (itstpezwpez) THEN
                run_tpez_wpez = .TRUE.
            END IF
        CASE DEFAULT
            CALL READ_waterbodyfile(hh)
        END SELECT

        WRITE(*,*) 'Doing waterbody: ',  trim(waterbody_names(hh))
        DO i = 1, number_of_schemes
            WRITE(*,*) 'Doing Scheme No. ', i

            ! *******************************************************************
            ! will use spray table in here --need to make spray table correct if changed by tpez
            ! gets the individual application scheme from the whole scheme table, non scenario specfic
            CALL set_chemical_applications(i)

            ! If running TPEZ, set the unchanging (with scheme) PARAMETERs here like spraydrift
            CALL set_tpez_spray(i)
            !*******************************************************************

            IF (is_batch_scenario(i)) THEN
                WRITE(*,'("Batch Scenario File: ", A200) ')  adjustl( scenario_batchfile(i))
                OPEN (Unit = BatchFileUnit, FILE=scenario_batchfile(i),STATUS='OLD', IOSTAT= iostatus )
                IF (iostatus /= 0) THEN
                    WRITE(*,*) "Can't find scenario batch file for scheme ",i
                    STOP
                END IF

                READ(BatchFileUnit,*) dummy  ! skip header
                end_of_file = .FALSE. ! reset the batch scenario READing
                error_on_READ = .FALSE.
            END IF

            kk=0 !index for scenario loop

            DO  ! scenario do loop
                ! Loop controled by either the number of files in the batch or by the number of scenarios READ in from input file
                kk=kk+1
                IF ( is_batch_scenario(i)) THEN
                    CALL READ_batch_scenarios(BatchFileUnit, end_of_file, error_on_READ)

                    ! *****ERROR and EOF CHECKING ON SCENARIO BATCH FILE ****************
                    IF (end_of_file) THEN
                        CLOSE(BatchFileUnit)
                        WRITE(*,*) 'end of batch scenario READ, exit after completing the run'
                        EXIT  !exit scenario loop
                    END IF

                    IF (error_on_READ) THEN
                        WRITE(*,*) 'bad scenario # ', kk
                        CYCLE
                    END IF
                    ! *****END ERROR CHECKING ON SCENARIO BATCH FILE ****************

                ELSE    ! *******use scenarios directly READ into input
                    IF (kk == number_of_scenarios(i) + 1) EXIT  !end of scenario list from gui inputs
                    CALL READ_scenario_file(i,kk, error)

                    IF (error) THEN
                        WRITE(*,*) 'exiting scenario loop due to scenario open problem'
                        CYCLE   ! error, try next scenario in scheme
                    END IF
                END IF

                CALL READ_Weatherfile ! this READs the new format weather file

                CALL INITL    ! initialize and ALLOCATIONS przm variables  (also sets application and drift)
                ! set scenario-dependent TPEZ drift here

                CALL Crop_Growth
                CALL hydrology_only
                CALL allocation_for_VVWM

                CALL tpez_drift_finalize !must be called after INITL because need to know total applications

                app_window_counter = 0  !use this to track app window to find medians

                DO jj = 0, app_window_span(i), app_window_step(i)
                    CALL reset_initial_masses

                    application_date= application_date_original + jj
                    app_window_counter = app_window_counter +1
                    ! makes a string that can be used for identifying output scheme#_scenario#_scenarioname
                    CALL make_run_id (i,kk, hh,jj)

                    ! "Rain Fast" Option
                    IF (is_adjust_for_rain) CALL adjust_application_dates_for_weather

                    CALL chem_transport_onfield
                    CALL groundwater

                    CALL VVWM

                    IF (run_tpez_wpez) THEN ! only do TPEZ WPEZ if its a pond run
                        CALL wpez     ! drift is same as pond noi need to recalculate
                        CALL tpez(i)  ! need to send in scheme number to find drift
                    END IF
                END DO

                CALL scenario_hydrolgy_summary

                ! process application date window into medians
                CALL calculate_medians(app_window_counter,run_tpez_wpez )

                CALL deallocate_scenario_PARAMETERs

            END DO  !END SCENARIO LOOP kk, begin at Line 125

            ! allocations are done in set_chmical_applications, need to deallocatte for next scheme
            CALL deallocate_application_PARAMETERs

            WRITE (message, '(A17,I3)') 'cpu time, scheme ', i
            CALL time_check(message)

        END DO  !End scheme Loop, i, begin at Line 98

        DEALLOCATE (spraytable )

    END DO ! End Waterbody/Watershed Loop, begin at Line 78

    !*******************************************
    CALL time_check('End program cpu time ')
    CALL SYSTEM_CLOCK(c_count, c_rate, c_max)
    clock_time = REAL(c_count)/REAL(c_rate)
    WRITE (*,*) 'Total clock time = ',clock_time  - clock_time_0

END PROGRAM PRZMVVWM
