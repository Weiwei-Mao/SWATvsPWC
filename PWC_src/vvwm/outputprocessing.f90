module outputprocessing

    contains
    

    subroutine output_processor(chem_index, output_unit, unit_number,unit_number_deg1,unit_number_deg2,&
                                summary_fiLEName, summary_fiLEName_deg1, summary_fiLEName_deg2, waterbody_name )

    use utilities
    use utilities_1, ONLY: pick_max, find_first_annual_dates
    use waterbody_PARAMETERs, ONLY: baseflow,SimTypeFlag, zero_depth, is_zero_depth, Afield
    
    use constants_and_variables, ONLY:  num_records, run_id, is_hed_files_made,is_add_return_frequency, additional_return_frequency, &
                                       num_years, startday, &
                                 gamma_1,        &
                                 gamma_2,        &
                                 fw1,            &
                                 fw2,            &
                                 aq1_store,      &   !beginning day after app concentration in water column
                                 aq2_store,      &
                                 aqconc_avg1,    &   !average daily concentration (after app)
                                 aqconc_avg2,    &
                                 daily_depth,    &
                                 runoff_total ,  &
                                 erosion_total,  &
                                 spray_total ,   &
                                 Daily_Avg_Runoff, Daily_avg_flow_out,  runoff_fraction, erosion_fraction, drift_fraction ,&
    k_burial, k_aer_aq, k_flow, k_hydro, k_photo, k_volatile,k_anaer_aq, gamma_1, gamma_2, gw_peak, post_bt_avg ,throughputs,simulation_avg, &
	is_waterbody_info_output, full_run_identification, applied_mass_sum_gram_per_cm2 , fraction_off_field
   
                         
    implicit none
    INTEGER,             intent(in)    :: chem_index
    INTEGER,             intent(in)    :: output_unit                                             !time series
    INTEGER,             intent(in)    :: unit_number,unit_number_deg1,unit_number_deg2            ! summary files
    CHARACTER(LEN= 500), intent(in)    :: summary_fiLEName, summary_fiLEName_deg1, summary_fiLEName_deg2
    !LOGICAL,             intent(inout) :: First_time_through
    CHARACTER(LEN= 20), intent(in)     :: waterbody_name
    
    
    CHARACTER(LEN=512) :: waterbody_outputfile
    
    !temporary PARAMETERs for esa, should make this more general in the future
    REAL:: return_frequency

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    REAL:: simulation_average
    
    REAL :: xxx  !local variable

 !   REAL(8),DIMENSION(num_records)::c1
 !   REAL(8),DIMENSION(num_records)::cavgw


    REAL,DIMENSION(num_records):: c4        !all 4-day averaged values
    REAL,DIMENSION(num_records):: c21
    REAL,DIMENSION(num_records):: c60
    REAL,DIMENSION(num_records):: c90
    REAL,DIMENSION(num_records):: c365
    REAL,DIMENSION(num_records):: benthic_c21
    
 !   REAL,DIMENSION(num_years):: onedayavg
    REAL,DIMENSION(num_years):: peak   
    REAL,DIMENSION(num_years):: benthic_peak  !now it 1-day average
    
    
    REAL,DIMENSION(num_years):: c1_max  !peak year to year daily average
    
    REAL,DIMENSION(num_years):: c4_max  !the peak 4-day average within the 365 days after application
    REAL,DIMENSION(num_years):: c21_max
    REAL,DIMENSION(num_years):: c60_max
    REAL, DIMENSION(num_years)::c90_max
    REAL,DIMENSION(num_years):: c365_max !the peak 365-day average within the 365 days after application
                                            !last year will be short depending on application date
                                        
    REAL,DIMENSION(num_years):: benthic_c21_max                         
                                        
    INTEGER :: i    
    INTEGER ::date_time(8)
    REAL :: convert                    !conversion factor kg/m3 

    REAL :: Total_Mass
    INTEGER :: YEAR,MONTH,DAY
    INTEGER :: eliminate_year
    INTEGER,DIMENSION(num_years) ::  first_annual_dates !array of yearly first dates (absolute days).
                                    ! First date is the caLENdar day of start of simulation 
    first_annual_dates= 0


 
if (is_waterbody_info_output) THEN
	select case (chem_index)
		
	    case (1)
		    waterbody_outputfile = trim(full_run_identification) // '_parent_'       // trim(waterbody_name) // '.out'
	    case (2)
		    waterbody_outputfile = trim(full_run_identification) // '_daughter_'      // trim(waterbody_name) // '.out'
	    case (3)
		    waterbody_outputfile = trim(full_run_identification) // '_granddaugter_' // trim(waterbody_name) // '.out'	
	    case default
		    waterbody_outputfile =trim(full_run_identification) // '_nada.out'	
        end select

    
	open (UNIT=output_unit,FILE= trim(waterbody_outputfile),  STATUS='unknown')
    
    

    
     WRITE(output_unit,'(A42)') 'Depth(m), Water Col(kg/m3), Benthic(kg/m3)'
 
    do i =1, num_records
        WRITE(output_unit,'(G12.4E3, "," ,ES12.4E3, "," ,ES12.4E3, "," ,ES12.4E3)')  daily_depth(i), aqconc_avg1(i), aqconc_avg2(i)
	END DO
	close (output_unit)
END IF


    !For Certain water bodies, users want to exclude concentrations below a certain level   
    if (is_zero_depth) THEN
          where (daily_depth < zero_depth) aqconc_avg1 = 0.0
    END IF

!
!    !***** EFED TIME SERIES ********************************************
!    WRITE(12,*) "col 1: Daily Depth (m)" 
!    WRITE(12,*) "col 2: Daily Average Aqueous Conc. in Water Columm, kg/m3 "
!    WRITE(12,*) "col 3: Daily Average Aqueous Conc. in Benthic Zone, kg/m3"
!    WRITE(12,*) "col 4: Daily Peak Aqueous Conc. in Water Column, kg/m3"
!    
!    WRITE(12,*)  startday, "= Start Day (number of days since 1900)"
!    do i =1, num_records
!         WRITE(12,'(G12.4E3, "," ,ES12.4E3, "," ,ES12.4E3, "," ,ES12.4E3)')  daily_depth(i), aqconc_avg1(i), aqconc_avg2(i), aq1_store(i) 
!    END DO
!    
!    !******HED TIME SERIES **************************************
!    if (SimTypeFlag /=2  .and. is_hed_files_made) THEN  !No need for CaLENdex and DEEM for Pond
!        !************************  DEEM Input File Header ********************************
!        WRITE(22,*) "'DEEM Input File" 
!        WRITE(22,'(" ''Performed on: ", i2,"/",i2,"/",i4,2x,"at " ,i2,":",i2) ') &
!                       date_time(2),date_time(3),date_time(1),date_time(5),date_time(6)
!        WRITE(22,*) "'Scenario: ", trim(scenario_id)," " ,trim( waterbodytext)
!        WRITE(22,*) "'List of Daily Average Aqueous Concentration in Water Columm, mg/L "
!        WRITE(22,*) "'Start data:"
!        
!        !************************  CALENDEX & DEEM DATA Input ********************************  
!        call get_date (startday, YEAR,MONTH,DAY)
!        eliminate_year = year
!        
!        do i =1, num_records
!             call get_date (startday+i-1, YEAR,MONTH,DAY)
!            if (year > eliminate_year) THEN 
!                WRITE(22,'(ES15.4)')  aqconc_avg1(i)*1000.0   !DEEM FILE DATA
!            !    WRITE(23,'(I2.2,"/",I2.2,"/", I4, "," ,ES15.4E3)')  month,day, year, aqconc_avg1(i)*1000.0   !CaLENdex File Data   
!            END IF    
!        END DO
!    END IF
!    !**************************************************************************
!END IF

    !Calculate chronic values *******************
    !The following returns the n-day running averages for each day in simulation
   
    call window_average(aqconc_avg1,4,num_records,c4)
    call window_average(aqconc_avg1,21,num_records,c21)
    call window_average(aqconc_avg1,60,num_records,c60)
    call window_average(aqconc_avg1,90,num_records,c90)  
    call window_average(aqconc_avg1,365,num_records,c365)

    Simulation_average = sum(aqconc_avg1)/num_records
    
    call window_average(aqconc_avg2,21,num_records,benthic_c21)

    call find_first_annual_dates (num_years, first_annual_dates )


	 
    call pick_max(num_years,num_records, first_annual_dates,aqconc_avg1,c1_max)     !NEW FIND DAILY AVERAGE CONCENTRATION RETURN   
    call pick_max(num_years,num_records, first_annual_dates,c4,c4_max)
    call pick_max(num_years,num_records, first_annual_dates,c21,c21_max)
    call pick_max(num_years,num_records, first_annual_dates,c60,c60_max)
    call pick_max(num_years,num_records, first_annual_dates,c90,c90_max)
    call pick_max(num_years,num_records, first_annual_dates,benthic_c21,benthic_c21_max)

    !treat the 365 day average somewhat differently:
    !In this case, we simply are calculating the average for the 365 day forward from the
    !day of application
 
    do concurrent (i=1:num_years-1) 
        c365_max(i) = c365(first_annual_dates(i)+365)  
    END DO
    
    c365_max(num_years) = c365(num_records)

    !****Calculate Acute values *******************
    call pick_max(num_years,num_records, first_annual_dates,aq1_store, peak)    !now using 1-d avg for acutes intead of peak
   ! call pick_max(num_years,num_records, first_annual_dates,aq2_store, benthic_peak)
    
    call pick_max(num_years,num_records, first_annual_dates,aqconc_avg2, benthic_peak) !new 4/5/2023 Changed to daily average for benthin acute
    

    !*****************************************

    convert = 1000000.
    peak               = peak*convert
    c1_max             = c1_max*convert
    c4_max             = c4_max*convert
    c21_max            = c21_max*convert
    c60_max            = c60_max*convert
    c90_max            = c90_max*convert
    c365_max           = c365_max*convert
    benthic_peak       = benthic_peak*convert
    benthic_c21_max    = benthic_c21_max*convert
    Simulation_average = Simulation_average*convert


 !   if (is_output_all== .false.) THEN  !append an output file for a reduced output batch run
       return_frequency = 10.0

 
        Total_Mass = runoff_total(chem_index) + erosion_total(chem_index) +  spray_total(chem_index)        !kg i think
        If (Total_Mass <= 0.0) THEN
            runoff_fraction  = 0.0
            erosion_fraction = 0.0
            drift_fraction   = 0.0
            fraction_off_field =0.0
        else
            runoff_fraction  = runoff_total(chem_index) /Total_Mass
            erosion_fraction = erosion_total(chem_index) /Total_Mass
            drift_fraction   =  spray_total(chem_index)/Total_Mass
            
            if (applied_mass_sum_gram_per_cm2 > 0.0) THEN
                  fraction_off_field = Total_Mass/(applied_mass_sum_gram_per_cm2*Afield*10.)  !applied mass is in kg/ha, afield is in m2
            else
                  fraction_off_field = 0.0
            endif
        END IF
        
    
       !WRITE(*,*) "total and fraction off field" , applied_mass_sum_gram_per_cm2*Afield*10., fraction_off_field
       !WRITE(*,*) 'Doing output process'  
       call calculate_effective_halflives()
       

       call WRITE_simple_batch_data(chem_index, unit_number,unit_number_deg1,unit_number_deg2,&
                                     summary_fiLEName, summary_fiLEName_deg1, summary_fiLEName_deg2, &
                                     return_frequency,num_years, peak,Simulation_average,c1_max,c4_max, &
                                     c21_max,c60_max,c90_max,c365_max,benthic_peak, benthic_c21_max )      


    end subroutine output_processor



!***************************************************************
!subroutine pick_max (num_years,num_records,bounds,c, output)
!!    !this subroutine choses the maximum values of subsets of the vector c
!!    !the subsets are defined by the vector "bounds"
!!    !maximum values of "c" are chosen from within the c indices defined by "bounds"
!!    !output is delivered in the vector "output"
!    implicit none
!    INTEGER, intent(in) :: num_records
!    INTEGER, intent(in) :: num_years
!    INTEGER, intent(in) :: bounds(num_years)
!    REAL, intent(in), DIMENSION(num_records) :: c
!    REAL, intent(out),DIMENSION(num_years) :: output
!
!    INTEGER :: i
!
!    !forall (i = 1: num_years-1) output(i) = maxval( c(bounds(i):bounds(i+1)-1) ) changed 2/5/2020
!    
!    
!    
!    do concurrent(i = 1: num_years-1)
!        output(i) = maxval( c(bounds(i):bounds(i+1)-1) )
!    END DO
!    
!    output(num_years)= maxval( c(bounds(num_years):num_records) )
!
!    
!end subroutine pick_max


!***************************************************************

!!*******************************************************************************
!subroutine find_first_annual_dates(num_years, first_annual_dates )
!   use constants_and_variables, ONLY: first_year, first_mon, first_day, startday
!   use utilities
!   implicit none
!   
!   INTEGER,intent(in) :: num_years
!   INTEGER,intent(out),DIMENSION(num_years) :: first_annual_dates
!   INTEGER i
!
!   do i = 1,num_years
!      first_annual_dates(i) =  jd(first_year+(i-1), first_mon,first_day )    
!   END DO
!
!   first_annual_dates = first_annual_dates - startday+1
!
!end subroutine find_first_annual_dates
!!********************************************************************************



 !subroutine WRITE_returnfrequency_data(return_frequency, unit_number,num_years, peak,Simulation_average,c1_max,c4_max,c21_max,c60_max,c90_max,c365_max,benthic_peak, benthic_c21_max   )
 !  use constants_and_variables, Only : version_number,Sediment_conversion_factor,fw2 
 !  use utilities_1, ONLY: Return_Frequency_Value
 !  
 !  implicit none
 ! 
 !    REAL, intent(in)                         :: return_frequency
 !    INTEGER, intent(in)                       :: unit_number
 !    INTEGER, intent(in)                       :: num_years
 !    REAL   , intent(in), DIMENSION(num_years) :: peak,c1_max,c4_max,c21_max,c60_max,c90_max,c365_max,benthic_peak, benthic_c21_max
 !    REAL, intent(in)                          :: Simulation_average
 !    
 !    
 !    REAL :: peak_out,c1_out, c4_out,c21_out,c60_out,c90_out,c365_out,benthic_peak_out,benthic_c21_out
 !    LOGICAL :: lowyearflag
 !    INTEGER :: date_time(8)
 !    INTEGER :: i
 !    CHARACTER(LEN=10) :: frequencystring 
 !
 ! 
 !    !**find values corresponding to  percentiles
 !    call Return_Frequency_Value(return_frequency, peak,           num_years, peak_out,         lowyearflag)
 !    
 !    call Return_Frequency_Value(return_frequency, c1_max,         num_years, c1_out,           lowyearflag)
 !    
 !    call Return_Frequency_Value(return_frequency, c4_max,         num_years, c4_out,           lowyearflag)
 !    call Return_Frequency_Value(return_frequency, c21_max,        num_years, c21_out,          lowyearflag)
 !    call Return_Frequency_Value(return_frequency, c60_max,        num_years, c60_out,          lowyearflag)
 !    call Return_Frequency_Value(return_frequency, c90_max,        num_years, c90_out,          lowyearflag)
 !    call Return_Frequency_Value(return_frequency, c365_max,       num_years, c365_out,         lowyearflag)
 !    call Return_Frequency_Value(return_frequency, benthic_peak,   num_years, benthic_peak_out, lowyearflag)
 !    call Return_Frequency_Value(return_frequency, benthic_c21_max,num_years, benthic_c21_out,  lowyearflag)
 !   
 !
 !   call date_and_time(VALUES = date_time)
 !   
 !   WRITE(unit_number,*) "Waterbody Output, PRZMVVWM Version ", Version_Number
 ! 
 !   WRITE(unit_number,*) 
 !   WRITE(unit_number,*)  '*******************************************'
 !   WRITE(unit_number,'("Performed on: ", i2,"/",i2,"/",i4,2x,"at " ,i2,":",i2) ') &
 !                  date_time(2),date_time(3),date_time(1),date_time(5),date_time(6)
 ! 
 !   WRITE(frequencystring, '(F4.1)')  return_frequency  !internal WRITE in order to right justify string
 ! 
 !
 ! 
 ! if (LowYearFlag)THEN
 !     WRITE(unit_number,  '(''* Insufficient years to calculate 1-in-'', A5, ''. Only maximums are reported here.'') '   ) adjustl(frequencystring)
 ! else 
 !      WRITE(unit_number,*)
 ! END IF
 !
 !! WRITE(unit_number, '(''Initial 1-in-'', A5 ,'' = '', G14.4E3,'' ppb'')'  )    adjustl(frequencystring), peak_out
 ! WRITE(unit_number, '(''1-d avg   1-in-'', A5 ,'' = '', G14.4E3,'' ppb'')'  )    adjustl(frequencystring),c1_out
 ! WRITE(unit_number, '(''365-d avg 1-in-'', A5 ,'' = '', G14.4E3,'' ppb'')'  )    adjustl(frequencystring),c365_out 
 ! WRITE(unit_number, '(''Simulation Avg       = '', G14.4E3,'' ppb'')'  )         simulation_average
 ! WRITE(unit_number, '(''4-d avg  1-in-'', A5 ,''  = '', G14.4E3,'' ppb'')'  )  adjustl(frequencystring),c4_out
 ! WRITE(unit_number, '(''21-d avg 1-in-'', A5 ,''  = '', G14.4E3,'' ppb'')'  )  adjustl(frequencystring), c21_out
 ! WRITE(unit_number, '(''60-d avg 1-in-'', A5 ,''  = '', G14.4E3,'' ppb'')'  )  adjustl(frequencystring), c60_out
 ! WRITE(unit_number, '(''90-d avg 1-in-'', A5 ,''  = '', G14.4E3,'' ppb'')'  )  adjustl(frequencystring),c90_out
 ! 
 ! 
 !   
 ! WRITE(unit_number, '(''Benthic Pore Water 1-d  avg 1-in-'', A5 ,''= '', G14.4E3,'' ppb'')'  ) adjustl(frequencystring),  benthic_peak_out
 ! WRITE(unit_number, '(''Benthic Pore Water 21-d avg 1-in-'', A5 ,''= '', G14.4E3,'' ppb'')'  ) adjustl(frequencystring),  benthic_c21_out
 !   
 ! WRITE(unit_number, '(''Benthic Conversion Factor             = '', G14.4E3,'' -Pore water (ug/L) to (total mass, ug)/(dry sed mass,kg)'')')&
 !       Sediment_conversion_factor*1000.
 ! 
 ! WRITE(unit_number, '(''Benthic Mass Fraction in Pore Water   = '', G14.4E3)'  )  fw2
 ! WRITE (unit_number,*)
 ! WRITE (unit_number,'(A)') 'YEAR    1-day     4-day      21-day     60-day     90-day   Yearly Avg Benthic Pk  Benthic 21-day'
 ! do I=1, num_years  
 !  WRITE(unit_number,'(I3,1x,8ES11.2E3)') i,c1_max(i),c4_max(i),c21_max(i),c60_max(i),c90_max(i),c365_max(i),benthic_peak(i),benthic_c21_max(i)
 ! END DO
 ! 
 !
 !end subroutine WRITE_returnfrequency_data
 




subroutine WRITE_simple_batch_data(chem_index,unit_number,unit_number_deg1,unit_number_deg2, &
                                   summary_fiLEName, summary_fiLEName_deg1, summary_fiLEName_deg2, &
                                   return_frequency,num_years, peak,Simulation_average,c1_max, &
                                   c4_max,c21_max,c60_max,c90_max,c365_max,benthic_peak, benthic_c21_max   )

use constants_and_variables, ONLY: run_id,Sediment_conversion_factor,fw2 ,&
    nchem,     runoff_fraction,erosion_fraction,drift_fraction,summary_outputfile, &
    effective_washout, effective_watercol_metab, effective_hydrolysis, effective_photolysis, effective_volatization, effective_total_deg1,&
    effective_burial, effective_benthic_metab, effective_benthic_hydrolysis, effective_total_deg2, &
    gw_peak, post_bt_avg ,throughputs,simulation_avg, fraction_off_field, family_name, app_window_counter, &
    hold_for_medians_wb, hold_for_medians_daughter,hold_for_medians_grandaughter, First_time_through_wb, working_directory

use waterbody_PARAMETERs, ONLY: FROC2

use utilities_1, ONLY: Return_Frequency_Value


    implicit none   
    INTEGER, intent(in)                       :: num_years
    REAL, intent(in)                          :: return_frequency
    REAL, intent(in), DIMENSION(num_years)    :: peak,c1_max,c4_max,c21_max,c60_max,c90_max,c365_max,benthic_peak, benthic_c21_max
    REAL, intent(in)                          :: Simulation_average
    !LOGICAL, intent(inout)                    :: First_time_through
    
    
    INTEGER, intent(in) :: unit_number,unit_number_deg1,unit_number_deg2
    CHARACTER(LEN= 500), intent(in):: summary_fiLEName, summary_fiLEName_deg1, summary_fiLEName_deg2
    
    
    INTEGER, intent(in) ::chem_index
    CHARACTER (LEN=457) :: header
    
    
    !****LOCAL*********************
    REAL      :: peak_out,c1_out, c4_out,c21_out,c60_out,c90_out,c365_out,benthic_peak_out,benthic_c21_out    
    LOGICAL   :: lowyearflag
    CHARACTER(LEN= 257) :: local_run_id
    
    If (First_time_through_wb(1) .AND. chem_index ==1 ) THEN
        header = 'Run Information                                                                  ,      1-d avg,    365-d avg,    Total avg,      4-d avg,     21-d avg,     60-d avg,      B 1-day,   B 21-d avg,    Off-Field,  Runoff Frac,   Erosn Frac,   Drift Frac,  col washout,    col metab,    col hydro,    col photo,    col volat,    col total,  ben sed rem,    ben metab,    ben hydro,    ben total,      gw_peak,  post_bt_avg,   throughput,   sim_avg_gw'

        
        Open(unit=unit_number,FILE=  (trim(working_directory) //trim(family_name) // "_" // trim(summary_fiLEName)),Status='unknown')  
            WRITE(unit_number, '(''Benthic Conversion Factor             = '', G14.4E3,'' -Pore water (ug/L) to (total mass, ug)/(dry sed mass,kg)'')') Sediment_conversion_factor(1)*1000.
            WRITE(unit_number, '(''OC sediment fraction                  = '', G14.4E3)') froc2
            WRITE(unit_number, '(A457)') header
        First_time_through_wb(1) = .FALSE.
    END IF             

        If (First_time_through_wb(2) .AND. chem_index ==2 ) THEN
        header = 'Run Information                                                                  ,      1-d avg,    365-d avg,    Total avg,      4-d avg,     21-d avg,     60-d avg,      B 1-day,   B 21-d avg,    Off-Field,  Runoff Frac,   Erosn Frac,   Drift Frac,  col washout,    col metab,    col hydro,    col photo,    col volat,    col total,  ben sed rem,    ben metab,    ben hydro,    ben total,      gw_peak,  post_bt_avg,   throughput,   sim_avg_gw'
        
            Open(unit=unit_number_deg1,FILE= (trim(working_directory) //trim(family_name) // "_" // trim(summary_fiLEName_deg1)),Status='unknown')  
            WRITE(unit_number_deg1, '(''Benthic Conversion Factor             = '', G14.4E3,'' -Pore water (ug/L) to (total mass, ug)/(dry sed mass,kg)'')') Sediment_conversion_factor(2)*1000.
            WRITE(unit_number_deg1, '(''OC sediment fraction                  = '', G14.4E3)') froc2
            WRITE(unit_number_deg1, '(A457)') header
            First_time_through_wb(2) = .FALSE.
        END IF
        
        If (First_time_through_wb(3) .AND. chem_index ==3 ) THEN
        header = 'Run Information                                                                  ,      1-d avg,    365-d avg,    Total avg,      4-d avg,     21-d avg,     60-d avg,      B 1-day,   B 21-d avg,    Off-Field,  Runoff Frac,   Erosn Frac,   Drift Frac,  col washout,    col metab,    col hydro,    col photo,    col volat,    col total,  ben sed rem,    ben metab,    ben hydro,    ben total,      gw_peak,  post_bt_avg,   throughput,   sim_avg_gw'
 
            Open(unit=unit_number_deg2,FILE= (trim(working_directory) //trim(family_name) // "_" // trim(summary_fiLEName_deg2)),Status='unknown')
            WRITE(unit_number_deg2, '(''Benthic Conversion Factor             = '', G14.4E3,'' -Pore water (ug/L) to (total mass, ug)/(dry sed mass,kg)'')') Sediment_conversion_factor(3)*1000.
            WRITE(unit_number_deg2, '(''OC sediment fraction                  = '', G14.4E3)') froc2
            WRITE(unit_number_deg2, '(A457)') header
            
            First_time_through_wb(3) = .FALSE.
            
        END IF
        

   
    
    !**find values corresponding to  percentiles
    call Return_Frequency_Value(return_frequency, peak,           num_years, peak_out,         lowyearflag)   
    call Return_Frequency_Value(return_frequency, c1_max,         num_years, c1_out,           lowyearflag)
    
    call Return_Frequency_Value(return_frequency, c4_max,         num_years, c4_out,           lowyearflag)
    call Return_Frequency_Value(return_frequency, c21_max,        num_years, c21_out,          lowyearflag)
    call Return_Frequency_Value(return_frequency, c60_max,        num_years, c60_out,          lowyearflag)
    !call Return_Frequency_Value(return_frequency, c90_max,        num_years, c90_out,          lowyearflag)
    call Return_Frequency_Value(return_frequency, c365_max,       num_years, c365_out,         lowyearflag)
    call Return_Frequency_Value(return_frequency, benthic_peak,   num_years, benthic_peak_out, lowyearflag)
    call Return_Frequency_Value(return_frequency, benthic_c21_max,num_years, benthic_c21_out,  lowyearflag)

	
    
    select case (chem_index)
    case (1)
        local_run_id = trim(run_id) // '_Parent'
        WRITE(unit_number,'(A80,1x,26(",", ES13.4E3))') (adjustl(local_run_id)), c1_out, c365_out , simulation_average, c4_out, c21_out,c60_out,benthic_peak_out, benthic_c21_out, fraction_off_field, runoff_fraction,erosion_fraction,drift_fraction, &
        effective_washout, effective_watercol_metab, effective_hydrolysis, effective_photolysis, effective_volatization, effective_total_deg1, effective_burial, effective_benthic_metab, effective_benthic_hydrolysis, effective_total_deg2, gw_peak(1), post_bt_avg(1) ,throughputs(1),simulation_avg(1)    !effective_total_deg2 does not mean degradate, means benthic

 !**capture data for median calculations here
       hold_for_medians_wb( app_window_counter,1 )= c1_out
       hold_for_medians_wb( app_window_counter,2 )= c365_out
       hold_for_medians_wb( app_window_counter,3 )= simulation_average
       hold_for_medians_wb( app_window_counter,4 )= c4_out
       hold_for_medians_wb( app_window_counter,5 )= c21_out
       hold_for_medians_wb( app_window_counter,6 )= c60_out
       hold_for_medians_wb( app_window_counter,7 )= benthic_peak_out
       hold_for_medians_wb( app_window_counter,8 )= benthic_c21_out
       hold_for_medians_wb( app_window_counter,9 )= post_bt_avg(1)
       hold_for_medians_wb( app_window_counter,10)= throughputs(1)
       hold_for_medians_wb( app_window_counter,11)=  gw_peak(1)
    case (2)
        local_run_id = trim(run_id) // '_deg1'
        WRITE(unit_number_deg1,'(A80,1x,26(",", ES13.4E3))') (adjustl(local_run_id)), c1_out, c365_out , simulation_average, c4_out, c21_out,c60_out,benthic_peak_out, benthic_c21_out,fraction_off_field,runoff_fraction,erosion_fraction,drift_fraction, &
        effective_washout, effective_watercol_metab, effective_hydrolysis, effective_photolysis, effective_volatization, effective_total_deg1, effective_burial, effective_benthic_metab, effective_benthic_hydrolysis, effective_total_deg2, gw_peak(2), post_bt_avg(2) ,throughputs(2) ,simulation_avg(2)

       hold_for_medians_daughter( app_window_counter ,1  )= c1_out
       hold_for_medians_daughter( app_window_counter ,2  )= c365_out
       hold_for_medians_daughter( app_window_counter ,3  )= simulation_average
       hold_for_medians_daughter( app_window_counter ,4  )= c4_out
       hold_for_medians_daughter( app_window_counter ,5  )= c21_out
       hold_for_medians_daughter( app_window_counter ,6  )= c60_out
       hold_for_medians_daughter( app_window_counter ,7  )= benthic_peak_out
       hold_for_medians_daughter( app_window_counter ,8  )= benthic_c21_out
       hold_for_medians_daughter( app_window_counter ,9  )= post_bt_avg(2)
       hold_for_medians_daughter( app_window_counter ,10 )= throughputs(2)  
       hold_for_medians_daughter( app_window_counter ,11 )=  gw_peak(2)
    case (3)
        local_run_id = trim(run_id)  // '_deg2'
        WRITE(unit_number_deg2,'(A80,1x,26(",", ES13.4E3))')(adjustl(local_run_id)), c1_out, c365_out , simulation_average, c4_out, c21_out,c60_out,benthic_peak_out, benthic_c21_out,fraction_off_field, runoff_fraction,erosion_fraction,drift_fraction, &
        effective_washout, effective_watercol_metab, effective_hydrolysis, effective_photolysis, effective_volatization, effective_total_deg1, effective_burial, effective_benthic_metab, effective_benthic_hydrolysis, effective_total_deg2, gw_peak(3), post_bt_avg(3) ,throughputs(3) ,simulation_avg(3)

        hold_for_medians_grandaughter(  app_window_counter,1   )= c1_out
        hold_for_medians_grandaughter(  app_window_counter,2   )= c365_out
        hold_for_medians_grandaughter(  app_window_counter,3   )= simulation_average
        hold_for_medians_grandaughter(  app_window_counter,4   )= c4_out
        hold_for_medians_grandaughter(  app_window_counter,5   )= c21_out
        hold_for_medians_grandaughter(  app_window_counter,6   )= c60_out
        hold_for_medians_grandaughter(  app_window_counter,7   )= benthic_peak_out
        hold_for_medians_grandaughter(  app_window_counter,8   )= benthic_c21_out
        hold_for_medians_grandaughter(  app_window_counter,9   )= post_bt_avg(3)
        hold_for_medians_grandaughter(  app_window_counter,10  )= throughputs(3)  
        hold_for_medians_grandaughter(  app_window_counter,11  )=  gw_peak(3)

        case default
    end select

     
end subroutine WRITE_simple_batch_data



















Subroutine calculate_effective_halflives()
!Effective compartment halflives averaged over simulation duration:
use constants_and_variables, ONLY: num_records, k_burial, k_aer_aq, k_flow, k_hydro, k_photo, k_volatile,k_anaer_aq, gamma_1, gamma_2, fw1, fw2, &
    effective_washout, effective_watercol_metab, effective_hydrolysis, effective_photolysis, effective_volatization, effective_total_deg1, effective_burial, effective_benthic_metab, effective_benthic_hydrolysis, effective_total_deg2

implicit none
REAL :: xxx

        xxx = sum(k_flow)
        if (xxx > 0.) THEN
            effective_washout= 0.69314/(xxx/num_records)/86400
        else
            effective_washout = 0.
        END IF                      
       
        xxx = sum(k_aer_aq)
        if (xxx > 0.) THEN
            effective_watercol_metab =  0.69314/(xxx/num_records)/86400
        else
            effective_watercol_metab = 0.0
        END IF
       
        xxx = sum(k_hydro*fw1)  
        if (xxx > 0.) THEN
            effective_hydrolysis = 0.69314/(xxx/num_records)/86400
        else
            effective_hydrolysis= 0.0
        END IF
       
        xxx =  sum(k_photo*fw1)
        if (xxx > 0.) THEN
            effective_photolysis =  0.69314/(xxx/num_records)/86400
        else
            effective_photolysis = 0.0
        END IF       
       
        xxx =  sum(k_volatile*fw1)
        if (xxx > 0.) THEN
           effective_volatization =  0.69314/(xxx/num_records)/86400
        else
           effective_volatization = 0.0
        END IF      
        
        xxx =  sum(gamma_1)
        if (xxx > 0.) THEN
            effective_total_deg1 =  0.69314/(xxx/num_records)/86400
        else
            effective_total_deg1 = 0.0
        END IF    
        
       
        xxx = sum(k_burial)
        if (xxx > 0.) THEN
             effective_burial =  0.69314/(xxx/num_records)/86400
        else
             effective_burial = 0.0
        END IF
       
        xxx = sum(k_anaer_aq)
        if (xxx > 0.) THEN
           effective_benthic_metab =  0.69314/(xxx/num_records)/86400
        else
           effective_benthic_metab = 0.0
        END IF
        
        xxx = sum(k_hydro*fw2)
        if (xxx > 0.) THEN
           effective_benthic_hydrolysis =  0.69314/(xxx/num_records)/86400
        else
           effective_benthic_hydrolysis = 0.0
        END IF 
       
       
        xxx = sum(gamma_2)
        if (xxx > 0.) THEN
            effective_total_deg2 = 0.69314/(xxx/num_records)/86400
        else
            effective_total_deg2 = 0.0
        END IF 


end Subroutine calculate_effective_halflives





     
     
     
     
end module outputprocessing