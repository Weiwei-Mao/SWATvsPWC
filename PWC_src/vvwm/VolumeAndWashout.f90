module volumeAndwashout

contains

   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subroutine constant_volume_calc
        !This subroutine calculates volume and washout rate
        use utilities

        use constants_and_variables, ONLY: num_records, flowthru_the_body,k_flow,daily_depth,      &  !output array of daily wash out rates (per second)
                                     volume1,     &  !for this case, v=constant, but keep array to be consistent with program
                                     Daily_avg_flow_out
        
        use waterbody_PARAMETERs, ONLY: area_waterbody,depth_0,flow_averaging
        
        implicit none
        REAL,DIMENSION(num_records):: q_avg        !Adjustable-day flow average
        REAL:: v_0                                 ! water body volume
                
        Daily_avg_flow_out = 0.0
               
        v_0 = area_waterbody*depth_0

        if (flow_averaging ==0) THEN
              q_avg = sum(flowthru_the_body)/num_records 
        else
             call window_average(flowthru_the_body,flow_averaging,num_records,q_avg) !calculate 30-day previous average  
        endif 

        k_flow = q_avg/v_0   !array of washout rates
        volume1= v_0
        daily_depth = depth_0    
        Daily_avg_flow_out = sum(q_avg)/num_records    !used for output CHARACTERization only
        
    end subroutine constant_volume_calc
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subroutine volume_calc
        !Full VVWM, This subroutine calculates volume and washout rate
        
        use constants_and_variables, ONLY: num_records, evap_m, precip_m, DELT_vvwm,minimum_depth,flowthru_the_body,&
            daily_depth,volume1,k_flow ,Daily_avg_flow_out
        
        
        use waterbody_PARAMETERs, ONLY: depth_0, depth_max,area_waterbody

        implicit none
        INTEGER:: day
        REAL:: v_0                                  !initial water body volume [m3]
        REAL:: v_max                                !maximum water body volume [m3]
        REAL:: v_min                                !minimum water body volume [m3]
        REAL:: v_previous
        REAL:: check
        REAL,DIMENSION(num_records)::vol_net
        REAL,DIMENSION(num_records)::evap_area
        REAL,DIMENSION(num_records)::precip_area
        
        
        INTEGER :: i
        
        Daily_avg_flow_out = 0.0
        
        v_0 = area_waterbody*depth_0
        v_max = area_waterbody*depth_max
        v_min = area_waterbody*minimum_depth
        k_flow = 0.        !sets all values of the array to zero
        v_previous = v_0

        precip_area = precip_m*area_waterbody /86400.    !m3/s
        evap_area = evap_m*area_waterbody /86400.    !m3/s, evap factor pfac now calculated in convert_weatherdata_for_VVWM  

        vol_net = (flowthru_the_body-evap_area+precip_area)*DELT_vvwm  !volume of water added in day; whole array operations
        
        
       
        do day = 1,num_records
             
            check = v_previous + vol_net(day)
            if (check > v_max) THEN
                volume1(day) = v_max
                k_flow(day) = (check-v_max)/DELT_vvwm/v_max   !day # and washout VOLUME
            ELSE IF (check < v_min) THEN
                volume1(day) = v_min
            else
                volume1(day) = check
            END IF
            v_previous = volume1(day)
          END DO

        
        
               
        Daily_avg_flow_out = sum(k_flow)*v_max/num_records  !used for output CHARACTERization only
        daily_depth = volume1/area_waterbody !whole array operation

        
        
        
    end subroutine volume_calc
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subroutine burial_calc(koc)
	!this is not REALly burial,mut a mass sediment mass balance onthe benthic zone

     use waterbody_PARAMETERs, ONLY: froc2
     use constants_and_variables, ONLY:   burial,capacity_2, & 
                                          k_burial                !output(kg/sec)
                              
        implicit none
        REAL,intent(in) :: koc
        REAL :: kd_sed_2  

        kd_sed_2 = KOC*FROC2*.001       !Kd of sediment  [m3/kg]
        k_burial=  burial* kd_sed_2/capacity_2

    end subroutine burial_calc
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



end module volumeAndwashout