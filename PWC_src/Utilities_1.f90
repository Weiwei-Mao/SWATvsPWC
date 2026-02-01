
	module utilities_1
    implicit none
    
	contains 
	
    function get_order(x) result(order)
    !returns an INTEGER array with indices idicating the sort order of original array x
       implicit none
       INTEGER, intent(in) :: x(:)
       INTEGER :: order(size(x))      ! Function result
    
       LOGICAL :: mask(size(x))
       INTEGER :: i

       mask = .true.
       do i = 1, size(x)
             order(i) = findloc( x, value=minval(x,dim=1,mask=mask), dim=1 )
             mask(order(i)) = .false.
       END DO  
    end function
    
    
    subroutine make_run_id (i,j, ii,mm)
    !makes a string that can be used for identifying output: Scheme#_Scenario#_ScenarioFiLEName (eg., 2_3_NDpumpkins)
    
    USE constants_and_variables, ONLY: run_id ,  scenario_id ,  full_run_identification, working_directory, family_name   
    USE waterbody_PARAMETERs, ONLY: waterbody_names
    implicit none
    INTEGER, intent(in) :: i,j, ii,mm
    CHARACTER(LEN=25) :: schemnumber, scenarionumber, appnumber
    INTEGER :: last_slash, last_dot,last_slash2, last_dot2
    
    CHARACTER(LEN=512) :: local_name
    
    !*****turn numbers into text*****
    WRITE(schemnumber,*) i
    WRITE(scenarionumber,*) j
    
    !last_slash = index(scenario_names(i,j), '\', .TRUE.)
    !last_dot   = index(scenario_names(i,j), '.', .TRUE.)
    
    last_slash2 = index(waterbody_names(ii), '\', .TRUE.)
    last_dot2   = index(waterbody_names(ii), '.', .TRUE.)
    
    IF (index(waterbody_names(ii), '.')==0) THEN
         local_name = waterbody_names(ii)
    else
         local_name= (waterbody_names(ii)((last_slash2+1):(last_dot2-1)))         
    END IF
    
    WRITE(appnumber, '(I4.4)') mm  !app window
    
    ! scheme_number_ScenarioName_WaterbodyName
    
    run_ID =  trim(adjustl(schemnumber)) // '_' //& 
    trim(scenario_id)  &
        // '_' // trim(adjustl(local_name)) // '_' // trim(adjustl(appnumber))
    
    !WRITE(*,'(A8, A256)') 'Run ID:', adjustl(run_ID )


    
    !run_ID =  trim(adjustl(schemnumber)) // '_' //&
    !trim(adjustl(scenario_names(i,j)((last_slash+1):(last_dot-1))))  &
    !    // '_' // trim(adjustl(local_name)) // '_' // trim(adjustl(appnumber))
    
    full_run_identification = trim(working_directory) // trim(family_name) // "_" // trim(run_id)
    

     
    end subroutine make_run_id 
    
    
    
    SUBROUTINE tridiagonal_solution(A,B,C,X,F,n)
        implicit none
        INTEGER,intent(in) :: n
        REAL,intent(in)    :: A(n), B(n),C(n) , F(n)
        REAL,intent(out)   :: X(n)
        
        REAL :: Y(n), ALPHA(n), BETA(n)
        INTEGER :: nu, i, j
        
        ALPHA(1) = B(1)
        BETA(1) = C(1)/ALPHA(1)
        Y(1) = F(1)/ALPHA(1)
        do I=2, N
            ALPHA(I) = B(I) - A(I)*BETA(I-1)
            BETA(I) = C(I)/ALPHA(I)
            Y(I) = (F(I)-A(I)*Y(I-1))/ALPHA(I)
        END DO
        

        
        
        X(N) = Y(N)
        NU=N-1

        do I=1, NU
            J=N-I
            X(J) = Y(J) - BETA(J)*X(J+1)          
        END DO
        

        
    end subroutine tridiagonal_solution
    

               
  !*****************************************************************************
   pure elemental INTEGER function jd (YEAR,MONTH,DAY)
     !calculate the days since 1/1/1900 given year,month, day, from Fliegel and van Flandern (1968)
    !Fliegel, H. F. and van Flandern, T. C. (1968). Communications of the ACM, Vol. 11, No. 10 (October, 1968). 

     implicit none
     INTEGER, intent(in) :: year,month,day

      JD= day-32075+1461*(year+4800+(month-14)/12)/4+367*(month-2-(month-14)/12*12) /12-3*((year+4900+(month-14)/12)/100)/4 -2415021

   end function jd
   !*****************************************************************************
   
   pure subroutine get_date (date1900, YEAR,MONTH,DAY)
 !computes THE GREGORIAN CALENDAR DATE (YEAR,MONTH,DAY) given days since 1900
   implicit none

   INTEGER,intent(out) :: YEAR,MONTH,DAY

   INTEGER,intent(in) :: date1900  !days since 1900
   INTEGER :: L,n,i,j

   L= 2483590 + date1900

   n= 4*L/146097

   L= L-(146097*n+3)/4
   I= 4000*(L+1)/1461001
   L= L-1461*I/4+31
   J= 80*L/2447

   day= L-2447*J/80
   L= J/11
   month = J+2-12*L
   year = 100*(N-49)+I+L

 !   YEAR= I
 !  MONTH= J
 !  DAY= K

   end subroutine get_date
      

     pure INTEGER function find_depth_node(n,depth,desired) 
     !Given an array "depth" of size "n" that is ordered from low to high values, this 
     !function will give the index of "depth" that is closest to the value "desired"
    
      implicit none
      INTEGER,intent(in)            :: n       !size of depth vector 
      REAL,DIMENSION(n), intent(in) :: depth   !vector holding incremental depths)
      REAL,intent(in)               :: desired !desired depth
      
      INTEGER :: i, index

      !WRITE(*,*) "CHECK for N herree:   n, depth, desired"
      !WRITE(*,*) n,depth,desired
      
      
      do i=1, n 
          index = i  !store value for the case where we go to the max n and i would be incremented another 1 value
          if (depth(i) > desired) exit
      END DO

      
      i = index     !if i falls out the loop above it will have a value of n+1, so we use index to capture REAL value
      if  (i==1) THEN 
          find_depth_node = 1
      ELSE IF (abs(depth(i) - desired) < abs (depth(i-1) - desired)) THEN
          find_depth_node = i
      else
          find_depth_node = i-1
      END IF

     end function find_depth_node  
           
     subroutine find_average_property(n, depth, target_depth, property, average)
     !weigted average, given vector of thicknesses and vector of correspnding property, Finds the average property
     !value of a target_depth
     use clock_variables
     implicit none
        INTEGER,intent(in) :: n                      !size of input vectors for depth and prperty
        REAL   ,intent(in) :: depth(n)          !vector of cumulative depths 
        REAL   ,intent(in) :: property(n)            !property of interest corresponding to thicknrss vector
        REAL   ,intent(in) :: target_depth           ! target depth is lower depth for averaging
        REAL   ,intent(out):: average                 !average property value from zero depth to target depth
        
        INTEGER :: i
        REAL ::  previous_depth, weighted_tally

        weighted_tally = 0.0
        previous_depth = 0.0
        
        do i = 1, n    
              if (depth(i) < target_depth) THEN
                  weighted_tally = weighted_tally+ property(i)*(depth(i) - previous_depth ) 
              elseif(depth(i) >= target_depth) THEN
      
                 weighted_tally = weighted_tally + property(i)*(target_depth -  previous_depth)

                 exit !exit do loop   
              END IF 
              previous_depth = depth(i)
        END DO
       average =  weighted_tally/target_depth     
       
       
     end subroutine
     

     
     subroutine find_medians(rows, columns, x , medians)
        !Find medians of each column in array x 
         INTEGER, intent(in) ::  rows, columns      
         REAL, intent(in),DIMENSION(:,:)    ::  x       !assumed shape array 
         REAL, intent(out),DIMENSION(:)     ::  medians
         
         REAL :: hold(rows) 
         
         INTEGER:: i,j
         INTEGER :: col
         REAL :: max
         
        do col = 1, columns !loop for each conc
           
            !arrange in accending order
            max = 0.0
            
            !x is assumed shape. also it is intent in, so need explicit bounds and local variable hold bcuz variable hold changes below          
            hold(1:rows) = x(1:rows, col) 
            
            do i=1,rows
               do j=i+1,rows
                  IF(hold(i)>hold(j))THEN
                     max=hold(i)
                     hold(i)=hold(j)
                     hold(j)=max
                  endif
               END DO
            END DO

            if(MOD(rows,2)==0)THEN
               medians(col) = (hold(rows/2)+hold(rows/2+1))/2.0
            else
               medians(col) = hold(rows/2 +1)
            endif
     END DO
     

     end subroutine find_medians
     
     
     
     subroutine pick_max (num_years,num_records,bounds,c, output)
!    !this subroutine choses the maximum values of subsets of the vector c
!    !the subsets are defined by the vector "bounds"
!    !maximum values of "c" are chosen from within the c indices defined by "bounds"
!    !output is delivered in the vector "output"
    implicit none
    INTEGER, intent(in) :: num_records
    INTEGER, intent(in) :: num_years
    INTEGER, intent(in) :: bounds(num_years)
    REAL, intent(in), DIMENSION(num_records) :: c
    REAL, intent(out),DIMENSION(num_years) :: output

    INTEGER :: i

    !forall (i = 1: num_years-1) output(i) = maxval( c(bounds(i):bounds(i+1)-1) ) changed 2/5/2020
    
    
    
    do concurrent(i = 1: num_years-1)
        output(i) = maxval( c(bounds(i):bounds(i+1)-1) )
    END DO
    
    output(num_years)= maxval( c(bounds(num_years):num_records) )

    
     end subroutine pick_max

     
     
     
 !***************************************************************
subroutine Return_Frequency_Value(returnfrequency, c_in, n, c_out, lowYearFlag)
    !CALCULATES THE Concentration at the given yearly return frequency
    implicit none
    
    REAL,intent(in) :: returnfrequency             !Example 1 in 10 years would be 10.0
    INTEGER,intent(in) :: n                        !number of items in list
    REAL, intent(in), DIMENSION(n):: c_in          !list of items
    
    REAL,intent(out):: c_out                       !output of 90th centile of peaks
    REAL:: f,DEC      
    INTEGER:: m    
    REAL,DIMENSION(n):: c_sorted
    LOGICAL, intent(out) :: LowYearFlag  !if n is less than 10, returns max value and LowYearFlag =1
    LowYearFlag = .false.

    call hpsort(n,c_sorted, c_in)  !returns a sorted array
    
    f = (1.0 -1.0/returnfrequency )*(n+1)
    m=int(f)
    DEC = f-m      
    
   if (n < returnfrequency)THEN
      c_out = c_sorted(n)
      LowYearFlag = .true.
   else 
    c_out = c_sorted(m)+DEC*(c_sorted(m+1)-c_sorted(m))
   END IF

end subroutine Return_Frequency_Value    
     
!****************************************************************
subroutine hpsort(n,ra,b)
!  from numerical recipes  (should be upgraded to new f90 routine)
    implicit none
    INTEGER,intent(in):: n
    REAL,intent(out),DIMENSION(n)::ra !ordered output array
    REAL,intent(in),DIMENSION(n):: b  !original unordered input array

    INTEGER i,ir,j,l
    REAL rra
    
    ra=b    ! this added to conserve original order

    if (n.lt.2) return

    l=n/2+1
    ir=n
10    continue
    if(l.gt.1)THEN 
    l=l-1
    rra=ra(l)
    else 
    rra=ra(ir)    
    ra(ir)=ra(1)
    ir=ir-1
    if(ir.eq.1)THEN 
    ra(1)=rra 
    return
    endif
    endif
    i=l 
    j=l+l
20    if(j.le.ir)THEN 
        if(j.lt.ir)THEN
            if(ra(j).lt.ra(j+1))j=j+1 
        endif
        if(rra.lt.ra(j))THEN 
            ra(i)=ra(j)
            i=j
            j=j+j
        else 
            j=ir+1
        endif
        goto 20
        endif
    ra(i)=rra
    goto 10
end subroutine hpsort
!*******************************************************************************     
     
!*******************************************************************************
subroutine find_first_annual_dates(num_years, first_annual_dates )
   use constants_and_variables, ONLY: first_year, first_mon, first_day, startday
   use utilities
   implicit none
   
   INTEGER,intent(in) :: num_years
   INTEGER,intent(out),DIMENSION(num_years) :: first_annual_dates
   INTEGER i

   do i = 1,num_years
      first_annual_dates(i) =  jd(first_year+(i-1), first_mon,first_day )    
   END DO

   first_annual_dates = first_annual_dates - startday+1

end subroutine find_first_annual_dates
!********************************************************************************


	subroutine find_in_table(row, column, table, tablerows,tablecolumns, output)

     !Given a REAL "table"  
     !with "tablerows" = number of rows of table
     !and "tablecolumns" = number of columns of table
     ! find the "output" in the table with input row index of INTEGER "row" and where
     ! the output is interpolated from the columns by the input REAL "column" value
     !In other words, rows are exact and columns are interpolated
    
	   INTEGER, intent(in) :: row      ! the method
	   REAL, intent(in)    :: column   !interpolate between column 
       
       INTEGER, intent(in) :: tablerows
       INTEGER, intent(in) :: tablecolumns      
       REAL, intent(in)    :: table(tablerows,tablecolumns)
	   REAL, intent(out)   :: output
	
	   INTEGER :: i
	   REAL    :: previous
	
	   !interpolate column value for values in row
	   previous = 0.0
	   do i = 1, tablecolumns
	   	 if (column == table(1, i)) THEN !exact match, get value and end
	   		 output =  table(row, i)
	   		 exit
	   	 elseif (table(1, i) < column ) THEN
	   		 previous = table(1, i)
	   	 else      !table(1, i))< column  !do interpolation and quit
	   		 output =  table(row, i-1) +   (table(row, i)-table(row, i-1)) *  (column - previous)/(table(1, i)- previous)
	   	!	 WRITE(*,'("            row = ",i2, " interpolate between columns ", i2, " and ", i2, ", fraction = ", g10.4 )') row, i-1, i,  (column - previous)/(table(1, i)- previous)
	   		 exit
		 END IF	 
       END DO


       
	
	end subroutine find_in_table





end module utilities_1
