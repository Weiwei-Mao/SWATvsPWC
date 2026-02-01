module clock_variables
!************************Clock Tinming  ***************************************    
REAL :: cputime_begin     !processor start reference time for begining of program
REAL :: time_1			  !cpu time at measured process
REAL :: time_2			  !cpu time at measured process
REAL :: time_3
REAL :: time_4
REAL :: time_5
REAL :: time_6

INTEGER :: c_count, c_rate, C_max  !used for system clock routine calls	


REAL :: clock_time_0   !clock start reference time for begining of program
REAL :: clock_time	   !clock time at process being measured

REAL :: Cumulative_cpu_1, Cumulative_cpu_2, Cumulative_cpu_3
REAL :: Cumulative_cpu_4 
REAL :: Cumulative_cpu_5  
REAL :: Cumulative_cpu_6 
REAL :: Cumulative_cpu_7 ,  Cumulative_cpu_8, Cumulative_cpu_9

contains
      subroutine time_check(message)
         implicit none                         
         CHARACTER(LEN = *) :: message

         CALL CPU_TIME (time_1)

         !Note to me: dont delete this output, delete the call to tine_check instead
         WRITE(*,*) message, time_1
         
      end subroutine 
      



end module clock_variables