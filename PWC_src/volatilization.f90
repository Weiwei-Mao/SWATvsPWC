module volatilization
    implicit none
    
         ! The wind measurements in the meteoroLOGICAL Daily Values File (*.dvf)
         ! were normalized to 10 meters, Open Flat Terrain.
         
         REAL,PARAMETER :: uWind_z0 = 0.03
         REAL,PARAMETER :: uWind_D  = 0.0
            
         !REAL, PARAMETER :: vonKarman = 0.4

         REAL, PARAMETER :: Minimum_Canopy_Height_cm = 5.0
         REAL, PARAMETER :: Minimum_Canopy_Height_m = Minimum_Canopy_Height_cm * 1.0E-02
    
    contains
    
 
    
    
    subroutine volatilization_setup(k)
     use  constants_and_Variables, ONLY: HEIGHT,uWind_Reference_Height, wind, CONDUC, CONDUC,UBT,air_TEMP, &
           is_temperature_simulated,Henry_unitless, DAIR, Height_stagnant_air_layer_cm, &
           thair_new,  dair,DGAIR, theta_sat
     
        implicit none
        REAL     :: CNDBDY(3)     ! CNDBDY: Boundary Layer's Conductance (cm day^-1) = 1/Rdb, (3 chemicals)
        INTEGER, intent(in) :: k ! chemical number
        REAL :: ATEMP(2),ZCH, z0, urh, uch,PWIND(2),TOTCR,d

    
        
        ! DGAIR now includes correction for air to bulk volume
        DGAIR  =  (THAIR_new**(10./3.)/theta_sat**2)*DAIR(K) * THAIR_new 

       !***********************************************************  
   

       ! Rdb:    Boundary layer resistance (day/cm) = d / Dair,
       ! d:      Height of the stagnant air layer above the soil
       ! Dair:   Molecular Diffusivity in air (cm^2 day^-1)
        CNDBDY(K) = DAIR(K) / Height_stagnant_air_layer_cm

       ! CONDUC: Canopy Conductance Including Boundary Layer's Conductance (cm day^-1)
       CONDUC(K) = CNDBDY(K)

       ! When canopy develops, resistance type approach is used to estimate the volatilization             
       ! flux and concentration retains in the canopy

       ! HEIGHT: Canopy height (cm)


       IF (HEIGHT .GT. Minimum_Canopy_Height_cm) THEN
           ZCH = HEIGHT/100.0     ! convert to meter
           
           
           IF (is_temperature_simulated) THEN
                 ATEMP(1)= UBT
           ELSE  
                ATEMP(1)= 15.0
           ENDIF 
           
           
           !IF (ITFLAG .EQ. 0) THEN
           !     ATEMP(1)= 15.0
           !ELSE
           !     ATEMP(1)= UBT
           !ENDIF
           
           ATEMP(2)= air_TEMP

           ! Let u_2 and u_1 be wind speeds measured at
           ! heights z_2 and z_1 respectively. THEN
           !
           !        u_2                    u_1
           ! ------------------  =  ------------------
           ! Ln((z_2-d_2)/z0_2)     Ln((z_1-d_1)/z0_1)
           !
           ! where
           !     u_i : wind speed at height z_i
           !     z_i : height at which the measurement was taken (m)
           !     d_i : zero plane displacement (m)
           !     z0_i: surface roughness LENgth or roughness height (m)
           !
           ! This equation assumes the atmosphere is neutrally stable,
           ! i.e., phi_m = 1, which implies psi_m = 0.
           !
           ! Given the wind speed at reference conditions (urh), compute
           ! the wind speed at the top of the canopy (uch). Assume the
           ! atmosphere is neutrally stable.
           !
           ! The wind speed at reference conditions (urh) is retrieved from the
           ! metereoLOGICAL file. The aerodynamic PARAMETERs for wind speed
           ! computations are set by the subroutine Get_Aerodynamic_PARAMETERs.
           ! In the absence of przm input file values, the routine assumes
           ! the conditions of the meteoroLOGICAL Daily Values File (*.dvf),
           ! i.e., Open Flat Terrain (used for MetereoLOGICAL Stations), and
           ! wind measurements normalized to 10 meters.
           ! See subroutines IniVar and Get_Aerodynamic_PARAMETERs.
           !   Wind_Reference_Height = 10.0
           !   Wind_z0 = 0.03
           !   Wind_D  = 0.0

           ! Computes zero displacement height, D (meter)
           ! and the roughness LENgth, Z0 (meter)
           Call Get_Crop_Params (zch, z0, d)

           ! urh: wind speed (meter/day) at reference height.
           !      units of WIND are cm/sec.
           !      1 cm/sec is equivaLENt to 864.0 meter/day
           !      subroutine canopy expects wind speed in meter/day
           ! uch: wind speed at the top of the canopy (zch)

           urh = WIND * 864.0
           uch = urh * Log((zch-d)/z0) /Log((uWind_Reference_Height-uWind_D)/uWind_z0)

           PWIND(1)= 0.0
           PWIND(2)= uch

           IF(Henry_unitless(K).GT.0.0.AND.URH.GT.0.0)THEN
               CALL Canopy(ATEMP,PWIND,ZCH,TOTCR )
               !this puts crcnc (which is array of 2 canopy resistances) into module for use in output ??? 
               ! CNDBDY: Boundary Layer's Conductance (cm/day)
               ! CONDUC: Canopy Conductance Including Boundary Layer's Conductance (cm/day)
               ! TOTCR:  Total canopy resistance (cm/day)
               CONDUC(K) = 1.0 / (1.0/CNDBDY(K) + TOTCR)
           ELSE
               TOTCR=0.0
           ENDIF
       ENDIF
    
    end subroutine volatilization_setup
    
    
    
    
    
    
    Subroutine Get_Crop_Params (Crop_Height, Z0, D)

      ! Given a height, compute the zero displacement
      ! height and the roughness LENgth
      !
      ! Crop_Height -- The canopy height (meter)
      ! D  -- The zero plane displacement height (meter)
      ! Z0 -- The roughness LENgth (meter).

      Implicit None

      REAL,  Intent(In) :: Crop_Height
      REAL, Intent(Out) :: D
      REAL, Intent(Out) :: Z0

      REAL :: logd, logz0

      If (Crop_Height > Minimum_Canopy_Height_m) THEN
         ! (Thibodeaux 1996) Regression for Zero plane displacement.
         ! height generated for 0.02 m < Crop_Height < 25 m
         ! Zero displacement height, m
         logd = 0.9793*Log10(Crop_Height) - 0.1536
         D    = 10.0 ** logd

         ! Roughness Height, m
         ! Valid for tall crops (Thibodeaux 1996).
         logz0 = 0.997*Log10(Crop_Height) - 0.883     ! copy 1/2; propagate changes
         Z0    = 10.0 ** logz0
      Else
         ! The crop height is less than the Minimum_Canopy_Height_m:
         ! the roughness height is set equal to the value
         ! of the regression for z0 evaluated at
         ! Crop_Height = Minimum_Canopy_Height_m
         logz0 = 0.997*Log10(Minimum_Canopy_Height_m) - 0.883  ! copy 2/2; propagate changes
         Z0    = 10.0 ** logz0
         D  = 0.0
      END IF

   End Subroutine Get_Crop_Params
    
    
    
    
    
     Subroutine Canopy (utemp, uwind, zch, totr)
      use constants_and_variables, ONLY: vonKarman
      Implicit None

      ! CH denotes canopy height
      ! RH denotes reference height
      !
      ! Utemp(2) - Air temperature (Celsius)
      !     Utemp(1): Air temperature at the soil surface
      !     Utemp(2): Ambient temperature
      !
      ! Uwind(2) - wind speed (meter/day)
      !     Uwind(1): wind speed at the soil surface
      !     Uwind(2): wind speed at the top of the canopy
      !
      ! ZCH  - canopy height (meter)  ( Zch > 0 )
      ! TOTR - total canopy resistance (cm/day)
      !
      ! CRC(2) - canopy resistance (cm/day)
      !     crc(1): total resistance in the lower half of the canopy
      !     crc(2): total resistance in the upper half of the canopy
      !
      ! History:
      ! * Wed Mar 24 15:29:51 EST 2004
      !   - translated to f90
      !   - total canopy resistance computed by integration,
      !     rather than a Riemann sum approximation.
      !   - The formulation of phi_h, phi_m, and psi_m updated as
      !     described in Thibodeaux. 1996. Environmental Chemodynamics:
      !     Movement of Chemicals in Air, Water, and Soil. Wiley.
      !     2nd Edition.
      !   - various bugs fixed; code cleaned.
      !
      ! * Modification date: 2/18/92 JAM
      !   To calculates the overall vertical transport resistance

      REAL,  Intent(in) :: zch
      REAL, Intent(out) :: totr
      REAL,  Intent(in) :: utemp(2)
      REAL,  Intent(in) :: uwind(2)
      REAL              :: crc(2)

      REAL :: gradw, gradt, RiNum, phi_h, phi_m, psi_m, ustar
      REAL :: d, z0, diffch, e2, c1, c2, ct, meanT, uch, zeta
      REAL :: temp0, RiMax

      REAL, PARAMETER :: Zero = 0.0

      ! If the Richardson number is "close" to Fuzzy, THEN
      ! the Richardson number is effectively equal to zero.
      REAL, PARAMETER :: Fuzzy = 0.003

      ! To convert from 1/m to 1/cm :
      ! value in 1/m is equivaLENt to (value * Im_to_Icm) 1/cm
      REAL, PARAMETER :: Im_to_Icm = 1.0e-02

      ! Temperature conversion: kelvin = Celsius + c2k
      REAL, PARAMETER :: c2k = 273.15

      ! Maximum value of the Richardson number
      REAL, PARAMETER :: Max_Richardson = 0.19

      ! g_grav: acceleration due to gravity. Express in m/day^2
      ! so that the Richardson number is DIMENSIONless.
      ! g_grav = 9.8 m/s^2
      !        = 9.8 m/s^2 * (86400 s/day)**2 = 7.31567E+10 m/day^2
      REAL, PARAMETER :: g_grav = 9.8 * 86400.0**2
      REAL, PARAMETER :: Pi = 3.14159265358979
      
      
      
      ! Gradients
      gradt = (utemp(2)-utemp(1)) / zch      ! kelvin/meter
      gradw = (uwind(2)-uwind(1)) / zch      ! 1/day
      meanT = Sum(utemp(1:2))/2 + c2k        ! mean Temperature, kelvin

      ! Computes Richardson number (RiNum) (DIMENSIONless).
      ! Louis J. Thibodeaux. 1996. Environmental Chemodynamics:
      ! Movement of Chemicals in Air, Water, and Soil. Wiley.
      ! 2nd Edition, p 373-375.
      ! Typically, -2.0 <= RiNum < 0.2, but for some of
      ! the PRZM scenarios RiNum was outside the nominal range.
      RiNum = g_grav / meanT * gradt / gradw**2

      ! Computes the DIMENSIONless height (zeta). See
      ! Louis J. Thibodeaux. 1996. Environmental Chemodynamics:
      ! Movement of Chemicals in Air, Water, and Soil. Wiley.
      ! 2nd Edition, p 379-381.
      If (RiNum < (-Fuzzy)) THEN
         ! Richardson number less than "fuzzy zero".
         zeta = RiNum

      ELSE IF (RiNum > Fuzzy) THEN
         ! Richardson number greater than "fuzzy zero".
         ! Make sure function for zeta is always evaluated with
         ! Richardson numbers < 0.2
         RiMax = Min (Max_Richardson, RiNum)
         zeta = RiMax / (1.0 - 5.0*RiMax)

      Else
         ! Abs(RiNum) <= Fuzzy, i.e., RiNum is equal to "Fuzzy Zero"
         zeta = Zero
      END IF

      ! phi_h: stability function for sensible heat (DIMENSIONless)
      ! phi_m: stability function for momentum (DIMENSIONless)
      ! psi_m: integrated momentum stability PARAMETER (DIMENSIONless)
      If (zeta < Zero) THEN
         phi_h = 1.0 / Sqrt(1.0 - 15.0*zeta)
         phi_m = Sqrt(phi_h)
         temp0 = (1.0+phi_m**2)/2.0 * ((1.0+phi_m)/2.0)**2
         psi_m = Pi/2.0 - 2.0*Atan(phi_m) + Log(temp0)
      Else
         ! zeta >= 0
         phi_h = 1.0 + 5.0*zeta
         phi_m = phi_h
         psi_m = -5.0 * zeta
      END IF

      ! Computes zero displacement height, D (meter)
      ! and the roughness LENgth, Z0 (meter)
      Call Get_Crop_Params (zch, z0, d)

      uch = Uwind(2)

      ! Compute friction velocity, USTAR (meter/day)
      ustar = vonKarman * uch / (Log((zch-d)/z0) - psi_m)

      ! Thermal eddy diffusivity (m^2/day) at canopy height
      diffch = ustar * vonKarman * (zch-d) / phi_h

      ! Compute canopy resistance.
      ! Let K(z) be the thermal eddy diffusivity at height z,
      !     K(z) = diffch * Exp[4*(z/ZCH - 1)]
      !
      ! The resistance at height z is 1/K(z).
      ! The resistance in the lower half of the canopy is:
      !
      !               z=zch/2   dz     e^2 (e^2-1) zch
      !     crc(1) = Integral [----] = ---------------
      !                z=0     K(z)       4 diffch
      !
      ! The resistance in the upper half of the canopy is:
      !
      !               z=zch     dz     (e^2-1) zch
      !     crc(2) = Integral [----] = -----------
      !              z=zch/2   K(z)     4 diffch
      !
      ! The total canopy resistance is:
      !     totr = crc(1) + crc(2), or
      !
      !             z=zch     dz     (e^4-1) zch
      !     totr = Integral [----] = -----------
      !              z=0     K(z)     4 diffch
      !
      ! Crc and totr are computed in day/meter.
      ! Im_to_Icm converts to day/cm.

      e2 = Exp(2.0)     ! e^2
      c2 = (e2 - 1) * zch / (4 * diffch) * Im_to_Icm
      c1 = e2 * c2
      ct = (e2**2 - 1) * zch / (4 * diffch) * Im_to_Icm

      crc(1) = c1
      crc(2) = c2
      totr   = ct


   End Subroutine Canopy
  
    
    
    
    
    
    
    
    end module volatilization
    
    
    !Module m_Wind
!   
!   Implicit None
!          ! The wind measurements in the meteoroLOGICAL Daily Values File (*.dvf)
!          ! were normalized to 10 meters, Open Flat Terrain.
!         
!         REAL,save      :: uWind_Reference_Height = 10.0  !(can be changed by input file)
!         REAL,PARAMETER :: uWind_z0 = 0.03
!         REAL,PARAMETER :: uWind_D  = 0.0
!               
!         
!         REAL, PARAMETER :: vonKarman = 0.4
!         REAL, Save      :: Height_stagnant_air_layer_cm = 0.5
!         REAL, PARAMETER :: Minimum_Canopy_Height_cm = 5.0
!         REAL, PARAMETER :: Minimum_Canopy_Height_m = Minimum_Canopy_Height_cm * 1.0E-02
!         
!         
!         
!!   Private
!!   Public :: Get_Crop_Params, inivar
!!   Public :: Get_Aerodynamic_PARAMETERs
!!
!!   ! Aerodynamic Environments, used by Get_Aerodynamic_PARAMETERs
!!   INTEGER, PARAMETER, Public :: t_Open_Flat_Terrain = 1
!!   INTEGER, PARAMETER, Public :: t_Class_A_Pan_Anemometer = 2
!!   INTEGER, PARAMETER, Public :: t_FAO_Reference_Short_Grass_Crop = 3
!!   INTEGER, PARAMETER, Public :: t_User = 4  ! Values READ from the input file.
!!
!!   ! The wind measurements in the meteoroLOGICAL Daily Values File (*.dvf)
!!   ! were normalized to 10 meters, Open Flat Terrain.
!!   !> Plan for the future: READ from the PRZM file the aerodynamic
!!   !>      environment of the metereoLOGICAL file (e.g., *.dvf).
!!
!!   INTEGER, Save, Public :: uWind_Environs = t_Open_Flat_Terrain
!!   REAL,    Save, Public :: uWind_Reference_Height, uWind_z0, uWind_D
!!
!!   ! Height of the stagnant air layer above the soil = 5 mm = 0.5 cm
!!   ! In the literature, generally denoted as "d".
!!   ! Value expressed in cm.
!!   REAL, Save, Public :: Height_stagnant_air_layer_cm = 0.5
!!
!!   ! Minimum Canopy Height = 0.05 m = 5 cm
!!   ! In the literature, Canopy Height is generally denoted as "Zch".
!!   ! Minimum_Canopy_Height_cm -- value in cm
!!   ! Minimum_Canopy_Height_c  -- value in meter
!!   ! 1 cm is equivaLENt to 1.0E-02 m
!!   REAL, PARAMETER, Public :: Minimum_Canopy_Height_cm = 5.0
!!   REAL, PARAMETER, Public :: Minimum_Canopy_Height_m = Minimum_Canopy_Height_cm * 1.0E-02
!!
!!   ! von Karman's constant, DIMENSIONless
!!   REAL, PARAMETER, Public :: vonKarman = 0.4
!!
!!   
!!Contains
!!
!! Subroutine IniVar ()
!!      !calculate wind PARAMETERs
!!      Implicit None
!!
!!      uWind_Environs = t_Open_Flat_Terrain
!!      Call Get_Aerodynamic_PARAMETERs(uWind_Environs, uWind_Reference_Height, uWind_z0, uWind_D)
!!   
!!
!!      ! Description taken from the przm manual.
!!      ! przm input file RECORD 31
!!      ! ZWIND: height of wind speed measurement above the soil surface (meter).
!!      ! The wind speed anemometer is usually fixed at 10 meters (30 feet)
!!      ! above the ground surface.  This height may differ at some weather
!!      ! stations such as at a class A station where the anemometer may be
!!      ! attached to the evaporation pan.  The correct value can be obtained
!!      ! from the meteoroLOGICAL data reports for the station whose data are
!!      ! in the simulation.
!!      !
!!      ! Thu Apr 08 16:54:02 EDT 2004 --
!!      ! o The variable ZWIND was replaced with uWind_Reference_Height
!!
!!      ! The wind measurements in the meteoroLOGICAL Daily Values File (*.dvf)
!!      ! were normalized to 10 meters, Open Flat Terrain.    
!!
!!   End Subroutine IniVar
!!
!!
!!
!!
!!   Subroutine Get_Aerodynamic_PARAMETERs (Wind_Env, Wind_Reference_Height, Wind_z0, Wind_D)
!!
!!      ! Returns Aerodynamic PARAMETERs for wind speed computations
!!      !         for the named environment.
!!      !
!!      ! Wind_Env -- (INTEGER) environment type: t_Open_Flat_Terrain or
!!      !             t_Class_A_Pan_Anemometer, t_FAO_Reference_Short_Grass_Crop
!!      !
!!      ! Wind_Reference_Height -- reference height (meter)
!!      ! Wind_z0 -- surface roughness LENgth or roughness height (meter)
!!      ! Wind_D  -- zero plane displacement (meter)
!!      !
!!      ! Let u_2 and u_1 be wind speeds measured at
!!      ! heights z_2 and z_1 respectively. THEN
!!      !
!!      !        u_2                    u_1
!!      ! ------------------  =  ------------------
!!      ! Ln((z_2-d_2)/z0_2)     Ln((z_1-d_1)/z0_1)
!!      !
!!      ! where
!!      !     u_i : wind speed at height z_i
!!      !     z_i : height at which the measurement was taken (m)
!!      !     d_i : zero plane displacement (m)
!!      !     z0_i: surface roughness LENgth or roughness height (m)
!!      !
!!      ! This equation assumes the atmosphere is neutrally stable,
!!      ! i.e., phi_m = 1, which implies psi_m = 0. See the PRZM
!!      ! manual and subroutine "Canopy" for more information.
!!
!!      Implicit None
!!
!!      INTEGER, Intent(In) :: Wind_Env
!!      REAL,   Intent(Out) :: Wind_Reference_Height
!!      REAL,   Intent(Out) :: Wind_z0
!!      REAL,   Intent(Out) :: Wind_D
!!
!!      Select Case (Wind_Env)
!!      Case(t_Open_Flat_Terrain)
!!         ! The wind measurements in the meteoroLOGICAL Daily Values File
!!         ! (*.dvf) were normalized to 10 meters, Open Flat Terrain (used
!!         ! for MetereoLOGICAL Stations):
!!         Wind_Reference_Height = 10.0
!!         Wind_z0 = 0.03
!!         Wind_D  = 0.0
!!
!!      Case(t_Class_A_Pan_Anemometer)
!!         Wind_Reference_Height = 0.6
!!         Wind_z0 = 0.01476
!!         Wind_D  = 0.08
!!
!!      Case(t_FAO_Reference_Short_Grass_Crop)
!!         Wind_Reference_Height = 2.0
!!         Wind_z0 = 0.01476
!!         Wind_D  = 0.08
!!
!!      Case(t_User)
!!         !> This is a plan for the future.
!!         !> User provides all PARAMETERs in the PRZM input file.
!!         !> Requires modification of the PRZM input file.
!!         Wind_Reference_Height = uWind_Reference_Height
!!         Wind_z0 = uWind_z0
!!         Wind_D  = uWind_D
!!
!!      Case Default
!!         ! Error
!!
!!      End Select
!!   End Subroutine Get_Aerodynamic_PARAMETERs
!
!Contains
!
!   Subroutine Get_Crop_Params (Crop_Height, Z0, D)
!
!       Given a height, compute the zero displacement
!       height and the roughness LENgth
!      
!       Crop_Height -- The canopy height (meter)
!       D  -- The zero plane displacement height (meter)
!       Z0 -- The roughness LENgth (meter).
!
!      Implicit None
!
!      REAL,  Intent(In) :: Crop_Height
!      REAL, Intent(Out) :: D
!      REAL, Intent(Out) :: Z0
!
!      REAL :: logd, logz0
!
!      If (Crop_Height > Minimum_Canopy_Height_m) THEN
!          (Thibodeaux 1996) Regression for Zero plane displacement.
!          height generated for 0.02 m < Crop_Height < 25 m
!          Zero displacement height, m
!         logd = 0.9793*Log10(Crop_Height) - 0.1536
!         D    = 10.0 ** logd
!
!         ! Roughness Height, m
!         ! Valid for tall crops (Thibodeaux 1996).
!         logz0 = 0.997*Log10(Crop_Height) - 0.883     ! copy 1/2; propagate changes
!         Z0    = 10.0 ** logz0
!      Else
!         ! The crop height is less than the Minimum_Canopy_Height_m:
!         ! the roughness height is set equal to the value
!         ! of the regression for z0 evaluated at
!         ! Crop_Height = Minimum_Canopy_Height_m
!         logz0 = 0.997*Log10(Minimum_Canopy_Height_m) - 0.883  ! copy 2/2; propagate changes
!         Z0    = 10.0 ** logz0
!         D  = 0.0
!      END IF
!
!   End Subroutine Get_Crop_Params
!
!
!End Module m_Wind