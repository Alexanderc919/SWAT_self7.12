 subroutine ICEm2(iceleb,ib)    !!fix_here
    
    !!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine predicts daily snom melt when the average air
!!    temperature exceeds 0 degrees Celcius

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name         |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    elevb(:,:)   |m             |elevation at center of band
!!    elevb_fr(:,:)|none          |fraction of subbasin area within elevation 
!!                                |band
!!    iida         |julian date   |day being simulated (current julian date)
!!    ihru         |none          |HRU number
!!    pcpband(:,:) |mm H2O        |precipitation for the day in band in HRU
!!    precipday    |mm H2O        |precipitation on the current day in the HRU
!!    sub_sftmp    |deg C         |Snowfall temperature
!!                                |Mean air temperature at which precipitation
!!                                |is equally likely to be rain as snow/freezing
!!                                |rain.
!!    sub_icefmn    |mm/deg C/day  |Minimum melt rate for ice during year (Dec.
!!                                |21) where deg C refers to the air temperature
!!    sub_icefmx    |mm/deg C/day  |Maximum melt rate for ice during year (June
!!                                |21) where deg C refers to the air temperature
!!                                |iceFMX and iceFMN allow the rate of ice melt
!!                                |to vary through the year. These parameters
!!                                |are accounting for the impact of soil
!!                                |temperature on snow melt.
!!    sub_CT    |deg C         |ICE melt base temperature
!!                                |Mean air temperature at which snow melt will 
!!                                |occur.
!!    sno_hru(:)   |mm H2O        |amount of water in snow in HRU on current day
!!    icecov1      |none          |1st shape parameter for ice cover equation
!!                                |This parameter is determined by solving the
!!                                |equation for 50% ice cover
!!    icecov2      |none          |2nd shape parameter for ice cover equation
!!                                |This parameter is determined by solving the
!!                                |equation for 95% ice cover
!!    icecovmx     |mm H2O        |Minimum ice water content that corresponds
!!                                |to 100% ice cover. If the ice water content
!!                                |is less than ICECOVMX, then a certain 
!!                                |percentage of the ground will be bare.
!!    iceeb(:,:)   |mm H2O        |ice water content in elevation band on 
!!                                |current day
!!    snotmp(:)    |deg C         |temperature of snow pack in HRU
!!    snotmpeb(:,:)|deg C         |temperature of snow pack in elevation band
!!    tavband(:,:) |deg C         |average temperature for the day in band in HRU
!!    sub_timp     |none          |Snow pack temperature lag factor (0-1)
!!                                |1 = no lag (snow pack temp=current day air
!!                                |temp) as the lag factor goes to zero, the
!!                                |snow pack's temperature will be less
!!                                |influenced by the current day's air 
!!                                |temperature
!!    tmpav(:)     |deg C         |average air temperature on current day for 
!!                                |HRU
!!    tmx(:)       |deg C         |maximum air temperature on current day for 
!!                                |HRU
!!    tmxband(:,:) |deg C         |maximum temperature for the day in band in HRU
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name         |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    precipday    |mm H2O        |amount of water in effective precipitation
!!                                |in HRU
!!    precipdt(:)  |mm H2O        |precipitation for the time step during day
!!    sno_hru(:)   |mm H2O        |amount of water in snow in HRU on current day
!!    iceeb(:,:)   |mm H2O        |ice water content in elevation band on 
!!                                |current day
!!    snofall      |mm H2O        |amount of precipitation falling as freezing 
!!                                |rain/snow on day in HRU
!!    icemlt       |mm H2O        |amount of water in ice melt for the day in 
!!                                |HRU
!!    snotmp(:)    |deg C         |temperature of snow pack in HRU
!!    snotmpeb(:,:)|deg C         |temperature of snow pack in elevation band
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ib          |none          |counter
!!    j           |none          |HRU number
!!    icefac       |
!!    iceleb       |mm H2O        |amount of ice melt in elevation band on 
!!                               |current day
!!    smp         |mm H2O        |precipitation on current day for HRU
!!    icecov      |none          |fraction of HRU area covered with snow
!!    sum         |mm H2O        |snow water content in HRU on current day
!!    yy          |none          |ratio of amount of current day's snow water
!!                               |content to the minimum amount needed to
!!                               |cover ground completely
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Real, Sin, Exp

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~
    
    use parm

      integer :: j, ib
      real :: sum, smp, icefac, iceleb
      real :: yy, icecov 

         j = 0
         j = ihru
         sum =0.
         smp =0.
         icemlt=0.
         isub = hru_sub(j)
          !!增加一个判断，当高程大于雪与冰川的分界线的程（Elesnoline）时才进行融冰计算，否则为0.
          if(elevb(1,isub) >= Elesnoline)then
           if (tmxband(ib,j) > sub_CT(ib,isub)) then
              icefac = 0.
              iceleb = 0.
    icefac = (sub_icefmx(ib,isub) + sub_icefmn(ib,isub)) / 2. + Sin((iida - 81) / 58.09) * (sub_icefmx(ib,isub) - sub_icefmn(ib,isub)) / 2. 
    iceleb = icefac * (((snotmpeb(ib,j) + tmxband(ib,j)) / 2.)  - sub_CT(ib,isub))     
    !&         + Sin((iida - 81) / 58.09) *                              
    !&          (sub_icefmx(ib,isub) - sub_icefmn(ib,isub)) / 2.    !! 365/2pi = 58.09
    !          iceleb = icefac * (((snotmpeb(ib,j) + tmxband(ib,j)) / 2.)  
    !&                                             - sub_CT(ib,isub))

              !! adjust for areal extent of snow cover
              if (iceeb(ib,j) < icecovmx) then
                yy = 0.
                icecov = 0.
                yy = iceeb(ib,j) / icecovmx
                icecov = yy / (yy + Exp(icecov1 - icecov2 * yy))
              else
                icecov = 1.
              endif
              iceleb = iceleb * icecov
              
            endif
          else
            iceleb = 0.
    endif
    
    
    end
    
    