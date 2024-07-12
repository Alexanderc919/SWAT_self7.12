      subroutine snom
      
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
!!    sub_smfmn    |mm/deg C/day  |Minimum melt rate for snow during year (Dec.
!!                                |21) where deg C refers to the air temperature
!!    sub_smfmx    |mm/deg C/day  |Maximum melt rate for snow during year (June
!!                                |21) where deg C refers to the air temperature
!!                                |SMFMX and SMFMN allow the rate of snow melt
!!                                |to vary through the year. These parameters
!!                                |are accounting for the impact of soil
!!                                |temperature on snow melt.
!!    sub_smtmp    |deg C         |Snow melt base temperature
!!                                |Mean air temperature at which snow melt will 
!!                                |occur.
!!    sno_hru(:)   |mm H2O        |amount of water in snow in HRU on current day
!!    snocov1      |none          |1st shape parameter for snow cover equation
!!                                |This parameter is determined by solving the
!!                                |equation for 50% snow cover
!!    snocov2      |none          |2nd shape parameter for snow cover equation
!!                                |This parameter is determined by solving the
!!                                |equation for 95% snow cover
!!    snocovmx     |mm H2O        |Minimum snow water content that corresponds
!!                                |to 100% snow cover. If the snow water content
!!                                |is less than SNOCOVMX, then a certain 
!!                                |percentage of the ground will be bare.
!!    snoeb(:,:)   |mm H2O        |snow water content in elevation band on 
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
!!    snoeb(:,:)   |mm H2O        |snow water content in elevation band on 
!!                                |current day
!!    snofall      |mm H2O        |amount of precipitation falling as freezing 
!!                                |rain/snow on day in HRU
!!    snomlt       |mm H2O        |amount of water in snow melt for the day in 
!!                                |HRU
!!    snotmp(:)    |deg C         |temperature of snow pack in HRU
!!    snotmpeb(:,:)|deg C         |temperature of snow pack in elevation band
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ib          |none          |counter
!!    j           |none          |HRU number
!!    smfac       |
!!    smleb       |mm H2O        |amount of snow melt in elevation band on 
!!                               |current day
!!    smp         |mm H2O        |precipitation on current day for HRU
!!    snocov      |none          |fraction of HRU area covered with snow
!!    sum         |mm H2O        |snow water content in HRU on current day
!!    xx          |none          |ratio of amount of current day's snow water
!!                               |content to the minimum amount needed to
!!                               |cover ground completely
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Real, Sin, Exp

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~


      use parm

      integer :: j, ib
      real :: sum, smp, smfac, smleb, iceleb !!fix_here
      real :: xx, snocov 

      j = 0
      j = ihru
      sum =0.
      smp =0.
      isub = hru_sub(j)

      
     
      if (elevb(1,isub) > 0. .and. elevb_fr(1,isub) > 0.) then 

!! elevation bands
        !! compute snow fall and melt for each elevation band
        do ib = 1, 10

          if (elevb_fr(ib,isub) < 0.) exit
          snotmpeb(ib,j) = snotmpeb(ib,j) * (1.-sub_timp(ib,isub)) +
     &                                tavband(ib,j) * sub_timp(ib,isub)

          !!if (tavband(ib,j) < sub_sftmp(ib,isub)) then


          if (tmxband(ib,j) < sub_sftmp(ib,isub).and.tmaband(ib,j) 
     &     < sub_sftmp_accu(ib,isub)) then

            !! compute snow fall if temperature is below sftmp
            snoeb(ib,j) = snoeb(ib,j) + pcpband(ib,j)
            snofall = snofall + pcpband(ib,j) * elevb_fr(ib,isub)

          else
          
          
          
          !!增加是否有降雨存在，并分别计算两种情况下的融雪量
           !if(precipday==0.)
          if (pcpband(ib,j)==0.)then
          
            if (tmxband(ib,j) >= sub_smtmp(ib,isub).and.tmaband(ib,j)>=
     &      sub_smtmp_accu(ib,isub)) then
              smfac = 0.
              smleb = 0.
              smfac = (sub_smfmx(ib,isub) + sub_smfmn(ib,isub)) / 2. +  
     &          Sin((iida - 81) / 58.09) *                              
     &          (sub_smfmx(ib,isub) - sub_smfmn(ib,isub)) / 2.    !! 365/2pi = 58.09
              smleb = smfac * (((snotmpeb(ib,j) + tmxband(ib,j)) / 2.)  
     &                                             - sub_smtmp(ib,isub))

              !! adjust for areal extent of snow cover
              if (snoeb(ib,j) < snocovmx) then
                xx = 0.
                snocov = 0.
                xx = snoeb(ib,j) / snocovmx
                snocov = xx / (xx + Exp(snocov1 - snocov2 * xx))
              else
                snocov = 1.
              endif
              smleb = smleb * snocov
              if (smleb < 0.) smleb = 0.
              if (smleb > snoeb(ib,j)) smleb = snoeb(ib,j)
              snoeb(ib,j) = snoeb(ib,j) - smleb
              snomlt = snomlt + smleb * elevb_fr(ib,isub)
             endif
          else if(pcpband(ib,j)>0.)then
              
              
              if (tmxband(ib,j) >= sub_smtmp(ib,isub).and.tmaband(ib,j)>=
     &      sub_smtmp_accu(ib,isub)) then
              smfac = 0.
              smleb = 0.
              smfac = (sub_smfmx(ib,isub) + sub_smfmn(ib,isub)) / 2. +  
     &          Sin((iida - 81) / 58.09) *                              
     &          (sub_smfmx(ib,isub) - sub_smfmn(ib,isub)) / 2.    !! 365/2pi = 58.09
              !!计算融雪时加进降雨带来的能量所产生的的融雪。
              smleb = smfac * (((snotmpeb(ib,j) + tmxband(ib,j)) / 2.)  
     &                - sub_smtmp(ib,isub))+  
     &    sub_crain(ib,isub)*pcpband(ib,j)*(tmpav(j)-sub_smtmp(ib,isub))

              !! adjust for areal extent of snow cover
              if (snoeb(ib,j) < snocovmx) then
                xx = 0.
                snocov = 0.
                xx = snoeb(ib,j) / snocovmx
                snocov = xx / (xx + Exp(snocov1 - snocov2 * xx))
              else
                snocov = 1.
              endif
              call ICEm2(iceleb,ib)  !!fix_here
              smleb = smleb * snocov!+iceleb
              if (smleb < 0.) smleb = 0.
              if (smleb > snoeb(ib,j)) smleb = snoeb(ib,j)
              snoeb(ib,j) = snoeb(ib,j) - smleb
              snomlt = snomlt + smleb * elevb_fr(ib,isub)
              
              
              
              
              snom_hru(2)=snom_hru(2)+(smleb-iceleb) * elevb_fr(ib,isub)!融雪
              snom_hru(3)=snom_hru(3)+iceleb * elevb_fr(ib,isub)!融冰
          if (snom_hru(2)<0)snom_hru(2)=0.
          if (snom_hru(3)<0)snom_hru(3)=0.
          
          
          
            endif
           end if
          endif
          !!改进部分结束
            !! compute snow melt if temperature is above smtmp
     !!       if (tmxband(ib,j) > sub_smtmp(ib,isub)) then
     !!         smfac = 0.
     !!         smleb = 0.
     !!         smfac = (sub_smfmx(ib,isub) + sub_smfmn(ib,isub)) / 2. +  
     !!&          Sin((iida - 81) / 58.09) *                              
     !!&          (sub_smfmx(ib,isub) - sub_smfmn(ib,isub)) / 2.    !! 365/2pi = 58.09
     !!         smleb = smfac * (((snotmpeb(ib,j) + tmxband(ib,j)) / 2.)  
     !!&                                             - sub_smtmp(ib,isub))
     !!
     !!         !! adjust for areal extent of snow cover
     !!         if (snoeb(ib,j) < snocovmx) then
     !!           xx = 0.
     !!           snocov = 0.
     !!           xx = snoeb(ib,j) / snocovmx
     !!           snocov = xx / (xx + Exp(snocov1 - snocov2 * xx))
     !!         else
     !!           snocov = 1.
     !!         endif
     !!         smleb = smleb * snocov
     !!         if (smleb < 0.) smleb = 0.
     !!         if (smleb > snoeb(ib,j)) smleb = snoeb(ib,j)
     !!         snoeb(ib,j) = snoeb(ib,j) - smleb
     !!         snomlt = snomlt + smleb * elevb_fr(ib,isub)
     !!       endif
     !!     endif
          sum = sum + snoeb(ib,j) * elevb_fr(ib,isub)
          smp = smp + pcpband(ib,j) * elevb_fr(ib,isub)
          
        end do
        


        !! add/sub aggregate snow fall and melt from effective precip 
        !! and snow cover
        precipday = smp + snomlt - snofall

        if (precipday < 0.) precipday = 0.
        if (nstep > 0) then
          do ii = 1, nstep
            precipdt(ii+1) = precipdt(ii+1) + (snomlt - snofall) / nstep
            if (precipdt(ii+1) < 0.) precipdt(ii+1) = 0.
          end do
        end if
        sno_hru(j) = sum

      else
          
!! no elevation bands
      call ICEm2(iceleb,ib)  !!fix_here
	ib = 1
          
        !! estimate snow pack temperature


      !! 更改的部分开始
        snotmp(j)=snotmp(j) * (1. - sub_timp(ib,isub)) + tmpav(j) *     
     &            sub_timp(ib,isub)
        if (tmx(j) < sub_sftmp(ib,isub).and.tmaband(ib,isub) 
     &     < sub_sftmp_accu(ib,isub)) then
        !!if (tmpav(j) <= sub_sftmp(ib,isub)) then
          !! calculate snow fall
          sno_hru(j) = sno_hru(j) + precipday
          snofall = precipday
          precipday = 0.
          precipdt = 0.
          endif
 
       
                if (pcpband(ib,isub)==0.)then
               
        if (tmx(j) > sub_smtmp(ib,isub) .and. tmaband(ib,isub)>=
     &      sub_smtmp_accu(ib,isub).and.sno_hru(j) > 0.) then
          !! adjust melt factor for time of year
          smfac = 0.
          snomlt = 0.
          smfac = (sub_smfmx(ib,isub) + sub_smfmn(ib,isub)) / 2. +      
     &       Sin((iida - 81) / 58.09) *                                 
     &       (sub_smfmx(ib,isub) - sub_smfmn(ib,isub)) / 2.    !! 365/2pi = 58.09
          snomlt = smfac * (((snotmp(j)+tmx(j))/2.)-sub_smtmp(ib,isub))

          !! adjust for areal extent of snow cover
          if (sno_hru(j) < snocovmx) then
            xx = 0.
            xx = sno_hru(j) / snocovmx
            snocov = xx / (xx + Exp(snocov1 - snocov2 * xx))
          else
            snocov = 1.
          endif
          snomlt = snomlt * snocov+iceleb
          if (snomlt < 0.) snomlt = 0.
          if (snomlt > sno_hru(j)) snomlt = sno_hru(j)
          sno_hru(j) = sno_hru(j) - snomlt
          precipday = precipday + snomlt
          if (nstep > 0) then
            do ii = 1, nstep
             precipdt(ii+1) = precipdt(ii+1) + snomlt / nstep
            end do
          end if
          if (precipday < 0.) precipday = 0.
        else
          snomlt = 0.
        end if
      else if (pcpband(ib,isub)>0.)then
      if (tmx(j) > sub_smtmp(ib,isub) .and. tmaband(ib,isub)>=
     &      sub_smtmp_accu(ib,isub).and.sno_hru(j) > 0.) then
          !! adjust melt factor for time of year
          smfac = 0.
          snomlt = 0.
          smfac = (sub_smfmx(ib,isub) + sub_smfmn(ib,isub)) / 2. +      
     &       Sin((iida - 81) / 58.09) *                                 
     &       (sub_smfmx(ib,isub) - sub_smfmn(ib,isub)) / 2.    !! 365/2pi = 58.09
          !!计算融雪时加进降雨带来的能量所产生的的融雪。
          snomlt = smfac * (((snotmp(j)+tmx(j))/2.)-sub_smtmp(ib,isub))+  
     & sub_crain(ib,isub)* precipday*(tmpav(j)-sub_smtmp(ib,isub))

          !! adjust for areal extent of snow cover
          if (sno_hru(j) < snocovmx) then
            xx = 0.
            xx = sno_hru(j) / snocovmx
            snocov = xx / (xx + Exp(snocov1 - snocov2 * xx))
          else
            snocov = 1.
          endif
          snomlt = snomlt * snocov!+iceleb
          
          
          
          snom_hru(2)=snom_hru(2)+(smleb-iceleb) * elevb_fr(ib,isub)!融雪
          snom_hru(3)=snom_hru(3)+iceleb * elevb_fr(ib,isub)!融冰
          
          
          
          
          
          if (snom_hru(2)<0)snom_hru(2)=0.
          if (snom_hru(3)<0)snom_hru(3)=0.
          if (snomlt < 0.) snomlt = 0.
          if (snomlt > sno_hru(j)) snomlt = sno_hru(j)
          sno_hru(j) = sno_hru(j) - snomlt
          precipday = precipday + snomlt
          if (nstep > 0) then
            do ii = 1, nstep
             precipdt(ii+1) = precipdt(ii+1) + snomlt / nstep
            end do
          end if
          if (precipday < 0.) precipday = 0.
        else
          snomlt = 0.
        end if
      endif
      !! 更改的部分结束。2018-1-19    
          
     !!   if (tmx(j) > sub_smtmp(ib,isub) .and. sno_hru(j) > 0.) then
     !!     !! adjust melt factor for time of year
     !!     smfac = 0.
     !!     snomlt = 0.
     !!     smfac = (sub_smfmx(ib,isub) + sub_smfmn(ib,isub)) / 2. +      
     !!&       Sin((iida - 81) / 58.09) *                                 
     !!&       (sub_smfmx(ib,isub) - sub_smfmn(ib,isub)) / 2.    !! 365/2pi = 58.09
     !!     snomlt = smfac * (((snotmp(j)+tmx(j))/2.)-sub_smtmp(ib,isub))
     !!
     !!     !! adjust for areal extent of snow cover
     !!     if (sno_hru(j) < snocovmx) then
     !!       xx = 0.
     !!       xx = sno_hru(j) / snocovmx
     !!       snocov = xx / (xx + Exp(snocov1 - snocov2 * xx))
     !!     else
     !!       snocov = 1.
     !!     endif
     !!     snomlt = snomlt * snocov
     !!     if (snomlt < 0.) snomlt = 0.
     !!     if (snomlt > sno_hru(j)) snomlt = sno_hru(j)
     !!     sno_hru(j) = sno_hru(j) - snomlt
     !!     precipday = precipday + snomlt
     !!     if (nstep > 0) then
     !!       do ii = 1, nstep
     !!        precipdt(ii+1) = precipdt(ii+1) + snomlt / nstep
     !!       end do
     !!     end if
     !!     if (precipday < 0.) precipday = 0.
     !!   else
     !!     snomlt = 0.
     !!   end if
      end if
       
      snom_hru(1)=snofall!下雪
      !IF (snofall>1)PRINT*,snofall

      return
      end