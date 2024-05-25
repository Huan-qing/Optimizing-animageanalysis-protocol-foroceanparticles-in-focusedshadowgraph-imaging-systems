function [hour_ctd,minute_ctd,second_ctd]=julian_conversion(julian_day)
    hour_ctd=floor(24*(julian_day-floor(julian_day)));
    hour_ctd_res=24*(julian_day-floor(julian_day));
    minute_ctd=floor(60*(hour_ctd_res-floor(hour_ctd_res)));
    minute_ctd_res=60*(hour_ctd_res-floor(hour_ctd_res));    
    second_ctd=round(60*(minute_ctd_res-minute_ctd));
    
    % if the second=60, it should be equal to 1 min
    if second_ctd==60
        second_ctd=0;
        minute_ctd=minute_ctd+1;
    end
    % the same for minute
    if minute_ctd==60
        minute_ctd=0;
        second_ctd=round(60*(minute_ctd_res-minute_ctd));
        hour_ctd=hour_ctd+1;
    end
    if hour_ctd==24
      hour_ctd=0;
    end
    time(:,1)=hour_ctd;
    time(:,2)=minute_ctd;
    time(:,3)=second_ctd;

