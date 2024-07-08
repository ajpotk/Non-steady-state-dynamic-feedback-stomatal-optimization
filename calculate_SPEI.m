function [SPEI] = calculate_SPEI(Date, P, T_a, R_0)

Year_begin = year(Date(1));
Year_end = year(Date(end));

N_Month = 12*(1 + Year_end - Year_begin);

PET_month = nan(1, N_Month); %[mm]
P_Month = nan(1, N_Month); %[mm]
CWD = nan(size(Date)); %climatic water deficit [mm]
for Year_local = Year_begin:Year_end
    
    %Thornthwaite equation (1948)
    is_Year = (year(Date) == Year_local);
    Date_local = Date(is_Year == 1);
    T_a_local = T_a(is_Year == 1);
    T_a_local(T_a_local < 0) = 0;
    R_0_local = R_0(is_Year == 1);
    P_local = 1e3 * P(is_Year == 1); %[mm]
    
    T_a_Monthly_mean = nan(1, 12);
    P_Monthly_sum = nan(1, 12); %[mm]
    Average_Day_Length = nan(1, 12);
    Days_in_Month = nan(1, 12);
    for j = 1:12
        is_Month = (month(Date_local) == j);
        T_a_Monthly_mean(j) = mean(T_a_local(is_Month == 1));
        P_Monthly_sum(j) = sum(P_local(is_Month == 1));
        Average_Day_Length(j) = 24 * sum(R_0_local(is_Month == 1) > 0)/sum(is_Month == 1); %[hr]
        Days_in_Month(j) = eomday(Year_local, j); %[day]
    end
    Heat_index = sum((T_a_Monthly_mean/5).^1.514);
    alpha = 6.75e-7 * Heat_index^3 - 7.71e-5 * Heat_index^2 + 1.792e-2 * Heat_index + 0.49239;
    PET_month(12*(Year_local - Year_begin) + (1:12)) = 16 * Average_Day_Length/12 .* Days_in_Month/30 .* (10*T_a_Monthly_mean/Heat_index).^alpha; %[mm]
    P_Month(12*(Year_local - Year_begin) + (1:12)) = P_Monthly_sum; %[mm]
    
    % store monthly climatic water deficit at the resolution of original environmental conditions
    for j = 1:12
        is_Year_and_Month = ((year(Date) == Year_local) + (month(Date) == j) == 2);
        CWD(is_Year_and_Month == 1) = P_Month(12*(Year_local - Year_begin) + j) - PET_month(12*(Year_local - Year_begin) + j);
    end
    
end

% fit log-logistic probability density function to monhtly climatic water deficit through the probability-weighted moments
CWD_Month = sort(P_Month - PET_month); %climatic water deficit [mm]
F = ((1:N_Month) - 0.35)/N_Month;
PWM_0 = mean(CWD_Month);
PWM_1 = sum((1 - F).*CWD_Month)/N_Month;
PWM_2 = sum((1 - F).^2.*CWD_Month)/N_Month;
b = (2*PWM_1 - PWM_0)/(6*PWM_1 - PWM_0 - 6*PWM_2); 
a = b*(PWM_0 - 2*PWM_1)/gamma(1+1/b)/gamma(1-1/b);
c = PWM_0 - a*gamma(1+1/b)*gamma(1-1/b);

% CDF for log-logistic
CDF_func = @(x) (1 + (a./(x-c)).^b).^(-1);

% SPEI for monhtly CWD
Prob_Month = 1 - CDF_func(CWD_Month);
W_Month = (-2*log(Prob_Month)).^0.5;
SPEI_Month = W_Month - (2.515517 + 0.802853*W_Month + 0.010328*W_Month.^2)./(1 + 1.432788*W_Month + 0.189269*W_Month.^2 + 0.001308*W_Month.^3);

% check that mean is 0 and std is 1
error_mean = abs(0 - mean(SPEI_Month));
error_std = abs(1 - std(SPEI_Month));
if (error_mean > 0.05) || (error_std > 0.05)
    error('ERROR: The mean and std of monthly SPEI are not 0 and 1, respectively!')
end

% SPEI for CWD at the resolution of original environmental conditions
Prob = 1 - CDF_func(CWD);
W = (-2*log(Prob)).^0.5;
SPEI = W - (2.515517 + 0.802853*W + 0.010328*W.^2)./(1 + 1.432788*W + 0.189269*W.^2 + 0.001308*W.^3);

end

