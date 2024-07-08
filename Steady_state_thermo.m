function [T_L_vect, VPD_L_vect, g_w_vect] = Steady_state_thermo( T_a_vect, E_vect, RH_vect, P_atm_vect, R_abs, g_H_a, emiss_L )

%% Inputs:
%   T_a_vect -- vector of air temperatures [C]
%   E_vect -- vector of transpiration rates [mol m-2 s-1]
%   RH_vect -- vector of relative humidities [-]
%   P_atm_vect -- vector of atmospheric pressure [kPa]
%   R_abs -- absorbed shortwave radiation [W m-2]
%   g_H_a -- boundary layer heat conductance [mol m-2 s-1]
%   emiss_L -- leaf emissivity [-]

%% Constants
c_p = 29.3; %specific heat capacity of iar [J mol-1 K-1]
sigma = 5.67e-8; %Stefan–Boltzman constant [W m-2 K-4]
m_w = 18e-3; %molar mass of water [kg mol-1]
lambda_E_0C = 2.5e6; %latent heat of vaporization at reference temperature of 0C [J kg-1]
dlambda_E_dT_L = -2365; %slope between latent heat of vaporization, lambda_E, and leaf temperature in degrees Celcius, T_L [J kg-1 C-1]

%% Steady-state radiation balance
lambda_E_func = @(T_L) lambda_E_0C + dlambda_E_dT_L * T_L; %latent heat of vaporization at a given leaf temperature [J kg-1]
Rad_balance_func = @(T_L, T_a, E) R_abs - 2*sigma*(emiss_L*(T_L+273.15).^4 - (T_a+273.15)^4) - 2*c_p*g_H_a*(T_L - T_a) - m_w*lambda_E_func(T_L)*E; %in [W m-2]
N_vect = length(E_vect);
T_L_vect = nan(1, N_vect);
for i = 1:N_vect
    Rad_balance_local_func = @(T_L) Rad_balance_func(T_L, T_a_vect(i), E_vect(i));
    T_L_vect(i) = fzero(Rad_balance_local_func, T_a_vect(i));
end

%% Determine VPD_L
e_L_vect = 0.61078 * exp(17.27 * T_L_vect ./ (T_L_vect + 237.3)); % saturated leaf vapor pressure in [kPa] -- Teten's equation
e_a = RH_vect .* 0.61078 .* exp(17.27 * T_a_vect ./ (T_a_vect + 237.3)); %air vapor pressure in [kPa] -- Teten's equation
VPD_L_vect = e_L_vect - e_a; %vapor pressure deficit in [kPa]
g_v_vect = E_vect .* P_atm_vect ./ VPD_L_vect; %total vapor conductance [mol m-2 s-1]
g_b = g_H_a; %boundary layer vapor conductance [mol m-2 s-1]
g_w_vect = 1./(1./g_v_vect - 1/g_b); %stomatal conductance to vapor [mol m-2 s-1]

end

