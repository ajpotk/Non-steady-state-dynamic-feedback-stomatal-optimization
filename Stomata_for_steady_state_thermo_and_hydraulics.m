function [outputs_stomata_steady_state] = ...
Stomata_for_steady_state_thermo_and_hydraulics(SWC_L, alpha, a_f_max, pi_L_star, beta, epsilon_L_max, pi_L_0_25C, ...
                                               k_L_max_25C, A, Q_10_k_L, psi_soil, ...
                                               T_a, RH, P_atm, c_a, ...
                                               R_abs, g_H_a, emiss_L, ...
                                               V_cmax25_func, ...
                                               varargin)

%% Inputs
% SWC_L -- leaf saturated water content [kg kg-1]
% alpha -- number of moles of carbon in mole of first product of photosynthesis [mol mol-1]
% a_f_max -- maximum apoplastic fraction [-]
% pi_L_star -- shape parameter for relationship between apoplastic fraction and osmotic potential [MPa]
% beta -- shape parameter for relationship between apoplastic fraction and osmotic potential [-]
% epsilon_L_max -- maximum elastic modulus [MPa]
% pi_L_0_25C -- osmotic potential at full hydration at standardized temperature of 25C [MPa]
% k_L_max_25C -- maximum soil-plant conductance per unit leaf area at 25C [mol m-2 s-1 MPa-1]
% A -- exponent for simple soil-plant conductance [MPa-1] -- k_L = k_L_max * exp(A*psi_L)
% Q_10_k_L -- Q10 value for the temperature-dependence of the maximum soil-plant conductance per unit leaf area [-]
% psi_soil -- soil water potential [MPa]
% T_a -- air temperature [C]
% RH -- relative humidity [-]
% P_atm -- atmospheric pressure [kPa]
% c_a -- air CO2 partial pressure [mol mol-1]
% R_abs -- absorbed shortwave radiation [W m-2]
% g_H_a -- boundary layer heat conductance [mol m-2 s-1]
% emiss_L -- leaf emissivity [-]
% V_cmax25_func -- maximum carboxylation capacity at 25C [mol m-2 s-1] expressed as function of leaf water potential

%% Default settings
PLC_max = 0.999;
dpsi_L_vect = 0.01; %[MPa]
R_abs_min = 10; %[W m-2]
Constant_chi_w = nan; %default to nan
Nonmonotonic_V_cmax_temperature_response = 0;
Electron_transport_limitation = 0;

%% Check ''varargin''
n_varargin = length(varargin);
n_varargin_name_and_value = floor(n_varargin/2);
if 2*n_varargin_name_and_value ~= n_varargin
    error('ERROR: varargin inputs should have an even total number of inputs: each name must has a value!')
end

for i = 1:n_varargin_name_and_value
   
    switch varargin{1+2*(i-1)}
        case 'Constant_chi_w'
            Constant_chi_w = varargin{2*i};
        case 'Nonmonotonic_V_cmax_temperature_response'
            Nonmonotonic_V_cmax_temperature_response = varargin{2*i};
            if isnumeric(Nonmonotonic_V_cmax_temperature_response)
                if (Nonmonotonic_V_cmax_temperature_response ~=0) && (Nonmonotonic_V_cmax_temperature_response ~= 1)
                    error('ERROR: ''Nonmonotonic_V_cmax_temperature_response'' must be a numeric value set to either 0 or 1!')
                end
            else
                error('ERROR: ''Nonmonotonic_V_cmax_temperature_response'' must be a numeric value set to either 0 or 1!')
            end
        case 'Electron_transport_limitation'
            Electron_transport_limitation = varargin{2*i};
            if isnumeric(Electron_transport_limitation)
                if (Electron_transport_limitation ~=0) && (Electron_transport_limitation ~= 1)
                    error('ERROR: ''Electron_transport_limitation'' must be a numeric value set to either 0 or 1!')
                end
            else
                error('ERROR: ''Electron_transport_limitation'' must be a numeric value set to either 0 or 1!')
            end
        otherwise
            if isstring(varargin{1+2*(i-1)}) || ischar(varargin{1+2*(i-1)})
                error(['ERROR: unsupported varargin name: ''', varargin{1+2*(i-1)}, '''', newline, 'Acceptable varargin names are ''Constant_chi_w'', ''Nonmonotonic_V_cmax_temperature_response'', and ''Electron_transport_limitation''.'])
            else
                error('ERROR: odd-numbered varargin inputs must be strings or characters!')
            end
    end
    
end

% determine discretization
psi_L_max = psi_soil;
psi_L_max = 0.1*ceil(psi_L_max/0.1); 
psi_L_min = log(1-PLC_max)/A;
psi_L_min = max(psi_L_min, psi_L_max-10);
psi_L_vect = psi_L_max:-abs(dpsi_L_vect):psi_L_min;
N_vect = length(psi_L_vect);

% determine T_L, VPD_L, and g_w from steady-state thermodynamics
T_a_vect = T_a * ones(1, N_vect);
RH_vect = RH * ones(1, N_vect);
P_atm_vect = P_atm * ones(1, N_vect);
c_a_vect = c_a * ones(1, N_vect);
pi_L_0_25C_vect = pi_L_0_25C .* ones(1, N_vect);
E_supply_vect = k_L_max_25C * Q_10_k_L.^((T_a_vect-25)/10) .* exp(A*psi_L_vect) .* (psi_soil - psi_L_vect); %supply curve
[T_L_vect, VPD_L_vect, g_w_supply_vect] = Steady_state_thermo(T_a_vect, E_supply_vect, RH_vect, P_atm_vect, R_abs, g_H_a, emiss_L);

% predict pressure-volume traits from previously estimated T_L and psi_L
[outputs_PV] = PV_from_psi_for_Trt( a_f_max, pi_L_star, beta, epsilon_L_max, ...
                                    pi_L_0_25C_vect, T_L_vect, psi_L_vect);

RWC_t_vect = outputs_PV.RWC_t;
RWC_s_vect = outputs_PV.RWC_s;
pi_L_vect = outputs_PV.pi_L;
pi_L_0_vect = outputs_PV.pi_L_0;
epsilon_L_0_vect = outputs_PV.epsilon_L_0;
a_f_vect = outputs_PV.a_f;
a_f_0_vect = outputs_PV.a_f_0;
da_f_dpi_L_vect = outputs_PV.da_f_dpi_L;
da_f_0dpi_L_0_vect = outputs_PV.da_f_0dpi_L_0;

% local temperature-dependent photosynthetic parameters
[outputs_photo_param] = Photosynthesis_parameters_temperature_response(T_L_vect, psi_L_vect, V_cmax25_func, ...
                                                                       'Nonmonotonic_V_cmax_temperature_response', Nonmonotonic_V_cmax_temperature_response, ...
                                                                       'Electron_transport_limitation', Electron_transport_limitation, ...
                                                                       'R_abs', R_abs);
Gamma_star_vect = outputs_photo_param.Gamma_star_vect;
K_m_vect = outputs_photo_param.K_m_vect;
V_cmax_vect = outputs_photo_param.V_cmax_vect;
R_d_vect = outputs_photo_param.R_d_vect;

% determine secondary estimate of g_w from previously estimated T_L and VPD_L
[outputs_stomata] = Stomata_for_RWC_and_Temperature(SWC_L, alpha, RWC_t_vect, RWC_s_vect, pi_L_vect, pi_L_0_vect, ...
                                                    a_f_vect, da_f_dpi_L_vect, a_f_0_vect, da_f_0dpi_L_0_vect, epsilon_L_0_vect, ...
                                                    VPD_L_vect, P_atm_vect, T_L_vect, V_cmax_vect, ...
                                                    R_d_vect, Gamma_star_vect, K_m_vect, c_a_vect, ...
                                                    'Constant_chi_w', Constant_chi_w);


g_w_demand_vect = outputs_stomata.g_w;
g_w_demand_vect(outputs_stomata.solved == 0) = nan;
chi_w_vect_demand = outputs_stomata.chi_w;


if Electron_transport_limitation == 1
    
    J_vect = outputs_photo_param.J_vect;
    
    [outputs_stomata_Electron_transport_limited] = Stomata_for_RWC_and_Temperature(SWC_L, alpha, RWC_t_vect, RWC_s_vect, pi_L_vect, pi_L_0_vect, ...
                                                                                   a_f_vect, da_f_dpi_L_vect, a_f_0_vect, da_f_0dpi_L_0_vect, epsilon_L_0_vect, ...
                                                                                   VPD_L_vect, P_atm_vect, T_L_vect, J_vect/4, ...
                                                                                   R_d_vect, Gamma_star_vect, 2*Gamma_star_vect, c_a_vect, ...
                                                                                   'Constant_chi_w', Constant_chi_w);
    g_w_demand_electron_transport_vect = outputs_stomata_Electron_transport_limited.g_w;
    g_w_demand_electron_transport_vect(outputs_stomata_Electron_transport_limited.solved == 0) = nan;
    if all(isnan(g_w_demand_electron_transport_vect)) && (R_abs < R_abs_min)
        g_w_demand_electron_transport_vect = zeros(1, N_vect);
    end
    chi_w_vect_3 = outputs_stomata_Electron_transport_limited.chi_w;
    
    g_w_demand_vect = min(g_w_demand_vect, g_w_demand_electron_transport_vect); 
    chi_w_vect_demand(g_w_demand_vect == g_w_demand_electron_transport_vect) = chi_w_vect_3(g_w_demand_vect == g_w_demand_electron_transport_vect);
    
end

E_demand_vect = g_w_demand_vect.*VPD_L_vect./P_atm_vect; %demand curve

% find intersection of two estimates of g_w
g_w_error_vect = g_w_supply_vect - g_w_demand_vect;
diff_sign_g_w_error_vect = diff(sign(g_w_error_vect));
diff_sign_g_w_error_vect(isnan(diff_sign_g_w_error_vect)) = 0;
not_zero_diff_sign_g_w_error_vect = (diff_sign_g_w_error_vect ~= 0);
if any(not_zero_diff_sign_g_w_error_vect)
    % interpolate where two estimates of g_w are equal
    convergence = 4; %convergence score -- 4 for complete convergence
    ind_LB = find(not_zero_diff_sign_g_w_error_vect, 1, 'first');
    ind_UB = ind_LB + 1;
    g_w_error_LB = g_w_error_vect(ind_LB);
    g_w_error_UB = g_w_error_vect(ind_UB);
    g_w_steady_state = g_w_supply_vect(ind_LB) + (g_w_supply_vect(ind_UB) - g_w_supply_vect(ind_LB)) * (0 - g_w_error_LB)/(g_w_error_UB - g_w_error_LB);
    g_w_steady_state = max(g_w_steady_state, 0);
    chi_w_steady_state = chi_w_vect_demand(ind_LB) + (chi_w_vect_demand(ind_UB) - chi_w_vect_demand(ind_LB)) * (0 - g_w_error_LB)/(g_w_error_UB - g_w_error_LB);
    VPD_L_steady_state = VPD_L_vect(ind_LB) + (VPD_L_vect(ind_UB) - VPD_L_vect(ind_LB)) * (0 - g_w_error_LB)/(g_w_error_UB - g_w_error_LB);
    T_L_steady_state = T_L_vect(ind_LB) + (T_L_vect(ind_UB) - T_L_vect(ind_LB)) * (0 - g_w_error_LB)/(g_w_error_UB - g_w_error_LB);
    E_steady_state = E_supply_vect(ind_LB) + (E_supply_vect(ind_UB) - E_supply_vect(ind_LB)) * (0 - g_w_error_LB)/(g_w_error_UB - g_w_error_LB);
    E_crit_steady_state = max(E_supply_vect);
    psi_L_steady_state = psi_L_vect(ind_LB) + (psi_L_vect(ind_UB) - psi_L_vect(ind_LB)) * (0 - g_w_error_LB)/(g_w_error_UB - g_w_error_LB);
    psi_L_crit_steady_state = psi_L_vect(E_supply_vect == E_crit_steady_state);
elseif ((max(g_w_demand_vect) < 0) && (min(g_w_supply_vect) < 0) && (max(g_w_supply_vect) > 0)) || (R_abs < R_abs_min)
    % set g_w = 0 and interpolate other traits
    convergence = 3; %convergence score
    ind_LB = find((g_w_supply_vect <= 0), 1, 'last');
    ind_UB = ind_LB + 1;
    g_w_LB = g_w_supply_vect(ind_LB);
    g_w_UB = g_w_supply_vect(ind_UB);
    g_w_steady_state = 0;
    chi_w_steady_state = chi_w_vect_demand(ind_LB) + (chi_w_vect_demand(ind_UB) - chi_w_vect_demand(ind_LB)) * (0 - g_w_LB)/(g_w_UB - g_w_LB);
    VPD_L_steady_state = VPD_L_vect(ind_LB) + (VPD_L_vect(ind_UB) - VPD_L_vect(ind_LB)) * (0 - g_w_LB)/(g_w_UB - g_w_LB);
    T_L_steady_state = T_L_vect(ind_LB) + (T_L_vect(ind_UB) - T_L_vect(ind_LB)) * (0 - g_w_LB)/(g_w_UB - g_w_LB);
    E_steady_state = 0;
    E_crit_steady_state = max(E_supply_vect);
    psi_L_steady_state = psi_L_vect(ind_LB) + (psi_L_vect(ind_UB) - psi_L_vect(ind_LB)) * (0 - g_w_LB)/(g_w_UB - g_w_LB);
    psi_L_crit_steady_state = psi_L_vect(E_supply_vect == E_crit_steady_state);
elseif min(abs(g_w_error_vect)) < 1e-2
    % partial convergence
    abs_g_w_error_vect = abs(g_w_error_vect);
    [~, ind_min_abs_error] = min(abs_g_w_error_vect);
    if E_supply_vect(ind_min_abs_error) >= 0
        % true partial convergence
        convergence = 2; %convergence score 
        g_w_steady_state = 0.5*(g_w_supply_vect(ind_min_abs_error) + g_w_demand_vect(ind_min_abs_error));
        chi_w_steady_state = chi_w_vect_demand(ind_min_abs_error);
        VPD_L_steady_state = VPD_L_vect(ind_min_abs_error);
        T_L_steady_state = T_L_vect(ind_min_abs_error);
        E_steady_state = E_supply_vect(ind_min_abs_error);
        E_crit_steady_state = max(E_supply_vect);
        psi_L_steady_state = psi_L_vect(ind_min_abs_error);
        psi_L_crit_steady_state = psi_L_vect(E_supply_vect == E_crit_steady_state);
    else
        % partial convergence but negative transpiration
        convergence = 1; %convergence score
        ind_LB = find((g_w_supply_vect <= 0), 1, 'last');
        ind_UB = ind_LB + 1;
        g_w_LB = g_w_supply_vect(ind_LB);
        g_w_UB = g_w_supply_vect(ind_UB);
        g_w_steady_state = 0;
        chi_w_steady_state = chi_w_vect_demand(ind_LB) + (chi_w_vect_demand(ind_UB) - chi_w_vect_demand(ind_LB)) * (0 - g_w_LB)/(g_w_UB - g_w_LB);
        VPD_L_steady_state = VPD_L_vect(ind_LB) + (VPD_L_vect(ind_UB) - VPD_L_vect(ind_LB)) * (0 - g_w_LB)/(g_w_UB - g_w_LB);
        T_L_steady_state = T_L_vect(ind_LB) + (T_L_vect(ind_UB) - T_L_vect(ind_LB)) * (0 - g_w_LB)/(g_w_UB - g_w_LB);
        E_steady_state = 0;
        E_crit_steady_state = max(E_supply_vect);
        psi_L_steady_state = psi_L_vect(ind_LB) + (psi_L_vect(ind_UB) - psi_L_vect(ind_LB)) * (0 - g_w_LB)/(g_w_UB - g_w_LB);
        psi_L_crit_steady_state = psi_L_vect(E_supply_vect == E_crit_steady_state);
    end
else
    convergence = 0; %convergence score
    ind_LB = find((g_w_supply_vect <= 0), 1, 'last');
    ind_UB = ind_LB + 1;
    g_w_LB = g_w_supply_vect(ind_LB);
    g_w_UB = g_w_supply_vect(ind_UB);
    g_w_steady_state = 0;
    chi_w_steady_state = chi_w_vect_demand(ind_LB) + (chi_w_vect_demand(ind_UB) - chi_w_vect_demand(ind_LB)) * (0 - g_w_LB)/(g_w_UB - g_w_LB);
    VPD_L_steady_state = VPD_L_vect(ind_LB) + (VPD_L_vect(ind_UB) - VPD_L_vect(ind_LB)) * (0 - g_w_LB)/(g_w_UB - g_w_LB);
    T_L_steady_state = T_L_vect(ind_LB) + (T_L_vect(ind_UB) - T_L_vect(ind_LB)) * (0 - g_w_LB)/(g_w_UB - g_w_LB);
    E_steady_state = 0;
    E_crit_steady_state = max(E_supply_vect);
    psi_L_steady_state = psi_L_vect(ind_LB) + (psi_L_vect(ind_UB) - psi_L_vect(ind_LB)) * (0 - g_w_LB)/(g_w_UB - g_w_LB);
    psi_L_crit_steady_state = psi_L_vect(E_supply_vect == E_crit_steady_state);
end

PLC_steady_state = 1 - exp(A*psi_L_steady_state);
g_c_steady_state = g_w_steady_state/1.6;
[outputs_photo_param_steady_state] = Photosynthesis_parameters_temperature_response(T_L_steady_state, psi_L_steady_state, V_cmax25_func, ...
                                                                       'Nonmonotonic_V_cmax_temperature_response', Nonmonotonic_V_cmax_temperature_response, ...
                                                                       'Electron_transport_limitation', Electron_transport_limitation, ...
                                                                       'R_abs', R_abs);
Gamma_star_steady_state = outputs_photo_param_steady_state.Gamma_star_vect;
K_m_steady_state = outputs_photo_param_steady_state.K_m_vect;
V_cmax_steady_state = outputs_photo_param_steady_state.V_cmax_vect;
R_d_steady_state = outputs_photo_param_steady_state.R_d_vect;
A_n_steady_state = (g_c_steady_state.*(c_a+K_m_steady_state)+V_cmax_steady_state-R_d_steady_state)/2 - (((g_c_steady_state.*(c_a-K_m_steady_state)-V_cmax_steady_state+R_d_steady_state)/2).^2 + g_c_steady_state.*(V_cmax_steady_state.*Gamma_star_steady_state+R_d_steady_state.*K_m_steady_state) + g_c_steady_state.^2.*c_a.*K_m_steady_state).^0.5;
A_c_limited_steady_state = 1;
J_steady_state = nan;
J_max_steady_state = nan;

if Electron_transport_limitation == 1
    
    J_steady_state = outputs_photo_param_steady_state.J_vect;
    J_max_steady_state = outputs_photo_param_steady_state.J_max_vect;
    A_j_steady_state = (g_c_steady_state.*(c_a+2*Gamma_star_steady_state)+J_steady_state/4-R_d_steady_state)/2 - (((g_c_steady_state.*(c_a-2*Gamma_star_steady_state)-J_steady_state/4+R_d_steady_state)/2).^2 + g_c_steady_state.*(J_steady_state/4.*Gamma_star_steady_state+R_d_steady_state.*2*Gamma_star_steady_state) + g_c_steady_state.^2.*c_a.*2*Gamma_star_steady_state).^0.5;
    A_n_steady_state = min(A_n_steady_state, A_j_steady_state); 
    A_c_limited_steady_state(A_n_steady_state == A_j_steady_state) = 0;
    
end
c_i_steady_state = c_a - A_n_steady_state/g_c_steady_state; 


%% Store outputs
outputs_stomata_steady_state.g_w_steady_state = g_w_steady_state;
outputs_stomata_steady_state.chi_w_steady_state = chi_w_steady_state;
outputs_stomata_steady_state.VPD_L_steady_state = VPD_L_steady_state;
outputs_stomata_steady_state.T_L_steady_state = T_L_steady_state;
outputs_stomata_steady_state.E_steady_state = E_steady_state;
outputs_stomata_steady_state.E_crit_steady_state = E_crit_steady_state;
outputs_stomata_steady_state.psi_L_steady_state = psi_L_steady_state;
outputs_stomata_steady_state.psi_L_crit_steady_state = psi_L_crit_steady_state;
outputs_stomata_steady_state.PLC_steady_state = PLC_steady_state;
outputs_stomata_steady_state.A_n_steady_state = A_n_steady_state;
outputs_stomata_steady_state.c_i_steady_state = c_i_steady_state;
outputs_stomata_steady_state.V_cmax_steady_state = V_cmax_steady_state;
outputs_stomata_steady_state.J_steady_state = J_steady_state;
outputs_stomata_steady_state.J_max_steady_state = J_max_steady_state;
outputs_stomata_steady_state.A_c_limited_steady_state = A_c_limited_steady_state;
outputs_stomata_steady_state.convergence = convergence;
outputs_stomata_steady_state.psi_L_vect = psi_L_vect;
outputs_stomata_steady_state.E_supply_vect = E_supply_vect;
outputs_stomata_steady_state.E_demand_vect = E_demand_vect;
outputs_stomata_steady_state.VPD_L_vect = VPD_L_vect;
outputs_stomata_steady_state.g_w_supply_vect = g_w_supply_vect;

end

