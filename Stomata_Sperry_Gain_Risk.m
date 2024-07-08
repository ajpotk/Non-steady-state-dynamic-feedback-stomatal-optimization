function [outputs_stomata_Sperry] = ...
         Stomata_Sperry_Gain_Risk(k_L_max_25C, A, Q_10_k_L, psi_soil, ...
                                  T_a, RH, P_atm, c_a, ...
                                  R_abs, g_H_a, emiss_L, ...
                                  V_cmax25_func, ...
                                  varargin)

%% Default settings
PLC_max = 0.999;
dpsi_L_vect = 0.01; %[MPa]
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
                error(['ERROR: unsupported varargin name: ''', varargin{1+2*(i-1)}, '''', newline, 'Acceptable varargin names are ''Nonmonotonic_V_cmax_temperature_response'' and ''Electron_transport_limitation''.'])
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
psi_L_vect = unique([psi_L_max:-abs(dpsi_L_vect):psi_L_min, psi_soil]);
N_vect = length(psi_L_vect);
T_a_vect = T_a * ones(1, N_vect);
E_vect = k_L_max_25C * Q_10_k_L.^((T_a_vect-25)/10) .* exp(A*psi_L_vect) .* (psi_soil - psi_L_vect);
E_crit = max(E_vect);
psi_L_crit = psi_L_vect(E_vect == E_crit);
is_allow = ((E_vect >= 0) + (psi_L_vect > psi_L_crit) == 2);
E_vect = E_vect(is_allow == 1);
psi_L_vect = psi_L_vect(is_allow == 1);
T_a_vect = T_a_vect(is_allow == 1);
N_vect = length(psi_L_vect);

RH_vect = RH * ones(1, N_vect);
P_atm_vect = P_atm * ones(1, N_vect);
c_a_vect = c_a * ones(1, N_vect);

% determine T_L, VPD_L, and g_w from steady-state thermodynamics
[T_L_vect, VPD_L_vect, g_w_vect] = Steady_state_thermo(T_a_vect, E_vect, RH_vect, P_atm_vect, R_abs, g_H_a, emiss_L);
g_w_vect(g_w_vect < 0) = nan;

% Sperry et al.'s (2017) xylem hydraulics
k_c_vect = -k_L_max_25C*Q_10_k_L.^((T_a_vect-25)/10).*exp(A*psi_L_vect).*(A*(psi_soil - psi_L_vect) - 1);
k_c_max = max(k_c_vect, [], 'omitnan');
dk_cdE_vect = -A*(A*(psi_L_vect - psi_soil) + 2)./(A*(psi_L_vect - psi_soil) + 1);

% local temperature-dependent photosynthetic parameters
[outputs_photo_param] = Photosynthesis_parameters_temperature_response(T_L_vect, psi_L_vect, V_cmax25_func, ...
                                                                       'Nonmonotonic_V_cmax_temperature_response', Nonmonotonic_V_cmax_temperature_response, ...
                                                                       'Electron_transport_limitation', Electron_transport_limitation, ...
                                                                       'R_abs', R_abs);
Gamma_star_vect = outputs_photo_param.Gamma_star_vect;
K_m_vect = outputs_photo_param.K_m_vect;
V_cmax_vect = outputs_photo_param.V_cmax_vect;
R_d_vect = outputs_photo_param.R_d_vect;

% solve for A_n
g_c_vect = g_w_vect/1.6;
A_n_vect = (g_c_vect.*(c_a_vect+K_m_vect)+V_cmax_vect-R_d_vect)/2 - (((g_c_vect.*(c_a_vect-K_m_vect)-V_cmax_vect+R_d_vect)/2).^2 + g_c_vect.*(V_cmax_vect.*Gamma_star_vect+R_d_vect.*K_m_vect) + g_c_vect.^2.*c_a_vect.*K_m_vect).^0.5;
A_c_limited_vect = ones(1, N_vect);

if Electron_transport_limitation == 1
    
    J_vect = outputs_photo_param.J_vect;
    J_max_vect = outputs_photo_param.J_max_vect;
    A_j_vect = (g_c_vect.*(c_a_vect+2*Gamma_star_vect)+J_vect/4-R_d_vect)/2 - (((g_c_vect.*(c_a_vect-2*Gamma_star_vect)-J_vect/4+R_d_vect)/2).^2 + g_c_vect.*(J_vect/4.*Gamma_star_vect+R_d_vect*2.*Gamma_star_vect) + g_c_vect.^2.*c_a_vect*2.*Gamma_star_vect).^0.5;
    A_n_vect = min(A_n_vect, A_j_vect); 
    A_c_limited_vect(A_n_vect == A_j_vect) = 0;
    
end

A_n_vect(isnan(psi_L_vect)) = nan;
A_n_max = max(A_n_vect, [], 'omitnan');
chi_w_vect = -A_n_max/k_c_max*dk_cdE_vect;

if A_n_max > 0
    
    % Sperry et al.'s (2017) objective
    Sperry_objective = A_n_vect + A_n_max*k_c_vect/k_c_max;
    Sperry_objective_max = max(Sperry_objective, [], 'omitnan');

    % find optimal g_w
    ind_opt = find((Sperry_objective == Sperry_objective_max), 1, 'first');

else
    
    % stomatal closure
    ind_opt = find((g_w_vect == 0), 1, 'first');
    
end

g_w_Sperry = g_w_vect(ind_opt);
VPD_L_Sperry = VPD_L_vect(ind_opt);
T_L_Sperry = T_L_vect(ind_opt);
E_Sperry = E_vect(ind_opt);
A_n_Sperry = A_n_vect(ind_opt);
chi_w_Sperry = chi_w_vect(ind_opt);
c_i_Sperry = c_a - 1.6*A_n_Sperry/g_w_Sperry; 
psi_L_Sperry = psi_L_vect(ind_opt);
PLC_Sperry = 1 - exp(A*psi_L_Sperry);
V_cmax_Sperry = V_cmax_vect(ind_opt);
if Electron_transport_limitation == 1
    J_Sperry = J_vect(ind_opt);
    J_max_Sperry = J_max_vect(ind_opt);
else
    J_Sperry = nan;
    J_max_Sperry = nan;
end
A_c_limited_Sperry = A_c_limited_vect(ind_opt);

%% Store outputs
outputs_stomata_Sperry.g_w_Sperry = g_w_Sperry;
outputs_stomata_Sperry.VPD_L_Sperry = VPD_L_Sperry;
outputs_stomata_Sperry.T_L_Sperry = T_L_Sperry;
outputs_stomata_Sperry.E_Sperry = E_Sperry;
outputs_stomata_Sperry.A_n_Sperry = A_n_Sperry;
outputs_stomata_Sperry.chi_w_Sperry = chi_w_Sperry;
outputs_stomata_Sperry.c_i_Sperry = c_i_Sperry;
outputs_stomata_Sperry.psi_L_Sperry = psi_L_Sperry;
outputs_stomata_Sperry.PLC_Sperry = PLC_Sperry;
outputs_stomata_Sperry.V_cmax_Sperry = V_cmax_Sperry;
outputs_stomata_Sperry.J_Sperry = J_Sperry;
outputs_stomata_Sperry.J_max_Sperry = J_max_Sperry;
outputs_stomata_Sperry.A_c_limited_Sperry = A_c_limited_Sperry;

end

