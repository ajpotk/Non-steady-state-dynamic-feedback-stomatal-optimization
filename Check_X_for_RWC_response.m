function [Check] = Check_X_for_RWC_response(X, pi_L_0_25C)

% % % % desired range of slope between psi_L and chi_w
% % % beta_0_min = -2; %minimum of -1.62 in Manzoni et al. (2011)
% % % beta_0_max = -0.05; %maximum of -0.13 in Manzoni et al. (2011)

% desired ratio of g_w at RWC_t = 0.85 to g_w at RWC_t = 1.00
g_w_ratio_85_to_100_RWC_t_max = 0.15;

% Constants
o_i = 1e-3 * 210; %standard atmospheric partial pressure of oxygen [mol mol-1]
R = 8.314; %universal gas constant [J mol-1 K-1]
R_d25_per_V_cmax25 = 0.015; %from Collatz et al. (1991)

% default environmental conditions
P_atm = 101.325; %atmospheric pressure [kPa]
VPD_L = 1; %leaf-to-air vapor pressure deficit [kPa]
c_a = 4e-4; %atmospheric CO2 conc

% realistic conditions for RWC_t, T_L, and pi_L_0
RWC_t = [0.85, 1]; %0.85:0.01:1.00;
T_L = 25; %code below is not supported when T_L is an array

% reshape connditions for their interactions
N_RWC_t = length(RWC_t);
N_T_L = length(T_L);
N_pi_L_0_25C = length(pi_L_0_25C);
N_tot = N_RWC_t * N_T_L * N_pi_L_0_25C;
[RWC_t, T_L, pi_L_0_25C] = meshgrid(RWC_t, T_L, pi_L_0_25C);

if N_T_L > 1
    error('ERROR: ''T_L'' must not be an array!')
elseif N_pi_L_0_25C > 1
    error('ERROR: ''pi_L_0_25C'' must not be an array!')
elseif N_pi_L_0_25C < 1
    error('ERROR: ''pi_L_0_25C'' must not be empty!')
end

% reshape as column vectors for functions below
RWC_t = reshape(RWC_t, N_tot, 1, 1); 
T_L = reshape(T_L, N_tot, 1, 1);
pi_L_0_25C = reshape(pi_L_0_25C, N_tot, 1, 1);
P_atm = repmat(P_atm, N_tot, 1);
VPD_L = repmat(VPD_L, N_tot, 1);
c_a = repmat(c_a, N_tot, 1);

% load other data
load_data_file_name = 'data_subset.mat';
data_loaded = 0;
i = 0;
while data_loaded == 0
    pause(0)
    i = i + 1;
    if i > 100
        break
    end
    try
       load(load_data_file_name, 'T_L_subset', 'V_cmax_subset')
       data_loaded = 1;
    catch ME
    end
end
if data_loaded == 0
    rethrow(ME)
end
% % % load(load_data_file_name, 'T_L_subset', 'V_cmax_subset')
% % % load(load_data_file_name, 'Trt_numb_subset', 'T_L_subset', ...
% % %          'psi_L_MD_subset', 'K_m_subset', 'Gamma_star_subset', ...
% % %          'V_cmax_subset', 'R_d_subset', 'c_a_subset', 'VPD_L_subset', ...
% % %          'P_atm_subset', 'g_w_subset')


% determine mean photosynthetic parameters
V_cmax25_func = @(psi_L_MD) 0; 
[outputs_photo_param] = Photosynthesis_parameters_temperature_response(T_L_subset, 0, V_cmax25_func);
V_cmax_per_V_cmax25_subset = outputs_photo_param.V_cmax_per_V_cmax25_vect;
V_cmax25 = mean(V_cmax_subset ./ V_cmax_per_V_cmax25_subset, 'omitnan');
V_cmax25_func = @(psi_L_MD) V_cmax25; 

% determine photosynthetic parameters for T_L
[outputs_photo_param] = Photosynthesis_parameters_temperature_response(T_L, 0, V_cmax25_func);
Gamma_star = outputs_photo_param.Gamma_star_vect;
K_m = outputs_photo_param.K_m_vect;
V_cmax = outputs_photo_param.V_cmax_vect;
R_d = outputs_photo_param.R_d_vect;

% get parameters
[a_f_max, pi_L_star, beta, epsilon_L_max, SWC_L, alpha] = X_to_parameters(X);

% determine pressure-volume traits
[outputs_PV] = PV_from_RWC_t_for_Trt( a_f_max, pi_L_star, beta, epsilon_L_max, ...
                                       pi_L_0_25C, T_L, RWC_t);

RWC_s = outputs_PV.RWC_s;
pi_L = outputs_PV.pi_L;
pi_L_0 = outputs_PV.pi_L_0;
epsilon_L_0 = outputs_PV.epsilon_L_0;
a_f = outputs_PV.a_f;
a_f_0 = outputs_PV.a_f_0;
da_f_dpi_L = outputs_PV.da_f_dpi_L;
da_f_0dpi_L_0 = outputs_PV.da_f_0dpi_L_0;
psi_L = outputs_PV.psi_L;

% determine stomatal conductance
[outputs_stomata] = Stomata_for_RWC_and_Temperature(SWC_L, alpha, RWC_t, RWC_s, pi_L, pi_L_0, ...
                                                    a_f, da_f_dpi_L, a_f_0, da_f_0dpi_L_0, epsilon_L_0, ...
                                                    VPD_L, P_atm, T_L, V_cmax, ...
                                                    R_d, Gamma_star, K_m, c_a);

g_w = outputs_stomata.g_w; 
solved = outputs_stomata.solved;

g_w_ratio_85_to_100_RWC_t = g_w(RWC_t == 0.85)/g_w(RWC_t == 1);
Check = (g_w_ratio_85_to_100_RWC_t <= g_w_ratio_85_to_100_RWC_t_max) && all(solved);

% % % chi_w = outputs_stomata.chi_w;
% % % chi_w(chi_w <= 0) = nan;
% % % p = polyfit(psi_L, log(chi_w), 1);
% % % beta_0 = p(1); %slope between psi_L and log(chi_w), based on Manzoni et al. (2011)
% % % % % % chi_w_ww = exp(p(2)); %well-watered chi_w, based on Manzoni et al. (2011)
% % % Check = all(chi_w > 0) && (beta_0 >= beta_0_min) && (beta_0 <= beta_0_max) && all(solved);

end

