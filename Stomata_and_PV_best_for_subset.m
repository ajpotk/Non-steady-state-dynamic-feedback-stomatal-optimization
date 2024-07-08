function [output_subset] = Stomata_and_PV_best_for_subset(x)

%% Stomata_and_PV_best_for_subset
% written by Aaron Potkay 2023

%% Inputs:
%   x(1) = a_f_max -- maximum apoplasm fraction at an osmotic potential of zero [m3 m-3]
%   x(2) = pi_L_star -- osmotic potential at which a_f = exp(-1)*a_f_max = 0.37*a_f_max [MPa]
%   x(3) = beta -- expontent in Weibull function describing apoplasm fraction, a_f [-]
%   x(4) = epsilon_L_max -- theoretical maximum elastic modulus [MPa]
%   x(5) = SWC_L -- saturated water content of leaf [kg kg-1]

%% Outputs:
%	RWC_s -- predicted symplast relative water content [m3 m-3]
%	pi_L -- predicted osmotic potential [MPa]
%	pi_L_0_25C -- predicted osmotic potential at full hydration at a reference temperature of 25C [MPa]
%	P_L -- predicted turgor [MPa]
%	a_f_0 -- predicted apoplasm fraction at full hydration [m3 m-3]
%	a_f -- predicted apoplasm fraction [m3 m-3]
%	RWC_t -- predicted total/bulk relative water content [m3 m-3]
%	epsilon_L_0 -- predicted elastic modulus at full hydration [MPa]
%	chi_w_0 -- predicted marginal water-use efficiency at full hydration [mol mol-1]
%	chi_w -- predicted marginal water-use efficiency [mol mol-1]
%   g_w -- stomatal conductance to H2O [mol m-2 s-2]
%   pi_L_0_25C_Trt_combo -- ''N_Trt_combo'' values of osmotic potential at full hydration at a reference temperature of 25C for each treamtment combo [MPa]

%% Retrieve parameters
[a_f_max, pi_L_star, beta, epsilon_L_max, SWC_L, alpha] = X_to_parameters(x);

%% Load other data
load_data_file_name = 'data_subset.mat';
load(load_data_file_name, 'Trt_numb_subset', 'T_L_subset', ...
     'psi_L_MD_subset', 'K_m_subset', 'Gamma_star_subset', ...
     'V_cmax_subset', 'R_d_subset', 'c_a_subset', 'VPD_L_subset', ...
     'P_atm_subset', 'g_w_subset')

%% Calculate predictions
N_data = length(Trt_numb_subset);
N_Trt_combo = max(Trt_numb_subset);
N_test = 8e2; %number of potential pi_L_0 values that will be tested

pi_L_0_25C_Trt_combo = nan(1,N_Trt_combo);
SSE_best_Trt_combo = nan(1,N_Trt_combo);

chi_w_subset_pred = nan(N_data, 1);
RWC_s_subset_pred = nan(N_data, 1);
RWC_s_tlp_subset_pred = nan(N_data, 1);
RWC_t_subset_pred = nan(N_data, 1);
pi_L_subset_pred = nan(N_data, 1);
pi_L_0_subset_pred = nan(N_data, 1);
P_L_subset_pred = nan(N_data, 1);
a_f_0_subset_pred = nan(N_data, 1);
a_f_subset_pred = nan(N_data, 1);
epsilon_L_0_subset_pred = nan(N_data, 1);
solved_subset_pred = nan(N_data, 1);
g_w_subset_pred = nan(N_data, 1);
pi_L_min_subset_pred = nan(N_data, 1);
pi_L_max_subset_pred = nan(N_data, 1);

NumWorkers = 4; %25
distcomp.feature( 'LocalUseMpiexec', false );
% Check if parallel processor is on
p = gcp('nocreate');
if isempty(p)
    try
        % Set up parallel processing if off
        local = parcluster('local');
        local.NumWorkers = NumWorkers;
        parpool('local', NumWorkers);
    catch ME
        disp(ME.message)
        error('ERROR: Could not establish parallel pool!')
    end
end

for i = unique(Trt_numb_subset')
    
    % divide data by treatment
    T_L_Trt = T_L_subset(Trt_numb_subset == i);
    psi_L_MD_Trt = psi_L_MD_subset(Trt_numb_subset == i);
    K_m_Trt = K_m_subset(Trt_numb_subset == i);
    Gamma_star_Trt = Gamma_star_subset(Trt_numb_subset == i);
    V_cmax_Trt = V_cmax_subset(Trt_numb_subset == i);
    R_d_Trt = R_d_subset(Trt_numb_subset == i);
    c_a_Trt = c_a_subset(Trt_numb_subset == i);
    VPD_L_Trt = VPD_L_subset(Trt_numb_subset == i);
    P_atm_Trt = P_atm_subset(Trt_numb_subset == i);
    g_w_Trt = g_w_subset(Trt_numb_subset == i);
    
    % generate potential pi_L_0_25C values that will be tested
    pi_L_0_25C_Trt_test = linspace(-1, -3, N_test); %starts at least-negative values
    
    % calculate SSE for each potential pi_L_0 value
    SSE = nan(1,N_test);
    parfor j = 1:N_test
        
        pi_L_0_25C = pi_L_0_25C_Trt_test(j);
        
        [output_Trt] = Stomata_and_PV_for_Trt(a_f_max, pi_L_star, beta, epsilon_L_max, SWC_L, alpha, ...
                                               pi_L_0_25C, psi_L_MD_Trt, VPD_L_Trt, ...
                                               P_atm_Trt, T_L_Trt, V_cmax_Trt, ...
                                               R_d_Trt, Gamma_star_Trt, K_m_Trt, ...
                                               c_a_Trt);
        
        g_w_TrT_pred = output_Trt.g_w;
        g_w_TrT_pred(output_Trt.solved == 0) = 9999;
        
        Check = Check_X_for_RWC_response(x, pi_L_0_25C);
        if Check == 0
            g_w_TrT_pred = 9999*ones(size(g_w_TrT_pred));
        end
        
        SSE(j) = sum((g_w_Trt - g_w_TrT_pred).^2);
        
    end
    
    % find best value of pi_L_0_25C with lowest SSE
    SSE_best = min(SSE);
    ind_best = find(SSE == SSE_best, 1, 'first');
    pi_L_0_25C = pi_L_0_25C_Trt_test(ind_best);
    pi_L_0_25C_Trt_combo(i) = pi_L_0_25C;
    SSE_best_Trt_combo(i) = SSE_best;
    
    % solve again with best value of pi_L_0
    [output_Trt] = Stomata_and_PV_for_Trt(a_f_max, pi_L_star, beta, epsilon_L_max, SWC_L, alpha, ...
                                          pi_L_0_25C, psi_L_MD_Trt, VPD_L_Trt, ...
                                          P_atm_Trt, T_L_Trt, V_cmax_Trt, ...
                                          R_d_Trt, Gamma_star_Trt, K_m_Trt, ...
                                          c_a_Trt);
    g_w_TrT_pred = output_Trt.g_w;
	g_w_TrT_pred(output_Trt.solved == 0) = 9999;
    
    % fill in vectors of predicted values with best value of pi_L_0
    chi_w_subset_pred(Trt_numb_subset == i) = output_Trt.chi_w;
    RWC_s_subset_pred(Trt_numb_subset == i) = output_Trt.RWC_s;
    RWC_s_tlp_subset_pred(Trt_numb_subset == i) = output_Trt.RWC_s_tlp;
    RWC_t_subset_pred(Trt_numb_subset == i) = output_Trt.RWC_t;
    pi_L_subset_pred(Trt_numb_subset == i) = output_Trt.pi_L;
    pi_L_0_subset_pred(Trt_numb_subset == i) = output_Trt.pi_L_0;
    P_L_subset_pred(Trt_numb_subset == i) = output_Trt.P_L;
    a_f_0_subset_pred(Trt_numb_subset == i) = output_Trt.a_f_0;
    a_f_subset_pred(Trt_numb_subset == i) = output_Trt.a_f;
    epsilon_L_0_subset_pred(Trt_numb_subset == i) = output_Trt.epsilon_L_0;
    solved_subset_pred(Trt_numb_subset == i) = output_Trt.solved;
    g_w_subset_pred(Trt_numb_subset == i) = g_w_TrT_pred;
    pi_L_min_subset_pred(Trt_numb_subset == i) = output_Trt.pi_L_min;
    pi_L_max_subset_pred(Trt_numb_subset == i) = output_Trt.pi_L_max;
    
end

pi_L_min_combo = min(pi_L_min_subset_pred, [], 'all');
pi_L_max_combo = min(pi_L_max_subset_pred, [], 'all');

% fill in output structure
output_subset.N_Trt_combo = N_Trt_combo;
output_subset.pi_L_min_combo = pi_L_min_combo;
output_subset.pi_L_max_combo = pi_L_max_combo;
output_subset.SSE_best_total_subset = sum(SSE_best_Trt_combo, 'omitnan');
output_subset.SSE_best_Trt_combo = SSE_best_Trt_combo;
output_subset.pi_L_0_25C_Trt_combo = pi_L_0_25C_Trt_combo;
output_subset.chi_w = chi_w_subset_pred;
output_subset.RWC_s = RWC_s_subset_pred;
output_subset.RWC_s_tlp = RWC_s_tlp_subset_pred;
output_subset.RWC_t = RWC_t_subset_pred;
output_subset.pi_L = pi_L_subset_pred;
output_subset.pi_L_0 = pi_L_0_subset_pred;
output_subset.P_L = P_L_subset_pred;
output_subset.a_f_0 = a_f_0_subset_pred;
output_subset.a_f = a_f_subset_pred;
output_subset.epsilon_L_0 = epsilon_L_0_subset_pred;
output_subset.solved = solved_subset_pred;
output_subset.g_w = g_w_subset_pred;


end

