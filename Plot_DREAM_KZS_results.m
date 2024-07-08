clear
clc
close all

%% Plot_DREAM_KZS_results
% written by Aaron Potkay 2023

% remove old files -- see ''Plot_g_w.m'' code
All_spec_fig_save_name = 'ALL_SPECIES_Fig_gw_best_layered';
All_spec_fig_save_name2 = 'ALL_SPECIES_Fig_gw_best_panel';
spec_ind_file_save_name = 'ALL_SPECIES-index.mat';
spec_ind_file_save_name2 = 'ALL_SPECIES-index2.mat';
All_spec_fig_save_name3 = 'ALL_SPECIES_Fig_violin_pi_L_0';
spec_ind_file_save_name3 = 'ALL_SPECIES-index-3.mat';
All_spec_fig_save_name4 = 'ALL_SPECIES_Fig_a_f-psi_L';
spec_ind_file_save_name4 = 'ALL_SPECIES-index-4.mat';
All_spec_fig_save_name5 = 'ALL_SPECIES_Fig_gw_Medlyn_USO_panel';
spec_ind_file_save_name5 = 'ALL_SPECIES-index-5.mat';
files = struct2cell(dir);
files = files(1,:);
is_All_spec_fig = cell2mat(cellfun(@(x) strcmp(x, [All_spec_fig_save_name, '.fig']), files, 'UniformOutput', 0));
is_All_spec_fig2 = cell2mat(cellfun(@(x) strcmp(x, [All_spec_fig_save_name2, '.fig']), files, 'UniformOutput', 0));
is_All_spec_fig3 = cell2mat(cellfun(@(x) strcmp(x, [All_spec_fig_save_name3, '.fig']), files, 'UniformOutput', 0));
is_All_spec_fig4 = cell2mat(cellfun(@(x) strcmp(x, [All_spec_fig_save_name4, '.fig']), files, 'UniformOutput', 0));
is_All_spec_fig5 = cell2mat(cellfun(@(x) strcmp(x, [All_spec_fig_save_name5, '.fig']), files, 'UniformOutput', 0));
is_spec_ind_file = cell2mat(cellfun(@(x) strcmp(x, spec_ind_file_save_name), files, 'UniformOutput', 0));
is_spec_ind_file2 = cell2mat(cellfun(@(x) strcmp(x, spec_ind_file_save_name2), files, 'UniformOutput', 0));
is_spec_ind_file3 = cell2mat(cellfun(@(x) strcmp(x, spec_ind_file_save_name3), files, 'UniformOutput', 0));
is_spec_ind_file4 = cell2mat(cellfun(@(x) strcmp(x, spec_ind_file_save_name4), files, 'UniformOutput', 0));
is_spec_ind_file5 = cell2mat(cellfun(@(x) strcmp(x, spec_ind_file_save_name5), files, 'UniformOutput', 0));

if any(is_All_spec_fig)
    delete([All_spec_fig_save_name, '.fig'])
end
if any(is_All_spec_fig2)
    delete([All_spec_fig_save_name2, '.fig'])
end
if any(is_All_spec_fig3)
    delete([All_spec_fig_save_name3, '.fig'])
end
if any(is_All_spec_fig4)
    delete([All_spec_fig_save_name4, '.fig'])
end
if any(is_All_spec_fig5)
    delete([All_spec_fig_save_name5, '.fig'])
end
if any(is_spec_ind_file)
    delete(spec_ind_file_save_name)
end
if any(is_spec_ind_file2)
    delete(spec_ind_file_save_name2)
end
if any(is_spec_ind_file3)
    delete(spec_ind_file_save_name3)
end
if any(is_spec_ind_file4)
    delete(spec_ind_file_save_name4)
end
if any(is_spec_ind_file5)
    delete(spec_ind_file_save_name5)
end

% retrieve observed data
[Var] = Download_B4W_data();
[Var] = Download_Sapes_data(Var); 
[Var] = Download_Lotschental_data(Var); 

% retrieve structured data
T_L = Var.T_L;
K_c = Var.K_c;
K_o = Var.K_o;
K_m = Var.K_m;
Gamma_star = Var.Gamma_star;
c_a = Var.c_a;
c_i = Var.c_i;
V_cmax = Var.V_cmax;
R_d = Var.R_d;
A_n = Var.A_n;
E = Var.E;
g_w = Var.g_w;
lambda = Var.lambda;
psi_L_MD = Var.psi_L_MD;
psi_L_PD = Var.psi_L_PD;
VPD_L = Var.VPD_L;
P_atm = Var.P_atm; 
Trt_numb = Var.Trt_numb;
Comb_Trt_name = Var.Comb_Trt_name;
SWC_L = Var.SWC_L;
species = Var.species;

%% Subset data by species
species_unique = unique(species);
N_species_unique = length(species_unique);
M_data = length(g_w);

% identify nans
is_not_nan = (~isnan(g_w)) .* ...
             (~isnan(T_L)) .* ...
             (~isnan(A_n)) .* ...
             (~isnan(lambda)) .* ...
             (~isnan(psi_L_MD)) .* ...
             (~isnan(V_cmax));

for i = [4, 1, 2, 3, 5] %1:N_species_unique %order changed for desired order of panels in some figures
    
    species_subset = species_unique{i}; %choice of species
    disp(['Plotting DREAM KZS  results for species #', num2str(i), ' of ' num2str(N_species_unique), ': ', species_subset])
    disp('    Subsetting data by species')
    is_species = cell2mat(cellfun(@(x) strcmp(x, species_subset), species, 'UniformOutput', 0)); 
    
    is_species_and_not_nan = is_species .* is_not_nan; 
    ind_is_species_and_not_nan = is_species_and_not_nan .* reshape(1:M_data, size(is_species));
    ind_is_species_and_not_nan = ind_is_species_and_not_nan(ind_is_species_and_not_nan > 0);
    M_data_species = length(ind_is_species_and_not_nan);

    % species-specific data
    T_L_subset = T_L(ind_is_species_and_not_nan);
    K_c_subset = K_c(ind_is_species_and_not_nan);
    K_o_subset = K_o(ind_is_species_and_not_nan);
    K_m_subset = K_m(ind_is_species_and_not_nan);
    Gamma_star_subset = Gamma_star(ind_is_species_and_not_nan);
    c_a_subset = c_a(ind_is_species_and_not_nan);
    c_i_subset = c_i(ind_is_species_and_not_nan);
    V_cmax_subset = V_cmax(ind_is_species_and_not_nan);
    R_d_subset = R_d(ind_is_species_and_not_nan);
    A_n_subset = A_n(ind_is_species_and_not_nan);
    E_subset = E(ind_is_species_and_not_nan);
    g_w_subset = g_w(ind_is_species_and_not_nan);
    lambda_subset = lambda(ind_is_species_and_not_nan);
    psi_L_MD_subset = psi_L_MD(ind_is_species_and_not_nan);
    psi_L_PD_subset = psi_L_PD(ind_is_species_and_not_nan);
    VPD_L_subset = VPD_L(ind_is_species_and_not_nan);
    P_atm_subset = P_atm(ind_is_species_and_not_nan); 
    Trt_numb_subset = Trt_numb(ind_is_species_and_not_nan);
    Comb_Trt_name_subset = Comb_Trt_name(ind_is_species_and_not_nan);
    SWC_L_subset = SWC_L(ind_is_species_and_not_nan);
    
    % save local .mat file for subset of empirical data
    load_data_file_name = 'data_subset.mat';
    save(load_data_file_name, ...
         'T_L_subset', 'K_c_subset', 'K_o_subset', 'K_m_subset', ...
         'Gamma_star_subset', 'c_a_subset', 'c_i_subset', ...
         'V_cmax_subset', 'R_d_subset', 'A_n_subset', 'E_subset',  ...
         'g_w_subset', 'lambda_subset', 'psi_L_MD_subset', ...
         'psi_L_PD_subset', 'VPD_L_subset', 'P_atm_subset', ...
         'Trt_numb_subset', 'Comb_Trt_name_subset',  'SWC_L_subset')
    
    %% Number of subsets
    N_Trt_combo = max(Trt_numb_subset);
    
    %% Load past DREAM results
    DREAM_results_file = ['DREAM_KZS_results_for_', species_subset, '.mat'];
    load(DREAM_results_file)
    
    chain_kzs_best = chain_kzs(1:end, 1:Nx);
    
    DREAM_Output_file = ['Output_file_for_', species_subset, '.txt'];
    fileID = fopen(DREAM_Output_file, 'r');
    formatSpec = '%f, ';
    formatSpec = repmat(formatSpec, 1, 10+2*N_Trt_combo);
    DREAM_Output = textscan(fileID, formatSpec, 'HeaderLines', 1);
    fclose(fileID);
    
    a_f_max_results = DREAM_Output{1};
    pi_L_star_results = DREAM_Output{2};
    beta_results = DREAM_Output{3};
    epsilon_L_max_results = DREAM_Output{4};
    SWC_L_results = DREAM_Output{5};
    alpha_results = DREAM_Output{6};
    SSE_best_total_subset_results = DREAM_Output{7};
    N_Trt_combo_results = DREAM_Output{8};
    N_results = length(a_f_max_results);
    SSE_best_Trt_combo_results = nan(N_results, N_Trt_combo);
    pi_L_0_25C_Trt_combo_results = nan(N_results, N_Trt_combo);
    for j = 1:N_Trt_combo
        SSE_best_Trt_combo_results(:,j) = DREAM_Output{j+8};
        pi_L_0_25C_Trt_combo_results(:,j) = DREAM_Output{j+8+N_Trt_combo};
    end
    pi_L_min_combo_results = DREAM_Output{9+2*N_Trt_combo};
    pi_L_max_combo_results = DREAM_Output{10+2*N_Trt_combo};
    
    N_best = N_results/200;
    
    %% Make plots
%     Fit_Medlyn_USO(A_n_subset, c_a_subset, g_w_subset, VPD_L_subset, T_L_subset, psi_L_MD_subset, V_cmax_subset, R_d_subset, Gamma_star_subset, K_m_subset, species_subset, Trt_numb_subset, Comb_Trt_name_subset, N_species_unique);
%     Plot_a_f_psi_relationship(a_f_max_results, pi_L_star_results, beta_results, epsilon_L_max_results, SWC_L_results, alpha_results, pi_L_0_25C_Trt_combo_results, psi_L_MD_subset, T_L_subset, SSE_best_total_subset_results, N_best, species_subset, Trt_numb_subset, Comb_Trt_name_subset, xmax, N_species_unique);
%     Plot_hist_and_covar(a_f_max_results, pi_L_star_results, beta_results, epsilon_L_max_results, SWC_L_results, alpha_results, pi_L_0_25C_Trt_combo_results, SSE_best_total_subset_results, N_best, xmin, xmax, Nx, species_subset);
%     Plot_g_w_layered(a_f_max_results, pi_L_star_results, beta_results, epsilon_L_max_results, SWC_L_results, alpha_results, pi_L_0_25C_Trt_combo_results, SSE_best_total_subset_results, N_best, species_subset, Trt_numb_subset, Comb_Trt_name_subset, g_w_subset, N_species_unique);
%     Plot_g_w_panel(a_f_max_results, pi_L_star_results, beta_results, epsilon_L_max_results, SWC_L_results, alpha_results, pi_L_0_25C_Trt_combo_results, SSE_best_total_subset_results, N_best, species_subset, Trt_numb_subset, Comb_Trt_name_subset, g_w_subset, N_species_unique);
%     Plot_g_w(a_f_max_results, pi_L_star_results, beta_results, epsilon_L_max_results, SWC_L_results, alpha_results, pi_L_0_25C_Trt_combo_results, SSE_best_total_subset_results, N_best, species_subset, Trt_numb_subset, Comb_Trt_name_subset, g_w_subset);
%     Plot_g_w_environment_response(a_f_max_results, pi_L_star_results, beta_results, epsilon_L_max_results, SWC_L_results, alpha_results, pi_L_0_25C_Trt_combo_results, SSE_best_total_subset_results, N_best, species_subset, T_L_subset, V_cmax_subset);
%     Plot_g_w_environment_response_and_data(pi_L_0_25C_Trt_combo_results, SSE_best_total_subset_results, a_f_max_results, pi_L_star_results, beta_results, epsilon_L_max_results, SWC_L_results, alpha_results, N_best, species_subset, Trt_numb_subset, Comb_Trt_name_subset, g_w_subset, psi_L_MD_subset, T_L_subset, VPD_L_subset, P_atm_subset, c_a_subset, V_cmax_subset);
%     Plot_violin_pi_L_0(a_f_max_results, pi_L_star_results, beta_results, epsilon_L_max_results, SWC_L_results, alpha_results, pi_L_0_25C_Trt_combo_results, SSE_best_total_subset_results, N_species_unique, N_best, species_subset, Trt_numb_subset, Comb_Trt_name_subset);
%     Plot_g_w_VPD_response(a_f_max_results, pi_L_star_results, beta_results, epsilon_L_max_results, SWC_L_results, alpha_results, pi_L_0_25C_Trt_combo_results, SSE_best_total_subset_results, N_best, species_subset, T_L_subset, psi_L_MD_subset, V_cmax_subset); 
%     plot_sgn_d2Xdotdgw2(N_best, Trt_numb_subset, SSE_best_total_subset_results, pi_L_0_25C_Trt_combo_results, a_f_max_results, pi_L_star_results, beta_results, epsilon_L_max_results, T_L_subset, psi_L_MD_subset);
    if i == 4
        % Ponderosa pine
%         Plot_g_w_VPD_response_combined(a_f_max_results, pi_L_star_results, beta_results, epsilon_L_max_results, SWC_L_results, alpha_results, pi_L_0_25C_Trt_combo_results, SSE_best_total_subset_results, N_best, species_subset, T_L_subset, psi_L_MD_subset, V_cmax_subset);
%         Plot_T_a_crit(a_f_max_results, pi_L_star_results, beta_results, epsilon_L_max_results, SWC_L_results, alpha_results, pi_L_0_25C_Trt_combo_results, SSE_best_total_subset_results, N_best, species_subset, T_L_subset, psi_L_MD_subset, V_cmax_subset);
    end
    if i == 2
        % Norway Spruce
        Run_and_Plot_Dynamic_Hydro(pi_L_0_25C_Trt_combo_results, SSE_best_total_subset_results, N_best, species_subset, a_f_max_results, pi_L_star_results, beta_results, epsilon_L_max_results, SWC_L_results, alpha_results, T_L_subset, V_cmax_subset, psi_L_MD_subset);
    end
    
    disp(['Completed plotting DREAM KZS results for species #', num2str(i), ' of ' num2str(N_species_unique), ': ', species_subset])
    
end
