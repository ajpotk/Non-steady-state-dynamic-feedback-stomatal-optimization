function [] = Plot_g_w_environment_response_and_data(pi_L_0_25C_Trt_combo_results, SSE_best_total_subset_results, ...
                                                     a_f_max_results, pi_L_star_results, beta_results, epsilon_L_max_results, SWC_L_results, alpha_results, ...
                                                     N_best, species_subset, Trt_numb_subset, Comb_Trt_name_subset, g_w_subset, psi_L_MD_subset, T_L_subset, ...
                                                     VPD_L_subset, P_atm_subset, c_a_subset, V_cmax_subset)

% Number of lines per plot
N_lines_psi_L_MD = 2;
N_lines_T_L = 3;
N_lines_pi_L_0_25C = 3;
P_min = 0.10;
P_max = 0.90;

% colors & markers
colors = [102, 102, 255;...
          124, 203, 161; ...
          240, 116, 110; ...
          253, 222, 156; ...
          002, 144, 153; ...
          220, 057, 120; ...
          006, 082, 117]/255;
markers = {'o'; 's'; 'd'; '^'; 'v'; '+'; 'p'};
colors_1 = colors(1,:); %user-specified
colors_2 = colors(6,:); %user-specified
colors_3 = colors(4,:); %user-specified
N_colors1 = 50; %user-specified
N_colors2 = 50; %user-specified
colors_exp1 = [1, 1, 1]; %user-specified
colors_exp2 = [1, 1, 1]; %user-specified
colors_map1 = zeros(N_colors1, 3);
colors_map2 = zeros(N_colors2, 3);
for i = 1:3
    colors_map1(:,i) = colors_1(i) + (colors_2(i) - colors_1(i)) * linspace(0,1,N_colors1)'.^colors_exp1(i);
    colors_map2(:,i) = colors_2(i) + (colors_3(i) - colors_2(i)) * linspace(0,1,N_colors2)'.^colors_exp2(i);
end
N_colors = N_colors1 + N_colors2;
colors_map = [colors_map1; colors_map2];

color_grey = 0.70*[1, 1, 1];

% % % figure
% % % hold on
% % % for i = 1:N_colors
% % %     area([0, 1], (1-(i-1)/N_colors)*[1, 1], 'FaceColor', colors_map(i,:), 'EdgeColor', 'none')
% % % end
% % % hold off
% % % ylim([0, 1])
% % % xlim([0, 1])

N_Trt_subset_combo = max(Trt_numb_subset);
Trt_subset_combo = cell(1, N_Trt_subset_combo);
for i = 1:N_Trt_subset_combo
    Trt_subset_combo_local = unique(Comb_Trt_name_subset(Trt_numb_subset == i));
    if ~isempty(Trt_subset_combo_local)
        Trt_subset_combo(i) = Trt_subset_combo_local;
    end
end

% specify names of species and treatments for plots
[~, Trt_subset_combo_plot] = Specify_species_and_Trt_names(species_subset, Trt_subset_combo);

% add tree-specific labels if necesary 
[~, ind_order] = sort(Trt_subset_combo_plot); %ordered alphabetically by legend entry
checked = zeros(1, N_Trt_subset_combo);
for i = ind_order
    if ~checked(i)
        is_Trt_subset_combo_plot_same = cell2mat(cellfun(@(x) strcmp(Trt_subset_combo_plot{i}, x), Trt_subset_combo_plot, 'UniformOutput', 0));
        if sum(is_Trt_subset_combo_plot_same) > 1
            ind_change = find(is_Trt_subset_combo_plot_same);
            k = 0;
            for j = ind_change
                k = k + 1;
                Trt_subset_combo_plot{j} = [Trt_subset_combo_plot{j}, ' - ', num2str(k)];
            end
        end
        checked(is_Trt_subset_combo_plot_same == 1) = 1;
    end
end

%% Regress V_cmax25 as function of psi_L_MD
V_cmax25_func = Regress_Vcmax25_as_func_psi(T_L_subset, V_cmax_subset, psi_L_MD_subset);

%% user-specified number of VPD bins with equal ranges of VPD
dVPD_L_subset = 1;
VPD_L_subset_max = max(VPD_L_subset);
VPD_L_subset_min = min(VPD_L_subset);
VPD_L_subset_max = dVPD_L_subset * ceil(VPD_L_subset_max/dVPD_L_subset);
VPD_L_subset_min = dVPD_L_subset * floor(VPD_L_subset_min/dVPD_L_subset);
N_VPD_L_bins = (VPD_L_subset_max - VPD_L_subset_min)/dVPD_L_subset;
VPD_L_bins = linspace(VPD_L_subset_min, VPD_L_subset_max, N_VPD_L_bins+1);
VPD_L_bins_min = VPD_L_bins(1:N_VPD_L_bins);
VPD_L_bins_max = VPD_L_bins(2:N_VPD_L_bins+1);
N_VPD_L_bins = 2; %only consider the first two bins

% % % % two VPD bins: high & low with equal amount of data points in each
% % % N_VPD_L_bins = 2;
% % % VPD_L_subset_max = max(VPD_L_subset);
% % % VPD_L_subset_min = min(VPD_L_subset);
% % % VPD_L_subset_med = median(VPD_L_subset);
% % % VPD_L_bins = [VPD_L_subset_min, VPD_L_subset_med, VPD_L_subset_max];
% % % VPD_L_bins_min = VPD_L_bins(1:N_VPD_L_bins);
% % % VPD_L_bins_max = VPD_L_bins(2:N_VPD_L_bins+1);

switch N_VPD_L_bins
    case 1
        VPD_L_bins_names = {''};
    case 2
        VPD_L_bins_names = {'Low VPD', 'High VPD'};
    case 3
        VPD_L_bins_names = {'Low VPD', 'Medium VDP', 'High VPD'};
    otherwise
        error('ERROR: Not enough cases consdered to name bins for VPD!!!')
end


ind_tile_psi_L_MD = 1:N_VPD_L_bins;
ind_tile_T_L = N_VPD_L_bins + (1:N_VPD_L_bins);
ind_tile_pi_L_0_25C = 2*N_VPD_L_bins + (1:N_VPD_L_bins);

g_w_subset_max = max(g_w_subset(VPD_L_subset <= VPD_L_bins_max(N_VPD_L_bins)));
dg_w_subset_max = 0.1;
g_w_subset_max = dg_w_subset_max * ceil(g_w_subset_max/dg_w_subset_max);

neg_psi_L_MD_subset_max_true = max(-psi_L_MD_subset(VPD_L_subset <= VPD_L_bins_max(N_VPD_L_bins)));
neg_psi_L_MD_subset_min_true = min(-psi_L_MD_subset(VPD_L_subset <= VPD_L_bins_max(N_VPD_L_bins)));
dneg_psi_L_MD_subset_max = 1;
dneg_psi_L_MD_subset_min = dneg_psi_L_MD_subset_max;
neg_psi_L_MD_subset_max = dneg_psi_L_MD_subset_max * ceil(neg_psi_L_MD_subset_max_true/dneg_psi_L_MD_subset_max);
neg_psi_L_MD_subset_min = dneg_psi_L_MD_subset_min * floor(neg_psi_L_MD_subset_min_true/dneg_psi_L_MD_subset_min);

T_L_subset_max_true = max(T_L_subset(VPD_L_subset <= VPD_L_bins_max(N_VPD_L_bins)));
T_L_subset_min_true = min(T_L_subset(VPD_L_subset <= VPD_L_bins_max(N_VPD_L_bins)));
dT_L_subset_max = 2;
dT_L_subset_min = dT_L_subset_max;
T_L_subset_max = dT_L_subset_max * ceil(T_L_subset_max_true/dT_L_subset_max);
T_L_subset_min = dT_L_subset_min * floor(T_L_subset_min_true/dT_L_subset_min);

%% Plot theoretical predictions
f_right_TL = 0.2; %fraction of figure on the left that the tiledlayout does not occupy
set_ind = 0;
set(0,'units','pixels')
Fig_size = [240, 40, 470, 500];
figure
set(gcf, 'Position', Fig_size)
tiledlayout(3, N_VPD_L_bins, 'TileSpacing', 'tight', 'OuterPosition', [f_right_TL, 0, 1-f_right_TL, 1]);

% best parameters for theoretical predictions
[~, ind_best_sub, ~, ~, ~, ~, ~, ~, ~, a_f_max, pi_L_star, beta, epsilon_L_max, SWC_L, alpha] = Choose_results(pi_L_0_25C_Trt_combo_results, SSE_best_total_subset_results, N_best, species_subset, a_f_max_results, pi_L_star_results, beta_results, epsilon_L_max_results, SWC_L_results, alpha_results); 
ind_best_sub = ind_best_sub(1);
pi_L_0_25C = mean(pi_L_0_25C_Trt_combo_results(ind_best_sub,:), 'omitnan'); 

% add theoretical lines to data
N_X = 2e2;
P_lines_psi_L_MD = P_min + (P_max - P_min) * linspace(0, 1, N_lines_psi_L_MD)';
P_lines_T_L = P_min + (P_max - P_min) * linspace(0, 1, N_lines_T_L)';
for i = 1:N_VPD_L_bins
    
    % mean VPD_L for bin
    VPD_L_bins_mean = (VPD_L_bins_min(i) + VPD_L_bins_max(i))/2;
    ind_bin = find((VPD_L_subset >= VPD_L_bins_min(i)) + (VPD_L_subset <= VPD_L_bins_max(i)) == 2);
    X_VPD_L = VPD_L_bins_mean * ones(1, N_X);
    
    % other bin-specific terms
    P_atm_bins = P_atm_subset(ind_bin);
    X_P_atm = median(P_atm_bins) * ones(1, N_X);
    c_a_bins = c_a_subset(ind_bin);
    X_c_a = median(c_a_bins) * ones(1, N_X);
    T_L_bins = T_L_subset(ind_bin);
    
    %leaf water potential as x-axis
    X_psi_L_MD = linspace(-neg_psi_L_MD_subset_max, 0, N_X);
    
    T_L_lines = quantile(T_L_bins, P_lines_psi_L_MD);
    norm_T_L_lines = (T_L_lines - T_L_subset_min_true)/(T_L_subset_max_true - T_L_subset_min_true);
    ind_colors_lines = ceil(N_colors * norm_T_L_lines);
    colors_lines = colors_map(ind_colors_lines, :);
    nexttile(ind_tile_psi_L_MD(i))
    for j = 1:N_lines_psi_L_MD
        
        % one line for each T_L
        X_T_L = T_L_lines(j) * ones(1, N_X);
        [outputs_photo_param] = Photosynthesis_parameters_temperature_response(X_T_L, X_psi_L_MD, V_cmax25_func);
        X_V_cmax = outputs_photo_param.V_cmax_vect;
        X_R_d = outputs_photo_param.R_d_vect;
        X_Gamma_star = outputs_photo_param.Gamma_star_vect;
        X_K_m = outputs_photo_param.K_m_vect;
        
        output = Stomata_and_PV_for_Trt(a_f_max, pi_L_star, beta, epsilon_L_max, SWC_L, alpha, ...
                                        pi_L_0_25C, X_psi_L_MD, X_VPD_L, ...
                                        X_P_atm, X_T_L, X_V_cmax, ...
                                        X_R_d, X_Gamma_star, X_K_m, ...
                                        X_c_a);
        X_g_w = output.g_w;
        X_g_w(output.solved == 0) = nan;
        X_g_w(X_g_w < 0) = 0;
        hold on
        plot(-X_psi_L_MD, X_g_w, '-', 'LineWidth', 2, 'Color', colors_lines(j,:));
        hold off
    end
    
    %leaf temperature as x-axis
    X_T_L = linspace(T_L_subset_min, T_L_subset_max, N_X);
    
    neg_psi_L_MD_lines = quantile(-psi_L_MD_subset(ind_bin), P_lines_T_L);
    norm_neg_psi_L_MD_lines = (neg_psi_L_MD_lines - neg_psi_L_MD_subset_min_true)/(neg_psi_L_MD_subset_max_true - neg_psi_L_MD_subset_min_true);
    ind_colors_lines = ceil(N_colors * norm_neg_psi_L_MD_lines);
    colors_lines = colors_map(ind_colors_lines, :);
    nexttile(ind_tile_T_L(i))
    for j = 1:N_lines_T_L
        
        % one line for each psi_L_MD
        X_psi_L_MD = -neg_psi_L_MD_lines(j) * ones(1, N_X);
        [outputs_photo_param] = Photosynthesis_parameters_temperature_response(X_T_L, X_psi_L_MD, V_cmax25_func);
        X_V_cmax = outputs_photo_param.V_cmax_vect;
        X_R_d = outputs_photo_param.R_d_vect;
        X_Gamma_star = outputs_photo_param.Gamma_star_vect;
        X_K_m = outputs_photo_param.K_m_vect;
        
        output = Stomata_and_PV_for_Trt(a_f_max, pi_L_star, beta, epsilon_L_max, SWC_L, alpha, ...
                                        pi_L_0_25C, X_psi_L_MD, X_VPD_L, ...
                                        X_P_atm, X_T_L, X_V_cmax, ...
                                        X_R_d, X_Gamma_star, X_K_m, ...
                                        X_c_a);
        X_g_w = output.g_w;
        X_g_w(output.solved == 0) = nan;
        X_g_w(X_g_w < 0) = 0;
        
        % additional line for each psi_L_MD without temperature response of chi_w
        X_chi_w = output.chi_w;
        X_chi_w(output.solved == 0) = nan;
        X_chi_w(X_chi_w < 0) = 0;
        X_chi_w_mean = mean(X_chi_w); 
        
        output = Stomata_and_PV_for_Trt(a_f_max, pi_L_star, beta, epsilon_L_max, SWC_L, alpha, ...
                                        pi_L_0_25C, X_psi_L_MD, X_VPD_L, ...
                                        X_P_atm, X_T_L, X_V_cmax, ...
                                        X_R_d, X_Gamma_star, X_K_m, ...
                                        X_c_a, ...
                                        'Constant_chi_w', X_chi_w_mean);
        X_g_w_Constant_chi_w = output.g_w;
        X_g_w_Constant_chi_w(output.solved == 0) = nan;
        X_g_w_Constant_chi_w(X_g_w_Constant_chi_w < 0) = 0;
        
        hold on
        plot(X_T_L, X_g_w, '-', 'LineWidth', 2, 'Color', colors_lines(j,:));
        plot(X_T_L, X_g_w_Constant_chi_w, ':', 'LineWidth', 1.2, 'Color', color_grey);
        hold off
        
    end
    
end

rel_alphab_y = 0.95;
rel_alphab_x = 0.06;
for i = 1:N_VPD_L_bins
    
    nexttile(ind_tile_psi_L_MD(i))
    ylim([0, g_w_subset_max])
    xlim([0, neg_psi_L_MD_subset_max]);
    caxis([T_L_subset_min_true, T_L_subset_max_true])
    set(gca, 'ytick', 0:dg_w_subset_max:g_w_subset_max)
    %%% title('Leaf water potential', 'FontSize', 10)
    if i == 1
        ylabel('{\itg_w} [mol\cdotm^{-2}\cdots^{-1}]')
        pos_ax = get(gca, 'OuterPosition');
        pos_annot = [0, pos_ax(2), pos_ax(1) + f_right_TL, pos_ax(4)];
        annotation('textbox', pos_annot, 'string', 'Leaf water potential', 'FontSize', 10, 'fontweight', 'bold', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'LineStyle', 'none')
    end
    xlabel('-{\it\psi_L} [MPa]')
    TITLE = VPD_L_bins_names{i};
    SUBTITLE = [num2str(VPD_L_bins_min(i)), ' \leq {\itD_L} \leq ', num2str(VPD_L_bins_max(i)) ' kPa'];
    title(TITLE, 'FontSize', 10)
    subtitle(SUBTITLE, 'FontSize', 8)
    LETTER = ['(', char(96 + i), ')'];
    text(0 + rel_alphab_x*(neg_psi_L_MD_subset_max - 0), 0 + rel_alphab_y*(g_w_subset_max - 0), ...
         LETTER, 'fontsize', 10, 'fontweight', 'bold', 'BackgroundColor', 'w')
    
    
    nexttile(ind_tile_T_L(i))
    ylim([0, g_w_subset_max])
    xlim([T_L_subset_min, T_L_subset_max])
    caxis([neg_psi_L_MD_subset_min_true, neg_psi_L_MD_subset_max_true])
    set(gca, 'ytick', 0:dg_w_subset_max:g_w_subset_max)
    %%% title('Leaf temperature', 'FontSize', 10)
    if i == 1
        ylabel('{\itg_w} [mol\cdotm^{-2}\cdots^{-1}]')
        pos_ax = get(gca, 'OuterPosition');
        pos_annot = [0, pos_ax(2), pos_ax(1) + f_right_TL, pos_ax(4)];
        annotation('textbox', pos_annot, 'string', 'Leaf temperature', 'FontSize', 10, 'fontweight', 'bold', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'LineStyle', 'none')
    end
    xlabel(['{\itT_L} [', char(176), 'C]'])
    LETTER = ['(', char(96 + N_VPD_L_bins + i), ')'];
    text(T_L_subset_min + rel_alphab_x*(T_L_subset_max - T_L_subset_min), 0 + rel_alphab_y*(g_w_subset_max - 0), ...
         LETTER, 'fontsize', 10, 'fontweight', 'bold', 'BackgroundColor', 'w')
    
end

%% add data to plot
for i = ind_order
    ind_plot_local = (Trt_numb_subset == i);
    if sum(ind_plot_local) > 0
        
        set_ind = set_ind + 1;
        
        g_w_subset_local = g_w_subset(ind_plot_local == 1);
        psi_L_MD_subset_local = psi_L_MD_subset(ind_plot_local == 1);
        T_L_subset_local = T_L_subset(ind_plot_local == 1);
        VPD_L_subset_local = VPD_L_subset(ind_plot_local == 1);
        
        %leaf water potential as x-axis
        for j = 1:N_VPD_L_bins
            
            ind_bin = find((VPD_L_subset_local >= VPD_L_bins_min(j)) + (VPD_L_subset_local <= VPD_L_bins_max(j)) == 2);
            nexttile(ind_tile_psi_L_MD(j))
            hold on
            scatter(-psi_L_MD_subset_local(ind_bin), g_w_subset_local(ind_bin), 36, T_L_subset_local(ind_bin), markers{i}, 'filled')
            hold off
            
        end
        
        
        %leaf temperature as x-axis
        for j = 1:N_VPD_L_bins
            
            ind_bin = find((VPD_L_subset_local >= VPD_L_bins_min(j)) + (VPD_L_subset_local <= VPD_L_bins_max(j)) == 2);
            nexttile(ind_tile_T_L(j))
            hold on
            scatter(T_L_subset_local(ind_bin), g_w_subset_local(ind_bin), 36, -psi_L_MD_subset_local(ind_bin), markers{i}, 'filled')
            hold off
            
        end
    end
end

colormap(colors_map)

%% theoretical response to reference osmotic potential at full hydration
X_pi_L_0_25C = pi_L_0_25C * linspace(0.8, 1.2, N_X); %range based on plasticity of the turgor loss point
P_lines_pi_L_0_25C = P_min + (P_max - P_min) * linspace(0, 1, N_lines_pi_L_0_25C)';

for i = 1:N_VPD_L_bins
    
    nexttile(ind_tile_pi_L_0_25C(i))
    ind_bin = find((VPD_L_subset >= VPD_L_bins_min(i)) + (VPD_L_subset <= VPD_L_bins_max(i)) == 2);
    P_atm = median(P_atm_subset(ind_bin));
    c_a = median(c_a_subset(ind_bin));
    VPD_L = (VPD_L_bins_min(i) + VPD_L_bins_max(i))/2;
    T_L = median(T_L_subset(ind_bin));
    neg_psi_L_MD_lines = quantile(-psi_L_MD_subset(ind_bin), P_lines_pi_L_0_25C);
    norm_neg_psi_L_MD_lines = (neg_psi_L_MD_lines - neg_psi_L_MD_subset_min_true)/(neg_psi_L_MD_subset_max_true - neg_psi_L_MD_subset_min_true);
    ind_colors_lines = ceil(N_colors * norm_neg_psi_L_MD_lines);
    colors_lines = colors_map(ind_colors_lines, :);
    
    for j = 1:N_lines_pi_L_0_25C
        
        X_g_w = nan(1, N_X);
        psi_L_MD = -neg_psi_L_MD_lines(j);
        
        [outputs_photo_param] = Photosynthesis_parameters_temperature_response(T_L, psi_L_MD, V_cmax25_func);
        V_cmax = outputs_photo_param.V_cmax_vect;
        R_d = outputs_photo_param.R_d_vect;
        Gamma_star = outputs_photo_param.Gamma_star_vect;
        K_m = outputs_photo_param.K_m_vect;

        for k = 1:N_X
            output = Stomata_and_PV_for_Trt(a_f_max, pi_L_star, beta, epsilon_L_max, SWC_L, alpha, ...
                                            X_pi_L_0_25C(k), psi_L_MD, VPD_L, ...
                                            P_atm, T_L, V_cmax, ...
                                            R_d, Gamma_star, K_m, ...
                                            c_a);

            if output.solved == 1
                X_g_w(k) = max(output.g_w, 0);
            end
        end

        hold on
        plot(-X_pi_L_0_25C, X_g_w, '-', 'LineWidth', 2, 'Color', colors_lines(j,:))
        hold off
        
        neg_pi_L_0_25C_min_local = min(-X_pi_L_0_25C(X_g_w <= g_w_subset_max));
        neg_pi_L_0_25C_max_local = max(-X_pi_L_0_25C(X_g_w >= 0));
        dneg_pi_L_0_25C_min = 0.1;
        dneg_pi_L_0_25C_max = dneg_pi_L_0_25C_min;
        neg_pi_L_0_25C_min_local = dneg_pi_L_0_25C_min * floor(neg_pi_L_0_25C_min_local/dneg_pi_L_0_25C_min);
        neg_pi_L_0_25C_max_local = dneg_pi_L_0_25C_max * ceil(neg_pi_L_0_25C_max_local/dneg_pi_L_0_25C_max);
        if (i == 1) && (j == 1)
            neg_pi_L_0_25C_min = neg_pi_L_0_25C_min_local;
            neg_pi_L_0_25C_max = neg_pi_L_0_25C_max_local; 
        else
            neg_pi_L_0_25C_min = min(neg_pi_L_0_25C_min, neg_pi_L_0_25C_min_local);
            neg_pi_L_0_25C_max = max(neg_pi_L_0_25C_max, neg_pi_L_0_25C_max_local); 
        end
        
    end
end

if isempty(neg_pi_L_0_25C_min) && isempty(neg_pi_L_0_25C_max)
    neg_pi_L_0_25C_min = 0;
    neg_pi_L_0_25C_max = 1;
end

for i = 1:N_VPD_L_bins
    
    nexttile(ind_tile_pi_L_0_25C(i))
    ylim([0, g_w_subset_max])
    xlim([neg_pi_L_0_25C_min, neg_pi_L_0_25C_max])
    set(gca, 'ytick', 0:dg_w_subset_max:g_w_subset_max)
    if i == 1
        ylabel('{\itg_w} [mol\cdotm^{-2}\cdots^{-1}]')
        pos_ax = get(gca, 'OuterPosition');
        pos_annot = [0, pos_ax(2), pos_ax(1) + f_right_TL, pos_ax(4)];
        annotation('textbox', pos_annot, 'string', 'Osmotic potential at full hydration', 'FontSize', 10, 'fontweight', 'bold', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'LineStyle', 'none')
    end
    xtick = get(gca, 'xtick');
    N_xtick = length(xtick);
    while N_xtick > 4
        N_xtick = floor(N_xtick/2);
        xtick = linspace(xtick(1), xtick(end), N_xtick); 
    end
    set(gca, 'xtick', xtick);
    xlabel(['-{\it\pi_L}_{,0,25', char(176), 'C} [MPa]'])
    LETTER = ['(', char(96 + 2*N_VPD_L_bins + i), ')'];
    text(neg_pi_L_0_25C_min + rel_alphab_x*(neg_pi_L_0_25C_max - neg_pi_L_0_25C_min), 0 + rel_alphab_y*(g_w_subset_max - 0), ...
         LETTER, 'fontsize', 10, 'fontweight', 'bold', 'BackgroundColor', 'w')
    box off
    
end

nexttile(ind_tile_psi_L_MD(N_VPD_L_bins))
cb1 = colorbar;
ylabel(cb1, ['{\itT_L} [', char(176), 'C]'], 'Rotation', 270, 'VerticalAlignment', 'bottom', 'FontSize', 8)

nexttile(ind_tile_T_L(N_VPD_L_bins))
cb2 = colorbar; 
ylabel(cb2, '-{\it\psi_L} [MPa]', 'Rotation', 270, 'VerticalAlignment', 'bottom', 'FontSize', 8)

nexttile(ind_tile_pi_L_0_25C(N_VPD_L_bins))
cb3 = colorbar;
ylabel(cb3, '-{\it\psi_L} [MPa]', 'Rotation', 270, 'VerticalAlignment', 'bottom', 'FontSize', 8)
caxis([neg_psi_L_MD_subset_min_true, neg_psi_L_MD_subset_max_true])

print(gcf, '-djpeg', '-r300', [species_subset, '_Fig_environmental_response_w_data.jpg'], '-painters', '-noui' )

end

