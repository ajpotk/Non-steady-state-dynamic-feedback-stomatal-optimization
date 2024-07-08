function [] = Plot_T_a_crit(a_f_max_results, pi_L_star_results, beta_results, epsilon_L_max_results, ...
                            SWC_L_results, alpha_results, pi_L_0_25C_Trt_combo_results, SSE_best_total_subset_results, ...
                            N_best, species_subset, T_L_subset, psi_L_MD_subset, V_cmax_subset)

%% Colors

% colors
color_red = [240, 168, 168]/255;
color_dark_red = [167, 0, 0]/255;
color_blue = [170, 240, 255]/255;
color_dark_blue = [0, 68, 170]/255;
color_green = [161, 207, 171]/255;
color_purple = [102, 102, 255]/255;
color_dark_purple = [102, 0, 102]/255;
color_grey = [0.7, 0.7, 0.7];
color_gold = [0.9, 0.7, 0.1]; 
face_alpha = 0.3;
color_mod = [0.6, 0.6, 0.6];

%% Constants, parameters, and conditions
emiss_L = 0.97; %leaf emissivity [-]
g_H_a = 3; %boundary layer heat conductance [mol m-2 s-1]
R_abs = 300; %absorbed shortwave radiation [W m-2]
RH = 0.5;
P_atm = 101.325; %atmospheric pressure [kPa]
c_a = 4e-4; %atmospheric CO2 concentration [mol mol-1]
k_L_max_25C = 7.22e-3; % maximum soil-plant conductance per unit leaf area at 25C [mol m-2 s-1 MPa-1] -- value based on Pinon Pine
A = 1.4; %exponent for simple soil-plant conductance [MPa-1] -- k_L = k_L_max * exp(A*psi_L) -- value based on Pinon Pine
Q_10_k_L = 1.5; %Q10 value for the temperature-dependence of the maximum soil-plant conductance per unit leaf area [-]
Nonmonotonic_V_cmax_temperature_response = 1;

% estimate air temperature and soil water potential for subset
outputs_Regress_Sapes_VC = Regress_Sapes_VC('Plot', 0, 'Q10_regression', 1);
T_a_subset = outputs_Regress_Sapes_VC.T_L_raw;
psi_soil_subset = outputs_Regress_Sapes_VC.psi_soil_raw;
psi_L_subset = outputs_Regress_Sapes_VC.psi_L_raw;
psi_L_crit_subset = outputs_Regress_Sapes_VC.psi_L_crit_raw;
psi_L_crit_subset_conf_interv = outputs_Regress_Sapes_VC.psi_L_crit_raw_conf_interv;
SM_crit_subset = psi_L_subset - psi_L_crit_subset;
SM_crit_subset_conf_interv = psi_L_subset - psi_L_crit_subset_conf_interv;

%% Best parameters for theoretical predictions
N_best = N_best/10;
N_best(N_best < 100) = 100;

[~, ind_best_sub, ~, ~, ~, ~, ~, ~, ~, a_f_max, pi_L_star, beta, epsilon_L_max, SWC_L, alpha] = Choose_results(pi_L_0_25C_Trt_combo_results, SSE_best_total_subset_results, N_best, species_subset, a_f_max_results, pi_L_star_results, beta_results, epsilon_L_max_results, SWC_L_results, alpha_results, 'Option', 2); 
ind_best_sub = ind_best_sub(1);
pi_L_0_25C = mean(pi_L_0_25C_Trt_combo_results(ind_best_sub,:), 'omitnan'); 


%% Regress V_cmax25 as function of psi_L_MD
V_cmax25_func = Regress_Vcmax25_as_func_psi(T_L_subset, V_cmax_subset, psi_L_MD_subset); 


%% Plot steady-state response between g_w, T_L, and VPD, varying VPD by changing T_a (with leaf water potential calculated from steady-state hydraulics)

Mat_data_save_file = 'Saved_T_a_crit_data.mat';

N_test = 3e2;
psi_soil_min = -3; %[MPa]
psi_soil_max = 0; %[MPa]
psi_soil_test = linspace(psi_soil_max, psi_soil_min, N_test); 
dpsi_soil_test = (psi_soil_max - psi_soil_min)/(N_test-1);

N_x_axis = 3e2;
T_a_min = 0; %[C]
T_a_max = 60; %[C]
T_a_max_plot = 50; %[C]
T_a_varied = linspace(T_a_min, T_a_max, N_x_axis);
dT_a_varied = (T_a_max - T_a_min)/(N_x_axis-1);

N_exp = 4;
g_H_a_exp = g_H_a * [1, 0.5, 1, 1];
R_abs_exp = R_abs * [1, 1, 0.5, 1];
RH_exp = RH * [1, 1, 1, 0.5];
color_exp = {color_purple; ...
             color_gold; ...
             color_blue; ...
             color_red};
Display_Name_exp = {'Control'; ...
                    '-50% {\itg_{H,a}}'; ...
                    '-50% {\itR_{abs}}'; ...
                    '-50% RH'};

% load past data if any
files = struct2cell(dir);
files = files(1,:);
is_Mat_data_save_file = cell2mat(cellfun(@(x) strcmp(x, Mat_data_save_file), files, 'UniformOutput', 0));

if any(is_Mat_data_save_file)

    load(Mat_data_save_file, 'T_a_crit_min_exp', 'T_a_crit_max_exp', 'SM_crit_store_exp', 'g_w_store_exp', 'convergence_store_exp', 'fail_convergence_exp')
    h_start = find(cell2mat(cellfun(@(x) isempty(x), T_a_crit_min_exp, 'UniformOutput', 0)), 1, 'first');
    
else
    
    T_a_crit_min_exp = cell(1, N_exp);
    T_a_crit_max_exp = cell(1, N_exp);
    SM_crit_store_exp = cell(1, N_exp);
    g_w_store_exp = cell(1, N_exp);
    convergence_store_exp = cell(1, N_exp);
    fail_convergence_exp = cell(1, N_exp);
    h_start = 1;
    
end

for h = h_start:N_exp
    
    g_H_a = g_H_a_exp(h);
    R_abs = R_abs_exp(h);
    RH = RH_exp(h);
    
    T_a_crit_min = nan(1, N_test);
    T_a_crit_max = nan(1, N_test);
    SM_crit_store = nan(N_x_axis, N_test);
    g_w_store = nan(N_x_axis, N_test);
    convergence_store = nan(N_x_axis, N_test);

    RH = RH(1) * ones(1, N_x_axis);
    P_atm = P_atm(1) * ones(1, N_x_axis);
    c_a = c_a(1) * ones(1, N_x_axis);

    progress_disp_inc = 0.01; 
    progress_disp = progress_disp_inc;

    for i = 1:N_test

        psi_soil = psi_soil_test(i); %soil water potential [MPa] -- used for steady-state predictions

        % vectors to store calculated data
        g_w_with_varied_T_a = nan(1, N_x_axis);
        VPD_L_with_varied_T_a = nan(1, N_x_axis);
        T_L_with_varied_T_a = nan(1, N_x_axis);
        E_with_varied_T_a = nan(1, N_x_axis);
        E_crit_with_varied_T_a = nan(1, N_x_axis);
        psi_L_with_varied_T_a = nan(1, N_x_axis);
        psi_L_crit_with_varied_T_a = nan(1, N_x_axis);
        convergence_with_varied_T_a = nan(1, N_x_axis);

        for j = 1:N_x_axis

            % display progress
            progress = (h - 1 + (i - 1 + j/N_x_axis)/N_test)/N_exp;
            if progress >= progress_disp
                disp(['        Calculating steady-state g_w for varied T_a: ', num2str(100*progress_disp_inc*floor(progress/progress_disp_inc), '%.0f'), '% done!'])
                progress_disp = progress_disp_inc * (1 + floor(progress/progress_disp_inc));
            end

            [outputs_stomata_steady_state] = ...
            Stomata_for_steady_state_thermo_and_hydraulics(SWC_L, alpha, a_f_max, pi_L_star, beta, epsilon_L_max, pi_L_0_25C, ...
                                                           k_L_max_25C, A, Q_10_k_L, psi_soil, ...
                                                           T_a_varied(j), RH(j), P_atm(j), c_a(j), ...
                                                           R_abs, g_H_a, emiss_L, ...
                                                           V_cmax25_func, ...
                                                           'Nonmonotonic_V_cmax_temperature_response', Nonmonotonic_V_cmax_temperature_response);

            g_w_with_varied_T_a(j) = outputs_stomata_steady_state.g_w_steady_state;
            VPD_L_with_varied_T_a(j) = outputs_stomata_steady_state.VPD_L_steady_state;
            T_L_with_varied_T_a(j) = outputs_stomata_steady_state.T_L_steady_state;
            E_with_varied_T_a(j) = outputs_stomata_steady_state.E_steady_state;
            E_crit_with_varied_T_a(j) = outputs_stomata_steady_state.E_crit_steady_state;
            psi_L_with_varied_T_a(j) = outputs_stomata_steady_state.psi_L_steady_state;
            psi_L_crit_with_varied_T_a(j) = outputs_stomata_steady_state.psi_L_crit_steady_state;
            convergence_with_varied_T_a(j) = outputs_stomata_steady_state.convergence;
            
            if  (psi_L_with_varied_T_a(j) < psi_L_crit_with_varied_T_a(j)) && (T_a_varied(j) < 40)
               error('!') 
            end
            
        end

        % determine critical air temperature
        SM_crit = psi_L_with_varied_T_a - psi_L_crit_with_varied_T_a; %safety margin from critical water potential [MPa]
        diff_sign_SM_crit = diff(sign(SM_crit));
        N_intersect = sum(diff_sign_SM_crit ~= 0);
        if N_intersect > 0

            N_intersect = min(N_intersect, 2);
            ind_LB = find((diff_sign_SM_crit ~= 0), N_intersect, 'first');
            ind_min_LB = ind_LB(1);
            ind_min_UB = ind_min_LB + 1;
            T_a_crit_min_LB = T_a_varied(ind_min_LB);
            T_a_crit_min_UB = T_a_varied(ind_min_UB);
            SM_crit_min_LB = SM_crit(ind_min_LB);
            SM_crit_min_UB = SM_crit(ind_min_UB);
            T_a_crit_min(i) = T_a_crit_min_LB + (T_a_crit_min_UB - T_a_crit_min_LB)*(0 - SM_crit_min_LB)/(SM_crit_min_UB - SM_crit_min_LB);

            if N_intersect == 2
                ind_max_LB = ind_LB(2);
                ind_max_UB = ind_max_LB + 1;  
                T_a_crit_max_LB = T_a_varied(ind_max_LB);
                T_a_crit_max_UB = T_a_varied(ind_max_UB);
                SM_crit_max_LB = SM_crit(ind_max_LB);
                SM_crit_max_UB = SM_crit(ind_max_UB);
                T_a_crit_max(i) = T_a_crit_max_LB + (T_a_crit_max_UB - T_a_crit_max_LB)*(0 - SM_crit_max_LB)/(SM_crit_max_UB - SM_crit_max_LB);
            elseif N_intersect == 1
                T_a_crit_max(i) = T_a_max;
            end

        end

        SM_crit_store(:,i) = SM_crit';
        g_w_store(:,i) = g_w_with_varied_T_a';
        convergence_store(:,i) = convergence_with_varied_T_a';
        
    end
    
    % solutions without full convergence are supercritical
    SM_crit_store(convergence_store < 4) = -999; %arbitarily negative value
    fail_convergence = zeros(size(convergence_store)); %binary for convergence failure
    fail_convergence(convergence_store < 4) = 1;
    
    % store values
    T_a_crit_min_exp{h} = T_a_crit_min;
    T_a_crit_max_exp{h} = T_a_crit_max;
    SM_crit_store_exp{h} = SM_crit_store;
    g_w_store_exp{h} = g_w_store;
    convergence_store_exp{h} = convergence_store;
    fail_convergence_exp{h} = fail_convergence;
    
    % save values
	save(Mat_data_save_file, 'T_a_crit_min_exp', 'T_a_crit_max_exp', 'SM_crit_store_exp', 'g_w_store_exp', 'convergence_store_exp', 'fail_convergence_exp')

end

%% Plots
figure
tiledlayout(2, 2)
set(gcf, 'position', [1, -200, 1100, 700])
plot_ind = 0;
plot_ind_leg = nan(1, N_exp);

for h = 1:N_exp
    
    % load values
    T_a_crit_min = T_a_crit_min_exp{h};
    T_a_crit_max = T_a_crit_max_exp{h};
    SM_crit_store = SM_crit_store_exp{h};
    g_w_store = g_w_store_exp{h};
    convergence_store = convergence_store_exp{h};
    fail_convergence = fail_convergence_exp{h};
    
    % contour plot of safety margin from critical
    if h == 1
        
        N_color = 50;
        color_map = [(linspace(1, color_purple(1), N_color)').^0.75, ...
                     (linspace(1, color_purple(2), N_color)').^0.75, ...
                     (linspace(1, color_purple(3), N_color)').^0.75];
        SM_min = -0.4;
        SM_crit_store(SM_crit_store < SM_min) = SM_min;
        SM_crit_store_neg = SM_crit_store; %only negative values
        SM_crit_store_neg(SM_crit_store_neg >= 0) = nan;
        SM_crit_store_neg(SM_crit_store_neg < SM_min) = SM_min;
        nexttile(1)
        box on
        set(gca,'Layer','top')
        hold on
        imagesc(-psi_soil_test, T_a_varied, -SM_crit_store_neg);
        contourf(-psi_soil_test, T_a_varied, -SM_crit_store_neg, 'LineColor', 'none');
        contour(-psi_soil_test, T_a_varied, SM_crit_store, SM_min:0.1:0, 'k--', 'linewidth', 1);
        contour(-psi_soil_test, T_a_varied, SM_crit_store, [0, 0], '--', 'EdgeColor', color_purple, 'linewidth', 4);
        hold off
        set(gca, 'ydir', 'normal')
        xlim([-psi_soil_max, -psi_soil_min])
        ylim([T_a_min, T_a_max_plot])
        ylabel(['{\itT_a} [', char(176), 'C]'])
        xlabel('-{\it\psi_{soil}} [MPa]')
        colormap(color_map)
        caxis([0, -SM_min])
        color_bar = colorbar;
        set(color_bar, 'Ticks', 0:0.1:-SM_min);
        c_ticklabels = get(color_bar, 'TickLabels');
        c_ticklabels{1} = ['\leq ', c_ticklabels{1}];
        c_ticklabels{end} = ['\geq ', c_ticklabels{end}];
        set(color_bar, 'TickLabels', c_ticklabels);
        c_label = ylabel(color_bar, '{\it\psi_{L,crit}} - {\it\psi_L} [MPa]', 'Rotation', 270, 'FontSize', 11);
        c_label.Position(1) = c_label.Position(1) + 1.3;


        % contour plot of stomatal conductance
        relative_g_w_store = 1e2*g_w_store/max(g_w_store, [], 'all'); %[%]
        nexttile(2)
        box on
        set(gca,'Layer','top')
        hold on
        contourf(-psi_soil_test, T_a_varied, relative_g_w_store, 'LineColor', 'none');
        contour(-psi_soil_test, T_a_varied, relative_g_w_store, 0:10:100, 'k--', 'linewidth', 1);
        contour(-psi_soil_test, T_a_varied, SM_crit_store, [0, 0], '--', 'EdgeColor', color_purple, 'linewidth', 4);
        xlim([-psi_soil_max, -psi_soil_min])
        ylim([T_a_min, T_a_max_plot])
        ylabel(['{\itT_a} [', char(176), 'C]'])
        xlabel('-{\it\psi_{soil}} [MPa]')
        colormap(color_map)
        caxis([0, 100])
        color_bar = colorbar;
        set(color_bar, 'Ticks', 0:20:100);
        c_ticklabels = get(color_bar, 'TickLabels');
        c_ticklabels = cellfun(@(x) [x, '%'], c_ticklabels, 'UniformOutput', 0);
        set(color_bar, 'TickLabels', c_ticklabels);
        c_label = ylabel(color_bar, 'Relative {\itg_w} [%]', 'Rotation', 270, 'FontSize', 11);
        c_label.Position(1) = c_label.Position(1) + 1.3;
    
    end

    % determine portions to plot
    not_nan = ~isnan(T_a_crit_min);
    ind = 1:N_test;
    N_not_nan = 0;
    ind_start_not_nan = nan(1, N_test);
    ind_end_not_nan = nan(1, N_test);
    accept_new_not_nan = 1;
    for i = 1:N_test
        if not_nan(i) && accept_new_not_nan
            N_not_nan = N_not_nan + 1;
            ind_start_not_nan(N_not_nan) = i;
            if any(~not_nan(i:end))
                ind_end_not_nan(N_not_nan) = find(((~not_nan) + (ind > i) == 2), 1, 'first') - 1;
            else
                ind_end_not_nan(N_not_nan) = N_test;
            end
            accept_new_not_nan = 0;
        elseif ~not_nan(i)
            accept_new_not_nan = 1;
        end
    end
    ind_start_not_nan = ind_start_not_nan(1:N_not_nan);
    ind_end_not_nan = ind_end_not_nan(1:N_not_nan);

    % plot supercritical temperature regions
    nexttile(3)
    box on
    set(gca,'Layer','top')
    hold on
    for i = 1:N_x_axis
        for j = 1:N_test
            if fail_convergence(i,j)
                fill(-psi_soil_test(j) + 0.5*dpsi_soil_test*[-1, 1, 1, -1], T_a_varied(i) + 0.5*dT_a_varied*[1, 1, -1, -1], ...
                     color_exp{h}, 'LineStyle', 'none', 'DisplayName', Display_Name_exp{h}, 'FaceAlpha', face_alpha);
            end
        end
    end
    for i = 1:N_not_nan
        if ind_start_not_nan(i) == ind_end_not_nan(i)
            if T_a_crit_min(ind_start_not_nan(i)) == T_a_crit_max(ind_start_not_nan(i))
                plot(-psi_soil_test(ind_start_not_nan(i)), T_a_crit_min(ind_start_not_nan(i)), 'o', 'LineStyle', 'none', 'MarkerFaceColor', [color_exp{h}, face_alpha], 'Color', [color_exp{h}, face_alpha]);
            else
                plot(-psi_soil_test(ind_start_not_nan(i))*[1, 1], [T_a_crit_min(ind_start_not_nan(i)), T_a_crit_max(ind_start_not_nan(i))], '-', 'LineWidth', 2, 'Color', [color_exp{h}, face_alpha]);
            end
        else
            plot_ind = plot_ind + 1;
            if isnan(plot_ind_leg(h))
                plot_ind_leg(h) = plot_ind;
            end
            Fill(plot_ind) = fill(-[psi_soil_test(ind_start_not_nan(i):ind_end_not_nan(i)), fliplr(psi_soil_test(ind_start_not_nan(i):ind_end_not_nan(i)))], ...
                                  [T_a_crit_max(ind_start_not_nan(i):ind_end_not_nan(i)), fliplr(T_a_crit_min(ind_start_not_nan(i):ind_end_not_nan(i)))], ...
                                  color_exp{h}, 'LineStyle', 'none', 'DisplayName', Display_Name_exp{h}, 'FaceAlpha', face_alpha);
            contour(-psi_soil_test, T_a_varied, SM_crit_store, [0, 0], '--', 'EdgeColor', color_exp{h}, 'linewidth', 4);
        end
    end
    hold off
    xlim([-psi_soil_max, -psi_soil_min])
    ylim([T_a_min, T_a_max_plot])
    ylabel(['{\itT_a} [', char(176), 'C]'])
    xlabel('-{\it\psi_{soil}} [MPa]')
    Leg = legend(Fill(plot_ind_leg(1:h)), 'Location', 'Layout');
    Leg.Layout.Tile = 'east';
    
    if h == N_exp
        text(mean(sum(-psi_soil_test.*fail_convergence)./sum(fail_convergence)), ...
             47, ... % % % mean(sum(T_a_varied'.*fail_convergence)./sum(fail_convergence)), ...
             'Unstable supercritical', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle')
        text(sum(-psi_soil_test.*(T_a_crit_max - T_a_crit_min), 'omitnan')/sum(T_a_crit_max - T_a_crit_min, 'omitnan'), ...
             sum(0.5*(T_a_crit_max.^2 - T_a_crit_min.^2), 'omitnan')/sum(T_a_crit_max - T_a_crit_min, 'omitnan'), ...
             ['Stable', newline, 'supercritical'], 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle')


    end
    
end

% add empirical data
hold on
Obs = plot(-psi_soil_subset, T_a_subset, 'o', 'MarkerSize', 4, 'LineStyle', 'none', 'Color', color_dark_blue, 'MarkerFaceColor', color_dark_blue, 'DisplayName', 'Observed');
hold off
Leg = legend([Fill(plot_ind_leg), Obs], 'Location', 'Layout');
Leg.Layout.Tile = 'east';

% plot of observed Safety margin from critical
conf_interv_tick_size = 0.05;
N_subset = length(SM_crit_subset);
is_subcrit = (SM_crit_subset_conf_interv(:,2) > 0); %positive safery margin is subcritical
is_supercrit = (SM_crit_subset_conf_interv(:,1) < 0); %negative safety margin is supercritival
is_maybe_crit = (~is_subcrit) .* (~is_supercrit);
nexttile(4)
box on
set(gca,'Layer','top')
hold on
plot([0, 3.5], [0, 0], '--', 'Color', color_purple, 'LineWidth', 2)
for i = 1:N_subset
    plot(-psi_soil_subset(i)*[1, 1], -SM_crit_subset_conf_interv(i,:), 'k', 'LineWidth', 0.5)
    plot(-psi_soil_subset(i) + conf_interv_tick_size*[-0.5, 0.5], -SM_crit_subset_conf_interv(i,1)*[1,1], 'k', 'LineWidth', 0.5)
    plot(-psi_soil_subset(i) + conf_interv_tick_size*[-0.5, 0.5], -SM_crit_subset_conf_interv(i,2)*[1,1], 'k', 'LineWidth', 0.5)
end
plot(-psi_soil_subset(is_subcrit == 1), -SM_crit_subset(is_subcrit == 1), 'o', 'LineStyle', 'none', 'MarkerFaceColor', color_dark_blue, 'Color', color_dark_blue)
plot(-psi_soil_subset(is_supercrit == 1), -SM_crit_subset(is_supercrit == 1), 'o', 'LineStyle', 'none', 'MarkerFaceColor', color_dark_red, 'Color', color_dark_red)
plot(-psi_soil_subset(is_maybe_crit == 1), -SM_crit_subset(is_maybe_crit == 1), 'o', 'LineStyle', 'none', 'MarkerFaceColor', color_purple, 'Color', color_purple)
hold off
ylabel('{\it\psi_{L,crit}} - {\it\psi_L} [MPa]')
xlabel('-{\it\psi_{soil}} [MPa]')
text(1, 0.05, ['{\color[rgb]{', num2str(color_dark_red) '}Supercritical}'], 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12)
text(1, -0.05, ['{\color[rgb]{', num2str(color_dark_blue) '}Subcritical}'], 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 12)
ylim([-1.5, 1.5])
xlim([0, 3.5])

% add letters & titles
Titles = {'Degree of supercriticalness'; ...
          'Stomatal conductance'; ...
          'Domain of supercriticalness'; ...
          'Observed degree of supercriticalness'};
rel_alphab_y = 1.02;
rel_alphab_x = 0.05;
for i = 1:4
    nexttile(i)
    xlim_min = min(xlim);
    xlim_max = max(xlim);
    ylim_min = min(ylim);
    ylim_max = max(ylim);
    LETTER = ['(', char(96 + i), ')'];
    text(xlim_min + rel_alphab_x*(xlim_max - xlim_min), ylim_min + rel_alphab_y*(ylim_max - ylim_min), ...
         [LETTER, ' ', Titles{i}], 'fontsize', 9, 'fontweight', 'bold', 'BackgroundColor', 'w', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left')
end

% save figure
Output_file = [species_subset, '_supercrit_domain.jpg'];
% print(gcf, '-djpeg', '-r300', Output_file, '-painters', '-noui' )
pause(0)

%% 2nd plot of supply vs demand
psi_soil_test = [0, -1, -1]; %[MPa]
T_a_test = [25, 25, 50]; %[C]
N_test = length(psi_soil_test);

RH = RH(1);
P_atm = P_atm(1);
c_a = c_a(1);

psi_L_crit_test = nan(1, N_test);

y_max = 0;
x_max = 0;

figure
set(gcf, 'position', [1, -200, 400, 700])
tiledlayout(N_test, 1, 'TileSpacing', 'none')
for i = 1:N_test
    
    psi_soil = psi_soil_test(i);
    T_a = T_a_test(i);
    
    [outputs_stomata_steady_state] = ...
    Stomata_for_steady_state_thermo_and_hydraulics(SWC_L, alpha, a_f_max, pi_L_star, beta, epsilon_L_max, pi_L_0_25C, ...
                                                   k_L_max_25C, A, Q_10_k_L, psi_soil, ...
                                                   T_a, RH, P_atm, c_a, ...
                                                   R_abs, g_H_a, emiss_L, ...
                                                   V_cmax25_func, ...
                                                   'Nonmonotonic_V_cmax_temperature_response', Nonmonotonic_V_cmax_temperature_response);
    
    E = outputs_stomata_steady_state.E_steady_state;
    psi_L = outputs_stomata_steady_state.psi_L_steady_state;
    psi_L_crit_test(i) = outputs_stomata_steady_state.psi_L_crit_steady_state;
    psi_L_vect = outputs_stomata_steady_state.psi_L_vect;
    E_supply_vect = outputs_stomata_steady_state.E_supply_vect;
    E_demand_vect = outputs_stomata_steady_state.E_demand_vect;
    E_crit = max(E_supply_vect); 
    VPD_L_vect = outputs_stomata_steady_state.VPD_L_vect;
    g_w_supply_vect = outputs_stomata_steady_state.g_w_supply_vect;
    convergence = outputs_stomata_steady_state.convergence;
    g_w_crit = g_w_supply_vect(E_supply_vect == E_crit);
    if i == 1
        g_w_max = min(g_w_crit, 0.03);
    end
    if convergence < 4
        E = nan;
        psi_L = nan;
        E_demand_vect = g_w_max*VPD_L_vect/P_atm;
    end
    
    nexttile(i)
    box on
    hold on
    P(1) = plot(-psi_L_vect, 1e3*E_supply_vect, '-', 'Color', color_dark_blue, 'LineWidth', 3, 'DisplayName', 'Supply');
    P(2) = plot(-psi_L_vect, 1e3*E_demand_vect, '-', 'Color', color_dark_red, 'LineWidth', 3, 'DisplayName', 'Demand');
    plot(-psi_L_crit_test(i)*[1, 1], 1e3*E_crit*[0, 1], 'k-', 'LineWidth', 0.9)
    plot(-psi_L_crit_test(i), 1e3*E_crit, 'ko', 'MarkerSize', 6, 'LineStyle', 'none', 'MarkerFaceColor', 'w')
    P(3) = plot(-psi_L, 1e3*E, 'd', 'MarkerSize', 10, 'LineStyle', 'none', 'MarkerFaceColor', color_dark_purple, 'MarkerEdgeColor', color_dark_purple, 'DisplayName', 'Solution');
    hold off
    
    if i == 1
       text(-psi_L_crit_test(i), 1e3*E_crit, '{\itE_{crit}}', ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'Color', 'k', 'FontSize', 9) 
       legend(P(1:3), 'Location', 'E')
       legend box on
    end
    
    y_max = max(y_max, max(ylim));
    x_max = max(x_max, max(xlim));
    
end

y_max = 3;

for i = 1:N_test
   
    nexttile(i)
    ylim([0, y_max])
    xlim([0, x_max])
    
    if i > 1
        ytick = get(gca, 'ytick');
        ytick(end) = [];
        set(gca, 'ytick', ytick)
    end
    
    if i < N_test
       set(gca, 'xticklabel', []); 
    end
    
    if i == 1
       text(-psi_L_crit_test(i) + 0.01*x_max, 0.01*y_max, '-{\it\psi_{L,crit}}', ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'Color', 'k', 'FontSize', 9) 
    end
    
    ylabel(['{\itE}', newline, '[mmol\cdotm^{-2}\cdots^{-1}]'])
    
end
xlabel('-{\it\psi_L} [MPa]')


% add letters & titles
Titles = {'Subcritical'; ...
          'Stable supercritical'; ...
          'Unstable supercritical'};

rel_alphab_y = 0.95;
rel_alphab_x = 0.05;
for i = 1:3
    nexttile(i)
    LETTER = ['(', char(96 + 4 + i), ')'];
    text(rel_alphab_x*x_max, rel_alphab_y*y_max, ...
         [LETTER, ' ', Titles{i}, newline, '({\it\psi_{soil}} = ', num2str(psi_soil_test(i)), ' MPa, {\itT_a} = ', num2str(T_a_test(i)), char(176), 'C)'], ...
         'fontsize', 9, 'fontweight', 'bold', 'BackgroundColor', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left')
end

% save figure
Output_file = [species_subset, '_supercrit_examples.jpg'];
% print(gcf, '-djpeg', '-r300', Output_file, '-painters', '-noui' )
pause(0)

end

