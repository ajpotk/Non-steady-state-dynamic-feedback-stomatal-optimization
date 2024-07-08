function [] = Plot_g_w_VPD_response_combined(a_f_max_results, pi_L_star_results, beta_results, epsilon_L_max_results, ...
                                             SWC_L_results, alpha_results, pi_L_0_25C_Trt_combo_results, SSE_best_total_subset_results, ...
                                             N_best, species_subset, T_L_subset, psi_L_MD_subset, V_cmax_subset)

%% Constants, parameters, and conditions
o_i = 1e-3 * 210; %standard atmospheric partial pressure of oxygen [mol mol-1]
R = 8.314; %universal gas constant [J mol-1 K-1]
R_d25_per_V_cmax25 = 0.015; %from Collatz et al. (1991)
emiss_L = 0.97; %leaf emissivity [-]
g_H_a = 3; %boundary layer heat conductance [mol m-2 s-1]
R_abs = 300; %absorbed shortwave radiation [W m-2]
RH = 0.5;
T_a = 24;
P_atm = 101.325; %atmospheric pressure [kPa]
c_a = 4e-4; %atmospheric CO2 concentration [mol mol-1]
RWC_t = 1; %leaf total relative water content [-] -- used for predictions at constant leaf temperature
psi_soil = 0; %soil water potential [MPa] -- used for steady-state predictions
k_L_max_25C = 7.22e-3; % maximum soil-plant conductance per unit leaf area at 25C [mol m-2 s-1 MPa-1] -- value based on Pinon Pine
A = 1.4; %exponent for simple soil-plant conductance [MPa-1] -- k_L = k_L_max * exp(A*psi_L) -- value based on Pinon Pine
Q_10_k_L = 1.5; %Q10 value for the temperature-dependence of the maximum soil-plant conductance per unit leaf area [-]
Nonmonotonic_V_cmax_temperature_response = 1;

%% Figure colors
colors = [102, 102, 255;...
          124, 203, 161; ...
          240, 116, 110; ...
          253, 222, 156; ...
          002, 144, 153; ...
          220, 057, 120; ...
          006, 082, 117]/255;

colors_1 = colors(1,:); %user-specified
colors_2 = colors(6,:); %user-specified
colors_3 = [0.65, 0.05, 0.20]; %user-specified
colors_4 = colors(4,:); %user-specified
N_colors1 = 50; %user-specified
N_colors2 = 50; %user-specified
colors_exp1 = [1, 1, 1]; %user-specified
colors_exp2 = [1, 1, 1]; %user-specified
colors_map1 = zeros(N_colors1, 3);
colors_map2 = zeros(N_colors2, 3);
for i = 1:3
    colors_map1(:,i) = colors_1(i) + (colors_2(i) - colors_1(i)) * linspace(0,1,N_colors1)'.^colors_exp1(i);
    colors_map2(:,i) = colors_3(i) + (colors_4(i) - colors_3(i)) * linspace(0,1,N_colors2)'.^colors_exp2(i);
end
colors_map1 = flipud(colors_map1); %flip order of colors
colors(1:3,:) = flipud(colors(1:3,:)); %flip colors so that warmer leaf temperatures are redder

color_grey = 0.70*[1, 1, 1];

%% Best parameters for theoretical predictions
N_best = N_best/10;
N_best(N_best < 100) = 100;

[~, ind_best_sub, ~, ~, ~, ~, ~, ~, ~, a_f_max, pi_L_star, beta, epsilon_L_max, SWC_L, alpha] = Choose_results(pi_L_0_25C_Trt_combo_results, SSE_best_total_subset_results, N_best, species_subset, a_f_max_results, pi_L_star_results, beta_results, epsilon_L_max_results, SWC_L_results, alpha_results, 'Option', 2); 
ind_best_sub = ind_best_sub(1);
pi_L_0_25C = mean(pi_L_0_25C_Trt_combo_results(ind_best_sub,:), 'omitnan'); 

%% Regress V_cmax25 as function of psi_L_MD
V_cmax25_func = Regress_Vcmax25_as_func_psi(T_L_subset, V_cmax_subset, psi_L_MD_subset); 

%% Figure: Instantaneous response to VPD at different leaf temperatures
VPD_L_min = 0.6; %[kPa]
VPD_L_max = 7; %[kPa]
VPD_L = VPD_L_min:0.01:VPD_L_max; 
T_L_lines = 27:-3:21;
N_x_axis = length(VPD_L);
N_lines_max = length(T_L_lines);

P_atm = P_atm(1) * ones(1, N_x_axis);
c_a = c_a(1) * ones(1, N_x_axis);
RWC_t = RWC_t(1) * ones(1, N_x_axis);

f_TLO_width_fraction = 0.70;
f_TLO_height_fraction = 0.80;
N_lines = 0;
DisplayName1 = cell(1, N_lines_max);

Fig1 = figure;
set(gcf, 'position', [100, 75, 1150, 460])
TLO = tiledlayout(1, 2, 'TileSpacing', 'tight', 'Unit', 'pixels');
pos_TLO = get(TLO, 'Position');
TLO_orig_width = pos_TLO(3);
TLO_orig_height = pos_TLO(4);
pos_TLO(3) = TLO_orig_width * f_TLO_width_fraction;
pos_TLO(4) = TLO_orig_height * f_TLO_height_fraction;
pos_TLO(2) = pos_TLO(2) + TLO_orig_height * (1 - f_TLO_height_fraction);
set(TLO, 'position', pos_TLO)

ax1 = nexttile(TLO, 1);
hold on
for i = 1:N_lines_max
    
    % leaf temperature for line
    T_L = T_L_lines(i) * ones(1, N_x_axis);
    
    % predict pressure-volume traits
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
    
    % local temperature-dependent photosynthetic parameters
    [outputs_photo_param] = Photosynthesis_parameters_temperature_response(T_L, psi_L, V_cmax25_func, ...
                                                                           'Nonmonotonic_V_cmax_temperature_response', Nonmonotonic_V_cmax_temperature_response);
    Gamma_star = outputs_photo_param.Gamma_star_vect;
    K_m = outputs_photo_param.K_m_vect;
    V_cmax = outputs_photo_param.V_cmax_vect;
    R_d = outputs_photo_param.R_d_vect;
    
    % predict stomatal conductance
    [outputs_stomata] = Stomata_for_RWC_and_Temperature(SWC_L, alpha, RWC_t, RWC_s, pi_L, pi_L_0, ...
                                                        a_f, da_f_dpi_L, a_f_0, da_f_0dpi_L_0, epsilon_L_0, ...
                                                        VPD_L, P_atm, T_L, V_cmax, ...
                                                        R_d, Gamma_star, K_m, c_a);


    g_w = outputs_stomata.g_w;
    g_w(outputs_stomata.solved == 0) = nan;
    
    if any(~isnan(g_w))
        
        N_lines = N_lines + 1;
        
        % calculate m = -dg_w/dln(VPD_L)/g_w_ref from negatively-sloped, primary VPD_L-g_w response
        g_w_polyfit = g_w;
        VPD_L_polyfit = VPD_L;
        log_VPD_L_polyfit = log(VPD_L_polyfit);
        is_not_nan = (~isnan(g_w_polyfit)) .* (~isnan(log_VPD_L_polyfit));
        dg_wdVPD_L_polyfit = [0, diff(g_w_polyfit)./diff(VPD_L_polyfit)];
        is_g_wdVPD_L_neg = (dg_wdVPD_L_polyfit < 0);
        is_both = is_not_nan .* is_g_wdVPD_L_neg; 
        g_w_polyfit = g_w_polyfit(is_both == 1);
        log_VPD_L_polyfit = log_VPD_L_polyfit(is_both == 1); 
        p = polyfit(log_VPD_L_polyfit, g_w_polyfit, 1);
        g_w_ref = p(2); %reference g_w at VPD_L = 1 kPa
        m = - p(1)/g_w_ref; %m = -dg_w/dln(VPD_L)/g_w_ref
        
        % % % DisplayName1{N_lines} = ['{\itT_L} = ', num2str(T_L_lines(i)), ' ', char(176), 'C, {\itm} = ', num2str(m, '%.2f')];
        DisplayName1{N_lines} = ['{\itT_L} = ', num2str(T_L_lines(i)), ' ', char(176), 'C'];
        
        plot(VPD_L, g_w, '-', 'linewidth', 1, 'Color', colors(N_lines,:), 'DisplayName', DisplayName1{N_lines});
        
    end
end
hold off
xlabel('{\itD_L} [kPa]')
ylabel('{\itg_w} [mol\cdotm^{-2}\cdots^{-1}]')
ylim_max = max(ylim);
ylim([0, ylim_max])
pause(0)


%% Plot steady-state response between g_w, T_L, and VPD, varying VPD by changing T_a (with leaf water potential calculated from steady-state hydraulics)
N_x_axis = 1e2;
T_a_min = 0;
T_a_max = 50;
T_a_varied = linspace(T_a_min, T_a_max, N_x_axis);

T_a = T_a(1) * ones(1, N_x_axis);
RH = RH(1) * ones(1, N_x_axis);
P_atm = P_atm(1) * ones(1, N_x_axis);
c_a = c_a(1) * ones(1, N_x_axis);
RWC_t = RWC_t(1) * ones(1, N_x_axis);

progress_disp_inc = 0.02; 
progress_disp = progress_disp_inc;

%% determine optimum temperature for photosynthesis 
T_L_opt = 0:0.01:50;
V_cmax_per_V_cmax25_varied = exp(65330*(T_L_opt-25)/298.15/R./(T_L_opt+273.15));
V_cmax_per_V_cmax25_varied = V_cmax_per_V_cmax25_varied .* (1 + exp((298.15*640 - 2e5)/R/298.15)) ./ (1 + exp(((T_L_opt+273.15)*640 - 2e5)/R./(T_L_opt+273.15))); %modify to include temperature optimum at 36.4C
T_L_opt = T_L_opt(V_cmax_per_V_cmax25_varied == max(V_cmax_per_V_cmax25_varied));

% vectors to store calculated data
g_w_with_varied_T_a = nan(1, N_x_axis);
VPD_L_with_varied_T_a = nan(1, N_x_axis);
T_L_with_varied_T_a = nan(1, N_x_axis);
E_with_varied_T_a = nan(1, N_x_axis);
E_crit_with_varied_T_a = nan(1, N_x_axis);
psi_L_with_varied_T_a = nan(1, N_x_axis);
psi_L_crit_with_varied_T_a = nan(1, N_x_axis);

for i = 1:N_x_axis
    
    % display progress
    progress = i/N_x_axis;
    if progress >= progress_disp
        disp(['        Calculating steady-state g_w for varied T_a: ', num2str(100*progress_disp_inc*floor(progress/progress_disp_inc), '%.0f'), '% done!'])
        progress_disp = progress_disp_inc * (1 + floor(progress/progress_disp_inc));
    end
    
    [outputs_stomata_steady_state] = ...
    Stomata_for_steady_state_thermo_and_hydraulics(SWC_L, alpha, a_f_max, pi_L_star, beta, epsilon_L_max, pi_L_0_25C, ...
                                                   k_L_max_25C, A, Q_10_k_L, psi_soil, ...
                                                   T_a_varied(i), RH(i), P_atm(i), c_a(i), ...
                                                   R_abs, g_H_a, emiss_L, ...
                                                   V_cmax25_func, ...
                                                   'Nonmonotonic_V_cmax_temperature_response', 1);

    g_w_with_varied_T_a(i) = outputs_stomata_steady_state.g_w_steady_state;
    VPD_L_with_varied_T_a(i) = outputs_stomata_steady_state.VPD_L_steady_state;
    T_L_with_varied_T_a(i) = outputs_stomata_steady_state.T_L_steady_state;
    E_with_varied_T_a(i) = outputs_stomata_steady_state.E_steady_state;
    E_crit_with_varied_T_a(i) = outputs_stomata_steady_state.E_crit_steady_state;
    psi_L_with_varied_T_a(i) = outputs_stomata_steady_state.psi_L_steady_state;
    psi_L_crit_with_varied_T_a(i) = outputs_stomata_steady_state.psi_L_crit_steady_state;
    
end

T_L_opt_lb = max(T_L_with_varied_T_a(T_L_with_varied_T_a <= T_L_opt));
T_L_opt_ub = min(T_L_with_varied_T_a(T_L_with_varied_T_a >= T_L_opt));
T_a_opt = T_a_varied(T_L_with_varied_T_a == T_L_opt_lb) + (T_a_varied(T_L_with_varied_T_a == T_L_opt_ub) - T_a_varied(T_L_with_varied_T_a == T_L_opt_lb))*(T_L_opt - T_L_opt_lb)/(T_L_opt_ub - T_L_opt_lb);

% calculate m = -dg_w/dln(VPD_L)/g_w_ref from negatively-sloped, primary VPD_L-g_w response
g_w_polyfit = g_w_with_varied_T_a;
VPD_L_polyfit = VPD_L_with_varied_T_a;
log_VPD_L_polyfit = log(VPD_L_polyfit);
is_not_nan = (~isnan(g_w_polyfit)) .* (~isnan(log_VPD_L_polyfit));
dg_wdVPD_L_polyfit = [0, diff(g_w_polyfit)./diff(VPD_L_polyfit)];
is_g_wdVPD_L_neg = (dg_wdVPD_L_polyfit < 0);
is_both = is_not_nan .* is_g_wdVPD_L_neg; 
g_w_polyfit = g_w_polyfit(is_both == 1);
log_VPD_L_polyfit = log_VPD_L_polyfit(is_both == 1); 
p = polyfit(log_VPD_L_polyfit, g_w_polyfit, 1);
g_w_ref = p(2); %reference g_w at V{D_L = 1 kPa
m = - p(1)/g_w_ref; %m = -dg_w/dln(VPD_L)/g_w_ref

% % % DisplayName2 = {['Varied {\itT_a} & constant RH', newline, '({\itsteady-state} xylem hydraulics),', newline, '{\itm} = ', num2str(m, '%.2f')]};
DisplayName2 = {'{\itSteady-state} heat & water'};
                 
pause(0)
%% Plot steady-state response between g_w, T_L, and VPD, varying VPD by changing RH (with leaf water potential calculated from steady-state hydraulics)
N_x_axis = 1e2;
RH_min = 0.1;
RH_max = 0.9;
RH_varied = linspace(RH_min, RH_max, N_x_axis);

progress_disp_inc = 0.02; 
progress_disp = progress_disp_inc;

% vectors to store calculated data
g_w_with_varied_RH = nan(1, N_x_axis);
VPD_L_with_varied_RH = nan(1, N_x_axis);
T_L_with_varied_RH = nan(1, N_x_axis);
E_with_varied_RH = nan(1, N_x_axis);

for i = 1:N_x_axis
    
    % display progress
    progress = i/N_x_axis;
    if progress >= progress_disp
        disp(['        Calculating steady-state g_w for varied RH: ', num2str(100*progress_disp_inc*floor(progress/progress_disp_inc), '%.0f'), '% done!'])
        progress_disp = progress_disp_inc * (1 + floor(progress/progress_disp_inc));
    end
    
    [outputs_stomata_steady_state] = ...
    Stomata_for_steady_state_thermo_and_hydraulics(SWC_L, alpha, a_f_max, pi_L_star, beta, epsilon_L_max, pi_L_0_25C, ...
                                                   k_L_max_25C, A, Q_10_k_L, psi_soil, ...
                                                   T_a(i), RH_varied(i), P_atm(i), c_a(i), ...
                                                   R_abs, g_H_a, emiss_L, ...
                                                   V_cmax25_func, ...
                                                   'Nonmonotonic_V_cmax_temperature_response', 1);

    g_w_with_varied_RH(i) = outputs_stomata_steady_state.g_w_steady_state;
    VPD_L_with_varied_RH(i) = outputs_stomata_steady_state.VPD_L_steady_state;
    T_L_with_varied_RH(i) = outputs_stomata_steady_state.T_L_steady_state;
    E_with_varied_RH(i) = outputs_stomata_steady_state.E_steady_state;
    E_crit_with_varied_RH(i) = outputs_stomata_steady_state.E_crit_steady_state;
    psi_L_with_varied_RH(i) = outputs_stomata_steady_state.psi_L_steady_state;
    psi_L_crit_with_varied_RH(i) = outputs_stomata_steady_state.psi_L_crit_steady_state;

end

% calculate m = -dg_w/dln(VPD_L)/g_w_ref from negatively-sloped, primary VPD_L-g_w response
g_w_polyfit = g_w_with_varied_RH;
VPD_L_polyfit = VPD_L_with_varied_RH;
log_VPD_L_polyfit = log(VPD_L_polyfit);
is_not_nan = (~isnan(g_w_polyfit)) .* (~isnan(log_VPD_L_polyfit));
dg_wdVPD_L_polyfit = [0, diff(g_w_polyfit)./diff(VPD_L_polyfit)];
is_g_wdVPD_L_neg = (dg_wdVPD_L_polyfit < 0);
is_both = is_not_nan .* is_g_wdVPD_L_neg; 
g_w_polyfit = g_w_polyfit(is_both == 1);
log_VPD_L_polyfit = log_VPD_L_polyfit(is_both == 1); 
p = polyfit(log_VPD_L_polyfit, g_w_polyfit, 1);
g_w_ref = p(2); %reference g_w at V{D_L = 1 kPa
m = - p(1)/g_w_ref; %m = -dg_w/dln(VPD_L)/g_w_ref

% % % DisplayName3 = {['{\itSteady-state} heat &' newline, 'water, {\itm} = ', num2str(m, '%.2f')]};
DisplayName3 = {'{\itSteady-state} heat & water'};

hold on
for i = 1:N_lines
    P1(i) = plot(-1, 0, '-', 'linewidth', 1, 'Color', colors(i,:), 'DisplayName', DisplayName1{i});
end
P3 = plot(-1, 0, '-', 'linewidth', 6, 'color', colors_map1(N_colors1/2, :), 'DisplayName', DisplayName3{1});
surface([VPD_L_with_varied_RH; VPD_L_with_varied_RH],[g_w_with_varied_RH; g_w_with_varied_RH], zeros(2, N_x_axis), [RH_varied; RH_varied], 'facecol', 'no', 'edgecol', 'interp', 'linewidth', 6, 'linestyle', '-')
hold off
ylim([0, ylim_max])
colormap(ax1, colors_map1)
caxis(ax1, [RH_min, RH_max])
title(['RH response', newline, 'with new model'], 'fontsize', 11)
xlim([VPD_L_min, 5.3])
set(gca, 'xtick', 1:5)

lgd_1 = legend(ax1, [P1(1:N_lines), P3], [DisplayName1(1:N_lines), DisplayName3], 'Location', 'northeast', 'fontsize', 9, 'Color', 'w', 'EdgeColor', 'k', 'box', 'on', 'Units', 'Pixels');
title(lgd_1, ['New model with varied RH'])
pos_lgd_1 = get(lgd_1, 'Position');

pause(0)
%% Plot steady-state response between g_w, T_L, and VPD, varying VPD by changing T_a (with constant leaf water potential of zero)
N_x_axis = 1e2;
T_a_min = 0;
T_a_max = 50;
T_a_varied = linspace(T_a_min, T_a_max, N_x_axis);

T_L_minus_T_a_vect = -10:0.01:10;
N_vect = length(T_L_minus_T_a_vect);

progress_disp_inc = 0.02; 
progress_disp = progress_disp_inc;

% vectors to store calculated data
g_w_with_varied_T_a_constant_RWC = nan(1, N_x_axis);
VPD_L_with_varied_T_a_constant_RWC = nan(1, N_x_axis);
T_L_with_varied_T_a_constant_RWC = nan(1, N_x_axis);
E_with_varied_T_a_constant_RWC = nan(1, N_x_axis);
chi_w_with_varied_T_a_constant_RWC = nan(1, N_x_axis);

for i = 1:N_x_axis
    
    % display progress
    progress = i/N_x_axis;
    if progress >= progress_disp
        disp(['        Calculating steady-state g_w for varied T_a (with constant psi_L): ', num2str(100*progress_disp_inc*floor(progress/progress_disp_inc), '%.0f'), '% done!'])
        progress_disp = progress_disp_inc * (1 + floor(progress/progress_disp_inc));
    end
    
    % determine g_w for the given T_a and RH
    T_L_vect = T_a_varied(i) + T_L_minus_T_a_vect;
    RH_vect = RH(1) * ones(1, N_vect);
    e_L_vect = 0.61078 * exp(17.27 * T_L_vect ./ (T_L_vect + 237.3)); % saturated leaf vapor pressure in [kPa] -- Teten's equation

    % local temperature-dependent photosynthetic parameters
    [outputs_photo_param] = Photosynthesis_parameters_temperature_response(T_L_vect, 0, V_cmax25_func, ...
                                                                           'Nonmonotonic_V_cmax_temperature_response', Nonmonotonic_V_cmax_temperature_response);
    Gamma_star_vect = outputs_photo_param.Gamma_star_vect;
    K_m_vect = outputs_photo_param.K_m_vect;
    V_cmax_vect = outputs_photo_param.V_cmax_vect;
    R_d_vect = outputs_photo_param.R_d_vect;

    % determine secondary estimate of g_w from T_L and VPD_L
    T_a_vect = T_a_varied(i) * ones(1, N_vect);
    e_a_vect = RH_vect .* 0.61078 .* exp(17.27 * T_a_vect ./ (T_a_vect + 237.3)); %air vapor pressure in [kPa] -- Teten's equation
    VPD_L_vect = e_L_vect - e_a_vect; %vapor pressure deficit in [kPa]
    P_atm_vect = P_atm(i) * ones(1, N_vect);
    c_a_vect = c_a(i) * ones(1, N_vect);
    pi_L_0_25C_vect = pi_L_0_25C * ones(1, N_vect);
    RWC_t_vect = RWC_t(i) * ones(1, N_vect);
    
	% predict pressure-volume traits from T_L and psi_L
    [outputs_PV] = PV_from_RWC_t_for_Trt( a_f_max, pi_L_star, beta, epsilon_L_max, ...
                                          pi_L_0_25C_vect, T_L_vect, RWC_t_vect);
    
    RWC_s_vect = outputs_PV.RWC_s;
    pi_L_vect = outputs_PV.pi_L;
    pi_L_0_vect = outputs_PV.pi_L_0;
    epsilon_L_0_vect = outputs_PV.epsilon_L_0;
    a_f_vect = outputs_PV.a_f;
    a_f_0_vect = outputs_PV.a_f_0;
    da_f_dpi_L_vect = outputs_PV.da_f_dpi_L;
    da_f_0dpi_L_0_vect = outputs_PV.da_f_0dpi_L_0;
    
    [outputs_stomata] = Stomata_for_RWC_and_Temperature(SWC_L, alpha, RWC_t_vect, RWC_s_vect, pi_L_vect, pi_L_0_vect, ...
                                                        a_f_vect, da_f_dpi_L_vect, a_f_0_vect, da_f_0dpi_L_0_vect, epsilon_L_0_vect, ...
                                                        VPD_L_vect, P_atm_vect, T_L_vect, V_cmax_vect, ...
                                                        R_d_vect, Gamma_star_vect, K_m_vect, c_a_vect);


    g_w_vect_2 = outputs_stomata.g_w;
    g_w_vect_2(outputs_stomata.solved == 0) = nan;
    E_vect = outputs_stomata.E;
    E_vect(outputs_stomata.solved == 0) = 0;
    chi_w_vect = outputs_stomata.chi_w;
    chi_w_vect(outputs_stomata.solved == 0) = nan;
    
	% determine T_L, VPD_L, and g_w from steady-state thermodynamics
    [~, ~, g_w_vect_1] = Steady_state_thermo(T_a_vect, E_vect, RH_vect, P_atm_vect, R_abs, g_H_a, emiss_L);
    g_w_vect_1(g_w_vect_1 < 0) = nan;
    
    % find intersection of two estimates of g_w
    g_w_error_vect = g_w_vect_1 - g_w_vect_2;
    diff_sign_g_w_error_vect = diff(sign(g_w_error_vect));
    diff_sign_g_w_error_vect(isnan(diff_sign_g_w_error_vect)) = 0;
    not_zero_diff_sign_g_w_error_vect = (diff_sign_g_w_error_vect ~= 0);
    if any(not_zero_diff_sign_g_w_error_vect)
        ind_LB = find(not_zero_diff_sign_g_w_error_vect, 1, 'first');
        ind_UB = ind_LB + 1;
        g_w_error_LB = g_w_error_vect(ind_LB);
        g_w_error_UB = g_w_error_vect(ind_UB);
        g_w_with_varied_T_a_constant_RWC(i) = g_w_vect_1(ind_LB) + (g_w_vect_1(ind_UB) - g_w_vect_1(ind_LB)) * (0 - g_w_error_LB)/(g_w_error_UB - g_w_error_LB);
        VPD_L_with_varied_T_a_constant_RWC(i) = VPD_L_vect(ind_LB) + (VPD_L_vect(ind_UB) - VPD_L_vect(ind_LB)) * (0 - g_w_error_LB)/(g_w_error_UB - g_w_error_LB);
        T_L_with_varied_T_a_constant_RWC(i) = T_L_vect(ind_LB) + (T_L_vect(ind_UB) - T_L_vect(ind_LB)) * (0 - g_w_error_LB)/(g_w_error_UB - g_w_error_LB);
        E_with_varied_T_a_constant_RWC(i) = E_vect(ind_LB) + (E_vect(ind_UB) - E_vect(ind_LB)) * (0 - g_w_error_LB)/(g_w_error_UB - g_w_error_LB);
        chi_w_with_varied_T_a_constant_RWC(i) = chi_w_vect(ind_LB) + (chi_w_vect(ind_UB) - chi_w_vect(ind_LB)) * (0 - g_w_error_LB)/(g_w_error_UB - g_w_error_LB);
    end
    
end

% calculate m = -dg_w/dln(VPD_L)/g_w_ref from negatively-sloped, primary VPD_L-g_w response
g_w_polyfit = g_w_with_varied_T_a_constant_RWC;
VPD_L_polyfit = VPD_L_with_varied_T_a_constant_RWC;
log_VPD_L_polyfit = log(VPD_L_polyfit);
is_not_nan = (~isnan(g_w_polyfit)) .* (~isnan(log_VPD_L_polyfit));
dg_wdVPD_L_polyfit = [0, diff(g_w_polyfit)./diff(VPD_L_polyfit)];
is_g_wdVPD_L_neg = (dg_wdVPD_L_polyfit < 0);
is_both = is_not_nan .* is_g_wdVPD_L_neg; 
g_w_polyfit = g_w_polyfit(is_both == 1);
log_VPD_L_polyfit = log_VPD_L_polyfit(is_both == 1); 
p = polyfit(log_VPD_L_polyfit, g_w_polyfit, 1);
g_w_ref = p(2); %reference g_w at V{D_L = 1 kPa
m = - p(1)/g_w_ref; %m = -dg_w/dln(VPD_L)/g_w_ref

% % % DisplayName4 = {['Varied {\itT_a} & constant RH', newline, '(fully hydrated), {\itm} = ', num2str(m, '%.2f')]};
DisplayName4 = {['Fully hydrated ({\itsteady-state}', newline, 'heat & fixed RWC)']};


pause(0)
%% Plot steady-state response between g_w, T_L, and VPD, varying VPD by changing T_a (with constant chi_w)

chi_w_with_varied_T_a_constant_RWC(chi_w_with_varied_T_a_constant_RWC <= 0) = nan;
chi_w_mean = mean(chi_w_with_varied_T_a_constant_RWC, 'omitnan');

N_x_axis = 1e2;
T_a_min = 0;
T_a_max = 50;
T_a_varied = linspace(T_a_min, T_a_max, N_x_axis);

% determine g_w for the given T_a and RH
PLC_max = 0.999;
psi_L_min = log(1-PLC_max)/A;
psi_L_vect = 0:-0.01:psi_L_min;
N_vect = length(psi_L_vect);
RH_vect = RH(1) * ones(1, N_vect);

progress_disp_inc = 0.02; 
progress_disp = progress_disp_inc;

% vectors to store calculated data
g_w_with_varied_T_a_constant_chi_w = nan(1, N_x_axis);
VPD_L_with_varied_T_a_constant_chi_w = nan(1, N_x_axis);
T_L_with_varied_T_a_constant_chi_w = nan(1, N_x_axis);
E_with_varied_T_a_constant_chi_w = nan(1, N_x_axis);

for i = 1:N_x_axis
    
    % display progress
    progress = i/N_x_axis;
    if progress >= progress_disp
        disp(['        Calculating steady-state g_w for varied T_a (with constant chi_w): ', num2str(100*progress_disp_inc*floor(progress/progress_disp_inc), '%.0f'), '% done!'])
        progress_disp = progress_disp_inc * (1 + floor(progress/progress_disp_inc));
    end

    [outputs_stomata_steady_state] = ...
    Stomata_for_steady_state_thermo_and_hydraulics(SWC_L, alpha, a_f_max, pi_L_star, beta, epsilon_L_max, pi_L_0_25C, ...
                                                   k_L_max_25C, A, Q_10_k_L, psi_soil, ...
                                                   T_a_varied(i), RH(i), P_atm(i), c_a(i), ...
                                                   R_abs, g_H_a, emiss_L, ...
                                                   V_cmax25_func, 'Constant_chi_w', chi_w_mean, ...
                                                   'Nonmonotonic_V_cmax_temperature_response', 1);
    
    g_w_with_varied_T_a_constant_chi_w(i) = outputs_stomata_steady_state.g_w_steady_state;
    VPD_L_with_varied_T_a_constant_chi_w(i) = outputs_stomata_steady_state.VPD_L_steady_state;
    T_L_with_varied_T_a_constant_chi_w(i) = outputs_stomata_steady_state.T_L_steady_state;
    E_with_varied_T_a_constant_chi_w(i) = outputs_stomata_steady_state.E_steady_state;
    
end

% calculate m = -dg_w/dln(VPD_L)/g_w_ref from negatively-sloped, primary VPD_L-g_w response
g_w_polyfit = g_w_with_varied_T_a_constant_chi_w;
VPD_L_polyfit = VPD_L_with_varied_T_a_constant_chi_w;
log_VPD_L_polyfit = log(VPD_L_polyfit);
is_not_nan = (~isnan(g_w_polyfit)) .* (~isnan(log_VPD_L_polyfit));
dg_wdVPD_L_polyfit = [0, diff(g_w_polyfit)./diff(VPD_L_polyfit)];
is_g_wdVPD_L_neg = (dg_wdVPD_L_polyfit < 0);
is_both = is_not_nan .* is_g_wdVPD_L_neg; 
g_w_polyfit = g_w_polyfit(is_both == 1);
log_VPD_L_polyfit = log_VPD_L_polyfit(is_both == 1); 
p = polyfit(log_VPD_L_polyfit, g_w_polyfit, 1);
g_w_ref = p(2); %reference g_w at V{D_L = 1 kPa
m = - p(1)/g_w_ref; %m = -dg_w/dln(VPD_L)/g_w_ref

% % % DisplayName5 = {['Varied {\itT_a} & constant RH', newline, '(constant {\it\chi_w}), {\itm} = ', num2str(m, '%.2f')]};
% % % DisplayName5 = {['Varied {\itT_a} & constant RH', newline, '(constant {\it\chi_w})']};
DisplayName5 = {['Cowan & Farquhar', newline, '(1977; constant {\it\chi_w})']};


pause(0)
%% Plot steady-state response between g_w, T_L, and VPD, varying VPD by changing T_a (by Sperry et al., 2017)

N_x_axis = 1e2;
T_a_min = 0;
T_a_max = 50;
T_a_varied = linspace(T_a_min, T_a_max, N_x_axis);

P_atm = P_atm(1) * ones(1, N_x_axis);
c_a = c_a(1) * ones(1, N_x_axis);
RH = RH(1) * ones(1, N_x_axis);

progress_disp_inc = 0.02; 
progress_disp = progress_disp_inc;

% vectors to store calculated data
g_w_with_varied_T_a_Sperry = nan(1, N_x_axis);
VPD_L_with_varied_T_a_Sperry = nan(1, N_x_axis);
T_L_with_varied_T_a_Sperry = nan(1, N_x_axis);
E_with_varied_T_a_Sperry = nan(1, N_x_axis);

for i = 1:N_x_axis
    
    % display progress
    progress = i/N_x_axis;
    if progress >= progress_disp
        disp(['        Calculating steady-state g_w for varied T_a (Sperry et al., 2017): ', num2str(100*progress_disp_inc*floor(progress/progress_disp_inc), '%.0f'), '% done!'])
        progress_disp = progress_disp_inc * (1 + floor(progress/progress_disp_inc));
    end
    
    [outputs_stomata_Sperry] = ...
    Stomata_Sperry_Gain_Risk(k_L_max_25C, A, Q_10_k_L, psi_soil, ...
                             T_a_varied(i), RH(i), P_atm(i), c_a(i), ...
                             R_abs, g_H_a, emiss_L, ...
                             V_cmax25_func, ...
                             'Nonmonotonic_V_cmax_temperature_response', 1);
	
	g_w_with_varied_T_a_Sperry(i) = outputs_stomata_Sperry.g_w_Sperry;
    VPD_L_with_varied_T_a_Sperry(i) = outputs_stomata_Sperry.VPD_L_Sperry;
    T_L_with_varied_T_a_Sperry(i) = outputs_stomata_Sperry.T_L_Sperry;
    E_with_varied_T_a_Sperry(i) = outputs_stomata_Sperry.E_Sperry;
    
end

% calculate m = -dg_w/dln(VPD_L)/g_w_ref from negatively-sloped, primary VPD_L-g_w response
g_w_polyfit = g_w_with_varied_T_a_Sperry;
VPD_L_polyfit = VPD_L_with_varied_T_a_Sperry;
log_VPD_L_polyfit = log(VPD_L_polyfit);
is_not_nan = (~isnan(g_w_polyfit)) .* (~isnan(log_VPD_L_polyfit));
dg_wdVPD_L_polyfit = [0, diff(g_w_polyfit)./diff(VPD_L_polyfit)];
is_g_wdVPD_L_neg = (dg_wdVPD_L_polyfit < 0);
is_both = is_not_nan .* is_g_wdVPD_L_neg; 
g_w_polyfit = g_w_polyfit(is_both == 1);
log_VPD_L_polyfit = log_VPD_L_polyfit(is_both == 1); 
p = polyfit(log_VPD_L_polyfit, g_w_polyfit, 1);
g_w_ref = p(2); %reference g_w at V{D_L = 1 kPa
m = - p(1)/g_w_ref; %m = -dg_w/dln(VPD_L)/g_w_ref

% % % DisplayName6 = {['Varied {\itT_a} & constant RH', newline, '(Sperry et al., 2017), {\itm} = ', num2str(m, '%.2f')]};
% % % DisplayName6 = {['Varied {\itT_a} & constant RH', newline, '(Gain-Risk)']};
DisplayName6 = {['Sperry et al.', newline, '(2017; Gain-Risk)']};


pause(0)
%% Plot simulated data with varied air temperature

ax2 = nexttile(TLO, 2);
hold on
P2 = plot(-1, 0, '-', 'linewidth', 6, 'color', colors_map2(N_colors2/2, :), 'DisplayName', DisplayName2{1});
P4 = plot(VPD_L_with_varied_T_a_constant_RWC, g_w_with_varied_T_a_constant_RWC, '-', 'linewidth', 1.2, 'color', colors_1, 'DisplayName', DisplayName4{1});
P5 = plot(VPD_L_with_varied_T_a_constant_chi_w, g_w_with_varied_T_a_constant_chi_w, ':', 'linewidth', 1.2, 'color', color_grey, 'DisplayName', DisplayName5{1});
P6 = plot(VPD_L_with_varied_T_a_Sperry(g_w_with_varied_T_a_Sperry > 0), g_w_with_varied_T_a_Sperry(g_w_with_varied_T_a_Sperry > 0), '--', 'linewidth', 1.2, 'color', color_grey, 'DisplayName', DisplayName6{1});
surface([VPD_L_with_varied_T_a; VPD_L_with_varied_T_a],[g_w_with_varied_T_a; g_w_with_varied_T_a], zeros(2, N_x_axis), [T_a_varied; T_a_varied], 'facecol', 'no', 'edgecol', 'interp', 'linewidth', 6, 'linestyle', '-');
hold off
ylim([0, ylim_max])
xlim([VPD_L_min, 5.3])
set(gca, 'xtick', 1:5)
xlabel('{\itD_L} [kPa]')
title(['Intermodel comparison of', newline, 'air temperature responses'], 'fontsize', 11)

% limits of air temperature for caxis
dT_a_min_plot = 5;
dT_a_max_plot = dT_a_min_plot;
T_a_min_plot = min(T_a_varied(~isnan(g_w_with_varied_T_a)));
T_a_max_plot = max(T_a_varied(~isnan(g_w_with_varied_T_a)));
T_a_min_plot = dT_a_min_plot * floor(T_a_min_plot/dT_a_min_plot);
T_a_max_plot = dT_a_max_plot * ceil(T_a_max_plot/dT_a_max_plot);

colormap(ax2, colors_map2)
caxis(ax2, [T_a_min_plot, T_a_max_plot])

% add colors bars
cb1 = colorbar(ax1, 'SouthOutside');
ylabel(cb1, 'RH', 'FontSize', 10)
set(cb1, 'ytick', RH_min:0.2:RH_max)
cb2 = colorbar(ax2, 'SouthOutside');
ylabel(cb2, ['{\itT_a} [', char(176), 'C]'], 'FontSize', 10)

rel_alphab_y = 1.08;
rel_alphab_x = 0.06;
for i = 1:2
    nexttile(i)
	LETTER = ['(', char(96 + i), ')'];
	text(VPD_L_min + rel_alphab_x*(max(xlim) - VPD_L_min), 0 + rel_alphab_y*(ylim_max - 0), ...
         LETTER, 'fontsize', 11, 'fontweight', 'bold', 'BackgroundColor', 'w') 
end


% legend
lgd_2 = legend(ax2, [P2, P4, P5, P6], [DisplayName2, DisplayName4, DisplayName5, DisplayName6], 'Location', 'northeastoutside', 'fontsize', 9, 'Color', 'w', 'EdgeColor', 'k', 'box', 'on', 'Units', 'Pixels');
title(lgd_2, 'Models with varied {\itT_a}')
pos_lgd_2 = get(lgd_2, 'Position');
% % % set(lgd_2, 'Position', pos_lgd_2)

% third axis
ax3 = axes('fontsize', 8);
pos_ax2 = get(ax2, 'position');
pos_ax3 = [pos_ax2(1) + 0.0*pos_ax2(3), ...
           pos_ax2(2) + 0.6*pos_ax2(4), ...
           0.25*pos_ax2(3), ...
           0.3*pos_ax2(4)];
set(ax3, 'position', pos_ax3)
hold on
plot([0, T_a_max_plot], [0, 0], 'k')
plot(T_a_varied , g_w_with_varied_T_a_constant_chi_w - g_w_with_varied_T_a, ':', 'linewidth', 1.2, 'color', color_grey)
plot(T_a_varied , g_w_with_varied_T_a_Sperry - g_w_with_varied_T_a, '--', 'linewidth', 1.2, 'color', color_grey)
hold off
xlabel(['{\itT_a} [', char(176), 'C]'], 'fontsize', 8)
ylabel('\Delta{\itg_w}', 'fontsize', 8)
xlim([0, T_a_max_plot])
set(ax3, 'xtick', 0:25:T_a_max_plot)
title(['Differences from', newline, 'new model'], 'fontsize', 8)

% save figure
Output_file1 = [species_subset, '_Fig_VPD_ss_response_gw_combined.jpg'];
print(gcf, '-djpeg', '-r300', Output_file1, '-painters', '-noui' )
exportgraphics(gcf, [species_subset, '_Fig_VPD_ss_response_gw_combined.pdf'], 'ContentType', 'vector')

%% Separate plot 
figure
tiledlayout(3, 1, 'TileSpacing', 'tight');

nexttile(1)
plot(T_a_varied, T_L_with_varied_T_a - T_a_varied, '-', 'linewidth', 6, 'color', colors_map2(N_colors2/2, :))
xlabel(['{\itT_a} [', char(176), 'C]'])
ylabel(['{\itT_L} - {\itT_a} [', char(176), 'C]'])
xlim([T_a_min_plot, T_a_max_plot])
dylim_max = 1;
dylim_min = dylim_max;
ylim_max = max(ylim);
ylim_max = dylim_max * ceil(ylim_max/dylim_max);
if min(T_L_with_varied_T_a - T_a_varied) > 0
    ylim_min = 0;
else
    ylim_min = min(ylim);
    ylim_min = dylim_min * floor(ylim_min/dylim_min);
end
ylim([ylim_min, ylim_max])
ylim([0, 3]) %ad hoc adjustment

hold on
plot(T_a_opt*[1, 1], [min(ylim), max(ylim)], '--', 'linewidth', 2, 'color', colors_3)
hold off

nexttile(2)
hold on
plot(T_a_varied, 1e3 * E_with_varied_T_a, '-', 'linewidth', 6, 'color', colors_map2(N_colors2/2, :))
plot(T_a_varied, 1e3 * E_crit_with_varied_T_a, '--', 'linewidth', 2, 'color', colors_1)
hold off
xlabel(['{\itT_a} [', char(176), 'C]'])
ylabel('{\itE} [mmol\cdotm^{-1}\cdots^{-1}]')
xlim([T_a_min_plot, T_a_max_plot])
dylim_max = 1;
ylim_max = max(ylim);
ylim_max = dylim_max * ceil(ylim_max/dylim_max);
ylim([0, ylim_max])

hold on
plot(T_a_opt*[1, 1], [min(ylim), max(ylim)], '--', 'linewidth', 2, 'color', colors_3)
hold off

nexttile(3)
hold on
plot(T_a_varied, -psi_L_with_varied_T_a, '-', 'linewidth', 6, 'color', colors_map2(N_colors2/2, :))
plot(T_a_varied, -psi_L_crit_with_varied_T_a, '--', 'linewidth', 2, 'color', colors_1)
hold off
xlabel(['{\itT_a} [', char(176), 'C]'])
ylabel('-{\it\psi_L} [MPa]')
xlim([T_a_min_plot, T_a_max_plot])
set(gca, 'yscale', 'log')
ylim([0.01, 10])
set(gca, 'ytick', 10.^(-2:1))

hold on
plot(T_a_opt*[1, 1], [min(ylim), max(ylim)], '--', 'linewidth', 2, 'color', colors_3)
hold off


rel_alphab_y = 0.95;
rel_alphab_x = 0.06;
for i = 1:3
    nexttile(i)
    xlim_min = min(xlim);
    xlim_max = max(xlim);
    ylim_min = min(ylim);
    ylim_max = max(ylim);
    LETTER = ['(', char(96 + i), ')'];
    text(xlim_min + rel_alphab_x*(xlim_max - xlim_min), ylim_min + rel_alphab_y*(ylim_max - ylim_min), ...
         LETTER, 'fontsize', 10, 'fontweight', 'bold', 'BackgroundColor', 'w')
    box off
end

% save figure
Output_file2 = [species_subset, '_Fig_VPD_ss_response_T_L_and_E.jpg'];
print(gcf, '-djpeg', '-r300', Output_file2, '-painters', '-noui' )

end

