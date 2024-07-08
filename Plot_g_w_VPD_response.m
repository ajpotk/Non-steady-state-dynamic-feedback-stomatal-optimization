function [] = Plot_g_w_VPD_response(a_f_max_results, pi_L_star_results, beta_results, epsilon_L_max_results, ...
                                    SWC_L_results, alpha_results, pi_L_0_25C_Trt_combo_results, SSE_best_total_subset_results, ...
                                    N_best, species_subset, T_L_subset, psi_L_MD_subset, V_cmax_subset)

%% Constants, parameters, and conditions
o_i = 1e-3 * 210; %standard atmospheric partial pressure of oxygen [mol mol-1]
R = 8.314; %universal gas constant [J mol-1 K-1]
R_d25_per_V_cmax25 = 0.015; %from Collatz et al. (1991)
emiss_L = 0.97; %leaf emissivity [-]
g_H_a = 2; %boundary layer heat conductance [mol m-2 s-1]
R_abs = 400; %absorbed shortwave radiation [W m-2]

%% Parameter estimates with lowest RMSE
[~, ind_best_sub, ~, ~, ~, ~, ~, ~, ~, a_f_max, pi_L_star, beta, epsilon_L_max, SWC_L, alpha] = Choose_results(pi_L_0_25C_Trt_combo_results, SSE_best_total_subset_results, N_best, species_subset, a_f_max_results, pi_L_star_results, beta_results, epsilon_L_max_results, SWC_L_results, alpha_results);

%% Default environmental conditions
default.P_atm = 101.325; %atmospheric pressure [kPa]
default.c_a = 4e-4; %atmospheric CO2 concentration [mol mol-1]
default.RWC_t = 1; %leaf total relative water content [-]
ind_lowest_SSE = ind_best_sub(1);
default.pi_L_0_25C = mean(pi_L_0_25C_Trt_combo_results(ind_lowest_SSE,:), 'omitnan'); %osmotic potential at full hydration [MPa] -- taken as average among treatments

%% Regress V_cmax25 as function of psi_L_MD
V_cmax25_func = Regress_Vcmax25_as_func_psi(T_L_subset, V_cmax_subset, psi_L_MD_subset);

%% Figure colors
colors = [102, 102, 255;...
          124, 203, 161; ...
          240, 116, 110; ...
          253, 222, 156; ...
          002, 144, 153; ...
          220, 057, 120; ...
          006, 082, 117]/255;

%% Figure: Instantaneous response to VPD at different leaf temperatures
VPD_L = 0.6:0.01:8; 
T_L_lines = 21:3:27;
N_x_axis = length(VPD_L);
N_lines = length(T_L_lines);

% use default parameters from ''default'' structure
P_atm = default.P_atm * ones(1, N_x_axis);
c_a = default.c_a * ones(1, N_x_axis);
RWC_t = default.RWC_t * ones(1, N_x_axis);
pi_L_0_25C = default.pi_L_0_25C * ones(1, N_x_axis);

figure
hold on
for i = 1:N_lines
    
    % leaf temperature for line
    T_L = T_L_lines(i) * ones(1, N_x_axis);
    
    % local temperature-dependent photosynthetic parameters
    [outputs_photo_param] = Photosynthesis_parameters_temperature_response(T_L, 0, V_cmax25_func);
    Gamma_star = outputs_photo_param.Gamma_star_vect;
    K_m = outputs_photo_param.K_m_vect;
    V_cmax = outputs_photo_param.V_cmax_vect;
    R_d = outputs_photo_param.R_d_vect;
    
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
    
    % predict stomatal conductance
    [outputs_stomata] = Stomata_for_RWC_and_Temperature(SWC_L, alpha, RWC_t, RWC_s, pi_L, pi_L_0, ...
                                                        a_f, da_f_dpi_L, a_f_0, da_f_0dpi_L_0, epsilon_L_0, ...
                                                        VPD_L, P_atm, T_L, V_cmax, ...
                                                        R_d, Gamma_star, K_m, c_a);


    g_w = outputs_stomata.g_w;
    g_w(outputs_stomata.solved == 0) = nan;
    
    % calculate m = -dg_w/dln(VPD_L)/g_w_ref from negatively-sloped, primary VPD_L-g_w response
    ind_VPD_sign_change = find(sign(diff(VPD_L)) == -1, 1, 'first') - 1; %where VPD begins to decrease along vector
    if isempty(ind_VPD_sign_change)
        ind_VPD_sign_change = N_x_axis;
    end
    g_w_polyfit = g_w(1:ind_VPD_sign_change); %only consider primary VPD_L-g_w response
    VPD_L_polyfit = VPD_L(1:ind_VPD_sign_change); %only consider primary VPD_L-g_w response
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
    
    P1(i) = plot(VPD_L, g_w, '.-', 'linewidth', 2, 'Color', colors(i,:), 'DisplayName', ['{\itT_L} = ', num2str(T_L_lines(i)), ' ', char(176), 'C, {\itm} = ', num2str(m, '%.2f')]);
    
end
hold off
xlabel('{\itD_L} [kPa]')
ylabel('{\itg_w} [mol\cdotm^{-2}\cdots^{-1}]')
legend(P1(N_lines:-1:1), 'Location', 'NE')
legend boxoff

Output_file = [species_subset, '_Fig_VPD_inst_response_gw_best.jpg'];
print(gcf, '-djpeg', '-r300', Output_file, '-painters', '-noui' )

pause(0)
%% Initializing properties of steady-state plots
g_w_max_plot = 0;
T_L_minus_T_a_max_plot = 0;
T_L_minus_T_a_min_plot = inf;
E_max_plot = 0;

%% Figure: Steady-state response between g_w, T_L, and VPD, varying VPD by changing T_a
RH_lines = 0.3:0.2:0.7;
N_lines = length(RH_lines);
N_x_axis = 50; %2e2;
T_a_min = 0;
T_a_max = 50;
T_a = linspace(T_a_min, T_a_max, N_x_axis);

% use default parameters from ''default'' structure
P_atm = default.P_atm * ones(1, N_x_axis);
c_a = default.c_a * ones(1, N_x_axis);
RWC_t = default.RWC_t * ones(1, N_x_axis);
pi_L_0_25C = default.pi_L_0_25C * ones(1, N_x_axis);

E_vect = 0:1e-4:3e-2; %0:1e-5:3e-2; %initial vector for transpiration
N_vect = length(E_vect);

progress_disp_inc = 0.02; 
progress_disp = progress_disp_inc;

set(0,'units','pixels')
Screen_size = get(0,'screensize');
Fig_size = Screen_size;
Fig_size(3) = 0.8*Fig_size(3);
Fig_size(4) = 0.9*Fig_size(4);

figure
set(gcf, 'Position', Fig_size)
tiledlayout(3,4);
for h = 1:N_lines
    
    RH = RH_lines(h);
    
    % determine g_w for the given T_a and RH
    g_w = nan(1, N_x_axis);
    VPD_L = nan(1, N_x_axis);
    T_L = nan(1, N_x_axis);
    E = nan(1, N_x_axis);
    RH_vect = RH * ones(1, N_vect);
    
    for i = 1:N_x_axis

        progress = (h-1)/N_lines + i/N_x_axis/N_lines;
        if progress >= progress_disp
            disp(['        Calculating steady-state g_w for varied T_a: ', num2str(100*progress_disp_inc*floor(progress/progress_disp_inc), '%.0f'), '% done!'])
            progress_disp = progress_disp_inc * (1 + floor(progress/progress_disp_inc));
        end

        % determine T_L, VPD_L, and g_w from steady-state thermodynamics
        T_a_vect = T_a(i) * ones(1, N_vect);
        P_atm_vect = P_atm(i) * ones(1, N_vect);
        c_a_vect = c_a(i) * ones(1, N_vect);
        RWC_t_vect = RWC_t(i) * ones(1, N_vect);
        pi_L_0_25C_vect = pi_L_0_25C(i) * ones(1, N_vect);
        [T_L_vect, VPD_L_vect, g_w_vect_1] = Steady_state_thermo( T_a_vect, E_vect, RH_vect, P_atm_vect, R_abs, g_H_a, emiss_L ); 
        g_w_vect_1(g_w_vect_1 < 0) = nan;
        
        % local temperature-dependent photosynthetic parameters
        [outputs_photo_param] = Photosynthesis_parameters_temperature_response(T_L_vect, 0, V_cmax25_func);
        Gamma_star_vect = outputs_photo_param.Gamma_star_vect;
        K_m_vect = outputs_photo_param.K_m_vect;
        V_cmax_vect = outputs_photo_param.V_cmax_vect;
        R_d_vect = outputs_photo_param.R_d_vect;

        % predict pressure-volume traits from previously estimated T_L
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

        % determine secondary estimate of g_w from previously estimated T_L and VPD_L
        [outputs_stomata] = Stomata_for_RWC_and_Temperature(SWC_L, alpha, RWC_t_vect, RWC_s_vect, pi_L_vect, pi_L_0_vect, ...
                                                            a_f_vect, da_f_dpi_L_vect, a_f_0_vect, da_f_0dpi_L_0_vect, epsilon_L_0_vect, ...
                                                            VPD_L_vect, P_atm_vect, T_L_vect, V_cmax_vect, ...
                                                            R_d_vect, Gamma_star_vect, K_m_vect, c_a_vect);


        g_w_vect_2 = outputs_stomata.g_w;
        g_w_vect_2(outputs_stomata.solved == 0) = nan;
        
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
            g_w(i) = g_w_vect_1(ind_LB) + (g_w_vect_1(ind_UB) - g_w_vect_1(ind_LB)) * (0 - g_w_error_LB)/(g_w_error_UB - g_w_error_LB);
            VPD_L(i) = VPD_L_vect(ind_LB) + (VPD_L_vect(ind_UB) - VPD_L_vect(ind_LB)) * (0 - g_w_error_LB)/(g_w_error_UB - g_w_error_LB);
            T_L(i) = T_L_vect(ind_LB) + (T_L_vect(ind_UB) - T_L_vect(ind_LB)) * (0 - g_w_error_LB)/(g_w_error_UB - g_w_error_LB);
            E(i) = E_vect(ind_LB) + (E_vect(ind_UB) - E_vect(ind_LB)) * (0 - g_w_error_LB)/(g_w_error_UB - g_w_error_LB);
        end

    end
    
    % calculate m = -dg_w/dln(VPD_L)/g_w_ref from negatively-sloped, primary VPD_L-g_w response
    ind_VPD_sign_change = find(sign(diff(VPD_L)) == -1, 1, 'first') - 1; %where VPD begins to decrease along vector
    if isempty(ind_VPD_sign_change)
        ind_VPD_sign_change = N_x_axis;
    end
    g_w_polyfit = g_w(1:ind_VPD_sign_change); %only consider primary VPD_L-g_w response
    VPD_L_polyfit = VPD_L(1:ind_VPD_sign_change); %only consider primary VPD_L-g_w response
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
    
    j = [3, 2, 1]; %ad hoc fix for desired colors
    j = j(h);
    
    nexttile(1,[2,2])
    hold on
    P2(h) = plot(VPD_L, g_w, '.-', 'linewidth', 2, 'Color', colors(j,:), 'DisplayName', ['RH = ', num2str(RH), ', {\itm} = ', num2str(m, '%.2f')]);
    hold off
    
    nexttile(9)
    hold on
    plot(T_a, T_L-T_a, '.-', 'linewidth', 2, 'Color', colors(j,:));
    hold off
    
    nexttile(10)
    hold on
    plot(T_a, 1e3*E, '.-', 'linewidth', 2, 'Color', colors(j,:));
    hold off
    
    pause(0)
    
end

nexttile(1,[2,2])
g_w_max_plot = max(g_w_max_plot, max(ylim));
title('varied {\itT_a}/constant RH')
xlabel('{\itD_L} [kPa]')
ylabel('{\itg_w} [mol\cdotm^{-2}\cdots^{-1}]')
legend(P2(N_lines:-1:1), 'Location', 'NE')
legend boxoff

nexttile(9)
T_L_minus_T_a_max_plot = max(T_L_minus_T_a_max_plot, max(ylim));
T_L_minus_T_a_min_plot = min(T_L_minus_T_a_min_plot, min(ylim));
xlabel(['{\itT_a} [', char(176), 'C]'])
ylabel(['{\itT_L} - {\itT_a} [', char(176), 'C]'])

nexttile(10)
E_max_plot = max(E_max_plot, max(ylim));
xlabel(['{\itT_a} [', char(176), 'C]'])
ylabel('{\itE} [mmol\cdotm^{-2}\cdots^{-1}]')

pause(0)

%% Figure: Steady-state response between g_w, T_L, and VPD, varying VPD by changing RH
T_a_lines = 15:10:35;
N_lines = length(RH_lines);
N_x_axis = 50; %2e2;
RH_min = 0;
RH_max = 1;
RH = linspace(RH_max, RH_min, N_x_axis); %RH must decline along vector for calculation of m, see below

% use default parameters from ''default'' structure
P_atm = default.P_atm * ones(1, N_x_axis);
c_a = default.c_a * ones(1, N_x_axis);
RWC_t = default.RWC_t * ones(1, N_x_axis);
pi_L_0_25C = default.pi_L_0_25C * ones(1, N_x_axis);

E_vect = 0:1e-4:3e-2; %0:1e-5:3e-2; %initial vector for transpiration
N_vect = length(E_vect);

progress_disp_inc = 0.02; 
progress_disp = progress_disp_inc;

for h = 1:N_lines
    
    T_a = T_a_lines(h);
    
    % determine g_w for the given T_a and RH
    g_w = nan(1, N_x_axis);
    VPD_L = nan(1, N_x_axis);
    T_L = nan(1, N_x_axis);
    E = nan(1, N_x_axis);
    T_a_vect = T_a * ones(1, N_vect);
    
    for i = 1:N_x_axis

        progress = (h-1)/N_lines + i/N_x_axis/N_lines;
        if progress >= progress_disp
            disp(['        Calculating steady-state g_w for varied RH: ', num2str(100*progress_disp_inc*floor(progress/progress_disp_inc), '%.0f'), '% done!'])
            progress_disp = progress_disp_inc * (1 + floor(progress/progress_disp_inc));
        end

        % determine T_L, VPD_L, and g_w from steady-state thermodynamics
        RH_vect = RH(i) * ones(1, N_vect);
        P_atm_vect = P_atm(i) * ones(1, N_vect);
        c_a_vect = c_a(i) * ones(1, N_vect);
        RWC_t_vect = RWC_t(i) * ones(1, N_vect);
        pi_L_0_25C_vect = pi_L_0_25C(i) * ones(1, N_vect);
        [T_L_vect, VPD_L_vect, g_w_vect_1] = Steady_state_thermo( T_a_vect, E_vect, RH_vect, P_atm_vect, R_abs, g_H_a, emiss_L ); 
        g_w_vect_1(g_w_vect_1 < 0) = nan;
        
        % local temperature-dependent photosynthetic parameters
        [outputs_photo_param] = Photosynthesis_parameters_temperature_response(T_L_vect, 0, V_cmax25_func);
        Gamma_star_vect = outputs_photo_param.Gamma_star_vect;
        K_m_vect = outputs_photo_param.K_m_vect;
        V_cmax_vect = outputs_photo_param.V_cmax_vect;
        R_d_vect = outputs_photo_param.R_d_vect;

        % predict pressure-volume traits from previously estimated T_L
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

        % determine secondary estimate of g_w from previously estimated T_L and VPD_L
        [outputs_stomata] = Stomata_for_RWC_and_Temperature(SWC_L, alpha, RWC_t_vect, RWC_s_vect, pi_L_vect, pi_L_0_vect, ...
                                                            a_f_vect, da_f_dpi_L_vect, a_f_0_vect, da_f_0dpi_L_0_vect, epsilon_L_0_vect, ...
                                                            VPD_L_vect, P_atm_vect, T_L_vect, V_cmax_vect, ...
                                                            R_d_vect, Gamma_star_vect, K_m_vect, c_a_vect);


        g_w_vect_2 = outputs_stomata.g_w;
        g_w_vect_2(outputs_stomata.solved == 0) = nan;
        
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
            g_w(i) = g_w_vect_1(ind_LB) + (g_w_vect_1(ind_UB) - g_w_vect_1(ind_LB)) * (0 - g_w_error_LB)/(g_w_error_UB - g_w_error_LB);
            VPD_L(i) = VPD_L_vect(ind_LB) + (VPD_L_vect(ind_UB) - VPD_L_vect(ind_LB)) * (0 - g_w_error_LB)/(g_w_error_UB - g_w_error_LB);
            T_L(i) = T_L_vect(ind_LB) + (T_L_vect(ind_UB) - T_L_vect(ind_LB)) * (0 - g_w_error_LB)/(g_w_error_UB - g_w_error_LB);
            E(i) = E_vect(ind_LB) + (E_vect(ind_UB) - E_vect(ind_LB)) * (0 - g_w_error_LB)/(g_w_error_UB - g_w_error_LB);
        end

    end

    % calculate m = -dg_w/dln(VPD_L)/g_w_ref from negatively-sloped, primary VPD_L-g_w response
    ind_VPD_sign_change = find(sign(diff(VPD_L)) == -1, 1, 'first') - 1; %where VPD begins to decrease along vector
    if isempty(ind_VPD_sign_change)
        ind_VPD_sign_change = N_x_axis;
    end
    g_w_polyfit = g_w(1:ind_VPD_sign_change); %only consider primary VPD_L-g_w response
    VPD_L_polyfit = VPD_L(1:ind_VPD_sign_change); %only consider primary VPD_L-g_w response
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
    
    nexttile(3,[2,2])
    hold on
    P3(h) = plot(VPD_L, g_w, '.-', 'linewidth', 2, 'Color', colors(h,:), 'DisplayName', ['{T_a} = ', num2str(T_a), ', {\itm} = ', num2str(m, '%.2f')]);
    hold off
    
    nexttile(11)
    hold on
    plot(RH, T_L-T_a, '.-', 'linewidth', 2, 'Color', colors(h,:));
    hold off
    
    nexttile(12)
    hold on
    plot(RH, 1e3*E, '.-', 'linewidth', 2, 'Color', colors(h,:));
    hold off
    
    pause(0)
    
end

nexttile(3,[2,2])
g_w_max_plot = max(g_w_max_plot, max(ylim));
title('varied RH/constant {\itT_a}')
xlabel('{\itD_L} [kPa]')
ylabel('{\itg_w} [mol\cdotm^{-2}\cdots^{-1}]')
legend(P3(N_lines:-1:1), 'Location', 'NE')
legend boxoff

nexttile(11)
T_L_minus_T_a_max_plot = max(T_L_minus_T_a_max_plot, max(ylim));
T_L_minus_T_a_min_plot = min(T_L_minus_T_a_min_plot, min(ylim));
xlabel('RH')
ylabel(['{\itT_L} - {\itT_a} [', char(176), 'C]'])

nexttile(12)
E_max_plot = max(E_max_plot, max(ylim));
xlabel('RH')
ylabel('{\itE} [mmol\cdotm^{-2}\cdots^{-1}]')

%% Additional figure adjustments
% set y-lims for steady-state plot
nexttile(1,[2,2])
ylim([0, g_w_max_plot])
nexttile(3,[2,2])
ylim([0, g_w_max_plot])
nexttile(9)
ylim([T_L_minus_T_a_min_plot, T_L_minus_T_a_max_plot])
xlim([T_a_min, T_a_max])
nexttile(10)
ylim([0, E_max_plot])
xlim([T_a_min, T_a_max])
nexttile(11)
ylim([T_L_minus_T_a_min_plot, T_L_minus_T_a_max_plot])
xlim([RH_min, RH_max])
nexttile(12)
ylim([0, E_max_plot])
xlim([RH_min, RH_max])

% vertical line to separate varied T_a from varied RH
L = annotation('line', 0.495*[1,1],[0.03, 0.97]);
set(L, 'Linewidth', 1);

% save figure
Output_file = [species_subset, '_Fig_VPD_ss_response_gw_best.jpg'];
print(gcf, '-djpeg', '-r300', Output_file, '-painters', '-noui' )


end

