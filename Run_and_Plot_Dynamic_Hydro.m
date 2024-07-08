function [] = Run_and_Plot_Dynamic_Hydro(pi_L_0_25C_Trt_combo_results, SSE_best_total_subset_results, N_best, species_subset, ...
                                         a_f_max_results, pi_L_star_results, beta_results, epsilon_L_max_results, SWC_L_results, alpha_results, ...
                                         T_L_subset, V_cmax_subset, psi_L_MD_subset)

save_file_extension = '_DE-Obe_x2-Vcmax_1.30e-3_k_L_0.2_A_0.80_Z_r_0.172_g_soil'; %'_CH-Dav';

%% PV traits
[~, ind_best_sub, ~, ~, ~, ~, ~, ~, ~, a_f_max, pi_L_star, beta, epsilon_L_max, SWC_L, alpha] = Choose_results(pi_L_0_25C_Trt_combo_results, SSE_best_total_subset_results, N_best, species_subset, a_f_max_results, pi_L_star_results, beta_results, epsilon_L_max_results, SWC_L_results, alpha_results); 
ind_best_sub = ind_best_sub(1);
pi_L_0_25C = mean(pi_L_0_25C_Trt_combo_results(ind_best_sub,:), 'omitnan'); 

%% Photosynthetic traits
mult_V_cmax25 = 2;
V_cmax25_func = Regress_Vcmax25_as_func_psi(T_L_subset, mult_V_cmax25*V_cmax_subset, psi_L_MD_subset, 'Mean_value', 1); 

%% Plant hydraulic conductance traits
k_L_max_25C = 1.30e-3; % maximum soil-plant conductance per unit leaf area at 25C [mol m-2 s-1 MPa-1]
A = 0.2; %exponent for simple soil-plant conductance [MPa-1] -- k_L = k_L_max * exp(A*psi_L)
Q_10_k_L = 1.5; %Q10 value for the temperature-dependence of the maximum soil-plant conductance per unit leaf area [-]

%% Other plant traits
emiss_L = 0.97; %leaf emissivity [-]
Z_r = 0.80; %effective root depth [m] 

%% Soil evaporation properties
g_soil = 0.172; %soil vapor conductance [mol m-2 s-1] -- per ground area basis
Standardized_moisture = 0;

for Sperry_model = 0:1

    %% Save file
    Mat_data_save_file = 'Saved_Dynamic_Hydro_Simulation_data';
    Mat_data_save_file = [Mat_data_save_file, save_file_extension]; 
    if Sperry_model == 1
        Mat_data_save_file = [Mat_data_save_file, '_Sperry'];
        if Standardized_moisture == 1
            Mat_data_save_file = [Mat_data_save_file, '_standardized_moisture'];
        end
    end
    Mat_data_save_file = [Mat_data_save_file, '.mat'];


    files = struct2cell(dir);
    files = files(1,:);
    is_Mat_data_save_file = cell2mat(cellfun(@(x) strcmp(x, Mat_data_save_file), files, 'UniformOutput', 0));

    if any(is_Mat_data_save_file)

        load(Mat_data_save_file)

    else

        %% Time series of environmental conditions
        [Date, P, T_a, RH, P_atm, R_0, g_H_a, c_a, GPP, ET, SPEI, dt] = Download_Flux_data;
        
        %% LAI
        [LAI] = Download_LAI_data(Date);
        A_n_flux = GPP./LAI; %observed canopy scale net carbon assimilation [mol m-2 s-1]
        E_flux = ET./LAI; %observed canopy scale transpiration [mol m-2 s-1]

        %% Soil properties
        sand = 40; %sand content [%] -- taken from https://soilgrids.org/ for location of flux site
        silt = 45; %silt content [%] -- taken from https://soilgrids.org/ for location of flux site
        clay = 15; %clay content [%] -- taken from https://soilgrids.org/ for location of flux site
        rho_soil_bulk = 1.4; %soil bulk density [g cm-3] -- taken from https://soilgrids.org/ for location of flux site
        [psi_soil_func, theta_max, ~, ~] = Pedotransfer_for_psi_as_func_theta(sand, silt, clay, rho_soil_bulk);
        theta_initial = theta_max;

        %% Dynamic simulation
        if (Sperry_model == 1) && (Standardized_moisture == 1)
            
            % load soil moisture from non-Sperry simulations
            Mat_data_save_file_0 = ['Saved_Dynamic_Hydro_Simulation_data', save_file_extension, '.mat']; 
            load(Mat_data_save_file_0, 'theta_store')
            
            [outputs_Simulate_Hydro] = ...
            Simulate_Hydro(SWC_L, alpha, a_f_max, pi_L_star, beta, epsilon_L_max, pi_L_0_25C, ...
                           k_L_max_25C, A, Q_10_k_L, V_cmax25_func, ...
                           theta_initial, theta_max, psi_soil_func, Z_r, LAI, emiss_L, ...
                           Date, P, T_a, RH, P_atm, R_0, g_H_a, c_a, dt, ...
                           'Sperry_model', Sperry_model, 'save_file_extension', save_file_extension, 'g_soil', g_soil, 'theta_store', theta_store); 

        else
            
            [outputs_Simulate_Hydro] = ...
            Simulate_Hydro(SWC_L, alpha, a_f_max, pi_L_star, beta, epsilon_L_max, pi_L_0_25C, ...
                           k_L_max_25C, A, Q_10_k_L, V_cmax25_func, ...
                           theta_initial, theta_max, psi_soil_func, Z_r, LAI, emiss_L, ...
                           Date, P, T_a, RH, P_atm, R_0, g_H_a, c_a, dt, ...
                           'Sperry_model', Sperry_model, 'save_file_extension', save_file_extension, 'g_soil', g_soil); 
            
        end
        
        theta_store = outputs_Simulate_Hydro.theta_store;
        g_w_store = outputs_Simulate_Hydro.g_w_store;
        chi_w_store = outputs_Simulate_Hydro.chi_w_store;
        E_store = outputs_Simulate_Hydro.E_store;
        A_n_store = outputs_Simulate_Hydro.A_n_store;
        c_i_store = outputs_Simulate_Hydro.c_i_store;
        V_cmax_store = outputs_Simulate_Hydro.V_cmax_store;
        J_store = outputs_Simulate_Hydro.J_store;
        J_max_store = outputs_Simulate_Hydro.J_max_store;
        A_c_limited_store = outputs_Simulate_Hydro.A_c_limited_store;
        PLC_store = outputs_Simulate_Hydro.PLC_store;
        psi_L_store = outputs_Simulate_Hydro.psi_L_store;
        psi_L_crit_store = outputs_Simulate_Hydro.psi_L_crit_store;
        E_soil_store = outputs_Simulate_Hydro.E_soil_store;
        
        %% Save data as .mat
        save(Mat_data_save_file, 'Date', 'P', 'T_a', 'RH', 'P_atm', 'R_0', 'g_H_a', 'c_a', 'SPEI', 'LAI', 'A_n_flux', 'E_flux', ...
             'psi_soil_func', 'theta_store', 'g_w_store', 'chi_w_store', 'E_store', 'A_n_store', 'c_i_store', ...
             'V_cmax_store', 'J_store', 'J_max_store', 'A_c_limited_store', 'PLC_store', 'psi_L_store', 'psi_L_crit_store', 'E_soil_store');

    end
    
end

%% Plot results

% colors
color_red = [240, 168, 168]/255;
color_dark_red = [167, 0, 0]/255;
color_blue = [170, 240, 255]/255;
color_dark_blue = [0, 68, 170]/255;
color_green = [161, 207, 171]/255;
color_purple = [184, 122, 184]/255;
color_dark_purple = [102, 0, 102]/255;
color_grey = [0.7, 0.7, 0.7];
color_gold = [0.9, 0.7, 0.1]; 

% load saved file for environmental values 
Sperry_model = 0;
Mat_data_save_file = 'Saved_Dynamic_Hydro_Simulation_data';
Mat_data_save_file = [Mat_data_save_file, save_file_extension]; 
if Sperry_model == 1
    Mat_data_save_file = [Mat_data_save_file, '_Sperry']; 
    if Standardized_moisture == 1
        Mat_data_save_file = [Mat_data_save_file, '_standardized_moisture'];
    end
end
Mat_data_save_file = [Mat_data_save_file, '.mat'];
load(Mat_data_save_file)

% plotted values
is_md = ((hour(Date) == 12) + (minute(Date) == 0) == 2);
Date_md = Date(is_md == 1);
T_a_md = T_a(is_md == 1);
e_a_sat = 0.61078 .* exp(17.27 * T_a ./ (T_a + 237.3)); %air vapor pressure [kPa] -- Teten's equation
VPD = e_a_sat.*(1 - RH); %vapor pressure deficit [kPa]
VPD_md = VPD(is_md == 1);
N_t = length(Date);
P_daily = 1e3*sum(reshape(P, 48, N_t/48)); %daily precipitation [mm]
SPEI_md = SPEI(is_md == 1);
LAI_md = LAI(is_md == 1);

% plot figure
figure
pos = get(gcf, 'position');
pos(3) = 0.70e3;
pos(4) = 0.50e3;
set(gcf, 'position', pos)
TLO = tiledlayout(4, 1, 'TileSpacing', 'none');
pos_TLO = get(TLO, 'InnerPosition');
fract = 6/7;
pos_TLO(3) = pos_TLO(3)*fract;
set(TLO, 'InnerPosition', pos_TLO);
colororder({'k','k'})

for Sperry_model = 0:1
    
    % load saved file
    Mat_data_save_file = 'Saved_Dynamic_Hydro_Simulation_data';
    Mat_data_save_file = [Mat_data_save_file, save_file_extension]; 
    if Sperry_model == 1
        Mat_data_save_file = [Mat_data_save_file, '_Sperry']; 
        if Standardized_moisture == 1
            Mat_data_save_file = [Mat_data_save_file, '_standardized_moisture'];
        end
    end
    Mat_data_save_file = [Mat_data_save_file, '.mat'];
    load(Mat_data_save_file)
    
    % plotted values
    A_n_flux_md = 1e6*A_n_flux(is_md == 1); %[umol m-2 s-1]
    A_n_flux_md(LAI_md < 0.6) = nan;
    A_n_md = 1e6*A_n_store(is_md == 1); %[umol m-2 s-1]
    E_md = 1e3*E_store(is_md == 1); %[mmol m-2 s-1]
    theta_md = theta_store(is_md == 1); 
    PLC_md = 1e2*PLC_store(is_md == 1); %[%]
    
    % lines and color modification
    line = '-';
    linewidth = 1.5;
    color_mod = [1, 1, 1];
    if Sperry_model == 1
        line = '-';
        linewidth = 0.9;
        color_mod = [0.6, 0.6, 0.6];
    end
    
    % plots
    nexttile(1)
    hold on
    plot(Date_md, theta_md, line, 'color', color_mod.*color_green, 'linewidth', linewidth)
    hold off
    
    nexttile(2)
    hold on
    if Sperry_model == 0
        plot(Date_md, A_n_flux_md, line, 'color', 'k', 'linewidth', 0.6)
    end
    plot(Date_md, A_n_md, line, 'color', color_mod.*color_red, 'linewidth', linewidth)
    hold off
    
    nexttile(3)
    hold on
    plot(Date_md, E_md, line, 'color', color_mod.*color_blue, 'linewidth', linewidth)
    hold off
    
    nexttile(4)
    hold on
    plot(Date_md, PLC_md, line, 'color', color_mod.*color_purple, 'linewidth', linewidth)
    hold off
    
    if Sperry_model == 0
        % plot PLC threshold for hydraulic failure for conifers
        hold on
        plot([min(Date_md), max(Date_md)], 50*[1, 1], ':', 'color', color_dark_red, 'linewidth', 1.2)
        plot([min(Date_md), max(Date_md)], 85*[1, 1], ':', 'color', 'b', 'linewidth', 1.2)
        hold off
        text(Date_md(1) + (Date_md(end) - Date_md(1))/6, 50, 'Hydraulic failure threshold', ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Color', color_dark_red)
    end
    
end

% modify plots
ymax_tile = {[1, nan], ...
             [nan], ...
             [nan], ...
             [nan]};
for i = 1:4
    nexttile(i)
    xlim([datetime(2017,1,1), datetime(2019,1,1)])
    ticklength = get(gca, 'ticklength');
    ticklength(1) = ticklength(1)/3;
    set(gca, 'ticklength', ticklength);
    if i < 4
      set(gca, 'xticklabel', []) 
    end
    N_yaxis = length(get(gca, 'yaxis'));
    for j = 1:N_yaxis
        if N_yaxis > 1
            if j == 1
                yyaxis left
            elseif j == 2
                yyaxis right
            end
        end
        ymin = min(0, min(ylim));
        if isnan(ymax_tile{i}(j))
            ymax = max(ylim);
        else
            ymax = ymax_tile{i}(j); 
        end
        N_tick = 5;
        dytick = ymax/(N_tick-1);
        dytick_exponent = floor(log10(dytick));
        dytick_significand = dytick/10^dytick_exponent;
        if ~isinteger(2*dytick_significand)
           dytick = 0.5*ceil(dytick_significand/0.5)*10^dytick_exponent;
           ymax = (N_tick-1)*dytick;
        end
        ylim([ymin, ymax])
        ytick = fliplr(ymax:-dytick:ymin);
        set(gca, 'ytick', ytick);
        yticklabel = num2cell(ytick);
        yticklabel = cellfun(@(x) num2str(x), yticklabel, 'UniformOutput', 0);
        yticklabel{1} = '';
        yticklabel{end} = '';
        set(gca, 'yticklabel', yticklabel);
        box on
        if (N_yaxis > 1) && (j == 2)
            yyaxis left
        end
    end
end

% more plots of supercritical days from new model
Sperry_model = 0;
Mat_data_save_file = 'Saved_Dynamic_Hydro_Simulation_data';
Mat_data_save_file = [Mat_data_save_file, save_file_extension]; 
if Sperry_model == 1
    Mat_data_save_file = [Mat_data_save_file, '_Sperry']; 
    if Standardized_moisture == 1
        Mat_data_save_file = [Mat_data_save_file, '_standardized_moisture'];
    end
end
Mat_data_save_file = [Mat_data_save_file, '.mat'];
load(Mat_data_save_file)

% additional plotting of supercritical days and low SPEI data
is_supercrit = (psi_L_store < psi_L_crit_store);
is_supercrit_daily = any(reshape(is_supercrit', 48, N_t/48));
Date_supercrit = Date_md(is_supercrit_daily == 1);
N_Date_supercrit = length(Date_supercrit);
is_drought = ((SPEI_md < -1) + (year(Date_md) == 2018) == 2);
N_drought = sum(diff(is_drought) == 1);
ind_begin_drought = find(diff(is_drought) == 1);
ind_end_drought = find(diff(is_drought) == -1);
Date_begin_drought = Date_md(ind_begin_drought);
Date_end_drought = Date_md(ind_end_drought);
for i = 1:4
    nexttile(i)
    hold on
    for j = 1:N_drought
        fill([Date_begin_drought(j), Date_end_drought(j), Date_end_drought(j), Date_begin_drought(j)], ...
             [max(ylim)*ones(1,2), min(ylim)*ones(1,2)], color_gold, 'FaceAlpha', 0.2, ...
             'EdgeColor', 'none')
    end
    for j = 1:N_Date_supercrit
        plot([Date_supercrit(j), Date_supercrit(j)], ylim, '-', 'color', color_grey, 'linewidth', 1.5)
    end
    hold off
    % move supercritical lines behind other lines
    h = get(gca, 'Children');
    N_h = length(h);
    set(gca, 'Children', [h(N_Date_supercrit+N_drought+1:N_h); h(1:N_Date_supercrit+N_drought)])
end

% identify 2018 European drought
color_mod = [0.4, 0.4, 0.4];
text(Date_begin_drought(1) + 0.5*(Date_end_drought(1) - Date_begin_drought(1)), 70, ['2018 European', newline, 'drought'], ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', color_mod.*color_gold, 'FontSize', 8)

% ylabels
color_mod = [0.6, 0.6, 0.6];
nexttile(1)
ylabel(['{\color[rgb]{', num2str(color_mod.*color_green) '}{\bf\it\theta}}', newline, '[m^{3}\cdotm^{-3}]'])
nexttile(2)
ylabel(['{\color[rgb]{', num2str(color_mod.*color_red) '}{\bf\itA_n}}', newline, '[\mumol\cdotm^{-2}\cdots^{-1}]'])
nexttile(3)
ylabel(['{\color[rgb]{', num2str(color_mod.*color_blue) '}{\bf\itE}}', newline, '[mmol\cdotm^{-2}\cdots^{-1}]'])
nexttile(4)
ylabel(['{\color[rgb]{', num2str(color_mod.*color_purple) '}\bfPLC}', newline, '[%]'])

% plot of precipitation
line = '-';
linewidth = 1.5;
nexttile(1)
yyaxis right
hold on
for i = 1:length(Date_md)
    plot([Date_md(i), Date_md(i)], [0, -P_daily(i)], line, 'color', color_blue, 'linewidth', linewidth)
end
hold off
ylabel(['{\color[rgb]{', num2str(color_mod.*color_blue) '}{\bf\itP}}', newline, '[mm]'])
ylim([-60, 0])
xlim([datetime(2017,1,1), datetime(2019,1,1)])
ytick = -60:15:0;
set(gca, 'ytick', ytick)
yticklabel = get(gca, 'yticklabel');
yticklabel = cellfun(@(x) num2str(-str2double(x)), yticklabel, 'UniformOutput', 0);
yticklabel{1} = '';
yticklabel{end} = '';
set(gca, 'yticklabel', yticklabel);
yyaxis left

% add letters
rel_alphab_y = 0.93;
rel_alphab_x = 0.015;
for i = 1:4
    nexttile(i)
    xlim_min = min(xlim);
    xlim_max = max(xlim);
    ylim_min = min(ylim);
    ylim_max = max(ylim);
    LETTER = ['(', char(96 + i), ')'];
    text(xlim_min + rel_alphab_x*(xlim_max - xlim_min), ylim_min + rel_alphab_y*(ylim_max - ylim_min), ...
         LETTER, 'fontsize', 9, 'fontweight', 'bold', 'BackgroundColor', 'w', 'VerticalAlignment', 'top')
end

% stars for example subcritical and supercritcial water potentials
is_supercrit = (psi_L_store < psi_L_crit_store);
[~, ind_select_supercrit] = max((psi_L_crit_store - psi_L_store) .* is_supercrit); %most supercritical value chosen as representative value
[~, ind_select_supercrit_md] = min(abs(Date_md - Date(ind_select_supercrit)));
Date_select_supercrit_md = Date_md(ind_select_supercrit_md);
PLC_select_supercrit = 1e2*PLC_store(ind_select_supercrit);
psi_L_md = psi_L_store(is_md == 1);
psi_L_crit_md = psi_L_crit_store(is_md == 1);
PLC_md = 1e2*PLC_store(is_md == 1); %[%]
is_subcrit_and_2017_md = ((psi_L_md > psi_L_crit_md)'+ (Date_md < datetime(2018, 1, 1)) == 2);
[PLC_select_subcrit_md, ind_select_subcrit_md] = max(PLC_md' .* is_subcrit_and_2017_md); %subcritical value with highest PLC chosen as representative value
Date_select_subcrit_md = Date_md(ind_select_subcrit_md);
ind_select_subcrit = find(Date == Date_select_subcrit_md, 1, 'first');
nexttile(4)
hold on
plot(Date_select_supercrit_md, PLC_select_supercrit, 'p', 'MarkerSize', 9, 'LineStyle', 'none', 'MarkerFaceColor', color_dark_red, 'MarkerEdgeColor', color_dark_red)
plot(Date_select_subcrit_md, PLC_select_subcrit_md, 'p', 'MarkerSize', 9, 'LineStyle', 'none', 'MarkerFaceColor', color_dark_blue, 'MarkerEdgeColor', color_dark_blue)
hold off

% plot of temperature
line = '-';
linewidth = 1.5;
nexttile(1)
pos_ax_1 = get(gca, 'position');
pos_ax_1(1) = pos_ax_1(1) + 0.005;
pos_ax_1(3) = pos_ax_1(3)*fract - 0.005;
h_ax_1 = axes('position', pos_ax_1);
set(h_ax_1, 'color', 'none')
hold on
plot(h_ax_1, Date_md, T_a_md, line, 'color', color_red, 'linewidth', linewidth)
hold off
ylim([-20, 30])
xlim([datetime(2017,1,1), datetime(2019,1,1)])
ylim_ax = get(h_ax_1, 'ylim');
ytick_ax = get(h_ax_1, 'ytick');
ticklength = get(h_ax_1, 'ticklength');
ticklength(1) = 0;
set(h_ax_1, 'ticklength', ticklength, 'yticklabel', '')
set(h_ax_1, 'XTick', [], 'XTickLabel', []);
set(h_ax_1, 'YTick', [], 'YTickLabel', []);
set(get(h_ax_1, 'XAxis'), 'Visible', 'off');
set(get(h_ax_1, 'YAxis'), 'Visible', 'off');

% create yaxis and ylabel for axis with temperature
pos_ax_2 = pos_ax_1;
pos_ax_2(1) = pos_ax_2(1) + pos_ax_2(3) + 0.075;
pos_ax_2(3) = 0;
h_ax_2 = axes('position', pos_ax_2);
yyaxis right
set(h_ax_2, 'color', 'none', 'ylim', ylim_ax, 'ytick', ytick_ax)
ylabel(h_ax_2, ['{\color[rgb]{', num2str(color_mod.*color_red) '}{\bf\itT_a}}', newline, '[', char(176), 'C]'])
yyaxis left
set(h_ax_2, 'ytick', [])

% save figure
Output_file = [species_subset, '_flux', save_file_extension, '.jpg'];
if Standardized_moisture == 1
    Output_file = [Output_file, '_standardized_moisture'];
end
Output_file = [Output_file, '.jpg']; 
print(gcf, '-djpeg', '-r300', Output_file, '-painters', '-noui' )



%% second figure to compare to GPP and ET

is_md = ((hour(Date) == 12) + (minute(Date) == 0) == 2);
Date_md = Date(is_md == 1);
GPP_flux_daily = 1800*sum(reshape(LAI.*A_n_flux, 48, length(Date)/48)); %[mol m-2]
ET_flux_daily = 1800*sum(reshape(LAI.*E_flux, 48, length(Date)/48)); %[mol m-2]

% plot figure
figure
tiledlayout(2, 1, 'TileSpacing', 'none')
pos = get(gcf, 'position');
pos(3) = 0.90e3;
pos(4) = 0.30e3;
set(gcf, 'position', pos)

% plot of GPP
nexttile(1)
box on
line = '-';
linewidth = 0.9;
color_mod = [0.6, 0.6, 0.6];
hold on
plot(Date_md, GPP_flux_daily, line, 'color', 'k', 'linewidth', linewidth)
hold off
xlim([datetime(2017,1,1), datetime(2019,1,1)])
ylabel(['{\color[rgb]{', num2str(color_mod.*color_red) '}{\bfDaily sum GPP}}' newline, '[mol\cdotm^{-2}]'])
set(gca, 'xticklabel', []) 

% plot of ET
nexttile(2)
box on
line = '-';
linewidth = 0.9;
color_mod = [0.6, 0.6, 0.6];
hold on
plot(Date_md, ET_flux_daily, line, 'color', 'k', 'linewidth', linewidth)
hold off
xlim([datetime(2017,1,1), datetime(2019,1,1)])
ylabel(['{\color[rgb]{', num2str(color_mod.*color_blue) '}{\bfDaily sum ET}}' newline, '[mol\cdotm^{-2}]'])

for Sperry_model = 0:1
    
    % load data
    Mat_data_save_file = 'Saved_Dynamic_Hydro_Simulation_data';
    Mat_data_save_file = [Mat_data_save_file, save_file_extension]; 
    if Sperry_model == 1
        Mat_data_save_file = [Mat_data_save_file, '_Sperry']; 
        if Standardized_moisture == 1
            Mat_data_save_file = [Mat_data_save_file, '_standardized_moisture'];
        end
    end
    Mat_data_save_file = [Mat_data_save_file, '.mat'];
    load(Mat_data_save_file)
    
    GPP_daily = 1800*sum(reshape(LAI.*A_n_store', 48, length(Date)/48)); %[mol m-2]
    ET_daily = 1800*sum(reshape(LAI.*E_store' + E_soil_store', 48, length(Date)/48)); %[mol m-2]
    
    switch Sperry_model
        case 0
            color_mod = [1, 1, 1];
            linewidth = 1.2;
        case 1
            color_mod = [0.6, 0.6, 0.6];
            linewidth = 0.9;
    end  
    
    nexttile(1)
    hold on
    plot(Date_md, GPP_daily, line, 'color', color_mod.*color_red, 'linewidth', linewidth)
    hold off
    
    nexttile(2)
    hold on
    plot(Date_md, ET_daily, line, 'color', color_mod.*color_blue, 'linewidth', linewidth)
    hold off
    
    
end


% add letters
rel_alphab_y = 0.93;
rel_alphab_x = 0.015;
for i = 1:2
    nexttile(i)
    ticklength = get(gca, 'ticklength');
    ticklength(1) = ticklength(1)/3;
    set(gca, 'ticklength', ticklength);
    xlim_min = min(xlim);
    xlim_max = max(xlim);
    ylim_min = min(ylim);
    ylim_max = max(ylim);
    LETTER = ['(', char(96 + i), ')'];
    text(xlim_min + rel_alphab_x*(xlim_max - xlim_min), ylim_min + rel_alphab_y*(ylim_max - ylim_min), ...
         LETTER, 'fontsize', 9, 'fontweight', 'bold', 'BackgroundColor', 'w', 'VerticalAlignment', 'top')
end

% save figure
Output_file = [species_subset, '_flux_GPP+ET_comparison_to_observations', save_file_extension, '.jpg'];
print(gcf, '-djpeg', '-r300', Output_file, '-painters', '-noui' )


%% third figure of psi_L vs E at subcritical and supercritical periods

% additional plant properties
k_LAI = 0.5; %extinction coefficient for Beer's Law [-]
fR_L_0 = 0.3; %fraction of soil-to_leaf hydraulic resistance in leaves [-] at reference LAI
LAI_0 = 2.5; %reference LAI [m2 m-2]
k_L_max_25C_0 = k_L_max_25C; %maximum soil-plant conductance per unit leaf area at 25C [mol m-2 s-1 MPa-1] at reference LAI
PLC_max = 0.9999;
dpsi_L_vect = 0.01; %[MPa]

% calculate LAI effect on k_L_max_25C
k_L_max_25C = k_L_max_25C_0./(fR_L_0 + (1-fR_L_0)*LAI/LAI_0);

% canopy-averaged absorbed radiation
R_abs = R_0.*(1-exp(-k_LAI*LAI))./LAI;
R_abs(LAI == 0) = 0; 

% plot figure
figure
pos = get(gcf, 'position');
pos(3) = 0.30e3;
pos(4) = 0.50e3;
set(gcf, 'position', pos)
tiledlayout(2, 1, 'TileSpacing', 'none')

% psi_L vs. E curve under subcritical conditions
ind_select = [ind_select_subcrit, ind_select_supercrit];
Title = {'Subcritical period'; 'Supercritical period'};
for Sperry_model = 0:1
    
    % load data
    Mat_data_save_file = 'Saved_Dynamic_Hydro_Simulation_data';
    Mat_data_save_file = [Mat_data_save_file, save_file_extension]; 
    if Sperry_model == 1
        Mat_data_save_file = [Mat_data_save_file, '_Sperry']; 
        if Standardized_moisture == 1
            Mat_data_save_file = [Mat_data_save_file, '_standardized_moisture'];
        end
    end
    Mat_data_save_file = [Mat_data_save_file, '.mat'];
    load(Mat_data_save_file)
    
    line = '-';
    linewidth = 1.5;
    color_mod = [0.6, 0.6, 0.6];
    
    switch Sperry_model
        case 0
            color_mod = [1, 1, 1];
            linewidth = 1.5;
        case 1
            color_mod = [0.6, 0.6, 0.6];
            linewidth = 1;
    end  
    
    for i = 1:2

        % determine discretization
        psi_soil = psi_soil_func(theta_store(ind_select(i)));
        psi_L_max = psi_soil;
        psi_L_max = 0.1*ceil(psi_L_max/0.1); 
        psi_L_min = log(1-PLC_max)/A;
        psi_L_vect = [psi_L_max:-abs(dpsi_L_vect):psi_L_min, -inf];
        PLC_vect = 1 - exp(A*psi_L_vect);

        % determine psi_L vs. E curve
        E_vect = k_L_max_25C(ind_select(i)) * Q_10_k_L.^((T_a(ind_select(i))-25)/10) .* exp(A*psi_L_vect) .* (psi_soil - psi_L_vect);
        E_vect(psi_L_vect == -inf) = 0;
        [E_crit, ind_crit] = max(E_vect);
        psi_L_crit = psi_L_vect(ind_crit);
        PLC_crit = 1 - exp(A*psi_L_crit);
        
        if Sperry_model == 0
        % if not the Sperry model, predict stomatal conductance by Sperry model for psi_soil from my model
            Nonmonotonic_V_cmax_temperature_response = 1;
            Electron_transport_limitation = 1;
            [outputs_stomata_Sperry] = ...
            Stomata_Sperry_Gain_Risk(k_L_max_25C(ind_select(i)), A, Q_10_k_L, psi_soil, ...
                                     T_a(ind_select(i)), RH(ind_select(i)), P_atm(ind_select(i)), c_a(ind_select(i)), ...
                                     R_abs(ind_select(i)), g_H_a(ind_select(i)), emiss_L, ...
                                     V_cmax25_func, ...
                                     'Nonmonotonic_V_cmax_temperature_response', Nonmonotonic_V_cmax_temperature_response, ...
                                     'Electron_transport_limitation', Electron_transport_limitation);
            E_Sperry = outputs_stomata_Sperry.E_Sperry;
            PLC_Sperry = outputs_stomata_Sperry.PLC_Sperry;
        end

        % plot psi_L vs. E curve
        switch Sperry_model
            case 0
                switch i
                    case 1
                        MarkerFaceColor = color_dark_blue;
                    case 2
                        MarkerFaceColor = color_dark_red;
                end   
                MarkerShape = 'p';
                MarkerSize = 14;
            case 1
                MarkerFaceColor = color_mod.*color_purple;
                MarkerShape = 's';
                MarkerSize = 5;
        end
        nexttile(i)
        box on
        hold on
        plot(1e2*PLC_vect, 1e3*E_vect, line, 'color', color_mod.*color_purple, 'linewidth', linewidth)
        plot(1e2*PLC_crit*[1, 1], 1e3*E_crit*[0, 1], 'k-', 'LineWidth', 0.9)
        plot(1e2*PLC_crit, 1e3*E_crit, 'ko', 'MarkerSize', 6, 'LineStyle', 'none', 'MarkerFaceColor', 'w')
        if Sperry_model == 0
            plot(1e2*PLC_Sperry, 1e3*E_Sperry, 's', 'MarkerSize', 7, 'LineStyle', 'none', 'MarkerFaceColor', color_mod.*color_purple, 'MarkerEdgeColor', color_mod.*color_purple)
        end
        plot(1e2*PLC_store(ind_select(i)), 1e3*E_store(ind_select(i)), MarkerShape, 'MarkerSize', MarkerSize, 'LineStyle', 'none', 'MarkerFaceColor', MarkerFaceColor, 'MarkerEdgeColor', MarkerFaceColor)
        hold off
        xlim([0, 100])
        set(gca, 'xtick', 0:25:100)
        ylim([0, 2.25])
        set(gca, 'ytick', 0.5:0.5:2)
        ylabel(['{\color[rgb]{', num2str(color_mod.*color_blue) '}{\bf\itE}}', newline, '[mmol\cdotm^{-2}\cdots^{-1}]'])
        if i == 2
           xlabel(['{\color[rgb]{', num2str(color_mod.*color_purple) '}{\bfPLC}} [%]'])
        elseif i == 1
            set(gca, 'xticklabel', []) 
            text(1e2*PLC_crit, 1e3*E_crit, '{\itE_{crit}}', ...
                'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'Color', 'k', 'FontSize', 9)
            text(1e2*PLC_crit + 0.01*max(xlim), 0.01*max(ylim), 'PLC({\it\psi_{L,crit}})', ...
                'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'Color', 'k', 'FontSize', 9)
            % arrows to denote supercritical and subcritical ranges
            arrow_start_y = 0.4;
            supercrit_arrow_start_x = 1e2*(1 + PLC_crit)/2;
            supercrit_arrow_length_x = 1e2 - supercrit_arrow_start_x;
            subcrit_arrow_start_x = 1e2*PLC_crit/2;
            subcrit_arrow_length_x = subcrit_arrow_start_x;
            hold on
            quiver(supercrit_arrow_start_x, arrow_start_y, supercrit_arrow_length_x, 0, 'color', 'k', 'MaxHeadSize', 0.2/supercrit_arrow_length_x)
            quiver(supercrit_arrow_start_x, arrow_start_y, -supercrit_arrow_length_x, 0, 'color', 'k', 'MaxHeadSize', 0.2/supercrit_arrow_length_x)
            quiver(subcrit_arrow_start_x, arrow_start_y, subcrit_arrow_length_x, 0, 'color', 'k', 'MaxHeadSize', 0.2/subcrit_arrow_length_x)
            quiver(subcrit_arrow_start_x, arrow_start_y, -subcrit_arrow_length_x, 0, 'color', 'k', 'MaxHeadSize', 0.2/subcrit_arrow_length_x)
            hold off
            text(supercrit_arrow_start_x, arrow_start_y, 'Supercritical', ...
                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Color', 'k', 'FontSize', 8)
            text(subcrit_arrow_start_x, arrow_start_y, 'Subcritical', ...
                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Color', 'k', 'FontSize', 8)
        end
        text(50, 2, [Title{i}, newline, num2str(month(Date(ind_select(i)))), '/', num2str(day(Date(ind_select(i)))), '/', num2str(year(Date(ind_select(i))))], ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', 'k', 'FontSize', 12)
    end
end

% add letters
rel_alphab_y = 0.97;
rel_alphab_x = 0.05;
for i = 1:2
    nexttile(i)
    xlim_min = min(xlim);
    xlim_max = max(xlim);
    ylim_min = min(ylim);
    ylim_max = max(ylim);
    LETTER = ['(', char(96 + 4 + i), ')'];
    text(xlim_min + rel_alphab_x*(xlim_max - xlim_min), ylim_min + rel_alphab_y*(ylim_max - ylim_min), ...
         LETTER, 'fontsize', 9, 'fontweight', 'bold', 'BackgroundColor', 'w', 'VerticalAlignment', 'top')
end

% save figure
Output_file = [species_subset, '_flux_E-PLC_supercrit', save_file_extension];
if Standardized_moisture == 1
    Output_file = [Output_file, '_standardized_moisture'];
end
Output_file = [Output_file, '.jpg']; 
print(gcf, '-djpeg', '-r300', Output_file, '-painters', '-noui' )

end

