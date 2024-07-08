function [] = Plot_g_w_environment_response(a_f_max_results, pi_L_star_results, beta_results, epsilon_L_max_results, ...
                                            SWC_L_results, alpha_results, pi_L_0_25C_Trt_combo_results, SSE_best_total_subset_results, ...
                                            N_best, species_subset, T_L_subset, V_cmax_subset)

%% Constants
o_i = 1e-3 * 210; %standard atmospheric partial pressure of oxygen [mol mol-1]
R = 8.314; %universal gas constant [J mol-1 K-1]
R_d25_per_V_cmax25 = 0.015; %from Collatz et al. (1991)

%% Parameter estimates with lowest RMSE
[~, ind_best_sub, ~, ~, ~, ~, ~, ~, ~, a_f_max, pi_L_star, beta, epsilon_L_max, SWC_L, alpha] = Choose_results(pi_L_0_25C_Trt_combo_results, SSE_best_total_subset_results, N_best, species_subset, a_f_max_results, pi_L_star_results, beta_results, epsilon_L_max_results, SWC_L_results, alpha_results);

%% Default environmental conditions
default.P_atm = 101.325; %atmospheric pressure [kPa]
default.VPD_L = 1; %leaf-to-air vapor pressure deficit [kPa]
default.c_a = 4e-4; %atmospheric CO2 concentration [mol mol-1]
default.T_L = 25; %leaf temperature [C]
default.RWC_t = 1; %leaf total relative water content [-]
ind_lowest_SSE = ind_best_sub(1);
default_pi_L_0_25C = mean(pi_L_0_25C_Trt_combo_results(ind_lowest_SSE,:), 'omitnan');
default_pi_L_0_25C_round = 0.1;
default_pi_L_0_25C = default_pi_L_0_25C_round*floor(default_pi_L_0_25C/default_pi_L_0_25C_round);
default.pi_L_0_25C = default_pi_L_0_25C; %osmotic potential at full hydration [MPa] -- taken as average among treatments

%% Photosynthetic parameters -- using mean ''V_cmax25'' and ''R_d25''
V_cmax25_func = @(psi_L_MD) 0;
[outputs_photo_param] = Photosynthesis_parameters_temperature_response(T_L_subset, 0, V_cmax25_func);
V_cmax_per_V_cmax25_subset = outputs_photo_param.V_cmax_per_V_cmax25_vect;
V_cmax25 = mean(V_cmax_subset ./ V_cmax_per_V_cmax25_subset, 'omitnan');
V_cmax25_func = @(psi_L_MD) V_cmax25;

%% Figure colors
colors = [102, 102, 255;...
          124, 203, 161; ...
          240, 116, 110; ...
          253, 222, 156; ...
          002, 144, 153; ...
          220, 057, 120; ...
          006, 082, 117]/255;

%% Figure
default_names = fieldnames(default);
N_default = length(default_names);
x_axis_name = {'RWC_t', 'pi_L_0_25C', 'T_L'};
lines_mod = {[1, 0.95, 0.9], [1.10, 1, 0.90], [35/default.T_L, 1, 15/default.T_L]}; %multiplier
N_response = length(x_axis_name);
x_axis_min = [0.85, (default.pi_L_0_25C - 1), 0]; %same x-axis within a row
x_axis_max = [1, min(0.5, default.pi_L_0_25C + 1), 40]; %same x-axis within a row
N_x_axis = 1e4;
x_axis_label_symbol = {'RWC_{\itt}', ['{\it\pi_L}_{,0,25', char(176), 'C}'], '{\itT_L}'}; %same xlabels within a row
x_axis_label_unit = {'', 'MPa', [char(176), 'C']}; %same xlabels within a row
x_axis_label_unit_bracket = cellfun(@(x) [' [', x, ']'], x_axis_label_unit, 'UniformOutput', 0); %add brackets
x_axis_label_unit_bracket(cell2mat(cellfun(@(x) strcmp('', x), x_axis_label_unit, 'UniformOutput', 0))) = {''}; %remove brackets from originally empty entries
x_axis_label = cellfun(@(x,y) [x, y], x_axis_label_symbol, x_axis_label_unit_bracket, 'UniformOutput', 0);

% information for lines per plot -- not affecting x-axis
lines_mod_names = repmat(x_axis_name, N_response, 1);
lines_mod_names(eye(N_response) == 1) = {[]};
lines_mod = repmat(lines_mod, N_response, 1);
lines_mod_label_symbol = repmat(x_axis_label_symbol, N_response, 1);
lines_mod_label_unit = repmat(x_axis_label_unit, N_response, 1);
for i = 1:N_response %row
    j = 0; %column
    while j < N_response
        j = j + 1; %column
        if isempty(lines_mod_names{i,j})
            lines_mod_names(i, j:N_response-1) = lines_mod_names(i, j+1:N_response);
            lines_mod(i, j:N_response-1) = lines_mod(i, j+1:N_response);
            lines_mod_label_symbol(i, j:N_response-1) = lines_mod_label_symbol(i, j+1:N_response);
            lines_mod_label_unit(i, j:N_response-1) = lines_mod_label_unit(i, j+1:N_response);
            break
        end
    end
end
lines_mod_names(:, N_response) = []; %field names for parameter being modified
lines_mod(:, N_response) = []; %multiplier for default value of parameter being modified 
lines_mod_label_symbol(:, N_response) = [];
lines_mod_label_unit(:, N_response) = [];
N_lines_mod = cell2mat(cellfun(@(x) length(x), lines_mod, 'UniformOutput', 0)); %number of lines per plot

% information for upper ylim of plots
g_w_max_include_for_ylim = 1; %maximum g_w considered for determinin the upper ylim in [mol m-2 s-1]
include_for_ylim = [1, 1; ...
                    0, 0; ...
                    1, 1]; %Boolean
size_include_for_ylim = size(include_for_ylim);
if ~all(size_include_for_ylim == [N_response,N_response-1])
   error(['ERROR: ''include_for_ylim'' must be a ', num2str(N_response), ' by ', num2str(N_response-1), ' array!']) 
end
g_w_max = g_w_max_include_for_ylim; %initiailizing upper ylim

% information for discarding excess trend lines
discard_excess = [0, 0; ...
                  1, 1; ...
                  0, 0]; %Boolean
size_discard_excess = size(discard_excess);
if ~all(size_discard_excess == [N_response,N_response-1])
   error(['ERROR: ''discard_excess'' must be a ', num2str(N_response), ' by ', num2str(N_response-1), ' array!']) 
end

set(0,'units','pixels')
Screen_size = get(0,'screensize');
Fig_size = Screen_size;
Fig_size(3) = 0.8*Fig_size(3);
Fig_size(4) = 0.9*Fig_size(4);

figure
set(gcf, 'Position', Fig_size)
TL = tiledlayout(N_response,N_response-1, 'TileSpacing', 'tight');
TL_title = [species_subset, newline, ...
            '{\ita_{f,max}} = ' num2str(a_f_max), ', {\it\pi_L^*} = ', num2str(pi_L_star), ', {\it\beta} = ', num2str(beta), ...
            newline, '{\it\epsilon_{L,max}} = ', num2str(epsilon_L_max), ', {\itSWC_L} = ', num2str(SWC_L), ', {\it\alpha} = ', num2str(alpha)];
title(TL, TL_title)
for i = 1:N_response %row
    for j = 1:(N_response-1) %column

        % tile index
        ind = j + (N_response-1)*(i-1);
        nexttile(ind)
        hold on
        for m = 1:N_lines_mod(i,j)

            clear local
            local.NAN = nan;
            % access default parameters -- set fields in ''local'' equal to the fields in ''default''
            for k = 1:N_default
                local = setfield(local, default_names{k}, getfield(default, default_names{k}));
                % modify ''local'' for specific line
                if strcmp(default_names{k}, lines_mod_names{i,j})
                    local_orig = getfield(default, default_names{k});
                    local_mod = lines_mod{i,j}(m);
                    local_modified = local_orig * local_mod;
                    local = setfield(local, default_names{k}, local_modified); 
                end
            end
            local = rmfield(local, 'NAN');

            % x-axis of local parameters -- determined by ith entry in ''x_axis_name''
            x_axis = linspace(x_axis_min(i), x_axis_max(i), N_x_axis);
            local = setfield(local, x_axis_name{i}, x_axis);

            % use local parameters from ''local'' structure
            P_atm = local.P_atm .* ones(1, N_x_axis);
            VPD_L = local.VPD_L .* ones(1, N_x_axis);
            c_a = local.c_a .* ones(1, N_x_axis);
            T_L = local.T_L .* ones(1, N_x_axis);
            RWC_t = local.RWC_t .* ones(1, N_x_axis);
            pi_L_0_25C = local.pi_L_0_25C .* ones(1, N_x_axis);

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
            RWC_t = outputs_PV.RWC_t;
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
            
            dchi_wdRWC_s = outputs_stomata.dchi_wdRWC_s;
            if any(dchi_wdRWC_s > 0)
                warning(['WARNING: positive dchi_w/dRWC_s in line ', num2str(m), ' for x-axis of ', x_axis_name{i}, ' with ', lines_mod_names{i,j}, ' = ', num2str(local_modified), ' ', lines_mod_label_unit{i,j}, '!!!'])
            end
            
            % discard excess trend if included
            if discard_excess(i,j) && any(~isnan(g_w))
                
                ind_begin = [];
                ind_end = [];
                
                searching_for_begin = 1; %Boolean
                
                for k = 1:N_x_axis
                    if searching_for_begin
                        % searching for beginning of non-nan's
                        if ~isnan(g_w(k))
                            ind_begin = [ind_begin, k];
                            searching_for_begin = 0;
                        end
                    else
                        % searching for end of non-nan's
                        if isnan(g_w(k))
                            ind_end = [ind_end, k-1];
                            searching_for_begin = 1;
                        elseif (k == N_x_axis) && (~isnan(g_w(k)))
                            ind_end = [ind_end, k];
                        end
                    end
                end
                
                if length(ind_begin) ~= length(ind_end)
                    error('ERROR: Unresolved code preventing equality of lengths of ''ind_begin'' and ''ind_end''!!!')
                end
                
                ind_length = ind_end+1-ind_begin;
                ind_begin = ind_begin(ind_length == max(ind_length));
                ind_end = ind_end(ind_length == max(ind_length));
                ind_g_w = 1:N_x_axis;
                g_w(ind_g_w < ind_begin) = nan;
                g_w(ind_g_w > ind_end) = nan;
                
            end
            
            plot(x_axis, g_w, '.-', 'linewidth', 2, 'Color', colors(m,:), 'DisplayName', [lines_mod_label_symbol{i,j}, ' = ', num2str(local_modified), ' ', lines_mod_label_unit{i,j}])

        end
        hold off

        % plot features
        ylabel('{\itg_w} [mol\cdotm^{-2}\cdots^{-1}]')
        xlabel(x_axis_label{i})
        xlim([x_axis_min(i), x_axis_max(i)])
        legend('Location', 'eastoutside')
        legend boxoff
        
        % find upper ylim of plot if included
        if include_for_ylim(i,j) && (max(ylim) <= g_w_max_include_for_ylim)
            g_w_max = max(g_w_max, max(ylim));
        end
        
    end
end

% apply upper ylim
if any(include_for_ylim)
    for i = 1:N_response %row
        for j = 1:(N_response-1) %column

            % tile index
            ind = j + (N_response-1)*(i-1);
            nexttile(ind)
            ylim([0, g_w_max])

        end
    end
end

Output_file = [species_subset, '_Fig_enviro_response_gw_best.jpg'];
print(gcf, '-djpeg', '-r300', Output_file, '-painters', '-noui' )


end

