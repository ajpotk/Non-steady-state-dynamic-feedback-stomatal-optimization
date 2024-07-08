function [] = Fit_Medlyn_USO(A_n_subset, c_a_subset, g_w_subset, VPD_L_subset, T_L_subset, psi_L_MD_subset, V_cmax_subset, R_d_subset, Gamma_star_subset, K_m_subset, species_subset, Trt_numb_subset, Comb_Trt_name_subset, N_species_unique)

include_enviro_factors = 1; %1 for including environmental conditions as fixed effects, and 0 for only treatment as a random effect

colors = [102, 102, 255;...
          124, 203, 161; ...
          240, 116, 110; ...
          253, 222, 156; ...
          002, 144, 153; ...
          220, 057, 120; ...
          006, 082, 117]/255;
markers = {'o'; 's'; 'd'; '^'; 'v'; '+'};

N_Trt_subset_combo = max(Trt_numb_subset);
N_Trt_numb_subset_actual = length(unique(Trt_numb_subset));
[N_colors, ~] = size(colors);
N_markers = length(markers);

if N_colors < N_Trt_numb_subset_actual
    error('Error: More colors must be set for plotting!')
end

if N_markers < N_Trt_numb_subset_actual
    error('Error: More markers must be set for plotting!')
end

Trt_subset_combo = cell(1, N_Trt_subset_combo);
for i = 1:N_Trt_subset_combo
    Trt_subset_combo_local = unique(Comb_Trt_name_subset(Trt_numb_subset == i));
    if ~isempty(Trt_subset_combo_local)
        Trt_subset_combo(i) = Trt_subset_combo_local;
    end
end

% specify names of species and treatments for plots
[species_subset_plot, Trt_subset_combo_plot] = Specify_species_and_Trt_names(species_subset, Trt_subset_combo);

%% Medlyn et al.'s (2011) USO Model
% point estimate of g1
g_c_subset = g_w_subset/1.6;
g1_subset = VPD_L_subset.^0.5 .* (g_c_subset.*c_a_subset./A_n_subset - 1);

if include_enviro_factors == 1
    % mixed effect modeling with environmental conditions as fixed effects and treatment as random effet
    if length(unique(c_a_subset)) > 1
        tbl = table(log(g1_subset), log(c_a_subset), log(T_L_subset + 273.15), psi_L_MD_subset, categorical(Trt_numb_subset), ...
                    'VariableNames', {'log_g1', 'log_c_a', 'log_T_L_K', 'psi_L_MD', 'Trt_numb'});
        lme = fitlme(tbl, 'log_g1 ~ 1 + log_c_a + log_T_L_K + psi_L_MD + (1 | Trt_numb)');
    else
        tbl = table(log(g1_subset), log(T_L_subset + 273.15), psi_L_MD_subset, categorical(Trt_numb_subset), ...
                    'VariableNames', {'log_g1', 'log_T_L_K', 'psi_L_MD', 'Trt_numb'});
        lme = fitlme(tbl, 'log_g1 ~ 1 + log_T_L_K + psi_L_MD + (1 | Trt_numb)');
    end
elseif include_enviro_factors == 0
    % mixed effect modeling with just treatment as random effet
    tbl = table(log(g1_subset), categorical(Trt_numb_subset), ...
                'VariableNames', {'log_g1', 'Trt_numb'});
    lme = fitlme(tbl, 'log_g1 ~ 1 + (1 | Trt_numb)');
else
    error('ERROR: ''include_enviro_factors'' must be either 0 or 1!')
end
disp(lme.Coefficients)
g1_subset_pred = exp(predict(lme));

% predict g_w
c_i_subset_pred = c_a_subset.*g1_subset_pred./(g1_subset_pred + VPD_L_subset.^0.5);
g_w_subset_pred = 1.6./c_a_subset.*(1 + g1_subset_pred./VPD_L_subset.^0.5).*(V_cmax_subset .* (c_i_subset_pred - Gamma_star_subset)./(c_i_subset_pred + K_m_subset) - R_d_subset);


% range of x- and y-axes
daxis = 0.1;
axis_min = min(g_w_subset);
axis_min = daxis * floor(axis_min/daxis); 
axis_max = max(g_w_subset);
axis_max = daxis * ceil(axis_max/daxis);
if isempty(axis_min)
    axis_min = 0;
end
if isempty(axis_max)
    axis_max = 1;
end

%% Plot with multiple species (one species added each time ''Plot_g_w.m'' is run) with one panel per species
All_spec_fig_save_name5 = 'ALL_SPECIES_Fig_gw_Medlyn_USO_panel';
spec_ind_file_save_name5 = 'ALL_SPECIES-index-5.mat';
max_DisplayNames = N_species_unique;

% see if file is already present
files = struct2cell(dir);
files = files(1,:);
is_All_spec_fig5 = cell2mat(cellfun(@(x) strcmp(x, [All_spec_fig_save_name5, '.fig']), files, 'UniformOutput', 0));
is_spec_ind_file5 = cell2mat(cellfun(@(x) strcmp(x, spec_ind_file_save_name5), files, 'UniformOutput', 0));

if any(is_All_spec_fig5) && any(is_spec_ind_file5)
    load(spec_ind_file_save_name5, 'spec_ind', 'TITLE_saved', 'SUBTITLE_saved');
    spec_ind = spec_ind + 1; %update the species-index (for plot marker colors)
    % open past figure if present
    openfig(All_spec_fig_save_name5);
    if N_colors < spec_ind
        error('Error: More colors must be set for plotting!')
    end
    if max_DisplayNames < spec_ind
        error('Error: ''max_DisplayNames'' must be increased for plotting more species!')
    end
else
    spec_ind = 1;
    TITLE_saved = cell(1, max_DisplayNames);
    SUBTITLE_saved = cell(1, max_DisplayNames);
    % create new figure if not present
    set(0,'units','pixels')
    Screen_size = get(0,'screensize');
    Fig_size = Screen_size;
    Fig_size(3) = Fig_size(3);
    Fig_size(4) = 0.55*Fig_size(4);
    figure
    set(gcf, 'Position', Fig_size)
    tiledlayout(1, max_DisplayNames, 'TileSpacing', 'tight');
    for i = 1:max_DisplayNames
        nexttile(i)
    end
end


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


nexttile(spec_ind)
set(gca, 'FontSize', 8)
hold on
plot([0, axis_max], [0, axis_max], 'k--', 'linewidth', 2)
set_ind = 0;
for i = ind_order
    ind_plot_local = (Trt_numb_subset == i);
    if sum(ind_plot_local) > 0
        set_ind = set_ind + 1;
        
        DisplayName = Trt_subset_combo_plot{i};
        
        plot_set(set_ind) = plot(g_w_subset(ind_plot_local == 1), g_w_subset_pred(ind_plot_local == 1), 'o', 'color', colors(set_ind,:), 'markerfacecolor', colors(set_ind,:), 'DisplayName', DisplayName, 'MarkerSize', 3);
        legend(plot_set, 'location', 'southoutside', 'box', 'off')
    end
end
hold off

if axis_max >= 0.4
    daxis_tick = 0.2;
else
    daxis_tick = 0.1;
end
axis_tick = 0:daxis_tick:axis_max;
set(gca, 'xtick', axis_tick, 'ytick', axis_tick)

xlim([0, axis_max])
ylim([0, axis_max])
xlabel('Observed {\itg_w} [mol\cdotm^{-2}\cdots^{-1}]', 'fontsize', 8)
if spec_ind == 1
    ylabel('Predicted {\itg_w} [mol\cdotm^{-2}\cdots^{-1}]', 'fontsize', 8)
end

rel_alphab_y = 0.97;
rel_alphab_x = 0.04;
LETTER = ['(', char(96+spec_ind), ')'];
text(0 + rel_alphab_x*(axis_max - 0), 0 + rel_alphab_y*(axis_max - 0), ...
     LETTER, 'fontsize', 10, 'fontweight', 'bold', 'BackgroundColor', 'w')

 
try
    mdl = fitlm(g_w_subset, g_w_subset_pred);
    p = coefTest(mdl);
    if p < 1e-3
       stars = '***';
    elseif p < 1e-2
       stars = '**';
    elseif p < 5e-2
       stars = '*';
    else
       stars = '';
    end
catch
    stars = '';
end
r = corrcoef(g_w_subset(~isnan(g_w_subset_pred)), g_w_subset_pred(~isnan(g_w_subset_pred)));
RMSE = mean((g_w_subset(~isnan(g_w_subset_pred)) - g_w_subset_pred(~isnan(g_w_subset_pred))).^2)^0.5;
TITLE = species_subset_plot; 
SUBTITLE = ['{\itr}^2 = ', num2str(r(1,2)^2, '%.2f'), stars, newline, 'RMSE = ', num2str(RMSE, '%.3f'), ' mol\cdotm^{-2}\cdots^{-1}', newline];
TITLE_saved{spec_ind} = TITLE;
SUBTITLE_saved{spec_ind} = SUBTITLE;
for i = 1:spec_ind
    nexttile(i)
    title(TITLE_saved{i}, 'fontsize', 10)
    subtitle(SUBTITLE_saved{i}, 'fontsize', 7)
end

% save figure, species-index (for plot marker colors), and results
save(spec_ind_file_save_name5, 'spec_ind', 'TITLE_saved', 'SUBTITLE_saved')
savefig(All_spec_fig_save_name5)
if include_enviro_factors == 0
    extra_text = '_constant_g1_among_Trt';
else
    extra_text = '';
end
% print(gcf, '-djpeg', '-r300', [All_spec_fig_save_name5, extra_text, '.jpg'], '-painters', '-noui' )
close(gcf)

end

