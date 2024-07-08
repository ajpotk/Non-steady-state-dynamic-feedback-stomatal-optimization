function [] = Plot_g_w_panel(a_f_max_results, pi_L_star_results, beta_results, epsilon_L_max_results, SWC_L_results, alpha_results, pi_L_0_25C_Trt_combo_results, SSE_best_total_subset_results, N_best, species_subset, Trt_numb_subset, Comb_Trt_name_subset, g_w_subset, N_species_unique)

%% Predictions of parameter estimates with lowest RMSE
[~, ~, ~, ~, ~, ~, ~, ~, ~, a_f_max_best_sub, pi_L_star_best_sub, beta_best_sub, epsilon_L_max_best_sub, SWC_L_best_sub, alpha_best_sub] = Choose_results(pi_L_0_25C_Trt_combo_results, SSE_best_total_subset_results, N_best, species_subset, a_f_max_results, pi_L_star_results, beta_results, epsilon_L_max_results, SWC_L_results, alpha_results); 
x = [a_f_max_best_sub, pi_L_star_best_sub, beta_best_sub, epsilon_L_max_best_sub, SWC_L_best_sub, alpha_best_sub]; 
[output_subset] = Stomata_and_PV_best_for_subset(x);

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

% range of x- and y-axes
daxis = 0.1;
axis_min = min([min(g_w_subset(output_subset.solved == 1)), min(output_subset.g_w(output_subset.solved == 1))]);
axis_min = daxis * floor(axis_min/daxis); 
axis_max = max([max(g_w_subset(output_subset.solved == 1)), max(output_subset.g_w(output_subset.solved == 1))]);
axis_max = daxis * ceil(axis_max/daxis);
if isempty(axis_min)
    axis_min = 0;
end
if isempty(axis_max)
    axis_max = 1;
end

%% Plot with multiple species (one species added each time ''Plot_g_w.m'' is run) with one panel per species
All_spec_fig_save_name2 = 'ALL_SPECIES_Fig_gw_best_panel';
spec_ind_file_save_name2 = 'ALL_SPECIES-index-2.mat';
max_DisplayNames = N_species_unique;

% see if file is already present
files = struct2cell(dir);
files = files(1,:);
is_All_spec_fig2 = cell2mat(cellfun(@(x) strcmp(x, [All_spec_fig_save_name2, '.fig']), files, 'UniformOutput', 0));
is_spec_ind_file2 = cell2mat(cellfun(@(x) strcmp(x, spec_ind_file_save_name2), files, 'UniformOutput', 0));

if any(is_All_spec_fig2) && any(is_spec_ind_file2)
    load(spec_ind_file_save_name2, 'spec_ind', 'TITLE_saved', 'SUBTITLE_saved');
    spec_ind = spec_ind + 1; %update the species-index (for plot marker colors)
    % open past figure if present
    openfig(All_spec_fig_save_name2);
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
    ind_plot_local = (output_subset.solved == 1) .* (Trt_numb_subset == i);
    if sum(ind_plot_local) > 0
        set_ind = set_ind + 1;
        
        DisplayName = Trt_subset_combo_plot{i};
        
        plot_set(set_ind) = plot(g_w_subset(ind_plot_local == 1), output_subset.g_w(ind_plot_local == 1), 'o', 'color', colors(set_ind,:), 'markerfacecolor', colors(set_ind,:), 'DisplayName', DisplayName, 'MarkerSize', 3);
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
    mdl = fitlm(g_w_subset(output_subset.solved == 1), output_subset.g_w(output_subset.solved == 1));
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
r = corrcoef(g_w_subset(output_subset.solved == 1), output_subset.g_w(output_subset.solved == 1));
RMSE = mean((g_w_subset(output_subset.solved == 1) - output_subset.g_w(output_subset.solved == 1)).^2)^0.5;
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
save(spec_ind_file_save_name2, 'spec_ind', 'TITLE_saved', 'SUBTITLE_saved')
savefig(All_spec_fig_save_name2)
print(gcf, '-djpeg', '-r300', [All_spec_fig_save_name2, '.jpg'], '-painters', '-noui' )
close(gcf)

end

