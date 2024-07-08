function [] = Plot_g_w_layered(a_f_max_results, pi_L_star_results, beta_results, epsilon_L_max_results, SWC_L_results, alpha_results, pi_L_0_25C_Trt_combo_results, SSE_best_total_subset_results, N_best, species_subset, Trt_numb_subset, Comb_Trt_name_subset, g_w_subset, N_species_unique)

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

%% Plot with multiple species (one species added each time ''Plot_g_w.m'' is run) with different colors for different species
All_spec_fig_save_name = 'ALL_SPECIES_Fig_gw_best_layered';
spec_ind_file_save_name = 'ALL_SPECIES-index.mat';
max_DisplayNames = N_species_unique;

% see if file is already present
files = struct2cell(dir);
files = files(1,:);
is_All_spec_fig = cell2mat(cellfun(@(x) strcmp(x, [All_spec_fig_save_name, '.fig']), files, 'UniformOutput', 0));
is_spec_ind_file = cell2mat(cellfun(@(x) strcmp(x, spec_ind_file_save_name), files, 'UniformOutput', 0));

if any(is_All_spec_fig) && any(is_spec_ind_file)
    % open past figure if present
    openfig(All_spec_fig_save_name);
    load(spec_ind_file_save_name, 'spec_ind', 'DisplayName_saved');
    spec_ind = spec_ind + 1; %update the species-index (for plot marker colors)
    if N_colors < spec_ind
        error('Error: More colors must be set for plotting!')
    end
    if max_DisplayNames < spec_ind
        error('Error: ''max_DisplayNames'' must be increased for plotting more species!')
    end
else
    % create new figure if not present
    figure;
    spec_ind = 1;
    DisplayName_saved = cell(1, max_DisplayNames);
    hold on
    plot([axis_min, axis_max], [axis_min, axis_max], 'k--', 'linewidth', 2)
    hold off
end

% determine statistics across all treatments 
mdl = fitlm(g_w_subset(output_subset.solved == 1), output_subset.g_w(output_subset.solved == 1));
try
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
DisplayName = [species_subset_plot, ', {\itr}^2 = ', num2str(r(1,2)^2, 2), stars, ', RMSE = ', num2str(RMSE, 2), ' mol\cdotm^{-2}\cdots^{-1}'];
DisplayName_saved{spec_ind} = DisplayName;

hold on
set_ind = 0;
for i = 1:N_Trt_subset_combo
    
    % plot by treatment with different markers for each treatment
    ind_plot_local = (output_subset.solved == 1) .* (Trt_numb_subset == i);
    if sum(ind_plot_local) > 0
        set_ind = set_ind + 1;
        plot(g_w_subset(ind_plot_local == 1), output_subset.g_w(ind_plot_local == 1), markers{set_ind}, 'color', colors(spec_ind,:), 'markerfacecolor', colors(spec_ind,:))
    end
    
end
plot_DisplayName = gobjects(1, spec_ind);
for i = 1:spec_ind
    plot_DisplayName(i) = plot(-1, -1, markers{1}, 'color', colors(i,:), 'markerfacecolor', colors(i,:), 'DisplayName', DisplayName_saved{i}); %arbitrary location outside of axis limits set below
end
hold off

legend(plot_DisplayName, 'location', 'NW')

xlim([axis_min, axis_max])
ylim([axis_min, axis_max])
xlabel('Observed {\itg_w} [mol\cdotm^{-2}\cdots^{-1}]')
ylabel('Predicted {\itg_w} [mol\cdotm^{-2}\cdots^{-1}]')

P_fig_legend = 25; %percent of y-axis taken up by letters denoting groupings
axis_max_new = axis_max/(1 - P_fig_legend/100);
ytick = get(gca, 'ytick');
dytick = diff(ytick);
dytick = dytick(1);
ytick = min(ytick):dytick:axis_max_new;
ylim([axis_min, axis_max_new])
set(gca, 'ytick', ytick);

% save figure, species-index (for plot marker colors), and results
save(spec_ind_file_save_name, 'spec_ind', 'DisplayName_saved')
savefig(All_spec_fig_save_name)
print(gcf, '-djpeg', '-r300', [All_spec_fig_save_name, '.jpg'], '-painters', '-noui' )
close(gcf)

end

