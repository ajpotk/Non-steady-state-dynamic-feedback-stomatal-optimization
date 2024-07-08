function [] = Plot_g_w(a_f_max_results, pi_L_star_results, beta_results, epsilon_L_max_results, SWC_L_results, alpha_results, pi_L_0_25C_Trt_combo_results, SSE_best_total_subset_results, N_best, species_subset, Trt_numb_subset, Comb_Trt_name_subset, g_w_subset)

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

%% Plot for just the single species with different colors for different treatments
figure
hold on
plot([axis_min, axis_max], [axis_min, axis_max], 'k--', 'linewidth', 2)
set_ind = 0;
for i = ind_order
    ind_plot_local = (output_subset.solved == 1) .* (Trt_numb_subset == i);
    if sum(ind_plot_local) > 0
        set_ind = set_ind + 1;
        
        mdl = fitlm(g_w_subset(ind_plot_local == 1), output_subset.g_w(ind_plot_local == 1));
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
        r = corrcoef(g_w_subset(ind_plot_local == 1), output_subset.g_w(ind_plot_local == 1));
        RMSE = mean((g_w_subset(ind_plot_local == 1) - output_subset.g_w(ind_plot_local == 1)).^2)^0.5;
        DisplayName = [Trt_subset_combo_plot{i}, ', {\itr}^2 = ', num2str(r(1,2)^2, 2), stars, ', RMSE = ', num2str(RMSE, 2), ' mol\cdotm^{-2}\cdots^{-1}'];
        
        plot_set(set_ind) = plot(g_w_subset(ind_plot_local == 1), output_subset.g_w(ind_plot_local == 1), 'o', 'color', colors(set_ind,:), 'markerfacecolor', colors(set_ind,:), 'DisplayName', DisplayName);
        legend(plot_set, 'location', 'NW')
    end
end
hold off

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
TITLE = [species_subset_plot, newline '{\itr}^2 = ', num2str(r(1,2)^2, 2), stars, newline, 'RMSE = ', num2str(RMSE, 2), ' mol\cdotm^{-2}\cdots^{-1}', newline];
title(TITLE)

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

print(gcf, '-djpeg', '-r300', [species_subset, '_Fig_gw_best.jpg'], '-painters', '-noui' )

end

