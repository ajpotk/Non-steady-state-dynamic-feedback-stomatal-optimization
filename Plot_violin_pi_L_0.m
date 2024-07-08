function [] = Plot_violin_pi_L_0(a_f_max_results, pi_L_star_results, beta_results, epsilon_L_max_results, SWC_L_results, alpha_results, pi_L_0_25C_Trt_combo_results, SSE_best_total_subset_results, N_species_unique, N_best, species_subset, Trt_numb_subset, Comb_Trt_name_subset)
%% Variation in pi_L_0_25C across treatments

colors = [102, 102, 255;...
          124, 203, 161; ...
          240, 116, 110; ...
          253, 222, 156; ...
          002, 144, 153; ...
          220, 057, 120; ...
          006, 082, 117]/255;
[N_colors, ~] = size(colors);

N_Trt_subset_combo = max(Trt_numb_subset);
Trt_subset_combo = cell(1, N_Trt_subset_combo);
for i = 1:N_Trt_subset_combo
    Trt_subset_combo_local = unique(Comb_Trt_name_subset(Trt_numb_subset == i));
    if ~isempty(Trt_subset_combo_local)
        Trt_subset_combo(i) = Trt_subset_combo_local;
    end
end

% specify names of species and treatments for plots
[species_subset_plot, Trt_subset_combo_plot] = Specify_species_and_Trt_names(species_subset, Trt_subset_combo);
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
Trt_subset_combo = Trt_subset_combo(ind_order);
Trt_subset_combo_plot = Trt_subset_combo_plot(ind_order);
pi_L_0_25C_Trt_combo_results = pi_L_0_25C_Trt_combo_results(:, ind_order);

is_not_empty = cell2mat(cellfun(@(x) ~isempty(x), Trt_subset_combo, 'UniformOutput', 0));
Trt_subset_combo_not_empty = Trt_subset_combo_plot(is_not_empty == 1);
pi_L_0_25C_Trt_combo_best_not_empty = pi_L_0_25C_Trt_combo_results(:, is_not_empty == 1);

[ind_best, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = Choose_results(pi_L_0_25C_Trt_combo_results, SSE_best_total_subset_results, N_best, species_subset, a_f_max_results, pi_L_star_results, beta_results, epsilon_L_max_results, SWC_L_results, alpha_results);
pi_L_0_25C_Trt_combo_best_not_empty = pi_L_0_25C_Trt_combo_best_not_empty(ind_best(1:N_best), :);

% mean, median, Q1, Q3, and Tukey's fence
Mean = mean(-pi_L_0_25C_Trt_combo_best_not_empty, 'omitnan');
Median = median(-pi_L_0_25C_Trt_combo_best_not_empty, 'omitnan');
Q1 = prctile(-pi_L_0_25C_Trt_combo_best_not_empty, 25);
Q3 = prctile(-pi_L_0_25C_Trt_combo_best_not_empty, 75);
IQR = Q3 - Q1;
Tukey_max = Q3 + 1.5*IQR;
Tukey_min = Q1 - 1.5*IQR;

cat_Trt_subset_combo = categorical(Trt_subset_combo_not_empty);
cat_Trt_subset_combo = reordercats(cat_Trt_subset_combo, Trt_subset_combo_not_empty);
N_Trt_subset_combo = length(cat_Trt_subset_combo);

% initial y-axes range of violin plot
y_max = 4; %%%max(max(Tukey_max), max(-pi_L_0_25C_Trt_combo_best_not_empty, [], 'all'));
y_min = 0; %%%max(min(min(Tukey_min), min(-pi_L_0_25C_Trt_combo_best_not_empty, [], 'all')), 0);
dy = 0.1;
y_max = dy*ceil(y_max/dy);
y_min = dy*floor(y_min/dy);
y = y_min:1e-3:y_max;
N_y = length(y);
ind = 1:N_y;

% determine normalization factor for PDF
max_PDF = 0;
for i = 1:N_Trt_subset_combo
    pd = fitdist(-pi_L_0_25C_Trt_combo_best_not_empty(:,i),'Kernel','Kernel','epanechnikov');
    PDF = pdf(pd, y);
    max_PDF = max(max_PDF, max(PDF));
end

%% Plot with multiple species (one species added each time ''Plot_g_w.m'' is run) with one panel per species
All_spec_fig_save_name3 = 'ALL_SPECIES_Fig_violin_pi_L_0';
spec_ind_file_save_name3 = 'ALL_SPECIES-index-3.mat';
max_DisplayNames = N_species_unique;

% see if file is already present
files = struct2cell(dir);
files = files(1,:);
is_All_spec_fig3 = cell2mat(cellfun(@(x) strcmp(x, [All_spec_fig_save_name3, '.fig']), files, 'UniformOutput', 0));
is_spec_ind_file3 = cell2mat(cellfun(@(x) strcmp(x, spec_ind_file_save_name3), files, 'UniformOutput', 0));

if any(is_All_spec_fig3) && any(is_spec_ind_file3)
    load(spec_ind_file_save_name3, 'spec_ind', 'TITLE_saved');
    spec_ind = spec_ind + 1; %update the species-index (for plot marker colors)
    % open past figure if present
    openfig(All_spec_fig_save_name3);
    if N_colors < spec_ind
        error('Error: More colors must be set for plotting!')
    end
    if max_DisplayNames < spec_ind
        error('Error: ''max_DisplayNames'' must be increased for plotting more species!')
    end
else
    spec_ind = 1;
    TITLE_saved = cell(1, max_DisplayNames);
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

nexttile(spec_ind)
hold on
for i = 1:N_Trt_subset_combo
    
    % Calculate Kernel density distribution
    pd = fitdist(-pi_L_0_25C_Trt_combo_best_not_empty(:,i),'Kernel','Kernel','epanechnikov');
    PDF = 0.8 * pdf(pd, y)/2/max_PDF;
    PDF(PDF == 0) = nan;
    
    % Plot Kernel density distribution
    isnan_PDF = isnan(PDF);
    N = 1; %number of continuous regions in PDF where positive
    while 1
        if N == 1
            ind_search_start = 1;
        else
            ind_search_start = ind_end+1; 
        end
        ind_start = find(((isnan_PDF == 0) + (ind >= ind_search_start) == 2), 1, 'first');
        if isempty(ind_start)
           break 
        end
        ind_end = find(((isnan_PDF == 1) + (ind >= ind_start) == 2), 1, 'first') - 1;
        if isempty(ind_end)
            ind_end = N_y;
        end
        N = N + 1;
        
        PDF_local = PDF(ind_start:ind_end);
        y_local = y(ind_start:ind_end);
        fill(i + [PDF_local, -fliplr(PDF_local)], [y_local, fliplr(y_local)], colors(i,:), 'LineStyle', 'none') 
        
        if ind_end == N_y
            break
        end
    end
    
    % Plot Tukey's fence, IQR, Median, and Mean
    plot(i*[1,1], [Tukey_min(i), Tukey_max(i)], '-', 'LineWidth', 1, 'Color', [0.5, 0.5, 0.7])
    plot(i*[1,1], [Q1(i), Q3(i)], '-', 'LineWidth', 4, 'Color', [0.5, 0.5, 0.7])
    plot(i, Median(i), 'o', 'MarkerFaceColor', 'w', 'LineWidth', 0.75, 'MarkerSize', 6, 'Color', [0.5, 0.5, 0.7])
    plot(i, Mean(i), 'k+', 'MarkerFaceColor', 'k', 'LineWidth', 1.5, 'MarkerSize', 6)
    
end
hold off
if spec_ind == 1
    ylabel('-{\it\pi_L}_{,0,25^{\circ}C} [MPa]')
end
xlim(0.5+ [0, N_Trt_subset_combo])
set(gca, 'xtick', 1:N_Trt_subset_combo)
set(gca, 'xticklabel', cat_Trt_subset_combo)
xtickangle(60)

% ANOVA to determine if groups are different
[grouping] = Group_by_ANOVA(-pi_L_0_25C_Trt_combo_best_not_empty);
P_fig_letters = 0; %percent of y-axis taken up by letters denoting groupings
y_max_new = y_max/(1-P_fig_letters/100);
% % % ytick = get(gca, 'ytick');
% % % dytick = diff(ytick);
% % % dytick = dytick(1);
dytick = 1;
ytick = y_min:dytick:y_max_new;
set(gca, 'ytick', ytick);
ylim([y_min, y_max_new])
for i = 1:N_Trt_subset_combo
    %%%text(i + 0.15, y_max + 0.3*(y_max_new - y_max), ['{\bf', grouping{i}, '}'], 'HorizontalAlignment', 'center')
    text(i + 0.2, y_min + 0.10*(y_max - y_min), ['{\bf', grouping{i}, '}'], 'HorizontalAlignment', 'center')
end

TITLE_saved{spec_ind} = species_subset_plot;
rel_alphab_y = 0.95;
rel_alphab_x = 0.06;
for i = 1:spec_ind
    nexttile(i)
    x_max = max(xlim);
    x_min = min(xlim);
    title(TITLE_saved{i}, 'fontsize', 10)
    LETTER = ['(', char(96 + i), ')'];
    text(x_min + rel_alphab_x*(x_max - x_min), y_min + rel_alphab_y*(y_max_new - y_min), ...
         LETTER, 'fontsize', 10, 'fontweight', 'bold', 'BackgroundColor', 'w')
end

% save figure, species-index (for plot marker colors), and results
save(spec_ind_file_save_name3, 'spec_ind', 'TITLE_saved')
savefig(All_spec_fig_save_name3)
print(gcf, '-djpeg', '-r300', [All_spec_fig_save_name3, '_for_', num2str(N_best), '_best.jpg'], '-painters', '-noui' )
close(gcf)


end

