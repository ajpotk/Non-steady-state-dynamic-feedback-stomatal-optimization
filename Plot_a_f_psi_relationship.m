function [] = Plot_a_f_psi_relationship(a_f_max_results, pi_L_star_results, beta_results, epsilon_L_max_results, SWC_L_results, alpha_results, pi_L_0_25C_Trt_combo_results, ...
                                        psi_L_MD_subset, T_L_subset, ...
                                        SSE_best_total_subset_results, N_best, species_subset, Trt_numb_subset, Comb_Trt_name_subset, xmax, N_species_unique)

P = 0.50; %P-value for amount of data shown for determining error bars
colors = [102, 102, 255;...
          124, 203, 161; ...
          240, 116, 110; ...
          253, 222, 156; ...
          002, 144, 153; ...
          220, 057, 120; ...
          006, 082, 117]/255;
x_lim = [0.5, 4];
y_lim = [0, xmax(1)];

N_Trt_subset_combo = max(Trt_numb_subset);
Trt_subset_combo = cell(1, N_Trt_subset_combo);
for i = 1:N_Trt_subset_combo
    Trt_subset_combo_local = unique(Comb_Trt_name_subset(Trt_numb_subset == i));
    if ~isempty(Trt_subset_combo_local)
        Trt_subset_combo(i) = Trt_subset_combo_local;
    end
end

% specify names of species and treatments for plots
[species_subset_plot, ~] = Specify_species_and_Trt_names(species_subset, Trt_subset_combo);

[~, ind_best_sub, N_best_sub, a_f_max_results_sub, pi_L_star_results_sub, beta_results_sub, epsilon_L_max_results_sub, ~, ~, ~, ~, ~, ~, ~, ~] = Choose_results(pi_L_0_25C_Trt_combo_results, SSE_best_total_subset_results, N_best, species_subset, a_f_max_results, pi_L_star_results, beta_results, epsilon_L_max_results, SWC_L_results, alpha_results);
pi_L_0_25C_results = mean(pi_L_0_25C_Trt_combo_results(ind_best_sub,:), 2, 'omitnan'); %mean among treatments

[psi_L_MD_subset, ind_sort] = sort(psi_L_MD_subset); %important to sort for plotting below
T_L_subset = T_L_subset(ind_sort); 

%% Determine a_f corresponding to the measured psi_L_MD for each set of parameters
N_psi = length(psi_L_MD_subset); 
a_f = nan(N_best_sub, N_psi);

for i = 1:N_best_sub
    
    [outputs_PV] = PV_from_psi_for_Trt( a_f_max_results_sub(i), pi_L_star_results_sub(i), beta_results_sub(i), epsilon_L_max_results_sub(i), ...
                                        pi_L_0_25C_results(i), T_L_subset, psi_L_MD_subset); 
    
    a_f(i,:) = outputs_PV.a_f;
    
end

% bin data by psi_L_MD
dpsi_L_MD_bin = 0.5;
psi_L_MD_bin_max = dpsi_L_MD_bin*ceil(max(psi_L_MD_subset)/dpsi_L_MD_bin);
psi_L_MD_bin_min = dpsi_L_MD_bin*floor(min(psi_L_MD_subset)/dpsi_L_MD_bin);
psi_L_MD_bin_lb = psi_L_MD_bin_min:dpsi_L_MD_bin:(psi_L_MD_bin_max-dpsi_L_MD_bin);
psi_L_MD_bin_ub = (psi_L_MD_bin_min+dpsi_L_MD_bin):dpsi_L_MD_bin:psi_L_MD_bin_max;
N_bin = length(psi_L_MD_bin_lb);

% determine the mean trend and error bars
% % % a_f_mean = median(a_f, 'omitnan'); 
% % % a_f_upper = quantile(a_f, 0.5+P/2);
% % % a_f_lower = quantile(a_f, 0.5-P/2);
a_f_bin_mean = nan(1,N_bin); 
a_f_bin_upper = nan(1,N_bin); 
a_f_bin_lower = nan(1,N_bin); 
for i = 1:N_bin
    
    a_f_local = a_f(:, (psi_L_MD_subset >= psi_L_MD_bin_lb(i)) + (psi_L_MD_subset < psi_L_MD_bin_ub(i)) == 2);
    [m,n] = size(a_f_local);
    a_f_local = reshape(a_f_local, 1, m*n); 
    a_f_bin_mean(i) = median(a_f_local, 'omitnan'); 
    a_f_bin_upper(i) = quantile(a_f_local, 0.5+P/2);
    a_f_bin_lower(i) = quantile(a_f_local, 0.5-P/2);
    
end

psi_L_MD_bin_plot = reshape([psi_L_MD_bin_lb; psi_L_MD_bin_ub], 1, 2*N_bin);
a_f_bin_mean_plot = reshape([a_f_bin_mean; a_f_bin_mean], 1, 2*N_bin);
a_f_bin_upper_plot = reshape([a_f_bin_upper; a_f_bin_upper], 1, 2*N_bin);
a_f_bin_lower_plot = reshape([a_f_bin_lower; a_f_bin_lower], 1, 2*N_bin);

%% Plot with multiple species (one species added each time ''Plot_g_w.m'' is run) with one panel per species
All_spec_fig_save_name4 = 'ALL_SPECIES_Fig_a_f-psi_L';
spec_ind_file_save_name4 = 'ALL_SPECIES-index-4.mat';
max_DisplayNames = N_species_unique;

% see if file is already present
files = struct2cell(dir);
files = files(1,:);
is_All_spec_fig4 = cell2mat(cellfun(@(x) strcmp(x, [All_spec_fig_save_name4, '.fig']), files, 'UniformOutput', 0));
is_spec_ind_file4 = cell2mat(cellfun(@(x) strcmp(x, spec_ind_file_save_name4), files, 'UniformOutput', 0));

if any(is_All_spec_fig4) && any(is_spec_ind_file4)
    load(spec_ind_file_save_name4, 'spec_ind', 'TITLE_saved');
    spec_ind = spec_ind + 1; %update the species-index (for plot marker colors)
    % open past figure if present
    openfig(All_spec_fig_save_name4);
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

Plot_density(repmat(-psi_L_MD_subset', N_best_sub, 1), a_f, 20, 20, 'xlim', x_lim, 'ylim', y_lim, 'color', colors(1,:), 'linewidth', 0.5);
hold on
% % % plot(-psi_L_MD_subset, a_f_bin_upper(ind_sort), ':', 'color', 'k', 'linewidth', 3)
% % % plot(-psi_L_MD_subset, a_f_bin_lower(ind_sort), ':', 'color', 'k', 'linewidth', 3)
% % % plot(-psi_L_MD_subset, a_f_bin_mean(ind_sort), '-', 'color', colors(3,:), 'linewidth', 4)
plot(-psi_L_MD_bin_plot, a_f_bin_upper_plot, ':', 'color', 'k', 'linewidth', 3)
plot(-psi_L_MD_bin_plot, a_f_bin_lower_plot, ':', 'color', 'k', 'linewidth', 3)
plot(-psi_L_MD_bin_plot, a_f_bin_mean_plot, '-', 'color', colors(3,:), 'linewidth', 4)
hold off
xlabel('-{\it\psi_L} [MPa]')
if spec_ind == 1
    ylabel('{\ita_f}')
end
xlim(x_lim)
ylim(y_lim)

y_min = y_lim(1);
y_max = y_lim(2);
TITLE_saved{spec_ind} = species_subset_plot;
rel_alphab_y = 0.95;
rel_alphab_x = 0.06;
for i = 1:spec_ind
    nexttile(i)
    x_max = max(xlim);
    x_min = min(xlim);
    title(TITLE_saved{i}, 'fontsize', 10)
    LETTER = ['(', char(96 + i), ')'];
    text(x_min + rel_alphab_x*(x_max - x_min), y_min + rel_alphab_y*(y_max - y_min), ...
         LETTER, 'fontsize', 10, 'fontweight', 'bold', 'BackgroundColor', 'w')
end

% save figure, species-index (for plot marker colors), and results
save(spec_ind_file_save_name4, 'spec_ind', 'TITLE_saved')
savefig(All_spec_fig_save_name4)
print(gcf, '-djpeg', '-r300', [All_spec_fig_save_name4, '_for_', num2str(N_best), '_best.jpg'], '-painters', '-noui' )
close(gcf)

end

