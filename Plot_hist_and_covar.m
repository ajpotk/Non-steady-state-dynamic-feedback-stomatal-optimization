function [] = Plot_hist_and_covar(a_f_max_results, pi_L_star_results, beta_results, epsilon_L_max_results, SWC_L_results, alpha_results, ...
                                  pi_L_0_25C_Trt_combo_results, SSE_best_total_subset_results, N_best, ...
                                  xmin, xmax, Nx, species_subset)

%% Density plots of ''N_best'' best parameters and determine most frequently found parameters
N_bins = 20;
color = [102, 102, 255]/255;
axis_label = {'{\ita_{f,max}}'; ...
              '{\it\pi_L^*}';
              '{\it\beta}'; ...
              '{{\it\epsilon_{L,max}}}'; ...
              'SWC_{\itL}'; ...
              '{\it\alpha}'};

[~, ~, ~, a_f_max_results_sub, pi_L_star_results_sub, beta_results_sub, epsilon_L_max_results_sub, SWC_L_results_sub, alpha_results_sub, ~, ~, ~, ~, ~, ~] = Choose_results(pi_L_0_25C_Trt_combo_results, SSE_best_total_subset_results, N_best, species_subset, a_f_max_results, pi_L_star_results, beta_results, epsilon_L_max_results, SWC_L_results, alpha_results);
results = [a_f_max_results_sub, pi_L_star_results_sub, beta_results_sub, epsilon_L_max_results_sub, SWC_L_results_sub, alpha_results_sub];


figure
tiledlayout(Nx, Nx, 'TileSpacing', 'tight')
% plot histograms on-diagonal
for j = 1:Nx
   %%%[xp,xx] = ksdensity(chain_kzs(end-N_best:end,j));
   %%%nexttile(j+(j-1)*Nx)
   %%%plot(xx, xp, 'b-', 'linewidth', 2);

   xbins = xmin(j) + (xmax(j) - xmin(j))*(0:N_bins)/N_bins;

   nexttile(j+(j-1)*Nx)
   h = histogram(results(:,j), 'BinEdges', xbins, 'EdgeColor', 'w', 'FaceColor', color, 'FaceAlpha', 1, 'LineWidth', 1);
   ylim([0, max(h.Values)])
   xlim([xmin(j), xmax(j)])
end

% plot covariation of best parameter estimates off-diagonal
for j = 1:Nx
   for k = 1:Nx
       if j < k
            
           mdl = fitlm(results(:,j), results(:,k));
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
           r = corrcoef(results(:,j), results(:,k));

           nexttile(j+(k-1)*Nx)
           Plot_density(results(:,j), results(:,k), 10, 10, 'xlim', [xmin(j), xmax(j)], 'ylim', [xmin(k), xmax(k)], 'color', color, 'linewidth', 0.5)
           
           title(['{\itr} = ', num2str(r(1,2), 2), stars])

           ylim([xmin(k), xmax(k)])
           xlim([xmin(j), xmax(j)])

       end
   end
end

% edit axes and ticks
for j = 1:Nx
   for k = 1:Nx
       if j <= k
           
           xtick_mod = 1;
           while 1
               if all([xmin(j), xmax(j)] >= 0)
                   dxtick = xtick_mod * 10^floor(log10(xmax(j)));
                   xtick = unique([xmin(j), fliplr(xmax(j):-dxtick:xmin(j))]);
                   xtick(xtick < xmin(j)) = [];
               elseif all([xmin(j), xmax(j)] <= 0)
                   dxtick = xtick_mod * 10^floor(log10(-xmin(j)));
                   xtick = unique([xmin(j):dxtick:xmax(j), xmax(j)]);
                   xtick(xtick > xmax(j)) = [];
               else 
                   error('Error: unplanned possibility!')
               end
               if length(xtick) <= 5
                  break 
               end
               xtick_mod = 2*xtick_mod;
           end
           
           ytick_mod = 1;
           while 1
               if all([xmin(k), xmax(k)] >= 0)
                   dytick = ytick_mod * 10^floor(log10(xmax(k)));
                   ytick = unique([xmin(k), fliplr(xmax(k):-dytick:xmin(k))]);
                   ytick(ytick < xmin(k)) = [];
               elseif all([xmin(k), xmax(k)] <= 0)
                   dytick = ytick_mod * 10^floor(log10(-xmin(k)));
                   ytick = unique([xmin(k):dytick:xmax(k), xmax(k)]);
                   ytick(ytick > xmax(k)) = [];
               else
                   error('Error: unplanned possibility!') 
               end
               if length(ytick) <= 5
                   break
               end
               ytick_mod = 2*ytick_mod;
           end
           
           nexttile(j+(k-1)*Nx)
           box off

           if (j == 1) && (k > 1)
               ylabel(axis_label{k})
           end

           if k == Nx
               xlabel(axis_label{j})
           end

           if j == k
               % remove 'tick' for histograms
              set(gca, 'ytick', '') 
           elseif (j < k) && (j > 1)
               % remove 'yticklabels' for covariation plots 
               set(gca, 'yticklabel', '')
           else
               set(gca, 'ytick', ytick)
           end
           
           set(gca, 'xtick', xtick)
           
           if k < Nx
               % remove 'xticklabels'
               set(gca, 'xticklabel', '')
           end

       end
   end
end

print(gcf, '-djpeg', '-r300', [species_subset, '_Fig_Hist_Covar_for_', num2str(N_best), '_best.jpg'], '-painters', '-noui' )

end

