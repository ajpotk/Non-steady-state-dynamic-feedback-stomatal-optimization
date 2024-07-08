function [] = plot_sgn_d2Xdotdgw2(N_best, Trt_numb_subset, SSE_best_total_subset_results, pi_L_0_25C_Trt_combo_results, a_f_max_results, pi_L_star_results, beta_results, epsilon_L_max_results, T_L_subset, psi_L_MD_subset)

[~, ind_best] = sort(SSE_best_total_subset_results);
ind_best = ind_best(1:N_best);

unique_Trt_numb_subset = unique(Trt_numb_subset');
N_Trt = length(unique_Trt_numb_subset);

sgn_d2pi_Ldotdgw2 = cell(1, N_Trt);
sgn_d2P_Ldotdgw2 = cell(1, N_Trt);
sgn_d2psi_Ldotdgw2 = cell(1, N_Trt);

for i = unique_Trt_numb_subset
    
    % divide data by treatment
    T_L_Trt = T_L_subset(Trt_numb_subset == i);
    psi_L_MD_Trt = psi_L_MD_subset(Trt_numb_subset == i);
    N_Trt_numb_subset_local = sum(Trt_numb_subset == i);
    sgn_d2pi_Ldotdgw2_local = nan(N_best, N_Trt_numb_subset_local);
    sgn_d2P_Ldotdgw2_local = nan(N_best, N_Trt_numb_subset_local);
    sgn_d2psi_Ldotdgw2_local = nan(N_best, N_Trt_numb_subset_local);
    
    for j = 1:N_best
        
        ind = ind_best(j);
        a_f_max = a_f_max_results(ind);
        pi_L_star = pi_L_star_results(ind);
        beta = beta_results(ind);
        epsilon_L_max = epsilon_L_max_results(ind);
        pi_L_0_25C = pi_L_0_25C_Trt_combo_results(ind, i);
        
        [outputs_PV] = PV_from_psi_for_Trt(a_f_max, pi_L_star, beta, epsilon_L_max, ...
                                            pi_L_0_25C, T_L_Trt, psi_L_MD_Trt);

        RWC_s = outputs_PV.RWC_s;
        RWC_t = outputs_PV.RWC_t;
        pi_L = outputs_PV.pi_L;
        pi_L_0 = outputs_PV.pi_L_0;
        epsilon_L = outputs_PV.epsilon_L;
        epsilon_L_0 = outputs_PV.epsilon_L_0;
        a_f = outputs_PV.a_f;
        a_f_0 = outputs_PV.a_f_0;
        da_f_dpi_L = outputs_PV.da_f_dpi_L;
        da_f_0dpi_L_0 = outputs_PV.da_f_0dpi_L_0;
        
        % osmotic potential as constraint
        sgn_d2pi_Ldotdgw2_local(j,:) = sign(1 + a_f.*pi_L_0./(RWC_t-a_f).*(da_f_0dpi_L_0./(1-a_f_0)-1./epsilon_L_0)./(1-pi_L_0./epsilon_L_0))';
        
        % turgor pressure as constraint
        sgn_d2P_Ldotdgw2_local(j,:) = -sign(da_f_dpi_L./RWC_s+a_f.*(da_f_0dpi_L_0./(1-a_f_0)-1./epsilon_L_0)./(1-pi_L_0./epsilon_L_0))';
        
        % total water potential as constraint
        sgn_d2psi_Ldotdgw2_local(j,:) = -sign(epsilon_L./RWC_s.*da_f_dpi_L./(1-a_f_0) - 1 + a_f./(1-a_f_0).*(epsilon_L - pi_L)./(epsilon_L_0 - pi_L_0).*(epsilon_L_0.*da_f_0dpi_L_0./(1-a_f_0) - 1));
        
    end
    
    sgn_d2pi_Ldotdgw2{i} = sgn_d2pi_Ldotdgw2_local;
    sgn_d2P_Ldotdgw2{i} = sgn_d2P_Ldotdgw2_local;
    sgn_d2psi_Ldotdgw2{i} = sgn_d2psi_Ldotdgw2_local;
    
end
sgn_d2pi_Ldotdgw2 = cellfun(@(x) reshape(x, 1, numel(x)), sgn_d2pi_Ldotdgw2, 'UniformOutput', 0);
sgn_d2pi_Ldotdgw2 = cell2mat(sgn_d2pi_Ldotdgw2);
sgn_d2P_Ldotdgw2 = cellfun(@(x) reshape(x, 1, numel(x)), sgn_d2P_Ldotdgw2, 'UniformOutput', 0);
sgn_d2P_Ldotdgw2 = cell2mat(sgn_d2P_Ldotdgw2);
sgn_d2psi_Ldotdgw2 = cellfun(@(x) reshape(x, 1, numel(x)), sgn_d2psi_Ldotdgw2, 'UniformOutput', 0);
sgn_d2psi_Ldotdgw2 = cell2mat(sgn_d2psi_Ldotdgw2);

perc_pos_sgn_d2pi_Ldotdgw2 = sum(sgn_d2pi_Ldotdgw2 == 1)/length(sgn_d2pi_Ldotdgw2);
perc_pos_sgn_d2P_Ldotdgw2 = sum(sgn_d2P_Ldotdgw2 == 1)/length(sgn_d2P_Ldotdgw2);
perc_pos_sgn_d2psi_Ldotdgw2 = sum(sgn_d2psi_Ldotdgw2 == 1)/length(sgn_d2psi_Ldotdgw2);

disp([num2str(100*perc_pos_sgn_d2pi_Ldotdgw2), '%, ', ...
      num2str(100*perc_pos_sgn_d2P_Ldotdgw2), '%, ', ...
      num2str(100*perc_pos_sgn_d2psi_Ldotdgw2), '%'])

end

