function [outputs_PV] = PV_from_psi_for_Trt( a_f_max, pi_L_star, beta, epsilon_L_max, ...
                                             pi_L_0_25C, T_L_Trt, psi_L_MD_Trt)

%% Fixed parameters
K = 0.1; %parameter controlling shape of relationship between turgor, P_L, and elastic modulus, epsilon_L[MPa]

%% Local functions
[a_f_func, da_f_dpi_L_func] = Local_a_f_functions(a_f_max, pi_L_star, beta);

%% Calculate predictions

pi_L_0 = pi_L_0_25C .* (273.15+T_L_Trt)/298.15;
a_f_0 = a_f_func(pi_L_0); 
da_f_0dpi_L_0 = da_f_dpi_L_func(pi_L_0);
epsilon_L_0 = epsilon_L_max*(1 - exp(pi_L_0/K));
RWC_s_tlp = 1 + pi_L_0./epsilon_L_0;

pi_L_func = @(x, pi_L_0) pi_L_0./x;
P_L_1_func = @(x, psi_L_MD_subset, pi_L_0) psi_L_MD_subset - pi_L_func(x, pi_L_0);
P_L_2_func = @(x, pi_L_0) K * log(x.^(epsilon_L_max/K) .* (exp(-pi_L_0/K) - 1) + 1);
diff_P_L_func = @(x, psi_L_MD_subset, pi_L_0) P_L_1_func(x, psi_L_MD_subset, pi_L_0) - P_L_2_func(x, pi_L_0);

[M, N] = size(psi_L_MD_Trt);
if min(M,N) ~= 1
    error('ERROR: Conditions must be in vector format!')
end

% solve for symplasm relative water content, RWC_s
RWC_s = nan(M, N);
solved_RWC_s = ones(M, N);
N_test = 5;
for i = 1:max(M,N)
    
    diff_P_L_local_func = @(x) diff_P_L_func(x, psi_L_MD_Trt(i), pi_L_0(i));
    
    % solve RWC_s by bisection
    RWC_s_test_max = 1;
    RWC_s_test_min = 0;
    j = 0;
    while 1
        
        j = j + 1;
        if j > 100
           solved_RWC_s(i) = 0;
           break 
        end
        
        RWC_s_test = linspace(RWC_s_test_min, RWC_s_test_max, N_test);
        diff_P_L_local = diff_P_L_local_func(RWC_s_test);
        
        if all(diff_P_L_local > 0) || all(diff_P_L_local < 0)
            solved_RWC_s(i) = 0;
            break 
        end
        
        if any(abs(diff_P_L_local) < 1e-3)
           break 
        end
        
        diff_sign_diff_P_L_local = diff(sign(diff_P_L_local));
        ind_best = find(diff_sign_diff_P_L_local, 1, 'last');
        RWC_s_test_best = RWC_s_test(ind_best + [0,1]);
        RWC_s_test_min = min(RWC_s_test_best);
        RWC_s_test_max = max(RWC_s_test_best);
        
    end
    
    if solved_RWC_s(i)
        [~, ind_best] = min(abs(diff_P_L_local));
        RWC_s(i) = RWC_s_test(ind_best);
    end
    
    
end

P_L = P_L_2_func(RWC_s, pi_L_0);
P_L(P_L < 0) = 0;
pi_L = pi_L_func(RWC_s, pi_L_0);
epsilon_L = epsilon_L_max*(1 - exp(-P_L/K));
a_f = a_f_func(pi_L);
da_f_dpi_L = da_f_dpi_L_func(pi_L);
RWC_t = (1-a_f_0).*RWC_s + a_f; 

%% Store outputs
outputs_PV.RWC_s = RWC_s;
outputs_PV.RWC_s_tlp = RWC_s_tlp;
outputs_PV.RWC_t = RWC_t;
outputs_PV.pi_L = pi_L;
outputs_PV.pi_L_0 = pi_L_0;
outputs_PV.P_L = P_L;
outputs_PV.epsilon_L = epsilon_L;
outputs_PV.epsilon_L_0 = epsilon_L_0;
outputs_PV.a_f = a_f;
outputs_PV.a_f_0 = a_f_0;
outputs_PV.da_f_dpi_L = da_f_dpi_L;
outputs_PV.da_f_0dpi_L_0 = da_f_0dpi_L_0;

end

