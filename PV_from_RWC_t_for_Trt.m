function [outputs_PV] = PV_from_RWC_t_for_Trt( a_f_max, pi_L_star, beta, epsilon_L_max, ...
                                               pi_L_0_25C, T_L_Trt, RWC_t_Trt)

%% Fixed parameters
K = 0.1; %parameter controlling shape of relationship between turgor, P_L, and elastic modulus, epsilon_L[MPa]

%% Local functions
[a_f_func, da_f_dpi_L_func] = Local_a_f_functions(a_f_max, pi_L_star, beta);

%% Calculate predictions
pi_L_0 = pi_L_0_25C .* (T_L_Trt + 273.15)/298.15; 
a_f_0 = a_f_func(pi_L_0);
da_f_0dpi_L_0 = da_f_dpi_L_func(pi_L_0);
epsilon_L_0 = epsilon_L_max*(1 - exp(pi_L_0/K));
RWC_s_tlp = 1 + pi_L_0./epsilon_L_0;

[M, N] = size(RWC_t_Trt);
if min(M,N) ~= 1
    error('ERROR: Conditions must be in vector format!')
end

% solve for symplasm relative water content, RWC_s
RWC_s = nan(M, N);
solved_RWC_s = ones(M, N);
N_test = 5;
diff_RWC_s_func = @(x, pi_L_0, RWC_t, a_f_0) x - (RWC_t - a_f_func(pi_L_0./x))./(1 - a_f_0);
for i = 1:max(M,N)

    diff_RWC_s_local_func = @(x) diff_RWC_s_func(x, pi_L_0(i), RWC_t_Trt(i), a_f_0(i));
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
        diff_RWC_s_local = diff_RWC_s_local_func(RWC_s_test);

        if all(diff_RWC_s_local > 0) || all(diff_RWC_s_local < 0)
            solved_RWC_s(i) = 0;
            break 
        end

        if any(abs(diff_RWC_s_local) < 1e-3)
           break 
        end

        diff_sign_diff_RWC_s_local = diff(sign(diff_RWC_s_local));
        ind_best = find(diff_sign_diff_RWC_s_local, 1, 'last');
        RWC_s_test_best = RWC_s_test(ind_best + [0,1]);
        RWC_s_test_min = min(RWC_s_test_best);
        RWC_s_test_max = max(RWC_s_test_best);

    end

    if solved_RWC_s(i)
        [~, ind_best] = min(abs(diff_RWC_s_local));
        RWC_s(i) = RWC_s_test(ind_best);
    end

end
pi_L = pi_L_0./RWC_s;
a_f = a_f_func(pi_L);
da_f_dpi_L = da_f_dpi_L_func(pi_L);

% determine turgor and total water potential
P_L_2_func = @(x, pi_L_0) K * log(x.^(epsilon_L_max/K) .* (exp(-pi_L_0/K) - 1) + 1);
P_L = P_L_2_func(RWC_s, pi_L_0);
psi_L = P_L + pi_L; 

%% Store outputs
outputs_PV.RWC_s = RWC_s;
outputs_PV.RWC_s_tlp = RWC_s_tlp;
outputs_PV.RWC_t = RWC_t_Trt;
outputs_PV.pi_L = pi_L;
outputs_PV.pi_L_0 = pi_L_0;
outputs_PV.epsilon_L_0 = epsilon_L_0;
outputs_PV.a_f = a_f;
outputs_PV.a_f_0 = a_f_0;
outputs_PV.da_f_dpi_L = da_f_dpi_L;
outputs_PV.da_f_0dpi_L_0 = da_f_0dpi_L_0;
outputs_PV.P_L = P_L;
outputs_PV.psi_L = psi_L;

end

