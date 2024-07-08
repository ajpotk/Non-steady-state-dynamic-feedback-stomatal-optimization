function [SSE_test] = Estimate_a_f_star_and_beta(RWC_t_subset, psi_L_MD_subset, RWC_t_tlp, pi_L_tlp)
% a_f_star, beta, pi_L_0, a_f_0, epsilon_L_0

a_f_0_test_max = 1;
da_f_0_test = 1e-2;
a_f_0_test = da_f_0_test:da_f_0_test:a_f_0_test_max;
N_a_f_0_test = length(a_f_0_test);

beta_test_max = 3;
dbeta_test = 1e-2;
beta_test = dbeta_test:dbeta_test:beta_test_max;
N_beta_test = length(beta_test);

% remove nan's
is_not_nan = ~isnan(RWC_t_subset) .* ~isnan(psi_L_MD_subset);
RWC_t_subset = RWC_t_subset(is_not_nan == 1);
psi_L_MD_subset = psi_L_MD_subset(is_not_nan == 1);
N_data_subset = sum(is_not_nan); 
N_data_below_tlp = sum(RWC_t_subset <= RWC_t_tlp);
N_data_above_tlp = sum(RWC_t_subset > RWC_t_tlp);
weight_subset = nan(N_data_subset, 1); %weights for calculating sum-square error
weight_subset(RWC_t_subset <= RWC_t_tlp) = 0.5*N_data_subset/N_data_below_tlp;
weight_subset(RWC_t_subset > RWC_t_tlp) = 0.5*N_data_subset/N_data_above_tlp;

diff_RWC_s_func = @(RWC_s, RWC_t, a_f_0, beta) RWC_s - (RWC_t - a_f_0*RWC_s.^beta)/(1-a_f_0);
options = optimset('display', 'off');

SSE_test = nan(N_a_f_0_test, N_beta_test); %sum-square error
pi_L_0_test = nan(N_a_f_0_test, N_beta_test);
epsilon_L_0_test = nan(N_a_f_0_test, N_beta_test);
a_f_tlp_test = nan(N_a_f_0_test, N_beta_test);

N_total_test = N_a_f_0_test * N_beta_test;
iter = 0;
dperc = 1;
perc_disp = dperc;

for i = 1:N_a_f_0_test
    
    a_f_0 = a_f_0_test(i); 
    
    for j = 1:N_beta_test
        
        iter = iter + 1;
        perc = 100*iter/N_total_test;
        if perc >= perc_disp
            disp(['Solving for sum-square error over range of a_f_0 and beta: ', num2str(perc), '% done'])
            perc_disp = dperc*(ceil(perc/dperc)+1);
        end
        
        beta = beta_test(j); 
        
        % solve for symplastic relative water content
        RWC_s_pred = nan(N_data_subset, 1);
        for k = 1:N_data_subset
            
            RWC_t = RWC_t_subset(k);
            diff_RWC_s_func_local = @(RWC_s) diff_RWC_s_func(RWC_s, RWC_t, a_f_0, beta);
            RWC_s_guess = 1;
            while (abs(imag(diff_RWC_s_func_local(RWC_s_guess))) > 0) || (abs(diff_RWC_s_func_local(RWC_s_guess)) == inf) || isnan(diff_RWC_s_func_local(RWC_s_guess))
                RWC_s_guess = RWC_s_guess - 1e-2;
                if RWC_s_guess < 0
                    break
                end
            end
            
            if RWC_s_guess >= 0
                [RWC_s_pred(k), ~, exitflag, ~] = fzero(diff_RWC_s_func_local, RWC_s_guess, options);
            else
                exitflag = -1;
            end
            
            if exitflag ~= 1
                %if fzero could not find the root, then solve again in
                %alternative way
                RWC_s = 0:1e-3:1;
                diff_RWC_s_local = diff_RWC_s_func_local(RWC_s);
                is_pass = (abs(diff_RWC_s_local) < inf) .* ~isnan(diff_RWC_s_local) .* (imag(diff_RWC_s_local) == 0);
                if sum(is_pass) == 0
                    % no solution
                    RWC_s_pred(k) = nan;
                else
                    RWC_s = RWC_s(is_pass == 1);
                    diff_RWC_s_local = diff_RWC_s_local(is_pass == 1);
                    if (max(diff_RWC_s_local) > 0) && (min(diff_RWC_s_local) < 0)
                        diff_sign_diff_RWC_s_local = diff(sign(diff_RWC_s_local));
                        ind_1 = find((diff_sign_diff_RWC_s_local ~= 0), 1, 'first');
                        ind_2 = ind_1 + 1;
                        RWC_s_pred(k) = RWC_s(ind_1) + (RWC_s(ind_2) - RWC_s(ind_1))*(0-diff_RWC_s_local(ind_1))/(diff_RWC_s_local(ind_2)-diff_RWC_s_local(ind_1));
                    else
                        % no solution
                        RWC_s_pred(k) = nan;
                    end
                end
            end
        end
        
        % solve for the apoplasm fraction
        a_f_pred = a_f_0*RWC_s_pred.^beta; 
        
        % solve for symplastic relative water content and apoplasm fraction at the turgor loss point
        diff_RWC_s_func_local = @(RWC_s) diff_RWC_s_func(RWC_s, RWC_t_tlp, a_f_0, beta);
        while (abs(imag(diff_RWC_s_func_local(RWC_s_guess))) > 0) || (abs(diff_RWC_s_func_local(RWC_s_guess)) == inf) || isnan(diff_RWC_s_func_local(RWC_s_guess))
            RWC_s_guess = RWC_s_guess - 0.01;
            if RWC_s_guess < 0
                break
            end
        end
        
        if RWC_s_guess >= 0
            [RWC_s_tlp, ~, exitflag, ~] = fzero(diff_RWC_s_func_local, RWC_s_guess, options);
        else
            exitflag = -1;
        end
        
        if exitflag ~= 1
            %if fzero could not find the root, then solve again in
            %alternative way
            RWC_s = 0:1e-3:1;
            diff_RWC_s_local = diff_RWC_s_func_local(RWC_s);
            is_pass = (abs(diff_RWC_s_local) < inf) .* ~isnan(diff_RWC_s_local) .* (imag(diff_RWC_s_local) == 0);
            if sum(is_pass) == 0
                % no solution
                RWC_s_tlp = nan;
            else
                RWC_s = RWC_s(is_pass == 1);
                diff_RWC_s_local = diff_RWC_s_local(is_pass == 1);
                if (max(diff_RWC_s_local) > 0) && (min(diff_RWC_s_local) < 0)
                    diff_sign_diff_RWC_s_local = diff(sign(diff_RWC_s_local));
                    ind_1 = find((diff_sign_diff_RWC_s_local ~= 0), 1, 'first');
                    ind_2 = ind_1 + 1;
                    RWC_s_tlp = RWC_s(ind_1) + (RWC_s(ind_2) - RWC_s(ind_1))*(0-diff_RWC_s_local(ind_1))/(diff_RWC_s_local(ind_2)-diff_RWC_s_local(ind_1));
                else
                    % no solution
                    RWC_s_tlp = nan;
                end
            end
        end
        a_f_tlp = a_f_0*RWC_s_tlp^beta; 
        
        % solve for osmotic potential at full hydration
        pi_L_0 = pi_L_tlp * RWC_s_tlp;
        
        % solve for the elastic modulus at full hydration
        epsilon_L_0 = -pi_L_0/(1 - a_f_0 - RWC_t_tlp + a_f_tlp);
        
        % predict total water potential
        psi_L_MD_pred = pi_L_0./RWC_s_pred + epsilon_L_0*max((RWC_t_subset-a_f_pred-RWC_t_tlp+a_f_tlp), 0);
        
        SSE_test(i,j) = sum(weight_subset .* (psi_L_MD_pred-psi_L_MD_subset).^2);
        pi_L_0_test(i,j) = pi_L_0;
        epsilon_L_0_test(i,j) = epsilon_L_0;
        a_f_tlp_test(i,j) = a_f_tlp;
        
    end
end

[beta_test, a_f_0_test] = meshgrid(beta_test, a_f_0_test);
           
% calculate ''a_f_star''
pi_L_star = -1; %reference osmotic potential [MPa]
a_f_star_test = a_f_0_test .* (pi_L_0_test/pi_L_star).^beta_test;

% log10 of SSE
log10_SSE_test = log10(SSE_test);

% find troughs of log10(SSE) surface
trough = zeros(N_a_f_0_test,N_beta_test);
for i = 2:(N_a_f_0_test-1)
    for j = 2:(N_beta_test-1)
        local = log10_SSE_test(i-1:i+1,j-1:j+1);
        if all(local(5) < local([2,8])) || all(local(5) < local([1,9])) || all(local(5) < local([4,6])) || all(local(5) < local([3,7]))
            trough(i,j) = 1;
        end
    end
end

% mask high SSE
log10_SSE_mask = ceil(min(log10_SSE_test, [], 'all'));
if abs(log10_SSE_mask - min(log10_SSE_test, [], 'all')) < 0.1
    log10_SSE_mask = log10_SSE_mask + 1;
end
log10_SSE_test(log10_SSE_test > log10_SSE_mask) = log10_SSE_mask; %hide larger SSE


% figure of log10(SSE)
figure
hold on
contourf(beta_test, a_f_0_test, log10_SSE_test, 'displayname', 'log_{10}(SSE)');
contour(beta_test, a_f_0_test, trough, [1,1], 'c', 'linewidth', 2, 'displayname', 'Local minimum')
hold off

c = colorbar;
c.Label.String = 'log_{10}(SSE)';
caxis([min(caxis), log10_SSE_mask])
c_ticks = get(c, 'Ticks');
c_ticks = unique([c_ticks, log10_SSE_mask]);
set(c, 'Ticks', c_ticks)
c_ticklabels = get(c, 'TickLabels');
c_ticklabels{end} = [c_ticklabels{end}, '+']; 
set(c, 'TickLabels', c_ticklabels)

xlabel('{\it\beta}')
ylabel('{\ita_f}_{,0}')
legend('location', 'se')

a_f_star_max_solve = max(a_f_0_test(~isnan(SSE_test)));
beta_max_solve = max(beta_test(~isnan(SSE_test)));
ylim([0, a_f_star_max_solve])
xlim([0, beta_max_solve])


% figure of elastic modulus
epsilon_L_0_mask = 100;
epsilon_L_0_masked = epsilon_L_0_test;
epsilon_L_0_masked(epsilon_L_0_masked > epsilon_L_0_mask) = epsilon_L_0_mask;
figure
hold on
contourf(beta_test, a_f_0_test, epsilon_L_0_masked);
contour(beta_test, a_f_0_test, trough, [1,1], 'c', 'linewidth', 2, 'displayname', 'Local minimum')
hold off
xlabel('{\it\beta}')
ylabel('{\ita_f}_{,0}')
c = colorbar;
c.Label.String = '{\it\epsilon_L}_{,0}';
caxis([min(caxis), epsilon_L_0_mask])
c_ticks = get(c, 'Ticks');
c_ticks = unique([c_ticks, epsilon_L_0_mask]);
set(c, 'Ticks', c_ticks)
c_ticklabels = get(c, 'TickLabels');
c_ticklabels{end} = [c_ticklabels{end}, '+']; 
set(c, 'TickLabels', c_ticklabels)
ylim([0, a_f_star_max_solve])
xlim([0, beta_max_solve])


% figure of osmotic potential at full hydration
figure
hold on
contourf(beta_test, a_f_0_test, pi_L_0_test);
contour(beta_test, a_f_0_test, trough, [1,1], 'c', 'linewidth', 2, 'displayname', 'Local minimum')
hold off
xlabel('{\it\beta}')
ylabel('{\ita_f}_{,0}')
c = colorbar;
c.Label.String = '{\it\pi_L}_{,0}';
ylim([0, a_f_star_max_solve])
xlim([0, beta_max_solve])


% figure of a_f_star
figure
hold on
contourf(beta_test, a_f_0_test, a_f_star_test);
contour(beta_test, a_f_0_test, trough, [1,1], 'c', 'linewidth', 2, 'displayname', 'Local minimum')
hold off
xlabel('{\it\beta}')
ylabel('{\ita_f}_{,0}')
c = colorbar;
c.Label.String = '{\ita_f}^*';
ylim([0, a_f_star_max_solve])
xlim([0, beta_max_solve])

pause(0)

end

