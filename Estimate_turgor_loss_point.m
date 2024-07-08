function [RWC_t_tlp, pi_L_tlp, pi_L_0_trad, a_f_trad] = Estimate_turgor_loss_point(RWC_t_subset, psi_L_MD_subset, varargin)

X_max = inf; %default value for ''X_max''
X_min = -inf; %default value for ''X_max''

for i = 1:ceil(length(varargin)/2)
    switch varargin{1 + 2*(i-1)}
        case 'X_max'
            X_max = varargin{2*i};
        case 'X_min'
            X_min = varargin{2*i};
        otherwise
            if isstring(varargin{1 + 2*(i-1)}) || ischar(varargin{1 + 2*(i-1)})
                error(['ERROR: ''Estimate_turgor_loss_point'' does not accept ''', varargin{1 + 2*(i-1)}, ''' as an input!'])
            else
                error('ERROR: Odd-numbered ''varargin'' inputs must be strings or characters!')
            end
    end
end

Y = -1./psi_L_MD_subset;
X = 1 - RWC_t_subset;
X(X > X_max) = nan; 
X(X < X_min) = nan; 

% remove nan's
is_not_nan = (~isnan(Y)) .* (~isnan(X));
Y = Y(is_not_nan == 1);
X = X(is_not_nan == 1);

[X, ind_sort] = sort(X);
Y = Y(ind_sort);

%% Determine turgor loss point
X_tlp_test = 0:1e-3:X(end - 2);
N_test = length(X_tlp_test);
X_tlp_max = 0.20; %physiological bounds for ''X''
X_tlp_min = 0.02; %physiological bounds for ''X''
X_tlp_default = 0.1;
SSE = nan(N_test, 1); %sum-square error
is_slope_neg_lin = ones(N_test, 1);
is_slope_neg_log = ones(N_test, 1);
for k = 1:N_test

    X_lin = X(X >= X_tlp_test(k));
    Y_lin = Y(X >= X_tlp_test(k));
    p_lin = polyfit(X_lin, Y_lin, 1);
    Y_lin_pred = polyval(p_lin, X_lin);

    X_log = X(X < X_tlp_test(k));
    Y_log = Y(X < X_tlp_test(k));
    Y_tlp = polyval(p_lin, X_tlp_test(k));
    b = (X_log-X_tlp_test(k))\log(Y_log/Y_tlp); 
    Y_log_pred = Y_tlp*exp(b*(X_log-X_tlp_test(k)));
    
    if p_lin(1) >= 0
        is_slope_neg_lin(k) = 0;
    end
    if b >= 0
        is_slope_neg_log(k) = 0;
    end

    SSE(k) = sum((Y_lin_pred - Y_lin).^2) + sum((Y_log_pred - Y_log).^2);

end
are_slopes_neg = is_slope_neg_lin .* is_slope_neg_log;
ind_best = find((SSE == min(SSE(are_slopes_neg == 1))), 1, 'first');
X_tlp = X_tlp_test(ind_best);

% plot fit
figure
hold on
plot(X,Y, 'bo', 'markerfacecolor', 'b') 


if (X_tlp >= X_tlp_min) && (X_tlp <= X_tlp_max) %based on realistic bounds for RWC_t_tlp from Bartlett et al. (2012)
    % if satisfied with estimate of ''X_tlp''
    X_lin = X(X >= X_tlp);
    Y_lin = Y(X >= X_tlp);
    p_lin = polyfit(X_lin, Y_lin, 1);
    
    X_log = X(X < X_tlp);
    Y_log = Y(X < X_tlp);
    Y_tlp = polyval(p_lin, X_tlp);
    b = (X_log-X_tlp)\log(Y_log/Y_tlp);
    
    Y_lin_pred = polyval(p_lin, X_lin);
    Y_lin_extrap = polyval(p_lin, X_log);
    Y_log_pred = Y_tlp*exp(b*(X_log-X_tlp));
    
    plot(X_lin, Y_lin_pred, 'r-', 'linewidth', 2)
    plot(X_log, Y_log_pred, 'r-', 'linewidth', 2)
    plot(X_log, Y_lin_extrap, 'k:', 'linewidth', 1.5)
    
    r = corrcoef([Y_lin; Y_log], [Y_lin_pred; Y_log_pred]);
    r2 = r(1,2)^2;
    title(['{\itr}^2 = ', num2str(r2, 2)])
    
else
    warning(['WARNING: ''X_tlp'' was estimated to be ', num2str(X_tlp), '. Resetting ''X_tlp'' to ', num2str(X_tlp_default), '.'])
    % otherwise assume valye of ''X_tlp''
    X_tlp = X_tlp_default;
    p_lin = polyfit(X((X >= X_tlp_min) .* (X <= X_tlp_max) == 1), Y((X >= X_tlp_min) .* (X <= X_tlp_max) == 1), 1); 
end
hold off
xlabel('1 - RWC_{\itT}')
ylabel('-1/{\it\psi_L} [MPa^{-1}]')

Y_lin_func = @(X) polyval(p_lin, X); 
Y_tlp = Y_lin_func(X_tlp);
Y_0 = Y_lin_func(0);

RWC_t_tlp = 1 - X_tlp;
pi_L_tlp = -1./Y_tlp; 
pi_L_0_trad = -1./Y_0;  %traditional estimate of ''pi_L_0'', ignoring varying apoplasm fraction


% estimate the apoplasm fraction through the traditional approach
X_a_f = fzero(Y_lin_func, X_tlp);
a_f_trad = 1 - X_a_f;

end

