function [V_cmax25_func] = Regress_Vcmax25_as_func_psi(T_L_subset, V_cmax_subset, psi_L_MD_subset, ...
                                                       varargin)

%% Constants
R = 8.314; %universal gas constant [J mol-1 K-1]

%% Default settings
Mean_value = 0;

%% Check ''varargin''
n_varargin = length(varargin);
n_varargin_name_and_value = floor(n_varargin/2);
if 2*n_varargin_name_and_value ~= n_varargin
    error('ERROR: varargin inputs should have an even total number of inputs: each name must has a value!')
end

for i = 1:n_varargin_name_and_value
   
    switch varargin{1+2*(i-1)}
        case 'Mean_value'
            Mean_value = varargin{2*i};
            if isnumeric(Mean_value)
                if (Mean_value ~=0) && (Mean_value ~= 1)
                    error('ERROR: ''Mean_value'' must be a numeric value set to either 0 or 1!')
                end
            else
                error('ERROR: ''Mean_value'' must be a numeric value set to either 0 or 1!')
            end
        otherwise
            if isstring(varargin{1+2*(i-1)}) || ischar(varargin{1+2*(i-1)})
                error(['ERROR: unsupported varargin name: ''', varargin{1+2*(i-1)}, '''', newline, 'Acceptable varargin names are ''Mean_value''.'])
            else
                error('ERROR: odd-numbered varargin inputs must be strings or characters!')
            end
    end
    
end

%% Regress V_cmax25 as function of psi_L_MD
V_cmax_per_V_cmax25_subset = exp(65330*(T_L_subset-25)/298.15/R./(T_L_subset+273.15));
V_cmax25_subset = V_cmax_subset ./ V_cmax_per_V_cmax25_subset;

if Mean_value == 1
    
    V_cmax25_func = @(x) mean(V_cmax25_subset) * ones(size(x));
    
else
    
    % V_cmax25_0 == B(1), S_f == B(2), psi_f == B(3);
    % % % modelfun_crrct =  @(B, x) B(1) .* (1 + exp(B(2)*B(3))) ./ (1 + exp(B(2)*(B(3)-x)));
    modelfun_crrct = @(B, x) ModelFunction(B, x); %see local function at bottom of script
    V_cmax25_0_guess = max(V_cmax25_subset);
    p = polyfit(psi_L_MD_subset, log(V_cmax25_subset), 1);
    S_f_guess = 2*p(1); %since d(ln(V_cmax25))/dpsi_L = S_f/2 at psi_L = psi_f
    S_f_guess(S_f_guess < 0) = 0;
    psi_f_guess = log((V_cmax25_subset/V_cmax25_0_guess - 1)./(1 - V_cmax25_subset/V_cmax25_0_guess.*exp(-S_f_guess*psi_L_MD_subset)))/S_f_guess; %by rearrangement of function
    psi_f_guess(abs(imag(psi_f_guess)) > 0) = nan;
    psi_f_guess(abs(real(psi_f_guess)) == inf) = nan;
    psi_f_guess = median(psi_f_guess, 'omitnan');
    psi_f_guess(psi_f_guess > 0) = 0;
    B_guess = [V_cmax25_0_guess, S_f_guess, psi_f_guess];
    options.Lower = [0, 0, -inf];
    options.Upper = [inf, inf, 0];
    B = nlinfit(psi_L_MD_subset, V_cmax25_subset, modelfun_crrct, B_guess, options);
    V_cmax25_func = @(x) modelfun_crrct(B, x); 

end

% % % x = min(psi_L_MD_subset):0.01:0;
% % % figure
% % % hold on
% % % plot(psi_L_MD_subset, V_cmax25_subset, 'bo', 'MarkerFaceColor', 'b')
% % % plot(x, modelfun_crrct(B_guess, x), 'r:', 'linewidth', 2)
% % % plot(x, V_cmax25_func(x), 'k-', 'linewidth', 3)
% % % hold off



%% Local function
    function a = ModelFunction(B, x)
        modelfun = @(B, x) B(1) .* (1 + exp(B(2)*B(3))) ./ (1 + exp(B(2)*(B(3)-x)));
        a = modelfun(B, x);
        a(isnan(a)) = 0;
        a(abs(a) == inf) = 0;
    end

end

