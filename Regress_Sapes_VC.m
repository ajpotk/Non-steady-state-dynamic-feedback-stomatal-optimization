function [outputs_Regress_Sapes_VC] = Regress_Sapes_VC(varargin)

%% Default settings
Plot = 1;
Q10_regression = 0;

%% Check ''varargin''
n_varargin = length(varargin);
n_varargin_name_and_value = floor(n_varargin/2);
if 2*n_varargin_name_and_value ~= n_varargin
    error('ERROR: varargin inputs should have an even total number of inputs: each name must has a value!')
end

for i = 1:n_varargin_name_and_value
   
    switch varargin{1+2*(i-1)}
        case 'Plot'
            Plot = varargin{2*i};
            if isnumeric(Plot)
                if (Plot ~=0) && (Plot ~= 1)
                    error('ERROR: ''Plot'' must be a numeric value set to either 0 or 1!')
                end
            else
                error('ERROR: ''Plot'' must be a numeric value set to either 0 or 1!')
            end
        case 'Q10_regression'
            Q10_regression = varargin{2*i};
            if isnumeric(Q10_regression)
                if (Q10_regression ~=0) && (Q10_regression ~= 1)
                    error('ERROR: ''Q10_regression'' must be a numeric value set to either 0 or 1!')
                end
            else
                error('ERROR: ''Q10_regression'' must be a numeric value set to either 0 or 1!')
            end
        otherwise
            if isstring(varargin{1+2*(i-1)}) || ischar(varargin{1+2*(i-1)})
                error(['ERROR: unsupported varargin name: ''', varargin{1+2*(i-1)}, '''', newline, 'Acceptable varargin names are ''Plot'' and ''Q10_regression''.'])
            else
                error('ERROR: odd-numbered varargin inputs must be strings or characters!')
            end
    end
    
end

%% Data location
home = pwd;

directory_data_file = 'C:\Users\potka002\Desktop\Research\Growth Optimizing Stomata\Test_GOH\Sapes & Sala 2021 Plant, Cell & Environment';
data_file = 'Master_Data.xlsx';
data_sheet = 'Mastersheet'; 

%% Download data
disp('Reading Sapes & Sala (2021) data')
cd(directory_data_file);
opts = detectImportOptions(data_file, 'Sheet', data_sheet);
raw_data = readtable(data_file, opts);
cd(home)

%% Collect specified variables from raw data
data_var = {'Trmmol', ...
            'Leaf_WP_MPa', ...
            'Soil_WP_Mpa_from_VWC', ...
            'Tleaf'};

data_multiplier = [1e-3, ... %mmol m-2s-1 to mol m-2 s-1
                   nan(1,3)];

raw_data =  table2cell(raw_data);
[M_data, ~] = size(raw_data);
N_data_var = length(data_var);
raw_data_var = opts.VariableNames;
data = cell(M_data, N_data_var);

% correct data by ''data_multiplier''
for i = 1:N_data_var
    is_data_var = cell2mat(cellfun(@(x) strcmp(x, data_var{i}), raw_data_var, 'UniformOutput', 0)); 
    data(:,i) = raw_data(:,is_data_var == 1); 
    if ~isnan(data_multiplier(i))
        data(:,i) = cellfun(@(x) x * data_multiplier(i), data(:,i), 'UniformOutput', 0);
    end
end

% remove nan's
is_nan = cell2mat(cellfun(@(x) isnan(x), data, 'UniformOutput', 0));
is_nan = any(is_nan, 2);
ind_not_nan = (1 - is_nan) .* (1:M_data)';
ind_not_nan(ind_not_nan == 0) = [];
data = cell2mat(data(ind_not_nan, :));

E = data(:, 1);
psi_L = data(:, 2);
psi_soil = data(:, 3);
psi_soil(psi_soil > 0) = 0;
T_L = data(:, 4);

% remove negative E
is_not_neg_E = (E >=  0);
is_soil_drier_than_leaf = (psi_soil >= psi_L);
is_all = is_not_neg_E .* is_soil_drier_than_leaf;
E = E(is_all == 1);
psi_L = psi_L(is_all == 1);
psi_soil = psi_soil(is_all == 1);
T_L = T_L(is_all == 1);

% raw data before binning
E_raw = E;
psi_L_raw = psi_L;
psi_soil_raw = psi_soil;
T_L_raw = T_L;

% bin data by psi_soil and psi_L
psi_bin = 0.05;
psi_soil_bin_ub = psi_bin*(ceil(min(psi_soil)/psi_bin)):psi_bin:0;
psi_soil_bin_lb = psi_soil_bin_ub - psi_bin;
psi_soil_bin_mid = (psi_soil_bin_ub + psi_soil_bin_lb)/2;
N_psi_soil_bin = length(psi_soil_bin_ub);
psi_L_bin_ub = psi_bin*(ceil(min(psi_L)/psi_bin)):psi_bin:0;
psi_L_bin_lb = psi_L_bin_ub - psi_bin;
psi_L_bin_mid = (psi_L_bin_ub + psi_L_bin_lb)/2;
N_psi_L_bin = length(psi_L_bin_ub);
[psi_L_bin_ub, psi_soil_bin_ub] = meshgrid(psi_L_bin_ub, psi_soil_bin_ub);
[psi_L_bin_lb, psi_soil_bin_lb] = meshgrid(psi_L_bin_lb, psi_soil_bin_lb);
[psi_L_bin_mid, psi_soil_bin_mid] = meshgrid(psi_L_bin_mid, psi_soil_bin_mid);
E_bin = nan(N_psi_soil_bin, N_psi_L_bin);
for i = 1:N_psi_soil_bin
    for j = 1:N_psi_L_bin
        
        ind_data_local = find((psi_soil >= psi_soil_bin_lb(i,j)) + (psi_soil < psi_soil_bin_ub(i,j)) + (psi_L >= psi_L_bin_lb(i,j)) + (psi_L < psi_L_bin_ub(i,j)) == 4);
        if ~isempty(ind_data_local)
            E_bin(i,j) = mean(E(ind_data_local));
        end
    end
end

psi_soil_bin_mid = reshape(psi_soil_bin_mid, N_psi_soil_bin*N_psi_L_bin, 1);
psi_L_bin_mid = reshape(psi_L_bin_mid, N_psi_soil_bin*N_psi_L_bin, 1);
E_bin = reshape(E_bin, N_psi_soil_bin*N_psi_L_bin, 1);

ind_not_nan = ~isnan(E_bin);
psi_soil = psi_soil_bin_mid(ind_not_nan == 1);
psi_L = psi_L_bin_mid(ind_not_nan == 1);
E = E_bin(ind_not_nan == 1);

% leaf-area specific soil-plant conductance
k_L = E ./ (psi_soil - psi_L);
k_L(k_L == inf) = nan;
k_L(k_L < 0) = nan;

% % % %% Attempt 1
% % % % for k_L = k_L_max / (1 + exp(-A*(psi - B)))
% % % % C(1) = k_L_max
% % % % C(2) = A
% % % % C(3) = B
% % % x = [psi_soil, psi_L];
% % % E_func = @(C, x) ModelFunction1(C, x);
% % % dE = 1e-20;
% % % log_E_func = @(C, x) log(E_func(C, x) + dE);
% % % lb = [0, 0, -inf];
% % % ub = [inf, inf, 0];
% % % C0 = [max(k_L)/2, 2, -2];
% % % C = lsqcurvefit(log_E_func, C0, x, log(E + dE), lb, ub);
% % % 
% % % SSE = sum((E - E_func(C, x)).^2);
% % % disp('MODEL 1')
% % % disp(['    SSE = ', num2str(SSE)])
% % % 
% % % 
% % % figure
% % % tiledlayout(1,2)
% % % 
% % % lims = [1e-20, 5e-3];
% % % 
% % % nexttile(1)
% % % hold on
% % % plot(lims, lims, 'k--')
% % % plot(E, E_func(C, x), 'bo', 'markerfacecolor', 'b', 'DisplayName', 'Single sigmoidal soil-plant conductance')
% % % hold off
% % % xlabel('Measured {\itE}')
% % % ylabel('Predicted {\itE}')
% % % set(gca, 'yScale', 'log')
% % % set(gca, 'xScale', 'log')
% % % xlim(lims)
% % % ylim(lims)
% % % 
% % % nexttile(2)
% % % hold on
% % % plot(-psi_L, E, 'ks', 'markerfacecolor', 'k', 'DisplayName', 'Data')
% % % plot(-psi_L, E_func(C, x), 'bo', 'markerfacecolor', 'b', 'DisplayName', 'Single sigmoidal soil-plant conductance')
% % % hold off
% % % xlabel('-{\it\psi_L}')
% % % ylabel('{\itE}')
% % % 
% % % %% Attempt 2
% % % % for 1/k_L = 1/k_soil + 1/k_plant
% % % % where k_soil = k_soil_max * exp(A*psi)
% % % % and k_plant = constant
% % % % C(1) = k_soil_max
% % % % C(2) = A
% % % % C(3) = k_plant
% % % x = [psi_soil, psi_L];
% % % E_func = @(C, x) ModelFunction2(C, x);
% % % dE = 1e-20;
% % % log_E_func = @(C, x) log(E_func(C, x) + dE);
% % % lb = [0, 0, 0];
% % % ub = [inf, inf, inf];
% % % C0 = [max(k_L), 1, max(k_L)];
% % % C = lsqcurvefit(log_E_func, C0, x, log(E + dE), lb, ub);
% % % 
% % % SSE = sum((E - E_func(C, x)).^2);
% % % disp('MODEL 2')
% % % disp(['    SSE = ', num2str(SSE)])
% % % 
% % % 
% % % nexttile(1)
% % % hold on
% % % plot(E, E_func(C, x), 'ro', 'markerfacecolor', 'r', 'DisplayName', 'Exponential soil conductance & constant plant conductance')
% % % hold off
% % % 
% % % nexttile(2)
% % % hold on
% % % plot(-psi_L, E_func(C, x), 'ro', 'markerfacecolor', 'r', 'DisplayName', 'Exponential soil conductance & constant plant conductance')
% % % hold off

if Q10_regression == 0
    
    %% Attempt 3
    % for k_L = k_L_max * exp(A*psi_L)
    [p, p_conf_interv] = regress(log(k_L(~isnan(k_L))), [psi_L(~isnan(k_L)), ones(sum(~isnan(k_L)), 1)]);
    A_conf_interv = p_conf_interv(1,:);
    k_L_max_conf_interv = exp(p_conf_interv(2,:));
    % % % p = polyfit(psi_L(~isnan(k_L)), log(k_L(~isnan(k_L))), 1);
    A = p(1);
    k_L_max = exp(p(2));
    C = [k_L_max, A];
    x = [psi_soil, psi_L];
    E_func = @(C, x) C(1)*exp(C(2)*x(:,2)).*(x(:,1)-x(:,2));

    if Plot == 1

        SSE = sum((E - E_func(C, x)).^2);
        disp('MODEL 3')
        disp(['    SSE = ', num2str(SSE)])

        figure
        tiledlayout(2,1)

        nexttile(1)
        hold on
        plot(E, E_func(C, x), 'mo', 'markerfacecolor', 'm', 'DisplayName', 'Single exponential soil-plant conductance')
        hold off
        ylim('Predicted {\itE}')
        xlim('Observed {\itE}')


        nexttile(2)
        hold on
        plot(-psi_L, E_func(C, x), 'mo', 'markerfacecolor', 'm', 'DisplayName', 'Single exponential soil-plant conductance')
        hold off
        ylim('Predicted {\itE}')
        xlim('Observed -{\it\psi_L}')

        % other figure
        X = -4.5:0.01:0;
        figure
        hold on
        plot(-psi_L, k_L, 'ks', 'markerfacecolor', 'k', 'DisplayName', 'Data')
        plot(-X, k_L_max*exp(A*X), '-m', 'linewidth', 2)
        hold off
        xlabel('-{\it\psi_L}')
        ylabel('{\itk_L}')

    end

elseif Q10_regression == 1

    %% Attempt 4
    % for k_L = k_L_max_25C * Q10^((T - 25)/10) * exp(A*psi_L)
    k_L_raw = E_raw./(psi_soil_raw - psi_L_raw);
    [p, p_conf_interv] = regress(log(k_L_raw), [psi_L_raw, ones(length(k_L_raw), 1), (T_L_raw-25)/10]);
    A = p(1);
    k_L_max = exp(p(2)); %value at 25C
    Q10 = exp(p(3)); 
    A_conf_interv = p_conf_interv(1,:);
    k_L_max_conf_interv = exp(p_conf_interv(2,:)); %value at 35C
    Q10_conf_interv = exp(p_conf_interv(3,:));
    
    % predict k_L
    log_k_L_pred = p(2) + A*psi_L_raw + p(3)*(T_L_raw-25)/10;
    std_A = diff(A_conf_interv)/1.96/2; 
    std_log_k_L_max = diff(p_conf_interv(2,:))/1.96/2;
    std_log_Q10 = diff(p_conf_interv(3,:))/1.96/2;
    std_log_k_L_pred = (std_log_k_L_max^2 + psi_L_raw.^2*std_A^2 + ((T_L_raw-25)/10).^2*std_log_Q10^2).^0.5;
    log_k_L_pred_conf_interv = log_k_L_pred + 1.96*std_log_k_L_pred.*[-1, 1];
    log_k_L_raw = log(k_L_raw);
    is_pred_okay = ((log_k_L_pred_conf_interv(:,1) < log_k_L_raw) + (log_k_L_pred_conf_interv(:,2) > log(k_L_raw)) == 2); %is within confidence intervals
    
    if Plot == 1
        
        N_pred = length(k_L_raw); 
        xy_min = floor(min([log_k_L_raw, log_k_L_pred_conf_interv(:,1)], [], 'all'));
        xy_max = ceil(max([log_k_L_raw, log_k_L_pred_conf_interv(:,2)], [], 'all'));
        tick_length = 0.1;
        
        figure
        box on
        hold on
        plot([xy_min, xy_max], [xy_min, xy_max], 'r--', 'LineWidth', 2)
        for i = 1:N_pred
            plot(log_k_L_raw(i)+tick_length*[1,-1], log_k_L_pred_conf_interv(i,1)*[1,1], '-k', 'LineWidth', 0.5)
            plot(log_k_L_raw(i)+tick_length*[1,-1], log_k_L_pred_conf_interv(i,2)*[1,1], '-k', 'LineWidth', 0.5)
            plot(log_k_L_raw(i)*[1,1], log_k_L_pred_conf_interv(i,1:2), '-k', 'LineWidth', 0.5)
        end
        plot(log_k_L_raw(is_pred_okay == 1), log_k_L_pred(is_pred_okay == 1), 'ko', 'LineStyle', 'none', 'MarkerFaceColor', 'k')
        plot(log_k_L_raw(is_pred_okay == 0), log_k_L_pred(is_pred_okay == 0), 'ko', 'MarkerFaceColor', 'w')
        hold off
        xlim([xy_min, xy_max])
        ylim([xy_min, xy_max])
        xlabel(['ln({\itk_{L,observed}})', newline, '[ln(mol\cdotm^{-2}\cdots^{-1}\cdotMPa^{-1})]'])
        ylabel(['ln({\itk_{L,predicted}})', newline, '[ln(mol\cdotm^{-2}\cdots^{-1}\cdotMPa^{-1})]'])
        
    end
    
end

psi_L_crit_raw = psi_soil_raw - 1/A;
psi_L_crit_raw_conf_interv = psi_soil_raw - 1./A_conf_interv;

%% outputs
outputs_Regress_Sapes_VC.E_raw = E_raw;
outputs_Regress_Sapes_VC.psi_L_raw = psi_L_raw;
outputs_Regress_Sapes_VC.psi_soil_raw = psi_soil_raw;
outputs_Regress_Sapes_VC.T_L_raw = T_L_raw;
outputs_Regress_Sapes_VC.A = A;
outputs_Regress_Sapes_VC.k_L_max = k_L_max;
outputs_Regress_Sapes_VC.A_conf_interv = A_conf_interv;
outputs_Regress_Sapes_VC.k_L_max_conf_interv = k_L_max_conf_interv;
outputs_Regress_Sapes_VC.psi_L_crit_raw = psi_L_crit_raw;
outputs_Regress_Sapes_VC.psi_L_crit_raw_conf_interv = psi_L_crit_raw_conf_interv;

% % % %% Local functions
% % %     function a = ModelFunction1(C, x)
% % %         modelfun1 = @(C, x) C(1)/C(2)*log((1 + exp(-C(2)*(x(:,1) - C(3))))./(1 + exp(-C(2)*(x(:,2) - C(3))))) + C(1)*(x(:,1) - x(:,2));
% % %         a = modelfun1(C, x);
% % %         a(a < 0) = 0;
% % %     end
% % % 
% % %     function a = ModelFunction2(C, x)
% % %         modelfun2 = @(C, x) C(1)/C(2)*exp(C(2)*x(:,1)) - C(3)/C(2)*lambertw(C(1)/C(3)*exp(C(2)*x(:,2) + C(1)/C(3)*exp(C(2)*x(:,1))));
% % %         a = modelfun2(C, x);
% % %         a(a < 0) = 0;
% % %     end

end

