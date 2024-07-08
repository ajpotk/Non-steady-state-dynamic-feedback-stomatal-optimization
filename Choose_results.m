function [ind_best, ind_best_sub, N_best_sub, ...
          a_f_max_results_sub, pi_L_star_results_sub, beta_results_sub, epsilon_L_max_results_sub, SWC_L_results_sub, alpha_results_sub, ...
          a_f_max_best_sub, pi_L_star_best_sub, beta_best_sub, epsilon_L_max_best_sub, SWC_L_best_sub, alpha_best_sub] = ...
          Choose_results(pi_L_0_25C_Trt_combo_results, SSE_best_total_subset_results, N_best, species_subset, ...
                         a_f_max_results, pi_L_star_results, beta_results, epsilon_L_max_results, SWC_L_results, alpha_results, ...
                         varargin)

%% First constrain results to the ''N_best'' resutls with the lowest SSE
[~, ind_best] = sort(SSE_best_total_subset_results);
ind_best = ind_best(1:N_best);

%% Option for choosing subset and best parameters
% Option = 1 -- subset based on data with pi_L_0_25C within a specified range (based on the standard deviation) of its mean value from DREAM
% Option = 2 -- no subsetting and best parameters taken as mode (if not uniformly distributed) of values from DREAM
% Option = 3 -- no subsetting and best parameters taken as best

Option = 3; %default setting -- may be changed through varargin below

%% Check ''varargin''
n_varargin = length(varargin);
n_varargin_name_and_value = floor(n_varargin/2);
if 2*n_varargin_name_and_value ~= n_varargin
    error('ERROR: varargin inputs should have an even total number of inputs: each name must has a value!')
end

for i = 1:n_varargin_name_and_value
   
    switch varargin{1+2*(i-1)}
        case 'Option'
            Option = varargin{2*i};
        otherwise
            if isstring(varargin{1+2*(i-1)}) || ischar(varargin{1+2*(i-1)})
                error(['ERROR: unsupported varargin name: ''', varargin{1+2*(i-1)}, '''', newline, 'Acceptable varagin names are ''Option'''])
            else
                error('ERROR: odd-numbered varargin inputs must be strings or characters!')
            end
    end
    
end


%% Choose results
switch Option
    case 1
        
        %% Find data within a specified range of the mean of ''pi_L_0_25C_Trt_combo_results''
        N_std = 0.25;
        pi_L_0_25C_Trt_combo_results = pi_L_0_25C_Trt_combo_results(ind_best,:);
        pi_L_0_25C_Trt_combo_mean = mean(pi_L_0_25C_Trt_combo_results);
        pi_L_0_25C_Trt_combo_std = std(pi_L_0_25C_Trt_combo_results);
        pi_L_0_25C_Trt_combo_upper = pi_L_0_25C_Trt_combo_mean + N_std*pi_L_0_25C_Trt_combo_std;
        pi_L_0_25C_Trt_combo_lower = pi_L_0_25C_Trt_combo_mean - N_std*pi_L_0_25C_Trt_combo_std;
        is_in_range = ((pi_L_0_25C_Trt_combo_results <= pi_L_0_25C_Trt_combo_upper) + (pi_L_0_25C_Trt_combo_results >= pi_L_0_25C_Trt_combo_lower) == 2);
        is_col_nan = all(isnan(pi_L_0_25C_Trt_combo_results)); %denotes if entire column is missing/nan
        for i = 1:length(is_col_nan)
           if is_col_nan(i)
               is_in_range(:,i) = 1; %ignore columns
           end
        end
        ind_best_sub = ind_best(all(is_in_range, 2));
        N_best_sub = length(ind_best_sub);

        %% Figure of PDFs of two data sets
        % first check if plot is already saved
        Output_file = [species_subset, '_Fig_SSE-PDF_for_', num2str(N_best), '_best.jpg'];
        files = struct2cell(dir);
        files = files(1,:);
        is_output_file = cell2mat(cellfun(@(x) strcmp(x, Output_file), files, 'UniformOutput', 0));
        if ~any(is_output_file)

            figure
            hold on
            histogram(SSE_best_total_subset_results(ind_best), 'Normalization', 'pdf', 'Displayname', [num2str(N_best), ' best data'])
            histogram(SSE_best_total_subset_results(ind_best_sub), 'Normalization', 'pdf', 'Displayname', ['Data with {\it\pi_L}_{0,25^{\circ}C} within \pm', num2str(N_std), ' std of mean {\it\pi_L}_{0,25^{\circ}C}'])
            hold off

            xlabel('SSE [mol^2\cdotm^{-4}\cdots^{-2}]')
            ylabel('PDF')
            legend('Location', 'NW')
            legend boxoff

            print(gcf, '-djpeg', '-r300', Output_file, '-painters', '-noui' )

        end

        %% Parameters in subset of data within the specified range of the mean of ''pi_L_0_25C_Trt_combo_results''
        a_f_max_results_sub = a_f_max_results(ind_best_sub);
        pi_L_star_results_sub = pi_L_star_results(ind_best_sub);
        beta_results_sub = beta_results(ind_best_sub);
        epsilon_L_max_results_sub = epsilon_L_max_results(ind_best_sub);
        SWC_L_results_sub = SWC_L_results(ind_best_sub);
        alpha_results_sub = alpha_results(ind_best_sub);

        %% Choose best parameters from subset of data within the specified range of the mean of ''pi_L_0_25C_Trt_combo_results''
        % choose 1st value, since already sorted by SSE
        a_f_max_best_sub = a_f_max_results_sub(1);
        pi_L_star_best_sub = pi_L_star_results_sub(1);
        beta_best_sub = beta_results_sub(1);
        epsilon_L_max_best_sub = epsilon_L_max_results_sub(1);
        SWC_L_best_sub = SWC_L_results_sub(1);
        alpha_best_sub = alpha_results_sub(1);

    case 2
        
        %% No subsetting
        ind_best_sub = ind_best;
        N_best_sub = N_best;
        a_f_max_results_sub = a_f_max_results(ind_best_sub);
        pi_L_star_results_sub = pi_L_star_results(ind_best_sub);
        beta_results_sub = beta_results(ind_best_sub);
        epsilon_L_max_results_sub = epsilon_L_max_results(ind_best_sub);
        SWC_L_results_sub = SWC_L_results(ind_best_sub);
        alpha_results_sub = alpha_results(ind_best_sub);
        
        %% If not uniformly distributed, then set best parameter value as mode of Kernel density
        results = [reshape(a_f_max_results_sub, N_best_sub, 1), ...
                   reshape(pi_L_star_results_sub, N_best_sub, 1), ...
                   reshape(beta_results_sub, N_best_sub, 1), ...
                   reshape(epsilon_L_max_results_sub, N_best_sub, 1), ...
                   reshape(SWC_L_results_sub, N_best_sub, 1), ...
                   reshape(alpha_results_sub, N_best_sub, 1)];
        
        [~, N_param] = size(results);
        results_best = nan(1, 5);
        for i = 1:N_param
            
            results_local = results(:,i);
            x = linspace(min(results_local), max(results_local), 5e2);
            pd = fitdist(results_local,'Kernel','Kernel','epanechnikov');
            PDF = pdf(pd, x);
            [~, ind_PDF_max] = sort(PDF, 'descend'); 
            results_best(i) = x(ind_PDF_max(1)); 
            
        end
        
        a_f_max_best_sub = results_best(1);
        pi_L_star_best_sub = results_best(2);
        beta_best_sub = results_best(3);
        epsilon_L_max_best_sub = results_best(4);
        SWC_L_best_sub = results_best(5);
        alpha_best_sub = results_best(6);
        
    case 3
        
        %% No subsetting
        ind_best_sub = ind_best;
        N_best_sub = N_best;
        a_f_max_results_sub = a_f_max_results(ind_best_sub);
        pi_L_star_results_sub = pi_L_star_results(ind_best_sub);
        beta_results_sub = beta_results(ind_best_sub);
        epsilon_L_max_results_sub = epsilon_L_max_results(ind_best_sub);
        SWC_L_results_sub = SWC_L_results(ind_best_sub);
        alpha_results_sub = alpha_results(ind_best_sub);
        
        %% Choose best parameters from subset of data
        % choose 1st value, since already sorted by SSE
        a_f_max_best_sub = a_f_max_results_sub(1);
        pi_L_star_best_sub = pi_L_star_results_sub(1);
        beta_best_sub = beta_results_sub(1);
        epsilon_L_max_best_sub = epsilon_L_max_results_sub(1);
        SWC_L_best_sub = SWC_L_results_sub(1);
        alpha_best_sub = alpha_results_sub(1);
        
end


end

