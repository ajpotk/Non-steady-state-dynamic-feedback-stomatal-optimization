function [Var] = Download_B4W_data(varargin)

%% Check varargin
[Var] = Check_varargin_for_download(varargin);

%% Check for past data .mat file
Mat_data_save_file = 'Saved_B4W_data.mat';

files = struct2cell(dir);
files = files(1,:);
is_Mat_data_save_file = cell2mat(cellfun(@(x) strcmp(x, Mat_data_save_file), files, 'UniformOutput', 0));

if any(is_Mat_data_save_file)
    
    load(Mat_data_save_file)
    
else
    %% Constants
    o_i = 1e-3 * 210; %standard atmospheric partial pressure of oxygen [mol mol-1]
    R = 8.314; %universal gas constant [J mol-1 K-1]

    %% Data location
    home = pwd;

    directory_data_file = 'C:\Users\potka002\Desktop\Research\Growth Optimizing Stomata\Test_GOH\B4W data-20220706T114512Z-001\B4W data';
    data_file = '2013 water potential and gas exchange.csv';

    %% Download data
    disp('Reading B4W data')
    cd(directory_data_file);
    opts = detectImportOptions(data_file);
    raw_data = readtable(data_file, opts);
    cd(home)

    %% Collect specified variables from raw data
    data_var = {'Species', ...
                'PlotID', ...
                'WarmingTrt', ...
                'heatTrt', ...
                'WaterTrt', ...
                'Photo', ... %µmol m-2 s-1
                'Pre_dawnWaterPotential', ... %MPa
                'Mid_dayWaterPotential',... %MPa
                'Cond', ... %mol m-2 s-1
                'Ci', ... %µmol mol-1
                'CO2S', ... %µmol mol-1
                'Trmmol', ... %mmol m-2 s-1
                'VpdL', ... %kPa
                'Tleaf'}; %C

    data_multiplier = [nan(1,5), ...
                       1e-6, ... %µmol m-2s-1 to mol m-2 s-1
                       nan(1,3), ...
                       1e-6 * ones(1,2), ... %µmol mol-1 to mol mol-1
                       1e-3, ...%mmol m-2s-1 to mol m-2 s-1
                       nan(1,2)];

    raw_data =  table2cell(raw_data);
    [M_data, N_raw_data_var] = size(raw_data);
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

    %% Calculate leaf temperature-dependent photosynthetic traits from Bernacchi et al. (2001)

    is_data_var = cell2mat(cellfun(@(x) strcmp(x, 'Tleaf'), data_var, 'UniformOutput', 0));
    T_L = cell2mat(data(:,is_data_var == 1));

    is_data_var = cell2mat(cellfun(@(x) strcmp(x, 'Photo'), data_var, 'UniformOutput', 0));
    A_n = cell2mat(data(:,is_data_var == 1));

    is_data_var = cell2mat(cellfun(@(x) strcmp(x, 'Ci'), data_var, 'UniformOutput', 0));
    c_i = cell2mat(data(:,is_data_var == 1));

    is_data_var = cell2mat(cellfun(@(x) strcmp(x, 'CO2S'), data_var, 'UniformOutput', 0));
    c_a = cell2mat(data(:,is_data_var == 1));

    is_data_var = cell2mat(cellfun(@(x) strcmp(x, 'Cond'), data_var, 'UniformOutput', 0));
    g_w = cell2mat(data(:,is_data_var == 1));

    is_data_var = cell2mat(cellfun(@(x) strcmp(x, 'Trmmol'), data_var, 'UniformOutput', 0));
    E = cell2mat(data(:,is_data_var == 1));

    is_data_var = cell2mat(cellfun(@(x) strcmp(x, 'VpdL'), data_var, 'UniformOutput', 0));
    VPD_L = cell2mat(data(:,is_data_var == 1));

    P_atm = g_w.*VPD_L./E; %[kPa]
    
    V_cmax25_func = @(psi_L_MD) 0;
    [outputs_photo_param] = Photosynthesis_parameters_temperature_response(T_L, 0, V_cmax25_func);
    
    Gamma_star = outputs_photo_param.Gamma_star_vect;
    K_c = outputs_photo_param.K_c_vect;
    K_o = outputs_photo_param.K_o_vect;
    K_m = outputs_photo_param.K_m_vect;

    % remove bad c_i
    c_i(c_i < Gamma_star) = nan;
    c_i(c_i > c_a) = nan;

    R_d_per_R_d25 = outputs_photo_param.R_d_per_R_d25_vect;
    V_cmax_per_V_cmax25 = outputs_photo_param.V_cmax_per_V_cmax25_vect;
    R_d25_per_V_cmax25 = 0.015; %from Collatz et al. (1991)

    %% Estimate V_cmax by one-point method
    V_cmax = A_n ./ ((c_i - Gamma_star)./(c_i + K_m) - R_d_per_R_d25./V_cmax_per_V_cmax25*R_d25_per_V_cmax25); %mol m-2 s-1

    V_cmax25 = V_cmax ./ V_cmax_per_V_cmax25;
    V_cmax(V_cmax25 < 0) = nan; %ignore data with unrealistic Vcmax25
    V_cmax25(V_cmax25 > 4e-4) = nan; %ignore data with unrealistic Vcmax25
    V_cmax(isnan(V_cmax25)) = nan; %ignore data with unrealistic Vcmax25
    R_d25 = R_d25_per_V_cmax25 * V_cmax25;
    R_d = R_d_per_R_d25 .* R_d25;

    %% Estimate k = dA_n/dc_i
    k = V_cmax .* (Gamma_star + K_m) ./ (c_i + K_m).^2;

    %% Estimate lambda
    lambda = 1.6*k./(g_w + 1.6*k).*A_n./E;
    lambda(lambda < 0) = nan; 

    %% Other data
    is_data_var = cell2mat(cellfun(@(x) strcmp(x, 'Species'), data_var, 'UniformOutput', 0));
    species = data(:,is_data_var == 1);

    is_data_var = cell2mat(cellfun(@(x) strcmp(x, 'Mid_dayWaterPotential'), data_var, 'UniformOutput', 0));
    psi_L_MD = cell2mat(data(:,is_data_var == 1));

    is_data_var = cell2mat(cellfun(@(x) strcmp(x, 'Pre_dawnWaterPotential'), data_var, 'UniformOutput', 0));
    psi_L_PD = cell2mat(data(:,is_data_var == 1));

    % remove outliers from (log10-transformed) lambda
    species_unique = unique(species);
    N_species_unique = length(species_unique);
    for i = 1:N_species_unique
        is_species = cell2mat(cellfun(@(x) strcmp(x, species_unique{i}), species, 'UniformOutput', 0));
        log10_lambda_Q1_species = prctile(log10(lambda(is_species == 1)), 25);
        log10_lambda_Q3_species = prctile(log10(lambda(is_species == 1)), 75);
        log10_lambda_IQR_species = log10_lambda_Q3_species - log10_lambda_Q1_species;
        log10_lambda_LB_species = log10_lambda_Q1_species - 1.5*log10_lambda_IQR_species;
        log10_lambda_UB_species = log10_lambda_Q3_species + 1.5*log10_lambda_IQR_species;
        lambda((is_species == 1) .* (log10(lambda) < log10_lambda_LB_species) == 1) = nan;
        lambda((is_species == 1) .* (log10(lambda) > log10_lambda_UB_species) == 1) = nan;
    end

    % remove zero g_w estimates
    g_w(g_w == 0) = nan;

    %% Treatment combinations
    is_data_var = cell2mat(cellfun(@(x) strcmp(x, 'heatTrt'), data_var, 'UniformOutput', 0));
    HeatTrt = data(:,is_data_var == 1);

    is_data_var = cell2mat(cellfun(@(x) strcmp(x, 'WaterTrt'), data_var, 'UniformOutput', 0));
    WaterTrt = data(:,is_data_var == 1);

    [Trt_numb, Comb_Trt_name] = Assign_Trt_number(species, HeatTrt, WaterTrt); 
    
    %% Saturated water content
    SWC_L = nan(M_data, 1);
    
    is_aceru = cell2mat(cellfun(@(x) strcmp(x, 'aceru'), species, 'UniformOutput', 0));
    is_queru = cell2mat(cellfun(@(x) strcmp(x, 'queru'), species, 'UniformOutput', 0));
    SWC_L(is_aceru == 1) = 2.6; %based on Belluau et al. (2021; Functional Ecology)
    SWC_L(is_queru == 1) = 2.5; %based on Belluau et al. (2021; Functional Ecology)
    
    
    %% Save data as .mat
    save(Mat_data_save_file, 'T_L', 'A_n', 'c_i', 'c_a', 'g_w', 'E', ...
         'VPD_L', 'P_atm', 'Gamma_star', 'K_c', 'K_o', 'K_m', 'V_cmax', ...
         'R_d', 'lambda', 'psi_L_MD', 'psi_L_PD', 'species', 'Trt_numb', ...
         'Comb_Trt_name', 'SWC_L');

end

%% Store Outputs
if isempty(Var)
    Var.T_L = T_L;
    Var.A_n = A_n;
    Var.c_i = c_i;
    Var.c_a = c_a; 
    Var.g_w = g_w;
    Var.E = E;
    Var.VPD_L = VPD_L; 
    Var.P_atm = P_atm; 
    Var.Gamma_star = Gamma_star;
    Var.K_c = K_c;
    Var.K_o = K_o;
    Var.K_m = K_m;
    Var.V_cmax = V_cmax;
    Var.R_d = R_d;
    Var.lambda = lambda;
    Var.psi_L_MD = psi_L_MD;
    Var.psi_L_PD = psi_L_PD;
    Var.species = species; 
    Var.Trt_numb = Trt_numb;
    Var.Comb_Trt_name = Comb_Trt_name; 
    Var.SWC_L = SWC_L;
else
    Var.T_L = [Var.T_L; T_L];
    Var.A_n = [Var.A_n; A_n];
    Var.c_i = [Var.c_i; c_i];
    Var.c_a = [Var.c_a; c_a]; 
    Var.g_w = [Var.g_w; g_w];
    Var.E = [Var.E; E];
    Var.VPD_L = [Var.VPD_L; VPD_L]; 
    Var.P_atm = [Var.P_atm; P_atm]; 
    Var.Gamma_star = [Var.Gamma_star; Gamma_star];
    Var.K_c = [Var.K_c; K_c];
    Var.K_o = [Var.K_o; K_o];
    Var.K_m = [Var.K_m; K_m];
    Var.V_cmax = [Var.V_cmax; V_cmax];
    Var.R_d = [Var.R_d; R_d];
    Var.lambda = [Var.lambda; lambda];
    Var.psi_L_MD = [Var.psi_L_MD; psi_L_MD];
    Var.psi_L_PD = [Var.psi_L_PD; psi_L_PD];
    Var.species = [Var.species; species]; 
    Var.Trt_numb = [Var.Trt_numb; Trt_numb];
    Var.Comb_Trt_name = [Var.Comb_Trt_name; Comb_Trt_name]; 
    Var.SWC_L = [Var.SWC_L; SWC_L];
end


end

