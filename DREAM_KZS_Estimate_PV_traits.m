clear
clc
close all

%% DREAM_Estimate_PV_traits_for_B4W
% written by Aaron Potkay 2023

[Var] = Download_B4W_data();
[Var] = Download_Sapes_data(Var); 
[Var] = Download_Lotschental_data(Var); 

% retrieve structured data
T_L = Var.T_L;
K_c = Var.K_c;
K_o = Var.K_o;
K_m = Var.K_m;
Gamma_star = Var.Gamma_star;
c_a = Var.c_a;
c_i = Var.c_i;
V_cmax = Var.V_cmax;
R_d = Var.R_d;
A_n = Var.A_n;
E = Var.E;
g_w = Var.g_w;
lambda = Var.lambda;
psi_L_MD = Var.psi_L_MD;
psi_L_PD = Var.psi_L_PD;
VPD_L = Var.VPD_L;
P_atm = Var.P_atm; 
Trt_numb = Var.Trt_numb;
Comb_Trt_name = Var.Comb_Trt_name;
SWC_L = Var.SWC_L;
species = Var.species;

%% Subset data by species
species_unique = unique(species);
N_species_unique = length(species_unique);
M_data = length(g_w);

% identify nans
is_not_nan = (~isnan(g_w)) .* ...
             (~isnan(T_L)) .* ...
             (~isnan(A_n)) .* ...
             (~isnan(lambda)) .* ...
             (~isnan(psi_L_MD)) .* ...
             (~isnan(V_cmax));

for i = 1:N_species_unique
    
    species_subset = species_unique{i}; %choice of species
    disp(['Beginning DREAM KZS for species #', num2str(i), ' of ' num2str(N_species_unique), ': ', species_subset])
    disp('    Subsetting data by species')
    is_species = cell2mat(cellfun(@(x) strcmp(x, species_subset), species, 'UniformOutput', 0)); 
    
    is_species_and_not_nan = is_species .* is_not_nan; 
    ind_is_species_and_not_nan = is_species_and_not_nan .* reshape(1:M_data, size(is_species));
    ind_is_species_and_not_nan = ind_is_species_and_not_nan(ind_is_species_and_not_nan > 0);
    M_data_species = length(ind_is_species_and_not_nan);

    % species-specific data
    T_L_subset = T_L(ind_is_species_and_not_nan);
    K_c_subset = K_c(ind_is_species_and_not_nan);
    K_o_subset = K_o(ind_is_species_and_not_nan);
    K_m_subset = K_m(ind_is_species_and_not_nan);
    Gamma_star_subset = Gamma_star(ind_is_species_and_not_nan);
    c_a_subset = c_a(ind_is_species_and_not_nan);
    c_i_subset = c_i(ind_is_species_and_not_nan);
    V_cmax_subset = V_cmax(ind_is_species_and_not_nan);
    R_d_subset = R_d(ind_is_species_and_not_nan);
    A_n_subset = A_n(ind_is_species_and_not_nan);
    E_subset = E(ind_is_species_and_not_nan);
    g_w_subset = g_w(ind_is_species_and_not_nan);
    lambda_subset = lambda(ind_is_species_and_not_nan);
    psi_L_MD_subset = psi_L_MD(ind_is_species_and_not_nan);
    psi_L_PD_subset = psi_L_PD(ind_is_species_and_not_nan);
    VPD_L_subset = VPD_L(ind_is_species_and_not_nan);
    P_atm_subset = P_atm(ind_is_species_and_not_nan); 
    Trt_numb_subset = Trt_numb(ind_is_species_and_not_nan);
    Comb_Trt_name_subset = Comb_Trt_name(ind_is_species_and_not_nan);
    SWC_L_subset = SWC_L(ind_is_species_and_not_nan);

    % save local .mat file for subset of empirical data
    load_data_file_name = 'data_subset.mat';
    save(load_data_file_name, ...
         'T_L_subset', 'K_c_subset', 'K_o_subset', 'K_m_subset', ...
         'Gamma_star_subset', 'c_a_subset', 'c_i_subset', ...
         'V_cmax_subset', 'R_d_subset', 'A_n_subset', 'E_subset',  ...
         'g_w_subset', 'lambda_subset', 'psi_L_MD_subset', ...
         'psi_L_PD_subset', 'VPD_L_subset', 'P_atm_subset', ...
         'Trt_numb_subset', 'Comb_Trt_name_subset',  'SWC_L_subset')
     
    %% Delete past Output file (if any; see "forwardmodel.m")
    Output_file = 'Output_text.txt';
    files = struct2cell(dir);
    files = files(1,:);
    is_output_file = cell2mat(cellfun(@(x) strcmp(x, Output_file), files, 'UniformOutput', 0));
    if any(is_output_file)
       delete(Output_file)
    end
    
    %% Run DREAM KZS
    
    DREAM_results_file = ['DREAM_KZS_results_for_', species_subset, '.mat'];
    files = struct2cell(dir);
    files = files(1,:);
    is_DREAM_results_file = cell2mat(cellfun(@(x) strcmp(x, DREAM_results_file), files, 'UniformOutput', 0));
    
    if any(is_DREAM_results_file)
        
        disp('    Loading previously-saved species-specific data')
        load(DREAM_results_file)
        
    else
    
        disp('    Running DREAM KZS for species-specific data')
        % Parameters:
        %   x(1) = a_f_max -- maximum apoplasm fraction at an osmotic potential of zero [m3 m-3]
        %   x(2) = pi_L_star -- osmotic potential at which a_f = exp(-1)*a_f_max = 0.37*a_f_max [MPa]
        %   x(3) = beta -- expontent in Weibull function describing apoplasm fraction, a_f [-]
        %   x(4) = epsilon_L_max -- theoretical maximum elastic modulus [MPa]
        %   x(5) = SWC_L -- saturated water content of leaf [kg kg-1]
        %   x(6) = alpha -- molar conversion constant from carbon assimilation to osmolytes [mol mol-1]

        N = 20;                         % Number of parallel chains in MCMC (default 4 or 20)
        T = 1e4;                        % Number of iterations (default 2000 or 5000)
        Nx = 6;                         % Dimension of model parameters (number of model parameters)
        Ny = M_data_species;            % Dimension of model responses (number of model outputs and number of observations)
        t1 = 100;                       % After which iteration the Kalman proposal is used (default 80 or 100)
        t2 = ceil(0.3*T);               % After which iteration the Kalman proposal is not used (default ceil(0.3*T))
        Ne = 200;                       % Number of archive samples for the Kalman proposal (default 100 or 300)
        p_k = 0.3;                      % The probability of using the Kalman proposal distribution (default 0.3)

        xmin = [0, -20,   0, 5,  1.5, 0];   % Lower bounds of the parameters
        xmax = [1, -1e-2, 4, 40, 3.5, 3];	% Upper bounds of the parameters
        range = [xmin' xmax'];              % Range of the parameters

        Obs = g_w_subset;           % The measurements
        sd = 0.35 * g_w_subset;     % The standard deviation of measurement errors (set to [] if unknown)

        Z  = dream_kzs(N,T,Nx,Ny,Obs,sd,range,p_k,t1,t2,Ne);
        fid_x = fopen('x.bin'); x_kzs = fread(fid_x,[Nx inf],'double'); fclose(fid_x); delete x.bin y.bin p.bin

        delete DREAM_KZS.mat  % The intermediate results saved when running dream_kzs

        chain_kzs = x_kzs';
        x_kzs =  permute(reshape(chain_kzs',[Nx,N,T]),[3,1,2]);

        save(DREAM_results_file)
        
        %% Rename Output file (see "forwardmodel.m")
        DREAM_Output_file = ['Output_file_for_', species_subset, '.txt'];
        movefile(Output_file, DREAM_Output_file)
        
    end
    
    disp(['Completed DREAM KZS for species #', num2str(i), ' of ' num2str(N_species_unique), ': ', species_subset])
    
end
