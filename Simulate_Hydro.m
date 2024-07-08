function [outputs_Simulate_Hydro] = ...
Simulate_Hydro(SWC_L, alpha, a_f_max, pi_L_star, beta, epsilon_L_max, pi_L_0_25C, ...
               k_L_max_25C, A, Q_10_k_L, V_cmax25_func, ...
               theta_initial, theta_max, psi_soil_func, Z_r, LAI, emiss_L, ...
               Date, P, T_a, RH, P_atm, R_0, g_H_a, c_a, dt, ...
               varargin)

% SWC_L -- leaf saturated water content [kg kg-1]
% alpha -- number of moles of carbon in mole of first product of photosynthesis [mol mol-1]
% a_f_max -- maximum apoplastic fraction [-]
% pi_L_star -- shape parameter for relationship between apoplastic fraction and osmotic potential [MPa]
% beta -- shape parameter for relationship between apoplastic fraction and osmotic potential [-]
% epsilon_L_max -- maximum elastic modulus [MPa]
% pi_L_0_25C -- osmotic potential at full hydration at standardized temperature of 25C [MPa]
% k_L_max_25C -- maximum soil-plant conductance per unit leaf area at 25C [mol m-2 s-1 MPa-1]
% A -- exponent for simple soil-plant conductance [MPa-1] -- k_L = k_L_max * exp(A*psi_L)
% Q_10_k_L -- Q10 value for the temperature-dependence of the maximum soil-plant conductance per unit leaf area [-]
% V_cmax25_func -- maximum carboxylation capacity at 25C [mol m-2 s-1] expressed as function of leaf water potential
% theta_initial -- initial volumetric soil water content [m3 m-3]
% theta_max -- maximum volumetric soil water content (porosity) [m3 m-3]
% psi_soil_func -- soil water potential [MPa] expressed as function of volumetric soil water content
% Z_r -- root depth [m]
% LAI -- time series of leaf area index [m2 m-2]
% P -- time series of Precipitation values [m s-1]
% T_a -- time series of air temperature [C]
% RH -- time series of relative humidity [-]
% P_atm -- time series of atmospheric pressure [kPa]
% c_a -- time series of air CO2 partial pressures [mol mol-1]
% R_0 -- time series of shortwave radiation above the canopy [W m-2]
% g_H_a -- time series of boundary layer heat conductance [mol m-2 s-1]
% emiss_L -- leaf emissivity [-]

%% Default settings
Nonmonotonic_V_cmax_temperature_response = 1;
Electron_transport_limitation = 1;
Sperry_model = 0;
save_file_extension = '';
theta_store = nan;
g_soil = 0; %soil vapor conductance [mol m-2 s-1] -- per ground area basis


%% Check ''varargin''
n_varargin = length(varargin);
n_varargin_name_and_value = floor(n_varargin/2);
if 2*n_varargin_name_and_value ~= n_varargin
    error('ERROR: varargin inputs should have an even total number of inputs: each name must has a value!')
end

for i = 1:n_varargin_name_and_value
   
    switch varargin{1+2*(i-1)}
        case 'Nonmonotonic_V_cmax_temperature_response'
            Nonmonotonic_V_cmax_temperature_response = varargin{2*i};
            if isnumeric(Nonmonotonic_V_cmax_temperature_response)
                if (Nonmonotonic_V_cmax_temperature_response ~=0) && (Nonmonotonic_V_cmax_temperature_response ~= 1)
                    error('ERROR: ''Nonmonotonic_V_cmax_temperature_response'' must be a numeric value set to either 0 or 1!')
                end
            else
                error('ERROR: ''Nonmonotonic_f0_temperature_response'' must be a numeric value set to either 0 or 1!')
            end
        case 'Electron_transport_limitation'
            Electron_transport_limitation = varargin{2*i};
            if isnumeric(Electron_transport_limitation)
                if (Electron_transport_limitation ~=0) && (Electron_transport_limitation ~= 1)
                    error('ERROR: ''Electron_transport_limitation'' must be a numeric value set to either 0 or 1!')
                end
            else
                error('ERROR: ''Electron_transport_limitation'' must be a numeric value set to either 0 or 1!')
            end
        case 'Sperry_model'
            Sperry_model = varargin{2*i};
            if isnumeric(Sperry_model)
                if (Sperry_model ~=0) && (Sperry_model ~= 1)
                    error('ERROR: ''Sperry_model'' must be a numeric value set to either 0 or 1!')
                end
            else
                error('ERROR: ''Sperry_model'' must be a numeric value set to either 0 or 1!')
            end
        case 'save_file_extension'
            save_file_extension = varargin{2*i};
            if ~ischar(save_file_extension)
                error('ERROR: ''save_file_extension'' must be a string!')
            end
        case 'g_soil'
            g_soil = varargin{2*i};
        case 'theta_store'
            theta_store = varargin{2*i};
        otherwise
            if isstring(varargin{1+2*(i-1)}) || ischar(varargin{1+2*(i-1)})
                error(['ERROR: unsupported varargin name: ''', varargin{1+2*(i-1)}, '''', newline, 'Acceptable varargin names are ''Nonmonotonic_V_cmax_temperature_response'', ''Electron_transport_limitation'', ''save_file_extension'', ''g_soil'', and ''theta_store''.'])
            else
                error('ERROR: odd-numbered varargin inputs must be strings or characters!')
            end
    end
    
end

%% Local save file
Mat_data_save_file = 'Saved_Simulate_hydro_data_temporary';
Mat_data_save_file = [Mat_data_save_file, save_file_extension]; 
if Sperry_model == 1
    Mat_data_save_file = [Mat_data_save_file, '_Sperry']; 
end
Mat_data_save_file = [Mat_data_save_file, '.mat'];

%% Local constants
R = 8.314; %universal gas constant [J mol-1 K-1]
m_w = 18e-3; %molar mass of water [kg mol-1]
rho_w = 997; %density of water [kg m-3]
w = 1e-6; %conversion factor from pascal to Megapascal [MPa Pa-1]
k_LAI = 0.5; %extinction coefficient for Beer's Law [-]
fR_L_0 = 0.3; %fraction of soil-to_leaf hydraulic resistance in leaves [-] at reference LAI
LAI_0 = 2.5; %reference LAI [m2 m-2]
k_L_max_25C_0 = k_L_max_25C; %maximum soil-plant conductance per unit leaf area at 25C [mol m-2 s-1 MPa-1] at reference LAI

% calculate LAI effect on k_L_max_25C
k_L_max_25C = k_L_max_25C_0./(fR_L_0 + (1-fR_L_0)*LAI/LAI_0);

% canopy-averaged absorbed radiation
R_abs = R_0.*(1-exp(-k_LAI*LAI))./LAI;
R_abs(LAI == 0) = 0; 

% calculalte saturated vapor pressure
e_a_sat = 0.61078 .* exp(17.27 * T_a ./ (T_a + 237.3)); %air vapor pressure [kPa] -- Teten's equation

% time discretization
N_t = length(Date);

% initial conditions
i_start = 1;
theta = theta_initial;
psi_soil = psi_soil_func(theta);

% storage vectors
if isnan(theta_store)
    outputs_Simulate_Hydro.theta_store = nan(1, N_t);
elseif length(theta_store) == N_t
    outputs_Simulate_Hydro.theta_store = theta_store;
else
    error('ERROR: ''theta_store'' must be a vector of ''N_t'' length!')
end
outputs_Simulate_Hydro.g_w_store = nan(1, N_t); 
outputs_Simulate_Hydro.chi_w_store = nan(1, N_t); 
outputs_Simulate_Hydro.E_store = nan(1, N_t); 
outputs_Simulate_Hydro.A_n_store = nan(1, N_t); 
outputs_Simulate_Hydro.c_i_store = nan(1, N_t); 
outputs_Simulate_Hydro.V_cmax_store = nan(1, N_t); 
outputs_Simulate_Hydro.J_store = nan(1, N_t);
outputs_Simulate_Hydro.J_max_store = nan(1, N_t); 
outputs_Simulate_Hydro.A_c_limited_store = nan(1, N_t); 
outputs_Simulate_Hydro.PLC_store = nan(1, N_t); 
outputs_Simulate_Hydro.psi_L_store = nan(1, N_t); 
outputs_Simulate_Hydro.psi_L_crit_store = nan(1, N_t); 
outputs_Simulate_Hydro.E_soil_store = nan(1, N_t); 

% progress tracking
dperc_disp = 1; %[%]
perc_disp = dperc_disp; %[%]

% temporary file saving
dperc_save = 5; %[%]
perc_save = dperc_save; %[%]

%% Check for previous temporary save file
files = struct2cell(dir);
files = files(1,:);
is_Mat_data_save_file = cell2mat(cellfun(@(x) strcmp(x, Mat_data_save_file), files, 'UniformOutput', 0));

if any(is_Mat_data_save_file)

    load(Mat_data_save_file, 'outputs_Simulate_Hydro')
    i_start = find(isnan(outputs_Simulate_Hydro.theta_store), 1, 'first');
    theta = outputs_Simulate_Hydro.theta_store(i_start - 1);
    perc = 100 * i_start/N_t;
    perc_disp = dperc_disp*(1 + floor(perc/dperc_disp));
    perc_save = dperc_save*(1 + floor(perc/dperc_save));
    
end

%% Dynamic simulations
for i = i_start:N_t
    
    % track progress
    perc = 100 * i/N_t;
    if perc > perc_disp
         if Sperry_model == 1
            disp(['Simulate_Hydro with Sperry: ' num2str(floor(perc)), '% done!'])
         else
             disp(['Simulate_Hydro: ' num2str(floor(perc)), '% done!'])
         end
        perc_disp = dperc_disp*(1 + floor(perc/dperc_disp));
    end
    
    % regularly save progress
    if perc > perc_save
        save(Mat_data_save_file, 'outputs_Simulate_Hydro')
        perc_save = dperc_save*(1 + floor(perc/dperc_save));
    end
    
    
    % simulate gas exchange
    if Sperry_model == 1
        
        [outputs_stomata_Sperry] = ...
        Stomata_Sperry_Gain_Risk(k_L_max_25C(i), A, Q_10_k_L, psi_soil, ...
                                 T_a(i), RH(i), P_atm(i), c_a(i), ...
                                 R_abs(i), g_H_a(i), emiss_L, ...
                                 V_cmax25_func, ...
                                 'Nonmonotonic_V_cmax_temperature_response', Nonmonotonic_V_cmax_temperature_response, ...
                                 'Electron_transport_limitation', Electron_transport_limitation);
        
        outputs_Simulate_Hydro.g_w_store(i) = outputs_stomata_Sperry.g_w_Sperry;
        outputs_Simulate_Hydro.chi_w_store(i) = outputs_stomata_Sperry.chi_w_Sperry;
        outputs_Simulate_Hydro.E_store(i) = outputs_stomata_Sperry.E_Sperry;
        outputs_Simulate_Hydro.A_n_store(i) = outputs_stomata_Sperry.A_n_Sperry;
        outputs_Simulate_Hydro.c_i_store(i) = outputs_stomata_Sperry.c_i_Sperry;
        outputs_Simulate_Hydro.V_cmax_store(i) = outputs_stomata_Sperry.V_cmax_Sperry;
        outputs_Simulate_Hydro.J_store(i) = outputs_stomata_Sperry.J_Sperry;
        outputs_Simulate_Hydro.J_max_store(i) = outputs_stomata_Sperry.J_max_Sperry;
        outputs_Simulate_Hydro.A_c_limited_store(i) = outputs_stomata_Sperry.A_c_limited_Sperry;
        outputs_Simulate_Hydro.PLC_store(i) = outputs_stomata_Sperry.PLC_Sperry;
        outputs_Simulate_Hydro.psi_L_store(i) = outputs_stomata_Sperry.psi_L_Sperry;
        
    else

        [outputs_stomata_steady_state] = ...
        Stomata_for_steady_state_thermo_and_hydraulics(SWC_L, alpha, a_f_max, pi_L_star, beta, epsilon_L_max, pi_L_0_25C, ...
                                                       k_L_max_25C(i), A, Q_10_k_L, psi_soil, ...
                                                       T_a(i), RH(i), P_atm(i), c_a(i), ...
                                                       R_abs(i), g_H_a(i), emiss_L, ...
                                                       V_cmax25_func, ...
                                                       'Nonmonotonic_V_cmax_temperature_response', Nonmonotonic_V_cmax_temperature_response, ...
                                                       'Electron_transport_limitation', Electron_transport_limitation);


        outputs_Simulate_Hydro.g_w_store(i) = outputs_stomata_steady_state.g_w_steady_state;
        outputs_Simulate_Hydro.chi_w_store(i) = outputs_stomata_steady_state.chi_w_steady_state;
        outputs_Simulate_Hydro.E_store(i) = outputs_stomata_steady_state.E_steady_state;
        outputs_Simulate_Hydro.A_n_store(i) = outputs_stomata_steady_state.A_n_steady_state;
        outputs_Simulate_Hydro.c_i_store(i) = outputs_stomata_steady_state.c_i_steady_state;
        outputs_Simulate_Hydro.V_cmax_store(i) = outputs_stomata_steady_state.V_cmax_steady_state;
        outputs_Simulate_Hydro.J_store(i) = outputs_stomata_steady_state.J_steady_state;
        outputs_Simulate_Hydro.J_max_store(i) = outputs_stomata_steady_state.J_max_steady_state;
        outputs_Simulate_Hydro.A_c_limited_store(i) = outputs_stomata_steady_state.A_c_limited_steady_state;
        outputs_Simulate_Hydro.PLC_store(i) = outputs_stomata_steady_state.PLC_steady_state;
        outputs_Simulate_Hydro.psi_L_store(i) = outputs_stomata_steady_state.psi_L_steady_state;
        outputs_Simulate_Hydro.psi_L_crit_store(i) = outputs_stomata_steady_state.psi_L_crit_steady_state;

    end
    
    if any(isnan([outputs_Simulate_Hydro.g_w_store(i), outputs_Simulate_Hydro.chi_w_store(i), ...
                  outputs_Simulate_Hydro.E_store(i), outputs_Simulate_Hydro.A_n_store(i), ...
                  outputs_Simulate_Hydro.PLC_store(i), outputs_Simulate_Hydro.psi_L_store(i)]))
        error('ERROR: ''NaN''s found!')
    end
    
    % calculate soil evaporation
    RH_soil = exp(m_w*psi_soil/w/rho_w/R/(T_a(i)+273.15));
    E_soil = g_soil*e_a_sat(i)/P_atm(i)*(RH_soil - RH(i)); %soil evporation [mol m-2 s-1] -- per ground area basis
    E_soil(E_soil < 0) = 0;
    
    if isnan(theta_store)
        % update soil moisture
        dtheta_dt = (P(i)/dt - m_w/rho_w*(E_soil + LAI(i)*outputs_Simulate_Hydro.E_store(i)))/Z_r;
        theta = theta + dtheta_dt*dt;
        theta = min(theta, theta_max); 
        if theta < 0
            error('ERROR: Volumetric soil water content cannot go below zero!') 
        end
    else
        theta = theta_store(i); 
    end
    psi_soil = psi_soil_func(theta);
    
    outputs_Simulate_Hydro.theta_store(i) = theta;
    outputs_Simulate_Hydro.E_soil_store(i) = E_soil;
    
end

end

