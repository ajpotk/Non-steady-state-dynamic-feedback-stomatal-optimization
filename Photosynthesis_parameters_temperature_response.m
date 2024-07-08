function [outputs_photo_param] = Photosynthesis_parameters_temperature_response(T_L_vect, psi_L_vect, V_cmax25_func, varargin)

% default settings
Nonmonotonic_V_cmax_temperature_response = 0;
Electron_transport_limitation = 0;
R_abs = inf; %absorbed shortwave radiation [W m-2]

% check ''varargin''
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
                error('ERROR: ''Nonmonotonic_V_cmax_temperature_response'' must be a numeric value set to either 0 or 1!')
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
        case 'R_abs'
            R_abs = varargin{2*i};
            if ~isnumeric(R_abs)
                error('ERROR: ''Electron_transport_limitation'' must be a numeric value set to either 0 or 1!')
            end
        otherwise
            if isstring(varargin{1+2*(i-1)}) || ischar(varargin{1+2*(i-1)})
                error(['ERROR: unsupported varargin name: ''', varargin{1+2*(i-1)}, '''', newline, 'Acceptable varargin names are ''Nonmonotonic_f0_temperature_response'', ''Electron_transport_limitation'', and ''R_abs''.'])
            else
                error('ERROR: odd-numbered varargin inputs must be strings or characters!')
            end
    end
    
end

% local parameters
o_i = 1e-3 * 210; %standard atmospheric partial pressure of oxygen [mol mol-1]
R = 8.314; %universal gas constant [J mol-1 K-1]
R_d25_per_V_cmax25 = 0.015; %from Collatz et al. (1991)
J_max25_per_V_cmax25 = 1.67; %from Medlyn et al. (2002)

% local temperature-dependent photosynthetic parameters -- from Bernacchi et al. (2001)
Gamma_star_vect = 1e-6 * 42.75*exp(37830*(T_L_vect-25)/298.15/R./(T_L_vect+273.15)); %[mol mol-1]
K_c_vect = 1e-6 * 404.9*exp(79430*(T_L_vect-25)/298.15/R./(T_L_vect+273.15)); %[mol mol-1]
K_o_vect = 1e-3 * 278.4*exp(36380*(T_L_vect-25)/298.15/R./(T_L_vect+273.15)); %[mol mol-1]
K_m_vect = K_c_vect.*(1 + o_i./K_o_vect); %[mol mol-1]

V_cmax25_vect = V_cmax25_func(psi_L_vect);
V_cmax_per_V_cmax25_vect = exp(65330*(T_L_vect-25)/298.15/R./(T_L_vect+273.15));
if Nonmonotonic_V_cmax_temperature_response == 1
    V_cmax_per_V_cmax25_vect = V_cmax_per_V_cmax25_vect .* (1 + exp((298.15*640 - 2e5)/R/298.15)) ./ (1 + exp(((T_L_vect+273.15)*640 - 2e5)/R./(T_L_vect+273.15))); %modify to include temperature optimum at 36.4C
end
V_cmax_vect = V_cmax_per_V_cmax25_vect .* V_cmax25_vect;

R_d25_vect = R_d25_per_V_cmax25 * V_cmax25_vect; 
R_d_per_R_d25_vect = exp(46390*(T_L_vect-25)/298.15/R./(T_L_vect+273.15));
R_d_vect = R_d_per_R_d25_vect .* R_d25_vect;

%store outputs
outputs_photo_param.Gamma_star_vect = Gamma_star_vect;
outputs_photo_param.K_c_vect = K_c_vect;
outputs_photo_param.K_o_vect = K_o_vect;
outputs_photo_param.K_m_vect = K_m_vect;
outputs_photo_param.V_cmax25_vect = V_cmax25_vect;
outputs_photo_param.V_cmax_per_V_cmax25_vect = V_cmax_per_V_cmax25_vect;
outputs_photo_param.V_cmax_vect = V_cmax_vect;
outputs_photo_param.R_d25_vect = R_d25_vect;
outputs_photo_param.R_d_per_R_d25_vect = R_d_per_R_d25_vect;
outputs_photo_param.R_d_vect = R_d_vect;

if Electron_transport_limitation == 1
    % more local temperature-dependent photosynthetic parameters -- from Bernacchi et al. (2003)
    J_max25_vect = J_max25_per_V_cmax25 * V_cmax25_vect; 
    J_max_per_J_max25_vect = exp(43900*(T_L_vect-25)/298.15/R./(T_L_vect+273.15));
    J_max_vect = J_max_per_J_max25_vect .* J_max25_vect; 
    
    PAR = 2.1e-6 * R_abs; %convert R_abs [W m-2] to PAR [mol m-2 s-1]
    J_vect = min(0.3*PAR, J_max_vect); %assuming a constant quantum yield of electron transport of 0.3 mol electrons mol−1 photon (Medlyn et al., 2002)
    
    %store more outputs
    outputs_photo_param.J_max25_vect = J_max25_vect;
    outputs_photo_param.J_max_per_J_max25_vect = J_max_per_J_max25_vect;
    outputs_photo_param.J_max_vect = J_max_vect;
    outputs_photo_param.PAR = PAR;
    outputs_photo_param.J_vect = J_vect;
end



end

