function [output_Trt] = Stomata_and_PV_for_Trt(a_f_max, pi_L_star, beta, epsilon_L_max, SWC_L, alpha, ...
                                               pi_L_0_25C, psi_L_MD_Trt, VPD_L_Trt, ...
                                               P_atm_Trt, T_L_Trt, V_cmax_Trt, ...
                                               R_d_Trt, Gamma_star_Trt, K_m_Trt, ...
                                               c_a_Trt, ...
                                               varargin)

%% Stomata_and_PV_for_Trt
% written by Aaron Potkay 2023

%% Inputs:
%   a_f_star -- apoplasm fraction at reference osmotic potential, pi_L_star [m3 m-3]
%   beta -- expontent between apoplasm fraction, a_f, and symplasm relative water content, RWC_s [-]
%   epsilon_L_max -- theoretical maximum elastic modulus [MPa]
%   SWC_L -- saturated water content of leaf [kg kg-1]
%   alpha -- molar conversion constant from carbon assimilation to osmolytes [mol mol-1]
%   pi_L_0_25C -- osmotic potential at full hydration at a reference temperature of 25C [MPa]

%% Outputs:
%	RWC_s -- predicted symplast relative water content [m3 m-3]
%	pi_L -- predicted osmotic potential [MPa]
%	P_L -- predicted turgor [MPa]
%	a_f_0 -- predicted apoplasm fraction at full hydration [m3 m-3]
%	a_f -- predicted apoplasm fraction [m3 m-3]
%	RWC_t -- predicted total/bulk relative water content [m3 m-3]
%	epsilon_L_0 -- predicted elastic modulus at full hydration [MPa]
%	chi_w -- predicted marginal water-use efficiency [mol mol-1]
%   g_w -- stomatal conductance to H2O [mol m-2 s-2]
%   solved -- Boolean operator denoting whether solution is realistic

%% User-specified settings
DFC = 1; % diffusion fixation coupling -- diffusion and photosynthetic fixation are coupled if DFC == 1 or decoupled if DFC == 0
Constant_chi_w = nan; %default to nan

%% Check ''varargin''
n_varargin = length(varargin);
n_varargin_name_and_value = floor(n_varargin/2);
if 2*n_varargin_name_and_value ~= n_varargin
    error('ERROR: varargin inputs should have an even total number of inputs: each name must has a value!')
end

for i = 1:n_varargin_name_and_value
   
    switch varargin{1+2*(i-1)}
        case 'DFC'
            if varargin{2*i} == 0
                DFC = 0;
            elseif varargin{2*i} == 1
                DFC = 1;
            else
                error('ERROR: ''DFC'' must be either 0 or 1!')
            end
        case 'Constant_chi_w'
            Constant_chi_w = varargin{2*i};
        otherwise
            if isstring(varargin{1+2*(i-1)}) || ischar(varargin{1+2*(i-1)})
                error(['ERROR: unsupported varargin name: ''', varargin{1+2*(i-1)}, '''', newline, 'Acceptable varagin names are ''xlim'', ''ylim'', ''color'', ''linecolor'', and ''linewidth'''])
            else
                error('ERROR: odd-numbered varargin inputs must be strings or characters!')
            end
    end
    
end

%% Calculate Predictions

[outputs_PV] = PV_from_psi_for_Trt( a_f_max, pi_L_star, beta, epsilon_L_max, ...
                                    pi_L_0_25C, T_L_Trt, psi_L_MD_Trt);

RWC_s = outputs_PV.RWC_s;
RWC_s_tlp = outputs_PV.RWC_s_tlp;
RWC_t = outputs_PV.RWC_t;
pi_L = outputs_PV.pi_L;
pi_L_0 = outputs_PV.pi_L_0;
P_L = outputs_PV.P_L;
epsilon_L = outputs_PV.epsilon_L;
epsilon_L_0 = outputs_PV.epsilon_L_0;
a_f = outputs_PV.a_f;
a_f_0 = outputs_PV.a_f_0;
da_f_dpi_L = outputs_PV.da_f_dpi_L;
da_f_0dpi_L_0 = outputs_PV.da_f_0dpi_L_0;

% range of pi_L
pi_L_min = min(pi_L, [], 'all');
pi_L_max = max(pi_L, [], 'all');

[outputs_stomata] = Stomata_for_RWC_and_Temperature(SWC_L, alpha, RWC_t, RWC_s, pi_L, pi_L_0, ...
                                                    a_f, da_f_dpi_L, a_f_0, da_f_0dpi_L_0, epsilon_L_0, ...
                                                    VPD_L_Trt, P_atm_Trt, T_L_Trt, V_cmax_Trt, ...
                                                    R_d_Trt, Gamma_star_Trt, K_m_Trt, c_a_Trt, ...
                                                    'DFC', DFC, 'Constant_chi_w', Constant_chi_w);

%% Store outputs
output_Trt.RWC_s = RWC_s;
output_Trt.RWC_s_tlp = RWC_s_tlp;
output_Trt.RWC_t = RWC_t;
output_Trt.pi_L = pi_L;
output_Trt.pi_L_0 = pi_L_0;
output_Trt.P_L = P_L;
output_Trt.a_f_0 = a_f_0;
output_Trt.a_f = a_f;
output_Trt.epsilon_L = epsilon_L;
output_Trt.epsilon_L_0 = epsilon_L_0;
output_Trt.chi_w = outputs_stomata.chi_w;
output_Trt.g_w = outputs_stomata.g_w;
output_Trt.solved = outputs_stomata.solved;
output_Trt.pi_L_min = pi_L_min;
output_Trt.pi_L_max = pi_L_max;
output_Trt.chi_w = outputs_stomata.chi_w;

end

