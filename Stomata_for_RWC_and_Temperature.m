function [outputs_stomata] = Stomata_for_RWC_and_Temperature(SWC_L, alpha, RWC_t, RWC_s, pi_L, pi_L_0, ...
                                                             a_f, da_f_dpi_L, a_f_0, da_f_0dpi_L_0, epsilon_L_0, ...
                                                             VPD_L_Trt, P_atm_Trt, T_L_Trt, V_cmax_Trt, ...
                                                             R_d_Trt, Gamma_star_Trt, K_m_Trt, c_a_Trt, ...
                                                             varargin)

%% Stomata_for_RWC_and_Temperature
% written by Aaron Potkay 2023

%% Inputs:
%   SWC_L -- saturated water content of leaf [kg kg-1]
%   alpha -- molar conversion constant from carbon assimilation to osmolytes [mol mol-1]
%   RWC_s -- predicted symplast relative water content [m3 m-3]
%	RWC_t -- predicted total/bulk relative water content [m3 m-3]
%	pi_L -- predicted osmotic potential [MPa]
%   pi_L_0 -- osmotic potential at full hydration [MPa]
%	a_f -- predicted apoplasm fraction [m3 m-3]
%   da_fdpi_L -- partial derivative of a_f with respect to pi_L [MPa-1]
%	a_f_0 -- predicted apoplasm fraction at full hydration [m3 m-3]
%   da_f_0dpi_L_0 -- partial derivative of a_f_0 with respect to pi_L_0 [MPa-1]
%	epsilon_L_0 -- predicted elastic modulus at full hydration [MPa]

%% Outputs:
%	chi_w -- predicted marginal water-use efficiency [mol mol-1]
%   g_w -- stomatal conductance to H2O [mol m-2 s-2]
%   solved -- Boolean operator denoting whether solution is realistic

%% User-specified settings
DFC = 1; % diffusion fixation coupling -- diffusion and photosynthetic fixation are coupled if DFC == 1 or decoupled if DFC == 0
Constant_chi_w = nan; %default to nan

%% Constants
R = 8.314; %universal gas constant [J mol-1 K-1]
m_w = 18e-3; %molar mass of water [kg mol-1]
rho_w = 997; %density of water [kg m-3]
w = 1e-6; %conversion factor from pascal to Megapascal [MPa Pa-1]
lambda_E_0C = 2.5e6; %latent heat of vaporization at reference temperature of 0C [J kg-1]
dlambda_E_dT_L = -2365; %slope between latent heat of vaporization, lambda_E, and leaf temperature in degrees Celcius, T_L [K kg-1 C-1]
c_w = 4181; %specific heat capacity of water [J kg-1 K-1]
c_DM = 2814; %specific heat capacity of leaf dry matter [J kg-1 K-1] -- mean of 7 species from Jayalakshmy & Philip (2010, "Thermophysical Properties of Plant Leaves and Their Influence on the Environment Temperature", Int J Thermophys)

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
                error(['ERROR: unsupported varargin name: ''', varargin{1+2*(i-1)}, '''', newline, 'Acceptable varagin names are ''DFC'' and ''Constant_chi_w''.'])
            else
                error('ERROR: odd-numbered varargin inputs must be strings or characters!')
            end
    end
    
end

%% Calculate Predictions

lambda_E = lambda_E_0C + dlambda_E_dT_L * T_L_Trt;

if isnan(Constant_chi_w)
    
    % marginal carbon cost of water -- leaf osmotic potential as constraint
    term1 = SWC_L.*RWC_t./(SWC_L.*RWC_t*c_w + c_DM) .* lambda_E./(273.15+T_L_Trt);
    Delta = epsilon_L_0./(epsilon_L_0 - pi_L_0).*(pi_L_0./(1 - a_f_0).*da_f_0dpi_L_0 - 1);
    term2 = 1./(1 + Delta.*a_f);
    chi_w = -(1-a_f)*alpha*m_w.*pi_L/w/rho_w/R./(273.15+T_L_Trt).*(term1 - term2);

    % check partial derivative of chi_w with respect to RWC_s
    term3 = SWC_L*c_DM./(SWC_L.*RWC_t*c_w + c_DM).^2 .* lambda_E./(273.15+T_L_Trt) .* (1 - a_f_0 - pi_L./RWC_s.*da_f_dpi_L);
    term4 = Delta./(1 + Delta.*a_f).^2 .* pi_L./RWC_s.*da_f_dpi_L;
    dchi_wdRWC_s = chi_w./RWC_s.*(pi_L./(1-a_f).*da_f_dpi_L - 1) ...
                   - (1-a_f)*alpha*m_w.*pi_L/w/rho_w/R./(273.15+T_L_Trt).*(term3 - term4);
               
else
    
    chi_w = Constant_chi_w * ones(1, length(T_L_Trt));
	dchi_wdRWC_s = zeros(1, length(T_L_Trt));
    
end

L = 1.6*chi_w.*VPD_L_Trt./P_atm_Trt;
if DFC == 1
    
    % solve for optimal g_w with diffusion-fixation coupling
    
    r = R_d_Trt./V_cmax_Trt;
    z = L./(K_m_Trt + Gamma_star_Trt);
    y = L./(c_a_Trt - Gamma_star_Trt);
    X = ((1-y.*r) - ((1-y.*r).^2 - (1-z.*(1-r)).*(1-y+y.^2.*r./z)).^0.5)./(1-z.*(1-r)); % X = (c_i - Gamma_star)/(c_a - Gamma_star)
    g_w = 1.6./(c_a_Trt - Gamma_star_Trt)./(1-X).*(V_cmax_Trt.*X./(X+y./z) - R_d_Trt);
    g_c = g_w/1.6; 
    A_n = 0.5*(g_c.*(c_a_Trt + K_m_Trt) + V_cmax_Trt - R_d_Trt) - ...
          (0.25*(g_c.*(c_a_Trt - K_m_Trt) + R_d_Trt - V_cmax_Trt).^2 + g_c.*(V_cmax_Trt.*Gamma_star_Trt + R_d_Trt.*K_m_Trt + g_c.*c_a_Trt.*K_m_Trt)).^0.5;
    
elseif DFC == 0
    
    % solve for optimal g_w with diffusion-fixation decoupling
    g_w = 1.6./L.*(V_cmax_Trt.*(c_a_Trt - Gamma_star_Trt - L)./(c_a_Trt + K_m_Trt - L) - R_d_Trt);
    g_c = g_w/1.6; 
    A_n = 0.5*(g_c.*(c_a_Trt + K_m_Trt) + V_cmax_Trt - R_d_Trt) - ...
          (0.25*(g_c.*(c_a_Trt - K_m_Trt) + R_d_Trt - V_cmax_Trt).^2 + g_c.*(V_cmax_Trt.*Gamma_star_Trt + R_d_Trt.*K_m_Trt + g_c.*c_a_Trt.*K_m_Trt)).^0.5;
    c_i = c_a_Trt - A_n./g_c;
    X = (c_i - Gamma_star_Trt)./(c_a_Trt - Gamma_star_Trt);
    
else
    
    error('ERROR: ''DFC'' may be only either 0 or 1!')
    
end
E = g_w.*VPD_L_Trt./P_atm_Trt;

% check viability of solution
solved_real_X = (imag(X) == 0); %checks for real solution for g_w
if isnan(Constant_chi_w)
    solved_a_f_0 = (a_f_0 < 1); %checks that predicted a_f_0 is realistic
    solved = solved_a_f_0 .* solved_real_X;
else
    solved = solved_real_X;
end


%% Store outputs
outputs_stomata.chi_w = chi_w;
outputs_stomata.dchi_wdRWC_s = dchi_wdRWC_s;
outputs_stomata.g_w = g_w;
outputs_stomata.A_n = A_n;
outputs_stomata.E = E;
outputs_stomata.solved = solved;

end

