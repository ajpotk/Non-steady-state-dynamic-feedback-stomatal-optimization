function [psi_soil_func, theta_max, a, n] = Pedotransfer_for_psi_as_func_theta(sand, silt, clay, rho_soil_bulk)
%% Input:
% sand -- sand content [%] (0-100%)
% silt -- silt content [%] (0-100%)
% clay -- clay content [%] (0-100%)
% rho_soil_bulk -- soil bulk density [g cm-3]
%% Output:
% psi_soil_func -- soil water potential [MPa] expressed as function of volumetric soil water content
% theta_max -- maximum volumetric soil water content [m3 m-3]


% check that sand + silt + clay = 100
if (sand + silt + clay) ~= 100
   error('ERROR: Sum of sand, silt, and clay content must equal 100%!') 
end

% Pedotransfer function for European forest soils from Teepe et al. (2003; J. Plant Nutr. Soil Sci.)
theta_max = 0.9786-0.36686*rho_soil_bulk;
if theta_max > 1
    error('ERROR: ''theta_max'' must not exceed 1!')
end
if theta_max <= 0
    error('ERROR: ''theta_max'' must be positive!')
end

a = exp(55.576 - 4.433*rho_soil_bulk - 0.002*silt^2 - 0.470*clay - 0.066*sand/rho_soil_bulk - 3.683*sand^0.5 - 0.0359*silt/rho_soil_bulk - 0.0016*sand^2 - 3.6916*silt^0.5 + 1.8643*log(sand) + 1.575*log(silt)); %[kPa-1]
a = max(a, 0.005); %a must be no smaller than the smallest value in Teepe et al. (2003)
a = 1e3 * a; %convert from [kPa-1] to [MPa-1]
n = 1 + exp(-2.8497 + 0.00027395*sand^2 + 0.01637*silt); %[-]
theta_min = 0;

psi_soil_func = @(theta) -((((theta - theta_min)/(theta_max - theta_min)).^(n/(1-n)) - 1).*(1/n))/a; %[MPa]
end

