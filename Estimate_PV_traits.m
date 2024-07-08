function [] = Estimate_PV_traits(species, Trt_numb, RWC_t, psi_L_MD, varargin)

X_max = inf; %default value for ''X_max''
X_min = -inf; %default value for ''X_max''

for i = 1:ceil(length(varargin)/2)
    switch varargin{1 + 2*(i-1)}
        case 'X_max'
            X_max = varargin{2*i};
        case 'X_min'
            X_min = varargin{2*i};
        otherwise
            if isstring(varargin{1 + 2*(i-1)}) || ischar(varargin{1 + 2*(i-1)})
                error(['ERROR: ''Estimate_PV_traits'' does not accept ''', varargin{1 + 2*(i-1)}, ''' as an input!'])
            else
                error('ERROR: Odd-numbered ''varargin'' inputs must be strings or characters!')
            end
    end
end

% first subset by species
species_unique = unique(species);
N_species_unique = length(species_unique);
for i = 1:N_species_unique
    is_species = cell2mat(cellfun(@(x) strcmp(x, species_unique{i}), species, 'UniformOutput', 0)); 
    
    % first subset by treatment
    Trt_numb_unique = unique(Trt_numb(is_species == 1));
    N_Trt_unique = length(Trt_numb_unique);
    for j = 1:N_Trt_unique
        
        is_Trt = (Trt_numb == j);
        is_species_and_Trt = is_species .* is_Trt;
        RWC_t_subset = RWC_t(is_species_and_Trt == 1);
        psi_L_MD_subset = psi_L_MD(is_species_and_Trt == 1);
        
        %% Determine turgor loss point
        [RWC_t_tlp, pi_L_tlp, pi_L_0_trad, a_f_trad] = Estimate_turgor_loss_point(RWC_t_subset, psi_L_MD_subset, 'X_max', X_max, 'X_min', X_min); 
        
        %% Plot SSE for ''a_f_star" and ''beta"
        %%%[SSE_test] = Estimate_a_f_star_and_beta(RWC_t_subset, psi_L_MD_subset, RWC_t_tlp, pi_L_tlp, pi_L_0_trad, a_f_trad);
        
    end
end

end

