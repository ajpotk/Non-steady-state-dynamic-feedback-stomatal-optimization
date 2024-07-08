function [Trt_numb, Comb_Trt_name] = Assign_Trt_number(species, varargin)

M_data = length(species);

if isempty(varargin)
    Trt_numb = ones(M_data, 1);
    Comb_Trt_name = repmat({'none'}, M_data, 1);
else
    
    Trt_numb = nan(M_data, 1);
    Comb_Trt_name = cell(M_data, 1);
    N_ax = length(varargin); %number of axes of variation

    % first subset by species
    species_unique = unique(species);
    N_species_unique = length(species_unique);
    for i = 1:N_species_unique

        is_species = cell2mat(cellfun(@(x) strcmp(x, species_unique{i}), species, 'UniformOutput', 0)); 
        ind_is_species = is_species .* reshape(1:M_data, size(is_species));
        ind_is_species = ind_is_species(ind_is_species > 0);
        M_data_subset = length(ind_is_species);
        
        % species-specific treatments
        Trt_subset = cell(1,N_ax);
        for j = 1:N_ax
            Trt_subset{j} = varargin{j}(ind_is_species);
        end

        % unique treatments and number of them
        Trt_subset_unique = cellfun(@(x) unique(x), Trt_subset, 'UniformOutput', 0);
        N_Trt_subset_unique = cell2mat(cellfun(@(x) length(x), Trt_subset_unique, 'UniformOutput', 0));
        N_Trt_subset_combo = prod(N_Trt_subset_unique); % number of possible combinations
        Trt_subset_combo = cell(N_Trt_subset_combo, N_ax);
        for j = 1:N_ax
            N_repmat = prod(N_Trt_subset_unique(1:j-1));
            N_repelem = prod(N_Trt_subset_unique(j+1:end));
            Trt_subset_combo(1:N_Trt_subset_combo,j) = repmat(repelem(Trt_subset_unique{j}, N_repelem, 1), N_repmat, 1);
        end
        
        % assign treatment identification number for subset
        Trt_numb_subset = nan(M_data_subset, 1);
        for j = 1:M_data_subset
            k = 0;
            while 1 
                k = k + 1;
                is_same = nan(1,N_ax);
                for l = 1:N_ax
                    is_same(l) = strcmp(Trt_subset{l}{j}, Trt_subset_combo{k,l});
                end
                
                if all(is_same)
                    Trt_numb_subset(j) = k;
                    break
                end
                
                if k == N_Trt_subset_combo
                    error('ERROR: Code could not assign Treatment number!')
                end
            end
        end
        
        % combined/total name for Trt_subset_combo
        Comb_Trt_name_subset_combo = cell(N_Trt_subset_combo, 1);
        Trt_subset_combo_alt = Trt_subset_combo;
        if N_ax > 1
            Trt_subset_combo_alt(:,2:end) = cellfun(@(x) [' - ', x], Trt_subset_combo_alt(:,2:end), 'UniformOutput', 0); 
        end
        for j = 1:N_Trt_subset_combo
            Comb_Trt_name_subset_combo{j} = cell2mat(Trt_subset_combo_alt(j,1:N_ax));
        end
        Comb_Trt_name_subset = Comb_Trt_name_subset_combo(Trt_numb_subset);
        
        % store ''Trt_numb_subset'' in ''Trt_numb''
        Trt_numb(ind_is_species) = Trt_numb_subset;
        Comb_Trt_name(ind_is_species) = Comb_Trt_name_subset;
        
    end
end
end

