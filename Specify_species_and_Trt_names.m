function [species_subset_plot, Trt_subset_combo_plot] = Specify_species_and_Trt_names(species_subset, Trt_subset_combo)

species_covar = {'Larix',   'European larch';...
                 'Picea',   'Norway spruce'; ...
                 'aceru',   'Red maple'; ...
                 'pinpond', 'Ponderosa pine'; ...
                 'queru',   'Red oak'};

Trt_covar = {'ambT - Ambient', ['+0', char(176), 'C - Ambient P']; ...
             'ambT - Shelter', ['+0', char(176), 'C - Shelter']; ...
             '1.7C - Ambient', ['+1.7', char(176), 'C - Ambient P']; ...
             '1.7C - Shelter', ['+1.7', char(176), 'C - Shelter']; ...
             '3.4C - Ambient', ['+3.4', char(176), 'C - Ambient P']; ...
             '3.4C - Shelter', ['+3.4', char(176), 'C - Shelter']; ...
             'CP',             'North Plateau population'; ... 
             'RMR',            'Rocky Mountain population'; ...
             'N13Ad_S1',       '1300 m ASL - Dry'; ...
             'N13Ad_S2',       '1300 m ASL - Dry'; ...    
             'N13Bd_L1',       '1300 m ASL - Dry'; ...
             'N13Bd_L2',       '1300 m ASL - Dry'; ...
             'N13WAd_L1',      '1300 m ASL - Wet'; ...
             'N13WAd_S1',      '1300 m ASL - Wet'; ... 
             'N13WAd_S2',      '1300 m ASL - Wet'; ...
             'N13WBd_L2',      '1300 m ASL - Wet'; ...
             'N13WBd_L3',      '1300 m ASL - Wet'; ...
             'N13WBd_S3',      '1300 m ASL - Wet'; ...
             'S22Ad_L1',       '2200 m ASL'; ...
             'S22Ad_L2',       '2200 m ASL'};

is_species_subset = cell2mat(cellfun(@(x) strcmp(species_subset, x), species_covar(:,1), 'UniformOutput', 0));
if any(is_species_subset)
    ind_species_subset = find(is_species_subset);
    species_subset_plot = species_covar{ind_species_subset, 2}; 
else
    species_subset_plot = ['??? ', species_subset, ' ???'];
end

N_Trt_subset_combo = length(Trt_subset_combo);
Trt_subset_combo_plot = cell(1, N_Trt_subset_combo);
for i = 1:N_Trt_subset_combo
    is_Trt_subset_combo = cell2mat(cellfun(@(x) strcmp(Trt_subset_combo{i}, x), Trt_covar(:,1), 'UniformOutput', 0));
    if any(is_Trt_subset_combo)
        ind_Trt_subset_combo = find(is_Trt_subset_combo);
        Trt_subset_combo_plot{i} = Trt_covar{ind_Trt_subset_combo, 2}; 
    else
        Trt_subset_combo_plot{i} = ['??? ', Trt_subset_combo{i}, ' ???'];
    end
end

end

