function [Var] = Download_Lotschental_data(varargin)

%% Figure colors
colors = [102, 102, 255;...
          124, 203, 161; ...
          240, 116, 110; ...
          253, 222, 156; ...
          002, 144, 153; ...
          220, 057, 120; ...
          006, 082, 117]/255;

%% Local variables
hour_day_start = 7; %[hour]
hour_day_end = 19; %[hour]
            
%% Check varargin
[Var] = Check_varargin_for_download(varargin);

%% Check for past data .mat file
Mat_data_save_file = 'Saved_Lotschental_data.mat';

files = struct2cell(dir);
files = files(1,:);
is_Mat_data_save_file = cell2mat(cellfun(@(x) strcmp(x, Mat_data_save_file), files, 'UniformOutput', 0));

if any(is_Mat_data_save_file)
    
    load(Mat_data_save_file)
    
else
    %% Constants
    o_i = 1e-3 * 210; %standard atmospheric partial pressure of oxygen [mol mol-1]
    R = 8.314; %universal gas constant [J mol-1 K-1]
    P_atm_0 = 101.325; %standard atmospheric pressure [kPa]
    g = 9.81; %gravitational acceleration [m s-2]
    gamma_T_a = 2.5e-3; %temperature lapse rate [K m-1]
    r_a = 287.04; %gas constant for dry air [J kg-1 K-1]

    %% Identification for trees
    % column 1: site/tree identification
    % column 2: site name
    % column 3: species
    Trt_name = {'N13Ad_S1',     'N13',  'Picea'; ...
                'N13Ad_S2',     'N13',  'Picea'; ...
                'N13Bd_L1',     'N13',  'Larix'; ...
                'N13Bd_L2',     'N13',  'Larix'; ...
                'N13WAd_L1',    'N13W', 'Larix'; ...
                'N13WAd_S1',    'N13W', 'Picea'; ...
                'N13WAd_S2',    'N13W', 'Picea'; ...
                'N13WBd_L2',    'N13W', 'Larix'; ...
                'N13WBd_L3',    'N13W', 'Larix'; ...
                'N13WBd_S3',    'N13W', 'Picea'; ...
                'S22Ad_L1',     'S22',  'Larix'; ...
                'S22Ad_L2',     'S22',  'Larix'};
    [N_Trt, ~] = size(Trt_name);
    species_unique = unique(Trt_name(:,3));
    N_species_unique = length(species_unique);
    max_Trt_numb_for_species = zeros(N_species_unique, 1);

    %% Data location
    home = pwd;

    mother_directory_data = 'C:\Users\potka002\Desktop\Research\Growth Optimizing Stomata\Test_GOH\Lotschental_data_2012-2015';

    % column 1: variable name, 
    % column 2: subfolder, 
    % column 3: begininng of text file, 
    % column 4: Boolean for whether text file name is followed by ''Trt_name'', 
    % column 5: ''formatSpec'' for textscan()
    % column 6: ''HeaderLines'' for textscan()
    % column 7: ''InputFormat'' for datetime() for date
    % column 8: ''InputFormat'' for datetime() for time of day
    daughter_directory_data = {'psi_L_MD',  '4. Leaf water potential',  'WPleaf_',              1,  '%s %s %s',                                 2,  'dd/MM/yyyy',   'HH:mm'; ...
                               'F_d',       '2. Sap flow',              'Fcrown_',              1,  '%s %s %s',                                 2,  'dd/MM/yyyy',   'HH:mm'; ...
                               'D',         '3. Dendrometer',           'D_stem_',              1,  '%s %s %s',                                 2,  'dd/MM/yyyy',   'HH:mm'; ...
                               'T_a',       '1. Environmental data',    'Temperature',          0,  '(%s %[^)] )%s %s %s %s %s %s %s %s %s',	1,  'MM/dd/yy',     'HH:mm:ss'; ...
                               'RH',        '1. Environmental data',    'Relative_humidity',    0,  '(%s %[^)] )%s %s %s %s %s %s %s %s %s',	1,  'MM/dd/yy',     'HH:mm:ss'};

    [N_daughter_data, ~] = size(daughter_directory_data);

    %% Download data by treatment
    % initialize data-storage structure
    Var = struct;
    for i = 1:N_daughter_data
        Var = setfield(Var, daughter_directory_data{i, 1}, []);
    end
    Var = setfield(Var, 'Trt_numb', []);
    Var = setfield(Var, 'Comb_Trt_name', []);
    Var = setfield(Var, 'species', []);
    Var = setfield(Var, 'Date', []);
    Var = setfield(Var, 'psi_L_PD', []);
    Var = setfield(Var, 'hour', []);


    disp('Downloading Lotschental data')
    temporary_save_name = 'Lotschental_temporary_data.mat';
    for i = 1:N_Trt

        disp(['Downloading Lotschental data for: ', Trt_name{i, 1}, ', #', num2str(i), '/', num2str(N_Trt)])
        Var_local = struct;

        for j = 1:N_daughter_data

            % change directoy to subfolder of data
            cd([mother_directory_data, '\', daughter_directory_data{j, 2}]);

            % read in raw data
            if daughter_directory_data{j, 4}
                file_name = [daughter_directory_data{j, 3}, Trt_name{i, 1}, '.txt'];
            else
                file_name = [daughter_directory_data{j, 3}, '.txt'];
            end

            disp(['     Downloading: ''', file_name, ''', ', num2str(j), '/', num2str(N_daughter_data)])
            fileID = fopen(file_name);
            raw_data = textscan(fileID, daughter_directory_data{j, 5}, 'HeaderLines', daughter_directory_data{j, 6});
            fclose(fileID);
            cd(home)

            % get time
            Date_raw = cellfun(@(x) datetime(x, 'InputFormat', daughter_directory_data{j,7}), raw_data{1}, 'UniformOutput', 0);
            Date_raw = cat(1, Date_raw{:}); %cell to array
            Timeofday_raw = cellfun(@(x) timeofday(datetime(x, 'InputFormat', daughter_directory_data{j,8})), raw_data{2}, 'UniformOutput', 0);
            Timeofday_raw = cat(1, Timeofday_raw{:}); %cell to array
            Datetime_raw = Date_raw + Timeofday_raw;

            % get measurements and convert to numbers
            raw_data_else = horzcat(raw_data{3:end});
            raw_data_else = cell2mat(cellfun(@(x) str2double(x), raw_data_else, 'UniformOutput', 0));
            [~, N_raw_data_else] = size(raw_data_else);

            if N_raw_data_else > 1

                % choose environmental conditions for correct site
                cd([mother_directory_data, '\', daughter_directory_data{j, 2}]);
                fileID = fopen(file_name);
                raw_Var_names = textscan(fileID, ['%s', repmat(' %s', 1, N_raw_data_else)], 1); %read just first line for variable names
                fclose(fileID);
                cd(home)

                raw_Var_names = raw_Var_names(2:end); %assuming first variable is ''Date''
                is_site = cell2mat(cellfun(@(x) strcmp(x, Trt_name{i, 2}), raw_Var_names, 'UniformOutput', 0));

                if ~any(is_site)
                    error('ERROR: Cannot match Lotschental data with proper site!')
                end

                ind_site = find(is_site, 1, 'first');
                raw_data_else = raw_data_else(:,ind_site);

            end

            % locally store data in ''Var_local'' structure
            Var_local = setfield(Var_local, daughter_directory_data{j, 1}, raw_data_else); 
            Var_local = setfield(Var_local, [daughter_directory_data{j, 1}, '_Datetime'], Datetime_raw); 

        end

        disp(['Collating Lotschental data for: ', Trt_name{i, 1}, ', #', num2str(i), '/', num2str(N_Trt)])
        for j = 2:N_daughter_data

            disp(['     Collating data arrays: 1 through ', num2str(j-1), ' with ', num2str(j)])

            if j == 2
                F1 = getfield(Var_local, daughter_directory_data{1, 1});
                Date1 = getfield(Var_local, [daughter_directory_data{1, 1}, '_Datetime']);
            else
                F1 = nan(N_collated, j-1);
                for k = 1:(j-1)
                    F1(1:N_collated, k) = getfield(Var_local, daughter_directory_data{k, 1});
                end
                Date1 = Date_collated;
            end
            F2 = getfield(Var_local, daughter_directory_data{j, 1});
            Date2 = getfield(Var_local, [daughter_directory_data{j, 1}, '_Datetime']);

            [ F1_collated, F2_collated, Date_collated ] = Collate_data( F1, Date1, F2, Date2, hours(1) );
            N_collated = length(Date_collated); %number of data points with sharred times

            for k = 1:(j-1)
                Var_local = setfield(Var_local, daughter_directory_data{k, 1}, F1_collated(:, k));
            end
            Var_local = setfield(Var_local, daughter_directory_data{j, 1}, F2_collated);

        end
        
        % identify leaf water potential as either midday or predawn by hour of day
        psi_L_PD_local = nan(size(Var_local.psi_L_MD));
        hour_local = hour(Date_collated);
        day_local = day(Date_collated);
        month_local = month(Date_collated);
        year_local = year(Date_collated);
        Date_local = datetime(year_local, month_local, day_local);
        Date_unique_local = unique(Date_local);
        N_Date_unique_local = length(Date_unique_local);
        for j = 1:N_Date_unique_local

            psi_L_MD_current_day_local = Var_local.psi_L_MD(Date_local == Date_unique_local(j));
            psi_L_MD_Datetime_current_day_hour_local = hour(Var_local.psi_L_MD_Datetime(Date_local == Date_unique_local(j)));

            psi_L_MD_current_day_morning_local = psi_L_MD_current_day_local(psi_L_MD_Datetime_current_day_hour_local <= hour_day_start);

            psi_L_MD_previous_day_local = Var_local.psi_L_MD(Date_local == (Date_unique_local(j) - day(1)));
            psi_L_MD_Datetime_previous_day_hour_local = hour(Var_local.psi_L_MD_Datetime(Date_local == (Date_unique_local(j) - day(1))));

            psi_L_MD_previous_day_night_local = psi_L_MD_previous_day_local(psi_L_MD_Datetime_previous_day_hour_local >= hour_day_end);

            psi_L_PD_current_day_local = [psi_L_MD_current_day_morning_local; psi_L_MD_previous_day_night_local]; 
            psi_L_PD_current_day_local = mean(psi_L_PD_current_day_local);

            psi_L_PD_local(Date_local == Date_unique_local(j)) = psi_L_PD_current_day_local;

        end

        % treatment identification number, treatment name, and species
        ind_species_unique = find(cell2mat(cellfun(@(x) strcmp(x, Trt_name{i, 3}), species_unique, 'UniformOutput', 0)), 1, 'first');
        max_Trt_numb_for_species(ind_species_unique) = max_Trt_numb_for_species(ind_species_unique) + 1;
        Trt_numb_local = max_Trt_numb_for_species(ind_species_unique)*ones(N_collated, 1);
        Comb_Trt_name_local = repmat(Trt_name(i, 1), N_collated, 1); 
        species_local = repmat(Trt_name(i, 3), N_collated, 1); 

        % store collated data in 'Var'
        for j = 1:N_daughter_data
            Var = setfield(Var, daughter_directory_data{j, 1}, [getfield(Var, daughter_directory_data{j, 1}); getfield(Var_local, daughter_directory_data{j, 1})]);
        end
        Var = setfield(Var, 'Trt_numb', [getfield(Var, 'Trt_numb'); Trt_numb_local]);
        Var = setfield(Var, 'Comb_Trt_name', [getfield(Var, 'Comb_Trt_name'); Comb_Trt_name_local]);
        Var = setfield(Var, 'species', [getfield(Var, 'species'); species_local]);
        Var = setfield(Var, 'Date', [getfield(Var, 'Date'); Date_collated]);
        Var = setfield(Var, 'psi_L_PD', [getfield(Var, 'psi_L_PD'); psi_L_PD_local]);
        Var = setfield(Var, 'hour', [getfield(Var, 'hour'); hour_local]);

        % save data temporarily as .mat
        save(temporary_save_name, 'Var')

    end

    %% local data variables
    psi_L_MD = Var.psi_L_MD;
    psi_L_PD = Var.psi_L_PD;
    F_d = Var.F_d;
    D = Var.D;
    T_a = Var.T_a;
    RH = Var.RH;
    Trt_numb = Var.Trt_numb;
    Comb_Trt_name = Var.Comb_Trt_name;
    species = Var.species;
    Date = Var.Date;
    hod = Var.hour;
    N_data = length(psi_L_MD);

    %% Download sapwood area
    daughter_directory_SWA = '5. Metadata';
    file_name_SWA = 'Metadata.txt';
    formatSpec_SWA = '%s %s %s %s %s %s %s %s %s %s %s %s %s %s (%s %[^)] %*[^\n]';
    HeaderLines_SWA = 1;

    % change directoy to subfolder of data
    cd([mother_directory_data, '\', daughter_directory_SWA]);
    fileID = fopen(file_name_SWA);
    raw_data = textscan(fileID, formatSpec_SWA, 'HeaderLines', HeaderLines_SWA);
    fclose(fileID);
    cd(home)

    Trt_name_metadata = raw_data{4};
    SWA_metadata = raw_data{13};
    SWA_metadata = cell2mat(cellfun(@(x) str2double(x), SWA_metadata, 'UniformOutput', 0));
    Elevation_metadata = raw_data{1};
    Elevation_metadata = cell2mat(cellfun(@(x) str2double(x(1:4)), Elevation_metadata, 'UniformOutput', 0));

    SWA = nan(N_data, 1);
    Elevation = nan(N_data, 1);
    for i = 1:N_Trt
        is_Trt = cell2mat(cellfun(@(x) strcmp(x, Trt_name{i, 1}), Trt_name_metadata, 'UniformOutput', 0));
        SWA(Trt_numb == i) = SWA_metadata(is_Trt == 1);
        Elevation(Trt_numb == i) = Elevation_metadata(is_Trt == 1);
    end

    % save data temporarily as .mat
    Var = setfield(Var, 'SWA', SWA);
    Var = setfield(Var, 'Elevation', Elevation);
    save(temporary_save_name, 'Var')

    %% Download V_cmax data
    daughter_directory_Aci = 'ACi\Loetschental_2022-07_ACi';
    file_name_Aci = 'ACi_Lotschental_selection.xlsx';
    sheet_name_Aci = 'Data';

    % change directoy to subfolder of data
    cd([mother_directory_data, '\', daughter_directory_Aci]);
    opts = detectImportOptions(file_name_Aci, 'Sheet', sheet_name_Aci);
    raw_data = readtable(file_name_Aci, opts);
    cd(home)

    Trt_name_Aci = table2cell(raw_data(:, 3));
    V_cmax_Aci = cell2mat(table2cell(raw_data(:, 5)));
    notes_Aci = table2cell(raw_data(:, 8));

    % remove 'D' from Trt_name_Aci to be comparable to names in other types of data
    ind_D_Trt_name_Aci = cellfun(@(x) strfind(x, 'D', 'ForceCellOutput', true), Trt_name_Aci, 'UniformOutput', 0);
    contains_D_Trt_name_Aci = cell2mat(cellfun(@(x) ~isempty(x{1}), ind_D_Trt_name_Aci, 'UniformOutput', 0));
    Trt_name_Aci(contains_D_Trt_name_Aci == 1) = cellfun(@(x,y) [x(1:y{1}-1), x(y{1}+1:end)], Trt_name_Aci(contains_D_Trt_name_Aci == 1), ind_D_Trt_name_Aci(contains_D_Trt_name_Aci == 1), 'UniformOutput', 0);

    % remove bad V_cmax measurements
    is_closed = cell2mat(cellfun(@(x) strcmp(x, 'Stomata appear closed'), notes_Aci, 'UniformOutput', 0));
    V_cmax_Aci(is_closed == 1) = nan;
    %%%is_unsat = cell2mat(cellfun(@(x) strcmp(x, 'Not reaching saturation'), notes_Aci, 'UniformOutput', 0));
    %%%V_cmax_Aci(is_unsat == 1) = nan;

    V_cmax23 = nan(N_data, 1);
    for i = 1:N_Trt

        Trt_name_alt = Trt_name{i, 1};

        % remove 'Ad' and/or 'Bd' from ''Trt_name_alt''to be comparable to A_ci data
        ind_Ad = strfind(Trt_name_alt, 'Ad');
        if ~isempty(ind_Ad)
            Trt_name_alt = [Trt_name_alt(1:ind_Ad-1), Trt_name_alt(ind_Ad+2:end)];
        end
        ind_Bd = strfind(Trt_name_alt, 'Bd');
        if ~isempty(ind_Bd)
            Trt_name_alt = [Trt_name_alt(1:ind_Bd-1), Trt_name_alt(ind_Bd+2:end)];
        end

        is_Trt_Aci = cell2mat(cellfun(@(x) strcmp(x, Trt_name_alt), Trt_name_Aci, 'UniformOutput', 0));
        V_cmax_Aci_local = V_cmax_Aci(is_Trt_Aci == 1);
        V_cmax_Aci_local = max(V_cmax_Aci_local(~isnan(V_cmax_Aci_local)));
        V_cmax_Aci_local(isempty(V_cmax_Aci_local)) = nan;

        is_Trt = cell2mat(cellfun(@(x) strcmp(x, Trt_name{i, 1}), Comb_Trt_name, 'UniformOutput', 0));
        V_cmax23(is_Trt == 1) = V_cmax_Aci_local;

    end

    % correct V_cmax for temperature, considering that V_cmac measurements were performed at 23C
    V_cmax25_func = @(psi_L_MD) 0;
    [outputs_photo_param_at_T_L_23] = Photosynthesis_parameters_temperature_response(23, 0, V_cmax25_func);
    V_cmax_per_V_cmax25_at_T_L_23 = outputs_photo_param_at_T_L_23.V_cmax_per_V_cmax25_vect;
    V_cmax25 = V_cmax23 ./ V_cmax_per_V_cmax25_at_T_L_23;
    [outputs_photo_param] = Photosynthesis_parameters_temperature_response(T_a, 0, V_cmax25_func);
    V_cmax_per_V_cmax25 = outputs_photo_param.V_cmax_per_V_cmax25_vect;
    V_cmax = V_cmax25 .* V_cmax_per_V_cmax25;

    % other photosynthetic traits
    Gamma_star = outputs_photo_param.Gamma_star_vect;
    K_c = outputs_photo_param.K_c_vect;
    K_o = outputs_photo_param.K_o_vect;
    K_m = outputs_photo_param.K_m_vect;

    % estimate dark respiration
    R_d_per_R_d25 = outputs_photo_param.R_d_per_R_d25_vect;
    R_d25_per_V_cmax25 = 0.015; %from Collatz et al. (1991)
    R_d25 = R_d25_per_V_cmax25 .* V_cmax25;
    R_d = R_d_per_R_d25 .* R_d25;

    % change units for photosynthetic variables
    V_cmax = V_cmax*1e-6; %convert from [mumol m-2 s-1] to [mol m-2 s-1]
    V_cmax25 = V_cmax25*1e-6; %convert from [mumol m-2 s-1] to [mol m-2 s-1]
    R_d = R_d*1e-6; %convert from [mumol m-2 s-1] to [mol m-2 s-1];

    % store V_cmax data
    Var = setfield(Var, 'V_cmax', V_cmax);
    Var = setfield(Var, 'V_cmax25', V_cmax25);
    Var = setfield(Var, 'Gamma_star', Gamma_star);
    Var = setfield(Var, 'K_c', K_c);
    Var = setfield(Var, 'K_o', K_o);
    Var = setfield(Var, 'K_m', K_m);
    Var = setfield(Var, 'R_d', R_d);
    save(temporary_save_name, 'Var')

    %% Convert units
    F_d = F_d/18/3600; %convert from [g h-1] to [mol s-1]
    RH = RH/100; %convert from percent to decimal

    %% Estimate stomatal conductance
    LA = nan(N_data, 1); %leaf area in [m2]
    is_Picea = cell2mat(cellfun(@(x) strcmp(x, 'Picea'), species, 'UniformOutput', 0));
    is_Larix = cell2mat(cellfun(@(x) strcmp(x, 'Larix'), species, 'UniformOutput', 0));
    LA(is_Picea == 1) = 0.457 * SWA(is_Picea == 1); %equation from Peters et al. (2018; Plant, Cell & Environment)
    LA(is_Larix == 1) = 0.532 * SWA(is_Larix == 1); %equation from Peters et al. (2018; Plant, Cell & Environment)

    E = F_d ./ LA; %transpiration per leaf area in [mol m-2 s-1]
    T_a_pot = T_a + gamma_T_a*Elevation; %potential temperature in [C]
    P_atm = P_atm_0 * exp(-g*Elevation/r_a./(T_a_pot+273.15)); %atmospheric pressure [kPa]
    VPD_L = (1-RH) * 0.61078 .* exp(17.27*T_a./(T_a + 237.3)); %vapor pressure deficit [kPa]
    VPD_L(VPD_L < 0.1) = nan;
    g_w = E .* P_atm ./ VPD_L; %stomatal conductance [mol m-2 s-1]

    % find maximum g_w for each treatment
    g_w_max_Trt = nan(N_data, 1);
    for i = 1:N_Trt
        is_Trt = cell2mat(cellfun(@(x) strcmp(x, Trt_name{i, 1}), Comb_Trt_name, 'UniformOutput', 0));
        g_w_max_Trt(is_Trt == 1) = max(g_w(is_Trt == 1));
    end

    Var = setfield(Var, 'E', E);
    Var = setfield(Var, 'P_atm', P_atm);
    Var = setfield(Var, 'VPD_L', VPD_L);
    Var = setfield(Var, 'g_w', g_w);
    Var = setfield(Var, 'g_w_max_Trt', g_w_max_Trt);
    save(temporary_save_name, 'Var')

    % solve for net carbon assimilation
    c_a = 400e-6 * ones(N_data, 1);
    c_i = 0.5*((c_a - K_m - 1.6./g_w.*(V_cmax - R_d)) + ((c_a - K_m - 1.6./g_w.*(V_cmax - R_d)).^2 + 4*1.6./g_w.*(V_cmax.*Gamma_star + R_d.*K_m) + 4*c_a.*K_m).^0.5);
    c_i(c_i < Gamma_star) = nan; % remove bad c_i
    c_i(c_i > c_a) = nan; % remove bad c_i
    A_n = V_cmax .* (c_i - Gamma_star)./(c_i + K_m) - R_d;
    k = V_cmax .* (Gamma_star + K_m) ./ (c_i + K_m).^2; % Estimate k = dA_n/dc_i
    lambda = 1.6*k./(g_w + 1.6*k).*A_n./E;
    lambda(lambda < 0) = nan;

    Var = setfield(Var, 'c_i', c_i);
    Var = setfield(Var, 'A_n', A_n);
    Var = setfield(Var, 'lambda', lambda);
    save(temporary_save_name, 'Var')

    % remove outliers from (log10-transformed) lambda
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

    % remove faulty & zero g_w estimates
    g_w(g_w_max_Trt < 0.1) = nan;
    g_w(g_w == 0) = nan;

    % only consider midday values (when A_n is carboxylation-limited)
    TimeofDay = timeofday(Date);
    g_w(TimeofDay > hours(15)) = nan;
    g_w(TimeofDay < hours(11)) = nan;

    %% Saturated water content
    SWC_L = nan(N_data, 1);
    SWC_L(is_Picea == 1) = 2.8; %based on Belluau et al. (2021; Functional Ecology)
    SWC_L(is_Larix == 1) = 3.0; %based on Belluau et al. (2021; Functional Ecology)
    
    
    %% Leaf area-specific plant soil hydraulic conductance
    Q_10_k_L = 1.5;
    E(isnan(g_w)) = nan;
    k_L = E ./ (psi_L_PD - psi_L_MD);
    k_L(k_L == inf) = nan;
    k_L(k_L <= 0) = nan;
    k_L((hod <= hour_day_start) + (hod >= hour_day_end) == 2) = nan;
    k_L_25C = k_L ./ (Q_10_k_L .^ ((T_a - 25)/10));

    species_unique = unique(species);
    N_species_unique = length(species_unique);
    k_L_25C_max_species_unique = nan(N_species_unique, 1);
    A_species_unique = nan(N_species_unique, 1);
    
    for i = 1:N_species_unique
        
        is_species_local = cell2mat(cellfun(@(x) strcmp(x, species_unique{i}), species, 'UniformOutput', 0));
        k_L_25C_local = k_L_25C(is_species_local == 1);
        psi_L_MD_local = psi_L_MD(is_species_local == 1);
        Trt_numb_local = Trt_numb(is_species_local == 1);
        Comb_Trt_name_local = Comb_Trt_name(is_species_local == 1);
        N_Trt_local = max(Trt_numb_local);
        
        % specify names of species and treatments for plots
        Comb_Trt_name_plot = cell(1, N_Trt_local);
        for j = 1:N_Trt_local
            Comb_Trt_name_plot(j) = unique(Comb_Trt_name_local(Trt_numb_local == j));
        end
        [species_plot, Trt_plot] = Specify_species_and_Trt_names(species_unique{i}, Comb_Trt_name_plot);
        
        % add tree-specific labels if necesary 
        [~, ind_order] = sort(Trt_plot); %ordered alphabetically by legend entry
        checked = zeros(1, N_Trt_local);
        for j = ind_order
            if ~checked(j)
                is_Trt_subset_combo_plot_same = cell2mat(cellfun(@(x) strcmp(Trt_plot{j}, x), Trt_plot, 'UniformOutput', 0));
                if sum(is_Trt_subset_combo_plot_same) > 1
                    ind_change = find(is_Trt_subset_combo_plot_same);
                    l = 0;
                    for k = ind_change
                        l = l + 1;
                        Trt_plot{k} = [Trt_plot{k}, ' - ', num2str(l)];
                    end
                end
                checked(is_Trt_subset_combo_plot_same == 1) = 1;
            end
        end
        
        figure
        title(species_plot)
        hold on
        for j = 1:N_Trt_local
            k_L_25C_local_Trt = k_L_25C_local(Trt_numb_local == j);
            if any(~isnan(k_L_25C_local_Trt))
                plot(-psi_L_MD_local(Trt_numb_local == j), k_L_25C_local_Trt, 'o', 'Color', colors(j,:), 'MarkerFaceColor', colors(j,:), 'DisplayName', Trt_plot{j})
            end
        end
        hold off
        set(gca, 'yscale', 'log')
        legend('location', 'eastoutside')
        xlabel('-{\it\psi_L} [MPa]')
        ylabel(['{\itk_L}_{25', char(176), '} [mol\cdotm^{-2}\cdots^{-1}\cdotMPa^{-1}]'])
        
        % for k_L_25C = k_L_25C_max * exp(A*psi_L)
        [p, S] = polyfit(psi_L_MD_local(~isnan(k_L_25C_local)), log(k_L_25C_local(~isnan(k_L_25C_local))), 1);
        A_local = p(1);
        k_L_25C_max_local = exp(p(2));
        
        % check significance
        mdl = fitlm(psi_L_MD_local(~isnan(k_L_25C_local)), log(k_L_25C_local(~isnan(k_L_25C_local))));
        p_val = coefTest(mdl);

        % plot exponential fit
        xlim_max = max(xlim);
        xlim([0, xlim_max]);
        x = 0:0.01:xlim_max;
        hold on
        plot(x, k_L_25C_max_local*exp(-A_local*x), 'k', 'LineWidth', 2, 'DisplayName', ['Exponential fit, p = ', num2str(p_val, 3)])
        hold off
        
        % 95 confidence interval
        [y_fit,delta] = polyval(p,-x,S);
        UB = exp(y_fit + 2*delta);
        LB = exp(y_fit - 2*delta);
        hold on
        fill([x, fliplr(x)], [UB, fliplr(LB)], 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.1, 'DisplayName', '95% Confidence interval')
        hold off
        
        if p_val > 0.05
            
            % plot average
            k_L_25C_average_local = mean(k_L_25C_local(~isnan(k_L_25C_local)));
            exponent = floor(log10(k_L_25C_average_local));
            significand = k_L_25C_average_local/10^exponent;
            hold on
            plot(x, k_L_25C_average_local*ones(size(x)), 'k--', 'LineWidth', 1.5, ...
                 'DisplayName', ['Mean {\itk_L}_{25', char(176), '}  = ', num2str(significand, 2), ' 	\times 10^{', num2str(exponent), '}', newline, 'mol\cdotmol^{-2}\cdots^{-1}\cdotMPa^{-1}'])
            hold off
            
            % plot most senesitive from mean as max
            A_max_local = -log(LB(end)/k_L_25C_average_local)/x(end);
            hold of
            plot(x, k_L_25C_average_local*exp(-A_max_local*x), 'k:', 'LineWidth', 1, ...
                 'DisplayName', ['{\itk_{L,max}}_{25', char(176), '} = mean {\itk_{L,max}}_{25', char(176), '},', newline, 'maximum possible vulnerability' newline, '({\itS_k} = ' num2str(A_max_local, 3), ' MPa^{-1})'])
            hold off
        end
        
        % save plot
        print(gcf, '-djpeg', '-r300', [species_unique{i}, '_k_L_25C_max.jpg'], '-painters', '-noui' )
        
        % store outputs
        k_L_25C_max_species_unique(i) = k_L_25C_max_local;
        A_species_unique(i) = A_local;
        
    end
    
    
    %% Save data as .mat
    save(Mat_data_save_file, 'T_a', 'A_n', 'c_i', 'c_a', 'g_w', 'E', ...
         'VPD_L', 'P_atm', 'Gamma_star', 'K_c', 'K_o', 'K_m', 'V_cmax', ...
         'R_d', 'lambda', 'psi_L_MD', 'psi_L_PD', 'species', 'Trt_numb', ...
         'Comb_Trt_name', 'SWC_L', 'k_L_25C', ...
         'species_unique', 'k_L_25C_max_species_unique', 'A_species_unique');

    
end


%% Store Outputs
if isempty(Var)
    Var.T_L = T_a;
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
    Var.T_L = [Var.T_L; T_a];
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

