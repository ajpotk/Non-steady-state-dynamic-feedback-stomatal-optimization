function [Date, P, T_a, RH, P_atm, R_0, g_H_a, c_a, GPP, ET, SPEI, dt] = Download_Flux_data()

%% Data location
home = pwd;

mother_directory_data = "C:\Users\potka002\Desktop\Research\Growth Optimizing Stomata\Test_GOH\Fluxnet";
file_name = "FLX_DE-Obe_FLUXNET2015_FULLSET_HH_2008-2020_beta-3.csv"; %"FLX_CH-Dav_FLUXNET2015_FULLSET_HH_1997-2020_beta-3.csv";

save_file_extension = '_DE-Obe'; %'_CH-Dav';

Mat_data_save_file = 'Saved_Flux_data';
Mat_data_save_file = [Mat_data_save_file, save_file_extension, '.mat']; 

files = struct2cell(dir);
files = files(1,:);
is_Mat_data_save_file = cell2mat(cellfun(@(x) strcmp(x, Mat_data_save_file), files, 'UniformOutput', 0));

if any(is_Mat_data_save_file)
    
    load(Mat_data_save_file)
    
else

    %% Local constants
    R = 8.314; %universal gas constant [J mol-1 K-1]
    m_w = 18e-3; %molar mass of water [kg mol-1]
    d_L = 1e-3; %characteristic leaf dimension [m] -- based on Niinements & Kull (1995) for Picea Abies
    C_v = 0.01; %turbulent transfer coefficient [m s-1/2]

    %% User-specified
    Date_begin = datetime(2017, 1, 1);
    Date_end = datetime(2020, 1, 1); 

    %% Set up import options and import the data
    opts = delimitedTextImportOptions("NumVariables", 247);

    % Specify range and delimiter
    opts.DataLines = [2, Inf];
    opts.Delimiter = ",";

    % Specify column names and types
    opts.VariableNames = ["TIMESTAMP_START", "TIMESTAMP_END", "TA_F_MDS", "TA_F_MDS_QC", "TA_ERA", "TA_F", "TA_F_QC", "SW_IN_POT", "SW_IN_F_MDS", "SW_IN_F_MDS_QC", "SW_IN_ERA", "SW_IN_F", "SW_IN_F_QC", "LW_IN_F_MDS", "LW_IN_F_MDS_QC", "LW_IN_ERA", "LW_IN_F", "LW_IN_F_QC", "LW_IN_JSB", "LW_IN_JSB_QC", "LW_IN_JSB_ERA", "LW_IN_JSB_F", "LW_IN_JSB_F_QC", "VPD_F_MDS", "VPD_F_MDS_QC", "VPD_ERA", "VPD_F", "VPD_F_QC", "PA", "PA_ERA", "PA_F", "PA_F_QC", "P", "P_ERA", "P_F", "P_F_QC", "WS", "WS_ERA", "WS_F", "WS_F_QC", "WD", "USTAR", "RH", "NETRAD", "PPFD_IN", "PPFD_DIF", "PPFD_OUT", "SW_DIF", "SW_OUT", "LW_OUT", "CO2_F_MDS", "CO2_F_MDS_QC", "TS_F_MDS_1", "TS_F_MDS_2", "TS_F_MDS_3", "TS_F_MDS_4", "TS_F_MDS_5", "TS_F_MDS_6", "TS_F_MDS_1_QC", "TS_F_MDS_2_QC", "TS_F_MDS_3_QC", "TS_F_MDS_4_QC", "TS_F_MDS_5_QC", "TS_F_MDS_6_QC", "SWC_F_MDS_1", "SWC_F_MDS_2", "SWC_F_MDS_3", "SWC_F_MDS_4", "SWC_F_MDS_5", "SWC_F_MDS_6", "SWC_F_MDS_1_QC", "SWC_F_MDS_2_QC", "SWC_F_MDS_3_QC", "SWC_F_MDS_4_QC", "SWC_F_MDS_5_QC", "SWC_F_MDS_6_QC", "G_F_MDS", "G_F_MDS_QC", "LE_F_MDS", "LE_F_MDS_QC", "LE_CORR", "LE_CORR_25", "LE_CORR_75", "LE_RANDUNC", "LE_RANDUNC_METHOD", "LE_RANDUNC_N", "LE_CORR_JOINTUNC", "H_F_MDS", "H_F_MDS_QC", "H_CORR", "H_CORR_25", "H_CORR_75", "H_RANDUNC", "H_RANDUNC_METHOD", "H_RANDUNC_N", "H_CORR_JOINTUNC", "EBC_CF_N", "EBC_CF_METHOD", "NIGHT", "NEE_CUT_REF", "NEE_VUT_REF", "NEE_CUT_REF_QC", "NEE_VUT_REF_QC", "NEE_CUT_REF_RANDUNC", "NEE_VUT_REF_RANDUNC", "NEE_CUT_REF_RANDUNC_METHOD", "NEE_VUT_REF_RANDUNC_METHOD", "NEE_CUT_REF_RANDUNC_N", "NEE_VUT_REF_RANDUNC_N", "NEE_CUT_REF_JOINTUNC", "NEE_VUT_REF_JOINTUNC", "NEE_CUT_USTAR50", "NEE_VUT_USTAR50", "NEE_CUT_USTAR50_QC", "NEE_VUT_USTAR50_QC", "NEE_CUT_USTAR50_RANDUNC", "NEE_VUT_USTAR50_RANDUNC", "NEE_CUT_USTAR50_RANDUNC_METHOD", "NEE_VUT_USTAR50_RANDUNC_METHOD", "NEE_CUT_USTAR50_RANDUNC_N", "NEE_VUT_USTAR50_RANDUNC_N", "NEE_CUT_USTAR50_JOINTUNC", "NEE_VUT_USTAR50_JOINTUNC", "NEE_CUT_MEAN", "NEE_VUT_MEAN", "NEE_CUT_MEAN_QC", "NEE_VUT_MEAN_QC", "NEE_CUT_SE", "NEE_VUT_SE", "NEE_CUT_05", "NEE_CUT_16", "NEE_CUT_25", "NEE_CUT_50", "NEE_CUT_75", "NEE_CUT_84", "NEE_CUT_95", "NEE_VUT_05", "NEE_VUT_16", "NEE_VUT_25", "NEE_VUT_50", "NEE_VUT_75", "NEE_VUT_84", "NEE_VUT_95", "NEE_CUT_05_QC", "NEE_CUT_16_QC", "NEE_CUT_25_QC", "NEE_CUT_50_QC", "NEE_CUT_75_QC", "NEE_CUT_84_QC", "NEE_CUT_95_QC", "NEE_VUT_05_QC", "NEE_VUT_16_QC", "NEE_VUT_25_QC", "NEE_VUT_50_QC", "NEE_VUT_75_QC", "NEE_VUT_84_QC", "NEE_VUT_95_QC", "RECO_NT_VUT_REF", "RECO_NT_VUT_USTAR50", "RECO_NT_VUT_MEAN", "RECO_NT_VUT_SE", "RECO_NT_VUT_05", "RECO_NT_VUT_16", "RECO_NT_VUT_25", "RECO_NT_VUT_50", "RECO_NT_VUT_75", "RECO_NT_VUT_84", "RECO_NT_VUT_95", "RECO_NT_CUT_REF", "RECO_NT_CUT_USTAR50", "RECO_NT_CUT_MEAN", "RECO_NT_CUT_SE", "RECO_NT_CUT_05", "RECO_NT_CUT_16", "RECO_NT_CUT_25", "RECO_NT_CUT_50", "RECO_NT_CUT_75", "RECO_NT_CUT_84", "RECO_NT_CUT_95", "GPP_NT_VUT_REF", "GPP_NT_VUT_USTAR50", "GPP_NT_VUT_MEAN", "GPP_NT_VUT_SE", "GPP_NT_VUT_05", "GPP_NT_VUT_16", "GPP_NT_VUT_25", "GPP_NT_VUT_50", "GPP_NT_VUT_75", "GPP_NT_VUT_84", "GPP_NT_VUT_95", "GPP_NT_CUT_REF", "GPP_NT_CUT_USTAR50", "GPP_NT_CUT_MEAN", "GPP_NT_CUT_SE", "GPP_NT_CUT_05", "GPP_NT_CUT_16", "GPP_NT_CUT_25", "GPP_NT_CUT_50", "GPP_NT_CUT_75", "GPP_NT_CUT_84", "GPP_NT_CUT_95", "RECO_DT_VUT_REF", "RECO_DT_VUT_USTAR50", "RECO_DT_VUT_MEAN", "RECO_DT_VUT_SE", "RECO_DT_VUT_05", "RECO_DT_VUT_16", "RECO_DT_VUT_25", "RECO_DT_VUT_50", "RECO_DT_VUT_75", "RECO_DT_VUT_84", "RECO_DT_VUT_95", "RECO_DT_CUT_REF", "RECO_DT_CUT_USTAR50", "RECO_DT_CUT_MEAN", "RECO_DT_CUT_SE", "RECO_DT_CUT_05", "RECO_DT_CUT_16", "RECO_DT_CUT_25", "RECO_DT_CUT_50", "RECO_DT_CUT_75", "RECO_DT_CUT_84", "RECO_DT_CUT_95", "GPP_DT_VUT_REF", "GPP_DT_VUT_USTAR50", "GPP_DT_VUT_MEAN", "GPP_DT_VUT_SE", "GPP_DT_VUT_05", "GPP_DT_VUT_16", "GPP_DT_VUT_25", "GPP_DT_VUT_50", "GPP_DT_VUT_75", "GPP_DT_VUT_84", "GPP_DT_VUT_95", "GPP_DT_CUT_REF", "GPP_DT_CUT_USTAR50", "GPP_DT_CUT_MEAN", "GPP_DT_CUT_SE", "GPP_DT_CUT_05", "GPP_DT_CUT_16", "GPP_DT_CUT_25", "GPP_DT_CUT_50", "GPP_DT_CUT_75", "GPP_DT_CUT_84", "GPP_DT_CUT_95", "RECO_SR", "RECO_SR_N"];
    opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

    % Specify file level properties
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";

    % Import the data
    cd(mother_directory_data)
    disp('Downloading Flux data')
    Flux_data = readtable(file_name, opts);
    clear opts
    cd(home)

    %% Choose data from selected years
    var_names = Flux_data.Properties.VariableNames;
    data_select = {'TIMESTAMP_START',   'Date',     1; ...
                   'TA_ERA',            'T_a',      1; ...
                   'SW_IN_ERA',         'R_0',      1; ...
                   'VPD_ERA',           'VPD',      1e-1; ... %[hPa] to [kPa]
                   'PA_ERA',            'P_atm',    1; ...
                   'P_ERA',             'P',        1e-3; ... %[mm] to [m]
                   'WS_ERA',            'u',        1; ...
                   'USTAR',             'u_star',   1;
                   'CO2_F_MDS',         'c_a',      1e-6; ... %[umol mol-1] to [mol mol-1]
                   'GPP_NT_VUT_REF',    'GPP',      1e-6; ... %[umol m-2 s-1] to [mol m-2 s-1]
                   'LE_CORR',           'LE',       1};
    [N_data_select, ~] = size(data_select);
    Var = struct;
    for i = 1:N_data_select

        ind_var_local = find(cell2mat(cellfun(@(x) strcmp(x, data_select{i,1}), var_names, 'UniformOutput', 0)), 1, 'first');
        data_local = table2array(Flux_data(:,ind_var_local));
        data_local(data_local == -9999) = nan;
        data_local = data_select{i,3} * data_local; %convert units
        if any(isnan(data_local))
            data_local = fillmissing(data_local, 'linear');
        end
        Var = setfield(Var, data_select{i,2}, data_local);

    end

    Date = Var.Date;

    Date = cellfun(@(x) num2str(x), num2cell(Date), 'UniformOutput', 0);
    Year = cell2mat(cellfun(@(x) str2num(x(1:4)), Date, 'UniformOutput', 0));
    Month = cell2mat(cellfun(@(x) str2num(x(5:6)), Date, 'UniformOutput', 0));
    Day = cell2mat(cellfun(@(x) str2num(x(7:8)), Date, 'UniformOutput', 0));
    Hour = cell2mat(cellfun(@(x) str2num(x(9:10)), Date, 'UniformOutput', 0));
    Minute = cell2mat(cellfun(@(x) str2num(x(11:12)), Date, 'UniformOutput', 0));
    Second = zeros(size(Date));
    Date = datetime(Year, Month, Day, Hour, Minute, Second);
    
    % conditions at all times
    T_a = Var.T_a;
    R_0 = Var.R_0;
    VPD = Var.VPD;
    P_atm = Var.P_atm;
    P = Var.P;
    u_star = Var.u_star;
    c_a = Var.c_a;
    GPP = Var.GPP;
    LE = Var.LE;
    
    % calculate SPEI
    SPEI = calculate_SPEI(Date, P, T_a, R_0);
    
    % conditions at desired times
    is_Date_select = ((Date >= Date_begin) + (Date < Date_end) == 2);
    Date = Date(is_Date_select == 1);
    T_a = T_a(is_Date_select == 1);
    R_0 = R_0(is_Date_select == 1);
    VPD = VPD(is_Date_select == 1);
    P_atm = P_atm(is_Date_select == 1);
    P = P(is_Date_select == 1);
    u_star = u_star(is_Date_select == 1);
    c_a = c_a(is_Date_select == 1);
    GPP = GPP(is_Date_select == 1);
    LE = LE(is_Date_select == 1);
    SPEI = SPEI(is_Date_select == 1);
    
    
    c_a = c_a + 400e-6; %ad hoc fix, because CO2 concentrations are too small
    
    e_a_sat = 0.61078 .* exp(17.27 * T_a ./ (T_a + 237.3)); %air vapor pressure in [kPa] -- Teten's equation
    RH = 1 - VPD./e_a_sat;
    RH(RH < 0) = 0;
    RH(RH > 1) = 1;

    g_H_a = (1e3 * P_atm)/R./(T_a + 273.5)*C_v.*(u_star/d_L).^0.5; %boundary layer heat conductance [mol m-2 s-1]
    
    lambda_E = 2.501e6 - 2.361e3*T_a; %latent heat of vaporization [J kg-1] -- expression from Allen et al. (1998; Crop evapotranspiration - Guidelines for computing crop water requirements. FAO Irrigation and drainage Paper 56. http://www.fao.org/3/x0490e/x0490e00.htm)
    ET = LE./lambda_E/m_w; %evapotranspiration [mol m-2 s-1]
    
    dt = seconds(mean(diff(Date))); %[s]
    
    %% Save data as .mat
    save(Mat_data_save_file, 'Date', 'P', 'T_a', 'RH', 'P_atm', 'R_0', 'g_H_a', 'c_a', 'GPP', 'ET', 'SPEI', 'dt');

end

end

