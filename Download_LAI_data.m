function [LAI] = Download_LAI_data(Date)

%% User-specified
Year_first = min(year(Date));
Year_last = max(year(Date));
LAT_select = 50.78666; %46.81533;
LONG_select = 13.72129; %9.85591;

save_file_extension = '_DE-Obe'; %'_CH-Dav';

Mat_data_save_file = 'Saved_LAI_data';
Mat_data_save_file = [Mat_data_save_file, save_file_extension, '.mat']; 

files = struct2cell(dir);
files = files(1,:);
is_Mat_data_save_file = cell2mat(cellfun(@(x) strcmp(x, Mat_data_save_file), files, 'UniformOutput', 0));

if any(is_Mat_data_save_file)
    
    load(Mat_data_save_file)
    
else

    %% Data location
    home = pwd;
    mother_directory_data = "C:\Users\potka002\Desktop\Research\Growth Optimizing Stomata\Test_GOH\Fluxnet\Sin_5km_bimonth";

    % Sinusoidal projection
    LAT_map_limited = linspace(90, -90, 4320)';
    LONG_map = linspace(-180, 180, 8640);
    LONG_map = LONG_map./cos(LAT_map_limited * pi/180);
    LONG_map(abs(LONG_map) > 180) = nan;
    LAT_map = repmat(LAT_map_limited, 1, 8640);

    error_select = ((LAT_map - LAT_select).^2 + (LONG_map - LONG_select).^2).^0.5;
    min_error = min(error_select, [], 'all');
    ind_selct = find(error_select == min_error, 1, 'first');
    LAT_found = LAT_map(ind_selct);
    LONG_found = LONG_map(ind_selct);

    ind_LAT_found = find(LAT_map_limited == LAT_found);
    LONG_map_limited = LONG_map(ind_LAT_found,:);
    ind_LONG_found = find(LONG_map_limited == LONG_found);

    N_bimonthly = 24*(Year_last + 1 - Year_first);
    LAI_bimonthly = nan(1, N_bimonthly);
    Date_bimonthly = NaT(1, N_bimonthly);
    ind_bimonthly = 0;
    cd(mother_directory_data)
    for Year = Year_first:Year_last

        Year_str = num2str(Year);

        for Month = 1:12

            Month_str = num2str(Month);
            len_Month_str = length(Month_str);
            while len_Month_str < 2
                Month_str = ['0', Month_str];
                len_Month_str = length(Month_str);
            end

            Days_in_Month = eomday(Year, Month);

            for Half = 1:2

                Half_str = ['0', num2str(Half)];
                file_name = ['SI_LAI_FPAR_CDR_Sin_5km_bimonth_', Year_str, Month_str, Half_str, '.tif'];
                disp(['Downloading LAI data from: ', file_name])
                TIF = Tiff(file_name,'r');
                TIF_data = double(read(TIF));
                TIF_data(TIF_data == 255) = nan;
                LAI_map = TIF_data(:,:,1)/10;

                ind_bimonthly = ind_bimonthly + 1;
                LAI_bimonthly(ind_bimonthly) = LAI_map(ind_LAT_found, ind_LONG_found);
                Day = (0.5*Half-0.25)*Days_in_Month;
                Hour = 24*(Day - floor(Day));
                Day = floor(Day);
                Date_bimonthly(ind_bimonthly) = datetime(Year, Month, Day, Hour, 0, 0);

            end
        end
    end
    cd(home)
    
    %% Save data as .mat
    save(Mat_data_save_file, 'Date_bimonthly', 'LAI_bimonthly');

end

%% Interpolate LAI for selected Dates
LAI = interp1(Date_bimonthly, LAI_bimonthly, Date, 'linear', 'extrap');
LAI(LAI < 0) = 0;

end

