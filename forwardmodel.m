function  [Y] = forwardmodel( x )
    
% Y -- predicted values for parameters, x
[output_subset] = Stomata_and_PV_best_for_subset(x);
Y = output_subset.g_w; 


% save results in other .txt file
Output_file = 'Output_text.txt';
files = struct2cell(dir);
files = files(1,:);
is_output_file = cell2mat(cellfun(@(x) strcmp(x, Output_file), files, 'UniformOutput', 0));

% open or create 'Output_file'
fileID = fopen(Output_file, 'a'); 

% Output header names
Output_names = {'a_f_max', ...
                'pi_L_star', ...
                'beta', ...
                'epsilon_L_max', ...
                'SWC_L', ...
                'alpha', ...
                'SSE_best_total_subset', ...
                'N_Trt_combo', ...
                'SSE_best_Trt_combo', ...
                'pi_L_0_25C_Trt_combo', ...
                'pi_L_min_combo', ...
                'pi_L_max_combo'};

if ~any(is_output_file)
    % if creating 'Output_file' for the first time, then create header
    header_text = cellfun(@(x) [x, ', '], Output_names, 'UniformOutput', 0);
    header_text = horzcat(header_text{:});
    fprintf(fileID, header_text);
end

N_Outputs = length(Output_names);
Outputs_new = cell(1, N_Outputs);
for i = 1:6
    Outputs_new{i} = x(i);
    Outputs_new{i} = [num2str(Outputs_new{i}), ', '];
end
for i = 7:N_Outputs
    Outputs_new_local = getfield(output_subset, Output_names{i});
    Outputs_new_local = num2cell(Outputs_new_local);
    Outputs_new_local = cellfun(@(x) [num2str(x), ', '], Outputs_new_local, 'UniformOutput', 0);
    Outputs_new_local = horzcat(Outputs_new_local{:});
    Outputs_new{i} = Outputs_new_local;
end
Outputs_new = horzcat(Outputs_new{:});

fprintf(fileID, ['\r\n', Outputs_new]);

fclose(fileID);

end



