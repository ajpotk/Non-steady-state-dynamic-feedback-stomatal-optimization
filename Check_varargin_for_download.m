function [Var] = Check_varargin_for_download(varargin)

Var = []; 

if ~isempty(varargin)
    varargin = varargin{1}; %because this varagin is another code's varargin

    if length(varargin) > 1
        %check that only one input was given
        error('ERROR: Only ''Var'' is accepted as a single optional input!')
    elseif length(varargin) == 1

        Var = varargin{1};

        % check that single input is a structure
        if ~isstruct(Var)
            error('ERROR: Input ''Var'' must be a structure!')
        end

        % check that single input contains any of the relevant fields
        varargin_field = fieldnames(Var);
        Var_field_opt = {'T_L'; ... %list of potetnial fields in ''varargin''
                         'A_n'; ...
                         'c_i'; ...
                         'c_a';  ...
                         'g_w'; ...
                         'E'; ...
                         'VPD_L'; ... 
                         'P_atm';  ...
                         'Gamma_star'; ...
                         'K_c'; ...
                         'K_o'; ...
                         'K_m'; ...
                         'V_cmax'; ...
                         'R_d'; ...
                         'lambda'; ...
                         'psi_L_MD'; ...
                         'psi_L_PD'; ...
                         'species'; ...
                         'Trt_numb'; ...
                         'Comb_Trt_name'; ...
                         'SWC_L'};

        N_Var_field_opt = length(Var_field_opt);
        Contains_Var_field_opt = zeros(N_Var_field_opt,1);
        for i = 1:N_Var_field_opt
            is_Var_field_opt_i = cell2mat(cellfun(@(x) strcmp(x, Var_field_opt{i}), varargin_field, 'UniformOutput', 0));
            if any(is_Var_field_opt_i)
                Contains_Var_field_opt(i) = 1; 
            end
        end

        if any(Contains_Var_field_opt)
            % if varargin contained any relevant fields, then first check length of each field
            length_field_opt = nan(N_Var_field_opt,1);
            for i = 1:N_Var_field_opt
                if Contains_Var_field_opt(i) == 1
                    length_field_opt(i) = length(getfield(Var, Var_field_opt{i}));
                end
            end
            length_field_opt = length_field_opt(~isnan(length_field_opt)); 
            if all(length_field_opt == length_field_opt(1))
                length_field = length_field_opt(1);
            else
                error('ERROR: Fields in ''Var'' structure must be the same lengths!')
            end

            % then second fill missing fields with 'none' or nan's
            for i = 1:N_Var_field_opt
                if Contains_Var_field_opt(i) == 0
                    if strcmp(Var_field_opt{i}, 'species') || strcmp(Var_field_opt{i}, 'Comb_Trt_name')
                        % if missing cell of text, then set to 'none'
                        Var = setfield(Var, Var_field_opt{i}, repmat({'none'},length_field,1));
                    else
                        % if missing numeric, then set to nan
                        Var = setfield(Var, Var_field_opt{i}, nan(length_field,1));
                    end
                end
            end

        else
            % if varargin did not contain any relevant fields, then error
            error(['ERROR: Input ''Var'' structure must contain at least one relevant field!', ...
                   newline, '       Relevant fields are:', newline, ...
                   cell2mat(cellfun(@(x) ['       ''', x, '''', newline], Var_field_opt, 'UniformOutput', 0)')])
        end

    end
end

end

