function X = prior(N,Nx,range)

xmin = range(:,1);
xmax = range(:,2);
X = nan(N,Nx);

%% Option for choosing subset and best parameters
% Option = 1 -- uniform distribution as used in the original DREAM code
% Option = 2 -- modified so that the 'prior' function only generates
% parameters that result in a realistic response of stomata to relative water content (for an assumed set of values for the osmotic potential)
Option = 2;

switch Option
    case 1
        
        for i = 1:N
            X(i,:) = unifrnd(xmin,xmax);
        end 
        
    case 2
        
        % Case 2-specific parameters
        N_test = 21;
        pi_L_0_25C_test = linspace(-1, -3, N_test);
        
        for i = 1:N

            j = 0;
            while 1

                j = j + 1;
                if j > 1e3
                    error('ERROR: Cannot find parameters that allow for realistic psi_L-chi_w response!')
                end

                X_local = unifrnd(xmin,xmax);
                Check = zeros(1, N_test);
                for j = 1:N_test
                    pi_L_0_25C_local = pi_L_0_25C_test(j);
                    Check_local = Check_X_for_RWC_response(X_local, pi_L_0_25C_local); 
                    Check(j) = Check_local; 
                end
                if any(Check)
                    break 
                end

            end

            X(i,:) = X_local;

        end    

end
end