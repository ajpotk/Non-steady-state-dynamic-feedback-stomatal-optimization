function [ F1_collated, F2_collated, Date_collated ] = Collate_data( F1, Date1, F2, Date2, varargin )

% ''Collate_data()'' finds data points in ''F1' and ''F2'' that share a
% common Date (i.e., ''Date1'' == ''Date2'')

diff_Date_allow = seconds(0); %default value

if ~isempty(varargin)
    if isduration(varargin{1})
        diff_Date_allow = varargin{1};
    else
        error('ERROR: First ''varargin'' entry must be ''duration'' class!')
    end
end

if ~isdatetime(Date1)
   error('ERROR: ''Date1'' must be a ''duration''!') 
end
if ~isdatetime(Date2)
   error('ERROR: ''Date2'' must be a ''duration''!') 
end

[N1, ~] = size(F1);
[N2, ~] = size(F2);
swap = 0;

if N2 < N1
    
    swap = 1;
    %swap F1 and F2 so that shorter array is in the position of F1
    F1_orig = F1;
    F2_orig = F2;
    Date1_orig = Date1;
    Date2_orig = Date2;
    N1_orig = N1;
    N2_orig = N2;
    
    F1 = F2_orig;
    F2 = F1_orig;
    Date1 = Date2_orig;
    Date2 = Date1_orig;
    N1 = N2_orig;
    N2 = N1_orig;
    
end

% sort data by time
[Date1, ind_sort1] = sort(Date1);
F1 = F1(ind_sort1, :);
[Date2, ind_sort2] = sort(Date2);
F2 = F2(ind_sort2, :);

Date1 = reshape(Date1, N1, 1);
Date2 = reshape(Date2, N2, 1);
diff_Date = Date1' - Date2;
[min_abs_diff_Date, ind_min_abs_diff_Date] = min(abs(diff_Date));
is_allow = (min_abs_diff_Date <= diff_Date_allow);

F1_collated = F1((is_allow == 1), :);
F2_collated = F2(ind_min_abs_diff_Date(is_allow == 1), :); 
Date_collated = Date1((is_allow == 1), :);

if swap
    %unswap
    F1_collated_swap = F1_collated;
    F2_collated_swap = F2_collated;
    
    F1_collated = F2_collated_swap;
    F2_collated = F1_collated_swap;
end
    
end

