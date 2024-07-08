function [grouping] = Group_by_ANOVA(data)

% sort data by mean
[~, ind_sort] = sort(mean(data));
data = data(:,ind_sort);

% create group names
[~, N_group_name] = size(data);
group_name = num2cell(1:N_group_name);
group_name = cellfun(@(x) num2str(x), group_name, 'UniformOutput', 0);
group_name_cat = categorical(group_name);
group_name_cat = reordercats(group_name_cat, group_name);

% perform ANOVA and build 'same' matrix -- if the two means of ith and jth 
% distributions are not siginificantlydifferent, then index (i,j) is 1,
% otherwise, the index is zero -- only works if data has already been
% sorted by their means
[~, ~, stats] = anova1(data, group_name_cat, 'off');
results = multcompare(stats, 'Display', 'off');
[N_results, ~] = size(results);
same = eye(N_group_name);
for i = 1:N_results
    if results(i,6) > 0.05
        same(results(i,1), results(i,2)) = 1;
        same(results(i,2), results(i,1)) = 1;
    end
end

% determine groups -- square formations of ones in the 'same' matrix are a
% grouping
ind_start = 1:N_group_name;
ind_end = nan(1, N_group_name);
for i = ind_start
   ind_end_local = N_group_name;
   while 1
       if all(same(i:ind_end_local, i:ind_end_local), 'all')
          break 
       end
       ind_end_local = ind_end_local - 1;
   end
   ind_end(i) = ind_end_local;
end

% remove non-unique groups
[~, ind_uniq] = unique(ind_end);
ind_start = ind_start(ind_uniq);
ind_end = ind_end(ind_uniq);

% assign grouping numbers to each mean
N_groups = length(ind_start); %number of groupings
grouping_number = cell(1, N_group_name); %numerical identification of groupings
for i = 1:N_groups
    for j = ind_start(i):ind_end(i)
        grouping_number{j} = [grouping_number{j}, i];
    end
end

% assign grouping letters to each mean
grouping = cell(1, N_group_name); %alphabetical identification of groupings
grouping_letter_opts = num2cell(96 + (1:N_groups));
grouping_letter_opts = cellfun(@(x) char(x), grouping_letter_opts, 'UniformOutput', 0); %option of letters to identify each grouping
for i = 1:N_group_name
    group_letters_local = grouping_letter_opts(grouping_number{i});
    group_letters_local = horzcat(group_letters_local{:});
    grouping{i} = group_letters_local;
end

% unsort data
[~, ind_unsort] = sort(ind_sort);
grouping = grouping(ind_unsort);

end

