%% Stats

% Reshape it all to column vectors and concatenate
allData = cellfun(@(x)x(:), layerGroupedArr, 'UniformOutput', false);
allData = vertcat(allData{:});

% Generate group indices for each set of data as column vectors and
% concatenate
groups = arrayfun(@(x, y)y * ones(numel(x{:}), 1), layerGroupedArr', 1:length(layerGroupedArr), 'UniformOutput', false);
groups = vertcat(groups{:});

[p,tbl,stats] = kruskalwallis(allData,groups); % violates assumption of normality so use kw instead of anova
[c,m,h,gnames] = multcompare(stats,"Alpha",0.05, "ctype","dunn-sidak");

pvalues = c(:,end);
groupnames = num2cell(c(:,1:2),2);
