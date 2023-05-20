%% Supplementary figure S5
% Run this after plot4EFG 

dTh = 0.03; % uA/day
initColor = 'black';
chronicColor = '#77CCFF';

x0 = s_new(3).allData.daysArr; % session dates converted to indices

data = cell2mat(s_new(3).allData.Ldata);
if useChargeUnit == 1
    data = data * 1e-6 * 167e-6 / 1e-9;
end

for i = 1:size(data,1)
    y_scatter = data(i,:);
    scatter(x0(1:end-1),y_scatter(1:end-1),12,scatterColor,...
    'filled','MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.3,'jitter','on','jitterAmount',0.5);
    hold on
    scatter(x0,y_scatter,12,scatterColor,...
    'filled','MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.3,'jitter','on','jitterAmount',0.5);  
    hold on
end

y0 = mean(data,'omitnan');

xpre0 = x0(1:39);
ypre0 = y0(1:39);
xpost0 = x0(40:end);
ypost0 = y0(40:end);
[~, gof1, curve_xpre, curve_ypre] = createExpFit(xpre0, ypre0);
[~, gof2, curve_xpost, curve_ypost] = createExpFit(xpost0, ypost0);

curve_x = [curve_xpre' curve_xpost'];
curve_y = [curve_ypre' curve_ypost'];

% Find indices of last N sessions pre-shock and last N sessions post-shock 
N = 4;
% Find boundary
boundaryidx = 40;

% Color data points for comparison
y_pre_arr = cell(N,1);
y_post_arr = cell(N,1);
hold on;
for i = 1:N
    y_pre = data(:,40-i);
    y_pre_arr{i} = y_pre;
    y_post = data(:,end-i+1);
    y_post_arr{i} = y_post;
end

yo = ylim;
yl = ylim + 1;
patchV = [0 0; 0 yl(2); x0(40) yl(2); x0(40) 0]; patchF = [1 2 3 4];
S1 = patch('Faces',patchF,'Vertices',patchV,'FaceColor',initColor,'FaceAlpha',.1);
S1.EdgeColor = 'none';
hold on;
patchV = [x0(40) 0; x0(40) yl(2); x0(end)+1 yl(2); x0(end)+1 0];
S3 = patch('Faces',patchF,'Vertices',patchV,'FaceColor',chronicColor,'FaceAlpha',.1);
S3.EdgeColor = 'none';

xlim([0 x0(end)+1])
ylim(yo);

xlabel('Days of behavioral testing')

% Add missing points to x and y pre-shock segments
[ax, ~] = plot2axes([curve_xpre' 99 100], [curve_ypre' 4.0833 4.0833], scatterColor, 'yscale', 0.167, 'yloc', 'left');
hold on;
plot(curve_xpost',curve_ypost','color',scatterColor,'LineWidth',1);
% plot(curve_xpost',curve_ypost','color',scatterColor,'LineWidth',1);
ylabel(ax(1), 'Threshold ({\mu}A)');
ylabel(ax(2), 'Threshold (nC/phase)');
xline(100,'LineWidth',1)

%% Plot groups
group1 = cell2mat(y_pre_arr); % n=12 contacts x 4 sessions
group2 = cell2mat(y_post_arr); % n=9 contacts x 4 sessions
group2 = group2(~isnan(group2)); 

figure;
scatter(ones(1,numel(group1)),group1')
hold on;
scatter(2*ones(1,numel(group2)),group2')
xlim([0 3])
%% Krusal-Wallis test to find difference between the two groups 
groups = [group1;group2];
groupid = [ones(1,numel(group1))' ; 2*ones(1,numel(group2))'];
[p,h,stats] = kruskalwallis(groups,groupid); 