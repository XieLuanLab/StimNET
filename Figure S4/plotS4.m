%% Supplementary figure S4

% due to stimulation misalignment with trigger signals, had to stop
% mid-frequency sweep and start another with remaining frequencies which
% explains multiple files for each date
% still 3 sessions though (1 session = 15 thresholds for 15 frequencies) 

d = dir('data');
d = d(~ismember({d.name},{'.','..'}));
s = cell(numel(d),2);
x = [];
y = [];

for i = 1:numel(d)
    file_i = dir(fullfile('data',d(i).name,'*.mat'));
    load(fullfile(file_i.folder,file_i.name))
    results = staircase.results_arr;
    idx = ~cellfun(@isempty, results);
    s{i,1} = staircase.frequencies(idx);
    s{i,2} = staircase.thresholds;
    x = [x staircase.frequencies(idx)];
    y = [y staircase.thresholds];
end

% Plot
xx = unique(x);
arrx = cell(numel(xx),1);
arry = cell(numel(xx),1);
for i = 1:numel(xx)
    idx = find(x == xx(i));
    arry{i} = y(idx);
    arrx{i} = repmat(xx(i),1,numel(idx));
end

mean_thresholds = mean(cell2mat(arry),2);
x0 = xx;
y0 = mean_thresholds';
[fitresult, gof, curve_x, curve_y] = createExpFit(x0, y0);

x2 = 0:0.1:105;
y2 = fitresult.a * exp(-fitresult.b * x2) + fitresult.c;
plot(x2,y2,'color',[1, 0, 0, 0.5],'LineWidth',2);
hold on;

cmap = parula(numel(xx)+10);
for i = 1:numel(xx)
    yvals = [arry{i}]; xvals = [arrx{i}];
    [GC,GR] = groupcounts(yvals');
    idx = find(GC > 1);
    % jitter to prevent overlap 
    if ~isempty(idx)
        idx2 = find(yvals == GR(idx));
        xvals(idx2(1)) = xvals(idx2(1)) - 0.25;
        xvals(idx2(2)) = xvals(idx2(2)) + 0.25;
    end
    scatter(xvals,yvals,18,'r','filled');
    hold on;
end
errorbar(x0,y0,y0'-min(cell2mat(arry),[],2),max(cell2mat(arry),[],2)-y0','x','Color','k')
xlim([0 105])
xlabel('Frequency (Hz)');
ylabel('Detection threshold ({\mu}A)')
box on
set(gca,'fontsize', 15)
set(gca,'linewidth',1.5)

% Red shading over 6 Hz
hold on;
patchF = [1 2 3 4];
patchV = [5.5 1; 5.5 4; 6.5 4; 6.5 1];
S2 = patch('Faces',patchF,'Vertices',patchV,'FaceColor','red','FaceAlpha',.2);
S2.EdgeColor = 'none';