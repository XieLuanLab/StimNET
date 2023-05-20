load('data\data.mat');
addpath('util\');

nANIMALS = length(s);
good_chs = cell(nANIMALS,1);
good_idxs = cell(nANIMALS,1);
all_date_arr = cell(nANIMALS,1);

for i = 1:nANIMALS
    animalID = s(i).animalID;
    data = s(i).data;
    thresholdArr = data.thresholdArr; % behavioral thresholds for a session 
    testedChArr = data.testedChArr; % channels tested in that session
    dateArr = data.dateArr; % session dates
    daysArr = data.daysArr; % session dates converted to days
    intanCh = data.intanCh; % Intan channel IDs of stim channels 
    all_date_arr{i} = dateArr;
end

layerAssignment()
groupThresholdsByLayer();

%% Plot all channel thresholds for all sessions
figure; 
plotThresholdsWFit(); % all sessions, plot 
splitInitAndChronicPhases();
sgtitle('All channel thresholds for all sessions')
%% Statistics for barplot
calculateStats();

%% Figure 4E 
% Barplot 
figure;
Lidx = 1:5;
nLAYERS = numel(Lidx);
mean_thresholds = zeros(nLAYERS,1);
med_thresholds = zeros(nLAYERS,1);
CI_thresholds = zeros(nLAYERS,2);
for L = 1:nLAYERS
    Ldata = layerGroupedArr{L};
    Ldata = Ldata(~isnan(Ldata)); % current 
    Ldata = Ldata * 1e-6 * 167e-6 / 1e-9; % nC/phase 
    [mu,CI] = calculateCI(Ldata,95);
    med_thresholds(L) = median(Ldata);
    mean_thresholds(L) = mu;
    CI_thresholds(L,:) = CI; 
end

h1 = bar(1:nLAYERS,mean_thresholds,'LineWidth',1); 
h1.FaceColor = 'flat';
h1.EdgeColor = 'white';
h1.CData = default_c(1:5,:);
h1.CData = default_c(1,:);
hold on;
set(gca,'fontsize', 12)
set(gca,'fontweight','bold')
set(gca,'linewidth',2)
sig_idx = pvalues < 0.05;
sigstar(groupnames(sig_idx),pvalues(sig_idx))

er = errorbar(1:nLAYERS,mean_thresholds,mean_thresholds-CI_thresholds(:,1),CI_thresholds(:,2)-mean_thresholds,'vertical');  

er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.LineWidth = 1.5;

xtickangle(90)
set(gca,'YAxisLocation','right')
ytickangle(90)
h=gca; h.XAxis.TickLength = [0 0]; h.YAxis.TickLength = [0 0];

xticklabels({'L1','L2/3','L4','L5','L6'});
ylabel('Detection threshold ({\mu}A)')
ylabel('Detection threshold (nC/phase)')

%% Figure 4F
% Violin plot for 5 animals (all session thresholds, L4-6 channels)

addpath('violin')
figure;
xarr = [];
yarr = [];
newIDs = {'Mouse1','Mouse2','Mouse3','Mouse4','Mouse5'};
for i = 1:5
    animalID = s_new(i).animalID;
    L46data = s_new(i).allData.Ldata(3:5);
    temp = cell2mat(L46data);
    L46data = temp(:);
    L46data = L46data(~isnan(L46data));
    y_data = L46data * 1e-6 * 167e-6 / 1e-9;
    x_data = repmat({animalID},length(y_data(:)),1);
    yarr = [yarr ; y_data(:)];
    xarr = [xarr; x_data];
end

violins = violinplot(yarr,xarr,'Bandwidth',0.1);

for i = 1:5
    violins(i).BoxColor = [0.2 0.2 0.2];
    violins(i).EdgeColor = [0.2 0.2 0.2];
    violins(i).ViolinAlpha{1} = 0.3;
end
colororder({'k','k'})

yyaxis left;
ylim([0 6])
ylabel('Thresholds (nC/phase)')
yyaxis right;
ylim([0 35])
ylh = ylabel('Thresholds ({\mu}A)');
set(ylh,'rotation',-90,'VerticalAlignment','bottom');

set(gca,'fontsize', 12)
set(gca,'fontweight','bold')
set(gca,'linewidth',1.5)

%% Figure 4F
% Longitudinal plot with error bars

figure;
mu_arr = {};
med_arr = {};
rsq = zeros(nANIMALS,1);
for i = 1:nANIMALS
    dateArr = s_new(i).allData.dateArr;
    daysArr = s_new(i).allData.daysArr;
    nSESSIONS = numel(dateArr);

    groupedData = cell2mat(s_new(i).allData.Ldata(3:5));
    mu_sig = zeros(1,nSESSIONS);
    med_sig = zeros(1,nSESSIONS);
    if i == 3
        nSESSIONS = 39; % hardcoded to only include pre-shock segment
    end
    for j = 1:nSESSIONS
        thresholdnC = groupedData(:,j)* 1e-6 * 167e-6 / 1e-9;
        [mu,CI] = calculateCI(thresholdnC,67);
        lowLim = CI(1); upLim = CI(2);
        if lowLim < 0
            lowLim = 0;
        end
        sd = std(thresholdnC,'omitnan');
        mu_sig(j) = mu;
        med_sig(j) = median(thresholdnC,'omitnan');
        
        e = errorbar(daysArr(j), mu, mu-lowLim, upLim-mu ,"LineStyle","none");
        e.Color = default_c(i,:);
        e.CapSize = 5;
        alpha = 0.4;
        set([e.Bar, e.Line], 'ColorType', 'truecoloralpha', 'ColorData', [e.Line.ColorData(1:3); 255*alpha])  
        hold on;
        scatter(daysArr(j), mu, 20, default_c(i,:),'filled','MarkerFaceAlpha',0.4);
        hold on;
    end

    if i == 3
        [fitresult, gof, curve_x, curve_y] = createExpFit(daysArr(1:39), mu_sig(1:39));
    else
        [fitresult, gof, curve_x, curve_y] = createExpFit(daysArr, mu_sig);
    end

    X = 0:curve_x(end);
    Y = fitresult(X);
    
    plot(X,Y,'color',default_c(i,:),'LineWidth',2)
end

% Make figure look pretty
xlim([0 230])
set(gca,'fontsize', 12)
set(gca,'fontweight','bold')
set(gca,'linewidth',2)
xlabel('Day of behavioral training')
ylabel('Threshold (nC/phase)')
yyaxis right
ylabel('Threshold ({\mu}A)')
ylim([0 36])
yticks(0:5:35)
yticklabels(string((0:5:35)))
yyrightlabel = get(gca,'YLabel');
set(yyrightlabel,'rotation',-90,'VerticalAlignment','baseline','HorizontalAlignment','center')

ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';

box on
ax = gca;
ax.LineWidth = 1.5;

% Legend
for i = 1:5
    h(i) = scatter(nan,nan,20,default_c(i,:),'filled','MarkerFaceAlpha',1);
end
lgd = legend(h,{'Mouse1','Mouse2','Mouse3','Mouse4','Mouse5'});
