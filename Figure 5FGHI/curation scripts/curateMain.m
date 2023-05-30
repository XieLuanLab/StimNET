load('startDates.mat');
addpath '\\10.129.151.108\xieluanlabs\xl_stimulation\Robin\Ripple32Ch'
addpath '\\10.129.151.108\xieluanlabs\xl_stimulation\Robin\Manuscript figures\Figure 5\Spike analysis'
load('\\10.129.151.108\xieluanlabs\xl_stimulation\Robin\Manuscript figures\Figure 5\Spike analysis\stimChannels.mat')
color_arr = get(gca,'colororder');
close;

%% Data curation
animalIDs = {'ICMS24','ICMS42','ICMS37','ICMS44','ICMS49'};
maxWeeks = 100;
weekArrs = cell(maxWeeks,5);
for i = 1:length(animalIDs)
    animalID = animalIDs{i};
    file = dir(strcat(animalID,'*.mat'));
    f = file(end).name;
    fprintf('Loading variables from %s...\n',file(end).name);
    load(f,'dates','wvf_date_arr','MADArr','segmentLens','ripple_ch');
    
    startDate = startDates.(animalID);
    curate;
    
    weekArr = convertDates2WPI(dates,i);
    column = groupWvfsByWeek(weekArr,curated_wvfs,maxWeeks);
    weekArrs(:,i) = column;

    curatedData.(animalID).dates = dates; 
    curatedData.(animalID).rippleChannels = ripple_ch;
    curatedData.(animalID).intanChannels = Ripple2Intan(ripple_ch);
end

save('curatedData.mat',"curatedData")

%% Plot spikes 
plot_weeks = [8 24 42; 9 29 40; 5 21 42; 5 21 37; 5 18 31]; 
channels2plot = [18,23,2,29,22]; 
totalPulses = [925000 1.9e6 106000 171946 10100];

ylimmax = 80;
ylimmin = -300;
tiledlayout(5,3,'TileSpacing','none');

spikeArr = cell(5,3);

for row = 1:5
    animalID = animalIDs{row};
    channel2plot = channels2plot(row);
    intanChannels = curatedData.(animalID).intanChannels;
    chIdx = find(intanChannels == channel2plot);
    plotW = plot_weeks(row,:);
    weekArr = weekArrs(:,row);
    for col = 1:3
        nexttile;
        x = (1:48)/3e1 + (col-1)*2;
        spikes = weekArr{plotW(col)}{chIdx};
        spikeArr{row,col} = spikes;
        h=stdshade(spikes',0.3, color_arr(row,:), x, 1, 1);

        xlim([x(1)-1 x(end)+1])
        ylim([ylimmin ylimmax]);

        hold on;
        axis off
        if row == 5 && col == 2
            xline(-1/3e1)
            line([0 20/3e1], [-200 -200])
            line([x(1)+25/3e1 x(1)+30/3e1], [0 0])
            line([x(1)+25/3e1 x(1)+30/3e1], [-200 -200])
        end
        % scale bars
        if row == 5 && col == 1
            xline(-1/3e1)
            line([0 30/3e1], [-200 -200])
        end
    end
end

t=tiledlayout(5,3,'TileSpacing','none');
% 
% for i = 1:5
%     for j = 1:3
%         nexttile;
%         plot(spikeArr{i,j}','k','LineWidth',0.1)
%     end
% end