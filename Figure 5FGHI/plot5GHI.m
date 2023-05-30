% Plot Figure 5GHI metrics

load('data/metrics.mat');
addpath('util')
color_arr = get(gca,'colororder');
plotColor = color_arr(2,:);
close all
for figID = 1:3
    figure  
    startDate = metrics.stage2StartDate;
    dates = metrics.dates; 

    idx1 = 1:length(dates);
    idx2 = find(cellfun('length',dates) == 11);
    idx3 = find(~contains(dates,'2023'));

    stage2Start = find(contains(dates,startDate));
    goodDatesIdx = intersect(intersect(idx1,idx2),idx3);
    dates = dates(goodDatesIdx); nDATES = numel(dates);
    p2p = metrics.p2p(goodDatesIdx); 
    noise = metrics.noise(goodDatesIdx);
    SNR = metrics.SNR(goodDatesIdx); 
    nCH = metrics.nCH;
    nActiveCh = metrics.nActiveCh(goodDatesIdx); 

    weekArr = convertDates2WPI(dates,2);

    [P2P_arr,noise_arr,SNR_arr] = groupDataByWeek(weekArr,p2p,SNR,noise);

    switch figID
        case 1; data = P2P_arr; ylbl = 'Peak-to-peak amplitude (\muV)'; titlestr = 'P2P voltage';
        case 2; data = SNR_arr; ylbl = 'SNR'; titlestr = 'SNR';
        case 3; data = noise_arr; ylbl = 'Noise (\muV)'; titlestr = 'Noise';
    end        

    plotMetricWPI(weekArr,data,plotColor);
    yl = get(gca, 'YLim');
    ylabel(ylbl);
    box on
    set(gca,'fontsize', 15);
    set(gca,'linewidth',1.5);
    xlim([0 281])
    limsy=get(gca,'YLim');
    ylim([0 limsy(2)*1.1])
end