function [P2P_arr,noise_arr,SNR_arr,amp_arr] = groupDataByWeek(weekArr,p2p,SNR,noise)
    
    uniqueWeeks = unique(weekArr);
    P2P_arr = cell(numel(uniqueWeeks),1);
    noise_arr = cell(numel(uniqueWeeks),1);
    SNR_arr = cell(numel(uniqueWeeks),1);
    amp_arr = cell(numel(uniqueWeeks),1);
    for i = 1:numel(uniqueWeeks)
        week = uniqueWeeks(i);
        rows = find(weekArr == week);
        P2P_arr{i} = [p2p{rows}];
        noise_arr{i} = [noise{rows}];
        SNR_arr{i} = [SNR{rows}];
    end

end