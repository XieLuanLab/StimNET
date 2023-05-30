function [wvf_ch_arr, p2p_ch_arr, preStimLength, MADval] = getThresholdCrossings(ns5_file_path,nev_file_path,ripple_ch)

stim_times_arr = readNEV_2(nev_file_path); 
if isempty(stim_times_arr)
    wvf_ch_arr = []; p2p_ch_arr = []; preStimLength = []; MADval = [];
    return
end
stimTimes = sort(cell2mat(stim_times_arr(:,2)));
firstStim = stimTimes(1)-1;
% First stimulation pulse time is variable.
% To remove potential confound of stimulation on spiking activity, try to
% record before first stimulation pulse if first stim pulse is after 60s.
% If before 60s, record up to 60s anyways. 
preStimLength = max(firstStim,60*3e4);

raw = readNS5_2(ns5_file_path, preStimLength);
if isempty(raw)
    wvf_ch_arr = []; p2p_ch_arr = []; preStimLength = []; MADval = [];
    return
end

% Remove common-mode signal 
cmr = median(raw,2);
raw = raw - cmr;
if size(raw,1) < 32
    wvf_ch_arr = []; p2p_ch_arr = []; preStimLength = []; MADval = [];
    return
end
raw = raw(ripple_ch,:);

%% Filter, remove movement artifacts, and blank 
copy = raw;
filtSig = bandpass(copy',[300 5000],3e4);
nSCH = size(copy,1);
filtSig(1:30,:) = zeros(30,nSCH);

% find artifacts using arbitrary thresholds
[r,~] = find(filtSig > 100 | filtSig < -1000);
artIdx = unique(r);

% blank 1 ms before and after artifact 
for i = 1:numel(artIdx)
    preArt = artIdx(i) - 30;
    postArt = artIdx(i) + 30;
    filtSig(preArt:postArt,:) = zeros(61,nSCH);
end

%% Get rms, MAD

rmsVal = rms(filtSig(1:firstStim,:));
MADval = 1.4826*mad(filtSig(1:firstStim,:));

%% Clear copy and raw to free up memory 
clear raw

%% Get threshold crossing events
spikeThresholds = rmsVal * 4.5; % threshold
MinPeakDistance = 45; % 1.5 ms
MaxPeakWidth = 25;

nSCH = size(filtSig,2);

wvf_ch_arr = cell(nSCH,1);
p2p_ch_arr = cell(nSCH,1);

for ch_idx = 1:nSCH
    seg = filtSig(:,ch_idx);
    threshold = spikeThresholds(ch_idx);
    segLength = length(filtSig);

    [pk, ts] = findpeaks(-seg,'MinPeakHeight',threshold,'MaxPeakWidth',MaxPeakWidth,...
        'MinPeakProminence',40,'MinPeakDistance',MinPeakDistance,'MinPeakWidth',5,'WidthReference','halfprom');
    
    ts = ts(pk < 300); % reject events whose peak is over 300 uA
    ts = ts(ts - 10 > 0 & ts + 37 < segLength); 
      
    wvf_arr = nan(numel(ts),48);
    p2p_arr = nan(numel(ts),1);
    for evt = 1:numel(ts)
        wvf = seg(ts(evt)-10:ts(evt)+37);
        prePeak = wvf(1:5); % if all zeros, indicates threshold crossing is due to post-blanking artifact
        last10 = wvf(end-9:end); 

        if max(wvf) < 200 && min(wvf) == -pk(evt) && any(prePeak) && any(last10)
            wvf_arr(evt,:) = wvf;
            p2p_arr(evt) = max(wvf) - wvf(11);
        end
    end
    wvf_ch_arr{ch_idx} = wvf_arr;
    if isempty(wvf_arr)
        p2p_ch_arr{ch_idx} = nan;
    else
        p2p_ch_arr{ch_idx} = p2p_arr;
    end
end

end