% PCA 
warning('off')
plotC = [0 0 0 0.25];
minWvf = 60;
w = wvf_date_arr;
nCH = numel(ripple_ch);
nDATES = numel(dates);
curated_nActiveCh = nan(nDATES,1);
curated_wvfs = cell(nDATES,1);
curated_p2p = cell(nDATES,1);
curated_noise = cell(nDATES,1);
curated_SNR = cell(nDATES,1);
curated_amp = cell(nDATES,1);
for j = 1:nDATES
    dateJ = w{j};
    if size(dateJ,1) == 1 || isempty(dateJ)
        curated_wvfs{j} = [];       
        continue;
    end
    dateCell = cell(nCH,1);
    dateP2P = zeros(nCH,1);
    dateNoise = zeros(nCH,1);
    dateSNR = zeros(nCH,1);
    dateAmp = zeros(nCH,1);
    dateMAD = MADArr{j};
    for k = 1:nCH
        chMAD = dateMAD(k);
        chK = dateJ{k};
        chK(any(isnan(chK), 2), :) = [];
        if size(chK,1) < 4
            dateCell{k} = [];
            dateP2P(k) = nan;
            dateNoise(k) = nan;
            dateSNR(k) = nan;
            dateAmp(k) = nan;
            continue
        end
        
        [cluster,amp,p2p,noise,snr,nClu] = getBestCluster(chK,minWvf,chMAD);
        dateCell{k} = cluster;
        dateP2P(k) = p2p;
        dateNoise(k) = noise;
        dateSNR(k) = snr;
        dateAmp(k) = amp;
    end
    curated_nActiveCh(j) = sum(~cellfun(@isempty,dateCell));
    curated_wvfs{j} = dateCell;
    curated_p2p{j} = dateP2P;
    curated_noise{j} = dateNoise;
    curated_SNR{j} = dateSNR;
    curated_amp{j} = dateAmp;
end

