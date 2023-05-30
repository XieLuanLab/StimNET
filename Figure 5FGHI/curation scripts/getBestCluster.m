function [cluster,amp,p2p,noise,snr,nClu] = getBestCluster(chK,minWvf,chMAD)

[~,score,~,~,~] = pca(chK-mean(chK));
X = score(:,1:2);
evaluation = evalclusters(X,"kmeans","silhouette","KList",1:2); % 2 works best 
clusterID = evaluation.OptimalY; % gives the cluster assignments with optimal K 
optK = max(evaluation.InspectedK);

centroids = nan(optK,48);
amp = zeros(optK,1);
p2p = zeros(optK,1);
noise = zeros(optK,1);
snr = zeros(optK,1);
nClu = zeros(optK,1);
for c = 1:optK
    cluster_k = chK(clusterID == c,:);
    centroid = mean(cluster_k,1);
    centroids(c,:) = centroid;
    amp(c) = abs(centroid(11));
    p2p(c) = max(centroid) - centroid(11);
    noise(c) = 1.4826*chMAD;
    snr(c) = amp(c) / noise(c);
    nClu(c) = size(cluster_k,1);
    % exclusion criteria
    if nClu(c) < minWvf 
        centroids(c,:) = nan(1,48);
        noise(c) = 0;
        p2p(c) = 0;
        snr(c) = 0;
        amp(c) = 0;
    end
end

%% Output

if sum(amp) > 0
    [~,winnerIdx] = max(amp);
    cluster = chK(clusterID(winnerIdx),:);
    amp = amp(winnerIdx);
    p2p = p2p(winnerIdx);
    noise = noise(winnerIdx);
    snr = snr(winnerIdx);
    nClu = nClu(winnerIdx);    
else
    cluster = [];
    amp = nan;
    p2p = nan;
    noise = nan;
    snr = nan;
    nClu = nan;
end


end