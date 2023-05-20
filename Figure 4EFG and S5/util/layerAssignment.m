% Get channels, their depths, and layer assignments 

close;
default_c = get(gca,'ColorOrder'); close;
default_colors = default_c(1:6,:);
animalIDs = {'ICMS24','ICMS42','ICMS37','ICMS44','ICMS49'}; 
nANIMALS = numel(animalIDs);

chList = cell(nANIMALS,1);
depthList = cell(nANIMALS,1);
for i = 1:nANIMALS
    animalID = animalIDs{i};
    data = s(i).data;
    [intan2depth,depth_per_ch,nELECS_OUT] = getDepthFromIntanCh(data.intanCh,animalID);
    chList{i} = intan2depth(:,1);
    depthList{i} = intan2depth(:,2);
end

%% Laminar values obtained from Allen brain atlas 
borderData = cell(5,1);
% Mouse1 laminar border value assignments in M2
l1 = 100; l23 = 550; l5 = 900; l6 = 1400;
M1 = [l1 l23 l23 l5 l6];
% S1 laminar border values
% Mouse2
l1 = 100; l23 = 300; l4 = 450; l5 = 750; l6 = 1200;
M2 = [l1 l23 l4 l5 l6];
% Mouse3
l1 = 100; l23 = 400; l4 = 500; l5 = 650; l6 = 1200;
M3 = [l1 l23 l4 l5 l6];
% Mouse4
l1 = 130; l23 = 600; l4 = 900; l5 = 1400; l6 = 1600;
M4 = [l1 l23 l4 l5 l6];
% Mouse5
l1 = 130; l23 = 390; l4 = 500; l5 = 1000; l6 = 1400;
M5 = [l1 l23 l4 l5 l6];

borderData{1} = M1;
borderData{2} = M2;
borderData{3} = M3;
borderData{4} = M4;
borderData{5} = M5;
%%
layerAssignments = cell(5,1);
layersOnly = cell(5,1);
for i = 1:5
    b = borderData{i};
    chs = chList{i};
    depths = depthList{i};
    layeridxs = zeros(numel(depths),1);
    for j = 1:numel(depths)
        depth = depths(j);
        ch = chs(j);
        if depth <= b(1)
            layeridxs(j) = 1;
        elseif (depth > b(1)) && (depth <= b(2))
            layeridxs(j) = 2;
        elseif (depth > b(2)) && (depth <= b(3))
            layeridxs(j) = 3;
        elseif (depth > b(3)) && (depth <= b(4))
            layeridxs(j) = 4;
        elseif (depth > b(4)) 
            layeridxs(j) = 5;
        end
    end
    layerAssignments{i} = [chs depths layeridxs];
    layersOnly{i} = layeridxs';
end
%% Clear variables
clear M1 M2 M3 M4 M5 nELECS_OUT j l1 l23 l4 l5 l6 borderData b ch chList chs depth depth_per_ch i animalID