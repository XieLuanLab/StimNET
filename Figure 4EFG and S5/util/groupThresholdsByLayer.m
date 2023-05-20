%% Create layer-dependent threshold cell arrays grouped by animalID and layer 
s_new = s;
animalGroup = cell(nANIMALS,1); % grouped by animalID
layerGroup = cell(5,1); % grouped by layer 
for i = 1:nANIMALS
    L1data = [];
    L23data = [];
    L4data = [];
    L5data = [];
    L6data = [];

    animalID = s(i).animalID;
    data = s(i).data;
    thresholdArr = data.thresholdArr; % behavioral thresholds for a session 
    testedChArr = data.testedChArr; % channels tested in that session
    dateArr = data.dateArr; % session dates
    daysArr = data.daysArr; % session dates converted to indices
    intanCh = data.intanCh; % Intan channel IDs of stim channels 
    
    layerData = layerAssignments{i};

    assert(isequal(intanCh,layerData(:,1)));
    for j = 1:size(layerData,1)
        thresholdData = thresholdArr(j,:);
        nanRemoved = thresholdData(~isnan(thresholdData));
        layerIdx = layerData(j,3);
        switch layerIdx
            case 1
                L1data = [L1data ; thresholdData];
                layerGroup{1} = [layerGroup{1} nanRemoved];
            case 2
                L23data = [L23data ;thresholdData];
                layerGroup{2} = [layerGroup{2} nanRemoved];
            case 3
                L4data = [L4data ;thresholdData];
                layerGroup{3} = [layerGroup{3} nanRemoved];
            case 4
                L5data = [L5data ;thresholdData];
                layerGroup{4} = [layerGroup{4} nanRemoved];
            case 5
                L6data = [L6data ;thresholdData];
                layerGroup{5} = [layerGroup{5} nanRemoved];
        end
        
    end

    Ldata = {L1data;L23data;L4data;L5data;L6data};
    
    % Save into new structure
    s_new(i).allData.thresholdArr = thresholdArr;
    s_new(i).allData.testedChArr = testedChArr;
    s_new(i).allData.dateArr = dateArr;
    s_new(i).allData.daysArr = daysArr;
    s_new(i).allData.intanCh = intanCh;
    s_new(i).allData.Ldata = Ldata;
end

%% Group by layer
layerGroupedArr = cell(5,1);
for i = 1:5
    layerData = s_new(i).allData.Ldata;
    for l = 1:5
        temp = layerData{l};
        layerGroupedArr{l} = [layerGroupedArr{l} temp(:)'];
        temp2 = layerGroupedArr{l};
        layerGroupedArr{l} = temp2(~isnan(temp2));
    end
end
%%
clear L1data L23data L4data L5data L6data layerIdx j good_ch nanRemoved
