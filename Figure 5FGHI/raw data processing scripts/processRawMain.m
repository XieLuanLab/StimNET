% Gets threshold crossing events and other metrics from raw ephys data 

addpath '\\10.129.151.108\xieluanlabs\xl_stimulation\Robin\Ripple32Ch'
addpath '\\10.129.151.108\xieluanlabs\xl_stimulation\Robin\Manuscript figures\Figure 5\Spike analysis'
addpath '\\10.129.151.108\xieluanlabs\xl_stimulation\Robin\lib\neuroshare'
load('\\10.129.151.108\xieluanlabs\xl_stimulation\Robin\Manuscript figures\Figure 5\Spike analysis\stimChannels.mat')

name = getenv('COMPUTERNAME'); 
if strcmp(name, 'BRC8065DGV513')
    data_dir = 'D:\ICMS';
else
    data_dir = 'E:\Robin\Behavioral';
end
d = dir(data_dir);
animalID = {'ICMS24','ICMS42','ICMS37','ICMS44','ICMS49'};

MADArrs = cell(5,1);
segmentLenArr = cell(5,1);
datesArr = cell(5,1);
bigWvfArr = cell(5,1);
P2PArrs = cell(5,1);

for i = 1:numel(animalID)
    fprintf('\nProcessing animalID %s...\n\n',animalID{i})
    animalFolder = fullfile(data_dir,animalID{i});
    % get Ripple stim channels 
    ripple_ch = intan2ripple(S.(animalID{i}));

    a = dir(animalFolder);
    a = a(~ismember({a(:).name},{'.','..','other','psychometric figures','staircase_figures','old','frequency','stage1Figures'}));
    dates = {a.name};
    [~, idx] = sort(datenum(dates, 'dd-mmm-yyyy'), 1, 'ascend');
    dates = dates(idx); datesArr{i} = dates;
    
    nDATES = numel(dates);
    wvf_date_arr = cell(nDATES,1);
    P2P_date_arr = cell(nDATES,1);
    segmentLens = nan(nDATES,1);
    MADArr = cell(nDATES,1);
    good_folder_idx = zeros(nDATES,1);
    
    for j = 1:nDATES
        date_folder = fullfile(animalFolder,dates{j});
        fns5 = dir(fullfile(date_folder,'*.ns5'));     
        fnev = dir(fullfile(date_folder,'*.nev'));    
        if isempty(fns5) || isempty(fnev)
            warning('File path for %s not found!',date_folder);
            continue;
        end
        if length({fns5.name}) > 1
            ns5_file_path = fullfile(date_folder,fns5(1).name);
            nev_file_path = fullfile(date_folder,fnev(1).name);
        else 
            ns5_file_path = fullfile(date_folder,fns5.name);
            nev_file_path = fullfile(date_folder,fnev.name);
        end
        
        fprintf('Processing folder %s...\n',dates{j});
        [wvf_ch_arr,p2p_ch_arr,preStimLength, MADval] = getThresholdCrossings(ns5_file_path,nev_file_path,ripple_ch);
        if isempty(MADval)
            warning('Ephys data for date folder %s is not suitable for analysis!\n',date_folder);
            wvf_date_arr{j} = nan;
            P2P_date_arr{j} = nan;
            MADArr{j} = nan;
        else
            wvf_date_arr{j} = wvf_ch_arr;
            P2P_date_arr{j} = p2p_ch_arr;
            good_folder_idx(j) = 1;
            segmentLens(j) = preStimLength/3e4;
            MADArr{j} = MADval;
        end
        MADArrs{i} = MADArr;
        segmentLenArr{i} = segmentLens;
        bigWvfArr{i} = wvf_date_arr;
        P2PArrs{i} = p2p_ch_arr;
    end
end
