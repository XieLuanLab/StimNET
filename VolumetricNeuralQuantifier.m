% VolumetricNeuralQuantifer analyzes sequences of volumetric image 2P scans
% to quantify the neural activation trends in response to neural
% stimulation during paired stimulation. 
%
% Usage:
%  Current version of program requires specific folder stucture and source
%  variable 'sourceFolders' listing files to process
%
% - Roy Lycke (rjl6@rice.edu)
% Released 1.0  Date: 04/29/2023

warning('off','all')
warning

flipFlag =0; % Flips the order of baseline and stim file, necissary for older datasets which had flipped experimental settings

% microns_per_pixel= 2.66;%2.232; %2.66;
probThresh = 0.75;
runDistance = 0;
elecPos = zeros(32,1); % Possible electrode position file
warning('off','MATLAB:MKDIR:DirectoryExists');

javaaddpath 'C:\Program Files\MATLAB\R2019a\java\mij.jar'
javaaddpath 'C:\Program Files\MATLAB\R2019a\java\ij.jar'


% Load supporting files
addpath(genpath('./rasterplot'));
addpath(genpath('./natsortfiles'));
addpath(genpath('./imshow3D'));
addpath(genpath('./bfmatlab'));
addpath(genpath('./NeuronAnalysis'));
addpath(genpath('./distinguishable_colors'));
addpath(genpath('./strel3d'));
addpath(genpath('./mij'));


% Start ImageJ interface
%    If you need to load the libraries look here: imagej.net/plugins/miji
% Update file location to your installation of FIJI
MIJ.start('C:\Users\XXXX\Desktop\fiji-win64\Fiji.app\plugins','C:\Users\XXXX\Desktop\fiji-win64\Fiji.app\plugins',0);



for curFile = 1:numel(sourceFolders)
    
    selpath = sourceFolders{curFile};
%     selpath = uigetdir('D:\','Select folder containing raw files');
    
    idcs = strfind(selpath,'\');
    pullPath = selpath(1:idcs(end)-1);
    selpathOut = strcat(pullPath,'\OUTPUT_',datestr(datetime('now'),'mm-dd-yyyy_HH-MM'),'_VER2July');
    keypath = strcat(pullPath,'\Key');
    segListing = dir(selpath);
    segListing(1:2) = []; % remove directory entries

    
    szFilePath = strcat(pullPath,'\ImgSize.mat');
    if(isfile(szFilePath))
        szVar = load(szFilePath);
        fov = szVar.imgSize;
        microns_per_pixel = fov/512; % all images collected at 512 resolution
    else
        microns_per_pixel = 2.18; % Assume zoom of 1 and FOV of 1120 um
    end
    

    
    % Check if electrode posistion file exisits
    elecPosFile = strcat(pullPath,'\ElecPos.mat');
    if(isfile(elecPosFile))
        runDistance=1;
        distData = load(elecPosFile);
        elecPos = distData.contactPos;
    end
    
    selpathOutOrig = selpathOut;
    fprintf('\nProcessing Dataset: ');
    disp(selpath)
    
    % Parse the vertical step distance from the file header. Could cause
    % issues if naming is shifted around too much
    firstFile = segListing(1).name;
    fragments = split(firstFile,'_');
    stepStr = split(fragments{3},'um');
    vertStep = str2double(stepStr{1});

    % Generate temporary data folder (accelerate image access)
    tempPath = strcat(selpathOut,'\TEMP');
    mkdir(tempPath)




    %% Extract listing of channels and currents used in experiment
    keyListing = dir(keypath);
    keyListing(1:2) = []; % remove directory entries
    keyNames = zeros(numel(keyListing),1);
    arrSize = 0;
    for i = 1:numel(keyListing)
        temp = keyListing(i).name;
        C = strsplit(temp,'_');
        C = strsplit(C{2},'.');
        keyNames(i) = str2double(C{1});
    end
    [~,keyInds] = sort(keyNames);

    

    
    
    tempMast = cell(numel(keyListing),1);
    tempPair = cell(numel(keyListing),1);
    for i = 1:numel(keyListing)
        currents_used = open(strcat(keypath,'\',keyListing(keyInds(i)).name));
        curKey = currents_used.replicant;
        curFlag = zeros(size(curKey,1),1);
        maxInd = 1;
        for k = 1:size(curKey)
            if(~ischar(curKey{k,1}))
                maxInd = 2;
            end
        end
        if(size(curKey,2)==4) % Paired stimulation
            temp=cell(size(curKey,1),2);
            temp(:,2)= curKey(:,4);
            for j = 1:size(curKey,1)
                chan = strcat(curKey{j,2},curKey{j,3});
                if(isequal(chan,'00'))
                    chan = '0';
                end
                temp{j,1}=chan;
            end
            curKey = temp;
            curFlag = ones(size(curKey,1),1);  
        elseif(maxInd==2)
            for k = 1:size(curKey)
                if(~ischar(curKey{k,1}))
                    curFlag(k)=1;
                end
            end
        end
        tempMast{i} =curKey(:,1:2);
        tempPair{i} = curFlag;
        arrSize = arrSize + numel(curFlag);
    end

    %Transfer key data to arrays
    masterKey = cell(arrSize,2);
    pairFlag = false(arrSize,1);
    c = 1;
    for i = 1:numel(keyListing)
        curPair = tempPair{i};
        curMast = tempMast{i};
        for j = 1:numel(curPair)
            pairFlag(c) = curPair(j);
            masterKey{c,1} = curMast{j,1};
            masterKey{c,2} = curMast{j,2};
            c = c+1;
        end
    end
    % Flatten the cell arrays and ensure only one layer deep
    c = 1;
    flatKey = cell(size(masterKey));
    for i = 1:size(masterKey,1)
        temp1 = masterKey{i,1};
        temp2 = masterKey{i,2};
        if(size(temp1,1)==1)
            flatKey{c,1}= temp1;
            if(isequal(temp2,'0'))  
                flatKey{c,2}=0;
            else
                flatKey{c,2}= temp2;
            end
            c = c+1;
        else
            for k = 1:size(temp1)
               flatKey{c,1}= temp1{k};
               flatKey{c,2}= temp2{k};
               c = c+1;
            end
        end
    end
    clear temp1 temp2

    % Pull data from temporally ordered key
    scanCurrent = cell2mat(flatKey(:,2)); 
    scanCurrent(scanCurrent==0.567) = 0.33; % Adjust current values to approximate ones
    scanCurrent(scanCurrent==0.533) = 0.16;
    scanChannel = flatKey(:,1);
    for i = 1:size(scanChannel,1)
        temp = scanChannel{i};
        if(isstring(temp) || ischar(temp))
            continue
        else
            scanChannel{i} = strcat(temp{1},temp{2});
        end
    end



    % provide numerical tags for unique groups
    [chNames,~,~]= unique(scanChannel);
    chNames(strcmp(chNames,'0'))=[];

    currNames = unique(scanCurrent);
    currNames(currNames==0)=[];


    baseScans = zeros(size(scanCurrent));
    baseScans(scanCurrent==0)=1;
    finalKey = flatKey(baseScans==0,:); % Remove baseline scans & Update scan key
    scanCurrent = scanCurrent(baseScans==0,:); 
    scanChannel = scanChannel(baseScans==0,:);

    % ensure first element in key is a string
    for i = 1:size(finalKey,1)
        if(~ischar(finalKey{i,1}))
            temp = finalKey{i,1};
            finalKey{i,1} = strcat(temp{1},temp{2});
        end
    end

    % Rearrange raw folders into temporal order - Key files must have Rep and
    % then an Underscore in their name
    rawScanID = zeros(numel(segListing),1);
    for i = 1:numel(segListing)
        filep = strsplit(segListing(i).name,'_');
        ind = strfind(filep,'Rep');% Find source file
        ind2 = strfind(filep,'REP');% Find source file
        ind3 = strfind(filep,'rep');% Find source file
        tInd=0;
        for j = 1:numel(ind)
            if(~isempty(ind{j}) && ind{j})
                tInd = j;
            end
            if(~isempty(ind2{j}) && ind2{j})
                tInd = j;
            end
            if(~isempty(ind3{j}) && ind3{j})
                tInd = j;
            end
        end

        repInfo = filep{tInd};
        temp = strsplit(repInfo,'-');
        if(contains(temp{2},'TIFF'))
            temp2 = strsplit(temp{2},'T');
            rawScanID(i) = str2double(temp2{1});
        else
            rawScanID(i) = str2double(temp{2});
        end
    end
    [~,order] = sort(rawScanID);
    segListing = segListing(order); % Now the seg listing should be ordered temporally


    %% Load raw data

    % Only open folders containing tiff converted files
    toKeep = false(numel(segListing),1);
    for i = 1:numel(segListing)
        if(contains(segListing(i).name,'TIFF'))
            toKeep(i)=1;
        end
    end
    segListing = segListing(toKeep); 





    % Iterate through all folders and compile a list of all scan files in order
    % of aquisition so they can be pulled via scan index
    SegFileListing = cell(numel(pairFlag),1);
    c = 1;
    for j = 1: numel(segListing)
        curFolder = strcat(selpath,'\',segListing(j).name);
        subFolderListing = dir(curFolder);
        subFolderListing(1:2) = []; % remove directory entries


        subListing = cell(numel(subFolderListing),1);
        for k = 1:numel(subFolderListing)
            subListing{k} = subFolderListing(k).name;
        end

        %Order all the scans in order of original collection 
        rawScanID = zeros(numel(subListing),1);
        for k = 1:numel(subListing)
            [~,fileName,~] = fileparts(subListing{k});
            filep = strsplit(fileName,'_');
            orderTag = filep{numel(filep)};
            temp = strsplit(orderTag,'-');
            rawScanID(k) = str2double(temp{2});
        end
        [~,order] = sort(rawScanID);
        subListing = subListing(order); % Now the seg listing should be ordered temporally
        for k = 1:numel(subListing)
            SegFileListing{c} = strcat(curFolder,'\',subListing{k});
            c = c + 1;
        end
    end

















    %% Load data - Data will be loaded by channel/current group
    c = 0;
    baseMean = zeros(numel(chNames),numel(currNames));
    baseSTD = zeros(numel(chNames),numel(currNames));
    stimMean = zeros(numel(chNames),numel(currNames));
    stimSTD = zeros(numel(chNames),numel(currNames));
    probabilityMap = cell(numel(chNames),numel(currNames));
    probabilityMapVol = cell(numel(chNames),numel(currNames));
    probabilityMapInhibit = cell(numel(chNames),numel(currNames));
    activeMap = cell(numel(chNames),numel(currNames));
    activeMapInhibit = cell(numel(chNames),numel(currNames));
    IntensityMap = cell(numel(chNames),numel(currNames));
    subtract_individual = cell(numel(chNames),numel(currNames));
    segment_individual = cell(numel(chNames),numel(currNames));
    IntensityMapInhibit = cell(numel(chNames),numel(currNames));
    baseMeanPure = zeros(size(scanCurrent));
    baseSTDPure = zeros(size(scanCurrent));


    % Load Initial scan volume to define electrode location
    opener = ij.io.Opener();
    opener.openImage(SegFileListing{1}).show();
    % perform gaussian blurring on the current volume
    MIJ.run('Measure Stack...', 'channels slices frames order=czt(default)');
    maxChan = max(MIJ.getColumn('Ch'));
    MIJ.run("Clear Results");
    if(maxChan==2)
       MIJ.run('Duplicate...', 'duplicate channels=2');
   elseif(maxChan==3)
       MIJ.run('Duplicate...', 'duplicate channels=3'); % only get the last channel in the scan 
    end
    
    % Push the data directly from ImageJ into Matlab
    frameBlockRaw = double(MIJ.getCurrentImage());

    

    
    % Load all scans to generage scanwide metrics
    runningMean = zeros(size(frameBlockRaw));
    runningM2 = zeros(size(frameBlockRaw));
    runningN = 0;
    
    meanBaseVol = zeros(size(frameBlockRaw));
    allMeans = zeros(numel(SegFileListing),1);
    allDiff = zeros(numel(SegFileListing)/2,1);
    allDiffSTD = zeros(numel(SegFileListing)/2,1);
    allDiffPos = zeros(numel(SegFileListing)/2,1);
    allDiffSTDPos = zeros(numel(SegFileListing)/2,1);
    allMax = zeros(numel(SegFileListing),1);
    matFileListing = cell(size(SegFileListing));
    rollFileListing = cell(size(SegFileListing));
    fprintf('Performing Initial Intensity Sweep... ');
    strCR = -1;
    pause(0.1);
    
    diffMeanInd = zeros(numel(SegFileListing)/2,1);
    diffSTDInd = zeros(numel(SegFileListing)/2,1);
    
    for i = 1:2:numel(SegFileListing)


        % Display processing status
        pVal = round(i/numel(SegFileListing)*100);     
        percentageOut = [num2str(pVal) '%%'];
        strOut = [percentageOut repmat(' ',1,10-length(percentageOut)-1)];
        fprintf([strCR strOut]);
        strCR = repmat('\b',1,length(strOut)-1);

        c = c+1;
        % Load paired baseliune scan volume
        opener = ij.io.Opener();
        if(flipFlag==0)
            opener.openImage(SegFileListing{i+1}).show();
        else
            opener.openImage(SegFileListing{i}).show();
        end
        if(maxChan==2) % only get the last channel in the scan 
            MIJ.run('Duplicate...', 'duplicate channels=2');
        elseif(maxChan==3)
            MIJ.run('Duplicate...', 'duplicate channels=3');
        end
        % Push the data directly from ImageJ into Matlab
        frameBlockB = double(MIJ.getCurrentImage());
                
        allMeans(i) = mean(frameBlockB,'all');
        allMax(i) = max(frameBlockB,[],'all');
        meanBaseVol = meanBaseVol + frameBlockB;

  
        delta = frameBlockB - runningMean;
        runningMean = runningMean + (delta / c);
        delta2 = frameBlockB - runningMean;
        runningM2 = runningM2 + delta .* delta2;
        
        
        
        opener = ij.io.Opener();
        if(flipFlag==0)
            opener.openImage(SegFileListing{i}).show();
        else
            opener.openImage(SegFileListing{i+1}).show();
        end
        if(maxChan==2) % only get the last channel in the scan 
            MIJ.run('Duplicate...', 'duplicate channels=2');
        elseif(maxChan==3)
            MIJ.run('Duplicate...', 'duplicate channels=3');
        end
        % Push the data directly from ImageJ into Matlab
        frameBlockS = double(MIJ.getCurrentImage());
                
        allMeans(i+1) = mean(frameBlockS,'all');
        allMax(i+1) = max(frameBlockS,[],'all');
        temp = frameBlockS - frameBlockB; % Excitation analysis
        % Calculate volume wide mean and std
        allVals = reshape(temp,1,[]);
        diffMeanInd(c) = mean(allVals);
        
        allDiff(c)= mean(temp,'all');
        allDiffSTD(c) = std(temp,0,'all');
        temp(temp<=0)=NaN;
        allDiffPos(c) = mean(temp,'all','omitnan');
        allDiffSTDPos(c) = std(temp,0,'all','omitnan');

        % Save 3D image volume into mat file for quicker retrieval in
        % the future
        combineFileName = strcat(tempPath,'\vol',num2str(i),'.mat');
        matFileListing{i} = combineFileName;
        save(combineFileName,'frameBlockB');  

        combineFileName = strcat(tempPath,'\vol',num2str(i+1),'.mat');
        matFileListing{i+1} = combineFileName;
        save(combineFileName,'frameBlockS');  

        try
            MIJ.run('Close'); 
        catch
        end

        try
            MIJ.closeAllWindows % Close figures
        catch
        end
    end
    fprintf('\n')        

    meanBaseVol = runningMean;
    varianceVol = runningM2./c;
%     sampleVarianceVol = runningM2./(c-1);
    stdBaseVolume = zeros(size(varianceVol));
    for l = 1:size(varianceVol,1)
        for m = 1:size(varianceVol,2)
            for n = 1:size(varianceVol,3)
                stdBaseVolume(l,m,n) = sqrt(varianceVol(l,m,n));
            end
        end
    end
    

%     diffMean = mean(diffMeanInd);
%     meanBaseVol = meanBaseVol/(numel(SegFileListing)/2);
%     maxRawInten = max(allMax) * 0.6; % Use 10% off of maximum for plotting figures
%     % Generate folder for unique imaging parametters
%     SourcePathOut = strcat(selpathOutOrig,'\SourceFiles');
%     mkdir(SourcePathOut)
%     stdThreshList = 0.5 :0.1: 2;


%% Calculate baseline volume properties and locations of vasculature
    % Plot and save the STD base Volume Plot
    curFig = figure(); % STD volume intensity
    imshow(squeeze(max(stdBaseVolume/(max(stdBaseVolume,[],'all')),[],3)))
    colormap(jet)
    colorbar
    title(strcat('Baseline STD volume : Max intensity = ',num2str(max(stdBaseVolume,[],'all'))));
    combineFileName = strcat(selpathOutOrig,'\BaselineSTDVol.fig');
    saveas(curFig,combineFileName);  
    combineFileName = strcat(selpathOutOrig,'\BaselineSTDVol.tiff');
    saveas(curFig,combineFileName);  
    close(curFig)
    combineFileName = strcat(selpathOutOrig,'\BaselineSTDVol.mat');
    save(combineFileName,'stdBaseVolume');
    
    
    curFig = figure(); % Mean Volume intensity
    set(0,'CurrentFigure',curFig)
    clf
    imshow(squeeze(max(meanBaseVol/(max(meanBaseVol,[],'all')),[],3)))
    colormap(jet)
    colorbar
    title(strcat('Baseline Mean volume : Max intensity = ',num2str(max(meanBaseVol,[],'all'))));
    combineFileName = strcat(selpathOutOrig,'\BaselineMeanVol.fig');
    saveas(curFig,combineFileName);  
    combineFileName = strcat(selpathOutOrig,'\BaselineMeanVol.tiff');
    saveas(curFig,combineFileName);  
    close(curFig)
    combineFileName = strcat(selpathOutOrig,'\BaselineMeanVol.mat');
    save(combineFileName,'meanBaseVol');
    
    
    % Locate regions likely to be vasculature ( low intensity, large, and
    % with low variance)
    overallMean = mean(meanBaseVol,'all');
    overallSTD = mean(stdBaseVolume,'all');
    vasculature = false(size(meanBaseVol));
    vasculature(meanBaseVol<(overallMean*0.45))=1; % low intensity regions
    vasculature(stdBaseVolume>(overallSTD*0.45))=0; % Low variability
    vasculature = bwareaopen(vasculature,500); % Larger regions
    vasculature = imclose(vasculature,strel(5));
    curFig = figure(); % Mean Volume intensity
    set(0,'CurrentFigure',curFig)
    clf
    imshow(squeeze(max(vasculature,[],3)))
    title('Baseline Vasculature volume');
    combineFileName = strcat(selpathOutOrig,'\VasculatureMap.fig');
    saveas(curFig,combineFileName);  
    combineFileName = strcat(selpathOutOrig,'\VasculatureMap.tiff');
    saveas(curFig,combineFileName);  
    close(curFig)
    combineFileName = strcat(selpathOutOrig,'\VasculatureMap.mat');
    save(combineFileName,'vasculature');
    
    
    
    
    
    
    
    
    





    
    
    
    %% Perform initial pass to define ROIs of neurons
    fprintf('Generating ROI volume master... ');
    strCR = -1;
    pause(0.1);
    ROIsumBlock = zeros(size(frameBlockRaw));
    for i = 2:2:size(matFileListing)
        % Display processing status
        pVal = round(i/numel(matFileListing)*100);     
        percentageOut = [num2str(pVal) '%%'];
        strOut = [percentageOut repmat(' ',1,10-length(percentageOut)-1)];
        fprintf([strCR strOut]);
        strCR = repmat('\b',1,length(strOut)-1);
        
        
        % Load Stimulation scan volume
        frameBlockRaw = load(matFileListing{i});
        frameBlockRaw = frameBlockRaw.frameBlockS;

        % Generate paired base scan from neighboring base scans
        neigRange = 3;
        % determine start point by limits of the available
        % baseline scans
        midInd = i-1;
        startInd = midInd - (neigRange-1); % canceled out the 2x and /2 due to index stepping and half range down
        if(startInd<1)
            % index is near start, set start to 1
            startInd=1;
        end
        if(startInd>(numel(matFileListing) - (2*neigRange)))
            % index is near start, set start to 1
            startInd=numel(matFileListing) - ((2*neigRange)+1);
        end
        frameBlockBase = zeros(size(frameBlockRaw));
        for neig = 1:neigRange
            offset = startInd + 2*(neig-1);
            frameBlockBaseCur = load(matFileListing{offset});
            frameBlockBaseCur = frameBlockBaseCur.frameBlockB;
            frameBlockBase = frameBlockBase + frameBlockBaseCur;
        end
        frameBlockBase = frameBlockBase./neigRange;
    
        baseMeanPure(i/2) = mean(frameBlockBase,'all'); % measure baseline trends
        baseSTDPure(i/2) = std(frameBlockBase,0,'all');
    
    
        curSubtract = frameBlockRaw - frameBlockBase;

        thresholdedVolume = false(size(frameBlockRaw));
        thresholdedVolume(curSubtract>(3.*stdBaseVolume))=1;

        tempBin = bwareaopen(thresholdedVolume,50);
        ROIsumBlock = ROIsumBlock + tempBin;
    end
    fprintf('\n')        
    
    
    curFig = figure(); % Mean Volume intensity
    set(0,'CurrentFigure',curFig)
    clf
    imshow(squeeze(max(ROIsumBlock/(max(ROIsumBlock,[],'all')),[],3)))
    colormap(jet)
    colorbar
    title(strcat('ROIsumBlock : Max intensity = ',num2str(max(ROIsumBlock,[],'all'))));
    combineFileName = strcat(selpathOutOrig,'\ROIsumBlock.fig');
    saveas(curFig,combineFileName);  
    combineFileName = strcat(selpathOutOrig,'\ROIsumBlock.tiff');
    saveas(curFig,combineFileName);  
    close(curFig)
    combineFileName = strcat(selpathOutOrig,'\ROIsumBlock.mat');
    save(combineFileName,'ROIsumBlock');
    
    

    %% Segment Neuron ROIs from stimulation averages
    MIJ.createImage(uint8(ROIsumBlock)) % Pass variable to ImageJ
    MIJ.run("Subtract Background...", "rolling=2 sliding disable stack"); % minimize the background
    
    %Segment ROIs (parameters may be effected by the number of replicants
    %run or how prolifent the activated regions are across trials.... may
    %have to revisit this calculation.
    MIJ.run("3D Iterative Thresholding", "min_vol_pix=50 max_vol_pix=10000 min_threshold=3 min_contrast=2 criteria_method=EDGES threshold_method=STEP segment_results=All value_method=2");

    % Transfer data from ImageJ back into Matlab
    mijOut = MIJ.getListImages();
    chkStr = '';
    for k = 1:numel(mijOut)
        chkStr = strcat(chkStr,mijOut(k).toCharArray');
    end
    segROIs = false(size(frameBlockBase));

    if(contains(chkStr,'draw')) % Check if ouput figure was generated
        MIJ.selectWindow("draw");

        MIJ.run('Measure Stack...', 'channels slices frames order=czt(default)');
        maxChan = max(MIJ.getColumn('Ch'));
        MIJ.run("Clear Results");
        if(isempty(maxChan))
            MIJ.selectWindow("draw");
            curChanVol = MIJ.getCurrentImage();
            segROIs(curChanVol>0) = 1;
        else
            for chanN = 1:maxChan % Process all channels in output segmentation
                MIJ.run('Duplicate...', strcat('duplicate channels=',num2str(chanN)));
                curChanVol = MIJ.getCurrentImage();
                segROIs(curChanVol>0) = 1;
            end
        end

    end
    MIJ.run('Close All'); 
    try
        MIJ.closeAllWindows(); % Close MIJI figures
    catch
    end           
    
  

    
    % Collect all regions after cleaning up ROI 
    segROIs = imopen(segROIs,strel3d(5));
    segROIs2 = bwareaopen(segROIs,50);
    regions = regionprops(bwconncomp(segROIs,18),'PixelList','Centroid');

    
    % Save ROI sum block and the resulting segmented volume for follow up
    % analysis/validation
    curFig = figure();
    zProjected = sum(ROIsumBlock,3);
    imshow(zProjected/max(zProjected,[],'all'))
    title('Summed Regions 3X above baseline')
    combineFileName = strcat(selpathOut,'\ROISource.fig');
    saveas(curFig,combineFileName);  
    combineFileName = strcat(selpathOut,'\ROISource.tif');
    saveas(curFig,combineFileName);
    close(curFig)
    combineFileName = strcat(selpathOut,'\ROISource.mat');
    save(combineFileName,'ROIsumBlock');  
 
    
    curFig = figure();
    zProjected = sum(ROIsumBlock,3);
    imshow(zProjected/max(zProjected,[],'all'))
    title('Z projected ROI segmentation')
    combineFileName = strcat(selpathOut,'\ROISeg.fig');
    saveas(curFig,combineFileName);  
    combineFileName = strcat(selpathOut,'\ROISeg.tif');
    saveas(curFig,combineFileName);
    close(curFig)
    combineFileName = strcat(selpathOut,'\ROISeg.mat');
    save(combineFileName,'segROIs'); 
    
    SourcePathOut = strcat(selpathOut,'\Source');
    mkdir(SourcePathOut)
    

    
    
    if(numel(regions)~=0) % only perform processing on the rest of this dataset if there exists ROIs to process
        %% collect average intensity response for each detected ROI across all stimulation trials    
        fprintf('Performing ROI Intensity Thresholding... ');
        strCR = -1;
        pause(0.1);

        AveRegionSubIntensity = zeros(numel(regions),numel(matFileListing)/2);
        AveRegionIntensity = zeros(numel(regions),numel(matFileListing)/2);
        AveRegionBaseIntensity = zeros(numel(regions),numel(matFileListing)/2);
        AveRegionIntensity2 = zeros(numel(regions),numel(matFileListing)/2);
        rasterDataAllPer = zeros(numel(regions),numel(matFileListing)/2);
        STDRegionSubIntensity = zeros(numel(regions),numel(matFileListing)/2);
        STDRegionIntensity = zeros(numel(regions),numel(matFileListing)/2);
        STDRegionBaseIntensity = zeros(numel(regions),numel(matFileListing)/2);
        rawRegionSubIntensity = zeros(numel(regions),numel(matFileListing)/2);
        STDrawRegionSubIntensity = zeros(numel(regions),numel(matFileListing)/2);
        for i = 2:2:size(matFileListing)

            % Display processing status
            pVal = round(i/numel(matFileListing)*100);     
            percentageOut = [num2str(pVal) '%%'];
            strOut = [percentageOut repmat(' ',1,10-length(percentageOut)-1)];
            fprintf([strCR strOut]);
            strCR = repmat('\b',1,length(strOut)-1);


            % Load Stimulation scan volume
            frameBlockRaw = load(matFileListing{i});
            frameBlockRaw = frameBlockRaw.frameBlockS;

            % Load paired baseliune scan volume
            frameBlockBaseO = load(matFileListing{i-1});
            frameBlockBaseO = frameBlockBaseO.frameBlockB;

            % Generate paired base scan from neighboring base scans
            neigRange = 3;
            % determine start point by limits of the available
            % baseline scans
            midInd = i-1;
            startInd = midInd - (neigRange-1); % canceled out the 2x and /2 due to index stepping and half range down
            if(startInd<1)
                % index is near start, set start to 1
                startInd=1;
            end
            if(startInd>(numel(matFileListing) - (2*neigRange)))
                % index is near start, set start to 1
                startInd=numel(matFileListing) - ((2*neigRange)+1);
            end
            frameBlockBase = zeros(size(frameBlockRaw));
            for neig = 1:neigRange
                offset = startInd + 2*(neig-1);
                frameBlockBaseCur = load(matFileListing{offset});
                frameBlockBaseCur = frameBlockBaseCur.frameBlockB;
                frameBlockBase = frameBlockBase + frameBlockBaseCur;
            end
            frameBlockBase = frameBlockBase./neigRange;

            curSubtract = frameBlockRaw - frameBlockBase;


            thresholdedVolume = false(size(frameBlockRaw));
            thresholdedVolume(curSubtract>(3.*stdBaseVolume))=1;

            tempBin = bwareaopen(thresholdedVolume,50);
    %         ROIsumBlock = ROIsumBlock + tempBin;


            % for each subtraction, measure the average intesnity of each ROI
            for k = 1:numel(regions)
                curPixels = regions(k).PixelList;
                curRegionVals = size(curPixels,1);
                curRegionCount = size(curPixels,1);
                curRegionIntensity = size(curPixels,1);
                curBaseIntensity = size(curPixels,1);
                for m = 1:size(curPixels,1)
                    curRegionCount(m) = tempBin(curPixels(m,2),curPixels(m,1),curPixels(m,3));
                    curRegionVals(m) = curSubtract(curPixels(m,2),curPixels(m,1),curPixels(m,3));
                    curRegionIntensity(m) = frameBlockRaw(curPixels(m,2),curPixels(m,1),curPixels(m,3));
                    curBaseIntensity(m) = frameBlockBaseO(curPixels(m,2),curPixels(m,1),curPixels(m,3));
                end   

                rasterDataAllPer(k,i/2) = mean(curRegionCount);
                AveRegionSubIntensity(k,i/2) = mean(curRegionVals);
                AveRegionIntensity(k,i/2) = mean(curRegionIntensity);
                AveRegionBaseIntensity(k,i/2) = mean(curBaseIntensity);

                rawRegionSubIntensity(k,i/2) = mean(curRegionIntensity-curBaseIntensity);
                STDrawRegionSubIntensity(k,i/2) = std(curRegionIntensity-curBaseIntensity);

                STDRegionSubIntensity(k,i/2) = std(curRegionVals);
                STDRegionIntensity(k,i/2) = std(curRegionIntensity);
                STDRegionBaseIntensity(k,i/2) = std(curBaseIntensity);
            end
        end
        fprintf('\n')        


        rasterDataAll = false(size(rasterDataAllPer));
        rasterDataAll(rasterDataAllPer>0.5)=1; % consider active if 50% of voxels over 3x STD threshold

    %     for i = size(rasterDataAll,2)
    %         temp = AveRegionIntensity(:,i);
    %         rasterDataAll(temp>regionThresh,i)=1;
    %     end

    %     %% Calcualte the threhsolds of activation for each ROI 
    %     % based on the baseline scan intensity standard deviation (3x)
    %     regionThresh = zeros(numel(regions),1);
    %     for k = 1:numel(regions)
    %         curPixels = regions(k).PixelList;
    %         curRegionVals = size(curPixels,1);
    %         for m = 1:size(curPixels,1)
    %             curRegionVals(m) = stdBaseVolume(curPixels(m,2),curPixels(m,1),curPixels(m,3));
    %         end   
    %         regionThresh(k) = mean(curRegionVals);
    %     end
    %     
    %         regionThresh2 = zeros(numel(regions),1);
    %     for k = 1:numel(regions)
    %         curPixels = regions(k).PixelList;
    %         curRegionVals = size(curPixels,1);
    %         for m = 1:size(curPixels,1)
    %             curRegionVals(m) = stdBaseVolume(curPixels(m,1),curPixels(m,2),curPixels(m,3));
    %         end   
    %         regionThresh2(k) = mean(curRegionVals);
    %     end

    %     %% Determine which ROIs are active during each stimulation scan session
    %     rasterDataAll = false(size(AveRegionIntensity));
    %     for i = size(rasterDataAll,2)
    %         temp = AveRegionIntensity(:,i);
    %         rasterDataAll(temp>regionThresh,i)=1;
    %     end
    %     
    %     rasterDataAll2 = false(size(AveRegionIntensity2));
    %     for i = size(rasterDataAll,2)
    %         temp = AveRegionIntensity2(:,i);
    %         rasterDataAll2(temp>regionThresh2,i)=1;
    %     end

    selpathOutIndROIs = strcat(selpathOut,'\IndividualROIs');
    mkdir(selpathOutIndROIs)
    

        %% Process each stimulation channel/current combination 
        % Determine the ROIs that are consistently which ROIs are consistently
        % above intensity threshold
        rasterData = false(numel(regions),numel(chNames)*numel(currNames));
        c=0;
        for chN = 1:numel(chNames)
            for cuN = 1: numel(currNames)
                c = c+1;
                curCh = chNames{chN};
                curCurr = currNames(cuN);

                %Specify valid inds for current combination
                chRasterMask = strcmp(scanChannel,curCh);
                currRasterMask = false(size(scanCurrent));
                currRasterMask(scanCurrent==curCurr)=1;
                currRasterMask(chRasterMask==0)=0;
                validInds = find(currRasterMask); 

                %pull data
                curRaster = rasterDataAll(:,validInds);
                ROIprobability = mean(curRaster,2);

                %Determine which ROIs were consistently active 
                tempRaster = ROIprobability>probThresh;
                rasterData(:,c)= tempRaster;

                % generate probability mask 
                probVol = zeros(size(frameBlockRaw));
                for k = 1:numel(regions)
                    if(ROIprobability(k)>0)
                        curPixels = regions(k).PixelList;
                        for m = 1:size(curPixels,1)
                            probVol(curPixels(m,2),curPixels(m,1),curPixels(m,3)) = ROIprobability(k);
                        end   
                    end
                end
                probabilityMap{chN,cuN} = probVol;  

                
                % generate individual trial activation maps
                for trialNum =1:size(curRaster,2)
                    activeVol = false(size(frameBlockRaw));
                    for k = 1:numel(regions)
                        if(curRaster(k,trialNum)>0)
                            curPixels = regions(k).PixelList;
                            for m = 1:size(curPixels,1)
                                activeVol(curPixels(m,2),curPixels(m,1),curPixels(m,3)) = 1;
                            end   
                        end
                    end
                    curFig = figure();
                    curView = squeeze(max(activeVol ,[],3));
                    imshow(curView)
                    daspect([1,1,1])
                    combineFileName = strcat(selpathOutIndROIs,'\','IndividualActivated_CH-',chNames{chN},'_CURR-',num2str(currNames(cuN)),'_trial',num2str(trialNum),'.fig');
                    saveas(curFig,combineFileName);  
                    combineFileName = strcat(selpathOutIndROIs,'\','IndividualActivated_CH-',chNames{chN},'_CURR-',num2str(currNames(cuN)),'_trial',num2str(trialNum),'.tiff');
                    saveas(curFig,combineFileName);  
                    close(curFig)
                    combineFileName = strcat(selpathOutIndROIs,'\','IndividualActivated_CH-',chNames{chN},'_CURR-',num2str(currNames(cuN)),'_trial',num2str(trialNum),'.tiff');
                    save(combineFileName,'activeVol');
                end
                
                
                % generate activation mask 
                activeVol = false(size(frameBlockRaw));
                for k = 1:numel(regions)
                    if(tempRaster(k)>0)
                        curPixels = regions(k).PixelList;
                        for m = 1:size(curPixels,1)
                            activeVol(curPixels(m,2),curPixels(m,1),curPixels(m,3)) = 1;
                        end   
                    end
                end
                activeMap{chN,cuN} = activeVol;

                
                
                ChanPathOut = strcat(SourcePathOut,'\',chNames{chN});
                mkdir(ChanPathOut)
                ChanCurPathOut = strcat(ChanPathOut,'\',num2str(currNames(cuN)));
                mkdir(ChanCurPathOut)
                baselineOut = strcat(ChanCurPathOut,'\Baseline');
                mkdir(baselineOut)
              StimOut = strcat(ChanCurPathOut,'\Stimulation');
                mkdir(StimOut)
                
                % Save source 
                for k = 1:numel(validInds)
                    % Load baseline scan volume
                    frameBlockBase = load(matFileListing{2*validInds(k)-1});
                    frameBlockBase = frameBlockBase.frameBlockB;
    
                    activityMap = figure();
                    imshow(squeeze(max(frameBlockBase,[],3)))
                    colormap('jet')
                    title(strcat('Source Baseline Images for : ',curCh,' ',num2str(curCurr)));
                    combineFileName = strcat(baselineOut,'\',curCh,'___',num2str(curCurr),'BaselineImages',num2str(k),'.fig');
                    saveas(activityMap,combineFileName);  
                    combineFileName = strcat(baselineOut,'\',curCh,'___',num2str(curCurr),'BaselineImages',num2str(k),'.tiff');
                    saveas(activityMap,combineFileName);  
                    close(activityMap)
                    ccombineFileName = strcat(baselineOut,'\',curCh,'___',num2str(curCurr),'BaselineImages',num2str(k),'.mat');
                    save(combineFileName,'frameBlockBase');

                    
                    % Load stimulation scan volume
                    frameBlockRaw = load(matFileListing{2*validInds(k)});
                    frameBlockRaw = frameBlockRaw.frameBlockS;
    
                    activityMap = figure();
                    imshow(squeeze(max(frameBlockRaw,[],3)))
                    colormap('jet')
                    title(strcat('Source Stimulation Images for : ',curCh,' ',num2str(curCurr)));
                    combineFileName = strcat(StimOut,'\',curCh,'___',num2str(curCurr),'StimulationImages',num2str(k),'.fig');
                    saveas(activityMap,combineFileName);  
                    combineFileName = strcat(StimOut,'\',curCh,'___',num2str(curCurr),'StimulationImages',num2str(k),'.tiff');
                    saveas(activityMap,combineFileName);  
                    close(activityMap)
                    combineFileName = strcat(StimOut,'\',curCh,'___',num2str(curCurr),'StimulationImages',num2str(k),'.mat');
                    save(combineFileName,'frameBlockRaw');
                    
                    
                    
                end
    
                %Plot subtraction of stim and baseline
                subtractionOut = strcat(ChanCurPathOut,'\Subtraction');
                mkdir(subtractionOut)
                for k = 1:numel(subtract_individual)
                    activityMap = figure();
                    imshow(squeeze(max(subtract_individual{k},[],3)))
                    colormap('jet')
                    title(strcat('Subtractions for : ',curCh,'_',num2str(curCurr),'_',num2str(k)));
                    combineFileName = strcat(subtractionOut,'\',curCh,'___',num2str(curCurr),'SubtractionImages',num2str(k),'.fig');
                    saveas(activityMap,combineFileName);  
                    combineFileName = strcat(subtractionOut,'\',curCh,'___',num2str(curCurr),'SubtractionImages',num2str(k),'.tiff');
                    saveas(activityMap,combineFileName);  
                    close(activityMap)
                    combineFileName = strcat(subtractionOut,'\',curCh,'___',num2str(curCurr),'SubtractionImages',num2str(k),'.tiff');
                    save(combineFileName,'subtract_individual{k}');

                end
    
    
                %Plot segmentation of stim and baseline
                subtractionOut = strcat(ChanCurPathOut,'\Segmentation');
                mkdir(subtractionOut)
                for k = 1:numel(segment_individual)
                    activityMap = figure();
                    imshow(squeeze(max(segment_individual{k},[],3)))
                    title(strcat('Segmentations for : ',curCh,'_',num2str(curCurr),'_',num2str(k)));
                    combineFileName = strcat(subtractionOut,'\',curCh,'___',num2str(curCurr),'SegmentationImages',num2str(k),'.fig');
                    saveas(activityMap,combineFileName);  
                    combineFileName = strcat(subtractionOut,'\',curCh,'___',num2str(curCurr),'SegmentationImages',num2str(k),'.tiff');
                    saveas(activityMap,combineFileName);  
                    close(activityMap)
    
                    activityMap = figure();
                    temp = double(squeeze(max(segment_individual{k},[],3)));
                    temp(temp==1)=0.9;
                    imshow(temp)
                    colormap('jet')
                    title(strcat('Segmentations for : ',curCh,'_',num2str(curCurr),'_',num2str(k)));
                    combineFileName = strcat(subtractionOut,'\',curCh,'___',num2str(curCurr),'RedSegmentationImages',num2str(k),'.fig');
                    saveas(activityMap,combineFileName);  
                    combineFileName = strcat(subtractionOut,'\',curCh,'___',num2str(curCurr),'RedSegmentationImages',num2str(k),'.tiff');
                    saveas(activityMap,combineFileName);  
                    close(activityMap)
                    combineFileName = strcat(subtractionOut,'\',curCh,'___',num2str(curCurr),'RedSegmentationImages',num2str(k),'.mat');
                    save(combineFileName,'segment_individual{k}');
                end
                
                % Plot 2D intensity heatmap plots (maxed, not summed)
                curView = squeeze(max(IntensityMap{chN,cuN},[],3));
                if(~isempty(curView))
                    curFig = figure();
                    set(0,'CurrentFigure',curFig)
                    clf
                    imshow(curView)
                    colormap(jet)
                    title(strcat('Average Activated Subtracted Image :  ',chNames{chN},'  :  ',num2str(currNames(cuN)),' uA'))
                    colorbar
                    combineFileName = strcat(ChanPathOUTPUT,'\',chNames{chN},'___',num2str(currNames(cuN)),'_MeanActiveSubract','.fig');
                    saveas(curFig,combineFileName);  
                    combineFileName = strcat(ChanPathOUTPUT,'\',chNames{chN},'___',num2str(currNames(cuN)),'_MeanActiveSubract','.tiff');
                    saveas(curFig,combineFileName);
                    close(curFig)
                    combineFileName = strcat(ChanPathOUTPUT,'\',chNames{chN},'___',num2str(currNames(cuN)),'_MeanActiveSubract','.mat');
                    save(combineFileName,'IntensityMap{chN,cuN}');
                end
                
                
            end
        end






       %% Plot average basescan intensity
        curFig = figure();
        set(0,'CurrentFigure',curFig)
        clf
        errorbar(baseMeanPure,baseSTDPure,'x')
        xlabel('Scan')
        ylabel('Average Intensity')
        combineFileName = strcat(selpathOut,'\BaseIntensity.fig');
        saveas(curFig,combineFileName);  
        combineFileName = strcat(selpathOut,'\BaseIntensity.tiff');
        saveas(curFig,combineFileName);
        close(curFig)









        %% Generate  plot of neurons acticated over course of experiment
        curFig = figure();
        set(0,'CurrentFigure',curFig)
        clf
        hold on
        title('Individual neuron regions - Activated')
        curRegion = false(size(frameBlockRaw));

        for i = 1:numel(regions)% highlight individual regions
            curPixels = regions(i).PixelList;
            for m = 1:size(curPixels,1)
                curRegion(curPixels(m,2),curPixels(m,1),curPixels(m,3))=1;
            end    
        end

        curView = squeeze(max(curRegion ,[],3));
        imshow(curView)
        daspect([1,1,1])
        axis tight
        camlight 
        lighting gouraud
        combineFileName = strcat(selpathOut,'\','RegionsActivated.fig');
        saveas(curFig,combineFileName);  
        combineFileName = strcat(selpathOut,'\','RegionsActivated.tiff');
        saveas(curFig,combineFileName);  
        close(curFig)


        % save segmentation volumes
        SegPathOut = strcat(selpathOutOrig,'\SegVolumes');
        mkdir(SegPathOut)
        for k = 1:numel(chNames)
            for j = 1:numel(currNames) 
                combineFileName = strcat(SegPathOut,'\','SegVolume_CH-',chNames{k},'_CURR-',num2str(currNames(j)),'.mat');
                segVol = activeMap{k,j};
                save(combineFileName,'segVol');  
            end
        end

        %% Plot probability masks for different channels and currents
        for k = 1:numel(chNames)
            % Generate figures for each channel then iterate through currents
            activityMap = figure();
            for j = 1:numel(currNames) 
                % Display activated regions
                if(~isempty(probabilityMap{k,j}))
                    subplot(numel(currNames),2,((j-1)*2)+1)
                    imshow(squeeze(max(probabilityMap{k,j},[],3)))
                    colormap('parula')
                    ax = gca;
                    set(ax,'XTick',[],'YTick',[]);
                    ylabel(strcat(num2str(currNames(j)), ' uA')) 
                    if(j==1)
                        title('Probability')
                    end

                    subplot(numel(currNames),2,((j-1)*2)+2)
                    imshow(squeeze(max(double(100*activeMap{k,j}),[],3)))
                    colormap('parula')
                    ax = gca;
                    set(ax,'XTick',[],'YTick',[]);
                    if(j==1)
                        title('Segmented')
                    end
                end
            end

            % Save activation summary plots
            figure(activityMap)
            activityMap.WindowState = 'maximized';
            sgtitle(strcat('MaxProjSubtracted CH: ',chNames{k}))
            combineFileName = strcat(selpathOut,'\','MaxProjSubtracted_CH-',chNames{k},'.fig');
            saveas(activityMap,combineFileName);  
            combineFileName = strcat(selpathOut,'\','MaxProjSubtracted_CH-',chNames{k},'.tiff');
            saveas(activityMap,combineFileName);  
            close(activityMap)
        end



        %% Process the raster data to find only regions with monotonically
        % increasing activity (pretty strightforward with 0 ands 1 only)
        monoRasterData = rasterData;
        for r = 1:size(rasterData,1)
            for cCh = 1:numel(chNames)
                mono = 1;
                prior = 0;
                for cCu = 1:numel(currNames)
                    ind =((cCh - 1)* numel(currNames)) + cCu;
                    if(rasterData(r,ind)<prior) % Check if current less than last
                        mono = 0; % monotonically increasing check failed
                    end
                    prior = rasterData(r,ind);
                end
                if(mono==0)
                    % if region was not monotonically increasing for
                    % channel/current combo remove it
                    indA =((cCh - 1)* numel(currNames)) + 1;
                    indB =((cCh - 1)* numel(currNames)) + cCu;
                    monoRasterData(r,indA:indB)=0;
                end
            end
        end








        %% Perform analysis on all regions and determine the intensity history for each ROI 
        % Process to both represent results as a function of
        % channel/current combo and over time.
        regionBaseValMeans = cell(size(rasterData));
        regionStimValMeans = cell(size(rasterData));
        regionActiveValMeans = cell(size(rasterData));
        regionBaseValSTDs = cell(size(rasterData));
        regionStimValSTDs = cell(size(rasterData));
        regionActiveValSTDs = cell(size(rasterData));
        differenceMeans = cell(size(rasterData));
        differenceSTDs = cell(size(rasterData));
        differenceActiveMeans = cell(size(rasterData));
        differenceActiveSTDs = cell(size(rasterData));

        selpathOutIntensity = strcat(selpathOut,'\IntensityAnalysis');
        mkdir(selpathOutIntensity)

        if(size(rasterData,1)>0) % process if there are detected regions
            c=0;
            fprintf('Performing Intensity Analysis... ');
            strCR = -1;
            pause(0.1);
            for chN = 1:numel(chNames) % Iterate through all raster groups.... again
                for cuN = 1: numel(currNames)
                    c = c+1;
                    ind =((chN - 1)* numel(currNames)) + cuN;

                    curCh = chNames{chN};
                    curCurr = currNames(cuN);

                    chRasterMask = strcmp(scanChannel,curCh);
                    currRasterMask = false(size(scanCurrent));
                    currRasterMask(scanCurrent==curCurr)=1;
                    currRasterMask(chRasterMask==0)=0;
                    validInds = find(currRasterMask);


                    % Display processing status
                    pVal = round(c/(numel(chNames)*numel(currNames))*100);     
                    percentageOut = [num2str(pVal) '%%'];
                    strOut = [percentageOut repmat(' ',1,10-length(percentageOut)-1)];
                    fprintf([strCR strOut]);
                    strCR = repmat('\b',1,length(strOut)-1);
                    pause(0.1);

                    baseSum = zeros(size(rasterData,1),size(validInds,1),1);
                    activeSum = zeros(size(rasterData,1),size(validInds,1),1);
                    baseSTD = zeros(size(rasterData,1),size(validInds,1),1);
                    activeSTD = zeros(size(rasterData,1),size(validInds,1),1);
                    differenceMean = zeros(size(rasterData,1),size(validInds,1),1);
                    differenceSTD = zeros(size(rasterData,1),size(validInds,1),1);
                    if(~isempty(validInds)) % Collect region intensities
                        for k = 1:size(validInds,1)       

                            % Load Stimulation scan volume
                            frameBlockRaw = load(matFileListing{2*validInds(k)});
                            frameBlockRaw = frameBlockRaw.frameBlockS;

                            % Load paired baseliune scan volume
                            frameBlockBase = load(matFileListing{2*validInds(k)-1});
                            frameBlockBase = frameBlockBase.frameBlockB;

                            % Run each region
                            for m = 1:size(rasterData,1)
                                % Load masks for each target ROI
                                curPixels = regions(m).PixelList;
                                curBase = zeros(size(curPixels,1),1);
                                curStim =  zeros(size(curPixels,1),1);
                                for px = 1:size(curPixels,1)
                                    curBase(px) = frameBlockBase(curPixels(px,2),curPixels(px,1),curPixels(px,3));
                                    curStim(px) = frameBlockRaw(curPixels(px,2),curPixels(px,1),curPixels(px,3));
                                end    


                                baseSum(m,k) = mean(curBase);
                                activeSum(m,k) = mean(curStim);
                                baseSTD(m,k) = std(curBase);
                                activeSTD(m,k) = std(curStim);

                                differenceMean(m,k) = mean(curStim-curBase);
                                differenceSTD(m,k) = std(curStim-curBase);

                            end
                        end
                        % Save data to raster structure
                        for m = 1:size(rasterData,1)

                            regionBaseValMeans{m,ind} = baseSum(m,:);
                            regionBaseValSTDs{m,ind} = baseSTD(m,:);

                            if(rasterData(m,ind))
                                regionActiveValMeans{m,ind} = activeSum(m,:);
                                regionActiveValSTDs{m,ind} = activeSTD(m,:);

                                differenceActiveMeans{m,ind} = differenceMean(m,:);
                                differenceActiveSTDs{m,ind} = differenceSTD(m,:);
                            else
                                regionStimValMeans{m,ind} = activeSum(m,:);
                                regionStimValSTDs{m,ind} = activeSTD(m,:);

                                differenceMeans{m,ind} = differenceMean(m,:);
                                differenceSTDs{m,ind} = differenceSTD(m,:);
                            end
                        end
                    end
                end
            end
            fprintf('\n')
        end 


        % Find the 100 regions that have the highest signal to noise ratio 
        if(size(rasterData,1)>100)
            snr = zeros(size(rasterData,1),1);
            for j = 1:size(rasterData,1)
                curBase = regionBaseValMeans(j,:);
                curBase = curBase(~cellfun('isempty',curBase));
                curBase = cell2mat(curBase);
                curActive = regionActiveValMeans(j,:);
                curActive = curActive(~cellfun('isempty',curActive));
                curActive = cell2mat(curActive);
                snr(j) = mean(curActive)/mean(curBase);
            end

            [~,largeInd] = maxk(snr,100);
            else
            largeInd = 1:size(rasterData,1);
        end
    %%
        lrgNum = 0;
        if(size(rasterData,1)>0) % process if there are detected regions
            for j = 1:size(rasterData,1)
                if(ismember(j,largeInd)) % Is roi one of the largest?
                    if(sum(rasterData(j,:))>0)
                        lrgNum = lrgNum +1;

                        curBase = regionBaseValMeans(j,:);
                        curBase = curBase(~cellfun('isempty',curBase));
                        curBase = cell2mat(curBase);
                        curStim = regionStimValMeans(j,:);
                        stimInd = ~cellfun('isempty',curStim);
                        curStim = curStim(~cellfun('isempty',curStim));
                        curStim = cell2mat(curStim);
                        curActive = regionActiveValMeans(j,:);
                        actInd = ~cellfun('isempty',curActive);
                        curActive = curActive(~cellfun('isempty',curActive));
                        curActive = cell2mat(curActive);
                        curDiff = differenceMeans(j,:);
                        curDiff = curDiff(~cellfun('isempty',curDiff));
                        curDiff = cell2mat(curDiff);
                        curDiffA = differenceActiveMeans(j,:);
                        curDiffA = curDiffA(~cellfun('isempty',curDiffA));
                        curDiffA = cell2mat(curDiffA);


                        % Get master for stimulation scenarios
                        masterCell = regionStimValMeans(j,:);
                        for k = 1:numel(masterCell)
                            temp = regionActiveValMeans(j,k);
                            if(~isempty(temp{1}))
                                masterCell(k) = temp;
                            end
                        end
                        cumlaInd = zeros(numel(masterCell),1);
                        cumlaInd(1) = 0;
                        cumlaInd(2) = numel(masterCell{1});
                        for k = 2:numel(masterCell)
                           cumlaInd(k+1) =  numel(masterCell{k}) + cumlaInd(k);
                        end

                        curBaseSTD = regionBaseValSTDs(j,:);
                        curBaseSTD = curBaseSTD(~cellfun('isempty',curBaseSTD));
                        curBaseSTD = cell2mat(curBaseSTD);
                        curStimSTD = regionStimValSTDs(j,:);                
                        curStimSTD = curStimSTD(~cellfun('isempty',curStimSTD));
                        curStimSTD = cell2mat(curStimSTD);
                        curActiveSTD = regionActiveValSTDs(j,:);
                        curActiveSTD = curActiveSTD(~cellfun('isempty',curActiveSTD));
                        curActiveSTD = cell2mat(curActiveSTD);


                        % Plot intensity analysis of Ro1
                        intMax = 1.1*max([curBase,curStim,curActive]) + max([curBaseSTD,curStimSTD,curActiveSTD]);
                        regFig = figure(); % Above is to set y limit 
                        set(0,'CurrentFigure',regFig)
                        clf
                        subplot(1,9,1:3)

                        errorbar(curBase,curBaseSTD,'x')
                        ylim([0 intMax])
                        title('Baseline')

                        % Stretch out inds to match data - this was a pain to
                        % code
                        stimInd2 = zeros(numel(curBase),1);
                        for k = 1:numel(actInd)
                            if(stimInd(k)==1)
                                stimInd2(cumlaInd(k)+1:cumlaInd(k+1))=cumlaInd(k)+1:cumlaInd(k+1);
                            end
                        end
                        stimInd2(stimInd2==0)=[];
                        actInd2 = zeros(numel(curBase),1);
                        for k = 1:numel(actInd)
                            if(actInd(k)==1)
                                actInd2(cumlaInd(k)+1:cumlaInd(k+1))=cumlaInd(k)+1:cumlaInd(k+1);
                            end
                        end
                        actInd2(actInd2==0)=[];

                        % Generate current ROI
                        curPixels = regions(j).PixelList;
                        temp = zeros(size(curPixels,1),1); % Find the average STD for this ROI
                        for px = 1:size(curPixels,1)
                            temp(px) = stdBaseVolume(curPixels(px,2),curPixels(px,1),curPixels(px,3));
                        end
                        roiThresh1 = mean(temp);
                        roiThresh2 = 2*roiThresh1;

                        % continue plotting intensity analysis of Ro1
                        subplot(1,9,4:6)
                        hold on
                        errorbar(stimInd2,curStim,curStimSTD,'x')
                        errorbar(actInd2,curActive,curActiveSTD,'o','Color','red')
                        ylim([0 intMax])
                        title('Stimulated')

                        subplot(1,9,7:9)
                        hold on
                        scatter(stimInd2,curDiff,'x')
                        scatter(actInd2,curDiffA,'o','MarkerEdgeColor','red')

                        ylim([1.1*min([curDiff curDiffA]) 1.1*max([roiThresh2,curDiff curDiffA])])
                        yline(roiThresh1,'--','1x STD');
                        yline(roiThresh2,'--','2x STD');
                        title('Difference')

                        regFig.WindowState = 'maximized';
                        combineFileName = strcat(selpathOutIntensity,'\RegionIntensities',num2str(lrgNum),'.fig');
                        saveas(regFig,combineFileName);  
                        combineFileName = strcat(selpathOutIntensity,'\RegionIntensities',num2str(lrgNum),'.tiff');
                        saveas(regFig,combineFileName); 
                        close(regFig)

                        % Plot individual region
                        curRegion = false(size(frameBlockRaw));
                        curPixels = regions(j).PixelList;
                        for m = 1:size(curPixels,1)
                            curRegion(curPixels(m,2),curPixels(m,1),curPixels(m,3))=1;
                        end    
                        regIntFig = figure();
                        set(0,'CurrentFigure',regIntFig)
                        clf
                        curView = squeeze(max(curRegion,[],3));
                        imshow(curView)

                        combineFileName = strcat(selpathOutIntensity,'\RegionPlot',num2str(lrgNum),'.fig');
                        saveas(regIntFig,combineFileName);  
                        combineFileName = strcat(selpathOutIntensity,'\RegionPlot',num2str(lrgNum),'.tiff');
                        saveas(regIntFig,combineFileName); 
                        close(regIntFig)
                    end
                end
            end
        end


        %% perform analysis of regions grouped
        curBase = regionBaseValMeans(j,:);
        curBase = curBase(~cellfun('isempty',curBase));
        curBase = cell2mat(curBase);
        % Process each channel current combination
        maxROIs = 40;
        c=0;
        grpSz = size(curBase,2)/(numel(chNames)*numel(currNames));
        for chN = 1:numel(chNames)
            for cuN = 1: numel(currNames)
                c = c+1;
                curCh = chNames{chN};
                curCurr = currNames(cuN);

    %             %Specify valid inds for current combination, which ROI
    %             %traces to show
                chRasterMask = strcmp(scanChannel,curCh);
                currRasterMask = false(size(scanCurrent));
                currRasterMask(scanCurrent==curCurr)=1;
                currRasterMask(chRasterMask==0)=0;
                validInds = find(currRasterMask); % This specifices the trials to highlight (temporal)



                stimFig = figure();
                subFig = figure();
                stimFigAndBase = figure();
                curRegion = false(size(frameBlockRaw));

                % Find valid active regions for channel/current
                curRaster = rasterData(:,c);
                roiInds = find(curRaster);

                if(~isempty(roiInds))
                    allActive = regionActiveValMeans(roiInds,:);
                    
                    %carfully remove empty cells
                    allActive = allActive(~cellfun('isempty',allActive));
        
                    % pad any active maps with 0s if we come across them...
                    maxSz = 0;
                    for k = 1:numel(allActive)
                        maxSz = max([maxSz, numel(allActive{k})]);
                    end
                    for k = 1:numel(allActive)
                        if(numel(allActive{k})<maxSz)
                            %pad array
                            newarr = zeros(maxSz,1);
                            newarr(1:numel(allActive{k})) = allActive{k};
                            allActive{k} = newarr';
                        end
                    end
                    allActive = cell2mat(allActive);
                    stimMod = 0.3*max(allActive,[],'all');

                    numToProc = min([numel(roiInds) maxROIs]);
                    
                    
                    % divide the active regions into thirds and around 15
                    % for each group if possible, if not, then take what we
                    % can get
                    grp1 = [];
                    grp2 = [];
                    grp3 = [];
                    HM = size(frameBlockRaw,3); % height of image stack
                    for v = 1:numel(roiInds)
                        curCent = regions(roiInds(v)).Centroid;
                        if(curCent(3)<HM/3)
                            %lower 3rd
                            if(numel(grp1)<16)
                                grp1 = [grp1,roiInds(v)];
                            end
                        elseif(curCent(3)<(2*HM/3))
                            % mid third
                            if(numel(grp2)<16)
                                grp2 = [grp2,roiInds(v)];
                            end
                        else
                            % upper third
                            if(numel(grp3)<16)
                                grp3 = [grp3,roiInds(v)];
                            end
                        end
                    end
                    roiInds2 = [grp1,grp2,grp3]; % new group of ROIs to display
                    for v = 1:numel(roiInds2)% % Iterate through ROIs active for current scenatio
                        j = roiInds2(v); % pull index

                        curBase = regionBaseValMeans(j,:);
                        curBase = curBase(~cellfun('isempty',curBase));
                        curBase = cell2mat(curBase);
                        curStim = regionStimValMeans(j,:);
                        stimInd = ~cellfun('isempty',curStim);
                        curStim = curStim(~cellfun('isempty',curStim));
                        curStim = cell2mat(curStim);
                        curActive = regionActiveValMeans(j,:);
                        actInd = ~cellfun('isempty',curActive);
                        curActive = curActive(~cellfun('isempty',curActive));
                        curActive = cell2mat(curActive);
                        curDiff = differenceMeans(j,:);
                        curDiff = curDiff(~cellfun('isempty',curDiff));
                        curDiff = cell2mat(curDiff);
                        curDiffA = differenceActiveMeans(j,:);
                        curDiffA = curDiffA(~cellfun('isempty',curDiffA));
                        curDiffA = cell2mat(curDiffA);


                        % Get master for stimulation scenarios
                        masterCell = regionStimValMeans(j,:);
                        for k = 1:numel(masterCell)
                            temp = regionActiveValMeans(j,k);
                            if(~isempty(temp{1}))
                                masterCell(k) = temp;
                            end
                        end
                        cumlaInd = zeros(numel(masterCell),1);
                        cumlaInd(1) = 0;
                        cumlaInd(2) = numel(masterCell{1});
                        for k = 2:numel(masterCell)
                           cumlaInd(k+1) =  numel(masterCell{k}) + cumlaInd(k);
                        end

                        curBaseSTD = regionBaseValSTDs(j,:);
                        curBaseSTD = curBaseSTD(~cellfun('isempty',curBaseSTD));
                        curBaseSTD = cell2mat(curBaseSTD);
                        curStimSTD = regionStimValSTDs(j,:);                
                        curStimSTD = curStimSTD(~cellfun('isempty',curStimSTD));
                        curStimSTD = cell2mat(curStimSTD);
                        curActiveSTD = regionActiveValSTDs(j,:);
                        curActiveSTD = curActiveSTD(~cellfun('isempty',curActiveSTD));
                        curActiveSTD = cell2mat(curActiveSTD);


                        % Stretch out inds to match data
                        stimInd2 = zeros(numel(curBase),1);
                        for k = 1:numel(actInd)
                            if(stimInd(k)==1)
                                stimInd2(cumlaInd(k)+1:cumlaInd(k+1))=cumlaInd(k)+1:cumlaInd(k+1);
                            end
                        end
                        stimInd2(stimInd2==0)=[];
                        actInd2 = zeros(numel(curBase),1);
                        for k = 1:numel(actInd)
                            if(actInd(k)==1)
                                actInd2(cumlaInd(k)+1:cumlaInd(k+1))=cumlaInd(k)+1:cumlaInd(k+1);
                            end
                        end
                        actInd2(actInd2==0)=[];

                        % Generate current ROI
                        curPixels = regions(j).PixelList;
                        temp = zeros(size(curPixels,1),1); % Find the average STD for this ROI
                        for px = 1:size(curPixels,1)
                            temp(px) = stdBaseVolume(curPixels(px,2),curPixels(px,1),curPixels(px,3));
                        end




                        % continue plotting intensity analysis of Ro1
                        figure(stimFig)
                        hold on
                        stimAll = zeros(numel(curStim)+numel(curActive),1);
                        stimAll(stimInd2)=curStim;
                        stimAll(actInd2)=curActive;
                        stimAll = stimAll+stimMod*(v-1);
                        plot(1:numel(stimAll),stimAll,'Color','black')


                        figure(stimFigAndBase)
                        hold on
                        stimBaseAll = zeros(2*numel(stimAll),1);
                        for k = 1:(numel(chNames)*numel(currNames))
                            stimSt = 2*grpSz*(k-1) + 1;
                            stimEn = 2*grpSz*(k-1) + grpSz;
                            baseSt = 2*grpSz*(k-1) + grpSz+1;
                            baseEn = 2*grpSz*(k-1) + 2*grpSz;
                            curBaseAdd = curBase+stimMod*(v-1);
                            stimBaseAll(stimSt:stimEn) = stimAll(grpSz*(k-1)+1:grpSz*(k));
                            stimBaseAll(baseSt:baseEn) = curBaseAdd(grpSz*(k-1)+1:grpSz*(k));
                        end
                        plot(1:numel(stimBaseAll),stimBaseAll,'Color','black')


                        figure(subFig) % add another ROI trace to 
                        hold on
                        subAll = zeros(numel(curDiff)+numel(curDiffA),1);
                        subAll(stimInd2)=curDiff;
                        subAll(actInd2)=curDiffA;
                        subAll = subAll+stimMod*(v-1);
                        plot(1:numel(subAll),subAll,'Color','black')


                        % Add individual region to volume for current
                        % channel/current block
                        curPixels = regions(j).PixelList;
                        for m = 1:size(curPixels,1)
                            curRegion(curPixels(m,2),curPixels(m,1),curPixels(m,3))=1;
                        end    

                    end


                    figMax = numel(roiInds2)*stimMod;

                    % Plot intensity of Channel/Current activated ROIs
                    figure(stimFig)
                    title('Stimulated Intensity')
                    ylabel('ROI Traces')
                    xlabel('Trials')
                    intense=0.2;
                    sInd = grpSz*(c-1);
                    eInd = grpSz*c;
                    patch([sInd eInd eInd sInd], [0 0 figMax figMax], 'red', 'FaceAlpha', intense, 'EdgeColor', 'none')
                    if(curCurr==10)% if max current highlight lower currents as well
                        for j = 1:3
                            intense2=0.2-(j*0.05);
                            sInd = grpSz*((c-j)-1);
                            eInd = grpSz*(c-j);
                            patch([sInd eInd eInd sInd], [0 0 figMax figMax], 'red', 'FaceAlpha', intense2, 'EdgeColor', 'none')
                        end
                    end
                    ylim([-stimMod stimMod*(numel(roiInds2)+1)])
                    set(gca,'YTick',[])
                    combineFileName = strcat(selpathOutIntensity,'\RegionIntensities_Ch',curCh, '_Curr',num2str(curCurr),'.fig');
                    saveas(stimFig,combineFileName);  
                    combineFileName = strcat(selpathOutIntensity,'\RegionIntensities_Ch',curCh, '_Curr',num2str(curCurr),'.tiff');
                    saveas(stimFig,combineFileName); 
                    close(stimFig)


                    % Plot intensity of Channel/Current activated ROIs inluding
                    % baselines interspersed between stim trial blocks
                    figure(stimFigAndBase)
                    title('Stimulated Intensity + Baseline')
                    ylabel('ROI Traces')
                    xlabel('Trials')
                    sInd = 2*grpSz*(c-1);
                    eInd = 2*grpSz*(c-1) + grpSz;
                    patch([sInd eInd eInd sInd], [0 0 figMax figMax], 'red', 'FaceAlpha', intense, 'EdgeColor', 'none')
                    if(curCurr==10)% if max current highlight lower currents as well
                        for j = 1:3
                            intense2=0.2-(j*0.05);
                            sInd = 2*grpSz*((c-j)-1);
                            eInd = 2*grpSz*((c-j)-1) + grpSz;
                            patch([sInd eInd eInd sInd], [0 0 figMax figMax], 'red', 'FaceAlpha', intense2, 'EdgeColor', 'none')
                        end
                    end
                    ylim([-stimMod stimMod*(numel(roiInds2)+1)])
                    set(gca,'YTick',[])

                    combineFileName = strcat(selpathOutIntensity,'\RegionIntensitiesBase_Ch',curCh, '_Curr',num2str(curCurr),'.fig');
                    saveas(stimFigAndBase,combineFileName);  
                    combineFileName = strcat(selpathOutIntensity,'\RegionIntensitiesBase_Ch',curCh, '_Curr',num2str(curCurr),'.tiff');
                    saveas(stimFigAndBase,combineFileName); 
                    close(stimFigAndBase)


                    % Plot difference between stim and baseline
                    figure(subFig)
                    title('Activation Intensity Difference')
                    ylabel('ROI Traces')
                    xlabel('Trials')
                    sInd = grpSz*(c-1);
                    eInd = grpSz*c;
                    patch([sInd eInd eInd sInd], [0 0 figMax figMax], 'red', 'FaceAlpha', intense, 'EdgeColor', 'none')
                    if(curCurr==10)% if max current highlight lower currents as well
                        for j = 1:3
                            intense2=0.2-(j*0.05);
                            sInd = grpSz*((c-j)-1);
                            eInd = grpSz*(c-j);
                            patch([sInd eInd eInd sInd], [0 0 figMax figMax], 'red', 'FaceAlpha', intense2, 'EdgeColor', 'none')
                        end
                    end
                    ylim([-stimMod stimMod*(numel(roiInds2)+1)])
                    set(gca,'YTick',[])

                    combineFileName = strcat(selpathOutIntensity,'\RegionDifferences_Ch',curCh, '_Curr',num2str(curCurr),'.fig');
                    saveas(subFig,combineFileName);  
                    combineFileName = strcat(selpathOutIntensity,'\RegionDifferenecs_Ch',curCh, '_Curr',num2str(curCurr),'.tiff');
                    saveas(subFig,combineFileName); 
                    close(subFig)



                    % Plot ROIs activated
                    if(sum(curRegion,'all')>2)
                        regIntFig = figure();
                        set(0,'CurrentFigure',regIntFig)
                        clf
        %                 curView = squeeze(max(curRegion,[],3));
        %                 imshow(curView)
                        sSurfaces = isosurface(curRegion,0.5);
                        xyz = sSurfaces.vertices;
                        tri = sSurfaces.faces;
                        outXYZ = SurfaceSmooth(xyz, tri, 2, 0.02, 10, 2, 0);
                        sSurfaces.vertices = outXYZ;
                        p2 = patch(sSurfaces,'FaceColor','red','EdgeColor','none','FaceAlpha',1);
                        isonormals(curRegion,p2);  
                        daspect([1,1,1]); view(3); axis tight; camlight; 
                        lighting gouraud % Format view

                        combineFileName = strcat(selpathOutIntensity,'\RegionPlot_Ch',curCh, '_Curr',num2str(curCurr),'.fig');
                        saveas(regIntFig,combineFileName);  
                        combineFileName = strcat(selpathOutIntensity,'\RegionPlot_Ch',curCh, '_Curr',num2str(curCurr),'.tiff');
                        saveas(regIntFig,combineFileName); 
                        close(regIntFig)
                        
                        combineFileName = strcat(selpathOutIntensity,'\RegionPlot_Ch',curCh, '_Curr',num2str(curCurr),'.mat');
                        save(combineFileName,'curRegion','roiInds2','grp1','grp2','grp3');
                    end
                end
            end
        end
        close all



    %         %% Chronologically order activation 
    %         numROIshown = 15;
    %         c=0;
    %         for chN = 1:numel(chNames)
    %             for cuN = 1: numel(currNames)
    %                 c = c+1;
    %                 curCh = chNames{chN};
    %                 curCurr = currNames(cuN);
    % 
    %                 %Specify valid inds for current combination, which ROI
    %                 %traces to show
    %                 chRasterMask = strcmp(scanChannel,curCh);
    %                 currRasterMask = false(size(scanCurrent));
    %                 currRasterMask(scanCurrent==curCurr)=1;
    %                 currRasterMask(chRasterMask==0)=0;
    %                 validInds = find(currRasterMask); % This specifices the trials to highlight (temporal)
    % 
    %                 
    %                 % Find ROIs that are active for given channel/current combo
    %                 curRaster = rasterData(:,c);
    %                 activeROIs = find(curRaster);
    %                 if(~isempty(activeROIs)) %only process datasets with active ROIs
    %                 
    %                     if(numel(activeROIs)>numROIshown)
    %                         % trim down # of ROI traces shown to meet limit
    %                         shownInds = activeROIs(1:numROIshown);
    %                     elseif(numel(activeROIs)<numROIshown)
    %                         % add in additional ROI trances to meet limit shown
    %                         allInds = 1:numel(regions);
    %                         allInds(activeROIs)=[];
    %                         toAdd = datasample(squeeze(allInds),(numROIshown-numel(activeROIs)),'Replace',false);
    %                         shownInds = [activeROIs',toAdd];
    %                         shownInds = sort(shownInds);
    %                     else
    %                         shownInds = activeROIs;
    %                     end
    % 
    % 
    %                     
    %                     
    %                     
    %                     %Pull and process subtracted ROI intensity data %%%%%%%%%
    %                     curIntensity = AveRegionSubIntensity(shownInds,:);
    % 
    %                     intMax = max(curIntensity,[],'all')*0.5; % set limits for spacing rows
    %                     figMax = numROIshown*intMax;
    %                     % space rows by adding to each consecutive row's values
    %                     dispIntensity = curIntensity;
    %                     for i = 2:numROIshown
    %                         dispIntensity(i,:) = curIntensity(i,:) + (i-1)*intMax;
    %                     end
    % 
    %                     regFig = figure();
    %                     hold on
    %                     for i = 1:numROIshown
    %                         if(ismember(shownInds(i),activeROIs))
    %                             % Active trace, plot in blue
    %                             plot(dispIntensity(i,:),'Color',[0.6000    0.8000    1.0000])
    %                         else
    %                             % Inactive trace, plot in grey
    %                             plot(dispIntensity(i,:),'Color',[0.8000    0.8000    0.8000])
    %                         end
    %                     end
    %                     for i = 1:numel(validInds)
    %                         validInd = validInds(i);
    %                         intense=0.5;
    %                         sInd = validInd - 0.5;
    %                         eInd = validInd + 0.5;
    %                         patch([sInd eInd eInd sInd], [0 0 figMax figMax], 'red', 'FaceAlpha', intense, 'EdgeColor', 'none')
    %                     end
    %                     title('ChronicSubIntensities')
    %                     xlabel('Imaging Trial')
    %                     ylabel('Average ROI Intensity')
    %                     regFig.WindowState = 'maximized';
    %                     combineFileName = strcat(selpathOutIntensity,'\ChronicSubIntensities--Ch',curCh,'--Curr',num2str(curCurr),'.fig');
    %                     saveas(regFig,combineFileName);  
    %                     combineFileName = strcat(selpathOutIntensity,'\ChronicSubIntensities--Ch',curCh,'--Curr',num2str(curCurr),'.tiff');
    %                     saveas(regFig,combineFileName); 
    %                     close(regFig)
    %                     
    %                     
    %                     regFig = figure();
    %                     hold on
    %                     for i = 1:numROIshown
    %                         if(ismember(shownInds(i),activeROIs))
    %                             % Active trace, plot in blue
    %                             plot(curIntensity(i,:),'Color',[0.6000    0.8000    1.0000])
    %                         else
    %                             % Inactive trace, plot in grey
    %                             plot(curIntensity(i,:),'Color',[0.8000    0.8000    0.8000])
    %                         end
    %                     end
    %                     for i = 1:numel(validInds)
    %                         validInd = validInds(i);
    %                         intense=0.5;
    %                         sInd = validInd - 0.5;
    %                         eInd = validInd + 0.5;
    %                         patch([sInd eInd eInd sInd], [0 0 intMax intMax], 'red', 'FaceAlpha', intense, 'EdgeColor', 'none')
    %                     end
    %                     title('ChronicSubIntensities')
    %                     xlabel('Imaging Trial')
    %                     ylabel('Average ROI Intensity')
    %                     regFig.WindowState = 'maximized';
    %                     combineFileName = strcat(selpathOutIntensity,'\ChronicSubIntensitiesStacked--Ch',curCh,'--Curr',num2str(curCurr),'.fig');
    %                     saveas(regFig,combineFileName);  
    %                     combineFileName = strcat(selpathOutIntensity,'\ChronicSubIntensitiesStacked--Ch',curCh,'--Curr',num2str(curCurr),'.tiff');
    %                     saveas(regFig,combineFileName); 
    %                     close(regFig)
    %                     
    %                     
    %                     
    %                     
    %                     
    %                     
    %                     %Pull and process raw ROI intensity data %%%%%%%%%
    %                     curIntensity = AveRegionIntensity(shownInds,:);
    % 
    %                     intMax = max(curIntensity,[],'all')*0.5; % set limits for spacing rows
    %                     figMax = numROIshown*intMax;
    %                     % space rows by adding to each consecutive row's values
    %                     dispIntensity = curIntensity;
    %                     for i = 2:numROIshown
    %                         dispIntensity(i,:) = curIntensity(i,:) + (i-1)*intMax;
    %                     end
    % 
    %                     regFig = figure();
    %                     hold on
    %                     for i = 1:numROIshown
    %                         if(ismember(shownInds(i),activeROIs))
    %                             % Active trace, plot in blue
    %                             plot(dispIntensity(i,:),'Color',[0.6000    0.8000    1.0000])
    %                         else
    %                             % Inactive trace, plot in grey
    %                             plot(dispIntensity(i,:),'Color',[0.8000    0.8000    0.8000])
    %                         end
    %                     end
    %                     for i = 1:numel(validInds)
    %                         validInd = validInds(i);
    %                         intense=0.5;
    %                         sInd = validInd - 0.5;
    %                         eInd = validInd + 0.5;
    %                         patch([sInd eInd eInd sInd], [0 0 figMax figMax], 'red', 'FaceAlpha', intense, 'EdgeColor', 'none')
    %                     end
    %                     title('ChronicIntensities')
    %                     xlabel('Imaging Trial')
    %                     ylabel('Average ROI Intensity')
    %                     regFig.WindowState = 'maximized';
    %                     combineFileName = strcat(selpathOutIntensity,'\ChronicIntensities--Ch',curCh,'--Curr',num2str(curCurr),'.fig');
    %                     saveas(regFig,combineFileName);  
    %                     combineFileName = strcat(selpathOutIntensity,'\ChronicIntensities--Ch',curCh,'--Curr',num2str(curCurr),'.tiff');
    %                     saveas(regFig,combineFileName); 
    %                     close(regFig)
    %                     
    %                     
    %                     regFig = figure();
    %                     hold on
    %                     for i = 1:numROIshown
    %                         if(ismember(shownInds(i),activeROIs))
    %                             % Active trace, plot in blue
    %                             plot(curIntensity(i,:),'Color',[0.6000    0.8000    1.0000])
    %                         else
    %                             % Inactive trace, plot in grey
    %                             plot(curIntensity(i,:),'Color',[0.8000    0.8000    0.8000])
    %                         end
    %                     end
    %                     for i = 1:numel(validInds)
    %                         validInd = validInds(i);
    %                         intense=0.5;
    %                         sInd = validInd - 0.5;
    %                         eInd = validInd + 0.5;
    %                         patch([sInd eInd eInd sInd], [0 0 intMax intMax], 'red', 'FaceAlpha', intense, 'EdgeColor', 'none')
    %                     end
    %                     title('ChronicIntensities')
    %                     xlabel('Imaging Trial')
    %                     ylabel('Average ROI Intensity')
    %                     regFig.WindowState = 'maximized';
    %                     combineFileName = strcat(selpathOutIntensity,'\ChronicIntensitiesStacked--Ch',curCh,'--Curr',num2str(curCurr),'.fig');
    %                     saveas(regFig,combineFileName);  
    %                     combineFileName = strcat(selpathOutIntensity,'\ChronicIntensitiesStacked--Ch',curCh,'--Curr',num2str(curCurr),'.tiff');
    %                     saveas(regFig,combineFileName); 
    %                     close(regFig)
    %                     
    %                     
    %                     
    %                     
    %                     
    %                     
    %                     %Pull and process raw ROI intensity data %%%%%%%%%
    %                     curIntensity = rawRegionSubIntensity(shownInds,:);
    % 
    %                     intMax = max(curIntensity,[],'all')*0.5; % set limits for spacing rows
    %                     figMax = numROIshown*intMax;
    %                     % space rows by adding to each consecutive row's values
    %                     dispIntensity = curIntensity;
    %                     for i = 2:numROIshown
    %                         dispIntensity(i,:) = curIntensity(i,:) + (i-1)*intMax;
    %                     end
    % 
    %                     regFig = figure();
    %                     hold on
    %                     for i = 1:numROIshown
    %                         if(ismember(shownInds(i),activeROIs))
    %                             % Active trace, plot in blue
    %                             plot(dispIntensity(i,:),'Color',[0.6000    0.8000    1.0000])
    %                         else
    %                             % Inactive trace, plot in grey
    %                             plot(dispIntensity(i,:),'Color',[0.8000    0.8000    0.8000])
    %                         end
    %                     end
    %                     for i = 1:numel(validInds)
    %                         validInd = validInds(i);
    %                         intense=0.5;
    %                         sInd = validInd - 0.5;
    %                         eInd = validInd + 0.5;
    %                         patch([sInd eInd eInd sInd], [0 0 figMax figMax], 'red', 'FaceAlpha', intense, 'EdgeColor', 'none')
    %                     end
    %                     title('RawSubChronicIntensities')
    %                     xlabel('Imaging Trial')
    %                     ylabel('Average ROI Intensity')
    %                     regFig.WindowState = 'maximized';
    %                     combineFileName = strcat(selpathOutIntensity,'\RawSubChronicIntensities--Ch',curCh,'--Curr',num2str(curCurr),'.fig');
    %                     saveas(regFig,combineFileName);  
    %                     combineFileName = strcat(selpathOutIntensity,'\RawSubChronicIntensities--Ch',curCh,'--Curr',num2str(curCurr),'.tiff');
    %                     saveas(regFig,combineFileName); 
    %                     close(regFig)
    %                     
    %                                         
    %                     regFig = figure();
    %                     hold on
    %                     for i = 1:numROIshown
    %                         if(ismember(shownInds(i),activeROIs))
    %                             % Active trace, plot in blue
    %                             plot(curIntensity(i,:),'Color',[0.6000    0.8000    1.0000])
    %                         else
    %                             % Inactive trace, plot in grey
    %                             plot(curIntensity(i,:),'Color',[0.8000    0.8000    0.8000])
    %                         end
    %                     end
    %                     for i = 1:numel(validInds)
    %                         validInd = validInds(i);
    %                         intense=0.5;
    %                         sInd = validInd - 0.5;
    %                         eInd = validInd + 0.5;
    %                         patch([sInd eInd eInd sInd], [0 0 intMax intMax], 'red', 'FaceAlpha', intense, 'EdgeColor', 'none')
    %                     end
    %                     title('RawSubChronicIntensities')
    %                     xlabel('Imaging Trial')
    %                     ylabel('Average ROI Intensity')
    %                     regFig.WindowState = 'maximized';
    %                     combineFileName = strcat(selpathOutIntensity,'\RawSubChronicIntensitiesStacked--Ch',curCh,'--Curr',num2str(curCurr),'.fig');
    %                     saveas(regFig,combineFileName);  
    %                     combineFileName = strcat(selpathOutIntensity,'\RawSubChronicIntensitiesStacked--Ch',curCh,'--Curr',num2str(curCurr),'.tiff');
    %                     saveas(regFig,combineFileName); 
    %                     close(regFig)
    %                     
    %                     
    %                     
    %                     
    %                     
    %                     
    %                     %Pull and process  average ROI intensity data %%%%%%%%%
    %                     curIntensity = AveRegionIntensity(shownInds,:);
    %                     baseIntensity = AveRegionBaseIntensity(shownInds,:);
    % 
    %                     %merge stim and baseline
    %                     allIntensity = zeros(size(curIntensity,1),2*size(curIntensity,2));
    %                     for i = 1:size(curIntensity,2)
    %                         allIntensity(:,2*i) = curIntensity(:,i);
    %                         allIntensity(:,(2*i)-1) = baseIntensity(:,i);
    %                     end
    % 
    %                     intMax = max(allIntensity,[],'all')*0.5; % set limits for spacing rows
    %                     figMax = numROIshown*intMax;
    %                     % space rows by adding to each consecutive row's values
    %                     dispIntensity = allIntensity;
    %                     for i = 2:numROIshown
    %                         dispIntensity(i,:) = allIntensity(i,:) + (i-1)*intMax;
    %                     end
    % 
    %                     regFig = figure();
    %                     hold on
    %                     for i = 1:numROIshown
    %                         if(ismember(shownInds(i),activeROIs))
    %                             % Active trace, plot in blue
    %                             plot(dispIntensity(i,:),'Color',[0.6000    0.8000    1.0000])
    %                         else
    %                             % Inactive trace, plot in grey
    %                             plot(dispIntensity(i,:),'Color',[0.8000    0.8000    0.8000])
    %                         end
    %                     end
    %                     for i = 1:numel(validInds)
    %                         validInd = 2*validInds(i);
    %                         intense=0.5;
    %                         sInd = validInd - 0.5;
    %                         eInd = validInd + 0.5;
    %                         patch([sInd eInd eInd sInd], [0 0 figMax figMax], 'red', 'FaceAlpha', intense, 'EdgeColor', 'none')
    %                     end
    %                     title('AllChronicIntensities')
    %                     xlabel('Imaging Trial')
    %                     ylabel('Average ROI Intensity')
    %                     regFig.WindowState = 'maximized';
    %                     combineFileName = strcat(selpathOutIntensity,'\AllChronicIntensities--Ch',curCh,'--Curr',num2str(curCurr),'.fig');
    %                     saveas(regFig,combineFileName);  
    %                     combineFileName = strcat(selpathOutIntensity,'\AllChronicIntensities--Ch',curCh,'--Curr',num2str(curCurr),'.tiff');
    %                     saveas(regFig,combineFileName); 
    %                     close(regFig)
    %                     
    %                     
    %                     regFig = figure();
    %                     hold on
    %                     for i = 1:numROIshown
    %                         if(ismember(shownInds(i),activeROIs))
    %                             % Active trace, plot in blue
    %                             plot(allIntensity(i,:),'Color',[0.6000    0.8000    1.0000])
    %                         else
    %                             % Inactive trace, plot in grey
    %                             plot(allIntensity(i,:),'Color',[0.8000    0.8000    0.8000])
    %                         end
    %                     end
    %                     for i = 1:numel(validInds)
    %                         validInd = 2*validInds(i);
    %                         intense=0.5;
    %                         sInd = validInd - 0.5;
    %                         eInd = validInd + 0.5;
    %                         patch([sInd eInd eInd sInd], [0 0 intMax intMax], 'red', 'FaceAlpha', intense, 'EdgeColor', 'none')
    %                     end
    %                     title('AllChronicIntensities')
    %                     xlabel('Imaging Trial')
    %                     ylabel('Average ROI Intensity')
    %                     regFig.WindowState = 'maximized';
    %                     combineFileName = strcat(selpathOutIntensity,'\AllChronicIntensitiesStacked--Ch',curCh,'--Curr',num2str(curCurr),'.fig');
    %                     saveas(regFig,combineFileName);  
    %                     combineFileName = strcat(selpathOutIntensity,'\AllChronicIntensitiesStacked--Ch',curCh,'--Curr',num2str(curCurr),'.tiff');
    %                     saveas(regFig,combineFileName); 
    %                     close(regFig)
    %                     
    %                     
    %                 end
    %             end
    %         end
    %         
    %       
    %         
    %         
    %         







            %% Plot the radial distane of actived neurons fo each channel/ current and density
            % 100 um bins will be used
            selpathOutDist = strcat(selpathOut,'\StimulationDistances');
            mkdir(selpathOutDist)
            contactsOrdered = [9,13,11,29,13,27,15,25,7,17,5,19,3,21,1,23,10,32,12,30,14,28,16,26,8,18,6,20,4,22,2,24];
            if(runDistance==1) % if the mapping file was located

                OverallDepth = zeros(numel(chNames),1); % Construct variables for overall depth analysis
                OverallmaxDist2D = zeros(numel(chNames),numel(currNames));
                OverallmaxDist3D = zeros(numel(chNames),numel(currNames));
                OverallAveDist2D = zeros(numel(chNames),numel(currNames));
                OverallAveDist3D = zeros(numel(chNames),numel(currNames));
                OverallDist2DErr = zeros(numel(chNames),numel(currNames));
                OverallDist3DErr = zeros(numel(chNames),numel(currNames));

                for i = 1:numel(chNames) % Step in reverse for for loop, allows for plotting of maximum axis
                    maxDists2D = zeros(numel(currNames),1);
                    maxDists3D = zeros(numel(currNames),1);
                    density2DFull = zeros(numel(currNames),10);
                    density3DFull = zeros(numel(currNames),10);

                    curCh = chNames{i};
                    max2D = 0;
                    max3D = 0;

                    % Calculate electrode position
                    curElecPos = [0,0,0];
                    for k = 1:numel(contactsOrdered)
                        if(isequal(num2str(contactsOrdered(k)),curCh))
                            % Current channel matches stimulating channel,
                            % pull position
                            curElecPos = elecPos(k,:);
                        end
                    end
                    OverallDepth(i) = curElecPos(3);

                    for j = numel(currNames):-1:1
                                            % For each detected neuron, calculate the distance from
                        % the electrode that is was, in both 2D and 3D space
                        curSegData = activeMap{i,j};
                        neuronRegs = regionprops(curSegData,'centroid');
                        centroids = cat(1,neuronRegs.Centroid);
                        if(numel(centroids)==0)
                            % If no centroids exist just skip this current level
                            continue
                        end
                        allDist2D = zeros(size(centroids,1),1);
                        allDist3D = zeros(size(centroids,1),1);

                        for k = 1:size(centroids,1) % Perform 2D & 3D distance calculations
                            allDist2D(k) = sqrt(((centroids(k,1)-curElecPos(1))*2.187)^2 + ((centroids(k,2)-curElecPos(2))*2.187)^2);
                            allDist3D(k) = sqrt(((centroids(k,1)-curElecPos(1))*2.187)^2 + ((centroids(k,2)-curElecPos(2))*2.187)^2  + ((centroids(k,3)-curElecPos(3))*vertStep)^2 );
    %                         allDist2D(k) = pdist([centroids(k,1:2);curElecPos(1:2)],'euclidean')*microns_per_pixel; % 2D distance
    %                         allDist3D(k) = pdist([centroids(k,:);curElecPos],'euclidean')*microns_per_pixel; % 3D distance
                        end

                        % Calculate max distance for channel/current
                        maxDists2D(j) = max(allDist2D);
                        maxDists3D(j) = max(allDist3D);
                        OverallmaxDist2D(i,j) = max(allDist2D);
                        OverallmaxDist3D(i,j) = max(allDist3D);
                        OverallAveDist2D(i,j) = mean(allDist2D);
                        OverallAveDist3D(i,j) = mean(allDist3D);
                        OverallDist2DErr(i,j) = std(allDist2D);
                        OverallDist3DErr(i,j) = std(allDist3D);

                        % Calculate density for channel/current
                        density2D = histcounts(allDist2D,0:100:1000);
                        density3D = histcounts(allDist3D,0:100:1000);
                        for k = 1:10
                            % calculate 2D density for area surveyed
                            density2DFull(j,k) = density2D(k)/((pi()*(100*k)^3)-(pi()*(100*(k-1))^3));

                            % calculate 3D density for area surveyed
                            density3DFull(j,k) = density3D(k)/(((4/3)*pi()*(100*k)^2)-((4/3)*pi()*(100*(k-1))^2));
                        end


                        % Plot Histogram of all neurons detected for
                        % channel/current
                        curFig = figure();
                        histogram(allDist2D,1:100:1000);
                        xlabel('Stimulation Distance')
                        ylabel('Number of neurons')
                        title(strcat('2D Stimulation counts-CH',curCh,'-CUR',num2str(currNames(j))))

                        combineFileName = strcat(selpathOutDist,'\StimulationHist2D-CH',curCh,'-CUR',num2str(currNames(j)),'.fig');
                        saveas(curFig,combineFileName);  
                        combineFileName = strcat(selpathOutDist,'\StimulationHist2D-CH',curCh,'-CUR',num2str(currNames(j)),'.tiff');
                        saveas(curFig,combineFileName);  
                        close(curFig)


                        curFig = figure();
                        histogram(allDist3D,1:100:1000);
                        xlabel('Stimulation Distance')
                        ylabel('Number of neurons')
                        title(strcat('3D Stimulation counts-CH',curCh,'-CUR',num2str(currNames(j))))

                        combineFileName = strcat(selpathOutDist,'\StimulationHist3D-CH',curCh,'-CUR',num2str(currNames(j)),'.fig');
                        saveas(curFig,combineFileName);  
                        combineFileName = strcat(selpathOutDist,'\StimulationHist3D-CH',curCh,'-CUR',num2str(currNames(j)),'.tiff');
                        saveas(curFig,combineFileName);  
                        close(curFig)




                        % Plot Histogram distribution of neuron distances                    

                        % Plot 2D individual channel+currenct response
                        bins = 0:100:1000;
                        [counts,bins] = hist(allDist2D,bins);
                        max2D = max([max2D,counts]);
                        curFig = figure();
                        set(0,'CurrentFigure',curFig)
                        clf
                        hold on 
                        barh(bins,counts)
                        xlabel('Count')
                        ylabel('Distance (um)')
                        xlim([0 (max2D+1)])
                        title(strcat('Stimulation Distance Histogram 2D: ',curCh,'_-_-',num2str(currNames(j))))

                        combineFileName = strcat(selpathOutDist,'\StimulationDistances2D',curCh,'--',num2str(currNames(j)),'.fig');
                        saveas(curFig,combineFileName);  
                        combineFileName = strcat(selpathOutDist,'\StimulationDistances2D',curCh,'--',num2str(currNames(j)),'.tiff');
                        saveas(curFig,combineFileName);  
                        close(curFig)




                        % Plot 3D individual channel+currenct response
                        bins = 0:100:1000;
                        [counts,bins] = hist(allDist3D,bins);
                        max3D = max([max3D,counts]);
                        curFig = figure();
                        set(0,'CurrentFigure',curFig)
                        clf
                        hold on 
                        barh(bins,counts)
                        xlabel('Count')
                        ylabel('Distance (um)')
                        xlim([0 (max3D+1)])
                        title(strcat('Stimulation Distance Histogram 3D: ',curCh,'_-_-',num2str(currNames(j))))

                        combineFileName = strcat(selpathOutDist,'\StimulationDistances3D',curCh,'--',num2str(currNames(j)),'.fig');
                        saveas(curFig,combineFileName);  
                        combineFileName = strcat(selpathOutDist,'\StimulationDistances3D',curCh,'--',num2str(currNames(j)),'.tiff');
                        saveas(curFig,combineFileName);  
                        close(curFig)
                    end

                    % Plot data for channel analyzed
                    % Maximum radial bouding radius/sphere
                    curFig = figure();
                    set(0,'CurrentFigure',curFig)
                    clf
                    hold on 
                    plot(currNames,maxDists2D,'-o')
                    xlabel('Currents')
                    xticks(currNames)
                    ylabel('Distance (um)')
                    title(strcat('Max 2D Stimulation Distances: Ch',curCh))

                    combineFileName = strcat(selpathOutDist,'\StimulationDistancesMax2D',curCh,'.fig');
                    saveas(curFig,combineFileName);  
                    combineFileName = strcat(selpathOutDist,'\StimulationDistancesMax2D',curCh,'.tiff');
                    saveas(curFig,combineFileName);  
                    close(curFig)


                    curFig = figure();
                    set(0,'CurrentFigure',curFig)
                    clf
                    hold on 
                    plot(currNames,maxDists3D,'-o')
                    xlabel('Currents')
                    xticks(currNames)
                    ylabel('Distance (um)')
                    title(strcat('Max 3D Stimulation Distances: Ch',curCh))

                    combineFileName = strcat(selpathOutDist,'\StimulationDistancesMax3D',curCh,'.fig');
                    saveas(curFig,combineFileName);  
                    combineFileName = strcat(selpathOutDist,'\StimulationDistancesMax3D',curCh,'.tiff');
                    saveas(curFig,combineFileName);  
                    close(curFig)


                    % Plot density of neurons activated in different regions  (One way of plotting the data)               
                    curFig = figure(); % 2D density
                    set(0,'CurrentFigure',curFig)
                    clf
                    hold on 
                    h = bar(currNames,density2DFull);
                    set(h, {'DisplayName'}, {'0-100','100-200','200-300','300-400','400-500','500-600','600-700','700-800','800-900','900-1000'}')
                    % Legend will show names for each color
                    legend() 
                    xlabel('Currents')
                    xticks(currNames)
                    ylabel('Density (Neurons/um^2)')
                    title(strcat('2D Stimulation Density: Ch',curCh))
                    curFig.WindowState = 'maximized';
                    combineFileName = strcat(selpathOutDist,'\StimulationDensity_2D',curCh,'.fig');
                    saveas(curFig,combineFileName);  
                    combineFileName = strcat(selpathOutDist,'\StimulationDensity_2D',curCh,'.tiff');
                    saveas(curFig,combineFileName);  
                    close(curFig)


                    curFig = figure(); % 3D density
                    set(0,'CurrentFigure',curFig)
                    clf
                    hold on 
                    h = bar(currNames,density3DFull);
                    set(h, {'DisplayName'}, {'0-100','100-200','200-300','300-400','400-500','500-600','600-700','700-800','800-900','900-1000'}')
                    % Legend will show names for each color
                    legend() 

                    xlabel('Currents')
                    xticks(currNames)
                    ylabel('Density (Neurons/um^3)')
                    title(strcat('3D Stimulation Density: Ch',curCh))
                    curFig.WindowState = 'maximized';
                    combineFileName = strcat(selpathOutDist,'\StimulationDensity_3D',curCh,'.fig');
                    saveas(curFig,combineFileName);  
                    combineFileName = strcat(selpathOutDist,'\StimulationDensity_3D',curCh,'.tiff');
                    saveas(curFig,combineFileName);  
                    close(curFig)

                end



                % order data in depth analysis... well by depth
                [OverallDepth,newOrder] = sort(OverallDepth);
                OverallmaxDist2D = OverallmaxDist2D(newOrder,:);
                OverallmaxDist3D = OverallmaxDist3D(newOrder,:);
                OverallAveDist2D = OverallAveDist2D(newOrder,:);
                OverallAveDist3D = OverallAveDist3D(newOrder,:);
                OverallDist2DErr = OverallDist2DErr(newOrder,:);
                OverallDist3DErr = OverallDist3DErr(newOrder,:);


                % Depth vs stimulation spread plots
                currNameStr = cell(numel(currNames),1);
                cmap = distinguishable_colors(numel(currNames)); % Channel Colors
                for curr = 1:numel(currNames)

                    % Plot stimulation 2D spread as a factor of depth and
                    % current, individual plots per current level
                    curFig = figure();
                    errorbar(OverallAveDist2D(:,curr),OverallDepth,OverallDist2DErr(:,curr),'-o','horizontal')
                    set(gca, 'YDir','reverse')
                    xlabel('2D Stimulation Spread (um)')
                    ylabel('Depth (um)')
                    title(strcat('2D Stimulation Spread Vs Depth: Current',num2str(currNames(curr))))
                    combineFileName = strcat(selpathOutDist,'\StimulationSpreadDepth_2D_',num2str(currNames(curr)),'uA.fig');
                    saveas(curFig,combineFileName);  
                    combineFileName = strcat(selpathOutDist,'\StimulationSpreadDepth_2D_',num2str(currNames(curr)),'uA.tiff');
                    saveas(curFig,combineFileName);  
                    close(curFig)


                    % Plot stimulation 3D spread as a factor of depth and
                    % current, individual plots per current level
                    curFig = figure();
                    errorbar(OverallAveDist3D(:,curr),OverallDepth,OverallDist3DErr(:,curr),'-o','horizontal')
                    set(gca, 'YDir','reverse')
                    xlabel('2D Stimulation Spread (um)')
                    ylabel('Depth (um)')
                    title(strcat('3D Stimulation Spread Vs Depth: Current',num2str(currNames(curr))))
                    combineFileName = strcat(selpathOutDist,'\StimulationSpreadDepth_3D_',num2str(currNames(curr)),'uA.fig');
                    saveas(curFig,combineFileName);  
                    combineFileName = strcat(selpathOutDist,'\StimulationSpreadDepth_3D_',num2str(currNames(curr)),'uA.tiff');
                    saveas(curFig,combineFileName);  
                    close(curFig)



                    % Plot max stimulation 2D spread as a factor of depth and
                    % current, individual plots per current level
                    curFig = figure();
                    plot(OverallmaxDist2D(:,curr),OverallDepth,'-o')
                    set(gca, 'YDir','reverse')
                    xlabel('2D Stimulation Spread (um)')
                    ylabel('Depth (um)')
                    title(strcat('Max 2D Stimulation Spread Vs Depth: Current',num2str(currNames(curr))))
                    combineFileName = strcat(selpathOutDist,'\MaxStimulationSpreadDepth_2D_',num2str(currNames(curr)),'uA.fig');
                    saveas(curFig,combineFileName);  
                    combineFileName = strcat(selpathOutDist,'\MaxStimulationSpreadDepth_2D_',num2str(currNames(curr)),'uA.tiff');
                    saveas(curFig,combineFileName);  
                    close(curFig)


                    % Plot max stimulation 3D spread as a factor of depth and
                    % current, individual plots per current level
                    curFig = figure();
                    plot(OverallmaxDist3D(:,curr),OverallDepth,'-o')
                    set(gca, 'YDir','reverse')
                    xlabel('2D Stimulation Spread (um)')
                    ylabel('Depth (um)')
                    title(strcat('Max 3D Stimulation Spread Vs Depth: Current',currNames(curr)))
                    combineFileName = strcat(selpathOutDist,'\MaxStimulationSpreadDepth_3D_',num2str(currNames(curr)),'uA.fig');
                    saveas(curFig,combineFileName);  
                    combineFileName = strcat(selpathOutDist,'\MaxStimulationSpreadDepth_3D_',num2str(currNames(curr)),'uA.tiff');
                    saveas(curFig,combineFileName);  
                    close(curFig)

                    currNameStr{curr} = num2str(currNames(curr));
                end



                % Plot stimulation 2D spread as a factor of depth and currents
                curFig = figure();
                hold on 
                for curr = 1:numel(currNames) % Plot traces for each current level
                    RGB = cmap(curr,:);
                    e = errorbar(OverallAveDist2D(:,curr),OverallDepth,OverallDist2DErr(:,curr),'-o','horizontal');
                    e.Color = RGB;
                end
                set(gca, 'YDir','reverse')
                legend(currNameStr) % Legend will show the color assigned for each current
                xlabel('2D Stimulation Spread (um)')
                ylabel('Depth (um)')
                title('2D Stimulation Spread Vs Depth')
                combineFileName = strcat(selpathOutDist,'\StimulationSpreadDepth_2D_Overall.fig');
                saveas(curFig,combineFileName);  
                combineFileName = strcat(selpathOutDist,'\StimulationSpreadDepth_2D_Overall.tiff');
                saveas(curFig,combineFileName);  
                close(curFig)


                % Plot stimulation 3D spread as a factor of depth and currents
                curFig = figure();
                hold on 
                for curr = 1:numel(currNames) % Plot traces for each current level
                    RGB = cmap(curr,:);
                    e = errorbar(OverallAveDist2D(:,curr),OverallDepth,OverallDist2DErr(:,curr),'-o','horizontal');
                    e.Color = RGB;
                end
                set(gca, 'YDir','reverse')
                legend(currNameStr) % Legend will show the color assigned for each current
                xlabel('3D Stimulation Spread (um)')
                ylabel('Depth (um)')
                title('3D Stimulation Spread Vs Depth')
                combineFileName = strcat(selpathOutDist,'\StimulationSpreadDepth_3D_Overall.fig');
                saveas(curFig,combineFileName);  
                combineFileName = strcat(selpathOutDist,'\StimulationSpreadDepth_3D_Overall.tiff');
                saveas(curFig,combineFileName);  
                close(curFig)


                % Plot maximum stimulation 2D spread as a factor of depth and currents
                curFig = figure();
                hold on 
                for curr = 1:numel(currNames) % Plot traces for each current level
                    RGB = cmap(curr,:);
                    e = plot(OverallmaxDist2D(:,curr),OverallDepth,'-o');
                    e.Color = RGB;
                end
                set(gca, 'YDir','reverse')
                legend(currNameStr) % Legend will show the color assigned for each current
                xlabel('2D Maximum Stimulation Spread (um)')
                ylabel('Depth (um)')
                title('2D Maximum Stimulation Spread Vs Depth')
                combineFileName = strcat(selpathOutDist,'\MaxStimulationSpreadDepth_2D_Overall.fig');
                saveas(curFig,combineFileName);  
                combineFileName = strcat(selpathOutDist,'\MaxStimulationSpreadDepth_2D_Overall.tiff');
                saveas(curFig,combineFileName);  
                close(curFig)


                % Plot maximum stimulation 3D spread as a factor of depth and currents
                curFig = figure();
                hold on 
                for curr = 1:numel(currNames) % Plot traces for each current level
                    RGB = cmap(curr,:);
                    e = plot(OverallmaxDist3D(:,curr),OverallDepth,'-o');
                    e.Color = RGB;
                end
                set(gca, 'YDir','reverse')
                legend(currNameStr) % Legend will show the color assigned for each current
                xlabel('3D Maximum Stimulation Spread (um)')
                ylabel('Depth (um)')
                title('3D Maximum Stimulation Spread Vs Depth')
                combineFileName = strcat(selpathOutDist,'\MaxStimulationSpreadDepth_3D_Overall.fig');
                saveas(curFig,combineFileName);  
                combineFileName = strcat(selpathOutDist,'\MaxStimulationSpreadDepth_3D_Overall.tiff');
                saveas(curFig,combineFileName);  
                close(curFig)


            end





            %% Plot the COM distance and volume of activation for each channel/current [just start off with activation]
            chMaxs2D = zeros(numel(chNames),1);
            chMins2D =  zeros(numel(chNames),1);
            chMeans2D =  zeros(numel(chNames),1);
            chSTD2D =  zeros(numel(chNames),1);        

            chMaxs3D = zeros(numel(chNames),1);
            chMins3D =  zeros(numel(chNames),1);
            chMeans3D =  zeros(numel(chNames),1);
            chSTD3D =  zeros(numel(chNames),1); 


            for i = 1:numel(chNames)
                curCh = chNames{i};


                % calculate distances  for current channel and all currents
                %     channelDistances = zeros(1,1);
                curMaxs = zeros(numel(currNames),1);
                curMins =  zeros(numel(currNames),1);
                curMeans =  zeros(numel(currNames),1);
                curSTD =  zeros(numel(currNames),1);
                tempDists = cell(numel(currNames),1);
                maxCOM = 0;

                for j = numel(currNames):-1:1

                    % Process each region for this current channel/current
                    curSegData = activeMap{i,j};

                    if(~isempty(curSegData)) % if volume isn't empty analyze centroids

                        actProps = regionprops(curSegData,'Centroid');
                        curDistAct = zeros(numel(actProps),1);
                        x = zeros(numel(actProps),1);
                        y = zeros(numel(actProps),1);
                        z = zeros(numel(actProps),1);
                        for k = 1:numel(actProps)
                            curCent = actProps(k).Centroid;
                            x(k) = curCent(1);
                            y(k) = curCent(2);
                            z(k) = curCent(3);
                        end

                        CoM = [(min(x)+max(x))/2,(min(y)+max(y))/2,(min(z)+max(z))/2];

                        for k = 1:numel(actProps)
                            % Calculate distance          
                            curCent = actProps(k).Centroid;
                            curDistAct(k) = sqrt(((CoM(1)-curCent(1))*2.187)^2 + ((CoM(2)-curCent(2))*2.187)^2  + ((CoM(3)-curCent(3))*2)^2 );

    %                        curDistAct(k) = pdist([CoM;curCent],'euclidean')*microns_per_pixel; % 3D distance
    %                        curDistAct(k) = sqrt(sum((CoM-curCent).^2))*microns_per_pixel;
                        end


                        if(numel(actProps)>0)
                            curMaxs(j) = max(curDistAct);
                            curMins(j) = min(curDistAct);
                            curMeans(j) = mean(curDistAct);
                            curSTD(j) = std(curDistAct);

                        end
                        tempDists{j} = curDistAct;


                        % Plot Histogram distribution of neuron distances
                        bins = 0:100:1000;
                        [counts,bins] = hist(curDistAct,bins);
                        maxCOM = max([maxCOM,counts]);
                        if(maxCOM<1)
                            maxCOM = 1;
                        end

                        % Plot individual channel+currenct response
                        curFig = figure();
                        set(0,'CurrentFigure',curFig)
                        clf
                        hold on 
                        barh(bins,counts)
                        xlabel('Count')
                        ylabel('Distance (um)')
                        xlim([0 maxCOM])
                        title(strcat('COM Stimulation Distance Histogram: ',curCh,'_-_-',num2str(currNames(j))))

                        combineFileName = strcat(selpathOutDist,'\StimulationDistancesCOM',curCh,'--',num2str(currNames(j)),'.fig');
                        saveas(curFig,combineFileName);  
                        combineFileName = strcat(selpathOutDist,'\StimulationDistancesCOM',curCh,'--',num2str(currNames(j)),'.tiff');
                        saveas(curFig,combineFileName);  
                        close(curFig)
                    end
                end
                channelDistances = cell2mat(tempDists);


                % Plot individual channel distance activation in 2D
                curFig = figure();
                set(0,'CurrentFigure',curFig)
                clf
                hold on
                errorbar(curMeans,curSTD,'x'); % Plot current varied responses
                plot(curMaxs,'bo')
                plot(curMins,'ro')
                xlabel('Current (uA)')
                xticklabels(num2str(currNames))
                xticks(1:numel(currNames))
                xlim([0.5 (numel(currNames)+0.5)])
                ylabel('Distance (um)')
                title(strcat('Stimulation Distance: ',curCh))

                legend('mean','max','min')
                combineFileName = strcat(selpathOutDist,'\StimulationDistancesCOM',curCh,'.fig');
                saveas(curFig,combineFileName);  
                combineFileName = strcat(selpathOutDist,'\StimulationDistancesCOM',curCh,'.tiff');
                saveas(curFig,combineFileName);  
                close(curFig)

                if(~isempty(channelDistances))
                    chMaxs2D(i) = max(channelDistances);
                    chMins2D(i) = min(channelDistances);
                    chMeans2D(i) = mean(channelDistances);
                    chSTD2D(i) = std(channelDistances);
                else
                    chMaxs2D(i) = 0;
                    chMins2D(i) = 0;
                    chMeans2D(i) = 0;
                    chSTD2D(i) = 0;
                end

            end


            % Plot overall channel reponses
            curFig = figure();
            set(0,'CurrentFigure',curFig)
            clf
            hold on
            errorbar(chMeans2D,chSTD2D,'x'); % Plot current varied responses
            plot(chMaxs2D,'bo')
            plot(chMins2D,'ro')

            xlabel('Channel')
            xticklabels(chNames)
            xticks(1:numel(chNames))
            xlim([0.5 (numel(chNames)+0.5)])
            ylabel('Distance (um)')

            legend('mean','max','min')
            combineFileName = strcat(selpathOutDist,'\StimulationDistancesOverallCOM.fig');
            saveas(curFig,combineFileName);  
            combineFileName = strcat(selpathOutDist,'\StimulationDistancesOverallCOM.tiff');
            saveas(curFig,combineFileName);  
            close(curFig)


            chMaxs2D = zeros(numel(chNames),1);
            chMins2D =  zeros(numel(chNames),1);
            chMeans2D =  zeros(numel(chNames),1);
            chSTD2D =  zeros(numel(chNames),1);
            for i = 1:numel(chNames)
                curCh = chNames{i};

                % calculate distances  for current channel and all currents
                %     channelDistances = zeros(1,1);
                curMaxs = zeros(numel(currNames),1);
                curMins =  zeros(numel(currNames),1);
                curMeans =  zeros(numel(currNames),1);
                curSTD =  zeros(numel(currNames),1);
                tempDists = cell(numel(currNames),1);
                for j = 1:numel(currNames)

                    % Process each region for this current channel/current
                    curSegData = activeMap{i,j};

                    if(~isempty(curSegData)) % if volume isn't empty analyze centroids

                        actProps = regionprops(curSegData,'Centroid');
                        curDistAct = zeros(numel(actProps),1);
                        x = zeros(numel(actProps),1);
                        y = zeros(numel(actProps),1);
                        z = zeros(numel(actProps),1);
                        for k = 1:numel(actProps)
                            curCent = actProps(k).Centroid;
                            x(k) = curCent(1);
                            y(k) = curCent(2);
                            z(k) = curCent(3);
                        end

                        CoM = [(min(x)+max(x))/2,(min(y)+max(y))/2,(min(z)+max(z))/2];

                        for k = 1:numel(actProps)
                            % Calculate distance          
                           curCent = actProps(k).Centroid;
                           curDistAct(k) = sqrt(sum((CoM-curCent).^2))*microns_per_pixel;
                        end

                        if(numel(actProps)>0)
                            curMaxs(j) = max(curDistAct);
                            curMins(j) = min(curDistAct);
                            curMeans(j) = mean(curDistAct);
                            curSTD(j) = std(curDistAct);
                        end
                        tempDists{j} = curDistAct;
                    end
                end
                channelDistances = cell2mat(tempDists);

                % Plot individual channel response
                curFig = figure();
                set(0,'CurrentFigure',curFig)
                clf
                hold on
                errorbar(curMeans,curSTD,'x'); % Plot current varied responses
                plot(curMaxs,'bo')
                plot(curMins,'ro')
                xlabel('Current (uA)')
                xticklabels(num2str(currNames))
                xticks(1:numel(currNames))
                xlim([0.5 (numel(currNames)+0.5)])
                ylabel('Distance (um)')
                title(strcat('Stimulation Distance: ',curCh))

                legend('mean','max','min')
                combineFileName = strcat(selpathOutDist,'\StimulationDistancesCOM',curCh,'.fig');
                saveas(curFig,combineFileName);  
                combineFileName = strcat(selpathOutDist,'\StimulationDistancesCOM',curCh,'.tiff');
                saveas(curFig,combineFileName);  
                close(curFig)

                if(~isempty(channelDistances))
                    chMaxs2D(i) = max(channelDistances);
                    chMins2D(i) = min(channelDistances);
                    chMeans2D(i) = mean(channelDistances);
                    chSTD2D(i) = std(channelDistances);
                else
                    chMaxs2D(i) = 0;
                    chMins2D(i) = 0;
                    chMeans2D(i) = 0;
                    chSTD2D(i) = 0;
                end
            end


            % Plot overall channel reponses
            curFig = figure();
            set(0,'CurrentFigure',curFig)
            clf
            hold on
            errorbar(chMeans2D,chSTD2D,'x'); % Plot current varied responses
            plot(chMaxs2D,'bo')
            plot(chMins2D,'ro')

            xlabel('Channel')
            xticklabels(chNames)
            xticks(1:numel(chNames))
            xlim([0.5 (numel(chNames)+0.5)])
            ylabel('Distance (um)')

            legend('mean','max','min')
            combineFileName = strcat(selpathOutDist,'\StimulationDistancesOverallCOM.fig');
            saveas(curFig,combineFileName);  
            combineFileName = strcat(selpathOutDist,'\StimulationDistancesOverallCOM.tiff');
            saveas(curFig,combineFileName);  
            close(curFig)






            %% Plot neuron specific responses to activation from different channel current combinations
            selpathOutRegDist = strcat(selpathOut,'\RegionDistances');
            mkdir(selpathOutRegDist)

            for i = 1:numel(chNames)
                % Construct framework for output of current sweep analysis
                curCh = chNames{i};

                Responses2D = zeros(10, numel(currNames));% Examine out to 1 mm from contact site
                Responses3D = zeros(20, numel(currNames));% Examine out to 2 mm from contact site

                for curReg = 1:numel(regions)


                    % Pull raster data for region/channel/ (all currents)
                    sInd = ((i-1)*numel(currNames))+ 1;
                    eInd = ((i-1)*numel(currNames))+ numel(currNames);
                    regChanRaster = rasterData(curReg,sInd:eInd);

                    % If region was active during channel stim, analyze it
                    if(sum(regChanRaster)>0)

                        % Calculate the 3D midpoint for the current region if interest
                        curPixels = regions(curReg).PixelList;
                        centroid(1) = mean(curPixels(:,2));
                        centroid(2) = mean(curPixels(:,1));
                        centroid(3) = mean(curPixels(:,3));

                        % Calculate distance to channel contact site
                         curElecPos = [0,0,0];
                        for k = 1:numel(contactsOrdered)
                            if(isequal(num2str(contactsOrdered(k)),curCh))
                                % Current channel matches stimulating channel,
                                % pull position
                                curElecPos = elecPos(k,:);
                            end
                        end



                         Dist2D = sqrt(((centroid(1)-curElecPos(1))*2.187)^2 + ((centroid(2)-curElecPos(2))*2.187)^2);
                         Dist3D = sqrt(((centroid(1)-curElecPos(1))*2.187)^2 + ((centroid(2)-curElecPos(2))*2.187)^2  + ((centroid(3)-curElecPos(3))*2)^2 );

    %                     Dist2D = pdist([centroid(1:2);curElecPos(1:2)],'euclidean')*microns_per_pixel; % 2D distance
                        ind2D = floor(Dist2D/100)+1;
    %                     Dist3D = pdist([centroid;curElecPos],'euclidean')*microns_per_pixel; % 3D distance
                        ind3D = floor(Dist3D/100)+1;

                        if(ind2D>10)
                            ind2D=10;
                        end
                        if(ind3D>20)
                            ind3D=20;
                        end

                        % Add raster to distance maps
                        Responses2D(ind2D,:) = Responses2D(ind2D,:) + regChanRaster;
                        Responses3D(ind3D,:) = Responses3D(ind3D,:) + regChanRaster;
                    end
                end


                % Plot response curve for specific neuron for this
                % chan/curr (Do I need to tie back to the raster???)
                for k = 1:10 % 2D analysis
                    curFig = figure();
                    response = Responses2D(k,:);
                    plot(1:numel(currNames),response)
                    xlabel('Current')
                    xticks(1:numel(currNames))
                    xticklabels(currNames)
                    ylabel('Activations')
                    title(strcat(num2str((k-1)*100),'-',num2str(k*100)))
                    combineFileName = strcat(selpathOutRegDist,'\RegionStimulation-CH',curCh,'-2D-Dist',num2str(k*100),'.fig');
                    saveas(curFig,combineFileName);  
                    combineFileName = strcat(selpathOutRegDist,'\RegionStimulation-CH',curCh,'-2D-Dist',num2str(k*100),'.tiff');
                    saveas(curFig,combineFileName);  
                    close(curFig)   
                end



                for k = 1:20 % 3D Analysis
                    curFig = figure();
                    response = Responses3D(k,:);
                    plot(1:numel(currNames),response)
                    xlabel('Current')
                    xticks(1:numel(currNames))
                    xticklabels(currNames)
                    ylabel('Distance (um)')
                    title(strcat(num2str((k-1)*100),'-',num2str(k*100)))
                    combineFileName = strcat(selpathOutRegDist,'\RegionStimulation-CH',curCh,'-3D-Dist',num2str(k*100),'.fig');
                    saveas(curFig,combineFileName);  
                    combineFileName = strcat(selpathOutRegDist,'\RegionStimulation-CH',curCh,'-3D-Dist',num2str(k*100),'.tiff');
                    saveas(curFig,combineFileName);  
                    close(curFig)   
                end


            end







            %% Plot overall Raster plot with currents segmented into channel groups
            N = size(rasterData,1);
            T = size(rasterData,2);
            [row,col] = find(rasterData);
            raster1D = zeros(numel(row),1);
            for i = 1:numel(row)
                raster1D(i) = (row(i)-1)*size(rasterData,2) + col(i);
            end
            if(~isempty(raster1D))
                rasterplot(raster1D,N,T);
                curFig = gcf;
                hold on
                title(strcat('Overall Raster - Channels & Currents Highlighted'))
                xlabel('Group') 
                ylabel('Active Neurons') 
                labelColor = zeros(numel(chNames),3);
                colorsCUR = distinguishable_colors(numel(chNames)); % Channel Colors
                intensLvl = 0.1:(0.7/numel(currNames)):0.8;
                        % Current Intensity Levels
                for k = 1:numel(chNames)% highlight individual current regions
                    for j = 1:numel(currNames)
                        validInd = ((k-1)*numel(currNames))+ j;
                        curColor = colorsCUR(k,:);
                        if(~strcmp(chNames(k),'0'))
                            intense = intensLvl(j);
                        else
                            intense=0.5;
                        end
                        labelColor(k,:) = curColor;
                        sInd = validInd - 0.5;
                        eInd = validInd + 0.5;
                        patch([sInd eInd eInd sInd], [0 0 2*N 2*N], curColor, 'FaceAlpha', intense, 'EdgeColor', 'none')
                    end
                end
                h = zeros(numel(chNames), 1); % Add legend to data
                    for i = 1:numel(chNames)
                      h(i) = plot(NaN,NaN,'gs','MarkerSize',10,'MarkerFaceColor',labelColor(i,:)); 
                    end
                 legend(h, chNames);
                xlim([0.5, size(rasterData,2)+0.5])
                curFig.WindowState = 'maximized';
                combineFileName = strcat(selpathOut,'\','Overall_Raster.fig');
                saveas(curFig,combineFileName);  
                combineFileName = strcat(selpathOut,'\','Overall_Raster.tiff');
                saveas(curFig,combineFileName);  
                close(curFig)
            end


            % Plot overall Raster plot with currents segmented into channel groups
            rasterActivation = sum(rasterData,1);
            x = 1:numel(rasterActivation);
            curFig = figure();
            set(0,'CurrentFigure',curFig)
            clf
            bar(x,rasterActivation)
            histLim = round(max(rasterActivation)*1.1);
            if(histLim==0)
                histLim=1;
            end
            hold on
            title(strcat('Overall Histogram - Channels & Currents Highlighted'))
            xlabel('Scan') 
            ylabel('Active Neuron Count') 
            labelColor = zeros(numel(chNames),3);
            colorsCUR = distinguishable_colors(numel(chNames)); % Channel Colors
            intensLvl = 0.1:(0.7/numel(currNames)):0.8;
            for k = 1:numel(chNames)% highlight channel & current regions
                for j = 1:numel(currNames)
                    validInd = ((k-1)*numel(currNames))+ j;
                    curColor = colorsCUR(k,:);
                    if(~strcmp(chNames(k),'0'))
                        intense = intensLvl(j);
                    else
                        intense=0.5;
                    end
                    labelColor(k,:) = curColor;
                    sInd = validInd - 0.5;
                    eInd = validInd + 0.5;
                    patch([sInd eInd eInd sInd], [0 0 2*N 2*N], curColor, 'FaceAlpha', intense, 'EdgeColor', 'none')
                end
            end
            h = zeros(numel(chNames), 1); % Add legend to data
            for i = 1:numel(chNames)
              h(i) = plot(NaN,NaN,'gs','MarkerSize',10,'MarkerFaceColor',labelColor(i,:)); 
            end
            legend(h, chNames);
            xlim([0.5, size(rasterData,2)+0.5])
            ylim([0 histLim])
            curFig.WindowState = 'maximized';
            combineFileName = strcat(selpathOut,'\','Overall_Histogram.fig');
            saveas(curFig,combineFileName);  
            combineFileName = strcat(selpathOut,'\','Overall_Histogram.tiff');
            saveas(curFig,combineFileName);  
            close(curFig);










            % Plot monotonically increasing Raster plot with currents segmented into channel groups
            selpathOutMono = strcat(selpathOut,'\Monotonic'); %create new folder
            mkdir(selpathOutMono)




            [row,col] = find(monoRasterData);
            raster1D = zeros(numel(row),1);
            for i = 1:numel(row)
                raster1D(i) = (row(i)-1)*size(monoRasterData,2) + col(i);
            end
            if(~isempty(raster1D))
                rasterplot(raster1D,N,T);
                curFig = gcf;
                hold on
                title(strcat('Overall Raster - Channels & Currents Highlighted'))
                xlabel('Group') 
                ylabel('Active Neurons') 
                labelColor = zeros(numel(chNames),3);
                colorsCUR = distinguishable_colors(numel(chNames)); % Channel Colors
                intensLvl = 0.1:(0.7/numel(currNames)):0.8;
                        % Current Intensity Levels
                for k = 1:numel(chNames)% highlight individual current regions
                    for j = 1:numel(currNames)
                        validInd = ((k-1)*numel(currNames))+ j;
                        curColor = colorsCUR(k,:);
                        if(~strcmp(chNames(k),'0'))
                            intense = intensLvl(j);
                        else
                            intense=0.5;
                        end
                        labelColor(k,:) = curColor;
                        sInd = validInd - 0.5;
                        eInd = validInd + 0.5;
                        patch([sInd eInd eInd sInd], [0 0 2*N 2*N], curColor, 'FaceAlpha', intense, 'EdgeColor', 'none')
                    end
                end
                h = zeros(numel(chNames), 1); % Add legend to data
                    for i = 1:numel(chNames)
                      h(i) = plot(NaN,NaN,'gs','MarkerSize',10,'MarkerFaceColor',labelColor(i,:)); 
                    end
                 legend(h, chNames);
                xlim([0.5, size(monoRasterData,2)+0.5])
                curFig.WindowState = 'maximized';
                combineFileName = strcat(selpathOutMono,'\','Overall_Raster.fig');
                saveas(curFig,combineFileName);  
                combineFileName = strcat(selpathOutMono,'\','Overall_Raster.tiff');
                saveas(curFig,combineFileName);  
                close(curFig)
            end


            %% Plot overall Raster plot with currents segmented into channel groups
            rasterActivation = sum(monoRasterData,1);
            x = 1:numel(rasterActivation);
            curFig = figure();
            set(0,'CurrentFigure',curFig)
            clf
            bar(x,rasterActivation)
            histLim = round(max(rasterActivation)*1.1);
            if(histLim==0)
                histLim=1;
            end
            hold on
            title(strcat('Overall Histogram - Channels & Currents Highlighted'))
            xlabel('Scan') 
            ylabel('Active Neuron Count') 
            labelColor = zeros(numel(chNames),3);
            colorsCUR = distinguishable_colors(numel(chNames)); % Channel Colors
            intensLvl = 0.1:(0.7/numel(currNames)):0.8;
            for k = 1:numel(chNames)% highlight channel & current regions
                for j = 1:numel(currNames)
                    validInd = ((k-1)*numel(currNames))+ j;
                    curColor = colorsCUR(k,:);
                    if(~strcmp(chNames(k),'0'))
                        intense = intensLvl(j);
                    else
                        intense=0.5;
                    end
                    labelColor(k,:) = curColor;
                    sInd = validInd - 0.5;
                    eInd = validInd + 0.5;
                    patch([sInd eInd eInd sInd], [0 0 2*N 2*N], curColor, 'FaceAlpha', intense, 'EdgeColor', 'none')
                end
            end
            h = zeros(numel(chNames), 1); % Add legend to data
            for i = 1:numel(chNames)
              h(i) = plot(NaN,NaN,'gs','MarkerSize',10,'MarkerFaceColor',labelColor(i,:)); 
            end
            legend(h, chNames);
            xlim([0.5, size(monoRasterData,2)+0.5])
            ylim([0 histLim])
            curFig.WindowState = 'maximized';
            combineFileName = strcat(selpathOutMono,'\','Overall_Histogram.fig');
            saveas(curFig,combineFileName);  
            combineFileName = strcat(selpathOutMono,'\','Overall_Histogram.tiff');
            saveas(curFig,combineFileName);  
            close(curFig);




            %% Calculate the trend of monotonically increasing activation
            % i.e. ratio of monotonically increasing ROIs for a given channel
            % verses those that are not consistently activated at higher
            % currents



            rasterSum = sum(rasterData,1);
            monoSum = sum(monoRasterData,1);
            roiCounts = reshape(rasterSum,numel(currNames),numel(chNames))';
            monoRoiCounts = reshape(monoSum,numel(currNames),numel(chNames))';
            roiMeans = mean(roiCounts,1);
            monoROImeans = mean(monoRoiCounts,1);
            roiSTD = std(roiCounts,0,1);
            monoROISTD = std(monoRoiCounts,0,1);

            curFig = figure();
            hold on
            for i = 1:size(roiCounts,1)
                % process each channel

                %calculate circle size
                s_all = roiCounts(i,:);
                s_mono = monoRoiCounts(i,:);
                s_allInd = find(roiCounts(i,:));
                s_monoInd = find(monoRoiCounts(i,:));
                s_all = s_all(s_allInd).*2;
                s_mono = s_mono(s_monoInd).*2;

                yPos = ones(numel(s_all),1)*i;
                yPos2 = ones(numel(s_mono),1)*i;
                scatter(s_allInd,yPos,s_all,'filled','MarkerFaceColor',[0.8000    0.8000    0.8000]);
                scatter(s_monoInd,yPos2,s_mono,'filled','MarkerFaceColor','blue');
            end
            xlim([0 5])
            ylim([0 size(roiCounts,1)+1])
            combineFileName = strcat(selpathOutMono,'\','MonotonicCompareChannels.fig');
            saveas(curFig,combineFileName);  
            combineFileName = strcat(selpathOutMono,'\','MonotonicCompareChannels.tiff');
            saveas(curFig,combineFileName);  
            close(curFig);


            % Plot bar plots of average channel responses
            curFig = figure();
            allMeans = [roiMeans;monoROImeans]';
            allErr = [roiSTD;monoROISTD]';
            h = bar(allMeans);           

            hold on
            xCnt = (get(h(1),'XData') + cell2mat(get(h,'XOffset'))).'; 
            errorbar(xCnt(:), allMeans(:), allErr(:), allErr(:), 'k', 'LineStyle','none')
            xlabel('Stimulation Current (\muA)')
            ylabel('Active Neural ROI counts')
            yyaxis right
            plot(1:size(allMeans,1),allMeans(:,2)./allMeans(:,1))
            combineFileName = strcat(selpathOutMono,'\','MonotonicTrend.fig');
            saveas(curFig,combineFileName);  
            combineFileName = strcat(selpathOutMono,'\','MonotonicTrend.tiff');
            saveas(curFig,combineFileName);  
            close(curFig);

            
            
            combineFileName = strcat(selpathOutMono,'\','MonotonicRegions.mat'); % save key of monotonic regions per channel
            save(combineFileName,'monoRasterData','chNames','currNames','roiCounts','monoRoiCounts'); 
            
    %         er = errorbar(x,data,errlow,errhigh);  % All ROI error bars  
    %         er.Color = [0 0 0];                            
    %         er.LineStyle = 'none';  
    % 
    %         % Monotonic ROI error bars
    %         er = errorbar(,monoROImeans,monoROImeans-monoROISTD,monoROImeans+monoROISTD);  % All ROI error bars  
    %         er.Color = [0 0 0];                            
    %         er.LineStyle = 'none';  



    %         
    %         steps = 1/(numel(chNames)+1);
    %         curFig = figure();
    %         hold on
    %         for i = 1:size(rasterData,2)
    %             % process each current (merge all channels)
    %             currPos = mod(i,numel(currNames));
    %             if(currPos==0)
    %                 currPos=4;
    %             end
    %             curRaster = rasterData(:,i);
    %             rasterInfo = find(curRaster);
    %             rOffset =  steps*floor(i/numel(currNames));
    %             rasterInfo = rasterInfo + rOffset;
    %             xPos = ones(size(rasterInfo))*currPos;
    %             scatter(xPos,rasterInfo,'filled','MarkerFaceColor',[0.8000    0.8000    0.8000],'MarkerFaceAlpha',0.3)
    %         end
    %         
    %         for i = 1:numel(currNames):size(monoRasterData,2) % Plot regions that monotonically increase
    %             % process each current (merge all channels)
    %             curRaster = monoRasterData(:,i:i+(numel(currNames)-1));
    %             [rasterInfo,currMono] = find(curRaster);
    %             rOffset =  steps*floor(i/numel(currNames));
    %             
    %             for k  =1:size(curRaster,1)
    %                 temp = curRaster(k,:);
    %                 inds = find(temp);
    %                 if(numel(inds)>1)
    %                     yPos = ones(numel(inds))*(k+rOffset);
    %                     plot(inds,yPos,'b-o')
    %                 end
    %             end
    %         end
    %         xlim([0.5 4.5])




       %% Calculate overlap and uniquely activated neurons between all channels at each current level
        overlapTable = zeros(numel(chNames),numel(chNames),numel(currNames));
        aOnlyTable = zeros(numel(chNames),numel(chNames),numel(currNames));
        bOnlyTable = zeros(numel(chNames),numel(chNames),numel(currNames));
        strTable = cell(numel(chNames),numel(chNames),numel(currNames));
        selpathOutShared = strcat(selpathOut,'\SharedActivation'); %create new folder
        mkdir(selpathOutShared)
        for i= 1:numel(currNames)
            for j=1:numel(chNames)
                indA =((j - 1)* numel(currNames)) + i; % pull activation data 
                rasterColA = rasterData(:,indA);

                for k=j+1:numel(chNames)
                    %Perform comparison of overlap by examining raster data
                    indB =((k - 1)* numel(currNames)) + i;
                    rasterColB = rasterData(:,indB);

                    overlap = rasterColA;
                    overlap(rasterColB==0)=0;
                    aOnly = rasterColA;
                    aOnly(rasterColB==1)=0;
                    bOnly = rasterColB;
                    bOnly(rasterColA==1)=0;
                    overlapTable(j,k,i) = sum(overlap);
                    aOnlyTable(j,k,i) = sum(aOnly);
                    bOnlyTable(j,k,i) = sum(bOnly);
                    strTable{j,k,i} = strcat(num2str(sum(aOnly)),'/',num2str(sum(overlap)),'/',num2str(sum(bOnly)));

                    % Generate output figure
                    curFig = figure();
                    mainTitle = strcat('Ch',chNames{j},' & Ch',chNames{k},' Curr',num2str(currNames(i)));
                    subTitle = strcat('Ch',chNames{j},'=',num2str(aOnlyTable(j,k,i)),' (Blue)         Ch',chNames{k},'=',num2str(bOnlyTable(j,k,i)),'(Red)        Shared=',num2str(overlapTable(j,k,i)),' (Yellow)'); 
                    title(sprintf(strcat(mainTitle,'\n',subTitle)))
                    ax = gca;
                    ax.TitleFontSizeMultiplier = 1;
                    % Background image
                    sharedSubVol = meanBaseVol;
                    sharedSub = squeeze(max(sharedSubVol,[],3))/max(sharedSubVol,[],'all');%maxRawInten;

                    ax1 = axes;
                    imagesc(sharedSub);
                    colormap(ax1,'gray');
                    daspect([1 1 1])
                    hold on

                    %Assign masked regions to volume
                    colorNeurons = zeros(size(meanBaseVol));
                    for reg = 1:numel(regions)% highlight individual regions
                        curPixels = regions(reg).PixelList;
                        curRegion = false(size(meanBaseVol));
                        for m = 1:size(curPixels,1)
                            curRegion(curPixels(m,2),curPixels(m,1),curPixels(m,3))=1;
                        end   

                        if(aOnly(reg)==1)
%                             A only Neurons - blue
                            colorNeurons(curRegion) = 0.25;

%                             sSurfaces = isosurface(curRegion,0.5);
%                             xyz = sSurfaces.vertices;
%                             tri = sSurfaces.faces;
%                             outXYZ = SurfaceSmooth(xyz, tri, 2, 0.02, 10, 2, 0);
%                             sSurfaces.vertices = outXYZ;
%                             p2 = patch(sSurfaces,'FaceColor','cyan','EdgeColor','none','FaceAlpha',1);
%                             isonormals(curRegion,p2);  

                        elseif(bOnly(reg)==1)
%                             B only Neurons - red
                            colorNeurons(curRegion) = 0.60;


%                             sSurfaces = isosurface(curRegion,0.5);
%                             xyz = sSurfaces.vertices;
%                             tri = sSurfaces.faces;
%                             outXYZ = SurfaceSmooth(xyz, tri, 2, 0.02, 10, 2, 0);
%                             sSurfaces.vertices = outXYZ;
%                             p2 = patch(sSurfaces,'FaceColor','red','EdgeColor','none','FaceAlpha',1);
%                             isonormals(curRegion,p2);  

                        elseif(overlap(reg)==1)
%                             Overlapping neurons - green
                            colorNeurons(curRegion) = 0.45;

%                             sSurfaces = isosurface(curRegion,0.5);
%                             xyz = sSurfaces.vertices;
%                             tri = sSurfaces.faces;
%                             outXYZ = SurfaceSmooth(xyz, tri, 2, 0.02, 10, 2, 0);
%                             sSurfaces.vertices = outXYZ;
%                             p2 = patch(sSurfaces,'FaceColor','yellow','EdgeColor','none','FaceAlpha',1);
%                             isonormals(curRegion,p2);  

                        end
                    end
%                     daspect([1,1,1]); view(3); axis tight; camlight; 
%                     lighting gouraud % Format view
                    colorNeurons = squeeze(max(colorNeurons,[],3));

                    ax2 = axes;
                    imagesc(ax2,colorNeurons,'alphadata',colorNeurons>0);
                    colormap(ax2,'jet');
                    caxis(ax2,[0 1]);
                    ax2.Visible = 'off';
                    linkprop([ax1 ax2],'Position');
                    set(ax1,'XTick',[], 'YTick', [])
                    set(ax,'XTick',[], 'YTick', [])
                    daspect([1 1 1])
                    view(2)

                    curFig.WindowState = 'maximized';
                    combineFileName = strcat(selpathOutShared,'\','SharedActivation_Ch',chNames{j},'_Ch',chNames{k},'_Cur',num2str(currNames(i)),'.fig');
                    saveas(curFig,combineFileName);  
                    combineFileName = strcat(selpathOutShared,'\','SharedActivation_Ch',chNames{j},'_Ch',chNames{k},'_Cur',num2str(currNames(i)),'.tiff');
                    saveas(curFig,combineFileName);  
                    close(curFig);

                end
            end
        end
        combineFileName = strcat(selpathOutShared,'\','SharedRegionData.mat'); % save key of shared activation analysis
        save(combineFileName,'overlapTable','aOnlyTable','bOnlyTable'); 

        chNamesSTR = chNames;
        for i = 1: numel(chNames)
            chNames{i} = erase(chNames{i},'-');
            chNamesSTR{i} = strcat('Ch ',chNames{i});
        end


        % Plot table of shared activations
        for i = 1:numel(currNames)
            T = cell2table(strTable(:,:,i),'VariableNames',chNamesSTR,'RowNames',chNamesSTR);

            % Get the table in string form.
            TString = evalc('disp(T)');
            % Use TeX Markup for bold formatting and underscores.
            TString = strrep(TString,'<strong>','\bf');
            TString = strrep(TString,'</strong>','\rm');
            TString = strrep(TString,'_','\_');
            % Get a fixed-width font.
            FixedWidth = get(0,'FixedWidthFontName');
            % Output the table using the annotation command.

            mainTitle = strcat('Activated neural Population comparison at  ',num2str(currNames(i)),' uA');
            subTitle = '[First channel only / shared / Second channel only]';


            annotation(gcf,'Textbox','String',TString,'Interpreter','Tex',...
            'FontName',FixedWidth,'Units','Normalized','Position',[0 0 1 1]);

            curFig = gcf;
            set(curFig,'Name',sprintf(strcat(mainTitle,'\n',subTitle)),'NumberTitle','off')
            curFig.Position =  [653   668   939   231];

            combineFileName = strcat(selpathOutShared,'\','SharedActivationTable_Cur',num2str(i),'.fig');
            saveas(curFig,combineFileName);  
            combineFileName = strcat(selpathOutShared,'\','SharedActivationTable_Cur',num2str(i),'.tiff');
            saveas(curFig,combineFileName);  
            close(curFig);
        end

        % save data from regional activity to analyze in follow up multisession analysis
        combineFileName = strcat(selpathOut,'\','RegionMetrics.mat'); 
        save(combineFileName,'regions','chNames','currNames','activeMap','rasterData'); 
    end      
end
% When finished loading and processing all source data close MIJI
try
    MIJ.exit % close the ImageJ interface
catch
    disp('Close Error');
end  
fprintf('Finished Processing \n\n');

