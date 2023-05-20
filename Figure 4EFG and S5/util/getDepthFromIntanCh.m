function [intan2depth,depth_per_ch,nELECS_OUT] = getDepthFromIntanCh(all_stim_ch,animalID)

    if strcmp(animalID,'ICMS37')
        load('channel maps/32LinearChMap_misaligned.mat','channel_map');
    else
        load('channel maps/32LinearChMap.mat','channel_map');
    end
    
    map = channel_map(:,2);

    nSCH = numel(all_stim_ch); % number of stim channels

    switch animalID
        case 'ICMS24'
        horz = 660; vert = 933;
        insertion_angle = rad2deg(atan(horz/vert));
        pitch = 60; % distance between electrodes in um
        nELECS_OUT = 12;
        offset = 0;
        depth_per_ch = pitch * cos(deg2rad(insertion_angle));
    
        case 'ICMS37'
        horz = 840; vert = 1454;
        insertion_angle = rad2deg(atan(horz/vert));
        pitch = 60; % distance between electrodes in um
        nELECS_OUT = 0;
        offset = 260;
        depth_per_ch = pitch * cos(deg2rad(insertion_angle));
    
        case 'ICMS42'
        horz = 810; vert = 1400;
        insertion_angle = rad2deg(atan(horz/vert));
        pitch = 60; % distance between electrodes in um
        nELECS_OUT = 17;
        offset = 0;
        depth_per_ch = pitch * cos(deg2rad(insertion_angle));

        case 'ICMS44'
        horz = 660; vert = 933;
        insertion_angle = rad2deg(atan(horz/vert));
        pitch = 60; % distance between electrodes in um
        nELECS_OUT = 9;
        offset = 0;
        depth_per_ch = pitch * cos(deg2rad(insertion_angle));

        case 'ICMS49'
        horz = 720; vert = 1247;
        insertion_angle = rad2deg(atan(horz/vert));
        pitch = 60; % distance between electrodes in um
        nELECS_OUT = 9;
        offset = 0;
        depth_per_ch = pitch * cos(deg2rad(insertion_angle));
            
    end

    intan2depth = zeros(nSCH,2);
    for ch = 1:nSCH
        Intan_ch = all_stim_ch(ch);
        % depth index ranges from 1-32 where 1 shallowest, 32 deepest
        depth_idx = find(map == Intan_ch);  
        % sanity check: depth_idx of nELECS_OUT + 1 should give depth of
        % depth_per_ch
        depth_um = round((depth_idx - nELECS_OUT) * depth_per_ch) + offset;
        intan2depth(ch,:) = [Intan_ch depth_um];
    end

end

