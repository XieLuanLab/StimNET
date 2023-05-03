%% ICMS detection threshold task stage 1 - premature turns and post-stim turns

fileName = 'raw.jsonable';
fileID = fopen(fullfile('data',fileName),'r');

raw = fread(fileID,inf); 
str = char(raw'); 

% Parse raw data for trial start, response period start, reward 
key_trial_start  = '"Trial start timestamp":';
key_trial_end = '"Trial end timestamp":';
idx_trial_start = strfind(str, key_trial_start);
idx_trial_end = strfind(str, key_trial_end);

key_response_start   = '"response_period": [[';
idx_response_start = strfind(str, key_response_start);

idx_response_period = strfind(str,'"response_period":');
rp = str(idx_response_period(1):idx_response_period(1)+20);
response_period = str2double(regexp(rp,'(?<="response_period":[^0-9]*)[0-9]*\.?[0-9]+', 'match')); 

idx_ITI = strfind(str,'"ITI":');
iti_str = str(idx_ITI(1):idx_ITI(1)+20);
ITI = str2double(regexp(iti_str,'(?<="ITI":[^0-9]*)[0-9]*\.?[0-9]+', 'match')); 

key_response_end1 = ']], "reward":';
idx_response_end1 = strfind(str, key_response_end1);
if isempty(idx_response_end1)
    key_response_end1 = ']], "rewardA":';
    idx_response_end1 = strfind(str, key_response_end1);
end

pattern = '"PSP":';
idx_PSP = strfind(str,'"PSP":');
psp_str = char(extractBetween(str(idx_PSP(1)+1:idx_PSP(1)+20),'[',']'));
stage_flag = 1; 
PSP = psp_str;
PSP_turn_str1 = '"wait_before_PSP_reset": ';   
PSP_turn_str2 = '], "stim_on":';
PSP_turn_idx1 = strfind(str, PSP_turn_str1);
PSP_turn_idx2 = strfind(str, PSP_turn_str2);

reward_times = []; 
no_reward_counter = 0;
trial_start_times = [];
response_start_times = [];
iti_threshold_crossings = [];
nTrials = numel(idx_trial_start);
for i = 1:nTrials
    trial_start = sscanf(str(idx_trial_start(i) + length(key_trial_start):idx_trial_end(i)), '%g', 1);
    response_start_parse = sscanf(str(idx_response_start(i)+ length(key_response_start):idx_response_end1(i)), '%g', 1);
    trial_start_times = [trial_start_times trial_start];
    response_start_times = [response_start_times trial_start+response_start_parse];

    key1 = '"reward": ';
    key2 = ', "ITI": [[';

    idx1 = strfind(str, key1);
    if isempty(idx1)
        key1 = '"rewardA": ';
        idx1 = strfind(str, key1);
    end
    idx2 = strfind(str, key2);
    s = str(idx1(i):idx1(i)+40);
    r = str2double(regexp(s,'(?<="reward[^0-9]*)[0-9]*\.?[0-9]+', 'match'));
    if isempty(r)
        no_reward_counter = no_reward_counter + 1;
    end
    reward_times = [reward_times trial_start + r];
    PSP_str = str(PSP_turn_idx1(i):PSP_turn_idx2(i));
    PSP_turn_time = regexp(PSP_str,'\d+\.?\d*','match');
    iti_threshold_crossings = [iti_threshold_crossings trial_start + str2double(PSP_turn_time(1:2:end))];
end

    
% Rasterplot and peri-event rate histogram
h = figure;
PRESTIM_WIN = response_period; % seconds
POSTSTIM_WIN = response_period;
ts = {};
nTrials = numel(idx_trial_start);
for i = 1:nTrials

        trial_start = trial_start_times(i);
        response_start = response_start_times(i);

        if i == 1
            prev_response_start = 0;
        else
            prev_response_start = response_start_times(i-1);
        end

        if i == nTrials
            next_trial_start = trial_start_times(i) + 99999;
        else
            next_trial_start = trial_start_times(i+1);
        end

        premature_t = iti_threshold_crossings(iti_threshold_crossings > prev_response_start & ...
        iti_threshold_crossings < response_start) - response_start;
        reward_t = reward_times(reward_times > response_start & reward_times < next_trial_start) - response_start;
        ts{i} = [premature_t reward_t];
end


preMature_x_vec = [];
preMature_y_vec = [];
postStim_x_vec = [];
postStim_y_vec = [];
for i = 1:nTrials
    timestamps_i = ts{i};
    if isempty(timestamps_i); continue; end
    preMature_x = timestamps_i(timestamps_i <= 0);
    postStim_x = timestamps_i(timestamps_i > 0);

    preMature_x_vec = [preMature_x_vec preMature_x];
    preMature_y_vec = [preMature_y_vec ones(1,length(preMature_x)) .* i];

    postStim_x_vec = [postStim_x_vec postStim_x];
    postStim_y_vec = [postStim_y_vec ones(1,length(postStim_x)) .* i];

end
scatter(preMature_x_vec,preMature_y_vec,8,'filled','Color','blue');
hold on;
scatter(postStim_x_vec,postStim_y_vec,8,'filled');

xlabel('Time (s)');
ylabel('Trial number');
ylim([0 nTrials+1])
xlim([-PRESTIM_WIN POSTSTIM_WIN])
xline(0,'color',[0.466,0.674,0.199],'Linewidth',3)
set(gca,'fontsize', 12)
set(gca,'fontweight','bold')
set(gca,'linewidth',1.5)