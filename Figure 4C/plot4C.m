% Plot psychometric curve

addpath('psignifit'); % https://github.com/wichmann-lab/psignifit

load('data\data.mat');
unpackStruct = @(s) cellfun(@(name) assignin('base',name,getfield(s,name)),fieldnames(s));
unpackStruct(data);

nSTIM_PARAMS = numel(currents);
correct_trials = zeros(nSTIM_PARAMS,1);
total_trials = zeros(nSTIM_PARAMS,1);
for i=1:size(trial_responses_arr,1)
    correct_trials(i) = sum(trial_responses_arr{i});
    total_trials(i) = length(trial_responses_arr{i});
end
input_data = [currents' correct_trials total_trials];

options = struct; 
options.sigmoidName = 'norm';   
options.expType = 'YesNo';
options.plotThresh = 'false';
result = psignifit(input_data,options);
plotPsych(result)
xlabel('Current ({\mu}A)','FontSize',14)
ylabel('Proportion detected','FontSize',14)

xlim([0 12])
ylim([-0.05 1.05])
box on
set(gca,'linewidth',1.5)
set(gca,'TickDir','in');
