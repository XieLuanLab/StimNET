% Plot staircase example for single channel

load('data\data.mat');
ch16_idx = find(staircase.Intan_stim_channel == 16);
results = staircase.results{ch16_idx};
n = 12;
trial_col = results(:,1);
current_col = results(:,4);
response_col = results(:,5);
rev_col = results(:,6);
p = plot(trial_col(1:n),current_col(1:n),'k','LineWidth',2);
hold on;
rev_acum = 0;
for i = 1:n  
    response = response_col(i);
    reversal = rev_col(i);
    if reversal == 1
        rev_acum = rev_acum + 1;
        rev_x = trial_col(i);
        rev_y = current_col(i);
        if response == 1
            t = text(rev_x-0.25,rev_y+0.45,sprintf('R%d',rev_acum),'FontSize',12);
        else
            t = text(rev_x-0.25,rev_y-0.4,sprintf('R%d',rev_acum),'FontSize',12);
        end
    end
    colors = [[1 1 1];[0 0 0]];
    color = colors(response+1,:);
    sctr = scatter(trial_col(i),current_col(i),150,'s','filled',...
        'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',color,'LineWidth',2);
end
ylim([0 n+1])
xlim([0 n+1])
xlabel('Trial number','FontSize',20);
ylabel('Current ({\mu}A)','FontSize',20);
yl = yline(1.5,'--r',{'Threshold=1.5{\mu}A'},'LineWidth',2,'FontSize',12);
yl.LabelHorizontalAlignment = 'left';
f=get(gca,'Children');
lgd = legend([f(4),f(6)],'ICMS detected','ICMS not detected');
lgd.FontWeight = 'bold';
set(gca,'fontsize', 12)
set(gca,'fontweight','bold')
set(gca,'linewidth',1.5)
