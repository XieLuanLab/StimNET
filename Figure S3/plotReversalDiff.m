openfig('Behavior S1.fig')
hSc=findobj(gcf,'Type','scatter');

r4 = hSc(1).YData;
r3 = hSc(2).YData;
r2 = hSc(3).YData;
r1 = hSc(4).YData;

d12 = abs(r1 - r2);
d23 = abs(r2 - r3);
d34 = abs(r3 - r4);
figure;
h = histogram([d12 d23 d34],'Normalization','probability');
h.LineWidth = 1;
xticks(1:4);
% xticklabels({'1','2','3','4'})

xlabel('Difference between reversals ({\mu}A)')
ylabel('Proportion')

box on
ylim([0 1])
set(gca,'fontsize', 15)
set(gca,'linewidth',1.5)