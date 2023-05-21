%% Plot voltage transients
[x1,y1] = importDTAfile('data/Mouse1_Jan09.dta');
[x2,y2] = importDTAfile('data/Mouse1_Jul07.dta');

% Convert to microseconds
x1 = 1e6 * x1;
x2 = 1e6 * x2;

plot(x1,y1,'LineWidth',2,'Color','k');
hold on;
plot(x2,y2,'LineWidth',2,'Color','r');
xlabel('Time ({\mu}s)')
ylabel('Potential (V)')

set(gca,'fontsize', 20)
set(gca,'fontweight','bold')
set(gca,'linewidth',2)
ylim([-2.2 2])
legend('Week 16','Week 42','','','Location','southeast')