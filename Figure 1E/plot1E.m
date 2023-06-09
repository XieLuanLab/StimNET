% Figure 1E

% Plot multiple VT curves from multiple files
% Define list of file names to read

files = {'10uA.DTA',...
    '20uA.DTA',...
    '30uA.DTA',...
    '40uA.DTA',...
    '50uA.DTA'};

nFILES = numel(files);
for i = 1:nFILES

    [T, Vf] = importDTAfile(fullfile('data',files{i}));
    plot(1e6*T,Vf,'LineWidth',2)
    hold on;
    
end
legend({'10 \muA','20 \muA','30 \muA','40 \muA','50 \muA'},'Location','southeast')

set(gca,'fontsize', 12)
set(gca,'fontweight','bold')
set(gca,'linewidth',1.5)

xlabel('Time (\mus)')
ylabel('Potential V vs. Ag|AgCl')
