% Plot Figure 5F spikes
% Text and scale bar added in Powerpoint

load('data/spikes.mat')
addpath('util')
color_arr = get(gca,'colororder');
close;
ylimmax = 80;
ylimmin = -300;
tiledlayout(5,3,'TileSpacing','none');

for row = 1:5
    for col = 1:3
        nexttile;
        x = (1:48)/3e1 + (col-1)*2;
        spikes = spikeArr{row,col};
        h=stdshade(spikes',0.3, color_arr(row,:), x, 1, 1);

        xlim([x(1)-1 x(end)+1])
        ylim([ylimmin ylimmax]);

        hold on;
        axis off
        if row == 5 && col == 2
            xline(-1/3e1)
            line([0 20/3e1], [-200 -200])
            line([x(1)+25/3e1 x(1)+30/3e1], [0 0])
            line([x(1)+25/3e1 x(1)+30/3e1], [-200 -200])
        end
        % scale bars
        if row == 5 && col == 1
            xline(-1/3e1)
            line([0 30/3e1], [-200 -200])
        end
    end
end