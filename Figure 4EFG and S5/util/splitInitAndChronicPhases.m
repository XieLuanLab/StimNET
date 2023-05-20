%% Find derivative that drops below threshold and split into two phases
dTh = 0.03; % uA/day threshold 
if useChargeUnit 
    dTh = dTh * 1e-6 * 167e-6 / 1e-9;
end
boundaries = cell(5,1);
timeToReachChronic = cell(5,1);
for i = 1:nANIMALS
subplot(2,3,i)
hold on;
initColor = 'black';
chronicColor = '#77CCFF';
fitColor = 'k';
if i == 3
    X = curve_fits{i,1};
    X1 = X(1:39);
    X2 = X(40:end);
    Y = curve_fits{i,2};
    Y1 = Y(1:39);
    Y2 = Y(40:end);
    y_idx1 = find(abs(gradient(Y1)) < dTh);
    y_idx2 = find(abs(gradient(Y2)) < dTh) + 39;
    boundary1 = y_idx1(1);
    boundary2 = y_idx2(1);
    xline(X(boundary1),'LineWidth',1)
    hold on;
    xline(X(39),'LineWidth',1)
    hold on;
    xline(X(boundary2),'LineWidth',1)
    boundary = [boundary1 boundary2];

    yo = ylim;
    yl = ylim + 1;
    patchV = [0 0; 0 yl(2); X(boundary1) yl(2); X(boundary1) 0]; patchF = [1 2 3 4];
    S1 = patch('Faces',patchF,'Vertices',patchV,'FaceColor',initColor,'FaceAlpha',.1);
    S1.EdgeColor = 'none';
    hold on;
    patchV = [X(boundary1) 0; X(boundary1) yl(2); X(39) yl(2); X(39) 0];
    S2 = patch('Faces',patchF,'Vertices',patchV,'FaceColor',chronicColor,'FaceAlpha',.2);
    S2.EdgeColor = 'none';
    patchV = [X(39) 0; X(39) yl(2); X(boundary2) yl(2); X(boundary2) 0];
    S3 = patch('Faces',patchF,'Vertices',patchV,'FaceColor',initColor,'FaceAlpha',.1);
    S3.EdgeColor = 'none';
    hold on;
    patchV = [X(boundary2) 0; X(boundary2) yl(2); X(end)+1 yl(2); X(end)+1 0];
    S4 = patch('Faces',patchF,'Vertices',patchV,'FaceColor',chronicColor,'FaceAlpha',.2);
    S4.EdgeColor = 'none';
else
    X = curve_fits{i,1};
    Y = curve_fits{i,2};
    y_idx = find(abs(diff(Y)) < dTh);
    boundary = y_idx(1);
    xline(X(boundary),'LineWidth',1)

    yo = ylim;
    yl = ylim + 1;
    patchV = [0 0; 0 yl(2); X(boundary) yl(2); X(boundary) 0];
    patchF = [1 2 3 4];
    S1 = patch('Faces',patchF,'Vertices',patchV,'FaceColor',initColor,'FaceAlpha',.1);
    S1.EdgeColor = 'none';
    hold on;
    patchV = [X(boundary) 0; X(boundary) yl(2); X(end)+1 yl(2); X(end)+1 0];
    patchF = [1 2 3 4];
    S2 = patch('Faces',patchF,'Vertices',patchV,'FaceColor',chronicColor,'FaceAlpha',.2);
    S2.EdgeColor = 'none';
    timeToReachChronic{i} = X(boundary);
end
plot(X,Y,'Color',fitColor,'LineWidth',2);
title(animalIDs{i})

boundaries{i} = boundary;
xlim([0 X(end)+1])
ylim(yo);
if i == 1
    xlabel('Day')
    ylabel('Threshold ({\mu}A)')


end
box on
ax = gca;
ax.LineWidth = 1;
ax.XAxis.TickLength = [0 0];
ax.YAxis.TickLength = [0 0];
end
