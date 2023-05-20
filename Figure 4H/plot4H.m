%% Plot ICMS thresholds 

% Read csv file
a = readcell('data.csv');
a = a(2:end,1:end-1);
nSOURCES = length(a);

% Modalities: 1-behavior, 2-imaging
threshold_modalities = a(:,5);
modalities = zeros(nSOURCES,1);
for i = 1:nSOURCES
    modality = threshold_modalities{i};
    if strcmp(modality,'Behavior')
        modalities(i) = 1;
    else
        modalities(i) = 2;
    end
end

% Lowest threshold
thresholds = [a{:,8}]';

% Charge per second
chargePerPhase = thresholds;
chargePerSec = [a{:,10}]';

figure
mkr = 'o';
mkr_fcolor = 'white';
mkr_ecolor = 'k';

mkr_size = 400;
lbl_fsize = 14;

letters = {'a','b','c','d','e','f','g','h','i','j','k','l','m','n'};
% Marker colors
mkr_fcolors = [255, 175, 175; 99, 167, 212; 99, 212, 139]/255;
x_shift = 0; y_shift = 0;
for i = 1:nSOURCES
    x = thresholds(i);
    y = chargePerSec(i);
    % our work
    if i == nSOURCES || i == nSOURCES - 1
        mkr = 'pentagram';
        mkr_size = 800;
    end

    mkr_ecolor = 'k';
    face_color = mkr_fcolors(modalities(i),:);
    mkr_falpha = 1;
    scatter(x,y,mkr_size,'Marker',mkr,'MarkerFaceColor',...
        face_color,'MarkerEdgeColor',mkr_ecolor,'MarkerFaceAlpha',mkr_falpha);
    if i < 12
        text(x+x_shift,y+y_shift,sprintf('%s',letters{i}),'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle','FontSize',lbl_fsize)
    elseif i == 12
        text(x+x_shift,y-1,'6 Hz','HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle','FontSize',lbl_fsize)
    else
        text(x+x_shift,y-10,'100 Hz','HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle','FontSize',lbl_fsize)
    end
    hold on;
    x_shift = 0; y_shift = 0;
end

xlabel('Charge/phase (nC/phase)')
ylabel('Charge/s (nC/phase/s)')

set(gca,'fontsize', 13)
set(gca,'fontweight','bold')
set(gca,'linewidth',1.5)

h(1) = scatter(nan,nan,800,'Marker','o','MarkerFaceColor',mkr_fcolors(1,:),'MarkerEdgeColor','k');
h(2) = scatter(nan,nan,800,'Marker','o','MarkerFaceColor',mkr_fcolors(2,:),'MarkerEdgeColor','k');
h(3) = scatter(nan,nan,600,'Marker','pentagram','MarkerFaceColor',mkr_fcolors(1,:),'MarkerEdgeColor','k');
legend(h,{'Behavioral','Imaging','This work'},'LineWidth',1.5,'Location','Southeast')

height = chargePerSec(end) - chargePerSec(end-1);
width = thresholds(end-1) - thresholds(end);
ourPts = [thresholds(end) chargePerSec(end-1) width height];
xx = [thresholds(end-1)*0.95 thresholds(end-1)*0.95 thresholds(end)*1.1 thresholds(end)*1.1];
yy = [chargePerSec(end)*0.9 chargePerSec(end-1)*1.1 chargePerSec(end-1)*1.1 chargePerSec(end)*0.9];
% Log scale
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlim([0.18 5])
ylim([1 2000])
box on
xtix = [0.2:0.2:1 2:4]; 
set(gca,'XTick',xtix);
set(gca,'XTickLabel',xtix)