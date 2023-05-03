% Plot multiple CV curves from multiple files
% Define list of file names to read
files = {'0pulses.DTA',...
        '2.5M.DTA',...
        '5M.DTA',...
        '7.5M.DTA',...
        '10M.DTA',...
        '15M.DTA',...
        '50M.DTA'};

nFILES = numel(files);
% Read File
for i = 1:nFILES
fid = fopen(fullfile('data',files{i}),'rt');
impData = textscan(fid, '%s', 'Delimiter', ' ');
fclose(fid);
impData = impData{1};
%Remove empty lines
impData(cell2mat(cellfun(@isempty,impData,'UniformOutput',0))) = [];
%Split data strings
strData = cellfun(@strsplit, impData, 'UniformOutput',0);
%Search for the keyword curve in the data
idxCurve = cellfun(@(X) any(cell2mat(strfind(X, 'CURVE'))), strData, 'UniformOutput',0);
%Return indicies of the keyword curve
idxCurve = find(~cell2mat(cellfun(@(X) X == 0,idxCurve,'UniformOutput',0)));
%Save the data in a new variable
for j = 1:length(idxCurve)-1
    eval(['Data.Curve' num2str(j) ' = strData(idxCurve(j)+6 : idxCurve(j+1)-1);']);
end
eval(['Data.Curve' num2str(j+1) ' = strData(idxCurve(j)+6 : end);']);

% Smooth 50M curves
if i == 7
    vs = [];
    ims = [];
    for k = [3 4 5 6]
        curve = eval(['Data.Curve' num2str(k)]);
        vs = [vs cell2mat(cellfun(@(X) str2double(X{3}), curve, 'UniformOutput', 0))];
        ims = [ims cell2mat(cellfun(@(X) str2double(X{4}), curve, 'UniformOutput', 0))];
    end
    v = mean(vs,2);
    im = smooth(mean(ims,2));
else
    curveNum = 3;
    curve = eval(['Data.Curve' num2str(curveNum)]);
    t = cell2mat(cellfun(@(X) str2double(X{2}), curve, 'UniformOutput', 0));
    v = cell2mat(cellfun(@(X) str2double(X{3}), curve, 'UniformOutput', 0));
    im = cell2mat(cellfun(@(X) str2double(X{4}), curve, 'UniformOutput', 0));
end

% Plot as two separate curves
vmax_idx = find(v == max(v));
v1 = v(1:vmax_idx);
i1 = im(1:vmax_idx);
v2 = v(vmax_idx+1:end);
i2 = im(vmax_idx+1:end);
% Need to connect two curves together
v1c = [v2(end);v1];
i1c = [i2(end);i1];
v2c = [v1c(end);v2];
i2c = [i1c(end);i2];

% Plot CV for individual file
linewidth = 3;
cmap = parula(nFILES + 3);
color_idx = i;
% Multiply current by 1e9 to get nA
y1 = i1c * 1e9;
y2 = i2c * 1e9;
plot(v1c,y1,'Color',cmap(color_idx,:),'LineWidth',linewidth);
hold on;
plot(v2c,y2,'Color',cmap(color_idx,:),'LineWidth',linewidth);
hold on;
end

h1 = zeros(nFILES,1);
for i = 1:nFILES
    h1(i) = plot(0,0,'Color',cmap(i,:),'LineWidth',linewidth);
    hold on;
end

lgdNames = {'0 pulses','2.5M pulses','5M pulses','7.5M pulses','10M pulses','15M pulses','50M pulses'};
legend(h1,lgdNames,'Interpreter','none','Location','southeast');
ylabel('Current (nA)')
xlabel('Potential V vs. Ag|AgCl')
xlim([-0.8,1]);
ylim([-150,100]);
set(gca,'fontsize', 12)
set(gca,'fontweight','bold')
set(gca,'linewidth',1.5)