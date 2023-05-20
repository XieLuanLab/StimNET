%% Plot data
% Look at L1-6 
useChargeUnit = 0;
Lidx = 1:5;
scatterColor = 'k';

curve_fits = cell(5,2);
scatterArr = cell(5,2);
gofArr = cell(5,1); % goodness of fit metrics (R^2)

for i = 1:nANIMALS
    data = s_new(i).allData;
    
    x0 = data.daysArr; % session dates converted to indices
    subplot(2,3,i);
    Ldata = data.Ldata(Lidx,:);
 
    b = [];
    for j = 1:numel(Lidx)
        layerData = Ldata{j};
        if useChargeUnit == 1
            layerData = layerData * 1e-6 * 167e-6 / 1e-9;
        end
        if isempty(layerData); continue; end
        for k = 1:size(layerData,1)
            y_scatter = layerData(k,:);
            scatter(x0(1:end-1),y_scatter(1:end-1),12,scatterColor,...
            'filled','MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.3,'jitter','on','jitterAmount',0.5);
            hold on
            scatter(x0,y_scatter,12,scatterColor,...
            'filled','MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.3,'jitter','off');  
            hold on  
        end
        b = [b; layerData];
    end

    y0 = median(b,'omitnan'); % take median across channels 

    % For Mouse3, split into 2 sections, pre-shock and post-shock
    if i == 3
        xpre0 = x0(1:39);
        ypre0 = y0(1:39);
        xpost0 = x0(40:end);
        ypost0 = y0(40:end);
        [~, gof1, curve_xpre, curve_ypre] = createExpFit(xpre0, ypre0); % fit to median threshold points  
        [fitresult, gof2, curve_xpost, curve_ypost] = createExpFit(xpost0, ypost0);
        gofArr{i} = [gof1.rsquare gof2.rsquare]; % R^2 for pre-shock and post-shock fits
        curve_x = [curve_xpre' curve_xpost'];
        curve_y = [curve_ypre' curve_ypost'];
    else
        [fitresult, gof, curve_x, curve_y] = createExpFit(x0, y0);
        gofArr{i} = gof.rsquare;
    end
    curve_fits{i,1} = curve_x;
    curve_fits{i,2} = curve_y;
    plot(curve_x,curve_y,'color',scatterColor,'LineWidth',2);
    hold on
    
    
end
