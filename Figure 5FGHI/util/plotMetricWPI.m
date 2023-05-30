function mdl = plotMetricWPI(weekArr,data,plotColor)

    uniqueWeeks = unique(weekArr);
    meanWeekData = cellfun(@(x) mean(x,2,'omitnan'), data,'UniformOutput',false);
    meanTrace = cell2mat(cellfun(@(x) mean(x,1,'omitnan'), meanWeekData,'UniformOutput',false));

    y = meanTrace(~isnan(meanTrace));
    x = uniqueWeeks(~isnan(meanTrace))*7;

    y_scatter = cell2mat(meanWeekData(~isnan(meanTrace))');

    xq = uniqueWeeks(1):max(uniqueWeeks);
    xq = xq * 7;
    yq = interp1(x,y,xq);
 
    stdshade(y_scatter',0.2,plotColor,x,1,1.5);
    hold on;
    for i = 1:size(meanWeekData,1)
        t = uniqueWeeks(i)*7;
        scatter(t,meanTrace(i),20,plotColor,'filled');
        hold on;
    end
    hold on;
    plot(xq,yq,'Color',plotColor);
    xlabel('Days post-implantation')
    hold on;
end