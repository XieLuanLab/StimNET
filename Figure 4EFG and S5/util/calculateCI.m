function [mu,CI] = calculateCI(data,alpha)
    SEM = std(data,'omitnan')/sqrt(length(data)); % Standard error
    width = (100 - alpha) / 100;
    left = width/2;
    right = 1-(width/2);
    ts = tinv([left  right],length(data)-1); % T-Score
    mu = mean(data,'omitnan');
    CI = mu + ts*SEM;       
    % https://www.mathworks.com/matlabcentral/answers/159417-how-to-calculate-the-confidence-interval
end