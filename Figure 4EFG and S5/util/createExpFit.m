function [fitresult, gof, curve_x, curve_y] = createExpFit(x0, y0)

[xData, yData] = prepareCurveData( x0, y0 );

% Set up fittype and options.
ft = fittype( 'a*exp(-b*x)+c', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.MaxFunEvals = 1000;
opts.MaxIter = 1000;
opts.Robust = 'LAR';
opts.StartPoint = [100 0.001 2];


% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

curve_x = xData;
curve_y = feval(fitresult,xData);
end

