function [fitresult, gof, curve_x, curve_y] = createExpFit(x, y)

[xData, yData] = prepareCurveData( x, y );

% Set up fittype and options.
ft = fittype( 'a*exp(-b*x)+c', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.MaxFunEvals = 1000;
opts.MaxIter = 1000;
opts.Robust = 'LAR';
opts.StartPoint = [30 0.005 2];
opts.StartPoint = [30 0.0001 2];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

curve_x = xData;
curve_y = feval(fitresult,xData);


