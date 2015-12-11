function [ figureHandle ] = plotPixelAndEstimation(truePixel, estimatedPixel, legendStr )
%function [ figureHandle ] = plotPixelAndEstimation(truePixel,
%estimatedPixel, legendStr )
%
%   Plot the two pixels passed as arguments truePixel, estimatedPixel, legendStr (optional) and
%   return the figure Handle figureHandle.
%   


figureHandle = figure;
set(gca,'FontSize',14)

plot(truePixel,'color',[0.7 0.7 0.7])
hold on
plot(estimatedPixel,'r')


xlim([1 length(truePixel)]);

if (nargin < 3)
    legend('true pixel','estimated pixel','Location','Best')
else
    legend(legendStr,'Location','Best')
end

xlabel('Bands','FontSize',14);
ylabel('Reflectance','FontSize',14);
grid on

