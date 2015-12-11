function [ handle ] = plotGPPredictionAndErrorBars(gpPrediction, gpVariance)
%[ handle ] = plotGPPredictionAndErrorBars(gpPrediction, gpVariance)
%
%   Creates figure with the gpMeanPrediction and +-1 std error bars.

    
    handle = figure;
    set(gca,'FontSize',14)
    t = 1:length(gpPrediction);
    errorbar(t,gpPrediction,3*sqrt(gpVariance),'color',[0.8 0.8 0.8]);
    hold on;
    plot(t,gpPrediction,'k');
    xlabel('Bands')
    ylabel('Radiance')
    title('GP Prediction and Error bars')
    legend('99.7% confidence region','gp prediction','Location','Best')
    xlim([1 t(end)]);
    grid on;
    
end

