function figureHandle = rocploter( h0Scores, h1Scores, plotStr, figureHandle, oinvertThreshold )
%
%Usage:
% figureHandle = rocploter( h0Scores, h1Scores, plotStr, figureHandle,
% oinvertThreshold )
%
%   This function plots the empirical ROC curve for the scores h0Scores and
%   h1Scores 
%
%   Variables: 
%
%   h0Scores    ->  name of the ascii file with scores for the H0 hypothesis
%   h1Scores    ->  name of the ascii file with scores for the H1 hypothesis
%   plotStr     ->  plot arguments e.g. '-sk'.
%   figureHandle -> to plot in a existing figure
%   oinvertThreshold -> to invert the inequality.
%
%  Author: Tales Imbiriba
%  Date: April 17, 2013
%
    if nargin < 5
        oinvertThreshold = false;
        if nargin < 4,
            figureHandle = [];
            if nargin < 3
                plotStr = [];
            end
        end
    end
    
    if isempty(plotStr)
        plotStr = '-k';
    end

%     h0Scores = load(h0Scores);
%     h1Scores = load(h1Scores);

    %h0Scores = h0Scores;
    %h1Scores = h1Scores;

    minScores = min([min(h0Scores), min(h1Scores)]);
    maxScores = max([max(h0Scores), max(h1Scores)]);


    thresholdStep = (maxScores-minScores)/100;

    
    t = minScores:thresholdStep:maxScores;
    %t = minScores:thresholdStep:maxScores+1;

    PFA = zeros(size(t));
    PD = zeros(size(t));

    if (oinvertThreshold)
        for i=1:length(t),
            PFA(i) =  length(h0Scores(h0Scores <t(i)));
            PD(i) = length(h1Scores(h1Scores <t(i)));
        end
    else
        for i=1:length(t),
            PFA(i) =  length(h0Scores(h0Scores >t(i)));
            PD(i) = length(h1Scores(h1Scores >t(i)));
        end
    end

    PFA = PFA/length(h0Scores);
    PD = PD/length(h1Scores);
    
    if isempty(figureHandle)
    figureHandle = figure;
    else
        figure(figureHandle);
        hold on;
    end
    set(gca,'FontSize',14)
    plot(PFA,PD, plotStr);
    %plot(PFA,PD);
    grid on;
    xlabel('Probability of False Alarm','fontsize',14);
    ylabel('Probability of Detection','fontsize',14);
    title('ROC Curve','fontsize',14);
    set(gca,'XTick',[0:0.1:1])
    set(gca,'YTick',[0:0.1:1])

end