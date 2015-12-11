%function rocploter( fileNameH0, fileNameH1, thresholdInfLimit, thresholdSupLimit, numberOfSteps )
function rocploter( fileNameH0, fileNameH1, bInvertIneq, colorString, lineWidth, numberOfSteps )
%
%Usage:
% rocploter( fileNameH0, fileNameH1, thresholdInfLimit, thresholdSupLimit,
% thresholdStep )
%
%  Variables: 
%
%  fileNameH0 -> name of the ascii file with scores for the H0 hypothesis
%
%  fileNameH1 -> name of the ascii file with scores for the H1 hypothesis
%
%  thresholdInfLimit -> inferior Threshold limit
%
%  thresholdSupLimit -> superior Threshold limit
%
%  thresholdStep -> step of variation in the interval [thresholdInfLimit,thresholdSupLimit]
%this last variable is optional and if not used will be set to
%thresholdStep = 0.0001;
%
%  Author: Tales Imbiriba
%  Date: April 17, 2013
%


if (nargin < 6) 
    numberOfSteps = 200;
end
if (nargin < 5)
    lineWidth = 1;
end
if (nargin < 4)
    colorString = '';
end

if (nargin < 3)
    bInvertIneq = 0;
end


    h0Scores = load(fileNameH0);
    h1Scores = load(fileNameH1);
    
    [l,c] = size(h0Scores);
    if(c>l)
        h0Scores = h0Scores';
    end
    [l,c] = size(h1Scores);
    if(c>l)
        h1Scores = h1Scores';
    end
    
    maxt = max([h0Scores;h1Scores]);
    mint = min([h0Scores;h1Scores]);

    %t = [maxt:0.00001:0.5];
    t = linspace(mint,maxt,numberOfSteps);

    PFA = zeros(size(t));
    PD = zeros(size(t));

    if bInvertIneq == 0
        for i=1:length(t),

            PFA(i) =  length(h0Scores(h0Scores > t(i)));
            PD(i) = length(h1Scores(h1Scores > t(i)));
        end
    else
        for i=1:length(t),

            PFA(i) =  length(h0Scores(h0Scores < t(i)));
            PD(i) = length(h1Scores(h1Scores < t(i)));
        end
        
    end

    PFA = PFA/length(h0Scores);
    PD = PD/length(h1Scores);
    plot(PFA,PD,colorString,'LineWidth',lineWidth);
   % grid
    %xlabel('PFA','fontsize',16);
    %ylabel('PD','fontsize',16);
    xlabel('Probability of False Alarm','fontSize',16)
    ylabel('Probability of Detection','fontsize',16)
    title('ROC Curve','fontsize',16);


end