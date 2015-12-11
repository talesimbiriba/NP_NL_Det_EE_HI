function [endmembers, nlPixIdx] = iterativeEndmemberEstAndNonlinDetect(Y,numOfEndmembers, PFA, maxNumOfRuns, relaxingDetectionFactor,useJavaProgressMonitor, realM)
%[endmembers, nlPixIdx] = iterativeEndmemberEstAndNonlinDetect(Y,numOfEndmembers, 
% PFA, maxNumOfRuns, relaxingDetectionFactor,useJavaProgressMonitor, realM)
%   This function performs 

oplotSimplex=1;

if nargin < 7
    oplotSimplex=0;
    if nargin <6
        useJavaProgressMonitor=1;
        if nargin <5
            relaxingDetectionFactor = 0.9;
            if nargin < 4
                maxNumOfRuns = 10;
                if nargin < 3
                    PFA = 0.01;
                    if nargin <2
                        disp('Usage Error!');
                        return;
                    end
                end
            end
        end
    end
end

Ytmp=Y;

count = 0;
maxStats = 1;
minStats = 0;

statsFig = figure;
simplexFig = figure;
% useJavaProgressMonitor=1;


% Estimate the endmembers using Minimum volume simplex
Mest = MVES(Ytmp,numOfEndmembers,0);
t = getThresholdForGivenPFAUsingLinearImageApproach(Y,Mest,PFA,numOfEndmembers,useJavaProgressMonitor);
relaxingFactorIncrement = (1-relaxingDetectionFactor)/maxNumOfRuns;

while ((maxStats-minStats)>0.05 &&  count < maxNumOfRuns)
    
    disp(['Iteration number ', num2str(count+1)]);
    
    % relaxing the detection threshold
    relaxedDetThreshold = relaxingDetectionFactor*t;
    
    % Control IDX
    controlIdx=[1:size(Ytmp,2)]';
    
    
    %%  2 - perform detection

    [stats,Aest,kpar] = computeRatioNLDetectiosStats(Ytmp,Mest,[],useJavaProgressMonitor);    
    
    % Detect nonlinear pixels with relaxed decision threshold
    nlPixelsIdx = find(stats<=relaxedDetThreshold);
    
    % plot statistics 
    figure(statsFig);clf;
    plot(stats)
    hold on
    plot(nlPixelsIdx,stats(nlPixelsIdx),'sr')
    
    
    controlIdx(nlPixelsIdx) = [];
    %plot3(tt(smallestStatsIdx,1),tt(smallestStatsIdx,2),tt(smallestStatsIdx,3),'om')
    
    %rocploter(stats(1:cNl),stats(cNl+1:end),pltStrgs(count+1,:),figHandle,'true')
    %cNl = cNl - length(find(smallestStatsIdx<Nl));
    %%  3 - exclude pixels detected as nonlinear and return to number 1
    if (oplotSimplex)
        clf(simplexFig)
        draw3dSimplex(realM,simplexFig);
        hold on;
        mm = realM'*Y;
        scatter3(mm(1,:),mm(2,:),mm(3,:),'.k')
        scatter3(mm(1,nlPixelsIdx),mm(2,nlPixelsIdx),mm(3,nlPixelsIdx),'ob')
        mm = realM'*Mest;
        scatter3(mm(1,:),mm(2,:),mm(3,:),'og','fill')
        pause(1);
    end
    
    Ytmp(:,nlPixelsIdx)=[];
    Mest = MVES(Ytmp,numOfEndmembers,0);
    
    
    
    maxStats = max(stats);
    minStats = min(stats);
    relaxingDetectionFactor = relaxingDetectionFactor + relaxingFactorIncrement;
    count = count +1;
    
    
     if(isempty(nlPixelsIdx))
        while (isempty(nlPixelsIdx) &&  count < maxNumOfRuns)          
            relaxedDetThreshold = relaxingDetectionFactor*t;            
            nlPixelsIdx = find(stats<=relaxedDetThreshold);
            % cannot discard more than 20% of the remaining samples
%             if length(nlPixelsIdx) > 0.2*size(Ytmp,2)
%                nlPixelsIdx = findSmallestNumbers(stats,0.2);
%             end
            %smallestStatsIdx = findSmallestNumbers(stats,p);
%             pause(0.1);
            clf(statsFig);
            plot(stats)
            hold on
            plot(nlPixelsIdx,stats(nlPixelsIdx),'sr')
            
            if (oplotSimplex)
                clf(simplexFig)
                draw3dSimplex(realM,simplexFig);
                hold on;
                mm = realM'*Y;
                scatter3(mm(1,:),mm(2,:),mm(3,:),'.k')
                scatter3(mm(1,nlPixelsIdx),mm(2,nlPixelsIdx),mm(3,nlPixelsIdx),'ob')
                mm = realM'*Mest;
                scatter3(mm(1,:),mm(2,:),mm(3,:),'og','fill')
                pause(1);
            end
            
            relaxingDetectionFactor = relaxingDetectionFactor + relaxingFactorIncrement;
            count = count +1;
        end
    end
    
    
end 

if (oplotSimplex)
    clf(simplexFig)
    draw3dSimplex(realM,simplexFig);
    hold on;
    mm = realM'*Y;
    scatter3(mm(1,:),mm(2,:),mm(3,:),'.k')
    scatter3(mm(1,nlPixelsIdx),mm(2,nlPixelsIdx),mm(3,nlPixelsIdx),'ob')
    mm = realM'*Mest;
    scatter3(mm(1,:),mm(2,:),mm(3,:),'og','fill')
    pause(1);
end

endmembers = Mest;
[ nlPixIdx, stats, kpar, threshold ] = GPLSRatioDetection(Y, Mest, PFA, useJavaProgressMonitor);

end

