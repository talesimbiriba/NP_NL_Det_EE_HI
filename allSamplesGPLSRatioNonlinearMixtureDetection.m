clear all;
close all;

% config. parallel computing
%myCluster = parcluster('local');
%nCores = myCluster.NumWorkers;
nCores = feature('numCores');
isOpen = matlabpool('size') > 0;
if (isOpen)
    matlabpool close;
end


numOfEndmembers = 3;
numOfSamples = 4000;        % 2*numOfSamples will be generated (numOfSamples 
                            %    linear samples and numOfSamples nonlinear ones)

% SNR in dB
SNR = 21;
abundanceVector = [0.3; 0.6; 0.1];

plotMarks = ['--k';'-sk';'-^k'];

nlD = [0.3 0.5 0.8];
% nlD = [0.5];
%plotMarks = ['-sk'];


nparCovFunc = 2; %type of kernel, and number of kernel parameters

ngammas = length(nlD);
dF = 4;    %decimation factor


rocFig = figure;  % GP/LS detector ROCs
rocFigLS = figure;  % Robust detector ROCs
for r=1:ngammas

    linearScores = zeros(numOfSamples,1);
    nonlinearScores = zeros(numOfSamples,1);
    
    linGPErrorNorms = zeros(numOfSamples,1);
    linLSErrorNorms = zeros(numOfSamples,1);
    nonlinGPErrorNorms = zeros(numOfSamples,1);
    nonlinLSErrorNorms = zeros(numOfSamples,1);
    
    
    linKernelPar = zeros(nparCovFunc+1,numOfSamples);
    nonlinKernelPar = zeros(nparCovFunc+1,numOfSamples);

    LSlinearScores = zeros(numOfSamples,1);
    LSnonlinearScores = zeros(numOfSamples,1);
    
    % generate "numOfSamples" samples using the LMM
    modelStr = 'linear';
    
    [Yl,M,a,noiseVar] = createDecimatedDataFromRealEndMembersSNR_NLD(numOfEndmembers,modelStr,numOfSamples,SNR,abundanceVector,nlD(r),dF);
    
    % generate "numOfSamples" samples using the GBM 
    modelStr = 'bilinear'; 
    [Ynl,~,~,~,gamma,k,spectra_names] = createDecimatedDataFromRealEndMembersSNR_NLD(numOfEndmembers,modelStr,numOfSamples,SNR,abundanceVector,nlD(r),dF);
    
    
    progressStepSize = max([1,floor(numOfSamples/100)]);

    % open matlabpool
    
    matlabpool('open', nCores)
    runOnAllJavaMonitorCP;
    ppm = ParforProgMon('Example: ', numOfSamples, progressStepSize, 300, 80);

    
    tic
    % making the tests for the samples
    parfor i=1:numOfSamples,
        
        [linearScores(i),linKernelPar(:,i),linGPErrorNorms(i),linLSErrorNorms(i)] = allSamplesGPLSRatioTest(Yl(:,i),M);
        [nonlinearScores(i), nonlinKernelPar(:,i),nonlinGPErrorNorms(i),nonlinLSErrorNorms(i)] = allSamplesGPLSRatioTest(Ynl(:,i),M);
        
        LSlinearScores(i) = robustTestForNonlinearMixtureDetection(Yl(:,i),M);
        LSnonlinearScores(i) = robustTestForNonlinearMixtureDetection(Ynl(:,i),M);
        
        % To compare a new detector with the detectors used here you should
        % call your detector here and save the liner and nonlinear test
        % statistics (scores).
        
        
        %parfor update
        if mod(i,progressStepSize)==0
            ppm.increment();
        end

    end
    ppm.delete();
    toc

    matlabpool('close')

    %rocploter(linearScores,nonlinearScores,'-k')
    rocploter(nonlinearScores,linearScores,plotMarks(r,:),rocFig);
    %rocploter(LSlinearScores,LSnonlinearScores,plotMarks(r,:),rocFigLs);
    rocploter(LSlinearScores,LSnonlinearScores,plotMarks(r,:),rocFigLS,false);
    %plotROCApproximation(M,a,mean([nonlinKernelPar linKernelPar],2),dF,'k',figHandle);
    if r==2
        FigComp = figure;
        rocploter(nonlinearScores,linearScores,'-sk',FigComp);
        rocploter(LSlinearScores,LSnonlinearScores,'-^k',FigComp,false);
    end
    
end
figure(rocFig)
legend(['NLD = ', num2str(nlD(1))],['NLD = ', num2str(nlD(2))],['NLD = ', num2str(nlD(3))],'location','best')
figure(rocFigLS)
legend(['NLD = ', num2str(nlD(1))],['NLD = ', num2str(nlD(2))],['NLD = ', num2str(nlD(3))],'location','best')
figure(FigComp)
legend('GP','LS')

%% Ploting Statistics Histograms

smax = max([linGPErrorNorms;nonlinGPErrorNorms]);
smin = min([linGPErrorNorms;nonlinGPErrorNorms]);
figure;
subplot(2,1,1);
set(gca,'FontSize',14)
hist(linGPErrorNorms,40);
xlim([smin, smax])
title('Squared Norm of the GP Fitting Error under H_0')
subplot(2,1,2);
set(gca,'FontSize',14)
hist(nonlinGPErrorNorms,40);
xlim([smin, smax])
title('Squared Norm of the GP Fitting Error under H_1')

smax = max([linLSErrorNorms;nonlinLSErrorNorms]);
smin = min([linLSErrorNorms;nonlinLSErrorNorms]);
figure;
subplot(2,1,1);
set(gca,'FontSize',14)
hist(linLSErrorNorms,40);
xlim([smin, smax])
title('Squared Norm of the LS Fitting Error under H_0')
subplot(2,1,2);
set(gca,'FontSize',14)
hist(nonlinLSErrorNorms,40);
xlim([smin, smax])
title('Squared Norm of the LS Fitting Error under H_1')

smax = 1;
smin = min([linearScores;nonlinearScores]);
figure
subplot(2,1,1);
set(gca,'FontSize',14)
hist(linearScores,40);
xlim([smin, smax])
title('Proposed Test Statistics under H_0')
subplot(2,1,2);
set(gca,'FontSize',14)
hist(nonlinearScores,40);
xlim([smin, smax])
title('Proposed Test Statistics under H_1')

% legend(['\gamma = ' num2str(gammas(1))],['\gamma = ' num2str(gammas(2))],['\gamma = ' num2str(gammas(3))],'Location','Best')
% title('Empirical ROCs for the GP Detector')
%legend('Empirical','Theory Approximation','Location','Best')
%title('Empirical and Approximated ROCs for the GP Detector')
