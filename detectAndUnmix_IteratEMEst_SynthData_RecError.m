clear all;
close all;


SNR = 21;                   % SNR in dB
oestimateM = 1;             % estimate the endmember matrix using the proposed model (1) or not (0)!
numOfEndmembers = 3;        % number of endmembers
numOfSamples = 100;        % number of pixels
%numOfSamples = 200; 

useJavaProgressMonitor = 1;     % Use java progress monitor in parfor
% useJavaProgressMonitor = 0;   % Use '0' to turn it off.

nlD = 0.5;  % degree of nonlinearity
dF = 6;     %decimation factor
PFA = 0.01; % Probability of False Alarm

detectionRelaxingFactor=0.9;   % relaxing factor for the detection threshold.
maxNumOfRuns =10;              % max number of iterations for iterativeEndmemberEstAndNonlinDetect()

proportions = 0.5;  % proportoin of nonlinear pixels 
Nl = (1-proportions)*numOfSamples;
Nnl = numOfSamples-Nl;

% SK-Hype parameters
C=100;  % reg = 1/C % used C=100 for the PNMM: check if we could use C=100 here too!
kernel='gaussian';
kbw=2;  % kernel bandwidth


%% Generate synthetic data using a linear and a nonlinear model  

modelStr = 'linear';
[Yl,M,Al,noiseVar,gamma,k,spectra_names] = createDecimatedDataFromRealEndMembersSNR_NLD(numOfEndmembers,modelStr,Nl,SNR,[],nlD,dF);

 
modelStr = 'ppnmm';   % this is a type of bilinear model
% modelStr = 'pnmm';   % (Ma)^{\xi}
[Yg,M,Ag,~,gamma,k,spectra_names] = createDecimatedDataFromRealEndMembersSNR_NLD(numOfEndmembers,modelStr,Nnl,SNR,[],nlD,dF);


Y=[Yl,Yg];      % pixels
A = [Al,Ag];    % true abundances

%% Estimate endmember and detect nonlineraly mixed pixels. 
if oestimateM
    [endmembers, nlPixIdx] = iterativeEndmemberEstAndNonlinDetect(Y,numOfEndmembers, PFA, maxNumOfRuns,detectionRelaxingFactor,useJavaProgressMonitor,M);
else
    endmembers = M;
     
    t = getThresholdForGivenPFAUsingLinearImageApproach(Y,M,PFA,numOfEndmembers,useJavaProgressMonitor);
    [stats,Aest,kpar] = computeRatioNLDetectiosStats(Y,M,[],useJavaProgressMonitor);
    %    
    % Detect nonlinear pixels with decision threshold
    nlPixIdx=zeros(numOfSamples,1);
    nlPixIdx(find(stats<=t)) = 1;
end


%% Other detectors! Here you can place a new detector to compare using the same conditions
nlIdx_ppnmm = ppnmmNLDetection(Y,M,PFA);
nlIdx_robustLS = robustLSDetection(Y,M,PFA);




%% Unmixing using linear and nonlinear models
[Y_est_lin, Aest_lin] = estimateLinearModel(Y, endmembers, 1);
[Aest_skp,~,Y_est_skp] = tskHype(Y,endmembers,kernel,C,kbw);
%[Aest_skp] = tKHype(Y,endmembers,kernel,C,kbw);
 

% [Y_est_lin, Aest_lin] = estimateLinearModel(Y, M, 1);
% [Aest_skp,~,Y_est_skp] = tskHype(Y,M);
Aest_det = zeros(size(A));
for i=1:numOfSamples,
    if nlPixIdx(i) ==1
        Aest_det(:,i) = Aest_skp(:,i);
%         Y_est_det(:,i) = Y_est_skp(:,i);        
    else      
        Aest_det(:,i) = Aest_lin(:,i);
%         Y_est_det(:,i) = Y_est_lin(:,i);
    end
end


A_ppnmm =zeros(size(A));
for i=1:numOfSamples,
    if nlIdx_ppnmm(i) ==1
        A_ppnmm(:,i) = Aest_skp(:,i);        
    else      
        A_ppnmm(:,i) = Aest_lin(:,i);        
    end
end

A_rLS =zeros(size(A));
for i=1:numOfSamples,
    if nlIdx_robustLS(i) ==1
        A_rLS(:,i) = Aest_skp(:,i);        
    else      
        A_rLS(:,i) = Aest_lin(:,i);        
    end
end





% Errors for the linearly mixed pixels
[ RMSE_det1LMM, STD_det1LMM] = RMSEAndSTDForMatrix(A(:,1:numOfSamples/2), Aest_det(:,1:numOfSamples/2));    % our method
[ RMSE_lin1LMM, STD_lin1LMM] = RMSEAndSTDForMatrix(A(:,1:numOfSamples/2), Aest_lin(:,1:numOfSamples/2));    % FCLS
[ RMSE_skp1LMM, STD_skp1LMM] = RMSEAndSTDForMatrix(A(:,1:numOfSamples/2), Aest_skp(:,1:numOfSamples/2));    % SKHype
[ RMSE_ppnmm1LMM, STD_ppnmm1LMM] = RMSEAndSTDForMatrix(A(:,1:numOfSamples/2), A_ppnmm(:,1:numOfSamples/2));  %  PPNMM Detect
[ RMSE_rLS1LMM, STD_rLS1LMM] = RMSEAndSTDForMatrix(A(:,1:numOfSamples/2), A_rLS(:,1:numOfSamples/2));      %  robust LS Detect

% Errors for the nonlinearly mixed pixels
[ RMSE_det1NL, STD_det1NL] = RMSEAndSTDForMatrix(A(:,numOfSamples/2 + 1:end), Aest_det(:,numOfSamples/2 + 1:end));  % Our
[ RMSE_lin1NL, STD_lin1NL] = RMSEAndSTDForMatrix(A(:,numOfSamples/2 + 1:end), Aest_lin(:,numOfSamples/2 + 1:end));  % FCLS
[ RMSE_skp1NL, STD_skp1NL] = RMSEAndSTDForMatrix(A(:,numOfSamples/2 + 1:end), Aest_skp(:,numOfSamples/2 + 1:end));  % SKHype
[ RMSE_ppnmm1NL, STD_ppnmm1NL] = RMSEAndSTDForMatrix(A(:,numOfSamples/2 + 1:end), A_ppnmm(:,numOfSamples/2 + 1:end));  % ppnmm detect
[ RMSE_rLS1NL, STD_rLS1NL] = RMSEAndSTDForMatrix(A(:,numOfSamples/2 + 1:end), A_rLS(:,numOfSamples/2 + 1:end));  % robust LS detect

% Errors for the full image
[RMSE_lin,STD_lin] = RMSEAndSTDForMatrix(A, Aest_lin);  % FCLS
[RMSE_skp,STD_skp] = RMSEAndSTDForMatrix(A, Aest_skp);  % SKHype
[RMSE_det,STD_det] = RMSEAndSTDForMatrix(A, Aest_det);  % Ours
[RMSE_ppnmm,STD_ppnmm] = RMSEAndSTDForMatrix(A, A_ppnmm);  % ppnmm detect
[RMSE_rLS,STD_rLS] = RMSEAndSTDForMatrix(A, A_rLS);  % robust LS detect


% Detection errors (our)
detError1LMM = 100*length(find(nlPixIdx(1:numOfSamples/2)==1))/(numOfSamples/2);    % linear part
detError1NL = 100*length(find(nlPixIdx(numOfSamples/2+1:end)==0))/(numOfSamples/2); % nonlinear part
detError = 100*(length(find(nlPixIdx(1:numOfSamples/2)==1)) + length(find(nlPixIdx(numOfSamples/2+1:end)==0)))/numOfSamples;    % full img.

% Detection errors (ppnmm)
detError1LMM_ppnmm = 100*length(find(nlIdx_ppnmm(1:numOfSamples/2)==1))/(numOfSamples/2);    % linear part
detError1NL_ppnmm = 100*length(find(nlIdx_ppnmm(numOfSamples/2+1:end)==0))/(numOfSamples/2); % nonlinear part
detError_ppnmm = (detError1LMM_ppnmm + detError1NL_ppnmm)/2;    % full img.

% Detection errors (rLS)
detError1LMM_rLS = 100*length(find(nlIdx_robustLS(1:numOfSamples/2)==1))/(numOfSamples/2);    % linear part
detError1NL_rLS = 100*length(find(nlIdx_robustLS(numOfSamples/2+1:end)==0))/(numOfSamples/2); % nonlinear part
detError_rLS = (detError1LMM_rLS + detError1NL_rLS)/2;    % full img.



disp('Model & FCLS & Sk-Hype & Ours & robust LS & PPNMM ')
disp(['LMM & ',num2str(RMSE_lin1LMM),' $\pm$ ', num2str(STD_lin1LMM),' & ',num2str(RMSE_skp1LMM),' $\pm$ ', num2str(STD_skp1LMM),' & ',num2str(RMSE_det1LMM),' $\pm$ ', num2str(STD_det1LMM),' $|$ ', num2str(detError1LMM),' & ',num2str(RMSE_rLS1LMM),' $\pm$ ', num2str(STD_rLS1LMM),' $|$ ', num2str(detError1LMM_rLS),' & ',num2str(RMSE_ppnmm1LMM),' $\pm$ ', num2str(STD_ppnmm1LMM),' $|$ ', num2str(detError1LMM_ppnmm),'\\' ])
disp(['NLM & ',num2str(RMSE_lin1NL),' $\pm$ ', num2str(STD_lin1NL),' & ',num2str(RMSE_skp1NL),' $\pm$ ', num2str(STD_skp1NL),' & ',num2str(RMSE_det1NL),' $\pm$ ', num2str(STD_det1NL),' $|$ ', num2str(detError1NL),' & ',num2str(RMSE_rLS1NL),' $\pm$ ', num2str(STD_rLS1NL),' $|$ ', num2str(detError1NL_rLS),' & ',num2str(RMSE_ppnmm1NL),' $\pm$ ', num2str(STD_ppnmm1NL),' $|$ ', num2str(detError1NL_ppnmm),'\\' ])
disp(['F.Img & ',num2str(RMSE_lin),' $\pm$ ', num2str(STD_lin),' & ',num2str(RMSE_skp),' $\pm$ ', num2str(STD_skp),' & ',num2str(RMSE_det),' $\pm$ ', num2str(STD_det),' $|$ ', num2str(detError),' & ',num2str(RMSE_rLS),' $\pm$ ', num2str(STD_rLS),' $|$ ', num2str(detError_rLS),' & ',num2str(RMSE_ppnmm),' $\pm$ ', num2str(STD_ppnmm),' $|$ ', num2str(detError_ppnmm),'\\' ])


