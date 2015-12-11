function [detStatistic, a_est, kernelPar, norm_els, norm_eKernel] = computeRatioNLDetectiosStats(r,M, kernelType, useJavaProgressMonitor,inferenceMethod)
%skhypeNLDetector run the skHype pro algorithm developed by Jie CHEN & Cedric
%Richard
%   [ detStatistic, a_est, kernelPar, norm_els, norm_eKernel ] = tskHype(r,M, kernelType, C, kbw)
%
%   r - is the LxN data matrix, where L is the number of bands and N is the
%   number of data samples. 
%   M - is the LxR endmember matrix, where R is the number of endmembers.
%   kernelType - is a string with the kernel to be used, 'gaussian' (default)
%   or 'polynomial'.
%   C - is the regularization parameter (default: C=100).
%   kbw - is the gaussian kernel bandwidth parameter (default: kbw = 2).
%   ousePureNLKernel - boolean "1" use pure nonlinear kernel (default and
%   "0" to use regular kernels.


% Jie CHEN & Cedric Richard
% chen@unice.fr, cedric.richard@unice.fr
% Laboratoire Lagrange, Universite de Nice Sophia-antipolis
% Nice, France
%
if nargin < 5 
    inferenceMethod = {@infExact};
    if nargin < 4
        useJavaProgressMonitor = 1;
        if nargin < 3
            kernelType = [];
            if nargin < 2
                error('Usage: [ a_est, d ] = skHype(r,M, kbw, C)')
            end
        end
    end
end

if isempty(kernelType)
    kernelType = 'gaussian';
end

[L,R] = size(M);
N = size(r,2);

a_est = zeros(R,N);
norm_els = zeros(N,1);
norm_eKernel = zeros(N,1);
detStatistic = zeros(N,1);

% ============  Parameters to tune ======================== 

% Regualrization parameter : 
% \mu in the paper = 1/C


%config. parallel computing
nCores = feature('numCores');
isOpen = matlabpool('size') > 0;
if (~isOpen)
    %matlabpool close;
    matlabpool('open', nCores)
end

if (useJavaProgressMonitor)
    progressStepSize = max([1,floor(N/100)]);

    %pctRunOnAll javaaddpath /home/tales/Work/matlab/HyperspectralCode/MultiThread/ParforProgMonv2/java
    runOnAllJavaMonitorCP;
    ppm = ParforProgMon('Compute Ratio Stats: ', N, progressStepSize, 300, 80);
else
    ppm =0;
    progressStepSize=0;
end

%parfor n = 1 : N
parfor n = 1 : N
%     if mod(n,100)==0, n, end


% these 2 lines were made just to avoid parfor warnings.
C=100;
KM=0;



% ============= Gaussian kernel calculation =================
    switch lower(kernelType)
        case {'gaussian'}
            % Gaussian kernel bandwidth : 

            %kbw = 2;
            covfunc = {@covSEisoU};  % for gaussian kernel
            %covfunc = {@covSEiso};  % for sqrExponential kernel
            hyp = estimateGPMUsingGPML( r(:,n), M , covfunc,1, inferenceMethod);
            kbw = exp(hyp.cov(1));
            C = 1/(exp(hyp.lik))^2;
            kernelPar(:,n) = [kbw;C];


            Q=eye(R)/kbw^2;
            MQM=M*Q*M';
            dMQM = diag(MQM);
            KM = exp(-0.5*(dMQM*ones(L,1)'+ones(L,1)*dMQM'-2*M*Q*M'));
        case {'gaussian2'}
            % Gaussian kernel bandwidth : 
            
            %kbw = 2;
            %covfunc = {@covSEiso};  % for gaussian kernel
            %covfunc = {@covSEiso};  % for sqrExponential kernel
            hyp = estimateGPMUsingGPML( r(:,n), M);

            kbw(1) = exp(2*hyp.cov(2));                                           % signal variance
            kbw(2) = exp(hyp.cov(1));                                 % characteristic length scale
            C = 1/(exp(hyp.lik))^2;


            Q=eye(R)/kbw(2)^2;
            MQM=M*Q*M';
            dMQM = diag(MQM);
            KM = kbw(1)^2*exp(-0.5*(dMQM*ones(L,1)'+ones(L,1)*dMQM'-2*M*Q*M'));
        case {'polynomial2'}
            % For using the polynomial proposed polynomial kernel, remove the
            % comment symbol:
            KM = (1+1/R^2*(M-0.5)*(M-0.5)').^2;
        case {'polynomial'}
            % For using the polynomial proposed polynomial kernel, remove the
            % comment symbol:
            KM = (1+M*M').^2;
        otherwise
            error(['Kernel function ', kernelType, ' is not defined or implemented!'])
    end

    MM =M*M';


    KM = (KM + KM')/2;
    MM = (MM + MM')/2;



    %[rl,a_est(:,n)] = estimateLinearModel(r(:,n),M,'true');
    [rl,a_est(:,n)] = estimateLinearModel(r(:,n),M,0);
    norm_els(n) = norm(r(:,n) - rl ,2)^2;

    Lc = chol(KM+ 1/C*eye(L), 'lower');
    Kinv_r = (Lc'\(Lc\r(:,n)));
    rnl = KM'*Kinv_r;
    norm_eKernel(n) = norm(r(:,n) - rnl,2)^2;
    %norm_eKernel(n) = norm(r(:,n) - KM'*z ,2)^2;
    %norm_eKernel(n) = norm(r(:,n) - KM'*((KM + 1/C*eye(L))\r(:,n)) ,2)^2;

     detStatistic(n) = 2*norm_eKernel(n)/ (norm_eKernel(n) + norm_els(n) );
    
    if (mod(n,progressStepSize)==0 && useJavaProgressMonitor==1)
        ppm.increment();
    end
     
end
%toc
%[RMSE, std] = ErrComput(a(:,1:N),a_est)



% sum-to-one normalization.
%a_est = a_est.*repmat(1./sum(a_est),R,1); 




end

