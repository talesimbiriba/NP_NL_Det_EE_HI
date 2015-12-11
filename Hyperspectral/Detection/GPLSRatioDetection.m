function [ detectionIdx, stats, kpar, t] = GPLSRatioDetection(Y, M, PFA,useJavaProgressMonitor)
%function [ detectionIdx, kpar ] = GPLSRatioDetection(Y, M, PFA)
%	Returns the detection index map 'detectionIdx', the detection 
%   statistics 'stats' and kernel parameters 'kpar' for the input 
%   matrix Y and endmembers M for a given PFA.


if nargin <4            
    useJavaProgressMonitor = 1;
end

N = size(Y,2);
R = size(M,2);

[stats,Al,kpar] = computeRatioNLDetectiosStats(Y,M,[],useJavaProgressMonitor);

t = getThresholdForGivenPFAUsingLinearImageApproach(Y,M,PFA,R,useJavaProgressMonitor);

% Detecting
detectionIdx = zeros(N,1);
for n=1:N,
    if (stats(n) < t)
        detectionIdx(n) = 1;
    end
end

end

