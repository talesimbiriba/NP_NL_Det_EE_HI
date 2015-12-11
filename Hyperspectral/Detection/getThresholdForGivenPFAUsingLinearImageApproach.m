function [ threshold ] = getThresholdForGivenPFAUsingLinearImageApproach( Y,M, PFA, R, useJavaProgressMonitor )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

[L,N] = size(Y);

%Aest = unconstrainedLeastSquaresEstimation(Y,M);
Aest = constrainedAbundanceEstimation(Y,M);
%Aest = hyperFcls(Y,M);
% Noise power estimation

Yest = M*Aest;
esqn=0;
for i=1:N,
    esqn = esqn + norm(Y(:,i)-Yest(:,i),2)^2;
end
esqn = esqn/N;

% Generate new Y.
Yest = Yest + randn(size(Y))*sqrt(esqn/L);

 %% Compute the GP detection statistics for the linear Y i.e. Yest

% compute the GP statistics for the new (linear) data


gpScores = computeRatioNLDetectiosStats(Yest,M, 'gaussian',useJavaProgressMonitor );

%% Fit a distribution to the statistics
l = length(PFA);
eta = zeros(l,1);
temp = gpScores(gpScores <=1);
distPar = betafit(temp);

for i=1:l,
    eta(i) = betainv(PFA(i),distPar(1),distPar(2));
end
threshold = eta;
end

