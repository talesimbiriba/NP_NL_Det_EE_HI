function [ d ] = robustLSDetection(Y,M,PFA)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here


[L,N] = size(Y);
R = size(M,2);
% determining threshold!
K = L-R+1;  % degrees of freedom 

eta = chi2inv(1-PFA,K);

Ym = mean(Y')';
ev = eig((Y - repmat(Ym,1,N))*(Y - repmat(Ym,1,N))'/N);
ev = sort(ev,'descend');
estNoisePower = mean(ev(R:end));
%estNoisePower = mean(ev(1:L-R+1));
% 
% RR=R+20;
% estNoisePower = mean(ev(1:L-RR+1));

d = zeros(N,1);

for n=1:N
    delta2 = robustTestForNonlinearMixtureDetection( Y(:,n), M );
    
    if (delta2/estNoisePower) > eta
        d(n) = 1;
    end
end

end

