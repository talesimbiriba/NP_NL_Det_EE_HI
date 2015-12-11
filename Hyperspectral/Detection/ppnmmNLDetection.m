function [ d, PFA_emp ] = ppnmmNLDetection(Y,M, PFA)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%   This code uses the PPNMM models to detect nonlinearly mixed pixesl [1].
%
%   Author: Tales. November 2015.
%
%   [1] @article{Altmann-2013-ID308,
% 	Author = {Altmann, Y. and Dobigeon, N. and Tourneret, J.-Y.},
% 	Journal = {{IEEE} Transactions on Image Processing},
% 	Number = {4},
% 	Pages = {1267-1276},
% 	Title = {Nonlinearity detection in hyperspectral images using a polynomial post-nonlinear mixing model.},
% 	Volume = {22},
% 	Year = {2013}}

[L N]=size(Y);
[~,R] = size(M);

% init var
sigma2_est_PPNMM = zeros(N,1);
al = zeros(R,N);
B = zeros(N,1);


% Unmixing procedurels
for i=1:N
    %i
    [al1 b compt]=subgradient_poly2b_u2(Y(:,i),M);
    al(:,i)=al1;
    B(i)=b;
    sigma2_est_PPNMM(i)=1/L*norm(Y(:,i)-gene_poly2(M,al(:,i),B(i)))^2;
end
% 
% save result_unmix_PPNMM_LMM.mat al B sigma2_est_PPNMM

% Detection
%load result_unmix_PPNMM_LMM.mat al B sigma2_est_PPNMM
%PFA=0.05;
for n=1:N
    d(n)=nonlinearity_detection_PPNMM(M,al(:,n),B(n),sigma2_est_PPNMM(n),PFA);
end
% PFA_emp=sum(d)/L

end

