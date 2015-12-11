function [ snorme_gp,Ygp ] = computeSquaredNormOfGPFittError(Y,M)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[L,N] = size(Y);
Ygp = zeros(L,N);
snorme_gp = zeros(N,1);
for i=1:N,
    hyp = estimateGPMUsingGPML( Y(:,i), M);
    covfunc = {@covSEiso};
    likfunc = @likGauss;
    Ygp(:,i) = gp(hyp, @infExact, @meanZero, covfunc, likfunc, M, Y(:,i), M);

    snorme_gp(i)  = norm(Y(:,i)-Ygp(:,i),2)^2;
end


end

