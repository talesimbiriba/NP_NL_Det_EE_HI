function [ hyp ] = estimateGPMUsingGPML( y, X ,covfunc, nparCovFunc, inferenceMethod )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if nargin < 5
    inferenceMethod = [];
    if nargin < 4
        covfunc = [];
    end
end

if isempty(covfunc)
    covfunc = {@covSEiso};
    hyp.cov = zeros(2,1);
else
    hyp.cov = zeros(nparCovFunc,1);
end

if isempty(inferenceMethod)
    inferenceMethod = {@infExact};
end

likfunc = @likGauss;
hyp.lik = log(0.1);


%hyp = minimize(hyp, @gp, -100, @infExact, @meanZero, covfunc, likfunc, X,y);
[~,hyp] = evalc('minimize(hyp, @gp, -100, inferenceMethod, @meanZero, covfunc, likfunc, X,y)');
%[~,hyp] = minimize(hyp, @gp, -100, inferenceMethod, @meanZero, covfunc, likfunc, X,y);
%[~,hyp] = evalc('minimize(hyp, @gp, -100, @infOptTales, @meanZero, covfunc, likfunc, X,y)');
%theta = [hyp.cov;exp(hyp.lik)^2];


end

