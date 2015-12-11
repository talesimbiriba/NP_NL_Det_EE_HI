function [ figureHandle ] = plotROCApproximation(M, a, kernelPar, dF, plotStr, figureHandle )
% [ figureHandle ] = plotROCApproximation(M, a, kernelPar, plotStr,
% figureHandle )
%
%   Detailed explanation goes here


if nargin < 6,
    figureHandle = [];
    if nargin < 5
        plotStr = [];
        if nargin <4
            dF =1;
            if nargin < 3
                error('Usage: [ figureHandle ] = plotROCApproximation(M, a, kernelPar, plotStr, figureHandle )')
            end
        end
    end
end
    
if isempty(plotStr)
    plotStr = '-k';
end


mp = kernelPar;


%load simDataGP_ROC

yl = M*a;

%[~, mp ] = allSamplesGPLSRatioTest(Yl(:,end),M)
[~,K] = covFuncCalc(mp,M',yl,'sqrExponential');
Xtr  = M'; Xtst = Xtr;
K1 = zeros(size(K));

L = length(yl);

K1 = zeros(size(K));
for i=1:L,
    for j=1:L,
        K1(i,j) = mp(1) * exp(- 0.5* ((Xtst(:,i) - Xtr(:,j))'*(Xtst(:,i) - Xtr(:,j)))/(mp(2)));
    end
end


H = eye(L) - K1*inv(K);
P = eye(L) - M*((M'*M)\M');

PSI = eye(L)*mp(3);

v = yl;

G = H*H';
rank(G)
mx = trace(G*PSI) + v'*G*v;
my = trace(P*PSI) + v'*P*v;
varx = 2*trace((G*PSI)*(G*PSI)) + 4*v'*G*PSI*G*v;
vary = 2*trace((P*PSI)*(P*PSI)) + 4*v'*P*PSI*P*v;
covxy = 2*trace(G*PSI*P*PSI) + 4*v'*G*PSI*P*v;

meanZeroT =(mx/my)*(1 - covxy/(mx*my) + vary/(my^2))
varZeroT = (mx^2/my^2)*(varx/(mx^2) + vary/(my^2) - 2*covxy/(mx*my) )

v= createDecimatedDataFromRealEndMembersSNR(3,'bilinear',1,1000000,a,3,dF);

% mx = mp(3)*trace(G) + v'*G*v;
% my = mp(3)*trace(P) + v'*P*v;
% varx = mp(3)*mp(3)*2*trace(G^2) + 4*mp(3)*v'*G*G*v;
% vary = mp(3)*mp(3)*2*trace(P^2) + 4*mp(3)*v'*P*v;
% covxy = 2*mp(3)*mp(3)*trace(G*K) + 4*mp(3)*v'*G*K*v;
mx = trace(G*PSI) + v'*G*v;
my = trace(P*PSI) + v'*P*v;
varx = 2*trace((G*PSI)*(G*PSI)) + 4*v'*G*PSI*G*v;
vary = 2*trace((P*PSI)*(P*PSI)) + 4*v'*P*PSI*P*v;
covxy = 2*trace(G*PSI*P*PSI) + 4*v'*G*PSI*P*v;


meanOneT =(mx/my)*(1 - covxy/(mx*my) + vary/(my^2))
varOneT = (mx^2/my^2)*(varx/(mx^2) + vary/(my^2) - 2*covxy/(mx*my) )


t = [0:0.01:1.5];
pfa = normcdf(  (t - meanZeroT)/sqrt(varZeroT));
pd = normcdf( (t-meanOneT)/sqrt(varOneT));

if isempty(figureHandle)
    figureHandle = figure;
else
    figure(figureHandle);
    hold on;
end
set(gca,'FontSize',14)
plot(pfa,pd, plotStr);
%plot(PFA,PD);
grid on;
xlabel('Probability of False Alarm','fontsize',14);
ylabel('Probability of Detection','fontsize',14);
title('ROC Curve','fontsize',14);
set(gca,'XTick',[0:0.1:1])
set(gca,'YTick',[0:0.1:1])

end

