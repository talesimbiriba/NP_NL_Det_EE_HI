function [ a_est, d, r_est ] = tskHype(r,M, kernelType, C, kbw, ousePureNLKernel)
%tskHype run the skHype pro algorithm developed by Jie CHEN & Cedric
%Richard
%   [ a_est, d ] = tskHype(r,M, kernelType, C, kbw)
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
if nargin < 6
    ousePureNLKernel = 0;
    if nargin < 5
        kbw =[];
        if nargin < 4
            C = [];
            if nargin < 3
                kernelType = [];
                if nargin < 2
                    error('Usage: [ a_est, d ] = skHype(r,M, kbw, C)')
                end
            end
        end
    end
end

[L,R] = size(M);
N = size(r,2);

% ============  Parameters to tune ======================== 
% Gaussian kernel bandwidth : 
if isempty(kbw)
    kbw = 2; % Cedric's original
    %kbw = 0.5;
    %kbw = 1.1;
%     kbw=0.8;
end
% Regualrization parameter : 
% \mu in the paper = 1/C
if isempty(C)
    C = 100;
end

if isempty(kernelType)
    kernelType = 'gaussian';
end


% ============= Gaussian kernel calculation =================
switch lower(kernelType)
    case {'gaussian'}
        Q=eye(R)/kbw^2;
        MQM=M*Q*M';
        dMQM = diag(MQM);
        KM = exp(-0.5*(dMQM*ones(L,1)'+ones(L,1)*dMQM'-2*M*Q*M'));
    case {'gaussian2'}
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

% [Qm,Dm]= eig(MM);
% Dm = diag(Dm);
% Dm(1:end-R) = 0;
% MM = Qm*diag(Qm)*Qm';

if (ousePureNLKernel)
    % minimum eigen value size considered nonzero
    minEVSize = 1e-5;

    % joint diagonalization following "C.W. Therrien Discrete random signals and 
    % statistical signal processing, pg. 60"
    [Qm,Dm]= eig(MM);
    idx = find(diag(Dm)>minEVSize);
    DmInvSqrt = zeros(size(Dm));
    for i=1:length(idx)
        DmInvSqrt(idx(i),idx(i)) = sqrt(1/Dm(idx(i),idx(i)));
    end
    MMt1 = (DmInvSqrt*Qm')*(Qm*Dm*Qm')*(DmInvSqrt*Qm')';
    %MMt1 = (sqrt(Dm)\Qm')*(Qm*Dm*Qm')*(sqrt(Dm)\Qm')';
    KMt1 = (DmInvSqrt*Qm')*(KM)*(DmInvSqrt*Qm')';
    %KMt1 = (sqrt(Dm)\Qm')*(KM)*(sqrt(Dm)\Qm')';
    [Qk1,Dk1] = eig(KMt1);
    KMDiag =  diag(Qk1'*KMt1*Qk1);
    MMDiag = diag(Qk1'*MMt1*Qk1);

    rho = min(KMDiag(KMDiag > minEVSize));

    KM = KM - rho*MM;
    etmp = eig(KM);
    etmp = min(etmp(etmp<0));
    if ~isempty(etmp)    
        KM = KM - min(etmp(etmp<0))*eye(size(KM));
    end
end


d = [0.5;0.5];

iter_total=0;
ttt=0;

v=0.0;
%tic

for n = 1 : N
 %   if mod(n,10)==0, n, end
    y = r(:,n);
    d = [0.5;0.5];
    K = d(1)*MM+d(2)*KM;
    Kt = [K+1/C*eye(L)];
    H = [Kt,    d(1)*M;
         d(1)*M',  d(1)*eye(R)];        
    H=(H+H')/2+0e-5*eye(L+R);
    f = -[y;zeros(R,1)];
    A = -[zeros(R,L),eye(R)];
    b = zeros(R,1);    

    z0 = qpas(H,f,A,b);
    
    z = z0(1:L);
    gam = z0(L+1:L+R);
    
    ht = (M'*z+gam)*d(1);

    
    for iter = 1 : 15                        % Maximum iteration 15
 
        
        d(1) = 1/(1+sqrt((d(2)^2*z'*KM*z)/(ht'*ht)));
        d(2)=1-d(1);

        K = d(1)*MM+d(2)*KM;
        Kt = [K+1/C*eye(L)];
        H = [Kt,    d(1) *M;
            d(1) *M',  d(1) *eye(R)];        
        H=(H+H')/2+0e-5*eye(L+R);
        z0 = qpas(H,f,A,b); 
        z = z0(1:L);
        gam = z0(L+1:L+R); 
        ht = (M'*z+gam)*d(1); 
        if (iter >1 & norm(d-dp)/R< 1e-4)    % Stop criterion
             break;
        end
        dp=d;
        
    end


  %   h = d(1)*(M'*beta+gam);
     a_est(:,n)=ht/sum(ht);
     r_est(:,n) = M*ht + d(2)*KM*z;
end
%toc
%[RMSE, std] = ErrComput(a(:,1:N),a_est)





end

