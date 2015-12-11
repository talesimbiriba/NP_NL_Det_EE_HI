function [ a_est ] = constrainedAbundanceEstimation(Y,M,mu)
%function [ a_est ] = constrainedAbundanceEstimation(Y,M)
%
% Returns the estimated abundances for the observations Y and  
if nargin < 3
    mu = 0.001;
end
    
N = size(Y,2);
R = size(M,2);
a_est= zeros(R,N);

H = M'*M + mu*eye(R);
L = -eye(R);
k = zeros(R,1);
A = ones(1,R);
b = 1;
%options = optimset('TolX',1e-18,'LargeScale','off');

for n=1:N
    f = (-Y(:,n)'*M)';
    a_est(:,n) = qpas(H,f,L,k,A,b);                 
    %a_est(:,n) = quadprog(H,f,L,k,[],[],[],[],[],options);
end

%a_est  = (hyperFcls(Y,M));

end

