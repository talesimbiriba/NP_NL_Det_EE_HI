function [f,K] = covFuncCalc(modelPar,X,y,kernelType)

npar = length(modelPar);
[R,L] = size(X);
if (L < R)
  % disp('Error L < R');
end

K = zeros(L,L);

tempConst = -1/(2*modelPar(2));
% for i=1:L,
%     for j=1:L,        
%         K(i,j) = modelPar(1) * exp(tempConst* ((X(:,i) - X(:,j))'*(X(:,i) - X(:,j))) );         
%     end
% end

% faster implementation than the one above.
% for i=1:L,    
%     for j=i+1:L,
%         K(i,j) = modelPar(1) * exp(tempConst* ((X(:,i) - X(:,j))'*(X(:,i) - X(:,j))));
%         K(j,i) = K(i,j);        
%     end
% end

if strcmp(kernelType,'sqrExponential')

%     %even faster than before. Computes only the upper diagonal matrix.
%     for i=1:L-1,
%         K(i,i+1:end) = modelPar(1) * exp(tempConst* sum((repmat(X(:,i),1,L-i) - X(:,i+1:end)).^2));
%     end
%     K = K + K';
% 
%     %K = K + eye(L)*modelPar(1);
%     K(1:L+1:L*L) = modelPar(1);
% 

    for i=1:L,    
        for j=i+1:L,
            K(i,j) = modelPar(1) * exp(tempConst* ((X(:,i) - X(:,j))'*(X(:,i) - X(:,j))));
            K(j,i) = K(i,j);        
        end
    end
    K(1:L+1:L*L) = modelPar(1);
    
    
    
    if (npar ==3)
        K = K + modelPar(3)*eye(L);
    else
       % K = K + 0.001*eye(L);
    end


    %using the cholesky decomposition to avoid numerical problems with det(K)
    Lc = chol(K,'lower');
    temp = sum(log(diag(Lc)));

    % the 1st minus bellow turns the maximization in a minimization problem!
    %f = - (-0.5 * y'*(K\y) -0.5*log(det(K)) - (L/2) * log(2*pi) );
    %f = - (-0.5 * y'*(Lc'\(Lc\y)) -temp - (L/2) * log(2*pi) );
    f = -0.5 * y'*(Lc'\(Lc\y));
    %f = - (-0.5 * y'*(K\y) -temp - (L/2) * log(2*pi) );
    
else
    error('Kernel type not implemented!');
end