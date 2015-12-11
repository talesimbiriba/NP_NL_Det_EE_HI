function [ Y, M, a, noiseVar,gamma,k,spectra_names] = createDecimatedDataFromRealEndMembersSNR_NLD( numOfEndmembers, modelStr, numOfSamples, SNR, abundanceVector, nlDegree, dF)
%
%   [ Y, M, a] = createDacimatedDataFromRealEndMembers( numOfEndmembers, modelStr,
%   numOfSamples, noiseVar, abundanceVector, gamma, dF)
%
%   return the "image" Y built with the "M" mixing matrix and abundances
%   "a"
%   
%
%   numOfEndmembers = number of endmembers
%   modelStr = 'linear','bilinear', ...  
%   numOfSamples = 100 (default)
%   noiseVar  =  0.01 (default)
%   a = vector of the same size as numOfEndmembers (default =
%   randomic)
%   gamma = 0.9 (default)
%   dF = 1; (default decimation factor - no decimation)
%
%   Example:
%   createDataFromRealEndMembers( 3, 'bilinear', 1,
%       0.0001, [], 0.9)
%

if nargin ~= 7
    error('Error: the function createDataFromRealEndMembers must have 7 input args.');
end

if (numOfEndmembers < 1 || numOfEndmembers >8 )
    disp('Error: numOfEndmembers must be in the interval [1,8]');
    return;
end

% if isempty(gamma)
%     gamma = 0.9;
% end

if isempty(numOfSamples)
    numOfSamples = 100;
end

aMat = false;
if isempty(abundanceVector)
    aMat = true;
    A = zeros(numOfEndmembers,numOfSamples);
%     A = rand(numOfEndmembers,numOfSamples);
%     for i=1:numOfSamples,
%         A(:,i) = A(:,i)/sum(A(:,i));
%     end
    count = 1;
    while count <=numOfSamples,
        temp = rand(numOfEndmembers -1,1);
        if sum(temp) < 1
            A(:,count) = [temp; 1-sum(temp)];
            count = count +1;
        end
    end
else
    a = abundanceVector;
    if (length(a) > numOfEndmembers)
        a = a(1:numOfEndmembers);
    end
    a = repmat(a,1,numOfSamples);
end

if size(abundanceVector,2) >1
    a = abundanceVector;
end


if isempty(dF)
    dF = 1;
end



% load JY Endmember Matrix
load MPlus;
M = MPlus(1:dF:end,:);

% load 8endMember_envi_spectralLibrary;
% loads 2 variables, M and spectra_names
%load endmembers.mat

%M = M(1:dF:end,:);

numOfBands = size(M,1);

Y=zeros(numOfBands,numOfSamples);
% 
% M = Mcell{1};
% if numOfEndmembers >1,
%     for i=2:numOfEndmembers,
%         M = [M,Mcell{i}];
%     end
% end

%M = [M(:,1),M(:,3:end)];
%M=M(:,4:end);
M = M(:,1:numOfEndmembers);
%spectra_names = spectra_names(1:numOfEndmembers)';
%spectra_names = spectra_names(4:numOfEndmembers)';
spectra_names = '';


if isempty(SNR)
    noiseVar = 0.001;
else
    noiseVar = noisePowerSNR(M, SNR);
end


%E = norm(mean(M,2),2)^2;
Ar = (1-nlDegree)/nlDegree;
k = sqrt(Ar/(1 + Ar));
%B = (1-k^2)*E;


if (strcmp(lower(modelStr),'linear'))
    if (aMat)
        Y = M*A + randn(numOfBands,numOfSamples)*sqrt(noiseVar);
        a =A;
    else
%         for i=1:numOfSamples,
%             Y(:,i) = M*a + randn(numOfBands,1)*sqrt(noiseVar);            
%         end
        Y = M*a + randn(numOfBands,numOfSamples)*sqrt(noiseVar);        
    end
    k = 1;
    nlDegree = 0;
    gamma=0;
end


if (strcmp(lower(modelStr),'bilinear'))
    if aMat
        a = A;
    end
    linPart = M*a;
       
    nonlinPart = zeros(numOfBands, numOfSamples);
    %compute the nonlinear contributions
    for n=1:numOfSamples,                        
        if numOfEndmembers >1,    
            count = 1;
%             for i=1:numOfEndmembers-1,
%                 for j=1+i:numOfEndmembers;
            for i=1:numOfEndmembers,
                for j=1:numOfEndmembers;

                    M_bilin(:,count) = M(:,i).*M(:,j);
                    a_bilin(count,1) = a(i,n)*a(j,n);
                    count = count+1;
                end
            end
        end
        nonlinPart(:,n) = (M_bilin*a_bilin);                                    
    end
    % Sample mean of the Energies
    El = trace(linPart'*linPart)/numOfSamples;
    Enl = trace(nonlinPart'*nonlinPart)/numOfSamples;
    Elnl = trace(linPart'*nonlinPart)/numOfSamples;
    % compute the mean K
    %K = ( -2*gamma*Elnl + sqrt( 4*(gamma^2)*(Elnl^2) -4*El*( (gamma^2)*Enl-El) ) )/ (2*El);
    B = (1-k^2)*El;
    gamma = (-2*k*Elnl + sqrt(4*(k^2)*Elnl^2 + 4*B*Enl))/(2*Enl);
    % generate the noisy signal.
    %Y = K*linPart + gamma*nonlinPart + randn(numOfBands,numOfSamples)*sqrt(noiseVar);    
    Y = k*linPart + gamma*nonlinPart + randn(numOfBands,numOfSamples)*sqrt(noiseVar);    
    % compute the nonlinearity degree.
    %nlDegree =  (2*K*gamma*Elnl + (gamma^2)*Enl)/( (K^2)*El + 2*K*gamma*Elnl + (gamma^2)*Enl );
end


if (strcmp(lower(modelStr),'ppnmm'))
    if aMat
        a= A;        
    end
    linPart = M*a;
    nonlinPart = linPart.*linPart;
    El = trace(linPart'*linPart)/numOfSamples;
    Enl = trace(nonlinPart'*nonlinPart)/numOfSamples;
    Elnl = trace(linPart'*nonlinPart)/numOfSamples;
    
    B = (1-k^2)*El;
    gamma = (-2*k*Elnl + sqrt(4*(k^2)*Elnl^2 + 4*B*Enl))/(2*Enl);
    Y = k*linPart + gamma*nonlinPart + randn(numOfBands,numOfSamples)*sqrt(noiseVar);    
    %nlDegree =  (2*K*gamma*Elnl + (gamma^2)*Enl)/( (K^2)*El + 2*K*gamma*Elnl + (gamma^2)*Enl );
end
   
if (strcmp(lower(modelStr),'pnmm'))
    if aMat
        a= A;        
    end
    p = 3;    
    linPart = M*a;
    nlPart = (M*a).^p;
    El = trace(linPart'*linPart)/numOfSamples;
    Enl = trace(nlPart'*nlPart)/numOfSamples;
    Elnl = trace(linPart'*nlPart)/numOfSamples;
    
    
    B = (1-k^2)*El;
    gamma = (-2*k*Elnl + sqrt(4*(k^2)*Elnl^2 + 4*B*Enl))/(2*Enl);
      
    Y = k*linPart + gamma*nlPart + randn(numOfBands,numOfSamples)*sqrt(noiseVar);
    
   
    %nld = trace(sconst*2*k*linPart'*nonlinPart + (sconst^2)*nonlinPart'*nonlinPart)/trace((k^2)*linPart'*linPart + sconst*2*k*linPart'*nonlinPart + (sconst^2)*nonlinPart'*nonlinPart);
    nld =  (2*k*gamma*Elnl + (gamma^2)*Enl)/( (k^2)*El + 2*k*gamma*Elnl + (gamma^2)*Enl );
    disp(['NLD = ',num2str(nld)])

end




end

