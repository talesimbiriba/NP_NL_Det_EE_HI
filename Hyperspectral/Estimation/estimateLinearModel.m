function [yest, aest] = estimateLinearModel(y, M, respectConstrains  )
% abundanceVector = estimateAbundancesUsingLeastSquares(yPixel,
% mixingMatrix, respctConstrains  )
% 
%   return the abundances "abundanceVector" using Least-Squares.
%
%   arguments: the pixel "yPixel", the Mixing Matix "mixingMatix" and
%   respectConstrains = {true, false} to make aest respect the sum-to-one and positivity constraints. 
%   

if nargin < 3
    respectConstrains = false;
    if nargin < 2
        error('The function estimateAbundancesUsingLeastSquares() must have at least 2 input arguments!');
    end
end

if (respectConstrains)
    % Check Cedric/ ohter papers;  
    aest = constrainedAbundanceEstimation(y,M);
    %normalizing sum(a) = 1;
    %aest = aest.*repmat(1./sum(aest),3,1); 
    
else
    aest = unconstrainedLeastSquaresEstimation(y,M);    
end


yest = M*aest;

end

