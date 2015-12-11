function delta2 = robustTestForNonlinearMixtureDetection( y, M )
%function delta = robustTestForNonlinearMixtureDetection( y, M )
%   
%   Return the l2 squared norm of the estimated error 
%   Inputs:
%           y: Lx1 target pixel
%           M: LxR mixing matrix

mR = M(:,end);
yr = y - mR;

respectConstrains = false;
[yr_est, c] = estimateLinearModelReducedEndmemberMatrix( y, M, respectConstrains);

er = yr-yr_est;

delta2 = er'*er;

end

