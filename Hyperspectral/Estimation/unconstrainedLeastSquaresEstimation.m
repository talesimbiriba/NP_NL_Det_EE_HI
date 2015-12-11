function a  = unconstrainedLeastSquaresEstimation( y, M )
%function a  = unconstrainedLeastSquaresEstimation( y, M )
%   
%   return the model parameters "a" using unconstrained Least-Squares
%
%   inputs: 
%           y: the Lx1 target pixel
%           M: the Mx1 mixing matrix

MTM = M'*M;
c = cond(MTM);
if (c > 1e3) 
    a  = (MTM + 0.0001*eye(size(MTM)))\M'*y;
else
    %a  = (M'*M)\M'*y;
    a  = (MTM)\M'*y;
end

end

