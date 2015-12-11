function [varargout] = likErf(hyp, y, mu, s2, inf, i)

% likErf - Error function or cumulative Gaussian likelihood function for binary
% classification or probit regression. The expression for the likelihood is 
%   likErf(t) = (1+erf(t/sqrt(2)))/2 = normcdf(t).
%
% Several modes are provided, for computing likelihoods, derivatives and moments
% respectively, see likFunctions.m for the details. In general, care is taken
% to avoid numerical issues when the arguments are extreme.
% 
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2013-02-05.
%
% See also LIKFUNCTIONS.M.

if nargin<3, varargout = {'0'}; return; end   % report number of hyperparameters
if nargin>1, y = sign(y); y(y==0) = 1; else y = 1; end % allow only +/- 1 values
if numel(y)==0, y = 1; end

if nargin<5                              % prediction mode if inf is not present
  y = y.*ones(size(mu));                                       % make y a vector
  s2zero = 1; if nargin>3, if norm(s2)>0, s2zero = 0; end, end         % s2==0 ?
  if s2zero                                         % log probability evaluation
    [p,lp] = cumGauss(y,mu);
  else                                                              % prediction
    lp = likErf(hyp, y, mu, s2, 'infEP'); p = exp(lp);
  end
  ymu = {}; ys2 = {};
  if nargout>1
    ymu = 2*p-1;                                                % first y moment
    if nargout>2
      ys2 = 4*p.*(1-p);                                        % second y moment
    end
  end
  varargout = {lp,ymu,ys2};
else                                                            % inference mode
  switch inf 
  case 'infLaplace'
    if nargin<6                                             % no derivative mode
      f = mu; yf = y.*f;                            % product latents and labels
      dlp = {}; d2lp = {}; d3lp = {};                         % return arguments
      [p,lp] = cumGauss(y,f);
      if nargout>1                                % derivative of log likelihood
        n_p = gauOverCumGauss(yf,p);
        dlp = y.*n_p;                             % derivative of log likelihood
        if nargout>2                          % 2nd derivative of log likelihood
          d2lp = -n_p.^2 - yf.*n_p;
          if nargout>3                        % 3rd derivative of log likelihood
            d3lp = 2*y.*n_p.^3 +3*f.*n_p.^2 +y.*(f.^2-1).*n_p; 
          end
        end
      end
      varargout = {lp,dlp,d2lp,d3lp};
    else                                                       % derivative mode
      varargout = {[],[],[]};                         % derivative w.r.t. hypers
    end

  case 'infEP'
    if nargin<6                                             % no derivative mode
      z = mu./sqrt(1+s2); dlZ = {}; d2lZ = {};
      [junk,lZ] = cumGauss(y,z);                             % log part function
      if numel(y)>0, z=z.*y; end
      if nargout>1
        if numel(y)==0, y=1; end
        n_p = gauOverCumGauss(z,exp(lZ));
        dlZ = y.*n_p./sqrt(1+s2);                      % 1st derivative wrt mean
        if nargout>2, d2lZ = -n_p.*(z+n_p)./(1+s2); end         % 2nd derivative
      end
      varargout = {lZ,dlZ,d2lZ};
    else                                                       % derivative mode
      varargout = {[]};                                     % deriv. wrt hyp.lik
    end
  
  case 'infVB'
    % naive variational lower bound based on asymptotical properties of lik
    % normcdf(t) -> -(t²-2dt+c)/2 for t->-oo (tight lower bound)
    % the bound has the form: (b+z/ga)*f - f.^2/(2*ga) - h(ga)/2
    d =  0.158482605320942;
    n = numel(s2); b = d*y.*ones(n,1); z = zeros(n,1);
    varargout = {b,z};
  end
end

function [p,lp] = cumGauss(y,f)
  if numel(y)>0, yf = y.*f; else yf = f; end     % product of latents and labels
  p = (1+erf(yf/sqrt(2)))/2;                                        % likelihood
  if nargout>1, lp = logphi(yf); end                            % log likelihood
  
function n_p = gauOverCumGauss(f,p)
  n_p = zeros(size(f));       % safely compute Gaussian over cumulative Gaussian
  ok = f>-5;                            % naive evaluation for large values of f
  n_p(ok) = (exp(-f(ok).^2/2)/sqrt(2*pi)) ./ p(ok); 

  bd = f<-6;                                      % tight upper bound evaluation
  n_p(bd) = sqrt(f(bd).^2/4+1)-f(bd)/2;

  interp = ~ok & ~bd;                % linearly interpolate between both of them
  tmp = f(interp);
  lam = -5-f(interp);
  n_p(interp) = (1-lam).*(exp(-tmp.^2/2)/sqrt(2*pi))./p(interp) + ...
                                                 lam .*(sqrt(tmp.^2/4+1)-tmp/2);