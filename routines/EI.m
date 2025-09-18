function [ret,sigma2] = EI(xstar, y, pred, dmodel)
% Calculates the expected improvement at xstar.
% 
% Inputs:
%       xstar - untried design point
%       y - values of objective function at training points x
%       pred - the prediction of xstar
%       dmodel - DACE model: a struct with the elements
%
%   Outputs:
%       ret - the expected improvement at xstar with training values x and y.
%
% Calculate the GP posterior
[mu,sigma2] = pred(dmodel, xstar);
ybest = max(y);

% The explicit formula of EI
if (sigma2 > 0)
    u = (mu - ybest)/sqrt(sigma2);
    ret = sqrt(sigma2)*(normpdf(u) + u*normcdf(u));
    % change to minimum for ga
    ret = -ret;
else
    ret = 0;
end

end