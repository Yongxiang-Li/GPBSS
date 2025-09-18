function ret=UCB(xstar, pred, dmodel,t)
% Calculates the upper confidence bound at xstar.
% 
% Inputs:
%       xstar - untried design point
%       y - values of objective function at training points x
%       lbx - the lower bound of the input space
%       ubx - the upper bound of the input space
%       pred - the prediction of xstar
%       dmodel - DACE model: a struct with the elements
%
%   Outputs:
%       ret - the upper confidence bound at xstar with training values x and y.
%
% Calculate the GP posterior
[mu,sigma2] = pred(dmodel, xstar);
delta = 0.25; d = size(xstar,2);
beta_t = 2*log(pi^2*t^2/(6*delta)) + 2*d*log(t);
ret = mu + sqrt(beta_t)*sqrt(sigma2);
% change to minimum for ga
ret = -ret;
end