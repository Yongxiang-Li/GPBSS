function [ Y ] = randGP(X, theta)
    N = size(X,1);
    C = corr_gauss0(theta, X);
    L = chol(C+1e-10*eye(N), 'lower');
    Z = randn(N,1);
    Y = L*Z;   
    if nargout==0, figure; plot(X, Y, '.'); end
end
