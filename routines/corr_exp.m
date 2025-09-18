function [R, D] = corr_exp(theta, gamma, X, Y)
    X = X.*theta';
    if nargin == 3
        XY = X*X';
        D2 = (diag(XY)+diag(XY)') - 2*XY;
        D = D2.^(gamma/2);
        R = exp(-D); 
    else
        Y = (Y.*theta')';
        D2 = ((X.^2)*ones(size(theta))+ones(size(theta'))*(Y.^2)) - 2*X*Y;
        D = D2.^(gamma/2);
        R = exp(-D);
    end
end