
function [ R, D ] = corr_gauss0( theta, X, Y )
    X = X.*theta';
    if nargin == 2
        XY = X*X';
        D = abs((diag(XY)+diag(XY)') - 2*XY);
        R = exp(-D); % C = exp(2*XY-(diag(XY)+diag(XY)'));
    else
        Y = (Y.*theta')';
        D = abs(((X.^2)*ones(size(theta))+ones(size(theta'))*(Y.^2)) - 2*X*Y);
        R = exp(-D); % C = exp(2*X*Y-((X.^2)*ones(size(theta))+ones(size(theta'))*(Y.^2)));
    end
end