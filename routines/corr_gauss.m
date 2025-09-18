function [ R ] = corr_gauss( theta, X, Y )
% 
    if nargin == 3
        D = theta * (X - Y');
    else
        Y = X;
        D = theta * (X - Y');
    end
    R = exp(-D.^2);
end