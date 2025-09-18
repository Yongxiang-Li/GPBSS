function [design, trail] = rand_GP( Nt, k, N, p, domain, beta, sigma, phi )
    
    beta = beta(:);   phi = phi(:);
    
    Sdata = domain*slhd(N, k, p);

    % test design sites
    if p == 1
        S = domain * linspace(0, 1, Nt)';
    elseif p == 2
        S = domain * linspace(0, 1, floor(sqrt(Nt)));
        [X1, X2] = meshgrid(S,S);
        S = [reshape(X1,numel(X1),1) reshape(X2,numel(X2),1)];
    else
        S = domain*lhsdesign(Nt, p);
    end
    level = struct('S', {[]}, 'Y', {[]});
    design = repmat(level, k, 1);
    trail = level;
    trail.S = S;

    X = [];
    for j = 1 : k
        design(j).S = Sdata(:,:,j);
        X = [X; Sdata(:,:,j)];
    end
    MU = regpoly0(X)*beta;
    R = corr_gauss_t(phi, X);
    Y = MU + sigma*mvnrnd(zeros(size(MU)),R)';
    dY = reshape(Y, N, k);
    for j = 1 : k
        design(j).Y = dY(:,j);
    end

    Rsx = corr_gauss_t(phi, S, X);
    Rxx = corr_gauss_t(phi, X);
    Rss = corr_gauss_t(phi, S);
    C = chol(Rxx)';
    IRs = Rsx / C';
    MU = regpoly0(S)*beta + IRs * (C \ (Y-regpoly0(X)*beta));
    R = (Rss - IRs * IRs');
    O = MU + sigma*mvnrnd(zeros(size(MU)),R)';
    trail.Y = O;

function [ C ] = corr_gauss_t( phi, X, Y ) % the same as corr_gauss_x
%COV Summary of this function goes here
%   Detailed explanation goes here
    X = X.*repmat(sqrt(phi'), size(X,1), 1);
    if nargin == 2
        Y = X';
        C = 0;
        for p = 1 : size(X,2)
            C = C - (X(:,p)-Y(p,:)).^2;
        end
        C = exp(C);
        C = C + eps*eye(size(C));
    else
        Y = (Y.*repmat(sqrt(phi'), size(Y,1), 1))';
        C = 0;
        for p = 1 : size(X,2)
            C = C - (X(:,p)-Y(p,:)).^2;
        end
        C = exp(C);
    end