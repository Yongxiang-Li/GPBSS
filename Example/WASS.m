function WASS = WASS(mu0, v0, mu1, v1)
    % mu0, v0: ground truth predictive mean and variance
    % mu1, v1: model predictive mean and variance
    epsv = 1e-12;
    v0 = max(v0(:), epsv * ones(size(v0(:))));
    v1 = max(v1(:), epsv * ones(size(v1(:))));

    sigma0 = sqrt(v0);
    sigma1 = sqrt(v1);
    dmu = mu1(:) - mu0(:);

    wass_vec = sqrt( dmu.^2 + (sigma0 - sigma1).^2 );
    WASS.mean = mean(wass_vec);
    WASS.all = wass_vec;
end
