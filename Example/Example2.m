addpath('../routines/');
clear;    close all;    rng('default');    dbstop if error

p = 2;
lob = [5*ones(p,1); 0.001];    upb = [20*ones(p,1); 0.5];
theta0 = (lob + upb)/2;

n = 900;    RMSE = 0;    S = 100;
filename ='Results_knotsnumber.mat';
for s = 1 : S
    X = lhsdesign(n, p);
    Y = peaks(6*X(:,1)-3, 6*X(:,2)-3) + randn(n,1)/10;
    Xt = lhsdesign(10*n, p);
    Yt = peaks(6*Xt(:,1)-3, 6*Xt(:,2)-3);
    
    GPBSSModels = [];    Ki = 10 : 5 : 30;    Kj = 10 : 5 : 30;
    for ki = Ki
        for kj = Kj
            tic;
            GPBSSModel = fit_GPBSS(X, Y, @regpoly1, @corr_gauss, init_Bspline(5,[ki kj]), lob, upb, theta0);
            GPBSSModel.Yhat = predict_GPBSS(GPBSSModel, Xt);
            GPBSSModel.time = toc;
            GPBSSModel.Nk = [GPBSSModel.bspline(1).p, GPBSSModel.bspline(2).p];
            GPBSSModel.RMSE = sqrt(mean((GPBSSModel.Yhat - Yt).^2));
            GPBSSModels = [GPBSSModels; GPBSSModel];
        end
    end
    RMSE = RMSE + reshape([GPBSSModels(:).RMSE], length(Ki), length(Kj));
end
RMSE = RMSE/S;
save(filename,'-v7.3');    
disp(RMSE);