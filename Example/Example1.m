addpath('../routines/');
clear;    close all;    rng('default');

p = 2;
lob = [1*ones(2,1); 0.001];    upb = [30*ones(2,1); 0.5];
theta0 = (lob + upb)/2;

results = []; % 用于存储结果的矩阵
models = {}; % 用于存储所有GPBSSModel的元胞数组
fields = {'RMSE', 'time','p1'};
filename ='Results_sequential.mat';

S = 10; % 重复次数
if ~exist(filename, 'file')
    cellResults = [];
    for n = 10000:10000:100000
        result_n = [];
        % 生成训练数据
        
        for repeat = 1:S
            rng(repeat);
            X = lhsdesign(n, p);
            Y = peaks(6*X(:,1)-3, 6*X(:,2)-3) + randn(n,1)/10;
            Xt = lhsdesign(2*n, p);
            Yt = peaks(6*Xt(:,1)-3, 6*Xt(:,2)-3);
            result_s = struct();
            tic;
            % 构建GPBSSModel_AIC模型
            GPBSSModel_AIC = fit_GPBSS_AIC(X, Y, @regpoly1, @corr_gauss, init_Bspline(5, 10:30), lob, upb, theta0);
            GPBSSModel_AIC.Yhat = predict_GPBSS(GPBSSModel_AIC, Xt);
            GPBSSModel_AIC.time = toc;
            GPBSSModel_AIC.p1 = [GPBSSModel_AIC.bspline(1).p, GPBSSModel_AIC.bspline(2).p];
            GPBSSModel_AIC.RMSE = sqrt(mean((GPBSSModel_AIC.Yhat - Yt).^2));

            tic;
            % 构建GPBSSModel_tAIC模型
            GPBSSModel_tAIC = fit_GPBSS_tAIC(X, Y, @regpoly1, @corr_gauss, init_Bspline(5, 10:30), lob, upb, theta0);
            GPBSSModel_tAIC.Yhat = predict_GPBSS(GPBSSModel_tAIC, Xt);
            GPBSSModel_tAIC.time = toc;
            GPBSSModel_tAIC.p1 = [GPBSSModel_tAIC.bspline(1).p, GPBSSModel_tAIC.bspline(2).p];
            GPBSSModel_tAIC.RMSE = sqrt(mean((GPBSSModel_tAIC.Yhat - Yt).^2));
            
            for i = 1 : length(fields)
                result_s.(fields{i}) = [GPBSSModel_AIC.(fields{i}); GPBSSModel_tAIC.(fields{i}); ]';
                result_s.n = n;
            end
            result_n = [result_n; result_s];

        end
        cellResults = [cellResults; result_n];
        resultsRMSE = reshape([result_n(:).RMSE],2, S)';
        resultsTime = reshape([result_n(:).time],2, S)';
        resultsp1 = reshape([result_n(:).p1],[], S)';
        resultsN = [result_n(:).n]';
        
        % 输出结果
        disp(['n: ' num2str(n) ', Model: AIC, Average RMSE: ' num2str( mean(resultsRMSE(:, 1))) ', Average Time: ' num2str(mean(resultsTime(:, 1))) ' seconds, p1: ' num2str(mean(resultsp1(:, 1:2), 1))]);
        disp(['n: ' num2str(n) ', Model: tAIC, Average RMSE: ' num2str(mean(resultsRMSE(:, 2))) ', Average Time: ' num2str(mean(resultsTime(:, 2))) ' seconds, p1: ' num2str(mean(resultsp1(:, 3:4), 1))]);

    end
    save(filename,'-v7.3');
end
% 保存结果到文件

