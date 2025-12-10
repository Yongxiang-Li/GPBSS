addpath('../routines/');
clear;    close all;    rng('default');

p = 2;
lob = [1*ones(2,1); 0.001];    upb = [30*ones(2,1); 0.5];
theta0 = (lob + upb)/2;

fields = {'RMSE', 'time'};
filename ='Example1_1.mat';

S = 10; % 重复次数
if ~exist(filename, 'file')
    cellResults = [];
    for n = 10000:10000:100000
        result_n = [];
        
        for repeat = 1:S
            rng(repeat);
            X = lhsdesign(n, p);
            Y = peaks(6*X(:,1)-3, 6*X(:,2)-3) + randn(n,1)/10;
            Xt = lhsdesign(2*n, p);
            Yt = peaks(6*Xt(:,1)-3, 6*Xt(:,2)-3);
            result_s = struct();
            tic;
            GPBSSModel_AIC = fit_GPBSS_AIC(X, Y, @regpoly1, @corr_gauss, init_Bspline(5, 10:30), lob, upb, theta0);
            GPBSSModel_AIC.Yhat = predict_GPBSS(GPBSSModel_AIC, Xt);
            GPBSSModel_AIC.time = toc;
            GPBSSModel_AIC.RMSE = sqrt(mean((GPBSSModel_AIC.Yhat - Yt).^2));

            tic;
            % 构建GPBSSModel_tAIC模型
            GPBSSModel_tAIC = fit_GPBSS_tAIC(X, Y, @regpoly1, @corr_gauss, init_Bspline(5, 10:30), lob, upb, theta0);
            GPBSSModel_tAIC.Yhat = predict_GPBSS(GPBSSModel_tAIC, Xt);
            GPBSSModel_tAIC.time = toc;
            GPBSSModel_tAIC.RMSE = sqrt(mean((GPBSSModel_tAIC.Yhat - Yt).^2));
            
            result_s = struct();
            result_s.RMSE = [GPBSSModel_AIC.RMSE; GPBSSModel_tAIC.RMSE]';
            result_s.time = [GPBSSModel_AIC.time; GPBSSModel_tAIC.time]';
            result_s.n    = n;
            result_s.AIC_time_total   = GPBSSModel_AIC.AIC_time_total;
            result_s.AIC_time_dim     = GPBSSModel_AIC.AIC_time_dim(:)';  
            result_s.tAIC_time_total  = GPBSSModel_tAIC.tAIC_time_total;
            
            result_n = [result_n; result_s];
        end
        
        cellResults = [cellResults; result_n];
        
        resultsRMSE = reshape([result_n(:).RMSE],2, S)';  
        resultsTime = reshape([result_n(:).time],2, S)'; 
        resultsN    = [result_n(:).n]'; 
        AIC_time_total_all  = [result_n(:).AIC_time_total]';   
        tAIC_time_total_all = [result_n(:).tAIC_time_total]';   
        AIC_time_dim_all   = vertcat(result_n(:).AIC_time_dim);
        disp(['n: ' num2str(n) ...
              ', Model: AIC, Average RMSE: ' num2str(mean(resultsRMSE(:, 1))) ...
              ', Average Time: ' num2str(mean(resultsTime(:, 1))) ' seconds, ' ]);
        disp(['    AIC knot selection: total = ' num2str(mean(AIC_time_total_all)) ' seconds, ' ...
              'per-dimension = [' num2str(mean(AIC_time_dim_all, 1)) '] seconds']);
        
        disp(['n: ' num2str(n) ...
              ', Model: tAIC, Average RMSE: ' num2str(mean(resultsRMSE(:, 2))) ...
              ', Average Time: ' num2str(mean(resultsTime(:, 2))) ' seconds, ']);
        disp(['    tAIC knot selection: total = ' num2str(mean(tAIC_time_total_all)) ' seconds ']);
        
    end
    
    % 保存所有结果
    save(filename,'cellResults','-v7.3');
end

