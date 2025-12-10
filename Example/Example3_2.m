clear; clc; dbstop if error; rng('default');
addpath('../routines/');  

filename = 'WASS_2_theta=5.mat';
fields = {'WASS', 'RMSE', 'time'};
models = {'GPBSS'};

S = 10;  
p = 2; N = 500; Nt = 1000;
theta_true = 5*ones(p,1); 
summary_table = table();  
Nk_list = {[4,4],[5,5],[6,6],[7,7]};
results_by_Nk = struct(); 
for i = 1:length(Nk_list)
    Nk = Nk_list{i};  results_all = struct();
    for m = 1:length(models);  results_all.(models{m}) = []; end
    for s = 1:S
        rng(s);
        X = sort(lhsdesign(N, p));    Xt = sort(lhsdesign(Nt, p));   
        Y = randGP(X,theta_true); 
        Kxx  = corr_gauss0(theta_true, X ) + sqrt(eps)*eye(N);
        Kxxt = corr_gauss0(theta_true, Xt, X);
        Kxtxt= corr_gauss0(theta_true, Xt) + sqrt(eps)*eye(Nt);
        L = chol(Kxx,'lower'); 
        alpha = L'\(L\Y); 
        mu_t = Kxxt * alpha;
        V = L'\(L\Kxxt'); 
        Sigma_t = Kxtxt - Kxxt * V; Sigma_t = (Sigma_t+Sigma_t')/2;
        S2_t = diag(Sigma_t);
        Yt =  mvnrnd(mu_t, Sigma_t)';
        lob = [0.1; 0.1; 0]; upb = [10;10;0.001]; theta0 = (lob + upb)/3;
        result_s = struct();
        tic;
        d = fit_GPBSS(X, Y, @regpoly0, @corr_gauss, init_Bspline(2,Nk), lob, upb, theta0);
        [d.Yhat,d.sigma2] = predict_GPBSS(d, Xt);
        result_s.time = toc;
        result_s.RMSE = sqrt(mean((d.Yhat - Yt).^2));
        result_s.WASS =  WASS(mu_t, S2_t, d.Yhat, d.sigma2).mean;
        results_all.GPBSS = [results_all.GPBSS; result_s];
    end
    for m = 1:length(models)
        model_name = models{m};
        R = results_all.(model_name);
        WASS_mean    = mean([R(:).WASS]);
        RMSE_mean  = mean([R(:).RMSE]);
        Time_mean  = mean([R(:).time]);
        fprintf('%-10s | WASS=%.6f | RMSE=%.6f | Time=%.2fs\n', ...
            model_name, WASS_mean, RMSE_mean, Time_mean);
        summary_table = [summary_table; table(Nk, {model_name}, WASS_mean, RMSE_mean, Time_mean, ...
    'VariableNames', {'Nk','Model','WASS','RMSE','Time'})];

    end
    results_by_Nk.(['Nk_', num2str(Nk(1)*Nk(2))]) = results_all;

end
save(filename, 'summary_table', 'results_by_Nk', '-v7.3');
figure; hold on;
Nk_list = [16,25,36,49];
model = 'GPBSS';  color = lines(1);  mu = nan(length(Nk_list),1);  sigma = nan(length(Nk_list),1);
for i = 1:length(Nk_list)
    Nk = Nk_list(i);  key = ['Nk_', num2str(Nk)];
    if isfield(results_by_Nk, key)
        result_struct = results_by_Nk.(key);
        if isfield(result_struct, model)
            res_array = result_struct.(model);  wass_vals = [res_array(:).WASS];
            mu(i) = mean(wass_vals, 'omitnan');  sigma(i) = std(wass_vals, 'omitnan');
        end
    end
end
Time_mu = Nk_list(:);  lower = mu - sigma;  upper = mu + sigma;
fill([Time_mu; flipud(Time_mu)],[lower; flipud(upper)],color, 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility','off');
plot(Time_mu, mu, '-o', 'LineWidth', 1,  'Color', color, 'DisplayName', model);
xlabel('Number of inducing point','Interpreter','latex','FontSize',12);
ylabel('Wasserstein Distance','Interpreter','latex','FontSize',12);
legend('Location', 'best');  xticks(Nk_list);  grid on;

colors = lines(length(models));
figure; hold on;
RMSE_mu = nan(length(Nk_list),1);   RMSE_std = nan(length(Nk_list),1);   Time_mu = nan(length(Nk_list),1);
for i = 1:length(Nk_list)
    Nk = Nk_list(i);
    key = ['Nk_', num2str(Nk)];
    if isfield(results_by_Nk, key)
        result_struct = results_by_Nk.(key);
        if isfield(result_struct, model)
            res_array = result_struct.(model);   RMSE_vals = [res_array(:).RMSE];  Time_vals = [res_array(:).time];
            RMSE_mu(i) = mean(RMSE_vals);    RMSE_std(i) = std(RMSE_vals);   Time_mu(i) = mean(Time_vals);
        end
    end
end
fill([Time_mu; flipud(Time_mu)], [RMSE_mu - RMSE_std; flipud(RMSE_mu + RMSE_std)],colors(m,:), 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(Time_mu, RMSE_mu, '-o', 'LineWidth', 1, 'Color', colors(m,:), 'DisplayName', model);

xlabel('Time (second)','Interpreter','latex','FontSize',12);
ylabel('RMSE','Interpreter','latex','FontSize',12);
legend('Location', 'best');
grid on;

