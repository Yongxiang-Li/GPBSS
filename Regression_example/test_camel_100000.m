clearvars -except methods
dbstop if error;
addpath('../routines/');
fields = {'RMSE', 'time'};
p = 2;  n = 100000;  m = 50000;
filename ='Sixhump_GPBSS_100000.mat';
S = 10;
if ~exist(filename, 'file')
    results_GPBSS = [];
    base_par = [5, 5];
    for i = 1:1:25
            new_par = base_par + i;
            result_n = [];
            for s = 1:S
                result_s = struct();
                rng(s);
                X = lhsdesign(n, p);    
                Y = camel6(X);
                Xt = lhsdesign(m, p);    Yt = camel6(Xt);
                % GPBSS
                lob = [0.1*ones(size(Xt,2),1); 0.001];
                upb = [20*ones(size(Xt,2),1); 0.1];
                theta0 = (lob + upb)/2;
                tic;
                dmodelGPBSS = fit_GPBSS(X, Y, @regpoly0, @corr_gauss, init_Bspline(4,new_par ), lob, upb, theta0);
                [dmodelGPBSS.Yhat,dmodelGPBSS.sigma2] = predict_GPBSS(dmodelGPBSS, Xt);
                dmodelGPBSS.time = toc;
                dmodelGPBSS.RMSE = sqrt(mean((dmodelGPBSS.Yhat - Yt).^2));
                for j = 1 : length(fields)
                    result_s.(fields{j}) = [dmodelGPBSS.(fields{j}) ]';
                    result_s.Nk = new_par;
                end
                result_n = [result_n; result_s];
            end
            results_GPBSS = [results_GPBSS; result_n];
            avg_RMSE = mean([result_n(:).RMSE]);    avg_time = mean([result_n(:).time]);  std_RMSE = std([result_n(:).RMSE]);
            disp(['Nk value ', num2str(new_par), ' completed.']);
            disp(['Average RMSE for Nk = ', num2str(new_par), ': ', num2str(avg_RMSE)]);
            disp(['std of RMSE for Nk = ', num2str(new_par), ': ', num2str(std_RMSE)]);
            disp(['Average time for Nk = ', num2str(new_par), ': ', num2str(avg_time), ' seconds']);
    end
    save(filename,'results_GPBSS','-v7.3');
end
Nk_GPBSS = [];  RMSE_GPBSS = [];  time_GPBSS = [];
for i = 1:numel(results_GPBSS)
    Nk_GPBSS = [Nk_GPBSS; results_GPBSS(i).Nk];  RMSE_GPBSS = [RMSE_GPBSS; results_GPBSS(i).RMSE];  time_GPBSS = [time_GPBSS; results_GPBSS(i).time];
end
% Calculate mean time and mean RMSE for each unique Nk for GPBSS
unique_Nk_GPBSS = unique(Nk_GPBSS, 'rows');
mean_time_GPBSS = zeros(size(unique_Nk_GPBSS, 1), 1);
mean_RMSE_GPBSS = zeros(size(unique_Nk_GPBSS, 1), 1);
std_RMSE_GPBSS = zeros(size(unique_Nk_GPBSS, 1), 1);

for i = 1:size(unique_Nk_GPBSS, 1)
    idx = ismember(Nk_GPBSS, unique_Nk_GPBSS(i, :), 'rows');
    mean_time_GPBSS(i) = mean(time_GPBSS(idx));   mean_RMSE_GPBSS(i) = mean(RMSE_GPBSS(idx));  std_RMSE_GPBSS(i) = std(RMSE_GPBSS(idx));
end
figure;
hold on;
colors = {[255/256 69/256 0/256]}; 
patch([mean_time_GPBSS; flipud(mean_time_GPBSS)], ...
      [mean_RMSE_GPBSS + std_RMSE_GPBSS; flipud(mean_RMSE_GPBSS - std_RMSE_GPBSS)], ...
      colors{1}, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
hLine_GPBSS = plot(mean_time_GPBSS, mean_RMSE_GPBSS, 'Color', colors{1}, 'LineWidth', 1);

xlabel('Time / s','Interpreter','latex');
ylabel('RMSE','Interpreter','latex');
legend(hLine_GPBSS, {'GPBSS'}, 'Interpreter','latex','Location','best');
xlim([-5,350]);  ylim([-2,17]);
