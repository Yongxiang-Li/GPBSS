clearvars -except methods
dbstop if error;

addpath('../routines/');    
data = readtable('output.csv');
fields = {'RMSE', 'time'};
n = 50000;
test_size = 50000;
rng('default');
filename ='Forest_GPBSS_50000.mat';
S = 10;
if ~exist(filename, 'file')
    results_GPBSS = [];
    base_par = [10, 10, 10];
    for i = 0:1:11
        new_par = base_par + i;
        result_n = [];
        for s = 1:S
            result_s = struct();
            rng(s);
            Xd = data{:, [1, 2, 4]};
            Yd = data{:, [3]};
            Xd = [1 1 1].*(Xd - min(Xd))./(max(Xd) - min(Xd));
            index = false(length(Yd),1); index(randperm(length(Yd),n)) = true;
            X = Xd(index,:); Y = Yd(index, end);
            Xt = Xd(~index,:); Yt = Yd(~index, end);
            test_indices = randperm(size(Xt, 1), test_size);
            Xt = Xt(test_indices, :);
            Yt = Yt(test_indices);
            lob = [1*ones(size(Xt,2),1); 0.001];
            upb = [20*ones(size(Xt,2),1); 0.1];
            theta0 = (lob + upb)/2;
            tic;
            dmodelGPBSS = fit_GPBSS(X, Y, @regpoly1, @corr_gauss, init_Bspline(4, new_par), lob, upb, theta0);
            dmodelGPBSS.Yhat = predict_GPBSS(dmodelGPBSS, Xt);
            dmodelGPBSS.time = toc;
            dmodelGPBSS.RMSE = sqrt(mean((dmodelGPBSS.Yhat - Yt).^2));
            for j = 1 : length(fields)
                result_s.(fields{j}) = [dmodelGPBSS.(fields{j}) ]';
                result_s.Nk = new_par;
            end
            result_n = [result_n; result_s];
        end
        results_GPBSS = [results_GPBSS; result_n];
        avg_RMSE = mean([result_n(:).RMSE]);
        avg_time = mean([result_n(:).time]);
        std_RMSE = std([result_n(:).RMSE]);
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
colors = {[189/256 214/256 159/256]}; 
edge_colors = {[79/256 139/256 87/256]}; 
patch([mean_time_GPBSS; flipud(mean_time_GPBSS)], ...
      [mean_RMSE_GPBSS + std_RMSE_GPBSS; flipud(mean_RMSE_GPBSS - std_RMSE_GPBSS)], ...
      colors{1}, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
hLine_GPBSS = plot(mean_time_GPBSS, mean_RMSE_GPBSS, 'Color', edge_colors{1}, 'LineWidth', 0.75);
