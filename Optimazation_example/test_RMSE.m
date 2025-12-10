clear;    dbstop if error;
addpath('..\routines\');    addpath('..\GPRR\')
% filename ='Hills_data.mat';
% OjbFun = @Fun_Hills;  
% X0 = generate_data(10000,2);    
% Y0 = OjbFun(X0,true);
% save(filename,'X0','Y0','-v7.3');
data = load('Hills_data.mat');
fields = {'RMSE', 'time'};
N = [100,500,1000,2000];
test_size = 5000;
rng('default');
filename ='Hills_RMSE_results.mat';
S = 100;
if ~exist(filename, 'file')
    cellResults = [];
    for n = N
        result_n = [];
        for s = 1:S
            result_s = struct();
            rng(s);
            Xd = data.X0;     Yd = data.Y0; 
            index = false(length(Yd),1); index(randperm(length(Yd),n)) = true;
            X = Xd(index,:); Y = Yd(index, end);
            Xt = Xd(~index,:); Yt = Yd(~index, end);
            test_indices = randperm(size(Xt, 1), test_size);
            Xt = Xt(test_indices, :);   Yt = Yt(test_indices);
            p = size(Xt,2);
            lb = [0.1*ones(1,p), 0.01];    ub = [20*ones(1,p),1];     theta0 = (lb + ub) /2;
            tic;
            dmodelGPBSS = fit_GPBSS(X, Y, @regpoly1, @corr_gauss, init_Bspline(4, [30,30]), lb, ub, theta0);
            dmodelGPBSS.Yhat = predict_GPBSS(dmodelGPBSS, Xt);
            dmodelGPBSS.time = toc;
            dmodelGPBSS.RMSE = sqrt(mean((dmodelGPBSS.Yhat - Yt).^2));
            fprintf('GPBSS - N: %d, Iteration: %d, RMSE: %.4f, Time: %.4f seconds\n', n, s, dmodelGPBSS.RMSE, dmodelGPBSS.time);

            tic;
            dmodelKriging = dacefit(X,Y, @regpoly0, @corr_gauss_s, theta0, lb, ub, 'dace');
            dmodelKriging.Yhat = predictor(dmodelKriging, Xt);
            dmodelKriging.time = toc;
            dmodelKriging.RMSE = sqrt(mean((dmodelKriging.Yhat - Yt).^2));
            fprintf('Kriging - N: %d, Iteration: %d, RMSE: %.4f, Time: %.4f seconds\n', n, s, dmodelKriging.RMSE, dmodelKriging.time);

            for i = 1 : length(fields)
                result_s.(fields{i}) = [dmodelGPBSS.(fields{i}) dmodelKriging.(fields{i})]';
                result_s.n = n;
            end
            result_n = [result_n; result_s];
        end
        cellResults = [cellResults; result_n];
    end
    save(filename,'-v7.3');
end

load('Hills_RMSE_results.mat')
gp_rmse = cell(1, length(N));
kriging_rmse = cell(1, length(N));

% 提取 GPBSS 和 Kriging 的 RMSE 数据
for n_index = 1:length(N)
    n = N(n_index);
    gp_rmse{n_index} = [];
    kriging_rmse{n_index} = [];
    for s = 1:S
        result = cellResults((n_index-1)*S + s);
        gp_rmse{n_index} = [gp_rmse{n_index}; result.RMSE(1)];
        kriging_rmse{n_index} = [kriging_rmse{n_index}; result.RMSE(2)];
    end
end

% 将cell数组转换为矩阵
gp_rmse_all = cell2mat(gp_rmse);
kriging_rmse_all = cell2mat(kriging_rmse);

figure;
h = boxplot(gp_rmse_all);
labelStrings = {'$100$', '$500$', '$1000$', '$2000$'};
set(gca, 'XTickLabel', labelStrings, 'TickLabelInterpreter', 'latex');
ylabel('RMSE', 'Interpreter', 'LaTeX');
xlabel('n','Interpreter', 'LaTeX');
ylim([0 5]);

figure;
h = boxplot(kriging_rmse_all);
labelStrings = {'$100$', '$500$', '$1000$', '$2000$'};
set(gca, 'XTickLabel', labelStrings, 'TickLabelInterpreter', 'latex');
ylabel('RMSE', 'Interpreter', 'LaTeX');
xlabel('n','Interpreter', 'LaTeX');
ylim([0 5]);


