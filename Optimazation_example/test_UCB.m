clear; clc;
addpath('..\routines\');
rng('default');

OjbFun = @Fun_Hills;
ACFun = @UCB;
p = 2;
numExperiments = 10;
maxIter = 800;
lb = [1 * ones(1, p), 0];
ub = [50 * ones(1, p),0.1];
theta0 = (lb + ub) / 2;
opts = optimoptions('ga', 'PopulationSize', 100, 'PlotFcn', @gaplotbestf);
OptACFun = @(func) ga(func, p, [], [], [], [], [0;0], [1;1], [], [], opts);

fitGP_GPBSS = @(X, Y, theta) fit_GPBSS(X, Y, @regpoly1, @corr_gauss, init_Bspline(4, [30 30]), lb, ub, theta0);
results_GPBSS = zeros(maxIter+1, numExperiments);
iterTimes_GPBSS = zeros(maxIter+1, numExperiments);
expTimes_GPBSS = zeros(1, numExperiments);

for expIdx = 1:numExperiments
    rng(expIdx);
    X0 = generate_data(100, p);
    Y0 = OjbFun(X0);
    
    expStartTime = tic;
    dmodelStartTime = tic;
    dmodel0 = fitGP_GPBSS(X0, Y0, theta0);
    dmodel0.time = toc(dmodelStartTime);

    [design, dmodels] = BO(dmodel0, OjbFun, fitGP_GPBSS, @predict_GPBSS, ACFun, OptACFun, [], maxIter);
    expTimes_GPBSS(expIdx) = toc(expStartTime);
    
    maxY = [dmodels(:).ybest]';
    results_GPBSS(:, expIdx) = maxY;
    for iterIdx = 1:maxIter+1
        iterTimes_GPBSS(iterIdx, expIdx) = dmodels(iterIdx).time;
    end
end

save(['result_GPBSSSO_UCB' datestr(today, 'mmdd') '.mat'], 'results_GPBSS', 'iterTimes_GPBSS', 'expTimes_GPBSS');


iterations = 1:maxIter+1; iterations = iterations(:);
mean_GPBSS = mean(results_GPBSS, 2);     std_GPBSS = std(results_GPBSS, 0, 2);
color_gpbss = [0.85, 0.33, 0.1];
fontname = 'Times New Roman';  fontsize = 10;
figure('Units','centimeters','Position',[0,0,8,5.5]); hold on;
fill([iterations; flipud(iterations)], [mean_GPBSS + std_GPBSS; flipud(mean_GPBSS - std_GPBSS)], color_gpbss, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
h_gpbss = plot(iterations, mean_GPBSS, '-', 'LineWidth', 1.2, 'Color', color_gpbss);
xlabel('Number of Iterations', 'Interpreter', 'latex', 'FontSize', fontsize);
ylabel('Objective Value', 'Interpreter', 'latex', 'FontSize', fontsize);
legend(h_gpbss, {'GPBSS'}, 'Interpreter', 'latex', 'FontSize', fontsize,  'Location', 'southeast');
xlim([1, maxIter+1]);  ylim([14, 20.5]);
grid on;  box on;
set(gca, 'FontName', fontname, 'FontSize', fontsize, 'LineWidth', 1);

figure('Units','centimeters','Position',[0,0,4,4]);
boxplot(expTimes_GPBSS(:),'Labels', {'GPBSS'}, 'Colors', [color_gpbss], 'Symbol', '', 'Widths', 0.5);
ylabel('Runtime (seconds)', 'Interpreter', 'latex', 'FontSize', fontsize);
set(gca, 'FontName', fontname, 'FontSize', fontsize, 'LineWidth', 1);
grid on; box on;

